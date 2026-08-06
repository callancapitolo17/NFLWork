"""Per-match quoting brain: consumes the tick's devigged ladder + fresh books,
maintains resting quotes through the gateway, enforces the cap stack via the
goal-grid ledger. All decisions are logged to the maker DB (measurement-first)."""
import datetime
import logging
import math

from unabated_edge import config
from unabated_edge.maker import ledger, store
from unabated_edge.venues.kalshi import yes_ask_from_book, no_ask_from_book
from kalshi_common.ev_calc import maker_fee_per_contract

log = logging.getLogger("unabated_edge")
_DT_MIN = datetime.datetime.min.replace(tzinfo=datetime.timezone.utc)


class MakerEngine:
    def __init__(self, gateway, state):
        self.gateway = gateway
        self.state = state
        self.last_success = {}      # sport -> last successful tick (watchdog)
        self._last_skip = {}        # (ticker, side) -> reason (dedup skip rows)
        self._halted = False
        self._fair_by_event = {}    # eid -> {line: p_over}, latest tick

    # ---------- gates ----------

    def _daily_halted(self, now):
        """Finding #75/F6: the halt is a LATCH, not a fresh-each-tick
        recompute. Once realized loss breaches the threshold this trading
        day, stay halted regardless of subsequent P&L recovery -- it only
        clears on a genuine ET trading-day roll (state.roll_day)."""
        self.state.roll_day(now)
        if self.state.settled_pnl_today <= -config.DAILY_LOSS_HALT_PCT * config.BANKROLL:
            self.state.latch_halt(now)
        self._halted = self.state.halt_latched
        return self._halted

    def _hard_stopped(self, now):
        """Finding #74/F2/F8: the realized component (settled_pnl_today)
        is always live -- it reacts to in-play games settling badly even
        though there's no in-play anchor feed to reprice them (runner.py
        stops calling on_match at kickoff, so _fair_by_event is frozen at
        the last pre-kickoff number for the rest of the game; that's the
        best available mark for STILL-OPEN positions, not invented live
        data). Settled fills are excluded from the unrealized mark
        (state.open_fills) so a position isn't counted twice -- once
        correctly via settled_pnl_today, once via a stale phantom mark."""
        realized = self.state.settled_pnl_today
        unreal = sum(ledger.mark_to_fair(self.state.open_fills(e),
                                          self._fair_by_event.get(e, {}))
                     for e in self.state.fills)
        if realized + unreal <= -config.HARD_STOP_DOLLARS:
            self._halted = True
            return True
        return False

    def _fill_burst(self, eid, now):
        times = self.state.fill_times.get(eid, [])
        return len([t for t in times if (now - t).total_seconds() < 60]) > config.FILL_BURST_N

    @staticmethod
    def _parse_mo(s):
        if not s:
            return None
        try:
            dt = datetime.datetime.fromisoformat(s)
        except (ValueError, TypeError):
            return None
        return dt.replace(tzinfo=datetime.timezone.utc) if dt.tzinfo is None else dt

    def _anchor_stale(self, ladder, start_utc, now):
        """Dead-feed guard on the anchor's own modifiedOn — kickoff-aware.

        Near kickoff the sharp total churns constantly, so a rung older than
        ANCHOR_STALE_SEC means the feed is lagging/dark -> pull. Far from
        kickoff (> ANCHOR_STALE_FARK_SEC out) the sharp total legitimately sits
        unchanged for hours; there the poll-success watchdog (MAX_STALENESS_SEC)
        is the real dead-feed guard, so we quote against Unabated's latest even
        if it's old. Always fail stale if NO rung carries a parseable timestamp
        (can't prove provenance)."""
        ages = [(now - dt).total_seconds()
                for dt in (self._parse_mo(r.get("modified_on")) for r in ladder.values())
                if dt is not None]
        if not ages:
            return True
        if start_utc is not None and (start_utc - now).total_seconds() > config.ANCHOR_STALE_FARK_SEC:
            return False
        return min(ages) > config.ANCHOR_STALE_SEC

    # ---------- quote math ----------

    def _margin_cents(self, fair_cents, alt):
        p = max(0.01, fair_cents / 100.0)
        m = max(maker_fee_per_contract(p) * 100 + config.PICKOFF_BUFFER_CENTS,
                config.ROI_MARGIN * fair_cents)
        if alt:
            m *= config.ALT_MARGIN_MULT
        return math.ceil(m)

    def _touch_join_price(self, ticker, side, fair_cents, alt, book, opp_ask_cents):
        """Crowd-touch join price (cents) or None. Entry: net edge within
        [threshold, MAX_MARGIN_CENTS] — same intent as the legacy crowd_tighter
        too-good-to-be-true guard, but fee-inclusive: effective floor here is
        approx fair − 5.5c (5c margin + 0.5c maker fee) vs legacy's fair − 5c.
        Hysteresis: an
        already-resting touch-join holds (queue position is capital) until its
        edge decays below TOUCH_JOIN_EXIT_EDGE_CENTS or turns suspicious. The
        touch is self-excluded: our own resting qty never counts as crowd."""
        bids = book["yes_bids"] if side == "yes" else book["no_bids"]
        cur = self.state.resting_for(ticker, side)

        def edge(price_cents):
            return fair_cents - price_cents - maker_fee_per_contract(price_cents / 100.0) * 100

        touch = None
        for p, q in bids:                       # best-first
            pc = int(round(p * 100))
            if cur and cur["price_cents"] == pc:
                q -= cur["count"]
            if q > 0.5:
                touch = pc
                break
        thr = config.TOUCH_JOIN_ALT_MIN_EDGE_CENTS if alt else config.TOUCH_JOIN_MIN_EDGE_CENTS
        if touch is not None and touch <= opp_ask_cents - 1 and \
                thr <= edge(touch) <= config.MAX_MARGIN_CENTS:
            return touch
        if cur and cur.get("mode") == "touch_join" and cur["price_cents"] <= opp_ask_cents - 1 and \
                config.TOUCH_JOIN_EXIT_EDGE_CENTS <= edge(cur["price_cents"]) <= config.MAX_MARGIN_CENTS:
            return cur["price_cents"]
        return None

    def _global_room(self, eid):
        committed = 0.0
        for other in self.state.events_with_exposure():
            if other != eid:
                committed += max(0.0, -ledger.worst_case(self.state.exposure_fills(other)))
        return config.GLOBAL_CAP_PCT * config.BANKROLL - committed

    def _desired(self, eid, ticker, line, rung, book, side):
        """((price_cents, count, mode), None) or (None, skip_reason)."""
        alt = bool(rung.get("alt"))
        if alt and not (config.ALT_OVERROUND_MIN <= (rung.get("overround") or 0) <= config.ALT_OVERROUND_MAX):
            return None, "alt_overround"
        fair = rung["p_over"] if side == "yes" else 1 - rung["p_over"]
        fair_cents = fair * 100
        opp_ask = yes_ask_from_book(book) if side == "yes" else no_ask_from_book(book)
        if opp_ask is None:
            return None, "no_crowd"
        opp_ask_cents = int(round(opp_ask * 100))
        legacy = math.floor(fair_cents + 1e-9) - self._margin_cents(fair_cents, alt)
        legacy = min(legacy, opp_ask_cents - 1)               # never cross
        join = self._touch_join_price(ticker, side, fair_cents, alt, book, opp_ask_cents)
        if join is not None and join > legacy:
            price, mode = join, "touch_join"
        else:
            price, mode = legacy, "quote"
            if price < fair_cents - config.MAX_MARGIN_CENTS - 1e-9:
                return None, "crowd_tighter"
        if price < 1:
            return None, "price_floor"
        pd = price / 100.0
        n = math.floor(config.MAX_QUOTE_PCT * config.BANKROLL / pd)
        if alt:
            n = math.floor(n * config.ALT_SIZE_MULT)
        n = min(n, config.MAKER_MAX_CONTRACTS)
        budget = min(config.MATCH_CAP_PCT * config.BANKROLL, self._global_room(eid))
        if budget <= 0:
            return None, "global_cap"
        n = min(n, ledger.max_contracts(
            self.state.exposure_fills(eid, exclude=(ticker, side)), line, side, pd, budget))
        if n <= 0:
            return None, "no_room"
        return (price, n, mode), None

    # ---------- order sync ----------

    def _sync(self, sport, eid, ticker, side, desired, reason, fair, alt, now):
        cur = self.state.resting_for(ticker, side)
        if desired is None:
            if cur:
                if self.gateway.cancel(cur["order_id"]):
                    self.state.on_cancel(ticker, side)
                    store.log_quote(now, sport, eid, ticker, side, "cancel", None, None,
                                    fair, None, alt, reason, cur["order_id"])
            elif self._last_skip.get((ticker, side)) != reason:
                store.log_quote(now, sport, eid, ticker, side, "skip", None, None,
                                fair, None, alt, reason, None)
            self._last_skip[(ticker, side)] = reason
            return
        self._last_skip.pop((ticker, side), None)
        price, count, mode = desired
        if cur and cur["price_cents"] == price:
            # hold: queue position is capital. Note `mode` is NOT updated here,
            # so a held touch_join keeps its hysteresis even if `desired`'s mode
            # changed at the same price (attribution drift only, no money impact).
            return
        action = "replace" if cur else "rest"
        if cur:
            if not self.gateway.cancel(cur["order_id"]):
                return                                  # keep local state; retry next tick
            self.state.on_cancel(ticker, side)
        coid = f"uemk-{ticker}-{side}-{price}-{int(now.timestamp())}"
        oid = self.gateway.place(ticker, side, price, count, coid)
        margin = round(fair - price / 100.0, 4)
        if oid is None:
            store.log_quote(now, sport, eid, ticker, side, "skip", price / 100.0, count,
                            fair, margin, alt, "place_failed", None)
            return
        self.state.on_place(ticker, side, oid, price, count, mode)
        store.log_quote(now, sport, eid, ticker, side, action, price / 100.0, count,
                        fair, margin, alt, mode, oid)

    # ---------- pulls ----------

    def _pull_match(self, eid, now, reason):
        for ticker in self.state.tickers_for(eid):
            info = self.state.ticker_info[ticker]
            for side in ("yes", "no"):
                cur = self.state.resting_for(ticker, side)
                if cur and self.gateway.cancel(cur["order_id"]):
                    self.state.on_cancel(ticker, side)
                    store.log_quote(now, info["sport"], eid, ticker, side, "cancel",
                                    None, None, None, None, None, reason, cur["order_id"])

    def pull_all(self, now, reason):
        for eid in list(self.state.events_with_quotes()):
            self._pull_match(eid, now, reason)

    def note_success(self, sport, now):
        self.last_success[sport] = now

    def watchdog(self, now):
        """Feed-staleness guard: if any adapter hasn't completed a successful
        tick within MAX_STALENESS_SEC, our fair is dark — pull everything.
        Called every main-loop iteration, including after tick exceptions."""
        for _sport, ts in self.last_success.items():
            if (now - ts).total_seconds() > config.MAX_STALENESS_SEC:
                if self.state.quotes_live():
                    log.warning("maker feed stale — pulling all quotes")
                self.pull_all(now, "feed_stale")
                return

    def sweep(self, sport, seen_eids, now):
        """Pull quotes on events that stopped appearing in the paired pre-kick
        set (kickoff passed, pairing lost, market closed)."""
        for eid in list(self.state.events_with_quotes()):
            tickers = self.state.tickers_for(eid)
            if tickers and self.state.ticker_info[tickers[0]]["sport"] == sport and eid not in seen_eids:
                self._pull_match(eid, now, "unpaired")

    def stats(self):
        committed = sum(max(0.0, -ledger.worst_case(self.state.exposure_fills(e)))
                        for e in self.state.events_with_exposure())
        # Same settled-fill exclusion as _hard_stopped (finding #74/F2/F8) --
        # this pnl number is what the heartbeat/operator watches, so it must
        # not double-count a settled leg's now-known value.
        unreal = sum(ledger.mark_to_fair(self.state.open_fills(e),
                                         self._fair_by_event.get(e, {}))
                     for e in self.state.fills)
        pnl = round(self.state.settled_pnl_today + unreal, 2)
        return {"quotes_live": self.state.quotes_live(),
                "worst_total": round(committed, 2), "pnl": pnl, "halted": self._halted}

    # ---------- main entry ----------

    def on_match(self, adapter, event_meta, kalshi_event, ladder, books, now):
        eid = event_meta.event_id
        sport = adapter.sport
        if self._daily_halted(now):
            return self.pull_all(now, "daily_halt")
        if self._hard_stopped(now):
            return self.pull_all(now, "hard_stop")
        if event_meta.start_utc is not None and \
                (event_meta.start_utc - now).total_seconds() < config.QUOTE_PULL_MIN * 60:
            return self._pull_match(eid, now, "pull_window")
        if now < self.state.cooloff_until.get(eid, _DT_MIN):
            return self._pull_match(eid, now, "cooloff")
        if self._fill_burst(eid, now):
            self.state.cooloff_until[eid] = now + datetime.timedelta(minutes=config.COOLOFF_MIN)
            log.warning("maker fill burst on event %s — %.0fmin cooloff", eid, config.COOLOFF_MIN)
            return self._pull_match(eid, now, "fill_burst")
        if not ladder:
            return self._pull_match(eid, now, "anchor_gone")
        self._fair_by_event[eid] = {ln: r["p_over"] for ln, r in ladder.items()}
        if self._anchor_stale(ladder, event_meta.start_utc, now):
            return self._pull_match(eid, now, "anchor_stale")
        baseline = self.state.position_baseline
        if baseline:
            for mk in kalshi_event.get("markets", []):
                t = mk.get("ticker")
                if t and abs(baseline.get(t, 0.0)) > 0.01 and t not in self.state.settled:
                    return self._pull_match(eid, now, "baseline_blocked")
        for mk in kalshi_event.get("markets", []):
            book = books.get(mk.get("ticker"))
            if book is None:
                continue
            ya, na = yes_ask_from_book(book), no_ask_from_book(book)
            if ya is not None and na is not None and ya + na < 1 - 2 * maker_fee_per_contract(0.5):
                log.warning("maker crossed book on %s — pulling match", mk.get("ticker"))
                return self._pull_match(eid, now, "crossed_book")
        for mk in kalshi_event.get("markets", []):
            if mk.get("strike_type") != "greater":
                continue
            try:
                line = float(mk.get("floor_strike"))
            except (TypeError, ValueError):
                continue
            rung, ticker = ladder.get(line), mk.get("ticker")
            if rung is None or not ticker:
                continue
            book = books.get(ticker)
            if book is None:
                continue                        # no book this tick: hold existing quotes
            self.state.register_ticker(ticker, sport, eid, line)
            for side in ("yes", "no"):
                fair = rung["p_over"] if side == "yes" else 1 - rung["p_over"]
                desired, reason = self._desired(eid, ticker, line, rung, book, side)
                self._sync(sport, eid, ticker, side, desired, reason, fair,
                           bool(rung.get("alt")), now)
        exp = self.state.exposure_fills(eid)
        store.log_ledger(now, sport, eid, ledger.worst_case(exp),
                         ledger.pnl_grid(exp) if exp else [], self.state.quotes_live_for(eid))
