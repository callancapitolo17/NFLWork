"""In-memory maker state + live-mode reconciliation against Kalshi.

Kalshi is the source of truth for fills/positions (MLB lesson): we poll
/portfolio/fills each tick and /portfolio/positions + /portfolio/settlements
each minute, and a divergence between local bookkeeping and Kalshi positions
is a tripwire, never silently reconciled."""
import datetime
import logging
from zoneinfo import ZoneInfo

from unabated_edge import config
from unabated_edge.maker import ledger, store
from kalshi_common import auth_client

log = logging.getLogger("unabated_edge")

_ET = ZoneInfo("America/New_York")


def trading_day(now: datetime.datetime) -> datetime.date:
    """The trading day a timestamp belongs to, rolled in US-Eastern at
    config.DAILY_ROLL_HOUR_ET (default 6am) rather than UTC midnight --
    see config.py for why (UTC midnight = 8pm ET, mid-slate)."""
    now_et = now.astimezone(_ET)
    shifted = now_et - datetime.timedelta(hours=config.DAILY_ROLL_HOUR_ET)
    return shifted.date()


def _fp(d, key):
    """Numeric field: prefer `<key>_fp` STRING decimal, fall back to bare int."""
    v = d.get(f"{key}_fp")
    if v is None:
        v = d.get(key)
    try:
        return float(v) if v is not None else None
    except (TypeError, ValueError):
        return None


def _price_dollars(d, side):
    v = d.get(f"{side}_price_dollars")
    if v is not None:
        return float(v)
    v = d.get(f"{side}_price")
    return v / 100.0 if v is not None else None


def _money(d, key):
    v = d.get(f"{key}_dollars")
    if v is not None:
        return float(v)
    v = d.get(f"{key}_fp")
    if v is None:
        v = d.get(key)
    if v is None:
        return 0.0
    if isinstance(v, str):
        return float(v)          # string decimals are dollars (e.g. fee_cost "2.168600")
    return float(v) / 100.0      # bare numbers are integer cents (e.g. revenue 16807)


class MakerState:
    def __init__(self):
        self.resting = {}          # (ticker, side) -> {"order_id","price_cents","count"}
        self.ticker_info = {}      # ticker -> {"sport","event_id","line"}
        self.our_orders = {}       # order_id -> (ticker, side)
        self.fills = {}            # event_id -> [(line, side, contracts, price_dollars)]
        self.fills_by_ticker = {}  # ticker -> net yes contracts (positions check)
        self.fill_times = {}       # event_id -> [fill datetimes] (burst tripwire)
        self.cooloff_until = {}    # event_id -> datetime
        self.settled = set()       # tickers already settled (skip in positions check)
        self.position_baseline = {}  # ticker -> Kalshi position adopted at startup (restart-with-inventory)
        self.settled_pnl_today = 0.0
        self._day = None
        self.halt_latched = False  # daily-loss halt latch (finding #75/F6): once tripped,
                                    # stays halted for the rest of this trading day even if
                                    # settled P&L recovers -- clears only on a genuine roll
        self._fills_min_ts = None
        self._done_trades = set()

    def roll_day(self, now):
        day = trading_day(now)
        if self._day != day:
            self._day = day
            self.settled_pnl_today = 0.0
            if self.halt_latched:                      # only touch the DB on an actual clear
                self.halt_latched = False
                store.set_halt_latch(day, False)
                log.info("maker daily halt latch cleared -- new trading day %s", day)

    def latch_halt(self, now):
        """Trip the daily-halt latch (idempotent). Persisted so a restart
        within the same halted trading day comes back up still halted --
        a deliberate same-day resume is an operator restart, not automatic."""
        if not self.halt_latched:                       # only touch the DB on an actual trip
            self.halt_latched = True
            store.set_halt_latch(self._day, True)
            log.warning("maker daily halt latched for trading day %s", self._day)

    def register_ticker(self, ticker, sport, event_id, line):
        self.ticker_info[ticker] = {"sport": sport, "event_id": event_id, "line": line}

    def resting_for(self, ticker, side):
        return self.resting.get((ticker, side))

    def on_place(self, ticker, side, order_id, price_cents, count, mode="quote"):
        self.resting[(ticker, side)] = {"order_id": order_id, "price_cents": price_cents,
                                        "count": count, "mode": mode}
        self.our_orders[order_id] = (ticker, side)

    def on_cancel(self, ticker, side):
        self.resting.pop((ticker, side), None)

    def tickers_for(self, eid):
        return [t for t, i in self.ticker_info.items() if i["event_id"] == eid]

    def events_with_quotes(self):
        return {self.ticker_info[t]["event_id"] for (t, _s) in self.resting if t in self.ticker_info}

    def events_with_exposure(self):
        return set(self.fills) | self.events_with_quotes()

    def quotes_live(self):
        return len(self.resting)

    def quotes_live_for(self, eid):
        return sum(1 for (t, _s) in self.resting
                   if self.ticker_info.get(t, {}).get("event_id") == eid)

    def settled_lines(self, eid):
        """Lines of this event whose ticker has already settled. Used to
        exclude settled fills from the unrealized mark (finding #74/F2/F8):
        a settled fill's true value is already captured exactly in
        settled_pnl_today, so marking it AGAIN against a stale pre-kickoff/
        frozen fair would double-count the same position."""
        return {info["line"] for t, info in self.ticker_info.items()
               if info["event_id"] == eid and t in self.settled}

    def open_fills(self, eid):
        """This event's fills, excluding any already-settled ticker's fill
        (see settled_lines) -- the set that should still be marked
        unrealized against the last-known fair."""
        settled = self.settled_lines(eid)
        return [f for f in self.fills.get(eid, []) if f[0] not in settled]

    def exposure_fills(self, eid, exclude=None):
        """Committed fills PLUS all resting quotes treated as filled, so
        simultaneous fills across rungs can never breach the match cap."""
        out = list(self.fills.get(eid, []))
        for (t, s), r in self.resting.items():
            if (t, s) == exclude:
                continue
            info = self.ticker_info.get(t)
            if info and info["event_id"] == eid:
                out.append((info["line"], s, r["count"], r["price_cents"] / 100.0))
        return out


def poll_fills(state: MakerState, now) -> list[str]:
    min_ts = state._fills_min_ts or int(now.timestamp()) - 3600
    status, body, _ = auth_client.api("GET", f"/portfolio/fills?limit=100&min_ts={min_ts}")
    if status != 200 or not isinstance(body, dict):
        log.warning("maker poll_fills failed: status=%s", status)
        return []
    new = []
    for f in body.get("fills") or []:
        tid, oid = f.get("trade_id"), f.get("order_id")
        if not tid or tid in state._done_trades or oid not in state.our_orders:
            continue
        ticker, placed_side = state.our_orders[oid]
        info = state.ticker_info.get(ticker)
        if info is None:
            continue
        # Trust OUR logical side (from our_orders), not the fill's label: a v2
        # "no" quote is placed as a sell-YES (ask), so the fill reports the YES
        # execution price. Derive the ledger cost on our logical side —
        # no-cost = 1 - yes_exec_price. Fall back to the side-labelled price.
        side = placed_side
        n = _fp(f, "count") or 0.0
        yes_px = _price_dollars(f, "yes")
        if yes_px is not None:
            price = yes_px if side == "yes" else round(1.0 - yes_px, 4)
        else:
            price = _price_dollars(f, side)
        if price is None or n <= 0:
            log.warning("maker fill %s unparseable payload keys=%s", tid, sorted(f))
            continue
        eid = info["event_id"]
        state.fills.setdefault(eid, []).append((info["line"], side, n, price))
        state.fills_by_ticker[ticker] = state.fills_by_ticker.get(ticker, 0.0) + (n if side == "yes" else -n)
        state.fill_times.setdefault(eid, []).append(now)
        state._done_trades.add(tid)
        cur = state.resting.get((ticker, placed_side))
        if cur and cur["order_id"] == oid:
            cur["count"] -= n
            if cur["count"] <= 1e-9:
                state.resting.pop((ticker, placed_side), None)
        worst = ledger.worst_case(state.fills[eid])
        store.log_fill(now, info["sport"], eid, oid, ticker, side, price, n,
                       _money(f, "fee"), worst, tid)
        log.info("maker FILL %s %s %.2f@%.2f worst_after=%.2f", ticker, side, n, price, worst)
        new.append(tid)
    # overlap the next window; trade_id dedup absorbs the re-delivery
    state._fills_min_ts = int(now.timestamp()) - 60
    return new


def poll_positions(state: MakerState, series_prefixes=()) -> bool:
    """True when Kalshi's net positions match baseline + our fills-derived book
    on every relevant ticker (in-series, locally-filled, or baselined; settled
    excluded). An in-series Kalshi position we can't account for is a mismatch —
    that is exactly the orphan shape a restart or a lost fill produces. API
    failure returns True — 'cannot verify' must not false-trip."""
    status, body, _ = auth_client.api("GET", "/portfolio/positions?limit=1000")
    if status != 200 or not isinstance(body, dict):
        log.warning("maker poll_positions failed: status=%s", status)
        return True
    kalshi_pos = {}
    for p in body.get("market_positions") or []:
        t = p.get("ticker") or ""
        if _in_series(t, series_prefixes) or t in state.fills_by_ticker or t in state.position_baseline:
            kalshi_pos[t] = _fp(p, "position") or 0.0
    for t in set(kalshi_pos) | set(state.fills_by_ticker) | set(state.position_baseline):
        if t in state.settled:
            continue
        if (t in state.position_baseline and t not in state.fills_by_ticker
                and abs(kalshi_pos.get(t, 0.0)) < 0.01):
            # baselined ticker we never traded this run went to zero on Kalshi
            # (settled or closed) — retire the baseline instead of tripping forever
            log.info("maker baseline retired for %s (position now flat)", t)
            state.position_baseline.pop(t)
            continue
        expected = state.position_baseline.get(t, 0.0) + state.fills_by_ticker.get(t, 0.0)
        actual = kalshi_pos.get(t, 0.0)
        if abs(actual - expected) > 0.01:
            log.warning("maker position mismatch %s: kalshi=%.2f local=%.2f (baseline=%.2f)",
                        t, actual, expected, state.position_baseline.get(t, 0.0))
            return False
    return True


def poll_settlements(state: MakerState, now):
    status, body, _ = auth_client.api("GET", "/portfolio/settlements?limit=100")
    if status != 200 or not isinstance(body, dict):
        return
    for s in body.get("settlements") or []:
        ticker = s.get("ticker")
        info = state.ticker_info.get(ticker)
        if info is None or ticker in state.settled:
            continue
        pnl = (_money(s, "revenue") - _money(s, "yes_total_cost")
               - _money(s, "no_total_cost") - _money(s, "fee_cost"))
        state.settled.add(ticker)
        state.fills_by_ticker.pop(ticker, None)
        state.roll_day(now)
        state.settled_pnl_today += pnl
        store.log_settlement(now, info["sport"], ticker, pnl)
        log.info("maker SETTLED %s pnl=%.2f day_total=%.2f", ticker, pnl, state.settled_pnl_today)


def _in_series(ticker, prefixes):
    return any(ticker.startswith(p) for p in prefixes)


def restore_halt_latch(state: MakerState, now):
    """Restart-time restore of the daily-halt latch (finding #75/F6): a
    fresh MakerState knows nothing about a previous run, so without this a
    restart during an already-halted trading day would silently un-halt.
    Only restores when the persisted latch belongs to the SAME trading day
    as `now` -- a latch from a prior trading day is stale and must not
    carry forward (that's the roll's job, not a restart's)."""
    state.roll_day(now)
    saved_day, saved_halted = store.get_halt_latch()
    if saved_halted and saved_day == state._day:
        state.halt_latched = True
        log.warning("maker restart: daily halt latch restored for trading day %s", state._day)


def startup_sync(state: MakerState, gateway, series_prefixes):
    """Live-mode startup reconciliation. A fresh MakerState knows nothing about
    a previous run: adopt existing in-series positions as a baseline (so the
    mismatch tripwire doesn't false-trip after a restart with inventory) and
    cancel any in-series resting orders — they are invisible to the fresh cap
    stack and must not stay live."""
    status, body, _ = auth_client.api("GET", "/portfolio/positions?limit=1000")
    if status == 200 and isinstance(body, dict):
        for p in body.get("market_positions") or []:
            t = p.get("ticker") or ""
            pos = _fp(p, "position") or 0.0
            if pos and _in_series(t, series_prefixes):
                state.position_baseline[t] = pos
                log.warning("maker startup: adopting existing position %s=%.2f as baseline", t, pos)
    else:
        log.warning("maker startup: positions fetch failed status=%s (baseline empty)", status)
    sweep_orphan_orders(state, gateway, series_prefixes)


def sweep_orphan_orders(state: MakerState, gateway, series_prefixes) -> int:
    """Cancel in-series resting orders Kalshi knows about but we don't — left
    by a previous run, or created by a POST that errored after landing (never
    retried, next tick gets a fresh client_order_id). Runs at startup and every
    recon cycle, bounding the orphan window to ~60s."""
    status, body, _ = auth_client.api("GET", "/portfolio/orders?status=resting&limit=1000")
    if status != 200 or not isinstance(body, dict):
        log.warning("maker orphan sweep: orders fetch failed status=%s", status)
        return 0
    n = 0
    for o in body.get("orders") or []:
        oid, t = o.get("order_id"), o.get("ticker") or ""
        if not oid or oid in state.our_orders or not _in_series(t, series_prefixes):
            continue
        log.warning("maker orphan order %s on %s — cancelling", oid, t)
        if gateway.cancel(oid):
            n += 1
    return n
