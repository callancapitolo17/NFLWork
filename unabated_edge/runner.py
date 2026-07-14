import argparse, datetime, time, signal, threading
from unabated_edge import config, feed, ev, sizing, mapping, storage
from unabated_edge.venues import kalshi
from unabated_edge.sports import registry
from unabated_edge.log_setup import setup_logging
from unabated_edge.risk import tipoff_ok, kill_switch_ok
from unabated_edge.maker import engine as maker_engine, gateway as maker_gateway, \
    state as maker_state, store as maker_store

log = setup_logging()
_running = threading.Event(); _running.set()
# candidates priced since the last heartbeat: >0 proves the full chain works
# (v2 parsed -> anchors unblurred -> ladder built -> Kalshi paired -> asks live);
# 0 with lines>0 exposes the silent-zero failure modes the raw line count hides.
# null_dropped: priceless snapshot rows skipped since the last heartbeat.
# books/trades: Kalshi-side capture rows written since the last heartbeat.
_stats = {"candidates": 0, "null_dropped": 0, "books": 0, "trades": 0}


def _market_fp(mk, key):
    """Numeric market field: live API sends `<key>_fp` as a STRING decimal
    (contracts are fractional, e.g. "2084.00" — verified 2026-07-10); fall back
    to a bare integer `<key>` for older payload shapes."""
    v = mk.get(f"{key}_fp")
    if v is None:
        v = mk.get(key)
    try:
        return float(v) if v is not None else None
    except (TypeError, ValueError):
        return None


def run_tick(adapter, state, kalshi_events, *, now, dry_run=True, book_fn, trades_fn=None, maker=None) -> list[dict]:
    if not kill_switch_ok():
        if maker is not None:
            maker.pull_all(now, "kill_switch")   # emergency stop must also clear resting orders
        return []
    flagged = []
    events = [e for e in state.events.values() if e.league_key == adapter.league_prefix]
    # snapshot only pre-kickoff events: the v2 pregame rows can freeze once a match
    # goes live, and frozen rows with later ts would masquerade as the closing line.
    # Stopping at kickoff makes "last snapshot = the close" true by construction.
    event_ids = {e.event_id for e in events
                 if e.start_utc is None or now < e.start_utc}
    # persist line snapshot for these events (CLV backbone)
    snap, dropped = [], 0
    for k, v in state.lines.items():
        if int(k.split("|")[0]) not in event_ids:
            continue
        price = feed.line_american_price(v)
        if price is None:                       # NULL price would poison CLV averages — skip + count
            dropped += 1
            continue
        snap.append({
            "ts": now,
            "event_id": int(k.split("|")[0]),
            "market_source_id": int(k.split("|")[2]),
            "bet_type": k.split("|")[3],
            "side": k.split("|")[1],
            "price": price,
            "points": v.get("points"),
            "modified_on": v.get("modifiedOn"),   # feed-side timestamp: staleness gate + CLV need it
        })
    # priceless lines (books listed but not quoting) are permanent in the v2 file —
    # report the drop count in the heartbeat, not per tick (would spam every 5s)
    _stats["null_dropped"] += dropped
    storage.snapshot_lines(adapter.sport, snap)
    flagged_sides = {}          # market_ticker -> set of flagged sides (crossed-book tripwire)
    seen_eids = set()
    for event_meta, kev in mapping.pair_events(adapter, events, kalshi_events):
        if event_meta.start_utc is not None and now >= event_meta.start_utc:
            continue            # in-play: Kalshi book/trade capture stops at kickoff too
        seen_eids.add(event_meta.event_id)
        # One book fetch per market per tick: persist the full depth (maker
        # design data — spread width, re-centering speed, where flow rests),
        # then reuse the same book for candidate asks below. Halves the REST
        # load vs the old two-fetches-per-rung ask path.
        books, book_rows = {}, []
        for mk in kev.get("markets", []):
            tkr = mk.get("ticker")
            if not tkr:
                continue
            book = book_fn(tkr)
            if book is None:
                continue
            books[tkr] = book
            yes_bid = book["yes_bids"][0] if book["yes_bids"] else (None, None)
            no_bid = book["no_bids"][0] if book["no_bids"] else (None, None)
            try:
                strike = float(mk.get("floor_strike"))
            except (TypeError, ValueError):
                strike = None
            book_rows.append({
                "ts": now, "market_ticker": tkr, "floor_strike": strike,
                "yes_bid": yes_bid[0], "yes_bid_qty": yes_bid[1],
                "no_bid": no_bid[0], "no_bid_qty": no_bid[1],
                "yes_ask": kalshi.yes_ask_from_book(book),
                "no_ask": kalshi.no_ask_from_book(book),
                "volume": _market_fp(mk, "volume"),
                "open_interest": _market_fp(mk, "open_interest"),
                "depth": book,
            })
        storage.snapshot_books(adapter.sport, book_rows)
        _stats["books"] += len(book_rows)
        if trades_fn is not None:
            for tkr in books:
                trades = trades_fn(tkr)
                storage.insert_trades(adapter.sport, trades, now)
                _stats["trades"] += len(trades)
        if maker is not None:
            maker.on_match(adapter, event_meta, kev,
                           adapter.fair_ladder(state, event_meta), books, now)
        if not tipoff_ok(event_meta.start_utc, config.KICKOFF_CUTOFF_MIN, now):
            continue
        for c in adapter.price_event(state, event_meta, kev):
            book = books.get(c.market_ticker)
            ask = kalshi.yes_ask_from_book(book) if c.side == "yes" else kalshi.no_ask_from_book(book)
            if ask is None:
                continue
            _stats["candidates"] += 1
            ev_d, ev_pct = ev.edge_for_yes(c.fair_prob, ask)
            storage.emit(
                "candidate_priced", adapter.sport,
                event_id=event_meta.event_id,
                label=c.label,
                side=c.side,
                fair=c.fair_prob,
                ask=ask,
                ev_pct=ev_pct,
                ev_dollars=ev_d,
                **(c.meta or {}),
            )
            # Gate on BOTH percentage and absolute EV: ev_pct=ev_d/ask inflates on
            # cheap longshots, so require a real per-contract edge too.
            if ev_pct < config.MIN_EV_PCT or ev_d < config.MIN_EV_DOLLARS:
                continue
            n = sizing.kelly_contracts(
                c.fair_prob, ask, config.BANKROLL,
                config.KELLY_FRACTION, config.BANKROLL * config.PER_MATCH_CAP_PCT,
            )
            row = {
                "ts": now,
                "sport": adapter.sport,
                "event_id": event_meta.event_id,
                "market_ticker": c.market_ticker,
                "outcome": c.label,
                "fair_prob": c.fair_prob,
                "yes_ask": ask,
                "ev_pct": ev_pct,
                "kelly_contracts": n,
                "dry_run": dry_run,
            }
            storage.log_flagged(row)
            flagged.append(row)
            sides = flagged_sides.setdefault(c.market_ticker, set())
            sides.add(c.side)
            if sides == {"yes", "no"}:
                # both sides of one rung "+EV" requires yes_ask+no_ask < 1-2*fee:
                # a crossed/stale orderbook, i.e. a data error — never a real double edge.
                log.warning("both sides of %s flagged in one tick — crossed/stale book, "
                            "treat as data error not edge", c.market_ticker)
            log.info(
                "EDGE %s %s %s fair=%.3f ask=%.2f ev=%.1f%% n=%d [DRY]",
                adapter.sport, c.label, c.side, c.fair_prob, ask, ev_pct * 100, n,
            )
    if maker is not None:
        maker.sweep(adapter.sport, seen_eids, now)
    return flagged


def main_loop(dry_run: bool):
    storage.init()
    kalshi.init()
    # v2 per-league polling: the v2 odds file carries UNBLURRED anchors anonymously
    # (no token needed), and the legacy changes/query delta feed does not carry
    # soccer at all — so each tick re-fetches the per-league file from the CDN.
    last_k = 0.0
    last_recon = 0.0
    kalshi_events = {}
    ticks = 0
    flagged_since_hb = 0
    hb_every = max(1, round(60 / config.V2_POLL_SEC))    # ~60s heartbeat
    trades_every = max(1, round(config.TRADES_POLL_SEC / config.V2_POLL_SEC))
    maker = None
    gw = maker_gateway.make_gateway(config.MAKER_MODE, config.MAKER_LIVE_ACK)
    if gw is not None:
        maker_store.init()
        maker = maker_engine.MakerEngine(gw, maker_state.MakerState())
        log.info("maker enabled mode=%s live=%s", config.MAKER_MODE, gw.is_live)
    while _running.is_set():
        try:
            now = datetime.datetime.now(datetime.timezone.utc)
            if time.time() - last_k > 30:
                kalshi_events = {a.sport: kalshi.list_events(a.kalshi_series()) for a in registry.ADAPTERS}
                last_k = time.time()
            # trades poll rides every Nth tick; lookback overlaps 2x so a slow
            # tick can't open a gap (dedup on trade_id absorbs the overlap)
            trades_fn = ((lambda t: kalshi.recent_trades(t, int(config.TRADES_POLL_SEC * 2)))
                         if ticks % trades_every == 0 else None)
            n_events = n_lines = 0
            for a in registry.ADAPTERS:
                state = feed.fetch_v2(a.league_id, a.league_prefix)
                n_events += len(state.events)
                n_lines += len(state.lines)
                flagged_since_hb += len(run_tick(a, state, kalshi_events.get(a.sport, []),
                                                 now=now, dry_run=dry_run,
                                                 book_fn=kalshi.get_book, trades_fn=trades_fn,
                                                 maker=maker))
                if maker:
                    maker.note_success(a.sport, now)
            storage.flush()
            if maker is not None and gw.is_live:
                maker_state.poll_fills(maker.state, now)
                if time.time() - last_recon > 60:
                    if not maker_state.poll_positions(maker.state):
                        maker.pull_all(now, "position_mismatch")
                    maker_state.poll_settlements(maker.state, now)
                    last_recon = time.time()
            ticks += 1
            if ticks % hb_every == 0:            # heartbeat: distinguishes "broken" from "no edges"
                log.info("heartbeat tick=%d events=%d lines=%d kalshi_events=%d candidates_recent=%d flagged_recent=%d null_dropped_recent=%d books_recent=%d trades_recent=%d maker=%s",
                         ticks, n_events, n_lines,
                         sum(len(v) for v in kalshi_events.values()), _stats["candidates"],
                         flagged_since_hb, _stats["null_dropped"], _stats["books"], _stats["trades"],
                         maker.stats() if maker else None)
                _stats["candidates"] = 0
                _stats["null_dropped"] = 0
                _stats["books"] = 0
                _stats["trades"] = 0
                flagged_since_hb = 0
        except Exception:
            log.exception("tick failed")
        if maker is not None:
            maker.watchdog(datetime.datetime.now(datetime.timezone.utc))
        time.sleep(config.V2_POLL_SEC)


def _stop(*_):
    _running.clear()


def cli():
    argparse.ArgumentParser().parse_args()
    signal.signal(signal.SIGINT, _stop)
    signal.signal(signal.SIGTERM, _stop)
    main_loop(dry_run=True)


if __name__ == "__main__":
    cli()
