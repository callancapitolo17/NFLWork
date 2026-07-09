import argparse, datetime, time, signal, threading
from unabated_edge import config, feed, ev, sizing, mapping, storage
from unabated_edge.venues import kalshi
from unabated_edge.sports import registry
from unabated_edge.log_setup import setup_logging
from unabated_edge.risk import tipoff_ok, kill_switch_ok

log = setup_logging()
_running = threading.Event(); _running.set()
# candidates priced since the last heartbeat: >0 proves the full chain works
# (v2 parsed -> anchors unblurred -> ladder built -> Kalshi paired -> asks live);
# 0 with lines>0 exposes the silent-zero failure modes the raw line count hides.
_stats = {"candidates": 0}


def run_tick(adapter, state, kalshi_events, *, now, dry_run=True, ask_fn) -> list[dict]:
    if not kill_switch_ok():
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
    if dropped:
        log.warning("%s: dropped %d line-snapshot rows with NULL price", adapter.sport, dropped)
    storage.snapshot_lines(adapter.sport, snap)
    flagged_sides = {}          # market_ticker -> set of flagged sides (crossed-book tripwire)
    for event_meta, kev in mapping.pair_events(adapter, events, kalshi_events):
        if not tipoff_ok(event_meta.start_utc, config.KICKOFF_CUTOFF_MIN, now):
            continue
        for c in adapter.price_event(state, event_meta, kev):
            ask = ask_fn(c.market_ticker, c.side)
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
    return flagged


def main_loop(dry_run: bool):
    storage.init()
    kalshi.init()
    # v2 per-league polling: the v2 odds file carries UNBLURRED anchors anonymously
    # (no token needed), and the legacy changes/query delta feed does not carry
    # soccer at all — so each tick re-fetches the per-league file from the CDN.
    last_k = 0.0
    kalshi_events = {}
    ticks = 0
    flagged_since_hb = 0
    hb_every = max(1, round(60 / config.V2_POLL_SEC))    # ~60s heartbeat
    ask_fn = lambda t, side: kalshi.best_yes_ask(t) if side == "yes" else kalshi.best_no_ask(t)
    while _running.is_set():
        try:
            now = datetime.datetime.now(datetime.timezone.utc)
            if time.time() - last_k > 30:
                kalshi_events = {a.sport: kalshi.list_events(a.kalshi_series()) for a in registry.ADAPTERS}
                last_k = time.time()
            n_events = n_lines = 0
            for a in registry.ADAPTERS:
                state = feed.fetch_v2(a.league_id, a.league_prefix)
                n_events += len(state.events)
                n_lines += len(state.lines)
                flagged_since_hb += len(run_tick(a, state, kalshi_events.get(a.sport, []),
                                                 now=now, dry_run=dry_run, ask_fn=ask_fn))
            storage.flush()
            ticks += 1
            if ticks % hb_every == 0:            # heartbeat: distinguishes "broken" from "no edges"
                log.info("heartbeat tick=%d events=%d lines=%d kalshi_events=%d candidates_recent=%d flagged_recent=%d",
                         ticks, n_events, n_lines,
                         sum(len(v) for v in kalshi_events.values()), _stats["candidates"], flagged_since_hb)
                _stats["candidates"] = 0
                flagged_since_hb = 0
        except Exception:
            log.exception("tick failed")
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
