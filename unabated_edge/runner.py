import argparse, datetime, time, signal, threading
from unabated_edge import config, feed, ev, sizing, mapping, storage
from unabated_edge.venues import kalshi
from unabated_edge.sports import registry
from unabated_edge.log_setup import setup_logging
from unabated_edge.risk import tipoff_ok, kill_switch_ok

log = setup_logging()
_running = threading.Event(); _running.set()


def run_tick(adapter, state, kalshi_events, *, now, dry_run=True, ask_fn) -> list[dict]:
    if not kill_switch_ok():
        return []
    flagged = []
    events = [e for e in state.events.values() if e.league_key == adapter.league_prefix]
    # persist line snapshot for these events (CLV backbone)
    snap = [
        {
            "ts": now,
            "event_id": int(k.split("|")[0]),
            "market_source_id": int(k.split("|")[2]),
            "bet_type": k.split("|")[3],
            "side": k.split("|")[1],
            "price": feed.line_american_price(v),
            "points": v.get("points"),
        }
        for k, v in state.lines.items()
        if int(k.split("|")[0]) in {e.event_id for e in events}
    ]
    storage.snapshot_lines(adapter.sport, snap)
    for p in mapping.pair_events(adapter, events, kalshi_events):
        fair = adapter.fair(state, p.event_meta)
        if fair is None or not mapping.validate(adapter, p, fair):
            continue
        if not tipoff_ok(p.event_meta.start_utc, config.KICKOFF_CUTOFF_MIN, now):
            continue
        for outcome in adapter.outcomes:
            ticker = p.outcome_tickers[outcome]
            ask = ask_fn(ticker)
            if ask is None:
                continue
            ev_d, ev_pct = ev.edge_for_yes(fair[outcome], ask)
            storage.emit(
                "candidate_priced", adapter.sport,
                event_id=p.event_meta.event_id,
                outcome=outcome,
                fair=fair[outcome],
                ask=ask,
                ev_pct=ev_pct,
            )
            if ev_pct < config.MIN_EV_PCT:
                continue
            n = sizing.kelly_contracts(
                fair[outcome], ask, config.BANKROLL,
                config.KELLY_FRACTION, config.BANKROLL * config.PER_MATCH_CAP_PCT,
            )
            row = {
                "ts": now,
                "sport": adapter.sport,
                "event_id": p.event_meta.event_id,
                "market_ticker": ticker,
                "outcome": outcome,
                "fair_prob": fair[outcome],
                "yes_ask": ask,
                "ev_pct": ev_pct,
                "kelly_contracts": n,
                "dry_run": dry_run,
            }
            storage.log_flagged(row)
            flagged.append(row)
            log.info(
                "EDGE %s %s fair=%.3f ask=%.2f ev=%.1f%% n=%d [DRY]",
                adapter.sport, outcome, fair[outcome], ask, ev_pct * 100, n,
            )
    return flagged


def main_loop(dry_run: bool):
    storage.init()
    if not config.UNABATED_TOKEN:
        log.warning(
            "UNABATED_AT_PROD is not set — feed will return blurred/anonymous lines "
            "and no real edges will be found until a token is set in .env"
        )
    kalshi.init()
    prefixes = registry.league_prefixes()
    state = feed.fetch_snapshot(prefixes)
    cursor = None
    last_k = 0.0
    kalshi_events = {}
    while _running.is_set():
        try:
            evs, cursor = feed.fetch_deltas(config.UNABATED_TOKEN, cursor)
            feed.apply_deltas(state, evs, prefixes)
            now = datetime.datetime.now(datetime.timezone.utc)
            if time.time() - last_k > 30:
                kalshi_events = {a.sport: kalshi.list_events(a.kalshi_series()) for a in registry.ADAPTERS}
                last_k = time.time()
            for a in registry.ADAPTERS:
                run_tick(a, state, kalshi_events.get(a.sport, []), now=now, dry_run=dry_run, ask_fn=kalshi.best_yes_ask)
            storage.flush()
        except Exception:
            log.exception("tick failed")
        time.sleep(2)


def _stop(*_):
    _running.clear()


def cli():
    argparse.ArgumentParser().parse_args()
    signal.signal(signal.SIGINT, _stop)
    signal.signal(signal.SIGTERM, _stop)
    main_loop(dry_run=True)


if __name__ == "__main__":
    cli()
