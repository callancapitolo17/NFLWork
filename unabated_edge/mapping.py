import logging

log = logging.getLogger("unabated_edge")
# Remember which (sport, event_id) we've already warned about so an unmatched
# event in a 2s tick loop is logged once, not thousands of times.
_warned_unmatched: set = set()


def pair_events(adapter, events_meta, kalshi_events) -> list[tuple]:
    """Match Unabated events to Kalshi events by canon team pair.

    Returns a list of (event_meta, kalshi_event) tuples where
    frozenset({canon(home), canon(away)}) == adapter.event_teams(kalshi_event).
    """
    # Index Kalshi events by their canonical team pair
    idx: dict[frozenset, dict] = {}
    for kev in kalshi_events:
        try:
            key = adapter.event_teams(kev)
        except Exception:
            continue
        if key:
            idx[key] = kev

    out = []
    for ev in events_meta:
        key = frozenset({adapter.canon_team(ev.home), adapter.canon_team(ev.away)})
        kev = idx.get(key)
        if kev is None:
            warn_key = (adapter.sport, getattr(ev, "event_id", None))
            if warn_key not in _warned_unmatched:
                _warned_unmatched.add(warn_key)
                log.info(
                    "unmatched %s event %s: %s vs %s (canon key %s) — no Kalshi market; "
                    "available Kalshi team-pair keys: %s",
                    adapter.sport, warn_key[1], ev.home, ev.away, sorted(key),
                    [sorted(k) for k in idx],
                )
            continue
        out.append((ev, kev))
    return out
