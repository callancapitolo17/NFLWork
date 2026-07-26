import logging

log = logging.getLogger("unabated_edge")
# Remember which (sport, event_id) we've already warned about so an unmatched
# event in a 2s tick loop is logged once, not thousands of times.
_warned_unmatched: set = set()


def pair_events(adapter, events_meta, kalshi_events) -> list[tuple]:
    """Match Unabated events to Kalshi events by canon team pair.

    Returns a list of (event_meta, kalshi_event) tuples where
    frozenset({canon(home), canon(away)}) == adapter.event_teams(kalshi_event).

    Fails closed on doubleheaders / ambiguous team-pairs: nothing downstream
    cross-checks game start times, so if a team-pair key maps to more than
    one Kalshi event, or more than one Unabated event shares that team-pair,
    there is no safe way to tell which Kalshi game belongs to which Unabated
    game. Rather than guess (which would quote a live market off the wrong
    game's anchor), that whole team-pair is excluded from pairing this tick
    and logged once as a WARNING. No quoting on doubleheaders is the
    accepted cost of v1.
    """
    # Index Kalshi events by their canonical team pair, keeping ALL events
    # per key (not last-write-wins) so a doubleheader collision is visible.
    kalshi_by_key: dict[frozenset, list[dict]] = {}
    for kev in kalshi_events:
        try:
            key = adapter.event_teams(kev)
        except Exception:
            continue
        if key:
            kalshi_by_key.setdefault(key, []).append(kev)

    # Same grouping on the Unabated side, so a doubleheader there is caught too.
    events_by_key: dict[frozenset, list] = {}
    for ev in events_meta:
        key = frozenset({adapter.canon_team(ev.home), adapter.canon_team(ev.away)})
        events_by_key.setdefault(key, []).append(ev)

    out = []
    for key, evs in events_by_key.items():
        kevs = kalshi_by_key.get(key, [])
        if len(evs) > 1 or len(kevs) > 1:
            warn_key = (adapter.sport, "ambiguous", key)
            if warn_key not in _warned_unmatched:
                _warned_unmatched.add(warn_key)
                log.warning(
                    "doubleheader/ambiguous team-pair %s: %d kalshi events, %d unabated "
                    "events — excluded (fail closed)",
                    sorted(key), len(kevs), len(evs),
                )
            continue
        if not kevs:
            ev = evs[0]
            warn_key = (adapter.sport, getattr(ev, "event_id", None))
            if warn_key not in _warned_unmatched:
                _warned_unmatched.add(warn_key)
                log.info(
                    "unmatched %s event %s: %s vs %s (canon key %s) — no Kalshi market; "
                    "available Kalshi team-pair keys: %s",
                    adapter.sport, warn_key[1], ev.home, ev.away, sorted(key),
                    [sorted(k) for k in kalshi_by_key],
                )
            continue
        out.append((evs[0], kevs[0]))
    return out
