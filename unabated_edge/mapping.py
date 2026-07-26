import logging
from zoneinfo import ZoneInfo

from unabated_edge.sports.base import SportAdapter

log = logging.getLogger("unabated_edge")
# Remember which (sport, event_id) we've already warned about so an unmatched
# event in a 2s tick loop is logged once, not thousands of times.
_warned_unmatched: set = set()

# Kalshi's MLB event_ticker carries the US-local game date. A 7:05pm ET game
# is 23:05Z (or after-midnight UTC for later starts), so comparing it against
# a naive UTC .date() would misdate evening games by a day -- exactly the
# games most likely to matter. Convert the Unabated side to this zone before
# taking .date() so both sides compare the same US-local calendar day.
_EASTERN = ZoneInfo("America/New_York")


def _adapter_is_date_aware(adapter) -> bool:
    """True when this adapter overrides SportAdapter.event_date (MLB) rather
    than inheriting the "no date support" default (soccer). Decided at the
    class level, not per-event, so a sport either always disambiguates by
    date or never does -- a single unparseable ticker on a date-aware sport
    fails closed by exclusion (see pair_events), it doesn't silently drop
    the whole sport back to legacy team-pair-only matching."""
    return type(adapter).event_date is not SportAdapter.event_date


def _et_date(event_meta):
    return event_meta.start_utc.astimezone(_EASTERN).date()


def pair_events(adapter, events_meta, kalshi_events) -> list[tuple]:
    """Match Unabated events to Kalshi events by canon team pair (and, for
    date-aware adapters, US-local game date).

    Returns a list of (event_meta, kalshi_event) tuples where
    frozenset({canon(home), canon(away)}) == adapter.event_teams(kalshi_event)
    (and, when date-aware, adapter.event_date(kalshi_event) ==
    the Unabated event's Eastern-time date).

    Fails closed on doubleheaders / ambiguous team-pairs: nothing downstream
    cross-checks game start times, so if a (team-pair[, date]) key maps to
    more than one Kalshi event, or more than one Unabated event shares that
    key, there is no safe way to tell which Kalshi game belongs to which
    Unabated game. Rather than guess (which would quote a live market off
    the wrong game's anchor), that whole key is excluded from pairing this
    tick and logged once as a WARNING. No quoting on true doubleheaders
    (same teams, same date) is the accepted cost -- date-aware pairing only
    fixes the far more common case of the same two teams meeting again on a
    later day in a multi-game series, which used to collide the same way.
    """
    date_aware = _adapter_is_date_aware(adapter)

    # Index Kalshi events by their (canonical team pair[, date]) key, keeping
    # ALL events per key (not last-write-wins) so a collision is visible.
    kalshi_by_key: dict = {}
    for kev in kalshi_events:
        try:
            team_key = adapter.event_teams(kev)
        except Exception:
            continue
        if not team_key:
            continue
        if date_aware:
            try:
                date_key = adapter.event_date(kev)
            except Exception:
                date_key = None
            if date_key is None:
                # Unparseable ticker on a date-aware sport: fail closed by
                # exclusion (never surfaces as a match), not by falling
                # back to blind team-pair-only pairing.
                continue
            full_key = (team_key, date_key)
        else:
            full_key = team_key
        kalshi_by_key.setdefault(full_key, []).append(kev)

    # Same grouping on the Unabated side, so a collision there is caught too.
    events_by_key: dict = {}
    for ev in events_meta:
        team_key = frozenset({adapter.canon_team(ev.home), adapter.canon_team(ev.away)})
        full_key = (team_key, _et_date(ev)) if date_aware else team_key
        events_by_key.setdefault(full_key, []).append(ev)

    out = []
    for full_key, evs in events_by_key.items():
        kevs = kalshi_by_key.get(full_key, [])
        team_key = full_key[0] if date_aware else full_key
        if len(evs) > 1 or len(kevs) > 1:
            warn_key = (adapter.sport, "ambiguous", full_key)
            if warn_key not in _warned_unmatched:
                _warned_unmatched.add(warn_key)
                if date_aware:
                    log.warning(
                        "doubleheader/ambiguous team-pair %s on %s: %d kalshi events, "
                        "%d unabated events — excluded (fail closed)",
                        sorted(team_key), full_key[1], len(kevs), len(evs),
                    )
                else:
                    log.warning(
                        "doubleheader/ambiguous team-pair %s: %d kalshi events, %d unabated "
                        "events — excluded (fail closed)",
                        sorted(team_key), len(kevs), len(evs),
                    )
            continue
        if not kevs:
            ev = evs[0]
            warn_key = (adapter.sport, getattr(ev, "event_id", None))
            if warn_key not in _warned_unmatched:
                _warned_unmatched.add(warn_key)
                log.info(
                    "unmatched %s event %s: %s vs %s (canon key %s) — no Kalshi market; "
                    "available Kalshi keys: %s",
                    adapter.sport, warn_key[1], ev.home, ev.away, sorted(team_key),
                    list(kalshi_by_key),
                )
            continue
        out.append((evs[0], kevs[0]))
    return out
