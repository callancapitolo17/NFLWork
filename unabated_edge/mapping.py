import logging
from dataclasses import dataclass

log = logging.getLogger("unabated_edge")
# Remember which (sport, event_id) we've already warned about so an unmatched
# event in a 2s tick loop is logged once, not thousands of times.
_warned_unmatched: set = set()


@dataclass(frozen=True)
class Pairing:
    event_meta: object
    outcome_tickers: dict      # outcome name -> kalshi market_ticker


def pair_events(adapter, events_meta, kalshi_events) -> list[Pairing]:
    # index kalshi events by canon team pair
    idx = {}
    for kev in kalshi_events:
        tickers = adapter.map_outcome_tickers(kev)
        named = {adapter.canon_team(n): t for n, t in tickers.pop("_named", [])}
        key = frozenset(named)
        idx[key] = (kev, named, tickers)
    out = []
    for ev in events_meta:
        key = frozenset({adapter.canon_team(ev.home), adapter.canon_team(ev.away)})
        hit = idx.get(key)
        if not hit:
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
        _, named, extra = hit
        ot = dict(extra)
        ot["home"] = named.get(adapter.canon_team(ev.home))
        ot["away"] = named.get(adapter.canon_team(ev.away))
        out.append(Pairing(ev, ot))
    return out


def validate(adapter, pairing: Pairing, fair: dict) -> bool:
    if set(fair) != set(adapter.outcomes):
        return False
    if not all(pairing.outcome_tickers.get(o) for o in adapter.outcomes):
        return False
    return abs(sum(fair.values()) - 1.0) < 1e-6
