from dataclasses import dataclass


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
