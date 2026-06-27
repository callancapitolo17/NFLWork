import datetime
from dataclasses import dataclass, field

@dataclass(frozen=True)
class EventMeta:
    event_id: int; league_key: str; start_utc: datetime.datetime
    home_id: int; away_id: int; home: str | None; away: str | None

@dataclass
class FeedState:
    books: dict = field(default_factory=dict)
    teams: dict = field(default_factory=dict)
    lines: dict = field(default_factory=dict)
    events: dict = field(default_factory=dict)
    cursor: int | None = None

def _dt(s):
    d = datetime.datetime.fromisoformat(s)
    if d.tzinfo is None: d = d.replace(tzinfo=datetime.timezone.utc)
    return d.astimezone(datetime.timezone.utc)

def _league_prefix(lk): return lk.split(":")[0]

def parse_snapshot(raw: dict, league_prefixes: set[str]) -> FeedState:
    st = FeedState()
    for s in raw.get("marketSources", []): st.books[s["id"]] = s.get("name")
    for tid, t in raw.get("teams", {}).items():
        st.teams[tid] = t.get("name") if isinstance(t, dict) else None
    for lk, events in raw.get("gameOddsEvents", {}).items():
        if _league_prefix(lk) not in league_prefixes: continue
        for ev in events: _ingest(st, ev, lk)
    return st

def _ingest(st, ev, lk):
    eid = ev["eventId"]
    et = ev.get("eventTeams", {})
    hid, aid = et.get("1", {}).get("id"), et.get("0", {}).get("id")
    st.events[eid] = EventMeta(eid, _league_prefix(lk), _dt(ev["eventStart"]),
                               hid, aid, st.teams.get(str(hid)), st.teams.get(str(aid)))
    for key, btmap in ev.get("gameOddsMarketSourcesLines", {}).items():
        si, ms, an = (p[2:] for p in key.split(":"))
        if an != "0": continue
        for bt, line in btmap.items():
            st.lines[f"{eid}|{si}|{ms}|{bt}"] = line

def events_for_league(st, league_prefix):
    return [e for e in st.events.values() if e.league_key == league_prefix]
