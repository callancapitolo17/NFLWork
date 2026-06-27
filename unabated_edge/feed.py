import datetime
from dataclasses import dataclass, field
import requests
from unabated_edge import config


def line_american_price(ln: dict):
    """American moneyline price from a raw Unabated line dict.
    Field name is americanPrice in the authed feed; some captures show
    'price'. Tolerant of both until Task 0 live recon pins it down."""
    return ln.get("americanPrice", ln.get("price"))

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
            if line.get("isBlurred") is True:
                continue
            st.lines[f"{eid}|{si}|{ms}|{bt}"] = line

def events_for_league(st, league_prefix):
    return [e for e in st.events.values() if e.league_key == league_prefix]

_HEADERS = {"accept":"application/json, text/plain, */*","origin":"https://unabated.com",
    "referer":"https://unabated.com/","user-agent":"Mozilla/5.0 (Macintosh; Intel Mac OS X 10_15_7) "
    "AppleWebKit/537.36 (KHTML, like Gecko) Chrome/147.0.0.0 Safari/537.36"}

def fetch_snapshot(league_prefixes, session=None) -> FeedState:
    r = (session or requests).get(config.UNABATED_SNAPSHOT_URL, headers=_HEADERS, timeout=30)
    r.raise_for_status()
    return parse_snapshot(r.json(), set(league_prefixes))

def fetch_deltas(token, cursor, session=None):
    s = session or requests
    if cursor is None:
        now = datetime.datetime.now(datetime.timezone.utc).replace(microsecond=0).isoformat().replace("+00:00","Z")
        url = f"{config.UNABATED_CHANGES_URL}?full_refresh_ISO={now}"
    else:
        url = f"{config.UNABATED_CHANGES_URL}/{cursor}"
    r = s.get(url, headers=_HEADERS, cookies={"unabated_at_prod": token}, timeout=30)
    r.raise_for_status(); body = r.json()
    ev = [mlc["gameOdds"]["gameOddsEvents"]
          for batch in body.get("results", []) for mlc in batch.get("marketLineChanges", [])
          if "gameOddsEvents" in mlc.get("gameOdds", {})]
    return ev, body.get("latestTimestamp")

def apply_deltas(st, event_dicts, league_prefixes):
    pref = set(league_prefixes)
    for evmap in event_dicts:
        for lk, events in evmap.items():
            if _league_prefix(lk) not in pref: continue
            for ev in events: _ingest(st, ev, lk)
