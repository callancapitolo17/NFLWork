import datetime
import logging
from dataclasses import dataclass, field
import requests
from unabated_edge import config, pricing

log = logging.getLogger("unabated_edge")


def line_american_price(ln: dict):
    """American moneyline price from a raw Unabated line dict.
    Field name is americanPrice in the authed feed; some captures show
    'price'. Tolerant of both until Task 0 live recon pins it down.

    Fail-closed shape check (Finding #76): American odds are never in
    (-100, 100). This is the choke point for that invariant — a malformed
    value (e.g. a stray decimal-odds price like 1.91) is rejected here
    (None, logged once) so it never reaches pricing.american_to_prob and
    silently produces a plausible-but-wrong probability. Callers
    (sports/totals.py::_side_prices, runner.py's line-snapshot loop) already
    treat None as "no price" and drop the side/rung."""
    px = ln.get("americanPrice", ln.get("price"))
    if px is None:
        return None
    try:
        out_of_range = -pricing.MIN_ABS_AMERICAN_ODDS < px < pricing.MIN_ABS_AMERICAN_ODDS
    except TypeError:
        log.warning("line_american_price: non-numeric price %r", px)
        return None
    if out_of_range:
        log.warning("line_american_price: rejecting out-of-range American odds %r "
                    "(must be <=-%s or >=%s)", px, pricing.MIN_ABS_AMERICAN_ODDS,
                    pricing.MIN_ABS_AMERICAN_ODDS)
        return None
    return px

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
        for ev in events: _safe_ingest(st, ev, lk)
    return st


def _safe_ingest(st, ev, lk):
    """Ingest one event, isolating malformed events so a single bad event
    (schema drift, missing key) can't drop the rest of its batch."""
    try:
        _ingest(st, ev, lk)
    except Exception:
        log.warning("skipped malformed event in %s (eventId=%s)", lk,
                    ev.get("eventId") if isinstance(ev, dict) else "?")

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


# ---------------------------------------------------------------------------
# v2 per-league feed (current app path; supersedes snapshot+changes for soccer).
# Shape: {odds: {"lg21:pt1:pregame": [row, ...], ...}, teams: {tid: {name}}, ...}
# Each row = one event x betType: {eventId, eventStart, eventName ("Away @ Home"),
# betTypeId, sides: {"si0:tid<away>": {"ms<book>": line}, "si1:tid<home>": {...}}}.
# Line dicts carry points/americanPrice/isBlurred/modifiedOn + alternateLines.
# ---------------------------------------------------------------------------

_V2_BET_TYPES = {1, 2, 3, 4}   # moneyline/spread/total/soccer-extra; skip props/corners/etc.


def parse_v2(raw: dict, league_prefix: str) -> FeedState:
    st = FeedState()
    for s in raw.get("marketSources", []) or []:
        if isinstance(s, dict):
            st.books[s["id"]] = s.get("name")
    for tid, t in (raw.get("teams", {}) or {}).items():
        st.teams[tid] = t.get("name") if isinstance(t, dict) else None
    for lk, rows in (raw.get("odds", {}) or {}).items():
        # full-game pregame only: "lg21:pt1:pregame"
        if _league_prefix(lk) != league_prefix or not lk.endswith(":pt1:pregame"):
            continue
        for row in rows:
            _ingest_v2_row(st, row, league_prefix)
    return st


def _ingest_v2_row(st, row, league_prefix):
    try:
        bt = row.get("betTypeId")
        if bt not in _V2_BET_TYPES:
            return
        eid = row["eventId"]
        sides = row.get("sides") or {}
        if eid not in st.events:
            hid = aid = None
            for sk in sides:
                si, _, tid = sk.partition(":")
                tid = tid[3:] if tid.startswith("tid") else tid
                if si == "si1": hid = tid
                elif si == "si0": aid = tid
            st.events[eid] = EventMeta(eid, league_prefix, _dt(row["eventStart"]),
                                       hid, aid, st.teams.get(hid), st.teams.get(aid))
        for sk, books in sides.items():
            si = sk.split(":")[0][2:]                     # "si0:tid2058" -> "0"
            for ms_key, line in (books or {}).items():
                if not isinstance(line, dict) or line.get("isBlurred") is True:
                    continue
                ms = ms_key[2:]                           # "ms7" -> "7"
                st.lines[f"{eid}|{si}|{ms}|bt{bt}"] = line
    except Exception:
        log.warning("skipped malformed v2 row (eventId=%s)",
                    row.get("eventId") if isinstance(row, dict) else "?")


def fetch_v2(league_id: int, league_prefix: str, session=None) -> FeedState:
    url = config.UNABATED_V2_LEAGUE_URL.format(league_id=league_id)
    r = (session or requests).get(url, headers=_HEADERS,
                                  params={"v": str(int(datetime.datetime.now(datetime.timezone.utc).timestamp() * 1000))},
                                  timeout=30)
    r.raise_for_status()
    return parse_v2(r.json(), league_prefix)

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
            for ev in events: _safe_ingest(st, ev, lk)
