# Unabated-Anchored Kalshi Edge Engine — Core + Soccer Adapter (Plan 1 of N)

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Build the **sport-agnostic `unabated_edge/` core** + a **`sports/soccer.py` adapter** that authenticates to Unabated, devigs the sharp anchor into fair probabilities, compares to Kalshi World Cup contracts, **logs +EV opportunities in dry-run**, and **persists the full line history + research firehose** for CLV/validation. Adding a future sport = adding one adapter file.

**Architecture:** A core that owns the feed, devig math, EV/fees, mapping, storage, and the Kalshi venue. Each sport is a small adapter implementing a `SportAdapter` interface (its Unabated league key, fair-outcome computation + devig arity, team-name canonicalization, and venue-ticker mapping). A `runner` iterates registered adapters each tick. Live order placement, exposure-cap gates, and additional sports/market-types are **later plans**.

**Tech Stack:** Python 3.14, `requests`, `duckdb`, `scipy`/`numpy`/`pandas`, `kalshi_common` (devig, fees, Kalshi auth). Tests via `pytest` against saved JSON fixtures.

## Global Constraints

- **Python 3.14.**
- **All timestamps `TIMESTAMPTZ` UTC.** Authenticated feed `eventStart` is TZ-aware (`+00:00`); never store naive timestamps.
- **No `print()`** — `logging` + rotating `bot.log`.
- **Secrets in gitignored `.env` only** (Unabated token, Kalshi keys); ship `.env.example`.
- **Reuse `kalshi_common`** for devig, fee math, Kalshi auth — never reimplement.
- **Dry-run is the only mode in Plan 1.** No `--live` flag, no order placement.
- **Sport-agnostic core.** The core must not hardcode `lg21`, "draw", or any soccer concept — those live in `sports/soccer.py`. Every storage table carries a `sport` column.
- **Anchor = `Sharp Book Price` (7), Circa (6/68).** Pinnacle absent.

---

## Task 0: Live recon — pin the two unknowns, capture fixtures

Discovery task (not TDD). Resolves the draw encoding + Kalshi World Cup ticker structure; outputs fixtures for later TDD tasks.

**Files:**
- Create: `unabated_edge/recon/capture_fixtures.py`, `unabated_edge/recon/FINDINGS.md`
- Create: `unabated_edge/tests/fixtures/unabated_lg21_authed.json`, `unabated_edge/tests/fixtures/kalshi_worldcup_markets.json`

- [ ] **Step 1: Capture an authenticated Unabated delta** spanning a near-kickoff World Cup match: poll `GET https://api-k.unabated.com/api/markets/changes/query?full_refresh_ISO=<nowISO>` with cookie `unabated_at_prod=<token from .env>`, union `results[].marketLineChanges[].gameOdds.gameOddsEvents` for keys starting `lg21:`, save to `unabated_lg21_authed.json`.
- [ ] **Step 2: Locate the draw price.** In the captured lg21 data, for one event find the 3-way draw: check (a) a `bt` slot beyond `bt1/2/3` with `points==null` and a distinct moneyline-shaped price, (b) a third side `si2`, (c) a separate `pt` key. Record exact location in `FINDINGS.md`, or "DRAW ABSENT → Poisson fallback."
- [ ] **Step 3: Discover Kalshi World Cup market structure.** Via `kalshi_common.auth_client.api("GET", ...)` query `/events?series_ticker=KXMENWORLDCUP` and search `/series` for the **match-winner** series (not just tournament futures). Capture: match series ticker, event-ticker format, how home/draw/away map to market tickers, the live-ask field. Save raw to `kalshi_worldcup_markets.json`.
- [ ] **Step 4: Document in `FINDINGS.md`** — the exact draw field, Kalshi match series ticker, event/market ticker patterns, and the orderbook ask field. This is the contract Tasks 4/6 implement.
- [ ] **Step 5: Commit** (no secrets):
```bash
cd unabated_edge && git add recon/ tests/fixtures/
git commit -m "recon(unabated-edge): capture lg21 draw + Kalshi WC market fixtures"
```

---

## Task 1: Core scaffold + config

**Files:**
- Create: `unabated_edge/__init__.py`, `unabated_edge/config.py`, `unabated_edge/.env.example`, `unabated_edge/tests/__init__.py`
- Test: `unabated_edge/tests/test_config.py`

**Interfaces:**
- Produces: `config._get(key, default)` + constants `UNABATED_SNAPSHOT_URL`, `UNABATED_CHANGES_URL`, `UNABATED_TOKEN`, `ANCHOR_SOURCE_IDS=[7,6,68]`, `SHARP_BOOK_PRICE_ID=7`, `BANKROLL`, `KELLY_FRACTION`, `MIN_EV_PCT`, `MAX_STALENESS_SEC`, `KICKOFF_CUTOFF_MIN`, `PER_MATCH_CAP_PCT`, `KALSHI_*`, `PKG_DIR`, `MARKET_DB_PATH`, `RESEARCH_DB_PATH`, `LOG_PATH`, `KILL_FILE`. **No soccer constants here.**

- [ ] **Step 1: Write the failing test**
```python
# unabated_edge/tests/test_config.py
from unabated_edge import config

def test_core_defaults():
    assert config.SHARP_BOOK_PRICE_ID == 7
    assert {7, 6, 68}.issubset(set(config.ANCHOR_SOURCE_IDS))
    assert config.BANKROLL == 1000.0
    assert config.UNABATED_SNAPSHOT_URL.startswith("https://content.unabated.com")
    assert not hasattr(config, "WC_LEAGUE_KEY_PREFIX")   # sport specifics live in adapters
```

- [ ] **Step 2: Run test to verify it fails**

Run: `python -m pytest unabated_edge/tests/test_config.py -v`
Expected: FAIL `ModuleNotFoundError: unabated_edge.config`

- [ ] **Step 3: Write minimal implementation**
```python
# unabated_edge/config.py
import os
from pathlib import Path
PKG_DIR = Path(__file__).resolve().parent

def _load_env(path: Path) -> dict[str, str]:
    env = {}
    if path.exists():
        for line in path.read_text().splitlines():
            line = line.strip()
            if line and not line.startswith("#") and "=" in line:
                k, v = line.split("=", 1)
                env[k.strip()] = v.strip().strip('"').strip("'")
    return env

_FILE_ENV = _load_env(PKG_DIR / ".env")
def _get(key, default=None): return os.environ.get(key, _FILE_ENV.get(key, default))

UNABATED_SNAPSHOT_URL = "https://content.unabated.com/markets/game-odds/b_gameodds.json"
UNABATED_CHANGES_URL = "https://api-k.unabated.com/api/markets/changes/query"
UNABATED_TOKEN = _get("UNABATED_AT_PROD")
SHARP_BOOK_PRICE_ID = 7
ANCHOR_SOURCE_IDS = [7, 6, 68]

BANKROLL = float(_get("BANKROLL", "1000.0"))
KELLY_FRACTION = float(_get("KELLY_FRACTION", "0.25"))
MIN_EV_PCT = float(_get("MIN_EV_PCT", "0.03"))
MAX_STALENESS_SEC = int(_get("MAX_STALENESS_SEC", "20"))
KICKOFF_CUTOFF_MIN = int(_get("KICKOFF_CUTOFF_MIN", "3"))
PER_MATCH_CAP_PCT = float(_get("PER_MATCH_CAP_PCT", "0.03"))

KALSHI_API_KEY_ID = _get("KALSHI_API_KEY_ID")
KALSHI_PRIVATE_KEY_PATH = _get("KALSHI_PRIVATE_KEY_PATH")
KALSHI_BASE_URL = _get("KALSHI_BASE_URL", "https://api.elections.kalshi.com/trade-api/v2")

MARKET_DB_PATH = PKG_DIR / "unabated_edge_market.duckdb"
RESEARCH_DB_PATH = PKG_DIR / "unabated_edge_research.duckdb"
LOG_PATH = PKG_DIR / "bot.log"
KILL_FILE = PKG_DIR / ".kill"
```
```python
# unabated_edge/__init__.py
```
```python
# unabated_edge/tests/__init__.py
```
```bash
# unabated_edge/.env.example
UNABATED_AT_PROD=paste_unabated_at_prod_cookie_value
KALSHI_API_KEY_ID=your_kalshi_api_key_id
KALSHI_PRIVATE_KEY_PATH=/absolute/path/to/kalshi_private_key.pem
BANKROLL=1000.0
```

- [ ] **Step 4: Run tests to verify they pass** — Run: `python -m pytest unabated_edge/tests/test_config.py -v` → PASS
- [ ] **Step 5: Commit**
```bash
git add unabated_edge/__init__.py unabated_edge/config.py unabated_edge/.env.example unabated_edge/tests/
git commit -m "feat(unabated-edge): core scaffold + sport-agnostic config"
```

---

## Task 2: Pricing primitives (devig, arity-parametric)

Core math the adapters call. No sport specifics.

**Files:** Create `unabated_edge/pricing.py`; Test `unabated_edge/tests/test_pricing.py`

**Interfaces:**
- Produces: `american_to_prob(odds)->float`; `devig(probs: list[float])->list[float]` (wraps `kalshi_common.fair_value._probit_devig_n`, works for 2..n outcomes).

- [ ] **Step 1: Write the failing test**
```python
# unabated_edge/tests/test_pricing.py
import pytest
from unabated_edge import pricing

def test_american_to_prob():
    assert pricing.american_to_prob(-200) == pytest.approx(2/3, abs=1e-6)
    assert pricing.american_to_prob(150) == pytest.approx(0.4, abs=1e-6)

def test_devig_three_way_sums_to_one():
    out = pricing.devig([0.50, 0.30, 0.28])
    assert sum(out) == pytest.approx(1.0, abs=1e-6)
    assert out[0] > out[1]

def test_devig_two_way():
    out = pricing.devig([0.55, 0.52])
    assert sum(out) == pytest.approx(1.0, abs=1e-6)
```

- [ ] **Step 2: Run → FAIL** `ModuleNotFoundError`
- [ ] **Step 3: Write minimal implementation**
```python
# unabated_edge/pricing.py
import sys
from pathlib import Path
sys.path.insert(0, str(Path(__file__).resolve().parent.parent))  # repo root for kalshi_common
from kalshi_common.fair_value import _probit_devig_n

def american_to_prob(odds: float) -> float:
    odds = float(odds)
    return (-odds) / (-odds + 100) if odds < 0 else 100 / (odds + 100)

def devig(probs: list[float]) -> list[float]:
    return _probit_devig_n(list(probs))
```
- [ ] **Step 4: Run → PASS**
- [ ] **Step 5: Commit** `feat(unabated-edge): devig + american-odds pricing primitives`

---

## Task 3: Unabated feed — snapshot parse (adapter-driven league filter)

**Files:** Create `unabated_edge/feed.py`; Test `unabated_edge/tests/test_feed_snapshot.py`

**Interfaces:**
- Produces:
  - `FeedState` dataclass: `books: dict[int,str]`, `teams: dict[str,str]`, `lines: dict[str,dict]` (key `"evId|si|ms|bt"`), `events: dict[int,EventMeta]`, `cursor: int|None`.
  - `EventMeta`: `(event_id, league_key, start_utc, home_id, away_id, home, away)`.
  - `parse_snapshot(raw: dict, league_prefixes: set[str]) -> FeedState` — keeps only events whose `lgNN:` key prefix is in `league_prefixes` (supplied by registered adapters).
  - `events_for_league(st, league_prefix) -> list[EventMeta]`.

- [ ] **Step 1: Write the failing test**
```python
# unabated_edge/tests/test_feed_snapshot.py
import datetime
from unabated_edge import feed

def _snap():
    return {"marketSources":[{"id":7,"name":"Sharp Book Price"},{"id":2,"name":"FanDuel"}],
            "teams":{"2063":{"name":"Argentina"},"2038":{"name":"Austria"}},
            "gameOddsEvents":{"lg21:pt1:pregame":[{
                "eventId":126539,"eventStart":"2026-06-22T17:00:00+00:00",
                "eventTeams":{"1":{"id":2063},"0":{"id":2038}},
                "gameOddsMarketSourcesLines":{
                    "si1:ms7:an0":{"bt1":{"marketSourceId":7,"points":None,"price":-150}},
                    "si0:ms7:an0":{"bt1":{"marketSourceId":7,"points":None,"price":130}}}}],
              "lg5:pt1:pregame":[{"eventId":1,"eventStart":"2026-06-22T17:00:00+00:00",
                "eventTeams":{"1":{"id":2063},"0":{"id":2038}},"gameOddsMarketSourcesLines":{}}]}}

def test_parse_filters_to_registered_leagues():
    st = feed.parse_snapshot(_snap(), {"lg21"})
    assert st.books[7] == "Sharp Book Price"
    assert 126539 in st.events and 1 not in st.events      # lg5 filtered out
    e = st.events[126539]
    assert e.home == "Argentina" and e.away == "Austria"
    assert e.start_utc == datetime.datetime(2026,6,22,17,0,tzinfo=datetime.timezone.utc)
    assert st.lines["126539|1|7|bt1"]["price"] == -150
```

- [ ] **Step 2: Run → FAIL** `ModuleNotFoundError`
- [ ] **Step 3: Write minimal implementation**
```python
# unabated_edge/feed.py
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
```
- [ ] **Step 4: Run → PASS**
- [ ] **Step 5: Commit** `feat(unabated-edge): snapshot parser with adapter-driven league filter`

---

## Task 4: Authenticated fetch + cursor delta merge

**Files:** Modify `unabated_edge/feed.py`; Test `unabated_edge/tests/test_feed_deltas.py`

**Interfaces:**
- Produces: `fetch_snapshot(league_prefixes, session=None) -> FeedState`; `fetch_deltas(token, cursor, session=None) -> tuple[list[dict], int]`; `apply_deltas(st, event_dicts, league_prefixes) -> None`.

- [ ] **Step 1: Write the failing test**
```python
# unabated_edge/tests/test_feed_deltas.py
from unabated_edge import feed
def test_apply_deltas_updates_line():
    st = feed.parse_snapshot({"marketSources":[{"id":7,"name":"S"}],"teams":{"1":{"name":"A"},"2":{"name":"B"}},
        "gameOddsEvents":{"lg21:pt1:pregame":[{"eventId":1,"eventStart":"2026-06-22T17:00:00+00:00",
            "eventTeams":{"1":{"id":1},"0":{"id":2}},
            "gameOddsMarketSourcesLines":{"si1:ms7:an0":{"bt1":{"price":-150}}}}]}}, {"lg21"})
    feed.apply_deltas(st, [{"lg21:pt1:pregame":[{"eventId":1,"eventStart":"2026-06-22T17:00:00+00:00",
        "eventTeams":{"1":{"id":1},"0":{"id":2}},
        "gameOddsMarketSourcesLines":{"si1:ms7:an0":{"bt1":{"price":-135}}}}]}], {"lg21"})
    assert st.lines["1|1|7|bt1"]["price"] == -135
```
- [ ] **Step 2: Run → FAIL** `AttributeError: apply_deltas`
- [ ] **Step 3: Write minimal implementation** (append to `feed.py`)
```python
import requests
from unabated_edge import config

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
```
- [ ] **Step 4: Run → PASS**
- [ ] **Step 5: Commit** `feat(unabated-edge): authenticated fetch + cursor delta merge`

---

## Task 5: SportAdapter interface + soccer adapter

The architectural heart: the contract every sport implements, plus soccer.

**Files:** Create `unabated_edge/sports/__init__.py`, `unabated_edge/sports/base.py`, `unabated_edge/sports/soccer.py`, `unabated_edge/sports/registry.py`; Test `unabated_edge/tests/test_soccer_adapter.py`

**Interfaces:**
- Produces:
  - `base.SportAdapter` (ABC) with: attribute `sport: str`, `league_prefix: str`, `outcomes: tuple[str,...]`; methods `fair(state, event_meta) -> dict[str,float] | None`, `canon_team(name: str) -> str`, `kalshi_series() -> str`, `map_outcome_tickers(kalshi_event: dict) -> dict[str,str]`.
  - `soccer.Soccer()` implementing it: `sport="soccer"`, `league_prefix="lg21"`, `outcomes=("home","draw","away")`.
  - `registry.ADAPTERS: list[SportAdapter]` (= `[Soccer()]`), `registry.league_prefixes() -> set[str]`.

- [ ] **Step 1: Write the failing test**
```python
# unabated_edge/tests/test_soccer_adapter.py
import datetime
from unabated_edge import feed
from unabated_edge.sports.soccer import Soccer

def _state_with_event():
    return feed.parse_snapshot({"marketSources":[{"id":7,"name":"S"}],
        "teams":{"1":{"name":"Korea Republic"},"2":{"name":"Mexico"}},
        "gameOddsEvents":{"lg21:pt1:pregame":[{"eventId":9,"eventStart":"2026-06-22T17:00:00+00:00",
            "eventTeams":{"1":{"id":1},"0":{"id":2}},
            "gameOddsMarketSourcesLines":{
                "si1:ms7:an0":{"bt1":{"price":-150},"bt4":{"price":300}},
                "si0:ms7:an0":{"bt1":{"price":400}}}}]}}, {"lg21"})

def test_canon_alias():
    assert Soccer().canon_team("Korea Republic") == "South Korea"

def test_fair_three_way_sums_to_one():
    s = Soccer(); st = _state_with_event()
    fair = s.fair(st, st.events[9])
    assert set(fair) == {"home","draw","away"}
    assert abs(sum(fair.values()) - 1.0) < 1e-6
```

- [ ] **Step 2: Run → FAIL** `ModuleNotFoundError`
- [ ] **Step 3: Write minimal implementation**
```python
# unabated_edge/sports/__init__.py
```
```python
# unabated_edge/sports/base.py
from abc import ABC, abstractmethod

class SportAdapter(ABC):
    sport: str
    league_prefix: str
    outcomes: tuple

    @abstractmethod
    def fair(self, state, event_meta) -> dict | None: ...
    @abstractmethod
    def canon_team(self, name: str) -> str: ...
    @abstractmethod
    def kalshi_series(self) -> str: ...
    @abstractmethod
    def map_outcome_tickers(self, kalshi_event: dict) -> dict: ...
```
```python
# unabated_edge/sports/soccer.py
from unabated_edge.sports.base import SportAdapter
from unabated_edge import pricing, config

_ALIASES = {"korea republic":"South Korea","usa":"United States","ir iran":"Iran"}
WC_MATCH_SERIES = "KXWCMATCH"          # REPLACE from Task 0 FINDINGS

class Soccer(SportAdapter):
    sport = "soccer"
    league_prefix = "lg21"
    outcomes = ("home", "draw", "away")

    def canon_team(self, name: str) -> str:
        return _ALIASES.get((name or "").strip().lower(), (name or "").strip())

    def kalshi_series(self) -> str:
        return WC_MATCH_SERIES

    def _anchor_ml(self, st, eid, side):
        for ms in config.ANCHOR_SOURCE_IDS:
            ln = st.lines.get(f"{eid}|{side}|{ms}|bt1")
            if ln and ln.get("price") is not None and ln.get("points") is None:
                return ln["price"]
        return None

    def _draw(self, st, eid):
        # FROM TASK 0 FINDINGS — e.g. draw lives in bt4 on side "1":
        ln = st.lines.get(f"{eid}|1|{config.SHARP_BOOK_PRICE_ID}|bt4")
        return ln["price"] if ln and ln.get("price") is not None else None

    def fair(self, st, ev) -> dict | None:
        h, a, d = self._anchor_ml(st, ev.event_id, "1"), self._anchor_ml(st, ev.event_id, "0"), self._draw(st, ev.event_id)
        if h is None or a is None or d is None:
            return None
        ph, pd, pa = (pricing.american_to_prob(x) for x in (h, d, a))
        dh, dd, da = pricing.devig([ph, pd, pa])
        return {"home": dh, "draw": dd, "away": da}

    def map_outcome_tickers(self, kalshi_event: dict) -> dict:
        out = {}
        for mk in kalshi_event.get("markets", []):
            sub = (mk.get("yes_sub_title") or "").lower()
            if "draw" in sub or "tie" in sub:
                out["draw"] = mk["ticker"]
            else:
                out.setdefault("_named", []).append((sub, mk["ticker"]))
        return out   # home/away resolved against canon team names in mapping (Task 7)
```
```python
# unabated_edge/sports/registry.py
from unabated_edge.sports.soccer import Soccer
ADAPTERS = [Soccer()]
def league_prefixes(): return {a.league_prefix for a in ADAPTERS}
def by_league(prefix): return next((a for a in ADAPTERS if a.league_prefix == prefix), None)
```
- [ ] **Step 4: Run → PASS**
- [ ] **Step 5: Commit** `feat(unabated-edge): SportAdapter interface + soccer adapter`
- [ ] **Step 6 (conditional — Task 0 found DRAW ABSENT): Poisson draw fallback.** Add `_poisson_fair(total, spread)` to `soccer.py`, wire `fair()` to use it when `_draw` is None; test sums-to-1 + monotone in spread. Commit separately.

---

## Task 6: Kalshi venue adapter

**Files:** Create `unabated_edge/venues/__init__.py`, `unabated_edge/venues/kalshi.py`; Test `unabated_edge/tests/test_kalshi_venue.py`

**Interfaces:**
- Produces: `init()` (calls `auth_client.configure(...)`); `list_events(series_ticker) -> list[dict]`; `best_yes_ask(market_ticker) -> float | None` (from `/markets/{t}/orderbook`, `yes_ask = 1 - best_no_bid`).

- [ ] **Step 1: Write the failing test**
```python
# unabated_edge/tests/test_kalshi_venue.py
from unabated_edge.venues import kalshi
def test_best_yes_ask_from_no_bids(monkeypatch):
    monkeypatch.setattr(kalshi.auth_client, "api",
        lambda m,p,*a,**k: (200, {"orderbook":{"no":[[30,100],[28,50]]}}, {}))
    assert kalshi.best_yes_ask("X") == 0.70    # 1 - 0.30
```
- [ ] **Step 2: Run → FAIL** `ModuleNotFoundError`
- [ ] **Step 3: Write minimal implementation**
```python
# unabated_edge/venues/__init__.py
```
```python
# unabated_edge/venues/kalshi.py
import sys
from pathlib import Path
sys.path.insert(0, str(Path(__file__).resolve().parent.parent.parent))
from kalshi_common import auth_client
from unabated_edge import config

def init():
    auth_client.configure(api_key_id=config.KALSHI_API_KEY_ID,
        private_key_path=config.KALSHI_PRIVATE_KEY_PATH, base_url=config.KALSHI_BASE_URL,
        project_root=str(Path(__file__).resolve().parent.parent.parent))

def list_events(series_ticker: str) -> list[dict]:
    status, body, _ = auth_client.api("GET", f"/events?series_ticker={series_ticker}&status=open&with_nested_markets=true")
    return body.get("events", []) if isinstance(body, dict) else []

def best_yes_ask(market_ticker: str) -> float | None:
    status, body, _ = auth_client.api("GET", f"/markets/{market_ticker}/orderbook")
    if status != 200 or not isinstance(body, dict): return None
    no = (body.get("orderbook") or {}).get("no") or []
    if not no: return None
    return round((100 - max(l[0] for l in no)) / 100.0, 2)
```
- [ ] **Step 4: Run → PASS**
- [ ] **Step 5: Commit** `feat(unabated-edge): Kalshi venue adapter (events + best YES ask)`

---

## Task 7: Generic mapping + fail-closed validation

**Files:** Create `unabated_edge/mapping.py`; Test `unabated_edge/tests/test_mapping.py`

**Interfaces:**
- Produces: `Pairing(event_meta, outcome_tickers: dict[str,str])`; `pair_events(adapter, events_meta, kalshi_events) -> list[Pairing]` (joins on canon team pair, resolves home/away tickers via canon names); `validate(adapter, pairing, fair) -> bool`.

- [ ] **Step 1: Write the failing test**
```python
# unabated_edge/tests/test_mapping.py
from unabated_edge import mapping
from unabated_edge.sports.soccer import Soccer

def test_validate_rejects_incomplete():
    a = Soccer()
    p = mapping.Pairing(None, {"home":"H","draw":"D"})       # missing away
    assert mapping.validate(a, p, {"home":.4,"draw":.3,"away":.3}) is False

def test_validate_accepts_complete():
    a = Soccer()
    p = mapping.Pairing(None, {"home":"H","draw":"D","away":"A"})
    assert mapping.validate(a, p, {"home":.4,"draw":.3,"away":.3}) is True
```
- [ ] **Step 2: Run → FAIL** `ModuleNotFoundError`
- [ ] **Step 3: Write minimal implementation**
```python
# unabated_edge/mapping.py
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
        if not hit: continue
        _, named, extra = hit
        ot = dict(extra)
        ot["home"] = named.get(adapter.canon_team(ev.home))
        ot["away"] = named.get(adapter.canon_team(ev.away))
        out.append(Pairing(ev, ot))
    return out

def validate(adapter, pairing: Pairing, fair: dict) -> bool:
    if set(fair) != set(adapter.outcomes): return False
    if not all(pairing.outcome_tickers.get(o) for o in adapter.outcomes): return False
    return abs(sum(fair.values()) - 1.0) < 1e-6
```
- [ ] **Step 4: Run → PASS**
- [ ] **Step 5: Commit** `feat(unabated-edge): generic Unabated↔Kalshi mapping + validation gate`

---

## Task 8: EV (net of fee)

**Files:** Create `unabated_edge/ev.py`; Test `unabated_edge/tests/test_ev.py`

**Interfaces:** `edge_for_yes(fair_prob, yes_ask) -> tuple[float,float]` (ev$/contract, ev% of stake), net of Kalshi taker fee.

- [ ] **Step 1: Write the failing test**
```python
# unabated_edge/tests/test_ev.py
from unabated_edge import ev
def test_positive_and_negative():
    assert ev.edge_for_yes(0.55, 0.45)[0] > 0
    assert ev.edge_for_yes(0.40, 0.50)[0] < 0
```
- [ ] **Step 2: Run → FAIL** `ModuleNotFoundError`
- [ ] **Step 3: Write minimal implementation**
```python
# unabated_edge/ev.py
import sys
from pathlib import Path
sys.path.insert(0, str(Path(__file__).resolve().parent.parent))
from kalshi_common.ev_calc import fee_per_contract

def edge_for_yes(fair_prob: float, yes_ask: float) -> tuple[float, float]:
    fee = fee_per_contract(yes_ask)
    ev_d = fair_prob * (1 - yes_ask) - (1 - fair_prob) * yes_ask - fee
    return ev_d, (ev_d / yes_ask if yes_ask > 0 else 0.0)
```
- [ ] **Step 4: Run → PASS**  — [ ] **Step 5: Commit** `feat(unabated-edge): net-of-fee EV for taking a YES offer`

---

## Task 9: Kelly sizing

**Files:** Create `unabated_edge/sizing.py`; Test `unabated_edge/tests/test_sizing.py`

**Interfaces:** `kelly_contracts(fair_prob, yes_ask, bankroll, kelly_fraction, per_match_cap_dollars) -> int`.

- [ ] **Step 1: Write the failing test**
```python
# unabated_edge/tests/test_sizing.py
from unabated_edge import sizing
def test_zero_when_no_edge():
    assert sizing.kelly_contracts(0.40, 0.50, 1000, 0.25, 30) == 0
def test_positive_and_capped():
    n = sizing.kelly_contracts(0.60, 0.45, 1000, 0.25, 30)
    assert n >= 1 and n * 0.45 <= 30 + 0.45
```
- [ ] **Step 2: Run → FAIL** `ModuleNotFoundError`
- [ ] **Step 3: Write minimal implementation**
```python
# unabated_edge/sizing.py
import math
def kelly_contracts(fair_prob, yes_ask, bankroll, kelly_fraction, per_match_cap_dollars) -> int:
    if not (0 < yes_ask < 1): return 0
    b = (1 - yes_ask) / yes_ask
    f = (b * fair_prob - (1 - fair_prob)) / b
    if f <= 0: return 0
    stake = min(bankroll * kelly_fraction * f, per_match_cap_dollars)
    return max(0, math.floor(stake / yes_ask))
```
- [ ] **Step 4: Run → PASS** — [ ] **Step 5: Commit** `feat(unabated-edge): fractional-Kelly binary sizing`

---

## Task 10: Storage — line history + research firehose + flagged edges

**Files:** Create `unabated_edge/storage.py`; Test `unabated_edge/tests/test_storage.py`

**Interfaces:**
- Produces:
  - `connect(path, read_only=False, retries=10)` (context manager, IOException backoff).
  - `init()` — creates `line_snapshots`, `flagged_edges` (market DB) + `research_events` (research DB); all tables include `sport VARCHAR`.
  - `snapshot_lines(sport, rows: list[dict])` — bulk insert current anchor lines (ts, sport, event_id, market_source_id, bet_type, side, price, points).
  - `log_flagged(row: dict)` — one flagged edge (incl. `sport`).
  - `emit(event_type, sport, **payload)` + `flush()` — buffered research firehose (never raises).

- [ ] **Step 1: Write the failing test**
```python
# unabated_edge/tests/test_storage.py
import datetime
from unabated_edge import storage, config
def test_snapshot_flag_and_research(tmp_path, monkeypatch):
    monkeypatch.setattr(config, "MARKET_DB_PATH", tmp_path/"m.duckdb")
    monkeypatch.setattr(config, "RESEARCH_DB_PATH", tmp_path/"r.duckdb")
    storage.init()
    now = datetime.datetime.now(datetime.timezone.utc)
    storage.snapshot_lines("soccer", [{"ts":now,"event_id":9,"market_source_id":7,"bet_type":"bt1","side":"1","price":-150.0,"points":None}])
    storage.log_flagged({"ts":now,"sport":"soccer","event_id":9,"market_ticker":"X","outcome":"home",
                         "fair_prob":0.6,"yes_ask":0.5,"ev_pct":0.1,"kelly_contracts":3,"dry_run":True})
    storage.emit("candidate_priced", "soccer", event_id=9, ev_pct=0.1); storage.flush()
    with storage.connect(config.MARKET_DB_PATH, read_only=True) as c:
        assert c.execute("SELECT count(*) FROM line_snapshots").fetchone()[0] == 1
        assert c.execute("SELECT count(*) FROM flagged_edges").fetchone()[0] == 1
    with storage.connect(config.RESEARCH_DB_PATH, read_only=True) as c:
        assert c.execute("SELECT count(*) FROM research_events").fetchone()[0] == 1
```
- [ ] **Step 2: Run → FAIL** `ModuleNotFoundError`
- [ ] **Step 3: Write minimal implementation**
```python
# unabated_edge/storage.py
import time, random, json, datetime, duckdb
from contextlib import contextmanager
from unabated_edge import config

_BUFFER = []

@contextmanager
def connect(path, read_only=False, retries=10):
    attempt = 0
    while True:
        try: con = duckdb.connect(str(path), read_only=read_only); break
        except duckdb.IOException:
            if attempt >= retries: raise
            time.sleep(0.05 * 2**attempt + random.random()*0.05); attempt += 1
    try: yield con
    finally: con.close()

def init():
    with connect(config.MARKET_DB_PATH) as c:
        c.execute("""CREATE TABLE IF NOT EXISTS line_snapshots(
            ts TIMESTAMPTZ, sport VARCHAR, event_id BIGINT, market_source_id INTEGER,
            bet_type VARCHAR, side VARCHAR, price DOUBLE, points DOUBLE)""")
        c.execute("""CREATE TABLE IF NOT EXISTS flagged_edges(
            ts TIMESTAMPTZ, sport VARCHAR, event_id BIGINT, market_ticker VARCHAR, outcome VARCHAR,
            fair_prob DOUBLE, yes_ask DOUBLE, ev_pct DOUBLE, kelly_contracts INTEGER, dry_run BOOLEAN)""")
    with connect(config.RESEARCH_DB_PATH) as c:
        c.execute("""CREATE TABLE IF NOT EXISTS research_events(
            ts TIMESTAMPTZ, event_type VARCHAR, sport VARCHAR, payload JSON)""")

def snapshot_lines(sport, rows):
    if not rows: return
    with connect(config.MARKET_DB_PATH) as c:
        c.executemany("INSERT INTO line_snapshots VALUES (?,?,?,?,?,?,?,?)",
            [[r["ts"], sport, r["event_id"], r["market_source_id"], r["bet_type"], r["side"], r["price"], r["points"]] for r in rows])

def log_flagged(r):
    with connect(config.MARKET_DB_PATH) as c:
        c.execute("INSERT INTO flagged_edges VALUES (?,?,?,?,?,?,?,?,?,?)",
            [r["ts"], r["sport"], r["event_id"], r["market_ticker"], r["outcome"],
             r["fair_prob"], r["yes_ask"], r["ev_pct"], r["kelly_contracts"], r["dry_run"]])

def emit(event_type, sport, **payload):
    try: _BUFFER.append((datetime.datetime.now(datetime.timezone.utc), event_type, sport, json.dumps(payload)))
    except Exception: pass

def flush():
    if not _BUFFER: return
    batch, _BUFFER[:] = list(_BUFFER), []
    try:
        with connect(config.RESEARCH_DB_PATH) as c:
            c.executemany("INSERT INTO research_events VALUES (?,?,?,?)", batch)
    except Exception:
        pass
```
- [ ] **Step 4: Run → PASS** — [ ] **Step 5: Commit** `feat(unabated-edge): DuckDB line history + research firehose + flagged edges`

---

## Task 11: Runner — dry-run loop over adapters + gates + logging + CLI

**Files:** Create `unabated_edge/runner.py`, `unabated_edge/log_setup.py`, `unabated_edge/risk.py`; Test `unabated_edge/tests/test_runner_tick.py`

**Interfaces:**
- Produces: `run_tick(adapter, state, kalshi_events, *, now, dry_run=True, ask_fn) -> list[dict]` (pure-ish; flagged rows, persists via storage + research emit); `main_loop(dry_run)`, `cli()` (no `--live`).

- [ ] **Step 1: Write the failing test**
```python
# unabated_edge/tests/test_runner_tick.py
import datetime
from unabated_edge import feed, runner, config, storage
from unabated_edge.sports.soccer import Soccer

def test_tick_flags_positive_edge(tmp_path, monkeypatch):
    monkeypatch.setattr(config, "MARKET_DB_PATH", tmp_path/"m.duckdb")
    monkeypatch.setattr(config, "RESEARCH_DB_PATH", tmp_path/"r.duckdb")
    storage.init()
    st = feed.parse_snapshot({"marketSources":[{"id":7,"name":"S"}],
        "teams":{"1":{"name":"Argentina"},"2":{"name":"Austria"}},
        "gameOddsEvents":{"lg21:pt1:pregame":[{"eventId":1,"eventStart":"2026-12-31T17:00:00+00:00",
            "eventTeams":{"1":{"id":1},"0":{"id":2}},
            "gameOddsMarketSourcesLines":{"si1:ms7:an0":{"bt1":{"price":-200},"bt4":{"price":300}},
                                          "si0:ms7:an0":{"bt1":{"price":400}}}}]}}, {"lg21"})
    kev = {"event_ticker":"E","markets":[
        {"ticker":"T-ARG","yes_sub_title":"Argentina"},
        {"ticker":"T-DRAW","yes_sub_title":"Draw"},
        {"ticker":"T-AUT","yes_sub_title":"Austria"}]}
    rows = runner.run_tick(Soccer(), st, [kev],
        now=datetime.datetime(2026,6,1,tzinfo=datetime.timezone.utc), dry_run=True,
        ask_fn=lambda t: 0.30)
    assert any(r["ev_pct"] > 0 for r in rows)
```
- [ ] **Step 2: Run → FAIL** `ModuleNotFoundError`
- [ ] **Step 3: Write minimal implementation**
```python
# unabated_edge/log_setup.py
import logging, sys
from logging.handlers import RotatingFileHandler
from unabated_edge import config
def setup_logging():
    log = logging.getLogger("unabated_edge")
    if any(getattr(h,"_m",False) for h in log.handlers): return log
    log.setLevel(logging.INFO)
    fh = RotatingFileHandler(config.LOG_PATH, maxBytes=10_000_000, backupCount=5); fh._m=True
    fh.setFormatter(logging.Formatter("%(asctime)s %(levelname)s %(name)s: %(message)s")); log.addHandler(fh)
    if sys.stderr.isatty():
        sh = logging.StreamHandler(); sh._m=True; sh.setFormatter(fh.formatter); log.addHandler(sh)
    return log
```
```python
# unabated_edge/risk.py   (copy exact signatures from kalshi_mlb_rfq/risk.py)
import datetime
from unabated_edge import config
def staleness_ok(generated_at, max_age_sec):
    now = datetime.datetime.now(datetime.timezone.utc)
    if generated_at.tzinfo is None: generated_at = generated_at.replace(tzinfo=datetime.timezone.utc)
    return (now - generated_at).total_seconds() <= max_age_sec
def tipoff_ok(commence_time, cancel_min, now=None):
    if commence_time is None: return True
    now = now or datetime.datetime.now(datetime.timezone.utc)
    return (commence_time - now).total_seconds() > cancel_min * 60
def kill_switch_ok(): return not config.KILL_FILE.exists()
```
```python
# unabated_edge/runner.py
import argparse, datetime, time, signal, threading
from unabated_edge import config, feed, ev, sizing, mapping, storage
from unabated_edge.venues import kalshi
from unabated_edge.sports import registry
from unabated_edge.log_setup import setup_logging
from unabated_edge.risk import staleness_ok, tipoff_ok, kill_switch_ok

log = setup_logging()
_running = threading.Event(); _running.set()

def run_tick(adapter, state, kalshi_events, *, now, dry_run=True, ask_fn) -> list[dict]:
    if not kill_switch_ok(): return []
    flagged = []
    events = [e for e in state.events.values() if e.league_key == adapter.league_prefix]
    # persist line snapshot for these events (CLV backbone)
    snap = [{"ts": now, "event_id": k.split("|")[0], "market_source_id": int(k.split("|")[2]),
             "bet_type": k.split("|")[3], "side": k.split("|")[1],
             "price": v.get("price"), "points": v.get("points")}
            for k, v in state.lines.items() if int(k.split("|")[0]) in {e.event_id for e in events}]
    storage.snapshot_lines(adapter.sport, snap)
    for p in mapping.pair_events(adapter, events, kalshi_events):
        fair = adapter.fair(state, p.event_meta)
        if fair is None or not mapping.validate(adapter, p, fair): continue
        if not tipoff_ok(p.event_meta.start_utc, config.KICKOFF_CUTOFF_MIN, now): continue
        for outcome in adapter.outcomes:
            ticker = p.outcome_tickers[outcome]; ask = ask_fn(ticker)
            if ask is None: continue
            ev_d, ev_pct = ev.edge_for_yes(fair[outcome], ask)
            storage.emit("candidate_priced", adapter.sport, event_id=p.event_meta.event_id,
                         outcome=outcome, fair=fair[outcome], ask=ask, ev_pct=ev_pct)
            if ev_pct < config.MIN_EV_PCT: continue
            n = sizing.kelly_contracts(fair[outcome], ask, config.BANKROLL,
                                       config.KELLY_FRACTION, config.BANKROLL*config.PER_MATCH_CAP_PCT)
            row = {"ts": now, "sport": adapter.sport, "event_id": p.event_meta.event_id,
                   "market_ticker": ticker, "outcome": outcome, "fair_prob": fair[outcome],
                   "yes_ask": ask, "ev_pct": ev_pct, "kelly_contracts": n, "dry_run": dry_run}
            storage.log_flagged(row); flagged.append(row)
            log.info("EDGE %s %s fair=%.3f ask=%.2f ev=%.1f%% n=%d [DRY]",
                     adapter.sport, outcome, fair[outcome], ask, ev_pct*100, n)
    return flagged

def main_loop(dry_run: bool):
    storage.init(); kalshi.init()
    prefixes = registry.league_prefixes()
    state = feed.fetch_snapshot(prefixes); cursor = None
    last_k = 0.0; kalshi_events = {}
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

def _stop(*_): _running.clear()
def cli():
    argparse.ArgumentParser().parse_args()
    signal.signal(signal.SIGINT, _stop); signal.signal(signal.SIGTERM, _stop)
    main_loop(dry_run=True)
if __name__ == "__main__": cli()
```
- [ ] **Step 4: Run → PASS** — Run: `python -m pytest unabated_edge/tests/test_runner_tick.py -v`
- [ ] **Step 5: Run full suite** — `python -m pytest unabated_edge/tests/ -v` → all green
- [ ] **Step 6: Commit** `feat(unabated-edge): dry-run runner over adapters + gates + persistence + CLI`

---

## Task 12: Docs, requirements, repo registration

**Files:** Create `unabated_edge/README.md`, `unabated_edge/requirements.txt`; Modify root `CLAUDE.md`, `.gitignore`

- [ ] **Step 1: `README.md`** — architecture (core + adapters), the Unabated token capture (copy-as-cURL → `.env UNABATED_AT_PROD`, 30-day expiry/refresh), `python -m unabated_edge.runner`, reading `flagged_edges`/`research_events` for CLV, **how to add a sport (write one `sports/<x>.py` adapter + register it)**, troubleshooting (401 → token expired).
- [ ] **Step 2: `requirements.txt`**
```
requests>=2.31
duckdb>=0.10.0
scipy>=1.10
numpy>=1.24.0
pandas>=2.0.0
pytest>=7.0.0
```
- [ ] **Step 3: root `CLAUDE.md`** — add `unabated_edge/` project-structure bullet (sport-agnostic Unabated-anchored Kalshi edge engine; soccer adapter first; dry-run + persistence; reuses `kalshi_common`).
- [ ] **Step 4: `.gitignore`** — ensure `unabated_edge/.env` + `unabated_edge/*.duckdb` covered.
- [ ] **Step 5: Commit** `docs(unabated-edge): README, requirements, repo registration`

---

## Out of scope (later plans)

- **Plan 2 — live execution:** `execution.py` order placement via `auth_client.api("POST","/portfolio/orders",...)`, full gate set (line-move recompute, exposure caps from a positions table, fill-ratio halt), `--live` flag, CLV report off `line_snapshots`.
- **Plan 3 — soccer breadth:** totals (`bt3`) + spreads (`bt2`) outcomes on the soccer adapter; tournament futures.
- **Breadth plans — new sports:** one `sports/<sport>.py` adapter each (US sports gated on measured CLV, since Unabated hands everyone their EV).

## Self-review notes

- Spec coverage: feed (T3/T4), devig (T2), adapter interface + soccer fair/draw (T5 + Task 0), Kalshi venue (T6), mapping/validation (T7), EV net fees (T8), Kelly (T9), **storage incl. line history + research firehose** (T10), dry-run runner + gates + persistence (T11), docs + "how to add a sport" (T12). The generalize-now + storage amendments are both reflected.
- Sport-agnostic core enforced: `config` test asserts no `WC_LEAGUE_KEY_PREFIX`; `lg21`/"draw"/`bt4` appear only in `sports/soccer.py`.
- Unknowns quarantined to Task 0 + the two soccer-adapter constants (`WC_MATCH_SERIES`, `_draw` location).
- Type consistency: `Pairing.outcome_tickers`, `adapter.outcomes`, `ask_fn`, `storage.emit/flush/log_flagged/snapshot_lines` names match across T5–T11.
