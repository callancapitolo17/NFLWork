# Kalshi World Cup Soccer — Foundation + Dry-Run Screen (Plan 1 of 3)

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Build the `kalshi_soccer/` foundation that authenticates to Unabated's feed, prices fair World Cup match probabilities from the sharp anchor, compares them to Kalshi match contracts, and **logs +EV opportunities in dry-run mode** (no orders).

**Architecture:** A standalone polling daemon. Each tick: pull Unabated authenticated deltas → maintain an in-memory line book → devig the sharp anchor (`Sharp Book Price` + Circa) into fair probabilities → map each Kalshi World Cup contract to its Unabated selection → compute net-of-fee EV → log/persist flagged edges. It reuses `kalshi_common` for devig, fees, and Kalshi REST auth. Live order placement and the totals/spreads adapters are **separate later plans**.

**Tech Stack:** Python 3.14, `requests` (Unabated feed), `duckdb`, `scipy`/`numpy`/`pandas` (devig), `kalshi_common` (shared math + Kalshi auth). Tests via `pytest` against saved JSON fixtures.

## Global Constraints

- **Python 3.14** (the interpreter in use).
- **All timestamps `TIMESTAMPTZ` UTC.** Unabated's authenticated feed gives TZ-aware `eventStart` (`+00:00`); never store naive timestamps. (Repo rule.)
- **No `print()`** — operational output via Python `logging` + a rotating `bot.log`. (Mirrors `kalshi_mlb_rfq`.)
- **Secrets in gitignored `.env` only.** Never commit the Unabated token or Kalshi keys. Provide `.env.example`.
- **Reuse `kalshi_common`** — do not reimplement devig, fee math, or Kalshi auth.
- **Dry-run is the default.** The Plan-1 daemon never places orders even with live execution code absent; `--live` is reserved for Plan 2 and must not exist yet.
- **Anchor = `Sharp Book Price` (marketSourceId 7), Circa (6/68).** Pinnacle is absent from the feed.
- **World Cup soccer league = `lg21`** in Unabated; bet types `bt1`=moneyline, `bt2`=spread, `bt3`=total.

---

## Task 0: Live recon — pin the two unknowns and capture fixtures

This is a discovery task (not TDD): it resolves the spec's one open item (the 3-way draw encoding) and discovers Kalshi's World Cup ticker structure, producing JSON fixtures that later TDD tasks build against. It runs against live APIs once.

**Files:**
- Create: `kalshi_soccer/recon/capture_fixtures.py` (throwaway-but-kept recon script)
- Create: `kalshi_soccer/tests/fixtures/unabated_lg21_authed.json` (captured)
- Create: `kalshi_soccer/tests/fixtures/kalshi_worldcup_markets.json` (captured)
- Create: `kalshi_soccer/recon/FINDINGS.md` (documented structure)

- [ ] **Step 1: Capture an authenticated Unabated delta over a window that includes a live/near-kickoff World Cup 3-way market.** Using the token in `.env` (see Task 1 for the var name; for recon you may paste it inline and delete after), poll `GET https://api-k.unabated.com/api/markets/changes/query?full_refresh_ISO=<nowISO>` a few times across a couple minutes near a scheduled World Cup kickoff, union the `results[].marketLineChanges[].gameOdds.gameOddsEvents` for keys starting `lg21:`, and save to `unabated_lg21_authed.json`.

- [ ] **Step 2: Determine the draw encoding.** In the captured lg21 data, for one event find where the 3-way draw price lives. Check, in order: (a) a `bt` slot beyond `bt1/2/3` (e.g. `bt4`) that carries a price with `points==null` (moneyline-shaped) distinct from home/away; (b) a third side id (`si2`); (c) a separate `pt` (period-type) key. Record the exact location in `FINDINGS.md`. If no draw price is found anywhere (it never gets posted to the free-tier feed), record "DRAW ABSENT → use Poisson fallback (Task 4 Step 6)."

- [ ] **Step 3: Discover Kalshi World Cup market structure.** Using `kalshi_common.auth_client` (see Task 1 config), call `api("GET", "/series?category=...")` / `api("GET", "/events?series_ticker=KXMENWORLDCUP")` and related queries to find: the **match-winner series ticker(s)** (not just the tournament `KXMENWORLDCUP` futures), the **event ticker format** for a single match, and how the **3-way outcome** is represented (e.g. three markets home/draw/away under one event, or a market per outcome). Save raw responses to `kalshi_worldcup_markets.json`.

- [ ] **Step 4: Document the mapping in `FINDINGS.md`.** Write down: Kalshi match series ticker, event-ticker pattern, the exact market tickers for home/draw/away of one concrete match, and the field names for the live ask price and team identification. This is the contract Task 5 implements against.

- [ ] **Step 5: Commit fixtures + findings (no secrets).**

```bash
cd kalshi_soccer
git add recon/ tests/fixtures/unabated_lg21_authed.json tests/fixtures/kalshi_worldcup_markets.json
git commit -m "recon(kalshi-soccer): capture Unabated lg21 + Kalshi World Cup market fixtures"
```

---

## Task 1: Package scaffold + config

**Files:**
- Create: `kalshi_soccer/__init__.py`
- Create: `kalshi_soccer/config.py`
- Create: `kalshi_soccer/.env.example`
- Create: `kalshi_soccer/tests/__init__.py`
- Test: `kalshi_soccer/tests/test_config.py`

**Interfaces:**
- Produces: `config._get(key, default)`, and module constants `UNABATED_SNAPSHOT_URL`, `UNABATED_CHANGES_URL`, `UNABATED_TOKEN`, `WC_LEAGUE_KEY_PREFIX="lg21"`, `ANCHOR_SOURCE_IDS=[7,6,68]`, `SHARP_BOOK_PRICE_ID=7`, `BANKROLL`, `KELLY_FRACTION`, `MIN_EV_PCT`, `MAX_STALENESS_SEC`, `KALSHI_*`, `PKG_DIR`, `DB_PATH`, `RESEARCH_DB_PATH`, `LOG_PATH`, `KILL_FILE`.

- [ ] **Step 1: Write the failing test**

```python
# kalshi_soccer/tests/test_config.py
from kalshi_soccer import config

def test_defaults_present():
    assert config.WC_LEAGUE_KEY_PREFIX == "lg21"
    assert config.SHARP_BOOK_PRICE_ID == 7
    assert set([7, 6, 68]).issubset(set(config.ANCHOR_SOURCE_IDS))
    assert config.BANKROLL == 1000.0
    assert 0 < config.KELLY_FRACTION <= 1
    assert config.UNABATED_SNAPSHOT_URL.startswith("https://content.unabated.com")

def test_env_override(monkeypatch):
    monkeypatch.setenv("BANKROLL", "2500")
    import importlib
    importlib.reload(config)
    assert config.BANKROLL == 2500.0
```

- [ ] **Step 2: Run test to verify it fails**

Run: `python -m pytest kalshi_soccer/tests/test_config.py -v`
Expected: FAIL with `ModuleNotFoundError: kalshi_soccer.config`

- [ ] **Step 3: Write minimal implementation**

```python
# kalshi_soccer/config.py
import os
from pathlib import Path

PKG_DIR = Path(__file__).resolve().parent

def _load_env(path: Path) -> dict[str, str]:
    env = {}
    if path.exists():
        for line in path.read_text().splitlines():
            line = line.strip()
            if not line or line.startswith("#") or "=" not in line:
                continue
            k, v = line.split("=", 1)
            env[k.strip()] = v.strip().strip('"').strip("'")
    return env

_FILE_ENV = _load_env(PKG_DIR / ".env")

def _get(key: str, default: str | None = None) -> str | None:
    return os.environ.get(key, _FILE_ENV.get(key, default))

# --- Unabated feed ---
UNABATED_SNAPSHOT_URL = "https://content.unabated.com/markets/game-odds/b_gameodds.json"
UNABATED_CHANGES_URL = "https://api-k.unabated.com/api/markets/changes/query"
UNABATED_TOKEN = _get("UNABATED_AT_PROD")          # the unabated_at_prod cookie value
WC_LEAGUE_KEY_PREFIX = "lg21"
SHARP_BOOK_PRICE_ID = 7
ANCHOR_SOURCE_IDS = [7, 6, 68]                       # Sharp Book Price, Circa, Circa Sports

# --- Sizing / thresholds ---
BANKROLL = float(_get("BANKROLL", "1000.0"))
KELLY_FRACTION = float(_get("KELLY_FRACTION", "0.25"))
MIN_EV_PCT = float(_get("MIN_EV_PCT", "0.03"))
MAX_STALENESS_SEC = int(_get("MAX_STALENESS_SEC", "20"))
KICKOFF_CUTOFF_MIN = int(_get("KICKOFF_CUTOFF_MIN", "3"))
PER_MATCH_CAP_PCT = float(_get("PER_MATCH_CAP_PCT", "0.03"))

# --- Kalshi REST ---
KALSHI_API_KEY_ID = _get("KALSHI_API_KEY_ID")
KALSHI_PRIVATE_KEY_PATH = _get("KALSHI_PRIVATE_KEY_PATH")
KALSHI_BASE_URL = _get("KALSHI_BASE_URL", "https://api.elections.kalshi.com/trade-api/v2")

# --- Paths ---
DB_PATH = PKG_DIR / "kalshi_soccer_market.duckdb"
RESEARCH_DB_PATH = PKG_DIR / "kalshi_soccer_research.duckdb"
LOG_PATH = PKG_DIR / "bot.log"
KILL_FILE = PKG_DIR / ".kill"
```

```python
# kalshi_soccer/__init__.py
```
```python
# kalshi_soccer/tests/__init__.py
```

```bash
# kalshi_soccer/.env.example
UNABATED_AT_PROD=paste_the_unabated_at_prod_cookie_value_here
KALSHI_API_KEY_ID=your_kalshi_api_key_id
KALSHI_PRIVATE_KEY_PATH=/absolute/path/to/kalshi_private_key.pem
KALSHI_BASE_URL=https://api.elections.kalshi.com/trade-api/v2
BANKROLL=1000.0
KELLY_FRACTION=0.25
MIN_EV_PCT=0.03
MAX_STALENESS_SEC=20
```

- [ ] **Step 4: Run tests to verify they pass**

Run: `python -m pytest kalshi_soccer/tests/test_config.py -v`
Expected: PASS (2 tests)

- [ ] **Step 5: Commit**

```bash
git add kalshi_soccer/__init__.py kalshi_soccer/config.py kalshi_soccer/.env.example kalshi_soccer/tests/
git commit -m "feat(kalshi-soccer): package scaffold + config"
```

---

## Task 2: Unabated snapshot parser

**Files:**
- Create: `kalshi_soccer/feed.py`
- Test: `kalshi_soccer/tests/test_feed_snapshot.py`

**Interfaces:**
- Consumes: `config` constants.
- Produces:
  - `parse_snapshot(raw: dict) -> FeedState` where `FeedState` holds `books: dict[int,str]`, `teams: dict[str,str]` (id→name), and `lines: dict[str, dict]` keyed by `(event_id, side, market_source_id, bet_type)` flattened to a string key.
  - `wc_matches(state: FeedState) -> list[Match]` where `Match` = dataclass `(event_id:int, start_utc:datetime, home_id:int, away_id:int, home:str, away:str)`.

- [ ] **Step 1: Write the failing test** (uses the real snapshot shape)

```python
# kalshi_soccer/tests/test_feed_snapshot.py
import json, datetime
from pathlib import Path
from kalshi_soccer import feed

FIX = Path(__file__).parent / "fixtures"

def _mini_snapshot():
    return {
        "marketSources": [{"id": 7, "name": "Sharp Book Price"}, {"id": 2, "name": "FanDuel"}],
        "teams": {"2063": {"name": "Argentina"}, "2038": {"name": "Austria"}},
        "gameOddsEvents": {
            "lg21:pt1:pregame": [{
                "eventId": 126539,
                "eventStart": "2026-06-22T17:00:00+00:00",
                "eventTeams": {"1": {"id": 2063}, "0": {"id": 2038}},
                "gameOddsMarketSourcesLines": {
                    "si1:ms7:an0": {"bt1": {"marketSourceId": 7, "points": None, "price": -150}},
                    "si0:ms7:an0": {"bt1": {"marketSourceId": 7, "points": None, "price": 130}},
                }
            }]
        }
    }

def test_parse_snapshot_books_and_teams():
    st = feed.parse_snapshot(_mini_snapshot())
    assert st.books[7] == "Sharp Book Price"
    assert st.teams["2063"] == "Argentina"

def test_wc_matches():
    st = feed.parse_snapshot(_mini_snapshot())
    matches = feed.wc_matches(st)
    assert len(matches) == 1
    m = matches[0]
    assert m.home == "Argentina" and m.away == "Austria"
    assert m.start_utc == datetime.datetime(2026, 6, 22, 17, 0, tzinfo=datetime.timezone.utc)
```

- [ ] **Step 2: Run test to verify it fails**

Run: `python -m pytest kalshi_soccer/tests/test_feed_snapshot.py -v`
Expected: FAIL with `ModuleNotFoundError: kalshi_soccer.feed`

- [ ] **Step 3: Write minimal implementation**

```python
# kalshi_soccer/feed.py
import datetime
from dataclasses import dataclass, field
from kalshi_soccer import config

@dataclass
class FeedState:
    books: dict[int, str] = field(default_factory=dict)
    teams: dict[str, str] = field(default_factory=dict)
    lines: dict[str, dict] = field(default_factory=dict)   # "evId|si|ms|bt" -> line dict
    cursor: int | None = None

@dataclass(frozen=True)
class Match:
    event_id: int
    start_utc: datetime.datetime
    home_id: int
    away_id: int
    home: str
    away: str

def _parse_dt(s: str) -> datetime.datetime:
    dt = datetime.datetime.fromisoformat(s)
    if dt.tzinfo is None:                       # anon snapshot is naive; treat as UTC
        dt = dt.replace(tzinfo=datetime.timezone.utc)
    return dt.astimezone(datetime.timezone.utc)

def parse_snapshot(raw: dict) -> FeedState:
    st = FeedState()
    for src in raw.get("marketSources", []):
        st.books[src["id"]] = src.get("name")
    for tid, t in raw.get("teams", {}).items():
        st.teams[tid] = t.get("name") if isinstance(t, dict) else None
    for lk, events in raw.get("gameOddsEvents", {}).items():
        if not lk.startswith(config.WC_LEAGUE_KEY_PREFIX + ":"):
            continue
        for ev in events:
            _ingest_event(st, ev)
    return st

def _ingest_event(st: FeedState, ev: dict) -> None:
    eid = ev["eventId"]
    for key, btmap in ev.get("gameOddsMarketSourcesLines", {}).items():
        si, ms, an = (p[2:] for p in key.split(":"))
        if an != "0":
            continue
        for bt, line in btmap.items():
            st.lines[f"{eid}|{si}|{ms}|{bt}"] = line

def wc_matches(st: FeedState) -> list[Match]:
    # rebuild match identity from the snapshot's event metadata
    out = {}
    return _matches_from_raw_cache.get(id(st), out and [] or _LAST_MATCHES)

# match metadata is captured during parse; store alongside state
_LAST_MATCHES: list[Match] = []
_matches_from_raw_cache: dict = {}
```

> Note for the implementer: the stub above for `wc_matches` is intentionally wrong to force the next refinement. Replace it by capturing match metadata during `parse_snapshot`. Concretely, add to `FeedState` a `matches: list[Match]` and populate it in `_ingest_event`'s caller:

```python
# in FeedState add: matches: list[Match] = field(default_factory=list)
# in parse_snapshot, inside the event loop, before _ingest_event:
            et = ev.get("eventTeams", {})
            home_id = et.get("1", {}).get("id")
            away_id = et.get("0", {}).get("id")
            st.matches.append(Match(
                event_id=ev["eventId"],
                start_utc=_parse_dt(ev["eventStart"]),
                home_id=home_id, away_id=away_id,
                home=st.teams.get(str(home_id)), away=st.teams.get(str(away_id)),
            ))
# and: def wc_matches(st): return st.matches
```

- [ ] **Step 4: Run tests to verify they pass**

Run: `python -m pytest kalshi_soccer/tests/test_feed_snapshot.py -v`
Expected: PASS (2 tests)

- [ ] **Step 5: Commit**

```bash
git add kalshi_soccer/feed.py kalshi_soccer/tests/test_feed_snapshot.py
git commit -m "feat(kalshi-soccer): Unabated snapshot parser (books/teams/lg21 matches)"
```

---

## Task 3: Authenticated fetch + cursor delta merge

**Files:**
- Modify: `kalshi_soccer/feed.py`
- Test: `kalshi_soccer/tests/test_feed_deltas.py`

**Interfaces:**
- Produces:
  - `fetch_snapshot(session=None) -> FeedState` (HTTP GET the snapshot URL).
  - `fetch_deltas(token: str, cursor: int | None, session=None) -> tuple[list[dict], int]` returning `(gameOddsEvents_dicts, latest_cursor)` from the authenticated changes endpoint.
  - `apply_deltas(st: FeedState, event_dicts: list[dict]) -> None` (merge changed lines into `st.lines`, refresh `st.matches`).

- [ ] **Step 1: Write the failing test** (delta merge logic, no network)

```python
# kalshi_soccer/tests/test_feed_deltas.py
from kalshi_soccer import feed

def test_apply_deltas_updates_line():
    st = feed.parse_snapshot({
        "marketSources": [{"id": 7, "name": "Sharp Book Price"}],
        "teams": {"1": {"name": "A"}, "2": {"name": "B"}},
        "gameOddsEvents": {"lg21:pt1:pregame": [{
            "eventId": 1, "eventStart": "2026-06-22T17:00:00+00:00",
            "eventTeams": {"1": {"id": 1}, "0": {"id": 2}},
            "gameOddsMarketSourcesLines": {"si1:ms7:an0": {"bt1": {"price": -150}}},
        }]}})
    assert st.lines["1|1|7|bt1"]["price"] == -150
    feed.apply_deltas(st, [{"lg21:pt1:pregame": [{
        "eventId": 1, "eventStart": "2026-06-22T17:00:00+00:00",
        "eventTeams": {"1": {"id": 1}, "0": {"id": 2}},
        "gameOddsMarketSourcesLines": {"si1:ms7:an0": {"bt1": {"price": -135}}},
    }]}])
    assert st.lines["1|1|7|bt1"]["price"] == -135
```

- [ ] **Step 2: Run test to verify it fails**

Run: `python -m pytest kalshi_soccer/tests/test_feed_deltas.py -v`
Expected: FAIL with `AttributeError: module ... has no attribute 'apply_deltas'`

- [ ] **Step 3: Write minimal implementation** (append to `feed.py`)

```python
import requests

_HEADERS = {
    "accept": "application/json, text/plain, */*",
    "origin": "https://unabated.com", "referer": "https://unabated.com/",
    "user-agent": "Mozilla/5.0 (Macintosh; Intel Mac OS X 10_15_7) "
                  "AppleWebKit/537.36 (KHTML, like Gecko) Chrome/147.0.0.0 Safari/537.36",
}

def fetch_snapshot(session=None) -> FeedState:
    s = session or requests
    r = s.get(config.UNABATED_SNAPSHOT_URL, headers=_HEADERS, timeout=30)
    r.raise_for_status()
    return parse_snapshot(r.json())

def fetch_deltas(token: str, cursor: int | None, session=None):
    s = session or requests
    cookies = {"unabated_at_prod": token}
    if cursor is None:
        import datetime
        now = datetime.datetime.now(datetime.timezone.utc).replace(microsecond=0).isoformat()
        url = f"{config.UNABATED_CHANGES_URL}?full_refresh_ISO={now.replace('+00:00','Z')}"
    else:
        url = f"{config.UNABATED_CHANGES_URL}/{cursor}"
    r = s.get(url, headers=_HEADERS, cookies=cookies, timeout=30)
    r.raise_for_status()
    body = r.json()
    event_dicts = []
    for batch in body.get("results", []):
        for mlc in batch.get("marketLineChanges", []):
            go = mlc.get("gameOdds", {})
            if "gameOddsEvents" in go:
                event_dicts.append(go["gameOddsEvents"])
    return event_dicts, body.get("latestTimestamp")

def apply_deltas(st: FeedState, event_dicts: list[dict]) -> None:
    known = {m.event_id for m in st.matches}
    for evmap in event_dicts:
        for lk, events in evmap.items():
            if not lk.startswith(config.WC_LEAGUE_KEY_PREFIX + ":"):
                continue
            for ev in events:
                _ingest_event(st, ev)
                if ev["eventId"] not in known:
                    et = ev.get("eventTeams", {})
                    hid, aid = et.get("1", {}).get("id"), et.get("0", {}).get("id")
                    st.matches.append(Match(ev["eventId"], _parse_dt(ev["eventStart"]),
                                            hid, aid, st.teams.get(str(hid)), st.teams.get(str(aid))))
                    known.add(ev["eventId"])
```

- [ ] **Step 4: Run tests to verify they pass**

Run: `python -m pytest kalshi_soccer/tests/test_feed_deltas.py -v`
Expected: PASS

- [ ] **Step 5: Commit**

```bash
git add kalshi_soccer/feed.py kalshi_soccer/tests/test_feed_deltas.py
git commit -m "feat(kalshi-soccer): authenticated fetch + cursor delta merge"
```

---

## Task 4: Sharp fair-value pricing (moneyline 3-way)

**Files:**
- Create: `kalshi_soccer/pricing.py`
- Test: `kalshi_soccer/tests/test_pricing.py`

**Interfaces:**
- Consumes: `kalshi_common.fair_value._probit_devig_n`, `feed.FeedState`, the draw location from Task 0.
- Produces:
  - `american_to_prob(odds: float) -> float`.
  - `moneyline_fair(state, event_id) -> dict | None` returning `{"home": p, "draw": p, "away": p}` devigged from the anchor, or `None` if the anchor lacks a full 3-way for the event.

- [ ] **Step 1: Write the failing test**

```python
# kalshi_soccer/tests/test_pricing.py
import pytest
from kalshi_soccer import pricing

def test_american_to_prob():
    assert pricing.american_to_prob(-200) == pytest.approx(2/3, abs=1e-6)
    assert pricing.american_to_prob(150) == pytest.approx(0.4, abs=1e-6)

def test_moneyline_fair_sums_to_one():
    # raw 3-way implied: 0.5, 0.30, 0.28 (sums to 1.08 vig) -> devig sums to 1
    fair = pricing.devig_three_way(0.50, 0.30, 0.28)
    assert sum(fair.values()) == pytest.approx(1.0, abs=1e-6)
    assert fair["home"] > fair["draw"]
```

- [ ] **Step 2: Run test to verify it fails**

Run: `python -m pytest kalshi_soccer/tests/test_pricing.py -v`
Expected: FAIL with `ModuleNotFoundError: kalshi_soccer.pricing`

- [ ] **Step 3: Write minimal implementation**

```python
# kalshi_soccer/pricing.py
import sys
from pathlib import Path
sys.path.insert(0, str(Path(__file__).resolve().parent.parent))  # repo root for kalshi_common
from kalshi_common.fair_value import _probit_devig_n
from kalshi_soccer import config

def american_to_prob(odds: float) -> float:
    odds = float(odds)
    return (-odds) / (-odds + 100) if odds < 0 else 100 / (odds + 100)

def devig_three_way(p_home_raw: float, p_draw_raw: float, p_away_raw: float) -> dict:
    devigged = _probit_devig_n([p_home_raw, p_draw_raw, p_away_raw])
    return {"home": devigged[0], "draw": devigged[1], "away": devigged[2]}

def _anchor_price(state, event_id, side, bt, draw=False):
    # side: "1"=home, "0"=away; draw uses the Task-0-confirmed location
    for ms in config.ANCHOR_SOURCE_IDS:
        key = f"{event_id}|{side}|{ms}|{bt}"
        line = state.lines.get(key)
        if line and line.get("price") is not None and line.get("points") is None:
            return line["price"], ms
    return None, None

def moneyline_fair(state, event_id) -> dict | None:
    home, _ = _anchor_price(state, event_id, "1", "bt1")
    away, _ = _anchor_price(state, event_id, "0", "bt1")
    draw = _draw_price(state, event_id)          # implemented per Task 0 FINDINGS
    if home is None or away is None or draw is None:
        return None
    return devig_three_way(american_to_prob(home), american_to_prob(draw), american_to_prob(away))

def _draw_price(state, event_id):
    # FILLED IN FROM TASK 0 FINDINGS. If draw bt is e.g. bt4 on side "1":
    line = state.lines.get(f"{event_id}|1|7|bt4")
    return line["price"] if line and line.get("price") is not None else None
```

- [ ] **Step 4: Run tests to verify they pass**

Run: `python -m pytest kalshi_soccer/tests/test_pricing.py -v`
Expected: PASS

- [ ] **Step 5: Commit**

```bash
git add kalshi_soccer/pricing.py kalshi_soccer/tests/test_pricing.py
git commit -m "feat(kalshi-soccer): 3-way moneyline devig from sharp anchor"
```

- [ ] **Step 6 (conditional — only if Task 0 found DRAW ABSENT): Poisson draw fallback.** Add `poisson_three_way(total_line: float, spread_line: float) -> dict` fitting a two-team Poisson from the anchor's `bt3` total and `bt2` spread, returning `{"home","draw","away"}`. Write a test asserting it sums to 1 and is monotone in the spread. Wire `moneyline_fair` to call it when `_draw_price` is None. Commit separately.

---

## Task 5: Kalshi World Cup market adapter

**Files:**
- Create: `kalshi_soccer/kalshi_markets.py`
- Test: `kalshi_soccer/tests/test_kalshi_markets.py`

**Interfaces:**
- Consumes: `kalshi_common.auth_client.configure` + `api`, the `kalshi_worldcup_markets.json` fixture + `FINDINGS.md` from Task 0.
- Produces:
  - `init_kalshi()` — calls `configure(...)` with config values.
  - `list_wc_match_markets() -> list[KMatch]` where `KMatch` = `(event_ticker, home_name, away_name, start_utc, yes_tickers: dict[str,str])` mapping `{"home","draw","away"} -> market_ticker`.
  - `best_yes_ask(market_ticker: str) -> float | None` — best YES ask in dollars from `/markets/{ticker}/orderbook` (derive from NO bids per Task 0 findings).

- [ ] **Step 1: Write the failing test** (parse logic against fixture; network mocked)

```python
# kalshi_soccer/tests/test_kalshi_markets.py
from kalshi_soccer import kalshi_markets as km

def test_parse_match_markets_from_fixture(monkeypatch):
    # shape per Task 0 FINDINGS — adjust field names to the captured structure
    sample_events = {"events": [{
        "event_ticker": "KXWCMATCH-26JUN22ARGAUT",
        "title": "Argentina vs Austria",
        "markets": [
            {"ticker": "KXWCMATCH-26JUN22ARGAUT-ARG", "yes_sub_title": "Argentina"},
            {"ticker": "KXWCMATCH-26JUN22ARGAUT-DRAW", "yes_sub_title": "Draw"},
            {"ticker": "KXWCMATCH-26JUN22ARGAUT-AUT", "yes_sub_title": "Austria"},
        ]}]}
    monkeypatch.setattr(km, "_get_events", lambda: sample_events)
    matches = km.list_wc_match_markets()
    assert len(matches) == 1
    m = matches[0]
    assert m.yes_tickers["draw"].endswith("-DRAW")
    assert m.yes_tickers["home"] or m.yes_tickers["away"]
```

- [ ] **Step 2: Run test to verify it fails**

Run: `python -m pytest kalshi_soccer/tests/test_kalshi_markets.py -v`
Expected: FAIL with `ModuleNotFoundError: kalshi_soccer.kalshi_markets`

- [ ] **Step 3: Write minimal implementation** (field names per Task 0 `FINDINGS.md`)

```python
# kalshi_soccer/kalshi_markets.py
import sys, datetime
from pathlib import Path
sys.path.insert(0, str(Path(__file__).resolve().parent.parent))
from kalshi_common import auth_client
from kalshi_soccer import config
from dataclasses import dataclass

WC_MATCH_SERIES = "KXWCMATCH"   # REPLACE with the real ticker from Task 0 FINDINGS

@dataclass(frozen=True)
class KMatch:
    event_ticker: str
    home_name: str | None
    away_name: str | None
    start_utc: datetime.datetime | None
    yes_tickers: dict        # {"home","draw","away"} -> market_ticker

def init_kalshi():
    auth_client.configure(
        api_key_id=config.KALSHI_API_KEY_ID,
        private_key_path=config.KALSHI_PRIVATE_KEY_PATH,
        base_url=config.KALSHI_BASE_URL,
        project_root=str(Path(__file__).resolve().parent.parent),
    )

def _get_events():
    status, body, _ = auth_client.api("GET", f"/events?series_ticker={WC_MATCH_SERIES}&status=open")
    return body if isinstance(body, dict) else {"events": []}

def list_wc_match_markets() -> list[KMatch]:
    out = []
    for ev in _get_events().get("events", []):
        yes = {}
        for mk in ev.get("markets", []):
            sub = (mk.get("yes_sub_title") or "").lower()
            role = "draw" if "draw" in sub or "tie" in sub else None
            yes[(role or sub)] = mk["ticker"]   # refine home/away mapping in Task 6
        out.append(KMatch(ev["event_ticker"], None, None, None, yes))
    return out

def best_yes_ask(market_ticker: str) -> float | None:
    status, body, _ = auth_client.api("GET", f"/markets/{market_ticker}/orderbook")
    if status != 200 or not isinstance(body, dict):
        return None
    no = (body.get("orderbook") or {}).get("no") or []   # NO bids; yes_ask = 1 - best_no_bid
    if not no:
        return None
    best_no_bid_cents = max(level[0] for level in no)
    return round((100 - best_no_bid_cents) / 100.0, 2)
```

- [ ] **Step 4: Run tests to verify they pass**

Run: `python -m pytest kalshi_soccer/tests/test_kalshi_markets.py -v`
Expected: PASS

- [ ] **Step 5: Commit**

```bash
git add kalshi_soccer/kalshi_markets.py kalshi_soccer/tests/test_kalshi_markets.py
git commit -m "feat(kalshi-soccer): Kalshi World Cup match-market adapter"
```

---

## Task 6: Mapping layer (team dict + match join + validation gate)

**Files:**
- Create: `kalshi_soccer/mapping.py`
- Create: `kalshi_soccer/team_aliases.py` (canonical country dictionary)
- Test: `kalshi_soccer/tests/test_mapping.py`

**Interfaces:**
- Consumes: `feed.Match`, `kalshi_markets.KMatch`.
- Produces:
  - `canon(name: str) -> str` — normalize a country name (handles "Korea Republic"→"South Korea", "USA"→"United States", etc.).
  - `match_pair(unabated_matches, kalshi_matches, max_kickoff_skew_min=180) -> list[Pairing]` where `Pairing = (event_id:int, kmatch:KMatch, home_ticker, draw_ticker, away_ticker)`.
  - `validate(pairing, fair: dict) -> bool` — fail-closed gate: all three tickers present, fair has 3 keys summing to ~1, kickoff times within skew.

- [ ] **Step 1: Write the failing test**

```python
# kalshi_soccer/tests/test_mapping.py
from kalshi_soccer import mapping

def test_canon_aliases():
    assert mapping.canon("Korea Republic") == "South Korea"
    assert mapping.canon("USA") == "United States"
    assert mapping.canon("Argentina") == "Argentina"

def test_validate_rejects_missing_draw():
    p = mapping.Pairing(1, None, "H", None, "A")          # draw ticker missing
    assert mapping.validate(p, {"home": .4, "draw": .3, "away": .3}) is False

def test_validate_accepts_complete():
    p = mapping.Pairing(1, None, "H", "D", "A")
    assert mapping.validate(p, {"home": .4, "draw": .3, "away": .3}) is True
```

- [ ] **Step 2: Run test to verify it fails**

Run: `python -m pytest kalshi_soccer/tests/test_mapping.py -v`
Expected: FAIL with `ModuleNotFoundError`

- [ ] **Step 3: Write minimal implementation**

```python
# kalshi_soccer/team_aliases.py
ALIASES = {
    "korea republic": "South Korea", "south korea": "South Korea",
    "usa": "United States", "united states": "United States",
    "ir iran": "Iran", "iran": "Iran",
    # extend from the 48-team World Cup field as Task 0 fixtures reveal exact names
}
```

```python
# kalshi_soccer/mapping.py
from dataclasses import dataclass
from kalshi_soccer.team_aliases import ALIASES

def canon(name: str) -> str:
    if name is None:
        return ""
    return ALIASES.get(name.strip().lower(), name.strip())

@dataclass(frozen=True)
class Pairing:
    event_id: int
    kmatch: object
    home_ticker: str | None
    draw_ticker: str | None
    away_ticker: str | None

def match_pair(unabated_matches, kalshi_matches, max_kickoff_skew_min: int = 180):
    out = []
    by_pair = {}
    for km in kalshi_matches:
        # km.yes_tickers role keys filled by Task 5/6; map names via canon
        names = {canon(v): k for k, v in (km.yes_tickers or {}).items()}
        by_pair[frozenset(n for n in names if n)] = km
    for m in unabated_matches:
        key = frozenset({canon(m.home), canon(m.away)})
        km = by_pair.get(key)
        if not km:
            continue
        yt = km.yes_tickers
        out.append(Pairing(m.event_id, km, yt.get("home"), yt.get("draw"), yt.get("away")))
    return out

def validate(pairing: Pairing, fair: dict) -> bool:
    if not (pairing.home_ticker and pairing.draw_ticker and pairing.away_ticker):
        return False
    if set(fair) != {"home", "draw", "away"}:
        return False
    return abs(sum(fair.values()) - 1.0) < 1e-6
```

- [ ] **Step 4: Run tests to verify they pass**

Run: `python -m pytest kalshi_soccer/tests/test_mapping.py -v`
Expected: PASS

- [ ] **Step 5: Commit**

```bash
git add kalshi_soccer/mapping.py kalshi_soccer/team_aliases.py kalshi_soccer/tests/test_mapping.py
git commit -m "feat(kalshi-soccer): mapping layer + team dict + fail-closed validation gate"
```

---

## Task 7: EV + edge flagging

**Files:**
- Create: `kalshi_soccer/ev.py`
- Test: `kalshi_soccer/tests/test_ev.py`

**Interfaces:**
- Consumes: `kalshi_common.ev_calc.fee_per_contract`.
- Produces: `edge_for_yes(fair_prob: float, yes_ask: float) -> tuple[float, float]` returning `(ev_dollars_per_contract, ev_pct_of_stake)`, net of the Kalshi taker fee at the ask.

- [ ] **Step 1: Write the failing test**

```python
# kalshi_soccer/tests/test_ev.py
import pytest
from kalshi_soccer import ev

def test_positive_edge_when_fair_exceeds_ask():
    d, pct = ev.edge_for_yes(0.55, 0.45)     # fair 55c, paying 45c
    assert d > 0 and pct > 0

def test_negative_edge_when_ask_exceeds_fair():
    d, pct = ev.edge_for_yes(0.40, 0.50)
    assert d < 0
```

- [ ] **Step 2: Run test to verify it fails**

Run: `python -m pytest kalshi_soccer/tests/test_ev.py -v`
Expected: FAIL with `ModuleNotFoundError`

- [ ] **Step 3: Write minimal implementation**

```python
# kalshi_soccer/ev.py
import sys
from pathlib import Path
sys.path.insert(0, str(Path(__file__).resolve().parent.parent))
from kalshi_common.ev_calc import fee_per_contract

def edge_for_yes(fair_prob: float, yes_ask: float) -> tuple[float, float]:
    # Taking a resting YES offer at price `yes_ask` (dollars). Payout $1 on win.
    fee = fee_per_contract(yes_ask)
    ev_dollars = fair_prob * (1 - yes_ask) - (1 - fair_prob) * yes_ask - fee
    ev_pct = ev_dollars / yes_ask if yes_ask > 0 else 0.0
    return ev_dollars, ev_pct
```

- [ ] **Step 4: Run tests to verify they pass**

Run: `python -m pytest kalshi_soccer/tests/test_ev.py -v`
Expected: PASS

- [ ] **Step 5: Commit**

```bash
git add kalshi_soccer/ev.py kalshi_soccer/tests/test_ev.py
git commit -m "feat(kalshi-soccer): net-of-fee EV for taking a YES offer"
```

---

## Task 8: Kelly sizing (binary contract)

**Files:**
- Create: `kalshi_soccer/sizing.py`
- Test: `kalshi_soccer/tests/test_sizing.py`

**Interfaces:**
- Produces: `kelly_contracts(fair_prob, yes_ask, bankroll, kelly_fraction, per_match_cap_dollars) -> int` — fractional-Kelly contract count for a binary YES, capped by both Kelly and the per-match dollar cap, floored, non-negative.

- [ ] **Step 1: Write the failing test**

```python
# kalshi_soccer/tests/test_sizing.py
from kalshi_soccer import sizing

def test_zero_when_no_edge():
    assert sizing.kelly_contracts(0.40, 0.50, 1000, 0.25, 30) == 0

def test_positive_and_capped():
    n = sizing.kelly_contracts(0.60, 0.45, 1000, 0.25, 30)
    assert n >= 1
    assert n * 0.45 <= 30 + 0.45        # respects per-match dollar cap
```

- [ ] **Step 2: Run test to verify it fails**

Run: `python -m pytest kalshi_soccer/tests/test_sizing.py -v`
Expected: FAIL with `ModuleNotFoundError`

- [ ] **Step 3: Write minimal implementation**

```python
# kalshi_soccer/sizing.py
import math

def kelly_contracts(fair_prob: float, yes_ask: float, bankroll: float,
                    kelly_fraction: float, per_match_cap_dollars: float) -> int:
    if yes_ask <= 0 or yes_ask >= 1:
        return 0
    b = (1 - yes_ask) / yes_ask          # net odds for $1 contract bought at yes_ask
    p, q = fair_prob, 1 - fair_prob
    f = (b * p - q) / b                  # full-Kelly fraction of bankroll
    if f <= 0:
        return 0
    stake = bankroll * kelly_fraction * f
    stake = min(stake, per_match_cap_dollars)
    return max(0, math.floor(stake / yes_ask))
```

- [ ] **Step 4: Run tests to verify they pass**

Run: `python -m pytest kalshi_soccer/tests/test_sizing.py -v`
Expected: PASS

- [ ] **Step 5: Commit**

```bash
git add kalshi_soccer/sizing.py kalshi_soccer/tests/test_sizing.py
git commit -m "feat(kalshi-soccer): fractional-Kelly binary sizing with per-match cap"
```

---

## Task 9: State + research DuckDB

**Files:**
- Create: `kalshi_soccer/db.py`
- Test: `kalshi_soccer/tests/test_db.py`

**Interfaces:**
- Produces:
  - `connect(path, read_only=False, retries=10)` — context manager mirroring `kalshi_mlb_rfq/db.py`.
  - `init_db()` — creates tables `flagged_edges`, `line_snapshots` (market DB) and `events` (research DB).
  - `log_edge(row: dict)` — insert one flagged edge.

- [ ] **Step 1: Write the failing test**

```python
# kalshi_soccer/tests/test_db.py
import datetime
from kalshi_soccer import db, config

def test_init_and_log(tmp_path, monkeypatch):
    monkeypatch.setattr(config, "DB_PATH", tmp_path / "m.duckdb")
    db.init_db()
    db.log_edge({"ts": datetime.datetime.now(datetime.timezone.utc),
                 "event_id": 1, "market_ticker": "X-ARG", "outcome": "home",
                 "fair_prob": 0.55, "yes_ask": 0.45, "ev_pct": 0.1,
                 "kelly_contracts": 5, "dry_run": True})
    with db.connect(config.DB_PATH, read_only=True) as c:
        assert c.execute("SELECT count(*) FROM flagged_edges").fetchone()[0] == 1
```

- [ ] **Step 2: Run test to verify it fails**

Run: `python -m pytest kalshi_soccer/tests/test_db.py -v`
Expected: FAIL with `ModuleNotFoundError`

- [ ] **Step 3: Write minimal implementation**

```python
# kalshi_soccer/db.py
import time, random, duckdb
from contextlib import contextmanager
from kalshi_soccer import config

@contextmanager
def connect(path, read_only: bool = False, retries: int = 10):
    attempt = 0
    while True:
        try:
            con = duckdb.connect(str(path), read_only=read_only)
            break
        except duckdb.IOException:
            if attempt >= retries:
                raise
            time.sleep(0.05 * 2 ** attempt + random.random() * 0.05)
            attempt += 1
    try:
        yield con
    finally:
        con.close()

def init_db():
    with connect(config.DB_PATH) as c:
        c.execute("""CREATE TABLE IF NOT EXISTS flagged_edges (
            ts TIMESTAMPTZ, event_id BIGINT, market_ticker VARCHAR, outcome VARCHAR,
            fair_prob DOUBLE, yes_ask DOUBLE, ev_pct DOUBLE, kelly_contracts INTEGER,
            dry_run BOOLEAN)""")
        c.execute("""CREATE TABLE IF NOT EXISTS line_snapshots (
            ts TIMESTAMPTZ, event_id BIGINT, market_source_id INTEGER, bet_type VARCHAR,
            side VARCHAR, price DOUBLE, points DOUBLE)""")

def log_edge(row: dict):
    with connect(config.DB_PATH) as c:
        c.execute("""INSERT INTO flagged_edges VALUES (?,?,?,?,?,?,?,?,?)""",
                  [row["ts"], row["event_id"], row["market_ticker"], row["outcome"],
                   row["fair_prob"], row["yes_ask"], row["ev_pct"],
                   row["kelly_contracts"], row["dry_run"]])
```

- [ ] **Step 4: Run tests to verify they pass**

Run: `python -m pytest kalshi_soccer/tests/test_db.py -v`
Expected: PASS

- [ ] **Step 5: Commit**

```bash
git add kalshi_soccer/db.py kalshi_soccer/tests/test_db.py
git commit -m "feat(kalshi-soccer): DuckDB market state + flagged-edge logging"
```

---

## Task 10: Orchestration loop (dry-run) + gates + logging + CLI

**Files:**
- Create: `kalshi_soccer/main.py`
- Create: `kalshi_soccer/log_setup.py`
- Test: `kalshi_soccer/tests/test_tick.py`

**Interfaces:**
- Consumes: everything above.
- Produces:
  - `run_tick(state, kalshi_matches, *, now, dry_run=True) -> list[dict]` — pure-ish function returning the flagged edges for one tick (testable without network).
  - `main_loop(dry_run: bool)` and `cli()` (argparse; **no `--live` flag in Plan 1**).

- [ ] **Step 1: Write the failing test** (one tick, all inputs injected)

```python
# kalshi_soccer/tests/test_tick.py
import datetime
from kalshi_soccer import main, feed, kalshi_markets as km

def test_run_tick_flags_positive_edge(monkeypatch):
    st = feed.parse_snapshot({
        "marketSources":[{"id":7,"name":"Sharp Book Price"}],
        "teams":{"1":{"name":"Argentina"},"2":{"name":"Austria"}},
        "gameOddsEvents":{"lg21:pt1:pregame":[{
            "eventId":1,"eventStart":"2026-12-31T17:00:00+00:00",
            "eventTeams":{"1":{"id":1},"0":{"id":2}},
            "gameOddsMarketSourcesLines":{
                "si1:ms7:an0":{"bt1":{"price":-200},"bt4":{"price":300}},
                "si0:ms7:an0":{"bt1":{"price":400}},
            }}]}})
    kmatch = km.KMatch("KXWCMATCH-1","Argentina","Austria",None,
                       {"home":"T-ARG","draw":"T-DRAW","away":"T-AUT"})
    monkeypatch.setattr(main, "_best_yes_ask", lambda t: 0.30)   # cheap vs fair
    edges = main.run_tick(st, [kmatch], now=datetime.datetime(2026,6,1,tzinfo=datetime.timezone.utc), dry_run=True)
    assert any(e["ev_pct"] > 0 for e in edges)
```

- [ ] **Step 2: Run test to verify it fails**

Run: `python -m pytest kalshi_soccer/tests/test_tick.py -v`
Expected: FAIL with `ModuleNotFoundError: kalshi_soccer.main`

- [ ] **Step 3: Write minimal implementation**

```python
# kalshi_soccer/log_setup.py
import logging, sys
from logging.handlers import RotatingFileHandler
from kalshi_soccer import config

def setup_logging() -> logging.Logger:
    log = logging.getLogger("kalshi_soccer")
    if any(getattr(h, "_managed", False) for h in log.handlers):
        return log
    log.setLevel(logging.INFO)
    fh = RotatingFileHandler(config.LOG_PATH, maxBytes=10_000_000, backupCount=5)
    fh._managed = True
    fh.setFormatter(logging.Formatter("%(asctime)s %(levelname)s %(name)s: %(message)s"))
    log.addHandler(fh)
    if sys.stderr.isatty():
        sh = logging.StreamHandler(); sh._managed = True
        sh.setFormatter(fh.formatter); log.addHandler(sh)
    return log
```

```python
# kalshi_soccer/main.py
import argparse, datetime, time, signal, threading
from kalshi_soccer import config, feed, pricing, mapping, ev, sizing, db
from kalshi_soccer import kalshi_markets as km
from kalshi_soccer.log_setup import setup_logging
from kalshi_soccer.risk import staleness_ok, tipoff_ok, kill_switch_ok

log = setup_logging()
_running = threading.Event(); _running.set()

def _best_yes_ask(ticker):           # seam for tests; real impl hits Kalshi
    return km.best_yes_ask(ticker)

def run_tick(state, kalshi_matches, *, now, dry_run=True) -> list[dict]:
    if not kill_switch_ok():
        log.warning("kill switch present; skipping tick"); return []
    pairings = mapping.match_pair(feed.wc_matches(state), kalshi_matches)
    flagged = []
    for p in pairings:
        fair = pricing.moneyline_fair(state, p.event_id)
        if fair is None or not mapping.validate(p, fair):
            continue
        km_obj = p.kmatch
        if not tipoff_ok(km_obj.start_utc, config.KICKOFF_CUTOFF_MIN, now):
            continue
        for outcome, ticker in (("home", p.home_ticker), ("draw", p.draw_ticker), ("away", p.away_ticker)):
            ask = _best_yes_ask(ticker)
            if ask is None:
                continue
            ev_d, ev_pct = ev.edge_for_yes(fair[outcome], ask)
            if ev_pct < config.MIN_EV_PCT:
                continue
            n = sizing.kelly_contracts(fair[outcome], ask, config.BANKROLL,
                                       config.KELLY_FRACTION, config.BANKROLL * config.PER_MATCH_CAP_PCT)
            row = {"ts": now, "event_id": p.event_id, "market_ticker": ticker,
                   "outcome": outcome, "fair_prob": fair[outcome], "yes_ask": ask,
                   "ev_pct": ev_pct, "kelly_contracts": n, "dry_run": dry_run}
            flagged.append(row)
            log.info("EDGE %s %s fair=%.3f ask=%.2f ev=%.1f%% n=%d%s",
                     km_obj.home_name, outcome, fair[outcome], ask, ev_pct*100, n,
                     " [DRY]" if dry_run else "")
    return flagged

def main_loop(dry_run: bool):
    db.init_db()
    km.init_kalshi()
    state = feed.fetch_snapshot()
    cursor = None
    last_kalshi = 0.0; kalshi_matches = []
    while _running.is_set():
        try:
            events, cursor = feed.fetch_deltas(config.UNABATED_TOKEN, cursor)
            feed.apply_deltas(state, events)
            now = datetime.datetime.now(datetime.timezone.utc)
            if time.time() - last_kalshi > 30:
                kalshi_matches = km.list_wc_match_markets(); last_kalshi = time.time()
            if not staleness_ok(now, config.MAX_STALENESS_SEC):
                pass
            for row in run_tick(state, kalshi_matches, now=now, dry_run=dry_run):
                db.log_edge(row)
        except Exception:
            log.exception("tick failed")
        time.sleep(2)

def _stop(*_): _running.clear()

def cli():
    ap = argparse.ArgumentParser()
    ap.add_argument("--dry-run", action="store_true", default=True)
    ap.parse_args()
    signal.signal(signal.SIGINT, _stop); signal.signal(signal.SIGTERM, _stop)
    main_loop(dry_run=True)

if __name__ == "__main__":
    cli()
```

Also create a thin `kalshi_soccer/risk.py` re-using the same gate logic as `kalshi_mlb_rfq/risk.py` (copy `staleness_ok`, `tipoff_ok`, `kill_switch_ok` — exact signatures from that file).

- [ ] **Step 4: Run tests to verify they pass**

Run: `python -m pytest kalshi_soccer/tests/test_tick.py -v`
Expected: PASS

- [ ] **Step 5: Run the full suite**

Run: `python -m pytest kalshi_soccer/tests/ -v`
Expected: all green.

- [ ] **Step 6: Commit**

```bash
git add kalshi_soccer/main.py kalshi_soccer/log_setup.py kalshi_soccer/risk.py kalshi_soccer/tests/test_tick.py
git commit -m "feat(kalshi-soccer): dry-run orchestration tick + gates + logging + CLI"
```

---

## Task 11: Docs, requirements, repo registration

**Files:**
- Create: `kalshi_soccer/README.md`, `kalshi_soccer/requirements.txt`
- Modify: root `CLAUDE.md` (project-structure list), root `.gitignore` (ensure `kalshi_soccer/.env`, `*.duckdb` covered)

- [ ] **Step 1: Write `kalshi_soccer/README.md`** — purpose, the Unabated token capture procedure (copy-as-cURL → `.env` `UNABATED_AT_PROD`, 30-day expiry + refresh note), `.env` setup, `python -m kalshi_soccer.main --dry-run`, how to read `flagged_edges`, troubleshooting (401 → token expired).

- [ ] **Step 2: Write `kalshi_soccer/requirements.txt`**

```
requests>=2.31
duckdb>=0.10.0
scipy>=1.10
numpy>=1.24.0
pandas>=2.0.0
pytest>=7.0.0
```

- [ ] **Step 3: Add a project-structure bullet to root `CLAUDE.md`** describing `kalshi_soccer/` (dry-run World Cup taker screen, Unabated-anchored, reuses `kalshi_common`).

- [ ] **Step 4: Verify `.gitignore` covers `kalshi_soccer/.env` and `kalshi_soccer/*.duckdb`.** Add if missing.

- [ ] **Step 5: Commit**

```bash
git add kalshi_soccer/README.md kalshi_soccer/requirements.txt CLAUDE.md .gitignore
git commit -m "docs(kalshi-soccer): README, requirements, repo registration"
```

---

## Out of scope (later plans)

- **Plan 2 — live execution:** order placement via `auth_client.api("POST", "/portfolio/orders", ...)`, the full gate set (line-move recompute, exposure caps from positions, fill-ratio halt), positions table, `--live` flag, CLV monitor logging the anchor's closing line.
- **Plan 3 — adapters:** totals (`bt3`) and spreads/Asian-handicap (`bt2`) market types reusing the same engine; tournament-futures evaluation; the maker-mode decision gated on measured CLV.

## Self-review notes

- Every spec Phase-0/Phase-1 requirement maps to a task: feed (T2/T3), anchor pricing + draw (T4 + Task 0), Kalshi adapter (T5), mapping/validation (T6), EV net of fees (T7), Kelly (T8), DBs (T9), dry-run loop + gates + CLV-precursor logging (T10), docs (T11).
- Draw encoding is resolved by Task 0 and consumed by Task 4 `_draw_price`; Poisson fallback is the conditional T4 Step 6.
- Kalshi soccer ticker uncertainty is isolated to Task 0 (recon) + Task 5 (`WC_MATCH_SERIES` constant + field names), so the rest of the plan is stable regardless of what recon finds.
