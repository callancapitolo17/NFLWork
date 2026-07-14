# WC Totals Maker Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Market maker for Kalshi `KXWCTOTAL` resting limit orders around the devigged Unabated anchor ladder, with an exact goal-grid risk ledger, live/shadow gateway, and full measurement logging.

**Architecture:** New `unabated_edge/maker/` module (ledger → store → gateway → state → engine) invoked from `runner.run_tick` on the existing 5s tick, per approved spec `docs/superpowers/specs/2026-07-10-wc-totals-maker-design.md`. Pure math isolated in `ledger.py`; all Kalshi REST behind `QuoteGateway`.

**Tech Stack:** Python 3.14, DuckDB, pytest, `kalshi_common.auth_client` / `ev_calc`.

## Global Constraints

- Branch: `worktree-wc-maker-depth-capture` (worktree at `/Users/callancapitolo/NFLWork/.claude/worktrees/wc-maker-depth-capture`). All commands run from that worktree root. NEVER merge to main, never push, never place a real order — the user gate happens after this plan completes.
- The live runner is currently executing from this worktree (pid varies). Editing files does not affect the running process; do NOT restart it, do NOT delete the worktree or its `.duckdb` files.
- Existing 55 tests must stay green after every task: `python3 -m pytest unabated_edge/tests/ -q`.
- All config via `unabated_edge/config.py` `_get()` pattern (env > .env > default). Percent caps are fractions of `BANKROLL`.
- `MAKER_MODE` defaults to `"off"` (spec deviation, safety: the capture-only runner must start without quoting; live is enabled in `.env` at the user gate with `MAKER_LIVE_ACK=1`).
- Staleness gate = feed-fetch watchdog (last successful tick age), NOT `modifiedOn` age — a quiet market's old `modifiedOn` is legitimate; a failing feed fetch is not (spec §3 deviation, documented).
- Ledger hypothetical includes ALL resting quotes of the match treated as filled (conservative spec tightening: simultaneous fills across rungs can never breach the match cap).
- Money convention: prices in tables/logs are DOLLARS (0.42); gateway payloads use integer CENTS. Contract counts are floats (Kalshi fp).
- Adverse-selection marks (+60s/+5min per fill) are NOT live columns: they are computed offline by joining `maker_fills.ts` to `line_snapshots` (the anchor fair is already captured every 5s). Same for spread-captured attribution — derivable from `maker_fills` + `maker_pnl`. Spec §6's mark columns are intentionally dropped from the live schema.

---

### Task 1: Risk ledger (`maker/ledger.py`)

**Files:**
- Create: `unabated_edge/maker/__init__.py` (empty)
- Create: `unabated_edge/maker/ledger.py`
- Test: `unabated_edge/tests/test_maker_ledger.py`

**Interfaces:**
- Produces: `pnl_grid(fills) -> list[float]` (len 11, g=0..10), `worst_case(fills) -> float` (0.0 when empty), `max_contracts(fills, line, side, price, budget) -> int`. A fill is the tuple `(line: float, side: "yes"|"no", contracts: float, price: float dollars)`; `side` is relative to the Over-market (YES = over, NO = under). `G_MAX = 10`.

- [ ] **Step 1: Write the failing tests**

```python
# unabated_edge/tests/test_maker_ledger.py
from unabated_edge.maker import ledger


def test_empty_book_is_flat():
    assert ledger.worst_case([]) == 0.0
    assert ledger.pnl_grid([]) == [0.0] * (ledger.G_MAX + 1)


def test_single_yes_fill_worst_case():
    # 40x YES Over-2.5 @ 0.42: loses 16.8 when g<=2, wins 23.2 when g>=3
    fills = [(2.5, "yes", 40, 0.42)]
    grid = ledger.pnl_grid(fills)
    assert round(grid[0], 2) == -16.8 and round(grid[2], 2) == -16.8
    assert round(grid[3], 2) == 23.2
    assert round(ledger.worst_case(fills), 2) == -16.8


def test_offset_under25_over15_worst_and_middle():
    # From the design conversation: 20x NO Over-2.5 @0.55 + 20x YES Over-1.5 @0.70
    fills = [(2.5, "no", 20, 0.55), (1.5, "yes", 20, 0.70)]
    grid = ledger.pnl_grid(fills)
    assert round(grid[0], 2) == -5.0          # low-goal worlds
    assert round(grid[2], 2) == 15.0          # the middle: both win at exactly 2
    assert round(grid[3], 2) == -5.0          # high-goal worlds
    assert round(ledger.worst_case(fills), 2) == -5.0


def test_grid_constant_above_top_rung():
    # No Kalshi rung above 5.5 -> pnl must be identical for g=6..10
    fills = [(5.5, "yes", 10, 0.05), (2.5, "no", 7, 0.6)]
    grid = ledger.pnl_grid(fills)
    assert len(set(round(v, 9) for v in grid[6:])) == 1


def test_max_contracts_shrinks_to_budget():
    # existing 40x YES Over-2.5 @0.42 (worst -16.8); candidate YES Over-1.5 @0.70,
    # budget 30 -> allowed n = floor((30-16.8)/0.70) = 18
    fills = [(2.5, "yes", 40, 0.42)]
    n = ledger.max_contracts(fills, 1.5, "yes", 0.70, 30.0)
    assert n == 18
    assert ledger.worst_case(fills + [(1.5, "yes", 18, 0.70)]) >= -30.0
    assert ledger.worst_case(fills + [(1.5, "yes", 19, 0.70)]) < -30.0


def test_max_contracts_opposite_direction_gets_room():
    # candidate NO Over-2.5 offsets the existing YES: plenty of room
    fills = [(2.5, "yes", 40, 0.42)]
    n = ledger.max_contracts(fills, 2.5, "no", 0.58, 30.0)
    assert n > 40


def test_max_contracts_empty_book():
    assert ledger.max_contracts([], 2.5, "yes", 0.42, 400.0) == int(400 / 0.42)


def test_max_contracts_zero_when_already_breached_or_no_budget():
    fills = [(2.5, "yes", 100, 0.42)]          # worst -42
    assert ledger.max_contracts(fills, 1.5, "yes", 0.70, 30.0) == 0
    assert ledger.max_contracts([], 2.5, "yes", 0.42, 0.0) == 0
```

- [ ] **Step 2: Run tests to verify they fail**

Run: `python3 -m pytest unabated_edge/tests/test_maker_ledger.py -q`
Expected: FAIL / error `ModuleNotFoundError: unabated_edge.maker`

- [ ] **Step 3: Implement**

```python
# unabated_edge/maker/__init__.py   (empty file)
```

```python
# unabated_edge/maker/ledger.py
"""Per-match risk ledger: exact worst-case P&L over the total-goals grid.

Every rung of a match settles on one integer (regulation-time total goals),
so the match portfolio's P&L is a function of that single number. Evaluating
it exactly on g=0..G_MAX replaces per-market caps and prices offsets, hedges,
and middles correctly with no special cases.

A fill is (line, side, contracts, price): side is relative to the Over
market — "yes" wins when g > line, "no" wins when g <= line. Prices in
dollars; P&L per contract is (1 - price) on a win, -price on a loss.
"""

G_MAX = 10   # "10 or more" bucket; P&L is constant above the top rung (tested)


def pnl(fills, g: int) -> float:
    total = 0.0
    for line, side, contracts, price in fills:
        win = (g > line) if side == "yes" else (g <= line)
        total += contracts * ((1.0 - price) if win else -price)
    return total


def pnl_grid(fills) -> list[float]:
    return [pnl(fills, g) for g in range(G_MAX + 1)]


def worst_case(fills) -> float:
    if not fills:
        return 0.0
    return min(pnl_grid(fills))


def max_contracts(fills, line: float, side: str, price: float, budget: float) -> int:
    """Largest integer n >= 0 such that adding (line, side, n, price) keeps
    worst_case >= -budget.

    Closed form: each goal outcome g imposes base_g + n*unit_g >= -budget and
    only outcomes where the candidate loses (unit_g < 0) bind. If the book is
    already beyond budget we never add exposure (fail closed)."""
    if budget <= 0:
        return 0
    base = pnl_grid(fills) if fills else [0.0] * (G_MAX + 1)
    if min(base) < -budget - 1e-9:
        return 0
    bound = None
    for g in range(G_MAX + 1):
        win = (g > line) if side == "yes" else (g <= line)
        unit = (1.0 - price) if win else -price
        if unit < 0:
            allowed = (budget + base[g]) / (-unit)
            bound = allowed if bound is None else min(bound, allowed)
    return int(bound + 1e-9) if bound is not None else 0
```

- [ ] **Step 4: Run tests to verify they pass**

Run: `python3 -m pytest unabated_edge/tests/test_maker_ledger.py -q`
Expected: `8 passed`

- [ ] **Step 5: Commit**

```bash
git add unabated_edge/maker/ unabated_edge/tests/test_maker_ledger.py
git commit -m "feat(unabated-edge): maker risk ledger — exact worst-case over the goal grid"
```

---

### Task 2: Config + maker DB (`config.py`, `maker/store.py`)

**Files:**
- Modify: `unabated_edge/config.py` (append after `PER_MATCH_CAP_PCT` line)
- Create: `unabated_edge/maker/store.py`
- Test: `unabated_edge/tests/test_maker_store.py`

**Interfaces:**
- Produces config names (module attrs on `unabated_edge.config`): `MAKER_MODE ("off")`, `MAKER_LIVE_ACK (None)`, `ROI_MARGIN (0.03)`, `PICKOFF_BUFFER_CENTS (1)`, `MAX_MARGIN_CENTS (5)`, `ALT_MARGIN_MULT (1.5)`, `ALT_SIZE_MULT (0.5)`, `ALT_OVERROUND_MIN (1.01)`, `ALT_OVERROUND_MAX (1.15)`, `QUOTE_PULL_MIN (3.0)`, `MAX_QUOTE_PCT (0.30)`, `MATCH_CAP_PCT (0.40)`, `GLOBAL_CAP_PCT (0.75)`, `DAILY_LOSS_HALT_PCT (0.40)`, `FILL_BURST_N (3)`, `COOLOFF_MIN (10.0)`, `MAKER_DB_PATH`.
- Produces store functions: `init()`, `log_quote(ts, sport, event_id, ticker, side, action, price, size, fair, margin, alt, reason, order_id)`, `log_fill(ts, sport, event_id, order_id, ticker, side, price, contracts, fee, worst_after, trade_id)`, `log_ledger(ts, sport, event_id, worst_case, grid, quotes_live)`, `log_settlement(ts, sport, ticker, pnl)`.

- [ ] **Step 1: Write the failing test**

```python
# unabated_edge/tests/test_maker_store.py
import datetime
from unabated_edge import config
from unabated_edge.maker import store
from unabated_edge.storage import connect


def test_maker_config_defaults():
    assert config.MAKER_MODE == "off"
    assert config.MATCH_CAP_PCT == 0.40
    assert config.GLOBAL_CAP_PCT == 0.75
    assert config.MAX_QUOTE_PCT == 0.30
    assert config.QUOTE_PULL_MIN == 3.0


def test_store_roundtrip(tmp_path, monkeypatch):
    monkeypatch.setattr(config, "MAKER_DB_PATH", tmp_path / "mk.duckdb")
    store.init()
    now = datetime.datetime.now(datetime.timezone.utc)
    store.log_quote(now, "soccer", 1, "T-O25", "yes", "rest", 0.40, 30, 0.42, 0.02, False, "quote", "oid1")
    store.log_fill(now, "soccer", 1, "oid1", "T-O25", "yes", 0.40, 5.0, 0.01, -2.0, "tr1")
    store.log_fill(now, "soccer", 1, "oid1", "T-O25", "yes", 0.40, 5.0, 0.01, -2.0, "tr1")  # dup trade_id
    store.log_ledger(now, "soccer", 1, -2.0, [-2.0] * 11, 4)
    store.log_settlement(now, "soccer", "T-O25", 3.5)
    with connect(config.MAKER_DB_PATH, read_only=True) as c:
        assert c.execute("SELECT count(*) FROM maker_quotes").fetchone()[0] == 1
        assert c.execute("SELECT count(*) FROM maker_fills").fetchone()[0] == 1   # deduped
        assert c.execute("SELECT count(*) FROM ledger_snapshots").fetchone()[0] == 1
        assert c.execute("SELECT settled_pnl FROM maker_pnl").fetchone()[0] == 3.5
```

- [ ] **Step 2: Run test to verify it fails**

Run: `python3 -m pytest unabated_edge/tests/test_maker_store.py -q`
Expected: FAIL (`MAKER_MODE` missing / no module `store`)

- [ ] **Step 3: Implement**

Append to `unabated_edge/config.py` (after the `PER_MATCH_CAP_PCT` line):

```python
# ---- maker (unabated_edge/maker/) — spec docs/superpowers/specs/2026-07-10-wc-totals-maker-design.md ----
MAKER_MODE = _get("MAKER_MODE", "off")            # off | shadow | live
MAKER_LIVE_ACK = _get("MAKER_LIVE_ACK")           # dead-man switch: must be "1" for live
ROI_MARGIN = float(_get("ROI_MARGIN", "0.03"))
PICKOFF_BUFFER_CENTS = int(_get("PICKOFF_BUFFER_CENTS", "1"))
MAX_MARGIN_CENTS = int(_get("MAX_MARGIN_CENTS", "5"))
ALT_MARGIN_MULT = float(_get("ALT_MARGIN_MULT", "1.5"))
ALT_SIZE_MULT = float(_get("ALT_SIZE_MULT", "0.5"))
ALT_OVERROUND_MIN = float(_get("ALT_OVERROUND_MIN", "1.01"))
ALT_OVERROUND_MAX = float(_get("ALT_OVERROUND_MAX", "1.15"))
QUOTE_PULL_MIN = float(_get("QUOTE_PULL_MIN", "3"))
MAX_QUOTE_PCT = float(_get("MAX_QUOTE_PCT", "0.30"))
MATCH_CAP_PCT = float(_get("MATCH_CAP_PCT", "0.40"))
GLOBAL_CAP_PCT = float(_get("GLOBAL_CAP_PCT", "0.75"))
DAILY_LOSS_HALT_PCT = float(_get("DAILY_LOSS_HALT_PCT", "0.40"))
FILL_BURST_N = int(_get("FILL_BURST_N", "3"))
COOLOFF_MIN = float(_get("COOLOFF_MIN", "10"))
MAKER_DB_PATH = PKG_DIR / "unabated_edge_maker.duckdb"
```

```python
# unabated_edge/maker/store.py
"""Maker sibling DB (unabated_edge_maker.duckdb): quote decisions, fills,
ledger snapshots, settlements. Separate file so maker state can be reset
without touching the capture history (three-DB discipline)."""
import json
import logging

from unabated_edge import config
from unabated_edge.storage import connect

log = logging.getLogger("unabated_edge")


def init():
    with connect(config.MAKER_DB_PATH) as c:
        c.execute("""CREATE TABLE IF NOT EXISTS maker_quotes(
            ts TIMESTAMPTZ, sport VARCHAR, event_id BIGINT, market_ticker VARCHAR,
            side VARCHAR, action VARCHAR, price DOUBLE, size DOUBLE,
            fair DOUBLE, margin DOUBLE, alt BOOLEAN, reason VARCHAR, order_id VARCHAR)""")
        c.execute("""CREATE TABLE IF NOT EXISTS maker_fills(
            trade_id VARCHAR PRIMARY KEY, ts TIMESTAMPTZ, sport VARCHAR, event_id BIGINT,
            order_id VARCHAR, market_ticker VARCHAR, side VARCHAR, price DOUBLE,
            contracts DOUBLE, fee DOUBLE, ledger_worst_after DOUBLE)""")
        c.execute("""CREATE TABLE IF NOT EXISTS ledger_snapshots(
            ts TIMESTAMPTZ, sport VARCHAR, event_id BIGINT, worst_case DOUBLE,
            pnl_grid JSON, quotes_live INTEGER)""")
        c.execute("""CREATE TABLE IF NOT EXISTS maker_pnl(
            market_ticker VARCHAR PRIMARY KEY, ts TIMESTAMPTZ, sport VARCHAR,
            settled_pnl DOUBLE)""")


def log_quote(ts, sport, event_id, ticker, side, action, price, size, fair, margin, alt, reason, order_id):
    with connect(config.MAKER_DB_PATH) as c:
        c.execute("INSERT INTO maker_quotes VALUES (?,?,?,?,?,?,?,?,?,?,?,?,?)",
                  [ts, sport, event_id, ticker, side, action, price, size, fair, margin, alt, reason, order_id])


def log_fill(ts, sport, event_id, order_id, ticker, side, price, contracts, fee, worst_after, trade_id):
    with connect(config.MAKER_DB_PATH) as c:
        c.execute("INSERT OR IGNORE INTO maker_fills VALUES (?,?,?,?,?,?,?,?,?,?,?)",
                  [trade_id, ts, sport, event_id, order_id, ticker, side, price, contracts, fee, worst_after])


def log_ledger(ts, sport, event_id, worst_case, grid, quotes_live):
    with connect(config.MAKER_DB_PATH) as c:
        c.execute("INSERT INTO ledger_snapshots VALUES (?,?,?,?,?,?)",
                  [ts, sport, event_id, worst_case, json.dumps(grid), quotes_live])


def log_settlement(ts, sport, ticker, pnl):
    with connect(config.MAKER_DB_PATH) as c:
        c.execute("INSERT OR IGNORE INTO maker_pnl VALUES (?,?,?,?)", [ticker, ts, sport, pnl])
```

- [ ] **Step 4: Run tests**

Run: `python3 -m pytest unabated_edge/tests/test_maker_store.py unabated_edge/tests/test_config.py -q`
Expected: all pass

- [ ] **Step 5: Commit**

```bash
git add unabated_edge/config.py unabated_edge/maker/store.py unabated_edge/tests/test_maker_store.py
git commit -m "feat(unabated-edge): maker config params + sibling maker DB (store.py)"
```

---

### Task 3: Quote gateway (`maker/gateway.py`)

**Files:**
- Create: `unabated_edge/maker/gateway.py`
- Test: `unabated_edge/tests/test_maker_gateway.py`

**Interfaces:**
- Produces: `QuoteGateway` ABC with `is_live: bool`, `place(ticker, side, price_cents, count, client_order_id) -> str | None`, `cancel(order_id) -> bool`; `ShadowGateway` (fabricated ids `shadow-N`); `LiveGateway` (Kalshi REST); `make_gateway(mode, ack) -> QuoteGateway | None` (None when mode "off"/empty; `SystemExit` when live without ack=="1" or unknown mode).

- [ ] **Step 1: Write the failing tests**

```python
# unabated_edge/tests/test_maker_gateway.py
import pytest
from unabated_edge.maker import gateway


def test_make_gateway_modes():
    assert gateway.make_gateway("off", None) is None
    assert isinstance(gateway.make_gateway("shadow", None), gateway.ShadowGateway)
    with pytest.raises(SystemExit):
        gateway.make_gateway("live", None)         # dead-man switch
    assert isinstance(gateway.make_gateway("live", "1"), gateway.LiveGateway)
    with pytest.raises(SystemExit):
        gateway.make_gateway("bogus", None)


def test_shadow_place_and_cancel():
    gw = gateway.ShadowGateway()
    oid = gw.place("T-O25", "yes", 40, 30, "coid")
    assert oid == "shadow-1" and not gw.is_live
    assert gw.cancel(oid) is True


def test_live_place_payload_yes_and_no(monkeypatch):
    calls = []
    def fake_api(method, path, body=None, timeout=30):
        calls.append((method, path, body))
        return 201, {"order": {"order_id": "ord-1"}}, {}
    monkeypatch.setattr(gateway.auth_client, "api", fake_api)
    gw = gateway.LiveGateway()
    assert gw.place("T-O25", "yes", 40, 30, "c1") == "ord-1"
    assert gw.place("T-O25", "no", 57, 12, "c2") == "ord-1"
    (m1, p1, b1), (m2, p2, b2) = calls
    assert (m1, p1) == ("POST", "/portfolio/orders")
    assert b1 == {"ticker": "T-O25", "client_order_id": "c1", "action": "buy",
                  "side": "yes", "type": "limit", "count": 30, "yes_price": 40}
    assert b2["no_price"] == 57 and "yes_price" not in b2


def test_live_place_failure_returns_none(monkeypatch):
    monkeypatch.setattr(gateway.auth_client, "api", lambda m, p, body=None, timeout=30: (500, "err", {}))
    assert gateway.LiveGateway().place("T", "yes", 40, 30, "c") is None


def test_live_cancel(monkeypatch):
    seen = {}
    def fake_api(method, path, body=None, timeout=30):
        seen["call"] = (method, path)
        return 200, {}, {}
    monkeypatch.setattr(gateway.auth_client, "api", fake_api)
    assert gateway.LiveGateway().cancel("ord-9") is True
    assert seen["call"] == ("DELETE", "/portfolio/orders/ord-9")
```

- [ ] **Step 2: Run to verify failure**

Run: `python3 -m pytest unabated_edge/tests/test_maker_gateway.py -q`
Expected: FAIL (no module `gateway`)

- [ ] **Step 3: Implement**

```python
# unabated_edge/maker/gateway.py
"""Order-placement seam. Engine talks only to QuoteGateway; swapping
shadow<->live changes no engine code (interface + implementations —
the kalshi_mlb_mm pattern)."""
import logging
from abc import ABC, abstractmethod

from kalshi_common import auth_client

log = logging.getLogger("unabated_edge")


class QuoteGateway(ABC):
    is_live = False

    @abstractmethod
    def place(self, ticker: str, side: str, price_cents: int, count: int,
              client_order_id: str) -> str | None: ...

    @abstractmethod
    def cancel(self, order_id: str) -> bool: ...


class ShadowGateway(QuoteGateway):
    """Places nothing; fabricates order ids. The engine logs every decision to
    the maker DB regardless of gateway, so shadow mode costs no code."""

    def __init__(self):
        self._n = 0

    def place(self, ticker, side, price_cents, count, client_order_id):
        self._n += 1
        return f"shadow-{self._n}"

    def cancel(self, order_id):
        return True


class LiveGateway(QuoteGateway):
    is_live = True

    def place(self, ticker, side, price_cents, count, client_order_id):
        body = {"ticker": ticker, "client_order_id": client_order_id,
                "action": "buy", "side": side, "type": "limit", "count": int(count)}
        body["yes_price" if side == "yes" else "no_price"] = int(price_cents)
        status, resp, _ = auth_client.api("POST", "/portfolio/orders", body)
        if status not in (200, 201) or not isinstance(resp, dict):
            log.warning("maker place failed %s %s %dc x%d: status=%s resp=%s",
                        ticker, side, price_cents, count, status, resp)
            return None
        return (resp.get("order") or {}).get("order_id")

    def cancel(self, order_id):
        status, _, _ = auth_client.api("DELETE", f"/portfolio/orders/{order_id}")
        if status not in (200, 204):
            log.warning("maker cancel failed %s: status=%s", order_id, status)
            return False
        return True


def make_gateway(mode: str | None, ack: str | None) -> QuoteGateway | None:
    if not mode or mode == "off":
        return None
    if mode == "shadow":
        return ShadowGateway()
    if mode == "live":
        if ack != "1":
            raise SystemExit("MAKER_MODE=live requires MAKER_LIVE_ACK=1 (dead-man switch)")
        return LiveGateway()
    raise SystemExit(f"unknown MAKER_MODE={mode!r}")
```

- [ ] **Step 4: Run tests** — `python3 -m pytest unabated_edge/tests/test_maker_gateway.py -q` → `5 passed`

- [ ] **Step 5: Commit**

```bash
git add unabated_edge/maker/gateway.py unabated_edge/tests/test_maker_gateway.py
git commit -m "feat(unabated-edge): maker QuoteGateway — shadow + live with MAKER_LIVE_ACK dead-man"
```

---

### Task 4: Maker state + reconciliation (`maker/state.py`)

**Files:**
- Create: `unabated_edge/maker/state.py`
- Test: `unabated_edge/tests/test_maker_state.py`

**Interfaces:**
- Produces class `MakerState` with attrs `resting: dict[(ticker,side) -> {"order_id","price_cents","count"}]`, `ticker_info: dict[ticker -> {"sport","event_id","line"}]`, `our_orders: dict[order_id -> (ticker,side)]`, `fills: dict[event_id -> list[fill tuple]]`, `fills_by_ticker: dict[ticker -> net yes contracts]`, `fill_times: dict[event_id -> list[datetime]]`, `cooloff_until: dict[event_id -> datetime]`, `settled: set[ticker]`, `settled_pnl_today: float`; methods `roll_day(now)`, `register_ticker(...)`, `resting_for(ticker, side)`, `on_place(...)`, `on_cancel(ticker, side)`, `tickers_for(eid)`, `events_with_quotes()`, `events_with_exposure()`, `quotes_live()`, `quotes_live_for(eid)`, `exposure_fills(eid, exclude=None)`.
- Produces module functions (live-mode polling, all defensive on fp/cents shapes): `poll_fills(state, now) -> list[trade_id]`, `poll_positions(state) -> bool`, `poll_settlements(state, now)`.
- Consumes: `maker.store.log_fill/log_settlement`, `maker.ledger.worst_case`, `kalshi_common.auth_client.api`.

- [ ] **Step 1: Write the failing tests**

```python
# unabated_edge/tests/test_maker_state.py
import datetime
from unabated_edge import config
from unabated_edge.maker import state as mstate

_NOW = datetime.datetime(2026, 7, 11, 12, 0, tzinfo=datetime.timezone.utc)


def _state_with_order():
    s = mstate.MakerState()
    s.register_ticker("T-O25", "soccer", 1, 2.5)
    s.on_place("T-O25", "yes", "ord-1", 40, 30)
    return s


def test_exposure_includes_resting_as_fills():
    s = _state_with_order()
    s.fills[1] = [(2.5, "no", 10, 0.58)]
    exp = s.exposure_fills(1)
    assert (2.5, "yes", 30, 0.40) in exp and (2.5, "no", 10, 0.58) in exp
    assert s.exposure_fills(1, exclude=("T-O25", "yes")) == [(2.5, "no", 10, 0.58)]


def test_poll_fills_parses_updates_and_dedups(tmp_path, monkeypatch):
    monkeypatch.setattr(config, "MAKER_DB_PATH", tmp_path / "mk.duckdb")
    from unabated_edge.maker import store
    store.init()
    s = _state_with_order()
    payload = {"fills": [
        {"trade_id": "tr1", "order_id": "ord-1", "ticker": "T-O25", "side": "yes",
         "count_fp": "12.5", "yes_price_dollars": "0.4000", "no_price_dollars": "0.6000"},
        {"trade_id": "tr2", "order_id": "someone-else", "ticker": "T-O25", "side": "yes",
         "count": 5, "yes_price": 40},
    ]}
    monkeypatch.setattr(mstate.auth_client, "api", lambda m, p, body=None, timeout=30: (200, payload, {}))
    new = mstate.poll_fills(s, _NOW)
    assert new == ["tr1"]                                  # not ours -> ignored
    assert s.fills[1] == [(2.5, "yes", 12.5, 0.40)]
    assert s.fills_by_ticker["T-O25"] == 12.5
    assert s.resting[("T-O25", "yes")]["count"] == 17.5    # 30 - 12.5
    assert mstate.poll_fills(s, _NOW) == []                # dedup on trade_id


def test_poll_fills_full_fill_clears_resting(tmp_path, monkeypatch):
    monkeypatch.setattr(config, "MAKER_DB_PATH", tmp_path / "mk.duckdb")
    from unabated_edge.maker import store
    store.init()
    s = _state_with_order()
    payload = {"fills": [{"trade_id": "trX", "order_id": "ord-1", "ticker": "T-O25",
                          "side": "yes", "count": 30, "yes_price": 40}]}
    monkeypatch.setattr(mstate.auth_client, "api", lambda m, p, body=None, timeout=30: (200, payload, {}))
    mstate.poll_fills(s, _NOW)
    assert ("T-O25", "yes") not in s.resting


def test_poll_positions_detects_mismatch(monkeypatch):
    s = _state_with_order()
    s.fills_by_ticker["T-O25"] = 12.5
    monkeypatch.setattr(mstate.auth_client, "api",
        lambda m, p, body=None, timeout=30: (200, {"market_positions": [{"ticker": "T-O25", "position_fp": "12.50"}]}, {}))
    assert mstate.poll_positions(s) is True
    monkeypatch.setattr(mstate.auth_client, "api",
        lambda m, p, body=None, timeout=30: (200, {"market_positions": [{"ticker": "T-O25", "position": 3}]}, {}))
    assert mstate.poll_positions(s) is False


def test_poll_settlements_updates_daily_pnl(tmp_path, monkeypatch):
    monkeypatch.setattr(config, "MAKER_DB_PATH", tmp_path / "mk.duckdb")
    from unabated_edge.maker import store
    store.init()
    s = _state_with_order()
    s.fills_by_ticker["T-O25"] = 12.5
    payload = {"settlements": [{"ticker": "T-O25", "revenue": 1250,        # cents
                                "yes_total_cost": 500, "no_total_cost": 0}]}
    monkeypatch.setattr(mstate.auth_client, "api", lambda m, p, body=None, timeout=30: (200, payload, {}))
    mstate.poll_settlements(s, _NOW)
    assert round(s.settled_pnl_today, 2) == 7.50
    assert "T-O25" in s.settled and "T-O25" not in s.fills_by_ticker
    mstate.poll_settlements(s, _NOW)                       # idempotent
    assert round(s.settled_pnl_today, 2) == 7.50


def test_roll_day_resets_pnl():
    s = mstate.MakerState()
    s.roll_day(_NOW)
    s.settled_pnl_today = -100.0
    s.roll_day(_NOW + datetime.timedelta(days=1))
    assert s.settled_pnl_today == 0.0
```

- [ ] **Step 2: Run to verify failure** — `python3 -m pytest unabated_edge/tests/test_maker_state.py -q` → FAIL (no module)

- [ ] **Step 3: Implement**

```python
# unabated_edge/maker/state.py
"""In-memory maker state + live-mode reconciliation against Kalshi.

Kalshi is the source of truth for fills/positions (MLB lesson): we poll
/portfolio/fills each tick and /portfolio/positions + /portfolio/settlements
each minute, and a divergence between local bookkeeping and Kalshi positions
is a tripwire, never silently reconciled."""
import datetime
import logging

from unabated_edge import config
from unabated_edge.maker import ledger, store
from kalshi_common import auth_client

log = logging.getLogger("unabated_edge")


def _fp(d, key):
    """Numeric field: prefer `<key>_fp` STRING decimal, fall back to bare int."""
    v = d.get(f"{key}_fp")
    if v is None:
        v = d.get(key)
    try:
        return float(v) if v is not None else None
    except (TypeError, ValueError):
        return None


def _price_dollars(d, side):
    v = d.get(f"{side}_price_dollars")
    if v is not None:
        return float(v)
    v = d.get(f"{side}_price")
    return v / 100.0 if v is not None else None


def _money(d, key):
    v = d.get(f"{key}_dollars")
    if v is not None:
        return float(v)
    v = _fp(d, key)
    return v / 100.0 if v is not None else 0.0   # legacy bare value = cents


class MakerState:
    def __init__(self):
        self.resting = {}          # (ticker, side) -> {"order_id","price_cents","count"}
        self.ticker_info = {}      # ticker -> {"sport","event_id","line"}
        self.our_orders = {}       # order_id -> (ticker, side)
        self.fills = {}            # event_id -> [(line, side, contracts, price_dollars)]
        self.fills_by_ticker = {}  # ticker -> net yes contracts (positions check)
        self.fill_times = {}       # event_id -> [fill datetimes] (burst tripwire)
        self.cooloff_until = {}    # event_id -> datetime
        self.settled = set()       # tickers already settled (skip in positions check)
        self.settled_pnl_today = 0.0
        self._day = None
        self._fills_min_ts = None
        self._done_trades = set()

    def roll_day(self, now):
        if self._day != now.date():
            self._day = now.date()
            self.settled_pnl_today = 0.0

    def register_ticker(self, ticker, sport, event_id, line):
        self.ticker_info[ticker] = {"sport": sport, "event_id": event_id, "line": line}

    def resting_for(self, ticker, side):
        return self.resting.get((ticker, side))

    def on_place(self, ticker, side, order_id, price_cents, count):
        self.resting[(ticker, side)] = {"order_id": order_id, "price_cents": price_cents, "count": count}
        self.our_orders[order_id] = (ticker, side)

    def on_cancel(self, ticker, side):
        self.resting.pop((ticker, side), None)

    def tickers_for(self, eid):
        return [t for t, i in self.ticker_info.items() if i["event_id"] == eid]

    def events_with_quotes(self):
        return {self.ticker_info[t]["event_id"] for (t, _s) in self.resting if t in self.ticker_info}

    def events_with_exposure(self):
        return set(self.fills) | self.events_with_quotes()

    def quotes_live(self):
        return len(self.resting)

    def quotes_live_for(self, eid):
        return sum(1 for (t, _s) in self.resting
                   if self.ticker_info.get(t, {}).get("event_id") == eid)

    def exposure_fills(self, eid, exclude=None):
        """Committed fills PLUS all resting quotes treated as filled, so
        simultaneous fills across rungs can never breach the match cap."""
        out = list(self.fills.get(eid, []))
        for (t, s), r in self.resting.items():
            if (t, s) == exclude:
                continue
            info = self.ticker_info.get(t)
            if info and info["event_id"] == eid:
                out.append((info["line"], s, r["count"], r["price_cents"] / 100.0))
        return out


def poll_fills(state: MakerState, now) -> list[str]:
    min_ts = state._fills_min_ts or int(now.timestamp()) - 3600
    status, body, _ = auth_client.api("GET", f"/portfolio/fills?limit=100&min_ts={min_ts}")
    if status != 200 or not isinstance(body, dict):
        log.warning("maker poll_fills failed: status=%s", status)
        return []
    new = []
    for f in body.get("fills") or []:
        tid, oid = f.get("trade_id"), f.get("order_id")
        if not tid or tid in state._done_trades or oid not in state.our_orders:
            continue
        ticker, placed_side = state.our_orders[oid]
        info = state.ticker_info.get(ticker)
        if info is None:
            continue
        side = f.get("side") or placed_side
        n = _fp(f, "count") or 0.0
        price = _price_dollars(f, side)
        if price is None or n <= 0:
            log.warning("maker fill %s unparseable payload keys=%s", tid, sorted(f))
            continue
        eid = info["event_id"]
        state.fills.setdefault(eid, []).append((info["line"], side, n, price))
        state.fills_by_ticker[ticker] = state.fills_by_ticker.get(ticker, 0.0) + (n if side == "yes" else -n)
        state.fill_times.setdefault(eid, []).append(now)
        state._done_trades.add(tid)
        cur = state.resting.get((ticker, placed_side))
        if cur and cur["order_id"] == oid:
            cur["count"] -= n
            if cur["count"] <= 1e-9:
                state.resting.pop((ticker, placed_side), None)
        worst = ledger.worst_case(state.fills[eid])
        store.log_fill(now, info["sport"], eid, oid, ticker, side, price, n,
                       _fp(f, "fee") or 0.0, worst, tid)
        log.info("maker FILL %s %s %.2f@%.2f worst_after=%.2f", ticker, side, n, price, worst)
        new.append(tid)
    # overlap the next window; trade_id dedup absorbs the re-delivery
    state._fills_min_ts = int(now.timestamp()) - 60
    return new


def poll_positions(state: MakerState) -> bool:
    """True when Kalshi's net positions match our fills-derived book on every
    ticker we quoted (settled tickers excluded). API failure returns True —
    'cannot verify' must not false-trip the tripwire."""
    if not state.fills_by_ticker:
        return True
    status, body, _ = auth_client.api("GET", "/portfolio/positions?limit=1000")
    if status != 200 or not isinstance(body, dict):
        log.warning("maker poll_positions failed: status=%s", status)
        return True
    kalshi_pos = {p.get("ticker"): (_fp(p, "position") or 0.0)
                  for p in body.get("market_positions") or []}
    for ticker, expected in state.fills_by_ticker.items():
        if ticker in state.settled:
            continue
        actual = kalshi_pos.get(ticker, 0.0)
        if abs(actual - expected) > 0.01:
            log.warning("maker position mismatch %s: kalshi=%.2f local=%.2f", ticker, actual, expected)
            return False
    return True


def poll_settlements(state: MakerState, now):
    status, body, _ = auth_client.api("GET", "/portfolio/settlements?limit=100")
    if status != 200 or not isinstance(body, dict):
        return
    for s in body.get("settlements") or []:
        ticker = s.get("ticker")
        info = state.ticker_info.get(ticker)
        if info is None or ticker in state.settled:
            continue
        pnl = _money(s, "revenue") - _money(s, "yes_total_cost") - _money(s, "no_total_cost")
        state.settled.add(ticker)
        state.fills_by_ticker.pop(ticker, None)
        state.roll_day(now)
        state.settled_pnl_today += pnl
        store.log_settlement(now, info["sport"], ticker, pnl)
        log.info("maker SETTLED %s pnl=%.2f day_total=%.2f", ticker, pnl, state.settled_pnl_today)
```

- [ ] **Step 4: Run tests** — `python3 -m pytest unabated_edge/tests/test_maker_state.py -q` → `6 passed`

- [ ] **Step 5: Commit**

```bash
git add unabated_edge/maker/state.py unabated_edge/tests/test_maker_state.py
git commit -m "feat(unabated-edge): maker state — resting book, fills/positions/settlements reconciliation"
```

---

### Task 5: Quoting engine (`maker/engine.py`)

**Files:**
- Create: `unabated_edge/maker/engine.py`
- Test: `unabated_edge/tests/test_maker_engine.py`

**Interfaces:**
- Produces class `MakerEngine(gateway, state)` with `on_match(adapter, event_meta, kalshi_event, ladder, books, now)`, `note_success(sport, now)`, `watchdog(now)`, `sweep(sport, seen_eids, now)`, `pull_all(now, reason)`, `stats() -> dict`.
- Consumes: ladder `{line: {"p_over","book","alt","overround"}}` (soccer adapter shape), `books {ticker: {"yes_bids","no_bids"}}`, `EventMeta(event_id, league_key, start_utc)`, kalshi_event markets `{"ticker","strike_type","floor_strike"}`, Task 1-4 interfaces, `kalshi_common.ev_calc.maker_fee_per_contract(price_dollars)`, `venues.kalshi.yes_ask_from_book/no_ask_from_book`.

- [ ] **Step 1: Write the failing tests**

```python
# unabated_edge/tests/test_maker_engine.py
import datetime
import pytest
from unabated_edge import config
from unabated_edge.feed import EventMeta
from unabated_edge.maker import engine as mengine, state as mstate
from unabated_edge.sports.soccer import Soccer

_NOW = datetime.datetime(2026, 7, 11, 12, 0, tzinfo=datetime.timezone.utc)
_KICK = _NOW + datetime.timedelta(hours=3)
_EM = EventMeta(event_id=1, league_key="lg21", start_utc=_KICK)
_KEV = {"title": "A vs B: Regulation Time Total Goals",
        "markets": [{"ticker": "T-O25", "strike_type": "greater", "floor_strike": 2.5}]}
_LADDER = {2.5: {"p_over": 0.42, "book": 7, "alt": False, "overround": 1.05}}
# wide crowd: 0.30 yes bid / 0.55 yes ask (no bid 0.45) -> our quotes fit inside
_BOOK = {"yes_bids": [(0.30, 500.0)], "no_bids": [(0.45, 500.0)]}
_BOOKS = {"T-O25": _BOOK}


class FakeGateway:
    is_live = False
    def __init__(self):
        self.placed, self.cancelled, self._n = [], [], 0
    def place(self, ticker, side, price_cents, count, client_order_id):
        self._n += 1
        self.placed.append((ticker, side, price_cents, count))
        return f"o-{self._n}"
    def cancel(self, order_id):
        self.cancelled.append(order_id)
        return True


@pytest.fixture
def eng(tmp_path, monkeypatch):
    monkeypatch.setattr(config, "MAKER_DB_PATH", tmp_path / "mk.duckdb")
    from unabated_edge.maker import store
    store.init()
    return mengine.MakerEngine(FakeGateway(), mstate.MakerState())


def test_rests_two_sides_at_fair_minus_margin(eng):
    eng.on_match(Soccer(), _EM, _KEV, _LADDER, _BOOKS, _NOW)
    sides = {(s, p) for _t, s, p, _n in eng.gateway.placed}
    # fee(0.42)*100+1 ≈ 1.44c, roi 3%*42 = 1.26c -> m=ceil(1.44)=2c -> yes 40c
    # no side: fair 58c, fee≈1.43c+1, roi 1.74c -> m=2c -> 56c, capped by ask? no ask=1-0.30=0.70 -> cap 69 -> 56c
    assert ("yes", 40) in sides and ("no", 56) in sides


def test_holds_same_price_no_churn(eng):
    eng.on_match(Soccer(), _EM, _KEV, _LADDER, _BOOKS, _NOW)
    n = len(eng.gateway.placed)
    eng.on_match(Soccer(), _EM, _KEV, _LADDER, _BOOKS, _NOW + datetime.timedelta(seconds=5))
    assert len(eng.gateway.placed) == n and eng.gateway.cancelled == []


def test_requotes_on_fair_move(eng):
    eng.on_match(Soccer(), _EM, _KEV, _LADDER, _BOOKS, _NOW)
    moved = {2.5: {**_LADDER[2.5], "p_over": 0.47}}
    eng.on_match(Soccer(), _EM, _KEV, moved, _BOOKS, _NOW + datetime.timedelta(seconds=5))
    assert len(eng.gateway.cancelled) == 2          # both sides replaced
    assert any(p == 45 for _t, s, p, _n in eng.gateway.placed if s == "yes")  # 47 - 2c margin


def test_never_crosses_the_ask(eng):
    tight = {"T-O25": {"yes_bids": [(0.30, 500.0)], "no_bids": [(0.59, 500.0)]}}  # yes ask 41c
    eng.on_match(Soccer(), _EM, _KEV, _LADDER, tight, _NOW)
    yes = [p for _t, s, p, _n in eng.gateway.placed if s == "yes"]
    assert yes == [40]                              # min(40, 41-1) = 40 still fine
    tighter = {"T-O25": {"yes_bids": [(0.30, 500.0)], "no_bids": [(0.64, 500.0)]}}  # yes ask 36c
    eng.on_match(Soccer(), _EM, _KEV, _LADDER, tighter, _NOW + datetime.timedelta(seconds=5))
    # 36-1=35 < fair 42 - MAX_MARGIN 5 = 37 -> crowd tighter than we'll quote: cancel, don't rest
    assert any(o for o in eng.gateway.cancelled)
    assert not any(p == 35 for _t, s, p, _n in eng.gateway.placed if s == "yes")


def test_alt_rung_gated_on_overround(eng):
    bad_alt = {2.5: {"p_over": 0.42, "book": 7, "alt": True, "overround": 1.30}}
    eng.on_match(Soccer(), _EM, _KEV, bad_alt, _BOOKS, _NOW)
    assert eng.gateway.placed == []


def test_pull_window_cancels_everything(eng):
    eng.on_match(Soccer(), _EM, _KEV, _LADDER, _BOOKS, _NOW)
    near_kick = _KICK - datetime.timedelta(minutes=2)
    eng.on_match(Soccer(), _EM, _KEV, _LADDER, _BOOKS, near_kick)
    assert len(eng.gateway.cancelled) == 2 and eng.state.quotes_live() == 0


def test_fill_burst_triggers_cooloff(eng):
    eng.state.fill_times[1] = [_NOW - datetime.timedelta(seconds=i) for i in range(config.FILL_BURST_N + 1)]
    eng.on_match(Soccer(), _EM, _KEV, _LADDER, _BOOKS, _NOW)
    assert eng.gateway.placed == []
    assert eng.state.cooloff_until[1] > _NOW


def test_ledger_caps_size(eng, monkeypatch):
    monkeypatch.setattr(config, "MATCH_CAP_PCT", 0.02)      # $20 cap
    eng.on_match(Soccer(), _EM, _KEV, _LADDER, _BOOKS, _NOW)
    # yes quote at 40c: alone it could take floor(20/0.40)=50, but the no-side
    # resting quote is also assumed filled; every placement stays within cap
    for _t, s, p, n in eng.gateway.placed:
        assert n * (p / 100.0) <= 300.0 + 1e-9              # never above MAX_QUOTE_PCT
    from unabated_edge.maker import ledger
    assert ledger.worst_case(eng.state.exposure_fills(1)) >= -20.0 - 1e-9


def test_watchdog_pulls_on_stale_feed(eng):
    eng.on_match(Soccer(), _EM, _KEV, _LADDER, _BOOKS, _NOW)
    eng.note_success("soccer", _NOW)
    eng.watchdog(_NOW + datetime.timedelta(seconds=config.MAX_STALENESS_SEC + 5))
    assert eng.state.quotes_live() == 0


def test_sweep_pulls_unpaired_events(eng):
    eng.on_match(Soccer(), _EM, _KEV, _LADDER, _BOOKS, _NOW)
    eng.sweep("soccer", set(), _NOW)                        # event 1 not seen this tick
    assert eng.state.quotes_live() == 0


def test_daily_halt_pulls_everything(eng, monkeypatch):
    eng.on_match(Soccer(), _EM, _KEV, _LADDER, _BOOKS, _NOW)
    eng.state.roll_day(_NOW)
    eng.state.settled_pnl_today = -config.DAILY_LOSS_HALT_PCT * config.BANKROLL - 1
    eng.on_match(Soccer(), _EM, _KEV, _LADDER, _BOOKS, _NOW + datetime.timedelta(seconds=5))
    assert eng.state.quotes_live() == 0
    assert eng.stats()["halted"] is True
```

- [ ] **Step 2: Run to verify failure** — `python3 -m pytest unabated_edge/tests/test_maker_engine.py -q` → FAIL (no module)

- [ ] **Step 3: Implement**

```python
# unabated_edge/maker/engine.py
"""Per-match quoting brain: consumes the tick's devigged ladder + fresh books,
maintains resting quotes through the gateway, enforces the cap stack via the
goal-grid ledger. All decisions are logged to the maker DB (measurement-first)."""
import datetime
import logging
import math

from unabated_edge import config
from unabated_edge.maker import ledger, store
from unabated_edge.venues.kalshi import yes_ask_from_book, no_ask_from_book
from kalshi_common.ev_calc import maker_fee_per_contract

log = logging.getLogger("unabated_edge")
_DT_MIN = datetime.datetime.min.replace(tzinfo=datetime.timezone.utc)


class MakerEngine:
    def __init__(self, gateway, state):
        self.gateway = gateway
        self.state = state
        self.last_success = {}      # sport -> last successful tick (watchdog)
        self._last_skip = {}        # (ticker, side) -> reason (dedup skip rows)
        self._halted = False

    # ---------- gates ----------

    def _daily_halted(self, now):
        self.state.roll_day(now)
        self._halted = self.state.settled_pnl_today <= -config.DAILY_LOSS_HALT_PCT * config.BANKROLL
        return self._halted

    def _fill_burst(self, eid, now):
        times = self.state.fill_times.get(eid, [])
        return len([t for t in times if (now - t).total_seconds() < 60]) > config.FILL_BURST_N

    # ---------- quote math ----------

    def _margin_cents(self, fair_cents, alt):
        p = max(0.01, fair_cents / 100.0)
        m = max(maker_fee_per_contract(p) * 100 + config.PICKOFF_BUFFER_CENTS,
                config.ROI_MARGIN * fair_cents)
        if alt:
            m *= config.ALT_MARGIN_MULT
        return math.ceil(m)

    def _global_room(self, eid):
        committed = 0.0
        for other in self.state.events_with_exposure():
            if other != eid:
                committed += max(0.0, -ledger.worst_case(self.state.exposure_fills(other)))
        return config.GLOBAL_CAP_PCT * config.BANKROLL - committed

    def _desired(self, eid, ticker, line, rung, book, side):
        """((price_cents, count), None) or (None, skip_reason)."""
        alt = bool(rung.get("alt"))
        if alt and not (config.ALT_OVERROUND_MIN <= (rung.get("overround") or 0) <= config.ALT_OVERROUND_MAX):
            return None, "alt_overround"
        fair = rung["p_over"] if side == "yes" else 1 - rung["p_over"]
        fair_cents = fair * 100
        price = math.floor(fair_cents + 1e-9) - self._margin_cents(fair_cents, alt)
        opp_ask = yes_ask_from_book(book) if side == "yes" else no_ask_from_book(book)
        if opp_ask is not None:
            price = min(price, int(round(opp_ask * 100)) - 1)     # never cross
        if price < 1:
            return None, "price_floor"
        if price < fair_cents - config.MAX_MARGIN_CENTS - 1e-9:
            return None, "crowd_tighter"
        pd = price / 100.0
        n = math.floor(config.MAX_QUOTE_PCT * config.BANKROLL / pd)
        if alt:
            n = math.floor(n * config.ALT_SIZE_MULT)
        budget = min(config.MATCH_CAP_PCT * config.BANKROLL, self._global_room(eid))
        if budget <= 0:
            return None, "global_cap"
        n = min(n, ledger.max_contracts(
            self.state.exposure_fills(eid, exclude=(ticker, side)), line, side, pd, budget))
        if n <= 0:
            return None, "no_room"
        return (price, n), None

    # ---------- order sync ----------

    def _sync(self, sport, eid, ticker, side, desired, reason, fair, alt, now):
        cur = self.state.resting_for(ticker, side)
        if desired is None:
            if cur:
                if self.gateway.cancel(cur["order_id"]):
                    self.state.on_cancel(ticker, side)
                    store.log_quote(now, sport, eid, ticker, side, "cancel", None, None,
                                    fair, None, alt, reason, cur["order_id"])
            elif self._last_skip.get((ticker, side)) != reason:
                store.log_quote(now, sport, eid, ticker, side, "skip", None, None,
                                fair, None, alt, reason, None)
            self._last_skip[(ticker, side)] = reason
            return
        self._last_skip.pop((ticker, side), None)
        price, count = desired
        if cur and cur["price_cents"] == price:
            return                                      # hold: queue position is capital
        action = "replace" if cur else "rest"
        if cur:
            if not self.gateway.cancel(cur["order_id"]):
                return                                  # keep local state; retry next tick
            self.state.on_cancel(ticker, side)
        coid = f"uemk-{ticker}-{side}-{price}-{int(now.timestamp())}"
        oid = self.gateway.place(ticker, side, price, count, coid)
        margin = round(fair - price / 100.0, 4)
        if oid is None:
            store.log_quote(now, sport, eid, ticker, side, "skip", price / 100.0, count,
                            fair, margin, alt, "place_failed", None)
            return
        self.state.on_place(ticker, side, oid, price, count)
        store.log_quote(now, sport, eid, ticker, side, action, price / 100.0, count,
                        fair, margin, alt, "quote", oid)

    # ---------- pulls ----------

    def _pull_match(self, eid, now, reason):
        for ticker in self.state.tickers_for(eid):
            info = self.state.ticker_info[ticker]
            for side in ("yes", "no"):
                cur = self.state.resting_for(ticker, side)
                if cur and self.gateway.cancel(cur["order_id"]):
                    self.state.on_cancel(ticker, side)
                    store.log_quote(now, info["sport"], eid, ticker, side, "cancel",
                                    None, None, None, None, None, reason, cur["order_id"])

    def pull_all(self, now, reason):
        for eid in list(self.state.events_with_quotes()):
            self._pull_match(eid, now, reason)

    def note_success(self, sport, now):
        self.last_success[sport] = now

    def watchdog(self, now):
        """Feed-staleness guard: if any adapter hasn't completed a successful
        tick within MAX_STALENESS_SEC, our fair is dark — pull everything.
        Called every main-loop iteration, including after tick exceptions."""
        for _sport, ts in self.last_success.items():
            if (now - ts).total_seconds() > config.MAX_STALENESS_SEC:
                if self.state.quotes_live():
                    log.warning("maker feed stale — pulling all quotes")
                self.pull_all(now, "feed_stale")
                return

    def sweep(self, sport, seen_eids, now):
        """Pull quotes on events that stopped appearing in the paired pre-kick
        set (kickoff passed, pairing lost, market closed)."""
        for eid in list(self.state.events_with_quotes()):
            tickers = self.state.tickers_for(eid)
            if tickers and self.state.ticker_info[tickers[0]]["sport"] == sport and eid not in seen_eids:
                self._pull_match(eid, now, "unpaired")

    def stats(self):
        committed = sum(max(0.0, -ledger.worst_case(self.state.exposure_fills(e)))
                        for e in self.state.events_with_exposure())
        return {"quotes_live": self.state.quotes_live(),
                "worst_total": round(committed, 2), "halted": self._halted}

    # ---------- main entry ----------

    def on_match(self, adapter, event_meta, kalshi_event, ladder, books, now):
        eid = event_meta.event_id
        sport = adapter.sport
        if self._daily_halted(now):
            return self.pull_all(now, "daily_halt")
        if event_meta.start_utc is not None and \
                (event_meta.start_utc - now).total_seconds() < config.QUOTE_PULL_MIN * 60:
            return self._pull_match(eid, now, "pull_window")
        if now < self.state.cooloff_until.get(eid, _DT_MIN):
            return self._pull_match(eid, now, "cooloff")
        if self._fill_burst(eid, now):
            self.state.cooloff_until[eid] = now + datetime.timedelta(minutes=config.COOLOFF_MIN)
            log.warning("maker fill burst on event %s — %.0fmin cooloff", eid, config.COOLOFF_MIN)
            return self._pull_match(eid, now, "fill_burst")
        if not ladder:
            return self._pull_match(eid, now, "anchor_gone")
        for mk in kalshi_event.get("markets", []):
            book = books.get(mk.get("ticker"))
            if book is None:
                continue
            ya, na = yes_ask_from_book(book), no_ask_from_book(book)
            if ya is not None and na is not None and ya + na < 1 - 2 * maker_fee_per_contract(0.5):
                log.warning("maker crossed book on %s — pulling match", mk.get("ticker"))
                return self._pull_match(eid, now, "crossed_book")
        for mk in kalshi_event.get("markets", []):
            if mk.get("strike_type") != "greater":
                continue
            try:
                line = float(mk.get("floor_strike"))
            except (TypeError, ValueError):
                continue
            rung, ticker = ladder.get(line), mk.get("ticker")
            if rung is None or not ticker:
                continue
            book = books.get(ticker)
            if book is None:
                continue                        # no book this tick: hold existing quotes
            self.state.register_ticker(ticker, sport, eid, line)
            for side in ("yes", "no"):
                fair = rung["p_over"] if side == "yes" else 1 - rung["p_over"]
                desired, reason = self._desired(eid, ticker, line, rung, book, side)
                self._sync(sport, eid, ticker, side, desired, reason, fair,
                           bool(rung.get("alt")), now)
        exp = self.state.exposure_fills(eid)
        store.log_ledger(now, sport, eid, ledger.worst_case(exp),
                         ledger.pnl_grid(exp) if exp else [], self.state.quotes_live_for(eid))
```

- [ ] **Step 4: Run tests** — `python3 -m pytest unabated_edge/tests/test_maker_engine.py -q` → `12 passed`. Adjust the two price-assertion tests ONLY if the fee constant makes the margin 1c different — recompute by hand first (`maker_fee_per_contract(0.42)` in a REPL), never loosen an assertion to `>`.

- [ ] **Step 5: Commit**

```bash
git add unabated_edge/maker/engine.py unabated_edge/tests/test_maker_engine.py
git commit -m "feat(unabated-edge): maker quoting engine — fair±margin, ledger caps, pulls, tripwires"
```

---

### Task 6: Runner + adapter integration

**Files:**
- Modify: `unabated_edge/sports/base.py` (add `fair_ladder` default)
- Modify: `unabated_edge/sports/soccer.py` (implement `fair_ladder`)
- Modify: `unabated_edge/runner.py` (maker wiring: run_tick param, main_loop lifecycle, watchdog, heartbeat)
- Test: `unabated_edge/tests/test_runner_tick.py` (extend), `unabated_edge/tests/test_soccer_adapter.py` (extend)

**Interfaces:**
- Produces: `SportAdapter.fair_ladder(state, event_meta) -> dict | None` (default None); `run_tick(..., maker=None)`; maker lifecycle in `main_loop` driven by `config.MAKER_MODE`.
- Consumes: `MakerEngine` API from Task 5, `make_gateway` from Task 3, polling functions from Task 4.

- [ ] **Step 1: Write the failing tests**

Append to `unabated_edge/tests/test_soccer_adapter.py`:

```python
def test_fair_ladder_exposes_anchor_ladder():
    from unabated_edge.tests.test_runner_tick import _state, _NOW  # reuse fixture helpers
    st = _state()
    em = next(iter(st.events.values()))
    ladder = Soccer().fair_ladder(st, em)
    assert 2.5 in ladder and 0 < ladder[2.5]["p_over"] < 1


def test_fair_ladder_none_when_no_anchor():
    from unabated_edge.tests.test_runner_tick import _state
    st = _state()
    st.lines.clear()
    em = next(iter(st.events.values()))
    assert Soccer().fair_ladder(st, em) is None
```

Append to `unabated_edge/tests/test_runner_tick.py`:

```python
class _FakeMaker:
    def __init__(self):
        self.calls, self.sweeps = [], []
    def on_match(self, adapter, event_meta, kev, ladder, books, now):
        self.calls.append((event_meta.event_id, ladder is not None, sorted(books)))
    def sweep(self, sport, seen, now):
        self.sweeps.append((sport, set(seen)))


def test_maker_hook_receives_ladder_and_books(tmp_path, monkeypatch):
    _init_dbs(tmp_path, monkeypatch)
    mk = _FakeMaker()
    runner.run_tick(Soccer(), _state(), [_KEV], now=_NOW, dry_run=True,
                    book_fn=lambda t: _BOOK, maker=mk)
    assert mk.calls == [(1, True, ["T-O25"])]
    assert mk.sweeps == [("soccer", {1})]


def test_maker_hook_skipped_when_none(tmp_path, monkeypatch):
    _init_dbs(tmp_path, monkeypatch)
    rows = runner.run_tick(Soccer(), _state(), [_KEV], now=_NOW, dry_run=True,
                           book_fn=lambda t: _BOOK)     # maker defaults to None
    assert any(rows)
```

- [ ] **Step 2: Run to verify failure** — `python3 -m pytest unabated_edge/tests/test_runner_tick.py unabated_edge/tests/test_soccer_adapter.py -q` → new tests FAIL (`fair_ladder` missing / unexpected kwarg `maker`)

- [ ] **Step 3: Implement**

`unabated_edge/sports/base.py` — add method to `SportAdapter` after `price_event`:

```python
    def fair_ladder(self, state, event_meta) -> dict | None:
        """Devigged {line: {p_over, book, alt, overround}} for the maker.
        None when this adapter/event has no anchored ladder this tick
        (adapters without maker support simply inherit this default)."""
        return None
```

`unabated_edge/sports/soccer.py` — add method to `Soccer` after `_anchor_totals`:

```python
    def fair_ladder(self, state, event_meta):
        return self._anchor_ladder(state, event_meta.event_id) or None
```

`unabated_edge/runner.py` changes:

1. Signature: `def run_tick(adapter, state, kalshi_events, *, now, dry_run=True, book_fn, trades_fn=None, maker=None) -> list[dict]:`
2. In the pairing loop, right after the trades-poll block and BEFORE the `tipoff_ok` check, add (and collect `seen_eids`):

```python
        if maker is not None:
            maker.on_match(adapter, event_meta, kev,
                           adapter.fair_ladder(state, event_meta), books, now)
```

   Initialize `seen_eids = set()` before the pairing loop; add `seen_eids.add(event_meta.event_id)` immediately AFTER the pre-kick `continue` (post-kick events are intentionally NOT added, so sweep pulls their quotes). At the end of `run_tick`, before `return flagged`:

```python
    if maker is not None:
        maker.sweep(adapter.sport, seen_eids, now)
```

   Note: `seen_eids.add(...)` goes AFTER the `now >= start_utc` continue, i.e. only pre-kick paired events count as seen.
3. `main_loop` — build the maker once before the loop (imports at top of file: `from unabated_edge.maker import engine as maker_engine, gateway as maker_gateway, state as maker_state, store as maker_store`):

```python
    maker = None
    gw = maker_gateway.make_gateway(config.MAKER_MODE, config.MAKER_LIVE_ACK)
    if gw is not None:
        maker_store.init()
        maker = maker_engine.MakerEngine(gw, maker_state.MakerState())
        log.info("maker enabled mode=%s live=%s", config.MAKER_MODE, gw.is_live)
```

4. Inside the per-adapter loop, after the `run_tick(...)` call: `if maker: maker.note_success(a.sport, now)`.
5. After `storage.flush()`, still inside the `try`:

```python
            if maker is not None and gw.is_live:
                maker_state.poll_fills(maker.state, now)
                if time.time() - last_recon > 60:
                    if not maker_state.poll_positions(maker.state):
                        maker.pull_all(now, "position_mismatch")
                    maker_state.poll_settlements(maker.state, now)
                    last_recon = time.time()
```

   (initialize `last_recon = 0.0` next to `last_k`.)
6. After the whole `try/except` block, before `time.sleep(...)` — runs even when the tick raised:

```python
        if maker is not None:
            maker.watchdog(datetime.datetime.now(datetime.timezone.utc))
```

7. Heartbeat: extend the log line with `" maker=%s"` and value `maker.stats() if maker else None`.

- [ ] **Step 4: Run the FULL suite** — `python3 -m pytest unabated_edge/tests/ -q`
Expected: all pass (55 existing + ~33 new). The existing `run_tick` tests pass unchanged because `maker` defaults to `None`.

- [ ] **Step 5: Commit**

```bash
git add unabated_edge/runner.py unabated_edge/sports/base.py unabated_edge/sports/soccer.py unabated_edge/tests/test_runner_tick.py unabated_edge/tests/test_soccer_adapter.py
git commit -m "feat(unabated-edge): wire maker into the tick loop — fair_ladder hook, watchdog, sweep, heartbeat"
```

---

### Task 7: Docs + full verification

**Files:**
- Modify: `unabated_edge/README.md` (maker section)
- Modify: `CLAUDE.md` (unabated_edge bullet: add maker)

- [ ] **Step 1: README** — add a `## Market maker (maker/)` section after "Architecture" covering: what it does (rests fair±margin quotes on all anchor rungs), the cap stack table (% of bankroll: quote 30% / match 40% / global 75% / daily halt 40%), MAKER_MODE off/shadow/live + MAKER_LIVE_ACK, the goal-grid ledger in two sentences, pull triggers (fair move, feed watchdog, 3-min pull window, cooloff, tripwires), the maker DB tables, and a runbook line (`MAKER_MODE=live MAKER_LIVE_ACK=1 python3 -m unabated_edge.runner`; kill file pulls quotes via watchdog on shutdown — verify `.kill` behavior note). Add all new config rows to the config table.
- [ ] **Step 2: CLAUDE.md** — extend the unabated_edge bullet: "…now includes an in-process market maker (`unabated_edge/maker/`) quoting KXWCTOTAL around the devigged anchor at fair±margin, goal-grid worst-case ledger caps (% of bankroll), live/shadow QuoteGateway with MAKER_LIVE_ACK dead-man; MAKER_MODE=off by default."
- [ ] **Step 3: Full suite + import smoke**

```bash
python3 -m pytest unabated_edge/tests/ -q
python3 -c "from unabated_edge import runner"    # import-time syntax check
```
Expected: all pass, clean import.

- [ ] **Step 4: Commit**

```bash
git add unabated_edge/README.md CLAUDE.md
git commit -m "docs(unabated-edge): maker section — architecture, cap stack, runbook"
```

---

## After the plan: NOT in scope for executors

The executive pre-merge review (`git diff main..HEAD`), the user merge/live gate, `.env` edit (`MAKER_MODE=live`, `MAKER_LIVE_ACK=1`), runner restart from main, and capture-DB folding are the orchestrating session's job — never an executor's.
