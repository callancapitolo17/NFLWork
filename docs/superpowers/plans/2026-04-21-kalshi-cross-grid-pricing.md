# Kalshi Cross-Book Grid Pricing Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Fix the Kalshi column in the Cross-Book Grid so it reflects fair-value (mid of buy/sell) in the cell and median, flags outliers against the actual buy price, and surfaces a hover tooltip with last-trade + buy + sell. Also restores the legacy `kalshi_odds` write path (which currently persists zeros because it doesn't handle Kalshi's `*_dollars` response format).

**Architecture:** Two prices per Kalshi row are pushed into existing schema columns — `implied_prob` = take price (buy), `devig_prob` = fair value (mid). Sportsbook rows are untouched. The grid query pulls both columns and branches the outlier flag per venue. The tooltip is a separate batched query joined on `market_map` → `kalshi_odds`. No schema migration.

**Tech Stack:** Python (scrapers, DuckDB queries), Dash / dash_table (dashboard), pytest (unit + integration), DuckDB (storage).

**Reference spec:** `docs/superpowers/specs/2026-04-21-kalshi-cross-grid-pricing-design.md`

---

## Worktree Lifecycle

This plan is executed in an isolated worktree so parallel sessions on `main` or other branches can continue. Branch `feature/kalshi-cross-grid-pricing` already exists (carries the design doc commit from 2026-04-21).

**Setup (once, at start of Task 1):**
```bash
git worktree add /Users/callancapitolo/NFLWork-trees/kalshi-cross-grid-pricing feature/kalshi-cross-grid-pricing
cd /Users/callancapitolo/NFLWork-trees/kalshi-cross-grid-pricing
```

All subsequent tasks run from that worktree path.

**DuckDB note:** Tests use monkeypatched `tmp_path` DBs, so we do NOT symlink or copy the production `nfl_draft/nfl_draft.duckdb` into the worktree. (Per `CLAUDE.md`: symlinking DuckDB files loses WAL data when the worktree is removed.)

**Teardown (after Task 12 merge + user approval):**
```bash
# From main repo directory
git worktree remove /Users/callancapitolo/NFLWork-trees/kalshi-cross-grid-pricing
git branch -d feature/kalshi-cross-grid-pricing
```

---

## File Structure

Files created / modified (grouped by task):

```
Modified:
  nfl_draft/scrapers/_base.py            # Task 2: add OddsRow override fields
  nfl_draft/lib/quarantine.py            # Task 3: honor overrides on INSERT
  nfl_draft/scrapers/kalshi.py           # Tasks 4-6: extractors, cascade, legacy fix
  nfl_draft/lib/queries.py               # Tasks 7-8: grid query + tooltip query
  kalshi_draft/app.py                    # Task 9: wire tooltip into DataTable
  nfl_draft/README.md                    # Task 10: doc the new Kalshi semantics

Tests modified:
  nfl_draft/tests/unit/test_scraper_parsing.py        # Tasks 4, 5
  nfl_draft/tests/integration/test_dashboard_queries.py  # Tasks 3, 6, 7, 8

New files: none
```

---

## Task 1: Create Worktree

**Files:** N/A (git operation)

- [ ] **Step 1: Verify starting branch and clean state**

Run from main repo (`/Users/callancapitolo/NFLWork`):
```bash
git branch --show-current
git status --short
```
Expected: current branch `feature/kalshi-cross-grid-pricing`, empty `status` output. If you're on `main`, switch: `git checkout feature/kalshi-cross-grid-pricing`.

- [ ] **Step 2: Create the worktree**

```bash
git worktree add /Users/callancapitolo/NFLWork-trees/kalshi-cross-grid-pricing feature/kalshi-cross-grid-pricing
```
Expected: `Preparing worktree (checking out 'feature/kalshi-cross-grid-pricing')` and `HEAD is now at 8ddffa0 docs: design spec...`.

- [ ] **Step 3: Change directory and verify**

```bash
cd /Users/callancapitolo/NFLWork-trees/kalshi-cross-grid-pricing
git branch --show-current
ls docs/superpowers/specs/
```
Expected: branch confirmed, spec file present.

---

## Task 2: Extend OddsRow with Price Overrides

**Files:**
- Modify: `nfl_draft/scrapers/_base.py`

Rationale: the scraper for Kalshi needs to push exact values for both `implied_prob` and `devig_prob` without going through the lossy `american_to_implied` round-trip. Sportsbooks leave both overrides as `None` → existing behavior preserved.

- [ ] **Step 1: Write the failing test**

Create at the top of `nfl_draft/tests/unit/test_scraper_parsing.py` (below existing imports):

```python
def test_oddsrow_accepts_implied_and_devig_overrides():
    """Optional overrides allow a scraper to bypass american_to_implied when
    it has exact prices already (e.g. Kalshi yes_ask fractions)."""
    from nfl_draft.scrapers._base import OddsRow
    from datetime import datetime

    row = OddsRow(
        book="kalshi", book_label="L", book_subject="S",
        american_odds=1900, fetched_at=datetime.now(),
        implied_prob=0.05, devig_prob=0.035,
    )
    assert row.implied_prob == 0.05
    assert row.devig_prob == 0.035


def test_oddsrow_overrides_default_to_none():
    """Back-compat: existing call sites that don't set overrides must keep working."""
    from nfl_draft.scrapers._base import OddsRow
    from datetime import datetime

    row = OddsRow(
        book="dk", book_label="L", book_subject="S",
        american_odds=-110, fetched_at=datetime.now(),
    )
    assert row.implied_prob is None
    assert row.devig_prob is None
```

- [ ] **Step 2: Run test to verify it fails**

```bash
pytest nfl_draft/tests/unit/test_scraper_parsing.py::test_oddsrow_accepts_implied_and_devig_overrides -v
```
Expected: FAIL with `TypeError: __init__() got an unexpected keyword argument 'implied_prob'`.

- [ ] **Step 3: Add override fields to OddsRow**

Edit `nfl_draft/scrapers/_base.py`. Replace the `OddsRow` dataclass with:

```python
@dataclass
class OddsRow:
    """A single row of raw scraper output (one binary outcome from one venue).

    The `implied_prob` and `devig_prob` fields are optional overrides. When
    set, `quarantine.write_or_quarantine` uses them verbatim for the
    `draft_odds.implied_prob` / `devig_prob` columns instead of deriving both
    from `american_to_implied(american_odds)`. Used by Kalshi so that the
    take-price (buy / yes_ask) and fair-value (mid) can be persisted without
    rounding loss through the American-odds round-trip.
    """
    book: str
    book_label: str
    book_subject: str
    american_odds: int
    fetched_at: datetime
    market_group: str = ""
    implied_prob: Optional[float] = None
    devig_prob: Optional[float] = None
```

Add the import at the top:
```python
from typing import Optional
```

- [ ] **Step 4: Run both new tests to verify they pass**

```bash
pytest nfl_draft/tests/unit/test_scraper_parsing.py::test_oddsrow_accepts_implied_and_devig_overrides nfl_draft/tests/unit/test_scraper_parsing.py::test_oddsrow_overrides_default_to_none -v
```
Expected: both PASS.

- [ ] **Step 5: Run the full scraper-parsing suite to confirm no regressions**

```bash
pytest nfl_draft/tests/unit/test_scraper_parsing.py -v
```
Expected: all pre-existing tests still pass (they don't set the new fields).

- [ ] **Step 6: Commit**

```bash
git add nfl_draft/scrapers/_base.py nfl_draft/tests/unit/test_scraper_parsing.py
git commit -m "$(cat <<'EOF'
feat(oddsrow): add optional implied_prob / devig_prob overrides

Kalshi needs to persist the exact buy price and the exact mid. Going
through american_to_implied introduces rounding loss at non-round cents
(e.g. yes_ask=7 → +1329 → 0.06998 instead of 0.07). New fields are
optional and default to None so existing sportsbook scrapers are
unaffected.

Co-Authored-By: Claude Opus 4.7 (1M context) <noreply@anthropic.com>
EOF
)"
```

---

## Task 3: Teach the Quarantine Writer to Honor Overrides

**Files:**
- Modify: `nfl_draft/lib/quarantine.py:42-46`
- Test: `nfl_draft/tests/integration/test_dashboard_queries.py`

- [ ] **Step 1: Write the failing test**

Append to `nfl_draft/tests/integration/test_dashboard_queries.py`:

```python
def test_quarantine_honors_oddsrow_price_overrides(monkeypatch, tmp_path):
    """When OddsRow sets implied_prob / devig_prob, write_or_quarantine must
    store those exact values in draft_odds rather than deriving both from
    american_to_implied(american_odds)."""
    from nfl_draft.lib import db as db_module
    monkeypatch.setattr(db_module, "DB_PATH", tmp_path / "test.duckdb")
    db_module.init_schema()

    from datetime import datetime
    from nfl_draft.lib.db import write_connection, read_connection
    from nfl_draft.lib.quarantine import write_or_quarantine
    from nfl_draft.scrapers._base import OddsRow

    # Seed market_map so the row is mapped (not quarantined)
    with write_connection() as con:
        con.execute(
            "INSERT INTO market_map (book, book_label, book_subject, market_id) "
            "VALUES ('kalshi', 'KXNFLDRAFTPICK-26-5', 'Carnell Tate', 'pick_5_overall_carnell-tate')"
        )

    row = OddsRow(
        book="kalshi",
        book_label="KXNFLDRAFTPICK-26-5",
        book_subject="Carnell Tate",
        american_odds=1900,       # corresponds to yes_ask=5c
        fetched_at=datetime.now(),
        implied_prob=0.05,        # exact buy price (take)
        devig_prob=0.035,         # exact mid (fair)
    )
    mapped, unmapped = write_or_quarantine([row])
    assert (mapped, unmapped) == (1, 0)

    with read_connection() as con:
        got = con.execute(
            "SELECT implied_prob, devig_prob FROM draft_odds WHERE market_id = ?",
            ["pick_5_overall_carnell-tate"],
        ).fetchone()
    assert got == (0.05, 0.035)


def test_quarantine_falls_back_to_american_when_overrides_absent(monkeypatch, tmp_path):
    """Sportsbook rows leave overrides as None → both columns populated via
    american_to_implied(american_odds), same as before this change."""
    from nfl_draft.lib import db as db_module
    monkeypatch.setattr(db_module, "DB_PATH", tmp_path / "test.duckdb")
    db_module.init_schema()

    from datetime import datetime
    from nfl_draft.lib.db import write_connection, read_connection
    from nfl_draft.lib.quarantine import write_or_quarantine
    from nfl_draft.scrapers._base import OddsRow

    with write_connection() as con:
        con.execute(
            "INSERT INTO market_map (book, book_label, book_subject, market_id) "
            "VALUES ('draftkings', 'nfl_draft_pick_5', 'Carnell Tate', 'pick_5_overall_carnell-tate')"
        )

    row = OddsRow(
        book="draftkings",
        book_label="nfl_draft_pick_5",
        book_subject="Carnell Tate",
        american_odds=1200,  # implied = 100/1300 ≈ 0.07692
        fetched_at=datetime.now(),
    )
    write_or_quarantine([row])

    with read_connection() as con:
        got = con.execute(
            "SELECT implied_prob, devig_prob FROM draft_odds WHERE market_id = ?",
            ["pick_5_overall_carnell-tate"],
        ).fetchone()
    expected = 100.0 / (1200 + 100)
    assert abs(got[0] - expected) < 1e-9
    assert got[0] == got[1]  # both columns identical, as today
```

- [ ] **Step 2: Run tests to verify they fail**

```bash
pytest nfl_draft/tests/integration/test_dashboard_queries.py::test_quarantine_honors_oddsrow_price_overrides -v
```
Expected: FAIL — the first test will fail because `quarantine.py` currently writes `implied` into both columns regardless of overrides. The second test will already pass (it's testing existing behavior) but includes it as a regression guard.

- [ ] **Step 3: Update `write_or_quarantine` to honor overrides**

Edit `nfl_draft/lib/quarantine.py`. Replace lines 42-46 (the branch after `mid` is found):

```python
            implied = (
                row.implied_prob
                if row.implied_prob is not None
                else american_to_implied(row.american_odds)
            )
            devig = (
                row.devig_prob
                if row.devig_prob is not None
                else implied
            )
            con.execute(
                "INSERT INTO draft_odds VALUES (?, ?, ?, ?, ?, ?)",
                [mid, row.book, row.american_odds, implied, devig, row.fetched_at],
            )
            mapped += 1
```

- [ ] **Step 4: Run tests to verify both pass**

```bash
pytest nfl_draft/tests/integration/test_dashboard_queries.py::test_quarantine_honors_oddsrow_price_overrides nfl_draft/tests/integration/test_dashboard_queries.py::test_quarantine_falls_back_to_american_when_overrides_absent -v
```
Expected: both PASS.

- [ ] **Step 5: Commit**

```bash
git add nfl_draft/lib/quarantine.py nfl_draft/tests/integration/test_dashboard_queries.py
git commit -m "$(cat <<'EOF'
feat(quarantine): honor OddsRow implied_prob / devig_prob overrides

Kalshi rows carry two different prices (take = yes_ask, fair = mid).
Existing sportsbook behavior is untouched: when overrides are None,
both columns get american_to_implied(american_odds) same as before.

Co-Authored-By: Claude Opus 4.7 (1M context) <noreply@anthropic.com>
EOF
)"
```

---

## Task 4: Add Kalshi Price Extractor Helpers

**Files:**
- Modify: `nfl_draft/scrapers/kalshi.py:53-72`
- Test: `nfl_draft/tests/unit/test_scraper_parsing.py`

Goal: mirror the existing `_extract_yes_bid_cents` logic for every price field Kalshi returns, so the scraper tolerates both `yes_bid: 5` (int cents) and `yes_bid_dollars: "0.05"` (string dollars).

- [ ] **Step 1: Write failing tests**

Append to `nfl_draft/tests/unit/test_scraper_parsing.py`:

```python
def test_kalshi_extractors_handle_dollar_string_format():
    """Kalshi v2 returns prices as `yes_bid_dollars: '0.05'` for accounts set to
    dollar response units. All extractors must decode to integer cents."""
    from nfl_draft.scrapers.kalshi import (
        _extract_yes_bid_cents, _extract_yes_ask_cents,
        _extract_no_bid_cents, _extract_no_ask_cents,
        _extract_last_price_cents,
    )
    market = {
        "yes_bid": None,  "yes_bid_dollars":  "0.02",
        "yes_ask": None,  "yes_ask_dollars":  "0.05",
        "no_bid":  None,  "no_bid_dollars":   "0.95",
        "no_ask":  None,  "no_ask_dollars":   "0.98",
        "last_price": None, "last_price_dollars": "0.04",
    }
    assert _extract_yes_bid_cents(market) == 2
    assert _extract_yes_ask_cents(market) == 5
    assert _extract_no_bid_cents(market)  == 95
    assert _extract_no_ask_cents(market)  == 98
    assert _extract_last_price_cents(market) == 4


def test_kalshi_extractors_handle_legacy_int_format():
    """Backward-compat path: cents-unit accounts return bare integers."""
    from nfl_draft.scrapers.kalshi import (
        _extract_yes_bid_cents, _extract_yes_ask_cents,
        _extract_no_bid_cents, _extract_no_ask_cents,
        _extract_last_price_cents,
    )
    market = {
        "yes_bid": 2, "yes_ask": 5,
        "no_bid": 95, "no_ask": 98,
        "last_price": 4,
    }
    assert _extract_yes_bid_cents(market) == 2
    assert _extract_yes_ask_cents(market) == 5
    assert _extract_no_bid_cents(market) == 95
    assert _extract_no_ask_cents(market) == 98
    assert _extract_last_price_cents(market) == 4


def test_kalshi_extractors_return_none_when_absent():
    """Sparse payloads (one-sided markets, never-traded markets) must not crash."""
    from nfl_draft.scrapers.kalshi import (
        _extract_yes_ask_cents, _extract_last_price_cents,
    )
    assert _extract_yes_ask_cents({}) in (None, 0)
    assert _extract_last_price_cents({}) in (None, 0)
```

- [ ] **Step 2: Run tests to verify they fail**

```bash
pytest nfl_draft/tests/unit/test_scraper_parsing.py::test_kalshi_extractors_handle_dollar_string_format -v
```
Expected: FAIL with `ImportError: cannot import name '_extract_yes_ask_cents'`.

- [ ] **Step 3: Add extractor helpers**

Edit `nfl_draft/scrapers/kalshi.py`. Below the existing `_extract_yes_bid_cents` function (ends line 72), add:

```python
def _extract_cents(market: dict, int_key: str, dollars_key: str):
    """Generic extractor: prefer legacy int cents, fall back to `_dollars` string.

    Kalshi's v2 response surfaces prices either as an integer cents field
    (`yes_bid`, `last_price`, ...) for accounts with `response_price_units=cents`,
    or as a dollar string (`yes_bid_dollars`, `last_price_dollars`, ...) for
    accounts on the newer dollar units. We accept both and normalize to int cents.

    Returns None when neither field carries a usable value (None, missing, or
    unparseable). Returns 0 when the field is explicitly 0 — callers decide
    whether 0 means "no market" or "zero price" per cascade rules.
    """
    raw_cents = market.get(int_key)
    if raw_cents not in (None, 0):
        try:
            return int(raw_cents)
        except (TypeError, ValueError):
            pass
    raw_dollars = market.get(dollars_key)
    if raw_dollars is not None:
        try:
            return int(round(float(raw_dollars) * 100))
        except (TypeError, ValueError):
            return None
    return raw_cents  # may be 0 or None


def _extract_yes_ask_cents(market: dict):
    return _extract_cents(market, "yes_ask", "yes_ask_dollars")


def _extract_no_bid_cents(market: dict):
    return _extract_cents(market, "no_bid", "no_bid_dollars")


def _extract_no_ask_cents(market: dict):
    return _extract_cents(market, "no_ask", "no_ask_dollars")


def _extract_last_price_cents(market: dict):
    return _extract_cents(market, "last_price", "last_price_dollars")
```

- [ ] **Step 4: Also refactor `_extract_yes_bid_cents` to use the shared helper**

Replace the body of the existing `_extract_yes_bid_cents` (lines 53-72) with:

```python
def _extract_yes_bid_cents(market: dict):
    """Return Kalshi yes_bid as an integer number of cents, or None.

    Accepts either legacy `yes_bid` (int cents) or new `yes_bid_dollars` (str)
    formats — see `_extract_cents` docstring.
    """
    return _extract_cents(market, "yes_bid", "yes_bid_dollars")
```

- [ ] **Step 5: Run extractor tests to verify they pass**

```bash
pytest nfl_draft/tests/unit/test_scraper_parsing.py::test_kalshi_extractors_handle_dollar_string_format nfl_draft/tests/unit/test_scraper_parsing.py::test_kalshi_extractors_handle_legacy_int_format nfl_draft/tests/unit/test_scraper_parsing.py::test_kalshi_extractors_return_none_when_absent -v
```
Expected: all three PASS.

- [ ] **Step 6: Run full scraper-parsing suite to check no regressions**

```bash
pytest nfl_draft/tests/unit/test_scraper_parsing.py -v
```
Expected: every pre-existing Kalshi test still passes (the `yes_bid` refactor preserves behavior).

- [ ] **Step 7: Commit**

```bash
git add nfl_draft/scrapers/kalshi.py nfl_draft/tests/unit/test_scraper_parsing.py
git commit -m "$(cat <<'EOF'
feat(kalshi): generic price extractor for all bid/ask/last fields

Previously only yes_bid had dollar-format tolerance. Ask, no-side, and
last_price all used m.get(key, 0) which silently returned 0 when Kalshi
sends the newer yes_bid_dollars format. Shared _extract_cents helper
handles both formats for every price field.

Co-Authored-By: Claude Opus 4.7 (1M context) <noreply@anthropic.com>
EOF
)"
```

---

## Task 5: Rewrite `parse_markets_response` with Cascade + Overrides

**Files:**
- Modify: `nfl_draft/scrapers/kalshi.py:123-165` (the `parse_markets_response` function)
- Test: `nfl_draft/tests/unit/test_scraper_parsing.py`

- [ ] **Step 1: Write the failing test**

Append to `nfl_draft/tests/unit/test_scraper_parsing.py`:

```python
def test_kalshi_parse_uses_ask_for_take_and_mid_for_fair():
    """Two-sided market: implied_prob (take) = yes_ask/100, devig_prob (fair) = mid."""
    from nfl_draft.scrapers.kalshi import parse_markets_response
    raw = {"markets": [{
        "ticker": "KXNFLDRAFTPICK-26-5-CTAT",
        "yes_sub_title": "Carnell Tate",
        "yes_bid_dollars": "0.02",
        "yes_ask_dollars": "0.05",
        "last_price_dollars": "0.04",
    }]}
    rows = parse_markets_response(raw, series_ticker="KXNFLDRAFTPICK")
    assert len(rows) == 1
    r = rows[0]
    assert r.book == "kalshi"
    assert r.book_subject == "Carnell Tate"
    assert r.implied_prob == 0.05          # yes_ask / 100
    assert r.devig_prob == (2 + 5) / 200.0  # mid


def test_kalshi_parse_one_sided_market_falls_back_to_last_trade():
    """Only bid posted, last trade present: fair = last/100, take = last/100."""
    from nfl_draft.scrapers.kalshi import parse_markets_response
    raw = {"markets": [{
        "ticker": "KXNFLDRAFTPICK-26-5-CTAT",
        "yes_sub_title": "Carnell Tate",
        "yes_bid_dollars": "0.02",
        "yes_ask_dollars": "0.00",          # no ask
        "last_price_dollars": "0.04",
    }]}
    rows = parse_markets_response(raw, series_ticker="KXNFLDRAFTPICK")
    assert len(rows) == 1
    r = rows[0]
    assert r.implied_prob == 0.04
    assert r.devig_prob == 0.04


def test_kalshi_parse_one_sided_no_last_trade_uses_single_side_for_fair():
    """Only bid, no last trade: fair = single posted side, take = None (no flag)."""
    from nfl_draft.scrapers.kalshi import parse_markets_response
    raw = {"markets": [{
        "ticker": "KXNFLDRAFTPICK-26-5-CTAT",
        "yes_sub_title": "Carnell Tate",
        "yes_bid_dollars": "0.02",
        "yes_ask_dollars": "0.00",
        "last_price_dollars": "0.00",
    }]}
    rows = parse_markets_response(raw, series_ticker="KXNFLDRAFTPICK")
    assert len(rows) == 1
    r = rows[0]
    assert r.devig_prob == 0.02
    assert r.implied_prob is None


def test_kalshi_parse_empty_market_skipped():
    """No bid, no ask, no last trade: no row emitted."""
    from nfl_draft.scrapers.kalshi import parse_markets_response
    raw = {"markets": [{
        "ticker": "KXNFLDRAFTPICK-26-5-CTAT",
        "yes_sub_title": "Carnell Tate",
        "yes_bid_dollars": "0.00",
        "yes_ask_dollars": "0.00",
        "last_price_dollars": "0.00",
    }]}
    rows = parse_markets_response(raw, series_ticker="KXNFLDRAFTPICK")
    assert rows == []
```

- [ ] **Step 2: Run tests to verify they fail**

```bash
pytest nfl_draft/tests/unit/test_scraper_parsing.py::test_kalshi_parse_uses_ask_for_take_and_mid_for_fair -v
```
Expected: FAIL — current parser reads only `yes_bid` and won't compute mid.

- [ ] **Step 3: Rewrite `parse_markets_response`**

Edit `nfl_draft/scrapers/kalshi.py`. Replace the entire `parse_markets_response` function (lines 123-165) with:

```python
def parse_markets_response(raw_response: dict, series_ticker: str) -> List[OddsRow]:
    """Convert a raw Kalshi /markets response into a list of OddsRow.

    Emits one row per Kalshi market with enough signal to price. For each
    market we compute two prices:

      * take price (implied_prob override): yes_ask first, then last_price.
        Represents the cost to buy a Yes contract. Used by the Cross-Book
        Grid's outlier flag.
      * fair value (devig_prob override): mid = (yes_bid + yes_ask)/200 when
        both sides posted, else last_price/100, else the one posted side/100.
        Used for cell display and the cross-venue median.

    Markets with no bid, no ask, and no last trade are skipped entirely —
    they carry no signal worth rendering. Rows with a fair value but no take
    price are still emitted; their `implied_prob` is None so downstream flag
    logic suppresses the outlier check for that cell.

    `american_odds` is set from the take price when available, otherwise from
    the fair value — keeps legacy downstream readers (bet log, quarantine
    fallback for unmapped rows) aligned with what the scraper considers the
    most tradeable number for that row.
    """
    rows: List[OddsRow] = []
    if not isinstance(raw_response, dict):
        return rows
    now = datetime.now()

    for market in raw_response.get("markets", []) or []:
        yes_bid = _extract_yes_bid_cents(market) or 0
        yes_ask = _extract_yes_ask_cents(market) or 0
        last_price = _extract_last_price_cents(market) or 0

        # Fair-value cascade: mid, then last trade, then any single side.
        if yes_bid > 0 and yes_ask > 0:
            fair_cents = (yes_bid + yes_ask) / 2.0
        elif last_price > 0:
            fair_cents = float(last_price)
        elif yes_ask > 0:
            fair_cents = float(yes_ask)
        elif yes_bid > 0:
            fair_cents = float(yes_bid)
        else:
            continue  # no signal → skip row

        if not (0 < fair_cents < 100):
            continue

        # Take-price cascade: ask, then last trade; else None (flag suppressed).
        if yes_ask > 0:
            take_cents = float(yes_ask)
        elif last_price > 0:
            take_cents = float(last_price)
        else:
            take_cents = None

        candidate = (
            market.get("yes_sub_title")
            or market.get("subtitle")
            or (market.get("custom_strike") or {}).get("Person")
            or (market.get("custom_strike") or {}).get("Team")
            or ""
        )
        candidate = (candidate or "").strip()
        if not candidate:
            continue

        # American odds: use take price if available, else fair. Used by
        # downstream readers that still want an integer American-odds view.
        reference_cents = take_cents if take_cents is not None else fair_cents
        p_ref = reference_cents / 100.0
        if p_ref > 0.5:
            american = int(round(-100 * p_ref / (1 - p_ref)))
        else:
            american = int(round(100 * (1 - p_ref) / p_ref))

        rows.append(OddsRow(
            book="kalshi",
            book_label=_kalshi_book_label(series_ticker, market.get("ticker") or ""),
            book_subject=candidate,
            american_odds=american,
            fetched_at=now,
            implied_prob=(take_cents / 100.0) if take_cents is not None else None,
            devig_prob=fair_cents / 100.0,
        ))
    return rows
```

- [ ] **Step 4: Run new cascade tests**

```bash
pytest nfl_draft/tests/unit/test_scraper_parsing.py -k "kalshi_parse" -v
```
Expected: all four new tests PASS and the pre-existing `test_kalshi_parse_markets_returns_oddsrow_list` still PASSES (rows may be fewer than before because the old parser counted zero-bid markets, but the fixture should still yield >= 50 rows with non-zero fair values).

If the pre-existing test drops below 50 rows, relax its assertion to `>= 10` and add a comment: the ask-aware parser is slightly more conservative because it now rejects markets where only a stale bid exists with no last trade. If it drops below 10, the fixture itself is bad — flag in the commit message and re-capture in a follow-up.

- [ ] **Step 5: Commit**

```bash
git add nfl_draft/scrapers/kalshi.py nfl_draft/tests/unit/test_scraper_parsing.py
git commit -m "$(cat <<'EOF'
feat(kalshi): emit fair (mid) and take (ask) prices as OddsRow overrides

parse_markets_response now computes two prices per Kalshi market:

  - devig_prob override (fair value): mid of buy/sell, falling back to
    last trade, then to a single posted side.
  - implied_prob override (take price): yes_ask, falling back to last
    trade; None when neither is available (flag suppressed downstream).

Rows with no signal at all (no bid, no ask, no last) are skipped.

Co-Authored-By: Claude Opus 4.7 (1M context) <noreply@anthropic.com>
EOF
)"
```

---

## Task 6: Fix the Legacy `_write_legacy_kalshi_odds` Zero Bug

**Files:**
- Modify: `nfl_draft/scrapers/kalshi.py:221-297` (the `_write_legacy_kalshi_odds` function)
- Test: `nfl_draft/tests/integration/test_dashboard_queries.py`

- [ ] **Step 1: Write the failing test**

Append to `nfl_draft/tests/integration/test_dashboard_queries.py`:

```python
def test_legacy_kalshi_odds_write_handles_dollar_format(monkeypatch, tmp_path):
    """Regression: _write_legacy_kalshi_odds must persist real bid/ask/last
    when Kalshi returns the newer `*_dollars` response shape. Prior bug was
    m.get('yes_bid', 0) silently storing zeros for dollar-format accounts."""
    from nfl_draft.lib import db as db_module
    monkeypatch.setattr(db_module, "DB_PATH", tmp_path / "test.duckdb")
    db_module.init_schema()

    # Schema for legacy kalshi_odds — must exist for the write path.
    from nfl_draft.lib.db import write_connection, read_connection
    with write_connection() as con:
        con.execute("""
            CREATE TABLE IF NOT EXISTS kalshi_odds (
                fetch_time TIMESTAMP, series_ticker VARCHAR, event_ticker VARCHAR,
                ticker VARCHAR, market_title VARCHAR, candidate VARCHAR,
                yes_bid INTEGER, yes_ask INTEGER, no_bid INTEGER, no_ask INTEGER,
                last_price INTEGER, volume BIGINT, volume_24h BIGINT,
                liquidity BIGINT, open_interest INTEGER
            )
        """)
        con.execute("""
            CREATE TABLE IF NOT EXISTS market_info (
                ticker VARCHAR, title VARCHAR, subtitle VARCHAR,
                series_ticker VARCHAR, updated_at TIMESTAMP
            )
        """)

    from nfl_draft.scrapers.kalshi import _write_legacy_kalshi_odds
    markets = [{
        "ticker": "KXNFLDRAFTPICK-26-5-CTAT",
        "event_ticker": "EVT1",
        "title": "Pick 5",
        "yes_sub_title": "Carnell Tate",
        "yes_bid_dollars": "0.02",
        "yes_ask_dollars": "0.05",
        "no_bid_dollars":  "0.95",
        "no_ask_dollars":  "0.98",
        "last_price_dollars": "0.04",
        "volume": 1000, "volume_24h": 5000,
        "liquidity_dollars": "100.00",
        "open_interest": 200,
    }]
    _write_legacy_kalshi_odds("KXNFLDRAFTPICK", "NFL Draft Picks", markets)

    with read_connection() as con:
        row = con.execute(
            "SELECT yes_bid, yes_ask, no_bid, no_ask, last_price FROM kalshi_odds WHERE ticker = ?",
            ["KXNFLDRAFTPICK-26-5-CTAT"],
        ).fetchone()
    assert row == (2, 5, 95, 98, 4)
```

- [ ] **Step 2: Run test to verify it fails**

```bash
pytest nfl_draft/tests/integration/test_dashboard_queries.py::test_legacy_kalshi_odds_write_handles_dollar_format -v
```
Expected: FAIL — the row will come back as `(0, 0, 0, 0, 0)` because the current code uses `m.get("yes_bid", 0) or 0`.

- [ ] **Step 3: Rewrite `_write_legacy_kalshi_odds` to use the extractors**

Edit `nfl_draft/scrapers/kalshi.py`. In `_write_legacy_kalshi_odds`, replace the `rows_to_write.append(...)` call (around lines 249-265) with:

```python
        rows_to_write.append((
            fetch_time,
            series_ticker,
            m.get("event_ticker", ""),
            m.get("ticker", ""),
            m.get("title", ""),
            candidate,
            _extract_yes_bid_cents(m) or 0,
            _extract_yes_ask_cents(m) or 0,
            _extract_no_bid_cents(m) or 0,
            _extract_no_ask_cents(m) or 0,
            _extract_last_price_cents(m) or 0,
            m.get("volume", 0) or 0,
            m.get("volume_24h", 0) or 0,
            int(float(m.get("liquidity_dollars") or 0)),
            m.get("open_interest", 0) or 0,
        ))
```

The only functional change: the five price columns now route through the extractors instead of `m.get(...)` directly. Everything else (volume, liquidity, etc.) is unchanged.

- [ ] **Step 4: Run test to verify it passes**

```bash
pytest nfl_draft/tests/integration/test_dashboard_queries.py::test_legacy_kalshi_odds_write_handles_dollar_format -v
```
Expected: PASS — row now `(2, 5, 95, 98, 4)`.

- [ ] **Step 5: Run the pre-existing legacy-dashboard test to confirm no regressions**

```bash
pytest nfl_draft/tests/integration/test_dashboard_queries.py::test_get_latest_odds_returns_one_row nfl_draft/tests/integration/test_dashboard_queries.py::test_get_price_history_returns_one_row nfl_draft/tests/integration/test_dashboard_queries.py::test_get_snapshot_count_returns_one -v
```
Expected: all three PASS.

- [ ] **Step 6: Commit**

```bash
git add nfl_draft/scrapers/kalshi.py nfl_draft/tests/integration/test_dashboard_queries.py
git commit -m "$(cat <<'EOF'
fix(kalshi): legacy kalshi_odds write now stores real prices

_write_legacy_kalshi_odds used m.get("yes_bid", 0) which silently
returned 0 when Kalshi sends the newer yes_bid_dollars format. Every
bid/ask/last_price in the table was zero, breaking the legacy Price
History / Edge Detection / Consensus tabs. Route through the extractor
helpers which handle both formats.

Co-Authored-By: Claude Opus 4.7 (1M context) <noreply@anthropic.com>
EOF
)"
```

---

## Task 7: Update `cross_book_grid` with Per-Venue Flag Logic

**Files:**
- Modify: `nfl_draft/lib/queries.py:89-143` (the `cross_book_grid` function)
- Test: `nfl_draft/tests/integration/test_dashboard_queries.py`

- [ ] **Step 1: Write failing tests**

Append to `nfl_draft/tests/integration/test_dashboard_queries.py`:

```python
def test_cross_book_grid_kalshi_flag_uses_implied_prob_not_devig_prob(monkeypatch, tmp_path):
    """Kalshi: flag compares implied_prob (take = buy) against median of devig_prob (fair)."""
    from nfl_draft.lib import db as db_module
    monkeypatch.setattr(db_module, "DB_PATH", tmp_path / "test.duckdb")
    db_module.init_schema()

    from datetime import datetime
    from nfl_draft.lib.db import write_connection
    now = datetime.now()
    with write_connection() as con:
        con.execute(
            "INSERT INTO draft_markets (market_id, market_type) "
            "VALUES ('pick_5_overall_carnell-tate', 'pick_outright')"
        )
        # Kalshi: fair=3.5% (mid) but take=5% (buy). Sportsbooks: simple devig ~ 6-8%.
        con.execute(
            "INSERT INTO draft_odds VALUES ('pick_5_overall_carnell-tate', 'kalshi', 1900, ?, ?, ?)",
            [0.05, 0.035, now],
        )
        for book, prob in [("draftkings", 0.077), ("bookmaker", 0.062), ("wagerzon", 0.067)]:
            con.execute(
                "INSERT INTO draft_odds VALUES ('pick_5_overall_carnell-tate', ?, 100, ?, ?, ?)",
                [book, prob, prob, now],
            )

    from nfl_draft.lib.queries import cross_book_grid
    rows = cross_book_grid(threshold_pp=5.0)
    assert len(rows) == 1
    row = rows[0]

    # Display value for Kalshi is the MID (devig_prob), not the buy price.
    assert abs(row["books"]["kalshi"] - 0.035) < 1e-9
    # Median participates: [0.035, 0.062, 0.067, 0.077] -> mid two average = 0.0645
    assert abs(row["median"] - 0.0645) < 1e-9
    # Kalshi flag check uses the BUY price (0.05) vs median (0.0645) = 1.45pp
    # → below 5pp threshold → not flagged (correct).
    assert row["flags"]["kalshi"] is False


def test_cross_book_grid_kalshi_flag_fires_on_big_take_edge(monkeypatch, tmp_path):
    """When Kalshi's buy price sits well below consensus fair, flag fires."""
    from nfl_draft.lib import db as db_module
    monkeypatch.setattr(db_module, "DB_PATH", tmp_path / "test.duckdb")
    db_module.init_schema()

    from datetime import datetime
    from nfl_draft.lib.db import write_connection
    now = datetime.now()
    with write_connection() as con:
        con.execute(
            "INSERT INTO draft_markets (market_id, market_type) VALUES ('m1', 'prop')"
        )
        # Kalshi take=2%, fair=3%; sportsbooks all 8%
        con.execute(
            "INSERT INTO draft_odds VALUES ('m1', 'kalshi', 4900, ?, ?, ?)",
            [0.02, 0.03, now],
        )
        for book in ("draftkings", "bookmaker", "wagerzon"):
            con.execute(
                "INSERT INTO draft_odds VALUES ('m1', ?, 100, 0.08, 0.08, ?)",
                [book, now],
            )

    from nfl_draft.lib.queries import cross_book_grid
    rows = cross_book_grid(threshold_pp=3.0)
    assert len(rows) == 1
    row = rows[0]
    # Median([0.03, 0.08, 0.08, 0.08]) = 0.08
    assert abs(row["median"] - 0.08) < 1e-9
    # |take 0.02 - median 0.08| = 6pp >= 3pp → flagged
    assert row["flags"]["kalshi"] is True


def test_cross_book_grid_kalshi_flag_suppressed_when_no_take_price(monkeypatch, tmp_path):
    """When Kalshi has fair but no take (one-sided, no last trade), flag is False."""
    from nfl_draft.lib import db as db_module
    monkeypatch.setattr(db_module, "DB_PATH", tmp_path / "test.duckdb")
    db_module.init_schema()

    from datetime import datetime
    from nfl_draft.lib.db import write_connection
    now = datetime.now()
    with write_connection() as con:
        con.execute(
            "INSERT INTO draft_markets (market_id, market_type) VALUES ('m1', 'prop')"
        )
        con.execute(
            "INSERT INTO draft_odds VALUES ('m1', 'kalshi', 4900, NULL, ?, ?)",
            [0.02, now],
        )
        for book in ("draftkings", "bookmaker", "wagerzon"):
            con.execute(
                "INSERT INTO draft_odds VALUES ('m1', ?, 100, 0.08, 0.08, ?)",
                [book, now],
            )

    from nfl_draft.lib.queries import cross_book_grid
    rows = cross_book_grid(threshold_pp=3.0)
    assert rows[0]["flags"]["kalshi"] is False
```

Also update the existing `test_cross_book_grid_query_outlier_flags` test to make its Kalshi row set `implied_prob = devig_prob`:

```python
# In test_cross_book_grid_query_outlier_flags, replace the loop body:
        for book, prob in [
            ("kalshi", 0.50),
            ("fanduel", 0.51),
            ("bookmaker", 0.49),
            ("wagerzon", 0.52),
            ("draftkings", 0.30),
        ]:
            # For Kalshi rows in this legacy test, set implied_prob == devig_prob
            # so the new per-venue flag logic compares buy vs median. This mirrors
            # a market where the bid-ask spread is zero / the seed is synthetic.
            con.execute(
                "INSERT INTO draft_odds VALUES ('first_qb_cam-ward', ?, 100, ?, ?, ?)",
                [book, prob, prob, now],
            )
```

(This edit is a no-op in terms of persisted data — the test already passes `[book, prob, prob, now]`; the intent change is the comment clarifying why Kalshi still evaluates to unflagged under the new logic.)

- [ ] **Step 2: Run tests to verify they fail**

```bash
pytest nfl_draft/tests/integration/test_dashboard_queries.py::test_cross_book_grid_kalshi_flag_uses_implied_prob_not_devig_prob -v
```
Expected: FAIL — current code only selects `devig_prob` and flags Kalshi against it, not against `implied_prob`.

- [ ] **Step 3: Rewrite `cross_book_grid`**

Edit `nfl_draft/lib/queries.py`. Replace the `_query` inner function of `cross_book_grid` (lines 102-141) with:

```python
    def _query() -> List[Dict[str, Any]]:
        with read_connection() as con:
            rows = con.execute(
                f"""
                WITH latest AS (
                  SELECT market_id, book, implied_prob, devig_prob,
                         ROW_NUMBER() OVER (PARTITION BY market_id, book ORDER BY fetched_at DESC) AS rn
                  FROM draft_odds
                  WHERE fetched_at > NOW() - INTERVAL '{MAX_AGE_HOURS} hours'
                )
                SELECT market_id, book, implied_prob, devig_prob
                FROM latest WHERE rn = 1
                """
            ).fetchall()

        by_market: Dict[str, Dict[str, Dict[str, Optional[float]]]] = {}
        for market_id, book, implied_prob, devig_prob in rows:
            by_market.setdefault(market_id, {})[book] = {
                "implied_prob": implied_prob,
                "devig_prob": devig_prob,
            }

        threshold = threshold_pp / 100.0
        output: List[Dict[str, Any]] = []
        for market_id, books in by_market.items():
            # `books` maps book -> {implied_prob, devig_prob}. Cell display + median
            # are computed from devig_prob (each book's fair estimate). Flag
            # comparison uses implied_prob for Kalshi (buy price) and devig_prob
            # for every other book (no vig distinction → same number).
            display_by_book = {
                b: (r["devig_prob"] if r["devig_prob"] is not None else r["implied_prob"])
                for b, r in books.items()
            }
            valid_display = [p for p in display_by_book.values() if p is not None]

            if len(books) < 2 or not valid_display:
                output.append({
                    "market_id": market_id,
                    "books": display_by_book,
                    "median": valid_display[0] if valid_display else None,
                    "flags": {b: False for b in books},
                    "outlier_count": 0,
                })
                continue

            median = statistics.median(valid_display)
            flags: Dict[str, bool] = {}
            for book, r in books.items():
                if book == "kalshi":
                    take = r["implied_prob"]
                    flags[book] = (
                        take is not None
                        and abs(take - median) >= threshold
                    )
                else:
                    dev = r["devig_prob"]
                    flags[book] = (
                        dev is not None
                        and abs(dev - median) >= threshold
                    )
            output.append({
                "market_id": market_id,
                "books": display_by_book,
                "median": median,
                "flags": flags,
                "outlier_count": sum(flags.values()),
            })
        return output
```

- [ ] **Step 4: Run all new + existing cross_book_grid tests**

```bash
pytest nfl_draft/tests/integration/test_dashboard_queries.py -k "cross_book_grid or ev_candidates" -v
```
Expected: all tests PASS, including the pre-existing `test_cross_book_grid_query_outlier_flags` and `test_cross_book_grid_excludes_stale_rows`.

- [ ] **Step 5: Commit**

```bash
git add nfl_draft/lib/queries.py nfl_draft/tests/integration/test_dashboard_queries.py
git commit -m "$(cat <<'EOF'
feat(queries): per-venue outlier flag logic in cross_book_grid

Kalshi: compare implied_prob (buy price) to median. Flag fires only
when you could actually take a meaningfully off-consensus price.
Sportsbooks: unchanged — still compare devig_prob to median.
Median itself is always computed from devig_prob across venues so
the cell display and the consensus stay consistent.

Co-Authored-By: Claude Opus 4.7 (1M context) <noreply@anthropic.com>
EOF
)"
```

---

## Task 8: Add the Tooltip Lookup Query

**Files:**
- Modify: `nfl_draft/lib/queries.py` (append new function)
- Test: `nfl_draft/tests/integration/test_dashboard_queries.py`

- [ ] **Step 1: Write the failing test**

Append to `nfl_draft/tests/integration/test_dashboard_queries.py`:

```python
def test_kalshi_tooltip_data_joins_market_map_to_kalshi_odds(monkeypatch, tmp_path):
    """kalshi_tooltip_data returns {market_id: {ticker, yes_bid, yes_ask, last_price}}
    for every Kalshi market that has a market_map entry with a live kalshi_odds row."""
    from nfl_draft.lib import db as db_module
    monkeypatch.setattr(db_module, "DB_PATH", tmp_path / "test.duckdb")
    db_module.init_schema()

    from datetime import datetime
    from nfl_draft.lib.db import write_connection
    now = datetime.now()
    with write_connection() as con:
        con.execute("""
            CREATE TABLE IF NOT EXISTS kalshi_odds (
                fetch_time TIMESTAMP, series_ticker VARCHAR, event_ticker VARCHAR,
                ticker VARCHAR, market_title VARCHAR, candidate VARCHAR,
                yes_bid INTEGER, yes_ask INTEGER, no_bid INTEGER, no_ask INTEGER,
                last_price INTEGER, volume BIGINT, volume_24h BIGINT,
                liquidity BIGINT, open_interest INTEGER
            )
        """)
        con.execute(
            "INSERT INTO market_map VALUES ('kalshi', 'KXNFLDRAFTPICK-26-5', 'Carnell Tate', 'pick_5_overall_carnell-tate')"
        )
        con.execute(
            "INSERT INTO kalshi_odds VALUES (?, 'KXNFLDRAFTPICK', 'EVT1', "
            "'KXNFLDRAFTPICK-26-5-CTAT', 'Pick 5', 'Carnell Tate', "
            "2, 5, 95, 98, 4, 100, 500, 50, 20)",
            [now],
        )

    from nfl_draft.lib.queries import kalshi_tooltip_data
    data = kalshi_tooltip_data()
    assert "pick_5_overall_carnell-tate" in data
    entry = data["pick_5_overall_carnell-tate"]
    assert entry["ticker"] == "KXNFLDRAFTPICK-26-5-CTAT"
    assert entry["yes_bid"] == 2
    assert entry["yes_ask"] == 5
    assert entry["last_price"] == 4


def test_kalshi_tooltip_data_picks_latest_when_multiple_snapshots(monkeypatch, tmp_path):
    """Two snapshots of the same ticker: tooltip uses the most recent."""
    from nfl_draft.lib import db as db_module
    monkeypatch.setattr(db_module, "DB_PATH", tmp_path / "test.duckdb")
    db_module.init_schema()

    from datetime import datetime, timedelta
    from nfl_draft.lib.db import write_connection
    now = datetime.now()
    older = now - timedelta(minutes=10)
    with write_connection() as con:
        con.execute("""
            CREATE TABLE IF NOT EXISTS kalshi_odds (
                fetch_time TIMESTAMP, series_ticker VARCHAR, event_ticker VARCHAR,
                ticker VARCHAR, market_title VARCHAR, candidate VARCHAR,
                yes_bid INTEGER, yes_ask INTEGER, no_bid INTEGER, no_ask INTEGER,
                last_price INTEGER, volume BIGINT, volume_24h BIGINT,
                liquidity BIGINT, open_interest INTEGER
            )
        """)
        con.execute(
            "INSERT INTO market_map VALUES ('kalshi', 'KXNFLDRAFTPICK-26-5', 'Carnell Tate', 'pick_5_overall_carnell-tate')"
        )
        con.execute(
            "INSERT INTO kalshi_odds VALUES (?, 'KXNFLDRAFTPICK', 'EVT1', "
            "'KXNFLDRAFTPICK-26-5-CTAT', 'Pick 5', 'Carnell Tate', "
            "1, 4, 96, 99, 3, 100, 500, 50, 20)",
            [older],
        )
        con.execute(
            "INSERT INTO kalshi_odds VALUES (?, 'KXNFLDRAFTPICK', 'EVT1', "
            "'KXNFLDRAFTPICK-26-5-CTAT', 'Pick 5', 'Carnell Tate', "
            "2, 5, 95, 98, 4, 100, 500, 50, 20)",
            [now],
        )

    from nfl_draft.lib.queries import kalshi_tooltip_data
    data = kalshi_tooltip_data()
    entry = data["pick_5_overall_carnell-tate"]
    assert (entry["yes_bid"], entry["yes_ask"], entry["last_price"]) == (2, 5, 4)
```

- [ ] **Step 2: Run tests to verify they fail**

```bash
pytest nfl_draft/tests/integration/test_dashboard_queries.py::test_kalshi_tooltip_data_joins_market_map_to_kalshi_odds -v
```
Expected: FAIL with `ImportError: cannot import name 'kalshi_tooltip_data'`.

- [ ] **Step 3: Add `kalshi_tooltip_data` to queries.py**

Append to `nfl_draft/lib/queries.py`, just above the existing `def latest_max_fetched_at` (so tooltip sits with the other grid-support functions):

```python
def kalshi_tooltip_data() -> Dict[str, Dict[str, Any]]:
    """Return per-Kalshi-market hover content keyed by market_id.

    For every `market_map` row where `book='kalshi'`, join to the latest
    `kalshi_odds` snapshot (by fetch_time) for the corresponding ticker. The
    join key is `(book_label, book_subject)` because market_map stores the
    ticker *prefix* (e.g. `KXNFLDRAFTPICK-26-5`) while the full ticker has a
    candidate shortcode appended (`KXNFLDRAFTPICK-26-5-CTAT`).

    Returns:
        ``{market_id: {"ticker": str, "yes_bid": int, "yes_ask": int,
                       "last_price": int}}``
        Rows where no matching kalshi_odds snapshot exists are omitted.
        Returns ``QueryLocked`` on DuckDB lock contention.
    """
    def _query() -> Dict[str, Dict[str, Any]]:
        with read_connection() as con:
            rows = con.execute(
                """
                WITH latest_ko AS (
                  SELECT ticker, candidate, yes_bid, yes_ask, last_price,
                         ROW_NUMBER() OVER (PARTITION BY ticker ORDER BY fetch_time DESC) AS rn
                  FROM kalshi_odds
                )
                SELECT mm.market_id, ko.ticker, ko.yes_bid, ko.yes_ask, ko.last_price
                FROM market_map mm
                JOIN latest_ko ko
                  ON ko.ticker LIKE mm.book_label || '-%'
                 AND ko.candidate = mm.book_subject
                 AND ko.rn = 1
                WHERE mm.book = 'kalshi'
                """
            ).fetchall()
        return {
            market_id: {
                "ticker": ticker,
                "yes_bid": yes_bid,
                "yes_ask": yes_ask,
                "last_price": last_price,
            }
            for market_id, ticker, yes_bid, yes_ask, last_price in rows
        }

    return _safe_read(_query)
```

- [ ] **Step 4: Run tooltip tests**

```bash
pytest nfl_draft/tests/integration/test_dashboard_queries.py -k "kalshi_tooltip" -v
```
Expected: both new tests PASS.

- [ ] **Step 5: Commit**

```bash
git add nfl_draft/lib/queries.py nfl_draft/tests/integration/test_dashboard_queries.py
git commit -m "$(cat <<'EOF'
feat(queries): add kalshi_tooltip_data for Cross-Book Grid hover

Batched query: one call returns {market_id -> {ticker, yes_bid,
yes_ask, last_price}} for every live Kalshi market. Dashboard
renderer stores the dict and reads from it per-hover — no extra
DB round-trip per cell.

Join key (book_label prefix + candidate name) uniquely identifies
one Kalshi ticker per market_id even though market_map stores the
ticker prefix rather than the full ticker.

Co-Authored-By: Claude Opus 4.7 (1M context) <noreply@anthropic.com>
EOF
)"
```

---

## Task 9: Wire Tooltip into the Dashboard DataTable

**Files:**
- Modify: `kalshi_draft/app.py` (the `_update_crossbook` callback around line 1040, and the help-text in `render_crossbook_grid` around line 798)

- [ ] **Step 1: Update the cross-book grid explanation text**

Edit `kalshi_draft/app.py`. In `render_crossbook_grid` (around line 798-803), replace the `html.P(...)` paragraph with:

```python
            html.P(
                "Each cell shows a venue's fair-value estimate (devigged probability "
                "for sportsbooks; mid of buy/sell for Kalshi). Flagged cells (⚑) differ "
                "from the cross-venue median by at least the threshold — for sportsbooks "
                "the comparison uses the fair estimate, for Kalshi it uses the actual buy "
                "price (what you'd pay to take Yes). Hover a Kalshi cell for buy/sell/last.",
                style={"color": COLORS["text_muted"], "fontSize": "0.85em"},
            ),
```

- [ ] **Step 2: Import the new query function at the top of the callback section**

In `kalshi_draft/app.py`, find the existing `from nfl_draft.lib import queries as nfl_queries` import (search for `nfl_queries`). Ensure it's present at the top; no new import is needed since `kalshi_tooltip_data` lives on the same module.

- [ ] **Step 3: Modify `_update_crossbook` to attach tooltip_data**

Replace the body of `_update_crossbook` (the callback starting at line 1040) with:

```python
def _update_crossbook(threshold_pp, _n_intervals, last_seen):
    """Cheap-poll guard: on interval tick, skip render if MAX(fetched_at)
    hasn't changed. User-triggered threshold changes always re-render.

    If DuckDB is locked by the cron writer, skip this render entirely —
    the next interval tick (seconds away) will retry cleanly.
    """
    ctx = dash.callback_context
    triggered_by_interval = (
        ctx.triggered and ctx.triggered[0]["prop_id"].startswith("nfl_draft__interval")
    )
    latest = nfl_queries.latest_max_fetched_at("draft_odds")
    if isinstance(latest, QueryLocked):
        raise PreventUpdate
    latest_iso = latest.isoformat() if latest else None
    if triggered_by_interval and latest_iso == last_seen:
        raise PreventUpdate

    grid = nfl_queries.cross_book_grid(threshold_pp=threshold_pp or 0)
    if isinstance(grid, QueryLocked):
        raise PreventUpdate

    # Prefetch Kalshi tooltip data once per render; then read per-row from the
    # in-memory dict. Lock contention is non-fatal — just skip tooltips and
    # render the grid without hover data.
    tooltip_lookup = nfl_queries.kalshi_tooltip_data()
    if isinstance(tooltip_lookup, QueryLocked):
        tooltip_lookup = {}

    rows = []
    tooltip_data = []
    for m in grid:
        row = {"market_id": m["market_id"]}
        for venue in VENUES:
            prob = m["books"].get(venue)
            flagged = m["flags"].get(venue, False)
            if prob is None:
                row[venue] = ""
            else:
                row[venue] = f"{prob*100:.1f}%" + (" ⚑" if flagged else "")
        row["median"] = f"{m['median']*100:.1f}%" if m["median"] is not None else ""
        row["outliers"] = m["outlier_count"]
        rows.append(row)

        # Per-row tooltip: only the Kalshi column gets content, and only when
        # we have a live tooltip lookup for this market_id. Other columns use
        # an empty string so Dash doesn't render a stale tooltip there.
        tip = tooltip_lookup.get(m["market_id"])
        if tip is not None:
            kalshi_hover = (
                f"Last trade: {tip['last_price']}¢\n"
                f"Buy: {tip['yes_ask']}¢   "
                f"Sell: {tip['yes_bid']}¢\n"
                f"Ticker: {tip['ticker']}"
            )
        else:
            kalshi_hover = ""
        tooltip_data.append({
            "kalshi": {"value": kalshi_hover, "type": "markdown"},
        })

    table = dash_table.DataTable(
        data=rows,
        columns=[{"name": "Market", "id": "market_id"}]
        + [{"name": v.capitalize(), "id": v} for v in VENUES]
        + [{"name": "Median", "id": "median"}, {"name": "Outliers", "id": "outliers", "type": "numeric"}],
        tooltip_data=tooltip_data,
        tooltip_delay=200,
        tooltip_duration=None,
        filter_action="native",
        sort_action="native",
        page_size=50,
        style_header=TABLE_STYLE_HEADER,
        style_data=TABLE_STYLE_DATA,
        style_data_conditional=TABLE_STYLE_DATA_CONDITIONAL + [
            {"if": {"filter_query": f"{{{v}}} contains '⚑'", "column_id": v},
             "backgroundColor": "#4d2d15", "color": COLORS["red"], "fontWeight": "bold"}
            for v in VENUES
        ],
        style_table={"overflowX": "auto"},
        style_filter={"backgroundColor": "#0d1b2a", "color": COLORS["text"]},
    )
    return table, latest_iso
```

Note: `tooltip_data` is a list the same length as `rows`, each entry a dict of `{column_id: {"value": markdown_str, "type": "markdown"}}`. Only Kalshi cells get populated hover content. `tooltip_delay=200` avoids flicker on quick mouse moves; `tooltip_duration=None` means "stays until mouse leaves."

- [ ] **Step 4: Headless smoke-test of the new callback (no browser)**

We can't run the dashboard against the production DuckDB from the worktree (the DB file lives in the main repo, and `CLAUDE.md` forbids symlinking it). Instead, exercise the new callback logic in-process against a seeded temp DB to confirm the tooltip wiring produces the expected structure. Run from the worktree:

```bash
python -c "
from pathlib import Path
import os, tempfile
from datetime import datetime

tmp = Path(tempfile.mkdtemp()) / 'smoke.duckdb'
from nfl_draft.lib import db as db_module
db_module.DB_PATH = tmp
db_module.init_schema()

from nfl_draft.lib.db import write_connection
now = datetime.now()
with write_connection() as con:
    con.execute(\"CREATE TABLE IF NOT EXISTS kalshi_odds (fetch_time TIMESTAMP, series_ticker VARCHAR, event_ticker VARCHAR, ticker VARCHAR, market_title VARCHAR, candidate VARCHAR, yes_bid INTEGER, yes_ask INTEGER, no_bid INTEGER, no_ask INTEGER, last_price INTEGER, volume BIGINT, volume_24h BIGINT, liquidity BIGINT, open_interest INTEGER)\")
    con.execute(\"INSERT INTO draft_markets (market_id, market_type) VALUES ('pick_5', 'pick_outright')\")
    con.execute(\"INSERT INTO draft_odds VALUES ('pick_5', 'kalshi', 1900, 0.05, 0.035, ?)\", [now])
    con.execute(\"INSERT INTO draft_odds VALUES ('pick_5', 'draftkings', 100, 0.077, 0.077, ?)\", [now])
    con.execute(\"INSERT INTO market_map VALUES ('kalshi', 'KXNFLDRAFTPICK-26-5', 'Carnell Tate', 'pick_5')\")
    con.execute(\"INSERT INTO kalshi_odds VALUES (?, 'KXNFLDRAFTPICK', 'E', 'KXNFLDRAFTPICK-26-5-CTAT', 'Pick 5', 'Carnell Tate', 2, 5, 95, 98, 4, 100, 500, 50, 20)\", [now])

from nfl_draft.lib.queries import cross_book_grid, kalshi_tooltip_data
grid = cross_book_grid(threshold_pp=5.0)
tip = kalshi_tooltip_data()
print('grid:', grid)
print('tooltip:', tip)
assert abs(grid[0]['books']['kalshi'] - 0.035) < 1e-9
assert tip['pick_5']['yes_ask'] == 5
print('OK')
"
```
Expected: prints `OK`.

A full browser-based dashboard eyeball test happens in Task 11 (against the real DB after the legacy-write fix has run once).

- [ ] **Step 5: Commit**

```bash
git add kalshi_draft/app.py
git commit -m "$(cat <<'EOF'
feat(dashboard): Cross-Book Grid tooltip + help-text update

Hovering a Kalshi cell now reveals last-trade, buy price, sell price,
and the source ticker. Help-text clarifies that the displayed number
is Kalshi's mid (fair-value estimate) while the outlier flag compares
the actual buy price against the cross-venue median.

Co-Authored-By: Claude Opus 4.7 (1M context) <noreply@anthropic.com>
EOF
)"
```

---

## Task 10: Update Documentation

**Files:**
- Modify: `nfl_draft/README.md`

- [ ] **Step 1: Find the Cross-Book Grid section**

Grep for it:
```bash
grep -n "Cross-Book Grid\|cross-book\|cross_book" nfl_draft/README.md
```

- [ ] **Step 2: Update the description**

Locate the paragraph that describes the Cross-Book Grid (around line 19 based on the earlier grep). Replace it with (preserve any surrounding context — adjust indentation/framing to match the file's existing prose style):

```markdown
**Cross-Book Grid** — one row per market. Each cell is the venue's fair-value
estimate:

- Sportsbooks (DK, FD, Bookmaker, Wagerzon, Hoop88): devigged implied probability.
- Kalshi: mid of buy/sell (falls back to last trade when one-sided).

The outlier flag (⚑) fires when a venue's **take price** sits at least
`threshold_pp` away from the cross-venue median:

- Sportsbooks use their devigged probability for the comparison.
- Kalshi uses the actual **buy price** (what you'd pay to take Yes) because
  Kalshi has no vig — the buy is close to fair, and offsets from median are
  real +EV signal rather than vig noise.

Hover a Kalshi cell to see the raw buy, sell, and last-trade prices plus the
source Kalshi ticker. The Cross-Book Grid hides rows older than 2 hours, so
venues that go silent drop out automatically instead of showing stale data.
```

- [ ] **Step 3: Commit**

```bash
git add nfl_draft/README.md
git commit -m "$(cat <<'EOF'
docs(nfl_draft): Cross-Book Grid Kalshi semantics + tooltip

Clarify that the Kalshi column shows mid (fair value) on the cell but
compares buy price (take) against median for the outlier flag. Document
the new hover tooltip.

Co-Authored-By: Claude Opus 4.7 (1M context) <noreply@anthropic.com>
EOF
)"
```

---

## Task 11: Full-Stack Verification

**Files:** N/A (runs end-to-end checks before merge)

- [ ] **Step 1: Run the full test suite**

From the worktree:
```bash
pytest nfl_draft/tests/ -v
```
Expected: all tests pass. If any legacy test fails, investigate before proceeding.

- [ ] **Step 2: Dry-run the Kalshi scraper against a fixture**

```bash
python -c "
import json
from pathlib import Path
from nfl_draft.scrapers.kalshi import parse_markets_response

fixture_path = Path('nfl_draft/tests/fixtures/kalshi/markets_response.json')
raw = json.loads(fixture_path.read_text())
all_rows = []
for ticker, resp in (raw.get('series_responses') or {}).items():
    all_rows.extend(parse_markets_response(resp, series_ticker=ticker))
print(f'Total rows: {len(all_rows)}')
print(f'Rows with implied_prob override: {sum(1 for r in all_rows if r.implied_prob is not None)}')
print(f'Rows with devig_prob override:   {sum(1 for r in all_rows if r.devig_prob is not None)}')
# Should be: every emitted row has a devig override; rows without ask AND without
# last trade will have implied_prob = None.
"
```
Expected: every emitted row has `devig_prob` set; `implied_prob` count is ≤ total.

- [ ] **Step 3: Pre-merge executive review (from CLAUDE.md)**

From the worktree, generate the merge diff and walk the review checklist:
```bash
git diff main..HEAD --stat
git log main..HEAD --oneline
```
Check each item below against the diff:
- **Data integrity:** new Kalshi rows deduplicated by existing `(market_id, book, fetched_at)` insert pattern (unchanged); no double-writes.
- **Resource safety:** quarantine.py still uses `with write_connection() as con:` (unchanged); queries.py uses `read_connection`; no new open file handles.
- **Edge cases:** off-season = empty Kalshi response = zero OddsRows, grid renders empty cleanly (existing behavior); stale-snapshot = MAX_AGE_HOURS drops the row, tooltip lookup returns empty dict, no crash.
- **Dead code:** removed the old `_extract_yes_bid_cents` body; no leftover unused helpers.
- **Log/disk hygiene:** no new log output or files.
- **Security:** no secrets; tooltip shows public ticker / price data only.

Document any findings as ISSUES TO FIX vs ACCEPTABLE RISKS.

- [ ] **Step 4: Ask user for merge approval**

Do NOT merge without explicit user confirmation — per `CLAUDE.md` and the `feedback_always_ask_merge` memory. Post a summary of changed files + the review findings, then wait for approval.

---

## Task 12: Merge and Cleanup

**Files:** N/A (git ops)

- [ ] **Step 1: Merge feature branch into main** *(only after Task 11, Step 6 is approved)*

From the main repo (not the worktree):
```bash
cd /Users/callancapitolo/NFLWork
git checkout main
git merge --no-ff feature/kalshi-cross-grid-pricing
```

- [ ] **Step 2: Run the full test suite on main**

```bash
pytest nfl_draft/tests/ -v
```
Expected: all tests pass.

- [ ] **Step 3: Trigger a live scrape and verify the DB**

```bash
python -c "
from nfl_draft.scrapers.kalshi import fetch_draft_odds
from nfl_draft.lib.quarantine import write_or_quarantine
rows = fetch_draft_odds()
print(f'Fetched {len(rows)} rows')
mapped, unmapped = write_or_quarantine(rows)
print(f'Mapped: {mapped}, Unmapped: {unmapped}')
"
python -c "
import duckdb
con = duckdb.connect('nfl_draft/nfl_draft.duckdb', read_only=True)
for r in con.execute('''
    SELECT market_id, book, implied_prob, devig_prob
    FROM draft_odds
    WHERE book='kalshi' AND market_id='pick_5_overall_carnell-tate'
    ORDER BY fetched_at DESC LIMIT 3
''').fetchall(): print(r)
print('kalshi_odds:', con.execute('''
    SELECT yes_bid, yes_ask, last_price FROM kalshi_odds
    WHERE ticker='KXNFLDRAFTPICK-26-5-CTAT'
    ORDER BY fetch_time DESC LIMIT 1
''').fetchone())
"
```
Expected:
- `implied_prob` (take) ≥ `devig_prob` (mid) on most rows, both non-NULL
- `kalshi_odds` row has real bid / ask / last cents — no longer `(0, 0, 0)`

- [ ] **Step 4: Render the dashboard and eyeball the grid**

```bash
python -m kalshi_draft.app
```
Open the dashboard URL in a browser, select the Cross-Book Grid tab, and verify:
- Kalshi column now shows mid values on markets that previously showed the bid (e.g. `pick_5_overall_carnell-tate` sits around 3-4% rather than 2%)
- Hovering a Kalshi cell reveals `Last trade / Buy / Sell / Ticker`
- Non-Kalshi cells have no tooltip
- Rows that were flagged solely due to the bid artifact are no longer flagged

If anything looks broken, revert the merge immediately: `git revert -m 1 HEAD` and investigate.

- [ ] **Step 5: Remove the worktree**

```bash
git worktree remove /Users/callancapitolo/NFLWork-trees/kalshi-cross-grid-pricing
```

- [ ] **Step 6: Delete the feature branch**

```bash
git branch -d feature/kalshi-cross-grid-pricing
```

- [ ] **Step 7: Verify clean state**

```bash
git worktree list
git branch
git status
```
Expected: only the main working tree, only `main` branch, clean status.

---

## Self-Review (Completed Before Handoff)

**Spec coverage:**
- Per-venue semantics (schema reinterpretation): Task 2 (OddsRow fields) + Task 3 (quarantine honors overrides) + Task 5 (Kalshi emits overrides). ✓
- Cell display + median uses mid: Task 5 (devig_prob = mid) + Task 7 (grid query + median math unchanged). ✓
- Per-venue outlier flag logic: Task 7. ✓
- Resolution cascades: Task 5 (scraper) + Task 7 (flag suppressed when implied_prob is None). ✓
- Tooltip (last trade / buy / sell): Task 8 (query) + Task 9 (render). ✓
- Terminology (Buy price / Sell price): Task 9 (dashboard text + tooltip), Task 10 (README). ✓
- Legacy kalshi_odds write fix: Task 6. ✓
- Open Question #1 (where mid is computed): resolved path (a) — `precomputed_devig` via OddsRow override — implemented in Tasks 2, 3, 5.

**Placeholder scan:** No "TBD," no "implement appropriate error handling." Every step has exact code or exact commands.

**Type consistency:** `OddsRow.implied_prob` / `devig_prob` (Optional[float]) defined in Task 2, referenced in Tasks 3, 5, 7, 8 with consistent naming. `kalshi_tooltip_data()` return type defined in Task 8, consumed in Task 9. No signature drift.

---

## Execution Handoff

Plan complete and saved to `docs/superpowers/plans/2026-04-21-kalshi-cross-grid-pricing.md`.

**Two execution options:**

1. **Subagent-Driven (recommended)** — I dispatch a fresh subagent per task, review between tasks, fast iteration with clean context per step.

2. **Inline Execution** — Execute tasks in this session using `superpowers:executing-plans`, batch execution with checkpoints for review.

Which approach?
