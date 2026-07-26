# Kalshi Single-Leg Anchor (issue #23) Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Use Kalshi's own real-time single-leg markets as an independent anchor — reject quotes whose combo fair is inconsistent with the live marginals (`corr_sanity`), and cancel resting quotes when a constituent single jumps after we quoted (`constituent_jump`).

**Architecture:** One new pure-math module `kalshi_mlb_mm/singles.py` owns everything about constituent singles: devigging a 2-way Kalshi book from top-of-book, per-leg marginals, the Fréchet/premium sanity test, and jump detection. `main.py` calls it from two places — the quote path (reusing the leg snapshot #17 already fetches, so **zero extra API calls per quote**) and the 10s risk sweep (one GET per *distinct* constituent ticker across all open quotes). The placement baseline for the jump check is `live_quotes.leg_prices_json`, already written by #17 — **no schema change**.

**Tech Stack:** Python 3, DuckDB, pytest, `kalshi_common.fair_value.devig_two_way` (probit), `kalshi_common.legset`, `kalshi_common.auth_client`.

## Global Constraints

- Branch: `worktree-mm-issue-23-singles-anchor`, worktree `.claude/worktrees/mm-issue-23-singles-anchor`, based on local `main` @ `060a9ce` (NOT `origin/main` — local is 3 merges ahead).
- Use `python3`, never `python` (no `python` on PATH).
- Full suite must stay green: `python3 -m pytest kalshi_mlb_mm/tests/ kalshi_common/tests/` — **1235 passing on main today**.
- Never `print()` in the bot — use the module `log` (Python `logging`).
- `research.emit(...)` can never raise into the trading loop; all new firehose calls go through it.
- Nothing in the risk sweep or quote path may raise: every new helper returns `None`/empty on bad input rather than throwing.
- Decision reasons must be **distinct strings** so the #14 report (which reads reason vocabularies from data, not a hardcoded list — verified: `report.py` has no skip-reason constants) picks them up automatically.
- Do NOT start/stop/restart the live maker. Do NOT merge without explicit user approval.

---

## Design decisions locked before coding

These were settled with the user across the scoping conversation. Do not re-litigate them mid-implementation.

1. **No ring buffer.** The breaker compares the placement baseline against the current read — two numbers. Last-N history buys only flicker-smoothing, and the asymmetry argues against it: a spurious cancel costs one re-postable quote, a missed cancel costs a pickoff. Cancel eagerly. Add a 2-sample confirmation later only if tests show real flicker.
2. **Placement baseline = `live_quotes.leg_prices_json`** (raw bid/ask per leg, written at quote time by #17). Devigged mid is derived from it. No new column, no migration.
3. **Item 4 is a shared-primitive refactor, not a second veto.** `_market_bid_ask` / `_leg_market_prices` / `_singles_moved` move into `singles.py`; `main.py` imports them under their existing private names so #17's tests and monkeypatches keep working unchanged.
4. **The confirm path keeps its live fetch.** The sweep's constituent prices are up to 10s old; a last look reading 10s-old prices is not a last look. Shared *code*, never shared *values*.
5. **Fréchet gates by default; the premium band ships log-only.** Fréchet is parameter-free and mathematically unarguable. A guessed `[0.5, 2.0]` band could kill legitimately correlated same-game combos, so `CORR_SANITY_PREMIUM_ENABLED` defaults `False` — it logs from day one and gets switched on once the firehose shows the real premium distribution.
6. **Degenerate price shapes skip the gate, they don't reject the quote.** A `yes_ask=100` / empty 0–100 book is *uninformative*, not alarming. `marginals_for_legs` returns `None`, the firehose records the miss, and the quote proceeds on book consensus alone (which is exactly the status quo before this ticket).
7. **An unreadable constituent during a sweep is NOT a cancel.** Cancelling on an API blip would flush the whole book on a transient 500. The fill-blocking backstop already exists: #17's confirm veto fails closed, so a pickoff on an unread constituent still cannot fill.
8. **`constituent_jump` is evaluated before `book_drift`** — it is the more specific and faster signal, so it wins the logged reason when both would fire.
9. **A `corr_sanity` rejection blocks the new quote only**; it does not cancel an existing resting quote on that combo. Resting quotes are the breaker's job.

---

## File structure

| File | Responsibility |
|---|---|
| `kalshi_mlb_mm/singles.py` (**create**) | Everything about Kalshi constituent singles: 2-way devig from top-of-book, per-leg marginals, Fréchet/premium sanity, jump detection, and the relocated fetch/compare helpers. Pure except `fetch_market_prices`, which is the only function that touches the API. |
| `kalshi_mlb_mm/main.py` (**modify**) | Import the relocated helpers; wire the `corr_sanity` gate into the quote path; wire the constituent-jump breaker into `_risk_sweep_tick`. |
| `kalshi_mlb_mm/config.py` (**modify**) | 5 new knobs with rationale comments. |
| `kalshi_mlb_mm/tests/test_singles.py` (**create**) | Unit tests for the pure math, incl. every degenerate price shape. |
| `kalshi_mlb_mm/tests/test_corr_sanity_gate.py` (**create**) | Quote-path integration: Fréchet reject, premium reject/log-only, degenerate pass-through, firehose event. |
| `kalshi_mlb_mm/tests/test_constituent_jump.py` (**create**) | Sweep integration: jump+quiet→cancel, jump+moving books→no cancel (pre-pivot), unconditional mode, quiet constituents→untouched, unreadable→untouched, API call bound. |
| `kalshi_mlb_mm/README.md` (**modify**) | Defense hierarchy items + config table rows + documented per-sweep API bound. |
| `CLAUDE.md` (**modify**) | One-line maker bullet update. |

---

### Task 1: `singles.py` — the constituent-singles module

**Files:**
- Create: `kalshi_mlb_mm/singles.py`
- Create: `kalshi_mlb_mm/tests/test_singles.py`
- Modify: `kalshi_mlb_mm/main.py:473-530` (delete the three relocated helpers, add the import)

**Interfaces:**
- Consumes: `kalshi_common.fair_value.devig_two_way`, `kalshi_common.auth_client.api`
- Produces:
  - `market_bid_ask(mkt: dict) -> tuple[float, float] | None`
  - `leg_market_prices(legs: list[dict]) -> dict | None`
  - `singles_moved(snapshot: dict, fresh: dict) -> bool`
  - `fetch_market_prices(tickers) -> dict[str, dict]`
  - `devigged_yes(yes_bid, yes_ask) -> float | None`
  - `marginals_for_legs(snapshot: dict, legs: list[dict]) -> list[float] | None`
  - `CorrSanity` dataclass with fields `ok, reason, baseline_independent, premium, frechet_lo, frechet_hi`
  - `corr_sanity(combo_fair, marginals, premium_min, premium_max) -> CorrSanity | None`
  - `jumped_tickers(baseline: dict, current: dict, threshold: float) -> dict[str, float]`

- [ ] **Step 1: Write the failing tests for the pure math**

Create `kalshi_mlb_mm/tests/test_singles.py`:

```python
"""Issue #23: Kalshi constituent-singles anchor — pure math.

Every function here must fail SOFT (return None / skip) on a degenerate
Kalshi book rather than raise: `yes_ask=100` (the divide-by-zero gotcha),
an empty 0-100 book, a crossed book, or a missing ticker.
"""
import math

from kalshi_mlb_mm import singles

SPREAD = "KXMLBSPREAD-25JUN271905TEXLAA-LAA2"
TOTAL = "KXMLBTOTAL-25JUN271905TEXLAA-9"
LEGS = [{"market_ticker": SPREAD, "event_ticker": "E", "side": "yes"},
        {"market_ticker": TOTAL, "event_ticker": "E", "side": "yes"}]
SNAP = {SPREAD: {"yes_bid": 0.45, "yes_ask": 0.47},
        TOTAL: {"yes_bid": 0.52, "yes_ask": 0.54}}


def test_devigged_yes_removes_the_spread_overround():
    """Buying YES costs the ask (0.47); buying NO costs 1 - bid (0.55).
    Those imply 1.02 of probability — the 2c spread. Devigging returns the
    fair YES inside the quoted band."""
    p = singles.devigged_yes(0.45, 0.47)
    assert p is not None
    assert 0.45 < p < 0.47


def test_devigged_yes_is_symmetric_about_a_balanced_book():
    assert abs(singles.devigged_yes(0.49, 0.51) - 0.50) < 1e-9


def test_devigged_yes_returns_none_on_degenerate_shapes():
    assert singles.devigged_yes(0.45, 1.00) is None    # yes_ask=100 gotcha
    assert singles.devigged_yes(0.0, 1.00) is None     # empty book
    assert singles.devigged_yes(0.0, 0.03) is None     # no_ask = 1.0
    assert singles.devigged_yes(0.60, 0.40) is None    # crossed
    assert singles.devigged_yes(0.50, 0.50) is None    # locked, no spread
    assert singles.devigged_yes(None, 0.47) is None
    assert singles.devigged_yes("x", 0.47) is None
    assert singles.devigged_yes(float("nan"), 0.47) is None


def test_marginals_follow_the_side_each_leg_takes():
    yes_side = singles.marginals_for_legs(SNAP, LEGS)
    no_legs = [dict(LEGS[0], side="no"), LEGS[1]]
    no_side = singles.marginals_for_legs(SNAP, no_legs)
    assert abs(yes_side[0] + no_side[0] - 1.0) < 1e-9
    assert abs(yes_side[1] - no_side[1]) < 1e-12


def test_marginals_fail_closed_on_missing_or_degenerate_leg():
    assert singles.marginals_for_legs({SPREAD: SNAP[SPREAD]}, LEGS) is None
    bad = {SPREAD: SNAP[SPREAD], TOTAL: {"yes_bid": 0.0, "yes_ask": 1.0}}
    assert singles.marginals_for_legs(bad, LEGS) is None


def test_corr_sanity_accepts_a_plausible_positive_correlation():
    """p1=0.50, p2=0.50 -> independent baseline 0.25, Frechet [0.0, 0.50].
    A same-game fair of 0.30 is a 1.2x premium: legitimate, must pass."""
    s = singles.corr_sanity(0.30, [0.50, 0.50], 0.5, 2.0)
    assert s.ok is True and s.reason is None
    assert abs(s.baseline_independent - 0.25) < 1e-12
    assert abs(s.premium - 1.2) < 1e-12
    assert (s.frechet_lo, s.frechet_hi) == (0.0, 0.50)


def test_corr_sanity_rejects_above_the_frechet_upper_bound():
    """A joint can never exceed its smallest marginal. 0.55 > min(0.50,0.60)."""
    s = singles.corr_sanity(0.55, [0.50, 0.60], 0.5, 2.0)
    assert s.ok is False and s.reason == "frechet"


def test_corr_sanity_rejects_below_the_frechet_lower_bound():
    """p1+p2-1 = 0.50 is the floor when both legs are likely."""
    s = singles.corr_sanity(0.40, [0.80, 0.70], 0.5, 2.0)
    assert s.ok is False and s.reason == "frechet"


def test_corr_sanity_flags_premium_band_only_inside_frechet():
    """0.13 sits inside Frechet [0, 0.5] but is a 0.52x premium on a 0.25
    baseline -> premium violation, not a Frechet one."""
    s = singles.corr_sanity(0.13, [0.50, 0.50], 0.6, 2.0)
    assert s.ok is False and s.reason == "premium"


def test_corr_sanity_returns_none_without_a_usable_anchor():
    assert singles.corr_sanity(0.3, [], 0.5, 2.0) is None
    assert singles.corr_sanity(None, [0.5, 0.5], 0.5, 2.0) is None
    assert singles.corr_sanity(float("nan"), [0.5, 0.5], 0.5, 2.0) is None


def test_jumped_tickers_reports_only_moves_past_the_threshold():
    cur = {SPREAD: {"yes_bid": 0.55, "yes_ask": 0.57},   # ~+0.10
           TOTAL: {"yes_bid": 0.525, "yes_ask": 0.545}}  # ~+0.005
    moved = singles.jumped_tickers(SNAP, cur, 0.03)
    assert set(moved) == {SPREAD}
    assert moved[SPREAD] > 0.03


def test_jumped_tickers_skips_unreadable_or_degenerate_sides():
    """No signal is NOT a jump — an API blip must never flush the book."""
    assert singles.jumped_tickers(SNAP, {}, 0.03) == {}
    degenerate = {SPREAD: {"yes_bid": 0.0, "yes_ask": 1.0},
                  TOTAL: SNAP[TOTAL]}
    assert singles.jumped_tickers(SNAP, degenerate, 0.03) == {}
```

- [ ] **Step 2: Run the tests to verify they fail**

Run: `python3 -m pytest kalshi_mlb_mm/tests/test_singles.py -v`
Expected: FAIL — `ModuleNotFoundError: No module named 'kalshi_mlb_mm.singles'`

- [ ] **Step 3: Create `kalshi_mlb_mm/singles.py`**

```python
"""Kalshi constituent single-leg markets as an independent anchor (issue #23).

Every leg of an MLB combo IS its own 2-way Kalshi market (KXMLBSPREAD-*,
KXMLBTOTAL-*, KXMLBGAME-*), and those books move in REAL TIME while our
sportsbook SGP scrapes refresh every ~150-165s. That makes them the only
pricing input we have that is both independent of the books and fast.

Two consumers, both in main.py:
  * quote path  -> corr_sanity(): is the book-consensus combo fair even
    consistent with the live marginals? (Frechet bounds + correlation premium)
  * risk sweep  -> jumped_tickers(): did a constituent move after we quoted?

Side effects: `leg_market_prices` and `fetch_market_prices` issue GET
/markets/{ticker}. Everything else is pure. No DB access.
"""
import math

from kalshi_common import auth_client, fair_value
from dataclasses import dataclass
import logging

log = logging.getLogger(__name__)


def market_bid_ask(mkt: dict) -> tuple[float, float] | None:
    """(yes_bid, yes_ask) in DOLLARS from either Kalshi market response shape
    — int cents (`yes_bid`) or string-dollar (`yes_bid_dollars`); see the
    kalshi_price_gotchas note on *_dollars vs int shapes."""
    bid, ask = mkt.get("yes_bid"), mkt.get("yes_ask")
    if bid is not None and ask is not None:
        return float(bid) / 100.0, float(ask) / 100.0
    bid, ask = mkt.get("yes_bid_dollars"), mkt.get("yes_ask_dollars")
    if bid is not None and ask is not None:
        return float(bid), float(ask)
    return None


def leg_market_prices(legs: list[dict]) -> dict | None:
    """Raw Kalshi odds for every leg of a combo: {leg market_ticker:
    {"yes_bid": $, "yes_ask": $}}. Each leg IS its own Kalshi singles market,
    so this is one fast GET per leg (~<1s for a 2-3 leg combo). Returns None
    if ANY leg can't be read (fail-safe: no partial baselines)."""
    try:
        out = {}
        for leg in legs:
            mt = str(leg.get("market_ticker") or "")
            if not mt:
                return None
            if mt in out:
                continue
            _status, body, _hdrs = auth_client.api("GET", f"/markets/{mt}")
            mkt = body.get("market") if isinstance(body, dict) else None
            prices = market_bid_ask(mkt) if isinstance(mkt, dict) else None
            if prices is None:
                return None
            out[mt] = {"yes_bid": prices[0], "yes_ask": prices[1]}
        return out or None
    except Exception as e:
        log.warning("[leg_snapshot] fetch failed: %s", e)
        return None


def fetch_market_prices(tickers) -> dict[str, dict]:
    """Best-effort {ticker: {"yes_bid": $, "yes_ask": $}} for a SET of tickers.

    Unlike `leg_market_prices` this does NOT fail closed: unreadable tickers
    are simply absent from the result. That is deliberate and is the risk
    sweep's contract — a transient 500 must not look like a price move and
    flush every resting quote. The fill-blocking backstop is #17's confirm
    veto, which does fail closed.

    Cost: exactly one GET per DISTINCT ticker. Callers dedup before calling.
    """
    out = {}
    for mt in tickers:
        if not mt:
            continue
        try:
            _status, body, _hdrs = auth_client.api("GET", f"/markets/{mt}")
            mkt = body.get("market") if isinstance(body, dict) else None
            prices = market_bid_ask(mkt) if isinstance(mkt, dict) else None
            if prices is not None:
                out[mt] = {"yes_bid": prices[0], "yes_ask": prices[1]}
        except Exception as e:
            log.warning("[constituent] fetch failed for %s: %s", mt, e)
    return out


def singles_moved(snapshot: dict, fresh: dict) -> bool:
    """True if ANY leg's yes_bid or yes_ask differs between the quote-time
    snapshot and the fresh read. Zero tolerance by design: Kalshi prices move
    in 1c ticks, so 'moved at all' == 'moved >= one tick'. A leg-set mismatch
    counts as moved (fail-safe)."""
    if set(snapshot) != set(fresh):
        return True
    for mt, snap in snapshot.items():
        cur = fresh[mt]
        if (abs(float(snap["yes_bid"]) - float(cur["yes_bid"])) > 1e-9
                or abs(float(snap["yes_ask"]) - float(cur["yes_ask"])) > 1e-9):
            return True
    return False


def devigged_yes(yes_bid, yes_ask) -> float | None:
    """Probit-devigged fair P(YES) for one 2-way Kalshi market, from its
    top-of-book. Buying YES costs `yes_ask`; buying NO costs `1 - yes_bid`
    (the NO ask). Those two implied probabilities sum to more than 1 by
    exactly the quoted spread — the 2-way overround `devig_two_way` removes.

    Returns None, never raises, on every degenerate shape: `yes_ask = 1.00`
    (the known divide-by-zero gotcha), `yes_bid = 0`, an empty 0-100 book, or
    a crossed/locked book with no spread to devig. Callers treat None as
    "this leg has no usable anchor".
    """
    try:
        bid, ask = float(yes_bid), float(yes_ask)
    except (TypeError, ValueError):
        return None
    if not (math.isfinite(bid) and math.isfinite(ask)):
        return None
    no_ask = 1.0 - bid
    if not (0.0 < ask < 1.0) or not (0.0 < no_ask < 1.0):
        return None
    if bid >= ask:
        return None
    devigged = fair_value.devig_two_way(1.0 / ask, 1.0 / no_ask)
    if devigged is None:
        return None
    fair_yes = float(devigged[0])
    return fair_yes if 0.0 < fair_yes < 1.0 else None


def marginals_for_legs(snapshot: dict, legs: list[dict]) -> list[float] | None:
    """Devigged probability of the side EACH LEG ACTUALLY TAKES, in leg order.

    A leg resting on the NO side contributes 1 - P(YES). Returns None if any
    leg is missing from the snapshot or its book is degenerate: the premium
    gate needs every marginal or it has no baseline at all.
    """
    if not legs or not snapshot:
        return None
    out = []
    for leg in legs:
        prices = snapshot.get(str(leg.get("market_ticker") or ""))
        if not prices:
            return None
        p_yes = devigged_yes(prices.get("yes_bid"), prices.get("yes_ask"))
        if p_yes is None:
            return None
        p = p_yes if str(leg.get("side") or "yes").lower() == "yes" else 1.0 - p_yes
        if not (0.0 < p < 1.0):
            return None
        out.append(p)
    return out


@dataclass(frozen=True)
class CorrSanity:
    """Verdict of the correlation-premium / Frechet check on one combo fair."""
    ok: bool
    reason: str | None            # None | "frechet" | "premium"
    baseline_independent: float
    premium: float
    frechet_lo: float
    frechet_hi: float


def corr_sanity(combo_fair, marginals: list[float],
                premium_min: float, premium_max: float) -> CorrSanity | None:
    """Is `combo_fair` consistent with the live single-leg marginals?

    Two independent tests, most-certain first:

    1. FRECHET (parameter-free, always correct). For ANY correlation
       structure whatsoever, a joint probability obeys
       `max(0, sum(p) - (n-1)) <= p_joint <= min(p)`. A violation is not a
       judgement call — it is arithmetic saying our fair cannot be a
       probability of these legs jointly.
    2. PREMIUM BAND (heuristic, tunable). `premium = combo_fair / prod(p)`
       is the correlation multiplier the books are implying. Same-game legs
       legitimately run well above 1 (run line + moneyline ~2x), so the band
       is deliberately wide and ships log-only.

    Returns None when there is nothing to check (no marginals, unusable
    fair) so the caller can log a miss and carry on.
    """
    if not marginals:
        return None
    try:
        fair = float(combo_fair)
    except (TypeError, ValueError):
        return None
    if not math.isfinite(fair) or not (0.0 < fair < 1.0):
        return None
    baseline = 1.0
    for p in marginals:
        baseline *= p
    if baseline <= 0.0:
        return None
    lo = max(0.0, sum(marginals) - (len(marginals) - 1))
    hi = min(marginals)
    premium = fair / baseline
    reason = None
    if not (lo <= fair <= hi):
        reason = "frechet"
    elif not (premium_min <= premium <= premium_max):
        reason = "premium"
    return CorrSanity(ok=reason is None, reason=reason,
                      baseline_independent=baseline, premium=premium,
                      frechet_lo=lo, frechet_hi=hi)


def jumped_tickers(baseline: dict, current: dict,
                   threshold: float) -> dict[str, float]:
    """{ticker: |delta devigged P(YES)|} for tickers that moved more than
    `threshold` between the two raw bid/ask snapshots.

    Tickers absent from either snapshot, or degenerate on either side, are
    SKIPPED — absence of signal is not a jump (see fetch_market_prices).
    Side is irrelevant: |delta P(YES)| == |delta P(NO)|.
    """
    out = {}
    if not baseline or not current:
        return out
    for mt, base in baseline.items():
        cur = current.get(mt)
        if not cur or not isinstance(base, dict):
            continue
        p_base = devigged_yes(base.get("yes_bid"), base.get("yes_ask"))
        p_cur = devigged_yes(cur.get("yes_bid"), cur.get("yes_ask"))
        if p_base is None or p_cur is None:
            continue
        delta = abs(p_cur - p_base)
        if delta > threshold:
            out[mt] = delta
    return out
```

- [ ] **Step 4: Run the tests to verify they pass**

Run: `python3 -m pytest kalshi_mlb_mm/tests/test_singles.py -v`
Expected: PASS (13 tests)

- [ ] **Step 5: Relocate the three helpers out of `main.py`**

Delete `_market_bid_ask`, `_leg_market_prices`, and `_singles_moved` from `main.py` (currently lines 473–530) and replace them with an import alias block placed with the other module imports at the top of `main.py`:

```python
# #23: the constituent-singles helpers moved to kalshi_mlb_mm/singles.py so the
# quote-time gate, the risk sweep, and #17's confirm veto all share ONE
# implementation of "read a Kalshi single and compare it". Imported under the
# original private names so existing call sites and test monkeypatches
# (main._leg_market_prices, main._singles_moved) keep working unchanged.
from kalshi_mlb_mm import singles
from kalshi_mlb_mm.singles import market_bid_ask as _market_bid_ask
from kalshi_mlb_mm.singles import leg_market_prices as _leg_market_prices
from kalshi_mlb_mm.singles import singles_moved as _singles_moved
```

- [ ] **Step 6: Verify the relocation broke nothing**

Run: `python3 -m pytest kalshi_mlb_mm/tests/test_confirm_singles_veto.py kalshi_mlb_mm/tests/test_singles.py -v`
Expected: PASS — all of #17's veto tests (including the three helper-level unit tests that reference `main._market_bid_ask`, `main._singles_moved`, `main._leg_market_prices`) plus the new module tests.

- [ ] **Step 7: Commit**

```bash
git add kalshi_mlb_mm/singles.py kalshi_mlb_mm/tests/test_singles.py kalshi_mlb_mm/main.py
git commit -m "feat(mm): singles.py — Kalshi constituent anchor math (#23)

Devig a 2-way Kalshi book from top-of-book, per-leg marginals honouring
each leg's side, Frechet/premium sanity, and jump detection. Relocates
#17's fetch/compare helpers here so the confirm veto, the quote gate and
the risk sweep share one implementation. All degenerate price shapes
(yes_ask=100, empty, crossed, locked) return None instead of raising."
```

---

### Task 2: `corr_sanity` quote-time gate + firehose logging (spec items 3 & 5)

**Files:**
- Modify: `kalshi_mlb_mm/config.py` (append after the `MIN_AGREEING_BOOKS` block, ~line 109)
- Modify: `kalshi_mlb_mm/main.py` (quote path, immediately after the `leg_snapshot` fetch at ~line 1115)
- Create: `kalshi_mlb_mm/tests/test_corr_sanity_gate.py`

**Interfaces:**
- Consumes: `singles.marginals_for_legs`, `singles.corr_sanity`, `CorrSanity` (Task 1)
- Produces: decision reasons `corr_sanity_frechet` / `corr_sanity_premium`; research event `corr_sanity_check`

- [ ] **Step 1: Write the failing integration tests**

Create `kalshi_mlb_mm/tests/test_corr_sanity_gate.py`:

```python
"""Issue #23 item 3: correlation-premium / Frechet gate at quote time.

The #20 dispersion gate asks whether the BOOKS agree with each other. This
asks whether their consensus is consistent with Kalshi's own live single-leg
prices — books can agree tightly and still be jointly wrong. The two gates
are independent and log distinct reasons.

Scaffolding mirrors test_confirm_singles_veto.py.
"""
import importlib
import json

import pandas as pd

from kalshi_common.leg_types import SPREAD_TOTAL_FAMILY
from mlb_sgp._shared import GameRef

EVT = "KXMLBGAME-25JUN271905TEXLAA"
SPREAD_TICKER = "KXMLBSPREAD-25JUN271905TEXLAA-LAA2"
TOTAL_TICKER = "KXMLBTOTAL-25JUN271905TEXLAA-9"
LEGS = [{"market_ticker": SPREAD_TICKER, "event_ticker": EVT, "side": "yes"},
        {"market_ticker": TOTAL_TICKER, "event_ticker": EVT, "side": "yes"}]
GREF = GameRef(game_id="game1", home_team="Los Angeles Angels",
               away_team="Texas Rangers", commence_time=None)

# 4 equal decimals -> probit devig gives exactly 0.25 per cell, per book.
BOOK_FAIR = 0.25


def _grid_df():
    rows = []
    for book in ("draftkings", "fanduel"):
        for cell in SPREAD_TOTAL_FAMILY:
            rows.append(dict(game_id="game1", combo=cell, period="FG",
                             bookmaker=book, sgp_decimal=3.8, fetch_time=None,
                             spread_line=-1.5, total_line=8.5))
    return pd.DataFrame(rows)


def _setup(monkeypatch, tmp_path, db_name, *, snapshot):
    """Quote path primed to reach the corr_sanity gate; `snapshot` is what
    the leg-price fetch returns (the live marginal anchor)."""
    import kalshi_mlb_mm.config as cfg
    import kalshi_mlb_mm.db as db
    from kalshi_mlb_mm import main
    monkeypatch.setattr(cfg, "DB_PATH", tmp_path / db_name)
    monkeypatch.setattr(cfg, "KILL_FILE", tmp_path / ".kill")
    importlib.reload(db)
    db.init_database()
    monkeypatch.setattr(main, "_SGP_ODDS", _grid_df())
    monkeypatch.setattr(main, "_resolve_game_for_legs", lambda gl: "game1")
    monkeypatch.setattr(main, "_game_ref", lambda gid: GREF)
    monkeypatch.setattr(main, "_ENGINE", None)
    monkeypatch.setattr(main, "_leg_market_prices", lambda legs: snapshot)
    return main, db, cfg


def _gate(main, cfg, monkeypatch, *, frechet=True, premium=False,
          pmin=0.5, pmax=2.0):
    monkeypatch.setattr(cfg, "CORR_SANITY_FRECHET_ENABLED", frechet)
    monkeypatch.setattr(cfg, "CORR_SANITY_PREMIUM_ENABLED", premium)
    monkeypatch.setattr(cfg, "CORR_PREMIUM_MIN", pmin)
    monkeypatch.setattr(cfg, "CORR_PREMIUM_MAX", pmax)


# Marginals ~0.46 and ~0.53 -> baseline ~0.244, Frechet hi = 0.46.
# BOOK_FAIR 0.25 sits inside Frechet and at premium ~1.02 -> passes cleanly.
HEALTHY = {SPREAD_TICKER: {"yes_bid": 0.45, "yes_ask": 0.47},
           TOTAL_TICKER: {"yes_bid": 0.52, "yes_ask": 0.54}}

# Both marginals ~0.10 -> Frechet hi ~0.10, so a 0.25 combo fair is
# arithmetically impossible: the joint exceeds its smallest leg.
FRECHET_VIOLATING = {SPREAD_TICKER: {"yes_bid": 0.09, "yes_ask": 0.11},
                     TOTAL_TICKER: {"yes_bid": 0.09, "yes_ask": 0.11}}

# Marginals ~0.90 each -> baseline ~0.81, Frechet [0.80, 0.90]. A 0.25 fair
# breaks Frechet too, so for a PREMIUM-only violation we need a fair inside
# Frechet but off-baseline: marginals ~0.70/0.71 give baseline ~0.497,
# Frechet [0.41, 0.70]; a 0.25 fair -> premium 0.50, inside Frechet? No.
# Use marginals ~0.55/0.52: baseline ~0.286, Frechet [0.07, 0.52];
# fair 0.25 -> premium ~0.874, so widen the band to trip it deliberately.
PREMIUM_TEST = {SPREAD_TICKER: {"yes_bid": 0.54, "yes_ask": 0.56},
                TOTAL_TICKER: {"yes_bid": 0.51, "yes_ask": 0.53}}

# yes_ask = 1.00 on one leg -> no usable anchor -> gate skipped, no crash.
DEGENERATE = {SPREAD_TICKER: {"yes_bid": 0.0, "yes_ask": 1.00},
              TOTAL_TICKER: {"yes_bid": 0.52, "yes_ask": 0.54}}


def _decisions(db):
    with db.connect(read_only=True) as con:
        return con.execute(
            "SELECT decision, reason FROM quote_decisions "
            "ORDER BY observed_at").fetchall()
```

The remaining tests in this file drive `main._discovery_tick` (or the
narrowest existing entry point the other quote-path tests use — mirror
`test_dispersion_gate.py`'s harness exactly) with an in-scope RFQ over
`LEGS`, and assert:

```python
def test_frechet_violation_rejects_the_quote(monkeypatch, tmp_path):
    """A combo fair above min(marginals) cannot be a joint probability."""
    # ... drive the tick with snapshot=FRECHET_VIOLATING, gate frechet=True
    # assert ("skipped", "corr_sanity_frechet") in _decisions(db)
    # assert the gateway received NO submit_quote call


def test_premium_band_logs_but_does_not_gate_by_default(monkeypatch, tmp_path):
    """Ships log-only: a band trip must still quote until we have data."""
    # ... snapshot=PREMIUM_TEST, pmin/pmax narrowed so premium trips,
    #     premium=False -> quote IS submitted, corr_sanity_check event has
    #     reason == "premium"


def test_premium_band_rejects_once_enabled(monkeypatch, tmp_path):
    """Same inputs, premium=True -> reason corr_sanity_premium, no quote."""


def test_healthy_marginals_pass_through(monkeypatch, tmp_path):
    """No violation -> quote submitted, event logged with reason None."""


def test_degenerate_marginals_skip_the_gate_without_crashing(monkeypatch, tmp_path):
    """yes_ask=100 on a leg: uninformative, not alarming. The quote proceeds
    on book consensus alone and the firehose records the miss."""
    # ... snapshot=DEGENERATE -> quote submitted, event marginals is None


def test_corr_sanity_check_event_carries_the_calibration_fields(monkeypatch, tmp_path):
    """Item 5: premium + marginals in the firehose per quote — the dataset a
    Phase-3 correlation model is calibrated on."""
    # payload has marginals, baseline_independent, premium, frechet_lo,
    # frechet_hi, combo_fair, n_legs
```

- [ ] **Step 2: Run to verify they fail**

Run: `python3 -m pytest kalshi_mlb_mm/tests/test_corr_sanity_gate.py -v`
Expected: FAIL — `AttributeError: module 'kalshi_mlb_mm.config' has no attribute 'CORR_SANITY_FRECHET_ENABLED'`

- [ ] **Step 3: Add the config knobs**

Append to `kalshi_mlb_mm/config.py` after the `MIN_AGREEING_BOOKS` block:

```python
# Correlation sanity vs Kalshi's own single-leg markets (issue #23, spec §13).
# Every leg of a combo is its own 2-way Kalshi market trading in REAL TIME,
# which makes the devigged marginals the one pricing input we have that is
# independent of the books AND fast. Two tests on the combo fair:
#
#   Frechet:  max(0, Sum(p) - (n-1)) <= combo_fair <= min(p)   [always true]
#   premium:  combo_fair / Prod(p)   in [CORR_PREMIUM_MIN, MAX] [heuristic]
#
# Frechet gates by default: it is parameter-free, and a violation means our
# fair is arithmetically not a joint probability of these legs. The premium
# band ships LOG-ONLY (enabled=False) because same-game legs legitimately run
# well above 1x (run line + moneyline ~2x) and a guessed band would decline
# real business; switch it on once the corr_sanity_check firehose shows the
# true premium distribution. This is a LEVEL check against an outside anchor
# and is deliberately independent of SIGMA_Z_MAX, which is a DISPERSION check
# among the books themselves — tightly-agreeing books can still be jointly
# wrong, and only this gate can see that.
CORR_SANITY_FRECHET_ENABLED = _get("CORR_SANITY_FRECHET_ENABLED", "true").lower() == "true"
CORR_SANITY_PREMIUM_ENABLED = _get("CORR_SANITY_PREMIUM_ENABLED", "false").lower() == "true"
CORR_PREMIUM_MIN = float(_get("CORR_PREMIUM_MIN", "0.5"))
CORR_PREMIUM_MAX = float(_get("CORR_PREMIUM_MAX", "2.0"))
```

Verify the boolean-parsing idiom matches the file's existing convention before committing; if `config.py` already has a `_get_bool` helper, use that instead.

- [ ] **Step 4: Wire the gate into the quote path**

In `main.py`, immediately after the existing `leg_snapshot` fail-closed block (~line 1120) and **before** `qid = gateway.submit_quote(...)`:

```python
        # #23 item 3: correlation sanity against Kalshi's live singles. The
        # leg snapshot we just fetched for #17's veto IS the marginal anchor,
        # so this costs ZERO extra API calls. Independent of the #20 gate:
        # that one asks whether the books agree with each other, this asks
        # whether their consensus is consistent with the real-time single-leg
        # prices. Degenerate books (yes_ask=100, empty) yield no marginals —
        # we log the miss and quote on book consensus alone, exactly as before
        # this ticket, rather than declining on missing information.
        marginals = singles.marginals_for_legs(leg_snapshot, legs)
        sanity = singles.corr_sanity(blended, marginals,
                                     config.CORR_PREMIUM_MIN,
                                     config.CORR_PREMIUM_MAX) if marginals else None
        research.emit("corr_sanity_check", rfq_id=rid, ticker=ticker,
                      payload=dict(game_id=game_id, combo_fair=blended,
                                   n_legs=len(legs), marginals=marginals,
                                   baseline_independent=(
                                       sanity.baseline_independent if sanity else None),
                                   premium=sanity.premium if sanity else None,
                                   frechet_lo=sanity.frechet_lo if sanity else None,
                                   frechet_hi=sanity.frechet_hi if sanity else None,
                                   reason=sanity.reason if sanity else None))
        gated = sanity is not None and (
            (sanity.reason == "frechet" and config.CORR_SANITY_FRECHET_ENABLED)
            or (sanity.reason == "premium" and config.CORR_SANITY_PREMIUM_ENABLED))
        if gated:
            _log_decision("skipped", rfq_id=rid, ticker=ticker, game_id=game_id,
                          reason=f"corr_sanity_{sanity.reason}",
                          book=book_med, blended=blended)
            continue
```

- [ ] **Step 5: Run the gate tests**

Run: `python3 -m pytest kalshi_mlb_mm/tests/test_corr_sanity_gate.py -v`
Expected: PASS

- [ ] **Step 6: Run the neighbouring gate suites for interference**

Run: `python3 -m pytest kalshi_mlb_mm/tests/test_dispersion_gate.py kalshi_mlb_mm/tests/test_uncertainty_margin.py kalshi_mlb_mm/tests/test_confirm_singles_veto.py -v`
Expected: PASS — #19's margin and #20's dispersion gate are untouched; reasons stay distinct.

- [ ] **Step 7: Commit**

```bash
git add kalshi_mlb_mm/config.py kalshi_mlb_mm/main.py kalshi_mlb_mm/tests/test_corr_sanity_gate.py
git commit -m "feat(mm): corr_sanity quote gate vs live Kalshi marginals (#23)

Frechet bounds + correlation premium against the devigged single-leg
markets, reusing #17's quote-time leg snapshot so it costs no extra API
calls. Frechet gates by default; the premium band ships log-only until
the corr_sanity_check firehose shows the real distribution. Distinct
reasons corr_sanity_frechet / corr_sanity_premium keep the #14 report's
funnel readable alongside #20's too_few_books / consensus_dispersion."
```

---

### Task 3: constituent-jump circuit breaker in the risk sweep (spec items 1 & 2)

**Files:**
- Modify: `kalshi_mlb_mm/config.py` (append after the Task 2 block)
- Modify: `kalshi_mlb_mm/main.py::_risk_sweep_tick` (~line 1502) + two new module-level helpers
- Create: `kalshi_mlb_mm/tests/test_constituent_jump.py`

**Interfaces:**
- Consumes: `singles.fetch_market_prices`, `singles.jumped_tickers` (Task 1); `_quote_game_ids`, `_resolve_game_for_legs`, `legset.parse_legs`, `legset.partition_by_game` (existing)
- Produces: sweep cancel reason `constituent_jump`; research event `constituent_jump`; helpers `_open_quote_constituent_tickers`, `_constituent_jumped_games`

- [ ] **Step 1: Write the failing sweep tests**

Create `kalshi_mlb_mm/tests/test_constituent_jump.py` following `test_confirm_singles_veto.py`'s `_setup` shape (seed `live_quotes` rows with `leg_prices_json` as the placement baseline, seed `seen_rfqs.legs_json`, monkeypatch `main._SGP_ODDS`, `main._resolve_game_for_legs`, `main._ENGINE`), then:

```python
def test_jump_with_quiet_books_cancels_quotes_on_that_game(monkeypatch, tmp_path):
    """THE case: a constituent moved 10c since we quoted while our book
    consensus barely moved -> we are the stale side -> pull."""
    # baseline snapshot 0.45/0.47; current 0.55/0.57; book fair unchanged
    # mode "book_quiet"
    # assert live_quotes.status == 'cancelled'
    # assert last quote_decisions row == ("sweep_cancel", "constituent_jump")
    # assert gateway.cancelled == ["q-grid"]


def test_jump_with_moving_books_does_not_fire_in_book_quiet_mode(monkeypatch, tmp_path):
    """Everything moved together -> not a pickoff signal. The existing
    book_drift breaker owns that case and keeps its own reason."""
    # current constituent jumped AND cur_med drifted past the quiet bound
    # assert reason != "constituent_jump"


def test_unconditional_mode_fires_regardless_of_book_movement(monkeypatch, tmp_path):
    """Post-#54 mode: the fair was live at placement, so ANY constituent
    jump during the resting window means the market moved after us."""
    # mode "unconditional", books moving -> still cancelled w/ constituent_jump


def test_quiet_constituents_leave_quotes_untouched(monkeypatch, tmp_path):
    """A sub-threshold wiggle is not a jump."""
    # current 0.455/0.475 (delta ~0.005 < 0.03) -> status stays 'open'


def test_unreadable_constituent_does_not_cancel(monkeypatch, tmp_path):
    """A transient API failure must not flush the book — absence of signal
    is not a jump. #17's confirm veto still fails closed on any fill."""
    # fetch_market_prices returns {} -> status stays 'open'


def test_jump_cancels_sibling_quotes_touching_the_same_game(monkeypatch, tmp_path):
    """Game-level cancel: quote B rests on a different market of the same
    game and never saw the jump itself, but its game is compromised."""
    # two open quotes on game1, only quote A's leg jumped -> BOTH cancelled


def test_missing_placement_baseline_is_skipped_not_cancelled(monkeypatch, tmp_path):
    """Pre-#17 rows have NULL leg_prices_json: no baseline, no jump verdict.
    The tipoff/staleness/book_drift gates still protect them."""


def test_api_calls_are_one_per_distinct_ticker_per_sweep(monkeypatch, tmp_path):
    """Documented rate bound: two quotes sharing a leg cost ONE GET for it."""
    # count auth_client.api calls across a sweep with 2 quotes / 3 distinct
    # tickers -> exactly 3 GETs


def test_sweep_makes_zero_calls_when_flat(monkeypatch, tmp_path):
    """No resting quotes -> no constituent polling at all."""
```

- [ ] **Step 2: Run to verify they fail**

Run: `python3 -m pytest kalshi_mlb_mm/tests/test_constituent_jump.py -v`
Expected: FAIL — `AttributeError: module 'kalshi_mlb_mm.config' has no attribute 'CONSTITUENT_JUMP_THRESHOLD'`

- [ ] **Step 3: Add the config knobs**

```python
# Constituent-jump circuit breaker (issue #23 items 1-2). Our books refresh
# every ~150-165s; the combo's constituent Kalshi singles trade in real time.
# If a constituent's devigged mid moves past this threshold AFTER we quoted,
# the market has moved and our resting quote is the stale side.
CONSTITUENT_JUMP_THRESHOLD = float(_get("CONSTITUENT_JUMP_THRESHOLD", "0.03"))
# Mode flag, flipped by the live-pricing pivot (#54) rather than a rewrite:
#   "book_quiet"    (default, TODAY) — also require our book consensus to have
#                   stayed within CONSTITUENT_BOOK_QUIET_MAX, which separates
#                   "we are stale" from "the whole market moved together".
#                   Quotes are priced off a <=150s cache, so a jump WITH the
#                   books moving is just our own data catching up.
#   "unconditional" (post-#54) — once every quote is priced from a live fetch,
#                   the fair was fresh at placement, so ANY constituent jump
#                   during the resting window means the market moved after us.
CONSTITUENT_JUMP_MODE = _get("CONSTITUENT_JUMP_MODE", "book_quiet")
CONSTITUENT_BOOK_QUIET_MAX = float(_get("CONSTITUENT_BOOK_QUIET_MAX", "0.01"))
```

- [ ] **Step 4: Add the two sweep helpers to `main.py`**

Place immediately above `_risk_sweep_tick`:

```python
def _open_quote_constituent_tickers(live_rows) -> set[str]:
    """Every DISTINCT constituent ticker across all open quotes.

    This set IS the sweep's API cost: one GET per element, per sweep. It is
    bounded by (open quotes x legs) and is usually far smaller, because
    quotes on the same game share legs. Zero open quotes -> zero calls.
    """
    tickers = set()
    for row in live_rows:
        leg_prices_json = row[6]
        if not leg_prices_json:
            continue
        try:
            tickers.update(str(t) for t in json.loads(leg_prices_json))
        except Exception:
            continue
    return tickers


def _constituent_jumped_games(live_rows, current_prices) -> tuple[set, dict]:
    """(game_ids whose constituents jumped since quote time, per-quote detail).

    The placement baseline is live_quotes.leg_prices_json — the raw Kalshi
    odds #17 already snapshots when the quote is submitted, so this needs no
    schema of its own. In "book_quiet" mode a jump only counts when our own
    book consensus stayed put; see CONSTITUENT_JUMP_MODE.

    Fail-safe: any row we cannot evaluate contributes NO jump signal (it is
    still covered by the tipoff / staleness / book-drift gates).
    """
    jumped_games, detail = set(), {}
    if not current_prices:
        return jumped_games, detail
    for qid, _gid, _ticker, book_fair_at_q, _rid, legs_json, leg_prices_json in live_rows:
        if not legs_json or not leg_prices_json:
            continue
        try:
            baseline = json.loads(leg_prices_json)
            raw_legs = json.loads(legs_json)
        except Exception:
            continue
        moves = singles.jumped_tickers(baseline, current_prices,
                                       config.CONSTITUENT_JUMP_THRESHOLD)
        if not moves:
            continue
        if config.CONSTITUENT_JUMP_MODE == "book_quiet":
            cur_med = _current_consensus_fair(legs_json)
            if cur_med is None or book_fair_at_q is None:
                continue                      # cannot establish quiet -> no signal
            if abs(float(cur_med) - float(book_fair_at_q)) >= config.CONSTITUENT_BOOK_QUIET_MAX:
                continue                      # market moved together, not a pickoff
        detail[qid] = moves
        jumped_games.update(_games_for_tickers(raw_legs, set(moves)))
    return jumped_games, detail


def _games_for_tickers(raw_legs: list[dict], tickers: set) -> set:
    """game_ids of the legs whose market_ticker is in `tickers`.

    CanonicalLeg carries the event_ticker but not the market_ticker, so the
    raw leg dicts map ticker -> event_ticker and legset.partition_by_game
    maps event_ticker -> the CanonicalLegs _resolve_game_for_legs consumes.
    """
    out = set()
    try:
        canon = legset.parse_legs(raw_legs)
        if not canon:
            return out
        by_event = legset.partition_by_game(canon)
        events = {str(l.get("event_ticker") or "") for l in raw_legs
                  if str(l.get("market_ticker") or "") in tickers}
        for event_ticker in events:
            game_legs = by_event.get(event_ticker)
            gid = _resolve_game_for_legs(game_legs) if game_legs else None
            if gid:
                out.add(gid)
    except Exception:
        return out
    return out
```

`_current_consensus_fair(legs_json)` is an extraction of the `router.combo_fair(...)` call already inlined in the sweep's book-drift branch — pull it out verbatim into a module-level helper returning `float | None` (swallowing exceptions exactly as the inline version does) so both the drift check and the quiet check use one implementation.

- [ ] **Step 5: Wire the breaker into `_risk_sweep_tick`**

Three edits inside `_risk_sweep_tick`:

(a) add `lq.leg_prices_json` as the 7th selected column:

```python
        live = con.execute(
            "SELECT lq.quote_id, lq.game_id, lq.combo_market_ticker, lq.book_fair, "
            "       lq.rfq_id, sr.legs_json, lq.leg_prices_json "
            "FROM live_quotes lq LEFT JOIN seen_rfqs sr ON lq.rfq_id = sr.rfq_id "
            "WHERE lq.status='open'").fetchall()
```

(b) before the per-quote loop, poll the constituents once:

```python
    # #23: real-time constituent poll. One GET per DISTINCT ticker across all
    # open quotes; skipped entirely when flat or when books are already stale
    # (everything is being cancelled anyway).
    jumped_games, jump_detail = set(), {}
    if live and not books_stale:
        current_prices = singles.fetch_market_prices(
            _open_quote_constituent_tickers(live))
        jumped_games, jump_detail = _constituent_jumped_games(live, current_prices)
        if jumped_games:
            research.emit("constituent_jump",
                          payload=dict(games=sorted(jumped_games),
                                       threshold=config.CONSTITUENT_JUMP_THRESHOLD,
                                       mode=config.CONSTITUENT_JUMP_MODE,
                                       moves={q: m for q, m in jump_detail.items()}))
```

(c) unpack the new column and evaluate the breaker **before** the book-drift
check (it is the more specific, faster signal, so it wins the logged reason
when both would fire). The loop header becomes:

```python
    for qid, game_id, ticker, book_fair_at_q, rid, legs_json, leg_prices_json in live:
```

and inside, immediately after the tipoff/`books_stale` block and before the
existing drift check:

```python
        if not cancel and jumped_games:
            gids = _quote_game_ids(legs_json)
            # gids is None only when the tipoff branch already cancelled this
            # quote, so a None here cannot fail open.
            if gids and (set(gids) & jumped_games):
                cancel, cancel_reason = True, "constituent_jump"
```

- [ ] **Step 6: Run the sweep tests**

Run: `python3 -m pytest kalshi_mlb_mm/tests/test_constituent_jump.py -v`
Expected: PASS

- [ ] **Step 7: Run every suite that touches the sweep**

Run: `python3 -m pytest kalshi_mlb_mm/tests/test_main_smoke.py kalshi_mlb_mm/tests/test_correctness_batch.py kalshi_mlb_mm/tests/test_open_quote_caps.py kalshi_mlb_mm/tests/test_per_game_exposure.py kalshi_mlb_mm/tests/test_n7_n8_n9_n10_n11_n12.py -v`
Expected: PASS — the added SELECT column changes the row tuple width, so any other consumer of that query must be updated in the same edit.

- [ ] **Step 8: Commit**

```bash
git add kalshi_mlb_mm/config.py kalshi_mlb_mm/main.py kalshi_mlb_mm/tests/test_constituent_jump.py
git commit -m "feat(mm): constituent-jump circuit breaker in the risk sweep (#23)

Polls each resting quote's Kalshi constituent singles every sweep (one GET
per DISTINCT ticker) and cancels every quote touching a game whose
constituent moved past CONSTITUENT_JUMP_THRESHOLD since placement. The
placement baseline is #17's leg_prices_json, so no schema change. Two modes
behind CONSTITUENT_JUMP_MODE: book_quiet (default today — a jump only counts
while our book fair stayed put) and unconditional (post-#54 live pricing).
An unreadable constituent is NOT a jump: a transient 500 must not flush the
resting book, and #17's confirm veto still fails closed on any fill."
```

---

### Task 4: Documentation, full suite, pre-merge review

**Files:**
- Modify: `kalshi_mlb_mm/README.md` (defense hierarchy + config table)
- Modify: `CLAUDE.md` (maker bullet)

- [ ] **Step 1: Insert two new defense-hierarchy entries in `README.md`**

Insert after the current item 2 (book-consensus gate) and renumber the rest:

> **3. Correlation sanity vs Kalshi's own singles (issue #23).** Every leg of a combo is its own real-time 2-way Kalshi market. Before submitting, the devigged marginals of all legs are read (the same snapshot #17 takes for its confirm veto — no extra API calls) and the combo fair is checked against them: it must respect the **Fréchet bounds** `max(0, Σp − (n−1)) ≤ fair ≤ min(p)` and its implied **correlation premium** `fair / Πp` must sit inside `[CORR_PREMIUM_MIN, CORR_PREMIUM_MAX]`. This is the only defense that can catch the books being *jointly* wrong — the #20 dispersion gate only sees whether they agree with **each other**, and tightly-agreeing books still get a thin margin under #19. Fréchet gates by default (`corr_sanity_frechet`); the premium band ships log-only (`CORR_SANITY_PREMIUM_ENABLED=false`) because same-game stacks legitimately price well above 1× — every quote emits a `corr_sanity_check` research event carrying marginals, baseline, premium and both bounds, which is both the tuning dataset for the band and the calibration dataset for a Phase-3 correlation model. Degenerate books (`yes_ask=100`, empty, crossed) yield no marginals: the miss is logged and the quote proceeds on book consensus alone rather than declining on missing information.

> **4. Constituent-jump circuit breaker (issue #23).** Our books refresh every ~150–165s; the constituent singles trade in real time. Every risk sweep (`RISK_SWEEP_SEC`, 10s) the bot polls the **distinct** constituent tickers of all resting quotes — **one GET per distinct ticker per sweep, bounded by (open quotes × legs) and zero when flat** — and compares each against the placement baseline stored in `live_quotes.leg_prices_json`. A devigged move past `CONSTITUENT_JUMP_THRESHOLD` cancels **every** resting quote touching that game (`constituent_jump`), not just the quote whose own leg moved. `CONSTITUENT_JUMP_MODE` selects the guard: `book_quiet` (default) additionally requires our book consensus to have stayed within `CONSTITUENT_BOOK_QUIET_MAX`, distinguishing "we are the stale ones" from "the whole market moved together"; `unconditional` drops that requirement and is the mode to flip once #54 prices every quote from a live fetch. An unreadable constituent is deliberately **not** treated as a jump — a transient API failure must not flush the resting book, and #17's confirm veto still fails closed on any resulting accept.

Also amend the existing item 2's closing sentence — it currently says the v1.1 correlation-premium gate "is documented in spec section 13 and **deferred**". Replace "and deferred" with "and is now implemented as defense item 3 (issue #23)".

- [ ] **Step 2: Add the config table rows in `README.md`**

Add next to the existing `SIGMA_Z_MAX` / `MIN_AGREEING_BOOKS` rows:

| Knob | Default | Meaning |
|---|---|---|
| `CORR_SANITY_FRECHET_ENABLED` | `true` | Reject quotes whose fair violates the Fréchet bounds implied by the live Kalshi marginals |
| `CORR_SANITY_PREMIUM_ENABLED` | `false` | Also reject on correlation premium outside the band — ships log-only until the firehose shows the real distribution |
| `CORR_PREMIUM_MIN` / `CORR_PREMIUM_MAX` | `0.5` / `2.0` | Band for `combo_fair / Π(marginals)` |
| `CONSTITUENT_JUMP_THRESHOLD` | `0.03` | Devigged move on a constituent single, since quote placement, that pulls every quote on that game |
| `CONSTITUENT_JUMP_MODE` | `book_quiet` | `book_quiet` (today) also requires quiet books; `unconditional` is the post-#54 mode |
| `CONSTITUENT_BOOK_QUIET_MAX` | `0.01` | How far book consensus may drift and still count as "quiet" |

- [ ] **Step 3: Update the maker bullet in `CLAUDE.md`**

Append to the `kalshi_mlb_mm/` bullet:

> Kalshi's own constituent single-leg markets are the bot's one real-time, book-independent anchor (issue #23): at quote time the devigged marginals gate the combo fair on Fréchet bounds (`corr_sanity_frechet`; the correlation-premium band ships log-only), reusing #17's leg snapshot for zero extra API calls; every 10s risk sweep re-polls the distinct constituent tickers of all resting quotes and cancels every quote touching a game whose constituent moved past `CONSTITUENT_JUMP_THRESHOLD` since placement (`constituent_jump`, one GET per distinct ticker per sweep). `CONSTITUENT_JUMP_MODE` flips from `book_quiet` to `unconditional` when #54 makes every quote live-priced.

- [ ] **Step 4: Run the FULL suite**

Run: `python3 -m pytest kalshi_mlb_mm/tests/ kalshi_common/tests/`
Expected: PASS — **≥1235 tests** (1235 on main + the new ones). Any failure here blocks the merge; do not proceed on a partial pass.

- [ ] **Step 5: Pre-merge executive review**

Run `git diff main..HEAD` and review against the repo's mandatory checklist — data integrity, resource safety (`db.connect` context managers, no leaked API sessions), edge cases (no open quotes, NULL `leg_prices_json`, unresolvable game, `yes_ask=100`, crossed book, cross-game combos), dead code, log/disk hygiene, secrets. Document findings as ISSUES TO FIX vs ACCEPTABLE RISKS, fix the issues, then re-run the full suite.

Specific things to verify in review, because they are the ways this change could silently hurt:
- The sweep's added SELECT column did not break another consumer of that row tuple.
- `research.emit` payloads are JSON-serialisable (`marginals` is a plain list of floats, `moves` a dict of str→float).
- The API call count per sweep really is one per distinct ticker — assert it in a test, do not eyeball it.
- No new `print()`.
- `corr_sanity` cannot divide by zero: every marginal is strictly inside (0,1) by construction in `marginals_for_legs`.

- [ ] **Step 6: Commit and request approval**

```bash
git add kalshi_mlb_mm/README.md CLAUDE.md
git commit -m "docs(mm): defense hierarchy + config for the Kalshi singles anchor (#23)"
```

Then present the diff summary and the review findings and **ask for explicit merge approval**. Do not merge. Remind the user that this touches the live quote path and the risk sweep, so a maker restart is needed after merge — and that there are already 3 merges pending a restart, with #42 as the scraper-side gate.

---

## Version control

- **Branch:** `worktree-mm-issue-23-singles-anchor` (created from local `main` @ `060a9ce`).
- **Worktree:** `.claude/worktrees/mm-issue-23-singles-anchor` — created before any file write, including this plan.
- **Commits:** one per task (4 total), each with its tests green.
- **Merge:** only after explicit user approval. Then `git worktree remove .claude/worktrees/mm-issue-23-singles-anchor` and `git branch -d worktree-mm-issue-23-singles-anchor`.
- **Push:** nothing is pushed without asking — local `main` is intentionally 3 merges ahead of `origin`.

## Documentation

- `kalshi_mlb_mm/README.md` — defense hierarchy items 3 & 4 (renumbering the rest), config table rows, the per-sweep API bound. Same merge as the code.
- `CLAUDE.md` — maker bullet.
- No new README needed; `kalshi_mlb_mm/` already has one.

---

## Self-review

**Spec coverage.** Item 1 (constituent snapshot + placement baseline) → Task 3 steps 4–5, with the ring buffer deliberately dropped per locked decision 1 and the baseline sourced from `leg_prices_json` per decision 2. Item 2 (jump breaker, both modes) → Task 3. Item 3 (Fréchet/premium gate) → Task 2. Item 4 (`singles_moved_since` fast path) → Task 1 step 5, satisfied as the shared-primitive refactor per decision 3, with the confirm path keeping its live fetch per decision 4. Item 5 (premium + marginals in the firehose) → Task 2 step 4 (`corr_sanity_check`). Acceptance criteria: jump+quiet→cancel, quiet→untouched, premium out of band→rejected, Fréchet violation→rejected, degenerate shapes→no crash, all covered; API bound documented in Task 3 step 4's docstring, README, and asserted by `test_api_calls_are_one_per_distinct_ticker_per_sweep`.

**Placeholder scan.** Task 2 step 1 and Task 3 step 1 give test *names, docstrings and assertions* but defer the harness body to "mirror `test_dispersion_gate.py` / `test_confirm_singles_veto.py`". That is a real gap — the implementer must read those two files first and copy their `_setup` shape. Everything else ships literal code.

**Type consistency.** `corr_sanity` returns `CorrSanity | None`; every call site guards on `sanity is not None` before touching `.reason`. `marginals_for_legs` returns `list[float] | None`; the call site passes it to `corr_sanity` only when truthy. `jumped_tickers` returns `dict[str, float]`, consumed as both a truthiness test and a `set(moves)` of tickers. `fetch_market_prices` returns a possibly-partial dict, and `jumped_tickers` skips absent tickers — consistent with decision 7. The sweep row tuple is 7-wide everywhere after Task 3 step 5(a).
