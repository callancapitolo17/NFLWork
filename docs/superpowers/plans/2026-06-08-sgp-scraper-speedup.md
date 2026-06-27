# MLB SGP Scraper Speedup Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Cut the SGP refresh cycle from ~160s to ≤60s for both Kalshi bots by parallelizing per-target pricing calls, replacing subprocess-per-cycle with an in-process `SGPService`, and TTL-caching structure fetches.

**Architecture:** Three layers change. (1) `mlb_sgp/draftkings.py` and `fanduel.py` gain target-level `ThreadPoolExecutor` parallelism (mirroring the existing pattern in `novig.py`/`prophetx.py`); all four orchestrators gain a `parallelism` parameter. (2) A new `kalshi_common/sgp_service.py` holds persistent per-book HTTP clients across cycles and runs the four book orchestrators concurrently with a per-book deadline; `kalshi_common/sgp_runner.py::sgp_cycle` gains a `service=` parameter selecting the in-process path (legacy subprocess path kept for rollback). (3) DK/FD orchestrators gain injected fetcher hooks so the service can TTL-cache event lists and selection-ID dictionaries; prices are never cached. The dashboard's CLI-shim path is untouched (it inherits the Phase-1 parallelism win only).

**Tech Stack:** Python 3 stdlib (`concurrent.futures`, `threading`), `curl_cffi` (existing), DuckDB, pytest.

**Spec:** `docs/superpowers/specs/2026-06-08-sgp-scraper-speedup-design.md`

---

## Conventions for every task

- **Worktree:** all edits happen in `/Users/callancapitolo/NFLWork/.claude/worktrees/sgp-scraper-speedup` (branch `worktree-sgp-scraper-speedup`). Run all commands from the worktree root.
- **Venvs are gitignored**, so the worktree has none. Use the main repo's venv *binaries* with the worktree as cwd — pytest resolves modules from cwd, so worktree code is what runs:
  - mlb_sgp tests: `/Users/callancapitolo/NFLWork/mlb_sgp/venv/bin/python -m pytest mlb_sgp/tests/... -q`
  - kalshi tests: `/Users/callancapitolo/NFLWork/kalshi_mlb_rfq/venv/bin/python -m pytest kalshi_mlb_rfq/tests/... -q`
- **Never touch live DuckDB files** (`kalshi_mlb_rfq/*.duckdb`, `Answer Keys/*.duckdb`) from the worktree except read-only. Tests use `tmp_path` DBs.
- **Commit after every task** with the message given in the task's final step.

### Known deviations from the spec (intentional, decided at planning)

1. `SGPService` lives in a **new file** `kalshi_common/sgp_service.py` (re-exported from `sgp_runner.py`), not inside `sgp_runner.py` — that file is already 392 lines with a different responsibility (Kalshi enumeration + DB I/O).
2. TTL structure-caching is implemented for **DK and FD only**. PX/NV's per-cycle structure fetches are one or two cheap calls (vs DK's 2MB-per-game payload); their big win is the persistent client, which they do get. Revisit if post-ship profiling disagrees. (Task 12 updates the spec table to match.)
3. `SGPService.refresh()` returns `dict[book → rows-or-None]` rather than a flat list — the per-book success/failure distinction is required for correct stale-row clearing (a failed book must keep its old rows, mirroring today's subprocess-crash behavior).

---

### Task 1: DraftKings target-level parallelism

**Files:**
- Modify: `mlb_sgp/draftkings.py`
- Test: `mlb_sgp/tests/test_dk_target_parallelism.py` (create)

DK currently prices its ~365 surviving targets in a sequential `for t in filtered_targets:` loop (`draftkings.py:282`) — only the 4 combos inside each target run concurrently. This task lifts the loop into a `ThreadPoolExecutor`, exactly mirroring the proven pattern in `novig.py:377`.

- [ ] **Step 1: Write the failing tests**

Create `mlb_sgp/tests/test_dk_target_parallelism.py`:

```python
"""Target-level parallelism tests for mlb_sgp/draftkings.py::price_sgps.

Strategy: monkeypatch the legacy helper functions on the real
scraper_draftkings_sgp module (price_sgps imports them lazily at call
time, so patching module attributes works), then drive price_sgps with
a many-target fixture and assert the complete row set comes back —
identical content to what the sequential loop produced, regardless of
completion order.
"""
import sys
import threading
import time
from datetime import datetime, timezone
from pathlib import Path
from unittest.mock import MagicMock

# price_sgps does `from scraper_draftkings_sgp import ...` (top-level
# module name) — only resolvable with mlb_sgp/ itself on sys.path.
MLB_SGP_DIR = Path(__file__).resolve().parents[1]
if str(MLB_SGP_DIR) not in sys.path:
    sys.path.insert(0, str(MLB_SGP_DIR))

import scraper_draftkings_sgp as legacy  # noqa: E402

from mlb_sgp import draftkings  # noqa: E402
from mlb_sgp._shared import TargetLine  # noqa: E402

CT = datetime(2026, 6, 8, 23, 0, tzinfo=timezone.utc)


def _mk_target(game_id, spread, total):
    return TargetLine(
        game_id=game_id, home_team="New York Yankees",
        away_team="Boston Red Sox", commence_time=CT,
        period="FG", spread=spread, total=total,
    )


def _wire_fixture(monkeypatch, n_targets, fail_index=None, call_sleep=0.0,
                  concurrency_tracker=None):
    """One game, n_targets distinct (spread, total) lines, all offered.

    fail_index: calculate_sgp raises for that target's total sel-ids —
    its 4 combos must vanish without affecting other targets.
    concurrency_tracker: dict with 'lock', 'live', 'max' — records the
    peak number of in-flight calculate_sgp calls.
    """
    spreads, totals, targets = {}, {}, []
    for i in range(n_targets):
        spread, total = -1.5 - i, 7.5 + i
        targets.append(_mk_target("g1", spread, total))
        spreads[("N", abs(spread), "1")] = [f"H{i}"]
        spreads[("P", abs(spread), "3")] = [f"A{i}"]
        totals[("O", total)] = [f"O{i}"]
        totals[("U", total)] = [f"U{i}"]
    sel_ids_all = {
        "fg": {"spreads": spreads, "totals": totals, "canonical": {"M"}},
        "f5": {"spreads": {}, "totals": {}},
    }
    monkeypatch.setattr(legacy, "fetch_dk_events",
                        lambda s: [{"dk_event_id": "e1"}])
    monkeypatch.setattr(legacy, "match_events",
                        lambda evs, tgt: [{"game_id": "g1", "dk_event_id": "e1"}])
    monkeypatch.setattr(legacy, "fetch_main_market_nums", lambda s, eid: {})
    monkeypatch.setattr(legacy, "fetch_selection_ids",
                        lambda s, eid, nums, v=False: sel_ids_all)
    monkeypatch.setattr(legacy, "_market_num", lambda sel: "M")

    fail_sels = {f"O{fail_index}", f"U{fail_index}"} if fail_index is not None else set()

    def fake_calculate(session, sp, to, verbose=False):
        if concurrency_tracker is not None:
            with concurrency_tracker["lock"]:
                concurrency_tracker["live"] += 1
                concurrency_tracker["max"] = max(
                    concurrency_tracker["max"], concurrency_tracker["live"])
        try:
            if call_sleep:
                time.sleep(call_sleep)
            if to in fail_sels:
                raise RuntimeError("boom")
            return {"trueOdds": 3.5, "displayOdds": "+250"}
        finally:
            if concurrency_tracker is not None:
                with concurrency_tracker["lock"]:
                    concurrency_tracker["live"] -= 1

    monkeypatch.setattr(legacy, "calculate_sgp", fake_calculate)
    return targets


def test_parallel_pricing_returns_full_row_set(monkeypatch):
    targets = _wire_fixture(monkeypatch, n_targets=12)
    rows = draftkings.price_sgps(targets, periods=("FG",),
                                 client=MagicMock(), parallelism=8)
    assert len(rows) == 12 * 4
    keys = {(r.combo, r.spread_line, r.total_line) for r in rows}
    assert ("Home Spread + Over", -1.5, 7.5) in keys
    assert ("Away Spread + Under", -12.5, 18.5) in keys
    assert all(r.sgp_decimal == 3.5 and r.sgp_american == 250 for r in rows)
    assert all(r.source == "draftkings_direct" for r in rows)
    assert all(r.bookmaker == "draftkings" for r in rows)


def test_one_target_failure_does_not_drop_others(monkeypatch):
    targets = _wire_fixture(monkeypatch, n_targets=6, fail_index=3)
    rows = draftkings.price_sgps(targets, periods=("FG",),
                                 client=MagicMock(), parallelism=4)
    # target 3's four combos vanish; the other 5 targets price fully
    assert len(rows) == 5 * 4
    assert not any(r.total_line == 7.5 + 3 for r in rows)


def test_targets_actually_run_concurrently(monkeypatch):
    tracker = {"lock": threading.Lock(), "live": 0, "max": 0}
    targets = _wire_fixture(monkeypatch, n_targets=8, call_sleep=0.02,
                            concurrency_tracker=tracker)
    rows = draftkings.price_sgps(targets, periods=("FG",),
                                 client=MagicMock(), parallelism=8)
    assert len(rows) == 8 * 4
    # 8 targets x 4 combos with 20ms calls: sequential targets would
    # never exceed 4 in flight. Parallel must.
    assert tracker["max"] > 4


def test_parallelism_default_reads_env(monkeypatch):
    monkeypatch.setenv("MLB_SGP_DK_PARALLELISM", "3")
    assert draftkings._resolve_parallelism(None) == 3
    assert draftkings._resolve_parallelism(7) == 7
    monkeypatch.delenv("MLB_SGP_DK_PARALLELISM")
    assert draftkings._resolve_parallelism(None) == 8  # shipped default
```

- [ ] **Step 2: Run the tests to verify they fail**

Run: `/Users/callancapitolo/NFLWork/mlb_sgp/venv/bin/python -m pytest mlb_sgp/tests/test_dk_target_parallelism.py -q`
Expected: FAIL — `TypeError: price_sgps() got an unexpected keyword argument 'parallelism'` (and `AttributeError` for `_resolve_parallelism`).

- [ ] **Step 3: Implement target-level parallelism in `draftkings.py`**

3a. Add `import os` to the imports block (after `from concurrent.futures import ...`).

3b. Add below the `SOURCE_LABEL_FALLBACK` constant:

```python
# Target-level parallelism. Env-overridable for ops tuning without a
# code edit; the shipped default comes from the Phase-0 probe
# (mlb_sgp/probe_concurrency.py). Total in-flight DK requests is
# parallelism x 4 (the per-target combo pool nests inside).
DK_TARGET_PARALLELISM_DEFAULT = 8


def _resolve_parallelism(parallelism: int | None) -> int:
    if parallelism is not None:
        return parallelism
    return int(os.environ.get("MLB_SGP_DK_PARALLELISM",
                              str(DK_TARGET_PARALLELISM_DEFAULT)))
```

3c. Change the `price_sgps` signature:

```python
def price_sgps(
    target_lines: list[TargetLine],
    periods: tuple[str, ...] = ("FG",),
    client: DraftKingsClient | None = None,
    verbose: bool = False,
    parallelism: int | None = None,
) -> list[PricedRow]:
```

3d. Replace the entire Phase-2 block — from the comment
`# ----- Phase 2: per target row, build and price 4 combos ----- #`
(line ~278) through the end of the `for t in filtered_targets:` loop
(the line before `return out`) — with a closure + pool. The closure body
is the existing loop body verbatim, with two mechanical changes:
`continue` → `return target_rows`, and `out.append(` → `target_rows.append(`:

```python
    # ----- Phase 2: per target row, build and price 4 combos ----- #
    # Each TargetLine is priced independently (its own sel-id lookups +
    # 4 calculateBets calls), so targets fan out on a thread pool. The
    # per-target combo pool (4-wide) nests inside, giving parallelism x 4
    # in-flight requests — the granularity the Phase-0 probe validates.
    # Pattern mirrors novig.py / prophetx.py (in production since the
    # PX/NV integration shipped).
    def _price_one_target(t: TargetLine) -> list[PricedRow]:
        target_rows: list[PricedRow] = []
        cache = per_game_cache.get(t.game_id)
        if cache is None:
            return target_rows

        period_key = t.period.lower()  # "FG" -> "fg" for legacy dict keys
        sel_ids_all = cache["sel_ids_all"]
        sel_ids = sel_ids_all.get(period_key, {"spreads": {}, "totals": {}})
        canonical = sel_ids.get("canonical", set())

        # Sign convention: home favored (spread < 0) means home is "N"
        # (negative-line side of the spread market), away is "P". Flip
        # when home is the underdog. ``spread`` itself is stored as
        # ``abs(t.spread)`` because DK selection IDs encode unsigned
        # magnitudes (the sign letter carries direction).
        if t.spread < 0:
            home_sign, away_sign = "N", "P"
        else:
            home_sign, away_sign = "P", "N"
        spread_abs = abs(t.spread)

        # _1 = home participant, _3 = away participant in DK sel_id suffix.
        home_spread_sels = sel_ids["spreads"].get((home_sign, spread_abs, "1")) or []
        away_spread_sels = sel_ids["spreads"].get((away_sign, spread_abs, "3")) or []
        over_sels = sel_ids["totals"].get(("O", t.total)) or []
        under_sels = sel_ids["totals"].get(("U", t.total)) or []

        # If any of the four legs is missing, the requested (spread,
        # total) is off-main. Fall back to integer-line derivation.
        if not (home_spread_sels and away_spread_sels and over_sels and under_sels):
            fallback = try_integer_fallback_dk(
                client.session, sel_ids, t.spread, t.total, canonical,
                verbose=verbose,
            )
            if fallback is None:
                return target_rows
            # 4 derived combos at once, all tagged _interpolated
            prefix = "" if t.period == "FG" else "F5 "
            fair_probs = fallback["fair_probs"]
            for key, base_combo in _FALLBACK_KEY_TO_COMBO.items():
                fair_p = fair_probs.get(key)
                if not fair_p or fair_p <= 0:
                    continue
                dec = 1.0 / fair_p
                target_rows.append(PricedRow(
                    game_id=t.game_id,
                    combo=prefix + base_combo,
                    period=t.period,
                    spread_line=t.spread,
                    total_line=t.total,
                    bookmaker=BOOK_NAME,
                    source=SOURCE_LABEL_FALLBACK,
                    sgp_decimal=round(dec, 4),
                    sgp_american=decimal_to_american(dec),
                    fetch_time=fetch_now,
                ))
            return target_rows

        # ----- Main path: price the 4 canonical combos in parallel ----- #
        prefix = "" if t.period == "FG" else "F5 "
        combos = (
            ("Home Spread + Over",  home_spread_sels, over_sels),
            ("Home Spread + Under", home_spread_sels, under_sels),
            ("Away Spread + Over",  away_spread_sels, over_sels),
            ("Away Spread + Under", away_spread_sels, under_sels),
        )

        priced_by_combo = _price_combos_parallel(
            client.session, combos, canonical,
            calculate_sgp, _market_num, verbose,
        )

        for combo_name, _sp_sels, _tot_sels in combos:
            sgp = priced_by_combo.get(combo_name)
            if not sgp:
                continue
            dec = sgp["trueOdds"]
            target_rows.append(PricedRow(
                game_id=t.game_id,
                combo=prefix + combo_name,
                period=t.period,
                spread_line=t.spread,
                total_line=t.total,
                bookmaker=BOOK_NAME,
                source=SOURCE_LABEL,
                sgp_decimal=round(dec, 4),
                sgp_american=decimal_to_american(dec),
                fetch_time=fetch_now,
            ))
        return target_rows

    n_workers = max(1, _resolve_parallelism(parallelism))
    with ThreadPoolExecutor(max_workers=n_workers) as pool:
        futures = [pool.submit(_price_one_target, t) for t in filtered_targets]
        for f in as_completed(futures):
            try:
                out.extend(f.result())
            except Exception as e:
                if verbose:
                    print(f"  dk target error: {e}", flush=True)

    return out
```

- [ ] **Step 4: Run the new tests + the existing DK suites**

Run: `/Users/callancapitolo/NFLWork/mlb_sgp/venv/bin/python -m pytest mlb_sgp/tests/test_dk_target_parallelism.py mlb_sgp/tests/test_draftkings_orchestrator.py -q`
Expected: ALL PASS.

- [ ] **Step 5: Commit**

```bash
git add mlb_sgp/draftkings.py mlb_sgp/tests/test_dk_target_parallelism.py
git commit -m "perf(sgp): DK target-level parallelism (env-tunable, default 8)

DK priced ~365 targets sequentially (~400ms each = ~146s/cycle, the
cycle's long pole). Targets now fan out on a thread pool mirroring the
proven novig.py/prophetx.py pattern; per-target combo pool unchanged."
```

---

### Task 2: FanDuel target-level parallelism

**Files:**
- Modify: `mlb_sgp/fanduel.py`
- Test: `mlb_sgp/tests/test_fd_target_parallelism.py` (create)

Same refactor as Task 1, FD shapes. FD is not the long pole (strict line
matching means few targets survive), so its default stays conservative.

- [ ] **Step 1: Write the failing tests**

Create `mlb_sgp/tests/test_fd_target_parallelism.py`:

```python
"""Target-level parallelism tests for mlb_sgp/fanduel.py::price_sgps."""
import sys
from datetime import datetime, timezone
from pathlib import Path
from unittest.mock import MagicMock

MLB_SGP_DIR = Path(__file__).resolve().parents[1]
if str(MLB_SGP_DIR) not in sys.path:
    sys.path.insert(0, str(MLB_SGP_DIR))

import scraper_fanduel_sgp as legacy  # noqa: E402

from mlb_sgp import fanduel  # noqa: E402
from mlb_sgp._shared import TargetLine  # noqa: E402

CT = datetime(2026, 6, 8, 23, 0, tzinfo=timezone.utc)


def _mk_target(game_id, spread, total):
    return TargetLine(
        game_id=game_id, home_team="New York Yankees",
        away_team="Boston Red Sox", commence_time=CT,
        period="FG", spread=spread, total=total,
    )


def _wire_fixture(monkeypatch, n_targets, fail_index=None):
    spreads, totals, targets = {}, {}, []
    for i in range(n_targets):
        spread, total = -1.5 - i, 7.5 + i
        targets.append(_mk_target("g1", spread, total))
        spreads[("home", spread)] = (f"SM{i}", f"SH{i}")
        spreads[("away", -spread)] = (f"SM{i}", f"SA{i}")
        totals[("O", total)] = (f"TM{i}", f"TO{i}")
        totals[("U", total)] = (f"TM{i}", f"TU{i}")
    runners = {"fg": {"spreads": spreads, "totals": totals},
               "f5": {"spreads": {}, "totals": {}}}
    monkeypatch.setattr(legacy, "fetch_fd_events", lambda s: [{"id": "e1"}])
    monkeypatch.setattr(legacy, "match_events", lambda evs, tgt: [
        {"game_id": "g1", "fd_event_id": "e1",
         "fd_home": "New York Yankees", "fd_away": "Boston Red Sox"}])
    monkeypatch.setattr(legacy, "fetch_event_runners",
                        lambda s, eid, h, a: runners)

    fail_sids = {f"TO{fail_index}", f"TU{fail_index}"} if fail_index is not None else set()

    def fake_price_combo(session, smid, ssid, tmid, tsid, verbose=False):
        if tsid in fail_sids:
            raise RuntimeError("boom")
        return {"decimal": 3.4, "american": 240}

    monkeypatch.setattr(legacy, "price_combo", fake_price_combo)
    return targets


def test_parallel_pricing_returns_full_row_set(monkeypatch):
    targets = _wire_fixture(monkeypatch, n_targets=10)
    rows = fanduel.price_sgps(targets, periods=("FG",),
                              client=MagicMock(), parallelism=6)
    assert len(rows) == 10 * 4
    assert all(r.sgp_decimal == 3.4 and r.sgp_american == 240 for r in rows)
    assert all(r.source == "fanduel_direct" for r in rows)


def test_one_target_failure_does_not_drop_others(monkeypatch):
    targets = _wire_fixture(monkeypatch, n_targets=5, fail_index=2)
    rows = fanduel.price_sgps(targets, periods=("FG",),
                              client=MagicMock(), parallelism=4)
    assert len(rows) == 4 * 4
    assert not any(r.total_line == 7.5 + 2 for r in rows)


def test_parallelism_default_reads_env(monkeypatch):
    monkeypatch.setenv("MLB_SGP_FD_PARALLELISM", "2")
    assert fanduel._resolve_parallelism(None) == 2
    assert fanduel._resolve_parallelism(9) == 9
    monkeypatch.delenv("MLB_SGP_FD_PARALLELISM")
    assert fanduel._resolve_parallelism(None) == 4  # conservative default
```

- [ ] **Step 2: Run to verify failure**

Run: `/Users/callancapitolo/NFLWork/mlb_sgp/venv/bin/python -m pytest mlb_sgp/tests/test_fd_target_parallelism.py -q`
Expected: FAIL — unexpected keyword `parallelism` / missing `_resolve_parallelism`.

- [ ] **Step 3: Implement in `fanduel.py`**

3a. Add `import os` to imports.

3b. Below `_FALLBACK_KEY_TO_COMBO`, add:

```python
# Target-level parallelism. FD is not the cycle's long pole (strict
# line matching filters most targets), so the default is conservative;
# env-overridable for ops tuning.
FD_TARGET_PARALLELISM_DEFAULT = 4


def _resolve_parallelism(parallelism: int | None) -> int:
    if parallelism is not None:
        return parallelism
    return int(os.environ.get("MLB_SGP_FD_PARALLELISM",
                              str(FD_TARGET_PARALLELISM_DEFAULT)))
```

3c. Add `parallelism: int | None = None` to the `price_sgps` signature
(after `verbose`).

3d. Replace the Phase-2 block (`# ----- Phase 2: per target row, price 4
combos ----- #` through the loop end, before `return out`) with the
closure + pool — same mechanical transformation as Task 1
(`continue` → `return target_rows`, `out.append(` → `target_rows.append(`),
loop body otherwise verbatim from the current file:

```python
    # ----- Phase 2: per target row, price 4 combos ----- #
    # Targets fan out on a thread pool (mirrors draftkings.py /
    # novig.py). Runners are pre-fetched in the filter loop above, so
    # workers only do implyBets calls — no shared-cache mutation races.
    def _price_one_target(t: TargetLine) -> list[PricedRow]:
        target_rows: list[PricedRow] = []
        game = matched_by_gid.get(t.game_id)
        if game is None:
            return target_rows

        period_key = t.period.lower()  # "FG" -> "fg"

        # Runners cache is guaranteed-populated by the filter loop
        # above (we only added a target if its game's runners were
        # fetched successfully).
        sel_ids_per_period = runners_cache[t.game_id]
        sel = sel_ids_per_period.get(period_key, {"spreads": {}, "totals": {}})

        # Sign convention: TargetLine.spread is the home-perspective
        # signed line, so the away leg is keyed by -spread.
        if not sel.get("spreads"):
            return target_rows
        home_line = t.spread
        away_line = -t.spread

        home_spread = sel["spreads"].get(("home", home_line))
        away_spread = sel["spreads"].get(("away", away_line))
        over = sel["totals"].get(("O", t.total))
        under = sel["totals"].get(("U", t.total))

        # If any of the four legs is missing, the requested (spread,
        # total) is off-main / off-alt. Try the integer-line fallback.
        if not (home_spread and away_spread and over and under):
            fallback = try_integer_fallback_fd(
                client.session, sel, home_line, away_line, t.total,
                verbose=verbose,
            )
            if fallback is None:
                return target_rows
            prefix = "" if t.period == "FG" else "F5 "
            fair_probs = fallback["fair_probs"]
            for key, base_combo in _FALLBACK_KEY_TO_COMBO.items():
                fair_p = fair_probs.get(key)
                if not fair_p or fair_p <= 0:
                    continue
                dec = 1.0 / fair_p
                target_rows.append(PricedRow(
                    game_id=t.game_id,
                    combo=prefix + base_combo,
                    period=t.period,
                    spread_line=t.spread,
                    total_line=t.total,
                    bookmaker=BOOK_NAME,
                    source=SOURCE_LABEL_FALLBACK,
                    sgp_decimal=round(dec, 4),
                    sgp_american=decimal_to_american(dec),
                    fetch_time=fetch_now,
                ))
            return target_rows

        # ----- Main path: price the 4 canonical combos in parallel ----- #
        prefix = "" if t.period == "FG" else "F5 "
        combos = (
            ("Home Spread + Over",  home_spread, over),
            ("Home Spread + Under", home_spread, under),
            ("Away Spread + Over",  away_spread, over),
            ("Away Spread + Under", away_spread, under),
        )

        priced_by_combo = _price_combos_parallel(
            client.session, combos, price_combo, verbose,
        )

        for combo_name, _sp_pair, _tot_pair in combos:
            result = priced_by_combo.get(combo_name)
            if not result:
                continue
            dec = float(result["decimal"])
            am = int(result["american"])
            target_rows.append(PricedRow(
                game_id=t.game_id,
                combo=prefix + combo_name,
                period=t.period,
                spread_line=t.spread,
                total_line=t.total,
                bookmaker=BOOK_NAME,
                source=SOURCE_LABEL,
                sgp_decimal=round(dec, 4),
                sgp_american=am,
                fetch_time=fetch_now,
            ))
        return target_rows

    n_workers = max(1, _resolve_parallelism(parallelism))
    with ThreadPoolExecutor(max_workers=n_workers) as pool:
        futures = [pool.submit(_price_one_target, t) for t in filtered_targets]
        for f in as_completed(futures):
            try:
                out.extend(f.result())
            except Exception as e:
                if verbose:
                    print(f"  fd target error: {e}", flush=True)

    return out
```

- [ ] **Step 4: Run new + existing FD tests**

Run: `/Users/callancapitolo/NFLWork/mlb_sgp/venv/bin/python -m pytest mlb_sgp/tests/test_fd_target_parallelism.py mlb_sgp/tests/test_fanduel_orchestrator.py -q`
Expected: ALL PASS.

- [ ] **Step 5: Commit**

```bash
git add mlb_sgp/fanduel.py mlb_sgp/tests/test_fd_target_parallelism.py
git commit -m "perf(sgp): FD target-level parallelism (env-tunable, default 4)"
```

---

### Task 3: ProphetX parallelism made tunable

**Files:**
- Modify: `mlb_sgp/prophetx.py` (constant at line ~81, signature, pool at line ~421)
- Test: append to `mlb_sgp/tests/test_prophetx_orchestrator.py`

PX already fans targets out 2-wide. This task only makes the width
env-overridable and call-time-injectable so the probe (Task 4) and the
post-probe tuning (Task 13) don't need code edits. **Default stays 2
until the probe says otherwise** — PX RFQs are a real market footprint.

- [ ] **Step 1: Write the failing test** (append to `test_prophetx_orchestrator.py`)

```python
def test_px_parallelism_resolution(monkeypatch):
    """PX pool width: explicit arg > env > module default (2)."""
    from mlb_sgp import prophetx
    monkeypatch.setenv("MLB_SGP_PX_PARALLELISM", "5")
    assert prophetx._resolve_parallelism(None) == 5
    assert prophetx._resolve_parallelism(3) == 3
    monkeypatch.delenv("MLB_SGP_PX_PARALLELISM")
    assert prophetx._resolve_parallelism(None) == 2


def test_px_price_sgps_accepts_parallelism_kwarg():
    """Contract: the kwarg exists and the empty-input early-exit honors it."""
    from mlb_sgp import prophetx
    assert prophetx.price_sgps([], parallelism=4) == []
```

- [ ] **Step 2: Run to verify failure**

Run: `/Users/callancapitolo/NFLWork/mlb_sgp/venv/bin/python -m pytest mlb_sgp/tests/test_prophetx_orchestrator.py -q -k parallelism`
Expected: FAIL — `_resolve_parallelism` missing.

- [ ] **Step 3: Implement**

In `prophetx.py`: ensure `import os` is present; replace
`PX_TARGET_PARALLELISM = 2` with:

```python
# PX RFQs are a real market footprint — keep the default at 2 until the
# Phase-0 probe (mlb_sgp/probe_concurrency.py) justifies more. Env-
# overridable so tuning needs no code edit.
PX_TARGET_PARALLELISM_DEFAULT = 2


def _resolve_parallelism(parallelism: int | None) -> int:
    if parallelism is not None:
        return parallelism
    return int(os.environ.get("MLB_SGP_PX_PARALLELISM",
                              str(PX_TARGET_PARALLELISM_DEFAULT)))
```

Add `parallelism: int | None = None` to `price_sgps` (after `verbose`),
and change the pool line `with ThreadPoolExecutor(max_workers=PX_TARGET_PARALLELISM) as pool:`
to:

```python
    with ThreadPoolExecutor(max_workers=max(1, _resolve_parallelism(parallelism))) as pool:
```

If any other code references the old `PX_TARGET_PARALLELISM` name, run
`grep -rn "PX_TARGET_PARALLELISM" mlb_sgp/ kalshi_mlb_rfq/ kalshi_mlb_mm/ kalshi_common/`
and update those references to `PX_TARGET_PARALLELISM_DEFAULT`.

- [ ] **Step 4: Run the full PX test file**

Run: `/Users/callancapitolo/NFLWork/mlb_sgp/venv/bin/python -m pytest mlb_sgp/tests/test_prophetx_orchestrator.py -q`
Expected: ALL PASS.

- [ ] **Step 5: Commit**

```bash
git add mlb_sgp/prophetx.py mlb_sgp/tests/test_prophetx_orchestrator.py
git commit -m "perf(sgp): make PX target parallelism env-tunable (default stays 2 pre-probe)"
```

---

### Task 4: Concurrency probe harness

**Files:**
- Create: `mlb_sgp/probe_concurrency.py`
- Test: `mlb_sgp/tests/test_probe_concurrency.py` (create)

The probe answers "how wide can we go per book?" by timing **real
`price_sgps` runs** on a small fixed target set at escalating
parallelism, measuring rows/sec and miss-rate per level. Reusing the
orchestrators (instead of raw request loops) means the probe exercises
exactly the request mix production will use. Budgets are hard-coded per
the spec: ≤150 DK calls, ≤40 PX RFQs, cooldown between levels, post-run
health check. **The script is manual-only — never invoked by the bots.**

- [ ] **Step 1: Write the failing tests for the pure helpers**

Create `mlb_sgp/tests/test_probe_concurrency.py`:

```python
"""Unit tests for the pure planning/summary helpers in probe_concurrency."""
from datetime import datetime, timezone

from mlb_sgp._shared import TargetLine
from mlb_sgp.probe_concurrency import (
    BOOK_PLANS,
    pick_probe_targets,
    summarize_level,
)

CT = datetime(2026, 6, 8, 23, 0, tzinfo=timezone.utc)


def _t(gid, spread=-1.5, total=8.5):
    return TargetLine(game_id=gid, home_team="H", away_team="A",
                      commence_time=CT, period="FG", spread=spread, total=total)


def test_pick_probe_targets_one_per_game_capped():
    targets = [_t("g1"), _t("g1", -2.5), _t("g2"), _t("g3"), _t("g3", -3.5)]
    picked = pick_probe_targets(targets, n_games=2)
    assert len(picked) == 2
    assert [p.game_id for p in picked] == ["g1", "g2"]


def test_pick_probe_targets_prefers_main_spread():
    targets = [_t("g1", -4.5), _t("g1", -1.5), _t("g1", -2.5)]
    picked = pick_probe_targets(targets, n_games=1)
    assert picked[0].spread == -1.5


def test_summarize_level_computes_rates():
    s = summarize_level(level=8, n_targets=6, n_rows=20, elapsed_sec=4.0)
    assert s["level"] == 8
    assert s["expected_rows"] == 24
    assert s["rows_per_sec"] == 5.0
    assert abs(s["miss_pct"] - (4 / 24 * 100)) < 1e-9


def test_book_plans_respect_budgets():
    # call budget = sum over levels of (targets_per_level * 4 combos)
    # + health check (4 calls). Spec: DK <= 150, PX <= 40.
    for book, plan in BOOK_PLANS.items():
        calls = sum(plan["targets_per_level"] * 4 for _ in plan["levels"]) + 4
        assert calls <= plan["budget"], f"{book} plan exceeds budget"
```

- [ ] **Step 2: Run to verify failure**

Run: `/Users/callancapitolo/NFLWork/mlb_sgp/venv/bin/python -m pytest mlb_sgp/tests/test_probe_concurrency.py -q`
Expected: FAIL — `ModuleNotFoundError: mlb_sgp.probe_concurrency`.

- [ ] **Step 3: Write `mlb_sgp/probe_concurrency.py`**

```python
"""Concurrency-ceiling probe for the MLB SGP scrapers.

Times real ``price_sgps`` runs on a small fixed target set at
escalating ``parallelism`` levels and reports rows/sec + miss-rate per
level. The knee (throughput stops scaling, or misses spike) is the
book's empirical ceiling; ship the last clean level.

MANUAL-ONLY tool — never invoked by the bots. Budgets per the speedup
spec (docs/superpowers/specs/2026-06-08-sgp-scraper-speedup-design.md):
DK <= 150 pricing calls, PX <= 40 RFQs, cooldown between levels, and a
post-run health check at parallelism 1. Run in the morning window
(~7-8am PT, no MLB games in progress).

Usage (from repo root):
  mlb_sgp/venv/bin/python -m mlb_sgp.probe_concurrency --book dk
  mlb_sgp/venv/bin/python -m mlb_sgp.probe_concurrency --book px \
      --db kalshi_mlb_rfq/kalshi_mlb_rfq_market.duckdb
"""
from __future__ import annotations
import argparse
import sys
import time
from pathlib import Path

_REPO_ROOT = Path(__file__).resolve().parents[1]
_MLB_SGP_DIR = Path(__file__).resolve().parent
for p in (str(_REPO_ROOT), str(_MLB_SGP_DIR)):
    if p not in sys.path:
        sys.path.insert(0, p)

from mlb_sgp._shared import TargetLine, load_target_lines  # noqa: E402

DEFAULT_DB = str(_REPO_ROOT / "kalshi_mlb_rfq" / "kalshi_mlb_rfq_market.duckdb")
COOLDOWN_SEC = 20

# Per-book ramp plans. calls-per-level = targets_per_level * 4 combos.
#   DK: 5 levels x 6 targets x 4 = 120 calls + 4 health = 124 <= 150
#   PX: 4 levels x 2 targets x 4 =  32 calls + 4 health =  36 <= 40
BOOK_PLANS = {
    "dk": {"levels": [2, 4, 8, 12, 16], "targets_per_level": 6, "budget": 150},
    "px": {"levels": [2, 3, 4, 6],      "targets_per_level": 2, "budget": 40},
}


def pick_probe_targets(targets: list[TargetLine], n_games: int) -> list[TargetLine]:
    """One target per game (first n_games games, input order), preferring
    the main-ish spread (smallest |spread|) so the book always offers it."""
    by_game: dict[str, list[TargetLine]] = {}
    order: list[str] = []
    for t in targets:
        if t.game_id not in by_game:
            order.append(t.game_id)
        by_game.setdefault(t.game_id, []).append(t)
    picked = []
    for gid in order[:n_games]:
        picked.append(min(by_game[gid], key=lambda t: abs(t.spread)))
    return picked


def summarize_level(level: int, n_targets: int, n_rows: int,
                    elapsed_sec: float) -> dict:
    expected = n_targets * 4
    return {
        "level": level,
        "n_targets": n_targets,
        "n_rows": n_rows,
        "expected_rows": expected,
        "elapsed_sec": round(elapsed_sec, 2),
        "rows_per_sec": round(n_rows / elapsed_sec, 2) if elapsed_sec > 0 else 0.0,
        "miss_pct": round((expected - n_rows) / expected * 100, 1) if expected else 0.0,
    }


def _book_module(book: str):
    if book == "dk":
        from mlb_sgp import draftkings
        return draftkings
    if book == "px":
        from mlb_sgp import prophetx
        return prophetx
    raise SystemExit(f"unknown book {book!r} (use dk|px)")


def _build_client(book: str):
    if book == "dk":
        from mlb_sgp.dk_client import DraftKingsClient
        return DraftKingsClient()
    from mlb_sgp.prophetx_client import ProphetXClient
    return ProphetXClient()


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--book", required=True, choices=("dk", "px"))
    ap.add_argument("--db", default=DEFAULT_DB,
                    help="DuckDB holding mlb_target_lines (read-only)")
    args = ap.parse_args()

    plan = BOOK_PLANS[args.book]
    all_targets = [t for t in load_target_lines(args.db) if t.period == "FG"]
    if not all_targets:
        print(f"no FG target lines in {args.db} — run while the bot has "
              f"populated mlb_target_lines (or pass --db).")
        return 2
    probe_targets = pick_probe_targets(all_targets, plan["targets_per_level"])
    print(f"book={args.book} probe_targets={len(probe_targets)} "
          f"games={[t.game_id[:8] for t in probe_targets]}")

    mod = _book_module(args.book)
    client = _build_client(args.book)  # ONE persistent session for the whole ramp

    results = []
    for level in plan["levels"]:
        t0 = time.monotonic()
        rows = mod.price_sgps(probe_targets, periods=("FG",),
                              client=client, parallelism=level)
        s = summarize_level(level, len(probe_targets), len(rows),
                            time.monotonic() - t0)
        results.append(s)
        print(f"  level={s['level']:<3} rows={s['n_rows']}/{s['expected_rows']} "
              f"elapsed={s['elapsed_sec']}s rows/sec={s['rows_per_sec']} "
              f"miss={s['miss_pct']}%")
        # Stop the ramp on degradation: any misses at this level when the
        # previous level was clean, or throughput regressing >20%.
        if len(results) >= 2:
            prev = results[-2]
            deg_miss = s["miss_pct"] > prev["miss_pct"]
            deg_tput = s["rows_per_sec"] < prev["rows_per_sec"] * 0.8
            if deg_miss or deg_tput:
                print(f"  DEGRADATION at level {level} "
                      f"(miss {prev['miss_pct']}->{s['miss_pct']}%, "
                      f"tput {prev['rows_per_sec']}->{s['rows_per_sec']}). "
                      f"Stopping ramp.")
                break
        time.sleep(COOLDOWN_SEC)

    # Post-probe health check: one target at parallelism 1 must price.
    hc = mod.price_sgps(probe_targets[:1], periods=("FG",),
                        client=client, parallelism=1)
    print(f"health check: {len(hc)}/4 rows "
          f"{'OK' if hc else '** FAILED — book may have tripped a defense **'}")

    clean = [r for r in results if r["miss_pct"] == results[0]["miss_pct"]]
    if clean:
        print(f"recommended parallelism: {clean[-1]['level']} "
              f"(last level matching baseline miss-rate)")
    return 0


if __name__ == "__main__":
    sys.exit(main())
```

- [ ] **Step 4: Run the unit tests**

Run: `/Users/callancapitolo/NFLWork/mlb_sgp/venv/bin/python -m pytest mlb_sgp/tests/test_probe_concurrency.py -q`
Expected: ALL PASS. (The network `main()` is exercised live in Task 13, not here.)

- [ ] **Step 5: Commit**

```bash
git add mlb_sgp/probe_concurrency.py mlb_sgp/tests/test_probe_concurrency.py
git commit -m "feat(sgp): concurrency-ceiling probe harness (manual-only, budgeted)

Times real price_sgps runs at escalating parallelism; rows/sec +
miss-rate per level; DK<=150-call / PX<=40-RFQ budgets, cooldowns,
post-run health check per the speedup spec."
```

---

### Task 5: `TTLCache` utility

**Files:**
- Modify: `mlb_sgp/_shared.py` (append)
- Test: append to `mlb_sgp/tests/test_shared.py`

- [ ] **Step 1: Write the failing tests** (append to `mlb_sgp/tests/test_shared.py`)

```python
def test_ttl_cache_caches_within_ttl_and_expires():
    from mlb_sgp._shared import TTLCache
    clock = {"t": 100.0}
    calls = {"n": 0}

    def fetch():
        calls["n"] += 1
        return f"v{calls['n']}"

    cache = TTLCache(ttl_sec=60, now_fn=lambda: clock["t"])
    assert cache.get_or_fetch("k", fetch) == "v1"
    clock["t"] = 130.0                      # 30s later: cached
    assert cache.get_or_fetch("k", fetch) == "v1"
    assert calls["n"] == 1
    clock["t"] = 161.0                      # 61s after store: expired
    assert cache.get_or_fetch("k", fetch) == "v2"
    assert calls["n"] == 2


def test_ttl_cache_keys_are_independent_and_clear_works():
    from mlb_sgp._shared import TTLCache
    cache = TTLCache(ttl_sec=60, now_fn=lambda: 0.0)
    assert cache.get_or_fetch("a", lambda: 1) == 1
    assert cache.get_or_fetch("b", lambda: 2) == 2
    assert cache.get_or_fetch("a", lambda: 99) == 1   # still cached
    cache.clear()
    assert cache.get_or_fetch("a", lambda: 99) == 99
```

- [ ] **Step 2: Run to verify failure**

Run: `/Users/callancapitolo/NFLWork/mlb_sgp/venv/bin/python -m pytest mlb_sgp/tests/test_shared.py -q -k ttl`
Expected: FAIL — `ImportError: cannot import name 'TTLCache'`.

- [ ] **Step 3: Implement** (append to `mlb_sgp/_shared.py`; add `import time` and `import threading` to its imports)

```python
class TTLCache:
    """Tiny per-key TTL cache for structure fetches (event lists,
    selection-id dictionaries). NEVER cache prices with this.

    Thread-safe for the lock around store access; concurrent misses on
    the same key may both fetch (last write wins) — acceptable for
    idempotent GETs, and in practice each book's structure fetches run
    single-threaded (the hoisting phase of price_sgps).
    """

    def __init__(self, ttl_sec: float, now_fn=time.monotonic):
        self.ttl_sec = ttl_sec
        self._now = now_fn
        self._lock = threading.Lock()
        self._store: dict = {}

    def get_or_fetch(self, key, fetch_fn):
        with self._lock:
            ent = self._store.get(key)
            if ent is not None and (self._now() - ent[0]) < self.ttl_sec:
                return ent[1]
        val = fetch_fn()
        with self._lock:
            self._store[key] = (self._now(), val)
        return val

    def clear(self):
        with self._lock:
            self._store.clear()
```

- [ ] **Step 4: Run the full shared suite**

Run: `/Users/callancapitolo/NFLWork/mlb_sgp/venv/bin/python -m pytest mlb_sgp/tests/test_shared.py -q`
Expected: ALL PASS.

- [ ] **Step 5: Commit**

```bash
git add mlb_sgp/_shared.py mlb_sgp/tests/test_shared.py
git commit -m "feat(sgp): TTLCache utility for structure fetches"
```

---

### Task 6: Fetcher hooks in DK + FD orchestrators

**Files:**
- Modify: `mlb_sgp/draftkings.py`, `mlb_sgp/fanduel.py`
- Test: append to `mlb_sgp/tests/test_dk_target_parallelism.py` and `test_fd_target_parallelism.py`

`price_sgps` gains an optional `fetchers` dict overriding the structure
fetches. Default behavior (no dict) is byte-identical to today — the
dashboard shims never pass it. `SGPService` (Task 7) passes cache-wrapped
versions.

- [ ] **Step 1: Write the failing tests**

Append to `mlb_sgp/tests/test_dk_target_parallelism.py`:

```python
def test_fetcher_hooks_override_structure_fetches(monkeypatch):
    """When `fetchers` is passed, price_sgps uses it instead of the
    legacy module functions — the seam SGPService caches through."""
    targets = _wire_fixture(monkeypatch, n_targets=2)
    # Sabotage the legacy fns: if the hooks are ignored, these blow up.
    monkeypatch.setattr(legacy, "fetch_dk_events",
                        lambda s: (_ for _ in ()).throw(AssertionError("hook bypassed")))
    monkeypatch.setattr(legacy, "fetch_main_market_nums",
                        lambda s, e: (_ for _ in ()).throw(AssertionError("hook bypassed")))
    monkeypatch.setattr(legacy, "fetch_selection_ids",
                        lambda s, e, n, v=False: (_ for _ in ()).throw(AssertionError("hook bypassed")))

    spreads = {("N", 1.5, "1"): ["H0"], ("P", 1.5, "3"): ["A0"],
               ("N", 2.5, "1"): ["H1"], ("P", 2.5, "3"): ["A1"]}
    totals = {("O", 7.5): ["O0"], ("U", 7.5): ["U0"],
              ("O", 8.5): ["O1"], ("U", 8.5): ["U1"]}
    sel_ids_all = {"fg": {"spreads": spreads, "totals": totals, "canonical": {"M"}},
                   "f5": {"spreads": {}, "totals": {}}}
    fetchers = {
        "fetch_dk_events": lambda session: [{"dk_event_id": "e1"}],
        "fetch_main_market_nums": lambda session, eid: {},
        "fetch_selection_ids": lambda session, eid, nums, verbose=False: sel_ids_all,
    }
    rows = draftkings.price_sgps(targets, periods=("FG",), client=MagicMock(),
                                 parallelism=2, fetchers=fetchers)
    assert len(rows) == 2 * 4
```

Append to `mlb_sgp/tests/test_fd_target_parallelism.py`:

```python
def test_fetcher_hooks_override_structure_fetches(monkeypatch):
    targets = _wire_fixture(monkeypatch, n_targets=2)
    monkeypatch.setattr(legacy, "fetch_fd_events",
                        lambda s: (_ for _ in ()).throw(AssertionError("hook bypassed")))
    monkeypatch.setattr(legacy, "fetch_event_runners",
                        lambda s, e, h, a: (_ for _ in ()).throw(AssertionError("hook bypassed")))

    spreads = {("home", -1.5): ("SM0", "SH0"), ("away", 1.5): ("SM0", "SA0"),
               ("home", -2.5): ("SM1", "SH1"), ("away", 2.5): ("SM1", "SA1")}
    totals = {("O", 7.5): ("TM0", "TO0"), ("U", 7.5): ("TM0", "TU0"),
              ("O", 8.5): ("TM1", "TO1"), ("U", 8.5): ("TM1", "TU1")}
    runners = {"fg": {"spreads": spreads, "totals": totals},
               "f5": {"spreads": {}, "totals": {}}}
    fetchers = {
        "fetch_fd_events": lambda session: [{"id": "e1"}],
        "fetch_event_runners": lambda session, eid, h, a: runners,
    }
    rows = fanduel.price_sgps(targets, periods=("FG",), client=MagicMock(),
                              parallelism=2, fetchers=fetchers)
    assert len(rows) == 2 * 4
```

- [ ] **Step 2: Run to verify failure**

Run: `/Users/callancapitolo/NFLWork/mlb_sgp/venv/bin/python -m pytest mlb_sgp/tests/test_dk_target_parallelism.py mlb_sgp/tests/test_fd_target_parallelism.py -q -k hook`
Expected: FAIL — unexpected keyword argument `fetchers`.

- [ ] **Step 3: Implement**

In `draftkings.py::price_sgps`: add `fetchers: dict | None = None` to
the signature (after `parallelism`). Immediately after the lazy
`from scraper_draftkings_sgp import (...)` block, add:

```python
    # Structure-fetch seam: SGPService injects TTL-cached wrappers here;
    # the dashboard shims pass nothing and get the legacy direct fetches.
    _f = fetchers or {}
    fetch_events_fn = _f.get("fetch_dk_events", fetch_dk_events)
    fetch_nums_fn = _f.get("fetch_main_market_nums", fetch_main_market_nums)
    fetch_selids_fn = _f.get("fetch_selection_ids", fetch_selection_ids)
```

Then replace the three call sites:
- `dk_events = fetch_dk_events(client.session)` → `dk_events = fetch_events_fn(client.session)`
- `main_nums = fetch_main_market_nums(client.session, game["dk_event_id"])` → `main_nums = fetch_nums_fn(client.session, game["dk_event_id"])`
- `sel_ids_all = fetch_selection_ids(client.session, game["dk_event_id"], main_nums, verbose,)` → `sel_ids_all = fetch_selids_fn(client.session, game["dk_event_id"], main_nums, verbose)`

In `fanduel.py::price_sgps`: add `fetchers: dict | None = None` (after
`parallelism`); after the lazy import block add:

```python
    _f = fetchers or {}
    fetch_events_fn = _f.get("fetch_fd_events", fetch_fd_events)
    fetch_runners_fn = _f.get("fetch_event_runners", fetch_event_runners)
```

Replace the two call sites:
- `fd_events = fetch_fd_events(client.session)` → `fd_events = fetch_events_fn(client.session)`
- `runners_cache[game_id] = fetch_event_runners(client.session, game["fd_event_id"], game["fd_home"], game["fd_away"],)` → `runners_cache[game_id] = fetch_runners_fn(client.session, game["fd_event_id"], game["fd_home"], game["fd_away"])`

- [ ] **Step 4: Run both parallelism test files + both orchestrator suites**

Run: `/Users/callancapitolo/NFLWork/mlb_sgp/venv/bin/python -m pytest mlb_sgp/tests/test_dk_target_parallelism.py mlb_sgp/tests/test_fd_target_parallelism.py mlb_sgp/tests/test_draftkings_orchestrator.py mlb_sgp/tests/test_fanduel_orchestrator.py -q`
Expected: ALL PASS.

- [ ] **Step 5: Commit**

```bash
git add mlb_sgp/draftkings.py mlb_sgp/fanduel.py \
  mlb_sgp/tests/test_dk_target_parallelism.py mlb_sgp/tests/test_fd_target_parallelism.py
git commit -m "feat(sgp): injectable structure-fetch hooks in DK/FD orchestrators

Default path unchanged (dashboard shims pass nothing). SGPService will
inject TTL-cached wrappers so event lists and DK's 2MB sel-id payload
stop being refetched every cycle."
```

---

### Task 7: `SGPService`

**Files:**
- Create: `kalshi_common/sgp_service.py`
- Modify: `kalshi_common/sgp_runner.py` (one re-export line)
- Test: `kalshi_mlb_rfq/tests/test_sgp_service.py` (create)

- [ ] **Step 1: Write the failing tests**

Create `kalshi_mlb_rfq/tests/test_sgp_service.py`:

```python
"""Tests for kalshi_common.sgp_service.SGPService.

All tests inject fake `runners` (the test seam) — no HTTP, no clients.
"""
import time
from datetime import datetime, timezone

from kalshi_common.sgp_service import SGPService
from mlb_sgp._shared import PricedRow, TargetLine

CT = datetime(2026, 6, 8, 23, 0, tzinfo=timezone.utc)
TARGETS = [TargetLine(game_id="g1", home_team="H", away_team="A",
                      commence_time=CT, period="FG", spread=-1.5, total=8.5)]


def _row(book):
    return PricedRow(game_id="g1", combo="Home Spread + Over", period="FG",
                     spread_line=-1.5, total_line=8.5, bookmaker=book,
                     source=f"{book}_direct", sgp_decimal=3.5,
                     sgp_american=250, fetch_time=CT)


def test_refresh_returns_rows_per_successful_book():
    svc = SGPService(books=("draftkings", "novig"), runners={
        "draftkings": lambda t: [_row("draftkings")],
        "novig": lambda t: [],
    })
    out = svc.refresh(TARGETS)
    assert set(out) == {"draftkings", "novig"}
    assert len(out["draftkings"]) == 1
    assert out["novig"] == []          # success with zero rows != failure


def test_failed_book_is_none_and_does_not_affect_others():
    def boom(t):
        raise RuntimeError("dead session")
    svc = SGPService(books=("draftkings", "fanduel"), runners={
        "draftkings": boom,
        "fanduel": lambda t: [_row("fanduel")],
    })
    out = svc.refresh(TARGETS)
    assert out["draftkings"] is None
    assert len(out["fanduel"]) == 1


def test_client_reset_after_three_consecutive_failures():
    def boom(t):
        raise RuntimeError("dead session")
    svc = SGPService(books=("draftkings",), runners={"draftkings": boom})
    sentinel = object()
    svc._state["draftkings"].client = sentinel
    svc.refresh(TARGETS)
    svc.refresh(TARGETS)
    assert svc._state["draftkings"].client is sentinel   # not yet
    svc.refresh(TARGETS)                                 # 3rd failure
    assert svc._state["draftkings"].client is None       # torn down


def test_success_resets_failure_counter():
    state = {"fail": True}

    def flaky(t):
        if state["fail"]:
            raise RuntimeError("x")
        return [_row("draftkings")]
    svc = SGPService(books=("draftkings",), runners={"draftkings": flaky})
    svc._state["draftkings"].client = object()
    svc.refresh(TARGETS); svc.refresh(TARGETS)       # 2 failures
    state["fail"] = False
    svc.refresh(TARGETS)                             # success
    state["fail"] = True
    svc.refresh(TARGETS)                             # 1 failure (not 3rd)
    assert svc._state["draftkings"].client is not None


def test_per_book_deadline_times_out_slow_book_only():
    def slow(t):
        time.sleep(2.0)
        return [_row("prophetx")]
    svc = SGPService(books=("prophetx", "novig"),
                     per_book_deadline_sec=0.3,
                     runners={"prophetx": slow,
                              "novig": lambda t: [_row("novig")]})
    t0 = time.monotonic()
    out = svc.refresh(TARGETS)
    assert out["prophetx"] is None
    assert len(out["novig"]) == 1
    assert time.monotonic() - t0 < 1.5   # did not wait the full 2s sleep


def test_min_refresh_skips_book_until_due():
    clock = {"t": 1000.0}
    calls = {"n": 0}

    def px_runner(t):
        calls["n"] += 1
        return [_row("prophetx")]
    svc = SGPService(books=("prophetx",),
                     min_refresh_sec={"prophetx": 120},
                     now_fn=lambda: clock["t"],
                     runners={"prophetx": px_runner})
    out1 = svc.refresh(TARGETS)
    assert "prophetx" in out1 and calls["n"] == 1
    clock["t"] += 30                      # 30s later: not due
    out2 = svc.refresh(TARGETS)
    assert "prophetx" not in out2 and calls["n"] == 1
    clock["t"] += 200                     # well past 120s: due again
    out3 = svc.refresh(TARGETS)
    assert "prophetx" in out3 and calls["n"] == 2
```

- [ ] **Step 2: Run to verify failure**

Run: `/Users/callancapitolo/NFLWork/kalshi_mlb_rfq/venv/bin/python -m pytest kalshi_mlb_rfq/tests/test_sgp_service.py -q`
Expected: FAIL — `ModuleNotFoundError: kalshi_common.sgp_service`.

- [ ] **Step 3: Write `kalshi_common/sgp_service.py`**

```python
"""In-process SGP pricing service shared by the Kalshi MLB bots.

Replaces the subprocess-per-cycle scraper model for the taker
(kalshi_mlb_rfq) and maker (kalshi_mlb_mm): one persistent HTTP client
per book held across cycles (no per-cycle TLS handshake), the four book
orchestrators run concurrently under a per-book deadline, and slow
structure fetches (event lists, DK's 2MB selection-id payload) are
TTL-cached. Prices are NEVER cached — every refresh() re-prices.

The dashboard keeps the CLI-shim subprocess path and never touches this.
"""
from __future__ import annotations
import logging
import sys
import time
from concurrent.futures import ThreadPoolExecutor
from pathlib import Path

from mlb_sgp._shared import PricedRow, TargetLine, TTLCache

log = logging.getLogger(__name__)

_REPO_ROOT = Path(__file__).resolve().parents[1]
_MLB_SGP_DIR = _REPO_ROOT / "mlb_sgp"

DEFAULT_BOOKS = ("draftkings", "fanduel", "prophetx", "novig")

# TTLs for structure fetches (market *structure*, never prices).
EVENTS_TTL_SEC = 300       # game list churns ~daily
STRUCTURE_TTL_SEC = 180    # sel-ids / runners churn when a book re-mains a line


class _BookState:
    __slots__ = ("client", "failures", "last_success", "caches")

    def __init__(self):
        self.client = None          # persistent per-book HTTP client
        self.failures = 0           # consecutive refresh failures
        self.last_success = None    # now_fn timestamp of last success
        self.caches = {}            # name -> TTLCache


class SGPService:
    """Holds persistent book clients; prices targets at all due books
    concurrently. See refresh() for the result contract."""

    MAX_FAILURES_BEFORE_REINIT = 3

    def __init__(
        self,
        books: tuple[str, ...] = DEFAULT_BOOKS,
        per_book_deadline_sec: float = 75.0,
        min_refresh_sec: dict[str, float] | None = None,
        runners: dict | None = None,   # test seam: {book: callable(targets)->rows}
        now_fn=time.monotonic,
    ):
        self.books = tuple(books)
        self.per_book_deadline_sec = per_book_deadline_sec
        self.min_refresh_sec = dict(min_refresh_sec or {})
        self._runners = runners
        self._now = now_fn
        self._state = {b: _BookState() for b in self.books}
        # The orchestrators lazily import the legacy scraper modules by
        # top-level name (`from scraper_draftkings_sgp import ...`),
        # which only resolves with mlb_sgp/ itself on sys.path — true
        # for CLI runs (cwd=mlb_sgp/) but not for the bots (cwd=repo
        # root). Make it resolvable here, once.
        if str(_MLB_SGP_DIR) not in sys.path:
            sys.path.insert(0, str(_MLB_SGP_DIR))

    # ------------------------------------------------------------------ #
    # Public API                                                          #
    # ------------------------------------------------------------------ #

    def refresh(self, targets: list[TargetLine]) -> dict[str, list[PricedRow] | None]:
        """Price `targets` at every due book concurrently.

        Returns {book: result} for books ATTEMPTED this call:
          list[PricedRow] (possibly empty) -> success
          None                             -> failure or deadline timeout
        Books skipped by min_refresh_sec are ABSENT from the dict —
        callers must not clear their previously-written rows.
        """
        due = [b for b in self.books if self._due(b)]
        results: dict[str, list[PricedRow] | None] = {}
        if not due:
            return results
        # Fresh pool per refresh: a hung book thread from a previous
        # cycle must not occupy a worker slot forever. shutdown(wait=
        # False) lets a still-running (timed-out) thread finish in the
        # background; its client gets torn down via the failure path.
        pool = ThreadPoolExecutor(max_workers=len(due),
                                  thread_name_prefix="sgp-book")
        try:
            futs = {b: pool.submit(self._run_book_safe, b, targets) for b in due}
            wall_deadline = time.monotonic() + self.per_book_deadline_sec
            for b, fut in futs.items():
                remaining = max(0.0, wall_deadline - time.monotonic())
                try:
                    rows = fut.result(timeout=remaining)
                except Exception:          # TimeoutError or runner crash
                    rows = None
                self._book_done(b, rows)
                results[b] = rows
        finally:
            pool.shutdown(wait=False, cancel_futures=True)
        return results

    def close(self):
        """Drop all persistent clients (sessions close on GC)."""
        for st in self._state.values():
            st.client = None

    # ------------------------------------------------------------------ #
    # Internals                                                           #
    # ------------------------------------------------------------------ #

    def _due(self, book: str) -> bool:
        min_s = self.min_refresh_sec.get(book, 0)
        st = self._state[book]
        if not min_s or st.last_success is None:
            return True
        return (self._now() - st.last_success) >= min_s

    def _book_done(self, book: str, rows) -> None:
        st = self._state[book]
        if rows is None:
            st.failures += 1
            log.warning("sgp_service: %s failed (consecutive=%d)",
                        book, st.failures)
            if st.failures >= self.MAX_FAILURES_BEFORE_REINIT:
                log.warning("sgp_service: %s client torn down for reinit", book)
                st.client = None
                st.failures = 0
        else:
            st.failures = 0
            st.last_success = self._now()

    def _run_book_safe(self, book: str, targets):
        """Runs in a worker thread. Returns rows or None — never raises
        (a raise would surface as a generic future error and lose the
        book attribution in logs)."""
        try:
            if self._runners is not None:
                return self._runners[book](targets)
            return self._run_book(book, targets)
        except Exception as e:
            log.warning("sgp_service: %s runner error: %s", book, e)
            return None

    def _run_book(self, book: str, targets):
        st = self._state[book]
        if book == "draftkings":
            from mlb_sgp import draftkings as mod
            if st.client is None:
                from mlb_sgp.dk_client import DraftKingsClient
                st.client = DraftKingsClient()
                st.caches = {"events": TTLCache(EVENTS_TTL_SEC),
                             "structure": TTLCache(STRUCTURE_TTL_SEC)}
            return mod.price_sgps(targets, periods=("FG",), client=st.client,
                                  fetchers=self._dk_fetchers(st))
        if book == "fanduel":
            from mlb_sgp import fanduel as mod
            if st.client is None:
                from mlb_sgp.fd_client import FanDuelClient
                st.client = FanDuelClient()
                st.caches = {"events": TTLCache(EVENTS_TTL_SEC),
                             "structure": TTLCache(STRUCTURE_TTL_SEC)}
            return mod.price_sgps(targets, periods=("FG",), client=st.client,
                                  fetchers=self._fd_fetchers(st))
        if book == "prophetx":
            from mlb_sgp import prophetx as mod
            if st.client is None:
                from mlb_sgp.prophetx_client import ProphetXClient
                st.client = ProphetXClient()
            return mod.price_sgps(targets, periods=("FG",), client=st.client)
        if book == "novig":
            from mlb_sgp import novig as mod
            if st.client is None:
                from mlb_sgp.novig_client import NovigClient
                st.client = NovigClient()
            return mod.price_sgps(targets, periods=("FG",), client=st.client)
        raise ValueError(f"unknown book {book!r}")

    @staticmethod
    def _dk_fetchers(st: _BookState) -> dict:
        import scraper_draftkings_sgp as legacy
        ev, struct = st.caches["events"], st.caches["structure"]
        return {
            "fetch_dk_events": lambda session: ev.get_or_fetch(
                "dk_events", lambda: legacy.fetch_dk_events(session)),
            "fetch_main_market_nums": lambda session, eid: struct.get_or_fetch(
                ("nums", eid), lambda: legacy.fetch_main_market_nums(session, eid)),
            "fetch_selection_ids": lambda session, eid, nums, verbose=False:
                struct.get_or_fetch(
                    ("selids", eid),
                    lambda: legacy.fetch_selection_ids(session, eid, nums, verbose)),
        }

    @staticmethod
    def _fd_fetchers(st: _BookState) -> dict:
        import scraper_fanduel_sgp as legacy
        ev, struct = st.caches["events"], st.caches["structure"]
        return {
            "fetch_fd_events": lambda session: ev.get_or_fetch(
                "fd_events", lambda: legacy.fetch_fd_events(session)),
            "fetch_event_runners": lambda session, eid, h, a: struct.get_or_fetch(
                ("runners", eid),
                lambda: legacy.fetch_event_runners(session, eid, h, a)),
        }
```

3b. In `kalshi_common/sgp_runner.py`, add below the existing imports:

```python
from kalshi_common.sgp_service import SGPService  # noqa: F401  (re-export)
```

- [ ] **Step 4: Run the tests**

Run: `/Users/callancapitolo/NFLWork/kalshi_mlb_rfq/venv/bin/python -m pytest kalshi_mlb_rfq/tests/test_sgp_service.py -q`
Expected: ALL PASS.

- [ ] **Step 5: Commit**

```bash
git add kalshi_common/sgp_service.py kalshi_common/sgp_runner.py \
  kalshi_mlb_rfq/tests/test_sgp_service.py
git commit -m "feat(sgp): in-process SGPService — persistent clients, per-book deadline, TTL caches

Holds one client per book across cycles (kills the per-cycle TLS
handshake + 2MB DK sel-id refetch), runs the 4 orchestrators
concurrently, isolates failures per book (3 strikes -> client reinit),
and supports per-book min_refresh_sec (the PX footprint dial, off by
default)."
```

---

### Task 8: Service-backed `sgp_cycle`

**Files:**
- Modify: `kalshi_common/sgp_runner.py` (the `sgp_cycle` function, line ~367)
- Test: append to `kalshi_mlb_rfq/tests/test_sgp_runner.py`

`sgp_cycle` gains `service=None`. With a service: in-process refresh,
write rows to the bot market DB, return per-book row counts (−1 =
failed book, old rows preserved). Without: the legacy subprocess path,
byte-identical (rollback hatch + anything else calling it).

- [ ] **Step 1: Write the failing tests** (append to `kalshi_mlb_rfq/tests/test_sgp_runner.py`)

```python
def test_sgp_cycle_service_path_writes_rows_and_preserves_failed_books(
        tmp_path, monkeypatch):
    """Service path: success -> stale source rows replaced; failure ->
    old rows preserved (mirrors today's subprocess-crash behavior)."""
    import duckdb
    from datetime import datetime, timezone
    from kalshi_common import sgp_runner
    from mlb_sgp import db as sgp_db
    from mlb_sgp._shared import PricedRow, TargetLine

    db_path = str(tmp_path / "market.duckdb")
    ct = datetime(2026, 6, 8, 23, 0, tzinfo=timezone.utc)
    target = TargetLine(game_id="g1", home_team="H", away_team="A",
                        commence_time=ct, period="FG", spread=-1.5, total=8.5)
    monkeypatch.setattr(sgp_runner, "enumerate_kalshi_targets", lambda: [target])

    # Pre-seed: one stale DK row (must be replaced) + one FD row (must
    # survive, because FD "fails" this cycle).
    stale_dk = PricedRow(game_id="gOLD", combo="Home Spread + Over",
                         period="FG", spread_line=-1.5, total_line=9.5,
                         bookmaker="draftkings", source="draftkings_direct",
                         sgp_decimal=9.99, sgp_american=899, fetch_time=ct)
    old_fd = PricedRow(game_id="gOLD", combo="Home Spread + Over",
                       period="FG", spread_line=-1.5, total_line=9.5,
                       bookmaker="fanduel", source="fanduel_direct",
                       sgp_decimal=8.88, sgp_american=788, fetch_time=ct)
    sgp_db.upsert_priced_rows([stale_dk, old_fd], db_path=db_path)

    fresh_dk = PricedRow(game_id="g1", combo="Home Spread + Over",
                         period="FG", spread_line=-1.5, total_line=8.5,
                         bookmaker="draftkings", source="draftkings_direct",
                         sgp_decimal=3.5, sgp_american=250, fetch_time=ct)

    class FakeService:
        def refresh(self, targets):
            assert targets == [target]
            return {"draftkings": [fresh_dk], "fanduel": None}

    counts = sgp_runner.sgp_cycle(bot_market_db=db_path, service=FakeService())
    assert counts == {"draftkings": 1, "fanduel": -1}

    con = duckdb.connect(db_path, read_only=True)
    try:
        rows = con.execute(
            "SELECT bookmaker, game_id, sgp_decimal FROM mlb_sgp_odds "
            "ORDER BY bookmaker").fetchall()
    finally:
        con.close()
    assert rows == [("draftkings", "g1", 3.5), ("fanduel", "gOLD", 8.88)]

    # Target lines were still written (tipoff gating reads them).
    con = duckdb.connect(db_path, read_only=True)
    try:
        n = con.execute("SELECT COUNT(*) FROM mlb_target_lines").fetchone()[0]
    finally:
        con.close()
    assert n == 1


def test_sgp_cycle_service_path_skipped_book_rows_untouched(tmp_path, monkeypatch):
    """A book absent from refresh() results (min_refresh skip) keeps its
    rows AND is absent from the returned counts."""
    from datetime import datetime, timezone
    from kalshi_common import sgp_runner
    from mlb_sgp import db as sgp_db
    from mlb_sgp._shared import PricedRow, TargetLine
    import duckdb

    db_path = str(tmp_path / "market.duckdb")
    ct = datetime(2026, 6, 8, 23, 0, tzinfo=timezone.utc)
    monkeypatch.setattr(sgp_runner, "enumerate_kalshi_targets", lambda: [
        TargetLine(game_id="g1", home_team="H", away_team="A",
                   commence_time=ct, period="FG", spread=-1.5, total=8.5)])
    px_row = PricedRow(game_id="g1", combo="Home Spread + Over", period="FG",
                       spread_line=-1.5, total_line=8.5, bookmaker="prophetx",
                       source="prophetx_direct", sgp_decimal=3.3,
                       sgp_american=230, fetch_time=ct)
    sgp_db.upsert_priced_rows([px_row], db_path=db_path)

    class FakeService:
        def refresh(self, targets):
            return {}   # everything skipped

    counts = sgp_runner.sgp_cycle(bot_market_db=db_path, service=FakeService())
    assert counts == {}
    con = duckdb.connect(db_path, read_only=True)
    try:
        n = con.execute("SELECT COUNT(*) FROM mlb_sgp_odds").fetchone()[0]
    finally:
        con.close()
    assert n == 1
```

- [ ] **Step 2: Run to verify failure**

Run: `/Users/callancapitolo/NFLWork/kalshi_mlb_rfq/venv/bin/python -m pytest kalshi_mlb_rfq/tests/test_sgp_runner.py -q -k service_path`
Expected: FAIL — `sgp_cycle() got an unexpected keyword argument 'service'`.

- [ ] **Step 3: Implement in `kalshi_common/sgp_runner.py`**

Replace the existing `sgp_cycle` function with:

```python
# Book name -> orchestrator module (source-label owner). Used by the
# service path to clear exactly the labels each book writes.
_BOOK_MODULES = {
    "draftkings": "mlb_sgp.draftkings",
    "fanduel": "mlb_sgp.fanduel",
    "prophetx": "mlb_sgp.prophetx",
    "novig": "mlb_sgp.novig",
}


def sgp_cycle(
    bot_market_db: str,
    scraper_dir: str | None = None,
    venv_python: str | None = None,
    timeout_sec: int | None = None,
    service=None,
) -> dict[str, int]:
    """One full SGP scrape tick.

    With `service` (an SGPService): in-process path —
      1. Enumerate Kalshi MVE -> list[TargetLine]
      2. Write mlb_target_lines (tipoff gating + debugging read it)
      3. service.refresh(targets); per successful book, clear that
         book's source labels and upsert its fresh rows. Failed books
         (None) keep their old rows — same outcome as a crashed
         subprocess today. Returns {book: row_count, failed: -1}.

    Without `service`: legacy subprocess path, unchanged. Returns
    {scraper_name: return_code}. Kept as the rollback hatch
    (scraper_dir / venv_python / timeout_sec are required then).
    """
    targets = enumerate_kalshi_targets()
    write_target_lines(targets, db_path=bot_market_db)

    if service is None:
        if not (scraper_dir and venv_python and timeout_sec):
            raise ValueError(
                "sgp_cycle: legacy subprocess path needs scraper_dir, "
                "venv_python and timeout_sec")
        return run_scrapers(
            scraper_dir=scraper_dir,
            scraper_names=SCRAPER_NAMES,
            venv_python=venv_python,
            timeout_sec=timeout_sec,
            env={
                "MLB_SGP_DB_PATH": bot_market_db,
                "MLB_SGP_PERIODS": "FG",
            },
        )

    import importlib
    from mlb_sgp import db as sgp_db

    results = service.refresh(targets)
    counts: dict[str, int] = {}
    for book, rows in results.items():
        if rows is None:
            counts[book] = -1
            continue
        mod = importlib.import_module(_BOOK_MODULES[book])
        sgp_db.clear_source(mod.SOURCE_LABEL, db_path=bot_market_db)
        fallback_label = getattr(mod, "SOURCE_LABEL_FALLBACK", None)
        if fallback_label:
            sgp_db.clear_source(fallback_label, db_path=bot_market_db)
        sgp_db.upsert_priced_rows(rows, db_path=bot_market_db)
        counts[book] = len(rows)
    return counts
```

- [ ] **Step 4: Run the full sgp_runner suite**

Run: `/Users/callancapitolo/NFLWork/kalshi_mlb_rfq/venv/bin/python -m pytest kalshi_mlb_rfq/tests/test_sgp_runner.py -q`
Expected: ALL PASS (legacy-path tests must still pass untouched).

- [ ] **Step 5: Commit**

```bash
git add kalshi_common/sgp_runner.py kalshi_mlb_rfq/tests/test_sgp_runner.py
git commit -m "feat(sgp): service-backed sgp_cycle (in-process path; subprocess path kept for rollback)"
```

---

### Task 9: Bot venv dependencies

**Files:**
- Modify: `kalshi_mlb_rfq/requirements.txt`, `kalshi_mlb_mm/requirements.txt`

The bots' venvs lack `curl_cffi` (verified 2026-06-08) because scrapers
used to run under `mlb_sgp/venv`. In-process pricing imports it inside
the bot process.

- [ ] **Step 1: Add the dependency to both requirements files**

Append to `kalshi_mlb_rfq/requirements.txt` and `kalshi_mlb_mm/requirements.txt`:

```
curl_cffi>=0.6      # in-process SGP scraping (SGPService) — matches mlb_sgp/venv
```

- [ ] **Step 2: Install into both live venvs** (additive — does not affect the running bots until restart)

```bash
/Users/callancapitolo/NFLWork/kalshi_mlb_rfq/venv/bin/pip install "curl_cffi>=0.6"
/Users/callancapitolo/NFLWork/kalshi_mlb_mm/venv/bin/pip install "curl_cffi>=0.6"
```

- [ ] **Step 3: Verify**

```bash
/Users/callancapitolo/NFLWork/kalshi_mlb_rfq/venv/bin/python -c "import curl_cffi; print('taker OK')"
/Users/callancapitolo/NFLWork/kalshi_mlb_mm/venv/bin/python -c "import curl_cffi; print('maker OK')"
```
Expected: `taker OK` / `maker OK`.

- [ ] **Step 4: Commit**

```bash
git add kalshi_mlb_rfq/requirements.txt kalshi_mlb_mm/requirements.txt
git commit -m "chore(bots): add curl_cffi for in-process SGP scraping"
```

---

### Task 10: Taker integration (`kalshi_mlb_rfq`)

**Files:**
- Modify: `kalshi_mlb_rfq/main.py` (two `sgp_cycle` call sites: warm-up ~line 1524, cadence ~line 1581)

No new unit tests — the service and cycle paths are tested in Tasks 7–8;
this is wiring. Verified by compile check here and live in Task 13.

- [ ] **Step 1: Construct the service at startup**

In `main_loop`, just before the warm-up block (the
`log.info("startup: warming SGP cache ...")` line), add:

```python
    from kalshi_common.sgp_service import SGPService
    sgp_service = SGPService(per_book_deadline_sec=config.SGP_SCRAPER_TIMEOUT_SEC)
```

- [ ] **Step 2: Switch both call sites to the service path**

Warm-up call (~line 1524) — replace:

```python
        rcs = sgp_runner.sgp_cycle(
            bot_market_db=str(config.BOT_MARKET_DB),
            scraper_dir=str(config.MLB_SGP_DIR),
            venv_python=str(config.MLB_SGP_DIR / "venv" / "bin" / "python"),
            timeout_sec=config.SGP_SCRAPER_TIMEOUT_SEC,
        )
```

with:

```python
        rcs = sgp_runner.sgp_cycle(
            bot_market_db=str(config.BOT_MARKET_DB),
            service=sgp_service,
        )
```

Cadence call (~line 1581): same replacement (the surrounding
`try/except` + `log.info("sgp_cycle: rcs=%s ...")` lines stay — `rcs`
is now per-book row counts instead of subprocess return codes, and the
`%s` format handles either).

- [ ] **Step 3: Compile check**

Run: `/Users/callancapitolo/NFLWork/kalshi_mlb_rfq/venv/bin/python -m py_compile kalshi_mlb_rfq/main.py && echo OK`
Expected: `OK`.

- [ ] **Step 4: Run the taker's full test suite**

Run: `/Users/callancapitolo/NFLWork/kalshi_mlb_rfq/venv/bin/python -m pytest kalshi_mlb_rfq/tests/ -q`
Expected: ALL PASS (one pre-existing broken env test is known — see
memory `kalshi_mlb_rfq_research_logging`; a failure that exists on
`main` too is not a regression from this work. Verify by running the
same test on `main` if anything fails.)

- [ ] **Step 5: Commit**

```bash
git add kalshi_mlb_rfq/main.py
git commit -m "feat(taker): SGP cycle via in-process SGPService (no more subprocess spawn)"
```

---

### Task 11: Maker integration (`kalshi_mlb_mm`)

**Files:**
- Modify: `kalshi_mlb_mm/main.py` (warm-up call ~line 822, cadence call ~line 867)

- [ ] **Step 1: Construct the service at startup**

In `main_loop`, just before the
`# synchronous warm-up: one SGP cycle` comment, add:

```python
    from kalshi_common.sgp_service import SGPService
    sgp_service = SGPService(per_book_deadline_sec=config.SGP_SCRAPER_TIMEOUT_SEC)
```

- [ ] **Step 2: Switch both call sites**

Replace both occurrences of:

```python
        rc = sgp_runner.sgp_cycle(bot_market_db=str(config.MARKET_DB),
                                   scraper_dir=str(config.MLB_SGP_DIR),
                                   venv_python=str(config.MLB_SGP_DIR / "venv" / "bin" / "python"),
                                   timeout_sec=config.SGP_SCRAPER_TIMEOUT_SEC)
```

with:

```python
        rc = sgp_runner.sgp_cycle(bot_market_db=str(config.MARKET_DB),
                                   service=sgp_service)
```

(Use `grep -n "sgp_runner.sgp_cycle" kalshi_mlb_mm/main.py` to find both;
exact indentation may differ between the two sites — preserve it.)

- [ ] **Step 3: Compile check + maker test suite**

Run: `/Users/callancapitolo/NFLWork/kalshi_mlb_mm/venv/bin/python -m py_compile kalshi_mlb_mm/main.py && /Users/callancapitolo/NFLWork/kalshi_mlb_mm/venv/bin/python -m pytest kalshi_mlb_mm/tests/ -q`
Expected: compile OK, tests PASS.

- [ ] **Step 4: Commit**

```bash
git add kalshi_mlb_mm/main.py
git commit -m "feat(maker): SGP cycle via in-process SGPService"
```

---

### Task 12: Documentation

**Files:**
- Modify: `mlb_sgp/README.md`, `kalshi_mlb_rfq/README.md`, `kalshi_mlb_mm/README.md`, `CLAUDE.md` (root), `docs/superpowers/specs/2026-06-08-sgp-scraper-speedup-design.md`

- [ ] **Step 1: `mlb_sgp/README.md`** — in the "Concurrency & Logs" section, add a "Target-level parallelism" subsection:

```markdown
### Target-level parallelism (2026-06)

Each orchestrator fans target lines out on a thread pool; the 4-combo
pool nests inside (total in-flight requests = parallelism × 4):

| Book | Default | Env override | Probed ceiling |
|---|---|---|---|
| DraftKings | 8 | `MLB_SGP_DK_PARALLELISM` | (fill from probe) |
| FanDuel | 4 | `MLB_SGP_FD_PARALLELISM` | not probed (not the long pole) |
| ProphetX | 2 | `MLB_SGP_PX_PARALLELISM` | (fill from probe) |
| Novig | 4 | n/a (module constant) | not probed |

`price_sgps(..., parallelism=N)` overrides both. Ceilings come from
`probe_concurrency.py` — a **manual-only** ramp harness (budgeted:
≤150 DK calls, ≤40 PX RFQs, cooldowns, post-run health check). Run it
in the morning (~7–8am PT, no games live):

    mlb_sgp/venv/bin/python -m mlb_sgp.probe_concurrency --book dk

The Kalshi bots no longer spawn these scrapers as subprocesses — they
price in-process via `kalshi_common/sgp_service.py::SGPService`
(persistent sessions + TTL-cached structure fetches). The dashboard
still uses the CLI shims, which behave exactly as before.
```

- [ ] **Step 2: `kalshi_mlb_rfq/README.md` and `kalshi_mlb_mm/README.md`** — find the SGP-cadence wording (`grep -n "sgp_cycle\|subprocess\|scraper" <README>`) and update: the bot prices SGPs in-process via `SGPService` (persistent per-book HTTP clients, concurrent books, per-book deadline = `SGP_SCRAPER_TIMEOUT_SEC`, TTL-cached structure fetches, failed book keeps prior rows). Mention the legacy subprocess path remains available by calling `sgp_cycle` without `service` (rollback).

- [ ] **Step 3: Root `CLAUDE.md`** — in the `kalshi_common/` bullet, extend the `sgp_runner` mention:

```
`sgp_runner` (SGP scrape orchestration + in-process `SGPService` — persistent book clients, both bots price in-process; dashboard still uses CLI shims)
```

- [ ] **Step 4: Spec touch-up** — in the spec's §3.5 table, mark the PX/NV row: structure caching deferred (their per-cycle structure fetches are 1–2 cheap calls; persistent clients cover them) and note `SGPService` lives in `kalshi_common/sgp_service.py` re-exported via `sgp_runner` (§3.4).

- [ ] **Step 5: Commit**

```bash
git add mlb_sgp/README.md kalshi_mlb_rfq/README.md kalshi_mlb_mm/README.md \
  CLAUDE.md docs/superpowers/specs/2026-06-08-sgp-scraper-speedup-design.md
git commit -m "docs(sgp): concurrency section, SGPService architecture, spec deviations"
```

---

### Task 13: Live verification & tuning (pre-merge gate)

**Files:**
- Modify: `mlb_sgp/draftkings.py` / `mlb_sgp/prophetx.py` (final defaults), `mlb_sgp/README.md` (probed ceilings)

No code until the measurements are in. **Steps 1–2 must run in the
morning window (~7–8am PT, no MLB games in progress)** per the spec.
Everything runs from the worktree against **copies** of live DBs —
never the live files.

- [ ] **Step 1: Run the DK probe**

```bash
cd /Users/callancapitolo/NFLWork/.claude/worktrees/sgp-scraper-speedup
cp /Users/callancapitolo/NFLWork/kalshi_mlb_rfq/kalshi_mlb_rfq_market.duckdb /tmp/probe_market.duckdb
/Users/callancapitolo/NFLWork/mlb_sgp/venv/bin/python -m mlb_sgp.probe_concurrency --book dk --db /tmp/probe_market.duckdb
```

Record the per-level table and the recommended parallelism. Confirm the
health check prints OK. (If `mlb_target_lines` is empty — bot idle /
off-season — note it and fall back to shipped defaults; do not invent
numbers.)

- [ ] **Step 2: Run the PX probe**

```bash
/Users/callancapitolo/NFLWork/mlb_sgp/venv/bin/python -m mlb_sgp.probe_concurrency --book px --db /tmp/probe_market.duckdb
```

Record the same. PX stops at ANY anomaly by design.

- [ ] **Step 3: Set final defaults from the probe**

Edit `DK_TARGET_PARALLELISM_DEFAULT` in `draftkings.py` and
`PX_TARGET_PARALLELISM_DEFAULT` in `prophetx.py` to the last clean
probe level (back off one level from any degradation point). Fill the
"Probed ceiling" column in `mlb_sgp/README.md`. Re-run both
parallelism test files (the env-default assertions in
`test_parallelism_default_reads_env` / `test_px_parallelism_resolution`
reference the shipped defaults — update those two assertions to the new
values in the same change).

- [ ] **Step 4: Supervised in-process refresh vs CLI-shim cross-check**

```bash
cp /Users/callancapitolo/NFLWork/kalshi_mlb_rfq/kalshi_mlb_rfq_market.duckdb /tmp/svc_market.duckdb
/Users/callancapitolo/NFLWork/kalshi_mlb_rfq/venv/bin/python - <<'EOF'
import time
from kalshi_common.sgp_service import SGPService
from mlb_sgp._shared import load_target_lines

targets = [t for t in load_target_lines("/tmp/svc_market.duckdb") if t.period == "FG"]
print(f"{len(targets)} targets")
svc = SGPService()
t0 = time.monotonic()
results = svc.refresh(targets)
elapsed = time.monotonic() - t0
for book, rows in results.items():
    print(f"{book}: {'FAILED' if rows is None else len(rows)} rows")
print(f"refresh wall-clock: {elapsed:.1f}s")
EOF
```

Acceptance: wall-clock **≤60s**; every book returns rows (not FAILED);
DK row count within ~10% of the recent subprocess runs (~1,400 on a
full slate — compare against `mlb_sgp/scraper_draftkings_sgp.py.runner.log`).
Spot-check ~20 DK `sgp_decimal` values against a concurrent CLI-shim
run on the same copied DB:
`MLB_SGP_DB_PATH=/tmp/svc_market.duckdb MLB_SGP_PERIODS=FG /Users/callancapitolo/NFLWork/mlb_sgp/venv/bin/python mlb_sgp/scraper_draftkings_sgp.py`
— differences should be live-line movement only (small, sparse), not
systematic.

- [ ] **Step 5: Golden regression + full unit suites**

```bash
/Users/callancapitolo/NFLWork/mlb_sgp/venv/bin/python -m pytest mlb_sgp/tests/ -q
/Users/callancapitolo/NFLWork/kalshi_mlb_rfq/venv/bin/python -m pytest kalshi_mlb_rfq/tests/ -q
/Users/callancapitolo/NFLWork/kalshi_mlb_mm/venv/bin/python -m pytest kalshi_mlb_mm/tests/ -q
```
Expected: ALL PASS (modulo the known pre-existing env-test failure noted in Task 10 — verify any failure also fails on `main` before attributing it).

- [ ] **Step 6: Dashboard parity**

Run the DK CLI shim against a **copy** of the dashboard DB and confirm
it still populates (proves the no-`fetchers`, no-`service` path is
untouched):

```bash
cp "/Users/callancapitolo/NFLWork/Answer Keys/mlb_mm.duckdb" /tmp/parity_mm.duckdb
MLB_SGP_DB_PATH=/tmp/parity_mm.duckdb /Users/callancapitolo/NFLWork/mlb_sgp/venv/bin/python mlb_sgp/scraper_draftkings_sgp.py
/Users/callancapitolo/NFLWork/mlb_sgp/venv/bin/python -c "
import duckdb
con = duckdb.connect('/tmp/parity_mm.duckdb', read_only=True)
print(con.execute(\"SELECT source, COUNT(*) FROM mlb_sgp_odds WHERE bookmaker='draftkings' GROUP BY source\").fetchall())
con.close()"
rm /tmp/parity_mm.duckdb /tmp/svc_market.duckdb /tmp/probe_market.duckdb
```

Expected: non-zero `draftkings_direct` count, no errors.

- [ ] **Step 7: Commit the tuning**

```bash
git add mlb_sgp/draftkings.py mlb_sgp/prophetx.py mlb_sgp/README.md \
  mlb_sgp/tests/test_dk_target_parallelism.py mlb_sgp/tests/test_prophetx_orchestrator.py
git commit -m "perf(sgp): set probed parallelism ceilings (DK=<N>, PX=<N>)

Probe results: <paste the per-level tables here>"
```

---

## Post-merge operational steps (after user approves the merge)

These happen on `main`, **not** in the worktree (per restart-gotchas
memory: bots must be restarted from the main repo cwd or `python -m`
loads the worktree package with no .env and wrong DBs):

1. Pre-merge executive engineer review of `git diff main..HEAD`
   (data integrity, resource safety, edge cases, dead code, log
   hygiene, secrets) → present ISSUES vs ACCEPTABLE RISKS → **explicit
   user approval to merge**.
2. Merge to `main`; `git worktree remove .claude/worktrees/sgp-scraper-speedup`
   and `git branch -d worktree-sgp-scraper-speedup`.
3. Restart both bots from `/Users/callancapitolo/NFLWork` (taker may
   need `kill -9` — SIGTERM is starved by the 640s rfq_refresh block;
   flag before killing per the kill-permission memory; the maker is a
   live-money MM — flag its restart too).
4. **1-hour soak:** watch both `bot.log`s — `sgp_cycle` wall-clock ≤60s,
   per-book counts comparable to pre-change (`rcs={'draftkings': ~1400,
   ...}` instead of return codes), no rising failure (-1) streaks.
5. Rollback hatch if a book misbehaves: revert the bot's call site to
   the legacy signature (subprocess path is untouched), or set the
   book's `MLB_SGP_*_PARALLELISM` env down and restart.
```
