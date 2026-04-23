# Asymmetric Outlier Flags Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Make the Cross-Book Grid / +EV Candidates flag non-Kalshi venues only when the book's YES implied probability is cheaper than the consensus fair (a bettable YES wager). Kalshi stays symmetric because both sides are purchasable.

**Architecture:** One-branch change inside `cross_book_grid()` — Kalshi keeps `abs(take - median) >= threshold`, every other book uses `(median - take) >= threshold`. `ev_candidates()` delegates to it so the filter propagates for free. UI copy, docstring, tests, and README are updated in lockstep.

**Tech Stack:** Python 3.11, DuckDB, Dash, pytest.

**Spec:** `docs/superpowers/specs/2026-04-23-asymmetric-outlier-flags-design.md`

---

## File Structure

Files created/modified and their responsibilities:

- **Modify** `nfl_draft/lib/queries.py` — change the flag condition inside `cross_book_grid()` and update its docstring. This is the only behavior change.
- **Modify** `nfl_draft/tests/unit/test_cross_book_grid.py` — repair the one test that relied on symmetric flagging; add two new tests covering both asymmetric cases (non-Kalshi above median not flagged, non-Kalshi below median flagged) and a Kalshi regression guard.
- **Modify** `nfl_draft/tests/integration/test_dashboard_queries.py` — `test_ev_candidates_ranks_by_abs_delta` currently expects an above-median non-Kalshi row to be an +EV candidate; that row will no longer flag. Adjust the fixture or assertions so both markets still produce bettable candidates.
- **Modify** `kalshi_draft/app.py` — replace the Cross-Book Grid blurb (lines 834-842) so users understand the asymmetric rule.
- **Modify** `nfl_draft/README.md` — replace the outlier-flag paragraph (lines 33-38) with the new rule and its rationale.

No new files.

---

### Task 1: Create worktree and feature branch

**Files:**
- None (infra only)

- [ ] **Step 1: Verify you are on `main` and the tree is clean**

Run: `git -C /Users/callancapitolo/NFLWork status --short && git -C /Users/callancapitolo/NFLWork branch --show-current`

Expected: empty output for `status --short` (or only the pre-existing untracked entries `.playwright-mcp/` and `docs/superpowers/plans/2026-04-22-bookmaker-no-popup.md`), and `main` for the branch.

- [ ] **Step 2: Create the worktree**

Run:
```bash
git -C /Users/callancapitolo/NFLWork worktree add \
  /Users/callancapitolo/NFLWork/.worktrees/asymmetric-flags \
  -b feature/asymmetric-outlier-flags
```

Expected: `Preparing worktree ...` and `HEAD is now at ...` messages.

- [ ] **Step 3: Verify the worktree**

Run: `git -C /Users/callancapitolo/NFLWork/.worktrees/asymmetric-flags branch --show-current`

Expected: `feature/asymmetric-outlier-flags`.

**From Task 2 onward, all file paths are rooted at `/Users/callancapitolo/NFLWork/.worktrees/asymmetric-flags/` — run every command from inside that worktree.**

---

### Task 2: Write the new failing tests (red)

**Files:**
- Modify: `nfl_draft/tests/unit/test_cross_book_grid.py`

Rationale: TDD — we add the asymmetric-rule tests first and prove they fail against the current symmetric implementation. We leave the existing symmetric-behavior test alone at this step; Task 4 will repair it once the implementation flips over.

- [ ] **Step 1: Append three new tests to the bottom of `nfl_draft/tests/unit/test_cross_book_grid.py`**

Append exactly this block after the last existing test:

```python
def test_non_kalshi_below_median_is_flagged(fresh_db):
    """Non-Kalshi venue whose YES implied_prob is well below the consensus
    median (bettable YES edge) must flag under the asymmetric rule.

    Median fair ~0.50; draftkings implied 0.30 -> median - take = 20pp > 10pp.
    """
    _insert_odds("m1", "kalshi",     100, 0.500, 0.500, datetime.now())
    _insert_odds("m1", "fanduel",    100, 0.500, 0.500, datetime.now())
    _insert_odds("m1", "draftkings", +233, 0.300, 0.305, datetime.now())
    rows = q.cross_book_grid(threshold_pp=10)
    m = rows[0]
    assert m["flags"]["draftkings"] is True


def test_non_kalshi_above_median_is_not_flagged(fresh_db):
    """Non-Kalshi venue whose YES implied_prob is well above the consensus
    median (book is overpricing YES — unactionable on YES-only sportsbook
    futures) must NOT flag under the asymmetric rule, even when the absolute
    delta exceeds the threshold.

    Median fair ~0.42; bookmaker implied 0.60 -> take - median = 18pp; because
    the book is above median and non-Kalshi, the flag stays False.
    """
    _insert_odds("m1", "kalshi",     100, 0.500, 0.420, datetime.now())
    _insert_odds("m1", "fanduel",    100, 0.500, 0.420, datetime.now())
    _insert_odds("m1", "draftkings", 100, 0.500, 0.420, datetime.now())
    _insert_odds("m1", "bookmaker",  -150, 0.600, 0.605, datetime.now())
    rows = q.cross_book_grid(threshold_pp=10)
    m = rows[0]
    assert m["flags"]["bookmaker"] is False


def test_kalshi_above_median_still_flags(fresh_db):
    """Regression guard: Kalshi must keep the symmetric rule. A Kalshi
    implied_prob well above the median (NO is underpriced -> bettable NO
    wager) must still flag.

    Median fair ~0.30; kalshi implied 0.55 -> take - median = 25pp; since
    Kalshi is symmetric, this flags.
    """
    _insert_odds("m1", "draftkings", +233, 0.300, 0.300, datetime.now())
    _insert_odds("m1", "fanduel",    +233, 0.300, 0.300, datetime.now())
    _insert_odds("m1", "kalshi",     -122, 0.550, 0.555, datetime.now())
    rows = q.cross_book_grid(threshold_pp=10)
    m = rows[0]
    assert m["flags"]["kalshi"] is True
```

- [ ] **Step 2: Run the new tests and confirm they fail as expected**

Run:
```bash
cd /Users/callancapitolo/NFLWork/.worktrees/asymmetric-flags && \
/Users/callancapitolo/NFLWork/kalshi_draft/venv/bin/pytest \
  nfl_draft/tests/unit/test_cross_book_grid.py::test_non_kalshi_below_median_is_flagged \
  nfl_draft/tests/unit/test_cross_book_grid.py::test_non_kalshi_above_median_is_not_flagged \
  nfl_draft/tests/unit/test_cross_book_grid.py::test_kalshi_above_median_still_flags -v
```

Expected output:
- `test_non_kalshi_below_median_is_flagged` PASSES (current symmetric rule also flags this case — `|0.30 - 0.50| = 20pp` triggers).
- `test_non_kalshi_above_median_is_not_flagged` FAILS with `assert True is False` — current rule flags the above-median case (bug the plan fixes).
- `test_kalshi_above_median_still_flags` PASSES (current symmetric rule already flags this case).

The one intentional failure (`test_non_kalshi_above_median_is_not_flagged`) is the red-phase target.

- [ ] **Step 3: Commit the red-phase tests**

```bash
cd /Users/callancapitolo/NFLWork/.worktrees/asymmetric-flags && \
git add nfl_draft/tests/unit/test_cross_book_grid.py && \
git commit -m "$(cat <<'EOF'
test(nfl_draft): add failing test for asymmetric outlier-flag rule

Adds three coverage tests for the asymmetric rule: non-Kalshi below-median
flags (bettable YES), non-Kalshi above-median does not flag (unactionable),
Kalshi above-median still flags (symmetric). The middle test is
intentionally red against the current symmetric implementation.

Co-Authored-By: Claude Opus 4.7 (1M context) <noreply@anthropic.com>
EOF
)"
```

Expected: one new commit on `feature/asymmetric-outlier-flags`.

---

### Task 3: Implement the asymmetric flag rule (green)

**Files:**
- Modify: `nfl_draft/lib/queries.py` (around lines 89-168)

- [ ] **Step 1: Replace the flag-assignment block in `cross_book_grid`**

In `nfl_draft/lib/queries.py`, find the loop:

```python
            flags: Dict[str, bool] = {}
            for book, r in books.items():
                # Uniform rule: flag if |raw take price - median fair| >= threshold.
                take = r["implied_prob"]
                flags[book] = (take is not None and abs(take - median) >= threshold)
```

Replace it with:

```python
            flags: Dict[str, bool] = {}
            for book, r in books.items():
                take = r["implied_prob"]
                if take is None:
                    flags[book] = False
                    continue
                # Asymmetric rule:
                #   * Kalshi is a two-sided binary (buy YES at yes_ask, buy NO at
                #     100 - yes_bid), so an edge in either direction is bettable.
                #     Keep the absolute-delta rule.
                #   * Every other book in this portal is a YES-only sportsbook
                #     futures market. If the book is pricing YES *above* the
                #     consensus fair, there's no way to take the other side;
                #     surfacing that flag is noise. Only flag when YES is cheap
                #     on this book (median - take >= threshold) -> bettable YES.
                if book == "kalshi":
                    flags[book] = abs(take - median) >= threshold
                else:
                    flags[book] = (median - take) >= threshold
```

- [ ] **Step 2: Rewrite the `cross_book_grid` docstring so the rule is explicit**

In the same file, replace the existing docstring of `cross_book_grid` (roughly lines 90-112) with:

```python
    """Return one row per market with all-venue prices, median, and outlier flags.

    Median is ``statistics.median(devig_prob across posting venues)`` — the
    true cross-book fair. Flags compare each venue's **raw take price**
    (``implied_prob``) to that median. A flagged cell is a direct bettable
    +EV signal.

    Flag rules (asymmetric by venue)
    --------------------------------
      * **Kalshi**: ``abs(implied_prob - median) >= threshold``. Both sides
        of a Kalshi binary are purchasable (buy YES at ``yes_ask``, buy NO
        at ``100 - yes_bid``), so an edge in either direction is bettable.
      * **All other books** (DK, FD, Bookmaker, Wagerzon, Hoop88, BetOnline):
        ``(median - implied_prob) >= threshold``. NFL-draft markets on these
        sportsbooks are YES-only futures — if a book prices YES above the
        consensus fair, the bettor cannot take the other side, so flagging
        an overpriced YES is noise. Only "YES is cheap on this book" fires.

    The signed delta exposed by ``ev_candidates``
    (``implied_prob - median``) still preserves direction:
      * delta > 0: price above fair -> YES is overpriced -> bet NO (Kalshi
        only; non-Kalshi books never reach this case because the flag
        suppresses it).
      * delta < 0: price below fair -> YES is underpriced -> bet YES.

    Rows with ``implied_prob`` NULL (one-sided Kalshi with no ask)
    participate in the median via ``devig_prob`` but can't be flagged.

    A market with only one posting venue has no median to compare against,
    so it gets no flags — still listed for completeness.

    Returns ``QueryLocked`` sentinel instead of raising on DuckDB lock
    contention; see module docstring.
    """
```

- [ ] **Step 3: Run the three new tests from Task 2 and confirm they all pass**

Run:
```bash
cd /Users/callancapitolo/NFLWork/.worktrees/asymmetric-flags && \
/Users/callancapitolo/NFLWork/kalshi_draft/venv/bin/pytest \
  nfl_draft/tests/unit/test_cross_book_grid.py::test_non_kalshi_below_median_is_flagged \
  nfl_draft/tests/unit/test_cross_book_grid.py::test_non_kalshi_above_median_is_not_flagged \
  nfl_draft/tests/unit/test_cross_book_grid.py::test_kalshi_above_median_still_flags -v
```

Expected: all three PASS.

- [ ] **Step 4: Commit the green-phase implementation**

```bash
cd /Users/callancapitolo/NFLWork/.worktrees/asymmetric-flags && \
git add nfl_draft/lib/queries.py && \
git commit -m "$(cat <<'EOF'
feat(nfl_draft): asymmetric outlier flag rule in cross_book_grid

Kalshi keeps the symmetric |take - median| >= threshold rule because both
sides of a binary are purchasable. Every other venue (DK/FD/Bookmaker/
Wagerzon/Hoop88/BetOnline) now only flags when median - take >= threshold,
i.e. YES is cheap on that book -> bettable YES wager. Above-median
non-Kalshi flags were unactionable on YES-only draft futures markets.

Co-Authored-By: Claude Opus 4.7 (1M context) <noreply@anthropic.com>
EOF
)"
```

Expected: one new commit.

---

### Task 4: Repair tests that relied on the old symmetric behavior

**Files:**
- Modify: `nfl_draft/tests/unit/test_cross_book_grid.py`
- Modify: `nfl_draft/tests/integration/test_dashboard_queries.py`

Two existing tests fail under the new rule because they encoded the symmetric contract. We repair both so the full suite is green.

- [ ] **Step 1: Run the full unit test file to see what's red**

Run:
```bash
cd /Users/callancapitolo/NFLWork/.worktrees/asymmetric-flags && \
/Users/callancapitolo/NFLWork/kalshi_draft/venv/bin/pytest \
  nfl_draft/tests/unit/test_cross_book_grid.py -v
```

Expected failures:
- `test_cross_book_grid_flag_uses_implied_vs_median_fair_all_books` — asserts `m["flags"]["bookmaker"] is True` when bookmaker is 11.5pp *above* median. Under the new rule this must be `False`.

All other tests in the file should pass.

- [ ] **Step 2: Repair `test_cross_book_grid_flag_uses_implied_vs_median_fair_all_books`**

In `nfl_draft/tests/unit/test_cross_book_grid.py`, replace the entire existing test (lines 42-60, docstring and body) with:

```python
def test_cross_book_grid_flag_uses_implied_vs_median_fair_all_books(fresh_db):
    """Flag fires when a book's implied_prob is >= threshold_pp below the
    consensus median (bettable YES edge). Kalshi keeps a symmetric rule;
    every other book only flags on the below-median side because their
    draft markets are YES-only futures."""
    # Median fair should be ~0.42; bookmaker implied 0.30 -> 12pp *below*.
    _insert_odds("m1", "draftkings", 100, 0.500, 0.420, datetime.now())
    _insert_odds("m1", "fanduel",    100, 0.500, 0.417, datetime.now())
    _insert_odds("m1", "bookmaker", +233, 0.300, 0.305, datetime.now())
    _insert_odds("m1", "kalshi",    +108, 0.480, 0.420, datetime.now())
    rows = q.cross_book_grid(threshold_pp=10)
    m = rows[0]
    # Median fair is ~0.420
    assert abs(m["median"] - 0.420) < 0.01
    # Bookmaker: median - take = 0.420 - 0.300 = 12pp >= 10pp -> flag
    assert m["flags"]["bookmaker"] is True
    # DK/FD: median - take = 0.420 - 0.500 = -8pp (book is above median)
    # -> non-Kalshi above-median -> no flag
    assert m["flags"]["draftkings"] is False
    assert m["flags"]["fanduel"] is False
    # Kalshi: |0.480 - 0.420| = 6pp < 10pp -> no flag
    assert m["flags"]["kalshi"] is False
```

- [ ] **Step 3: Run the unit test file again and confirm all green**

Run:
```bash
cd /Users/callancapitolo/NFLWork/.worktrees/asymmetric-flags && \
/Users/callancapitolo/NFLWork/kalshi_draft/venv/bin/pytest \
  nfl_draft/tests/unit/test_cross_book_grid.py -v
```

Expected: every test in the file passes.

- [ ] **Step 4: Run the integration suite to see what's red there**

Run:
```bash
cd /Users/callancapitolo/NFLWork/.worktrees/asymmetric-flags && \
/Users/callancapitolo/NFLWork/kalshi_draft/venv/bin/pytest \
  nfl_draft/tests/integration/test_dashboard_queries.py -v
```

Expected failures:
- `test_ev_candidates_ranks_by_abs_delta` — fixture inserts `("bm", 0.65)` on market `m2` with others at 0.50. Under the new rule, `bm` is a non-Kalshi book above the median, so it no longer flags. The test expects two +EV rows (one per market) and breaks when only `m1` produces a row.

`test_cross_book_grid_query_outlier_flags` should still pass (its DK row at 0.30 is below the median, which flags under both rules).

- [ ] **Step 5: Repair `test_ev_candidates_ranks_by_abs_delta`**

In `nfl_draft/tests/integration/test_dashboard_queries.py`, replace the existing test (lines 88-108) with:

```python
def test_ev_candidates_ranks_by_abs_delta(monkeypatch, tmp_path):
    from nfl_draft.lib import db as db_module
    monkeypatch.setattr(db_module, "DB_PATH", tmp_path / "test.duckdb")
    db_module.init_schema()
    from datetime import datetime
    from nfl_draft.lib.db import write_connection
    now = datetime.now()
    with write_connection() as con:
        con.execute("INSERT INTO draft_markets (market_id, market_type) VALUES ('m1', 'prop')")
        con.execute("INSERT INTO draft_markets (market_id, market_type) VALUES ('m2', 'prop')")
        # m1: sportsbook (bm) 25pp below median -> bettable YES -> flags.
        for book, prob in [("kalshi", 0.50), ("dk", 0.50), ("fd", 0.50), ("bm", 0.25)]:
            con.execute("INSERT INTO draft_odds VALUES ('m1', ?, 100, ?, ?, ?)", [book, prob, prob, now])
        # m2: kalshi 15pp above median -> bettable NO on the binary -> flags
        # (Kalshi remains symmetric). Non-Kalshi books here are at median, so
        # no sportsbook-side contribution.
        for book, prob in [("kalshi", 0.65), ("dk", 0.50), ("fd", 0.50), ("bm", 0.50)]:
            con.execute("INSERT INTO draft_odds VALUES ('m2', ?, 100, ?, ?, ?)", [book, prob, prob, now])
    from nfl_draft.lib.queries import ev_candidates
    rows = ev_candidates(threshold_pp=10)
    assert len(rows) == 2
    # Ranked by |delta| desc: m1 (|−0.25|) then m2 (|+0.15|).
    assert rows[0]["market_id"] == "m1"
    assert rows[0]["book"] == "bm"
    assert rows[0]["delta"] == pytest.approx(-0.25)
    assert rows[1]["market_id"] == "m2"
    assert rows[1]["book"] == "kalshi"
    assert rows[1]["delta"] == pytest.approx(0.15)
```

- [ ] **Step 6: Run the integration file and confirm green**

Run:
```bash
cd /Users/callancapitolo/NFLWork/.worktrees/asymmetric-flags && \
/Users/callancapitolo/NFLWork/kalshi_draft/venv/bin/pytest \
  nfl_draft/tests/integration/test_dashboard_queries.py -v
```

Expected: every test in the file passes.

- [ ] **Step 7: Run the full nfl_draft test suite as a final sanity check**

Run:
```bash
cd /Users/callancapitolo/NFLWork/.worktrees/asymmetric-flags && \
/Users/callancapitolo/NFLWork/kalshi_draft/venv/bin/pytest nfl_draft/tests/ -v
```

Expected: all tests pass. No new warnings beyond what main already produces.

- [ ] **Step 8: Commit the test repairs**

```bash
cd /Users/callancapitolo/NFLWork/.worktrees/asymmetric-flags && \
git add nfl_draft/tests/unit/test_cross_book_grid.py \
         nfl_draft/tests/integration/test_dashboard_queries.py && \
git commit -m "$(cat <<'EOF'
test(nfl_draft): update existing flag tests for asymmetric rule

The unit test for uniform flag semantics and the +EV ordering integration
test both encoded the old symmetric contract (non-Kalshi above-median
flagged). Rewrite them so the non-Kalshi cases use below-median fixtures
(bettable YES edges), and use Kalshi for the above-median case in the
ordering test (since Kalshi stays symmetric).

Co-Authored-By: Claude Opus 4.7 (1M context) <noreply@anthropic.com>
EOF
)"
```

Expected: one new commit.

---

### Task 5: Update the Cross-Book Grid UI copy

**Files:**
- Modify: `kalshi_draft/app.py` (around lines 834-842)

- [ ] **Step 1: Replace the grid blurb `html.P(...)` text**

In `kalshi_draft/app.py`, find this `html.P` block inside `render_crossbook_grid`:

```python
            html.P(
                "Each cell shows the bettable price — American odds for "
                "sportsbooks, cents for Kalshi. The median column is the "
                "cross-venue devigged fair (the true probability with "
                "vig stripped). A flagged cell (⚑) means the book's "
                "raw take price differs from the fair median by at least "
                "the threshold — i.e. direct +EV. Positive delta → "
                "price is above fair → bet NO; negative → bet YES.",
                style={"color": COLORS["text_muted"], "fontSize": "0.85em"},
            ),
```

Replace it with:

```python
            html.P(
                "Each cell shows the bettable price — American odds for "
                "sportsbooks, cents for Kalshi. The median column is the "
                "cross-venue devigged fair (the true probability with "
                "vig stripped). A flagged cell (⚑) means there is a "
                "bettable +EV edge: for sportsbooks (DK, FD, Bookmaker, "
                "Wagerzon, Hoop88, BetOnline), the book's YES price is "
                "cheaper than the fair median by at least the threshold "
                "(bet YES). For Kalshi, the take price differs from fair "
                "in either direction by at least the threshold (bet YES "
                "if underpriced, NO if overpriced). Sportsbook futures "
                "are YES-only, so an overpriced YES there is noise and "
                "does not flag.",
                style={"color": COLORS["text_muted"], "fontSize": "0.85em"},
            ),
```

- [ ] **Step 2: Smoke-check that Python parses the module**

Run:
```bash
cd /Users/callancapitolo/NFLWork/.worktrees/asymmetric-flags && \
/Users/callancapitolo/NFLWork/kalshi_draft/venv/bin/python -c "import kalshi_draft.app"
```

Expected: no output, exit code 0. (If import has side effects that launch the server, stop after the import succeeds — do not leave a dashboard running. If import fails, fix the copy change.)

- [ ] **Step 3: Commit the UI copy change**

```bash
cd /Users/callancapitolo/NFLWork/.worktrees/asymmetric-flags && \
git add kalshi_draft/app.py && \
git commit -m "$(cat <<'EOF'
chore(app): explain asymmetric flag rule in grid blurb

Updates the Cross-Book Grid description so users know why non-Kalshi
venues only flag below-median (YES-only sportsbook futures) while Kalshi
flags both directions (bettable binary).

Co-Authored-By: Claude Opus 4.7 (1M context) <noreply@anthropic.com>
EOF
)"
```

---

### Task 6: Update the README

**Files:**
- Modify: `nfl_draft/README.md` (around lines 33-46)

- [ ] **Step 1: Replace the outlier-flag paragraph**

In `nfl_draft/README.md`, find this block:

```markdown
The outlier flag (⚑) fires when a venue's **raw take price** (what
you'd actually pay) differs from the median fair by at least
`threshold_pp`. A flagged cell is a direct +EV signal:

- Delta > 0 (price above fair) → YES is overpriced → bet NO.
- Delta < 0 (price below fair) → YES is underpriced → bet YES.

Hover a Kalshi cell to see raw buy / sell / last-trade prices plus the
source ticker.

**Note on sportsbook vig.** Every sportsbook cell sits ~1–3pp above
the median fair by construction (that's the vig tax). The default
threshold of 10pp filters that out; avoid dropping below ~5pp or you'll
flag normal vig as "edges."
```

Replace it with:

```markdown
The outlier flag (⚑) fires on a bettable +EV edge. The rule is
asymmetric by venue:

- **Kalshi** (two-sided binary, both YES and NO purchasable): flag fires
  whenever `|take − median_fair| ≥ threshold_pp`. Delta < 0 → buy YES;
  delta > 0 → buy NO.
- **Sportsbooks** (DK, FD, Bookmaker, Wagerzon, Hoop88, BetOnline —
  YES-only futures): flag fires only when `median_fair − take ≥
  threshold_pp` (YES is cheap on this book → bet YES). An overpriced
  YES on a sportsbook is unactionable (you can't take the other side
  on a draft futures market), so those cases are intentionally not
  flagged even when `|delta|` exceeds the threshold.

Hover a Kalshi cell to see raw buy / sell / last-trade prices plus the
source ticker.

**Note on sportsbook vig.** Every sportsbook cell sits ~1–3pp above
the median fair by construction (that's the vig tax). Under the
asymmetric rule the sportsbook flag only looks at the "below" side, so
vig cannot cause a sportsbook cell to flag — but keep the threshold at
~10pp anyway so a minor edge (or a stale line) from only one book
doesn't fire on normal noise.
```

- [ ] **Step 2: Commit the README update**

```bash
cd /Users/callancapitolo/NFLWork/.worktrees/asymmetric-flags && \
git add nfl_draft/README.md && \
git commit -m "$(cat <<'EOF'
docs(nfl_draft): document asymmetric outlier-flag rule

Kalshi stays symmetric (bettable binary). Sportsbooks now only flag
when YES is cheap on the book (bettable YES wager); above-median flags
were unactionable on YES-only draft futures.

Co-Authored-By: Claude Opus 4.7 (1M context) <noreply@anthropic.com>
EOF
)"
```

---

### Task 7: Pre-merge executive review, approval, merge, cleanup

**Files:**
- None beyond the diff

- [ ] **Step 1: Run the full nfl_draft test suite one more time**

Run:
```bash
cd /Users/callancapitolo/NFLWork/.worktrees/asymmetric-flags && \
/Users/callancapitolo/NFLWork/kalshi_draft/venv/bin/pytest nfl_draft/tests/ -v
```

Expected: all tests pass, no errors.

- [ ] **Step 2: Produce the full diff against main for executive review**

Run:
```bash
cd /Users/callancapitolo/NFLWork/.worktrees/asymmetric-flags && \
git log main..HEAD --oneline && \
git diff main..HEAD --stat && \
git diff main..HEAD
```

Review against the CLAUDE.md pre-merge checklist:
- **Data integrity**: No schema or ingest changes — only read-path flag logic.
- **Resource safety**: No new DB connections introduced.
- **Edge cases**: `take is None` still short-circuits to `flags[book] = False` (explicitly, via the early `continue`). Single-venue markets still skip the flag loop entirely.
- **Dead code**: No unused flags/functions/imports introduced. The old comment "Uniform rule:" is gone.
- **Security**: No secrets, no logging changes.

- [ ] **Step 3: Ask the user to approve the merge**

Report the commit list and diff summary to the user and ask for explicit merge approval. Do NOT merge without a "yes".

- [ ] **Step 4: After approval, merge to main via fast-forward**

```bash
git -C /Users/callancapitolo/NFLWork checkout main && \
git -C /Users/callancapitolo/NFLWork merge --ff-only feature/asymmetric-outlier-flags && \
git -C /Users/callancapitolo/NFLWork log --oneline -6
```

Expected: fast-forward succeeds; the worktree commits now show on main.

- [ ] **Step 5: Remove the worktree and branch**

```bash
git -C /Users/callancapitolo/NFLWork worktree remove /Users/callancapitolo/NFLWork/.worktrees/asymmetric-flags && \
git -C /Users/callancapitolo/NFLWork branch -d feature/asymmetric-outlier-flags && \
git -C /Users/callancapitolo/NFLWork worktree list
```

Expected: worktree entry gone; branch deleted.

- [ ] **Step 6: Final sanity run from main**

```bash
cd /Users/callancapitolo/NFLWork && \
/Users/callancapitolo/NFLWork/kalshi_draft/venv/bin/pytest nfl_draft/tests/ -v
```

Expected: all tests pass from the merged main tree.
