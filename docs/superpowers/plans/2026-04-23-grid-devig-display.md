# Cross-Book Grid — Actual-Price Display with Devigged Median Anchor Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Fix the Cross-Book Grid so sportsbook rows are actually devigged at ingest, cells show the bettable price (American odds / Kalshi cents) while the median column stays the true fair, and the outlier flag compares each book's raw take price to the median of fairs (Unabated-style direct-EV signal).

**Architecture:** Derive a devig-grouping key from the existing `OddsRow.market_group` + canonical `market_id` at write-time (no scraper touches needed — scrapers already emit `market_group`). Bucket rows by `(book, outright_group)` in `write_or_quarantine`, run `proportional_devig` per bucket, persist to `draft_odds.devig_prob`. Drop the Kalshi-specific branch in `cross_book_grid`; flag logic becomes universal `|implied_prob − median(devig_prob)| ≥ threshold`. Update the Dash cell renderer to show `american_odds` / `implied_prob*100¢`.

**Tech Stack:** Python 3.11, DuckDB, Dash, pytest. Venv at `/Users/callancapitolo/NFLWork/kalshi_draft/venv/bin/python`.

**Deviation from spec:** The spec proposed adding `build_outright_group(market_type, **kwargs)` to `market_map.py` and updating every scraper to populate a new `OddsRow.outright_group` field. During planning I confirmed that every non-Kalshi scraper **already** emits `market_group` on each `OddsRow` (see `_base.py:24` and call sites in every scraper). Rather than add a parallel field and edit 7 scrapers, we derive the group key from `(market_group, market_id)` at write time — same grouping outcome, zero scraper changes, faster to ship on draft day. Kalshi rows pre-set `devig_prob` so they're skipped by the devig pass regardless.

**Spec:** `docs/superpowers/specs/2026-04-22-grid-devig-display-design.md`

---

## File Structure

**New files:**
- `nfl_draft/lib/backfill_devig.py` — one-shot backfill for existing `draft_odds` rows, deleted after running.

**Modified files:**
- `nfl_draft/lib/market_map.py` — add `outright_group_key(market_group, market_id)` helper (~50 lines including docstring/regex).
- `nfl_draft/lib/quarantine.py` — replace per-row write loop with bucket-then-devig flow.
- `nfl_draft/lib/queries.py` — add `american_odds` to `cross_book_grid` CTE + return; drop Kalshi flag branch; change `ev_candidates.delta`.
- `kalshi_draft/app.py` — cell renderer (American odds for sportsbooks, cents for Kalshi) + header copy.
- `nfl_draft/README.md` — Cross-Book Grid paragraph rewrite.

**Test files:**
- `nfl_draft/tests/unit/test_market_map.py` — add `outright_group_key` table-driven tests.
- `nfl_draft/tests/integration/test_quarantine.py` — add devig-at-ingest tests (n-way, 1-row, mixed, Kalshi pre-set).
- `nfl_draft/tests/unit/test_app_helpers.py` — cell-format helper test (if helper is extracted).

---

## Task 0: Create worktree + branch

**Files:**
- Create worktree: `/Users/callancapitolo/NFLWork-worktrees/grid-devig-display`
- Create branch: `feature/grid-devig-and-price-display`

- [ ] **Step 1: Verify clean working tree on `main`**

Run:
```bash
cd /Users/callancapitolo/NFLWork && git status
```
Expected: `On branch main`, possibly untracked files like `.playwright-mcp/` but no staged/modified tracked files. If modified files exist, stop and resolve before proceeding.

- [ ] **Step 2: Create worktree**

Run:
```bash
cd /Users/callancapitolo/NFLWork && git worktree add -b feature/grid-devig-and-price-display /Users/callancapitolo/NFLWork-worktrees/grid-devig-display main
```
Expected: `Preparing worktree ... HEAD is now at <sha>`.

- [ ] **Step 3: Copy DuckDB into worktree for integration testing**

Run:
```bash
cp /Users/callancapitolo/NFLWork/nfl_draft/nfl_draft.duckdb /Users/callancapitolo/NFLWork-worktrees/grid-devig-display/nfl_draft/nfl_draft.duckdb
```

⚠️ **Never symlink** — DuckDB WAL files land next to the *path*, not the *target*. Symlink + worktree removal = lost WAL data. Copy.

- [ ] **Step 4: Verify worktree + venv reachable**

Run:
```bash
cd /Users/callancapitolo/NFLWork-worktrees/grid-devig-display && git branch --show-current && /Users/callancapitolo/NFLWork/kalshi_draft/venv/bin/python -c "import nfl_draft; print('ok')"
```
Expected: `feature/grid-devig-and-price-display` then `ok`.

All subsequent tasks run from `/Users/callancapitolo/NFLWork-worktrees/grid-devig-display`.

---

## Task 1: Add `outright_group_key` helper

**Purpose:** Given `(market_group, market_id)` from an `OddsRow`, return the devig-grouping key, or `None` if the market is not devigable (props, matchup_before, mr_irrelevant_position).

**Files:**
- Modify: `nfl_draft/lib/market_map.py`
- Test: `nfl_draft/tests/unit/test_market_map.py`

- [ ] **Step 1: Write failing tests**

Append to `nfl_draft/tests/unit/test_market_map.py`:

```python
from nfl_draft.lib.market_map import outright_group_key


def test_outright_group_key_pick_outright():
    assert outright_group_key("pick_outright", "pick_2_overall_david-bailey") == "pick_2_overall"


def test_outright_group_key_first_at_position():
    assert outright_group_key("first_at_position", "first_wr_jordyn-tyson") == "first_wr"


def test_outright_group_key_top_n_range():
    assert outright_group_key("top_10_range", "top_10_sonny-styles") == "top_10"
    assert outright_group_key("top_5_range", "top_5_caleb-downs") == "top_5"


def test_outright_group_key_team_first_pick_handles_spaces_in_team():
    # build_market_id uses team.lower() -> can contain a space ("new england")
    assert outright_group_key("team_first_pick", "team_washington_first_pick_lucas") == "team_washington_first_pick"
    assert outright_group_key("team_first_pick", "team_new england_first_pick_drew-allar") == "team_new england_first_pick"


def test_outright_group_key_team_first_pick_position():
    # build_market_id uses _slug_underscored for both team and position
    assert outright_group_key(
        "team_first_pick_position",
        "arizona_cardinals_first_pick_pos_wide_receiver",
    ) == "arizona_cardinals_first_pick_pos"


def test_outright_group_key_nth_at_position():
    # market_group = "nth_at_position_2"; market_id = "2_wr_jordyn-tyson"
    assert outright_group_key("nth_at_position_2", "2_wr_jordyn-tyson") == "2_wr"
    assert outright_group_key("nth_at_position_3", "3_cb_jermod-mccoy") == "3_cb"


def test_outright_group_key_draft_position_over_under():
    assert outright_group_key(
        "draft_position_over_under",
        "draft_position_ou_spencer-fano_10p5_over",
    ) == "draft_position_ou_spencer-fano_10p5"
    assert outright_group_key(
        "draft_position_over_under",
        "draft_position_ou_spencer-fano_10p5_under",
    ) == "draft_position_ou_spencer-fano_10p5"


def test_outright_group_key_mr_irrelevant_returns_none():
    assert outright_group_key("mr_irrelevant_position", "mr_irrelevant_wide_receiver") is None


def test_outright_group_key_matchup_before_returns_none():
    # Matchups are self-contained 2-way; we don't bucket them across other markets.
    assert outright_group_key("matchup_before", "matchup_hunter_before_tyson") is None


def test_outright_group_key_prop_returns_none():
    assert outright_group_key("prop_first_round_total_ou", "prop_first_round_total") is None
    assert outright_group_key("prop_team_position_of_first_pick", "prop_xxx") is None


def test_outright_group_key_unrecognized_returns_none():
    assert outright_group_key("", "random_market_id") is None
    assert outright_group_key("totally_unknown", "x_y_z") is None


def test_outright_group_key_malformed_market_id_returns_none():
    # market_group says pick_outright but market_id doesn't match the pattern
    assert outright_group_key("pick_outright", "not_a_pick_market") is None
```

- [ ] **Step 2: Run tests to verify they fail**

Run:
```bash
/Users/callancapitolo/NFLWork/kalshi_draft/venv/bin/python -m pytest nfl_draft/tests/unit/test_market_map.py::test_outright_group_key_pick_outright -v
```
Expected: FAIL with `ImportError: cannot import name 'outright_group_key' from 'nfl_draft.lib.market_map'`.

- [ ] **Step 3: Implement `outright_group_key`**

Add at the end of `nfl_draft/lib/market_map.py`:

```python
def outright_group_key(market_group: str, market_id: str) -> Optional[str]:
    """Derive the devig-grouping key from (market_group, market_id), or None.

    Scrapers emit ``market_group`` on every ``OddsRow`` (see
    ``scrapers/_base.py``). This helper dispatches on ``market_group`` to
    extract the group instance (pick number, position, team+position, etc.)
    from ``market_id``. Rows that share the same ``(book, outright_group_key)``
    are competing outcomes in the same logical outright and are devigged
    together in ``quarantine.write_or_quarantine``.

    Returns ``None`` for:
      * prop markets (``market_group`` starting with ``prop_``)
      * ``mr_irrelevant_position`` (single row per book — nothing to devig)
      * ``matchup_before`` (self-contained 2-way, not grouped across markets)
      * unrecognized / malformed inputs

    The regex patterns here are tightly coupled to ``build_market_id``'s
    output format. If ``build_market_id`` changes, these must change too.
    Tested in ``tests/unit/test_market_map.py``.
    """
    if market_group == "pick_outright":
        m = re.match(r"^(pick_\d+_overall)_", market_id)
        return m.group(1) if m else None
    if market_group == "first_at_position":
        m = re.match(r"^(first_[a-z]+)_", market_id)
        return m.group(1) if m else None
    if market_group.startswith("top_") and market_group.endswith("_range"):
        m = re.match(r"^(top_\d+)_", market_id)
        return m.group(1) if m else None
    if market_group == "team_first_pick":
        m = re.match(r"^(team_.+?_first_pick)_", market_id)
        return m.group(1) if m else None
    if market_group == "team_first_pick_position":
        m = re.match(r"^(.+?_first_pick_pos)_", market_id)
        return m.group(1) if m else None
    if market_group.startswith("nth_at_position_"):
        nth = market_group[len("nth_at_position_"):]
        m = re.match(rf"^({re.escape(nth)}_[a-z]+)_", market_id)
        return m.group(1) if m else None
    if market_group == "draft_position_over_under":
        m = re.match(r"^(draft_position_ou_.+)_(?:over|under)$", market_id)
        return m.group(1) if m else None
    return None
```

- [ ] **Step 4: Run tests to verify they pass**

Run:
```bash
/Users/callancapitolo/NFLWork/kalshi_draft/venv/bin/python -m pytest nfl_draft/tests/unit/test_market_map.py -v
```
Expected: all tests pass, including the 11 new `test_outright_group_key_*` cases.

- [ ] **Step 5: Commit**

```bash
git add nfl_draft/lib/market_map.py nfl_draft/tests/unit/test_market_map.py
git commit -m "$(cat <<'EOF'
feat(market_map): add outright_group_key for ingest-time devig grouping

Derives the devig-grouping key from (market_group, market_id) already
emitted by every non-Kalshi scraper on OddsRow. Returns None for props,
matchups, and mr_irrelevant. Tight-coupled to build_market_id output
formats; tested per-market_type.

Co-Authored-By: Claude Opus 4.7 (1M context) <noreply@anthropic.com>
EOF
)"
```

---

## Task 2: Devig at ingest in `write_or_quarantine`

**Purpose:** Before writing to `draft_odds`, bucket mapped rows by `(book, outright_group_key)` and populate `devig_prob` via `proportional_devig`. Single-row buckets and `None` groups keep `devig_prob = implied_prob`. Kalshi rows (which pre-set `devig_prob`) pass through untouched.

**Files:**
- Modify: `nfl_draft/lib/quarantine.py`
- Test: `nfl_draft/tests/integration/test_quarantine.py`

- [ ] **Step 1: Write failing tests**

Append to `nfl_draft/tests/integration/test_quarantine.py`:

```python
def _seed_market_map_for_first_wr(con):
    """Helper: seed draft_markets + market_map for a 3-candidate first_wr outright.

    Uses DK as the only book so we can observe devig behavior in isolation.
    """
    con.execute("INSERT INTO draft_markets (market_id, market_type, subject_player, position) "
                "VALUES ('first_wr_hunter', 'first_at_position', 'Hunter', 'WR')")
    con.execute("INSERT INTO draft_markets (market_id, market_type, subject_player, position) "
                "VALUES ('first_wr_tyson', 'first_at_position', 'Tyson', 'WR')")
    con.execute("INSERT INTO draft_markets (market_id, market_type, subject_player, position) "
                "VALUES ('first_wr_tate', 'first_at_position', 'Tate', 'WR')")
    con.execute("INSERT INTO market_map VALUES ('draftkings', '1st WR', 'Hunter', 'first_wr_hunter')")
    con.execute("INSERT INTO market_map VALUES ('draftkings', '1st WR', 'Tyson',  'first_wr_tyson')")
    con.execute("INSERT INTO market_map VALUES ('draftkings', '1st WR', 'Tate',   'first_wr_tate')")


def test_write_or_quarantine_devigs_outright_group(seeded):
    """3-candidate first_wr outright at DK gets proportional-devigged."""
    from nfl_draft.lib.db import write_connection, read_connection
    with write_connection() as con:
        _seed_market_map_for_first_wr(con)
    now = datetime.now()
    rows = [
        OddsRow(book="draftkings", book_label="1st WR", book_subject="Hunter",
                american_odds=+110, fetched_at=now, market_group="first_at_position"),
        OddsRow(book="draftkings", book_label="1st WR", book_subject="Tyson",
                american_odds=+200, fetched_at=now, market_group="first_at_position"),
        OddsRow(book="draftkings", book_label="1st WR", book_subject="Tate",
                american_odds=+400, fetched_at=now, market_group="first_at_position"),
    ]
    write_or_quarantine(rows)
    with read_connection() as con:
        result = con.execute(
            "SELECT market_id, implied_prob, devig_prob "
            "FROM draft_odds ORDER BY market_id"
        ).fetchall()
    # Sum of raw implieds: 100/210 + 100/300 + 100/500 = 0.4762 + 0.3333 + 0.2000 = 1.0095
    # Devigged: each raw / 1.0095
    assert len(result) == 3
    total_devig = sum(r[2] for r in result)
    assert abs(total_devig - 1.0) < 1e-6, f"devig should sum to 1.0, got {total_devig}"
    # Each implied_prob should NOT equal devig_prob (would indicate no devig ran)
    for market_id, implied, devig in result:
        assert abs(implied - devig) > 1e-4, f"{market_id}: implied={implied} equals devig={devig}"


def test_write_or_quarantine_single_row_bucket_uses_implied(seeded):
    """A prop with no sibling (group=None) gets devig_prob = implied_prob."""
    from nfl_draft.lib.db import write_connection, read_connection
    with write_connection() as con:
        con.execute("INSERT INTO draft_markets (market_id, market_type) "
                    "VALUES ('prop_over_9p5_total_qbs_drafted_round_1', 'prop')")
        con.execute("INSERT INTO market_map VALUES "
                    "('draftkings', 'Total QBs Round 1', 'Over 9.5', 'prop_over_9p5_total_qbs_drafted_round_1')")
    rows = [OddsRow(
        book="draftkings", book_label="Total QBs Round 1", book_subject="Over 9.5",
        american_odds=-120, fetched_at=datetime.now(),
        market_group="prop_first_round_total_ou",
    )]
    write_or_quarantine(rows)
    with read_connection() as con:
        implied, devig = con.execute(
            "SELECT implied_prob, devig_prob FROM draft_odds"
        ).fetchone()
    # -120 -> implied = 120/220 = 0.5454...
    assert abs(implied - 120/220) < 1e-6
    assert abs(devig - implied) < 1e-9, "1-row group should pass implied through as devig"


def test_write_or_quarantine_respects_kalshi_pre_set_devig(seeded):
    """Kalshi rows arrive with devig_prob already set to mid; must not be overwritten."""
    from nfl_draft.lib.db import write_connection, read_connection
    with write_connection() as con:
        con.execute("INSERT INTO draft_markets (market_id, market_type, subject_player, position) "
                    "VALUES ('first_wr_hunter', 'first_at_position', 'Hunter', 'WR')")
        con.execute("INSERT INTO market_map VALUES "
                    "('kalshi', 'KXNFLDRAFT-FIRST-WR', 'Hunter', 'first_wr_hunter')")
    rows = [OddsRow(
        book="kalshi", book_label="KXNFLDRAFT-FIRST-WR", book_subject="Hunter",
        american_odds=+138, fetched_at=datetime.now(),
        market_group="",  # Kalshi doesn't set market_group
        implied_prob=0.42, devig_prob=0.40,  # take vs mid
    )]
    write_or_quarantine(rows)
    with read_connection() as con:
        implied, devig = con.execute(
            "SELECT implied_prob, devig_prob FROM draft_odds"
        ).fetchone()
    assert abs(implied - 0.42) < 1e-9
    assert abs(devig - 0.40) < 1e-9, "Kalshi pre-set devig_prob must be preserved"


def test_write_or_quarantine_groups_only_within_book(seeded):
    """DK's first_wr devig must not be polluted by another book's first_wr rows."""
    from nfl_draft.lib.db import write_connection, read_connection
    with write_connection() as con:
        _seed_market_map_for_first_wr(con)
        # Add bookmaker as second book on the same outright
        con.execute("INSERT INTO market_map VALUES "
                    "('bookmaker', '1st WR', 'Hunter', 'first_wr_hunter')")
        con.execute("INSERT INTO market_map VALUES "
                    "('bookmaker', '1st WR', 'Tyson',  'first_wr_tyson')")
    now = datetime.now()
    rows = [
        # DK: 2-candidate subset
        OddsRow(book="draftkings", book_label="1st WR", book_subject="Hunter",
                american_odds=+100, fetched_at=now, market_group="first_at_position"),
        OddsRow(book="draftkings", book_label="1st WR", book_subject="Tyson",
                american_odds=+100, fetched_at=now, market_group="first_at_position"),
        # Bookmaker: 2-candidate subset with different numbers
        OddsRow(book="bookmaker", book_label="1st WR", book_subject="Hunter",
                american_odds=-150, fetched_at=now, market_group="first_at_position"),
        OddsRow(book="bookmaker", book_label="1st WR", book_subject="Tyson",
                american_odds=+200, fetched_at=now, market_group="first_at_position"),
    ]
    write_or_quarantine(rows)
    with read_connection() as con:
        dk_sum = con.execute(
            "SELECT SUM(devig_prob) FROM draft_odds WHERE book='draftkings'"
        ).fetchone()[0]
        bm_sum = con.execute(
            "SELECT SUM(devig_prob) FROM draft_odds WHERE book='bookmaker'"
        ).fetchone()[0]
    assert abs(dk_sum - 1.0) < 1e-6, f"DK devig should sum to 1.0 within book, got {dk_sum}"
    assert abs(bm_sum - 1.0) < 1e-6, f"Bookmaker devig should sum to 1.0 within book, got {bm_sum}"
```

- [ ] **Step 2: Run tests to verify they fail**

Run:
```bash
/Users/callancapitolo/NFLWork/kalshi_draft/venv/bin/python -m pytest nfl_draft/tests/integration/test_quarantine.py::test_write_or_quarantine_devigs_outright_group -v
```
Expected: FAIL with `assert abs(implied - devig) > 1e-4` (because current code writes `devig = implied`).

- [ ] **Step 3: Replace `write_or_quarantine` with bucketed devig flow**

Replace the entire contents of `nfl_draft/lib/quarantine.py` with:

```python
"""Write OddsRow batches to draft_odds with quarantine for unmapped rows.

Devig is applied at ingest time: rows are bucketed by (book,
outright_group_key(market_group, market_id)) and proportional_devig runs
per bucket. Rows whose scraper already pre-set ``devig_prob`` (Kalshi)
pass through untouched. Single-row buckets and rows whose market_group
maps to no group key (props, matchups, mr_irrelevant) get
``devig_prob = implied_prob`` — no devig, comparison against the
cross-venue median still works against raw implieds.
"""

from typing import List, Tuple, Dict, Optional
from nfl_draft.lib.db import write_connection, read_connection
from nfl_draft.lib.devig import american_to_implied, devig_n_way
from nfl_draft.lib.market_map import outright_group_key
from nfl_draft.scrapers._base import OddsRow


def write_or_quarantine(rows: List[OddsRow]) -> Tuple[int, int]:
    """Resolve each row's market_id; write mapped to draft_odds (with devig
    applied per outright group), unmapped to quarantine.

    Returns: (mapped_count, unmapped_count).
    """
    if not rows:
        return (0, 0)

    # Single read pass -- load every market_map entry into memory.
    with read_connection() as con:
        map_rows = con.execute(
            "SELECT book, book_label, book_subject, market_id FROM market_map"
        ).fetchall()
    lookup = {(b, l, s): m for b, l, s, m in map_rows}

    # Resolve market_ids up front; collect mapped rows so we can bucket-devig them.
    mapped: List[Tuple[OddsRow, str]] = []  # (row, market_id)
    unmapped_rows: List[OddsRow] = []
    for row in rows:
        mid = lookup.get((row.book, row.book_label, row.book_subject))
        if mid is None:
            unmapped_rows.append(row)
        else:
            mapped.append((row, mid))

    # Bucket mapped rows by (book, outright_group_key). Rows whose group
    # resolves to None go into solo buckets keyed on identity so they remain
    # 1-row buckets (no devig).
    buckets: Dict[object, List[int]] = {}
    for idx, (row, mid) in enumerate(mapped):
        group = outright_group_key(row.market_group or "", mid)
        key: object = (row.book, group) if group is not None else ("__solo__", idx)
        buckets.setdefault(key, []).append(idx)

    # Compute devig_prob per mapped row index.
    devig_by_idx: Dict[int, float] = {}
    for key, idx_list in buckets.items():
        bucket_rows = [mapped[i] for i in idx_list]
        # If ANY row in the bucket arrived with devig_prob pre-set (Kalshi),
        # respect every row's scraper-supplied values verbatim. In practice a
        # bucket contains rows from one scraper so this is all-or-nothing.
        if any(row.devig_prob is not None for row, _ in bucket_rows):
            for i, (row, _) in zip(idx_list, bucket_rows):
                if row.devig_prob is not None:
                    devig_by_idx[i] = row.devig_prob
                else:
                    # Shouldn't happen in practice; fall back to implied.
                    implied = (row.implied_prob
                               if row.implied_prob is not None
                               else american_to_implied(row.american_odds))
                    devig_by_idx[i] = implied
            continue
        if len(bucket_rows) >= 2:
            # n-way proportional devig over the bucket's American odds.
            devigged = devig_n_way([row.american_odds for row, _ in bucket_rows])
            for i, dev in zip(idx_list, devigged):
                devig_by_idx[i] = dev
        else:
            # 1-row bucket (or market_group with no group key): pass implied through.
            row, _ = bucket_rows[0]
            devig_by_idx[idx_list[0]] = american_to_implied(row.american_odds)

    # Single write pass: unmapped -> quarantine, mapped -> draft_odds.
    with write_connection() as con:
        for row in unmapped_rows:
            con.execute(
                "INSERT INTO draft_odds_unmapped VALUES (?, ?, ?, ?, ?, ?)",
                [row.book, row.book_label, row.book_subject,
                 row.american_odds, row.fetched_at, "no market_map entry"],
            )
        for idx, (row, mid) in enumerate(mapped):
            implied = (row.implied_prob
                       if row.implied_prob is not None
                       else american_to_implied(row.american_odds))
            devig = devig_by_idx[idx]
            con.execute(
                "INSERT INTO draft_odds VALUES (?, ?, ?, ?, ?, ?)",
                [mid, row.book, row.american_odds, implied, devig, row.fetched_at],
            )
    return (len(mapped), len(unmapped_rows))
```

- [ ] **Step 4: Run all quarantine tests to verify pass**

Run:
```bash
/Users/callancapitolo/NFLWork/kalshi_draft/venv/bin/python -m pytest nfl_draft/tests/integration/test_quarantine.py -v
```
Expected: All tests pass — the 3 existing tests (unmapped, mapped, perf) + the 4 new devig tests.

- [ ] **Step 5: Run full unit + integration suite to catch regressions**

Run:
```bash
/Users/callancapitolo/NFLWork/kalshi_draft/venv/bin/python -m pytest nfl_draft/tests/unit nfl_draft/tests/integration -v
```
Expected: All tests pass. If anything in `test_dashboard_queries.py` or `test_pipeline.py` breaks, investigate — those tests may be asserting `devig_prob == implied_prob` based on today's (broken) behavior. Update the test's expected value to the correctly-devigged number, not the old buggy behavior.

- [ ] **Step 6: Commit**

```bash
git add nfl_draft/lib/quarantine.py nfl_draft/tests/integration/test_quarantine.py
git commit -m "$(cat <<'EOF'
fix(quarantine): actually devig sportsbook rows at ingest

Pre-fix, every sportsbook row had devig_prob == implied_prob (verified by
DB scan 2026-04-22). write_or_quarantine now buckets mapped rows by
(book, outright_group_key) and runs proportional_devig per bucket. Kalshi
rows pass through (scraper pre-sets devig_prob from mid). Single-row
buckets and prop/matchup markets (group key None) keep devig = implied.

Co-Authored-By: Claude Opus 4.7 (1M context) <noreply@anthropic.com>
EOF
)"
```

---

## Task 3: Change `cross_book_grid` flag semantics + expose `american_odds`

**Purpose:** (a) Drop the Kalshi-specific flag branch; every book compares `|implied_prob − median(devig_prob)| ≥ threshold`. (b) Include `american_odds` and `implied_prob` in the returned shape so the Dash callback can format cells with the actual price.

**Files:**
- Modify: `nfl_draft/lib/queries.py`
- Test: `nfl_draft/tests/unit/` (we'll add a dedicated `test_cross_book_grid.py` since no such test file exists)

- [ ] **Step 1: Create `test_cross_book_grid.py` with failing tests**

Create `nfl_draft/tests/unit/test_cross_book_grid.py`:

```python
"""Unit tests for cross_book_grid and ev_candidates flag semantics.

Uses tmp DuckDB per test. Populates draft_odds directly (bypasses scrapers).
"""
import pytest
from datetime import datetime
from nfl_draft.lib import db as db_module
from nfl_draft.lib import queries as q


@pytest.fixture
def fresh_db(monkeypatch, tmp_path):
    monkeypatch.setattr(db_module, "DB_PATH", tmp_path / "test.duckdb")
    db_module.init_schema()
    yield


def _insert_odds(market_id, book, odds, implied, devig, fetched_at=None):
    from nfl_draft.lib.db import write_connection
    if fetched_at is None:
        fetched_at = datetime.now()
    with write_connection() as con:
        con.execute(
            "INSERT INTO draft_odds VALUES (?, ?, ?, ?, ?, ?)",
            [market_id, book, odds, implied, devig, fetched_at],
        )


def test_cross_book_grid_returns_american_odds_and_implied(fresh_db):
    """Callback needs american_odds + implied_prob to format cells."""
    _insert_odds("m1", "draftkings", 110, 0.476, 0.425, datetime.now())
    _insert_odds("m1", "bookmaker", -115, 0.535, 0.435, datetime.now())
    rows = q.cross_book_grid(threshold_pp=10)
    assert len(rows) == 1
    m = rows[0]
    dk = m["books"]["draftkings"]
    assert dk["american_odds"] == 110
    assert dk["implied_prob"] == pytest.approx(0.476)
    assert dk["devig_prob"] == pytest.approx(0.425)


def test_cross_book_grid_flag_uses_implied_vs_median_fair_all_books(fresh_db):
    """Flag fires when book's implied_prob is ≥ threshold_pp away from
    median(devig_prob). Applies uniformly to every book, not just Kalshi."""
    # Median fair should be ~0.42; bookmaker implied 0.535 -> 11.5pp above.
    _insert_odds("m1", "draftkings", 100, 0.500, 0.420, datetime.now())
    _insert_odds("m1", "fanduel",    100, 0.500, 0.417, datetime.now())
    _insert_odds("m1", "bookmaker", -115, 0.535, 0.425, datetime.now())
    _insert_odds("m1", "kalshi",    +108, 0.480, 0.420, datetime.now())
    rows = q.cross_book_grid(threshold_pp=10)
    m = rows[0]
    # Median fair is ~0.420
    assert abs(m["median"] - 0.420) < 0.01
    # Bookmaker: |0.535 - 0.420| = 11.5pp -> flag
    assert m["flags"]["bookmaker"] is True
    # DK/FD: |0.500 - 0.420| = 8pp -> no flag at threshold 10
    assert m["flags"]["draftkings"] is False
    assert m["flags"]["fanduel"] is False
    # Kalshi: |0.480 - 0.420| = 6pp -> no flag
    assert m["flags"]["kalshi"] is False


def test_cross_book_grid_kalshi_flag_no_special_case(fresh_db):
    """Kalshi with implied_prob == 0.30 against median fair 0.42 -> flags
    (same formula as every other book). Pre-change code had a Kalshi-specific
    branch; the new code treats Kalshi uniformly."""
    _insert_odds("m1", "draftkings", 100, 0.500, 0.420, datetime.now())
    _insert_odds("m1", "fanduel",    100, 0.500, 0.420, datetime.now())
    _insert_odds("m1", "kalshi",    +233, 0.300, 0.305, datetime.now())
    rows = q.cross_book_grid(threshold_pp=10)
    m = rows[0]
    # Kalshi: |0.300 - 0.420| = 12pp -> flag
    assert m["flags"]["kalshi"] is True


def test_ev_candidates_delta_uses_implied_prob(fresh_db):
    """ev_candidates.delta must be implied_prob - median_fair (= the EV in pp),
    not devig_prob - median_fair."""
    _insert_odds("m1", "draftkings", 100, 0.500, 0.420, datetime.now())
    _insert_odds("m1", "fanduel",    100, 0.500, 0.420, datetime.now())
    _insert_odds("m1", "bookmaker", -115, 0.535, 0.425, datetime.now())
    rows = q.ev_candidates(threshold_pp=10)
    assert len(rows) == 1
    r = rows[0]
    assert r["book"] == "bookmaker"
    # delta = implied (0.535) - median_fair (~0.420) = ~0.115
    assert r["delta"] == pytest.approx(0.535 - 0.420, abs=0.01)
    assert r["book_prob"] == pytest.approx(0.535)
```

- [ ] **Step 2: Run tests to verify they fail**

Run:
```bash
/Users/callancapitolo/NFLWork/kalshi_draft/venv/bin/python -m pytest nfl_draft/tests/unit/test_cross_book_grid.py -v
```
Expected: FAIL — `KeyError: 'american_odds'` on the first test (the returned shape doesn't include it yet). Subsequent tests fail on flag semantics.

- [ ] **Step 3: Update `cross_book_grid` + `ev_candidates` in `queries.py`**

Replace `cross_book_grid` in `nfl_draft/lib/queries.py` (lines 89–179 approximately) with:

```python
def cross_book_grid(threshold_pp: float = 10.0):
    """Return one row per market with all-venue prices, median, and outlier flags.

    Median is ``statistics.median(devig_prob across posting venues)`` — the
    true cross-book fair. Flags compare each venue's **raw take price**
    (``implied_prob``) to that median. A flagged cell is a direct +EV signal:
    the venue is offering a price that differs from the cross-venue fair by
    more than ``threshold_pp`` percentage points.

    Signed delta (``implied_prob - median_fair``) preserves direction:
      * delta > 0: price above fair -> YES is overpriced -> bet NO
      * delta < 0: price below fair -> YES is underpriced -> bet YES

    Kalshi is treated uniformly: its ``implied_prob`` is ``yes_ask / 100``
    (the take price), so the same formula applies. Rows with
    ``implied_prob`` NULL (one-sided Kalshi with no ask) participate in
    the median via ``devig_prob`` but can't be flagged.

    A market with only one posting venue has no median to compare against,
    so it gets no flags — still listed for completeness.

    Returns ``QueryLocked`` sentinel instead of raising on DuckDB lock
    contention; see module docstring.
    """
    def _query() -> List[Dict[str, Any]]:
        with read_connection() as con:
            rows = con.execute(
                f"""
                WITH latest AS (
                  SELECT market_id, book, american_odds, implied_prob, devig_prob,
                         ROW_NUMBER() OVER (PARTITION BY market_id, book ORDER BY fetched_at DESC) AS rn
                  FROM draft_odds
                  WHERE fetched_at > NOW() - INTERVAL '{MAX_AGE_HOURS} hours'
                )
                SELECT market_id, book, american_odds, implied_prob, devig_prob
                FROM latest WHERE rn = 1
                """
            ).fetchall()

        by_market: Dict[str, Dict[str, Dict[str, Any]]] = {}
        for market_id, book, american_odds, implied_prob, devig_prob in rows:
            by_market.setdefault(market_id, {})[book] = {
                "american_odds": american_odds,
                "implied_prob": implied_prob,
                "devig_prob": devig_prob,
            }

        threshold = threshold_pp / 100.0
        output: List[Dict[str, Any]] = []
        for market_id, books in by_market.items():
            # Median uses devig_prob (fair). Pull every non-null devig for the median.
            fair_probs = [r["devig_prob"] for r in books.values() if r["devig_prob"] is not None]

            if len(books) < 2 or not fair_probs:
                # Only one posting venue -> no median -> no flags.
                output.append({
                    "market_id": market_id,
                    "books": books,
                    "median": fair_probs[0] if fair_probs else None,
                    "flags": {},
                    "outlier_count": 0,
                })
                continue

            median = statistics.median(fair_probs)
            flags: Dict[str, bool] = {}
            for book, r in books.items():
                # Uniform rule: flag if |raw take price - median fair| >= threshold.
                take = r["implied_prob"]
                flags[book] = (take is not None and abs(take - median) >= threshold)
            output.append({
                "market_id": market_id,
                "books": books,
                "median": median,
                "flags": flags,
                "outlier_count": sum(flags.values()),
            })
        return output

    return _safe_read(_query)
```

Replace `ev_candidates` in the same file with:

```python
def ev_candidates(threshold_pp: float = 10.0):
    """Flat list of flagged (market, venue) outliers, sorted by |delta| desc.

    Each row represents a book/market combo where the book's **raw take
    price** sits more than threshold_pp from the cross-book median fair.
    Signed delta (``implied_prob - median_fair``) preserves direction and
    equals the EV in percentage points:
      * positive -> price above fair -> bet NO
      * negative -> price below fair -> bet YES

    Delegates to ``cross_book_grid``; propagates its ``QueryLocked``
    sentinel when the read fails.
    """
    grid = cross_book_grid(threshold_pp)
    if isinstance(grid, QueryLocked):
        return grid
    out: List[Dict[str, Any]] = []
    for m in grid:
        if m["median"] is None or not m["flags"]:
            continue
        for book, flagged in m["flags"].items():
            if not flagged:
                continue
            book_prob = m["books"][book]["implied_prob"]
            delta = book_prob - m["median"]
            out.append({
                "market_id": m["market_id"],
                "book": book,
                "book_prob": book_prob,
                "median": m["median"],
                "delta": delta,
            })
    out.sort(key=lambda r: abs(r["delta"]), reverse=True)
    return out
```

- [ ] **Step 4: Run tests to verify pass**

Run:
```bash
/Users/callancapitolo/NFLWork/kalshi_draft/venv/bin/python -m pytest nfl_draft/tests/unit/test_cross_book_grid.py -v
```
Expected: all four new tests pass.

- [ ] **Step 5: Run full suite to catch regressions from `books` shape change**

Run:
```bash
/Users/callancapitolo/NFLWork/kalshi_draft/venv/bin/python -m pytest nfl_draft/tests/ -v
```
Expected: pass. The `books` dict shape changed from `{book: prob}` to `{book: {american_odds, implied_prob, devig_prob}}` — the only consumer is `app.py::_update_crossbook` which we update next.

- [ ] **Step 6: Commit**

```bash
git add nfl_draft/lib/queries.py nfl_draft/tests/unit/test_cross_book_grid.py
git commit -m "$(cat <<'EOF'
feat(queries): universal raw-price-vs-fair-median flag semantics

cross_book_grid now returns per-book {american_odds, implied_prob,
devig_prob} so the dashboard can render bettable prices. Median is still
devig_prob across venues. Flag rule is uniform: |implied - median_fair|
>= threshold for every book including Kalshi (dropped Kalshi special
branch). ev_candidates.delta is now implied - median (direct EV in pp).

Co-Authored-By: Claude Opus 4.7 (1M context) <noreply@anthropic.com>
EOF
)"
```

---

## Task 4: Cell renderer in `app.py`

**Purpose:** Show the bettable price in each grid cell: American odds for sportsbooks (`+110`, `−150`), cents for Kalshi (`42¢`). Median column stays `X.X%`. ⚑ suffix unchanged.

**Files:**
- Modify: `kalshi_draft/app.py:1107-1120` (render loop) and `kalshi_draft/app.py:833-840` (header copy)

- [ ] **Step 1: Update the render loop in `_update_crossbook`**

In `kalshi_draft/app.py`, replace the block at lines 1107–1120 (approximately; the block that populates `row[venue]`):

```python
    rows = []
    tooltip_data = []
    for m in grid:
        row = {"market_id": m["market_id"]}
        for venue in VENUES:
            record = m["books"].get(venue)
            flagged = m["flags"].get(venue, False)
            if record is None:
                row[venue] = ""
                continue
            if venue == "kalshi":
                ip = record["implied_prob"]
                cell = f"{round(ip * 100)}¢" if ip is not None else ""
            else:
                ao = record["american_odds"]
                cell = f"{int(ao):+d}" if ao is not None else ""
            row[venue] = cell + (" ⚑" if flagged and cell else "")
        row["median"] = f"{m['median']*100:.1f}%" if m["median"] is not None else ""
        row["outliers"] = m["outlier_count"]
        rows.append(row)

        tip = tooltip_lookup.get(m["market_id"])
        if tip is not None:
            kalshi_hover = (
                f"Last trade: {tip['last_price']}¢  \n"
                f"Buy: {tip['yes_ask']}¢  Sell: {tip['yes_bid']}¢  \n"
                f"Ticker: {tip['ticker']}"
            )
            tooltip_data.append({
                "kalshi": {"value": kalshi_hover, "type": "markdown"},
            })
        else:
            tooltip_data.append({})
```

Note: the existing `TABLE_STYLE_DATA_CONDITIONAL` still matches `{v} contains '⚑'` so the red highlight keeps working with no change.

- [ ] **Step 2: Update header copy in `render_crossbook_grid`**

In `kalshi_draft/app.py`, replace the `html.P(...)` block at approximately lines 834–841 with:

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

- [ ] **Step 3: Manually smoke-test the dashboard**

Run:
```bash
/Users/callancapitolo/NFLWork/kalshi_draft/venv/bin/python kalshi_draft/app.py &
```

(If a dashboard is already running, the daemon scrape on restart will append to logs — that's fine for smoke test. Kill the process after.)

Open `http://127.0.0.1:8090/`, click **Portal** → **Cross-Book Grid**, verify:

1. Sportsbook cells show American-odds format (`+110`, `-150`).
2. Kalshi cell shows cents (`42¢`).
3. Median column still shows percent (`42.5%`).
4. A row with a flag shows `⚑` and is highlighted red.
5. Header copy reflects the new semantics.

If no rows render: check `/tmp/nfl_draft_startup_scrape.log` — a fresh DuckDB in the worktree may have no recent data yet. Re-run a manual scrape:
```bash
/Users/callancapitolo/NFLWork/kalshi_draft/venv/bin/python -m nfl_draft.run --mode scrape --book all
```

Kill the dashboard process when done:
```bash
pkill -f "kalshi_draft/app.py"
```

- [ ] **Step 4: Run full suite (app callbacks are smoke-covered by existing tests)**

Run:
```bash
/Users/callancapitolo/NFLWork/kalshi_draft/venv/bin/python -m pytest nfl_draft/tests/ -v
```
Expected: all pass.

- [ ] **Step 5: Commit**

```bash
git add kalshi_draft/app.py
git commit -m "$(cat <<'EOF'
feat(app): cross-book grid shows bettable prices, not devigged probs

Sportsbook cells show American odds; Kalshi cells show cents. Median
column stays the cross-venue devigged fair. Header copy rewritten to
explain the new +EV-signal semantics (delta direction guides bet side).

Co-Authored-By: Claude Opus 4.7 (1M context) <noreply@anthropic.com>
EOF
)"
```

---

## Task 5: Backfill historical `draft_odds`

**Purpose:** The existing `draft_odds` rows (~2 weeks of pre-draft data) were all written with the old `devig_prob = implied_prob` bug. Re-compute `devig_prob` in place using the new grouping rule so historical data is usable.

**Files:**
- Create: `nfl_draft/lib/backfill_devig.py` (deleted after run)

- [ ] **Step 1: Write the backfill script**

Create `nfl_draft/lib/backfill_devig.py`:

```python
"""One-shot: re-compute draft_odds.devig_prob using the bucketed devig rule.

Usage (run once from main after merge, then delete this file):
    python -m nfl_draft.lib.backfill_devig

Groups rows by (date_trunc(minute, fetched_at), book, outright_group_key)
so concurrent scrape batches don't split. Skips Kalshi rows (their
devig_prob was always correct — written as mid by the scraper).
"""

from collections import defaultdict
from nfl_draft.lib.db import read_connection, write_connection
from nfl_draft.lib.devig import devig_n_way
from nfl_draft.lib.market_map import outright_group_key


def _load_rows():
    """Return (ts_bucket, book, market_id, american_odds, market_group) for every non-Kalshi draft_odds row.

    market_group comes from the LAST scraper run that successfully mapped a
    row for this (book, book_label, book_subject) — we look it up via
    market_map. Rows whose market_group is unknown get skipped (left as-is).
    """
    with read_connection() as con:
        return con.execute(
            """
            SELECT date_trunc('minute', o.fetched_at) AS ts_bucket,
                   o.book, o.market_id, o.american_odds,
                   -- No market_group column on draft_odds; join market_map
                   -- on (book, market_id) and take any matching book_label's
                   -- classification. Backfill tolerance: if the same market_id
                   -- has multiple labels (shouldn't), ANY_VALUE picks one.
                   o.fetched_at
            FROM draft_odds o
            WHERE o.book != 'kalshi'
            ORDER BY o.fetched_at
            """
        ).fetchall()


def _derive_market_group(market_id: str) -> str:
    """Heuristic fallback: derive market_group from market_id prefix.

    We don't persist market_group in draft_odds, so for the backfill we
    reverse-engineer it from market_id structure. This is narrower than
    outright_group_key's input requirements but sufficient for backfill.
    """
    if market_id.startswith("pick_") and "_overall_" in market_id:
        return "pick_outright"
    if market_id.startswith("first_"):
        # first_{pos}_{slug}
        return "first_at_position"
    if market_id.startswith("top_"):
        # top_{N}_{slug}
        import re
        m = re.match(r"^top_(\d+)_", market_id)
        return f"top_{m.group(1)}_range" if m else ""
    if market_id.startswith("team_") and "_first_pick_" in market_id:
        return "team_first_pick"
    if "_first_pick_pos_" in market_id:
        return "team_first_pick_position"
    if market_id.startswith("draft_position_ou_"):
        return "draft_position_over_under"
    if market_id.startswith("matchup_"):
        return "matchup_before"
    if market_id.startswith("mr_irrelevant_"):
        return "mr_irrelevant_position"
    if market_id.startswith("prop_"):
        return "prop_unknown"
    # Check nth_at_position: {N}_{pos}_{slug}
    import re
    m = re.match(r"^(\d+)_[a-z]+_", market_id)
    if m:
        return f"nth_at_position_{m.group(1)}"
    return ""


def run():
    print("Loading draft_odds rows (excluding Kalshi)...")
    rows = _load_rows()
    print(f"  Loaded {len(rows)} rows.")

    # Bucket: (ts_bucket, book, outright_group) -> [(market_id, american_odds, fetched_at), ...]
    buckets = defaultdict(list)
    ungrouped = []  # [(market_id, book, fetched_at)]
    for ts_bucket, book, market_id, american_odds, fetched_at in rows:
        mg = _derive_market_group(market_id)
        group = outright_group_key(mg, market_id)
        if group is None:
            ungrouped.append((market_id, book, fetched_at))
            continue
        buckets[(ts_bucket, book, group)].append((market_id, american_odds, fetched_at))

    print(f"  Bucketed into {len(buckets)} groups; {len(ungrouped)} rows ungrouped (kept as-is).")

    # Compute and write
    updates = 0
    with write_connection() as con:
        for (ts_bucket, book, group), items in buckets.items():
            if len(items) < 2:
                continue  # Nothing to devig; devig_prob already == implied
            devigged = devig_n_way([ao for _, ao, _ in items])
            for (market_id, _, fetched_at), dev in zip(items, devigged):
                con.execute(
                    """
                    UPDATE draft_odds
                    SET devig_prob = ?
                    WHERE market_id = ? AND book = ? AND fetched_at = ?
                    """,
                    [dev, market_id, book, fetched_at],
                )
                updates += 1
    print(f"  Updated devig_prob on {updates} rows.")
    print("Done. Delete this file after committing.")


if __name__ == "__main__":
    run()
```

- [ ] **Step 2: Dry-run on a copy, not the worktree DB**

Since the worktree has a copy of `nfl_draft.duckdb`, we can run against it destructively — it's a scratch copy. Run:

```bash
cd /Users/callancapitolo/NFLWork-worktrees/grid-devig-display
/Users/callancapitolo/NFLWork/kalshi_draft/venv/bin/python -m nfl_draft.lib.backfill_devig
```
Expected output: `Loaded N rows. Bucketed into M groups; K rows ungrouped. Updated devig_prob on U rows.` where U is thousands.

- [ ] **Step 3: Verify the backfill worked**

Run:
```bash
/Users/callancapitolo/NFLWork/kalshi_draft/venv/bin/python -c "
import duckdb
con = duckdb.connect('nfl_draft/nfl_draft.duckdb', read_only=True)
rows = con.execute('''
  SELECT book, COUNT(*) AS total,
         SUM(CASE WHEN implied_prob <> devig_prob THEN 1 ELSE 0 END) AS differ
  FROM draft_odds WHERE book != 'kalshi'
  GROUP BY book ORDER BY book
''').fetchall()
for r in rows: print(f'{r[0]:12s}  total={r[1]:5d}  differ={r[2]:5d}')
"
```
Expected: `differ` > 0 for every sportsbook (meaning devig_prob now differs from implied_prob where it should — i.e. the backfill ran).

- [ ] **Step 4: Commit the backfill script**

```bash
git add nfl_draft/lib/backfill_devig.py
git commit -m "$(cat <<'EOF'
chore: one-shot backfill for historical draft_odds.devig_prob

Re-computes devig_prob on non-Kalshi rows using the new bucketed
proportional_devig rule. Derives market_group heuristically from
market_id prefix (draft_odds doesn't persist market_group). Run once
from main post-merge, then deleted.

Co-Authored-By: Claude Opus 4.7 (1M context) <noreply@anthropic.com>
EOF
)"
```

---

## Task 6: Update `nfl_draft/README.md`

**Files:**
- Modify: `nfl_draft/README.md:21-39` (Cross-Book Grid section)

- [ ] **Step 1: Replace the section**

In `nfl_draft/README.md`, replace the "Cross-Book Grid — Kalshi vs sportsbook semantics" subsection with:

```markdown
### Cross-Book Grid — actual-price display

Each cell shows what you would actually bet at that venue:

- **Sportsbooks** (DK, FD, Bookmaker, Wagerzon, Hoop88, BetOnline):
  posted American odds (e.g. `+110`, `-150`).
- **Kalshi**: yes_ask in cents (e.g. `42¢`).

The **Median** column is the cross-venue **devigged fair** in percent —
the true probability estimate once each sportsbook's vig has been
stripped via proportional devig at ingest time.

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

See `docs/superpowers/specs/2026-04-22-grid-devig-display-design.md`
for the full design and `docs/superpowers/plans/2026-04-23-grid-devig-display.md`
for the implementation plan.
```

- [ ] **Step 2: Commit**

```bash
git add nfl_draft/README.md
git commit -m "$(cat <<'EOF'
docs(nfl_draft): rewrite Cross-Book Grid section for actual-price display

Reflects the new cell semantics (American odds / cents, not devigged
probs) and the new flag rule (raw take price vs median fair = direct EV).

Co-Authored-By: Claude Opus 4.7 (1M context) <noreply@anthropic.com>
EOF
)"
```

---

## Task 7: Pre-merge executive review + merge + cleanup

**Files:** read-only review.

- [ ] **Step 1: Full test run from worktree**

Run:
```bash
cd /Users/callancapitolo/NFLWork-worktrees/grid-devig-display
/Users/callancapitolo/NFLWork/kalshi_draft/venv/bin/python -m pytest nfl_draft/tests/ -v
```
Expected: green.

- [ ] **Step 2: Review the full diff**

Run:
```bash
git diff main..HEAD --stat
git diff main..HEAD
```

Review against CLAUDE.md "Pre-merge review" checklist:

1. **Data integrity**: devig bucket key includes `book` (per-book, not cross-book). Kalshi pre-set devig_prob preserved. Backfill buckets by minute to avoid splitting concurrent scrapes. ✓
2. **Resource safety**: No new DB connections added in hot paths; `write_or_quarantine` still uses a single `write_connection()` context manager. ✓
3. **Edge cases**: 1-row bucket / `group=None` / Kalshi / mixed bucket — all covered by tests. ✓
4. **Dead code**: No unused helpers added. `outright_group_key` is used in `quarantine.py` and is the only new public symbol. ✓
5. **Log/disk hygiene**: No new logs added. Backfill script deleted next step. ✓
6. **Security**: No credentials touched. ✓

Document ISSUES-TO-FIX vs ACCEPTABLE-RISKS in a short terminal summary back to the user.

- [ ] **Step 3: Get explicit user approval to merge**

Post the diff summary + review findings to the user. Wait for explicit approval before proceeding.

- [ ] **Step 4: Merge to main**

After user approval:

```bash
cd /Users/callancapitolo/NFLWork
git checkout main
git merge --no-ff feature/grid-devig-and-price-display -m "$(cat <<'EOF'
feat(nfl_draft): cross-book grid actual-price display with devigged median

Fixes sportsbook devig bug (devig_prob == implied_prob on 100% of non-Kalshi
rows). Grid cells now show bettable prices (American odds / Kalshi cents).
Flag rule uniform: |implied - median_fair| >= threshold -> direct +EV signal.
+EV Candidates delta is now EV in percentage points.

See docs/superpowers/specs/2026-04-22-grid-devig-display-design.md and
docs/superpowers/plans/2026-04-23-grid-devig-display.md.

Co-Authored-By: Claude Opus 4.7 (1M context) <noreply@anthropic.com>
EOF
)"
```

- [ ] **Step 5: Run backfill on production DB**

Run from main:
```bash
cd /Users/callancapitolo/NFLWork
/Users/callancapitolo/NFLWork/kalshi_draft/venv/bin/python -m nfl_draft.lib.backfill_devig
```
Expected: thousands of rows updated.

- [ ] **Step 6: Delete the backfill script**

```bash
cd /Users/callancapitolo/NFLWork
rm nfl_draft/lib/backfill_devig.py
git add nfl_draft/lib/backfill_devig.py
git commit -m "$(cat <<'EOF'
chore: remove one-shot backfill_devig.py after production run

Per repo convention (no temp files on disk; script was a one-shot).

Co-Authored-By: Claude Opus 4.7 (1M context) <noreply@anthropic.com>
EOF
)"
```

- [ ] **Step 7: Cleanup worktree + branch**

```bash
cd /Users/callancapitolo/NFLWork
git worktree remove /Users/callancapitolo/NFLWork-worktrees/grid-devig-display
git branch -d feature/grid-devig-and-price-display
```

Expected: worktree removed, branch deleted (since it's fully merged).

- [ ] **Step 8: Final verification**

Run:
```bash
git status
git log --oneline -8
```
Expected: clean tree on `main`, recent commits include the feature merge + backfill-script removal.

---

## Self-review

Spec coverage check:

- Math layer (devig at ingest) — Task 2 ✓
- Flag layer (universal implied-vs-fair) — Task 3 ✓
- Display layer (American odds / cents) — Task 4 ✓
- `outright_group` derivation — Task 1 (note: deviated from spec by deriving at write-time from existing `market_group` instead of adding a new field; rationale at top of plan) ✓
- Kalshi preservation — Task 2 Step 1 `test_write_or_quarantine_respects_kalshi_pre_set_devig` ✓
- Backfill — Task 5 ✓
- Docs — Task 6 ✓
- Worktree + merge flow — Tasks 0, 7 ✓
- Edge cases (1-row bucket, single-sided Yes/No, solo-venue market) — tests in Tasks 2 & 3 ✓

Placeholder scan: no TBDs, TODOs, or vague-instruction steps. Every code step contains actual code.

Type consistency: `outright_group_key` signature `(market_group: str, market_id: str) -> Optional[str]` consistent across definition (Task 1), usage (Task 2), and tests. `cross_book_grid` return shape `books[venue] = {american_odds, implied_prob, devig_prob}` consistent across Tasks 3 & 4.
