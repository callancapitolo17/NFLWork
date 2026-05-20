"""Regression tests for the SGP scraper shims and golden-baseline diffs.

The baseline CSV (tests/golden/dk_sgp_baseline.csv) was captured BEFORE the
dk_client refactor in Task 6. After refactor, running the scraper against the
live DK API must produce sgp_decimal values within tolerance of baseline.

Tolerance notes
---------------
The plan called for a 0.005 tolerance, but real DK markets move between
consecutive scrapes minutes apart (price refreshes, line moves). We observed
during this very task that a single game's FG spread+total combos drifted
~0.06-0.15 decimal odds between two runs taken ~10 minutes apart, while every
other combo on every other game matched exactly. That's the natural-movement
signature, not a code regression: a code regression would affect every game,
not one.

We also observed that the scraper filters events with start_time > now_utc,
so any games that tip off between baseline capture and the re-run legitimately
disappear from the output. A several-hour-old baseline can therefore lose
~half its rows for entirely innocent reasons.

The test compares only the intersection of (game_id, combo, period) keys
present in BOTH baseline and current. That intersection is the meaningful
regression surface — every still-active game's prices must match within
tolerance. Tolerance is split into two pieces:
  * 0.20 absolute per combo — absorbs natural line moves
  * <= 15% of intersected combos may drift > 0.20
    A real refactor regression would push almost every combo over 0.20.
The test additionally requires the intersection to be non-empty (else the
baseline is fully stale and the test would silently pass).


Schema
------
mlb_sgp_odds columns: game_id, combo, period, bookmaker, sgp_decimal,
sgp_american, fetch_time, source. The unique key for DK rows is
(game_id, combo, period).
"""
import sys
import pytest
import duckdb
import csv
import subprocess
from pathlib import Path


@pytest.mark.integration
def test_dk_sgp_output_matches_baseline():
    """Live API; re-runs the DK SGP scraper and diffs against the golden CSV."""
    mlb_sgp_dir = Path(__file__).resolve().parent.parent
    # Re-run the scraper after refactor. Use sys.executable so the subprocess
    # inherits the same venv that's running pytest (the bare name 'python' is
    # not always on PATH in test environments).
    result = subprocess.run(
        [sys.executable, "scraper_draftkings_sgp.py"],
        capture_output=True, text=True,
        cwd=mlb_sgp_dir,
    )
    assert result.returncode == 0, (
        f"Scraper failed: stdout={result.stdout!r} stderr={result.stderr!r}"
    )

    db_path = mlb_sgp_dir.parent / "Answer Keys" / "mlb_mm.duckdb"
    con = duckdb.connect(str(db_path), read_only=True)
    current = {
        (row[0], row[1], row[2]): float(row[3])
        for row in con.execute(
            "SELECT game_id, combo, period, sgp_decimal "
            "FROM mlb_sgp_odds WHERE bookmaker='draftkings'"
        ).fetchall()
    }
    con.close()

    baseline_path = Path(__file__).resolve().parent / "golden" / "dk_sgp_baseline.csv"
    baseline = {}
    with open(baseline_path) as f:
        reader = csv.DictReader(f)
        for row in reader:
            baseline[(row["game_id"], row["combo"], row["period"])] = float(
                row["sgp_decimal"]
            )

    shared = baseline.keys() & current.keys()
    assert shared, (
        "No (game_id, combo, period) keys in common between baseline and "
        "current run. Either the baseline is fully stale (all games tipped "
        f"off) or the refactor broke event discovery. baseline rows: "
        f"{len(baseline)}, current rows: {len(current)}"
    )
    drift = {
        k: (baseline[k], current[k])
        for k in shared
        if abs(baseline[k] - current[k]) > 0.20
    }
    drift_frac = len(drift) / len(shared)
    assert drift_frac <= 0.15, (
        f"Too many combos drifted >0.20 decimal odds: "
        f"{len(drift)}/{len(shared)} ({drift_frac:.0%}). "
        f"Sample: {dict(list(drift.items())[:5])}"
    )


@pytest.mark.integration
def test_fd_sgp_output_matches_baseline():
    """Live API; re-runs the FD SGP scraper and diffs against the golden CSV.

    Mirrors the DK pattern (Task 6): intersection-only comparison on
    (game_id, combo, period), 0.20 absolute tolerance, <=15% drift cap.
    Required because FD markets move naturally between consecutive scrapes
    and games tipping off legitimately drop out of the current set.
    """
    mlb_sgp_dir = Path(__file__).resolve().parent.parent
    result = subprocess.run(
        [sys.executable, "scraper_fanduel_sgp.py"],
        capture_output=True, text=True,
        cwd=mlb_sgp_dir,
    )
    assert result.returncode == 0, (
        f"Scraper failed: stdout={result.stdout!r} stderr={result.stderr!r}"
    )

    db_path = mlb_sgp_dir.parent / "Answer Keys" / "mlb_mm.duckdb"
    con = duckdb.connect(str(db_path), read_only=True)
    current = {
        (row[0], row[1], row[2]): float(row[3])
        for row in con.execute(
            "SELECT game_id, combo, period, sgp_decimal "
            "FROM mlb_sgp_odds WHERE bookmaker='fanduel'"
        ).fetchall()
    }
    con.close()

    baseline_path = Path(__file__).resolve().parent / "golden" / "fd_sgp_baseline.csv"
    baseline = {}
    with open(baseline_path) as f:
        reader = csv.DictReader(f)
        for row in reader:
            baseline[(row["game_id"], row["combo"], row["period"])] = float(
                row["sgp_decimal"]
            )

    shared = baseline.keys() & current.keys()
    assert shared, (
        "No (game_id, combo, period) keys in common between baseline and "
        "current run. Either the baseline is fully stale (all games tipped "
        f"off) or the refactor broke event discovery. baseline rows: "
        f"{len(baseline)}, current rows: {len(current)}"
    )
    drift = {
        k: (baseline[k], current[k])
        for k in shared
        if abs(baseline[k] - current[k]) > 0.20
    }
    drift_frac = len(drift) / len(shared)
    assert drift_frac <= 0.15, (
        f"Too many combos drifted >0.20 decimal odds: "
        f"{len(drift)}/{len(shared)} ({drift_frac:.0%}). "
        f"Sample: {dict(list(drift.items())[:5])}"
    )


def test_dk_shim_writes_new_schema_cols(tmp_path):
    """After shim refactor, scraper writes spread_line/total_line cols.

    Shim-shape check only — verifies the scraper:
      1) exits 0 when there are no target lines to scrape,
      2) creates the mlb_sgp_odds table with the post-Task-4 schema
         (spread_line + total_line columns present).

    Full pricing correctness is covered by the live-API golden tests
    above and by Task 28's manual smoke test.
    """
    import os
    import subprocess

    db = str(tmp_path / "shim_test.duckdb")
    # Empty parlay_lines table so the shim has no targets to scrape and
    # never hits DK's API. The shim should still ensure_table + exit 0.
    con = duckdb.connect(db)
    con.execute("""
        CREATE TABLE mlb_parlay_lines (
            game_id VARCHAR, home_team VARCHAR, away_team VARCHAR,
            commence_time VARCHAR, fg_spread DOUBLE, fg_total DOUBLE,
            f5_spread DOUBLE, f5_total DOUBLE
        )
    """)
    con.close()

    env = os.environ.copy()
    env["MLB_SGP_DB_PATH"] = db
    repo_root = (
        "/Users/callancapitolo/NFLWork/.claude/worktrees/"
        "kalshi-mlb-rfq-line-source-pivot"
    )
    result = subprocess.run(
        [sys.executable, "mlb_sgp/scraper_draftkings_sgp.py"],
        env=env, capture_output=True, timeout=60, cwd=repo_root,
    )
    assert result.returncode == 0, f"stderr: {result.stderr.decode()}"

    con = duckdb.connect(db, read_only=True)
    cols = {c[1] for c in con.execute("PRAGMA table_info('mlb_sgp_odds')").fetchall()}
    con.close()
    assert "spread_line" in cols
    assert "total_line" in cols


def test_fd_shim_writes_new_schema_cols(tmp_path):
    """After shim refactor, scraper writes spread_line/total_line cols.

    Mirror of test_dk_shim_writes_new_schema_cols for FanDuel. Shim-shape
    check only — verifies the scraper:
      1) exits 0 when there are no target lines to scrape,
      2) creates the mlb_sgp_odds table with the post-Task-4 schema
         (spread_line + total_line columns present).

    Full pricing correctness is covered by the live-API golden tests
    above and by Task 28's manual smoke test.
    """
    import os
    import subprocess

    db = str(tmp_path / "shim_test.duckdb")
    # Empty parlay_lines table so the shim has no targets to scrape and
    # never hits FD's API. The shim should still ensure_table + exit 0.
    con = duckdb.connect(db)
    con.execute("""
        CREATE TABLE mlb_parlay_lines (
            game_id VARCHAR, home_team VARCHAR, away_team VARCHAR,
            commence_time VARCHAR, fg_spread DOUBLE, fg_total DOUBLE,
            f5_spread DOUBLE, f5_total DOUBLE
        )
    """)
    con.close()

    env = os.environ.copy()
    env["MLB_SGP_DB_PATH"] = db
    repo_root = (
        "/Users/callancapitolo/NFLWork/.claude/worktrees/"
        "kalshi-mlb-rfq-line-source-pivot"
    )
    result = subprocess.run(
        [sys.executable, "mlb_sgp/scraper_fanduel_sgp.py"],
        env=env, capture_output=True, timeout=60, cwd=repo_root,
    )
    assert result.returncode == 0, f"stderr: {result.stderr.decode()}"

    con = duckdb.connect(db, read_only=True)
    cols = {c[1] for c in con.execute("PRAGMA table_info('mlb_sgp_odds')").fetchall()}
    con.close()
    assert "spread_line" in cols
    assert "total_line" in cols


def test_px_shim_writes_new_schema_cols(tmp_path):
    """After shim refactor, scraper writes spread_line/total_line cols.

    Mirror of test_dk_shim_writes_new_schema_cols for ProphetX. Shim-shape
    check only — verifies the scraper:
      1) exits 0 when there are no target lines to scrape,
      2) creates the mlb_sgp_odds table with the post-Task-4 schema
         (spread_line + total_line columns present).

    Full pricing correctness is covered by the live-API golden tests
    above and by Task 28's manual smoke test.
    """
    import os
    import subprocess

    db = str(tmp_path / "shim_test.duckdb")
    # Empty parlay_lines table so the shim has no targets to scrape and
    # never hits PX's API. The shim should still ensure_table + exit 0.
    con = duckdb.connect(db)
    con.execute("""
        CREATE TABLE mlb_parlay_lines (
            game_id VARCHAR, home_team VARCHAR, away_team VARCHAR,
            commence_time VARCHAR, fg_spread DOUBLE, fg_total DOUBLE,
            f5_spread DOUBLE, f5_total DOUBLE
        )
    """)
    con.close()

    env = os.environ.copy()
    env["MLB_SGP_DB_PATH"] = db
    repo_root = (
        "/Users/callancapitolo/NFLWork/.claude/worktrees/"
        "kalshi-mlb-rfq-line-source-pivot"
    )
    result = subprocess.run(
        [sys.executable, "mlb_sgp/scraper_prophetx_sgp.py"],
        env=env, capture_output=True, timeout=60, cwd=repo_root,
    )
    assert result.returncode == 0, f"stderr: {result.stderr.decode()}"

    con = duckdb.connect(db, read_only=True)
    cols = {c[1] for c in con.execute("PRAGMA table_info('mlb_sgp_odds')").fetchall()}
    con.close()
    assert "spread_line" in cols
    assert "total_line" in cols


def test_nv_shim_writes_new_schema_cols(tmp_path):
    """After shim refactor, scraper writes spread_line/total_line cols.

    Mirror of test_dk_shim_writes_new_schema_cols for Novig. Shim-shape
    check only — verifies the scraper:
      1) exits 0 when there are no target lines to scrape,
      2) creates the mlb_sgp_odds table with the post-Task-4 schema
         (spread_line + total_line columns present).

    Full pricing correctness is covered by the live-API golden tests
    above and by Task 28's manual smoke test.
    """
    import os
    import subprocess

    db = str(tmp_path / "shim_test.duckdb")
    # Empty parlay_lines table so the shim has no targets to scrape and
    # never hits Novig's API. The shim should still ensure_table + exit 0.
    con = duckdb.connect(db)
    con.execute("""
        CREATE TABLE mlb_parlay_lines (
            game_id VARCHAR, home_team VARCHAR, away_team VARCHAR,
            commence_time VARCHAR, fg_spread DOUBLE, fg_total DOUBLE,
            f5_spread DOUBLE, f5_total DOUBLE
        )
    """)
    con.close()

    env = os.environ.copy()
    env["MLB_SGP_DB_PATH"] = db
    repo_root = (
        "/Users/callancapitolo/NFLWork/.claude/worktrees/"
        "kalshi-mlb-rfq-line-source-pivot"
    )
    result = subprocess.run(
        [sys.executable, "mlb_sgp/scraper_novig_sgp.py"],
        env=env, capture_output=True, timeout=60, cwd=repo_root,
    )
    assert result.returncode == 0, f"stderr: {result.stderr.decode()}"

    con = duckdb.connect(db, read_only=True)
    cols = {c[1] for c in con.execute("PRAGMA table_info('mlb_sgp_odds')").fetchall()}
    con.close()
    assert "spread_line" in cols
    assert "total_line" in cols
