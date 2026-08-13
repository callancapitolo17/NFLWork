"""Shim-shape regression tests for the SGP scraper CLI entry points.

Each shim must, given a target table with no rows to scrape:
  1) exit 0 without touching the book's API, and
  2) create ``mlb_sgp_odds`` with the post-Task-4 schema (spread_line +
     total_line present).

Pricing correctness lives in ``test_golden_baselines.py`` (frozen payloads
replayed offline) and in the per-book fixture tests.
"""
import sys
import duckdb
import subprocess
from pathlib import Path


# The two live-API "golden" tests that used to live here (DK and FD) were
# REMOVED by issue #42, not merely regenerated.
#
# They shelled out to the real scraper against the live book, wrote into the
# SHARED production DuckDB (Answer Keys/mlb_mm.duckdb) from a test run, and
# then diffed against a CSV baseline keyed by game_id. Every game_id in that
# baseline expires when the game starts, so the tests failed with "no keys in
# common" on any normal day and regenerating the CSV would have bought exactly
# one more day. They also could not run offline.
#
# Their replacement is mlb_sgp/tests/test_golden_baselines.py: payloads are
# captured once into a frozen fixture and replayed offline, pinning the same
# thing these pinned (the orchestrator's PricedRow output) for all six books
# instead of two. See mlb_sgp/tests/golden_replay.py.

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
    # Derive the repo root from this file. It used to be a hardcoded
    # absolute path into a worktree that has long since been deleted,
    # so every one of these four tests died with FileNotFoundError
    # regardless of the code under test (issue #42).
    repo_root = str(Path(__file__).resolve().parents[2])
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
    # Derive the repo root from this file. It used to be a hardcoded
    # absolute path into a worktree that has long since been deleted,
    # so every one of these four tests died with FileNotFoundError
    # regardless of the code under test (issue #42).
    repo_root = str(Path(__file__).resolve().parents[2])
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
    # Derive the repo root from this file. It used to be a hardcoded
    # absolute path into a worktree that has long since been deleted,
    # so every one of these four tests died with FileNotFoundError
    # regardless of the code under test (issue #42).
    repo_root = str(Path(__file__).resolve().parents[2])
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
    # Derive the repo root from this file. It used to be a hardcoded
    # absolute path into a worktree that has long since been deleted,
    # so every one of these four tests died with FileNotFoundError
    # regardless of the code under test (issue #42).
    repo_root = str(Path(__file__).resolve().parents[2])
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
