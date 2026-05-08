"""
Manual integration test for the integer-line derivation feature.

Run this after a slate refresh to verify the new fallback path triggers
on real data. Not part of the regular pytest suite — run explicitly:

    /Users/callancapitolo/NFLWork/mlb_sgp/venv/bin/python -m mlb_sgp.tests.test_integration_integer_line

Exits 0 if a slate has integer-line games AND _interpolated rows were written.
Exits 1 if integer-line games exist but no _interpolated rows appeared
(indicates the fallback isn't triggering or all bounds checks are failing).
Exits 2 if no integer-line games on slate (test inconclusive).
"""
import sys
from pathlib import Path

import duckdb

# Resolve repo root from this file's location: <repo>/mlb_sgp/tests/test_integration_integer_line.py
_REPO_ROOT = Path(__file__).resolve().parents[2]
DB_PATH = str(_REPO_ROOT / "Answer Keys" / "mlb_mm.duckdb")


def main() -> int:
    con = duckdb.connect(DB_PATH, read_only=True)

    # Are there any integer-line games on the current slate?
    integer_games = con.execute("""
        SELECT game_id, fg_total, f5_total
        FROM mlb_parlay_lines
        WHERE fg_total = ROUND(fg_total) OR f5_total = ROUND(f5_total)
    """).fetchall()

    if not integer_games:
        print("No integer-line games on current slate — test inconclusive.")
        return 2

    print(f"Found {len(integer_games)} integer-line game(s) on slate.")

    # Did any book write _interpolated rows for those games?
    interpolated = con.execute("""
        SELECT bookmaker, source, COUNT(*) AS n
        FROM mlb_sgp_odds
        WHERE source LIKE '%_interpolated'
        GROUP BY bookmaker, source
        ORDER BY bookmaker
    """).fetchall()

    if not interpolated:
        print("FAIL: integer-line games exist but no _interpolated rows in mlb_sgp_odds.")
        print("Investigation: check scraper logs for bounds-check WARN lines.")
        return 1

    print("PASS: _interpolated rows present.")
    for row in interpolated:
        print(f"  {row[0]:12s} {row[1]:30s} n={row[2]}")

    # Sum-to-one invariant check
    invariants = con.execute("""
        SELECT bookmaker, game_id, period,
               COUNT(*) AS n_combos,
               ROUND(SUM(1.0/sgp_decimal), 4) AS sum_implied
        FROM mlb_sgp_odds
        WHERE source LIKE '%_interpolated'
        GROUP BY bookmaker, game_id, period
        HAVING n_combos != 4 OR sum_implied < 0.97 OR sum_implied > 1.03
    """).fetchall()

    if invariants:
        print(f"FAIL: {len(invariants)} interpolated game/period(s) violate sum-to-one or n_combos invariant:")
        for row in invariants:
            print(f"  {row}")
        return 1

    print("Sum-to-one invariant holds for all interpolated game/periods.")
    return 0


if __name__ == "__main__":
    sys.exit(main())
