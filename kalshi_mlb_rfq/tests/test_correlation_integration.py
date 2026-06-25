"""Integration test: correlated same-game second bet is down-sized.

Strategy
--------
We exercise Tasks 2–5 end-to-end in memory:
  - Seed `_SGP_ODDS_CACHE` with a complete 4-cell grid (all four Home/Away ×
    Over/Under combos) for game "g1" at spread -1.5 / total 8.5 from two
    books.  That gives `_grid_lookup` enough rows to call `devig_book` (the
    probit path requires >=4 cells per book) and return a median probability.
  - Call `_book_implied_cov(game_id, new_region, new_price, pos_region,
    pos_price)` directly to compute the return-space covariance between a new
    "Home Spread + Over" bet and an already-held identical position.
  - Assert that the covariance is strictly positive (same-direction combos are
    positively correlated).
  - Assert that `kelly.kelly_size_combo(…, existing_positions=[held])` returns
    strictly fewer contracts than the baseline (no held position).

Why not wire `_load_existing_positions_book` against a live DuckDB?
  That function reads `positions` + `combo_cache`, parses legs_json, and
  calls `_leg_dict_to_typed` — requiring a complete schema and populated rows.
  Setting that up in a single test adds substantial fixture weight for a check
  that is already unit-tested by Tasks 2 and 3.  The core invariant (correlated
  add is down-sized) is demonstrated fully by the direct `_book_implied_cov` +
  `kelly_size_combo` path, which exercises the real book-implied covariance
  engine end-to-end.  We do include a lightweight smoke test that confirms
  `_load_existing_positions_book` returns `[]` for a game with no held
  positions, so the DB-read path is exercised without building a full fixture.
"""

import pandas as pd
import pytest
import duckdb

import kalshi_mlb_rfq.main as m
import kalshi_mlb_rfq.config as cfg
from kalshi_mlb_rfq import correlation, kelly


# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------

def _seed_grid() -> pd.DataFrame:
    """Minimal 4-cell grid per book for game 'g1' at spread -1.5, total 8.5.

    All four Home/Away × Over/Under combos are needed so that `devig_book`
    can use the probit 4-cell path rather than the single-side heuristic.
    We set symmetric decimal odds (3.3 ≈ 30.3% implied each side) so the
    devigged fair for any single quadrant lands around 0.25.
    """
    rows = []
    cells = [
        ("Home Spread + Over",  3.3),
        ("Home Spread + Under", 3.3),
        ("Away Spread + Over",  3.3),
        ("Away Spread + Under", 3.3),
    ]
    for book in ("draftkings", "fanduel"):
        for (label, dec) in cells:
            rows.append({
                "game_id":    "g1",
                "combo":      label,
                "bookmaker":  book,
                "sgp_decimal": dec,
                "spread_line": -1.5,
                "total_line":   8.5,
                "period":      "FG",
            })
    return pd.DataFrame(rows)


# ---------------------------------------------------------------------------
# Tests
# ---------------------------------------------------------------------------

class TestCorrelatedAddIsDownsized:
    """Book-implied covariance engine correctly penalises correlated positions."""

    def test_book_implied_cov_is_positive_for_same_direction_combo(self, monkeypatch):
        """Two identical Home+Over regions must have strictly positive cov_return."""
        monkeypatch.setattr(cfg, "USE_MODEL", False)
        monkeypatch.setattr(cfg, "MIN_BOOK_COUNT_FOR_BLEND", 2)
        monkeypatch.setattr(m, "_SGP_ODDS_CACHE", _seed_grid())

        region = correlation.ComboRegion("home", -1.5, "over", 8.5)
        price = 0.30

        cov = m._book_implied_cov("g1", region, price, region, price)
        assert cov > 0, (
            f"Same-direction combos must have positive cov_return, got {cov}"
        )

    def test_correlated_second_add_is_strictly_smaller(self, monkeypatch):
        """Kelly contract count for a correlated new bet must be < independent baseline.

        This is the core invariant: after holding a strongly-correlated position
        in the same game-combo, the sizing engine must reduce new exposure.
        """
        monkeypatch.setattr(cfg, "USE_MODEL", False)
        monkeypatch.setattr(cfg, "MIN_BOOK_COUNT_FOR_BLEND", 2)
        monkeypatch.setattr(m, "_SGP_ODDS_CACHE", _seed_grid())

        region = correlation.ComboRegion("home", -1.5, "over", 8.5)
        price = 0.30
        blended_fair = 0.45   # strongly +EV at price 0.30
        bankroll = 1_000.0
        kf = 0.25

        # Compute real cov_return from the book-implied grid.
        cov = m._book_implied_cov("g1", region, price, region, price)

        # Baseline: no existing position.
        base = kelly.kelly_size_combo(
            outcome_vec=None,
            blended_fair=blended_fair,
            existing_positions=[],
            effective_price=price,
            bankroll=bankroll,
            kelly_fraction=kf,
        )

        # With a 40-contract held position in the same correlated combo.
        held_position = [{"cov_return": cov, "contracts": 40, "effective_price": price}]
        corr = kelly.kelly_size_combo(
            outcome_vec=None,
            blended_fair=blended_fair,
            existing_positions=held_position,
            effective_price=price,
            bankroll=bankroll,
            kelly_fraction=kf,
        )

        assert base > 0, "Baseline Kelly must be positive (sanity check on +EV parameters)"
        assert corr < base, (
            f"Correlated second add should be smaller than baseline: "
            f"corr={corr} vs base={base}, cov_return={cov:.4f}"
        )

    def test_independent_position_does_not_downsize(self, monkeypatch):
        """A held position with zero covariance (cov_return=0) must not reduce sizing."""
        monkeypatch.setattr(cfg, "USE_MODEL", False)
        monkeypatch.setattr(m, "_SGP_ODDS_CACHE", _seed_grid())

        price = 0.30
        blended_fair = 0.45
        bankroll = 1_000.0
        kf = 0.25

        base = kelly.kelly_size_combo(
            outcome_vec=None,
            blended_fair=blended_fair,
            existing_positions=[],
            effective_price=price,
            bankroll=bankroll,
            kelly_fraction=kf,
        )

        indep_position = [{"cov_return": 0.0, "contracts": 40, "effective_price": price}]
        indep = kelly.kelly_size_combo(
            outcome_vec=None,
            blended_fair=blended_fair,
            existing_positions=indep_position,
            effective_price=price,
            bankroll=bankroll,
            kelly_fraction=kf,
        )

        assert indep == base, (
            f"Independent position (cov=0) must not change sizing: "
            f"indep={indep} vs base={base}"
        )

    def test_load_existing_positions_book_empty_when_no_held(self, monkeypatch, tmp_path):
        """Smoke: `_load_existing_positions_book` returns [] when positions table is empty."""
        monkeypatch.setattr(cfg, "USE_MODEL", False)
        monkeypatch.setattr(m, "_SGP_ODDS_CACHE", _seed_grid())

        # Build a minimal DuckDB with the required tables but no held positions.
        bot_db = str(tmp_path / "bot.duckdb")
        con = duckdb.connect(bot_db)
        con.execute("""
            CREATE TABLE positions (
                combo_market_ticker VARCHAR,
                game_id VARCHAR,
                net_contracts DOUBLE,
                weighted_price DOUBLE,
                side VARCHAR
            )
        """)
        con.execute("""
            CREATE TABLE combo_cache (
                combo_market_ticker VARCHAR,
                legs_json VARCHAR,
                game_id VARCHAR
            )
        """)
        con.close()

        # `db.connect` reads `DB_PATH` from `kalshi_mlb_rfq.config` via `db.py`.
        monkeypatch.setattr(m.db, "DB_PATH", tmp_path / "bot.duckdb")

        region = correlation.ComboRegion("home", -1.5, "over", 8.5)
        result = m._load_existing_positions_book("g1", region, 0.30, "yes")
        assert result == [], f"Expected [] for empty positions table, got {result}"
