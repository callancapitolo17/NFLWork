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
    not devig_prob - median_fair. Non-Kalshi flag requires the book to be
    *below* the median (bettable YES edge)."""
    _insert_odds("m1", "draftkings", 100, 0.500, 0.420, datetime.now())
    _insert_odds("m1", "fanduel",    100, 0.500, 0.420, datetime.now())
    _insert_odds("m1", "bookmaker", +233, 0.300, 0.305, datetime.now())
    rows = q.ev_candidates(threshold_pp=10)
    assert len(rows) == 1
    r = rows[0]
    assert r["book"] == "bookmaker"
    # delta = implied (0.300) - median_fair (~0.420) = ~-0.120
    assert r["delta"] == pytest.approx(0.300 - 0.420, abs=0.01)
    assert r["book_prob"] == pytest.approx(0.300)


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
