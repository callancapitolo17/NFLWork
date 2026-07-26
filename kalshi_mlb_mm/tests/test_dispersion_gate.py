"""Issue #20: consensus gate = z-space dispersion threshold (no outlier band).

Gate: quote only if n_books >= MIN_AGREEING_BOOKS AND stddev of the book
fairs in z-space (norm.ppf, ddof=1) <= SIGMA_Z_MAX. Fair = median of ALL
books — no outlier removal, no second median. Tested against BOTH
implementations (main._consensus_filter oracle, router.consensus production),
which must agree on shared cases.

The tail tests are the ticket's must-fail-first proof: under the old absolute
±2¢ band a book 2¢ away at p≈0.08 (25% relative disagreement) was tolerated.
"""
import pytest

from kalshi_mlb_mm import main, config, router


# ---------------------------------------------------------------------------
# Tail behavior (must FAIL on the old ±2¢ band code)
# ---------------------------------------------------------------------------

def test_tail_three_books_one_2c_away_declined(monkeypatch):
    """p≈0.08: 8¢/8¢/6¢ — a 2¢ gap is a 25% relative fight at this price.
    Old band kept all three; the dispersion gate declines the combo."""
    monkeypatch.setattr(config, "MIN_AGREEING_BOOKS", 2)
    out = main._consensus_filter({"fd": 0.08, "mgm": 0.08, "px": 0.06})
    assert out == {}, f"longshot dispersion must decline the quote, got {out}"


def test_tail_two_books_2c_apart_declined(monkeypatch):
    """p≈0.08 with 2 books 2¢ apart: old band passed both (each 1¢ from
    median); z-space dispersion declines."""
    monkeypatch.setattr(config, "MIN_AGREEING_BOOKS", 2)
    out = main._consensus_filter({"fd": 0.08, "mgm": 0.06})
    assert out == {}


# ---------------------------------------------------------------------------
# Continuity at mid prices (same accept/decline as the old band)
# ---------------------------------------------------------------------------

def test_mid_two_books_3c_gap_kept(monkeypatch):
    """p≈0.50, 3¢ gap: old band kept both (each 1.5¢ from median); new gate
    keeps both too (sigma_z ≈ 0.053 <= 0.07)."""
    monkeypatch.setattr(config, "MIN_AGREEING_BOOKS", 2)
    out = main._consensus_filter({"fd": 0.50, "mgm": 0.47})
    assert set(out) == {"fd", "mgm"}


def test_mid_two_books_5c_gap_declined(monkeypatch):
    """p≈0.50, 5¢ gap: old band declined (each 2.5¢ from median); new gate
    declines too (sigma_z ≈ 0.089 > 0.07)."""
    monkeypatch.setattr(config, "MIN_AGREEING_BOOKS", 2)
    out = main._consensus_filter({"fd": 0.50, "mgm": 0.45})
    assert out == {}


def test_mid_three_books_tight_cluster_kept(monkeypatch):
    monkeypatch.setattr(config, "MIN_AGREEING_BOOKS", 2)
    out = main._consensus_filter({"fd": 0.42, "mgm": 0.43, "px": 0.44})
    assert len(out) == 3


# ---------------------------------------------------------------------------
# Degenerate cases
# ---------------------------------------------------------------------------

def test_single_book_declined_by_count(monkeypatch):
    """One book has sigma_z = 0 by construction — the count check must refuse
    it BEFORE dispersion is even considered."""
    monkeypatch.setattr(config, "MIN_AGREEING_BOOKS", 2)
    assert main._consensus_filter({"fd": 0.43}) == {}


def test_two_identical_books_kept(monkeypatch):
    monkeypatch.setattr(config, "MIN_AGREEING_BOOKS", 2)
    out = main._consensus_filter({"fd": 0.43, "mgm": 0.43})
    assert set(out) == {"fd", "mgm"}


def test_no_outlier_rescue(monkeypatch):
    """DELIBERATE behavior change vs the old band: one loud dissenter now
    kills the whole combo's quote (it may be the informed book) instead of
    being outvoted by the pack."""
    monkeypatch.setattr(config, "MIN_AGREEING_BOOKS", 2)
    out = main._consensus_filter({"dk": 0.42, "fd": 0.43, "px": 0.43, "nv": 0.55})
    assert out == {}


# ---------------------------------------------------------------------------
# Production path: router.consensus — same gate, explicit reasons
# ---------------------------------------------------------------------------

def test_router_consensus_ok_returns_median_of_all():
    cons, reason = router.consensus({"fd": 0.42, "mgm": 0.43, "px": 0.44},
                                    min_books=2, sigma_z_max=0.07)
    assert reason == "ok"
    assert cons.fair == pytest.approx(0.43)
    assert cons.n_books == 3
    assert cons.sigma_pts == pytest.approx(0.01, abs=1e-9)
    assert cons.sigma_z > 0


def test_router_consensus_too_few_books():
    cons, reason = router.consensus({"fd": 0.43}, min_books=2, sigma_z_max=0.07)
    assert cons is None and reason == "too_few_books"


def test_router_consensus_dispersion_reason():
    cons, reason = router.consensus({"fd": 0.08, "mgm": 0.06},
                                    min_books=2, sigma_z_max=0.07)
    assert cons is None and reason == "consensus_dispersion"


@pytest.mark.parametrize("book_fairs", [
    {"fd": 0.08, "mgm": 0.08, "px": 0.06},
    {"fd": 0.08, "mgm": 0.06},
    {"fd": 0.50, "mgm": 0.47},
    {"fd": 0.50, "mgm": 0.45},
    {"fd": 0.42, "mgm": 0.43, "px": 0.44},
    {"fd": 0.43},
    {"dk": 0.42, "fd": 0.43, "px": 0.43, "nv": 0.55},
])
def test_oracle_and_production_agree(book_fairs, monkeypatch):
    """main._consensus_filter (test oracle) and router.consensus (production)
    must make identical quote/no-quote decisions on the same inputs."""
    monkeypatch.setattr(config, "MIN_AGREEING_BOOKS", 2)
    oracle_quotes = main._consensus_filter(book_fairs) != {}
    cons, _ = router.consensus(book_fairs, min_books=2,
                               sigma_z_max=config.SIGMA_Z_MAX)
    assert oracle_quotes == (cons is not None)
    if cons is not None:
        import statistics
        assert cons.fair == pytest.approx(statistics.median(book_fairs.values()))


def test_consensus_detail_lockstep():
    """consensus_detail returns (fair, ALL books) when the gate passes."""
    det = router.consensus_detail({"fd": 0.42, "mgm": 0.43, "px": 0.44},
                                  min_books=2, sigma_z_max=0.07)
    assert det is not None
    fair, books = det
    assert fair == pytest.approx(0.43)
    assert books == ["fd", "mgm", "px"]
    assert router.consensus_detail({"fd": 0.08, "mgm": 0.06},
                                   min_books=2, sigma_z_max=0.07) is None
