"""Tests for the consensus gate oracle (main._consensus_filter).

Issue #20: z-space dispersion gate. The filter takes a dict of
{book: devigged_fair} and returns ALL books when they agree (sample stddev of
norm.ppf(fair) <= SIGMA_Z_MAX and count >= MIN_AGREEING_BOOKS), or {} when
they don't — there is no outlier removal. Deeper coverage (tail behavior,
continuity with the old ±2¢ band, oracle-vs-production agreement) lives in
test_dispersion_gate.py.
"""
from kalshi_mlb_mm import main, config


def test_consensus_three_books_agree_pass():
    """Three tightly-clustered books → gate passes, ALL kept."""
    out = main._consensus_filter({"dk": 0.42, "fd": 0.43, "px": 0.43})
    assert len(out) == 3
    assert set(out) == {"dk", "fd", "px"}


def test_consensus_one_outlier_declines_all():
    """Issue #20 behavior change: a loud dissenter (0.55 vs the 0.43 pack) is
    no longer outvoted — it may be the informed book, so the whole combo
    declines."""
    out = main._consensus_filter({"dk": 0.42, "fd": 0.43, "px": 0.43, "nv": 0.55})
    assert out == {}


def test_consensus_scattered_books_skip(monkeypatch):
    """Books telling wildly different stories (0.30 / 0.43 / 0.43 / 0.60) →
    dispersion far above SIGMA_Z_MAX → empty dict."""
    monkeypatch.setattr(config, "MIN_AGREEING_BOOKS", 3)
    out = main._consensus_filter({"dk": 0.30, "fd": 0.43, "px": 0.43, "nv": 0.60})
    assert out == {}


def test_consensus_too_few_books_skip(monkeypatch):
    """Below MIN_AGREEING_BOOKS at the count check → empty dict.
    Pinned to threshold=3: only 2 books supplied < 3 → skip before the
    dispersion check."""
    monkeypatch.setattr(config, "MIN_AGREEING_BOOKS", 3)
    out = main._consensus_filter({"dk": 0.43, "fd": 0.43})
    assert out == {}


def test_consensus_dispersion_threshold(monkeypatch):
    """A tight cluster passes; adding real disagreement declines. Pinned to
    threshold=3 so the count logic is deterministic."""
    monkeypatch.setattr(config, "MIN_AGREEING_BOOKS", 3)
    assert config.SIGMA_Z_MAX == 0.07  # sanity: default dispersion cap
    # {0.42, 0.43, 0.44}: sigma_z ~ 0.025 <= 0.07 → all kept.
    out = main._consensus_filter({"dk": 0.42, "fd": 0.43, "px": 0.44})
    assert len(out) == 3
    # {0.42, 0.43, 0.50}: sigma_z ~ 0.11 > 0.07 → decline everything.
    out = main._consensus_filter({"dk": 0.42, "fd": 0.43, "px": 0.50})
    assert out == {}, f"dispersion above cap must decline, got {out}"


def test_consensus_two_book_default(monkeypatch):
    """Production default (MIN_AGREEING_BOOKS=2): two agreeing books PASS
    (the 2026-06-18 change that unblocks quote volume)."""
    monkeypatch.setattr(config, "MIN_AGREEING_BOOKS", 2)
    # Two books in exact agreement → kept (would have been skipped at 3).
    out = main._consensus_filter({"dk": 0.43, "fd": 0.43})
    assert set(out) == {"dk", "fd"}
    # A 2-book set in real disagreement → still skipped.
    out = main._consensus_filter({"dk": 0.30, "fd": 0.43})
    assert out == {}
    # One book alone → below threshold → skip.
    out = main._consensus_filter({"dk": 0.43})
    assert out == {}
