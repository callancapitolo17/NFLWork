"""Issue #86: F5-winner tie conversion math (fail-first).

Kalshi's KXMLBF5 team markets are UNCONDITIONAL — "If the score at the end
of the first 5 innings is tied, the 'Tie' strike will resolve to Yes and
all team strikes will resolve to No" (live rules_secondary, 2026-08-12).
Books' F5 ML is push-refund 2-way (recon: DK/FD/MGM/NV implied sums
1.01-1.07), i.e. P(team wins F5 | not tied) — a CONDITIONAL number.

Because "team leads after 5" is a subset of "not tied", the conversion is
an exact identity, valid for any combo containing an F5-winner leg:

    P(combo) = P(combo | not tied) x (1 - P(tie))

with the MARGINAL P(tie) — no correlation assumption, even when the combo
also carries F5 totals.
"""
import pytest

from kalshi_common import fair_value


# --------------------------------------------------------------------- #
# f5_winner_fair_from_book: conditional -> unconditional                 #
# --------------------------------------------------------------------- #

def test_push_2way_even_odds_hand_case():
    # Issue #86 spec hand case: even-odds conditional 0.50, P(tie)=0.13
    out = fair_value.f5_winner_fair_from_book(0.5, 0.13, "push_2way")
    assert out == pytest.approx(0.435)


def test_push_2way_live_fd_case():
    # Live recon 2026-08-12 KC@LAD: FD 2-way devig LAD 0.6236, FD 3-way
    # devigged tie 0.1423 -> 0.5349 unconditional (FD's own direct 3-way
    # LAD devig was 0.5364 — the identity closes to 0.15pt).
    out = fair_value.f5_winner_fair_from_book(0.6236, 0.1423, "push_2way")
    assert out == pytest.approx(0.6236 * (1.0 - 0.1423))


def test_three_way_unconditional_is_identity():
    # A book that priced the team leg from a true 3-way market already
    # quotes unconditional space — no adjustment, p_tie unused.
    out = fair_value.f5_winner_fair_from_book(0.5364, 0.1423,
                                              "three_way_unconditional")
    assert out == pytest.approx(0.5364)
    # ... even when no tie anchor exists at all.
    out = fair_value.f5_winner_fair_from_book(0.5364, None,
                                              "three_way_unconditional")
    assert out == pytest.approx(0.5364)


def test_unknown_semantics_fails_closed():
    assert fair_value.f5_winner_fair_from_book(0.5, 0.13, "draw_no_bet") is None
    assert fair_value.f5_winner_fair_from_book(0.5, 0.13, None) is None
    assert fair_value.f5_winner_fair_from_book(0.5, 0.13, "") is None


def test_push_2way_without_tie_anchor_fails_closed():
    assert fair_value.f5_winner_fair_from_book(0.5, None, "push_2way") is None


def test_push_2way_insane_tie_prob_fails_closed():
    # 0 / negative / >= the 0.5 sanity ceiling (a plausible F5 tie mass is
    # ~0.05-0.30; 0.5+ means swapped outcomes in a parse, not baseball).
    for bad in (0.0, -0.1, 0.5, 0.7, 1.0, 1.5):
        assert fair_value.f5_winner_fair_from_book(0.5, bad, "push_2way") is None


def test_invalid_conditional_fair_fails_closed():
    for bad in (None, 0.0, 1.0, -0.2, 1.7, "x"):
        assert fair_value.f5_winner_fair_from_book(bad, 0.13, "push_2way") is None


# --------------------------------------------------------------------- #
# tie_prob_from_three_way: devigged P(tie) off a 3-way result market     #
# --------------------------------------------------------------------- #

def test_symmetric_three_way_is_one_third():
    # Equal decimals must devig to exactly 1/3 each (probit preserves
    # symmetry), independent of the vig level.
    out = fair_value.tie_prob_from_three_way(3.0, 3.0, 3.0)
    assert out == pytest.approx(1.0 / 3.0)
    out = fair_value.tie_prob_from_three_way(2.7, 2.7, 2.7)
    assert out == pytest.approx(1.0 / 3.0)


def test_fd_live_case_lands_in_plausible_band():
    # FD "First 5 Innings Result" KC@LAD 2026-08-12: KC 2.88 / Tie 6.5 /
    # LAD 1.7246 (raw sum 1.0809; proportional devig gives tie 0.1423).
    out = fair_value.tie_prob_from_three_way(1.72463768115942, 6.5, 2.88)
    assert out is not None
    assert 0.10 < out < 0.17


def test_higher_tie_decimal_means_lower_tie_prob():
    lo = fair_value.tie_prob_from_three_way(1.72, 8.0, 2.88)
    hi = fair_value.tie_prob_from_three_way(1.72, 5.0, 2.88)
    assert lo is not None and hi is not None
    assert lo < hi


def test_overround_envelope_fails_closed():
    # Sum of implieds must sit in a sane vig envelope: a crossed book
    # (sum < 1) or absurd overround is a feed bug, never devigged.
    assert fair_value.tie_prob_from_three_way(5.0, 9.0, 5.0) is None   # ~0.51
    assert fair_value.tie_prob_from_three_way(1.5, 1.5, 1.5) is None   # 2.0


def test_invalid_decimals_fail_closed():
    assert fair_value.tie_prob_from_three_way(None, 6.5, 2.88) is None
    assert fair_value.tie_prob_from_three_way(1.72, 0.0, 2.88) is None
    assert fair_value.tie_prob_from_three_way(1.72, 6.5, 1.0) is None
    assert fair_value.tie_prob_from_three_way("x", 6.5, 2.88) is None
