"""Issue #23: Kalshi constituent-singles anchor — pure math.

Every leg of a combo IS its own 2-way Kalshi market, trading in real time
while our book scrapes refresh every ~150-165s. These functions turn those
books into (a) devigged marginals for the quote-time sanity gate and (b) a
jump signal for the risk sweep.

Every function must fail SOFT on a degenerate Kalshi book rather than raise:
`yes_ask=100` (the divide-by-zero gotcha), an empty 0-100 book, a crossed or
locked book, or a ticker missing from a snapshot.
"""
from kalshi_mlb_mm import singles

SPREAD = "KXMLBSPREAD-25JUN271905TEXLAA-LAA2"
TOTAL = "KXMLBTOTAL-25JUN271905TEXLAA-9"
LEGS = [{"market_ticker": SPREAD, "event_ticker": "E", "side": "yes"},
        {"market_ticker": TOTAL, "event_ticker": "E", "side": "yes"}]
SNAP = {SPREAD: {"yes_bid": 0.45, "yes_ask": 0.47},
        TOTAL: {"yes_bid": 0.52, "yes_ask": 0.54}}


# --------------------------------------------------------------------------- #
# devigged_yes                                                                #
# --------------------------------------------------------------------------- #

def test_devigged_yes_removes_the_spread_overround():
    """Buying YES costs the ask (0.47); buying NO costs 1 - bid (0.55). Those
    imply 1.02 of probability — the 2c spread. Devigging lands the fair inside
    the quoted band."""
    p = singles.devigged_yes(0.45, 0.47)
    assert p is not None
    assert 0.45 < p < 0.47


def test_devigged_yes_is_symmetric_about_a_balanced_book():
    assert abs(singles.devigged_yes(0.49, 0.51) - 0.50) < 1e-9


def test_devigged_yes_returns_none_on_degenerate_shapes():
    assert singles.devigged_yes(0.45, 1.00) is None    # yes_ask=100 gotcha
    assert singles.devigged_yes(0.0, 1.00) is None     # empty 0-100 book
    assert singles.devigged_yes(0.0, 0.03) is None     # no_ask = 1.0
    assert singles.devigged_yes(0.60, 0.40) is None    # crossed
    assert singles.devigged_yes(0.50, 0.50) is None    # locked, no spread
    assert singles.devigged_yes(None, 0.47) is None
    assert singles.devigged_yes("x", 0.47) is None
    assert singles.devigged_yes(float("nan"), 0.47) is None
    assert singles.devigged_yes(0.45, float("inf")) is None


# --------------------------------------------------------------------------- #
# marginals_for_legs                                                          #
# --------------------------------------------------------------------------- #

def test_marginals_follow_the_side_each_leg_takes():
    yes_side = singles.marginals_for_legs(SNAP, LEGS)
    no_legs = [dict(LEGS[0], side="no"), LEGS[1]]
    no_side = singles.marginals_for_legs(SNAP, no_legs)
    assert abs(yes_side[0] + no_side[0] - 1.0) < 1e-9
    assert abs(yes_side[1] - no_side[1]) < 1e-12


def test_marginals_fail_closed_on_missing_or_degenerate_leg():
    assert singles.marginals_for_legs({SPREAD: SNAP[SPREAD]}, LEGS) is None
    bad = {SPREAD: SNAP[SPREAD], TOTAL: {"yes_bid": 0.0, "yes_ask": 1.0}}
    assert singles.marginals_for_legs(bad, LEGS) is None
    assert singles.marginals_for_legs({}, LEGS) is None
    assert singles.marginals_for_legs(SNAP, []) is None


# --------------------------------------------------------------------------- #
# corr_sanity                                                                 #
# --------------------------------------------------------------------------- #

def test_corr_sanity_accepts_a_plausible_positive_correlation():
    """p1=p2=0.50 -> independent baseline 0.25, Frechet [0.0, 0.50]. A
    same-game fair of 0.30 is a 1.2x premium: legitimate, must pass."""
    s = singles.corr_sanity(0.30, [0.50, 0.50], 0.5, 2.0)
    assert s.ok is True and s.reason is None
    assert abs(s.baseline_independent - 0.25) < 1e-12
    assert abs(s.premium - 1.2) < 1e-12
    assert (s.frechet_lo, s.frechet_hi) == (0.0, 0.50)


def test_corr_sanity_rejects_above_the_frechet_upper_bound():
    """A joint can never exceed its smallest marginal. 0.55 > min(0.50, 0.60)."""
    s = singles.corr_sanity(0.55, [0.50, 0.60], 0.5, 2.0)
    assert s.ok is False and s.reason == "frechet"


def test_corr_sanity_rejects_below_the_frechet_lower_bound():
    """p1 + p2 - 1 = 0.50 is the floor when both legs are likely."""
    s = singles.corr_sanity(0.40, [0.80, 0.70], 0.5, 2.0)
    assert s.ok is False and s.reason == "frechet"


def test_corr_sanity_flags_premium_band_only_inside_frechet():
    """0.13 sits inside Frechet [0, 0.5] but is a 0.52x premium on a 0.25
    baseline -> premium violation, not a Frechet one."""
    s = singles.corr_sanity(0.13, [0.50, 0.50], 0.6, 2.0)
    assert s.ok is False and s.reason == "premium"


def test_corr_sanity_checks_frechet_before_premium():
    """Both violated -> the arithmetic verdict wins; premium is heuristic."""
    s = singles.corr_sanity(0.55, [0.50, 0.60], 0.99, 1.01)
    assert s.reason == "frechet"


def test_corr_sanity_handles_three_legs():
    """Frechet generalises: lo = max(0, sum - (n-1))."""
    s = singles.corr_sanity(0.10, [0.50, 0.50, 0.50], 0.5, 2.0)
    assert abs(s.baseline_independent - 0.125) < 1e-12
    assert s.frechet_lo == 0.0 and s.frechet_hi == 0.50
    assert s.ok is True


def test_corr_sanity_returns_none_without_a_usable_anchor():
    assert singles.corr_sanity(0.3, [], 0.5, 2.0) is None
    assert singles.corr_sanity(0.3, None, 0.5, 2.0) is None
    assert singles.corr_sanity(None, [0.5, 0.5], 0.5, 2.0) is None
    assert singles.corr_sanity(float("nan"), [0.5, 0.5], 0.5, 2.0) is None
    assert singles.corr_sanity(0.0, [0.5, 0.5], 0.5, 2.0) is None
    assert singles.corr_sanity("x", [0.5, 0.5], 0.5, 2.0) is None


# --------------------------------------------------------------------------- #
# jumped_tickers                                                              #
# --------------------------------------------------------------------------- #

def test_jumped_tickers_reports_only_moves_past_the_threshold():
    current = {SPREAD: {"yes_bid": 0.55, "yes_ask": 0.57},    # ~ +0.10
               TOTAL: {"yes_bid": 0.525, "yes_ask": 0.545}}   # ~ +0.005
    moved = singles.jumped_tickers(SNAP, current, 0.03)
    assert set(moved) == {SPREAD}
    assert moved[SPREAD] > 0.03


def test_jumped_tickers_is_side_agnostic():
    """|delta P(YES)| == |delta P(NO)|, so a downward move counts too."""
    current = {SPREAD: {"yes_bid": 0.35, "yes_ask": 0.37}, TOTAL: SNAP[TOTAL]}
    assert set(singles.jumped_tickers(SNAP, current, 0.03)) == {SPREAD}


def test_jumped_tickers_skips_unreadable_or_degenerate_sides():
    """Absence of signal is NOT a jump — an API blip must never flush the
    resting book (#17's confirm veto is the fail-closed backstop)."""
    assert singles.jumped_tickers(SNAP, {}, 0.03) == {}
    assert singles.jumped_tickers({}, SNAP, 0.03) == {}
    degenerate = {SPREAD: {"yes_bid": 0.0, "yes_ask": 1.0}, TOTAL: SNAP[TOTAL]}
    assert singles.jumped_tickers(SNAP, degenerate, 0.03) == {}
    partial = {TOTAL: SNAP[TOTAL]}                 # SPREAD missing entirely
    assert singles.jumped_tickers(SNAP, partial, 0.03) == {}


def test_jumped_tickers_ignores_a_move_exactly_at_the_threshold():
    """Strictly greater-than, so a threshold of 0 does not fire on noise-free
    equality."""
    assert singles.jumped_tickers(SNAP, dict(SNAP), 0.0) == {}


# --------------------------------------------------------------------------- #
# fetch_market_prices — the only API-touching function                        #
# --------------------------------------------------------------------------- #

def test_fetch_market_prices_is_one_get_per_distinct_ticker(monkeypatch):
    from kalshi_common import auth_client
    calls = []

    def fake_api(method, path, *a, **k):
        calls.append(path)
        return 200, {"market": {"yes_bid": 45, "yes_ask": 47}}, None

    monkeypatch.setattr(auth_client, "api", fake_api)
    out = singles.fetch_market_prices({SPREAD, TOTAL})
    assert len(calls) == 2
    assert set(out) == {SPREAD, TOTAL}
    assert out[SPREAD] == {"yes_bid": 0.45, "yes_ask": 0.47}


def test_fetch_market_prices_returns_partial_on_failure(monkeypatch):
    """Best-effort by design: unreadable tickers are absent, not fatal."""
    from kalshi_common import auth_client

    def fake_api(method, path, *a, **k):
        if SPREAD in path:
            return 200, {"market": {"yes_bid": 45, "yes_ask": 47}}, None
        raise RuntimeError("boom")

    monkeypatch.setattr(auth_client, "api", fake_api)
    out = singles.fetch_market_prices({SPREAD, TOTAL})
    assert set(out) == {SPREAD}


def test_fetch_market_prices_makes_no_calls_for_an_empty_set(monkeypatch):
    from kalshi_common import auth_client
    monkeypatch.setattr(auth_client, "api",
                        lambda *a, **k: (_ for _ in ()).throw(
                            AssertionError("must not call the API")))
    assert singles.fetch_market_prices(set()) == {}
