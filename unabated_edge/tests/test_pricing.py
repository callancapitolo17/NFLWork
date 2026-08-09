import pytest
from unabated_edge import pricing

def test_american_to_prob():
    assert pricing.american_to_prob(-200) == pytest.approx(2/3, abs=1e-6)
    assert pricing.american_to_prob(150) == pytest.approx(0.4, abs=1e-6)

def test_devig_three_way_sums_to_one():
    out = pricing.devig([0.50, 0.30, 0.28])
    assert sum(out) == pytest.approx(1.0, abs=1e-6)
    assert out[0] > out[1]

def test_devig_two_way():
    out = pricing.devig([0.55, 0.52])
    assert sum(out) == pytest.approx(1.0, abs=1e-6)


# ---------- overround feed-integrity gate (issue #73) ----------

_BAND = {"band_min": 1.005, "band_max": 1.20}

def test_overround_reject_crossed_pair():
    """Sum < 1 is arithmetically not a two-way market — hard reject regardless of band."""
    assert pricing.overround_reject(0.82, **_BAND) == "anchor_crossed"
    assert pricing.overround_reject(0.999999, **_BAND) == "anchor_crossed"

def test_overround_reject_out_of_band():
    assert pricing.overround_reject(1.35, **_BAND) == "anchor_overround"
    assert pricing.overround_reject(1.201, **_BAND) == "anchor_overround"
    assert pricing.overround_reject(1.002, **_BAND) == "anchor_overround"  # in [1, band_min): no real book quotes ~zero vig

def test_overround_reject_passes_sane_vig():
    assert pricing.overround_reject(1.05, **_BAND) is None
    assert pricing.overround_reject(1.005, **_BAND) is None   # band edges inclusive
    assert pricing.overround_reject(1.20, **_BAND) is None

def test_overround_reject_band_params_keyword_only():
    """Band knobs must be passed explicitly by the caller (no env/config reads
    inside the helper) — pin the keyword-only signature."""
    import inspect
    sig = inspect.signature(pricing.overround_reject)
    assert sig.parameters["band_min"].kind == inspect.Parameter.KEYWORD_ONLY
    assert sig.parameters["band_max"].kind == inspect.Parameter.KEYWORD_ONLY
