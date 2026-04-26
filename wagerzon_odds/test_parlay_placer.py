# wagerzon_odds/test_parlay_placer.py
"""Unit tests for parlay_placer. All Wagerzon HTTP calls are mocked."""
import pytest
from unittest.mock import patch, MagicMock
from parlay_placer import Leg, ParlaySpec, encode_sel, encode_detail_data


@pytest.fixture(autouse=True)
def reset_session():
    """Clear cached session before/after each test to avoid leakage."""
    import parlay_placer
    parlay_placer._clear_session()
    yield
    parlay_placer._clear_session()


def test_leg_dataclass_basic():
    leg = Leg(idgm=5632938, play=0, points=-1.5, odds=117, pitcher=0)
    assert leg.idgm == 5632938
    assert leg.play == 0


def test_parlay_spec_dataclass_basic():
    spec = ParlaySpec(
        parlay_hash="abc123",
        legs=[
            Leg(idgm=5632938, play=1, points=-1.5, odds=117),
            Leg(idgm=5632938, play=2, points=-7.5, odds=105),
        ],
        amount=15.0,
        expected_win=12.50,
        expected_risk=15.0,
    )
    assert len(spec.legs) == 2


def test_encode_sel_two_legs():
    legs = [
        Leg(idgm=5632938, play=1, points=-1.5, odds=117),
        Leg(idgm=5632938, play=2, points=-7.5, odds=105),
    ]
    assert encode_sel(legs) == "1_5632938_-1.5_117,2_5632938_-7.5_105"


def test_encode_sel_negative_odds():
    legs = [Leg(idgm=5632938, play=5, points=0, odds=-140)]
    assert encode_sel(legs) == "5_5632938_0_-140"


def test_encode_detail_data_shape():
    legs = [Leg(idgm=5632938, play=1, points=-1.5, odds=117)]
    out = encode_detail_data(legs, amount=15.0)
    # detail_data is a list of dicts ready for json.dumps
    assert out[0]["IdGame"] == 5632938
    assert out[0]["Play"] == 1
    assert out[0]["Amount"] == "15"
    assert out[0]["Points"]["selected"] is True


def test_encode_detail_data_non_integer_amount():
    """Kelly recommendations may be fractional dollars (e.g. $15.75)."""
    legs = [Leg(idgm=5632938, play=1, points=-1.5, odds=117)]
    out = encode_detail_data(legs, amount=15.75)
    assert out[0]["Amount"] == "15.75"


def test_get_session_caches_within_call(monkeypatch):
    """Two calls to _get_session in the same place_parlays() invocation
    should reuse one session, not log in twice."""
    import parlay_placer
    monkeypatch.setenv("WAGERZON_USERNAME", "user")
    monkeypatch.setenv("WAGERZON_PASSWORD", "pass")

    with patch("parlay_placer.requests.Session") as MockSession:
        mock_session = MagicMock()
        # Simulate already-authenticated GET (URL contains NewSchedule)
        response = MagicMock(url="https://backend.wagerzon.com/wager/NewSchedule.aspx")
        mock_session.get.return_value = response
        MockSession.return_value = mock_session

        s1 = parlay_placer._get_session()
        s2 = parlay_placer._get_session()
        assert s1 is s2
        # Only one Session() instantiation
        assert MockSession.call_count == 1
