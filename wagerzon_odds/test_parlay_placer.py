# wagerzon_odds/test_parlay_placer.py
"""Unit tests for parlay_placer. All Wagerzon HTTP calls are mocked."""
import pytest
from parlay_placer import Leg, ParlaySpec, encode_sel, encode_detail_data


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
