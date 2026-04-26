import sys
from pathlib import Path
sys.path.insert(0, str(Path(__file__).resolve().parents[3] / "wagerzon_odds"))

import pytest
from unittest.mock import MagicMock
import parlay_pricer


def test_get_combined_parlay_price_builds_sel_with_per_leg_idgm():
    """Cross-game parlay must encode each leg's own idgm in the sel string."""
    legs = [
        {"idgm": 100001, "play": 1, "points": "-1.5", "odds": 110},   # game A home spread
        {"idgm": 100001, "play": 2, "points": "9.5", "odds": -105},   # game A over
        {"idgm": 100002, "play": 0, "points": "+1.5", "odds": -120},  # game B away spread
        {"idgm": 100002, "play": 3, "points": "7.5", "odds": -110},   # game B under
    ]
    expected_sel_parts = [
        "1_100001_-1.5_110",
        "2_100001_9.5_-105",
        "0_100002_+1.5_-120",
        "3_100002_7.5_-110",
    ]

    fake_session = MagicMock()
    fake_response = MagicMock()
    fake_response.json.return_value = {
        "result": {
            "details": [{"Win": 21000, "Risk": 1000}],
        }
    }
    fake_response.raise_for_status.return_value = None
    fake_session.post.return_value = fake_response

    result = parlay_pricer.get_combined_parlay_price(fake_session, legs, amount=1000)

    call_args = fake_session.post.call_args
    sel_value = call_args.kwargs["data"]["sel"]
    for part in expected_sel_parts:
        assert part in sel_value, f"Missing leg in sel: {part}"
    assert result is not None
    assert result["win"] == 21000
    assert result["decimal"] == pytest.approx(22.0, rel=0.01)


def test_get_combined_parlay_price_returns_none_on_error():
    legs = [
        {"idgm": 100001, "play": 1, "points": "-1.5", "odds": 110},
        {"idgm": 100002, "play": 0, "points": "+1.5", "odds": -120},
    ]
    fake_session = MagicMock()
    fake_response = MagicMock()
    fake_response.json.return_value = {
        "result": {
            "details": [{"Win": 0, "Risk": 0}],
            "ErrorMsgKey": "MAXPARLAYRISKEXCEED",
        }
    }
    fake_response.raise_for_status.return_value = None
    fake_session.post.return_value = fake_response

    result = parlay_pricer.get_combined_parlay_price(fake_session, legs, amount=10000)
    assert result is None
