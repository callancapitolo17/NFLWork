"""Tests for mlb_sgp/_shared.py — shared dataclasses + helpers."""
import dataclasses
from datetime import datetime, timezone

import pytest

from mlb_sgp._shared import TargetLine, PricedRow


def test_target_line_immutable():
    t = TargetLine(
        game_id="abc-123",
        home_team="New York Yankees",
        away_team="Boston Red Sox",
        commence_time=datetime(2026, 5, 13, 23, 0, tzinfo=timezone.utc),
        period="FG",
        spread=-1.5,
        total=8.5,
    )
    with pytest.raises(dataclasses.FrozenInstanceError):
        t.spread = -2.5  # type: ignore[misc]


def test_priced_row_fields():
    p = PricedRow(
        game_id="abc-123",
        combo="Home Spread + Over",
        period="FG",
        spread_line=-1.5,
        total_line=8.5,
        bookmaker="draftkings",
        source="draftkings_direct",
        sgp_decimal=2.85,
        sgp_american=185,
        fetch_time=datetime.now(timezone.utc),
    )
    assert p.combo == "Home Spread + Over"
    assert p.spread_line == -1.5
    assert p.sgp_american == 185


def test_priced_row_immutable():
    p = PricedRow(
        game_id="abc-123",
        combo="Home Spread + Over",
        period="FG",
        spread_line=-1.5,
        total_line=8.5,
        bookmaker="draftkings",
        source="draftkings_direct",
        sgp_decimal=2.85,
        sgp_american=185,
        fetch_time=datetime.now(timezone.utc),
    )
    with pytest.raises(dataclasses.FrozenInstanceError):
        p.sgp_decimal = 3.0  # type: ignore[misc]
