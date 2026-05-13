"""Tests for mlb_sgp/_shared.py — shared dataclasses + helpers."""
import dataclasses
from datetime import datetime, timezone

import pytest

from mlb_sgp._shared import (
    TargetLine,
    PricedRow,
    decimal_to_american,
    american_to_decimal,
    _utc_bucket,
)


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


def test_decimal_to_american_favorite():
    # Standard formula: dec >= 2.0 → +((dec-1)*100); dec < 2.0 → -100/(dec-1)
    assert decimal_to_american(1.5) == -200  # -100 / 0.5 = -200
    assert decimal_to_american(2.0) == 100   # +(1.0 * 100)
    assert decimal_to_american(2.5) == 150   # +(1.5 * 100)
    assert decimal_to_american(3.0) == 200


def test_american_to_decimal_round_trip():
    for am in [-300, -150, 100, 150, 300]:
        dec = american_to_decimal(am)
        back = decimal_to_american(dec)
        assert back == am, f"round trip {am} → {dec} → {back}"


def test_utc_bucket_isolates_hour():
    from datetime import datetime, timezone
    t = datetime(2026, 5, 13, 23, 17, 42, tzinfo=timezone.utc)
    assert _utc_bucket(t) == "2026-05-13T23"


def test_utc_bucket_handles_iso_string():
    assert _utc_bucket("2026-05-13T23:17:42Z") == "2026-05-13T23"
    assert _utc_bucket("2026-05-13T23:17:42+00:00") == "2026-05-13T23"


def test_utc_bucket_handles_empty():
    # Existing scrapers pass `lines.get("commence_time", "")` and expect ""
    # back so the team-only fallback matcher fires. Match that contract.
    assert _utc_bucket("") == ""
    assert _utc_bucket(None) == ""


def test_utc_bucket_naive_datetime_assumed_utc():
    t = datetime(2026, 5, 13, 23, 17, 42)  # no tzinfo
    assert _utc_bucket(t) == "2026-05-13T23"


def test_utc_bucket_converts_non_utc_tz():
    from datetime import timedelta
    eastern = timezone(timedelta(hours=-4))  # EDT-ish
    # 23:00 EDT = 03:00 UTC next day
    t = datetime(2026, 5, 13, 23, 0, tzinfo=eastern)
    assert _utc_bucket(t) == "2026-05-14T03"
