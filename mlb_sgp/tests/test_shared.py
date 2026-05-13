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


import duckdb
import tempfile
from pathlib import Path


def test_load_target_lines_prefers_mlb_target_lines(tmp_path):
    """When both tables exist, mlb_target_lines wins."""
    from mlb_sgp._shared import load_target_lines
    db = str(tmp_path / "t.duckdb")
    con = duckdb.connect(db)
    con.execute("""
        CREATE TABLE mlb_target_lines (
            game_id VARCHAR, home_team VARCHAR, away_team VARCHAR,
            commence_time TIMESTAMP, period VARCHAR,
            spread DOUBLE, total DOUBLE, written_at TIMESTAMP
        )
    """)
    con.execute("""
        INSERT INTO mlb_target_lines VALUES
            ('g1', 'NYY', 'BOS', '2026-05-13 23:00', 'FG', -1.5, 8.5, NOW()),
            ('g1', 'NYY', 'BOS', '2026-05-13 23:00', 'FG', -2.5, 8.5, NOW()),
            ('g1', 'NYY', 'BOS', '2026-05-13 23:00', 'FG', -1.5, 9.5, NOW())
    """)
    con.execute("""
        CREATE TABLE mlb_parlay_lines (
            game_id VARCHAR, home_team VARCHAR, away_team VARCHAR,
            commence_time VARCHAR, fg_spread DOUBLE, fg_total DOUBLE,
            f5_spread DOUBLE, f5_total DOUBLE
        )
    """)
    con.execute("INSERT INTO mlb_parlay_lines VALUES ('g1','NYY','BOS','2026-05-13T23:00Z',-3.5,7.5,NULL,NULL)")
    con.close()

    rows = load_target_lines(db_path=db)
    assert len(rows) == 3, "multi-row from mlb_target_lines"
    spreads = sorted({r.spread for r in rows})
    assert -3.5 not in spreads, "should NOT pull from mlb_parlay_lines"
    assert all(r.period == "FG" for r in rows)


def test_load_target_lines_legacy_fallback(tmp_path):
    """When only mlb_parlay_lines exists, emit FG + F5 rows from each game."""
    from mlb_sgp._shared import load_target_lines
    db = str(tmp_path / "t.duckdb")
    con = duckdb.connect(db)
    con.execute("""
        CREATE TABLE mlb_parlay_lines (
            game_id VARCHAR, home_team VARCHAR, away_team VARCHAR,
            commence_time VARCHAR, fg_spread DOUBLE, fg_total DOUBLE,
            f5_spread DOUBLE, f5_total DOUBLE
        )
    """)
    con.execute("""
        INSERT INTO mlb_parlay_lines VALUES
            ('g1', 'NYY', 'BOS', '2026-05-13T23:00:00+00:00', -1.5, 8.5, -0.5, 4.5),
            ('g2', 'LAD', 'SF',  '2026-05-14T02:00:00+00:00', -2.5, 7.5, NULL, NULL)
    """)
    con.close()

    rows = load_target_lines(db_path=db)
    fg = [r for r in rows if r.period == "FG"]
    f5 = [r for r in rows if r.period == "F5"]
    assert len(fg) == 2, "two FG rows (one per game)"
    assert len(f5) == 1, "one F5 row (g2 had NULL F5)"
    assert {r.game_id for r in fg} == {"g1", "g2"}
    assert {r.game_id for r in f5} == {"g1"}


def test_load_target_lines_empty_db(tmp_path):
    from mlb_sgp._shared import load_target_lines
    db = str(tmp_path / "empty.duckdb")
    con = duckdb.connect(db)
    con.close()
    assert load_target_lines(db_path=db) == []
