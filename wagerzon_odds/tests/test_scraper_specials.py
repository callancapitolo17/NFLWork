"""Tests for wagerzon_odds.scraper_specials"""
import json
from pathlib import Path

# pytest discovers this via conftest path; for now just import directly
import sys
sys.path.insert(0, str(Path(__file__).parent.parent))
from scraper_specials import parse_specials_json, extract_team_from_htm, extract_prop_type


# Fixture: today's recon JSON saved to disk during reconnaissance
RECON_FULL = Path(__file__).parent.parent / "recon_specials_full.json"


def test_extract_team_triple_play():
    assert extract_team_from_htm("DODGERS TRIPLE-PLAY (SCR 1ST, 1H & GM)") == "DODGERS"

def test_extract_team_grand_slam_multiword():
    assert extract_team_from_htm("WHITE SOX GRAND-SLAM (SCR 1ST, 1H, GM & SCR U4½)") == "WHITE SOX"

def test_extract_team_no_prop_type_returns_none():
    assert extract_team_from_htm("MARLINS, ANGELS & YANKEES ALL WIN") is None

def test_extract_team_game_level_g_s_returns_none():
    # CHC-SDG G-S is game-level (cross-team) — no single subject team
    assert extract_team_from_htm("CHC-SDG G-S (GM O7½, 1H U4½, HITS U15½, 1ST HR 2R)") is None

def test_extract_team_grd_slm_alias():
    # Wagerzon abbreviated grand-slam variant — same leg grammar as GRAND-SLAM
    assert extract_team_from_htm("ASTROS GRD-SLM (SCR 1ST, 1H, GM & SCR U3½)") == "ASTROS"

def test_extract_team_trple_play_typo():
    # Wagerzon misspells TRIPLE-PLAY as TRPLE-PLAY in some props
    assert extract_team_from_htm("ORIOLES TRPLE-PLAY (SCR 1ST, 1H & GM)") == "ORIOLES"

def test_extract_prop_type_canonicalizes_grd_slm():
    assert extract_prop_type("YANKEES GRD-SLM (SCR 1ST, 1H, GM & SCR O4½)") == "GRAND-SLAM"

def test_extract_prop_type_canonicalizes_trple_play():
    assert extract_prop_type("DODGERS TRPLE-PLAY (SCR 1ST, 1H & GM)") == "TRIPLE-PLAY"

def test_extract_prop_type_returns_canonical_for_canonical_input():
    assert extract_prop_type("DODGERS TRIPLE-PLAY (SCR 1ST, 1H & GM)") == "TRIPLE-PLAY"
    assert extract_prop_type("RANGERS GRAND-SLAM (SCR 1ST, 1H, GM & SCR U3½)") == "GRAND-SLAM"

def test_parse_specials_json_canonicalizes_alias():
    """GRD-SLM rows should land in the DB as prop_type='GRAND-SLAM' so the
    pricer's `WHERE prop_type IN ('TRIPLE-PLAY','GRAND-SLAM')` picks them up."""
    raw = {"result": {"listLeagues": [[{"Description": "MLB - SPECIALS", "Games": [
        {"htm": "ASTROS GRD-SLM (SCR 1ST, 1H, GM & SCR U3½)",
         "hnum": 757042, "idgm": 5635935,
         "gmdt": "20260427", "gmtm": "19:40:00",
         "GameLines": [{"odds": "365", "oddsh": "+365"}]},
    ]}]]}}
    rows = parse_specials_json(raw, sport="mlb", league_id=4899)
    assert len(rows) == 1
    assert rows[0]["prop_type"] == "GRAND-SLAM"  # canonicalized, not raw GRD-SLM
    assert rows[0]["team"] == "ASTROS"
    assert rows[0]["odds"] == 365


def test_parse_specials_json_today_fixture():
    """Run against today's saved JSON. Should yield the 14 priceable rows."""
    if not RECON_FULL.exists():
        import pytest
        pytest.skip(f"{RECON_FULL} not present; run scraper_specials.py once first")
    raw = json.loads(RECON_FULL.read_text())
    rows = parse_specials_json(raw, sport="mlb", league_id=4899)
    # Expect at least 14 priceable props (8 TP + 6 single-team GS)
    priceable = [r for r in rows if r["prop_type"] in ("TRIPLE-PLAY", "GRAND-SLAM") and r["team"]]
    assert len(priceable) >= 14, f"Got {len(priceable)} priceable, expected >=14"
    # Required fields present
    for r in priceable:
        assert r["scraped_at"] is not None
        assert r["sport"] == "mlb"
        assert r["league_id"] == 4899
        assert isinstance(r["odds"], int)
        assert r["description"]
        assert r["prop_type"] in ("TRIPLE-PLAY", "GRAND-SLAM")
        assert r["team"]


def test_odds_parsed_as_signed_int():
    raw = {"result": {"listLeagues": [[{"Description": "MLB - SPECIALS", "Games": [
        {"htm": "DODGERS TRIPLE-PLAY (SCR 1ST, 1H & GM)",
         "hnum": 757041, "idgm": 5635934,
         "gmdt": "20260427", "gmtm": "19:40:00",
         "GameLines": [{"odds": "110", "oddsh": "-110"}]},
    ]}]]}}
    rows = parse_specials_json(raw, sport="mlb", league_id=4899)
    dodgers = [r for r in rows if r["team"] == "DODGERS"][0]
    assert dodgers["odds"] == -110  # signed


def test_no_specials_section_returns_empty():
    raw = {"result": {"listLeagues": [[{"Description": "MLB - REGULAR", "Games": [
        {"htm": "Yankees", "hnum": 1, "idgm": 1,
         "gmdt": "20260427", "gmtm": "19:40:00",
         "GameLines": [{"odds": "100", "oddsh": "+100"}]},
    ]}]]}}
    rows = parse_specials_json(raw, sport="mlb", league_id=4899)
    assert rows == []
