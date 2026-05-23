"""Unit tests for validate_idgmtyp_shape().

The validator is warn-only — it doesn't drop records, just surfaces
that the response shape has drifted. These tests pin down:
  - Known idgmtyps with valid shape pass without warning.
  - Known idgmtyps with all-null probe fields warn (drift signal).
  - Unknown idgmtyps warn — but only ONCE per scrape (dedupe).
  - Empty GameLines is handled defensively (no crash).
"""
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))
import scraper_v2  # noqa: E402


def _child(idgmtyp, line_fields: dict, vtm: str = "TEST CHILD"):
    return {
        "idgmtyp": idgmtyp,
        "vtm": vtm,
        "GameLines": [line_fields],
    }


# ---- Known idgmtyps pass when probe fields populated ----------------------

def test_idgmtyp_10_full_line_passes_with_spread():
    child = _child(10, {"vsprdt": "-1.5"})
    assert scraper_v2.validate_idgmtyp_shape(child) is True


def test_idgmtyp_10_full_line_passes_with_total():
    child = _child(10, {"ovt": "8.5"})
    assert scraper_v2.validate_idgmtyp_shape(child) is True


def test_idgmtyp_19_hits_total_passes_with_ovt():
    child = _child(19, {"ovt": "13.5"}, vtm="GUARDIANS TOTAL HITS")
    assert scraper_v2.validate_idgmtyp_shape(child) is True


def test_idgmtyp_31_oddeven_passes_with_voddst():
    child = _child(31, {"voddst": "-110"}, vtm="TOTAL RUNS ODD")
    assert scraper_v2.validate_idgmtyp_shape(child) is True


def test_idgmtyp_35_team_total_passes_with_ovt():
    child = _child(35, {"ovt": "4.5"}, vtm="GUARDIANS TEAM TOTAL")
    assert scraper_v2.validate_idgmtyp_shape(child) is True


# ---- Drift signal: all probe fields null --------------------------------

def test_idgmtyp_35_warns_when_all_probes_null(capsys):
    child = _child(35, {"vsprdt": "-1.5"}, vtm="GUARDIANS TEAM TOTAL")
    # ovt is None → drift warning.
    assert scraper_v2.validate_idgmtyp_shape(child) is False
    log = capsys.readouterr().out
    assert "WARNING" in log
    assert "idgmtyp=35" in log
    assert "shape may have drifted" in log


def test_idgmtyp_31_warns_when_no_ml_fields(capsys):
    child = _child(31, {"ovt": "8.5"}, vtm="TOTAL RUNS ODD")
    # ML probe fields all None → drift warning.
    assert scraper_v2.validate_idgmtyp_shape(child) is False
    log = capsys.readouterr().out
    assert "idgmtyp=31" in log


# ---- Unknown idgmtyp + dedupe ---------------------------------------------

def test_unknown_idgmtyp_warns_once_and_returns_false(capsys):
    warned = set()
    c1 = _child(999, {"ovt": "5.5"}, vtm="MYSTERY MARKET A")
    c2 = _child(999, {"ovt": "6.5"}, vtm="MYSTERY MARKET B")
    assert scraper_v2.validate_idgmtyp_shape(c1, warned) is False
    assert scraper_v2.validate_idgmtyp_shape(c2, warned) is False
    log = capsys.readouterr().out
    # Should warn for the FIRST occurrence and stay silent for the second.
    assert log.count("unknown WZ idgmtyp=999") == 1
    assert 999 in warned


def test_different_unknown_idgmtyps_each_warn(capsys):
    warned = set()
    scraper_v2.validate_idgmtyp_shape(_child(998, {"ovt": "5.5"}), warned)
    scraper_v2.validate_idgmtyp_shape(_child(999, {"ovt": "5.5"}), warned)
    log = capsys.readouterr().out
    assert "idgmtyp=998" in log
    assert "idgmtyp=999" in log


def test_unknown_idgmtyp_without_dedupe_set_still_warns(capsys):
    """Backward-compat: caller may omit `warned_unknown` (e.g., one-shot
    validation in a test). Should still print the warning."""
    assert scraper_v2.validate_idgmtyp_shape(
        _child(7777, {"ovt": "5.5"})
    ) is False
    log = capsys.readouterr().out
    assert "idgmtyp=7777" in log


# ---- Defensive paths -------------------------------------------------------

def test_empty_gamelines_returns_false_no_crash():
    child = {"idgmtyp": 10, "vtm": "EMPTY", "GameLines": []}
    assert scraper_v2.validate_idgmtyp_shape(child) is False


def test_missing_gamelines_key_returns_false_no_crash():
    child = {"idgmtyp": 10, "vtm": "NO LINES KEY"}
    assert scraper_v2.validate_idgmtyp_shape(child) is False
