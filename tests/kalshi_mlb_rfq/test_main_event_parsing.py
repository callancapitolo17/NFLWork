from kalshi_mlb_rfq import main


def test_parse_3letter_both():
    """DET @ ATL — both 3-letter codes."""
    assert main._parse_event_suffix("26APR291915DETATL") == ("DET", "ATL")


def test_parse_2letter_away_3letter_home():
    """KC @ ATH (post-rebrand) — 2-letter away, 3-letter home.

    Regression: the old fixed-width slice grabbed the trailing '0' from '40'
    as part of away_code, producing '0KC'. With ATH missing from the dict it
    also dropped the home → game silently skipped.
    """
    assert main._parse_event_suffix("26APR292140KCATH") == ("KC", "ATH")


def test_parse_2letter_away_3letter_home_az():
    """AZ @ MIL — same regression class as KCATH (rebranded code + 2-letter)."""
    assert main._parse_event_suffix("26APR291940AZMIL") == ("AZ", "MIL")


def test_parse_3letter_away_2letter_home():
    """CWS @ SD — disambiguation matters: 'CWSSD' must NOT parse as ('CW','SSD')."""
    assert main._parse_event_suffix("26MAY012140CWSSD") == ("CWS", "SD")


def test_parse_2letter_both():
    """SF @ TB — both 2-letter codes."""
    assert main._parse_event_suffix("26MAY021810SFTB") == ("SF", "TB")


def test_parse_unknown_codes_returns_none():
    """Unknown codes should return (None, None) so the caller drops the event."""
    assert main._parse_event_suffix("26APR292140XXXYYY") == (None, None)


def test_parse_too_short_returns_none():
    """Suffix shorter than date+team-block fails cleanly, no IndexError."""
    assert main._parse_event_suffix("26APR29") == (None, None)
    assert main._parse_event_suffix("") == (None, None)


def test_home_code_from_event_ticker_2letter_home():
    """SPREAD event_ticker with 2-letter home (TB) — old code returned 'BSF', should be 'TB'."""
    assert main._home_code_from_event_ticker("KXMLBSPREAD-26MAY021810SFTB") == "TB"


def test_home_code_from_event_ticker_3letter_home():
    assert main._home_code_from_event_ticker("KXMLBSPREAD-26APR292140KCATH") == "ATH"


def test_parse_doubleheader_game_1():
    """Kalshi appends G1/G2 to both games of a doubleheader.

    Live 2026-09-01: KXMLBGAME-26SEP041410DETCLEG1. The old fixed grammar
    read 'EG1' / 'G1' as the home code and rejected the event entirely.
    """
    assert main._parse_event_suffix("26SEP041410DETCLEG1") == ("DET", "CLE")


def test_parse_doubleheader_game_2():
    assert main._parse_event_suffix("26SEP041915DETCLEG2") == ("DET", "CLE")


def test_parse_doubleheader_with_2letter_home():
    # KC @ SD game 2 — the marker is stripped before the 3-then-2 home probe,
    # so it cannot steal characters from a short home code.
    assert main._parse_event_suffix("26SEP042140KCSDG2") == ("KC", "SD")


def test_home_code_from_doubleheader_event_ticker():
    """Regression: an unparseable home code made every spread leg type as the
    AWAY team, so CLE -3.5 was stored as +3.5 — a wrong line, not a decline."""
    assert main._home_code_from_event_ticker(
        "KXMLBSPREAD-26SEP041410DETCLEG1") == "CLE"
