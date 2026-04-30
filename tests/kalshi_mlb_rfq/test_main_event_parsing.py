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
