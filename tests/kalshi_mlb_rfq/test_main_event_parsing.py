from datetime import datetime, timezone

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


class TestResolveGameIdPicksTheRightGame:
    """The team pair narrows; the suffix's first pitch identifies.

    The cache's 24h horizon holds both games of a doubleheader AND the next
    game of a series, so returning the first team-pair hit returned an
    arbitrary one of them (#95: fairs off by 0.05-0.11).
    """

    CLE, DET = "Cleveland Guardians", "Detroit Tigers"

    def _cache(self, monkeypatch, rows):
        cache = {gid: {"home_team": self.CLE, "away_team": self.DET,
                       "commence_time": ct, "fg_lines": []}
                 for gid, ct in rows}
        monkeypatch.setattr(main, "_PARLAY_LINES_CACHE", cache)

    def _now(self, monkeypatch, when):
        class _DT(datetime):
            @classmethod
            def now(cls, tz=None):
                return when
        monkeypatch.setattr(main, "datetime", _DT)

    def test_each_doubleheader_game_resolves_to_its_own_id(self, monkeypatch):
        self._cache(monkeypatch, [
            ("dh-1", datetime(2026, 9, 4, 18, 10, tzinfo=timezone.utc)),
            ("dh-2", datetime(2026, 9, 4, 23, 15, tzinfo=timezone.utc)),
        ])
        self._now(monkeypatch, datetime(2026, 9, 4, 12, 0, tzinfo=timezone.utc))
        g1 = main.parse_suffix_start_utc("26SEP041410DETCLEG1")
        g2 = main.parse_suffix_start_utc("26SEP041915DETCLEG2")
        assert main._resolve_game_id("CLE", "DET", g1) == "dh-1"
        assert main._resolve_game_id("CLE", "DET", g2) == "dh-2"

    def test_indistinguishable_rows_fail_closed(self, monkeypatch):
        self._cache(monkeypatch, [
            ("a", datetime(2026, 9, 4, 18, 10, tzinfo=timezone.utc)),
            ("b", datetime(2026, 9, 4, 18, 25, tzinfo=timezone.utc)),
        ])
        self._now(monkeypatch, datetime(2026, 9, 4, 12, 0, tzinfo=timezone.utc))
        start = main.parse_suffix_start_utc("26SEP041410DETCLEG1")
        assert main._resolve_game_id("CLE", "DET", start) is None

    def test_an_unreadable_start_declines(self, monkeypatch):
        self._cache(monkeypatch, [
            ("only", datetime(2026, 9, 4, 18, 10, tzinfo=timezone.utc))])
        self._now(monkeypatch, datetime(2026, 9, 4, 12, 0, tzinfo=timezone.utc))
        assert main._resolve_game_id("CLE", "DET", None) is None

    def test_ambiguity_warns_once_until_the_cache_refreshes(self, monkeypatch,
                                                            caplog):
        """Pre-merge review fix: enumeration runs every RFQ_REFRESH_SEC, so an
        unconditional warning was ~2,880 lines/day for one stuck pair."""
        import logging
        self._cache(monkeypatch, [
            ("a", datetime(2026, 9, 4, 18, 10, tzinfo=timezone.utc)),
            ("b", datetime(2026, 9, 4, 18, 25, tzinfo=timezone.utc)),
        ])
        self._now(monkeypatch, datetime(2026, 9, 4, 12, 0, tzinfo=timezone.utc))
        main._AMBIGUOUS_WARNED.clear()
        start = main.parse_suffix_start_utc("26SEP041410DETCLEG1")
        with caplog.at_level(logging.WARNING, logger="kalshi_mlb_rfq"):
            for _ in range(5):
                assert main._resolve_game_id("CLE", "DET", start) is None
        assert sum("resolve_game_ambiguous" in r.getMessage()
                   for r in caplog.records) == 1
        # A cache refresh re-arms it.
        main._AMBIGUOUS_WARNED.clear()
        with caplog.at_level(logging.WARNING, logger="kalshi_mlb_rfq"):
            main._resolve_game_id("CLE", "DET", start)
        assert sum("resolve_game_ambiguous" in r.getMessage()
                   for r in caplog.records) == 2
