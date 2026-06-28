import datetime
import logging
from unabated_edge import mapping, feed
from unabated_edge.sports.soccer import Soccer


def _ev(home, away, event_id=1):
    return feed.EventMeta(
        event_id=event_id, league_key="lg21",
        start_utc=datetime.datetime(2026, 6, 28, 18, 0, tzinfo=datetime.timezone.utc),
        home_id=10, away_id=11, home=home, away=away,
    )


def test_validate_rejects_incomplete():
    a = Soccer()
    p = mapping.Pairing(None, {"home": "H", "draw": "D"})       # missing away
    assert mapping.validate(a, p, {"home": .4, "draw": .3, "away": .3}) is False


def test_validate_accepts_complete():
    a = Soccer()
    p = mapping.Pairing(None, {"home": "H", "draw": "D", "away": "A"})
    assert mapping.validate(a, p, {"home": .4, "draw": .3, "away": .3}) is True


def test_pair_events_non_aliased_teams():
    """pair_events must match non-aliased teams regardless of case from Kalshi feed."""
    # Unabated-side event: Argentina (home) vs Austria (away), title-cased
    ev = feed.EventMeta(
        event_id=1, league_key="lg21",
        start_utc=datetime.datetime(2026, 6, 28, 18, 0, tzinfo=datetime.timezone.utc),
        home_id=10, away_id=11, home="Argentina", away="Austria",
    )
    # Kalshi-side event: markets have title-cased yes_sub_title (as the API sends them)
    kalshi_event = {
        "markets": [
            {"ticker": "WC-ARG", "yes_sub_title": "Argentina"},
            {"ticker": "WC-AUT", "yes_sub_title": "Austria"},
            {"ticker": "WC-DRW", "yes_sub_title": "Draw"},
        ]
    }
    result = mapping.pair_events(Soccer(), [ev], [kalshi_event])
    assert len(result) == 1, f"Expected 1 Pairing, got {len(result)}"
    ot = result[0].outcome_tickers
    assert ot["home"] == "WC-ARG", f"home ticker wrong: {ot}"
    assert ot["away"] == "WC-AUT", f"away ticker wrong: {ot}"
    assert ot["draw"] == "WC-DRW", f"draw ticker wrong: {ot}"


def test_pair_events_matches_across_divergent_spellings():
    """Aliased Unabated names must pair with plain-English Kalshi names."""
    ev = _ev("Côte d'Ivoire", "Czechia")
    kalshi_event = {
        "markets": [
            {"ticker": "WC-CIV", "yes_sub_title": "Ivory Coast"},
            {"ticker": "WC-CZE", "yes_sub_title": "Czech Republic"},
            {"ticker": "WC-DRW", "yes_sub_title": "Draw"},
        ]
    }
    result = mapping.pair_events(Soccer(), [ev], [kalshi_event])
    assert len(result) == 1
    assert result[0].outcome_tickers["home"] == "WC-CIV"
    assert result[0].outcome_tickers["away"] == "WC-CZE"


def test_pair_events_logs_unmatched_once(caplog):
    """An unmatched event is logged once, not re-logged on the next tick (dedup)."""
    mapping._warned_unmatched.clear()
    ev = _ev("Nowhere United", "Parts Unknown", event_id=777)
    kalshi_events = []                                   # nothing to match against
    with caplog.at_level(logging.INFO, logger="unabated_edge"):
        mapping.pair_events(Soccer(), [ev], kalshi_events)
        mapping.pair_events(Soccer(), [ev], kalshi_events)   # second tick, same event
    unmatched_logs = [r for r in caplog.records if "unmatched" in r.getMessage()]
    assert len(unmatched_logs) == 1, f"expected 1 dedup'd log, got {len(unmatched_logs)}"
    assert "Nowhere United" in unmatched_logs[0].getMessage()
    mapping._warned_unmatched.clear()
