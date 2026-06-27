import datetime
from unabated_edge import mapping, feed
from unabated_edge.sports.soccer import Soccer


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
