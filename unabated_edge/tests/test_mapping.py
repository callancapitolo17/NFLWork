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


def _total_kev(home, away, ticker="T-O25"):
    """Minimal KXWCTOTAL-style Kalshi event with one Over market."""
    return {
        "title": f"{home} vs {away}: Regulation Time Total Goals",
        "markets": [
            {"ticker": ticker, "strike_type": "greater", "floor_strike": 2.5,
             "yes_sub_title": "Reg Time: Over 2.5 goals scored"},
        ]
    }


def test_pair_events_matches_by_team_pair():
    """pair_events returns (event_meta, kalshi_event) tuple when team pairs agree."""
    ev = _ev("Argentina", "Austria", event_id=1)
    kev = _total_kev("Argentina", "Austria")
    result = mapping.pair_events(Soccer(), [ev], [kev])
    assert len(result) == 1
    em, matched_kev = result[0]
    assert em.home == "Argentina"
    assert matched_kev is kev


def test_pair_events_matches_across_divergent_spellings():
    """Aliased Unabated names must pair with plain-English Kalshi title names."""
    ev = _ev("Côte d'Ivoire", "Czechia")
    kev = _total_kev("Ivory Coast", "Czech Republic")
    result = mapping.pair_events(Soccer(), [ev], [kev])
    assert len(result) == 1


def test_pair_events_no_match_returns_empty():
    ev = _ev("Nowhere United", "Parts Unknown", event_id=99)
    result = mapping.pair_events(Soccer(), [ev], [])
    assert result == []


def test_pair_events_multiple_events():
    """Multiple events matched independently."""
    ev1 = _ev("Argentina", "Austria", event_id=1)
    ev2 = _ev("Colombia", "Ghana", event_id=2)
    kev1 = _total_kev("Argentina", "Austria", "T-ARG")
    kev2 = _total_kev("Colombia", "Ghana", "T-COL")
    result = mapping.pair_events(Soccer(), [ev1, ev2], [kev1, kev2])
    assert len(result) == 2


def test_pair_events_logs_unmatched_once(caplog):
    """An unmatched event is logged once, not re-logged on the next tick (dedup)."""
    mapping._warned_unmatched.clear()
    ev = _ev("Nowhere United", "Parts Unknown", event_id=777)
    kalshi_events = []
    with caplog.at_level(logging.INFO, logger="unabated_edge"):
        mapping.pair_events(Soccer(), [ev], kalshi_events)
        mapping.pair_events(Soccer(), [ev], kalshi_events)   # second tick, same event
    unmatched_logs = [r for r in caplog.records if "unmatched" in r.getMessage()]
    assert len(unmatched_logs) == 1, f"expected 1 dedup'd log, got {len(unmatched_logs)}"
    assert "Nowhere United" in unmatched_logs[0].getMessage()
    mapping._warned_unmatched.clear()
