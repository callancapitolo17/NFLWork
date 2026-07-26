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


def test_pair_events_excludes_doubleheader_both_sides_ambiguous(caplog):
    """Two Kalshi events AND two Unabated events share the same team pair
    (a doubleheader) -> that team-pair returns NO pairs at all (fail closed,
    since nothing downstream can tell which Kalshi game is which Unabated
    game). A separate, unambiguous pair in the same tick still matches."""
    mapping._warned_unmatched.clear()
    ev1 = _ev("Yankees", "Phillies", event_id=1)
    ev2 = _ev("Yankees", "Phillies", event_id=2)
    ev3 = _ev("Argentina", "Austria", event_id=3)
    kev1 = _total_kev("Yankees", "Phillies", "T-YP1")
    kev2 = _total_kev("Yankees", "Phillies", "T-YP2")
    kev3 = _total_kev("Argentina", "Austria", "T-ARG")
    with caplog.at_level(logging.WARNING, logger="unabated_edge"):
        result = mapping.pair_events(Soccer(), [ev1, ev2, ev3], [kev1, kev2, kev3])
    assert len(result) == 1
    em, kev = result[0]
    assert em.event_id == 3 and kev is kev3
    warn_logs = [r for r in caplog.records if "doubleheader" in r.getMessage()]
    assert len(warn_logs) == 1
    assert "2 kalshi events, 2 unabated events" in warn_logs[0].getMessage()
    mapping._warned_unmatched.clear()


def test_pair_events_excludes_duplicate_unabated_side_only(caplog):
    """Two Unabated events share a team pair even though only one Kalshi
    event exists for it -> still excluded (ambiguous on the Unabated side
    alone is enough; guessing which Unabated event the lone Kalshi game
    belongs to is exactly the unsafe case)."""
    mapping._warned_unmatched.clear()
    ev1 = _ev("Yankees", "Phillies", event_id=1)
    ev2 = _ev("Yankees", "Phillies", event_id=2)
    kev1 = _total_kev("Yankees", "Phillies", "T-YP1")
    with caplog.at_level(logging.WARNING, logger="unabated_edge"):
        result = mapping.pair_events(Soccer(), [ev1, ev2], [kev1])
    assert result == []
    warn_logs = [r for r in caplog.records if "doubleheader" in r.getMessage()]
    assert len(warn_logs) == 1
    assert "1 kalshi events, 2 unabated events" in warn_logs[0].getMessage()
    mapping._warned_unmatched.clear()


def test_pair_events_excludes_duplicate_kalshi_side_only(caplog):
    """Two Kalshi events share a team pair even though only one Unabated
    event exists for it (the failure mode actually observed live: a
    doubleheader's second Kalshi game used to silently vanish via
    last-write-wins and the survivor got paired to the wrong game) ->
    excluded, not silently resolved to either one."""
    mapping._warned_unmatched.clear()
    ev1 = _ev("Yankees", "Phillies", event_id=1)
    kev1 = _total_kev("Yankees", "Phillies", "T-YP1")
    kev2 = _total_kev("Yankees", "Phillies", "T-YP2")
    with caplog.at_level(logging.WARNING, logger="unabated_edge"):
        result = mapping.pair_events(Soccer(), [ev1], [kev1, kev2])
    assert result == []
    warn_logs = [r for r in caplog.records if "doubleheader" in r.getMessage()]
    assert len(warn_logs) == 1
    assert "2 kalshi events, 1 unabated events" in warn_logs[0].getMessage()
    mapping._warned_unmatched.clear()
