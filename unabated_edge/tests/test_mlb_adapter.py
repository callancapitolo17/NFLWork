from unabated_edge.sports.mlb import Mlb
from unabated_edge.sports import registry
from unabated_edge.tests.conftest import build_bt3_state


def test_canon_team_city_and_nickname_collapse():
    m = Mlb()
    assert m.canon_team("New York Yankees") == "Yankees"
    assert m.canon_team("Yankees") == "Yankees"
    assert m.canon_team("St. Louis Cardinals") == "Cardinals"
    assert m.canon_team("Athletics") == "Athletics"


def test_event_teams_parses_kalshi_ticker():
    """Kalshi MLB titles are city-based with short disambiguator suffixes
    (Task 1 recon: "New York Y vs Philadelphia: Total Runs"), not
    nickname-based — too fragile to parse as free text. event_teams instead
    reads the two club codes out of event_ticker, e.g.
    "KXMLBTOTAL-26JUL251805NYYPHI" (NYY=Yankees, PHI=Phillies)."""
    m = Mlb()
    ev = {"event_ticker": "KXMLBTOTAL-26JUL251805NYYPHI"}
    assert m.event_teams(ev) == frozenset({"Yankees", "Phillies"})


def test_event_teams_parses_mixed_length_club_codes():
    """LAA (3-letter) + SF (2-letter) exercises the 3+2 split, not just the
    3+3 case above — Task 1 recon: "Los Angeles A vs San Francisco" (Angels
    @ Giants)."""
    m = Mlb()
    ev = {"event_ticker": "KXMLBTOTAL-26JUL251840LAASF"}
    assert m.event_teams(ev) == frozenset({"Angels", "Giants"})


def test_event_teams_parses_az_diamondbacks_code():
    """Regression: live pairing smoke against today's slate found Kalshi
    uses "AZ" for Arizona, not "ARI" as first guessed from convention alone
    — a real code-map gap, not a doubleheader/date artifact."""
    m = Mlb()
    ev = {"event_ticker": "KXMLBTOTAL-26JUL261335AZWSH"}
    assert m.event_teams(ev) == frozenset({"Diamondbacks", "Nationals"})


def test_event_teams_malformed_ticker_fails_closed():
    m = Mlb()
    # Ticker shape matches but codes aren't in the known club-code table.
    assert m.event_teams({"event_ticker": "KXMLBTOTAL-26JUL251805XXXYYY"}) == frozenset()
    # Doesn't match the ticker shape at all.
    assert m.event_teams({"event_ticker": "GARBAGE"}) == frozenset()
    # No event_ticker key.
    assert m.event_teams({}) == frozenset()


def test_mlb_registered():
    assert any(a.sport == "mlb" for a in registry.ADAPTERS)
    assert len({a.league_prefix for a in registry.ADAPTERS}) == len(registry.ADAPTERS)


def test_price_event_prices_matching_rung_skips_others():
    """Inherited TotalsLadderAdapter path with an MLB-typical line: anchor
    8.5 both sides -> yes+no candidates at the 8.5 rung; a Kalshi market at
    a different floor_strike (no anchor line there) is skipped."""
    m = Mlb()
    st = build_bt3_state("lg5", "Philadelphia Phillies", "New York Yankees",
                          eid=108549, ms=7, over_price=-119, under_price=104, line=8.5)
    ev_meta = st.events[108549]
    kalshi_event = {
        "event_ticker": "KXMLBTOTAL-26JUL251805NYYPHI",
        "title": "New York Y vs Philadelphia: Total Runs",
        "markets": [
            {"ticker": "KXMLBTOTAL-26JUL251805NYYPHI-8", "strike_type": "greater", "floor_strike": 8.5},
            {"ticker": "KXMLBTOTAL-26JUL251805NYYPHI-9", "strike_type": "greater", "floor_strike": 9.5},
        ],
    }
    candidates = m.price_event(st, ev_meta, kalshi_event)
    assert len(candidates) == 2  # only the 8.5 rung; 9.5 has no anchor line -> skipped
    assert {c.market_ticker for c in candidates} == {"KXMLBTOTAL-26JUL251805NYYPHI-8"}
    sides = {c.side for c in candidates}
    assert sides == {"yes", "no"}
    labels = {c.label for c in candidates}
    assert labels == {"over_8.5", "under_8.5"}
    total_prob = sum(c.fair_prob for c in candidates)
    assert abs(total_prob - 1.0) < 1e-5
