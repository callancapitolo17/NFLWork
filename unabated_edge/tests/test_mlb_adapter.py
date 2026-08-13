import datetime

from unabated_edge.sports.mlb import Mlb, _MLB_CLUB_CODES
from unabated_edge.sports import registry
from unabated_edge.tests.conftest import build_bt3_state


def test_event_date_parses_ticker():
    """Task 4b: date-aware pairing needs the game date, not just the club
    codes, out of event_ticker's date block."""
    m = Mlb()
    assert m.event_date({"event_ticker": "KXMLBTOTAL-26JUL251805NYYPHI"}) == datetime.date(2026, 7, 25)
    assert m.event_date({"event_ticker": "KXMLBTOTAL-26JUL261335CHCPIT"}) == datetime.date(2026, 7, 26)


def test_event_date_fails_closed_on_garbage():
    m = Mlb()
    assert m.event_date({"event_ticker": "GARBAGE"}) is None
    assert m.event_date({}) is None


def test_event_date_fails_closed_on_unknown_month_not_locale_dependent():
    """event_date must use an explicit month dict, not strptime's %b (which
    is locale-dependent and would silently fail-closed on every ticker in a
    non-English locale, zeroing out date-aware pairing). A ticker whose date
    block doesn't match any known 3-letter English month abbreviation must
    fail closed (None), and mixed-case month text must not match the
    ticker regex at all (also None) -- neither path should raise."""
    m = Mlb()
    # Shape matches (3 uppercase letters) but "XXX" isn't a real month.
    assert m.event_date({"event_ticker": "KXMLBTOTAL-26XXX251805NYYPHI"}) is None
    # Lowercase month text doesn't match _TICKER_RE's uppercase-only group.
    assert m.event_date({"event_ticker": "KXMLBTOTAL-26jul251805NYYPHI"}) is None


def test_import_enforces_no_2letter_code_is_prefix_of_3letter_code():
    """mlb.py asserts this invariant at module load (a 2-letter code being a
    prefix of a 3-letter one would make the ticker away/home split
    ambiguous — see the assertion's comment in sports/mlb.py). This test
    both proves the module already imported clean (the assert didn't fire
    above) and re-checks the same condition directly against the live
    table, so a future code addition that reintroduces the collision fails
    loudly here instead of silently mis-parsing tickers."""
    assert not any(
        c2 != c3 and c3.startswith(c2)
        for c2 in _MLB_CLUB_CODES for c3 in _MLB_CLUB_CODES
        if len(c2) == 2 and len(c3) == 3
    )


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


# ---------- overround feed-integrity gate (issue #73), shared via TotalsLadderAdapter ----------

def test_crossed_main_rung_rejected_for_mlb():
    """The #73 feed-integrity gate lives in the shared TotalsLadderAdapter
    base, not in Soccer -- prove it fires for the Mlb adapter too. A crossed
    pair (raw implied sum < 1, one side stale mid-refresh) must be refused
    before devig, marked rejected, never priced."""
    m = Mlb()
    st = build_bt3_state("lg5", "Philadelphia Phillies", "New York Yankees",
                          eid=108549, ms=7, over_price=150, under_price=150, line=8.5)  # 0.4+0.4=0.8
    ladder = m._anchor_ladder(st, 108549)
    assert ladder == {} or ladder.get(8.5, {}).get("reject") == "anchor_crossed"
    assert m._anchor_totals(st, 108549) == {}
    kalshi_event = {
        "event_ticker": "KXMLBTOTAL-26JUL251805NYYPHI",
        "title": "New York Y vs Philadelphia: Total Runs",
        "markets": [{"ticker": "KXMLBTOTAL-26JUL251805NYYPHI-8",
                     "strike_type": "greater", "floor_strike": 8.5}],
    }
    assert m.price_event(st, st.events[108549], kalshi_event) == []


def test_blown_vig_main_rung_rejected_for_mlb():
    """Same gate, the other rejection reason: implausible vig (sum > 1.20)
    on an MLB main rung is refused with the distinct band reason."""
    m = Mlb()
    st = build_bt3_state("lg5", "Philadelphia Phillies", "New York Yankees",
                          eid=108549, ms=7, over_price=-200, under_price=-200, line=8.5)  # 0.667*2
    # ANCHOR_SOURCE_IDS falls through past book 7 (only book with data, and its
    # only rung is rejected) to books 6/68 (no data at all) -> {}, same as an
    # all-rejected book on the soccer path.
    assert m._anchor_totals(st, 108549) == {}
    kalshi_event = {
        "event_ticker": "KXMLBTOTAL-26JUL251805NYYPHI",
        "title": "New York Y vs Philadelphia: Total Runs",
        "markets": [{"ticker": "KXMLBTOTAL-26JUL251805NYYPHI-8",
                     "strike_type": "greater", "floor_strike": 8.5}],
    }
    assert m.price_event(st, st.events[108549], kalshi_event) == []
