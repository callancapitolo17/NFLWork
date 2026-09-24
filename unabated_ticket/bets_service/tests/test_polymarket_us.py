"""Polymarket US parser on tests/fixtures/bets/polymarket_us_account.json, and
PolymarketUSSource's signing, paging, event cache and never-partial rules against a fake session."""
import base64
import copy
import json
from datetime import datetime, timezone
from pathlib import Path

import pytest
from cryptography.hazmat.primitives import serialization
from cryptography.hazmat.primitives.asymmetric import ed25519

from unabated_ticket.bets_service.sources import polymarket_us
from unabated_ticket.bets_service.sources.polymarket_us import (
    REASON_COMBO, REASON_NO_FILLS, PolymarketUSSource, fill_of, normalize_polymarket_us, parse_iso, parse_market_type, team_name_of)

FIXTURE_PATH = Path(__file__).parents[2] / "tests" / "fixtures" / "bets" / "polymarket_us_account.json"
FETCHED_AT = "2026-09-23T23:00:00Z"
CONTRACT_KEYS = {"id", "source", "venue", "league", "eventStart", "eventDate", "awayTeam", "homeTeam",
                 "awayKey", "homeKey", "rotation", "betType", "period", "side", "points", "price", "stake",
                 "toWin", "contracts", "placedAt", "status", "closedAt", "isParlayLeg", "parlayId",
                 "legIndex", "legCount", "approx", "unmatchable", "sourceFetchedAt", "venueIds", "raw"}
F5_SLUG = "asc-mlb-tor-bal-2026-09-23-f5-neg-1pt5"


@pytest.fixture
def account() -> dict:
    return json.loads(FIXTURE_PATH.read_text())


@pytest.fixture
def records(account) -> list[dict]:
    return normalize_polymarket_us(account["positions"], account["activities"], account["events"], FETCHED_AT)


def by_id(records: list[dict], record_id: str) -> dict:
    found = [record for record in records if record["id"] == record_id]
    assert len(found) == 1, f"record {record_id} missing"
    return found[0]


# ---- parser: the live pull ------------------------------------------------------------

def test_open_nfl_total_held_no_is_the_under_at_its_own_price(records):
    record = by_id(records, "polymarket_us:tsc-nfl-min-tb-2026-09-27-total-46pt5:no")
    assert (record["status"], record["closedAt"]) == ("open", None)
    assert (record["league"], record["betType"], record["period"], record["side"], record["points"]) == \
        ("nfl", "total", "FG", "under", 46.5)
    assert (record["awayTeam"], record["homeTeam"]) == ("Minnesota Vikings", "Tampa Bay Buccaneers")
    assert (record["awayKey"], record["homeKey"]) == (None, None)  # the panel resolves keys
    # Three fills at YES 0.545 = NO 0.455; the open size comes from the position (1.3), not the fills.
    assert (record["contracts"], record["price"], record["stake"], record["toWin"]) == (1.3, 120, 0.59, 0.71)
    assert (record["eventStart"], record["eventDate"]) == ("2026-09-27T20:05:00Z", "2026-09-27")
    assert record["placedAt"] == "2026-09-15T20:22:51Z"
    assert record["venueIds"] == {"marketSlug": "tsc-nfl-min-tb-2026-09-27-total-46pt5", "eventSlug": "nfl-min-tb-2026-09-27"}
    assert (record["source"], record["venue"], record["sourceFetchedAt"]) == ("polymarket_us_api", "polymarket_us", FETCHED_AT)
    assert (record["approx"], record["unmatchable"]) == ([], None)


def test_settled_moneylines_follow_the_resolution_side(records):
    lost_yes = by_id(records, "polymarket_us:aec-wnba-conn-wsh-2026-09-22:yes")
    assert (lost_yes["status"], lost_yes["side"], lost_yes["awayTeam"], lost_yes["homeTeam"]) == \
        ("lost", "away", "Connecticut", "Washington")
    assert (lost_yes["contracts"], lost_yes["price"], lost_yes["stake"], lost_yes["toWin"]) == (600, 1329, 42, 558)
    assert lost_yes["closedAt"] == "2026-09-23T01:34:16Z"  # the settlement's time
    lost_no = by_id(records, "polymarket_us:aec-wnba-ind-tor-2026-08-18:no")
    assert (lost_no["status"], lost_no["side"], lost_no["homeTeam"], lost_no["stake"]) == ("lost", "home", "Toronto", 72)
    won_yes = by_id(records, "polymarket_us:aec-wnba-por-sea-2026-08-14:yes")
    assert (won_yes["status"], won_yes["side"], won_yes["price"], won_yes["toWin"]) == ("won", "away", 133, 28.5)
    # 10 PM Eastern on the 14th is 02:00 UTC on the 15th: eventDate is the Eastern date.
    assert (won_yes["eventStart"], won_yes["eventDate"]) == ("2026-08-15T02:00:00Z", "2026-08-14")


def test_spread_held_no_is_the_other_team_with_its_own_signed_number(records):
    record = by_id(records, "polymarket_us:asc-nfl-bal-min-2026-08-22-pos-2pt5:no")
    assert (record["status"], record["betType"], record["side"], record["points"]) == ("lost", "spread", "home", -2.5)
    assert (record["awayTeam"], record["homeTeam"]) == ("Baltimore Ravens", "Minnesota Vikings")
    assert (record["contracts"], record["price"], record["stake"]) == (7.6, -116, 4.08)


def test_a_side_sold_back_to_zero_is_closed_at_its_last_fill(records):
    record = by_id(records, "polymarket_us:aec-wnba-dal-sea-2026-09-23:no")
    assert (record["status"], record["contracts"], record["stake"], record["toWin"]) == ("closed", 0, 0, 0)
    assert record["price"] == 400  # the buys' VWAP (NO at 0.20); the sells at 0.21 do not move it
    assert (record["closedAt"], record["raw"]["fillCount"], record["raw"]["netContracts"]) == ("2026-09-23T05:06:21Z", 10, 0)
    assert (record["side"], record["homeTeam"]) == ("home", "Seattle")


def test_a_combo_is_one_unmatchable_record_with_its_legs_in_raw(records):
    record = by_id(records, "polymarket_us:caoc-97039773d2e5d3cb:yes")
    assert (record["status"], record["unmatchable"]) == ("won", REASON_COMBO)
    assert (record["league"], record["betType"], record["side"], record["isParlayLeg"]) == (None, "other", None, False)
    assert (record["contracts"], record["price"], record["stake"], record["toWin"]) == (193.57, 920, 18.97, 174.6)
    legs = record["raw"]["comboLegs"]
    assert len(legs) == 3
    assert legs[0] == {"slug": "aec-nfl-gb-min-2026-09-13", "eventSlug": "nfl-gb-min-2026-09-13", "outcome": "Vikings",
                       "outcomeSide": "OUTCOME_SIDE_NO", "eventStartTime": "2026-09-13T20:25:00Z",
                       "state": "COMBO_LEG_STATE_WON"}


# ---- parser: hand-written shapes ------------------------------------------------------

def test_first_five_run_line_whose_yes_side_is_the_favourite_and_a_busted_trade_is_skipped(records):
    record = by_id(records, f"polymarket_us:{F5_SLUG}:yes")
    assert (record["status"], record["league"], record["betType"], record["period"]) == ("open", "mlb", "spread", "F5")
    assert (record["side"], record["points"], record["awayTeam"]) == ("away", -1.5, "Toronto Blue Jays")
    # The busted 1,000 at 0.90 would drag the VWAP to ~0.895 (-852) if it counted.
    assert (record["contracts"], record["price"], record["raw"]["fillCount"]) == (10, 150, 1)


def test_a_team_total_fails_closed_with_the_market_type(records):
    [record] = [r for r in records if r["raw"]["marketType"] == "football_team_first_half_total"]
    assert (record["status"], record["unmatchable"]) == ("open", "team total (football_team_first_half_total)")
    assert (record["league"], record["betType"], record["contracts"]) == (None, "other", 5)


def test_a_neutral_settlement_is_unknown_and_closes_at_the_settlement(records):
    [record] = [r for r in records if r["raw"]["resolutionSide"] == "POSITION_RESOLUTION_SIDE_NEUTRAL"]
    assert (record["status"], record["closedAt"], record["contracts"]) == ("unknown", "2026-09-25T04:00:00Z", 20)
    assert (record["league"], record["betType"], record["side"], record["awayTeam"]) == \
        ("nfl", "moneyline", "away", "Atlanta Falcons")


def test_soccer_fails_closed_and_a_deposit_is_not_a_bet(records):
    [record] = [r for r in records if r["raw"]["marketType"] == "soccer_team_full_game_total"]
    assert record["unmatchable"] == "market type not supported (soccer_team_full_game_total)"
    # No position and no settlement while the fills still hold contracts: open (the Kalshi rule —
    # a false open undersizes the next bet, a false settle would oversize it).
    assert (record["status"], record["contracts"], record["closedAt"]) == ("open", 3, None)
    assert len(records) == 11  # 7 live (slug, side) groups + 4 hand-written; the deposit makes none


def test_every_record_carries_the_contract_keys(records):
    for record in records:
        assert set(record) == CONTRACT_KEYS, record["id"]
        assert record["status"] in {"open", "won", "lost", "closed", "unknown"}
        assert (record["status"] == "open") == (record["closedAt"] is None), record["id"]


def test_an_open_position_with_no_fill_in_the_history_is_listed_unpriced(account):
    positions = dict(account["positions"])
    positions["aec-nfl-kc-mia-2026-09-27"] = {
        "netPositionDecimal": "-12.0000", "updateTime": "2026-08-01T00:00:00Z",
        "marketMetadata": {"slug": "aec-nfl-kc-mia-2026-09-27", "title": "KC Chiefs vs MIA Dolphins",
                           "outcome": "Dolphins", "eventSlug": "nfl-kc-mia-2026-09-27"}}
    records = normalize_polymarket_us(positions, account["activities"], account["events"], FETCHED_AT)
    record = by_id(records, "polymarket_us:aec-nfl-kc-mia-2026-09-27:no")
    assert (record["status"], record["contracts"], record["unmatchable"]) == ("open", 12, REASON_NO_FILLS)
    assert (record["price"], record["stake"], record["placedAt"], record["closedAt"]) == (None, None, None, None)
    assert record["venueIds"] == {"marketSlug": "aec-nfl-kc-mia-2026-09-27", "eventSlug": "nfl-kc-mia-2026-09-27"}


def test_a_game_without_its_event_payload_fails_closed_naming_it(account):
    events = dict(account["events"])
    del events["nfl-min-tb-2026-09-27"]
    records = normalize_polymarket_us(account["positions"], account["activities"], events, FETCHED_AT)
    record = by_id(records, "polymarket_us:tsc-nfl-min-tb-2026-09-27-total-46pt5:no")
    assert record["unmatchable"] == "unreadable Polymarket US event (no event payload for nfl-min-tb-2026-09-27)"
    assert record["status"] == "open"  # still listed, just not matched


def test_an_away_home_that_disagrees_with_a_side_fails_closed(account):
    events = copy.deepcopy(account["events"])
    events["nfl-bal-min-2026-08-22"]["teams"].reverse()  # the sides say BAL away, MIN home
    records = normalize_polymarket_us(account["positions"], account["activities"], events, FETCHED_AT)
    record = by_id(records, "polymarket_us:asc-nfl-bal-min-2026-08-22-pos-2pt5:no")
    assert record["unmatchable"].endswith("a side's away/home disagrees with the event)")


def test_an_undocumented_trade_state_fails_the_poll(account):
    trade = copy.deepcopy(account["activities"][0]["trade"])
    trade["state"] = "TRADE_STATE_SOMETHING_NEW"
    with pytest.raises(RuntimeError, match="TRADE_STATE_SOMETHING_NEW"):
        fill_of(trade)


def test_an_own_order_without_its_side_fails_the_poll(account):
    trade = copy.deepcopy(account["activities"][0]["trade"])
    del trade["aggressor"]["outcomeSide"]
    with pytest.raises(RuntimeError, match="own order without outcomeSide/action"):
        fill_of(trade)


@pytest.mark.parametrize("market_type, expected", [
    ("baseball_team_full_game_winner", {"betType": "moneyline", "period": "FG", "sport": "baseball"}),
    ("baseball_team_first_five_total", {"betType": "total", "period": "F5", "sport": "baseball"}),
    ("football_game_first_quarter_total", {"betType": "total", "period": "1Q", "sport": "football"}),
    ("football_team_second_half_spread", {"betType": "spread", "period": "2H", "sport": "football"}),
    ("hockey_team_full_game_winner", {"betType": "moneyline", "period": "FG", "sport": "hockey"}),
    ("football_team_points_full_game_total", "market type not supported (football_team_points_full_game_total)"),
    ("basketball_team_first_five_spread", "market type not supported (basketball_team_first_five_spread)"),
    ("soccer_team_full_time_winner", "market type not supported (soccer_team_full_time_winner)"),
    (None, "market type not supported (none)"),
])
def test_parse_market_type(market_type, expected):
    assert parse_market_type(market_type) == expected


def test_team_spelling_per_league():
    liberty = {"name": "Flames", "safeName": "Liberty", "alias": "Flames"}
    ottawa = {"name": "Senators", "safeName": "Ottawa Senators"}
    falcons = {"name": "Atlanta Falcons", "safeName": "ATL Falcons"}
    assert (team_name_of(liberty, "cfb"), team_name_of(ottawa, "nhl"), team_name_of(falcons, "nfl")) == \
        ("Liberty", "Ottawa Senators", "Atlanta Falcons")


def test_parse_iso_reads_nanoseconds():
    assert parse_iso("2026-09-15T20:22:51.500805479Z") == datetime(2026, 9, 15, 20, 22, 51, 500805, tzinfo=timezone.utc)
    assert parse_iso("2026-09-27T20:05:00Z") == datetime(2026, 9, 27, 20, 5, tzinfo=timezone.utc)
    assert parse_iso("not a time") is None


# ---- network half ---------------------------------------------------------------------

NOW = datetime(2026, 9, 23, 23, 0, tzinfo=timezone.utc).timestamp()
KEY_ID = "00000000-0000-4000-8000-000000000000"


def make_secret() -> tuple[str, ed25519.Ed25519PublicKey]:
    """A Polymarket-shaped secret: base64 of the 32-byte seed followed by the public key."""
    private = ed25519.Ed25519PrivateKey.generate()
    seed = private.private_bytes(serialization.Encoding.Raw, serialization.PrivateFormat.Raw,
                                 serialization.NoEncryption())
    public = private.public_key()
    public_raw = public.public_bytes(serialization.Encoding.Raw, serialization.PublicFormat.Raw)
    return base64.b64encode(seed + public_raw).decode(), public


class FakeResponse:
    def __init__(self, status_code: int, body: object):
        self.status_code = status_code
        self._body = body
        self.text = json.dumps(body)

    def json(self) -> object:
        return self._body


class FakeSession:
    """Answers the two signed account GETs page by page and the public event lookup."""

    def __init__(self, account: dict, position_pages: list[dict] | None = None,
                 activity_pages: list[list[dict]] | None = None, event_status: int = 200):
        self.account = account
        self.position_pages = position_pages or [account["positions"]]
        self.activity_pages = activity_pages or [account["activities"]]
        self.event_status = event_status
        self.status_overrides: dict[str, int] = {}
        self.calls: list[tuple[str, dict, dict]] = []

    def _page(self, pages: list, key: str, params: dict) -> dict:
        index = int(params.get("cursor") or 0)
        last = index == len(pages) - 1
        return {key: pages[index], "nextCursor": "" if last else str(index + 1), "eof": last}

    def get(self, url, params=None, headers=None, timeout=None):
        params = params or {}
        self.calls.append((url, dict(params), dict(headers or {})))
        path = url.split(".us", 1)[1]
        if path in self.status_overrides:
            return FakeResponse(self.status_overrides[path], {"message": "Invalid API key signature"})
        if url.startswith(polymarket_us.API_BASE_URL) and path == polymarket_us.POSITIONS_PATH:
            return FakeResponse(200, self._page(self.position_pages, "positions", params))
        if url.startswith(polymarket_us.API_BASE_URL) and path == polymarket_us.ACTIVITIES_PATH:
            return FakeResponse(200, self._page(self.activity_pages, "activities", params))
        if url.startswith(polymarket_us.GATEWAY_BASE_URL) and path == polymarket_us.EVENTS_PATH:
            if self.event_status != 200:
                return FakeResponse(self.event_status, {"message": "boom"})
            event = self.account["events"].get(params["slug"])
            return FakeResponse(200, {"events": [event] if event else []})
        raise AssertionError(f"unexpected GET {url}")

    def calls_to(self, base: str, path: str) -> list[tuple[str, dict, dict]]:
        return [call for call in self.calls if call[0] == f"{base}{path}"]


def make_source(session: FakeSession, secret: str, **kwargs) -> PolymarketUSSource:
    return PolymarketUSSource(key_id=KEY_ID, secret_key=secret, history_days=31, poll_sec=60,
                              session_factory=lambda: session, clock=lambda: NOW, **kwargs)


def test_fetch_signs_the_path_without_the_query_and_matches_the_pure_parser(account):
    secret, public = make_secret()
    session = FakeSession(account)
    records = make_source(session, secret).fetch()
    expected = normalize_polymarket_us(account["positions"], account["activities"], account["events"], None)
    strip = lambda rows: sorted((dict(row, sourceFetchedAt=None) for row in rows), key=lambda row: row["id"])  # noqa: E731
    assert strip(records) == strip(expected)
    for path in (polymarket_us.POSITIONS_PATH, polymarket_us.ACTIVITIES_PATH):
        [(_url, params, headers)] = session.calls_to(polymarket_us.API_BASE_URL, path)
        assert params == {"limit": polymarket_us.PAGE_LIMIT}
        assert headers["X-PM-Access-Key"] == KEY_ID
        assert headers["X-PM-Timestamp"] == str(int(NOW * 1000))
        # Verifies over timestamp + GET + path; a query string in the message would fail here.
        public.verify(base64.b64decode(headers["X-PM-Signature"]), f"{headers['X-PM-Timestamp']}GET{path}".encode())
    event_calls = session.calls_to(polymarket_us.GATEWAY_BASE_URL, polymarket_us.EVENTS_PATH)
    assert all("X-PM-Access-Key" not in headers for _url, _params, headers in event_calls)  # public: no key
    by_slug = {params["slug"]: params for _url, params, _headers in event_calls}
    assert by_slug["nfl-min-tb-2026-09-27"] == {"slug": "nfl-min-tb-2026-09-27", "limit": 1,
                                                "sportsMarketTypes": "football_team_full_game_total"}


def test_events_are_looked_up_once_and_a_failed_lookup_retries_next_poll(account):
    secret, _public = make_secret()
    session = FakeSession(account, event_status=500)
    source = make_source(session, secret)
    first = source.fetch()  # every game fails closed, the poll itself succeeds
    assert by_id(first, "polymarket_us:tsc-nfl-min-tb-2026-09-27-total-46pt5:no")["unmatchable"].startswith(
        "unreadable Polymarket US event")
    n_events = len({params["slug"] for _url, params, _headers in
                    session.calls_to(polymarket_us.GATEWAY_BASE_URL, polymarket_us.EVENTS_PATH)})
    session.event_status = 200
    second = source.fetch()
    assert by_id(second, "polymarket_us:tsc-nfl-min-tb-2026-09-27-total-46pt5:no")["unmatchable"] is None
    lookups_after_second = len(session.calls_to(polymarket_us.GATEWAY_BASE_URL, polymarket_us.EVENTS_PATH))
    source.fetch()
    assert len(session.calls_to(polymarket_us.GATEWAY_BASE_URL, polymarket_us.EVENTS_PATH)) == lookups_after_second
    assert lookups_after_second == 2 * n_events


def test_positions_and_activities_follow_the_cursor_to_eof(account):
    secret, _public = make_secret()
    slugs = sorted(account["positions"])
    position_pages = [{slugs[0]: account["positions"][slugs[0]]},
                      {slug: account["positions"][slug] for slug in slugs[1:]}]
    session = FakeSession(account, position_pages=position_pages,
                          activity_pages=[account["activities"][:10], account["activities"][10:]])
    records = make_source(session, secret).fetch()
    assert len(session.calls_to(polymarket_us.API_BASE_URL, polymarket_us.POSITIONS_PATH)) == 2
    assert [params.get("cursor") for _url, params, _headers in
            session.calls_to(polymarket_us.API_BASE_URL, polymarket_us.ACTIVITIES_PATH)] == [None, "1"]
    assert len(records) == 11


def test_activity_paging_stops_once_a_page_ends_before_the_history_window(account):
    secret, _public = make_secret()
    old = {"type": "ACTIVITY_TYPE_ACCOUNT_DEPOSIT",
           "accountBalanceChange": {"createTime": "2026-07-01T00:00:00Z", "amount": {"value": "0"}}}
    session = FakeSession(account, activity_pages=[account["activities"] + [old], [old]])
    make_source(session, secret).fetch()
    assert len(session.calls_to(polymarket_us.API_BASE_URL, polymarket_us.ACTIVITIES_PATH)) == 1


def old_page_split(account: dict, slug: str) -> list[list[dict]]:
    """Page 1: every activity but `slug`'s trades, ending on one older than the window;
    page 2: `slug`'s trades."""
    old = {"type": "ACTIVITY_TYPE_ACCOUNT_DEPOSIT",
           "accountBalanceChange": {"createTime": "2026-07-01T00:00:00Z", "amount": {"value": "0"}}}
    is_slug_trade = lambda a: a["type"] == "ACTIVITY_TYPE_TRADE" and a["trade"]["marketSlug"] == slug  # noqa: E731
    return [[a for a in account["activities"] if not is_slug_trade(a)] + [old],
            [a for a in account["activities"] if is_slug_trade(a)]]


@pytest.mark.parametrize("slug, record_id, status", [
    # Open: the position needs its fills for a price.
    ("tsc-nfl-min-tb-2026-09-27-total-46pt5", "polymarket_us:tsc-nfl-min-tb-2026-09-27-total-46pt5:no", "open"),
    # Just settled: the settlement on page 1 needs its fills, or the store's open row never closes.
    ("aec-wnba-conn-wsh-2026-09-22", "polymarket_us:aec-wnba-conn-wsh-2026-09-22:yes", "lost"),
])
def test_paging_goes_past_the_window_until_open_and_settling_positions_have_their_fills(account, slug, record_id, status):
    secret, _public = make_secret()
    session = FakeSession(account, activity_pages=old_page_split(account, slug))
    records = make_source(session, secret).fetch()
    assert len(session.calls_to(polymarket_us.API_BASE_URL, polymarket_us.ACTIVITIES_PATH)) == 2
    record = by_id(records, record_id)
    assert (record["status"], record["unmatchable"]) == (status, None)
    assert record["price"] is not None


def test_a_cursor_that_never_ends_fails_the_poll(account, monkeypatch):
    secret, _public = make_secret()
    monkeypatch.setattr(polymarket_us, "MAX_PAGES", 3)
    session = FakeSession(account, activity_pages=[account["activities"][:1]] * 5)
    with pytest.raises(RuntimeError, match="no eof after 3 pages"):
        make_source(session, secret).fetch()


def test_a_refused_signature_fails_the_poll_with_the_status(account):
    secret, _public = make_secret()
    session = FakeSession(account)
    session.status_overrides[polymarket_us.POSITIONS_PATH] = 401
    with pytest.raises(RuntimeError, match="GET /v1/portfolio/positions: expected HTTP 200, got 401"):
        make_source(session, secret).fetch()


def test_missing_or_malformed_keys_name_the_fix(account):
    session = FakeSession(account)
    with pytest.raises(RuntimeError, match="POLYMARKET_US_KEY_ID and POLYMARKET_US_SECRET_KEY"):
        PolymarketUSSource(key_id="", secret_key="", session_factory=lambda: session, clock=lambda: NOW).fetch()
    with pytest.raises(RuntimeError, match="not a base64 Ed25519 key"):
        make_source(session, "not base64 !!").fetch()


def test_source_if_configured_needs_both_keys(monkeypatch):
    monkeypatch.setattr(polymarket_us.config, "POLYMARKET_US_KEY_ID", KEY_ID)
    monkeypatch.setattr(polymarket_us.config, "POLYMARKET_US_SECRET_KEY", None)
    assert polymarket_us.source_if_configured() is None
    monkeypatch.setattr(polymarket_us.config, "POLYMARKET_US_SECRET_KEY", make_secret()[0])
    source = polymarket_us.source_if_configured()
    assert isinstance(source, PolymarketUSSource)
    assert (source.name, source.poll_sec) == ("polymarket_us", 60.0)
