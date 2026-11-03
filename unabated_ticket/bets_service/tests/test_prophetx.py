"""ProphetX (2026-10-31): how a game's markets read as board lines
(prophetx_markets.py), the bet source's records (sources/prophetx.py), the
client's session and order body (prophetx_client.py), and the Ticket's Place
desk (prophetx_place.py) against a scripted client and runner quote. Shapes are
built from the docs' v4 schemas (no live capture yet). No network."""
import json
from datetime import datetime, timedelta, timezone

import pytest

from unabated_ticket.bets_service import prophetx_markets
from unabated_ticket.bets_service.prophetx_client import ProphetXBlocked, ProphetXClient, ProphetXError, fill_or_kill_body
from unabated_ticket.bets_service.prophetx_place import (ManualProphetXDesk, settle_unconfirmed_place_orders,
                                                         teams_mismatch,
                                                         validate_prophetx_place_request)
from unabated_ticket.bets_service.sources.prophetx import ProphetXSource, normalize_prophetx
from unabated_ticket.bets_service.store import BetsStore

NOW = datetime(2026, 11, 1, 15, 0, tzinfo=timezone.utc)
PX_EVENT_ID = 5501
EVENT = {"event_id": PX_EVENT_ID, "name": "Buffalo Bills at Kansas City Chiefs", "scheduled": "2026-11-01T21:25:00Z",
         "status": "not_started", "tournament_id": 31, "tournament_name": "NFL",
         "competitors": [{"id": 11, "name": "Kansas City Chiefs", "side": "home"},
                         {"id": 12, "name": "Buffalo Bills", "side": "away"}]}


def offer(strike_id, price, quantity, competitor_id=None, name="", strike=None, display_strike=None, outcome_id=1):
    return {"outcome_id": outcome_id, "competitor_id": competitor_id, "name": name, "price": price, "strike": strike,
            "display_strike": display_strike, "strike_id": strike_id, "quantity": quantity}


MONEYLINE = {"id": 900, "name": "Moneyline", "type": "moneyline", "selections": [
    [offer("a" * 32, 150, 300, competitor_id=12, name="Buffalo Bills")],
    [offer("b" * 32, -160, 500, competitor_id=11, name="Kansas City Chiefs")],
]}
SPREAD = {"id": 901, "name": "Spread", "type": "spread", "selections": [
    [offer("c" * 32, -110, 200, 12, "Buffalo Bills", 3.5, "+3.5"), offer("c" * 32, -115, 900, 12, "Buffalo Bills", 3.5, "+3.5")],
    [offer("d" * 32, -108, 400, 11, "Kansas City Chiefs", -3.5, "-3.5")],
], "market_strikes": [{"id": 902, "name": "Spread", "type": "spread", "selections": [
    [offer("e" * 32, 120, 50, 12, "Buffalo Bills", 2.5, "+2.5")],
    [offer("f" * 32, -140, 50, 11, "Kansas City Chiefs", -2.5, "-2.5")],
]}]}
TOTAL = {"id": 903, "name": "Total Points", "type": "total", "selections": [
    [offer("g" * 32, -105, 100, None, "Over 47.5", 47.5)],
    [offer("h" * 32, -115, 100, None, "Under 47.5", 47.5)],
]}
FIRST_HALF_SPREAD = {"id": 904, "name": "1st Half Spread", "type": "spread", "selections": [
    [offer("i" * 32, -110, 100, 12, "Buffalo Bills", 1.5, "+1.5")],
    [offer("j" * 32, -110, 100, 11, "Kansas City Chiefs", -1.5, "-1.5")],
]}
QUARTER = {"id": 905, "name": "1st Quarter Spread", "type": "spread", "selections": MONEYLINE["selections"]}
MARKETS = [MONEYLINE, SPREAD, TOTAL, FIRST_HALF_SPREAD, QUARTER]


# ---- markets ---------------------------------------------------------------------------


def test_every_game_market_side_reads_as_a_board_line():
    lines = prophetx_markets.lines_of_event(EVENT, MARKETS)
    by_strike = {line["strike_id"]: (line["betType"], line["period"], line["side"], line["points"]) for line in lines}
    assert by_strike == {
        "a" * 32: ("moneyline", "FG", "away", None), "b" * 32: ("moneyline", "FG", "home", None),
        "c" * 32: ("spread", "FG", "away", 3.5), "d" * 32: ("spread", "FG", "home", -3.5),
        "e" * 32: ("spread", "FG", "away", 2.5), "f" * 32: ("spread", "FG", "home", -2.5),
        "g" * 32: ("total", "FG", "over", 47.5), "h" * 32: ("total", "FG", "under", 47.5),
        "i" * 32: ("spread", "1H", "away", 1.5), "j" * 32: ("spread", "1H", "home", -1.5),
    }, "the quarter market is left out, never read as a full game"
    assert {line["awayTeam"] for line in lines} == {"Buffalo Bills"}


def test_a_main_number_listed_again_under_its_alt_numbers_is_one_line():
    repeated = {**SPREAD, "market_strikes": SPREAD["market_strikes"] + [{**SPREAD, "market_strikes": None}]}
    strikes = [line["strike_id"] for line in prophetx_markets.lines_of_event(EVENT, [repeated])]
    assert sorted(strikes) == ["c" * 32, "d" * 32, "e" * 32, "f" * 32]


def test_a_spread_whose_sides_carry_one_unsigned_number_does_not_read():
    one_number = {"id": 1, "name": "Spread", "type": "spread", "selections": [
        [offer("x" * 32, -110, 10, 12, "Buffalo Bills", 3.5)], [offer("y" * 32, -110, 10, 11, "Kansas City Chiefs", 3.5)]]}
    assert prophetx_markets.lines_of_event(EVENT, [one_number]) == []


def test_a_game_without_home_and_away_or_in_another_league_is_named():
    assert "home and away" in prophetx_markets.lines_of_event(
        {**EVENT, "competitors": [{"id": 1, "name": "A"}, {"id": 2, "name": "B"}]}, MARKETS)
    assert "league not supported" in prophetx_markets.lines_of_event(
        {**EVENT, "tournament_id": 7, "tournament_name": "Bundesliga"}, MARKETS)


def test_the_best_offer_is_the_highest_price_with_the_money_resting_at_it():
    assert prophetx_markets.best_offer_of(SPREAD["selections"][0]) == {"price": -110, "depth": 200}
    assert prophetx_markets.best_offer_of([offer("z", -110, 0)]) is None
    assert prophetx_markets.american_to_probability(-110) == pytest.approx(0.5238095)


# ---- bet source ------------------------------------------------------------------------


def order(order_id, strike_id, filled, status="open", winning="tbd", open_quantity=0.0, price=-110, **extra):
    return {"order_id": order_id, "external_id": f"ext-{order_id}", "sport_event_id": PX_EVENT_ID, "market_id": 901,
            "strike_id": strike_id, "price": price, "quantity": filled + open_quantity,
            "original_quantity": filled + open_quantity, "filled_quantity": filled, "open_quantity": open_quantity,
            "status": status, "matching_status": "fully_matched" if not open_quantity else "partially_matched",
            "winning_status": winning, "created_at": "2026-11-01T14:00:00Z", "updated_at": "2026-11-01T14:00:05Z",
            "settled_at": None, **extra}


def normalized(orders, trades=(), known=None):
    lines = {PX_EVENT_ID: prophetx_markets.lines_of_event(EVENT, MARKETS)}
    return {record["id"]: record for record in normalize_prophetx(list(orders), list(trades), {PX_EVENT_ID: EVENT},
                                                                  lines, known or {}, "2026-11-01T15:00:00Z")}


def test_a_matched_order_is_a_bet_on_its_line_priced_off_its_fills():
    trades = [{"order_id": "o1", "price": -110, "quantity": 110, "fee": 0.5},
              {"order_id": "o1", "price": 100, "quantity": 100, "fee": 0.5}]
    record = normalized([order("o1", "c" * 32, 210)], trades)["prophetx:o1"]
    assert (record["venue"], record["league"], record["betType"], record["period"], record["side"], record["points"]) == \
        ("prophetx", "nfl", "spread", "FG", "away", 3.5)
    assert (record["awayTeam"], record["homeTeam"], record["eventDate"]) == ("Buffalo Bills", "Kansas City Chiefs", "2026-11-01")
    # $110 at -110 wins $100, $100 at +100 wins $100: $210 to win $200, and $1 of fees: $211 to net $199.
    assert (record["stake"], record["toWin"], record["price"], record["status"]) == (211, 199, -106, "open")
    assert record["raw"]["fees"] == 1.0 and record["unmatchable"] is None


def test_a_resting_remainder_is_its_own_unmatchable_record():
    records = normalized([order("o2", "g" * 32, 40, open_quantity=60, price=-105)])
    assert records["prophetx:o2"]["stake"] == 40 and records["prophetx:o2"]["toWin"] == pytest.approx(38.1)
    resting = records["prophetx:o2:order"]
    assert (resting["stake"], resting["status"], resting["approx"]) == (60, "open", ["prophetx_unmatched_order"])
    assert resting["unmatchable"] == "unmatched order (resting $60.0 at -105)"


@pytest.mark.parametrize("status, winning, expected", [
    ("closed", "profit", "won"), ("closed", "loss", "lost"), ("settled", "push", "push"),
    ("manually_settled", "manually_lost", "lost"), ("closed", "no_result", "void"), ("closed", "tbd", "unknown"),
    ("canceled", "tbd", "open"), ("wiped", "loss", "lost"), ("void", "profit", "void"), ("invalid", "tbd", "void"),
])
def test_status_grades_the_matched_part(status, winning, expected):
    assert normalized([order("o3", "d" * 32, 50, status=status, winning=winning)])["prophetx:o3"]["status"] == expected


def test_an_order_with_nothing_matched_and_nothing_resting_is_no_record():
    assert normalized([order("o4", "d" * 32, 0, status="canceled")]) == {}


def test_a_game_prophetx_no_longer_lists_keeps_its_stored_line():
    stored = normalized([order("o5", "h" * 32, 30)])["prophetx:o5"]
    settled = order("o5", "h" * 32, 30, status="closed", winning="loss", settled_at="2026-11-02T01:00:00Z")
    records = {r["id"]: r for r in normalize_prophetx([settled], [], {}, {}, {"prophetx:o5": stored}, None)}
    assert (records["prophetx:o5"]["side"], records["prophetx:o5"]["points"], records["prophetx:o5"]["status"]) == \
        ("under", 47.5, "lost")
    assert records["prophetx:o5"]["closedAt"] == "2026-11-02T01:00:00Z"
    unknown = normalize_prophetx([settled], [], {}, {}, {}, None)[0]
    assert unknown["unmatchable"] == f"ProphetX does not list game {PX_EVENT_ID} any more"


def test_a_remainder_that_stops_resting_closes_its_stored_open_copy_as_void():
    stored = normalized([order("o7", "g" * 32, 40, open_quantity=60, price=-105)])
    matched_all = order("o7", "g" * 32, 100, status="closed", winning="tbd", updated_at="2026-11-01T14:30:00Z")
    records = normalized([matched_all], known=stored)
    ended = records["prophetx:o7:order"]
    assert (ended["status"], ended["stake"], ended["toWin"], ended["closedAt"]) == ("void", 0.0, 0.0, "2026-11-01T14:30:00Z")
    assert records["prophetx:o7"]["stake"] == 100
    # Once closed, it is not sent again.
    assert "prophetx:o7:order" not in normalized([matched_all], known=records)


class ScriptedReads:
    """`orders` answer the polled window (the one ending now); `older[k]` the k-th
    slice before it (k = 1 the newest), every other slice none."""
    def __init__(self, orders=None, listed=True, older=None, window_days=31):
        self.event_reads = []
        self.market_reads = []
        self.orders = orders if orders is not None else [order("o6", "b" * 32, 80, price=-160)]
        self.listed = listed
        self.older = older or {}
        self.window_sec = window_days * 86400
        self.history_reads = []
        self.trade_reads = []

    def _slice(self, start, end):
        # 0 = the polled window; k = the k-th older slice.
        return round((NOW.timestamp() - end) / self.window_sec)

    def order_history(self, start, end):
        assert end - start <= self.window_sec
        self.history_reads.append(self._slice(start, end))
        k = self._slice(start, end)
        return self.orders if k == 0 else self.older.get(k, [])

    def trades(self, start, end):
        self.trade_reads.append(self._slice(start, end))
        return []

    def sport_events(self, event_ids=None):
        self.event_reads.append(event_ids)
        return [EVENT] if self.listed else []

    def markets(self, event_id):
        self.market_reads.append(event_id)
        return MARKETS

    def balance(self):
        return 200, {"data": {"balance": 412.5, "matched_order_balance": 80.0}}


def test_the_source_reads_each_game_once_and_reads_the_balance():
    client = ScriptedReads()
    source = ProphetXSource(client, dict, history_days=31, clock=lambda: NOW.timestamp())
    [first] = source.fetch()
    source.fetch()
    assert first["id"] == "prophetx:o6" and first["side"] == "home" and first["betType"] == "moneyline"
    assert client.event_reads == [[PX_EVENT_ID]]
    assert source.balance() == {"cash": 412.5, "inBets": 80.0, "available": 412.5, "creditLimit": None}


def test_orders_older_than_the_window_are_read_once_per_boot_back_to_the_account_start():
    old = order("o1", "c" * 32, 9)
    client = ScriptedReads(older={1: [], 2: [old]})
    source = ProphetXSource(client, dict, history_days=31, clock=lambda: NOW.timestamp())
    first = {record["id"] for record in source.fetch()}
    # Slices 1 (empty), 2 (an order), then 3 and 4 empty in a row end the walk.
    assert client.history_reads == [1, 2, 3, 4, 0] and client.trade_reads == [1, 2, 3, 4, 0]
    second = {record["id"] for record in source.fetch()}
    assert first == second == {"prophetx:o6", "prophetx:o1"}
    assert client.history_reads == [1, 2, 3, 4, 0, 0]


def test_markets_are_read_again_only_for_a_strike_they_lack_after_the_ttl():
    now = [NOW.timestamp()]
    client = ScriptedReads()
    source = ProphetXSource(client, dict, history_days=31, clock=lambda: now[0])
    source.fetch()
    now[0] += 3600
    source.fetch()
    assert client.market_reads == [PX_EVENT_ID]
    client.orders = client.orders + [order("o8", "z" * 32, 10)]
    source.fetch()
    assert client.market_reads == [PX_EVENT_ID, PX_EVENT_ID]


def test_a_game_prophetx_no_longer_lists_is_asked_for_again_only_hourly():
    now = [NOW.timestamp()]
    client = ScriptedReads(listed=False)
    source = ProphetXSource(client, dict, history_days=31, clock=lambda: now[0])
    source.fetch()
    now[0] += 60
    source.fetch()
    assert client.event_reads == [[PX_EVENT_ID]] and client.market_reads == []
    now[0] += 3600
    source.fetch()
    assert client.event_reads == [[PX_EVENT_ID], [PX_EVENT_ID]]


# ---- client ----------------------------------------------------------------------------


def test_the_client_logs_in_once_per_token_and_again_after_a_401():
    sent = []
    answers = [(200, {"data": {"access_token": "t1"}}), (401, {}), (200, {"data": {"access_token": "t2"}}),
               (200, {"data": {"balance": 1}})]

    def transport(method, url, headers, body):
        sent.append((method, url, headers.get("Authorization"), json.loads(body) if body else None))
        status, reply = answers.pop(0)
        return status, json.dumps(reply).encode()

    client = ProphetXClient("ak", "sk", "https://px.test/partner", transport=transport, clock=lambda: 1000.0, min_send_interval_sec=0)
    assert client.balance() == (200, {"data": {"balance": 1}})
    assert [(method, url, auth) for method, url, auth, _ in sent] == [
        ("POST", "https://px.test/partner/auth/login", None), ("GET", "https://px.test/partner/v4/mm/get_balance", "Bearer t1"),
        ("POST", "https://px.test/partner/auth/login", None), ("GET", "https://px.test/partner/v4/mm/get_balance", "Bearer t2")]
    assert sent[0][3] == {"access_key": "ak", "secret_key": "sk"}


def scripted_transport(answers, sent):
    def transport(method, url, headers, body):
        sent.append((method, url.rsplit("/partner", 1)[-1], headers.get("Authorization"), json.loads(body) if body else None))
        status, reply = answers.pop(0)
        return status, reply if isinstance(reply, bytes) else json.dumps(reply).encode()
    return transport


def test_an_expired_token_is_refreshed_not_logged_in_again():
    now = [1000.0]
    sent = []
    answers = [(200, {"data": {"access_token": "t1", "refresh_token": "r1"}}), (200, {"data": {"balance": 1}}),
               (200, {"data": {"access_token": "t2"}}), (200, {"data": {"balance": 2}})]
    client = ProphetXClient("ak", "sk", "https://px.test/partner", transport=scripted_transport(answers, sent),
                            clock=lambda: now[0], min_send_interval_sec=0)
    client.balance()
    now[0] += 600
    assert client.balance() == (200, {"data": {"balance": 2}})
    assert [(path, auth, body) for _, path, auth, body in sent] == [
        ("/auth/login", None, {"access_key": "ak", "secret_key": "sk"}), ("/v4/mm/get_balance", "Bearer t1", None),
        ("/auth/refresh", None, {"refresh_token": "r1"}), ("/v4/mm/get_balance", "Bearer t2", None)]


def test_a_refused_refresh_falls_back_to_one_login():
    now = [1000.0]
    sent = []
    answers = [(200, {"data": {"access_token": "t1", "refresh_token": "r1"}}), (200, {"data": {}}),
               (401, {"error": "expired"}), (200, {"data": {"access_token": "t3", "refresh_token": "r3"}}),
               (200, {"data": {}})]
    client = ProphetXClient("ak", "sk", "https://px.test/partner", transport=scripted_transport(answers, sent),
                            clock=lambda: now[0], min_send_interval_sec=0)
    client.balance()
    now[0] += 600
    client.balance()
    assert [path for _, path, _, _ in sent] == ["/auth/login", "/v4/mm/get_balance", "/auth/refresh", "/auth/login",
                                               "/v4/mm/get_balance"]


def test_a_cloudfront_block_stops_every_call_for_an_hour_then_doubles():
    now = [1000.0]
    sent = []
    page = b"<HTML><H1>403 ERROR</H1>Request blocked. Generated by cloudfront</HTML>"
    answers = [(403, page), (403, page), (200, {"data": {"access_token": "t1"}}), (200, {"data": {}})]
    client = ProphetXClient("ak", "sk", "https://px.test/partner", transport=scripted_transport(answers, sent),
                            clock=lambda: now[0], min_send_interval_sec=0)
    with pytest.raises(ProphetXBlocked, match="CloudFront refused this server"):
        client.balance()
    now[0] += 3599
    with pytest.raises(ProphetXBlocked):
        client.balance()
    assert len(sent) == 1 and client.blocked_until == 1000.0 + 3600
    now[0] += 1
    with pytest.raises(ProphetXBlocked):
        client.balance()
    assert len(sent) == 2 and client.blocked_until == now[0] + 7200
    now[0] += 7200
    assert client.balance() == (200, {"data": {}})
    assert client.blocked_until == 0.0


def test_a_refused_login_backs_off_instead_of_retrying_each_poll():
    sent = []
    answers = [(403, {"error": "session_num_exceed", "message": "only 20 sessions"})]
    client = ProphetXClient("ak", "sk", "https://px.test/partner", transport=scripted_transport(answers, sent),
                            clock=lambda: 1000.0, min_send_interval_sec=0)
    with pytest.raises(ProphetXBlocked, match="session_num_exceed"):
        client.balance()
    with pytest.raises(ProphetXBlocked):
        client.balance()
    assert len(sent) == 1


def test_a_refused_login_raises_with_what_came_back():
    client = ProphetXClient("ak", "sk", transport=lambda *_: (401, b'{"message":"bad key"}'), min_send_interval_sec=0)
    with pytest.raises(ProphetXError, match="login refused: HTTP 401"):
        client.balance()


def test_the_order_body_is_one_fill_or_kill_buy_in_dollars():
    assert fill_or_kill_body("c" * 32, -110, 25.456, "xman-1") == {
        "strike_id": "c" * 32, "price": -110, "quantity": 25.46, "external_id": "xman-1", "order_strategy": "fillOrKill"}
    with pytest.raises(ValueError):
        fill_or_kill_body("c" * 32, -99, 10, "x")


# ---- Place desk ------------------------------------------------------------------------

EVENT_ID = 123742
LINE = {"league": "nfl", "betType": "spread", "period": "FG", "side": "away", "points": 3.5,
        "awayTeam": "Buffalo Bills", "homeTeam": "Kansas City Chiefs", "awayTeamId": "1", "homeTeamId": "2"}
QUOTE = {"ok": True, "eventId": EVENT_ID, "eventStart": "2026-11-01T21:25:00.000Z", "label": "Bills +3.5 FG",
         "price": -110, "line": LINE, "lineStrikeId": None, "strikeIds": ["e" * 32]}
BODY = {"clickId": "abc12345-0001", "eventId": EVENT_ID, "periodTypeId": 1, "betTypeId": 2, "sideIndex": 0,
        "points": 3.5, "stake": 50, "ticketAllIn": 0.5238, "prophetxReadAt": "2026-11-01T14:59:40Z",
        "kalshiReadAt": None, "novigReadAt": None, "polymarketReadAt": None}


class ScriptedProphetX:
    """ProphetXClient's calls, scripted: the order is taken and reads fully matched."""

    def __init__(self, submit=None, order_reads=None, events=None, markets=None):
        self.submits = []
        self.submit = submit
        self.order_reads = list(order_reads or [])
        self.history_windows = []
        self.events = [EVENT] if events is None else events
        self.markets_list = MARKETS if markets is None else markets

    def strikes(self, strike_ids):
        return [{"event_id": PX_EVENT_ID, "strike_id": strike_id} for strike_id in strike_ids]

    def sport_events(self, event_ids=None):
        return self.events

    def markets(self, event_id):
        return self.markets_list

    def submit_fill_or_kill(self, strike_id, price, stake, external_id):
        self.submits.append((strike_id, price, stake, external_id))
        if isinstance(self.submit, Exception):
            raise self.submit
        return self.submit or (200, {"data": {"success": True, "order": {
            "order_id": "px-1", "status": "open", "filled_quantity": 0, "open_quantity": stake}}})

    def order_history(self, created_from, created_to, max_pages=50):
        """Place's order reads: each scripted read is a list of history rows,
        or (200, {"data": order}) for the one order px-1, or an exception."""
        self.history_windows.append((created_from, created_to))
        if self.order_reads:
            read = self.order_reads.pop(0)
            if isinstance(read, Exception):
                raise read
            if isinstance(read, tuple):
                return [{"order_id": "px-1", **read[1]["data"]}]
            return read
        stake = self.submits[-1][2]
        return [{"order_id": "px-0", "status": "settled", "filled_quantity": 9, "open_quantity": 0},
                {"order_id": "px-1", "status": "open", "matching_status": "fully_matched", "filled_quantity": stake,
                 "open_quantity": 0, "fee": 0.25}]


@pytest.fixture
def store(tmp_path):
    held = BetsStore(tmp_path / "bets.duckdb", 0)
    held.log_source_run("prophetx", NOW - timedelta(seconds=40), NOW - timedelta(seconds=30), True, None, 3)
    yield held
    held.close()


# The runner's team index as /team_keys.json resolves names: ProphetX's spelling -> Unabated key.
TEAM_KEYS = {"Buffalo Bills": "nfl:1", "Kansas City Chiefs": "nfl:2", "Miami Dolphins": "nfl:3"}


def team_keys_reader(keys=None):
    keys = TEAM_KEYS if keys is None else keys

    def read(league, away, home):
        return {"ok": True, "awayKey": keys.get(away), "homeKey": keys.get(home)}
    return read


def desk_with(store, client, quote=None, now=NOW, team_keys=None):
    return ManualProphetXDesk(client, store, lambda _line: dict(QUOTE if quote is None else quote),
                              team_keys or team_keys_reader(), now=lambda: now,
                              sleep=lambda _sec: None, clock=iter(range(0, 1000)).__next__)


def request(**changes):
    return validate_prophetx_place_request({**BODY, **changes})


def test_validate_names_what_is_wrong():
    assert validate_prophetx_place_request({**BODY, "prophetxReadAt": None}).startswith("prophetxReadAt")
    assert validate_prophetx_place_request({**BODY, "stake": 0.5}).startswith("stake")
    assert validate_prophetx_place_request({**BODY, "betTypeId": 1})["points"] is None
    assert validate_prophetx_place_request({**BODY, "novigReadAt": "soon"}).startswith("novigReadAt")
    assert request()["stake"] == 50 and request()["kalshi_read_at"] is None


def test_a_click_finds_the_strike_by_the_games_ids_and_sends_one_fill_or_kill_at_the_best_price(store):
    client = ScriptedProphetX()
    result = desk_with(store, client).place(request())
    assert client.submits == [("c" * 32, -110, 50, "xman-abc12345-0001")]
    assert result["status"] == "filled" and result["fee"] == 0.25
    assert result["message"] == "Bet $50.25 on Bills +3.5 at -111 to win $45.20 (limit -110 on ProphetX, fee $0.25)."
    row = store._con.execute("SELECT status, prophetx_order_id, count_filled, cost_filled, ticker FROM manual_orders "
                             "WHERE client_order_id = 'xman-abc12345-0001'").fetchone()
    assert row == ("filled", row[1], 5000, 50.25, "c" * 32) and row[1]


def test_without_strike_ids_the_game_is_found_by_its_teams_and_start(store):
    client = ScriptedProphetX(events=[EVENT, {**EVENT, "event_id": 1, "scheduled": "2026-11-02T21:25:00Z"}])
    assert desk_with(store, client, {**QUOTE, "strikeIds": []}).place(request())["status"] == "filled"
    assert client.submits[0][0] == "c" * 32


def test_the_same_click_is_never_sent_twice(store):
    client = ScriptedProphetX()
    desk_with(store, client).place(request())
    again = desk_with(store, client).place(request(prophetxReadAt="2026-11-01T15:30:00Z"))
    assert again["status"] == "refused" and "already sent" in again["message"]
    assert len(client.submits) == 1


@pytest.mark.parametrize("quote_changes, client_changes, words", [
    ({"eventStart": "2026-11-01T14:00:00Z"}, {}, "has started"),
    ({}, {"events": [{**EVENT, "status": "live"}]}, "has started"),
    ({}, {"events": []}, "no longer lists"),
    ({"line": {**LINE, "awayTeam": "Kansas City Chiefs", "homeTeam": "Buffalo Bills", "awayTeamId": "2", "homeTeamId": "1"}},
     {}, "other way round"),
    ({"line": {**LINE, "points": 4.5}}, {}, "0 ProphetX strikes"),
    ({"line": {**LINE, "betType": "total", "side": "over", "points": 47.5, "awayTeam": "Miami Dolphins", "awayTeamId": "3"}},
     {}, "does not read as"),
    ({"line": {**LINE, "homeTeamId": None}}, {}, "names no team ids"),
    ({"lineStrikeId": "e" * 32}, {}, "is not the one"),
    ({"line": {**LINE, "betType": "moneyline", "period": "1H", "points": None}}, {}, "tied period"),
    ({"ok": False, "detail": "the board holds 0 ProphetX lines"}, {}, "the board holds 0"),
    ({}, {"markets": [{**SPREAD, "selections": [[offer("c" * 32, -110, 0, 12, "Buffalo Bills", 3.5, "+3.5")],
                                                SPREAD["selections"][1]]}]}, "nothing rests"),
])
def test_refusals_send_nothing(store, quote_changes, client_changes, words):
    client = ScriptedProphetX(**client_changes)
    result = desk_with(store, client, {**QUOTE, **quote_changes}).place(request())
    assert result["status"] == "refused" and words in result["message"], result["message"]
    assert client.submits == []


def test_a_price_more_than_a_tenth_of_a_cent_worse_refuses(store):
    client = ScriptedProphetX()
    result = desk_with(store, client).place(request(ticketAllIn=0.51))
    assert result["status"] == "refused" and "-110 now" in result["message"]


def test_holding_the_other_side_refuses(store):
    store.upsert_bets([{"id": "prophetx:held", "venue": "prophetx", "source": "prophetx_api", "status": "open",
                        "unmatchable": None, "venueIds": {"strikeId": "d" * 32}, "placedAt": "2026-11-01T12:00:00Z"}], NOW)
    result = desk_with(store, ScriptedProphetX()).place(request())
    assert result["status"] == "refused" and "other side" in result["message"]


def test_an_open_order_matching_has_not_reached_is_read_again_not_called_unfilled(store):
    unreached = (200, {"data": {"success": True, "order": {
        "order_id": "px-1", "status": "open", "filled_quantity": 0, "open_quantity": 0}}})
    client = ScriptedProphetX(submit=unreached)
    result = desk_with(store, client).place(request())
    assert result["status"] == "filled" and result["filled"] == 50


def test_a_stale_or_failing_prophetx_read_refuses(store):
    store.log_source_run("prophetx", NOW - timedelta(seconds=5), NOW, False, "HTTP 500", 0)
    result = desk_with(store, ScriptedProphetX()).place(request())
    assert result["status"] == "refused" and "read is failing" in result["message"]


def test_a_part_fill_reads_partial_and_nothing_filled_reads_unfilled(store):
    partial = ScriptedProphetX(order_reads=[(200, {"data": {"status": "canceled", "filled_quantity": 20, "open_quantity": 0}})])
    assert desk_with(store, partial).place(request())["status"] == "partial"
    none = ScriptedProphetX(order_reads=[(200, {"data": {"status": "canceled", "filled_quantity": 0, "open_quantity": 0}})])
    result = desk_with(store, none).place(request(clickId="abc12345-0002", prophetxReadAt="2026-11-01T15:00:01Z"))
    assert result["status"] == "unfilled" and "Nothing rests" in result["message"]


# 2026-11-03: two orders read "inactive" with money open for the whole window
# while ProphetX had matched them in full (its cancel answered order_already_filled).
STILL_MATCHING = (200, {"data": {"order_id": "px-1", "status": "inactive", "matching_status": "unmatched",
                                 "filled_quantity": 0, "open_quantity": 50}})
HISTORY_FILLED = {"order_id": "px-1", "external_id": "xman-abc12345-0001", "status": "open",
                  "matching_status": "fully_matched", "filled_quantity": 50, "open_quantity": 0, "price": -110}


def manual_row(store, click="xman-abc12345-0001"):
    return store._con.execute("SELECT status, prophetx_order_id, count_filled, cost_filled FROM manual_orders "
                              "WHERE client_order_id = ?", [click]).fetchone()


def test_an_order_not_read_finished_is_left_unknown_never_cancelled_and_says_so(store):
    client = ScriptedProphetX(order_reads=[STILL_MATCHING] * 10)
    result = desk_with(store, client).place(request())
    assert result["status"] == "unknown" and "still matching" in result["message"]
    assert "resting" not in result["message"] and not hasattr(client, "cancel_order")
    assert manual_row(store)[:2] == ("unknown", "px-1")


def test_the_order_history_settles_an_unknown_place_order_at_its_fills(store):
    desk_with(store, ScriptedProphetX(order_reads=[STILL_MATCHING] * 10)).place(request())
    trades = [{"order_id": "px-1", "quantity": 50, "price": -105, "fee": 0.25, "status": "open"}]
    settle_unconfirmed_place_orders(store, [HISTORY_FILLED], trades, NOW + timedelta(seconds=30))
    assert manual_row(store) == ("filled", "px-1", 5000, 50.25)
    # Settled rows are not looked for again.
    assert store.unconfirmed_prophetx_place_orders(NOW + timedelta(minutes=1)) == []


def test_the_history_leaves_an_order_still_matching_or_unlisted_unknown(store):
    desk_with(store, ScriptedProphetX(order_reads=[STILL_MATCHING] * 10)).place(request())
    settle_unconfirmed_place_orders(store, [], [], NOW + timedelta(seconds=30))
    settle_unconfirmed_place_orders(store, [STILL_MATCHING[1]["data"]], [], NOW + timedelta(seconds=30))
    assert manual_row(store)[0] == "unknown"
    # A day on, the row is no longer looked for.
    assert store.unconfirmed_prophetx_place_orders(NOW + timedelta(hours=25)) == []


def test_a_lost_submit_answer_is_settled_by_its_external_id(store):
    desk_with(store, ScriptedProphetX(submit=OSError("timed out"))).place(request())
    assert manual_row(store)[:2] == ("unknown", None)
    killed = {**HISTORY_FILLED, "status": "canceled", "matching_status": "partially_matched", "filled_quantity": 20,
              "open_quantity": 0}
    settle_unconfirmed_place_orders(store, [killed], [], NOW + timedelta(seconds=30))
    assert manual_row(store) == ("partial", "px-1", 2000, 20.0)


def test_each_source_poll_hands_its_order_history_to_the_place_settler():
    class History:
        def order_history(self, start, end):
            return [HISTORY_FILLED]

        def trades(self, start, end):
            return []

        def sport_events(self, event_ids=None):
            return []

    seen = []
    source = ProphetXSource(History(), dict, clock=lambda: NOW.timestamp(),
                            settle_place_orders=lambda orders, trades: seen.append((orders, trades)))
    source._older = ([], [])
    source.fetch()
    assert seen == [([HISTORY_FILLED], [])]


def test_place_reads_its_order_off_the_order_history_until_it_is_listed(store):
    # 2026-11-03: get_order 404'd for 8 s+ after a send; the history listed the order within 3 s.
    client = ScriptedProphetX(order_reads=[[], OSError("timed out")])
    result = desk_with(store, client).place(request())
    assert result["status"] == "filled" and result["filled"] == 50 and result["message"].startswith("Bet $50.25")
    assert len(client.history_windows) == 3
    start, end = client.history_windows[0]
    assert start <= NOW.timestamp() - 60 and end >= NOW.timestamp()


def test_a_refusal_or_no_answer_is_named(store):
    refused = ScriptedProphetX(submit=(400, {"message": "invalid_price"}))
    result = desk_with(store, refused).place(request())
    assert result["status"] == "rejected" and "invalid_price" in result["message"]
    silent = ScriptedProphetX(submit=OSError("timed out"))
    result = desk_with(store, silent).place(request(clickId="abc12345-0003", prophetxReadAt="2026-11-01T15:00:00Z"))
    assert result["status"] == "unknown" and "did not answer" in result["message"]


def test_teams_match_by_unabated_id_not_spelling():
    # Cal's 2026-11-03 refusal: ProphetX's "Sam Houston State" is the board's "Sam Houston St".
    sam = {**EVENT, "competitors": [{"id": 21, "name": "Alabama", "side": "home"},
                                    {"id": 22, "name": "Sam Houston State", "side": "away"}]}
    line = {**LINE, "league": "cbb", "awayTeam": "Sam Houston St", "homeTeam": "Alabama", "awayTeamId": "7", "homeTeamId": "8"}
    asked = []

    def read(league, away, home):
        asked.append((league, away, home))
        return {"ok": True, "awayKey": "cbb:7", "homeKey": "cbb:8"}
    assert teams_mismatch(sam, line, read) is None
    assert asked == [("cbb", "Sam Houston State", "Alabama")]
    assert teams_mismatch(sam, {**line, "awayTeamId": "8", "homeTeamId": "7"}, read).endswith("the other way round from the board")


def test_without_strike_ids_a_differently_spelled_game_is_still_found(store):
    spelled = {**EVENT, "competitors": [{"id": 11, "name": "KC Chiefs", "side": "home"},
                                        {"id": 12, "name": "Buffalo Bills", "side": "away"}]}
    client = ScriptedProphetX(events=[spelled])
    keys = team_keys_reader({**TEAM_KEYS, "KC Chiefs": "nfl:2"})
    assert desk_with(store, client, {**QUOTE, "strikeIds": []}, team_keys=keys).place(request())["status"] == "filled"


def test_a_name_nothing_resolves_or_a_runner_that_cannot_answer_refuses(store):
    unknown = desk_with(store, ScriptedProphetX(), team_keys=team_keys_reader({"Kansas City Chiefs": "nfl:2"})).place(request())
    assert unknown["status"] == "refused" and "Buffalo Bills resolves to no nfl team" in unknown["message"]
    down = desk_with(store, ScriptedProphetX(),
                     team_keys=lambda *_names: {"ok": False, "reason": "runner", "detail": "runner timed out"}).place(request())
    assert down["status"] == "refused" and "could not be checked" in down["message"] and "runner timed out" in down["message"]
    found_by_teams = desk_with(store, ScriptedProphetX(), {**QUOTE, "strikeIds": []},
                               team_keys=lambda *_names: {"ok": False, "detail": "runner timed out"}).place(request())
    assert found_by_teams["status"] == "refused" and "runner timed out" in found_by_teams["message"]
    assert store._con.execute("SELECT count(*) FROM manual_orders").fetchone() == (0,)


def test_place_opens_while_the_server_is_blocked_and_off_on_a_bad_key(store):
    from unabated_ticket.bets_service import service

    class Blocked:
        access_key = "ak123456"

        def __init__(self, error):
            self.error = error

        def balance(self):
            raise self.error

    blocked = ProphetXSource(Blocked(ProphetXBlocked("ProphetX's CloudFront refused this server", 5000.0, True)), dict)
    desk, reason = service.open_prophetx_place_desk(store, [blocked])
    assert desk is not None and reason == ""
    for error in (ProphetXError("ProphetX login refused: HTTP 401"),
                  ProphetXBlocked("ProphetX login refused: HTTP 403 session_num_exceed", 5000.0, False)):
        desk, reason = service.open_prophetx_place_desk(store, [ProphetXSource(Blocked(error), dict)])
        assert desk is None and "could not be used" in reason


def test_a_not_found_event_list_reads_as_no_games_not_an_error():
    sent = []
    answers = [(200, {"data": {"access_token": "t1"}}),
               (404, {"error": "data_not_found", "message": "data not found"})]
    client = ProphetXClient("ak", "sk", "https://px.test/partner", transport=scripted_transport(answers, sent),
                            clock=lambda: 1000.0, min_send_interval_sec=0)
    assert client.sport_events([50016501]) == []


def test_a_batch_that_lists_nothing_is_asked_game_by_game_once_an_hour():
    now = [NOW.timestamp()]
    ended = order("o9", "y" * 32, 10, sport_event_id=1700008779)

    class OneEnded(ScriptedReads):
        def sport_events(self, event_ids=None):
            self.event_reads.append(event_ids)
            return [EVENT] if event_ids == [PX_EVENT_ID] else []

    client = OneEnded(orders=[order("o6", "b" * 32, 80, price=-160), ended])
    source = ProphetXSource(client, dict, history_days=31, clock=lambda: now[0])
    records = {record["id"]: record for record in source.fetch()}
    assert records["prophetx:o6"]["betType"] == "moneyline"
    now[0] += 60
    source.fetch()
    assert client.event_reads == [sorted([PX_EVENT_ID, 1700008779]), [PX_EVENT_ID], [1700008779]]


def test_sends_are_spaced_a_second_apart_and_a_rate_limited_read_goes_once_more():
    mono = [0.0]
    slept = []

    def sleep(seconds):
        slept.append(round(seconds, 2))
        mono[0] += seconds

    sent = []
    answers = [(200, {"data": {"access_token": "t1"}}),
               (429, {"error": "rate_limit_reached", "message": "please send only 1 requests/second"}),
               (200, {"data": {"balance": 1}})]
    client = ProphetXClient("ak", "sk", "https://px.test/partner", transport=scripted_transport(answers, sent),
                            clock=lambda: 1000.0, monotonic=lambda: mono[0], sleep=sleep)
    assert client.balance() == (200, {"data": {"balance": 1}})
    assert [path for _, path, _, _ in sent] == ["/auth/login", "/v4/mm/get_balance", "/v4/mm/get_balance"]
    # login at 0, balance 1.1 s later, the 429 waits 2 s, which already spaces the resend.
    assert slept == [1.1, 2.0]


def test_a_rate_limited_order_is_not_sent_again():
    sent = []
    answers = [(200, {"data": {"access_token": "t1"}}), (429, {"error": "rate_limit_reached"})]
    client = ProphetXClient("ak", "sk", "https://px.test/partner", transport=scripted_transport(answers, sent),
                            clock=lambda: 1000.0, min_send_interval_sec=0)
    status, _reply = client.submit_fill_or_kill("c" * 32, -110, 10, "xman-1")
    assert status == 429 and len(sent) == 2
