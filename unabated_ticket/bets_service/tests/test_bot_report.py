"""GET /bot.json (bot_report.py): the Bet Tracker's Bot tab, read from taker.duckdb."""
import json
import threading
import urllib.error
import urllib.request
from datetime import datetime, timedelta, timezone
from http.server import ThreadingHTTPServer

import pytest

from unabated_ticket.bets_service import bot_report, config, service
from unabated_ticket.bets_service.store import BetsStore
from unabated_ticket.bets_service.taker_store import TakerStore

NOW = datetime(2026, 11, 4, 20, 0, tzinfo=timezone.utc)
KALSHI_ID = "ubot-" + "a" * 32
KALSHI_OLD_ID = "ubot-" + "b" * 32
NOVIG_ID = "0b3f5f20-9c1b-4d6e-8f3a-5c2d1e7b4a90"


def kalshi_decision(reason: str, row_key: str, order: dict | None = None, venue: str | None = "kalshi",
                    live_edge: float = 0.05) -> dict:
    decision = {"rowKey": row_key, "league": "NFL", "eventId": 7, "eventStart": "2026-11-05T01:00:00.000Z",
                "decidedAt": NOW.isoformat(), "reason": reason, "ticker": "KX-1", "side": "yes", "unabatedPrice": -110,
                "liveAllIn": 0.52, "liveEdge": live_edge, "contracts": 10, "eventName": "Bills @ Dolphins",
                "period": "FG", "betType": "spread", "sideLabel": "Bills", "points": -3.5}
    if venue is not None:
        decision["venue"] = venue
    if order is not None:
        decision["order"] = order
    return decision


def log(store: TakerStore, decision: dict, decided_at: datetime) -> None:
    store.append_decisions([{"decision": decision, "decided_at": decided_at, "event_start": NOW + timedelta(hours=5)}],
                           decided_at)


@pytest.fixture
def store(tmp_path):
    taker_store = TakerStore(tmp_path / "taker.duckdb")
    # A filled Kalshi order today, its decision in the log.
    taker_store.insert_order({"client_order_id": KALSHI_ID, "placed_at": NOW - timedelta(hours=1), "row_key": "k1",
                              "event_id": 7, "ticker": "KX-1", "contract_side": "yes", "ask_price": 0.5,
                              "all_in_price": 0.52, "count_requested": 10, "event_start": NOW + timedelta(hours=5),
                              "open_cost": 5.2, "subaccount": 1})
    taker_store.update_order(KALSHI_ID, {"status": "partial", "count_filled": 4.0, "cost_filled": 2.08})
    log(taker_store, kalshi_decision("bet", "k1", {"clientOrderId": KALSHI_ID, "status": "partial"}, live_edge=0.06),
        NOW - timedelta(hours=1))
    # An older, settled Kalshi order (a week ago), and an old decision with no venue (before the Novig bot).
    taker_store.insert_order({"client_order_id": KALSHI_OLD_ID, "placed_at": NOW - timedelta(days=7), "row_key": "k0",
                              "event_id": 6, "ticker": "KX-0", "contract_side": "no", "ask_price": 0.4,
                              "all_in_price": 0.42, "count_requested": 5, "event_start": NOW - timedelta(days=6),
                              "open_cost": 2.1, "subaccount": 1})
    taker_store.update_order(KALSHI_OLD_ID, {"status": "filled", "count_filled": 5.0, "cost_filled": 2.1,
                                             "settled_at": NOW - timedelta(days=6), "open_cost": 0.0})
    log(taker_store, kalshi_decision("edge_short", "k0", venue=None), NOW - timedelta(days=7))
    # Today: two lines short of edge, one refused at the ceiling.
    log(taker_store, kalshi_decision("edge_short", "k2"), NOW - timedelta(minutes=30))
    log(taker_store, kalshi_decision("edge_short", "k3"), NOW - timedelta(minutes=20))
    log(taker_store, kalshi_decision("edge_short", "k3"), NOW - timedelta(minutes=10))  # same line again: one line
    log(taker_store, kalshi_decision("bet", "k4", {"clientOrderId": "ubot-" + "c" * 32, "status": "ceiling"}),
        NOW - timedelta(minutes=5))
    # A Novig order with no fills read yet.
    taker_store.insert_novig_order({
        "client_id": NOVIG_ID, "placed_at": NOW - timedelta(minutes=2), "row_key": "n1", "event_id": 8,
        "event_start": NOW + timedelta(days=1), "outcome_id": "o", "market_id": "m", "price": 0.465,
        "qty_requested": 2150, "status": "sending", "open_cost": 9.9975, "league": "cfb", "bet_type": "spread",
        "period": "FG", "side": "home", "points": -35.5, "away_team": None, "home_team": None,
        "away_team_id": None, "home_team_id": None})
    log(taker_store, kalshi_decision("not_best", "n2", venue="novig"), NOW - timedelta(minutes=3))
    yield taker_store
    taker_store.close()


def test_activity_counts_orders_and_lines_in_the_window(store):
    activity = store.bot_activity(NOW - timedelta(days=1), NOW)
    kalshi = activity["kalshi"]
    assert (kalshi["orders"], kalshi["filledOrders"]) == (1, 1)
    assert kalshi["fillRate"] == pytest.approx(0.4)
    assert kalshi["staked"] == pytest.approx(2.08)
    assert kalshi["avgEdgeAtFill"] == pytest.approx(0.06)
    reasons = {(row["reason"], row["orderStatus"]): row["lines"] for row in kalshi["reasons"]}
    assert reasons == {("edge_short", None): 2, ("bet", "partial"): 1, ("bet", "ceiling"): 1}
    novig = activity["novig"]
    # No fills read yet: the fill rate is unknown, never 0%.
    assert (novig["orders"], novig["fillRate"], novig["staked"]) == (1, None, 0.0)
    assert novig["reasons"] == [{"reason": "not_best", "orderStatus": None, "lines": 1}]


def test_all_time_takes_old_decisions_with_no_venue_as_kalshi(store):
    kalshi = store.bot_activity(None, None)["kalshi"]
    assert kalshi["orders"] == 2
    assert kalshi["fillRate"] == pytest.approx(9 / 15)
    assert {(row["reason"], row["orderStatus"]): row["lines"] for row in kalshi["reasons"]}[("edge_short", None)] == 3


def test_open_orders_and_recent_orders_name_the_line(store):
    assert store.bot_open_orders()["kalshi"]["openOrders"] == 1
    assert store.bot_open_orders()["novig"]["openOrders"] == 1
    recent = store.recent_bot_orders(10)
    assert [order["orderId"] for order in recent] == [NOVIG_ID, KALSHI_ID, KALSHI_OLD_ID]
    kalshi = recent[1]
    assert (kalshi["eventName"], kalshi["sideLabel"], kalshi["points"], kalshi["liveEdge"]) == (
        "Bills @ Dolphins", "Bills", -3.5, 0.06)
    assert recent[0]["eventName"] is None  # its decision is not in the log


def test_desk_status_off_stopped_running(tmp_path):
    stop = tmp_path / "taker.stop"
    assert bot_report.desk_status(False, "no key", (stop,)) == {"status": "off", "detail": "no key"}
    assert bot_report.desk_status(True, "", (stop,))["status"] == "running"
    stop.write_text("")
    assert bot_report.desk_status(True, "", (stop,)) == {"status": "stopped", "detail": f"stop file {stop} is present"}


def test_runner_pass_reads_the_last_pass_or_says_why_not():
    body = {"desk": "armed", "novigDesk": "not asked yet", "rules": {"novigKellyBankroll": 2000},
            "pass": {"finishedAt": "2026-11-04T19:59:50.000Z", "skipped": None, "novigSkipped": "the Novig poll is stale",
                     "rows": 40}}
    read = bot_report.runner_pass(lambda: (200, json.dumps(body).encode()))
    assert (read["finishedAt"], read["novigSkipped"], read["novigKellyBankroll"], read["error"]) == (
        "2026-11-04T19:59:50.000Z", "the Novig poll is stale", 2000, None)
    down = bot_report.runner_pass(lambda: (502, json.dumps({"error": "server runner unreachable"}).encode()))
    assert (down["finishedAt"], down["error"]) == (None, "server runner unreachable")


def test_parse_window():
    start, end = bot_report.parse_window({"start": ["2026-11-04T08:00:00Z"], "end": ["2026-11-05T08:00:00Z"]})
    assert end - start == timedelta(days=1)
    assert bot_report.parse_window({}) == (None, None)
    assert "timezone" in bot_report.parse_window({"start": ["2026-11-04T08:00:00"]})
    assert "after" in bot_report.parse_window({"start": ["2026-11-05T08:00:00Z"], "end": ["2026-11-04T08:00:00Z"]})


def test_get_bot_json_over_http(store, tmp_path, monkeypatch):
    monkeypatch.setattr(config, "TAKER_STOP_PATH", tmp_path / "taker.stop")
    monkeypatch.setattr(config, "NOVIG_TAKER_STOP_PATH", tmp_path / "novig_taker.stop")
    bets_store = BetsStore(tmp_path / "bets.duckdb", 7)
    handler = service.make_handler(bets_store, 0.0, ["kalshi"], runner_url="http://127.0.0.1:9", runner_timeout_sec=0.5,
                                   taker_store=store, novig_taker_desk_off_reason="no bot Novig key is configured")
    server = ThreadingHTTPServer(("127.0.0.1", 0), handler)
    threading.Thread(target=server.serve_forever, daemon=True).start()
    base = f"http://127.0.0.1:{server.server_address[1]}"
    try:
        with urllib.request.urlopen(f"{base}/bot.json?start=2026-11-03T20:00:00Z", timeout=5) as response:
            body = json.loads(response.read())
        with pytest.raises(urllib.error.HTTPError) as refused:
            urllib.request.urlopen(f"{base}/bot.json?start=yesterday", timeout=5)
    finally:
        server.shutdown()
        server.server_close()
        bets_store.close()
    assert refused.value.code == 400
    kalshi, novig = body["bots"]
    assert (kalshi["status"], novig["status"], novig["detail"]) == ("off", "off", "no bot Novig key is configured")
    assert kalshi["capDollars"] == 4500 and novig["capDollars"] == 850
    assert kalshi["openCost"] == pytest.approx(5.2)
    assert kalshi["activity"]["orders"] == 1
    assert body["runner"]["error"]  # no runner on port 9
    assert body["recentOrders"][0]["placedAt"].endswith("Z")
