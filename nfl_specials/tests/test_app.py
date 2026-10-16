import http.client
import json
import threading
from datetime import date, datetime, timedelta, timezone
from http.server import ThreadingHTTPServer

import pytest

from nfl_specials import app as app_module
from nfl_specials.app import SpecialsApp
from nfl_specials.board import Board, FectaLine
from nfl_specials.books import BookGame
from nfl_specials.pricing import BookFair
from nfl_specials.special_parser import parse_fecta
from nfl_specials.store import Store
from nfl_specials.wz import WzSpecial

GAME = BookGame("evt", home="SEA", away="LAC", game_start_time="2099-01-01T00:00:00Z")
SPECIAL = WzSpecial(wz_game_id=6008634, rotation=777033, description="CHARGERS TRIFECTA (1Q, 1H & GM)",
                    wz_american=1600, week_date=date(2098, 12, 31))


@pytest.fixture
def app(tmp_path, monkeypatch):
    monkeypatch.setattr(app_module.wz, "account_labels", lambda: ["Wagerzon"])
    placements = []

    def fake_place(account, special, risk):
        placements.append(risk)
        return {"status": "placed", "ticket_number": f"T{len(placements)}", "balance_after": None}

    monkeypatch.setattr(app_module.wz, "place_fecta", fake_place)
    specials_app = SpecialsApp(Store(tmp_path / "state.duckdb"))
    line = FectaLine(SPECIAL, parse_fecta(SPECIAL.description), "priced", game=GAME,
                     book_fairs={"FanDuel": BookFair("FanDuel", 0.081, None, 1.4, 18)})
    specials_app.publish(Board(datetime.now(timezone.utc), datetime.now(timezone.utc), [line], {}))
    specials_app._balances = {"Wagerzon": 500.0}
    specials_app.placements = placements
    return specials_app


def _bet(risk=25):
    return {"wz_game_id": SPECIAL.wz_game_id, "wz_american": 1600, "risk": risk, "account": "Wagerzon"}


def test_a_second_click_on_the_same_special_is_refused(app):
    assert app.place(_bet())["status"] == "placed"
    with pytest.raises(ValueError, match="already placed"):
        app.place(_bet())
    assert app.placements == [25.0]


def test_placement_lowers_the_budget(app):
    app.place(_bet(100))
    assert app.budget("Wagerzon") == 400.0


def test_unknown_balance_recommends_nothing(app):
    app._balances = {}
    payload = app.board_payload("Wagerzon")
    assert payload["available_balance"] is None
    assert all(line["recommended_stake"] == 0 for line in payload["board"]["lines"])


def test_under_the_minimum_or_a_changed_price_is_refused(app):
    with pytest.raises(ValueError, match="limits"):
        app.place(_bet(15))
    with pytest.raises(ValueError, match="limits"):
        app.place(_bet(300))
    with pytest.raises(ValueError, match="refresh the page"):
        app.place({**_bet(), "wz_american": 1500})
    assert app.placements == []


def test_payload_carries_structured_legs_and_book_progress(app):
    app._board.book_progress = {"FanDuel": {"done": 1, "total": 1}}
    payload = app.board_payload("Wagerzon")
    assert payload["board"]["book_progress"] == {"FanDuel": {"done": 1, "total": 1}}
    assert payload["board"]["lines"][0]["legs"][0] == {"kind": "win", "period": "Q1", "line": None,
                                                       "text": "LAC wins Q1"}


def test_a_placed_special_is_topped_up_to_kelly_not_doubled(app):
    app.store.save_settings(35_000.0, 0.25)
    before = app.board_payload("Wagerzon")["board"]["lines"][0]["recommended_stake"]
    assert before > 25
    app.place(_bet(25))
    line = app.board_payload("Wagerzon")["board"]["lines"][0]
    assert line["recommended_stake"] == pytest.approx(before - 25, abs=1)

def test_a_top_up_cannot_pass_the_per_special_maximum(app, monkeypatch):
    app.place(_bet(200))
    monkeypatch.setattr(app_module, "DUPLICATE_PLACEMENT_WINDOW", timedelta(seconds=0))
    with pytest.raises(ValueError, match="at most \\$50 more"):
        app.place(_bet(60))
    assert app.place(_bet(50))["status"] == "placed"
    assert app.placements == [200.0, 50.0]


def test_the_tailnet_name_is_allowed_only_on_the_tailscale_serve_port():
    hosts = app_module.allowed_hosts(8096, ("mlb-stack.tail87ab2d.ts.net",))
    assert hosts == ("127.0.0.1:8096", "localhost:8096", "mlb-stack.tail87ab2d.ts.net:8443")


def test_without_tailnet_names_only_loopback_is_allowed():
    assert app_module.allowed_hosts(8096, ()) == ("127.0.0.1:8096", "localhost:8096")


def test_the_page_answers_its_names_and_refuses_any_other_host(app):
    server = ThreadingHTTPServer(("127.0.0.1", 0), app_module.make_handler(app, 0, ()))
    port = server.server_address[1]
    server.RequestHandlerClass = app_module.make_handler(app, port, ("mlb-stack.tail87ab2d.ts.net",))
    threading.Thread(target=server.serve_forever, daemon=True).start()
    try:
        answers = {}
        for host in (f"127.0.0.1:{port}", "mlb-stack.tail87ab2d.ts.net:8443", "mlb-stack.tail87ab2d.ts.net",
                     "evil.example"):
            connection = http.client.HTTPConnection("127.0.0.1", port, timeout=5)
            connection.request("GET", "/api/board", headers={"Host": host})
            response = connection.getresponse()
            body = json.loads(response.read())
            answers[host] = (response.status, body.get("error"))
            connection.close()
    finally:
        server.shutdown()
        server.server_close()
    assert answers == {f"127.0.0.1:{port}": (200, None),
                       "mlb-stack.tail87ab2d.ts.net:8443": (200, None),
                       "mlb-stack.tail87ab2d.ts.net": (403, "foreign Host header"),
                       "evil.example": (403, "foreign Host header")}
