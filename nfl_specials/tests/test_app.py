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
from nfl_specials.pricing import BookFair, ScoresFirstShare
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
    free_play_flags = []

    # What Wagerzon holds: a bet draws on the pot it was sent for, unless a
    # test sets charges_cash to have Wagerzon book a free play as cash.
    wagerzon = {"available": 500.0, "free_play": 300.0, "charges_cash": False}

    def fake_place(account, special, risk, use_free_play=False):
        placements.append(risk)
        free_play_flags.append(use_free_play)
        pot = "free_play" if use_free_play and not wagerzon["charges_cash"] else "available"
        wagerzon[pot] -= risk
        return {"status": "placed", "ticket_number": f"T{len(placements)}", "balance_after": None}

    monkeypatch.setattr(app_module.wz, "place_fecta", fake_place)
    monkeypatch.setattr(app_module.wz, "read_account_balances",
                        lambda account: app_module.wz.WzBalances(wagerzon["available"], wagerzon["free_play"]))
    specials_app = SpecialsApp(Store(tmp_path / "state.duckdb"))
    line = FectaLine(SPECIAL, parse_fecta(SPECIAL.description), "priced", game=GAME,
                     book_fairs={"FanDuel": BookFair("FanDuel", 0.081, None, 1.4, 18)})
    specials_app.publish(Board(datetime.now(timezone.utc), datetime.now(timezone.utc), [line], {}))
    specials_app._balances = {"Wagerzon": 500.0}
    specials_app._free_play = {"Wagerzon": 300.0}
    specials_app.placements = placements
    specials_app.free_play_flags = free_play_flags
    specials_app.wagerzon = wagerzon
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


def test_payload_carries_a_draftkings_down_superfectas_books_and_fixed_share(app):
    special = WzSpecial(wz_game_id=2, rotation=2, description="CHARGERS SUPERFECTA (SCR 1ST, 1Q, 1H & GM)",
                        wz_american=1600, week_date=date(2098, 12, 31))
    line = FectaLine(special, parse_fecta(special.description), "priced", game=GAME,
                     book_fairs={"FanDuel": BookFair("FanDuel", 0.12, None, 1.3, 8),
                                 "BetMGM": BookFair("BetMGM", 0.10, None, 1.3, 8)},
                     sf_share=ScoresFirstShare("fixed", 0.776, 0, reason="DraftKings down: error: HTTP 403"))
    app.publish(Board(datetime.now(timezone.utc), datetime.now(timezone.utc), [line], {}))
    payload = app.board_payload("Wagerzon")["board"]["lines"][0]
    assert {key: payload[key] for key in ("is_superfecta", "books", "sf_share", "fair_prob")} == {
        "is_superfecta": True,
        "books": {"FanDuel": {"fair_prob": 0.12, "fair_american": 733, "sgp_american": None, "overround": 1.3,
                              "n_cells": 8, "reason": None},
                  "BetMGM": {"fair_prob": 0.10, "fair_american": 900, "sgp_american": None, "overround": 1.3,
                             "n_cells": 8, "reason": None}},
        "sf_share": {"share": 0.776, "book": "fixed"},
        "fair_prob": pytest.approx(0.10 * 0.776),
    }


def test_a_placed_special_is_topped_up_to_kelly_not_doubled(app):
    app.store.save_settings(35_000.0, 0.25, False)
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


def test_free_play_fits_stakes_to_the_free_play_balance(app):
    app.save_settings({"bankroll": 35_000.0, "kelly_fraction": 0.25, "use_free_play": True})
    app._free_play = {"Wagerzon": 60.0}
    payload = app.board_payload("Wagerzon")
    line = payload["board"]["lines"][0]
    assert payload["settings"]["use_free_play"] is True
    assert (payload["available_balance"], payload["cash_available"], payload["free_play_balance"]) == (60.0, 500.0, 60.0)
    assert 55 <= line["recommended_stake"] <= 60   # Kelly is far above it, so the stake fills the pot
    # fair 0.081 at +1600: $1 of free play is worth 0.081 x 16
    assert line["free_play_value"] == pytest.approx(0.081 * 16)


def test_a_free_play_bet_spends_the_free_play_balance_not_cash(app):
    app.save_settings({"bankroll": 1000.0, "kelly_fraction": 0.25, "use_free_play": True})
    result = app.place({**_bet(100), "free_play": True})
    assert (result["status"], result["free_play"], result["warning"]) == ("placed", True, None)
    assert app.free_play_flags == [True]
    assert (app.budget("Wagerzon", use_free_play=True), app.budget("Wagerzon")) == (200.0, 500.0)
    assert app.store.placements_since(datetime.now(timezone.utc) - timedelta(days=1))[0]["free_play"] is True


def test_a_free_play_wagerzon_books_as_cash_is_recorded_as_cash_and_flagged(app):
    app.save_settings({"bankroll": 1000.0, "kelly_fraction": 0.25, "use_free_play": True})
    app.wagerzon["charges_cash"] = True
    result = app.place({**_bet(100), "free_play": True})
    assert (result["status"], result["free_play"]) == ("placed", False)
    assert result["warning"] == "Wagerzon charged this to cash, not free play"
    assert (app.budget("Wagerzon", use_free_play=True), app.budget("Wagerzon")) == (300.0, 400.0)
    assert app.store.placements_since(datetime.now(timezone.utc) - timedelta(days=1))[0]["free_play"] is False


def test_an_unread_balance_after_a_free_play_lowers_the_pot_and_says_so(app, monkeypatch):
    app.save_settings({"bankroll": 1000.0, "kelly_fraction": 0.25, "use_free_play": True})
    monkeypatch.setattr(app_module.wz, "read_account_balances", lambda account: app_module.wz.WzBalances(None, None))
    result = app.place({**_bet(100), "free_play": True})
    assert (result["free_play"], result["warning"]) == (True, "couldn't confirm which balance Wagerzon charged — check the account")
    assert app.budget("Wagerzon", use_free_play=True) == 200.0


def test_refreshed_balances_fill_each_accounts_two_pots(app, monkeypatch):
    monkeypatch.setattr(app_module.wz, "read_balances", lambda: {
        "Wagerzon": app_module.wz.WzBalances(500.0, 60.0), "WagerzonC": app_module.wz.WzBalances(900.0, None)})
    app._read_balances()
    assert [app.budget(a, f) for a in ("Wagerzon", "WagerzonC") for f in (False, True)] == [500.0, 60.0, 900.0, None]


def test_the_cash_view_still_shows_the_free_play_balance(app):
    payload = app.board_payload("Wagerzon")
    assert (payload["available_balance"], payload["cash_available"], payload["free_play_balance"]) == (500.0, 500.0, 300.0)


def test_free_play_over_its_balance_or_against_the_setting_is_refused(app):
    with pytest.raises(ValueError, match="refresh the page"):
        app.place({**_bet(25), "free_play": True})
    app.save_settings({"bankroll": 1000.0, "kelly_fraction": 0.25, "use_free_play": True})
    with pytest.raises(ValueError, match="refresh the page"):
        app.place(_bet(25))
    app._free_play = {"Wagerzon": 100.0}
    with pytest.raises(ValueError, match="\\$100.00 Wagerzon shows free play"):
        app.place({**_bet(120), "free_play": True})
    app._free_play = {}
    with pytest.raises(ValueError, match="free play balance unknown"):
        app.place({**_bet(25), "free_play": True})
    assert app.placements == []


def test_settings_keep_free_play_when_the_page_omits_it(app):
    app.save_settings({"bankroll": 1000.0, "kelly_fraction": 0.25, "use_free_play": True})
    assert app.save_settings({"bankroll": 2000.0, "kelly_fraction": 0.25})["use_free_play"] is True
    with pytest.raises(ValueError, match="use_free_play"):
        app.save_settings({"bankroll": 2000.0, "kelly_fraction": 0.25, "use_free_play": 1})


def test_the_tailnet_name_is_allowed_only_on_the_tailscale_serve_port():
    hosts = app_module.allowed_hosts(8096, ("mlb-stack.tail87ab2d.ts.net",))
    assert hosts == ("127.0.0.1:8096", "localhost:8096", "mlb-stack.tail87ab2d.ts.net:8443")


def test_without_tailnet_names_only_loopback_is_allowed():
    assert app_module.allowed_hosts(8096, ()) == ("127.0.0.1:8096", "localhost:8096")


@pytest.fixture
def page_port(app):
    server = ThreadingHTTPServer(("127.0.0.1", 0), app_module.make_handler(app, 0, ()))
    port = server.server_address[1]
    server.RequestHandlerClass = app_module.make_handler(app, port, ("mlb-stack.tail87ab2d.ts.net",))
    threading.Thread(target=server.serve_forever, daemon=True).start()
    yield port
    server.shutdown()
    server.server_close()


def _get(port, path, host):
    connection = http.client.HTTPConnection("127.0.0.1", port, timeout=5)
    connection.request("GET", path, headers={"Host": host})
    response = connection.getresponse()
    body = response.read()
    connection.close()
    return response, body


def test_the_page_answers_its_names_and_refuses_any_other_host(page_port):
    port = page_port
    answers = {}
    for host in (f"127.0.0.1:{port}", "mlb-stack.tail87ab2d.ts.net:8443", "mlb-stack.tail87ab2d.ts.net",
                 "evil.example"):
        response, body = _get(port, "/api/board", host)
        answers[host] = (response.status, json.loads(body).get("error"))
    assert answers == {f"127.0.0.1:{port}": (200, None),
                       "mlb-stack.tail87ab2d.ts.net:8443": (200, None),
                       "mlb-stack.tail87ab2d.ts.net": (403, "foreign Host header"),
                       "evil.example": (403, "foreign Host header")}


def test_no_other_site_can_frame_the_page(page_port):
    response, body = _get(page_port, "/", "mlb-stack.tail87ab2d.ts.net:8443")
    assert response.status == 200 and body.startswith(b"<!")
    assert {name: response.getheader(name) for name in app_module.PAGE_SECURITY_HEADERS} == {
        "Content-Security-Policy": "frame-ancestors 'none'; base-uri 'none'",
        "X-Content-Type-Options": "nosniff",
        "Referrer-Policy": "no-referrer"}
