from datetime import date, datetime, timezone

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
    with pytest.raises(ValueError, match="at least"):
        app.place(_bet(15))
    with pytest.raises(ValueError, match="refresh the page"):
        app.place({**_bet(), "wz_american": 1500})
    assert app.placements == []
