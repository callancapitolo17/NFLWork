import json
import threading
from datetime import datetime, timedelta, timezone
from http.client import HTTPConnection
from http.server import ThreadingHTTPServer

import duckdb
import pytest

from sgp_compare.app import CompareApp, make_handler
from sgp_compare.legs import CanonicalLeg, Game, OfferedLeg, ParlayQuote, ParlayTier
from sgp_compare.pricing import price_at_venue
from sgp_compare.store import Store
from sgp_compare.venues.base import VenueTransportError

GAME = Game("g1", "Houston Texans", "Indianapolis Colts", datetime.now(timezone.utc) + timedelta(days=1))


def leg(leg_id, opposite, market, odds=-110.0, opposite_odds=-110.0):
    return OfferedLeg(leg_id, "g1", market, "Game Lines", leg_id, CanonicalLeg("FG", "total", None, "over", 1.0),
                      odds, opposite_odds, 100.0, opposite, market, True, {"lineId": leg_id})


class FakeVenue:
    """Prices any parlay at a fixed decimal per combination; records requests."""
    name, display_name = "fake", "Fake"

    def __init__(self, cell_prices=None, fail_cells=False):
        self.cell_prices = cell_prices or {}
        self.fail_cells = fail_cells
        self.requests = []
        self.legs = [leg("a", "a2", "m1"), leg("a2", "a", "m1"), leg("b", "b2", "m2"), leg("b2", "b", "m2")]

    def list_games(self):
        return [GAME]

    def list_legs(self, game, fresh=False):
        return self.legs

    def price_parlay(self, game, legs):
        key = tuple(l.leg_id for l in legs)
        self.requests.append(key)
        if self.fail_cells and key != ("a", "b"):
            raise VenueTransportError("boom")
        price = self.cell_prices.get(key)
        return None if price is None else ParlayQuote((ParlayTier(price, 50.0),), {}, {"key": list(key)})


ALL_CELLS = {("a", "b"): 3.0, ("a", "b2"): 4.5, ("a2", "b"): 4.5, ("a2", "b2"): 3.0}


def test_two_legs_price_all_four_combinations():
    venue = FakeVenue(ALL_CELLS)
    result = price_at_venue(venue, GAME, ["a", "b"])
    assert result.status == "priced"
    assert venue.requests == [("a", "b"), ("a", "b2"), ("a2", "b"), ("a2", "b2")]
    assert result.report.correlation_fair > 0 and result.report.correlation_raw > 0


def test_declined_parlay_and_missing_legs_and_transport_errors():
    assert price_at_venue(FakeVenue({}), GAME, ["a", "b"]).status == "declined"
    assert price_at_venue(FakeVenue(ALL_CELLS), GAME, ["a", "zzz"]).status == "missing_legs"
    failing = price_at_venue(FakeVenue(ALL_CELLS, fail_cells=True), GAME, ["a", "b"])
    assert failing.status == "priced" and failing.report.correlation_fair is None


def test_two_legs_from_one_market_are_refused():
    with pytest.raises(ValueError, match="same market"):
        price_at_venue(FakeVenue(ALL_CELLS), GAME, ["a", "a2"])


@pytest.fixture
def server(tmp_path):
    store = Store(tmp_path / "q.duckdb")
    app = CompareApp({"fake": FakeVenue(ALL_CELLS)}, store)
    httpd = ThreadingHTTPServer(("127.0.0.1", 0), None)
    port = httpd.server_address[1]
    httpd.RequestHandlerClass = make_handler(app, port)
    thread = threading.Thread(target=httpd.serve_forever, daemon=True)
    thread.start()
    yield port, store
    httpd.shutdown()
    httpd.server_close()


def request(port, method, path, body=None, host=None, content_type="application/json"):
    conn = HTTPConnection("127.0.0.1", port, timeout=5)
    headers = {"Host": host or f"127.0.0.1:{port}"}
    if body is not None:
        headers["Content-Type"] = content_type
    conn.request(method, path, body=json.dumps(body) if body is not None else None, headers=headers)
    response = conn.getresponse()
    return response.status, json.loads(response.read() or b"{}")


def test_price_route_appends_one_row(server):
    port, store = server
    status, body = request(port, "POST", "/api/price", {"venue": "fake", "game_id": "g1", "leg_ids": ["a", "b"]})
    assert status == 200 and body["results"][0]["status"] == "priced"
    with duckdb.connect(str(store.db_path)) as con:
        rows = con.execute("SELECT venue, n_legs, correlation_raw FROM parlay_quotes").fetchall()
    assert len(rows) == 1 and rows[0][:2] == ("fake", 2)


def test_foreign_host_and_non_json_posts_are_refused(server):
    port, _ = server
    assert request(port, "GET", "/api/games", host="evil.example:80")[0] == 403
    assert request(port, "POST", "/api/price", {"x": 1}, content_type="text/plain")[0] == 415
    assert request(port, "POST", "/api/price", {"venue": "nope", "game_id": "g1", "leg_ids": []})[0] == 400
