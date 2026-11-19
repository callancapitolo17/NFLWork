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


def leg(leg_id, opposite, market, line=43.5, odds=-110.0, opposite_odds=-110.0):
    return OfferedLeg(leg_id, market, "Game Lines", f"side {leg_id}", CanonicalLeg("FG", "total", None, "over", line),
                      odds, opposite_odds, 100.0, opposite, market, True, {"lineId": leg_id})


TWO_MARKETS = [leg("a", "a2", "m1"), leg("a2", "a", "m1"), leg("b", "b2", "m2"), leg("b2", "b", "m2")]
THREE_MARKETS = TWO_MARKETS + [leg("c", "c2", "m3"), leg("c2", "c", "m3")]
# Two 50% legs that hit together 30% of the time, 10% parlay vig: true factor +20%.
ALL_CELLS = {("a", "b"): 1 / 0.33, ("a", "b2"): 1 / 0.22, ("a2", "b"): 1 / 0.22, ("a2", "b2"): 1 / 0.33}


class FakeVenue:
    """Prices any parlay at a fixed decimal per combination; records requests."""
    name, display_name = "fake", "Fake"

    def __init__(self, cell_prices, legs=TWO_MARKETS, fail_on=(), list_error=None):
        self.cell_prices = cell_prices
        self.legs = legs
        self.fail_on = set(fail_on)
        self.list_error = list_error
        self.requests = []

    def list_games(self):
        return [GAME]

    def list_legs(self, game, fresh=False):
        if self.list_error:
            raise self.list_error
        return self.legs

    def price_parlay(self, game, legs):
        key = tuple(l.leg_id for l in legs)
        self.requests.append(key)
        if key in self.fail_on:
            raise VenueTransportError("HTTP 403")
        price = self.cell_prices.get(key)
        return None if price is None else ParlayQuote((ParlayTier(price, 50.0),), {}, {"key": list(key)})


def test_two_legs_price_all_four_combinations_in_order():
    venue = FakeVenue(ALL_CELLS)
    result = price_at_venue(venue, GAME, ["a", "b"])
    assert result.status == "priced"
    assert venue.requests == [("a", "b"), ("a", "b2"), ("a2", "b"), ("a2", "b2")]
    assert result.partition_decimal_odds == [ALL_CELLS[k] for k in venue.requests]
    assert result.report.correlation_fair == pytest.approx(0.20)
    assert result.report.correlation_raw == pytest.approx(0.33 / (110 / 210) ** 2 - 1)


def test_three_legs_price_eight_combinations():
    cells = {}
    for a in ("a", "a2"):
        for b in ("b", "b2"):
            for c in ("c", "c2"):
                cells[(a, b, c)] = 7.2
    venue = FakeVenue(cells, legs=THREE_MARKETS)
    result = price_at_venue(venue, GAME, ["a", "b", "c"])
    assert len(venue.requests) == 8
    assert result.report.correlation_fair == pytest.approx(0.0)


def test_four_legs_skip_the_partition_with_a_reason():
    legs = THREE_MARKETS + [leg("d", "d2", "m4"), leg("d2", "d", "m4")]
    venue = FakeVenue({("a", "b", "c", "d"): 20.0}, legs=legs)
    result = price_at_venue(venue, GAME, ["a", "b", "c", "d"])
    assert venue.requests == [("a", "b", "c", "d")]
    assert result.report.correlation_fair is None
    assert "only done up to 3 legs" in result.report.fair_unavailable_reason


def test_a_failed_combination_stops_the_partition_and_says_it_failed():
    venue = FakeVenue(ALL_CELLS, fail_on=[("a", "b2")])
    result = price_at_venue(venue, GAME, ["a", "b"])
    assert result.status == "priced"
    assert venue.requests == [("a", "b"), ("a", "b2")]
    assert result.report.fair_unavailable_reason == "pricing the combination side a + side b2 failed: HTTP 403"
    assert result.report.correlation_raw is not None


def test_a_declined_combination_stops_the_partition_and_names_it():
    cells = dict(ALL_CELLS)
    del cells[("a2", "b")]
    venue = FakeVenue(cells)
    result = price_at_venue(venue, GAME, ["a", "b"])
    assert len(venue.requests) == 3
    assert result.report.fair_unavailable_reason == "Fake would not price the combination side a2 + side b"


def test_a_whole_number_line_skips_the_partition():
    legs = [leg("a", "a2", "m1", line=44.0), leg("a2", "a", "m1", line=44.0), *TWO_MARKETS[2:]]
    venue = FakeVenue(ALL_CELLS, legs=legs)
    result = price_at_venue(venue, GAME, ["a", "b"])
    assert venue.requests == [("a", "b")]
    assert "is a whole number, so it can push" in result.report.fair_unavailable_reason


def test_declined_missing_and_failed_reads_keep_the_requested_legs():
    declined = price_at_venue(FakeVenue({}), GAME, ["a", "b"])
    assert (declined.status, declined.detail) == ("declined", "Fake would not price this combination")
    missing = price_at_venue(FakeVenue(ALL_CELLS), GAME, ["a", "zzz"])
    assert missing.status == "missing_legs" and missing.requested_leg_ids == ["a", "zzz"]
    assert missing.detail.startswith("1 of 2 leg(s) no longer offered")
    failed = price_at_venue(FakeVenue(ALL_CELLS, list_error=VenueTransportError("HTTP 500")), GAME, ["a", "b"])
    assert (failed.status, failed.detail, failed.requested_leg_ids) == ("error", "HTTP 500", ["a", "b"])


def test_bad_leg_sets_are_refused_before_any_request():
    venue = FakeVenue(ALL_CELLS)
    with pytest.raises(ValueError, match="same market"):
        price_at_venue(venue, GAME, ["a", "a2"])
    with pytest.raises(ValueError, match=">= 2 legs"):
        price_at_venue(venue, GAME, ["a"])
    with pytest.raises(ValueError, match="picked twice"):
        price_at_venue(venue, GAME, ["a", "a"])
    assert venue.requests == []


# ---------------------------------------------------------------- the web app

@pytest.fixture
def server(tmp_path):
    store = Store(tmp_path / "q.duckdb")
    app = CompareApp({"fake": FakeVenue(ALL_CELLS)}, store)
    httpd = ThreadingHTTPServer(("127.0.0.1", 0), None)
    port = httpd.server_address[1]
    httpd.RequestHandlerClass = make_handler(app, port)
    threading.Thread(target=httpd.serve_forever, daemon=True).start()
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
    raw = response.read()
    return response.status, (json.loads(raw) if response.getheader("Content-Type") == "application/json" else raw)


def test_price_route_appends_one_full_row(server):
    port, store = server
    status, body = request(port, "POST", "/api/price", {"venue": "fake", "game_id": "g1", "leg_ids": ["a", "b"]})
    assert status == 200 and body["store_error"] is None
    assert body["results"][0]["report"]["correlation_fair"] == pytest.approx(0.20)
    with duckdb.connect(str(store.db_path)) as con:
        rows = con.execute("""
            SELECT venue, status, venue_game_id, game_start_time = ?, n_legs, requested_leg_ids_json, legs_json,
                   best_decimal_odds, partition_json, correlation_fair, partition_overround, fair_unavailable_reason
            FROM parlay_quotes""", [GAME.game_start_time]).fetchall()
    assert len(rows) == 1
    venue, status_, game_id, start_matches, n_legs, requested, legs_json, best, partition, fair, overround, reason = rows[0]
    assert (venue, status_, game_id, start_matches, n_legs, json.loads(requested)) == ("fake", "priced", "g1", True, 2, ["a", "b"])
    stored_legs = json.loads(legs_json)
    assert [l["leg_id"] for l in stored_legs] == ["a", "b"]
    assert stored_legs[0]["canonical"] == {"period": "FG", "stat": "total", "subject": None, "side": "over", "line": 43.5}
    assert stored_legs[0]["opposite_leg_id"] == "a2" and stored_legs[0]["venue_ref"] == {"lineId": "a"}
    assert best == pytest.approx(1 / 0.33) and len(json.loads(partition)) == 4
    assert (fair, overround, reason) == (pytest.approx(0.20), pytest.approx(0.10), None)


def test_a_locked_store_still_returns_the_price(server, monkeypatch):
    port, store = server

    def locked(*args, **kwargs):
        raise duckdb.IOException("Could not set lock on file")
    monkeypatch.setattr(store, "append_quote", locked)
    status, body = request(port, "POST", "/api/price", {"venue": "fake", "game_id": "g1", "leg_ids": ["a", "b"]})
    assert status == 200 and body["results"][0]["status"] == "priced"
    assert body["store_error"].startswith("not saved: IOException")


def test_read_routes(server):
    port, _ = server
    status, games = request(port, "GET", "/api/games")
    assert status == 200 and games["by_venue"]["fake"]["games"][0]["game_id"] == "g1"
    status, legs = request(port, "GET", "/api/legs?venue=fake&game_id=g1")
    assert status == 200 and [l["leg_id"] for l in legs["legs"]] == ["a", "a2", "b", "b2"]
    status, page = request(port, "GET", "/")
    assert status == 200 and b"Correlation factor" in page
    assert request(port, "GET", "/nope")[0] == 404


def test_refusals(server):
    port, _ = server
    assert request(port, "GET", "/api/games", host="evil.example:80")[0] == 403
    assert request(port, "POST", "/api/price", {"x": 1}, content_type="text/plain")[0] == 415
    assert request(port, "POST", "/api/price", {"venue": "nope", "game_id": "g1", "leg_ids": []})[0] == 400
    assert request(port, "POST", "/api/price", ["not", "an", "object"])[0] == 400
    assert request(port, "POST", "/api/price", {"venue": "fake", "game_id": "started", "leg_ids": ["a", "b"]})[0] == 400
