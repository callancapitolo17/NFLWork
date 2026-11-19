"""Local web page: build an NFL same-game parlay, price it at each venue, see
how correlated each venue thinks the legs are.

    sgp_compare/run.sh        ->  http://127.0.0.1:8097

Routes (bound to loopback; every verb refuses a Host header other than the
loopback names, so a page whose DNS flips to 127.0.0.1 cannot use it; POSTs
must be application/json):
    GET  /                                   the page (static/index.html)
    GET  /api/games                          upcoming NFL games per venue
    GET  /api/legs?venue=V&game_id=G         the game's legs at venue V
    POST /api/price {venue, game_id, leg_ids}  price the parlay; APPENDS one row
                                             to parlay_quotes (store.py)
    GET  /api/history                        the 50 newest stored quotes
Nothing here places a bet: every venue call is a read or a price request.
"""
from __future__ import annotations

import json
import logging
import os
from http.server import BaseHTTPRequestHandler, ThreadingHTTPServer
from pathlib import Path
from urllib.parse import parse_qs, urlparse

from sgp_compare.legs import Game, OfferedLeg
from sgp_compare.pricing import price_at_venue
from sgp_compare.store import Store
from sgp_compare.venues import build_venues
from sgp_compare.venues.base import Venue, VenueTransportError

log = logging.getLogger("sgp_compare.app")

APP_HOST = "127.0.0.1"
APP_PORT = int(os.environ.get("SGP_COMPARE_PORT", "8097"))
INDEX_HTML = Path(__file__).resolve().parent / "static" / "index.html"
MAX_REQUEST_BYTES = 64 * 1024
PAGE_SECURITY_HEADERS = {
    "Content-Security-Policy": "default-src 'self'; script-src 'unsafe-inline'; style-src 'unsafe-inline'",
    "X-Content-Type-Options": "nosniff",
    "Referrer-Policy": "no-referrer",
}


class CompareApp:
    def __init__(self, venues: dict[str, Venue], store: Store) -> None:
        self.venues = venues
        self.store = store

    def games_payload(self) -> dict:
        per_venue = {}
        for name, venue in self.venues.items():
            try:
                per_venue[name] = {"games": [_game_json(g) for g in venue.list_games()], "error": None}
            except VenueTransportError as exc:
                per_venue[name] = {"games": [], "error": str(exc)}
        return {"venues": {name: v.display_name for name, v in self.venues.items()}, "by_venue": per_venue}

    def legs_payload(self, venue_name: str, game_id: str) -> dict:
        venue, game = self._venue_and_game(venue_name, game_id)
        return {"game": _game_json(game), "legs": [_offered_leg_json(leg) for leg in venue.list_legs(game)]}

    def price(self, body: dict) -> dict:
        venue, game = self._venue_and_game(str(body.get("venue") or ""), str(body.get("game_id") or ""))
        leg_ids = body.get("leg_ids")
        if not isinstance(leg_ids, list) or not all(isinstance(x, str) for x in leg_ids):
            raise ValueError("leg_ids must be a list of strings")
        result = price_at_venue(venue, game, leg_ids)
        self.store.append_quote(game, len(leg_ids), result)
        return {"game": _game_json(game), "results": [result.to_json()]}

    def _venue_and_game(self, venue_name: str, game_id: str) -> tuple[Venue, Game]:
        venue = self.venues.get(venue_name)
        if venue is None:
            raise ValueError(f"unknown venue {venue_name!r}; expected one of {sorted(self.venues)}")
        game = next((g for g in venue.list_games() if g.venue_game_id == game_id), None)
        if game is None:
            raise ValueError(f"game {game_id!r} is not an upcoming NFL game at {venue.display_name}")
        return venue, game


def _game_json(game: Game) -> dict:
    return {"game_id": game.venue_game_id, "label": game.label, "home_team": game.home_team,
            "away_team": game.away_team, "game_start_time": game.game_start_time.isoformat()}


def _offered_leg_json(leg: OfferedLeg) -> dict:
    return {"leg_id": leg.leg_id, "market_label": leg.market_label, "category": leg.category,
            "side_label": leg.side_label, "american_odds": leg.american_odds,
            "opposite_american_odds": leg.opposite_american_odds, "available_stake": leg.available_stake,
            "market_key": leg.market_key, "is_main_line": leg.is_main_line,
            "canonical": leg.canonical.__dict__ if leg.canonical else None}


def allowed_hosts(port: int) -> tuple[str, ...]:
    return (f"127.0.0.1:{port}", f"localhost:{port}")


def make_handler(app: CompareApp, port: int):
    class Handler(BaseHTTPRequestHandler):
        def do_GET(self) -> None:  # noqa: N802 — http.server's name
            if self._refused_foreign_host():
                return
            url = urlparse(self.path)
            query = parse_qs(url.query)
            try:
                if url.path == "/":
                    self._send_page()
                elif url.path == "/api/games":
                    self._send_json(200, app.games_payload())
                elif url.path == "/api/legs":
                    self._send_json(200, app.legs_payload((query.get("venue") or [""])[0],
                                                          (query.get("game_id") or [""])[0]))
                elif url.path == "/api/history":
                    self._send_json(200, {"quotes": app.store.recent_quotes()})
                else:
                    self._send_json(404, {"error": "not found"})
            except ValueError as exc:
                self._send_json(400, {"error": str(exc)})
            except VenueTransportError as exc:
                self._send_json(502, {"error": str(exc)})
            except Exception as exc:
                log.exception("GET %s failed", self.path)
                self._send_json(500, {"error": f"{type(exc).__name__}: {exc}"[:300]})

        def do_POST(self) -> None:  # noqa: N802
            if self._refused_foreign_host():
                return
            if (self.headers.get("Content-Type") or "").split(";")[0].strip() != "application/json":
                self._send_json(415, {"error": "Content-Type must be application/json"})
                return
            try:
                length = int(self.headers.get("Content-Length") or 0)
                if length > MAX_REQUEST_BYTES:
                    raise ValueError(f"request body over {MAX_REQUEST_BYTES} bytes")
                body = json.loads(self.rfile.read(length) or b"{}")
                if self.path == "/api/price":
                    self._send_json(200, app.price(body))
                else:
                    self._send_json(404, {"error": "not found"})
            except ValueError as exc:
                self._send_json(400, {"error": str(exc)})
            except VenueTransportError as exc:
                self._send_json(502, {"error": str(exc)})
            except Exception as exc:
                log.exception("POST %s failed", self.path)
                self._send_json(500, {"error": f"{type(exc).__name__}: {exc}"[:300]})

        def _refused_foreign_host(self) -> bool:
            if self.headers.get("Host") in allowed_hosts(port):
                return False
            self._send_json(403, {"error": "foreign Host header"})
            return True

        def _send_page(self) -> None:
            body = INDEX_HTML.read_bytes()
            self.send_response(200)
            self.send_header("Content-Type", "text/html; charset=utf-8")
            self.send_header("Content-Length", str(len(body)))
            for name, value in PAGE_SECURITY_HEADERS.items():
                self.send_header(name, value)
            self.end_headers()
            self.wfile.write(body)

        def _send_json(self, status: int, payload: dict) -> None:
            body = json.dumps(payload, default=str).encode()
            self.send_response(status)
            self.send_header("Content-Type", "application/json")
            self.send_header("Content-Length", str(len(body)))
            self.end_headers()
            self.wfile.write(body)

        def log_message(self, format: str, *args: object) -> None:  # noqa: A002
            log.debug("%s " + format, self.address_string(), *args)

    return Handler


def main() -> None:
    logging.basicConfig(level=logging.INFO, format="%(asctime)s %(levelname)s %(name)s: %(message)s")
    app = CompareApp(build_venues(), Store())
    server = ThreadingHTTPServer((APP_HOST, APP_PORT), make_handler(app, APP_PORT))
    log.info("SGP compare on http://%s:%d (venues: %s; quotes -> %s)", APP_HOST, APP_PORT,
             ", ".join(app.venues), app.store.db_path)
    try:
        server.serve_forever()
    finally:
        server.server_close()


if __name__ == "__main__":
    main()
