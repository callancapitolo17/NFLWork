"""Local web app for Wagerzon NFL trifectas/superfectas: fair price, edge,
recommended stake, and a Place button.

    python -m nfl_specials.app        (or nfl_specials/run.sh)
    -> http://127.0.0.1:8096

Routes (loopback only; every verb refuses a foreign Host header, and POSTs
must be application/json, which no other site can send here without a
preflight this server never answers):
    GET  /                       the page (static/index.html)
    GET  /api/board?account=X    latest board, sized at the current settings
                                 and fitted to account X's Wagerzon balance
    POST /api/refresh            start a refresh now (no-op while one is running)
    POST /api/settings           {bankroll, kelly_fraction}
    POST /api/place              {wz_game_id, wz_american, risk, account}
Refreshes run only when asked (Refresh button; user decision 2026-10-02).
Side effects: fecta_quotes APPEND (each refresh), placed_fectas APPEND (each
placement attempt), settings UPSERT — all in config.STATE_DB_PATH. A Place
click submits a REAL Wagerzon bet (wz.place_fecta).
"""
from __future__ import annotations

import json
import logging
import threading
from datetime import datetime, timedelta, timezone
from http.server import BaseHTTPRequestHandler, ThreadingHTTPServer
from urllib.parse import parse_qs, urlparse

from nfl_specials import config, wz
from nfl_specials.board import Board, refresh_board, size_board
from nfl_specials.books import parse_start
from nfl_specials.pricing import decimal_to_american
from nfl_specials.store import Store

log = logging.getLogger("nfl_specials.app")

INDEX_HTML = config.PACKAGE_DIR / "static" / "index.html"
PLACEMENT_HISTORY_DAYS = 7
# A second placement on the same special within this window is refused: a
# slow Wagerzon answer invites a second click, and each click is a real bet.
DUPLICATE_PLACEMENT_WINDOW = timedelta(minutes=2)


def _american(prob: float | None) -> int | None:
    if prob is None or not 0.0 < prob < 1.0:
        return None
    return decimal_to_american(1.0 / prob)


class SpecialsApp:
    """Holds the latest published board and runs refreshes one at a time."""

    def __init__(self, store: Store) -> None:
        self.store = store
        self._lock = threading.Lock()
        self._board: Board | None = None
        self._refreshing = False
        self._last_error: str | None = None
        self._place_lock = threading.Lock()   # one Wagerzon submission at a time
        # Wagerzon account -> 'available' (the budget), read when a refresh
        # starts and ends, lowered after each placement. None = unknown.
        self._balances: dict[str, float | None] = {}

    # --- refresh -----------------------------------------------------------
    def publish(self, board: Board) -> None:
        with self._lock:
            self._board = board

    def start_refresh(self) -> bool:
        """Start a refresh in the background; False if one is already running."""
        with self._lock:
            if self._refreshing:
                return False
            self._refreshing = True
        threading.Thread(target=self._run_refresh, name="refresh", daemon=True).start()
        return True

    def _read_balances(self) -> None:
        balances = wz.available_balances()
        with self._lock:
            self._balances = balances

    def _run_refresh(self) -> None:
        try:
            self._read_balances()
            refresh_board(self.store, self.publish)
            self._read_balances()   # bets may have settled while it priced
            self._last_error = None
        except Exception as exc:
            log.exception("refresh failed")
            self._last_error = f"{type(exc).__name__}: {exc}"[:300]
        finally:
            with self._lock:
                self._refreshing = False

    # --- read --------------------------------------------------------------
    def budget(self, account: str) -> float | None:
        """Wagerzon 'available' for `account`, floored at 0; None = unknown."""
        with self._lock:
            available = self._balances.get(account)
        return None if available is None else max(0.0, available)

    def board_payload(self, account: str | None) -> dict:
        with self._lock:
            board, refreshing = self._board, self._refreshing
        settings = self.store.settings()
        accounts = wz.account_labels()
        account = account if account in accounts else (accounts[0] if accounts else None)
        since = datetime.now(timezone.utc) - timedelta(days=PLACEMENT_HISTORY_DAYS)
        placements = self.store.placements_since(since)
        payload = {
            "refreshing": refreshing,
            "last_error": self._last_error,
            "settings": settings,
            "accounts": accounts,
            "account": account,
            "min_stake": config.WZ_MIN_STAKE,
            "max_stake": config.WZ_MAX_STAKE,
            "available_balance": None,
            "placements": [_placement_json(p) for p in placements],
            "board": None,
        }
        if board is None:
            return payload
        budget = self.budget(account)
        payload["available_balance"] = budget
        # Unknown balance: recommend nothing rather than uncapped Kelly.
        sized = size_board(board, settings["bankroll"], settings["kelly_fraction"],
                           0.0 if budget is None else budget,
                           placed_ids=_placed_ids(placements, account))
        payload["board"] = {
            "started_at": board.started_at.isoformat(),
            "finished_at": board.finished_at.isoformat() if board.finished_at else None,
            "progress_done": board.progress_done,
            "progress_total": board.progress_total,
            "book_status": board.book_status,
            "book_progress": board.book_progress,
            "lines": [_line_json(line, sizing) for line, sizing in zip(board.lines, sized)],
        }
        return payload

    # --- write -------------------------------------------------------------
    def save_settings(self, body: dict) -> dict:
        bankroll, kelly_fraction = body.get("bankroll"), body.get("kelly_fraction")
        if not isinstance(bankroll, (int, float)) or bankroll <= 0:
            raise ValueError(f"bankroll must be a positive number, got {bankroll!r}")
        if not isinstance(kelly_fraction, (int, float)) or not 0 < kelly_fraction <= 1:
            raise ValueError(f"kelly_fraction must be in (0, 1], got {kelly_fraction!r}")
        self.store.save_settings(float(bankroll), float(kelly_fraction))
        return self.store.settings()

    def place(self, body: dict) -> dict:
        """Submit one real Wagerzon bet after re-checking it against the board."""
        wz_game_id, risk = body.get("wz_game_id"), body.get("risk")
        shown_american, account = body.get("wz_american"), body.get("account")
        if not isinstance(wz_game_id, int) or not isinstance(shown_american, int):
            raise ValueError("wz_game_id and wz_american must be whole numbers")
        if not isinstance(risk, (int, float)) or not config.WZ_MIN_STAKE <= risk <= config.WZ_MAX_STAKE:
            raise ValueError(f"risk must be ${config.WZ_MIN_STAKE:.0f}-${config.WZ_MAX_STAKE:.0f} "
                             f"(Wagerzon's limits), got {risk!r}")
        if account not in wz.account_labels():
            raise ValueError(f"unknown Wagerzon account {account!r}")
        with self._lock:
            board = self._board
        line = next((l for l in (board.lines if board else []) if l.special.wz_game_id == wz_game_id), None)
        if line is None or line.game is None:
            raise ValueError(f"special {wz_game_id!r} is not on the current board")
        if line.special.wz_american != shown_american:
            raise ValueError(f"price on the board is {line.special.wz_american:+d}, page showed "
                             f"{shown_american:+d} — refresh the page")
        if parse_start(line.game.game_start_time) <= datetime.now(timezone.utc):
            raise ValueError("game has started")
        settings = self.store.settings()
        budget = self.budget(account)
        if budget is not None and risk > budget:
            raise ValueError(f"risk ${risk:.2f} is more than the ${budget:.2f} Wagerzon shows available")
        sizing = size_board(board, settings["bankroll"], settings["kelly_fraction"],
                            budget)[board.lines.index(line)]

        with self._place_lock:
            # Checked inside the lock so a double click waits for the first
            # placement to land, then sees it.
            since = datetime.now(timezone.utc) - DUPLICATE_PLACEMENT_WINDOW
            recent = self.store.last_placed(account, wz_game_id, since)
            if recent is not None:
                raise ValueError(f"${recent['risk']:.0f} already placed on this special at "
                                 f"{recent['placed_at']:%H:%M:%S} (ticket {recent['ticket_number']}); "
                                 f"wait 2 minutes to add more")
            result = wz.place_fecta(account, line.special, float(risk))
            if result.get("status") == "placed":
                self._lower_balance(account, float(risk), result.get("balance_after"))
            self.store.record_placement({
                "placed_at": datetime.now(timezone.utc), "account": account,
                "wz_game_id": line.special.wz_game_id, "rotation": line.special.rotation,
                "description": line.special.description, "game_start_time": line.game.game_start_time,
                "wz_american": line.special.wz_american, "risk": float(risk),
                "fair_prob": sizing.fair_prob, "ev": sizing.ev,
                "status": result.get("status"), "ticket_number": result.get("ticket_number"),
                "error": result.get("error_msg"),
            })
        return {"status": result.get("status"), "ticket_number": result.get("ticket_number"),
                "error": result.get("error_msg"), "balance_after": result.get("balance_after")}

    def _lower_balance(self, account: str, risk: float, balance_after: float | None) -> None:
        """Keep the budget right between refreshes: Wagerzon's reported
        balance after the bet, else the old balance less the risk."""
        with self._lock:
            before = self._balances.get(account)
            if balance_after is not None:
                self._balances[account] = float(balance_after)
            elif before is not None:
                self._balances[account] = before - risk


def _line_json(line, sizing) -> dict:
    books = {}
    for name, fair in line.book_fairs.items():
        books[name] = {"fair_prob": fair.fair_prob, "fair_american": _american(fair.fair_prob),
                       "sgp_american": decimal_to_american(fair.sgp_decimal) if fair.sgp_decimal else None,
                       "overround": fair.overround, "n_cells": fair.n_cells, "reason": fair.reason}
    fecta = line.fecta
    return {
        "wz_game_id": line.special.wz_game_id,
        "rotation": line.special.rotation,
        "description": line.special.description,
        "wz_american": line.special.wz_american,
        "team": fecta.team if fecta else None,
        "prop_type": fecta.prop_type if fecta else None,
        "legs": [{"kind": leg.kind, "period": leg.period, "line": leg.line, "text": leg.describe(fecta.team)}
                 for leg in fecta.legs] if fecta else [],
        "home": line.game.home if line.game else None,
        "away": line.game.away if line.game else None,
        "game_start_time": line.game.game_start_time if line.game else None,
        "status": line.status,
        "note": line.note,
        "books": books,
        "is_superfecta": line.is_superfecta(),
        "fair_prob": sizing.fair_prob,
        "fair_american": _american(sizing.fair_prob),
        "ev": sizing.ev,
        "kelly_stake": round(sizing.kelly_stake, 2),
        "recommended_stake": sizing.recommended_stake,
        "yields_to": sizing.yields_to,
    }


def _placed_ids(placements: list[dict], account: str | None) -> frozenset[int]:
    """wz_game_ids `account` has a placed bet on, among `placements`."""
    return frozenset(p["wz_game_id"] for p in placements
                     if p["account"] == account and p["status"] == "placed")


def _placement_json(row: dict) -> dict:
    return {key: (value.isoformat() if isinstance(value, datetime) else value) for key, value in row.items()}


def allowed_hosts(port: int) -> tuple[str, ...]:
    return (f"127.0.0.1:{port}", f"localhost:{port}")


def make_handler(app: SpecialsApp, port: int):
    class Handler(BaseHTTPRequestHandler):
        def do_GET(self) -> None:  # noqa: N802 — http.server's name
            if self._refused_foreign_host():
                return
            if self.path == "/":
                body = INDEX_HTML.read_bytes()
                self.send_response(200)
                self.send_header("Content-Type", "text/html; charset=utf-8")
                self.send_header("Content-Length", str(len(body)))
                self.end_headers()
                self.wfile.write(body)
            elif self.path.split("?")[0] == "/api/board":
                query = parse_qs(urlparse(self.path).query)
                self._send_json(200, app.board_payload((query.get("account") or [None])[0]))
            else:
                self._send_json(404, {"error": "not found"})

        def do_POST(self) -> None:  # noqa: N802
            if self._refused_foreign_host():
                return
            if (self.headers.get("Content-Type") or "").split(";")[0].strip() != "application/json":
                self._send_json(415, {"error": "Content-Type must be application/json"})
                return
            try:
                length = int(self.headers.get("Content-Length") or 0)
                body = json.loads(self.rfile.read(length) or b"{}")
                if self.path == "/api/refresh":
                    self._send_json(200, {"started": app.start_refresh()})
                elif self.path == "/api/settings":
                    self._send_json(200, {"settings": app.save_settings(body)})
                elif self.path == "/api/place":
                    self._send_json(200, app.place(body))
                else:
                    self._send_json(404, {"error": "not found"})
            except ValueError as exc:
                self._send_json(400, {"error": str(exc)})
            except Exception as exc:
                log.exception("POST %s failed", self.path)
                self._send_json(500, {"error": f"{type(exc).__name__}: {exc}"[:300]})

        def _refused_foreign_host(self) -> bool:
            # A page whose DNS flips to 127.0.0.1 is same-origin with this
            # server, so CORS alone would not stop it from placing bets.
            if self.headers.get("Host") in allowed_hosts(port):
                return False
            self._send_json(403, {"error": "foreign Host header"})
            return True

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
    app = SpecialsApp(Store())
    server = ThreadingHTTPServer((config.APP_HOST, config.APP_PORT), make_handler(app, config.APP_PORT))
    log.info("NFL fecta pricer on http://%s:%d (click Refresh to price)", config.APP_HOST, config.APP_PORT)
    try:
        server.serve_forever()
    finally:
        server.server_close()


if __name__ == "__main__":
    main()
