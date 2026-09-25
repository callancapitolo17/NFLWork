"""Entry point: the local bets service the Unabated Ticket panel polls (#114).

    python -m unabated_ticket.bets_service.service      (or bets_service/run.sh)

Inputs:  each registered Source (sources/kalshi.py; sources/betonline.py when
         its cookie file exists; sources/novig.py when its token file exists;
         sources/bfa.py and sources/wagerzon.py when their logins are configured;
         sources/polymarket_us.py when its API key is configured)
         on its own poll_sec.
Outputs: HTTP on 127.0.0.1:8094 (loopback only, no auth):
           GET /bets.json[?days=N]  {generatedAt, sources: {name: {fetchedAt, ok,
                                    error, count}}, bets: [records open + settled
                                    within N days (default RETENTION_DAYS=30)],
                                    crosswalk: [team_crosswalk rows, newest first],
                                    pins: [bet_pins rows, newest first]}
           GET /health              {ok, generatedAt, uptimeSec, sources}
           POST /crosswalk.json     body {rows: [{venue, league, venueTeamKey,
                                    unabatedTeamId, venueTeamName?, unabatedTeamName?,
                                    learnedFrom?}]} -> {ok, learned, conflicts, crosswalk}
                                    (#118 step 4: the panel's lessons from id joins;
                                    Content-Type must be application/json — a web page
                                    cannot send that cross-origin without a preflight
                                    this server never answers, so no site can write here)
           DELETE /crosswalk.json   -> {ok, cleared, crosswalk: []}
           POST /pins.json          body {pin: {betId, venue, league, eventId, eventStart?,
                                    awayTeamId?, homeTeamId?, awayTeamName?, homeTeamName?},
                                    crosswalk: [at most 2 rows, the crosswalk shape]}
                                    -> {ok, pins, crosswalk}; 404 when no bet has that id,
                                    400 when the pin's venue is not the bet's
                                    (Cal's manual attach: the pin plus the team names it
                                    teaches, which REPLACE a held key)
           DELETE /pins.json?betId= -> {ok, removedPin, removedRows, pins, crosswalk}
                                    (Undo: the pin and the rows it taught)
         Every verb refuses a request whose Host header is not the loopback
         name the service is serving on (403): a page at evil.example whose
         DNS flips to 127.0.0.1 is SAME-ORIGIN with this server, so neither
         CORS nor the JSON Content-Type guard applies to it.
Side effects: UPSERTs records into bets.duckdb::bets and APPENDs a row to
bets.duckdb::source_runs per poll (see store.py); INSERTs into / DELETEs from
bets.duckdb::team_crosswalk on the two crosswalk routes; UPSERTs / DELETEs
bets.duckdb::bet_pins and the crosswalk rows a pin taught on the two pin
routes; rotating log at
bets_service.log. A poll that raises writes a failed source_runs row and
leaves `bets` untouched — a dark source never blanks the list.
"""
import json
import logging
import signal
import threading
import time
from datetime import datetime, timezone
from http.server import BaseHTTPRequestHandler, ThreadingHTTPServer
from urllib.parse import parse_qs, urlparse

from unabated_ticket.bets_service import config
from unabated_ticket.bets_service.log_setup import setup_logging
from unabated_ticket.bets_service.sources import Source
from unabated_ticket.bets_service.sources.betonline import source_if_configured as betonline_source_if_configured
from unabated_ticket.bets_service.sources.bfa import source_if_configured as bfa_source_if_configured
from unabated_ticket.bets_service.sources.kalshi import KalshiSource
from unabated_ticket.bets_service.sources.novig import source_if_connected as novig_source_if_connected
from unabated_ticket.bets_service.sources.polymarket_us import source_if_configured as polymarket_us_source_if_configured
from unabated_ticket.bets_service.sources.wagerzon import source_if_configured as wagerzon_source_if_configured
from unabated_ticket.bets_service.store import BetsStore

log = logging.getLogger(__name__)

POLL_TICK_SEC = 1.0
MAX_DAYS = 3650
# A crosswalk POST is a few rows per poll; anything near this is not the panel.
MAX_BODY_BYTES = 1024 * 1024
MAX_CROSSWALK_ROWS_PER_POST = 1000
# The only Host headers a request can carry (DNS rebinding sends its own name).
LOOPBACK_HOST_NAMES = ("127.0.0.1", "localhost")
HTTP_DEFAULT_PORT = 80
CROSSWALK_REQUIRED_FIELDS = ("venue", "league", "venueTeamKey", "unabatedTeamId")
CROSSWALK_OPTIONAL_FIELDS = ("venueTeamName", "unabatedTeamName", "learnedFrom")
# A bet names at most two teams, so one attach teaches at most two rows.
MAX_CROSSWALK_ROWS_PER_PIN = 2
PIN_REQUIRED_FIELDS = ("betId", "venue", "league", "eventId")
PIN_OPTIONAL_ID_FIELDS = ("awayTeamId", "homeTeamId")
PIN_OPTIONAL_TEXT_FIELDS = ("eventStart", "awayTeamName", "homeTeamName")


def _now() -> datetime:
    return datetime.now(timezone.utc)


def _iso(value: datetime) -> str:
    return value.strftime("%Y-%m-%dT%H:%M:%SZ")


def allowed_hosts(port: int) -> tuple[str, ...]:
    """The Host values a request to this port can carry. A browser omits the
    scheme's default port, so on 80 the bare names are valid too."""
    with_port = tuple(f"{name}:{port}" for name in LOOPBACK_HOST_NAMES)
    return with_port + LOOPBACK_HOST_NAMES if port == HTTP_DEFAULT_PORT else with_port


def host_allowed(host_header: str | None, port: int) -> bool:
    """Whether a request's Host names this loopback service. A browser sends
    the name the page was loaded from, so a rebound evil.example does not
    match — the one check that survives the attacker being same-origin."""
    return (host_header or "").strip().lower() in allowed_hosts(port)


def run_source_once(source: Source, store: BetsStore) -> bool:
    """One poll of one source: upsert its records and log the run. Never raises."""
    started_at = _now()
    try:
        records = source.fetch()
        written = store.upsert_bets(records, started_at)
        if written:  # unchanged records are not rewritten, so this is what moved
            log.info("source %s: %d of %d record(s) changed", source.name, written, len(records))
        store.log_source_run(source.name, started_at, _now(), True, None, len(records))
        return True
    except Exception as error:  # a source failure is a failed run, not a dead service
        log.warning("source %s failed: %s: %s", source.name, type(error).__name__, error)
        store.log_source_run(source.name, started_at, _now(), False,
                             f"{type(error).__name__}: {error}", 0)
        return False


def poll_loop(sources: list[Source], store: BetsStore, stop: threading.Event) -> None:
    """Run each source on its own cadence until `stop`. A store write that
    raises (disk full, a locked file) is logged and retried next poll — the
    thread must outlive it, or the HTTP side would serve ageing data with no
    failed run to show for it."""
    next_due = {source.name: 0.0 for source in sources}
    while not stop.is_set():
        now = time.monotonic()
        for source in sources:
            if now >= next_due[source.name]:
                try:
                    run_source_once(source, store)
                except Exception:  # noqa: BLE001 — the poll thread is the service's heartbeat
                    log.exception("store write failed for source %s; retrying next poll", source.name)
                next_due[source.name] = time.monotonic() + source.poll_sec
        stop.wait(POLL_TICK_SEC)


# A registered source that has not finished a poll yet (the first Kalshi poll
# takes ~1-2 min: one throttled GET per market and per event). Without this
# entry the panel would read "no source configured" — the text for venues
# with no source at all (#115-#117).
NO_POLL_YET = {"fetchedAt": None, "ok": False, "error": "no completed poll yet", "count": 0}


def source_status(store: BetsStore, source_names: list[str]) -> dict[str, dict]:
    status = store.source_status()
    for name in source_names:
        status.setdefault(name, dict(NO_POLL_YET))
    return status


def bets_payload(store: BetsStore, days: int, source_names: list[str] = ()) -> dict:
    now = _now()
    return {"generatedAt": _iso(now), "sources": source_status(store, list(source_names)),
            "bets": store.load_bets(days, now), "crosswalk": store.load_crosswalk(), "pins": store.load_pins()}


def _non_empty_string(value: object) -> bool:
    return isinstance(value, str) and value != ""


def _id_as_string(value: object) -> object:
    """Unabated's ids are numeric; the store keeps them as strings. Anything
    else is returned as it came, for the caller's type check to name."""
    if isinstance(value, int) and not isinstance(value, bool):
        return str(value)
    return value


def _is_iso_timestamp(value: str) -> bool:
    try:
        datetime.fromisoformat(value.replace("Z", "+00:00"))
    except ValueError:
        return False
    return True


def validate_pin_request(body: object) -> tuple[dict, list[dict]] | str:
    """(pin, crosswalk rows) of a POST /pins.json body, or an error message.
    Every crosswalk row must be the pin's own venue and league: an attach
    teaches only what that bet's venue calls that game's teams."""
    if not isinstance(body, dict) or not isinstance(body.get("pin"), dict):
        return "body must be an object with a `pin` object"
    raw = body["pin"]
    pin: dict = {}
    for field in PIN_REQUIRED_FIELDS:
        value = _id_as_string(raw.get(field)) if field == "eventId" else raw.get(field)
        if not _non_empty_string(value):
            return f"pin.{field} must be a non-empty string"
        pin[field] = value
    for field in PIN_OPTIONAL_ID_FIELDS:
        value = _id_as_string(raw.get(field))
        if value is not None and not _non_empty_string(value):
            return f"pin.{field} must be a non-empty string, a number or null"
        pin[field] = value
    for field in PIN_OPTIONAL_TEXT_FIELDS:
        value = raw.get(field)
        if value is not None and not isinstance(value, str):
            return f"pin.{field} must be a string or null"
        pin[field] = value
    if pin["eventStart"] is not None and not _is_iso_timestamp(pin["eventStart"]):
        return f"pin.eventStart must be an ISO timestamp or null, got {pin['eventStart']!r}"
    raw_rows = body.get("crosswalk", [])
    if not isinstance(raw_rows, list):
        return "crosswalk must be an array"
    if len(raw_rows) > MAX_CROSSWALK_ROWS_PER_PIN:
        return f"at most {MAX_CROSSWALK_ROWS_PER_PIN} crosswalk rows per pin, got {len(raw_rows)}"
    rows = validate_crosswalk_rows({"rows": raw_rows})
    if isinstance(rows, str):
        return f"crosswalk: {rows}"
    for index, row in enumerate(rows):
        if (row["venue"], row["league"]) != (pin["venue"], pin["league"]):
            return (f"crosswalk[{index}] is {row['venue']}/{row['league']}, "
                    f"the pin is {pin['venue']}/{pin['league']}")
    return pin, rows


def validate_crosswalk_rows(body: object) -> list[dict] | str:
    """The rows of a POST /crosswalk.json body in the store's shape, or an
    error message naming the first bad row. `unabatedTeamId` may arrive as an
    int (Unabated's ids are numeric) and is stored as its string."""
    if not isinstance(body, dict) or not isinstance(body.get("rows"), list):
        return "body must be an object with a `rows` array"
    raw_rows = body["rows"]
    if len(raw_rows) > MAX_CROSSWALK_ROWS_PER_POST:
        return f"at most {MAX_CROSSWALK_ROWS_PER_POST} rows per request, got {len(raw_rows)}"
    rows: list[dict] = []
    for index, raw in enumerate(raw_rows):
        if not isinstance(raw, dict):
            return f"rows[{index}] must be an object"
        row = {"unabatedTeamId": _id_as_string(raw.get("unabatedTeamId"))}
        for field in CROSSWALK_REQUIRED_FIELDS:
            value = row.get(field, raw.get(field))
            if not _non_empty_string(value):
                return f"rows[{index}].{field} must be a non-empty string"
            row[field] = value
        for field in CROSSWALK_OPTIONAL_FIELDS:
            value = raw.get(field)
            if value is not None and not isinstance(value, str):
                return f"rows[{index}].{field} must be a string or null"
            row[field] = value
        rows.append(row)
    return rows


def health_payload(store: BetsStore, started_at: float, source_names: list[str] = ()) -> dict:
    return {"ok": True, "generatedAt": _iso(_now()), "uptimeSec": int(time.monotonic() - started_at),
            "sources": source_status(store, list(source_names))}


def parse_days(query: str) -> int | str:
    """`days` from a query string, or an error message."""
    raw = parse_qs(query).get("days")
    if not raw:
        return config.RETENTION_DAYS
    try:
        days = int(raw[0])
    except ValueError:
        return f"days must be an integer, got {raw[0]!r}"
    if not 0 <= days <= MAX_DAYS:
        return f"days must be between 0 and {MAX_DAYS}, got {days}"
    return days


def make_handler(store: BetsStore, started_at: float,
                 source_names: list[str] = ()) -> type[BaseHTTPRequestHandler]:
    names = list(source_names)

    class BetsHandler(BaseHTTPRequestHandler):
        def do_GET(self) -> None:  # noqa: N802 — http.server's name
            if self._refused_foreign_host():
                return
            url = urlparse(self.path)
            if url.path == "/health":
                self._send_json(200, health_payload(store, started_at, names))
                return
            if url.path == "/bets.json":
                days = parse_days(url.query)
                if isinstance(days, str):
                    self._send_json(400, {"error": days})
                    return
                self._send_json(200, bets_payload(store, days, names))
                return
            self._send_json(404, {"error": f"no route for {url.path}"})

        def do_POST(self) -> None:  # noqa: N802 — http.server's name
            if self._refused_foreign_host():
                return
            url = urlparse(self.path)
            if url.path not in ("/crosswalk.json", "/pins.json"):
                self._send_json(404, {"error": f"no route for POST {url.path}"})
                return
            body = self._read_json_body()
            if isinstance(body, tuple):
                self._send_json(*body)
                return
            if url.path == "/pins.json":
                self._post_pin(body)
                return
            rows = validate_crosswalk_rows(body)
            if isinstance(rows, str):
                self._send_json(400, {"error": rows})
                return
            result = store.learn_crosswalk(rows, _now())
            self._send_json(200, {"ok": True, "learned": result["learned"], "conflicts": result["conflicts"],
                                  "crosswalk": store.load_crosswalk()})

        def _post_pin(self, body: object) -> None:
            request = validate_pin_request(body)
            if isinstance(request, str):
                self._send_json(400, {"error": request})
                return
            pin, rows = request
            venue = store.bet_venue(pin["betId"])
            if venue is None:
                self._send_json(404, {"error": f"no bet with id {pin['betId']!r}"})
                return
            # The rows carry the pin's venue, so the pin must carry the bet's.
            if venue != pin["venue"]:
                self._send_json(400, {"error": f"bet {pin['betId']!r} is {venue}, the pin says {pin['venue']}"})
                return
            store.pin_bet(pin, rows, _now())
            self._send_json(200, {"ok": True, "pins": store.load_pins(), "crosswalk": store.load_crosswalk()})

        def do_DELETE(self) -> None:  # noqa: N802 — http.server's name
            if self._refused_foreign_host():
                return
            url = urlparse(self.path)
            if url.path == "/pins.json":
                self._delete_pin(url.query)
                return
            if url.path != "/crosswalk.json":
                self._send_json(404, {"error": f"no route for DELETE {url.path}"})
                return
            cleared = store.clear_crosswalk()
            self._send_json(200, {"ok": True, "cleared": cleared, "crosswalk": []})

        def _delete_pin(self, query: str) -> None:
            bet_ids = parse_qs(query).get("betId")
            if not bet_ids or not _non_empty_string(bet_ids[0]):
                self._send_json(400, {"error": "betId query parameter required"})
                return
            result = store.unpin_bet(bet_ids[0])
            self._send_json(200, {"ok": True, **result, "pins": store.load_pins(), "crosswalk": store.load_crosswalk()})

        # DNS-rebinding guard, first thing on every verb: a page whose name
        # resolves to 127.0.0.1 is same-origin with this server, so nothing
        # else here stops it reading every bet or clearing the crosswalk. Its
        # Host is still its own name. The port comes from the socket, not
        # config, so a service on a non-default port guards itself.
        def _refused_foreign_host(self) -> bool:
            port = self.server.server_address[1]
            host = self.headers.get("Host")
            if host_allowed(host, port):
                return False
            self._send_json(403, {"error": f"Host must be one of {list(allowed_hosts(port))}, got {host!r}"})
            return True

        # The parsed JSON body, or a (status, error payload) tuple to send. The
        # Content-Type check is the write guard: without CORS headers here a
        # browser page can only reach this server with a "simple" request
        # (form or text/plain), never application/json, and the extension
        # page is exempt from CORS for its host permission.
        def _read_json_body(self) -> object | tuple[int, dict]:
            content_type = self.headers.get("Content-Type", "")
            if not content_type.split(";")[0].strip().lower() == "application/json":
                return 415, {"error": "Content-Type must be application/json"}
            try:
                length = int(self.headers.get("Content-Length", ""))
            except ValueError:
                return 411, {"error": "Content-Length required"}
            if length < 0 or length > MAX_BODY_BYTES:
                return 413, {"error": f"body must be at most {MAX_BODY_BYTES} bytes, got {length}"}
            try:
                return json.loads(self.rfile.read(length))
            except ValueError as error:
                return 400, {"error": f"body is not JSON: {error}"}

        def _send_json(self, status: int, payload: dict) -> None:
            body = json.dumps(payload).encode()
            self.send_response(status)
            self.send_header("Content-Type", "application/json")
            self.send_header("Content-Length", str(len(body)))
            self.send_header("Cache-Control", "no-store")
            self.end_headers()
            self.wfile.write(body)

        def log_message(self, format: str, *args: object) -> None:  # noqa: A002
            log.debug("http %s", format % args)

    return BetsHandler


def serve(sources: list[Source], store: BetsStore, host: str, port: int) -> None:
    """Run the poll thread and the HTTP server until SIGINT/SIGTERM."""
    stop = threading.Event()
    started_at = time.monotonic()
    source_names = [source.name for source in sources]
    server = ThreadingHTTPServer((host, port), make_handler(store, started_at, source_names))
    server.daemon_threads = True
    poller = threading.Thread(target=poll_loop, args=(sources, store, stop), name="poll", daemon=True)

    def request_stop(_signum, _frame) -> None:
        log.info("stopping")
        stop.set()
        threading.Thread(target=server.shutdown, daemon=True).start()

    signal.signal(signal.SIGINT, request_stop)
    signal.signal(signal.SIGTERM, request_stop)
    poller.start()
    log.info("bets service on http://%s:%d (sources: %s)", host, port, ", ".join(source_names))
    try:
        server.serve_forever()
    finally:
        server.server_close()
        stop.set()
        poller.join(timeout=5)
        store.close()


def main() -> None:
    setup_logging()
    store = BetsStore(config.DB_PATH, config.SOURCE_RUNS_RETENTION_DAYS)
    sources: list[Source] = [KalshiSource()]
    for optional in (betonline_source_if_configured(), novig_source_if_connected(), bfa_source_if_configured(),
                     wagerzon_source_if_configured(), polymarket_us_source_if_configured()):
        if optional is not None:
            sources.append(optional)
    serve(sources, store, config.BIND_HOST, config.PORT)


if __name__ == "__main__":
    main()
