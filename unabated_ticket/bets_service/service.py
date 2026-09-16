"""Entry point: the local bets service the Unabated Ticket panel polls (#114).

    python -m unabated_ticket.bets_service.service      (or bets_service/run.sh)

Inputs:  each registered Source (sources/kalshi.py; sources/betonline.py when
         its cookie file exists; sources/novig.py when its token file exists) on
         its own poll_sec.
Outputs: HTTP on 127.0.0.1:8094 (loopback only, no auth):
           GET /bets.json[?days=N]  {generatedAt, sources: {name: {fetchedAt, ok,
                                    error, count}}, bets: [records open + settled
                                    within N days (default RETENTION_DAYS=30)],
                                    crosswalk: [team_crosswalk rows, newest first]}
           GET /health              {ok, generatedAt, uptimeSec, sources}
           POST /crosswalk.json     body {rows: [{venue, league, venueTeamKey,
                                    unabatedTeamId, venueTeamName?, unabatedTeamName?,
                                    learnedFrom?}]} -> {ok, learned, conflicts, crosswalk}
                                    (#118 step 4: the panel's lessons from id joins;
                                    Content-Type must be application/json — a web page
                                    cannot send that cross-origin without a preflight
                                    this server never answers, so no site can write here)
           DELETE /crosswalk.json   -> {ok, cleared, crosswalk: []}
Side effects: UPSERTs records into bets.duckdb::bets and APPENDs a row to
bets.duckdb::source_runs per poll (see store.py); INSERTs into / DELETEs from
bets.duckdb::team_crosswalk on the two crosswalk routes; rotating log at
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
from unabated_ticket.bets_service.sources.kalshi import KalshiSource
from unabated_ticket.bets_service.sources.novig import source_if_connected as novig_source_if_connected
from unabated_ticket.bets_service.store import BetsStore

log = logging.getLogger(__name__)

POLL_TICK_SEC = 1.0
MAX_DAYS = 3650
# A crosswalk POST is a few rows per poll; anything near this is not the panel.
MAX_BODY_BYTES = 1024 * 1024
MAX_CROSSWALK_ROWS_PER_POST = 1000
CROSSWALK_REQUIRED_FIELDS = ("venue", "league", "venueTeamKey", "unabatedTeamId")
CROSSWALK_OPTIONAL_FIELDS = ("venueTeamName", "unabatedTeamName", "learnedFrom")


def _now() -> datetime:
    return datetime.now(timezone.utc)


def _iso(value: datetime) -> str:
    return value.strftime("%Y-%m-%dT%H:%M:%SZ")


def run_source_once(source: Source, store: BetsStore) -> bool:
    """One poll of one source: upsert its records and log the run. Never raises."""
    started_at = _now()
    try:
        records = source.fetch()
        store.upsert_bets(records, started_at)
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
            "bets": store.load_bets(days, now), "crosswalk": store.load_crosswalk()}


def _non_empty_string(value: object) -> bool:
    return isinstance(value, str) and value != ""


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
        team_id = raw.get("unabatedTeamId")
        if isinstance(team_id, int) and not isinstance(team_id, bool):
            team_id = str(team_id)
        row = {"unabatedTeamId": team_id}
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
            url = urlparse(self.path)
            if url.path != "/crosswalk.json":
                self._send_json(404, {"error": f"no route for POST {url.path}"})
                return
            body = self._read_json_body()
            if isinstance(body, tuple):
                self._send_json(*body)
                return
            rows = validate_crosswalk_rows(body)
            if isinstance(rows, str):
                self._send_json(400, {"error": rows})
                return
            result = store.learn_crosswalk(rows, _now())
            self._send_json(200, {"ok": True, "learned": result["learned"], "conflicts": result["conflicts"],
                                  "crosswalk": store.load_crosswalk()})

        def do_DELETE(self) -> None:  # noqa: N802 — http.server's name
            url = urlparse(self.path)
            if url.path != "/crosswalk.json":
                self._send_json(404, {"error": f"no route for DELETE {url.path}"})
                return
            cleared = store.clear_crosswalk()
            self._send_json(200, {"ok": True, "cleared": cleared, "crosswalk": []})

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
    store = BetsStore(config.DB_PATH)
    sources: list[Source] = [KalshiSource()]
    for optional in (betonline_source_if_configured(), novig_source_if_connected()):
        if optional is not None:
            sources.append(optional)
    serve(sources, store, config.BIND_HOST, config.PORT)


if __name__ == "__main__":
    main()
