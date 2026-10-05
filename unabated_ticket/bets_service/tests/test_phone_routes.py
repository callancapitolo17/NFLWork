"""The phone page's routes on the bets service (phone page plan step 2): the
static files from a fixed allowlist (anything else, traversal included, is a
404), the Host allowlist in front of them, and GET /edges.json proxied to a
stub runner — passed through on a 200, a 502 naming the runner on a refused
connection, a timeout or a non-200 answer."""
import http.client
import json
import socket
import threading
import time
from http.server import BaseHTTPRequestHandler, ThreadingHTTPServer

import pytest

from unabated_ticket.bets_service import service
from unabated_ticket.bets_service.store import BetsStore

RUNNER_BODY = {"generatedAt": "2026-10-03T12:00:00.000Z", "items": [], "total": 0}
RUNNER_TIMEOUT_SEC = 0.3
SLOW_RUNNER_DELAY_SEC = 1.5


class StubRunnerHandler(BaseHTTPRequestHandler):
    """GET /edges.json as the runner would: 200 + RUNNER_BODY, or the server's `mode`."""

    def do_GET(self) -> None:  # noqa: N802 — http.server's name
        mode = self.server.mode
        self.server.hosts_seen.append(self.headers.get("Host"))
        if mode == "slow":
            time.sleep(SLOW_RUNNER_DELAY_SEC)
        status = 500 if mode == "error" else 200
        body = json.dumps({"error": "boom"} if mode == "error" else RUNNER_BODY).encode()
        try:
            self.send_response(status)
            self.send_header("Content-Type", "application/json")
            self.send_header("Content-Length", str(len(body)))
            self.end_headers()
            self.wfile.write(body)
        except OSError:  # the proxy gave up on a slow reply
            pass

    def log_message(self, format: str, *args: object) -> None:  # noqa: A002
        pass


def start_server(server: ThreadingHTTPServer) -> threading.Thread:
    server.daemon_threads = True
    thread = threading.Thread(target=server.serve_forever, daemon=True)
    thread.start()
    return thread


@pytest.fixture
def stub_runner():
    server = ThreadingHTTPServer(("127.0.0.1", 0), StubRunnerHandler)
    server.mode = "ok"
    server.hosts_seen = []
    start_server(server)
    yield server
    server.shutdown()
    server.server_close()


@pytest.fixture
def store(tmp_path):
    bets_store = BetsStore(tmp_path / "bets.duckdb", 7)
    yield bets_store
    bets_store.close()


def serve_bets(store: BetsStore, runner_url: str):
    server = ThreadingHTTPServer(("127.0.0.1", 0), service.make_handler(
        store, 0.0, ["kalshi"], runner_url=runner_url, runner_timeout_sec=RUNNER_TIMEOUT_SEC))
    start_server(server)
    return server


@pytest.fixture
def bets_port(store, stub_runner):
    server = serve_bets(store, f"http://127.0.0.1:{stub_runner.server_address[1]}")
    yield server.server_address[1]
    server.shutdown()
    server.server_close()


def get(port: int, path: str, host: str | None = None) -> tuple[int, dict, bytes]:
    """A raw GET, so the path reaches the server exactly as written (no client normalising `..`)."""
    connection = http.client.HTTPConnection("127.0.0.1", port, timeout=5)
    connection.putrequest("GET", path, skip_host=True)
    connection.putheader("Host", host or f"127.0.0.1:{port}")
    connection.endheaders()
    response = connection.getresponse()
    body = response.read()
    headers = {name.lower(): value for name, value in response.getheaders()}
    connection.close()
    return response.status, headers, body


def unused_port() -> int:
    with socket.socket() as probe:
        probe.bind(("127.0.0.1", 0))
        return probe.getsockname()[1]


# ---- static files ----------------------------------------------------------------------

@pytest.mark.parametrize("path, content_type", [
    ("/", "text/html; charset=utf-8"),
    ("/index.html", "text/html; charset=utf-8"),
    ("/phone.css", "text/css; charset=utf-8"),
    ("/phone.js", "text/javascript; charset=utf-8"),
    ("/phoneview.js", "text/javascript; charset=utf-8"),
    ("/ext/kelly.js", "text/javascript; charset=utf-8"),
    ("/ext/edgerows.js", "text/javascript; charset=utf-8"),
    ("/tracker", "text/html; charset=utf-8"),
    ("/tracker/", "text/html; charset=utf-8"),
    ("/tracker/tracker.css", "text/css; charset=utf-8"),
    ("/tracker/tracker.js", "text/javascript; charset=utf-8"),
    ("/tracker/trackerstats.js", "text/javascript; charset=utf-8"),
])
def test_every_allowlisted_file_is_served_whole_with_its_type_and_no_cache(bets_port, path, content_type):
    status, headers, body = get(bets_port, path)
    assert status == 200
    assert headers["content-type"] == content_type
    assert headers["cache-control"] == "no-store"
    assert headers["x-content-type-options"] == "nosniff"
    assert "default-src 'self'" in headers["content-security-policy"]
    assert body == service.STATIC_FILES[path].read_bytes()


def test_the_page_loads_every_extension_module_it_is_served_and_nothing_else():
    """index.html's <script src="ext/..."> list and the allowlist must agree:
    a module missing from the allowlist would 404 on the phone, one missing
    from the page would be served for nothing."""
    page = (service.PHONE_DIR / "index.html").read_text()
    loaded = [chunk.split('"', 1)[0] for chunk in page.split('<script src="ext/')[1:]]
    assert loaded == list(service.PHONE_EXTENSION_MODULES)


def test_the_tracker_page_loads_only_allowlisted_files_by_absolute_path():
    """The tracker is served at both /tracker and /tracker/, so its assets must
    be absolute paths (a relative one would resolve to / from the first) and
    each must be on the allowlist."""
    page = (service.TRACKER_DIR / "index.html").read_text()
    assets = [chunk.split('"', 1)[0] for marker in ('<script src="', 'rel="stylesheet" href="')
              for chunk in page.split(marker)[1:]]
    assert sorted(assets) == ["/tracker/tracker.css", "/tracker/tracker.js", "/tracker/trackerstats.js"]
    assert all(asset in service.STATIC_FILES for asset in assets)


@pytest.mark.parametrize("path", [
    "/tracker/../bets_service/service.py",
    "/tracker/index.html",
    "/../bets_service/service.py",
    "/ext/../bets_service/config.py",
    "/ext/..%2fbets_service%2fconfig.py",
    "/ext/%2e%2e/bets_service/.env",
    "/ext/panel.js",
    "/ext/manifest.json",
    "/phone/../../bets.duckdb",
    "//etc/passwd",
    "/phone.css/",
    "/PHONE.CSS",
    "/bets.duckdb",
])
def test_anything_off_the_allowlist_is_a_404_traversal_included(bets_port, path):
    status, headers, body = get(bets_port, path)
    assert status == 404
    assert headers["content-type"] == "application/json"
    assert json.loads(body)["error"].startswith("no route for ")


@pytest.mark.parametrize("path", ["/", "/phone.js", "/ext/kelly.js", "/edges.json"])
def test_a_foreign_host_is_refused_before_any_file_or_the_runner(bets_port, stub_runner, path):
    status, _headers, body = get(bets_port, path, host="evil.example:8094")
    assert status == 403
    assert "Host must be one of" in json.loads(body)["error"]
    assert stub_runner.hosts_seen == []


# ---- /edges.json -------------------------------------------------------------------------

def test_edges_json_passes_the_runner_body_through(bets_port, stub_runner):
    status, headers, body = get(bets_port, "/edges.json")
    assert status == 200
    assert headers["content-type"] == "application/json"
    assert headers["cache-control"] == "no-store"
    assert json.loads(body) == RUNNER_BODY
    # The runner's own Host allowlist accepts what the proxy sends.
    assert stub_runner.hosts_seen == [f"127.0.0.1:{stub_runner.server_address[1]}"]


def test_edges_json_is_a_502_naming_the_runner_when_it_is_down(store):
    runner_url = f"http://127.0.0.1:{unused_port()}"
    server = serve_bets(store, runner_url)
    try:
        status, headers, body = get(server.server_address[1], "/edges.json")
    finally:
        server.shutdown()
        server.server_close()
    assert status == 502
    assert headers["content-type"] == "application/json"
    reply = json.loads(body)
    assert reply["runnerUrl"] == runner_url
    assert reply["error"].startswith(f"server runner at {runner_url} unreachable (")
    assert "node unabated_ticket/server/runner.js" in reply["error"]


def test_edges_json_is_a_502_when_the_runner_is_slower_than_the_timeout(bets_port, stub_runner):
    stub_runner.mode = "slow"
    started = time.monotonic()
    status, _headers, body = get(bets_port, "/edges.json")
    assert time.monotonic() - started < SLOW_RUNNER_DELAY_SEC
    assert status == 502
    assert "unreachable" in json.loads(body)["error"]


def test_edges_json_is_a_502_naming_the_runner_status_on_a_runner_error(bets_port, stub_runner):
    stub_runner.mode = "error"
    status, _headers, body = get(bets_port, "/edges.json")
    assert status == 502
    runner_url = f"http://127.0.0.1:{stub_runner.server_address[1]}"
    assert json.loads(body)["error"].startswith(f"server runner at {runner_url} answered HTTP 500 for /edges.json")

