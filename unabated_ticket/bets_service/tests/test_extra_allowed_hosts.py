"""BETS_EXTRA_ALLOWED_HOSTS (phone page plan step 4): the tailnet name that
`tailscale serve` forwards as the browser's Host passes the #125 allowlist,
bare and with :443, while every other name is still a 403; a bad value stops
the service at import instead of running with a misread allowlist."""
import http.client
import os
import subprocess
import sys
import threading
from http.server import ThreadingHTTPServer
from pathlib import Path

import pytest

from unabated_ticket.bets_service import config, service
from unabated_ticket.bets_service.store import BetsStore

TAILNET_NAME = "mlb-stack.tail1234.ts.net"
REPO_ROOT = Path(__file__).resolve().parents[3]


@pytest.mark.parametrize("raw, expected", [
    (None, ()),
    ("", ()),
    ("  ,  ", ()),
    (TAILNET_NAME, (TAILNET_NAME,)),
    (f" {TAILNET_NAME} , other-vm.tail1234.ts.net ", (TAILNET_NAME, "other-vm.tail1234.ts.net")),
    (f"{TAILNET_NAME},{TAILNET_NAME}", (TAILNET_NAME,)),
    ("100.64.0.7", ("100.64.0.7",)),  # the tailnet IP is a valid name too
])
def test_parse_extra_allowed_hosts_accepts_bare_lowercase_names(raw, expected):
    assert config.parse_extra_allowed_hosts(raw) == expected


@pytest.mark.parametrize("raw", [
    "https://mlb-stack.tail1234.ts.net",   # scheme
    "mlb-stack.tail1234.ts.net/",          # path
    "mlb-stack.tail1234.ts.net:443",       # port (added by the service)
    "*.tail1234.ts.net",                   # wildcard
    "MLB-stack.tail1234.ts.net",           # not lowercase
    "-bad.ts.net",                         # label starts with a hyphen
    "bad..ts.net",                         # empty label
    "mlb stack.ts.net",                    # space
    f"{TAILNET_NAME},evil example",        # one bad entry fails the whole value
])
def test_parse_extra_allowed_hosts_refuses_anything_but_a_host_name(raw):
    with pytest.raises(ValueError, match="BETS_EXTRA_ALLOWED_HOSTS"):
        config.parse_extra_allowed_hosts(raw)


def test_a_bad_value_stops_the_service_at_import():
    env = {**os.environ, "BETS_EXTRA_ALLOWED_HOSTS": "https://evil.example"}
    result = subprocess.run([sys.executable, "-c", "import unabated_ticket.bets_service.config"],
                            cwd=REPO_ROOT, env=env, capture_output=True, text=True, check=False)
    assert result.returncode != 0
    assert "BETS_EXTRA_ALLOWED_HOSTS" in result.stderr
    assert "'https://evil.example'" in result.stderr


def test_host_allowed_with_the_tailnet_name():
    extra = (TAILNET_NAME,)
    assert service.host_allowed(TAILNET_NAME, 8094, extra) is True
    assert service.host_allowed(f"{TAILNET_NAME}:443", 8094, extra) is True
    assert service.host_allowed(TAILNET_NAME.upper(), 8094, extra) is True
    assert service.host_allowed("127.0.0.1:8094", 8094, extra) is True  # loopback still passes
    # The proxy speaks HTTPS on 443; the name on the service's own port is not what it sends.
    assert service.host_allowed(f"{TAILNET_NAME}:8094", 8094, extra) is False
    assert service.host_allowed(f"{TAILNET_NAME}:80", 8094, extra) is False
    # Look-alikes: a sub- or super-domain is another name.
    assert service.host_allowed(f"evil.{TAILNET_NAME}", 8094, extra) is False
    assert service.host_allowed(f"{TAILNET_NAME}.evil.example", 8094, extra) is False
    assert service.host_allowed("evil.example", 8094, extra) is False
    # Without the setting the tailnet name is foreign.
    assert service.host_allowed(TAILNET_NAME, 8094) is False


@pytest.fixture
def tailnet_server(tmp_path):
    store = BetsStore(tmp_path / "bets.duckdb", 7)
    server = ThreadingHTTPServer(("127.0.0.1", 0), service.make_handler(
        store, 0.0, ["kalshi"], extra_hosts=(TAILNET_NAME,)))
    thread = threading.Thread(target=server.serve_forever, daemon=True)
    thread.start()
    yield server.server_address[1]
    server.shutdown()
    server.server_close()
    store.close()


def get_status(port: int, path: str, host: str) -> tuple[int, bytes]:
    connection = http.client.HTTPConnection("127.0.0.1", port, timeout=5)
    try:
        connection.putrequest("GET", path, skip_host=True)
        connection.putheader("Host", host)
        connection.endheaders()
        response = connection.getresponse()
        return response.status, response.read()
    finally:
        connection.close()


@pytest.mark.parametrize("host", [TAILNET_NAME, f"{TAILNET_NAME}:443"])
def test_http_serves_the_tailnet_name_tailscale_serve_forwards(tailnet_server, host):
    for path in ("/health", "/bets.json", "/"):
        status, _body = get_status(tailnet_server, path, host)
        assert status == 200, f"{host} {path}"


@pytest.mark.parametrize("host", ["evil.example", f"evil.{TAILNET_NAME}", f"{TAILNET_NAME}:8094"])
def test_http_still_refuses_every_other_host(tailnet_server, host):
    status, body = get_status(tailnet_server, "/bets.json", host)
    assert status == 403
    assert b"Host must be one of" in body
