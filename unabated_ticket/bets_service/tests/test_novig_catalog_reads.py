"""Keyed Novig catalog reads for the runner (novig_catalog_reads.py): only the
four catalog GETs pass validation, each is signed with the key and sent as a
GET, Novig's answer (a 429 too) comes back as is, and POST /novig_read.json
serves them loopback-only (503 without the bot key, 502 when Novig is silent)."""
import json
import threading
import urllib.error
import urllib.request
from http.server import ThreadingHTTPServer

import pytest
from unabated_ticket.bets_service import service
from unabated_ticket.bets_service.novig_catalog_reads import (
    NovigCatalogReads,
    pooled_transport,
    validate_read_request,
)
from unabated_ticket.bets_service.novig_orders import NovigOrderClient
from unabated_ticket.bets_service.store import BetsStore
from unabated_ticket.bets_service.tests.test_novig_orders import KEYPAIRS

MARKET = "01a1758e-f87f-7602-a07f-33183f9f8960"
EVENT = "01a16b43-def9-7992-9e3b-95b2b04f0756"


@pytest.mark.parametrize("body, expected", [
    ({"path": "/v3/catalog/markets", "query": {"league": "NBA,NCAAB", "eventStatus": "OPEN_PREGAME", "limit": "5000", "after": "c1"}},
     ("/v3/catalog/markets", [("league", "NBA,NCAAB"), ("eventStatus", "OPEN_PREGAME"), ("limit", "5000"), ("after", "c1")])),
    ({"path": f"/v3/catalog/markets/{MARKET}"}, (f"/v3/catalog/markets/{MARKET}", [])),
    ({"path": f"/v3/catalog/markets/{MARKET}/book", "query": None}, (f"/v3/catalog/markets/{MARKET}/book", [])),
    ({"path": f"/v3/catalog/events/{EVENT}"}, (f"/v3/catalog/events/{EVENT}", [])),
])
def test_the_four_catalog_reads_pass(body, expected):
    assert validate_read_request(body) == expected


@pytest.mark.parametrize("body, refusal", [
    ("x", "JSON object"),
    ({"path": "/v3/orders"}, "path must be"),
    ({"path": "/v3/portfolio/positions"}, "path must be"),
    ({"path": f"/v3/catalog/markets/{MARKET}/../../orders"}, "path must be"),
    ({"path": f"/v3/public/catalog/markets/{MARKET}"}, "path must be"),
    ({"path": "/v3/catalog/markets/a?b"}, "path must be"),
    ({"path": f"/v3/catalog/events/{EVENT}\n"}, "path must be"),
    ({"path": "/v3/catalog/markets/a%2Fb"}, "path must be"),
    ({"path": "/v3/catalog/markets/.."}, "path must be"),
    ({"path": f"/v3/catalog/markets/{MARKET}", "query": {"league": "NBA"}}, "only /v3/catalog/markets takes a query"),
    ({"path": "/v3/catalog/markets", "query": {"account": "x"}}, "is not one of"),
    ({"path": "/v3/catalog/markets", "query": {"league": 5}}, "non-empty string"),
    ({"path": "/v3/catalog/markets", "query": ["league"]}, "query must be an object"),
])
def test_anything_else_is_refused_before_signing(body, refusal):
    assert refusal in validate_read_request(body)


def _reads(tmp_path, replies):
    key_path = tmp_path / "bot.pem"
    key_path.write_text(KEYPAIRS["ed25519-test-1"]["private_key_pkcs8_pem"])
    sent = []

    def transport(method, url, headers, body):
        sent.append({"method": method, "url": url, "headers": headers, "body": body})
        status, reply = replies.pop(0)
        return status, json.dumps(reply).encode()

    client = NovigOrderClient("4365ed85-test", key_path, "https://api.novig.com", transport=transport)
    return NovigCatalogReads(client), sent


def test_a_read_is_a_signed_get_and_novigs_answer_comes_back_as_is(tmp_path):
    reads, sent = _reads(tmp_path, [(200, {"items": [], "next": None}), (429, {"code": "RATE_LIMIT_EXCEEDED"})])
    assert reads.read("/v3/catalog/markets", [("league", "NBA"), ("eventStatus", "OPEN_PREGAME")]) == {
        "status": 200, "body": {"items": [], "next": None}}
    assert reads.read(f"/v3/catalog/events/{EVENT}", []) == {"status": 429, "body": {"code": "RATE_LIMIT_EXCEEDED"}}
    assert [request["method"] for request in sent] == ["GET", "GET"]
    assert sent[0]["url"] == "https://api.novig.com/v3/catalog/markets?league=NBA&eventStatus=OPEN_PREGAME"
    assert sent[0]["headers"]["Novig-Key-Id"] == "4365ed85-test" and sent[0]["headers"]["Novig-Signature"]
    assert sent[0]["body"] is None


def test_the_pooled_transport_runs_reads_at_once_on_separate_sessions():
    made, inside, both_in = [], [], threading.Barrier(2, timeout=5)

    def make_transport():
        session = len(made)
        made.append(session)

        def send(method, url, headers, body):
            inside.append(session)
            both_in.wait()
            return 200, b"{}"

        return send

    send = pooled_transport(make_transport, size=2)
    threads = [threading.Thread(target=send, args=("GET", "u", {}, None)) for _ in range(2)]
    for thread in threads:
        thread.start()
    for thread in threads:
        thread.join(5)
    assert sorted(inside) == [0, 1], "two reads at once went to two sessions"


def test_a_read_with_every_session_busy_is_refused_after_the_wait():
    release = threading.Event()

    def make_transport():
        def send(method, url, headers, body):
            release.wait(5)
            return 200, b"{}"
        return send

    send = pooled_transport(make_transport, size=1, wait_sec=0.05)
    busy = threading.Thread(target=send, args=("GET", "u", {}, None))
    busy.start()
    try:
        with pytest.raises(OSError, match="all 1 Novig read sessions busy for 0.05 s"):
            send("GET", "u", {}, None)
    finally:
        release.set()
        busy.join(5)


def _serve(handler):
    server = ThreadingHTTPServer(("127.0.0.1", 0), handler)
    threading.Thread(target=server.serve_forever, daemon=True).start()
    return server, f"http://127.0.0.1:{server.server_address[1]}"


def _post(url, payload, host=None):
    headers = {"Content-Type": "application/json"}
    if host:
        headers["Host"] = host
    request = urllib.request.Request(url, data=json.dumps(payload).encode(), method="POST", headers=headers)
    try:
        with urllib.request.urlopen(request, timeout=5) as response:
            return response.status, json.loads(response.read())
    except urllib.error.HTTPError as error:
        return error.code, json.loads(error.read())


@pytest.fixture
def bets_store(tmp_path):
    store = BetsStore(tmp_path / "bets.duckdb", 7)
    yield store
    store.close()


class _Silent:
    def read(self, path, query):
        raise OSError("ConnectTimeout: no answer")


def test_post_novig_read_serves_keyed_reads_loopback_only(bets_store, tmp_path):
    reads, sent = _reads(tmp_path, [(200, {"marketId": MARKET, "status": "OPEN"}), (429, {"code": "RATE_LIMIT_EXCEEDED"})])
    server, base = _serve(service.make_handler(bets_store, 0.0, ["kalshi"], novig_catalog_reads=reads, extra_hosts=("mlb-stack.test",)))
    try:
        market = _post(f"{base}/novig_read.json", {"path": f"/v3/catalog/markets/{MARKET}"})
        limited = _post(f"{base}/novig_read.json", {"path": f"/v3/catalog/events/{EVENT}"})
        bad = _post(f"{base}/novig_read.json", {"path": "/v3/orders"})
        foreign = _post(f"{base}/novig_read.json", {"path": f"/v3/catalog/markets/{MARKET}"}, host="mlb-stack.test")
    finally:
        server.shutdown()
        server.server_close()
    assert market == (200, {"status": 200, "body": {"marketId": MARKET, "status": "OPEN"}})
    assert limited == (200, {"status": 429, "body": {"code": "RATE_LIMIT_EXCEEDED"}}), "Novig's 429 passes through for the runner to name"
    assert bad[0] == 400 and "path must be" in bad[1]["error"]
    assert foreign[0] == 403, "the tailnet name never reaches the bot key"
    assert len(sent) == 2, "the refused bodies were never signed or sent"


def test_post_novig_read_is_503_without_the_bot_key_and_502_when_novig_is_silent(bets_store):
    for reads, expected_status, expected_text in ((None, 503, "NOVIG_TAKER_KEY_ID"), (_Silent(), 502, "Novig did not answer")):
        server, base = _serve(service.make_handler(bets_store, 0.0, ["kalshi"], novig_catalog_reads=reads))
        try:
            status, body = _post(f"{base}/novig_read.json", {"path": f"/v3/catalog/markets/{MARKET}"})
        finally:
            server.shutdown()
            server.server_close()
        assert status == expected_status and expected_text in body["error"]
