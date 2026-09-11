"""KalshiSource against a fake API, the DuckDB store's keep-previous-on-failure
rule, and the /bets.json + /health contract over a real HTTP server."""
import json
import threading
import urllib.request
from datetime import datetime, timedelta, timezone
from http.server import ThreadingHTTPServer
from urllib.parse import parse_qs, urlparse

import pytest

from unabated_ticket.bets_service import service
from unabated_ticket.bets_service.sources.kalshi import KalshiSource
from unabated_ticket.bets_service.store import BetsStore
from unabated_ticket.bets_service.tests.conftest import fill, position

NE_FILL = fill("KXNFLGAME-26SEP20PITNE-NE", "yes", 482.42, 0.64, "2026-09-11T15:05:04.400105Z", trade_id="t1")
NE_FILL_2 = fill("KXNFLGAME-26SEP20PITNE-NE", "yes", 15.28, 0.64, "2026-09-11T15:05:04.400105Z", trade_id="t2")
NE_POSITION = position("KXNFLGAME-26SEP20PITNE-NE", 550, 352)


class FakeKalshiApi:
    """Answers the signed GETs the source makes; records every path it saw."""

    def __init__(self, fixture: dict, fills: list[dict], positions: list[dict]):
        self.fixture = fixture
        self.fills = fills
        self.positions = positions
        self.paths: list[str] = []
        self.fail_paths: set[str] = set()

    def __call__(self, method: str, path: str) -> tuple[int, object, dict]:
        assert method == "GET"
        self.paths.append(path)
        url = urlparse(path)
        if url.path in self.fail_paths:
            return 503, "unavailable", {}
        if url.path == "/portfolio/fills":
            min_ts = parse_qs(url.query).get("min_ts")
            fills = self.fills
            if min_ts:
                fills = [f for f in fills if datetime.fromisoformat(
                    f["created_time"].replace("Z", "+00:00")).timestamp() >= int(min_ts[0])]
            return 200, {"fills": fills, "cursor": ""}, {}
        if url.path == "/portfolio/positions":
            return 200, {"market_positions": self.positions, "cursor": ""}, {}
        if url.path.startswith("/markets/"):
            market = self.fixture["markets"].get(url.path[len("/markets/"):])
            return (200, {"market": market}, {}) if market else (404, {"message": "not found"}, {})
        if url.path.startswith("/events/"):
            event = self.fixture["events"].get(url.path[len("/events/"):])
            return (200, {"event": event}, {}) if event else (404, {"message": "not found"}, {})
        raise AssertionError(f"unexpected path {path}")

    def count(self, prefix: str) -> int:
        return sum(1 for path in self.paths if path.startswith(prefix))


def make_source(api, clock, poll_sec=60, reconcile_sec=3600):
    return KalshiSource(api=api, poll_sec=poll_sec, reconcile_sec=reconcile_sec,
                        configure_auth=False, lookup_gap_sec=0, clock=clock)


# ---- KalshiSource -------------------------------------------------------------------

def test_first_fetch_is_a_full_pull_then_incremental_with_overlap(kalshi_fixture):
    api = FakeKalshiApi(kalshi_fixture, [NE_FILL, NE_FILL_2], [NE_POSITION])
    now = [1_800_000_000.0]
    source = make_source(api, clock=lambda: now[0])
    [record] = source.fetch()
    assert record["id"] == "kalshi:KXNFLGAME-26SEP20PITNE-NE:yes"
    assert (record["stake"], record["price"], record["contracts"]) == (352, -178, 550)
    assert api.paths[0] == "/portfolio/fills?limit=200"  # full pull, no min_ts
    now[0] += 60
    source.fetch()
    incremental = [p for p in api.paths if p.startswith("/portfolio/fills?min_ts=")]
    assert incremental == [f"/portfolio/fills?min_ts={int(1_800_000_000 - 60)}&limit=200"]
    # Market and event were looked up exactly once across both polls.
    assert api.count("/markets/") == 1
    assert api.count("/events/") == 1


def test_fills_are_deduped_on_trade_id_and_a_missing_id_fails_loudly(kalshi_fixture):
    api = FakeKalshiApi(kalshi_fixture, [NE_FILL, NE_FILL], [NE_POSITION])
    now = [1_800_000_000.0]
    source = make_source(api, clock=lambda: now[0])
    [record] = source.fetch()
    assert record["raw"]["fillCount"] == 1
    api.fills = [fill("KXNFLGAME-26SEP20PITNE-NE", "yes", 1, 0.64, "2026-09-11T15:05:04Z")]
    now[0] += 3600  # next full pull sees the id-less fill
    with pytest.raises(RuntimeError, match="trade_id"):
        source.fetch()


def test_reconcile_re_pulls_fills_and_re_reads_unsettled_markets(kalshi_fixture):
    api = FakeKalshiApi(kalshi_fixture, [NE_FILL], [NE_POSITION])
    now = [1_800_000_000.0]
    source = make_source(api, clock=lambda: now[0], reconcile_sec=3600)
    source.fetch()
    assert source.fetch()[0]["status"] == "open"
    # Kalshi settles the market: the position vanishes and the result lands.
    api.fixture["markets"]["KXNFLGAME-26SEP20PITNE-NE"] = {
        **kalshi_fixture["markets"]["KXNFLGAME-26SEP20PITNE-NE"], "result": "yes", "status": "finalized"}
    api.positions = []
    now[0] += 3600
    [record] = source.fetch()
    assert api.count("/portfolio/fills?limit=200") == 2  # second full pull
    assert api.count("/markets/") == 2  # re-read because it had no result
    assert record["status"] == "won"
    assert record["closedAt"] == "2026-09-20T20:00:00Z"


def test_a_failed_market_lookup_fails_closed_for_that_record_only(kalshi_fixture):
    other = fill("KXNCAAFTOTAL-26SEP12RICEND-60", "yes", 500, 0.32, "2026-09-11T20:08:40Z", trade_id="t9")
    api = FakeKalshiApi(kalshi_fixture, [NE_FILL, other], [NE_POSITION])
    api.fail_paths.add("/markets/KXNFLGAME-26SEP20PITNE-NE")
    source = make_source(api, clock=lambda: 1_800_000_000.0)
    records = {record["id"]: record for record in source.fetch()}
    assert records["kalshi:KXNFLGAME-26SEP20PITNE-NE:yes"]["unmatchable"] == (
        "unreadable Kalshi market (no market payload for KXNFLGAME-26SEP20PITNE-NE)")
    assert records["kalshi:KXNCAAFTOTAL-26SEP12RICEND-60:yes"]["unmatchable"] is None


def test_a_failed_fills_pull_raises(kalshi_fixture):
    api = FakeKalshiApi(kalshi_fixture, [NE_FILL], [NE_POSITION])
    api.fail_paths.add("/portfolio/fills")
    with pytest.raises(RuntimeError, match="expected 200"):
        make_source(api, clock=lambda: 1_800_000_000.0).fetch()


# ---- store + service loop ------------------------------------------------------------

class ScriptedSource:
    name = "kalshi"
    poll_sec = 60

    def __init__(self, outcomes: list):
        self.outcomes = list(outcomes)

    def fetch(self) -> list[dict]:
        outcome = self.outcomes.pop(0)
        if isinstance(outcome, Exception):
            raise outcome
        return outcome


def record(record_id: str, status: str, closed_at: str | None, placed_at: str = "2026-09-11T15:05:04Z") -> dict:
    return {"id": record_id, "source": "kalshi_api", "venue": "kalshi", "league": "nfl", "status": status,
            "placedAt": placed_at, "closedAt": closed_at, "eventDate": "2026-09-20",
            "sourceFetchedAt": "2026-09-11T21:00:00Z"}


@pytest.fixture
def store(tmp_path):
    bets_store = BetsStore(tmp_path / "bets.duckdb")
    yield bets_store
    bets_store.close()


def test_failed_poll_keeps_previous_records_and_logs_a_failed_run(store):
    source = ScriptedSource([[record("kalshi:a:yes", "open", None)], RuntimeError("Kalshi 503")])
    assert service.run_source_once(source, store) is True
    assert service.run_source_once(source, store) is False
    payload = service.bets_payload(store, days=30)
    assert [bet["id"] for bet in payload["bets"]] == ["kalshi:a:yes"]
    status = payload["sources"]["kalshi"]
    assert status["ok"] is False
    assert status["error"] == "RuntimeError: Kalshi 503"
    assert status["count"] == 1
    assert status["fetchedAt"] is not None  # the last SUCCESSFUL run
    runs = store._con.execute("SELECT ok, error, n_records FROM source_runs ORDER BY started_at").fetchall()
    assert runs == [(True, None, 1), (False, "RuntimeError: Kalshi 503", 0)]


def test_upsert_replaces_a_record_by_id(store):
    source = ScriptedSource([[record("kalshi:a:yes", "open", None)],
                             [record("kalshi:a:yes", "won", "2026-09-12T00:00:00Z")]])
    service.run_source_once(source, store)
    service.run_source_once(source, store)
    [bet] = service.bets_payload(store, days=30)["bets"]
    assert bet["status"] == "won"
    assert store._con.execute("SELECT count(*) FROM bets").fetchone()[0] == 1


def test_bets_json_window_keeps_open_and_recently_closed(store):
    now = datetime.now(timezone.utc)
    recent = (now - timedelta(days=3)).strftime("%Y-%m-%dT%H:%M:%SZ")
    old = (now - timedelta(days=40)).strftime("%Y-%m-%dT%H:%M:%SZ")
    source = ScriptedSource([[
        record("kalshi:open:yes", "open", None),
        record("kalshi:recent:yes", "lost", recent),
        record("kalshi:old:yes", "won", old),
        record("kalshi:undated:yes", "closed", None),
    ]])
    service.run_source_once(source, store)
    ids = lambda days: sorted(bet["id"] for bet in service.bets_payload(store, days)["bets"])  # noqa: E731
    assert ids(30) == ["kalshi:open:yes", "kalshi:recent:yes", "kalshi:undated:yes"]
    assert ids(1) == ["kalshi:open:yes", "kalshi:undated:yes"]
    assert ids(60) == ["kalshi:old:yes", "kalshi:open:yes", "kalshi:recent:yes", "kalshi:undated:yes"]


def test_parse_days():
    assert service.parse_days("") == 30
    assert service.parse_days("days=7") == 7
    assert service.parse_days("days=x").startswith("days must be an integer")
    assert service.parse_days("days=-1").startswith("days must be between")


# ---- HTTP ------------------------------------------------------------------------------

@pytest.fixture
def http_server(store):
    server = ThreadingHTTPServer(("127.0.0.1", 0), service.make_handler(store, 0.0))
    thread = threading.Thread(target=server.serve_forever, daemon=True)
    thread.start()
    yield f"http://127.0.0.1:{server.server_address[1]}"
    server.shutdown()
    server.server_close()


def get_json(url: str) -> tuple[int, dict]:
    try:
        with urllib.request.urlopen(url, timeout=5) as response:
            return response.status, json.loads(response.read())
    except urllib.error.HTTPError as error:
        return error.code, json.loads(error.read())


def test_http_bets_json_and_health_shape(store, http_server):
    service.run_source_once(ScriptedSource([[record("kalshi:a:yes", "open", None)]]), store)
    status, payload = get_json(f"{http_server}/bets.json")
    assert status == 200
    assert set(payload) == {"generatedAt", "sources", "bets"}
    assert set(payload["sources"]["kalshi"]) == {"fetchedAt", "ok", "error", "count"}
    assert payload["sources"]["kalshi"]["ok"] is True
    assert payload["bets"][0]["id"] == "kalshi:a:yes"
    status, payload = get_json(f"{http_server}/bets.json?days=0")
    assert status == 200 and payload["bets"][0]["id"] == "kalshi:a:yes"  # open bets ignore the window
    status, payload = get_json(f"{http_server}/bets.json?days=abc")
    assert status == 400 and "days" in payload["error"]
    status, payload = get_json(f"{http_server}/health")
    assert status == 200 and payload["ok"] is True and "kalshi" in payload["sources"]
    status, _payload = get_json(f"{http_server}/nope")
    assert status == 404


def test_http_before_any_poll_serves_an_empty_list(http_server):
    status, payload = get_json(f"{http_server}/bets.json")
    assert status == 200
    assert payload == {"generatedAt": payload["generatedAt"], "sources": {}, "bets": []}
