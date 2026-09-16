"""KalshiSource against a fake API, the DuckDB store's keep-previous-on-failure
rule, and the /bets.json + /health contract over a real HTTP server."""
import json
import threading
import time
import urllib.request
from datetime import datetime, timedelta, timezone
from http.server import ThreadingHTTPServer
from pathlib import Path
from urllib.parse import parse_qs, urlparse

import pytest

from unabated_ticket.bets_service import config, service
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


class FailingStore:
    """A store whose every write raises — what a full disk looks like to the loop."""

    def __init__(self):
        self.attempts = 0

    def upsert_bets(self, records, seen_at):
        self.attempts += 1
        raise RuntimeError("disk full")

    def log_source_run(self, *args):
        raise RuntimeError("disk full")


def test_poll_loop_survives_a_store_write_error_and_retries():
    source = ScriptedSource([[record("kalshi:a:yes", "open", None)]] * 3)
    failing = FailingStore()
    stop = threading.Event()
    thread = threading.Thread(target=service.poll_loop, args=([source], failing, stop), daemon=True)
    original_tick_sec = service.POLL_TICK_SEC
    service.POLL_TICK_SEC = 0.01
    source.poll_sec = 0.02
    try:
        thread.start()
        deadline = time.monotonic() + 2
        while failing.attempts < 2 and time.monotonic() < deadline:
            time.sleep(0.01)
    finally:
        stop.set()
        thread.join(timeout=2)
        service.POLL_TICK_SEC = original_tick_sec
    assert failing.attempts >= 2  # the thread outlived the first failure
    assert not thread.is_alive()


def test_registered_sources_report_no_completed_poll_yet(store):
    payload = service.bets_payload(store, days=30, source_names=["kalshi"])
    assert payload["sources"] == {"kalshi": service.NO_POLL_YET}
    service.run_source_once(ScriptedSource([[record("kalshi:a:yes", "open", None)]]), store)
    assert service.bets_payload(store, days=30, source_names=["kalshi"])["sources"]["kalshi"]["ok"] is True


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


# ---- team crosswalk (#118 step 4) -------------------------------------------------------

def crosswalk_row(venue: str, venue_team_key: str, unabated_team_id: str, **extra) -> dict:
    return {"venue": venue, "league": "cfb", "venueTeamKey": venue_team_key, "unabatedTeamId": unabated_team_id,
            "venueTeamName": venue_team_key, "unabatedTeamName": f"team {unabated_team_id}",
            "learnedFrom": "kalshi:x:yes on board event 123742", **extra}


def test_crosswalk_learns_once_per_venue_team_and_never_rewrites_a_held_id(store):
    now = datetime(2026, 9, 15, 19, 0, tzinfo=timezone.utc)
    first = store.learn_crosswalk([crosswalk_row("kalshi", "Wazzu (venue spelling)", "717"),
                                   crosswalk_row("novig", "nv-wsu", "717")], now)
    assert first == {"learned": 2, "conflicts": []}
    # The same rows again: no duplicates, nothing changes, no conflict.
    again = store.learn_crosswalk([crosswalk_row("kalshi", "Wazzu (venue spelling)", "717", learnedFrom="later bet")], now + timedelta(hours=1))
    assert again == {"learned": 0, "conflicts": []}
    # A different id for a held venue team is refused and reported, the held row untouched.
    refused = store.learn_crosswalk([crosswalk_row("kalshi", "Wazzu (venue spelling)", "1")], now + timedelta(hours=2))
    assert refused == {"learned": 0, "conflicts": [{"venue": "kalshi", "league": "cfb", "venueTeamKey": "Wazzu (venue spelling)",
                                                    "held": "717", "proposed": "1"}]}
    rows = store.load_crosswalk()
    assert [(row["venue"], row["venueTeamKey"], row["unabatedTeamId"], row["learnedAt"], row["learnedFrom"]) for row in rows] == [
        ("kalshi", "Wazzu (venue spelling)", "717", "2026-09-15T19:00:00Z", "kalshi:x:yes on board event 123742"),
        ("novig", "nv-wsu", "717", "2026-09-15T19:00:00Z", "kalshi:x:yes on board event 123742"),
    ]
    assert store._con.execute("SELECT count(*) FROM team_crosswalk").fetchone()[0] == 2
    assert set(rows[0]) == {"venue", "league", "venueTeamKey", "venueTeamName", "unabatedTeamId", "unabatedTeamName", "learnedFrom", "learnedAt"}
    # Served with every bets payload; cleared on request.
    assert service.bets_payload(store, days=30)["crosswalk"] == rows
    assert store.clear_crosswalk() == 2
    assert store.load_crosswalk() == []
    assert store.clear_crosswalk() == 0


def test_validate_crosswalk_rows_names_the_first_bad_row():
    good = crosswalk_row("kalshi", "Wazzu", "717")
    assert service.validate_crosswalk_rows({"rows": [good]}) == [good]
    # A numeric Unabated id is stored as its string; optional fields default to None.
    numeric = {"venue": "novig", "league": "cfb", "venueTeamKey": "nv-1", "unabatedTeamId": 717}
    assert service.validate_crosswalk_rows({"rows": [numeric]}) == [
        {**numeric, "unabatedTeamId": "717", "venueTeamName": None, "unabatedTeamName": None, "learnedFrom": None}]
    assert service.validate_crosswalk_rows([]) == "body must be an object with a `rows` array"
    assert service.validate_crosswalk_rows({"rows": [1]}) == "rows[0] must be an object"
    assert service.validate_crosswalk_rows({"rows": [good, {**good, "venue": ""}]}) == "rows[1].venue must be a non-empty string"
    assert service.validate_crosswalk_rows({"rows": [{**good, "unabatedTeamId": True}]}) == "rows[0].unabatedTeamId must be a non-empty string"
    assert service.validate_crosswalk_rows({"rows": [{**good, "learnedFrom": 3}]}) == "rows[0].learnedFrom must be a string or null"
    assert service.validate_crosswalk_rows({"rows": [good] * 1001}).startswith("at most 1000 rows")


# ---- HTTP ------------------------------------------------------------------------------

@pytest.fixture
def http_server(store):
    server = ThreadingHTTPServer(("127.0.0.1", 0), service.make_handler(store, 0.0, ["kalshi"]))
    thread = threading.Thread(target=server.serve_forever, daemon=True)
    thread.start()
    yield f"http://127.0.0.1:{server.server_address[1]}"
    server.shutdown()
    server.server_close()


def get_json(url: str) -> tuple[int, dict]:
    return request_json("GET", url)


def request_json(method: str, url: str, body: bytes | None = None, content_type: str | None = None) -> tuple[int, dict]:
    headers = {"Content-Type": content_type} if content_type else {}
    request = urllib.request.Request(url, data=body, method=method, headers=headers)
    try:
        with urllib.request.urlopen(request, timeout=5) as response:
            return response.status, json.loads(response.read())
    except urllib.error.HTTPError as error:
        return error.code, json.loads(error.read())


def test_http_bets_json_and_health_shape(store, http_server):
    service.run_source_once(ScriptedSource([[record("kalshi:a:yes", "open", None)]]), store)
    status, payload = get_json(f"{http_server}/bets.json")
    assert status == 200
    assert set(payload) == {"generatedAt", "sources", "bets", "crosswalk"}
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


def test_http_before_any_poll_serves_an_empty_list_with_the_source_pending(http_server):
    status, payload = get_json(f"{http_server}/bets.json")
    assert status == 200
    assert payload == {"generatedAt": payload["generatedAt"], "sources": {"kalshi": service.NO_POLL_YET}, "bets": [], "crosswalk": []}
    status, payload = get_json(f"{http_server}/health")
    assert status == 200 and payload["sources"] == {"kalshi": service.NO_POLL_YET}


def test_http_crosswalk_post_learns_and_delete_clears(store, http_server):
    rows = [crosswalk_row("kalshi", "Wazzu (venue spelling)", "717"), crosswalk_row("novig", "nv-duq", 927)]
    body = json.dumps({"rows": rows}).encode()
    status, reply = request_json("POST", f"{http_server}/crosswalk.json", body, "application/json; charset=utf-8")
    assert status == 200
    assert (reply["ok"], reply["learned"], reply["conflicts"]) == (True, 2, [])
    assert [(row["venue"], row["unabatedTeamId"]) for row in reply["crosswalk"]] == [("kalshi", "717"), ("novig", "927")]
    # The table rides along with every bets payload.
    status, payload = get_json(f"{http_server}/bets.json")
    assert status == 200 and len(payload["crosswalk"]) == 2
    # A conflicting re-learn is refused, reported, and the held row survives.
    conflict = json.dumps({"rows": [crosswalk_row("kalshi", "Wazzu (venue spelling)", "1")]}).encode()
    status, reply = request_json("POST", f"{http_server}/crosswalk.json", conflict, "application/json")
    assert status == 200 and reply["learned"] == 0 and reply["conflicts"][0]["held"] == "717"
    assert store.load_crosswalk()[0]["unabatedTeamId"] == "717" or store.load_crosswalk()[1]["unabatedTeamId"] == "717"
    status, reply = request_json("DELETE", f"{http_server}/crosswalk.json")
    assert status == 200 and reply == {"ok": True, "cleared": 2, "crosswalk": []}
    assert get_json(f"{http_server}/bets.json")[1]["crosswalk"] == []


def test_http_crosswalk_refuses_non_json_writes_bad_bodies_and_other_paths(http_server):
    body = json.dumps({"rows": [crosswalk_row("kalshi", "Wazzu", "717")]}).encode()
    # The CSRF guard: a web page's cross-origin POST can only be form- or text-encoded.
    status, reply = request_json("POST", f"{http_server}/crosswalk.json", body, "text/plain")
    assert status == 415 and "application/json" in reply["error"]
    status, reply = request_json("POST", f"{http_server}/crosswalk.json", body, "application/x-www-form-urlencoded")
    assert status == 415
    status, reply = request_json("POST", f"{http_server}/crosswalk.json", b"{not json", "application/json")
    assert status == 400 and reply["error"].startswith("body is not JSON")
    status, reply = request_json("POST", f"{http_server}/crosswalk.json", b'{"rows": [{"venue": "kalshi"}]}', "application/json")
    assert status == 400 and reply["error"] == "rows[0].league must be a non-empty string"
    status, reply = request_json("POST", f"{http_server}/bets.json", body, "application/json")
    assert status == 404
    status, reply = request_json("DELETE", f"{http_server}/bets.json")
    assert status == 404
    assert get_json(f"{http_server}/bets.json")[1]["crosswalk"] == []


@pytest.mark.parametrize("raw_root, expected", [
    ("/u/NFLWork", "/u/NFLWork"),
    ("/u/NFLWork/.worktrees/feature-x", "/u/NFLWork"),
    # Claude desktop's layout: without this the service read no Kalshi credentials from a worktree.
    ("/u/NFLWork/.claude/worktrees/bet-venue-ids", "/u/NFLWork"),
])
def test_main_checkout_root_finds_the_main_checkout_from_either_worktree_layout(raw_root: str, expected: str) -> None:
    assert config.main_checkout_root(Path(raw_root)) == Path(expected)
