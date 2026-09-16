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
from unabated_ticket.bets_service.store import UPSERT_ROWS_PER_STATEMENT, BetsStore
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


SOURCE_RUNS_RETENTION_DAYS = 7


@pytest.fixture
def store(tmp_path):
    bets_store = BetsStore(tmp_path / "bets.duckdb", SOURCE_RUNS_RETENTION_DAYS)
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


def test_upsert_writes_only_records_whose_content_changed(store):
    first = record("kalshi:a:yes", "open", None)
    seen = datetime(2026, 9, 15, 19, 0, tzinfo=timezone.utc)
    assert store.upsert_bets([first], seen) == 1
    # The same content again, a minute later: nothing to write. Kalshi re-sends
    # every record it has ever seen on every poll, and a DuckDB update is
    # copy-on-write — the WAL is what an identical rewrite costs (#125).
    assert store.upsert_bets([dict(first)], seen + timedelta(minutes=1)) == 0
    # `sourceFetchedAt` is the poll's own clock and changes every poll: were it
    # compared, no row would ever be skipped.
    restamped = {**first, "sourceFetchedAt": "2026-09-15T19:01:00Z"}
    assert store.upsert_bets([restamped], seen + timedelta(minutes=1)) == 0
    # A real change is written, and only it.
    settled = record("kalshi:a:yes", "won", "2026-09-15T20:00:00Z")
    assert store.upsert_bets([settled, record("kalshi:b:no", "open", None)], seen + timedelta(minutes=2)) == 2
    assert store.upsert_bets([settled, record("kalshi:b:no", "open", None)], seen + timedelta(minutes=3)) == 0
    assert store._con.execute("SELECT count(*) FROM bets").fetchone()[0] == 2
    [won] = [bet for bet in service.bets_payload(store, days=30)["bets"] if bet["id"] == "kalshi:a:yes"]
    assert won["status"] == "won"


def test_an_unchanged_record_serves_the_stamp_of_the_poll_that_last_changed_it(store):
    # The served sourceFetchedAt is the stored one — the last poll that CHANGED
    # the record — never the latest poll: the panel breaks merge ties on it,
    # and a record the source no longer returns must not read as fresh.
    seen = datetime(2026, 9, 15, 19, 0, tzinfo=timezone.utc)
    first = {**record("kalshi:a:yes", "open", None), "sourceFetchedAt": "2026-09-15T19:00:00Z"}
    store.upsert_bets([first], seen)
    store.upsert_bets([{**first, "sourceFetchedAt": "2026-09-15T19:01:00Z"}], seen + timedelta(minutes=1))
    [bet] = service.bets_payload(store, days=30, source_names=["kalshi"])["bets"]
    assert bet["sourceFetchedAt"] == "2026-09-15T19:00:00Z"
    settled = {**record("kalshi:a:yes", "won", "2026-09-15T20:00:00Z"), "sourceFetchedAt": "2026-09-15T20:02:00Z"}
    store.upsert_bets([settled], seen + timedelta(hours=1))
    [bet] = service.bets_payload(store, days=30, source_names=["kalshi"])["bets"]
    assert (bet["status"], bet["sourceFetchedAt"]) == ("won", "2026-09-15T20:02:00Z")


def test_a_poll_larger_than_one_upsert_statement_is_chunked(store):
    seen = datetime(2026, 9, 15, 19, 0, tzinfo=timezone.utc)
    records = [record(f"novig:{n}:yes", "open", None) for n in range(UPSERT_ROWS_PER_STATEMENT * 2 + 7)]
    assert store.upsert_bets(records, seen) == len(records)
    assert store.upsert_bets(records, seen + timedelta(minutes=1)) == 0
    records[-1] = record(records[-1]["id"], "won", "2026-09-15T20:00:00Z")
    assert store.upsert_bets(records, seen + timedelta(minutes=2)) == 1
    assert store._con.execute("SELECT count(*) FROM bets").fetchone()[0] == len(records)


def test_a_duplicated_id_in_one_poll_keeps_the_last_record_and_warns(store, caplog):
    # A multi-row INSERT ... ON CONFLICT silently keeps the FIRST row of a key
    # repeated in the statement; the old per-row executemany kept the last.
    seen = datetime(2026, 9, 15, 19, 0, tzinfo=timezone.utc)
    with caplog.at_level("WARNING"):
        assert store.upsert_bets([record("kalshi:a:yes", "open", None),
                                  record("kalshi:b:yes", "open", None),
                                  record("kalshi:a:yes", "won", "2026-09-15T20:00:00Z")], seen) == 2
    assert store._con.execute("SELECT status FROM bets WHERE id = 'kalshi:a:yes'").fetchone() == ("won",)
    assert "1 record id(s) repeated in one poll" in caplog.text and "kalshi:a:yes" in caplog.text
    assert store.upsert_bets([], seen + timedelta(minutes=1)) == 0  # an empty poll writes nothing
    # The count is the statement's own, not a scan keyed on seen_at: two polls
    # sharing a clock reading (two sources in one microsecond) do not add up.
    assert store.upsert_bets([record("kalshi:b:yes", "open", None)], seen) == 0


def test_rows_from_before_the_content_hash_column_are_rewritten_once(store):
    seen = datetime(2026, 9, 15, 19, 0, tzinfo=timezone.utc)
    store.upsert_bets([record("kalshi:a:yes", "open", None)], seen)
    store._con.execute("UPDATE bets SET content_hash = NULL")  # what a pre-#125 row looks like
    assert store.upsert_bets([record("kalshi:a:yes", "open", None)], seen + timedelta(minutes=1)) == 1
    assert store.upsert_bets([record("kalshi:a:yes", "open", None)], seen + timedelta(minutes=2)) == 0
    # An exact text compare: a value that Python would call equal is a change.
    assert store.upsert_bets([{**record("kalshi:a:yes", "open", None), "contracts": 3}], seen) == 1
    assert store.upsert_bets([{**record("kalshi:a:yes", "open", None), "contracts": 3.0}], seen) == 1


def log_run(store: BetsStore, source: str, at: datetime, ok: bool = True) -> None:
    store.log_source_run(source, at, at + timedelta(seconds=1), ok, None if ok else "RuntimeError: 503", 1 if ok else 0)


def count_runs(store: BetsStore) -> int:
    return store._con.execute("SELECT count(*) FROM source_runs").fetchone()[0]


def test_source_runs_are_pruned_to_the_retention_window_and_no_more_than_hourly(store):
    start = datetime(2026, 9, 1, 12, 0, tzinfo=timezone.utc)
    log_run(store, "kalshi", start)                                       # first run prunes; nothing is old yet
    log_run(store, "kalshi", start + timedelta(hours=2))
    log_run(store, "kalshi", start + timedelta(days=SOURCE_RUNS_RETENTION_DAYS, hours=3))  # prunes the first two
    assert count_runs(store) == 1
    # The freshness read is unchanged by the prune: the surviving run is the latest.
    status = store.source_status()["kalshi"]
    assert (status["ok"], status["count"], status["fetchedAt"]) == (True, 1, "2026-09-08T15:00:00Z")
    # At most one prune per SOURCE_RUNS_PRUNE_INTERVAL_SEC: a row old enough to
    # go survives a run logged inside the interval, and goes on the next one.
    ancient = start - timedelta(days=30)
    store._con.execute("INSERT INTO source_runs (source, started_at, finished_at, ok, error, n_records) "
                       "VALUES (?, ?, ?, ?, ?, ?)", ["kalshi", ancient, ancient, True, None, 1])
    log_run(store, "kalshi", start + timedelta(days=SOURCE_RUNS_RETENTION_DAYS, hours=3, minutes=30))
    assert count_runs(store) == 3
    log_run(store, "kalshi", start + timedelta(days=SOURCE_RUNS_RETENTION_DAYS, hours=5))
    assert count_runs(store) == 3  # the two recent runs plus this one; the ancient row is gone
    assert store._con.execute("SELECT count(*) FROM source_runs WHERE started_at < ?",
                              [start]).fetchone()[0] == 0


def test_prune_keeps_every_source_its_latest_run_and_latest_successful_run(store):
    # A source failing for longer than the window must keep showing its last
    # good read; a source no longer polled must keep its last run at all.
    start = datetime(2026, 9, 1, 12, 0, tzinfo=timezone.utc)
    log_run(store, "novig", start)                          # the last time Novig succeeded
    log_run(store, "betonline", start + timedelta(hours=1))  # BetOnline's last run ever
    for day in range(1, 2 * SOURCE_RUNS_RETENTION_DAYS + 1):
        log_run(store, "novig", start + timedelta(days=day), ok=False)
        log_run(store, "kalshi", start + timedelta(days=day))
    status = store.source_status()
    assert (status["novig"]["ok"], status["novig"]["fetchedAt"], status["novig"]["count"]) == (False, "2026-09-01T12:00:00Z", 1)
    assert status["novig"]["error"] == "RuntimeError: 503"
    assert (status["betonline"]["ok"], status["betonline"]["fetchedAt"]) == (True, "2026-09-01T13:00:00Z")
    assert status["kalshi"]["fetchedAt"] == "2026-09-15T12:00:00Z"
    # And everything else past the window is gone.
    assert store._con.execute("SELECT count(*) FROM source_runs WHERE started_at < ?",
                              [start + timedelta(days=SOURCE_RUNS_RETENTION_DAYS)]).fetchone()[0] == 2


def test_a_zero_retention_window_disables_the_prune_and_a_negative_one_is_refused(tmp_path):
    with pytest.raises(ValueError, match="0 \\(no prune\\) or more"):
        BetsStore(tmp_path / "bad.duckdb", -1)
    unpruned = BetsStore(tmp_path / "bets.duckdb", 0)
    try:
        start = datetime(2026, 9, 1, 12, 0, tzinfo=timezone.utc)
        log_run(unpruned, "kalshi", start)
        log_run(unpruned, "kalshi", start + timedelta(days=400))
        assert count_runs(unpruned) == 2
    finally:
        unpruned.close()


def test_a_failing_prune_is_logged_and_never_fails_the_poll(store, caplog):
    start = datetime(2026, 9, 1, 12, 0, tzinfo=timezone.utc)
    log_run(store, "kalshi", start)
    store._con.execute("ALTER TABLE source_runs RENAME TO source_runs_gone")  # the DELETE now raises
    try:
        with caplog.at_level("ERROR"):
            source = ScriptedSource([[record("kalshi:a:yes", "open", None)]])
            # The INSERT would fail too; exercise the prune path directly instead.
            assert store._prune_source_runs_locked(start + timedelta(hours=2)) == 0
    finally:
        store._con.execute("ALTER TABLE source_runs_gone RENAME TO source_runs")
    assert "prune of rows older than" in caplog.text
    assert service.run_source_once(source, store) is True


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


def request_json(method: str, url: str, body: bytes | None = None, content_type: str | None = None,
                 host: str | None = None) -> tuple[int, dict]:
    headers = {"Content-Type": content_type} if content_type else {}
    if host:
        headers["Host"] = host  # what a DNS-rebound page sends: its own name
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


def test_http_refuses_a_foreign_host_on_every_verb(store, http_server):
    """DNS rebinding: evil.example resolving to 127.0.0.1 is SAME-ORIGIN with
    this server, so CORS and the JSON Content-Type guard do not apply to it.
    The Host it sends is still its own name."""
    service.run_source_once(ScriptedSource([[record("kalshi:a:yes", "open", None)]]), store)
    store.learn_crosswalk([crosswalk_row("kalshi", "Wazzu", "717")], datetime.now(timezone.utc))
    body = json.dumps({"rows": [crosswalk_row("novig", "nv-1", "927")]}).encode()
    for method, path, payload, content_type in [
            ("GET", "/bets.json", None, None), ("GET", "/health", None, None),
            ("POST", "/crosswalk.json", body, "application/json"),
            ("DELETE", "/crosswalk.json", None, None)]:
        status, reply = request_json(method, f"{http_server}{path}", payload, content_type, host="evil.example")
        assert status == 403, f"{method} {path}"
        assert "Host must be one of" in reply["error"]
    # Nothing leaked and nothing was written or cleared.
    assert len(store.load_crosswalk()) == 1
    # The loopback names the service is actually serving on still pass.
    port = urlparse(http_server).port
    for host in [f"127.0.0.1:{port}", f"localhost:{port}", f"LOCALHOST:{port}"]:
        status, payload = request_json("GET", f"{http_server}/bets.json", host=host)
        assert (status, payload["bets"][0]["id"]) == (200, "kalshi:a:yes"), host


def test_host_allowed():
    assert service.host_allowed("127.0.0.1:8094", 8094) is True
    assert service.host_allowed("localhost:8094", 8094) is True
    assert service.host_allowed("evil.example", 8094) is False
    assert service.host_allowed("evil.example:8094", 8094) is False
    assert service.host_allowed("127.0.0.1:8095", 8094) is False  # another service's port
    assert service.host_allowed("127.0.0.1", 8094) is False       # port 80, not ours
    assert service.host_allowed(None, 8094) is False
    # On the scheme's default port a browser omits it from Host.
    assert service.host_allowed("127.0.0.1", 80) is True
    assert service.host_allowed("localhost:80", 80) is True
    assert service.host_allowed("evil.example", 80) is False


@pytest.mark.parametrize("raw_root, expected", [
    ("/u/NFLWork", "/u/NFLWork"),
    ("/u/NFLWork/.worktrees/feature-x", "/u/NFLWork"),
    # Claude desktop's layout: without this the service read no Kalshi credentials from a worktree.
    ("/u/NFLWork/.claude/worktrees/bet-venue-ids", "/u/NFLWork"),
])
def test_main_checkout_root_finds_the_main_checkout_from_either_worktree_layout(raw_root: str, expected: str) -> None:
    assert config.main_checkout_root(Path(raw_root)) == Path(expected)
