"""Dry-run parlay quoter (#283): scope rule, pricer client, decision math,
decision log, the feed's door, the leg table and the decision loop."""
import io
import json
import urllib.error
from datetime import datetime, timedelta, timezone

import pytest

import queue
import time

from kalshi_parlay_quoter import decide, feed, legtable, pricer, scope
from kalshi_parlay_quoter.db import DecisionLog
from kalshi_parlay_quoter.main import Quoter, stored_in_sample


def leg(ticker, side="yes"):
    return {"market_ticker": ticker, "event_ticker": ticker.rsplit("-", 1)[0],
            "side": side}


NFL_SPREAD = "KXNFLSPREAD-26NOV09KCDEN-KC7"
NFL_TOTAL_SAME = "KXNFLTOTAL-26NOV09KCDEN-45"
NBA_ML = "KXNBAGAME-26NOV08LALBOS-LAL"
CBB_TOTAL = "KXNCAAMBTOTAL-26NOV08DUKEUNC-150"
CFB_SPREAD = "KXNCAAFSPREAD-26NOV08OSUMICH-OSU3"
MLB_TOTAL = "KXMLBTOTAL-26APR011905NYYBOS-9"


def rfq(rid="r1", legs=None, cost="500.00"):
    m = {"id": rid, "market_ticker": "KXMVECROSSCATEGORY-X-" + rid,
         "creator_id": "c1", "contracts_fp": "0.00",
         "target_cost_dollars": cost, "created_ts": "2026-11-07T18:00:00Z"}
    if legs is not None:
        m["mve_selected_legs"] = legs
    return m


def priced(ticker, side, fair, lo=None, hi=None, sigma=0.01, ok=True,
           reason=None, kickoff_ms=1_800_000_000_000):
    return {"market_ticker": ticker, "side": side, "ok": ok, "reason": reason,
            "fair": fair, "fair_lo": lo if lo is not None else fair,
            "fair_hi": hi if hi is not None else fair,
            "bracketed": lo is not None, "sigma": sigma, "c": 0.2, "tau": 1.1,
            "line": {"league": "NFL", "market": "spread", "number": -7.0},
            "kickoff_ms": kickoff_ms, "fair_age_ms": 1200,
            "kalshi_status": "active", "kalshi_yes_bid": 0.48,
            "kalshi_yes_ask": 0.52, "kalshi_no_ask": 0.52}


def fixed_pricer(fairs):
    calls = []

    def fn(rid, legs):
        calls.append((rid, legs))
        return pricer.PriceResult(
            [priced(t, s, fairs[t]) for t, s in legs], 123, 4.0)
    fn.calls = calls
    return fn


# -- scope -------------------------------------------------------------------

@pytest.mark.parametrize("ticker,league,market", [
    (NFL_SPREAD, "NFL", "spread"), (NBA_ML, "NBA", "ml"),
    (CBB_TOTAL, "CBB", "total"), (CFB_SPREAD, "CFB", "spread"),
    (MLB_TOTAL, "MLB", "total"), ("KXNFLGAME-26NOV09KCDEN-KC", "NFL", "ml"),
])
def test_classify_full_game_series(ticker, league, market):
    info = scope.classify_leg(leg(ticker))
    assert (info.league, info.market) == (league, market)


@pytest.mark.parametrize("ticker", [
    "KXNCAAF1HTOTAL-26NOV08OSUMICH-24",   # first half
    "KXMLBF5TOTAL-26APR011905NYYBOS-4",   # first five
    "KXNBATEAMTOTAL-26NOV08LALBOS-LAL110",  # team total prop
    "KXNHLGAME-26NOV08BOSTOR-BOS",        # league Unabated fairs are not wired for
    "KXNFLSPREAD",                        # malformed
])
def test_out_of_scope_legs(ticker):
    assert scope.classify_leg(leg(ticker)).league is None


def test_cross_game_cross_league_in_scope_any_leg_count():
    legs = [leg(NFL_SPREAD), leg(NBA_ML, "no"), leg(CBB_TOTAL), leg(CFB_SPREAD),
            leg(MLB_TOTAL)]
    v = scope.evaluate(legs)
    assert v.in_scope and v.n_games == 5 and not v.same_game
    assert v.leagues == ["NFL", "NBA", "CBB", "CFB", "MLB"]


def test_same_game_declined():
    v = scope.evaluate([leg(NFL_SPREAD), leg(NFL_TOTAL_SAME)])
    assert (v.in_scope, v.reason, v.same_game) == (False, "same_game", True)
    assert v.detail == "NFL:26NOV09KCDEN"


def test_other_declines():
    assert scope.evaluate(None).reason == "no_legs"
    assert scope.evaluate([leg(NFL_SPREAD)]).reason == "single_leg"
    v = scope.evaluate([leg(NFL_SPREAD), leg("KXNCAAF1HTOTAL-26NOV08OSUMICH-24")])
    assert (v.reason, v.detail) == ("unsupported_leg", "KXNCAAF1HTOTAL")


# -- pricer client -----------------------------------------------------------

class _Resp(io.BytesIO):
    def __enter__(self):
        return self

    def __exit__(self, *a):
        return False


def test_pricer_builds_get_url_and_parses():
    seen = {}
    body = {"rfq_id": "r1", "priced_at_ms": 9,
            "legs": [priced(NFL_SPREAD, "yes", 0.5), priced(NBA_ML, "no", 0.4)]}

    def opener(req, timeout):
        seen.update(url=req.full_url, method=req.get_method(), timeout=timeout)
        return _Resp(json.dumps(body).encode())

    res = pricer.price("http://127.0.0.1:8095", "r1",
                       [(NFL_SPREAD, "yes"), (NBA_ML, "no")],
                       timeout=1.5, opener=opener)
    assert seen["method"] == "GET" and seen["timeout"] == 1.5
    assert seen["url"] == ("http://127.0.0.1:8095/rfq_price.json?rfq_id=r1"
                           f"&legs={NFL_SPREAD}:yes,{NBA_ML}:no")
    assert res.priced_at_ms == 9 and len(res.legs) == 2


@pytest.mark.parametrize("make,err", [
    (lambda r, timeout: (_ for _ in ()).throw(urllib.error.HTTPError(
        r.full_url, 400, "bad", {}, io.BytesIO(b""))), "http_400"),
    (lambda r, timeout: (_ for _ in ()).throw(urllib.error.URLError("refused")), "transport"),
    (lambda r, timeout: _Resp(b"not json"), "bad_json"),
    (lambda r, timeout: _Resp(json.dumps({"legs": []}).encode()), "leg_count_mismatch"),
    (lambda r, timeout: _Resp(json.dumps({"legs": [priced(NBA_ML, "yes", .5)]}).encode()),
     "leg_order_mismatch"),
])
def test_pricer_errors(make, err):
    with pytest.raises(pricer.PricerError, match=err):
        pricer.price("http://x", "r1", [(NFL_SPREAD, "yes")], timeout=1,
                     opener=make)


# -- decision ----------------------------------------------------------------

def test_would_quote_parlay_fair_is_product():
    fn = fixed_pricer({NFL_SPREAD: 0.5, NBA_ML: 0.4})
    row, legs = decide.decide(rfq(legs=[leg(NFL_SPREAD), leg(NBA_ML, "no")]),
                              [leg(NFL_SPREAD), leg(NBA_ML, "no")], fn,
                              now_ms=1_800_000_000_000 - 3_600_000)
    assert row["decision"] == "would_quote" and row["reason"] is None
    assert row["parlay_fair"] == pytest.approx(0.2)
    assert row["hours_to_kickoff"] == pytest.approx(1.0)
    assert row["legs_kalshi_ask_product"] == pytest.approx(0.52 * 0.52)
    assert row["target_cost_usd"] == 500.0 and row["leagues"] == "NFL,NBA"
    assert fn.calls == [("r1", [(NFL_SPREAD, "yes"), (NBA_ML, "no")])]
    assert [l["fair"] for l in legs] == [0.5, 0.4]


def test_bracketed_bounds_multiply():
    def fn(rid, ls):
        return pricer.PriceResult([priced(NFL_SPREAD, "yes", 0.5, 0.45, 0.55),
                                   priced(NBA_ML, "yes", 0.4)], 1, 1.0)
    row, _ = decide.decide(rfq(), [leg(NFL_SPREAD), leg(NBA_ML)], fn)
    assert row["parlay_fair_lo"] == pytest.approx(0.18)
    assert row["parlay_fair_hi"] == pytest.approx(0.22)


def test_out_of_scope_never_calls_pricer():
    fn = fixed_pricer({})
    row, legs = decide.decide(rfq(), [leg(NFL_SPREAD), leg(NFL_TOTAL_SAME)], fn)
    assert (row["decision"], row["reason"]) == ("decline", "same_game")
    assert fn.calls == [] and len(legs) == 2


def test_unpriced_leg_declines_with_its_reason():
    def fn(rid, ls):
        return pricer.PriceResult(
            [priced(NFL_SPREAD, "yes", 0.5),
             priced(NBA_ML, "yes", None, ok=False, reason="past_ladder")], 1, 1.0)
    row, legs = decide.decide(rfq(), [leg(NFL_SPREAD), leg(NBA_ML)], fn)
    assert (row["decision"], row["reason"]) == ("decline", "leg_past_ladder")
    assert row["parlay_fair"] is None and legs[0]["fair"] == 0.5


def test_pricer_error_declines():
    def fn(rid, ls):
        raise pricer.PricerError("transport: refused")
    row, _ = decide.decide(rfq(), [leg(NFL_SPREAD), leg(NBA_ML)], fn)
    assert (row["reason"], row["reason_detail"]) == ("pricer_error", "transport: refused")


# -- decision log -------------------------------------------------------------


def test_existing_db_gains_detail_column(tmp_path):
    import duckdb
    from kalshi_parlay_quoter import db as dbmod
    path = tmp_path / "old.duckdb"
    con = duckdb.connect(str(path))
    con.execute(dbmod.SCHEMA_SQL.split("-- Added after")[0])  # the #288 schema
    assert "detail" not in [r[0] for r in con.execute("DESCRIBE rfq_legs").fetchall()]
    con.close()
    dlog = DecisionLog(path)
    row, legs = decide.decide(rfq(), [leg(NFL_SPREAD), leg(NBA_ML)],
                              fixed_pricer({NFL_SPREAD: .5, NBA_ML: .4}))
    legs[0]["detail"] = "on rung"
    dlog.add(row, legs)
    dlog.flush()
    assert dlog.query("SELECT detail FROM rfq_legs WHERE leg_idx=0")[0][0] == "on rung"
    dlog.close()


def test_url_timeout_maps_to_pricer_timeout():
    def opener(req, timeout):
        raise urllib.error.URLError(TimeoutError("timed out"))
    with pytest.raises(pricer.PricerTimeout):
        pricer.price("http://x", "r1", [(NFL_SPREAD, "yes")], timeout=1, opener=opener)


def test_flush_writes_more_than_one_chunk(tmp_path):
    from kalshi_parlay_quoter import db as dbmod
    dlog = DecisionLog(tmp_path / "q.duckdb")
    fn = fixed_pricer({NFL_SPREAD: .5, NBA_ML: .4})
    n = dbmod.CHUNK_ROWS + 7
    for i in range(n):
        dlog.add(*decide.decide(rfq(f"r{i}"), [leg(NFL_SPREAD), leg(NBA_ML)], fn))
    dlog.add(*decide.decide(rfq("r0"), [leg(NFL_SPREAD), leg(NBA_ML)], fn))  # dup ignored
    assert dlog.flush() == n + 1
    assert dlog.query("SELECT count(*) FROM rfq_decisions")[0][0] == n
    assert dlog.query("SELECT count(*) FROM rfq_legs")[0][0] == 2 * n
    dlog.close()


class _BrokenBody:
    def __init__(self, exc):
        self.exc = exc

    def __enter__(self):
        return self

    def __exit__(self, *a):
        return False

    def read(self):
        raise self.exc


def test_mid_read_failures_map_to_pricer_errors():
    import http.client
    with pytest.raises(pricer.PricerTimeout):
        pricer.price("http://x", "r1", [(NFL_SPREAD, "yes")], timeout=1,
                     opener=lambda r, timeout: _BrokenBody(TimeoutError("read timed out")))
    with pytest.raises(pricer.PricerError, match="IncompleteRead"):
        pricer.price("http://x", "r1", [(NFL_SPREAD, "yes")], timeout=1,
                     opener=lambda r, timeout: _BrokenBody(http.client.IncompleteRead(b"")))


def test_timeout_logs_pricer_error_with_timeout_detail():
    def fn(rid, ls):
        raise pricer.PricerTimeout("timeout")
    row, _ = decide.decide(rfq(), [leg(NFL_SPREAD), leg(NBA_ML)], fn)
    assert (row["decision"], row["reason"], row["reason_detail"]) == (
        "decline", "pricer_error", "timeout")


def test_reused_price_does_not_outlive_kickoff():
    kickoff = 1_800_000_000_000
    fn = fixed_pricer({NFL_SPREAD: .5, NBA_ML: .4})
    row, legs = decide.decide(rfq(), [leg(NFL_SPREAD), leg(NBA_ML)], fn,
                              now_ms=kickoff + 1)
    assert (row["decision"], row["reason"]) == ("decline", "leg_started")
    assert legs[0]["detail"] == "started after it was priced"


def test_in_chunk_duplicate_keeps_first_row_and_columns_line_up(tmp_path):
    dlog = DecisionLog(tmp_path / "q.duckdb")
    fn = fixed_pricer({NFL_SPREAD: .5, NBA_ML: .4})
    first = decide.decide(rfq("dup", cost="500.00"), [leg(NFL_SPREAD), leg(NBA_ML)], fn)
    second = decide.decide(rfq("dup", cost="999.00"), [leg(NFL_SPREAD), leg(NBA_ML)], fn)
    dlog.add(*first)
    dlog.add(*second)                                 # same chunk, same rfq_id
    dlog.flush()
    assert dlog.query("SELECT target_cost_usd, decision, parlay_fair, n_legs, leagues "
                      "FROM rfq_decisions") == [(500.0, "would_quote", 0.2, 2, "NFL,NBA")]
    assert dlog.query("SELECT leg_idx, market_ticker, side, fair, kalshi_no_ask, ok "
                      "FROM rfq_legs ORDER BY leg_idx") == [
        (0, NFL_SPREAD, "yes", 0.5, 0.52, True), (1, NBA_ML, "yes", 0.4, 0.52, True)]
    dlog.close()


# -- fairs client ---------------------------------------------------------------

def test_fairs_builds_its_own_url_and_refuses_an_oversized_chunk():
    seen = {}

    def opener(req, timeout):
        seen.update(url=req.full_url, timeout=timeout)
        return _Resp(json.dumps({"priced_at_ms": 7, "legs": [priced(NFL_SPREAD, "no", 0.5)]}).encode())

    res = pricer.fairs("http://127.0.0.1:8095", [(NFL_SPREAD, "no")], timeout=15, opener=opener)
    assert seen == {"url": f"http://127.0.0.1:8095/rfq_fairs.json?legs={NFL_SPREAD}:no", "timeout": 15}
    assert res.priced_at_ms == 7
    with pytest.raises(ValueError, match="at most 200"):
        pricer.fairs("http://x", [(NFL_SPREAD, "yes")] * 201, timeout=1, opener=opener)


def test_an_answer_for_the_other_side_is_refused():
    opener = lambda r, timeout: _Resp(json.dumps({"legs": [priced(NFL_SPREAD, "yes", .5)]}).encode())
    with pytest.raises(pricer.PricerError, match="leg_order_mismatch"):
        pricer.fairs("http://x", [(NFL_SPREAD, "no")], timeout=1, opener=opener)


# -- feed door ------------------------------------------------------------------

def frame(msg, kind="rfq_created"):
    return json.dumps({"type": kind, "sid": 1, "msg": msg}, separators=(",", ":")).encode()


@pytest.mark.parametrize("raw,key", [
    (frame(rfq("a", [leg(NFL_SPREAD), leg(NBA_ML)])), "kept"),
    (frame(rfq("a", [leg(NFL_SPREAD), leg(NBA_ML)], cost="49.99")), "door_below_floor"),
    (frame(rfq("a", [leg("KXNHLGAME-26NOV08BOSTOR-BOS"), leg("KXEPLGAME-26NOV08ARSCHE-ARS")])), "door_other_league"),
    (frame(rfq("a", [leg(NFL_SPREAD), leg("KXNHLGAME-26NOV08BOSTOR-BOS")])), "door_unsupported_leg"),
    (frame(rfq("a", [leg(NFL_SPREAD), leg("KXNFL1HSPREAD-26NOV09KCDEN-KC4")])), "door_unsupported_leg"),
    (frame(rfq("a", [leg(NFL_SPREAD), leg(NFL_TOTAL_SAME)])), "door_same_game"),
    (frame(rfq("a", [leg(NFL_SPREAD)])), "door_single_leg"),
    (frame(rfq("a", [leg(NFL_SPREAD), leg(NBA_ML, side="maybe")])), "door_bad_side"),
    (frame(rfq("a")), "door_other_league"),                      # legless
    # Shapes the runner's TICKER_RE refuses must not reach the leg table.
    (frame(rfq("a", [leg(NFL_SPREAD), leg("KXNBAGAME-26NOV08LALBOS-LAL-G2")])), "door_unsupported_leg"),
    (frame(rfq("a", [leg(NFL_SPREAD), leg("KXNBAGAME-26NOV08LALBOS")])), "door_unsupported_leg"),
    (frame(rfq("a", [leg(NFL_SPREAD), leg("KXNBAGAME-26nov08LALBOS-LAL")])), "door_unsupported_leg"),
    (frame(rfq("a", [leg(NFL_SPREAD), leg("KXNBAGAME-26NOV08LALBOS-LAL ")])), "door_unsupported_leg"),
    (frame({"id": "a", "market_ticker": "x"}, "rfq_deleted"), "frames_deleted"),
    (frame({"channel": "communications", "sid": 1}, "subscribed"), "frames_other"),
    (b"not json", "malformed"),
    (b'{"type":"rfq_created","msg":"KXNFLSPREAD-X-Y' + b'"' + b'-"}', "malformed"),
    (frame({"market_ticker": "x", "mve_selected_legs": [leg(NFL_SPREAD), leg(NBA_ML)]}), "malformed"),  # no id
])
def test_door_keys(raw, key):
    assert feed.door(raw, 50)[0] == key


def test_door_contract_sized_rfqs_meet_the_floor():
    small = rfq("small", [leg(NFL_SPREAD), leg(NBA_ML)], cost="0")
    small["contracts_fp"] = "20.00"                 # pays at most $20 < $50
    big = dict(small, id="big", contracts_fp="205.00")
    assert feed.door(frame(small), 50)[0] == "door_below_floor"
    key, kept = feed.door(frame(big), 50)
    assert (key, kept["id"]) == ("kept", "big")
    assert feed.door(frame(small), 0)[0] == "kept", "a zero floor keeps every size"
    sizeless = dict(small, contracts_fp="0.00")
    unreadable = dict(small, contracts_fp="lots")
    assert [feed.door(frame(one), 50)[0] for one in (sizeless, unreadable)] == ["door_below_floor"] * 2


class FakeConnection:
    """Serves `frames`, then raises `then` forever."""
    def __init__(self, frames, then=TimeoutError):
        self.frames, self.then, self.sent, self.closed = list(frames), then, [], False

    def send(self, text):
        self.sent.append(json.loads(text))

    def recv(self, timeout):
        if self.frames:
            return self.frames.pop(0)
        raise self.then()

    def close(self):
        self.closed = True


def _wait_for(predicate, what, seconds=3.0):
    import time
    deadline = time.monotonic() + seconds
    while time.monotonic() < deadline:
        if predicate():
            return
        time.sleep(0.01)
    raise AssertionError(f"timed out waiting for {what}")


def test_feed_subscribes_queues_kept_rfqs_and_counts_the_rest():
    ack = frame({"channel": "communications", "sid": 1}, "subscribed")
    conn = FakeConnection([ack,
                           frame(rfq("r1", [leg(NFL_SPREAD), leg(NBA_ML)])),
                           frame(rfq("r2", [leg(NFL_SPREAD), leg(NBA_ML)], cost="10")),
                           frame({"id": "r1"}, "rfq_deleted")])
    f = feed.RFQFeed(lambda: conn, min_cost_usd=50, queue_max=10, idle_timeout_sec=30)
    f.start()
    _wait_for(lambda: f.rfqs.qsize() == 1 and not conn.frames, "the frames to be read")
    assert conn.sent == [{"id": 1, "cmd": "subscribe", "params": {"channels": ["communications"]}}]
    kept, read_at = f.rfqs.get_nowait()
    assert kept["id"] == "r1" and read_at.tzinfo is timezone.utc
    assert f.take_counts() == {"frames_other": 1, "kept": 1, "door_below_floor": 1, "frames_deleted": 1}
    assert f.take_counts() == {}, "counts are per call"
    f.stop()
    assert conn.closed


def test_feed_full_queue_drops_and_counts_instead_of_blocking():
    f = feed.RFQFeed(lambda: None, min_cost_usd=50, queue_max=1, idle_timeout_sec=30)
    f.handle_frame(frame(rfq("r1", [leg(NFL_SPREAD), leg(NBA_ML)])))
    f.handle_frame(frame(rfq("r2", [leg(NFL_SPREAD), leg(NBA_ML)])))
    assert f.rfqs.qsize() == 1
    assert f.take_counts() == {"kept": 1, "queue_dropped": 1}


def test_feed_reconnects_after_an_error_frame_and_stays_closed_while_paused(monkeypatch):
    monkeypatch.setattr(feed, "RECONNECT_BASE_SEC", 0.01)
    made = []

    def connect():
        first = not made
        conn = FakeConnection([frame({"code": 6, "msg": "buffer overflow"}, "error")] if first else
                              [frame({"channel": "communications"}, "subscribed")])
        made.append(conn)
        return conn

    f = feed.RFQFeed(connect, min_cost_usd=50, queue_max=10, idle_timeout_sec=30)
    f.start()
    _wait_for(lambda: len(made) == 2 and f.connected, "the reconnect")
    assert made[0].closed
    assert f.take_counts()["reconnects"] == 1
    f.set_paused(True)
    _wait_for(lambda: made[1].closed, "the paused feed to close its socket", seconds=8)
    assert len(made) == 2, "paused: no new connection"
    f.stop()


def test_feed_silent_socket_is_presumed_dead(monkeypatch):
    monkeypatch.setattr(feed, "RECONNECT_BASE_SEC", 0.01)
    clock = [0.0]
    made = []

    def connect():
        made.append(FakeConnection([frame({"channel": "communications"}, "subscribed")]))
        return made[-1]

    f = feed.RFQFeed(connect, min_cost_usd=50, queue_max=10, idle_timeout_sec=30, clock=lambda: clock[0])
    f.start()
    _wait_for(lambda: f.connected, "the first connection")
    clock[0] = 31.0                                   # every recv now times out past the idle limit
    _wait_for(lambda: len(made) >= 2, "the reconnect")
    f.stop()


# -- leg table ------------------------------------------------------------------

def test_leg_table_serves_true_ages_and_refuses_old_or_missing_fairs():
    clock = [100.0]
    table = legtable.LegTable(max_age_sec=15, leg_ttl_sec=600, clock=lambda: clock[0])
    keys = [(NFL_SPREAD, "yes"), (NBA_ML, "no")]
    table.want(keys)
    assert table.missing(keys) == keys
    with pytest.raises(pricer.LegsNotReady, match="2 of 2 legs"):
        table.price("r1", keys)
    table.store(keys[:1], [priced(NFL_SPREAD, "yes", .5)], 1000, sent_at=100.0)
    clock[0] = 104.0
    table.store(keys[1:], [priced(NBA_ML, "no", .4)], 5000, sent_at=104.0)
    clock[0] = 106.0
    res = table.price("r1", keys)
    assert [l["cache_age_ms"] for l in res.legs] == [6000.0, 2000.0]
    assert [l["fair_age_ms"] for l in res.legs] == [7200.0, 3200.0], "runner's 1200 ms plus the reuse"
    assert [l["priced_at_ms"] for l in res.legs] == [1000, 5000]
    assert (res.priced_at_ms, res.elapsed_ms) == (1000, 0.0), "the RFQ takes its oldest leg's stamp"
    clock[0] = 115.5                                  # the first fair is now 15.5 s old
    assert table.missing(keys) == keys[:1]
    with pytest.raises(pricer.LegsNotReady, match=NFL_SPREAD):
        table.price("r1", keys)
    assert isinstance(pricer.LegsNotReady(""), pricer.PricerError)


def test_leg_table_forgets_legs_nobody_asked_for_lately():
    clock = [0.0]
    table = legtable.LegTable(max_age_sec=15, leg_ttl_sec=600, clock=lambda: clock[0])
    table.want([(NFL_SPREAD, "yes")])
    table.store([(NFL_SPREAD, "yes")], [priced(NFL_SPREAD, "yes", .5)], 1, sent_at=0.0)
    clock[0] = 500.0
    table.want([(NBA_ML, "no")])
    clock[0] = 601.0
    assert table.wanted() == [(NBA_ML, "no")] and len(table) == 1
    # A refresh that was in flight when the leg expired must not bring it back.
    table.store([(NFL_SPREAD, "yes")], [priced(NFL_SPREAD, "yes", .5)], 1, sent_at=601.0)
    assert table.missing([(NFL_SPREAD, "yes")]) == [(NFL_SPREAD, "yes")]


def test_refresher_prices_every_wanted_leg_in_chunks_and_a_failed_chunk_ends_the_pass():
    table = legtable.LegTable(max_age_sec=15, leg_ttl_sec=600)
    keys = [(f"KXNBATOTAL-26NOV08LALBOS-{n}", "yes") for n in range(1, 451)]
    table.want(keys)
    calls = []

    def fetch(chunk):
        calls.append(len(chunk))
        if len(calls) == 2 and fail[0]:
            raise pricer.PricerTimeout("timeout")
        return pricer.PriceResult([priced(t, s, .5) for t, s in chunk], 9, 1.0)

    fail = [True]
    refresher = legtable.Refresher(table, fetch, interval_sec=2)
    assert refresher.refresh_once() is False
    assert calls == [200, 200] and len(table.missing(keys)) == 250
    assert refresher.take_stats()["table_refresh_failed"] == 1
    fail[0] = False
    calls.clear()
    generation = table.generation
    assert refresher.refresh_once() is True
    assert calls == [200, 200, 50] and table.missing(keys) == []
    assert table.generation > generation
    assert refresher.take_stats()["table_refresh_failed"] == 0


# -- decision loop --------------------------------------------------------------

def put(rfqs, one, age_sec=0.0):
    """Queue an RFQ the way the feed does: with the time its frame was read."""
    rfqs.put((one, datetime.now(timezone.utc) - timedelta(seconds=age_sec)))


def table_pricer(fairs, **leg_overrides):
    """A table fetch: board-only legs (no Kalshi fields)."""
    def fetch(chunk):
        legs = []
        for ticker, side in chunk:
            one = priced(ticker, side, fairs[ticker], **leg_overrides.get(ticker, {}))
            one.update(kalshi_status=None, kalshi_yes_bid=None, kalshi_yes_ask=None, kalshi_no_ask=None)
            legs.append(one)
        return pricer.PriceResult(legs, 123, 1.0)
    return fetch


def make_quoter(tmp_path, fairs, *, confirm_fn=None, clock=None, **kwargs):
    dlog = DecisionLog(tmp_path / "q.duckdb")
    table = legtable.LegTable(max_age_sec=15, leg_ttl_sec=600, **({"clock": clock} if clock else {}))
    refresher = legtable.Refresher(table, table_pricer(fairs), interval_sec=2,
                                   **({"clock": clock} if clock else {}))
    rfqs = queue.Queue()
    kwargs.setdefault("stored_sample_pct", 100)
    kwargs.setdefault("confirms_per_sec", 0)
    quoter = Quoter(rfqs, table, dlog, confirm_fn=confirm_fn or fixed_pricer(fairs),
                    **kwargs, **({"clock": clock} if clock else {}))
    return quoter, rfqs, refresher, dlog


FAIRS = {NFL_SPREAD: 0.5, NBA_ML: 0.4, CBB_TOTAL: 0.5, CFB_SPREAD: 0.5}


def test_new_legs_wait_one_refresh_then_price_from_memory(tmp_path):
    quoter, rfqs, refresher, dlog = make_quoter(tmp_path, FAIRS)
    put(rfqs, rfq("r1", [leg(NFL_SPREAD), leg(NBA_ML)]))
    assert quoter.step() == 1
    assert quoter.take_counts()["gauge_pending"] == 1, "parked: no fair for its legs yet"
    refresher.refresh_once()
    quoter.step()
    put(rfqs, rfq("r2", [leg(NBA_ML), leg(NFL_SPREAD, side="no")]))   # NO side is its own leg
    quoter.step()
    assert quoter.take_counts()["gauge_pending"] == 1
    refresher.refresh_once()
    put(rfqs, rfq("r3", [leg(NFL_SPREAD), leg(NBA_ML)]))              # both legs known: no wait
    quoter.step()
    dlog.flush()
    rows = dlog.query("SELECT rfq_id, decision, parlay_fair, confirmed, legs_kalshi_ask_product, queue_ms >= 0 "
                      "FROM rfq_decisions ORDER BY rfq_id")
    assert rows == [("r1", "would_quote", 0.2, False, None, True), ("r2", "would_quote", 0.2, False, None, True),
                    ("r3", "would_quote", 0.2, False, None, True)]
    ages = dlog.query("SELECT count(*), min(cache_age_ms) >= 0, count(kalshi_yes_ask) FROM rfq_legs")
    assert ages == [(6, True, 0)], "table legs carry their age and no Kalshi price"
    dlog.close()


def test_a_leg_never_priced_declines_pricer_error_after_the_wait(tmp_path):
    clock = [0.0]
    quoter, rfqs, refresher, dlog = make_quoter(tmp_path, FAIRS, clock=lambda: clock[0])
    put(rfqs, rfq("r1", [leg(NFL_SPREAD), leg(NBA_ML)]))
    quoter.step()
    clock[0] = 4.9
    quoter.step()
    assert quoter.take_counts().get("pricer_error") is None, "still inside the 5 s wait"
    clock[0] = 5.0
    quoter.step()
    counts = quoter.take_counts()
    assert (counts["pricer_error"], counts["gauge_pending"]) == (1, 0)
    dlog.flush()
    reason, detail = dlog.query("SELECT reason, reason_detail FROM rfq_decisions")[0]
    assert reason == "pricer_error" and "2 of 2 legs" in detail
    dlog.close()


def test_a_full_waiting_room_logs_pricer_backlog(tmp_path):
    quoter, rfqs, refresher, dlog = make_quoter(tmp_path, FAIRS, max_backlog=1)
    put(rfqs, rfq("r1", [leg(NFL_SPREAD), leg(NBA_ML)]))
    put(rfqs, rfq("r2", [leg(CBB_TOTAL), leg(CFB_SPREAD)]))
    quoter.step()
    counts = quoter.take_counts()
    assert (counts["pricer_backlog"], counts["gauge_pending"]) == (1, 1)
    dlog.close()


def test_a_declined_leg_declines_the_parlay_and_is_not_confirmed(tmp_path):
    dlog = DecisionLog(tmp_path / "q.duckdb")
    table = legtable.LegTable(max_age_sec=15, leg_ttl_sec=600)
    fetch = table_pricer(FAIRS, **{NBA_ML: {"ok": False, "reason": "no_event"}})
    confirm = fixed_pricer(FAIRS)
    rfqs = queue.Queue()
    quoter = Quoter(rfqs, table, dlog, confirm_fn=confirm, confirms_per_sec=1000, stored_sample_pct=100)
    table.want([(NFL_SPREAD, "yes"), (NBA_ML, "yes")])
    legtable.Refresher(table, fetch, interval_sec=2).refresh_once()
    put(rfqs, rfq("r1", [leg(NFL_SPREAD), leg(NBA_ML)]))
    quoter.step()
    quoter.shutdown()
    assert quoter.take_counts()["leg_no_event"] == 1 and confirm.calls == []
    dlog.close()


def test_would_quotes_are_confirmed_up_to_the_rate_and_stored_in_full(tmp_path):
    clock = [0.0]
    confirm = fixed_pricer(FAIRS)
    quoter, rfqs, refresher, dlog = make_quoter(
        tmp_path, FAIRS, confirm_fn=confirm, clock=lambda: clock[0],
        confirms_per_sec=1, stored_sample_pct=0)
    quoter._table.want([(NFL_SPREAD, "yes"), (NBA_ML, "yes")])
    refresher.refresh_once()
    for rid in ("r1", "r2", "r3"):
        put(rfqs, rfq(rid, [leg(NFL_SPREAD), leg(NBA_ML)]))
    quoter.step()                                     # one confirm slot this second
    clock[0] = 1.0
    put(rfqs, rfq("r4", [leg(NFL_SPREAD), leg(NBA_ML)]))
    quoter.step()
    quoter.shutdown()
    assert [rid for rid, _ in confirm.calls] == ["r1", "r4"]
    counts = quoter.take_counts()
    assert (counts["would_quote"], counts["would_quote_unconfirmed"], counts["stored"]) == (2, 2, 2)
    dlog.flush()
    rows = dlog.query("SELECT rfq_id, decision, confirmed, round(legs_kalshi_ask_product, 4), pricer_ms "
                      "FROM rfq_decisions ORDER BY rfq_id")
    assert rows == [("r1", "would_quote", True, 0.2704, 4.0), ("r4", "would_quote", True, 0.2704, 4.0)]
    assert dlog.query("SELECT count(*), count(kalshi_yes_ask), count(priced_at_ms) FROM rfq_legs") == [(4, 4, 4)]
    dlog.close()


def test_a_confirm_that_declines_is_stored_as_the_decline(tmp_path):
    def confirm(rfq_id, legs):
        return pricer.PriceResult(
            [priced(t, s, FAIRS[t], ok=t != NBA_ML, reason=None if t != NBA_ML else "market_mismatch")
             for t, s in legs], 5, 2.0)
    quoter, rfqs, refresher, dlog = make_quoter(tmp_path, FAIRS, confirm_fn=confirm,
                                                confirms_per_sec=1000, stored_sample_pct=0)
    quoter._table.want([(NFL_SPREAD, "yes"), (NBA_ML, "yes")])
    refresher.refresh_once()
    put(rfqs, rfq("r1", [leg(NFL_SPREAD), leg(NBA_ML)]))
    quoter.step()
    quoter.shutdown()
    dlog.flush()
    assert dlog.query("SELECT decision, reason, confirmed FROM rfq_decisions") == [
        ("decline", "leg_market_mismatch", True)]
    assert quoter.take_counts()["confirm_leg_market_mismatch"] == 1
    dlog.close()


def test_unconfirmed_rows_are_stored_by_a_fixed_sample_of_ids(tmp_path):
    ids = [f"rfq-{n}" for n in range(20_000)]
    share = sum(stored_in_sample(i, 1) for i in ids) / len(ids)
    assert 0.007 < share < 0.013
    assert all(stored_in_sample(i, 100) for i in ids[:50])
    assert not any(stored_in_sample(i, 0) for i in ids[:50])
    assert [stored_in_sample(i, 1) for i in ids[:500]] == [stored_in_sample(i, 1) for i in ids[:500]]
    quoter, rfqs, refresher, dlog = make_quoter(tmp_path, FAIRS, stored_sample_pct=0)
    quoter._table.want([(NFL_SPREAD, "yes"), (NBA_ML, "yes")])
    refresher.refresh_once()
    put(rfqs, rfq("r1", [leg(NFL_SPREAD), leg(NBA_ML)]))
    quoter.step()
    counts = quoter.take_counts()
    assert counts["would_quote_unconfirmed"] == 1 and "stored" not in counts
    dlog.flush()
    assert dlog.query("SELECT count(*) FROM rfq_decisions")[0][0] == 0, "counted, not stored"
    dlog.close()


def test_one_bad_rfq_does_not_stop_the_loop(tmp_path):
    quoter, rfqs, refresher, dlog = make_quoter(tmp_path, FAIRS)
    quoter._table.want([(NFL_SPREAD, "yes"), (NBA_ML, "yes")])
    refresher.refresh_once()
    put(rfqs, {"id": "bad", "mve_selected_legs": [None, 7]})
    put(rfqs, rfq("r1", [leg(NFL_SPREAD), leg(NBA_ML)]))
    assert quoter.step() == 2
    counts = quoter.take_counts()
    assert (counts["handle_crashed"], counts["would_quote_unconfirmed"]) == (1, 1)
    dlog.close()


def test_minute_counts_round_trip_and_old_files_gain_the_new_columns(tmp_path):
    import duckdb
    from kalshi_parlay_quoter import db as dbmod
    path = tmp_path / "old.duckdb"
    con = duckdb.connect(str(path))
    con.execute(dbmod.SCHEMA_SQL.split("-- 2026-11-13")[0])      # the schema before this change
    assert "confirmed" not in [r[0] for r in con.execute("DESCRIBE rfq_decisions").fetchall()]
    con.close()
    dlog = DecisionLog(path)
    minute = datetime(2026, 11, 13, 2, 0, tzinfo=timezone.utc)
    dlog.add_counts(minute, {"kept": 475, "frames_deleted": 0, "would_quote": 3})
    dlog.add(*decide.decide(rfq(), [leg(NFL_SPREAD), leg(NBA_ML)], fixed_pricer(FAIRS)))
    dlog.flush()
    dlog.checkpoint()
    assert dlog.query("SELECT key, n FROM minute_counts ORDER BY key") == [("kept", 475), ("would_quote", 3)]
    assert dlog.query("SELECT confirmed FROM rfq_decisions") == [(False,)]
    dlog.close()


def test_a_failed_confirm_backs_off_and_counts_apart_from_first_pass_errors(tmp_path):
    clock = [0.0]
    calls = []

    def confirm(rfq_id, legs):
        calls.append(rfq_id)
        raise pricer.PricerTimeout("timeout")

    quoter, rfqs, refresher, dlog = make_quoter(tmp_path, FAIRS, confirm_fn=confirm, clock=lambda: clock[0],
                                                confirms_per_sec=1, confirm_backoff_sec=30, stored_sample_pct=0)
    quoter._table.want([(NFL_SPREAD, "yes"), (NBA_ML, "yes")])
    refresher.refresh_once()
    for second, rid in ((0.0, "r1"), (5.0, "r2"), (29.0, "r3"), (31.0, "r4")):
        clock[0] = second
        refresher.refresh_once()
        put(rfqs, rfq(rid, [leg(NFL_SPREAD), leg(NBA_ML)]))
        quoter.step()
        _wait_for(lambda: quoter._confirms_in_flight == 0, "the confirm to finish")
    quoter.shutdown()
    assert calls == ["r1", "r4"], "no confirm for 30 s after one the runner failed"
    counts = quoter.take_counts()
    assert (counts["confirm_pricer_error"], counts["would_quote_unconfirmed"]) == (2, 2)
    assert "pricer_error" not in counts
    dlog.flush()
    assert dlog.query("SELECT rfq_id, reason, reason_detail, confirmed FROM rfq_decisions ORDER BY 1") == [
        ("r1", "pricer_error", "timeout", True), ("r4", "pricer_error", "timeout", True)]
    dlog.close()


def test_queue_wait_is_in_queue_ms_and_an_rfq_too_old_is_not_decided(tmp_path):
    quoter, rfqs, refresher, dlog = make_quoter(tmp_path, FAIRS, max_rfq_age_sec=15)
    quoter._table.want([(NFL_SPREAD, "yes"), (NBA_ML, "yes")])
    refresher.refresh_once()
    put(rfqs, rfq("waited", [leg(NFL_SPREAD), leg(NBA_ML)]), age_sec=2)
    put(rfqs, rfq("gone", [leg(NFL_SPREAD), leg(NBA_ML)]), age_sec=16)
    assert quoter.step() == 2
    counts = quoter.take_counts()
    assert (counts["would_quote_unconfirmed"], counts["expired_in_queue"]) == (1, 1)
    dlog.flush()
    (rid, queue_ms, lag_sec), = dlog.query(
        "SELECT rfq_id, queue_ms, epoch(now()) - epoch(seen_at) FROM rfq_decisions")
    assert rid == "waited" and 2000 <= queue_ms < 4000 and 2 <= lag_sec < 6
    dlog.close()


def test_pause_drops_what_was_waiting_and_counts_it(tmp_path):
    quoter, rfqs, refresher, dlog = make_quoter(tmp_path, FAIRS)
    put(rfqs, rfq("parked", [leg(NFL_SPREAD), leg(NBA_ML)]))
    quoter.step()                                     # parked: no fair yet
    put(rfqs, rfq("queued", [leg(NFL_SPREAD), leg(NBA_ML)]))
    quoter.drop_waiting()
    refresher.refresh_once()
    quoter.step()
    counts = quoter.take_counts()
    assert counts == {"dropped_on_pause": 2, "gauge_pending": 0, "gauge_queue": 0, "gauge_table_legs": 2}
    dlog.close()


def test_a_paused_refresher_asks_the_runner_nothing_and_resumes():
    table = legtable.LegTable(max_age_sec=15, leg_ttl_sec=600)
    keys = [(f"KXNBATOTAL-26NOV08LALBOS-{n}", "yes") for n in range(1, 451)]
    table.want(keys)
    calls = []

    def fetch(chunk):
        calls.append(len(chunk))
        if len(calls) == 1:
            refresher.set_paused(True)       # paused mid-pass
        return pricer.PriceResult([priced(t, s, .5) for t, s in chunk], 9, 1.0)

    refresher = legtable.Refresher(table, fetch, interval_sec=0.05)
    refresher.refresh_once()
    assert calls == [200], "a pause ends the pass at the next chunk"
    refresher.start()
    time.sleep(0.4)
    assert calls == [200], "the running thread asks nothing while paused"
    refresher.set_paused(False)
    deadline = time.monotonic() + 3
    while len(calls) < 3 and time.monotonic() < deadline:
        time.sleep(0.02)
    refresher.stop()
    assert len(calls) >= 3, "and resumes when the pause lifts"


def test_new_legs_are_priced_first_and_between_passes():
    table = legtable.LegTable(max_age_sec=15, leg_ttl_sec=600)
    old = [(f"KXNBATOTAL-26NOV08LALBOS-{n}", "yes") for n in range(1, 301)]
    table.want(old)
    calls = []

    def fetch(chunk):
        calls.append(list(chunk))
        return pricer.PriceResult([priced(t, s, .5) for t, s in chunk], 9, 1.0)

    refresher = legtable.Refresher(table, fetch, interval_sec=2)
    refresher.refresh_once()
    new = [(NFL_SPREAD, "yes"), (NBA_ML, "no")]
    table.want(new)
    assert table.unpriced() == new
    calls.clear()
    assert refresher.price_new_legs() is True
    assert calls == [new] and table.missing(new) == [] and table.unpriced() == []
    table.want([(CBB_TOTAL, "yes")])
    calls.clear()
    refresher.refresh_once()
    assert calls[0][0] == (CBB_TOTAL, "yes"), "a full pass prices the leg with no fair first"
    assert sum(map(len, calls)) == 303


def test_a_leg_the_runner_rejects_is_isolated_and_the_rest_still_price(tmp_path):
    bad = ("KXNBAGAME-26NOV08LALBOS-BAD", "yes")
    clock = [1000.0]
    table = legtable.LegTable(max_age_sec=15, leg_ttl_sec=600, clock=lambda: clock[0])
    good = [(f"KXNBATOTAL-26NOV08LALBOS-{n}", "yes") for n in range(1, 451)]
    keys = good[:100] + [bad] + good[100:]
    table.want(keys)
    calls = []

    def fetch(chunk):
        calls.append(len(chunk))
        if bad in chunk:
            raise pricer.PricerRejected("http_400")
        return pricer.PriceResult([priced(t, s, .5) for t, s in chunk], 9, 1.0)

    refresher = legtable.Refresher(table, fetch, interval_sec=2, clock=lambda: clock[0])
    assert refresher.refresh_once() is True
    assert table.missing(keys) == [], "every good leg priced, the bad one answered as a decline"
    assert len(calls) <= 3 + 2 * 8, "halving, not one call per leg"
    assert refresher.take_stats()["table_legs_rejected"] == 1
    declined = table.price("r1", [bad, good[0]]).legs[0]
    assert (declined["ok"], declined["reason"], declined["detail"]) == (False, "runner_rejected", "http_400")
    calls.clear()
    refresher.refresh_once()
    assert calls == [200, 200, 50], "a rejected leg is not asked again"
    other = (NFL_SPREAD, "yes")
    table.want([other])
    table.store([other], [priced(*other, .5)], 9, clock[0])
    row, _ = decide.decide(rfq("r1", [leg(bad[0]), leg(NFL_SPREAD)]), [leg(bad[0]), leg(NFL_SPREAD)], table.price)
    assert (row["decision"], row["reason"]) == ("decline", "leg_runner_rejected")
    clock[0] += legtable.REJECT_RETRY_SEC + 1
    table.want(keys)
    calls.clear()
    refresher.refresh_once()
    assert sum(calls) > 450, "after 5 minutes the rejected leg is asked again"
    assert refresher.take_stats()["table_legs_rejected"] == 1, "and rejected again"


def test_4xx_and_misaligned_answers_are_rejections_and_5xx_is_not():
    def http(code):
        return lambda r, timeout: (_ for _ in ()).throw(urllib.error.HTTPError(r.full_url, code, "x", {}, io.BytesIO(b"")))
    with pytest.raises(pricer.PricerRejected, match="http_400"):
        pricer.fairs("http://x", [(NFL_SPREAD, "yes")], timeout=1, opener=http(400))
    for code in (403, 404, 500):     # about the runner, not these legs
        with pytest.raises(pricer.PricerError, match=f"http_{code}") as raised:
            pricer.fairs("http://x", [(NFL_SPREAD, "yes")], timeout=1, opener=http(code))
        assert not isinstance(raised.value, pricer.PricerRejected)


def test_minute_counts_add_up_across_a_restart_inside_the_minute(tmp_path):
    minute = datetime(2026, 11, 13, 2, 0, tzinfo=timezone.utc)
    for kept, gauge in ((100, 7), (40, 0)):
        dlog = DecisionLog(tmp_path / "q.duckdb")
        dlog.add_counts(minute, {"kept": kept, "gauge_pending": gauge})
        dlog.close()
    dlog = DecisionLog(tmp_path / "q.duckdb")
    assert dlog.query("SELECT key, n FROM minute_counts ORDER BY key") == [("gauge_pending", 0), ("kept", 140)]
    dlog.close()
