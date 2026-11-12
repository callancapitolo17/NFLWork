"""Dry-run parlay quoter (#283): scope rule, pricer client, decision math,
decision log, and the loop's one-decision-per-RFQ accounting."""
import io
import json
import urllib.error
from datetime import datetime, timezone

import pytest

from kalshi_parlay_quoter import decide, pricer, scope
from kalshi_parlay_quoter.db import DecisionLog
from kalshi_parlay_quoter.main import Quoter
from kalshi_mlb_mm.ws_rfq_source import _door_verdict


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


def test_ws_door_widened_without_changing_mlb_default():
    m = {"mve_selected_legs": [leg(NFL_SPREAD), leg(MLB_TOTAL)],
         "target_cost_dollars": "500"}
    assert _door_verdict(m, max_legs=3, min_cost_usd=250) == "non_mlb"
    assert _door_verdict(m, max_legs=3, min_cost_usd=250,
                         leg_prefixes=scope.LEAGUE_PREFIXES) == "keep"
    m["mve_selected_legs"].append(leg("KXNHLGAME-26NOV08BOSTOR-BOS"))
    assert _door_verdict(m, max_legs=10**6, min_cost_usd=0,
                         leg_prefixes=scope.LEAGUE_PREFIXES) == "non_mlb"


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


# -- decision log + loop -------------------------------------------------------

class FakeSource:
    def __init__(self, rfqs, markets=None):
        self.rfqs = rfqs
        self.markets = markets or {}
        self.ingestion_drops = {"non_mlb": 0, "too_many_legs": 0, "below_floor": 0}
        self.malformed_count = 0
        self._mode = "ws"

    def poll(self):
        return list(self.rfqs)

    def get_market(self, ticker):
        return self.markets.get(ticker)


def test_every_rfq_gets_one_decision_row(tmp_path):
    dlog = DecisionLog(tmp_path / "q.duckdb")
    legless = rfq("r3")
    src = FakeSource([
        rfq("r1", [leg(NFL_SPREAD), leg(NBA_ML)]),
        rfq("r2", [leg(NFL_SPREAD), leg(NFL_TOTAL_SAME)]),
        legless,
        rfq("r4", [leg("KXNHLGAME-26NOV08BOSTOR-BOS"), leg(NBA_ML)]),  # REST fallback leak
        rfq("r5"),  # legless, market lookup fails
    ], markets={legless["market_ticker"]: {
        "mve_selected_legs": [leg(CBB_TOTAL), leg(CFB_SPREAD)]}})
    q = Quoter(src, dlog, price_fn=fixed_pricer(
        {NFL_SPREAD: 0.5, NBA_ML: 0.4, CBB_TOTAL: 0.5, CFB_SPREAD: 0.5}), workers=1)
    assert q.tick() == 5
    assert q.tick() == 0          # already decided: no second row
    q.shutdown()
    src.ingestion_drops["non_mlb"] = 7
    src.ingestion_drops["below_floor"] = 3
    minute = datetime(2026, 11, 7, 18, 0, tzinfo=timezone.utc)
    dlog.add_door(q.door_row(minute))
    dlog.flush()
    rows = dict((r[0], (r[1], r[2])) for r in dlog.query(
        "SELECT rfq_id, decision, reason FROM rfq_decisions"))
    assert rows == {"r1": ("would_quote", None), "r2": ("decline", "same_game"),
                    "r3": ("would_quote", None), "r5": ("decline", "no_legs")}
    # Leg rows only for RFQs that reached the pricer (r1, r3).
    assert dlog.query("SELECT count(*) FROM rfq_legs")[0][0] == 4
    assert dlog.query("SELECT reason_detail FROM rfq_decisions WHERE rfq_id='r5'"
                      )[0][0] == "market_fetch_failed"
    door = dlog.query("SELECT logged, other_league, below_floor, feed_mode "
                      "FROM door_counts")
    assert door == [(4, 8, 3, "ws")]
    assert dlog.query("SELECT count(*) FROM rfq_decisions WHERE queue_ms >= 0")[0][0] == 4
    line = dlog.query("SELECT line->>'league' FROM rfq_legs "
                      "WHERE rfq_id='r1' AND leg_idx=0")[0][0]
    assert line == "NFL"
    dlog.close()


def test_caching_pricer_prices_each_leg_once_per_ttl():
    t = [0.0]
    inner = fixed_pricer({NFL_SPREAD: 0.5, NBA_ML: 0.4, CBB_TOTAL: 0.3})
    cp = pricer.CachingPricer(inner, ttl_sec=3, backoff_sec=5, clock=lambda: t[0])
    a = cp("r1", [(NFL_SPREAD, "yes"), (NBA_ML, "no")])
    b = cp("r2", [(NBA_ML, "no"), (CBB_TOTAL, "yes")])
    assert [l["fair"] for l in a.legs] == [0.5, 0.4]
    assert [l["market_ticker"] for l in b.legs] == [NBA_ML, CBB_TOTAL]
    assert inner.calls[1] == ("r2", [(CBB_TOTAL, "yes")])   # only the new leg
    cp("r3", [(NFL_SPREAD, "yes"), (NBA_ML, "no")])
    assert len(inner.calls) == 2 and cp.cache_hits == 3
    t[0] = 3.5                                              # expired
    cp("r4", [(NFL_SPREAD, "yes")])
    assert inner.calls[-1] == ("r4", [(NFL_SPREAD, "yes")])


def test_caching_pricer_backs_off_after_failure():
    t = [0.0]
    calls = []

    def failing(rid, legs):
        calls.append(rid)
        raise pricer.PricerError("http_404")
    cp = pricer.CachingPricer(failing, ttl_sec=3, backoff_sec=5, clock=lambda: t[0])
    with pytest.raises(pricer.PricerError, match="http_404"):
        cp("r1", [(NFL_SPREAD, "yes")])
    with pytest.raises(pricer.PricerError, match="backoff after http_404"):
        cp("r2", [(NFL_SPREAD, "yes")])
    t[0] = 6
    with pytest.raises(pricer.PricerError, match="http_404"):
        cp("r3", [(NFL_SPREAD, "yes")])
    assert calls == ["r1", "r3"]


def test_full_backlog_logs_pricer_backlog_without_pricing(tmp_path):
    dlog = DecisionLog(tmp_path / "q.duckdb")
    fn = fixed_pricer({NFL_SPREAD: 0.5, NBA_ML: 0.4})
    src = FakeSource([rfq("r1", [leg(NFL_SPREAD), leg(NBA_ML)]),
                      rfq("r2", [leg(NFL_SPREAD), leg(NFL_TOTAL_SAME)])])
    q = Quoter(src, dlog, price_fn=fn, workers=1, max_backlog=0)
    assert q.tick() == 2
    q.shutdown()
    dlog.flush()
    rows = dict(dlog.query("SELECT rfq_id, reason FROM rfq_decisions"))
    assert rows == {"r1": "pricer_backlog", "r2": "same_game"}
    assert fn.calls == [] and q.door_row(datetime.now(timezone.utc))["backlog"] == 0
    dlog.close()


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


def test_contract_sized_rfqs_meet_the_floor(tmp_path):
    dlog = DecisionLog(tmp_path / "q.duckdb")
    small = rfq("small", [leg(NFL_SPREAD), leg(NBA_ML)], cost="0")
    small["contracts_fp"] = "20.00"           # pays at most $20 < $50
    big = rfq("big", [leg(NFL_SPREAD), leg(NBA_ML)], cost="0")
    big["contracts_fp"] = "205.00"
    q = Quoter(FakeSource([small, big]), dlog, workers=1, min_cost_usd=50,
               price_fn=fixed_pricer({NFL_SPREAD: .5, NBA_ML: .4}))
    assert q.tick() == 1
    q.shutdown()
    dlog.flush()
    assert [r[0] for r in dlog.query("SELECT rfq_id FROM rfq_decisions")] == ["big"]
    assert q.door_row(datetime.now(timezone.utc))["below_floor"] == 1
    dlog.close()


def test_timeout_does_not_back_off():
    t = [0.0]
    calls = []

    def slow(rid, legs):
        calls.append(rid)
        raise pricer.PricerTimeout("timeout")
    cp = pricer.CachingPricer(slow, ttl_sec=10, backoff_sec=5, clock=lambda: t[0])
    for rid in ("r1", "r2"):
        with pytest.raises(pricer.PricerTimeout):
            cp(rid, [(NFL_SPREAD, "yes")])
    assert calls == ["r1", "r2"]


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
