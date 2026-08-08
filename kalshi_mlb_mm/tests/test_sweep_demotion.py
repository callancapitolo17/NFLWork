"""Issue #57 — demote the sweep: the quote path must run entirely without
sweep rows, and the post-fill cooldown must clear via a TARGETED on-demand
fetch instead of waiting for a full-slate sweep pass.

Written FIRST, failing on pre-#57 code, per repo rule. What each block pins:

1. Cooldown: `in_cooldown_awaiting_refresh` keys on the OnDemandEngine's
   completed post-fill landing for every game the combo touches — sweep rows
   can neither clear it nor are they needed; the gate lazily re-triggers the
   targeted fetch (restart / store-prune self-heal), and the confirm tick
   fires it eagerly the moment a fill is recorded.
2. `same_price_block` compares the filled fair against the LIVE fair the
   targeted fetch produced (with zero sweep rows in sight).
3. Discovery tick has no sweep-freshness top gate: an empty `_SGP_ODDS`
   prices RFQs normally via the live feed.
4. Risk sweep: no `books_stale` mass-cancel; the constituent-jump breaker
   runs (and fires) with zero sweep rows; the "can't trust our fair" pull is
   re-keyed on #38 live-fetch health (`books_unhealthy` via
   BookHealthAlerter.consensus_dark()).
5. `_current_consensus_fair` (book-drift breaker) prices from live engine
   fairs, not sweep rows.
6. Zombie-gate proof: retired symbols are gone from the source (repo rule:
   delete dead code, git is the archive — mirrors #56's deleted-knob test).
7. Config: sweep slowed to 300s; #37 alerts count only the on_demand path
   so a slow/failing sweep can never trip them (path-filter behavior itself
   is covered in kalshi_common/tests/test_book_health_alerts.py).

Scaffolding style follows test_live_pricing.py; clock control is injected
(engine ages / now_fn), never slept.
"""
import importlib
import inspect
import json
from datetime import datetime, timedelta, timezone
from pathlib import Path

import pandas as pd

from kalshi_common import legset
from kalshi_common.book_health import BookHealthAlerter
from mlb_sgp._shared import GameRef, OnDemandBookResult

from kalshi_mlb_mm.tests.conftest import leg

# Same-game spread+total (grid route); LAA is home in 25JUN271905TEXLAA.
GRID_LEGS = [
    leg("KXMLBSPREAD-25JUN271905TEXLAA-LAA2", "yes"),
    leg("KXMLBTOTAL-25JUN271905TEXLAA-9", "yes"),
]
GREF = GameRef(game_id="game1", home_team="Los Angeles Angels",
               away_team="Texas Rangers", commence_time=None)
TICKER = "COMBO-SD"
LIVE_FAIRS = {"draftkings": 0.30, "fanduel": 0.32}   # consensus median 0.31
LIVE_FAIR = 0.31


def _hash():
    return legset.leg_set_hash(legset.parse_legs(GRID_LEGS))


class Eng:
    """Duck-typed OnDemandEngine: canned live fairs + a configurable
    completed-landing age (the #57 accessor), records ensure/refetch calls."""

    def __init__(self, fairs=None, completed_age=None):
        self.fairs = dict(fairs) if fairs else None
        self.completed_age = completed_age    # sec since completed landing
        self.ensure_calls = []
        self.refetch_calls = []

    def lookup(self, h):
        return dict(self.fairs) if self.fairs else None

    def lookup_results(self, h):
        if not self.fairs:
            return None
        return {b: OnDemandBookResult(book=b, fair=f, route="partition",
                                      n_cells_priced=4, latency_sec=1.0)
                for b, f in self.fairs.items()}

    def landed_at(self, h):
        return 1000.0 if self.fairs else None

    def landed_empty(self, h):
        return False

    def result_age_sec(self, h):
        return 2.0 if self.fairs else None

    def dropped_books(self, h):
        return ()

    def ensure_fetch(self, h, game, legs):
        self.ensure_calls.append((h, game, tuple(legs)))
        return True

    def refetch_now(self, jobs, deadline_sec):
        self.refetch_calls.append((list(jobs), deadline_sec))
        return self.fairs is not None

    def completed_fetch_age_sec(self, h):
        return self.completed_age


class Src:
    def poll(self):
        return [{"id": "r-sd", "market_ticker": TICKER, "contracts": 3}]

    def get_market(self, t):
        return {}


class NoQuoteGW:
    def submit_quote(self, *a):
        raise AssertionError("must not submit in these tests")


def _env(monkeypatch, tmp_path, db_name, *, engine, sgp_odds=None,
         fill_fair=None, filled_secs_ago=120):
    """Sweep-free discovery/confirm env. `sgp_odds` defaults to None — the
    whole point of #57 is that the quote path never needs it. `fill_fair`
    seeds a prior fill + expired cooldown floor on TICKER."""
    import kalshi_mlb_mm.config as cfg
    import kalshi_mlb_mm.db as db
    import kalshi_mlb_mm.risk as risk
    from kalshi_mlb_mm import main

    monkeypatch.setattr(cfg, "DB_PATH", tmp_path / db_name)
    monkeypatch.setattr(cfg, "KILL_FILE", tmp_path / ".kill")
    importlib.reload(db)
    db.init_database()

    monkeypatch.setattr(main, "_SGP_ODDS", sgp_odds)
    monkeypatch.setattr(main, "_today_fills", lambda: [])
    monkeypatch.setattr(main, "_today_fills_by_game", lambda: [])
    monkeypatch.setattr(main, "_resolve_game_for_legs", lambda gl: "game1")
    monkeypatch.setattr(main, "_commence_time", lambda gid: None)
    monkeypatch.setattr(main, "_game_ref", lambda gid: GREF)
    monkeypatch.setattr(risk, "tipoff_ok", lambda ct, min_: True)
    monkeypatch.setattr(main, "_ENGINE", engine)
    monkeypatch.setattr(main, "_SCOPE_CACHE", {TICKER: (True, None, GRID_LEGS)})
    monkeypatch.setattr(main, "_OD_RESULT_EMITTED", {})
    monkeypatch.setattr(main, "_PREV_BOOK_FAIR", {})
    emitted = []
    monkeypatch.setattr(main.research, "emit",
                        lambda ev, **kw: emitted.append((ev, kw)))

    filled_at = datetime.now(timezone.utc) - timedelta(seconds=filled_secs_ago)
    if fill_fair is not None:
        with db.connect() as con:
            con.execute(
                "INSERT INTO fills (fill_id, quote_id, rfq_id, "
                "combo_market_ticker, game_id, side_held, contracts, price, "
                "fee, model_fair_at_quote, book_fair_at_quote, "
                "blended_fair_at_quote, fair_at_confirm, realized_pnl, "
                "filled_at, reconciled) VALUES (?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?)",
                ["fill-sd", "q-sd-old", "r-sd-prior", TICKER, "game1", "yes",
                 1, 0.25, 0.0, None, fill_fair, fill_fair, fill_fair, None,
                 filled_at, True])
            con.execute("INSERT OR REPLACE INTO fill_games (fill_id, game_id) "
                        "VALUES (?, ?)", ["fill-sd", "game1"])
            con.execute(
                "INSERT OR REPLACE INTO combo_cooldown "
                "(combo_market_ticker, cooled_until) VALUES (?, ?)",
                [TICKER, datetime.now(timezone.utc) - timedelta(seconds=60)])
    return main, db, emitted, filled_at


def _last_decision(db, rfq_id="r-sd"):
    with db.connect(read_only=True) as con:
        return con.execute(
            "SELECT decision, reason FROM quote_decisions "
            "WHERE rfq_id=? ORDER BY observed_at DESC LIMIT 1",
            [rfq_id]).fetchone()


# --------------------------------------------------------------------------- #
# 1. Cooldown clears via the targeted live fetch — sweep-free                 #
# --------------------------------------------------------------------------- #

def test_awaiting_refresh_until_targeted_fetch_lands(monkeypatch, tmp_path):
    """Floor expired, live fairs are even being served (quorum partial), but
    no COMPLETED post-fill flight exists → stay cooled AND lazily re-trigger
    the targeted fetch so the gate clears itself without any sweep."""
    eng = Eng(fairs=LIVE_FAIRS, completed_age=None)
    main, db, _, _ = _env(monkeypatch, tmp_path, "sd1.duckdb", engine=eng,
                          fill_fair=0.55)
    main._discovery_tick(Src(), NoQuoteGW(), dry_run=True)
    assert _last_decision(db) == ("skipped", "in_cooldown_awaiting_refresh")
    assert eng.ensure_calls and eng.ensure_calls[0][0] == _hash(), \
        "the cooldown gate must (re-)trigger the targeted fetch itself"
    assert eng.ensure_calls[0][1] == GREF


def test_pre_fill_landing_does_not_clear_cooldown(monkeypatch, tmp_path):
    """A flight that completed BEFORE the fill is exactly the data that
    priced the picked-off fill — it must not clear the gate."""
    eng = Eng(fairs=LIVE_FAIRS, completed_age=300.0)     # fill was 120s ago
    main, db, _, _ = _env(monkeypatch, tmp_path, "sd2.duckdb", engine=eng,
                          fill_fair=0.55)
    main._discovery_tick(Src(), NoQuoteGW(), dry_run=True)
    assert _last_decision(db) == ("skipped", "in_cooldown_awaiting_refresh")
    assert eng.ensure_calls, "pre-fill landing must re-trigger a fresh fetch"


def test_post_fill_landing_and_moved_fair_quotes(monkeypatch, tmp_path):
    """Targeted fetch landed after the fill and the live fair moved off the
    filled fair → quote, with ZERO sweep rows anywhere."""
    eng = Eng(fairs=LIVE_FAIRS, completed_age=10.0)
    main, db, _, _ = _env(monkeypatch, tmp_path, "sd3.duckdb", engine=eng,
                          fill_fair=0.55)
    main._discovery_tick(Src(), NoQuoteGW(), dry_run=True)
    assert _last_decision(db) == ("dry_run_quote", None)


def test_same_price_block_compares_live_fair(monkeypatch, tmp_path):
    """Post-fill fetch landed but its LIVE consensus still sits on the
    filled fair → same_price_block. Sweep rows absent, so the only possible
    comparison source is the targeted fetch's fair."""
    eng = Eng(fairs={"draftkings": LIVE_FAIR, "fanduel": LIVE_FAIR},
              completed_age=10.0)
    main, db, _, _ = _env(monkeypatch, tmp_path, "sd4.duckdb", engine=eng,
                          fill_fair=LIVE_FAIR)
    main._discovery_tick(Src(), NoQuoteGW(), dry_run=True)
    assert _last_decision(db) == ("skipped", "same_price_block")


def test_sweep_rows_alone_cannot_clear_cooldown(monkeypatch, tmp_path):
    """The old clearing condition — post-fill sweep rows in _SGP_ODDS — must
    no longer open the gate: only a completed targeted fetch counts."""
    post_fill_rows = pd.DataFrame(
        {"game_id": ["game1"], "combo": ["c"], "period": ["FG"],
         "bookmaker": ["dk"], "sgp_decimal": [2.0],
         "fetch_time": [pd.Timestamp(datetime.now(timezone.utc))],
         "spread_line": [-1.5], "total_line": [8.5]})
    eng = Eng(fairs=LIVE_FAIRS, completed_age=None)      # no targeted landing
    main, db, _, _ = _env(monkeypatch, tmp_path, "sd5.duckdb", engine=eng,
                          sgp_odds=post_fill_rows, fill_fair=0.55)
    main._discovery_tick(Src(), NoQuoteGW(), dry_run=True)
    assert _last_decision(db) == ("skipped", "in_cooldown_awaiting_refresh")


def test_confirm_fill_triggers_targeted_fetch(monkeypatch, tmp_path):
    """The moment a fill is recorded, the confirm tick fires the targeted
    post-fill fetch for the combo's game leg set (research-tagged), so the
    cooldown usually clears right at the 60s floor. Also pins that the whole
    confirm path (live re-fetch last look included) runs with _SGP_ODDS=None."""
    eng = Eng(fairs=LIVE_FAIRS)
    main, db, emitted, _ = _env(monkeypatch, tmp_path, "sd6.duckdb", engine=eng)
    now = datetime.now(timezone.utc)
    snapshot = {"L": {"yes_bid": 0.5, "yes_ask": 0.52}}
    with db.connect() as con:
        con.execute(
            "INSERT INTO live_quotes (quote_id, rfq_id, combo_market_ticker, "
            "game_id, yes_bid, no_bid, model_fair, book_fair, blended_fair, "
            "status, submitted_at, closed_at, leg_prices_json) "
            "VALUES (?,?,?,?,?,?,?,?,?,?,?,?,?)",
            ["q-sd", "r-sd", TICKER, "game1", 0.20, 0.75, None, LIVE_FAIR,
             LIVE_FAIR, "open", now, None, json.dumps(snapshot)])
        con.execute(
            "INSERT OR REPLACE INTO seen_rfqs (rfq_id, market_ticker, in_scope, "
            "game_id, legs_json, first_seen_at, last_decision, creator_id) "
            "VALUES (?,?,?,?,?,?,?,?)",
            ["r-sd", TICKER, True, "game1", json.dumps(GRID_LEGS), now,
             "quoted", ""])
    monkeypatch.setattr(main, "_leg_market_prices", lambda legs: snapshot)
    monkeypatch.setattr(main.notify, "fill", lambda *a, **k: None)
    monkeypatch.setattr(
        main.auth_client, "api",
        lambda *a, **k: (200, {"quote": {"status": "accepted",
                                         "accepted_side": "no",
                                         "contracts_fp": "1.00"}}, None))

    class GW:
        def confirm(self, qid):
            return True

    main._confirm_tick(GW(), dry_run=False)
    with db.connect(read_only=True) as con:
        n_fills = con.execute("SELECT COUNT(*) FROM fills").fetchone()[0]
    assert n_fills == 1, "harness must reach a recorded fill"
    assert eng.ensure_calls and eng.ensure_calls[0][0] == _hash(), \
        "recording a fill must immediately queue the targeted post-fill fetch"
    post_fill_reqs = [kw for ev, kw in emitted if ev == "on_demand_requested"
                      and kw["payload"].get("trigger") == "post_fill"]
    assert post_fill_reqs, "post-fill fetch must be research-tagged"


# --------------------------------------------------------------------------- #
# 3. Discovery tick has no sweep-freshness top gate                           #
# --------------------------------------------------------------------------- #

def test_discovery_tick_prices_with_empty_sweep_frame(monkeypatch, tmp_path):
    """A fresh, never-filled combo prices from the live feed while _SGP_ODDS
    is an EMPTY frame (mid-gap between two 300s sweep passes)."""
    eng = Eng(fairs=LIVE_FAIRS)
    main, db, emitted, _ = _env(monkeypatch, tmp_path, "sd7.duckdb", engine=eng,
                                sgp_odds=pd.DataFrame())
    main._discovery_tick(Src(), NoQuoteGW(), dry_run=True)
    assert _last_decision(db) == ("dry_run_quote", None)
    payloads = [kw["payload"] for ev, kw in emitted if ev == "quote_priced"]
    assert payloads and abs(payloads[-1]["blended_fair"] - LIVE_FAIR) < 1e-9


# --------------------------------------------------------------------------- #
# 4. Risk sweep: no books_stale; jump breaker sweep-free; health-dark pull    #
# --------------------------------------------------------------------------- #

_SPREAD_T = "KXMLBSPREAD-25JUN271905TEXLAA-LAA2"
_TOTAL_T = "KXMLBTOTAL-25JUN271905TEXLAA-9"
_BASELINE = {_SPREAD_T: {"yes_bid": 0.45, "yes_ask": 0.47},
             _TOTAL_T: {"yes_bid": 0.52, "yes_ask": 0.54}}
_QUIET = {_SPREAD_T: {"yes_bid": 0.455, "yes_ask": 0.475},
          _TOTAL_T: {"yes_bid": 0.52, "yes_ask": 0.54}}
_JUMPED = {_SPREAD_T: {"yes_bid": 0.55, "yes_ask": 0.57},
           _TOTAL_T: {"yes_bid": 0.52, "yes_ask": 0.54}}


def _sweep_env(monkeypatch, tmp_path, db_name, *, current):
    """One open quote with a placement baseline; _SGP_ODDS is None (the #57
    steady state between slow sweep passes); constituent poll stubbed."""
    import kalshi_mlb_mm.config as cfg
    import kalshi_mlb_mm.db as db
    from kalshi_mlb_mm import main, singles

    monkeypatch.setattr(cfg, "DB_PATH", tmp_path / db_name)
    monkeypatch.setattr(cfg, "KILL_FILE", tmp_path / ".kill")
    importlib.reload(db)
    db.init_database()
    monkeypatch.setattr(main, "_SGP_ODDS", None)
    monkeypatch.setattr(main, "_ENGINE", None)
    monkeypatch.setattr(main, "_resolve_game_for_legs", lambda gl: "game1")
    monkeypatch.setattr(main, "_commence_time",
                        lambda gid: datetime(2099, 1, 1, tzinfo=timezone.utc))
    polled = []

    def fake_poll(tickers, budget_sec=None):
        polled.append(tuple(tickers))
        return current

    monkeypatch.setattr(singles, "fetch_market_prices", fake_poll)
    now = datetime.now(timezone.utc)
    with db.connect() as con:
        con.execute(
            "INSERT INTO live_quotes (quote_id, rfq_id, combo_market_ticker, "
            "game_id, yes_bid, no_bid, model_fair, book_fair, blended_fair, "
            "status, submitted_at, closed_at, leg_prices_json) "
            "VALUES (?,?,?,?,?,?,?,?,?,?,?,?,?)",
            ["q1", "r-q1", "COMBO-q1", "game1", 0.20, 0.72, None, 0.25, 0.25,
             "open", now, None, json.dumps(_BASELINE)])
        con.execute(
            "INSERT OR REPLACE INTO seen_rfqs (rfq_id, market_ticker, in_scope, "
            "game_id, legs_json, first_seen_at, last_decision, creator_id) "
            "VALUES (?,?,?,?,?,?,?,?)",
            ["r-q1", "COMBO-q1", True, "game1", json.dumps(GRID_LEGS), now,
             "quoted", ""])
    return main, db, polled


def _quote_state(db, qid="q1"):
    with db.connect(read_only=True) as con:
        status = con.execute(
            "SELECT status FROM live_quotes WHERE quote_id=?", [qid]).fetchone()[0]
        row = con.execute(
            "SELECT decision, reason FROM quote_decisions WHERE quote_id=? "
            "ORDER BY observed_at DESC LIMIT 1", [qid]).fetchone()
    return status, row


class _GW:
    def __init__(self):
        self.cancelled = []

    def cancel(self, qid):
        self.cancelled.append(qid)
        return True


def test_no_books_stale_mass_cancel_without_sweep_rows(monkeypatch, tmp_path):
    """_SGP_ODDS=None is the NORMAL state under a 300s sweep — quiet
    constituents must leave the resting quote alone (pre-#57 this cancelled
    everything with books_stale ~40% of every cycle)."""
    main, db, polled = _sweep_env(monkeypatch, tmp_path, "sd8.duckdb",
                                  current=_QUIET)
    gw = _GW()
    main._risk_sweep_tick(gw)
    status, _ = _quote_state(db)
    assert status == "open" and gw.cancelled == []
    assert polled, "constituent poll must run without sweep rows — it is " \
                   "the real-time protection that replaces books_stale"


def test_constituent_jump_fires_without_sweep_rows(monkeypatch, tmp_path):
    """The jump breaker is the fast pull signal — it must work sweep-free."""
    main, db, _ = _sweep_env(monkeypatch, tmp_path, "sd9.duckdb",
                             current=_JUMPED)
    gw = _GW()
    main._risk_sweep_tick(gw)
    status, decision = _quote_state(db)
    assert status == "cancelled" and gw.cancelled == ["q1"]
    assert decision == ("sweep_cancel", "constituent_jump")


def _dark_alerter():
    """Real alerter driven dark: 2 known books, both on a 3-failure streak
    on the on_demand path → healthy < MIN_AGREEING_BOOKS."""
    a = BookHealthAlerter(label="T", streak_threshold=3, min_healthy_books=2,
                          paths=("on_demand",), notifier=lambda t, m: None,
                          notify_async=False)
    for book in ("draftkings", "fanduel"):
        for _ in range(3):
            a.observe(book=book, path="on_demand", outcome="transport_error")
    return a


def _healthy_alerter():
    a = BookHealthAlerter(label="T", streak_threshold=3, min_healthy_books=2,
                          paths=("on_demand",), notifier=lambda t, m: None,
                          notify_async=False)
    for book in ("draftkings", "fanduel"):
        a.observe(book=book, path="on_demand", outcome="ok")
    return a


def test_health_dark_pulls_open_quotes(monkeypatch, tmp_path):
    """#38 live-fetch health replaces data-age as the 'we can no longer trust
    our fair' pull: consensus-dark cancels every resting quote."""
    main, db, _ = _sweep_env(monkeypatch, tmp_path, "sd10.duckdb",
                             current=_QUIET)
    gw = _GW()
    main._risk_sweep_tick(gw, book_health=_dark_alerter())
    status, decision = _quote_state(db)
    assert status == "cancelled" and gw.cancelled == ["q1"]
    assert decision == ("sweep_cancel", "books_unhealthy")


def test_health_ok_no_pull(monkeypatch, tmp_path):
    main, db, _ = _sweep_env(monkeypatch, tmp_path, "sd11.duckdb",
                             current=_QUIET)
    gw = _GW()
    main._risk_sweep_tick(gw, book_health=_healthy_alerter())
    status, _ = _quote_state(db)
    assert status == "open" and gw.cancelled == []


def test_book_health_param_is_keyword_only():
    from kalshi_mlb_mm import main
    param = inspect.signature(main._risk_sweep_tick).parameters["book_health"]
    assert param.kind is inspect.Parameter.KEYWORD_ONLY
    assert param.default is None, "fail-safe default: no alerter, no pull"


def test_consensus_dark_accessor():
    """The risk sweep's read of Rule B state: live computation, no dispatch
    needed, first-observation precondition intact."""
    assert _dark_alerter().consensus_dark() is True
    assert _healthy_alerter().consensus_dark() is False
    fresh = BookHealthAlerter(label="T", streak_threshold=1,
                              min_healthy_books=2, paths=("on_demand",),
                              notifier=lambda t, m: None, notify_async=False)
    assert fresh.consensus_dark() is False, "no observations → unknown, not dark"
    fresh.observe(book="draftkings", path="on_demand", outcome="error")
    assert fresh.consensus_dark() is False, \
        "one known book < min_healthy_books → precondition blocks a false dark"


# --------------------------------------------------------------------------- #
# 5. Book-drift breaker prices from live fairs                                #
# --------------------------------------------------------------------------- #

def test_current_consensus_fair_uses_live_fairs(monkeypatch):
    from kalshi_mlb_mm import main
    monkeypatch.setattr(main, "_SGP_ODDS", None)
    monkeypatch.setattr(main, "_ENGINE", Eng(fairs=LIVE_FAIRS))
    monkeypatch.setattr(main, "_resolve_game_for_legs", lambda gl: "game1")
    fair = main._current_consensus_fair(json.dumps(GRID_LEGS))
    assert fair is not None and abs(fair - LIVE_FAIR) < 1e-9, \
        "drift breaker must see live fairs or it goes blind between sweeps"


# --------------------------------------------------------------------------- #
# 6. Zombie-gate proof + engine accessor unit                                 #
# --------------------------------------------------------------------------- #

def test_retired_gates_are_deleted_not_flagged():
    from kalshi_mlb_mm import config, main
    src_main = Path(inspect.getsourcefile(main)).read_text()
    assert "books_stale" not in src_main, \
        "books_stale must be DELETED (repo rule), not left as a zombie gate"
    assert "_post_fill_refresh_landed" not in src_main, \
        "the sweep-row cooldown check must be deleted with its gate"
    src_cfg = Path(inspect.getsourcefile(config)).read_text()
    assert "MAX_BOOK_STALENESS_SEC" not in src_cfg, \
        "the sweep-staleness knob has no consumer left post-#57"


def test_engine_completed_fetch_age_sec():
    """Real-engine unit for the #57 accessor: age of the last COMPLETED
    landing with >=1 book result; None for in-flight, never-fetched, or
    zero-book flights. Clock injected via now_fn."""
    from kalshi_mlb_mm.on_demand import OnDemandEngine
    now = [100.0]
    results = {"draftkings": OnDemandBookResult(
        book="draftkings", fair=0.3, route="partition", n_cells_priced=4,
        latency_sec=0.1)}

    class Svc:
        books = ("draftkings",)
        on_demand_deadline_sec = 10.0

        def price_on_demand(self, book, game, legs):
            return results.get(book)

    eng = OnDemandEngine(Svc(), now_fn=lambda: now[0], autostart=False)
    h = _hash()
    assert eng.completed_fetch_age_sec(h) is None          # never fetched
    eng.ensure_fetch(h, GREF, legset.parse_legs(GRID_LEGS))
    assert eng.completed_fetch_age_sec(h) is None          # in flight
    eng._drain_once()
    assert eng.completed_fetch_age_sec(h) == 0.0           # just completed
    now[0] += 42.0
    assert eng.completed_fetch_age_sec(h) == 42.0
    # Zero-book completion (every book declined) is NOT a priced landing.
    results.clear()
    eng.ensure_fetch(h, GREF, legset.parse_legs(GRID_LEGS))
    eng._drain_once()
    assert eng.completed_fetch_age_sec(h) is None


# --------------------------------------------------------------------------- #
# 7. Config: slow sweep, alerts keyed off it                                  #
# --------------------------------------------------------------------------- #

def test_sweep_cadence_default_slowed_to_300():
    from kalshi_mlb_mm import config
    assert config.SGP_REFRESH_SEC == 300
    assert config.STRUCTURE_WARM_SEC < config.SGP_REFRESH_SEC, \
        "#50 structure warming must stay on its own HOT cadence"


def test_book_alerts_count_only_the_on_demand_path():
    """The #57 flip BOOK_ALERT_PATHS promised in config: a slow (or dead)
    sweep can never trip #37's streak alerts once only on_demand counts."""
    from kalshi_mlb_mm import config
    assert config.BOOK_ALERT_PATHS == ("on_demand",)
