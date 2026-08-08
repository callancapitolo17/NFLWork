"""Issue #23 items 1-2: constituent-jump circuit breaker in the risk sweep.

Our books refresh every ~150-165s; a combo's constituent Kalshi singles trade
in real time. If one moves past CONSTITUENT_JUMP_THRESHOLD after we quoted,
the market has left our resting price behind and we pull every quote touching
that game.

The placement baseline is live_quotes.leg_prices_json — the raw Kalshi odds
#17 already snapshots at submission — so this needs no schema of its own.
Scaffolding follows test_confirm_singles_veto.py.
"""
import importlib
import json
from datetime import datetime, timezone

import pandas as pd

from kalshi_common.leg_types import SPREAD_TOTAL_FAMILY
from kalshi_mlb_mm.tests.conftest import leg


EVT = "25JUN271905TEXLAA"
SPREAD_T = "KXMLBSPREAD-25JUN271905TEXLAA-LAA2"
TOTAL_T = "KXMLBTOTAL-25JUN271905TEXLAA-9"
LEGS = [leg(SPREAD_T, "yes"),
        leg(TOTAL_T, "yes")]
# Second quote on the SAME game, resting on the total only.
SIBLING_LEGS = [leg(TOTAL_T, "yes")]

BASELINE = {SPREAD_T: {"yes_bid": 0.45, "yes_ask": 0.47},
            TOTAL_T: {"yes_bid": 0.52, "yes_ask": 0.54}}
SIBLING_BASELINE = {TOTAL_T: {"yes_bid": 0.52, "yes_ask": 0.54}}
# Spread jumped ~10c; total essentially unchanged.
JUMPED = {SPREAD_T: {"yes_bid": 0.55, "yes_ask": 0.57},
          TOTAL_T: {"yes_bid": 0.52, "yes_ask": 0.54}}
QUIET = {SPREAD_T: {"yes_bid": 0.455, "yes_ask": 0.475},
         TOTAL_T: {"yes_bid": 0.52, "yes_ask": 0.54}}

BOOK_FAIR_AT_QUOTE = 0.25


def _grid_df():
    """A healthy book cache that prices the combo at exactly 0.25 per cell."""
    rows = []
    for book in ("draftkings", "fanduel"):
        for cell in SPREAD_TOTAL_FAMILY:
            rows.append(dict(game_id="game1", combo=cell, period="FG",
                             bookmaker=book, sgp_decimal=3.8, fetch_time=None,
                             spread_line=-1.5, total_line=8.5))
    return pd.DataFrame(rows)


def _setup(monkeypatch, tmp_path, db_name, *, current, quotes=None,
           book_fair_now=None):
    """Seed open quote(s) with a placement baseline and stub the live poll.

    `current` is what the constituent poll returns this sweep; `book_fair_now`
    overrides the recomputed book consensus (None = leave the real cache fair,
    which equals BOOK_FAIR_AT_QUOTE, i.e. quiet books).
    """
    import kalshi_mlb_mm.config as cfg
    import kalshi_mlb_mm.db as db
    from kalshi_mlb_mm import main, singles
    monkeypatch.setattr(cfg, "DB_PATH", tmp_path / db_name)
    monkeypatch.setattr(cfg, "KILL_FILE", tmp_path / ".kill")
    importlib.reload(db)
    db.init_database()
    monkeypatch.setattr(main, "_SGP_ODDS", _grid_df())
    monkeypatch.setattr(main, "_resolve_game_for_legs", lambda gl: "game1")
    monkeypatch.setattr(main, "_ENGINE", None)
    monkeypatch.setattr(main, "_commence_time",
                        lambda gid: datetime(2099, 1, 1, tzinfo=timezone.utc))
    monkeypatch.setattr(cfg, "CONSTITUENT_JUMP_THRESHOLD", 0.03)
    # current=None leaves the real fetch in place (for the API-bound test).
    if current is not None:
        monkeypatch.setattr(singles, "fetch_market_prices",
                            lambda tickers, budget_sec=None: current)
    if book_fair_now is not None:
        monkeypatch.setattr(main, "_current_consensus_fair",
                            lambda legs_json: book_fair_now)
    now = datetime.now(timezone.utc)
    quotes = quotes or [("q1", LEGS, BASELINE)]
    with db.connect() as con:
        for qid, legs, baseline in quotes:
            con.execute(
                "INSERT INTO live_quotes (quote_id, rfq_id, combo_market_ticker, "
                "game_id, yes_bid, no_bid, model_fair, book_fair, blended_fair, "
                "status, submitted_at, closed_at, leg_prices_json) "
                "VALUES (?,?,?,?,?,?,?,?,?,?,?,?,?)",
                [qid, f"r-{qid}", f"COMBO-{qid}", "game1", 0.20, 0.72, None,
                 BOOK_FAIR_AT_QUOTE, BOOK_FAIR_AT_QUOTE, "open", now, None,
                 json.dumps(baseline) if baseline else None])
            con.execute(
                "INSERT OR REPLACE INTO seen_rfqs (rfq_id, market_ticker, in_scope, "
                "game_id, legs_json, first_seen_at, last_decision, creator_id) "
                "VALUES (?,?,?,?,?,?,?,?)",
                [f"r-{qid}", f"COMBO-{qid}", True, "game1", json.dumps(legs),
                 now, "quoted", ""])
    return main, db


def _state(db, qid="q1"):
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


# --------------------------------------------------------------------------- #
# The core case                                                               #
# --------------------------------------------------------------------------- #

def test_jump_with_quiet_books_cancels_the_quote(monkeypatch, tmp_path):
    """THE case: a constituent moved ~10c since we quoted while our book
    consensus barely moved -> we are the stale side -> pull."""
    main, db = _setup(monkeypatch, tmp_path, "cj_fire.duckdb", current=JUMPED)
    gw = _GW()
    main._risk_sweep_tick(gw)
    status, decision = _state(db)
    assert status == "cancelled"
    assert decision == ("sweep_cancel", "constituent_jump")
    assert gw.cancelled == ["q1"]


def test_quiet_constituents_leave_the_quote_open(monkeypatch, tmp_path):
    """A sub-threshold wiggle (~0.005) is not a jump."""
    main, db = _setup(monkeypatch, tmp_path, "cj_quiet.duckdb", current=QUIET)
    gw = _GW()
    main._risk_sweep_tick(gw)
    status, _ = _state(db)
    assert status == "open"
    assert gw.cancelled == []


# --------------------------------------------------------------------------- #
# Unconditional since #54                                                     #
# --------------------------------------------------------------------------- #

def test_jump_fires_regardless_of_book_movement(monkeypatch, tmp_path):
    """#54 live-only: the fair was live at placement, so ANY constituent jump
    during the resting window means the market moved after us — even when the
    books moved too (the pre-#54 book_quiet guard excused exactly this case
    for cache-priced quotes, and was deleted with cache pricing)."""
    main, db = _setup(monkeypatch, tmp_path, "cj_uncond.duckdb", current=JUMPED,
                      book_fair_now=0.30)          # 5c from the 0.25 at quote
    gw = _GW()
    main._risk_sweep_tick(gw)
    status, decision = _state(db)
    assert status == "cancelled"
    assert decision == ("sweep_cancel", "constituent_jump")


# --------------------------------------------------------------------------- #
# Fail-safe behaviour                                                         #
# --------------------------------------------------------------------------- #

def test_unreadable_constituents_do_not_cancel(monkeypatch, tmp_path):
    """A transient API failure must NOT flush the resting book — absence of
    signal is not a jump. #17's confirm veto still fails closed on any fill."""
    main, db = _setup(monkeypatch, tmp_path, "cj_unread.duckdb", current={})
    gw = _GW()
    main._risk_sweep_tick(gw)
    status, _ = _state(db)
    assert status == "open"
    assert gw.cancelled == []


def test_missing_placement_baseline_is_skipped_not_cancelled(monkeypatch, tmp_path):
    """A pre-#17 row has NULL leg_prices_json: no baseline, no jump verdict.
    The tipoff / staleness / book-drift gates still cover it."""
    main, db = _setup(monkeypatch, tmp_path, "cj_nobase.duckdb", current=JUMPED,
                      quotes=[("q1", LEGS, None)])
    gw = _GW()
    main._risk_sweep_tick(gw)
    status, _ = _state(db)
    assert status == "open"


# --------------------------------------------------------------------------- #
# Game-level propagation                                                      #
# --------------------------------------------------------------------------- #

def test_jump_cancels_sibling_quotes_touching_the_same_game(monkeypatch, tmp_path):
    """Quote 2 rests on a different market of the same game and never saw the
    jump itself — but its game is compromised, so it goes too."""
    main, db = _setup(monkeypatch, tmp_path, "cj_sibling.duckdb", current=JUMPED,
                      quotes=[("q1", LEGS, BASELINE),
                              ("q2", SIBLING_LEGS, SIBLING_BASELINE)])
    gw = _GW()
    main._risk_sweep_tick(gw)
    assert _state(db, "q1")[0] == "cancelled"
    assert _state(db, "q2")[0] == "cancelled"
    assert _state(db, "q2")[1] == ("sweep_cancel", "constituent_jump")
    assert set(gw.cancelled) == {"q1", "q2"}


# --------------------------------------------------------------------------- #
# Documented rate load                                                        #
# --------------------------------------------------------------------------- #

def test_api_calls_are_one_per_distinct_ticker_per_sweep(monkeypatch, tmp_path):
    """The documented bound, asserted rather than eyeballed: two quotes sharing
    the total leg cost ONE GET for it, not two."""
    from kalshi_common import auth_client
    # current=None -> the REAL fetch_market_prices runs, so the GETs are real.
    main, db = _setup(monkeypatch, tmp_path, "cj_bound.duckdb", current=None,
                      quotes=[("q1", LEGS, BASELINE),
                              ("q2", SIBLING_LEGS, SIBLING_BASELINE)])
    calls = []

    def fake_api(method, path, *a, **k):
        calls.append(path)
        return 200, {"market": {"yes_bid": 45, "yes_ask": 47}}, None

    monkeypatch.setattr(auth_client, "api", fake_api)
    main._risk_sweep_tick(_GW())
    market_calls = [c for c in calls if c.startswith("/markets/")]
    assert len(market_calls) == 2, f"expected 2 distinct GETs, got {market_calls}"
    assert set(market_calls) == {f"/markets/{SPREAD_T}", f"/markets/{TOTAL_T}"}


def test_poll_order_rotates_so_the_tail_is_not_starved(monkeypatch):
    """With a budget, a large book cannot be fully polled in one sweep. The
    start point must advance, or the same trailing tickers never get read."""
    from kalshi_mlb_mm import main
    monkeypatch.setattr(main, "_CONSTITUENT_POLL_CURSOR", 0)
    tickers = {"A", "B", "C", "D"}
    assert main._rotated_poll_order(tickers) == ["A", "B", "C", "D"]
    monkeypatch.setattr(main, "_CONSTITUENT_POLL_CURSOR", 2)
    assert main._rotated_poll_order(tickers) == ["C", "D", "A", "B"]
    monkeypatch.setattr(main, "_CONSTITUENT_POLL_CURSOR", 6)   # wraps
    assert main._rotated_poll_order(tickers) == ["C", "D", "A", "B"]
    assert main._rotated_poll_order(set()) == []


def test_sweep_advances_the_poll_cursor(monkeypatch, tmp_path):
    """Each sweep starts where the last one left off."""
    main, db = _setup(monkeypatch, tmp_path, "cj_cursor.duckdb", current=QUIET)
    monkeypatch.setattr(main, "_CONSTITUENT_POLL_CURSOR", 0)
    main._risk_sweep_tick(_GW())
    assert main._CONSTITUENT_POLL_CURSOR > 0


def test_sweep_makes_no_constituent_calls_when_flat(monkeypatch, tmp_path):
    """No resting quotes -> no constituent polling at all."""
    import kalshi_mlb_mm.config as cfg
    import kalshi_mlb_mm.db as db
    from kalshi_mlb_mm import main, singles
    monkeypatch.setattr(cfg, "DB_PATH", tmp_path / "cj_flat.duckdb")
    monkeypatch.setattr(cfg, "KILL_FILE", tmp_path / ".kill")
    importlib.reload(db)
    db.init_database()
    monkeypatch.setattr(main, "_SGP_ODDS", _grid_df())

    def must_not_call(tickers):
        raise AssertionError("flat book must not poll constituents")

    monkeypatch.setattr(singles, "fetch_market_prices", must_not_call)
    main._risk_sweep_tick(_GW())        # no exception == no polling


def test_unhealthy_books_skip_the_poll_and_cancel_anyway(monkeypatch, tmp_path):
    """#57: the "can't trust our fair" pull keys on #38 live-fetch health,
    not sweep-row age. When the alerter is dark everything is cancelled
    regardless, so spending API calls on constituents would be pure waste.
    (An EMPTY sweep frame alone must cancel nothing — test_sweep_demotion.)"""
    from kalshi_common.book_health import BookHealthAlerter
    from kalshi_mlb_mm import singles
    main, db = _setup(monkeypatch, tmp_path, "cj_stale.duckdb", current=JUMPED)

    def must_not_call(tickers, budget_sec=None):
        raise AssertionError("dark books must not poll constituents")

    monkeypatch.setattr(singles, "fetch_market_prices", must_not_call)
    monkeypatch.setattr(main, "_SGP_ODDS", pd.DataFrame())
    dark = BookHealthAlerter(label="T", streak_threshold=1,
                             min_healthy_books=2, paths=("on_demand",),
                             notifier=lambda t, m: None, notify_async=False)
    for book in ("draftkings", "fanduel"):
        dark.observe(book=book, path="on_demand", outcome="transport_error")
    gw = _GW()
    main._risk_sweep_tick(gw, book_health=dark)
    status, decision = _state(db)
    assert status == "cancelled"
    assert decision == ("sweep_cancel", "books_unhealthy")


# --------------------------------------------------------------------------- #
# Firehose                                                                    #
# --------------------------------------------------------------------------- #

def test_constituent_jump_event_carries_games_and_moves(monkeypatch, tmp_path):
    main, db = _setup(monkeypatch, tmp_path, "cj_evt.duckdb", current=JUMPED)
    events = []
    monkeypatch.setattr(main.research, "emit",
                        lambda ev, **kw: events.append((ev, kw)))
    main._risk_sweep_tick(_GW())
    jumps = [kw for ev, kw in events if ev == "constituent_jump"]
    assert len(jumps) == 1
    payload = jumps[0]["payload"]
    assert payload["games"] == ["game1"]
    assert payload["threshold"] == 0.03
    assert SPREAD_T in payload["moves"]["q1"]
    assert TOTAL_T not in payload["moves"]["q1"]
