"""R-1 (issue #22): per-game and daily cap gates must count OPEN quotes'
worst-case exposure, not just fills. Up to MAX_OPEN_QUOTES resting quotes on
different combos of one game could previously fill in a single burst to many
multiples of the per-game / daily caps because those gates read `fills` only.

Gate tests mirror the N7 harness (test_n7_n8_n9_n10_n11_n12.py): mocked
_SGP_ODDS / scope cache / router fair, fresh DuckDB per test, a Source that
serves one RFQ, and a GW that records submits.
"""
import importlib
from datetime import datetime, timezone

import pandas as pd

_EVT = "KXMLBGAME-25JUN271905TEXLAA"
_LEGS = [{"market_ticker": "KXMLBSPREAD-25JUN271905TEXLAA-LAA2",
          "event_ticker": _EVT, "side": "yes"},
         {"market_ticker": "KXMLBTOTAL-25JUN271905TEXLAA-9",
          "event_ticker": _EVT, "side": "yes"}]


def _setup(monkeypatch, tmp_path, db_name, candidate_ticker, candidate_game):
    """Common gate-test harness. Returns the reloaded db module."""
    import kalshi_mlb_mm.config as cfg
    import kalshi_mlb_mm.db as db
    import kalshi_mlb_mm.risk as risk
    import kalshi_mlb_mm.router as router_mod
    from kalshi_mlb_mm import main

    monkeypatch.setattr(cfg, "DB_PATH", tmp_path / db_name)
    monkeypatch.setattr(cfg, "KILL_FILE", tmp_path / ".kill")
    monkeypatch.setattr(cfg, "BANKROLL", 500.0)
    monkeypatch.setattr(cfg, "MAX_FILL_EXPOSURE_PCT", 0.10)   # per-fill cap $50
    monkeypatch.setattr(cfg, "MAX_GAME_EXPOSURE_PCT", 0.10)   # per-game cap $50
    monkeypatch.setattr(cfg, "DAILY_EXPOSURE_CAP_PCT", 0.75)  # daily cap $375
    importlib.reload(db)
    db.init_database()

    monkeypatch.setattr(main, "_SGP_ODDS",
                        pd.DataFrame({"game_id": [candidate_game], "combo": ["c"],
                                      "period": ["FG"], "bookmaker": ["dk"],
                                      "sgp_decimal": [2.0], "fetch_time": [None],
                                      "spread_line": [-1.5], "total_line": [8.5]}))
    monkeypatch.setattr(risk, "tipoff_ok", lambda ct, m: True)
    monkeypatch.setattr(router_mod, "combo_fair", lambda *a, **k: 0.55)
    monkeypatch.setattr(main, "_commence_time", lambda gid: None)
    monkeypatch.setattr(main, "_PREV_BOOK_FAIR", {})
    monkeypatch.setattr(main, "_SCOPE_CACHE",
                        {candidate_ticker: (True, candidate_game, _LEGS)})
    monkeypatch.setattr(main, "_resolve_game_for_legs", lambda gl: candidate_game)
    monkeypatch.setattr(main, "_leg_market_prices",
                        lambda legs: {"L": {"yes_bid": 0.5, "yes_ask": 0.52}})
    return db


def _insert_open_quote(con, quote_id, ticker, primary_game, worst_exposure_usd,
                       games=None):
    """One OPEN live_quotes row + its quote_games attribution rows."""
    con.execute(
        "INSERT INTO live_quotes (quote_id, rfq_id, combo_market_ticker, "
        "game_id, yes_bid, no_bid, model_fair, book_fair, blended_fair, "
        "status, submitted_at, closed_at, worst_exposure_usd) "
        "VALUES (?,?,?,?,?,?,?,?,?,?,?,?,?)",
        [quote_id, "r_" + quote_id, ticker, primary_game, 0.50, 0.43, None,
         0.55, 0.55, "open", datetime.now(timezone.utc), None,
         worst_exposure_usd])
    for g in (games if games is not None else [primary_game]):
        con.execute("INSERT OR REPLACE INTO quote_games (quote_id, game_id) "
                    "VALUES (?, ?)", [quote_id, g])


class _GW:
    def __init__(self):
        self.submits = []

    def submit_quote(self, *a):
        self.submits.append(a)
        return "qid-new"


class _Src:
    def __init__(self, rfq):
        self._rfq = rfq

    def poll(self):
        return [self._rfq]

    def get_market(self, t):
        return {}


def _last_decision(db, rfq_id):
    with db.connect(read_only=True) as con:
        return con.execute(
            "SELECT decision, reason FROM quote_decisions "
            "WHERE rfq_id=? ORDER BY observed_at DESC LIMIT 1",
            [rfq_id]).fetchone()


# ---------------------------------------------------------------------------
# Acceptance 1: open quotes summing to the PER-GAME cap block the next quote
# on that game — with the fills table completely empty.
# ---------------------------------------------------------------------------
def test_open_quotes_reaching_per_game_cap_block_next_quote(monkeypatch, tmp_path):
    from kalshi_mlb_mm import main
    db = _setup(monkeypatch, tmp_path, "oq_pg.duckdb", "COMBO-PG-NEW", "GAME-PG")
    # Two open quotes on other combos of the SAME game: $30 + $20 = $50 = cap.
    with db.connect() as con:
        _insert_open_quote(con, "oq-pg-1", "COMBO-PG-1", "GAME-PG", 30.0)
        _insert_open_quote(con, "oq-pg-2", "COMBO-PG-2", "GAME-PG", 20.0)

    gw = _GW()
    main._discovery_tick(_Src({"id": "r-pg-new", "market_ticker": "COMBO-PG-NEW",
                               "contracts": 3}), gw, dry_run=False)

    assert gw.submits == [], (
        "per-game cap must count open quotes' worst-case exposure")
    d = _last_decision(db, "r-pg-new")
    assert d == ("skipped", "per_game_cap"), f"expected per_game_cap, got {d}"


# ---------------------------------------------------------------------------
# Acceptance 2: open quotes across DIFFERENT games summing to the DAILY cap
# block a quote on a fresh game (per-game headroom for that game is clean).
# ---------------------------------------------------------------------------
def test_open_quotes_reaching_daily_cap_block_next_quote(monkeypatch, tmp_path):
    from kalshi_mlb_mm import main
    db = _setup(monkeypatch, tmp_path, "oq_daily.duckdb", "COMBO-D-NEW", "GAME-D-NEW")
    # 8 open quotes at $50 on 8 OTHER games = $400 >= $375 daily cap.
    with db.connect() as con:
        for i in range(8):
            _insert_open_quote(con, f"oq-d-{i}", f"COMBO-D-{i}", f"GAME-D-{i}", 50.0)

    gw = _GW()
    main._discovery_tick(_Src({"id": "r-d-new", "market_ticker": "COMBO-D-NEW",
                               "contracts": 3}), gw, dry_run=False)

    assert gw.submits == [], (
        "daily cap must count open quotes' worst-case exposure")
    d = _last_decision(db, "r-d-new")
    assert d == ("skipped", "daily_cap"), f"expected daily_cap, got {d}"


# ---------------------------------------------------------------------------
# Acceptance 3: a CROSS-GAME open quote counts against EVERY game it touches —
# a new quote on the secondary game must see no headroom.
# ---------------------------------------------------------------------------
def test_cross_game_open_quote_counts_against_secondary_game(monkeypatch, tmp_path):
    from kalshi_mlb_mm import main
    db = _setup(monkeypatch, tmp_path, "oq_xg.duckdb", "COMBO-XG-NEW", "GAME-XG-B")
    # One open cross-game quote (primary GAME-XG-A) touching A and B at $50.
    with db.connect() as con:
        _insert_open_quote(con, "oq-xg", "COMBO-XG", "GAME-XG-A", 50.0,
                           games=["GAME-XG-A", "GAME-XG-B"])

    gw = _GW()
    main._discovery_tick(_Src({"id": "r-xg-new", "market_ticker": "COMBO-XG-NEW",
                               "contracts": 3}), gw, dry_run=False)

    assert gw.submits == [], (
        "secondary game of a cross-game open quote must be blocked at the cap")
    d = _last_decision(db, "r-xg-new")
    assert d == ("skipped", "per_game_cap"), f"expected per_game_cap, got {d}"


# ---------------------------------------------------------------------------
# Headroom sanity: open exposure BELOW the caps must still quote (the gates
# must not over-block).
# ---------------------------------------------------------------------------
def test_open_quotes_below_caps_still_quote(monkeypatch, tmp_path):
    from kalshi_mlb_mm import main
    db = _setup(monkeypatch, tmp_path, "oq_ok.duckdb", "COMBO-OK-NEW", "GAME-OK")
    with db.connect() as con:
        _insert_open_quote(con, "oq-ok-1", "COMBO-OK-1", "GAME-OK", 20.0)

    gw = _GW()
    main._discovery_tick(_Src({"id": "r-ok-new", "market_ticker": "COMBO-OK-NEW",
                               "contracts": 3}), gw, dry_run=False)

    assert len(gw.submits) == 1, "under-cap open exposure must not block quoting"
    d = _last_decision(db, "r-ok-new")
    assert d[0] == "quoted", f"expected quoted, got {d}"


# ---------------------------------------------------------------------------
# Conservative fallbacks in the exposure readers (helper-level).
# ---------------------------------------------------------------------------
def test_own_open_quote_excluded_when_requoting_same_rfq(monkeypatch, tmp_path):
    """Replace safety: the RFQ's OWN open quote is superseded on re-quote
    (hysteresis keeps it, or Kalshi auto-cancels it when the replacement
    lands), so it must NOT count against the caps — otherwise every replace
    self-blocks once its game nears the cap."""
    import kalshi_mlb_mm.config as cfg
    import kalshi_mlb_mm.pricing as pricing_mod
    from kalshi_mlb_mm.pricing import Quote
    from kalshi_mlb_mm import main
    db = _setup(monkeypatch, tmp_path, "oq_self.duckdb", "COMBO-SELF", "GAME-SELF")
    # The PRE-EXISTING per-combo gate (N7) also counts the own open quote in
    # worst_inflight and would fire first at the default $50 cap; issue #22
    # keeps N7 as-is, so give it headroom (same as test_quote_replace) to pin
    # the per-game/daily exclusion specifically.
    monkeypatch.setattr(cfg, "MAX_COMBO_EXPOSURE_USD", 200.0)
    # Existing open quote for THIS rfq at the full per-game cap ($50). Its
    # prices differ from the new quote by >> hysteresis → replace path.
    with db.connect() as con:
        con.execute(
            "INSERT INTO live_quotes (quote_id, rfq_id, combo_market_ticker, "
            "game_id, yes_bid, no_bid, model_fair, book_fair, blended_fair, "
            "status, submitted_at, closed_at, worst_exposure_usd) "
            "VALUES (?,?,?,?,?,?,?,?,?,?,?,?,?)",
            ["oq-self-old", "r-self", "COMBO-SELF", "GAME-SELF", 0.500, 0.430,
             None, 0.55, 0.55, "open", datetime.now(timezone.utc), None, 50.0])
        con.execute("INSERT OR REPLACE INTO quote_games (quote_id, game_id) "
                    "VALUES (?, ?)", ["oq-self-old", "GAME-SELF"])
    monkeypatch.setattr(pricing_mod, "quote",
                        lambda fair, roi: Quote(yes_bid=0.520, no_bid=0.410))

    gw = _GW()
    main._discovery_tick(_Src({"id": "r-self", "market_ticker": "COMBO-SELF",
                               "contracts": 3}), gw, dry_run=False)

    assert len(gw.submits) == 1, "re-quote of own RFQ must not self-block on caps"
    with db.connect(read_only=True) as con:
        rows = dict(con.execute(
            "SELECT quote_id, status FROM live_quotes").fetchall())
    assert rows == {"oq-self-old": "replaced", "qid-new": "open"}, rows


def test_null_worst_exposure_counts_at_per_fill_cap(monkeypatch, tmp_path):
    """A pre-migration open row (worst_exposure_usd NULL) counts at
    max_fill_exposure_usd() — N8 convention — never as zero."""
    import kalshi_mlb_mm.config as cfg
    import kalshi_mlb_mm.db as db
    from kalshi_mlb_mm import main
    monkeypatch.setattr(cfg, "DB_PATH", tmp_path / "oq_null.duckdb")
    monkeypatch.setattr(cfg, "BANKROLL", 500.0)
    monkeypatch.setattr(cfg, "MAX_FILL_EXPOSURE_PCT", 0.10)
    importlib.reload(db)
    db.init_database()
    with db.connect() as con:
        _insert_open_quote(con, "oq-null", "COMBO-NULL", "GAME-NULL", None)

    daily_rows = main._open_quote_exposure_rows()
    assert [r["price"] for r in daily_rows] == [50.0]
    game_rows = main._open_quote_exposure_by_game()
    assert game_rows == [{"game_id": "GAME-NULL", "price": 50.0}]


def test_open_quote_without_quote_games_falls_back_to_primary(monkeypatch, tmp_path):
    """Defense-in-depth: an open row with NO quote_games mapping still counts
    against its primary game_id (LEFT JOIN + COALESCE), never disappears."""
    import kalshi_mlb_mm.config as cfg
    import kalshi_mlb_mm.db as db
    from kalshi_mlb_mm import main
    monkeypatch.setattr(cfg, "DB_PATH", tmp_path / "oq_orphan.duckdb")
    importlib.reload(db)
    db.init_database()
    with db.connect() as con:
        con.execute(
            "INSERT INTO live_quotes (quote_id, rfq_id, combo_market_ticker, "
            "game_id, yes_bid, no_bid, model_fair, book_fair, blended_fair, "
            "status, submitted_at, closed_at, worst_exposure_usd) "
            "VALUES (?,?,?,?,?,?,?,?,?,?,?,?,?)",
            ["oq-orphan", "r-orphan", "COMBO-ORPHAN", "GAME-ORPHAN", 0.50, 0.43,
             None, 0.55, 0.55, "open", datetime.now(timezone.utc), None, 25.0])

    game_rows = main._open_quote_exposure_by_game()
    assert game_rows == [{"game_id": "GAME-ORPHAN", "price": 25.0}]


def test_closed_quotes_do_not_count(monkeypatch, tmp_path):
    """Only status='open' rows contribute — filled/cancelled/replaced quotes'
    exposure is owned by the fills path (or gone)."""
    import kalshi_mlb_mm.config as cfg
    import kalshi_mlb_mm.db as db
    from kalshi_mlb_mm import main
    monkeypatch.setattr(cfg, "DB_PATH", tmp_path / "oq_closed.duckdb")
    importlib.reload(db)
    db.init_database()
    now = datetime.now(timezone.utc)
    with db.connect() as con:
        for status in ("filled", "cancelled", "replaced", "voided"):
            con.execute(
                "INSERT INTO live_quotes (quote_id, rfq_id, combo_market_ticker, "
                "game_id, yes_bid, no_bid, model_fair, book_fair, blended_fair, "
                "status, submitted_at, closed_at, worst_exposure_usd) "
                "VALUES (?,?,?,?,?,?,?,?,?,?,?,?,?)",
                [f"oq-{status}", f"r-{status}", "COMBO-C", "GAME-C", 0.50, 0.43,
                 None, 0.55, 0.55, status, now, now, 50.0])
            con.execute("INSERT OR REPLACE INTO quote_games (quote_id, game_id) "
                        "VALUES (?, ?)", [f"oq-{status}", "GAME-C"])

    assert main._open_quote_exposure_rows() == []
    assert main._open_quote_exposure_by_game() == []
