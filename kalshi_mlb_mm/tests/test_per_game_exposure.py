"""Per-game exposure ledger (fill_games) — cross-game combos must count their
full stake against EVERY game they touch, without double-counting the daily cap
or P&L. See main._today_fills_by_game / _fill_game_ids and the fill_games table.
"""
import importlib
from datetime import datetime, timezone


def _fresh_db(tmp_path, monkeypatch, name="pge.duckdb"):
    import kalshi_mlb_mm.config as cfg
    monkeypatch.setattr(cfg, "DB_PATH", tmp_path / name)
    import kalshi_mlb_mm.db as db
    importlib.reload(db)
    db.init_database()
    return db


def _insert_fill(con, fill_id, game_id, price, contracts, games):
    """One fills row (one-per-combo) + one fill_games row per game it touches."""
    con.execute(
        "INSERT INTO fills (fill_id, quote_id, rfq_id, combo_market_ticker, "
        "game_id, side_held, contracts, price, fee, filled_at, reconciled) "
        "VALUES (?,?,?,?,?,?,?,?,?,?,?)",
        [fill_id, "q_" + fill_id, "r_" + fill_id, "TICK-" + fill_id, game_id,
         "yes", contracts, price, 0.0, datetime.now(timezone.utc), True])
    for g in games:
        con.execute("INSERT OR REPLACE INTO fill_games (fill_id, game_id) VALUES (?, ?)",
                    [fill_id, g])


def test_fill_games_table_exists(tmp_path, monkeypatch):
    db = _fresh_db(tmp_path, monkeypatch)
    with db.connect(read_only=True) as con:
        # querying the table must not raise (empty is fine)
        n = con.execute("SELECT COUNT(*) FROM fill_games").fetchone()[0]
    assert n == 0


def test_cross_game_exposure_counts_against_both_games(tmp_path, monkeypatch):
    db = _fresh_db(tmp_path, monkeypatch)
    from kalshi_mlb_mm import main
    with db.connect() as con:
        # Cross-game combo: exposure 0.30*10 = 3.0, touches GAME_A and GAME_B.
        _insert_fill(con, "X", "GAME_A", 0.30, 10, ["GAME_A", "GAME_B"])
        # Single-game combo on GAME_A: exposure 0.20*5 = 1.0.
        _insert_fill(con, "Y", "GAME_A", 0.20, 5, ["GAME_A"])

    rows = main._today_fills_by_game()
    a = sum(r["price"] for r in rows if r["game_id"] == "GAME_A")
    b = sum(r["price"] for r in rows if r["game_id"] == "GAME_B")
    # GAME_A sees both combos (3.0 + 1.0); GAME_B sees the cross-game combo (3.0).
    assert a == 4.0
    assert b == 3.0  # <-- the fix: B's exposure as a SECONDARY leg is now visible


def test_daily_cap_not_double_counted(tmp_path, monkeypatch):
    """`_today_fills` (daily cap source) stays one-row-per-combo, so the cross-game
    combo's stake is counted ONCE for the daily cap even though it fans out to two
    games for the per-game cap."""
    db = _fresh_db(tmp_path, monkeypatch)
    from kalshi_mlb_mm import main
    with db.connect() as con:
        _insert_fill(con, "X", "GAME_A", 0.30, 10, ["GAME_A", "GAME_B"])  # 3.0
        _insert_fill(con, "Y", "GAME_A", 0.20, 5, ["GAME_A"])             # 1.0
    daily = sum(r["price"] for r in main._today_fills())
    assert daily == 4.0  # 3.0 + 1.0, NOT 7.0


def test_per_game_cap_blocks_secondary_game(tmp_path, monkeypatch):
    """After a cross-game (A,B) fill, a game that was only ever the 2nd leg (B)
    must be blocked once its attributed exposure crosses the cap."""
    import kalshi_mlb_mm.risk as risk
    db = _fresh_db(tmp_path, monkeypatch)
    from kalshi_mlb_mm import main
    with db.connect() as con:
        _insert_fill(con, "X", "GAME_A", 0.50, 10, ["GAME_A", "GAME_B"])  # 5.0 each

    rows = main._today_fills_by_game()
    # bankroll 100 * pct 0.02 = $2 cap. B's attributed exposure is $5 → blocked.
    assert risk.per_game_cap_ok("GAME_B", rows, bankroll=100.0, pct=0.02) is False
    # Same combo, larger cap ($10) → allowed.
    assert risk.per_game_cap_ok("GAME_B", rows, bankroll=100.0, pct=0.10) is True


def test_fill_game_ids_failsafe_to_primary(tmp_path, monkeypatch):
    """_fill_game_ids ALWAYS includes the primary game_id, even when legs can't
    be parsed/resolved — so a recorded fill is never lost from the cap."""
    _fresh_db(tmp_path, monkeypatch)
    from kalshi_mlb_mm import main
    # Empty legs → degrade to primary only.
    assert main._fill_game_ids([], "PRIMARY") == ["PRIMARY"]
    assert main._fill_game_ids(None, "PRIMARY") == ["PRIMARY"]


def test_fill_game_ids_collects_and_dedupes_games(tmp_path, monkeypatch):
    """_fill_game_ids returns every resolved game plus the primary, deduped/sorted.
    legset is mocked so this pins the collect/dedupe logic, not ticker parsing."""
    _fresh_db(tmp_path, monkeypatch)
    from kalshi_mlb_mm import main
    # Fake a 2-game canon; resolver maps each game-leg-group to an id.
    monkeypatch.setattr(main.legset, "parse_legs", lambda legs: ["canon"])
    monkeypatch.setattr(main.legset, "partition_by_game",
                        lambda canon: {"A": ["legA"], "B": ["legB"]})
    monkeypatch.setattr(main, "_resolve_game_for_legs",
                        lambda gl: "GAME_A" if gl == ["legA"] else "GAME_B")
    # Primary GAME_A overlaps a resolved id → must dedupe to two entries, sorted.
    assert main._fill_game_ids([{"x": 1}], "GAME_A") == ["GAME_A", "GAME_B"]


def test_fill_game_ids_skips_unresolvable_game_but_keeps_primary(tmp_path, monkeypatch):
    """If one game resolves to None, it's dropped, but the primary is always kept."""
    _fresh_db(tmp_path, monkeypatch)
    from kalshi_mlb_mm import main
    monkeypatch.setattr(main.legset, "parse_legs", lambda legs: ["canon"])
    monkeypatch.setattr(main.legset, "partition_by_game",
                        lambda canon: {"A": ["legA"], "B": ["legB"]})
    monkeypatch.setattr(main, "_resolve_game_for_legs",
                        lambda gl: "GAME_A" if gl == ["legA"] else None)
    assert main._fill_game_ids([{"x": 1}], "GAME_A") == ["GAME_A"]


def test_fill_game_ids_never_raises_degrades_to_primary(tmp_path, monkeypatch):
    """Even if leg parsing throws unexpectedly, _fill_game_ids must NOT raise —
    it degrades to primary-only so the fill is never lost from the per-game cap
    and the exception can't escape into the fill-write path."""
    _fresh_db(tmp_path, monkeypatch)
    from kalshi_mlb_mm import main

    def _boom(_legs):
        raise ValueError("bad leg shape")

    monkeypatch.setattr(main.legset, "parse_legs", _boom)
    assert main._fill_game_ids([{"x": 1}], "PRIMARY") == ["PRIMARY"]


def test_orphan_fill_invisible_to_per_game_but_caught_by_daily(tmp_path, monkeypatch):
    """Defense-in-depth: a fills row with NO fill_games rows (should never happen
    given the fail-safe) is invisible to the per-game cap (INNER JOIN) but is
    STILL counted by the daily cap (which reads `fills` directly) — so a combo is
    never fully lost from exposure control."""
    db = _fresh_db(tmp_path, monkeypatch)
    from kalshi_mlb_mm import main
    with db.connect() as con:
        # Deliberately insert a fills row WITHOUT any fill_games mapping.
        con.execute(
            "INSERT INTO fills (fill_id, quote_id, rfq_id, combo_market_ticker, "
            "game_id, side_held, contracts, price, fee, filled_at, reconciled) "
            "VALUES (?,?,?,?,?,?,?,?,?,?,?)",
            ["orphan", "q", "r", "TICK", "GAME_A", "yes", 10, 0.40, 0.0,
             datetime.now(timezone.utc), True])
    per_game = main._today_fills_by_game()
    daily = main._today_fills()
    assert per_game == []                      # INNER JOIN drops the unmapped fill
    assert sum(r["price"] for r in daily) == 4.0  # daily cap still sees it (0.40*10)
