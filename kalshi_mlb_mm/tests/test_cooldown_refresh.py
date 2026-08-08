"""P-7 (issue #21, re-pointed by #57): the post-fill cooldown must outlast
the data refresh, and same-price re-quotes must be blocked until consensus
fair actually moves.

The 60s COMBO_COOLDOWN_SEC floor expires while the counterparty who picked
us off can immediately post a new RFQ on the same combo. Two gates:

- `in_cooldown_awaiting_refresh`: floor passed, but no TARGETED on-demand
  fetch has completed post-fill for every game the combo touches (#57 — the
  old rule waited for a full sweep pass; sweep rows are now irrelevant, see
  test_sweep_demotion.py for the re-point regressions).
- `same_price_block`: books re-asked, but consensus fair is still within
  QUOTE_HYSTERESIS of the fair the fill transacted against.

Gate tests mirror the N7 harness (mocked scope cache / router fair, fresh
DuckDB per test); _SGP_ODDS stays None throughout — the cooldown lifecycle
must run entirely sweep-free. The unit block covers
`_post_fill_live_refresh_landed`'s cross-game and fail-closed edges.
"""
import importlib
from datetime import datetime, timedelta, timezone

from kalshi_common import legset

from kalshi_mlb_mm.tests.conftest import FakeLiveEngine, leg

_EVT = "KXMLBGAME-25JUN271905TEXLAA"
_LEGS = [{"market_ticker": "KXMLBSPREAD-25JUN271905TEXLAA-LAA2",
          "event_ticker": _EVT, "side": "yes"},
         {"market_ticker": "KXMLBTOTAL-25JUN271905TEXLAA-9",
          "event_ticker": _EVT, "side": "yes"}]

_TICKER = "COMBO-CD"
_GAME = "GAME-CD"
_FILL_FAIR = 0.55


def _cooldown_env(monkeypatch, tmp_path, db_name, *, completed_age,
                  cooled_until, fair_now):
    """Discovery-tick env with one prior fill on _TICKER (fair_at_confirm =
    _FILL_FAIR, filled 120s ago) and its combo_cooldown row. `completed_age`
    is the engine's post-fill-fetch landing age (None = no completed
    targeted fetch yet). Returns (db, main, src, gw, filled_at)."""
    import kalshi_mlb_mm.config as cfg
    import kalshi_mlb_mm.db as db
    import kalshi_mlb_mm.risk as risk
    import kalshi_mlb_mm.router as router_mod
    from kalshi_mlb_mm import main

    monkeypatch.setattr(cfg, "DB_PATH", tmp_path / db_name)
    monkeypatch.setattr(cfg, "KILL_FILE", tmp_path / ".kill")
    importlib.reload(db)
    db.init_database()

    # Sweep-free on purpose (#57): the cooldown lifecycle must never need it.
    monkeypatch.setattr(main, "_SGP_ODDS", None)
    monkeypatch.setattr(risk, "tipoff_ok", lambda ct, m: True)
    monkeypatch.setattr(router_mod, "combo_fair_detail",
                        lambda *a, **k: (router_mod.ComboFair(fair_now, 0.0, 1), "ok"))
    monkeypatch.setattr(main, "_commence_time", lambda gid: None)
    monkeypatch.setattr(main, "_PREV_BOOK_FAIR", {})
    monkeypatch.setattr(main, "_SCOPE_CACHE", {_TICKER: (True, _GAME, _LEGS)})
    monkeypatch.setattr(main, "_resolve_game_for_legs", lambda gl: _GAME)
    from mlb_sgp._shared import GameRef
    monkeypatch.setattr(main, "_game_ref",
                        lambda gid: GameRef(game_id=_GAME, home_team="H",
                                            away_team="A", commence_time=None))
    monkeypatch.setattr(main, "_leg_market_prices",
                        lambda legs: {"L": {"yes_bid": 0.5, "yes_ask": 0.52}})
    monkeypatch.setattr(main, "_ENGINE",
                        FakeLiveEngine(completed_age=completed_age))

    filled_at = datetime.now(timezone.utc) - timedelta(seconds=120)
    with db.connect() as con:
        con.execute(
            "INSERT INTO fills (fill_id, quote_id, rfq_id, combo_market_ticker, "
            "game_id, side_held, contracts, price, fee, model_fair_at_quote, "
            "book_fair_at_quote, blended_fair_at_quote, fair_at_confirm, "
            "realized_pnl, filled_at, reconciled) "
            "VALUES (?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?)",
            ["fill-cd", "qid-cd", "r-cd-prior", _TICKER, _GAME, "yes", 1, 0.50,
             0.0, None, _FILL_FAIR, _FILL_FAIR, _FILL_FAIR, None, filled_at, True])
        con.execute("INSERT OR REPLACE INTO fill_games (fill_id, game_id) "
                    "VALUES (?, ?)", ["fill-cd", _GAME])
        con.execute(
            "INSERT OR REPLACE INTO combo_cooldown "
            "(combo_market_ticker, cooled_until) VALUES (?, ?)",
            [_TICKER, cooled_until])

    class Src:
        def poll(self):
            return [{"id": "r-cd-new", "market_ticker": _TICKER, "contracts": 3}]

        def get_market(self, t):
            return {}

    class GW:
        def __init__(self):
            self.submits = []

        def submit_quote(self, *a):
            self.submits.append(a)
            return "qid-cd-new"

    return db, main, Src(), GW(), filled_at


def _last_decision(db, rfq_id="r-cd-new"):
    with db.connect(read_only=True) as con:
        return con.execute(
            "SELECT decision, reason FROM quote_decisions "
            "WHERE rfq_id=? ORDER BY observed_at DESC LIMIT 1",
            [rfq_id]).fetchone()


def test_floor_active_still_skips_in_cooldown(monkeypatch, tmp_path):
    """Pins existing behavior: within the COMBO_COOLDOWN_SEC floor the combo
    is skipped with the original in_cooldown reason."""
    now = datetime.now(timezone.utc)
    db, main, src, gw, _ = _cooldown_env(
        monkeypatch, tmp_path, "cd_floor.duckdb",
        completed_age=10.0,
        cooled_until=now + timedelta(seconds=60),   # floor still active
        fair_now=0.58)
    main._discovery_tick(src, gw, dry_run=False)
    assert gw.submits == []
    assert _last_decision(db) == ("skipped", "in_cooldown")


def test_floor_passed_but_no_targeted_fetch_skips(monkeypatch, tmp_path):
    """THE regression test for issue #21: cooldown floor expired, but no
    targeted post-fill fetch has completed → the books were never re-asked
    after the picked-off fill → must stay cooled (and re-trigger the fetch)."""
    now = datetime.now(timezone.utc)
    db, main, src, gw, _ = _cooldown_env(
        monkeypatch, tmp_path, "cd_stale.duckdb",
        completed_age=None,                          # no post-fill landing
        cooled_until=now - timedelta(seconds=60),    # floor expired
        fair_now=0.58)
    main._discovery_tick(src, gw, dry_run=False)
    assert gw.submits == [], (
        "must not re-quote before a targeted post-fill fetch lands")
    assert _last_decision(db) == ("skipped", "in_cooldown_awaiting_refresh")
    assert main._ENGINE.ensure_calls, (
        "the gate must lazily re-trigger the targeted fetch")


def test_post_fill_fetch_and_fair_moved_allows_quote(monkeypatch, tmp_path):
    """Targeted fetch completed after the fill AND consensus fair moved →
    quote."""
    now = datetime.now(timezone.utc)
    db, main, src, gw, _ = _cooldown_env(
        monkeypatch, tmp_path, "cd_moved.duckdb",
        completed_age=10.0,                          # landed 10s ago, post-fill
        cooled_until=now - timedelta(seconds=60),
        fair_now=0.58)                               # moved 3¢
    main._discovery_tick(src, gw, dry_run=False)
    assert len(gw.submits) == 1, "refreshed + moved fair must quote"
    assert _last_decision(db)[0] == "quoted"


def test_post_fill_fetch_but_fair_unchanged_same_price_block(monkeypatch, tmp_path):
    """Books re-asked after the fill but consensus fair is within
    QUOTE_HYSTERESIS of the fill's fair → same_price_block."""
    now = datetime.now(timezone.utc)
    db, main, src, gw, _ = _cooldown_env(
        monkeypatch, tmp_path, "cd_same.duckdb",
        completed_age=10.0,
        cooled_until=now - timedelta(seconds=60),
        fair_now=_FILL_FAIR + 0.002)                # within ε=0.005
    main._discovery_tick(src, gw, dry_run=False)
    assert gw.submits == [], (
        "unchanged fair must not be re-quoted at the just-picked-off price")
    assert _last_decision(db) == ("skipped", "same_price_block")


def test_no_cooldown_row_means_no_new_gates(monkeypatch, tmp_path):
    """A combo that never filled (no combo_cooldown row) is untouched by the
    refresh/same-price gates even with no targeted fetch and unchanged fair."""
    now = datetime.now(timezone.utc)
    db, main, src, gw, _ = _cooldown_env(
        monkeypatch, tmp_path, "cd_none.duckdb",
        completed_age=None,
        cooled_until=now - timedelta(seconds=60),
        fair_now=_FILL_FAIR)
    # Remove the cooldown row AND the fill so this combo looks never-filled.
    with db.connect() as con:
        con.execute("DELETE FROM combo_cooldown WHERE combo_market_ticker=?",
                    [_TICKER])
        con.execute("DELETE FROM fills WHERE combo_market_ticker=?", [_TICKER])
        con.execute("DELETE FROM fill_games WHERE fill_id='fill-cd'")
    main._discovery_tick(src, gw, dry_run=False)
    assert len(gw.submits) == 1
    assert _last_decision(db)[0] == "quoted"


# ---------------------------------------------------------------------------
# _post_fill_live_refresh_landed unit coverage (cross-game + fail-closed +
# timezone edges), against per-hash landing ages.
# ---------------------------------------------------------------------------
_G1_LEGS = [leg("KXMLBSPREAD-25JUN271905TEXLAA-LAA2", "yes"),
            leg("KXMLBTOTAL-25JUN271905TEXLAA-9", "yes")]
_G2_LEGS = [leg("KXMLBTOTAL-25JUN271910NYYBOS-9", "yes")]


class _GateEng:
    """Engine stub with per-hash completed-landing ages."""

    def __init__(self, ages):
        self.ages = dict(ages)

    def completed_fetch_age_sec(self, h):
        return self.ages.get(h)


def _by_game(*leg_lists):
    canon = legset.parse_legs([l for ls in leg_lists for l in ls])
    return legset.partition_by_game(canon)


def _hash_of(leg_list):
    return legset.leg_set_hash(legset.parse_legs(leg_list))


def _check(monkeypatch, by_game, filled_at, ages):
    from kalshi_mlb_mm import main
    monkeypatch.setattr(main, "_ENGINE", _GateEng(ages))
    return main._post_fill_live_refresh_landed(by_game, filled_at)


def test_gate_post_fill_landing_true(monkeypatch):
    fill = datetime.now(timezone.utc) - timedelta(seconds=120)
    bg = _by_game(_G1_LEGS)
    assert _check(monkeypatch, bg, fill, {_hash_of(_G1_LEGS): 10.0}) is True


def test_gate_pre_fill_landing_false(monkeypatch):
    """A flight that completed BEFORE the fill is the picked-off data."""
    fill = datetime.now(timezone.utc) - timedelta(seconds=120)
    bg = _by_game(_G1_LEGS)
    assert _check(monkeypatch, bg, fill, {_hash_of(_G1_LEGS): 300.0}) is False


def test_gate_cross_game_requires_every_game(monkeypatch):
    fill = datetime.now(timezone.utc) - timedelta(seconds=120)
    bg = _by_game(_G1_LEGS, _G2_LEGS)
    assert len(bg) == 2, "harness must produce two games"
    one_landed = {_hash_of(_G1_LEGS): 10.0}          # g2 missing
    assert _check(monkeypatch, bg, fill, one_landed) is False
    both = {_hash_of(_G1_LEGS): 10.0, _hash_of(_G2_LEGS): 5.0}
    assert _check(monkeypatch, bg, fill, both) is True


def test_gate_naive_filled_at_is_session_local(monkeypatch):
    """fills.filled_at is TIMESTAMPTZ (aware), but a naive value must be
    treated as session-LOCAL (risk._now_matching convention), not UTC."""
    naive_local_fill = datetime.now() - timedelta(seconds=120)
    bg = _by_game(_G1_LEGS)
    assert _check(monkeypatch, bg, naive_local_fill,
                  {_hash_of(_G1_LEGS): 10.0}) is True
    assert _check(monkeypatch, bg, naive_local_fill,
                  {_hash_of(_G1_LEGS): 300.0}) is False


def test_gate_fails_closed(monkeypatch):
    from kalshi_mlb_mm import main
    fill = datetime.now(timezone.utc)
    bg = _by_game(_G1_LEGS)
    # No engine at all.
    monkeypatch.setattr(main, "_ENGINE", None)
    assert main._post_fill_live_refresh_landed(bg, fill) is False
    # No landing for the hash / no filled_at.
    assert _check(monkeypatch, bg, fill, {}) is False
    assert _check(monkeypatch, bg, None, {_hash_of(_G1_LEGS): 1.0}) is False
