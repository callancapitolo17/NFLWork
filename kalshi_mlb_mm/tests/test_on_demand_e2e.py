"""Phase 2 end-to-end loop invariants: REAL OnDemandEngine + fake SGPService
driven through main._discovery_tick — pending → landed → quoted → feed
re-quote on movement → feed stops when the RFQ leaves the poll."""
import importlib
import threading


from kalshi_common import legset
from kalshi_mlb_mm.on_demand import OnDemandEngine, QUOTE_FRESH_SEC
from mlb_sgp._shared import GameRef, OnDemandBookResult
from kalshi_mlb_mm.tests.conftest import leg


EVT = "25JUN271905TEXLAA"
OD_LEGS = [
    leg("KXMLBSPREAD-25JUN271905TEXLAA-LAA2", "yes"),
    leg("KXMLBTOTAL-25JUN271905TEXLAA-9", "yes"),
    leg("KXMLBGAME-25JUN271905TEXLAA-LAA", "yes"),
]
GREF = GameRef(game_id="game1", home_team="Los Angeles Angels",
               away_team="Texas Rangers", commence_time=None)


class FakeClock:
    def __init__(self):
        self.t = 1000.0
    def __call__(self):
        return self.t


class FakeService:
    """Controllable fair; counts fetches."""
    books = ("draftkings", "fanduel")

    def __init__(self):
        self.fair = 0.20
        self.fetches = 0
        self.lock = threading.Lock()

    def price_on_demand(self, book, game, legs):
        with self.lock:
            self.fetches += 1
        return OnDemandBookResult(book=book, fair=self.fair, route="partition",
                                  n_cells_priced=8, latency_sec=0.5)


class GW:
    def __init__(self):
        self.submits = []
    def submit_quote(self, rid, yb, nb):
        self.submits.append((rid, yb, nb))
        return f"q{len(self.submits)}"


class Src:
    def __init__(self):
        self.open_rfqs = [{"id": "r-od", "market_ticker": "COMBO-OD", "contracts": 1}]
    def poll(self):
        return list(self.open_rfqs)
    def get_market(self, t):
        return {}


def _setup(monkeypatch, tmp_path):
    import kalshi_mlb_mm.config as cfg
    import kalshi_mlb_mm.db as db
    import kalshi_mlb_mm.risk as risk
    from kalshi_mlb_mm import main
    monkeypatch.setattr(cfg, "DB_PATH", tmp_path / "e2e.duckdb")
    monkeypatch.setattr(cfg, "KILL_FILE", tmp_path / ".kill")
    # An open quote's worst-case exposure counts against the per-combo cap
    # (N7); raise it so the hysteresis REPLACE path is reachable (same
    # workaround as test_main_smoke's dedup test).
    monkeypatch.setattr(cfg, "MAX_COMBO_EXPOSURE_USD", 500.0)
    importlib.reload(db)
    db.init_database()
    clock = FakeClock()
    svc = FakeService()
    eng = OnDemandEngine(svc, now_fn=clock, autostart=False)
    monkeypatch.setattr(main, "_today_fills", lambda: [])
    monkeypatch.setattr(main, "_today_fills_by_game", lambda: [])
    monkeypatch.setattr(main, "_resolve_game_for_legs", lambda gl: "game1")
    monkeypatch.setattr(main, "_commence_time", lambda gid: None)
    monkeypatch.setattr(main, "_game_ref", lambda gid: GREF)
    monkeypatch.setattr(risk, "tipoff_ok", lambda ct, min_: True)
    monkeypatch.setattr(main, "_ENGINE", eng)
    monkeypatch.setattr(main, "_SCOPE_CACHE", {"COMBO-OD": (True, None, OD_LEGS)})
    monkeypatch.setattr(main, "_OD_RESULT_EMITTED", {})
    monkeypatch.setattr(main, "_PREV_BOOK_FAIR", {})   # isolate circuit breaker
    monkeypatch.setattr(main, "_leg_market_prices",
                        lambda legs: {"L": {"yes_bid": 0.5, "yes_ask": 0.52}})
    return main, db, clock, svc, eng


def test_full_feed_lifecycle(monkeypatch, tmp_path):
    main, db, clock, svc, eng = _setup(monkeypatch, tmp_path)
    src, gw = Src(), GW()

    # Tick 1: on_demand, nothing fresh -> pending skip, fetch queued.
    main._discovery_tick(src, gw, dry_run=False)
    assert gw.submits == []
    assert eng._queue_len() == 1

    # Fetch lands (worker drain), next tick quotes.
    eng._drain_once()
    assert svc.fetches == 2                       # one per book, shared fetch
    main._discovery_tick(src, gw, dry_run=False)
    assert len(gw.submits) == 1                   # quoted from live result

    # Within freshness + fair unchanged: hysteresis leaves the quote alone.
    clock.t += 2
    main._discovery_tick(src, gw, dry_run=False)
    assert len(gw.submits) == 1

    # Result ages out while the RFQ is still open -> the poll re-feeds.
    # Fair moves a little (beyond hysteresis, inside the 0.03 circuit-breaker
    # band) -> the landing REPLACES the resting quote.
    clock.t += QUOTE_FRESH_SEC
    svc.fair = 0.215
    main._discovery_tick(src, gw, dry_run=False)  # pending again + enqueue
    assert eng._queue_len() == 1
    assert len(gw.submits) == 1
    eng._drain_once()
    main._discovery_tick(src, gw, dry_run=False)  # landing re-quotes: replace
    assert len(gw.submits) == 2
    with db.connect(read_only=True) as con:
        statuses = dict(con.execute(
            "SELECT quote_id, status FROM live_quotes").fetchall())
    assert statuses["q1"] == "replaced" and statuses["q2"] == "open"

    # Fair JUMPS (>= BOOK_MOVE_CB_THRESHOLD) -> circuit breaker pulls the
    # resting quote instead of re-quoting (existing defense, now fed by
    # on-demand fairs).
    clock.t += QUOTE_FRESH_SEC + 1
    svc.fair = 0.35
    main._discovery_tick(src, gw, dry_run=False)  # pending + enqueue
    eng._drain_once()
    main._discovery_tick(src, gw, dry_run=False)  # CB fires on the landing
    assert len(gw.submits) == 2                   # no new quote
    with db.connect(read_only=True) as con:
        st_q2 = con.execute(
            "SELECT status FROM live_quotes WHERE quote_id='q2'").fetchone()[0]
    assert st_q2 == "cancelled"

    # RFQ leaves the poll -> feed stops: no new fetches ever again.
    src.open_rfqs = []
    clock.t += QUOTE_FRESH_SEC + 1
    fetches_before = svc.fetches
    for _ in range(5):
        clock.t += QUOTE_FRESH_SEC + 1
        main._discovery_tick(src, gw, dry_run=False)
        eng._drain_once()
    assert svc.fetches == fetches_before
    assert eng._queue_len() == 0


def test_mixed_grid_plus_on_demand_combo_prices_product(monkeypatch, tmp_path):
    main, db, clock, svc, eng = _setup(monkeypatch, tmp_path)
    # Combo: game A 2-leg grid + game B 3-leg on-demand.
    legs_a = [leg("KXMLBSPREAD-25JUN271905TEXLAA-LAA2", "yes"),
              leg("KXMLBTOTAL-25JUN271905TEXLAA-9", "yes")]
    legs_b = [leg("KXMLBSPREAD-25JUN272005SDLAD-LAD2", "yes"),
              leg("KXMLBTOTAL-25JUN272005SDLAD-8", "yes"),
              leg("KXMLBGAME-25JUN272005SDLAD-LAD", "yes")]
    rows = []
    for book in ("draftkings", "fanduel"):
        for combo, dec in (("Home Spread + Over", 3.10), ("Home Spread + Under", 4.20),
                           ("Away Spread + Over", 3.60), ("Away Spread + Under", 3.05)):
            rows.append(dict(game_id="gA", combo=combo, period="FG",
                             bookmaker=book, sgp_decimal=dec,
                             spread_line=-1.5, total_line=8.5))
    monkeypatch.setattr(main, "_resolve_game_for_legs",
                        lambda gl: "gA" if gl[0].game_id == EVT else "gB")
    monkeypatch.setattr(main, "_SCOPE_CACHE",
                        {"COMBO-MIX": (True, None, legs_a + legs_b)})

    class SrcMix:
        def poll(self):
            return [{"id": "r-mix", "market_ticker": "COMBO-MIX", "contracts": 1}]
        def get_market(self, t):
            return {}

    gw = GW()
    # #54 live-only: BOTH games ride the engine — no cached grid rows
    # exist to price from.
    svc.fair = 0.30           # 0.30 * 0.30 = 0.09, above MIN_FAIR_PROB
    main._discovery_tick(SrcMix(), gw, dry_run=False)     # pending (both games)
    assert gw.submits == []
    assert eng._queue_len() == 2, "one live job per game, grid included"
    eng._drain_once()
    eng._drain_once()
    main._discovery_tick(SrcMix(), gw, dry_run=False)     # both games priceable
    assert len(gw.submits) == 1
    # quoted fair = live consensus(gA) * live consensus(gB) = 0.30 * 0.30
    with db.connect(read_only=True) as con:
        blended = con.execute("SELECT blended_fair FROM live_quotes").fetchone()[0]
    assert abs(blended - 0.09) < 1e-9
