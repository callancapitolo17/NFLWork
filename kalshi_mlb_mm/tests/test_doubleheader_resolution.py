"""Game identity when the team pair is not unique.

Two Kalshi events can carry the same team pair: a doubleheader's two games
(Kalshi appends G1/G2 — KXMLBGAME-26SEP041410DETCLEG1 / ...1915DETCLEG2, live
2026-09-01) share it on the SAME day, and a series' consecutive games share it
a day apart. The Odds API ``events`` endpoint returns a multi-day window, so
both shapes are in ``mlb_target_lines`` at once.

``WHERE home_team=? AND away_team=? LIMIT 1`` answers either with whichever row
DuckDB hands back first — a wrong game, which #95 measured at 0.05-0.11 in
probability. So the pair NARROWS and the event suffix's own first pitch
IDENTIFIES, with ambiguity failing closed.

The invariant this file pins:
  * each game of a pair -> its OWN Odds-API id, never its twin's
  * two rows still indistinguishable inside the tolerance -> decline
  * leg-surface identity (full suffix) -> unchanged, still keeps them apart
"""
from datetime import datetime, timedelta, timezone

import pytest

from kalshi_common import legset

G1 = "26SEP041410DETCLEG1"      # 14:10 ET -> 18:10 UTC
G2 = "26SEP041915DETCLEG2"      # 19:15 ET -> 23:15 UTC
PLAIN = "26AUG271910MILNYM"
SINGLE = "26SEP051410DETCLE"    # same pairing, no doubleheader


def _total_leg(suffix: str):
    return legset.parse_leg({"event_ticker": f"KXMLBTOTAL-{suffix}",
                             "market_ticker": f"KXMLBTOTAL-{suffix}-9",
                             "side": "yes"})


def _market_db(tmp_path, monkeypatch, rows):
    """A MARKET_DB holding the given (game_id, home, away, commence_utc) rows.

    commence_time is NAIVE UTC here because that is what write_target_lines
    stores — the comparison only works if both sides agree on that.
    """
    import duckdb

    from kalshi_mlb_mm import config, main
    path = tmp_path / "market.duckdb"
    con = duckdb.connect(str(path))
    con.execute("CREATE TABLE mlb_target_lines (game_id VARCHAR, "
                "home_team VARCHAR, away_team VARCHAR, commence_time TIMESTAMP)")
    for row in rows:
        con.execute("INSERT INTO mlb_target_lines VALUES (?, ?, ?, ?)", list(row))
    con.close()
    monkeypatch.setattr(config, "MARKET_DB", path)
    main._RESOLVE_CACHE.clear()
    return main


CLE, DET = "Cleveland Guardians", "Detroit Tigers"


@pytest.fixture
def doubleheader_db(tmp_path, monkeypatch):
    return _market_db(tmp_path, monkeypatch, [
        ("dh-game-1", CLE, DET, datetime(2026, 9, 4, 18, 10)),
        ("dh-game-2", CLE, DET, datetime(2026, 9, 4, 23, 15)),
    ])


class TestStartTimeResolution:
    def test_an_ordinary_game_still_resolves(self, tmp_path, monkeypatch):
        main = _market_db(tmp_path, monkeypatch, [
            ("odds-api-id", CLE, DET, datetime(2026, 9, 5, 18, 10))])
        assert (main._resolve_game_for_legs_uncached([_total_leg(SINGLE)])
                == "odds-api-id")

    def test_each_doubleheader_game_resolves_to_its_own_id(self,
                                                           doubleheader_db):
        # The whole point: NOT the same id, and each one its own — this is what
        # unblocks quoting a combo that touches a doubleheader.
        main = doubleheader_db
        assert main._resolve_game_for_legs_uncached([_total_leg(G1)]) == "dh-game-1"
        assert main._resolve_game_for_legs_uncached([_total_leg(G2)]) == "dh-game-2"

    def test_the_two_games_never_share_an_id(self, doubleheader_db):
        main = doubleheader_db
        first = main._resolve_game_for_legs_uncached([_total_leg(G1)])
        second = main._resolve_game_for_legs_uncached([_total_leg(G2)])
        assert first is not None and second is not None and first != second

    def test_a_next_day_series_game_does_not_steal_the_id(self, tmp_path,
                                                          monkeypatch):
        """The live 2026-09-01 shape: same pair, one day apart.

        Resolving tonight's Kalshi event to tomorrow's row is a wrong game_id
        under the exposure ledgers AND a first pitch ~24h late under the
        tipoff gate.
        """
        main = _market_db(tmp_path, monkeypatch, [
            ("tonight", CLE, DET, datetime(2026, 9, 5, 18, 10)),
            ("tomorrow", CLE, DET, datetime(2026, 9, 6, 18, 10)),
        ])
        assert (main._resolve_game_for_legs_uncached([_total_leg(SINGLE)])
                == "tonight")

    def test_indistinguishable_rows_fail_closed(self, tmp_path, monkeypatch):
        main = _market_db(tmp_path, monkeypatch, [
            ("a", CLE, DET, datetime(2026, 9, 4, 18, 10)),
            ("b", CLE, DET, datetime(2026, 9, 4, 18, 25)),
        ])
        assert main._resolve_game_for_legs_uncached([_total_leg(G1)]) is None

    def test_a_game_the_schedule_does_not_carry_declines(self,
                                                         doubleheader_db):
        # Absent is a decline, not a nearest-neighbour match.
        main = doubleheader_db
        assert main._resolve_game_for_legs_uncached([_total_leg(PLAIN)]) is None

    def test_the_cache_keeps_the_two_games_apart(self, doubleheader_db):
        # _RESOLVE_CACHE keys on CanonicalLeg.game_id, which KEEPS the G1/G2
        # marker — if it did not, the first lookup would poison the second.
        main = doubleheader_db
        assert main._resolve_game_for_legs([_total_leg(G1)]) == "dh-game-1"
        assert main._resolve_game_for_legs([_total_leg(G2)]) == "dh-game-2"

    def test_leg_surface_identity_still_separates_them(self):
        assert _total_leg(G1).game_id != _total_leg(G2).game_id


class TestRfiDoubleheaderExclusion:
    """kalshi_rfi quotes one market per game off team-name book matching, so
    it drops doubleheaders outright. Its same-ET-day COUNT heuristic misses a
    pair whose first game already started and delisted; the suffix marker
    does not."""

    def _game(self, suffix, **kw):
        from kalshi_rfi.discovery import RfiGame, parse_suffix_start_utc
        return RfiGame(ticker=f"KXMLBRFI-{suffix}", suffix=suffix,
                       home_team=CLE,
                       away_team=DET,
                       commence_utc=parse_suffix_start_utc(suffix),
                       yes_bid_cents=40, yes_ask_cents=45, status="open",
                       **kw)

    def test_a_lone_surviving_doubleheader_game_is_still_dropped(self):
        from kalshi_rfi.discovery import drop_doubleheaders
        assert drop_doubleheaders([self._game(G2)]) == []

    def test_both_games_are_dropped(self):
        from kalshi_rfi.discovery import drop_doubleheaders
        assert drop_doubleheaders([self._game(G1), self._game(G2)]) == []

    def test_an_ordinary_game_is_kept(self):
        from kalshi_rfi.discovery import drop_doubleheaders
        kept = drop_doubleheaders([self._game(PLAIN)])
        assert [g.suffix for g in kept] == [PLAIN]


class TestAmbiguityIsWarnedOnce:
    """Resolution FAILURES are deliberately not cached (a transient DB lock
    must not stick until the next refresh), so the warning needs its own
    guard — the discovery tick re-enters every open RFQ every 2s, and evening
    peak runs ~100 RFQs/min."""

    def test_repeated_ambiguous_lookups_warn_once(self, tmp_path, monkeypatch,
                                                  caplog):
        import logging
        main = _market_db(tmp_path, monkeypatch, [
            ("a", CLE, DET, datetime(2026, 9, 4, 18, 10)),
            ("b", CLE, DET, datetime(2026, 9, 4, 18, 25)),
        ])
        main._AMBIGUOUS_WARNED.clear()
        with caplog.at_level(logging.WARNING, logger="kalshi_mlb_mm"):
            for _ in range(5):
                assert main._resolve_game_for_legs_uncached(
                    [_total_leg(G1)]) is None
        warnings = [r for r in caplog.records
                    if "resolve_game_ambiguous" in r.getMessage()]
        assert len(warnings) == 1
        main._AMBIGUOUS_WARNED.clear()

    def test_a_target_line_refresh_re_arms_the_warning(self, tmp_path,
                                                       monkeypatch):
        # The ambiguity may be gone after the refresh; if it is not, the next
        # cycle says so again rather than going quiet forever.
        from kalshi_common import sgp_runner
        main = _market_db(tmp_path, monkeypatch, [
            ("a", CLE, DET, datetime(2026, 9, 4, 18, 10)),
            ("b", CLE, DET, datetime(2026, 9, 4, 18, 25)),
        ])
        main._AMBIGUOUS_WARNED.clear()
        main._resolve_game_for_legs_uncached([_total_leg(G1)])
        assert main._AMBIGUOUS_WARNED
        monkeypatch.setattr(sgp_runner, "target_line_cycle", lambda **kw: [])
        main._target_line_tick()
        assert main._AMBIGUOUS_WARNED == set()


class TestTipoffGateFailsClosedOnAnyUnreadableClock:
    """Fix from the pre-merge review: the discovery tick used to drop None
    clocks and take min() of the rest, so a 2-game combo passed on its OTHER
    game's clock. The sweep's _quote_first_pitches already failed closed;
    the tick now does too."""

    def test_quote_first_pitches_is_none_if_any_game_is_unreadable(self):
        import json
        from kalshi_mlb_mm import main
        legs = [{"event_ticker": f"KXMLBTOTAL-{G1}",
                 "market_ticker": f"KXMLBTOTAL-{G1}-9", "side": "yes"},
                {"event_ticker": "KXMLBTOTAL-garbage",
                 "market_ticker": "KXMLBTOTAL-garbage-9", "side": "yes"}]
        assert main._quote_first_pitches(json.dumps(legs)) is None

    def test_quote_first_pitches_reads_both_games(self):
        import json
        from kalshi_mlb_mm import main
        legs = [{"event_ticker": f"KXMLBTOTAL-{G1}",
                 "market_ticker": f"KXMLBTOTAL-{G1}-9", "side": "yes"},
                {"event_ticker": f"KXMLBTOTAL-{PLAIN}",
                 "market_ticker": f"KXMLBTOTAL-{PLAIN}-9", "side": "yes"}]
        starts = main._quote_first_pitches(json.dumps(legs))
        assert starts is not None and len(starts) == 2
        assert all(s.tzinfo is not None for s in starts)


class TestDiscoveryTickTipoffFailsClosed:
    """The tick used to take min() of the clocks it COULD read, so a 2-game
    combo passed on its other game's clock. Now any unreadable clock skips."""

    _EVT_A = f"KXMLBGAME-{PLAIN}"
    _EVT_B = f"KXMLBGAME-{G1}"
    _LEGS = [{"market_ticker": f"KXMLBTOTAL-{PLAIN}-9",
              "event_ticker": f"KXMLBTOTAL-{PLAIN}", "side": "yes"},
             {"market_ticker": f"KXMLBTOTAL-{G1}-9",
              "event_ticker": f"KXMLBTOTAL-{G1}", "side": "yes"}]

    def _setup(self, monkeypatch, tmp_path, first_pitch):
        import importlib
        import kalshi_mlb_mm.config as cfg
        import kalshi_mlb_mm.db as db
        import kalshi_mlb_mm.risk as risk
        import kalshi_mlb_mm.router as router_mod
        from kalshi_mlb_mm import main
        from kalshi_mlb_mm.tests.conftest import FakeLiveEngine

        monkeypatch.setattr(cfg, "DB_PATH", tmp_path / "tick.duckdb")
        monkeypatch.setattr(cfg, "KILL_FILE", tmp_path / ".kill")
        importlib.reload(db)
        db.init_database()
        # tipoff_ok is waved through: the point is that the gate must never
        # REACH it with a clock it could not read.
        monkeypatch.setattr(risk, "tipoff_ok", lambda ct, m: True)
        monkeypatch.setattr(router_mod, "combo_fair_detail",
                            lambda *a, **k: (router_mod.ComboFair(0.55, 0.0, 1), "ok"))
        monkeypatch.setattr(main, "_first_pitch_utc", first_pitch)
        monkeypatch.setattr(main, "_PREV_BOOK_FAIR", {})
        monkeypatch.setattr(main, "_SCOPE_CACHE",
                            {"COMBO-XG": (True, "gA", self._LEGS)})
        monkeypatch.setattr(main, "_resolve_game_for_legs",
                            lambda gl: {PLAIN: "gA", G1: "gB"}.get(gl[0].game_id))
        monkeypatch.setattr(main, "_leg_market_prices",
                            lambda legs: {"L": {"yes_bid": 0.5, "yes_ask": 0.52}})
        monkeypatch.setattr(main, "_ENGINE", FakeLiveEngine())
        return main

    class _GW:
        def __init__(self):
            self.submits = []

        def submit_quote(self, *a):
            self.submits.append(a)
            return "qid-new"

    class _Src:
        def poll(self):
            return [{"id": "r-xg", "market_ticker": "COMBO-XG", "contracts": 1}]

        def get_market(self, t):
            return {}

    def test_one_unreadable_clock_skips_the_whole_combo(self, monkeypatch,
                                                        tmp_path):
        soon = datetime.now(timezone.utc) + timedelta(hours=1)
        main = self._setup(monkeypatch, tmp_path,
                           lambda gl: soon if gl[0].game_id == PLAIN else None)
        gw = self._GW()
        main._discovery_tick(self._Src(), gw, dry_run=False)
        assert gw.submits == [], "a combo with an unreadable clock was quoted"

    def test_both_clocks_readable_still_quotes(self, monkeypatch, tmp_path):
        # Control: the same harness quotes when every clock reads.
        soon = datetime.now(timezone.utc) + timedelta(hours=1)
        main = self._setup(monkeypatch, tmp_path, lambda gl: soon)
        gw = self._GW()
        main._discovery_tick(self._Src(), gw, dry_run=False)
        assert len(gw.submits) == 1
