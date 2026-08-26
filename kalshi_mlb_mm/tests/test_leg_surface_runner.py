"""Structure-route pass, pass accounting and the DuckDB mirror (issue #96).

Two failures these guard against. First, a DARK book must not read as a book
that simply offers none of these alt rungs — #95 showed rungs legitimately
vanish far from first pitch, so the two look identical in the row count and
only the outcome tells them apart. Second, a FAILED pass must not publish an
empty slice: that would blank a book's prices on a transient blip instead of
letting them age out under the staleness gate.
"""
from datetime import datetime, timezone

import pytest

from kalshi_common.legset import CanonicalLeg
from kalshi_common.sgp_service import StructureLegOdds
from kalshi_mlb_mm.leg_surface import db, structure
from kalshi_mlb_mm.leg_surface.runner import SurfaceIngest
from kalshi_mlb_mm.leg_surface.slate import SurfaceGame

FETCHED_AT = datetime(2026, 8, 26, 1, 30, tzinfo=timezone.utc)


def make_game():
    gid = "26AUG252138CLELAA"
    return SurfaceGame(
        game_id=gid, home_team="Los Angeles Angels",
        away_team="Cleveland Guardians",
        start_utc=datetime(2026, 8, 26, 1, 38),
        legs=(CanonicalLeg(gid, "total", 8.5, "over"),
              CanonicalLeg(gid, "total", 8.5, "under"),
              CanonicalLeg(gid, "total", 12.5, "over"),
              CanonicalLeg(gid, "total", 12.5, "under")))


class FakeService:
    """One structure fetch per game, then local resolution — the real seam."""

    def __init__(self, odds=None, outcome="ok"):
        self.odds = odds if odds is not None else {}
        self.outcome = outcome
        self.calls = []

    def structure_leg_odds(self, book, game, legs):
        self.calls.append((book, game.game_id, len(legs)))
        return StructureLegOdds(
            book=book, fetched_at=FETCHED_AT if self.outcome == "ok" else None,
            odds=dict(self.odds), outcome=self.outcome)


class TestStructurePass:
    def test_one_fetch_serves_the_whole_ladder(self):
        # The cost model of the whole epic: work is O(games), not O(legs).
        service = FakeService(odds={0: (1.91, 1.91), 1: (1.91, 1.91)})
        rows, _counts, priced, _dark = structure.run_pass(
            "betmgm", service, [make_game()], band_min=1.005, band_max=1.20,
            max_req_per_sec=0)
        assert len(service.calls) == 1
        assert len(rows) == 2 and priced == 1

    def test_a_missing_alt_rung_does_not_unprice_the_others(self):
        # resolve_legs is all-or-nothing over the list it is handed, so legs
        # resolve one at a time. If that ever regressed, this test goes to
        # zero rows instead of two.
        service = FakeService(odds={0: (1.91, 1.91), 1: (1.91, 1.91)})
        rows, counts, _priced, _dark = structure.run_pass(
            "betmgm", service, [make_game()], band_min=1.005, band_max=1.20,
            max_req_per_sec=0)
        assert len(rows) == 2
        assert counts.unresolved == 1        # the 12.5 rung

    def test_one_resolved_leg_devigs_against_its_opposite_price(self):
        # The book resolved only the OVER leg but published both prices.
        service = FakeService(odds={0: (1.91, 1.91)})
        rows, counts, _priced, _dark = structure.run_pass(
            "betmgm", service, [make_game()], band_min=1.005, band_max=1.20,
            max_req_per_sec=0)
        assert {r.side for r in rows} == {"over", "under"}
        assert counts.one_sided == 0

    def test_a_genuinely_one_sided_rung_is_excluded(self):
        service = FakeService(odds={0: (1.91, None)})
        rows, counts, _priced, _dark = structure.run_pass(
            "betmgm", service, [make_game()], band_min=1.005, band_max=1.20,
            max_req_per_sec=0)
        assert rows == [] and counts.one_sided == 1

    def test_rows_carry_the_structure_fetch_time(self):
        service = FakeService(odds={0: (1.91, 1.91)})
        rows, _counts, _priced, _dark = structure.run_pass(
            "betmgm", service, [make_game()], band_min=1.005, band_max=1.20,
            max_req_per_sec=0)
        assert all(r.built_at == FETCHED_AT for r in rows)

    @pytest.mark.parametrize("outcome", ["no_event", "transport_error"])
    def test_a_dark_book_is_one_unmatched_game_not_a_ladder_of_misses(
            self, outcome):
        # Charging it per rung would make one dark book look like a coverage
        # collapse across every line in the slate.
        service = FakeService(outcome=outcome)
        rows, counts, priced, dark = structure.run_pass(
            "betmgm", service, [make_game()], band_min=1.005, band_max=1.20,
            max_req_per_sec=0)
        assert rows == [] and priced == 0
        assert dark == (1 if outcome == "transport_error" else 0)
        assert counts.game_unmatched == 1 and counts.unresolved == 0


class TestPassAccounting:
    def make_ingest(self, tmp_path, monkeypatch, service):
        from kalshi_mlb_mm import config
        monkeypatch.setattr(config, "SURFACE_DB",
                            tmp_path / "surface.duckdb")
        db.init_database()
        ingest = SurfaceIngest(service=service, books_structure=("betmgm",),
                               books_singles=())
        ingest._slate = [make_game()]
        return ingest

    def test_a_successful_pass_publishes_and_mirrors(self, tmp_path,
                                                     monkeypatch):
        service = FakeService(odds={0: (1.91, 1.91)})
        ingest = self.make_ingest(tmp_path, monkeypatch, service)
        result = ingest.structure_pass("betmgm")
        assert result.legs_written == 2 and result.rungs_priced == 1
        assert ingest.surface.row_count() == 2
        with db.connect(read_only=True) as con:
            assert con.execute(
                "SELECT COUNT(*) FROM mlb_leg_surface").fetchone()[0] == 2
            assert con.execute(
                "SELECT COUNT(*) FROM surface_refresh_log").fetchone()[0] == 1

    def test_a_failed_pass_keeps_the_previous_rows(self, tmp_path,
                                                   monkeypatch):
        # A blip must not blank a book. The rows age out under #99's gate
        # instead, which is a decline the operator can count.
        service = FakeService(odds={0: (1.91, 1.91)})
        ingest = self.make_ingest(tmp_path, monkeypatch, service)
        ingest.structure_pass("betmgm")

        class Boom:
            def structure_leg_odds(self, *a, **kw):
                raise RuntimeError("book exploded")

        ingest._service = Boom()
        result = ingest.structure_pass("betmgm")
        assert result.error_class == "RuntimeError"
        assert ingest.surface.row_count() == 2       # untouched
        with db.connect(read_only=True) as con:
            assert con.execute(
                "SELECT COUNT(*) FROM mlb_leg_surface").fetchone()[0] == 2

    def test_no_duplicate_keys_after_repeated_passes(self, tmp_path,
                                                     monkeypatch):
        # Acceptance criterion. Uniqueness is structural (DELETE the slice,
        # re-insert it) rather than a constraint, because `line` is
        # legitimately NULL for moneyline and DuckDB PKs reject NULL.
        service = FakeService(odds={0: (1.91, 1.91)})
        ingest = self.make_ingest(tmp_path, monkeypatch, service)
        for _ in range(3):
            ingest.structure_pass("betmgm")
        with db.connect(read_only=True) as con:
            dupes = con.execute(
                "SELECT COUNT(*) FROM (SELECT book, game_id, period, "
                "market_type, line, side FROM mlb_leg_surface "
                "GROUP BY ALL HAVING COUNT(*) > 1)").fetchone()[0]
            total = con.execute(
                "SELECT COUNT(*) FROM mlb_leg_surface").fetchone()[0]
        assert dupes == 0 and total == 2

    def test_exclusion_reasons_are_countable_in_sql(self, tmp_path,
                                                    monkeypatch):
        service = FakeService(odds={0: (2.10, 2.10)})     # crossed
        ingest = self.make_ingest(tmp_path, monkeypatch, service)
        ingest.structure_pass("betmgm")
        with db.connect(read_only=True) as con:
            crossed, unresolved = con.execute(
                "SELECT SUM(n_crossed), SUM(n_unresolved) "
                "FROM surface_refresh_log").fetchone()
        assert crossed == 1 and unresolved == 1

    def test_refresh_log_prunes_to_its_retention_window(self, tmp_path,
                                                        monkeypatch):
        service = FakeService(odds={0: (1.91, 1.91)})
        ingest = self.make_ingest(tmp_path, monkeypatch, service)
        ingest.structure_pass("betmgm")
        with db.connect() as con:
            con.execute("UPDATE surface_refresh_log "
                        "SET started_at = started_at - INTERVAL 48 HOUR")
        assert db.prune_refresh_log(retention_hours=24) == 1


def test_slate_refresh_keeps_the_previous_slate_on_an_empty_result(
        tmp_path, monkeypatch):
    # An empty Kalshi response is far more often a blip than an off-day, and
    # dropping the slate blanks every book on the next pass.
    from kalshi_mlb_mm import config
    from kalshi_mlb_mm.leg_surface import runner as runner_mod
    monkeypatch.setattr(config, "SURFACE_DB", tmp_path / "surface.duckdb")
    db.init_database()
    ingest = SurfaceIngest(service=FakeService())
    ingest._slate = [make_game()]
    monkeypatch.setattr(runner_mod.slate, "discover_slate",
                        lambda **kw: [])
    assert ingest.refresh_slate() == 0
    assert len(ingest.current_slate()) == 1


def test_worker_honours_its_cadence(tmp_path, monkeypatch):
    """Regression: the first live run had every worker spinning at full speed.

    ``threading.Event.wait()`` returns IMMEDIATELY when the event is set, so
    sleeping on a set "running" flag is a no-op — FanDuel's 45s singles scrape
    was firing every 2.5s, i.e. ~18x the intended book traffic. The loop
    sleeps on a STOP event (clear while running) instead.
    """
    import threading
    import time

    from kalshi_mlb_mm import config
    monkeypatch.setattr(config, "SURFACE_DB", tmp_path / "surface.duckdb")
    db.init_database()

    ingest = SurfaceIngest(service=FakeService())
    ingest._slate = [make_game()]
    ingest._slate_ready.set()
    starts = []
    monkeypatch.setattr(ingest, "structure_pass",
                        lambda book: starts.append(time.monotonic()))

    ingest._stop.clear()
    worker = threading.Thread(
        target=ingest._book_loop, args=("betmgm", "structure", 0.25),
        daemon=True)
    worker.start()
    time.sleep(1.0)
    ingest._stop.set()
    worker.join(timeout=2)

    assert len(starts) <= 6            # a spinning loop lands in the hundreds
    gaps = [b - a for a, b in zip(starts, starts[1:])]
    assert all(gap >= 0.2 for gap in gaps), gaps


class TestRouteOwnershipIsExclusive:
    """Regression: the first live run wrote 44 duplicate keys.

    FanDuel's structure route was pricing its FG totals even though its
    singles route owns them, so the same surface key landed from two slices —
    duplicate rows in the mirror, and a coin flip between two prices in
    memory. Ownership is now read from ONE place by both routes.
    """

    def make_game(self):
        gid = "26AUG252138CLELAA"
        return SurfaceGame(
            gid, "Los Angeles Angels", "Cleveland Guardians",
            datetime(2026, 8, 26, 1, 38),
            (CanonicalLeg(gid, "total", 8.5, "over"),
             CanonicalLeg(gid, "total", 8.5, "under"),
             CanonicalLeg(gid, "spread", -1.5, "home"),
             CanonicalLeg(gid, "spread", -1.5, "away")))

    def test_structure_route_skips_rungs_the_singles_route_owns(self):
        service = FakeService(odds={i: (1.91, 1.91) for i in range(4)})
        rows, counts, _priced, _dark = structure.run_pass(
            "fanduel", service, [self.make_game()],
            skip_markets={("total", "FG")}, band_min=1.005, band_max=1.20,
            max_req_per_sec=0)
        assert {r.market_type for r in rows} == {"spread"}
        # Skipped, not missed: this route was never asked for those rungs.
        assert counts.unresolved == 0

    def test_fanduel_is_the_only_book_with_a_split(self):
        ingest = SurfaceIngest()
        assert ingest.singles_markets("fanduel") == {("total", "FG"),
                                                     ("total", "F5")}
        # BetMGM runs structure only, so it must skip nothing.
        assert ingest.singles_markets("betmgm") == set()

    def test_both_fanduel_routes_together_write_no_duplicate_key(
            self, tmp_path, monkeypatch):
        """The end-to-end form of the bug: run BOTH FD routes, check the
        mirror. This is the query that caught it on live data."""
        from kalshi_mlb_mm import config
        from kalshi_mlb_mm.leg_surface import singles as singles_mod
        monkeypatch.setattr(config, "SURFACE_DB", tmp_path / "surface.duckdb")
        db.init_database()

        game = self.make_game()
        ingest = SurfaceIngest(
            service=FakeService(odds={i: (1.91, 1.91) for i in range(4)}),
            books_structure=("fanduel",), books_singles=("fanduel",))
        ingest._slate = [game]
        monkeypatch.setattr(singles_mod, "_scrape", lambda book: [{
            "game_id": "fd1", "game_start_time": game.start_utc,
            "home_team": game.home_team, "away_team": game.away_team,
            "period": "FG", "fetch_time": datetime.now(timezone.utc),
            "total": 8.5, "over_price": -110, "under_price": -110,
            "home_ml": None, "away_ml": None, "home_spread": -1.5,
            "home_spread_price": -110, "away_spread": 1.5,
            "away_spread_price": -110}])

        ingest.structure_pass("fanduel")
        ingest.singles_pass("fanduel")
        with db.connect(read_only=True) as con:
            dupes = con.execute(
                "SELECT COUNT(*) FROM (SELECT book, game_id, period, "
                "market_type, line, side FROM mlb_leg_surface "
                "GROUP BY ALL HAVING COUNT(*) > 1)").fetchone()[0]
        assert dupes == 0
        # And the split actually happened: spreads from structure, totals
        # from singles, one row each.
        leg_total = CanonicalLeg(game.game_id, "total", 8.5, "over")
        leg_spread = CanonicalLeg(game.game_id, "spread", -1.5, "home")
        assert ingest.surface.get("fanduel", leg_total).route == "singles"
        assert ingest.surface.get("fanduel", leg_spread).route == "structure"


def test_a_fully_dark_book_keeps_its_previous_rows(tmp_path, monkeypatch):
    """Regression: a Caesars transport blip blanked its whole slice.

    Every game transport-failed, which is not an exception, so the pass
    "succeeded" with zero rows and published an empty slice — losing a book
    instantly on a wobble. A dark pass now publishes nothing and the rows age
    out under #99's gate instead, which is a countable decline.
    """
    from kalshi_mlb_mm import config
    monkeypatch.setattr(config, "SURFACE_DB", tmp_path / "surface.duckdb")
    db.init_database()
    ingest = SurfaceIngest(service=FakeService(odds={0: (1.91, 1.91)}),
                           books_structure=("caesars",), books_singles=())
    ingest._slate = [make_game()]
    ingest.structure_pass("caesars")
    assert ingest.surface.row_count() == 2

    ingest._service = FakeService(outcome="transport_error")
    result = ingest.structure_pass("caesars")
    assert result.error_class == "book_dark"
    assert ingest.surface.row_count() == 2


def test_a_book_that_simply_stops_offering_a_rung_does_lose_it(tmp_path,
                                                               monkeypatch):
    # The mirror image of the test above: 'no_event' is the book answering,
    # so its rows MUST clear rather than rest at their last price forever.
    from kalshi_mlb_mm import config
    monkeypatch.setattr(config, "SURFACE_DB", tmp_path / "surface.duckdb")
    db.init_database()
    ingest = SurfaceIngest(service=FakeService(odds={0: (1.91, 1.91)}),
                           books_structure=("betmgm",), books_singles=())
    ingest._slate = [make_game()]
    ingest.structure_pass("betmgm")

    ingest._service = FakeService(outcome="no_event")
    result = ingest.structure_pass("betmgm")
    assert result.error_class is None
    assert ingest.surface.row_count() == 0


def test_concurrent_workers_share_one_sgp_service(tmp_path, monkeypatch):
    """Four structure workers wake together on the first slate.

    A bare `if self._service is None: build()` builds FOUR services — three
    orphaned with their per-book HTTP clients never closed, and which one wins
    decided by whichever thread assigns last. Double-checked under a lock.
    """
    import threading

    from kalshi_mlb_mm import config
    monkeypatch.setattr(config, "SURFACE_DB", tmp_path / "surface.duckdb")
    db.init_database()

    import time as _time

    built = []
    barrier = threading.Barrier(4)
    ingest = SurfaceIngest(books_structure=("fanduel", "betmgm", "novig",
                                            "caesars"), books_singles=())

    def slow_build():
        # Widen the window a real SGPService construction opens (six HTTP
        # clients) so an unlocked check-then-build would reliably lose.
        _time.sleep(0.05)
        service = FakeService()
        built.append(service)
        return service

    def racer():
        barrier.wait(timeout=5)      # all four enter _ensure_service together
        ingest._ensure_service()

    monkeypatch.setattr(ingest, "_build_service", slow_build)
    threads = [threading.Thread(target=racer) for _ in range(4)]
    for t in threads:
        t.start()
    for t in threads:
        t.join(timeout=10)

    # The barrier proves all four raced; the lock proves only one built.
    assert len(built) == 1
    assert ingest._service is built[0]
