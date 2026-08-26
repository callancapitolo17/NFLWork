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
            odds=dict(self.odds), n_legs_requested=len(legs),
            outcome=self.outcome)


class TestStructurePass:
    def test_one_fetch_serves_the_whole_ladder(self):
        # The cost model of the whole epic: work is O(games), not O(legs).
        service = FakeService(odds={0: (1.91, 1.91), 1: (1.91, 1.91)})
        rows, _counts, priced = structure.run_pass(
            "betmgm", service, [make_game()], band_min=1.005, band_max=1.20,
            max_req_per_sec=0)
        assert len(service.calls) == 1
        assert len(rows) == 2 and priced == 1

    def test_a_missing_alt_rung_does_not_unprice_the_others(self):
        # resolve_legs is all-or-nothing over the list it is handed, so legs
        # resolve one at a time. If that ever regressed, this test goes to
        # zero rows instead of two.
        service = FakeService(odds={0: (1.91, 1.91), 1: (1.91, 1.91)})
        rows, counts, _ = structure.run_pass(
            "betmgm", service, [make_game()], band_min=1.005, band_max=1.20,
            max_req_per_sec=0)
        assert len(rows) == 2
        assert counts.unresolved == 1        # the 12.5 rung

    def test_one_resolved_leg_devigs_against_its_opposite_price(self):
        # The book resolved only the OVER leg but published both prices.
        service = FakeService(odds={0: (1.91, 1.91)})
        rows, counts, _ = structure.run_pass(
            "betmgm", service, [make_game()], band_min=1.005, band_max=1.20,
            max_req_per_sec=0)
        assert {r.side for r in rows} == {"over", "under"}
        assert counts.one_sided == 0

    def test_a_genuinely_one_sided_rung_is_excluded(self):
        service = FakeService(odds={0: (1.91, None)})
        rows, counts, _ = structure.run_pass(
            "betmgm", service, [make_game()], band_min=1.005, band_max=1.20,
            max_req_per_sec=0)
        assert rows == [] and counts.one_sided == 1

    def test_rows_carry_the_structure_fetch_time(self):
        service = FakeService(odds={0: (1.91, 1.91)})
        rows, _, _ = structure.run_pass(
            "betmgm", service, [make_game()], band_min=1.005, band_max=1.20,
            max_req_per_sec=0)
        assert all(r.built_at == FETCHED_AT for r in rows)

    @pytest.mark.parametrize("outcome", ["no_event", "transport_error"])
    def test_a_dark_book_is_one_unmatched_game_not_a_ladder_of_misses(
            self, outcome):
        # Charging it per rung would make one dark book look like a coverage
        # collapse across every line in the slate.
        service = FakeService(outcome=outcome)
        rows, counts, priced = structure.run_pass(
            "betmgm", service, [make_game()], band_min=1.005, band_max=1.20,
            max_req_per_sec=0)
        assert rows == [] and priced == 0
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
