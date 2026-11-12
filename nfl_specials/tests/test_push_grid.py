"""Whole-number 1H lines: the push grid fitted from the partitions either
side, Wagerzon's push payout, and the effective fair. No network."""
from datetime import datetime, timezone

import pytest
from scipy.stats import norm

from nfl_specials import push_grid
from nfl_specials.books import BookGame, LegMarket, Outcome
from nfl_specials.halves_board import HalfBookFair, HalfLine, price_half_at_book
from nfl_specials.pricing import american_to_decimal
from nfl_specials.wz_halves import CORNERS, WzHalfGame, WzParlayPrice

A, P, H, O, U = push_grid.AWAY, push_grid.PUSH, push_grid.HOME, push_grid.OVER, push_grid.UNDER
# Away -6, total 26 (both even, so both can land together).
TRUE_GRID = {(A, O): 0.20, (A, P): 0.02, (A, U): 0.18,
             (P, O): 0.03, (P, P): 0.01, (P, U): 0.03,
             (H, O): 0.22, (H, P): 0.02, (H, U): 0.29}
BOOK_GAME = BookGame("evt", "NYG", "WAS", "2099-01-01T18:00:00Z")


def _game(**overrides) -> WzHalfGame:
    fields = dict(league="NFL", wz_game_id=7, away="WAS", home="NYG",
                  game_start_time=datetime(2099, 1, 1, 18, tzinfo=timezone.utc),
                  away_spread=-6.0, away_spread_american=-110, home_spread_american=-110,
                  total=26.0, over_american=-110, under_american=-110)
    fields.update(overrides)
    return WzHalfGame(**fields)


class GridBook:
    """Prices any 1H alt line off a true outcome grid, shifted by `vig_z` in
    probit space (what probit devig removes exactly)."""
    name = "DraftKings"

    def __init__(self, grid: dict, game: WzHalfGame, vig_z: float = 0.1, missing_total: float | None = None):
        self.grid, self.game, self.vig_z, self.missing_total = grid, game, vig_z, missing_total

    def leg_market(self, book_game, role, leg):
        if leg.kind == "total" and leg.line == self.missing_total:
            return None
        sides = ("away", "home") if leg.kind == "spread" else ("over", "under")
        outcomes = tuple(Outcome((leg.kind, side, leg.line), side) for side in sides)
        return LegMarket.single(outcomes[0], outcomes)

    def price(self, book_game, refs):
        (_, spread_side, spread_line), (_, total_side, total_line) = refs
        spreads = push_grid.spread_band(spread_side, spread_line, self.game.away_spread)
        totals = push_grid.total_band(total_side, total_line, self.game.total)
        prob = sum(p for (s, t), p in self.grid.items() if s in spreads and t in totals)
        return 1.0 / norm.cdf(norm.ppf(prob) + self.vig_z)


def _partitions(grid: dict, game: WzHalfGame, noise: dict | None = None) -> dict:
    """Exact devigged corner fairs per bracket pair, with optional per-(pair, corner) noise."""
    noise = noise or {}
    out = {}
    for pair in push_grid.corner_bracket_pairs(game):
        fairs = {}
        for corner in CORNERS:
            spread_side, total_side = corner.split("_")
            spreads = push_grid.spread_band(spread_side, pair[0], game.away_spread)
            totals = push_grid.total_band(total_side, pair[1], game.total)
            fairs[corner] = sum(p for (s, t), p in grid.items() if s in spreads and t in totals)
            fairs[corner] += noise.get((pair, corner), 0.0)
        out[pair] = fairs
    return out


def test_brackets_and_bands():
    assert push_grid.bracket_lines(-6.0) == (-6.5, -5.5) and push_grid.bracket_lines(26.5) == (26.5,)
    # Away -6.5 needs a 7-point win: only Wagerzon's away cover. Away -5.5 wins on a 6-point win too.
    assert push_grid.spread_band(A, -6.5, -6.0) == {A} and push_grid.spread_band(A, -5.5, -6.0) == {A, P}
    assert push_grid.spread_band(H, -6.5, -6.0) == {P, H}
    assert push_grid.total_band(O, 25.5, 26.0) == {O, P} and push_grid.total_band(U, 26.5, 26.0) == {P, U}


def test_fit_recovers_the_grid_from_both_whole_numbers():
    fit = push_grid.fit_push_grid(_partitions(TRUE_GRID, _game()), -6.0, 26.0)
    assert fit.reason is None
    assert fit.grid.cells == pytest.approx(TRUE_GRID, abs=1e-9)


def test_one_whole_number_leg_fits_a_three_by_two_grid():
    grid = {(A, O): 0.22, (A, U): 0.18, (P, O): 0.04, (P, U): 0.03, (H, O): 0.24, (H, U): 0.29}
    game = _game(total=26.5)
    fit = push_grid.fit_push_grid(_partitions(grid, game), -6.0, 26.5)
    assert fit.grid.cells == pytest.approx(grid, abs=1e-9)


def test_odd_and_even_lines_cannot_both_push():
    assert not push_grid.double_push_possible(-3.0, 24.0)
    assert push_grid.double_push_possible(-3.0, 27.0) and push_grid.double_push_possible(-6.0, 26.0)
    grid = dict(TRUE_GRID)
    grid[(A, O)] += grid.pop((P, P))
    fit = push_grid.fit_push_grid(_partitions(grid, _game(away_spread=-3.0)), -3.0, 26.0)
    assert (P, P) not in fit.grid.cells and fit.grid.cells == pytest.approx(grid, abs=1e-9)


def test_a_negative_push_slice_skips_the_game():
    # Prices that agree with each other but imply a 6-point win with the over is below zero.
    grid = dict(TRUE_GRID)
    grid[(A, P)], grid[(A, O)] = -0.02, grid[(A, O)] + grid[(A, P)] + 0.02
    fit = push_grid.fit_push_grid(_partitions(grid, _game()), -6.0, 26.0)
    assert fit.grid is None and fit.reason == "implied away/push chance -0.020 is negative"


def test_partitions_that_disagree_skip_the_game():
    # One whole-number leg: two partitions, six unknowns, eight equations. Shifting
    # mass between two cells the other partition pins cannot be fitted.
    grid = {(A, O): 0.22, (A, U): 0.18, (P, O): 0.04, (P, U): 0.03, (H, O): 0.24, (H, U): 0.29}
    game = _game(total=26.5)
    noise = {((-6.5, 26.5), "away_over"): 0.05, ((-6.5, 26.5), "away_under"): -0.05}
    fit = push_grid.fit_push_grid(_partitions(grid, game, noise), -6.0, 26.5)
    assert fit.grid is None and fit.reason.startswith("nearby lines disagree")


def test_rounding_sized_noise_is_clamped_not_skipped():
    grid = dict(TRUE_GRID)
    grid[(A, P)], grid[(A, O)] = 0.0, grid[(A, O)] + grid[(A, P)]
    noise = {((-6.5, 25.5), "away_over"): -0.002, ((-6.5, 25.5), "away_under"): 0.002}
    fit = push_grid.fit_push_grid(_partitions(grid, _game(), noise), -6.0, 26.0)
    assert fit.grid is not None and min(fit.grid.cells.values()) >= 0.0
    assert sum(fit.grid.cells.values()) == pytest.approx(1.0)


def test_an_implausible_push_chance_skips_the_game():
    grid = {(A, O): 0.15, (A, U): 0.15, (P, O): 0.15, (P, U): 0.15, (H, O): 0.2, (H, U): 0.2}
    fit = push_grid.fit_push_grid(_partitions(grid, _game(total=26.5)), -6.0, 26.5)
    assert fit.grid is None and "spread push" in fit.reason
    grid = {(A, O): 0.2, (A, P): 0.1, (A, U): 0.2, (H, O): 0.2, (H, P): 0.1, (H, U): 0.2}
    fit = push_grid.fit_push_grid(_partitions(grid, _game(away_spread=-6.5)), -6.5, 26.0)
    assert fit.grid is None and "total push" in fit.reason


def test_positive_pickem_and_half_spread_bands():
    assert push_grid.spread_band(A, 2.5, 3.0) == {A} and push_grid.spread_band(A, 3.5, 3.0) == {A, P}
    assert push_grid.spread_band(H, -0.5, 0.0) == {P, H} and push_grid.spread_band(A, 0.5, 0.0) == {A, P}
    for away_spread, total in ((3.0, 27.0), (0.0, 26.0), (-6.5, 26.0)):
        grid = dict(TRUE_GRID) if push_grid.is_whole_number(away_spread) else {
            (A, O): 0.23, (A, P): 0.02, (A, U): 0.21, (H, O): 0.22, (H, P): 0.03, (H, U): 0.29}
        game = _game(away_spread=away_spread, total=total)
        fit = push_grid.fit_push_grid(_partitions(grid, game), away_spread, total)
        assert fit.grid.cells == pytest.approx(grid, abs=1e-9), away_spread


def test_wagerzon_pays_the_other_leg_alone_on_a_push():
    game = _game(away_spread_american=-120, home_spread_american=110, over_american=100, under_american=-130)
    grid = push_grid.PushGrid(TRUE_GRID)
    expected = {
        "away_over": 0.02 * american_to_decimal(-120) + 0.03 * american_to_decimal(100) + 0.01,
        "away_under": 0.02 * american_to_decimal(-120) + 0.03 * american_to_decimal(-130) + 0.01,
        "home_over": 0.02 * american_to_decimal(110) + 0.03 * american_to_decimal(100) + 0.01,
        "home_under": 0.02 * american_to_decimal(110) + 0.03 * american_to_decimal(-130) + 0.01,
    }
    assert {c: push_grid.push_return(grid, c, game) for c in CORNERS} == pytest.approx(expected)
    assert push_grid.corner_outcomes(grid, "home_under") == (0.29, 0.02, 0.03, 0.01)


def test_book_prices_the_four_bracket_partitions_and_the_fair_is_effective():
    game = _game()
    fair = price_half_at_book(GridBook(TRUE_GRID, game), BOOK_GAME, game)
    assert fair.reason is None
    assert fair.corner_fairs["away_over"] == pytest.approx(0.20, abs=1e-6)
    line = HalfLine(game, "priced", book_fairs={"DraftKings": fair, "BetMGM": fair},
                    wz_prices={c: WzParlayPrice(3.6, 500) for c in CORNERS})
    push_back = 0.02 * american_to_decimal(-110) + 0.03 * american_to_decimal(-110) + 0.01
    assert line.fair_prob("away_over") == pytest.approx(0.20 + push_back / 3.6, abs=1e-6)
    # EV of the effective binary at 3.6 equals the real bet's EV.
    assert line.fair_prob("away_over") * 3.6 - 1 == pytest.approx(0.20 * 3.6 + push_back - 1, abs=1e-6)


def test_a_missing_bracket_line_prices_nothing_and_no_wz_price_no_fair():
    game = _game()
    fair = price_half_at_book(GridBook(TRUE_GRID, game, missing_total=26.5), BOOK_GAME, game)
    assert fair.corner_fairs is None and fair.reason == "no SGP 1H total at 26.5"
    priced = HalfBookFair("DraftKings", {c: 0.2 for c in CORNERS}, 1.2, push_returns={c: 0.05 for c in CORNERS})
    assert priced.effective_fair("away_over", None) is None


def test_a_three_way_half_point_spread_adds_the_tie_cells():
    """A ±0.5 spread read off a 3-way 1H winner wins on the team OR the tie:
    the corner is the sum of two cells, not the last one priced."""
    three_way = {("away", "over"): 0.25, ("away", "under"): 0.20, ("tie", "over"): 0.05,
                 ("tie", "under"): 0.08, ("home", "over"): 0.18, ("home", "under"): 0.24}

    class ThreeWayBook:
        name = "FanDuel"

        def leg_market(self, book_game, role, leg):
            if leg.kind == "spread":
                away, tie, home = (Outcome(x, x) for x in ("away", "tie", "home"))
                return LegMarket(winning=frozenset([away]), group=(away, tie, home))
            over, under = Outcome("over", "over"), Outcome("under", "under")
            return LegMarket.single(over, (over, under))

        def price(self, book_game, refs):
            return 1.0 / norm.cdf(norm.ppf(three_way[refs]) + 0.1)

    fair = price_half_at_book(ThreeWayBook(), BOOK_GAME, _game(away_spread=-0.5, total=26.5))
    assert fair.corner_fairs == pytest.approx({"away_over": 0.25, "away_under": 0.20,
                                               "home_over": 0.23, "home_under": 0.32}, abs=1e-6)
    # +0.5: the away team wins on its own outcome OR the tie.
    class PlusHalfBook(ThreeWayBook):
        def leg_market(self, book_game, role, leg):
            if leg.kind == "spread":
                away, tie, home = (Outcome(x, x) for x in ("away", "tie", "home"))
                return LegMarket(winning=frozenset([away, tie]), group=(away, tie, home))
            return super().leg_market(book_game, role, leg)

    fair = price_half_at_book(PlusHalfBook(), BOOK_GAME, _game(away_spread=0.5, total=26.5))
    assert fair.corner_fairs == pytest.approx({"away_over": 0.30, "away_under": 0.28,
                                               "home_over": 0.18, "home_under": 0.24}, abs=1e-6)


def test_push_return_is_stored_and_an_old_file_gains_the_column(tmp_path):
    import duckdb
    from nfl_specials.store import Store

    db = tmp_path / "old.duckdb"
    con = duckdb.connect(str(db))
    con.execute("""CREATE TABLE half_parlay_quotes (quoted_at TIMESTAMPTZ, league VARCHAR, wz_game_id BIGINT,
                   away_team VARCHAR, home_team VARCHAR, game_start_time TIMESTAMPTZ, away_spread DOUBLE,
                   total DOUBLE, corner VARCHAR, book VARCHAR, fair_prob DOUBLE, overround DOUBLE,
                   wz_decimal DOUBLE, reason VARCHAR)""")
    con.close()
    store = Store(db_path=db)
    game = _game()
    fair = price_half_at_book(GridBook(TRUE_GRID, game), BOOK_GAME, game)
    from nfl_specials.halves_board import _quote_rows
    store.append_half_quotes(_quote_rows(HalfLine(game, "priced"), datetime.now(timezone.utc), "DraftKings", fair=fair))
    con = duckdb.connect(str(db))
    stored = dict(con.execute("SELECT corner, push_return FROM half_parlay_quotes").fetchall())
    con.close()
    assert stored == pytest.approx(fair.push_returns)


def test_real_betmgm_brackets_fit_within_rounding():
    """LA @ ARI 1H, LA -4.5, total 20, BetMGM 2026-11-12 (0.05 price steps):
    the two direct reads of each push cell agree, so the fit holds."""
    from kalshi_common.fair_value import _probit_devig_n as probit_devig_n
    decimals = {(-4.5, 19.5): (2.70, 4.40, 3.90, 3.20), (-4.5, 20.5): (3.25, 3.50, 4.60, 2.85)}
    partitions = {pair: dict(zip(CORNERS, (float(p) for p in probit_devig_n([1 / d for d in prices]))))
                  for pair, prices in decimals.items()}
    fit = push_grid.fit_push_grid(partitions, -4.5, 20.0)
    assert fit.reason is None and fit.residual < 0.002
    assert fit.grid.prob(A, P) == pytest.approx(0.056, abs=0.002)
    assert fit.grid.prob(H, P) == pytest.approx(0.036, abs=0.002)


def test_a_whole_number_line_needs_two_books():
    game = _game()
    one = price_half_at_book(GridBook(TRUE_GRID, game), BOOK_GAME, game)
    prices = {c: WzParlayPrice(3.6, 500) for c in CORNERS}
    alone = HalfLine(game, "priced", book_fairs={"DraftKings": one}, wz_prices=prices)
    assert all(alone.fair_prob(c) is None for c in CORNERS)
    lower = HalfBookFair("BetMGM", {c: f - 0.01 for c, f in one.corner_fairs.items()}, 1.2,
                         push_returns=one.push_returns)
    both = HalfLine(game, "priced", book_fairs={"DraftKings": one, "BetMGM": lower}, wz_prices=prices)
    assert both.fair_prob("away_over") == pytest.approx(lower.effective_fair("away_over", 3.6))
    half = _game(away_spread=-6.5, total=26.5)
    assert HalfLine(half, "priced", book_fairs={"DraftKings": HalfBookFair("DraftKings", {c: 0.25 for c in CORNERS}, 1.2)},
                    wz_prices=prices).fair_prob("away_over") == 0.25


def test_a_held_whole_number_parlay_pays_pushes_at_its_own_price(monkeypatch):
    from nfl_specials import halves_board
    from nfl_specials.halves_board import HalfBoard, HeldHalfParlay, size_halves
    game = _game()
    fair = price_half_at_book(GridBook(TRUE_GRID, game), BOOK_GAME, game)
    line = HalfLine(game, "priced", book_fairs={"DraftKings": fair, "BetMGM": fair},
                    wz_prices={c: WzParlayPrice(5.0, 500) for c in CORNERS})
    groups_seen = []
    real = halves_board.budgeted_stakes
    monkeypatch.setattr(halves_board, "budgeted_stakes",
                        lambda groups, *args: groups_seen.append(groups) or real(groups, *args))
    board = HalfBoard(datetime.now(timezone.utc), None, [line], {})
    size_halves(board, 5000, 0.25, None, held=(HeldHalfParlay(game.wz_game_id, "away_over", 50.0, 300),))
    held = [h for group in groups_seen[0] for h in group.held]
    # Placed at +300 (4.0): its pushes pay against 4.0, not today's 5.0.
    assert [h.fair_prob for h in held] == pytest.approx([0.20 + fair.push_returns["away_over"] / 4.0])
