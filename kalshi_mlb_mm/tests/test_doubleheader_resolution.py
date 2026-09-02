"""Doubleheader identity across the two ways a game gets resolved.

Kalshi appends G1/G2 to both games of a doubleheader
(KXMLBGAME-26SEP041410DETCLEG1 / ...1915DETCLEG2, live 2026-09-01). Teaching
the suffix parser that grammar un-drops those games, which is the point — but
it also un-drops them for the OLDER resolvers that match on team names alone.
Those must keep declining: a doubleheader is exactly the case where
``WHERE home_team=? AND away_team=? LIMIT 1`` returns the wrong game, and #95
measured the resulting fairs off by 0.05-0.11.

So the invariant this file pins is a split one:
  * leg-surface identity (full suffix)  -> both games priceable, kept apart
  * team-name resolution                -> both games declined
"""
import pytest

from kalshi_common import legset
from kalshi_common.leg_types import game_number_from_suffix

G1 = "26SEP041410DETCLEG1"
G2 = "26SEP041915DETCLEG2"
PLAIN = "26AUG271910MILNYM"


def _total_leg(suffix: str):
    return legset.parse_leg({"event_ticker": f"KXMLBTOTAL-{suffix}",
                             "market_ticker": f"KXMLBTOTAL-{suffix}-9",
                             "side": "yes"})


SINGLE = "26SEP051410DETCLE"     # same pairing, no doubleheader


@pytest.fixture
def target_lines_db(tmp_path, monkeypatch):
    """A MARKET_DB holding exactly one DET @ CLE row.

    The positive control matters more than the negative one: without it the
    doubleheader assertions would pass on an absent database and prove
    nothing.
    """
    import duckdb

    from kalshi_mlb_mm import config, main
    path = tmp_path / "market.duckdb"
    con = duckdb.connect(str(path))
    con.execute("CREATE TABLE mlb_target_lines "
                "(game_id VARCHAR, home_team VARCHAR, away_team VARCHAR)")
    con.execute("INSERT INTO mlb_target_lines VALUES "
                "('odds-api-id', 'Cleveland Guardians', 'Detroit Tigers')")
    con.close()
    monkeypatch.setattr(config, "MARKET_DB", path)
    return main


class TestTeamNameResolutionFailsClosed:
    def test_an_ordinary_game_still_resolves(self, target_lines_db):
        main = target_lines_db
        assert (main._resolve_game_for_legs_uncached([_total_leg(SINGLE)])
                == "odds-api-id")

    def test_a_doubleheader_leg_resolves_to_no_game(self, target_lines_db):
        # Both games match that same single row on teams alone, so LIMIT 1
        # would hand each of them the other's line. Decline instead.
        main = target_lines_db
        assert main._resolve_game_for_legs_uncached([_total_leg(G1)]) is None
        assert main._resolve_game_for_legs_uncached([_total_leg(G2)]) is None

    def test_the_guard_is_the_game_number_not_a_parse_failure(self):
        # The codes DO parse now; the decline is a deliberate risk gate, so it
        # must key on the marker rather than on _parse_event_suffix failing.
        from kalshi_common.leg_types import _parse_event_suffix
        assert _parse_event_suffix(G1) == ("DET", "CLE")
        assert game_number_from_suffix(G1) == 1


class TestRfiDoubleheaderExclusion:
    """kalshi_rfi quotes one market per game off team-name book matching, so
    it drops doubleheaders outright. Its same-ET-day COUNT heuristic misses a
    pair whose first game already started and delisted; the suffix marker
    does not."""

    def _game(self, suffix, **kw):
        from kalshi_rfi.discovery import RfiGame, parse_suffix_start_utc
        return RfiGame(ticker=f"KXMLBRFI-{suffix}", suffix=suffix,
                       home_team="Cleveland Guardians",
                       away_team="Detroit Tigers",
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
