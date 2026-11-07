"""Scope rule (#283): an RFQ is in scope when every leg is a full-game
moneyline, spread or total in a league Unabated prices, and no two legs share
a game. No leg-count cap (Cal, 2026-11-07).

Kalshi tickers: series-EVENTCODE-STRIKE, e.g.
    KXNFLSPREAD-26NOV09KCDEN-KC7   series KXNFLSPREAD, game 26NOV09KCDEN
    KXMLBTOTAL-25JUN271905NYYBOS-9 (MLB codes carry the first-pitch time,
                                    so doubleheader games stay distinct)
Period series (KXNCAAF1HTOTAL, KXMLBF5...) and props (KXNBATEAMTOTAL) are not
in SERIES, so they fall out as unsupported_leg.
"""
from dataclasses import dataclass

# Kalshi series root -> league label. Order matters only for readability.
LEAGUES: dict[str, str] = {
    "KXNFL": "NFL",
    "KXNCAAF": "CFB",
    "KXNBA": "NBA",
    "KXNCAAMB": "CBB",
    "KXMLB": "MLB",
}
MARKETS: dict[str, str] = {"GAME": "ml", "SPREAD": "spread", "TOTAL": "total"}

# Exact full-game series names, e.g. "KXNFLSPREAD" -> ("NFL", "spread").
SERIES: dict[str, tuple[str, str]] = {
    root + suffix: (league, market)
    for root, league in LEAGUES.items()
    for suffix, market in MARKETS.items()
}
# Door prefixes for the WS source: anything that could be one of SERIES.
LEAGUE_PREFIXES: tuple[str, ...] = tuple(LEAGUES)
# Decline reasons evaluate() can return (no pricing happened).
SCOPE_REASONS = frozenset({"no_legs", "unsupported_leg", "bad_side",
                           "single_leg", "same_game"})


@dataclass(frozen=True)
class LegInfo:
    market_ticker: str
    side: str
    series: str
    league: str | None     # None when the leg is out of scope
    market: str | None     # ml | spread | total
    game_key: str | None   # league:EVENTCODE


@dataclass(frozen=True)
class ScopeVerdict:
    in_scope: bool
    reason: str | None          # None when in scope
    detail: str | None
    legs: list[LegInfo]
    leagues: list[str]          # distinct, in leg order
    n_games: int
    same_game: bool


def classify_leg(leg: dict) -> LegInfo:
    ticker = str(leg.get("market_ticker") or "")
    side = str(leg.get("side") or "").lower()
    parts = ticker.split("-")
    series = parts[0]
    hit = SERIES.get(series)
    if hit is None or len(parts) < 2 or not parts[1]:
        return LegInfo(ticker, side, series, None, None, None)
    league, market = hit
    return LegInfo(ticker, side, series, league, market, f"{league}:{parts[1]}")


def in_league(legs: list[dict]) -> bool:
    """Cheap door check: every leg's ticker starts with a priced league root.
    Mirrors the WS door so the REST fallback (unfiltered) is counted the
    same way instead of writing a row per foreign RFQ."""
    return all(str(l.get("market_ticker") or "").startswith(LEAGUE_PREFIXES)
               for l in legs)


def evaluate(legs: list[dict] | None) -> ScopeVerdict:
    if not legs:
        return ScopeVerdict(False, "no_legs", None, [], [], 0, False)
    infos = [classify_leg(l) for l in legs]
    leagues = list(dict.fromkeys(i.league for i in infos if i.league))
    games = [i.game_key for i in infos if i.game_key]
    n_games = len(set(games))
    same_game = len(games) != n_games

    def verdict(reason, detail=None):
        return ScopeVerdict(reason is None, reason, detail, infos, leagues,
                            n_games, same_game)

    bad = [i for i in infos if i.league is None]
    if bad:
        return verdict("unsupported_leg", ",".join(sorted({i.series for i in bad})))
    bad_side = [i for i in infos if i.side not in ("yes", "no")]
    if bad_side:
        return verdict("bad_side", bad_side[0].market_ticker)
    if len(infos) < 2:
        return verdict("single_leg")
    if same_game:
        dup = sorted({g for g in games if games.count(g) > 1})
        return verdict("same_game", ",".join(dup))
    return verdict(None)
