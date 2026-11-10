"""How a ProphetX game and its markets read as board lines: shared by the
ProphetX bet source (sources/prophetx.py: which line an order bought) and the
Place desk (prophetx_place.py: which strike a board line is, and what it costs
now), so a Place bet reads back as the line it was placed on.

ProphetX shapes these rules rest on (docs.prophetx.co, 2026-10-31; the
Market Data guide and the v4 schemas; NOT a live capture):
  sport event  {event_id, name, scheduled, status, tournament_id,
               tournament_name, competitors [{id, name, display_name,
               abbreviation, side: "home" | "away"}]}
  market       {id, name, type moneyline | spread | total | sup_moneyline,
               status, strike, selections, market_strikes [market...]}: a
               spread or total's other numbers sit under market_strikes, each
               shaped like a market.
  selections   a list per side of the market (one per team, or Over and
               Under), each a ladder of offers {outcome_id, competitor_id,
               name, display_name, price (AMERICAN), strike, display_strike,
               strike_id, quantity, value}. The ladder's best price is the one
               a buyer takes (mlb_sgp/prophetx.py reads the top entry as the
               standalone price); `quantity` is the money resting there.
  periods      only in the market's name ("1st Half Spread").
  sub_type     live read 2026-11-10 (Kent State @ Bowling Green): the game
               total is "total"; team totals "team_total_points"; player
               props "player_total_receiving_yards" etc., which also carry a
               player_id and are typed "total" with names like "Winn Sharp
               Receiving Yards" (no word "player"). Spread and moneyline
               sub_types are unread, so props are dropped by what marks them
               (game_line_market), not by a list of game sub_types.

Unverified until the first live reads, and so failing closed: whether a spread
offer's `strike` is its own side's number (+3.5 on the dog) or the market's
one number; both sides must read as opposite numbers, else the market is
unreadable. Pure: no network, no clock.
"""
import math
import re

# ProphetX tournament (by id, then by name) -> the board's league path.
# Ids from the docs' get_sport_events example (31 NFL, 132 NBA, 234 NHL, 109 MLB).
LEAGUE_OF_TOURNAMENT_ID = {31: "nfl", 132: "nba", 234: "nhl", 109: "mlb"}
LEAGUE_OF_TOURNAMENT_NAME = {
    "nfl": "nfl", "ncaaf": "cfb", "ncaa football": "cfb", "college football": "cfb", "ncaa fbs": "cfb",
    "nba": "nba", "ncaab": "cbb", "ncaa basketball": "cbb", "college basketball": "cbb", "ncaa men's basketball": "cbb",
    "wnba": "wnba", "mlb": "mlb", "major league baseball": "mlb", "nhl": "nhl",
}
BET_TYPES = {"moneyline": "moneyline", "spread": "spread", "total": "total"}
HALF_RE = re.compile(r"\b(1st|first)\s+half\b|\b1h\b", re.IGNORECASE)
FIRST_FIVE_RE = re.compile(r"\b(1st|first)\s+(5|five)\b|\bf5\b", re.IGNORECASE)
# Periods the board has no line for: never read as a full-game line.
OTHER_PERIOD_RE = re.compile(
    r"\b(2nd|second|3rd|third|4th|fourth)\s+(half|quarter|period)\b|\b(1st|first)\s+(quarter|period|inning)\b"
    r"|\bquarter\b|\binning\b|\bovertime\b|\bteam total\b|\bplayer\b", re.IGNORECASE)
TEAM_SIDES = ("away", "home")
# sub_type words that mark a market as not the game's own line (a player prop, a team total).
NOT_GAME_LINE_SUB_TYPE_WORDS = ("player", "team")


def _text(value: object) -> str:
    return value.strip() if isinstance(value, str) else ""


def number_of(value: object) -> float | None:
    """A finite number off a ProphetX field (number or numeric text, "+150"
    included), else None. Shared by the source and the Place desk so both read
    a field the same way."""
    if isinstance(value, bool) or value is None or value == "":
        return None
    try:
        number = float(str(value).replace("+", "")) if isinstance(value, str) else float(value)
    except (TypeError, ValueError):
        return None
    return number if math.isfinite(number) else None


def league_of(event: dict) -> str | None:
    """The board's league path of a ProphetX game, or None (not covered)."""
    league = LEAGUE_OF_TOURNAMENT_ID.get(event.get("tournament_id"))
    return league or LEAGUE_OF_TOURNAMENT_NAME.get(_text(event.get("tournament_name")).lower())


def team_name_of(competitor: dict) -> str | None:
    return _text(competitor.get("name")) or _text(competitor.get("display_name")) or None


def teams_of(event: dict) -> dict | str:
    """{away: competitor, home: competitor} by each competitor's own `side`,
    or what was wrong."""
    competitors = event.get("competitors")
    if not isinstance(competitors, list) or len(competitors) != 2:
        return f"ProphetX game {event.get('event_id')} has {len(competitors) if isinstance(competitors, list) else 'no'} competitors"
    by_side = {_text(c.get("side")).lower(): c for c in competitors if isinstance(c, dict)}
    if set(by_side) != set(TEAM_SIDES):
        return f"ProphetX game {event.get('event_id')} names no home and away side ({sorted(by_side)})"
    if not team_name_of(by_side["away"]) or not team_name_of(by_side["home"]):
        return f"ProphetX game {event.get('event_id')} has a team with no name"
    return {"away": by_side["away"], "home": by_side["home"]}


def period_of(market: dict, league: str) -> str | None:
    """"FG" / "1H" / "F5" off the market's name, or None for a period the
    board has no line for."""
    name = " ".join(_text(market.get(key)) for key in ("name", "group_name", "category_name"))
    if FIRST_FIVE_RE.search(name):
        return "F5"
    if OTHER_PERIOD_RE.search(name):
        return None
    if HALF_RE.search(name):
        return "F5" if league == "mlb" else "1H"
    return "FG"


def game_line_market(market: dict) -> bool:
    """False for a player prop or a team total, which ProphetX types "total"
    like the game total and names without the word "player": a player_id, or
    a sub_type naming a player or a team, marks one."""
    if market.get("player_id") not in (None, "", 0):
        return False
    sub_type = _text(market.get("sub_type")).lower()
    return not any(word in sub_type for word in NOT_GAME_LINE_SUB_TYPE_WORDS)


def sub_markets_of(market: dict) -> list[dict]:
    """A market and every alt-number market nested under it."""
    nested = market.get("market_strikes")
    return [market] + [sub for sub in (nested if isinstance(nested, list) else []) if isinstance(sub, dict)]


def side_ladders_of(market: dict) -> list[list[dict]]:
    """The market's two sides' offer ladders (each a list of offers), or []
    when it does not have exactly two."""
    selections = market.get("selections")
    if not isinstance(selections, list):
        return []
    ladders = [[offer for offer in side if isinstance(offer, dict)] for side in selections if isinstance(side, list)]
    ladders = [ladder for ladder in ladders if ladder]
    return ladders if len(ladders) == 2 else []


def strike_id_of(ladder: list[dict]) -> str | None:
    """The one strike id a side's ladder buys, or None when it names none or two."""
    ids = {offer.get("strike_id") for offer in ladder if isinstance(offer.get("strike_id"), str) and offer.get("strike_id")}
    return ids.pop() if len(ids) == 1 else None


def _signed_points(offer: dict) -> float | None:
    """The offer's own number: display_strike when it carries a sign, else strike."""
    display = _text(offer.get("display_strike"))
    if display[:1] in ("+", "-"):
        return number_of(display)
    return number_of(offer.get("strike"))


def _team_slot(offer: dict, teams: dict) -> str | None:
    """'away' / 'home': the team an offer backs, by competitor id, then by name."""
    competitor_id = offer.get("competitor_id")
    if competitor_id is not None:
        for slot in TEAM_SIDES:
            if str(teams[slot].get("id")) == str(competitor_id):
                return slot
    name = (_text(offer.get("name")) or _text(offer.get("display_name"))).lower()
    hits = [slot for slot in TEAM_SIDES
            if any(name == candidate.lower() or name.startswith(f"{candidate.lower()} ")
                   for candidate in (_text(teams[slot].get(key)) for key in ("name", "display_name", "abbreviation"))
                   if candidate)]
    return hits[0] if len(hits) == 1 else None


def sides_of_market(market: dict, bet_type: str, teams: dict) -> list[dict] | str:
    """Both sides of one (sub-)market as board sides: [{side, points,
    strike_id, ladder}], or what made it unreadable."""
    ladders = side_ladders_of(market)
    if not ladders:
        return f"market {market.get('id')} has no two sides"
    sides = []
    for ladder in ladders:
        strike_id = strike_id_of(ladder)
        if strike_id is None:
            return f"a side of market {market.get('id')} names no single strike id"
        head = ladder[0]
        if bet_type == "total":
            word = (_text(head.get("name")) or _text(head.get("display_name"))).lower()
            side = "over" if word.startswith("over") else "under" if word.startswith("under") else None
            points = number_of(head.get("strike"))
            if side is None or points is None or points <= 0:
                return f"total market {market.get('id')} has a side that reads as neither Over nor Under at a number"
        else:
            side = _team_slot(head, teams)
            points = None if bet_type == "moneyline" else _signed_points(head)
            if side is None:
                return f"a side of market {market.get('id')} names neither team"
            if bet_type == "spread" and points is None:
                return f"a spread side of market {market.get('id')} has no number"
        sides.append({"side": side, "points": points, "strike_id": strike_id, "ladder": ladder})
    if sides[0]["side"] == sides[1]["side"]:
        return f"both sides of market {market.get('id')} read as {sides[0]['side']}"
    if bet_type == "spread" and sides[0]["points"] != -sides[1]["points"]:
        return (f"spread market {market.get('id')}'s sides read {sides[0]['points']} and {sides[1]['points']}, "
                "not one number each way")
    if bet_type == "total" and sides[0]["points"] != sides[1]["points"]:
        return f"total market {market.get('id')}'s sides read different numbers"
    return sides


def lines_of_event(event: dict, markets: list[dict]) -> list[dict] | str:
    """Every side of every game market of one game as a board line: [{league,
    betType, period, side, points, strike_id, otherStrikeId (the same pair's
    other side), marketId, ladder, awayTeam, homeTeam}], or what made the game unreadable. Markets whose period or
    sides do not read, and player props and team totals (game_line_market), are left out, never guessed."""
    league = league_of(event)
    if league is None:
        return f"league not supported ({event.get('tournament_name') or event.get('tournament_id')})"
    teams = teams_of(event)
    if isinstance(teams, str):
        return teams
    lines = []
    seen_strikes: set[str] = set()
    for market in markets:
        bet_type = BET_TYPES.get(_text(market.get("type")).lower())
        if bet_type is None or not game_line_market(market):
            continue
        period = period_of(market, league)
        if period is None:
            continue
        for sub in sub_markets_of(market):
            if not game_line_market(sub):
                continue
            sides = sides_of_market(sub, bet_type, teams)
            if isinstance(sides, str):
                continue
            for side, other in zip(sides, reversed(sides)):
                # A main number listed again under market_strikes is one line, not two.
                if side["strike_id"] in seen_strikes:
                    continue
                seen_strikes.add(side["strike_id"])
                lines.append({"league": league, "betType": bet_type, "period": period, "side": side["side"],
                              "points": side["points"], "strike_id": side["strike_id"], "otherStrikeId": other["strike_id"],
                              "marketId": sub.get("id"),
                              "marketName": _text(market.get("name")), "ladder": side["ladder"],
                              "awayTeam": team_name_of(teams["away"]), "homeTeam": team_name_of(teams["home"])})
    return lines


def best_offer_of(ladder: list[dict]) -> dict | None:
    """{price, depth}: the best (highest) American price a buyer takes off a
    side's ladder and the money resting at it, or None when nothing rests."""
    offers = [(offer.get("price"), number_of(offer.get("quantity"))) for offer in ladder]
    offers = [(price, depth) for price, depth in offers
              if isinstance(price, (int, float)) and not isinstance(price, bool) and abs(price) >= 100
              and depth is not None and depth > 0]
    if not offers:
        return None
    best = max(price for price, _depth in offers)
    return {"price": int(best), "depth": round(sum(depth for price, depth in offers if price == best), 2)}


def american_to_probability(american: float) -> float:
    """The implied probability of an American price (no vig removed)."""
    return 100 / (american + 100) if american > 0 else -american / (-american + 100)


def profit_at(stake: float, american: float) -> float:
    """What `stake` wins at an American price."""
    return stake * american / 100 if american > 0 else stake * 100 / -american
