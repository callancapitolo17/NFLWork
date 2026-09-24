"""Polymarket US bet source: the account's positions + activity history -> normalised bet records.

Polymarket US is the CFTC exchange (the app, api.polymarket.us), not the international
polymarket.com; Unabated lists the two as separate books, so the venue key is
`polymarket_us`.

Two halves, like sources/kalshi.py:
  normalize_polymarket_us(...)  pure: our own fills (from the activity history's trades),
                                open positions, settlements and one public event payload per
                                game -> one record per (market slug, contract side). Pinned by
                                tests/test_polymarket_us.py on
                                tests/fixtures/bets/polymarket_us_account.json.
  PolymarketUSSource            the network half (Source protocol): every poll signs two GETs
                                with the account's API key (positions; activities back to
                                HISTORY_DAYS, and further until every open or just-settled
                                position's fills are in) and makes one cached public GET per game.

Inputs:  POLYMARKET_US_KEY_ID / POLYMARKET_US_SECRET_KEY (created at polymarket.us/developer;
         config reads them from the environment, bets_service/.env, kalshi_draft/.env, then
         bet_logger/.env in the main checkout). Every request to api.polymarket.us carries
         X-PM-Access-Key, X-PM-Timestamp (ms) and X-PM-Signature = base64 Ed25519 signature of
         timestamp + method + path — the path WITHOUT its query string (a signed query string is
         refused 401, live 2026-09-23).
         GET api.polymarket.us/v1/portfolio/positions   {positions: {slug: position}, nextCursor, eof}
         GET api.polymarket.us/v1/portfolio/activities  {activities: [...], nextCursor, eof},
             newest first
         GET gateway.polymarket.us/v1/events?slug=<eventSlug>&limit=1&sportsMarketTypes=<type>
             public, no key — the event (teams, startTime) with only the markets of the traded
             type (~10 KB; the unfiltered event is up to 4.7 MB)
Outputs: list of records. Side effects: none on disk; the event cache lives in memory. Raises on
         a failed positions or activities pull, a trade in a state the exchange has not
         documented, or an own order without its side, so the service records a failed run and
         keeps the previous records (never a partial list). A failed event lookup only makes
         that game's records unmatchable for one poll.

Wire grammar (live pull of 2026-09-23: 1 open position, 26 trades, 5 settlements):
  One instrument per market: the YES contract. `price` on a trade is always the YES price.
  Buying NO is selling YES; the own order says which: `outcomeSide` OUTCOME_SIDE_YES|NO,
  `action` ORDER_ACTION_BUY|SELL (intent BUY_LONG / SELL_LONG / BUY_SHORT / SELL_SHORT).
  A trade carries both orders; ours is `aggressor` when `isAggressor`, else `passive` (the
  other is the counterparty's and is never read). Quantities are fractional decimal strings.
  position.netPositionDecimal: positive = YES held, negative = NO held.
  positionResolution.side: POSITION_RESOLUTION_SIDE_LONG = YES won, _SHORT = NO won,
  _NEUTRAL = neither (a tie settles at $0.50 on a winner market; unobserved, so status
  "unknown"); beforePosition is the size that settled.
  A trade's `market` is the market object: sportsMarketType
  "<sport>_(team|game)_<period>_<winner|spread|total>", `line`, `gameStartTime`, and
  marketSides [{long, description, team {id, name, alias, safeName, league, ordering}}]:
    winner   both sides name a team; the YES (long) side is either one.
    spread   each side's description is its own signed number ("+2.50" / "-2.50"); `line`
             may carry either sign, and the YES side may be the favourite (MLB F5 "-1.50").
    total    sides "Over" (YES) / "Under" with no team. A total whose sides name a team is a
             TEAM total ("GB over 6.5 1H", football_team_first_half_total) — unmatchable, as is
             football_team_points_full_game_total.
  Combos (slugs "caoc-…") are parlays: comboLegDetails on the trade, no `market`. Listed as one
  unmatchable record with its legs in raw (plan: combos fail closed).
  Event teams come in away, home order (checked against Kalshi's MINTB / DALSEA / TORBAL
  tickers and the `ordering` field on five leagues); a side's `ordering` that disagrees fails
  the record closed. Team spelling per league: `name` ("Minnesota Vikings", WNBA "Dallas"),
  but `safeName` for college ("Liberty", where `name` is the nickname "Flames") and the NHL
  ("Ottawa Senators", where `name` is "Senators").

Record conventions (plan § Kalshi specifics, same shape): one record per (slug, contract side)
seen in our fills; the positions endpoint is the truth for the open size; a settlement decides
won/lost; a side sold back to zero is closed; fills still holding contracts with neither a
position nor a settlement read open (the Kalshi rule: a false open undersizes the next bet, a
false settle would oversize it). An open position whose fills are nowhere in the account's whole
history is still listed, unmatchable and unpriced. Price is the VWAP of our buys on that side, in
the side's own dollars (NO = 1 - YES price), fees excluded; stake = contracts x price. awayKey /
homeKey stay None — the extension resolves them through teams.js.
"""
import base64
import logging
import re
import time
from collections.abc import Callable
from datetime import datetime, timezone

import requests
from cryptography.hazmat.primitives.asymmetric import ed25519

from unabated_ticket.bets_service import config
from unabated_ticket.bets_service.normalize import (
    EASTERN, cents_to_american, json_clean, round_cents, to_number, utc_now_iso)

log = logging.getLogger(__name__)

SOURCE = "polymarket_us_api"
VENUE = "polymarket_us"
API_BASE_URL = "https://api.polymarket.us"
GATEWAY_BASE_URL = "https://gateway.polymarket.us"
POSITIONS_PATH = "/v1/portfolio/positions"
ACTIVITIES_PATH = "/v1/portfolio/activities"
EVENTS_PATH = "/v1/events"
PAGE_LIMIT = 100
# A cursor that never reaches eof would otherwise loop forever; 200 pages of 100 is far past
# any history window this service reads.
MAX_PAGES = 200
HTTP_TIMEOUT_SEC = 20
# The secret key is base64 of the 32-byte Ed25519 seed followed by the public key.
ED25519_SEED_BYTES = 32

ACTIVITY_TRADE = "ACTIVITY_TYPE_TRADE"
ACTIVITY_RESOLUTION = "ACTIVITY_TYPE_POSITION_RESOLUTION"
# Every trade state the docs list: counted while the trade stands, skipped once the
# exchange voided it (busted) or the clearinghouse refused it (rejected).
COUNTED_TRADE_STATES = {
    "TRADE_STATE_NEW", "TRADE_STATE_CLEARED", "TRADE_STATE_INFLIGHT", "TRADE_STATE_PENDING_RISK",
    "TRADE_STATE_PENDING_CLEARED", "TRADE_STATE_CLEARING_ACKNOWLEDGED", "TRADE_STATE_RETRY_REQUEST",
}
REVERSED_TRADE_STATES = {"TRADE_STATE_BUSTED", "TRADE_STATE_REJECTED"}
CONTRACT_SIDES = {"OUTCOME_SIDE_YES": "yes", "OUTCOME_SIDE_NO": "no"}
ORDER_ACTIONS = {"ORDER_ACTION_BUY": "buy", "ORDER_ACTION_SELL": "sell"}
# The contract side that won; NEUTRAL (neither) is absent on purpose -> status "unknown".
RESOLUTION_WINNERS = {"POSITION_RESOLUTION_SIDE_LONG": "yes", "POSITION_RESOLUTION_SIDE_SHORT": "no"}

# Polymarket US league slug (team.league) -> the extension's league path (feed.LEAGUES).
LEAGUES = {"nfl": "nfl", "cfb": "cfb", "nba": "nba", "cbb": "cbb", "wnba": "wnba", "mlb": "mlb", "nhl": "nhl"}
# The team field that spells the team the way Unabated does (module docstring); "name" otherwise.
# cbb follows cfb (college) — no cbb event was listed on 2026-09-23 to check it against.
TEAM_NAME_FIELDS = {"cfb": "safeName", "cbb": "safeName", "nhl": "safeName"}
DEFAULT_TEAM_NAME_FIELD = "name"
SPORTS = {"football", "basketball", "baseball", "hockey"}
PERIODS = {"full_game": "FG", "first_half": "1H", "second_half": "2H", "first_quarter": "1Q",
           "second_quarter": "2Q", "third_quarter": "3Q", "fourth_quarter": "4Q", "first_five": "F5"}
BET_TYPES = {"winner": "moneyline", "spread": "spread", "total": "total"}
MARKET_TYPE_RE = re.compile(r"^(?P<sport>[a-z]+)_(?:team|game)_(?P<period>[a-z_]+?)_(?P<kind>winner|spread|total)$")
TOTAL_SIDES = {"over": "over", "under": "under"}

REASON_COMBO = "combo (parlay): its legs are not matched"
REASON_NO_FILLS = "open position with no fill in the account's activity history"


# ---- pure normaliser --------------------------------------------------------------

def parse_iso(value: object) -> datetime | None:
    """"2026-09-15T20:22:51.500805479Z" (nanoseconds) -> aware UTC datetime, or None."""
    if not isinstance(value, str) or not value:
        return None
    match = re.match(r"^(\d{4}-\d{2}-\d{2}T\d{2}:\d{2}:\d{2})(?:\.(\d+))?(Z|[+-]\d{2}:\d{2})?$", value)
    if not match:
        return None
    base, fraction, zone = match.groups()
    micros = (fraction or "0")[:6].ljust(6, "0")
    offset = "+00:00" if zone in (None, "Z") else zone
    return datetime.fromisoformat(f"{base}.{micros}{offset}").astimezone(timezone.utc)


def iso_utc(moment: datetime | None) -> str | None:
    return None if moment is None else moment.strftime("%Y-%m-%dT%H:%M:%SZ")


def activity_time(activity: dict) -> datetime | None:
    """When an activity happened: a trade's createTime, a settlement's updateTime, a balance
    change's createTime (the order the endpoint sorts on, newest first)."""
    for key in ("trade", "positionResolution", "accountBalanceChange"):
        body = activity.get(key)
        if isinstance(body, dict):
            return parse_iso(body.get("createTime") or body.get("updateTime"))
    return None


def own_order_of(trade: dict) -> dict:
    """Our order on a trade (the other one is the counterparty's)."""
    order = trade.get("aggressor") if trade.get("isAggressor") else trade.get("passive")
    if not isinstance(order, dict) or order.get("outcomeSide") not in CONTRACT_SIDES \
            or order.get("action") not in ORDER_ACTIONS:
        keys = sorted(order.keys()) if isinstance(order, dict) else type(order).__name__
        raise RuntimeError(f"Polymarket US trade {trade.get('id')} on {trade.get('marketSlug')}: own order "
                           f"without outcomeSide/action (keys {keys}); capture it and extend "
                           "sources/polymarket_us.py")
    return order


def fill_of(trade: dict) -> dict | None:
    """Our side of one trade, or None when the exchange reversed it."""
    state = trade.get("state")
    if state in REVERSED_TRADE_STATES:
        return None
    if state not in COUNTED_TRADE_STATES:
        raise RuntimeError(f"Polymarket US trade {trade.get('id')} on {trade.get('marketSlug')} is in state "
                           f"{state!r}, which the docs do not list; decide whether it stands before counting it")
    order = own_order_of(trade)
    contract_side = CONTRACT_SIDES[order["outcomeSide"]]
    yes_price = to_number((trade.get("price") or {}).get("value"))
    metadata = order.get("marketMetadata") or {}
    return {
        "tradeId": trade.get("id"),
        "slug": trade.get("marketSlug"),
        "eventSlug": metadata.get("eventSlug") or None,
        "contractSide": contract_side,
        "action": ORDER_ACTIONS[order["action"]],
        "contracts": to_number(trade.get("qtyDecimal")),
        "sidePrice": yes_price if contract_side == "yes" else 1 - yes_price,
        "createTime": trade.get("createTime"),
        "market": trade.get("market") if isinstance(trade.get("market"), dict) else None,
        "comboLegs": trade.get("comboLegDetails") or [],
    }


def aggregate_fills(fills: list[dict]) -> dict:
    """VWAP in side dollars over buys (sells only when there were no buys), net contracts and
    the first / last fill time, for one (slug, contract side); all None / 0 with no fills."""
    if not fills:
        return {"vwapDollars": None, "netContracts": 0.0, "firstFillAt": None, "lastFillAt": None, "fillCount": 0}
    ordered = sorted(fills, key=lambda fill: parse_iso(fill["createTime"]) or datetime.min.replace(tzinfo=timezone.utc))
    bought = bought_cost = sold = sold_cost = 0.0
    for fill in ordered:
        if fill["action"] == "sell":
            sold += fill["contracts"]
            sold_cost += fill["contracts"] * fill["sidePrice"]
        else:
            bought += fill["contracts"]
            bought_cost += fill["contracts"] * fill["sidePrice"]
    if bought > 0:
        vwap = bought_cost / bought
    elif sold > 0:
        vwap = sold_cost / sold
    else:
        vwap = None
    return {"vwapDollars": vwap, "netContracts": bought - sold, "firstFillAt": ordered[0]["createTime"],
            "lastFillAt": ordered[-1]["createTime"], "fillCount": len(ordered)}


def held_side_of(net_position: float) -> str | None:
    if net_position > 0:
        return "yes"
    if net_position < 0:
        return "no"
    return None


def net_of(position: dict | None) -> float:
    return to_number(position.get("netPositionDecimal")) if isinstance(position, dict) else 0.0


def status_of(contract_side: str, position: dict | None, resolution: dict | None,
              net_contracts: float) -> str:
    """A settlement decides first, then the open position; with neither, a side sold back to
    zero is closed and one still holding contracts is open (module docstring)."""
    if resolution:
        if held_side_of(net_of(resolution.get("beforePosition"))) != contract_side:
            return "closed"
        winner = RESOLUTION_WINNERS.get(resolution.get("side"))
        if winner is None:
            return "unknown"
        return "won" if winner == contract_side else "lost"
    open_side = held_side_of(net_of(position))
    if open_side is not None:
        return "open" if open_side == contract_side else "closed"
    return "closed" if net_contracts <= 0 else "open"


def contracts_of(status: str, position: dict | None, resolution: dict | None, net_contracts: float) -> float:
    if status == "open" and held_side_of(net_of(position)) is not None:
        return abs(net_of(position))
    if resolution and status != "closed":
        return abs(net_of(resolution.get("beforePosition")))
    return max(net_contracts, 0.0)


def closed_at_of(status: str, resolution: dict | None, last_fill_at: str | None) -> str | None:
    if status == "open":
        return None
    if resolution and status != "closed":
        return iso_utc(parse_iso(resolution.get("updateTime")))
    return iso_utc(parse_iso(last_fill_at))


def parse_market_type(market_type: object) -> dict | str:
    """"baseball_team_first_five_spread" -> {betType, period, sport}, or the unmatchable reason."""
    match = MARKET_TYPE_RE.match(market_type) if isinstance(market_type, str) else None
    if not match or match["sport"] not in SPORTS or match["period"] not in PERIODS:
        return f"market type not supported ({market_type or 'none'})"
    if match["period"] == "first_five" and match["sport"] != "baseball":
        return f"market type not supported ({market_type})"
    return {"betType": BET_TYPES[match["kind"]], "period": PERIODS[match["period"]], "sport": match["sport"]}


def team_name_of(team: dict, league: str) -> str | None:
    name = team.get(TEAM_NAME_FIELDS.get(league, DEFAULT_TEAM_NAME_FIELD))
    return name if isinstance(name, str) and name.strip() else None


def event_teams_of(event: dict) -> dict | str:
    """{league, away, home} from the event's two teams (away first), or the unmatchable reason."""
    teams = event.get("teams")
    if not isinstance(teams, list) or len(teams) != 2 or not all(isinstance(team, dict) for team in teams):
        return f"unreadable Polymarket US event ({event.get('slug')}: expected 2 teams)"
    venue_leagues = {team.get("league") for team in teams}
    if len(venue_leagues) != 1:
        return f"unreadable Polymarket US event ({event.get('slug')}: teams in leagues {sorted(map(str, venue_leagues))})"
    venue_league = venue_leagues.pop()
    league = LEAGUES.get(venue_league)
    if league is None:
        return f"league not supported ({venue_league})"
    away, home = teams
    names = (team_name_of(away, league), team_name_of(home, league))
    if None in names:
        return f"unreadable Polymarket US event ({event.get('slug')}: a team without a name)"
    return {"league": league, "awayId": away.get("id"), "homeId": home.get("id"),
            "awayTeam": names[0], "homeTeam": names[1]}


def market_sides_of(market: dict, contract_side: str) -> tuple[dict, dict] | str:
    """(held side, other side) of a two-sided market; YES is the `long` side."""
    sides = market.get("marketSides")
    if not isinstance(sides, list) or len(sides) != 2 or {bool(side.get("long")) for side in sides} != {True, False}:
        return f"unreadable Polymarket US market ({market.get('slug')}: expected one YES and one NO side)"
    yes_side, no_side = sorted(sides, key=lambda side: not side.get("long"))
    return (yes_side, no_side) if contract_side == "yes" else (no_side, yes_side)


def team_slot_of(side: dict, teams: dict) -> str | None:
    """"away" / "home" for a side whose team is one of the event's two, else None."""
    team = side.get("team") if isinstance(side.get("team"), dict) else {}
    team_id = team.get("id", side.get("teamId"))
    if team_id is not None and team_id == teams["awayId"]:
        return "away"
    if team_id is not None and team_id == teams["homeId"]:
        return "home"
    return None


def ordering_disagrees(sides: tuple[dict, dict], teams: dict) -> bool:
    """A side whose team says away/home against the event's list order."""
    for side in sides:
        ordering = (side.get("team") or {}).get("ordering") if isinstance(side.get("team"), dict) else None
        if ordering in ("away", "home") and team_slot_of(side, teams) not in (None, ordering):
            return True
    return False


def parse_signed_number(text: object) -> float | None:
    try:
        return float(str(text).strip())
    except (TypeError, ValueError):
        return None


def side_and_points_of(market: dict, spec: dict, contract_side: str, teams: dict) -> dict | str:
    """{side, points} in the held side's own terms, or the unmatchable reason."""
    slug = market.get("slug")
    sides = market_sides_of(market, contract_side)
    if isinstance(sides, str):
        return sides
    held, other = sides
    if spec["betType"] == "total":
        if any(side.get("team") or side.get("teamId") for side in sides):
            return f"team total ({market.get('sportsMarketType')})"
        side = TOTAL_SIDES.get(str(held.get("description") or "").strip().lower())
        line = parse_signed_number(market.get("line"))
        if side is None or line is None or line <= 0:
            return f"unreadable Polymarket US market ({slug}: total side {held.get('description')!r}, line {market.get('line')!r})"
        return {"side": side, "points": line}
    if ordering_disagrees(sides, teams):
        return f"unreadable Polymarket US market ({slug}: a side's away/home disagrees with the event)"
    held_slot, other_slot = team_slot_of(held, teams), team_slot_of(other, teams)
    if held_slot is None or other_slot is None or held_slot == other_slot:
        return f"unreadable Polymarket US market ({slug}: sides are not the event's two teams)"
    if spec["betType"] == "moneyline":
        return {"side": held_slot, "points": None}
    points, other_points = parse_signed_number(held.get("description")), parse_signed_number(other.get("description"))
    line = parse_signed_number(market.get("line"))
    if points is None or other_points is None or line is None or points != -other_points or abs(points) != abs(line):
        return (f"unreadable Polymarket US market ({slug}: spread sides {held.get('description')!r} / "
                f"{other.get('description')!r} against line {market.get('line')!r})")
    return {"side": held_slot, "points": points}


def game_fields_of(market: dict | None, event: dict | None, contract_side: str, event_slug: str | None) -> dict | str:
    """The league / teams / type / side fields of a single-market record, or the reason it
    fails closed."""
    if market is None:
        return "unreadable Polymarket US market (no market payload on its trades)"
    spec = parse_market_type(market.get("sportsMarketType"))
    if isinstance(spec, str):
        return spec
    if event is None:
        return f"unreadable Polymarket US event (no event payload for {event_slug})"
    teams = event_teams_of(event)
    if isinstance(teams, str):
        return teams
    taken = side_and_points_of(market, spec, contract_side, teams)
    if isinstance(taken, str):
        return taken
    event_start = parse_iso(event.get("startTime")) or parse_iso(market.get("gameStartTime"))
    return {
        "league": teams["league"],
        "eventStart": iso_utc(event_start),
        "eventDate": None if event_start is None else event_start.astimezone(EASTERN).strftime("%Y-%m-%d"),
        "awayTeam": teams["awayTeam"], "homeTeam": teams["homeTeam"],
        "betType": spec["betType"], "period": spec["period"],
        "side": taken["side"], "points": taken["points"],
    }


def _unmatchable(record: dict, reason: str) -> dict:
    record["unmatchable"] = reason
    return record


def combo_legs_of(legs: list[dict]) -> list[dict]:
    """The legs of a combo, trimmed to what says what was bet (raw only, never matched)."""
    return [{"slug": leg.get("slug"), "eventSlug": leg.get("eventSlug"), "outcome": leg.get("outcome"),
             "outcomeSide": leg.get("outcomeSide"), "eventStartTime": leg.get("eventStartTime"),
             "state": leg.get("state")} for leg in legs if isinstance(leg, dict)]


def latest_resolutions(activities: list[dict]) -> dict[str, dict]:
    """slug -> its latest settlement."""
    latest: dict[str, dict] = {}
    for activity in activities:
        resolution = activity.get("positionResolution")
        if activity.get("type") != ACTIVITY_RESOLUTION or not isinstance(resolution, dict):
            continue
        slug = resolution.get("marketSlug")
        held = latest.get(slug)
        if held is None or (parse_iso(resolution.get("updateTime")) or datetime.min.replace(tzinfo=timezone.utc)) \
                > (parse_iso(held.get("updateTime")) or datetime.min.replace(tzinfo=timezone.utc)):
            latest[slug] = resolution
    return latest


def fills_of(activities: list[dict]) -> list[dict]:
    fills = []
    for activity in activities:
        if activity.get("type") != ACTIVITY_TRADE or not isinstance(activity.get("trade"), dict):
            continue
        fill = fill_of(activity["trade"])
        if fill is not None:
            fills.append(fill)
    return fills


def normalize_polymarket_us(positions: dict[str, dict], activities: list[dict], events: dict[str, dict],
                            fetched_at: str | None) -> list[dict]:
    """One record per (slug, contract side) seen in our fills. `positions` is keyed by market
    slug, `events` by event slug; a game missing its event fails closed with a reason. Pure."""
    resolutions = latest_resolutions(activities)
    groups: dict[tuple[str, str], list[dict]] = {}
    for fill in fills_of(activities):
        groups.setdefault((fill["slug"], fill["contractSide"]), []).append(fill)
    for slug, position in positions.items():
        held = held_side_of(net_of(position))
        if held is not None:
            groups.setdefault((slug, held), [])
    records = []
    for (slug, contract_side), group in groups.items():
        agg = aggregate_fills(group)
        position = positions.get(slug)
        resolution = resolutions.get(slug)
        status = status_of(contract_side, position, resolution, agg["netContracts"])
        contracts = contracts_of(status, position, resolution, agg["netContracts"])
        entry = agg["vwapDollars"]
        market = next((fill["market"] for fill in reversed(group) if fill["market"]), None)
        event_slug = next((fill["eventSlug"] for fill in group if fill["eventSlug"]), None) \
            or ((position or {}).get("marketMetadata") or {}).get("eventSlug") or None
        combo_legs = next((fill["comboLegs"] for fill in group if fill["comboLegs"]), None) \
            or (resolution or {}).get("comboLegDetails") or []
        record = {
            "id": f"{VENUE}:{slug}:{contract_side}",
            "source": SOURCE,
            "venue": VENUE,
            "league": None, "eventStart": None, "eventDate": None,
            "awayTeam": None, "homeTeam": None, "awayKey": None, "homeKey": None,
            "rotation": None, "betType": "other", "period": None, "side": None, "points": None,
            "price": None if entry is None else cents_to_american(entry * 100),
            "stake": None if entry is None else round_cents(contracts * entry),
            "toWin": None if entry is None else round_cents(contracts * (1 - entry)),
            "contracts": contracts,
            "placedAt": iso_utc(parse_iso(agg["firstFillAt"])),
            "status": status,
            "closedAt": closed_at_of(status, resolution, agg["lastFillAt"]),
            "isParlayLeg": False, "parlayId": None, "legIndex": None, "legCount": None,
            "approx": [],
            "unmatchable": None,
            "sourceFetchedAt": fetched_at,
            "venueIds": {"marketSlug": slug, "eventSlug": event_slug},
            "raw": {
                "slug": slug, "eventSlug": event_slug, "contractSide": contract_side,
                "fillCount": agg["fillCount"], "vwapDollars": entry, "netContracts": agg["netContracts"],
                "netPosition": net_of(position) if position else None,
                "resolutionSide": resolution.get("side") if resolution else None,
                "marketType": market.get("sportsMarketType") if market else None,
                "line": market.get("line") if market else None,
                "marketQuestion": market.get("question") if market else None,
                "positionTitle": ((position or {}).get("marketMetadata") or {}).get("title"),
                "positionOutcome": ((position or {}).get("marketMetadata") or {}).get("outcome"),
            },
        }
        if not group:
            records.append(json_clean(_unmatchable(record, REASON_NO_FILLS)))
            continue
        if combo_legs:
            record["raw"]["comboLegs"] = combo_legs_of(combo_legs)
            records.append(json_clean(_unmatchable(record, REASON_COMBO)))
            continue
        game = game_fields_of(market, events.get(event_slug) if event_slug else None, contract_side, event_slug)
        if isinstance(game, str):
            records.append(json_clean(_unmatchable(record, game)))
            continue
        record.update(game)
        records.append(json_clean(record))
    return records


# ---- network half -----------------------------------------------------------------

def new_session() -> requests.Session:
    return requests.Session()


def signer_of(secret_key: str) -> ed25519.Ed25519PrivateKey:
    try:
        seed = base64.b64decode(secret_key, validate=True)[:ED25519_SEED_BYTES]
        return ed25519.Ed25519PrivateKey.from_private_bytes(seed)
    except ValueError as error:
        raise RuntimeError("POLYMARKET_US_SECRET_KEY is not a base64 Ed25519 key "
                           f"(the Secret Key polymarket.us/developer showed once): {error}") from None


def json_or_raise(response: object, what: str) -> dict:
    if response.status_code != 200:
        raise RuntimeError(f"Polymarket US {what}: expected HTTP 200, got {response.status_code}: "
                           f"{response.text[:200]}")
    try:
        body = response.json()
    except ValueError:
        raise RuntimeError(f"Polymarket US {what}: HTTP 200 without a JSON body") from None
    if not isinstance(body, dict):
        raise RuntimeError(f"Polymarket US {what}: expected a JSON object, got {type(body).__name__}")
    return body


class PolymarketUSSource:
    """Source protocol implementation for the Polymarket US account (signed GETs only).

    `session_factory()` builds an HTTP session with .get() (requests by default; tests inject a
    fake). The key is read from config here and nowhere else; it never enters a record or a log.
    """

    name = VENUE

    def __init__(self, key_id: str | None = None, secret_key: str | None = None,
                 history_days: int | None = None, poll_sec: float | None = None,
                 session_factory: Callable[[], object] = new_session,
                 clock: Callable[[], float] = time.time):
        self.poll_sec = poll_sec if poll_sec is not None else config.POLYMARKET_US_POLL_SEC
        self._key_id = key_id if key_id is not None else config.POLYMARKET_US_KEY_ID
        self._secret_key = secret_key if secret_key is not None else config.POLYMARKET_US_SECRET_KEY
        self._history_days = history_days if history_days is not None else config.POLYMARKET_US_HISTORY_DAYS
        self._session = session_factory()
        self._clock = clock
        self._signer: ed25519.Ed25519PrivateKey | None = None
        self._events: dict[str, dict] = {}

    # -- signed account GETs ------------------------------------------------------------

    def _signed_headers(self, path: str) -> dict:
        if not self._key_id or not self._secret_key:
            raise RuntimeError("Polymarket US API key missing: set POLYMARKET_US_KEY_ID and "
                               "POLYMARKET_US_SECRET_KEY in bet_logger/.env")
        if self._signer is None:
            self._signer = signer_of(self._secret_key)
        timestamp = str(int(self._clock() * 1000))
        signature = self._signer.sign(f"{timestamp}GET{path}".encode())
        return {"X-PM-Access-Key": self._key_id, "X-PM-Timestamp": timestamp,
                "X-PM-Signature": base64.b64encode(signature).decode(), "Accept": "application/json"}

    def _get_signed_page(self, path: str, cursor: str | None) -> dict:
        params = {"limit": PAGE_LIMIT}
        if cursor:
            params["cursor"] = cursor
        response = self._session.get(f"{API_BASE_URL}{path}", params=params, headers=self._signed_headers(path),
                                     timeout=HTTP_TIMEOUT_SEC)
        return json_or_raise(response, f"GET {path}")

    def _fetch_positions(self) -> dict[str, dict]:
        positions: dict[str, dict] = {}
        cursor = None
        for _page in range(MAX_PAGES):
            body = self._get_signed_page(POSITIONS_PATH, cursor)
            page = body.get("positions")
            if not isinstance(page, dict):
                raise RuntimeError(f"Polymarket US positions: expected a slug -> position map, got "
                                   f"{type(page).__name__}")
            positions.update(page)
            cursor = body.get("nextCursor")
            if body.get("eof") or not cursor:
                return positions
        raise RuntimeError(f"Polymarket US positions: no eof after {MAX_PAGES} pages")

    def _fetch_activities(self, positions: dict[str, dict]) -> list[dict]:
        """Newest first. Paging stops once a page ends before the history_days cutoff AND every
        open position and every settlement read so far has a fill in hand: a bet placed weeks
        before its game keeps its price while open, and its settlement still finds its fills
        (without them the store's open row would never close)."""
        cutoff = datetime.fromtimestamp(self._clock() - self._history_days * 86400, tz=timezone.utc)
        needs_fill = {slug for slug, position in positions.items() if held_side_of(net_of(position)) is not None}
        has_fill: set[str] = set()
        activities: list[dict] = []
        cursor = None
        for _page in range(MAX_PAGES):
            body = self._get_signed_page(ACTIVITIES_PATH, cursor)
            page = body.get("activities")
            if not isinstance(page, list):
                raise RuntimeError(f"Polymarket US activities: expected a list, got {type(page).__name__}")
            activities.extend(page)
            for activity in page:
                if activity.get("type") == ACTIVITY_TRADE:
                    has_fill.add((activity.get("trade") or {}).get("marketSlug"))
                elif activity.get("type") == ACTIVITY_RESOLUTION:
                    needs_fill.add((activity.get("positionResolution") or {}).get("marketSlug"))
            cursor = body.get("nextCursor")
            oldest = activity_time(page[-1]) if page else None
            past_window = oldest is not None and oldest < cutoff
            if body.get("eof") or not cursor or (past_window and needs_fill <= has_fill):
                return activities
        raise RuntimeError(f"Polymarket US activities: no eof after {MAX_PAGES} pages")

    # -- public event lookups (cached) -----------------------------------------------------

    def _lookup_event(self, event_slug: str, market_type: str) -> dict | None:
        """One public GET, filtered to the traded market's type so the event comes back small
        (and at all: a filter that matches none of its markets returns no event). None (logged)
        on any failure, so one bad game fails closed instead of failing the poll."""
        params = {"slug": event_slug, "limit": 1, "sportsMarketTypes": market_type}
        try:
            response = self._session.get(f"{GATEWAY_BASE_URL}{EVENTS_PATH}", params=params,
                                         headers={"Accept": "application/json"}, timeout=HTTP_TIMEOUT_SEC)
            events = json_or_raise(response, f"event {event_slug}").get("events")
        except (RuntimeError, requests.RequestException) as error:
            log.warning("polymarket_us: event lookup %s failed: %s", event_slug, error)
            return None
        matches = [event for event in events or [] if isinstance(event, dict) and event.get("slug") == event_slug]
        if len(matches) != 1:
            log.warning("polymarket_us: event lookup %s returned %d matching events", event_slug, len(matches))
            return None
        return matches[0]

    def _ensure_events(self, fills: list[dict]) -> None:
        wanted = {}
        for fill in fills:
            market_type = (fill["market"] or {}).get("sportsMarketType")
            if fill["eventSlug"] and market_type and fill["eventSlug"] not in self._events:
                wanted.setdefault(fill["eventSlug"], market_type)
        for event_slug, market_type in sorted(wanted.items()):
            event = self._lookup_event(event_slug, market_type)
            if event is not None:
                self._events[event_slug] = event

    # -- Source protocol --------------------------------------------------------------

    def fetch(self) -> list[dict]:
        positions = self._fetch_positions()
        activities = self._fetch_activities(positions)
        self._ensure_events(fills_of(activities))
        records = normalize_polymarket_us(positions, activities, self._events, utc_now_iso())
        n_open = sum(1 for record in records if record["status"] == "open")
        n_unmatchable = sum(1 for record in records if record["unmatchable"])
        log.info("polymarket_us: %d positions + %d activities -> %d records (%d open, %d unmatchable)",
                 len(positions), len(activities), len(records), n_open, n_unmatchable)
        return records


def source_if_configured() -> PolymarketUSSource | None:
    """The source when the API key is configured (the BFA / Wagerzon pattern); a missing key
    is logged once with the fix."""
    if not (config.POLYMARKET_US_KEY_ID and config.POLYMARKET_US_SECRET_KEY):
        log.warning("polymarket_us: no POLYMARKET_US_KEY_ID / POLYMARKET_US_SECRET_KEY (bet_logger/.env in "
                    "the main checkout) — Polymarket US source not registered")
        return None
    return PolymarketUSSource()
