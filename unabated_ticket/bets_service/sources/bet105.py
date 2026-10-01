"""Bet105 bet source: the open bets the extension reads from the account's own Chrome
-> normalised bet records (2026-09-29).

Bet105 (bet105.ag) is a LinePros white-label behind Cloudflare, which challenges any
request that is not the browser session's own — so unlike every other venue the network
half is not here. The Unabated Ticket extension, running in Cal's logged-in Chrome, makes
the same `getHistory` call the site's My Plays page makes (state 0 = open) on both LinePros
feeds and POSTs the answers to this service's `POST /bet105.json` (extension/bet105.js and
panel.js; service.py). This module is the pure half plus the two service rules:

  normalize_bet105(feeds, fetched_at)  every bet group of both feeds -> records, one per
                                       leg of a parlay; pinned by tests/test_bet105.py on
                                       tests/fixtures/bets/bet105_history.json.
  validate_push(body)                  the POST body's shape, or what is wrong with it.
  closed_by_absence(stored, pushed_ids, now)
                                       Bet105 shows no settled list (My Plays only ever asks
                                       for state 0), so an open record a complete push no
                                       longer carries is marked `closed` with no result —
                                       never guessed won or lost.

Inputs:  {fetchedAt, feeds: {prematch: [betGroup], live: [betGroup]}} — the `betGroups` of
         POST /__bff/__partner-{feed}/betLobbyV2/logic/ {a: "getHistory", state: "0"} as
         captured from the site 2026-09-29 (three open NFL first-half totals; the live
         feed answered []). Or {error} when the extension could not read the account.
Outputs: records. Side effects: none (the service stores them).

Bet group (one ticket): betGroupId (the native id; the app keys `${feed}:${betGroupId}`),
  ticketNumber, acceptTime (epoch s, fractional), risk / toWin (USD), betType (0 =
  straight; the app: Object.keys(BET_TYPES)[betType] || "straight"), state (0 = pending,
  -1 = void: the app's BET_HISTORY_STATE_* constants; nothing else observed), result,
  isWin, finalOdds (decimal, the ticket's), componentBets[] (the legs).
Leg: team1 / team2 in the feed's event order — away, home (bet105_odds/scraper.py reads the
  same eventData as [away, home, start]); leagueName ("NFL"), leagueId, sportId; periodId
  (the coefficient feed's period keys: m = full game, h1 / h2 halves, f5 = first five
  innings, q1-q4); marketId, the LinePros wager type (3 Money Line, 5 Total, 6 Spread,
  7 / 8 Team Total 1 / 2, 1 = 1X2, the rest props — the site's eventsMetadata wagertypes);
  marketStyleId (0 on every observed leg; a styled market is not a game line); key (the
  line: the total, or the AWAY spread — the odds scraper reads the spread `r` as the away
  number and negates it for home); subKey (the side: "1" the first side — Over, or the
  away team — "2" the second: the odds scraper's ML dict keys and its [over, under]
  arrays); finalOdds (decimal); eventStartTime (epoch s); eventId; description ("" on every
  observed leg: the selection lives in the ids, not in text); state.
Verified on the capture: three first-half totals, side code "1". Unobserved, pinned only by
hand-written fixture rows: an Under, a spread's sign, a moneyline, a parlay, the live feed.
Everything else fails closed with the reason and the codes: a market other than 3 / 5 / 6,
a styled market, a side code other than 1 / 2, a period or league name the tables do not
know (a college league name is added once it is seen, never guessed).
"""
from datetime import datetime, timezone

from unabated_ticket.bets_service.normalize import EASTERN, js_round, json_clean, parse_iso_ms, round_cents

SOURCE = "bet105_extension"
VENUE = "bet105"
FEEDS = ("prematch", "live")

MARKET_MONEYLINE = "3"
MARKET_TOTAL = "5"
MARKET_SPREAD = "6"
GAME_MARKETS = {MARKET_MONEYLINE: "moneyline", MARKET_TOTAL: "total", MARKET_SPREAD: "spread"}
# The other wager types the site's metadata names, for the reason a leg fails closed.
OTHER_MARKET_NAMES = {"1": "1X2", "2": "Double Chance", "4": "Draw No Bet", "7": "Team Total 1",
                      "8": "Team Total 2", "30": "Outright Winner"}
TOTAL_SIDES = {"1": "over", "2": "under"}
TEAM_SIDES = {"1": "away", "2": "home"}
PERIODS = {"m": "FG", "h1": "1H", "h2": "2H", "f5": "F5", "q1": "1Q", "q2": "2Q", "q3": "3Q", "q4": "4Q"}
LEAGUES = {"NFL": "nfl", "NBA": "nba", "MLB": "mlb", "NHL": "nhl", "WNBA": "wnba"}
STATUS_BY_STATE = {0: "open", -1: "void"}
BET_TYPE_STRAIGHT = 0
PLAIN_MARKET_STYLE = 0
MAX_GROUPS_PER_PUSH = 2000
REASON_CLOSED_BY_ABSENCE = "left the open list (Bet105 shows no settled bets)"
REASON_NO_LEGS = "bet carries no legs"


# ---- pure parser --------------------------------------------------------------------

def american_from_decimal(decimal_odds: object) -> int | None:
    """1.9524 -> -105, 2.5 -> +150; None below evens-of-nothing (<= 1) or not a number."""
    if isinstance(decimal_odds, bool) or not isinstance(decimal_odds, (int, float)) or decimal_odds <= 1:
        return None
    if decimal_odds >= 2:
        return js_round((decimal_odds - 1) * 100)
    return -js_round(100 / (decimal_odds - 1))


def _number(value: object) -> float | None:
    if value is None or isinstance(value, bool):
        return None
    try:
        return float(value)
    except (TypeError, ValueError):
        return None


def _money(value: object) -> float | None:
    number = _number(value)
    return None if number is None else round_cents(number)


def moment_of_epoch(seconds: object) -> datetime | None:
    number = _number(seconds)
    if number is None or number <= 0:
        return None
    return datetime.fromtimestamp(number, timezone.utc)


def iso_utc(moment: datetime | None) -> str | None:
    return None if moment is None else moment.strftime("%Y-%m-%dT%H:%M:%SZ")


def eastern_date_of(moment: datetime) -> str:
    return moment.astimezone(EASTERN).strftime("%Y-%m-%d")


def native_id_of(group: dict) -> str:
    value = group.get("betGroupId")
    if value in (None, "", 0, "0"):
        raise RuntimeError(f"Bet105 bet group carries no betGroupId; keys seen: {sorted(group.keys())}")
    return str(value)


def status_of(state: object) -> str:
    number = _number(state)
    return STATUS_BY_STATE.get(int(number), "unknown") if number is not None and number.is_integer() else "unknown"


def _text(value: object) -> str:
    return value.strip() if isinstance(value, str) else ""


def parse_leg(leg: dict) -> dict | str:
    """One componentBet -> {league, period, betType, side, points, price, awayTeam,
    homeTeam, eventStart}, or the reason it does not read as a game line."""
    market = str(leg.get("marketId") if leg.get("marketId") is not None else "")
    style = _number(leg.get("marketStyleId")) or 0
    if style != PLAIN_MARKET_STYLE and str(int(style)) != market:
        return f"styled market (marketStyleId {int(style)}, market {market}) is not a game line"
    bet_type = GAME_MARKETS.get(market)
    if bet_type is None:
        return f"market {market or 'blank'} ({OTHER_MARKET_NAMES.get(market, 'unknown wager type')}) is not a game line"
    period = PERIODS.get(_text(leg.get("periodId")).lower())
    if period is None:
        return f"unknown period code ({leg.get('periodId')!r})"
    league_name = _text(leg.get("leagueName"))
    league = LEAGUES.get(league_name)
    if league is None:
        return f"league not supported (leagueName {league_name or 'blank'!r})"
    away, home = _text(leg.get("team1")), _text(leg.get("team2"))
    if not away or not home:
        return "leg names no teams"
    price = american_from_decimal(leg.get("finalOdds"))
    if price is None:
        return f"no price (finalOdds {leg.get('finalOdds')!r})"
    side_code = _text(leg.get("subKey"))
    parsed = {"league": league, "period": period, "betType": bet_type, "price": price,
              "awayTeam": away, "homeTeam": home, "eventStart": moment_of_epoch(leg.get("eventStartTime")),
              "points": None}
    if bet_type == "moneyline":
        side = TEAM_SIDES.get(side_code)
        if side is None:
            return f"unknown moneyline side code ({side_code or 'blank'})"
        parsed["side"] = side
        return parsed
    line = _number(leg.get("key"))
    if line is None:
        return f"{bet_type} carries no line (key {leg.get('key')!r})"
    if bet_type == "total":
        side = TOTAL_SIDES.get(side_code)
        if side is None:
            return f"unknown total side code ({side_code or 'blank'})"
        parsed.update(side=side, points=line)
        return parsed
    side = TEAM_SIDES.get(side_code)
    if side is None:
        return f"unknown spread side code ({side_code or 'blank'})"
    # The key is the away number; the home side gets its negation (the odds scraper's rule).
    parsed.update(side=side, points=line if side == "away" else -line)
    return parsed


def _base_record(group: dict, feed: str, native_id: str, fetched_at: str | None) -> dict:
    return {
        "id": f"{VENUE}:{feed}:{native_id}",
        "source": SOURCE,
        "venue": VENUE,
        "league": None, "eventStart": None, "eventDate": None,
        "awayTeam": None, "homeTeam": None, "awayKey": None, "homeKey": None,
        "rotation": None, "betType": "other", "period": None, "side": None, "points": None,
        "price": None,
        "stake": _money(group.get("risk")),
        "toWin": _money(group.get("toWin")),
        "contracts": None,
        "placedAt": iso_utc(moment_of_epoch(group.get("acceptTime"))),
        "status": status_of(group.get("state")),
        "closedAt": None,
        "isParlayLeg": False, "parlayId": None, "legIndex": None, "legCount": None,
        "approx": [],
        "unmatchable": None,
        "sourceFetchedAt": fetched_at,
        "raw": {
            "nativeId": native_id,
            "feed": feed,
            "betGroupId": group.get("betGroupId"),
            "ticketNumber": group.get("ticketNumber"),
            "betType": group.get("betType"),
            "state": group.get("state"),
            "result": group.get("result"),
            "isWin": group.get("isWin"),
            "finalOdds": group.get("finalOdds"),
            "risk": group.get("risk"),
            "toWin": group.get("toWin"),
            "acceptTime": group.get("acceptTime"),
        },
    }


def _leg_record(record: dict, leg: dict) -> dict:
    """The leg's own fields into raw (so an unmatched row shows its codes), then the
    parsed line — or the reason it failed."""
    record["raw"].update({
        "betId": leg.get("betId"), "description": leg.get("description"), "team1": leg.get("team1"),
        "team2": leg.get("team2"), "leagueName": leg.get("leagueName"), "leagueId": leg.get("leagueId"),
        "sportId": leg.get("sportId"), "periodId": leg.get("periodId"), "marketId": leg.get("marketId"),
        "marketStyleId": leg.get("marketStyleId"), "key": leg.get("key"), "subKey": leg.get("subKey"),
        "eventId": leg.get("eventId"), "eventStartTime": leg.get("eventStartTime"),
        "legFinalOdds": leg.get("finalOdds"), "legState": leg.get("state"),
    })
    parsed = parse_leg(leg)
    if isinstance(parsed, str):
        record["unmatchable"] = parsed
        return record
    event_start = parsed["eventStart"]
    record.update({
        "league": parsed["league"],
        "eventStart": iso_utc(event_start),
        "eventDate": None if event_start is None else eastern_date_of(event_start),
        "awayTeam": parsed["awayTeam"], "homeTeam": parsed["homeTeam"],
        "betType": parsed["betType"], "period": parsed["period"], "side": parsed["side"],
        "points": parsed["points"], "price": parsed["price"],
    })
    return record


def normalize_group(group: dict, feed: str, fetched_at: str | None) -> list[dict]:
    """One bet group -> one record, or one per leg of a parlay. Pure."""
    native_id = native_id_of(group)
    legs = group.get("componentBets") or []
    base = _base_record(group, feed, native_id, fetched_at)
    if not legs:
        base["unmatchable"] = REASON_NO_LEGS
        return [json_clean(base)]
    if len(legs) == 1 and _number(group.get("betType")) == BET_TYPE_STRAIGHT:
        return [json_clean(_leg_record(base, legs[0]))]
    parlay_price = american_from_decimal(group.get("finalOdds"))
    records = []
    for index, leg in enumerate(legs):
        record = _base_record(group, feed, native_id, fetched_at)
        record["id"] = f"{base['id']}:leg{index}"
        record["raw"]["parlayPrice"] = parlay_price
        record.update({"isParlayLeg": True, "parlayId": base["id"], "legIndex": index, "legCount": len(legs)})
        records.append(json_clean(_leg_record(record, leg)))
    return records


def _groups_of(feed_groups: object) -> list[dict]:
    """The site answers a list; the older LinePros app read a dict keyed by id."""
    if isinstance(feed_groups, dict):
        return list(feed_groups.values())
    return list(feed_groups or [])


def normalize_bet105(feeds: dict, fetched_at: str | None) -> list[dict]:
    """Every bet group of both feeds -> records. Raises only when a group has no
    betGroupId — the one thing the store cannot work without."""
    records: list[dict] = []
    for feed in FEEDS:
        for group in _groups_of(feeds.get(feed)):
            records.extend(normalize_group(group, feed, fetched_at))
    return records


# ---- the push -------------------------------------------------------------------------

def validate_push(body: object) -> dict | str:
    """{fetchedAt, feeds: {prematch: [...], live: [...]}} or {error} of a
    POST /bet105.json body, or what is wrong with it. Both feeds must be present:
    a push is the account's whole open list, or it is an error — never half."""
    if not isinstance(body, dict):
        return f"body must be an object, got {type(body).__name__}"
    if "error" in body:
        if not isinstance(body["error"], str) or not body["error"].strip():
            return "error must be a non-empty string"
        return {"error": body["error"].strip()}
    fetched_at = body.get("fetchedAt")
    if not isinstance(fetched_at, str) or parse_iso_ms(fetched_at) is None:
        return f"fetchedAt must be an ISO time, got {fetched_at!r}"
    feeds = body.get("feeds")
    if not isinstance(feeds, dict) or set(feeds) != set(FEEDS):
        return f"feeds must be an object with exactly {list(FEEDS)}, got {sorted(feeds) if isinstance(feeds, dict) else type(feeds).__name__}"
    total = 0
    clean: dict[str, list[dict]] = {}
    for feed in FEEDS:
        groups = feeds[feed]
        if not isinstance(groups, (list, dict)):
            return f"feeds.{feed} must be a list, got {type(groups).__name__}"
        groups = _groups_of(groups)
        for index, group in enumerate(groups):
            if not isinstance(group, dict):
                return f"feeds.{feed}[{index}] must be an object, got {type(group).__name__}"
            if group.get("betGroupId") in (None, "", 0, "0"):
                return f"feeds.{feed}[{index}] carries no betGroupId; keys seen: {sorted(group.keys())}"
        total += len(groups)
        clean[feed] = groups
    if total > MAX_GROUPS_PER_PUSH:
        return f"at most {MAX_GROUPS_PER_PUSH} bet groups per push, got {total}"
    return {"fetchedAt": fetched_at, "feeds": clean}


def closed_by_absence(stored_records: list[dict], pushed_ids: set[str], now_iso: str) -> list[dict]:
    """The stored open Bet105 records a complete push no longer lists, as closed
    copies (no result: the venue shows none). Pure."""
    closed = []
    for record in stored_records:
        if record.get("venue") != VENUE or record.get("status") != "open" or record.get("id") in pushed_ids:
            continue
        copy = json_clean({**record, "raw": {**record.get("raw", {}), "closedReason": REASON_CLOSED_BY_ABSENCE}})
        copy.update({"status": "closed", "closedAt": now_iso, "sourceFetchedAt": now_iso})
        closed.append(copy)
    return closed
