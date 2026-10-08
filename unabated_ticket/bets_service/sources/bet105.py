"""Bet105 bet source: the open and settled bets the extension reads from the account's
own Chrome -> normalised bet records (open 2026-09-29, settled 2026-10-07).

Bet105 (bet105.ag) is a LinePros white-label behind Cloudflare, which challenges any
request that is not the browser session's own — so unlike every other venue the network
half is not here. The Unabated Ticket extension, running in Cal's logged-in Chrome, makes
the site's own calls and POSTs the answers to this service's `POST /bet105.json`
(extension/bet105.js and panel.js; service.py): `getHistory` (state 0 = open) on both
LinePros feeds, and the My Bets page's `wagers/search` (every wager the account has, graded
or not). This module is the pure half plus the service rules:

  normalize_bet105(feeds, fetched_at)  every open bet group of both feeds -> records, one
                                       per leg of a parlay.
  normalize_settled(wagers, fetched_at, known_feeds)
                                       every GRADED wager of wagers/search -> records with
                                       the venue's result, plus the reasons for the ones it
                                       could not read. Pending wagers are left to the open
                                       read (the venue also keeps a few never-graded March
                                       2026 wagers Pending, which getHistory does not list).
  feeds_by_native_id(records)          the feed each known ticket's record sits under
                                       (stored, or open in this push) — normalize_settled's
                                       known_feeds.
  merge_settled(open_records, settled_records)
                                       a ticket the settled read carries is settled: its
                                       open records give way.
  validate_push(body)                  the POST body's shape, or what is wrong with it.
  closed_by_absence(stored, pushed_ids, now)
                                       the fallback only: an open record neither read
                                       carries any more is marked `closed` with no result —
                                       never guessed won or lost.
Pinned by tests/test_bet105.py on tests/fixtures/bets/bet105_history.json.

Inputs:  {fetchedAt, feeds: {prematch: [betGroup], live: [betGroup]}, settled: [wager]}, or
         {error} when the extension could not read the account.
Outputs: records. Side effects: none (the service stores them).

Open read — POST /__bff/__partner-{feed}/betLobbyV2/logic/ {a: "getHistory", state: "0"},
captured 2026-09-29 (three open NFL first-half totals; the live feed answered []).
Bet group (one ticket): betGroupId (the native id; the app keys `${feed}:${betGroupId}`),
  ticketNumber, acceptTime (epoch s, fractional), risk / toWin (USD), isFreePlay, betType
  (0 = straight; the app: Object.keys(BET_TYPES)[betType] || "straight"), state (0 = pending,
  -1 = void: the app's BET_HISTORY_STATE_* constants; nothing else observed), result,
  isWin, finalOdds (decimal, the ticket's), componentBets[] (the legs).
Leg: team1 / team2 in the feed's event order — away, home (bet105_odds/scraper.py reads the
  same eventData as [away, home, start]); leagueName ("NFL"), leagueId, sportId; periodId
  (m = full game, h1 / h2 halves, f5 = first five innings, s1 / s4 the 1st / 4th quarter —
  seen on the settled read; q1-q4 the coefficient feed's quarter keys, unseen on a bet);
  marketId, the LinePros wager type (3 Money Line, 5 Total, 6 Spread, 7 / 8 Team Total 1 /
  2, 1 = 1X2, the rest props — the site's eventsMetadata wagertypes); marketStyleId (0 on
  every observed leg; a styled market is not a game line); key (the line: the total, or the
  AWAY spread — the odds scraper reads the spread `r` as the away number and negates it for
  home); subKey (the side: "1" the first side — Over, or the away team — "2" the second);
  finalOdds (decimal); eventStartTime (epoch s); eventId; description ("" on every observed
  leg: the selection lives in the ids, not in text); state.

Settled read — POST /__bff/api/wagers/search {} with the same X-Broker-CSRF header (403
{code: "CSRF_FAILED"} without it), captured 2026-10-07: 230 wagers placed 2025-09-06 to
2026-09-28, all productCode "PreMatch", category "SPORTS". The extension sends each wager
cut to the fields read here (extension/bet105.js settledWagerOf).
Wager: ticketNumber (a string: getHistory's betGroupId — so the record id is the open
  record's; one sequence across products, 230 of 230 unique), wagerId (getHistory's
  ticketNumber; the number My Bets shows), productCode ("PreMatch" -> the prematch feed;
  the feed comes first from the ticket's own open record, under whichever feed it sits),
  wagerStatus ("Win" 117, "Loss" 104, "Push" 4,
  "Pending" 5), placeTime / gradeTime (ISO UTC; gradeTime is the settle time, set on every
  graded wager), risk, toWin, result (the money the bet made: toWin on a win — to the tenth
  of a cent — minus risk on a loss, 0 on a push, 0 on a lost free play), isFreePlay,
  isCashout, wagerDetails (the text My Bets shows), properties {odds (the ticket's
  decimal), fmtOdds, grades (one "W" / "L" / "P" per leg), teaserName, fixedParlayName,
  legs}.
Settled leg: the open leg's ids under other names — league (the leagueName), team1 / team2,
  periodId, marketId (a number), side (1 / 2, a number: subKey), figure (key: the total, or
  the AWAY spread — checked on six real spreads against the bet's own text, both sides),
  fmtOdds (the leg's American price, a signed string on every leg; `odds` is decimal on
  newer legs and the American string on older ones, so it is not read), startTime
  ("2026-10-04 17:00:00+00:00"), legId (getHistory's betId), eventId. No marketStyleId.
Status: the site's own My Bets code (its Te(): isCashout first, then PENDING, CANCEL and
  NO_ACTION, else by result) and the 225 graded wagers: Win -> won, Loss -> lost, Push ->
  push, CANCEL / NO_ACTION -> void, isCashout -> closed (sold before settlement: the store's
  "no known P&L" status). A free play stakes nothing (the BFA convention), so a lost one
  costs 0, as its result says. A won record's toWin is the venue's `result`, what it paid.
  With those, won -> toWin, lost -> -stake, push / void -> 0 reproduces `result` on every
  graded wager captured. The record also carries `pnl` = `result`, the venue's own P&L
  (the field the tracker counts first, as for Kalshi), so a cash-out books its result.

Verified on the captures: totals of both sides, both spread sides, moneylines of both sides,
1st / 4th quarters, a 3-leg parlay (a $0 free bet), a free play, pushes; NFL, MLB and the
college leagues. Unobserved, pinned only by hand-written fixture rows: a void, a cash-out,
an open parlay, the live feed, a productCode other than PreMatch (it settles under its open
record's feed; with no open record ever seen it is skipped with the reason). Everything else fails closed with the reason
and the codes: a market other than 3 / 5 / 6, a styled market, a side code other than
1 / 2, a period or league name the tables do not know (a league name is added once it is
seen on a bet, never guessed; soccer stays out).
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
PERIODS = {"m": "FG", "h1": "1H", "h2": "2H", "f5": "F5", "s1": "1Q", "s4": "4Q",
           "q1": "1Q", "q2": "2Q", "q3": "3Q", "q4": "4Q"}
# LinePros league names as bets carry them. "College Football - FCS", "College Basketball
# Extra" and "NBA Preseason" are the venue's own extra leagues for the same sport.
LEAGUES = {"NFL": "nfl", "NBA": "nba", "MLB": "mlb", "NHL": "nhl", "WNBA": "wnba",
           "College Football": "cfb", "College Football - FCS": "cfb",
           "College Basketball": "cbb", "College Basketball Extra": "cbb", "NBA Preseason": "nba"}
STATUS_BY_STATE = {0: "open", -1: "void"}
BET_TYPE_STRAIGHT = 0
PLAIN_MARKET_STYLE = 0
MAX_GROUPS_PER_PUSH = 2000
REASON_CLOSED_BY_ABSENCE = "left the open list and the settled list does not carry it"
REASON_NO_LEGS = "bet carries no legs"

# The settled read (wagers/search).
SETTLED_STATUS_BY_WAGER_STATUS = {"WIN": "won", "LOSS": "lost", "PUSH": "push",
                                  "CANCEL": "void", "NO_ACTION": "void"}
WAGER_STATUS_PENDING = "PENDING"
STATUS_CASHED_OUT = "closed"
FEED_BY_PRODUCT = {"PreMatch": "prematch"}


# ---- pure parser --------------------------------------------------------------------

def american_from_decimal(decimal_odds: object) -> int | None:
    """1.9524 -> -105, 2.5 -> +150; None below evens-of-nothing (<= 1) or not a number."""
    if isinstance(decimal_odds, bool) or not isinstance(decimal_odds, (int, float)) or decimal_odds <= 1:
        return None
    if decimal_odds >= 2:
        return js_round((decimal_odds - 1) * 100)
    return -js_round(100 / (decimal_odds - 1))


def american_from_text(text: object) -> int | None:
    """"-117" -> -117, "+192" -> 192 (the settled leg's fmtOdds); None for anything else."""
    if not isinstance(text, str):
        return None
    stripped = text.strip()
    digits = stripped[1:] if stripped[:1] in ("+", "-") else stripped
    if not digits.isdigit():
        return None
    value = int(stripped)
    return value if abs(value) >= 100 else None


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


def stake_of(risk: object, is_free_play: object) -> float | None:
    """A free play stakes nothing — the venue settles a lost one at 0 (the BFA convention)."""
    return 0 if is_free_play is True else _money(risk)


def moment_of_epoch(seconds: object) -> datetime | None:
    number = _number(seconds)
    if number is None or number <= 0:
        return None
    return datetime.fromtimestamp(number, timezone.utc)


def moment_of_iso(text: object) -> datetime | None:
    """"2026-10-06T01:56:19Z" or "2026-10-06 00:15:00+00:00" -> UTC; None when it is not
    an ISO time carrying its zone (a zone is never guessed)."""
    if not isinstance(text, str) or not text.strip():
        return None
    try:
        moment = datetime.fromisoformat(text.strip())
    except ValueError:
        return None
    if moment.tzinfo is None:
        return None
    return moment.astimezone(timezone.utc)


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


def _code(value: object) -> str:
    """A code the feeds send as a string ("5", "1") or a number (5, 1) -> its text."""
    if isinstance(value, bool) or value is None:
        return ""
    if isinstance(value, float) and value.is_integer():
        return str(int(value))
    return value.strip() if isinstance(value, str) else str(value)


def parse_game_line(*, market: str, period_code: object, league_name: str, away: str, home: str,
                    price: int | None, price_text: str, side_code: str, line: float | None,
                    line_text: str, event_start: datetime | None) -> dict | str:
    """One leg's selection, already read off its feed -> {league, period, betType, side,
    points, price, awayTeam, homeTeam, eventStart}, or the reason it is not a game line.
    `price_text` / `line_text` name the feed's own field in a reason."""
    bet_type = GAME_MARKETS.get(market)
    if bet_type is None:
        return f"market {market or 'blank'} ({OTHER_MARKET_NAMES.get(market, 'unknown wager type')}) is not a game line"
    period = PERIODS.get(_text(period_code).lower())
    if period is None:
        return f"unknown period code ({period_code!r})"
    league = LEAGUES.get(league_name)
    if league is None:
        return f"league not supported ({league_name or 'blank'!r})"
    if not away or not home:
        return "leg names no teams"
    if price is None:
        return f"no price ({price_text})"
    parsed = {"league": league, "period": period, "betType": bet_type, "price": price,
              "awayTeam": away, "homeTeam": home, "eventStart": event_start, "points": None}
    if bet_type == "moneyline":
        side = TEAM_SIDES.get(side_code)
        if side is None:
            return f"unknown moneyline side code ({side_code or 'blank'})"
        parsed["side"] = side
        return parsed
    if line is None:
        return f"{bet_type} carries no line ({line_text})"
    if bet_type == "total":
        side = TOTAL_SIDES.get(side_code)
        if side is None:
            return f"unknown total side code ({side_code or 'blank'})"
        parsed.update(side=side, points=line)
        return parsed
    side = TEAM_SIDES.get(side_code)
    if side is None:
        return f"unknown spread side code ({side_code or 'blank'})"
    # The line is the away number; the home side gets its negation (the odds scraper's rule).
    parsed.update(side=side, points=line if side == "away" else -line)
    return parsed


def parse_leg(leg: dict) -> dict | str:
    """One open componentBet (getHistory) -> parse_game_line's result."""
    market = _code(leg.get("marketId"))
    style = _number(leg.get("marketStyleId")) or 0
    if style != PLAIN_MARKET_STYLE and str(int(style)) != market:
        return f"styled market (marketStyleId {int(style)}, market {market}) is not a game line"
    return parse_game_line(
        market=market, period_code=leg.get("periodId"), league_name=_text(leg.get("leagueName")),
        away=_text(leg.get("team1")), home=_text(leg.get("team2")),
        price=american_from_decimal(leg.get("finalOdds")), price_text=f"finalOdds {leg.get('finalOdds')!r}",
        side_code=_code(leg.get("subKey")), line=_number(leg.get("key")), line_text=f"key {leg.get('key')!r}",
        event_start=moment_of_epoch(leg.get("eventStartTime")))


def parse_settled_leg(leg: dict) -> dict | str:
    """One settled leg (wagers/search) -> parse_game_line's result. This feed carries no
    marketStyleId, so the market id alone decides."""
    return parse_game_line(
        market=_code(leg.get("marketId")), period_code=leg.get("periodId"), league_name=_text(leg.get("league")),
        away=_text(leg.get("team1")), home=_text(leg.get("team2")),
        price=american_from_text(leg.get("fmtOdds")), price_text=f"fmtOdds {leg.get('fmtOdds')!r}",
        side_code=_code(leg.get("side")), line=_number(leg.get("figure")), line_text=f"figure {leg.get('figure')!r}",
        event_start=moment_of_iso(leg.get("startTime")))


def _empty_record(record_id: str, fetched_at: str | None) -> dict:
    return {
        "id": record_id,
        "source": SOURCE,
        "venue": VENUE,
        "league": None, "eventStart": None, "eventDate": None,
        "awayTeam": None, "homeTeam": None, "awayKey": None, "homeKey": None,
        "rotation": None, "betType": "other", "period": None, "side": None, "points": None,
        "price": None, "stake": None, "toWin": None, "contracts": None,
        "placedAt": None, "status": None, "closedAt": None,
        "isParlayLeg": False, "parlayId": None, "legIndex": None, "legCount": None,
        "approx": [],
        "unmatchable": None,
        "sourceFetchedAt": fetched_at,
        "raw": {},
    }


def _base_record(group: dict, feed: str, native_id: str, fetched_at: str | None) -> dict:
    record = _empty_record(f"{VENUE}:{feed}:{native_id}", fetched_at)
    record.update({
        "stake": stake_of(group.get("risk"), group.get("isFreePlay")),
        "toWin": _money(group.get("toWin")),
        "placedAt": iso_utc(moment_of_epoch(group.get("acceptTime"))),
        "status": status_of(group.get("state")),
        "raw": {
            "nativeId": native_id,
            "feed": feed,
            "betGroupId": group.get("betGroupId"),
            "ticketNumber": group.get("ticketNumber"),
            "betType": group.get("betType"),
            "state": group.get("state"),
            "result": group.get("result"),
            "isWin": group.get("isWin"),
            "isFreePlay": group.get("isFreePlay"),
            "finalOdds": group.get("finalOdds"),
            "risk": group.get("risk"),
            "toWin": group.get("toWin"),
            "acceptTime": group.get("acceptTime"),
        },
    })
    return record


def _apply_parsed_leg(record: dict, parsed: dict | str) -> dict:
    """The parsed line onto the record — or the reason it failed."""
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
    return _apply_parsed_leg(record, parse_leg(leg))


def _as_parlay_leg(record: dict, ticket_id: str, index: int, leg_count: int, parlay_price: int | None) -> dict:
    record["id"] = f"{ticket_id}:leg{index}"
    record["raw"]["parlayPrice"] = parlay_price
    record.update({"isParlayLeg": True, "parlayId": ticket_id, "legIndex": index, "legCount": leg_count})
    return record


def normalize_group(group: dict, feed: str, fetched_at: str | None) -> list[dict]:
    """One open bet group -> one record, or one per leg of a parlay. Pure."""
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
        record = _as_parlay_leg(_base_record(group, feed, native_id, fetched_at), base["id"], index, len(legs), parlay_price)
        records.append(json_clean(_leg_record(record, leg)))
    return records


def _groups_of(feed_groups: object) -> list[dict]:
    """The site answers a list; the older LinePros app read a dict keyed by id."""
    if isinstance(feed_groups, dict):
        return list(feed_groups.values())
    return list(feed_groups or [])


def normalize_bet105(feeds: dict, fetched_at: str | None) -> list[dict]:
    """Every open bet group of both feeds -> records. Raises only when a group has no
    betGroupId — the one thing the store cannot work without."""
    records: list[dict] = []
    for feed in FEEDS:
        for group in _groups_of(feeds.get(feed)):
            records.extend(normalize_group(group, feed, fetched_at))
    return records


# ---- the settled read -----------------------------------------------------------------

def settled_status_of(wager: dict) -> str | None:
    """won / lost / push / void / closed (cashed out), or None: pending, or a status the
    site's code and the captures do not name."""
    if wager.get("isCashout") is True:
        return STATUS_CASHED_OUT
    return SETTLED_STATUS_BY_WAGER_STATUS.get(_text(wager.get("wagerStatus")).upper())


def _settled_base_record(wager: dict, feed: str, native_id: str, status: str, fetched_at: str | None) -> dict:
    properties = wager.get("properties") or {}
    paid = wager.get("result") if status == "won" else wager.get("toWin")
    record = _empty_record(f"{VENUE}:{feed}:{native_id}", fetched_at)
    record.update({
        "stake": stake_of(wager.get("risk"), wager.get("isFreePlay")),
        "toWin": _money(paid),
        "pnl": _money(wager.get("result")),
        "placedAt": iso_utc(moment_of_iso(wager.get("placeTime"))),
        "status": status,
        "closedAt": iso_utc(moment_of_iso(wager.get("gradeTime"))),
        "raw": {
            "nativeId": native_id,
            "feed": feed,
            "read": "settled",
            "ticketNumber": wager.get("ticketNumber"),
            "wagerId": wager.get("wagerId"),
            "productCode": wager.get("productCode"),
            "wagerStatus": wager.get("wagerStatus"),
            "placeTime": wager.get("placeTime"),
            "gradeTime": wager.get("gradeTime"),
            "risk": wager.get("risk"),
            "toWin": wager.get("toWin"),
            "result": wager.get("result"),
            "isFreePlay": wager.get("isFreePlay"),
            "isCashout": wager.get("isCashout"),
            "wagerDetails": wager.get("wagerDetails"),
            "odds": properties.get("odds"),
            "grades": properties.get("grades"),
            "teaserName": properties.get("teaserName"),
            "fixedParlayName": properties.get("fixedParlayName"),
        },
    })
    return record


def _settled_leg_record(record: dict, leg: dict, grade: object) -> dict:
    record["raw"].update({
        "legId": leg.get("legId"), "eventId": leg.get("eventId"), "sportId": leg.get("sportId"),
        "leagueId": leg.get("leagueId"), "league": leg.get("league"), "team1": leg.get("team1"),
        "team2": leg.get("team2"), "periodId": leg.get("periodId"), "period": leg.get("period"),
        "marketId": leg.get("marketId"), "market": leg.get("market"), "side": leg.get("side"),
        "figure": leg.get("figure"), "fmtOdds": leg.get("fmtOdds"), "startTime": leg.get("startTime"),
        "legGrade": grade,
    })
    return _apply_parsed_leg(record, parse_settled_leg(leg))


def normalize_settled_wager(wager: dict, feed: str, status: str, fetched_at: str | None) -> list[dict]:
    """One graded wager -> one record, or one per leg of a parlay, every leg on the
    ticket's stake, payout and status. Pure."""
    native_id = str(wager["ticketNumber"]).strip()
    properties = wager.get("properties") or {}
    legs = properties.get("legs") or []
    grades = properties.get("grades") if isinstance(properties.get("grades"), list) else []
    base = _settled_base_record(wager, feed, native_id, status, fetched_at)
    if not legs:
        base["unmatchable"] = REASON_NO_LEGS
        return [json_clean(base)]
    if len(legs) == 1:
        return [json_clean(_settled_leg_record(base, legs[0], grades[0] if grades else None))]
    parlay_price = american_from_decimal(properties.get("odds"))
    records = []
    for index, leg in enumerate(legs):
        record = _as_parlay_leg(_settled_base_record(wager, feed, native_id, status, fetched_at),
                                base["id"], index, len(legs), parlay_price)
        grade = grades[index] if index < len(grades) else None
        records.append(json_clean(_settled_leg_record(record, leg, grade)))
    return records


def feeds_by_native_id(records: list[dict]) -> dict[str, str]:
    """nativeId -> the feed its Bet105 record sits under, for every record given (the
    stored ones, then this push's open ones, which win). Ticket numbers are one sequence
    across both feeds, so the native id alone names the ticket."""
    feeds: dict[str, str] = {}
    for record in records:
        raw = record.get("raw") or {}
        if record.get("venue") == VENUE and raw.get("nativeId") and raw.get("feed") in FEEDS:
            feeds[str(raw["nativeId"])] = raw["feed"]
    return feeds


def normalize_settled(wagers: list[dict], fetched_at: str | None,
                      known_feeds: dict[str, str]) -> tuple[list[dict], list[str]]:
    """Every graded wager -> records; plus, for each graded-looking wager it cannot read,
    the reason (the service logs them). A wager settles under the feed its ticket's record
    already sits under (`known_feeds`, feeds_by_native_id), else its productCode's — so it
    lands on its open record's id. Pending wagers are skipped silently: the open read owns
    open bets."""
    records: list[dict] = []
    skipped: list[str] = []
    for wager in wagers:
        status = settled_status_of(wager)
        ticket = str(wager.get("ticketNumber")).strip()
        if status is None:
            if _text(wager.get("wagerStatus")).upper() != WAGER_STATUS_PENDING:
                skipped.append(f"ticket {ticket}: unknown wagerStatus {wager.get('wagerStatus')!r}")
            continue
        feed = known_feeds.get(ticket) or FEED_BY_PRODUCT.get(_text(wager.get("productCode")))
        if feed is None:
            skipped.append(f"ticket {ticket}: unknown productCode {wager.get('productCode')!r} and no open record to take the feed from")
            continue
        records.extend(normalize_settled_wager(wager, feed, status, fetched_at))
    return records, skipped


def merge_settled(open_records: list[dict], settled_records: list[dict]) -> list[dict]:
    """The settled read's records over the open read's: a ticket the venue has graded is
    settled even while getHistory still lists it, so every open record of it gives way."""
    settled_tickets = {record["raw"]["nativeId"] for record in settled_records}
    still_open = [record for record in open_records if record["raw"]["nativeId"] not in settled_tickets]
    return still_open + settled_records


# ---- the push -------------------------------------------------------------------------

def _validate_feeds(feeds: object) -> dict | str:
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
    return clean


def _validate_settled(settled: object) -> list[dict] | str:
    if not isinstance(settled, list):
        return f"settled must be a list, got {type(settled).__name__}"
    for index, wager in enumerate(settled):
        if not isinstance(wager, dict):
            return f"settled[{index}] must be an object, got {type(wager).__name__}"
        ticket = wager.get("ticketNumber")
        if isinstance(ticket, bool) or not isinstance(ticket, (str, int)) or str(ticket).strip() in ("", "0"):
            return f"settled[{index}] carries no ticketNumber; keys seen: {sorted(wager.keys())}"
    return settled


def validate_push(body: object) -> dict | str:
    """{fetchedAt, feeds: {prematch: [...], live: [...]}, settled: [...]} or {error} of a
    POST /bet105.json body, or what is wrong with it. Every read must be present: a push
    is the account's whole open list and its whole wager list, or it is an error — never
    half (a half push would close real bets)."""
    if not isinstance(body, dict):
        return f"body must be an object, got {type(body).__name__}"
    if "error" in body:
        if not isinstance(body["error"], str) or not body["error"].strip():
            return "error must be a non-empty string"
        return {"error": body["error"].strip()}
    fetched_at = body.get("fetchedAt")
    if not isinstance(fetched_at, str) or parse_iso_ms(fetched_at) is None:
        return f"fetchedAt must be an ISO time, got {fetched_at!r}"
    feeds = _validate_feeds(body.get("feeds"))
    if isinstance(feeds, str):
        return feeds
    settled = _validate_settled(body.get("settled"))
    if isinstance(settled, str):
        return settled
    return {"fetchedAt": fetched_at, "feeds": feeds, "settled": settled}


def closed_by_absence(stored_records: list[dict], pushed_ids: set[str], now_iso: str) -> list[dict]:
    """The stored open Bet105 records a complete push no longer lists — in neither the
    open nor the settled read — as closed copies (no result: the venue shows none). Pure."""
    closed = []
    for record in stored_records:
        if record.get("venue") != VENUE or record.get("status") != "open" or record.get("id") in pushed_ids:
            continue
        copy = json_clean({**record, "raw": {**record.get("raw", {}), "closedReason": REASON_CLOSED_BY_ABSENCE}})
        copy.update({"status": "closed", "closedAt": now_iso, "sourceFetchedAt": now_iso})
        closed.append(copy)
    return closed
