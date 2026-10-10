"""DraftKings bet source: the open and settled bets the extension reads from the account's
own Chrome -> normalised bet records (2026-10-10).

DraftKings sits behind Akamai's bot checks and asks for a code on a new login, so — like
Bet105 — the network half is not here. The Unabated Ticket extension, running in Cal's
logged-in Chrome, mints My Bets' token, asks DraftKings' bets socket for every Open and
Settled bet page and POSTs them to this service's `POST /draftkings.json`
(extension/draftkings.js and panel.js; service.py). This module is the pure half plus the
service rules:

  normalize_draftkings(push)          every bet of both lists -> records, one per leg of a
                                      parlay; a settled bet over its open copy.
  validate_push(body)                 the POST body's shape, or what is wrong with it.
  closed_by_absence(stored, ids, now) an open record neither list carries any more is
                                      marked `closed` with no result — never guessed.
Pinned by tests/test_draftkings.py on tests/fixtures/bets/draftkings_bets.json.

Inputs:  {fetchedAt, open: [bet], settled: [bet], events: {eventId: event}}, or {error}
         when the extension could not read the account.
Outputs: records. Side effects: none (the service stores them).

Shapes — Cal's My Bets HAR of 2026-10-10 (three open CFB alternate spreads) and the page's
bundle (@draftkings/dk-my-bets-list 3.11.0), which names every field and status:
Bet: betId (the native id), receiptId, type ("Single" seen; parlays carry several
  selections), status ("Unsettled" while open), settlementStatus (the bundle's table:
  Open, Won, Lost, Cancelled, CashOut, PartialCashOut, Draw, NonRunner, HalfWon, HalfLost,
  WonDeadHeat, Placed, PladedDeadHeat, None), numberOfBets (> 1 = a round robin / system
  bet), displayOdds ("+566"; DraftKings writes a minus as U+2212), placementDate /
  settlementDate (ISO UTC), stake, potentialReturns (stake included: $75 at +566 ->
  499.5), returns (what it paid, stake included), plus the extension's freeBetAmount
  (bonus.freeBetAmount) and combinationCount.
Selection: eventId, marketId, displayOdds, selectionDisplayName ("Temple -17.5"),
  marketDisplayName ("Spread Alternate"), participants [{id, name}] (the picked team;
  empty on a total), nestedSelectionCount (an SGP group inside a parlay).
Event: eventStartDate (ISO UTC), homeTeamName, awayTeamName, leagueId, participants
  [{id, name, venueRole: "Home" | "Away"}].

Verified on the capture: open single alternate spreads of both venue roles' shape (all
three were the away or home favourite). Pinned only by hand-written fixture rows: every
settled status, a total, a moneyline, a first-half line, a parlay, a free bet. Fails closed
with the reason on anything else: a market name the tables do not know (team totals, props,
quarters' names not seen), a league id not listed (added once seen on a bet, never
guessed), a selection whose team is neither side of its event, an SGP group in a parlay, a
round robin.
"""
import re
from datetime import datetime, timezone

from unabated_ticket.bets_service.normalize import EASTERN, json_clean, parse_iso_ms, round_cents

SOURCE = "draftkings_extension"
VENUE = "draftkings"

# DraftKings league ids (nfl_specials/dk_book.py: 88808; mlb_sgp: 84240; the capture: 87637
# for college football; the rest are DraftKings' long-standing ids for these leagues).
LEAGUES = {"88808": "nfl", "87637": "cfb", "84240": "mlb", "42648": "nba", "92483": "cbb",
           "42133": "nhl", "94682": "wnba"}
# A market name is a period phrase (anywhere, or none = full game) plus one of these.
SPREAD_MARKETS = {"spread", "spread alternate", "alternate spread", "point spread", "run line",
                  "alternate run line", "run line alternate", "puck line", "alternate puck line",
                  "puck line alternate"}
TOTAL_MARKETS = {"total", "total alternate", "alternate total", "total points", "total runs",
                 "total goals", "alternate total runs", "alternate total points"}
MONEYLINE_MARKETS = {"moneyline", "money line"}
PERIOD_PHRASES = (("1st half", "1H"), ("2nd half", "2H"), ("1st quarter", "1Q"), ("2nd quarter", "2Q"),
                  ("3rd quarter", "3Q"), ("4th quarter", "4Q"), ("1st 5 innings", "F5"))
STATUS_BY_SETTLEMENT = {"Open": "open", "Won": "won", "Lost": "lost", "Draw": "push",
                        "Cancelled": "void", "NonRunner": "void", "CashOut": "closed",
                        "PartialCashOut": "closed"}
# Settled with a payout but no plain result (dead heats, half wins): kept as closed with
# the venue's P&L, never forced into won or lost.
DEAD_HEAT_STATUSES = {"HalfWon", "HalfLost", "WonDeadHeat", "Placed", "PladedDeadHeat"}
LIST_NAMES = ("open", "settled")
MAX_BETS_PER_PUSH = 2000
REASON_CLOSED_BY_ABSENCE = "left the open list and the settled list does not carry it"
REASON_NO_SELECTIONS = "bet carries no selections"
REASON_ROUND_ROBIN = "round robin / system bet (several bets on one ticket)"
REASON_SGP_GROUP = "SGP group inside a parlay"
UNICODE_MINUS = "−"
TOTAL_SELECTION = re.compile(r"^(over|under)\s+(\d+(?:\.\d+)?)$", re.IGNORECASE)
SPREAD_SELECTION = re.compile(r"^(.+?)\s+([+-]\d+(?:\.\d+)?)$")


# ---- small readers ------------------------------------------------------------------

def _text(value: object) -> str:
    return value.strip() if isinstance(value, str) else ""


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


def american_of(display_odds: object) -> int | None:
    """"+566" -> 566, "−110" (U+2212) or "-110" -> -110, "EVEN" -> 100; None otherwise."""
    text = _text(display_odds).replace(UNICODE_MINUS, "-")
    if text.upper() == "EVEN":
        return 100
    digits = text[1:] if text[:1] in ("+", "-") else text
    if not digits.isdigit():
        return None
    value = int(text)
    return value if abs(value) >= 100 else None


def moment_of(text: object) -> datetime | None:
    """ISO UTC with its zone -> datetime; None when absent or zoneless (never guessed)."""
    if not _text(text):
        return None
    try:
        moment = datetime.fromisoformat(_text(text).replace("Z", "+00:00"))
    except ValueError:
        return None
    return None if moment.tzinfo is None else moment.astimezone(timezone.utc)


def iso_utc(moment: datetime | None) -> str | None:
    return None if moment is None else moment.strftime("%Y-%m-%dT%H:%M:%SZ")


def _signed(text: str) -> float:
    return float(text.replace(UNICODE_MINUS, "-"))


# ---- one selection --------------------------------------------------------------------

def market_of(market_name: str) -> tuple[str, str] | str:
    """"1st Half Spread" -> ("spread", "1H"); "Total Alternate" -> ("total", "FG"); or the
    reason it is not a game line."""
    lowered = re.sub(r"\s+", " ", market_name.lower()).strip()
    period = "FG"
    for phrase, label in PERIOD_PHRASES:
        if phrase in lowered:
            period = label
            lowered = lowered.replace(phrase, " ")
            break
    base = re.sub(r"\s+", " ", lowered.replace("-", " ")).strip()
    if base in SPREAD_MARKETS:
        return "spread", period
    if base in TOTAL_MARKETS:
        return "total", period
    if base in MONEYLINE_MARKETS:
        return "moneyline", period
    return f"market {market_name or 'blank'!r} is not a game line"


def team_side_of(selection: dict, event: dict, picked_name: str) -> str | None:
    """"away" / "home" for the picked team: by participant id against the event's venue
    roles, else by name against the event's team names."""
    roles = {str(p.get("id")): _text(p.get("venueRole")).lower() for p in event.get("participants") or []}
    for participant in selection.get("participants") or []:
        role = roles.get(str(participant.get("id")))
        if role in ("away", "home"):
            return role
    if picked_name and picked_name == _text(event.get("awayTeamName")):
        return "away"
    if picked_name and picked_name == _text(event.get("homeTeamName")):
        return "home"
    return None


def parse_selection(selection: dict, event: dict | None) -> dict | str:
    """One selection and its event -> {league, period, betType, side, points, price,
    awayTeam, homeTeam, eventStart}, or the reason it is not a game line."""
    if selection.get("nestedSelectionCount"):
        return REASON_SGP_GROUP
    if not event:
        return f"event {selection.get('eventId')!r} missing from the push"
    league = LEAGUES.get(_text(str(event.get("leagueId") or "")))
    if league is None:
        return f"league id not supported ({event.get('leagueId')!r})"
    market = market_of(_text(selection.get("marketDisplayName")))
    if isinstance(market, str):
        return market
    bet_type, period = market
    away, home = _text(event.get("awayTeamName")), _text(event.get("homeTeamName"))
    if not away or not home:
        return "event names no teams"
    price = american_of(selection.get("displayOdds"))
    if price is None:
        return f"no price (displayOdds {selection.get('displayOdds')!r})"
    name = _text(selection.get("selectionDisplayName")).replace(UNICODE_MINUS, "-")
    parsed = {"league": league, "period": period, "betType": bet_type, "price": price, "awayTeam": away,
              "homeTeam": home, "eventStart": moment_of(event.get("eventStartDate")), "points": None}
    if bet_type == "total":
        match = TOTAL_SELECTION.match(name)
        if not match:
            return f"total selection {name!r} is not Over / Under a number"
        parsed.update(side=match.group(1).lower(), points=float(match.group(2)))
        return parsed
    picked_name = name
    if bet_type == "spread":
        match = SPREAD_SELECTION.match(name)
        if not match:
            return f"spread selection {name!r} carries no signed number"
        picked_name = match.group(1).strip()
        parsed["points"] = _signed(match.group(2))
    side = team_side_of(selection, event, picked_name)
    if side is None:
        return f"selection {name!r} is neither {away!r} nor {home!r}"
    parsed["side"] = side
    return parsed


# ---- one bet --------------------------------------------------------------------------

def native_id_of(bet: dict) -> str:
    value = bet.get("betId")
    if value in (None, "", 0, "0") or isinstance(value, bool):
        raise RuntimeError(f"DraftKings bet carries no betId; keys seen: {sorted(bet.keys())}")
    return str(value).strip()


def status_of(bet: dict) -> str:
    settlement = _text(bet.get("settlementStatus"))
    if settlement in DEAD_HEAT_STATUSES:
        return "closed"
    return STATUS_BY_SETTLEMENT.get(settlement, "unknown")


def stake_of(bet: dict) -> float | None:
    """A free bet stakes nothing of Cal's (the BFA / Bet105 convention)."""
    stake = _money(bet.get("stake"))
    free = _number(bet.get("freeBetAmount")) or 0
    if stake is None:
        return None
    return round_cents(max(stake - free, 0))


def to_win_of(bet: dict, status: str) -> float | None:
    """What the bet wins over its stake: the payout's profit once won, else the potential."""
    stake = _money(bet.get("stake"))
    paid = _money(bet.get("returns")) if status == "won" else _money(bet.get("potentialReturns"))
    if stake is None or paid is None:
        return None
    if (_number(bet.get("freeBetAmount")) or 0) > 0:
        # A free bet's returns carry no stake back.
        return paid
    return round_cents(paid - stake)


def pnl_of(bet: dict, status: str) -> float | None:
    """The venue's own result once settled (returns less Cal's stake); None while open."""
    if status in ("open", "unknown"):
        return None
    returns, stake = _money(bet.get("returns")), stake_of(bet)
    if returns is None or stake is None:
        return None
    return round_cents(returns - stake)


def _empty_record(record_id: str, fetched_at: str | None) -> dict:
    return {
        "id": record_id, "source": SOURCE, "venue": VENUE,
        "league": None, "eventStart": None, "eventDate": None,
        "awayTeam": None, "homeTeam": None, "awayKey": None, "homeKey": None,
        "rotation": None, "betType": "other", "period": None, "side": None, "points": None,
        "price": None, "stake": None, "toWin": None, "contracts": None,
        "placedAt": None, "status": None, "closedAt": None,
        "isParlayLeg": False, "parlayId": None, "legIndex": None, "legCount": None,
        "approx": [], "unmatchable": None, "sourceFetchedAt": fetched_at, "raw": {},
    }


def _base_record(bet: dict, native_id: str, read: str, fetched_at: str | None) -> dict:
    status = status_of(bet)
    record = _empty_record(f"{VENUE}:{native_id}", fetched_at)
    record.update({
        "stake": stake_of(bet),
        "toWin": to_win_of(bet, status),
        "pnl": pnl_of(bet, status),
        "placedAt": iso_utc(moment_of(bet.get("placementDate"))),
        "status": status,
        "closedAt": None if status == "open" else iso_utc(moment_of(bet.get("settlementDate"))),
        "raw": {
            "nativeId": native_id, "read": read, "receiptId": bet.get("receiptId"), "type": bet.get("type"),
            "betStatus": bet.get("status"), "settlementStatus": bet.get("settlementStatus"),
            "numberOfBets": bet.get("numberOfBets"), "displayOdds": bet.get("displayOdds"),
            "stake": bet.get("stake"), "potentialReturns": bet.get("potentialReturns"),
            "returns": bet.get("returns"), "freeBetAmount": bet.get("freeBetAmount"),
            "placementDate": bet.get("placementDate"), "settlementDate": bet.get("settlementDate"),
        },
    })
    return record


def _leg_record(record: dict, selection: dict, event: dict | None) -> dict:
    """The selection's own fields into raw (so an unmatched row shows them), then the
    parsed line — or the reason it failed."""
    record["raw"].update({
        "selectionId": selection.get("selectionId"), "eventId": selection.get("eventId"),
        "marketId": selection.get("marketId"), "selectionDisplayName": selection.get("selectionDisplayName"),
        "marketDisplayName": selection.get("marketDisplayName"), "legDisplayOdds": selection.get("displayOdds"),
        "legSettlementStatus": selection.get("settlementStatus"),
        "leagueId": (event or {}).get("leagueId"), "eventDisplayName": (event or {}).get("eventDisplayName"),
    })
    parsed = parse_selection(selection, event)
    if isinstance(parsed, str):
        record["unmatchable"] = parsed
        return record
    event_start = parsed["eventStart"]
    record.update({
        "league": parsed["league"], "eventStart": iso_utc(event_start),
        "eventDate": None if event_start is None else event_start.astimezone(EASTERN).strftime("%Y-%m-%d"),
        "awayTeam": parsed["awayTeam"], "homeTeam": parsed["homeTeam"], "betType": parsed["betType"],
        "period": parsed["period"], "side": parsed["side"], "points": parsed["points"], "price": parsed["price"],
    })
    return record


def normalize_bet(bet: dict, events: dict, read: str, fetched_at: str | None) -> list[dict]:
    """One bet -> one record, or one per leg of a parlay (every leg on the ticket's stake,
    payout and status, like Bet105). Pure."""
    native_id = native_id_of(bet)
    selections = bet.get("selections") or []
    base = _base_record(bet, native_id, read, fetched_at)
    if not selections:
        base["unmatchable"] = REASON_NO_SELECTIONS
        return [json_clean(base)]
    if (_number(bet.get("numberOfBets")) or 1) > 1 or (_number(bet.get("combinationCount")) or 0) > 0:
        base["unmatchable"] = REASON_ROUND_ROBIN
        return [json_clean(base)]
    if len(selections) == 1:
        return [json_clean(_leg_record(base, selections[0], events.get(str(selections[0].get("eventId")))))]
    parlay_price = american_of(bet.get("displayOdds"))
    records = []
    for index, selection in enumerate(selections):
        record = _base_record(bet, native_id, read, fetched_at)
        record["id"] = f"{base['id']}:leg{index}"
        record["raw"]["parlayPrice"] = parlay_price
        record.update({"isParlayLeg": True, "parlayId": base["id"], "legIndex": index, "legCount": len(selections)})
        records.append(json_clean(_leg_record(record, selection, events.get(str(selection.get("eventId"))))))
    return records


def normalize_draftkings(push: dict) -> tuple[list[dict], list[str]]:
    """Both lists -> records, plus the reasons for settled bets it would not store (the
    service logs them). A bet the settled list carries wins over its open copy; a settled
    bet with a status the bundle's table does not name is skipped, never guessed."""
    events = push["events"]
    fetched_at = push["fetchedAt"]
    settled_records: list[dict] = []
    skipped: list[str] = []
    for bet in push["settled"]:
        if status_of(bet) in ("open", "unknown"):
            skipped.append(f"bet {native_id_of(bet)}: settlementStatus {bet.get('settlementStatus')!r} in the settled list")
            continue
        settled_records.extend(normalize_bet(bet, events, "settled", fetched_at))
    settled_ids = {record["raw"]["nativeId"] for record in settled_records}
    open_records: list[dict] = []
    for bet in push["open"]:
        if native_id_of(bet) in settled_ids:
            continue
        open_records.extend(normalize_bet(bet, events, "open", fetched_at))
    return open_records + settled_records, skipped


# ---- the push -------------------------------------------------------------------------

def _validate_bets(name: str, bets: object) -> str | None:
    if not isinstance(bets, list):
        return f"{name} must be a list, got {type(bets).__name__}"
    for index, bet in enumerate(bets):
        if not isinstance(bet, dict):
            return f"{name}[{index}] must be an object, got {type(bet).__name__}"
        bet_id = bet.get("betId")
        if isinstance(bet_id, bool) or not isinstance(bet_id, (str, int)) or str(bet_id).strip() in ("", "0"):
            return f"{name}[{index}] carries no betId; keys seen: {sorted(bet.keys())}"
        if not isinstance(bet.get("selections", []), list):
            return f"{name}[{index}].selections must be a list"
    return None


def validate_push(body: object) -> dict | str:
    """{fetchedAt, open, settled, events} or {error} of a POST /draftkings.json body, or
    what is wrong with it. Both lists must be present: a push is the account's whole open
    list and its settled window, or it is an error — never half."""
    if not isinstance(body, dict):
        return f"body must be an object, got {type(body).__name__}"
    if "error" in body:
        if not isinstance(body["error"], str) or not body["error"].strip():
            return "error must be a non-empty string"
        return {"error": body["error"].strip()}
    fetched_at = body.get("fetchedAt")
    if not isinstance(fetched_at, str) or parse_iso_ms(fetched_at) is None:
        return f"fetchedAt must be an ISO time, got {fetched_at!r}"
    for name in LIST_NAMES:
        problem = _validate_bets(name, body.get(name))
        if problem:
            return problem
    if len(body["open"]) + len(body["settled"]) > MAX_BETS_PER_PUSH:
        return f"at most {MAX_BETS_PER_PUSH} bets per push, got {len(body['open']) + len(body['settled'])}"
    events = body.get("events")
    if not isinstance(events, dict) or not all(isinstance(event, dict) for event in events.values()):
        return f"events must be an object of objects, got {type(events).__name__}"
    return {"fetchedAt": fetched_at, "open": body["open"], "settled": body["settled"], "events": events}


def closed_by_absence(stored_records: list[dict], pushed_ids: set[str], now_iso: str) -> list[dict]:
    """The stored open DraftKings records a complete push no longer lists — in neither list —
    as closed copies (no result: none is known). Pure."""
    closed = []
    for record in stored_records:
        if record.get("venue") != VENUE or record.get("status") != "open" or record.get("id") in pushed_ids:
            continue
        copy = json_clean({**record, "raw": {**record.get("raw", {}), "closedReason": REASON_CLOSED_BY_ABSENCE}})
        copy.update({"status": "closed", "closedAt": now_iso, "sourceFetchedAt": now_iso})
        closed.append(copy)
    return closed
