"""Novig bet source (issue #116; rebuilt 2026-09-22 on the Portfolio REST
feed): the account's Portfolio -> normalised bet records, polled from the
bets service with no browser involved.

Why REST: on 2026-09-22 Novig moved its web app from app.novig.us to
novig.com and put its Hasura GraphQL behind a query allowlist the same hour
("query is not allowed", even for the old app's own queries). The new app's
Portfolio screen reads a REST feed on the app's own Auth0 token, and so does
this source — the same client id and audience NovigAuth already mints for,
no allowlist in the path.

    GET {NOVIG_REST_URL}/portfolio/{active|settled}
        ?currency=CASH&sort={active|settled}_recency&limit=N[&cursor=...]
        Authorization: Bearer <access token>        (NovigAuth)
        novig-client-capabilities: card_stack,player_props   (what the app sends)
    -> {"items": [card, ...], "nextCursor": str | null}   cursor-paged to the end

Cards (the app's *PortfolioCardDTO shapes, live-captured 2026-09-22 and
pinned in tests/fixtures/bets/novig_portfolio.json):
  straight  one order: {orderId, state, outcome, price, amounts, qty,
            eventHeader, metadata.createdAt, sortAt}; cardId is the market id
  union     several orders on one market — the same order-line shape under
            stateBlocks[].{home,away}.lines[] (active) and settledItems[].line
            (settled), sharing the card's eventHeader and marketLabel
  parlay    {state, price, amountCopy, payoutCopy, legCount, sgpGroups[].legs[],
            singleLegs[]}; a leg is {outcomeId, marketId, eventHeader, outcome}
            and carries no price of its own; cardId is the parlay id

Novig facts these rules rest on (validated 2026-09-22 against the 476 orders
the GraphQL-era normaliser had already recorded):
  * the outcome on a line is the side HELD: a lay already reads as the other
    outcome, so there is no isBid to flip (476/476 agreed on side; the old
    lay flip had priced one order at the wrong side's 0.70 where the feed's
    own 0.305 is what was paid)
  * price is a 0-1 probability of that side; amounts.cost is the dollars
    risked; payoutCopy is the gross payout with thousands commas, which is
    the contract count (a contract pays $1); qty is what is still RESTING,
    in hundredths of a contract
  * one order can appear twice: its matched part (Matched / Win / Loss) and a
    Canceled line for the unmatched remainder under the same orderId. The
    matched line is the bet; a Canceled line is a bet only when the order
    never matched at all (void, $0)
  * outcome.index 0 = home / Over; competitor.id is the team; subtitle is
    "Moneyline" / "Spread" / "Total", prefixed "1st Half " for the first half
    or "First 5 " for MLB's first five innings; anything else is a prop — not
    a game market
  * every team carries `unabatedId`, Unabated's own team id: it rides in
    awayTeamVenue / homeTeamVenue and bets.js keys the record on it before
    falling back to names (#118)
  * states seen live: Matched, Win, Loss, Canceled, and Settled for a
    FRACTIONAL settlement (outcome.status "0.72", "0.50" — a first-five tie
    paid at half, a partial), graded by what it paid against its cost, with
    toWin the actual P&L; Draw / Push, Cash Out and Unmatched / Pending are
    mapped by name; anything else is "unknown" and logged — never silently
    open

Inputs:  the GETs above. Outputs: list of records (contract in
docs/2026-09-11-issue-114-bet-history-plan.md; ids "novig:<order id>" and
"novig:<parlay id>:<leg index>" as before, so bets.duckdb history continues).
Side effects: NovigAuth may rewrite the token file on rotation; nothing else
on disk. Raises on any failed request, malformed page or an unended list, so
the service records a failed run and keeps the previous records.
"""
import logging
import re
from collections.abc import Callable
from datetime import datetime, timezone
from pathlib import Path
from zoneinfo import ZoneInfo

from unabated_ticket.bets_service import config
from unabated_ticket.bets_service.normalize import js_round, json_clean, utc_now_iso
from unabated_ticket.bets_service.sources.novig_auth import NovigAuth, NovigAuthError

log = logging.getLogger(__name__)

SOURCE = "novig_rest"
VENUE = "novig"
EASTERN = ZoneInfo("America/New_York")

# (tab, sort) of the two Portfolio lists the app reads; both are paged to the end.
PORTFOLIO_TABS = (("active", "active_recency"), ("settled", "settled_recency"))
PAGE_LIMIT = 50
MAX_PAGES = 200
CLIENT_CAPABILITIES_HEADER = ("novig-client-capabilities", "card_stack,player_props")

# Novig eventHeader.league -> the feed.LEAGUES path the scanner uses.
LEAGUES = {
    "NFL": "nfl", "NCAAF": "cfb", "NBA": "nba", "NCAAB": "cbb", "WNBA": "wnba", "MLB": "mlb", "NHL": "nhl",
    "MLS": "soccer", "EPL": "soccer", "Bundesliga": "soccer", "Serie A": "soccer", "La Liga": "soccer",
    "Ligue 1": "soccer", "Champions League": "soccer", "Europa League": "soccer", "FIFA Club World Cup": "soccer",
}
# outcome.subtitle of a game market: an optional period prefix + the bet type.
GAME_MARKET_RE = re.compile(r"^(1st Half |First 5 )?(Moneyline|Spread|Total)$")
PERIOD_BY_PREFIX = {None: "FG", "1st Half ": "1H", "First 5 ": "F5"}
BET_TYPE_BY_SUBTITLE = {"Moneyline": "moneyline", "Spread": "spread", "Total": "total"}
TRAILING_NUMBER_RE = re.compile(r"([+-]?\d+(?:\.\d+)?)\s*$")
LEADING_NUMBER_RE = re.compile(r"^\s*([+-]?\d+(?:\.\d+)?)")
REASON_NOT_GAME = "not a game market"
APPROX_UNMATCHED = "novig_order_unmatched"
QTY_PER_CONTRACT = 100

# Card / line `state` (case- and space-insensitive) -> the contract's status.
STATUS_BY_STATE = {
    "matched": "open", "win": "won", "loss": "lost", "draw": "push", "push": "push",
    "canceled": "void", "cancelled": "void", "cashout": "closed", "cashedout": "closed",
    "unmatched": "open", "pending": "open", "open": "open",
}
# States of an order that has nothing matched yet: open on its resting size, flagged.
RESTING_STATES = {"unmatched", "pending", "open"}
CANCELED_STATES = {"canceled", "cancelled"}
# A fractional settlement: graded on the money, not the state name.
FRACTIONAL_STATE = "settled"


# ---- helpers ------------------------------------------------------------------------

def _to_number(value: object) -> float | None:
    """A finite number from a number or a string ("1,279.25" — the copy fields
    carry thousands commas), else None."""
    if isinstance(value, bool) or value is None:
        return None
    try:
        number = float(str(value).replace(",", "")) if isinstance(value, str) else float(value)
    except (TypeError, ValueError):
        return None
    return number if number == number and number not in (float("inf"), float("-inf")) else None


def _iso_utc(value: object) -> str | None:
    """Novig sends "...Z" already; re-emit at millisecond precision like the other sources."""
    if not isinstance(value, str):
        return None
    try:
        parsed = datetime.fromisoformat(value.replace("Z", "+00:00"))
    except ValueError:
        return None
    if parsed.tzinfo is None:
        parsed = parsed.replace(tzinfo=timezone.utc)
    parsed = parsed.astimezone(timezone.utc)
    return parsed.strftime("%Y-%m-%dT%H:%M:%S.") + f"{parsed.microsecond // 1000:03d}Z"


def _eastern_date_of(iso: str) -> str | None:
    try:
        parsed = datetime.fromisoformat(iso.replace("Z", "+00:00"))
    except ValueError:
        return None
    return parsed.astimezone(EASTERN).strftime("%Y-%m-%d")


def probability_to_american(probability: float | None) -> int | None:
    if probability is None or not (0 < probability < 1):
        return None
    if probability >= 0.5:
        return -js_round(probability / (1 - probability) * 100)
    return js_round((1 - probability) / probability * 100)


def _round_cents(dollars: float) -> float:
    return js_round(dollars * 100) / 100


def _state_key(state: object) -> str:
    return re.sub(r"[\s_]", "", state).lower() if isinstance(state, str) else ""


def _status_of(state: object, cost: float | None = None, payout: float | None = None) -> str:
    """The contract's status for a card / line state; "unknown" for a state
    this module has not seen (normalize_novig logs those once per poll)."""
    key = _state_key(state)
    if key == FRACTIONAL_STATE:
        if cost is None or payout is None:
            return "unknown"
        if abs(payout - cost) < 0.005:
            return "push"
        return "won" if payout > cost else "lost"
    return STATUS_BY_STATE.get(key, "unknown")


def _venue_id(value: object) -> str | None:
    return value if isinstance(value, str) and value != "" else None


def _venue_team(team: dict | None) -> dict | None:
    """The team as Novig names it, plus Novig's copy of Unabated's team id."""
    if not isinstance(team, dict):
        return None
    unabated_id = team.get("unabatedId")
    return {
        "id": _venue_id(team.get("id")), "name": team.get("name") or None,
        "shortName": team.get("shortName") or None, "symbol": team.get("symbol") or None,
        "unabatedId": str(unabated_id) if isinstance(unabated_id, int) and not isinstance(unabated_id, bool) else None,
    }


def _team_name(team: dict | None) -> str | None:
    if not isinstance(team, dict):
        return None
    return team.get("name") or team.get("shortName") or team.get("symbol") or None


def _trailing_number(text: object) -> float | None:
    match = TRAILING_NUMBER_RE.search(text.strip()) if isinstance(text, str) else None
    return _to_number(match.group(1)) if match else None


def _leading_number(text: object) -> float | None:
    match = LEADING_NUMBER_RE.search(text) if isinstance(text, str) else None
    return _to_number(match.group(1)) if match else None


# ---- side / points ------------------------------------------------------------------

def _team_side(outcome: dict, header: dict) -> str | None:
    """"away" / "home" for the team an outcome names: its competitor id
    against the event's teams first, Novig's index convention (0 = home)
    second. None when neither resolves — a wrong side flags the wrong line."""
    competitor_id = _venue_id((outcome.get("competitor") or {}).get("id"))
    home_id = _venue_id((header.get("homeTeam") or {}).get("id"))
    away_id = _venue_id((header.get("awayTeam") or {}).get("id"))
    if competitor_id and competitor_id == home_id:
        return "home"
    if competitor_id and competitor_id == away_id:
        return "away"
    index = outcome.get("index")
    if index == 0:
        return "home"
    if index == 1:
        return "away"
    return None


def _total_side(outcome: dict) -> str | None:
    title = (outcome.get("title") or "").strip().lower()
    label = (outcome.get("displayLabel") or "").strip().lower()
    if title.startswith("over") or label.startswith("o ") or label.startswith("over"):
        return "over"
    if title.startswith("under") or label.startswith("u ") or label.startswith("under"):
        return "under"
    index = outcome.get("index")
    return "over" if index == 0 else "under" if index == 1 else None


def _outcome_number(outcome: dict, market_label: str | None) -> float | None:
    """The line's number: the outcome title's trailing number ("ATL -1.5",
    "Under 58.5"), else its display label ("O 7.5", "-1.5"), else the union
    card's market label ("7.5 Total") — union lines title just "Over"."""
    number = _trailing_number(outcome.get("title"))
    if number is None:
        number = _trailing_number(outcome.get("displayLabel"))
    if number is None:
        number = _leading_number(market_label)
    return number


def _side_and_points(bet_type: str, outcome: dict, header: dict, market_label: str | None) -> dict | str:
    """{side, points} of the side held, or an error string."""
    outcome_id = outcome.get("outcomeId") or "?"
    if bet_type == "total":
        side = _total_side(outcome)
        points = _outcome_number(outcome, market_label)
        if side is None:
            return f"total outcome {outcome_id} names neither over nor under"
        if points is None:
            return f"total outcome {outcome_id} has no number"
        return {"side": side, "points": points}
    side = _team_side(outcome, header)
    if side is None:
        return f"outcome {outcome_id} names no team"
    if bet_type == "moneyline":
        return {"side": side, "points": None}
    points = _outcome_number(outcome, market_label)
    if points is None:
        return f"spread outcome {outcome_id} has no number"
    return {"side": side, "points": points}


# ---- records ------------------------------------------------------------------------

def _unmatchable(base: dict, reason: str) -> dict:
    base.update({
        "league": None, "eventStart": None, "eventDate": None, "awayTeam": None, "homeTeam": None,
        "awayKey": None, "homeKey": None, "betType": "other", "period": None, "side": None, "points": None,
        "unmatchable": reason,
    })
    return base


def _game_fields(outcome: dict, header: dict, market_label: str | None) -> dict:
    """league / event / teams / market of a line or a parlay leg, or {unmatchable}."""
    if header.get("eventType") not in (None, "Game"):
        return {"unmatchable": REASON_NOT_GAME}
    match = GAME_MARKET_RE.match(outcome.get("subtitle") or "")
    if not match:
        return {"unmatchable": REASON_NOT_GAME}
    home, away = header.get("homeTeam"), header.get("awayTeam")
    if not isinstance(home, dict) or not isinstance(away, dict):
        return {"unmatchable": f"unreadable Novig card (no teams on event {header.get('eventId') or '?'})"}
    novig_league = header.get("league")
    league = LEAGUES.get(novig_league) if isinstance(novig_league, str) else None
    if not league:
        return {"unmatchable": f"league not supported ({novig_league or 'unknown'})"}
    bet_type = BET_TYPE_BY_SUBTITLE[match.group(2)]
    taken = _side_and_points(bet_type, outcome, header, market_label)
    if isinstance(taken, str):
        return {"unmatchable": f"unreadable Novig card ({taken})"}
    event_start = _iso_utc(header.get("scheduledStart"))
    period = PERIOD_BY_PREFIX[match.group(1)]
    # Novig labels MLB's first five innings "First 5"; an MLB "1st Half" is the same thing.
    if period == "1H" and league == "mlb":
        period = "F5"
    return {
        "league": league, "eventStart": event_start, "eventDate": _eastern_date_of(event_start) if event_start else None,
        "awayTeam": _team_name(away), "homeTeam": _team_name(home), "awayKey": None, "homeKey": None,
        "betType": bet_type, "period": period,
        "side": taken["side"], "points": taken["points"], "unmatchable": None,
    }


def _venue_fields(market_id: object, outcome: dict, header: dict) -> dict:
    return {
        "venueIds": {"marketId": _venue_id(market_id), "outcomeId": _venue_id(outcome.get("outcomeId")),
                     "eventId": _venue_id(header.get("eventId")), "gameId": None},
        "awayTeamVenue": _venue_team(header.get("awayTeam")), "homeTeamVenue": _venue_team(header.get("homeTeam")),
    }


def _market_title(outcome: dict, header: dict) -> str:
    return " · ".join(part for part in (header.get("description"), outcome.get("subtitle"), outcome.get("title")) if part)


def order_lines_of(card: dict) -> list[tuple[dict, dict, str | None]]:
    """(line, eventHeader, marketLabel) for every order line a straight or
    union card carries — a straight card is its own single line."""
    header = card.get("eventHeader") or {}
    if card.get("type") == "straight":
        return [(card, header, None)]
    if card.get("type") != "union":
        return []
    label = card.get("marketLabel")
    # A union's lines carry no sortAt of their own; the card's (its latest
    # settlement) is the closest thing to when each line closed.
    sort_at = card.get("sortAt")
    lines = []
    for block in card.get("stateBlocks") or []:
        for side in ("home", "away"):
            for line in ((block.get(side) or {}).get("lines") or []):
                lines.append(({"sortAt": sort_at, **line}, header, label))
    for item in card.get("settledItems") or []:
        if isinstance(item, dict) and isinstance(item.get("line"), dict):
            lines.append(({"sortAt": sort_at, **item["line"]}, header, label))
    return lines


def normalize_order_line(line: dict, header: dict, market_id: object, market_label: str | None,
                         fetched_at: str | None) -> dict:
    """One order line -> one record. `market_id` is the card's id (a straight
    or union card's cardId is its market id)."""
    outcome = line.get("outcome") or {}
    amounts = line.get("amounts") or {}
    state_key = _state_key(line.get("state"))
    probability = _to_number(line.get("price"))
    cost = _to_number(amounts.get("cost"))
    # payoutCopy is the gross payout: contracts x $1 (or what a fractional
    # settlement actually paid); absent on a Canceled line.
    payout = _to_number(amounts.get("payoutCopy"))
    if payout is None and cost is not None and probability:
        payout = cost / probability
    status = _status_of(line.get("state"), cost, payout)
    resting = state_key in RESTING_STATES
    # A Canceled line is the unmatched remainder (or a never-matched order):
    # nothing is at risk on it whatever its resting cost reads.
    if status == "void":
        stake, contracts, to_win = 0.0, 0.0, 0.0
    else:
        stake = _round_cents(cost) if cost is not None else None
        contracts = _round_cents(cost / probability) if cost is not None and probability else None
        to_win = _round_cents(payout - cost) if payout is not None and cost is not None else None
    placed_at = _iso_utc((line.get("metadata") or {}).get("createdAt"))
    resting_qty = _to_number(line.get("qty"))
    base = {
        "id": f"novig:{line.get('orderId')}",
        "source": SOURCE,
        "venue": VENUE,
        "rotation": None,
        "price": probability_to_american(probability),
        "stake": stake,
        "toWin": to_win,
        "contracts": contracts,
        "placedAt": placed_at,
        "status": status,
        "closedAt": None if status == "open" else _iso_utc(line.get("sortAt")) or placed_at,
        "isParlayLeg": False, "parlayId": None, "legIndex": None, "legCount": None,
        "approx": [APPROX_UNMATCHED] if resting else [],
        "sourceFetchedAt": fetched_at,
        **_venue_fields(market_id, outcome, header),
        "raw": {
            "orderId": line.get("orderId"), "state": line.get("state"), "probability": probability, "cost": cost,
            "payout": payout, "restingQty": resting_qty, "qtyUnit": QTY_PER_CONTRACT,
            "outcomeIndex": outcome.get("index"), "outcomeTitle": outcome.get("title") or None,
            "outcomeSubtitle": outcome.get("subtitle") or None, "outcomeStatus": outcome.get("status") or None,
            "novigLeague": header.get("league"), "marketTitle": _market_title(outcome, header),
        },
    }
    game = _game_fields(outcome, header, market_label)
    if game.get("unmatchable"):
        return _unmatchable(base, game["unmatchable"])
    base.update(game)
    return base


def parlay_legs_of(card: dict) -> list[dict]:
    """Every leg of a parlay card in display order: the same-game groups' legs, then the singles."""
    legs = []
    for group in card.get("sgpGroups") or []:
        legs.extend(leg for leg in (group.get("legs") or []) if isinstance(leg, dict))
    legs.extend(leg for leg in (card.get("singleLegs") or []) if isinstance(leg, dict))
    return legs


def normalize_parlay(card: dict, fetched_at: str | None) -> list[dict]:
    legs = parlay_legs_of(card)
    status = _status_of(card.get("state"))
    wager = _to_number(card.get("amountCopy"))
    payout = _to_number(card.get("payoutCopy"))
    parlay_price = _to_number(card.get("price"))
    placed_at = _iso_utc((card.get("metadata") or {}).get("createdAt"))
    closed_at = None if status == "open" else _iso_utc(card.get("sortAt")) or placed_at
    records = []
    for leg_index, leg in enumerate(legs):
        outcome = leg.get("outcome") or {}
        header = leg.get("eventHeader") or {}
        base = {
            "id": f"novig:{card.get('cardId')}:{leg_index}",
            "source": SOURCE,
            "venue": VENUE,
            "rotation": None,
            "price": None,  # the feed prices the parlay, not its legs
            "stake": _round_cents(wager) if wager is not None else None,
            "toWin": _round_cents(payout - wager) if payout is not None and wager is not None else None,
            "contracts": None,
            "placedAt": placed_at,
            "status": status,
            "closedAt": closed_at,
            "isParlayLeg": True, "parlayId": f"novig:{card.get('cardId')}", "legIndex": leg_index, "legCount": len(legs),
            "approx": [],
            "sourceFetchedAt": fetched_at,
            **_venue_fields(leg.get("marketId"), outcome, header),
            "raw": {
                "parlayId": card.get("cardId"), "state": card.get("state"), "parlayPrice": parlay_price,
                "wager": wager, "payout": payout, "isSgp": card.get("isSgp"),
                "outcomeIndex": outcome.get("index"), "outcomeTitle": outcome.get("title") or None,
                "outcomeSubtitle": outcome.get("subtitle") or None, "outcomeStatus": outcome.get("status") or None,
                "novigLeague": header.get("league"), "marketTitle": _market_title(outcome, header),
            },
        }
        game = _game_fields(outcome, header, None)
        if game.get("unmatchable"):
            records.append(_unmatchable(base, game["unmatchable"]))
        else:
            base.update(game)
            records.append(base)
    return records


def normalize_novig(cards: list[dict], fetched_at: str | None) -> list[dict]:
    """Every record the Portfolio cards describe (active and settled lists
    together). An order seen as both its matched part and a Canceled
    remainder keeps the matched part; two lines of the same kind keep the
    later one."""
    matched: dict[str, dict] = {}
    canceled_only: dict[str, dict] = {}
    parlay_records: dict[str, dict] = {}
    for card in cards:
        if not isinstance(card, dict):
            continue
        if card.get("type") == "parlay":
            if not card.get("cardId"):
                continue
            for record in normalize_parlay(card, fetched_at):
                parlay_records[record["id"]] = record
            continue
        for line, header, market_label in order_lines_of(card):
            if not line.get("orderId"):
                continue
            record = normalize_order_line(line, header, card.get("cardId"), market_label, fetched_at)
            target = canceled_only if _state_key(line.get("state")) in CANCELED_STATES else matched
            target[record["id"]] = record
    for record_id, record in canceled_only.items():
        matched.setdefault(record_id, record)
    records = [*matched.values(), *parlay_records.values()]
    unknown_states = sorted({repr(record["raw"].get("state")) for record in records if record["status"] == "unknown"})
    if unknown_states:
        log.warning("novig: %d record(s) in a state this source does not know -> status 'unknown': %s",
                    sum(record["status"] == "unknown" for record in records), ", ".join(unknown_states))
    return [json_clean(record) for record in records]


# ---- network half -------------------------------------------------------------------

HttpGet = Callable[[str, str], dict]


def _make_session():
    """curl_cffi with Chrome impersonation, the transport the anonymous Novig
    SGP scraper already uses through Cloudflare; plain requests as a fallback."""
    try:
        from curl_cffi import requests as cffi_requests
        session = cffi_requests.Session(impersonate="chrome")
    except ImportError:  # pragma: no cover
        import requests
        session = requests.Session()
    return session


def portfolio_url(base_url: str, tab: str, sort: str, limit: int, cursor: str | None) -> str:
    url = f"{base_url}/portfolio/{tab}?currency=CASH&sort={sort}&limit={limit}"
    return f"{url}&cursor={cursor}" if cursor else url


class NovigSource:
    """Source protocol implementation for a Novig account (read-only GETs).

    `http_get(url, bearer) -> body` is injectable for tests; the default GETs
    through curl_cffi. `auth` is a NovigAuth over the token file.
    """

    name = "novig"

    def __init__(self, auth: NovigAuth | None = None, http_get: HttpGet | None = None,
                 poll_sec: float | None = None, base_url: str | None = None):
        self.poll_sec = poll_sec if poll_sec is not None else config.NOVIG_POLL_SEC
        self._base_url = (base_url or config.NOVIG_REST_URL).rstrip("/")
        self._auth = auth or NovigAuth(config.NOVIG_TOKEN_PATH)
        self._http_get = http_get or self._get_json
        self._session = None

    def _get_json(self, url: str, bearer: str) -> dict:
        if self._session is None:
            self._session = _make_session()
        response = self._session.get(url, headers={"Authorization": f"Bearer {bearer}", "Accept": "application/json",
                                                   CLIENT_CAPABILITIES_HEADER[0]: CLIENT_CAPABILITIES_HEADER[1]},
                                     timeout=20)
        if response.status_code != 200:
            raise RuntimeError(f"novig portfolio HTTP {response.status_code}: {response.text[:200]}")
        return response.json()

    def _page_to_end(self, tab: str, sort: str, bearer: str) -> list[dict]:
        cards: list[dict] = []
        cursor = None
        for _page in range(MAX_PAGES):
            body = self._http_get(portfolio_url(self._base_url, tab, sort, PAGE_LIMIT, cursor), bearer)
            items = body.get("items") if isinstance(body, dict) else None
            if not isinstance(items, list):
                raise RuntimeError(f"novig portfolio/{tab}: expected an `items` list, got {type(items).__name__}")
            cards.extend(items)
            cursor = body.get("nextCursor")
            if not cursor:
                return cards
        raise RuntimeError(f"novig portfolio/{tab} did not end within {MAX_PAGES} pages; refusing a partial pull")

    def fetch(self) -> list[dict]:
        try:
            access = self._auth.token()
        except NovigAuthError as error:
            raise RuntimeError(f"novig auth: {error}") from error
        cards: list[dict] = []
        counts = {}
        for tab, sort in PORTFOLIO_TABS:
            tab_cards = self._page_to_end(tab, sort, access.token)
            counts[tab] = len(tab_cards)
            cards.extend(tab_cards)
        records = normalize_novig(cards, utc_now_iso())
        log.info("novig: %d active + %d settled cards -> %d records", counts["active"], counts["settled"], len(records))
        return records


def source_if_connected(token_path: Path | None = None) -> NovigSource | None:
    """A NovigSource when `novig_auth connect` has been run (the token file
    exists), else None so the service leaves the venue "no source configured"."""
    path = token_path or config.NOVIG_TOKEN_PATH
    if not Path(path).exists():
        return None
    return NovigSource(auth=NovigAuth(path))
