"""Novig bet source (issue #116): the account's orders + parlays -> normalised
bet records, polled from the bets service with no browser involved.

Two halves:
  normalize_novig(...)  pure port of extension/novig_bets.js normalizeNovig;
                        the node tests on tests/fixtures/bets/novig_bets.json
                        pin its output and tests/test_parity_novig.py holds
                        the two byte-equivalent.
  NovigSource           the network half (Source protocol). Every poll: a
                        bearer from NovigAuth (the service's own refresh
                        token), then the `order` and `parlay` tables of
                        Novig's Hasura GraphQL, filtered to this trader and to
                        rows that are open or changed within the retention
                        window, paginated to the end. Selections are the
                        fragments the app's own Portfolio cards read
                        (OrderStateCard / OrderDescription / EventData /
                        ParlayStateCard), so the rows have the shape
                        novig_bets.js was written against.

Inputs:  POST https://api.novig.us/v1/graphql with Authorization: Bearer.
Outputs: list of records. Side effects: NovigAuth may rewrite the token file on
         rotation; nothing else on disk. Raises on any failed request or a
         GraphQL error so the service records a failed run and keeps the
         previous records.

Novig facts (bundle 2026-09-11, live-verified the same day — see the
fixture's "provenance" and novig_bets.js): outcome index 0 = HOME / OVER;
`price` is a 0-1 probability, one contract pays $1 and `qty` is in
HUNDREDTHS of a contract; `isBid` false LAYS the outcome (the other side at
1 - p, a spread's number negated); matched size = originalQty - qty; `_1H`
market types are the first five innings on MLB; outcome.status may be a
fractional settlement string; parlay.status is Title-case on the wire.
"""
import logging
import re
import time
from collections.abc import Callable
from datetime import datetime, timedelta, timezone
from pathlib import Path
from zoneinfo import ZoneInfo

from unabated_ticket.bets_service import config
from unabated_ticket.bets_service.normalize import js_round, json_clean, utc_now_iso
from unabated_ticket.bets_service.sources.novig_auth import NovigAuth, NovigAuthError

log = logging.getLogger(__name__)

SOURCE = "novig_page"
VENUE = "novig"
EASTERN = ZoneInfo("America/New_York")

LEAGUES = {
    "NFL": "nfl", "NCAAF": "cfb", "NBA": "nba", "NCAAB": "cbb", "WNBA": "wnba", "MLB": "mlb", "NHL": "nhl",
    "MLS": "soccer", "EPL": "soccer", "Bundesliga": "soccer", "Serie A": "soccer", "La Liga": "soccer",
    "Ligue 1": "soccer", "Champions League": "soccer", "Europa League": "soccer", "FIFA Club World Cup": "soccer",
}
MARKET_TYPES = {
    "MONEY": ("moneyline", False), "SPREAD": ("spread", False), "TOTAL": ("total", False),
    "MONEY_1H": ("moneyline", True), "SPREAD_1H": ("spread", True), "TOTAL_1H": ("total", True),
}
REASON_NOT_GAME = "not a game market"
APPROX_UNMATCHED = "novig_order_unmatched"
APPROX_PENDING = "novig_order_pending"
DESCRIPTION_NUMBER_RE = re.compile(r"([+-]?\d+(?:\.\d+)?)\s*$")
# Novig's qty unit: 100 = one $1-payout contract (live-verified 2026-09-11).
QTY_PER_CONTRACT = 100


# ---- pure normaliser (port of novig_bets.js) -----------------------------------------

def _to_number(value: object) -> float | None:
    """novig_bets.js toNumber: a finite number or null (NOT bets.js's 0)."""
    if isinstance(value, bool) or value is None:
        return None
    try:
        number = float(value)
    except (TypeError, ValueError):
        return None
    return number if number == number and number not in (float("inf"), float("-inf")) else None


def _iso_utc(value: object) -> str | None:
    """Date.parse + toISOString: millisecond precision, trailing Z."""
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


def _team_name(team: dict | None) -> str | None:
    if not team:
        return None
    return team.get("name") or team.get("short_name") or team.get("symbol") or None


def _matchup_of(market: dict) -> str | None:
    game = (market.get("event") or {}).get("game")
    return f"{_team_name(game.get('awayTeam'))} @ {_team_name(game.get('homeTeam'))}" if game else None


def _team_index_of(outcome: dict, game: dict) -> int | None:
    """0 = away / 1 = home: competitor symbol against the game's teams first,
    Novig's index convention (0 = home) second; None when neither resolves."""
    competitor = outcome.get("competitor") or {}
    symbol = competitor.get("symbol")
    home, away = game.get("homeTeam"), game.get("awayTeam")
    if symbol and home and away:
        if symbol == home.get("symbol"):
            return 1
        if symbol == away.get("symbol"):
            return 0
    if outcome.get("index") == 0:
        return 1
    if outcome.get("index") == 1:
        return 0
    return None


def _description_number(description: object) -> float | None:
    match = DESCRIPTION_NUMBER_RE.search(description.strip()) if isinstance(description, str) else None
    return _to_number(match.group(1)) if match else None


def _outcome_side(bet_type: str, outcome: dict, market: dict, game: dict) -> dict | str:
    """{side, points} of the OUTCOME (what a bid backs), or an error string."""
    if bet_type == "total":
        description = outcome["description"].strip().lower() if isinstance(outcome.get("description"), str) else ""
        if description.startswith("over"):
            side = "over"
        elif description.startswith("under"):
            side = "under"
        else:
            side = {0: "over", 1: "under"}.get(outcome.get("index"))
        points = _description_number(outcome.get("description"))
        if points is None:
            points = _to_number(market.get("strike"))
        if not side:
            return f"total outcome {outcome.get('id')} names neither over nor under"
        if points is None:
            return f"total outcome {outcome.get('id')} has no number"
        return {"side": side, "points": points}
    team_index = _team_index_of(outcome, game)
    if team_index is None:
        return f"outcome {outcome.get('id')} names no team"
    side = "away" if team_index == 0 else "home"
    if bet_type == "moneyline":
        return {"side": side, "points": None}
    points = _description_number(outcome.get("description"))
    if points is None:
        strike = _to_number(market.get("strike"))
        if strike is not None:
            points = strike if team_index == 1 else -strike
    if points is None:
        return f"spread outcome {outcome.get('id')} has no number"
    return {"side": side, "points": points}


def _lay_side(bet_type: str, backed: dict) -> dict:
    flip = {"away": "home", "home": "away", "over": "under", "under": "over"}
    return {"side": flip[backed["side"]], "points": -backed["points"] if bet_type == "spread" else backed["points"]}


def _grade_settled(outcome_status: object, is_bid: bool) -> str:
    """WIN/LOSS/PUSH, or a fractional payout per contract as a string ("0.50"
    = a first-five tie): 1 / 0 / 0.5 are won / lost / push, anything else a
    partial settlement the contract has no word for (unknown, never open)."""
    if outcome_status == "PUSH":
        return "push"
    if outcome_status == "WIN":
        return "won" if is_bid else "lost"
    if outcome_status == "LOSS":
        return "lost" if is_bid else "won"
    fraction = _to_number(outcome_status)
    if fraction is None:
        return "unknown"
    paid = fraction if is_bid else 1 - fraction
    return {1.0: "won", 0.0: "lost", 0.5: "push"}.get(paid, "unknown")


def _parse_ms(value: object) -> float:
    if not isinstance(value, str):
        return float("nan")
    try:
        return datetime.fromisoformat(value.replace("Z", "+00:00")).timestamp()
    except ValueError:
        return float("nan")


def _order_status(order: dict, market: dict, outcome: dict, matched_qty: float) -> str:
    fills = order.get("fills") if isinstance(order.get("fills"), list) else []
    if order.get("status") == "REJECTED":
        return "void"
    if _to_number(order.get("qty")) == 0 and fills and all(fill.get("isWash") for fill in fills):
        return "closed"
    cash_outs = market.get("cash_out_requests")
    cash_out = cash_outs[0] if isinstance(cash_outs, list) and cash_outs else None
    if cash_out and cash_out.get("status") == "APPROVED" and _parse_ms(cash_out.get("created_at")) > _parse_ms(order.get("created_at")):
        return "closed"
    if market.get("status") == "SETTLED":
        if order.get("status") == "CANCELED" and matched_qty <= 0:
            return "void"
        return _grade_settled(outcome.get("status"), order.get("isBid") is True)
    if order.get("status") == "CANCELED":
        return "open" if matched_qty > 0 else "void"
    if order.get("status") in ("OPEN", "FILLED", "PENDING"):
        return "open"
    return "unknown"


def _parlay_status(status: object) -> str:
    upper = status.upper() if isinstance(status, str) else ""
    return {"FILLED": "open", "WIN": "won", "LOSS": "lost", "PUSH": "push", "UNFILLED": "void"}.get(upper, "unknown")


def _unmatchable(base: dict, reason: str) -> dict:
    base.update({
        "league": None, "eventStart": None, "eventDate": None, "awayTeam": None, "homeTeam": None,
        "awayKey": None, "homeKey": None, "betType": "other", "period": None, "side": None, "points": None,
        "unmatchable": reason,
    })
    return base


def _game_fields(market: dict, outcome: dict, is_bid: bool) -> dict:
    event = market.get("event") or {}
    game = event.get("game") or None
    if event.get("type") == "FUTURE" or market.get("player"):
        return {"unmatchable": REASON_NOT_GAME}
    spec = MARKET_TYPES.get(market.get("type"))
    if not spec:
        return {"unmatchable": REASON_NOT_GAME}
    if not game or not game.get("awayTeam") or not game.get("homeTeam"):
        return {"unmatchable": f"unreadable Novig order (no teams on event {event.get('id') or '?'})"}
    novig_league = game.get("league") or event.get("league") or None
    league = LEAGUES.get(novig_league)
    if not league:
        return {"unmatchable": f"league not supported ({novig_league or 'unknown'})"}
    bet_type, half = spec
    backed = _outcome_side(bet_type, outcome, market, game)
    if isinstance(backed, str):
        return {"unmatchable": f"unreadable Novig order ({backed})"}
    taken = backed if is_bid else _lay_side(bet_type, backed)
    event_start = _iso_utc(game.get("scheduled_start") or event.get("scheduled_start"))
    return {
        "league": league, "eventStart": event_start,
        "eventDate": _eastern_date_of(event_start) if event_start else None,
        "awayTeam": _team_name(game.get("awayTeam")), "homeTeam": _team_name(game.get("homeTeam")),
        "awayKey": None, "homeKey": None, "betType": bet_type,
        "period": ("F5" if league == "mlb" else "1H") if half else "FG",
        "side": taken["side"], "points": taken["points"], "unmatchable": None,
    }


def _market_title(market: dict, outcome: dict) -> str:
    return " · ".join(part for part in (_matchup_of(market), market.get("type"), outcome.get("description")) if part)


def normalize_order(order: dict, read_at: str | None) -> dict:
    market = order.get("market") or {}
    outcome = order.get("outcome") or {}
    original_qty = (_to_number(order.get("originalQty")) or 0.0) / QTY_PER_CONTRACT
    remaining_qty = (_to_number(order.get("qty")) or 0.0) / QTY_PER_CONTRACT
    matched_qty = max(original_qty - remaining_qty, 0.0)
    is_bid = order.get("isBid") is True
    outcome_price = _to_number(order.get("price"))
    taken_price = None if outcome_price is None else (outcome_price if is_bid else 1 - outcome_price)
    status = _order_status(order, market, outcome, matched_qty)
    approx = []
    if status == "open" and matched_qty == 0:
        approx.append(APPROX_PENDING if order.get("status") == "PENDING" else APPROX_UNMATCHED)
    sized_qty = original_qty if status == "open" and matched_qty == 0 else matched_qty
    fills = order.get("fills") if isinstance(order.get("fills"), list) else []
    game = (market.get("event") or {}).get("game") or None
    base = {
        "id": f"novig:{order['id']}",
        "source": SOURCE,
        "venue": VENUE,
        "rotation": None,
        "price": None if taken_price is None else probability_to_american(taken_price),
        "stake": None if taken_price is None else _round_cents(sized_qty * taken_price),
        "toWin": None if taken_price is None else _round_cents(sized_qty * (1 - taken_price)),
        "contracts": sized_qty,
        "placedAt": _iso_utc(order.get("created_at")),
        "status": status,
        "closedAt": None if status == "open" else (_iso_utc(order.get("updated_at")) or _iso_utc(order.get("created_at"))),
        "isParlayLeg": False, "parlayId": None, "legIndex": None, "legCount": None,
        "approx": approx,
        "sourceFetchedAt": read_at,
        "raw": {
            "orderId": order["id"], "orderStatus": order.get("status"), "isBid": is_bid,
            "outcomePrice": outcome_price, "originalQty": original_qty, "remainingQty": remaining_qty, "qtyUnit": QTY_PER_CONTRACT,
            "fillCount": len(fills), "marketType": market.get("type") or None, "marketStatus": market.get("status") or None,
            "strike": _to_number(market.get("strike")), "outcomeIndex": outcome.get("index"),
            "outcomeDescription": outcome.get("description") or None, "outcomeStatus": outcome.get("status") or None,
            "novigLeague": game.get("league") if game else None,
            "marketTitle": _market_title(market, outcome),
        },
    }
    fields = _game_fields(market, outcome, is_bid)
    if fields.get("unmatchable"):
        return _unmatchable(base, fields["unmatchable"])
    base.update(fields)
    return base


def normalize_parlay(parlay: dict, read_at: str | None) -> list[dict]:
    legs = parlay.get("legs") if isinstance(parlay.get("legs"), list) else []
    status = _parlay_status(parlay.get("status"))
    wager = _to_number(parlay.get("wager"))
    parlay_price = _to_number(parlay.get("price"))
    records = []
    for leg_index, leg in enumerate(legs):
        outcome = leg.get("outcome") or {}
        market = outcome.get("market") or {}
        leg_price = _to_number(leg.get("price"))
        game = (market.get("event") or {}).get("game") or None
        base = {
            "id": f"novig:{parlay['id']}:{leg_index}",
            "source": SOURCE,
            "venue": VENUE,
            "rotation": None,
            "price": None if leg_price is None else probability_to_american(leg_price),
            "stake": wager,
            "toWin": None if wager is None or not (parlay_price and parlay_price > 0) else _round_cents(wager / parlay_price - wager),
            "contracts": None,
            "placedAt": _iso_utc(parlay.get("created_at")),
            "status": status,
            "closedAt": None if status == "open" else (_iso_utc(parlay.get("updated_at")) or _iso_utc(parlay.get("created_at"))),
            "isParlayLeg": True, "parlayId": f"novig:{parlay['id']}", "legIndex": leg_index, "legCount": len(legs),
            "approx": [],
            "sourceFetchedAt": read_at,
            "raw": {
                "parlayId": parlay["id"], "parlayStatus": parlay.get("status"), "parlayPrice": parlay_price, "legPrice": leg_price,
                "marketType": market.get("type") or None, "strike": _to_number(market.get("strike")), "outcomeIndex": outcome.get("index"),
                "outcomeDescription": outcome.get("description") or None, "outcomeStatus": outcome.get("status") or None,
                "novigLeague": game.get("league") if game else None,
                "marketTitle": _market_title(market, outcome),
            },
        }
        fields = _game_fields(market, outcome, True)
        records.append(_unmatchable(base, fields["unmatchable"]) if fields.get("unmatchable") else {**base, **fields})
    return records


def normalize_novig(orders: list[dict], parlays: list[dict], read_at: str | None) -> list[dict]:
    """Every record the rows describe; orders and parlay legs dedupe on their
    native id (a row seen twice keeps the later one). Pure."""
    by_id: dict[str, dict] = {}
    for order in orders:
        if not order or not order.get("id"):
            continue
        record = normalize_order(order, read_at)
        by_id[record["id"]] = record
    for parlay in parlays:
        if not parlay or not parlay.get("id"):
            continue
        for record in normalize_parlay(parlay, read_at):
            by_id[record["id"]] = record
    return [json_clean(record) for record in by_id.values()]


# ---- network half -----------------------------------------------------------------

PAGE_LIMIT = 100
MAX_PAGES = 50  # 5,000 rows per list per poll — far above any account's live window

# The fields the app's Portfolio cards read (OrderStateCard_Frag +
# OrderDescription_Frag + EventData_Frag), spelled out so the row shape is
# visible here and matches tests/fixtures/bets/novig_bets.json.
TEAM_FIELDS = "id name symbol short_name"
EVENT_FIELDS = f"""id type status scheduled_start league description
  game {{ id league sport scheduled_start home_score away_score
    awayTeam {{ {TEAM_FIELDS} }} homeTeam {{ {TEAM_FIELDS} }} }}"""
ORDER_FIELDS = f"""id timeToLiveMs qty originalQty created_at updated_at status price isBid tif currency
  fills {{ id created_at cost qty isWash isTaker }}
  market {{ id type strike status re_settled_at
    cash_out_requests(limit: 1, order_by: {{ created_at: desc }}) {{ id created_at status }}
    player {{ id full_name }} competitor {{ id name symbol }}
    event {{ {EVENT_FIELDS} }} }}
  outcome {{ id index description status competitor {{ id name symbol }} }}"""
PARLAY_FIELDS = f"""id created_at updated_at status wager price unboosted_price
  wallet {{ currency }}
  legs {{ price outcome {{ id index description status competitor {{ id name symbol }}
    market {{ id type strike player {{ id }} event {{ {EVENT_FIELDS} }} }} }} }}"""

USER_QUERY = """query BetsService_User($auth_id: String!) {
  user(where: { auth_id: { _eq: $auth_id } }) { id trader_id }
}"""
ORDERS_QUERY = f"""query BetsService_Orders($trader_id: uuid!, $since: timestamptz!, $limit: Int!, $offset: Int!) {{
  order(where: {{ trader: {{ id: {{ _eq: $trader_id }} }},
                  _or: [{{ market: {{ status: {{ _neq: "SETTLED" }} }} }}, {{ updated_at: {{ _gte: $since }} }}] }},
        order_by: {{ created_at: desc }}, limit: $limit, offset: $offset) {{ {ORDER_FIELDS} }}
}}"""
PARLAYS_QUERY = f"""query BetsService_Parlays($trader_id: uuid!, $since: timestamptz!, $limit: Int!, $offset: Int!) {{
  parlay(where: {{ trader: {{ id: {{ _eq: $trader_id }} }},
                   _or: [{{ status: {{ _eq: "FILLED" }} }}, {{ updated_at: {{ _gte: $since }} }}] }},
         order_by: {{ created_at: desc }}, limit: $limit, offset: $offset) {{ {PARLAY_FIELDS} }}
}}"""

GraphqlCall = Callable[[str, dict, str], dict]


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


class NovigSource:
    """Source protocol implementation for a Novig account (read-only queries).

    `graphql(query, variables, bearer) -> data` is injectable for tests; the
    default posts through curl_cffi. `auth` is a NovigAuth over the token file.
    """

    name = "novig"

    def __init__(self, auth: NovigAuth | None = None, graphql: GraphqlCall | None = None,
                 poll_sec: float | None = None, retention_days: int | None = None,
                 clock: Callable[[], float] = time.time):
        self.poll_sec = poll_sec if poll_sec is not None else config.NOVIG_POLL_SEC
        self._retention_days = retention_days if retention_days is not None else config.RETENTION_DAYS
        self._auth = auth or NovigAuth(config.NOVIG_TOKEN_PATH)
        self._graphql = graphql or self._post_graphql
        self._clock = clock
        self._session = None
        self._trader_id: str | None = None

    def _post_graphql(self, query: str, variables: dict, bearer: str) -> dict:
        if self._session is None:
            self._session = _make_session()
        response = self._session.post(config.NOVIG_GRAPHQL_URL, json={"query": query, "variables": variables},
                                      headers={"Authorization": f"Bearer {bearer}", "Content-Type": "application/json"},
                                      timeout=20)
        if response.status_code != 200:
            raise RuntimeError(f"novig graphql HTTP {response.status_code}: {response.text[:200]}")
        body = response.json()
        if body.get("errors"):
            raise RuntimeError(f"novig graphql error: {body['errors'][0].get('message', body['errors'][0])}")
        return body.get("data") or {}

    def _resolve_trader_id(self, bearer: str, auth_id: str) -> str:
        if self._trader_id:
            return self._trader_id
        users = self._graphql(USER_QUERY, {"auth_id": auth_id}, bearer).get("user") or []
        if len(users) != 1 or not users[0].get("trader_id"):
            raise RuntimeError(f"expected exactly one Novig user with a trader_id for {auth_id}, got {len(users)}")
        self._trader_id = users[0]["trader_id"]
        return self._trader_id

    def _paginate(self, query: str, key: str, variables: dict, bearer: str) -> list[dict]:
        rows: list[dict] = []
        for page in range(MAX_PAGES):
            data = self._graphql(query, {**variables, "limit": PAGE_LIMIT, "offset": page * PAGE_LIMIT}, bearer)
            batch = data.get(key)
            if not isinstance(batch, list):
                raise RuntimeError(f"novig graphql: expected a `{key}` list, got {type(batch).__name__}")
            rows.extend(batch)
            if len(batch) < PAGE_LIMIT:
                return rows
        raise RuntimeError(f"novig `{key}` did not end within {MAX_PAGES} pages; refusing a partial pull")

    def fetch(self) -> list[dict]:
        try:
            access = self._auth.token()
        except NovigAuthError as error:
            raise RuntimeError(f"novig auth: {error}") from error
        trader_id = self._resolve_trader_id(access.token, access.auth_id)
        since = (datetime.now(timezone.utc) - timedelta(days=self._retention_days)).strftime("%Y-%m-%dT%H:%M:%SZ")
        variables = {"trader_id": trader_id, "since": since}
        orders = self._paginate(ORDERS_QUERY, "order", variables, access.token)
        parlays = self._paginate(PARLAYS_QUERY, "parlay", variables, access.token)
        records = normalize_novig(orders, parlays, utc_now_iso())
        log.info("novig: %d orders, %d parlays -> %d records", len(orders), len(parlays), len(records))
        return records


def source_if_connected(token_path: Path | None = None) -> NovigSource | None:
    """The source when a token file exists, else None (logged) — an unconnected
    Novig must read "no source configured", not a permanently failing row."""
    path = token_path or config.NOVIG_TOKEN_PATH
    if not path.exists():
        log.info("novig: no token at %s — not registered (run novig_auth connect)", path)
        return None
    return NovigSource()
