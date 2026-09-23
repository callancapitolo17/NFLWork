"""Wagerzon (the C account) bet source: HistoryHelper weeks + the open-bets helper -> normalised bet records.

Two halves, like sources/betonline.py:
  normalize_wagerzon(wagers, fetched_at)   pure parser of HistoryHelper wager rows; the pytest on
                                           tests/fixtures/bets/wagerzon_history.json pins it.
  WagerzonSource                           the network half (Source protocol): an ASP.NET form
                                           login whose session cookie lives in memory, then
                                           `HistoryHelper.aspx?week=N` for the last HISTORY_WEEKS
                                           Mon-Sun weeks plus `OpenBetsHelper.aspx`, every wager
                                           normalised (pending bets kept).

Inputs:  WAGERZONC_USERNAME / WAGERZONC_PASSWORD, else WAGERZON_USERNAME / WAGERZON_PASSWORD
         (config: the environment, bets_service/.env, kalshi_draft/.env, then bet_logger/.env in
         the main checkout — the C account's login has sat in the primary slot since
         2026-06-26, bet_logger/scraper_wagerzon.py);
         GET backend.wagerzon.com/wager/HistoryHelper.aspx?week=N   {result: {details: [{Date,
             wager: [...]}], StartDate, EndDate, ErrorMsg, ...}}
         GET backend.wagerzon.com/wager/OpenBetsHelper.aspx          {result: [...]} (the helper
             the OpenBets.aspx widget loads; `{"result": []}` with no open bet on 2026-09-22, so a
             row's shape is unobserved — a row that is not a HistoryHelper-shaped wager fails the
             poll loudly, naming its keys, rather than being dropped).
Outputs: list of records. Side effects: none on disk. Raises on a failed login, a helper that
         answers HTML or a redirect after a fresh login (the session died), an error message in a
         week's body, or an open-bets row of unknown shape, so the service records a failed run
         and keeps the previous records (never a partial list).

Why not import bet_logger/scraper_wagerzon.py: it imports Google Sheets and dotenv at module
level, skips pending bets and prints. The login, endpoint and grammar are copied from it and
must stay in step.

Wager grammar (HistoryHelper live 2026-09-22 — props only that week — plus every
bet_logger run since 2026-04 in bet_logger/logs; the parser fails closed on anything else,
listing the record as unmatchable with the reason and the raw description):
  IdWager        219337542 — the stable native id (TicketNumber carries the same value; 0 on a
                 transaction row). WagerOrTrans "WAGER" | "TRAN" (a transfer, not a bet: skipped).
  Result         "" = open; WIN / LOSE / PUSH / CANCELLED / VOID / NO ACTION / NO BET.
  RiskAmount, WinAmount, WinLoss   USD as strings.
  PlacedDate + PlacedTime          "09/20/2026" + "04:09 PM " on the site's Eastern clock.
  HeaderDesc     "STRAIGHT BET" | "PARLAY (2 TEAMS)" | "NFL WEEK 2 - SPECIALS" | "DST Straight|ID:…"
                 — informational; the leg count is len(details).
  details[]      one per leg: DetailDesc, DetailResult, GameDate + GameTime ("09/21/2026",
                 "08:15 PM", Eastern — the event start), IdSport (NFL / CFB / NBA / CBK / WNBA /
                 MLB / NHL / SOC name the league; PROP / RBL / DST are props; XFER a transfer).
  DetailDesc     "[967] TOTAL o7½-120 (CHI CUBS vrs TB RAYS)<BR>( JAMESON TAILLON - R / SHANE MCCLANAHAN - L )"
                     total: rotation, o/u + number, price (EV = +100), "(AWAY vrs HOME)" away first,
                     then the pitchers bracket (ignored)
                 "[969] 1H ARI DBACKS -½+115<BR>( … pitchers … )"   spread: own team only, the period
                     token BEFORE the name, no opponent anywhere — placed by rotation parity
                 "TOTAL o5EV (1H HOU ASTROS GM#2 vrs 1H BAL ORIOLES GM#2)"   a doubleheader's game
                     number after the name (kept in raw.gameNumber)
                 "[212] TEXAS TECH +168"                             moneyline (scraper_wagerzon.py)
                 "( NYM vs COL Has Been Postponed. NO Action )"      a postponed leg: fails closed
                 "[777045] CARDINALS SUPERFECTA (…) +1770"           a prop: IdSport says so first
                 ½ ¼ ¾ fractions; HTML tags separate the brackets.
  MLB "1H" is the first five innings (Novig's rule in sources/novig.py): period F5.
eventStart is the leg's own GameDate + GameTime (Eastern -> UTC) and eventDate its GameDate, so
the matcher uses the 30-minute time rule; a settled bet's closedAt is its event start (no settle
time is served; the placed time is the fallback). A spread or moneyline names one team, placed
by rotation parity (odd = away, approx side_from_rotation_parity) unless a "(AWAY vrs HOME)"
bracket follows it (not seen live; honoured when present).
"""
import logging
import re
import threading
from collections.abc import Callable
from datetime import datetime, timezone

import requests

from unabated_ticket.bets_service import config
from unabated_ticket.bets_service.normalize import EASTERN, json_clean, round_cents, utc_now_iso
from unabated_ticket.bets_service.sources.betonline import (
    APPROX_SIDE_PARITY, american_from_payout, parse_points, side_from_rotation)
from unabated_ticket.bets_service.sources.bfa import NUMBER, PRICE, VERSUS_RE, iso_utc, parse_price

log = logging.getLogger(__name__)

# Login and endpoints: copied from bet_logger/scraper_wagerzon.py (keep in step).
BASE_URL = "https://backend.wagerzon.com"
HISTORY_URL = f"{BASE_URL}/wager/HistoryHelper.aspx"
OPEN_BETS_URL = f"{BASE_URL}/wager/OpenBetsHelper.aspx"
USER_AGENT = ("Mozilla/5.0 (Macintosh; Intel Mac OS X 10_15_7) "
              "AppleWebKit/537.36 (KHTML, like Gecko) Chrome/145.0.0.0 Safari/537.36")
XHR_HEADERS = {"X-Requested-With": "XMLHttpRequest", "Accept": "application/json"}
# Posted back exactly as served: ASP.NET validates the form against them.
ASPNET_HIDDEN_FIELDS = ("__VIEWSTATE", "__VIEWSTATEGENERATOR", "__EVENTVALIDATION",
                        "__EVENTTARGET", "__EVENTARGUMENT")
HTTP_TIMEOUT_SEC = 20

SOURCE = "wagerzon_api"
VENUE = "wagerzon"
# The site renders every timestamp in Eastern (Monday Night Football reads 08:15 PM).
WAGERZON_TZ = EASTERN
WAGER_ROW_MARKER = "WAGER"
LEAGUES = {"NFL": "nfl", "CFB": "cfb", "NBA": "nba", "CBK": "cbb", "WNBA": "wnba", "MLB": "mlb",
           "NHL": "nhl", "SOC": "soccer"}
PROP_SPORT_CODES = {"PROP", "RBL", "DST"}
STATUS_MAP = {"": "open", "win": "won", "lose": "lost", "push": "push", "cancelled": "void",
              "canceled": "void", "void": "void", "no action": "void", "no bet": "void"}
PERIODS = {"1H": "1H", "2H": "2H", "1Q": "1Q", "2Q": "2Q", "3Q": "3Q", "4Q": "4Q"}
POSTPONED_MARKER = "POSTPONED"
REASON_NO_LEGS = "wager carries no legs"

HTML_TAG_RE = re.compile(r"<[^>]+>")
ROTATION_RE = re.compile(r"^\[(?P<rotation>\d+)\]\s*(?P<body>.*)$")
# Zero or more brackets after the price: the pitchers, or "(AWAY vrs HOME)".
BRACKETS = r"(?P<brackets>(?:\s*\([^()]*\))*)"
TEAM = r"(?P<team>[^()]+?)"
TOTAL_RE = re.compile(rf"^TOTAL\s+(?P<direction>[ou])(?P<points>{NUMBER})(?P<price>{PRICE})"
                      rf"\s*\((?P<context>[^()]*)\){BRACKETS}\s*$", re.IGNORECASE)
SPREAD_RE = re.compile(rf"^{TEAM}\s+(?P<points>[+-]{NUMBER}|PK)(?P<price>{PRICE}){BRACKETS}\s*$", re.IGNORECASE)
MONEYLINE_RE = re.compile(rf"^{TEAM}\s+(?P<price>{PRICE}){BRACKETS}\s*$", re.IGNORECASE)
BRACKET_RE = re.compile(r"\(([^()]*)\)")
# "1H ARI DBACKS GM#2": the period token before the name, the doubleheader game after it.
TEAM_TOKENS_RE = re.compile(r"^(?:(?P<period>\d[HQ])\s+)?(?P<name>.+?)(?:\s+GM#(?P<game>\d))?$", re.IGNORECASE)


# ---- pure parser --------------------------------------------------------------------

def clean_description(raw: str) -> str:
    """HTML tags become spaces, whitespace collapses: one line to match against."""
    return re.sub(r"\s+", " ", HTML_TAG_RE.sub(" ", raw)).strip()


def parse_team(text: str) -> dict | str:
    """"1H ARI DBACKS GM#2" -> {name, period, gameNumber}; an unknown period token -> reason."""
    match = TEAM_TOKENS_RE.match(text.strip())
    if not match:
        return f"no team name ({text.strip()[:60]})"
    period = None
    if match.group("period"):
        period = PERIODS.get(match.group("period").upper())
        if period is None:
            return f"unknown period ({match.group('period')})"
    game = match.group("game")
    return {"name": match.group("name").strip(), "period": period, "gameNumber": int(game) if game else None}


def parse_team_pair(context: str) -> dict | str:
    """"1H HOU ASTROS GM#2 vrs 1H BAL ORIOLES GM#2" -> {awayTeam, homeTeam, period, gameNumber}."""
    parts = VERSUS_RE.split(context.strip())
    if len(parts) != 2:
        return f"no 'AWAY vrs HOME' bracket ({context.strip()[:60]})"
    away, home = parse_team(parts[0]), parse_team(parts[1])
    if isinstance(away, str):
        return away
    if isinstance(home, str):
        return home
    periods = {team["period"] for team in (away, home) if team["period"]}
    if len(periods) > 1:
        return f"periods differ ({away['period']} vs {home['period']})"
    games = {team["gameNumber"] for team in (away, home) if team["gameNumber"]}
    if len(games) > 1:
        return f"doubleheader games differ ({away['gameNumber']} vs {home['gameNumber']})"
    return {"awayTeam": away["name"], "homeTeam": home["name"],
            "period": periods.pop() if periods else None, "gameNumber": games.pop() if games else None}


def teams_bracket_of(brackets: str) -> str | None:
    """The "(AWAY vrs HOME)" bracket among those after a price, if any (the pitchers
    bracket names no teams)."""
    for content in BRACKET_RE.findall(brackets or ""):
        if VERSUS_RE.search(content):
            return content
    return None


def _one_team_leg(rotation: int | None, bet_type: str, team_text: str, points: float | None,
                  price: int, brackets: str) -> dict | str:
    team = parse_team(team_text)
    if isinstance(team, str):
        return team
    leg = {"rotation": rotation, "betType": bet_type, "points": points, "price": price,
           "period": team["period"], "gameNumber": team["gameNumber"]}
    context = teams_bracket_of(brackets)
    if context is None:
        if rotation is None:
            return f"no rotation to place the side of {team['name']}"
        side = side_from_rotation(rotation)
        away, home = (team["name"], None) if side == "away" else (None, team["name"])
        leg.update(side=side, awayTeam=away, homeTeam=home, sideFromParity=True)
        return leg
    pair = parse_team_pair(context)
    if isinstance(pair, str):
        return pair
    if team["name"].upper() == pair["awayTeam"].upper():
        side = "away"
    elif team["name"].upper() == pair["homeTeam"].upper():
        side = "home"
    else:
        return f"team {team['name']} is not in its bracket ({context.strip()[:60]})"
    leg.update(side=side, awayTeam=pair["awayTeam"], homeTeam=pair["homeTeam"],
               period=team["period"] or pair["period"], gameNumber=team["gameNumber"] or pair["gameNumber"],
               sideFromParity=False)
    return leg


def parse_leg(raw_description: str) -> dict | str:
    """A DetailDesc -> {rotation, betType, side, points, price, period (None = full game),
    awayTeam, homeTeam, gameNumber, sideFromParity}, or the reason it does not parse."""
    text = clean_description(raw_description)
    rotation = None
    body = text
    rotation_match = ROTATION_RE.match(text)
    if rotation_match:
        rotation = int(rotation_match.group("rotation"))
        body = rotation_match.group("body").strip()
    if POSTPONED_MARKER in body.upper():
        return f"postponed ({body[:60]})"
    total = TOTAL_RE.match(body)
    if total:
        pair = parse_team_pair(total.group("context"))
        if isinstance(pair, str):
            return pair
        return {"rotation": rotation, "betType": "total",
                "side": "over" if total.group("direction").lower() == "o" else "under",
                "points": parse_points(total.group("points")), "price": parse_price(total.group("price")),
                "period": pair["period"], "awayTeam": pair["awayTeam"], "homeTeam": pair["homeTeam"],
                "gameNumber": pair["gameNumber"], "sideFromParity": False}
    spread = SPREAD_RE.match(body)
    if spread:
        return _one_team_leg(rotation, "spread", spread.group("team"), parse_points(spread.group("points")),
                             parse_price(spread.group("price")), spread.group("brackets"))
    moneyline = MONEYLINE_RE.match(body)
    if moneyline:
        return _one_team_leg(rotation, "moneyline", moneyline.group("team"), None,
                             parse_price(moneyline.group("price")), moneyline.group("brackets"))
    return f"unrecognised selection ({body[:60]})"


def league_of(sport_code: object) -> str | None:
    return LEAGUES.get(str(sport_code or "").strip().upper())


def period_for(league: str, period: str | None) -> str:
    """The leg's period; Wagerzon's MLB "1H" is the first five innings."""
    if period is None:
        return "FG"
    if period == "1H" and league == "mlb":
        return "F5"
    return period


def parse_eastern(date_text: object, time_text: object) -> datetime | None:
    """"09/21/2026" + "08:15 PM " on the site's Eastern clock -> aware UTC; blank -> None."""
    if not isinstance(date_text, str) or not isinstance(time_text, str):
        return None
    try:
        naive = datetime.strptime(f"{date_text.strip()} {time_text.strip()}", "%m/%d/%Y %I:%M %p")
    except ValueError:
        return None
    return naive.replace(tzinfo=WAGERZON_TZ).astimezone(timezone.utc)


def eastern_date_of(date_text: object) -> str | None:
    if not isinstance(date_text, str):
        return None
    try:
        return datetime.strptime(date_text.strip(), "%m/%d/%Y").strftime("%Y-%m-%d")
    except ValueError:
        return None


def status_of(result: object) -> str:
    return STATUS_MAP.get(str(result or "").strip().lower(), "unknown")


def native_id_of(wager: dict) -> str:
    value = wager.get("IdWager")
    if value in (None, "", 0, "0"):
        raise RuntimeError(f"Wagerzon wager carries no IdWager; keys seen: {sorted(wager.keys())}")
    return str(value)


def _money(value: object) -> float | None:
    if value is None or isinstance(value, bool):
        return None
    try:
        return round_cents(float(value))
    except (TypeError, ValueError):
        return None


def _base_record(wager: dict, native_id: str, fetched_at: str | None) -> dict:
    status = status_of(wager.get("Result"))
    placed_at = parse_eastern(wager.get("PlacedDate"), wager.get("PlacedTime"))
    return {
        "id": f"{VENUE}:{native_id}",
        "source": SOURCE,
        "venue": VENUE,
        "league": None, "eventStart": None, "eventDate": None,
        "awayTeam": None, "homeTeam": None, "awayKey": None, "homeKey": None,
        "rotation": None, "betType": "other", "period": None, "side": None, "points": None,
        "price": None,
        "stake": _money(wager.get("RiskAmount")),
        "toWin": _money(wager.get("WinAmount")),
        "contracts": None,
        "placedAt": iso_utc(placed_at),
        "status": status,
        # No settle time is served: an open bet has none, a settled one closes at
        # its event start once the leg is read (_apply_leg), else its placed time.
        "closedAt": None if status == "open" else iso_utc(placed_at),
        "isParlayLeg": False, "parlayId": None, "legIndex": None, "legCount": None,
        "approx": [],
        "unmatchable": None,
        "sourceFetchedAt": fetched_at,
        "raw": {
            "nativeId": native_id,
            "headerDesc": wager.get("HeaderDesc"),
            "result": wager.get("Result"),
            "riskAmount": wager.get("RiskAmount"),
            "winAmount": wager.get("WinAmount"),
            "winLoss": wager.get("WinLoss"),
            "placedDate": wager.get("PlacedDate"),
            "placedTime": wager.get("PlacedTime"),
            "wagerType": wager.get("WagerType"),
            "ifBetWagerType": wager.get("IfBetWagerType"),
        },
    }


def _describe_leg(record: dict, leg_row: dict) -> dict:
    record["raw"].update({
        "description": leg_row.get("DetailDesc"),
        "legResult": leg_row.get("DetailResult"),
        "sportCode": leg_row.get("IdSport"),
        "gameDate": leg_row.get("GameDate"),
        "gameTime": leg_row.get("GameTime"),
        "idGame": leg_row.get("IdGame"),
    })
    return record


def _apply_leg(record: dict, league: str, leg: dict, leg_row: dict) -> dict:
    event_start = parse_eastern(leg_row.get("GameDate"), leg_row.get("GameTime"))
    record.update({
        "league": league,
        "eventStart": iso_utc(event_start),
        "eventDate": eastern_date_of(leg_row.get("GameDate")),
        "rotation": leg["rotation"], "betType": leg["betType"], "period": period_for(league, leg["period"]),
        "side": leg["side"], "points": leg["points"], "price": leg["price"],
        "awayTeam": leg["awayTeam"], "homeTeam": leg["homeTeam"],
        "approx": [APPROX_SIDE_PARITY] if leg["sideFromParity"] else [],
    })
    if record["status"] != "open" and event_start is not None:
        record["closedAt"] = iso_utc(event_start)
    record["raw"]["gameNumber"] = leg["gameNumber"]
    return record


def _unmatchable(record: dict, reason: str) -> dict:
    record["unmatchable"] = reason
    return record


def _leg_record(record: dict, leg_row: dict) -> dict:
    """The record for one leg row: league from IdSport, then the description."""
    _describe_leg(record, leg_row)
    sport_code = str(leg_row.get("IdSport") or "").strip().upper()
    if sport_code in PROP_SPORT_CODES:
        return _unmatchable(record, f"not a game market (IdSport {sport_code})")
    league = league_of(sport_code)
    if league is None:
        return _unmatchable(record, f"league not supported (IdSport {sport_code or 'blank'})")
    leg = parse_leg(str(leg_row.get("DetailDesc") or ""))
    if isinstance(leg, str):
        return _unmatchable(record, leg)
    return _apply_leg(record, league, leg, leg_row)


def normalize_wager(wager: dict, fetched_at: str | None) -> list[dict]:
    """One HistoryHelper wager -> one record, or one per leg of a parlay; a transaction
    row -> nothing. Pure."""
    if wager.get("WagerOrTrans") != WAGER_ROW_MARKER:
        return []
    native_id = native_id_of(wager)
    legs = wager.get("details") or []
    base = _base_record(wager, native_id, fetched_at)
    if not legs:
        return [json_clean(_unmatchable(base, REASON_NO_LEGS))]
    if len(legs) == 1:
        return [json_clean(_leg_record(base, legs[0]))]
    parlay_price = american_from_payout(base["stake"] or 0, base["toWin"] or 0)
    records = []
    for index, leg_row in enumerate(legs):
        record = _base_record(wager, native_id, fetched_at)
        record["id"] = f"{base['id']}:leg{index}"
        record["raw"]["parlayPrice"] = parlay_price
        record.update({"isParlayLeg": True, "parlayId": base["id"], "legIndex": index, "legCount": len(legs)})
        records.append(json_clean(_leg_record(record, leg_row)))
    return records


def normalize_wagerzon(wagers: list[dict], fetched_at: str | None) -> list[dict]:
    """Every wager row -> records (pending bets KEPT; a parlay is one record per leg;
    transaction rows dropped). Raises only when a wager has no IdWager — the one
    thing the store cannot work without."""
    records: list[dict] = []
    for wager in wagers:
        records.extend(normalize_wager(wager, fetched_at))
    return records


def wagers_of_history(result: dict) -> list[dict]:
    """The wager rows of one HistoryHelper week, in the order served (a day at a time)."""
    return [wager for day in result.get("details") or [] for wager in day.get("wager") or []]


def open_wagers_of(body: object) -> list[dict]:
    """The OpenBetsHelper rows when they are HistoryHelper-shaped wagers; an unknown
    shape fails loudly with its keys (the shape is unobserved — module docstring)."""
    rows = body.get("result") if isinstance(body, dict) else None
    if not isinstance(rows, list):
        keys = sorted(body.keys()) if isinstance(body, dict) else type(body).__name__
        raise RuntimeError(f"Wagerzon OpenBetsHelper: expected {{result: [...]}}, got {keys}")
    for row in rows:
        if not isinstance(row, dict) or "IdWager" not in row or not isinstance(row.get("details"), list):
            keys = sorted(row.keys()) if isinstance(row, dict) else type(row).__name__
            raise RuntimeError(f"Wagerzon OpenBetsHelper row in an unknown shape (keys {keys}); "
                               "capture it and extend sources/wagerzon.py")
    return rows


# ---- network half -------------------------------------------------------------------

def new_session() -> requests.Session:
    session = requests.Session()
    session.headers["User-Agent"] = USER_AGENT
    return session


def json_body_of(response: object) -> object | None:
    """The response's JSON, or None when the site answered anything else (a redirect
    to the login page, an HTML page): the session is dead."""
    if response.status_code != 200:
        return None
    if "json" not in response.headers.get("Content-Type", "").lower():
        return None
    try:
        return response.json()
    except ValueError:
        return None


class WagerzonSource:
    """Source protocol implementation for the Wagerzon C account (helper GETs only).

    `session_factory()` builds an HTTP session with .get()/.post() and a cookie jar
    (requests by default; tests inject a fake). The logged-in session is kept
    across polls and replaced when a helper stops answering JSON.
    """

    name = VENUE

    def __init__(self, username: str | None = None, password: str | None = None,
                 history_weeks: int | None = None, poll_sec: float | None = None,
                 session_factory: Callable[[], object] = new_session):
        self.poll_sec = poll_sec if poll_sec is not None else config.WAGERZON_POLL_SEC
        self._username = username if username is not None else config.WAGERZON_USERNAME
        self._password = password if password is not None else config.WAGERZON_PASSWORD
        self._history_weeks = history_weeks if history_weeks is not None else config.WAGERZON_HISTORY_WEEKS
        self._session_factory = session_factory
        self._session: object | None = None
        self._session_lock = threading.Lock()

    # -- login ------------------------------------------------------------------------

    def _login(self) -> object:
        """ASP.NET form login: the hidden fields are posted back exactly as served."""
        if not self._username or not self._password:
            raise RuntimeError("Wagerzon credentials missing: set WAGERZONC_USERNAME and WAGERZONC_PASSWORD "
                               "(or WAGERZON_USERNAME / WAGERZON_PASSWORD) in bet_logger/.env")
        session = self._session_factory()
        login_page = session.get(BASE_URL, timeout=HTTP_TIMEOUT_SEC)
        if login_page.status_code != 200:
            raise RuntimeError(f"Wagerzon login page: HTTP {login_page.status_code}")
        fields = {}
        for name in ASPNET_HIDDEN_FIELDS:
            match = re.search(rf'(?:name|id)="{name}"[^>]*value="([^"]*)"', login_page.text)
            if match:
                fields[name] = match.group(1)
        if "__VIEWSTATE" not in fields:
            raise RuntimeError("Wagerzon login page carries no __VIEWSTATE — its form changed")
        fields.update({"Account": self._username, "Password": self._password, "BtnSubmit": ""})
        submitted = session.post(BASE_URL, data=fields, timeout=HTTP_TIMEOUT_SEC)
        if submitted.status_code != 200:
            raise RuntimeError(f"Wagerzon login failed: HTTP {submitted.status_code}")
        log.info("wagerzon: logged in")
        return session

    def _get_json(self, url: str, params: dict | None) -> object:
        """A helper's JSON body; a dead session (HTML or a redirect back) is replaced by
        one fresh login, after which a second non-JSON answer is the poll's failure."""
        with self._session_lock:
            for attempt in (1, 2):
                if self._session is None:
                    self._session = self._login()
                response = self._session.get(url, params=params, headers=XHR_HEADERS,
                                             timeout=HTTP_TIMEOUT_SEC, allow_redirects=False)
                body = json_body_of(response)
                if body is not None:
                    return body
                if attempt == 1:
                    log.info("wagerzon: %s answered HTTP %d without JSON; logging in again",
                             url.rsplit("/", 1)[-1], response.status_code)
                    self._session = None
            raise RuntimeError(f"Wagerzon {url.rsplit('/', 1)[-1]} returned no JSON after a fresh login "
                               f"(HTTP {response.status_code}); the login may be refused")

    # -- history ----------------------------------------------------------------------

    def _fetch_history(self) -> list[dict]:
        wagers: list[dict] = []
        for week in range(self._history_weeks):
            body = self._get_json(HISTORY_URL, {"week": week})
            result = body.get("result") if isinstance(body, dict) else None
            if not isinstance(result, dict):
                raise RuntimeError(f"Wagerzon HistoryHelper week {week}: no result object in the body")
            if result.get("ErrorMsg"):
                raise RuntimeError(f"Wagerzon HistoryHelper week {week}: {result['ErrorMsg']}")
            wagers.extend(wagers_of_history(result))
        return wagers

    def _fetch_open_bets(self) -> list[dict]:
        return open_wagers_of(self._get_json(OPEN_BETS_URL, None))

    # -- Source protocol --------------------------------------------------------------

    def fetch(self) -> list[dict]:
        history = self._fetch_history()
        seen = {str(wager.get("IdWager")) for wager in history}
        # The history already lists a pending bet under its game day (Result ""); the
        # open-bets helper adds only what the weeks did not carry.
        open_only = [wager for wager in self._fetch_open_bets() if str(wager.get("IdWager")) not in seen]
        records = normalize_wagerzon(history + open_only, utc_now_iso())
        n_open = sum(1 for record in records if record["status"] == "open")
        n_unmatchable = sum(1 for record in records if record["unmatchable"])
        log.info("wagerzon: %d history rows + %d open-only rows -> %d records (%d open, %d unmatchable)",
                 len(history), len(open_only), len(records), n_open, n_unmatchable)
        return records


def source_if_configured() -> WagerzonSource | None:
    """The source when the login is configured (the BetOnline pattern); missing
    credentials are logged once with the fix."""
    if not (config.WAGERZON_USERNAME and config.WAGERZON_PASSWORD):
        log.warning("wagerzon: no WAGERZONC_* or WAGERZON_* login (bet_logger/.env in the main checkout) "
                    "— Wagerzon source not registered")
        return None
    return WagerzonSource()
