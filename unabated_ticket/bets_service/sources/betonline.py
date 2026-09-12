"""BetOnline bet source: the account's bet-history report -> normalised bet records (#115).

Two halves, like sources/kalshi.py:
  normalize_betonline(rows, fetched_at)   pure parser of the report's rows; the
                                          pytest on tests/fixtures/bets/betonline_history.json
                                          pins it.
  BetOnlineSource                         the network half (Source protocol): every
                                          poll refreshes the Keycloak access token
                                          only when it is about to expire, pulls the
                                          paged report for the retention window and
                                          normalises every row (pending bets kept —
                                          the sheet scraper skips them).

Inputs:  the Keycloak refresh token (`krefresh` cookie) in
         bet_logger/recon_betonline_cookies.json of the MAIN checkout (written by
         bet_logger/recon_betonline.py; config.BETONLINE_COOKIES_PATH), the
         Cloudflare cookies next to it, and BetOnline's report endpoint
         `POST api.betonline.ag/report/api/report/get-bet-history`.
Outputs: list of records. Side effects: the ROTATED refresh token is written back
         to the cookie file atomically (temp file + os.replace) — that file is shared
         with bet_logger/scraper_betonline.py and its LaunchAgent, see "Token
         rotation" below. In-memory access-token cache only otherwise. Raises on a
         failed token refresh or any non-200 report page so the service records a
         failed run and keeps the previous records (never a partial list).

Why not import bet_logger/scraper_betonline.py: it imports Google Sheets at module
level (not in this service's venv), writes the cookie file non-atomically and its
fetch_bet_history() returns a PARTIAL list on a failed page, which the Source
protocol forbids. The endpoint, header set and cookie-file layout are copied from it
and must stay in step.

Token rotation (the race the issue names): every refresh rotates `krefresh`, and
two processes refreshing the same parent token trip Keycloak's reuse detection.
Rules: (1) this source refreshes only when its cached access token is within
REFRESH_MARGIN_SEC of expiry; (2) read -> refresh -> write runs under an exclusive
flock on `<cookie file>.lock`, which bet_logger/scraper_betonline.py takes too, so
the loser of a race waits, then re-reads and rotates the WINNER's token instead of
the dead parent; (3) it RE-READS the cookie file inside the lock, so a rotation the
LaunchAgent made while the service was down is picked up; (4) it writes the
rotated token with os.replace. The LaunchAgent stays loaded and unchanged — it is
the keep-alive for when the service is not running.

Report grammar (from the issue and bet_logger/scraper_betonline.py's regexes; no
saved payload exists, so the first live pull must confirm it — the parser fails
closed on anything it does not recognise, listing the row as unmatchable with the
reason and the raw description):
  Description  "[Desktop|Mobile] - NFL - 465 Chicago Bears +3.5 -110[ - 1st Half]"
               <league> - <rotation> <team> <points> <price>, totals as
               "<team> over 44½ -110", moneyline "<team> +140"; a trailing
               " - <period>" names a non-full-game period; a trailing " for ..."
               tail is ignored. Same Game Parlay descriptions repeat the league
               ("NFL - NFL - ") and list the legs on one line; the separator is
               unverified — ", " / " | " / "; " are tried and EVERY piece must parse,
               else the ticket is one unmatchable record.
  WagerType    Spread | Total | Money Line | Same Game Parlay (the authority for a
               straight bet's type; a description that parses as another type is
               unmatchable).
  WagerStatus  Pending | Won | Lost | Push | Cancelled.
  Risk, ToWin  USD.  Date  the placed time (ISO; a naive value is read as UTC).
The report carries NO game date and names only the bet's own team, so
eventStart / eventDate are null and the side comes from rotation parity (odd =
away, even = home, the US rotation convention) — both recorded in `approx`
("game_date_unknown", "side_from_rotation_parity"); the matcher keys on the
rotation number. The report carries no settle time either: a settled bet's
closedAt is its placed time (a lower bound; it only decides when the bet leaves
the 30-day window, never a match).
"""
import fcntl
import json
import logging
import os
import re
import threading
import time
from collections.abc import Callable, Iterator
from contextlib import contextmanager
from datetime import datetime, timedelta, timezone
from pathlib import Path

import requests

from unabated_ticket.bets_service import config
from unabated_ticket.bets_service.normalize import json_clean, round_cents, utc_now_iso

log = logging.getLogger(__name__)

# Endpoints and headers: copied from bet_logger/scraper_betonline.py (keep in step).
TOKEN_URL = "https://api.betonline.ag/api/auth/realms/betonline/protocol/openid-connect/token"
BET_HISTORY_URL = "https://api.betonline.ag/report/api/report/get-bet-history"
CLIENT_ID = "betonline-web"
BASE_HEADERS = {
    "Accept": "application/json, text/plain, */*",
    "Origin": "https://www.betonline.ag",
    "Referer": "https://www.betonline.ag/",
    "User-Agent": (
        "Mozilla/5.0 (Macintosh; Intel Mac OS X 10_15_7) "
        "AppleWebKit/537.36 (KHTML, like Gecko) Chrome/131.0.0.0 Safari/537.36"
    ),
    "gsetting": "bolnasite",
    "contests": "na",
    "gmt-offset": "-8",
    "utc-offset": "480",
}
COOKIE_DOMAIN_MARKER = "betonline.ag"
REFRESH_TOKEN_COOKIE = "krefresh"
PAGE_SIZE = 100
MAX_PAGES = 50
HTTP_TIMEOUT_SEC = 15
# Refresh the access token only this close to its expiry (issue #115 rule).
REFRESH_MARGIN_SEC = 60
# When the token response carries no expires_in, assume Keycloak's usual 5 min.
DEFAULT_ACCESS_TTL_SEC = 300

# The report's own wager id. The field name is unverified until the first live
# pull: a row carrying none of these fails the poll loudly (never a hash of
# mutable fields — the store upserts on the id).
ID_FIELD_CANDIDATES = ("Id", "TicketNumber", "WagerNumber", "TicketId", "WagerId", "BetId")

SOURCE = "betonline_api"
VENUE = "betonline"
DEVICE_PREFIX_RE = re.compile(r"^(?:Desktop|Mobile)\s*-\s*", re.IGNORECASE)
LEAGUE_TOKENS = {
    "NFL": "nfl", "NCAAF": "cfb", "CFB": "cfb", "COLLEGE FOOTBALL": "cfb",
    "NBA": "nba", "NCAAB": "cbb", "NCAAM": "cbb", "CBB": "cbb", "COLLEGE BASKETBALL": "cbb",
    "WNBA": "wnba", "MLB": "mlb", "NHL": "nhl", "HOCKEY": "nhl", "SOCCER": "soccer",
}
# A sport prefix names two leagues the scanner keeps apart; fail closed.
SPORT_ONLY_TOKENS = {"FOOTBALL", "BASKETBALL", "BASEBALL"}
PERIODS = {
    "1st half": "1H", "first half": "1H", "2nd half": "2H", "second half": "2H",
    "1st quarter": "1Q", "2nd quarter": "2Q", "3rd quarter": "3Q", "4th quarter": "4Q",
    "1st 5 innings": "F5", "first 5 innings": "F5", "1st five innings": "F5",
    "first five innings": "F5", "1st inning": "I1", "first inning": "I1",
}
STATUS_MAP = {
    "pending": "open", "won": "won", "lost": "lost", "push": "push",
    "cancelled": "void", "canceled": "void", "void": "void",
}
WAGER_TYPE_MAP = {"spread": "spread", "total": "total", "money line": "moneyline",
                  "moneyline": "moneyline"}
PARLAY_WAGER_TYPES = {"same game parlay", "parlay"}
PARLAY_LEG_SEPARATORS = (" | ", "; ", ", ")
FRACTIONS = {"½": ".5", "¼": ".25", "¾": ".75"}
APPROX_DATE_UNKNOWN = "game_date_unknown"
APPROX_SIDE_PARITY = "side_from_rotation_parity"

NUMBER = r"(?:\d+(?:\.\d+)?[½¼¾]?|[½¼¾])"
PRICE = r"[+-]\d{3,}"
# A team name never carries a signed number or a leg separator, so a moneyline
# regex cannot swallow "Bears +3.5 -110, <another leg>" as one team.
TEAM = r"(?P<team>(?:(?![+-]\d)[^,;|])+?)"
LEG_RE = re.compile(r"^(?P<rotation>\d+)\s+(?P<body>.+)$")
TOTAL_RE = re.compile(rf"^{TEAM}\s+(?P<direction>over|under)\s+(?P<points>{NUMBER})\s+(?P<price>{PRICE})$",
                      re.IGNORECASE)
SPREAD_RE = re.compile(rf"^{TEAM}\s+(?P<points>[+-]{NUMBER}|pk|pick)\s+(?P<price>{PRICE})$", re.IGNORECASE)
MONEYLINE_RE = re.compile(rf"^{TEAM}\s+(?P<price>{PRICE})$")
FOR_TAIL_RE = re.compile(r"\s+for\s+.*$", re.IGNORECASE)


# ---- pure parser --------------------------------------------------------------------

def parse_points(text: str) -> float:
    lowered = text.lower()
    if lowered in ("pk", "pick"):
        return 0.0
    for glyph, decimal in FRACTIONS.items():
        text = text.replace(glyph, decimal)
    if text.startswith("."):
        text = "0" + text
    if text.startswith(("+.", "-.")):
        text = text[0] + "0" + text[1:]
    return float(text)


def split_period(text: str) -> tuple[str, str | None, str | None]:
    """"<body> - 1st Half" -> (body, "1H", None); no tail -> FG; an unknown tail
    -> (body, None, reason)."""
    head, separator, tail = text.rpartition(" - ")
    if not separator:
        return text, "FG", None
    period = PERIODS.get(tail.strip().lower())
    if period is None:
        return head, None, f"unknown period ({tail.strip()})"
    return head, period, None


def parse_leg(text: str) -> dict | str:
    """"465 Chicago Bears +3.5 -110 - 1st Half" -> {rotation, team, betType, side,
    points, price, period}, or the reason it does not parse."""
    leg_match = LEG_RE.match(text.strip())
    if not leg_match:
        return f"no rotation number ({text.strip()[:60]})"
    rotation = int(leg_match.group("rotation"))
    body, period, period_error = split_period(FOR_TAIL_RE.sub("", leg_match.group("body")).strip())
    if period_error:
        return period_error
    parsed = {"rotation": rotation, "period": period}
    total = TOTAL_RE.match(body)
    if total:
        parsed.update(team=total.group("team"), betType="total", side=total.group("direction").lower(),
                      points=parse_points(total.group("points")), price=int(total.group("price")))
        return parsed
    spread = SPREAD_RE.match(body)
    if spread:
        parsed.update(team=spread.group("team"), betType="spread", side=side_from_rotation(rotation),
                      points=parse_points(spread.group("points")), price=int(spread.group("price")))
        return parsed
    moneyline = MONEYLINE_RE.match(body)
    if moneyline:
        parsed.update(team=moneyline.group("team"), betType="moneyline", side=side_from_rotation(rotation),
                      points=None, price=int(moneyline.group("price")))
        return parsed
    return f"unrecognised selection ({body[:60]})"


def side_from_rotation(rotation: int) -> str:
    """US rotation convention: the away (visiting) team carries the odd number."""
    return "away" if rotation % 2 == 1 else "home"


def split_league(description: str) -> tuple[str | None, str, str | None]:
    """"Desktop - NFL - 465 ..." -> ("nfl", "465 ...", None); a sport-only or
    unknown prefix -> (None, rest, reason)."""
    text = DEVICE_PREFIX_RE.sub("", description.strip())
    token, separator, rest = text.partition(" - ")
    if not separator:
        return None, text, "no league prefix"
    upper = token.strip().upper()
    league = LEAGUE_TOKENS.get(upper)
    if league:
        return league, rest.strip(), None
    if upper in SPORT_ONLY_TOKENS:
        return None, rest.strip(), f"league not determined from sport prefix ({token.strip()})"
    return None, rest.strip(), f"unknown league prefix ({token.strip()})"


def split_parlay_legs(rest: str) -> list[dict] | str:
    """Legs of a Same Game Parlay description (after the league prefix). A second
    repeated league token is dropped. Every piece must parse, else the reason."""
    token, separator, tail = rest.partition(" - ")
    if separator and token.strip().upper() in LEAGUE_TOKENS:
        rest = tail.strip()
    for candidate in PARLAY_LEG_SEPARATORS:
        pieces = [piece.strip() for piece in rest.split(candidate)]
        if len(pieces) < 2:
            continue
        legs = [parse_leg(piece) for piece in pieces]
        if all(isinstance(leg, dict) for leg in legs):
            return legs
    single = parse_leg(rest)
    if isinstance(single, dict):
        return [single]
    return f"parlay legs not parsed ({rest[:80]})"


def native_id_of(row: dict) -> str:
    for key in ID_FIELD_CANDIDATES:
        value = row.get(key)
        if value not in (None, ""):
            return str(value)
    raise RuntimeError(
        f"BetOnline report row carries none of the id fields {ID_FIELD_CANDIDATES}; "
        f"keys seen: {sorted(row.keys())}")


def parse_placed_at(value: object) -> str | None:
    if not isinstance(value, str) or not value:
        return None
    try:
        parsed = datetime.fromisoformat(value.replace("Z", "+00:00"))
    except ValueError:
        return None
    if parsed.tzinfo is None:
        parsed = parsed.replace(tzinfo=timezone.utc)
    return parsed.astimezone(timezone.utc).strftime("%Y-%m-%dT%H:%M:%SZ")


def status_of(wager_status: object) -> str:
    return STATUS_MAP.get(str(wager_status or "").strip().lower(), "unknown")


def american_from_payout(risk: float, to_win: float) -> int | None:
    if risk <= 0 or to_win <= 0:
        return None
    ratio = to_win / risk
    return round(ratio * 100) if ratio >= 1 else -round(100 / ratio)


def _money(value: object) -> float | None:
    if value is None or isinstance(value, bool):
        return None
    try:
        return round_cents(float(value))
    except (TypeError, ValueError):
        return None


def _base_record(row: dict, native_id: str, fetched_at: str | None) -> dict:
    status = status_of(row.get("WagerStatus"))
    placed_at = parse_placed_at(row.get("Date"))
    return {
        "id": f"{VENUE}:{native_id}",
        "source": SOURCE,
        "venue": VENUE,
        "league": None, "eventStart": None, "eventDate": None,
        "awayTeam": None, "homeTeam": None, "awayKey": None, "homeKey": None,
        "rotation": None, "betType": "other", "period": None, "side": None, "points": None,
        "price": None,
        "stake": _money(row.get("Risk")),
        "toWin": _money(row.get("ToWin")),
        "contracts": None,
        "placedAt": placed_at,
        "status": status,
        # The report carries no settle time: the placed time is the lower bound.
        "closedAt": None if status == "open" else placed_at,
        "isParlayLeg": False, "parlayId": None, "legIndex": None, "legCount": None,
        "approx": [],
        "unmatchable": None,
        "sourceFetchedAt": fetched_at,
        "raw": {
            "nativeId": native_id,
            "description": row.get("Description"),
            "wagerType": row.get("WagerType"),
            "wagerStatus": row.get("WagerStatus"),
            "risk": row.get("Risk"),
            "toWin": row.get("ToWin"),
            "date": row.get("Date"),
        },
    }


def _apply_leg(record: dict, league: str, leg: dict) -> dict:
    record.update({
        "league": league, "rotation": leg["rotation"], "betType": leg["betType"],
        "period": leg["period"], "side": leg["side"], "points": leg["points"], "price": leg["price"],
        "approx": [APPROX_DATE_UNKNOWN] + ([APPROX_SIDE_PARITY] if leg["betType"] != "total" else []),
    })
    if side_from_rotation(leg["rotation"]) == "away":
        record["awayTeam"] = leg["team"]
    else:
        record["homeTeam"] = leg["team"]
    record["raw"]["team"] = leg["team"]
    return record


def _unmatchable(record: dict, reason: str) -> dict:
    record["unmatchable"] = reason
    return record


def normalize_row(row: dict, fetched_at: str | None) -> list[dict]:
    """One report row -> one record, or one record per leg of a parlay. Pure."""
    native_id = native_id_of(row)
    base = _base_record(row, native_id, fetched_at)
    league, rest, league_error = split_league(str(row.get("Description") or ""))
    wager_type = str(row.get("WagerType") or "").strip().lower()
    if wager_type in PARLAY_WAGER_TYPES:
        parlay_price = american_from_payout(base["stake"] or 0, base["toWin"] or 0)
        base["raw"]["parlayPrice"] = parlay_price
        if league_error:
            return [json_clean(_unmatchable(base, league_error))]
        legs = split_parlay_legs(rest)
        if isinstance(legs, str):
            return [json_clean(_unmatchable(base, legs))]
        records = []
        for index, leg in enumerate(legs):
            record = _base_record(row, native_id, fetched_at)
            record["id"] = f"{VENUE}:{native_id}:leg{index}"
            record["raw"]["parlayPrice"] = parlay_price
            record.update({"isParlayLeg": True, "parlayId": base["id"], "legIndex": index,
                           "legCount": len(legs)})
            records.append(json_clean(_apply_leg(record, league, leg)))
        return records
    if league_error:
        return [json_clean(_unmatchable(base, league_error))]
    leg = parse_leg(rest)
    if isinstance(leg, str):
        return [json_clean(_unmatchable(base, leg))]
    expected = WAGER_TYPE_MAP.get(wager_type)
    if expected is None:
        return [json_clean(_unmatchable(base, f"unknown wager type ({row.get('WagerType')})"))]
    if expected != leg["betType"]:
        return [json_clean(_unmatchable(
            base, f"wager type {row.get('WagerType')} but the selection reads as {leg['betType']}"))]
    return [json_clean(_apply_leg(base, league, leg))]


def normalize_betonline(rows: list[dict], fetched_at: str | None) -> list[dict]:
    """Every report row -> records (pending bets KEPT; a parlay is one record per
    leg). Raises only when a row has no native id — the one thing the store
    cannot work without."""
    records: list[dict] = []
    for row in rows:
        records.extend(normalize_row(row, fetched_at))
    return records


# ---- network half -------------------------------------------------------------------

def read_cookies(path: Path) -> list[dict]:
    if not path.exists():
        raise RuntimeError(f"BetOnline cookie file missing: {path} — run bet_logger/recon_betonline.py")
    cookies = json.loads(path.read_text())
    if not isinstance(cookies, list):
        raise RuntimeError(f"BetOnline cookie file {path} is not a JSON list")
    return cookies


def write_cookies_atomic(path: Path, cookies: list[dict]) -> None:
    """Temp file in the same directory, then os.replace: a reader never sees a
    half-written token file."""
    temp_path = path.with_name(path.name + ".tmp")
    temp_path.write_text(json.dumps(cookies, indent=2))
    os.replace(temp_path, path)


def refresh_token_of(cookies: list[dict]) -> str:
    for cookie in cookies:
        if cookie.get("name") == REFRESH_TOKEN_COOKIE and cookie.get("value"):
            return cookie["value"]
    raise RuntimeError(f"no {REFRESH_TOKEN_COOKIE} cookie in the BetOnline cookie file — "
                       "run bet_logger/recon_betonline.py")


def session_with_cookies(cookies: list[dict]) -> requests.Session:
    session = requests.Session()
    for cookie in cookies:
        if COOKIE_DOMAIN_MARKER in cookie.get("domain", ""):
            session.cookies.set(cookie["name"], cookie["value"], domain=cookie["domain"],
                                path=cookie.get("path", "/"))
    return session


def api_headers(access_token: str, now: datetime) -> dict:
    ms = int(now.timestamp() * 1000)
    return {
        "Authorization": f"Bearer {access_token}",
        "Content-Type": "application/json",
        **BASE_HEADERS,
        "actual-time": str(ms),
        "iso-time": now.strftime("%Y-%m-%dT%H:%M:%S.") + f"{ms % 1000:03d}Z",
        "utc-time": now.strftime("%a, %d %b %Y %H:%M:%S GMT"),
    }


@contextmanager
def cookie_file_lock(path: Path) -> Iterator[None]:
    """Exclusive flock on `<cookie file>.lock` for the read -> refresh -> write
    sequence. bet_logger/scraper_betonline.py takes the same lock (same path
    rule), so the two never rotate the same refresh token."""
    lock_path = path.with_name(path.name + ".lock")
    with open(lock_path, "w") as lock_file:
        fcntl.flock(lock_file, fcntl.LOCK_EX)
        try:
            yield
        finally:
            fcntl.flock(lock_file, fcntl.LOCK_UN)


class BetOnlineSource:
    """Source protocol implementation for the BetOnline account (report POSTs only).

    `session_factory(cookies)` builds the HTTP session (requests by default; tests
    inject a fake with .post()). The access token lives in memory with its expiry;
    the rotated refresh token goes back to the cookie file.
    """

    name = VENUE

    def __init__(self, cookies_path: Path | None = None, poll_sec: float | None = None,
                 history_days: int | None = None,
                 session_factory: Callable[[list[dict]], object] = session_with_cookies,
                 clock: Callable[[], float] = time.time,
                 local_today: Callable[[], datetime] = datetime.now):
        self.poll_sec = poll_sec if poll_sec is not None else config.BETONLINE_POLL_SEC
        self._cookies_path = Path(cookies_path or config.BETONLINE_COOKIES_PATH)
        self._history_days = history_days if history_days is not None else config.BETONLINE_HISTORY_DAYS
        self._session_factory = session_factory
        self._clock = clock
        self._local_today = local_today
        self._access_token: str | None = None
        self._access_expires_at = 0.0
        self._refresh_lock = threading.Lock()

    # -- token ------------------------------------------------------------------------

    def _access_token_fresh(self) -> bool:
        return (self._access_token is not None
                and self._access_expires_at - self._clock() > REFRESH_MARGIN_SEC)

    def _refresh_access_token(self) -> tuple[str, object]:
        """Under the cookie-file lock: re-read the file, exchange the refresh
        token, persist the rotated one. Returns (access_token, session)."""
        with cookie_file_lock(self._cookies_path):
            return self._refresh_access_token_locked()

    def _refresh_access_token_locked(self) -> tuple[str, object]:
        cookies = read_cookies(self._cookies_path)
        refresh_token = refresh_token_of(cookies)
        session = self._session_factory(cookies)
        response = session.post(
            TOKEN_URL,
            data={"grant_type": "refresh_token", "refresh_token": refresh_token, "client_id": CLIENT_ID},
            headers={"Content-Type": "application/x-www-form-urlencoded",
                     **{key: BASE_HEADERS[key] for key in ("User-Agent", "Origin", "Referer")}},
            timeout=HTTP_TIMEOUT_SEC,
        )
        if response.status_code != 200:
            error_code = None
            try:
                error_code = response.json().get("error")
            except Exception:  # noqa: BLE001 — a non-JSON error body is fine to ignore
                pass
            raise RuntimeError(
                f"BetOnline token refresh failed (HTTP {response.status_code}, {error_code}); "
                "if the refresh token expired, run bet_logger/recon_betonline.py")
        token = response.json()
        access_token = token["access_token"]
        ttl = float(token.get("expires_in") or DEFAULT_ACCESS_TTL_SEC)
        rotated = token.get("refresh_token")
        if rotated and rotated != refresh_token:
            for cookie in cookies:
                if cookie.get("name") == REFRESH_TOKEN_COOKIE:
                    cookie["value"] = rotated
            write_cookies_atomic(self._cookies_path, cookies)
            log.info("betonline: refresh token rotated and saved")
        self._access_token = access_token
        self._access_expires_at = self._clock() + ttl
        return access_token, session

    def _authenticated_session(self) -> tuple[str, object]:
        with self._refresh_lock:
            if self._access_token_fresh():
                return self._access_token, self._session_factory(read_cookies(self._cookies_path))
            return self._refresh_access_token()

    # -- report -----------------------------------------------------------------------

    def _report_window(self) -> tuple[str, str]:
        # Local dates: BetOnline returns 0 rows when EndDate is "tomorrow" in its
        # own timezone (bet_logger/CLAUDE.md), so never compute this in UTC.
        today = self._local_today()
        start = today - timedelta(days=self._history_days)
        return start.strftime("%Y-%m-%d"), today.strftime("%Y-%m-%d")

    def _fetch_report(self, session: object, access_token: str) -> list[dict]:
        start_date, end_date = self._report_window()
        headers = api_headers(access_token, datetime.now(timezone.utc))
        rows: list[dict] = []
        for page in range(MAX_PAGES):
            response = session.post(
                BET_HISTORY_URL, headers=headers,
                json={
                    "Id": None, "StartDate": f"{start_date}T00:00:00.000Z",
                    "EndDate": f"{end_date}T00:00:00.000Z", "Status": None, "Product": None,
                    "WagerType": None, "FreePlayFlag": None, "StartPosition": page * PAGE_SIZE,
                    "TotalPerPage": PAGE_SIZE, "IsDailyFigureReport": False,
                },
                timeout=HTTP_TIMEOUT_SEC,
            )
            if response.status_code != 200:
                raise RuntimeError(f"BetOnline report page {page} failed: HTTP {response.status_code}")
            body = response.json()
            page_rows = body.get("Data") or []
            total_rows = int(body.get("TotalRows") or 0)
            rows.extend(page_rows)
            if not page_rows or len(rows) >= total_rows:
                return rows
        raise RuntimeError(f"BetOnline report did not end within {MAX_PAGES} pages")

    # -- Source protocol --------------------------------------------------------------

    def fetch(self) -> list[dict]:
        access_token, session = self._authenticated_session()
        rows = self._fetch_report(session, access_token)
        records = normalize_betonline(rows, utc_now_iso())
        n_open = sum(1 for record in records if record["status"] == "open")
        n_unmatchable = sum(1 for record in records if record["unmatchable"])
        log.info("betonline: %d report rows -> %d records (%d open, %d unmatchable)",
                 len(rows), len(records), n_open, n_unmatchable)
        return records


def source_if_configured() -> BetOnlineSource | None:
    """The source when the recon cookie file exists (the Novig pattern: a venue
    with no credentials reads "no source configured" in the panel rather than
    failing every poll); a missing file is logged once with the fix."""
    if not config.BETONLINE_COOKIES_PATH.exists():
        log.warning("betonline: no cookie file at %s — run bet_logger/recon_betonline.py "
                    "--interactive to enable the BetOnline source", config.BETONLINE_COOKIES_PATH)
        return None
    return BetOnlineSource()
