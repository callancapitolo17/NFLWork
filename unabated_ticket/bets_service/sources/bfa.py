"""BFA (Betfastaction) bet source: the account's open bets + GetPlayerHistory -> normalised bet records.

Two halves, like sources/betonline.py:
  normalize_open_bets(wagers, fetched_at)  pure parser of GetPlayerOpenBets rows — the OPEN bets,
                                           the ones the panel flags — pinned by the pytest on
                                           tests/fixtures/bets/bfa_history.json ("openBets").
  normalize_bfa(wagers, fetched_at)        pure parser of GetPlayerHistory wagers: what settles an
                                           open record once it leaves the open list (a store row
                                           stays open until a later poll says otherwise).
  BFASource                                the network half (Source protocol): a Keycloak PASSWORD
                                           login (PKCE) this process owns — access and refresh
                                           token in memory only, the access token refreshed within
                                           REFRESH_MARGIN_SEC of expiry and a failed refresh
                                           replaced by a new login — then the open bets, then the
                                           paged history over the retention window; an open-bets
                                           record replaces the history's copy of the same wager.

Inputs:  BFA_USERNAME / BFA_PASSWORD (config: the environment, bets_service/.env,
         kalshi_draft/.env, then bet_logger/.env in the main checkout — the sheet scraper's
         own credentials, no copying);
         GET api.bfagaming.com/history/api/GetPlayerOpenBets?playerId   (player id from the JWT)
         GET api.bfagaming.com/history/api/GetPlayerHistory
             ?playerId&startDate&endDate&page&recordsByPage
Outputs: list of records. Side effects: none on disk. bet_logger/recon_bfa_auth.json is
         neither read nor written: the weekly LaunchAgent (bet_logger/scraper_bfa.py) rotates
         that refresh token, and two processes rotating one token trip Keycloak's reuse
         detection (why BetOnline needs a file lock) — owning a session avoids the race.
         Raises on a failed login, a failed refresh that a new login cannot replace, or any
         non-200 history page, so the service records a failed run and keeps the previous
         records (never a partial list).

Why not import bet_logger/scraper_bfa.py: it imports Google Sheets at module level (not in
this venv), rotates the shared token file and skips pending bets. The endpoints, headers and
the description grammar are copied from it and bet_logger/recon_bfa.py and must stay in step.

Open bet shape (GetPlayerOpenBets, captured 2026-02-24 with two open bets; a list, one object
per wager):
  idWager        343243731 — the same id the history later carries (ticketNumber repeats it).
  headerDescription  STRAIGHT BET | PARLAY (2 TEAMS) | …;  riskAmount, winAmount  USD;  result 255.
  betDetails[]   one per leg: idSport (CBB / CFB / NFL / NBA / MLB / NHL …, the league the
                 history never names — so an open college bet IS placed), gameDateTime (the
                 event start), idGame, detailDescription
                 "CBB - Alternative Lines <br> [1674] TOTAL u68½+110 \r(NEW MEXICO 1H vrs NEVADA 1H) [Sport:Basketball, League:NCAA]"
                 — the market group, then the history's own leg grammar, then a sport suffix.
  placedDate, gameDateTime   the Pacific wall-clock PLUS 7 HOURS, whatever the season: that
                 wager, placed between two history rows stamped 15:32 and 15:33 PST, reads 22:33,
                 and its 8 PM PT tip reads 03:00. True UTC only under daylight time, so the 7 hours
                 are subtracted and the wall-clock localised (OPEN_BETS_CLOCK_OFFSET).
  GetPlayerOpenBetsWithOpenSpot (if-bets awaiting a leg) answered [] and is not read.

History wager grammar (live pull 2026-09-22, 16 wagers; the parser fails closed on anything
else, listing the record as unmatchable with the reason and the raw description):
  id            354669433 — the wager id, the stable native id.
  type          STRAIGHT BET | STRAIGHT BET (FP) (free play: risk 0) | PARLAY (2 TEAMS) |
                4 TEAM TEASERS. A parlay or teaser's `description` is its FIRST leg and
                `picks[]` the rest (each with a description and its own result); one record
                per leg like BetOnline, every leg on the ticket's status.
  description   "[1340] TOTAL u24EV \\r(ARIZONA 1H vrs BYU 1H)"   total: rotation, o/u + number,
                    price (EV = +100), "(AWAY vrs HOME)" away first, "1H" on the names = 1H;
                    a baseball total carries the pitchers in a second bracket, ignored
                "[308945] ILLINOIS ST -3-110"                    spread: its own team only
                "[477] GREEN BAY PACKERS +8-110 (B+6)"           teaser leg: "(B+6)" = bought
                    points, the number as shown is the side's own
                "[1306551] GRAMBLING 1H +168"                    moneyline (scraper_bfa.py;
                    not in the pull)
                "TOTAL u38½-115 (SAINT LOUIS 1H TEAM PTS vrs …)"  team total -> betType "other",
                    unmatchable
                ½ ¼ ¾ fractions; "\\r" before the bracket.
  result        "" / PENDING / OPEN = open; WIN / LOSE / PUSH / CANCELLED / VOID.
  risk, win     USD (`amount` is the net P&L and is not read).
  placedDate, lastModification, settledDate
                naive wall-clock on the account's own clock — Pacific: every straight bet in
                the pull carries its scheduled kickoff in settledDate (09:00 on the Kent State
                @ Ohio State and Buffalo @ Penn State 1H totals, both noon-ET games) and its
                grading time in lastModification (the Rice @ Notre Dame 1H total graded at
                halftime, 12:30 -> 14:03). So a straight bet's settledDate is its eventStart
                (flagged approx event_start_from_settled_date; refused outside placedAt
                −1 day .. +60 days, where a .NET default date would otherwise parse) and
                lastModification its closedAt. A parlay or teaser carries ONE settledDate for
                the whole ticket, so its legs get no start (game_date_unknown: the matcher
                windows placedAt and keys on the rotation, as for BetOnline).
A history description names no sport: the league comes from bet_logger/utils.py parse_sport
(its nickname scan resolves the pro leagues only), and a settled college game — most of this
account — is "league unknown" and unmatchable, never guessed as CFB or CBB (hand-off rule,
2026-09-23). Only settled bets are affected: an open bet's league comes from its idSport.
A spread or moneyline names one team, placed by rotation parity (odd = away, approx
side_from_rotation_parity) — the convention the pull's totals follow too: every over carries
an odd rotation, every under an even one. BFA's rotations are Unabated's own numbers (the
2026-09-22 pull's 308945 / 308959 and 461-481 were all on that week's board), except that a
first-half leg's rotation is the game's with a "1" PREPENDED — 1340 for game 340, 1306551
for 306551, 1670 for 670 — so a 1H leg's rotation is served with that digit stripped (the
number as written stays in raw.rotationAsWritten) and a one-team 1H bet can match by rotation.
"""
import base64
import hashlib
import json
import logging
import re
import secrets
import threading
import time
from collections.abc import Callable
from datetime import datetime, timedelta, timezone
from urllib.parse import parse_qs, urlparse
from zoneinfo import ZoneInfo

import requests

from unabated_ticket.bets_service import config
from unabated_ticket.bets_service.normalize import EASTERN, json_clean, round_cents, utc_now_iso
from unabated_ticket.bets_service.sources.betonline import (
    APPROX_DATE_UNKNOWN, APPROX_SIDE_PARITY, SHEET_LABEL_TO_LEAGUE, american_from_payout, parse_points,
    parse_sport, side_from_rotation)

log = logging.getLogger(__name__)

# Endpoints and headers: copied from bet_logger/recon_bfa.py + scraper_bfa.py (keep in step).
KEYCLOAK_BASE = "https://auth.bfagaming.com/realms/players_realm/protocol/openid-connect"
AUTH_URL = f"{KEYCLOAK_BASE}/auth"
TOKEN_URL = f"{KEYCLOAK_BASE}/token"
CLIENT_ID = "bfagaming"
REDIRECT_URI = "https://bfagaming.com/"
HISTORY_URL = "https://api.bfagaming.com/history/api/GetPlayerHistory"
OPEN_BETS_URL = "https://api.bfagaming.com/history/api/GetPlayerOpenBets"
USER_AGENT = ("Mozilla/5.0 (Macintosh; Intel Mac OS X 10_15_7) "
              "AppleWebKit/537.36 (KHTML, like Gecko) Chrome/145.0.0.0 Safari/537.36")
BASE_HEADERS = {"Accept": "application/json", "Origin": "https://bfagaming.com",
                "Referer": "https://bfagaming.com/", "User-Agent": USER_AGENT}
RECORDS_PER_PAGE = 100
MAX_PAGES = 50
HTTP_TIMEOUT_SEC = 20
# Refresh the access token only this close to its expiry (the BetOnline rule).
REFRESH_MARGIN_SEC = 60
# When the token response carries no expires_in, assume Keycloak's usual 5 min.
DEFAULT_ACCESS_TTL_SEC = 300
LOGIN_FORM_ACTION_RE = re.compile(r'action="([^"]+)"')
INVALID_CREDENTIALS_MARKER = "Invalid username or password"

SOURCE = "bfa_api"
VENUE = "bfa"
# The account renders every history timestamp on its own clock (module docstring).
BFA_TZ = ZoneInfo("America/Los_Angeles")
# GetPlayerOpenBets stamps the Pacific wall-clock plus 7 hours (module docstring).
OPEN_BETS_CLOCK_OFFSET = timedelta(hours=7)
# An open leg's idSport (the codes the account's line profile lists); the rest are props.
OPEN_BET_LEAGUES = {"CBB": "cbb", "CFB": "cfb", "NFL": "nfl", "NBA": "nba", "WNBA": "wnba", "MLB": "mlb",
                    "NHL": "nhl", "SOC": "soccer"}
OPEN_BET_PROP_CODES = {"PROP", "TNT", "MU", "ESOC"}
REASON_NO_OPEN_LEGS = "open bet carries no legs"
# "… [Sport:Basketball, League:NCAA]" closes an open bet's description.
SPORT_SUFFIX_RE = re.compile(r"\s*\[Sport:[^\]]*\]\s*$", re.IGNORECASE)
HTML_TAG_RE = re.compile(r"<[^>]+>")
APPROX_START_FROM_SETTLED = "event_start_from_settled_date"
EVENT_START_DAYS_BEFORE_PLACED = 1
EVENT_START_DAYS_AFTER_PLACED = 60
REASON_LEAGUE_UNKNOWN = "league unknown (BFA names no sport; a college game cannot be placed)"
REASON_TEAM_TOTAL = "team total"
FIRST_HALF_ROTATION_PREFIX = "1"
# A game rotation has at least three digits, so a prefixed one has at least four.
FIRST_HALF_ROTATION_MIN_DIGITS = 4
TEAM_TOTAL_MARKER = "TEAM PTS"

STATUS_MAP = {"": "open", "pending": "open", "open": "open", "win": "won", "lose": "lost", "push": "push",
              "cancelled": "void", "canceled": "void", "void": "void"}
# "PARLAY (2 TEAMS)" / "4 TEAM TEASERS": the leg count the ticket declares.
MULTI_LEG_TYPE_RE = re.compile(r"PARLAY\s*\((?P<parlay>\d+)\s*TEAMS?\)|(?P<teaser>\d+)\s*TEAM\s*TEASERS?", re.IGNORECASE)
PERIODS = {"1H": "1H", "2H": "2H", "1Q": "1Q", "2Q": "2Q", "3Q": "3Q", "4Q": "4Q"}

ROTATION_RE = re.compile(r"^\[(?P<rotation>\d+)\]\s*(?P<body>.*)$")
NUMBER = r"(?:\d+(?:\.\d+)?[½¼¾]?|[½¼¾])"
PRICE = r"(?:[+-]\d{3,}|EV)"
# "(AWAY vrs HOME)" after the price; anything after that bracket (a baseball
# total's pitchers bracket) is ignored. Without a bracket nothing may follow.
BRACKET = r"(?:\s*\((?P<context>[^()]*)\)(?P<extra>.*))?"
TOTAL_RE = re.compile(rf"^TOTAL\s+(?P<direction>[ou])(?P<points>{NUMBER})(?P<price>{PRICE}){BRACKET}$", re.IGNORECASE)
# A team name never carries a bracket: a prop like "GIANTS SUPERFECTA (SCR 1ST, 1Q, 1H
# & GM +7½) +1935" must not read as a Giants moneyline.
TEAM = r"(?P<team>[^()]+?)"
SPREAD_RE = re.compile(rf"^{TEAM}\s+(?P<points>[+-]{NUMBER}|PK)(?P<price>{PRICE}){BRACKET}$", re.IGNORECASE)
MONEYLINE_RE = re.compile(rf"^{TEAM}\s+(?P<price>{PRICE}){BRACKET}$", re.IGNORECASE)
# A teaser leg's bought points, "(B+6)": the number shown already includes them.
BOUGHT_POINTS_RE = re.compile(r"\s*\(B[+-]\d+(?:\.\d+)?[½¼¾]?\)")
PERIOD_TOKEN_RE = re.compile(r"^(?P<name>.+?)\s+(?P<token>\d[HQ])$", re.IGNORECASE)
VERSUS_RE = re.compile(r"\s+vrs\s+", re.IGNORECASE)


# ---- pure parser --------------------------------------------------------------------

def normalize_description(raw: str) -> str:
    """One line, single spaces, the teaser's bought-points bracket dropped."""
    text = raw.replace("\r", " ").replace("\n", " ")
    text = BOUGHT_POINTS_RE.sub("", text)
    return re.sub(r"\s+", " ", text).strip()


def parse_price(text: str) -> int:
    return 100 if text.upper() == "EV" else int(text)


def leg_segment_of(raw: str) -> str:
    """The "[rotation] …" segment of a description that wraps it in HTML — an open
    bet's "CBB - Alternative Lines <br> [1674] TOTAL … [Sport:…, League:…]" — else
    the whole text. A history description has no wrapping and passes through."""
    text = SPORT_SUFFIX_RE.sub("", raw)
    for segment in HTML_TAG_RE.split(text):
        if ROTATION_RE.match(normalize_description(segment)):
            return segment
    return text


def period_for(league: str, period: str | None) -> str:
    """The leg's period; an MLB "1H" is the first five innings (Novig's rule)."""
    if period is None:
        return "FG"
    if period == "1H" and league == "mlb":
        return "F5"
    return period


def split_period(name: str) -> tuple[str, str | None] | str:
    """"ARIZONA 1H" -> ("ARIZONA", "1H"); no token -> (name, None); an unknown
    token -> the reason."""
    match = PERIOD_TOKEN_RE.match(name.strip())
    if not match:
        return name.strip(), None
    period = PERIODS.get(match.group("token").upper())
    if period is None:
        return f"unknown period ({match.group('token')})"
    return match.group("name").strip(), period


def parse_teams_bracket(context: str) -> dict | str:
    """"ARIZONA 1H vrs BYU 1H" -> {awayTeam, homeTeam, period}, or the reason."""
    if TEAM_TOTAL_MARKER in context.upper():
        return REASON_TEAM_TOTAL
    parts = VERSUS_RE.split(context.strip())
    if len(parts) != 2:
        return f"no 'AWAY vrs HOME' bracket ({context.strip()[:60]})"
    away, home = split_period(parts[0]), split_period(parts[1])
    if isinstance(away, str):
        return away
    if isinstance(home, str):
        return home
    periods = {period for _, period in (away, home) if period}
    if len(periods) > 1:
        return f"periods differ ({away[1]} vs {home[1]})"
    return {"awayTeam": away[0], "homeTeam": home[0], "period": periods.pop() if periods else "FG"}


def _one_team_leg(rotation: int, bet_type: str, team_text: str, points: float | None, price: int,
                  context: str | None) -> dict | str:
    """A spread or moneyline: its own team, placed by rotation parity — or by the
    "(AWAY vrs HOME)" bracket when the description carries one (not seen live)."""
    split = split_period(team_text)
    if isinstance(split, str):
        return split
    team, period = split
    leg = {"rotation": rotation, "betType": bet_type, "points": points, "price": price}
    if context is None:
        side = side_from_rotation(rotation)
        away, home = (team, None) if side == "away" else (None, team)
        leg.update(side=side, period=period or "FG", awayTeam=away, homeTeam=home, sideFromParity=True)
        return leg
    bracket = parse_teams_bracket(context)
    if isinstance(bracket, str):
        return bracket
    if team.upper() == bracket["awayTeam"].upper():
        side = "away"
    elif team.upper() == bracket["homeTeam"].upper():
        side = "home"
    else:
        return f"team {team} is not in its bracket ({context.strip()[:60]})"
    leg.update(side=side, period=period or bracket["period"], awayTeam=bracket["awayTeam"],
               homeTeam=bracket["homeTeam"], sideFromParity=False)
    return leg


def game_rotation_of(rotation: int, period: str | None) -> int:
    """A first-half leg's rotation as written -> the game's (module docstring)."""
    written = str(rotation)
    if period == "1H" and len(written) >= FIRST_HALF_ROTATION_MIN_DIGITS and written.startswith(FIRST_HALF_ROTATION_PREFIX):
        return int(written[len(FIRST_HALF_ROTATION_PREFIX):])
    return rotation


def _with_game_rotation(leg: dict) -> dict:
    leg["rotationAsWritten"] = leg["rotation"]
    leg["rotation"] = game_rotation_of(leg["rotation"], leg["period"])
    return leg


def parse_leg(raw_description: str) -> dict | str:
    """"[1340] TOTAL u24EV \\r(ARIZONA 1H vrs BYU 1H)" -> {rotation (the game's), rotationAsWritten,
    betType, side, points, price, period, awayTeam, homeTeam, sideFromParity}, or the reason it
    does not parse."""
    leg = _parse_leg_as_written(raw_description)
    return _with_game_rotation(leg) if isinstance(leg, dict) else leg


def _parse_leg_as_written(raw_description: str) -> dict | str:
    text = normalize_description(leg_segment_of(raw_description))
    rotation_match = ROTATION_RE.match(text)
    if not rotation_match:
        return f"no '[rotation]' prefix ({text[:60]})"
    rotation = int(rotation_match.group("rotation"))
    body = rotation_match.group("body").strip()
    total = TOTAL_RE.match(body)
    if total:
        if total.group("context") is None:
            return f"total names no teams ({body[:60]})"
        bracket = parse_teams_bracket(total.group("context"))
        if isinstance(bracket, str):
            return bracket
        return {"rotation": rotation, "betType": "total",
                "side": "over" if total.group("direction").lower() == "o" else "under",
                "points": parse_points(total.group("points")), "price": parse_price(total.group("price")),
                "period": bracket["period"], "awayTeam": bracket["awayTeam"], "homeTeam": bracket["homeTeam"],
                "sideFromParity": False}
    spread = SPREAD_RE.match(body)
    if spread:
        return _one_team_leg(rotation, "spread", spread.group("team"), parse_points(spread.group("points")),
                             parse_price(spread.group("price")), spread.group("context"))
    moneyline = MONEYLINE_RE.match(body)
    if moneyline:
        return _one_team_leg(rotation, "moneyline", moneyline.group("team"), None,
                             parse_price(moneyline.group("price")), moneyline.group("context"))
    return f"unrecognised selection ({body[:60]})"


def league_of(description: str) -> str | None:
    """feed.LEAGUES path for the description, or None: the sheet resolver knows the
    pro leagues by nickname (and baseball by its pitchers bracket), nothing else."""
    return SHEET_LABEL_TO_LEAGUE.get(parse_sport(normalize_description(description)))


def parse_account_time(value: object) -> datetime | None:
    """A naive BFA timestamp on the account's clock -> aware UTC; unparseable -> None."""
    if not isinstance(value, str) or not value:
        return None
    try:
        parsed = datetime.fromisoformat(value)
    except ValueError:
        return None
    if parsed.tzinfo is None:
        parsed = parsed.replace(tzinfo=BFA_TZ)
    return parsed.astimezone(timezone.utc)


def parse_open_bets_time(value: object) -> datetime | None:
    """A GetPlayerOpenBets timestamp (Pacific wall-clock + 7 h) -> aware UTC."""
    if not isinstance(value, str) or not value:
        return None
    try:
        stamped = datetime.fromisoformat(value)
    except ValueError:
        return None
    if stamped.tzinfo is not None:
        return stamped.astimezone(timezone.utc)
    pacific_wall_clock = stamped - OPEN_BETS_CLOCK_OFFSET
    return pacific_wall_clock.replace(tzinfo=BFA_TZ).astimezone(timezone.utc)


def iso_utc(moment: datetime | None) -> str | None:
    return None if moment is None else moment.strftime("%Y-%m-%dT%H:%M:%SZ")


def eastern_date_of(moment: datetime) -> str:
    return moment.astimezone(EASTERN).strftime("%Y-%m-%d")


def event_start_of(settled_date: object, placed_at: datetime | None) -> datetime | None:
    """A straight bet's settledDate as its kickoff (module docstring), only when it
    sits within a day before and 60 days after the bet was placed."""
    start = parse_account_time(settled_date)
    if start is None or placed_at is None:
        return None
    if start < placed_at - timedelta(days=EVENT_START_DAYS_BEFORE_PLACED):
        return None
    if start > placed_at + timedelta(days=EVENT_START_DAYS_AFTER_PLACED):
        return None
    return start


def status_of(result: object) -> str:
    return STATUS_MAP.get(str(result or "").strip().lower(), "unknown")


def native_id_of(wager: dict) -> str:
    value = wager.get("id")
    if value in (None, ""):
        raise RuntimeError(f"BFA wager carries no id; keys seen: {sorted(wager.keys())}")
    return str(value)


def _money(value: object) -> float | None:
    if value is None or isinstance(value, bool):
        return None
    try:
        return round_cents(float(value))
    except (TypeError, ValueError):
        return None


def _base_record(wager: dict, native_id: str, fetched_at: str | None) -> dict:
    status = status_of(wager.get("result"))
    placed_at = parse_account_time(wager.get("placedDate"))
    graded_at = parse_account_time(wager.get("lastModification"))
    closed_at = None if status == "open" else (graded_at or placed_at)
    return {
        "id": f"{VENUE}:{native_id}",
        "source": SOURCE,
        "venue": VENUE,
        "league": None, "eventStart": None, "eventDate": None,
        "awayTeam": None, "homeTeam": None, "awayKey": None, "homeKey": None,
        "rotation": None, "betType": "other", "period": None, "side": None, "points": None,
        "price": None,
        "stake": _money(wager.get("risk")),
        "toWin": _money(wager.get("win")),
        "contracts": None,
        "placedAt": iso_utc(placed_at),
        "status": status,
        "closedAt": iso_utc(closed_at),
        "isParlayLeg": False, "parlayId": None, "legIndex": None, "legCount": None,
        "approx": [],
        "unmatchable": None,
        "sourceFetchedAt": fetched_at,
        "raw": {
            "nativeId": native_id,
            "description": wager.get("description"),
            "type": wager.get("type"),
            "result": wager.get("result"),
            "risk": wager.get("risk"),
            "win": wager.get("win"),
            "amount": wager.get("amount"),
            "placedDate": wager.get("placedDate"),
            "lastModification": wager.get("lastModification"),
            "settledDate": wager.get("settledDate"),
        },
    }


def _apply_leg(record: dict, league: str, leg: dict, event_start: datetime | None,
               start_approx: str | None = None) -> dict:
    """`start_approx` names why the start is weaker than it looks (the history's
    settledDate); an open bet's gameDateTime needs no flag."""
    approx = []
    if event_start is None:
        approx.append(APPROX_DATE_UNKNOWN)
    elif start_approx:
        approx.append(start_approx)
    if leg["sideFromParity"]:
        approx.append(APPROX_SIDE_PARITY)
    record.update({
        "league": league,
        "eventStart": iso_utc(event_start),
        "eventDate": None if event_start is None else eastern_date_of(event_start),
        "rotation": leg["rotation"], "betType": leg["betType"], "period": period_for(league, leg["period"]),
        "side": leg["side"], "points": leg["points"], "price": leg["price"],
        "awayTeam": leg["awayTeam"], "homeTeam": leg["homeTeam"],
        "approx": approx,
    })
    record["raw"]["rotationAsWritten"] = leg["rotationAsWritten"]
    return record


def _unmatchable(record: dict, reason: str) -> dict:
    record["unmatchable"] = reason
    return record


def _leg_texts(wager: dict) -> list[str]:
    """A parlay or teaser's legs: the ticket's own description first, then each pick's."""
    picks = wager.get("picks") or []
    return [str(wager.get("description") or "")] + [str(pick.get("description") or "") for pick in picks]


def _normalize_multi_leg(wager: dict, native_id: str, fetched_at: str | None, declared_legs: int) -> list[dict]:
    base = _base_record(wager, native_id, fetched_at)
    parlay_price = american_from_payout(base["stake"] or 0, base["toWin"] or 0)
    base["raw"]["parlayPrice"] = parlay_price
    leg_texts = _leg_texts(wager)
    if len(leg_texts) != declared_legs:
        return [json_clean(_unmatchable(
            base, f"{wager.get('type')} names {declared_legs} legs but carries {len(leg_texts)}"))]
    legs = [parse_leg(text) for text in leg_texts]
    failed = next((leg for leg in legs if isinstance(leg, str)), None)
    if failed is not None:
        return [json_clean(_unmatchable(base, f"leg not parsed: {failed}"))]
    picks = wager.get("picks") or []
    records = []
    for index, (text, leg) in enumerate(zip(leg_texts, legs)):
        record = _base_record(wager, native_id, fetched_at)
        record["id"] = f"{base['id']}:leg{index}"
        record["raw"]["parlayPrice"] = parlay_price
        record["raw"]["legDescription"] = text
        # The ticket's own leg has no result of its own; each pick carries one.
        record["raw"]["legResult"] = None if index == 0 else picks[index - 1].get("result")
        record.update({"isParlayLeg": True, "parlayId": base["id"], "legIndex": index, "legCount": len(legs)})
        league = league_of(text)
        if league is None:
            records.append(json_clean(_unmatchable(record, REASON_LEAGUE_UNKNOWN)))
            continue
        records.append(json_clean(_apply_leg(record, league, leg, None)))
    return records


def normalize_wager(wager: dict, fetched_at: str | None) -> list[dict]:
    """One history wager -> one record, or one record per leg of a parlay/teaser. Pure."""
    native_id = native_id_of(wager)
    multi_leg = MULTI_LEG_TYPE_RE.search(str(wager.get("type") or ""))
    if multi_leg:
        return _normalize_multi_leg(wager, native_id, fetched_at,
                                    int(multi_leg.group("parlay") or multi_leg.group("teaser")))
    base = _base_record(wager, native_id, fetched_at)
    description = str(wager.get("description") or "")
    leg = parse_leg(description)
    if isinstance(leg, str):
        return [json_clean(_unmatchable(base, leg))]
    league = league_of(description)
    if league is None:
        return [json_clean(_unmatchable(base, REASON_LEAGUE_UNKNOWN))]
    event_start = event_start_of(wager.get("settledDate"), parse_account_time(wager.get("placedDate")))
    return [json_clean(_apply_leg(base, league, leg, event_start, APPROX_START_FROM_SETTLED))]


def normalize_bfa(wagers: list[dict], fetched_at: str | None) -> list[dict]:
    """Every history wager -> records (pending bets KEPT; a parlay or teaser is one
    record per leg). Raises only when a wager has no id — the one thing the store
    cannot work without."""
    records: list[dict] = []
    for wager in wagers:
        records.extend(normalize_wager(wager, fetched_at))
    return records


# ---- open bets ----------------------------------------------------------------------

def open_native_id_of(wager: dict) -> str:
    value = wager.get("idWager")
    if value in (None, "", 0, "0"):
        raise RuntimeError(f"BFA open bet carries no idWager; keys seen: {sorted(wager.keys())}")
    return str(value)


def open_bet_league_of(sport_code: object, description: str) -> tuple[str | None, str | None]:
    """(league, reason): the leg's idSport first, the nickname scan for a code the
    table does not know; a prop code or no league at all is the reason."""
    code = str(sport_code or "").strip().upper()
    if code in OPEN_BET_PROP_CODES:
        return None, f"not a game market (idSport {code})"
    league = OPEN_BET_LEAGUES.get(code) or league_of(description)
    if league is None:
        return None, f"league not supported (idSport {code or 'blank'})"
    return league, None


def _open_base_record(wager: dict, native_id: str, fetched_at: str | None) -> dict:
    placed_at = parse_open_bets_time(wager.get("placedDate"))
    return {
        "id": f"{VENUE}:{native_id}",
        "source": SOURCE,
        "venue": VENUE,
        "league": None, "eventStart": None, "eventDate": None,
        "awayTeam": None, "homeTeam": None, "awayKey": None, "homeKey": None,
        "rotation": None, "betType": "other", "period": None, "side": None, "points": None,
        "price": None,
        "stake": _money(wager.get("riskAmount")),
        "toWin": _money(wager.get("winAmount")),
        "contracts": None,
        "placedAt": iso_utc(placed_at),
        "status": "open",
        "closedAt": None,
        "isParlayLeg": False, "parlayId": None, "legIndex": None, "legCount": None,
        "approx": [],
        "unmatchable": None,
        "sourceFetchedAt": fetched_at,
        "raw": {
            "nativeId": native_id,
            "openBet": True,
            "headerDescription": wager.get("headerDescription"),
            "riskAmount": wager.get("riskAmount"),
            "winAmount": wager.get("winAmount"),
            "placedDate": wager.get("placedDate"),
            "wagerType": wager.get("wagerType"),
            "ifBetWagerType": wager.get("ifBetWagerType"),
        },
    }


def _open_leg_record(record: dict, leg_row: dict) -> dict:
    description = str(leg_row.get("detailDescription") or "")
    record["raw"].update({
        "description": description,
        "idSport": leg_row.get("idSport"),
        "gameDateTime": leg_row.get("gameDateTime"),
        "idGame": leg_row.get("idGame"),
    })
    league, reason = open_bet_league_of(leg_row.get("idSport"), description)
    if reason:
        return _unmatchable(record, reason)
    leg = parse_leg(description)
    if isinstance(leg, str):
        return _unmatchable(record, leg)
    return _apply_leg(record, league, leg, parse_open_bets_time(leg_row.get("gameDateTime")))


def normalize_open_wager(wager: dict, fetched_at: str | None) -> list[dict]:
    """One GetPlayerOpenBets wager -> one open record, or one per leg. Pure."""
    native_id = open_native_id_of(wager)
    legs = wager.get("betDetails") or []
    base = _open_base_record(wager, native_id, fetched_at)
    if not legs:
        return [json_clean(_unmatchable(base, REASON_NO_OPEN_LEGS))]
    if len(legs) == 1:
        return [json_clean(_open_leg_record(base, legs[0]))]
    parlay_price = american_from_payout(base["stake"] or 0, base["toWin"] or 0)
    records = []
    for index, leg_row in enumerate(legs):
        record = _open_base_record(wager, native_id, fetched_at)
        record["id"] = f"{base['id']}:leg{index}"
        record["raw"]["parlayPrice"] = parlay_price
        record.update({"isParlayLeg": True, "parlayId": base["id"], "legIndex": index, "legCount": len(legs)})
        records.append(json_clean(_open_leg_record(record, leg_row)))
    return records


def normalize_open_bets(wagers: list[dict], fetched_at: str | None) -> list[dict]:
    """Every open wager -> records. Raises only when one has no idWager."""
    records: list[dict] = []
    for wager in wagers:
        records.extend(normalize_open_wager(wager, fetched_at))
    return records


def merge_open_over_history(history_records: list[dict], open_records: list[dict]) -> list[dict]:
    """One list, an open-bets record replacing the history's copy of the same id:
    the history lists a pending wager too, with less (no league code, no start)."""
    by_id = {record["id"]: record for record in history_records}
    for record in open_records:
        by_id[record["id"]] = record
    return list(by_id.values())


# ---- network half -------------------------------------------------------------------

def pkce_pair() -> tuple[str, str]:
    verifier = secrets.token_urlsafe(32)
    challenge = base64.urlsafe_b64encode(hashlib.sha256(verifier.encode()).digest()).rstrip(b"=").decode()
    return verifier, challenge


def jwt_payload(token: str) -> dict:
    parts = token.split(".")
    if len(parts) != 3:
        return {}
    padded = parts[1] + "=" * (-len(parts[1]) % 4)
    try:
        return json.loads(base64.urlsafe_b64decode(padded))
    except ValueError:
        return {}


def new_session() -> requests.Session:
    session = requests.Session()
    session.headers["User-Agent"] = USER_AGENT
    return session


class BFASource:
    """Source protocol implementation for the BFA account (history GETs only).

    `session_factory()` builds an HTTP session with .get()/.post() (requests by
    default; tests inject a fake). The Keycloak session — access token with its
    expiry, refresh token, player id — lives in memory and nowhere else.
    """

    name = VENUE

    def __init__(self, username: str | None = None, password: str | None = None,
                 history_days: int | None = None, poll_sec: float | None = None,
                 session_factory: Callable[[], object] = new_session,
                 clock: Callable[[], float] = time.time):
        self.poll_sec = poll_sec if poll_sec is not None else config.BFA_POLL_SEC
        self._username = username if username is not None else config.BFA_USERNAME
        self._password = password if password is not None else config.BFA_PASSWORD
        self._history_days = history_days if history_days is not None else config.BFA_HISTORY_DAYS
        self._session_factory = session_factory
        self._clock = clock
        self._access_token: str | None = None
        self._access_expires_at = 0.0
        self._refresh_token: str | None = None
        self._player_id: str | None = None
        self._token_lock = threading.Lock()

    # -- Keycloak session -------------------------------------------------------------

    def _take_tokens(self, token: dict, how: str) -> None:
        access_token = token.get("access_token")
        if not access_token:
            raise RuntimeError(f"BFA {how} returned no access token")
        player_id = jwt_payload(access_token).get("player_id")
        if not player_id:
            raise RuntimeError(f"BFA {how} access token carries no player_id")
        self._access_token = access_token
        self._player_id = str(player_id)
        self._access_expires_at = self._clock() + float(token.get("expires_in") or DEFAULT_ACCESS_TTL_SEC)
        self._refresh_token = token.get("refresh_token") or self._refresh_token

    def _login(self) -> None:
        """Keycloak password login with PKCE: auth page -> form POST -> 302 with the
        code -> token exchange. Nothing about it is persisted."""
        if not self._username or not self._password:
            raise RuntimeError("BFA credentials missing: set BFA_USERNAME and BFA_PASSWORD "
                               "(bet_logger/.env in the main checkout)")
        session = self._session_factory()
        verifier, challenge = pkce_pair()
        login_page = session.get(
            AUTH_URL,
            params={"client_id": CLIENT_ID, "redirect_uri": REDIRECT_URI, "response_type": "code",
                    "scope": "openid", "code_challenge": challenge, "code_challenge_method": "S256"},
            timeout=HTTP_TIMEOUT_SEC)
        if login_page.status_code != 200:
            raise RuntimeError(f"BFA login page: HTTP {login_page.status_code}")
        action = LOGIN_FORM_ACTION_RE.search(login_page.text)
        if not action:
            raise RuntimeError("BFA login page carries no form action")
        submitted = session.post(action.group(1).replace("&amp;", "&"),
                                 data={"username": self._username, "password": self._password},
                                 allow_redirects=False, timeout=HTTP_TIMEOUT_SEC)
        if submitted.status_code not in (302, 303):
            detail = "invalid credentials" if INVALID_CREDENTIALS_MARKER in submitted.text else "no redirect"
            raise RuntimeError(f"BFA login failed: HTTP {submitted.status_code} ({detail})")
        code = parse_qs(urlparse(submitted.headers.get("Location", "")).query).get("code", [None])[0]
        if not code:
            raise RuntimeError("BFA login redirect carries no auth code")
        exchanged = session.post(
            TOKEN_URL,
            data={"grant_type": "authorization_code", "client_id": CLIENT_ID, "code": code,
                  "redirect_uri": REDIRECT_URI, "code_verifier": verifier},
            timeout=HTTP_TIMEOUT_SEC)
        if exchanged.status_code != 200:
            raise RuntimeError(f"BFA token exchange failed: HTTP {exchanged.status_code}")
        self._take_tokens(exchanged.json(), "login")
        log.info("bfa: logged in")

    def _refresh(self) -> bool:
        """Exchange the in-memory refresh token; False when Keycloak refuses it (the
        caller logs in again — nothing on disk to invalidate)."""
        session = self._session_factory()
        response = session.post(
            TOKEN_URL,
            data={"grant_type": "refresh_token", "refresh_token": self._refresh_token, "client_id": CLIENT_ID},
            timeout=HTTP_TIMEOUT_SEC)
        if response.status_code != 200:
            log.info("bfa: token refresh refused (HTTP %d); logging in again", response.status_code)
            return False
        self._take_tokens(response.json(), "refresh")
        return True

    def _ensure_access_token(self) -> None:
        with self._token_lock:
            if self._access_token is not None and self._access_expires_at - self._clock() > REFRESH_MARGIN_SEC:
                return
            if self._refresh_token is not None and self._refresh():
                return
            self._login()

    # -- history ----------------------------------------------------------------------

    def _history_window(self) -> tuple[str, str]:
        """(startDate, endDate) as YYYY-MM-DD on the account's clock: the retention
        window back, through tomorrow so a bet placed tonight is inside it."""
        today = datetime.fromtimestamp(self._clock(), BFA_TZ).date()
        return (today - timedelta(days=self._history_days)).isoformat(), (today + timedelta(days=1)).isoformat()

    def _fetch_open_bets(self) -> list[dict]:
        session = self._session_factory()
        response = session.get(OPEN_BETS_URL, headers={"Authorization": f"Bearer {self._access_token}", **BASE_HEADERS},
                               params={"playerId": self._player_id}, timeout=HTTP_TIMEOUT_SEC)
        if response.status_code != 200:
            raise RuntimeError(f"BFA open bets failed: HTTP {response.status_code}")
        body = response.json()
        if not isinstance(body, list):
            raise RuntimeError(f"BFA open bets: expected a list, got {type(body).__name__}")
        return body

    def _fetch_history(self) -> list[dict]:
        start_date, end_date = self._history_window()
        session = self._session_factory()
        headers = {"Authorization": f"Bearer {self._access_token}", **BASE_HEADERS}
        wagers: list[dict] = []
        for page in range(MAX_PAGES):
            response = session.get(
                HISTORY_URL, headers=headers, timeout=HTTP_TIMEOUT_SEC,
                params={"playerId": self._player_id, "startDate": start_date, "endDate": end_date,
                        "page": page, "recordsByPage": RECORDS_PER_PAGE})
            if response.status_code != 200:
                raise RuntimeError(f"BFA history page {page} failed: HTTP {response.status_code}")
            body = response.json()
            page_wagers = body.get("wagers") or []
            wagers.extend(page_wagers)
            # totalRecords counts transactions too (18 against 16 wagers live), so an
            # empty page is the other end condition.
            if not page_wagers or len(wagers) >= int(body.get("totalRecords") or 0):
                return wagers
        raise RuntimeError(f"BFA history did not end within {MAX_PAGES} pages")

    # -- Source protocol --------------------------------------------------------------

    def fetch(self) -> list[dict]:
        self._ensure_access_token()
        fetched_at = utc_now_iso()
        open_wagers = self._fetch_open_bets()
        history_wagers = self._fetch_history()
        records = merge_open_over_history(normalize_bfa(history_wagers, fetched_at),
                                          normalize_open_bets(open_wagers, fetched_at))
        n_open = sum(1 for record in records if record["status"] == "open")
        n_unmatchable = sum(1 for record in records if record["unmatchable"])
        log.info("bfa: %d open + %d history wagers -> %d records (%d open, %d unmatchable)",
                 len(open_wagers), len(history_wagers), len(records), n_open, n_unmatchable)
        return records


def source_if_configured() -> BFASource | None:
    """The source when the login is configured (the BetOnline pattern: a venue with
    no credentials reads "no source configured" in the panel rather than failing
    every poll); missing credentials are logged once with the fix."""
    if not (config.BFA_USERNAME and config.BFA_PASSWORD):
        log.warning("bfa: BFA_USERNAME / BFA_PASSWORD not set (bet_logger/.env in the main checkout) "
                    "— BFA source not registered")
        return None
    return BFASource()
