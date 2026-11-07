"""Underdog bet source: the account's prediction-market entries -> normalised bet records,
read on the service's own Underdog login (sources/underdog_auth.py).

Underdog has no public API for its prediction markets (UDX's participant API is
approval-gated), so this reads the same private endpoints Underdog's web app reads
for its Live and Results pages. Cal chose this knowing Underdog's terms ban scripted
access (2026-11-07): reads only, his own account, one poll a minute.

Two halves, like sources/prophetx.py:
  normalize_underdog(pages, fetched_at)  pure: the slips of every page read -> records.
                                         Pinned by tests/test_underdog.py on
                                         tests/fixtures/bets/underdog_slips.json.
  UnderdogSource                         the network half (Source protocol).

Endpoints (captured from Cal's logged-in Chrome, 2026-11-07; only the bearer token is
required, no location headers — tested both ways):
  GET api.underdogfantasy.com/v10/user/active_entry_slips             every open slip
  GET api.underdogfantasy.com/v10/user/settled_entry_slips?page=N     15 per page, newest
                                                                      created first; meta
                                                                      {count, items, next, page}
Each answer is {data: {entry_slips, games, teams, over_under_lines, over_under_options,
over_unders, appearances, ...}, meta}: the slips plus every game, team, line and option
they name, so each page is read on its own lookups.

Wire grammar (one live $1 prediction slip, 2026-11-07, plus the app's NFL board):
  slip    {id, product ("prediction_markets"; "fantasy" = pick'em, skipped), type
          ("Prediction"), status ("active" | "settled"), fee (stake, dollars),
          payout (dollars paid at settlement), created_at, payout_at,
          prediction {contract_quantity (= dollars paid if it wins), trade_cost,
          total_fees}, selection_groups [{match_id, match_type "Game", selections
          [{option_id, payout_multiplier, result ("pending")}]}]}
  option  (over_under_options) {id, over_under_line_id, choice ("home" | "away" on a
          moneyline or spread, "higher" | "lower" on a total), choice_display
          ("49ers to win", "CIN -3.5", "Higher"), american_price, decimal_price}
  line    (over_under_lines) {id, over_under_id, provider_id "udx", stat_value (the
          total, or the spread's size unsigned: "3.5")}
  over_under {id, title ("LV @ SF Moneyline", "CIN @ ATL Spread", "CIN @ ATL Total
          Points O/U"), appearance_stat {stat: "moneyline" | "spread" | "points"}}
  game    {id, sport_id ("NFL"), away_team_id, home_team_id, scheduled_at (UTC)}
  team    {id, name ("San Francisco 49ers")}
Unseen until a slip settles or is cashed out: the grade fields. A settled slip is graded
off its money, which the slip states whatever the words: payout 0 lost, payout equal to
contract_quantity won, payout equal to the stake void (refunded); any other payout is a
cash-out, `closed` with what it paid in raw (the Polymarket US rule for a side sold back).

Record conventions: one record per slip, id "underdog:<slip id>"; stake = fee (what left
the wallet), toWin = contract_quantity - fee, price = the American price of that ratio
(Underdog's 1.23x on $1 -> -435). A slip of more than one pick (a combo) is one
unmatchable record with its picks in raw. awayKey / homeKey stay None — the extension
resolves the names through teams.js.

Outputs: records. Side effects: network reads only. Raises on a failed or malformed read
(the service keeps the previous records and logs the failed run).
"""
import logging
import re
import time
from collections.abc import Callable
from datetime import datetime, timezone

import requests

from unabated_ticket.bets_service import config
from unabated_ticket.bets_service.normalize import (
    EASTERN, json_clean, mark_parse_failed, round_cents, to_number, utc_now_iso)
from unabated_ticket.bets_service.sources.betonline import american_from_payout
from unabated_ticket.bets_service.sources.underdog_auth import UnderdogAuth

log = logging.getLogger(__name__)

SOURCE = "underdog_api"
VENUE = "underdog"
API_BASE_URL = "https://api.underdogfantasy.com"
ACTIVE_PATH = "/v10/user/active_entry_slips"
SETTLED_PATH = "/v10/user/settled_entry_slips"
# The web app's own build stamp; any value was accepted on 2026-11-07, so it is only
# read once per boot to look like the current client.
CLIENT_VERSION_URL = "https://app.underdogsports.com/client-version"
FALLBACK_CLIENT_VERSION = "20261106182101"
HTTP_TIMEOUT_SEC = 20
# 1,726 settled pick'em slips at 15 a page is 116 pages; a cursor that never ends stops here.
MAX_SETTLED_PAGES = 150
PRODUCT_PREDICTION = "prediction_markets"
# Underdog sport_id -> the extension's league path (feed.LEAGUES).
LEAGUES = {"NFL": "nfl", "CFB": "cfb", "NBA": "nba", "CBB": "cbb", "WNBA": "wnba", "MLB": "mlb", "NHL": "nhl"}
BET_TYPES = {"moneyline": "moneyline", "spread": "spread", "points": "total"}
TEAM_CHOICES = {"away", "home"}
TOTAL_CHOICES = {"higher": "over", "lower": "under"}
# A full-game market's title; a half or quarter market would name its period, and none
# has been seen, so one fails closed rather than being filed as the full game.
PERIOD_WORDS_RE = re.compile(r"\b(1H|2H|1Q|2Q|3Q|4Q|half|quarter|period|inning|F5)\b", re.IGNORECASE)
SIGNED_NUMBER_RE = re.compile(r"([+-]\d+(?:\.\d+)?)\s*$")
# A payout within a cent of a reference amount is that amount.
MONEY_TOLERANCE = 0.011


def iso_utc_of(value: object) -> str | None:
    if not isinstance(value, str) or not value.strip():
        return None
    try:
        parsed = datetime.fromisoformat(value.strip().replace("Z", "+00:00"))
    except ValueError:
        return None
    if parsed.tzinfo is None:
        parsed = parsed.replace(tzinfo=timezone.utc)
    return parsed.astimezone(timezone.utc).strftime("%Y-%m-%dT%H:%M:%SZ")


def eastern_date_of(iso: str | None) -> str | None:
    if not iso:
        return None
    return datetime.fromisoformat(iso.replace("Z", "+00:00")).astimezone(EASTERN).date().isoformat()


def _by_id(items: object) -> dict:
    """Underdog sends each lookup as a list (slip answers) or an id-keyed dict (board answers)."""
    if isinstance(items, dict):
        items = list(items.values())
    if not isinstance(items, list):
        return {}
    return {item["id"]: item for item in items if isinstance(item, dict) and item.get("id") is not None}


def selections_of(slip: dict) -> list[tuple[dict, dict]]:
    """(selection group, selection) pairs of a slip, in order."""
    pairs = []
    for group in slip.get("selection_groups") or []:
        for selection in group.get("selections") or []:
            pairs.append((group, selection))
    return pairs


def money_of(slip: dict) -> tuple[float, float | None]:
    """(stake, dollars paid if it wins) of a slip: the fee, and the contracts it holds
    ($1 each), else the fee times the slip's own payout multiplier."""
    stake = round_cents(to_number(slip.get("fee")))
    prediction = slip.get("prediction") if isinstance(slip.get("prediction"), dict) else {}
    contracts = to_number(prediction.get("contract_quantity"))
    if contracts > 0:
        return stake, round_cents(contracts)
    multiplier = to_number(slip.get("initial_max_payout_multiplier"))
    return stake, round_cents(stake * multiplier) if multiplier > 0 else None


def status_of(slip: dict, stake: float, win_payout: float | None) -> str:
    """The record status of a slip, graded off its money once settled (module docstring)."""
    status = str(slip.get("status") or "").lower()
    if status == "active":
        return "open"
    if status in ("cancelled", "canceled", "refunded", "voided"):
        return "void"
    if status != "settled":
        return "unknown"
    paid = to_number(slip.get("payout"))
    if paid <= MONEY_TOLERANCE:
        return "lost"
    if win_payout is not None and abs(paid - win_payout) <= MONEY_TOLERANCE:
        return "won"
    if abs(paid - stake) <= MONEY_TOLERANCE:
        return "void"
    return "closed"


def line_of(selection: dict, group: dict, lookups: dict) -> dict | str:
    """The board line a single pick is on, or why it has none."""
    option = lookups["options"].get(selection.get("option_id"))
    if option is None:
        return f"option {selection.get('option_id')} not in the answer"
    line = lookups["lines"].get(option.get("over_under_line_id"))
    over_under = lookups["over_unders"].get(line.get("over_under_id")) if line else None
    if line is None or over_under is None:
        return f"line {option.get('over_under_line_id')} not in the answer"
    stat = (over_under.get("appearance_stat") or {}).get("stat")
    bet_type = BET_TYPES.get(stat)
    if bet_type is None:
        return f"market {stat!r} ({over_under.get('title')}) is not a game line"
    title = str(over_under.get("title") or "")
    if PERIOD_WORDS_RE.search(title):
        return f"market {title!r} names a period; only full-game lines are read"
    game = lookups["games"].get(group.get("match_id"))
    if game is None or group.get("match_type") != "Game":
        return f"game {group.get('match_id')} not in the answer"
    league = LEAGUES.get(str(game.get("sport_id") or ""))
    if league is None:
        return f"league {game.get('sport_id')!r} not read"
    away = lookups["teams"].get(game.get("away_team_id"))
    home = lookups["teams"].get(game.get("home_team_id"))
    if not away or not home or not away.get("name") or not home.get("name"):
        return f"game {game.get('id')}: a team is not in the answer"
    side_points = side_and_points_of(bet_type, option, line)
    if isinstance(side_points, str):
        return side_points
    start = iso_utc_of(game.get("scheduled_at"))
    return {"league": league, "eventStart": start, "eventDate": eastern_date_of(start),
            "awayTeam": away["name"], "homeTeam": home["name"], "betType": bet_type, "period": "FG",
            **side_points}


def side_and_points_of(bet_type: str, option: dict, line: dict) -> dict | str:
    choice = str(option.get("choice") or "").lower()
    if bet_type == "total":
        side = TOTAL_CHOICES.get(choice)
        points = to_number(line.get("stat_value"))
        if side is None or points <= 0:
            return f"total pick {choice!r} at {line.get('stat_value')!r} does not read"
        return {"side": side, "points": points}
    if choice not in TEAM_CHOICES:
        return f"{bet_type} pick {choice!r} names no team side"
    if bet_type == "moneyline":
        return {"side": choice, "points": None}
    # The spread's sign is only in the pick's text ("CIN -3.5"); the line holds its size.
    signed = SIGNED_NUMBER_RE.search(str(option.get("choice_display") or ""))
    size = to_number(line.get("stat_value"))
    if signed is None or abs(abs(float(signed.group(1))) - size) > 1e-9:
        return f"spread pick {option.get('choice_display')!r} does not match the line {line.get('stat_value')!r}"
    return {"side": choice, "points": float(signed.group(1))}


def lookups_of(data: dict) -> dict:
    return {"options": _by_id(data.get("over_under_options")), "lines": _by_id(data.get("over_under_lines")),
            "over_unders": _by_id(data.get("over_unders")), "games": _by_id(data.get("games")),
            "teams": _by_id(data.get("teams"))}


GAME_FIELDS = ("league", "eventStart", "eventDate", "awayTeam", "homeTeam", "betType", "period", "side", "points")


def record_of(slip: dict, lookups: dict, fetched_at: str | None) -> dict:
    stake, win_payout = money_of(slip)
    to_win = round_cents(win_payout - stake) if win_payout is not None else None
    status = status_of(slip, stake, win_payout)
    picks = selections_of(slip)
    prediction = slip.get("prediction") if isinstance(slip.get("prediction"), dict) else {}
    raw = {"slipId": slip.get("id"), "type": slip.get("type"), "slipStatus": slip.get("status"),
           "displayResult": slip.get("display_result"), "fee": slip.get("fee"), "payout": slip.get("payout"),
           "contractQuantity": prediction.get("contract_quantity"), "tradeCost": prediction.get("trade_cost"),
           "totalFees": prediction.get("total_fees"), "freeEntry": slip.get("free_entry"),
           "usedBonusFunds": slip.get("used_bonus_funds"),
           "picks": [{"optionId": selection.get("option_id"), "matchId": group.get("match_id"),
                      "multiplier": selection.get("payout_multiplier"), "result": selection.get("result")}
                     for group, selection in picks]}
    first = picks[0] if picks else ({}, {})
    # The game's own sport, kept for a slip that does not place (the Tracker's league column).
    raw["sportId"] = (lookups["games"].get(first[0].get("match_id")) or {}).get("sport_id")
    record = {"id": f"{VENUE}:{slip.get('id')}", "source": SOURCE, "venue": VENUE,
              "awayKey": None, "homeKey": None, "rotation": None, "contracts": win_payout,
              "isParlayLeg": False, "parlayId": None, "legIndex": None, "legCount": None,
              "sourceFetchedAt": fetched_at, "placedAt": iso_utc_of(slip.get("created_at")),
              "venueIds": {"slipId": slip.get("id"), "optionId": first[1].get("option_id"),
                           "matchId": first[0].get("match_id")},
              "price": american_from_payout(stake, to_win) if to_win else None,
              "stake": stake, "toWin": to_win, "status": status, "approx": [], "unmatchable": None,
              "closedAt": None if status == "open" else iso_utc_of(slip.get("payout_at")) or fetched_at,
              "raw": raw}
    record.update({field: None for field in GAME_FIELDS})
    if len(picks) != 1:
        record["betType"] = "other"
        record["unmatchable"] = f"combo of {len(picks)} picks: its picks are not matched"
        return json_clean(record)
    line = line_of(first[1], first[0], lookups)
    if isinstance(line, dict):
        record.update(line)
        return json_clean(record)
    record["betType"] = "other"
    if line.startswith(("market", "league")):
        # A prop, a period market or a league the scanner does not carry: a deliberate exclusion.
        record["unmatchable"] = line
        return json_clean(record)
    return json_clean(mark_parse_failed(record, line))


def normalize_underdog(pages: list[dict], fetched_at: str | None) -> list[dict]:
    """The records of every prediction slip on the pages read (active first, then
    settled). A slip on two pages (it settled between the reads) keeps its last copy,
    the settled one. Pick'em slips (product "fantasy") are skipped. Pure."""
    by_id: dict[str, dict] = {}
    for page in pages:
        data = page.get("data") if isinstance(page, dict) else None
        if not isinstance(data, dict):
            raise ValueError(f"underdog: answer without a data object: {str(page)[:200]}")
        lookups = lookups_of(data)
        for slip in data.get("entry_slips") or []:
            if not isinstance(slip, dict) or slip.get("product") != PRODUCT_PREDICTION or not slip.get("id"):
                continue
            by_id[str(slip["id"])] = record_of(slip, lookups, fetched_at)
    return list(by_id.values())


# ---- network half -----------------------------------------------------------------

def _oldest_created_ms(page: dict) -> float | None:
    """Epoch ms of the oldest slip created on a settled page (they come newest first)."""
    slips = (page.get("data") or {}).get("entry_slips") or []
    stamps = [iso_utc_of(slip.get("created_at")) for slip in slips if isinstance(slip, dict)]
    stamps = [stamp for stamp in stamps if stamp]
    if not stamps:
        return None
    return datetime.fromisoformat(min(stamps).replace("Z", "+00:00")).timestamp() * 1000


class UnderdogSource:
    """Polls the account's Underdog prediction slips on the service's own login.

    Each poll reads every active slip, then settled pages newest first until a page
    holds only slips created more than `history_days` ago AND every slip the store
    still holds as open has been seen (a slip settles late, and its page is ordered by
    when it was created, not settled) — at most MAX_SETTLED_PAGES pages.
    """

    name = VENUE

    def __init__(self, auth: UnderdogAuth, known_open_ids: Callable[[], set[str]] = set,
                 poll_sec: float = config.UNDERDOG_POLL_SEC, history_days: int = config.UNDERDOG_HISTORY_DAYS,
                 session: requests.Session | None = None, clock: Callable[[], float] = time.time):
        """`known_open_ids`: the slip ids of the stored Underdog records still open."""
        self.poll_sec = poll_sec
        self._auth = auth
        self._known_open_ids = known_open_ids
        self._history_days = history_days
        self._session = session or requests.Session()
        self._clock = clock
        self._client_version: str | None = None

    def _headers(self) -> dict:
        if self._client_version is None:
            self._client_version = self._read_client_version()
        return {"Accept": "application/json", "Client-Type": "web", "Client-Version": self._client_version,
                "Authorization": f"Bearer {self._auth.token().token}"}

    def _read_client_version(self) -> str:
        try:
            reply = self._session.get(CLIENT_VERSION_URL, timeout=HTTP_TIMEOUT_SEC)
            version = reply.json().get("clientVersion") if reply.status_code == 200 else None
        except (requests.RequestException, ValueError):
            version = None
        return str(version) if version else FALLBACK_CLIENT_VERSION

    def _get(self, path: str, params: dict | None = None) -> dict:
        reply = self._session.get(API_BASE_URL + path, params=params, headers=self._headers(), timeout=HTTP_TIMEOUT_SEC)
        if reply.status_code == 401:
            # A token Underdog refused before its expiry: refresh once and retry.
            self._auth.forget()
            reply = self._session.get(API_BASE_URL + path, params=params, headers=self._headers(),
                                      timeout=HTTP_TIMEOUT_SEC)
        if reply.status_code != 200:
            raise RuntimeError(f"underdog: GET {path} {params or ''} answered HTTP {reply.status_code}: {reply.text[:200]}")
        try:
            payload = reply.json()
        except ValueError as error:
            raise RuntimeError(f"underdog: GET {path} answered non-JSON: {reply.text[:200]}") from error
        if not isinstance(payload, dict):
            raise RuntimeError(f"underdog: GET {path} answered {type(payload).__name__}, expected an object")
        return payload

    def fetch(self) -> list[dict]:
        pages = [self._get(ACTIVE_PATH)]
        seen = {str(slip.get("id")) for slip in (pages[0].get("data") or {}).get("entry_slips") or []
                if isinstance(slip, dict)}
        still_open = set(self._known_open_ids()) - seen
        cutoff_ms = (self._clock() - self._history_days * 86400) * 1000
        page_number = 1
        while page_number <= MAX_SETTLED_PAGES:
            page = self._get(SETTLED_PATH, {"page": page_number})
            pages.append(page)
            still_open -= {str(slip.get("id")) for slip in (page.get("data") or {}).get("entry_slips") or []
                           if isinstance(slip, dict)}
            oldest = _oldest_created_ms(page)
            next_page = (page.get("meta") or {}).get("next")
            if not next_page or oldest is None:
                break
            if oldest < cutoff_ms and not still_open:
                break
            page_number = int(next_page)
        else:
            log.warning("underdog: stopped after %d settled pages; %d open slip(s) not found",
                        MAX_SETTLED_PAGES, len(still_open))
        return normalize_underdog(pages, utc_now_iso())


def source_if_connected(known_open_ids: Callable[[], set[str]] = set) -> UnderdogSource | None:
    """An UnderdogSource when `underdog_auth connect` has been run (the token file
    exists), else None so the service leaves the venue "no source configured"."""
    if not config.UNDERDOG_TOKEN_PATH.exists():
        log.warning("underdog: no token at %s — run `python -m unabated_ticket.bets_service.sources.underdog_auth "
                    "connect`; Underdog source not registered", config.UNDERDOG_TOKEN_PATH)
        return None
    return UnderdogSource(UnderdogAuth(config.UNDERDOG_TOKEN_PATH), known_open_ids)
