"""Underdog bet source: the account's prediction-market entries -> normalised bet records,
read on the service's own Underdog login (sources/underdog_auth.py).

Underdog has no public API for its prediction markets (UDX's participant API is
approval-gated), so this reads the same private endpoints Underdog's web app reads
for its Live and Results pages. Cal chose this knowing Underdog's terms ban scripted
access (2026-11-07): reads only, his own account, one poll a minute.

Two halves, like sources/prophetx.py:
  normalize_underdog(pages, fetched_at)  pure: the slips of every page read -> records.
                                         Pinned by tests/test_underdog.py.
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
Unseen until a slip settles or is cashed out: the grade fields. Every slip the active list
returns is open, whatever its status word. A settled slip is graded off its money, which
the slip states whatever the words: payout 0 lost, payout equal to contract_quantity won,
payout equal to the stake void (refunded); any other payout is a cash-out, `closed` with
`pnl` = payout - stake (the Tracker's closed-bet money); a settled slip with no payout field
is `unknown`. A slip the store holds open that neither list returns, once the settled walk
has passed the time it was placed, is `closed` with no result (the Bet105 rule).

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
    EASTERN, json_clean, mark_parse_failed, parse_iso_ms, round_cents, to_number, utc_now_iso)
from unabated_ticket.bets_service.sources.betonline import american_from_payout
from unabated_ticket.bets_service.sources.underdog_auth import UnderdogAuth

log = logging.getLogger(__name__)

SOURCE = "underdog_api"
VENUE = "underdog"
API_BASE_URL = "https://api.underdogfantasy.com"
ACTIVE_PATH = "/v10/user/active_entry_slips"
SETTLED_PATH = "/v10/user/settled_entry_slips"
# The web app's build stamp on 2026-11-07; any value was accepted, so it is a constant.
CLIENT_VERSION = "20261106182101"
HTTP_TIMEOUT_SEC = 20
# A poll reads at most this many settled pages (15 slips each), so one slip the walk
# cannot place never holds the shared poll lane for long.
MAX_SETTLED_PAGES = 20
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
    """(stake, dollars paid if it wins) of a slip: the fee, and the contracts it holds ($1 each)."""
    stake = round_cents(to_number(slip.get("fee")))
    prediction = slip.get("prediction") if isinstance(slip.get("prediction"), dict) else {}
    contracts = to_number(prediction.get("contract_quantity"))
    return stake, round_cents(contracts) if contracts > 0 else None


def status_of(slip: dict, stake: float, win_payout: float | None, on_active_list: bool) -> str:
    """The record status of a slip, graded off its money once settled (module docstring).
    Anything the active list returns is open, whatever its status word."""
    if on_active_list:
        return "open"
    if str(slip.get("status") or "").lower() != "settled" or slip.get("payout") is None:
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


def _stored_game(stored: dict | None) -> dict | None:
    """A stored record's game fields, when it had a line."""
    if not stored or stored.get("unmatchable") or stored.get("betType") not in BET_TYPES.values():
        return None
    return {field: stored.get(field) for field in GAME_FIELDS}


def record_of(slip: dict, lookups: dict, stored: dict | None, on_active_list: bool, fetched_at: str | None) -> dict:
    """One slip's record. `stored`: its stored record, if any — its game fields stand in
    when this answer no longer carries the game, and its closedAt when the slip names no
    payout time (never the poll clock twice)."""
    stake, win_payout = money_of(slip)
    to_win = round_cents(win_payout - stake) if win_payout is not None else None
    status = status_of(slip, stake, win_payout, on_active_list)
    picks = selections_of(slip)
    prediction = slip.get("prediction") if isinstance(slip.get("prediction"), dict) else {}
    first = picks[0] if picks else ({}, {})
    raw = {"slipId": slip.get("id"), "slipStatus": slip.get("status"), "fee": slip.get("fee"),
           "payout": slip.get("payout"), "contractQuantity": prediction.get("contract_quantity"),
           "openContractQuantity": prediction.get("open_contract_quantity"),
           # The game's own sport, kept for a slip that does not place (the Tracker's league column).
           "sportId": (lookups["games"].get(first[0].get("match_id")) or {}).get("sport_id"),
           "picks": [{"optionId": selection.get("option_id"), "matchId": group.get("match_id"),
                      "multiplier": selection.get("payout_multiplier"), "result": selection.get("result")}
                     for group, selection in picks]}
    closed_at = None
    if status != "open":
        closed_at = iso_utc_of(slip.get("payout_at")) or (stored or {}).get("closedAt") or fetched_at
    record = {"id": f"{VENUE}:{slip.get('id')}", "source": SOURCE, "venue": VENUE,
              "awayKey": None, "homeKey": None, "rotation": None, "contracts": win_payout,
              "isParlayLeg": False, "parlayId": None, "legIndex": None, "legCount": None,
              "sourceFetchedAt": fetched_at, "placedAt": iso_utc_of(slip.get("created_at")),
              "venueIds": {"slipId": slip.get("id")},
              "price": american_from_payout(stake, to_win) if to_win else None,
              "stake": stake, "toWin": to_win, "status": status, "approx": [], "unmatchable": None,
              "closedAt": closed_at, "raw": raw}
    if status == "closed":
        # A cash-out: the Tracker counts a closed bet's money only off `pnl`.
        record["pnl"] = round_cents(to_number(slip.get("payout")) - stake)
    record.update({field: None for field in GAME_FIELDS})
    if len(picks) != 1:
        record["betType"] = "other"
        record["unmatchable"] = f"combo of {len(picks)} picks: its picks are not matched"
        return json_clean(record)
    line = line_of(first[1], first[0], lookups)
    if isinstance(line, str) and _stored_game(stored) is not None:
        line = _stored_game(stored)
    if isinstance(line, dict):
        record.update(line)
        return json_clean(record)
    record["betType"] = "other"
    if line.startswith(("market", "league")):
        # A prop, a period market or a league the scanner does not carry: a deliberate exclusion.
        record["unmatchable"] = line
        return json_clean(record)
    return json_clean(mark_parse_failed(record, line))


def slips_of(page: object) -> tuple[list[dict], dict]:
    """(slips, lookups) of one answer; raises on an answer that is not a slip list."""
    data = page.get("data") if isinstance(page, dict) else None
    if not isinstance(data, dict) or not isinstance(data.get("entry_slips"), list):
        raise ValueError(f"underdog: answer without data.entry_slips: {str(page)[:200]}")
    return [slip for slip in data["entry_slips"] if isinstance(slip, dict)], lookups_of(data)


def normalize_underdog(active_page: dict, settled_pages: list[dict], stored: dict[str, dict],
                       fetched_at: str | None) -> list[dict]:
    """The records of every prediction slip on the pages read. A slip on both lists (it
    settled between the reads) keeps its settled copy. Pick'em slips (product "fantasy")
    are skipped. `stored`: the stored Underdog records by id. Pure."""
    by_id: dict[str, dict] = {}
    for page, on_active_list in [(active_page, True)] + [(page, False) for page in settled_pages]:
        slips, lookups = slips_of(page)
        for slip in slips:
            if slip.get("product") != PRODUCT_PREDICTION or not slip.get("id"):
                continue
            record_id = f"{VENUE}:{slip['id']}"
            by_id[record_id] = record_of(slip, lookups, stored.get(record_id), on_active_list, fetched_at)
    return list(by_id.values())


def closed_by_absence(stored: dict[str, dict], seen_ids: set[str], placed_before_ms: float | None,
                      fetched_at: str) -> list[dict]:
    """Stored open records neither list returned, once the settled walk has gone past the
    time each was placed (so it would have been seen): `closed`, no result — never
    guessed won or lost (the Bet105 rule). Pure."""
    closed = []
    for record_id, record in stored.items():
        if record.get("status") != "open" or record_id in seen_ids:
            continue
        placed_ms = parse_iso_ms(record.get("placedAt") or "")
        if placed_before_ms is None or placed_ms is None or placed_ms <= placed_before_ms:
            continue
        ended = dict(record)
        ended.update({"status": "closed", "closedAt": fetched_at, "sourceFetchedAt": fetched_at,
                      "raw": {**(record.get("raw") or {}), "closedByAbsence": True}})
        closed.append(json_clean(ended))
    return closed


# ---- network half -----------------------------------------------------------------

def _oldest_created_ms(slips: list[dict]) -> int | None:
    """Epoch ms of the oldest slip created on a settled page (they come newest first)."""
    stamps = [parse_iso_ms(slip["created_at"]) for slip in slips if isinstance(slip.get("created_at"), str)]
    stamps = [stamp for stamp in stamps if stamp is not None]
    return min(stamps) if stamps else None


class UnderdogSource:
    """Polls the account's Underdog prediction slips on the service's own login.

    Each poll reads every active slip and the newest settled page, then further
    settled pages (newest made first) while any slip the store still holds open was
    placed before the oldest slip read so far — a slip settles late, and the list is
    ordered by when it was made. The first poll after a boot also reads back
    `history_days`. A stored open slip the walk passes without seeing is closed with
    no result. At most MAX_SETTLED_PAGES pages a poll.
    """

    name = VENUE

    def __init__(self, auth: UnderdogAuth, known_records: Callable[[], dict[str, dict]] = dict,
                 poll_sec: float = config.UNDERDOG_POLL_SEC, history_days: int = config.UNDERDOG_HISTORY_DAYS,
                 session: requests.Session | None = None, clock: Callable[[], float] = time.time):
        """`known_records`: the stored Underdog records by id (store.load_records_of_venue)."""
        self.poll_sec = poll_sec
        self._auth = auth
        self._known_records = known_records
        self._history_days = history_days
        self._session = session or requests.Session()
        self._clock = clock
        self._read_history = False

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
            return reply.json()
        except ValueError as error:
            raise RuntimeError(f"underdog: GET {path} answered non-JSON: {reply.text[:200]}") from error

    def _headers(self) -> dict:
        return {"Accept": "application/json", "Client-Type": "web", "Client-Version": CLIENT_VERSION,
                "Authorization": f"Bearer {self._auth.token().token}"}

    def fetch(self) -> list[dict]:
        stored = self._known_records()
        active_page = self._get(ACTIVE_PATH)
        active_slips, _ = slips_of(active_page)
        seen_ids = {f"{VENUE}:{slip.get('id')}" for slip in active_slips}
        # Stored open slips not on the active list: placed time by id, until a settled page shows them.
        open_placed_ms = {record_id: parse_iso_ms(record.get("placedAt") or "") for record_id, record in stored.items()
                          if record.get("status") == "open" and record_id not in seen_ids}
        history_ms = None if self._read_history else (self._clock() - self._history_days * 86400) * 1000
        settled_pages: list[dict] = []
        oldest_read_ms = None
        complete = False
        for page_number in range(1, MAX_SETTLED_PAGES + 1):
            page = self._get(SETTLED_PATH, {"page": page_number})
            slips, _ = slips_of(page)
            settled_pages.append(page)
            seen_ids |= {f"{VENUE}:{slip.get('id')}" for slip in slips}
            oldest = _oldest_created_ms(slips)
            oldest_read_ms = oldest if oldest is not None else oldest_read_ms
            if not (page.get("meta") or {}).get("next") or oldest is None:
                complete = True
                break
            wanted = [stamp for record_id, stamp in open_placed_ms.items()
                      if record_id not in seen_ids and stamp is not None and stamp <= oldest]
            if not wanted and (history_ms is None or oldest < history_ms):
                break
        else:
            log.warning("underdog: stopped after %d settled pages", MAX_SETTLED_PAGES)
        self._read_history = True
        fetched_at = utc_now_iso()
        records = normalize_underdog(active_page, settled_pages, stored, fetched_at)
        # Every slip ever made was read when the list ended; else only those placed after the oldest read.
        passed_ms = -1 if complete else oldest_read_ms
        return records + closed_by_absence(stored, seen_ids, passed_ms, fetched_at)


def source_if_connected(known_records: Callable[[], dict[str, dict]] = dict) -> UnderdogSource | None:
    """An UnderdogSource when `underdog_auth connect` has been run (the token file
    exists), else None so the service leaves the venue "no source configured"."""
    if not config.UNDERDOG_TOKEN_PATH.exists():
        log.warning("underdog: no token at %s — run `python -m unabated_ticket.bets_service.sources.underdog_auth "
                    "connect`; Underdog source not registered", config.UNDERDOG_TOKEN_PATH)
        return None
    return UnderdogSource(UnderdogAuth(config.UNDERDOG_TOKEN_PATH), known_records)
