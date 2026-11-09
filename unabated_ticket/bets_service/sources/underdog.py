"""Underdog bet source: the account's prediction-market slips -> normalised bet records,
PUSHED by the extension (POST /underdog.json), which reads them on the Underdog web
app's own login in Cal's Chrome (extension/underdog.js, underdogpage.js).

Underdog has no public API for its prediction markets (UDX's participant API is
approval-gated), so the extension reads the same private endpoints the web app reads
for its Live and Results pages, with the headers the app itself just sent. Cal chose
this knowing Underdog's terms ban scripted access (2026-11-07): reads only, his own
account. Not from the VM: Underdog's login server answers a script with a Cloudflare
bot wall (2026-11-09, refresh grants too), and getting past one is not built.

Pure: validate_push + records_of_push turn one push into records; pinned by
tests/test_underdog.py.

Endpoints (captured from Cal's logged-in Chrome, 2026-11-07 / 11-09):
  GET api.underdogfantasy.com/v10/user/active_entry_slips             every open slip
  GET api.underdogfantasy.com/v10/user/settled_entry_slips?page=N     15 per page, newest
                                                                      created first; meta
                                                                      {count, items, next, page}
Each answer is {data: {entry_slips, games, teams, over_under_lines, over_under_options,
over_unders, appearances, ...}, meta}: the slips plus every game, team, line and option
they name, so each page is read on its own lookups. A settled page is ~600 KB whole, so
the extension sends only the prediction slips and the lookups they name.

Wire grammar (Cal's $1 49ers prediction slip, open 2026-11-07 and settled 2026-11-09):
  slip    {id, product ("prediction_markets"; "fantasy" = pick'em, skipped), type
          ("Prediction"), status ("active" | "settled"), display_result ("won"), fee
          (stake, dollars), payout (dollars paid at settlement: "1.23" on the win),
          created_at, payout_at, prediction {contract_quantity (= dollars paid if it
          wins), open_contract_quantity, trade_cost, total_fees}, selection_groups
          [{match_id, match_type "Game", selections [{option_id, payout_multiplier,
          result ("pending" | "won")}]}]}
  option  (over_under_options) {id, over_under_line_id, choice ("home" | "away" on a
          moneyline or spread, "higher" | "lower" on a total), choice_display
          ("49ers to win", "CIN -3.5", "Higher"), american_price, decimal_price}
  line    (over_under_lines) {id, over_under_id, provider_id "udx", stat_value (the
          total, or the spread's size unsigned: "3.5"; null on a moneyline)}
  over_under {id, title ("LV @ SF Moneyline", "CIN @ ATL Spread", "CIN @ ATL Total
          Points O/U"), appearance_stat {stat: "moneyline" | "spread" | "points"}}
  game    {id, sport_id ("NFL"), away_team_id, home_team_id, scheduled_at (UTC)}
  team    {id, name ("San Francisco 49ers")}
Every slip the active list returns is open, whatever its status word. A settled slip is
graded off its money, which the slip states whatever the words: payout 0 lost, payout
equal to contract_quantity won, payout equal to the stake void (refunded); any other
payout is a cash-out, `closed` with `pnl` = payout - stake (the Tracker's closed-bet
money); a settled slip with no payout field is `unknown`. A slip the store holds open
that neither list returns, once the settled pages have passed the time it was placed,
is `closed` with no result (the Bet105 rule).

Record conventions: one record per slip, id "underdog:<slip id>"; stake = fee (what left
the wallet), toWin = contract_quantity - fee, price = the American price of that ratio
(Underdog's 1.23x on $1 -> -435). A slip of more than one pick (a combo) is one
unmatchable record with its picks in raw. awayKey / homeKey stay None — the extension
resolves the names through teams.js.

Outputs: records. Side effects: none (the service stores them).
"""
import re
from datetime import datetime, timezone

from unabated_ticket.bets_service.normalize import (
    EASTERN, json_clean, mark_parse_failed, parse_iso_ms, round_cents, to_number, utc_now_iso)
from unabated_ticket.bets_service.sources.betonline import american_from_payout

SOURCE = "underdog_api"
VENUE = "underdog"
# A push carries at most this many settled pages (15 slips each; the extension's own
# cap, extension/underdog.js MAX_SETTLED_PAGES) and slips.
MAX_SETTLED_PAGES = 20
MAX_SLIPS_PER_PUSH = 500
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


# ---- one push -----------------------------------------------------------------------

def _validate_page(name: str, page: object) -> str | None:
    if not isinstance(page, dict) or not isinstance(page.get("data"), dict):
        return f"{name} must be an answer object with data"
    data = page["data"]
    if not isinstance(data.get("entry_slips"), list) or not all(isinstance(slip, dict) for slip in data["entry_slips"]):
        return f"{name}.data.entry_slips must be a list of objects"
    for key in ("over_under_options", "over_under_lines", "over_unders", "games", "teams"):
        if key in data and not isinstance(data[key], list):
            return f"{name}.data.{key} must be a list, got {type(data[key]).__name__}"
    return None


def validate_push(body: object) -> dict | str:
    """{fetchedAt, active, settled, complete} or {error} of a POST /underdog.json body, or
    what is wrong with it. A push is the whole active list plus the settled pages read from
    page 1, or it is an error — never half (a half push would close open slips)."""
    if not isinstance(body, dict):
        return f"body must be an object, got {type(body).__name__}"
    if "error" in body:
        if not isinstance(body["error"], str) or not body["error"].strip():
            return "error must be a non-empty string"
        return {"error": body["error"].strip()}
    fetched_at = body.get("fetchedAt")
    if not isinstance(fetched_at, str) or parse_iso_ms(fetched_at) is None:
        return f"fetchedAt must be an ISO time, got {fetched_at!r}"
    problem = _validate_page("active", body.get("active"))
    if problem:
        return problem
    settled = body.get("settled")
    if not isinstance(settled, list) or not 1 <= len(settled) <= MAX_SETTLED_PAGES:
        return f"settled must be a list of 1 to {MAX_SETTLED_PAGES} pages"
    for index, page in enumerate(settled):
        problem = _validate_page(f"settled[{index}]", page)
        if problem:
            return problem
    if not isinstance(body.get("complete"), bool):
        return f"complete must be true or false, got {body.get('complete')!r}"
    slip_count = sum(len(page["data"]["entry_slips"]) for page in [body["active"]] + settled)
    if slip_count > MAX_SLIPS_PER_PUSH:
        return f"at most {MAX_SLIPS_PER_PUSH} slips per push, got {slip_count}"
    return {"fetchedAt": fetched_at, "active": body["active"], "settled": settled, "complete": body["complete"]}


def _oldest_created_ms(slips: list[dict]) -> int | None:
    stamps = [parse_iso_ms(slip["created_at"]) for slip in slips if isinstance(slip.get("created_at"), str)]
    stamps = [stamp for stamp in stamps if stamp is not None]
    return min(stamps) if stamps else None


def records_of_push(push: dict, stored: dict[str, dict]) -> tuple[list[dict], list[dict]]:
    """(records, closed) of a validated push. `stored`: the stored Underdog records by id.
    `complete` (the settled list ended) closes every stored open slip neither list
    returned; otherwise only those placed after the oldest slip the pages carry. Pure."""
    fetched_at = push["fetchedAt"]
    records = normalize_underdog(push["active"], push["settled"], stored, fetched_at)
    pages = [push["active"]] + push["settled"]
    seen_ids = {f"{VENUE}:{slip.get('id')}" for page in pages for slip in slips_of(page)[0]}
    settled_slips = [slip for page in push["settled"] for slip in slips_of(page)[0]]
    passed_ms = -1 if push["complete"] else _oldest_created_ms(settled_slips)
    return records, closed_by_absence(stored, seen_ids, passed_ms, fetched_at)
