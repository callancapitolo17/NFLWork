"""Kalshi bet source: the account's fills + positions -> normalised bet records.

Three parts:
  normalize_kalshi(...)  pure port of extension/bets.js normalizeKalshi; the
                         node tests on tests/fixtures/bets/kalshi_fixture.json
                         pin its output and tests/test_parity.py holds the two
                         byte-equivalent.
  apply_net_position_pnl Python only: each ticker's P&L off its NET position,
                         after fees (`pnl`, `fees`); the Bet Tracker counts it.
  KalshiSource           the network half (Source protocol). Every poll: fills
                         since the last poll minus a 60 s overlap (trade_id
                         dedupe), unsettled positions, and one cached public
                         GET per market and per event. A full fills re-pull
                         once an hour is the reconcile; it also re-reads every
                         cached market that had no result yet, because
                         `market.result` is the only thing on a market payload
                         that changes (settlement) and it decides won/lost.

Inputs:  Kalshi REST (signed GETs through kalshi_common.auth_client — read only).
Outputs: list of records. Side effects: none on disk; in-memory fill/market/
         event caches only. Raises on a failed fills or positions pull so the
         service records a failed run and keeps the previous records.

Record conventions (plan § Kalshi specifics): one record per (ticker, side)
seen in fills; the positions endpoint is the truth for the open size (signed
position_fp, negative = NO); position_fp 0 with total_traded_dollars > 0 is
closed; fills supply the VWAP entry price and the first fill time. awayKey /
homeKey are left None — the extension resolves them through teams.js.

Subaccounts: Kalshi keeps a separate balance and set of positions per
subaccount (0 = the primary account, 1-63 numbered; the automated taker trades
its own, config.KALSHI_TAKER_SUBACCOUNT). Positions read without a subaccount
are the primary's only, so every subaccount read is pulled and normalised on
its own: fills, positions and settlements with `subaccount=N`. A numbered
subaccount's ids end in ":s<N>" (record_id_of), so the bot's bet and Cal's on
the same market stay two records; the primary's ids are unchanged.

Startup skip: a restart empties every in-memory cache, and re-reading each
market ever held (its pre-cutoff fills, the market, the event; one throttled
GET each) took ~16 minutes on 2026-10-17. A market whose stored records are
all final (FINAL_STATUSES) AND were built from a complete read (every record
stamped `raw.readComplete` == READ_VERSION: its historical fills read, its
market and event payloads present) never changes, so the first fetch skips
it: no historical read, no lookup, no record re-sent; the stored copy stands
(the store never deletes). A row from a partial read, or from an older
READ_VERSION, is never skipped, so it is rebuilt at the next restart. A
skipped market is read again once it shows an open position or a fill newer
than its stored `raw.lastFillAt` (a re-entry while the service was down
included).
"""
import logging
import threading
import time
from collections.abc import Callable, Iterator

from kalshi_common import auth_client
from unabated_ticket.bets_service import config
from unabated_ticket.bets_service.normalize import (
    cents_to_american, json_clean, parse_iso_ms, round_cents, to_number, utc_now_iso)
from unabated_ticket.bets_service.sources import kalshi_ticker as tk

log = logging.getLogger(__name__)

PAGE_LIMIT = 200
FILLS_PATH = "/portfolio/fills"
# Kalshi moves fills older than its history cutoff (2026-08-07 when first
# checked) off /portfolio/fills; the opening trades of older positions live
# only here, and the cutoff advances, so every full pull reads both.
HISTORICAL_FILLS_PATH = "/historical/fills"
# Contract counts are fractional strings ("1168.36"); below this a net
# position is flat.
FLAT_POSITION_EPSILON = 1e-6
POSITIONS_PATH = "/portfolio/positions?settlement_status=unsettled"
# Every market the account settled since Kalshi's history cutoff; a bet placed
# before the cutoff and settled after it (the Sep-15 starting-QB markets,
# bought in March-May) has fills only under /historical.
SETTLEMENTS_PATH = "/portfolio/settlements"
PRIMARY_SUBACCOUNT = 0
MAX_SUBACCOUNT = 63
# Record statuses that no later fill or settlement changes for placing or
# sizing; a market whose stored records are all one of these is skipped at
# startup (module docstring, "Startup skip").
FINAL_STATUSES = frozenset({"won", "lost", "void", "closed"})
# Stamped into raw.readComplete on a record built from a complete read. Bump it
# when a change to how records are built must reach rows already stored as
# final: every stored row then fails the match and is rebuilt at the next restart.
READ_VERSION = 3
# A record whose market or event payload could not be read (_game_fields,
# _lookup failures) carries this prefix in `unmatchable`.
UNREADABLE_REASON_PREFIX = "unreadable"


def subaccount_suffix(subaccount: int) -> str:
    """What a numbered subaccount adds to a record id; nothing for the primary."""
    return "" if subaccount == PRIMARY_SUBACCOUNT else f":s{subaccount}"


def record_id_of(ticker: str, contract_side: str, subaccount: int = PRIMARY_SUBACCOUNT) -> str:
    """The bets service's id for a Kalshi bet: kalshi:<ticker>:<side>, plus
    ":s<N>" when it is held in numbered subaccount N."""
    return f"kalshi:{ticker}:{contract_side}{subaccount_suffix(subaccount)}"


def subaccount_of_record_id(record_id: str) -> int:
    """The subaccount a Kalshi record id names (record_id_of's inverse); 0
    for the primary's three-part ids."""
    parts = record_id.split(":")
    if len(parts) == 4 and parts[3].startswith("s") and parts[3][1:].isdigit():
        return int(parts[3][1:])
    if len(parts) == 3:
        return PRIMARY_SUBACCOUNT
    raise ValueError(f"expected a Kalshi record id kalshi:<ticker>:<side>[:s<N>], got {record_id!r}")


# ---- pure normaliser (port target) ------------------------------------------------

def _fill_price_cents(fill: dict) -> float:
    key = "yes_price_dollars" if fill.get("side") == "yes" else "no_price_dollars"
    return to_number(fill.get(key)) * 100


def aggregate_fills(fills: list[dict]) -> dict:
    """VWAP in cents over buys (sells only when there were no buys), net
    contracts and the first fill time, for one (ticker, side)."""
    ordered = sorted(fills, key=lambda fill: parse_iso_ms(fill.get("created_time") or "") or 0)
    bought = bought_cost = sold = sold_cost = 0.0
    for fill in ordered:
        count = to_number(fill.get("count_fp"))
        cents = _fill_price_cents(fill)
        if fill.get("action") == "sell":
            sold += count
            sold_cost += count * cents
        else:
            bought += count
            bought_cost += count * cents
    if bought > 0:
        vwap_cents = bought_cost / bought
    elif sold > 0:
        vwap_cents = sold_cost / sold
    else:
        vwap_cents = None
    return {"vwapCents": vwap_cents, "netContracts": bought - sold,
            "firstFillAt": ordered[0].get("created_time"), "fillCount": len(ordered)}


def status_of(contract_side: str, market: dict | None, position: dict | None,
              net_contracts: float) -> str:
    result = market.get("result") if market else None
    if position:
        position_fp = to_number(position.get("position_fp"))
        if position_fp == 0 and to_number(position.get("total_traded_dollars")) > 0:
            return "closed"
        open_side = "yes" if position_fp > 0 else "no" if position_fp < 0 else None
        if open_side != contract_side:
            return "closed"
    elif net_contracts <= 0:
        return "closed"
    if result in ("yes", "no"):
        return "won" if result == contract_side else "lost"
    if result == "void":
        return "void"
    if result:
        return "unknown"
    return "open" if position or net_contracts > 0 else "unknown"


def closed_at_of(status: str, position: dict | None, market: dict | None,
                 last_fill_at: str | None) -> str | None:
    if status == "open":
        return None
    if status == "closed":
        if position and position.get("last_updated_ts") is not None:
            return position["last_updated_ts"]
        return last_fill_at
    if market and market.get("expected_expiration_time"):
        return market["expected_expiration_time"]
    return last_fill_at


def _unmatchable(base: dict, reason: str) -> dict:
    base.update({
        "league": None, "eventStart": None, "eventDate": None, "awayTeam": None, "homeTeam": None,
        "awayKey": None, "homeKey": None, "betType": "other", "period": None, "side": None,
        "points": None, "approx": [], "unmatchable": reason,
    })
    return base


def _game_fields(ticker: str, contract_side: str, market: dict | None, event: dict | None,
                 spec: dict) -> dict | str:
    """The league/teams/side fields of a game-series record, or the unmatchable reason."""
    if not market or not event:
        return f"unreadable Kalshi market (no {'event' if market else 'market'} payload for {ticker})"
    suffix = tk.parse_event_suffix(market["event_ticker"])
    event_teams = tk.parse_event_teams(event)
    if not suffix or not event_teams:
        what = "event title" if suffix else "event suffix"
        return f"unreadable Kalshi market ({what} on {market['event_ticker']})"
    strike = tk.parse_strike(ticker, market, spec["betType"], spec.get("fixed_points"),
                             event_teams["codes"])
    if isinstance(strike, str):
        return f"unreadable Kalshi market ({strike})"
    contract = tk.side_of_contract(spec["betType"], strike, contract_side, spec["league"])
    return {
        "league": spec["league"], "eventStart": suffix["eventStart"], "eventDate": suffix["eventDate"],
        "awayTeam": event_teams["awayTeam"], "homeTeam": event_teams["homeTeam"],
        "awayKey": None, "homeKey": None, "betType": spec["betType"], "period": spec["period"],
        "side": contract["side"], "points": contract["points"], "approx": contract["approx"],
        "unmatchable": None,
    }


def fill_on_its_open_side(fill: dict, position: dict | None) -> dict:
    """A sell of one side while the account is net long the OTHER side is a
    buy of that other side. Kalshi's v2 order API trades the YES leg only, so
    a NO buy (the taker bot's, unabated_ticket/bets_service/taker_orders.py)
    is placed as a YES sell and its fill may be labelled either way; the
    position's sign says which it was. Every fill carries both prices, so
    only the labels change. A sell of a side that IS held stays a sell."""
    if fill.get("action") != "sell" or not position:
        return fill
    position_fp = to_number(position.get("position_fp"))
    sold = fill.get("side")
    if sold == "yes" and position_fp < 0:
        bought = "no"
    elif sold == "no" and position_fp > 0:
        bought = "yes"
    else:
        return fill
    return {**fill, "side": bought, "outcome_side": bought, "action": "buy"}


def combo_leg_event_tickers(market: dict | None) -> list[str] | None:
    """A combo market's legs' event tickers (KXNFLSPREAD-26OCT19WASSF, ...), which
    name each leg's league; None on a market with no legs. Kalshi's GET
    /markets/<ticker> carries mve_selected_legs on finalized combos too (checked
    2026-10-20 on three KXMVECROSSCATEGORY tickers, 2026-09-09 to 2026-10-20)."""
    legs = (market or {}).get("mve_selected_legs")
    if not isinstance(legs, list) or not legs:
        return None
    return [str(leg.get("event_ticker") or "") for leg in legs]


def normalize_kalshi(fills: list[dict], positions: list[dict], markets: dict[str, dict],
                     events: dict[str, dict], fetched_at: str | None) -> list[dict]:
    """One record per (ticker, side) seen in `fills`. `markets` / `events` are
    keyed by ticker / event_ticker; a game-series market missing either fails
    closed with a specific reason. Pure."""
    positions_by_ticker = {position["ticker"]: position for position in positions}
    groups: dict[str, list[dict]] = {}
    for fill in fills:
        fill = fill_on_its_open_side(fill, positions_by_ticker.get(fill["ticker"]))
        groups.setdefault(f"{fill['ticker']}:{fill['side']}", []).append(fill)
    records = []
    for key, group in groups.items():
        ticker = group[0]["ticker"]
        contract_side = group[0]["side"]
        series = tk.series_of(ticker)
        market = markets.get(ticker)
        position = positions_by_ticker.get(ticker)
        agg = aggregate_fills(group)
        status = status_of(contract_side, market, position, agg["netContracts"])
        if position and status == "open":
            contracts = abs(to_number(position.get("position_fp")))
        else:
            contracts = max(agg["netContracts"], 0)
        entry_dollars = None if agg["vwapCents"] is None else agg["vwapCents"] / 100
        last_fill_at = max(fill.get("created_time") or "" for fill in group)
        event = events.get(market["event_ticker"]) if market else None
        base = {
            "id": f"kalshi:{key}",
            "source": "kalshi_api",
            "venue": "kalshi",
            "rotation": None,
            "price": None if agg["vwapCents"] is None else cents_to_american(agg["vwapCents"]),
            "stake": None if entry_dollars is None else round_cents(contracts * entry_dollars),
            "toWin": None if entry_dollars is None else round_cents(contracts * (1 - entry_dollars)),
            "contracts": contracts,
            "placedAt": agg["firstFillAt"],
            "status": status,
            "closedAt": closed_at_of(status, position, market, last_fill_at),
            "isParlayLeg": False, "parlayId": None, "legIndex": None, "legCount": None,
            "sourceFetchedAt": fetched_at,
            # Kalshi's own ids (#118): the event ticker's suffix is what
            # Unabated's Kalshi alt rungs carry inside sourceKey.
            "venueIds": {"marketTicker": ticker, "eventTicker": market.get("event_ticker") if market else None},
            "raw": {
                "ticker": ticker, "series": series, "contractSide": contract_side,
                "fillCount": agg["fillCount"], "vwapCents": agg["vwapCents"],
                "positionFp": to_number(position.get("position_fp")) if position else None,
                "totalTradedDollars": to_number(position.get("total_traded_dollars")) if position else None,
                "marketTitle": market.get("title") if market else None,
                "marketStatus": market.get("status") if market else None,
                "marketResult": market.get("result") if market else None,
                "eventTitle": event.get("title") if event else None,
                "legEventTickers": combo_leg_event_tickers(market),
            },
        }
        spec = tk.GAME_SERIES.get(series)
        if not spec:
            reason = tk.REASON_NOT_GAME if series in tk.NON_GAME_SERIES else tk.REASON_UNKNOWN_SERIES
            records.append(json_clean(_unmatchable(base, reason)))
            continue
        game = _game_fields(ticker, contract_side, market, event, spec)
        if isinstance(game, str):
            records.append(json_clean(_unmatchable(base, game)))
            continue
        base.update(game)
        records.append(json_clean(base))
    return records


# ---- net-position P&L --------------------------------------------------------------
#
# Kalshi holds ONE net position per market: YES and NO contracts of the same
# market cancel. A fill's `side` is the outcome it adds (a "sell no at 98c"
# while holding YES closes YES at 2c), checked against Kalshi's own
# realized_pnl_dollars to the cent on 2026-10-06. normalize_kalshi keeps one
# record per (ticker, side) for the panel; apply_net_position_pnl then prices
# each ticker off its net position so a closing trade on the other side no
# longer leaves the opening record paying out its full original size.

def _own_side_price(position_sign: int, yes_price: float) -> float:
    """What one contract of the held side cost or sold for: YES at the yes
    price, NO at one minus it."""
    return yes_price if position_sign > 0 else 1 - yes_price


def _sign(value: float) -> int:
    return 1 if value > 0 else -1


def net_position_pnl(fills: list[dict], market_result: str | None) -> dict:
    """Realized P&L of one ticker's fills, average-cost, settled at the market
    result when it has one. Pure.

    Returns {entrySide, positionYes (signed: + YES, - NO), realized (dollars,
    before fees), fees (dollars, Kalshi's fee_cost summed), final (no open
    position left: settled, or flat), heldContracts / heldCost (what was still
    held when the market settled or voided; 0 when flat)}."""
    ordered = sorted(fills, key=lambda fill: (parse_iso_ms(fill.get("created_time") or "") or 0,
                                              fill.get("trade_id") or ""))
    position_yes = cost = realized = fees = 0.0
    for fill in ordered:
        fees += to_number(fill.get("fee_cost"))
        count = to_number(fill.get("count_fp"))
        yes_price = to_number(fill.get("yes_price_dollars"))
        delta = count if fill.get("side") == "yes" else -count
        if abs(position_yes) < FLAT_POSITION_EPSILON or _sign(delta) == _sign(position_yes):
            cost += count * _own_side_price(_sign(delta), yes_price)
            position_yes += delta
            continue
        closing = min(count, abs(position_yes))
        average_cost = cost / abs(position_yes)
        realized += closing * (_own_side_price(_sign(position_yes), yes_price) - average_cost)
        cost -= closing * average_cost
        position_yes += _sign(delta) * closing
        remainder = count - closing
        if remainder > FLAT_POSITION_EPSILON:
            position_yes = _sign(delta) * remainder
            cost = remainder * _own_side_price(_sign(delta), yes_price)
    holding = abs(position_yes) >= FLAT_POSITION_EPSILON
    held_contracts = held_cost = 0.0
    if holding and market_result in ("yes", "no", "void"):
        held_contracts, held_cost = abs(position_yes), cost
    if holding and market_result in ("yes", "no"):
        held_side = "yes" if position_yes > 0 else "no"
        payout = abs(position_yes) if held_side == market_result else 0.0
        realized += payout - cost
        holding = False
    elif holding and market_result == "void":
        holding = False  # a voided market refunds the cost of what was held
    return {"entrySide": ordered[0].get("side") if ordered else None, "positionYes": position_yes,
            "realized": realized, "fees": fees, "final": not holding,
            "heldContracts": held_contracts, "heldCost": held_cost}


def grade_settled_hold(entry: dict, net: dict, market: dict | None, ticker_fills: list[dict]) -> None:
    """A record normalize_kalshi called "closed" although its position was
    still held when the market settled gets the hold's size, cost and result.
    Once a market settles, Kalshi drops it from the unsettled positions, so
    fill_on_its_open_side can no longer relabel a bot NO buy that Kalshi
    reported as a NO *sell* (seen 2026-10-16, KXNCAAFSPREAD-26OCT16CLMBPENN:
    one `no` sell fill, no buys, market result `no`); aggregate_fills then
    counts it as a sale of nothing and the record read closed with 0 contracts
    and a $0 stake. net_position_pnl, which goes by `side` alone, already
    priced it as the NO hold it was. Mutates `entry`.

    Only a ticker with no buy fills at all is regraded: nothing can be closed
    out without first being bought, so its sells were the opening trades. A
    ticker with buys may hold a real close-out labelled a sell of the held
    side (the fixture's KXNEXTTEAMNFL-26MCROSBY-LV `yes` sells), which
    `side` alone would misread as a hold."""
    if entry.get("status") != "closed" or net["heldContracts"] < FLAT_POSITION_EPSILON:
        return
    if any(fill.get("action") != "sell" for fill in ticker_fills):
        return
    held_side = "yes" if net["positionYes"] > 0 else "no"
    if held_side != entry["raw"]["contractSide"]:
        return
    result = market.get("result") if market else None
    status = "void" if result == "void" else "won" if result == held_side else "lost"
    contracts = net["heldContracts"]
    entry.update({
        "status": status, "contracts": contracts,
        "stake": round_cents(net["heldCost"]), "toWin": round_cents(contracts - net["heldCost"]),
        "closedAt": closed_at_of(status, None, market, entry.get("closedAt")),
    })


def apply_net_position_pnl(records: list[dict], fills: list[dict], markets: dict[str, dict],
                           positions: list[dict]) -> list[dict]:
    """Prices each ticker off its net position (net_position_pnl). The record
    of the side that opened the position gets `pnl` (dollars after fees, once
    final; None while a position is open) and `fees`; the other side's record,
    if any, gets `mergedInto` = that record's id and status "closed", so it
    never counts on its own. A position that flipped to the other side and is
    still open is left as normalize_kalshi built it (logged). Mutates and
    returns `records`.

    The fills are relabelled by fill_on_its_open_side against the same
    `positions` normalize_kalshi used, so a bot NO buy reported as a YES sell
    prices the NO record it was filed under (unrelabelled, it named a YES
    record that does not exist and failed the whole poll)."""
    positions_by_ticker = {position["ticker"]: position for position in positions}
    fills_by_ticker: dict[str, list[dict]] = {}
    for fill in fills:
        fill = fill_on_its_open_side(fill, positions_by_ticker.get(fill["ticker"]))
        fills_by_ticker.setdefault(fill["ticker"], []).append(fill)
    records_by_id = {record["id"]: record for record in records}
    for ticker, ticker_fills in fills_by_ticker.items():
        market = markets.get(ticker)
        net = net_position_pnl(ticker_fills, market.get("result") if market else None)
        entry_id = f"kalshi:{ticker}:{net['entrySide']}"
        other_side = "no" if net["entrySide"] == "yes" else "yes"
        other_id = f"kalshi:{ticker}:{other_side}"
        open_side = "yes" if net["positionYes"] > 0 else "no"
        if not net["final"] and open_side != net["entrySide"]:
            log.warning("kalshi %s: position flipped to %s and is open; left unmerged", ticker, open_side)
            continue
        entry = records_by_id[entry_id]
        entry["pnl"] = round_cents(net["realized"] - net["fees"]) if net["final"] else None
        entry["fees"] = round_cents(net["fees"])
        grade_settled_hold(entry, net, market, ticker_fills)
        other = records_by_id.get(other_id)
        if other is not None:
            other.update({"status": "closed", "mergedInto": entry_id, "pnl": None, "fees": 0})
    return records


def tag_subaccount(records: list[dict], subaccount: int) -> list[dict]:
    """One numbered subaccount's records (normalize_kalshi +
    apply_net_position_pnl over its own fills and positions) with ":s<N>" on
    `id` and `mergedInto` and `raw.subaccount` = N. The primary's records come
    back unchanged. Pure: returns new dicts."""
    if subaccount == PRIMARY_SUBACCOUNT:
        return records
    suffix = subaccount_suffix(subaccount)
    tagged = []
    for record in records:
        copy = {**record, "id": record["id"] + suffix, "raw": {**record["raw"], "subaccount": subaccount}}
        if record.get("mergedInto"):
            copy["mergedInto"] = record["mergedInto"] + suffix
        tagged.append(copy)
    return tagged


# ---- network half -----------------------------------------------------------------

ApiCall = Callable[[str, str], tuple[int, object, dict]]


# Integer cents (verified 2026-10-13): `balance` is cash, open positions left out;
# `portfolio_value` is what the open positions are worth now.
BALANCE_PATH = "/portfolio/balance"


def balance_of_kalshi(body: dict) -> dict:
    cash_cents = body.get("balance")
    positions_cents = body.get("portfolio_value")
    for field, value in (("balance", cash_cents), ("portfolio_value", positions_cents)):
        if isinstance(value, bool) or not isinstance(value, (int, float)):
            raise RuntimeError(f"Kalshi {BALANCE_PATH}: expected integer cents in {field}, got {value!r}")
    cash = cash_cents / 100
    return {"cash": cash, "inBets": positions_cents / 100, "available": cash, "creditLimit": None}


def _paginate(api: ApiCall, base_path: str, items_key: str) -> Iterator[dict]:
    """Walk a Kalshi list endpoint by cursor. Raises on any non-200 page so a
    partial pull is never mistaken for a complete one."""
    separator = "&" if "?" in base_path else "?"
    cursor = None
    while True:
        path = f"{base_path}{separator}limit={PAGE_LIMIT}"
        if cursor:
            path += f"&cursor={cursor}"
        status, body, _headers = api("GET", path)
        if status != 200 or not isinstance(body, dict):
            raise RuntimeError(f"GET {base_path} expected 200 with a JSON body, got {status}")
        yield from body.get(items_key) or []
        cursor = body.get("cursor")
        if not cursor:
            return


def _hold_fills(held: dict[str, dict], fills) -> int:
    """Keep each fill in `held` by trade_id; returns how many were new."""
    n_new = 0
    for fill in fills:
        trade_id = fill.get("trade_id")
        if not trade_id:
            raise RuntimeError(f"fill without trade_id on {fill.get('ticker')}: cannot dedupe")
        if trade_id not in held:
            n_new += 1
        held[trade_id] = fill
    return n_new


def _own_subaccount_fills(fills: Iterator[dict], subaccount: int) -> Iterator[dict]:
    """The fills whose `subaccount_number` is `subaccount` or unknown. A fill
    Kalshi returned under `?subaccount=N` but labelled for another subaccount
    would be filed against the wrong positions (a bot NO buy would read as a
    closed YES in the primary); it is dropped and logged instead."""
    for fill in fills:
        labelled = fill.get("subaccount_number")
        if labelled is None or labelled == subaccount:
            yield fill
            continue
        log.warning("kalshi fill %s on %s labelled subaccount %s came back under subaccount %d; dropped",
                    fill.get("trade_id"), fill.get("ticker"), labelled, subaccount)


def _validated_subaccounts(subaccounts: tuple[int, ...]) -> tuple[int, ...]:
    """The subaccounts to read: the primary first, then distinct numbered ones."""
    if not subaccounts or subaccounts[0] != PRIMARY_SUBACCOUNT:
        raise ValueError(f"expected the primary subaccount 0 first, got {subaccounts}")
    numbered = subaccounts[1:]
    if len(set(numbered)) != len(numbered) or any(not 1 <= number <= MAX_SUBACCOUNT for number in numbered):
        raise ValueError(f"expected distinct numbered subaccounts 1-{MAX_SUBACCOUNT} after 0, got {subaccounts}")
    return subaccounts


class KalshiSource:
    """Source protocol implementation for the Kalshi account (read-only GETs).

    `api` is kalshi_common.auth_client.api by default; tests inject a fake.
    `configure_auth=True` calls auth_client.configure() from config — the only
    place credentials are read; they never enter a record or a log line.
    `subaccounts` are the ones read, the primary (0) first: the service adds
    the bot's (config.KALSHI_TAKER_SUBACCOUNT). The key must be unrestricted
    (Cal's main key) to read more than one.
    `known_tickers(subaccount)` returns the tickers the store already holds
    that subaccount's Kalshi records for (BetsStore.load_kalshi_tickers);
    their pre-cutoff fills are read.
    `final_tickers(subaccount)` returns {ticker: stored lastFillAt} for the
    tickers whose stored records of that subaccount are all final and complete
    (BetsStore.load_kalshi_final_tickers); read once, on the first fetch, they
    are skipped (module docstring).
    """

    name = "kalshi"

    def __init__(self, api: ApiCall | None = None, poll_sec: float | None = None,
                 reconcile_sec: float | None = None, configure_auth: bool = True,
                 lookup_gap_sec: float | None = None, clock: Callable[[], float] = time.time,
                 known_tickers: Callable[[int], set[str]] = lambda _subaccount: set(),
                 final_tickers: Callable[[int], dict[str, str]] = lambda _subaccount: {},
                 subaccounts: tuple[int, ...] = (PRIMARY_SUBACCOUNT,)):
        self.poll_sec = poll_sec if poll_sec is not None else config.KALSHI_POLL_SEC
        self._reconcile_sec = reconcile_sec if reconcile_sec is not None else config.KALSHI_RECONCILE_SEC
        self._lookup_gap_sec = lookup_gap_sec if lookup_gap_sec is not None else config.KALSHI_LOOKUP_GAP_SEC
        self._clock = clock
        self._known_tickers = known_tickers
        self._final_tickers = final_tickers
        self._subaccounts = _validated_subaccounts(tuple(subaccounts))
        # Per subaccount: tickers whose /historical/fills have been read (pre-cutoff fills never change).
        self._historical_read: dict[int, set[str]] = {number: set() for number in self._subaccounts}
        # Per subaccount: stored-final tickers not read this process -> their stored
        # last fill (epoch ms); set on the first fetch.
        self._skipped: dict[int, dict[str, int]] | None = None
        if configure_auth:
            if not config.KALSHI_API_KEY_ID or not config.KALSHI_PRIVATE_KEY_PATH:
                raise RuntimeError("KALSHI_API_KEY_ID / KALSHI_PRIVATE_KEY_PATH not set "
                                   "(bets_service/.env, kalshi_draft/.env or the environment)")
            auth_client.configure(config.KALSHI_API_KEY_ID, config.KALSHI_PRIVATE_KEY_PATH,
                                  config.KALSHI_BASE_URL, config.PROJECT_ROOT)
        self._api = api or auth_client.api
        # Per subaccount: its fills by trade_id.
        self._fills_by_trade_id: dict[int, dict[str, dict]] = {number: {} for number in self._subaccounts}
        self._markets: dict[str, dict] = {}
        self._events: dict[str, dict] = {}
        self._last_poll_ts: float | None = None
        self._last_full_pull_ts: float | None = None
        self._last_lookup_at = 0.0
        self._lookup_lock = threading.Lock()

    # -- fills ---------------------------------------------------------------------

    def _full_pull_due(self, now: float) -> bool:
        return (self._last_full_pull_ts is None
                or now - self._last_full_pull_ts >= self._reconcile_sec)

    def _pull_fills(self, now: float) -> None:
        """Every subaccount's fills; the poll clock moves only once all of them were read."""
        full = self._full_pull_due(now)
        for subaccount in self._subaccounts:
            self._pull_subaccount_fills(subaccount, full)
        if full:
            self._last_full_pull_ts = now
            self._refresh_unsettled_markets()
        self._last_poll_ts = now

    def _pull_subaccount_fills(self, subaccount: int, full: bool) -> None:
        path = f"{FILLS_PATH}?subaccount={subaccount}"
        if not full:
            min_ts = int(self._last_poll_ts) - config.KALSHI_FILLS_OVERLAP_SEC
            path += f"&min_ts={min_ts}"
        held = self._fills_by_trade_id[subaccount]
        n_new = _hold_fills(held, _own_subaccount_fills(_paginate(self._api, path, "fills"), subaccount))
        if full:
            n_new += self._hold_historical_fills(subaccount)
        log.info("kalshi fills, subaccount %d: %s pull, %d new, %d held", subaccount,
                 "full" if full else "incremental", n_new, len(held))

    def _hold_historical_fills(self, subaccount: int) -> int:
        """The subaccount's pre-cutoff fills of markets its bets touch: one in
        its recent fills, one the store holds a record of, or one it settled
        since the cutoff (/portfolio/settlements). Read one market at a time
        (/historical/fills?ticker=), each once per process: the full history
        (~6,900 fills on 2026-10-06, mostly the MLB bots' older trades) is ~33
        pages, and on Cal's Mac a page of it hung for 34 minutes and reset.
        A market whose read fails is logged and retried on the next full pull;
        it never fails the poll."""
        held = self._fills_by_trade_id[subaccount]
        settled = _paginate(self._api, f"{SETTLEMENTS_PATH}?subaccount={subaccount}", "settlements")
        wanted = ({fill["ticker"] for fill in held.values()} | self._known_tickers(subaccount)
                  | {settlement["ticker"] for settlement in settled})
        return self._read_historical_fills(subaccount, wanted - set(self._skipped_in(subaccount)))

    def _read_historical_fills(self, subaccount: int, tickers: set[str]) -> int:
        """/historical/fills of each ticker not yet read this process."""
        held = self._fills_by_trade_id[subaccount]
        already_read = self._historical_read[subaccount]
        n_new = 0
        for ticker in sorted(tickers - already_read):
            time.sleep(self._lookup_gap_sec)
            path = f"{HISTORICAL_FILLS_PATH}?ticker={ticker}&subaccount={subaccount}"
            try:
                fills = list(_own_subaccount_fills(_paginate(self._api, path, "fills"), subaccount))
            except (RuntimeError, OSError) as error:
                log.warning("kalshi historical fills for %s (subaccount %d) failed, retried next full pull: %s",
                            ticker, subaccount, error)
                continue
            n_new += _hold_fills(held, fills)
            already_read.add(ticker)
        return n_new

    # -- startup skip (module docstring) -------------------------------------------

    def _skipped_in(self, subaccount: int) -> dict[str, int]:
        return self._skipped[subaccount] if self._skipped is not None else {}

    def _load_skipped(self) -> None:
        """On the first fetch: the stored-final tickers of every subaccount. A
        stored lastFillAt that does not parse keeps its ticker out (read in full)."""
        if self._skipped is not None:
            return
        skipped: dict[int, dict[str, int]] = {}
        for subaccount in self._subaccounts:
            last_fills = {ticker: parse_iso_ms(last_fill_at or "")
                          for ticker, last_fill_at in self._final_tickers(subaccount).items()}
            skipped[subaccount] = {ticker: ms for ticker, ms in last_fills.items() if ms is not None}
            log.info("kalshi subaccount %d: %d stored-final markets skipped at startup",
                     subaccount, len(skipped[subaccount]))
        self._skipped = skipped

    def _unskip_active(self, subaccount: int, positions: list[dict]) -> None:
        """A skipped ticker with an open position, or a fill newer than its
        stored last fill, is live again: it is read in full from now on (its
        pre-cutoff fills now; a failed read retries next full pull and stamps
        its records incomplete, so they are never skipped at a restart)."""
        skipped = self._skipped_in(subaccount)
        if not skipped:
            return
        open_tickers = {position["ticker"] for position in positions
                        if abs(to_number(position.get("position_fp"))) > FLAT_POSITION_EPSILON}
        newer_fill_tickers = {fill["ticker"] for fill in self._fills_by_trade_id[subaccount].values()
                              if fill["ticker"] in skipped
                              and (parse_iso_ms(fill.get("created_time") or "") or 0) > skipped[fill["ticker"]]}
        active = set(skipped) & (open_tickers | newer_fill_tickers)
        if not active:
            return
        log.info("kalshi subaccount %d: %d skipped markets active again: %s", subaccount, len(active), sorted(active))
        for ticker in active:
            del skipped[ticker]
        self._read_historical_fills(subaccount, active)

    def _stamp_read_state(self, records: list[dict], fills: list[dict], subaccount: int) -> list[dict]:
        """raw.lastFillAt (the market's newest fill, either side) and
        raw.readComplete (READ_VERSION when its historical fills were read, its
        market payload cached and no lookup left it unreadable, else None) on
        each record."""
        last_fill_by_ticker: dict[str, str] = {}
        for fill in fills:
            created_time = fill.get("created_time") or ""
            if created_time > last_fill_by_ticker.get(fill["ticker"], ""):
                last_fill_by_ticker[fill["ticker"]] = created_time
        historical_read = self._historical_read[subaccount]
        stamped = []
        for record in records:
            ticker = record["raw"]["ticker"]
            unreadable = str(record.get("unmatchable") or "").startswith(UNREADABLE_REASON_PREFIX)
            complete = ticker in historical_read and ticker in self._markets and not unreadable
            stamped.append({**record, "raw": {**record["raw"], "lastFillAt": last_fill_by_ticker.get(ticker),
                                              "readComplete": READ_VERSION if complete else None}})
        return stamped

    def _live_fills(self, subaccount: int) -> list[dict]:
        skipped = self._skipped_in(subaccount)
        return [fill for fill in self._fills_by_trade_id[subaccount].values() if fill["ticker"] not in skipped]

    def _pull_positions(self, subaccount: int) -> list[dict]:
        return list(_paginate(self._api, f"{POSITIONS_PATH}&subaccount={subaccount}", "market_positions"))

    # -- public market / event lookups (cached; throttled) ----------------------------

    def _lookup(self, path: str, payload_key: str) -> dict | None:
        """One throttled GET; None (logged) on a non-200 so a single bad
        ticker fails closed as an unreadable record instead of failing the poll."""
        with self._lookup_lock:
            wait = self._lookup_gap_sec - (self._clock() - self._last_lookup_at)
            if wait > 0:
                time.sleep(wait)
            status, body, _headers = self._api("GET", path)
            self._last_lookup_at = self._clock()
        if status != 200 or not isinstance(body, dict) or not isinstance(body.get(payload_key), dict):
            log.warning("kalshi lookup GET %s failed: status=%s", path, status)
            return None
        return body[payload_key]

    def _ensure_market_and_event(self, ticker: str) -> None:
        if ticker not in self._markets:
            # A market settled before Kalshi's history cutoff may answer only
            # under /historical (its fills do).
            market = (self._lookup(f"/markets/{ticker}", "market")
                      or self._lookup(f"/historical/markets/{ticker}", "market"))
            if market is None:
                return
            self._markets[ticker] = market
        event_ticker = self._markets[ticker].get("event_ticker")
        if event_ticker and event_ticker not in self._events:
            event = self._lookup(f"/events/{event_ticker}", "event")
            if event is not None:
                self._events[event_ticker] = event

    def _refresh_unsettled_markets(self) -> None:
        unsettled = [ticker for ticker, market in self._markets.items() if not market.get("result")]
        for ticker in unsettled:
            market = self._lookup(f"/markets/{ticker}", "market")
            if market is not None:
                self._markets[ticker] = market
        if unsettled:
            log.info("kalshi reconcile: re-read %d unsettled markets", len(unsettled))

    # -- Source protocol -------------------------------------------------------------

    def fetch(self) -> list[dict]:
        now = self._clock()
        self._load_skipped()
        self._pull_fills(now)
        positions = {subaccount: self._pull_positions(subaccount) for subaccount in self._subaccounts}
        for subaccount in self._subaccounts:
            self._unskip_active(subaccount, positions[subaccount])
        live_fills = {subaccount: self._live_fills(subaccount) for subaccount in self._subaccounts}
        tickers = {fill["ticker"] for fills in live_fills.values() for fill in fills}
        for ticker in sorted(tickers):
            self._ensure_market_and_event(ticker)
        fetched_at = utc_now_iso()
        records = []
        for subaccount in self._subaccounts:
            fills = live_fills[subaccount]
            subaccount_records = normalize_kalshi(fills, positions[subaccount], self._markets, self._events, fetched_at)
            subaccount_records = apply_net_position_pnl(subaccount_records, fills, self._markets, positions[subaccount])
            subaccount_records = self._stamp_read_state(subaccount_records, fills, subaccount)
            records.extend(json_clean(record) for record in tag_subaccount(subaccount_records, subaccount))
            log.info("kalshi subaccount %d: %d fills, %d positions -> %d records", subaccount, len(fills),
                     len(positions[subaccount]), len(subaccount_records))
        return records

    def balance(self) -> dict:
        """The account balance for the tracker's Balance by book card (balances.py)."""
        status, body, _headers = self._api("GET", BALANCE_PATH)
        if status != 200 or not isinstance(body, dict):
            raise RuntimeError(f"GET {BALANCE_PATH} expected 200 with a JSON body, got {status}")
        return balance_of_kalshi(body)
