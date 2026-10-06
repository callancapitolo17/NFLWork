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


def normalize_kalshi(fills: list[dict], positions: list[dict], markets: dict[str, dict],
                     events: dict[str, dict], fetched_at: str | None) -> list[dict]:
    """One record per (ticker, side) seen in `fills`. `markets` / `events` are
    keyed by ticker / event_ticker; a game-series market missing either fails
    closed with a specific reason. Pure."""
    positions_by_ticker = {position["ticker"]: position for position in positions}
    groups: dict[str, list[dict]] = {}
    for fill in fills:
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
    position left: settled, or flat)}."""
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
    if holding and market_result in ("yes", "no"):
        held_side = "yes" if position_yes > 0 else "no"
        payout = abs(position_yes) if held_side == market_result else 0.0
        realized += payout - cost
        holding = False
    elif holding and market_result == "void":
        holding = False  # a voided market refunds the cost of what was held
    return {"entrySide": ordered[0].get("side") if ordered else None, "positionYes": position_yes,
            "realized": realized, "fees": fees, "final": not holding}


def apply_net_position_pnl(records: list[dict], fills: list[dict],
                           markets: dict[str, dict]) -> list[dict]:
    """Prices each ticker off its net position (net_position_pnl). The record
    of the side that opened the position gets `pnl` (dollars after fees, once
    final; None while a position is open) and `fees`; the other side's record,
    if any, gets `mergedInto` = that record's id and status "closed", so it
    never counts on its own. A position that flipped to the other side and is
    still open is left as normalize_kalshi built it (logged). Mutates and
    returns `records`."""
    fills_by_ticker: dict[str, list[dict]] = {}
    for fill in fills:
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
        other = records_by_id.get(other_id)
        if other is not None:
            other.update({"status": "closed", "mergedInto": entry_id, "pnl": None, "fees": 0})
    return records


# ---- network half -----------------------------------------------------------------

ApiCall = Callable[[str, str], tuple[int, object, dict]]


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


class KalshiSource:
    """Source protocol implementation for the Kalshi account (read-only GETs).

    `api` is kalshi_common.auth_client.api by default; tests inject a fake.
    `configure_auth=True` calls auth_client.configure() from config — the only
    place credentials are read; they never enter a record or a log line.
    """

    name = "kalshi"

    def __init__(self, api: ApiCall | None = None, poll_sec: float | None = None,
                 reconcile_sec: float | None = None, configure_auth: bool = True,
                 lookup_gap_sec: float | None = None, clock: Callable[[], float] = time.time):
        self.poll_sec = poll_sec if poll_sec is not None else config.KALSHI_POLL_SEC
        self._reconcile_sec = reconcile_sec if reconcile_sec is not None else config.KALSHI_RECONCILE_SEC
        self._lookup_gap_sec = lookup_gap_sec if lookup_gap_sec is not None else config.KALSHI_LOOKUP_GAP_SEC
        self._clock = clock
        if configure_auth:
            if not config.KALSHI_API_KEY_ID or not config.KALSHI_PRIVATE_KEY_PATH:
                raise RuntimeError("KALSHI_API_KEY_ID / KALSHI_PRIVATE_KEY_PATH not set "
                                   "(bets_service/.env, kalshi_draft/.env or the environment)")
            auth_client.configure(config.KALSHI_API_KEY_ID, config.KALSHI_PRIVATE_KEY_PATH,
                                  config.KALSHI_BASE_URL, config.PROJECT_ROOT)
        self._api = api or auth_client.api
        self._fills_by_trade_id: dict[str, dict] = {}
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
        full = self._full_pull_due(now)
        path = FILLS_PATH
        if not full:
            min_ts = int(self._last_poll_ts) - config.KALSHI_FILLS_OVERLAP_SEC
            path = f"{FILLS_PATH}?min_ts={min_ts}"
        n_new = 0
        paths = [path, HISTORICAL_FILLS_PATH] if full else [path]
        for fill in (fill for page_path in paths for fill in _paginate(self._api, page_path, "fills")):
            trade_id = fill.get("trade_id")
            if not trade_id:
                raise RuntimeError(f"fill without trade_id on {fill.get('ticker')}: cannot dedupe")
            if trade_id not in self._fills_by_trade_id:
                n_new += 1
            self._fills_by_trade_id[trade_id] = fill
        log.info("kalshi fills: %s pull, %d new, %d held", "full" if full else "incremental",
                 n_new, len(self._fills_by_trade_id))
        if full:
            self._last_full_pull_ts = now
            self._refresh_unsettled_markets()
        self._last_poll_ts = now

    def _pull_positions(self) -> list[dict]:
        return list(_paginate(self._api, POSITIONS_PATH, "market_positions"))

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
        self._pull_fills(now)
        positions = self._pull_positions()
        fills = list(self._fills_by_trade_id.values())
        for ticker in sorted({fill["ticker"] for fill in fills}):
            self._ensure_market_and_event(ticker)
        records = normalize_kalshi(fills, positions, self._markets, self._events, utc_now_iso())
        records = [json_clean(record) for record in apply_net_position_pnl(records, fills, self._markets)]
        log.info("kalshi: %d fills, %d positions -> %d records", len(fills), len(positions),
                 len(records))
        return records
