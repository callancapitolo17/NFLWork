"""In-memory maker state + Kalshi reconciliation (fills, settlements, orphans).

Kalshi is the source of truth for fills and positions (MLB lesson, mirrored
from unabated_edge/maker/state.py): /portfolio/fills is polled every cycle,
/portfolio/settlements on its own cadence, and startup cancels any in-series
resting order a previous run left behind.

Exposure model (one-sided YES bids, worst case = cost):
  game_filled_cost(ticker)  — open (unsettled) fill cost for one game
  daily_committed(now)      — today's fill cost + ALL resting order cost
The engine's caps subtract the resting order under replacement themselves
(see engine.decide docstring).
"""
import datetime
import logging
from zoneinfo import ZoneInfo

from kalshi_common import auth_client

from kalshi_rfi import config, storage

log = logging.getLogger(__name__)

_ET = ZoneInfo("America/New_York")
SERIES_PREFIX = "KXMLBRFI-"


def trading_day(now: datetime.datetime) -> datetime.date:
    """Trading day rolled at DAILY_ROLL_HOUR_ET US-Eastern, not UTC midnight
    (UTC midnight is 8pm ET — mid-slate)."""
    now_et = now.astimezone(_ET)
    shifted = now_et - datetime.timedelta(hours=config.DAILY_ROLL_HOUR_ET)
    return shifted.date()


def _fp(d, key):
    v = d.get(f"{key}_fp")
    if v is None:
        v = d.get(key)
    try:
        return float(v) if v is not None else None
    except (TypeError, ValueError):
        return None


def _price_dollars(d, side):
    v = d.get(f"{side}_price_dollars")
    if v is not None:
        return float(v)
    v = d.get(f"{side}_price")
    return v / 100.0 if v is not None else None


def _money(d, key):
    v = d.get(f"{key}_dollars")
    if v is not None:
        return float(v)
    v = d.get(f"{key}_fp")
    if v is None:
        v = d.get(key)
    if v is None:
        return 0.0
    if isinstance(v, str):
        return float(v)          # string decimals are dollars
    return float(v) / 100.0      # bare numbers are integer cents


class RfiState:
    def __init__(self):
        self.resting = {}       # ticker -> {"order_id","price_cents","count"}
        self.our_orders = {}    # order_id -> ticker
        self.fills = {}         # ticker -> [(count, price_cents, day)]
        self.settled = set()    # tickers already settled
        self._fills_min_ts = None
        self._done_trades = set()
        self._order_seq = 0

    def next_client_order_id(self, now: datetime.datetime) -> str:
        self._order_seq += 1
        return f"rfi-{int(now.timestamp())}-{self._order_seq}"

    def on_place(self, ticker, order_id, price_cents, count):
        self.resting[ticker] = {"order_id": order_id,
                                "price_cents": price_cents, "count": count}
        self.our_orders[order_id] = ticker

    def on_cancel(self, ticker):
        self.resting.pop(ticker, None)

    def resting_cost_usd(self, exclude_ticker: str | None = None) -> float:
        return sum(r["price_cents"] * r["count"] / 100.0
                   for t, r in self.resting.items() if t != exclude_ticker)

    def game_filled_cost_usd(self, ticker: str) -> float:
        if ticker in self.settled:
            return 0.0
        return sum(n * p / 100.0 for (n, p, _day) in self.fills.get(ticker, []))

    def daily_filled_cost_usd(self, now: datetime.datetime) -> float:
        """Cost of today's fills, SETTLED INCLUDED — settling a market frees
        per-game headroom but must not free the daily budget (the daily cap
        bounds money put at risk per day, not money currently at risk)."""
        day = trading_day(now)
        return sum(n * p / 100.0
                   for fills in self.fills.values()
                   for (n, p, fill_day) in fills if fill_day == day)

    def daily_committed_usd(self, now: datetime.datetime,
                            exclude_ticker: str | None = None) -> float:
        return (self.daily_filled_cost_usd(now)
                + self.resting_cost_usd(exclude_ticker))

    def hydrate_from_db(self, con):
        """Startup restore from the DB, so a restart cannot forget filled
        exposure (per-game + daily caps), lose unsettled fills from
        settlement matching, or orphan fills that landed on a PREVIOUS
        run's orders while the daemon was down. Resting orders are NOT
        restored — the live orphan sweep cancels them instead, but their
        order_id→ticker map IS restored so their fills still attribute."""
        self.settled |= storage.load_settled_tickers(con)
        for order_id, ticker in storage.load_recent_placed_orders(con):
            self.our_orders.setdefault(order_id, ticker)
        for trade_id, ticker, price_cents, count, ts in \
                storage.load_recent_fills(con):
            if trade_id in self._done_trades:
                continue
            self._done_trades.add(trade_id)
            if ts.tzinfo is None:
                ts = ts.replace(tzinfo=datetime.timezone.utc)
            self.fills.setdefault(ticker, []).append(
                (float(count), int(price_cents), trading_day(ts)))
        if self.fills or self.our_orders:
            log.info("hydrated %d fills / %d prior orders from DB",
                     sum(len(v) for v in self.fills.values()),
                     len(self.our_orders))


FIRST_FILL_POLL_BACKFILL_SEC = 24 * 3600   # cover an overnight restart gap


def poll_fills(state: RfiState, con, now: datetime.datetime) -> int:
    """Ingest new fills on our orders. Returns the number ingested.
    Never raises: a network failure is a skipped poll, not a dead daemon."""
    min_ts = (state._fills_min_ts
              or int(now.timestamp()) - FIRST_FILL_POLL_BACKFILL_SEC)
    try:
        status, body, _ = auth_client.api(
            "GET", f"/portfolio/fills?limit=100&min_ts={min_ts}")
    except Exception as e:
        log.warning("poll_fills network failure: %s", e)
        return 0
    if status != 200 or not isinstance(body, dict):
        log.warning("poll_fills failed: status=%s", status)
        return 0
    n_new = 0
    for f in body.get("fills") or []:
        tid, oid = f.get("trade_id"), f.get("order_id")
        if not tid or tid in state._done_trades or oid not in state.our_orders:
            continue
        ticker = state.our_orders[oid]
        n = _fp(f, "count") or 0.0
        # We only ever place YES bids, so the YES execution price IS our price.
        yes_px = _price_dollars(f, "yes")
        price_cents = round(yes_px * 100) if yes_px is not None else None
        if price_cents is None or n <= 0:
            log.warning("fill %s unparseable payload keys=%s", tid, sorted(f))
            continue
        state._done_trades.add(tid)
        state.fills.setdefault(ticker, []).append(
            (n, price_cents, trading_day(now)))
        cur = state.resting.get(ticker)
        if cur and cur["order_id"] == oid:
            cur["count"] -= n
            if cur["count"] <= 1e-9:
                state.resting.pop(ticker, None)
        storage.log_fill(con, tid, oid, ticker, price_cents, n,
                         _money(f, "fee"))
        log.info("FILL %s yes %.0f@%dc", ticker, n, price_cents)
        n_new += 1
    # Overlap the next poll window; trade_id dedup absorbs re-delivery.
    state._fills_min_ts = int(now.timestamp()) - 60
    return n_new


def poll_settlements(state: RfiState, con):
    try:
        status, body, _ = auth_client.api(
            "GET", "/portfolio/settlements?limit=100")
    except Exception as e:
        log.warning("poll_settlements network failure: %s", e)
        return
    if status != 200 or not isinstance(body, dict):
        return
    for s in body.get("settlements") or []:
        ticker = s.get("ticker") or ""
        if (not ticker.startswith(SERIES_PREFIX) or ticker in state.settled
                or ticker not in state.fills):
            continue
        pnl = (_money(s, "revenue") - _money(s, "yes_total_cost")
               - _money(s, "no_total_cost") - _money(s, "fee_cost"))
        state.settled.add(ticker)
        storage.log_settlement(con, ticker, pnl)
        log.info("SETTLED %s pnl=%.2f", ticker, pnl)


def sweep_orphan_orders(state: RfiState, gateway) -> int:
    """Cancel in-series resting orders Kalshi knows about but we don't —
    left by a previous run or a POST that errored after landing. Runs at
    startup AND every ORPHAN_SWEEP_SEC (finding 2: a lost place-response
    orphan must not rest until an operator restart). NOTE: this cancels ANY
    of the account's KXMLBRFI orders, including manually placed ones — the
    bot owns the series while it runs (documented in README)."""
    try:
        status, body, _ = auth_client.api(
            "GET", "/portfolio/orders?status=resting&limit=1000")
    except Exception as e:
        log.warning("orphan sweep network failure: %s", e)
        return 0
    if status != 200 or not isinstance(body, dict):
        log.warning("orphan sweep: orders fetch failed status=%s", status)
        return 0
    n = 0
    for o in body.get("orders") or []:
        oid, t = o.get("order_id"), o.get("ticker") or ""
        if not oid or oid in state.our_orders or not t.startswith(SERIES_PREFIX):
            continue
        log.warning("orphan order %s on %s — cancelling", oid, t)
        if gateway.cancel(oid):
            n += 1
    return n
