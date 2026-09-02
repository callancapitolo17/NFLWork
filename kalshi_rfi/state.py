"""In-memory maker state + Kalshi reconciliation (fills, settlements, orphans).

Kalshi is the source of truth for fills and positions (MLB lesson, mirrored
from unabated_edge/maker/state.py): /portfolio/fills is polled every cycle,
/portfolio/settlements on its own cadence, and reconcile_resting_orders()
diffs Kalshi's resting orders against local state every cycle.

Exposure model (one-sided YES bids, worst case = cost):
  game_filled_cost(ticker)  — open (unsettled) fill cost for one game
  game_exposure_usd(ticker) — that plus every resting order on the game
  daily_committed(now)      — today's fill cost + ALL resting order cost
The engine's caps subtract the resting order under replacement themselves
(see engine.decide docstring).

Resting orders are keyed by ORDER ID, not ticker: a ticker-keyed dict
silently dropped the first order id whenever a second order landed on the
same game, orphaning a live order from every cancel, cap and sweep (bug 3,
observed 2026-08-27). The bot still intends exactly one order per game —
multiple entries mean drift, which reconcile_resting_orders() cleans up.
"""
import datetime
import logging
import time
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
        # order_id -> {"ticker","price_cents","count","exchange_index",
        #              "placed_mono"}
        self.resting = {}
        self.our_orders = {}    # order_id -> ticker
        self.fills = {}         # ticker -> [(count, price_cents, day)]
        self.settled = set()    # tickers already settled
        self.last_fill_mono = {}  # ticker -> time.monotonic() of last fill
        self.last_fill_poll_mono = 0.0
        self._fills_min_ts = None
        self._done_trades = set()
        self._order_seq = 0

    def next_client_order_id(self, now: datetime.datetime) -> str:
        self._order_seq += 1
        return f"rfi-{int(now.timestamp())}-{self._order_seq}"

    def on_place(self, ticker, order_id, price_cents, count,
                 exchange_index=None):
        self.resting[order_id] = {"ticker": ticker,
                                  "price_cents": price_cents,
                                  "count": count,
                                  "exchange_index": exchange_index,
                                  "placed_mono": time.monotonic()}
        self.our_orders[order_id] = ticker

    def on_cancel(self, order_id):
        self.resting.pop(order_id, None)

    def resting_for(self, ticker: str) -> list[tuple[str, dict]]:
        """[(order_id, record)] for one game — normally 0 or 1 entries."""
        return [(oid, r) for oid, r in self.resting.items()
                if r["ticker"] == ticker]

    def resting_tickers(self) -> list[str]:
        return sorted({r["ticker"] for r in self.resting.values()})

    def resting_cost_usd(self, exclude_ticker: str | None = None) -> float:
        return sum(r["price_cents"] * r["count"] / 100.0
                   for r in self.resting.values()
                   if r["ticker"] != exclude_ticker)

    def game_resting_cost_usd(self, ticker: str) -> float:
        return sum(r["price_cents"] * r["count"] / 100.0
                   for r in self.resting.values() if r["ticker"] == ticker)

    def game_filled_cost_usd(self, ticker: str) -> float:
        if ticker in self.settled:
            return 0.0
        return sum(n * p / 100.0 for (n, p, _day) in self.fills.get(ticker, []))

    def game_exposure_usd(self, ticker: str) -> float:
        """Worst-case dollars already committed to one game: unsettled fills
        plus every order resting on it. This is what the per-game cap must
        be checked against immediately before a placement (bug 1)."""
        return self.game_filled_cost_usd(ticker) + self.game_resting_cost_usd(ticker)

    def in_fill_cooldown(self, ticker: str, cooldown_sec: float) -> bool:
        """True while a game is cooling off after a fill. At a 1c margin the
        jump guard re-quotes every ~30s, which turned one game into 5 fills
        in 60 seconds; standing down after each fill breaks that loop and
        gives the fill poll time to register what we just bought."""
        if cooldown_sec <= 0:
            return False
        last = self.last_fill_mono.get(ticker)
        return last is not None and time.monotonic() - last < cooldown_sec

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
        restored — the live reconcile sweep cancels them instead, but their
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
        state.last_fill_mono[ticker] = time.monotonic()
        cur = state.resting.get(oid)
        if cur is not None:
            cur["count"] -= n
            if cur["count"] <= 1e-9:
                state.resting.pop(oid, None)
        storage.log_fill(con, tid, oid, ticker, price_cents, n,
                         _money(f, "fee"))
        log.info("FILL %s yes %.0f@%dc", ticker, n, price_cents)
        n_new += 1
    # Overlap the next poll window; trade_id dedup absorbs re-delivery.
    state._fills_min_ts = int(now.timestamp()) - 60
    state.last_fill_poll_mono = time.monotonic()
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


def _remaining_count(order: dict) -> float | None:
    """Contracts still resting on a Kalshi order row, or None if the payload
    doesn't say. Never guess: an unparseable count must leave local sizing
    alone rather than silently free cap headroom."""
    for key in ("remaining_count", "count"):
        v = _fp(order, key)
        if v is not None and v > 0:
            return v
    return None


RECONCILE_MIN_AGE_SEC = 5.0


def reconcile_resting_orders(state: RfiState, gateway, con=None) -> dict:
    """Diff Kalshi's resting orders against local state and converge.

    Local tracking demonstrably drifts (a cancel that 404s on the wrong
    exchange shard leaves an order live while we drop it; a place POST whose
    response was lost leaves one we never knew about), so every cycle:

      * any in-series order resting on Kalshi that we are NOT tracking is
        cancelled — the bot owns KXMLBRFI while it runs;
      * any locally-tracked order Kalshi is not resting is dropped (it
        filled, expired, or was cancelled — poll_fills owns the accounting);
      * matched orders take Kalshi's remaining count and exchange shard,
        so the caps size against reality and later cancels route correctly.

    Returns a counts dict for logging. Never raises: a failed fetch holds
    all local state (an empty listing must never read as "everything is
    gone" and blank our exposure).
    """
    try:
        status, body, _ = auth_client.api(
            "GET", "/portfolio/orders?status=resting&limit=1000")
    except Exception as e:
        log.warning("reconcile network failure: %s", e)
        return {}
    if status != 200 or not isinstance(body, dict):
        log.warning("reconcile: orders fetch failed status=%s", status)
        return {}

    live_ids = set()
    n_cancelled = n_adjusted = 0
    for o in body.get("orders") or []:
        oid, ticker = o.get("order_id"), o.get("ticker") or ""
        if not oid or not ticker.startswith(SERIES_PREFIX):
            continue
        live_ids.add(oid)
        exch = o.get("exchange_index")
        rec = state.resting.get(oid)
        if rec is None:
            log.warning("untracked resting order %s on %s — cancelling",
                        oid, ticker)
            if gateway.cancel(oid, ticker=ticker, exchange_index=exch):
                n_cancelled += 1
                if con is not None:
                    # The audit trail is the whole point: an order the bot
                    # didn't know about must still show up in `orders`.
                    storage.log_order(con, ticker, oid, "cancel", None,
                                      _remaining_count(o), "live",
                                      "reconcile_untracked")
            continue
        if exch is not None and rec.get("exchange_index") != exch:
            rec["exchange_index"] = exch
        remaining = _remaining_count(o)
        if remaining is not None and abs(remaining - rec["count"]) > 1e-9:
            log.info("reconcile %s: count %.0f -> %.0f (Kalshi)", oid,
                     rec["count"], remaining)
            rec["count"] = remaining
            n_adjusted += 1

    n_dropped = 0
    now_mono = time.monotonic()
    for oid, rec in list(state.resting.items()):
        if oid in live_ids:
            continue
        # A just-placed order may not be listed yet; only drop settled state.
        if now_mono - rec.get("placed_mono", 0.0) < RECONCILE_MIN_AGE_SEC:
            continue
        log.info("reconcile: %s on %s no longer resting — dropping",
                 oid, rec["ticker"])
        state.resting.pop(oid, None)
        n_dropped += 1

    if n_cancelled or n_dropped or n_adjusted:
        log.info("reconcile: cancelled=%d dropped=%d adjusted=%d",
                 n_cancelled, n_dropped, n_adjusted)
    return {"cancelled": n_cancelled, "dropped": n_dropped,
            "adjusted": n_adjusted}
