"""In-memory maker state + live-mode reconciliation against Kalshi.

Kalshi is the source of truth for fills/positions (MLB lesson): we poll
/portfolio/fills each tick and /portfolio/positions + /portfolio/settlements
each minute, and a divergence between local bookkeeping and Kalshi positions
is a tripwire, never silently reconciled."""
import datetime
import logging

from unabated_edge import config
from unabated_edge.maker import ledger, store
from kalshi_common import auth_client

log = logging.getLogger("unabated_edge")


def _fp(d, key):
    """Numeric field: prefer `<key>_fp` STRING decimal, fall back to bare int."""
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
    v = _fp(d, key)
    return v / 100.0 if v is not None else 0.0   # legacy bare value = cents


class MakerState:
    def __init__(self):
        self.resting = {}          # (ticker, side) -> {"order_id","price_cents","count"}
        self.ticker_info = {}      # ticker -> {"sport","event_id","line"}
        self.our_orders = {}       # order_id -> (ticker, side)
        self.fills = {}            # event_id -> [(line, side, contracts, price_dollars)]
        self.fills_by_ticker = {}  # ticker -> net yes contracts (positions check)
        self.fill_times = {}       # event_id -> [fill datetimes] (burst tripwire)
        self.cooloff_until = {}    # event_id -> datetime
        self.settled = set()       # tickers already settled (skip in positions check)
        self.settled_pnl_today = 0.0
        self._day = None
        self._fills_min_ts = None
        self._done_trades = set()

    def roll_day(self, now):
        if self._day != now.date():
            self._day = now.date()
            self.settled_pnl_today = 0.0

    def register_ticker(self, ticker, sport, event_id, line):
        self.ticker_info[ticker] = {"sport": sport, "event_id": event_id, "line": line}

    def resting_for(self, ticker, side):
        return self.resting.get((ticker, side))

    def on_place(self, ticker, side, order_id, price_cents, count):
        self.resting[(ticker, side)] = {"order_id": order_id, "price_cents": price_cents, "count": count}
        self.our_orders[order_id] = (ticker, side)

    def on_cancel(self, ticker, side):
        self.resting.pop((ticker, side), None)

    def tickers_for(self, eid):
        return [t for t, i in self.ticker_info.items() if i["event_id"] == eid]

    def events_with_quotes(self):
        return {self.ticker_info[t]["event_id"] for (t, _s) in self.resting if t in self.ticker_info}

    def events_with_exposure(self):
        return set(self.fills) | self.events_with_quotes()

    def quotes_live(self):
        return len(self.resting)

    def quotes_live_for(self, eid):
        return sum(1 for (t, _s) in self.resting
                   if self.ticker_info.get(t, {}).get("event_id") == eid)

    def exposure_fills(self, eid, exclude=None):
        """Committed fills PLUS all resting quotes treated as filled, so
        simultaneous fills across rungs can never breach the match cap."""
        out = list(self.fills.get(eid, []))
        for (t, s), r in self.resting.items():
            if (t, s) == exclude:
                continue
            info = self.ticker_info.get(t)
            if info and info["event_id"] == eid:
                out.append((info["line"], s, r["count"], r["price_cents"] / 100.0))
        return out


def poll_fills(state: MakerState, now) -> list[str]:
    min_ts = state._fills_min_ts or int(now.timestamp()) - 3600
    status, body, _ = auth_client.api("GET", f"/portfolio/fills?limit=100&min_ts={min_ts}")
    if status != 200 or not isinstance(body, dict):
        log.warning("maker poll_fills failed: status=%s", status)
        return []
    new = []
    for f in body.get("fills") or []:
        tid, oid = f.get("trade_id"), f.get("order_id")
        if not tid or tid in state._done_trades or oid not in state.our_orders:
            continue
        ticker, placed_side = state.our_orders[oid]
        info = state.ticker_info.get(ticker)
        if info is None:
            continue
        side = f.get("side") or placed_side
        n = _fp(f, "count") or 0.0
        price = _price_dollars(f, side)
        if price is None or n <= 0:
            log.warning("maker fill %s unparseable payload keys=%s", tid, sorted(f))
            continue
        eid = info["event_id"]
        state.fills.setdefault(eid, []).append((info["line"], side, n, price))
        state.fills_by_ticker[ticker] = state.fills_by_ticker.get(ticker, 0.0) + (n if side == "yes" else -n)
        state.fill_times.setdefault(eid, []).append(now)
        state._done_trades.add(tid)
        cur = state.resting.get((ticker, placed_side))
        if cur and cur["order_id"] == oid:
            cur["count"] -= n
            if cur["count"] <= 1e-9:
                state.resting.pop((ticker, placed_side), None)
        worst = ledger.worst_case(state.fills[eid])
        store.log_fill(now, info["sport"], eid, oid, ticker, side, price, n,
                       _fp(f, "fee") or 0.0, worst, tid)
        log.info("maker FILL %s %s %.2f@%.2f worst_after=%.2f", ticker, side, n, price, worst)
        new.append(tid)
    # overlap the next window; trade_id dedup absorbs the re-delivery
    state._fills_min_ts = int(now.timestamp()) - 60
    return new


def poll_positions(state: MakerState) -> bool:
    """True when Kalshi's net positions match our fills-derived book on every
    ticker we quoted (settled tickers excluded). API failure returns True —
    'cannot verify' must not false-trip the tripwire."""
    if not state.fills_by_ticker:
        return True
    status, body, _ = auth_client.api("GET", "/portfolio/positions?limit=1000")
    if status != 200 or not isinstance(body, dict):
        log.warning("maker poll_positions failed: status=%s", status)
        return True
    kalshi_pos = {p.get("ticker"): (_fp(p, "position") or 0.0)
                  for p in body.get("market_positions") or []}
    for ticker, expected in state.fills_by_ticker.items():
        if ticker in state.settled:
            continue
        actual = kalshi_pos.get(ticker, 0.0)
        if abs(actual - expected) > 0.01:
            log.warning("maker position mismatch %s: kalshi=%.2f local=%.2f", ticker, actual, expected)
            return False
    return True


def poll_settlements(state: MakerState, now):
    status, body, _ = auth_client.api("GET", "/portfolio/settlements?limit=100")
    if status != 200 or not isinstance(body, dict):
        return
    for s in body.get("settlements") or []:
        ticker = s.get("ticker")
        info = state.ticker_info.get(ticker)
        if info is None or ticker in state.settled:
            continue
        pnl = _money(s, "revenue") - _money(s, "yes_total_cost") - _money(s, "no_total_cost")
        state.settled.add(ticker)
        state.fills_by_ticker.pop(ticker, None)
        state.roll_day(now)
        state.settled_pnl_today += pnl
        store.log_settlement(now, info["sport"], ticker, pnl)
        log.info("maker SETTLED %s pnl=%.2f day_total=%.2f", ticker, pnl, state.settled_pnl_today)
