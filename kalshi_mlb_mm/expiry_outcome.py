"""Expired-quote outcome labeler — were we outpriced, or did nobody trade?

For every quote that expires unfilled, answers the margin-tuning question:
did the combo market print a trade while our quote was resting (a competing
maker won the RFQ — we were outpriced) or did no trade print at all (the
creator never executed with anyone — cutting margin would donate edge)?

Method (live-verified 2026-08-13): the public trade tape
`GET /markets/trades?ticker={combo_market_ticker}` returns trade prints with
`created_time` (UTC) and `taker_side`. Quirk: on MVE combo markets
yes_price/no_price/count come back None — timestamps + taker_side are the
usable fields, which is all the labeler needs.

Labels (per quote, decided once, after the grace window has fully elapsed):
  - competitor_traded_in_window  trade printed in [submitted_at, closed_at]
  - traded_shortly_after         first trade after closed_at landed within
                                 EXPIRY_OUTCOME_GRACE_SEC
  - no_trade                     neither

Inputs:  kalshi_mlb_mm.duckdb::live_quotes (terminal, unlabeled rows);
         Kalshi public trade tape (one GET per distinct ticker per sweep).
Outputs: INSERTs kalshi_mlb_mm.duckdb::quote_expiry_outcomes (one row per
         quote — the dedupe marker AND the queryable label store); one
         `quote_expiry_outcome` research event per quote (buffered, flushed
         by the main loop's research.flush()).
Safety:  zero writes to the quote path; never raises into the trading loop —
         per-ticker failures are logged and retried next sweep.
"""
import logging
from datetime import datetime, timedelta, timezone

from kalshi_common import auth_client
from kalshi_mlb_mm import config, db, research

log = logging.getLogger("kalshi_mlb_mm.expiry_outcome")

LABEL_COMPETITOR_TRADED = "competitor_traded_in_window"
LABEL_TRADED_SHORTLY_AFTER = "traded_shortly_after"
LABEL_NO_TRADE = "no_trade"

# One page covers the window: MVE combo markets print trades rarely (most
# expiries are exactly the no-trade case this module measures), so a single
# max-size page is effectively the full tape. If a tape ever fills the page
# AND its oldest print postdates a quote's window start, the label is
# computed from what we have and flagged tape_truncated (loud log + payload).
TRADES_FETCH_LIMIT = 1000

# The sweep shares the trading loop's single thread — bounded timeout so a
# hung tape call can't starve the 2s confirm loop (same rationale as
# settlement.API_TIMEOUT_SEC).
API_TIMEOUT_SEC = 10


def to_utc(dt: datetime) -> datetime:
    """UTC view of a DuckDB timestamp.

    TIMESTAMPTZ columns come back tz-aware in the session timezone
    (America/Los_Angeles on this machine); a naive value can only be a
    legacy pre-migration row, which was written in local time — both cases
    are exactly what astimezone() assumes.
    """
    return dt.astimezone(timezone.utc)


def parse_trade_time(created_time: str) -> datetime | None:
    """Aware UTC datetime from a tape print's created_time, else None."""
    try:
        return datetime.fromisoformat(str(created_time).replace("Z", "+00:00"))
    except (TypeError, ValueError):
        return None


def label_quote(submitted_at_utc: datetime, closed_at_utc: datetime,
                trade_times_utc: list[datetime],
                grace_sec: float) -> tuple[str, int, float | None]:
    """Pure labeling rule. Returns (label, n_trades_in_window,
    first_trade_after_close_sec).

    Window bounds are inclusive: a print stamped exactly at closed_at raced
    our expiry and still means a competitor executed while we were quotable.
    first_trade_after_close_sec is measured to the first print strictly
    after closed_at (None if the tape has none), regardless of label — it is
    the "how close behind the money were we" diagnostic.
    """
    n_in_window = sum(1 for t in trade_times_utc
                      if submitted_at_utc <= t <= closed_at_utc)
    after_close = [t for t in trade_times_utc if t > closed_at_utc]
    first_after_sec = (min(after_close) - closed_at_utc).total_seconds() \
        if after_close else None
    if n_in_window > 0:
        return LABEL_COMPETITOR_TRADED, n_in_window, first_after_sec
    if first_after_sec is not None and first_after_sec <= grace_sec:
        return LABEL_TRADED_SHORTLY_AFTER, 0, first_after_sec
    return LABEL_NO_TRADE, 0, first_after_sec


def _labelable_statuses() -> tuple[str, ...]:
    """'expired' always; 'cancelled' (quotes WE pulled) only behind the knob —
    pulls are our own risk decisions, so by default they'd muddy the headline
    "did the market want this price" ratio."""
    if config.EXPIRY_OUTCOME_INCLUDE_CANCELLED:
        return ("expired", "cancelled")
    return ("expired",)


def _pending_quotes() -> list[tuple]:
    """Unlabeled terminal quotes whose grace window has fully elapsed.

    Waiting until closed_at + grace is what makes the traded_shortly_after /
    no_trade split decidable in one pass — labeling earlier would need a
    re-visit, and each quote is labeled exactly once.
    """
    statuses = _labelable_statuses()
    grace_cutoff = (datetime.now(timezone.utc)
                    - timedelta(seconds=config.EXPIRY_OUTCOME_GRACE_SEC))
    placeholders = ",".join("?" for _ in statuses)
    with db.connect(read_only=True) as con:
        return con.execute(
            f"SELECT q.quote_id, q.combo_market_ticker, q.status, "
            f"q.submitted_at, q.closed_at "
            f"FROM live_quotes q "
            f"LEFT JOIN quote_expiry_outcomes o ON o.quote_id = q.quote_id "
            f"WHERE q.status IN ({placeholders}) "
            f"AND q.closed_at IS NOT NULL AND q.closed_at <= ? "
            f"AND o.quote_id IS NULL "
            f"ORDER BY q.closed_at",
            [*statuses, grace_cutoff]).fetchall()


def _fetch_trade_times(ticker: str) -> tuple[list[datetime], bool] | None:
    """Tape prints for one combo market → (UTC times, page_full flag),
    or None on any API failure (retried next sweep)."""
    status_code, body, _ = auth_client.api(
        "GET", f"/markets/trades?ticker={ticker}&limit={TRADES_FETCH_LIMIT}",
        timeout=API_TIMEOUT_SEC)
    if status_code != 200 or not isinstance(body, dict):
        log.warning("[expiry_tape_fetch_failed] ticker=%s status=%s — will "
                    "retry next sweep", ticker, status_code)
        return None
    trades = body.get("trades")
    if not isinstance(trades, list):
        log.warning("[expiry_tape_bad_shape] ticker=%s body_keys=%s — will "
                    "retry next sweep", ticker, list(body.keys()))
        return None
    times = []
    for trade in trades:
        parsed = parse_trade_time(trade.get("created_time")) \
            if isinstance(trade, dict) else None
        if parsed is None:
            log.warning("[expiry_tape_bad_trade] ticker=%s trade=%r — "
                        "skipping this print", ticker, trade)
            continue
        times.append(parsed.astimezone(timezone.utc))
    return times, len(trades) >= TRADES_FETCH_LIMIT


def _compute_outcome(quote_row: tuple, trade_times_utc: list[datetime],
                     page_full: bool) -> dict:
    """Pure label computation for one quote (no I/O)."""
    quote_id, ticker, quote_status, submitted_at, closed_at = quote_row
    submitted_utc = to_utc(submitted_at)
    closed_utc = to_utc(closed_at)
    label, n_in_window, first_after_sec = label_quote(
        submitted_utc, closed_utc, trade_times_utc,
        config.EXPIRY_OUTCOME_GRACE_SEC)
    tape_truncated = bool(
        page_full and trade_times_utc
        and min(trade_times_utc) > submitted_utc)
    if tape_truncated:
        log.warning("[expiry_tape_truncated] ticker=%s quote_id=%s — full "
                    "%d-print page starts after the quote window; label %r "
                    "computed from a partial tape", ticker, quote_id,
                    TRADES_FETCH_LIMIT, label)
    return dict(quote_id=quote_id, ticker=ticker, quote_status=quote_status,
                label=label, n_trades_in_window=n_in_window,
                first_trade_after_close_sec=first_after_sec,
                window_start=submitted_utc, window_end=closed_utc,
                n_tape_trades=len(trade_times_utc),
                tape_truncated=tape_truncated)


def _record_outcomes(outcomes: list[dict]) -> None:
    """Write one ticker's labels in ONE connection/transaction, then emit.

    One db.connect() per label was the [[duckdb-connect-cost]] anti-pattern
    (~17ms/open); batching per ticker keeps the sweep's DB cost flat.
    """
    now_utc = datetime.now(timezone.utc)
    with db.connect() as con:
        con.executemany(
            "INSERT OR REPLACE INTO quote_expiry_outcomes (quote_id, "
            "combo_market_ticker, quote_status, label, n_trades_in_window, "
            "first_trade_after_close_sec, window_start, window_end, "
            "recorded_at) VALUES (?,?,?,?,?,?,?,?,?)",
            [[o["quote_id"], o["ticker"], o["quote_status"], o["label"],
              o["n_trades_in_window"], o["first_trade_after_close_sec"],
              o["window_start"], o["window_end"], now_utc]
             for o in outcomes])
    for o in outcomes:
        research.emit("quote_expiry_outcome", ticker=o["ticker"],
                      quote_id=o["quote_id"],
                      payload=dict(label=o["label"],
                                   quote_status=o["quote_status"],
                                   n_trades_in_window=o["n_trades_in_window"],
                                   first_trade_after_close_sec=o[
                                       "first_trade_after_close_sec"],
                                   window_start=o["window_start"].isoformat(),
                                   window_end=o["window_end"].isoformat(),
                                   n_tape_trades=o["n_tape_trades"],
                                   tape_truncated=o["tape_truncated"]))
        log.info("[quote_expiry_outcome] quote_id=%s ticker=%s label=%s "
                 "n_in_window=%d first_after_close_sec=%s",
                 o["quote_id"], o["ticker"], o["label"],
                 o["n_trades_in_window"], o["first_trade_after_close_sec"])


def expiry_outcome_sweep_tick() -> None:
    """One sweep pass. Idempotent; API-failure tolerant; never raises.

    One tape fetch per distinct ticker; tickers beyond the per-sweep cap
    wait for the next sweep (quotes/day is tens, so the cap only matters
    on a backlog after downtime).
    """
    try:
        pending = _pending_quotes()
    except Exception as exc:
        log.warning("[expiry_sweep_error] listing unlabeled quotes: %s — "
                    "will retry next sweep", exc)
        return
    if not pending:
        return
    quotes_by_ticker: dict[str, list[tuple]] = {}
    for row in pending:
        quotes_by_ticker.setdefault(row[1], []).append(row)
    tickers = list(quotes_by_ticker)
    cap = config.EXPIRY_OUTCOME_MAX_TICKERS_PER_SWEEP
    if len(tickers) > cap:
        log.info("[expiry_sweep_capped] %d tickers pending, labeling first "
                 "%d this sweep", len(tickers), cap)
        tickers = tickers[:cap]
    for ticker in tickers:
        try:
            fetched = _fetch_trade_times(ticker)
            if fetched is None:
                continue
            trade_times_utc, page_full = fetched
            _record_outcomes([
                _compute_outcome(quote_row, trade_times_utc, page_full)
                for quote_row in quotes_by_ticker[ticker]])
        except Exception as exc:
            log.warning("[expiry_sweep_error] ticker=%s: %s — will retry "
                        "next sweep", ticker, exc)
