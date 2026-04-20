"""Kalshi adapter: calls legacy kalshi_draft/fetcher.py for raw market data,
returns List[OddsRow] for the new nfl_draft pipeline.

The legacy fetcher's side-effect writes (kalshi_odds, draft_series,
market_info, positions, resting_orders) happen as a byproduct of the calls
we make here - this keeps the Kalshi dashboard + MM bot fed from the same
code path that the new portal uses.

Trade tape (fetch_trades):
------------------------
Kalshi's /markets/trades endpoint IGNORES the `series_ticker=` query param
(verified live 2026-04-17) - it only filters on `ticker=<specific-market>`.
So the trade poller:
  1. Discovers all NFL draft series (reusing discover_draft_series).
  2. For each series, pulls the distinct open-market tickers from
     kalshi_odds (cheap - that table is populated by fetch_draft_odds on
     every scrape tick).
  3. For each ticker, calls /markets/trades?ticker=...&min_ts=<epoch>&limit=100,
     paginating by cursor.
  4. INSERT OR IGNOREs by trade_id PK (natural dedup).
  5. Advances kalshi_poll_state.last_traded_at per series to max(traded_at).

Note: Kalshi's `min_ts` param is EPOCH SECONDS, not ISO. We convert the
local `last_traded_at` cursor to epoch-seconds at the boundary.
"""

import sys
import traceback
from datetime import datetime, timezone
from pathlib import Path
from typing import List

# Ensure both the repo root AND kalshi_draft/ itself are importable - the
# legacy fetcher uses bare `from auth import ...` and `from db import ...`.
_REPO_ROOT = Path(__file__).resolve().parent.parent.parent
_KALSHI_DIR = _REPO_ROOT / "kalshi_draft"
for p in (_REPO_ROOT, _KALSHI_DIR):
    sp = str(p)
    if sp not in sys.path:
        sys.path.insert(0, sp)

from kalshi_draft import fetcher as legacy_fetcher  # noqa: E402
from kalshi_draft.auth import public_request  # noqa: E402

from nfl_draft.lib.db import (  # noqa: E402
    write_connection,
    read_connection,
    utc_iso_to_local,
)
from nfl_draft.scrapers._base import OddsRow, TradeRow


def _extract_yes_bid_cents(market: dict):
    """Return Kalshi yes_bid as an integer number of cents, or None.

    Kalshi's v2 response includes both legacy (int cents) and new
    (`_dollars` string) keys depending on the account's `response_price_units`.
    We accept either and normalize to cents.
    """
    raw_cents = market.get("yes_bid")
    if raw_cents not in (None, 0):
        try:
            return int(raw_cents)
        except (TypeError, ValueError):
            pass
    raw_dollars = market.get("yes_bid_dollars")
    if raw_dollars is not None:
        try:
            return int(round(float(raw_dollars) * 100))
        except (TypeError, ValueError):
            return None
    return raw_cents  # may be 0 or None


def _kalshi_book_label(series_ticker: str, ticker: str) -> str:
    """Return the MARKET_MAP key 'book_label' for a Kalshi market.

    For series where pick_number / range_high varies per market
    (KXNFLDRAFTPICK, KXNFLDRAFTTOP), the same player can appear at multiple
    pick numbers -- if we used bare series_ticker as book_label, all of those
    rows would collapse onto a single MARKET_MAP key and the resolver would
    incorrectly map, e.g., "Ty Simpson @ pick 3" to pick_10_ty-simpson. So we
    scope the label to the ticker's middle-segment prefix (e.g.
    'KXNFLDRAFTPICK-26-10'), which uniquely identifies the market group.

    For series where the market-type is fixed per series (e.g. KXNFLDRAFT1,
    KXNFLDRAFTQB) we keep series_ticker as-is -- collisions can't happen.
    """
    if series_ticker in ("KXNFLDRAFTPICK", "KXNFLDRAFTTOP"):
        # Ticker shape: KXNFLDRAFTPICK-26-10-TSIM -> label 'KXNFLDRAFTPICK-26-10'
        parts = (ticker or "").split("-")
        if len(parts) >= 3:
            return "-".join(parts[:3])
    return series_ticker


def parse_markets_response(raw_response: dict, series_ticker: str) -> List[OddsRow]:
    """Convert a raw Kalshi /markets response into a list of OddsRow.

    One row per market with a usable yes_bid. Skips markets with zero/empty
    bids (no live market price) and markets missing a candidate name.
    """
    rows: List[OddsRow] = []
    if not isinstance(raw_response, dict):
        return rows
    now = datetime.now()

    for market in raw_response.get("markets", []) or []:
        yes_bid_cents = _extract_yes_bid_cents(market)
        if yes_bid_cents is None or yes_bid_cents <= 0:
            continue
        p = yes_bid_cents / 100.0
        if p <= 0 or p >= 1:
            continue
        # Convert implied probability to American odds
        if p > 0.5:
            american = int(round(-100 * p / (1 - p)))
        else:
            american = int(round(100 * (1 - p) / p))

        candidate = (
            market.get("yes_sub_title")
            or market.get("subtitle")
            or (market.get("custom_strike") or {}).get("Person")
            or (market.get("custom_strike") or {}).get("Team")
            or ""
        )
        candidate = (candidate or "").strip()
        if not candidate:
            continue

        rows.append(OddsRow(
            book="kalshi",
            book_label=_kalshi_book_label(series_ticker, market.get("ticker") or ""),
            book_subject=candidate,
            american_odds=american,
            fetched_at=now,
        ))
    return rows


def fetch_draft_odds() -> List[OddsRow]:
    """Discover all NFL Draft series, fetch their open markets, return rows.

    Pagination
    ----------
    Kalshi's /markets endpoint caps each page at 100 markets. KXNFLDRAFTPICK
    alone carries ~649 open markets (32 teams * ~20 candidates each across
    multiple pick tiers), so we MUST walk the cursor — otherwise ~85% of
    markets silently drop and downstream portal rows go stale.

    Side effect: calling the legacy fetcher writes kalshi_odds + draft_series
    + market_info into nfl_draft.duckdb (the legacy fetcher has been
    repointed to the new DB in Task 3-9).
    """
    rows: List[OddsRow] = []
    series_list = legacy_fetcher.discover_draft_series()
    for series in series_list:
        ticker = series["series_ticker"]
        # Legacy side-effect write (fills kalshi_odds via the run pipeline).
        # discover + fetch_markets_for_series together are what the legacy
        # dashboard depends on - we mirror that call pattern so nothing
        # downstream regresses.
        legacy_fetcher.fetch_markets_for_series(ticker)
        # Walk all pages for this series. Previously we hit only the first
        # 100 markets per series; KXNFLDRAFTPICK alone has ~649 so ~85% were
        # dropped, causing draft_odds for "missing" players to go stale from
        # an earlier scrape that happened to include them. Cursor pattern
        # mirrors _enumerate_fallback lower in this file.
        cursor = None
        while True:
            path = f"/markets?series_ticker={ticker}&status=open&limit=100"
            if cursor:
                path += f"&cursor={cursor}"
            raw = public_request(path)
            if not raw:
                break
            rows.extend(parse_markets_response(raw, ticker))
            cursor = raw.get("cursor")
            if not cursor or not raw.get("markets"):
                break
    return rows


def _local_to_epoch_seconds(local_ts: datetime) -> int:
    """Convert naive local datetime -> UTC epoch seconds (for Kalshi min_ts)."""
    aware_local = local_ts.astimezone()  # interpret naive as system local
    return int(aware_local.astimezone(timezone.utc).timestamp())


def parse_trades_response(raw: dict) -> List[TradeRow]:
    """Convert a raw Kalshi /markets/trades response into TradeRow list.

    Kalshi's v2 response shape (as of 2026-04-17) carries:
      trade_id, ticker, taker_side, yes_price_dollars (str), count_fp (str),
      created_time (ISO-8601 UTC with Z or offset).

    We also tolerate the older `yes_price` (int cents) / `count` (int) /
    `side` fields in case an account is toggled to `response_price_units=cents`
    or the API evolves. Rows missing critical fields are skipped.
    """
    rows: List[TradeRow] = []
    if not isinstance(raw, dict):
        return rows
    now = datetime.now()

    for trade in raw.get("trades", []) or []:
        trade_id = trade.get("trade_id")
        ticker = trade.get("ticker")
        if not trade_id or not ticker:
            continue

        # Side: prefer taker_side (Kalshi v2); fall back to generic "side".
        side_raw = trade.get("taker_side") or trade.get("side") or "yes"
        side = str(side_raw).lower()

        # Price: normalize to integer cents, accepting either shape.
        if "yes_price" in trade and trade["yes_price"] is not None:
            try:
                price_cents = int(trade["yes_price"])
            except (TypeError, ValueError):
                continue
        elif "yes_price_dollars" in trade and trade["yes_price_dollars"] is not None:
            try:
                price_cents = int(round(float(trade["yes_price_dollars"]) * 100))
            except (TypeError, ValueError):
                continue
        else:
            continue

        # Count: accept `count` (int) or `count_fp` (string/float, new API).
        # Kalshi uses `count_fp` for fractional-contract markets; for binary
        # YES/NO it's always an integer-valued float ("33.00"). We round
        # defensively so partial contracts don't silently truncate.
        if "count" in trade and trade["count"] is not None:
            try:
                count = int(trade["count"])
            except (TypeError, ValueError):
                continue
        elif "count_fp" in trade and trade["count_fp"] is not None:
            try:
                count = int(round(float(trade["count_fp"])))
            except (TypeError, ValueError):
                continue
        else:
            continue

        created = trade.get("created_time")
        if not created:
            continue
        try:
            traded_at = utc_iso_to_local(created)
        except (TypeError, ValueError):
            continue

        rows.append(TradeRow(
            trade_id=str(trade_id),
            ticker=str(ticker),
            side=side,
            price_cents=price_cents,
            count=count,
            traded_at=traded_at,
            fetched_at=now,
        ))
    return rows


def _tickers_for_series(con, series_ticker: str) -> List[str]:
    """Return distinct market tickers we've seen for this series.

    We use kalshi_odds as the ticker source: fetch_draft_odds writes every
    open market it sees to that table on every scrape tick, so it's a fresh
    enumeration of tickers worth polling for trades. Falls back to a live
    /markets query if the table is empty (e.g. first-run, fresh DB).
    """
    try:
        rows = con.execute(
            "SELECT DISTINCT ticker FROM kalshi_odds WHERE series_ticker = ?",
            [series_ticker],
        ).fetchall()
        tickers = [r[0] for r in rows if r[0]]
    except Exception:
        tickers = []

    if tickers:
        return tickers

    # Fallback: enumerate open markets live. Paginated via cursor.
    tickers = []
    cursor = None
    while True:
        path = f"/markets?series_ticker={series_ticker}&status=open&limit=100"
        if cursor:
            path += f"&cursor={cursor}"
        data = public_request(path)
        if not data:
            break
        for m in data.get("markets", []) or []:
            t = m.get("ticker")
            if t:
                tickers.append(t)
        cursor = data.get("cursor")
        if not cursor or not data.get("markets"):
            break
    return tickers


def fetch_trades() -> List[TradeRow]:
    """Poll Kalshi /markets/trades for each NFL draft series since last cursor.

    Returns the list of NEWLY-ingested TradeRows (already written to
    kalshi_trades via INSERT OR IGNORE). Callers can len() this to report
    ingestion count; empty list is a normal outcome during quiet periods.
    """
    ingested: List[TradeRow] = []

    try:
        series_list = legacy_fetcher.discover_draft_series()
    except Exception as e:
        print(f"  [trades] discover_draft_series failed: {e}")
        return ingested

    # Pull cursor state once and hold ticker enumeration inside a read conn.
    with read_connection() as rcon:
        cursor_by_series = {}
        try:
            cursor_rows = rcon.execute(
                "SELECT series_ticker, last_traded_at FROM kalshi_poll_state"
            ).fetchall()
            cursor_by_series = {r[0]: r[1] for r in cursor_rows}
        except Exception:
            # Table missing (shouldn't happen post-init_schema); treat as empty.
            cursor_by_series = {}

        tickers_by_series = {}
        for s in series_list:
            st = s["series_ticker"]
            tickers_by_series[st] = _tickers_for_series(rcon, st)

    now_local = datetime.now()

    for s in series_list:
        series_ticker = s["series_ticker"]
        last_local = cursor_by_series.get(series_ticker)
        # Kalshi min_ts is EPOCH SECONDS (not ISO) - verified live.
        min_ts = _local_to_epoch_seconds(last_local) if last_local else None

        max_traded_at = last_local  # track advance

        for ticker in tickers_by_series.get(series_ticker, []):
            # Per-ticker try/except so one bad ticker can't break the poll.
            try:
                cursor = None
                while True:
                    path = f"/markets/trades?ticker={ticker}&limit=100"
                    if min_ts is not None:
                        path += f"&min_ts={min_ts}"
                    if cursor:
                        path += f"&cursor={cursor}"

                    raw = public_request(path)
                    if not raw:
                        break

                    batch = parse_trades_response(raw)
                    if batch:
                        with write_connection() as wcon:
                            for t in batch:
                                wcon.execute(
                                    "INSERT OR IGNORE INTO kalshi_trades "
                                    "VALUES (?, ?, ?, ?, ?, ?, ?, ?)",
                                    [
                                        t.trade_id, t.ticker, t.side,
                                        t.price_cents, t.count,
                                        t.count * t.price_cents * 0.01,
                                        t.traded_at, t.fetched_at,
                                    ],
                                )
                        ingested.extend(batch)
                        # Track max traded_at for cursor advance.
                        batch_max = max(t.traded_at for t in batch)
                        if max_traded_at is None or batch_max > max_traded_at:
                            max_traded_at = batch_max

                    cursor = raw.get("cursor")
                    if not cursor or not raw.get("trades"):
                        break
            except Exception as e:
                print(f"  [trades] ticker={ticker} failed: {e}")
                traceback.print_exc()
                continue

        # Update poll_state for this series, only if we actually saw trades
        # (idempotent: don't rewrite the cursor when the window was empty).
        try:
            with write_connection() as wcon:
                if max_traded_at is not None and max_traded_at != last_local:
                    wcon.execute(
                        "INSERT INTO kalshi_poll_state "
                        "(series_ticker, last_traded_at, last_polled_at) "
                        "VALUES (?, ?, ?) "
                        "ON CONFLICT (series_ticker) DO UPDATE SET "
                        "last_traded_at=excluded.last_traded_at, "
                        "last_polled_at=excluded.last_polled_at",
                        [series_ticker, max_traded_at, now_local],
                    )
                else:
                    # Still record that we polled (for observability), but
                    # don't regress the traded_at cursor.
                    wcon.execute(
                        "INSERT INTO kalshi_poll_state "
                        "(series_ticker, last_traded_at, last_polled_at) "
                        "VALUES (?, ?, ?) "
                        "ON CONFLICT (series_ticker) DO UPDATE SET "
                        "last_polled_at=excluded.last_polled_at",
                        [series_ticker, last_local, now_local],
                    )
        except Exception as e:
            print(f"  [trades] poll_state update failed for {series_ticker}: {e}")

    return ingested
