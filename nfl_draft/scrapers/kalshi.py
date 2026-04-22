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


def _extract_cents(market: dict, int_key: str, dollars_key: str):
    """Generic extractor: prefer legacy int cents, fall back to `_dollars` string.

    Kalshi's v2 response surfaces prices either as an integer cents field
    (`yes_bid`, `last_price`, ...) for accounts with `response_price_units=cents`,
    or as a dollar string (`yes_bid_dollars`, `last_price_dollars`, ...) for
    accounts on the newer dollar units. We accept both and normalize to int cents.

    Returns None when neither field carries a usable value (None, missing, or
    unparseable). Returns 0 when the field is explicitly 0 — callers decide
    whether 0 means "no market" or "zero price" per cascade rules.
    """
    raw_cents = market.get(int_key)
    if raw_cents not in (None, 0):
        try:
            return int(raw_cents)
        except (TypeError, ValueError):
            pass
    raw_dollars = market.get(dollars_key)
    if raw_dollars is not None:
        try:
            return int(round(float(raw_dollars) * 100))
        except (TypeError, ValueError):
            return None
    return raw_cents  # may be 0 or None


def _extract_yes_bid_cents(market: dict):
    """Return Kalshi yes_bid as an integer number of cents, or None.

    Accepts either legacy `yes_bid` (int cents) or new `yes_bid_dollars` (str)
    formats — see `_extract_cents` docstring.
    """
    return _extract_cents(market, "yes_bid", "yes_bid_dollars")


def _extract_yes_ask_cents(market: dict):
    """Return Kalshi yes_ask as an integer number of cents, or None."""
    return _extract_cents(market, "yes_ask", "yes_ask_dollars")


def _extract_no_bid_cents(market: dict):
    """Return Kalshi no_bid as an integer number of cents, or None."""
    return _extract_cents(market, "no_bid", "no_bid_dollars")


def _extract_no_ask_cents(market: dict):
    """Return Kalshi no_ask as an integer number of cents, or None."""
    return _extract_cents(market, "no_ask", "no_ask_dollars")


def _extract_last_price_cents(market: dict):
    """Return Kalshi last_price as an integer number of cents, or None."""
    return _extract_cents(market, "last_price", "last_price_dollars")


def _extract_candidate(market: dict) -> str:
    """Return the human-readable candidate name (stripped), or ``""``.

    Both the portal write path (``parse_markets_response`` -> OddsRow ->
    market_map) and the legacy write path (``_write_legacy_kalshi_odds`` ->
    kalshi_odds) must derive identical candidate strings so that the
    Cross-Book Grid tooltip join ``kalshi_odds.candidate = market_map.book_subject``
    succeeds. Prior to 2026-04-22 these two paths had slightly different
    fallback orders and whitespace handling; harmonizing here eliminates a
    latent join-miss risk even though no in-prod mismatches exist today.

    Fallback order, most specific first:
      1. ``yes_sub_title`` -- primary source (player name on the YES contract).
      2. ``subtitle`` -- older Kalshi response shape.
      3. ``custom_strike.Person`` -- structured candidate tags.
      4. ``custom_strike.Team`` -- team-based markets.
      5. ``no_sub_title`` -- last resort (NO-contract label).
    """
    candidate = (
        market.get("yes_sub_title")
        or market.get("subtitle")
        or (market.get("custom_strike") or {}).get("Person")
        or (market.get("custom_strike") or {}).get("Team")
        or market.get("no_sub_title")
        or ""
    )
    return (candidate or "").strip()


# Position series carry a tier suffix (P1 / P2 / P3 / ...) in the ticker's
# 2nd segment - one Kalshi market per (player, tier). Without scoping the
# book_label by that tier, all three "Mauigoa is the Nth OL" markets would
# collapse onto the same MARKET_MAP key (first_ol_francis-mauigoa) and the
# Cross-Book Grid's latest-per-(market,book) query would non-deterministically
# pick any of them. See _kalshi_book_label for the scoping rule.
_POSITION_SERIES = {
    "KXNFLDRAFTQB", "KXNFLDRAFTRB", "KXNFLDRAFTWR", "KXNFLDRAFTTE",
    "KXNFLDRAFTOL", "KXNFLDRAFTDB", "KXNFLDRAFTLB", "KXNFLDRAFTDT",
    "KXNFLDRAFTEDGE",
}


def _kalshi_book_label(series_ticker: str, ticker: str) -> str:
    """Return the MARKET_MAP key 'book_label' for a Kalshi market.

    For series where pick_number / range_high varies per market
    (KXNFLDRAFTPICK, KXNFLDRAFTTOP), the same player can appear at multiple
    pick numbers -- if we used bare series_ticker as book_label, all of those
    rows would collapse onto a single MARKET_MAP key and the resolver would
    incorrectly map, e.g., "Ty Simpson @ pick 3" to pick_10_ty-simpson. So we
    scope the label to the ticker's middle-segment prefix (e.g.
    'KXNFLDRAFTPICK-26-10'), which uniquely identifies the market group.

    Position series (KXNFLDRAFTQB / RB / WR / TE / OL / DB / LB / DT / EDGE)
    have the same problem at a different segment: Kalshi posts multiple tier
    markets per player (P1 = "1st at position", P2 = "2nd at position", ...).
    We scope the label to 'SERIES-P<N>' (e.g. 'KXNFLDRAFTOL-26P1') so the
    three tiers get distinct book_labels. Downstream, the MARKET_MAP builder
    only maps P1 to first_<pos>_<player>; P2+ fall through to quarantine
    (correct -- there's no canonical market_type for 2nd/3rd-at-position yet).

    For series where the market-type is fixed per series (e.g. KXNFLDRAFT1)
    we keep series_ticker as-is -- collisions can't happen.
    """
    if series_ticker in ("KXNFLDRAFTPICK", "KXNFLDRAFTTOP"):
        # Ticker shape: KXNFLDRAFTPICK-26-10-TSIM -> label 'KXNFLDRAFTPICK-26-10'
        parts = (ticker or "").split("-")
        if len(parts) >= 3:
            return "-".join(parts[:3])
    if series_ticker in _POSITION_SERIES:
        # Ticker shape: KXNFLDRAFTOL-26P1-FMAU -> label 'KXNFLDRAFTOL-26P1'
        parts = (ticker or "").split("-")
        if len(parts) >= 2:
            return f"{series_ticker}-{parts[1]}"
    if series_ticker == "KXNFLDRAFTTEAM":
        # Ticker shape: KXNFLDRAFTTEAM-26OCOO-WAS (series-playerTag-teamCode).
        # Scope label by player tag so 32 teams * ~16 players = ~512 distinct
        # (label, subject) pairs survive MARKET_MAP dedup.
        parts = (ticker or "").split("-")
        if len(parts) >= 2:
            return f"{series_ticker}-{parts[1]}"
    return series_ticker


def parse_markets_response(raw_response: dict, series_ticker: str) -> List[OddsRow]:
    """Convert a raw Kalshi /markets response into a list of OddsRow.

    Emits one row per Kalshi market with enough signal to price. For each
    market we compute two prices:

      * take price (implied_prob override): yes_ask first, then last_price.
        Represents the cost to buy a Yes contract. Used by the Cross-Book
        Grid's outlier flag.
      * fair value (devig_prob override): mid = (yes_bid + yes_ask)/200 when
        both sides posted, else last_price/100, else the one posted side/100.
        Used for cell display and the cross-venue median.

    Markets with no bid, no ask, and no last trade are skipped entirely --
    they carry no signal worth rendering. Rows with a fair value but no take
    price are still emitted; their ``implied_prob`` is None so downstream flag
    logic suppresses the outlier check for that cell.

    ``american_odds`` is set from the take price when available, otherwise from
    the fair value -- keeps legacy downstream readers (bet log, quarantine
    fallback for unmapped rows) aligned with what the scraper considers the
    most tradeable number for that row.
    """
    rows: List[OddsRow] = []
    if not isinstance(raw_response, dict):
        return rows
    now = datetime.now()

    for market in raw_response.get("markets", []) or []:
        yes_bid = _extract_yes_bid_cents(market) or 0
        yes_ask = _extract_yes_ask_cents(market) or 0
        last_price = _extract_last_price_cents(market) or 0

        # Fair-value cascade: mid, then last trade, then any single side.
        if yes_bid > 0 and yes_ask > 0:
            fair_cents = (yes_bid + yes_ask) / 2.0
        elif last_price > 0:
            fair_cents = float(last_price)
        elif yes_ask > 0:
            fair_cents = float(yes_ask)
        elif yes_bid > 0:
            fair_cents = float(yes_bid)
        else:
            continue  # no signal -> skip row

        if not (0 < fair_cents < 100):
            continue

        # Take-price cascade: ask, then last trade; else None (flag suppressed).
        if yes_ask > 0:
            take_cents = float(yes_ask)
        elif last_price > 0:
            take_cents = float(last_price)
        else:
            take_cents = None

        candidate = _extract_candidate(market)
        if not candidate:
            continue

        # American odds: use take price if available, else fair. Used by
        # downstream readers that still want an integer American-odds view.
        reference_cents = take_cents if take_cents is not None else fair_cents
        p_ref = reference_cents / 100.0
        if p_ref > 0.5:
            american = int(round(-100 * p_ref / (1 - p_ref)))
        else:
            american = int(round(100 * (1 - p_ref) / p_ref))

        rows.append(OddsRow(
            book="kalshi",
            book_label=_kalshi_book_label(series_ticker, market.get("ticker") or ""),
            book_subject=candidate,
            american_odds=american,
            fetched_at=now,
            implied_prob=(take_cents / 100.0) if take_cents is not None else None,
            devig_prob=fair_cents / 100.0,
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

    Legacy kalshi_odds write
    ------------------------
    This function also writes a fresh snapshot to the legacy `kalshi_odds`
    table (feeds the Price History / Edge Detection / Consensus dashboard
    tabs). Historically the adapter relied on a side effect of
    `legacy_fetcher.fetch_markets_for_series`, but that helper only FETCHES
    — the writes lived in `fetcher.run()`'s `save_odds_snapshot` +
    `cache_market_info`, which are only invoked under `__main__`. As a
    result `kalshi_odds` hadn't been updated since 2026-02-19 (the last
    time the legacy CLI was run). We now port the write logic directly
    into the adapter so it runs on every scrape tick.
    """
    rows: List[OddsRow] = []
    series_list = legacy_fetcher.discover_draft_series()
    for series in series_list:
        ticker = series["series_ticker"]
        # Collect every page of open markets for this series. Merge all
        # raw market dicts so the legacy kalshi_odds write below sees
        # everything, and pass each page through the OddsRow parser for
        # the portal.
        cursor = None
        all_raw_markets: List[dict] = []
        while True:
            path = f"/markets?series_ticker={ticker}&status=open&limit=100"
            if cursor:
                path += f"&cursor={cursor}"
            raw = public_request(path)
            if not raw:
                break
            rows.extend(parse_markets_response(raw, ticker))
            page_markets = raw.get("markets") or []
            all_raw_markets.extend(page_markets)
            cursor = raw.get("cursor")
            if not cursor or not page_markets:
                break

        # Legacy kalshi_odds write. See function docstring for why this
        # lives in the adapter rather than in legacy_fetcher.
        if all_raw_markets:
            _write_legacy_kalshi_odds(ticker, series.get("title", ""), all_raw_markets)
    return rows


def _write_legacy_kalshi_odds(series_ticker: str, series_title: str, markets: list) -> None:
    """Persist a fresh kalshi_odds + market_info + draft_series snapshot.

    Restores the write path that `kalshi_draft/fetcher.py::run` used to do
    at the end of every fetch (`save_odds_snapshot` + `cache_market_info`).
    Columns stay identical to the legacy schema (see DESCRIBE kalshi_odds
    in the production DB) so existing dashboard queries
    (get_latest_odds / get_price_history / get_snapshot_count) keep working
    unchanged.

    market_info + draft_series have no PK in the production DB, so plain
    INSERT is used instead of INSERT OR REPLACE — DuckDB's binder rejects
    ON CONFLICT without a constraint target. The tables grow ~1 row per
    fetch-tick per ticker, which is cheap (few hundred rows/day).
    """
    fetch_time = datetime.now()
    rows_to_write = []
    market_info_rows = []
    for m in markets:
        candidate = _extract_candidate(m)
        if not candidate:
            continue
        rows_to_write.append((
            fetch_time,
            series_ticker,
            m.get("event_ticker", ""),
            m.get("ticker", ""),
            m.get("title", ""),
            candidate,
            _extract_yes_bid_cents(m) or 0,
            _extract_yes_ask_cents(m) or 0,
            _extract_no_bid_cents(m) or 0,
            _extract_no_ask_cents(m) or 0,
            _extract_last_price_cents(m) or 0,
            m.get("volume", 0) or 0,
            m.get("volume_24h", 0) or 0,
            int(float(m.get("liquidity_dollars") or 0)),
            m.get("open_interest", 0) or 0,
        ))
        market_info_rows.append((
            m.get("ticker", ""),
            m.get("title", ""),
            candidate,
            series_ticker,
            fetch_time,
        ))

    if not rows_to_write:
        return

    with write_connection() as con:
        con.executemany(
            "INSERT INTO kalshi_odds VALUES (?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?)",
            rows_to_write,
        )
        con.executemany(
            "INSERT INTO market_info (ticker, title, subtitle, series_ticker, updated_at) "
            "VALUES (?, ?, ?, ?, ?)",
            market_info_rows,
        )
        # Also touch draft_series so the legacy dashboard sees a fresh title.
        # Classification reuses the legacy categorizer for parity.
        try:
            category = legacy_fetcher.classify_series(series_ticker, series_title)
            con.execute(
                "INSERT INTO draft_series VALUES (?, ?, ?, ?)",
                [series_ticker, series_title, category, fetch_time],
            )
        except Exception:
            # draft_series is non-critical for the portal; swallow quietly.
            pass


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
