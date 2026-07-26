-- MLB capture readout queries for unabated_edge_market.duckdb
--
-- Run against a READ-ONLY COPY of the live DB, never the live file directly
-- while the bot is running (WAL-hang caveat from the WC maker). Copy target:
-- a scratch dir OUTSIDE the repo tree — e.g. $CLAUDE_JOB_DIR/tmp or a
-- user-chosen path. NEVER copy a .duckdb into the repo (see CLAUDE.md rule 5
-- on symlinking/copying DuckDB files) and never `cp` into /tmp (shared,
-- unmanaged). Example:
--   cp unabated_edge/unabated_edge_market.duckdb "$CLAUDE_JOB_DIR/tmp/ro.duckdb"
--   duckdb -readonly "$CLAUDE_JOB_DIR/tmp/ro.duckdb" < unabated_edge/mlb_readout.sql
--
-- Schema notes (unabated_edge/storage.py CREATE TABLE, checked 2026-07-25):
--   book_snapshots: ts, sport, market_ticker, floor_strike, yes_bid,
--     yes_bid_qty, no_bid, no_bid_qty, yes_ask, no_ask, volume,
--     open_interest, depth (JSON). One row per (market_ticker, tick) — this
--     is already top-of-book, not a full ladder, so yes_ask/yes_bid are
--     used directly (no "best_" prefix).
--   IMPORTANT: yes_ask/no_ask have NO matching qty column because Kalshi's
--   YES/NO books are complementary: yes_ask = 1 - no_bid and
--   no_ask = 1 - yes_bid (see venues/kalshi.py yes_ask_from_book /
--   no_ask_from_book). So the depth backing the yes-ask touch is
--   no_bid_qty, and the depth backing the no-ask touch is yes_bid_qty.
--   kalshi_trades: trade_id (PK), sport, market_ticker, created_time
--     (raw ISO string), yes_price, count, taker_side, captured_ts
--     (TIMESTAMPTZ, when we recorded the trade — used for the ASOF join
--     since it's the same tz-aware clock as book_snapshots.ts).
--   line_snapshots: ts, sport, event_id, market_source_id, bet_type, side,
--     price, points, modified_on. No market_ticker column, so it cannot be
--     joined to book_snapshots/kalshi_trades directly (no shared key) —
--     not used below.

-- 1. Crowd spread distribution per rung (is this a WC-like 1c book?)
SELECT market_ticker,
       round(avg(yes_ask - yes_bid), 4) AS avg_spread,
       round(quantile_cont(yes_ask - yes_bid, 0.5), 4) AS med_spread,
       count(*) AS snaps
FROM book_snapshots
WHERE sport = 'mlb' AND yes_ask IS NOT NULL AND yes_bid IS NOT NULL
GROUP BY market_ticker
ORDER BY med_spread DESC;

-- 2. Depth at touch (queue we'd be joining). No yes_ask_qty/no_ask_qty
-- column exists because yes_ask = 1 - no_bid and no_ask = 1 - yes_bid, so
-- the ask-side depth is read off the complementary bid's qty column.
SELECT market_ticker,
       round(avg(yes_bid_qty), 0) AS avg_yes_bid_depth,
       round(avg(no_bid_qty), 0) AS avg_yes_ask_depth
FROM book_snapshots
WHERE sport = 'mlb'
GROUP BY market_ticker
ORDER BY avg_yes_bid_depth DESC;

-- 3. Share of trades within 1c of mid (the WC kill-shot number was 99.99%)
WITH t AS (
  SELECT tr.market_ticker, tr.yes_price,
         (bs.yes_bid + bs.yes_ask) / 2 AS mid
  FROM kalshi_trades tr
  ASOF JOIN book_snapshots bs
    ON tr.market_ticker = bs.market_ticker AND tr.captured_ts >= bs.ts
  WHERE tr.sport = 'mlb' AND bs.sport = 'mlb'
)
SELECT count(*) AS trades,
       round(100.0 * avg(CASE WHEN abs(yes_price - mid) <= 0.01 THEN 1 ELSE 0 END), 2)
         AS pct_within_1c
FROM t
WHERE mid IS NOT NULL;

-- 4. Anchor coverage: distinct rungs Kalshi lists that we've snapshotted
-- (line_snapshots has no market_ticker to join against, so this counts
-- snapshot coverage, not a fair-vs-listed comparison — see schema note above)
SELECT sport, count(DISTINCT market_ticker) AS rungs_snapshotted
FROM book_snapshots WHERE sport = 'mlb' GROUP BY sport;

-- 5. Retention prune (disk guard; run when the market DB grows past ~2 GB)
-- DELETE FROM book_snapshots WHERE sport='mlb' AND ts          < now() - INTERVAL 14 DAY;
-- DELETE FROM kalshi_trades  WHERE sport='mlb' AND captured_ts < now() - INTERVAL 14 DAY;
