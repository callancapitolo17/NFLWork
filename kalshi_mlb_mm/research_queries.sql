-- Research firehose analysis queries for the Kalshi MLB MM (maker) bot.
--
-- The research DB (kalshi_mlb_mm_research.duckdb) is a separate sibling of the
-- trading DB. Run from the kalshi_mlb_mm/ directory, e.g.:
--   duckdb < research_queries.sql
-- or open an interactive shell and ATTACH both DBs:
--   duckdb -cmd "ATTACH 'kalshi_mlb_mm.duckdb' AS state (READ_ONLY); ATTACH 'kalshi_mlb_mm_research.duckdb' AS research (READ_ONLY);"
--
-- Payload fields are JSON. Use json_extract_string(payload, 'key') to read
-- them, NOT the ->> operator: in DuckDB 1.4.x, ->> can try to cast the whole
-- payload to a number when the same query both extracts a field in SELECT and
-- compares one in WHERE over rows whose JSON holds numeric values.
-- json_extract_string always returns VARCHAR and sidesteps that.
-- NOTE: it returns VARCHAR even for numeric payload fields (blended_fair,
-- current_rate, ...). To FILTER or compare them numerically, cast explicitly:
--   WHERE CAST(json_extract_string(payload,'agreeing_count') AS INTEGER) > 2

ATTACH 'kalshi_mlb_mm.duckdb'          AS state    (READ_ONLY);
ATTACH 'kalshi_mlb_mm_research.duckdb' AS research (READ_ONLY);

-- ---------------------------------------------------------------------------
-- 1) FIRST-FILL VERIFICATION — what fields did Kalshi actually return?
--    Read the payload to confirm field names (accepted_side, contracts, etc.)
--    match our assumptions before we rely on them in production code.
-- ---------------------------------------------------------------------------
SELECT json_extract_string(payload, 'response_keys') AS keys,
       payload
FROM research.events
WHERE event_type = 'accept_observed'
ORDER BY ts DESC LIMIT 1;

-- ---------------------------------------------------------------------------
-- 2) PER-CREATOR FILL COUNTS AND SKIP-RATE (last 24h).
--    Which counterparties are farming us? Which are well-behaved?
-- ---------------------------------------------------------------------------
WITH per_creator AS (
    SELECT json_extract_string(r.payload, 'rfq_raw') AS rfq_raw,
           f.rfq_id
    FROM research.events r
    LEFT JOIN research.events f ON f.rfq_id = r.rfq_id
                                AND f.event_type = 'fill_recorded'
    WHERE r.event_type = 'rfq_received'
      AND r.ts > now() - INTERVAL 24 HOUR
)
SELECT json_extract_string(rfq_raw, 'creator_user_id') AS creator_id,
       COUNT(*) AS rfqs_seen,
       COUNT(rfq_id) AS rfqs_filled
FROM per_creator
GROUP BY 1 ORDER BY rfqs_filled DESC;

-- ---------------------------------------------------------------------------
-- 3) VOID-RATE BREAKDOWN BY REASON (last 1h).
--    Which void paths are dominating? High voided_no_fresh_books = books stale.
--    High voided_last_look = our fair is drifting between quote and accept.
-- ---------------------------------------------------------------------------
SELECT json_extract_string(payload, 'reason') AS reason,
       COUNT(*) AS n
FROM research.events
WHERE event_type = 'decision'
  AND json_extract_string(payload, 'decision') LIKE 'voided_%'
  AND ts > now() - INTERVAL 1 HOUR
GROUP BY 1 ORDER BY 2 DESC;

-- ---------------------------------------------------------------------------
-- 4) RECONCILE MISMATCH RATE + OUTCOMES (last 24h).
--    Phantom fills, max-age fallbacks, or mismatches all point to Kalshi
--    API or timing issues that deserve investigation.
-- ---------------------------------------------------------------------------
SELECT json_extract_string(payload, 'outcome') AS outcome,
       COUNT(*) AS n
FROM research.events
WHERE event_type = 'reconcile_done'
  AND ts > now() - INTERVAL 24 HOUR
GROUP BY 1 ORDER BY 2 DESC;

-- ---------------------------------------------------------------------------
-- 5) PER-BOOK CONSENSUS PARTICIPATION RATE (last 24h).
--    Which books are frequently missing from the agreeing set?
--    A book that drops out often may be systematically stale.
-- ---------------------------------------------------------------------------
SELECT key AS book, COUNT(*) AS n_with_quote
FROM research.events,
     LATERAL (
         SELECT UNNEST(json_keys(json_extract(payload, '$.book_fairs'))) AS key
     )
WHERE event_type = 'quote_priced'
  AND ts > now() - INTERVAL 24 HOUR
GROUP BY 1 ORDER BY 2 DESC;

-- ---------------------------------------------------------------------------
-- 6) QUOTE → ACCEPT → FILL LATENCY P50/P95.
--    How long between when we price a quote and when Kalshi accepts it?
--    How long from accept to fill confirmation? Helps tune loop intervals.
-- ---------------------------------------------------------------------------
WITH latency_joined AS (
    SELECT qp.ts AS quoted_at, ao.ts AS accepted_at, fr.ts AS filled_at
    FROM research.events qp
    JOIN research.events ao ON ao.event_type = 'accept_observed'
                            AND ao.quote_id = qp.quote_id
    JOIN research.events fr ON fr.event_type = 'fill_recorded'
                            AND fr.quote_id = qp.quote_id
    WHERE qp.event_type = 'quote_priced'
)
SELECT
    percentile_cont(0.50) WITHIN GROUP (ORDER BY EXTRACT(EPOCH FROM (accepted_at - quoted_at)) * 1000) AS quote_to_accept_p50_ms,
    percentile_cont(0.95) WITHIN GROUP (ORDER BY EXTRACT(EPOCH FROM (accepted_at - quoted_at)) * 1000) AS quote_to_accept_p95_ms,
    percentile_cont(0.50) WITHIN GROUP (ORDER BY EXTRACT(EPOCH FROM (filled_at - accepted_at)) * 1000) AS accept_to_fill_p50_ms
FROM latency_joined;

-- ---------------------------------------------------------------------------
-- 7) P&L ATTRIBUTION PER SETTLED FILL (issue #12).
--    Decomposes each settled fill's PER-CONTRACT P&L into:
--      quoted_margin_per_ct      = fair(held side, at quote) - price paid
--                                  → the edge we THOUGHT we were capturing
--      fair_drift_per_ct         = fair moved between quote and confirm
--                                  → adverse selection our last look tolerated
--      settlement_vs_fair_per_ct = binary outcome (0/1) vs fair at confirm
--                                  → noise per fill that converges to
--                                    residual fair error over many fills
--      fee                       = maker fee per contract
--    Identity: (margin + drift + outcome_vs_fair - fee) * contracts
--              == realized_pnl, so identity_check should be ~0 for every row.
--    Uses state.fills + state.settlements (written by the settlement sweep).
-- ---------------------------------------------------------------------------
SELECT f.fill_id,
       f.combo_market_ticker,
       f.filled_at,
       f.side_held,
       f.contracts,
       f.price,
       f.fee,
       s.result,
       CASE WHEN f.side_held = 'yes' THEN f.blended_fair_at_quote
            ELSE 1 - f.blended_fair_at_quote END               AS fair_side_at_quote,
       CASE WHEN f.side_held = 'yes' THEN f.fair_at_confirm
            ELSE 1 - f.fair_at_confirm END                     AS fair_side_at_confirm,
       CASE WHEN s.result = f.side_held THEN 1.0 ELSE 0.0 END  AS outcome,
       fair_side_at_quote - f.price                            AS quoted_margin_per_ct,
       fair_side_at_confirm - fair_side_at_quote               AS fair_drift_per_ct,
       outcome - fair_side_at_confirm                          AS settlement_vs_fair_per_ct,
       f.realized_pnl,
       (quoted_margin_per_ct + fair_drift_per_ct
        + settlement_vs_fair_per_ct - f.fee) * f.contracts
         - f.realized_pnl                                      AS identity_check
FROM state.fills f
JOIN state.settlements s ON s.combo_market_ticker = f.combo_market_ticker
WHERE f.realized_pnl IS NOT NULL
ORDER BY f.filled_at;

-- ---------------------------------------------------------------------------
-- 8) THE MEASUREMENT-PHASE HEADLINE: does the quoted margin survive?
--    Aggregate of query 7 — average quoted margin vs average realized P&L
--    per contract. realized < quoted by more than noise = adverse selection
--    is eating the margin.
-- ---------------------------------------------------------------------------
SELECT COUNT(*)                                                   AS settled_fills,
       SUM(f.contracts)                                           AS contracts,
       AVG(CASE WHEN f.side_held = 'yes' THEN f.blended_fair_at_quote
                ELSE 1 - f.blended_fair_at_quote END - f.price)   AS avg_quoted_margin_per_ct,
       SUM(f.realized_pnl) / NULLIF(SUM(f.contracts), 0)          AS realized_pnl_per_ct,
       SUM(f.realized_pnl)                                        AS total_realized_pnl
FROM state.fills f
WHERE f.realized_pnl IS NOT NULL;

-- ---------------------------------------------------------------------------
-- 9) LIVE-PRICING HEALTH (issue #54): daily veto/void rate. In live mode the
--    confirm last look re-checks a SECONDS-old number, so the void rate
--    should collapse vs the cache era. quote_priced carries pricing_mode, so
--    the daily mode mix is joined in — if voids do not collapse once the mix
--    goes live, live pricing is not actually working.
-- ---------------------------------------------------------------------------
WITH daily_decisions AS (
    SELECT CAST(observed_at AS DATE)                                AS day,
           SUM(CASE WHEN decision LIKE 'voided_%' THEN 1 ELSE 0 END) AS voided,
           SUM(CASE WHEN decision = 'confirmed'   THEN 1 ELSE 0 END) AS confirmed
    FROM state.quote_decisions
    GROUP BY 1
), daily_modes AS (
    SELECT CAST(ts AS DATE)                                          AS day,
           json_extract_string(payload, 'pricing_mode')              AS pricing_mode,
           COUNT(*)                                                  AS quotes_priced
    FROM research.events
    WHERE event_type = 'quote_priced'
    GROUP BY 1, 2
)
SELECT d.day,
       m.pricing_mode,
       m.quotes_priced,
       d.voided,
       d.confirmed,
       d.voided * 1.0 / NULLIF(d.voided + d.confirmed, 0) AS void_rate
FROM daily_decisions d
LEFT JOIN daily_modes m USING (day)
ORDER BY d.day, m.pricing_mode;

-- ---------------------------------------------------------------------------
-- 10) SHADOW ROLLOUT EVIDENCE (issue #54): cached-vs-live fair pairs. The
--     decision dataset for flipping QUOTE_PRICING_MODE to 'live': a small,
--     centered cached_minus_live distribution means the cache was not lying;
--     a fat or skewed one is the adverse-selection window live pricing
--     closes. Novig caveat: NV rows inherit its open sweep-price broadcast
--     defect — check per-book live fairs in the live_games payload before
--     leaning on NV-heavy pairs.
-- ---------------------------------------------------------------------------
SELECT ts,
       ticker,
       CAST(json_extract_string(payload, 'cached_fair') AS DOUBLE)  AS cached_fair,
       CAST(json_extract_string(payload, 'live_fair')   AS DOUBLE)  AS live_fair,
       cached_fair - live_fair                                      AS cached_minus_live,
       json_extract_string(payload, 'cached_reason')                AS cached_reason,
       json_extract_string(payload, 'live_reason')                  AS live_reason,
       json_extract_string(payload, 'live_games')                   AS live_games
FROM research.events
WHERE event_type = 'shadow_fair_comparison'
ORDER BY ts DESC;
