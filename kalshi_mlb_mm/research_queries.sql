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
WITH joined AS (
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
FROM joined;
