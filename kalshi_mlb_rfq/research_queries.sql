-- Research firehose analysis queries (Approach 2b).
--
-- The research DB (kalshi_mlb_rfq_research.duckdb) is a separate sibling of the
-- trading DB. Run from the kalshi_mlb_rfq/ directory, e.g.:
--   duckdb < research_queries.sql
-- or open an interactive shell and ATTACH both DBs:
--   duckdb -cmd "ATTACH 'kalshi_mlb_rfq.duckdb' AS state (READ_ONLY); ATTACH 'kalshi_mlb_rfq_research.duckdb' AS research (READ_ONLY);"
--
-- Payload fields are JSON. Use json_extract_string(payload, 'key') to read
-- them, NOT the ->> operator: in DuckDB 1.4.x, ->> can try to cast the whole
-- payload to a number when the same query both extracts a field in SELECT and
-- compares one in WHERE over rows whose JSON holds numeric values.
-- json_extract_string always returns VARCHAR and sidesteps that.
--
-- Join keys:
--   * candidate_evaluated / kelly_sized / rfq_submit_failed are pre-RFQ, so
--     combo_ticker is NULL — join them by json_extract_string(payload,'leg_set_hash').
--   * quote_priced / gate_evaluated / walk_diagnosed / position_snapshot are
--     post-RFQ — they carry combo_ticker, rfq_id, quote_id.
--   * state.live_rfqs maps leg_set_hash <-> rfq_id <-> combo_market_ticker.

ATTACH 'kalshi_mlb_rfq.duckdb'          AS state    (READ_ONLY);
ATTACH 'kalshi_mlb_rfq_research.duckdb' AS research (READ_ONLY);

-- ---------------------------------------------------------------------------
-- 1) RFQ JOURNEY — replay every event for one combo, in time order.
--    Set :lsh to a leg_set_hash (from candidate_evaluated or state.live_rfqs).
--    Captures both the pre-RFQ events (matched by leg_set_hash in payload) and
--    the post-RFQ events (matched by combo_ticker via the live_rfqs mapping).
-- ---------------------------------------------------------------------------
-- SELECT e.ts, e.event_type, e.combo_ticker, e.payload
-- FROM research.events e
-- WHERE json_extract_string(e.payload, 'leg_set_hash') = :lsh
--    OR e.combo_ticker IN (
--        SELECT combo_market_ticker FROM state.live_rfqs WHERE leg_set_hash = :lsh)
-- ORDER BY e.ts;

-- ---------------------------------------------------------------------------
-- 2) MISSED EDGES — candidates we priced but did NOT submit (last 200).
--    "Where did we leave the surface on the table, and why?"
-- ---------------------------------------------------------------------------
SELECT ts,
       game_id,
       json_extract_string(payload, 'leg_set_hash') AS leg_set_hash,
       json_extract_string(payload, 'outcome')      AS outcome,
       json_extract_string(payload, 'model_fair')   AS model_fair,
       json_extract_string(payload, 'blended_fair') AS blended_fair,
       json_extract_string(payload, 'n_books')      AS n_books
FROM research.events
WHERE event_type = 'candidate_evaluated'
  AND json_extract_string(payload, 'outcome') <> 'submitted'
ORDER BY ts DESC
LIMIT 200;

-- ---------------------------------------------------------------------------
-- 3) GATE BREAKDOWN — which gate rejected accepts over the last 24h.
-- ---------------------------------------------------------------------------
SELECT json_extract_string(payload, 'decision') AS decision,
       count(*)                                 AS n
FROM research.events
WHERE event_type = 'gate_evaluated'
  AND ts > now() - INTERVAL 1 DAY
GROUP BY 1
ORDER BY n DESC;
