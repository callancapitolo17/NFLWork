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
-- The surface DB only exists once the leg-surface ingest (#96) has run at
-- least once. If you have never run it, comment this line out along with
-- queries 12-16 — a READ_ONLY ATTACH of a missing file aborts the script.
ATTACH 'kalshi_mlb_mm_surface.duckdb'  AS surface  (READ_ONLY);

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
-- 9) LIVE-PRICING HEALTH (issue #54): daily veto/void rate. Every quote is
--    priced from a post-RFQ live fetch, so the confirm last look re-checks a
--    SECONDS-old number — the void rate should collapse vs the pre-#54
--    cache-era days in this same table. If it does not, live pricing is not
--    actually working (books timing out, engine wedged, warming broken).
--    live_fetch_timeout / live_too_few_books daily counts ride along as the
--    "why aren't we quoting" companions.
-- ---------------------------------------------------------------------------
SELECT CAST(observed_at AS DATE)                                       AS day,
       SUM(CASE WHEN decision LIKE 'voided_%' THEN 1 ELSE 0 END)        AS voided,
       SUM(CASE WHEN decision = 'confirmed'   THEN 1 ELSE 0 END)        AS confirmed,
       SUM(CASE WHEN decision LIKE 'voided_%' THEN 1 ELSE 0 END) * 1.0
         / NULLIF(SUM(CASE WHEN decision LIKE 'voided_%' THEN 1 ELSE 0 END)
                  + SUM(CASE WHEN decision = 'confirmed' THEN 1 ELSE 0 END), 0)
                                                                        AS void_rate,
       SUM(CASE WHEN reason = 'live_fetch_timeout'  THEN 1 ELSE 0 END)  AS live_fetch_timeouts,
       SUM(CASE WHEN reason = 'live_too_few_books'  THEN 1 ELSE 0 END)  AS live_too_few_books
FROM state.quote_decisions
GROUP BY 1
ORDER BY 1;

-- ---------------------------------------------------------------------------
-- 10) LIVE FETCH TRACE (issue #54 acceptance): per-quote per-game live books,
--     latencies, and result age at pricing time — proof every quoted fair
--     traces to a fetch initiated after its RFQ landed (join the
--     on_demand_requested / on_demand_result timestamps by leg_set_hash for
--     the full chain). Novig caveat: NV fairs inherit its open sweep-price
--     broadcast defect — check per-book values before leaning on NV-heavy
--     consensus.
-- ---------------------------------------------------------------------------
SELECT ts,
       ticker,
       json_extract_string(payload, 'blended_fair')  AS blended_fair,
       json_extract_string(payload, 'live_games')    AS live_games
FROM research.events
WHERE event_type = 'quote_priced'
ORDER BY ts DESC;

-- ---------------------------------------------------------------------------
-- 11) EXPIRED-QUOTE OUTCOMES — the margin-tuning ratio. Of quotes that
--     expired unfilled, how many saw the combo market trade during our
--     resting window (competitor_traded_in_window = a competing maker won
--     the RFQ, we're being outpriced) vs no trade at all (no_trade = the
--     creator never executed with anyone — cutting margin donates edge)?
--     traded_shortly_after is the near-miss bucket in between.
--     Uses state.quote_expiry_outcomes (written by the expiry-outcome sweep).
-- ---------------------------------------------------------------------------
SELECT CAST(window_end AS DATE)  AS day,
       label,
       COUNT(*)                  AS quotes,
       ROUND(COUNT(*) * 1.0 / SUM(COUNT(*)) OVER (
           PARTITION BY CAST(window_end AS DATE)), 3) AS share_of_day
FROM state.quote_expiry_outcomes
GROUP BY 1, 2
ORDER BY 1, 2;


-- ---------------------------------------------------------------------------
-- 12) LEG SURFACE: ACHIEVED CADENCE per book (issue #96 acceptance).
--     "built_at freshness for each book stays within its target cadence" is
--     really two numbers: how long a pass TAKES and how often it STARTS. A
--     pass that overruns its cadence pushes every row past #99's age gate, so
--     compare med_gap against SURFACE_CADENCE_*_SEC and med_dur against it.
--     A gap far BELOW the cadence means the sleep is broken, not that the
--     book is fast — that bug shipped once (Event.wait on a set flag).
-- ---------------------------------------------------------------------------
SELECT book,
       route,
       COUNT(*)                         AS passes,
       ROUND(MEDIAN(duration_sec), 1)   AS med_dur_sec,
       ROUND(MEDIAN(gap), 1)            AS med_gap_sec,
       ROUND(MIN(gap), 1)               AS min_gap_sec,
       ROUND(MAX(duration_sec), 1)      AS worst_dur_sec
FROM (SELECT book, route, duration_sec,
             epoch(started_at - LAG(started_at) OVER (
                 PARTITION BY book, route ORDER BY started_at)) AS gap
      FROM surface.surface_refresh_log)
GROUP BY 1, 2
ORDER BY 1, 2;

-- ---------------------------------------------------------------------------
-- 13) LEG SURFACE: EXCLUSIONS BY REASON per book. crossed/overround should be
--     ~0 (they are feed-integrity guards, not filters). unresolved is EXPECTED
--     and time-of-day dependent — books shrink their ladders far from first
--     pitch — so it is not a health signal. A non-zero game_ambiguous means a
--     doubleheader was declined fail-closed, which is the correct outcome and
--     worth knowing about.
-- ---------------------------------------------------------------------------
SELECT book,
       route,
       SUM(rungs_priced)      AS priced,
       SUM(n_crossed)         AS crossed,
       SUM(n_overround)       AS overround,
       SUM(n_one_sided)       AS one_sided,
       SUM(n_unresolved)      AS unresolved,
       SUM(n_game_unmatched)  AS game_unmatched,
       SUM(n_game_ambiguous)  AS game_ambiguous,
       COUNT(*) FILTER (WHERE error_class IS NOT NULL) AS failed_passes
FROM surface.surface_refresh_log
GROUP BY 1, 2
ORDER BY 1, 2;

-- ---------------------------------------------------------------------------
-- 14) LEG SURFACE: BOOKS PER LEG — the quorum question. #98 needs
--     MIN_AGREEING_BOOKS (2) per leg, so the mass at n_books = 1 is the share
--     of legs that will decline no matter how fresh the surface is.
-- ---------------------------------------------------------------------------
SELECT n_books, COUNT(*) AS legs
FROM (SELECT game_id, period, market_type, line, side,
             COUNT(DISTINCT book) AS n_books
      FROM surface.mlb_leg_surface
      GROUP BY ALL)
GROUP BY 1
ORDER BY 1;

-- ---------------------------------------------------------------------------
-- 15) LEG SURFACE: DUPLICATE KEYS — must be zero (issue #96 acceptance).
--     Uniqueness is structural (each flush replaces one (book, route) slice)
--     rather than a constraint, because `line` is legitimately NULL for
--     moneyline and DuckDB PRIMARY KEYs reject NULL. This is the check that
--     caught FanDuel's two routes both pricing its FG totals.
-- ---------------------------------------------------------------------------
SELECT book, game_id, period, market_type, line, side, COUNT(*) AS rows_
FROM surface.mlb_leg_surface
GROUP BY ALL
HAVING COUNT(*) > 1
ORDER BY 1, 2;

-- ---------------------------------------------------------------------------
-- 16) LEG SURFACE: ROW AGE right now, per book. This is exactly what the
--     SURFACE_MAX_AGE_SEC gate (#99) sees. DraftKings is expected to sit at
--     30-90s (its slate scrape alone is ~21-28s), so the 30s gate removes it
--     from consensus — deliberately, and counted in query 20.
-- ---------------------------------------------------------------------------
SELECT book,
       route,
       COUNT(*)                                       AS legs,
       COUNT(DISTINCT game_id)                        AS games,
       ROUND(MEDIAN(epoch(now() - built_at)), 1)      AS med_age_sec,
       ROUND(MAX(epoch(now() - built_at)), 1)         AS max_age_sec
FROM surface.mlb_leg_surface
GROUP BY 1, 2
ORDER BY 1, 2;

-- ---------------------------------------------------------------------------
-- 17) #98 ROUTING SPLIT — how many quotes came from the surface vs the live
--     engine, and how many were MIXED. This is the throughput claim: a
--     surface-only quote cost zero book requests. Cross-check against the
--     book_requests_per_day query in kalshi_common/fetch_health_queries.sql —
--     post-#98, on_demand rows should track the same-game half only.
-- ---------------------------------------------------------------------------
SELECT CASE
         WHEN json_extract_string(payload, 'surface_games') IS NULL
           THEN 'live_only'
         WHEN json_extract_string(payload, 'live_games') IS NULL
           THEN 'surface_only'
         ELSE 'mixed'
       END                                            AS pricing_source,
       COUNT(*)                                       AS quotes,
       ROUND(AVG(CAST(json_extract_string(payload, 'blended_fair')
                      AS DOUBLE)), 4)                 AS avg_fair
FROM research.events
WHERE event_type = 'quote_priced'
GROUP BY 1
ORDER BY 2 DESC;

-- ---------------------------------------------------------------------------
-- 18) THE THRESHOLD QUERY (#99) — the age of the surface rows that ACTUALLY
--     backed a quote, per book. Query 16 shows what the surface HOLDS, this
--     shows what got USED, and query 20 shows what the gate REFUSED.
--     Re-read the percentiles per book before changing SURFACE_MAX_AGE_SEC or
--     any ingest cadence: a gate below a book's p95 removes that book from
--     consensus most of the time. Post-#99 every age here is by construction
--     <= SURFACE_MAX_AGE_SEC.
-- ---------------------------------------------------------------------------
WITH used AS (
    SELECT e.ts,
           game.key                                        AS leg_set_hash,
           bk.key                                          AS book,
           CAST(json_extract_string(bk.value, 'age_sec') AS DOUBLE) AS age_sec,
           json_extract_string(bk.value, 'route')          AS route
    FROM research.events                                    AS e,
         LATERAL json_each(json_extract(e.payload, '$.surface_games')) AS game,
         LATERAL json_each(json_extract(game.value, '$.books'))        AS bk
    WHERE e.event_type = 'quote_priced'
      AND json_extract_string(e.payload, 'surface_games') IS NOT NULL
)
SELECT book,
       route,
       COUNT(*)                          AS legs_used,
       ROUND(MEDIAN(age_sec), 1)         AS p50_age_sec,
       ROUND(QUANTILE_CONT(age_sec, 0.95), 1) AS p95_age_sec,
       ROUND(MAX(age_sec), 1)            AS max_age_sec
FROM used
GROUP BY 1, 2
ORDER BY 5 DESC;

-- ---------------------------------------------------------------------------
-- 19) #98/#99 DECLINE ATTRIBUTION — a thin CACHED surface vs a thin LIVE
--     fetch, and WHY the surface was thin. These are the same #20 gate on
--     different inputs and they call for opposite fixes:
--       surface_too_few_books    -> ingest COVERAGE (the matrix in config.py)
--       surface_stale            -> ingest CADENCE, or a book gone dark (#99)
--       surface_dispersion       -> the fresh books genuinely disagree
--       surface_constituent_moved-> Kalshi moved after the row was built (#99)
--       live_* / on_demand_*     -> book latency or the on-demand path
-- ---------------------------------------------------------------------------
SELECT reason, COUNT(*) AS n
FROM state.quote_decisions
WHERE reason IN ('surface_too_few_books', 'surface_stale',
                 'surface_dispersion', 'surface_constituent_moved',
                 'live_too_few_books', 'consensus_dispersion',
                 'live_fetch_timeout', 'on_demand_pending')
  AND observed_at >= now() - INTERVAL 24 HOUR
GROUP BY 1
ORDER BY 2 DESC;

-- ---------------------------------------------------------------------------
-- 20) #99 AGE-GATE EXCLUSIONS, per book. The companion to query 18: that one
--     measures the rows that PRICED a quote, this one measures the rows the
--     gate refused. DraftKings is EXPECTED at ~100% excluded (its slate
--     scrape alone is 21-28s against a 30s gate) and that is by design —
--     announced by a startup WARNING and counted here. ANY OTHER BOOK
--     drifting up is the signal that the ingest cadence and
--     SURFACE_MAX_AGE_SEC have fallen out of step.
--     Rows are cumulative-per-window; sum them, do not diff them.
-- ---------------------------------------------------------------------------
WITH windows AS (
    SELECT bk.key                                             AS book,
           CAST(json_extract_string(bk.value, 'used') AS BIGINT)    AS used,
           CAST(json_extract_string(bk.value, 'excluded_by_age')
                AS BIGINT)                                    AS excluded
    FROM research.events                                       AS e,
         LATERAL json_each(json_extract(e.payload, '$.books'))  AS bk
    WHERE e.event_type = 'surface_age_summary'
      AND e.ts >= now() - INTERVAL 24 HOUR
)
SELECT book,
       SUM(used)                                              AS legs_used,
       SUM(excluded)                                          AS legs_excluded,
       ROUND(100.0 * SUM(excluded)
             / NULLIF(SUM(used) + SUM(excluded), 0), 1)       AS pct_excluded
FROM windows
GROUP BY 1
ORDER BY 4 DESC;

-- ---------------------------------------------------------------------------
-- 21) #99 CONSTITUENT VETO — is the second guard earning its place, and is
--     SURFACE_CONSTITUENT_MOVE_THRESHOLD (0.03, inherited from #23) right?
--     Read three things:
--       * the verdict mix — a large `no_baseline` share means the tape is not
--         seeing these legs often enough for the veto to bite at all;
--       * the delta percentiles among 'ok' verdicts — if p95 sits far below
--         the threshold the guard is nearly inert and could tighten;
--       * `moved` volume — this is edge we are declining to quote, so it
--         should be a small minority, not the common case.
-- ---------------------------------------------------------------------------
WITH legs AS (
    SELECT json_extract_string(l.value, 'verdict')                  AS verdict,
           CAST(json_extract_string(l.value, 'delta') AS DOUBLE)     AS delta,
           CAST(json_extract_string(l.value, 'oldest_row_age_sec')
                AS DOUBLE)                                          AS row_age_sec,
           CAST(json_extract_string(l.value, 'baseline_age_sec')
                AS DOUBLE)                                          AS baseline_age_sec
    FROM research.events                                            AS e,
         LATERAL json_each(json_extract(e.payload, '$.legs'))       AS l
    WHERE e.event_type = 'surface_constituent_check'
      AND e.ts >= now() - INTERVAL 24 HOUR
)
SELECT verdict,
       COUNT(*)                                       AS n,
       ROUND(100.0 * COUNT(*) / SUM(COUNT(*)) OVER (), 1) AS pct,
       ROUND(MEDIAN(delta), 4)                        AS p50_delta,
       ROUND(QUANTILE_CONT(delta, 0.95), 4)           AS p95_delta,
       ROUND(MEDIAN(row_age_sec), 1)                  AS p50_row_age_sec,
       ROUND(MEDIAN(baseline_age_sec), 1)             AS p50_baseline_lead_sec
FROM legs
GROUP BY 1
ORDER BY 2 DESC;
