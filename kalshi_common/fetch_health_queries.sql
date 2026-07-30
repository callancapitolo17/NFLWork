-- Ready-made reads over sgp_fetch_health (issue #38).
--
-- Run against a bot's MARKET DB (the same file mlb_sgp_odds lives in):
--     kalshi_mlb_mm/kalshi_mlb_mm_market.duckdb      (maker)
--     kalshi_mlb_rfq/kalshi_mlb_rfq_market.duckdb    (taker)
--
--   duckdb kalshi_mlb_mm/kalshi_mlb_mm_market.duckdb -readonly
--
-- Queries are delimited by `-- name: <id>` headers; tests select by name.
--
-- ####################################################################
-- ##  UNITS WARNING — ALWAYS GROUP BY path, NEVER AVERAGE ACROSS IT  ##
-- ####################################################################
-- counters_json.targets_attempted and .prices_returned carry DIFFERENT
-- UNITS on the two paths (see mlb_sgp/_shared.py::COUNTER_NAMES):
--
--   sweep      targets_attempted = TARGET LINES priced (each fans out to
--                                  ~4 combo price calls)
--              prices_returned   = PricedRows produced
--   on_demand  targets_attempted = PRICE CALLS issued (2^N partition cells)
--              prices_returned   = those calls that returned a decimal
--
-- The sweep's ratio is ~1:4 BY CONSTRUCTION and on-demand's is 1:1. A
-- figure pooled across both paths is not a weaker signal — it is a wrong
-- one. Same for duration: a 40s sweep and a 2s live fetch share no scale.


-- name: outcome_mix_24h
-- Per book x path: last-24h outcome mix and latency spread.
-- This is the "is the feed alive?" readout — start here.
SELECT
    book,
    path,
    count(*)                                              AS n_fetches,
    count(*) FILTER (WHERE outcome = 'ok')                AS n_ok,
    count(*) FILTER (WHERE outcome = 'empty')             AS n_empty,
    count(*) FILTER (WHERE outcome = 'transport_error')   AS n_transport_error,
    count(*) FILTER (WHERE outcome = 'timeout')           AS n_timeout,
    count(*) FILTER (WHERE outcome = 'error')             AS n_error,
    round(100.0 * count(*) FILTER (WHERE outcome = 'ok') / count(*), 1)
                                                          AS pct_ok,
    round(quantile_cont(duration_sec, 0.50), 2)           AS p50_duration_sec,
    round(quantile_cont(duration_sec, 0.95), 2)           AS p95_duration_sec,
    max(fetched_at)                                       AS last_fetch_at
FROM sgp_fetch_health
WHERE fetched_at >= now() - INTERVAL 24 HOUR
GROUP BY book, path
ORDER BY book, path;


-- name: current_failure_streak
-- How many consecutive fetches a book's most recent run of failures spans,
-- per path. This is what #37 alerts on: one 403 is noise, forty in a row is
-- a dead book. A book with a fresh 'ok' reports streak 0.
WITH ordered AS (
    SELECT book, path, fetched_at, outcome, error_class,
           row_number() OVER (PARTITION BY book, path
                              ORDER BY fetched_at DESC) AS rn
    FROM sgp_fetch_health
    WHERE fetched_at >= now() - INTERVAL 7 DAY
),
last_ok AS (
    SELECT book, path, min(rn) AS rn_last_ok
    FROM ordered
    WHERE outcome IN ('ok', 'empty')
    GROUP BY book, path
)
SELECT
    o.book,
    o.path,
    coalesce(l.rn_last_ok - 1, max(o.rn))       AS failing_streak,
    max(o.fetched_at) FILTER (WHERE o.rn = 1)   AS last_fetch_at,
    arg_max(o.error_class, -o.rn)               AS newest_error_class
FROM ordered o
LEFT JOIN last_ok l USING (book, path)
GROUP BY o.book, o.path, l.rn_last_ok
HAVING failing_streak > 0
ORDER BY failing_streak DESC;


-- name: parse_health_24h
-- The #35 counters, kept strictly inside one path (see the units warning).
-- A book whose legs_resolved collapses while transport stays clean is the
-- FanDuel alt-total regression: healthy wire, blind parser.
SELECT
    book,
    path,
    count(*)                                                  AS n_fetches,
    sum(CAST(json_extract(counters_json, '$.events_seen')     AS BIGINT)) AS events_seen,
    sum(CAST(json_extract(counters_json, '$.events_matched')  AS BIGINT)) AS events_matched,
    sum(CAST(json_extract(counters_json, '$.legs_attempted')  AS BIGINT)) AS legs_attempted,
    sum(CAST(json_extract(counters_json, '$.legs_resolved')   AS BIGINT)) AS legs_resolved,
    sum(CAST(json_extract(counters_json, '$.targets_attempted') AS BIGINT)) AS targets_attempted,
    sum(CAST(json_extract(counters_json, '$.prices_returned') AS BIGINT)) AS prices_returned,
    sum(CAST(json_extract(counters_json, '$.parse_failures')  AS BIGINT)) AS parse_failures,
    sum(CAST(json_extract(counters_json, '$.sanity_drops')    AS BIGINT)) AS sanity_drops
FROM sgp_fetch_health
WHERE fetched_at >= now() - INTERVAL 24 HOUR
  AND counters_json IS NOT NULL
GROUP BY book, path
ORDER BY book, path;


-- name: retention_footprint
-- Sanity check on the 30-day prune: oldest row, row count, days covered.
SELECT
    count(*)                                         AS n_rows,
    min(fetched_at)                                  AS oldest,
    max(fetched_at)                                  AS newest,
    date_diff('day', min(fetched_at), max(fetched_at)) AS days_covered
FROM sgp_fetch_health;
