# Kalshi parlay quoter: phase 1, dry run (#283)

Part of epic #281. Watches every Kalshi RFQ, decides whether it is a parlay
we would quote, prices it from Unabated fairs, and logs the decision. **It
sends no quotes and no orders**; nothing in this package can.

## How it works

1. **Feed** (`feed.py`). One thread reads Kalshi's `communications`
   websocket with the `websockets` library: ~9,900 frames a second, half
   `rfq_created`, half `rfq_deleted` (2026-11-13). The library is the point:
   `websocket-client` spent ~76 us a frame, kept up with ~4,200 frames a
   second at 93% of a core, and Kalshi cut the socket every 15-20 s.
2. **Door** (`feed.door`, `scope.py`). A byte check drops any frame naming
   none of our series before it is decoded. The rest must pass the scope
   rule: every leg a full-game moneyline (`…GAME`), spread (`…SPREAD`) or
   total (`…TOTAL`) in NFL (`KXNFL`), CFB (`KXNCAAF`), NBA (`KXNBA`), CBB
   (`KXNCAAMB`) or MLB (`KXMLB`); at least two legs; no two legs on the same
   game; and the size floor ($50, or 50 contracts when sized in contracts).
   Period markets (1H, F5) and props fall out as `unsupported_leg`. About
   475 RFQs a second pass. The rest are counted, never stored.
3. **Leg table** (`legtable.py`). Every leg seen in the last 10 minutes is
   held in memory with its fair from `GET 127.0.0.1:8095/rfq_fairs.json`
   (200 legs a call, board only, no Kalshi read), and an RFQ is priced from
   that table, so the work grows with distinct legs (thousands), not with
   RFQs. **A leg is priced when it changes, not on a clock**: its fair moves
   only when Unabated's snapshot for its league lands on the runner, about
   once a minute. Every 2 s the quoter reads each league's snapshot time off
   the runner's `/health` (answered by the runner's main thread, at no cost
   to the thread its taker runs on) and re-prices only the legs of a league
   whose time moved, and any leg priced over 90 s ago. The first version
   re-priced every leg every 2 s: 30 times the runner calls for the same
   numbers, on the taker's thread (2026-11-13, the taker froze 40 minutes
   after the quoter went live).
   Per leg: Unabated fair for that side, `fair_lo`/`fair_hi` bracket, `sigma`
   (with `c`, `tau`), the resolved Unabated line, kickoff and fair age. A leg
   seen for the first time is priced within ~0.25 s; its RFQ waits for
   that, at most 5 s. A leg the runner refuses (HTTP
   400) is isolated by halving its chunk, declines `runner_rejected`, and is
   asked again after 5 minutes (so a ticker that stays bad counts in
   `table_legs_rejected` 12 times an hour).
4. **Decision** (`decide.py`). `would_quote` when every leg is priced (parlay
   fair = product of leg fairs, since legs are on different games), otherwise
   `decline` with a reason.
5. **Confirm (off by default since 2026-11-13).** With
   `PARLAY_CONFIRMS_PER_SEC` above 0, that many would-quote parlays a second
   are re-priced through `GET /rfq_price.json`, where the runner also reads Kalshi's market (status,
   strike), book and start for each leg. Those rows are `confirmed` and carry
   Kalshi's prices. The cap exists because those reads share Kalshi's 20
   reads/s limit with the taker bot. The others are counted
   `would_quote_unconfirmed`. A confirm that fails with `pricer_error` stops
   confirms for 30 s. On the VM most confirms timed out (16 confirmed, 288
   timed out in 80 minutes) and their Kalshi reads run on the taker's
   thread, so they stay off until they read the runner's Kalshi book stream.

Python never sends keys to the runner; the runner never sees an RFQ beyond
its legs. The margin formula is #284; a confirmed row logs everything it
needs.

## Decision log (`kalshi_parlay_quoter.duckdb`)

| table | one row per | notes |
|---|---|---|
| `minute_counts` | minute and key | every frame and every decision, stored or not. Totals come from here, never from counting rows |
| `rfq_decisions` | stored RFQ | every `confirmed` parlay, plus 0.5% of the rest (a fixed hash of the RFQ id). Size, leagues, leg count, decision, reason, parlay fair and bracket, max leg sigma, product of Kalshi's single-leg asks (confirmed only), hours to kickoff, `queue_ms` (seen to decided) |
| `rfq_legs` | leg of a stored RFQ | what the pricer returned (`line` is JSON, `detail` is the runner's note), `priced_at_ms` and `cache_age_ms`. Kalshi columns are null unless the RFQ is confirmed |

`minute_counts` keys:

- Frames: `frames_deleted`, `frames_other`, `malformed` (an `rfq_created`
  frame that could not be read), `reconnects`.
- Door: `door_other_league` (no leg in our series; legless frames too),
  `door_unsupported_leg`, `door_same_game`, `door_single_leg`,
  `door_bad_side`, `door_no_legs`, `door_below_floor`, `kept`,
  `queue_dropped`.
- Decisions, one per kept RFQ: `would_quote` (confirmed),
  `would_quote_unconfirmed`, `pricer_error`, `pricer_backlog`,
  `leg_<runner reason>` (`leg_past_ladder`, `leg_started`, `leg_stale_feed`,
  `leg_no_event`, `leg_runner_rejected`), `confirm_<reason>` (a would-quote
  the Kalshi-checked confirm declined: `confirm_leg_market_mismatch`,
  `confirm_leg_kalshi_unread`, `confirm_pricer_error`, …),
  `expired_in_queue` (waited over 15 s before it was decided),
  `dropped_on_pause`, `handle_crashed`, `confirm_crashed`.
- `stored`: rows written to `rfq_decisions`. It is on top of the decision
  keys, not one of them.
- Health: `table_refresh_failed`, `table_legs_rejected`, `table_legs_priced`
  (leg prices asked of the runner: about one per table leg a minute),
  `gauge_table_legs`,
  `gauge_table_refresh_ms`, `gauge_pending`, `gauge_queue` (gauges are the
  value at the minute's end, not a count).

`rfq_created` frames a minute = `malformed` + every `door_*` key + `kept` +
`queue_dropped`. `kept` = the decision keys, less RFQs still waiting at the
minute's end (`gauge_pending`, `gauge_queue`, and up to 2 confirms in
flight, which count in the minute they finish). A restart inside a minute adds
to that minute's counts.

Rows from before 2026-11-13 are a census (every in-scope RFQ, each priced
with Kalshi reads) and have `confirmed` NULL. Do not compare row counts
across that date; use `minute_counts` after it.

```sql
-- the 24 h check
SELECT key, sum(n) FROM minute_counts
WHERE minute > now() - INTERVAL 24 HOUR AND key NOT LIKE 'gauge_%'
GROUP BY ALL ORDER BY 2 DESC;
-- confirmed would-quotes against Kalshi's own leg prices
SELECT count(*), median(parlay_fair / legs_kalshi_ask_product)
FROM rfq_decisions WHERE confirmed AND decision = 'would_quote'
  AND seen_at > now() - INTERVAL 24 HOUR;
```

## Freshness

A fair is priced on only while all three hold:

- the runner's snapshot times were read under 15 s ago
  (`PARLAY_TABLE_MAX_AGE_SEC`). A runner that is down stops every fair
  within 15 s. A runner whose board thread is stuck keeps answering
  `/health` ok for 15 s more (it says not ok after 15 s without a post from
  that thread), so that case takes about 30-35 s;
- the leg was priced on its league's current snapshot, or the quoter saw
  that snapshot change under 15 s ago and the re-price is on its way. The
  quoter sees a change up to ~7 s after it lands (`/health` is posted every
  5 s and read every 2 s). A league the runner would decline (stale or
  failed snapshot, not loaded, a restarted runner) is left out of the
  snapshot times, which counts as a change;
- the leg was priced under 105 s ago (`PARLAY_TABLE_MAX_REUSE_SEC` + 15).

So a fair can be one snapshot behind the runner's board for about 22 s. What
a snapshot time does not show is caught by the 90 s re-price and can lag
that long: time to kickoff inside `sigma`/`tau`, and the runner's tail
measurement. Every leg row records `cache_age_ms` (now often 30-60 s, since
an unchanged fair is reused), and `fair_age_ms` is the age of Unabated's
snapshot plus that
reuse. A leg whose kickoff has passed since it was priced logs `started`.
With the runner down, RFQs decline `pricer_error`.

## Running

On the VM it sits behind the `parlay_quoter` compose profile. To ship it
with every `deploy.sh` run (the one deploy route), set
`COMPOSE_PROFILES=parlay_quoter` in `unabated_ticket/deploy/.env`; deploy.sh
then builds it, restarts it and checks it is running. Remove the line (and
`docker compose stop parlay_quoter`) to turn it off.

```
cd ~/NFLWork/unabated_ticket/deploy
docker compose logs -f parlay_quoter          # one counts line a minute
touch ~/unabated-data/parlay_quoter.stop      # close the feed, stop all runner calls, drop waiting RFQs (rm to resume)
```

Load, measured locally 2026-11-13 against a fake feed at 9,900 frames a
second with 475 kept RFQs a second and 9,400 table legs: 31% of one core in
total, no queue, no reconnects (no TLS; a live read-only run of the reader
alone cost 42% of a core with TLS), and 9,360 leg prices asked of the runner
a minute (`table_legs_priced`; the every-2-s version asked ~280,000). The log grows ~190 MB a day at 3.4
stored decisions a second. Rows are inserted in multi-row statements
(DuckDB's `executemany` goes row by row and was the first VM run's CPU sink)
and the file is checkpointed each minute.

It signs with the same read-only Kalshi key as the bets service
(`KALSHI_API_KEY_ID`, read from `bets_service/.env` and nothing else from
that file) and writes only `/data/kalshi_parlay_quoter.duckdb`.

Locally: `KALSHI_API_KEY_ID=… KALSHI_PRIVATE_KEY_PATH=… python -m kalshi_parlay_quoter.main`.

| env | default | |
|---|---|---|
| `PARLAY_QUOTER_DB_PATH` | package dir | decision log |
| `UNABATED_RUNNER_URL` | `http://127.0.0.1:8095` | leg pricers |
| `KALSHI_WS_URL` | Kalshi production | RFQ feed |
| `PARLAY_MIN_RFQ_TARGET_COST_USD` | 50 | size floor at the door; 0 keeps every size |
| `PARLAY_TABLE_REFRESH_SEC` | 2 | how often the runner's snapshot times are read off `/health` |
| `PARLAY_TABLE_HEALTH_TIMEOUT_SEC` | 3 | one `/health` read |
| `PARLAY_TABLE_MAX_REUSE_SEC` | 90 | a leg is re-priced at the latest this long after it was priced |
| `PARLAY_TABLE_MAX_AGE_SEC` | 15 | how old the snapshot times may be, and how long an old fair is served after its league changed |
| `PARLAY_TABLE_LEG_TTL_SEC` | 600 | how long a leg stays in the table after it was last asked |
| `PARLAY_TABLE_FETCH_TIMEOUT_SEC` | 15 | one refresh call |
| `PARLAY_PENDING_MAX_WAIT_SEC` | 5 | wait for a leg's first fair, then `pricer_error` |
| `PARLAY_MAX_RFQ_AGE_SEC` | 15 | an RFQ that waited longer in the queue is counted `expired_in_queue`, not decided |
| `PARLAY_MAX_BACKLOG` | 2000 | RFQs waiting that way; beyond it `pricer_backlog` |
| `PARLAY_CONFIRMS_PER_SEC` | 0 | Kalshi-checked confirms a second; 0 is off |
| `PARLAY_CONFIRM_WORKERS` | 2 | confirms in flight |
| `PARLAY_CONFIRM_TIMEOUT_SEC` | 3 | slower = `confirm_pricer_error` (detail `timeout`) |
| `PARLAY_CONFIRM_BACKOFF_SEC` | 30 | no confirms for this long after one fails with `pricer_error` |
| `PARLAY_STORED_SAMPLE_PCT` | 0.5 | share of unconfirmed RFQs stored as rows |
| `PARLAY_FEED_QUEUE_MAX` | 20000 | kept RFQs waiting to be decided; beyond it `queue_dropped` |
| `PARLAY_FEED_IDLE_TIMEOUT_SEC` | 30 | reconnect after this long with no frame |
| `PARLAY_FLUSH_SEC` | 1 | how often rows are written |
| `PARLAY_KALSHI_ENV_FILE` | unset | file to read `KALSHI_API_KEY_ID` from when it is not in the environment |
| `KALSHI_BASE_URL` | Kalshi production | the signer requires it; no REST call is made |
| `PARLAY_STOP_PATH` | package dir `.stop` | pause file |

Tests: `python -m pytest kalshi_parlay_quoter/tests`.
