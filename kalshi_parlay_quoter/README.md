# Kalshi parlay quoter: phase 1, dry run (#283)

Part of epic #281. Watches every Kalshi RFQ, decides whether it is a parlay
we would quote, prices it through the runner, and logs the decision. **It
sends no quotes and no orders**; nothing in this package can.

## How it works

1. **Feed.** `kalshi_mlb_mm.ws_rfq_source.WebSocketRFQSource` (the MLB
   maker's battle-tested WS mirror with loud REST fallback), with its door
   widened to the leagues below and no leg cap.
2. **Scope rule** (`scope.py`). Every leg is a full-game moneyline
   (`…GAME`), spread (`…SPREAD`) or total (`…TOTAL`) series in NFL
   (`KXNFL`), CFB (`KXNCAAF`), NBA (`KXNBA`), CBB (`KXNCAAMB`) or MLB
   (`KXMLB`); at least two legs; no two legs on the same game. Period markets
   (1H, F5) and props are declined as `unsupported_leg`.
3. **Pricing.** In-scope RFQs go to the runner (#282):
   `GET 127.0.0.1:8095/rfq_price.json?rfq_id=…&legs=TICKER:side,…`. Per leg:
   Unabated fair for that side, `fair_lo`/`fair_hi` bracket, `sigma` (with
   `c`, `tau`), the resolved Unabated line, kickoff, fair age and Kalshi's
   live price. Python never sends keys to the runner; the runner never sees
   an RFQ beyond its legs.
4. **Decision.** `would_quote` when every leg is priced (parlay fair =
   product of leg fairs, since legs are on different games), otherwise
   `decline` with a reason. The margin formula is #284; everything it needs
   is logged so margins can be replayed over this data.

## Decision log (`kalshi_parlay_quoter.duckdb`)

| table | one row per | notes |
|---|---|---|
| `rfq_decisions` | in-league RFQ id | size, leagues, leg count, same-game flag, decision, reason, parlay fair and bracket, max leg sigma, product of Kalshi's single-leg asks, hours to kickoff, pricer latency, `queue_ms` (wait for a free pricer worker) |
| `rfq_legs` | leg of an RFQ that reached the pricer | everything the pricer returned for the leg (`line` is JSON, `detail` is the runner's note) |
| `door_counts` | minute | RFQs that got no row: `other_league`, `below_floor` (size under `PARLAY_MIN_RFQ_TARGET_COST_USD`, default $50), `too_many_legs` (always 0, no cap), `malformed`; plus `logged`, the pricer `backlog` and the feed mode |

Every RFQ seen is either a row in `rfq_decisions` or a count in
`door_counts`. Decline reasons: `no_legs`, `unsupported_leg`, `bad_side`,
`single_leg`, `same_game`, `pricer_error`, `pricer_backlog`, and `leg_<runner reason>`
(`leg_past_ladder`, `leg_started`, `leg_stale_feed`, …).

```sql
-- the 24 h check
SELECT decision, reason, count(*) FROM rfq_decisions
WHERE seen_at > now() - INTERVAL 24 HOUR GROUP BY ALL ORDER BY 3 DESC;
SELECT sum(logged), sum(other_league), sum(below_floor) FROM door_counts
WHERE minute > now() - INTERVAL 24 HOUR;
```

## Running

On the VM it sits behind the `parlay_quoter` compose profile. To ship it
with every `deploy.sh` run (the one deploy route), set
`COMPOSE_PROFILES=parlay_quoter` in `unabated_ticket/deploy/.env`; deploy.sh
then builds it, restarts it and checks it is running. Remove the line (and
`docker compose stop parlay_quoter`) to turn it off.

```
cd ~/NFLWork/unabated_ticket/deploy
docker compose logs -f parlay_quoter          # one summary line a minute
touch ~/unabated-data/parlay_quoter.stop      # pause the loop (rm to resume);
                                              # the WS feed keeps running
```

Load (measured on the VM, 2026-11-07): ~100 in-league RFQs a second at the
$50 floor. Each leg is priced at most once per `PARLAY_PRICE_CACHE_SEC`; RFQs
beyond `PARLAY_MAX_BACKLOG` waiting for a worker log as `pricer_backlog`; a
runner failure pauses calls for `PARLAY_PRICER_BACKOFF_SEC`. Leg rows are
written only for RFQs that reached the pricer.

It signs with the same read-only Kalshi key as the bets service
(`KALSHI_API_KEY_ID`, read from `bets_service/.env` and nothing else from
that file) and writes only `/data/kalshi_parlay_quoter.duckdb`.

Locally: `KALSHI_API_KEY_ID=… KALSHI_PRIVATE_KEY_PATH=… python -m kalshi_parlay_quoter.main`.

| env | default | |
|---|---|---|
| `PARLAY_QUOTER_DB_PATH` | package dir | decision log |
| `UNABATED_RUNNER_URL` | `http://127.0.0.1:8095` | leg pricer |
| `PARLAY_PRICER_TIMEOUT_SEC` | 3 | slower = `pricer_error` (detail `timeout`, no backoff) |
| `PARLAY_PRICER_WORKERS` | 8 | concurrent pricer calls |
| `PARLAY_MIN_RFQ_TARGET_COST_USD` | 50 | 0 logs every size as a row |
| `PARLAY_MAX_MARKET_FETCHES_PER_SEC` | 5 | legless WS frames need a REST market fetch |
| `PARLAY_PRICE_CACHE_SEC` | 10 | reuse a leg's price this long |
| `PARLAY_MAX_BACKLOG` | 200 | queue cap; beyond it `pricer_backlog` |
| `PARLAY_PRICER_BACKOFF_SEC` | 5 | pause after a runner failure |
| `PARLAY_STOP_PATH` | package dir `.stop` | pause file |

Tests: `python -m pytest kalshi_parlay_quoter/tests`.
