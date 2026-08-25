# kalshi_rfi — one-sided Kalshi YRFI maker

Standalone daemon that rests **YES (YRFI) bids** on Kalshi's `KXMLBRFI`
("run in the first inning") markets at book-consensus fair minus a margin.
One-sided by design (user decision 2026-08-25): retail flow buys NRFI (NO),
which matches against resting YES bids — we take the other side at a
discount. We never quote NO and never cross the spread.

## How it prices

A `KXMLBRFI` market is exactly a 1st-inning total at 0.5 (issue #87). Each
game's fair P(YRFI) comes from a **live** fetch at the 4 books with working
period-aware I1 hooks — DraftKings, FanDuel, BetMGM, Novig — through the
shared `SGPService.price_on_demand` on a single `CanonicalLeg(..., "I1")`.
A lone leg takes the **opt-in** n==1 fast path added for this bot (exact
2-cell probit devig of the structure's own odds via `devig_partition`, so
crossed/degenerate pairs fail the [1.0, 1.25] overround gate; zero SGP price
calls). The fast path is enabled per `SGPService` instance
(`single_leg_structure_fair=True`) and this bot pairs it with
`structure_ttl_sec=0.0` so **every fair refresh re-fetches the book's odds
from the wire** — never the structure TTL cache. Both knobs default off/420s,
so the MM/taker SGP pipeline is completely unaffected. DK lacks structure
odds and its SGP endpoint refuses 1-selection sets, so in practice ~3 books
price (live-verified 2026-08-25: FD/MGM/Novig priced, σ_z 0.052).

Consensus gate mirrors the MM's issue #20 semantics: ≥ `RFI_MIN_BOOKS` (2)
books and sample stddev of probit-transformed fairs ≤ `RFI_SIGMA_Z_MAX`
(0.07), median fair, no outlier removal.

## Quote lifecycle (per game, per 30s cycle)

1. Board refresh: one `GET /markets?series_ticker=KXMLBRFI&status=open`
   gives every game's touch. First pitch is parsed from the ticker suffix
   (US/Eastern → UTC) — **never** `close_time` (it is first pitch + 72h).
2. Fair refresh on TTL: 300s normally, 120s inside the last hour
   (`RFI_FAIR_REFRESH_FAR_SEC` / `_NEAR_SEC` / `RFI_NEAR_WINDOW_MIN`).
3. Jump guard between refreshes: if Kalshi's mid moved ≥
   `RFI_KALSHI_JUMP_CENTS` (3¢) since the fair was fetched, cancel and
   refetch (single-market version of the MM's constituent_jump).
4. Desired bid = `floor(fair·100) − RFI_MARGIN_CENTS` (3¢), clamped below
   the ask (always maker), refused unless post-maker-fee edge ≥
   `RFI_MIN_EDGE_CENTS` (2¢).
5. Size: worst case of a YES bid = its cost, capped at
   `RFI_PER_GAME_CAP_USD` ($10) per game and `RFI_DAILY_CAP_USD` ($100)
   per ET trading day (fills + all resting orders; settled fills still
   count against the day).
6. Cancel on: consensus decline, edge gone, caps, market gone/suspended,
   and unconditionally at first pitch − `RFI_PULL_BEFORE_START_SEC` (120s).
   Shutdown cancels everything (no unmanaged GTC orders).
7. Doubleheaders are fail-closed excluded (book event-matchers key on team
   names and can pick the wrong game of a same-day pair).

Fills/settlements/orphans follow the unabated_edge maker pattern: Kalshi is
the source of truth (`/portfolio/fills` each cycle, `/portfolio/settlements`
each 600s, startup + in-series orphan-order cancel sweep).

## Running

```bash
./kalshi_rfi/run.sh                          # shadow (default): no orders, full dataset
RFI_MODE=live RFI_LIVE_ACK=1 ./kalshi_rfi/run.sh   # live (dead-man switch required)
```

Credentials: `KALSHI_API_KEY_ID` / `KALSHI_PRIVATE_KEY_PATH` via env or
`kalshi_rfi/.env` (same names as the other bots). `RFI_MODE=off` refuses to
start. Stop cleanly by creating `kalshi_rfi/.kill` or Ctrl-C — both cancel
all resting orders on the way out.

## Data (kalshi_rfi/kalshi_rfi.duckdb — this bot is the only writer)

- `snapshots` — one row per cycle × in-horizon game: per-book fairs (JSON),
  σ_z, consensus, Kalshi touch, decision + reason, our quote. **This is the
  research dataset**; shadow mode fills it identically to live.
- `orders` — every place/cancel with reason (shadow included).
- `fills` — trade_id-deduped executions; `settlements` — realized P&L.

Useful first queries: fills vs. snapshots consensus at fill time (adverse
selection), realized YRFI rate vs. consensus (calibration), reason counts
(why we weren't quoting).

## Troubleshooting

- `consensus DECLINED books=1` — only one book priced. DK declining is
  normal (see above); check FD/MGM/Novig auth/board if others drop.
- No games quoted — check the horizon (`RFI_QUOTE_HORIZON_HOURS`, 12h) and
  that markets are `active`; the bot skips games > 1h past start.
- `RFI_MODE=live requires RFI_LIVE_ACK=1` — the dead-man switch, on purpose.
