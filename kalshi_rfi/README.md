# kalshi_rfi — one-sided Kalshi YRFI maker

Standalone daemon that rests **YES (YRFI) bids** on Kalshi's `KXMLBRFI`
("run in the first inning") markets at book-consensus fair minus a margin.
One-sided by design (user decision 2026-08-25): retail flow buys NRFI (NO),
which matches against resting YES bids — we take the other side at a
discount. We never quote NO and never cross the spread.

## How it prices

A `KXMLBRFI` market is exactly a 1st-inning total at 0.5 (issue #87). Each
game's fair P(YRFI) comes from a **live** fetch at the 5 books with working
period-aware I1 hooks — DraftKings, FanDuel, BetMGM, Novig, Caesars —
through the shared `SGPService.price_on_demand` on a single
`CanonicalLeg(..., "I1")`.
A lone leg takes the **opt-in** n==1 fast path added for this bot (exact
2-cell probit devig of the structure's own odds via `devig_partition`, so
crossed/degenerate pairs fail the [1.0, 1.25] overround gate; zero SGP price
calls). The fast path is enabled per `SGPService` instance
(`single_leg_structure_fair=True`) and this bot pairs it with
`structure_ttl_sec=0.0` so **every fair refresh re-fetches the book's odds
from the wire** — never the structure TTL cache. Both knobs default off/420s,
so the MM/taker SGP pipeline is completely unaffected. DK never contributes:
its SGP structure carries selection ids but no odds, so the n==1 fast path
cannot run there, and its fallback price route is dead too — both documented
`calculateBets` paths returned Akamai 403 on 2026-08-25 (issue #39
recurring; DK's READ endpoints stay green, so only pricing is blocked).
In practice **4 books price**: FD/MGM/Novig/CZR (live-verified 2026-08-25
over 6 games — Caesars priced 6/6, four-book σ_z 0.017–0.030).

Consensus gate mirrors the MM's issue #20 semantics: ≥ `RFI_MIN_BOOKS` (2)
books and sample stddev of probit-transformed fairs ≤ `RFI_SIGMA_Z_MAX`
(0.07), median fair, no outlier removal.

## Quote lifecycle (per game, per 30s cycle)

1. Board refresh: one `GET /markets?series_ticker=KXMLBRFI&status=open`
   gives every game's touch. First pitch is parsed from the ticker suffix
   (US/Eastern → UTC) — **never** `close_time` (it is first pitch + 72h).
2. Fair refresh on TTL: 300s normally, 60s inside the last hour
   (`RFI_FAIR_REFRESH_FAR_SEC` / `_NEAR_SEC` / `RFI_NEAR_WINDOW_MIN`); each
   fetch is wall-bounded by `RFI_FAIR_FETCH_WALL_SEC` (15s) — a hung book
   counts as not-priced rather than freezing the loop.
3. Jump guard between refreshes (`engine.kalshi_book_moved`): bid and ask
   compared **per side** against the touch cached at fair-fetch time —
   ≥ `RFI_KALSHI_JUMP_CENTS` (3¢) on either side cancels and refetches. A
   best bid equal to our own resting price is ignored (our order pins the
   bid side; a mid-based guard both self-triggered on placement and needed
   a 6¢ ask move). Residual blind spot: depth evaporating *below* our bid
   is invisible at top-of-book.
4. Desired bid = `floor(fair·100) − RFI_MARGIN_CENTS`, clamped below
   the ask (always maker), refused unless post-maker-fee edge ≥
   `RFI_MIN_EDGE_CENTS` (2¢).
5. Size: worst case of a YES bid = its cost, capped at
   `RFI_PER_GAME_CAP_USD` per game and `RFI_DAILY_CAP_USD` per ET trading
   day (fills + all resting orders; settled fills still count against the
   day). Three layers, because one fill poll is not enough:
   - **every** placement re-polls fills and re-sizes first — fresh quotes
     included. (This is where the cap broke on 2026-08-27: an order that
     filled *completely* is dropped from state, so the next cycle took the
     fresh-quote path, which re-polled only when replacing an order. One
     game took 5 fills / 45 contracts (~$24) in ~60s against a $5 cap.)
     The extra poll is throttled to one per 2s across the slate — a
     seconds-old poll still closes the bug; the stale one was a full cycle
     old.
   - a **hard backstop** immediately before the wire (`main._fit_to_caps`)
     refuses or shrinks any order that would push the game's — or the
     day's — worst case past its cap, counting unsettled fills *and* every
     order resting on that game. It re-reads state and does not trust the
     count the engine sized at the top of the cycle.
   - a **post-fill cooldown** (`RFI_POST_FILL_COOLDOWN_SEC`, 60s) stands
     the game down after each fill. At a 1¢ margin the jump guard
     cancels/refetches every ~30s, so without it a filling game re-quotes
     into its own fill faster than the caps can see it.
   A same-price order still downsizes when headroom shrank.
6. Cancel on: consensus decline, edge gone, caps, market gone/suspended,
   and unconditionally at first pitch − `RFI_PULL_BEFORE_START_SEC`
   (600s — late scratches land T−10min to T−2min and books lag them; that
   window's fills are adverse selection, not volume).
   Shutdown cancels everything (no unmanaged GTC orders), and a network
   error during shutdown skips to the next ticker rather than stranding
   the rest.
7. Doubleheaders are fail-closed excluded (book event-matchers key on team
   names and can pick the wrong game of a same-day pair).

## Order tracking and reconciliation

Resting orders are keyed by **order id**, not ticker. A ticker-keyed map
dropped the first order id the moment a second order landed on the same
game, orphaning a live order from every cancel, cap and sweep by
construction (observed 2026-08-27). The bot still intends exactly one order
per game; more than one means state drifted, and `_cancel` pulls all of
them.

Because local tracking demonstrably drifts, **every cycle** diffs Kalshi's
`/portfolio/orders?status=resting` against local state
(`state.reconcile_resting_orders`, replacing the old startup-plus-300s
orphan sweep):

- an in-series order resting on Kalshi that we don't track is cancelled
  (the bot owns KXMLBRFI while it runs — see the warning below);
- a locally-tracked order Kalshi isn't resting is dropped (a placement
  younger than 5s is left alone: the listing can lag the POST);
- matched orders take Kalshi's remaining count and exchange shard, so caps
  size against reality and later cancels route correctly.

A failed listing fetch holds all local state — an empty listing must never
read as "everything is gone".

### Exchange sharding

Kalshi split its exchange into shards around 2026-08-24 and MLB lives on
`exchange_index` **3** (NFL/NBA are still 0). The index is read from each
market payload and threaded through placement (a body field — it also skips
the auto-routing lookup) and cancellation (a query param). It is never
hardcoded: markets can move shards, and other series are elsewhere.

This matters because `DELETE /portfolio/events/orders/{id}` carries no
ticker in its path — unrouted, it hits shard 0 and returns 404 for a
baseball order. Combined with the "explicit 404 = already gone" rule
(added for the MM phantom-open-quotes bug), the bot believed every cancel
succeeded and left 13 real orders resting while shutdown reported clean.
A 404 is now **verified** against the resting listing before it is
believed: still resting → the cancel failed and local state is held;
verifiably absent → treated as cancelled; listing unavailable → fails
closed and retries next cycle.

Fills/settlements follow the unabated_edge maker pattern: Kalshi is
the source of truth (`/portfolio/fills` each cycle with a 24h first-poll
backfill, `/portfolio/settlements` each 600s). Restart
hydration rebuilds filled exposure and the prior run's order→ticker map
from the DB, so fills that landed while the daemon was down still attribute.
A transient network failure anywhere in a cycle holds all state and retries
next cycle — it never crashes the daemon.

**The bot owns the KXMLBRFI series while running**: reconciliation
cancels ANY of the account's resting orders in the series, including
manually placed ones. Don't hand-trade KXMLBRFI on the same account while
the bot is live.

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
  normal (see above); check FD/MGM/Novig/CZR auth/board if others drop.
  Caesars needs a freshly minted AWS-WAF token; a mint failure shows as a
  clean per-cycle decline, and the other three still reach quorum.
- No games quoted — check the horizon (`RFI_QUOTE_HORIZON_HOURS`, 12h) and
  that markets are `active`; the bot skips games > 1h past start.
- `RFI_MODE=live requires RFI_LIVE_ACK=1` — the dead-man switch, on purpose.
