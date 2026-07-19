# unabated_edge

Sport-agnostic, Unabated-anchored Kalshi edge-detection engine. The taker/flagging path runs dry (no orders); the in-process market maker (`maker/`, see below) places real resting orders when `MAKER_MODE=live`. Soccer (World Cup **totals**) is the first adapter.

---

## How soccer pricing works

Soccer prices the **regulation-time total** (Over/Under goals), not the moneyline.

**Why totals only:**
- Unabated's sharp anchor books (IDs 7, 6, 68) quote `bt3` (total) directly for WC matches.
- They do NOT quote a moneyline (`bt1`) for soccer in the authenticated feed.
- We do **not** build a soccer model. The sharp anchor price IS truth; devig it and compare to Kalshi. If the anchor doesn't quote it, we don't price it.

**Pricing flow:**
1. For each Unabated event, find the first anchor book with at least one complete rung — Over (`side=0|bt3`) and Under (`side=1|bt3`) at the same line. Its main quote **plus its `alternateLines` ladder** (both sides) form the anchor ladder, e.g. `{1.5: p, 2.0: p, 2.5: p, …}`. Same-book only — never mixes books' vig.
2. Devig each rung's two-way quote (probit) → `p_over` per line.
3. Match the Kalshi series `KXWCTOTAL` — each event lists Over-ladder markets with `strike_type="greater"`, `floor_strike=X.5`.
4. For every Kalshi rung whose `floor_strike` appears in the anchor ladder, emit two candidates: **Over = buy YES**, **Under = buy NO** on that market.
5. Fail closed at every step: missing anchor side, no common line, or no matching Kalshi rung → skip, no flag (never interpolate between lines).

---

## Architecture

```
unabated_edge/
  config.py          # constants + .env loader (BANKROLL, KELLY_FRACTION, paths …)
  feed.py            # Unabated feeds: v2 per-league polling (current) + legacy snapshot/deltas
  pricing.py         # american_to_prob, devig (probit, wraps kalshi_common)
  ev.py              # edge_for_yes — net of Kalshi taker fees
  sizing.py          # kelly_contracts — fractional Kelly with per-match cap
  mapping.py         # pair_events — match Unabated↔Kalshi events by canon team pair
  storage.py         # DuckDB writes: line_snapshots, flagged_edges, research_events
  risk.py            # tipoff_ok, kill_switch_ok
  runner.py          # run_tick, main_loop, cli() entry point
  log_setup.py       # rotating file log + stderr (when tty)

  sports/
    base.py          # SportAdapter ABC + Candidate dataclass
    soccer.py        # Soccer adapter (lg21, KXWCTOTAL, totals-only)
    registry.py      # ADAPTERS list + league_prefixes() / by_league()

  venues/
    kalshi.py        # list_events, best_yes_ask, best_no_ask (via kalshi_common.auth_client)
```

**Generic core, per-sport adapters.** Everything inside `feed.py`, `pricing.py`, `ev.py`, `sizing.py`, `mapping.py`, `storage.py`, and `runner.py` is sport-agnostic. All soccer-specific knowledge lives in `sports/soccer.py`. Adding a sport = one adapter file + one registry line.

**Data flow per tick (every `V2_POLL_SEC`, default 5s):**
1. `feed.fetch_v2(league_id, league_prefix)` re-fetches Unabated's **v2 per-league odds file** (`content.unabated.com/markets/v2/league/<id>/odds.json`). Anchors are **unblurred anonymously** in this file — no token needed — and each line carries its full `alternateLines` ladder. (The legacy `changes/query` delta feed does not carry soccer at all; the legacy snapshot's anchors are blurred. Both legacy functions remain in `feed.py` for future US-league adapters.)
2. `mapping.pair_events` matches Unabated events to open Kalshi events by canonical team-pair (parsed from the Kalshi event title).
3. `adapter.price_event()` calls `_anchor_totals` to devig the bt3 over/under ladder from the first complete anchor book, then emits candidates at every Kalshi rung matching a ladder line by `floor_strike`.
4. **Kalshi microstructure capture:** for every paired pre-kickoff event, each market's full orderbook is fetched **once** per tick (`venues.kalshi.get_book`) and written to `book_snapshots` — top-of-book bid/ask/size columns plus the full depth ladder as JSON, and the market's `volume_fp`/`open_interest_fp`. Every `TRADES_POLL_SEC` (default 30s) the executed-trades tape is polled per market into `kalshi_trades` (PK `trade_id`, overlapping poll windows dedup via `INSERT OR IGNORE`; a >100-trade burst inside one window loses the excess). Both stop at kickoff, same close semantics as `line_snapshots`. This is the maker-design dataset: spread width, re-centering speed after anchor moves, and where flow trades.
5. For each `Candidate` (Over=YES, Under=NO), `ev.edge_for_yes` computes net EV after Kalshi taker fee. Asks are derived from the same fetched book (`yes_ask_from_book` / `no_ask_from_book`), so pricing costs no extra REST calls.
6. Candidates above both `MIN_EV_PCT` and `MIN_EV_DOLLARS` are Kelly-sized and written to `flagged_edges`. Line snapshots (with the feed's `modified_on`) are written every tick **until kickoff** — so the last snapshot per event is the closing line by construction. Every priced candidate goes to the research firehose with rung provenance (`book`, `alt`, `overround`) so alt-line trustworthiness stays auditable. If both sides of one rung flag at once, that's a crossed/stale book and is logged as a data error, not a double edge.

**Databases (both gitignored, pkg-local):**
- `unabated_edge_market.duckdb` — `line_snapshots`, `flagged_edges`, `book_snapshots` (Kalshi bid/ask/depth per rung per tick), `kalshi_trades` (executed-trades tape)
- `unabated_edge_research.duckdb` — `research_events` (full candidate lifecycle)

**Shared dependency:** `kalshi_common/` (repo root) provides probit devigging (`fair_value._probit_devig_n`), fee math (`ev_calc.fee_per_contract`), and the authenticated Kalshi REST client (`auth_client`).

---

## Market maker (maker/)

`unabated_edge/maker/` rests two-sided limit orders (fair ± margin) on every
rung the anchor ladder prices, in-process off the same 5s tick that already
fetched the anchor and the Kalshi book — no separate daemon, no duplicated
polling (unlike `kalshi_mlb_mm/`, which is standalone). It earns the
bid/ask spread from retail flow; the taker dry-run found no taker edge on
these markets (1c spreads sitting on the anchor), so the maker is the play.
Spec: `docs/superpowers/specs/2026-07-10-wc-totals-maker-design.md`.

```
unabated_edge/maker/
  engine.py   # per-match quote decisions: margins, never-cross, alt gating, pulls, cap enforcement
  ledger.py   # pure math: fills -> pnl(g) over goal-grid g=0..10 -> worst case
  gateway.py  # QuoteGateway ABC; LiveGateway (REST orders) / ShadowGateway (record only)
  state.py    # resting-order book + Kalshi fills/positions/settlements reconciliation
  store.py    # maker_quotes / maker_fills / ledger_snapshots / maker_pnl writes
```

**Quoting.** For each rung, margin = `max(maker_fee_per_contract + PICKOFF_BUFFER_CENTS, ROI_MARGIN × price)`
(alt rungs multiply margin by `ALT_MARGIN_MULT`). The bid is capped one tick
below the opposing best ask (never cross); if that cap would sit tighter
than `fair − MAX_MARGIN_CENTS`, the rung is skipped — the crowd is already
tighter than we're willing to quote. Alt (non-main) rungs only quote inside
an `[ALT_OVERROUND_MIN, ALT_OVERROUND_MAX]` vig band and get smaller size
(`ALT_SIZE_MULT`).

**Touch-join pricing.** When the crowd's best same-side bid offers better
entry than our resting fair−margin quote, we join it instead (the opposing
ask still caps the join one tick below it — never cross). Entry condition:
net edge (fair − join_price − maker_fee) ∈ [TOUCH_JOIN_MIN_EDGE_CENTS, MAX_MARGIN_CENTS];
alt rungs use TOUCH_JOIN_ALT_MIN_EDGE_CENTS (default 1.5c). The calculation
self-excludes our own resting order from the touch. Hysteresis: hold join order
while edge ≥ TOUCH_JOIN_EXIT_EDGE_CENTS (0.25c), exit if edge < 0.25c or >5c
(suspicious distortion). Config dials: `TOUCH_JOIN_MIN_EDGE_CENTS` (1.0c main),
`TOUCH_JOIN_ALT_MIN_EDGE_CENTS` (1.5c), `TOUCH_JOIN_EXIT_EDGE_CENTS` (0.25c).
Joins are attributed `touch_join` in `maker_quotes` vs `quote` for standard-margin
orders. Deferred protections: spec `docs/superpowers/specs/2026-07-18-wc-touchjoin-pricing-design.md`,
GitHub issues #1 (mid-jump circuit breaker) and #2 (markout/toxicity brake).

**The goal-grid ledger** (`ledger.py`) treats every fill on a match as
settling on one integer — the regulation-time total goals `g` — and computes
exact worst-case P&L by evaluating `pnl(g)` over `g = 0..10` and taking the
min, instead of per-market caps. Every candidate quote is sized (or
dropped) so a full hypothetical fill keeps that worst case inside the caps
below. **Deviation from spec:** the ledger's exposure snapshot
(`state.exposure_fills`) folds in *resting* (not-yet-filled) quotes as if
already filled, so simultaneous fills across rungs can never breach a cap —
the reported "worst case" is a hypothetical ceiling, not realized exposure.

**Cap stack** (percent of `BANKROLL`; all four are separate env vars):

| cap | param | default | at $1,000 |
|---|---|---|---|
| per resting order | `MAX_QUOTE_PCT` | 0.30 | $300 (ledger usually binds first) |
| per match worst case | `MATCH_CAP_PCT` | 0.40 | $400 |
| global Σ worst cases (all matches) | `GLOBAL_CAP_PCT` | 0.75 | $750 |
| daily realized-loss halt | `DAILY_LOSS_HALT_PCT` | 0.40 | $400 |

A second same-day match's budget is squeezed by whatever the first match's
worst case has already committed against the global cap.

**Modes — `MAKER_MODE=off|shadow|live`** (`config.py`, `.env`-overridable).
**Deviation from spec:** the default is `off` (spec drafted `live` as the
default; shipped safer). `shadow` builds a `ShadowGateway` that logs every
quote decision to `maker_quotes` but places nothing — for testing future
leagues, not this v1 path. `live` builds a `LiveGateway` that places/cancels
real orders via `kalshi_common.auth_client`, and refuses to start unless
`MAKER_LIVE_ACK=1` is also set (`gateway.make_gateway` raises `SystemExit`
otherwise) — a dead-man switch so `MAKER_MODE=live` can never go live by typo.

**Pull triggers** (checked per match, per tick, each cancels every resting
quote on that match unless noted as global):

| trigger | scope | action |
|---|---|---|
| recomputed quote price differs from the resting price (fair moved ≥ 1c, or the opposing ask moved the never-cross cap) | that rung | cancel/replace that rung |
| feed watchdog: no successful tick for a sport in `MAX_STALENESS_SEC` (20s) | all matches | pull everything (`watchdog()`, run every main-loop iteration) |
| within `QUOTE_PULL_MIN` (3 min) of kickoff | that match | pull; inventory rides to settlement |
| `COOLOFF_MIN` (10 min) after a fill burst | that match | pull + hold off requoting |
| fill burst: > `FILL_BURST_N` (3) fills in 60s | that match | pull + start cooloff |
| crossed/impossible book (`yes_ask + no_ask < 1 − 2·fee`) | that match | pull (data error) |
| unpaired / kickoff passed / market closed | that match | pull (`sweep()`) |
| Kalshi position mismatch vs local fills | all matches | pull everything (live mode only — reconciliation polling runs only when `MAKER_MODE=live`) |
| daily loss halt (§cap stack) | all matches | pull everything |
| hard $-stop: realized + mark-to-anchor unrealized P&L ≤ `-HARD_STOP_DOLLARS` ($50) | all matches | pull everything (`hard_stop`) — checked every tick, right after the daily halt |
| anchor ladder disappears for a match | that match | pull that match (`anchor_gone`) |
| anchor frozen: even the freshest rung's `modifiedOn` is older than `ANCHOR_STALE_SEC` (180s), or no rung has a parseable timestamp (fail-safe) | that match | pull that match (`anchor_stale`) |
| Kalshi position on this event's ticker still carries a startup baseline (restart-with-inventory not yet settled) | that match | pull/skip that match (`baseline_blocked`) — no fresh quotes rest on top of unreconciled inventory |
| no opposing crowd bid on a side (`yes_ask`/`no_ask` is `None`) | that side | skip, don't quote (`no_crowd`) — never quote unconstrained into an empty book |
| graceful shutdown (SIGINT/SIGTERM or kill file exit) | all matches | pull ALL resting quotes (`shutdown`/`kill_switch`) — global |
| `.kill` file present | all matches | pull everything — **deviation from spec's "watchdog" framing**: the check is at the top of `run_tick` itself (same 5s tick cadence), not inside the `watchdog()` staleness method |

**Data model** — sibling DB `unabated_edge_maker.duckdb` (`maker/store.py::init`),
kept separate from the capture DBs so maker state can reset independently:

| table | row = |
|---|---|
| `maker_quotes` | every quote decision (rest/replace/cancel/skip), with price, size, fair, margin, reason |
| `maker_fills` | every real fill (PK `trade_id`), price, contracts, fee, ledger worst-case after |
| `ledger_snapshots` | per match per tick: worst case, full `pnl_grid` JSON, quotes live |
| `maker_pnl` | per market (PK `market_ticker`) settled P&L |

**Reconciliation.** A fresh `MakerState` knows nothing about a previous run
(crash/restart) or a POST that landed but errored before returning an
`order_id` (`auth_client` never retries POSTs, so the next tick gets a new
`client_order_id`). In `MAKER_MODE=live`, startup adopts any existing
in-series Kalshi positions as a baseline (`state.position_baseline`) —
so the mismatch tripwire doesn't false-trip on inventory carried over from
before — and cancels any in-series resting orders left over from before
(`state.startup_sync`). Every 60s recon cycle repeats the order sweep
(`state.sweep_orphan_orders`), bounding the orphaned-order window from a
landed-but-errored POST to ~60s. Position mismatches beyond
baseline + local fills still pull all quotes (`position_mismatch`, see the
pull-triggers table below). Operational caveat: a mid-match restart rebuilds positions and orders but not the per-match goal-grid ledger — pre-restart fills don't count toward the match cap until settlement, so don't restart mid-match while holding inventory (flatten or wait for settlement first).

**Runbook:**

```bash
MAKER_MODE=live MAKER_LIVE_ACK=1 python3 -m unabated_edge.runner
```

`MAKER_MODE=shadow` (no ack needed) exercises the same code path without
placing orders. Creating the existing `.kill` file does **not** stop the
process loop — it makes `run_tick` skip pricing and, every tick (5s), call
`maker.pull_all(reason="kill_switch")` to cancel every resting quote, for as
long as the file exists. Graceful shutdown — `SIGINT`/`SIGTERM` (handlers
clear `_running`) or the `.kill` file — pulls every resting maker quote
before the process exits: `main_loop` calls `maker.pull_all(reason="shutdown")`
right after its `while _running.is_set():` loop ends. Only killing the
process directly (`kill -9`) or a crash skips this cleanup and leaves any
resting quotes live on Kalshi's book until the next `poll_positions`
reconciliation on restart.

**Operational runbook (live) — tuition run:**

- **$300 account note:** `BANKROLL` stays `1000` for the cap-stack math (§cap
  stack percentages are computed off it), but the leashes that actually bind
  at this account size are `MAKER_MAX_CONTRACTS=2` (every order is ~$1) and
  `HARD_STOP_DOLLARS=50`; the %-based caps (`MAX_QUOTE_PCT`/`MATCH_CAP_PCT`/
  `GLOBAL_CAP_PCT`/`DAILY_LOSS_HALT_PCT`) are outer ceilings that never bind at
  this size. **Never raise `MAKER_MAX_CONTRACTS` without first checking it
  against the real account balance.**
- **Single writer:** stop the capture runner before starting the maker (the
  maker captures too — running two `unabated_edge.runner` processes means two
  writers on the same DuckDB files). A `.runner.lock` PID file at
  `unabated_edge/.runner.lock` enforces one runner: a second `main_loop` call
  refuses to start (`SystemExit`) while the PID in the lock file is alive, and
  automatically reclaims a stale lock left by a crashed/`kill -9`'d process.
- **Launch inline, never in `.env`:** start live mode with
  `MAKER_MODE=live MAKER_LIVE_ACK=1 python3 -m unabated_edge.runner` from the
  worktree cwd — as an inline env var, not baked into `unabated_edge/.env` —
  so a bare restart (no env vars) always returns to `MAKER_MODE=off`.
- **Safe stop:** `touch unabated_edge/.kill`, wait one 5s tick so the kill
  switch pulls all resting quotes (`kill_switch`), THEN send `Ctrl-C` /
  `SIGTERM`. Never `kill -9` while quotes are live — it skips the shutdown
  pull and can leave resting orders live on Kalshi until the next
  `poll_positions` reconciliation.
- **No mid-match restart while holding inventory** (see the Reconciliation
  caveat above) — a restart rebuilds positions/orders but not the per-match
  goal-grid ledger, so pre-restart fills don't count toward the match cap
  until settlement.

---

## Authentication

### Unabated token (NOT needed for soccer)

The v2 per-league odds file the engine polls is **anonymous — anchors are unblurred without any login**, so the soccer path runs with no Unabated credentials at all.

The token below is only needed for the **legacy** changes endpoint (`/api/markets/changes/query`, used by future US-league adapters): a premium-account JWT passed as a cookie named `unabated_at_prod`, valid ~30 days.

**Capture steps:**
1. Log in to [unabated.com](https://unabated.com) in a browser.
2. Open DevTools → Network tab. Trigger a line update or navigate to any odds page.
3. Find a `changes/query` request. Right-click → Copy → Copy as cURL.
4. Extract the value of the `unabated_at_prod` cookie from the cURL command.
5. Store it in `unabated_edge/.env`:

```
UNABATED_AT_PROD=<paste token here>
```

**Refresh:** the token expires after ~30 days. When it does, the deltas endpoint returns HTTP 401. Recapture by repeating the steps above and updating `.env`.

### Kalshi credentials

Stored in `.env` alongside the Unabated token:

```
KALSHI_API_KEY_ID=<your key id>
KALSHI_PRIVATE_KEY_PATH=<absolute path to PEM file>
# optional — defaults to production URL
KALSHI_BASE_URL=https://api.elections.kalshi.com/trade-api/v2
```

---

## Configuration

All tuneable constants live in `config.py` and can be overridden via `.env` or environment variables:

| Variable | Default | Description |
|---|---|---|
| `BANKROLL` | `1000.0` | Total bankroll in USD |
| `KELLY_FRACTION` | `0.25` | Fractional Kelly multiplier |
| `MIN_EV_PCT` | `0.03` | Minimum net EV % to flag (3 %) |
| `MIN_EV_DOLLARS` | `0.02` | Minimum absolute per-contract EV (gates with `MIN_EV_PCT` so cheap longshots can't clear on % alone) |
| `UNABATED_V2_POLL_SEC` | `5` | Seconds between re-fetches of each league's v2 odds file (also the tick cadence) |
| `KALSHI_TRADES_POLL_SEC` | `30` | Seconds between executed-trades tape polls (rides every Nth main tick) |
| `MAX_STALENESS_SEC` | `20` | Taker path: reserved, not enforced. **Maker path: enforced** — `maker.watchdog()` pulls all quotes if a sport's feed hasn't ticked successfully within this many seconds |
| `KICKOFF_CUTOFF_MIN` | `3` | Stop flagging this many minutes before kickoff |
| `PER_MATCH_CAP_PCT` | `0.03` | Max fraction of bankroll per match (3 %) — taker sizing only |
| `MAKER_MODE` | `off` | `off` \| `shadow` \| `live` — see [Market maker](#market-maker-maker) |
| `MAKER_LIVE_ACK` | unset | Must be `"1"` for `MAKER_MODE=live` to start (dead-man switch) |
| `ROI_MARGIN` | `0.03` | Margin floor as a fraction of price (3 %), before the fee+buffer floor |
| `PICKOFF_BUFFER_CENTS` | `1` | Cents added above `maker_fee_per_contract` for the margin floor |
| `MAX_MARGIN_CENTS` | `5` | Skip a rung if never-cross would tighten our bid beyond `fair − this` |
| `ALT_MARGIN_MULT` | `1.5` | Margin multiplier on alt (non-main) rungs |
| `ALT_SIZE_MULT` | `0.5` | Size multiplier on alt rungs |
| `ALT_OVERROUND_MIN` | `1.01` | Alt rung quotes only if the rung's devig overround is ≥ this |
| `ALT_OVERROUND_MAX` | `1.15` | Alt rung quotes only if the rung's devig overround is ≤ this |
| `QUOTE_PULL_MIN` | `3` | Minutes before kickoff to pull a match's quotes (inventory rides) |
| `MAX_QUOTE_PCT` | `0.30` | Max fraction of bankroll per resting order (ledger cap usually binds first) |
| `MATCH_CAP_PCT` | `0.40` | Max ledger worst-case per match, as a fraction of bankroll |
| `GLOBAL_CAP_PCT` | `0.75` | Max Σ ledger worst-case across all matches, as a fraction of bankroll |
| `DAILY_LOSS_HALT_PCT` | `0.40` | Realized settled loss (fraction of bankroll) that halts quoting for the day |
| `FILL_BURST_N` | `3` | More than this many fills on one match in 60s trips the fill-burst tripwire |
| `COOLOFF_MIN` | `10` | Minutes a match stays pulled after a fill-burst trip |
| `MAKER_MAX_CONTRACTS` | `2` | Hard per-quote contract ceiling (tuition-run leash — binds before the %-based ledger caps at this account size) |
| `HARD_STOP_DOLLARS` | `50` | Cumulative realized + mark-to-anchor unrealized loss (dollars) that halts the maker for the day, separate from `DAILY_LOSS_HALT_PCT` |
| `ANCHOR_STALE_SEC` | `180` | Pull a match if even its freshest anchor rung's `modifiedOn` is older than this — frozen-feed guard |

---

## Running

```bash
# install deps (from repo root or unabated_edge/)
pip install -r unabated_edge/requirements.txt

# run (taker path is dry-run; maker places orders only with MAKER_MODE=live + MAKER_LIVE_ACK=1)
python3 -m unabated_edge.runner
```

The bot re-fetches each adapter's v2 league odds file every `V2_POLL_SEC` (default 5s), refreshes Kalshi event listings every 30 seconds, and logs flagged edges to `unabated_edge/bot.log` and to the two DuckDB files. A heartbeat line (~every 60s) reports event/line/market counts **and candidates priced** — `candidates_recent=0` while `lines>0` means the pricing chain is broken (re-blurred anchors, feed shape drift, or no Kalshi pairings), which a raw line count alone would hide.

**Kill switch:** create `unabated_edge/.kill` to pause pricing/flagging (and, if the maker is enabled, pull its resting quotes) every tick for as long as the file exists — this does **not** exit the process. `SIGINT`/`SIGTERM` exit the loop and, if the maker is enabled, pull all resting quotes on the way out — see [Market maker](#market-maker-maker) for details.

---

## Querying results

```python
import duckdb

# Edges flagged this session
con = duckdb.connect("unabated_edge/unabated_edge_market.duckdb", read_only=True)
print(con.execute("SELECT * FROM flagged_edges ORDER BY ts DESC LIMIT 20").df())

# CLV backbone — compare fair_prob at flag time to closing line
print(con.execute("""
    SELECT sport, event_id, market_source_id, bet_type, side,
           MIN(price) as open_price, MAX(ts) as last_ts
    FROM line_snapshots
    GROUP BY 1,2,3,4,5
    ORDER BY last_ts DESC
""").df())

# Research firehose (every candidate priced, not just flagged ones)
res = duckdb.connect("unabated_edge/unabated_edge_research.duckdb", read_only=True)
print(res.execute("SELECT * FROM research_events ORDER BY ts DESC LIMIT 50").df())

# Kalshi microstructure (maker design data): spread + depth per rung over time
print(con.execute("""
    SELECT market_ticker, floor_strike,
           avg(yes_ask - yes_bid) AS avg_spread,
           avg(yes_bid_qty) AS avg_top_depth,
           max(volume) AS volume
    FROM book_snapshots GROUP BY 1,2 ORDER BY 1
""").df())

# Where flow trades + taker direction (adverse-selection input)
print(con.execute("""
    SELECT market_ticker, taker_side, count(*) AS n, sum(count) AS contracts
    FROM kalshi_trades GROUP BY 1,2 ORDER BY contracts DESC
""").df())
```

Both DBs have a `sport` column so multi-sport rows never collide.

---

## How to add a sport

1. **Write `sports/<sport>.py`** implementing `SportAdapter`:

```python
from unabated_edge.sports.base import SportAdapter, Candidate
from unabated_edge import pricing, config
from unabated_edge.feed import line_american_price

class MyNewSport(SportAdapter):
    sport = "my_sport"           # unique string key
    league_prefix = "lgXX"       # Unabated league prefix

    def canon_team(self, name: str) -> str:
        return name.strip()

    def kalshi_series(self) -> str:
        return "KXMYSERIES"

    def event_teams(self, kalshi_event: dict) -> frozenset:
        # parse two team names from kalshi_event["title"] (or markets)
        ...

    def price_event(self, state, event_meta, kalshi_event) -> list[Candidate]:
        # anchor from state.lines, devig, find Kalshi market, return Candidates
        # return [] to fail closed
        ...
```

2. **Register it in `sports/registry.py`:**

```python
from unabated_edge.sports.soccer import Soccer
from unabated_edge.sports.my_sport import MyNewSport

ADAPTERS = [Soccer(), MyNewSport()]
```

That is all. The engine picks up every adapter in `ADAPTERS` automatically each tick.

**Rule:** only price markets that the sharp anchor quotes directly. Do not build a model to derive prices for markets the anchor doesn't quote.

---

## Troubleshooting

**`401` from the deltas endpoint**
Token expired. Recapture `unabated_at_prod` from a logged-in browser session and update `unabated_edge/.env`.

**No edges flagged / `_anchor_total` returns None**
Anchor books may not be quoting bt3 for these events yet, or both Over and Under must appear at the same line value. Check `state.lines` keys for `|bt3|` entries. If every event returns None, the token is likely blurred (logged-out) — re-authenticate. The heartbeat log line shows `events`/`lines`/`kalshi_events` counts so you can tell "broken" (zeros) from "healthy, no edges today".

**No Kalshi market matched (fail closed)**
The anchor quoted a line (e.g. 2.5) but `KXWCTOTAL` has no market with `floor_strike=2.5`. This is expected when Kalshi's rung ladder doesn't cover that line — nothing to trade against.

**Import errors / missing `kalshi_common`**
Run from the repo root so that `kalshi_common/` is on the path, or install in editable mode. The engine inserts the repo root into `sys.path` at import time.
