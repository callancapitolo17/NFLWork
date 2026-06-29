# unabated_edge

Sport-agnostic, Unabated-anchored Kalshi edge-detection engine. Runs dry only — no order placement. Soccer (World Cup **totals**) is the first adapter.

---

## How soccer pricing works

Soccer prices the **regulation-time total** (Over/Under goals), not the moneyline.

**Why totals only:**
- Unabated's sharp anchor books (IDs 7, 6, 68) quote `bt3` (total) directly for WC matches.
- They do NOT quote a moneyline (`bt1`) for soccer in the authenticated feed.
- We do **not** build a soccer model. The sharp anchor price IS truth; devig it and compare to Kalshi. If the anchor doesn't quote it, we don't price it.

**Pricing flow:**
1. For each Unabated event, find the first anchor book that has BOTH Over (`side=0|bt3`) and Under (`side=1|bt3`) at the same line value.
2. Devig the two-way anchor quote (probit) → `p_over`.
3. Match the Kalshi series `KXWCTOTAL` — each event lists Over-ladder markets with `strike_type="greater"`, `floor_strike=X.5`.
4. Find the Kalshi market whose `floor_strike` equals the anchor line.
5. Emit two candidates: **Over = buy YES** on that market, **Under = buy NO** on that market.
6. Fail closed at every step: missing anchor side, mismatched line values, or no matching Kalshi rung → return `[]`, no flag.

---

## Architecture

```
unabated_edge/
  config.py          # constants + .env loader (BANKROLL, KELLY_FRACTION, paths …)
  feed.py            # Unabated public feed: snapshot + cursor-delta polling
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

**Data flow per tick:**
1. `feed.py` polls the Unabated changes endpoint for sharp-book line updates.
2. `mapping.pair_events` matches Unabated events to open Kalshi events by canonical team-pair (parsed from the Kalshi event title).
3. `adapter.price_event()` calls `_anchor_total` to devig the bt3 over/under from the first complete anchor book, then finds the matching Kalshi rung by `floor_strike`.
4. For each `Candidate` (Over=YES, Under=NO), `ev.edge_for_yes` computes net EV after Kalshi taker fee. `ask_fn(ticker, side)` fetches `best_yes_ask` for YES or `best_no_ask` for NO.
5. Candidates above both `MIN_EV_PCT` and `MIN_EV_DOLLARS` are Kelly-sized and written to `flagged_edges`. All line snapshots are written every tick for CLV tracking.

**Databases (both gitignored, pkg-local):**
- `unabated_edge_market.duckdb` — `line_snapshots`, `flagged_edges`
- `unabated_edge_research.duckdb` — `research_events` (full candidate lifecycle)

**Shared dependency:** `kalshi_common/` (repo root) provides probit devigging (`fair_value._probit_devig_n`), fee math (`ev_calc.fee_per_contract`), and the authenticated Kalshi REST client (`auth_client`).

---

## Authentication

### Unabated token

The Unabated changes endpoint (`/api/markets/changes/query`) requires a premium-account JWT passed as a cookie named `unabated_at_prod`. The token is valid for approximately 30 days.

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
| `MAX_STALENESS_SEC` | `20` | Reserved for Plan 2 live-execution staleness gate — not yet enforced |
| `KICKOFF_CUTOFF_MIN` | `3` | Stop flagging this many minutes before kickoff |
| `PER_MATCH_CAP_PCT` | `0.03` | Max fraction of bankroll per match (3 %) |

---

## Running

```bash
# install deps (from repo root or unabated_edge/)
pip install -r unabated_edge/requirements.txt

# run (always dry-run — no orders are placed)
python3 -m unabated_edge.runner
```

The bot polls the Unabated feed every ~2 seconds, refreshes Kalshi event listings every 30 seconds, and logs flagged edges to `unabated_edge/bot.log` and to the two DuckDB files.

**Kill switch:** create `unabated_edge/.kill` to stop the loop gracefully. `SIGINT`/`SIGTERM` also stop it.

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
