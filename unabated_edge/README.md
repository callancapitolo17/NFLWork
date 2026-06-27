# unabated_edge

Sport-agnostic, Unabated-anchored Kalshi edge-detection engine. Runs dry only — no order placement. Soccer (World Cup moneylines) is the first adapter.

---

## Architecture

```
unabated_edge/
  config.py          # constants + .env loader (BANKROLL, KELLY_FRACTION, paths …)
  feed.py            # Unabated public feed: snapshot + cursor-delta polling
  pricing.py         # american_to_prob, devig (probit, wraps kalshi_common)
  ev.py              # edge_for_yes — net of Kalshi taker fees
  sizing.py          # kelly_contracts — fractional Kelly with per-match cap
  mapping.py         # pair_events, validate (fail-closed)
  storage.py         # DuckDB writes: line_snapshots, flagged_edges, research_events
  risk.py            # staleness_ok, tipoff_ok, kill_switch_ok
  runner.py          # run_tick, main_loop, cli() entry point
  log_setup.py       # rotating file log + stderr (when tty)

  sports/
    base.py          # SportAdapter ABC
    soccer.py        # Soccer adapter (lg21, home/draw/away, WC series)
    registry.py      # ADAPTERS list + league_prefixes() / by_league()

  venues/
    kalshi.py        # list_events, best_yes_ask (via kalshi_common.auth_client)
```

**Generic core, per-sport adapters.** Everything inside `feed.py`, `pricing.py`, `ev.py`, `sizing.py`, `mapping.py`, `storage.py`, and `runner.py` is sport-agnostic. All soccer-specific knowledge lives in `sports/soccer.py`. Adding a new sport is a single adapter file + one registry line (see below).

**Data flow per tick:**
1. `feed.py` polls the Unabated changes endpoint for sharp-book line updates.
2. `runner.run_tick` calls `adapter.fair()` to devig the sharp anchor prices into true probabilities.
3. `mapping.pair_events` matches Unabated events to open Kalshi markets by canon team name.
4. `ev.edge_for_yes` computes net EV (after Kalshi taker fee) for each outcome's YES market.
5. Candidates above `MIN_EV_PCT` are Kelly-sized and written to `flagged_edges` (MARKET_DB) and the research firehose (RESEARCH_DB). All line snapshots are written to `line_snapshots` every tick for CLV tracking.

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
| `MAX_STALENESS_SEC` | `20` | Reserved for Plan 2 live-execution staleness gate — not yet enforced |
| `KICKOFF_CUTOFF_MIN` | `3` | Stop flagging this many minutes before kickoff |
| `PER_MATCH_CAP_PCT` | `0.03` | Max fraction of bankroll per match (3 %) |

---

## Running

```bash
# install deps (from repo root or unabated_edge/)
pip install -r unabated_edge/requirements.txt

# run (always dry-run — no orders are placed)
python -m unabated_edge.runner
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
from unabated_edge.sports.base import SportAdapter
from unabated_edge import pricing

class MyNewSport(SportAdapter):
    sport = "my_sport"           # unique string key
    league_prefix = "lgXX"       # Unabated league prefix (e.g. "lg21" = soccer WC)
    outcomes = ("home", "away")  # tuple of outcome names

    def fair(self, state, event_meta) -> dict | None:
        # look up sharp anchor prices from state.lines, devig, return {outcome: prob}
        ...

    def canon_team(self, name: str) -> str:
        # normalize team name (lowercase, strip aliases, etc.)
        return name.strip()

    def kalshi_series(self) -> str:
        # Kalshi series_ticker for this sport's markets
        return "KXMYSERIES"

    def map_outcome_tickers(self, kalshi_event: dict) -> dict:
        # return {outcome_name: market_ticker} from Kalshi event dict
        ...
```

2. **Register it in `sports/registry.py`:**

```python
from unabated_edge.sports.soccer import Soccer
from unabated_edge.sports.my_sport import MyNewSport

ADAPTERS = [Soccer(), MyNewSport()]
```

That is all. The engine picks up every adapter in `ADAPTERS` automatically each tick.

The `feed.py` snapshot filters by `league_prefixes()`, so only the relevant leagues are ingested. The `storage.py` tables use the `sport` column to separate rows by adapter.

---

## Troubleshooting

**`401` from the deltas endpoint**
Token expired. Recapture `unabated_at_prod` from a logged-in browser session and update `unabated_edge/.env`.

**`fair()` returns `None` / no edges flagged**
Sharp-book prices may not be published for this sport or event yet. Check `state.lines` keys and verify the `league_prefix` matches the Unabated feed (e.g. `"lg21"` for World Cup). If all three anchor source IDs (7, 6, 68) return `None` for every event, the token is likely blurred (logged-out) — re-authenticate.

**Soccer draw bt4 location / `WC_MATCH_SERIES` value**
`sports/soccer.py` has two constants marked `# REPLACE from Task 0 FINDINGS`:
- `WC_MATCH_SERIES = "KXWCMATCH"` — Kalshi series ticker for WC match markets. Verify against `GET /series` or live Kalshi event listings.
- `_draw()` — reads `bt4` on side `"1"` from the sharp-book source. Confirm by live recon against the Unabated feed during a World Cup match.

Until these are verified against live data, the soccer adapter may produce no pairings.

**Import errors / missing `kalshi_common`**
Run from the repo root so that `kalshi_common/` is on the path, or install in editable mode. The engine inserts the repo root into `sys.path` at import time.
