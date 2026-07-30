# Kalshi MLB Bots Monitor

A read-only web dashboard that monitors **both** Kalshi MLB bots on one screen:

- **Maker** (`kalshi_mlb_mm`) — provides liquidity by quoting other people's RFQs.
- **Taker** (`kalshi_mlb_rfq`) — sends RFQs to take liquidity.

It answers four questions at a glance: **how many RFQs, how many fills, what's
outstanding, and why am I not getting filled.** It never writes to any bot DB and
imports no bot code, so it is safe to run against the live bots.

## Run

```bash
kalshi_mlb_monitor/run.sh
# → http://127.0.0.1:8092
```

Options (env vars):

| Var | Default | Meaning |
|-----|---------|---------|
| `KALSHI_MLB_MONITOR_PORT` | `8092` | Port to serve on |
| `KALSHI_MLB_MONITOR_ROOT` | `~/NFLWork` | Repo root holding the two bot dirs (where the live DBs are) |

The monitor always reads the **live** DBs under `~/NFLWork`, even if launched from a
git worktree — set `KALSHI_MLB_MONITOR_ROOT` only if your checkout lives elsewhere.

Requires `dash`, `plotly`, `duckdb`, `pandas` (already in the project's Python env).

## What each tab shows

Pick a bot (Maker / Taker) and a time window (1h / 24h / 7d / all) in the header.
The status strip shows whether the bot is running, how fresh the data is, dry-run
vs live, and a **STALE** banner if a bot is down.

| Tab | Purpose |
|-----|---------|
| **Overview** | Headline KPIs, a conversion **funnel** (Seen → In-scope → Quoted → Filled for the maker; Sent → Evaluated → Filled for the taker), and **SGP feed health** per book × path. |
| **Why Not Filled ★** | The core view. Decision/reason breakdown (bar + share table with a plain-language legend) and a stacked time-series so you see *when* a reason spikes. |
| **Fills & P&L** | Recent fills table + cumulative fills/stake. Friendly empty-state when a bot has 0 fills. |
| **Positions & Exposure** | Open positions with exposure, working-order status counts (flags **orphaned** `open` orders when the bot is down), and working-order detail. |
| **Adverse Selection** | Maker: per-fill fair drift (`fair_at_confirm − fair_at_quote`) vs the quoted 5% margin — the v1 "does the margin survive adverse selection?" question. Taker: accept vs walk vs halt over time. |

## SGP feed health (Overview)

Per book × path: 🟢 healthy · 🟡 degraded (1+ consecutive failed fetches) ·
🔴 DEAD (`BOOK_ALERT_STREAK`, default 3), with last-fetch age, 24h
%-answered, and the newest `error_class` (`BookTransportError:events:403`
tells DK's Akamai block apart from Novig's DNS death).

- Reads `sgp_fetch_health` from each bot's **market DB** (`bots.py::market_db`)
  through the shipped queries in `kalshi_common/fetch_health_queries.sql` — the
  panel never re-derives the health predicate, so it cannot disagree with what
  the bot alerts on (issue #37).
- **`empty` is healthy**, not a failure: the book answered and had no markets
  (off-day, thin slate). The column reads "% answered" for that reason.
- A book that died days ago still appears (the streak query looks back 7 days,
  the mix only 24h) — the longest outage must not be the most invisible one.
- No `sgp_fetch_health` table (a bot that hasn't run since #38) renders as
  "no history yet", never as a lock error.

## How it works

- **Lock-safe reads.** Every query opens a short-lived `read_only=True` DuckDB
  handle. DuckDB serves the last committed checkpoint, so the dashboard reads fine
  while the live maker holds the write lock. On a rare transient lock the query
  returns a `LOCKED` sentinel and the page keeps its last render.
- **Per-bot adapter** (`bots.py`) maps each bot's table/column names to logical
  fields, so the render code in `app.py` stays generic across the two schemas.
- **Reason vocabularies are read from the data**, never hardcoded — new decision
  codes appear automatically. `bots.py::REASON_GLOSS` adds human descriptions.
- **Naive-local timestamps.** The dashboard operates in naive local wall-clock:
  the maker's tables are `TIMESTAMPTZ` (migrated 2026-07) and `queries._read`
  converts them to naive-local on read; the taker's tables are still naive-local
  natively. Time-window cutoffs are computed in Python and passed as parameters.

## Files

| File | Role |
|------|------|
| `app.py` | Dash layout + callbacks (entrypoint, `python -m kalshi_mlb_monitor.app`) |
| `queries.py` | Lock-safe read-only DuckDB query helpers |
| `bots.py` | Per-bot adapter config + reason-code glossary |
| `run.sh` | Launcher |

## Troubleshooting

- **Everything shows "—" / empty:** check the data root points at the live DBs
  (`KALSHI_MLB_MONITOR_ROOT`) and that the bot DBs exist under it.
- **Taker tab looks frozen:** the taker bot is likely down — the status strip will
  say STALE. Its data is the last historical state until you restart the bot.
- **Port in use:** set `KALSHI_MLB_MONITOR_PORT` to a free port.
