# Wagerzon Multi-Account Support — MLB Correlated Parlay Dashboard

**Date:** 2026-05-01
**Branch:** `feature/wagerzon-multi-account-dashboard`
**Worktree:** `.worktrees/wagerzon-multi-account-dashboard`

## Goal

Display available balances for all configured Wagerzon accounts on the MLB correlated parlay dashboard, and let the user pick which account each placement goes to via a single global selector.

## Context

Today the MLB correlated parlay dashboard (`Answer Keys/MLB Dashboard/`) automates parlay placement to Wagerzon via `wagerzon_odds/parlay_placer.py`, but is hardwired to one account: `WAGERZON_USERNAME` / `WAGERZON_PASSWORD` from `bet_logger/.env`. A second account (`WagerzonJ`) already exists in env and is scraped by `bet_logger`, and a third (`WagerzonC`) is being created in a separate worktree. The dashboard has no awareness of either secondary account and no concept of available balance.

## Decisions (from brainstorming)

| # | Decision |
|---|---|
| Q1 | Global account selector (one dropdown sets account for all placements until changed). |
| Q2 | Show balances for **all** discovered accounts always. |
| Q3 | Three accounts at launch: `Wagerzon`, `WagerzonJ`, `WagerzonC`. Design supports N. |
| Q4 | Refresh balances on dashboard load + after every successful placement. Manual refresh button as a small affordance. No background polling. |
| Q5 | Soft-warn (not block) on insufficient available balance for the selected account. Warning is transient — recomputes when balance, account, or risk changes. |
| Q6 | Record `account` label on every placed bet for per-account P&L attribution. |
| Q7 | Default account on dashboard load = last-used (persisted across sessions). |
| Q8 | On balance-fetch failure: show last-known balance with `stale Xs ago` tag; allow placement; turn pill red if stale > 10 minutes. |
| Q9 | Scope is **only** the MLB correlated parlay dashboard placement flow. |
| —  | "Available balance" specifically (cash minus open exposure / pending), not raw cash. |

## Architecture overview

**New modules** (in `wagerzon_odds/`):

| File | Responsibility |
|---|---|
| `wagerzon_auth.py` | Login + cached session per `(username, password)`. Extracted from current `parlay_placer.py` so placement and balance share one auth path. |
| `wagerzon_accounts.py` | Account registry. Discovers accounts from env vars matching `WAGERZON{SUFFIX}_USERNAME` / `WAGERZON{SUFFIX}_PASSWORD`. Returns ordered list of `WagerzonAccount`. |
| `wagerzon_balance.py` | `fetch_available_balance(account)` — wraps the Wagerzon balance endpoint. Returns `BalanceSnapshot`. |

**Modified files:**

- `wagerzon_odds/parlay_placer.py` — accepts an explicit `account: WagerzonAccount`; no env reads at module top.
- `Answer Keys/MLB Dashboard/mlb_dashboard_server.py` — three new endpoints, one modified.
- `Answer Keys/MLB Dashboard/mlb_dashboard.R` — header pills + selector + transient warnings.

### Data flow on dashboard load

```
dashboard.R loads
  → JS: GET /api/wagerzon/balances
       → server: registry.list_accounts()
       → balance.fetch_all(accounts)        (parallel via ThreadPoolExecutor)
       → returns [{label, available, cash, fetched_at, error, stale_seconds}, ...]
  → JS: GET /api/wagerzon/last-used → "WagerzonJ"
  → header renders pills + dropdown selected to last-used
```

### Data flow on bet placement

```
User clicks Place on parlay
  → JS: POST /api/place-parlay {parlay: ..., account: "WagerzonJ"}
  → server: registry resolves "WagerzonJ" → WagerzonAccount
       → parlay_placer.place_parlays(parlay, account)
       → on success: INSERT placed_parlays (..., account="WagerzonJ")
       → balance.fetch_available_balance(account)
       → response includes balance_after snapshot
  → JS: pill for WagerzonJ flashes + updates
```

## Account registry

`wagerzon_odds/wagerzon_accounts.py`:

```python
@dataclass(frozen=True)
class WagerzonAccount:
    label: str        # "Wagerzon", "WagerzonJ", "WagerzonC"
    suffix: str       # "", "J", "C"
    username: str
    password: str

def list_accounts() -> list[WagerzonAccount]:
    """Discovers accounts from env. Always returns primary first (if present),
    then alphabetical by suffix."""
```

**Discovery rule:**

- Scan env for `WAGERZON_USERNAME` (primary, suffix `""`, label `"Wagerzon"`) and any `WAGERZON{X}_USERNAME` where `{X}` is a single uppercase letter.
- For each match with both `_USERNAME` and `_PASSWORD` set, build a `WagerzonAccount` with `label = "Wagerzon" + suffix`.
- Skip pairs missing one half (log a warning).
- Order: primary first, then alphabetical by suffix.

**Why suffix discovery:** matches the existing convention (`WAGERZONJ_*`), zero config files, adding WagerzonC = two env vars + Flask restart.

## Balance fetcher

`wagerzon_odds/wagerzon_balance.py`:

```python
@dataclass
class BalanceSnapshot:
    label: str
    available: float | None      # None when error is set
    cash: float | None           # optional, may be None even on success
    fetched_at: datetime         # UTC
    error: str | None            # one of: timeout, auth_failed, wz_error, None

def fetch_available_balance(account: WagerzonAccount) -> BalanceSnapshot: ...

def fetch_all(accounts: list[WagerzonAccount]) -> list[BalanceSnapshot]:
    """Parallel fan-out; one slow account doesn't block the others.
    Per-account fetch timeout = 5 seconds."""
```

**Endpoint discovery deferred to implementation:** Memory has notes on `ConfirmWagerHelper` (pricing) and `PostWagerMultipleHelper` (placement) but nothing on a balance endpoint. The implementation plan must include a "discover the WZ balance endpoint" task before writing this module — likely via browser devtools while logged into the WZ web UI.

## Auth helper

`wagerzon_odds/wagerzon_auth.py`:

```python
def get_session(account: WagerzonAccount) -> requests.Session:
    """Returns a logged-in session. Caches per-account in-memory for the
    lifetime of the Flask process; transparently re-logins on 401."""
```

In-memory cache only. Flask process restart = fresh logins. No cookie persistence to disk (avoids stale-cookie debugging headaches).

## Dashboard UI

### Header strip (added near top of `mlb_dashboard.R`)

```
┌────────────────────────────────────────────────────────────────────────┐
│  Wagerzon: $1,245.32   WagerzonJ: $890.10 ⚠   WagerzonC: $2,010.55  ↻ │
│  Placing on:  [ WagerzonJ ▾ ]                                           │
└────────────────────────────────────────────────────────────────────────┘
```

- **Pills** — one per account from the registry. Each shows label + available balance. The selected-account pill gets a subtle highlight (bold outline).
- **Stale tag** — if last fetch had `error`, pill shows last-known number with subtitle `stale Xs ago` and ⚠ glyph. Pill flips red if stale > 10 minutes.
- **Refresh button (↻)** — fetches all balances on click.
- **Selector** — single dropdown listing every label from the registry. On change, immediately persists via `POST /api/wagerzon/last-used`.

### Per-parlay placement area

- After a parlay's `risk` is computed, JS compares it against the selected account's `available` balance.
- If `risk > available`: render transient warning beside the Place button — `"⚠ insufficient on WagerzonJ ($890 < $1,000 risk)"`. Place button stays enabled (Q5 → soft warn). Warning recomputes when balance, account, or risk changes.
- If `available` is `null` (no successful fetch ever for that account): warning is **suppressed** — we can't compute it and don't want a misleading message. The pill itself already shows the stale/error state.

### Placement outcomes

- **Success:** the account's pill flashes + updates with `balance_after`. Toast: `"Placed on WagerzonJ — ticket #12345"`.
- **Failure** (`BALANCEEXCEED`, price moved, network error): existing error toast plus account label. No automatic account switch.

## Backend / Flask endpoints

### New

```
GET  /api/wagerzon/balances
       200 {"balances": [
           {"label": "Wagerzon",  "available": 1245.32, "cash": 1300.00,
            "fetched_at": "2026-05-01T15:22:11Z", "error": null,
            "stale_seconds": 0},
           {"label": "WagerzonJ", "available": 890.10,  "cash": 890.10,
            "fetched_at": "2026-05-01T15:22:11Z", "error": null,
            "stale_seconds": 0},
           {"label": "WagerzonC", "available": null,    "cash": null,
            "fetched_at": "2026-05-01T15:21:42Z", "error": "auth_failed",
            "stale_seconds": 29}
         ]}
       Always 200; per-account errors live in the snapshot, not HTTP status.

GET  /api/wagerzon/last-used
       200 {"label": "WagerzonJ"}     // null if never set or label removed

POST /api/wagerzon/last-used
       Body: {"label": "WagerzonJ"}
       200 {"ok": true}
       400 if label not in current registry
```

### Modified

```
POST /api/place-parlay
       Body adds field: "account": "WagerzonJ"   // required
       → server resolves account from registry
       → parlay_placer.place_parlays(parlay, account)
       → on success: INSERT placed_parlays with account label
       → balance.fetch_available_balance(account)
       Response gains: "balance_after": {label, available, cash, fetched_at, error}
       400 if account label missing or not in registry
```

### Server-side caches

- `_balance_cache: dict[str, BalanceSnapshot]` — last attempt per label. Powers the "stale" tag (`stale_seconds = now - fetched_at`, computed at response time). Updated by `/api/wagerzon/balances` and after successful placement.
- `_session_cache` (inside `wagerzon_auth`): one logged-in `requests.Session` per account label.
- Both caches are in-memory only; Flask restart = clean slate.

### Concurrency

- `fetch_all` parallelizes balance fetches via `ThreadPoolExecutor`; one slow account doesn't block the others.
- Per-account fetch timeout = 5 seconds; on timeout the snapshot has `error="timeout"` and the cached prior snapshot serves as the "stale" value.

## Persistence & data model

### Last-used account

New table in `mlb_dashboard.duckdb`:

```sql
CREATE TABLE IF NOT EXISTS dashboard_settings (
    key   TEXT PRIMARY KEY,
    value TEXT NOT NULL,
    updated_at TIMESTAMP DEFAULT CURRENT_TIMESTAMP
);
```

- Last-used account: `key = 'wagerzon_last_used'`, `value = '<label>'`.
- Generic `key/value` shape so future dashboard prefs (sort order, default sport, etc.) reuse it without schema churn.
- **Read on load:** `SELECT value FROM dashboard_settings WHERE key='wagerzon_last_used'`. If missing, or label no longer in registry, fall back to registry primary.
- **Write on selector change:** upsert via `INSERT ... ON CONFLICT (key) DO UPDATE SET value = excluded.value, updated_at = CURRENT_TIMESTAMP`.

### Per-account bet attribution

```sql
ALTER TABLE placed_parlays ADD COLUMN account TEXT;
```

- New rows: written with the resolved account label.
- Old rows: stay `NULL`. Backfill not safe — historical bets all went to primary, but encoding that as a value masks the absence of explicit recording. Downstream P&L queries should treat `NULL` as "primary" with a SQL comment, OR exclude pre-feature rows.

### Migration

Both schema changes (`CREATE TABLE dashboard_settings`, `ALTER TABLE placed_parlays ADD COLUMN account`) run idempotently at server startup — `IF NOT EXISTS` for the table, try-except on the `ALTER` (DuckDB raises if column exists).

### Not tracked (deferred YAGNI)

Balance-before / balance-after snapshots per bet. bet_logger already reconciles WZ history per account. Easy to add a column later if a query demands it.

### bet_logger interaction (no work needed)

`bet_logger/scraper_wagerzon.py` already scrapes per-account (primary → Sheet1, `--account j` → Shared tab). When the dashboard places a bet on WagerzonJ or WagerzonC via this feature, the next scheduled bet_logger run picks up that ticket from the appropriate account's WZ history and lands it in the right sheet tab automatically. No bet_logger code changes required, but `--account c` support in `bet_logger/scraper_wagerzon.py` (and a `run_all_scrapers.sh` entry for it) is a prerequisite — being delivered by the parallel WagerzonC worktree, not this one.

## Failure modes & edge cases

### Balance-fetch failures

| Failure | Behavior |
|---|---|
| Network timeout (5s) | Snapshot `error="timeout"`. UI shows last-known number with `stale Xs ago`. Placement still allowed. |
| Auth 401 | `wagerzon_auth.get_session()` re-logins transparently and retries once. Re-login fail → `error="auth_failed"`; same stale-display behavior. |
| WZ 5xx / outage | `error="wz_error"`. Same stale-display. Pill flips red if stale > 10 minutes. |
| Account env misconfigured (one of pair missing) | Account omitted from registry; never appears in UI. Warning logged at startup. |
| Last-used label no longer in registry | Falls back to primary; selector lands on primary; user re-picks. |

### Placement failures

| Failure | Behavior |
|---|---|
| `BALANCEEXCEED` from WZ despite pre-warn | Existing handling; toast `"Rejected by WagerzonJ: balance exceeded"`. Server triggers fresh balance fetch for that account. |
| Price moved | Existing handling; toast labels account. No balance change. |
| Network error mid-placement | Existing `orphaned` status; toast labels account. Server triggers balance refresh to reveal whether the bet landed. |
| Account label not in registry (env changed mid-session) | 400 with clear message. UI re-fetches `/api/wagerzon/balances` to re-render selector. |

### UI edge cases

- **Zero accounts discovered:** selector reads `"No Wagerzon accounts configured"`, Place button disabled with that tooltip. Server logs warning at startup.
- **One account discovered:** selector renders for visual consistency, locked to that label. Header still shows the pill.
- **Account added at runtime:** out of scope — registry is read at startup. Restart Flask to pick up new accounts. Documented in README.

### Concurrency

- Flask placement is single-threaded already; two placements in quick succession serialize naturally; balance pill updates between them.
- Multi-tab: last-used is per-DB so a second tab inherits the same selection. Selector changes don't push to other tabs (no SSE). Acceptable wart.

## Out of scope

- Account selection in `/api/place-bet` (single-bet flow is tracking-only today — no WZ call to redirect).
- NFL Draft portal (`nfl_draft/`) account selection — same modules will be reusable when wanted.
- Background balance polling (rejected in Q4).
- Per-account auto-suggestion or auto-switching on insufficient balance (rejected in Q5).
- Backfill of historical `placed_parlays.account` values.
- Balance-before/after snapshots per bet.
- Hot-reload of accounts without Flask restart.

## Version control plan

- **Branch:** `feature/wagerzon-multi-account-dashboard`
- **Worktree:** `.worktrees/wagerzon-multi-account-dashboard` (already created)
- **DuckDB caveat:** do NOT symlink `mlb_dashboard.duckdb` into the worktree (CLAUDE.md rule — DuckDB WAL data lost on worktree removal). Either copy the DB into the worktree for local testing, or test the schema migration after merging to main.
- **Commit structure** (logical):
  1. `feat(wagerzon): add account registry + auth/balance modules` — new files only, no integration
  2. `refactor(wagerzon): parlay_placer takes explicit account parameter` — `parlay_placer.py` only
  3. `feat(mlb-dashboard): balance display, account selector, persistence` — server + R + JS + DB migration
  4. `docs: README + CLAUDE.md updates for multi-account flow`
- **Pre-merge:** executive engineer review of `git diff main..HEAD` per CLAUDE.md (data integrity, resource safety, edge cases, dead code, log/disk hygiene, secrets). Run dashboard end-to-end on at least 2 accounts before merging. **Explicit user approval required before merge.**

## Documentation updates (same merge to main)

- `wagerzon_odds/README.md` — new "Multi-account support" section describing the registry, env-var convention, balance fetcher, parlay_placer's `account` parameter.
- `wagerzon_odds/CLAUDE.md` — if it exists; else add a note in `bet_logger/CLAUDE.md` cross-referencing.
- `Answer Keys/MLB Dashboard/CLAUDE.md` (or parent `Answer Keys/CLAUDE.md`) — note new endpoints, `dashboard_settings` table, `placed_parlays.account` column, and how to add a new account.
- `bet_logger/.env.example` — already includes `WAGERZONC_USERNAME` / `WAGERZONC_PASSWORD` (added by the parallel WagerzonC worktree). No change needed.
- Root `CLAUDE.md` — likely no change.

## Auto-memory updates after ship

- Add `wagerzon_multi_account.md` reference memory: registry pattern, where balance lives, how to add an account.
- Update navigator/Wagerzon notes if the auth refactor changes anything for future scrapers.
