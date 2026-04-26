# MLB Correlated Parlay Auto-Placement (Wagerzon) — Design Spec

**Author:** brainstorm session 2026-04-26
**Status:** Draft pending user review
**Scope:** Wagerzon-only. MLB correlated parlays only. One-click placement from the MLB answer key dashboard.

---

## 1. Mission & Problem

The MLB answer key dashboard surfaces correlated-parlay +EV opportunities in `mlb.duckdb::mlb_parlay_opportunities`. Today these are placed manually: the user reads the parlay, switches to Wagerzon, rebuilds the legs by clicking through the schedule UI, fills the amount, and submits.

A prior automation attempt for the CBB dashboard (`bet_placer/navigator_wagerzon.py`) uses Playwright to scrape the schedule DOM and click odds buttons. It works but is unsatisfactory:

- Brittle DOM scoring against text patterns (`+11-110`, `o169½-110`) breaks when Wagerzon adjusts UI.
- Fuzzy team-name matching with mascot stripping is fragile.
- Browser launch + DOM stabilization adds ~20–30s per bet.
- User must still manually click "Place Bet" at the end (defeats most of the value).
- Singleton-lock cleanup, persistent profiles, and Cloudflare quirks add operational overhead.

**The opportunity:** Wagerzon exposes `PostWagerMultipleHelper.aspx` — a clean REST endpoint that accepts a JSON array of wagers and returns confirmation. Recon (2026-04-26, see `wagerzon_odds/recon_place_parlay.json`) confirmed:

- No CSRF / `__VIEWSTATE` / antiforgery token requirements
- Standard session cookies (already obtained by `wagerzon_odds/parlay_pricer.py`)
- Body is the same `sel` / `detailData` encoding the pricer already builds
- Native bulk support — N parlays in one POST

This means the entire DOM-clicking middle phase can be deleted. Selection, pricing, and placement all become REST calls.

---

## 2. Scope

### In scope (this PR)
- A new Python module `wagerzon_odds/parlay_placer.py` that places parlays via the Wagerzon REST API.
- A new `/api/place-parlay` endpoint on the MLB dashboard server.
- A "Place" button on each parlay row in the MLB dashboard's parlay table.
- A "Dry run" toggle in the dashboard that exercises everything except the final POST.
- `placed_parlays` and `placement_orphans` tables in `mlb_dashboard.duckdb`.
- Failure surfacing: rejected/price-moved/auth/network errors visible per-row.

### Out of scope (explicit non-goals — listed as follow-ups in §9)
- Bulk approve UI ("Place selected" across many parlays)
- Autonomous cron placement above an edge threshold
- "Force place at current price" override button
- Sheets dual-write at placement time (the existing `bet_logger/scraper_wagerzon.py` covers Sheets via Wagerzon's HistoryHelper feed)
- Other books (DK, FD, ProphetX, Novig — they're read-only fair-value sources)
- Any change to the parlay pricer or `mlb_parlay_opportunities` schema

---

## 3. Architecture (Option A — Pure API)

### Components at a glance

```
┌──────────────────────────┐
│  MLB Dashboard (R/Shiny) │   ← user clicks "Place" on a parlay row
└────────────┬─────────────┘
             │ POST /api/place-parlay {parlay_hash}
             ▼
┌──────────────────────────────────────────┐
│  mlb_dashboard_server.py                 │
│  /api/place-parlay handler               │
│   1. Idempotency check (placed_parlays)  │
│   2. Load row from mlb_parlay_opps       │
│   3. Build ParlaySpec                    │
│   4. parlay_placer.place_parlays([spec]) │
│   5. Write result to placed_parlays      │
└────────────┬─────────────────────────────┘
             │
             ▼
┌──────────────────────────────────────────┐
│  wagerzon_odds/parlay_placer.py          │
│   • _api_login() (reused from pricer)    │
│   • ConfirmWagerHelper (price preflight) │
│   • Zero-tolerance drift check           │
│   • PostWagerMultipleHelper (the POST)   │
│   • Parse Confirm/IDWT/TicketNumber      │
└────────────┬─────────────────────────────┘
             │
             ▼
   backend.wagerzon.com
   (REST, single session, no browser)
```

### Wagerzon endpoints (all already understood)

| Endpoint | Purpose | Used by |
|---|---|---|
| `Default.aspx` (POST form) | Login → session cookies | `_api_login()` |
| `ConfirmWagerHelper.aspx` | Returns exact `Win/Risk` for a proposed parlay (preview, no placement) | Price-drift check |
| `PostWagerMultipleHelper.aspx` | **Places** a wager. Returns `Confirm: true`, `IDWT`, `TicketNumber` | Final placement |

### Why in-process (not subprocess like CBB)

CBB's `/api/place-bet` spawns `placer.py` as a child because Playwright owns a stateful browser. With pure-API, placement runs in the dashboard server's request handler thread. Faster (~1–2s vs 20–30s), simpler, no profile/SingletonLock plumbing.

---

## 4. Components

### 4.1 `wagerzon_odds/parlay_placer.py` (new)

```python
@dataclass
class ParlaySpec:
    parlay_hash: str
    legs: list[Leg]              # each Leg: idgm, play, points, odds, pitcher
    amount: float                # = recommended_size from mlb_parlay_opportunities
    expected_win: float          # from dashboard's wz_odds × amount
    expected_risk: float         # = amount

@dataclass
class PlacementResult:
    parlay_hash: str
    status: str                  # placed | price_moved | rejected | auth_error
                                 # | network_error | would_place (dry_run)
    ticket_number: str | None
    idwt: int | None
    actual_win: float | None
    actual_risk: float | None
    error_msg: str               # empty on success
    error_msg_key: str           # Wagerzon's error key, when applicable
    raw_response: str | None     # for orphan / forensics

def place_parlays(specs: list[ParlaySpec], dry_run: bool = False) -> list[PlacementResult]:
    ...
```

**Internals:**
1. Login if no warm session (reuse session-cache pattern from `parlay_pricer.py`).
2. For each spec: call `ConfirmWagerHelper` → assert `abs(returned_win - expected_win) <= $0.01` → if any spec fails, mark *that one* as `price_moved` and exclude from the placement batch.
3. Remaining specs: build `postWagerRequests` JSON array → POST to `PostWagerMultipleHelper.aspx`.
4. Parse response array (one `WagerPostResult` per spec). Map per-spec `ErrorMsgKey` → status.
5. If `dry_run=True`: skip step 3, return `would_place` results with the verified Win/Risk.

**Imports/reuses from `parlay_pricer.py`:** `_api_login`, leg encoding helpers (play codes 0/1=spread, 2/3=O/U, 4/5=ML), period handling for FG vs F5.

### 4.2 `Answer Keys/MLB Dashboard/mlb_dashboard_server.py` (modify)

Add endpoint `POST /api/place-parlay`:

```
Request:  {parlay_hash: str, dry_run: bool = false}
Response: {status: str, ticket_number: str|null, error_msg: str|null}
```

Handler (live mode):
1. Check `placed_parlays` for an existing row with this `parlay_hash` and status in `("placing", "placed")` → return that row (idempotent).
2. Insert `placed_parlays` row with `status="placing"`.
3. Load parlay row from `mlb.duckdb::mlb_parlay_opportunities` (read-only attach).
4. Build `ParlaySpec` from row fields (legs, amount, expected payout).
5. `parlay_placer.place_parlays([spec])` → `PlacementResult`.
6. Update `placed_parlays` row to terminal status. On orphan condition (Wagerzon confirmed but local write fails), additionally insert into `placement_orphans`.
7. Return JSON to dashboard.

Handler (dry_run=true mode): the placement is a no-op against `placed_parlays`. Skip steps 1, 2, 6 entirely; load the parlay (step 3), build the spec (step 4), call `place_parlays([spec], dry_run=True)`, and return the result. The `would_place` status is purely a transient API response — it never lands in the database. This means dry-running a parlay does not block a subsequent live placement.

### 4.3 `Answer Keys/MLB Dashboard/mlb_dashboard.R` (modify)

In `create_parlays_table()` (currently around `mlb_dashboard.R:246-330`):
- Add a "Place" column rendering one of:
  - A button labeled "Place" if no row in `placed_parlays` for this hash.
  - Status text + ticket number if already placed (e.g. `placed · #212147131`).
  - Error pill if failure (e.g. `price_moved`, `rejected: insufficient balance`). Row stays visible.
- Add a "Dry run" toggle (checkbox) above the table. When on, the Place button text becomes "Dry run" and the API call passes `dry_run=true`.

Status update: the dashboard already polls / refreshes the parlay table. Reuse that mechanism rather than introducing new live-update plumbing.

---

## 5. Data Model

Both tables live in `mlb_dashboard.duckdb` (NOT `mlb.duckdb` — placement state is dashboard state). Connection guard: every connection uses `on.exit(dbDisconnect(...))` per repo CLAUDE.md.

### 5.1 `placed_parlays`

```sql
CREATE TABLE IF NOT EXISTS placed_parlays (
  parlay_hash       VARCHAR PRIMARY KEY,   -- from mlb_parlay_opportunities
  status            VARCHAR NOT NULL,
  combo             VARCHAR,
  game_id           VARCHAR,
  game_time         TIMESTAMP,

  -- intent
  recommended_size  DOUBLE,
  expected_odds     INTEGER,
  expected_win      DOUBLE,

  -- result
  actual_size       DOUBLE,
  actual_win        DOUBLE,
  ticket_number     VARCHAR,
  idwt              BIGINT,

  -- forensics
  legs_json         VARCHAR,
  error_msg         VARCHAR,
  error_msg_key     VARCHAR,

  placed_at         TIMESTAMP DEFAULT CURRENT_TIMESTAMP,
  updated_at        TIMESTAMP DEFAULT CURRENT_TIMESTAMP
);
```

**Status state machine:**
```
(no row) → placing → placed                          (success)
                   → price_moved | rejected | auth_error
                     | network_error | orphaned      (terminal failure)
```

`would_place` (dry-run) is intentionally not in this diagram — it is returned by the API but never written to `placed_parlays`. See §4.2 for handler behavior in dry_run mode.

The `parlay_hash` PK gives idempotency on accidental double-clicks: a second concurrent insert fails. We do not retry within the handler (per §6).

If the existing `placed_parlays` schema in `mlb_dashboard.duckdb` is missing any of these columns, add them via `ALTER TABLE ... ADD COLUMN IF NOT EXISTS`. No destructive migration.

### 5.2 `placement_orphans`

Only written when Wagerzon confirms placement but our local DB write fails. Rare; documents real-money commitments that have no local audit row.

```sql
CREATE TABLE IF NOT EXISTS placement_orphans (
  idwt          BIGINT PRIMARY KEY,
  ticket_number VARCHAR,
  parlay_hash   VARCHAR,
  raw_response  VARCHAR,            -- full JSON for forensics
  error         VARCHAR,            -- the local DB write error
  created_at    TIMESTAMP DEFAULT CURRENT_TIMESTAMP
);
```

Reconciliation is manual: cross-reference against `bet_logger`'s Wagerzon HistoryHelper feed.

### 5.3 What does *not* change
- `mlb.duckdb::mlb_parlay_opportunities` — read-only for this feature.
- `bet_logger/` — Wagerzon's HistoryHelper picks up auto-placed parlays naturally; auto-placed bets appear in the bet log identically to manually-placed ones.

---

## 6. Failure Handling

Five failure classes, each terminal in this PR (no retry, no override).

### 6.1 Price drift (the user's hard constraint)
ConfirmWagerHelper returns the exact `Win/Risk` Wagerzon will honor. We compute expected payout as `recommended_size × decimal(wz_odds)`. If `abs(returned_win - expected_win) > $0.01`, abort.

- **Status:** `price_moved`
- **Error column:** `Expected $30.00, Wagerzon offered $28.50`
- **UI:** row stays visible. The dashboard's edge filter will drop it on next refresh if the new price makes it sub-threshold.

### 6.2 Wagerzon API rejection
Response carries `ErrorMsgKey` / `ErrorMsg` / `ErrorCode` populated. Map known keys to user-facing statuses:
- `insufficient_funds` → `rejected: insufficient balance`
- `bet_too_large` / limit-related → `rejected: exceeds limit`
- line-pulled / event-started → `rejected: line pulled`
- Anything else → `rejected: <ErrorMsgKey>` (raw passthrough so we discover what else exists)

### 6.3 Session/auth failure
Detect when Wagerzon returns an HTML login page instead of JSON (response `Content-Type` is not `application/json`, OR response body parses as HTML). On detection: re-run `_api_login()` to obtain a fresh session, then retry the **same** request (ConfirmWagerHelper preflight or PostWagerMultipleHelper, whichever triggered the auth failure) **once**. If the retry also returns HTML, status is `auth_error`, no further retries. The re-login is bounded to one attempt because a persistent auth failure means credentials or account state changed — silently looping won't fix it.

Why this is safe even on the placement POST: ConfirmWagerHelper's preflight runs first; if the session was already dead, the failure surfaces there and we re-login + retry the *preflight*, not the placement. If the session dies between preflight and placement (extremely rare), the retry on the placement POST is acceptable because Wagerzon's earlier preflight rejected nothing — but to be defensive, the retry on PostWagerMultipleHelper *also* re-runs ConfirmWagerHelper before the second placement attempt.

### 6.4 Network failure
`requests.exceptions.Timeout` / `ConnectionError` → status `network_error`, error column shows the exception class. **No automatic retry**: silently retrying network failures is how you accidentally double-place.

**Important ambiguity:** a network error on the PostWagerMultipleHelper request does NOT tell you whether the bet went through. Wagerzon may have received and processed the request, then the response was lost in transit. Before manually clicking Place again on the same parlay, verify against Wagerzon's history (and against `bet_logger/scraper_wagerzon.py`'s feed) that the bet was not placed. The error column should make this explicit so the user doesn't auto-retry.

### 6.5 Local DB write after Wagerzon confirmed
The danger case: bet is placed but `placed_parlays` write fails. Wrap the post-placement DB write in try/except. On failure, write to `placement_orphans` with `IDWT`, `TicketNumber`, raw response, and the local DB error. Log loudly to stderr.

---

## 7. Testing & Rollout

### 7.1 Pre-flight
1. **Recon JSON scrubbed** — `wagerzon_odds/recon_place_parlay.json` had the password in plaintext (`confirmPassword` field and login form body); replaced both occurrences with `REDACTED` on 2026-04-26. The file is gitignored (`.gitignore:38: recon_*.json`).
2. **F5 encoding verification** — the recon captured a full-game parlay (`Period: 0`, `IdGameType: 18`). The MLB pricer also produces F5 parlays. Confirm via `wagerzon_odds/parlay_pricer.py` how F5 is encoded for ConfirmWagerHelper; lift the same encoding into placement. (Existing pricer already handles F5, so this is a code-reading task, not a discovery task.)

### 7.2 Test stages

**Stage 1 — unit tests (mocked).** A `wagerzon_odds/test_parlay_placer.py` with one fixture-based test per scenario. Mocks `requests.Session.post`. Cases:
- success response → `status=placed`, ticket number parsed
- `ErrorMsgKey: "insufficient_funds"` → `status=rejected`, error_msg populated
- ConfirmWagerHelper returns `Win` differing > $0.01 → `status=price_moved`, no POST attempted
- HTML login page returned → one retry → success
- HTML login page returned → retry also fails → `status=auth_error`

**Stage 2 — dry-run via dashboard.** Toggle "Dry run" on. Click Place on real parlays. Verify:
- ConfirmWagerHelper is called against real Wagerzon
- Price-drift check fires correctly when a parlay's odds have shifted
- UI status column updates to `would_place` with the verified Win/Risk
- No actual bet is placed (verify against Wagerzon's history)

**Stage 3 — first real bet via dashboard.** Toggle "Dry run" off. Place ONE parlay at the **$15 minimum**. Verify the full round-trip:
- Wagerzon shows the bet in their history
- `placed_parlays` row exists with correct `ticket_number`, `idwt`, `actual_win`
- `bet_logger/scraper_wagerzon.py` picks it up on next run and pushes to Sheets
- Dashboard renders the row as `placed · #<ticket>`

**Important constraint (per user):** dry-run validates the placement *logic*, but the only real acceptance test is clicking Place on the real dashboard for a real $15 parlay. Don't merge without Stage 3.

### 7.3 Rollout sequence

1. Build on a feature branch in a worktree (per §8).
2. Stages 1 → 2 → 3 in order. No skipping.
3. After Stage 3 passes: pre-merge review per repo CLAUDE.md (data integrity, resource safety, edge cases, logging, security). Document findings.
4. Explicit user approval → merge to `main`. Never merge without approval.
5. Clean up worktree + delete branch.

---

## 8. Version Control & Documentation

### Branch & worktree
- Branch: `feature/mlb-parlay-auto-place`
- Worktree created via `/worktree` before any code changes.
- DuckDB rule: never symlink databases into the worktree (per repo CLAUDE.md). Test from `main` after merge for any DuckDB-touching validation, or `cp` the DB into the worktree.

### Files created / modified
| Path | Change |
|---|---|
| `wagerzon_odds/parlay_placer.py` | NEW — REST client for placement |
| `wagerzon_odds/test_parlay_placer.py` | NEW — unit tests with mocked responses |
| `wagerzon_odds/recon_place_parlay.py` | NEW — already created (commit alongside) |
| `Answer Keys/MLB Dashboard/mlb_dashboard_server.py` | MODIFY — add `/api/place-parlay` |
| `Answer Keys/MLB Dashboard/mlb_dashboard.R` | MODIFY — add Place column + Dry run toggle |
| `wagerzon_odds/README.md` | UPDATE — document parlay_placer + recon script |
| `Answer Keys/MLB Dashboard/README.md` (if exists, else create) | UPDATE — Place button + dry-run toggle docs |

### Commit hygiene
- Source code only. The recon JSON stays gitignored.
- Documentation updates ship in the same merge to `main` as the feature.
- No skipped hooks. Pre-merge review documented before approval.

---

## 9. Follow-up: hooks for options 2 + 3

These are not built in this PR but the architecture supports them with minimal work:

### Bulk approve
- `parlay_placer.place_parlays()` already takes a list; the Wagerzon endpoint is natively bulk.
- `/api/place-parlay` extends to accept `{parlay_hashes: [...]}`. Server filters out hashes already in `placed_parlays` (idempotency), packs the rest into one `PostWagerMultipleHelper` call.
- Dashboard adds a checkbox column + "Place selected" button.

### Autonomous cron
- New script `Answer Keys/MLB Dashboard/place_parlays_cron.R` (or `.py`) imported on cron schedule.
- Query: `SELECT * FROM mlb_parlay_opportunities op LEFT JOIN placed_parlays pp ON op.parlay_hash = pp.parlay_hash WHERE pp.parlay_hash IS NULL AND op.edge_pct > <threshold>`.
- Calls `parlay_placer.place_parlays()` directly. Same code path as one-click.
- Add an "auto-place on/off" master switch and a `min_edge_pct` config in the dashboard before enabling.

### Override / force-place
- Add an "override" column to the placement endpoint that, when true, skips the §6.1 price-drift abort and accepts whatever Wagerzon currently quotes.
- Useful only after several weeks of validation showing how often legitimate price moves are still profitable.

---

## 10. Risks & open questions

### Known risks
- **Password in placement body**: `confirmPassword` is sent in every `PostWagerMultipleHelper` POST. Stored in `WAGERZON_PASSWORD` env var. Don't log full request bodies. Recon JSON has been scrubbed.
- **Single-session correctness**: if multiple `/api/place-parlay` requests arrive concurrently (race), they could both pass the idempotency check before either inserts. Mitigation: the Flask dashboard server processes requests serially per-thread; if multi-threaded, wrap the idempotency check + insert in a transaction with a unique constraint on `parlay_hash` (DuckDB enforces PK).
- **Orphan window**: between Wagerzon's confirmation and our local DB write is a few-millisecond window where the bet exists at Wagerzon but not locally. Mitigated by §6.5 (orphans table). Acceptable.
- **F5 encoding untested in placement**: the recon was a full-game parlay. Existing pricer encodes F5 for ConfirmWagerHelper; verify it works for PostWagerMultipleHelper in Stage 2 of testing.

### Open questions (to surface during implementation)
- **Does `mlb_parlay_opportunities` store `idgm` per parlay?** Wagerzon's internal game ID is required to build the `sel` and `detailData` payloads (`play_idgm_points_odds`). The pricer (`wagerzon_odds/parlay_pricer.py`) already calls ConfirmWagerHelper, so it must obtain `idgm` somehow — either it's already in `mlb_parlay_opportunities`, or it's joined from another table at pricing time. The placer must use the same path. **Resolve in implementation step 1**: read the pricer to determine where `idgm` comes from, then either (a) lift it into `ParlaySpec` from the same source, or (b) add `idgm` to `mlb_parlay_opportunities` if it's currently being computed and discarded.
- Does the existing `placed_parlays` schema in `mlb_dashboard.duckdb` already have all the columns in §5.1, or do we need `ALTER TABLE`? Resolve in implementation step 1.
- Are there Wagerzon `ErrorMsgKey` values we don't yet know about that need a friendlier mapping in §6.2? Surfaced organically by the raw passthrough; revisit after some real-world placement runs.
- Does the dashboard's polling cadence give a smooth-enough UX for status updates, or do we need a more responsive update path? Revisit if the placement flow feels laggy.

---

## Appendix A: Wagerzon `PostWagerMultipleHelper` request shape

Captured 2026-04-26 (`wagerzon_odds/recon_place_parlay.json`). Body is `application/x-www-form-urlencoded`; `postWagerRequests` is URL-encoded JSON.

```
POST https://backend.wagerzon.com/wager/PostWagerMultipleHelper.aspx
Content-Type: application/x-www-form-urlencoded
Cookie: <session>

postWagerRequests=[
  {
    "WT": "1",
    "open": 0,
    "IDWT": "0",
    "sel": "5_5632938_0_-140,3_5632938_7.5_-130",
    "sameAmount": false,
    "amountType": "0",
    "detailData": "[{\"Amount\":\"15\",\"RiskWin\":0,...}]",
    "confirmPassword": "<WAGERZON_PASSWORD>",
    "sameAmountNumber": "15",
    "useFreePlayAmount": false,
    "roundRobinCombinations": ""
  }
]
```

`sel` encoding (per leg, comma-separated): `play_idgm_points_odds`
- play: 0=away spread, 1=home spread, 2=over, 3=under, 4=away ML, 5=home ML
- points: line value (negative for favorite spread / over)
- odds: integer American (e.g. `117` = +117, `-130` = -130)

## Appendix B: Wagerzon `PostWagerMultipleHelper` response shape

```json
{
  "result": [{
    "WagerPostResult": {
      "details": [<one entry per leg>],
      "IDWT": 438173,
      "TicketNumber": "212147131",
      "Risk": 15.0,
      "Win": 30.0,
      "WagerType": 1,
      "WagerTypeDesc": "PARLAY (2 TEAMS)",
      "Confirm": true,
      "ErrorCode": {},
      "ErrorMsg": "",
      "ErrorMsgKey": "",
      "ErrorMsgParams": ""
    },
    "UpdateLines": null
  }]
}
```

Success indicator: `Confirm: true` AND `ErrorMsgKey == ""`. Failure: `ErrorMsgKey` populated.
