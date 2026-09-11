# Unabated Ticket — bet history core: flag lines and games already bet (#114)

**Date:** 2026-09-11 · **Branch:** `feature/unabated-bet-history` · **Worktree:** `.worktrees/unabated-bet-history`
**Extends:** `unabated_ticket/` (#111 ticket, #112 Edges + alerts; #113 alt lines is still awaiting merge)
**Split (user decision 2026-09-11):** this ticket is the **core** — matcher, local bets service with **Kalshi**, panel UI. Each other venue is its own ticket plugging into the record contract defined here: **#115 BetOnline**, **#116 Novig**, **#117 ProphetX**.
**Pikkit is out** (user confirmed 2026-09-11 it does not auto-sync; it only syncs when the phone app is opened, so it cannot see a bet placed a minute ago). Research behind that call is in the appendix.
**Status:** plan only — no code until the user agrees.

## Goal

When a line is in the Ticket tab or the Edges tab, show whether you already
have a bet on it, on the other side of it, or on the same game — and what that
bet was. This ticket ships that end-to-end for Kalshi and leaves a contract
the other venues implement.

## Facts that shape the design (read from the repo, 2026-09-11)

| Area | What exists today | Consequence |
|---|---|---|
| Kalshi | `kalshi_common/auth_client.api()` (signed REST, `configure()` pattern); `/portfolio/fills?limit=100&min_ts=` polled by `kalshi_rfi/state.py` and `unabated_edge/maker/state.py`; `/portfolio/positions`; `kalshi_common/leg_types.parse_suffix_start_utc` + `_parse_event_suffix` for **MLB** event suffixes (`26SEP041410DETCLE`, `G1`/`G2` doubleheader marker) | Reuse the client and the MLB suffix parser as-is. NFL/CFB/NBA/NHL ticker grammar is **not** in the repo — step 0 reads it off the user's own fills. |
| Extension | `page.js` (MAIN world) → `window.postMessage` → `content.js` (ISOLATED) → `chrome.storage.local`; the panel re-renders on `storage.onChanged`; `feed.describeLine` and the captured ticket both carry `league`, `awayTeam`, `homeTeam`, `eventStart`, `betType`, `period`, `sideIndex`, `points`, `rotation` | The matcher consumes the ticket/edge-row shape unchanged; venue content scripts (#116/#117) publish to storage with the same pattern. |
| Ports | 8083 (MLB dashboard), 8090 (draft portal), 8092 (bots monitor), 8093 (worktree dashboard) are taken | Bets service on **8094**. |

## Architecture

```
unabated_ticket/
  bets_service/                      NEW — local Python service, 127.0.0.1:8094
    service.py                         entry point: ThreadingHTTPServer + one poll thread; GET /bets.json, GET /health
    store.py                           DuckDB bets.duckdb: `bets` (upsert on bet id), `source_runs` (per-source freshness log)
    normalize.py                       shared record helpers (american <-> cents, side/points conventions); python-tested
    sources/__init__.py                Source protocol: name, poll_sec, fetch() -> list[record]; the plug-in point for #115/#116/#117
    sources/kalshi.py                  auth_client.configure(); fills + positions -> records
    sources/kalshi_ticker.py           series -> betType/period; event suffix -> league/teams/start (MLB from leg_types, others from step 0)
    .env.example                       KALSHI_API_KEY_ID, KALSHI_PRIVATE_KEY_PATH, BETS_SERVICE_PORT=8094
    run.sh
    tests/                             pytest on the recorded fixtures
  extension/
    bets.js                            NEW pure: record shape, game match, tiers, merge, retention (node-tested)
    teams.js                           NEW pure: per-league team-name -> key tables (Unabated names, Kalshi codes)
    panel.html/js/css                  EDIT: Bets tab, Ticket banner, Edges row marker, hide-bet filter, per-source freshness, unmatched list
    manifest.json                      EDIT: host_permissions http://127.0.0.1:8094/*
  tests/
    bets.test.js, teams.test.js        NEW node tests
    fixtures/bets/                     kalshi_fills.json, kalshi_positions.json (real, anonymised: ids scrubbed, stakes kept)
```

**Why a local service.** Kalshi's private key must not live in the extension.
The service is the only place that signs requests; the panel polls
`http://127.0.0.1:8094/bets.json` (host permission for that origin; extension
pages with a host permission are exempt from CORS and mixed-content rules). It
binds 127.0.0.1 only, has no auth, and serves nothing but the user's own bet
list. Same pattern as the dashboards (`kalshi_mlb_monitor/run.sh`): start it
by hand or from a LaunchAgent; the panel says loudly when it is unreachable.
#115 and #117 (API route) add `sources/*.py` modules behind the same protocol;
#116/#117 content-script routes bypass the service and write storage directly.

**Service loop.** One thread; each registered source runs on its own
`poll_sec` (Kalshi 60 s: fills since the last cursor with a 60 s `min_ts`
overlap, dedupe on `trade_id`, plus positions). Each pass writes a
`source_runs` row (`source, started_at, ok, error, n_records`) whether it
succeeded or not, and a failed pass **keeps the previous records** and marks
the source stale — same rule as the leg surface: a dark source must never
blank the list. `/bets.json`:

```json
{ "generatedAt": "...",
  "sources": { "kalshi": {"fetchedAt": "...", "ok": true, "error": null, "count": 12} },
  "bets": [ <records, open + settled within RETENTION_DAYS=30> ] }
```

## Normalised bet record — the contract #115/#116/#117 implement

```
{
  id:            "kalshi:KXNFLGAME-26SEP14DALPHI-PHI:yes"          source-native, stable across refreshes
  source:        "kalshi_api" | "betonline_api" | "novig_api" | "novig_page" | "prophetx_api" | "prophetx_page"
  venue:         "kalshi" | "novig" | "betonline" | "prophetx"
  league:        "nfl" | "cfb" | "nba" | "cbb" | "wnba" | "mlb" | "nhl" | "soccer" | null    (feed.LEAGUES paths)
  eventStart:    ISO UTC | null            (Kalshi suffix; null when the venue does not give it)
  awayTeam, homeTeam:  raw names as the venue wrote them | null
  awayKey, homeKey:    teams.js keys ("nfl:phi") | null when unresolved
  rotation:      number | null            (only some venues carry it)
  betType:       "moneyline" | "spread" | "total" | "other"
  period:        "FG" | "1H" | "2H" | "1Q".."4Q" | "F5" | "I1"
  side:          "away" | "home" | "over" | "under" | null
  points:        number | null            (the side's own number, as the ticket shows it — no home-perspective sign)
  price:         American int             (Kalshi: from the VWAP fill price in cents)
  stake:         USD risked               (Kalshi: contracts × price)
  toWin:         USD | null
  contracts:     number | null            (exchanges)
  placedAt:      ISO UTC                  (Kalshi: first fill)
  status:        "open" | "won" | "lost" | "push" | "void" | "closed" | "unknown"   ("closed" = position sold to 0)
  isParlayLeg:   bool, parlayId: string | null, legIndex, legCount
  approx:        [ "kalshi_no_side_includes_tie", "game_date_unknown", ... ]   reasons the match is weaker than it looks
  sourceFetchedAt: ISO
  raw:           { trimmed venue fields } (for the unmatched list and debugging, never for matching)
}
```

One source per venue, so **no cross-source dedupe**; within a source, records
dedupe on their native id. Whatever a venue does not provide stays `null`,
never guessed.

**Kalshi specifics.** A bet is a *position*, not a ticket: aggregate fills by
`(ticker, side)` (VWAP price, summed contracts, first `created_time`), net
buys against sells, and cross-check the net against `/portfolio/positions`
(Kalshi is the source of truth — the bots' lesson). Net 0 ⇒ `status: closed`
(kept, never matched). Series → `betType`/`period`: `KXMLBGAME` ml,
`KXMLBSPREAD`, `KXMLBTOTAL`, `KXMLBF5*`, `KXMLBRFI` (I1 total at 0.5) are
known; other sports' series come from step 0 and the parser **fails closed**
on an unknown series (`league: null`, listed as unmatched with the raw
ticker). Side: YES on a team market = that team; **NO on a team market = the
other team with `approx: ["kalshi_no_side_includes_tie"]`** for sports with
ties (NFL/soccer; MLB has none). Spread/total strikes come from the market
ticker suffix, converted once, in the service, to the side's own number (same
sign convention `kalshi_common.legset` uses; a test per sign). Event, league
and start come from the **event-ticker suffix, never team names** (the
`close_time` note; `G1`/`G2` keeps doubleheaders apart).

## Matching (`bets.js`, pure, node-tested)

**Input line** is the existing ticket / edge-row shape (`league`, `awayTeam`,
`homeTeam`, `eventStartMs`, `betType`, `period`, `sideIndex`, `points`,
`rotation`). Side 0 = away/Over, 1 = home/Under, as everywhere else.

**Game match** (prerequisite for every tier), in order:

1. `league` equal, and
2. both `awayKey`/`homeKey` equal (order-insensitive: a source may swap them)
   **or** `rotation` equals the row's away/home rotation, and
3. time: when the bet has `eventStart`, |Δ| ≤ 30 min (the doubleheader rule
   from `leg_types.SCHEDULE_START_TOLERANCE_MIN`); when it does not, the
   event starts between `placedAt` and `placedAt + 8 days` and the bet is
   flagged `game_date_unknown`. Two candidate events for one bet (a series,
   a doubleheader with no time) ⇒ **no match, listed under "unmatched:
   ambiguous game"** — never a guess.

**Tiers**, per bet (and per parlay leg, labelled "(parlay leg)"), strongest first:

| Tier | Rule | Shown as |
|---|---|---|
| `same_line` | same `betType`, `period`, side, and `points` (ML: points null on both) | "You bet this: Eagles -3.5 -110 · $300 @ Kalshi · Sep 10 2:15 PM" |
| `same_side` | same `betType`, `period`, side; different `points` | "You have Eagles -3.5 -110 (this is -4.5)" |
| `opposite` | same `betType`, `period`; complementary side (away↔home, over↔under), any points | red: "You are on the OTHER side: Cowboys +3.5 -105 · $200 @ Kalshi" (adds "at a different number" when points differ) |
| `same_game` | game matched, any other market/period | "You have a bet on this game: Under 40.5 · $150 @ Kalshi" |

Kalshi NO-side bets match as the other team with the tie caveat in the label
("NO Eagles ≈ Cowboys or tie"). Closed positions and settled bets are
excluded from tiers but stay in the Bets tab list.

**Team keys (`teams.js`).** Per league: Unabated's `teams` names (from the
snapshot) and Kalshi codes (`_MLB_CODE_TO_TEAM` ported for MLB; other sports
from step 0). Unresolvable name ⇒ key `null` ⇒ "unmatched: team not
recognised (<name>)" with the raw text so the table can be extended. Venue
tickets add their own spellings to the same tables. CFB is the risk
(hundreds of teams); start with the FBS names from the Unabated CFB snapshot
and grow from the unmatched list.

## UI

- **Ticket tab** — a banner between the event line and the stake, one line
  per match (strongest first, at most 5 then "+N more"), `opposite` in the
  red `.bad` style. Nothing shown when no match; "bet sources unavailable" in
  the existing warning strip when every source is stale or absent.
- **Edges tab** — a badge per row: `BET` (same_line/same_side), `OTHER SIDE`
  (red), `GAME`. A **"Hide lines I've bet"** checkbox hides `same_line` +
  `same_side` rows only — `opposite` and `same_game` stay visible because
  they are warnings. Alerts skip lines that would be hidden (same setting;
  default off like the filter).
- **Bets tab** (third tab) — (1) per-venue freshness table, one row per
  registered source (`fetchedAt`, count, error; green/amber/red at 5 min /
  60 min; the service itself "unreachable since …" in red; a venue with no
  source yet shows "no source configured", so #115–#117 are visibly pending,
  not silently missing); (2) open bets list (venue, bet, stake, price,
  placed); (3) **unmatched** list with the reason per bet (team not
  recognised, ambiguous game, no event on the board yet, league not on the
  scanner, unknown Kalshi series). Nothing is dropped silently.
- **Header line** under the tabs: "bets: 14 open · kalshi 20 s · betonline —
  · novig — · prophetx —".

## Storage and retention

- `chrome.storage.local`: `betsService` (last `/bets.json` payload +
  `fetchedAt`), `betsSettings` (`serviceUrl` default
  `http://127.0.0.1:8094`, `hideBet`); venue content scripts (#116/#117) add
  `betsNovig` / `betsProphetxPage`. A record is ~1 KB; the extension keeps
  **open bets + settled within 30 days** and prunes on load, so a heavy
  month stays under 1 MB of the 10 MB quota.
- Service DuckDB `unabated_ticket/bets_service/bets.duckdb` (gitignored by
  `*.duckdb`): every record ever seen, unbounded — a few rows a day, and the
  future CLV computation needs the full history. `/bets.json` serves the
  30-day window by default (`?days=`).
- The panel polls the service every 30 s while open (the scanner's
  visibility pause applies), never in the service worker.

## Phases — one commit each

0. **Recon** — done 2026-09-11 (see "Recon findings"); the remaining step
   is saving the anonymised fixture (one fill + one market/event payload per
   in-scope series). Commit:
   `recon(unabated_ticket): Kalshi bet-history fixtures for #114`.
1. **`bets.js` + `teams.js` + node tests** — record shape, team keys, game
   match, tiers, merge, retention, on the phase-0 fixture. Pure files only.
2. **Bets service** — `bets_service/` with the Source protocol + Kalshi,
   DuckDB store, `/bets.json`, `/health`, pytest on fixtures, `.env.example`,
   `run.sh`, README section. Verified against the live account from the
   worktree (the service reads no repo DuckDB, so no seeding is needed).
3. **Panel integration** — poll, Bets tab, Ticket banner, Edges badge, hide
   filter + alert suppression, freshness, unmatched list, manifest host
   permission. Scripted-Chromium check (below).
4. **Docs + review** — README, CLAUDE.md bullet, pre-merge review, explicit
   approval, merge, worktree + branch cleanup.

#115/#116/#117 start after phase 2 lands (they need the Source protocol and
the fixture layout), each on its own branch.

## Tests

- **Node** (`node --test "unabated_ticket/tests/*.test.js"`): `bets.test.js`
  — every tier on real fixtures (a Kalshi NFL spread against the NFL v2 slice
  event 125807; an MLB doubleheader with `G1`/`G2`; a Kalshi NO-side ML with
  the tie caveat; a sold-to-zero position excluded; an ambiguous-game bet in
  unmatched; native-id dedupe across two consecutive service payloads;
  retention pruning; a record with `league: null` listed under "unknown
  series"). `teams.test.js` — every Unabated NFL/MLB name resolves; unknown
  name ⇒ null.
- **Python** (`pytest unabated_ticket/bets_service/tests`): fills →
  positions aggregation (VWAP, netting, closed), ticker parsing per sport in
  the fixture, both spread signs and both total directions, unknown series
  fails closed, a source failure keeps previous records and writes a failed
  `source_runs` row, `/bets.json` shape.
- **Scripted Chromium** (Playwright in `mlb_sgp/venv`, the #112/#113 harness
  recipe from memory: `--load-extension`, feeds routed to fixtures, fake grid
  page with `__reactFiber$` props, arrow-function-only `evaluate`): route
  `http://127.0.0.1:8094/bets.json` to a fixture; assert the Ticket banner
  text for each tier, the Edges badges, that the hide filter removes exactly
  the `same_line`/`same_side` rows, the freshness colours with a stale
  source, "no source configured" rows, and the unreachable-service message.
  Kept in the scratchpad, not committed (as before).
- **Manual**: service running, a real Kalshi position, click the matching
  line on Unabated, see the banner.

## Version control

- Branch `feature/unabated-bet-history` from `main` (9472414), worktree
  `.worktrees/unabated-bet-history` — both created 2026-09-11 before this
  file. One commit per phase; this plan is its own first commit.
- **#113 (`feature/unabated-alt-edges`) is unmerged and edits `feed.js`,
  `page.js`, `panel.js` heavily.** Phases 0–2 live in new files and cannot
  conflict; phase 3 touches `panel.js`/`panel.html`/`manifest.json` and is
  rebased onto `main` after #113 merges (or #113 onto this — decide before
  phase 3 starts).
- Pre-merge review of `git diff main..HEAD` against the CLAUDE.md checklist
  (data integrity: native-id dedupe and the keep-previous-on-failure rule;
  resource safety: DuckDB connection per pass, closed on error; edge cases:
  empty history, first run, service down, unknown series; dead code; log
  rotation for the service log; secrets: the Kalshi key never appears in
  `/bets.json`, logs, or storage). Findings as ISSUES TO FIX vs ACCEPTABLE
  RISKS, then **explicit user approval before merge**; no push without it.
- After merge: `git worktree remove .worktrees/unabated-bet-history` and
  `git branch -d feature/unabated-bet-history`.

## Documentation (same merge)

- `unabated_ticket/README.md`: a "Bets" section — what is matched and how,
  the four tiers, the service (install, `.env`, `run.sh`, port 8094,
  LaunchAgent optional), the Source protocol and how a venue ticket adds one,
  the unmatched list and how to extend `teams.js`, troubleshooting (service
  unreachable, unknown series, unmatched team).
- `CLAUDE.md`: extend the Unabated Ticket bullet with the bets service, port,
  tier names, and the venue-ticket split.
- This doc gets the phase-0 "Recon findings" section appended.

## Risks and open questions

**Risks**

1. **Kalshi non-MLB ticker grammar** is unverified; the parser fails closed
   rather than guessing a sport, and every such record is visible in the
   Bets tab.
2. **Team-name resolution** is the main source of false negatives (CFB,
   "Athletics", soccer). Every miss is visible in the unmatched list; tables
   grow from real misses. False positives are worse: same team pair on
   different dates (a series) is rejected as ambiguous, not matched.
3. **Wrong-side matches from perspective bugs.** `points` is always the
   side's own number; the Kalshi spread sign is converted once with a test
   per sign. The `opposite` tier is where a sign error would show as a wrong
   red banner — fixtures cover both spread signs and both total directions.
4. **Service unreachable ≠ no bets.** The panel keeps the last payload,
   marks it stale, and says so; a failed source keeps its previous records.
5. **Until #115–#117 land, three venues show "no source configured"** —
   the feature is honest about partial coverage rather than implying a
   clean board.
6. **#113 rebase** (see Version control).

**Open questions for the user**

1. Settled bets: stored, matched only while open. Any objection?
(The Kalshi sports question is answered by the recon section below: NFL +
CFB game markets, no hand-traded MLB.)

## Recon findings — Kalshi (2026-09-11, read-only pull of the account's fills and positions)

Pulled with `kalshi_draft/auth.py` (the helpers `kalshi_mm/analyze_performance.py`
already uses to compute CLV from fills + settlements — that script is the
closest existing code and its `paginate_auth`, fill-cost and P&L helpers are
reused by the service). 196 fills since 2026-07-18, 18 unsettled positions.

**Series traded by hand (taker fills) that the matcher must cover:**

| Series | Fills | Market | Ticker example |
|---|---|---|---|
| `KXNCAAF1HTOTAL` | 37 | CFB 1st-half total | `KXNCAAF1HTOTAL-26SEP05MRSHPSU-29` |
| `KXNCAAFTOTAL` | 26 | CFB total | `KXNCAAFTOTAL-26SEP05LEHGTWN-52` |
| `KXNCAAFSPREAD` | 9 | CFB spread | `KXNCAAFSPREAD-26SEP05MOSUTXAM-TXAM39` |
| `KXNFLGAME` | 5 | NFL moneyline | `KXNFLGAME-26SEP20PITNE-NE` |
| `KXNCAAFGAME` | 2 | CFB moneyline | `KXNCAAFGAME-26SEP05ECUALA-ECU` |

No MLB game markets by hand; `KXMLBRFI` fills are the RFI bot's (maker),
`KXMVECROSSCATEGORY` are the bots' combos. Futures/props (`KXJOINCLUB`,
`KXWCAWARD`, `KXNEXTTEAMNFL`, `KXNFLOROTY`, `KXSTARTINGQBWEEK1`,
`KXWCGOALLEADER`) are recorded as `betType: "other"` with unmatched reason
"not a game market" — shown, never matched. So step 0's grammar question is
answered: **NFL + CFB game markets**, and the MLB parser stays for the bots'
tickers.

**Football ticker grammar** (differs from MLB): `<SERIES>-<YYMMMDD><AWAY><HOME>-<STRIKE>`
— **date only, no HHMM** (MLB carries `HHMM`), so no start time comes from
the suffix; the matcher uses the suffix date (ET) against the Unabated
event's local date, ±1 day guard. Team codes are 2–4 letters (`MOSU`, `TXAM`,
`LEH`, `GTWN`, `PSU`, `NE`) and the away/home split is NOT derivable from the
code block alone — resolve through the public event, which is one cached
GET per event:

- `GET /events/{event_ticker}` → `title: "Missouri St. vs Texas A&M: Spread"`,
  `sub_title: "MOSU vs TXAM (Sep 5)"` — first team is away, second is home
  (same order as MLB). Names are Kalshi's short forms ("Missouri St.",
  "Penn St.", "New England"), so `teams.js` maps Kalshi short names ↔
  Unabated names; codes are never matched on.
- `GET /markets/{ticker}` gives the strike semantics:
  - spread `-TXAM39` → `strike_type: greater`, `floor_strike: 38.5`, title
    "Texas A&M wins by over 38.5 points" ⇒ **YES = TXAM −38.5, NO = MOSU
    +38.5** (integer in the ticker = floor_strike + 0.5; read
    `floor_strike`, never decode the integer);
  - total `-52` → `floor_strike: 51.5`, "Over 51.5 points scored" ⇒ YES =
    Over 51.5, NO = Under 51.5; `KXNCAAF1HTOTAL` identical with `period: "1H"`;
  - moneyline `-NE` → `strike_type: structured`, "New England wins",
    event `mutually_exclusive: true` ⇒ YES = NE ML; NO = PIT **or tie** →
    `approx: ["kalshi_no_side_includes_tie"]` (NFL/CFB overtime makes ties
    rare but the market's NO includes them).
  - `expected_expiration_time` is an end-of-game estimate and `close_time`
    is later still — neither is kickoff (the `close_time` note holds).

**Fill / position payload facts:** fills carry `count_fp` (fractional
contracts, e.g. `12.07`), `side` yes/no, `action` buy/sell,
`yes_price_dollars`/`no_price_dollars`, `is_taker`, `fee_cost`, `ts`,
`created_time`, `exchange_index`. Positions carry **signed** `position_fp`
(`-500.00` on `KXNCAAF1HTOTAL-26SEP12OKLAMICH-23` = 500 NO contracts),
`market_exposure_dollars`, `total_traded_dollars`; `position_fp 0` with
`total_traded_dollars > 0` is a closed position (`status: "closed"`). Stake
= |position| × entry price from the VWAP of that side's buy fills; the
positions endpoint is the truth for the open size, fills only supply price
and time.

**Consequences for the plan above:** `sources/kalshi_ticker.py` becomes a
thin series map (`GAME` → ml, `SPREAD`, `TOTAL`, `1HTOTAL`; `KXMLB*` via
`leg_types`) plus a cached event/market lookup through
`kalshi_draft.auth.public_request` (already throttled to ~1.6 rps) — no
hand-built code table. The Kalshi fixture for phase 0 is one fill and one
market/event payload per series above, anonymised. A full fills pull is a
single page today (196 < 200), so the service does `min_ts` overlap polling
as planned and a full re-pull once an hour as the reconcile.

## Appendix — research: Pikkit versus direct pulls (2026-09-11)

Outcome: **Pikkit dropped** (user confirmed it does not auto-sync). Kept
here so the reasoning is not re-derived.

- **Coverage:** BookSync lists Kalshi, Novig, ProphetX (auto-syncing
  prediction markets) and BetOnline (Pro-only offshore); all four venues in
  scope under one vendor. ([supported list](https://pikkit.com/resources/sportsbooks))
- **Sync trigger — the deciding fact:** "BookSync runs in the background when
  you open the Pikkit app"; "new bets sync when you open the app or take
  certain actions inside it"; backgrounding pauses it. No server-side
  schedule. Pending bets do sync. ([BookSync](https://pikkit.com/booksync),
  [setup guide](https://pikkit.com/blog/how-to-connect-sportsbooks-to-pikkit))
  A bet placed a minute ago is invisible until the phone app runs, which
  defeats the "other side" warning at bet time.
- **Export:** Pro "Bet Transaction Exporting" CSV; columns unpublished.
  ([Pro](https://pikkit.com/pro)) **API:** none. The web app talks JSON to
  `prod-website.pikkit.app` (the repo's `mlb_sgp/pikkit_common.py` already
  intercepts `/betslip`), so a passive content script was feasible, but the
  sync trigger makes it moot. `hoop88_correlation/.pikkit_session.json`
  holds no auth (analytics + UI prefs only).
- **Pikkit's one real advantage** was a single normalised schema (one parser,
  one team table for all books). That saves parsing effort, not
  correctness, and is now spread across #115–#117.
- **Direct routes found** (details in the venue tickets): Kalshi portfolio
  API (have keys); BetOnline `bet_logger` scraper (have code, token likely
  dead since 2026-08-24); ProphetX Trading API `GET /v4/mm/get-order-history`
  / `get-trades` behind an approval (`marketmaking@prophetexchange.com`),
  else content script; Novig NBX API v2 `/emm/orders|fills|positions/all`
  behind "request your client ID and secret from Novig", else content script.
