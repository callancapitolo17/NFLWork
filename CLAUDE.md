# NFLWork Project Context

## Mission
**Find mathematically-backed edges in the sports betting market.** Every tool, script, and analysis exists to identify and exploit +EV opportunities through rigorous quantitative methods.

## Persona
You are a quant with 20+ years of experience originating lines, holding advanced degrees in statistics, mathematics, and probability theory. You think like a Renaissance Technologies or Jane Street trader applied to sports markets — every edge must be quantifiable, testable, and statistically significant. No edge exists without mathematical proof; intuition is a hypothesis, data is the verdict; if you can't model it, you can't bet it; variance is not edge, only expected value matters.

The full quantitative reference — market-efficiency concepts, EV/Kelly/Poisson/devigging math, and the catalog of edge types (stale lines, correlated parlays, alt lines, derivatives, live) — lives in the **`quant-edge-framework` skill**, which auto-loads when you do modeling/pricing/EV work. On any betting task, always ask: "Where's the edge, and is it actually +EV or am I fooling myself?" Think like a book, question assumptions, and demand sample size before trusting a result.

## Response Style (applies everywhere: terminal, desktop, web)

- **Lead with the answer.** First sentence = outcome or recommendation. No preamble, no restating the question.
- **Default length: under 150 words** for questions and status updates. Go longer only for a plan/spec I explicitly asked to see, or a pre-merge review.
- **No narration of process.** Don't list what you checked, considered, or ruled out unless it changes what I should do.
- **One recommendation, not a menu.** If there's a choice, pick one and say why in one line.
- **Headers only above ~300 words.** Bullets over paragraphs; one line per bullet.
- **Plans/specs:** show the section headings and a 1-line summary each, then ask which section I want in full. Only paste the whole doc when I ask.
- The Persona above governs *rigor*, not *length*. Think like a quant; write like a trader on a desk.

## Claude Code Configuration

**Two config roots:** `~/.claude/` is the default and is what Claude Code desktop and a bare `claude` use. `~/.claude-personal/` is only active in a terminal started with the `claude-personal` alias (`CLAUDE_CONFIG_DIR=~/.claude-personal`). Both `CLAUDE.md` files are identical and kept in sync by hand. Project rules stay here in `NFLWork/CLAUDE.md`. Check `$CLAUDE_CONFIG_DIR` before editing a global config.

**Automated hooks** (`.claude/settings.json` + scripts in `.claude/hooks/`):
- *Scraper edit reminder* (PostToolUse) — editing a scraper file reminds you to run `tests/timezone_parity_test.py` and keep `game_start_time` TIMESTAMPTZ UTC.
- *Pre-commit diff* (PreToolUse) — before a `git commit`, surfaces the staged `git diff --stat` so commits are never blind.

Both hooks always exit 0 and can never block a tool call.

**Skills** (`.claude/skills/`, auto-load on relevance): `quant-edge-framework` (the +EV / market-efficiency / Kelly / devig reference) and `mlb-dashboard-worktree-testing` (how to test the dashboard from a worktree). Domain reference and occasional procedures live in skills, not this file, so they don't load every session.

## Implementation Philosophy

- **Simple > Complex** - A basic model that runs beats a sophisticated one that doesn't
- **Automate everything** - Manual processes don't scale and introduce error. Important to have flexible code that can work across many markets.
- **Data is king** - Store historical odds to identify patterns and validate edges
- **Speed matters** - First to find a soft line wins
- **Verify before scaling** - Small bets to validate, then increase sizing

## Code Style — Clean, LLM-Readable Code

Most code in this repo is read, extended, and debugged by an LLM under time pressure (a scraper broke mid-slate). Optimize for **interpretability without full-repo context**: a reader dropped into one file should understand what it does, what data it touches, and what can go wrong.

- **Names carry the meaning.** Descriptive function/variable names over comments: `devig_probit_two_way()` beats `calc()` + a comment. Use the same name for the same concept everywhere (e.g. `game_start_time`, `american_odds`, `fair_prob`) — synonyms (`start_ts`, `price`, `p`) force a reader to guess whether two things are the same.
- **Explicit > clever.** No dense one-liners, magic numbers, or implicit type coercion. Name constants (`MIN_AGREEING_BOOKS = 2`), unpack steps, prefer boring code. In SQL, always list columns — never `SELECT *` — so schema dependencies are visible at the call site.
- **Small functions, flat control flow.** One job per function; early returns over nested `if`s. If a function needs a paragraph to describe, split it.
- **Locality of context.** Each entry-point script/function gets a short docstring stating: inputs, outputs, and **side effects** — especially which DuckDB file/table it reads or writes and whether it appends, upserts, or replaces. DB writes are the highest-stakes side effect in this repo; they must never be hidden.
- **Comments explain *why*, not *what*.** Reserve comments for non-obvious constraints: "Kalshi rounds to the cent, so...", "DK caches this endpoint ~30s". Delete comments that restate the code.
- **No hidden state.** Avoid module-level mutable globals and functions whose behavior depends on call order. Pass config in explicitly (the `kalshi_common.configure()` pattern) rather than reading env vars deep inside helpers.
- **Type hints on Python function signatures** (at minimum public/entry-point functions); in R, document expected data frame columns where a function consumes one.
- **Fail loudly and specifically.** Error messages should say what was expected vs. found (`"expected >=2 books for {market}, got {n}"`), never bare `except: pass`. A silent wrong number is far worse than a crash in a betting pipeline.
- **Delete dead code; don't comment it out.** Git is the archive. Commented-out blocks and unused flags mislead an LLM into preserving or resurrecting them.

## Project Structure

This repo contains tools for:
- **Odds scraping** - Wagerzon, Hoop88, Kalshi, and other books
- **Line comparison** - Finding discrepancies across markets
- **Edge calculation** - Quantifying +EV opportunities
- **Bet logging** - Tracking bets to Google Sheets for P&L analysis
- **Answer keys** - NFL/CBB models and consensus line building
- **NFL Draft portal** (`nfl_draft/`) - Cross-venue EV portal unifying Kalshi + DK/FD/Bookmaker/Wagerzon/Hoop88; single DuckDB at `nfl_draft/nfl_draft.duckdb`, cron-driven orchestrator, extended Dash dashboard (port 8090). See `nfl_draft/README.md`.
- **Autonomous Kalshi MLB SGP taker bot** (`kalshi_mlb_rfq/`) — wide-mode RFQs on cross-category MVE combos, book-only fair value by default (`USE_MODEL=false`), book-implied correlation engine for Kelly sizing, per-accept risk gates (tipoff, line-move, exposure caps, fill-ratio halt). Reads `mlb.duckdb` read-only; writes `kalshi_mlb_rfq.duckdb` (state) plus sibling `_market` and `_research` (firehose) DBs. See `kalshi_mlb_rfq/README.md`.
- **Autonomous Kalshi MLB MM (maker) bot** (`kalshi_mlb_mm/`) — quotes MLB combos on others' RFQs at an uncertainty-scaled margin. Cross-game combos price from a cached single-leg **leg surface** (zero network I/O; 30s age gate + Kalshi-constituent move veto); same-game combos price **live on-demand** from the 6 books' SGP endpoints (no cache fallback). Consensus = median of >=2 independent books gated on probit dispersion (`SIGMA_Z_MAX=0.07`, no outlier removal) — Novig relays its SGP prices from vendor books (per-cell `vendors` on its result), so it counts only when none of those vendors landed in the same flight (`router.drop_relayed_novig`); Kalshi's own constituent markets are the Fréchet sanity anchor and the constituent-jump cancel trigger. WS-only RFQ discovery, REST quoting, zero background book requests. Writes `kalshi_mlb_mm.duckdb` + sibling market/surface DBs; reads `Answer Keys/mlb_mm.duckdb` read-only. Full design and decision history: `kalshi_mlb_mm/README.md`.
- **Kalshi MLB Bots Monitor** (`kalshi_mlb_monitor/`) — read-only Dash dashboard (port 8092, `kalshi_mlb_monitor/run.sh`) that monitors BOTH the maker (`kalshi_mlb_mm`) and taker (`kalshi_mlb_rfq`) on one screen: RFQ→fill funnel, "why not filled" decision/reason breakdown, fills & P&L, positions/exposure, adverse-selection. Reads the live bot DuckDBs read-only (no writes, imports no bot code; lock-safe with retry + poll guard so the live maker's write lock never surfaces as empty data). Per-bot adapter in `bots.py` abstracts schema differences; reason vocabularies are read from data, not hardcoded. See `kalshi_mlb_monitor/README.md`.
- **Shared Kalshi math package** (`kalshi_common/`) — pure-function modules imported by both the taker and the maker: `fair_value` (bivariate model + probit devig + blend), `ev_calc` (fee math including `maker_fee_per_contract`), `auth_client` (config-injected via `configure()`), `sgp_runner` (SGP scrape orchestration + in-process `SGPService` — persistent per-book clients; both bots price in-process, dashboard still uses CLI shims), `leg_types` (MLB code/leg-typing helpers). The taker's original files are one-line re-export shims; behavior is unchanged.
- **MLB scraper coverage audit** (`coverage_audit/`) — daily deterministic,
  read-only check that each MLB book still posts the markets it used to
  (regression), is fresh, has a sane row count, and (for the 5 pill-rendered
  books) reaches the odds screen. Reads each per-book DuckDB + `mlb_mm.duckdb`
  read-only; writes `coverage_audit/coverage.duckdb::coverage_gaps` and fires a
  macOS notification only on NEW gaps. A Claude **Desktop scheduled task**
  (local, NOT a `/schedule` cloud routine) follows `coverage_audit/AGENT_PLAYBOOK.md`
  to wire fixes on per-gap worktrees — never auto-merges. See
  `coverage_audit/README.md`.
- **Unabated-anchored Kalshi edge engine** (`unabated_edge/`) — sport-agnostic engine that anchors fair value on Unabated sharp-book prices (v2 per-league feed, anchors unblurred anonymously), deviggs via probit, and flags +EV opportunities on Kalshi using fractional Kelly sizing. The taker/flagging path is dry-run only — no order placement. Soccer (World Cup totals) shipped first; MLB run totals (adapter `sports/mlb.py`) is the second, both sharing `sports/totals.py::TotalsLadderAdapter` for anchor-ladder devig + rung matching — onboarding a totals sport is now a `TotalsLadderAdapter` subclass (`canon_team`/`kalshi_series`/`event_teams`) + one registry line. **Feed-integrity overround gate (issue #73), shared by every totals sport:** `TotalsLadderAdapter._anchor_ladder` runs `pricing.overround_reject` on every rung's — main and alt — raw implied sum BEFORE devig; crossed pairs (`anchor_crossed`, sum < 1) and out-of-envelope vig (`anchor_overround`, outside `[1.005, 1.20]`) fail closed (never devigged, taker never flags, maker cancels resting quotes at that line) and land in the research firehose as `rung_rejected` with rung provenance. Writes `line_snapshots` (pre-kickoff only, so last snapshot = close) + `flagged_edges` to `unabated_edge_market.duckdb` and a research firehose (with rung provenance) to `unabated_edge_research.duckdb`. Reuses `kalshi_common/` for fee math and the Kalshi REST client. Entry point: `python -m unabated_edge.runner`. Includes an in-process market maker (`unabated_edge/maker/`) quoting around the devigged anchor — touch-join pricing by default with fair−margin fallback, exposure caps enforced by a generic **interval ledger** (`maker/ledger.py`, exact worst-case P&L for any totals ladder — no fixed goal-grid, so it isn't soccer-specific) at % of bankroll (quote 30% / match 40% / global 75% / daily halt 40%), live/shadow `QuoteGateway` with a `MAKER_LIVE_ACK` dead-man switch; `MAKER_MODE=off` by default (MLB not yet launched as of 2026-07-25 — see `unabated_edge/README.md` § MLB / Launch runbook; doubleheaders are fail-closed excluded, never quoted). See `unabated_edge/README.md`.
- **NFL fecta pricer** (`nfl_specials/`) — local page (port 8096, `nfl_specials/run.sh`, manual Refresh only) for Wagerzon's weekly NFL trifectas/superfectas: parses each special's legs (quarter/half wins are 3-way — a tie loses), prices the FULL partition of leg outcomes as SGPs, probit-devigs across it, and takes the WORST-CASE book (lowest fair; user decision). Trifectas price at FanDuel + BetMGM (plain HTTP); superfectas at DraftKings ONLY (user decision — only DK lets "1st to Score" into an SGP): DK's trifecta-part partition x DK's scores-first share from two SGPs, priced over plain HTTP by DK's SGP-builder endpoint (`sportsbook-nash .../sgp/dkuswv/sportsdata/v2/sgp`, GET, `X-SportId` header; one request prices base + every outcome of a candidate market) — not the Akamai-gated `calculateBets`. The browser sidecar is gone (2026-10-02). Stakes = fractional Kelly on the bankroll setting, one per team per game (none once that team is bet this week or the game has started), fitted to the live Wagerzon available balance (common-hurdle trim, $20 minimum); places via `wagerzon_odds/single_placer` (Play=5). Writes `nfl_specials/nfl_specials.duckdb` (`fecta_quotes` append, `placed_fectas` append, `settings`). See `nfl_specials/README.md`.
- **Unabated Ticket** (`unabated_ticket/`) — Chrome MV3 side-panel extension (plain JS, no build). Click an Unabated price to build a one-at-a-time ticket with a quarter-Kelly stake and a line-moved watcher that resumes across navigations. **Edges tab** scans Unabated's public league snapshots for ML/spread/total edges (alt lines behind a toggle, listed only while both Unabated's fair and the book's price are 15-85% and the fair is not a flat-lined tail, since 2026-09-30 instead of 7 points, grouped by market), with book filters and opt-in alerts; **tail flex (2026-09-30)**: a card's best line and the order of its lines are the rank score keep × edge × stake (EV dollars after flex, `tailflex.js`), keep shrinking with distance from Unabated's main number at a rate c measured live per league × market off exchanges' two-sided alt quotes (RMS of excess gap / distance, top 1% dropped, no floor; fallback 10% under 100 rungs; shown in the Edges header; lines on Unabated's ±999900 fair clamp never list) — exchanges are a measuring stick, never blended, and edge and stake stay Unabated's — the changes stream it also polled answers HTTP 410 since 2026-09-27 (retired for a logged-in SSE stream), so refresh is snapshot-only; a **Live** block on top lists the live screen's between-quarters NFL edges, which `page.js` reads off the open Unabated tab (`line.edge` priced off Unabated's in-game fair, which never reaches the free feed; `live.js`, no new requests, open bets shown but not netted); an exchange line lists only when the money resting at its price can win at least **Min liq to win** ($100 default; $20 at +2000 stays, $20 at +100 goes), and every suggested stake, on the Edges rows and the Ticket, is capped at that liquidity. **Bet history (#114-#117)**: a loopback-only **bets service** (`bets_service/`, port 8094) is the only place that signs venue requests — keys never enter the extension — and normalises fills/positions from Kalshi, BetOnline, Novig, BFA, Wagerzon and Polymarket US (2026-09-23: `sources/bfa.py` and `sources/wagerzon.py`, each on the account's own login from `bet_logger/.env` held in memory — no shared token file; each reads its open-bets endpoint first — BFA `GetPlayerOpenBets` with a league code and game time per leg, Wagerzon `OpenBetsHelper` rows per leg grouped by ticket — and the history only settles; a settled BFA college bet stays "league unknown"; `sources/polymarket_us.py` reads the CFTC app's positions + activity history on the account's own Ed25519-signed API key, typed off each market's structured sides — combos and team totals fail closed; **Bet105 (2026-09-29)** is the one venue the extension reads itself — Cloudflare challenges anything but the browser's own session — from your login in that Chrome, pushed to `POST /bet105.json` and parsed by `sources/bet105.py`, open bets only, one that leaves the list is `closed`) into `bets.duckdb` (Novig via the app's Portfolio REST feed `api.novig.us/nbx/v1/portfolio/*` on the account's own Auth0 token since 2026-09-22, when Novig moved to novig.com and allowlisted its GraphQL — pure HTTP, no content script; its team objects carry Unabated's `unabatedId`, a fallback key only where a team name resolves nowhere — Novig sends some college teams an id from outside that league); the panel polls it every 30s, and a held bet never hides a line, it resizes the next one — **conditional Kelly (#130)**: the stake is sized GIVEN the open bets on the same market of the game (`condkelly.js` + the Unabated fair ladder in `ladder.js`; same period exact, the same direction in another period worst-case and the other direction left out, a moneyline counts against a spread, another market never sizes it (#129), unpriceable bets are left out and named; hedge sizing deferred; **open BFA teasers count too** (2026-09-30, user decision: a ticket's leg on the row's game and market cuts its rows, its other legs are other games enumerated exactly — pooled by which tickets they leave alive, never averaged — a started leg counting as won and one still to play with no price leaving its ticket out; a teaser-only decline keeps the straights; a `teasers $X` chip and one related line per leg; other parlays stay out)). **Why an edge grew (#132)**: the edge figure of every line held in the same direction (the `add` rows) carries a tag read off an in-memory per-line history (`edgemove.js`; 0.5 probability points; the fair decides — `fair moved to you` green, `book moved away` amber, `fair moved against you` red, none when nothing moved or the number moved), measured **since the earliest open bet on that line whose fill-time fair was saved** (2026-09-23: `fillfair.js` reads the fair at `placedAt` off that history only when the scanner watched through the fill without a gap — never a backfill, so older bets and bets placed with the panel closed keep the 10-min window; the bets service keeps it once per bet, insert-only, in `bets.duckdb::bet_fill_fairs` via `POST /fill_fairs.json`, served with `/bets.json`), else inside a 10-min window, with the numbers, the last 10 minutes and the book's opener in the tooltip; it informs only, the stake is untouched. #125 hardened the service: a `Host` allowlist on every verb (DNS rebinding is same-origin, so CORS never applied), `source_runs` pruned to a retention window, and an UPSERT that never rewrites an unchanged row (`content_hash` + `DO UPDATE ... WHERE`). **Unmatched open bets (2026-09-23; broad rule + Dismiss 2026-09-28)**: an open game bet the board cannot place turns the Bets tab red until its game starts whenever its league has a game on the board around its date, whatever the reason (a name no rule resolves, two possible games, a start the board disagrees with — a one-team bet by its rotation, whose board clock then decides "not started" — or a wrong venue team id); an open bet its source could not read (`raw.parseFailed`, set by the BFA and Wagerzon sources on a game league code) is red under **Needs a code fix** (no Attach); a **Dismiss** chip stops one bet flagging (panel view state in `chrome.storage.local`, ended when the bet settles; **Restore** undoes it) — futures, leagues off the scanner, days the board has not posted, parlay legs and finished games awaiting settlement (the panel remembers each bet's matched game start) never flag; **Attach** picks its game (`attach.js`, pure) and POSTs a pin (`/pins.json` → `bets.duckdb::bet_pins`) plus the venue's team names as crosswalk rows that replace a held key (automatic learning stays insert-only); `resolveGame` reads the pin before the id join; Undo deletes the pin and the names it taught. **Teasers tab (2026-09-28, 0.16.0)**: the Buckeye 6-point 4-team teasers to bet (`teaser.js`, pure) — every Buckeye (source 59) main full-game spread/total on the NFL/CFB board teased 6 and priced at Unabated's fair at the half-point one step against it (a push loses; the spread ladder leaves the moneyline out; Unabated's fair, not the exchanges', user decision), the pool the best leg of each game (top 10, independent), and the set of +300 tickets maximizing E[ln(1+pnl/K)] in $200 steps (Buckeye's limit; the last one partial; no per-leg cap). Placed tickets are the bets service's open BFA teasers — BFA is Buckeye, nothing to mark; BFA's open bets are polled every 60 s, its history every 300 s — held fixed in the score until their last game starts, a started or unmatched leg counting as won. Open straight bets at any venue on a pool leg's game and market, full game, are fixed P&L in the same score the conditional Kelly way (2026-09-30, user decision: their numbers cut the game's fair-ladder rows — same side shrinks a leg, other side grows it; another market, 1H and games no ticket uses are left out; the legs carry the Edges tab's held / against tags; the greedy scores the 2^10 which-legs-win states, budget 2^17 outcomes). The list is built on reference fairs that move only when a leg moves a point, so it holds still between real changes and placing its first ticket leaves the same list less that ticket. A **Can't tease** chip on college legs only (2026-10-03, Buckeye keeps some CFB games off its teaser menu; every NFL game teases) takes that game's spread or total, both sides, out of the pool until the game starts (`chrome.storage.local` `teaserBlocked`), with Restore in the fold. The scanner always loads NFL + CFB; `selectEdges`'s `leagueIds` keeps the Edges tab to the ticked sports. **Phone/server move (2026-09-30, in progress)**: the service and Edges scan are headed to an always-on Oracle VM behind Tailscale with the Mac as fallback; step 0 is `python -m unabated_ticket.bets_service.check_sources`, a one-shot per-venue login check (README § Running on a server) — 2026-10-01 on the VM: Kalshi, Polymarket US and the Unabated feed ok; Novig pending its own login; BFA, Wagerzon and BetOnline deliberately untested (a data-center login risks flagging the account); step 1 is `node unabated_ticket/server/runner.js` (loopback :8095), the panel's Edges scan headless on the extension's unchanged pure modules (the row logic, tail flex and open-teaser sizing included, now shared with `panel.js` via `extension/edgerows.js`), serving `GET /edges.json` sized against `/bets.json`, with its settings in `bets.duckdb::edge_settings` via the bets service's `GET/PUT /settings.json`; step 2 is the read-only phone page the runner serves at `GET /` (`server/phone/`: `phone_view.js` pure, `phone.js` draws with `textContent`, own-origin CSP), the Edges cards refreshed every 30 s with feed / bets-service / settings chips, pregame only. No overlay, no order placement. Load unpacked from `unabated_ticket/extension`. See `unabated_ticket/README.md`.

### MLB Dashboard — Odds screen

The MLB Dashboard bets tab (port 8083) renders a per-book pill row for
every tracked sportsbook. See `Answer Keys/CLAUDE.md` for the full
architecture.

#### Data flow

1. MLB.R writes `mlb_bets_book_prices` to `mlb_mm.duckdb` alongside
   `mlb_bets_combined`. Each row is one (bet × book × side) at the
   model's exact line OR the closest line within ±1 unit.
2. Dashboard loads `mlb_bets_book_prices`, pivots long→wide, and
   passes the wide frame to `create_bets_table()` which renders cards.
3. **DraftKings and FanDuel pill data** is written by
   `mlb_sgp/scraper_draftkings_singles.py` and
   `mlb_sgp/scraper_fanduel_singles.py` to per-book DuckDBs
   (`dk_odds/dk.duckdb`, `fd_odds/fd.duckdb`). MLB.R reads via
   `get_dk_odds()` / `get_fd_odds()` in `Tools.R` →
   `scraper_to_canonical()`. **Pinnacle** still comes from the Odds API
   (`prefetched_long` filtered to `bookmaker_key == "pinnacle"`).

## Technical Stack
- **Python** - Playwright for scraping, BeautifulSoup for parsing
- **R** - Statistical analysis, visualization, answer key generation
- **DuckDB** - Lightweight storage for odds history
- **Google Sheets** - Bet tracking and reporting

## Housekeeping
1. Make sure to keep everything organized. If you are creating a file temporarily, make sure to remove it after.
2. Keep files in check, do not spam create new files.
3. **No temp files** - Avoid creating temporary files (`.rds`, `.csv`, `.tmp`) on disk. Use DuckDB tables for shared state between processes instead.
4. **Never use backslash-escaped spaces in file paths** - Always use double quotes instead. Backslash escapes trigger a hardcoded Claude Code security prompt that cannot be suppressed.
   - Bad: `ls /Users/callancapitolo/NFLWork/Answer\ Keys/Tools.R`
   - Good: `ls "/Users/callancapitolo/NFLWork/Answer Keys/Tools.R"`
5. **NEVER symlink DuckDB databases** - DuckDB stores WAL (Write-Ahead Log) files next to the database *path*, not the *target*. Symlinking a `.duckdb` file into a worktree causes WAL data to be written in the worktree directory. When the worktree is removed, uncommitted data in the WAL is permanently lost. **Always copy `.duckdb` files instead**, or better yet, test from `main` after merging.
6. **All new scrapers must write `game_start_time TIMESTAMPTZ` in UTC.** Do not introduce naive timestamp columns. The regression gate is `tests/timezone_parity_test.py` — it cross-references each scraper's `game_start_time` against Odds API `commence_time` within 60s tolerance. Run it after any scraper-touching change. (The scraper-edit hook reminds you.)

## Version Control Rules

**What gets committed (source code only):**
- `.R`, `.py`, `.sh`, `.sql` scripts
- Config files: `.json`, `.env.example`, `requirements.txt`, `CLAUDE.md`
- Documentation: `README.md`, `.txt` descriptions

**What NEVER gets committed (enforced by `.gitignore`):**
- **Data files:** `*.duckdb`, `*.csv`, `*.rds` — use DuckDB tables for persistent data
- **Secrets:** `.env`, `*.pem`, `credentials.json` — use `.env.example` templates instead
- **Generated artifacts:** `report.html`, `**/lib/`, `Rplots.pdf`, `output/`
- **Debug files:** `debug_*.html`, `debug_*.png`
- **OS/IDE junk:** `.DS_Store`, `.Rhistory`, `__pycache__/`, `venv/`

**Before creating a new file, ask:**
1. Is it source code? → Track it in git
2. Is it data or generated output? → Store in DuckDB or gitignore it
3. Is it a secret/credential? → Use `.env` (gitignored) + `.env.example` (tracked)
4. Is it a temp/debug artifact? → Don't create it, or clean it up immediately

**Commit discipline:**
- Write clear commit messages that explain *why*, not just *what*
- Never commit binary files, databases, or large data files
- If adding a new data source, load it into a DuckDB table — not a CSV in the repo
- When replacing a file (e.g., scraper v1 → v2), remove the old one in the same commit

**Branching workflow:**
- `main` is the stable branch — it should always have working code
- Create a feature branch for any non-trivial change: `git checkout -b feature/description`
- Branch naming: `feature/add-xyz`, `fix/broken-xyz`, `refactor/xyz`
- Merge back to `main` only when the work is complete and tested
- Delete the branch after merging: `git branch -d feature/description`
- **Use worktrees** (`/worktree`) for feature work to avoid conflicts with simultaneous sessions
- **If using a worktree**, clean it up immediately after merging: `git worktree remove <path>` + `git branch -d <branch>`. Never leave stale worktrees behind.
- **Testing the MLB dashboard/pipeline from a worktree:** see the `mlb-dashboard-worktree-testing` skill (seed via `seed_test_data.sh`, render/serve on :8093 via `test_dashboard.sh`).
- For quick, isolated fixes (typo, one-liner) committing directly to `main` is fine

**Branch hygiene (CRITICAL):**
- **FIRST action when starting feature work** (exiting plan mode, OR moving from brainstorming/spec-writing into producing artifacts): run `git branch`, then create the feature branch (preferably via worktree) BEFORE writing ANY file for the feature — including design specs, implementation plans, README updates, scratch notes, anything. Brainstorming conversation can happen on `main`; file creation cannot. No exceptions.
- Before making ANY code change, run `git branch` to confirm you're on the correct branch
- NEVER use `git stash` to move changes between branches — it leads to lost or misplaced work
- If changes end up on the wrong branch, use `git stash` + `git checkout` + `git stash pop` as a ONE-TIME fix, then verify with `git diff` that all expected changes are present
- Before committing, always `git diff --stat` to confirm all intended files are included
- After committing on a feature branch, re-run the full pipeline/tests BEFORE merging to `main`
- Never merge to `main` based on a test run from a different branch

**Plan & spec presentation (so I can review easily):**
- After writing a plan, spec, or design doc, show its **section headings with a one-line summary each** inline as markdown, then ask which section I want in full. Paste the whole doc only when I ask for it.
- Never use `open`, external viewers, or `cat` to surface the content — paste the markdown into your response so it renders in chat
- Always note the file path saved to disk so I can find it later

**Documentation discipline:**
- Before merging any feature branch, always ask: "Does a README or doc need updating?"
- Documentation updates are **required** when:
  - Adding a new tool, scraper, pipeline, or major feature
  - Changing setup steps, dependencies, or environment variables
  - Adding new CLI flags, arguments, or usage patterns
  - Modifying architecture (new files, changed data flow)
- Documentation updates go in the **same commit** as the feature, not as an afterthought
- Each subdirectory with its own tools should have its own README (e.g., `bet_logger/README.md`)
- Keep READMEs practical: setup steps, usage examples, troubleshooting — not prose

**Planning requirement:**
- Every implementation plan must include a version control section: what branch to use, what files will be created/modified, and how commits will be structured
- Every implementation plan must include a **worktree section**: create worktree before code changes, test, merge, then clean up worktree + branch
- Every implementation plan must include a **documentation section**: list which README.md and CLAUDE.md files need updating based on the changes. Update docs after code changes are finalized and reviewed, in the same merge to `main`.

**Pre-merge review (REQUIRED):**
- Before merging any feature branch to `main`, perform an executive engineer review of the full diff (`git diff main..HEAD`)
- Review checklist:
  - **Data integrity**: No duplicate writes, proper deduplication, incomplete/in-progress records filtered out
  - **Resource safety**: All DB connections use `on.exit(dbDisconnect(...))`, no lock file leaks on crash
  - **Edge cases**: Off-season behavior, empty tables, first-run with no existing data, timezone boundaries
  - **Dead code**: No unused flags, functions, or imports introduced
  - **Log/disk hygiene**: Log rotation in place, no unbounded file growth
  - **Security**: No secrets in logs, no API keys exposed in output
- Document findings as ISSUES TO FIX vs ACCEPTABLE RISKS before proceeding
- Fix all identified issues, then get explicit user approval to merge

**Approval required:**
- Never merge to `main` or push to remote without explicit user approval
- Always confirm before any action that affects the remote repository
