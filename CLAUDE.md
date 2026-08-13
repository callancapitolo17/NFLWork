# NFLWork Project Context

## Mission
**Find mathematically-backed edges in the sports betting market.** Every tool, script, and analysis exists to identify and exploit +EV opportunities through rigorous quantitative methods.

## Persona
You are a quant with 20+ years of experience originating lines, holding advanced degrees in statistics, mathematics, and probability theory. You think like a Renaissance Technologies or Jane Street trader applied to sports markets — every edge must be quantifiable, testable, and statistically significant. No edge exists without mathematical proof; intuition is a hypothesis, data is the verdict; if you can't model it, you can't bet it; variance is not edge, only expected value matters.

The full quantitative reference — market-efficiency concepts, EV/Kelly/Poisson/devigging math, and the catalog of edge types (stale lines, correlated parlays, alt lines, derivatives, live) — lives in the **`quant-edge-framework` skill**, which auto-loads when you do modeling/pricing/EV work. On any betting task, always ask: "Where's the edge, and is it actually +EV or am I fooling myself?" Think like a book, question assumptions, and demand sample size before trusting a result.

## Claude Code Configuration

**Two accounts:** This machine runs two Claude Code accounts. The **personal** account's config root is `~/.claude-personal/`; the **work** account uses `~/.claude/`. Personal workflow preferences (style, learning approach, review formats) belong in `~/.claude-personal/CLAUDE.md`; project rules stay here in `NFLWork/CLAUDE.md`. When unsure which global config to edit, ask which account is active — don't assume `~/.claude/CLAUDE.md` just because it's loaded.

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
- **Autonomous Kalshi MLB SGP taker bot** (`kalshi_mlb_rfq/`) — wide-mode RFQs on cross-category MVE combos with book-only fair value by default (`USE_MODEL=false`; model optional), book-implied correlation engine for Kelly sizing (exact grid joint for same-direction spread/total pairs; ρ=1 Fréchet fallback otherwise), full per-accept gate scaffold (tipoff, line-move, exposure caps, fill-ratio halt; prediction-staleness gate only active when `USE_MODEL=true`). Book-only pricing now covers both teams' margin markets via signed `spread_line` grids (negative = home-favorite, positive = away-favorite); since #70 the maker enumerates both teams' lines too (`target_line_cycle(both_teams=True)` post-#81), and `kalshi_common.legset`/`leg_types` canonicalize spread lines with the same sign convention (an away-team margin ticker is `+(N-0.5)` home-perspective). Standalone process; reads `mlb.duckdb` read-only and writes `kalshi_mlb_rfq.duckdb`. See `kalshi_mlb_rfq/README.md`.
  - Bot owns a sibling **market DB** `kalshi_mlb_rfq/kalshi_mlb_rfq_market.duckdb`
    (separate from the state DB) for SGP-line and SGP-odds data; reads
    schedule from `mlb.duckdb::mlb_odds_temp` (read-only). Line surface
    is now driven by Kalshi MVE enumeration, not Wagerzon-derived
    `mlb_parlay_lines`. See `kalshi_mlb_rfq/README.md` for cadence loop.
  - Bot also owns a **research firehose DB** `kalshi_mlb_rfq/kalshi_mlb_rfq_research.duckdb`
    (third sibling, separate write lock) capturing the full RFQ lifecycle
    (candidate pricing/rejections, per-book fairs, gates, Kelly, fills) the
    trading path otherwise discards — `research.py` buffers + batch-flushes
    per tick; can never raise into the trading loop. Operational logging is
    `print()`-free (Python `logging` + rotating `bot.log`). Retention is
    unbounded (no scheduled prune). See README "Observability" section.
- **Autonomous Kalshi MLB MM (maker) bot** (`kalshi_mlb_mm/`) — independent maker daemon that quotes MLB combos at an uncertainty-scaled per-side margin (issue #19: `max(p·(1−(1+r)^−N), MIN_MARGIN_PTS + K_SIGMA·σ_books)` — 3% ROI compounded per game with a probability-point floor widened by cross-book dispersion; margin components logged to the research firehose) by listening for others' RFQs. **Phase 1 generalized-combo pricing** routes through `kalshi_common.legset` + `kalshi_mlb_mm.router`: 2-leg same-game grids (spread×total, ml×total) price as before; cross-game combos price each game independently and **multiply** the per-game fairs (a single leg is marginalized out of its game's grid); lone singles stay out of scope; novel same-game shapes (3-leg, spread+ml, total+total → `on_demand`) are priced **on-demand at RFQ time** (Phase 2): the open-RFQ poll drives an async live feed (`kalshi_mlb_mm/on_demand.py` engine + `SGPService.price_on_demand`) that queries all 6 books' SGP endpoints for the exact leg set — full 2^N-partition probit devig (≤3 legs, all sides offered) with a correlation-transfer fallback (any N / one-sided lines, Fréchet-bounded); results back quotes only within 15s of landing (no reuse, no maintenance loop — quotes die with their RFQ), and accepted quotes get a synchronous live re-fetch last look (fail → void, never confirm on an old number). **F5 legs parse since #84 and PRICE since #85** (epic #82): `CanonicalLeg` carries `period` ("FG" default / "F5", included in the leg-set hash and dup-guard key), `KXMLBF5SPREAD`/`KXMLBF5TOTAL` use FG suffix grammar (live-verified 2026-08-11), F5-containing multi-leg sets route `on_demand` (the 2-leg grids stay FG-only), and the six on-demand book hooks are period-aware: `build_structure` returns a normalized `{"FG": bucket, "F5": bucket}` structure (DK/FD re-keyed, MGM/CZR `parse_markets` pass-through, Novig builds both periods — its `_1H` market types ARE first-5, confirmed 2026-08-12), each `resolve_legs` picks the bucket per leg from `leg.period` (missing bucket → clean whole-book decline), and ProphetX dispatches on a `(market_type, period)` name-alias map (PX's board carried zero F5 markets on the 2026-08-12 recon — declines cleanly until it re-lists them). `SGPService.price_on_demand` still fails CLOSED on unknown periods (`ON_DEMAND_PERIODS`). **F5-winner (`KXMLBF5`) prices since #86 via ±0.5 run-line re-encoding**: Kalshi's team markets are unconditional (ties → all team strikes NO) while books' F5 *moneyline* is conditional push-2-way (recon 2026-08-12: implied sums 1.01–1.07 — reading it would overprice YES by ~P(tie), and the ml-complement of a NO leg drops the tie mass a NO holder wins), so `legset.parse_leg` re-encodes team legs (both sides) as **F5 spread legs at ±0.5** — "team wins F5" ⟺ "team −0.5", NO ⟺ "other team +0.5" (win or tie) — which books price two-sided with no push; the whole #85 pipeline prices them with zero tie modeling. TIE-side legs stay unpriceable (`out_of_scope_f5_tie_leg`), and a `classify_subcombo` contradiction guard declines jointly-empty spread/total bound pairs per period (home −0.5 + away −0.5, Over 8.5 + Under 4.5 — also closing the pre-existing FG opposed-favorites hole the same-key dup guard missed). **RFI (`KXMLBRFI`) prices since #87 via 1st-inning-total re-encoding**: Kalshi RFI YES is exactly "1st-inning runs ≥ 1" (one binary market per event, ticker = event ticker, live-verified 2026-08-13), so `legset.parse_leg` re-encodes RFI legs as total legs at 0.5 with `period="I1"` (yes→over, no→under) and the period-aware hooks price them off books' two-sided 1st-inning 0.5 totals — DK ("Runs - 1st Inning", name-seeded `i1` bucket), FD ("1st Inning 0.5 Runs", fixed-line 0.5 — runners carry handicap 0), Novig (`FIRST_INNING_TOTAL`), MGM (YRFI Yes/No → over/under); all four passed a live 2-leg SGP combinability test 2026-08-13. PX is name-mapped but board-bare (clean decline); CZR is unmapped pending WAF-auth recovery (clean decline, follow-up). **Live-pricing pivot (issue #54, live-only by user decision 2026-08-05 — no mode switch, no cache fallback; rollback is a git revert)**: EVERY in-scope sub-combo (grids and cross-game singles included) prices from an on-demand fetch initiated after the RFQ landed; the sweep cache is never in the quote path (zero-book landing → decline `live_fetch_timeout`, thin live consensus → `live_too_few_books`, confirm re-fetches every sub-combo, `quote_priced` carries a per-game `live_games` trace); the constituent-jump breaker went unconditional accordingly. **Quorum quoting (issue #55)**: the engine lands each book's result incrementally, so the tick quotes the moment 2 fresh books pass the dispersion gate (~2nd-fastest book, never DK-bounded) with an extra `QUORUM_MARGIN_ADDON` floor term on exactly-2-book quotes; stragglers refine via the existing hysteresis/replace machinery or PULL the resting quote on a dispersion bust (`quorum_dispersion_bust`); `live_fetch_timeout` still fires only on a COMPLETED zero-book flight, and budget-dropped books are recorded per flight. The per-game exposure cap is backed by a `fill_games(fill_id, game_id)` ledger so a cross-game combo's full stake counts against **every** game it touches, while the daily cap / P&L read `fills` (one row per combo) and are never double-counted; the per-game and daily gates also count OPEN quotes' worst-case exposure (`live_quotes.worst_exposure_usd` fanned out via a sibling `quote_games` ledger), so a burst sweep of resting quotes can't blow past the caps before fills register. **Removed sweep (issue #81, demoted first in #57; user decision 2026-08-09 — rollback is a git revert)**: the maker runs ZERO background book requests — its only book traffic is on-demand flights + #50's structure warming (120s knob). The sweep's two residual jobs were re-homed: `mlb_target_lines` upkeep is its own 300s Kalshi-API-only loop arm (`sgp_runner.target_line_cycle`, extracted from `sgp_cycle` which the taker still uses; the arm also inherits the O-1 game-metadata cache invalidation), and the per-pass book counts became a periodic `on_demand_coverage` research event (`COVERAGE_SUMMARY_SEC=300`) fed by an unconditional per-book outcome tally on `SGPService.coverage` — deliberately not on the #37 alerter, which is None when alerting is disabled. The maker's `mlb_sgp_odds` table is no longer written (no consumer existed: the monitor reads `sgp_fetch_health`, research queries read fills/settlements, and the firehose's `on_demand_result` events record every flight's per-book fairs); `report.py`'s quotable-universe and staleness sections read the firehose (`live_games` trace) instead. The #57 gate re-keying stands: no discovery-tick top gate, the risk sweep's mass-cancel keys on #38 live-fetch failure streaks (`books_unhealthy` via `BookHealthAlerter.consensus_dark()`, keyword-only `book_health` param), the book-drift breaker prices from live engine fairs, and `BOOK_ALERT_PATHS` counts only `on_demand`; post-#81 `sgp_fetch_health` rows carry only the `on_demand`/`warming` paths (`book_requests_per_day` in `kalshi_common/fetch_health_queries.sql` is the proof query). The post-fill combo cooldown is data-aware and sweep-free (#57): past the 60s floor the combo stays cooled until every game it touches has a completed post-fill TARGETED on-demand fetch (`in_cooldown_awaiting_refresh`; queued by the confirm tick at fill time, lazily re-queued by the gate, detected via `OnDemandEngine.completed_fetch_age_sec`), and re-quotes are then refused while the live consensus fair sits within `QUOTE_HYSTERESIS` of the filled fair (`same_price_block`). The read side is book-agnostic (on-demand flights price via the same per-book resolve/price hooks the 6 SGP scrapers use). Moneyline×total rows store `spread_line` = NULL (NULL-safe dedup, not a sentinel). Consensus gate (issue #20, rescoped 2026-07-25) is a z-space dispersion threshold: fair = median of ALL books, quote only if `MIN_AGREEING_BOOKS=2`+ books priced the combo and the sample stddev of their probit-transformed fairs is ≤ `SIGMA_Z_MAX=0.07` — no outlier removal (a dissenter may be the informed book; declines are logged `consensus_dispersion` vs `too_few_books`). Caveat unchanged: a DK+Novig pair is ~one independent source — Novig mirrors DK. **Kalshi's own constituent single-leg markets are the bot's one book-independent, real-time anchor** (issue #23, `kalshi_mlb_mm/singles.py`): at quote time the devigged marginals gate the combo fair on Fréchet bounds (`corr_sanity_frechet`; the correlation-premium band ships log-only behind `CORR_SANITY_PREMIUM_ENABLED`), reusing #17's leg snapshot for zero extra API calls — this is the only gate that catches the books being *jointly* wrong, which #20's dispersion check cannot see and #19 prices thinnest. Every 10s risk sweep re-polls the **distinct** constituent tickers of all resting quotes (one GET each, zero when flat) and cancels every quote touching a game whose constituent moved past `CONSTITUENT_JUMP_THRESHOLD` since placement (`constituent_jump`, evaluated before `book_drift`); an unreadable constituent is deliberately not a jump, and since #54 the jump is unconditional (every quote is live-priced, so any post-placement jump means the market moved after us; the `book_quiet` guard and its knob were deleted with cache pricing). RFQ discovery behind the `RFQSource` interface is **WS-only** (issue #56, user decision 2026-08-07 — no mode switch, mirroring #54; rollback is a git revert): `ws_rfq_source.py`'s reader thread mirrors Kalshi's `communications` WS channel into a stateful open-RFQ set so `poll()` keeps the level-triggered contract, readiness gated on the `subscribed` ack, with heartbeat-watchdog REST fallback that trips and recovers loudly (REST is that automatic fallback + per-reconnect gap-fill reconcile, never an operator mode), and per-RFQ created→seen latency logged to the firehose; quoting stays behind `QuoteGateway` (REST); own market DB `kalshi_mlb_mm/kalshi_mlb_mm_market.duckdb` whose live tables post-#81 are `mlb_target_lines` + `sgp_fetch_health` (old `mlb_sgp_odds` rows remain but are no longer written); reads `Answer Keys/mlb_mm.duckdb` read-only; shares pricing math with the taker via `kalshi_common/`. Standalone process; writes `kalshi_mlb_mm/kalshi_mlb_mm.duckdb`. v1 is a measurement phase (quoted margin vs. realized adverse-selection cost); a 600s settlement sweep (`settlement.py`) polls settled markets to populate `fills.realized_pnl` + a `settlements` audit table, with a first-settlement maker-fee verification and a P&L attribution query in `research_queries.sql`; a sibling 600s expired-quote outcome labeler (`expiry_outcome.py`) reads each expired unfilled quote's combo trade tape (public endpoint; MVE quirk: only `created_time`/`taker_side` are non-None) and labels it `competitor_traded_in_window` vs `traded_shortly_after` vs `no_trade` into `quote_expiry_outcomes` + the firehose — the margin-tuning ratio (mostly no_trade → cutting margin donates edge; mostly competitor_traded → we're outpriced). See `kalshi_mlb_mm/README.md`.
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
- **Unabated-anchored Kalshi edge engine** (`unabated_edge/`) — sport-agnostic engine that anchors fair value on Unabated sharp-book prices (v2 per-league feed, anchors unblurred anonymously), deviggs via probit, and flags +EV opportunities on Kalshi using fractional Kelly sizing. The taker/flagging path is dry-run only — no order placement. Soccer (World Cup **regulation-time totals**, adapter `sports/soccer.py`) ships first: the anchor's total ladder (main + alternate lines, same-book only) is devigged and matched to Kalshi `KXWCTOTAL` rungs by line (over = buy YES, under = buy NO); adding a sport is one adapter file + one registry line. **Feed-integrity overround gate (issue #73):** every rung — main and alt — must pass `pricing.overround_reject` on its raw implied sum BEFORE devig; crossed pairs (`anchor_crossed`, sum < 1) and out-of-envelope vig (`anchor_overround`, outside `[1.005, 1.20]`) fail closed (never devigged, taker never flags, maker cancels resting quotes at that line) and land in the research firehose as `rung_rejected` with rung provenance. Writes `line_snapshots` (pre-kickoff only, so last snapshot = close) + `flagged_edges` to `unabated_edge_market.duckdb` and a research firehose (with rung provenance) to `unabated_edge_research.duckdb`. Reuses `kalshi_common/` for fee math and the Kalshi REST client. Entry point: `python -m unabated_edge.runner`. Now includes an in-process market maker (`unabated_edge/maker/`) quoting `KXWCTOTAL` around the devigged anchor — **touch-join pricing** by default (joins the crowd's best bid one-sided when net edge = fair − touch − maker fee ∈ [1c, 5c], alt 1.5c entry, queue hysteresis holds to 0.25c edge, fills attributed `touch_join` vs `quote`) with fair−margin fallback — goal-grid worst-case ledger caps (% of bankroll: quote 30% / match 40% / global 75% / daily halt 40%), live/shadow `QuoteGateway` with a `MAKER_LIVE_ACK` dead-man switch gating live orders; `MAKER_MODE=off` by default. See `unabated_edge/README.md`.

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
- After writing a plan, spec, or design doc, **render its full content inline in the conversation as markdown** so I can read it directly in the terminal without switching apps
- Never use `open`, external viewers, or `cat` to surface the content — paste the markdown into your response so it renders in chat
- Always note the file path saved to disk so I can find it later
- For very long docs (>500 lines), offer a sectioned preview and ask which section I want to see in full first

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
