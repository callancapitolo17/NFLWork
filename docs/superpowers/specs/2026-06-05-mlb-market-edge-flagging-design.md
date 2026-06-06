# MLB Dashboard — Dual-Source Edge Flagging (Model OR Market)

## Review Pack

**What we're building** — Today the MLB dashboard "bets" tab only surfaces a bet
if our *model* thinks it's +EV (model EV ≥ 2%). We're adding a second, independent
edge signal: if a sportsbook's price is **out of line with the market consensus**
(a book's offered price beats the devigged consensus fair by ≥ 2%), that bet gets
flagged too — even when the model is neutral. A bet lights up if **model says edge
OR market says edge**. Each card shows a MODEL / MARKET / BOTH badge so you can see
which signal (or both) fired.

**Key decisions**

1. **Consensus median as the yardstick, not Pinnacle.** We judge each book's price
   against the *median devigged fair across all books quoting that line*, not against
   a single sharp book. Rejected "Pinnacle only" because coverage is too low (many
   lines lack Pinnacle on both sides); the user's intent was "flag it if a book is
   off [the market]," which a consensus captures directly.

2. **Expand the candidate universe (Approach B), not just badge existing cards
   (Approach A).** A market-only edge — model neutral, but a book badly mispriced —
   never reaches the tab today because the 2% model filter discards it before book
   prices are attached. We add a parallel path that finds these. Rejected badge-only
   because it finds zero *new* edges; it's cosmetic.

3. **Isolated module, model path untouched.** The market path is a new
   `find_market_edges()` function reading the canonical per-book odds that `MLB.R`
   already assembles. We merge its output with the model bets by `bet_row_id`. The
   model and its 2% filter are not modified. Rejected reworking the model's internal
   `prediction_set` filter because it's riskier and entangles two concerns.

4. **Leave-one-out consensus; only 1 other book required.** Each book's price is judged
   against the median devigged fair of the **other** books at that line — never against
   a yardstick that includes itself. This is essential because Wagerzon is the softest
   book: if Wagerzon were included in its own consensus, its good price would drag the
   median toward itself and hide its own edge (worst exactly when few books quote the
   line). With leave-one-out we only need the edge book **+ 1 other** (total 2). Rejected
   the original "≥ 3 books, plain median" rule — it suppressed the soft-book edges this
   feature exists to find.

5. **Same 2% threshold for both signals.** Keeps the bar consistent and the mental
   model simple; trivially raised later if market flags prove noisy.

**Risks / push back here**

- **Soft-book consensus can be systematically biased.** Our tracked books skew
  recreational (Wagerzon, Hoop88, BFA, Bet105). If they're collectively shaded the
  same way (e.g. toward popular teams), the consensus median inherits that bias and
  the "edge" is illusory. The model is our defense — when MODEL and MARKET disagree,
  trust the model. The BOTH badge marks the cases where they agree (highest conviction).
- **Volume.** Expanding the universe could noticeably lengthen the bets list. The 2%
  threshold should contain it, but you may want to raise the market threshold to 3%
  after seeing real volume.
- **Single-other-book false positives.** With the floor at 1 other book, a Wagerzon
  price can flag against just one comparison book. If *that* book is also mispriced,
  the "edge" is illusory. Accepted for v1 because Wagerzon edges are the point and
  thin-line coverage matters more than the occasional false flag; the model EV shown
  alongside is the sanity check. Raise the floor to 2 others if false flags become noise.

**Worth understanding** (opt-in)

- **Devigging needs a matched pair.** To strip the vig from a book you need *both*
  sides of the *same* market at the *same* line from that *same* book. Mixing lines or
  books gives a meaningless number. In R terms: think of it like you can't `merge()`
  two data frames on a key that doesn't actually align — the devig "key" is
  (book, market, line), and both sides must share it.
- **`bet_row_id` is a canonical, book-agnostic join key.** It's an md5 hash of the
  *bet identity* — `(id, base_market_type, period, normalized_line, bet_on)` — where
  `base_market_type` strips the `alternate_` prefix (so `"totals"` and `"alternate_totals"`
  collapse) and `period` is F5/FG/etc. It deliberately excludes the book *and* the raw
  market string, because the model uses verbose strings (`"totals_1st_5_innings"`) while
  the canonical book odds only carry the stripped form (`"totals"` + `period="F5"`) — the
  model→canonical mapping is lossy, so only the canonical identity is reproducible by
  *both* paths. Hashing the identity (not the raw string) is what lets us full-join the
  model and market frames and lets the book-prices expansion attach every book's pills to
  one card. (Like a primary key in R you join two tables on.)

---

## Spec Body

### Problem

The bets tab (`mlb_dashboard.R::create_bets_table()`) renders only rows from
`mlb_bets_combined`, which `MLB.R` filters to `ev >= EV_THRESHOLD` (2%), where `ev`
is **model probability vs the book's offered price**. The model is the sole source of
truth. The dashboard already devigs each book's two-sided price for a display-only
FAIR toggle, but those devigged numbers never drive flagging.

Consequence: a bet where the model is neutral but a sportsbook is clearly mispriced
relative to the rest of the market — a genuine soft-line / CLV edge, the core of the
project mission — is invisible. It's filtered out at the 2% model gate before book
prices are ever attached.

### Goal

Flag a bet on the tab if **model EV ≥ 2% OR market EV ≥ 2%**, where market EV measures
how far the best available book price beats the devigged consensus fair. Show on each
card which signal fired (MODEL / MARKET / BOTH) and the underlying numbers.

### Non-goals

- Changing the model, its probabilities, or its 2% filter.
- Surfacing market edges on lines no book quotes on both sides (can't devig).
- Any change to `expand_bets_to_book_prices()` internals.

### Architecture

```
                       (unchanged model path)
  prediction_set ──► format_bets_table (ev>=2%) ──► all_bets ──┐
                                                               │
  book_odds_by_book (canonical per-book odds, already built) ──┤
        │                                                      │
        └──► find_market_edges()  ── market_bets ──────────────┤
             (NEW module: market_edge.R)                       │
                                                               ▼
                                            full_join on bet_row_id
                                            → edge_source, model_ev, market_ev
                                                               │
                                                               ▼
                                            mlb_bets_combined (+ new columns)
                                                               │
                                            expand_bets_to_book_prices() (unchanged)
                                                               │
                                                               ▼
                                            mlb_bets_book_prices  +  dashboard badge
```

#### Component 1 — `market_edge.R` (new file, `Answer Keys/MLB Answer Key/`)

**Signature**

```r
find_market_edges(book_odds_by_book,
                  threshold   = 0.02,
                  min_others  = 1,          # leave-one-out: ≥1 OTHER book to compare against
                  devig_fn    = devig_american)
```

- **Input** `book_odds_by_book`: the named list of canonical per-book frames already
  built at `MLB.R:923`, each with columns `game_id, market, period, side, line,
  american_odds, fetch_time` (and game start time). Reused as-is — no new data source.
- **Output**: a long data frame, one row per flagged market edge, in a schema that
  maps cleanly onto `mlb_bets_combined`:

  | column | meaning |
  |---|---|
  | `bet_row_id` | `compute_bet_row_id(id, market_type, period, line, bet_on)` — the shared helper, identical to the model path |
  | `id` / `game_id` | game id |
  | `market_type` (`"totals"`/`"spreads"`/`"h2h"`), `period`, `line`, `bet_on` | canonical bet identity |
  | `market` | set to `market_type` so downstream `expand_bets_to_book_prices()` (which calls `.derive_market_type(market)`) resolves correctly; `period` is carried explicitly |
  | `home_team`, `away_team`, `pt_start_time` | joined from `.game_id_lookup` / `mlb_odds` so the card header + future-game filter work |
  | `bookmaker_key` | the book offering the edge price (best LOO-EV book) |
  | `odds` | that book's American price |
  | `prob` | the edge book's leave-one-out consensus fair (treated as "true" prob) |
  | `ev` | market EV = `compute_ev(loo_consensus_fair, edge_book_price)` |
  | `n_books` | how many books quote both sides at this line (yardstick uses `n_books − 1`) |
  | `edge_source` | `"market"` |

  After flagging, market rows are collapsed to the **single best-EV line** per
  `(id, market_type, period, side)` (mirrors the model's one-rec-per-side dedup, keeps
  volume sane).

**Algorithm**

1. Stack all books' canonical rows into one long frame, tagged with `bookmaker_key`.
2. Drop rows older than the staleness cutoff (reuse `BOOK_STALENESS_CUTOFF_MIN`).
3. Group by `(game_id, market, period, line, side)` and its opposite side. For each
   **book** that quotes **both** sides at that exact line, devig the pair with
   `devig_fn` → that book's fair prob for `side`. Let `n_books` = how many books produced
   a fair here.
4. **Leave-one-out, per book:** for each book `b`, build its yardstick from the *other*
   books only — `loo_fair_b = median(fair of all books except b)`. Require at least
   `min_others` other books (i.e. `n_books − 1 >= min_others`), else `b` can't be judged.
5. For each book `b`, `ev_b = compute_ev(loo_fair_b, price_b)`. Take the best (max-EV)
   book. If `best_ev >= threshold`, emit a row with that book as `bookmaker_key`, its
   price as `odds`, and `prob = its loo_fair`.
6. Compute `bet_row_id` with the identical hash recipe used in `MLB.R` so ids align.

**Like-for-like is the integrity rule of the whole signal.** The consensus median may
only pool fair probabilities from books that match on **every** identity field:
`game_id`, `market` (totals vs spreads vs moneyline), `period` (F5 vs Full Game vs
F3/F7), `line`, and `side`. A book at a *different* line (Over 8.5 vs Over 9.5), a
*different* period (F5 vs FG), or only quoting *one* side never enters that median.
Consequences baked in:

- The check runs **independently at every distinct line** the books offer — so alt
  lines are covered — but lines are **never blended** (8.5 is never used as evidence
  about 9.5).
- Lines are matched on a **normalized** value (`8.5 == 8.50`) so float/format
  differences don't split a true match.
- A book quoting only one side contributes to *neither* the devig nor `n_books`.
- The edge book is **never** in its own yardstick (leave-one-out), so a soft book can't
  hide its own edge by dragging the consensus toward itself.
- This is stricter than the display grid's ±3-unit "closest line" tolerance, which is
  display-only and never feeds the edge math.

**Why isolated:** the function has one job (find out-of-line books), one input
(canonical odds), one output (a long frame). It can be unit-tested with a hand-built
`book_odds_by_book` fixture and known devig outputs, with zero dependence on the model
or DuckDB.

#### Component 2 — merge in `MLB.R`

**Ordering change.** `book_odds_by_book` is currently built at ~line 923, *after* the
correlation-Kelly call (line 881). To let market bets get the same correlation
adjustment and to merge before that call, the `.game_id_lookup` / `prefetched_long` /
`book_odds_by_book` construction block (currently lines ~905–936) **moves up** to
immediately after the model `all_bets_combined` is finalized (after line 860) and the
old block is deleted from its current spot. The merge then sits between that block and
the existing `adjust_kelly_for_correlation()` call.

Sequence after the move:

1. Model `all_bets_combined` finalized (lines 831–860), now tagged
   `edge_source = "model"`, `model_ev = ev`, `market_ev = NA_real_`. Its `bet_row_id`
   comes from the shared `compute_bet_row_id()` (canonical recipe — see Component 1).
2. Build `book_odds_by_book` (moved block).
3. `market_bets <- find_market_edges(book_odds_by_book)`, tagged
   `edge_source = "market"`, `market_ev = ev`, `model_ev = NA_real_`, with `bet_size`
   computed via `kelly_stake()` (see below). `bet_row_id` from the same shared helper.
4. **`full_join(model, market, by = "bet_row_id")`** — now valid because both sides use
   the canonical recipe. Resolve each row:
   - **model only** → keep as is (`edge_source = "model"`).
   - **market only** → keep the market row (`edge_source = "market"`).
   - **both** → `edge_source = "both"`; keep the **model's** book, `odds`, `prob`,
     `bet_size` (model is the primary pricer); attach `market_ev`. Set
     `ev = pmax(model_ev, market_ev, na.rm = TRUE)` for ranking only.
5. Drop any exact duplicate `bet_row_id` (defensive), `arrange(desc(ev))`. **Do not**
   re-run the line-collapsing model dedup — it already ran on the model set, and
   re-running it would drop legitimately distinct market-only lines.
6. `adjust_kelly_for_correlation()` (existing call) now runs on the merged frame.
7. `expand_bets_to_book_prices()` runs unchanged — it keys on `bet_row_id`, shared by
   both paths, so every card (model, market, or both) gets its full per-book pill row.

**Kelly sizing for market bets:** `bet_size = kelly_stake(ev = market_ev,
book_prob = implied_prob(edge_book_price), bankroll, kelly_mult)` — the existing helper
(`kelly_stake(ev, book_prob, bankroll, kelly_mult)`), treating the LOO consensus as the
true probability via `market_ev`. `bankroll` and `kelly_mult` are already in scope at
the merge point. No new sizing code.

#### Component 3 — dashboard badge (`mlb_dashboard.R`)

- The loader (`SELECT * FROM mlb_bets_combined`) now also carries `edge_source`,
  `model_ev`, `market_ev`. No schema migration needed beyond these columns existing
  on the table that `MLB.R` writes.
- In `create_bets_table()` / the hero-strip renderer: add a small badge driven by
  `edge_source` (MODEL / MARKET / BOTH). When both EVs are present, render both
  (e.g. `Model +3.1% · Market +2.4%`); otherwise render the single driving EV.
- Badge styling lives with the existing card CSS; no new layout system.

### Data flow summary

1. `MLB.R` builds model bets (unchanged) and `book_odds_by_book` (unchanged).
2. `find_market_edges()` produces `market_bets` from `book_odds_by_book`.
3. Full-join → `mlb_bets_combined` gains `edge_source`, `model_ev`, `market_ev`.
4. `expand_bets_to_book_prices()` (unchanged) writes `mlb_bets_book_prices`.
5. Dashboard reads both tables, renders cards with the new badge.

### Error handling & edge cases

- **Only the edge book quotes the line / one-sided quotes:** can't build a leave-one-out
  yardstick (no other book), so skip silently — no market flag.
- **Stale quotes:** dropped by the existing 30-min cutoff before devig.
- **Devig failure** (`NA`/zero odds): `devig_american()` already returns `NA`; those
  books are excluded from the consensus, count toward neither `n_books` nor the median.
- **Off-season / empty odds:** `find_market_edges()` returns a 0-row frame with the
  correct columns; the full-join degenerates to the model set; tab behaves as today.
- **Hash drift:** the `bet_row_id` recipe is duplicated in two places — extract it to a
  shared helper (`compute_bet_row_id()`) so the model and market paths can never diverge.
- **Line representation:** group on a normalized/rounded line so `2.5` and `2.50` match;
  devig only within an exact (book, line) pair.

### Testing

- **Unit (`market_edge` test):** hand-built `book_odds_by_book` fixture with known
  American odds → assert the leave-one-out fair, `n_books`, the flagged book, and EV.
  Include: a clear edge, a sub-threshold near-miss, a **2-book line (edge book + 1 other)
  that MUST flag** (validates the floor and that leave-one-out excludes the edge book),
  a **1-book line that must be skipped** (no other book to compare against), a one-sided
  book (excluded from consensus), and a stale row (dropped).
- **Integration:** run the `MLB.R` pipeline on copied live DBs; confirm
  `mlb_bets_combined` gains the new columns, market-only bets appear, BOTH bets are
  labeled, and `mlb_bets_book_prices` still populates for the new rows.
- **Render check:** render the dashboard against copied live DBs and visually confirm
  the MODEL / MARKET / BOTH badge and dual-EV display (per the project's
  verify-by-rendering rule — static diff review is not sufficient).
- **Regression:** the model-only set is a strict subset of the new set; assert no
  previously-flagged model bet disappears.

### Version control

- **Branch / worktree:** `worktree-mlb-market-edge-flagging` (already created at
  `.claude/worktrees/mlb-market-edge-flagging`). All work happens here; merge to `main`
  only after the render check + explicit approval; remove the worktree + branch after merge.
- **DuckDB caution:** do **not** symlink any `.duckdb` into the worktree (WAL-loss
  hazard per CLAUDE.md). Test by *copying* live DBs into the worktree, or test from
  `main` after merge.
- **Files created:** `Answer Keys/MLB Answer Key/market_edge.R`; `Answer Keys/tests/test_market_edge.R`.
- **Files modified:**
  - `Answer Keys/MLB Answer Key/odds_screen.R` — add the shared `compute_bet_row_id()`
    helper (sits beside `.derive_market_type` / `.derive_period`, which it uses).
  - `Answer Keys/MLB Answer Key/MLB.R` — replace the inline `bet_row_id` hash with the
    shared helper; move the `book_odds_by_book` block up; insert the market merge.
  - `Answer Keys/MLB Dashboard/mlb_dashboard.R` — badge + dual-EV + `.edge-badge` CSS.
- **Commits:** (1) `compute_bet_row_id()` helper + test; (2) `market_edge.R` module +
  unit tests; (3) MLB.R merge (move block + full_join); (4) dashboard badge + CSS;
  (5) docs.

### Documentation

- Update `Answer Keys/CLAUDE.md` ("Odds screen" / bets-tab section) to document the
  dual-source flagging, the `edge_source`/`model_ev`/`market_ev` columns, and the
  `find_market_edges()` module.
- Note the new columns on `mlb_bets_combined` and the consensus/guard parameters.
- Docs land in the **same merge** to `main` as the code.
