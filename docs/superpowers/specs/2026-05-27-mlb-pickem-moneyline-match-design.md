# MLB Dashboard — Pick'em (0-line spread) ↔ Moneyline Matching

## Review Pack

**What we're building**
On the MLB bets-tab odds screen, a bet that is a spread at **line 0** (a
pick'em — e.g. Wagerzon "PHI 0" in the first 3 innings) is economically a
draw-no-bet moneyline: you back a team to lead, a tie refunds. Today the
per-book comparison cells match it to each book's *run line* (±0.5), which is
a different bet (no push), so the comparison is apples-to-oranges and books
like DraftKings show "—". This change makes each book's cell show that book's
true pick'em price, derived from its 2-way or 3-way period **winner** market.
Display-only — the model's bet generation is untouched.

**Key decisions**
1. **Trigger = spread/alt-spread with `line == 0` exactly** (any period).
   Rejected `|line| ≤ 0.5` because a ±0.5 run line has no push and is a
   genuinely different bet — only a 0-handicap equals a moneyline.
2. **Derive the price inside the matcher** (Approach A), via one pure helper
   `derive_pickem_american()`, rather than synthesizing fake "spread-0" rows
   upstream (Approach B). Keeps real and derived data separate in the frame
   and the DNB math in one tested place.
3. **Source priority per book: 2-way winner → 3-way winner (devig, drop tie,
   renormalize) → true 0-handicap → "—".** The 3-way path (collapse a
   Home/Tie/Away market to draw-no-bet) is what lets most books produce a
   pick'em number instead of "—".
4. **Display-only.** The model already prices the WZ pick'em correctly as a
   draw-no-bet (verified in `compare_alts_to_samples`); we do not add
   per-book bet generation/sizing. The derived prices populate comparison
   cells only.
5. **v1 coverage: DraftKings + FanDuel.** Both books get period-winner
   capture wired end-to-end so the F3 (and F5/F7) pick'em cards show real
   numbers for the two books most likely to be acted on. The 3-way path
   ships ready and is exercised against Wagerzon's F5 3-way *if* its
   scraper output already carries home/away/tie (to confirm in planning).
   Other books (Pinnacle/Bet105/BKM/Hoop88) show "—" on pick'em cards until
   wired — the matcher/helper are book-agnostic, so each addition is just
   scraper capture.

**Risks / push back here**
- **Coverage is the value lever, and it's incremental.** v1 wires DK + FD;
  the 3-way path ships but is only exercised in v1 if WZ's F5 3-way scraper
  output is already consumable. Pinnacle/Bet105/BKM/Hoop88 cells stay "—" on
  pick'em cards until each book's winner market is wired (fast-follow).
- **FD winner-market availability is a planning verification.** FD currently
  flows a run line to the F3 PK card; whether FD also posts a capturable F3
  2-way winner (vs. only the run line) is unconfirmed. If FD posts no F3/F5/
  F7 winner at all, FD will stay "—" in v1 despite being in scope, and we'd
  decide in planning whether to defer FD or surface their main-game ML
  another way.
- **Schema add.** `mlb_bets_book_prices` gains a `derived_fair_odds` column.
  It's additive and NULL for all non-pick'em rows, but it's still a schema
  change the dashboard loader must tolerate on first deploy.
- **"True 0-handicap" books are rare.** Almost every book expresses the
  pick'em as a moneyline, so the 0-handicap branch may never fire in
  practice. Kept in the priority list for correctness, but it's low-value.

**Worth understanding** (opt-in)
- **Draw-no-bet as conditional probability.** `p_dnb = p_win / (p_win +
  p_loss)` is just renormalizing after dropping the tie — the same move as
  `prop.table()` on a subset in R. A 2-way winner market is *already* in this
  form (no tie outcome), which is why devigging it directly gives the DNB. A
  3-way market carries the tie, so you devig all three then drop it.
- **Why the matcher, not the cell, computes FAIR here.** Normal cells store a
  raw quote and the cell devigs the pick/opposite pair on the fly. A
  3-way-derived price has no single raw pair the cell could devig correctly,
  so the matcher precomputes both the raw-implied DNB and the devigged DNB and
  stores them. This is the kind of "push the computation to where the inputs
  live" boundary decision that keeps the cell renderer dumb and testable.

---

## Problem

The bets-tab card for a first-3-innings "PHI 0" pick'em renders a per-book
comparison grid. The matcher (`odds_screen.R::expand_bets_to_book_prices`)
treats the bet as a spread at line 0 and, via the ±3.0-unit closest-line
picker, matches each book's nearest **run line** (e.g. FD −0.5, DK nothing).

Two problems:
1. A run line at ±0.5 is **not** the same bet as a 0-handicap. `-0.5` makes a
   tie a loss; `+0.5` makes a tie a win; the 0-handicap **pushes** on a tie.
   So the comparison comingles different payout structures.
2. DraftKings posts no F3 run line *and* no F3 0-handicap, so its cell is
   "—" — even though DK posts a 2-way "1st 3 Innings" winner (COL +140 / LA
   −180) that **is** the pick'em equivalent.

The model side is correct: `compare_alts_to_samples` (Tools.R:5082-5094)
prices a 0-line spread by excluding ties (`margins != -home_spread` with
`home_spread == 0`) and computing `P(lead | no tie)` — identical to its 2-way
moneyline branch (Tools.R:5223-5225). So the fix is purely on the
display/matching side.

## Architecture (Approach A)

Three new/changed pieces, all in the dashboard display path:

### 1. `derive_pickem_american()` — pure helper (`odds_screen.R`)

Signature (conceptual):

```
derive_pickem_american(home_raw, away_raw, tie_raw = NA) ->
  list(home_raw_dnb, away_raw_dnb, home_fair_dnb, away_fair_dnb)   # American odds
```

- **2-way source** (`tie_raw` is NA): a 2-way winner already excludes ties.
  - `raw_dnb` = the raw odds unchanged (`home_raw`, `away_raw`).
  - `fair_dnb` = probit 2-way devig (`Tools.R::devig_american`) → probs →
    American.
- **3-way source** (`tie_raw` present):
  - `fair_dnb`: `devig_american_3way(c(away_raw, home_raw, tie_raw))` →
    `p_away, p_home, p_tie`; drop tie → `p_home/(p_home+p_away)`,
    `p_away/(p_home+p_away)` → American.
  - `raw_dnb`: same drop-tie-renormalize on the **raw** implied probs (no
    devig) → American.

No I/O, no DB. Unit-tested against hand-computed numbers, including the
DK ARI/PHI example and a synthetic 3-way.

### 2. Canonical representation of period winners

So the matcher can find a book's winner market, scrapers must surface it in
the canonical book frame (`scraper_to_canonical` → `normalize_book_odds_frame`):

- **2-way winner** → existing `h2h` shape: two rows (home_team, away_team),
  `line = NA`, `market = h2h`, tagged with the period (e.g. market name
  `h2h_1st_3_innings` so `.derive_period` yields `F3`).
- **3-way winner** → new `h2h_3way` shape: three rows (home_team, away_team,
  `"Tie"`), `line = NA`, period-tagged. `scraper_to_canonical` gains a 3-way
  branch that emits the tie row; `normalize_book_odds_frame` carries
  `market = h2h_3way`. (Wagerzon's `h2h_3way_1st_5_innings` is the candidate
  first consumer — confirm during planning that the WZ scraper output carries
  home/away/tie odds, not just the pricing-path representation.)

### 3. Matcher rule (`expand_bets_to_book_prices`)

When the bet is `market_type ∈ {spreads, alternate_spreads}` **and**
`line == 0`:

- The related market set for that bet becomes the period **winner** markets
  (`h2h`, `h2h_3way`) — *not* the run-line union.
- For each book, gather its winner rows for `game_id` + `period`:
  - 2-way: home/away `h2h` rows → `derive_pickem_american(home, away)`.
  - 3-way: home/away/Tie `h2h_3way` rows → `derive_pickem_american(home,
    away, tie)`.
  - true 0-handicap: a `spreads`/`alternate_spreads` row at exactly `line 0`
    (rare) → use its odds directly via the 2-way path.
  - none of the above → emit no row (cell shows "—").
- Emit a `pick` row and an `opposite` row with:
  - `american_odds` = `*_raw_dnb` (RAW display),
  - `derived_fair_odds` = `*_fair_dnb` (FAIR display),
  - `is_exact_line = TRUE`, `line_quoted = 0`.

This branch is gated strictly on `line == 0`, so **every other card — main
spreads, totals, alt-spreads including the DK FG fix — is unaffected.**

### 4. Schema + cell display

- `mlb_bets_book_prices` gains `derived_fair_odds DOUBLE` (NULL for normal
  rows). MLB.R Phase 8 `CREATE TABLE` updated.
- `book_cell.R::render_book_cell`: when `derived_fair_odds` is non-NULL, show
  `american_odds` on the RAW toggle and `derived_fair_odds` on the FAIR
  toggle, skipping the on-the-fly pair devig. Normal rows behave exactly as
  today.
- A derived pick'em cell shows no line tag (it is the exact pick'em).

## Data flow

```
DK scraper: classify "1st 3 Innings" as a 2-way period winner (was dropped)
  → dk.duckdb mlb_odds row with period=F3, moneyline odds
get_dk_odds() / scraper_to_canonical()
  → canonical h2h rows (market=h2h_1st_3_innings → period F3, line NA)
Wagerzon scraper (existing h2h_3way_1st_5_innings)
  → scraper_to_canonical 3-way branch → h2h_3way rows (home/away/Tie)
MLB.R Phase 7b: expand_bets_to_book_prices(all_bets_combined, book_odds_by_book)
  → for a line==0 spread bet, derive_pickem_american per book
  → mlb_bets_book_prices rows with american_odds (raw DNB) + derived_fair_odds
Dashboard: create_bets_table → render_book_cell shows raw/fair DNB per toggle
```

## Edge cases & error handling

- **Book has only run lines** (no winner, no 0-handicap): cell "—". No
  fallback to the run line (the user's explicit call — avoids the
  apples-to-oranges).
- **3-way devig fails / degenerate odds** (e.g. one outcome missing,
  `uniroot` non-convergence): helper returns NA for that book; cell "—".
- **2-way winner with a missing side**: NA → "—" (mirrors the existing
  paired-side guard philosophy).
- **Non-zero spreads**: untouched — they never enter this branch.
- **`derived_fair_odds` NULL on legacy/first-deploy rows**: cell falls back to
  existing devig path; safe.

## Testing

- `derive_pickem_american()` unit tests (R): 2-way (hand-checked devig),
  3-way (hand-checked 3-way devig + drop-tie), NA/degenerate inputs.
- `expand_bets_to_book_prices` tests: a `line==0` spread bet matches a book's
  2-way `h2h` winner and a book's `h2h_3way` winner; a book with only run
  lines yields no row; a non-zero spread bet is unchanged (regression).
- DK scraper parser test: "1st 3 Innings" classifies as a 2-way period
  winner and yields a moneyline row at period F3.
- Manual render check against copied live DBs per the project's
  verify-by-rendering rule before merge.

## Documentation impact

- `Answer Keys/CLAUDE.md` — odds-screen section: document the pick'em→winner
  matching rule and the `derived_fair_odds` column.
- `mlb_sgp/README.md` (or DK scraper docs) — note that DK period winner
  markets ("1st N Innings") are now captured as moneylines.

## Out of scope (v1)

- Per-book bet generation/sizing on derived pick'em prices (display-only).
- Wiring books beyond DraftKings + FanDuel (and opportunistically Wagerzon's
  existing F5 3-way). Pinnacle/Bet105/BKM/Hoop88 are incremental fast-follow.
- The separate DK FG alt-spread `abs(line)` bucket fix (its own branch,
  already implemented and verified).
