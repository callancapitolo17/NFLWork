# Wagerzon MLB F5 3-Way Market — Design Spec

**Author:** Claude (with Callan)
**Date:** 2026-05-03
**Branch:** `feature/wagerzon-mlb-3way`
**Status:** Draft — pending review

---

## Goal

Match Wagerzon's **F5 3-way moneyline** market (`MLB - 1ST 5 INN WINNER (3-WAY)`) end-to-end: scrape the feed, devig the three-way prices, compare against historical sample margin distribution, and surface +EV bets on the MLB dashboard's straight-bets tab. Coexist with the existing 2-way F5 ML matching path (don't replace it).

## Why this exists

Wagerzon offers two distinct F5 ML markets:

1. **2-way + push refund** (existing): home/away with tie refund. Lives in `idgmtyp=15` children of `idgmtyp=10` FG parents. Currently scraped → `market = "spreads_h1"` (carries spread + total + 2-way ML on one row) → matched by `compare_moneylines_to_wagerzon` with `non_push <- margins[margins != 0]` exclusion semantics.
2. **3-way (home/away/draw, no refund)** (new): explicit draw price, treated as a real outcome. Lives in `lg=1280` as `idgmtyp=29` parent games. Currently **scraped page is fetched but the data is silently discarded** because the scraper has no branch for `idgmtyp=29`.

The 3-way market typically carries higher vig (~9.7% in the recon snapshot vs ~4-5% on the 2-way) but also exposes the draw as a separately bettable outcome. F5 ties are common enough (~16-18% historical) that the draw side can be +EV on a slate-by-slate basis. Today we book exactly zero of these bets.

## Recon evidence

Live capture from `https://backend.wagerzon.com/wager/NewScheduleHelper.aspx?WT=0&lg=1280` on 2026-05-03 (auth: primary Wagerzon account):

```json
{
  "result": {
    "listLeagues": [[{
      "Description": "MLB - 1ST 5 INN WINNER (3-WAY)",
      "Games": [
        {
          "vtm": "1H CLE GUARDIANS 3WAY",
          "htm": "1H ATHLETICS 3WAY",
          "gpd": "GAME",
          "idgmtyp": 29,
          "idlg": 1280,
          "idspt": "MLB",
          "GameLines": [{
            "voddst": "105",      // away ML
            "hoddst": "130",      // home ML
            "vspoddst": "475",    // ← DRAW PRICE (3rd outcome)
            "voddsh": "+105",     // away ML (display variant)
            "hoddsh": "+130",     // home ML (display variant)
            "vsph":   "+475"      // draw (display variant)
            // ...standard empty spread/total fields
          }]
        },
        // ...5 more games + occasional EmptyGame placeholders
      ]
    }]]
  }
}
```

**Six games today, all with the same shape.** Sanity check on the first game's vig:
- p_away_implied = 100 / (100 + 105) = 0.488
- p_home_implied = 100 / (100 + 130) = 0.435
- p_draw_implied = 100 / (100 + 475) = 0.174
- Sum = 1.097 → 9.7% vig
- Devigged: p_away = 0.444, p_home = 0.397, **p_draw = 0.158**

That ~16% draw probability matches historical F5 tie frequency in MLB.

## Architecture

```
wagerzon_odds/scraper_v2.py
  └── parse_odds()
        ├── (existing) iterate listLeagues
        ├── (existing) parent loop: idgmtyp 10/15/19/25/30/31/35/44/47/66
        └── NEW: standalone parent loop branch for idgmtyp == 29
              └── parse_3way_line()  (NEW helper)
                    Emits row with market="h2h_3way_1st_5_innings", period="f5",
                    away_ml/home_ml/draw_ml populated, all spread/total fields NULL

wagerzon_odds/wagerzon.duckdb
  └── mlb_odds (and cbb_odds, nfl_odds, college_baseball_odds for symmetry)
        New column: draw_ml INTEGER (NULL for non-3-way rows)
        Schema upgrade via init_database()'s idempotent
        ALTER TABLE ... ADD COLUMN IF NOT EXISTS

Answer Keys/Tools.R::get_wagerzon_odds()
  └── (existing) for-row dispatcher emits per-market_type records
        NEW branch: if row$market starts with "h2h_3way_", emit one record
        with market_type="h2h_3way", odds_away/odds_home/odds_draw populated

Answer Keys/Tools.R::compare_alts_to_samples()
  └── (existing) filter regex + suffix dispatcher + per-row branch chain
        EXTEND filter regex: include grepl("^h2h_3way_", market)
        EXTEND suffix dispatcher: add 1st_5_innings$ → "f5"
        NEW per-row branch: ^h2h_3way_ →
          - p_home = sum(margin > 0) / N
          - p_away = sum(margin < 0) / N
          - p_draw = sum(margin == 0) / N    (no exclusion)
          - american_prob_3way(odds_away, odds_home, odds_draw) → devigged probs
          - Emit up to 3 bets (Home/Away/Tie) with EV gate

Answer Keys/Tools.R::american_prob_3way()  (NEW helper)
  └── Sum-and-normalize devig: each implied prob / sum of all three

Answer Keys/MLB Answer Key/MLB.R
  └── (existing) compare_alts_to_samples called per book — no change
        The new 3-way market flows through the existing call site
        because it's just another offshore-odds market type.

mlb_bets_combined (downstream output)
  └── New rows: market="h2h_3way_1st_5_innings", bookmaker_key="wagerzon",
        bet_on ∈ {"<home_team>", "<away_team>", "Tie"}
        Dashboard market-filter built dynamically — picks up automatically.
```

## Component responsibilities

### 1. Scraper: `wagerzon_odds/scraper_v2.py`

**Responsibility:** Pull the lg=1280 league out of the existing combined response and emit one row per game in the standard `mlb_odds` shape, with the new `draw_ml` field populated.

**New code (~50 lines):**
- `init_database()`: add `draw_ml INTEGER` to the CREATE TABLE statement; add `ALTER TABLE {table_name} ADD COLUMN IF NOT EXISTS draw_ml INTEGER` immediately after for each affected sport table (mlb_odds, cbb_odds, nfl_odds, college_baseball_odds).
- New parser helper `parse_3way_line(line, game_id, period, market, base) → dict | None` mirroring `parse_moneyline_only` but extracting `voddst`/`hoddst`/`vspoddst` and writing them into `away_ml`/`home_ml`/`draw_ml`.
- New parent-loop branch in `parse_odds()` for `idgmtyp == 29`. Lives alongside the existing parent-game loop (NOT inside the `GameChilds` loop — `idgmtyp=29` is a parent, not a child). Strip the trailing " 3WAY" from `vtm`/`htm` before team-name resolution. Skip rows where `voddst` is empty (Wagerzon posts placeholder games with no prices).
- The `mlb_odds` INSERT uses named columns (`scraper_v2.py:642-643`), so the new `draw_ml` column gets included automatically once added to the `columns` list.

**Why a separate `parse_3way_line` (vs reusing `parse_game_line`):** semantic clarity. `parse_game_line` returns spread + total + 2-way ML; `parse_3way_line` returns 3-way ML only. Mixing them would force `parse_game_line` to grow a `vspoddst` field that's only ever populated in 3-way contexts. Cleaner to keep them separate.

### 2. Reader: `Answer Keys/Tools.R::get_wagerzon_odds()`

**Responsibility:** When iterating raw `mlb_odds` rows, detect 3-way markets and emit a single record per game (not the spread+total+ML triplet that 2-way markets produce).

**New code (~15 lines):** Add a branch in the per-row dispatcher (around `Tools.R:3290`) that runs **before** the existing `if (!is.na(row$away_spread))` spread branch:

```r
# 3-way moneyline (Wagerzon F5 3-way market): emit one record with all 3 prices
if (!is.na(row$draw_ml) && grepl("^h2h_3way_", row$market)) {
  three_way_rec <- c(base, list(
    market = row$market,
    market_type = "h2h_3way",
    line = NA_real_,
    odds_away = row$away_ml,
    odds_home = row$home_ml,
    odds_draw = row$draw_ml,
    odds_over = NA_integer_,
    odds_under = NA_integer_
  ))
  result_list[[length(result_list) + 1]] <- three_way_rec
  next   # don't fall through to spread/total/h2h record emission
}
```

**Edge case:** If a future sport's scraper writes `draw_ml` on a non-3-way row (defensive), the `grepl("^h2h_3way_", row$market)` check prevents misclassification.

### 3. Devig helper: `Answer Keys/Tools.R::american_prob_3way()`

**Responsibility:** Convert three American odds → three devigged probabilities summing to 1.0.

**New code (~15 lines, near existing `american_prob`):**

```r
american_prob_3way <- function(odds_away, odds_home, odds_draw) {
  # Convert each leg to implied prob, then normalize so the three sum to 1.
  implied <- function(odds) {
    if (is.na(odds) || odds == 0) return(NA_real_)
    if (odds > 0) 100 / (odds + 100)
    else (-odds) / (-odds + 100)
  }
  p_a <- implied(odds_away)
  p_h <- implied(odds_home)
  p_d <- implied(odds_draw)
  if (any(is.na(c(p_a, p_h, p_d)))) {
    return(list(p_away = NA_real_, p_home = NA_real_, p_draw = NA_real_))
  }
  total <- p_a + p_h + p_d
  list(p_away = p_a / total, p_home = p_h / total, p_draw = p_d / total)
}
```

**Why a new helper instead of extending `american_prob`:** keeps the existing 2-way helper's signature stable. `american_prob` is used in dozens of call sites; changing its signature or adding optional args risks subtle bugs.

### 4. Matcher: `Answer Keys/Tools.R::compare_alts_to_samples()`

**Responsibility:** Add a new `else if` branch to the per-row chain that prices 3-way markets against the sample margin distribution.

**Three small edits:**

(a) Filter regex extension (around `Tools.R:4525`):
```r
alt_odds <- offshore_odds %>%
  filter(
    grepl("^alternate_", market) |
    grepl("^team_totals_", market) |
    grepl("1st_3_innings", market) |
    grepl("1st_7_innings", market) |
    market == "odd_even_runs" |
    grepl("^h2h_3way_", market)            # NEW
  )
```

(b) Suffix dispatcher entry (around `Tools.R:4555-4566`):
```r
suffix <- if (grepl("1st_3_innings$", row$market)) {
  "f3"
} else if (grepl("1st_5_innings$", row$market)) {   # NEW (only for h2h_3way today; F5 mains go via Odds API)
  "f5"
} else if (grepl("1st_7_innings$", row$market)) {
  "f7"
} else if (row$market == "odd_even_runs") {
  "fg"
} else {
  sub(".*_", "", row$market)
}
```

(c) New `else if` branch in the per-row chain (immediately after the `odd_even_runs` branch, before the for-loop's closing `}`):
```r
} else if (grepl("^h2h_3way_", row$market) &&
           !is.na(row$odds_home) && !is.na(row$odds_away) && !is.na(row$odds_draw)) {
  # 3-way moneyline (Wagerzon F5 3-way). Tie is its OWN outcome, not a refund.
  col_name <- paste0(margin_col, "_", period)   # game_home_margin_period_F5
  if (!col_name %in% names(sample_df)) next
  margins <- sample_df[[col_name]]
  margins <- margins[!is.na(margins)]
  if (length(margins) == 0) next

  p_home <- sum(margins > 0) / length(margins)
  p_away <- sum(margins < 0) / length(margins)
  p_draw <- sum(margins == 0) / length(margins)

  probs <- american_prob_3way(row$odds_away, row$odds_home, row$odds_draw)
  if (any(is.na(unlist(probs)))) next

  home_ev <- compute_ev(p_home, probs$p_home)
  away_ev <- compute_ev(p_away, probs$p_away)
  draw_ev <- compute_ev(p_draw, probs$p_draw)
  home_size <- kelly_stake(home_ev, probs$p_home, bankroll, kelly_mult)
  away_size <- kelly_stake(away_ev, probs$p_away, bankroll, kelly_mult)
  draw_size <- kelly_stake(draw_ev, probs$p_draw, bankroll, kelly_mult)

  if (home_ev >= ev_threshold) {
    all_bets[[length(all_bets) + 1]] <- tibble(
      id = game_id, home_team = row$home_team, away_team = row$away_team,
      pt_start_time = pt_start_time, bookmaker_key = book_key,
      market = row$market, bet_on = row$home_team,
      line = NA_real_, bet_size = home_size, ev = home_ev,
      odds = row$odds_home, prob = p_home
    )
  }
  if (away_ev >= ev_threshold) {
    all_bets[[length(all_bets) + 1]] <- tibble(
      id = game_id, home_team = row$home_team, away_team = row$away_team,
      pt_start_time = pt_start_time, bookmaker_key = book_key,
      market = row$market, bet_on = row$away_team,
      line = NA_real_, bet_size = away_size, ev = away_ev,
      odds = row$odds_away, prob = p_away
    )
  }
  if (draw_ev >= ev_threshold) {
    all_bets[[length(all_bets) + 1]] <- tibble(
      id = game_id, home_team = row$home_team, away_team = row$away_team,
      pt_start_time = pt_start_time, bookmaker_key = book_key,
      market = row$market, bet_on = "Tie",
      line = NA_real_, bet_size = draw_size, ev = draw_ev,
      odds = row$odds_draw, prob = p_draw
    )
  }
```

### 5. Schema additions

**`wagerzon.duckdb`:** `mlb_odds`, `cbb_odds`, `nfl_odds`, `college_baseball_odds` each get `draw_ml INTEGER` via idempotent ALTER. NULL for all existing rows + every non-3-way row going forward.

**`mlb_bets_combined`:** Tibble shape unchanged. New `market` value `h2h_3way_1st_5_innings`. New possible `bet_on` value `"Tie"`. Dashboard market filter is dynamic — no JS change needed.

## Data flow (end-to-end)

```
1. Wagerzon scraper run
   └── scraper_v2.py fetches NewScheduleHelper.aspx?lg=…,1280,… (already in URL)
   └── parse_odds() detects idgmtyp=29 leagues
   └── For each non-empty game in that league:
         strip " 3WAY" suffix → resolve teams
         emit row: market=h2h_3way_1st_5_innings, draw_ml=vspoddst
   └── INSERT into mlb_odds (named columns; new draw_ml populated only on 3-way rows)

2. MLB.R pipeline run
   └── get_wagerzon_odds("mlb") loads raw mlb_odds rows
   └── for each row, dispatcher detects market=h2h_3way_… → emits one h2h_3way record
   └── compare_alts_to_samples() processes the row:
         filter regex catches it
         suffix dispatcher resolves "1st_5_innings" → "f5" → "F5"
         per-row branch: read game_home_margin_period_F5, devig, EV-gate, emit bets
   └── Up to 3 bets per game (home/away/tie) added to mlb_bets_combined

3. Dashboard
   └── Reads mlb_bets_combined, renders Bets tab
   └── Market filter dropdown picks up "h2h_3way_1st_5_innings" automatically
   └── User sees the new market alongside existing h2h_1st_5_innings (2-way), can place via existing wagerzon auto-place flow
```

## Testing strategy

### Unit tests (`Answer Keys/tests/test_compare_alts_to_samples.R` — extend existing file)

- **Test 1 (devig):** `american_prob_3way(+105, +130, +475)` returns `(0.444, 0.397, 0.158)` ± tolerance, sums to 1.0
- **Test 2 (devig NA):** any NA odds → all three NA
- **Test 3 (3-way fair-prob):** synthetic samples with 700 home wins / 200 away wins / 100 ties → assert `compare_alts_to_samples` emits home + tie bets at correct probs (away is -EV at any reasonable price)
- **Test 4 (3-way NA odds):** missing draw price → no bets emitted
- **Test 5 (3-way real-shape fixture):** emulate the recon JSON's exact field shape → end-to-end through `get_wagerzon_odds` → `compare_alts_to_samples` → assert correct row count and probs

### Scraper tests (`wagerzon_odds/tests/` — new file)

- **Test (parse):** feed a fixture mirroring the recon JSON → assert `parse_3way_line` returns expected row shape with `draw_ml` populated, `away_spread`/`total`/etc. NULL
- **Test (skip empty):** feed a JSON with `EmptyGame: true` → assert no row emitted

### Schema migration test

- **Test (migration idempotent):** create a `wagerzon.duckdb` with the OLD 18-column schema, run `init_database()` once → confirm `draw_ml` column added; run again → confirm no error and column count unchanged

### End-to-end smoke (manual, not committed)

- Copy live DBs into worktree (NEVER symlink)
- Run `python run.py mlb`
- Confirm `mlb_bets_combined` contains rows with `market = h2h_3way_1st_5_innings` from `wagerzon` (when edges exist)
- Hand-verify one row: prob × profit calc matches `compute_ev` output

## Risks & mitigations

| Risk | Likelihood | Mitigation |
|---|---|---|
| Wagerzon changes the field name (`vspoddst` → something else) | Low | Spec includes recon evidence; if it changes, scraper fails loudly (NULL `draw_ml`) and the matcher's `is.na(odds_draw)` guard prevents bad bets |
| F5 ties are rare in our sample, so `p_draw` is poorly estimated | Real but acceptable | Sample size is 12,719 historical games; ~16% are F5 ties, so ~2K observations — not noise-free but reasonable. Larger sample mitigates |
| Vig is ~10% on 3-way (vs ~4-5% on 2-way) so edges are rarer | Real but acceptable | Same EV threshold (0.02) as 2-way; sparse + accurate is better than dense + noisy |
| Wagerzon settles ties as "no action" instead of paying the draw price (bug in our model) | Verify before placing | Test by placing a small live bet and observing settlement; do NOT scale until confirmed |
| The new `draw_ml` column breaks a downstream reader I missed | Low | Verified all readers use named-column access (`SELECT *` then `row$col`); INSERT uses named columns |
| Idempotent ALTER fires on every scraper start (microseconds) | Negligible | Acceptable cost for zero-ops upgrade |

## Out of scope

- F3 / F7 3-way (recon only confirmed F5 — would need fresh recon to verify those leagues exist)
- Replacing the existing 2-way matching path
- A new dashboard tab — the existing Bets tab handles new market types via dynamic filter
- Schema migration tooling beyond the inline ALTER (no external migration script)
- 3-way support in `compare_moneylines_to_wagerzon` (the F5 mains pricer goes through Odds API, not offshore — and the Odds API doesn't expose Wagerzon's 3-way market). All 3-way pricing goes through `compare_alts_to_samples`.

## Version control plan

- **Branch:** `feature/wagerzon-mlb-3way` (already created via worktree)
- **Worktree:** `/Users/callancapitolo/NFLWork/.worktrees/wagerzon-mlb-3way`
- **Commits:** ~7 small commits along the TDD sequence (one per logical task: schema, scraper parser, reader branch, devig helper, matcher branch, end-to-end test, docs)
- **Merge:** `--no-ff` to main only after explicit user approval
- **Cleanup:** `git worktree remove` + `git branch -d` after merge

## Documentation plan

- `Answer Keys/CLAUDE.md` — add 3-way market to MLB pipeline bullet
- `Answer Keys/MLB Dashboard/README.md` — add `h2h_3way_1st_5_innings` to Markets section
- `wagerzon_odds/CLAUDE.md` — note `draw_ml` column + `idgmtyp=29` handling pattern

## Open questions

None — all four open questions resolved via prior conversation:
1. ✅ Market name: `h2h_3way_1st_5_innings`
2. ✅ Draw label: "Tie"
3. ✅ Don't gate to Wagerzon (leave matcher generic — book prefix not enforced)
4. ✅ Schema migration: inline idempotent ALTER (Option B from the brainstorm)
