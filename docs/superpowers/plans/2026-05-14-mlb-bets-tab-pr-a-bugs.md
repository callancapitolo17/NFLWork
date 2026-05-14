# MLB Bets Tab — PR A (bugs) Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Fix three visible bugs on the MLB Dashboard bets-tab V8 card layout: game times rendering in UTC, Bet105/BFA alt-total quotes missing from per-book grids, and the spread-opposite cell wrongly amber-tagged on every spread bet.

**Architecture:** Two of the three fixes (alt-totals join, opposite-line tagging) are isolated mutations inside `Answer Keys/MLB Answer Key/odds_screen.R` and are TDD-friendly using the existing `testthat` suite at `Answer Keys/tests/test_odds_screen.R`. The third (timezone) is a render-time fix in `Answer Keys/MLB Dashboard/mlb_dashboard.R` that mirrors the legacy parlay-table pattern (`mlb_dashboard.R:515-524`) — emit ISO+Z, let the browser convert with `toLocaleString()`. No new files; no new dependencies.

**Tech Stack:** R 4.x, testthat, dplyr, tibble (existing); plain DOM JavaScript for the timezone fix (existing pattern).

**Worktree:** Before Task 1, create a fresh worktree off `main` named `mlb-bets-tab-pr-a-bugs` (e.g. via `EnterWorktree({name: "mlb-bets-tab-pr-a-bugs"})` if you have the Superpowers worktree skill, or `git worktree add .claude/worktrees/mlb-bets-tab-pr-a-bugs -b worktree-mlb-bets-tab-pr-a-bugs main` otherwise). All paths in this plan assume that worktree as the working directory. The spec at `docs/superpowers/specs/2026-05-14-mlb-bets-tab-followups-design.md` should already be on `main`.

---

## File Structure

```
Answer Keys/
├── MLB Answer Key/
│   └── odds_screen.R                         (modify: 2 fixes — bugs 2 + 3)
├── MLB Dashboard/
│   └── mlb_dashboard.R                       (modify: 1 fix — bug 1)
├── tests/
│   └── test_odds_screen.R                    (modify: add tests for bugs 2 + 3)
└── CLAUDE.md                                 (modify: pitfall note — _fg suffix)
```

Each file's responsibility:

- **`odds_screen.R`** — pure helpers for the bets-tab odds-screen pipeline. We touch `.derive_period()`, `.derive_market_type()`, and `expand_bets_to_book_prices()`. No new functions; no signature changes.
- **`mlb_dashboard.R`** — V8 card layout R/JS rendering. We change one `mutate()` line in `create_bets_table()` and add one `<script>` block to the bets-tab page template.
- **`test_odds_screen.R`** — existing testthat file. We add 3 new `test_that(...)` blocks. Existing tests must continue to pass.
- **`CLAUDE.md`** (Answer Keys) — docs. Add one bullet to the Pitfalls section noting the dual-suffix convention.

---

## Task 1: Bug 3 — spread-opposite cell wrongly tagged as alt-line

**Why first.** Smallest, most contained change. Pure logic fix inside one function with no UI ramifications. Lays the foundation for the testthat workflow we'll reuse in Task 2.

**Files:**
- Modify: `Answer Keys/MLB Answer Key/odds_screen.R:151-249` (`expand_bets_to_book_prices`)
- Test:   `Answer Keys/tests/test_odds_screen.R` (append new `test_that` block)

- [ ] **Step 1: Write the failing test**

Append to `Answer Keys/tests/test_odds_screen.R`:

```r
test_that("spread bet: opposite-side cell on the same line is exact (not alt)", {
  # Pick: ATH home -0.5 at +120
  # Book: ATH -0.5 (+120) / STL +0.5 (-140)
  # Expected: BOTH the pick row AND the opposite row should report
  # is_exact_line = TRUE because both teams are on the bet's line.
  bets <- make_bet_row(
    market      = "spreads_1st_5_innings",
    line        = -0.5,
    bet_on      = "Athletics",
    market_type = "spreads"
  ) %>% mutate(home_team = "Athletics", away_team = "St. Louis Cardinals")

  book_odds <- list(hoop88 = bind_rows(
    book_row("g1", "spreads", "F5", "Athletics",            -0.5, +120),
    book_row("g1", "spreads", "F5", "St. Louis Cardinals",  +0.5, -140)
  ))

  out <- expand_bets_to_book_prices(bets, book_odds)
  expect_equal(nrow(out), 2)
  expect_setequal(out$side, c("pick", "opposite"))
  expect_true(all(out$is_exact_line),
              info = "both pick and opposite should be exact when book is on the bet's line")
})

test_that("spread bet: opposite-side cell on a different line is alt", {
  # Pick: ATH home -0.5
  # Book: ATH 0 (-117) / STL 0 (+102) -- different line
  # Expected: opposite row reports is_exact_line = FALSE (correctly alt).
  bets <- make_bet_row(
    market      = "spreads_1st_5_innings",
    line        = -0.5,
    bet_on      = "Athletics",
    market_type = "spreads"
  ) %>% mutate(home_team = "Athletics", away_team = "St. Louis Cardinals")

  book_odds <- list(bet105 = bind_rows(
    book_row("g1", "spreads", "F5", "Athletics",             0.0, -117),
    book_row("g1", "spreads", "F5", "St. Louis Cardinals",   0.0, +102)
  ))

  out <- expand_bets_to_book_prices(bets, book_odds)
  expect_equal(nrow(out), 2)
  expect_false(any(out$is_exact_line),
               info = "both pick and opposite are on line 0, not the bet's -0.5")
})
```

- [ ] **Step 2: Run the new tests; confirm the first one fails, second passes**

Run from `Answer Keys/`:
```bash
cd "/Users/callancapitolo/NFLWork/.claude/worktrees/mlb-bets-tab-followups/Answer Keys" && \
  Rscript -e 'testthat::test_file("tests/test_odds_screen.R")'
```

Expected output: the new "spread bet: opposite-side cell on the same line is exact" test fails with `expect_true(all(out$is_exact_line))` reporting one of the two rows has `is_exact_line = FALSE`. The "different line is alt" test passes (it's testing a case the old code already handled correctly — included as a regression guard for the fix).

- [ ] **Step 3: Apply the fix to `expand_bets_to_book_prices`**

Open `Answer Keys/MLB Answer Key/odds_screen.R`. Locate the slot loop around line 199-209:

```r
      for (slot in names(side_labels)) {
        side_value <- side_labels[[slot]]
        if (is.na(side_value)) next

        candidates <- book_frame %>%
          filter(game_id == bet$game_id,
                 market  == bet$market_type,
                 period  == bet$period,
                 side    == side_value)

        chosen <- .pick_closest_line(candidates, bet$line, bet$bet_on)
```

Replace the `chosen <- ...` line (and add the conditional just above it) with:

```r
        # Spread/alt_spread: away-side candidates carry the negated line
        # (see scraper_to_canonical()). The opposite slot must compare
        # against -bet$line, not bet$line, or every same-line book gets
        # flagged as a line mismatch. Totals share one line across both
        # sides; moneyline has line=NA.
        is_spread_market <- bet$market_type %in% c("spreads", "alternate_spreads")
        effective_line <- if (slot == "opposite" && is_spread_market && !is.na(bet$line)) {
          -bet$line
        } else {
          bet$line
        }
        chosen <- .pick_closest_line(candidates, effective_line, bet$bet_on)
```

- [ ] **Step 4: Run the tests again; confirm all green**

Run:
```bash
cd "/Users/callancapitolo/NFLWork/.claude/worktrees/mlb-bets-tab-followups/Answer Keys" && \
  Rscript -e 'testthat::test_file("tests/test_odds_screen.R")'
```

Expected: every `test_that` block passes (the two new ones plus all pre-existing ones — verify "0 failed").

- [ ] **Step 5: Commit**

```bash
cd "/Users/callancapitolo/NFLWork/.claude/worktrees/mlb-bets-tab-followups" && \
  git add "Answer Keys/MLB Answer Key/odds_screen.R" "Answer Keys/tests/test_odds_screen.R" && \
  git commit -m "$(cat <<'EOF'
fix(bets-tab): correct is_exact_line for spread opposite-side cells

The opposite slot in expand_bets_to_book_prices was passing bet$line
(pick-side line, e.g. -0.5) into .pick_closest_line, but the away-side
candidates are emitted with negated lines by scraper_to_canonical (e.g.
+0.5). Result: every same-line book got flagged as a line mismatch on the
opposite row, painting the V8 card grid amber across the entire
opposite-side row.

Pass an effective_line per slot: -bet$line for spread/alt_spread
opposite slot, bet$line everywhere else. Totals (Over/Under share one
line) and moneyline (line=NA) are unaffected.

Co-Authored-By: Claude Opus 4.7 (1M context) <noreply@anthropic.com>
EOF
)"
```

---

## Task 2: Bug 2 — Bet105/BFA alt-total quotes missing from cards

**Why second.** Same file as Task 1 (`odds_screen.R`); reuses the testthat workflow we just exercised. Verification requires re-running the MLB pipeline (Task 4 covers that for all three bugs at once).

**Files:**
- Modify: `Answer Keys/MLB Answer Key/odds_screen.R:31-42` (`.derive_period`, `.derive_market_type`)
- Test:   `Answer Keys/tests/test_odds_screen.R` (append new `test_that` block)

- [ ] **Step 1: Write the failing tests**

Append to `Answer Keys/tests/test_odds_screen.R`:

```r
test_that(".derive_period recognizes _fg/_f3/_f5/_f7 alt-suffix convention", {
  # Internal helpers are unexported. Source the file fresh and call the
  # underscore-prefixed names directly (they're top-level in the same env).
  expect_equal(odds_screen:::.derive_period("alternate_totals_fg"),  "FG")
  expect_equal(odds_screen:::.derive_period("alternate_spreads_f3"), "F3")
  expect_equal(odds_screen:::.derive_period("alternate_totals_f5"),  "F5")
  expect_equal(odds_screen:::.derive_period("alternate_spreads_f7"), "F7")
  # Existing convention still works
  expect_equal(odds_screen:::.derive_period("totals_1st_5_innings"), "F5")
  expect_equal(odds_screen:::.derive_period("totals"),               "FG")
})

test_that(".derive_market_type strips _fg/_f3/_f5/_f7 alt-suffix", {
  expect_equal(odds_screen:::.derive_market_type("alternate_totals_fg"),  "alternate_totals")
  expect_equal(odds_screen:::.derive_market_type("alternate_spreads_f3"), "alternate_spreads")
  expect_equal(odds_screen:::.derive_market_type("alternate_totals_f5"),  "alternate_totals")
  expect_equal(odds_screen:::.derive_market_type("alternate_spreads_f7"), "alternate_spreads")
  # Existing convention still works
  expect_equal(odds_screen:::.derive_market_type("totals_1st_5_innings"), "totals")
  expect_equal(odds_screen:::.derive_market_type("totals"),               "totals")
})

test_that("alt-total bet (e.g. Royals Under 6.5) joins to a book with un-suffixed alt market", {
  # Mimics the real failure mode: bet from compare_alts_to_samples uses
  # 'alternate_totals_fg'; Bet105/BFA scraper writes 'alternate_totals'.
  # After the fix both sides canonicalize to (alternate_totals, FG) and
  # the join produces both pick + opposite rows.
  bets <- make_bet_row(
    market      = "alternate_totals_fg",
    line        = 6.5,
    bet_on      = "Under",
    market_type = "alternate_totals"
  )
  # Force re-derivation by removing market_type so the inside-helper
  # path is exercised the same way mlb_bets_combined would feed it.
  bets$market_type <- NULL
  bets$period <- NULL

  book_odds <- list(bet105 = bind_rows(
    book_row("g1", "alternate_totals", "FG", "Over",  6.5, +200),
    book_row("g1", "alternate_totals", "FG", "Under", 6.5, -240)
  ))

  out <- expand_bets_to_book_prices(bets, book_odds)
  expect_equal(nrow(out), 2)
  expect_setequal(out$side, c("pick", "opposite"))
  expect_true(all(out$is_exact_line))
  expect_equal(out$american_odds[out$side == "pick"], -240)
})
```

> Note on the `odds_screen:::` qualifier: `odds_screen.R` is sourced as a flat script (not a package), so the `:::` form will not work. Replace `odds_screen:::.derive_period(...)` with the bare `.derive_period(...)` call — the `source()` at the top of the test file puts these helpers in the same environment as the test. (If you're reading this as the implementer: just remove the `odds_screen:::` prefix in the test code above before pasting.)

Apply the note inline. The tests should call `.derive_period(...)` and `.derive_market_type(...)` directly, no namespace prefix.

- [ ] **Step 2: Run the new tests; confirm both helper tests fail and the integration test fails**

Run:
```bash
cd "/Users/callancapitolo/NFLWork/.claude/worktrees/mlb-bets-tab-followups/Answer Keys" && \
  Rscript -e 'testthat::test_file("tests/test_odds_screen.R")'
```

Expected: the two `.derive_*` tests fail (existing functions don't recognize `_fg/_f3/_f5/_f7`); the integration test fails because the join can't match.

- [ ] **Step 3: Apply the fix to `.derive_period` and `.derive_market_type`**

Open `Answer Keys/MLB Answer Key/odds_screen.R`. Replace lines 31-42 (the two helper functions) with:

```r
.derive_period <- function(market_name) {
  case_when(
    str_detect(market_name, "_1st_3_innings$") ~ "F3",
    str_detect(market_name, "_1st_5_innings$") ~ "F5",
    str_detect(market_name, "_1st_7_innings$") ~ "F7",
    # Alt-market convention used by compare_alts_to_samples: <market>_<period>
    # (e.g. alternate_totals_fg, alternate_spreads_f5). Recognize all four.
    str_detect(market_name, "_f3$")            ~ "F3",
    str_detect(market_name, "_f5$")            ~ "F5",
    str_detect(market_name, "_f7$")            ~ "F7",
    str_detect(market_name, "_fg$")            ~ "FG",
    TRUE                                       ~ "FG"
  )
}

.derive_market_type <- function(market_name) {
  market_name %>%
    str_replace("_1st_[357]_innings$", "") %>%
    str_replace("_(fg|f3|f5|f7)$", "")
}
```

- [ ] **Step 4: Run the tests again; confirm all green**

Run:
```bash
cd "/Users/callancapitolo/NFLWork/.claude/worktrees/mlb-bets-tab-followups/Answer Keys" && \
  Rscript -e 'testthat::test_file("tests/test_odds_screen.R")'
```

Expected: every `test_that` block in the file passes (existing + the 3 new ones from Tasks 1 + 2). Verify "0 failed."

- [ ] **Step 5: Commit**

```bash
cd "/Users/callancapitolo/NFLWork/.claude/worktrees/mlb-bets-tab-followups" && \
  git add "Answer Keys/MLB Answer Key/odds_screen.R" "Answer Keys/tests/test_odds_screen.R" && \
  git commit -m "$(cat <<'EOF'
fix(bets-tab): canonicalize _fg/_f3/_f5/_f7 alt-market suffix in odds_screen

mlb_bets_combined uses alternate_totals_fg / alternate_spreads_f5 for
alt-line bets (the suffix convention compare_alts_to_samples emits).
Some scrapers (Bet105, BFA) write the un-suffixed alternate_totals /
alternate_spreads. The join in expand_bets_to_book_prices required
exact equality on market_type, so Bet105 alt-totals never joined and
no Bet105 pill rendered on alt-total cards. WZ slipped through only
because it writes alternate_totals_fg literally in its raw DB.

Extend .derive_period() and .derive_market_type() to recognize the
_fg/_f3/_f5/_f7 alt convention alongside the existing _1st_X_innings
recognition. Both sides now canonicalize to (alternate_totals, FG) and
the join matches.

Co-Authored-By: Claude Opus 4.7 (1M context) <noreply@anthropic.com>
EOF
)"
```

---

## Task 3: Bug 1 — game times in viewer-local timezone

**Why third.** UI-only change in a different file (`mlb_dashboard.R`); no existing test scaffolding for the dashboard rendering layer, so verification is manual via dashboard reload (covered in Task 4).

**Files:**
- Modify: `Answer Keys/MLB Dashboard/mlb_dashboard.R:1420-1437` (the `mutate()` block in `create_bets_table()`)
- Modify: `Answer Keys/MLB Dashboard/mlb_dashboard.R` (one-time `<script>` block in the bets-tab page template — see Step 3 for placement)

- [ ] **Step 1: Replace the server-side `tipoff` format with a UTC ISO marker**

Open `Answer Keys/MLB Dashboard/mlb_dashboard.R`. Locate the `mutate()` inside `create_bets_table()` around line 1420-1437:

```r
  table_data <- all_bets %>%
    mutate(
      bet_hash       = pmap_chr(list(id, market, bet_on, line), generate_bet_hash),
      is_placed      = bet_hash %in% placed_hashes,
      matchup        = paste(away_team, "@", home_team),
      tipoff         = ifelse(is.na(pt_start_time), "",
                              format(pt_start_time, "%a %I:%M %p")),
      ev_pct         = ev * 100,
      ...
```

Change the `tipoff` line to emit an ISO+Z UTC string (the browser will convert):

```r
      # tipoff carries a UTC ISO string; client-side JS in the bets-tab
      # template walks every [data-tipoff] element and renders it via
      # toLocaleString() in the viewer's local timezone (matches the
      # parlay-table approach at mlb_dashboard.R:515-524). Server-side
      # format() printed UTC clock time because pt_start_time loses its
      # America/Los_Angeles tag round-tripping through DuckDB.
      tipoff         = ifelse(is.na(pt_start_time), "",
                              format(lubridate::with_tz(pt_start_time, "UTC"),
                                     "%Y-%m-%dT%H:%M:%SZ")),
```

- [ ] **Step 2: Update `render_bet_card` to wrap tipoff in a `data-tipoff` span**

Open `Answer Keys/MLB Dashboard/mlb_dashboard.R`. Locate `render_bet_card()` around line 1252 — specifically the `matchup_line` HTML template at lines 1269-1277:

```r
  matchup_line <- sprintf(
    '<div class="matchup-line">
       <span class="primary">%s</span>
       <span class="sep">&middot;</span>
       <span class="secondary">%s</span>
     </div>',
    htmltools::htmlEscape(matchup),
    htmltools::htmlEscape(tipoff)
  )
```

Replace it with:

```r
  # tipoff is now a UTC ISO string (e.g. "2026-05-15T02:06:00Z"). The
  # bets-tab page-level <script> walks every [data-tipoff] element on
  # DOMContentLoaded and renders it via Date.toLocaleString() in the
  # viewer's local timezone. Pre-render with &nbsp; so layout doesn't
  # shift before JS runs. Empty tipoff (NA pt_start_time) renders as
  # an empty separator block.
  matchup_line <- if (nzchar(tipoff)) {
    sprintf(
      '<div class="matchup-line">
         <span class="primary">%s</span>
         <span class="sep">&middot;</span>
         <span class="secondary" data-tipoff="%s">&nbsp;</span>
       </div>',
      htmltools::htmlEscape(matchup),
      htmltools::htmlEscape(tipoff)
    )
  } else {
    sprintf(
      '<div class="matchup-line">
         <span class="primary">%s</span>
       </div>',
      htmltools::htmlEscape(matchup)
    )
  }
```

The `tipoff` parameter signature on `render_bet_card` doesn't change — it still receives a string. We're just changing what *kind* of string it is (UTC ISO instead of pre-formatted local time) and how it's rendered into HTML.

- [ ] **Step 3: Add the client-side renderer to the bets-tab page**

Find the existing inline `<script>` block at the bottom of the bets-tab page template (the same template that renders the cards — search for an existing `<script>` that runs on DOMContentLoaded inside `mlb_dashboard.R`). Add this initializer alongside the other DOMContentLoaded code:

```html
<script>
  // Render UTC ISO tipoffs in the viewer's local timezone.
  // Mirrors the parlay table renderer at mlb_dashboard.R:515-524.
  function renderTipoffs(root) {
    (root || document).querySelectorAll('[data-tipoff]').forEach(function(el) {
      var iso = el.getAttribute('data-tipoff');
      if (!iso) return;
      var d = new Date(iso);
      if (isNaN(d.getTime())) return;
      el.textContent = d.toLocaleString(undefined, {
        weekday: 'short', hour: 'numeric', minute: '2-digit'
      });
    });
  }
  document.addEventListener('DOMContentLoaded', function() { renderTipoffs(); });
  // Re-render after refresh swaps in new card HTML
  if (window.MLB_AFTER_REFRESH) {
    window.MLB_AFTER_REFRESH.push(renderTipoffs);
  } else {
    window.MLB_AFTER_REFRESH = [renderTipoffs];
  }
</script>
```

If a `MLB_AFTER_REFRESH` hook doesn't exist in this codebase, the second half of the snippet is a no-op (the array variable just sits there). The dashboard's existing `Refresh` button reloads the page, which re-fires DOMContentLoaded, so the initial render covers the common case.

- [ ] **Step 4: Spot-check by reloading the dashboard**

Start the dashboard locally (or rely on the user's running instance) and hard-reload the bets tab. Visually confirm one card's tipoff reads in PT — e.g. a 7:07 PM ET game shows as `Thu 4:07 PM` (PT). If the value still shows UTC clock time, the JS didn't run — open the browser devtools Console and look for `data-tipoff` attributes on the cards and any errors.

- [ ] **Step 5: Commit**

```bash
cd "/Users/callancapitolo/NFLWork/.claude/worktrees/mlb-bets-tab-followups" && \
  git add "Answer Keys/MLB Dashboard/mlb_dashboard.R" && \
  git commit -m "$(cat <<'EOF'
fix(bets-tab): render game times in viewer-local timezone

The V8 card hero strip formatted pt_start_time server-side, but the
column round-trips through DuckDB (mlb_bets_combined) which strips the
America/Los_Angeles tz tag set by format_bets_table. R's DuckDB driver
hands it back as UTC POSIXct, so format() printed UTC clock time
(e.g. a 4:06 PM PT game showed as Thu 11:06 PM).

Match the legacy parlay-table approach (mlb_dashboard.R:515-524):
emit pt_start_time as a UTC ISO string in a data-tipoff attribute,
and add a small DOMContentLoaded script that walks the cards and
renders each via Date.toLocaleString() in the viewer's tz. Works for
any viewer's OS timezone with no server config.

Co-Authored-By: Claude Opus 4.7 (1M context) <noreply@anthropic.com>
EOF
)"
```

---

## Task 4: Pipeline rerun + dashboard verification

**Why fourth.** Bug 2's fix only becomes visible after the MLB pipeline rebuilds `mlb_bets_book_prices` with the new join logic. Bug 3 is visible immediately on dashboard reload (no pipeline needed, but a fresh load is the right verification). Bug 1 is visible immediately. Verifying all three together saves a context switch.

**Files:** none (verification only).

- [ ] **Step 1: Re-run the MLB pipeline so `mlb_bets_book_prices` rebuilds with the new join**

Run from the worktree root:
```bash
cd "/Users/callancapitolo/NFLWork/.claude/worktrees/mlb-bets-tab-followups/Answer Keys" && \
  python run.py mlb 2>&1 | tail -40
```

Expected: pipeline completes without errors. Look in the tail of stdout for the line counts written to `mlb_bets_combined` and `mlb_bets_book_prices`. The number of rows in `mlb_bets_book_prices` should be larger than before (Bet105 alt-total rows now joining).

> Pipeline runs against the worktree's pricing code but writes to the shared `Answer Keys/mlb_mm.duckdb`. That's fine for verification (the dashboard reads from the same file), but be aware the user's other terminal sessions also see the rebuilt data.

- [ ] **Step 2: Confirm Bet105 alt-total rows now exist in the per-book table**

Run:
```bash
cd "/Users/callancapitolo/NFLWork/.claude/worktrees/mlb-bets-tab-followups" && \
  python3 -c "
import duckdb
con = duckdb.connect('Answer Keys/mlb_mm.duckdb', read_only=True)
df = con.execute('''
  SELECT bookmaker, market, period, COUNT(*) AS n
  FROM mlb_bets_book_prices
  WHERE market LIKE 'alternate_%'
  GROUP BY 1,2,3 ORDER BY 1,2,3
''').fetchdf()
print(df.to_string())
"
```

Expected: rows for `bookmaker = 'bet105' AND market = 'alternate_totals'` (and ideally `alternate_spreads`) appear. Before the fix, only `wagerzon` appeared in this query.

- [ ] **Step 3: Hard-reload the dashboard and visually verify all three fixes**

Open the dashboard at http://localhost:8083 (or wherever it's running for you). On the bets tab:

1. **Tipoff** — pick any card; the time should read in PT (e.g. `Thu 4:06 PM`), not UTC.
2. **Bet105 alt pills** — find an alt-total card (e.g. a Royals Under N.5 bet). Confirm Bet105 now shows a price pill alongside Wagerzon. If Bet105 doesn't post on this particular game it'll still be empty; pick a card where Bet105 had a quote.
3. **No spurious opposite-row tags** — pick a spread bet (e.g. ATH -0.5). The opposite row (STL +0.5) cells should NOT have the amber `+0.5` line tag for books that are on the same line. Only books on a genuinely different line (B105 in your earlier screenshot was on `0`) should still show the tag.

If any of the three fail visually, stop and dig into the corresponding task. If all pass, continue to Step 4.

- [ ] **Step 4: Update the Answer Keys CLAUDE.md pitfalls section**

Open `Answer Keys/CLAUDE.md`. Find the "Common Pitfalls" section. Add this bullet immediately after the existing pitfalls list:

```markdown
9. **Alt-market suffix conventions** — `compare_alts_to_samples` writes
   `alternate_totals_fg / alternate_spreads_f5` (suffixed). Some scrapers
   (Bet105, BFA) write the un-suffixed `alternate_totals / alternate_spreads`.
   Both conventions are now canonicalized in `MLB Answer Key/odds_screen.R`'s
   `.derive_period()` and `.derive_market_type()` so the join in
   `expand_bets_to_book_prices` matches across either form. If you add a
   new scraper that uses a third convention, extend those helpers.
```

- [ ] **Step 5: Commit the doc update**

```bash
cd "/Users/callancapitolo/NFLWork/.claude/worktrees/mlb-bets-tab-followups" && \
  git add "Answer Keys/CLAUDE.md" && \
  git commit -m "$(cat <<'EOF'
docs(answer-keys): note dual alt-suffix convention in odds_screen pitfalls

Records the rule that .derive_period() and .derive_market_type() now
recognize both _fg/_f3/_f5/_f7 (compare_alts_to_samples) and
_1st_X_innings (Odds API) suffix conventions, so a future scraper that
uses a third convention is a known extension point.

Co-Authored-By: Claude Opus 4.7 (1M context) <noreply@anthropic.com>
EOF
)"
```

- [ ] **Step 6: Summarize the branch state for the user**

Run:
```bash
cd "/Users/callancapitolo/NFLWork/.claude/worktrees/mlb-bets-tab-followups" && \
  git log --oneline main..HEAD && \
  git diff --stat main..HEAD
```

Expected: 4 commits ahead of main (1 spec + 3 fix commits + 1 doc commit = 5 actually). Diff stat shows `odds_screen.R`, `test_odds_screen.R`, `mlb_dashboard.R`, `Answer Keys/CLAUDE.md`, and the spec doc.

Hand off to the user with: "PR A is ready. 4 fix commits + 1 spec commit on `worktree-mlb-bets-tab-followups`. All testthat tests pass; bets tab visually verified for all three bugs. Ready to merge to main when you give the word."

---

## Out-of-band notes

- **Do NOT merge to main.** The user policy in `CLAUDE.md` requires explicit approval before merging or pushing. After Task 4 completes, stop and wait.
- **Don't write PR B / PR C plans yet.** They wait until PR A merges so the screenshots in those plans match the post-fix dashboard.
- **The test file `test_odds_screen.R` lives in `Answer Keys/tests/`, not `Answer Keys/MLB Answer Key/tests/`.** The `make_bet_row` and `book_row` helpers are defined inline at the top of that file; reuse them.
- **`pt_start_time` round-trip rule** is already documented in `Answer Keys/CLAUDE.md` (Common Pitfalls #7). The doc update we add in Task 4 covers a *different* class of pitfall (suffix canonicalization), not the timezone one — it's worth re-reading #7 if you want context on why Bug 1 happens at all.
