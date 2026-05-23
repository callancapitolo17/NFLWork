# MLB Parlay Tab — Books-Agree Count

## Review Pack

**What we're building**

A new "Agree k/n" pill appended to the end of the existing Books cell on the
MLB dashboard's Parlay tab. It counts how many of the five reference fair-
probability voices — DK, FD, ProphetX, Novig, and the model — agree on the
price, where "agree" means "within ±1pp of the median of the five." Display
only. No filter, no color coding, no auto-bet gate. It's a caution signal you
read alongside Edge %.

**Key decisions**

1. **Metric: `k_within_1pp`, median-anchored, model included in the count.**
   Alternatives considered: trimmed range (continuous, harder to set a
   threshold), MAD (collapses to ~0 too often at n=5 and ignores cluster
   size), splitting the model into a separate "model-in-cluster" check
   (statistically cleaner but more UI complexity). The "count of voices
   within ±1pp of the median" wins on interpretability and matches the
   sentence you actually say out loud: *"4 of 5 books agree."*
2. **Band width X = 1pp.** Tight enough to require real agreement. Permissive
   of small calibration noise (sub-percentage-point devig differences across
   books) but not of meaningful disagreement.
3. **Display: B3 — suffix "summary pill" at the end of the Books strip.**
   Same pill shape as the book pills so it integrates with the strip
   visually, but a distinct color (blue tint) so the eye can find it. Label
   text is `Agree k/n` (e.g. `Agree 4/5`, `Agree 3/4`) — the word
   "Agree" lives in the pill so the number is self-explanatory.
4. **Computed at dashboard render time in R, not in the pricer.** All five
   fair-prob columns already live on each parlay row
   (`dk_fair_prob`, `fd_fair_prob`, `px_fair_prob`, `nv_fair_prob`,
   `model_prob_raw`). No schema change to `mlb_parlay_opportunities`, no
   Python pricer change. Cost: a few microseconds per row on each refresh.
5. **NA books reduce the denominator, don't penalize the row.** If a book
   doesn't quote a parlay (NA fair prob), compute k out of however many
   actually quoted and display "k/n" (e.g. `Agree 3/4`). Treating NA as
   "disagree" would over-penalize markets where one book is offline.

**Risks / push back here**

- **Model-as-outlier (scenario H) is invisible.** With "include model" in
  the count, a row where the model is the lone dissenter from a tight book
  cluster reads as `Agree 4/5` — same as a row where a *book* is the lone
  dissenter. You've explicitly accepted this tradeoff; flagging once more
  because it's the one case the metric provably can't distinguish from a
  good-cluster row. Easy mitigation later: layer on a separate
  `model_in_cluster` boolean.
- **No filter / no color in v1.** Display-only means you read the value and
  decide manually. Risk: you forget to glance at it and bet a noisy row at
  the same size as a clean one. Mitigation path: layer a filter input
  ("Min books agree") after a few days of watching live values, once you
  know what threshold actually feels right.
- **Reference set composition is hardcoded at 5.** If a new SGP book scraper
  gets added (or one is dropped), `n` shifts. The code path uses whichever
  fair-prob columns are non-NA so it adapts mechanically, but the UI label
  drifts with it. Worth a code comment noting the dependency so future-you
  isn't surprised.

**Worth understanding** *(opt-in)*

1. **Median-anchored vs mutual-agreement counts.** "Books within ±1pp of the
   *median*" is a one-pass count. "Largest subset where every pair is within
   1pp of every other" is a combinatorial search. They produce the same
   number on every scenario we examined except smooth uniform spread (e.g.
   `42, 44, 46, 48, 50`), where median-anchored gives 3 and mutual gives 2.
   We chose median-anchored because the implementation is a one-liner and
   the difference on that one edge case is academic at n=5.
2. **Why `median` instead of `mean` as the anchor.** A single broken book
   in the cluster can't drag the median far (the median ignores extreme
   values by construction). If we anchored on the mean, the band would shift
   toward the outlier and miss the actual cluster. Same intuition as why
   MAD (median absolute deviation) beats Mean Absolute Deviation in robust
   statistics.

---

## Design body

### Metric definition

Per parlay row, compute `k` and `n`:

```r
compute_k_within <- function(probs, band_pp = 1) {
  probs <- probs[!is.na(probs)]
  n <- length(probs)
  if (n < 2) return(list(k = NA_integer_, n = n))
  med <- median(probs)
  list(k = sum(abs(probs - med) <= band_pp), n = n)
}
```

The input vector is
`c(row$dk_fair_prob, row$fd_fair_prob, row$px_fair_prob, row$nv_fair_prob, row$model_prob_raw)`.
Rows with fewer than two non-NA probs return `NA / n` and the pill is
omitted entirely (rather than showing a meaningless `Agree 1/1`).

### UI: suffix "Agree" summary pill

`render_books_strip()` in `mlb_dashboard.R` already builds the existing
M / DK / FD / PX / NV / Cons pill row. Extend it to accept `k` and `n`
arguments and append one final pill at the end of the strip:

- **Shape**: identical to the existing book pills (same padding, radius,
  font, monospace number).
- **Color**: distinct from the book pills so the eye finds it — light blue
  text (`#8ab4f8`) on a translucent blue background (`#1f6feb22`),
  font weight 700.
- **Spacing**: `margin-left: 8px` to visually separate it from the last
  book pill.
- **Label**: `Agree k/n` (e.g. `Agree 4/5`, `Agree 3/4`).
- **Empty case**: when `n < 2`, omit the pill entirely — no placeholder.

### Files touched

- `Answer Keys/MLB Dashboard/mlb_dashboard.R`
  - Add `compute_k_within()` helper (near the existing parlay-table
    helpers, around the current `apply_combo_residuals()` block).
  - Update `render_books_strip()` to accept `k` and `n` and append the
    suffix pill when both are defined.
  - Update the Books column's `cell` function inside
    `create_parlays_table()` so it computes `k`/`n` from the row and
    passes them into `render_books_strip()`.

No Python pricer changes. No DB schema changes. No new files.

### Worktree, branch, and merge

Already on worktree `parlay-agree-count` at
`/Users/callancapitolo/NFLWork/.claude/worktrees/parlay-agree-count`
on branch `worktree-parlay-agree-count`. Implementation and live-dashboard
verification happen here. Merge to `main` only after explicit approval per
project policy.

### Documentation

After implementation:

- Add one paragraph to `Answer Keys/CLAUDE.md` describing the new pill —
  what the count means, how it's computed, and the explicit "model is
  counted as a voice" tradeoff. Place it near the existing "MLB Dashboard
  — Odds screen" section so anyone scanning that file finds it.
- No README updates needed (no new entrypoint, no new setup steps).
- No memory entry at this stage — the design is small enough that the
  CLAUDE.md mention is enough institutional memory.
