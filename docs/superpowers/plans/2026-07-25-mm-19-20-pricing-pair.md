# Plan: #19 uncertainty-scaled margin + #20 z-space dispersion gate (Session F pricing pair)

## Review Pack

**What we're building** — Two coupled pricing fixes for the Kalshi MLB maker.
#19: the per-side margin becomes `max(ROI part, MIN_MARGIN_PTS + K_SIGMA·σ_books)`
with `(1+r)^N` compounding for N-game combos, so longshot and multi-game quotes
stop carrying near-zero absolute cushion. #20 (rescoped 2026-07-25): the
consensus gate stops tossing outlier books against an absolute ±2¢ band and
instead declines to quote when cross-book dispersion in z-space exceeds
`SIGMA_Z_MAX` — fair = median of ALL books, no outlier removal.

**Key decisions**
1. **Dispersion gate replaces outlier band** (user decision 2026-07-25): an
   outlier book may be the informed one; deleting it and quoting the stale
   majority is adverse selection by construction. False declines are ~free.
2. **`(1+r)^N` compounding** over additive per-game points: matches how every
   sportsbook builds multi-leg parlay hold (per-leg vig multiplies); reduces to
   current behavior exactly at N=1.
3. **σ propagation for cross-game combos**: relative variances add —
   `σ_combo = fair_combo · sqrt(Σ (σ_g/f_g)²)`. Standard error propagation for
   a product of independent estimates.
4. **Margin components go to the research firehose** (`quote_priced` payload),
   NOT into `quote_decisions` — the #14 report's column semantics
   (blended_fair/yes_bid/no_bid) stay untouched.
5. **Knobs on first principles** (`MIN_MARGIN_PTS=0.01`, `K_SIGMA=1.0`,
   `SIGMA_Z_MAX=0.07`): #13 markout data will never exist; calibration comes
   later from settlement P&L (#12) and the #14 demand curve.

**Risks / push back here** — Wider longshot/multi-game quotes may mean fewer
fills at the tails (0 fills to date already); this is accepted — an unfilled
wide quote is free, a filled thin tail quote is the leak. The single-outlier
case now kills the whole combo's quote instead of being outvoted (accepted:
at today's 2–3 live books the outvote rarely existed anyway).

**Worth understanding** — Error propagation for products (relative errors add
in quadrature, like R's `sqrt(sum((s/x)^2))`); sample vs population stddev
(`ddof=1`, R's default `sd()`).

## Design

### #19 — pricing.py
- `quote(fair, target_roi, *, sigma_pts=0.0, n_games=1, min_margin_pts=0.0,
  k_sigma=0.0) -> Quote | None`. Defaults reproduce current behavior exactly
  (N=1, floor 0): old `p/(1+r)` ≡ `p − p·(1−(1+r)^{-1})`.
- Per side (p = side win-prob; σ and floor are side-symmetric):
  - `roi_pts = p · (1 − (1+target_roi)^(−n_games))`
  - `floor_pts = min_margin_pts + k_sigma · sigma_pts`
  - `margin_pts = max(roi_pts, floor_pts)`; target all-in cost `raw = p − margin_pts`
  - fee iteration + grid floor-down + step-down guard unchanged, now
    guaranteeing `bid + fee ≤ p − margin_pts` (realized margin ≥ target).
- `Quote` gains frozen optional component fields (defaults None) for firehose
  logging: `sigma_pts, n_games, floor_pts, roi_pts_yes, roi_pts_no,
  margin_pts_yes, margin_pts_no`. Constructor calls with just bids still work.
- Existing guards kept: `0 < fair < 1`, both bids > 0, `yes_bid + no_bid < 1`.

### #20 — router.py + main.py
- `router.Consensus` frozen dataclass: `fair, sigma_pts, sigma_z, n_books`.
- `router.consensus(book_fairs, min_books, sigma_z_max) ->
  tuple[Consensus | None, str]` with reason ∈ {"ok", "too_few_books",
  "consensus_dispersion"}:
  1. `n < min_books` → (None, "too_few_books").
  2. `fair = median(all)`; `z_i = norm.ppf(clip(f_i, 1e−6, 1−1e−6))`;
     `sigma_z = stdev(z_i, ddof=1)`; `sigma_pts = stdev(f_i, ddof=1)`.
  3. `sigma_z > sigma_z_max` → (None, "consensus_dispersion").
  4. Else (Consensus(...), "ok"). No outlier removal; single median.
- `consensus_detail` updated in lockstep (returns fair + ALL books on pass).
- `subcombo_fair` → returns `(Consensus | None, reason)`;
  `combo_fair_detail(...) -> tuple[ComboFair | None, reason]` where
  `ComboFair = (fair, sigma_pts, n_games)`; cross-game: fairs multiply,
  σ via relative-variance propagation. `combo_fair(...)` stays float|None
  (thin wrapper) for the confirm last-look and risk-sweep call sites.
- `main._consensus_filter`: rewritten as the dispersion-gate oracle — returns
  all books when `n ≥ MIN_AGREEING_BOOKS` and `σ_z ≤ SIGMA_Z_MAX`, else {}.
- Discovery tick: `combo_fair_detail` → distinct skip reasons
  (`too_few_books` / `consensus_dispersion`; other failures stay `no_fair`);
  threads `sigma_pts`/`n_games` + config knobs into `pricing.quote`; adds
  margin components to the `quote_priced` firehose payload.

### config.py
- REMOVE `BOOK_CONSENSUS_BAND`. ADD `SIGMA_Z_MAX=0.07`, `MIN_MARGIN_PTS=0.01`,
  `K_SIGMA=1.0`, each with a rationale comment (continuity derivation for
  0.07; favorite-longshot/options-MM rationale for the floor).

## Tests (TDD: new tests written and shown failing BEFORE implementation)
- **#20 must-fail-first**: 3 books at p≈0.08 with one 2¢ away → declined
  (`consensus_dispersion`); 2 books 2¢ apart at p≈0.08 → declined. Both PASS
  under current code's band → must fail red first.
- **#20 continuity**: 2- and 3-book sets near p=0.50 (3¢-gap pass, 5¢-gap
  fail) behave identically to the old band.
- **#20 degenerate**: 1 book → too_few_books; 2 identical books → quoted;
  oracle (`_consensus_filter`) and production (`router.consensus`) agree on
  shared cases; 2–3-book realism (FD/BetMGM/ProphetX only).
- **#19**: at p=0.10, σ=0.008 → margin ≥ 0.018 (hand-computed floor); at
  p=0.50, σ≈0, N=1 → bids identical to current code; 3-game margin > 1-game
  at same p; property sweep over p×σ×N grid: bids ≥ 0 or None,
  yes+no < 1, realized margin ≥ max(roi_pts, floor_pts).
- Full suite: `python -m pytest kalshi_mlb_mm/tests/ kalshi_common/tests/`.

## Version control / worktree
- Worktree `mm-19-20-pricing-pair`, branch `worktree-mm-19-20-pricing-pair`
  off LOCAL main `45fb8f0` (verified; reset from origin-main base).
- Commits: (1) #19 pricing + router threading + tests; (2) #20 gate + tests;
  (3) docs (README + this plan). Push branch, draft PR closing #19 + #20.
- NO merge without explicit user approval. No pushes to main. Live maker
  restart is user-owned (single post-merge reminder required — live pricing
  path changes).

## Documentation
- `kalshi_mlb_mm/README.md`: rewrite pricing section (formula + components +
  interplay: moderate σ widens margin, σ_z past threshold kills the quote);
  knob table: add SIGMA_Z_MAX/MIN_MARGIN_PTS/K_SIGMA, remove
  BOOK_CONSENSUS_BAND, reword MIN_AGREEING_BOOKS ("books with a fair", no
  longer "band survivors").
- Root `CLAUDE.md` maker bullet: consensus described as z-space dispersion
  gate (one line).
