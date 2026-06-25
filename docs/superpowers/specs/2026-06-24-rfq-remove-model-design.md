# Kalshi MLB RFQ taker — book-only pricing + book-implied correlation

## Review Pack

**What we're building** — Remove the Monte-Carlo model from the Kalshi MLB RFQ
taker entirely (behind a `USE_MODEL` flag, default `false`) and replace its two
jobs with book data:

1. **Pricing** → fair value = `median(books)` (the model leaves the median).
2. **Risk / correlation sizing** → a **book-implied correlation engine** built
   from the per-game spread×total SGP grid the bot already scrapes, feeding
   Kelly's existing-positions branch so correlated same-game combos can't be
   overbet.

This fixes the fill collapse (the 600s model-staleness gate increasingly
rejected quotes as the external R model's refresh cadence degraded from ~7 min
in May to ~23 min on June 3, when 95% of quotes were declined) *and* preserves
the one thing the model was genuinely good at — correlation-aware sizing — using
the books instead of the sim.

**Key decisions**

1. **Book-only via a flag, not a rip-out.** `USE_MODEL` (default `false`); model
   code stays dormant. Reversible in one env var; lets us A/B later; lowest risk.
   Rejected full deletion.
2. **Keep the ≥2-book floor for pricing.** Fair = `median(books)` with the
   existing `MIN_BOOK_COUNT_FOR_BLEND=2` gate. Removing the model doesn't change
   which combos have book data, so coverage is unchanged. Rejected ≥3 (cuts
   volume) and "require a sharp book" (more code).
3. **Drop the prediction-staleness gate when the flag is off.** The model is no
   longer a pricing input, so its freshness is irrelevant. Live books stay
   guarded by `MAX_BOOK_STALENESS_SEC`.
4. **Correlation from the book SGP grid, not the model or a blunt cap.** Compute
   pairwise combo correlation exactly from grid lookups + interval algebra, with
   a conservative ρ=1 fallback when the needed grid cell isn't quoted by ≥2
   books. Rejected: model-for-correlation (keeps the stale dependency we're
   removing), ρ=1-everywhere / one-position-per-game (too blunt — the user
   explicitly wants a real correlation estimate), and on-demand union SGP
   scraping (needs scraper changes + accept-path latency).

**Risks / push back here**

- **Grid coverage of correlation.** The exact method covers spread/total legs
  (the core combos). Cross-category MVE legs not in the spread×total grid
  (player props, etc.) fall back to ρ=1 — safe but less sizing-efficient. If a
  lot of your held positions are cross-category, much of the "real correlation"
  benefit degrades to the conservative cap. Worth confirming the live position
  mix is mostly spread/total.
- **Activating dormant code.** Kelly's existing-positions branch is currently
  dormant ("revisit before relying on it"). This change makes it live for the
  first time, driven by book-implied covariance. New live sizing math → needs
  careful tests and a small-stake validation before ramp.
- **Pricing off two soft books.** With N=2 allowed and no model anchor, a combo
  quoted by only two soft books is priced off their average. EV gate + per-game
  cap limit damage; you accepted this for volume.
- **Model regression left unaddressed.** This routes around the degraded R
  pipeline rather than fixing it (out of scope).

**Worth understanding** (opt-in)

- **The SGP grid is already a joint distribution.** Each SGP price like
  `P(home covers −1.5 AND over 8.5)` is a *joint* probability of two correlated
  legs. The whole grid is a coarse picture of how the game's margin and total
  move together — which is exactly what you need for correlation, and why no
  separate model is required. (In R terms: it's an empirical 2-D contingency
  table over (margin-region, total-region).)
- **Why correlation reduces to a grid lookup.** Two combos are both
  step-functions of the same two numbers (margin M, total T). The event "both
  hit" is just "M and T land in the *tighter* of the two regions" — which is
  itself another spread×total combo, i.e. another grid cell. So the joint you
  need is usually a value you already scraped.

---

## Design body

### Background: why fills collapsed

The taker accepts maker quotes that are +EV vs its fair value, today
`median(model_fair, book₁, book₂, …)`. A per-accept gate rejects the quote if the
model's samples are older than `MAX_PREDICTION_STALENESS_SEC` (600s).

Evidence (`kalshi_mlb_rfq_research.duckdb`, `quote_log`):

- `declined_stale_predictions` is the #1 decline reason (53,446 in `quote_log`;
  66% of research-firehose gate evaluations).
- Its share tracks the model's refresh cadence degrading:

  | Day | model median refresh gap | % refreshes > 10 min | accepts |
  |-----|--------------------------|----------------------|---------|
  | May 19 | 7.0 min | 11% | high |
  | May 30 | 7.1 min | 32% | 33 |
  | May 31 | 15.2 min | 71% | 18 |
  | Jun 02 | 18.9 min | 76% | 17 |
  | Jun 03 | 23.4 min | 95% | 0 |

- A prior fix set `MAX_STALENESS_SEC=900` in `.env`, but the code reads
  `MAX_PREDICTION_STALENESS_SEC` — a different name — so it was a silent no-op.

### Part 1 — Book-only pricing (`USE_MODEL` flag)

The model feeds the trading path in six places. `USE_MODEL=false` changes each:

| # | Site (`main.py`) | Today (`USE_MODEL=true`) | Book-only (`USE_MODEL=false`) |
|---|---|---|---|
| 1 | `_fresh_blended_fair` (~450) — quote-eval pricing | `median(model, books)`; drop if no samples | `median(books)`; no sample load |
| 2 | candidate loop (~1452) — RFQ creation | `continue` if no samples; `blend(model, books)` | `median(books)`; never skip on samples |
| 3 | staleness gate (~486) | `declined_stale_predictions` if model > 600s | gate **skipped** |
| 4 | `_kelly_size_for_quote` (~620) | `return 0` if no samples; `outcome_vec` from samples | book-implied correlation Kelly (Part 2) |
| 5 | `_kelly_size_for_candidate` (~717) | same as #4 | same as #4 |
| 6 | `_refresh_caches` (~153) — loads `mlb_game_samples` + meta | always loads | skip sample/meta load |

`_load_book_fairs` already returns `{}` unless ≥2 books price the exact
(game, spread, total) tuple, so dropping the model yields `median(books)` over
the same candidate set. New helper:

```python
def _book_only_fair(book_fairs: dict[str, float]) -> float | None:
    vals = [v for v in book_fairs.values() if v is not None]
    return statistics.median(vals) if vals else None
```

`config.py`:
```python
USE_MODEL = _get("USE_MODEL", "false").lower() in ("1", "true", "yes")
```

### Part 2 — Book-implied correlation engine

**Purpose.** When sizing a new combo on a game where we already hold positions,
estimate the correlation between the new combo and each held combo, so Kelly
down-sizes correlated additions instead of treating them as independent.

**The covariance Kelly needs.** Kelly's existing-positions branch
(`kelly.kelly_size_combo`) uses, for each held position, the covariance of
*returns* between the new bet and that position:

```
cov_returns(new, pos) = [ P(new ∩ pos) − P(new)·P(pos) ] / (price_new · price_pos)
```

`P(new)` and `P(pos)` are the book-only fairs we already compute. The only new
quantity is the **joint** `P(new ∩ pos)`.

**Computing the joint from the grid.** Every combo is `(M-region) ∩ (T-region)`
where M = home margin, T = total. For two combos, the joint event is the
intersection of their M-regions and T-regions — itself a spread×total region:

```
A = (home −1.5, Over 8.5)   → M > 1.5, T > 8.5
B = (home −2.5, Over 9.5)   → M > 2.5, T > 9.5
A ∩ B                        → M > 2.5, T > 9.5  = grid cell (home −2.5, Over 9.5)
P(A ∩ B) = book SGP fair of (home −2.5, Over 9.5)        ← direct grid lookup
```

- **Same-direction pairs** (both Over, both Under, nested spreads): the joint is
  the SGP fair of the (tighter-spread, tighter-total) cell — one lookup.
- **Crossing pairs** (one Over + one Under, or opposing spread sides): the joint
  is computed by inclusion-exclusion over grid cells, e.g.
  `P(M>s, t₁<T<t₂) = P(M>s, T>t₁) − P(M>s, T>t₂)`, each term a grid lookup. A
  combo on incompatible regions (e.g. Over 9.5 ∩ Under 8.5) has joint 0.
- All lookups go through the same `MIN_BOOK_COUNT_FOR_BLEND=2`, `median`-across-
  books path as `_load_book_fairs`, so the joint is a book *consensus*.

**Fallback (ρ=1, conservative).** If the needed grid cell isn't priced by ≥2
books, or either leg is outside the spread×total grid (cross-category MVE legs:
player props, etc.), use the comonotonic bound `P(A∩B) = min(P(A), P(B))`. This
is the maximum possible joint given the marginals → maximum positive covariance
→ maximum down-sizing. Safe by construction; only sizing-efficiency is lost.

**New module** `kalshi_mlb_rfq/correlation.py` (pure functions, unit-testable):

```python
def joint_prob(legs_a, legs_b, grid_lookup) -> float | None:
    """P(A ∩ B) from the spread×total grid via region intersection +
    inclusion-exclusion. None if not resolvable from the grid."""

def cov_returns(p_a, p_b, p_joint, price_a, price_b) -> float:
    """Return-space covariance for Kelly."""
```

`grid_lookup(spread_line, total_line, quadrant) -> float | None` wraps the
existing `_SGP_ODDS_CACHE` median-across-books read.

**Kelly refactor.** `kelly.kelly_size_combo`'s existing-positions branch
currently derives covariance from model sample paths (`outcome_vec`). Replace
that inner computation with `cov_returns(...)` built from book-implied joints;
`mu_new`/`var_new`/`base_frac` already come from `blended_fair` and are
unchanged. The function signature drops its reliance on `outcome_vec` for the
covariance term (kept as an ignored/optional arg for `USE_MODEL=true` back-compat,
or branched on the flag). Existing positions are loaded from the `positions`
table (leg-sets via `combo_cache`); each position's `P(pos)` and `price_pos` are
recomputed from current book fairs for consistency with the new bet.

**Backstops retained.** `per_game_cap_ok` (per-game dollar cap) stays as an
independent second limit, so even if correlation is under-estimated the game's
aggregate exposure is bounded.

### Testing

`correlation.py` (pure unit tests — no DB):
- Same-direction joint = tighter-cell lookup (mocked grid).
- Crossing joint via inclusion-exclusion; incompatible regions → 0.
- Missing grid cell or non-grid leg → ρ=1 fallback `min(P(A),P(B))`.
- `cov_returns` sign/magnitude on hand-computed cases.

Kelly:
- Two perfectly-correlated held + new combos size like ~one bet (not double).
- Uncorrelated (independent grid joint = P(A)·P(B)) reproduces independent Kelly.
- `existing_positions=[]` → single-bet Kelly off `blended_fair`.

Pricing / gates:
- `USE_MODEL=false`: `_all_per_accept_gates_pass` never returns
  `declined_stale_predictions` even with `_SAMPLES_META_GENERATED_AT=None`.
- `USE_MODEL=false`: `_fresh_blended_fair` = `median(books)` without samples.
- Regression: `USE_MODEL=true` path + existing tests still pass (audit tests
  that assume the staleness gate runs by default, since the runtime default
  flips to `false`).

Pre-go-live verification (no live trading): full `tests/` suite green, then a
single-cycle / dry-run smoke from the worktree confirming the bot enumerates
candidates, prices book-only, computes correlation without raising, reaches the
accept path, and shows no `declined_stale_predictions` in `quote_log`. **User
approval required before turning live**, with small-stake validation of the
newly-activated correlation sizing before any ramp.

### Version control

- Worktree/branch `worktree-rfq-remove-model` off `main`.
- Files: `kalshi_mlb_rfq/config.py`, `kalshi_mlb_rfq/main.py`,
  `kalshi_mlb_rfq/correlation.py` (new), `kalshi_mlb_rfq/kelly.py`,
  `kalshi_mlb_rfq/tests/` (new tests), `kalshi_mlb_rfq/.env.example`
  (fix dead `MAX_STALENESS_SEC` → `MAX_PREDICTION_STALENESS_SEC`; add
  `USE_MODEL`), `kalshi_mlb_rfq/README.md`, root `CLAUDE.md`.
- Also fix the live `.env` `MAX_STALENESS_SEC` no-op as part of go-live config.
- Commits: (1) `USE_MODEL` flag + book-only pricing + gate/cache changes + tests;
  (2) `correlation.py` + Kelly refactor + tests; (3) docs. Pre-merge executive
  review of `git diff main..HEAD`; merge only after user approval.

### Documentation

- `kalshi_mlb_rfq/README.md`: `USE_MODEL`, book-only pricing, the book-implied
  correlation engine + ρ=1 fallback, staleness gate now model-only.
- Root `CLAUDE.md`: taker bullet → book-only pricing default + book-implied
  correlation sizing.
- `.env.example`: `USE_MODEL=false` + corrected staleness var name.

### Out of scope

- Fixing the R model pipeline cadence / whole-number-Over calibration leak.
- On-demand union SGP scraping for correlation (grid + ρ=1 fallback instead).
- Any change to line enumeration, RFQ submission, or accept mechanics.
