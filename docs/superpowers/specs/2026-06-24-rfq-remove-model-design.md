# Kalshi MLB RFQ taker — remove the model from pricing (book-only mode)

## Review Pack

**What we're building** — A `USE_MODEL` config flag (default `false`) for the
Kalshi MLB RFQ taker bot. When off, the bot prices and sizes every combo purely
from live book consensus (median of the sportsbooks that quote it), and the
model-prediction staleness gate is skipped entirely. This restores fills, which
collapsed because the bot's fixed 600s model-staleness gate increasingly
rejected quotes as the external R model pipeline's refresh cadence degraded from
~7 min (May) to ~23 min (June). The flag is reversible: set `USE_MODEL=true` and
the bot behaves exactly as it does today.

**Key decisions**

1. **Book-only via a config flag, not a rip-out.** Decision: add `USE_MODEL`
   (default `false`), leave the model code in place but dormant. Rejected: full
   deletion of `model_fair`/sample-loading. Why: reversible in one env var, lets
   us A/B the model later if its cadence/calibration is fixed, and is far lower
   risk to ship tonight.
2. **Keep the ≥2-book floor unchanged.** Decision: fair value = `median(books)`
   with the existing `MIN_BOOK_COUNT_FOR_BLEND=2` gate. Rejected: raise to ≥3
   (cuts fillable volume) or require a sharp book present (more code). Why:
   removing the model doesn't change which combos have book data, so the
   candidate set is essentially unchanged; DK/FD price SGPs with correlation
   built in, so 2 sharp books is a reasonable consensus.
3. **Drop the prediction-staleness gate when the flag is off.** Decision: skip
   `declined_stale_predictions` entirely in book-only mode. Why: the model is no
   longer a pricing input, so its freshness is irrelevant. Live book inputs are
   still protected by the separate `MAX_BOOK_STALENESS_SEC` gate.
4. **Kelly sizes each combo independently in book-only mode.** Decision: pass
   `existing_positions=[]` so Kelly uses the single-bet path off `blended_fair`.
   Why: the model's only unique job — cross-position correlation sizing — is
   already dormant in v1 (no positions branch is "revisit before relying on
   it"). The `per_game_cap_ok` dollar cap remains the backstop against
   over-loading a single game.

**Risks / push back here**

- **Pricing off two soft books.** With the model gone and N=2 allowed, a combo
  quoted by only two soft/correlated books is priced off their average with no
  sharp anchor. You chose to accept this for volume; the EV gate (`MIN_EV_PCT`)
  and per-game cap limit the damage, but it is a real exposure. If realized P&L
  on N=2 combos looks bad, the cheapest mitigation is decision #2's "require a
  sharp book present" variant.
- **Losing model correlation sizing.** Independent Kelly across correlated
  same-game combos can stack exposure faster than full-Kelly-with-correlation
  would allow. Backstopped by `per_game_cap_ok`, but if you ever ramp bankroll,
  revisit this before relying on the cap alone.
- **The model regression is left unaddressed.** This change routes *around* the
  degraded R pipeline rather than fixing it. The pipeline slowdown (and the
  separate whole-number-Over calibration leak) still exist; book-only mode just
  no longer depends on them. Fixing the model is out of scope here.

**Worth understanding** (opt-in)

- **Feature flag.** `USE_MODEL` is a boolean read from the environment that
  switches a code path at runtime — like an `if (config$use_model)` branch in R
  that you toggle without editing code. The value of a flag (vs. deleting code)
  is that the old behavior stays one setting away, so you can compare the two
  empirically instead of arguing about them.
- **Median as a robust consensus.** `median()` of book fairs ignores how far an
  outlier is, only its rank — so one stale/buggy book can't drag the price the
  way a `mean()` would. At exactly two sources median equals mean; the
  robustness only kicks in at N≥3. This is why decision #2's floor matters.

---

## Design body

### Background: why fills collapsed

The taker sends RFQs on Kalshi MVE combos and accepts maker quotes that are
+EV versus its fair value. Fair value today is
`median(model_fair, book₁, book₂, …)`. Before any quote is accepted, a per-accept
gate checks that the model's samples are fresh:

```
gen_at = latest mlb_samples_meta.generated_at
if age(gen_at) > MAX_PREDICTION_STALENESS_SEC (600s):  decline "declined_stale_predictions"
```

Evidence from `kalshi_mlb_rfq_research.duckdb` and `quote_log`:

- `declined_stale_predictions` is the #1 decline reason overall (53,446 in
  `quote_log`; 66% of all gate evaluations in the research firehose).
- Its share climbs in lockstep with the model's refresh cadence degrading:

  | Day | model median refresh gap | % refreshes > 10 min | accepts |
  |-----|--------------------------|----------------------|---------|
  | May 19 | 7.0 min | 11% | high |
  | May 30 | 7.1 min | 32% | 33 |
  | May 31 | 15.2 min | 71% | 18 |
  | Jun 02 | 18.9 min | 76% | 17 |
  | Jun 03 | 23.4 min | 95% | 0 |

- A prior fix attempt set `MAX_STALENESS_SEC=900` in `.env`, but the code reads
  `MAX_PREDICTION_STALENESS_SEC` — a different name — so the override was a
  silent no-op. Effective threshold stayed 600s.

The model samples are pre-game Monte-Carlo simulations that barely change over
10–25 minutes pre-game; the 600s threshold was an over-tight proxy coupled to an
assumption ("model refreshes < 10 min") that is now false. Rather than re-tune
the threshold, the user chose to remove the model from pricing entirely.

### What "remove the model" means precisely

The model feeds the trading path in exactly six places. `USE_MODEL=false`
changes each as follows:

| # | Site (`main.py`) | Today (`USE_MODEL=true`) | Book-only (`USE_MODEL=false`) |
|---|---|---|---|
| 1 | `_fresh_blended_fair` (~450) — quote-eval pricing | load samples; drop if none; `median(model, books)` | `median(books)`; no sample load |
| 2 | candidate loop (~1452) — RFQ creation | `continue` (skip game) if no samples; `blend(model, books)` | price `median(books)`; never skip on samples |
| 3 | staleness gate (~486) | `declined_stale_predictions` if model > 600s | gate **skipped** |
| 4 | `_kelly_size_for_quote` (~620) | `return 0` if no samples; `outcome_vec` from samples | single-bet Kelly off `blended_fair`; `existing=[]` |
| 5 | `_kelly_size_for_candidate` (~717) | same as #4 | same as #4 |
| 6 | `_refresh_caches` (~153) — loads `mlb_game_samples` + `mlb_samples_meta` | always loads from `mlb.duckdb` | skip sample/meta load — decouples from model DB |

Why each is safe without the model:

- **Fair value (#1, #2):** `_load_book_fairs` already returns `{}` unless
  `MIN_BOOK_COUNT_FOR_BLEND` (=2) books price the exact (game, spread, total)
  tuple. So today a candidate is only priced when ≥2 books quote it — the model
  is the 3rd+ voice in the median, not a coverage requirement. Dropping it
  yields `median(books)` over the same candidate set.
- **Staleness (#3):** the gate exists to detect a dead model pipeline. With the
  model out of pricing, it is irrelevant. The `MAX_BOOK_STALENESS_SEC` gate (and
  the `_load_book_fairs` freshness implicit in the SGP cache cadence) still
  guard the live inputs.
- **Kelly (#4, #5):** `kelly.kelly_size_combo` derives `mu_new`/`var_new`/
  `base_frac` entirely from `blended_fair`. The `outcome_vec` argument is used
  *only* in the existing-positions correlation branch, which the code documents
  as "Dormant in v1 (no positions yet); revisit before relying on it." Passing
  `existing_positions=[]` takes the single-bet path; correctness of the live
  path is unchanged.
- **Cache (#6):** skipping the sample/meta load when the flag is off means the
  bot no longer reads `mlb.duckdb` at all for pricing, removing the external
  dependency that caused the regression. (Schedule/commence-time still come from
  the bot's own market DB, unaffected.)

### Implementation sketch

`config.py`:
```python
USE_MODEL = _get("USE_MODEL", "false").lower() in ("1", "true", "yes")
```

New helper in `main.py`:
```python
def _book_only_fair(book_fairs: dict[str, float]) -> float | None:
    """Median of book fairs, or None if empty. The >=2 floor is enforced
    upstream in _load_book_fairs."""
    vals = [v for v in book_fairs.values() if v is not None]
    return statistics.median(vals) if vals else None
```

`_fresh_blended_fair` and the candidate loop branch on `config.USE_MODEL`: when
true, today's path; when false, skip sample loading and use `_book_only_fair`.
The staleness gate is wrapped in `if config.USE_MODEL:`. Both Kelly functions
branch: when false, skip the `samples`/`outcome_vec`/`existing` work and call
`kelly_size_combo(outcome_vec=np.array([]), existing_positions=[], …)`.
`_refresh_caches` skips the `mlb_game_samples` / `mlb_samples_meta` queries when
the flag is off (samples cache stays empty, meta stays `None`).

### Testing

New tests (pytest, alongside existing `tests/`):

- `_book_only_fair` returns the median of supplied book fairs; `None` on empty.
- With `USE_MODEL=false`, `_all_per_accept_gates_pass` does **not** return
  `declined_stale_predictions` even when `_SAMPLES_META_GENERATED_AT` is `None`
  or far in the past.
- With `USE_MODEL=false`, `_fresh_blended_fair` returns `median(books)` and does
  not require `_SAMPLES_CACHE` to be populated.
- Book-only Kelly: `kelly_size_combo` with `existing_positions=[]` returns a
  positive size for a +EV book-only fair and 0 for a -EV one.
- Regression: with `USE_MODEL=true`, all existing model-path tests still pass
  (audit existing tests for any that assume the staleness gate runs by default,
  since the runtime default flips to `false`).

Verification before go-live (no live trading): run the full `tests/` suite, then
a dry-run / single-cycle smoke of the bot from the worktree to confirm it
enumerates candidates, prices them book-only, and reaches the accept path
without raising — and that `quote_log` decisions no longer show
`declined_stale_predictions`. User approval required before turning it live.

### Version control

- Worktree `worktree-rfq-remove-model`, branch `worktree-rfq-remove-model`, off
  `main`.
- Files modified: `kalshi_mlb_rfq/config.py`, `kalshi_mlb_rfq/main.py`,
  `kalshi_mlb_rfq/tests/` (new test file), `kalshi_mlb_rfq/.env.example`
  (rename the dead `MAX_STALENESS_SEC` line → `MAX_PREDICTION_STALENESS_SEC` and
  add `USE_MODEL`), `kalshi_mlb_rfq/README.md`, root `CLAUDE.md`.
- Also fix the live `.env` (not committed) `MAX_STALENESS_SEC` no-op line as part
  of go-live config.
- Commits: (1) flag + pricing/gate/Kelly/cache changes + tests, (2) docs. Single
  feature branch; pre-merge executive review of `git diff main..HEAD`; merge only
  after user approval.

### Documentation

- `kalshi_mlb_rfq/README.md`: document `USE_MODEL`, book-only pricing path, and
  that the staleness gate is model-only.
- Root `CLAUDE.md`: update the taker bullet to note book-only pricing is the
  default (model optional, off by default).
- `.env.example`: `USE_MODEL=false` + the corrected staleness var name.

### Out of scope

- Fixing the R model pipeline cadence / whole-number-Over calibration leak.
- Re-enabling correlation-aware Kelly without the model.
- Any change to line enumeration, RFQ submission, or accept mechanics.
