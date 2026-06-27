# RFQ taker — price both teams' margin markets from books

**Date:** 2026-06-27
**Branch:** `worktree-rfq-remove-model` (existing)
**Status:** Design — awaiting review

---

## Review Pack

**What we're building** — In book-only mode the taker currently prices only
the *home* team's run-line correctly. Kalshi lists a separate "Team wins by
over N−0.5 runs" margin market for **both** teams, and the away-team ones get
routed to the wrong book-cache cell (the +1.5 complement instead of the away
−1.5 margin), mispricing ~half of all spread combos by up to ~30pp. This fix
teaches the book cache to hold *both* teams' margin grids and routes every leg
to the correct cell, so all spread×total combos price correctly from books
alone — no model.

**Key decisions**

1. **Select the grid by the line's *sign*, reusing the existing schema** —
   rather than add new cache columns/labels. A negative `spread_line` (−1.5) is
   the home-favorite grid; the positive mirror (+1.5) is the away-favorite grid.
   *Rejected:* a new `favorite_team` column + new combo labels — more churn,
   touches devig + storage, no extra expressiveness. The sign already encodes it.

2. **Fix lives in enumeration + leg-routing, not the scrapers** — all four
   scrapers (DK/FD/PX/NV) were verified sign-agnostic at fetch time; they fetch
   whatever signed line they're handed. The bug is that the *enumeration* only
   ever emits negative target lines. *Rejected:* rewriting scrapers — unnecessary.

3. **No model, anywhere** — the away-margin event was the one thing the model
   priced that books couldn't. We close that gap in the books instead of keeping
   a narrow model. *Rejected:* hybrid (books + a margin-only model) — user
   decision: "the model sucks, we won't use it for now."

4. **Book coverage is the ceiling; thin tails drop gracefully** — a combo no
   book quotes (e.g. away −3.5 longshot) is skipped, never guessed. With no
   model, an un-sourceable price is one we shouldn't bet.

5. **Hold go-live until per-book coverage is measured** — a one-shot coverage
   report (what % of the Kalshi spread×total ladder each book actually returns)
   gates the go-live decision.

**Risks / push back here**

- **Coverage of the dog-side alt ladder is unknown.** away −2.5/−3.5 (a dog
  winning by 3+/4+) may be sparsely quoted. If coverage is poor, the "price
  everything" win is partial in practice. The coverage report makes this
  visible *before* we commit — but if the number is disappointing, you may
  prefer to trade only the lines with solid coverage.
- **Devig validity on the away-favorite grid.** Probit devig assumes a clean
  complementary 4-cell grid. Dog-side alt grids can have wider books / missing
  cells; we keep the existing "need all 4 cells, ≥2 books" floor, which means
  sparse grids self-exclude (safe, but reduces coverage further).
- **This is the last correctness blocker before book-only go-live**, but the
  separate *strategic* risk (book-only makes our fair a public function →
  adverse selection) is unchanged by this fix and still warrants a small,
  measured live start.

**Worth understanding** (opt-in)

- **Complementary devig grid** — the four cells {home covers, away covers} ×
  {over, under} are devigged together so the spread dimension and total
  dimension stay internally consistent (probabilities reconcile). It's like
  normalizing a 2×2 contingency table so margins agree, rather than treating
  four odds independently. The key realization in this design: a single runline
  has *one* such grid, but Kalshi exposes *two* margin markets (one per team)
  that need *two* grids.

---

## Problem & root cause

Kalshi spread markets are **per-team margin markets**, confirmed live:

```
KXMLBSPREAD-...-PHI2   "Philadelphia wins by over 1.5 runs"   = PHI −1.5
KXMLBSPREAD-...-NYM2   "New York M wins by over 1.5 runs"     = NYM −1.5
```

These two are **not** complements — both lose on any 1-run game
(P(PHI−1.5) + P(NYM−1.5) ≈ 0.70, gap = 1-run-margin mass).

The book SGP cache, however, only stores **one** runline grid per game, because
`_fetch_kalshi_spread_lines` (`kalshi_common/sgp_runner.py`) hardcodes the
covering team to `"home"` and computes `line = -(n-0.5)` (always negative),
then dedups by `|line|`. So the cache holds:

```
spread_line = -1.5  →  Home Spread = home −1.5  |  Away Spread = away +1.5 (complement)
```

There is **no cell for "away wins by 2+"** (away −1.5). When the bot trades an
away-team margin leg, `_combo_region_from_legs` routes it to the `Away Spread`
cell (away +1.5, ~80%) when the true market is away −1.5 (~15%). On `main` the
Monte-Carlo model masked this (it priced the real leg and was blended in); in
book-only mode the wrong cell is the unchallenged price.

**Verified facts (code-level):**
- `main` `_load_book_fairs` hardcodes `combo="Home Spread + Over"`; the model
  (`fair_value.model_fair`) was the real per-leg pricer (`main.py:467-479`).
- `mlb_target_lines` contains only negative spreads (`-1.5,-2.5,-3.5`).
- DK (`scraper_draftkings_sgp.py:188`), PX (`scraper_prophetx_sgp.py:419`),
  FD (`scraper_fanduel_sgp.py:505` `(side, signed_line)` key), and NV
  (`scraper_novig_sgp.py:398` `strike == target_spread`) all fetch by signed
  line — none assume home-favorite at fetch time.

---

## Design

### Data model — one grid per team, selected by line sign

Reuse the existing cache schema. `spread_line` is **home-perspective signed**:

```
spread_line  -3.5 -2.5 -1.5      +1.5 +2.5 +3.5
             └─ home-favorite ─┘  └─ away-favorite ─┘
  Home cell  home −1.5 …          home +1.5 …   (home dog covers)
  Away cell  away +1.5 …          away −1.5 …   (away favorite covers = away wins by N+)
```

### Leg → grid+cell routing (the full ladder)

```
Kalshi market          spread_line   cell          book leg
─────────────────────────────────────────────────────────────
home wins by N+  YES  →  -(N-0.5)     Home Spread   home -(N-0.5)
home wins by N+  NO   →  -(N-0.5)     Away Spread   away +(N-0.5)   (complement)
away wins by N+  YES  →  +(N-0.5)     Away Spread   away -(N-0.5)   ← new
away wins by N+  NO   →  +(N-0.5)     Home Spread   home +(N-0.5)   ← new
```

Two rules:
- **Grid (sign):** `spread_line = -(line_n-0.5) if team_is_home else +(line_n-0.5)`
  — set by *whose* margin market the ticker is.
- **Cell (within grid):** existing rule
  `spread_side = "home" if (team_is_home == (side=="yes")) else "away"`.

### Component changes

1. **`kalshi_common/sgp_runner.py::_fetch_kalshi_spread_lines`** — emit the
   home-perspective *signed* line per Kalshi spread market and tag the covering
   team: `line = -(n-0.5) if team_chars == home_code else +(n-0.5)`, append
   `(line, who)`. Drop the hardcoded `"home"` and the dedup-by-`|line|` (it
   collapsed the two teams). `enumerate_kalshi_targets` keeps `_who` only for
   the line value (already signed); it can ignore the tag. Result:
   `mlb_target_lines` now contains both `±(n-0.5)` per game.

2. **`kalshi_mlb_rfq/main.py::_combo_region_from_legs`** — change
   `spread_line = -(spread_leg.line_n - 0.5)` to sign by team:
   `spread_line = -(line_n-0.5) if spread_leg.team_is_home else +(line_n-0.5)`.
   `spread_side` logic unchanged. *(This subsumes the Task-9 quadrant fix; the
   cell label is still built from `spread_side`/`total_side`.)*

3. **`kalshi_mlb_rfq/main.py::_spread_line_from_legs`** — apply the same
   signed-by-team derivation so the candidate loop's `spread_line` matches the
   region's.

4. **`_load_book_fairs`** — no change. Already filters by `spread_line` and
   devigs the quadrant label; it finds the `+1.5` grid once populated.

5. **`kalshi_mlb_rfq/correlation.py` grid_lookup / joint-cell logic** — review
   so same-game joint lookups read the spread cell from the correct (signed)
   grid. The joint of two same-game combos now must agree on which grid each
   leg's spread sits in. Adjust `_grid_lookup` / `joint_prob` cell selection to
   key on the signed `spread_line` already present in `ComboRegion`.

6. **Scrapers (DK/FD/PX/NV)** — **no logic change.** Verify only:
   - FD `price_combo` derives `away_line = -home_line` (symmetric) so a `+1.5`
     target fetches away −1.5 / home +1.5. *(Implementation-time check.)*
   - All four now receive `±` targets; confirm non-empty grids come back live.

### Coverage semantics

A combo is priced iff ≥2 books return a full 4-cell grid for its
`(spread_line, total_line)`. Thin-tail grids (sparse dog-side alts) self-exclude
via the existing floor and are **skipped, not guessed**. Emit a research/log
event (`rejected_no_book_data` / a new `rejected_thin_grid`) so drops are visible.

---

## Testing

**Unit (`tests/test_book_only_pricing.py` / `test_correlation.py`):**
- The 4-row routing table: each `(team_is_home, side, line_n)` → expected
  `(spread_line sign, spread_side)`. Explicitly assert away-YES → `+`,`away`
  and away-NO → `+`,`home`.
- `_combo_region_from_legs` signs the line by team (home → −, away → +).

**Integration (`test_book_only_pricing.py`):**
- Seed `_SGP_ODDS_CACHE` with **both** grids (`±1.5`) at asymmetric odds.
  Assert an away-margin leg prices to the `+1.5` away cell (≈ away −1.5 prob),
  NOT the `−1.5` away +1.5 complement. This is the regression that pins the bug.
- Keep existing suite green (no `USE_MODEL=true` path needed; model is gone in
  book-only, but the flag-guarded tests should still pass).

**Live coverage report (one-shot script, not committed as a scraper):**
- For today's slate: per book, fraction of the Kalshi spread×total ladder
  (both `±` grids, all N) that returns a full 4-cell grid. Output a table.
  This is the **go-live gate** — we read the real coverage number before
  trusting book-only on both grids.

---

## Version control & docs

- **Branch:** continue on existing `worktree-rfq-remove-model` (this worktree).
- **Commits:** (a) enumeration signed lines; (b) leg-routing sign fix + tests;
  (c) correlation grid_lookup; (d) coverage-report script (kept under
  `kalshi_mlb_rfq/` tools, gitignored output). Each with tests green.
- **Docs (same merge):** `kalshi_mlb_rfq/README.md` — book-only now prices both
  teams' margin grids via signed `spread_line`; the root `CLAUDE.md` taker
  bullet (book-only correlation engine line).
- **Pre-merge:** executive-engineer review of `git diff main..HEAD`; **no merge
  to main without explicit approval.** Worktree removed after merge.

## Open implementation checks (resolve during build, not blockers)

- FD `price_combo` symmetric `away_line` derivation (confirm or fix).
- `correlation.joint_prob` Fréchet/same-direction logic still valid when the
  two legs sit in *different* signed grids (cross-grid joints).
- Whether `enumerate_2leg` should suppress redundant same-event candidates
  (e.g. home-NO vs away-YES that map to overlapping cells) — the earlier
  "same-N collision" finding; fold in or defer explicitly.
