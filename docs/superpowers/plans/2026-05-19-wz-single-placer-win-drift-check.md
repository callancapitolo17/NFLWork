# Wagerzon Single Placer — Switch Drift Check from `Odds` to `Win`

**Date:** 2026-05-19
**Branch:** `worktree-fix-wz-single-win-drift-check`
**Author:** Callan + Claude (Opus 4.7)
**Status:** Draft — awaiting approval

---

## Review Pack

**What we're building.** Make `single_placer.place_single` verify Wagerzon's
accepted price by comparing the preflight-returned `Win` (dollars) against an
`expected_win` value sourced from the dashboard, with $0.01 tolerance. Stop
requiring an `Odds` field in WZ's response (WZ omits it for some single-bet
markets like alt-spreads and F-period totals). Round-trip the verified `Win`
from `/api/wz-quote-single` through the client into `/api/place-bet` so the
placer's drift check is `Win`-against-`Win` from the same endpoint — the
exact pattern that's been working in production for parlays for months.

**Key decisions.**

1. **Drift check on `Win`, not `Odds`.** WZ's `ConfirmWagerHelper` reliably
   returns `Win`; `Odds` is per-market-conditional. The parlay placer has
   only ever read `Win`/`Risk`. Rejected alternative: derive `Odds` from
   `Win` + `Risk`. Functionally equivalent math, but a worse fit to the
   existing protocol — adds an indirection where none is needed.
2. **`expected_win` source: dashboard-supplied verified Win first,
   computed math second.** When the user has clicked into Risk and the
   verify-with-wz round-trip succeeded, the dashboard has the authoritative
   WZ `Win` value. We send it as `expected_win`. When no verified Win is
   available (user clicked Place without ever editing Risk), the placer
   computes `expected_win` from `wz_odds_at_place × actual_size` using the
   standard American-odds formula. Mirrors `_resolve_amount_and_win` at
   `mlb_dashboard_server.py:1554-1582` for parlays.
3. **Tolerance: $0.01, same as parlays.** `DRIFT_TOLERANCE_USD` is already
   defined in `parlay_placer.py`. We will import it from there (single
   source of truth) rather than re-declaring. If the math-fallback path
   produces a Win that disagrees with WZ's Win by more than 1 cent due to
   WZ rounding to whole dollars, we'll see it on the first real placement
   and revisit. Rejected alternative: $1.00 tolerance to absorb potential
   WZ rounding. Strictness is intentional — the goal is to catch any real
   drift.
4. **`DRIFT_TOLERANCE_AMERICAN_ODDS = 1` constant goes away.** With the
   Odds-based drift check removed, the constant has no caller. We delete
   it along with the dead `_price_moved(expected, actual)` helper that
   accepted American-odds integers. A new `_price_moved_by_win(expected,
   actual)` helper accepts dollar values for the new check.
5. **The deployment note in `Answer Keys/CLAUDE.md` gets updated.**
   It currently says "endpoint + WT value are inferred from parlay placer
   recon, NOT verified end-to-end against a live WZ account for singles."
   After this fix lands and the user successfully places a real single,
   that caveat goes away.

**Risks / push back here.**

- **Stale `expected_win` on the button dataset.** If the user verifies at
  Risk=$205, then changes Risk to $300 *without* re-verifying, then clicks
  Place — the button would carry the old verified Win for $205. The placer
  would compare WZ's preflight Win for $300 against $205's verified Win and
  reject as `price_moved`. **Mitigation:** clear `data-expected-win`
  whenever Risk changes via edit, reset, or escape. Already a natural
  hook point — the existing `setOverride(card, ...)` calls and the keydown
  Escape handler at `mlb_dashboard.R:5818-5829` are exactly where we wire
  the `delete btn.dataset.expectedWin` calls in. The `verifyWithWz`
  success branch is the only place we *set* it.
- **The fallback math may not perfectly match WZ.** If WZ rounds Win to
  whole dollars in the preflight response, the math-fallback expected_win
  (precise to the cent) will trip the $0.01 drift check. We don't have
  evidence either way yet. **Mitigation:** make `_price_moved_by_win`'s
  error message include both values so a real placement that fails this
  way is debuggable in one glance.
- **JS browser caching.** Dashboard JS is inlined into HTML by R, served
  by Flask. After deploying, the user's browser may still have the old
  `placeBet` and `verifyWithWz` cached. **Mitigation:** hard reload (Cmd+
  Shift+R) is sufficient — no server-side cache busting needed because
  the HTML is regenerated per request.
- **Bet shapes that genuinely don't return `Win` either.** Possible but
  unlikely — parlay_placer assumes `details[0]["Win"]` exists (line 201
  uses subscript, not `.get`, and parlays place successfully today). If
  we hit a single shape where Win is also absent, we'll see it as a
  `missing_win` rejection with a clear key and can add it as a known case.

**Worth understanding.**

- **American odds → decimal payout** (R analogue: think of converting a
  betting line to a multiplier). For positive odds `+215`: $1 risk wins
  $2.15, so decimal odds = 1 + 215/100 = 3.15. For negative odds `-110`:
  $110 risk wins $100, so decimal odds = 1 + 100/110 = 1.909... The
  reverse mapping (Win → Odds) we discussed earlier is exact and bijective
  for fixed Risk; we just don't need it anymore now that the drift check
  lives in `Win` space.
- **What "drift check" actually means.** Imagine you're at a market stall.
  You point at a price tag ("$5.00"), the seller writes it down on a
  ticket and reaches for the item. Drift check = re-confirming the
  $5.00 on the ticket matches the tag before money changes hands. WZ's
  `ConfirmWagerHelper` preflight is the ticket-write step; comparing the
  preflight's Win to the dashboard's verified Win is the re-confirm.
- **Why parlays do this on `Win` instead of `Odds`.** A parlay's combined
  American odds are awkward — multiplying decimal odds and converting back
  is lossy. Win is the direct dollar payout the user is owed; it's the
  cleanest verifiable scalar. The parlay protocol just happens to also be
  the right protocol for singles, which we hadn't noticed until now.

---

## Problem statement

`single_placer.place_single` rejects all single-bet placements where WZ's
`ConfirmWagerHelper` preflight returns a non-empty `details` array but
omits the `"Odds"` key on `details[0]`. We've now observed this on at
least two unrelated market shapes (F7 totals and FG alt-spreads); the
pattern appears to be that WZ does not echo `Odds` on some single-bet
shapes despite returning a valid `Win`. The placer's strict-`Odds`
requirement is at fault — the parlay placer has worked correctly on the
same endpoint for months by reading `Win` instead.

The user-visible symptom is the red toast `Status: rejected — Wagerzon
preflight details missing Odds field` on otherwise valid placement
attempts.

## Approach

Refactor the single-bet placement path so it verifies the accepted price
via `Win` (dollars), exactly mirroring `parlay_placer`. Source the
expected `Win` from the dashboard (verified via `/api/wz-quote-single`,
round-tripped on the place button's dataset) with a math fallback when
no verified value is available.

### Architecture (proposed)

```
                    /api/wz-quote-single        single_pricer.get_single_price()
   user edits      ───────────────────────▶    ─────────────────────────────▶  WZ ConfirmWagerHelper (WT=0)
   Risk field                                                                       │
                                                                                    │ details[0].Win = $441
                                                                                    ▼
                                                          dashboard JS displays "$441 ✓ wz"
                                                          AND stashes 441 on btn.dataset.expectedWin   ◀── NEW
                       
   user clicks      /api/place-bet body now includes
   Place Bet       ───────────────────────▶ expected_win=441                                ◀── NEW
                                                            │
                                                            │ single_placer.place_single
                                                            ▼
                                                      WZ ConfirmWagerHelper (WT=0)
                                                            │
                                                            │ details[0].Win = $441
                                                            ▼
                                                      drift_ok( wz_win=441, expected_win=441, tol=$0.01 ) ✓
                                                            │
                                                            ▼
                                                      proceed to PostWagerMultipleHelper
                                                      (the real bet placement)
```

The Odds field becomes informational only — read if present, but the drift
check is on Win.

### Pseudocode of the new placer block

```python
# wagerzon_odds/single_placer.py — replaces lines 266-276
first_detail = details[0]
wz_win = first_detail.get("Win")
if wz_win is None:
    # WZ should always return Win — parlay placer relies on it. If absent,
    # something genuinely novel happened; refuse to place.
    return _rejected(
        msg="Wagerzon preflight details missing Win field",
        key="missing_win",
    )

expected_win = bet.get("expected_win")
if expected_win is None:
    # User clicked Place without verifying via /api/wz-quote-single.
    # Compute the expected Win from the user's stake and the dashboard-
    # quoted American odds. Mirrors _resolve_amount_and_win for parlays.
    expected_win = _compute_expected_win(
        odds=int(bet["wz_odds_at_place"]),
        risk=float(bet["actual_size"]),
    )

if abs(float(wz_win) - float(expected_win)) > DRIFT_TOLERANCE_USD:
    return _price_moved_by_win(
        expected=float(expected_win),
        actual=float(wz_win),
    )
```

And the helpers:

```python
def _compute_expected_win(odds: int, risk: float) -> float:
    """American-odds → expected Win for a single-leg bet at `risk` stake.
    Matches the math in mlb_dashboard_server.py::_resolve_amount_and_win's
    fallback path so the placer and dashboard agree by construction."""
    decimal = (odds / 100 + 1) if odds > 0 else (100 / -odds + 1)
    return round(risk * (decimal - 1), 2)


def _price_moved_by_win(expected: float, actual: float) -> dict:
    return {
        "status": "price_moved",
        "ticket_number": None,
        "balance_after": None,
        "error_msg": f"Expected ${expected:.2f}, Wagerzon offered ${actual:.2f}",
        "error_msg_key": "drift",
    }
```

---

## Change inventory

### `wagerzon_odds/single_placer.py`

| Line(s) | Action | Description |
|---|---|---|
| 79-82 | DELETE | Remove `DRIFT_TOLERANCE_AMERICAN_ODDS = 1` constant. |
| 151-158 | DELETE | Remove `_price_moved(expected: int, actual: int)` (American-odds version). |
| (new, near top) | ADD | Import `DRIFT_TOLERANCE_USD` from `parlay_placer` so we have one source of truth. |
| (new, near helpers) | ADD | `_compute_expected_win(odds, risk) -> float`. |
| (new, near helpers) | ADD | `_price_moved_by_win(expected, actual) -> dict`. |
| 266-276 | REPLACE | Switch the drift check to `Win`-based per the pseudocode above. |

### `Answer Keys/MLB Dashboard/mlb_dashboard.R`

| Location | Action | Description |
|---|---|---|
| `verifyWithWz` success branch (~line 5771) | EDIT | On `j.win != null`: set `btn.dataset.expectedWin = j.win`. On null: `delete btn.dataset.expectedWin`. |
| `placeBet` body builder (~line 4144) | EDIT | Add `expected_win: data.expectedWin ? parseFloat(data.expectedWin) : null`. |
| Risk reset handlers (~line 5791-5800, 5818-5829) | EDIT | When the user resets or escapes a risk edit, `delete btn.dataset.expectedWin` so we don't send a stale value. |
| Risk-edit-commit handler (`commitEdit` / wherever the Risk override is set) | EDIT | When Risk changes and verification hasn't yet succeeded for the new amount, clear `data-expected-win`. Already implicit if `verifyWithWz` is called on every Risk change — but explicit `delete` on edit-start prevents a brief window where the stale value is sent. |

### `wagerzon_odds/single_pricer.py`

No changes required. The pricer already returns `{win, current_wz_odds, ...}`
and the existing JS reads `j.win`. No new field needed.

### `Answer Keys/MLB Dashboard/mlb_dashboard_server.py`

No changes required. `/api/place-bet` passes the request body straight
through to `single_placer.place_single`, which already accepts a generic
`bet` dict. The new `expected_win` field is picked up by the placer via
`bet.get("expected_win")`.

### `Answer Keys/CLAUDE.md`

Update the deployment note under "WZ single-bet auto-placer":

```
Deployment note: place_single uses the recon-confirmed
PostWagerMultipleHelper.aspx endpoint with WT=0 for singles. Price-match
check is Win-on-Win against an expected_win value supplied by the
dashboard (verified via /api/wz-quote-single) with a $0.01 tolerance —
identical to the parlay placer. Verified end-to-end against a live WZ
account on YYYY-MM-DD via the bets-tab Place button on a CLE -2.5 alt
spread (the bet that originally exposed the Odds-vs-Win bug).
```

(Date filled in once the placement actually succeeds.)

### `wagerzon_odds/CLAUDE.md`

No changes required. The "Quick map" already describes single_placer
generically; the drift-check field choice was never documented at this
level.

---

## Test plan

### New unit tests in `Answer Keys/MLB Dashboard/tests/test_place_bet_dispatch.py`

Or a new file `tests/test_single_placer_win_drift.py` — TBD by what
fits the existing layout best (existing test_place_bet_dispatch.py uses
mocked sessions, which is the right shape; will extend it).

1. **`test_placer_uses_expected_win_when_supplied`** — WZ response has
   `Win=441`; bet payload has `expected_win=441`. → placement proceeds
   to Step 2.
2. **`test_placer_falls_back_to_math_when_expected_win_absent`** — WZ
   response has `Win=440.75`; bet payload has no `expected_win` but has
   `wz_odds_at_place=215`, `actual_size=205`. Math computes
   `_compute_expected_win(215, 205) = 440.75`. → placement proceeds.
3. **`test_placer_returns_price_moved_on_win_drift`** — WZ response has
   `Win=410`; bet payload has `expected_win=441`. → returns
   `price_moved` with error message naming both values.
4. **`test_placer_returns_missing_win_when_win_absent`** — WZ response
   has `details[0]` with no `Win`. → returns `rejected` with
   `error_msg_key="missing_win"`.
5. **`test_placer_no_longer_requires_odds_field`** — WZ response has
   `Win=441` but **no `Odds`** (the case from the user's screenshot).
   With `expected_win=441` in payload, placement proceeds to Step 2.
   This is the regression test for today's bug.
6. **`test_compute_expected_win_positive_odds`** — `(215, 205) → 440.75`.
7. **`test_compute_expected_win_negative_odds`** — `(-110, 100) → 90.91`.
8. **`test_compute_expected_win_even_money`** — `(100, 100) → 100`. Also
   `(-100, 100) → 100`.

### Existing tests to verify still pass

- All 14 tests in `tests/test_place_bet_dispatch.py` (bets pass `idgm`/
  `play` directly and don't depend on the resolver).
- All 14 tests in `tests/test_resolve_wagerzon_play_idgm.py` (from
  yesterday's resolver fix).
- All other dashboard tests (44 total) — should remain green.

### Manual end-to-end test

After the worktree is merged and the dashboard restarted:

1. Open the bets tab, find any WZ pick.
2. Edit Risk to a small amount (e.g. $1 or $5).
3. Observe `✓ wz` checkmark appears in the To Win column with WZ's
   verified Win value.
4. Click Place Bet.
5. **Expected:** "Placed at wagerzon #TICKET_NUMBER" success toast. Bet
   appears in WZ ticket history at the displayed risk + odds.
6. If the bet is rejected, the error message should clearly identify
   either `missing_win` (WZ schema surprise — unlikely) or `drift`
   (price genuinely moved, with both Win values in the message). Not
   `missing_odds` — that error key is being removed.

---

## Version control

- **Branch:** `worktree-fix-wz-single-win-drift-check` (already created
  via EnterWorktree).
- **Worktree dir:** `.claude/worktrees/fix-wz-single-win-drift-check/`.
- **Commit structure:** one commit covering all code changes + tests +
  docs update, since the changes are interdependent (the placer fix
  requires the JS round-trip to actually be useful end-to-end, and
  vice-versa). Single atomic merge to `main` after approval.
- **Suggested commit message:**
  ```
  fix(wagerzon): single_placer drift-checks on Win, not Odds

  Mirrors parlay_placer's Win-on-Win drift check. Round-trips the
  verified Win from /api/wz-quote-single through the place button's
  dataset into /api/place-bet so the placer can compare against WZ's
  preflight Win with $0.01 tolerance — the same protocol that's
  worked in parlay placement for months.

  Removes the strict `Odds` requirement that rejected valid placements
  when WZ omitted that field (observed on F-period totals and FG
  alt-spreads). `Odds` becomes informational; the drift check uses Win,
  which WZ returns reliably.

  Tests: 8 new cases including the regression for "Win present, Odds
  absent" (the user's CLE -2.5 alt-spread bug).
  ```

## Worktree lifecycle

1. Worktree already exists at
   `.claude/worktrees/fix-wz-single-win-drift-check/` on branch
   `worktree-fix-wz-single-win-drift-check`.
2. Implement all changes in the worktree.
3. Run dashboard test suite from the worktree to confirm green.
4. Show diff vs `main`.
5. Get explicit user approval to merge.
6. Merge `--no-ff` into local `main` from the main checkout.
7. Re-run dashboard test suite on `main` post-merge.
8. `ExitWorktree action=remove`.
9. Restart the dashboard (Flask process holds Python imports in memory).
10. User performs the manual end-to-end test.
11. On success, follow up with the `Answer Keys/CLAUDE.md` deployment
    note update (date-filled).

## Documentation

- `Answer Keys/CLAUDE.md` — update the WZ single-bet placer deployment
  note as described above. Done in the same commit as the code changes
  so the doc lands atomically with the behavior change.
- `wagerzon_odds/CLAUDE.md` — no change needed.
- This plan document persists in `docs/superpowers/plans/` as the durable
  record of the design decision (the artifact Callan called out was
  missing from the resolver fix yesterday).

## Rollback plan

If the manual end-to-end test fails in an unexpected way (e.g., `Win`
field also absent for some bet shape, or WZ rounds Win in a way that
breaks the $0.01 tolerance consistently):

1. Don't push to origin until placement actually works.
2. If we've already merged to local `main`, revert with
   `git revert <merge-commit>` — clean, no force-push needed since not
   pushed.
3. Use the captured error message (`missing_win` or `drift` with both
   values) to inform the next attempt. Likely fixes would be widening
   the tolerance, or computing `expected_win` differently.

---

## Open questions (would like Callan's call before implementing)

1. **JS placement of `expected_win` clearing on Risk edit.** The cleanest
   hook is inside `startEdit(valueEl)` at the moment the user opens a
   Risk cell for editing — clear there, then `verifyWithWz` re-sets on
   commit if successful. Sound right?
2. **Test file location.** Add the 8 new tests to existing
   `test_place_bet_dispatch.py`, or new file
   `test_single_placer_win_drift.py`? My slight preference: new file,
   because the existing one is dispatch-focused (which book gets which
   path) while these tests are placer-internals-focused. But happy to
   bundle if you'd rather have fewer files.
3. **Date in the Answer Keys/CLAUDE.md update.** Should I write the
   doc update with `YYYY-MM-DD` as a placeholder and fill in after the
   manual test passes? Or fill in `2026-05-19` now and let it be aspirational
   if the test reveals more work?
