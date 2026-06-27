## Task 7 Report: Generic mapping + fail-closed validation

### Status
DONE

### TDD Evidence

**RED:** Wrote `unabated_edge/tests/test_mapping.py` first. Running it immediately produced:
```
ImportError: cannot import name 'mapping' from 'unabated_edge'
```
(1 error, 0 collected — confirmed RED before writing implementation.)

**GREEN:** Created `unabated_edge/mapping.py` with `Pairing`, `pair_events`, and `validate`. Both tests passed immediately:
```
test_validate_rejects_incomplete PASSED
test_validate_accepts_complete   PASSED
2 passed in 1.34s
```

### Files Changed
- `unabated_edge/mapping.py` — created (35 lines)
- `unabated_edge/tests/test_mapping.py` — created (16 lines)

### Commit
`b4c7a5d` — `feat(unabated-edge): generic Unabated↔Kalshi mapping + validation gate`

### Full Suite
11/11 passed across all test files:
- test_config, test_feed_deltas, test_feed_snapshot, test_kalshi_venue, test_mapping (×2), test_pricing (×3), test_soccer_adapter (×2)

### Self-Review

**`validate` is fail-closed:**
- Returns `False` if `set(fair) != set(adapter.outcomes)` (wrong outcome keys)
- Returns `False` if any outcome ticker is missing or falsy (`not all(...)`)
- Returns `False` if `|sum(fair) - 1.0| >= 1e-6` (probabilities don't sum to 1)
- Only returns `True` when all three conditions pass

**`_named` resolution works correctly:**
- `adapter.map_outcome_tickers(kev)` returns a dict that may have a `"_named"` key containing `[(name, ticker), ...]`
- `tickers.pop("_named", [])` removes the sentinel and iterates it
- `named = {adapter.canon_team(n): t ...}` builds a canon-name → ticker lookup
- The frozenset key is built from `named` keys (the canon team names), so Unabated events are matched by canonical team pair
- `ot["home"]` and `ot["away"]` are resolved by looking up `adapter.canon_team(ev.home)` and `adapter.canon_team(ev.away)` in that `named` dict

**Sport-agnostic:**
- No soccer literals anywhere in `mapping.py`
- All sport-specific logic flows through the `adapter` argument: `adapter.canon_team`, `adapter.map_outcome_tickers`, `adapter.outcomes`
- Works for any `SportAdapter` subclass

### Concerns
None.

---

## Fix: case-mismatch bug in pair_events team-name join

### Bug
`map_outcome_tickers()` called `.lower()` on `yes_sub_title` before appending to `_named`. This caused `canon_team("argentina") == "argentina"` on the Kalshi side, while the Unabated feed sent `"Argentina"` → `canon_team("Argentina") == "Argentina"`. The frozenset keys never matched, so `pair_events` silently returned `[]` for all non-aliased teams.

### Test added (RED → GREEN)

**Test command:**
```
python3 -m pytest unabated_edge/tests/test_mapping.py unabated_edge/tests/test_soccer_adapter.py -v
```

**Before fix (FAILING):**
```
FAILED unabated_edge/tests/test_mapping.py::test_pair_events_non_aliased_teams
AssertionError: Expected 1 Pairing, got 0
```

**After fix (PASSING):**
```
unabated_edge/tests/test_mapping.py::test_validate_rejects_incomplete PASSED
unabated_edge/tests/test_mapping.py::test_validate_accepts_complete   PASSED
unabated_edge/tests/test_mapping.py::test_pair_events_non_aliased_teams PASSED
unabated_edge/tests/test_soccer_adapter.py::test_canon_alias          PASSED
unabated_edge/tests/test_soccer_adapter.py::test_fair_three_way_sums_to_one PASSED
5 passed in 1.21s
```

Full suite: 12/12 passed (no regressions).

### Files changed
- `unabated_edge/sports/soccer.py` — preserve raw-stripped name in `_named`; lowercase only a local var for draw/tie detection
- `unabated_edge/tests/test_mapping.py` — added `test_pair_events_non_aliased_teams` (Argentina vs Austria)

### Commit
`6bd2bf0` — `fix(unabated-edge): stop lowercasing Kalshi team names before canon join`
