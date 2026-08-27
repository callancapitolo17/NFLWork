# Issue #98 — router reads the leg surface for cross-game combos

**Branch:** `feature/mm-surface-router`
**Worktree:** `.claude/worktrees/issue-98-surface-router` (created before any file)
**Base:** `main` @ e2791ba (#96 merged). Epic #94. #99 (age gate) comes AFTER.

---

## 1. The one rule

A game group with **exactly one leg** prices from the in-memory `LegSurface`.
Everything else is untouched.

| per-game group | route today | route after #98 |
|---|---|---|
| 1 leg (only possible inside a cross-game combo) | live on-demand fetch | **leg surface, zero I/O** |
| 2-leg grid (spread×total, ml×total) | live on-demand fetch | live on-demand fetch |
| `on_demand` (3-leg, F5 mix, …) | live on-demand fetch | live on-demand fetch |
| `unpriceable` | decline | decline |

`classify_subcombo(gl) == "single"` is exactly "1 leg for that game", and it
already runs the F5-TIE guard, the duplicate-market guard and the #86
contradiction guard *before* the `n == 1` branch — so a lone TIE leg still
classifies `unpriceable`, not `single`. Lone single-leg RFQs are out of scope at
`_priceable` (`len(canon) < 2`), so a `single` group only ever appears inside a
multi-game combo. That is why "1 leg" and "cross-game" are the same predicate,
and why same-game never touches the surface.

## 2. Files changed

### `kalshi_mlb_mm/router.py`

- New public predicate:
  ```python
  def routes_to_surface(game_legs) -> bool:
      return legset.classify_subcombo(game_legs) == "single"
  ```
  One definition, read by `subcombo_consensus`, by `combo_fair_detail`, and by
  the discovery tick's feed block — so the pricer and the fetch-queuer can never
  disagree about which games cost network.
- `subcombo_consensus(..., surface_fairs=None)`: under `live_routing=True`, a
  `routes_to_surface` group with a non-None `surface_fairs` reads
  `surface_fairs(game_legs) -> {book: fair}` and runs the **unchanged** #20
  gate. A failing gate returns `surface_too_few_books` / `surface_dispersion` —
  renamed only so the decline counts stay readable next to `live_too_few_books`
  (constraint 3: no silent behaviour, everything countable).
- `combo_fair_detail(..., surface_fairs=None)`: calls `resolve_game` **only for
  groups that are not surface-routed**.

  *Why bypass it (issue note 1):* under `live_routing=True` the resolved
  `game_id` is never read — `subcombo_consensus` keys the live path on
  `leg_set_hash`, not `game_id` — so for a surface group `resolve_game` is a
  gate with no consumer, and `"unresolved_game"` would be a pure false decline.
  The surface is keyed on `CanonicalLeg.game_id` (the Kalshi suffix), which the
  legs already carry. `resolve_game` stays for every live-routed group, byte for
  byte.

  *What is NOT bypassed:* `main._discovery_tick` keeps its own `no_game` decline
  on `game_ids_list`. That resolve feeds the per-game exposure cap, the tipoff
  gate and the `fill_games` ledger — it is a **risk** gate, not a pricing one,
  and it stays. (Its team-name `LIMIT 1` doubleheader bug is real and explicitly
  out of scope here.)
- `combo_fair(...)` gains the same pass-through.

`surface_fairs is None` reproduces today's routing exactly — a `single` group
falls back to the on-demand path. That is the off switch.

### `kalshi_mlb_mm/config.py`

- `SURFACE_ENABLED` (see the open question in §6).

### `kalshi_mlb_mm/main.py`

1. Module global `_SURFACE_INGEST` / `_SURFACE` alongside `_ENGINE`.
2. `_surface_fairs(game_legs) -> dict[book, float]` — the injected lookup.
   Reads `LegSurface.book_fairs(leg)` (in memory). **Never** opens
   `kalshi_mlb_mm_surface.duckdb`: a connect is ~17ms and has caused three
   incidents in hot loops. `book_fairs` already collapses FanDuel's two routes
   to one `"fanduel"` key — not re-split, because two routes at one book are one
   opinion and splitting them would let FD satisfy `MIN_AGREEING_BOOKS` alone.
3. **Discovery feed block** (the acceptance-critical change): a group where
   `router.routes_to_surface(gl)` and the surface is live is skipped entirely —
   no `ensure_fetch`, no `od_pending`, no `landed_empty`. This is what makes a
   cross-game RFQ cost zero outbound requests. A surface **miss** is not a live
   fallback: the book simply is not in that leg's `book_fairs`, so it is
   excluded from consensus, exactly as the issue specifies.
4. **Confirm tick**: surface groups are skipped when building `refetch_jobs`,
   *before* `_resolve_game_for_legs` / `_game_ref` — otherwise an unresolvable
   game would void a fill we can price. A pure cross-game combo ends with an
   empty `refetch_jobs` and `refetch_ok = True`, and `combo_fair` re-prices it
   from the surface's current rows. The last look itself (`risk.last_look_ok`,
   the #17 singles veto) is untouched.
5. **Post-fill cooldown** (`_post_fill_live_refresh_landed`): a surface group
   has no flight to wait on, so the analogous condition is
   `>= MIN_AGREEING_BOOKS` books holding a row for that leg with
   `built_at > filled_at` — "the books were re-asked after the pick-off". Same
   fail-closed posture. `_ensure_post_fill_fetches` skips surface groups (there
   is nothing to queue; the ingest loop is already refreshing them).
6. **Risk sweep** `_current_consensus_fair` and the drift/constituent-jump path:
   pass the same lookup through. The breakers themselves are unchanged.
7. **Research**: `quote_priced` gains an additive `surface_games` key parallel
   to `live_games` — per game, per book: `fair`, `route`, and `age_sec`
   (`now - built_at`). This is the loud, countable placeholder constraint 3
   asks for, and it is precisely the dataset #99 needs to pick
   `SURFACE_MAX_AGE_SEC`. `live_games` is untouched so `report.py` keeps reading.
8. `main_loop`: construct `SurfaceIngest`, `.start()` it, `.stop()` in the
   `finally` — one construction, mirroring `OnDemandEngine`.

### Tests — `kalshi_mlb_mm/tests/test_surface_routing.py` (new)

- Cross-game 1+1: `_ENGINE` replaced by a counting fake → **`ensure_fetch`
  called 0 times**, quote priced. (Acceptance 1, asserted on a counter.)
- Same-game 2-leg grid: fetch count unchanged, surface never consulted — a
  populated surface is present but must not be read. (Acceptance 2, mirroring
  `test_live_pricing.py`'s "cache present but never consulted".)
- Mixed 1+2: exactly one `ensure_fetch`, for game B only; game A from surface.
  (Acceptance 3.)
- Surface with 1 book → `surface_too_few_books`; wide dispersion →
  `surface_dispersion`. Gate semantics identical to live.
- FanDuel present on both routes counts as **one** book.
- `surface_fairs=None` → routing byte-identical to today (regression oracle).

Plus the full existing suite, which pins same-game (`test_live_pricing.py`,
`test_router_on_demand.py`, `test_quorum_quoting.py`, `test_cooldown_refresh.py`).

## 3. Live verification (not optional — #96 shipped six bugs, four live-only)

1. `python -m kalshi_mlb_mm.leg_surface` for ~10 min against real books; confirm
   `mlb_leg_surface` fills and per-book pass logs are clean.
2. A read-only harness script (scratchpad, deleted after) that: starts a real
   `SurfaceIngest`, pulls real open `KXMLB*` RFQ legs from Kalshi, and prices
   them through `router.combo_fair_detail` with a counting `_ENGINE` — proving
   zero outbound book requests on real cross-game leg sets and reporting the
   surface hit rate and `built_at` age distribution.
3. Bot itself stays **off** (died 2026-08-19) unless you say otherwise.

## 4. Version control

- Commits: (1) router routing + predicate, (2) main wiring + config, (3) tests,
  (4) docs. `git diff --stat` before each.
- Pre-merge review against the CLAUDE.md checklist on `git diff main..HEAD`,
  findings presented as ISSUES vs ACCEPTABLE RISKS, then **explicit approval**
  before merging to `main`.
- After merge: `git worktree remove` + `git branch -d`.
- No `.duckdb` copied or symlinked into the worktree; the surface DB is written
  fresh by the ingest loop in the worktree's own directory.

## 5. Documentation (same merge)

- `kalshi_mlb_mm/README.md` — "Leg surface" status flips from *shipped dark* to
  *wired for cross-game*; new routing table; the `SURFACE_ENABLED` knob; the
  `surface_games` research key; the note that same-game is unchanged.
- `CLAUDE.md` — the maker blurb currently says live-only is the ONLY pricing
  path and that the surface is dark until #98. Both need the qualification.
- `kalshi_mlb_mm/research_queries.sql` — a surface-age / surface-vs-live decline
  query for #99.

## 6. Open question

`SURFACE_ENABLED` default. See the question asked alongside this plan.
