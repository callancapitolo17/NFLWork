# Kalshi MLB RFQ — Two-RFQ Per-Side Kelly Sizing

## Review Pack

**What we're building**
Replace the bot's hardcoded `target_cost_dollars=$1.00` per RFQ with **two
RFQs per combo, each sized via `target_cost_dollars` = (Kelly count × worst-
acceptable price) for its specific side**. Each RFQ is internally tagged
with `intended_side` ("yes" or "no"). At quote-evaluation time, if the
maker's +EV side matches the RFQ's intent, accept; otherwise decline.
Result: every fill lands at-or-near Kelly's count for whichever side
filled, and we naturally pick up more contracts when the maker quotes
better than our worst-acceptable estimate.

**Key decisions**

1. **Two RFQs per combo, side-tagged.** The fundamental constraint of
   Kalshi's RFQ is one size knob per RFQ. Any single-RFQ design is a
   mathematical compromise on at least one side. Two RFQs make the design
   symmetric and correct.

2. **`target_cost_dollars` over `contracts`**. With two RFQs and
   decline-on-mismatch, the opposite-side over-bet concern (which would
   have favored `contracts`) is gone. `target_cost_dollars` then wins on
   the intended side: it captures more contracts when makers quote
   better than our worst-acceptable estimate, while matching exactly when
   they quote at the cliff. Mild over-bet tail risk only at extreme low
   maker quotes (sub-$0.05 — rare fat-finger territory); accepting that
   tail for v1.

3. **Decline-on-side-mismatch**. If RFQ-A is intended for YES side but
   the maker quotes us NO-side EV, we decline (let it walk). The opposite
   side RFQ-B handles that scenario.

4. **Skip cheap RFQs**. If a side's Kelly count is 0 (no acceptable price
   gives 5% EV after fees on that side), don't send that side's RFQ.
   Combos with 0 on both sides skip entirely.

5. **API quota is not actually binding.** Today the bot caps at
   `MAX_LIVE_RFQS=80`. The `_enumerate_and_score_all_games` regularly
   returns 200+ candidates, so the 80-cap is already the bottleneck.
   Halving from "200 considered, 80 fit" to "200 considered, 40 fit"
   moves us from one ceiling to another; coverage doesn't actually drop.

**Risks / push back here**

- **Maker fingerprinting.** Sophisticated MMs will detect "this user
  always sends paired RFQs." They could quote tighter to extract more,
  or quote less competitively because they know we'll take one side
  regardless. Mitigation: randomize the order, add small delays, and
  monitor whether per-maker fill ratios degrade. Belt-and-suspenders: if
  fingerprint cost is real, we can later split the two RFQs across
  different time windows.
- **Schema migration on a live DB.** The bot's `kalshi_mlb_rfq.duckdb` has
  130+ existing fills. The `live_rfqs` ALTER ADD COLUMN is non-destructive
  (NULLs for old rows), but worth backing up before merging.
- **First-cycle behavior.** Bot startup will need to handle "what does
  `intended_side` mean for RFQs created before the deploy?" — likely
  treat NULL as "no constraint" (legacy behavior) so we don't accidentally
  decline mid-flight pre-deploy RFQs.
- **Kelly form conservatism (~10% underbet).** The bot's `kelly.py` uses
  variance form (mu/var); ~10% under-bets vs classical `(p-c)/(1-c)` on
  single bets but generalizes to correlated portfolios. Acceptable for v1;
  worth refactoring later to use classical when no correlated positions
  are involved.

**Worth understanding**

- *Why two RFQs vs. one.* In RFQ markets, the requester must commit
  size before price discovery. Each side of a contract has its own
  Kelly-optimal count, and those counts can differ by 10× at extreme
  fairs. One size knob can't represent two different optimal counts.
  Two RFQs decouples the constraint.
- *Why `target_cost_dollars` instead of `contracts`.* In a two-RFQ design
  with decline-on-mismatch, the opposite-side over-bet risk vanishes — so
  the argument for locking contract count (the main case for `contracts`)
  no longer applies. On the intended side, `target_cost_dollars` captures
  more contracts when the maker quotes better than worst-case, getting us
  closer to ideal Kelly. The mild over-bet tail at penny-level maker
  quotes (sub-$0.05) is rare enough to defer mitigation to v2.
- *Why the Kelly count differs between sides.* Kelly = expected return /
  variance. At fair=0.40, NO has higher win probability (60%), so lower
  variance per contract, so Kelly says you can put more bankroll at risk
  → bigger contract count on NO. Each side's Kelly is computed
  independently from its own (p, ask, fee) triple.

---

## Spec

### Current state (post rebase, commit `2451c0a`)

Already in place from the WIP merge:

- `_outcome_vec_for_legs(samples, typed_legs, side="yes")` — side-aware ✓
- `_load_existing_positions_for_game(game_id, samples)` — uses stored
  `positions.side` ✓
- `_kelly_size_for_quote(quote, fair, side)` — side-aware ✓
- `mint_and_create_rfq(c, target_cost_dollars=N)` — accepts dollar budget
  (will switch to `contracts`) ✓
- `rfq_client.create_rfq` — supports both `contracts` and
  `target_cost_dollars` (mutually exclusive) ✓
- `_refresh_rfqs` skips Kelly==0 candidates ✓ (but only YES side today)

**YES-only legacy from WIP commit (to be reworked):**
- `_worst_acceptable_yes_ask(fair, ev_floor)` — only computes YES side
- `_kelly_size_for_candidate(...)` returns `(int, float)` — only YES
- `_enumerate_and_score_all_games` returns `kelly_sizes` keyed only by YES
- `_refresh_rfqs` sends one RFQ per candidate, sized for YES
- Schema columns `kelly_contracts_at_submit`, `estimated_yes_ask_at_submit`,
  `target_cost_dollars_at_submit` are YES-biased; never reached live DB

### Target design

**`_worst_acceptable_ask(blended_fair, side, ev_floor_pct) -> float`**

Symmetric variant of the YES helper:

```python
def _worst_acceptable_ask(blended_fair, side, ev_floor_pct):
    if side == "yes":
        ev_fn = lambda ask: ev_calc.post_fee_ev_buy_yes(blended_fair, 1 - ask)[1]
        hi = max(0.0001, blended_fair - 1e-6)
    else:
        ev_fn = lambda ask: ev_calc.post_fee_ev_buy_no(blended_fair, 1 - ask)[1]
        hi = max(0.0001, (1 - blended_fair) - 1e-6)
    # binary search in (0.0001, hi)
    return lo if final_ev >= ev_floor_pct else 0.0
```

**`_kelly_size_for_candidate(game_id, typed_legs, samples, blended_fair) -> tuple[int, int, float, float]`**

Returns `(kelly_yes_n, kelly_no_n, worst_yes_ask, worst_no_ask)`.

For each side:
- Compute `worst_<side>_ask` via the symmetric helper above.
- Build `outcome_vec` using `_outcome_vec_for_legs(samples, typed_legs, side=side)`.
- Use side-appropriate `effective_price` and `blended_fair` (flips to
  `1 - fair` for NO).
- Existing positions loaded once via `_load_existing_positions_for_game`
  (already side-aware via stored `positions.side`).

**`_refresh_rfqs` issues up to two RFQs per candidate**:

```python
for c in target:
    yes_n, no_n, yes_ask, no_ask = kelly_sizes[c.leg_set_hash]
    sides_to_send = []
    if yes_n > 0:
        sides_to_send.append(("yes", yes_n, yes_ask))
    if no_n  > 0:
        sides_to_send.append(("no",  no_n,  no_ask))
    for side, n, ask in sides_to_send:
        # skip if we already have a live RFQ for this (leg_set_hash, side)
        if (c.leg_set_hash, side) in live_pairs:
            continue
        target_cost = round(n * ask, 2)
        if target_cost < 0.01:
            continue
        rid, combo_ticker = mint_and_create_rfq(
            c, target_cost_dollars=target_cost, intended_side=side)
        # insert into live_rfqs with intended_side + Kelly state
```

**`_evaluate_quote` side-match gate**:

```python
# Read intended_side from live_rfqs.
intended_side = ...  # from live_rfqs row
# Existing chosen-side logic still picks the +EV side.
if chosen_side != intended_side:
    _log_quote_decision(quote, fair, "declined_side_mismatch",
                         post_fee_ev=chosen_ev_pct, diag=diag)
    return
# Else: proceed with existing accept logic.
```

**Schema additions** (replacing the WIP-commit's YES-biased columns):

```sql
ALTER TABLE live_rfqs ADD COLUMN IF NOT EXISTS intended_side                VARCHAR;
ALTER TABLE live_rfqs ADD COLUMN IF NOT EXISTS kelly_yes_n_at_submit        INTEGER;
ALTER TABLE live_rfqs ADD COLUMN IF NOT EXISTS kelly_no_n_at_submit         INTEGER;
ALTER TABLE live_rfqs ADD COLUMN IF NOT EXISTS worst_yes_ask_at_submit      DOUBLE;
ALTER TABLE live_rfqs ADD COLUMN IF NOT EXISTS worst_no_ask_at_submit       DOUBLE;
ALTER TABLE live_rfqs ADD COLUMN IF NOT EXISTS target_cost_dollars_at_submit DOUBLE;
```

Drop the WIP-commit additions (`kelly_contracts_at_submit`,
`estimated_yes_ask_at_submit`, `target_cost_dollars_at_submit`) from
both `SCHEMA_SQL` and `MIGRATE_SQL` since they never reached a live DB.

Update `combo_cooldown` to be (leg_set_hash, side) PK — already done in
v2 NO-side migration ✓.

### Dedup

Live RFQs are now keyed by `(leg_set_hash, intended_side)`. Update
`live_hashes` lookup in `_refresh_rfqs` to use this composite key:

```python
live_pairs = {(h, s): rid for rid, h, s in con.execute(
    "SELECT rfq_id, leg_set_hash, intended_side "
    "FROM live_rfqs WHERE status='open'"
).fetchall()}
```

### Test plan

- Unit test `_worst_acceptable_ask` for both sides across `fair ∈ {0.1,
  0.295, 0.40, 0.50, 0.55, 0.70, 0.85}` × `ev_floor=0.05`.
- Unit test `_kelly_size_for_candidate` returns the right 4-tuple shape
  for representative fairs.
- Integration test: spin up bot with `--dry-run`, observe one full
  enumerate→score→refresh cycle. Verify `live_rfqs` has up to 2 rows
  per `leg_set_hash` with correct `intended_side` + Kelly counts.
- Manual sanity-check: dump `live_rfqs` after dry-run, eyeball that
  `kelly_yes_n × worst_yes_ask` and `kelly_no_n × worst_no_ask` produce
  reasonable dollar exposures for the bankroll.

### Version control + worktree + docs

- **Branch:** `worktree-rfq-kelly-create-time-sizing` (current). Already
  rebased atop main `b642b41`. The WIP commit `2451c0a` retains
  YES-only sizing post-merge.
- **Commits to add:**
  1. `refactor(kalshi-rfq): symmetric _worst_acceptable_ask + _kelly_size_for_candidate (yes/no)`
     — pure helper refactor.
  2. `feat(kalshi-rfq): two-RFQ per-side Kelly sizing via contracts` —
     wire-up + schema + accept-side gate.
  3. `docs(kalshi-rfq): README sizing section + .env.example`
- **Squash the WIP `2451c0a` into commit 1.** The WIP-commit message is
  misleading post-merge. Rewrite to "refactor(kalshi-rfq): symmetric
  Kelly helpers and migration cleanup".
- **Worktree cleanup:** after merge to main + bot restart, run
  `git worktree remove .claude/worktrees/rfq-kelly-create-time-sizing` +
  `git branch -d worktree-rfq-kelly-create-time-sizing`.
- **Docs updates:**
  - `kalshi_mlb_rfq/README.md` — new "Sizing" section documenting the
    two-RFQ design + the `_worst_acceptable_ask` estimator + the
    side-match accept gate.
  - `kalshi_mlb_rfq/.env.example` — add `KELLY_CREATE_EV_FLOOR_PCT` knob.
  - `kalshi_mlb_rfq/CLAUDE.md` (if exists) — update sizing architecture
    note.
  - MEMORY auto-update post-live (project memory entry for two-RFQ
    sizing once verified in production).

### Rollout

1. **Backup**: `cp kalshi_mlb_rfq.duckdb kalshi_mlb_rfq.duckdb.bak.<ts>`.
2. **Stop bot**: `kill $(pgrep -f "kalshi_mlb_rfq.main")`.
3. **Merge worktree branch to `main`**.
4. **Restart bot**: `./venv/bin/python -u -m kalshi_mlb_rfq.main >> bot.log 2>&1 &`.
5. **Watch first cycle**: look for `rfq_refresh: add=N skipped_kelly_zero=M`
   lines. Should see N > 0 within a few cycles. Each "add" line should
   correspond to a `(combo, side)` pair.
6. **After 10 min**: query `live_rfqs` for new audit columns; spot-check
   that 2 rows per combo where both sides have positive Kelly.
7. **Halt path**: if anything looks off, `kill` bot, `git revert` the
   merge commit on main, restart from previous state. Existing fills
   stay in DB; new audit columns are NULL on legacy rows.

### Out of scope (deliberately deferred)

- **Empirical estimator** using historical quote_log data (replaces
  worst-acceptable-ask with a calibrated median maker quote). Defer
  to v2 once we have 4–8 weeks of quote data.
- **Spread-conditional sizing** (use maker's spread to predict edge).
  Promising but unproven; needs data first.
- **Per-maker calibration** (some makers tighter than others). Same.
- **Cancel-and-reissue probe pattern** for rare big-mispricing
  opportunities. Complex, latency-sensitive; defer.
- **Kelly form unification** (classical formula when no correlated
  positions). Small refactor for ~10% size bump; worth doing but
  separately from this change.
- **Resize-in-flight** when blended_fair shifts between create and
  quote arrival. Defer until we have data on how often this matters.
