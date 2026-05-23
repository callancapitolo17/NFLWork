# Kalshi RFQ — Enable NO-side Accepts on REST

## Review Pack

**What we're building**
Today the bot only evaluates and accepts the YES side of LP quotes (via the inverted-semantics `accepted_side="no"` → BUY YES). This change adds symmetric NO-side evaluation: per quote we compute EV on both sides, take the higher-EV side if it passes the gate, and record the fill against the correct direction. The REST `accept_quote` endpoint already supports both sides — this is a side-selection change in the bot, not a protocol change. Stays on $1 diagnostic RFQs throughout; real-size sizing waits for the FIX migration.

**Key decisions**

1. **Symmetric EV gates do all the policing — no opposing-side block.** Each side is evaluated independently against current fair after fees. Three regimes fall out of the math:
   - **Same-LP both sides +EV:** mathematically impossible after fees. The LP's spread forces `yes_ask + no_ask > 1`, plus fees forces `yes_ask + no_ask + fees > 1`. For both gates to pass we'd need `1 > yes_ask + no_ask + fees` — contradiction. If it ever does fire it means the fee model or fair is broken; treat as defensive assert (decline both, alert).
   - **Cross-LP both sides +EV in same poll cycle:** real arb. Each gate clearing implies `yes_ask_A + no_ask_B + fees < 1` — combined cost < $1 guaranteed payout. Take both.
   - **Across-time opposing fills:** allowed. The incremental decision when adding NO to a held YES is `E[YES + NO] - E[YES alone] = 1 - p_current - b - fee_b > 0`, which is exactly the NO-side gate condition. So if NO passes its gate against current fair, adding it is strictly better than keeping YES alone. *Alternative rejected:* hard-block once we hold one side. Would forgo the cross-LP arb and refuse beneficial hedges when fair drifts against us.

2. **Schema: add `side` column to `positions` and `combo_cooldown` primary keys.** *Alternative rejected:* track direction implicitly via signed `net_contracts`. The two sides have different `weighted_price` semantics (`yes_ask` vs `no_ask`) so a single signed row can't carry both. Without this change, an upsert into `positions` for the opposite side corrupts the existing row's weighted_price. Backward-compat: existing rows backfill to `side='yes'` (correct — bot has only taken YES historically).

3. **Per-side cooldown.** `combo_cooldown` PK becomes `(leg_set_hash, side)`. After we accept YES on combo X, NO on combo X remains eligible (treated as a genuinely different position in our book).

4. **Diagnostic for hedge formation, not prevention.** When we accept a side opposite to a held position on the same combo, write a `hedge_added` flag to `quote_log` with `original_side`, `original_price`, `current_fair`, `new_price`, and the projected combined net P&L. Lets us monitor whether across-time hedges are happening and whether they're profitable in aggregate, without blocking them.

5. **Stay at `target_cost_dollars=$1.00` on RFQ creation for both sides.** *Alternative rejected:* attempt real-size NO accepts now. Real-size sizing on either side inherits the upstream price-guess problem the FIX migration is designed to fix. Better to gather walk-diagnostic data on both sides at $1 through the FIX build; real size lands when FIX does.

**Risks / push back here**

- **Across-time hedge can be cumulatively -EV even when each leg was +EV at take-time.** Mechanism: take YES at T1 against fair=`p1`, then take NO at T2 against fair=`p2 < p1`. Each leg cleared its gate, but if `p1 - p2 < 2 * fees` the combined position is a small *cumulative* loss relative to having taken neither leg. The forward-looking decision at T2 is still correct (NO improves the held YES position's expected P&L), but the original T1 YES was, in retrospect, a worse bet than we thought. We accept this because (a) the YES was sunk cost by T2, (b) refusing the NO would leave more expected loss on the table than accepting it. The `hedge_added` diagnostic will tell us how often this happens and whether it nets to a real loss pattern. **If after 1-2 weeks of data the hedge-pattern P&L is meaningfully negative, we revisit and consider a "block if fair drift < threshold" rule.**
- **Cross-LP profitable hedge sizing.** When both LPs give +EV on opposite sides in the same poll, we'll take both at $1 each — fine while diagnostic. Once FIX lands and sizing is real, the hedge is genuinely profitable and we *want* to scale both, but Kelly is currently single-position; we'd need to think about how Kelly handles "I'm about to hold both sides" sizing. **Out of scope for this REST patch — flag for the FIX migration plan.**
- **Position reconciliation between accept and DB write.** `get_position_contracts` returns a signed integer. The fill-record path assumes the side it just accepted. If the API call lies or there's a race, we record a position direction that doesn't match Kalshi's truth. Mitigation: cross-check by looking at the position sign (negative = long NO) and warning if it disagrees with the side we sent.
- **Schema migration is one-way.** Adding `side` to `combo_cooldown` and `positions` PKs is a DB change. If we have to roll back, existing rows are still readable but writes from a downgraded binary would clash. Mitigation: keep migration trivial (ALTER TABLE ADD COLUMN with default), backfill, then add to PK.

**Worth understanding**

1. **The LP-spread invariant: `yes_ask + no_ask > 1`.** A Kalshi binary pays $1 to exactly one side. For an LP to make money on a two-sided quote, they must pay less for both sides combined than they collect at settlement: `yes_bid + no_bid < 1`. That forces `yes_ask + no_ask = (1 - no_bid) + (1 - yes_bid) > 1`. After fees, even more so. This is why **same-LP** both sides +EV can't happen — the asks already sum to more than the $1 you'd get back at settlement.

2. **Why cross-LP can still be arb.** Different LPs are independent: nothing forces `yes_ask_A + no_ask_B > 1`. If LP_A is bidding aggressively on NO (so `yes_ask_A` is low) and LP_B is bidding aggressively on YES (so `no_ask_B` is low), their combined asks can dip below $1, and we can lift both sides for a guaranteed profit. Our independent EV gates passing on both sides is *equivalent* to that hedge being profitable after fees — the gates self-police.

3. **Sunk cost reasoning for across-time hedges.** Once we hold YES from a prior cycle, that capital is committed. The forward-looking question is: "given I hold YES, does buying NO now improve expected P&L?" The math says yes whenever NO passes its standalone EV gate against current fair. In R-stats terms, it's analogous to refitting a regression with new data — you don't unwind earlier decisions, you just keep updating with the best current estimate. The original YES might have been a bad bet in retrospect, but that's a separate issue from whether to add NO now.

4. **R analogue for what `accepted_side` does.** Think of `accepted_side` as the LP's offer that we're lifting, not the side we end up holding. In R terms, if `quote <- list(yes_bid = 0.45, no_bid = 0.55)`, calling `accept_quote(side = "no")` means "I'm taking the LP's NO bid" — i.e., the LP buys NO from us at $0.55, which means **we sell NO and keep YES**. The "inversion" in the code comment isn't a bug; it's that the API names the LP's side of the transaction, not ours.

---

## Current state (where we are)

`_evaluate_quote` (`kalshi_mlb_rfq/main.py:767-844`):

```python
no_bid = float(quote.get("no_bid_dollars") or 0)
ev_dollars, ev_pct = ev_calc.post_fee_ev_buy_yes(fair, no_bid)        # YES only
if ev_pct < config.MIN_EV_PCT: ...

contracts = _kelly_size_for_quote(quote, fair)                         # YES-hardcoded
...
resp, err_body = rfq_client.accept_quote(quote["id"], accepted_side="no")   # hardcoded → BUY YES
...
con.execute("INSERT INTO fills (..., side, ...) VALUES (..., 'yes', ...)") # hardcoded
```

The asymmetric pieces:
- EV: only `post_fee_ev_buy_yes` is called; `post_fee_ev_buy_no` already exists in `ev_calc.py:29` but is unused.
- Kelly: `_kelly_size_for_quote` (line 607) computes `effective_price = yes_ask + fee` and uses `outcome_vec = mask.astype(int).values` (combo-hits mask). For NO side, effective_price uses `no_ask = 1 - yes_bid` and outcome flips to `1 - mask`.
- Accept: side hardcoded.
- Fills + positions writes: side hardcoded `'yes'`; `weighted_price` reflects YES price.
- Cooldown: keyed on `leg_set_hash` only; opposite side blocked by accident.
- Inverse gate (`risk.inverse_combo_ok`): checks for *mirror-leg* positions (Yankees ML vs Mets ML), not for *opposite-side* positions on the same combo. Different concept; orthogonal to this change.

## Architecture changes

```
                ┌──────────────────────────────────────────────┐
                │              _evaluate_quote                 │
                │                                              │
                │   ev_yes = post_fee_ev_buy_yes(fair, no_bid) │
                │   ev_no  = post_fee_ev_buy_no (fair, yes_bid)│
                │                                              │
                │   yes_ok = ev_yes_pct >= MIN_EV_PCT          │
                │   no_ok  = ev_no_pct  >= MIN_EV_PCT          │
                │                                              │
                │   if yes_ok and no_ok:                       │
                │     ── defensive math-invariant check ──     │
                │     log MATH_INVARIANT_BROKEN, decline both  │
                │     (same-LP both +EV impossible after fees) │
                │                                              │
                │   chosen = "yes" if yes_ok else              │
                │            "no"  if no_ok  else None         │
                │   if chosen is None: decline_ev              │
                │                                              │
                │   ── hedge-formation diagnostic (non-block) ─│
                │   if positions[combo, opposite(chosen)]:     │
                │     diag.hedge_added = True                  │
                │     diag.original_side / price / fair        │
                │                                              │
                │   contracts = _kelly_size_for_quote(         │
                │       quote, fair, side=chosen)              │
                │   accepted_side = "no" if chosen=="yes"      │
                │                  else "yes"                  │
                │   accept_quote(quote_id, accepted_side)      │
                │                                              │
                │   record fill with side=chosen,              │
                │   price=yes_ask if "yes" else no_ask         │
                └──────────────────────────────────────────────┘
```

## Phase-by-phase

### Phase 1 — Schema migration (~½ day)

`kalshi_mlb_rfq/db.py`:

- Add `side VARCHAR NOT NULL DEFAULT 'yes'` to `positions` and `combo_cooldown`.
- Migrate primary keys:
  - `positions`: PK becomes `(combo_market_ticker, side)`.
  - `combo_cooldown`: PK becomes `(leg_set_hash, side)`.
- Existing rows backfill to `side='yes'` (correct historical interpretation — bot has only taken YES).
- DuckDB migration: `ALTER TABLE ... ADD COLUMN side VARCHAR NOT NULL DEFAULT 'yes'`. Then drop+recreate PK constraint (DuckDB allows this via temp-table swap). Wrap in a `_migrate_v2_side_columns()` function that runs on startup, idempotent.

**Why a migration script, not a fresh schema?** The bot has live state (open RFQs, cooldown entries, position rows) that we don't want to lose on the restart that loads the new code.

### Phase 2 — Symmetric EV + side selection (~½ day)

`kalshi_mlb_rfq/main.py::_evaluate_quote`:

```python
no_bid  = float(quote.get("no_bid_dollars")  or 0)
yes_bid = float(quote.get("yes_bid_dollars") or 0)

ev_yes_$, ev_yes_pct = ev_calc.post_fee_ev_buy_yes(fair, no_bid)
ev_no_$,  ev_no_pct  = ev_calc.post_fee_ev_buy_no (fair, yes_bid)

yes_ok = ev_yes_pct >= config.MIN_EV_PCT and no_bid > 0
no_ok  = ev_no_pct  >= config.MIN_EV_PCT and yes_bid > 0

# Defensive math-invariant check. yes_ask + no_ask + fees > 1 always holds
# for any LP making money on the spread, so both gates passing at once
# means the fee model or fair is broken. Should never fire.
if yes_ok and no_ok:
    log.warning("MATH_INVARIANT_BROKEN: both sides +EV from one LP "
                "yes_ask=%.4f no_ask=%.4f fair=%.4f", 1-no_bid, 1-yes_bid, fair)
    _log_quote_decision(quote, fair, "declined_math_invariant",
                        diag={**diag, "ev_yes_pct": ev_yes_pct,
                              "ev_no_pct": ev_no_pct})
    return

if not (yes_ok or no_ok):
    _log_quote_decision(quote, fair, "declined_ev",
                        diag={**diag, "ev_yes_pct": ev_yes_pct,
                              "ev_no_pct": ev_no_pct})
    return

chosen = "yes" if yes_ok else "no"
diag["chosen_side"] = chosen
diag["ev_yes_pct"]  = ev_yes_pct
diag["ev_no_pct"]   = ev_no_pct
```

Then thread `chosen` through Kelly and accept calls below.

### Phase 3 — Side-aware Kelly + accept + fill record (~½ day)

`_kelly_size_for_quote(quote, fair, side="yes")`:
- For `side="yes"`: existing logic (`effective_price = yes_ask + fee`, `outcome_vec = mask`).
- For `side="no"`: `effective_price = no_ask + fee` where `no_ask = 1 - yes_bid`; `outcome_vec = 1 - mask` (we win when the combo misses).
- Existing-positions handling: for each position's `side`, its `outcome_vec` is the combo-mask or its complement accordingly. (Today all positions are YES; once we accept NO fills the loop must respect each row's side.)

`accept_quote` call:
- `accepted_side = "no" if chosen == "yes" else "yes"` (the LP-side inversion).
- Comment updated to make this mapping obvious.

Fill record (`INSERT INTO fills`):
- `side` = `chosen` ("yes" or "no").
- `price_dollars` = `yes_ask` if chosen=="yes" else `no_ask`.
- `fee_dollars` = `ev_calc.fee_per_contract(effective_price)` for the chosen side.
- `actual` contract count: `get_position_contracts` returns signed int; cross-check that `sign(actual) == +1 for chosen=='yes' else -1`. If sign disagrees, log a `position_direction_mismatch` warning but still write the row (the API is authoritative).

Position upsert (`INSERT INTO positions`):
- Insert with `side=chosen`. PK collision triggers the weighted-price update path — but only against the same-side row. Different-side row is independent.

### Phase 4 — Hedge-formation diagnostic (~¼ day, non-blocking)

Not a gate — a logger. When we're about to accept a side that's opposite a held position on the same combo, write diag fields into `quote_log` so we can monitor whether across-time hedges are happening and whether they're net profitable. Does NOT block the accept.

```python
opposite = "no" if chosen == "yes" else "yes"
with db.connect(read_only=True) as con:
    held = con.execute(
        "SELECT net_contracts, weighted_price FROM positions "
        "WHERE combo_market_ticker=? AND side=? AND net_contracts > 0 LIMIT 1",
        [combo_market_ticker, opposite]
    ).fetchone()
if held:
    held_n, held_price = held
    new_price = (1 - no_bid) if chosen == "yes" else (1 - yes_bid)
    # Projected combined net P&L at settlement (deterministic given both sides
    # held): payout = min(held_n, new_n) of $1 each. For the matched contracts
    # the net is $1 - (held_price + new_price + fees). Surfaces what the hedge
    # actually locks in.
    diag["hedge_added"] = True
    diag["hedge_original_side"]   = opposite
    diag["hedge_original_price"]  = float(held_price)
    diag["hedge_new_price"]       = new_price
    diag["hedge_current_fair"]    = fair
    diag["hedge_projected_net"]   = 1 - held_price - new_price  # before fees
```

Add `hedge_*` columns to `quote_log` schema (Phase 1 migration handles this in the same migration function).

**Why a diagnostic instead of a gate:** the math says any individual +EV accept is forward-looking-correct even when it creates a hedge. The risk we want to *measure* (not prevent) is the cumulative-loss case when fair drift was less than `2 * fees`. After 1-2 weeks of data, if the diagnostic shows hedges netting to a real loss pattern, we revisit with a "block when projected_net < threshold" rule.

### Phase 5 — Cooldown side-awareness (~¼ day)

`risk.cooldown_ok(leg_set_hash, side, cooldown_map)` — change signature to take `side`. `cooldown_map` keyed by `(leg_set_hash, side)`. Callsite in `_all_per_accept_gates_pass` needs the chosen side too, so this gate also moves to the post-side-selection block in `_evaluate_quote`.

`INSERT INTO combo_cooldown` after fill writes `(leg_set_hash, game_id, side, cooled_until, reason)`.

### Phase 6 — Tests + dry-run validation (~½ day)

Unit tests (`tests/kalshi_mlb_rfq/`):

- `test_ev_calc.py`: symmetric `post_fee_ev_buy_no` already exists; add a test that both sides return zero EV when prices are symmetric around fair.
- `test_main_side_selection.py` (new): only-YES eligible → chosen="yes", `accepted_side="no"`.
- `test_main_side_selection.py`: only-NO eligible → chosen="no", `accepted_side="yes"`.
- `test_main_side_selection.py`: neither eligible → "declined_ev".
- `test_main_side_selection.py`: forced both-+EV (mocked quote with `yes_bid + no_bid > 1`, impossible IRL) → "declined_math_invariant", warning logged.
- `test_kelly_no_side.py` (new): NO-side Kelly returns 0 when fair > yes_bid (NO is -EV); returns >0 when fair < yes_bid (NO is +EV).
- `test_main_hedge_diag.py` (new): when held YES exists and incoming NO passes gate, fill proceeds AND `quote_log.hedge_added=True` with correct `hedge_projected_net`. Without held YES, `hedge_added` is absent.

Dry-run validation: start bot with `--dry-run`. Verify in `quote_log`:
- New columns `ev_yes_pct`, `ev_no_pct`, `chosen_side`, `hedge_*` populate correctly.
- `declined_math_invariant` count should be 0 across the dry-run window. Any non-zero count means the fee model or fair-value pipeline has drifted — investigate before merging.
- Existing YES-only quotes still result in "declined_dry_run" with chosen_side="yes".
- Find a quote where the model fair is below yes_bid (NO is +EV) — expect chosen_side="no" and the decline reason still "declined_dry_run".

Then a small live cutover with $1 diagnostic RFQs only:
- Monitor for the first NO-side fill. Expect: fills.side="no", price_dollars≈no_ask, get_position_contracts returns negative.
- Cross-check `positions` table has a row with side="no".

### Phase 7 — Documentation + merge

- `kalshi_mlb_rfq/README.md` — rewrite "Accept semantics" section to cover symmetric side selection, the LP-side inversion, the math-invariant guard, the hedge-formation diagnostic, and the new diag columns. Note that real-size sizing still awaits FIX.
- `kalshi_mlb_rfq/CLAUDE.md` if it exists — note that RFQ evaluation is now symmetric and any future leg-aware logic must respect both sides.
- Memory entry under "Reference": "Kalshi RFQ NO-side accepts" — semantics, schema changes, position-direction convention.

Merge to `main` with `--no-ff` (matches recent merge style). Restart bot. Confirm first NO-side walk-diagnostic row in `quote_log`.

## Risks (engineering specifics)

- **DuckDB PK migration**: ALTER TABLE on PRIMARY KEY may require the temp-table-swap dance in DuckDB. If the migration breaks, the bot fails to start. Mitigation: write the migration as a function with explicit before/after row-count assertions; smoke-test against a copy of the production state DB before merging.
- **Race on direction at hedge-diagnostic lookup time**: between the diag's `positions` read and the `accept_quote` call, a parallel fill on the opposite side could land — meaning we'd miss tagging this fill as `hedge_added`. Single-threaded today via `ACCEPT_LOCK`, but flag for future multi-threading. Since the diagnostic is informational only, a missed tag isn't a correctness issue, just a data-quality one.
- **`weighted_price` semantics on combined-side positions**: with opposite-side fills now allowed (both within poll cycles via cross-LP arb AND across cycles via fair drift), we'll routinely have two `positions` rows on the same combo with different sides. Queries that downstream consumers (dashboard, P&L) run against `positions` need to be reviewed; will mostly want to SUM or aggregate by side. Audit `_kelly_size_for_quote`'s existing-positions loop in particular — it joins `combo_cache` to `positions` and may need an explicit `side` predicate.
- **Kelly cross-leg correlation when holding both sides**: `kelly_size_combo`'s existing-positions branch computes correlation across `outcome_vec`s. If we hold side="yes" on combo X with outcome_vec=`mask`, then evaluating side="no" on the same combo gives outcome_vec=`1 - mask` — perfect negative correlation, which the Kelly math will treat as a hedge. The math probably works; should be verified in unit tests.
- **No backward compatibility with prior cooldown rows**: existing rows post-migration are stamped `side='yes'`. If we ever held NO positions before this change (we didn't), they'd be lost. We didn't, so safe.

## Version control

- **Branch:** `worktree-kalshi-rfq-no-side` (this worktree)
- **Worktree:** `~/NFLWork/.claude/worktrees/kalshi-rfq-no-side`
- **Commit plan** (rough — 6-8 commits):
  1. `db: add side column + PK update to positions/combo_cooldown; add hedge_* cols to quote_log (migration)`
  2. `ev: symmetric side selection in _evaluate_quote + math-invariant guard + diag columns`
  3. `kelly: make _kelly_size_for_quote side-aware (including existing-positions side predicate)`
  4. `risk: cooldown side-aware`
  5. `diag: hedge_added logging when fill creates an opposite-side combined position`
  6. `fills+positions: write actual side and price`
  7. `tests: side-selection, NO Kelly, math-invariant guard, hedge diagnostic`
  8. `docs: README accept-semantics rewrite for symmetric sides`

## Worktree lifecycle

- This plan file lives in `kalshi-rfq-no-side` worktree already.
- After explicit user approval per phase: merge `worktree-kalshi-rfq-no-side` → `main` with `--no-ff`.
- Remove worktree and delete branch immediately after final merge.
- Per CLAUDE.md: never use `git stash` to move work; never amend.

## Documentation updates (in the merge commit)

- `kalshi_mlb_rfq/README.md` — Accept semantics rewrite (new chosen-side decision tree, LP-side inversion table, math-invariant guard, hedge-formation diagnostic, diag columns).
- `kalshi_mlb_rfq/.env.example` — no env changes (sizing knob remains $1).
- Memory: add "Kalshi RFQ NO-side support" under Reference.

## Open questions for you before I start

1. **Migration safety** — OK to run the schema migration against the live state DB on first start of the new binary? The migration is idempotent and additive; existing rows backfill to `side='yes'`. Alternative: stop bot, manually back up `kalshi_mlb_rfq.duckdb`, run migration offline, restart.
2. **Walk-diagnostic SLO** — how long do you want to watch dry-run + small-live before merging? My default would be: 24 hours dry-run, then 48 hours live at $1, then full review of `quote_log` for NO-side activity AND hedge-diagnostic P&L pattern before declaring done.
3. **Hedge P&L review threshold** — at what cumulative `hedge_projected_net` loss would you want the diagnostic to escalate to a real block-gate? I'd suggest reviewing after 2 weeks of live data and deciding then, rather than picking a threshold now.
