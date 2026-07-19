# WC Maker Touch-Join Pricing — Design

## Review Pack

**What we're building** — The WC totals maker currently rests quotes ~3c behind
the crowd's 1c spread and has zero fills in 4 days live. This change makes the
engine quote *at* the crowd's best bid — one-sided, only when our devigged
Unabated fair says that price carries ≥1c of edge after maker fees — so our
orders sit where 99.99% of the flow actually trades. Ships tonight for the
ESP-ARG final (July 19); it is also the intended default pricing for the
next-sport maker.

**Key decisions**
1. **Touch-join is default behavior, no enable flag.** Rejected: an opt-in
   `TOUCH_JOIN_ENABLED` config. You chose the long-term path as the default;
   rollback is `git revert` + restart, not a flag.
2. **Overlay, not rewrite.** The legacy fair−margin quote remains the fallback
   whenever the join condition fails. Rejected: replacing legacy quoting
   entirely — bigger diff through tested live code on a one-day deadline.
3. **Queue hysteresis holds resting touch orders while edge ≥ 0.25c.**
   Rejected: requoting on every fair wiggle — queue position is the asset
   (Polymarket's best community bot cut cancellations 95% for this reason).
4. **Mid-jump circuit breaker and markout/toxicity brake deferred** to
   GitHub issues #1 and #2 (your scope call). Tonight's news protection is the
   existing trigger stack (anchor staleness, kickoff pull, fill burst,
   watchdog).
5. **Sizing and caps unchanged** — your uncapped %-caps stand; the goal-grid
   worst-case ledger remains the loss bound.

**Risks / push back here**
- **Fills are now expected, at your uncapped sizing (~$150/side/rung).** The
  replay found no adverse convergence when anchor and crowd disagree, but the
  bet is that our devig is right to within ~1c. If the crowd's shade is
  informed rather than retail bias, we pay for that lesson at full size.
- **No crowd-keyed news pull tonight** (deferred issue #1): a mid-match-eve
  news gap that outruns the 5s poll is bounded only by the ledger's worst case.
- **The fill-burst brake (3 fills/60s → 10-min cooloff) may throttle a healthy
  fill stream** on match day. Kept deliberately as the remaining brake.

**Worth understanding** — *Queue position*: Kalshi fills orders at the same
price first-come-first-served, like rows in an R data frame processed in
insertion order. Cancel/replace deletes your row and appends a new one at the
bottom. That's why the design holds orders through small fair moves: an early
join can be worth more than a slightly better price.

---

## Design body

### Pricing rule (`maker/engine.py::_desired`, per side)

Today: `bid = min(floor(fair_c) − margin, opp_ask − 1)`; skip if the crowd is
already tighter than fair − MAX_MARGIN_CENTS (`crowd_tighter`).

New, evaluated after the legacy price:

1. **Crowd touch** = best same-side bid with our own resting order's quantity
   subtracted at its level. If removing our qty empties the level, the touch is
   the next level down (self-exclusion — we never treat our own order as the
   crowd, so we can't chase ourselves when the real crowd retreats).
2. **Join condition**: `net_edge = fair_c − touch − maker_fee_cents(touch)`
   (NO side mirrored) with thresholds
   `TOUCH_JOIN_MIN_EDGE_CENTS = 1.0` (main rungs) /
   `TOUCH_JOIN_ALT_MIN_EDGE_CENTS = 1.5` (alt rungs), **and** the anchor is
   fresh per the existing kickoff-aware staleness gate. If met, the desired
   price is the touch (which is strictly below the opposing ask, so the
   never-cross invariant holds by construction).
3. **Fallback**: join condition not met → legacy fair−margin price, unchanged
   semantics (including `crowd_tighter` / `no_crowd` skips).
4. **Hysteresis**: if we are resting at a price that was a touch-join and
   `net_edge at our price ≥ TOUCH_JOIN_EXIT_EDGE_CENTS = 0.25`, hold the order
   even if fair moved — *unless* the crowd has moved to a new level where the
   entry condition holds (then cancel/replace to the new touch) or our price
   would violate never-cross. Below the exit edge → revert to legacy price.

One-sided by construction: the join condition holding on both sides
simultaneously would require the book to be ~2c+ wide around fair, which the
capture shows essentially never happens pre-kickoff.

### What deliberately does not change

Sizing (`MAX_QUOTE_PCT` / match / global %-caps, `ALT_SIZE_MULT`), the
goal-grid worst-case ledger (resting quotes counted as filled), every existing
pull trigger, reconciliation/startup sync, the v2 order gateway, shadow mode.

### Config (`config.py`)

`TOUCH_JOIN_MIN_EDGE_CENTS` (1.0), `TOUCH_JOIN_ALT_MIN_EDGE_CENTS` (1.5),
`TOUCH_JOIN_EXIT_EDGE_CENTS` (0.25). Dials only — no on/off switch.

### Observability

Touch-join placements write `maker_quotes` with reason `touch_join` (legacy
placements keep `quote`), preserving per-fill attribution for the post-final
review. Research firehose unchanged.

### Expected behavior (from the 3-day capture replay)

Joins concentrate on the **NO side of ESPARG-2/-3** (persistent ~1–1.5c
public-overs bias vs the anchor; net edge ≥1c on 16% of NO tick-rungs vs 6%
YES). Match-day queues historically turn over ~1–1.5× with taker-YES flow —
the flow that fills NO bids. Front-of-queue "improve" is not attempted:
spreads are 1c ~97% of the time.

### Testing

Unit tests: join condition (both sides, fee-inclusive), alt threshold,
self-exclusion (our order at/inside the touch), hysteresis hold / exit /
crowd-moved requote, never-cross under join, fallback parity with current
behavior when the join condition never fires. Full existing suite (118 tests)
must stay green.

### Version control & worktree

All work on the existing `worktree-wc-maker-depth-capture` branch in its
existing worktree (the maker exists only there; the live process runs from
this cwd — no branch switch, no new worktree). Spec committed first, then
implementation as its own commit. No merge to main tonight; post-WC the branch
merges with the rest of the maker (existing pending-merge state).

### Rollout

1. Tests green → flag in chat → clean SIGTERM of live runner (pid 55049).
2. Relaunch inline from this worktree cwd with the same env as the current
   live launch (`MAKER_MODE=live MAKER_LIVE_ACK=1`, caps as user configured).
3. Verify: heartbeat resumes, `maker_quotes` shows `touch_join` rests on the
   replay-predicted sides (NO on ESPARG-2/-3), prices equal the live touch.
4. Rollback: `git revert` the implementation commit + restart.

### Documentation

`unabated_edge/README.md` maker section gains a touch-join subsection (pricing
rule, dials, hysteresis) in the implementation commit. Memory update after
ship. Deferred work tracked in GitHub issues #1 (mid-jump breaker) and
#2 (markout brake).
