# Unabated Ticket: every +EV leg in the teaser pool (2026-10-16)

Approved by Cal 2026-10-16. Branch `feature/teaser-pool-all-ev-legs`.

## Problem

`unabated_ticket/extension/teaser.js` keeps only the 10 highest legs
(`POOL_SIZE = 10`). The cap existed only for speed.

## Approach

1. **Pool rule.** Every game's best leg whose win chance is strictly above
   `BREAK_EVEN_WIN` (70.7%, the +300 break-even for one leg) enters the pool.
   A leg at or below 70.7% stays out, even if it could make a +EV ticket with
   strong partners (Cal's decision: the individual leg must be +EV). No count
   cap; `POOL_SIZE` goes. Pairing is unchanged: every 4-leg combination of the
   pool is a candidate and the greedy adds the $200 ticket that most raises
   E[ln(1 + pnl / K)], exactly as before.
2. **Outcome budget.** `MAX_OUTCOMES` goes from 2^17 to 2^20, so held
   straights stop shedding +EV legs (on the 9/27 board, 15 straights made
   69,984 outcomes at 10 legs, so an 11th leg already overflowed).
   `structureWithinBudget`'s shedding stays as the fallback; a leg shed for
   the budget says so in the Legs list.
3. **Bug fix.** `legStatesOf` used a `Uint16Array`: a 17th pool leg's bit was
   lost, so its tickets never won. It moves to `Uint32Array`. A pool leg is
   its own game (at least 2 rows), so within 2^20 outcomes the pool has at
   most 20 legs and every `1 << leg` stays inside 32 bits.
4. **Scoring pass.** `bestTicketAt` scanned every unused ticket against all
   2^N leg states, about C(N,4) x 2^N per greedy step (~40 s at 17 legs,
   hours at 24, and in the server runner that blocks the live Kalshi taker on
   the same event loop). It becomes one sum-over-supersets pass over
   `gainByState` per step (N x 2^N); a ticket's growth is then
   loseTotal + superSum[mask]. Same tickets as the old scan, including the
   near-tie rule (`GROWTH_TIE_TOLERANCE`, the earlier ticket wins).

## Out of scope

The $200 ticket limit, the Kelly multiplier, panel and phone layout.

## Done when

- A 20-leg test board pools every leg above 70.7% and none below.
- On every existing test board the new scoring picks the same tickets as the
  old scan (a one-off run over all 57 builds before the pool rule changed;
  the committed test compares the two pick by pick on five boards).
- README and CLAUDE.md updated.

## Version control

- Branch `feature/teaser-pool-all-ev-legs` in worktree
  `.claude/worktrees/modest-newton-931339`.
- Files: `unabated_ticket/extension/teaser.js`,
  `unabated_ticket/extension/teaserview.js` (the too-few-legs wording),
  `unabated_ticket/extension/manifest.json` (0.30.0),
  `unabated_ticket/tests/teaser.test.js`,
  `unabated_ticket/tests/teaserview.test.js`, `unabated_ticket/README.md`,
  `docs/2026-09-27-unabated-ticket-teasers-plan.md`, root `CLAUDE.md`, this
  plan.
- Code, tests and docs in one commit; the pre-merge review's fixes in a second.

## Worktree

Built and tested in the worktree; merge to `main` only on Cal's yes, after
the pre-merge review; then remove the worktree and delete the branch.

## Documentation

- `unabated_ticket/README.md` § Teasers: pool rule, budget, scoring cost.
- `docs/2026-09-27-unabated-ticket-teasers-plan.md` § 3: pool rule.
- Root `CLAUDE.md`: the Teasers tab line ("top 10", "2^10", "2^17").
