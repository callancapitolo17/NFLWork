# Unabated Ticket — conditional Kelly against held bets (plan, issue #130)

Spec: GitHub issue #130 (agreed line by line). #129 (spread with total) is out
of scope and stays a grey `game` line. The bets service is untouched.

## 1. Goal

Today every line is sized alone and then adjusted by subtracting dollars
(`add $X`, `still $X against`). Dollars at different prices and numbers are
not comparable, and a moneyline held against a spread is ignored. Replace the
subtraction with conditional Kelly: held bets are fixed, the new stake `x`
maximizes `sum(prob * ln(1 + pnl / K))` over outcome rows, `K = bankroll *
multiplier`. With no related bets the answer is `kellyStakeFromEdge`, to the
cent.

## 2. Math — new `extension/condkelly.js`

Pure, IIFE + `module.exports` like `kelly.js`. No feed or bet-record
knowledge: it sees cuts, directions, dollars and probabilities.

```
solveStake({ kellyBankroll, candidate, held, ladders }) -> { stake, reason }
  candidate  { group, cut, direction: "above"|"below", netOdds, prob }
  held       [{ group, cut, direction, stake, toWin }]
  ladders    { [group]: [[cut, probAbove], ...] }   half-point cuts only
```

- **Group** = one period of one market axis (`FG`, `1H`, ...). Every bet passed
  in is already on the candidate's axis (total or margin).
- **Cuts per bet.** Half-point number: one cut. Whole number `k`: win side of
  `k+0.5`, lose side of `k-0.5`, push between (the whole-number rung's own
  fair is conditional on no push and is never read).
- **Candidate overrides its own cut** with `prob`, so a held bet on the same
  cut needs no ladder. Whole-number candidate: push mass from the ladder's
  `k±0.5`, win = `prob * (1 - push)`.
- **Rows** = gaps between a group's sorted cuts; `prob` = difference of
  `probAbove`. A negative slice down to -0.005 clamps to 0 (rows renormalised);
  beyond that `reason: "ladder not monotone"`.
- **Cross-period** (more than one group): inside the objective, sort each
  group's rows by P&L (the candidate's group depends on `x`) and pair them by
  cumulative probability (comonotonic = worst case the fairs allow).
- **Solver.** Upper bound: the largest `x` keeping `K + pnl > 0` on every row.
  Return 0 when `f(tiny x) <= f(0)`. Otherwise golden-section on `[0, hi]`
  (the objective is a min of concave functions, so one hill). Held bets that
  can already lose `K` or more: `reason: "held risk exceeds the Kelly bankroll"`.
- A candidate with no edge of its own is never solved: the stake is $0, as
  today (hedge sizing is deferred, section 11).
- A missing ladder cut throws (expected vs found) — the caller filters first.

## 3. Ladder — new `extension/ladder.js`

Pure. Input is `scannerState.lines`.

- `groupLinesByEvent(lines)` — one pass, `Map(eventId -> lines)`. Panel caches
  it per feed change (same pattern as `boardLinesCache`).
- `buildLadder(eventLines, { periodTypeId, axis, sport })` -> `Map(cut ->
  probAbove)`. Half-point rungs only. Per rung: median over every line with a
  `bacr` of `americanToProb(bacr)`, oriented to P(above cut).
  - Total: cut = points; side 0 (Over) = above, side 1 (Under) = 1 - p.
  - Margin (M = away - home): away at `a` -> cut `-a`, above; home at `h` ->
    cut `h`, 1 - p. Moneyline lines feed cuts `+0.5` (away) and `-0.5` (home),
    except soccer (three-way).
- `probAt(ladder, cut)` -> `{ prob }` or `{ reason }`: no rung, or the rung's
  fair equals a neighbouring rung's (Unabated flat-lines tails). The
  `(-0.5, +0.5)` pair is exempt.
- `periodTypeIdOf(period)` — inverse of `feed.PERIODS`; `F5` / `I1` have no id
  and fall out as "no fair".

## 4. Bets — `extension/bets.js`

- `tierOf`: moneyline vs spread in the same period tiers by team (`same_side`
  / `opposite`, never `same_line`). Same market in another period becomes
  `related_same` / `related_opposite` (names from the stale
  `feature/unabated-related-period-bets`, 87e77d6; that branch is not merged).
  Different market stays `same_game`.
- `TIER_RANK` and `TIER_LABELS` gain the two tiers; tags read `same side` /
  `other side` (the bet text already starts with its period).
- `labelOf`: the `now <number>` suffix only when bet type and period both
  match the row (a moneyline has no number to have moved from).
- New `positionOf(bet, line)` on every match: `{ axis, period, cut, direction,
  stake, toWin }` or `{ reason }`. Record-level guards live here: no stake /
  payout, parlay leg, bet type other, Kalshi NO moneyline
  (`kalshi_no_side_includes_tie`), soccer moneyline, quarter line, side not
  resolvable. Team side comes from `sideIndexByTeamKeys ?? ByVenueId ??
  ByRotation` (Unabated's frame, not the venue's). Same helper shapes the row
  itself (`linePosition`).
- `exposureOf` and the `exposure` field of `annotateRows` are deleted: dollars
  now come from the bets that are in the math (section 5).

## 5. Advice and words — `extension/betsview.js`

- `stakeAdvice({ line, price, edgePct, bankroll, multiplier, matches,
  ladderOf })` replaces `stakeAdvice(stake, exposure)`. Steps: standalone
  stake -> split matches into in-math (same axis, position known, fair known)
  and not-sized (reason) -> `condkelly.solveStake` -> advice:

```
{ kind: "none" | "sized" | "declined",
  bet, alone, verb: "add" | "bet", held, against, reason,
  matches: [...match, inMath, note] }
```

  - `none`: nothing in the math, `bet = alone` (today's number exactly).
  - `declined`: whole calc refused (`ladder not monotone`, ...): `bet = alone`
    and every would-be bet is greyed with the reason.
- `suggestedBetAmount(advice)` -> `advice.bet` (0 when unsizable). Filter,
  alerts, Ticket and Copy all read it.
- `stakeAdviceWords` -> `{ verb, bet, alone }`; `stakeAdviceLine` -> `add
  $188.32, $183.00 alone`. `have / target / net / note` are removed.
- `badges(flag)` replaces `badgeText` / `badgeKind`: `held $X` and
  `against $X`, both when both exist; bare kind from the tier when no dollars;
  `game` otherwise.
- `relatedLines`: in-math -> coloured tag; not in math -> grey tag + grey text
  (`game · not sized`, or `no fair at 36.5`). In-math first, cap of three
  stays in the panel.

## 6. Panel — `panel.js`, `panel.html`, `panel.css`, `manifest.json`

- `panel.html`: load `ladder.js` and `condkelly.js` before `panel.js`.
- Ladder cache: `linesByEventCache` + per-(event, period, axis) ladders, reset
  in the scanner `onChange` next to `boardLinesCache`.
- `withBetFlags`: build advice per row; `exposureDollars` reads
  `advice.held + advice.against`.
- `fillStakeCell`: `add|bet $stake` + one small `$X alone` line (only when
  sized and different). `at-size` class when the bet is $0.
- `relatedBlock` / `renderBetBanner`: grey class for not-sized lines.
- `renderStakeExposure` + `actedStake`: the same three pieces as the rail.
- Alerts: skip rows whose conditional stake is $0 (replaces the `at_size`
  check); `meetsMinStake` unchanged in shape.
- CSS: `.related-line.not-sized`, `.bet-match.not-sized`, second badge
  spacing. Manifest 0.8.0 -> 0.9.0.

## 7. Tests (Kelly bankroll 8000)

`tests/condkelly.test.js`, `tests/ladder.test.js` new; `bets.test.js`,
`betsview.test.js` updated. All checked in a scratch prototype already:

| Case | Expect |
|---|---|
| no held bets | equals `kellyStakeFromEdge` to the cent (297.60) |
| same line, held $181.50 at +125 | 116.10 |
| same line, held $400 | 0 |
| UConn: hold U47.5 $181.50 +203, new U52.5 +125 | 121.28 |
| UConn: hold U52.5 $181.50 +125, new U47.5 +203 | 24.10 |
| New Mexico @ Oklahoma | 295.99 |
| Lions @ Bills | 188.32; unchanged with Bills -8.5 $472 held (betsview level) |
| Chargers ML vs held +3.5 $400 | 0 |
| cross-period 1H Over $300 -> FG Over | 408.39 |

Guards, one test each: no rung, flat rung, `±0.5` exemption, missing stake,
parlay leg, Kalshi NO moneyline, soccer moneyline, quarter line, unknown
period, non-monotone beyond 0.005 (declined) and within it (clamped),
whole-number push row, held risk above K.

## 8. Version control

- Branch `feature/unabated-ticket-conditional-kelly` (this worktree's branch,
  renamed before any file was written).
- Commits: (1) this plan; (2) `condkelly.js` + `ladder.js` + tests;
  (3) `bets.js` tiers/positions + tests; (4) `betsview.js` + `panel.*` +
  manifest + tests; (5) README + `CLAUDE.md`.
- New: `extension/condkelly.js`, `extension/ladder.js`,
  `tests/condkelly.test.js`, `tests/ladder.test.js`, this plan.
  Modified: `bets.js`, `betsview.js`, `panel.js`, `panel.html`, `panel.css`,
  `manifest.json`, `tests/bets.test.js`, `tests/betsview.test.js`,
  `eslint.config.js` (new globals, if it lists them), `unabated_ticket/README.md`,
  root `CLAUDE.md`.
- No merge, no push without an explicit yes. `./unabated_ticket/check.sh` on
  the branch, then the pre-merge review, then ask.

## 9. Worktree

Work happens in `.claude/worktrees/suspicious-sanderson-e1c108`. No DuckDB
files are touched or copied. After an approved merge: `git worktree remove`
the path and `git branch -d feature/unabated-ticket-conditional-kelly`. Also
after merge, with a yes: delete the stale `feature/unabated-related-period-bets`
(the issue retires it).

## 10. Documentation

Same branch, after the code is final: `unabated_ticket/README.md` — Stake
(conditional Kelly, rail wording), Bets (tiers, in-math vs grey,
badges, guards), Tests (counts + the two new suites), design-decisions log
entry reversing 2026-09-12 "show it, don't size off it"; root `CLAUDE.md`
Unabated Ticket blurb; manifest description if the wording changes.

## 11. Known limits

- **Hedge sizing is deferred** (user, 2026-09-19): no `hedge only` label and
  no stake on a price with no edge of its own, even when a held bet on the
  other side would make the math want one. Such a line reads $0, as today.

- `bacr` is a whole American price: stakes are good to about ±$10.
- A cross-period hedge gets no credit (accepted in the issue).
- Re-ranking a card's best line by conditional stake is out of scope: cards
  still rank by the standalone stake.
