# Unabated Ticket — Buckeye teasers (plan)

A Teasers tab that builds the set of Buckeye 6-point teasers to bet right now,
each ticket at or under Buckeye's $200 limit, sized with the panel's Kelly
settings. Unabated's own teaser calculator prices one hand-built ticket at a
typed-in price; this prices every Buckeye leg on the board and spreads the
Kelly amount over as many $200 tickets as it is worth.

User decisions (2026-09-27): Buckeye only; legs priced off Unabated's fair
(not the exchange median, although at teaser numbers Unabated sits up to 3
points above Kalshi/Novig/Polymarket/ProphetX); $200 maximum per ticket.

## 1. Goal

Open the tab, see "bet these N tickets, $200 each, in this order", place them
at Buckeye, tick each one placed. The next visit sizes around the tickets
already placed instead of suggesting them again.

## 2. Legs — new `extension/teaser.js` (pure)

- Source: Buckeye (source 59) main full-game spreads and totals, NFL and CFB,
  on the board (`statusId` 1), game not started, line age within the Edges
  "Max line age h" setting (same meaning: a dead feed keeps old numbers).
- Teased number (`TEASER_POINTS = 6`): spread points + 6; Over - 6; Under + 6.
- Cut and direction, the `condkelly.js` convention: away spread at `a` wins
  above `-a`, home at `h` below `h`, Over above, Under below.
- Fair from `ladder.buildLadder` over the event's lines of THAT bet type only.
  Spreads leave out the moneyline: Unabated's moneyline fair and its alt
  spreads disagree by about 2.5 points near zero (SEA @ WAS 77.1% vs 79.5%),
  which put -0.1% on "Seattle by exactly 1".
- Half-point number: win = P(above) or 1 - P(above), no push. Whole number
  `k`: push = P(above k-0.5) - P(above k+0.5), win from the far rung.
- A leg with no rung, a flat rung (`ladder.probAbove` reasons) or a push below
  -0.005 is listed grey with its reason and never used.

## 3. Tickets and grading

- Pool: the best leg per game (one leg per game, so legs are independent),
  top `POOL_SIZE = 8` by win chance. 3^8 = 6,561 joint outcomes.
- Candidate tickets: every combination of the leg counts Buckeye pays
  (`BUCKEYE_PAYOUTS`; 4 teams +300 known, others pending).
- Grading per outcome: any leg lost -> -1; ties by `BUCKEYE_TIE_RULE`
  (pending; e.g. "reduce": a 4-teamer with one push pays the 3-team price).
- Ticket EV and "fair" (the full-ticket payout that makes EV zero with the
  other outcomes graded as usual; closed form, EV is linear in it).

## 4. Portfolio — Kelly under the $200 cap

- Objective: E[ln(1 + pnl / K)], `K = bankroll * multiplier` (the panel's
  settings; same objective `condkelly.js` uses for held bets).
- Greedy in $200 steps (`TICKET_MAX_STAKE = 200`): add the ticket that most
  raises the objective, repeat until none does; the last ticket may be
  partial (golden section on 0-200). No ticket twice.
- Placed tickets are fixed P&L in the objective. A game with a placed leg at
  another number (Buckeye's line moved) gets rows between ALL its cuts, as in
  `condkelly.js`, so the two legs are the same game, not independent. A
  placed total and a new spread on one game cannot be combined: that game
  leaves the pool.
- Today's board (prototype, $30k x 0.25, 4-teamers only): 12 tickets, $2,275,
  expected +$711 (31%), wins money 52%, all lose 14%; each leg in 4-8
  tickets; solved in 53 ms.

## 5. Placed tickets — `chrome.storage.local.teaserPlaced`

`[{ id, placedAt, stake, legs: [{ eventId, leagueId, betTypeId, sideIndex,
teasedPoints, label, eventStart }] }]`. Written by the Placed button, removed
by Undo, dropped once every leg's game has started. Nothing leaves the
browser (Buckeye bets are not in the bets service).

## 6. Scanner — football always loaded

`scanner.start` gets the Edges leagues plus NFL and CFB; `feed.selectEdges`
takes a `leagueIds` filter so the Edges tab and its alerts still show only the
sports you ticked.

## 7. Panel — `panel.html`, `panel.js`, `panel.css`, `manifest.json`

- Tab "Teasers" after Bets, count = tickets to bet.
- Summary line: tickets, dollars, expected profit, chance to make money,
  chance every ticket loses.
- Ticket cards in greedy order: legs, stake, EV, Copy, Placed. Placed
  tickets collapse into a "Placed" section with Undo.
- Legs list: teased number, from Buckeye's line and price, win and push
  chance, "in N tickets"; grey rows with the reason for unpriced legs.
- Recomputed on every scanner update while the tab is visible.
- Manifest 0.12.1 -> 0.13.0.

## 8. Tests — `tests/teaser.test.js` (+ `feed.test.js`)

Synthetic ladders, no live snapshots: away/home/Over/Under teased numbers;
whole-number push; the moneyline left out of the spread ladder; no rung,
flat, non-monotone; grading (no push, one push reduces, all push refunds);
probabilities sum to 1; fair closed form; one ticket under the cap equals
the single-bet Kelly stake; the greedy never exceeds $200 per ticket nor
repeats a ticket; a placed ticket lowers the next suggestion; a moved line
on a placed game shares rows; `selectEdges` league filter.

## 9. Version control

- Branch `feature/unabated-ticket-teasers`, worktree
  `.claude/worktrees/buckeye-teasers` (from local `main`, which is 10
  commits ahead of `origin/main`).
- Commits: (1) this plan; (2) `teaser.js` + tests; (3) scanner leagues +
  `selectEdges` filter + test; (4) panel tab, CSS, manifest; (5) README +
  `CLAUDE.md`.
- New: `extension/teaser.js`, `tests/teaser.test.js`, this plan. Modified:
  `feed.js`, `panel.js`, `panel.html`, `panel.css`, `manifest.json`,
  `tests/feed.test.js`, `eslint.config.js` (globals, if listed),
  `unabated_ticket/README.md`, root `CLAUDE.md`.
- `./unabated_ticket/check.sh` on the branch, pre-merge review, then ask. No
  merge, no push without an explicit yes.

## 10. Worktree

Work happens in `.claude/worktrees/buckeye-teasers`. No DuckDB files are
touched. After an approved merge into local `main`: `git worktree remove` the
path and `git branch -d feature/unabated-ticket-teasers`.

## 11. Documentation

Same branch, after the code is final: `unabated_ticket/README.md` — Teasers
section (legs, grading, portfolio, placed tickets, limits), Tests counts,
decisions log (Unabated fair over exchanges; $200 cap; one leg per game);
root `CLAUDE.md` Unabated Ticket blurb; manifest description.

## 12. Known limits

- Unabated's fair at teaser numbers runs up to 3 points above the exchanges;
  EVs and stakes are as good as that fair.
- One leg per game: same-game spread + total correlation is not modelled.
- Football, full game, 6 points only. Basketball teasers wait for Buckeye's
  basketball table.
- Straight bets held elsewhere on the same games are not in the math.
- Placed tickets live in this browser only.
