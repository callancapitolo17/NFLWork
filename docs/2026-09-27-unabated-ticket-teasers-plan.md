# Unabated Ticket — Buckeye teasers (plan)

A Teasers tab that builds the set of Buckeye 6-point teasers to bet right now,
each ticket at or under Buckeye's $200 limit, sized with the panel's Kelly
settings. Unabated's own teaser calculator prices one hand-built ticket at a
typed-in price; this prices every Buckeye leg on the board and spreads the
Kelly amount over as many $200 tickets as it is worth.

User decisions: Buckeye only; legs priced off Unabated's fair (not the
exchange median, although at teaser numbers Unabated sits up to 3 points above
Kalshi/Novig/Polymarket/ProphetX); $200 maximum per ticket; 4-team tickets
only (+300 is Buckeye's best teaser price); a push counts as a loss
(conservative, Buckeye's tie rule unconfirmed); the panel's bankroll and Kelly
multiplier; a juiced Buckeye line teases like a -110 one (2026-09-27/28).

## 1. Goal

Open the tab, see "bet these N tickets, $200 each, in this order", place them
at Buckeye, tick each one placed. The next visit sizes around the tickets
already placed instead of suggesting them again.

## 2. Legs — new `extension/teaser.js` (pure)

- Source: Buckeye (source 59) main full-game spreads and totals, NFL and CFB,
  on the board (`statusId` 1), game not started, line age within the Edges
  "Max line age h" setting (same meaning: a dead feed keeps old numbers).
  The Buckeye price is ignored: -120 and +100 tease the same.
- Teased number (`TEASER_POINTS = 6`): spread points + 6; Over - 6; Under + 6.
- Push = loss, so a leg wins exactly when the half-point one step against
  you wins: Bills -1 wins when -1.5 does, Browns +8 when +7.5 does, Over 38
  when Over 38.5 does. The leg's fair is Unabated's price at that half-point,
  read off `ladder.buildLadder`: Bills -1 -> -1.5 at -285 -> 74.0%. One
  lookup per leg, no push math.
- Cut and direction, the `condkelly.js` convention: away spread at `a` wins
  above `-a`, home at `h` below `h`, Over above, Under below.
- The ladder is built from the event's lines of THAT bet type only. Spreads
  leave out the moneyline: Unabated's moneyline fair and its alt spreads
  disagree by about 2.5 points near zero (SEA @ WAS 77.1% vs 79.5%).
- A leg with no rung or a flat rung (`ladder.probAbove` reasons) is listed
  grey with its reason and never used.

## 3. Tickets and grading

- Pool: the best leg per game (one leg per game, so legs are independent),
  top `POOL_SIZE = 10` by win chance. Every leg is win-or-lose, so 2^10 =
  1,024 joint outcomes (26 ms on the 9/27 board; 12 legs took 221 ms and
  added nothing).
- Candidate tickets: every 4-leg combination of the pool (210).
- Grading: all four legs win -> +3 per dollar (+300); anything else -> -1.
- Ticket EV = P(all four win) x 4 - 1.

## 4. Portfolio — Kelly under the $200 cap

- Objective: E[ln(1 + pnl / K)], `K = bankroll * multiplier` (the panel's
  settings; same objective `condkelly.js` uses for held bets).
- Greedy in $200 steps (`TICKET_MAX_STAKE = 200`): add the ticket that most
  raises the objective, repeat until none does; the last ticket may be
  partial (golden section on 0-200, whole dollars). No ticket twice.
- Placed tickets are fixed P&L in the objective until their LAST game starts.
  New tickets never use a started game, but a placed ticket whose 10:00 legs
  are live still rides on its 1:05 leg, so each placed leg keeps the win
  chance saved when it was placed and that is used once its game starts.
- A game with a placed leg at another number (Buckeye moved the line) gets
  rows between all its cuts, as in `condkelly.js`, so the two legs are the
  same game, not independent. A game with a placed leg only offers new legs
  on the same market (spread or total).
- The list holds still: it is rebuilt only when a Buckeye number moves, a
  game joins or leaves, a ticket is placed or undone, a setting changes, or
  a leg in the list moves `REBUILD_FAIR_MOVE = 1` point or more since the
  last build. Smaller ticks update the EVs in place. Why: many 4-leg combos
  are near-tied, so rebuilding on every tick gave lists sharing 1 of 7
  tickets at 9:06 and 9:37 on 9/27 (no leg moved over 0.51 points) that are
  worth the same at the same fairs (+$326 vs +$324 expected). A frozen list
  alone would miss news: the 49ers leg down 5 points halves the 9:37 list's
  value after risk ($141 -> $71).
- Splitting finer does not cut the swing (measured 9/27): 92-97% of the
  variance is how many legs cover. The Kelly multiplier is the variance dial.
- 9/27 board at $20k x 0.25, push = loss: 7 tickets, $1,292, expected +$324.

## 5. Placed tickets — `chrome.storage.local.teaserPlaced`

`[{ id, placedAt, stake, legs: [{ key, eventId, leagueId, betTypeId,
sideIndex, basePoints, teasedPoints, label, eventStartMs, win }] }]`. `win`
is the leg's fair when placed. Written by the Placed button, removed by
Undo, dropped once every leg's game has started. Nothing leaves the browser
(Buckeye bets are not in the bets service).

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
- Legs list: teased number, from Buckeye's line, win chance, line age,
  "in N tickets"; grey rows with the reason for unpriced legs.
- Empty states: feed loading; fewer than 4 priced games on Buckeye's board;
  no ticket worth betting.
- Manifest 0.12.1 -> 0.13.0.

## 8. Tests — `tests/teaser.test.js` (+ `feed.test.js`)

Synthetic ladders, no live snapshots: away/home/Over/Under teased numbers;
a whole number reads the half-point one step against you; the moneyline left
out of the spread ladder; no rung, flat; grading (all four win, anything
else loses); one ticket under the cap equals the single-bet Kelly stake; the
greedy never exceeds $200 per ticket nor repeats a ticket; a placed ticket
lowers the next suggestion; a placed ticket with a started game still counts
with its saved win chance; a moved line on a placed game shares rows; the
list is not rebuilt when only fairs move; `selectEdges` league filter.

## 9. Version control

- Branch `feature/unabated-ticket-teasers`, worktree
  `.claude/worktrees/buckeye-teasers` (from local `main`, which is 10
  commits ahead of `origin/main`).
- Commits: (1) this plan; (2) `teaser.js` + tests; (3) scanner leagues +
  `selectEdges` filter + test; (4) panel tab, CSS, manifest; (5) README +
  `CLAUDE.md`.
- New: `extension/teaser.js`, `tests/teaser.test.js`, this plan. Modified:
  `feed.js`, `panel.js`, `panel.html`, `panel.css`, `manifest.json`,
  `tests/feed.test.js`, `unabated_ticket/README.md`, root `CLAUDE.md`.
- `./unabated_ticket/check.sh` on the branch, pre-merge review, then ask. No
  merge, no push without an explicit yes.

## 10. Worktree

Work happens in `.claude/worktrees/buckeye-teasers`. No DuckDB files are
touched. After an approved merge into local `main`: `git worktree remove` the
path and `git branch -d feature/unabated-ticket-teasers`.

## 11. Documentation

Same branch, after the code is final: `unabated_ticket/README.md` — Teasers
section (legs, grading, portfolio, placed tickets, limits), Tests counts,
decisions log (Unabated fair over exchanges; $200 cap; 4-team only; push =
loss; one leg per game); root `CLAUDE.md` Unabated Ticket blurb; manifest
description.

## 12. Known limits

- Unabated's fair at teaser numbers runs up to 3 points above the exchanges;
  EVs and stakes are as good as that fair.
- Push = loss understates a ticket if Buckeye really drops a push to the
  3-team price (9/27: +$329 vs +$418 expected on the same board).
- One leg per game: same-game spread + total correlation is not modelled.
- Football, full game, 6 points only.
- Straight bets held elsewhere on the same games are not in the math.
- Placed tickets live in this browser only.
