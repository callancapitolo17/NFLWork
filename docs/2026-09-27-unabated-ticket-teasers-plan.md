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
multiplier; a juiced Buckeye line teases like a -110 one; BFA is Buckeye, so
placed tickets are read from the bets service's BFA open bets, nothing to
mark (2026-09-27/28).

## 1. Goal

Open the tab, see "bet these N tickets, $200 each, in this order", place them
at BFA. The bets service's BFA pull brings them back as open tickets — nothing
to click — and the next build sizes around them instead of suggesting them
again.

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
- Placed tickets = every open BFA teaser, from the tool or not (section 5):
  fixed P&L in the objective, +toWin when every leg wins and -stake
  otherwise (push = loss), until the ticket's LAST game starts. New tickets
  never use a started game, but a placed ticket whose 10:00 legs are live
  still rides on its 1:05 leg, so a started leg — and a leg no board game
  matches — counts as still alive. Conservative: the ticket keeps its full
  weight on the legs still to play, and no fair has to be saved.
- A game with a placed leg at another number (Buckeye moved the line) gets
  rows between all its cuts, as in `condkelly.js`, so the two legs are the
  same game, not independent. A game with a placed leg only offers new legs
  on the same market (spread or total).
- The list holds still: it is rebuilt only when a Buckeye number moves, a
  game joins or leaves, an open BFA teaser appears or closes, a setting
  changes, or
  a leg in the list moves `REBUILD_FAIR_MOVE = 1` point or more since the
  last build. Smaller ticks update the EVs in place. Why: many 4-leg combos
  are near-tied, so rebuilding on every tick gave lists sharing 1 of 7
  tickets at 9:06 and 9:37 on 9/27 (no leg moved over 0.51 points) that are
  worth the same at the same fairs (+$326 vs +$324 expected). A frozen list
  alone would miss news: the 49ers leg down 5 points halves the 9:37 list's
  value after risk ($141 -> $71).
- Splitting finer does not cut the swing (measured 9/27): 92-97% of the
  variance is how many legs cover. The Kelly multiplier is the variance dial.
- No per-leg cap (user declined 2026-09-28, shown "no leg on more than half
  the tickets": -6% value after risk). The same score stops the set and
  limits each leg: every extra ticket on a leg is scored against the
  tickets already riding on it. Each leg lands near its "target", quarter
  Kelly as a straight bet at -241 (the per-leg price inside +300): 9/27
  49ers $892 vs $1,089, Lions $400 vs $354. Filling targets alone (no log
  score) stacked the top pair on 4 of 7 tickets and was worth $135 vs $141.
- 9/27 board at $20k x 0.25, push = loss: 7 tickets, $1,292, expected +$324.

## 5. Placed tickets — BFA open bets from the bets service

- The panel already polls `/bets.json` every 30 s. An open BFA teaser arrives
  as one record per leg (`isParlayLeg`, `parlayId`, `legCount`, the ticket's
  `stake` / `toWin` on every leg, `raw.headerDescription` "4 TEAM TEASERS",
  the leg's teased `points`, `rotation` and `eventStart`). Sunday's 8 tickets
  came through this way, legs at the teased numbers.
- Legs are grouped by `parlayId` into tickets and joined to board games with
  `bets.js`'s matcher (rotation + start; BFA rotations are Unabated's own).
  `bets.js` gains `teaserLegPositionOf(bet, line)`: `positionOf` without the
  parlay-leg and stake guards (a leg's dollars are the ticket's).
- BFA is pulled every 300 s today (open bets + history in one poll), so a
  ticket could take 5.5 min to show and be suggested again meanwhile. The
  bets service splits it: open bets every 60 s, history every 300 s
  (`sources/bfa.py`, `config.py`). The tab shows how old the last BFA pull
  is; the bets service down -> red banner, the list stays up.

## 6. Scanner — football always loaded

`scanner.start` gets the Edges leagues plus NFL and CFB; `feed.selectEdges`
takes a `leagueIds` filter so the Edges tab and its alerts still show only the
sports you ticked.

## 7. Panel — `panel.html`, `panel.js`, `panel.css`, `manifest.json`

- Tab "Teasers" after Bets, count = tickets to bet.
- Summary line: tickets, dollars, expected profit, chance to make money,
  chance every ticket loses.
- Ticket cards in greedy order: legs, stake, EV, Copy. Open BFA teasers
  are listed read-only under "Open at BFA" with the age of the last BFA
  pull. Change from the approved mockup: frame 2's Placed button and Undo
  are gone — BFA says what is placed.
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
greedy never exceeds $200 per ticket nor repeats a ticket; open BFA teaser
legs group into tickets and join their games by rotation; a placed ticket
lowers the next suggestion; a started or unmatched leg counts as alive; a
moved line on a placed game shares rows; the list is not rebuilt when only
fairs move; `selectEdges` league filter. pytest: BFA open bets every 60 s,
history every 300 s.

## 9. Version control

- Branch `feature/unabated-ticket-teasers`, worktree
  `.claude/worktrees/buckeye-teasers` (from local `main`, which is 10
  commits ahead of `origin/main`).
- Commits: (1) this plan; (2) `teaser.js` + `bets.js` leg position + tests;
  (3) scanner leagues + `selectEdges` filter + test; (4) BFA open bets every
  60 s + pytest; (5) panel tab, CSS, manifest; (6) README + `CLAUDE.md`.
- New: `extension/teaser.js`, `tests/teaser.test.js`, this plan. Modified:
  `feed.js`, `bets.js`, `panel.js`, `panel.html`, `panel.css`,
  `manifest.json`, `tests/feed.test.js`, `tests/bets.test.js`,
  `bets_service/sources/bfa.py`, `bets_service/config.py`,
  `bets_service/tests/test_bfa.py`, `unabated_ticket/README.md`, root
  `CLAUDE.md`. The bets service must be restarted after merge.
- `./unabated_ticket/check.sh` on the branch, pre-merge review, then ask. No
  merge, no push without an explicit yes.

## 10. Worktree

Work happens in `.claude/worktrees/buckeye-teasers`. No DuckDB files are
touched. After an approved merge into local `main`: `git worktree remove` the
path and `git branch -d feature/unabated-ticket-teasers`.

## 11. Documentation

Same branch, after the code is final: `unabated_ticket/README.md` — Teasers
section (legs, grading, portfolio explained with the per-leg targets and a
step-by-step table like 9/27's, placed tickets from BFA, limits), the
bets service's BFA cadence, Tests counts, decisions log (Unabated fair over
exchanges; $200 cap; 4-team only; push = loss; one leg per game; BFA is
Buckeye, no marking); root `CLAUDE.md` Unabated Ticket blurb; manifest
description.

## 12. Known limits

- Unabated's fair at teaser numbers runs up to 3 points above the exchanges;
  EVs and stakes are as good as that fair.
- Push = loss understates a ticket if Buckeye really drops a push to the
  3-team price (9/27: +$329 vs +$418 expected on the same board).
- One leg per game: same-game spread + total correlation is not modelled.
- Football, full game, 6 points only.
- Straight bets on the same games are in the math since section 13, but only
  on the leg's own market and full game: a total held under a spread leg, a
  1H bet, and a bet on a game no new ticket can use are left out.
- The Edges tab sizes a straight bet with the open BFA teasers on its game
  since section 14; other parlays are still left out.
- A ticket placed in the last ~90 s may not be in the tab yet (BFA every
  60 s + panel poll 30 s); the pull age is on screen.

## 13. Straight bets on the same games (2026-09-30)

Cal's last item before merge: size the teasers around the straight bets he
already holds. User decisions 2026-09-30: "bets already placed" means open
straight bets (spreads, totals, moneylines) at every venue the bets service
reads, on the same games; same method as conditional Kelly (#130) — no
correlation number, the game's fair ladder split into results; design below
approved after the 9/27 walk-through.

Measured on 9/27 (16:37 UTC board, K $5,000): 15 open straights on 9 of the
10 pool games ($4,542) cut the set from 7 tickets / $1,292 to 4 / $676, and
the set's value after risk from -$8 to +$62 on top of the straights. Same
side shrinks a leg (Lions and Patriots out, 49ers $892 -> $476, Seahawks
$800 -> $400); the other side grows it (Steelers +9 $0 -> $200 with Bengals
-3.5 held).

- **Which bets** — `teaser.heldStraights(records, boardLines, {now,
  ladderOf})`: every open record that is not a parlay leg, joined to its
  game by `bets.annotateRows` (the Edges tab's matcher: pins, venue ids,
  names, rotation), placed by `bets.positionOf` (its reasons: bet type, no
  stake, no side, quarter line, Kalshi NO + tie, three-way), full game only
  (else "1H bet"), priced at `condkelly.cutsNeeded` off the teaser ladder of
  its market (spread-only margin ladder, so a moneyline reads the ±0.5
  rungs; 28 of 28 present on 9/27). Started games drop out; no fair at a
  needed number -> left out, named.
- **Which count** — only on a game a pool leg rides on, on that leg's
  market (axis). Other market, other period, games outside the pool: left
  out, as on the Edges tab (#129). A straight never adds a game, never
  changes which leg a game offers, never restricts a game's market.
- **How they count** — the game's factor gets the straight's cuts (a whole
  number both sides of its push); its P&L per row is +toWin / 0 / -stake
  (condkelly's win/push/lose rule, `resultAt` exported for it), added to the
  fixed P&L like an open teaser's. The score, the greedy, the $200 cap and
  the one-partial rule are unchanged. A straight whose cut makes its game's
  ladder non-monotone is left out with the reason (one bad rung must not
  blank the list). Straights and open teasers that can already lose K ->
  no tickets, with that reason (today's rule for open teasers).
- **Speed** — straights multiply a game's rows (9/27: 1,024 -> 69,984
  outcomes; the unchanged greedy took 412 ms in node). The greedy groups
  outcomes by which pool legs win (a 10-bit state, 1,024 states): a
  ticket's growth is a sum over states, not outcomes, and the per-combo
  outcome lists go. `MAX_OUTCOMES` 2^14 -> 2^17. Past it the weakest pool
  leg on a game no open teaser rides on is shed, and its game's straights
  with it.
- **The list holds still** — the straights in the math and their reference
  fairs join `inputKey`: a straight placed or closed on a pool game, or its
  fair moving `REBUILD_FAIR_MOVE`, rebuilds; one on another game does not.
- **Summary** — the tickets' and open teasers' figures stay teaser-only (a
  straight's P&L is not in Expected / Makes money / All lose);
  `summary.straights` {count, stake, games} feeds the note line.
- **Legs list** — a pool row carries `held $X` (straights on the leg's
  direction) and `against $Y` (the other), the Edges tab's tags; a bare
  `game` when the game has straights but none counts. Mockup approved
  before the panel is touched.
- **Tests** (`tests/teaser.test.js`) — heldStraights positions (away, home,
  Over, Under, moneyline, whole number = two cuts) and reasons (1H, no
  fair, parlay leg, started); same side lowers a leg's dollars, other side
  raises them; other market and non-pool games change nothing; a push row
  pays 0; rebuild on a placed straight / a 1-point move, not on another
  game or a 0.5-point move; held risk >= K -> no tickets; a non-monotone
  straight is left out and the list stands; the grouped greedy matches the
  existing tests unchanged; the summary leaves straights out.
- **Version control** — same branch `feature/unabated-ticket-teasers`
  (local main merged in first: c709fe34, manifest 0.15.0 since main is
  0.14.0). Commits: (1) this section; (2) `teaser.js` grouped greedy +
  straights, `condkelly.js` export, tests; (3) panel tags + note, CSS if
  needed; (4) README + root `CLAUDE.md`. Files: `extension/teaser.js`,
  `extension/condkelly.js`, `extension/panel.js`, `tests/teaser.test.js`,
  `tests/condkelly.test.js` (export), `unabated_ticket/README.md`,
  `CLAUDE.md`, this plan.
- **Worktree** — `.claude/worktrees/buckeye-teasers` as before; no DuckDB
  touched. After an approved merge: `git worktree remove` + `git branch -d`,
  restart the bets service (launchd `kickstart`, BFA open bets every 60 s)
  and reload the unpacked extension.
- **Documentation** — README Teasers: a "Straight bets on the same games"
  paragraph with the 9/27 numbers, the Limits line, Tests counts, a
  decisions-log entry; root `CLAUDE.md` Unabated Ticket blurb.

## 14. Open teasers on the Edges tab (2026-09-30)

The one-way gap section 12 listed: the Edges tab (and the Ticket tab and
alerts) sized a straight bet without the open teasers on its game, because
conditional Kelly (#130) leaves every parlay leg out. User decisions
2026-09-30, after the 9/27 walk-through below and the mockup: count them
exactly (a ticket's other legs enumerated — not averaged, not counted as
won); open BFA teasers only; a `teasers $X` chip.

Measured on 9/27 (16:46 UTC board, right after Cal's 8 teasers, K $5,000):
all 14 Novig/Kalshi/ProphetX rows that said bet something on the same side
as a teaser leg go to $0 — Colts +7.5 −270 at Novig (+2.3%, 3 legs on Colts
+7.5) $316 → $0: the teasers already carry those sides at about Kelly size
at a better price (a leg inside +300 is −241). The other side grows, the
standing lift rule: Chargers +6.5 +120 (4 legs on Bills −1) $137 → $426,
Commanders +6.5 $381 → $743. Shortcuts on the Chargers row: other legs
averaged $492, counted as won $1,104.

- **Which tickets** — `teaser.openTeasers`' tickets (the open BFA teasers
  the Teasers tab reads), in play, with a live leg on the row's game. In the
  math when that leg sits on the row's market (spread / moneyline vs
  total), with a fair at its number on the row's ladder; a 1H row takes a
  full-game leg on its own direction (worst-case pairing, as a straight)
  and leaves the other direction out. Else left out and named: `game · not
  sized`, `other period · not sized`, `no fair at …`. Other parlays stay
  `parlay leg`. A leg on another game counts as won once its game has
  started (the board's clock, else BFA's): BFA closes a teaser within
  minutes of a losing leg's game ending (9/27: the five Seahawks tickets at
  20:30 UTC, ~10 min after that game, with the 49ers leg still to play), so
  an open ticket's finished legs won; only a leg in progress is assumed. A
  leg still to play with no board game (CFB still loading, a failed league)
  or no fair is unknown, not won: the ticket is left out, `other leg not
  priced: …` (pre-merge review: counted as won, it sized $1,301 where the
  loaded board says $361).
- **How they count** — `condkelly.solveStake` gains `tickets` [{group, cut,
  direction, stake, toWin, others: [{factor, cut, direction}]}] and
  `factors` {factor: [[cut, P(above)], ...]}. The ticket's leg cuts the
  game's rows like a held bet (a half-point: a teaser push loses); its other
  legs sit on other games, independent of this one: their rows are
  enumerated (each game cut at every number an other leg needs) and grouped
  by which tickets they leave alive (legs shared across tickets make those
  joint), and the score is the probability-weighted sum of today's growth
  (per period worst-case pairing unchanged) over those groups. No ticket:
  one group, today's number to the cent. Budget: 2^16 joint other-game
  states, past it the calc is declined with the reason (9/27: at most 128
  states, 20 groups per row). A decline the teasers alone cause (budget, an
  other game's ladder out of line, their risk past K) re-solves on the
  straights and names the reason on the teasers, so adding teasers never
  leaves the row worse sized than ignoring them.
- **On the row** (mockup approved 2026-09-30) — chip `teasers $X`: the
  stakes of the tickets in the math, green (`held` style) on the row's
  direction, red (`against`) on the other, its tooltip each ticket's legs.
  Related bets: one line per leg number — "3 teasers on Colts +7.5 · $600 ·
  BFA · other legs win 40–44%" ("rides on this leg alone" when no other leg
  is still to play) — tagged like a straight (this line / same side / other
  side); the ticket's own leg records leave the list. The rail is
  unchanged; the verb is `add` when a straight or a teaser in the math is
  on the row's direction. Same on the Ticket tab, alerts and the Min
  suggested bet filter; the "by my exposure" sort counts teaser dollars.
- **Tests** — `tests/condkelly.test.js`: no ticket unchanged; a ticket with
  no other leg = a half-point straight; one and two other legs match a
  brute-force growth grid; shared other legs are joint; same side lowers,
  other side raises; budget and non-monotone declines; a whole-number
  ticket cut throws. `tests/betsview.test.js`: chips, related lines, verb,
  leg records dropped, other market / other period / no fair notes, a
  decline zeroes the teaser dollars, a ticket on another game ignored.
- **Version control** — same branch `feature/unabated-ticket-teasers`.
  Commits: (1) this section; (2) `condkelly.js` tickets + tests; (3)
  `betsview.js` + `teaser.js` leg reason + tests; (4) `panel.js` chip,
  related line, exposure sort; (5) README + root `CLAUDE.md`. Files:
  `extension/condkelly.js`, `extension/betsview.js`, `extension/teaser.js`,
  `extension/panel.js`, `tests/condkelly.test.js`, `tests/betsview.test.js`,
  `unabated_ticket/README.md`, `CLAUDE.md`, this plan.
- **Worktree** — `.claude/worktrees/buckeye-teasers`; no DuckDB touched.
  After an approved merge: `git worktree remove` + `git branch -d`, restart
  the bets service (launchd `kickstart`), reload the unpacked extension.
  Local main took the tail-flex work (0.15.1) while this was built; it was
  merged into the branch before the merge to main, manifest 0.16.0.
- **Documentation** — README Stake section: an "Open teasers count too"
  paragraph with the 9/27 numbers; the Teasers Limits line; Tests counts;
  a decisions-log entry; root `CLAUDE.md` Unabated Ticket blurb.
