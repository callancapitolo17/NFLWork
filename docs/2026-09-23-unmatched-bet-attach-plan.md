# Unmatched open bets: flag, attach, learn — design

2026-09-23. Status: approved by Cal (visual mockup approved the same day), building on `feature/unmatched-bet-attach`.
Replaces the "board history + settled-bet CLV" draft (d358ba1): Cal decided the point is catching open bets that match no game, not CLV. Settled bets, a board recorder and closes are out of scope.

## 1. Goal

An open bet the panel cannot match to a game is left out when the panel sizes the next bet on that game (conditional Kelly, #130). Today that bet only gets a red left edge inside the Bets tab. Make it impossible to miss, let Cal attach it to its game by hand, and learn from the attach so the next bet like it matches on its own.

## 2. Facts this rests on (verified 2026-09-23)

- Matcher (`extension/bets.js`): a bet's game is decided by its Kalshi / Novig venue id when the board carries it, else by the name rule — same league, then team pair or rotation with every known team in the game, then start within 30 min (a date-only bet: within a day). Two board events is "ambiguous game", never a guess.
- `unmatchedReasons` lists every open bet with no game and why: `not a game market` and other `unmatchable` records, `league not on the scanner`, `team not recognised (…)`, `ambiguous game`, `no event on the board yet`, and the id-join variants.
- The team crosswalk (`bets.duckdb::team_crosswalk`, #118) maps (venue, league, venue team) to an Unabated team id. It is applied before names and learned automatically only from id joins. The store is insert-only: a held key is never rewritten.
- The service is `bets.duckdb`'s only writer. Writes are guarded by the Host allowlist and a JSON content type.
- Started games stay on the board: 9 MLB events 0 to 6 h past their start still listed pregame rows at 23:16 UTC. So a bet on a game in progress keeps matching; a game drops off only once it is over.
- Live check at 16:00 PT: 3 of 145 open bets unmatched, all Kalshi futures (`KXNFLOROTY`, `KXJOINCLUB`) with nothing to attach to.

## 3. What turns the tab red

An open, unmatched bet **needs a game** when all three hold:

1. It is a game bet: not `unmatchable` (futures, props, unreadable Kalshi markets) and its league is on the scanner.
2. Its game has not started by what the bet knows: `eventStart` in the future, or `eventDate` today or later (Eastern), or no date at all.
3. The miss is fixable: `team not recognised`, `ambiguous game`, or **start time differs** — the same team pair is on the board within 12 h of the bet's start at another time (a new diagnosis: today it reads "no event on the board yet").

Everything else unmatched is **not on the board**: futures and props, leagues the scanner does not cover, games not posted yet, games already over and awaiting settlement. Those are listed, folded, and never flag. 12 h is below the gap between two games of one MLB series (16 h or more) and above any doubleheader gap.

## 4. Attach

1. **Pick the game.** Board events in the bet's league whose Eastern date is within a day of the bet's date (no date: from the day it was placed to 14 days later). Events where one of the bet's teams already resolves, or its rotation matches, come first and say so; the rest sort by start. Typing searches every event of the league by team name. A bet with no league searches every league.
2. **Confirm the teams.** Each team the venue names is lined up with the event's team on the same side; **Swap** flips it (neutral sites, a one-team bet on the other side). Each name is tagged `known` (already resolves to that team, nothing written), `learn` (resolves to nothing) or `fix` (resolves to another team, shown as "was …"). A "Your bet" line restates the bet on the picked game so a wrong side is visible before saving.
3. **Attach and learn** writes, through the service, one pin (bet id → event id) and the `learn` / `fix` crosswalk rows, marked as taught by that pin.
4. **Undo** deletes the pin and the crosswalk rows it taught. A row it overwrote is not restored; an id join can relearn it.

Rules: a manual row overwrites a held automatic one for the same key (Cal is the authority). Automatic learning stays insert-only, so it never overwrites a manual row. Pins are never pruned; they are a few a week and the history of manual matches.

## 5. Components

### 5.1 Service (`bets_service/`)

- `bet_pins(bet_id PK, venue, league, event_id, event_start, away_team_id, home_team_id, away_team_name, home_team_name, pinned_at)`, upsert on `bet_id`.
- `team_crosswalk` gains `pinned_bet_id VARCHAR` (null on automatic rows).
- `POST /pins.json {pin, crosswalk: [rows]}`: 404 when the bet id is unknown; upserts the pin and replaces the rows, in one transaction. `DELETE /pins.json?betId=`: removes the pin and its rows. Both reply `{ok, pins, crosswalk}` and keep the crosswalk routes' guards.
- `/bets.json` adds `pins`.

### 5.2 Matcher (`extension/bets.js`)

- `applyPins(records, pins)` puts `pin` on each pinned record (and the pin's league on a record that has none).
- `resolveGame` tries the pin first: its event on the board is the game (`join: "pin"`); a pinned event that left the board falls through to the usual rules.
- `unmatchedReasons(bets, lines, now)` adds `needsGame` per the rule in section 3 and the "start time differs" diagnosis.
- Automatic crosswalk learning skips pinned bets (the attach already taught them).

### 5.3 Picker logic (`extension/attach.js`, new, pure)

`attachCandidates(bet, lines, {query, now})`, `attachPlan(bet, event, {swapped, lines})` and `pinRequest(bet, event, plan)`. No DOM, node-tested.

### 5.4 Panel

- Bets tab label turns red with a red count of bets needing a game; the header line says "N not matched to a game".
- Bets tab: a red banner, a "Needs a game" section above Open (hidden when empty) with Attach on each row, and a folded "Not on the board" list below Open. Unmatched open rows keep their red edge; pinned rows carry an `attached` tag and Undo.
- The attach panel opens under its row (the two steps of section 4) and survives the 30 s poll re-render, search focus included.
- Crosswalk rows taught by an attach read "attached by you".

## 6. Not in scope

Settled bets, a board recorder, CLV (decided 2026-09-23). A pin on a settled bet is kept but not used.

## 7. Version control, tests, documentation

- Branch `feature/unmatched-bet-attach`, worktree `.claude/worktrees/board-history-clv`.
- Commits: (1) this document; (2) service pins, routes, pytest; (3) matcher pin join, flag rule, `attach.js`, node tests; (4) panel UI and CSS, docs.
- Docs in the same commits: `unabated_ticket/README.md` (service routes; Bets tab), root `CLAUDE.md` Unabated Ticket bullet.
- `./unabated_ticket/check.sh` green before any diff is shown. Pre-merge review per `CLAUDE.md`; Cal merges. After merging, Cal restarts the bets service and reloads the extension.

## 8. Risks

- A wrong attach teaches a wrong name. Step 2 shows every row it will write and the bet restated on the game; Undo removes both.
- Records whose venue id changes between polls would lose their pin. Record ids are the venue-native bet ids the store already upserts on.
- The panel keeps working against an older service: no `pins` in the payload means no pins, and Attach reports the failed POST.
