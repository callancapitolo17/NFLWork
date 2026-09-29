# Unabated Ticket: live edges between quarters

## Goal
List Unabated's between-quarters NFL live edges in the Edges tab, sized and
alerted like pregame edges, and make a click on a live price build a full
ticket (fair included).

## What the live screen does (read from Unabated's app code, 2026-09-27)
- The in-game fair is a per-game "fair price set": `status` ready / expired /
  unavailable, `producedUtc`, `checkpointType`, `periodNumber`, and fair
  probabilities per market. It arrives over the logged-in SSE stream
  (`POST <data API>/subscriptions`, then `EventSource`). The free v2 file
  carries none: 0 of ~5,700 live lines had a fair (NFL, MLB, soccer), including
  at the end of Q3 of SNF while the screen showed edges.
- The app enables it for NFL, straight markets, live mode only. It is "ready"
  at a checkpoint (between quarters) and "expired" once play resumes.
- The screen computes each live line's edge in the browser:
  `line.edge = {edge: EV %, linePrice: book price, simulatedPrice: fair
  American at that line's points}`, alt rungs included. The row gets
  `inGameFairPrice = {status, producedUtc, checkpointType, lines}`.
- The polling stream the scanner uses (`api-k.unabated.com/api/markets/changes/query`)
  now returns HTTP 410 ("Polling change subscriptions have been retired").

## Design
Data path: `page.js` (already running on the Unabated tab) → `content.js` →
`chrome.storage.local.liveEdges` → panel. No new network requests.

1. **Live reader (page.js).** Once a second, while the grid holds live rows:
   take every line with a finite `edge.edge` and `statusId` 1 from rows whose
   `inGameFairPrice.status` is `ready` (skip Unabated's own column ms49 and
   final games). Post `live_edges` = {at, league, games (eventId, checkpoint,
   fair produced time, status), rows}. Post an empty list the moment no fair is
   ready, so the panel clears when play resumes. Post on change plus a 5 s
   heartbeat.
2. **Relay (content.js).** Accept `live_edges`; store per league.
3. **Live section (panel).** Above the pregame list, only while the screen
   shows a live game. During a break: rows in the existing row style (side,
   book, price, fair, edge, stake, held-bet tags) with a `live` tag and the
   checkpoint instead of the kickoff time. In play: one muted line ("edges
   return at the next break"). Heartbeat older than 5 s: amber line, stakes
   hidden (the numbers may be gone).
4. **Filters.** The existing Books, Bets, min edge, min stake, min liquidity and
   include-alts settings apply. No new settings.
5. **Alerts.** The existing "Notify at or above X%" setting. One notification
   per game, checkpoint, book, side and line; the next break alerts again.
6. **Ticket on a live click.** Today a live click gets the edge but no fair
   (the ticket reads `bacr`, empty live). Fall back to `edge.simulatedPrice`,
   and label the fair "live, end of Q1".
7. **Row click.** Highlights the cell on the screen (existing locate), matched
   to the live grid row.
8. **Sizing.** Quarter-Kelly off the live edge (existing `kelly.js`). Open bets
   on the same game show on the row but are not netted: conditional Kelly
   needs a live fair ladder. The fair set carries per-market probabilities, so
   that is a follow-up.
9. **Retired stream.** Stop polling the 410 endpoint. Pregame edges keep
   refreshing from snapshots (1-5 min), and the "changes HTTP 410" banner goes.

## Not in scope
- Live edges with no live screen open (subscribing to Unabated's logged-in
  stream from the extension is more intrusive; see the ToS note in memory).
- Other sports: the app enables the in-game fair for NFL only today. The
  reader is league-agnostic and picks up any league Unabated turns on.
- Moving pregame refresh onto the SSE stream.

## Verification
- Unit tests: live reader over a fixture grid (ready, expired, final game,
  ms49, alt rungs, exchange liquidity); panel Live section states; alert keys.
- Playwright harness (existing recipe): fake live grid with fiber props, panel
  shows the Live section, clears on expiry, ticket carries the live fair.
- Live: MNF 2026-09-28 (Eagles @ Bears, 17:15 PT). At the end of Q1 compare
  five panel rows with the screen's cells (edge % and fair).

## Version control
Branch `feature/unabated-live-edges` off local `main` (d1c88097). Commits:
1. plan doc
2. scanner: drop the retired changes poll
3. page.js live reader + content.js relay + tests
4. panel Live section + alerts + tests
5. ticket live fair + locate on live rows
6. docs + manifest 0.13.0

## Worktree
`.claude/worktrees/unabated-live-edges`. Tests run there. You load the unpacked
extension from the worktree for the MNF check. After your approval: pre-merge
review, merge to `main`, remove the worktree and branch.

## Documentation
- `unabated_ticket/README.md`: Live section, in-game fair facts, retired stream.
- `NFLWork/CLAUDE.md`: one sentence on live edges in the Unabated Ticket bullet.
