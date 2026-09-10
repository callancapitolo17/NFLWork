# Unabated Ticket — edge scanner + alerts plan

**Date:** 2026-09-10 · **Branch:** `feature/unabated-edge-scanner` · **Worktree:** `.worktrees/unabated-edge-scanner`
**Extends:** `unabated_ticket/` (merged 2026-09-09, issue #111)

## Goal

Two features on one scanner: an **Edges** tab in the side panel listing every
positive-edge line across the leagues you follow, and a **Chrome notification**
when a new one crosses your alert threshold. Covers all leagues at once, not
just the one on screen, with no extra Unabated tabs.

## Data (verified 2026-09-10, no login required)

| Feed | URL | Size | Cadence | Carries |
|---|---|---|---|---|
| Snapshot | `content.unabated.com/markets/v2/league/{lg}/odds.json` | CFB 8 MB, NFL 3 MB, MLB 2 MB gzip; ETag + 304 supported | regenerates ~every 27 s | every row: `sides[sideKey][ms<book>]` line with `points, price, sourcePrice, sourceFormat, bacr, ge, liquidity, statusId, alternateLines`; plus `teams`, `marketSources` (id→name) |
| Updates | `api-k.unabated.com/api/markets/changes/query/{cursor}` | ~300 KB compressed per 10 s | batches every ~1.3 s; cursor = `latestTimestamp` | `gameOddsEvents["lg:pt:pregame"][].gameOddsMarketSourcesLines["si:ms:an"]["bt<N>"]` with the same price fields + `ge`, `bacr`, `sideKey`, `bestAltPoints/Price`, `bage`, `babacr` |

`ge` is Unabated's edge as a fraction (0.0096 = +0.96%), the same number the
screen shows and the ticket sizes from. `bacr` is the fair at the book's
points. The Unabated Line itself (`ms49`) is blurred in anonymous updates but
book lines keep `ge`/`bacr`; inside the extension the fetch carries the login
cookie anyway.

Keys: `lg1` NFL, `lg2` CFB, `lg5` MLB. `pt1` full game, `pt2` 1H, `pt3` 2H,
`pt4–7` quarters. `bt1` moneyline, `bt2` spread, `bt3` total (others are props).
Side 0 = away / Over, side 1 = home / Under. `eventStart` is UTC.

## Architecture

```
unabated_ticket/extension/
  feed.js        pure: parse v2 snapshot + changes batches into one line map   (new, node-tested)
  scanner.js     poll loop, state, edge selection, alert dedupe                (new)
  panel.html/js  + "Edges" tab, alert + scanner settings                      (edit)
  page.js        + publish enabled books / bet types from the page's userSettings (edit)
  background.js  + notification click -> focus Unabated tab                   (edit)
  manifest.json  + host_permissions content.unabated.com, api-k.unabated.com; "notifications", "tabs"
tests/feed.test.js + tests/fixtures/{v2_slice.json, changes_slice.json}  (a few KB each)
```

**Where it runs:** in the side panel page. It is open whenever you're betting,
`setInterval` is reliable there, and it stops the moment the panel closes, so
nothing polls Unabated in the background. (MV3 alarms are capped at 30 s and
would miss update batches; an offscreen document is more machinery than this
needs.)

**Loop:** on panel open, fetch the snapshot for each enabled league (ETag
cached, so a reopen within ~27 s costs nothing), then poll updates every 10 s
from the last `latestTimestamp`, applying each line into
`state[league][eventId][betType][period][sideKey][book]`. Re-fetch the
snapshot every 10 min as a resync. Pause when `document.hidden` and the panel
is not visible.

**Books filter:** your Unabated book selection lives in the page's
`context.userSettings.gameOdds` (and bet types in `oddsFilterContext`).
`page.js` publishes the enabled book ids + selected bet types to storage
whenever a Unabated tab is open; the scanner filters on them. If none has been
published yet the Edges tab says so and shows all books.

**Edge selection:** line on board (`statusId 1`), book enabled, bet type in
1/2/3, period in enabled set (default full game), `ge >= minEdge`, event not
started. Stake per line from the existing `kellyStakeFromEdge`.

## Edges tab

Rows sorted by edge (toggle: stake). Each row: side label (same wording as the
ticket), matchup, book, price `-111 · 52.5¢`, edge, stake, liquidity for
exchanges, minutes to start. Header shows leagues covered, lines scanned, last
update age, and the books filter in effect. Settings: leagues (NFL/CFB/MLB
default all on), periods (FG default), minimum edge to list (default 1.0%).

**Row click:** focuses the Unabated tab (opening `/{league}/odds` if that
league isn't showing), then `page.js` scrolls the grid to the row
(`api.ensureNodeVisible` by grid key pattern) and flashes the cell. You fire
the one-click there, so the deeplink is a real click and never popup-blocked.
If the row is hidden by your bet-type filter the panel says so.

## Alerts

`chrome.notifications` when a line first satisfies `ge >= alertMinEdge`
(default 2.0%). Dedupe key = `(marketId, sideKey, book, points, price)`; a
line re-alerts only if its price improves. Per-event cooldown 5 min so a
moving game doesn't spam. Clicking the notification focuses the Unabated tab
on that league and scrolls to the row (same path as row click). Toggle in
settings; off by default until you've watched the list for a session.

## Phases / commits

1. `feed.js` parsers + fixtures + node tests (snapshot and update slices
   captured 2026-09-10; assert `ge`, `bacr`, `sourcePrice`, side keys, alt
   fields, and that an update overwrites a snapshot line).
2. `scanner.js` loop + Edges tab (list, settings, header stats).
3. Books/bet-types publish from `page.js`; filter wired.
4. Row click → focus tab + scroll + flash.
5. Alerts.
6. README + CLAUDE.md line, pre-merge review, merge.

## Testing

- `node --test unabated_ticket/tests` for `kelly.js` and `feed.js`.
- Scripted Chromium harness (same recipe as the ticket): route the two feed
  URLs to fixture slices, open `panel.html`, assert the Edges tab lists the
  expected rows and that a second update batch changes an edge.
- Manual: open the panel with NFL + CFB enabled, compare three listed edges
  against the Unabated screen; click a row and confirm the scroll/flash; set
  alert threshold low and confirm one notification per line.

## Etiquette / limits

~300 KB per 10 s while the panel is open, snapshot only on open and every
10 min with ETag. Nothing runs when the panel is closed. Same endpoints the
page itself calls, at a far lower rate than its 0.6 s poll.

## Version control

Branch `feature/unabated-edge-scanner`, worktree `.worktrees/unabated-edge-scanner`.
One commit per phase. Pre-merge review, explicit approval, merge to `main`,
`git worktree remove`, `git branch -d`.

## Documentation

`unabated_ticket/README.md`: Edges tab, alerts, settings, data sources, the
books-filter dependency on an open tab. `CLAUDE.md`: extend the Unabated
Ticket bullet.

## Risks

- Unabated could gate `ge`/`bacr` behind login on the public feeds; the
  extension fetches with the session cookie, so it degrades only if they gate
  by referer/origin. Fallback is computing edge from `ms49` `bacr`, which the
  snapshot also carries.
- Update-stream schema differs from the snapshot (nested `bt<N>`, `sideKey`
  on the line); both parsers are fixture-tested.
- Row click across leagues navigates the Unabated tab; if you were mid-scan on
  another league it moves. Acceptable, and the Edges tab is the new home.
