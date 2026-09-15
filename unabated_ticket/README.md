# Unabated Ticket

Chrome extension (Manifest V3, plain JS, no build step). Two tabs in one
side panel:

- **Ticket** — click a price on the Unabated odds screen and the panel shows
  the bet (side, points, book price, Unabated fair price, edge) and the
  quarter-Kelly stake. The panel stays open when the sportsbook tab opens,
  so the stake is in view while you place the bet. Issue #111; plan in
  `docs/2026-09-08-unabated-ticket-extension-plan.md`.
- **Edges** — every positive-edge moneyline / spread / total across every
  team sport Unabated prices (NFL, CFB, NBA, CBB, WNBA, MLB, NHL and ~20
  soccer leagues) at once, read from Unabated's public market feeds while
  the panel is open, with a stake per line, a click that jumps to the row on the
  Unabated tab, and optional Chrome notifications when a new line crosses
  your alert threshold. Issue #112; plan in
  `docs/2026-09-10-unabated-edge-scanner-plan.md`. Alternate spreads and
  totals (every rung Unabated prices, not just the main number) list too
  behind an **Include alt lines** toggle — issue #113, see *Alt lines*.

One ticket at a time. No bet tracking, no overlay on the Unabated page, no
rounding of the stake.

## Install (load unpacked)

1. Chrome → `chrome://extensions` → turn on **Developer mode** (top right).
2. **Load unpacked** → pick `NFLWork/unabated_ticket/extension`.
3. Pin "Unabated Ticket" from the puzzle-piece menu (optional). Clicking the
   toolbar icon opens the side panel; a captured price also opens it.
4. Open `https://tools.unabated.com/cfb/odds` (premium login) and click a price.

After editing any file under `extension/`, press the reload icon on the
extension's card. The service worker re-injects `page.js` and `content.js`
into Unabated tabs that are already open (the old copy retires itself), so
the tab does not need a reload; if the Ticket tab ever shows "No Unabated
tab is running the capture script", reload the tab.

## How capture works

- `page.js` runs in the page's own JS world (`"world": "MAIN"`). The ticket
  is read from React fiber props (`marketLine`, `sideIndex`, `context`) and
  the AG Grid row (`node.data`) attached to the clicked
  `.odds-cell-action-shell`; those are invisible to an isolated content
  script. Capture runs on **pointerdown** (capture phase on `document`, so
  before any Unabated handler) with `click` as a fallback, deduped per cell —
  Unabated's one-click betting can open the book's deeplink on mouse-down and
  may navigate the tab away before a `click` ever fires. Unabated's own
  handler runs untouched afterwards. If the deeplink replaces the Unabated
  tab, the ticket is already stored and the panel shows it with
  "Not watching the line" (no tab left to watch).
- Fast path: fiber props. Fallback: `data-marketline-id` on the shell plus a
  `forEachNode` scan of every row's `sides`. If both fail the panel says
  **Could not read this cell** with the reason; it never shows a stake it
  cannot back.
- Fields: `price` = `americanPrice` (exchanges have only `price`), `fair` =
  `marketLine.bacr` (Unabated's no-vig price at that book's points),
  side 0 = away / Over, side 1 = home / Under.
- The ticket also carries `watch: {gridKey, sideKey, bookKey}` — the row id
  and `sides["si<n>:tid<id>"]["ms<book>"]` path used to find the same line
  again. Not in the plan's contract; needed by the watcher.
- `content.js` (isolated world) writes the ticket to `chrome.storage.local`
  directly. `background.js` is off the hot path — it only sets the
  panel-open-on-click behavior and handles edge-notification clicks — so a
  sleeping or crashed service worker can't stop a click from reaching the
  panel.
- Every 5 s `page.js` re-reads the same book line through the grid API. If
  price or points moved, the panel shows **Line moved**, re-sizes off the
  new price and fair, and keeps the captured line for comparison. Off the
  board shows in red. If the Unabated tab is closed or navigated away the
  panel says **Not watching the line**.
- A click on an **alternate-line cell** works the same way: its `marketLine`
  is one of the main line's `alternateLines`, so the ticket carries
  `watch.altPoints` and the watcher re-finds that rung by points inside
  the ladder (a rung the book pulls shows as off the board). The screen
  computes `edge` for main lines only, so an alt ticket is sized from the
  feed's `ge` on the object — the same number.

## Stake

Mode B of the Kelly sheet (`extension/kelly.js`): Unabated already publishes
the edge (EV per $1) for every line, so the stake is sized straight from it.

```
b      = decimal(american book price) - 1
full   = edge / b                 0 if edge <= 0
stake  = bankroll * full * multiplier        not rounded
```

Under the stake the panel shows **To win** (profit at the book's American
price) and **Payout** (stake plus profit), the way an exchange order slip
does. The panel shows Unabated's edge % as-is. The fair price (`bacr`) is
shown for information only. Prices print as American plus prediction-market cents
(implied probability), e.g. `-111 · 52.5¢`; on exchanges the cents use the
exchange's exact `sourcePrice` so they match Unabated's screen, while the
stake uses the American price because that is what Unabated's edge was
computed from.

Settings (bankroll, Kelly multiplier) sit at the bottom of the panel and
persist in `chrome.storage.local`. Defaults 30000 and 0.25.

Copy puts one line on the clipboard:
`Seattle Mariners -133 · 57.0¢ @ Novig | fair -139 · 58.2¢ | edge +1.89% | stake $188.55 | to win $141.77 | payout $330.32 | Texas Rangers @ Seattle Mariners · MLB`.

## Edges tab

### Data

Two public Unabated feeds (no login needed; verified 2026-09-10). The panel
fetches them only while it is open and stops the moment it closes; nothing
runs in the service worker.

| Feed | URL | What it carries |
|---|---|---|
| Snapshot | `content.unabated.com/markets/v2/league/{id}/odds.json?t=<30 s bucket>` (29 team-sport leagues, `feed.LEAGUES`; ~18 MB gzip in total, CFB alone 9.7 MB, regenerated ~every 27 s). The query is a cache buster: CloudFront hands gzip clients the bare URL from an edge cache that was 6.5 h old on 2026-09-10 | every row's `sides[side][ms<book>]` line: `points, americanPrice, sourcePrice, sourceFormat, bacr, ge, liquidity, statusId, sequenceNumber`, plus its `alternateLines[]` (same fields per rung); `teams`; `marketSources` |
| Changes | `api-k.unabated.com/api/markets/changes/query[/{cursor}]` (~300 KB per 10 s) | the same fields per changed line under `gameOddsEvents[lg:pt:pregame][].gameOddsMarketSourcesLines[si:ms:an][bt]`, plus `sideKey` |

`ge` is Unabated's edge as a fraction (0.0296 = +2.96%), the same number the
Ticket tab sizes from. `bacr` is the fair at the book's points.

Loop (`extension/scanner.js`): snapshot per enabled league on open (four
downloads at a time, each league listed the moment it lands, CFB last),
then the changes stream every 10 s, which covers all leagues in one call.
**The anonymous stream is incomplete**: measured 2026-09-10 over 3 min it
delivered 69 of the 191 NFL line changes the snapshot recorded, and the
misses were Kalshi (42), Caesars (45), ProphetX, Polymarket and Underdog —
the books that carry the edges. So each league's snapshot is re-downloaded
on its own cadence by compressed size (≤2 MB every 60 s, ≤5 MB every 2 min,
larger every 5 min; NFL is ~3 MB, CFB 9.7 MB), plus a full resync every
10 min. Expect roughly 6–8 MB/min with every sport on. The first cursor is derived from the snapshot's
`Last-Modified` (cursor = nanoseconds since 2021-01-06, kept as a string —
it is above 2⁵³) so nothing between the build and the first poll is lost; a
full page (7 batches) is followed immediately; a cursor the server rejects
(`resultCode: Failed`, older than ~3 min) resyncs from snapshots; hiding the
panel pauses polling and a pause over 2 min resyncs on return.

Parsing (`extension/feed.js`, node-tested on real slices under
`tests/fixtures/`):

- A line is keyed `(marketId, book, sideKey)`. Snapshot game rows are unique
  per event, period and bet type, but the changes stream tags other markets
  of the same event with the same `bt` key (team totals under `bt3`), so for
  updates event + bet type is **not** a key; only `marketId` tells them apart.
- An update is applied only when its `sequenceNumber` is newer than the line
  held. The stream replays old lines, and a snapshot can already be ahead of
  a batch.
- Only pregame moneyline/spread/total rows (`pt*:pregame:bt{1,2,3}:e*`).
  Props, team totals (`bt4`) and live rows are skipped.
- Books list when `isActive` and enabled for game odds in `marketSources`.
  `statusId` there is **not** a liveness flag: on 2026-09-11 Caesars,
  Bet365, Fliff, Bet105, BetOnline and Underdog Prediction Market carried
  `statusId` 2 or 3 with lines changed minutes earlier (the first rule,
  `statusId == 1`, hid them from the Books filter). Dead feeds (Matchbook,
  `isActive` false) carry lines like +5900 at -1.5 with a 3336% "edge";
  the max-line-age gate catches the rest.
- Each spread/total line's `alternateLines[]` expands into alt lines keyed
  `(marketId, book, sideKey, points)` — `<main key>:alt<points>` — flagged
  `isAlt` with `mainPoints` = the parent line's points. Measured on the live
  NFL file 2026-09-11 (38k alts, 1.4k with `ge` ≥ 1% at live books): the
  parent's `marketId` and `ms<id>` are the key because an alt's own
  `marketId` can be null (Fanatics) and its `marketSourceId` can name
  another book (Sports Interaction mirrors BetMGM's 4); `stn` is the
  market's standard number, not the book's main (Hard Rock: `stn` 47.5 on
  a 48.0 main); ladders can hold `null` entries; an alt on the main line's
  own points is dropped (same bet twice). Moneylines have no alts.

### Alt lines

Off by default. Tick **Include alt lines** in the filter box and every
book's alternate spreads and totals join the list under the same gates as
main lines (board, book, bet type, period, edge, start, line age) plus two
of their own — most alt "edges" are deep longshots (live 2026-09-11 the
median NFL alt edge sat 13 points off the number at +400 and up; -18.5 at
+800 for +3.6% and a few-dollar stake is typical) where Unabated's fair is
extrapolated:

- **Max pts from main** (default 7): distance from the book's *current*
  main-line points. 0 = no limit.
- **Min liquidity $** (default 100): for lines that report liquidity, i.e.
  exchanges (Kalshi's median alt depth was $129, Novig's $250); books with
  no figure pass. 0 = no limit.

An alt sitting on the main line's current number is hidden (it would be
the same bet twice; when a main line moves onto an alt's number via the
stream, that alt hides until the next snapshot replaces the ladder). A
main line the stream takes off the board leaves its ladder listed until
that refresh too. On the grid, an alt cell's book is read by object
identity in the row's ladders, then the column, before the line's own
`marketSourceId` (Sports Interaction's alts carry BetMGM's id).
Rows carry an `alt` badge and say `alt of -2.5` (the book's main number);
the header counts alts apart (`3,120 lines (+40,278 alts)`).

**Freshness.** The anonymous changes stream carries **no alt updates**
(2,603 keys in a live page, all `an0`, none with `alternateLines`; only a
`bestAlt*` summary rides on the main line), so alts are exactly as fresh as
the league's last snapshot refresh — 60 s / 2 min / 5 min by file size —
and a main line that moves between refreshes leaves its ladder stale until
the next one. Every alt's `modifiedOn` is the feed's `0001-01-01T00:00:00`
sentinel; its `sequenceNumber` is the change time in epoch ms (on 10,182
main lines it trailed `modifiedOn` by a median 1.2 s), so **Max line age**
applies to alts through that (`feed.lineChangedMs`). Caveat: on ~1% of
main lines the sequence ran far ahead of `modifiedOn`, so an alt's age can
read younger than the price really is.

### Grouped by market

On by default (**Group by market** in the filter box). With alts on, one
soft ladder yields 4–10 rows that all say one thing, so the list shows one
**card per (game, period, bet type, side)** — a +EV opinion is directional,
so the two sides of a market are two cards. The card's head is the side and
its best edge; under it the market, matchup and start; then the **best
line**, which is the highest Kelly stake (stake = edge / (decimal − 1)
already taxes longshots), so a -110 main line at +5% outranks a +944 rung at
+6%; then `▸ 2 books · 7 lines (+6)`, which opens the other books and rungs.
Every line is clickable (locate) as before. The count badge counts cards,
sort orders cards through their best line, and `feed.groupEdges` (pure,
node-tested) does the grouping; the panel passes the stake as the rank.
Turning the toggle off gives the flat list.

### What is listed

A line shows when: it is on the board (`statusId 1`), its book is allowed
(see filter), the bet type and period are enabled, `ge` is at or above the
minimum edge, the game has not started, and the book changed the line
within **Max line age** (default 168 h). That last one matters: on
2026-09-10 the biggest "edges" on the live feed were a 96-day-old Buckeye
-110 on a 44.5 total (+36.67%) and a 24-day-old SouthPoint Oklahoma +2 —
dead feeds at books Unabated still flags active. Every row prints the
line's age ("line 12d old") so a genuinely stale line at a live book, which
is a real edge, can be told from a dead one. Rows carry the same wording as
the ticket, the price as American plus cents (exchange cents from
`sourcePrice`), liquidity for exchanges, time to start, and the stake from
`kellyStakeFromEdge` with the panel's bankroll and multiplier. Sort by edge,
stake or start time. Settings (sports, periods, bets, books, minimum edge,
max line age, sort) persist in `chrome.storage.local` under `edges` (sports
as league ids; a sport checkbox toggles all of its leagues).

Leagues come from `feed.LEAGUES` (ids probed 1–70 on 2026-09-10, labels
read off team names). Tennis (ATP 9, WTA 10) and combat (22) are not listed:
their sides key on people and their bet types are not 1/2/3. The row-click
screen path is verified for nfl/cfb/mlb only; nba, cbb, nhl, wnba and
soccer follow the site's nav and, for a soccer league, the odds screen must
have that league selected for the row to be found.

The header shows leagues loaded, lines held, when the newest snapshot was
built (its `Last-Modified`; minutes or more means a stale edge copy),
stream age, and the filter in effect. A league that fails to load is named in a red banner while the rest
keep working; if every league fails the tab says **feed unavailable** rather
than showing an empty list.

**Books and bets.** The filter box has a **Books** dropdown (every live
book in the feed, multi-select, with "My Unabated selection" / "All live" /
"None") and **Bets** checkboxes (moneyline / spread / total). Until you
tick a book yourself the list follows the selection `page.js` reads off
your open Unabated odds tab (`context.userSettings.gameOdds`, published
every 10 s, stored as `booksFilter`); the summary says which source is in
effect, and the header line under the status explains it. Click that
header line to see the selection as read from the tab and the raw fields
it came from (a diagnostic; on 2026-09-10 the `isUnavailable` flag gave 33
books where the screen showed ~12, so the flag may still need adjusting —
your own ticks always win). Settings persist under `edges` (`bookIds` null
= follow Unabated).

**Row click.** Focuses the Unabated tab showing that league (navigates an
existing Unabated tab, or opens one, when none does), then `page.js` finds
the row through the grid API, scrolls it into view and outlines the price
cell for 2.5 s. You click the price there yourself, so the book's deeplink
is a real gesture and never popup-blocked. If the row is hidden by your
bet-type or period filter the panel says so. An **alt row** click first
expands the grid row's Alts (`node.setExpanded(true)`), then finds the
cell by its fiber props — points, side, book and, when the cell's row data
carries them, event and bet type — retrying for 2 s while the alt cells
mount. If no such cell renders it says so and outlines the main-line cell
instead, so the row is still found. The expand-and-match path is verified
against the scripted grid only (see Tests); the real screen's Alts row was
not reachable from the harness, so the first real click is the check.

### Alerts

Off by default. Turn on **Notify on new edges at or above N%** (default
2.0%). The first pass after enabling — or after changing the threshold or
the leagues — baselines every line already there without pinging. After
that: one Chrome notification per line the first time it crosses the
threshold, again only if its price improves (dedupe key = market, book,
side, points — so alt lines, when included, are deduped per rung, and
turning alts or their gates on re-baselines first), and at most one per
event per 5 min. With **Group by market** on, the unit is the card instead:
one notification per (game, market, side) about its best line, again only
when the card's best line improves by the card's own ranking — a higher
stake, so a pulled main line that leaves a +944 rung as "best" is not news
however its edge % compares — with `N books · M lines` in the body;
toggling grouping re-baselines. The alert card is built from the lines at
or above the *alert* threshold, so its best line and counts can differ
from the card on screen (built at the list threshold). In the flat
list a ladder with several rungs over the threshold pings once per 5 min
per rung until each has fired; the notification title says `(alt of -2.5)`
so an alt is never mistaken for the main line. Title is the bet and
book, body the edge, stake, matchup and time to start. Clicking the
notification runs the same jump-to-row path as a row click. No alerts
fire while the panel is closed. The alert log lives in `chrome.storage.local`
(`alertLog`, 24 h).

### Etiquette

While the panel is open: ~300 KB per 10 s on the changes stream and the
per-league snapshot refreshes above (6–8 MB/min with all sports on, ~1
MB/min with just football/baseball/basketball/hockey). Same endpoints the
page itself calls, at a far lower rate than its 0.6 s poll; nothing runs
when the panel is closed.

## Tests

```bash
node --test "unabated_ticket/tests/*.test.js"
```

`kelly.test.js` checks the sheet's worked example (-400 at +12.5% edge,
bankroll 30000, quarter Kelly = $3,750), the Seattle -133 / +1.89% case,
zero or negative edge → $0, that nothing rounds, and that exchange cents
use `sourcePrice`. `feed.test.js` parses the fixture slices (NFL event
125807, captured 2026-09-10; its `alternateLines` are real rungs from the
2026-09-11 file for the same event, one `null` entry added): `ge`/`bacr`/
`sourcePrice`/liquidity/side keys, the live-book flag, edge selection and
sorting, cursor extraction, that an update overwrites a snapshot line only
with a newer sequence number, and the alt path — keys, `mainPoints`,
`includeAlts` off by default, the distance / liquidity / same-number gates
following a moved main line, age through `sequenceNumber`, a pulled
ladder disappearing on the next parse, and `groupEdges` (card keys, book
and line counts, best by edge vs by stake, cards following their best).
`scanner.test.js` drives the loop with an injected fetch: cursor from
`Last-Modified`, poll, rejected-cursor resync, per-league failure, resume,
and a league switch while a snapshot is still downloading.

End-to-end without a real login: Playwright (in `mlb_sgp/venv`) with the
ms-playwright Chromium, `--load-extension`, the two feed URLs routed to the
fixtures and `tools.unabated.com` to a page that fakes the grid's React
fiber props; open `chrome-extension://<id>/panel.html` and read its text.
The panel page's CSP forbids string eval, so every `evaluate` /
`wait_for_function` must be an arrow-function string. The #113 run (14
checks, 2026-09-11) drove: alts off by default, the toggle listing 6 alt
rows under the default gates and 16 with them off, the distance gate,
settings persistence, an alt row click expanding the fake row (a
`setExpanded` that mounts the ladder's shells) and outlining the right
alt cell, a main row click unchanged, an alt-cell ticket with
`watch.altPoints` staying quiet for a watch tick then going off the board
when the rung was pulled, and a main-cell ticket unchanged. The grouped
run (20 checks) adds: 4 cards for 9 lines, the Bears card's best line being
the -110 main over the +944 rung, the badge counting cards, the expander,
a nested row click locating, and a card staying open across a re-render.

Manual checklist after loading unpacked: click a best-line price and a
book-column price, then a moneyline, a spread and a total; confirm side
label, points sign, price, fair and stake; change bankroll and watch the
stake move; wait for a line change and see the warning; close the Unabated
tab and see "Not watching".

## Troubleshooting

- **Panel empty / "Click a price"**: nothing captured yet. Click a price
  again.
- **No Unabated fair for this line**: Unabated has no edge for that line
  (common on lopsided moneylines and exchange-only lines), so there is
  nothing to size from. Not a bug; pick a line that shows an edge %.
- **Could not read this cell**: Unabated changed prop or class names. Check
  `.odds-cell-action-shell` still exists and the fiber props still carry
  `marketLine` / `sideIndex` (see the DOM notes in the plan doc); the
  reason text names which lookup failed.
- **Panel did not open on click**: Chrome only auto-opens the side panel
  with a user gesture attached; click the toolbar icon once, it stays open.
- **"Not watching the line"**: the Unabated tab is closed, navigated away,
  or the row left the grid (filter change). Re-click the price.
- **Cannot size: Unabated has no edge at the new line**: the line moved to
  points Unabated has not priced yet. Wait a tick or re-click.
- **Edges: "feed unavailable for CFB (HTTP 403)"**: Unabated blocked or moved
  the public feed for that league. The other leagues keep updating; check
  the URL in `scanner.js` against what the odds page requests.
- **Edges: "changes cursor rejected; resyncing"**: the laptop slept or the
  panel was hidden for minutes; the scanner reloads the snapshots by
  itself. Persistent repeats mean the cursor format changed
  (`feed.cursorFromDate`).
- **"No Unabated tab is running the capture script"** (Ticket warning or
  Edges header): open an odds tab, or reload the one you have. The worker
  re-injects on an extension reload, but a tab opened before the extension
  was first installed, or one where the re-injection failed, needs a reload.
- **Edges header shows "(N stale: …)"**: those leagues' snapshot builds are
  older than 15 min — a CDN edge copy the cache buster did not get past, or
  an off-season league whose file stopped regenerating (Serie B, WBC).
- **Edges: "books: all N live (Unabated selection unreadable …)"** with an
  Unabated tab open: click that line for the raw fields; `userSettings.gameOdds`
  moved. Tick your books in the Books dropdown meanwhile.
- **Edges row click: "row is not on the grid"**: the Unabated tab's own
  bet-type or period filter hides that row, or the game left the board.
- **Edges alt row click: "row expanded but no … cell at N rendered"**: the
  row was found (its main cell is outlined) but no cell for that book at
  that number mounted within 2 s — Unabated's Alts section renders
  differently from what `page.js` expects (`node.setExpanded`), or the
  book's ladder no longer carries that rung. Open the row's Alts by hand
  and look for the number; if it is there, the expand path needs the real
  DOM (`altCellShellFor`).
- **Edges alt rows look stale**: they only refresh with the league snapshot
  (60 s–5 min); the stream never carries alts. The row's line age comes
  from the alt's `sequenceNumber`.
- **No notifications**: they only fire while the panel is open and the
  toggle is on; the first pass after enabling is silent by design. Check
  Chrome's notification permission for the extension in System Settings.

## Design decisions log (moved from the root CLAUDE.md, 2026-09-15)

History of design decisions that used to live in `NFLWork/CLAUDE.md`. The sections above are the maintained reference; this log records *why* each choice was made and when, with issue numbers.

**Unabated Ticket** (`unabated_ticket/`) — Chrome MV3 side-panel extension (plain JS, no build). A capture-phase click on an Unabated odds-screen price reads the React fiber / AG Grid row (`page.js`, MAIN world), builds a one-at-a-time ticket (side, points, book price, `bacr` fair) in `chrome.storage.local`, computes the unrounded quarter-Kelly stake (`kelly.js`, node-tested), and re-reads the line every 5 s for a line-moved warning. **Edges tab (issue #112)**: while the panel is open, `scanner.js` reads Unabated's public feeds — the v2 league snapshot (`content.unabated.com/markets/v2/league/{lg}/odds.json?t=<30s bucket>` — the query busts a CloudFront edge cache that served gzip clients a 6.5 h-old copy on 2026-09-10 — on open + per-league refresh, 4 at a time, for the 29 team-sport leagues in `feed.LEAGUES`: NFL/CFB/NBA/CBB/WNBA/MLB/NHL + soccer; tennis and combat key sides on people and are out) and the changes stream (`api-k.unabated.com/api/markets/changes/query/{cursor}`, every 10 s; cursor = ns since 2021-01-06, kept as a string; **incomplete anonymously** — 69 of 191 NFL line changes in 3 min on 2026-09-10, exchange moves mostly missing — so each league's snapshot re-downloads on a size-tiered cadence, 60 s / 2 min / 5 min) — through `feed.js` (pure, fixture-tested: lines keyed `(marketId, book, sideKey)` because the changes stream tags other markets of an event with the same `bt` key, updates applied only on a newer `sequenceNumber`, books listed when `isActive` and enabled for game odds — `statusId` is not liveness, Caesars/Underdog carry 2 while live) and lists every ML/spread/total with Unabated's `ge` at or above the minimum and a `modifiedOn` within the max line age (default 168 h — dead feeds at "active" books carried 96-day-old lines with +36% "edges" on 2026-09-10; each row prints its line age), sized with `kellyStakeFromEdge`, filtered by a Books multi-select dropdown + Bets checkboxes in the panel (books default to the selection `page.js` publishes from `userSettings.gameOdds`; the user's own ticks win). Row click / notification click store a `locate` request (`locate.js`) that focuses the Unabated tab and has `page.js` scroll to the row and outline the cell — the price click stays the user's. Alerts are Chrome notifications, off by default, baselined on enable, deduped per line with a 5-min per-event cooldown. **Alt lines (issue #113)**: each snapshot line's `alternateLines[]` expands into lines keyed `(marketId, book, sideKey, points)` under the main line (`isAlt`, `mainPoints` = the book's own main number — NOT the feed's `stn`, which is the market's standard number), listed only behind an **Include alt lines** toggle (off by default) with two alt-only gates, **Max pts from main** (7) and **Min liquidity** ($100, exchanges only), on top of every main-line gate; an alt sitting on its main line's current number is hidden. Freshness caveat measured 2026-09-11: the anonymous changes stream carries NO alt updates and every alt's `modifiedOn` is the `0001-01-01` sentinel, so alts refresh only with the per-league snapshot and their age reads from `sequenceNumber` (the change time in epoch ms). An alt row click expands the grid row's Alts (`setExpanded`) and finds the cell by fiber props (points/side/book/event), failing loudly with the main cell outlined; a ticket captured on an alt cell carries `watch.altPoints` so the watcher re-reads the same rung. **Group by market** (on by default): `feed.groupEdges` collapses the list to one card per (game, period, bet type, side) — a +EV opinion is directional — whose best line is the highest Kelly stake (taxes longshots, so a -110 main at +5% beats a +944 rung at +6%), with `N books · M lines` behind an expander; alerts then key on the card and re-fire only when its best edge improves. No overlay, no tracking, no order placement. Load unpacked from `unabated_ticket/extension`.

