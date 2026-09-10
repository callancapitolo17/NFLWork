# Unabated Ticket

Chrome extension (Manifest V3, plain JS, no build step). Two tabs in one
side panel:

- **Ticket** — click a price on the Unabated odds screen and the panel shows
  the bet (side, points, book price, Unabated fair price, edge) and the
  quarter-Kelly stake. The panel stays open when the sportsbook tab opens,
  so the stake is in view while you place the bet. Issue #111; plan in
  `docs/2026-09-08-unabated-ticket-extension-plan.md`.
- **Edges** — every positive-edge moneyline / spread / total across NFL, CFB
  and MLB at once, read from Unabated's public market feeds while the panel
  is open, with a stake per line, a click that jumps to the row on the
  Unabated tab, and optional Chrome notifications when a new line crosses
  your alert threshold. Issue #112; plan in
  `docs/2026-09-10-unabated-edge-scanner-plan.md`.

One ticket at a time. No bet tracking, no overlay on the Unabated page, no
rounding of the stake.

## Install (load unpacked)

1. Chrome → `chrome://extensions` → turn on **Developer mode** (top right).
2. **Load unpacked** → pick `NFLWork/unabated_ticket/extension`.
3. Pin "Unabated Ticket" from the puzzle-piece menu (optional). Clicking the
   toolbar icon opens the side panel; a captured price also opens it.
4. Open `https://tools.unabated.com/cfb/odds` (premium login) and click a price.

After editing any file under `extension/`, press the reload icon on the
extension's card and reload the Unabated tab.

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

## Stake

Mode B of the Kelly sheet (`extension/kelly.js`): Unabated already publishes
the edge (EV per $1) for every line, so the stake is sized straight from it.

```
b      = decimal(american book price) - 1
full   = edge / b                 0 if edge <= 0
stake  = bankroll * full * multiplier        not rounded
```

The panel shows Unabated's edge % as-is. The fair price (`bacr`) is shown for
information only. Prices print as American plus prediction-market cents
(implied probability), e.g. `-111 · 52.5¢`; on exchanges the cents use the
exchange's exact `sourcePrice` so they match Unabated's screen, while the
stake uses the American price because that is what Unabated's edge was
computed from.

Settings (bankroll, Kelly multiplier) sit at the bottom of the panel and
persist in `chrome.storage.local`. Defaults 30000 and 0.25.

Copy puts one line on the clipboard:
`Seattle Mariners -133 · 57.0¢ @ Novig | fair -139 · 58.2¢ | edge +1.89% | stake $188.55 | Texas Rangers @ Seattle Mariners · MLB`.

## Edges tab

### Data

Two public Unabated feeds (no login needed; verified 2026-09-10). The panel
fetches them only while it is open and stops the moment it closes; nothing
runs in the service worker.

| Feed | URL | What it carries |
|---|---|---|
| Snapshot | `content.unabated.com/markets/v2/league/{1,2,5}/odds.json` (NFL, CFB, MLB; 2–8 MB gzip, regenerated ~every 27 s) | every row's `sides[side][ms<book>]` line: `points, americanPrice, sourcePrice, sourceFormat, bacr, ge, liquidity, statusId, sequenceNumber`; `teams`; `marketSources` |
| Changes | `api-k.unabated.com/api/markets/changes/query[/{cursor}]` (~300 KB per 10 s) | the same fields per changed line under `gameOddsEvents[lg:pt:pregame][].gameOddsMarketSourcesLines[si:ms:an][bt]`, plus `sideKey` |

`ge` is Unabated's edge as a fraction (0.0296 = +2.96%), the same number the
Ticket tab sizes from. `bacr` is the fair at the book's points.

Loop (`extension/scanner.js`): snapshot per enabled league on open and every
10 min (browser-cache revalidation, so a quick reopen is a 304), then the
changes stream every 10 s. The first cursor is derived from the snapshot's
`Last-Modified` (cursor = nanoseconds since 2021-01-06, kept as a string —
it is above 2⁵³) so nothing between the build and the first poll is lost; a
full page (7 batches) is followed immediately; a cursor the server rejects
(`resultCode: Failed`, older than ~3 min) resyncs from snapshots; hiding the
panel pauses polling and a pause over 2 min resyncs on return.

Parsing (`extension/feed.js`, node-tested on real slices under
`tests/fixtures/`):

- A line is keyed `(marketId, book, sideKey)`. Event + bet type is **not** a
  key: team totals reuse `bt3` under the same event.
- An update is applied only when its `sequenceNumber` is newer than the line
  held. The stream replays old lines, and a snapshot can already be ahead of
  a batch.
- Only pregame moneyline/spread/total rows (`pt*:pregame:bt{1,2,3}:e*`).
  Props, team totals (`bt4`) and live rows are skipped.
- Books list only when `isActive && statusId == 1` in `marketSources` —
  what the odds screen itself shows. Dead feeds (Matchbook, pool books)
  carry lines like +5900 at -1.5 with a 3336% "edge".

### What is listed

A line shows when: it is on the board (`statusId 1`), its book is allowed
(see filter), the bet type and period are enabled, `ge` is at or above the
minimum edge, and the game has not started. Rows carry the same wording as
the ticket, the price as American plus cents (exchange cents from
`sourcePrice`), liquidity for exchanges, time to start, and the stake from
`kellyStakeFromEdge` with the panel's bankroll and multiplier. Sort by edge,
stake or start time. Settings (leagues, periods, minimum edge, sort) persist
in `chrome.storage.local` under `edges`.

The header shows leagues loaded, lines held, update age, and the filter in
effect. A league that fails to load is named in a red banner while the rest
keep working; if every league fails the tab says **feed unavailable** rather
than showing an empty list.

**Books filter.** Your Unabated book selection lives in the page's
`context.userSettings.gameOdds` (entries with `isUnavailable === false`) and
the selected bet types under localStorage `oddsFilterContext:preferences`.
`page.js` publishes both every 10 s while an Unabated odds tab is open;
`content.js` stores them as `booksFilter`. Until one has been published the
header says **no books filter yet: showing all live books**. A failed read
keeps the last good filter and reports the error. The bet-type shape was
never seen logged in (it needs a premium session), so an unreadable one is
reported in the header and the list falls back to ML/spread/total — check
the tab's console line `[unabated-ticket] books filter` once and fix
`selectedBetTypeIdsOf` in `page.js` if the ids are somewhere else.

**Row click.** Focuses the Unabated tab showing that league (navigates an
existing Unabated tab, or opens one, when none does), then `page.js` finds
the row through the grid API, scrolls it into view and outlines the price
cell for 2.5 s. You click the price there yourself, so the book's deeplink
is a real gesture and never popup-blocked. If the row is hidden by your
bet-type or period filter the panel says so.

### Alerts

Off by default. Turn on **Notify on new edges at or above N%** (default
2.0%). The first pass after enabling — or after changing the threshold or
the leagues — baselines every line already there without pinging. After
that: one Chrome notification per line the first time it crosses the
threshold, again only if its price improves (dedupe key = market, book,
side, points), and at most one per event per 5 min. Title is the bet and
book, body the edge, stake, matchup and time to start. Clicking the
notification runs the same jump-to-row path as a row click. No alerts
fire while the panel is closed. The alert log lives in `chrome.storage.local`
(`alertLog`, 24 h).

### Etiquette

While the panel is open: ~300 KB per 10 s on the changes stream, a snapshot
per league on open and every 10 min. Same endpoints the page itself calls,
at a far lower rate than its 0.6 s poll.

## Tests

```bash
node --test "unabated_ticket/tests/*.test.js"
```

`kelly.test.js` checks the sheet's worked example (-400 at +12.5% edge,
bankroll 30000, quarter Kelly = $3,750), the Seattle -133 / +1.89% case,
zero or negative edge → $0, that nothing rounds, and that exchange cents
use `sourcePrice`. `feed.test.js` parses the fixture slices (NFL event
125807, captured 2026-09-10): `ge`/`bacr`/`sourcePrice`/liquidity/side keys,
the live-book flag, edge selection and sorting, cursor extraction, and that
an update overwrites a snapshot line only with a newer sequence number.
`scanner.test.js` drives the loop with an injected fetch: cursor from
`Last-Modified`, poll, rejected-cursor resync, per-league failure, resume.

End-to-end without a real login: Playwright (in `mlb_sgp/venv`) with the
ms-playwright Chromium, `--load-extension`, the two feed URLs routed to the
fixtures and `tools.unabated.com` to a page that fakes the grid's React
fiber props; open `chrome-extension://<id>/panel.html` and read its text.

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
- **Edges: "no books filter yet"** with an Unabated tab open: the tab's
  console line `[unabated-ticket] books filter` shows the read error;
  `userSettings.gameOdds` moved. Until fixed the list shows all live books.
- **Edges row click: "row is not on the grid"**: the Unabated tab's own
  bet-type or period filter hides that row, or the game left the board.
- **No notifications**: they only fire while the panel is open and the
  toggle is on; the first pass after enabling is silent by design. Check
  Chrome's notification permission for the extension in System Settings.
