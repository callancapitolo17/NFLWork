# Unabated Ticket

Chrome extension (Manifest V3, plain JS, no build step). Three tabs in one
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
- **Bets** — your own open bets (Kalshi today; BetOnline, Novig and ProphetX
  are separate tickets) read from a local service, matched against the
  board so the Ticket tab says when you already have this line, the other
  side of it, or a bet on the game, and the Edges list flags the same. Issue
  #114; plan in `docs/2026-09-11-issue-114-bet-history-plan.md`. See *Bets*.

One ticket at a time. No overlay on the Unabated page, no rounding of the
stake, no order placement.

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
  tab, the ticket is already stored; the watcher resumes when an Unabated
  tab loads again (Back, a reload), see below.
- Fast path: fiber props. Fallback: `data-marketline-id` on the shell plus a
  `forEachNode` scan of every row's `sides`. If both fail the panel says
  **Could not read this cell** with the reason; it never shows a stake it
  cannot back.
- Fields: `price` = `americanPrice` (exchanges have only `price`), `fair` =
  `marketLine.bacr` (Unabated's no-vig price at that book's points),
  side 0 = away / Over, side 1 = home / Under.
- Edge, three places in order: the clicked object (`edge.edge`, the screen's
  computed %, else the feed's `ge` fraction), then the row's own
  `sides[side][ms<book>]` entry at the same points (the cell's prop can be a
  copy without `ge` — live 2026-09-12 a Novig main line the Edges tab listed
  carried neither on the cell, twice, after a tab reload), then in the
  panel the Edges feed's copy of the same line, matched on game, bet type,
  period, side, book and points — never on the feed key, since an alt rung's
  cell object can lack `marketId`; two feed lines at that number are two
  markets (a team total the changes stream tagged `bt3`) and the panel
  refuses rather than guess — **only at the same price**; an edge is for
  one price. A ticket
  sized from the feed says so in the warning strip with the feed copy's
  age. When none of the three has an edge the panel shows **No Unabated
  fair for this line** with the cell's fields (`noEdgeDetail`) and what the
  feed holds instead.
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
- **The watcher outlives the page.** `page.js` dies with every navigation
  (one-click betting leaving the tab, Back, a tab reload or discard, an
  extension reload's takeover) while the ticket stays in storage. On load
  `page.js` posts `resume_request`; `content.js` answers with the stored
  ticket (and offers it unasked once on its own load, since the two scripts
  load in no guaranteed order) and `page.js` rebuilds the watcher from the
  ticket's identity with no grid API — `watchedRowNode` finds the grid from
  the DOM and the row by event, bet type, period, side, book and number.
  Until 0.6.6 nothing re-attached, so every ticket after such an event read
  "Not watching the line" under a live heartbeat with nothing to re-click
  for. A tab whose grid has not mounted yet reports "no odds grid
  reachable on the page" for a tick or two, then reads. Not resumed: a tab on another
  league (a second tab would otherwise post "wrong league" every 5 s), and
  a ticket whose `eventStart` (naive UTC on the grid row) has passed, and
  an offer older than the capture the tab is already watching. Two same-league tabs may both
  watch; `content.js` drops a failed read while a good one from the last
  7.5 s stands, so a tab on another game date cannot flap the banner.
  Harness: `tests/page_rows.test.js` (resume section).
- A click on an **alternate-line cell** works the same way: its `marketLine`
  is one of the main line's `alternateLines`, so the ticket carries
  `watch.altPoints` and the watcher re-finds that rung by points inside
  the ladder (a rung the book pulls shows as off the board). The screen
  computes `edge` for main lines only, so an alt ticket is sized from the
  feed's `ge` on the object — the same number. A cell in the expanded Alts
  section sits on its own grid row (an AG Grid child: `detail`, `level`
  > 0, or a parent carrying data) whose entry for the book is the rung
  itself, so both capture and the watcher resolve the market's rows
  through one ranking (`rankedMarketNodes`, over every AG Grid on the
  page — an open Alts section mounts its own): a row whose entry for the
  book IS the line being resolved (capture: the clicked object or its
  number, on the entry or in its ladder; the watcher: the captured
  number) beats every shape signal, then a top-level row beats a child, a
  row carrying the market's `bestLines` beats one without, and a row
  carrying this book's `alternateLines` ladder — counted by rung, not
  array length — beats a lone rung. The line-identity rank exists because
  Unabated's alt-lines views list a market's rungs as sibling TOP-LEVEL
  rows sharing one grid key (live 2026-09-12, CFB UAPB@ALCN: eight Over
  rows 33.5 .. 56.5 with `bestLines` and no ladder; NFL CHI@CAR: rung rows
  next to a main row whose 25-rung ladder held the clicked -2.5, the -2.5
  row's own entry carrying an alternateLines array with no rung in it),
  which tie on every shape signal — grid order picked the 33.5 row for a
  click on 56.5, the rungless array counted as a ladder and tied the -2.5
  row with the main row, the Row trace fired on every such capture, and
  the watcher's keyed lookup could answer with any sibling. Capture
  classifies, prices and watches against the best-ranked row that carries
  the book (`ladderRowFor`): on the NFL layout the main row, the ticket an
  alt line; on the CFB layout the rung's own row. The watcher is anchored
  on the captured NUMBER, never on which row answers: each tick trusts
  the grid-key lookup only for a row of the same shape carrying that
  number (its entry at it, or its ladder's rung — one rule,
  `lineOnEntry`), otherwise re-ranks requiring the book's entry as capture
  did, and reads the price at that number; a number gone from the entry
  and its ladder reads off the board at the captured price. So a main
  line whose number moves while the ladder keeps the old number reports
  the price at the captured number, not "moved to" the new one — you bet
  a number. `tests/page_rows.test.js` loads the real `page.js` into a vm
  sandbox (it exposes its row-resolution internals only when the sandbox
  sets `__unabatedTicketExposeInternals`) and replays both live grids, an
  expanded Alts section and a plain main line; run it against an older
  `page.js` with `PAGE_JS=<path>` to watch it reproduce the tie verbatim.
  Picking a child is
  legitimate when a book prices no main line, so the shape is compared,
  never forced. Measured against the
  Alts row the rung read as a main line and the watcher followed the real
  main number (live 2026-09-12: Alabama A&M -5.5 +264 became "now -102 at
  +1.5"); picking "any row of the market carrying a ladder" then followed
  the LOWEST rung (Under 19.5 +265 became "now +2242 at +2.5"). The ticket
  carries `rowResolution` (the candidate rows, the pick, the script
  build) and each watch tick names the row it read and its shape; the
  panel prints that **Row trace** under the warning only on a concrete
  doubt — the watcher reading a row of a different shape than capture
  picked, an alt ticket whose watched number moved (it is re-found by
  number, so it cannot), or two rows tied at the best rank. An ordinary
  line move never shows it. Send that line with a screenshot if a capture
  ever follows the wrong rung again.

## Panel layout

The panel is a fixed header over one scrolling pane per tab. The header holds
the tab bar, the bets header line, and (on the Edges tab) the filter toolbar;
everything else scrolls inside its own tab. The header stops at 60% of the
panel and scrolls itself past that, so the filter drawer and the settings
block — which live in it — stay reachable in a short window instead of
squeezing the pane to nothing. **This is what keeps your place in
the Edges list**: the three tabs used to share the document's scroller, so
hiding one collapsed the scroll height and Chrome clamped `scrollTop` to 0 —
every capture (which brings the Ticket tab forward on its own) sent the list
back to the top. Each pane now scrolls on its own, `showTab` remembers and
restores each pane's offset, and `renderEdges` preserves it across the
scanner's rebuild every few seconds. The row you last clicked keeps a tint, and
the Ticket tab shows a **← Back to edges** link to it.

The filter controls sit behind one chip that states the filter in words
(`Football · FG · Moneyline/Spread/Total · 12 books · ≥1.0%`) and opens the
drawer in place; sort and minimum edge stay out on the toolbar. Bankroll,
Kelly multiplier and the bets service URL live behind the **⚙** at the right
of the tab bar.

An Edges row is two columns: the pick, market, matchup and the book's line on
the left, and a right rail carrying the **edge %** and the **stake**, so both
line up in one column down the list. Edge magnitude also reads as colour in
three tiers (≥4%, 2–4%, under 2%) on the figure and on the row's left stripe,
and the time to first pitch warms to amber inside 12 hours and red inside 2.

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

The figure shown is the number to act on, with the verb on it: `bet $281`
when nothing is held on the market, `add $121` when a position is already
down, and under it what you hold and what full size is
(`$120 held · full size $241`). Against a position on the other side it reads
`bet $192` over `$120 on the other side · net $72 on this side` (the Ticket
spells it `$120 already on the other side`). The Ticket's label says
the same thing ("Bet" / "Add to your position" / "Already at full size"). A
line that cannot be sized keeps a `—`, never a computed-looking `$0`.

Settings (bankroll, Kelly multiplier, bets service URL) sit behind the ⚙ in
the tab bar and persist in `chrome.storage.local`. Defaults 30000 and 0.25.

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
The card is its best line — it carries that line's own number, and a rung
behind the expander names its own when it differs. Every line is clickable
(locate) as before. The count badge counts cards,
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

## Bets

Two kinds of source feed the flags:

| Venue | Source | How it refreshes |
|---|---|---|
| Kalshi | `bets_service/sources/kalshi.py` (local service, signed REST) | every 60 s while the service runs |
| Novig | `bets_service/sources/novig.py` (local service, the account's own Auth0 refresh token, #116) | every 60 s while the service runs; no tab needed |
| Novig (fallback) | `extension/novig_page.js` + `novig_content.js` (content scripts on `app.novig.us`) | whenever the Novig tab's Portfolio screen fetches its lists — open it to refresh |
| BetOnline | — (#115) | shows "no source configured" |
| ProphetX | — (#117) | shows "no source configured" |

Start the service (next section), keep the panel open. Every 30 s while the
panel is visible it fetches `http://127.0.0.1:8094/bets.json` (never from
the service worker), resolves each record's teams through `teams.js`, dedupes
on the venue's native id against what it already holds, and keeps open bets
plus settled ones from the last 30 days in `chrome.storage.local`
(`betsService`; `betsSettings` holds the service URL; `betsNovig` is what
the Novig content script wrote).
A poll that fails keeps the last records and says so; nothing is ever
blanked. A poll that succeeds is authoritative for every venue whose source
reports `ok`: a stored record of that venue the payload no longer lists is
dropped (a reset service DB, a purged fill), so a stale position cannot flag
lines forever; records of a venue whose source failed stay as they were.

**What is matched.** A bet matches a line when the league is the same, the
two teams resolve to the same pair (either order) or the rotation number
matches, and the time agrees: within 30 min when the venue gives a start
time (Kalshi MLB tickers, every Novig order), else the bet's Eastern date
within a day of the line's (Kalshi football tickers carry the date only). A bet that two board
events accept (a series, a doubleheader without a time) is **never**
guessed — it lands in the unmatched list as "ambiguous game". Only open bets
match; settled and closed positions stay in the list but never flag a line.

**Four tiers**, strongest first (`bets.js`, node-tested):

| Tier | Meaning | Ticket banner | Edges badge |
|---|---|---|---|
| `same_line` | same market, period, side and number | "You bet this: Eagles -3.5 -110 · $300 @ Kalshi · Sep 10 2:15 PM" | `held $300` |
| `same_side` | same market, period, side; different number | "You have Eagles -3.5 -110 (this is -4.5)" | `held $300` |
| `opposite` | same market and period, the other side | red: "You are on the OTHER side: Cowboys +3.5 -105 · $200 @ Kalshi" ("at a different number" when the points differ) | `against $200` (red) |
| `same_game` | same game, any other market or period | "You have a bet on this game: Under 40.5 · $150 @ Kalshi" | `game` |

**A bet you hold never hides a line — it changes the size of the next one.**
The edge still being there after you bet it is information (add, or at
least know the market has not moved against you). Per line, `held` is the
dollars risked on the same direction (`same_line` + `same_side`; a different
number is the same opinion) and `against` the dollars on the other side;
both are dollars risked at every venue, so they compare directly with the
Kelly stake (`bets.exposureOf`, `betsview.stakeAdvice`):

The wording is the same three numbers in the same order everywhere (user
choice 2026-09-11): **wagered** what you hold on this side (or "against" when
it is on the other side), **target** the Kelly stake, **bet** the number to
act on. `betsview.stakeAdviceWords` builds it once for the row, the Ticket
block and the Copy text:

| You hold | Stake column | Ticket stake block |
|---|---|---|
| nothing | `$500` | — |
| $300 same side, Kelly $500 | `bet $200` over "wagered $300 → target $500" | "wagered $300 → target $500, bet **$200**" |
| $600 same side, Kelly $520 | `bet $0` over "wagered $600 → target $520" (muted) | "wagered $600 → target $520, bet **$0**" |
| $200 other side, Kelly $500 | `bet $500` over "wagered $200 against → target $500 (net $300 on this side)" | red "wagered $200 against → target $500, bet **$500** (net $300 on this side)" |
| same game only | `$500` | — (the banner still lists the bet) |

Held and against rows (and cards) carry a dim line naming the position:
"you hold Texas A&M -38.5 -110 · $300 · Kalshi · Sep 10 2:15 PM". A `game`
badge is a plain marker — another market on the game does not change how
this line is sized.

A Kalshi NO on a team market is the other team **or a tie** (NFL/CFB/soccer);
it matches as that team and the label says so ("NO Eagles ≈ Cowboys or
tie"). Kalshi first-5 and RFI markets map to the `F5` / `I1` periods.

**Where it shows.**

- *Header line* under the tabs, always: "bets: 14 open · kalshi 20 s ·
  betonline — · novig — · prophetx —" — open bets known to the panel and the
  age of each venue's last successful pull (a dash = no source yet). Red
  when the service is unreachable.
- *Ticket tab*: a banner between the matchup and the stake, one line per
  matching bet, strongest first, at most 5 then "+N more"; nothing when no
  bet matches; under the Kelly stake, the held / add / other-side block
  above. The warning strip adds "Bet sources unavailable" when no venue has
  reported in the last hour (the flags may then be missing).
- *Edges tab*: the badge, the **Related bets** block and the sized stake on
  each row, or on each card from its best line. The block is labelled and
  ruled (red when a position is against you); a bet on that very line prints
  only what differs from the row — venue, its entry price, when — because the
  row already states the pick, while another market or the other side names
  itself. The Ticket shows the same matches with the full sentence. Sort **by my exposure** puts held and
  against lines first. Nothing is filtered; alerts skip only lines you
  already hold at size (nothing to act on) and fire as before otherwise.
- *Bets tab*: opens on **total at risk** across open bets, with the venue
  count and, when a venue reported a bet without a stake, how many are not in
  that total. Then one line per venue (a dot for freshness, what it holds, how
  old the last pull is) rather than a table: last pull green under 5 min, amber
  under 60, red past that or on a failed poll with its error; venues with
  no source yet read "no source configured", a source still on its first
  poll reads "no completed poll yet"; a page-sourced venue past the hour
  says how to refresh — "open app.novig.us and its Portfolio screen in a
  tab to refresh", or "Novig tab is open — open its Portfolio screen to
  refresh" when the tab was seen in the last 5 min; the service itself shows
  "unreachable since …" in red with the last records still listed), the
  open bets (venue, bet, stake, placed — each with a green left edge when
  the board matched it, red when it did not, so a problem bet is visible in
  the open list too), and the **unmatched** list — every open bet no board line matches, with why: team
  not recognised (the raw name, so `teams.js` can grow), ambiguous game, no
  event on the board yet, league not on the scanner, not a game market
  (futures, the bots' combos), unknown Kalshi series.

### Novig source (service)

Novig's official NBX API is a separate "Liquidity Provider" account with a
$30k minimum deposit, so it cannot see bets placed in the retail app. The
service instead logs in **as the retail account once** and polls the same
Hasura GraphQL the app uses (issue #116):

```bash
/Users/callancapitolo/NFLWork/kalshi_draft/venv/bin/python3 -m unabated_ticket.bets_service.sources.novig_auth connect
# or, to avoid the copy/paste race on the callback:
/Users/callancapitolo/NFLWork/mlb_sgp/venv/bin/python3 -m unabated_ticket.bets_service.sources.novig_auth connect --browser
```

- **Connect** runs Auth0's PKCE authorization-code flow against the app's
  public client (no secret exists). You log in yourself; the script only
  needs the URL you land on (`https://app.novig.us/?code=…&state=…` —
  copy it at once, the app strips it as it loads; `--browser` opens a
  Playwright window that intercepts the callback so nothing can be lost).
  The refresh token goes to `NOVIG_TOKEN_PATH`
  (`bets_service/novig_token.json`, gitignored, mode 0600) and the source
  registers on the next service start.
- **Why its own token.** The web app keeps a *rotating* refresh token in
  localStorage, and Auth0 revokes the whole chain when a rotated token is
  reused — copying the app's token would log the app out and kill the
  poller. A separate login is a separate chain; both live side by side.
- **Polling** (`sources/novig.py`, every 60 s): refresh the 30-min access
  token inside a 2-min margin (a rotated refresh token is rewritten at
  once), resolve the trader once from the JWT `sub` (`user.auth_id` →
  `trader_id`, the app's own `AppEntry_Query` chain), then read the `order`
  and `parlay` tables filtered to that trader and to rows that are open or
  changed within the retention window, 100 per page to a short page. The
  selections are the fields the app's Portfolio cards read, so the rows have
  the shape `novig_bets.js` was written against; `normalize_novig()` is a
  port of it and `tests/test_parity_novig.py` holds the two byte-equivalent.
  Transport is `curl_cffi` with Chrome impersonation, the session the
  anonymous Novig SGP scraper already gets through Cloudflare with.
- **Renewal.** Auth0 chains usually cap at about 30 days; when the refresh
  fails the poll fails loudly (the Bets tab row goes red with the error and
  the previous records stay) and you run `connect` again.

### Novig source (content scripts, fallback)

The content-script route reads the same bets without any token, while a
Novig tab has the Portfolio screen open. It stays as the fallback when the
service is not connected; the panel shows the service row when both report.
Two content scripts run on `app.novig.us` and send **zero** requests of their own:

- `novig_page.js` (MAIN world, `document_start`) wraps `window.fetch` before
  the app bundle captures it and mirrors the responses the app fetches for
  its Portfolio screen — the Apollo queries `ActivePortfolioOrders_Query`,
  `SettledPortfolioOrders_Query` and `ParlayPortfolioQuery` POSTed to
  `api.novig.us/v1/graphql` (subscriptions run over a WebSocket and are not
  mirrored; the app refetches these lists over fetch after a fill, which is
  caught). No token is read; the clone of a watched response goes to the
  isolated world by `window.postMessage`.
- `novig_content.js` (ISOLATED world, after `novig_bets.js`) keeps the pages
  each list returned in this tab (keyed on operation + where-clause, offset 0
  restarts a list), normalises them with `novig_bets.js` and writes
  `chrome.storage.local.betsNovig = {bets, readAt, url, error, complete,
  pageSeenAt}`. The panel merges `bets` on every change: a **complete** read
  (all three operations seen in this tab, every list's last page short of
  its limit) is authoritative for the venue — a stored Novig record it no
  longer lists is dropped; an incomplete read (a list still has pages the
  app has not loaded, or an Active-only refetch from another screen) only
  adds.

`novig_bets.js` (pure, node-tested on
`tests/fixtures/bets/novig_bets.json`) applies the app's own rules, read off
its bundle on 2026-09-11: outcome index 0 is the **home** team or **Over**;
`price` is a 0–1 probability and one contract pays $1; `isBid: true` backs
the outcome and `false` lays it, so a lay at `p` becomes the other side at
`1 - p` (a spread's number negated); the matched size is `originalQty -
qty`; `MONEY`/`SPREAD`/`TOTAL` (+ `_1H`, which on MLB is the first five
innings → period `F5`, else `1H`) are game markets, everything else
(props, futures, `1X2`, team totals) is `other` / "not a game market";
leagues map NFL→nfl, NCAAF→cfb, NBA/NCAAB/WNBA/MLB/NHL and the soccer
leagues, the rest fail closed as "league not supported (…)". Status follows
the app's card: settled by the outcome's WIN/LOSS/PUSH against `isBid`, a
cancel with no fill is `void`, a wash or an approved cash-out is `closed`, a
cancel with fills stays `open` on the matched part. A resting or pending
order is `open` on its full size with `approx: ["novig_order_unmatched"]`
(or `_pending`) — it is a bet you are trying to place. Parlays give one
record per leg (`id` `novig:<parlay>:<leg>`, stake = the parlay's wager).
Team keys are left null; the panel resolves `game.awayTeam.name` /
`homeTeam.name` through `teams.js` (the runtime index below — Novig's full
names are Unabated's for every team seen so far, with six aliases).
Novig sends no rotation number, so the team pair is the ONLY way a Novig
bet reaches a line: one unrecognised name blocks the bet outright, and
it lands in the Bets tab's unmatched list naming the spelling to add.

Caveat: the fixture's shapes come from the bundle's operation documents and
the fragments the cards read the `market` / `outcome` / `fills` JSON blobs
through, not from a logged-in capture. If a live blob differs, the record
lands in the unmatched list with a specific "unreadable Novig order (…)"
reason and its `raw` fields — nothing is guessed.

## Bets service

A local Python service that turns the user's own bet history into normalised
records the panel can match against a line (issue #114; plan in
`docs/2026-09-11-issue-114-bet-history-plan.md`). It is the only place that
signs Kalshi requests — the private key never enters the extension. Read-only
GETs; no order placement.

```bash
./unabated_ticket/bets_service/run.sh        # http://127.0.0.1:8094
```

- **Install**: nothing beyond the `kalshi_draft/venv` (duckdb, cryptography);
  `run.sh` uses it when present, else `python3`. Launch from the repo root.
- **Credentials**: `KALSHI_API_KEY_ID` + `KALSHI_PRIVATE_KEY_PATH`, read from
  the environment, then `unabated_ticket/bets_service/.env`, then the bots'
  `kalshi_draft/.env` in the main checkout — so with the bots configured no
  new file is needed. `.env.example` lists every knob (port, retention window,
  Kalshi cadence, log level). Never commit `.env`.
- **Endpoints** (loopback only, no auth): `GET /bets.json[?days=N]` →
  `{generatedAt, sources: {kalshi: {fetchedAt, ok, error, count}}, bets: [...]}`
  with open bets plus settled/closed ones within `N` days (default 30);
  `GET /health` → `{ok, uptimeSec, sources}`.
- **Kalshi source** (`sources/kalshi.py`): every 60 s pulls fills since the
  last poll with a 60 s overlap (deduped on `trade_id`) and unsettled
  positions, plus one cached public GET per market and per event; a full
  fills re-pull once an hour is the reconcile, and it re-reads every cached
  market without a result yet (settlement is the one thing on a market payload
  that changes). Records are one per (ticker, side): positions are the truth
  for the open size, fills give the VWAP entry price and first fill time,
  `market.result` gives won/lost. Team keys are left `null` — the panel fills
  them with `bets.resolveTeamKeys()` so the team table lives only in
  `teams.js`. `normalize_kalshi()` is a port of `extension/bets.js`
  `normalizeKalshi()`; `tests/test_parity.py` holds the two byte-equivalent.
- **Store** (`store.py`, `bets.duckdb`, gitignored): `bets` upserts on the
  record id and is never pruned (the CLV work needs the history);
  `source_runs` appends one row per poll. A failed poll writes a failed
  `source_runs` row and touches nothing else, so a dark source keeps serving
  its previous records; a store write that raises (disk full) is logged and
  retried next poll, never killing the poll thread. Until a source's first
  poll completes (Kalshi: ~1–2 min, one throttled GET per market and event)
  `/bets.json` lists it as `{ok: false, error: "no completed poll yet"}`.
  Log: `bets_service.log` (rotating, 10 MB × 3).
- **Adding a venue** (#115 BetOnline, #117 ProphetX; Novig is
  `sources/novig.py` above): a module in `bets_service/sources/` with
  `name`, `poll_sec` and `fetch() -> list[record]` (the `Source` protocol in
  `sources/__init__.py`), registered in `service.main()`. `fetch()` returns
  every record the venue knows and raises on failure — never a partial list.
  Records follow the contract in the plan (`id` = `"<venue>:<native id>"`,
  `side`/`points` in the side's own number, raw team names, keys `null`). A
  content-script venue instead writes `{bets<Venue>: {bets, readAt, url,
  error, complete}}` to `chrome.storage.local` and the panel merges it with
  `betsview.mergePageSource`.

## Tests

```bash
node --test "unabated_ticket/tests/*.test.js"
/Users/callancapitolo/NFLWork/kalshi_draft/venv/bin/python3 -m pytest unabated_ticket/bets_service/tests
```

`betsview.test.js` covers the panel's bet-history presentation helpers:
freshness colours at the 5 / 60 min bounds, the per-venue rows (unconfigured,
failed poll, never fetched), the service status texts, the "sources
unavailable" rule, the header line, the banner's 5-line cut, the badge
text and kind (held / against / game), `stakeAdvice` (none / add / at size /
reverse with the net), the position lines, the stored + fresh merge (newest
per id, a venue's ok pull authoritative, keys filled, old settled pruned),
the ticket → line shape, and settings sanitising.

`novig_bets.test.js` runs the Novig normaliser on
`fixtures/bets/novig_bets.json`: page bookkeeping (complete vs a full last
page, offset 0 restarting a list, parlay lists keyed on their where-clause),
a matched moneyline bid (index 0 = home, probability → American, stake =
contracts × price), a spread bid and the lay of the same outcome (other
team, negated number, `1 - p`), totals (partial fill sized on the matched
part, a resting order open with the caveat), MLB `_1H` → `F5`, NCAAF → cfb,
the strike fallback in the home perspective, props / unsupported leagues /
blobs without teams failing closed with a reason, the settled grades (bid on
WIN, lay on LOSS, PUSH, void cancel, wash closed), cancel-with-fills /
cash-out / REJECTED / PENDING, parlay legs, and native-id dedupe.
`bets.test.js` then matches those Novig records against the NFL slice: the
moneyline bid (same_line / opposite with Novig labels), the spread bid and
its lay on both signs, a resting Under and a parlay leg flagging the game,
settled orders never matching, and the unmatched reasons. On the Python
side `test_normalize_novig.py` restates those cases for the port,
`test_parity_novig.py` holds it byte-equivalent to the JS, and
`test_novig_source.py` covers the auth (refresh once then cache to the
margin, rotation persisted with 0600, missing/broken token file, callback
state/error parsing) and the source (trader resolved once, pagination to a
short page, no user / auth failure / GraphQL error all raise).

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
and a league switch while a snapshot is still downloading. `bets.test.js`
pins the Kalshi normaliser and the matcher on
`fixtures/bets/kalshi_fixture.json`; the pytest suite covers the service
(fills → positions aggregation, both spread signs and total directions, the
NO-moneyline tie caveat, unknown series failing closed, a failed poll keeping
the previous records, the `/bets.json` shape and `?days=` window, and node
parity on the same fixture).

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
The #114 run (31 checks, 2026-09-11) routes `127.0.0.1:8094/bets.json` to a
payload built from the Kalshi fixture plus synthetic CHI@CAR bets and
drove: the header line and open count, the Bets tab rows (kalshi green,
a red Novig row with its error, two "no source configured"), the open and
unmatched lists with their reasons, the Ticket banner for every tier
(moneyline both sides, NO with the tie caveat, spread same-side and
other-side-at-a-different-number, full-game total as same_game, 1H total
same_line), `held $N` / `against $N` / `game` badges with their position
lines and sized stakes ("bet $200" over "wagered $300 → target $500") on cards and
rows, the exposure sort, the ticket's held / other-side / at-size block,
settings and payload persistence, a stale source turning
the row red and raising the Ticket warning while the banner keeps the last
bets, the service going away (red header, "unreachable since", records
kept) and a reload restoring the stored records. A second script pointed
the panel at the running service: 29 real open positions listed, the bots'
combos as "not a game market", no console errors.

Manual checklist after loading unpacked: click a best-line price and a
book-column price, then a moneyline, a spread and a total; confirm side
label, points sign, price, fair and stake; change bankroll and watch the
stake move; wait for a line change and see the warning; close the Unabated
tab and see "Not watching". With the bets service running and a real
Kalshi position: click that line and see the banner; the other side of it
in red.

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
- **"Not watching the line"**: no Unabated tab is showing the line — the
  tab is closed, on another league, or the row left the grid (filter
  change, game off the board). A tab that merely reloaded or navigated
  and came back resumes on its own within ~5 s; if the message stays with
  the tab open on the right league, read the reason after the colon.
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
- **Bets: "bets service unreachable since …"** (red header, Bets tab
  banner): nothing is listening on the service URL. Start it with
  `./unabated_ticket/bets_service/run.sh` from the repo root and check
  `bets_service.log`; the panel keeps the last records it fetched and
  retries every 30 s. A URL on another port needs a matching
  `host_permissions` entry in `manifest.json` (only `127.0.0.1:8094` ships).
- **Bets: a venue row is red with an error**: that source's last poll
  failed (the text is the exception); the records shown are from its last
  good pull. Kalshi: expired or missing `KALSHI_API_KEY_ID` /
  `KALSHI_PRIVATE_KEY_PATH` (see the service's `.env.example`).
- **Bets: "Bet sources unavailable" on the Ticket tab**: no venue has
  reported in the last hour — the service is down or every source is
  failing — so a missing flag means nothing. Fix the service, not the bet.
- **Bets: unmatched "unknown Kalshi series"**: a game market whose series
  is not in `bets.js` `GAME_SERIES` (and its Python twin in
  `bets_service/sources/kalshi_ticker.py`); the raw ticker is in the list.
  Add the series with its league, bet type and period, with a fixture test.
- **Bets: unmatched "team not recognised (Name)"**: the venue's spelling
  resolves to no team, or to more than one, in the league's index. The
  index is not hand-written: every league snapshot the scanner parses
  carries Unabated's team list (id, name, abbreviation), `teams.js`
  registers it and the panel persists it (`teamsIndex`), so a name resolves
  once its league has been scanned this session or a previous one. Keys are
  `<league>:<Unabated team id>`, the ids the board's own lines carry, so
  the board side never name-matches at all. A venue spelling resolves by
  exact normalised name (St. = State), then a hand alias (`ALIASES` in
  `teams.js`: "UAlbany" → "Albany", "North Carolina State" → "NC State",
  "Southern Mississippi" → "Southern Miss", "Prairie View A&M" → "Prairie
  View", "Louisiana" → "UL Lafayette", "Southeastern Louisiana" → "SE
  Louisiana"), then the name minus a leading
  code token ("PIT Steelers"), then a UNIQUE word-boundary containment
  ("Steelers", "New England", "Middle Tennessee" → "Middle Tennessee State",
  "Grambling St." → "Grambling" only because the leftover is an
  institutional suffix). Two candidates is null, never a guess. If a
  spelling keeps failing, add one `ALIASES` row with a `teams.test.js` case.
  Resist turning an alias into a rule: "A&M" as an institutional suffix
  would map a bet on Texas A&M to Texas on any week Unabated lists Texas
  and not Texas A&M, since the index holds only the teams currently
  playing. Measured 2026-09-12 against 111 open bets, every one of the 16
  cross-spelling resolutions the rules make was correct, and the only
  failures were three names no rule can derive.
- **Bets: unmatched "ambiguous game"**: two board events accept the bet
  (a doubleheader or series without a start time on the venue side). The
  panel refuses to guess; the bet still counts in the header.
- **Bets: unmatched "no event on the board yet" / "league not on the
  scanner"**: the game is not in any loaded league snapshot — untick fewer
  sports in the Edges controls, or wait for Unabated to list it.
