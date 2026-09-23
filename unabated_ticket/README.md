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
into Unabated tabs that are already open, so the tab does not need a reload.
Both halves hand over rather than stand aside: the new `page.js` posts a
`takeover` and the old one retires; the new `content.js` calls the retire hook
its predecessor published on the isolated world's `window` (Chrome keys that
world on the extension id, which a reload does not change, so both copies see
it). Until 2026-09-15 `content.js` guarded on a boolean instead, which made the
*new* copy return and left the dead one — whose `chrome.storage` writes throw
"Extension context invalidated" — holding the page's messages; the tab did
need a reload, and this paragraph said otherwise. If the Ticket tab still shows
"No Unabated tab is running the capture script", reload the tab (a tab the
re-injection could not reach still needs one).

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
- **"Unabated changed" (the load-time shape check).** Every fiber read in
  `page.js` already reports its own failure where it is used — a failed
  capture prints the cell's actual keys, the watcher posts its error, the
  book-selection read posts the `userSettings` shape — but all of those only
  fire once you click. So `page.js` runs one check per load, on the odds screen
  only (a path segment `odds`; other Unabated pages may render a grid with no
  price cells at all): on a grid that has rendered `.ag-row`s, at least one
  `.odds-cell-action-shell` must carry React props with `marketLine` +
  `sideIndex` on them. A failing check is re-run once 3 s later and published
  only if it still fails, since AG Grid can mount a row a frame before React
  mounts its price cells. Zero is the new-bundle
  state and the panel shows one red banner above the tabs naming what it
  found (no price cells at all, or cells whose props no longer say what a
  ticket needs, with the shapes seen up the fiber). A board with **no rows**
  is not checked at all — a quiet slate, a non-odds page and a grid still
  loading are indistinguishable from each other, so the panel stays silent
  rather than guess; the check retries for 20 s waiting for rows and then
  gives up, leaving the stored verdict on `checking`. It is deliberately one
  count and not a per-prop matrix: the point-of-use errors above already name
  every other shape. The verdict is re-sent (not re-checked) with every 10 s
  heartbeat, since the extension-reload path injects `page.js` and
  `content.js` in two round trips and a single load-time post can land before
  `content.js` listens; `checking` is never re-sent. It re-runs on every
  `page.js` load (so every navigation),
  and the banner is gated on the capture script still heartbeating, so a
  verdict never outlives the tab that produced it. One verdict is stored, not
  one per tab: with two Unabated tabs open, the last one to load wins, so a
  second tab landing on a page with no rows can clear a real banner. The
  point-of-use errors still fire on the next click.
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
  panel says **Not watching the line**. The moved line is written to
  `watchStatus.current`, **not onto the ticket**: the watcher used to
  read-modify-write the whole `ticket` to hang `current` off it, and its
  `capturedAt` guard tested the ticket it had *read*, so a click landing
  inside that `get`→`set` window was written and then reverted ~1 ms later by
  the old ticket coming back — a stale ticket with a dead watcher, reading as
  "my click did not register". The two writers now touch different keys, so a
  capture always stands; every reader matches `watchStatus.capturedAt` against
  the ticket's, which is what drops a late tick for a previous capture. A
  failed read carries the last moved line forward, so the stake keeps sizing
  off it under **Not watching the line** instead of reverting to the captured
  price the market already left.
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
  line move never shows it. It is collapsed to its reason ("Row trace: two
  grid rows tied for this market"); open it for the candidates and the pick, and
  send that with a screenshot if a capture ever follows the wrong rung again.

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
drawer in place; sort and minimum edge stay out on the toolbar. Bankroll and
the Kelly multiplier sit at the foot of the Ticket tab under **Sizing**, where
the stake they size is; the **⚙** at the right of the tab bar jumps there from
any tab. The filter drawer also has **Min suggested bet $**: zero is off;
otherwise it hides a row unless its current actionable bet or top-up meets the
amount, and alerts skip the rows it hides. The bets service URL is on the Bets tab, with the venues it feeds.

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

On an exchange line (Kalshi, Novig — any line Unabated marks `sourceFormat 4`,
a probability) the panel also prints the **limit order** the stake means,
under the dollar figure: `466 contracts @ 53¢ · $246.98`
(`kelly.contractOrder`). A contract costs its price in whole cents and pays
$1, so the limit price is the exchange's exact probability rounded to the
cent — what you type into Kalshi — and the count is `floor(stake / price)` in
integer cents, never rounded up past Kelly; the cost shown is what the
contracts actually spend, so the leftover (under one contract) is visible
rather than hidden. The count divides the number to act on, so a top-up shows
the top-up's contracts, and a stake under one contract reads
`under 1 contract @ 53¢` rather than a zero. Sportsbook lines show no
contract row. Fees are not folded in: the stake and Unabated's edge are both
pre-fee, so the real outlay on Kalshi is a touch above the cost shown.

That is the stake with nothing held. **With open bets on the same market of
the game, the stake is sized GIVEN them — conditional Kelly (issue #130,
`extension/condkelly.js`)** — instead of sizing alone and subtracting dollars,
which compares dollars at different prices and numbers:

```
K      = bankroll * multiplier      the scale the held bets were sized on
choose x >= 0 to maximize  sum over outcome rows of  prob * ln(1 + pnl / K)
```

- **The new bet's win chance** comes from Unabated's edge, as above:
  `p = (1 + edge) / decimal`. **A held bet's** comes from Unabated's fair
  today at its number: the median `bacr` across books at that half-point rung
  (`extension/ladder.js`, snapshot lines only — the changes stream files team
  totals under the game total's bet type). Its payoff is `toWin` / `-stake`
  from the record. No fees are added: Unabated's prices already include them.
- **Same market, same period: exact.** Spreads, moneylines and every alt
  number are cuts on the margin (away minus home: an away bet at `a` wins above
  `-a`, a home bet at `h` wins below `h`, a moneyline is the `±0.5` cut);
  totals are cuts on the total. Outcome rows are the gaps between the cuts,
  their chances differences along the ladder. One procedure covers the same
  line, another number on the same side, the other side, middles, and a
  moneyline against a spread. A whole-number line gets a push row from the
  `k-0.5` and `k+0.5` rungs (the whole-number rung's own fair is conditional on
  no push and is never read). A held bet on the row's own number takes the
  row's chance, so it needs no rung.
- **Same market, another period, same direction (1H Over with FG Over):
  worst case.** Unabated says nothing about how two periods move together and
  no correlation is estimated (user decision 2026-09-18), so each period's
  rows are sorted by P&L and paired by cumulative probability — bad with bad,
  good with good. Hold 1H Over $300 at +110 (55%), new FG Over +120 (50%), K
  $8,000: $667 alone, $408 worst case (the measured NFL link 0.69 would give
  $530; $408 keeps 95% of the growth).
- **Same market, another period, the OTHER direction (1H Under with FG Over):
  not in the math** — grey `other period · not sized`, the stake is the
  standalone one (user decision 2026-09-19). It is a hedge in the real world,
  and the worst-case pairing would size it as though both bets won together:
  hold 1H Under $300 at +110, same FG Over: $407, where the measured link wants
  $802 ($407 keeps 76% of the growth, $667 keeps 97%). Left out is "no credit
  for the hedge" without the penalty.
- **Another market (a spread against a total): not in the math** — it stays a
  grey `game` line. Open question in #129.
- **Left out and named, never guessed:** no rung at the bet's number or a rung
  whose fair repeats a neighbour's (Unabated flat-lines deep tails; the two
  moneyline cuts are exempt), no stake on the record, parlay legs, a Kalshi NO
  moneyline (also wins on a tie), a soccer three-way moneyline, quarter lines,
  a period Unabated has no ladder for (`F5`, `I1`). A ladder that crosses by
  more than half a point of probability, held bets that can already lose `K`,
  a whole-number row without its two rungs, or an edge so large on a heavy
  favourite that `p >= 1`, decline the whole calc and the standalone stake
  stands.
- The search is a bounded golden section (the score is a single hill), capped
  so `K + pnl` stays positive, and `$0` when the slope at zero is not positive.
  With nothing held the panel uses `kellyStakeFromEdge` directly, so the
  number is unchanged to the cent. `bacr` is a whole American price, so
  conditional stakes are good to about ±$10.
- **Not sized: hedges.** A price with no edge of its own is `$0` even when a
  held bet on the other side would make the math want some (deferred, user
  decision 2026-09-19).

Measured on four live cards, 2026-09-17/18, K $8,000:

| Card | Before | Now |
|---|---|---|
| New Mexico @ Oklahoma: Over 53.5 +217, hold Under 39.5 $214 and Under 48.5 $6.20 | bet $183, "still $37 against" | bet $296 |
| Lions @ Bills: Over 61.5 +213, hold Over 61.5 $270 and Under 51.5 $413 | at full size | add $188 |
| UConn @ Southern Miss: Under 52.5 +125, hold Under 47.5 $181.50 | add $116 | add $121 |
| Chargers @ Bills: Chargers ML +213, hold Chargers +3.5 $400 | bet $189 | add $0 |

The figure shown is the number to act on, with the verb on it — `add` when
anything in the math is held on this direction, `bet` otherwise — and under
it one small line, what the stake would be alone: `add $188.32` over
`$270.05 alone`. The Ticket's label says the same thing ("Bet" / "Add to your
position" / "Already at full size") and its small line adds the position:
`held $270 · against $413 · $270.05 alone`. A line that cannot be sized keeps
a `—`, never a computed-looking `$0`.

Settings (bankroll, Kelly multiplier) sit at the foot of the Ticket tab under
**Sizing**, reachable from any tab via the ⚙, and persist in
`chrome.storage.local`. Defaults 30000 and 0.25. The bets service URL is on
the Bets tab, under the venue strip it feeds.

Copy puts one line on the clipboard:
`Seattle Mariners -133 · 57.0¢ @ Novig | fair -139 · 58.2¢ | edge +1.89% | stake $188.55 | to win $141.77 | payout $330.32 | Texas Rangers @ Seattle Mariners · MLB`.
When held bets changed the number the stake reads `stake $188.32 (add $188.32, $270.05 alone)`.
On an exchange line the contract order follows the stake: `stake $188.55 | 330 contracts @ 57¢ | to win ...`.

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
- **Venue ids on the rungs** (#118 step 2, measured on the live NFL and CFB
  files 2026-09-15). Alt lines keep the book's own `sourceKey` /
  `sourceData` strings. Kalshi rungs carry `sourceKey`
  `Y-KXNCAAFSPREAD-26SEP19DUQWSU-WSU36` (contract side + market ticker; the
  other side reads `N-…`; 12,318 of 12,318 rungs fit that shape), Novig
  rungs `sourceData` = the Novig outcome id (a UUID), ProphetX a 32-hex id.
  Main lines never carry either field, and the changes stream has neither,
  so they refresh with the snapshot only. Each event collects them in
  `event.venueIds`: `kalshiEventSuffixes` (e.g. `["26SEP19DUQWSU"]`, kept
  whole — Kalshi's team codes are not Unabated's abbreviations),
  `kalshiContracts` and `novigOutcomes` (id → `{lineKey, mainKey,
  points, sideIndex}` — `points` / `sideIndex` are the contract's own strike
  and Unabated side, fixed for the id; `lineKey` is the listed line at that number — the main line
  when a Kalshi or Novig rung sits on the main number (priced or not), else
  null for an unpriced rung; the map rebuilds only with the snapshot while
  the stream can move a main line, so a join checks `points` against the
  line's current number). An id
  of any other shape is "no id", never an error. Coverage that day: 95 of
  319 NFL/CFB/WNBA events had Kalshi ids, 98 Novig ids. `describeLine` hands
  every row its event's map by reference (`row.venueIds`), which is how the
  bet matcher joins on them (Bets → What is matched).

### Alt lines

Off by default. Tick **Include alt lines** in the filter box and every
book's alternate spreads and totals join the list under the same gates as
main lines (board, book, bet type, period, edge, start, line age) plus
one of their own — most alt "edges" are deep longshots (live 2026-09-11 the
median NFL alt edge sat 13 points off the number at +400 and up; -18.5 at
+800 for +3.6% and a few-dollar stake is typical) where Unabated's fair is
extrapolated:

- **Max pts from main** (default 7): distance from the book's *current*
  main-line points. 0 = no limit.

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
`kellyStakeFromEdge` with the panel's bankroll and multiplier, never more
than the line's resting liquidity (the rail then reads "all $17 liq").
There is no liquidity floor (removed 2026-09-23): a flat dollar floor hid
longshots whose whole Kelly stake is small, and before then it gated alts
only, so a Novig Portland Fire +809 main moneyline with $17 behind it
listed with "add $32.58". Instead the capped stake is the bet you can
place, and **Min suggested bet** hides the ones too small to bother with,
at any odds. Liquidity is read as dollars you can stake (not verified
against payout) and comes from the league snapshot only (the changes
stream carries none), so it can lag a price move by up to a snapshot cycle. Sort by edge,
stake or start time. Settings (sports, periods, bets, books, minimum edge,
minimum suggested bet, max line age, sort) persist in `chrome.storage.local` under `edges` (sports
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
book in the feed, multi-select, with "Default books" / "My Unabated
selection" / "All live" / "None") and **Bets** checkboxes (moneyline /
spread / total). Until you tick a book yourself the list uses the **default
books** (`DEFAULT_BOOK_NAMES` in `panel.js`, the user's list of 2026-09-15):
Bet105, BetOnline, BetOnline Direct, Bookmaker, Bookmaker-Internal, Buckeye,
Kalshi, Novig, NoVig-Internal, Poly US Ing, Polymarket, Polymarket US,
Prophet Exchange, Underdog Prediction Market — matched by name, because four
of them never appear in the anonymous feed their ids could be read from; a
book the feed does not list that day ticks nothing. "My Unabated selection"
switches the list to the selection `page.js` reads off
your open Unabated odds tab (`context.userSettings.gameOdds`, published
every 10 s, stored as `booksFilter`); the summary says which source is in
effect, and the header line under the status explains it. Click that
header line to see the selection as read from the tab and the raw fields
it came from (a diagnostic; on 2026-09-10 the `isUnavailable` flag gave 33
books where the screen showed ~12, so the flag may still need adjusting —
your own ticks always win). Settings persist under `edges` (`bookIds`
absent = default books, null = follow Unabated, an array = your ticks; an
install that stored null before the default books existed keeps following
Unabated).

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
(`betsService`, which also carries the service's team crosswalk;
`betsSettings` holds the service URL; `betsNovig` is what the Novig content
script wrote).
A poll that fails keeps the last records and says so; nothing is ever
blanked. A poll that succeeds is authoritative for every venue whose source
reports `ok`: a stored record of that venue the payload no longer lists is
dropped (a reset service DB, a purged fill), so a stale position cannot flag
lines forever; records of a venue whose source failed stay as they were.

**What is matched.** First by **venue id** (#118 step 3), exact or nothing:
a Kalshi bet joins the board event whose Kalshi rungs carry its event-ticker
suffix (`venueIds.eventTicker` after the first "-", the WHOLE string —
`26SEP19DUQWSU`, never split into team codes; a moneyline joins through the
suffix its event's spread/total rungs carry, so only when Kalshi lists a
spread or total ladder for that game — else it takes the name rule), a Novig spread/total bet the
event whose Novig rungs carry its `outcomeId`. The event must be in the bet's
league. No team name has to resolve. An id join is final: the name rule is
not consulted, so team names that point at another event cannot win, and a
bet's team key the joined game does not have is ignored for its side. One
keyed team is enough to place a bet: its own team's key names that side, the
other team's key the opposite one. An id
on two board events is "ambiguous game (Kalshi event … on 2 board events)";
an id on no board event of its league (Kalshi lists ladders on a third of CFB events,
Novig moneylines have no rungs at all, a number the ladder dropped) falls
through to the name rule below. Malformed or missing ids are no id, never an
error. BetOnline records carry no ids. The join decides the GAME; the tier
still reads the bet's own market, side and number against the row's CURRENT
line (so a -35.5 bet on a line the changes stream has since moved to -36.5
is `same_side`), and the id map's `lineKey` / `mainKey` are never read — the
map is as old as the last snapshot, and a Novig lay's outcome id names the
side the bet is against. All three surfaces — the Edges rows, the Ticket banner
and the unmatched list — decide the game on the whole board (one row per
event), so a bet whose id event has no listed edge row never flags a listed
row its team names happen to fit. A spread bet whose team names do not resolve is
placed on its side off its own contract in the map (Kalshi `Y-`/`N-` market
ticker, Novig outcome): the same side at the contract's strike, the other
side at the negated one (a Kalshi NO, a Novig lay). A Kalshi moneyline with
unresolved names has no contract on a rung, so it matches as `same_game` —
until the crosswalk below has keyed its teams.

**Team crosswalk (#118 step 4).** Every id join is also a lesson: the bet's
venue names both teams (Novig by its own team id, `awayTeamVenue.id`; Kalshi
has no team ids, so its event-title name — "PIT Steelers" — is the key) and
the joined board event carries both Unabated team ids. The panel learns
`(venue, league, venue team) → Unabated team id` for both teams
(`bets.learnCrosswalk`), sends the new rows to the bets service in one
`POST /crosswalk.json` after every poll and every board update, and the
service keeps them in `bets.duckdb::team_crosswalk` with `learned_from`
(the bet id and board event) and `learned_at`, serving the whole table with
every `/bets.json`. Records then resolve their keys through the crosswalk
BEFORE `teams.js` sees the names (`bets.resolveTeamKeys(records,
crosswalk)`), so a later bet of that venue on either team matches even
where the venue lists no ladder for the game (nothing to id-join) — and a
Kalshi moneyline whose names resolve nowhere gets its side. Fail-closed:
learned only from a join on exactly ONE board event whose row carries both
Unabated ids, from a bet naming both venue teams; a bet whose own contract
on the joined row (a Kalshi `Y-`/`N-` ticker or a Novig outcome at the bet's
number) sits on a different Unabated side than the bet's venue side is a
conflict — the venue's away/home is not Unabated's for that game, the one
orientation check that needs no name; a venue name that already
resolves to a DIFFERENT team than the joined event's is a conflict — neither
side is learned and the panel logs why once (`console.warn`) — because a
swapped away/home (a neutral site) would otherwise write a wrong row; a
held venue team is never rewritten as another id (the service refuses it
too and reports it in the POST reply). A wrong row could still only match a
game where the opponent and the start time also agree. The *Bets tab* lists
the table under **Team crosswalk** — "Wazzu → Washington State · Novig ·
CFB · learned Sep 15 3:00 PM", the source bet in the tooltip — with a
**Clear** button (two clicks: the first arms it as "Clear 68 rows?" for 6 s,
the second deletes — no dialog) that `DELETE`s the table on the service and
re-keys every record from its names alone; the board teaches the rows again
as bets join. Live 2026-09-15 (60 open records, 358 board events): 68 rows
would be learned (54 Novig, 14 Kalshi), 0 conflicts, and re-learning
against the table taught nothing new; every spelling learned that day
already resolved by name (step 1's eventName spellings), so the payoff is
the next venue spelling `teams.js` cannot key, not a rescue today.

Otherwise a bet matches a line when the league is the same, the
two teams resolve to the same pair (either order) or the rotation number
matches, and the time agrees: within 30 min when the venue gives a start
time (Kalshi MLB tickers, every Novig order), else the bet's Eastern date
within a day of the line's (Kalshi football tickers carry the date only),
else — for a venue that gives no game date at all (BetOnline's report,
flagged `approx: game_date_unknown`) — an event starting between 12 h before
and 14 days after the bet was placed, keyed on the rotation number; a
rotation match also needs every recognised team name in the row's game (the
same rotation comes round the next week, inside a dateless bet's 14-day
window), and a team name
`teams.js` does not know still matches by rotation and takes its side from
the row's away/home rotations. A bet that two board events accept (a series,
a doubleheader without a time, two weeks of the same rotation) is **never**
guessed — it lands in the unmatched list as "ambiguous game". Only open bets
match; settled and closed positions stay in the list but never flag a line.

**Six tiers**, strongest first (`bets.js`, node-tested). Spreads and
moneylines are one market here (the margin), totals the other:

| Tier | Meaning | Tag | In the stake's math |
|---|---|---|---|
| `same_line` | same bet type, period, side and number | `this line` | yes, exact |
| `same_side` | same market and period, same direction: another number, or a moneyline against a spread on the same team | `same side` | yes, exact |
| `opposite` | same market and period, the other direction | `other side` (red) | yes, exact |
| `related_same` | same market and direction, another period (1H total on an FG total row) | `same side` | yes, worst case |
| `related_opposite` | same market, other direction, another period | `other period · not sized` (grey) | never: a real-world hedge the worst case would penalise |
| `same_game` | another market of the game (a spread on a total row) | `game · not sized` (grey) | never (#129) |

Every match also carries its **position** on the market axis
(`bets.positionOf`: axis, period, cut, direction, stake, toWin) or the reason
it has none; the side comes from the bet's team keys, then its venue contract,
then its rotation — Unabated's frame, never the venue's away/home.

**A bet you hold never hides a line — it changes the size of the next one.**
The edge still being there after you bet it is information (add, or at
least know the market has not moved against you). The size is conditional
Kelly against the bets that are in the math (see **Stake**;
`betsview.stakeAdvice`). Per line, `held` is the dollars in the math on the
row's direction and `against` the dollars on the other one.

| You hold | Badge | Stake column | Ticket stake block |
|---|---|---|---|
| nothing | — | `bet $500.00` | "Bet" · **$500.00** |
| the same over $270 and two unders $413 at another number | `held $270` `against $413` | `add $188.32` over "$270.05 alone" | "Add to your position" · **$188.32** · "held $270 · against $413 · $270.05 alone" |
| the team's +3.5 $400, row is its moneyline | `held $400` | `add $0.00` over "$188.92 alone" (muted) | "Already at full size" · **$0.00** |
| a bet with no fair at its number, or a parlay leg | `held` / `against`, bare | `bet $500.00` | "Bet" · **$500.00** |
| another market only | `game` | `bet $500.00` | "Bet" · **$500.00** |

Both badges show when both exist. To win and Payout describe the number shown
above them, so a top-up prices the top-up and a $0 line shows no payout at
all; the Copy line carries the same figure. Rows and cards carry a labelled
**Related bets** block, one line per position, bets in the math first: a tag
and the bet itself, "Chattanooga -5.5 +138 · 42.0¢ · $168 · Kalshi", plus
`· now -6.5` when the line has moved off the number you bet. A **coloured**
tag (`this line` / `same side` / `other side`) means the bet is in the math; a
**grey** tag with grey text means it is not, and the tag says why —
`game · not sized` for another market, `other period · not sized` for the
other direction in another period, `no fair at 36.5`, `parlay leg`,
`no stake on the record`, `Kalshi NO also wins on a tie`,
`three-way moneyline`, `quarter line`, `side not resolved`, or the reason the
whole calc was declined (`ladder not monotone`). Capped at three with
"+N more on this game". The Ticket banner shows the same lines and tags (five,
then "+N more"). No placed-at: it never told you which bet was which, and it
was a third of the line (user decision 2026-09-14). The Bets tab still shows
it, where the bets are the subject.

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
  bet matches; under the stake, the position and the standalone size
  (`held $270 · against $413 · $270.05 alone`). The warning strip adds "Bet sources unavailable" when no venue has
  reported in the last hour (the flags may then be missing).
- *Edges tab*: the badge, the **Related bets** block and the sized stake on
  each row, or on each card from its best line. The block is labelled and
  ruled (red when a position is against you); a bet on that very line prints
  only what differs from the row — venue, its entry price, when — because the
  row already states the pick, while another market or the other side names
  itself. The Ticket shows the same matches with the full sentence. Sort **by my exposure** puts held and
  against lines first. Nothing is filtered; alerts skip only lines whose
  stake against what you hold is $0 (nothing to act on) and fire as before
  otherwise. The Min suggested bet filter reads the same number.
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
  (futures, the bots' combos), unknown Kalshi series. A bet with a venue id
  names both tiers that missed: "by id: Kalshi event 26SEP19DUQWSU not on
  any board ladder; by name: team not recognised (…)" ("Novig outcome" for
  Novig; "only on another league's board ladder" when the id sits on an
  event of a different league; no id tier for a Novig moneyline or a league
  off the scanner) — and, last, the **Team crosswalk** list with its Clear
  button (above).

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
names are Unabated's, or the long form of its `eventName`, for every team
seen so far, with four college aliases).
Novig sends no rotation number, so the team pair is the ONLY way a Novig
bet reaches a line: one unrecognised name blocks the bet outright, and
it lands in the Bets tab's unmatched list naming the spelling to add.

Every record (orders and each parlay leg, matchable or not) also keeps
Novig's own ids (#118 step 2): `venueIds: {marketId, outcomeId, eventId,
gameId}` and `awayTeamVenue` / `homeTeamVenue` `{id, name, shortName,
symbol}` — each a string as Novig sent it, else null. The outcome id is
what Unabated's Novig rungs carry as `sourceData` (see Edges → Parsing);
on a lay it is the outcome laid, the side taken stays in `side`. Novig's
`symbol` is not Unabated's abbreviation (`UTC` vs `CHT`, `EKU` vs `EKY`):
stored, never a key. The matcher joins on `outcomeId` (Bets → What is
matched). Live
2026-09-15: all 41 open Novig game bets sat on events with Novig rungs,
and 29 had a rung at their own number and side; the other 12 were 2
moneylines (moneylines have no rungs) and 10 numbers the ladder no longer
listed.

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
  `{generatedAt, sources: {kalshi: {...}, betonline: {fetchedAt, ok, error, count}}, bets: [...], crosswalk: [...]}`
  with open bets plus settled/closed ones within `N` days (default 30) and
  the whole team crosswalk (newest first); `GET /health` → `{ok, uptimeSec,
  sources}`; `POST /crosswalk.json` with `{rows: [{venue, league,
  venueTeamKey, unabatedTeamId, venueTeamName?, unabatedTeamName?,
  learnedFrom?}]}` (at most 1000 rows, 1 MiB) → `{ok, learned, conflicts,
  crosswalk}` — INSERT only, a held key with another id is returned in
  `conflicts` and never rewritten; `DELETE /crosswalk.json` → `{ok, cleared,
  crosswalk: []}`. The two write routes require `Content-Type:
  application/json` (415 otherwise): the service sends no CORS headers, so a
  web page can only reach it with a "simple" cross-origin request (a form or
  `text/plain` POST, which is refused) and never with JSON or DELETE (both
  need a preflight the service never answers), while the extension page is
  exempt from CORS for its `127.0.0.1:8094` host permission — nothing new
  in the manifest.
- **Host rule** (every verb, #125): a request whose `Host` header is not
  `127.0.0.1:<port>` or `localhost:<port>` is refused with 403 before any
  route runs. This is the guard the CORS and Content-Type rules above cannot
  be: a page at `evil.example` whose DNS flips to `127.0.0.1` (rebinding) is
  **same-origin** with this service, so no CORS applies to it at all and it
  could otherwise `GET /bets.json` — every venue's open positions, stakes and
  venue ids — or `DELETE /crosswalk.json`. A browser sends the name the page
  was loaded from, so the rebound request's Host is still `evil.example`.
  Chrome prompts on public→loopback; Safari and Firefox do not. The port
  comes from the listening socket, so a non-default `BETS_SERVICE_PORT`
  guards itself (on port 80 the bare names pass too — browsers omit the
  default port).
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
  Records carry `venueIds: {marketTicker, eventTicker}` (#118; `eventTicker`
  null when the market payload is missing — never derived from the ticker).
  The event ticker's suffix is the string Unabated's Kalshi rungs carry
  (join it whole; CFB codes are not Unabated abbreviations). Live
  2026-09-15: all 7 open Kalshi game bets' suffixes were on the board, each
  on one event and the same one the team names matched; 6 had their exact
  contract as a rung (the 7th a moneyline). The step-3 audit the same day
  (every open record through the matcher against the live NFL + CFB
  snapshots, 288 board events): Kalshi 7/7 game bets matched by name alone
  and 7/7 by id + name, every one joined by id, and across the 2,590 rows
  those bets matched the tier and match set were identical both ways; Novig
  (re-run on the branch's service, which sends outcome ids) 39/41 by name
  and 39/41 by id + name (the 2 misses are WNBA, not in the audited
  snapshots), 29 joined by outcome id, 0 of 11,012 rows different;
  BetOnline 3/3, 0 of 1,190 rows different. Where names resolve the id join
  changes nothing; it pays off on names teams.js cannot key.
- **BetOnline source** (`sources/betonline.py`, #115): every 300 s pulls the
  account's paged bet-history report (`POST api.betonline.ag/report/api/report/get-bet-history`,
  pure HTTP — the same endpoint `bet_logger/scraper_betonline.py` uses) for the
  last 31 days and keeps **pending** bets. The Keycloak refresh token is the
  `krefresh` cookie in `bet_logger/recon_betonline_cookies.json` (main
  checkout; `BETS_BETONLINE_COOKIES_PATH`), written by `bet_logger/recon_betonline.py`
  — run it with `--interactive` and log in by hand when the poll reports
  "token refresh failed" (the token dies after 3 days unused). The access
  token is refreshed only within 60 s of expiry, under an exclusive lock on
  `recon_betonline_cookies.json.lock` shared with the bet_logger scraper and
  its LaunchAgent, and the rotated refresh token is written back atomically.
  Record ids are `betonline:<TicketNumber>-<WagerNumber>`; the league comes
  from `bet_logger/utils.py parse_sport` (the report names the SPORT —
  "FOOTBALL" — never the league); a total names both teams, a spread or
  moneyline only its own team, placed by rotation parity (odd = away, `approx:
  side_from_rotation_parity`); the report carries no game date or settle time,
  so `eventStart`/`eventDate` are null and a settled bet's `closedAt` is its
  placed time. Unknown periods, sports outside the scanner and parlays whose
  legs do not parse fail closed as unmatchable with the reason. Same Game
  Parlay rows have not been seen live yet; their leg grammar is a guess the
  parser refuses rather than misreads.
- **Store** (`store.py`, `bets.duckdb`, gitignored): `bets` upserts on the
  record id and is never pruned (the CLV work needs the history), but only
  rows whose content actually CHANGED are written (#125): a source re-sends
  every record it has ever seen on every poll (Kalshi: ~111 records every 60 s, of
  which ~109 are identical) and a DuckDB update is copy-on-write, so an
  identical rewrite appends a new version and leaves the old one as garbage —
  that is what grew the WAL to 9.8 MB against a 4.7 MB file holding 500 rows.
  The compare key is a `content_hash` column (sha256 of the record with
  `sourceFetchedAt` removed — that field is the poll's own clock and would
  defeat the skip; rows from before the column are rewritten once), checked
  by the UPSERT's own `DO UPDATE … WHERE` — no read-back; measured 0 WAL
  bytes over 20 identical polls of the live 528 records, ~74 ms a poll. So a
  record's served `sourceFetchedAt`, like `last_seen_at`, is the last poll
  that CHANGED it, not the last poll that saw it — deliberately: the panel
  breaks merge ties on it (`bets.js dedupeByNativeId`), and restamping it
  with the latest poll would also mark records the source no longer returns
  (Novig and BetOnline pull a rolling window; `bets` serves open rows
  forever) as fresh. No column records when a venue last returned a record.
  `source_runs` appends one row per poll and is **pruned to
  `BETS_SOURCE_RUNS_RETENTION_DAYS` (default 7; 0 disables)** — it grows
  ~2,000 rows/day and `source_status()` full-scans it twice under the store
  lock on every `/bets.json` and `/health`. Every source always keeps its
  latest run and its latest successful run whatever their age (a source
  failing for weeks still shows its last good read). The prune runs on a
  source poll, at most hourly, and can never fail the poll (an error is
  logged and retried an hour later). `team_crosswalk` (#118 step 4,
  primary key `(venue, league, venue_team_key)`, plus `venue_team_name`,
  `unabated_team_id`, `unabated_team_name`, `learned_from`, `learned_at`)
  holds what the panel learned, INSERT-only and cleared on DELETE. A failed poll writes a failed
  `source_runs` row and touches nothing else, so a dark source keeps serving
  its previous records; a store write that raises (disk full) is logged and
  retried next poll, never killing the poll thread. Until a source's first
  poll completes (Kalshi: ~1–2 min, one throttled GET per market and event)
  `/bets.json` lists it as `{ok: false, error: "no completed poll yet"}`.
  Log: `bets_service.log` (rotating, 10 MB × 3).
- **Adding a venue** (#117 ProphetX; BetOnline and Novig are
  `sources/betonline.py` / `sources/novig.py` above): a module in
  `bets_service/sources/` with `name`, `poll_sec` and `fetch() -> list[record]`
  (the `Source` protocol in `sources/__init__.py`), registered in
  `service.main()`. `fetch()` returns every record the venue knows and raises
  on failure — never a partial list. Records follow the contract in the plan
  (`id` = `"<venue>:<native id>"`, `side`/`points` in the side's own number,
  raw team names, keys `null`). A venue with no game date sets `eventStart`
  and `eventDate` null, carries `rotation` and `approx: ["game_date_unknown"]`,
  and the matcher windows on `placedAt` (see Bets). A content-script venue
  instead writes `{bets<Venue>: {bets, readAt, url, error, complete}}` to
  `chrome.storage.local` and the panel merges it with
  `betsview.mergePageSource`.

## Tests

One command runs everything and exits non-zero if any part fails:

```bash
./unabated_ticket/check.sh
```

It runs, in order, ESLint over `extension/` and `tests/` (`npm run lint`),
the node suite (`npm test` = `node --test tests/*.test.js`, 240 tests) and
the bets service's pytest suite (109 tests, on the `kalshi_draft/venv`
python from the main checkout, resolved the way `bets_service/run.sh`
does, else `python3`). All three run even when an earlier one fails, so one
run shows every failure. ESLint comes from `unabated_ticket/package.json`
— run `npm install` once per checkout (it lands in the gitignored
`node_modules/`; `package-lock.json` is tracked so everyone gets the same
eslint). The config (`eslint.config.js`) errors only on real defects —
undefined or unused variables, unreachable code, dead assignments, a
rethrow that drops its cause — and has no style rules; it knows the
modules are dual-loaded (plain `<script>` publishing `globalThis.UnabatedX`
in the panel, `require()` in tests), so no file needs a disable comment.
The extension itself has no build step and stays loaded unpacked.

`betsview.test.js` covers the panel's bet-history presentation helpers:
freshness colours at the 5 / 60 min bounds, the per-venue rows (unconfigured,
failed poll, never fetched), the service status texts, the "sources
unavailable" rule, the header line, the banner's 5-line cut, the stored +
fresh merge (newest per id, a venue's ok pull authoritative, keys filled, old
settled pruned), the ticket → line shape, settings sanitising, and the #130
stake advice through the real matcher on a Kelly bankroll of $8,000: nothing
held equals `kellyStakeFromEdge` exactly; Lions @ Bills adds $188.32 with
`held $270` + `against $413` badges, and a held Bills -8.5 leaves it
unchanged as a grey `game · not sized` line; the Chargers moneyline with the
+3.5 held adds $0; a 1H over held sizes the FG over at $408.39 while a 1H
under is left out (`other period · not sized`, the standalone $666.67); the same
line held needs no ladder ($116.10, then $0 once over size); every guard (no
rung, flat rung, unknown period, parlay leg, no stake) leaves the bet out
with its reason and a bare badge; a non-monotone ladder and a whole-number
row without its two rungs decline the calc; a price with no edge is $0.

`condkelly.test.js` pins the solver on the cards measured 2026-09-17/18 (K
$8,000): the standalone stake to the cent at four prices, same line $116.10
and $0, UConn $121.28 (easier number held) and $24.10 (harder number), New
Mexico @ Oklahoma $295.99, Lions @ Bills $188.32, Chargers $0, the
cross-period worst case $408.39 against $666.67 alone and no hedge credit
across periods, the whole-number push row, the 0.005 clamp and the decline
past it, held risk above K, and loud failures on a missing rung, a quarter
line and a bad probability. `ladder.test.js` covers the fair ladder: totals
and margin cuts (away `-a`, home `h`, moneyline `±0.5`), the median across
books with both sides counted, the flat-tail guard and its moneyline-pair
exemption, no rung / whole number / other period, changes-stream lines
ignored (`fromSnapshot`), soccer moneylines feeding no cut, and the period
name → id map.

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
cash-out / REJECTED / PENDING, parlay legs, native-id dedupe, and the
venue ids (#118: kept on every record and parlay leg, a missing or
non-string id null). `bets.test.js` joins synthetic Kalshi and Novig records
with unrecognisable team names on `fixtures/v2_venue_ids_slice.json`'s ids
(#118 step 3): the Kalshi spread YES and NO and the moneyline through the
suffix, the Novig bid and its lay on the outcome id, a main line moved off
the id's number by the changes stream (`same_side`), a suffix on two events
(ambiguous), id vs team names pointing at different events (id wins),
malformed / missing / other-league ids, the two-tier unmatched wording, and
BetOnline staying on the name path. Its crosswalk cases (#118 step 4): a
Kalshi join teaching both event-title names with `learnedFrom`, a Novig join
keying on Novig's team ids (name fallback when the ids are missing, nothing
when a side has no name), nothing learned from a name match / an ambiguous
suffix / a settled bet / a row missing an Unabated id, a name resolving to
another team refused as a conflict (and a held row with another id), a
later moneyline on a no-ladder board matching through the crosswalk alone,
`rekeyRecords` taking the keys back on clear, malformed served rows
skipped, and `venueTeamOf`'s id → name → record-name order.
`betsview.test.js` adds the payload crosswalk keying a record before its
names in both merges and the Bets-tab rows; the pytest suite covers the
store (learn once, no duplicates, a different id refused with what is
held, served with the payload, cleared), the row validator naming the first
bad row, and over HTTP the POST / conflict / DELETE round trip plus the
415 on non-JSON writes, 400 on bad bodies and 404 on other paths.
`bets.test.js` then matches those Novig records against the NFL slice: the
moneyline bid (same_line / opposite with Novig labels), the spread bid and
its lay on both signs, a resting Under and a parlay leg flagging the game,
settled orders never matching, and the unmatched reasons. On the Python
side `test_normalize_novig.py` restates those cases for the port,
`test_parity_novig.py` holds it byte-equivalent to the JS (on the fixture,
and on extra rows with malformed or missing ids fed to node on stdin), and
`test_novig_source.py` covers the auth (refresh once then cache to the
margin, rotation persisted with 0600, missing/broken token file, callback
state/error parsing) and the source (trader resolved once, pagination to a
short page, no user / auth failure / GraphQL error all raise).

The #130 smoke run (2026-09-19, a scratch Playwright harness, not committed)
loaded the unpacked extension with the NFL slice at a future kickoff and four
synthetic CHI@CAR bets: both badges, `add $95.90` over `$85.36 alone`, a grey
`no fair at 36.5` tag, a Bears moneyline sizing the Bears -20.5 row as the
same side, grey `game · not sized` lines on the spread rows, and the Ticket
tab showing the same $104.87 as its Edges row, with no console errors.

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
ladder disappearing on the next parse, the rung venue ids on
`fixtures/v2_venue_ids_slice.json` (real CFB rows of 2026-09-15 with the
ids intact: per-event suffixes and id maps, the Novig rung on its main
number pointing at the main line, malformed ids ignored, the changes
stream leaving them alone), and `groupEdges` (card keys, book
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
same_line), `held $N` / `against $N` / `game` badges with their Related
bets block and sized stakes ("add $200" over "$300 held · full size $500") on
cards and rows, the exposure sort, the ticket's held / other-side / at-size block,
settings and payload persistence, a stale source turning
the row red and raising the Ticket warning while the banner keeps the last
bets, the service going away (red header, "unreachable since", records
kept) and a reload restoring the stored records. A second script pointed
the panel at the running service: 29 real open positions listed, the bots'
combos as "not a game market", no console errors.

Pre-merge checklist for this directory, on top of the repo-wide review in
`CLAUDE.md`: `./unabated_ticket/check.sh` passes on the branch (lint, node
suite, pytest — this is the one command the review runs; a failing or
skipped check blocks the merge), then the manual pass below after loading
the branch unpacked.

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
- **Red banner above the tabs, "Unabated changed: …"**: the load-time shape
  check found a grid with rows whose price cells no longer carry what a
  ticket is read from — a new Unabated bundle. Nothing will capture until
  `page.js` is updated for it. The banner names what was found: *no price
  cells* means the `.odds-cell-action-shell` class moved (capture, the
  watcher and locate all find their cell by it); *N price cells … but none
  carries the React props* means the props moved, and the detail prints the
  props shapes up the fiber so the new names are right there. Absence of the
  banner is not a pass on a board with no rows — nothing is checked there.
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
  carries Unabated's team list (id, name, abbreviation) AND, on each game
  row's `eventName`, a second, long-form spelling of both teams ("Prairie
  View A&M Panthers - PRV @ Baylor Bears - BAY" where the list says
  "Prairie View"; NFL/NBA/NHL "Ravens Baltimore BAL", MLB "Dodgers Los
  Angeles" — `feed.teamSpellingsFromEventName`, #118). `teams.js` registers
  both spellings under the same team id and the panel persists them
  (`teamsIndex`, one entry per team with `name` and `eventName`), so a name
  resolves once its league has been scanned this session or a previous
  one. Keys are `<league>:<Unabated team id>`, the ids the board's own
  lines carry, so the board side never name-matches at all. A venue
  spelling resolves by exact normalised name (St. = State) against either
  spelling, then a hand alias (`ALIASES` in `teams.js`: "UAlbany" →
  "Albany", "North Carolina State" → "NC State", "Southern Mississippi" →
  "Southern Miss", "Louisiana" → "UL Lafayette"), then the name minus a
  leading code token ("PIT Steelers"), then a UNIQUE word-boundary
  containment ("Steelers", "New England", "Middle Tennessee" → "Middle
  Tennessee State", "Prairie View A&M" → "Prairie View A&M Panthers",
  "Grambling St." → "Grambling" only because the leftover is an
  institutional suffix). Two candidates is null, never a guess. If a
  spelling keeps failing, add one `ALIASES` row with a `teams.test.js` case.
  Resist turning an alias into a rule: "A&M" as an institutional suffix
  would map a bet on Texas A&M to Texas on any week Unabated lists Texas
  and not Texas A&M, since the index holds only the teams currently
  playing. Measured 2026-09-12 against all 338 distinct (venue, league,
  name) pairs in 408 bet records with the day's NFL/CFB/NBA/MLB/WNBA
  snapshots registered: 58 cross-spelling resolutions, every one correct,
  and 4 unresolved — Kalshi MLB "A's" and "Chicago WS", Novig "Arkansas
  Pine Bluff" and "Texas A&M Commerce" (Unabated hyphenates both). The
  eventName spellings retired the "Prairie View A&M" and "Southeastern
  Louisiana" aliases; the four kept were each re-checked against the same
  snapshot (the long forms are "Albany Great Danes" and "Southern Miss
  Golden Eagles", and "North Carolina State" hits both the Wolfpack and
  "North Carolina" + State).
- **Bets: unmatched "ambiguous game"**: two board events accept the bet
  (a doubleheader or series without a start time on the venue side), or —
  "ambiguous game (Kalshi event … on 2 board events)" — the bet's venue id
  sits on two board events. The panel refuses to guess; the bet still counts
  in the header.
- **Bets: unmatched "by id: … not on any board ladder; by name: …"**: the
  bet carries a Kalshi / Novig id no board rung carries (the venue lists no
  ladder for the game, or the ladder dropped that number), and the name rule
  then missed for the reason after "by name". Fix the name reason; the id
  joins by itself once Unabated lists the ladder.
- **Bets: unmatched "no event on the board yet" / "league not on the
  scanner"**: the game is not in any loaded league snapshot — untick fewer
  sports in the Edges controls, or wait for Unabated to list it.

## Design decisions log (moved from the root CLAUDE.md, 2026-09-15)

History of design decisions that used to live in `NFLWork/CLAUDE.md`. The sections above are the maintained reference; this log records *why* each choice was made and when, with issue numbers.

**Unabated Ticket** (`unabated_ticket/`) — Chrome MV3 side-panel extension (plain JS, no build). A capture-phase click on an Unabated odds-screen price reads the React fiber / AG Grid row (`page.js`, MAIN world), builds a one-at-a-time ticket (side, points, book price, `bacr` fair) in `chrome.storage.local`, computes the unrounded quarter-Kelly stake (`kelly.js`, node-tested), and re-reads the line every 5 s for a line-moved warning — **the watcher resumes from the stored ticket** whenever `page.js` loads (it dies with every navigation: one-click betting leaving the tab, Back, reload, discard, extension-reload takeover; `content.js` hands the stored ticket back on `resume_request` and once on its own load, and `page.js` rebuilds the watch by identity with no grid API — before 0.6.6 every such ticket read "Not watching the line" under a live heartbeat). **Edges tab (issue #112)**: while the panel is open, `scanner.js` reads Unabated's public feeds — the v2 league snapshot (`content.unabated.com/markets/v2/league/{lg}/odds.json?t=<30s bucket>` — the query busts a CloudFront edge cache that served gzip clients a 6.5 h-old copy on 2026-09-10 — on open + per-league refresh, 4 at a time, for the 29 team-sport leagues in `feed.LEAGUES`: NFL/CFB/NBA/CBB/WNBA/MLB/NHL + soccer; tennis and combat key sides on people and are out) and the changes stream (`api-k.unabated.com/api/markets/changes/query/{cursor}`, every 10 s; cursor = ns since 2021-01-06, kept as a string; **incomplete anonymously** — 69 of 191 NFL line changes in 3 min on 2026-09-10, exchange moves mostly missing — so each league's snapshot re-downloads on a size-tiered cadence, 60 s / 2 min / 5 min) — through `feed.js` (pure, fixture-tested: lines keyed `(marketId, book, sideKey)` because the changes stream tags other markets of an event with the same `bt` key, updates applied only on a newer `sequenceNumber`, books listed when `isActive` and enabled for game odds — `statusId` is not liveness, Caesars/Underdog carry 2 while live) and lists every ML/spread/total with Unabated's `ge` at or above the minimum and a `modifiedOn` within the max line age (default 168 h — dead feeds at "active" books carried 96-day-old lines with +36% "edges" on 2026-09-10; each row prints its line age), sized with `kellyStakeFromEdge`, filtered by a Books multi-select dropdown + Bets checkboxes in the panel (books default to the selection `page.js` publishes from `userSettings.gameOdds`; the user's own ticks win). Row click / notification click store a `locate` request (`locate.js`) that focuses the Unabated tab and has `page.js` scroll to the row and outline the cell — the price click stays the user's. Alerts are Chrome notifications, off by default, baselined on enable, deduped per line with a 5-min per-event cooldown. **Alt lines (issue #113)**: each snapshot line's `alternateLines[]` expands into lines keyed `(marketId, book, sideKey, points)` under the main line (`isAlt`, `mainPoints` = the book's own main number — NOT the feed's `stn`, which is the market's standard number), listed only behind an **Include alt lines** toggle (off by default) with one alt-only gate, **Max pts from main** (7), on top of every main-line gate (there is no liquidity floor since 2026-09-23: every suggested stake is capped at the line's liquidity and **Min suggested bet** judges the capped bet); an alt sitting on its main line's current number is hidden. Freshness caveat measured 2026-09-11: the anonymous changes stream carries NO alt updates and every alt's `modifiedOn` is the `0001-01-01` sentinel, so alts refresh only with the per-league snapshot and their age reads from `sequenceNumber` (the change time in epoch ms). An alt row click expands the grid row's Alts (`setExpanded`) and finds the cell by fiber props (points/side/book/event), failing loudly with the main cell outlined; a ticket captured on an alt cell carries `watch.altPoints` so the watcher re-reads the same rung. **Group by market** (on by default): `feed.groupEdges` collapses the list to one card per (game, period, bet type, side) — a +EV opinion is directional — whose best line is the highest Kelly stake (taxes longshots, so a -110 main at +5% beats a +944 rung at +6%), with `N books · M lines` behind an expander; alerts then key on the card and re-fire only when its best edge improves. **Bet history flags (issue #114, core; merged 2026-09-11)**: a local **bets service** (`unabated_ticket/bets_service/`, `run.sh`, `127.0.0.1:8094` loopback-only, no auth) is the only place that signs Kalshi requests — the private key never enters the extension — and turns the account's fills + unsettled positions into normalised bet records (one per `(ticker, side)`; positions are the truth for open size, fills the VWAP entry price; `sources/kalshi.py`, `sources/kalshi_ticker.py` a series map + cached public `/markets` and `/events` GETs, football suffixes carry the Eastern DATE only, MLB `HHMM` too; unknown series and futures fail closed as "unmatchable", shown never matched), UPSERTs them on native id into `bets_service/bets.duckdb::bets` (never pruned — future CLV) with a `source_runs` row per poll, and serves `GET /bets.json[?days=30]` (`{generatedAt, sources: {kalshi: {fetchedAt, ok, error, count}}, bets}`) + `/health`; a failed poll keeps the previous records (a dark source never blanks the list), a source on its first poll reads `ok:false, "no completed poll yet"`. The panel polls it every 30 s while visible — **never from the service worker** — into `chrome.storage.local` `betsService` (payload + records, open + settled ≤30 days, ~1 KB each) and `betsSettings` (`serviceUrl`), resolving team keys through `teams.js` — NOT a hand table (user decision 2026-09-11): every league snapshot's own team list (id, name, abbreviation) is registered at runtime and persisted as `teamsIndex`, keys are `<league>:<Unabated team id>` (the ids the board's lines carry), venue spellings resolve by exact normalised name, then a three-row `ALIASES` list, then a leading-code strip ("PIT Steelers"), then a UNIQUE word-boundary containment (the query-extends-team direction only for State/University/College, so "Southern Mississippi" never becomes "Southern") — and treating a venue's `ok` pull as authoritative for that venue's records. `bets.js` (pure, node-tested on `tests/fixtures/bets/kalshi_fixture.json`; `normalize_kalshi` in Python is held byte-equivalent by `test_parity.py`) matches a bet to a line on league + team pair (either order) or rotation + time (≤30 min with a start; the Eastern date ±1 day without), refuses to guess when two board events fit ("ambiguous game"), and ranks four tiers: `same_line` (market, period, side and number), `same_side` (different number), `opposite` (the other side — red), `same_game` (any other market/period); Kalshi NO on a team market = the other team **or a tie** (NFL/CFB/soccer, `approx: kalshi_no_side_includes_tie`). Surfaces: a **Ticket banner** (strongest first, 5 then "+N more"), **Edges badges + sized stakes** — a held bet NEVER hides a line, it changes the size of the next one (user decision 2026-09-11, replacing a hide toggle): `bets.exposureOf` sums dollars risked on the same direction (`held`, same_line + same_side) and the other side (`against`), `betsview.stakeAdvice` turns the row's Kelly stake into one wording everywhere (user choice 2026-09-11, `stakeAdviceWords`): "wagered $350 → target $600, bet $250" (against positions read "wagered $200 against …, bet $500 (net $300 on this side)", full size reads "bet $0"), rows and cards carry a `held $N` / `against $N` (red) / `game` badge plus a dim "you hold …" position line, a **by my exposure** sort, and alerts skip only lines already held at size, a **Bets tab** (per-venue freshness green <5 min / amber <60 / red, open bets, and an **unmatched** list with a reason per bet — nothing dropped silently), and a **header line** under the tabs ("bets: N open · kalshi 20 s · betonline — …"). Venue split: **#115 BetOnline (merged 2026-09-15)** polls the account's bet-history report through the bets service on the bet_logger Keycloak refresh token (`sources/betonline.py`; registered when the recon cookie file exists; refresh only within 60 s of expiry under a file lock shared with `bet_logger/scraper_betonline.py`, pending bets kept, league via `bet_logger/utils.parse_sport`, no game date → `approx: game_date_unknown` and the matcher windows `placedAt` −12 h/+14 d keyed on the rotation, with `tierByRotation` for team names `teams.js` cannot key); **#116 Novig, #117 ProphetX** each add a `Source` (`sources/__init__.py` protocol: `name`, `poll_sec`, `fetch() -> list[record]`, raise on failure, never partial) or a content script writing storage directly; until they land those venues read "no source configured". **Novig (issue #116)**: the official NBX API is a separate "Liquidity Provider" account ($30k minimum deposit; docs.novig.com/lp-onboarding) and cannot see retail-app bets, so the primary source is `bets_service/sources/novig.py` — a one-time `novig_auth connect` (Auth0 PKCE against the app's public client, the user logs in, refresh token saved to gitignored `novig_token.json`; its OWN chain, never the app's rotating localStorage token) then a 60 s poll of the app's Hasura `order`/`parlay` tables for the trader resolved from the JWT `sub`, through `curl_cffi`; `normalize_novig` is held byte-equivalent to `novig_bets.js` by `test_parity_novig.py`. The content-script mirror is the token-free fallback: `novig_page.js` (MAIN world, `document_start`, wraps `window.fetch` before the app bundle captures it) mirrors the three Portfolio queries the `app.novig.us` app fetches for itself (`ActivePortfolioOrders_Query` / `SettledPortfolioOrders_Query` / `ParlayPortfolioQuery` → `api.novig.us/v1/graphql`; zero requests of its own, no token read) and `novig_content.js` normalises them via `novig_bets.js` (outcome index 0 = home/Over, `price` a 0–1 probability, `isBid:false` LAYS the outcome = the other side at `1 - p`, `_1H` on MLB = `F5`, props/futures/unsupported leagues fail closed) into `chrome.storage.local.betsNovig`, which the panel merges (a complete read is authoritative for the venue; a stale one says "open app.novig.us … to refresh"); shapes are bundle-derived, not a live capture — a divergent blob lands in the unmatched list as "unreadable Novig order (…)". No overlay, no order placement. Load unpacked from `unabated_ticket/extension`.

**Conditional Kelly against held bets (issue #130, 2026-09-19).** The panel sized every line alone and adjusted by subtracting dollars (`add $X`, `still $X against`); dollars at different prices and numbers are not comparable and a moneyline held against a spread was ignored (four live cards, K $8,000: New Mexico @ Oklahoma bet $183 → $296, Lions @ Bills "at full size" → add $188, UConn add $116 → $121, Chargers ML bet $189 → add $0). Now `betsview.stakeAdvice` sizes the row given the open bets on its market (`condkelly.js`, `ladder.js`): K = bankroll × multiplier (full Kelly on the whole bankroll then a quarter is wrong — it says add $500 to a line already held at full size $667); the new bet's chance from Unabated's edge (decision 2026-09-09 stands), a held bet's from the median `bacr` at its half-point rung today; no fees added (Unabated's prices include them, user 2026-09-17). Same market and period is exact off the ladder; another period is the worst case the fairs allow (comonotonic pairing) because no correlation is estimated (user decision 2026-09-18: no per-sport history, use what Unabated provides) — **this reverses the 2026-09-12 "show it, don't size off it" decision** on the unmerged `feature/unabated-related-period-bets`, whose `related_same` / `related_opposite` tier names were taken and which is retired; **only the same direction is sized across periods** — the logic review (2026-09-19) showed the literal rule sizes a cross-period hedge as if both bets won together (hold 1H Under $300, new FG Over: $407 vs $667 alone vs $802 under the measured 0.69 link), so the other direction in another period is left out and named `other period · not sized` (user decision 2026-09-19); another market is never in the math (worst case gave $35 vs $188 on the Bills card and independent differs from ignored by $1; open in #129). Guards are left out and named, never guessed. Found while building: the changes stream ADDS lines for known events and files team totals under the game total's `bt3`, so `feed.js` marks snapshot lines `fromSnapshot` and only those feed a ladder. `bets.exposureOf`, `badgeText` / `badgeKind` and the `add` / `at_size` / `reverse` advice kinds are gone. **Hedge sizing is deferred** (user, 2026-09-19): a price with no edge of its own stays $0 and there is no `hedge only` label. Plan: `docs/2026-09-18-unabated-ticket-conditional-kelly-plan.md`.
