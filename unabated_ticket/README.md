# Unabated Ticket

Chrome extension (Manifest V3, plain JS, no build step). Four tabs in one
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
- **Bets** — your own open bets (Kalshi, Novig, BetOnline, BFA, Wagerzon,
  Polymarket US) read from a local service, matched against the
  board so the Ticket tab says when you already have this line, the other
  side of it, or a bet on the game, and the Edges list flags the same. Issue
  #114; plan in `docs/2026-09-11-issue-114-bet-history-plan.md`. See *Bets*.
- **Teasers** — the Buckeye 6-point 4-team teasers to bet right now, each at
  or under Buckeye's $200 limit, priced off Unabated's fair and sized with
  Kelly over the whole set, around the teasers already open at BFA (BFA is
  Buckeye). Plan in `docs/2026-09-27-unabated-ticket-teasers-plan.md`. See
  *Teasers*.

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
  markets and the panel refuses rather than guess — **only at the same price**; an edge is for
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
(`Football · FG · Moneyline/Spread/Total · 12 books · ≥2.5%`) and opens the
drawer in place; sort and minimum edge (labelled `min edge [2.5] %`, default 2.5%) stay out on the toolbar. Bankroll and
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
a probability) the panel also prints the **order** the stake means, under the
dollar figure: `1,127 contracts @ 23.2¢ · $261.48` (`kelly.contractOrder`).
The price is Unabated's exact number for the book, taken as the all-in cost of
one contract; the count is `floor(stake / price)`, never rounded up past
Kelly; the cost shown is what the contracts actually spend, so the leftover
(under one contract) is visible rather than hidden. **Kalshi's price on
Unabated already carries Kalshi's fee** (a 22¢ ask shows as 23.2¢ = +331,
because 22¢ + 7% × 0.22 × 0.78 = 23.2¢), so the edge is net of the fee, the
count needs no fee model, and the cost matches Kalshi's own Cost line
(contracts × ask + fee) to the cent; the number in Kalshi's limit-price box is
its own ask, not Unabated's figure. The count divides the number to act on, so
a top-up shows the top-up's contracts, and a stake under one contract reads
`under 1 contract @ 23.2¢` rather than a zero. Sportsbook lines show no
contract row.

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
  (`extension/ladder.js`). Its payoff is `toWin` / `-stake`
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
  moneyline cuts are exempt), no stake on the record, parlay legs (an open
  BFA teaser counts as its ticket, below), a Kalshi NO
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
- **Open teasers count too (2026-09-30, teasers plan section 14).** An open
  BFA teaser (the tickets the Teasers tab reads, `teaser.openTeasers`) with a
  leg on the row's game and market is held like a bet whose payout also
  needs its other legs: the leg cuts the game's rows at its half-point (a
  teaser push loses), and the other legs — other games, independent, each
  at Unabated's fair — decide whether the ticket is still alive. Under the
  log objective a ticket cannot be replaced by its average or by "the other
  legs win", so `condkelly.solveStake` enumerates the other games (each cut
  at every number a leg needs), pools the outcomes by which tickets they
  leave alive (legs shared across tickets are joint) and sums the growth
  over those pools; at most 2^16 other-game outcomes. A leg on another game
  counts as won once its game has started (the board's clock, else BFA's)
  — BFA closes a teaser within minutes of a losing leg's game ending (9/27:
  the five Seahawks tickets at 20:30 UTC, ~10 min after that game, the
  49ers leg still to play), so an open ticket's finished legs won; one in
  progress is assumed to. A leg still to play that is not priced — no board
  game yet (CFB still loading, a league that failed) or no fair at its
  number — is unknown, not won: its ticket is left out, `other leg not
  priced: …`. A decline the teasers alone cause (their other games past the
  budget or out of line, their risk past `K`) does not throw the straights
  out: they size the row and the teasers say why they do not. The leg
  itself follows a straight's rules: another market is `game · not sized`,
  the other direction in another period is left out, a leg with no rung on
  the row's ladder is `no fair at …`. A ticket's leg on the row's game but
  the other market (a total under a spread row) is treated as an
  independent game, like every spread–total pair here. Other parlays stay
  out as `parlay leg`. Cal's 8 teasers on 9/27 (16:46
  UTC board, K $5,000): Colts +7.5 −270 at Novig (+2.3%, three legs on Colts
  +7.5) $315.90 alone → `add $0` — the teasers already carry that side at
  about Kelly size at a better price (a leg inside +300 is −241); all 14
  Novig/Kalshi/ProphetX rows that said bet something on the same side as a
  teaser leg went to $0. The other side grows (the lift above): Chargers
  +6.5 +120 against four Bills −1 legs $137.48 → $425.65 (averaging the
  other legs would say $492, counting them as won $1,104); Commanders +6.5
  $381 → $743. On the row a `teasers $600` chip (green on the row's
  direction, red on the other; the tickets in its tooltip) sits beside
  `held` / `against`, the related bets list one line per leg — "3 teasers
  on Indianapolis Colts +7.5 · $600 · BFA · other legs win 40–44%" ("rides
  on this leg alone" when no other leg is still to play) — instead of each
  leg record, and the verb is `add` when a teaser in the math is on the
  row's direction (so the edge-move tag shows there). The same advice sizes
  alerts, the Min suggested bet filter and the Ticket, whose position line
  reads `teasers $600` / `teasers against $800`; the "by my exposure" sort
  counts teaser stakes.

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

On an exchange line the number never passes what is resting at the price
(the feed's liquidity, on the Edges row and on the Ticket when the feed
holds the line at the ticket's price): `add $17 · all $17 liq · $71.06
alone`. With nothing resting the Ticket reads "Nothing resting at this
price" over `$0`.

Before acting on an `add`, read the tag next to the edge (Edges tab → [Why
an edge grew](#why-an-edge-grew), issue #132): `fair moved to you` is the
sharps agreeing, `book moved away` is the book ahead of a fair that has not
answered yet, `fair moved against you` is the market turning against the
side while the edge still reads bigger. The number itself is unchanged by
the tag.

Settings (bankroll, Kelly multiplier) sit at the foot of the Ticket tab under
**Sizing**, reachable from any tab via the ⚙, and persist in
`chrome.storage.local`. Defaults 30000 and 0.25. The bets service URL is on
the Bets tab, under the venue strip it feeds.

Copy puts one line on the clipboard:
`Seattle Mariners -133 · 57.0¢ @ Novig | fair -139 · 58.2¢ | edge +1.89% | stake $188.55 | to win $141.77 | payout $330.32 | Texas Rangers @ Seattle Mariners · MLB`.
When held bets changed the number the stake reads `stake $188.32 (add $188.32, $270.05 alone)`.
On an exchange line the contract order follows the stake: `stake $261.69 | 1127 contracts @ 23.2¢ | to win ...`.

## Edges tab

### Data

The public league snapshots (no login needed). The panel fetches them only
while it is open and stops the moment it closes; nothing runs in the service
worker. The changes stream below was the second feed until **2026-09-27**,
when every path under `api-k.unabated.com/api/markets/changes` started
answering HTTP 410 ("Polling change subscriptions have been retired. Use SSE
or WebSocket subscriptions."); its replacement (`POST <data API>/subscriptions`
then `EventSource`) runs on the logged-in session, so the scanner no longer
polls anything but snapshots.

| Feed | URL | What it carries |
|---|---|---|
| Snapshot | `content.unabated.com/markets/v2/league/{id}/odds.json?t=<30 s bucket>` (29 team-sport leagues, `feed.LEAGUES`; ~18 MB gzip in total, CFB alone 9.7 MB, regenerated ~every 27 s). The query is a cache buster: CloudFront hands gzip clients the bare URL from an edge cache that was 6.5 h old on 2026-09-10 | every row's `sides[side][ms<book>]` line: `points, americanPrice, sourcePrice, sourceFormat, bacr, ge, liquidity, statusId, sequenceNumber`, plus its `alternateLines[]` (same fields per rung); `teams`; `marketSources` |
| Changes (retired 2026-09-27, HTTP 410) | `api-k.unabated.com/api/markets/changes/query[/{cursor}]` (~300 KB per 10 s) | the same fields per changed line under `gameOddsEvents[lg:pt:pregame][].gameOddsMarketSourcesLines[si:ms:an][bt]`, plus `sideKey` |

`ge` is Unabated's edge as a fraction (0.0296 = +2.96%), the same number the
Ticket tab sizes from. `bacr` is the fair at the book's points.

Loop (`extension/scanner.js`): snapshot per enabled league on open (four
downloads at a time, each league listed the moment it lands, CFB last), then
each league's snapshot re-downloads on its own cadence by compressed size
(≤2 MB every 60 s, ≤5 MB every 2 min, larger every 5 min; NFL is ~3 MB, CFB
9.7 MB), plus a full resync every 10 min. Expect roughly 6–8 MB/min with
every sport on. A 10 s tick checks which leagues are due; hiding the panel
pauses it, and a pause over 2 min resyncs on return. Before 2026-09-27 the
changes stream topped this up every 10 s, but even then it was incomplete
anonymously (2026-09-10, 3 min: 69 of 191 NFL line changes; Kalshi, Caesars,
ProphetX, Polymarket and Underdog mostly missing), which is why the
per-league cadence exists. Its parser was deleted after the 410; tests move
a line with `tests/helpers.js::moveLine`.

Parsing (`extension/feed.js`, node-tested on real slices under
`tests/fixtures/`):

- A line is keyed `(marketId, book, sideKey)`. Snapshot game rows are unique
  per event, period and bet type; a refresh replaces the league's lines whole.
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
  Main lines never carry either field; they refresh with the snapshot. Each event collects them in
  `event.venueIds`: `kalshiEventSuffixes` (e.g. `["26SEP19DUQWSU"]`, kept
  whole — Kalshi's team codes are not Unabated's abbreviations),
  `kalshiContracts` and `novigOutcomes` (id → `{lineKey, mainKey,
  points, sideIndex}` — `points` / `sideIndex` are the contract's own strike
  and Unabated side, fixed for the id; `lineKey` is the listed line at that number — the main line
  when a Kalshi or Novig rung sits on the main number (priced or not), else
  null for an unpriced rung; a join checks `points` against the line's
  current number). An id
  of any other shape is "no id", never an error. Coverage that day: 95 of
  319 NFL/CFB/WNBA events had Kalshi ids, 98 Novig ids. `describeLine` hands
  every row its event's map by reference (`row.venueIds`), which is how the
  bet matcher joins on them (Bets → What is matched).

### Alt lines

Off by default. Tick **Include alt lines** in the filter box and every
book's alternate spreads and totals join the list under the same gates as
main lines (board, book, bet type, period, edge, start, line age, liquidity).
The depth cap is in probability, not points: an alt lists only while
BOTH Unabated's fair and the book's own price put its side between **15% and
85%** (about +567 / −567; `feed.altWithinDepthCap`, a constant in code, no
setting), so it means the same depth in every sport. A rung whose Unabated
fair equals the next rung's on the same side is a flat-lined tail, not a
fair, and never lists (`feed.flatFairRungKeys`). It replaced the 7-point **Max pts from
main** setting on 2026-09-30. Inside the cap a deep rung is ranked down by
the tail flex (below), not hidden.

An alt sitting on the main line's current number is hidden (it would be
the same bet twice). On the grid, an alt cell's book is read by object
identity in the row's ladders, then the column, before the line's own
`marketSourceId` (Sports Interaction's alts carry BetMGM's id).
Rows carry an `alt` badge and say `alt of -2.5` (the book's main number);
the header counts alts apart (`3,120 lines (+40,278 alts)`).

**Freshness.** Alts are exactly as fresh as the league's last snapshot
refresh — 60 s / 2 min / 5 min by file size. Every alt's `modifiedOn` is the feed's `0001-01-01T00:00:00`
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
line**, which is the highest **tail-flex rank score** — EV dollars after
flex, below — so a -14.5 at +270 for 13.9% outranks a -27.5 at +1300 for
17.3%; then `▸ 2 books · 7 lines (+6)`, which opens the other books and rungs.
The card is its best line — it carries that line's own number, and a rung
behind the expander names its own when it differs. Every line is clickable
(locate) as before. The count badge counts cards,
sort orders cards through their best line, and `feed.groupEdges` (pure,
node-tested) does the grouping; the panel passes the rank score as the rank.
Turning the toggle off gives the flat list.

### Tail flex (rank score)

Unabated's fair is solid at the main number and flexes in the tails (the
host on Unabated Live, 2026-09-11, 49:26: the fuselage is solid, the wings
flex). The panel prices that flex to **rank** lines only — the edge and the
stake on every row stay Unabated's (`kelly.kellyStakeFromEdge` on `ge`,
conditional Kelly untouched). Per line (`extension/tailflex.js`, pure,
node-tested), with Unabated fair p = (1 + edge) / decimal:

- `dist_sd = |qnorm(p) − qnorm(fair of the same side at Unabated's main
  number)|` (Unabated's own line, book 49; even money when it has none; 0 on
  a moneyline);
- `sigma_z = sqrt(DZ_MAIN² + (c × dist_sd)²)`, `DZ_MAIN` = one cent at the
  median in probit;
- `sigma_e = decimal × dnorm(z) × sigma_z`, `keep = edge² / (edge² + sigma_e²)`
  (Baker & McHale 2013);
- **rank score = keep × edge × stake** — EV dollars after flex, the stake
  being the row's own Kelly stake. It picks a card's best line and orders
  the lines inside it. It is not shown: the row shows the edge and the bet.
  The "by stake" sort still sorts cards by their best line's stake.

**c is measured live** per (league, period, bet type) on every scanner
update from exchanges (`hasLiquidity` books) quoting BOTH sides of a rung at
the mirrored number, fresh (on the board, inside **Max line age**) with an
implied sum in [0.995, 1.20]: probit-devig the pair, take the probit gap to
Unabated's fair, remove the RMS gap the exchanges show at Unabated's main
number, divide by `dist_sd` for rungs ≥ 0.3 SD out, drop the top 1% (stale
quotes), take the RMS — c is a standard deviation, so not the median. Fewer
than 100 such rungs → `c = 10%`. No floor on a measured c (user, 2026-09-30). Exchanges are a measuring stick only; no
fair is blended. The Edges header shows the c in use for every spread/total
market on the list (`tail flex: NFL spr 7.6% · CFB spr 9.3% · MLB tot 10%`).
Measured c, NFL + CFB full game, 2026-09-30: NFL spr 7.6% / tot 5.6%, CFB
spr 9.3% / tot 9.2% (~100 ms per measure over 194k lines).

A line whose fair is Unabated's ±999900 placeholder (a clamp, not a fair)
never lists, never measures, and never prices a rung of the conditional-Kelly
ladder (`feed.isClampedFair`; a held bet there is left out and named).

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
**Min liq to win $** (default 100) hides an exchange line, main or alt,
when the money resting at its price would win less than that: $20 at
+2000 wins $400 and stays, $20 at +100 wins $20 and goes (Cal,
2026-09-23). It replaced a flat $100 stake floor that gated alts only and
hid longshots whose whole Kelly bet is small; a Novig Portland Fire +809
main moneyline with $17 behind it listed with "add $32.58" under it, and
now lists (it can win $137) as "add $17". Books with no figure pass, 0 =
off. Liquidity is read as dollars you can stake (not verified against
payout) and comes from the league snapshot, so it can lag a price move by up to a snapshot cycle. Sort by edge,
stake or start time. Settings (sports, periods, bets, books, minimum edge,
minimum suggested bet, max line age, min liq to win, sort) persist in `chrome.storage.local` under `edges` (sports
as league ids; a sport checkbox toggles all of its leagues).

Leagues come from `feed.LEAGUES` (ids probed 1–70 on 2026-09-10, labels
read off team names). Tennis (ATP 9, WTA 10) and combat (22) are not listed:
their sides key on people and their bet types are not 1/2/3. The row-click
screen path is verified for nfl/cfb/mlb only; nba, cbb, nhl, wnba and
soccer follow the site's nav and, for a soccer league, the odds screen must
have that league selected for the row to be found.

The header shows leagues loaded, lines held, when the newest snapshot was
built (its `Last-Modified`; minutes or more means a stale edge copy),
and the filter in effect. A league that fails to load is named in a red banner while the rest
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

### Why an edge grew

Next to the edge figure of a line you already hold in the same direction —
the rows, card lines and Ticket whose stake reads `add $X` (conditional
Kelly, Stake above) — one small tag says what moved since the panel last
saw the line, inside a ten-minute window (`extension/edgemove.js`, issue
#132), or since your first fill on the line once its fair was saved (*Since
your first fill* below). A line you do not hold, or hold only on the other side, shows no
tag (user decision, 2026-09-23): the tag exists for the adverse selection
of *adding* to a position, where the reasons an edge can grow call for
opposite actions, and a first bet is not a top-up. The history behind the
tag is still recorded for every line, so a line bet later is tagged at
once. **The fair decides**; the price only refines the reading.

| Unabated's fair (`bacr`, in probability) | The book's price on the side | Tag | Meaning |
|---|---|---|---|
| moved toward the side | anything | `fair moved to you` (green) | the sharps agree; the price is a bonus — closing-line value |
| unchanged | got better | `book moved away` (amber) | the book is shading against the side and the fair has not answered yet: Unabated's fair is ~1–2 min behind the book (v2 rebuild ~1/min + ~40 s ingest lag, #126), so the edge spikes for a refresh or two. Wait one snapshot; if the fair holds and the price is still there, it is a stale soft line |
| moved against the side | anything | `fair moved against you` (red) | the market is moving against the side and the book is ahead of the fair. The edge can still read *bigger* because the price improved by more than the fair fell — adding on it is adverse selection: on an exchange a great price resting right after a move usually belongs to someone who knows the line moved |
| unchanged | unchanged or worse | none | |

The rule, in probability (`kelly.americanToProb` on the fair; the exchange's
exact `sourcePrice` on the price, so a Kalshi cent or a Novig half-cent is
measured as itself): a move under 0.5 points is unchanged — `bacr` is a
whole American price, so a one-point fair change near even money (-110 →
-111, 0.23 pts) is rounding, and the fair must move about three American
points there to register (at +200 a five-point move does). The comparison is
against the newest observation at or before ten minutes ago, else the line's
first sighting; a line seen once, or unchanged for ten minutes, has no tag,
and a book that only shortened never earns one. **A line whose number moved
reads as first seen**, not as a move: a price at 48.5 is not comparable to
one at 47.5 (1,466 of 2,720 live NFL spread/total lines sat on a different
number than they opened, measured 2026-09-22).

Under the tag, one small line names the mover: `fair 33.7% → 35.6%` when
the fair decided (green or red), `price +125 → +141` when it was the book
(amber). The tooltip carries everything: `fair 33.7% → 35.6% · price +199 →
+215 · moved 2m ago · opened +185`. "ago" is when the panel first *saw*
the move: every observation is a snapshot, up to one refresh interval
(60 s / 2 min / 5 min by league) after the book actually moved. `opened` is the book's own opening price (`openerPrice`, on every main
line and no alt rung, per book — #126, measured 2026-09-22); when the book
has moved its number since open it reads `opened -120 at -3`, because the
price is not comparable across numbers. The opener's fair is not in the feed,
so it is context, never a tag.

The history lives in the scanner's memory only (`scanner.getHistory()`,
handed to the panel with every update): it starts empty when the panel
opens or the leagues change, so the first refresh shows no tags, and a line
the snapshot no longer lists is forgotten. Alerts carry the tag's words in
the notification body (Alerts below), so a re-fire on an improved edge says
which improvement it was. **The stake never changes for the tag** — it
informs; the number to act on stays what conditional Kelly computes. A hold
rule ("no add for one snapshot after an amber, never on a red") is a
separate decision once the tag has been watched on live cards.

**Since your first fill** (user decisions 2026-09-23). On a line you hold in
the same direction, the tag measures from your EARLIEST open straight bet on
that very line — same period, bet type, side and number, at any venue; a
parlay leg is not a position on the line — whose fill-time fair was saved,
instead of from ten minutes ago: the fair then
against the fair now, and — on the bet's own book only — that bet's fill
price against the price now (Kalshi's exact VWAP in cents, Novig's exact
probability, a sportsbook's American price), through the same four cases at
the same 0.5-point threshold (`edgemove.classifyMove`, the one rule both
paths call). On another book (a Kalshi row for a Novig bet, any row for a
BFA or Wagerzon bet) the price is not compared — the gap between two books
is not a book moving away — and the fair decides alone; the fair itself is
the same across books. A saved fair is used only on the market it was read
from (the `marketId` its line key starts with, one per side of a bet type,
shared by every book), so a bet the matcher later puts on another game
never carries it there. The
detail line names the fill — `since your +208 bet: fair 36.0% → 34.7%`, or
`since your +208 bet: price +208 → +223` on amber — and the tooltip carries
the since-fill numbers, the last ten minutes as the plain tag would read
them (an amber there is the fair lag itself: the book moved and Unabated's
fair, ~1–2 min behind, has not answered — wait a snapshot) and the opener.
Nothing moved since the fill is no tag, whatever the last ten minutes did.
A held line with no saved fair at its number — including bets held only at
another number or in another period — keeps the ten-minute tag above,
unchanged. In the row's Related bets block and the Ticket banner, each bet
with a saved fair adds `fair then 36.0%`.

How a fill's fair is saved (`extension/fillfair.js`): after every bets poll
and every scanner update, each open bet without one is read off the
scanner's history — the newest observation of its line at or before
`placedAt`. A Kalshi NO that
also wins on a tie is refused — no line's fair is its payoff — and a
placement up to a minute past the panel's clock (a venue clock running
ahead) waits a pass, while one further out is refused as a time-zone error.
The bet's own venue's line first (Kalshi, Novig, BetOnline and
ProphetX are on the board), else the book at that number that changed its
line most recently: Unabated's fair is the same at every book on a number
(measured 2026-09-23: NFL 5,796 of 5,796 multi-book rungs identical, MLB
1,394 of 1,394; CFB 95.3%, the rest mostly alt rungs untouched for hours
that still carry an older fair, median gap 0.04 points). A fair is saved
only when the history is a gap-free record of the fill: the scanner has
watched since before `placedAt` (`status.observingSince` — set when the
first load after the panel opens lands, cleared when the panel hides, set
again only once the catch-up after it reappears has landed, and restarted
after more than two minutes without a successful snapshot),
the line's own league has loaded without a gap since before it too
(`status.leagueObservingSince`: one league's snapshot can keep failing
while the others land, and its exchange lines move on the snapshot only —
a load more than one refresh interval plus two minutes after the previous
one starts that league's run over), and the line had been seen by then,
with a fair. A bet on a market the board never carries (a prop, a team
total, Kalshi's F5) is refused at once. The panel POSTs each row on its own
to the bets service, which keeps the first capture per bet for good
(`bet_fill_fairs`, Bets service below); one row per request so a row the
service refuses never takes a good one down with it, and the served rows
are merged into the held ones rather than swapped in (a saved row never
changes, so a poll that left before a save landed cannot drop it). **No backfill** (user decision): a
bet placed before the panel was watching — every bet open when this shipped,
one placed with the panel hidden or closed — gets no saved fair and keeps
the ten-minute tag. A new bet reaches the panel well inside the ten minutes
the history keeps: the service polls Kalshi, Novig and BFA's open bets every
60 s and BetOnline, BFA's history and Wagerzon every 300 s, and the panel
polls the service every 30 s. Record quirks that shape the reading: a Kalshi record is one per
contract, with the FIRST fill's time and the VWAP of every fill, so its
baseline is the fair at the first fill against the average price; a Novig
order's time is its placement, resting or not; and Unabated's Kalshi price
includes the taker fee (a 22¢ ask shows 23.2¢) while the fill is the bare
contract price, so on Kalshi the price leg reads up to ~1.75¢ worse than it
is — it can hide an amber, never raise one. The stake never changes for the
tag.

### Live edges (between quarters)

Unabated's live screen prices NFL lines off an **in-game fair** that exists
only at a break: a per-game "fair price set" (`ready` at the end of a
quarter, `expired` once play resumes) that reaches the screen over the
logged-in stream and never the free v2 file (0 of ~5,700 live lines carried
a fair on 2026-09-27, including at the end of Q3 of SNF while the screen
showed edges). The app enables it for NFL, straight markets, live mode only.
Read from Unabated's app bundle, 2026-09-27: the screen sets
`row.inGameFairPrice = {status, producedUtc, checkpointType, lines}` on a live
row and `line.edge = {edge: EV %, linePrice, simulatedPrice}` on every line
of it, alt rungs included, `simulatedPrice` being the fair American price at
that line's points.

So the panel does not price live lines itself: `page.js` reads what the open
live screen computed. Once a second it walks the grid's live rows (skipping
final games and Unabated's own column) and posts `live_edges` — the games on
screen with their fair status, and every on-board line of a game whose fair
is `ready` with its edge and fair — on change, else every 5 s. `content.js`
stores it per league as `liveEdges:<league>`, and `live.js` (pure,
node-tested) turns it into rows for a **Live** block above the pregame list
(the pregame list gets a *Pregame* label while the block shows):

- **At a break**: rows in the usual style with a `live` tag, the checkpoint
  (`end of Q1`, `halftime`) instead of the kickoff, and `fair −124 · changed
  12s ago` on the book line; the Edges tab shows `N live` next to its count.
- **In play**: one muted line, `in play · edges return at the next break`.
- **Tab not reading** (no post for 10 s from a visible tab, 75 s from a
  hidden one): an amber line, `bring the Unabated live tab to the front`,
  and stakes read `—`, because the numbers may be gone. After 2 min without
  a post the block disappears; a tab that closes says so on `pagehide` and
  its rows go at once.

A hidden tab is the normal case (the bet is placed in another tab), and
Chrome runs a hidden tab's timers about once a minute after five minutes in
the background. So `page.js` also scans on the grid's own `modelUpdated` /
`cellValueChanged` / `rowDataUpdated` events (at most every 250 ms), which
the app's stream keeps firing in a hidden tab: a break's edges post when
they land, not on the next throttled tick. Every post carries `visible` and
the tab's `instanceId`; a pregame tab of the same league never blanks a
live tab's game while that tab is still posting (`content.js`).

The Edges filters apply as they do to pregame (Books, Bets, Periods, min
edge, min suggested bet, min liquidity to win, alt lines and their distance);
there are no live-only settings. Stakes are the line's own quarter-Kelly off
the live edge: open bets on the game show on the row as `live · not sized`
but are never netted, since conditional Kelly needs a fair ladder the screen
does not expose (the fair set carries per-market probabilities; a
follow-up). Clicking a live row outlines the cell on the screen (`locate`
with `live: true`, so `page.js` looks only among live rows). Alerts use the
same threshold, one per line per break (the key carries the fair's
`producedUtc`), and nothing alerts while the tab is not reading.

A price clicked on a live row builds a normal ticket whose fair is
`edge.simulatedPrice` (live lines carry no `bacr`), labelled `Unabated live
fair, end of Q1`, and whose watcher stays on the live row (`watch.live`).
The screen shows one league at a time, so live edges exist only for the
league of an open live tab; nothing is fetched to get them.

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
book, body the edge, the why-it-grew tag when the line is one you hold in
the same direction and something moved (`fair moved to you` / `book moved
away` / `fair moved against you`, Why an edge grew above — `… since your
+208 bet` when it reads from a saved fill fair), stake, matchup
and time to start. Clicking the
notification runs the same jump-to-row path as a row click. No alerts
fire while the panel is closed. The alert log lives in `chrome.storage.local`
(`alertLog`, 24 h).

### Etiquette

While the panel is open: the per-league snapshot refreshes above (6–8
MB/min with all sports on, ~1 MB/min with just
football/baseball/basketball/hockey), the same files the page itself loads;
nothing runs when the panel is closed. Live edges add no requests at all:
they are read off the page the user already has open.

## Bets

Two kinds of source feed the flags:

| Venue | Source | How it refreshes |
|---|---|---|
| Kalshi | `bets_service/sources/kalshi.py` (local service, signed REST) | every 60 s while the service runs |
| Novig | `bets_service/sources/novig.py` (local service, the app's Portfolio REST feed on the account's own Auth0 refresh token, #116) | every 60 s while the service runs; no tab needed |
| BetOnline | — (#115) | shows "no source configured" |
| BFA (Betfastaction) | `bets_service/sources/bfa.py` (local service, the account's own Keycloak password login from `bet_logger/.env`; 2026-09-23) | open bets every 60 s, history every 300 s, while the service runs |
| Wagerzon (the C account) | `bets_service/sources/wagerzon.py` (local service, the site's form login from `bet_logger/.env`; 2026-09-23) | every 300 s while the service runs |
| Polymarket US (the CFTC app) | `bets_service/sources/polymarket_us.py` (local service, the account's own API key from `bet_logger/.env`, Ed25519-signed; 2026-09-23) | every 60 s while the service runs |
| Bet105 | the panel itself (`extension/bet105.js`, your own login at app.bet105.ag in this Chrome — Cloudflare challenges anything else) → `POST /bet105.json`, parsed by `bets_service/sources/bet105.py`; 2026-09-29 | every 5 min while the panel is open; open bets only |
| ProphetX | — (#117) | shows "no source configured" |

Start the service (next section), keep the panel open. Every 30 s while the
panel is visible it fetches `http://127.0.0.1:8094/bets.json` (never from
the service worker), resolves each record's teams through `teams.js`, dedupes
on the venue's native id against what it already holds, and keeps open bets
plus settled ones from the last 30 days in `chrome.storage.local`
(`betsService`, which also carries the service's team crosswalk and the
saved fill fairs; `betsSettings` holds the service URL).
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
line (so a -35.5 bet on a line that has since moved to -36.5 is
`same_side`), and the id map's `lineKey` / `mainKey` are never read — a
Novig lay's outcome id names the
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

**Unmatched open bets: flag, attach, learn (2026-09-23; plan
`docs/2026-09-23-unmatched-bet-attach-plan.md`).** An open bet the board
cannot place is left out when the panel sizes the next bet on its game, so
it is flagged, attached by hand and learned from:
- *The flag.* `bets.unmatchedReasons(bets, lines, now)` marks each unmatched
  open bet `attachable` (a game bet — not `unmatchable` — whose league is on
  the board; a bet with no league counts when any league is) and
  `needsGame`: attachable, its game not started by what the bet knows
  (the start of the board event it last matched, remembered by the panel
  from `bets.matchedStarts` in `betsService.knownStarts` — so a BetOnline
  bet with no date, or a Kalshi football bet with only one, reads as
  started once its finished game leaves the board, not red until the
  venue settles it; else `eventStart` in the future, else an `eventDate`
  today or later in Eastern time, else no date at all), not a parlay or
  teaser leg (a leg sizes nothing either way), the board listing at least
  one game of its league on the bet's own Eastern day (a bet with no date:
  its placed window — so an MLB bet on tomorrow's game, placed before
  tomorrow's slate is posted, stays grey until it is) — **whatever the reason** for the
  miss (the broad rule, Cal 2026-09-28: the old rule flagged only named
  causes, and a wrong Novig team id that made its game read "not posted
  yet" stayed grey). The reason explains the miss: `team not recognised`,
  `ambiguous game`, `no event on the board yet`, or **start time differs** (the
  bet's game on the board — its team pair, or, for a bet that names one
  team (BFA, Wagerzon and BetOnline spreads and moneylines), its rotation
  on a row where that team fits — not started yet, within 12 h of the
  bet's start at another time; 12 h is under the gap between two games of
  one MLB series, and a doubleheader's game 1 in progress is left out so a
  game-2 bet is never pointed at it). When the bet's rotation is on that
  board event, the board's clock decides whether the game has started — the
  bet's own is the one in doubt (BFA's read 7 h early on 2026-09-26, a
  start before the bet was placed); a match by team pair alone keeps the
  bet's own clock, so a finished doubleheader game 1 never flags against
  game 2. Bets needing a game turn the Bets tab label red
  with a red count, add "N not matched to a game" to the header line, and
  put a red banner at the top of the Bets tab. Futures and props, leagues
  off the scanner, days the board has not posted, parlay legs and games over but not yet settled
  never flag; they fold into **Not on the board**. Measured on the live
  board of 2026-09-28 (437 events, 39 open bets): 6 unmatched, all futures
  and props, 0 red. Started games keep their
  pregame rows on the board until they end (9 MLB events 0–6 h past their
  start were still listed, 2026-09-23), so a bet in play stays matched.
- *Needs a code fix* (2026-09-26). An open bet its source could not read —
  the venue's own league code named a game league and the selection did not
  parse, so the source marked it `raw.parseFailed` (BFA and Wagerzon, see
  Bets service) — is `needsFix`: red from the poll that serves it, board or
  no board, until the parser is fixed and the service restarted, or the bet
  settles. It has its own **Needs a code fix** block, banner and "N needs a
  code fix" in the header line, and no Attach: the record names no game,
  team or line. Deliberate exclusions — props, team totals, a league not
  supported, a postponed leg — carry no marker and stay grey. A record with
  no parsed pick reads the venue's own text (`raw.description`, markup tags
  as " · ") where it used to read its id.
- *Attach* (`attach.js`, pure, node-tested). Step 1 lists the board's games
  in the bet's league whose Eastern date is within a day of the bet's (no
  date: from the day it was placed to 14 days later), games where one of
  the bet's teams already resolves or its rotation matches first; two
  letters in the search box search every game of the league at any date.
  Step 2 lines each team the venue names up with the picked game's team on
  the same side (**Swap** flips it: a neutral site, or a one-team bet on the
  other side), tagged `known` (already resolves there, nothing written),
  `learn` (resolves to nothing) or `fix` (resolves to another team, shown
  as "was …"), and restates the bet on the game ("Your bet: Abilene
  Christian +7.5") so a wrong side shows before saving. **Attach and learn**
  POSTs `/pins.json`: the pin (bet → board event) and the learn / fix names
  as crosswalk rows marked with the pin, which replace a held key — a
  manual lesson outranks an id join, and automatic learning stays
  insert-only, so it never overwrites one.
- *Dismiss* (2026-09-28). Every red row — Needs a game or Needs a code
  fix — carries a **Dismiss** chip (beside Attach, or alone): for a game
  Unabated never lists, or a bet Cal has decided to live with. The bet
  leaves the red count, banner and header, and moves into Not on the board
  tagged **dismissed**, keeping its reason and Attach; **Restore** flags it
  again. It still does not size the next bet. Dismissals are panel view
  state — the bet ids in `chrome.storage.local` (`betsService.dismissed`),
  passed to `bets.unmatchedReasons(…, {dismissedIds})` — and
  `betsview.keepDismissedOpen` drops a dismissal once its bet settles,
  closes or leaves the store, so it never hides a later bet.
- *After.* `bets.applyPins` puts each served pin on its record (`pin`; a
  record with no league takes the pin's, marked `leagueFromPin`, and gives
  it back on Undo) and `resolveGame` reads the pin before the id join and
  the name rule while the pinned event is on the board. The open row reads
  **attached** with **Undo** (`DELETE /pins.json?betId=`: the pin and the
  names it taught go; a row it replaced is not restored), and the taught
  names appear under Team crosswalk as "attached by you". Automatic
  learning skips pinned bets. The panel keeps the pins with the records in
  `betsService.pins`.

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
  poll reads "no completed poll yet"; the service itself shows
  "unreachable since …" in red with the last records still listed), the
  open bets (venue, bet, stake, placed — each with a green left edge when
  the board matched it, red when it did not, so a problem bet is visible in
  the open list too; an attached bet reads **attached** with **Undo**; a
  line under the list says how many open game bets match, or "Every open
  game bet matches a game on the board").
  Unmatched open bets are listed with why in two places (above): **Needs a
  game**, right under the money with an **Attach** button on each row (a bet
  its source could not read sits under it in **Needs a code fix**, with the
  parser's reason and no Attach), and
  the folded **Not on the board** under the open list (Attach there too when
  the bet is a game bet). The reasons: team
  not recognised (the raw name, so `teams.js` can grow), ambiguous game,
  start time differs, no event on the board yet, league not on the scanner,
  not a game market (futures, the bots' combos), unknown Kalshi series. A bet with a venue id
  names both tiers that missed: "by id: Kalshi event 26SEP19DUQWSU not on
  any board ladder; by name: team not recognised (…)" ("Novig outcome" for
  Novig; "only on another league's board ladder" when the id sits on an
  event of a different league; no id tier for a Novig moneyline or a league
  off the scanner) — and, last, the **Team crosswalk** list with its Clear
  button (above).

### Novig source

Novig's official NBX API is a separate "Liquidity Provider" account with a
$30k minimum deposit, so it cannot see bets placed in the retail app. The
service instead logs in **as the retail account once** and reads the same
Portfolio feed the app reads (issue #116):

```bash
/Users/callancapitolo/NFLWork/kalshi_draft/venv/bin/python3 -m unabated_ticket.bets_service.sources.novig_auth connect
# or, to avoid the copy/paste race on the callback:
/Users/callancapitolo/NFLWork/mlb_sgp/venv/bin/python3 -m unabated_ticket.bets_service.sources.novig_auth connect --browser
```

- **Connect** runs Auth0's PKCE authorization-code flow against the app's
  public client (no secret exists). You log in yourself; the script only
  needs the URL you land on (`https://novig.com/?code=…&state=…` — copy it
  at once, the app strips it as it loads; `--browser` opens a Playwright
  window that intercepts the callback so nothing can be lost). The refresh
  token goes to `NOVIG_TOKEN_PATH` (`bets_service/novig_token.json`,
  gitignored, mode 0600) and the source registers on the next service start.
- **Why its own token.** The web app keeps a *rotating* refresh token in
  localStorage, and Auth0 revokes the whole chain when a rotated token is
  reused — copying the app's token would log the app out and kill the
  poller. A separate login is a separate chain; both live side by side.
- **Polling** (`sources/novig.py`, every 60 s): refresh the access token
  inside a 2-min margin (a rotated refresh token is rewritten at once), then
  `GET https://api.novig.us/nbx/v1/portfolio/{active,settled}?currency=CASH&sort=<tab>_recency&limit=50`
  with the app's `novig-client-capabilities: card_stack,player_props`
  header, following `nextCursor` to the end of both lists (a list that
  does not end within 200 pages, a page without `items`, or any non-200
  fails the poll loudly and keeps the previous records). Transport is
  `curl_cffi` with Chrome impersonation, the session the anonymous Novig
  SGP scraper already gets through Cloudflare with.
- **Why REST (2026-09-22).** Novig moved its web app from `app.novig.us`
  to `novig.com` and put its Hasura GraphQL behind a query allowlist the
  same hour — every hand-written query, and the old app's own Portfolio
  queries, now answer `query is not allowed`. The new app's Portfolio screen
  is a REST feed on the same Auth0 client and audience the service already
  mints for, so nothing in the path can tell the poller from the app. The
  earlier content-script mirror of `app.novig.us` (three extension files)
  went with the old app.
- **Renewal.** Auth0 chains usually cap at about 30 days; when the refresh
  fails the poll fails loudly (the Bets tab row goes red with the error and
  the previous records stay) and you run `connect` again.

`normalize_novig()` (pytest on `tests/fixtures/bets/novig_portfolio.json`,
live cards captured 2026-09-22 and validated against the 476 orders the
GraphQL-era normaliser had recorded — side agreed on every one) applies the
feed's own shapes:

- A **straight** card is one order; a **union** card is several orders on
  one market, each an order line under `stateBlocks[].{home,away}.lines[]`
  (active) or `settledItems[].line` (settled), sharing the card's
  `eventHeader` and `marketLabel`; a **parlay** card lists its legs under
  `sgpGroups[].legs[]` and `singleLegs[]`, priced as a whole. The card id
  of a straight or union card is the market id.
- The outcome on a line is the side **held** — a lay already reads as the
  other outcome, so there is no `isBid` to flip (the old lay flip had
  priced one order at the wrong side's 0.70; the feed's own 0.305 is what
  was paid). `price` is that side's 0–1 probability, `amounts.cost` the
  dollars risked, `payoutCopy` the gross payout (thousands commas), so
  contracts = cost / price and toWin = payout − cost.
- One order can appear twice: its matched part (`Matched` / `Win` /
  `Loss`) and a `Canceled` line for the unmatched remainder, under the same
  `orderId`. The matched line is the bet; a `Canceled` line is a bet only
  when the order never matched at all — then `void` with $0 at risk. A
  `Settled` state is a **fractional** settlement (`outcome.status` `"0.72"`,
  `"0.50"` — a first-five tie paid at half, a partial) and is graded on
  what it paid against its cost, toWin being the actual P&L. `Draw`/`Push`,
  `Cash Out` and `Unmatched`/`Pending` (open on the resting size with
  `approx: ["novig_order_unmatched"]`) map by name; any other state is
  `unknown` and logged, never silently open.
- `outcome.subtitle` names the market: `Moneyline` / `Spread` / `Total`
  are full game, `1st Half …` is `1H`, `First 5 …` is `F5` (MLB), anything
  else (`Passing Yards`, …) is "not a game market"; the number is the
  outcome title's trailing number (`ATL -1.5`, `Under 58.5`), else its
  display label (`O 7.5`), else the union's `marketLabel` (`7.5 Total`).
  The team is `competitor.id` against the event's teams (index 0 = home /
  Over second). Leagues map NFL→nfl, NCAAF→cfb, NBA/NCAAB/WNBA/MLB/NHL and
  the soccer leagues; the rest fail closed as "league not supported (…)".
  Parlays give one record per leg (`id` `novig:<parlay>:<leg>`, stake =
  the parlay's wager, no leg price).
- Every team object carries **`unabatedId`**, Unabated's own team id, which
  rides in `awayTeamVenue` / `homeTeamVenue` (with Novig's `id`, `name`,
  `shortName`, `symbol`). In `bets.js` it keys a side only when the team's
  name resolves nowhere and the id is a team of that league: Novig now sends
  it for most CFB teams too, and some carry an id that is no CFB team
  (2026-09-26: Jackson State 276 where the board's is 777; UTEP, Prairie
  View A&M, Eastern Washington too), so the name goes first.
  Novig's `symbol` is not Unabated's abbreviation (`UTC` vs `CHT`): stored,
  never a key. `venueIds: {marketId, outcomeId, eventId}` are the venue's
  own; the outcome id is what Unabated's Novig rungs carry as `sourceData`
  (see Edges → Parsing) and the matcher joins on it (Bets → What is
  matched). Novig sends no rotation number, so a Novig bet reaches a line
  by its outcome id or its team pair only.

## Teasers

Buckeye 6-point teasers on NFL and CFB: the set of 4-team tickets to bet right
now, each at or under Buckeye's $200 limit, sized with the panel's Kelly
settings. Unabated's own teaser calculator prices one hand-built ticket at a
typed-in price; this prices every Buckeye leg on the board and spreads the
Kelly amount over as many $200 tickets as it is worth. Plan:
`docs/2026-09-27-unabated-ticket-teasers-plan.md`. Code: `extension/teaser.js`
(pure, node-tested) and the Teasers tab in `panel.js`.

**Legs.** Every Buckeye (Unabated source 59) main full-game spread and total on
the board, game not started, line changed within the Edges tab's **Max line
age h** (a dead feed keeps old numbers). Teased 6 points: a spread's own number
+6, an Over −6, an Under +6. Buckeye's price is ignored: −120 and +100 tease
the same. **A push loses** (Buckeye's tie rule is unconfirmed), so a leg wins
exactly when the half-point one step against it wins: Bills −1 (from −7) is
priced at −1.5, Browns +8 at +7.5, Over 38 at 38.5. The win chance is
Unabated's fair at that half-point, read off the game's ladder (`ladder.js`)
built from that bet type's lines alone — the spread ladder leaves the
moneyline out, since Unabated's moneyline fair and its alt spreads disagree by
about 2.5 points near zero (SEA @ WAS 77.1% vs 79.5%). A leg with no rung at
its number, or a rung that repeats its neighbour's fair (Unabated flat-lines
deep tails), is listed grey with the reason and never used.

**Can't tease (college only, 2026-10-03).** Buckeye keeps some college games
off its teaser menu; every NFL game teases. A CFB leg in the pool or above
break-even carries a **Can't tease** chip: it takes that game's spread or
total, both sides (`teaser.marketKeyOf`), out of the pool, and the game's other
market can still step in. The list rebuilds when the market held a pool leg;
otherwise only the Legs list changes. The marked market sits in the
fold (`can't tease`, `spread blocked` / `total blocked`) with **Restore**;
restoring on the same fairs brings the same list back. Marks live in
`chrome.storage.local` (`teaserBlocked`: market key → the game's start) and
each ends when its game starts.

**Tickets.** The pool is each game's best leg — one leg per game, so legs are
independent — the top 10 by win chance (2^10 = 1,024 joint outcomes; 12 legs
added nothing). Every 4-leg combination of the pool (210) is a candidate. All
four legs win: +3 per dollar (+300); anything else loses the stake. A leg
breaks even at 70.7% (0.25^¼).

**The set.** One score decides everything: E[ln(1 + P&L / K)], K = bankroll ×
Kelly multiplier (the ⚙ settings; the objective conditional Kelly uses for
held bets). Starting from the open BFA teasers and the straight bets held on
the same games, the tab adds the $200 ticket
that raises the score most, again and again, until no $200 ticket raises it;
that ticket's best stake below $200 (whole dollars) is the last, partial one.
No ticket twice, and no per-leg cap (user decision 2026-09-28: "no leg on more
than half the tickets" cost 6% of the value after risk). Open teasers are part
of the set: an open 4-team ticket whose four legs are all live is never
offered again, and once an open 4-team ticket under $200 is in play (the
list's partial, placed — a 2-team teaser does not count) no second partial
is; the empty state says so, with the partial it held back. So placing the
list in order leaves the rest of it, and placing all of it leaves nothing.
The joint outcomes stay under 131,072: open teasers on many games outside the
pool, or straight bets splitting many games into more rows, cost the pool its
weakest free legs (and the straights on their games) rather than the list.
Every ticket wins exactly where its four legs do, so the search scores the
2^10 states of "which pool legs win", not every outcome: 9/27 with its 15
straights (69,984 outcomes) builds in ~0.2 s. Two legs with the same fair
(Colts +7.5 and Titans +8.5 are both 74.7% on 9/27) make tickets that tie
exactly; a tie goes to the first combination in the pool's order. Two ways
to see what it does, on Sunday's board (2026-09-27 16:37 UTC, $20,000 at
quarter Kelly, K = $5,000, before counting straights):

*Each leg lands near its target.* A leg inside a +300 4-teamer carries −241
(the 70.7% break-even); its target is quarter Kelly as a straight bet at
−241, and the set puts about that much on it:

| Leg | Win | Target | In the set |
|---|---|---|---|
| 49ers −1.5 | 77.1% | $1,089 | $892 |
| Seahawks −2.5 | 76.0% | $911 | $800 |
| Colts +7.5 | 74.7% | $684 | $600 |
| Titans +8.5 | 74.7% | $684 | $600 |
| Bills −1 | 74.0% | $571 | $492 |
| Browns +8 | 73.9% | $548 | $600 |
| Broncos +7.5 | 73.9% | $548 | $492 |
| Lions −1 | 72.8% | $354 | $400 |
| Patriots +9 | 72.3% | $277 | $292 |
| Steelers +9 | 70.7% | $0 | $0 |

*Every extra ticket on a leg is scored against the tickets already riding on
it.* Step by step, the ticket bet and the best one the other way on the 49ers,
with what each adds to the value after risk:

| Step | Bet | EV | Adds | Best the other way | EV | Adds |
|---|---|---|---|---|---|---|
| 1 | 49ers / Seahawks / Colts / Titans | 30.7% | $48.0 | Seahawks / Colts / Titans / Bills | 25.6% | $38.0 |
| 2 | 49ers / Seahawks / Bills / Browns | 28.2% | $34.8 | Seahawks / Bills / Browns / Broncos | 22.9% | $29.3 |
| 3 | 49ers / Colts / Broncos / Lions | 23.8% | $23.2 | Colts / Titans / Bills / Broncos | 22.0% | $19.8 |
| 4 | Seahawks / Titans / Bills / Broncos | 24.2% | $16.9 | 49ers / Seahawks / Titans / Broncos | 29.3% | $15.5 |
| 5 | 49ers / Titans / Browns / Patriots | 23.0% | $12.1 | Colts / Titans / Browns / Patriots | 19.2% | $8.8 |
| 6 | Seahawks / Colts / Browns / Lions | 22.1% | $3.5 | 49ers / Seahawks / Colts / Patriots | 26.5% | $1.2 |
| 7 | 49ers / Bills / Broncos / Patriots, $92 | 21.9% | $2.5 | Colts / Bills / Broncos / Patriots | 18.1% | $1.7 |

At step 4 the better ticket by EV (29.3%, with the 49ers) adds less than a
24.2% one without them: the 49ers already ride on three tickets. The set: 7
tickets, $1,292, +$324 expected (25%), worth $141 after risk, 47% to make
money, 22% that every ticket loses. Splitting finer does not cut the swing
(measured 9/27): 92–97% of the variance is how many legs cover. The Kelly
multiplier is the variance dial.

**Placed tickets come from BFA — nothing to mark.** BFA is Buckeye, so every
open BFA teaser the bets service reports, from this tab or not, is a placed
ticket. Its legs arrive one record per leg (`isParlayLeg`, `parlayId`, the
ticket's stake and to-win on each, `raw.headerDescription` "4 TEAM TEASERS",
the leg's teased number, rotation and start), are grouped by `parlayId` and
joined to their games through `bets.js`'s matcher (rotation + start within
30 min; BFA's rotations are Unabated's), else by rotation alone within 12 h
of BFA's start — a leg the matcher misses (BFA's clock off, as it once was
by 7 h; a team name resolving elsewhere) would otherwise count as won and
its ticket come back into the list; past 12 h it could be next week's game
with the same rotation. `bets.teaserLegPositionOf` places a leg on its
game's axis without the parlay-leg and stake guards straight bets have. A
placed ticket is fixed P&L in the score — +to win when every leg covers,
−stake otherwise, a push losing — until its last game starts. A started leg,
a leg no board game matches and a leg with no fair count as won: the ticket
keeps its full weight on the legs still to play, and no fair has to be saved.
A game with a placed leg offers new legs only on that leg's market (spread or
total); when Buckeye has moved the number since, the placed and the new leg
share the game's rows — the game is one outcome cut at both numbers, as in
conditional Kelly — rather than counting as two independent legs. The bets
service polls BFA's open bets every 60 s, so a ticket shows within about 90
s of being placed; the tab says how old the last BFA pull is.

**Straight bets on the same games count too (2026-09-30).** Every open
straight bet the bets service reads — any venue, spread, total or moneyline,
not a parlay leg — that sits on a pool leg's game and market, full game, is
fixed P&L in the score the way conditional Kelly sizes a held bet: its
numbers cut the game's rows (a whole number both sides of its push) and it
pays +to win, 0 on a push or −stake in each, priced off the same ladder as
the leg (a moneyline at the spread ladder's ±0.5 rungs). No correlation
number: where the two numbers sit on the game's fair ladder is the link.
`teaser.heldStraights` joins each bet to its game with `bets.js`'s matcher
(the Edges tab's: pins, venue ids, names, rotation). A straight on the leg's
own side means fewer tickets on that leg (its losses get deeper), one on the
other side means more. On 9/27 Cal held 15 straights on 9 of the 10 pool
games ($4,542): counted, the set goes from 7 tickets / $1,292 to 4 / $676,
and the set's value after risk on top of the straights from −$8 to +$62 —
the Lions and Patriots legs drop out, the 49ers ($892 → $476) and Seahawks
($800 → $400) halve, and Steelers +9 comes in ($0 → $200) with Bengals −3.5
held. Left out, as on the Edges tab: a bet on the game's other market (a
total under a spread leg), a 1H bet, a bet on a game no ticket uses, and one
with no fair at its number; a bet whose rung is out of line with its game's
ladder is left out rather than failing the build. A straight never adds a
game or changes the leg a game offers. If the open teasers and straights on
these games could already lose K, there are no tickets and the tab says why.
Their own P&L is in none of the summary's figures; the note under it says
how much straight money the set is sized around.

**The list holds still.** Many 4-leg combinations are near-tied: rebuilding
on every tick gave lists at 9:06 and 9:37 on 9/27 that shared 1 of 7 tickets
(no leg had moved more than half a point) and were worth the same (+$326 vs
+$324). So the list is built on each leg's fair as of the build, replaced
only when that leg moves 1 point or more, and it is rebuilt only when what it
is built on changes: the pool, a pool leg's number, an open BFA teaser, a
straight bet placed or closed on a pool game, K, or a 1-point move (a
straight's own number counts). Smaller ticks move the win chances, EVs and summary in place.
Built on the same fairs, the list after you place its first ticket is the
same list less that ticket. The reference fairs are kept in
`chrome.storage.local` (`teaserRefs`, rewritten on each rebuild), so a
reopened panel builds the list it last showed; and after a start or a
scanner restart (ticking a sport on the Edges tab) the tab waits for both
NFL and CFB (or their errors) before building — a list off NFL alone while
CFB, loaded last, is still coming would not be the list. A frozen list alone would miss news: the 49ers
leg down 5 points halves the 9:37 list's value after risk ($141 → $71) —
that move rebuilds it.

**The tab.** The count on the tab is the tickets to bet. A summary (tickets
and dollars to bet, expected profit, the chance to make money and the chance
every ticket loses — with open BFA teasers, what is placed and the expected
profit of all of them together), then **Open at BFA** (each open teaser
read-only, its legs at their current fairs or why they count as won, and the
age of the last BFA pull), the ticket cards in the order to bet them (legs,
stake, EV, **Copy**: `[rotation] side` per leg, BFA's own form) and the legs
(each game's best leg: teased number, Buckeye's line, win chance, line age,
how many tickets use it, and on a pool leg the Edges tab's tags for the
straights on its game — `held $X` on the leg's side, `against $Y` on the
other, a bare `game` when none of them counts, every bet and why in the
tooltip; the break-even line, which a pool leg can sit
below; the other games below break-even, on an open teaser's other market or
with no fair behind a fold, the markets marked can't tease first). An open ticket out of the math says why: every
game started, or no leg priced on the board. The Placed figure in the
summary is the open tickets still in the math; the Open at BFA line says
both when they differ. Empty states: the
board loading, fewer than 4 priced games, no ticket worth betting. The bets
service down puts a red banner up and the list stays (an open teaser may be
missing from it). The scanner always loads NFL and CFB for this tab; the
Edges tab and its alerts still list only the sports you tick
(`feed.selectEdges`'s `leagueIds`).

**Limits.** Unabated's fair at teaser numbers runs up to 3 points above the
exchanges (Kalshi, Novig, Polymarket, ProphetX); EVs and stakes are as good as
that fair (user decision to price off Unabated anyway, 2026-09-27). Push = loss
understates a ticket if Buckeye really drops a push to the 3-team price (9/27:
+$329 vs +$418 expected on the same board). Same-game spread + total
correlation is not modelled (one leg per game). Football, full game, 6 points,
4 teams only. A straight bet counts only on its leg's market and the full
game: Unabated gives no spread–total or full-game–1H link, so a total held
under a spread leg and a 1H bet size nothing. The other way round, the Edges
tab sizes a straight bet with the open teasers on its game (Stake, "Open
teasers count too"); its other parlays are still left out. A ticket placed
in the last ~90 s may not be in the tab yet.

## Bets service

A local Python service that turns the user's own bet history into normalised
records the panel can match against a line (issue #114; plan in
`docs/2026-09-11-issue-114-bet-history-plan.md`). It is the only place that
signs Kalshi requests — the private key never enters the extension. Read-only
GETs; no order placement.

It runs as a launchd agent (`bets_service/com.nflwork.bets-service.plist`):
it starts at login and restarts itself if it dies. A copy started from a
Claude session dies with the app (an app update killed one 2026-09-30).

```bash
# install once (or again after editing the plist)
cp unabated_ticket/bets_service/com.nflwork.bets-service.plist ~/Library/LaunchAgents/
launchctl bootstrap gui/$(id -u) ~/Library/LaunchAgents/com.nflwork.bets-service.plist
# restart onto new code (after a merge that touches the service)
launchctl kickstart -k gui/$(id -u)/com.nflwork.bets-service
# stop and remove
launchctl bootout gui/$(id -u) ~/Library/LaunchAgents/com.nflwork.bets-service.plist
```

- **Never start `run.sh` by hand while the agent is loaded**: the second copy
  fails on the `bets.duckdb` lock. The agent runs the main checkout's
  `run.sh`; a worktree test copy has its own `bets.duckdb` but needs another
  `BETS_SERVICE_PORT`.
- **Logs**: `bets_service.log` (rotating); crash tracebacks go to
  `bets_service.stderr.log`. After a hand kill, the first relaunch can fail
  with "Conflicting lock" while the old process lets go; launchd retries 30 s
  later.
- **Install**: nothing beyond the `kalshi_draft/venv` (duckdb, cryptography);
  `run.sh` uses it when present, else `python3`. Launch from the repo root.
- **Credentials**: `KALSHI_API_KEY_ID` + `KALSHI_PRIVATE_KEY_PATH`, read from
  the environment, then `unabated_ticket/bets_service/.env`, then the bots'
  `kalshi_draft/.env` in the main checkout, then `bet_logger/.env` there (the
  sheet scrapers' book logins: `BFA_USERNAME` / `BFA_PASSWORD`, `WAGERZONC_*`
  else `WAGERZON_*`, and the Polymarket US API key `POLYMARKET_US_KEY_ID` /
  `POLYMARKET_US_SECRET_KEY` from polymarket.us/developer) — so with the
  bots and the sheet scrapers configured no new file is needed. `.env.example` lists every knob (port, retention window,
  Kalshi cadence, log level). Never commit `.env`.
- **Endpoints** (loopback only, no auth): `GET /bets.json[?days=N]` →
  `{generatedAt, sources: {kalshi: {...}, betonline: {fetchedAt, ok, error, count}}, bets: [...], crosswalk: [...], pins: [...], fillFairs: [...]}`
  with open bets plus settled/closed ones within `N` days (default 30), the
  whole team crosswalk and every pin (newest first), and the saved fill
  fairs of those bets (oldest placement first); `GET /health` → `{ok, uptimeSec,
  sources}`; `POST /crosswalk.json` with `{rows: [{venue, league,
  venueTeamKey, unabatedTeamId, venueTeamName?, unabatedTeamName?,
  learnedFrom?}]}` (at most 1000 rows, 1 MiB) → `{ok, learned, conflicts,
  crosswalk}` — INSERT only, a held key with another id is returned in
  `conflicts` and never rewritten; `DELETE /crosswalk.json` → `{ok, cleared,
  crosswalk: []}`. **Pins** (the Bets tab's Attach, 2026-09-23): `POST
  /pins.json` with `{pin: {betId, venue, league, eventId, eventStart?,
  awayTeamId?, homeTeamId?, awayTeamName?, homeTeamName?}, crosswalk: [at
  most 2 rows of the pin's own venue and league]}` → `{ok, pins, crosswalk}`
  (404 when no stored bet has that id, 400 when the pin's venue is not the
  stored bet's or `eventStart` is not an ISO timestamp) saves the bet → board event pin and
  the team names the attach teaches, in one transaction; those rows REPLACE a
  held key (Cal is the authority) and carry `pinned_bet_id`, and a re-attach
  first drops the rows the earlier attach taught. `DELETE
  /pins.json?betId=` → `{ok, removedPin, removedRows, pins, crosswalk}` is
  Undo: the pin and exactly the rows it taught (a row it replaced is not
  restored). `POST /fill_fairs.json` with `{rows: [{betId, lineKey,
  points, fairAmerican, fairObservedAt, placedAt}]}` (at most 1000 rows) →
  `{ok, saved, fillFairs}` — the fair a new bet's line had when it was
  placed (*Since your first fill* above), INSERT only: a bet that has one
  keeps it, so `saved` counts only new bets. A row is refused with a 400
  that names it unless `fairAmerican` is a whole American price, both times
  are ISO and the fair was observed at or before `placedAt`. `POST
  /closing_fairs.json` with `{rows: [{betId, lineKey, points, fairAmerican,
  fairObservedAt, eventStart}]}` (at most 1000) → `{ok, saved}` — the server
  runner's pregame read of each open bet's line (Bet Tracker, CLV); a row
  observed at or after `eventStart` is a 400, and a row older than the stored
  one is ignored. `POST /exclusions.json` with `{betIds: [...], excluded:
  true|false}` (at most 50) → `{ok, changed, exclusions}` — the Bet
  Tracker's Remove / Restore; a bet not at BFA or Wagerzon is a 400, an
  unknown id a 404. `/bets.json` also serves `closingFairs` (the window's
  bets) and `exclusions` (all). `POST
  /bet105.json` with `{fetchedAt, feeds: {prematch: [betGroup], live:
  [betGroup]}}` (both feeds, at most 2000 groups) → `{ok, count, closed}`, or
  `{error}` → `{ok, recorded: "error"}` — the panel's read of Bet105 (the Bet105
  source bullet), stored as that source's run. **Settings** (the server
  runner, Running on a server below): `GET /settings.json` → `{settings:
  {bankroll, multiplier, leagues, periods, betTypes, bookMode, bookIds,
  minEdgePct, minStake, maxLineAgeHours, minLiquidityToWin, includeAlts,
  sortBy, groupByMarket}, updatedAt}`, every field `null`
  until set — null is the panel's default, read from `extension/edgerows.js`
  by the runner, never copied into Python; `PUT /settings.json` with
  `{settings: {<some of those>: value or null}}` → `{ok, settings,
  updatedAt}` sets the fields it names, null resets one, and a field it
  leaves out keeps its value. `bookMode` is `default` (the default books),
  `all` (every live book) or `custom`, and `bookIds` is a list exactly when
  it is `custom`. An unknown field or a value the panel's inputs would refuse
  (bankroll ≤ 0, multiplier outside (0, 1], no period, a bet type outside
  1–3, …) is a 400 naming it. The write routes require `Content-Type:
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
  default port). `BETS_EXTRA_ALLOWED_HOSTS` (comma-separated, default empty)
  adds names, each bare and with `:443`: on the server, `tailscale serve`
  forwards the browser's `Host` (the VM's tailnet name) to the loopback
  port. Entries must be bare lowercase host names — a scheme, port, path or
  wildcard stops the service at startup.
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
- **BFA source** (`sources/bfa.py`, 2026-09-23): logs in to
  Betfastaction's Keycloak with the password in `bet_logger/.env`
  (`BFA_USERNAME` / `BFA_PASSWORD`) and reads the **open bets every 60 s**
  (`BFA_POLL_SEC`; BFA is Buckeye, and a teaser placed there must reach the
  Teasers tab before it is suggested again — this plus the panel's 30 s poll
  is about 90 s) and the history every 300 s (`BFA_HISTORY_POLL_SEC`).
  Between history pulls a poll serves the open list plus the last pull's
  settled bets, so the venue's count holds still; the history's pending
  copies stay out, so a bet that has just left the open list keeps its
  open-list record (league, start) until the next pull settles it. The open
  bets are
  `GET api.bfagaming.com/history/api/GetPlayerOpenBets?playerId`, one object
  per wager with `betDetails[]` per leg carrying `idSport` (the league: `CBB`,
  `CFB`, `NFL`, …, so an open college bet IS placed), `gameDateTime` (the
  start) and the leg's description wrapped as `CBB - Alternative Lines <br>
  [1674] TOTAL u68½+110 (NEW MEXICO 1H vrs NEVADA 1H) [Sport:…, League:…]`
  (September's suffix is two brackets, `[Sport:Football][League:NCAA]`; both
  are stripped); its timestamps are naive on the account's **Pacific** clock,
  the same as the history's (checked live 2026-09-26: a bet stamped 12:42
  first showed on the 12:43 PT poll, and its `16:00` kickoff is the board's
  23:00 UTC start). The history is
  `GET …/GetPlayerHistory` for the last 31 days through tomorrow, which
  settles an open record once it leaves the open list (a store row stays open
  until a poll says otherwise); an open-bets record replaces the history's
  pending copy of the same wager. The session is this process's
  own and lives in memory: the access token is refreshed within 60 s of
  expiry and a refused refresh is replaced by a new login;
  `bet_logger/recon_bfa_auth.json` is never read or written (the weekly
  LaunchAgent rotates that token, and two rotators trip Keycloak's reuse
  detection — the reason BetOnline needs its file lock). Record ids are
  `bfa:<wager id>`, `:legN` per leg of a parlay or teaser (the ticket's
  `description` is leg 0, `picks[]` the rest, every leg on the ticket's
  status). Grammar from the live pull of 2026-09-22 (16 wagers) and
  `bet_logger/scraper_bfa.py`: `[rotation] TOTAL o30½-110 (RICE 1H vrs NOTRE DAME 1H)`
  (away first, `1H` on the names = 1H, `EV` = +100, a baseball total's
  pitchers bracket ignored); `[rotation] ILLINOIS ST -3-110` and
  `GRAMBLING 1H +168` (own team only, placed by rotation parity — the pull's
  totals follow the same convention, every over odd, every under even —
  `approx: side_from_rotation_parity`); a teaser leg's `(B+6)` dropped. BFA's
  rotations are Unabated's own (the pull's 308945 / 308959 and 461–481 were
  all on that week's board) except that a **first-half leg's rotation is the
  game's with a "1" prepended** (1340 for game 340, 1306551 for 306551), so a
  1H leg is served with that digit stripped (`raw.rotationAsWritten` keeps it)
  and a one-team 1H bet can match by rotation. The
  history description names no sport, so there the league is
  `bet_logger/utils.py parse_sport`'s nickname scan — pro leagues only — and a
  **settled college game is unmatchable as "league unknown"**, never guessed
  as CFB or CBB (hand-off rule); only open bets are flagged and those carry
  their `idSport`, so nothing that matters is lost.
  Timestamps are on the account's Pacific clock: a straight bet's
  `settledDate` was its scheduled kickoff on every settled row (noon-ET games
  read 09:00), so it is served as `eventStart` (`approx:
  event_start_from_settled_date`; refused outside placed −1 d / +60 d, where a
  .NET default date would parse), and `lastModification` — the grading time —
  as `closedAt`; a parlay or teaser has one `settledDate` for the ticket, so
  its legs are dateless (`game_date_unknown`, matched by rotation). Team
  totals and unparsed selections fail closed with the reason. A parlay or
  teaser the history cannot read — a leg count off its type (`4 TEAM TEASERS
  names 4 legs but carries 3`), a leg that does not parse, or a type it does
  not know that carries `picks` — still comes back one record per leg,
  `:leg0` up to the declared count, each unmatchable with the reason and on
  the ticket's status: the open list stored the ticket under those ids, and
  a lone `bfa:<id>` would leave them open forever (the store never closes a
  record a poll does not name). An OPEN leg whose `idSport` is a game code (`CBB
  CFB NFL NBA WNBA MLB NHL SOC`) but whose description does not parse is
  also marked `raw.parseFailed` (`normalize.mark_parse_failed`), which the
  panel flags as needing a code fix; a team total, a prop code and a code
  only the nickname scan resolves are never marked, nor is the history
  (it names no league). What an OPEN straight bet's `settledDate` holds is
  unobserved (the pull had none pending): a null or placeholder falls to the
  dateless window, so nothing is lost either way.
- **Wagerzon source** (`sources/wagerzon.py`, 2026-09-23, the C account): every
  300 s pulls `HistoryHelper.aspx?week=N` for the last 6 Mon–Sun weeks (the
  30-day window with margin) plus `OpenBetsHelper.aspx`, on an ASP.NET form
  login (`WAGERZONC_USERNAME` / `WAGERZONC_PASSWORD`, else `WAGERZON_*`, from
  `bet_logger/.env` — the C account's login has sat in the primary slot since
  2026-06-26) whose session cookie lives in memory; a helper answering HTML or
  a redirect instead of JSON means the session died and is replaced by one
  fresh login (a second miss fails the poll). Pending bets are kept: the
  helper is read first — `{result: [row, …]}`, one row PER LEG grouped by
  `TicketNumber`, each with `RiskAmount`, `WinAmount`, `PlacedDate` and
  `GameDateTime` (`M/D/YYYY h:mm:ss A`, Eastern), `IdSport`, `IdGame`,
  `DetailDescription`, `RotationNumbers`, `GameDescription` — the fields the
  site's own `OpenBetsTable` widget reads off the `ui_c` bundle; the values
  are unobserved (`{"result": []}` live, no open bet that night), and a row
  without `TicketNumber` and `DetailDescription` fails the poll naming its
  keys rather than being dropped. The history lists a pending bet too, under
  its game day with an empty `Result`, and settles it once it leaves the
  helper; an open-bets record replaces the history's copy of the same ticket.
  Record ids are `wagerzon:<TicketNumber>` (= `IdWager`), `:legN` per leg of a
  parlay (`details[]` is the leg list; each leg carries its own `IdSport`,
  `GameDate` + `GameTime` and `DetailResult`), transfer rows (`WagerOrTrans`
  `TRAN`) are not bets. The league is the leg's `IdSport` (`NFL CFB NBA CBK
  WNBA MLB NHL SOC`; `PROP` / `RBL` / `DST` are props and fail closed as
  "not a game market"), the start its `GameDate` + `GameTime` on the site's
  Eastern clock, so the matcher uses the 30-minute rule and a settled bet's
  `closedAt` is its event start. Grammar from every bet_logger run since
  2026-04: `[967] TOTAL o7½-120 (CHI CUBS vrs TB RAYS)<BR>( pitchers )` (away
  first, the pitchers bracket ignored), `[969] 1H ARI DBACKS -½+115<BR>( … )`
  (the period token BEFORE the name, no opponent anywhere — placed by rotation
  parity, `approx: side_from_rotation_parity`), `GM#1` / `GM#2` after a name on
  a doubleheader (kept in `raw.gameNumber`), `EV` = +100; MLB `1H` is the first
  five innings (`F5`, Novig's rule). A postponed leg ("( NYM vs COL Has Been
  Postponed. NO Action )"), an unsupported sport and an unparsed selection fail
  closed per leg with the reason; an unparsed selection on a game `IdSport`,
  open or in the history, is also marked `raw.parseFailed` (the panel's
  code-fix flag), a postponed leg or a prop never.
- **Polymarket US source** (`sources/polymarket_us.py`, 2026-09-23; the CFTC
  app at polymarket.us, not the international polymarket.com — Unabated lists
  them as two books, so the venue key is `polymarket_us`): every 60 s two
  signed GETs on `api.polymarket.us` — `/v1/portfolio/positions` (a slug →
  position map, `netPositionDecimal` positive = YES held, negative = NO) and
  `/v1/portfolio/activities` (newest first, paged until a page ends before the
  31-day window AND the fills in hand add up to the size of every open
  position and every settlement inside the window, so a bet placed weeks
  before its game keeps its full entry price and its settlement still closes
  the store's open row; a settlement older than the window adds nothing) —
  with the account's own API key (`X-PM-Access-Key`,
  `X-PM-Timestamp` in ms, `X-PM-Signature` = base64 Ed25519 over timestamp +
  `GET` + the path **without** its query string; a signed query is refused
  401). A record is one (market slug, contract side) of our own fills: each
  trade carries both orders and ours is `aggressor` when `isAggressor`, else
  `passive`; its `intent` (BUY_LONG / SELL_LONG = YES opened / reduced,
  BUY_SHORT / SELL_SHORT = NO) decides the side and buy/sell, so a "buy NO"
  booked against a held YES nets the YES; `price` is always the
  YES price, so a NO fill costs `1 − price`. Price is the VWAP of our buys on
  that side (fees excluded), the open size comes from the position, and a
  `positionResolution` settles it (`side` LONG = YES won, SHORT = NO won,
  NEUTRAL — a tie at $0.50 — is `unknown`); a side sold back to zero is
  `closed`, and fills still holding contracts with neither a position nor a
  settlement are `unknown` (a held bet is only ever read off the positions
  endpoint); an open position with no fill
  anywhere in the history is listed unmatchable and unpriced; busted and clearinghouse-rejected trades never count, and a trade
  state the docs do not list fails the poll. The trade's own `market` object
  gives the type — `sportsMarketType` `<sport>_(team|game)_<period>_<winner|spread|total>`
  for football, basketball, baseball and hockey, periods FG / 1H / 2H / 1Q–4Q
  / F5 — and each side's own number (`+2.50` / `-2.50`; the YES side can be
  the favourite). The teams and start come from one cached public GET per game
  (`gateway.polymarket.us/v1/events?slug=…&sportsMarketTypes=<the traded type>`,
  ~10 KB; the unfiltered event is up to 4.7 MB), away first (checked against
  Kalshi tickers and the sides' `ordering` field, which must agree or the
  record fails closed), spelled `safeName` for CFB, CBB and the NHL ("Liberty",
  "Ottawa Senators") and `name` elsewhere. Team totals, player props, soccer,
  unsupported leagues and combos (`caoc-…` parlays, listed as one record with
  their legs in `raw.comboLegs`) fail closed with the reason.
- **Bet105 source** (`sources/bet105.py` + `extension/bet105.js`, 2026-09-29):
  the one PUSHED venue. Bet105 (a LinePros white-label) sits behind Cloudflare,
  which challenges any request that is not the browser session's own, so the
  panel reads the account from your own login in this Chrome and the service
  only parses and stores. Every 5 min while the panel is visible: `GET /__bff/api/customers`
  (401 = not logged in; its `csrfToken` goes in `X-Broker-CSRF`), then the
  site's own `getHistory` (`{a: "getHistory", state: "0"}` = open bets) on both
  LinePros feeds, `/__bff/__partner-prematch/betLobbyV2/logic/` and
  `…/__partner-live/…`, POSTed as `{fetchedAt, feeds: {prematch, live}}` to
  `POST /bet105.json`; a read that fails is POSTed as `{error}` so the Bets tab
  row turns red with the fix (log in at app.bet105.ag in this Chrome). Never a
  half push: both feeds or an error. Record ids are `bet105:<feed>:<betGroupId>`
  (`:legN` per leg of a parlay, on the ticket's stake and price). A leg names
  the selection by id, not text (`description` is empty): `marketId` 5 total /
  6 spread / 3 moneyline (the LinePros wager types; 7 / 8 team totals, 1 = 1X2
  and props fail closed), `key` the line — the total, or the AWAY spread number,
  negated for the home side — and `subKey` the side, `1` = Over or the away
  team, `2` = Under or home, the same ids `bet105_odds/scraper.py` reads off the
  odds feed (`team1` / `team2` are away / home there too); `periodId` `m` = FG,
  `h1` / `h2`, `f5`, `q1`–`q4`; `leagueName` (`NFL`, the pro leagues; a college
  name fails closed until seen); `eventStartTime` epoch seconds; `finalOdds`
  decimal. Verified on the 2026-09-29 capture: three NFL first-half totals,
  side `1`. **Unobserved**, pinned by hand-written fixture rows only: an
  Under, a spread's sign, a moneyline, a parlay, the live feed. Bet105 shows no
  settled list (My Plays only ever asks for state 0), so an open record a
  complete push no longer carries is marked `closed` with no result — never
  guessed won or lost — and the store keeps it. No credentials and no poll in
  the service; `service.PUSHED_SOURCES` lists the venue so the panel reads "no
  completed poll yet" until the first push.
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
  `unabated_team_id`, `unabated_team_name`, `learned_from`, `learned_at`,
  `pinned_bet_id`) holds what the panel learned, INSERT-only from id joins
  and cleared on DELETE; an attach's rows replace a held key and name their
  pin. `bet_pins` (primary key `bet_id`, plus `venue`, `league`, `event_id`,
  `event_start`, both team ids and names, `pinned_at`) holds every attach,
  never pruned (a few a week).
  `bet_fill_fairs` (primary key `bet_id`, plus `line_key`, `points`,
  `fair_american` — Unabated's `bacr`, a whole American price —
  `fair_observed_at`, `placed_at`, `captured_at`) holds one fill-time fair
  per bet, INSERT-only with `ON CONFLICT DO NOTHING` (the first capture wins;
  a request repeating a bet keeps its first row and logs it), never pruned
  (~14 bets a day) and never backfilled; `/bets.json` serves the rows of the
  bets in its window only, so the table's growth never reaches the poll.
  `bet_closing_fairs` (primary key `bet_id`, plus `line_key`, `points`,
  `fair_american`, `fair_observed_at`, `event_start`, `updated_at`) holds one
  closing fair per bet, UPSERTed only forward in time (`DO UPDATE ... WHERE
  excluded.fair_observed_at >` the stored one), never pruned, served for the
  window's bets. `bet_exclusions` (primary key `bet_id`, plus `venue`,
  `excluded_at`) holds the bets removed from the Bet Tracker; Restore deletes
  the row, the bet in `bets` is never touched.
  `edge_settings` holds at most one row (`settings_id` = 1, checked) of the
  settings above, one explicit column each (`bankroll`, `kelly_multiplier`,
  `league_ids` / `period_type_ids` / `bet_type_ids` / `book_ids` as
  `INTEGER[]`, `book_mode`, `min_edge_pct`, `min_stake`,
  `max_line_age_hours`, `min_liquidity_to_win`, `include_alts`,
  `sort_by`, `group_by_market`, `updated_at`), UPSERTed
  whole by each PUT; NULL = the default. A failed poll writes a failed
  `source_runs` row and touches nothing else, so a dark source keeps serving
  its previous records; a store write that raises (disk full) is logged and
  retried next poll, never killing the poll thread. Until a source's first
  poll completes (Kalshi: ~1–2 min, one throttled GET per market and event)
  `/bets.json` lists it as `{ok: false, error: "no completed poll yet"}`.
  Log: `bets_service.log` (rotating, 10 MB × 3).
- **Adding a venue** (#117 ProphetX; BetOnline, Novig, BFA, Wagerzon and
  Polymarket US are `sources/betonline.py` / `sources/novig.py` /
  `sources/bfa.py` / `sources/wagerzon.py` / `sources/polymarket_us.py`
  above): a module in
  `bets_service/sources/` with `name`, `poll_sec` and `fetch() -> list[record]`
  (the `Source` protocol in `sources/__init__.py`), registered in
  `service.main()`. `fetch()` returns every record the venue knows and raises
  on failure — never a partial list. A venue no script can reach (Cloudflare
  challenging anything but the browser's own session — Bet105) is read by the
  panel instead and PUSHED to a route of its own (`POST /bet105.json`), the
  parser still a module in `sources/` and the venue listed in
  `service.PUSHED_SOURCES`. Records follow the contract in the plan
  (`id` = `"<venue>:<native id>"`, `side`/`points` in the side's own number,
  raw team names, keys `null`). A venue with no game date sets `eventStart`
  and `eventDate` null, carries `rotation` and `approx: ["game_date_unknown"]`,
  and the matcher windows on `placedAt` (see Bets). A content-script venue
  instead writes `{bets<Venue>: {bets, readAt, url, error, complete}}` to
  `chrome.storage.local` and the panel merges it with
  `betsview.mergePageSource`.
- **Capturing a venue's raw rows** (`bets_service/recon_capture.py`): a
  parser is built and re-verified on a live pull, and the cloud sessions that
  write parsers cannot reach the books, so this standalone script (only
  `requests`) runs on the Mac and its output is attached to the project
  thread. Today it covers BFA (Keycloak password login from `bet_logger/.env`,
  every `GetPlayerHistory` wager of the last 30 days through tomorrow, raw,
  with counts per result and type) and the Wagerzon C account (form login;
  `HistoryHelper` weeks 0 and 1; the `OpenBets.aspx` page plus every
  open/pending helper it or its scripts name, fetched raw). `cd` to the repo
  root and run `bet_logger/venv/bin/python3 unabated_ticket/bets_service/recon_capture.py`;
  files land in `~/Downloads/bets_recon/` (`--out DIR`). It writes no token,
  cookie or password, and leaves `bet_logger/recon_bfa_auth.json` alone.

## Running on a server (Oracle)

Work in progress toward using the panel from a phone with the Mac closed
(plan agreed 2026-09-30): the bets service and the Edges scan move to an
always-on Oracle Cloud VM, reached privately over Tailscale, with the Mac as
the fallback for any book that refuses a data-center login. **Step 0** is the
login check, **step 1** the headless Edges runner, **step 2** the phone
page, **step 4** the deploy (Docker + `tailscale serve`; all below); the Mac
relay (3) is still to come. Plan: `phone_page_plan.md` in the project files.

### Step 0: does each book accept a login from the server?

1. **Use the existing VM if there is one.** The MLB cloud-migration work
   (local branch `worktree-cloud-migration-mlb`, `deploy/cloud/ORACLE_SETUP.md`)
   created an A1.Flex VM `mlb-stack` with the key `~/.ssh/oracle_mlb.key`. It
   takes the whole Always Free A1 allowance (4 OCPU / 24 GB), so a second free
   VM is not possible; run this on that one. Only if it no longer exists:
   Oracle Cloud > Compute > Instances > Create, Ubuntu, shape
   `VM.Standard.A1.Flex`, a US region, your SSH public key, only SSH open.
2. **Get the code on it.** The repo is private, so add a read-only deploy key:
   on the VM `ssh-keygen -t ed25519 -f ~/.ssh/nflwork -N ""`, paste
   `~/.ssh/nflwork.pub` into GitHub > repo Settings > Deploy keys (read
   access only), then
   `GIT_SSH_COMMAND="ssh -i ~/.ssh/nflwork" git clone git@github.com:callancapitolo17/NFLWork.git`.
3. **Install.**
   `sudo apt install -y python3-venv && cd NFLWork && python3 -m venv venv && venv/bin/pip install -r unabated_ticket/bets_service/requirements.txt`
4. **Credentials: API keys only.** Never `scp` `bet_logger/.env`: it holds
   the BFA and Wagerzon passwords, which stay off the server until you decide
   to move those books (a data-center login may flag the account). Copy the
   Kalshi key ID and the `.pem` its `KALSHI_PRIVATE_KEY_PATH` names (fix that
   path on the VM), and put only the Polymarket US key lines in
   `unabated_ticket/bets_service/.env` on the VM. Then `chmod 600` each.
   Result on the VM, 2026-10-01: Kalshi, Polymarket US and the Unabated feed
   `ok`; Novig awaiting its own login; BFA, Wagerzon and BetOnline untested.
   - **Novig: don't copy the token.** Auth0 revokes the whole chain when a
     rotated refresh token is reused, so the VM needs its own login:
     `venv/bin/python -m unabated_ticket.bets_service.sources.novig_auth connect`
     (the paste flow; open the printed URL on any device, log in, paste back
     the URL you land on).
   - **BetOnline (only once its risk is accepted): move the cookie file, don't copy it.** Every refresh
     rotates its token and Keycloak kills a reused one. Stop the Mac's bets
     service and BetOnline scraper, `scp` `bet_logger/recon_betonline_cookies.json`
     over, run the check, then `scp` it back and restart the Mac side.
5. **Run it** from the repo root:
   `venv/bin/python -m unabated_ticket.bets_service.check_sources`

It prints one line per venue: `ok` with a record count, `not configured`, or
`FAILED` with the error, and a last line for Unabated's public NFL feed
(what the Edges scan needs). Bet contents are never printed. Exit code 1 when
anything failed. A book that fails here with a block or a 403 stays on the
Mac; one that is `ok` can move.

**Result on the VM (2026-10-01):** Kalshi ok (221 records; 341 s on a first
read, its per-market lookups are paced 0.6 s apart), Polymarket US ok, the
Unabated feed ok. Novig awaits its own login; BFA, Wagerzon and BetOnline were
not tested (see step 4). **2026-10-05:** Cal decided those three books don't
mind a data-center login, so they move to the VM (`deploy/README.md` step 10). The VM has no `python3-venv`; running the check in a
`python:3.12-slim` container with the repo mounted works without it.

### Step 1: the Edges list, headless (`server/runner.js`)

A Node process runs the panel's own Edges scan with no browser and serves
the list as JSON for the phone page. It `require`s the extension's pure
modules unchanged — `scanner.js` (the same loop and cadence: the settings'
leagues plus NFL and CFB, which the panel always loads for its Teasers tab,
each league again every 60 s / 2 min / 5 min by snapshot size, a full
resync every 10 min, never paused), `edgerows.js` (the list's selection,
conditional-Kelly sizing against your open bets and open BFA teasers, the
stake capped at liquidity, the Min liq to win / min edge / max line age /
books / alts toggle, the sort, the market cards with their best line picked
by the tail-flex rank, the stake rail's words, the held / against / teasers
badges and the edge-move tag), `tailflex.js` (c measured off the same
snapshot), `teaser.js` (the open teasers joined to the board, as the panel's
`openTeasersNow`), `betsview.js`, `bets.js`, `teams.js`, `ladder.js`,
`condkelly.js`, `kelly.js`, `edgemove.js`, `fillfair.js` — so `/edges.json`
is what the panel's Edges tab shows for the same settings, bets and
snapshot. Pregame only: the Live block reads the logged-in Unabated
screen, which a server does not have.

```bash
# Node 18+ (Ubuntu 22.04's apt node is v12 — install 20 or 22 from NodeSource or nvm). No npm install.
node unabated_ticket/server/runner.js       # http://127.0.0.1:8095/edges.json
```

| Env | Default | |
|---|---|---|
| `UNABATED_RUNNER_HOST` | `127.0.0.1` | bind address; stays loopback on the server too (step 4) |
| `UNABATED_RUNNER_PORT` | `8095` | |
| `BETS_SERVICE_URL` | `http://127.0.0.1:8094` | where `/bets.json` and `/settings.json` come from |

- **Reads**: Unabated's public league snapshots (no login, ~6–8 MB/min with
  every league on, as the panel); the bets service's `/bets.json` every 30 s
  (the panel's cadence), applied exactly as the panel applies it
  (`edgeRows.applyBetsPayload`); `/settings.json` every 10 s. A league change
  in the settings restarts the scan, as the panel's checkboxes do (NFL and
  CFB stay loaded either way).
- **Settings**: one row in `bets.duckdb::edge_settings`, edited with `PUT
  /settings.json` on the bets service (Bets service above), e.g.
  `curl -X PUT -H 'Content-Type: application/json' -d '{"settings":{"bankroll":8000,"minEdgePct":2}}' http://127.0.0.1:8094/settings.json`.
  A field never set is the panel's default (`edgerows.DEFAULT_STAKE_SETTINGS`
  / `DEFAULT_EDGE_SETTINGS`: $30,000 at 0.25 Kelly, 2.5% edge, 168 h, $100 to
  win, the default books, no alts, cards). The panel's own settings stay in
  `chrome.storage` — the two are not synced. `bookMode: "all"` is every live
  book (the panel's "My Unabated selection" has no Unabated tab to follow
  here). While the service cannot be read the runner keeps the last settings
  it read (the defaults before the first) and says so.
- **Serves**: `GET /edges.json` → `{generatedAt, settings {source, error,
  okAt, updatedAt, stake, edges}, scanner {phase, error, leaguesLoaded,
  leagueErrors, lineCount, altLineCount, snapshotBuiltAt, staleLeagues, …},
  betsService {okAt, error, unreachableSince, openBets, sources}, books
  {mode, ids, names, liveCount}, tailFlex, grouped, unit, total, maxShown,
  items}` — the top 200 cards (`{key, sideName, …, bookCount, lineCount,
  best, others}`, best the highest tail-flex rank, others in rank order) or
  lines; `tailFlex` is the panel's header line ("tail flex: NFL spr 7.1% ·
  …"). Each line carries the feed fields the panel renders plus `edgeTier`,
  `stake` (standalone, capped at liquidity), `rankScore`, `advice` (`kind`,
  `bet`, `alone`, `verb`, `held`, `against`, `teasers {held, against}`,
  `cappedAt`), `rail`
  (`{text: "add $237.25", note: "$437.25 alone", atSize}`), `badges`,
  `related` (the bets on the game with their tags and `fair then`) and
  `move` (`{kind, label, detail, sinceFill}` or null). The full shape is
  documented at the top of `server/edges_payload.js`. `GET /health` →
  `{ok, uptimeSec, scanner, betsService, settings, closingFairs}`. Any verb but GET is
  405; a request whose `Host` is not `127.0.0.1:<port>`, `localhost:<port>`
  or the bound address with the port is 403 (the bets service's #125 rule).
- **Closing fairs** (2026-10-05): it also loads the leagues of open straight
  bets (an Edges list never shows a league its settings leave out; a league
  added this way stays loaded until the runner restarts, since every change
  of leagues restarts the scan; a soccer bet loads every soccer league, as
  the bet names no competition) and every
  30 s POSTs `/closing_fairs.json` with the rows that changed — each open
  bet's Unabated fair on its own line while its game is still to start
  (`server/closefair.js`; Bet Tracker, CLV).
- **Side effects**: none on disk; its one write to the bets service is that
  POST. Feed, line history and bets live in memory and rebuild on restart (so
  the edge-move tag starts empty, as when the panel opens). Logs state changes
  to stdout/stderr. Unlike the panel it does not POST fill fairs or crosswalk
  lessons — the Mac panel keeps doing that.

### Step 2: the phone page (`server/phone/`)

A phone-sized, read-only page (Edges, Bets and Settings tabs and a ticket
sheet) served by the bets service, so the page, its reads and its one write
share one origin. You bet in the book's own app; nothing on the page places
anything.

```bash
venv/bin/python -m unabated_ticket.bets_service.service   # or bets_service/run.sh  (:8094)
node unabated_ticket/server/runner.js                     # must be running too    (:8095)
# then open http://127.0.0.1:8094/
```

- **Routes on the bets service**: `GET /` serves `server/phone/index.html`;
  its CSS/JS and the extension's pure modules (under `/ext/`) come from a
  fixed allowlist of paths (`service.STATIC_FILES`; anything else, traversal
  included, is a 404), `Cache-Control: no-store`, a same-origin CSP.
  `GET /edges.json` proxies the runner (`UNABATED_RUNNER_URL`, default
  `http://127.0.0.1:8095`; `UNABATED_RUNNER_TIMEOUT_SEC`, default 5) and
  answers 502 `{error, runnerUrl}` when it is down, slow or not 200. The
  Host allowlist covers every route, so the page opens only on the machine
  itself unless `BETS_EXTRA_ALLOWED_HOSTS` lists the tailnet name (step 4).
- **Edges**: the runner's cards as the panel draws them: badges, side, market,
  game and kickoff, book price with cents, line age, edge in its tier colour,
  the stake rail and its note, the move tag, related bets, the tail-flex line.
  Tap a card for the ticket: price, fair, edge, the number to act on with the
  panel's label, position line, contracts (exchange lines), to win / payout,
  every related bet and the card's other lines. Freshness and errors (edges or
  bets service unreachable since …, scanner errors, settings fallback) sit at
  the top.
- **Bets**: money at risk, Needs a game / Needs a code fix / Not on the board
  (the runner now adds `betsService.unmatched`, `bets.unmatchedReasons` over
  its board, to `/edges.json`), venue freshness, open bets. Attach and
  Dismiss stay on the desktop panel.
- **Settings**: every `/settings.json` field, its default shown when unset;
  Save sends only what changed (`PUT`, JSON), a 400's text is shown, and each
  set field has Reset (sends null). The runner picks a change up within 10 s.
- **Polling**: `/edges.json` every 15 s and `/bets.json` every 30 s, only while
  the page is visible; at once on becoming visible and after a save.
- **Reuse**: `phoneview.js` (pure, tested in `tests/phoneview.test.js`) formats
  on top of the extension's `kelly.js`, `feed.js`, `teams.js`, `bets.js`,
  `ladder.js`, `condkelly.js`, `betsview.js`, `edgemove.js`, `fillfair.js`,
  `tailflex.js` and `edgerows.js`, loaded unchanged as plain scripts.

### Bet Tracker (`server/tracker/`)

Daily P&L, CLV and bet analysis, served by the bets service at `/tracker`
(on the VM: `https://<vm>.<tailnet>.ts.net/tracker`). It reads `GET
/bets.json?days=3650` (every bet the service has stored, plus the saved fill
fairs, closing fairs and removed bets) every 60 s while visible. Its one write
is the Bets page's Remove / Restore.

- **Overview**: net P&L with the record, ROI, handle, open risk; expected
  P&L at Unabated's fair and actual vs expected with its z-score; CLV;
  cumulative
  actual vs expected chart with daily bars; a 6-week calendar heatmap; daily
  results; P&L by venue; open bets with price, fair, edge and CLV.
- **Analysis**: filter by venue, league and type (straight, parlay, teaser),
  group by venue, league, market, period, type, odds, edge at fill, timing
  (hours placed before start), weekday or stake. Each group shows ROI with its
  95% interval (proven only when the interval clears zero), expected ROI,
  actual vs expected, z, CLV and how often the close was beaten. Calibration of fair win probability against the
  actual win rate (straights with a fair), and a searchable bet log with each
  bet's closing fair and CLV.
- **Bets**: every BFA and Wagerzon bet, newest first, searchable, filtered to
  counted or removed. **Remove** marks a bet that isn't yours (`POST
  /exclusions.json` → `bets.duckdb::bet_exclusions`; the bet itself stays in
  `bets`) and it leaves every number on Overview and Analysis — a parlay or
  teaser goes as a whole; **Restore** counts it again. Only those two books
  (Cal, 2026-10-05); the service refuses any other venue's bet.
- **CLV** (2026-10-05): the expected return at Unabated's **closing** fair,
  closing prob × the ticket's actual payout − 1, stake-weighted, plus the
  share of bets that beat the close. The close is read by the server runner
  (`server/closefair.js`): every 30 s it POSTs each open straight bet's fair
  on its own line and number (the bet's venue first, else the book that moved
  last, an alt rung when the market moved off the number) while the game is
  still to start, and the service keeps the newest reading per bet
  (`bets.duckdb::bet_closing_fairs`). A reading more than 30 min before the
  start (the runner was down) is not a close, and a league whose snapshot
  build is over 15 min old (a stale CDN copy) gives none. History starts the day the
  runner first ran: no backfill. No close for parlays, teasers, props,
  futures, Kalshi combos, leagues Unabated does not list, or a number the
  board stops carrying. **The runner must be running at kickoff** — on the
  Mac, its launchd agent:
  ```bash
  cp unabated_ticket/server/com.nflwork.edges-runner.plist ~/Library/LaunchAgents/
  launchctl bootstrap gui/$(id -u) ~/Library/LaunchAgents/com.nflwork.edges-runner.plist
  launchctl kickstart -k gui/$(id -u)/com.nflwork.edges-runner      # restart after a pull
  curl -s http://127.0.0.1:8095/health                              # closingFairs {okAt, error, rowsSent}
  ```
  It needs `node` 18+ on Homebrew's PATH (`/opt/homebrew/bin`); logs go to
  `unabated_ticket/server/edges_runner.log`.
- **Header**: range (7D, 30D, 90D, YTD, All), $ / units with the unit size
  (default $100; range, units and unit size are remembered in the browser),
  and how many venues' last poll succeeded.
- **Rules** (`trackerstats.js`, pure, tested in `tests/trackerstats.test.js`):
  P&L lands on the **Pacific** day of the record's `closedAt`, which each
  venue fills differently: Kalshi the market's expiration, Polymarket US its
  resolution, Novig the ticket's settle time, BFA its grade time, Wagerzon the
  game's start, and BetOnline the time it was **placed** (its report has no
  settle time), so a BetOnline bet lands on the day you placed it. Kalshi's
  multivariate combos (`KXMVECROSSCATEGORY`, the MLB bots' RFQ fills on the
  same account) are their own type, "Kalshi combo", and count in every total
  by default (Cal, 2026-10-05); switching off the header's **Bot combos**
  toggle hides them. Won pays
  `toWin`, lost costs `stake`, push and void are 0. Open bets are exposure,
  not P&L; a Kalshi position sold before settlement (`closed`) and a bet whose
  result is gone (`unknown`, e.g. Bet105 once it leaves the open list) have no
  known P&L and are counted as "without a result". A parlay or teaser is one
  ticket (its legs carry the ticket's stake and status). Edge = fair
  probability × the ticket's actual payout − 1, the fair being the fill fair
  when one was saved, else the closing fair; expected P&L, z and calibration
  cover only bets with one, and the rest group under "No fair saved". A
  removed bet is in none of the numbers. Kalshi fees are not in Kalshi's `toWin`, so its P&L is
  before fees.
- **History** is whatever this service's `bets.duckdb` holds: sources pull
  about 31 days back, and records are kept from then on. A service started on
  a fresh database (the VM) starts its history there.
- **CSP**: the page loads only its own files by absolute path (it is served at
  both `/tracker` and `/tracker/`) and sets per-element styles through the
  CSSOM, since the service's CSP blocks inline style attributes.

#### On your phone, served from the Mac (current setup, 2026-10-05)

The Mac is the primary host for now: it has every book's login, the full
`bets.duckdb` and the Bet105 extension (the VM cannot pass BetMGM/ProphetX
geo checks or mint Cloudflare cookies). `tailscale serve` puts the Mac's
bets service on your tailnet, the same way the VM does, so the phone opens
the tracker (and the Edges page at `/`) from anywhere while the Mac is awake.

1. Install Tailscale on the Mac (Mac App Store or tailscale.com) and log in
   with the account your phone uses. In the admin console's DNS page,
   MagicDNS and **HTTPS Certificates** must be on.
2. Serve the bets service (the CLI ships inside the app):
   ```bash
   TS=/Applications/Tailscale.app/Contents/MacOS/Tailscale
   $TS serve --bg 8094          # https://<mac>.<tailnet>.ts.net/ -> 127.0.0.1:8094, survives restarts
   $TS serve status             # prints the https:// name
   ```
3. Let the service accept that name (the #125 Host guard 403s any other):
   add `BETS_EXTRA_ALLOWED_HOSTS=<mac>.<tailnet>.ts.net` (lowercase, no
   `https://`, no port) to `unabated_ticket/bets_service/.env`, then
   `launchctl kickstart -k gui/$(id -u)/com.nflwork.bets-service`.
4. On the phone (Tailscale app on, same account) open
   `https://<mac>.<tailnet>.ts.net/tracker`.
5. Keep the Mac reachable: on the charger, System Settings > Battery >
   Options > "Prevent automatic sleeping when the display is off" (or
   `sudo pmset -c sleep 0`). A closed lid still sleeps a laptop unless an
   external display is attached.

Nothing listens publicly: the service stays bound to 127.0.0.1 and only
devices on your tailnet reach the name. `$TS serve reset` takes it off.

### Step 4: deploy (`deploy/`)

Docker Compose on the VM runs both processes, each still bound to
127.0.0.1 (`network_mode: host`, nothing published); `tailscale serve` on
the host proxies `https://<vm>.<tailnet>.ts.net` to `127.0.0.1:8094`, so
only devices on your tailnet reach the page. Setup, update, logs, stop and
the "nothing public listens" check: `deploy/README.md`.

- **Bets service** (`python:3.12-slim`, wheels only) and **runner**
  (`node:22-slim`), `restart: unless-stopped`, json-file logs 10 MB × 3.
- **Repo mounted read-only** at `/app`, so an update is `deploy/deploy.sh`
  (pull `main`, rebuild, recreate, then PASS/FAIL checks of `/health`,
  `/` and `/edges.json` with the tailnet `Host`, and a 403 for a foreign one).
- **Writable data dir** (`UNABATED_DATA_DIR` → `/data`): `bets.duckdb`
  (`BETS_DB_PATH`), `novig_token.json` (rotates), `logs/`
  (`BETS_SERVICE_LOG_PATH`; `BETS_SERVICE_LOG_CONSOLE=1` also logs to
  stderr for `docker compose logs`). Off the server all three keep their
  defaults next to the package.
- **Secrets are files, never in an image**: the Kalshi `.pem` mounted
  read-only, `bets_service/.env` (Kalshi key id, Polymarket US keys) inside
  the read-only repo mount, with the BFA and Wagerzon logins since
  2026-10-05; BetOnline's rotating token file is moved (never copied) into
  the data dir (`deploy/README.md` step 10).
- **`BETS_EXTRA_ALLOWED_HOSTS`** in `deploy/.env` is the tailnet name the
  #125 Host guard must accept (Bets service above).

## Tests

One command runs everything and exits non-zero if any part fails:

```bash
./unabated_ticket/check.sh
```

It runs, in order, ESLint over `extension/`, `server/` and `tests/` (`npm run lint`),
the node suite (`npm test` = `node --test tests/*.test.js`, 428 tests) and
the bets service's pytest suite (388 tests, on the `kalshi_draft/venv`
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

`edgerows.test.js` pins the Edges list's rows on the real NFL slice an hour
before kickoff — the panel's defaults, the book modes, quarter Kelly with
nothing held and each row's tail-flex rank score, the gates, the sorts, the
cards (best line by rank), the tail-flex header, a held Bears -2.5 that
turns the same line into `add $237.25` (`$437.25 alone`) and the Panthers
side into `bet $328.86`, and an open $400 BFA teaser on the Panthers +3.5
that puts the Panthers lines at size — and the edge-move tag only on a held
line. `runner.test.js` drives `server/runner.js` with the slice behind an
injected fetch and a scripted bets service: `/edges.json` equals what
`edgerows.js` lists for the same state and settings (with and without the
teaser), a flat list and the settings row's minimum edge and sort, the
defaults and named errors while the service is down, NFL and CFB always
scanned, a league change restarting the scan, and the HTTP Host / verb /
path guards. `bets_service/tests/test_edge_settings.py`
covers `edge_settings` and `GET/PUT /settings.json` (validation, partial
updates, the Content-Type and Host guards).

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
exemption, no rung / whole number / other period, soccer moneylines feeding no cut, and the period
name → id map.

`edgemove.test.js` pins the why-it-grew tag (#132): fair up with the price
flat, better or worse is green; price better with the fair flat is amber;
fair against the side with the price better (the edge reading bigger) or
flat is red; amber turns red once the fair follows the book, with "ago"
counted from the book's move; a book that only shortened, a sub-threshold
one-point fair change (-150 → -151) and a first sighting are no tag; a
whole-point change (-150 → -160) and a Novig half-cent are moves; a number
move resets; a move older than the window expires; a line first seen inside
the window compares to its first sighting; a missing fair decides on the
price; pruning keeps one baseline; bad input fails loudly; `classifyMove`
reads the same four cases off any two observations and is exactly what
`edgeMove` returns over its window. `scanner.test.js`
adds the history's plumbing on the fixtures: one snapshot observation per
line on start and a clean slate on restart; a re-downloaded snapshot
records a fair move and an alt rung's price move and forgets the lines it
no longer lists. Its
`observingSince` cases: set when the first load lands, cleared by
`pause()`, set again when the catch-up after `resume()` lands (not when
`resume()` ran — every request takes 5 s on the test clock), a pass that
lands while hidden starting no run, and a 110 s outage keeping the run
while a 140 s one restarts it; `leagueObservingSince` keeps CFB's run
through one missed refresh and starts it over after five minutes of failed
loads while NFL's run and the global one go on; `leagueLoadedAt` records when
each league's last snapshot landed (a failed load keeps it, `start()` clears
it) — the Teasers tab redoes its board pass only when NFL's or CFB's changes.

`fillfair.test.js` pins the since-your-first-fill capture and reading on
the 49ers card that motivated it (Seahawks @ 49ers, a 49ers -14.5 alt held
on Novig at +208): the fair the bet's own line showed at the fill (+178,
not the +188 it shows now); the own venue first, and for a venue off the
board (BFA) the most recently changed book at the number; an own-venue line
first seen after the fill falling back to a book watched through it; every
refusal with its reason (before the panel was watching, first seen after,
no fair, a Kalshi NO that also wins on a tie, a prop or Kalshi F5 the board
never lists, placed in the future, no placed time, and a fill inside a gap
in its own league's snapshots); a bet with no line yet, and one a venue
clock 30 s ahead placed "after now", left for the next pass; the other side of the number never
read; saved, settled and unmatchable bets skipped. On the display side: a
reply that lacks a held row keeps it, for the bets still held; the earliest
open `same_line` straight bet with a saved fair is the baseline (a
same-side bet at another number is not, nor a parlay leg however early, nor
a fair read off another market); the card reads red — the price improved to
+223 while the fair fell 36.0% → 34.7% — amber with the fair holding, green
with it rising, and exactly `edgemove.classifyMove`'s answer; on a Kalshi
row the Novig fill's price is not compared (+230 with the fair flat is no
tag) while the fair still decides; a record with no price decides on the
fair alone; and the fill price basis per venue (Kalshi's VWAP cents,
Novig's probability, a sportsbook's American price).

`bets.test.js` joins synthetic Kalshi and Novig records
with unrecognisable team names on `fixtures/v2_venue_ids_slice.json`'s ids
(#118 step 3): the Kalshi spread YES and NO and the moneyline through the
suffix, the Novig bid and its lay on the outcome id, a main line moved off
the id's number (`same_side`), a suffix on two events
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
`bets.test.js` then matches Novig records as the service emits them
against the NFL slice: a moneyline on Carolina (same_line / opposite with
the Novig label, keyed on the venue's `unabatedId`s), `resolveTeamKeys`
keying on `unabatedId` only where the name resolves nowhere (never over a
name that resolves, never an id that is no team of the league) with a
learned crosswalk row still winning and `rekeyRecords` keeping it, a
resting order and a parlay leg
flagging the game, settled and void records never matching, and the
source's unmatched reasons passing through. For attach (2026-09-23),
`bets.test.js` covers a pin matching a bet whose name resolves nowhere, the
pinned event leaving the board, Undo, a pin lending its league, and the flag
rule (a name not recognised before and after kickoff, futures and leagues
off the board, date-only and dateless bets, ambiguous vs not posted, start
time differs vs the next day's game); `attach.test.js` covers step 1's
window, ranking, search and cap, step 2's known / learn / fix with Swap, a
one-team bet by name and by rotation, a doubleheader that teaches nothing,
the pin body and the labels; the pytest suite covers a pin replacing a held
row, Undo removing only what it taught, a re-attach, rollback, the request
validator, and the `/pins.json` round trip with its 404 / 400 / 415 and the
Host rule. On the Python side
`test_normalize_novig.py` runs the normaliser on
`fixtures/bets/novig_portfolio.json` (live cards, 2026-09-22): a matched
straight spread (the side held at its own price, cost / payout / contracts,
`unabatedId` as a string, cardId = market id), a team without an
`unabatedId`, union lines as one record each with the number from the
market label, an order seen matched and as a Canceled remainder kept once,
never-matched cancels void at $0, the feed pricing the side held (no lay
flip), settled wins and losses with `closedAt`, `1st Half` → `1H` and
`First 5` → `F5`, a fractional `Settled` line graded on what it paid,
props / unsupported leagues / unreadable cards failing closed with a
reason, parlay legs sharing the wager with no leg price, unknown states
logged as `unknown`, resting states flagged, and the line / leg walkers;
`test_novig_source.py` covers the auth (refresh once then cache to the
margin, rotation persisted with 0600, missing/broken token file, callback
state/error parsing) and the source (both Portfolio lists cursor-paged to
the end with the bearer on every GET, the URL grammar, and a page without
`items` / an unended list / an HTTP error / an auth failure all raising).

`teaser.test.js` builds synthetic boards (fair rungs whose `bacr` lands on
round win chances) and pins the Teasers rules: the teased number and cut of
an away / home spread, an Over and an Under; a whole number read at the
half-point one step against the bet (the rung the other way carries another
fair, so reading it would show); the moneyline left out of the spread
ladder; no rung and a flat rung grey with their reasons; started games,
lines off the board, stale lines and other books never teased; grading (the
one ticket's EV, expected profit, chance to make money and chance every
ticket loses); one ticket under the cap sized at the single-bet Kelly stake
in whole dollars; the greedy never above $200 a ticket nor taking a ticket
twice; the pool as each game's best leg, top 10; fewer than 4 priced games;
open BFA teaser legs grouping into tickets and joining their games by
rotation (a BFA parlay, a settled teaser and another venue's leg ignored);
an open ticket lowering the next suggestion on the legs it shares ($200 → $84
with $200 open on three of its legs); an open ticket's own four legs never
offered again, so placing the top ticket leaves the list less it (a board
with four standout legs, where it used to come back at the top); the whole
list open leaving nothing (one partial ticket per set); open teasers on
eight games outside the pool shedding one pool leg to stay at 2^17
outcomes rather than going dark; a leg whose BFA start is 3 h off the board's
joining by rotation so its ticket stays taken (13 h off: counted as won); a
5-team teaser with a leg counted as won blocking no 4-leg combination; the
partial rule ($63 alone, $45 with a 2-team $110 teaser open, none and $32
held back with a 4-team $150 ticket open); a reopened panel's saved
reference fairs rebuilding the list it showed where a 0.08-point tick
changes the list on live fairs; a started leg and a leg no board
game matches counted as won, and a ticket with none still to play out of the
math; a game with an open leg offering only that market; an open leg at
another number sharing its game's rows with the new leg; the list holding
still on a half-point tick (the EV moving in place) and rebuilt on a 1-point
move, on a new K and on a pool leg's number (not on +2.5 → +3, which a
push-loses leg wins on alike); and the Legs list's standings. Straight bets
held on the same games: `heldStraights` placing an away spread, a home whole
number (two cuts), an Over and a moneyline (the +0.5 rung) and naming a 1H
bet, no fair and no stake (parlay legs, settled bets and started games
dropped); a same-side straight lowering the dollars on its leg ($284 alone)
and an other-side one raising them; another market, another period and a
game outside the pool changing nothing (the row says why); a whole number
splitting its game at both sides of the push; straights that can lose K
giving no tickets with the reason; a rung out of line with its game's
ladder leaving that straight out while the list stands; the list rebuilt on
a straight placed on a pool game and on its 1-point move, not on one off the
pool nor a sub-point tick; the summary leaving the straights' P&L out and
counting them for the note, with the row's held / against dollars; and the
outcome budget shedding the two weakest legs with their games' straights.
`condkelly.test.js` adds `resultAt` (win, push, loss either side).
`bets.test.js` adds `teaserLegPositionOf` (the leg at its own teased number,
no stake guard, the parlay-leg guard kept for straight sizing),
`feed.test.js` the `leagueIds` filter, and `test_bfa.py` the cadence: open
bets every poll, the history every 300 s (not at 299), the open list plus the
last pull's settled bets in between, and a failed history pull retried on
the next poll.

Open teasers in the Edges stake (2026-09-30): `condkelly.test.js` checks a
ticket with no other leg sizing like a straight at its half-point, one other
leg and two tickets on one other game (joint) against two on separate games
(independent) each matching the log growth written out by hand and topped
cent by cent (and not the averaged ticket), Cal's 9/27 Colts +7.5 ($315.90 →
$0) and Chargers +6.5 ($137.48 → $425.65) on that board's fairs, and the
declines (17 other games past the 2^16 budget, an other game's ladder rising
with the cut) and loud failures (a whole-number ticket leg, a missing
other-game rung). `betsview.test.js` runs the same two rows through
`stakeAdvice` with `teaser.openTeasers`-shaped tickets: the `teasers $600` /
`teasers $800` chips with their tickets in the tooltip, the leg records
dropped from the matches, one related line per leg sorted among the
straights by tier and dollars, `add` with only teasers held, "rides on this
leg alone" with its other legs started, another market / the other
direction in another period / no rung each left out and named, an other
leg still to play with no board game or no fair leaving its ticket out
(and counting as won once BFA's clock says its game began), a decline the
teasers alone cause (six tickets on 18 other games) leaving the straight
sizing the row (add $0), a declined calc zeroing the teaser dollars, and a
teaser on another game changing nothing. `teaser.test.js` adds each leg's
`started` (the board's clock, else BFA's) and an empty result with no open
teaser.

The Edges-with-teasers smoke run (2026-09-30, a scratch Playwright harness,
not committed) loaded the unpacked extension with Sunday's saved NFL board
shifted to a future kickoff, Cal's 8 teasers built by the real BFA
open-bets normaliser and his straights of that morning: the Novig rows read
`teasers $600` / `add $0.00` / `$315.90 alone` on Colts +7.5, `held $400` ·
`against $445` · `teasers $800` / `add $425.65` on Chargers +6.5 and
`against $1040` · `teasers $860` on Commanders +6.5 (capped at its $415 of
liquidity), with each teasers line's tickets on hover, as in the approved
mockup.

The Teasers smoke run (2026-09-28, a scratch Playwright harness, not
committed) loaded the unpacked extension with Sunday's saved NFL board
shifted to a future kickoff and a bets payload built by the real BFA
open-bets normaliser: with nothing open the tab showed the plan's 7 tickets,
$1,292, +$324 (25%), 47% to make money; with ticket #1 and an outside $60
ticket (one leg on no board game, "counted as won") open, "Bet 6 more",
both tickets under Open at BFA with the pull age, and the legs with their
ticket and open counts, in light and dark.

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
sorting, and the alt path — keys, `mainPoints`,
`includeAlts` off by default, the distance / liquidity / same-number gates
following a moved main line, age through `sequenceNumber`, a pulled
ladder disappearing on the next parse, the rung venue ids on
`fixtures/v2_venue_ids_slice.json` (real CFB rows of 2026-09-15 with the
ids intact: per-event suffixes and id maps, the Novig rung on its main
number pointing at the main line, malformed ids ignored, a re-parse
replacing them), and `groupEdges` (card keys, book
and line counts, best by edge vs by stake, cards following their best).
`scanner.test.js` drives the loop with an injected fetch: snapshot-only
refresh (nothing goes to the retired stream), per-league failure, resume,
and a league switch while a snapshot is still downloading. `bets.test.js`
pins the Kalshi normaliser and the matcher on
`fixtures/bets/kalshi_fixture.json`; the pytest suite covers the service
(fills → positions aggregation, both spread signs and total directions, the
NO-moneyline tie caveat, unknown series failing closed, a failed poll keeping
the previous records, the `/bets.json` shape and `?days=` window, and node
parity on the same fixture) and the fill fairs: insert-only with the first
capture winning (a later capture and a repeat inside one request are
no-ops), served only for the bets in the window, the validator naming the
first bad row with what it expected and found, and `POST /fill_fairs.json`
behind the same JSON and Host guards as the crosswalk.

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
- **A held line still shows the ten-minute tag, not `since your … bet`**:
  no bet on that line has a saved fill fair. Bets open before 0.12.0 and
  bets placed while the panel was hidden or closed never get one (no
  backfill, by decision). For a new bet, the panel console logs each
  decision — `fill fairs: 1 captured` then `fill fairs: <bet id> saved`, or
  `not saved: … placed before the panel was watching` / `its line was first
  seen after the bet` / `no Unabated fair on its line at the fill` / `the
  board lists no market of its bet type and period`. A
  `service write failed: HTTP 404` means the bets service predates the
  route: restart it (`launchctl kickstart -k gui/$(id -u)/com.nflwork.bets-service`)
  and reload the extension; the
  captured row is retried every minute until the panel closes.
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
  (60 s–5 min). The row's line age comes
  from the alt's `sequenceNumber`.
- **No notifications**: they only fire while the panel is open and the
  toggle is on; the first pass after enabling is silent by design. Check
  Chrome's notification permission for the extension in System Settings.
- **Bets: "bets service unreachable since …"** (red header, Bets tab
  banner): nothing is listening on the service URL. Check the agent with
  `launchctl print gui/$(id -u)/com.nflwork.bets-service` (not loaded: install
  it, Bets service above) and read `bets_service.stderr.log` then
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
  spelling keeps failing, **Attach** the bet (the Bets tab's Needs a game):
  the attach teaches that venue's spelling through the crosswalk. Add an
  `ALIASES` row with a `teams.test.js` case only for a spelling every venue
  shares.
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
  in the header. **Attach** it to the right game: nothing is learned when
  both names already resolve, the pin alone decides the game.
- **Bets: "start time differs (bet …, board …)"**: the bet's game — the same
  two teams, or its rotation for a bet naming one team — is on the board
  within 12 h of the bet's start at another time: a rescheduled game or a
  venue clock read in the wrong zone. Attach it; if every bet of one venue
  shows it, fix that source's clock (whole hours off: suspect the clock
  first — BFA moved its open-bets clock once, and a DST slip is 1 h).
- **Bets: "Needs a code fix"**: the bets service could not read an open bet
  on a game league — the venue changed its text (BFA's two-bracket suffix,
  2026-09-26) or sent a form the parser never saw. The row shows the venue's
  own text and the parser's reason. Teach `bets_service/sources/<venue>.py`
  the new form with a pytest case and restart the service; the flag clears on
  its next poll. Attach cannot help: the record has no game, team or line. An
  unreadable bet that sits grey under Not on the board instead means the
  running service predates the marker (0.12.1): restart it.
- **Bets: "Attach failed: HTTP 404: no route for POST /pins.json"**: the
  bets service predates attach (0.11.0). Restart it
  (`launchctl kickstart -k gui/$(id -u)/com.nflwork.bets-service`). "HTTP 404: no bet with id …"
  means the service has not stored that record yet; wait for its next poll.
- **Bets: unmatched "by id: … not on any board ladder; by name: …"**: the
  bet carries a Kalshi / Novig id no board rung carries (the venue lists no
  ladder for the game, or the ladder dropped that number), and the name rule
  then missed for the reason after "by name". Fix the name reason; the id
  joins by itself once Unabated lists the ladder.
- **Bets: unmatched "no event on the board yet" / "league not on the
  scanner"**: the game is not in any loaded league snapshot — untick fewer
  sports in the Edges controls, or wait for Unabated to list it.
- **Teasers: a ticket just placed at BFA is suggested again**: it reaches
  the tab within about 90 s (the bets service reads BFA's open bets every
  60 s, the panel polls the service every 30 s); "Open at BFA" says how old
  the last BFA pull is. Much older means the service is down (red banner) or
  still on the old 300 s cadence — restart it after updating. If the ticket
  is listed there but a leg reads "no board game; counted as won", that leg
  did not join its game (BFA's start more than 12 h off the board's, or the
  game off Buckeye's board), so the ticket is not recognised as the list's:
  attach the leg to its game on the Bets tab.
- **Teasers: "Fewer than 4 games with a priced Buckeye leg"**: Buckeye's
  NFL/CFB board is thin (Monday night, or before the week's lines post), or
  its lines are older than **Max line age h** in the Edges controls.

## Design decisions log (moved from the root CLAUDE.md, 2026-09-15)

History of design decisions that used to live in `NFLWork/CLAUDE.md`. The sections above are the maintained reference; this log records *why* each choice was made and when, with issue numbers.

**2026-10-03 — Can't tease on college legs (0.17.0).** Some teaser legs the
tab suggested could not be bet: Buckeye keeps some college games off its
teaser menu, while every NFL game teases (user). Cal marks one by hand — a
chip on CFB legs only, mockup approved — rather than a rule, since no rule for
which games Buckeye drops is known. A mark covers the game's spread or total,
both sides, not the whole game: a bettable leg is never dropped by mistake,
and a game shut off entirely costs a second click. Marks end when the game
starts, so nothing carries into the next week. Checked on the 2026-10-03
board: marking Norfolk State +8.5 (in 4 tickets) rebuilt the list
without it, and Restore brought back the identical list.

**2026-10-02 — Bet105 reads time out after 30 s (0.16.2).** Bet105 stopped
pushing at 17:27 PT on 2026-10-01 while every other venue stayed fresh. Its
fetches had no timeout, so one request Bet105 never answered held the poll's
busy flag and every later 5-min read was skipped until the panel reloaded.
Each request now gives up after 30 s and the Bets tab shows "Bet105 did not
answer within 30s"; the next read runs 5 min later.

**2026-10-02 — Min edge default 2.5%, labelled on the toolbar (0.16.1).** The
minimum-edge box sat on the toolbar with no visible label, so Cal did not know
it was there. It now reads `min edge [x] %`, and its default moved from 1.0% to
2.5% (user). A value already saved in the panel is kept.

**2026-09-30 — The Edges tab sizes against the open teasers on a game.** The
reverse of the entry below, the last gap before the Teasers merge: a
straight on a game with open BFA teasers was sized as if they did not exist
(#130 leaves parlay legs out). A ticket's P&L on a row of the game also
depends on its other legs, and under the log objective that gamble cannot be
replaced by its average, so the choice was between the exact joint (enumerate
the other games, pool by which tickets survive) and a shortcut. Walked
through on Cal's 9/27 teasers (16:46 UTC board, K $5,000): every same-side
row goes to $0 whichever way (Colts +7.5 −270 $316 → $0 with three Colts
legs), so the method only matters on the other side — Chargers +6.5 against
four Bills −1 legs: exact $426, other legs averaged $492, other legs counted
as won $1,104. Cal's call: exact; open BFA teasers only (the Teasers tab's
own read; other parlays stay `parlay leg`); a `teasers $X` chip after the
approved mockup. A leg on another game that started counts as won — BFA
closed the five Seahawks tickets ten minutes after that game ended, so an
open ticket's finished legs won. The pre-merge review caught two
departures, both fixed before merge: a leg still to play with no board game
yet (CFB loading) counted as won too ($361 → $1,301 on its board) — now it
leaves its ticket out; and a decline the teasers alone cause threw out the
straights — now they still size the row. The other-side lift applies as for
straights (Commanders +6.5 $381 → $743), per the standing decision.

**2026-09-30 — Teasers sized around straight bets on the same games.** Cal's
last item before the Teasers merge ("incorporate previously placed bets").
Asked what he meant, he picked open straight bets — spreads, totals,
moneylines at any venue — on the games the teasers ride on (not every open
bet on any game: $14.6k open on NFL/CFB against a $5,000 Kelly bankroll
would leave the tab nothing to offer; no teasers or parlays are held outside
BFA). Same method as conditional Kelly (#130), no correlation number: the
straight's numbers cut its game's fair-ladder rows and its P&L is fixed in
the score, so the same side shrinks a leg and the other side grows it (the
#130 other-side lift Cal kept three times). Decided by walking 9/27 through:
the 49ers −14.5 straights ($414) made a 49ers ticket's loss −$614 instead of
−$200 and flipped the second pick to a ticket without the 49ers; all 15
straights cut the set from $1,292 to $676 and lifted its value after risk
from −$8 to +$62. Other market and other period left out, as on the Edges
tab (#129: Unabated gives no spread–total or cross-period link; all of Cal's
1H bets are totals anyway). Built: straights multiply a game's rows, so the
greedy now scores the 2^10 "which legs win" states instead of every outcome
(9/27: 69,984 outcomes, 164 ms vs 412 ms) and the budget went from 2^14 to
2^17; exact ties go to the first combination in order. The screen (held /
against / game tags on the legs, the note under the summary) was approved
from a render of the 9/27 board.

**2026-09-30 — Alt cap on price too, 15-85%; flat tails out; EV line gone
(0.15.1).** On the first live run of 0.15.0 the list showed Under 2.5 1H at
Kalshi +2242 (4.3c) as a "554% edge" in two CFB games: Unabated's fair for
that total was flat at +258 (27.9%) on every rung from 2.5 to 16.5, inside
the 10-90% fair cap. The user: it should be out anyway because it's 4.3c.
So the cap now applies to the book's price as well as the fair, raised to
15-85% (user), and a rung whose fair equals its neighbour's is dropped from
the list and from the tail-flex measurement (ladder.js already dropped it
from sizing). The `EV $` line on each row was removed as repetitive next to
the edge and the bet (user); it still ranks lines.

**2026-09-30 — Tail flex ranks alt lines; the 7-point cap becomes a 10-90% fair cap (0.15.0).**
The card's best line was the highest Kelly stake, which takes Unabated's
deep-rung fairs at face value; the Unabated host's own pick (-12 at 4% over
-24 at 10%, main -8.5) only holds if tail fairs are discounted. The user
kept Unabated as the truth for edge and stake and rejected blending its
fair with the exchanges ("I don't want to blend"), so the discount is
used only to rank: rank score = keep × edge × stake (EV dollars after
flex), keep = edge² / (edge² + sigma_e²), sigma_e growing with the line's
distance from the main number at a rate c measured live off two-sided
exchange alt quotes (fallback 10% under 100 rungs). The video's check
holds: -24 wins below c ≈ 5.5%, -12 above. The cap was first removed
outright; on one NFL+CFB board that moved suggested stake on best lines
priced -200 or shorter from $1.4k to $6.4k (that bucket ran -16.8% on 13
settled bets; +180 to +399, where the profit is, ran +19.6% on 281) and
added 79 longshot cards at +400 or longer with no track record. The user
kept a cap but in probability so it fits every sport (points mean
different things in an MLB total and a CFB spread; a share of the total
breaks on spreads; a max book price would hide the soft books' mispricing):
an alt lists only while Unabated's fair is 10-90%. c is the RMS of the per-rung
ratios with the top 1% dropped: the first build took the median, which
runs about a third low for a standard deviation and floored a quarter of
rungs at zero (NFL spreads 2.7% median vs 7.1% RMS — the median flipped the
video's pick). Removed: the **Max pts from
main** setting (stored values are ignored); `feed.altPassesGates` and the
live block now gate alts on the fair instead. Not shown on cards: the uncertainty
and the discounted edge; the header shows c per market and each row its
`EV $` rank score. No floor on a measured c (user decision). Card alerts
re-fire when the best line's rank score improves AND the line itself
changed (another line, or its price or edge moved) — the score also moves
with every re-measure of c, which is not news. Lines on Unabated's ±999900
fair clamp no longer list: removing the cap had surfaced Vanderbilt +46 at
-1800, "5.55%", a $2,497 quarter-Kelly stake, ranked #2 overall.

**2026-09-29 — BFA: a ticket the history cannot read still settles its
legs.** The open list stores a BFA parlay or teaser as one record per leg
(`bfa:<id>:leg0` …), and once it is graded only the history can settle those
rows. When the history could not read the ticket — its type names a leg count
it does not carry, or a leg does not parse — it wrote one unmatchable
`bfa:<id>` instead, which settled nothing: the leg rows stayed open for good
(the store is upsert-only and serves every open row), on the Bets tab and in
the Teasers tab's open tickets. Now such a ticket comes back one record per
leg up to the count its type declares (the open list has one row per leg),
each unmatchable with the reason and on the ticket's status: never a lone
record, never a partial ticket. A type the grammar does not know that carries
picks (an if-bet, say) takes the same path with the legs it carries; it used
to read as a straight bet on its first leg. While the ticket is open, the open
list's legs replace these one for one, which also removes the extra
unmatchable record listed beside them. Nothing to clean up: on 2026-09-29 the
live store's 61 BFA records (12 graded teasers and parlays) held no such
ticket.

**2026-09-29 — Bet105 is read by the panel, not the service.** Bet105 is a
LinePros white-label (the odds scraper already reads its socket) behind
Cloudflare, which challenges any request that is not the browser session's
own, so no script route exists and none is attempted. The panel, running in
Cal's own logged-in Chrome, makes the site's own `getHistory` call and pushes
the answer to the service (`POST /bet105.json`); the service parses and stores
it like every other venue so pins, attach and freshness work unchanged. Open
bets only (the site shows no settled list); a bet a complete push no longer
lists is marked `closed` with no result. The selection ids (`marketId`,
`key`, `subKey`, `periodId`) are read the way `bet105_odds/scraper.py` reads
the same feed, and only the captured shape (NFL first-half totals, side `1`)
is verified; the rest is pinned by hand-written fixture rows and fails closed
on anything unseen. Rejected: a Pikkit-mediated read (its terms forbid it and
its data is only as fresh as the phone app's last sync).

**2026-09-28 — Buckeye teasers (the Teasers tab, 0.16.0).** Cal bets Buckeye
6-point 4-team teasers at BFA and wanted the set built for him. Decided with
him 2026-09-27/28 (plan: `docs/2026-09-27-unabated-ticket-teasers-plan.md`):
Buckeye only; legs priced off **Unabated's fair, not the exchanges** —
decided right after being shown that at teaser numbers Unabated sits up to 3
points above Kalshi / Novig / Polymarket / ProphetX, for consistency with
every other edge in the panel; **$200 a ticket** (Buckeye's limit), hence a
set of tickets rather than one; **4 teams only** (+300 is Buckeye's best
teaser price); **a push loses** (conservative while Buckeye's tie rule is
unconfirmed), which makes every leg win-or-lose at the half-point one step
against it; **one leg per game** (independent legs); the panel's own
bankroll and multiplier; **BFA is Buckeye**, so placed tickets are read from
the bets service's BFA open bets and there is no Placed button to press —
the approved mockup's Placed and Undo are gone; **no per-leg cap** (declined
after seeing "no leg on more than half the tickets" cost 6% of the value
after risk; Sunday's Seahawks leg on 5 of 8 tickets sank 5 of them while the
other 7 legs covered); the list rebuilt only on real changes or a 1-point
leg move (near-tied combinations reshuffled a rebuild-every-tick list for no
value). Built: the list is built on reference fairs that move only by a
point, so a rebuild on any other trigger — an open teaser appearing — does
not reshuffle it, and the list after placing its first ticket is the same
list less that ticket. The greedy is the prototype's: a full $200 ticket
while one still raises the score, then one partial ticket at its own best
stake (on 9/27 its sixth $200 ticket is past that ticket's own best, $129,
yet the set is worth $141 after risk against $140 for stopping there at
$129). The pre-merge review found a placed ticket coming back: open
teasers entered the score only as P&L, so on a board with four standout legs
the ticket just placed was the best $200 ticket again, and following the list
stacked one combination four times. Now an open ticket's own four legs are
never offered again, and once a ticket under $200 is open no partial one is
(else placing the list's partial last ticket brought a smaller one, $42,
$25, $13, $5). Past 2^14 joint outcomes the pool sheds free legs instead of
the tab going dark. A second review found the same bug through a leg that
missed the board join (BFA's clock off, a name resolving elsewhere): counted
as won, it left its ticket unrecognised — now such a leg joins by rotation
within 12 h. It also found the frozen list lost on a panel reopen or a
scanner restart (now the reference fairs persist and the tab waits for both
football boards), a 5-team or 2-team teaser tripping the combination and
partial rules (now 4-team only), and the Open at BFA wording and totals. BFA's open bets moved
from every 300 s to every 60 s (history stays at 300 s) so a placed ticket
shows within ~90 s. The scanner always loads NFL and CFB; `selectEdges`
gained `leagueIds` so the Edges tab still lists only the ticked sports.

**2026-09-26 — Novig's `unabatedId` no longer beats the team name.** An
open Novig over on Southern @ Jackson State never matched its game. Novig
sends Jackson State `unabatedId` 276, which is no CFB team (the board's
Jackson State is 777), and `resolveTeamKeys` put the venue id ahead of the
name, so the bet looked for a team that is not on the board. Measured
against that day's board: 5 open CFB sides carried an id that is no CFB
team (Jackson State, UTEP, Prairie View A&M, Eastern Washington), and
every other Novig side (38 CFB, 38 NFL) agreed with its name. The id never
rescued a name that failed. Now the name goes first, and the id keys a side
only when the name resolves nowhere and the id is a team of the bet's
league (`teams.hasTeam`). A learned crosswalk row still wins over both.

**2026-09-26 — Loud when an open bet cannot be read, or its clock is off.**
The Eastern Illinois bet (below) exposed two ways an open game bet fell off
the board with the Bets tab quiet; Cal had to notice it himself. (1) A
source's parse failure read exactly like a deliberate exclusion: both were
`unmatchable`, grey under Not on the board. The sources now mark
`raw.parseFailed` where the venue's own league code named a game league and
the selection did not parse, and the panel flags those red in a new **Needs
a code fix** block (no Attach: the record has no game, team or line) until
the parser is fixed or the bet settles. A marker from the source rather than
string-matching reasons in JS (Cal's preference): only the source knows
which of its reasons are on purpose (BFA: a team total; Wagerzon: a
postponed leg). Wired for BFA open bets and Wagerzon legs, the two
description grammars with a per-leg league code; BetOnline reads its league
off team nicknames and Kalshi, Novig and Polymarket read structured
payloads, so they stay unmarked until one breaks the same way (a
`mark_parse_failed` call). Rows with no parsed pick now read the venue's
text instead of their id. (2) BFA and BetOnline spreads and moneylines name
only their own team, and **start time differs** needed both team keys, so
the bet with its clock 7 h off read "no event on the board yet". The rule
now also takes the bet's rotation on a not-started board event within 12 h
where its named team fits (`knownTeamFits` — BetOnline's reused rotations
are a week apart, outside the window anyway). A rotation match also lets
the board's clock say "not started": that bet's own start, 9 AM PT, was
before it was even placed, so the old started-check would have kept it grey
too; a match by team pair alone keeps the bet's clock, so a finished
doubleheader game 1 never flags against game 2. Cal approved the mockup
before the build.

**2026-09-26 — BFA open bets: a two-bracket suffix, and the open list is on
Pacific time.** The first live BFA open bet (Eastern Illinois +35½, CFB) never
matched, for two reasons. Its description ended `[Sport:Football][League:NCAA]`,
two brackets, where the February capture had one (`[Sport:Basketball,
League:NCAA]`). The suffix stripper only knew the one-bracket form, so the leg
failed to parse ("unrecognised selection"). Behind that, the open list's
timestamps are the account's Pacific wall clock, not the "+7 h" clock the
February capture showed (that response was dated 18:50 PST but listed a bet
placed "22:33"). The 7 h subtraction put the 4 PM PT kickoff at 9 AM PT, so
even with the suffix fixed the bet would have missed the board's 23:00 UTC
start. Both halves of the source now read one clock (`parse_account_time`).
BFA has moved this clock once already; if it moves again, bets will miss the
board by whole hours.

**2026-09-23 — Why an edge grew: since your first fill (#132 follow-up).**
The trigger was a Novig card: 49ers -14.5, an alt of -8.5, held for $414 over
three fills at +208, +217 and +212, then showing +223 and `add $77.25` at
+12.15%. Before adding, Cal wanted to see whether Unabated's fair had risen
with the position or fallen while the price improved (the picked-off
pattern) — which the tag's ten-minute window cannot see across fills hours
apart. The tag on a line held in the same direction now reads from the
earliest open bet on that line whose fill-time fair was saved: the panel
reads the fair off its in-memory history at `placedAt` and the bets service
keeps it once per bet (`bet_fill_fairs`, `POST /fill_fairs.json`). User
decisions: (1) display only on held lines — the existing `add` gate; (2) save
for every new bet, never backfill — nobody knows at fill time which line will
pop back up, and a fair the panel did not observe is not saved, so bets open
at ship time and bets placed with the panel hidden or closed keep the
ten-minute tag; (3) the stake is untouched — a hold rule ("no add on red")
stays a separate decision; (4) no knobs. My calls: the four cases stay one
rule (`edgemove.classifyMove`, called by both paths); the price leg compares
the fill's own price on edgemove's basis (exact cents for exchanges), and only
on the bet's own book — across books it would read a venue spread as the book
moving away; the
scanner's `observingSince` starts when the catch-up after a resume has
landed rather than when `resume()` is called (until then the history is the
pre-pause board) and restarts after a two-minute gap with no successful
poll; a fallback to another book's line at the number, because the fair is
the same across books there (NFL and MLB 100% of multi-book rungs, CFB 95%,
measured 2026-09-23); and `/bets.json` serves only the fairs of the bets in
its window, since the table is never pruned. Known and accepted: Kalshi's
price on Unabated includes the taker fee while the fill does not, so on
Kalshi the price leg reads up to ~1.75¢ worse — it can hide an amber, never
raise one; a bet placed just before the panel hid and reported after it
reappears gets no fair (the rule reads the last resume, not each hidden
interval).

**2026-09-23 — Polymarket US source in the bets service.** Cal's account is on
Polymarket US (the CFTC app), not polymarket.com, so the public
wallet-address data API that would have served the international site does
not apply. Chosen over borrowing the app's session: the documented API key
from polymarket.us/developer, which Cal created and keeps in
`bet_logger/.env`; the service signs with it and nothing else sees it. No
title parsing: the live pull showed every trade carries its full market
object (type, line, both sides with their own signed number, team ids and
away/home), so the type and side come from structured fields and only the
team names and start come from a public event lookup, filtered to the traded
market type because the whole event runs to 4.7 MB. Verified on the live
pull of 2026-09-23 (1 open position, 26 trades, 5 settlements, one 3-leg
combo): all seven records' stakes and results agree with the venue's own
realized P&L. Combos fail closed as one record for now (plan decision), though
their legs carry slug, side and state. **Unobserved:** a NEUTRAL settlement
(served `unknown`), a busted trade, a CBB or NBA event (spelling rule
assumed from CFB and WNBA). Poll 60 s and the 31-day window are constants.

**2026-09-23 — Wagerzon (the C account) source in the bets service.** Second of
the Bet Logger books Cal asked for. Same shape as BFA: a form login this
process owns in memory (`sources/wagerzon.py`), history plus the open-bets
helper every 300 s, pending bets kept. The live capture that night held only
props (NFL specials, RBL touchdown props) and one transfer, so the game-bet
grammar — MLB totals with the pitchers bracket, first-five run lines with the
`1H` token before the name and no opponent, `GM#1`/`GM#2` doubleheaders,
`EV` prices, postponed legs — was taken from every bet_logger run since
2026-04 (`bet_logger/logs`) and pinned by a hand-written fixture week; the
matcher gets a real start time here (`GameDate` + `GameTime`, Eastern) and
the league from `IdSport`, so unlike BFA nothing depends on a nickname scan.
**Unobserved:** `OpenBetsHelper.aspx` answered `{"result": []}` — a row's
shape is unknown, so a row that is not a HistoryHelper-shaped wager fails the
poll naming its keys (loud, never a silent drop); the history already lists a
pending bet with an empty `Result`, so the helper may add nothing. Considered
and rejected: skipping the helper until a shape is seen (an open bet the
history did not carry would go unflagged with no trace). Later that day, on
Cal's "we just need the open bets right": the row shape was read off the
site's `ui_c` bundle (`OpenBetsTable.loadOpenBets` groups the helper's rows by
`TicketNumber` and reads the fields listed in the source bullet), so the
helper is now the primary source of open bets and the history only settles them.

**2026-09-23 — BFA (Betfastaction) source in the bets service.** Cal asked for
the Bet Logger books on the Bets tab: BFA first, the Wagerzon C account next,
Polymarket after a question on which Polymarket. The source owns its own
Keycloak password session in memory (`sources/bfa.py`) rather than sharing
`bet_logger/recon_bfa_auth.json`: a shared refresh token needs a file lock
(BetOnline, #115) and still dies when the LaunchAgent and the service race a
rotation. Built on a live capture (`bets_service/recon_capture.py`, run on the
Mac because the cloud thread cannot reach the book): 16 wagers over 30 days —
9 straight, 3 free-play, 3 four-team teasers, 1 two-team parlay, none pending.
Findings that shaped the parser: `settledDate` is the scheduled kickoff, on the
account's Pacific clock, on every settled straight bet (noon-ET games read
09:00), so it serves as the event start; the rotation-parity convention holds
(every over odd, every under even); the description names no sport and
`parse_sport` cannot tell CFB from CBB, so college bets — the bulk of the
account (1H totals) — land in the unmatched list as "league unknown", the
hand-off's rule (never guess). Later that day, on Cal's "we just need the
open bets right": `bet_logger/recon_bfa_api.json` (the February browser
capture) holds `GetPlayerOpenBets`, whose legs carry `idSport` and
`gameDateTime` — so an OPEN college bet is placed after all, and the source
reads that endpoint first with the history only settling records; its clock
is the Pacific wall-clock plus 7 hours (two February rows), true UTC only in
summer. The college question now touches settled bets only. Measured the same
evening against the live board: every BFA NFL (15) and CFB (21) spelling from
the pull keys through `teams.js`, every Wagerzon MLB spelling but "ARI DBACKS"
(now an `ALIASES` row), and both books' rotations are Unabated's numbers — BFA
prepends a "1" on first-half legs, which the source strips. CBB spellings are
untested until Unabated lists the teams (2 today). Credentials resolve
from `bet_logger/.env` as `config.py`'s fourth lookup place, so no password is
copied. Poll 300 s and the 31-day window are constants, not settings.

**2026-09-23 — Liquidity: a floor on what the resting money can win, and no
stake above it.** Min liquidity ($100) gated alts only, so a Novig Portland
Fire +809 main moneyline with $17 resting listed with "add $32.58". A flat
stake floor on every line was rejected: it hides longshots whose whole Kelly
bet is small. Cal's rule instead: $20 on +2000 (wins $400) is worth a look,
$20 on +100 (wins $20) is not, so **Min liq to win** ($100) gates main lines
and alts on liquidity × the price's win per dollar, and every suggested
stake (rows, cards, sort, alerts, the Ticket) is capped at the line's
liquidity. Unverified: that Unabated's liquidity is stake dollars, not payout.

**2026-09-22 — Why an edge grew: the fair decides (#132).** The panel sized on
the current edge every refresh and could not say why it was what it was;
on a held line the two reasons an edge grows call for opposite actions
(fair moved to you: add; the book moved away while Unabated's fair, ~1–2
min behind, has not answered: wait). The scanner now keeps a per-line
history in memory — one observation per snapshot and per stream update,
distinct values only, a number move resets the line — and `edgemove.js`
reads a tag off it inside a 10-min window at a 0.5-point threshold in
probability. Spec addition the same day: four cases, the fair's direction
picks the tag (`fair moved to you` green, `book moved away` amber, `fair
moved against you` red — the adverse-selection case, red whatever the price
did — none otherwise); "both moved, name the larger" was dropped. Openers
(`openerPrice` / `openerPoints`, #126) are parsed for the first time and
shown as tooltip context only: measured live, per book on every main line,
none on alt rungs, and 54% of NFL spread/total lines sit on another number
than they opened. Acting on the tag (a hold rule) is out of scope; the stake
is untouched. Plan: `docs/2026-09-22-unabated-ticket-edge-source-tag-plan.md`.
Narrowed 2026-09-23 (user decision): the tag is shown only next to lines
held in the same direction (the `add` rows, card lines and Ticket) — it
exists for the adverse selection of adding, and a first bet is not a
top-up — while the history is still kept for every line.

**2026-09-22 — Novig source rebuilt on the Portfolio REST feed (#116).** Novig
moved its web app from `app.novig.us` to `novig.com` and put its Hasura
GraphQL behind a query allowlist the same hour; the service's hand-written
`order`/`parlay` queries, and the old app's own Portfolio queries, all answer
`query is not allowed`. The new app reads `GET api.novig.us/nbx/v1/portfolio/{active,settled}`
on the same Auth0 client and audience the service already mints for, so
`sources/novig.py` now does the same — pure HTTP, cursor-paged, no allowlist
in the path. Considered and rejected: a capture-and-replay of the app's
allowlisted GraphQL documents (an arms race Novig had just shown it was
running, and moot once the app itself left GraphQL) and keeping the
content-script mirror as the primary source (the Portfolio screen pages 15
at a time; 91 open positions / $16k at risk would need seven "load more"
taps per refresh, and a position outside the window would go stale rather
than missing — worse for conditional Kelly). The mirror (`novig_page.js`,
`novig_content.js`, `novig_bets.js`) and the parity test went with the old
app. The REST shape validated against the 476 GraphQL-era records: side
agreed on all of them; the old normaliser's lay flip had priced one order at
the wrong side's 0.70 where the feed's 0.305 was paid. The Novig paragraph
below is the pre-rebuild text, kept as history.

**Unabated Ticket** (`unabated_ticket/`) — Chrome MV3 side-panel extension (plain JS, no build). A capture-phase click on an Unabated odds-screen price reads the React fiber / AG Grid row (`page.js`, MAIN world), builds a one-at-a-time ticket (side, points, book price, `bacr` fair) in `chrome.storage.local`, computes the unrounded quarter-Kelly stake (`kelly.js`, node-tested), and re-reads the line every 5 s for a line-moved warning — **the watcher resumes from the stored ticket** whenever `page.js` loads (it dies with every navigation: one-click betting leaving the tab, Back, reload, discard, extension-reload takeover; `content.js` hands the stored ticket back on `resume_request` and once on its own load, and `page.js` rebuilds the watch by identity with no grid API — before 0.6.6 every such ticket read "Not watching the line" under a live heartbeat). **Edges tab (issue #112)**: while the panel is open, `scanner.js` reads Unabated's public feeds — the v2 league snapshot (`content.unabated.com/markets/v2/league/{lg}/odds.json?t=<30s bucket>` — the query busts a CloudFront edge cache that served gzip clients a 6.5 h-old copy on 2026-09-10 — on open + per-league refresh, 4 at a time, for the 29 team-sport leagues in `feed.LEAGUES`: NFL/CFB/NBA/CBB/WNBA/MLB/NHL + soccer; tennis and combat key sides on people and are out) and the changes stream (`api-k.unabated.com/api/markets/changes/query/{cursor}`, every 10 s; cursor = ns since 2021-01-06, kept as a string; **incomplete anonymously** — 69 of 191 NFL line changes in 3 min on 2026-09-10, exchange moves mostly missing — so each league's snapshot re-downloads on a size-tiered cadence, 60 s / 2 min / 5 min) — through `feed.js` (pure, fixture-tested: lines keyed `(marketId, book, sideKey)` because the changes stream tags other markets of an event with the same `bt` key, updates applied only on a newer `sequenceNumber`, books listed when `isActive` and enabled for game odds — `statusId` is not liveness, Caesars/Underdog carry 2 while live) and lists every ML/spread/total with Unabated's `ge` at or above the minimum and a `modifiedOn` within the max line age (default 168 h — dead feeds at "active" books carried 96-day-old lines with +36% "edges" on 2026-09-10; each row prints its line age), sized with `kellyStakeFromEdge`, filtered by a Books multi-select dropdown + Bets checkboxes in the panel (books default to the selection `page.js` publishes from `userSettings.gameOdds`; the user's own ticks win). Row click / notification click store a `locate` request (`locate.js`) that focuses the Unabated tab and has `page.js` scroll to the row and outline the cell — the price click stays the user's. Alerts are Chrome notifications, off by default, baselined on enable, deduped per line with a 5-min per-event cooldown. **Alt lines (issue #113)**: each snapshot line's `alternateLines[]` expands into lines keyed `(marketId, book, sideKey, points)` under the main line (`isAlt`, `mainPoints` = the book's own main number — NOT the feed's `stn`, which is the market's standard number), listed only behind an **Include alt lines** toggle (off by default) with two alt-only gates, **Max pts from main** (7) and **Min liquidity** ($100, exchanges only), on top of every main-line gate; an alt sitting on its main line's current number is hidden. Freshness caveat measured 2026-09-11: the anonymous changes stream carries NO alt updates and every alt's `modifiedOn` is the `0001-01-01` sentinel, so alts refresh only with the per-league snapshot and their age reads from `sequenceNumber` (the change time in epoch ms). An alt row click expands the grid row's Alts (`setExpanded`) and finds the cell by fiber props (points/side/book/event), failing loudly with the main cell outlined; a ticket captured on an alt cell carries `watch.altPoints` so the watcher re-reads the same rung. **Group by market** (on by default): `feed.groupEdges` collapses the list to one card per (game, period, bet type, side) — a +EV opinion is directional — whose best line is the highest Kelly stake (taxes longshots, so a -110 main at +5% beats a +944 rung at +6%), with `N books · M lines` behind an expander; alerts then key on the card and re-fire only when its best edge improves. **Bet history flags (issue #114, core; merged 2026-09-11)**: a local **bets service** (`unabated_ticket/bets_service/`, `run.sh`, `127.0.0.1:8094` loopback-only, no auth) is the only place that signs Kalshi requests — the private key never enters the extension — and turns the account's fills + unsettled positions into normalised bet records (one per `(ticker, side)`; positions are the truth for open size, fills the VWAP entry price; `sources/kalshi.py`, `sources/kalshi_ticker.py` a series map + cached public `/markets` and `/events` GETs, football suffixes carry the Eastern DATE only, MLB `HHMM` too; unknown series and futures fail closed as "unmatchable", shown never matched), UPSERTs them on native id into `bets_service/bets.duckdb::bets` (never pruned — future CLV) with a `source_runs` row per poll, and serves `GET /bets.json[?days=30]` (`{generatedAt, sources: {kalshi: {fetchedAt, ok, error, count}}, bets}`) + `/health`; a failed poll keeps the previous records (a dark source never blanks the list), a source on its first poll reads `ok:false, "no completed poll yet"`. The panel polls it every 30 s while visible — **never from the service worker** — into `chrome.storage.local` `betsService` (payload + records, open + settled ≤30 days, ~1 KB each) and `betsSettings` (`serviceUrl`), resolving team keys through `teams.js` — NOT a hand table (user decision 2026-09-11): every league snapshot's own team list (id, name, abbreviation) is registered at runtime and persisted as `teamsIndex`, keys are `<league>:<Unabated team id>` (the ids the board's lines carry), venue spellings resolve by exact normalised name, then a three-row `ALIASES` list, then a leading-code strip ("PIT Steelers"), then a UNIQUE word-boundary containment (the query-extends-team direction only for State/University/College, so "Southern Mississippi" never becomes "Southern") — and treating a venue's `ok` pull as authoritative for that venue's records. `bets.js` (pure, node-tested on `tests/fixtures/bets/kalshi_fixture.json`; `normalize_kalshi` in Python is held byte-equivalent by `test_parity.py`) matches a bet to a line on league + team pair (either order) or rotation + time (≤30 min with a start; the Eastern date ±1 day without), refuses to guess when two board events fit ("ambiguous game"), and ranks four tiers: `same_line` (market, period, side and number), `same_side` (different number), `opposite` (the other side — red), `same_game` (any other market/period); Kalshi NO on a team market = the other team **or a tie** (NFL/CFB/soccer, `approx: kalshi_no_side_includes_tie`). Surfaces: a **Ticket banner** (strongest first, 5 then "+N more"), **Edges badges + sized stakes** — a held bet NEVER hides a line, it changes the size of the next one (user decision 2026-09-11, replacing a hide toggle): `bets.exposureOf` sums dollars risked on the same direction (`held`, same_line + same_side) and the other side (`against`), `betsview.stakeAdvice` turns the row's Kelly stake into one wording everywhere (user choice 2026-09-11, `stakeAdviceWords`): "wagered $350 → target $600, bet $250" (against positions read "wagered $200 against …, bet $500 (net $300 on this side)", full size reads "bet $0"), rows and cards carry a `held $N` / `against $N` (red) / `game` badge plus a dim "you hold …" position line, a **by my exposure** sort, and alerts skip only lines already held at size, a **Bets tab** (per-venue freshness green <5 min / amber <60 / red, open bets, and an **unmatched** list with a reason per bet — nothing dropped silently), and a **header line** under the tabs ("bets: N open · kalshi 20 s · betonline — …"). Venue split: **#115 BetOnline (merged 2026-09-15)** polls the account's bet-history report through the bets service on the bet_logger Keycloak refresh token (`sources/betonline.py`; registered when the recon cookie file exists; refresh only within 60 s of expiry under a file lock shared with `bet_logger/scraper_betonline.py`, pending bets kept, league via `bet_logger/utils.parse_sport`, no game date → `approx: game_date_unknown` and the matcher windows `placedAt` −12 h/+14 d keyed on the rotation, with `tierByRotation` for team names `teams.js` cannot key); **#116 Novig, #117 ProphetX** each add a `Source` (`sources/__init__.py` protocol: `name`, `poll_sec`, `fetch() -> list[record]`, raise on failure, never partial) or a content script writing storage directly; until they land those venues read "no source configured". **Novig (issue #116)**: the official NBX API is a separate "Liquidity Provider" account ($30k minimum deposit; docs.novig.com/lp-onboarding) and cannot see retail-app bets, so the primary source is `bets_service/sources/novig.py` — a one-time `novig_auth connect` (Auth0 PKCE against the app's public client, the user logs in, refresh token saved to gitignored `novig_token.json`; its OWN chain, never the app's rotating localStorage token) then a 60 s poll of the app's Hasura `order`/`parlay` tables for the trader resolved from the JWT `sub`, through `curl_cffi`; `normalize_novig` is held byte-equivalent to `novig_bets.js` by `test_parity_novig.py`. The content-script mirror is the token-free fallback: `novig_page.js` (MAIN world, `document_start`, wraps `window.fetch` before the app bundle captures it) mirrors the three Portfolio queries the `app.novig.us` app fetches for itself (`ActivePortfolioOrders_Query` / `SettledPortfolioOrders_Query` / `ParlayPortfolioQuery` → `api.novig.us/v1/graphql`; zero requests of its own, no token read) and `novig_content.js` normalises them via `novig_bets.js` (outcome index 0 = home/Over, `price` a 0–1 probability, `isBid:false` LAYS the outcome = the other side at `1 - p`, `_1H` on MLB = `F5`, props/futures/unsupported leagues fail closed) into `chrome.storage.local.betsNovig`, which the panel merges (a complete read is authoritative for the venue; a stale one says "open app.novig.us … to refresh"); shapes are bundle-derived, not a live capture — a divergent blob lands in the unmatched list as "unreadable Novig order (…)". No overlay, no order placement. Load unpacked from `unabated_ticket/extension`.

**Conditional Kelly against held bets (issue #130, 2026-09-19).** The panel sized every line alone and adjusted by subtracting dollars (`add $X`, `still $X against`); dollars at different prices and numbers are not comparable and a moneyline held against a spread was ignored (four live cards, K $8,000: New Mexico @ Oklahoma bet $183 → $296, Lions @ Bills "at full size" → add $188, UConn add $116 → $121, Chargers ML bet $189 → add $0). Now `betsview.stakeAdvice` sizes the row given the open bets on its market (`condkelly.js`, `ladder.js`): K = bankroll × multiplier (full Kelly on the whole bankroll then a quarter is wrong — it says add $500 to a line already held at full size $667); the new bet's chance from Unabated's edge (decision 2026-09-09 stands), a held bet's from the median `bacr` at its half-point rung today; no fees added (Unabated's prices include them, user 2026-09-17). Same market and period is exact off the ladder; another period is the worst case the fairs allow (comonotonic pairing) because no correlation is estimated (user decision 2026-09-18: no per-sport history, use what Unabated provides) — **this reverses the 2026-09-12 "show it, don't size off it" decision** on the unmerged `feature/unabated-related-period-bets`, whose `related_same` / `related_opposite` tier names were taken and which is retired; **only the same direction is sized across periods** — the logic review (2026-09-19) showed the literal rule sizes a cross-period hedge as if both bets won together (hold 1H Under $300, new FG Over: $407 vs $667 alone vs $802 under the measured 0.69 link), so the other direction in another period is left out and named `other period · not sized` (user decision 2026-09-19); another market is never in the math (worst case gave $35 vs $188 on the Bills card and independent differs from ignored by $1; open in #129). Guards are left out and named, never guessed. Found while building: the changes stream ADDS lines for known events and files team totals under the game total's `bt3`, so `feed.js` marks snapshot lines `fromSnapshot` and only those feed a ladder. `bets.exposureOf`, `badgeText` / `badgeKind` and the `add` / `at_size` / `reverse` advice kinds are gone. **Hedge sizing is deferred** (user, 2026-09-19): a price with no edge of its own stays $0 and there is no `hedge only` label. Plan: `docs/2026-09-18-unabated-ticket-conditional-kelly-plan.md`.
