# Unabated Ticket

Chrome extension (Manifest V3, plain JS, no build step). Click a price on the
Unabated odds screen and a side panel shows the bet — side, points, book
price, Unabated fair price, edge — and the quarter-Kelly stake. The panel
stays open when the sportsbook tab opens, so the stake is in view while you
place the bet. Issue #111; plan in
`docs/2026-09-08-unabated-ticket-extension-plan.md`.

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
  script. A capture-phase click listener runs first; Unabated's own
  deep-link handler runs untouched afterwards.
- Fast path: fiber props. Fallback: `data-marketline-id` on the shell plus a
  `forEachNode` scan of every row's `sides`. If both fail the panel says
  **Could not read this cell** with the reason; it never shows a stake it
  cannot back.
- Fields: `price` = `americanPrice` (exchanges have only `price`), `fair` =
  `marketLine.bacr` (Unabated's no-vig price at that book's points),
  side 0 = away / Over, side 1 = home / Under.
- `content.js` forwards the ticket to `background.js`, which writes
  `chrome.storage.session` (cleared when Chrome closes) and opens the panel.
- Every 5 s `page.js` re-reads the same book line through the grid API. If
  price or points moved, the panel shows **Line moved**, re-sizes off the
  new price and fair, and keeps the captured line for comparison. Off the
  board shows in red. If the Unabated tab is closed or navigated away the
  panel says **Not watching the line**.

## Stake

Mode A of the Kelly sheet (`extension/kelly.js`):

```
p_fair = prob(fair)               fair is already no-vig
b      = decimal(price) - 1
full   = (p_fair * b - (1 - p_fair)) / b     0 if no edge
stake  = bankroll * full * multiplier        not rounded
```

Settings (bankroll, Kelly multiplier) sit at the bottom of the panel and
persist in `chrome.storage.local`. Defaults 30000 and 0.25.

Copy puts one line on the clipboard:
`LSU Tigers -22.5 -110 @ DraftKings | fair -106 | edge +1.23% | stake $412.50 | LT @ LSU`.

## Tests

```bash
node --test unabated_ticket/tests/kelly.test.js
```

Checks the sheet's worked example (-400 vs fair -900, bankroll 30000,
quarter Kelly = $3,750), negative edge → $0, and that nothing rounds.

Manual checklist after loading unpacked: click a best-line price and a
book-column price, then a moneyline, a spread and a total; confirm side
label, points sign, price, fair and stake; change bankroll and watch the
stake move; wait for a line change and see the warning; close the Unabated
tab and see "Not watching".

## Troubleshooting

- **Panel empty / "Click a price"**: nothing captured yet, or Chrome was
  restarted (session storage clears). Click a price again.
- **Could not read this cell**: Unabated changed prop or class names. Check
  `.odds-cell-action-shell` still exists and the fiber props still carry
  `marketLine` / `sideIndex` (see the DOM notes in the plan doc); the
  reason text names which lookup failed.
- **Panel did not open on click**: Chrome only auto-opens the side panel
  with a user gesture attached; click the toolbar icon once, it stays open.
- **"Not watching the line"**: the Unabated tab is closed, navigated away,
  or the row left the grid (filter change). Re-click the price.
- **Cannot size: no Unabated fair at the new line**: the line moved to
  points Unabated has not priced yet (`bacr` missing). Wait a tick or re-click.
