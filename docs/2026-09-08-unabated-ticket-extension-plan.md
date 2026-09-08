# Unabated Ticket — Chrome extension plan

**Date:** 2026-09-08 · **Branch:** `feature/unabated-ticket` · **Worktree:** `.worktrees/unabated-ticket`

## Problem

Clicking a price on the Unabated odds screen deep-links to the sportsbook in a
new tab. By the time the book loads, the bet (side, line, price) is out of
sight, and the stake still has to be retyped into the Kelly sheet. The
extension remembers the bet and computes the stake at the moment of the click.

## Scope

**In:** one ticket at a time, shown in a Chrome side panel that stays open when
the book tab takes focus. Stake from the Kelly sheet's Mode A. Line-moved
warning while the ticket is open. Copy button. Settings: bankroll, Kelly
multiplier.

**Out (user decisions 2026-09-08):** bet tracking (Pikkit), on-screen overlay,
placed tags, recent-tickets list, exposure-aware Kelly, voice, phone.

## Architecture

Chrome extension, Manifest V3, plain JS, no build step. Lives in
`unabated_ticket/` (new top-level directory).

```
unabated_ticket/
  README.md
  extension/
    manifest.json        MV3; host permission tools.unabated.com; sidePanel + storage
    content.js           runs on tools.unabated.com/*/odds; captures clicks, watches the line
    kelly.js             pure functions: american->decimal, fair prob, stake (shared)
    background.js        opens the side panel on toolbar click; relays messages
    panel.html / panel.js / panel.css   the ticket UI
  tests/
    kelly.test.js        node --test; checks against the Kelly sheet's worked examples
```

### Data flow

1. `content.js` installs one capture-phase click listener on `document`.
   A click whose target is inside `.odds-cell-action-shell` reads the cell's
   React props (`marketLine`, `sideIndex`, `betType`) and the AG Grid row
   (`node.data`) and builds a **ticket**. Unabated's own handler runs
   untouched afterwards.
2. Ticket goes to `chrome.storage.session` (cleared when Chrome closes).
   The panel listens to `storage.onChanged` and re-renders.
3. `content.js` keeps a 5s watcher on the captured line: it re-reads the same
   row via the grid API and, if the book's price or points changed, writes
   `current` next to `captured` on the ticket. Panel shows the warning and the
   re-sized stake.
4. Settings live in `chrome.storage.local` and persist.

### Ticket contract

```
{
  capturedAt: 1788882208262,
  league: "cfb", eventId: 123629, eventStart: "2026-09-12T23:30:00" (UTC, naive),
  eventName: "Louisiana Tech Bulldogs - LT @ LSU Tigers - LSU",
  betType: "Moneyline" | "Spread" | "Total", sideIndex: 0|1, sideLabel: "LSU" | "Under",
  points: null | 55.5,
  book: { id: 89, name: "Novig" },
  price: -49900,           american, book's price
  fair: -5854,             american, Unabated Line at the book's points (marketLine.bacr)
  edgePct: -1.48,          Unabated's own number, display only
  current: { price, points, seenAt } | null
}
```

Field sources are documented in memory `unabated_odds_screen_dom.md`.
Exchanges (Novig, Kalshi, ProphetX) have no `americanPrice`; use `price`.

### Kelly (Mode A of the sheet)

```
dec_book  = american_to_decimal(price)
p_fair    = american_to_prob(fair)          # fair is already no-vig
b         = dec_book - 1
full      = (p_fair * b - (1 - p_fair)) / b # 0 if <= 0
stake     = bankroll * full * multiplier               # shown to the dollar
```

Shown: stake (large), full-Kelly dollars (small), edge %, fair. Negative edge
shows stake $0 and a muted "no edge" line. Defaults: bankroll 30000,
multiplier 0.25 (from the Kelly Calculator sheet). No rounding (user decision).

### Side labels

Moneyline: team name from `eventTeams[sideIndex]` resolved through
`context.fullOddsData.teams`. Spread: team plus signed points
(`-22.5`). Total: `Over`/`Under` plus points. Rotation number shown small.

### Fragility and fallback

Reading React fiber props is the fast path. If the prop names change,
`content.js` falls back to `data-marketline-id` on the clicked shell and a
`api.forEachNode` scan of `node.data.sides[*][ms<id>]` for the matching
`marketLineId`; the grid API is reached from any cell's fiber. If both fail,
the panel shows "could not read this cell" rather than a wrong stake.

## Phases

1. **Capture + panel.** Click on Unabated shows the ticket with stake.
   Settings editable. Copy button.
2. **Line watcher.** Warning and re-sized stake when the captured line moves.
3. **README + CLAUDE.md line, pre-merge review, merge.**

## Testing

- `node --test unabated_ticket/tests/kelly.test.js` for `kelly.js`: the sheet's worked
  example (-400 vs fair -900, bankroll 30000, quarter Kelly) must give $3,750.
- Manual: load unpacked, open `tools.unabated.com/cfb/odds`, click a best-line
  price and a book-column price, a moneyline, a spread and a total; confirm
  side label, points sign, price, fair and stake; edit bankroll and see the
  stake move; wait for a line change and see the warning.

## Version control

- Branch `feature/unabated-ticket` on worktree `.worktrees/unabated-ticket`.
- Commits: (1) plan doc, (2) extension skeleton + kelly + tests, (3) capture +
  panel, (4) line watcher, (5) README + CLAUDE.md.
- Pre-merge review of `git diff main..HEAD`, then explicit approval, then
  merge to `main`, `git worktree remove .worktrees/unabated-ticket`,
  `git branch -d feature/unabated-ticket`.

## Documentation

- `unabated_ticket/README.md`: install (load unpacked), how capture works,
  settings, troubleshooting (panel empty, "could not read this cell").
- `CLAUDE.md` project structure: one bullet for `unabated_ticket/`.

## Risks

- Unabated front-end deploys can rename props or classes. Fallback path above;
  worst case is a small selector fix.
- The Unabated tab must stay open for the line watcher. If it is closed the
  panel shows the captured line with "not watching".
- Extension is unpacked, so Chrome shows a developer-mode banner once per
  launch. Acceptable for a personal tool.
