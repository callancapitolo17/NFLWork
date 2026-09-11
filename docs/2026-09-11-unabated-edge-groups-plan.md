# Unabated Ticket — Edges grouped by market (plan)

Follows #113 (alt lines). With alts on, one soft ladder produces 4–10 rows
that all say one thing ("Novig's Idaho spread ladder is soft vs consensus").
The Edges tab gains a **Group by market** view: one card per
(game, period, bet type, side) showing the best line, with the rest behind
an expander.

## Design

- **Unit** = (eventId, periodTypeId, betTypeId, sideIndex). A +EV opinion is
  directional, so the two sides of a market are two cards.
- **Best line** = highest Kelly stake (stake = edge / (decimal − 1) already
  taxes longshots), edge as the tie-break. The card prints the best line's
  book, rung, price, edge and stake, plus `N books · M lines`; clicking the
  chevron reveals the other lines, each clickable like today. Clicking any
  line runs the existing locate.
- **Sort** (edge / stake / start) applies to cards through their best line.
- **Alerts** key on the card: one notification per (game, market, side) the
  first time its best line crosses the threshold, again only when the best
  edge improves; the per-event 5-min cooldown is unchanged.
- **Ungrouped** stays available (toggle off), byte-for-byte the #113 list.
- Pure grouping in `feed.js` (`groupEdges(rows, rankOf)`), node-tested;
  the panel only renders.

## Version control

Branch `feature/unabated-edge-groups` off `feature/unabated-alt-edges`,
worktree `.worktrees/unabated-edge-groups`. Merge after #113 (or together).
Commits: (1) `feed.groupEdges` + tests, (2) panel view + alerts + harness,
(3) docs. Pre-merge review, explicit approval, worktree + branch removed.

## Files

`unabated_ticket/extension/{feed.js, panel.js, panel.html, panel.css,
manifest.json}`, `unabated_ticket/tests/feed.test.js`,
`unabated_ticket/README.md`, `CLAUDE.md` (Unabated Ticket bullet).
