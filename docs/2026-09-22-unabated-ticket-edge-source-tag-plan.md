# Unabated Ticket — why an edge grew: fair moved to you vs book moved away (plan, issue #132)

Spec: GitHub issue #132 and its 2026-09-22 comment (four cases, the fair
decides). Context: #130 (conditional Kelly, merged 2026-09-19 as 8b04947 —
the `add $X` the tag qualifies) and #126 (the unread `openerPrice` /
`openerPoints` fields). The stake is never changed by this work: the tag
informs, the number to act on stays what #130 computes.

## 1. Goal

An edge on a held line grows and the row says `add $X`. Two causes, opposite
actions:

| Fair (`bacr`, in probability) | Book price on the side | Tag | Meaning |
|---|---|---|---|
| moved toward the side | anything | `fair moved to you` (green) | the sharps agree; the price is a bonus |
| unchanged | got better | `book moved away` (amber) | the book is shading against the side and the fair has not answered yet (Unabated's fair is ~1–2 min behind: v2 rebuild ~1/min + ~40 s ingest, #126); wait one snapshot — if the fair holds and the price is still there it is a stale soft line |
| moved against the side | anything | `fair moved against you` (red) | the market is moving against the side and the book is ahead of the fair; the edge can still read bigger because the price improved by more than the fair fell |
| unchanged | unchanged | no tag | |

The red row is the adverse-selection case the issue exists for: adding on
every "improvement" puts the most money on right as the market turns. The
panel keeps no previous value today (`feed.applyChanges` and the per-league
snapshot merge replace a line in place), so it cannot tell the cases apart.
This adds a small per-line history in memory and one tag next to every edge
figure.

## 2. Math — new `extension/edgemove.js`

Pure, IIFE + `module.exports` like `kelly.js` (`globalThis.UnabatedEdgeMove`).
Requires `kelly.js` for `americanToProb` / `bookProbOf`.

```
observe(history, line, { at, source })     source: "snapshot" | "stream"
forget(history, key)
edgeMove(entries, now) -> { kind, fairDelta, priceDelta, sinceMs, source, from, to }
```

- **History** = `history[line.key]` → `[{at, source, points, price,
  sourceFormat, sourcePrice, bacr}]`, newest last. `observe` appends only
  when the line is new or `price` / `sourcePrice` / `bacr` differ from the
  last entry, then prunes: entries older than `at − WINDOW` are dropped
  except the newest of them (the baseline the window compares against).
- **A number move resets the history** (not in the issue; found today: 1,466
  of 2,720 live NFL spread/total lines sit on a different number than they
  opened). A price at 48.5 is not comparable to one at 47.5, so an
  observation whose `points` differ from the last entry starts the key
  afresh and the line reads as first seen. Alt rungs are keyed by points and
  never change number.
- **Window** `EDGE_MOVE_WINDOW_MS = 10 min`: covers the fair lag plus two
  snapshots on the slowest tier (5 min), so "wait a snapshot" can be
  followed while the tag is still up. **Threshold**
  `MOVE_THRESHOLD = 0.005` (0.5 points of probability, the issue's start)
  with a `1e-9` epsilon so a Novig half-cent (`0.565 − 0.56` is
  `0.00499999…` in floats) counts as moved.
- **Reference** = the newest entry at or before `now − WINDOW`, else the
  oldest entry (a line first seen inside the window compares against its
  first sighting). Reference = the current entry means nothing moved in the
  window → `kind: "none"`. One entry (first sighting, or unchanged since a
  number move) → `none`.
- **Deltas in probability**, signed so positive is edge-increasing for the
  row's side: `fairDelta = prob(to.bacr) − prob(from.bacr)` (null when either
  fair is missing); `priceDelta = bookProb(from) − bookProb(to)` (the book
  lengthened the side's odds), where `bookProb` is `kelly.bookProbOf` so an
  exchange's exact `sourcePrice` is compared, not its whole-American rounding.
- **Kind — the fair decides** (issue comment 2026-09-22): `fair_to_you`
  when `fairDelta ≥ T`, `fair_against` when `fairDelta ≤ −T`, both whatever
  the price did; a fair inside the threshold (or unknown at either end) is
  `book_away` when `priceDelta ≥ T`, else `none`. A book that only
  shortened never earns a tag. So the amber tag turns red when the fair
  follows the book (the collapse the issue describes), with
  `fair 35.6% → 33.7%` in the tooltip.
- **`sinceMs` / `source`**: the first entry after the reference — when the
  move was first observed and whether that was a snapshot or a stream tick.
  For Kalshi/Novig and every alt rung that is always a snapshot, up to one
  refresh interval after the book moved.
- Bad input (no price, non-array entries) throws with expected vs found.

## 3. Feed — `extension/feed.js`

- `applyChanges` also returns `appliedKeys` (the keys it replaced or added)
  so the scanner can record stream observations without re-diffing 100k lines.
- `normalizeSnapshotLine` reads `openerPrice` and `openerPoints` (#126).
  Measured live 2026-09-22 on the NFL file: every one of 2,720 spread/total
  and 1,393 moneyline lines carries an opener, and it is **per book** (on
  all 136 sides priced by ≥3 books the openers differed across books); alt
  rungs carry none (0 of 40,499). The changes stream has no opener, so
  `applyChanges` carries both over from the held line, as it does
  `liquidity`. `describeLine` exposes them.
- The opener is tooltip context only — a longer window for free — never a
  tag: its fair is not in the feed, so it can only say the book moved since
  open. Shown as `opened -120` when the number is unchanged, `opened -120 at
  -3` when it is not (the price is not comparable then), nothing when absent.

## 4. Scanner — `extension/scanner.js`

- Owns `history` (reset by `start()`), records with `now()` as `at`:
  every line of a league's snapshot after `mergeInto` (`source:
  "snapshot"`), every `appliedKeys` line after `applyChanges` (`"stream"`).
  `mergeInto` reports the keys it dropped; those not re-listed are forgotten
  so a pulled line does not keep its history for the session.
- Handed to the panel as the third argument of `onChange(status, state,
  history)` and via `getHistory()`. In memory only, never storage.

## 5. Panel — `panel.js`, `panel.html`, `panel.css`, `manifest.json`

- `moveTag(row)` → `<span class="tag move-fair|move-book|move-against">`
  reading `fair moved to you` / `book moved away` / `fair moved against
  you`, with the numbers in `title`:
  `fair 33.7% → 35.6% · price +199 → +215 · moved 2m ago (snapshot) ·
  opened +185`. Nothing for `none`.
- Rows and cards (`rowParts`): the tag sits in the rail under the edge
  figure, with one small line under it naming the mover (`fair 33.7% →
  35.6%` for green / red, `price +125 → +141` for amber; user choice after
  the mock, 2026-09-22). Card expander lines (`renderGroupLine`): same tag next to the
  edge, each book moves on its own. Ticket (`renderTicket`): the tag next to
  the edge fact, from the feed's copy of the ticket's line only when its
  price equals the price being sized (the same one-price rule `pricedLine`
  uses).
- Alerts (`notifyEdge`): the tag's words go into the notification body
  after the edge (`+3.2% edge · book moved away · stake $120 · …`), so a
  re-fire on an improved edge (#112) says which improvement it was.
- CSS: `.tag.move-fair` green (the `held` treatment), `.tag.move-book`
  amber (`--warn-bg` / `--warn-fg`), `.tag.move-against` red (the `against`
  treatment); rail alignment. `panel.html` loads
  `edgemove.js` before `panel.js`. Manifest 0.9.0 → 0.10.0.

## 6. Tests

`tests/edgemove.test.js` (new): fair up / price flat → `fair_to_you`; price
better / fair flat → `book_away`; fair up + price better and fair up + price
worse → still `fair_to_you`; fair against + price better (the red row, with
the edge bigger) and fair against + price flat → `fair_against`; neither
(sub-threshold) → `none`; a book that only shortened → `none`; first
sighting → `none`; an alt rung with snapshot-only history reports `source:
"snapshot"`; a whole-point fair change (-150 → -160, 1.5 pts) vs a one-point
one (-150 → -151, 0.16 pts); a Novig half-cent counts; a number move resets;
a move older than the window → `none`; a missing fair decides on price
alone; pruning keeps one baseline.

`tests/scanner.test.js`: every line has one snapshot entry after `start`; the
fixture's spread move (-14 → -3.5, a number move) resets to one stream entry;
a scripted stream tick with a longer price on the same number → `book_away`
from `stream`; a re-downloaded snapshot with a moved `bacr` → `fair_to_you`
from `snapshot`, and the same for an alt rung's price; a line the new
snapshot no longer lists loses its history; `start` clears it.

`tests/feed.test.js`: `appliedKeys`; openers parsed on main lines, null on
alts, carried over on a stream replace.

## 7. Version control

- Branch `feature/unabated-ticket-edge-source-tag` (this worktree's branch,
  renamed from `claude/vibrant-bartik-eedc58` before this file was written).
- Commits: (1) this plan; (2) `edgemove.js` + tests, `feed.js`
  (`appliedKeys`, openers) + tests; (3) scanner history + tests; (4)
  `panel.js` / `panel.html` / `panel.css` / manifest; (5) README + root
  `CLAUDE.md`.
- New: `extension/edgemove.js`, `tests/edgemove.test.js`, this plan.
  Modified: `feed.js`, `scanner.js`, `panel.js`, `panel.html`, `panel.css`,
  `manifest.json`, `eslint.config.js` (module list comment only),
  `tests/feed.test.js`, `tests/scanner.test.js`, `unabated_ticket/README.md`,
  root `CLAUDE.md`.
- `./unabated_ticket/check.sh` green on the branch, then the pre-merge
  review, then ask. No merge, no push without an explicit yes.

## 8. Worktree

Work happens in `.claude/worktrees/vibrant-bartik-eedc58`. No DuckDB file is
touched or copied. After an approved merge: `git worktree remove` the path
and `git branch -d feature/unabated-ticket-edge-source-tag`.

## 9. Documentation

Same branch, after the code is final: `unabated_ticket/README.md` — a new
**Edges → Why an edge grew** subsection (the two cases and why they matter,
the window and threshold, what "ago" means per source, number moves, the
opener line, "the stake never changes for it"), one line in **Stake** after
the #130 paragraph linking there before an `add`, the alert body in
**Alerts**, test counts and the new suite in **Tests**, a dated entry in the
design decisions log. Root `CLAUDE.md`: one clause in the Unabated Ticket
blurb.

## 10. Known limits

- `bacr` is a whole American price: near even money the fair must move
  about three American points to clear 0.5 probability points (-110 → -112
  is 0.45 pts, unchanged; -110 → -113 is 0.67, moved). At +200 a five-point
  move clears it.
- History lives only while the panel is open and starts empty on every open
  and every league change; the first refresh after opening shows no tags.
- A number move within the window reads as first seen, not as a move.
- Verified by the node suite and fixtures only; no live panel session in
  this worktree. A hold rule ("no add for one snapshot after an amber, never
  on a red") is a separate decision once the tag has been watched on live
  cards.
