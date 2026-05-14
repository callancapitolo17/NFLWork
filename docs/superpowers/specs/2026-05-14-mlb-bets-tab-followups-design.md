# MLB Bets Tab — V8 Follow-ups (bug fixes + display/placement features)

**Date:** 2026-05-14
**Author:** brainstorming session
**Status:** spec for review
**Branch:** `worktree-mlb-bets-tab-followups`

---

## Review Pack

**What we're building.** Three bug fixes and three feature additions to the
MLB Dashboard bets tab (V8 card layout). Bugs: game times rendering in UTC
instead of viewer-local, Bet105 (and BFA) alt-total quotes not appearing in
the per-book grid, and the spread-opposite cell wrongly tagged as a
mismatched line. Features: a manual Kelly calculator widget below the
bankroll/Kelly settings, click-to-edit Risk on the bet card with
WZ-verified to-win on commit, and a per-cell devigged-odds toggle (RAW /
FAIR) on each card's grid.

**Key decisions** (genuine choice points):

1. **Timezone rendering — client-side `toLocaleString()`** (rejected:
   server-side `with_tz()`). Reason: matches the existing parlay-table
   pattern at `mlb_dashboard.R:515-524`, works for any viewer's tz without
   server config, no extra `pt_start_time` re-tagging gymnastics around
   the DuckDB round-trip that strips the `America/Los_Angeles` tag.

2. **Editable risk amount validation — let WZ be the source of truth**
   (rejected: client-side round-to-integer). Reason: WZ's
   `ConfirmWagerHelper` will reject invalid amounts with a specific
   `ErrorMsg` we can surface verbatim. Don't bake assumptions about WZ
   rules into the client; this also exposes the latent bug where the
   existing flow already passes raw float `bet_size` to WZ unchanged.

3. **WZ verified to-win — preflight on every risk commit** (rejected:
   trust local math). Reason: the same `ConfirmWagerHelper(RiskWin=2)`
   call we already use in `parlay_pricer.py` returns the *actual* `Win`
   value WZ would credit *and* the current odds — eliminates the
   `price_moved` failure class entirely. ~100 lines of code, half-day
   build, leverages existing pattern.

4. **Per-cell devig — toggle, default FAIR, American odds, probit math**
   (rejected: stacked %, inline %, RAW default, percentage display).
   Reason: matches user's manual workflow (devig sharps → check fair),
   keeps cell density unchanged when scanning, reuses the same
   `devig_american` (probit z-shift) function the model already uses for
   the `FAIR` chip — no methodology drift between the model's fair and
   the per-book fair.

**Risks / push back here:**

- **Scope** — six items in one spec. Recommended split into three PRs
  (bugs / display features / placement). Want different bundling?
- **Default view = FAIR** is a behavioral change. After this lands,
  every cell on every card defaults to showing devigged American odds
  instead of raw. Anyone glancing at the tab will need to learn that
  toggle exists. OK?
- **WZ verification latency** adds ~200–500ms between pressing Enter on
  the risk field and seeing "to win." We surface a "verifying…" spinner.
  Acceptable, or do you want optimistic local rendering with a quiet
  background reconciliation?
- **Off-line book devig** (e.g. B105 quoting line `0` when your bet is
  `-0.5`) is computed against the book's own line, not your bet's. We
  show it with the existing amber background + line tag — but the
  number genuinely answers a different question than the cell next to
  it. You said this was OK; flagging in case you change your mind once
  you see it on real data.
- **Bookmaker has no alt markets in its raw scraper output** (verified
  in `bookmaker_odds/bookmaker.duckdb` — only base spreads, no
  `alternate_*`). Even after fixing the alt-totals join, Bookmaker pills
  won't appear on alt-line cards. That's a scraper coverage gap, not a
  bug in this work. Worth a follow-up to extend the scraper.

**Worth understanding** (R analogues, opt-in):

- **Probit devig (z-shift).** Like normalizing a vector that should sum
  to 1 — but instead of dividing by the total (`p / sum(p)`, the
  multiplicative method), you convert each prob to its z-score with R's
  `qnorm()`, shift all z's by a constant `c` until the resulting probs
  sum to 100%, then convert back with `pnorm()`. For 2 sides the math
  is closed-form: `c = -(z1 + z2) / 2`. The whole point of the shift is
  that it preserves the *relative information* the bookmaker's line
  encoded — which is more accurate than naively pro-rating, especially
  for lopsided lines. The function lives in `Tools.R::devig_american`.
- **DOM `contentEditable`.** The browser's "make this text typeable on
  click" attribute. We use it for click-to-edit Risk because we want
  the value to *read as plain text* until interacted with — a regular
  `<input>` would always look like a form field. R analogue: think of
  it as the inverse of `shiny::textOutput`; instead of "always-text,"
  it's "text-by-default, input-on-click." The trade-off is you lose the
  built-in keyboard handling that `<input type=number>` gives you (we
  add `keydown` handlers for Enter/Escape).
- **Preflight pattern (ConfirmWagerHelper).** Call an endpoint that
  *would* place the bet, with a flag (`RiskWin=2`) that says "don't
  actually do it; tell me what *would* happen." Returns the price,
  balance check, and any error WZ would raise. R analogue: `dryrun=TRUE`
  in many DB wrappers, or `dbplyr::show_query()` returning the SQL that
  would execute. We already use this pattern for parlay pricing
  (`parlay_pricer.py`); this spec extends it to singles.

---

## Scope

This spec covers the MLB Dashboard bets tab (port 8083, V8 card layout
landed 2026-05-13). Six items:

| # | Item | Type | Files touched | Suggested PR |
|---|---|---|---|---|
| 1 | Game times in viewer-local tz | bug | `mlb_dashboard.R` | A |
| 2 | Bet105/BFA alt-total quotes appearing on cards | bug | `Answer Keys/MLB Answer Key/odds_screen.R` | A |
| 3 | Spread-opposite cell wrongly tagged as alt-line | bug | `Answer Keys/MLB Answer Key/odds_screen.R` | A |
| 4 | Kelly Calculator widget | feature | `mlb_dashboard.R` | B |
| 5 | Per-cell devig toggle (RAW / FAIR) | feature | `mlb_dashboard.R`, `book_cell.R`, `odds_screen.R` | B |
| 6 | Editable Risk + WZ verified to-win | feature | `mlb_dashboard.R`, `mlb_dashboard_server.py`, `wagerzon_odds/single_pricer.py` (new) | C |

**Out of scope.** Bookmaker scraper extension to cover alt markets;
parlay-tab parity for the same features (parlays already have
ConfirmWagerHelper preflight via `parlay_pricer.py`); the same
treatment for the trifecta tab; CBB dashboard parity.

---

## Bug fixes (PR A)

### Bug 1 — game times in viewer-local tz

**Symptom.** The hero-strip subtitle on each V8 card shows tipoff in UTC
clock time. A 4:06 PM PT game reads as `Thu 11:06 PM`.

**Root cause.** `mlb_dashboard.R:1425-1426` formats `pt_start_time`
server-side:
```r
tipoff = ifelse(is.na(pt_start_time), "",
                format(pt_start_time, "%a %I:%M %p"))
```
`pt_start_time` is written by `format_bets_table()` in `Tools.R:1148` as
`with_tz(..., "America/Los_Angeles")`, but it round-trips through DuckDB
(`mlb_bets_combined`) which stores `TIMESTAMP` without a tz tag. When R
reads it back via `dbGetQuery`, it lands as UTC POSIXct. `format()` then
prints UTC clock time.

**Fix.** Match the legacy parlay-table approach
(`mlb_dashboard.R:515-524`):
1. In `create_bets_table()`, emit the time as a UTC ISO string with `Z`
   suffix into a `data-tipoff` attribute on each card, *not* as a
   pre-formatted `tipoff` string.
2. Add a small `<script>` block to the page that walks every
   `[data-tipoff]` element and replaces its visible text with
   `new Date(t).toLocaleString(undefined, {weekday:'short', hour:'numeric', minute:'2-digit'})`.
3. The browser converts to whatever the viewer's OS timezone is — works
   for the user (PT) and any other viewer.

**Files.** `Answer Keys/MLB Dashboard/mlb_dashboard.R` only. No server,
no R-side tz logic.

### Bug 2 — Bet105/BFA alt-total quotes don't appear on cards

**Symptom.** A bet flagged as `alternate_totals_fg` (e.g. Royals Under
6.5) shows the WZ pill on the card, but Bet105 (and BFA) show no pill,
even though both books have `alternate_totals` rows for the game in
their raw scraper output.

**Root cause.** Market-name mismatch in the join performed by
`expand_bets_to_book_prices()` in `Answer Keys/MLB Answer Key/odds_screen.R`.

- `mlb_bets_combined.market` for the bet is `"alternate_totals_fg"`
  (the suffix convention used by `compare_alts_to_samples`).
- The Bet105 raw market is `"alternate_totals"` (no `_fg`). After
  `map_scraper_markets_mlb` and `scraper_to_canonical`, it lands as
  `market = "alternate_totals"`, `period = "FG"`.
- The join condition is `market == bet$market_type`. After
  `.derive_market_type()` (which today only strips `_1st_X_innings`),
  `bet$market_type` is still `"alternate_totals_fg"` (the `_fg` suffix
  isn't recognized). So `"alternate_totals" != "alternate_totals_fg"` →
  no match → no row emitted → no pill rendered.

**Why WZ works.** WZ writes `alternate_totals_fg` literally in its raw
DB, so its book frame matches the bet by accident. Bet105 and BFA do not.

**Fix.** Extend `.derive_period()` and `.derive_market_type()` in
`odds_screen.R` to also recognize the `_fg / _f3 / _f5 / _f7` suffix
convention, in addition to the existing `_1st_X_innings` recognition.
After the change, both `alternate_totals_fg` and `alternate_totals`
canonicalize to `(market = "alternate_totals", period = "FG")`. The
join then matches correctly. Same fix unlocks `alternate_spreads` from
Bet105 and any other book that uses the un-suffixed alt naming.

**Files.** `Answer Keys/MLB Answer Key/odds_screen.R` only.

**Verification.** Re-run the MLB pipeline; query `mlb_bets_book_prices`
should now contain rows where `bookmaker = 'bet105' AND market =
'alternate_totals'`. Reload the dashboard; Bet105 alt-totals pills
should render alongside WZ on alt cards.

### Bug 3 — spread-opposite cell wrongly tagged as alt-line

**Symptom.** On a spread bet (e.g. ATH -0.5 / STL +0.5), every cell on
the *opposite* row (STL +0.5) renders with the amber `+0.5` line tag,
implying every book is on a different line — even when they're on the
same line.

**Root cause.** In `expand_bets_to_book_prices()`,
`.pick_closest_line(candidates, bet$line, bet$bet_on)` is called with
`bet$line` for both the pick and opposite slots. But in
`scraper_to_canonical()`, the opposite (away) side of a spread is
emitted with the *negated* line (so a book offering ATH -0.5 / STL +0.5
ends up with `line = -0.5` for ATH and `line = +0.5` for STL).

So when matching the opposite slot:
```
.dist          = abs(+0.5 − (−0.5)) = 1.0
is_exact_line  = (1.0 < 1e-9)        = FALSE  ← bug
```
`render_book_cell` then renders the amber tag on every opposite cell,
even when the book is on the same `0.5` line you bet.

**Why FanDuel slips through.** On the example card, FD doesn't post F5
spreads — it's matched on F5 ML where `line = NA` on both sides. The
NA-NA branch in `.pick_closest_line` returns `is_exact_line = TRUE`, so
FD's cell appears tag-less by accident, masking the bug for one book.

**Fix.** Inside the slot loop in `expand_bets_to_book_prices()`,
compute an `effective_line` per slot before calling
`.pick_closest_line`:
```r
is_spread_market <- bet$market_type %in% c("spreads", "alternate_spreads")
effective_line <- if (slot == "opposite" && is_spread_market && !is.na(bet$line)) {
  -bet$line
} else {
  bet$line
}
chosen <- .pick_closest_line(candidates, effective_line, bet$bet_on)
```
For totals (Over/Under share the same line number) and moneyline
(line = NA), behavior is unchanged.

**Files.** `Answer Keys/MLB Answer Key/odds_screen.R` only.

---

## Display features (PR B)

### Feature 4 — Kelly Calculator widget

**Goal.** A standalone manual Kelly tool living in the bets tab. Type
American odds + a fair %, get a recommended risk size based on the
existing dashboard `Bankroll` × `Kelly Fraction`. Read-only — does not
mutate any bet card.

**Placement.** Own strip directly below the existing Bankroll/Kelly
settings, above the filter bar. Always visible while scanning bets.

**UI shape.**
```
┌─ Kelly Calculator ──────────────────────────────────────────────┐
│  Odds: [ +120  ]   Fair: [ 52.5% ]   →   Risk: $76              │
│  edge +7.0%   ·   Kelly 3.05%   ·   breakeven 45.5%             │
└─────────────────────────────────────────────────────────────────┘
```

- Two text inputs: `Odds` (American: `+120`, `-110`), `Fair`
  (accepts `52.5%`, `0.525`, *or* American fair odds `-110` /
  `+120` — auto-detect by sign and presence of `%`).
- Live output: `Risk` (computed) — prominent green; `to win`
  (computed); detail row showing edge / Kelly% / breakeven.
- Reads the existing dashboard `Bankroll` + `Kelly Fraction` inputs
  live (no separate inputs to keep in sync).
- Negative-EV behavior: if fair < breakeven, Risk → `$0` (greyed),
  edge chip → red.
- Invalid input (non-numeric Odds / Fair) → red border on the
  offending field; outputs blank.

**Math.** Reuse the existing `calculateKellyBet(prob, americanOdds,
bankroll, kellyMult)` JS function at `mlb_dashboard.R:3746`. Kelly% =
`max(0, (b·p − q) / b)` with `b = decimal_odds − 1`. Size = `min(bankroll,
bankroll × kellyMult × kelly%)`.

**Files.** `Answer Keys/MLB Dashboard/mlb_dashboard.R` only — new HTML
block in the bets-tab template, ~50 lines of inline JS reusing
`calculateKellyBet` and the bankroll/kelly-mult inputs.

**Visualization.** Mockup at
`/Users/callancapitolo/.claude-personal/jobs/fafc045f/kelly_widget_mockup.html`
served at http://localhost:8765/kelly_widget_mockup.html during the
brainstorm session.

### Feature 5 — Per-cell devig toggle (RAW / FAIR)

**Goal.** Each bet card's per-book grid gains a small `RAW` / `FAIR`
toggle in the corner. **Default is FAIR** — every cell shows the
book's two-sided quote devigged with probit math, rendered as American
odds. Click `RAW` to flip back to the original book quotes.

**Devig math.** Probit (additive z-shift), 2-way closed-form
`c = -(z1 + z2) / 2`. Same function as `Tools.R::devig_american`.
Implemented in JS for client-side rendering — math copied from the R
function, ~30 lines.

**Per-cell behavior.**
- For each cell, devig is computed from the book's *own* two-sided
  quote on the *book's own* line (not the bet's line). So WZ's
  `+120 / -140` (sum 103.8% vig) → `+129 / -129`.
- Off-line cells (book on a different line than the bet) still get
  the existing amber background + tiny line tag (`0`, `8.5`, etc.) —
  same visual cue used today on RAW view. The devig math runs on the
  book's line, the line tag tells you it's not your bet's line.
- Empty cells (no quote) show `—` in both views.
- Pick book cell keeps the green-bordered pick treatment; in FAIR
  view the displayed odds are the pick book's devigged price (so the
  "pick" is also visible as fair).

**Toggle UI.** Small two-button segmented control (`RAW` | `FAIR`)
in the top-right of each card's grid container. Toggle state is
*per-card* (independent for each card on the page) — no global
preference. Toggle defaults to `FAIR` on every card every render.

**Files.**
- `Answer Keys/MLB Dashboard/mlb_dashboard.R` — toggle markup, JS
  flip handler, JS probit-devig helper.
- `Answer Keys/MLB Dashboard/book_cell.R` — `render_book_cell()`
  emits both `<span class="raw">` and `<span class="fair">` per
  cell; CSS hides one or the other based on a class on the grid
  container.
- `Answer Keys/MLB Answer Key/odds_screen.R` — `expand_bets_to_book_prices()`
  needs to surface the *opposite-side American odds* per cell so the
  client can compute the devig pair. Today only one side's odds are
  in the wide row per book; we need both. Add
  `<book>_opposite_american_odds` columns alongside the existing
  `<book>_american_odds`.

**Visualization.** v2 mockup at
http://localhost:60167 (brainstorm session
`/Users/callancapitolo/.claude-personal/jobs/fafc045f/.superpowers/brainstorm/83700-1778783805/`)
shows the FAIR-default toggle, real probit-devigged values, and
off-line cell treatment.

---

## Placement feature (PR C)

### Feature 6 — Editable Risk + WZ verified to-win

**Goal.** The hero-strip `Risk $XX` value on each card becomes
click-to-edit. Type any number, press Enter / blur, and the dashboard
fires a WZ preflight quote that returns the *actual* to-win and the
*actual* current odds — surfacing line drift before you click Place
Bet. Eliminates the `price_moved` failure class for WZ-targeted bets.

**UX.**
1. Hover `$31` Risk value → dashed outline.
2. Click → becomes a green-bordered editable text field with current
   value selected.
3. Type a new amount (decimals allowed; no client-side rounding).
4. Press Enter / blur → commit:
   - Override stored locally; value renders amber + small `↶` snap-back
     button appears.
   - "To Win $XX" is replaced by a small spinner.
   - Fire `POST /api/wz-quote-single` with `{bet_hash, amount, account}`.
   - On response (~200-500ms): swap spinner for verified to-win.
     - If WZ accepts at displayed odds: `To Win $30` (verified ✓).
     - If odds drifted: surface inline `WZ now -125 (was +120) → $20`.
     - If WZ rejects: surface `ErrorMsg` verbatim near the field.
5. `↶` button → snap back to model's recommended size; clears verified
   state.
6. Escape mid-edit → cancel commit, no server call.
7. Page refresh → override lost (no persistence; user accepted in
   brainstorm).

**Verification badge.** WZ rows show a small "verified" indicator on
the To Win value (subtle green check or `✓ wz`) once a quote completes.
Non-WZ rows still show local math (no API to verify against).

**Server endpoint.** New `POST /api/wz-quote-single` in
`Answer Keys/MLB Dashboard/mlb_dashboard_server.py`:
```
Body:    {bet_hash, amount, account}
Returns: {win, current_wz_odds, error_msg, error_msg_key}
```
Looks up the bet in `mlb_bets_combined` to get `idgm`, `play`, `points`,
`american_odds`. Calls `wagerzon_odds.single_pricer.get_single_price(session, bet, amount)`
(new module, mirrors `parlay_pricer.get_parlay_price`). Returns the
parsed response. ~50 lines.

**New module.** `wagerzon_odds/single_pricer.py`. Mirror of
`parlay_pricer.py::get_parlay_price()` but with `WT="0"` (single, not
parlay) and a single-leg `sel` string. Same `RiskWin="2"` flag to skip
balance validation. ~30 lines, structured to match the parlay file's
conventions for future maintainability.

**Existing-flow fix bundled here.** Today, `mlb_dashboard.R:3803-3804`
sends both `actual_size` and `kelly_bet` from the same `data-size`
attribute. With editable risk this would lose the model's recommended
number. Split into two attributes:
- `data-size` — the editable amount, drives `actual_size`
- `data-model-size` — immutable Kelly result, drives `kelly_bet` /
  `recommended_size`

This preserves CLV-grade data integrity: future analysis can ask "when
did I override the model and how did I do?"

**Files.**
- `wagerzon_odds/single_pricer.py` (new) — the preflight helper.
- `Answer Keys/MLB Dashboard/mlb_dashboard_server.py` — new endpoint.
- `Answer Keys/MLB Dashboard/mlb_dashboard.R` — click-to-edit JS,
  snap-back button, verified-to-win swap, two-attribute split for
  `data-size` / `data-model-size`.

---

## Implementation order

Three sequential PRs, each on its own branch off main.

- **PR A — Bugs (1, 2, 3).** Smallest changes, all visible, all
  independent. Each fix isolated to one file. Bundled because they all
  ship the same week and reviewing 3 small fixes in one PR is faster
  than three separate PRs.
- **PR B — Display features (4, 5).** Kelly calculator + per-cell
  devig toggle. Both display-only (no server / no place-flow). The
  probit devig math used in feature 5 informs no other PR; the Kelly
  math reused from existing code. Independent of PR A but easier to
  review after the bugs land so the screenshots match expectations.
- **PR C — Editable risk + WZ verified to-win (6).** Larger surface
  area: new module, new server endpoint, click-to-edit JS, two-attr
  split in placement plumbing. Earns its own PR for review attention.

Each PR gets its own implementation plan written via the
`superpowers:writing-plans` skill, including version-control,
worktree, testing, and documentation sections per project policy
(`CLAUDE.md`).

## Documentation updates required

Same merge as the code:

- **`Answer Keys/CLAUDE.md`** — Add an entry to the "MLB Dashboard —
  Odds screen" section describing the per-cell devig toggle (default
  FAIR, probit math, file touched). Update the "Pitfalls" section to
  note the `_fg` suffix canonicalization in `odds_screen.R`.
- **`CLAUDE.md` (root)** — No update needed; this is internal to the
  MLB Dashboard.
- **`Answer Keys/MLB Dashboard/README.md`** — Brief mention of Kelly
  Calculator widget in the feature list, and the WZ verified-quote
  flow under "placing bets."
- **`wagerzon_odds/CLAUDE.md`** — Add `single_pricer.py` to the
  "Quick map" list, mirroring the `parlay_pricer.py` entry.

## Testing requirements

- **PR A bugs:** Re-run pipeline and load the dashboard. Visual
  verification of (1) tipoff in PT clock, (2) Bet105 alt-total pills
  appearing on a known alt-total card, (3) opposite-side cells on a
  spread bet rendering without amber line tags when the book is on
  the same line.
- **PR B display features:** Kelly calc — manual cases for known
  inputs (e.g. `+100` at `55%` should output `Kelly 10%, edge +5%`).
  Devig toggle — manual cases against `Tools.R::devig_american` for a
  few known book quotes; assert JS and R produce the same American
  fair odds within 0.5 cents. R-side existing tests in
  `Answer Keys/MLB Answer Key/odds_screen.R` tests for the new
  `<book>_opposite_american_odds` column.
- **PR C placement:** New unit test for `single_pricer.get_single_price`
  with mocked WZ session response (mirror `test_parlay_placer.py`
  pattern). Manual end-to-end: click into Risk on a real WZ-pickable
  bet, type `1`, press Enter, confirm spinner → `$X to win` swap, then
  do not actually place. Then place a $1 bet end-to-end to confirm
  the existing placement flow still works after the `data-size` /
  `data-model-size` split.

## Open questions (none expected — all resolved in brainstorming)

- ~~Timezone fix approach~~ → JS `toLocaleString()`
- ~~Editable risk validation~~ → trust WZ
- ~~Verified quote scope~~ → bundle into editable-risk PR
- ~~Devig display variant~~ → toggle (C)
- ~~Devig default view~~ → FAIR
- ~~Devig number format~~ → American odds, not %
- ~~Off-line cell treatment~~ → amber background + line tag (existing pattern)
- ~~Probit vs multiplicative~~ → probit (matches existing codebase)
