# Issue #86 — KXMLBF5 (F5 winner) tie-adjusted pricing: implementation plan

**Date:** 2026-08-13 · **Branch:** `worktree-issue-86-f5-winner-tie` · **Epic:** #82 · Depends: #84, #85 (merged)

Recon findings (posted on the issue): Kalshi F5-winner team markets are
**unconditional, ties lose** ("all team strikes will resolve to No"); every
book that prices F5 ML quotes a **conditional push-2-way** number
(implied sums 1.01–1.07); FanDuel and BetMGM both carry a **3-way F5
result market** in the same structures the on-demand flight already
fetches — a free P(tie) anchor; the anchors agree with each other and
with Kalshi's own TIE market within ~0.8pt.

## Core math (exact, not approximate)

"Team leads after 5" is a subset of "not tied", so for ANY combo C
containing an F5-winner leg:

    P(C) = P(C | not tied) × (1 − P(tie))

with the *marginal* P(tie) — no correlation assumption, even when C also
contains F5 totals. The books' SGP prices for such combos live entirely
in conditional space (a pushed leg refunds), so one multiplication at
the end converts.

## Design

1. **Conversion function** — `kalshi_common/fair_value.py::f5_winner_fair_from_book(conditional_fair, p_tie, semantics)`.
   Pure; semantics `"push_2way"` → `fair × (1 − p_tie)`;
   `"three_way_unconditional"` → identity (future books); anything else /
   invalid p_tie → `None` (fail closed). Unit-tested against hand cases.
2. **Semantics label + tie anchor live in the book adapters** (per issue
   spec). Each on-demand hook dict gains `f5_ml_semantics` (constant:
   `"push_2way"` at DK/FD/MGM/NV/CZR; PX has no F5 board) and a
   `tie_prob(structure)` hook: FD deviggs its "First 5 Innings Result"
   runners, MGM its "1st 5 innings - 1X2" options (both already inside
   the fetched structures — zero extra wire calls); DK/NV/CZR/PX → None.
3. **`OnDemandBookResult`** gains `f5_semantics` and `p_tie_book`
   (default None, additive) — set by `_price_on_demand` when the leg set
   contains an F5-winner team leg.
4. **Conversion happens in `OnDemandEngine.lookup`** — the single choke
   point every consumer (quote tick, confirm last-look, book-drift
   breaker) already reads. Per flight: if the leg set has an F5-winner
   leg, `p_tie` = median of landed `p_tie_book` anchors; every book fair
   converts via the shared function; **no anchor landed → `lookup`
   returns `{}`** (declines as `too_few_books`; firehose shows why).
   Raw results in `lookup_results` stay unconverted (research reads).
   This runs BEFORE `router.consensus`, satisfying the σ_z requirement.
5. **Guard lifts** (all three, coordinated):
   - `legset.classify_subcombo`: F5-winner **team** legs route
     `on_demand`; **TIE-side legs stay unpriceable** (no book can price
     them). home+away same game already dies on the dup guard
     (`("ml","F5",None)` collision).
   - `sgp_service.price_on_demand` period backstop: allow (ml, F5) team
     sides through; keep an explicit zero-wire decline for
     `side in ("tie","not_tie")`.
   - `main.py` scope label: F5-winner team legs become in-scope;
     TIE-side legs keep a distinct `out_of_scope_f5_tie_leg` label.
6. **Book-side name fixes** so legs actually resolve:
   - DK: `fetch_main_market_nums` also accepts bare `"1st 5 Innings"`
     as the F5 moneyline market (live board name; fails closed if its
     selections aren't 0ML-form).
   - CZR: classify `"1st 5 innings money line"` / `"1st 5 innings -
     money line"` → ("F5","moneyline") (name inferred from CZR's posted
     1st-3/1st-5 family; fails closed if wrong).
   - FD/MGM/NV already resolve (ml, F5); PX declines cleanly.
7. **Firehose provenance**: `on_demand_result` per-book entries gain
   `semantics` + `p_tie_book`; `_live_games_detail` (rides
   `quote_priced`) gains per-game `p_tie_used` + converted fairs.
8. **Not in v1**: Kalshi's own TIE market as a P(tie) fallback anchor
   (extra GET + orderbook guards). Book anchors are free and present
   whenever FD/MGM post F5 boards; if firehose data later shows
   anchor-missing declines, add it as a follow-up. The #23 Fréchet gate
   already cross-checks converted fairs against Kalshi's unconditional
   team marginal at quote time.

## Tests (fail-first)

- `fair_value`: even-odds push-2way with p_tie 0.13 → 0.435; the live FD
  case 0.6236/0.1423 → 0.5349; identity for 3-way; None for bad
  inputs/unknown semantics.
- `legset`: F5-winner 2-leg combo → `on_demand` (fails today:
  `unpriceable`); TIE-leg combos stay `unpriceable`; home+away dup.
- `sgp_service`: (ml,F5) team legs reach hooks (fails today); tie-side
  backstop; FD/MGM `tie_prob` devig from structure fixtures; result
  fields ride back.
- `mlb_sgp`: FD Result-market capture, MGM 1X2 parse, CZR/DK name
  classification.
- `kalshi_mlb_mm`: engine lookup conversion (median anchor, no-anchor →
  `{}`); mixed-semantics dispersion test (converted fairs pass the gate,
  unconverted would inflate σ_z); scope labels; research payload fields.

## Workflow

- Worktree `issue-86-f5-winner-tie` (already on it), no writes outside.
- Commits: (1) conversion math + fail-first tests; (2) book adapters
  (FD/MGM tie3, DK/CZR names); (3) guard lifts + engine conversion +
  firehose; (4) docs (README + CLAUDE.md blurb).
- Full suite `python -m pytest kalshi_mlb_mm/tests/ kalshi_common/tests/
  mlb_sgp/tests/` before proposing merge; executive review of
  `git diff main..HEAD`; **no merge without user approval**; bots are
  never touched (post-merge restart reminder — `_SCOPE_CACHE` pins
  scope verdicts until restart).

## Docs

- `kalshi_mlb_mm/README.md`: F5-winner pricing section (semantics table,
  conversion, anchor policy, decline behavior).
- Root `CLAUDE.md`: update the #86 sentence in the maker paragraph.
