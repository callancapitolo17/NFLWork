# Novig Vendor Guard for the Maker's SGP Consensus

**Date:** 2026-10-02
**Status:** Approved 2026-10-02 and built; awaiting merge approval
**Branch:** `feature/novig-vendor-consensus-guard` (cut from
`fix/novig-sgp-rest-discovery`, which made Novig price again and has since
merged to `main` as f23b01fe)

---

## 1. Problem

The maker (`kalshi_mlb_mm`) quotes a same-game combo when at least
`MIN_AGREEING_BOOKS` (2) books land and their devigged fairs sit within
`SIGMA_Z_MAX` (0.07) in probit space. Novig does not originate SGP prices: its
parlay endpoint relays a vendor book's price and names that vendor on every leg
of the response (`legs[].vendor`). So `{FanDuel, Novig-relaying-FanDuel}`
passes both gates while being one opinion.

The quorum is not the only damage. In a 3-book set a relayed copy also:

- **double-weights its vendor in the median**, and
- **shrinks the dispersion estimate**: the sample stddev of `{a, a, b}` is
  `0.58·|a−b|` against `0.71·|a−b|` for `{a, b}`. That loosens the dispersion
  gate and thins the `K_SIGMA·σ` margin term — a copy manufactures agreement.

## 2. Measurement (2026-10-02, maker's live path)

Harness: `SGPService.price_on_demand` at draftkings / fanduel / novig / betmgm
/ caesars in parallel per combo (as `OnDemandEngine._fetch_combo` does), with
each book's `price_selection_set`, `resolve_legs` and `NovigClient.submit_parlay`
wrapped to record every partition cell's decimal and Novig's per-leg vendor.
4 games on the 2026-10-03 slate (CWS@CLE, ATL@LAD, NYY@TB, SD@MIL), 14 shapes
per game (FG ML × total at three lines, FG spread ±1.5 × total, F5 spread ±0.5
× F5 total, FG ML × F5 total, FG total × 1st-inning 0.5, FG ML × 1st-inning,
two 3-leg sets, spread × ML), run twice 20 minutes apart (17:51 and 18:12 PT):
112 combo observations. No 403s; Novig's only declines were 400 "cannot
price".

| Fact | Both passes |
|---|---|
| Routing unit | **per cell**: each cell of a partition names its own vendor; every leg of one cell names the same vendor |
| Vendor share of 471 priced Novig cells | DraftKings 60%, FanDuel 29%, BetMGM 11%, Caesars <1% |
| Novig grids mixing 2-3 vendors | 91 of 104 (88%); 12 of the 13 single-vendor grids were DraftKings |
| Same cell, same vendor 20 minutes later | 162 of 231 (70%) — routing moves with prices, so no static shape→vendor map works |
| Shading (Novig cell ÷ vendor's own cell, same instant) | FanDuel cells 0.990 (IQR 0.983-0.994); BetMGM cells 1.000 |
| Routed vendor was the cheaper of FD/MGM | 109 of 139 cells (78%) — Novig usually takes the cheapest vendor quote |
| Novig vs a book it relayed in that combo: within σ_z 0.07 | **96 of 99 (97%)**, median σ_z 0.016 (~0.8 prob pts) |
| Independent baseline, FanDuel vs BetMGM: within σ_z 0.07 | 58 of 70 (83%), median σ_z 0.039 |
| DraftKings' own SGP endpoint | priced 0 of 112 (dark since #102) — Novig is our only read of DK |

Side fact: Novig's exchange singles are its own book, not relayed — CLE ML
traded at `available` 0.585 on Novig's market tree while the DK-relayed leg in
the same parlay response said 0.571. The leg surface (cross-game singles) is
therefore out of scope.

## 3. Rule

> **Novig's fair counts in a flight's consensus only if every cell it priced
> named a vendor and none of those vendors priced the same combo in that
> flight.** Otherwise Novig is dropped and its vendors speak for themselves.

One fixed rule, no knob. It keeps Novig exactly where it adds information:
as our only read of DraftKings (always), and of FanDuel/BetMGM alt lines our
own resolvers cannot price (Novig relayed FanDuel on FG totals 5.5/7.5 where
our FanDuel path declined).

Effect on the 112 measured combos (books landing inside the 10 s budget):

| Policy | Quotable | Too few books | Dispersion |
|---|---|---|---|
| Today (Novig always counts) | 94 | 8 | 10 |
| **This rule** | **78** | 21 | 13 |
| Never count Novig (config-only: drop it from `SGP_BOOKS`) | 58 | 42 | 12 |

Under the rule: 13 combos lose a quorum that was Novig + its own vendor, 6 now
fail dispersion (Novig was bridging a FanDuel-vs-BetMGM disagreement — e.g.
ATL@LAD ML × 8.5: FD 0.304, MGM 0.349, Novig 0.331 quoted at the median), 3
now pass (Novig was the dissenter, all Route B transfer fairs), and the 75
that quote both ways move a median 0.2 prob pts (max 1.0). Novig stays in 29
of the 104 flights it priced.

Rejected: *never count Novig* — simplest (no code), but drops 36 of 94
quotable combos instead of 16 and loses the only DraftKings read we have.
Rejected: *count Novig + vendor as one book* — the routing is per cell, so a
grid is usually 2-3 books at once; there is no single vendor to merge with.

## 4. Code

1. `mlb_sgp/novig_client.py` — `_parse_parlay_response` keeps each leg's
   `vendor` (pure parser).
2. `mlb_sgp/novig.py` — `price_selection_set(relayed=...)` records the
   vendors of each priced call into a per-call collector passed in by the
   caller (explicit, no client state); a priced leg naming no vendor is
   recorded as `"unknown"`.
3. `mlb_sgp/_shared.py` — `RelayedVendors` (the thread-safe per-call
   collector: partition cells price concurrently), `UNKNOWN_VENDOR`, and
   `OnDemandBookResult.vendors: tuple[str, ...] = ()`: the books whose prices
   this result relays, as our book keys; empty for a book that prices its own
   combos.
4. `kalshi_common/sgp_service.py` — `_price_on_demand` creates the collector,
   hands it to Novig's price hook, and stamps `vendors` on the result.
5. `kalshi_mlb_mm/router.py` — `drop_relayed_novig(book_fairs,
   novig_vendors)`: the rule, pure, on plain `{book: fair}` + a vendor tuple
   so both callers below share it.
6. `kalshi_mlb_mm/on_demand.py` — `lookup()` returns the fairs through
   `drop_relayed_novig`, so the quote, the confirm last look and the
   risk-sweep drift check all read through the one rule.
7. `kalshi_mlb_mm/main.py` — `on_demand_result` / `quote_priced` research
   payloads carry each book's `vendors` (and `live_games` the
   `consensus_books` that counted); `_on_demand_fill_info` prices from
   `lookup()` so research matches the quote path.
8. `kalshi_mlb_mm/report.py` — `universe_stats` (the daily report's
   approximation of the live gate) counts books per landing through the same
   rule when the payload carries vendors.

Tests: parser keeps vendors; Novig on-demand result carries the union over all
priced cells, other books carry `()`; the rule's cases (vendor in flight →
dropped; vendors absent → kept; unknown vendor → dropped; no Novig →
unchanged; Novig alone → kept or dropped); `lookup()` and `universe_stats`
apply it.

## 5. Not covered

- **Taker (`kalshi_mlb_rfq`)** — same flaw: `_load_book_fairs` gates on
  `MIN_BOOK_COUNT_FOR_BLEND=2` and medians sweep rows that carry no vendor.
  Dormant since 2026-06-29; the fix needs the vendor on `mlb_sgp_odds` rows.
  Do it when the taker is revived; a one-line warning goes next to
  `MIN_BOOK_COUNT_FOR_BLEND`.
- **Dashboard parlay blend** (`Answer Keys/mlb_correlated_parlay.R`) — medians
  model + up to 6 books including Novig; same double weight, no count gate.
- **Book-health rule B** counts Novig as a healthy book; it is a fetch-health
  alert, not a pricing gate. Accepted.

## 6. Version control, worktree, docs

- Worktree `.claude/worktrees/nice-rhodes-ec25c3`, branch
  `feature/novig-vendor-consensus-guard`; local `main` merged in before the
  pre-merge review.
- Commits, each carrying its own docs: (1) vendor capture (`mlb_sgp/` +
  `sgp_service` + result field, tests, `mlb_sgp/README.md`); (2) the rule in
  `router` + `on_demand.lookup` + research payloads + report, tests, maker
  README/config, taker warning, root `CLAUDE.md`, this spec.
- Docs in the same merge: `kalshi_mlb_mm/config.py` comment above
  `MIN_AGREEING_BOOKS` (replace "noted follow-up" with the rule),
  `kalshi_mlb_mm/README.md` (consensus section + decision log),
  `mlb_sgp/README.md` Novig section (per-cell routing facts),
  `kalshi_mlb_rfq/config.py` warning, root `CLAUDE.md` maker bullet.
- Pre-merge review per `CLAUDE.md`; merge, worktree removal and branch
  deletion only on explicit approval.
