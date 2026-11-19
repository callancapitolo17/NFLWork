# SGP Compare

Build an NFL same-game parlay, price it at each venue, and see each venue's
**Correlation factor**: how much more often the venue thinks the legs hit
together than if they had nothing to do with each other.

Research tool for Cal's own knowledge (2026-11-18): no fair price, no edge, no
bets. Every venue call is a read or a price request.

## Run

```bash
pip install curl_cffi duckdb scipy
sgp_compare/run.sh            # http://127.0.0.1:8097 (Mac first; VM container later)
```

Pick a game, search or browse the markets (tick **alt lines** for every
number), tap one side per market, press **Price**.

## What each number means

| Column | Meaning |
|---|---|
| Parlay price | The venue's best offer, with every price tier and the most it takes at each (ProphetX: bigger size, worse price). |
| Legs as separate bets | The legs' single prices multiplied together; "no vig" uses each single devigged against its other side (kalshi_common probit; none when the two sides are crossed or more than 25% apart). |
| Correlation factor, with vig | parlay chance / product of the single chances - 1, straight off takeable prices. Mixes correlation with the extra margin a venue charges on parlays. |
| Correlation factor, vig removed | From the venue's parlay price for every combination of the legs' outcomes (2^N price requests, 2-3 legs only): each combination's share of the total, parlay's share / product of each leg's share - 1. One price list on both sides of the ratio, so a margin spread evenly over the combinations cancels exactly. Blank, with the reason, when a combination is declined or fails (the first one stops the rest), a leg has no other side, a line is a whole number or a half/quarter moneyline (it can push or tie), the combinations add to more than 100% + 25% per leg (a stale price), or there are 4+ legs. ProphetX won't put a TD "NO" in a parlay, so a TD leg never gets this number. |
| Parlay vig | Sum of every combination's chance - 1. |
| Inside parlay | ProphetX folds its correlation charge into each leg's price inside the parlay. |

## Layout

| File | Job |
|---|---|
| `legs.py` | `Game`, `CanonicalLeg` (venue-neutral: period, stat, subject, side, line), `OfferedLeg`, `ParlayQuote`. |
| `venues/base.py` | The `Venue` interface: `list_games`, `list_legs`, `price_parlay`. |
| `venues/prophetx.py` | ProphetX: anonymous public endpoints over a curl_cffi Chrome session; game lines, halves, quarters, team totals, ~290 player props per game. |
| `venues/draftkings.py` | DraftKings: the SGP endpoints `nfl_specials/dk_book.py` already uses (no login); every SGP market of the game, ~3,200 legs incl. player props and milestone ladders. |
| `venues/__init__.py` | `build_venues()`: the venues the page uses. Adding a venue = a module + one line here. |
| `correlation.py` | Pure math: the Correlation factor with and without vig (singles devig from `kalshi_common.fair_value.two_way_fair`). |
| `pricing.py` | `price_at_venue`: fresh legs, the parlay price, the 2^N partition. |
| `store.py` | DuckDB `$SGP_COMPARE_DB` (default `sgp_compare/sgp_compare.duckdb`), table `parlay_quotes`, APPEND one row per venue per Price click: requested leg ids, legs whole (canonical key, venue ids, both sides' prices), tiers, partition prices, both factors with their reasons, raw venue response. Bad requests (400) write nothing; a store failure still shows the price with "not saved". |
| `app.py`, `static/index.html` | Loopback page and JSON API (Host allowlist, JSON-only POSTs). |

## Venues

| Venue | Status |
|---|---|
| ProphetX | Live (PR 1). Parlay = `POST /parlay/public/api/v1/user/request` with stake 1, a price request only. |
| Novig, Kalshi, Polymarket US | Next, one thread each (Novig relays DK/FD/MGM/CZR prices; Kalshi needs an RFQ). |
| DraftKings | Live. Parlay = DK's SGP-builder GET (`.../sgp/dkuswv/sportsdata/v2/sgp`): base legs + the last leg's market, answering every side of that market, so the 2^N combinations take 2^(N-1) requests (answers kept 10 s), at most one per 500 ms like DK's own widget. One price, no stated limit, no per-leg prices. Singles are DK's listed prices: the builder prices a single a bit worse (HOU ML -750 vs -700 listed, 2026-11-19). DK refuses data-center IPs: run on the Mac, or set `SGP_COMPARE_DK_PROXY` to the Mac's `nfl_specials` dk_relay (`http://<mac tailscale ip>:8097`). |
| FanDuel, BetMGM, Caesars | After the prediction markets. |

Cross-venue pricing (one parlay priced at every venue) matches legs by
`CanonicalLeg`; it lands with the second venue.

## Tests

```bash
python3 -m pytest -q sgp_compare/tests
```

Fixtures are trimmed copies of the 2026-11-18 ProphetX recon and the 2026-11-19
DraftKings recon (both IND at HOU, read-only, from Cal's Mac).

## Troubleshooting

- **Error in the game list, or HTTP 502 from /api/legs**: ProphetX unreachable, answered
  non-200, or sent a shape the parser does not expect; the message says which. The cloud
  sandbox cannot reach ProphetX; run on the Mac.
- **"no longer offered"**: the line moved between picking and pricing; reload the game.
- **"would not price this combination"**: ProphetX returned no offers (it refuses some combinations).
  At DraftKings: DK dropped a leg (`selectionsNotMapped`) or left the last leg's market out of
  `compatibleMarkets`; it never says "no" outright.
- **DraftKings "expected HTTP 200, got 403"**: the request left from a data-center IP. Run on the Mac,
  or point `SGP_COMPARE_DK_PROXY` at the dk_relay. **404 on SGP price**: DK moved its widget's API
  (see `nfl_specials/README.md` troubleshooting).
- **Vig-removed factor blank**: the reason is shown under it.
