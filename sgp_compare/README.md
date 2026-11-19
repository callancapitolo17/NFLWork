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
| `venues/fanduel.py` | FanDuel: the event list, nine event-page tabs merged (SGP markets only; any tab failing fails the read), and `implyBets` for the parlay price, on at most 3 reused sessions; headers and FanDuel's PerimeterX token come from `mlb_sgp/scraper_fanduel_sgp.py`. Main and alternate lines of one bet share a `market_key`, so the alt-lines toggle groups them. |
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
| FanDuel | Parlay = `implyBets` (`pricePolicy=SUGGESTED`), a bet-slip price request only; one price with FanDuel's max stake. Ladder props ("5+ Receptions") and TD scorers have no other side, so no vig-removed number. Moneyline + spread on one team is refused. |
| DraftKings, BetMGM, Caesars | After the prediction markets. |

Cross-venue pricing (one parlay priced at every venue) matches legs by
`CanonicalLeg`; it lands with the second venue.

## Tests

```bash
python3 -m pytest -q sgp_compare/tests
```

Fixtures are trimmed copies of the 2026-11-18 ProphetX recon and the 2026-11-19 FanDuel recon (IND at HOU).

## Troubleshooting

- **Error in the game list, or HTTP 502 from /api/legs**: the venue was unreachable, answered
  non-200, or sent a shape the parser does not expect; the message says which. The cloud
  sandbox cannot reach the venues; run on the Mac.
- **"no longer offered"**: the line moved between picking and pricing; reload the game.
- **"would not price this combination"**: the venue declined it (ProphetX returned no offers;
  FanDuel answered `INVALID_SGM_COMBINATION`, e.g. a moneyline and spread on one team).
- **"FanDuel moved ... while pricing"**: a main line changed between the leg read and the
  parlay request; reload the game.
- **FanDuel errors with HTTP 400/403 everywhere, or "no SGP entry"**: FanDuel's PerimeterX token (`FD_PX_CONTEXT` in
  `mlb_sgp/scraper_fanduel_sgp.py`) has expired; refresh it there (it is shared with the MLB scrapers).
  FanDuel also needs a US IP.
- **Vig-removed factor blank**: the reason is shown under it.
