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
| Legs as separate bets | The legs' single prices multiplied together; "no vig" uses each single devigged against its other side. |
| Correlation factor, with vig | parlay chance / product of the single chances - 1, straight off takeable prices. Mixes correlation with the extra margin a venue charges on parlays. |
| Correlation factor, vig removed | Same, with the parlay devigged across every combination of its legs' outcomes (2^N price requests, 2-3 legs only) and each single devigged (probit). Blank, with the reason, when a combination is declined, a leg has no usable other side (none, crossed, or wider than 25%), or there are 4+ legs. |
| Parlay vig | Sum of every combination's chance - 1. |
| Inside parlay | ProphetX folds its correlation charge into each leg's price inside the parlay. |

## Layout

| File | Job |
|---|---|
| `legs.py` | `Game`, `CanonicalLeg` (venue-neutral: period, stat, subject, side, line), `OfferedLeg`, `ParlayQuote`. |
| `venues/base.py` | The `Venue` interface: `list_games`, `list_legs`, `price_parlay`. |
| `venues/prophetx.py` | ProphetX: anonymous public endpoints over a curl_cffi Chrome session; game lines, halves, quarters, team totals, ~290 player props per game. |
| `venues/__init__.py` | `build_venues()`: the venues the page uses. Adding a venue = a module + one line here. |
| `correlation.py` | Pure math: probit devig, correlation factor with and without vig. |
| `pricing.py` | `price_at_venue`: fresh legs, the parlay price, the 2^N partition. |
| `store.py` | DuckDB `$SGP_COMPARE_DB` (default `sgp_compare/sgp_compare.duckdb`), table `parlay_quotes`, APPEND one row per venue per Price click, raw venue response included. |
| `app.py`, `static/index.html` | Loopback page and JSON API (Host allowlist, JSON-only POSTs). |

## Venues

| Venue | Status |
|---|---|
| ProphetX | Live (PR 1). Parlay = `POST /parlay/public/api/v1/user/request` with stake 1, a price request only. |
| Novig, Kalshi, Polymarket US | Next, one thread each (Novig relays DK/FD/MGM/CZR prices; Kalshi needs an RFQ). |
| DraftKings, FanDuel, BetMGM, Caesars | After the prediction markets. |

Cross-venue pricing (one parlay priced at every venue) matches legs by
`CanonicalLeg`; it lands with the second venue.

## Tests

```bash
python3 -m pytest -q sgp_compare/tests
```

Fixtures are trimmed copies of the 2026-11-18 ProphetX recon (IND at HOU).

## Troubleshooting

- **HTTP 502 on a ProphetX read**: ProphetX unreachable or answered non-200; the
  message carries the status and body. The cloud sandbox cannot reach ProphetX; run on the Mac.
- **"would not price this combination"**: ProphetX returned no offers (it refuses some combinations).
- **Vig-removed factor blank**: the reason is shown under it.
