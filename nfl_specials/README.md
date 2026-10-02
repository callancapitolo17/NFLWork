# NFL fecta pricer (Wagerzon trifectas and superfectas)

Prices Wagerzon's weekly NFL **trifectas** ("team wins 1Q, 1H and the game",
or a spread version) and **superfectas** (scores first + the three) off
sportsbook same-game-parlay prices, shows the edge and a Kelly stake, and
places the bet on Wagerzon from the page.

```
Wagerzon "NFL WEEK n - SPECIALS" ──▶ parse legs ──▶ FanDuel / BetMGM (HTTP) ─┐
                                                 └▶ DraftKings (dk_price_sidecar) ┴▶ fair, EV, stake ──▶ page :8096 ──▶ Place
```

## Run

```bash
dk_price_sidecar/run.sh     # own terminal; needed for superfectas (DK only)
nfl_specials/run.sh         # http://127.0.0.1:8096
```

The board re-prices every 15 min; Refresh forces one now. FanDuel and BetMGM
price every trifecta in ~90 s; DraftKings goes last, paced 1.5 s per call, and
fills in as it goes.

## How a special is priced

1. **Legs.** `special_parser.py` reads `SEAHAWKS TRIFECTA (1Q -½, 1H -4½ & GM -7½)`.
   Quarter and half wins are **3-way**: a tied quarter or half loses. With
   integer scores that is the same bet as -0.5, so either market works.
2. **Partition.** At each book every combination of the legs' outcomes is
   priced as an SGP (8 cells when every leg is two-way; 18 for a trifecta
   off 3-way period markets; 36 for a DK superfecta).
3. **Devig.** Probit devig across the whole partition (`kalshi_common`'s
   n-way probit). The implied sum is the book's real SGP hold — 1.24-1.40
   measured on 2026-10-02, far above compounded single-leg vig — so pricing
   the special's own SGP and stripping single-leg vig would overstate every
   edge. A partition with a declined cell or an implied sum outside
   `[1, 1 + 0.25 x legs]` gives that book no fair.
4. **Consensus** = mean of the books that priced the full partition.
5. **Stake** = bankroll x Kelly fraction x full Kelly at Wagerzon's price.
   Per (game, team) only the best special by expected log growth keeps a
   stake ("overlaps" on the rest): a team's fectas win together.

| Book | Trifecta | Superfecta | How |
|---|---|---|---|
| FanDuel | yes | no — "Team to Score First" is not SGP-eligible | 3-way period winners + ML / spreads, `implyBets` |
| BetMGM | yes | no — no first-to-score market | ±0.5 period spreads + ML / spreads, `tv2Picks` |
| DraftKings | yes | yes ("1st to Score") | 3-way period markets via the sidecar's real Chrome |

## Placing

**Place submits a real Wagerzon bet.** The server re-checks the special is on
the board at the price the page showed and that the game has not started, then
`wagerzon_odds/single_placer.place_single` previews it (`ConfirmWagerHelper`,
win must match to $0.01 or it refuses) and submits it. A special is a one-sided
prop on the "home" slot, so it goes in as Play=5 with no points. Every attempt
is logged, whatever Wagerzon answered. Wagerzon's minimum online wager is above
$10 (a $25 preview passes).

## Files

| File | Job |
|---|---|
| `wz.py` | scrape the week's fectas; place one |
| `special_parser.py` | description -> legs (anything unrecognized is refused, never guessed) |
| `books.py` | adapter interface: a leg = an exhaustive outcome group |
| `fd_book.py`, `mgm_book.py`, `dk_book.py` | per-book market lookup + SGP price |
| `pricing.py` | partition devig, EV, Kelly |
| `board.py` | one refresh; read-time sizing |
| `store.py` | DuckDB state |
| `app.py`, `static/index.html` | loopback web page + JSON API |

## Data — `nfl_specials/nfl_specials.duckdb`

| Table | Write | Row |
|---|---|---|
| `fecta_quotes` | APPEND every refresh | one (special, book) fair or the reason there is none |
| `placed_fectas` | APPEND every Place click | the attempt and Wagerzon's answer |
| `settings` | UPSERT | bankroll, Kelly fraction |

## Troubleshooting

- **DraftKings chip says sidecar not running** — start `dk_price_sidecar/run.sh` from
  your own terminal (a sandboxed process cannot draw Chrome's window and its
  price calls hang).
- **DraftKings chip shows HTTP 502 "Failed to fetch"** — DK denied the page. The
  sidecar reloads and retries once; if it keeps happening, check its `/health`
  (`blocks`, `reloads`, `calls_since_reload`).
- **A superfecta says "no SGP market for '... scores first'"** at FD/BetMGM — expected;
  only DK prices superfectas.
- **Nothing on the board** — Wagerzon posts the specials Thursday-ish; the league is
  found by name (`NFL WEEK <n> - SPECIALS`).

## Not done yet

- **Cloud.** DK only prices to a real Google Chrome; Playwright's Linux Chromium under
  Xvfb was denied (2026-10-02). Running remotely needs a host with real Chrome, untested.
- Correlation across games is ignored (cross-game NFL correlation measured ~0).
