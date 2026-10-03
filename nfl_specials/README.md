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

The board prices only when you click **Refresh** (no auto-refresh). FanDuel and
BetMGM price the trifectas in ~90 s; DraftKings then prices the superfectas
(~20 calls each, paced 1.5 s, page reloaded every 5 calls), ~7 more minutes,
filling in as it goes.

## How a special is priced

0. **Game.** Every book's game must be the special's: kickoff inside the
   Wagerzon week (Thursday to Monday night around the specials' date) and, at
   each book, the same two teams within 12 hours of the board's kickoff —
   books list future weeks, so otherwise a started or missing game would be
   priced off next week's.
1. **Legs.** `special_parser.py` reads `SEAHAWKS TRIFECTA (1Q -½, 1H -4½ & GM -7½)`.
   Quarter and half wins are **3-way**: a tied quarter or half loses. With
   integer scores that is the same bet as -0.5, so either market works.
2. **Partition.** At each book every combination of the legs' outcomes is
   priced as an SGP (8 cells when every leg is two-way; 12-18 for a trifecta
   off 3-way period markets).
3. **Devig.** Probit devig across the whole partition (`kalshi_common`'s
   n-way probit). The implied sum is the book's real SGP hold — 1.24-1.40
   measured on 2026-10-02, far above compounded single-leg vig — so pricing
   the special's own SGP and stripping single-leg vig would overstate every
   edge. A partition with a declined cell or an implied sum outside
   `[1, 1 + 0.25 x legs]` gives that book no fair.
4. **Fair = worst case**: the lowest probability among the books that priced
   the full partition (a special is only ever backed; user decision
   2026-10-02).
   **Superfectas price at DraftKings only** (user decision 2026-10-02: only DK
   lets "scores first" into an SGP): P(trifecta part) x P(scores first |
   trifecta part). The first factor is DK's own partition fair for the other
   three legs (12-18 cells); the second comes from TWO DK SGPs — (team scores
   first + the trifecta legs) and (opponent scores first + the same legs) —
   whose vig cancels in the ratio. DK's shares on 2026-10-02 were 0.83-0.85;
   history (2011-2025, 1,970 team-games that won 1Q, 1H and the game) says
   0.893 ± 0.007, so DK's number is the conservative one. ~20 DK calls per
   superfecta instead of a 36-cell partition.
5. **Stake** = bankroll x Kelly fraction x full Kelly at Wagerzon's price.
   Per (game, team) only the best special by expected log growth keeps a
   stake ("overlaps" on the rest): a team's fectas win together. The stakes
   are then **fitted to the Wagerzon balance available** (read live from
   Wagerzon when a refresh starts and ends, so pending bets and the week's
   results already count; an unreadable balance recommends nothing): if they
   add up past it, every bet's marginal log growth must clear one common
   hurdle, raised until they fit — weaker edges shrink first and drop to
   zero. A stake under Wagerzon's $20 minimum for specials (measured: $15-19
   rejected) is dropped and the rest re-fit. The page shows the uncapped
   Kelly stake next to any trimmed one.

| Book | Trifecta | Superfecta | How |
|---|---|---|---|
| FanDuel | yes | no — "Team to Score First" is not SGP-eligible | 3-way period winners + ML / spreads, `implyBets` |
| BetMGM | yes | no — no first-to-score market | ±0.5 period spreads + ML / spreads, `tv2Picks` |
| DraftKings | not used (calls are scarce) | yes — the only book ("1st to Score") | via the sidecar's real Chrome |

## Placing

**Place submits a real Wagerzon bet.** The server re-checks the special is on
the board at the price the page showed, that the game has not started, that
the stake is at least $20 and within the available balance, and that the same
special was not placed in the last 2 minutes (a slow answer invites a second
click; the button is also disabled while a placement is in flight), then
`wagerzon_odds/single_placer.place_single` previews it (`ConfirmWagerHelper`,
win must match to $0.01 or it refuses) and submits it. A special is a one-sided
prop on the "home" slot, so it goes in as Play=5 with no points. Every attempt
is logged, whatever Wagerzon answered; a placed bet lowers the page's available
balance right away. Wagerzon's minimum online wager on a special is $20.

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
  sidecar reloads every 5 calls and retries once after a denial; if it keeps
  happening, check its `/health` (`blocks`, `reloads`, `calls_since_reload`).
- **Restarting the sidecar fails with "Opening in existing browser session"** — an
  old sidecar (or its Chrome) still holds `~/.dk_price_sidecar/profile`; stop it first.
- **A superfecta says "no SGP market for '... scores first'"** at FD/BetMGM — expected;
  only DK prices superfectas.
- **Nothing on the board** — Wagerzon posts the specials Thursday-ish; every league
  named `NFL WEEK <n> - SPECIALS` is read (two can be up at once).
- **"Available: unknown"** — the Wagerzon balance read failed, so no stakes are
  recommended; click Refresh.

## Not done yet

- **Cloud.** DK only prices to a real Google Chrome; Playwright's Linux Chromium under
  Xvfb was denied (2026-10-02). Running remotely needs a host with real Chrome, untested.
- Correlation across games is ignored (cross-game NFL correlation measured ~0).
