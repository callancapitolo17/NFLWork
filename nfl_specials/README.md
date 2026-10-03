# NFL fecta pricer (Wagerzon trifectas and superfectas)

Prices Wagerzon's weekly NFL **trifectas** ("team wins 1Q, 1H and the game",
or a spread version) and **superfectas** (scores first + the three) off
sportsbook same-game-parlay prices, shows the edge and a Kelly stake, and
places the bet on Wagerzon from the page.

```
Wagerzon "NFL WEEK n - SPECIALS" ──▶ parse legs ──▶ FanDuel / BetMGM (HTTP) ─┐
                                                 └▶ DraftKings (HTTP) ────────┴▶ fair, EV, stake ──▶ page :8096 ──▶ Place
```

## Run

```bash
nfl_specials/run.sh         # http://127.0.0.1:8096
```

The board prices only when you click **Refresh** (no auto-refresh). FanDuel and
BetMGM price the trifectas in ~90 s; DraftKings then prices the superfectas
(7-10 requests each, 0.5 s apart), ~40 s for six, filling in as it goes.

## The page

- **Top bar:** Wagerzon account (hidden with one account), bankroll and Kelly
  fraction behind one button (saved on change), Refresh.
- **Tiles:** Wagerzon available, the recommended total (with its share of
  available) and the expected profit of those stakes at the worst-case fair.
- **Progress:** while a refresh runs the page polls every 3 s and shows done /
  total overall and per book (FD / MGM / DK); it stops polling when the refresh
  ends. Afterwards only a book that failed shows.
- **Board:** grouped by game in kickoff order (started games and specials with
  no game last), or one list by EV ("Best EV"); filter All / Trifectas /
  Superfectas and "+EV only" (remembered in the browser). Each row: legs as
  chips, every book's fair with the worst case marked (a superfecta shows DK's
  1Q/1H/GM part, the team-scores-first %, and DK's own price incl. vig), WZ
  price, fair, EV (shaded by size), and the stake box prefilled with the
  recommended stake, "to win", and the uncapped Kelly when trimmed.
- **Place:** a confirm sheet (risk, to win, price, fair, EV, account), then the
  button shows "Placing…" until Wagerzon answers; the ticket or the error shows
  on the row. **Placed this week** lists only bets Wagerzon took (refused
  attempts stay in `placed_fectas`).

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
   0.893 ± 0.007, so DK's number is the conservative one. 7-10 DK requests
   per superfecta: each prices every outcome of the cell's last leg.
5. **Stake** = bankroll x Kelly fraction x full Kelly at Wagerzon's price.
   Per (game, team) only the best special by expected log growth keeps a
   stake ("overlaps" on the rest): a team's fectas win together. The stakes
   are then **fitted to the Wagerzon balance available** (read live from
   Wagerzon when a refresh starts and ends, so pending bets and the week's
   results already count; an unreadable balance recommends nothing): if they
   add up past it, every bet's marginal log growth must clear one common
   hurdle, raised until they fit — weaker edges shrink first and drop to
   zero. A stake under Wagerzon's $20 minimum for specials (measured: $15-19
   rejected) is dropped and the rest re-fit; each stake is capped at the
   $250 maximum, and what a capped bet can't take goes to the next edge. The page shows the uncapped
   Kelly stake next to any trimmed one. No stake (and no budget) goes to a
   special whose game has started, or to any special of a team the account
   already has a placed bet on this week (that bet is already in Wagerzon's
   available; a second stake would double one position).

| Book | Trifecta | Superfecta | How |
|---|---|---|---|
| FanDuel | yes | no — "Team to Score First" is not SGP-eligible | 3-way period winners + ML / spreads, `implyBets` |
| BetMGM | yes | no — no first-to-score market | ±0.5 period spreads + ML / spreads, `tv2Picks` |
| DraftKings | not used | yes — the only book ("1st to Score") | its SGP builder's price endpoint (below) |

**DraftKings, without a browser (2026-10-02).** DK's betslip price call
(`calculateBets`) answers only a real Chrome (issue #102). DK's SGP *builder*
widget prices through a different, ungated call:

```
GET https://sportsbook-nash.draftkings.com/sites/US-WV-SB/api/sportscontent/sgp/dkuswv/sportsdata/v2/sgp
    ?eventId=<event>&selections=<base ids, comma>&marketCandidates=<market id>&oddsStyle=american
    header X-SportId: 3        (404 without it)
```

It returns the base SGP's `trueOdds` and, under `compatibleMarkets`, the
`trueOdds` of base + each outcome of the candidate market — so one request
prices every partition cell that differs only in its last leg, and both
sides of "scores first". Verified against `calculateBets` the same minute:
10/10 cells identical, and all six Week 5 superfectas' fairs identical to 4
decimals. **If DK cannot combine a base leg it drops it and prices the rest**
(`selectionsNotMapped`); the adapter treats that as a decline, never as a
price.

## Placing

**Place submits a real Wagerzon bet.** The server re-checks the special is on
the board at the price the page showed, that the game has not started, that
the stake is at least $20 and within the available balance, and that the same
special was not placed in the last 2 minutes (a slow answer invites a second
click; the button is also disabled while a placement is in flight), then
`wagerzon_odds/single_placer.place_single` previews it (`ConfirmWagerHelper`,
win must match to $0.01 or it refuses) and submits it. A special is a PROP game,
so it goes in with its rotation number as the Play and no points — what
Wagerzon's own site sends. (Play=5 passes the preview but the submit fails
"Couldn't find Game Line".) Every attempt
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

- **DraftKings chip shows "DK SGP price ... returned HTTP 404"** — the request lost
  its `X-SportId` header, or DK moved its SGP widget's API. Re-read the widget
  (`dk-same-game-parlay/<version>/samegameparlayweb-external.js`, linked from the
  sportsbook page's `sameGameParlayConfig`) for `getYourbetRequestParameters`.
- **A superfecta says "no SGP market for '... scores first'"** at FD/BetMGM — expected;
  only DK prices superfectas.
- **Nothing on the board** — Wagerzon posts the specials Thursday-ish; every league
  named `NFL WEEK <n> - SPECIALS` is read (two can be up at once).
- **"Available: unknown"** — the Wagerzon balance read failed, so no stakes are
  recommended; click Refresh.

## Not done yet

- **Cloud.** Every book is plain HTTP now, but only the home IP has been tested;
  Wagerzon logins from a data-center IP are untested.
- Correlation across games is ignored (cross-game NFL correlation measured ~0).
