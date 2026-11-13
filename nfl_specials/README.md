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

On the phone or any tailnet device: **https://mlb-stack.tail87ab2d.ts.net:8443/**
(the Oracle VM; see Running on the VM). For development on the Mac:

```bash
nfl_specials/run.sh         # http://127.0.0.1:8096
```

The board prices only when you click **Refresh** (no auto-refresh). The three
books price side by side, one thread each (`board.py::_price_lane`), so the
refresh takes about as long as the slowest book. Measured 2026-10-02 (17
specials): BetMGM done in 0.5 min, FanDuel in 0.9 min; DraftKings' six
superfectas (7-10 requests each, 0.5 s apart) took ~40 s over plain HTTP.

## Running on the VM

The page runs as the `fectas` container of the VM's Docker stack
(`unabated_ticket/deploy/`, next to the bets service), bound to the VM's
loopback and published on the tailnet by `tailscale serve` at port 8443.

```
phone ──tailnet──▶ VM :8443 ──▶ fectas (:8096) ──▶ FanDuel, BetMGM, Wagerzon (direct)
                                              └──▶ Mac dk_relay (:8097) ──▶ DraftKings
```

- **DraftKings goes through the Mac.** DraftKings answers "Access Denied" to
  Oracle's data-center IPs, so the container sends its DraftKings requests to
  `dk_relay.py` on the Mac over Tailscale (`NFL_SPECIALS_DK_PROXY` in the
  VM's `deploy/.env`). **Mac asleep or off, or abroad (DraftKings refuses
  non-US IPs): the DraftKings chip errors and superfectas fall back to
  FanDuel + BetMGM with a fixed scores-first share (Pricing step 4).**
- **Wagerzon login:** only the `WAGERZON*_USERNAME` / `_PASSWORD` lines of the
  VM's `unabated_ticket/bets_service/.env` (`NFL_SPECIALS_WAGERZON_ENV`), the
  same login the bets service uses.
- **State:** `~/unabated-data/nfl_specials.duckdb` on the VM.
- **Update and setup:** `unabated_ticket/deploy/README.md` (`deploy.sh`, step 12).
- Use one page at a time: the Mac's `run.sh` copy has its own database, so it
  does not see bets placed from the VM's page.

The relay on the Mac runs as a launchd agent (starts at login, restarts if it
dies, waits for Tailscale):

```bash
# install once (after editing the plist, run the bootout line first)
cp nfl_specials/com.nflwork.dk-relay.plist ~/Library/LaunchAgents/
launchctl bootstrap gui/$(id -u) ~/Library/LaunchAgents/com.nflwork.dk-relay.plist
# restart onto new code
launchctl kickstart -k gui/$(id -u)/com.nflwork.dk-relay
# stop and remove
launchctl bootout gui/$(id -u) ~/Library/LaunchAgents/com.nflwork.dk-relay.plist
```

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
  chips, every book's fair with the worst case marked (a superfecta shows its
  one DK fair, or with DK down its trifecta part's FD and MGM fairs and a
  `share 0.78 fixed, DK down` chip), WZ
  price, fair, EV (shaded by size), and the stake box prefilled with the
  recommended stake and its stake figures (see the stake column below).
- **Place:** a confirm sheet (risk, to win, price, fair, EV, account), then the
  button shows "Placing…" until Wagerzon answers; the ticket or the error shows
  on the row. **Placed this week** lists only bets Wagerzon took (refused
  attempts stay in `placed_fectas`).

## How a special is priced

0. **Game.** Every book's game must be the special's: kickoff inside the
   Wagerzon week (Thursday to Monday night around the first Sunday on or after
   the specials' date) and, at
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
   **DraftKings down** (it failed to open this refresh; user decision
   2026-10-22): the trifecta part prices at FanDuel + BetMGM exactly as a
   trifecta (worst case of the two) and the share is a fixed **0.776**
   (`board.DK_DOWN_SCORES_FIRST_SHARE`), the lowest DK ever quoted (122
   quotes in `fecta_quotes.sf_share`: 0.776-0.883, median
   0.842). Each superfecta appends one `fecta_quotes` row with `book =
   'fixed'`, `quote_kind = 'scores_first_share'`, `sf_share = 0.776` and the
   DK error in `reason`. DK failing mid-refresh does not fall back: those
   superfectas stay unpriced until the next Refresh.
5. **Stake** = bankroll x Kelly fraction x full Kelly at Wagerzon's price,
   sized **per (game, team) jointly** (user decision 2026-10-15, replacing
   "one special per team"): a team's fectas are nested — the stricter one
   winning means the looser one won (GM -7½ => GM; a superfecta => its
   trifecta part) — so the group is one chain where the k-th likeliest bet
   wins only if every likelier one did (exact when nested, the maximum
   correlation the fairs allow otherwise). Each stake is the Kelly stake
   GIVEN the others and the bets **already placed** on that team this week:
   a second special gets money only when it still adds log growth on top of
   the first, a placed bet shrinks the next stake instead of blocking it, and
   the same special is topped up to Kelly (never past Wagerzon's $250 per
   special). A placed bet whose special the board cannot price still blocks
   its team. A special the joint Kelly leaves at zero yields to the team's
   biggest position. Opposite teams are separate groups. The stakes are then
   **fitted to the Wagerzon balance available** (read live from
   Wagerzon when a refresh starts and ends, so pending bets and the week's
   results already count; an unreadable balance recommends nothing): if they
   add up past it, every bet's marginal log growth must clear one common
   hurdle, raised until they fit — weaker edges shrink first and drop to
   zero. A stake under
   Wagerzon's $20 minimum for specials (measured: $15-19
   rejected) is dropped and the rest re-fit; each stake is capped at the
   $250 maximum, and what a capped bet can't take goes to the next edge. No
   stake (and no budget) goes to a special whose game has started.

   Each row's stake column shows, in a
   fixed order, the Recommended stake, the stake with no balance limit (same
   sizing, $20 minimum and $250 cap kept) and the uncapped joint Kelly stake,
   even when they match or are $0; what is already bet on that special as a
   green chip; and at most one plain reason when something holds the stake
   back (amber when it is $0 for a reason worth a look, grey when it is only
   capped, trimmed, full at $250 or waiting for the balance). A line the $20
   minimum zeroes (`under_min_stake` on the payload) shows no reason when
   its team has other specials, rather than blaming a sibling; a lone line
   still says its Kelly is under the minimum. A
   team's specials in one game sit together; with more than one, a strip on
   top shows the team's total bet and recommended stake and its one Kelly
   figure, and the rows drop their own Kelly line.

| Book | Trifecta | Superfecta | How |
|---|---|---|---|
| FanDuel | yes | no — "Team to Score First" is not SGP-eligible | 3-way period winners + ML / spreads, `implyBets` |
| BetMGM | yes | no — no first-to-score market | ±0.5 period spreads + ML / spreads, `tv2Picks` |
| DraftKings | not used | yes — the only book ("1st to Score") | its SGP builder's price endpoint (below) |

With DraftKings down, FanDuel and BetMGM price a superfecta's trifecta part
(Pricing step 4).

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
the stake is at least $20 and within the available balance, that the account's
tickets on this special this week plus this stake stay within Wagerzon's $250
per special (user, 2026-10-15: the maximum is per special, not per ticket), and
that the same special was not placed in the last 2 minutes (a slow answer invites a second
click; the button is also disabled while a placement is in flight), then
`wagerzon_odds/single_placer.place_single` previews it (`ConfirmWagerHelper`,
win must match to $0.01 or it refuses) and submits it. A special is a PROP game,
so it goes in with its rotation number as the Play and no points — what
Wagerzon's own site sends. (Play=5 passes the preview but the submit fails
"Couldn't find Game Line".) The server keeps one Wagerzon login per account; when
the preview comes back logged out (an expired login), it logs in again and retries
once (nothing was submitted yet). A logout at the submit is never retried: that bet
may have gone through (2026-11-08). Every attempt
is logged, whatever Wagerzon answered; a placed bet lowers the page's available
balance right away. Wagerzon's minimum online wager on a special is $20.

### Free play (Cal, 2026-11-12)

**Use free play** in Settings (saved as `settings.use_free_play`) bets the
Wagerzon account's free play balance instead of cash:

- **Balance:** `FreePlayAmount` off the same `PlayerInfoHelper` read as the cash
  `RealAvailBalance` (`wagerzon_balance.BalanceSnapshot.free_play`), read when a
  refresh starts and ends, lowered by each free play bet.
- **Stakes:** the same Kelly stakes on the same cash edges, fitted to the free play
  balance instead of the cash one ($20 minimum and $250 per special unchanged) —
  the rule Cal chose for BFA free play (option A, 2026-11-09). The Available tile
  becomes Free play; Expected profit becomes free play value, stake x fair x net
  odds (a free play's win pays only the to-win and a loss costs no cash). A held
  free play bet is sized around as if it were cash (conservative: it assumes a
  loss the free play never costs).
- **Place:** the request carries `free_play`, refused unless it matches the saved
  setting or when the stake is over the free play balance (or the balance is
  unknown), checked inside the one-at-a-time placement lock; `single_placer` then
  sends `useFreePlayAmount` true on the preview and the submit.
- **Which pot it drew on:** after a placed free play the server reads both balances
  again. Free play fell by the stake: recorded `free_play` TRUE. Cash fell instead:
  recorded FALSE and the page warns "Wagerzon charged this to cash". Neither, or
  the read failed: recorded TRUE (as sent) with a "couldn't confirm" warning.
  Placed this week tags free play rows. `placed_fectas.ev` stays the cash EV.
- **Unverified until the first live click:** that Wagerzon takes free play on a
  special, and that its preview quotes the same to-win (a free play's stake is not
  returned on a win; if the preview's Win differs, the placer refuses it as a moved
  price and nothing is submitted).

## 1H parlays (Cal, 2026-11-12)

The **1H Parlays** tab prices Wagerzon's same-game **1st-half spread + 1st-half
total** parlays for the NFL and college football. Wagerzon prices them as if the
two legs were independent, minus a flat ~0.57% (measured 2026-11-12 on every
pairing), so the pairings that win together (favorite + over, underdog + under)
can pay more than their true chance. Quarter lines are refused in any Wagerzon
parlay (`WAGERTYPENOALLOWGAME`), and 2nd halves are posted only at halftime, so
only 1st halves are read.

1. **Wagerzon read** (`wz_halves.py`): the catalog leagues "NFL - 1ST HALVES" and
   "NCAA FOOTBALL - 1ST HALVES", one request each (a combined `lg=` answers only
   the first sport). Each 1H line is its own Wagerzon game (`idgmtyp` 15).
2. **Whole-number lines** are priced off the half-point lines either side and
   need two books; see "1H parlays on whole-number lines" below.
3. **Fair**: the game's four corners (away or home spread x over or under) are
   the four cells of ONE 2 x 2 partition at each book, at Wagerzon's exact
   numbers, probit-devigged together (`halves_board.price_half_at_book`, the
   fectas' `pricing.price_partition`). Fair = the lowest per corner among the
   books that priced the partition.

   | Book | 1H spread | 1H total |
   |---|---|---|
   | FanDuel | "1st Half Spread" or "1st Half Alternate Spread" | "1st Half Total" (main line only) |
   | BetMGM | "1st half spread" (one market per line) | "1st half totals" (one market per line) |
   | DraftKings | "Spread 1st Half" or "Spread Alternate - 1st Half" | "Total 1st Half" or "Total Alternate - 1st Half" |

   College lists by school name: both of Wagerzon's names must be the book's
   names or their start followed by a mascot (`MEMPHIS` = "Memphis Tigers"; `ST`
   reads `STATE`; never "Michigan State" for `MICHIGAN` or "Texas Tech" for
   `TEXAS`), in their home/away slots, with kickoffs within 30 minutes and one
   such listing only; anything else is "this game is not listed". A game a book
   names differently (`UL LAFAYETTE` vs "Louisiana") is not priced there. On
   2026-11-12 BetMGM declined every college cell (its college 1H markets are
   suspended until near kickoff) and DraftKings posts 1H alternates only a day or
   two out, so many college corners price at one or two books.
4. **Wagerzon price**: every corner of a priced game is previewed with
   `ConfirmWagerHelper` RiskWin=2 (nothing placed) at $500, else $100
   (`parlay_pricer.get_combined_parlay_price`). Sizing uses the whole American
   odds at or under that decimal (`halves_board.conservative_american`).
5. **Stake** (`halves_board.size_half_board` -> `size_halves`): quarter Kelly on
   the page's bankroll. **NFL** parlays come **after the fectas**: their budget
   is the cash available less the fectas' recommended stakes (cash only; with
   Use free play on, fectas draw on the free play balance and parlays keep the
   cash). **College** parlays are sized **by slate** (Cal, 2026-11-13): a slate
   is the college games still to kick off whose first halves overlap its first
   game's (`cfb_slates`: kickoff less than `CFB_FIRST_HALF_WINDOW` = 2 h after
   the slate's first kickoff; not chained, so a day of 1-2 h gaps still
   splits), and each slate splits the whole cash available, because a 1H
   parlay settles at halftime, before the next slate. The only thing held back
   is the recommended stake of a fecta on an NFL game kicking off before the
   slate's halftime (never on a Saturday; Thursday/Friday night or late-season
   Saturday NFL). A later slate's stakes are only final once the earlier slate
   settles and Refresh reads the new balance, and the page's Recommended tile
   adds every slate up, so it can show more than the cash. An NFL parlay sits in one group with
   that spread team's fectas in the game (the fectas' `board._team_key`): the
   fectas' placed bets and recommended stakes are fixed bets the parlay is sized
   around, as one nested chain, the maximum correlation the fairs allow, so a
   parlay beside a fecta gets less, never more. A team's over and under corners
   cannot both win, which the chain also over-counts (cautious). Placed 1H
   parlays (`placed_half_parlays`, counting the ones that may have gone through)
   shrink the next stake on that side; a placed bet the board cannot price
   blocks its side. $20 minimum, $250 per parlay a week
   (`config.WZ_HALF_PARLAY_MAX_STAKE`, a mistyped-stake guard, not a measured
   Wagerzon limit).
6. **Place**: the server re-checks the corner is priced at the price the page
   showed, the game has not started, the stake is whole dollars, $20 to $250
   with the week's earlier tickets on that parlay, and within the cash
   available, and nothing went to Wagerzon on that parlay in the last 2 minutes.
   Then `wz_halves.place_half_parlay` re-reads Wagerzon's 1H lines (a changed
   number or price is `line_moved`), runs `parlay_placer`'s ConfirmWagerHelper
   preflight at the stake (a to-win under the board's price by more than
   Wagerzon's whole-dollar rounding plus the board price's own, $0.50 + 0.005 x
   stake, is `price_moved`; a failure before the submit is `rejected`) and sends ONE PostWagerMultipleHelper submit on the chosen
   account. It does not use `parlay_placer.place_parlays`, which re-submits
   after a logged-out reply to the submit. A submit that fails, answers logged
   out or cannot be read is `unknown` / `network_error`: it may be on Wagerzon,
   so it blocks a retry, counts as held and shows "?" under placed. Every
   attempt past the server's own checks is recorded (those refusals are not). A failed 1H refresh clears the 1H board, so its old fairs
   cannot be placed.

Refresh prices the fectas, then the 1H board (about 70 s for ~70 games on
2026-11-12), one after the other so each book's pacing holds.

## Files

| File | Job |
|---|---|
| `wz.py` | scrape the week's fectas; place one |
| `wz_halves.py` | read Wagerzon's NFL and college 1H lines, preview every pairing, place a 1H parlay |
| `halves_board.py` | one 1H refresh (book matching, 2 x 2 partition fairs); NFL sizing after the fectas, college by slate |
| `special_parser.py` | description -> legs (anything unrecognized is refused, never guessed) |
| `books.py` | adapter interface: a leg = an exhaustive outcome group |
| `fd_book.py`, `mgm_book.py`, `dk_book.py` | per-book market lookup + SGP price |
| `pricing.py` | partition devig, EV, Kelly |
| `board.py` | one refresh; read-time sizing |
| `store.py` | DuckDB state |
| `app.py`, `static/index.html` | loopback web page + JSON API |
| `dk_relay.py`, `com.nflwork.dk-relay.plist` | DraftKings-only CONNECT relay on the Mac for the VM, and its launchd agent |
| `requirements.txt` | third-party packages for the VM image (`unabated_ticket/deploy/nfl_specials.Dockerfile`) |

## Data: `nfl_specials/nfl_specials.duckdb`

| Table | Write | Row |
|---|---|---|
| `fecta_quotes` | APPEND every refresh | one (special, book) fair or the reason there is none |
| `placed_fectas` | APPEND every Place click | the attempt and Wagerzon's answer; `free_play` TRUE for a free play bet |
| `settings` | UPSERT | bankroll, Kelly fraction, use free play |
| `half_parlay_quotes` | APPEND every refresh | one (1H game, book, corner) fair, or Wagerzon's parlay price (`book = 'Wagerzon'`, `wz_decimal`), or the reason there is none |
| `placed_half_parlays` | APPEND every 1H Place past the server's checks | the attempt and Wagerzon's answer, `wz_decimal` the price quoted at the stake |

## Troubleshooting

- **DraftKings chip shows "DK SGP price ... returned HTTP 404"** — the request lost
  its `X-SportId` header, or DK moved its SGP widget's API. Re-read the widget
  (`dk-same-game-parlay/<version>/samegameparlayweb-external.js`, linked from the
  sportsbook page's `sameGameParlayConfig`) for `getYourbetRequestParameters`.
- **A superfecta says "no SGP market for '... scores first'"** at FD/BetMGM — expected;
  only DK prices the scores-first leg (FD/BetMGM price only the trifecta part, with DK down).
- **Nothing on the board** — Wagerzon posts the specials Thursday-ish; every league
  named `NFL WEEK <n> - SPECIALS` is read (two can be up at once).
- **On the VM the DraftKings chip shows a curl error** ("Failed to connect to ...
  port 8097", or a timeout): the relay is not reachable. Check the Mac is awake and
  on Tailscale, `launchctl print gui/$(id -u)/com.nflwork.dk-relay` shows it running,
  and `nfl_specials/dk_relay.log` has a start line. "HTTP 403" on the chip means the
  request reached DraftKings without the relay (`NFL_SPECIALS_DK_PROXY` unset); "CONNECT
  tunnel failed, response 403" means the relay refused the host: a new DraftKings host
  belongs in `DK_HOSTS`.
- **"Available: unknown"** — the Wagerzon balance read failed, so no stakes are
  recommended; click Refresh.

## Not done yet

- **A live scores-first share with DraftKings down.** The fallback uses the
  fixed 0.776; FanDuel's standalone first-to-score (independence, ~0.50)
  never bets.
- Correlation across games is ignored (cross-game NFL correlation measured ~0).

## 1H parlays on whole-number lines (`push_grid.py`)

A 1H spread or total on a whole number (spread 6, total 26) can push. Wagerzon
then pays the other leg alone at its straight price; both pushing refunds.
Pricing is books only (Cal, 2026-11-12):

1. At each book, price the 2 x 2 corner partition at the half-point lines
   either side (5.5 / 6.5 x 25.5 / 26.5: four partitions, or two with one
   whole-number leg), probit-devigged per partition.
2. `fit_push_grid` least-squares fits the cover/push/lose x over/push/under
   grid those bands imply. A both-push cell is left out when the spread and
   total differ in parity (a 1H margin and total are both even or both odd).
3. Refused (the book prices nothing for the game): a cell below -0.005
   (smaller negatives are rounding, clamped to 0), a spread push over 20% or a
   total push over 15% (NFL 1H 2019-2025 maxes: margin 3 at 12.3%, total 20 at
   8.6%), or partitions that disagree with the fit by over 0.01 (real BetMGM
   brackets fit within 0.001). One book alone can still move mass between cells
   no other partition pins, which only a second book catches.
4. Each corner's fair is the **effective** fair at Wagerzon's parlay decimal D:
   P(both win) + push return / D, where push return = P(spread wins, total
   pushes) x spread decimal + P(spread pushes, total wins) x total decimal +
   P(both push). Same EV as the real bet; Kelly on it is the smaller,
   cautious stake (the binary swings more than the real payout). Worst case
   across books as for half-point lines.
5. A whole-number game needs two books priced (`MIN_BOOKS_WHOLE_NUMBER`,
   `halves_board.py`); with one, it is listed under "not priced" with the book
   named, stakes nothing, and Place refuses it (as it refuses any corner with
   no fair).

`half_parlay_quotes.fair_prob` is the both-win chance and `push_return` the
push return, so the effective fair can be rebuilt from the table;
`placed_half_parlays.fair_prob` is the effective fair the stake was sized on.
A held whole-number parlay counts its pushes at the price it was placed at.

Unverified: the push rule (other leg at its straight price, both push
refunds) is as Wagerzon was described, not read off a graded ticket. Check
the first graded whole-number push against it.

Coverage (probe 2026-11-12): DraftKings posts 1H alt lines about 2 days out
(NFL) and 1 day out (college); BetMGM college 1H legs price only near
kickoff; FanDuel has no 1H alt totals in the SGP tab, so it prices only a
whole-number spread with a half-point total.
