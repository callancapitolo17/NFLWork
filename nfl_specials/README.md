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
  VM's `deploy/.env`). **Mac asleep or off: the DraftKings chip errors and
  superfectas stay unpriced; trifectas still price.**
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
  one DK fair), WZ
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
| `dk_relay.py`, `com.nflwork.dk-relay.plist` | DraftKings-only CONNECT relay on the Mac for the VM, and its launchd agent |
| `requirements.txt` | third-party packages for the VM image (`unabated_ticket/deploy/nfl_specials.Dockerfile`) |

## Data: `nfl_specials/nfl_specials.duckdb`

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

- **Superfectas with the Mac off.** DraftKings needs a home IP; nothing at home
  but the Mac carries the relay.
- Correlation across games is ignored (cross-game NFL correlation measured ~0).
