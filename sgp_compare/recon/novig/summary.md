# Novig NFL recon (2026-11-19, read-only)

Game: IND @ HOU, 2026-11-20 01:15 UTC (soonest OPEN_PREGAME card on the NFL page; 29 game cards on the page).

## Parlay endpoint moved
- `POST https://api.novig.us/nbx/v1/parlay/request/unauthenticated` (what `mlb_sgp/scraper_novig_sgp.py` and `novig_client.py` still call) now answers **404 "Cannot POST"** (`legacy_v1_parlay_404.json`). Every existing Novig SGP caller is broken.
- The novig.com bundle's `requestAnonymousQuote` now calls **`POST https://api.novig.us/rfq/v2/quote/unauthenticated`** with body `{"outcomeIds": ["<outcome uuid>", ...]}` (base `RFQ_API_BASE = https://api.novig.us/rfq`, 10 s client timeout).
- Reply: `{"price": "0.271", "expiresAt": <ms>, "serverNow": <ms>}`. One price only: **no offers list, no per-leg prices, no `vendor`**. `price` = cost per $1 payout (implied prob, vig included); decimal = 1/price. Quote lives 15-30 s.
- Decline: **HTTP 400, plain text `Invalid parlay combination.`** (both sides of one total).
- The app wraps the reply as `maxCollateral = 250000 * price` (`USER_MARKET_MAXIMUM = 25e4`); execute is a separate authenticated `/execute`, never called here.

## Parlay results (`parlay_<name>.json` = request + raw response; `parlay_results.txt` = one line each)
| case | legs | HTTP | price | price / product of legs' `available` |
|---|---|---|---|---|
| a | HOU -12.5 + o43.5 | 200 | 0.271 | 1.053 |
| b | HOU ML + o43.5 | 200 | 0.448 | 1.014 |
| c | HOU -12.5 + Nico Collins o70.5 rec yds | 200 | 0.262 | 1.028 |
| d | Collins o70.5 rec + Stroud o257.5 pass yds | 200 | 0.343 | 1.269 |
| e | J. Taylor anytime TD + HOU -12.5 | 200 | 0.154 | 0.724 |
| f | o21.5 1H + HOU -12.5 | 200 | 0.288 | 1.108 |
| f2 | HOU -6.5 1H + o43.5 | 200 | 0.295 | 1.031 |
| g | o43.5 + u43.5 | 400 | Invalid parlay combination. | - |
| h | spread + total + rec yds + anytime TD | 200 | 0.072 | 1.263 |
| i | home/over (= case a), home/under, away/over, away/under | 200 x4 | 0.271 / 0.225 / 0.251 / 0.260 | sum **1.007** |
| j | Collins o70.5 rec yds PROTECTED + HOU -12.5 | 200 | 0.298 | 1.066 |

- Every leg type tried was accepted: spread, ML, total, 1H spread/total, player yards, anytime TD, PROTECTED prop, 4 legs. Only the same-market opposite pair declined.
- The (a) four-corner partition sums to 1.007, so the SGP overround is ~0.7% on a 2x2, about the singles' 0.5% per market. Correlation signs look sane (e: an IND RB TD against HOU -12.5 is priced well under independence).
- Vendor mix: unknowable from v2; the per-leg `vendor` the MLB Novig-relay guard (`router.drop_relayed_novig`) depends on is gone.

## Markets (`markets_ind_hou.json`, `market_types.txt`)
- 927 markets, 32 types, all OPEN. Mains: MONEY 1, SPREAD 28 (alt ladder), TOTAL 25, TEAM_TOTAL 4, SPREAD_1H 2, TOTAL_1H 3, MONEY_1H 1. **No quarter markets.**
- Props (player on `market.player`, incl. `player_competitors[].competitor.symbol` for the team): RECEIVING_YARDS 188, RUSHING_YARDS 87, RECEPTIONS 80, TOUCHDOWNS 38 (Over 0.5 = anytime TD), PASSING_YARDS 37, LONGEST_RECEPTION 31, FIRST_TOUCHDOWN_SCORER 30 (Yes/No), PASSING_AND_RUSHING_YARDS 24, RUSHING_AND_RECEIVING_YARDS 16, LONGEST_RUSH 13, PASSING_TOUCHDOWNS 8, RUSHING_ATTEMPTS 7, PASSING_ATTEMPTS 6, PASSING_COMPLETIONS 4, LONGEST_COMPLETION 3, FIELD_GOALS_MADE 3, INTERCEPTIONS_THROWN 2, KICKING_POINTS 2, TACKLES_ASSISTS 1, plus `_PROTECTED` twins of rec/rush/receptions/TD/pass yds/pass TD (separate outcome ids, priced a bit higher, e.g. Collins o70.5 0.565 vs 0.515; meaning not confirmed, likely void-if-inactive protection).
- Every market has `strike` (MONEY's is 0.0). `market.competitor` is set only on TEAM_TOTAL; spread/ML team is on `outcome.competitor`. Outcomes are always 2 (Over/Under, team/team, Yes/No).
- Price fields: `outcome.available` = best price to buy that outcome now (implied prob; main pairs sum 1.005), **null when nothing rests** (about 10% of prop outcomes, e.g. most anytime-TD Unders). `outcome.altAvailable` null on all 1,854 outcomes. `market.is_consensus` = the main line per (type, player) (1 per game market, 1-2 per player per prop type, 30/38 TD). `market.volume` = traded dollars (nonzero on 274/927, max 615k). `market_locks` non-empty on 9 (a lock id + created_at; e.g. Will Mallory TD, whose outcomes then have null `available`). `market_detail` always null. `tournament_markets` (1,082) is just a list of `{type}`; `scoreboard_markets` is the settlement MONEY market with outcome status TBD.

## Transport
- NFL page and EventMarkets_Query both 200 with the committed query (no refresh needed). ~17 parlay quotes at ~1 s spacing: no 429/403, no Cloudflare challenge.
