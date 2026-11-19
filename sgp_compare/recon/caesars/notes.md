# Caesars NFL recon (2026-11-19 UTC, state nj, logged out)

Read only. Captured with `mlb_sgp.caesars_client.CaesarsClient` (node-minted AWS-WAF token).

## Schedule feed
- `GET /sb/v4/sports/americanfootball/tabs` lists the tabs. There is no "Games" tab; the default is `HIGHLIGHTS%7CFeatured`.
- Working path: `v4/sports/americanfootball/tabs/HIGHLIGHTS%7CFeatured`. Any unknown tab id (e.g. `SCHEDULE%7CGames`) silently returns the same Featured body, so a 200 does not prove a tab exists.
- NFL competition: name `NFL`, id `007d7c61-07a7-4e18-bb40-15104b6eac92`. NCAAF (`b7eda1b3-0170-4510-9616-1bce561d7aa7`) is in the same feed.
- Featured listed 12 NFL games. `SCHEDULE%7CTD%20Scorer%20%F0%9F%99%8C` listed 13. The extra game is Monday night (Bengals at Commanders), not Saints at Bears; see Follow-up.
- Tab events carry no markets. Teams in `event.metadata.homeTeamName` / `awayTeamName`, start in `startTime` (UTC), `started` flag as in MLB. `metadata.sportsbookFilterName` = "Week 11".
- `v4/sports`, `v4/sports/americanfootball` and `v3/sports/americanfootball/tabs` are 404.
- `tabs.json` is the Featured body trimmed to the NFL competition.

## Event picked
Indianapolis Colts at Houston Texans, id `faa98a96-ed6b-e0a9-94c7-71f707d9e752`, start 2026-11-20T01:15:00Z. `GET /sb/v4/events/{id}` saved as `event.json` (70 markets in `keyMarketGroups[].markets[]`).

## Market templates (templateName, market type, metadata.period, has market.line, count)
- `|Alt Receiving Yards|`, standard-market-template, none, no line, 12
- `|Alt Receptions|`, standard-market-template, none, no line, 12
- `|Alt Rushing Yards|`, standard-market-template, none, no line, 6
- `|Team Total Points|`, over-under, none, line, 4
- `|Alt Passing Touchdowns|`, standard-market-template, none, no line, 2
- `|Alt Passing Yards|`, standard-market-template, none, no line, 2
- `|Alt Rushing + Receiving Yards|`, standard-market-template, none, no line, 2
- `|Correct Score|`, standard-market-template, none, no line, 2
- `|Match Spread|`, two-way-handicap, MATCH_INC_OT, line, 2
- `|Team| |First Score Of Game Will Be|`, standard-market-template, none, no line, 2
- `|Total Points|`, over-under, MATCH_INC_OT, line, 2
- `|Winning Margin|`, standard-market-template, none, no line, 2
- `|1st Half Betting|` (3-way), `|1st Half Money Line|`, FIRST_HALF, no line, 1 each
- `|1st Half Spread|` (two-way-handicap), `|1st Half Total Points|` (over-under), FIRST_HALF, line, 1 each
- `|1st Quarter Betting|` (3-way), `|1st Quarter Money Line|`, FIRST_QUARTER, no line, 1 each
- `|1st Quarter Spread|`, `|1st Quarter Total Points|`, FIRST_QUARTER, line, 1 each
- `|Money Line|`, standard-market-template, MATCH_INC_OT, no line, 1
- 1 each, no period, no line: `|1st Scoring Play|`, `|Anytime Touchdown Scorer|`, `|Defensive or Special Team Touchdown Scorer|`, `|Field Goal| |Specials|`, `|First Touchdown Scorer|`, `|Odds Boost|`, `|Player To Score 2 or More Touchdowns|`, `|Player To Score 3 or More Touchdowns|`, `|Player| |Specials|`, `|Scoring| |Specials|`, `|Will There Be Overtime?|`

Names are pipe-wrapped (`|Spread|`, `|Houston Texans|`). Main Spread and Total Points each appear twice: one plain, one also carrying `movingLines.linePrices[]` (the alt ladder, each entry `{line, selections:[{selectionType, price}]}`).

## Selections
- `type`: spread/ML `home`/`away`; totals `over`/`under`; 3-way `draw`; specials `generic`/`none`. Counts: over 298, generic 131, none 101, home 56, away 52, under 8, draw 2.
- Spread `market.line` is the HOME line (-12.5 = Texans -12.5), same as the MLB run line; the away leg sends the negated line.
- Player yardage/reception props are over-only alt ladders with no `market.line`: market `|Nico Collins - Receiving Yards|`, selections typed `over` named `|30+|`, `|40+|` ... `|140+|` (not in order). The threshold is only in the name. Player, team and ids in `market.metadata` (`player`, `playerId`, `teamName`, `team` HOME/AWAY, `marketTypeCode` e.g. `player-alt-receiving-yards`).
- Anytime TD: market `|Player to Score a Touchdown|`, one selection per player named `|David Montgomery|`, `type` = the player's side (`home`/`away`). No "no" side.
- Team totals: `|Houston Texans Total Points|`, over-under with `line`, no period.

## Periods
`market.metadata.period`: `MATCH_INC_OT` (main spread, total, ML), `FIRST_HALF`, `FIRST_QUARTER`. Props carry no period. No 2H or Q2-Q4 markets in this event.

## Limits
- Event JSON: every market has `betLimits {liabilityLimit, maxMultipleBet, maxWin, maxBet}`, apparently cents (main spread maxBet 5,000,000 = $50k; receiving-yards ladder maxBet 150,000, maxWin 500,000; one market has maxMultipleBet 0).
- `/v2/bets/details` carries no stake limit fields (only `potentialReturns: 0`, `bettable`, `delay`).

## Combo pricing (`POST /sb/v2/bets/details`, empty `combinationSelections`)
The parlay is `parlays[0]` with `type: "BYO"`, `externalBetType: "SGP2"`, `typeName: "Same Game Parlay"`, price `{decimal, fractional, american}` (american is a string, "+245").

- `parlay_a_spread_total`: Texans -12.5 (1.9434) + Over 43.5 (1.8547). Product 3.604, SGP 3.45 (+245), ratio 0.957.
- `parlay_b_spread_recyds`: Texans -12.5 (1.9434) + Nico Collins 60+ rec yds (1.45872). Product 2.835, SGP 2.89 (+189), ratio 1.019.
- `parlay_c_ml_anytd`: Texans ML (1.13333) + David Montgomery anytime TD (1.7874). Product 2.026, SGP 1.97 (-103), ratio 0.972.
- `parlay_d1_spread_alone`, `parlay_d2_total_alone`: HTTP 200, `parlays: []`, the leg price only in `singles[]`.
- `parlay_e_invalid_over_under`: Over + Under 43.5. HTTP 200, top-level `errors: []` and `bettable: true`; `parlays[0]` has no `price`, `errors[0].details.message` "Two or more selections needed for an SGP", reason "(2) Cannot be combined", both legs `included: false`. `price_combo` already returns None here.

Each response file holds the exact request body, status and raw response.

Surprise: combo b priced above the independent product even though a Texans cover and Collins 60+ should correlate positively. Check more prop-x-side combos before trusting them.

## Follow-up (2026-11-19 ~02:55 UTC, same event, not started)

Files `parlay2_<name>.json` (request, status, raw response). "Product" = the leg decimals multiplied.

### Alt lines
- The alt ladders live on the SAME market as the main line: `|Spread|` id `a444e69d...` and `|Total Points|` id `92987c11...` are listed twice in `event.json`, once plain and once with `movingLines`. Spread ladder -24.5 .. -1.5 (home perspective), total ladder 32.5 .. 56.5, both in 1-point steps.
- An alt leg is the main selectionId + marketId with `line` set to the alt number. The response's `singles[].legs[].requestedLine` echoes the alt line; `singles[].line` still shows the market's main line (-12.5 / 43.5), so read `requestedLine`, not `line`.
- Caesars reprices server side: sending -7.5 with the -12.5 price (`parlay2_ctrl_*`) still returns the -7.5 ladder price 1.54 and the same parlay 2.79. The price we send is ignored.
- `alt_home_m7.5_alone`: single 1.54 (ladder price), no parlay.
- `alt_home_m7.5_over`: Texans -7.5 + Over 43.5 = 2.79 (+179), product 2.856, ratio 0.977. Main -12.5 + Over = 3.45 (`main_home_over`, unchanged from the first capture). So the alt line is priced, not mapped back to the main.
- `alt_away_p7.5_alone`: Colts +7.5 sent as the away selection with line +7.5 (the ladder's -7.5 entry, away side): single 2.33 with requestedLine 7.5.
- `alt_over_47.5_alone`: single 2.27, requestedLine 47.5. `alt_over_47.5_home_main`: + Texans -12.5 = 4.20 (+320), product 4.412, ratio 0.952.

### Other combos
| file | legs | product | SGP | ratio |
|---|---|---|---|---|
| hou_tt_over_home_main | Texans TT o28.5 + Texans -12.5 | 3.743 | 2.50 (+150) | 0.668 |
| h1_spread_total | Texans 1H -7 + 1H Over 22.5 | 3.729 | 3.10 (+210) | 0.831 |
| two_atd | Montgomery + Collins anytime TD | 3.474 | 3.70 (+270) | 1.065 |
| collins60_80 | Collins 60+ and 80+ | 3.049 | none | |
| ml_collins60 | Texans ML + Collins 60+ | 1.653 | 1.71 (-141) | 1.034 |
| away_main_collins60 | Colts +12.5 + Collins 60+ | 2.738 | 2.94 (+194) | 1.074 |
| home_main_collins60 | Texans -12.5 + Collins 60+ | 2.835 | 2.90 (+190) | 1.023 |

- Two thresholds of one player ladder are refused: `parlays[0]` has no price, reason "(1) Too closely related to selections in your parlay". Two players from the same anytime-TD market combine fine.
- 1H spread + 1H total combine (3.10). The 1H total main line was 22.5 here.
- Every combo with a player prop leg prices ABOVE the product (1.02 to 1.07), even same-team positive pairs (two Texans TD scorers, Texans cover + Collins). The direction between them is sensible (Colts +12.5 + Collins is longer than Texans -12.5 + Collins), so the SGP engine seems to use a fairer prop price than the posted ladder price; the posted prop price carries more margin than the SGP's. A ratio over 1 does not mean negative correlation here.
- Texans -12.5 + Collins 60+ moved 2.89 -> 2.90 between captures (prices drift).

### Tabs coverage
- Featured (`HIGHLIGHTS%7CFeatured`): 12 NFL events, Thu 2026-11-20 Colts at Texans through Sun night Vikings at 49ers. It DOES include Saints at Bears (`865f91aa-b6a5-ddd7-3a19-122d67e64d3f`); the first note was wrong.
- TD Scorer (`SCHEDULE%7CTD%20Scorer%20%F0%9F%99%8C`): the same 12 plus Monday night Bengals at Commanders (`c1c9f6ff-ceb2-df25-9304-81b21a9dab71`, 2026-11-24T01:15Z), 13 total.
- Neither tab is proven complete; TD Scorer covered the whole Thu-Mon week here and Featured stopped at Sunday. Reading both and taking the union of NFL event ids is the safe listing.
- Featured ids/starts: faa98a96 Colts@Texans 11-20 01:15Z; c8c55c7a Dolphins@Bills, b6cc2370 Cardinals@Chiefs, 347b047a Buccaneers@Lions, a464d90c Steelers@Eagles, 865f91aa Saints@Bears, 08fc0e68 Titans@Cowboys, a7321fbe Ravens@Panthers all 11-22 18:00Z; daed573a Jets@Chargers 21:05Z; c3a79427 Jaguars@Giants, ae99d872 Raiders@Broncos 21:25Z; 06c5b44d Vikings@49ers 11-23 01:20Z.
