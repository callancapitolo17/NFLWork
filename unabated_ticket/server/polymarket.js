// Unabated Ticket server — the Ticket's Place button at Polymarket US (the
// CFTC exchange, api.polymarket.us; Cal 2026-10-21, option A): which
// Polymarket US contract a board line buys and what buying it costs now.
// The runner reads the venue (GET /polymarket_quote.json); the bets service
// (bets_service/polymarket_place.py) checks the price and sends the order.
//
// Unabated's feed names no Polymarket US market (main lines and moneylines
// carry no rung id at any venue), so the contract is found in Polymarket
// US's own list, as Place finds Novig and Kalshi moneylines (Cal 2026-10-19:
// "found in each venue's own list, refuse on any mismatch"):
//   1. the game: Polymarket US's events starting within START_TOLERANCE_MS of
//      the board's start (bets.js's own game rule), filtered to the line's
//      market type; the one whose two teams resolve (teams.js + the team crosswalk, exactly as Cal's
//      Polymarket US bets are matched to the board) to the board game's away
//      and home, in that order;
//   2. the market: among that event's markets of the type, the one market and
//      contract side whose own structured sides read as the row's side and
//      number (the bets source's grammar, sources/polymarket_us.py:
//      YES is the `long` side; a spread side's description is its own signed
//      number; a total's YES is Over). Zero or two of either refuses.
// A first-half (or any non-full-game) moneyline refuses: whether Polymarket
// US voids a tied half or settles it at $0.50 is unchecked, and a contract
// bought at Unabated's tie-no-bet fair that pays half on a tie is mispriced
// by about P(tie) / 2 (the same reason Place refuses them at Kalshi and Novig).
//
// Prices: Polymarket US trades ONE instrument per market, the YES contract;
// its book's bids and offers are YES prices. YES costs the best offer; NO
// costs 1 minus the best bid. Takers pay THETA x p x (1 - p) per contract
// (the exchange's fee curve, CFTC filings; allInOf) — counted in the all-in
// price the bets service holds against the Ticket's.
//
// Pure: no network, no clock. Reused by runner.js only.

const BOOK_NAME = "Polymarket US";
const VENUE = "polymarket_us";
const BET_TYPE_MONEYLINE = 1;
const BET_TYPE_SPREAD = 2;
const BET_TYPE_TOTAL = 3;
const PERIOD_FULL_GAME = 1;
const SIDE_AWAY_OR_OVER = 0;
// The board's league path -> Polymarket US's sport word in sportsMarketType.
const SPORT_OF_LEAGUE = { nfl: "football", cfb: "football", nba: "basketball", cbb: "basketball", wnba: "basketball", mlb: "baseball", nhl: "hockey" };
// The board's period (feed.js PERIODS) -> Polymarket US's period words.
const PERIOD_WORDS = { FG: "full_game", "1H": "first_half", "2H": "second_half", "1Q": "first_quarter", "2Q": "second_quarter", "3Q": "third_quarter", "4Q": "fourth_quarter" };
const KIND_OF_BET_TYPE = { [BET_TYPE_MONEYLINE]: "winner", [BET_TYPE_SPREAD]: "spread", [BET_TYPE_TOTAL]: "total" };
// Polymarket US's team field that spells a team the way Unabated does
// (sources/polymarket_us.py TEAM_NAME_FIELDS): safeName for college and the
// NHL, name elsewhere.
const TEAM_NAME_FIELDS = { cfb: "safeName", cbb: "safeName", nhl: "safeName" };
// Polymarket US's events are asked for, and must start, within this of the
// board's start (bets.js's own game rule).
const START_TOLERANCE_MS = 30 * 60 * 1000;
// Events per page asked for; a page this full may be cut, so the quote refuses
// rather than miss the game or a second copy of it.
const EVENTS_PAGE_LIMIT = 100;
const MARKET_STATE_OPEN = "MARKET_STATE_OPEN";
// Polymarket US's taker fee coefficient: fee per contract = THETA x p x (1 - p).
// 0.06 is the exchange-wide coefficient of the July 2026 CFTC filing; a
// sports-only 0.03 is reported second-hand, so the higher one is used:
// overcounting the fee can only refuse a click, never overpay. The order's
// own commission is read back and stored, so the first real fills show which.
const TAKER_FEE_THETA = 0.06;

function isPolymarketUsBook(book) {
  return Boolean(book) && book.name === BOOK_NAME;
}

// "football_team_full_game_spread" for a board row, or {reason, detail}.
function marketTypeOf(row) {
  const sport = SPORT_OF_LEAGUE[row.league];
  const period = PERIOD_WORDS[row.period];
  const kind = KIND_OF_BET_TYPE[row.betTypeId];
  if (!sport) return { reason: "no_market", detail: `Polymarket US Place does not cover league ${row.league}` };
  if (!period || !kind) return { reason: "no_market", detail: `no Polymarket US market type for ${row.period} bet type ${row.betTypeId}` };
  if (row.betTypeId === BET_TYPE_MONEYLINE && row.periodTypeId !== PERIOD_FULL_GAME) {
    return { reason: "no_market", detail: `a ${row.period} moneyline: whether Polymarket US voids a tied period is unchecked` };
  }
  return { type: `${sport}_team_${period}_${kind}` };
}

function teamNameOf(team, league) {
  const name = team ? team[TEAM_NAME_FIELDS[league] || "name"] : null;
  return typeof name === "string" && name.trim() ? name : null;
}

// {awayTeam, homeTeam, awayId, homeId} of an event in the row's league, away
// first (Polymarket US lists away, home), or null.
function eventTeamsOf(event, league) {
  const teams = event && Array.isArray(event.teams) ? event.teams : null;
  if (!teams || teams.length !== 2 || !teams.every((team) => team && typeof team === "object")) return null;
  if (!teams.every((team) => team.league === league)) return null;
  const [away, home] = teams;
  const awayTeam = teamNameOf(away, league);
  const homeTeam = teamNameOf(home, league);
  if (!awayTeam || !homeTeam) return null;
  return { awayTeam, homeTeam, awayId: away.id, homeId: home.id };
}

// The events that are the row's game: same away and home team (after the
// crosswalk) and a start within START_TOLERANCE_MS. Returns them with their teams.
//   keysOf     (record) -> {awayKey, homeKey}: bets.resolveTeamKeys over the held crosswalk
//   boardKeys  {away, home}: the row's own team keys (teams.keyOf of its team ids)
function matchingEvents(events, row, keysOf, boardKeys) {
  const startMs = Date.parse(row.eventStart);
  const found = [];
  for (const event of Array.isArray(events) ? events : []) {
    const teams = eventTeamsOf(event, row.league);
    if (!teams) continue;
    const eventStartMs = Date.parse(event.startTime);
    if (!Number.isFinite(eventStartMs) || !Number.isFinite(startMs) || Math.abs(eventStartMs - startMs) > START_TOLERANCE_MS) continue;
    const keys = keysOf({ venue: VENUE, league: row.league, awayTeam: teams.awayTeam, homeTeam: teams.homeTeam });
    if (!keys.awayKey || !keys.homeKey) continue;
    if (keys.awayKey === boardKeys.away && keys.homeKey === boardKeys.home) found.push({ event, teams });
  }
  return found;
}

function signedNumber(text) {
  if (text == null || String(text).trim() === "") return null;
  const value = Number(String(text).trim());
  return Number.isFinite(value) ? value : null;
}

// [held side, other side] of a two-sided market for `contractSide`, or null.
function sidesOf(market, contractSide) {
  const sides = Array.isArray(market.marketSides) ? market.marketSides : null;
  if (!sides || sides.length !== 2) return null;
  const yes = sides.filter((side) => side && side.long === true);
  const no = sides.filter((side) => side && side.long !== true);
  if (yes.length !== 1 || no.length !== 1) return null;
  return contractSide === "yes" ? [yes[0], no[0]] : [no[0], yes[0]];
}

function teamIdOf(side) {
  if (side.team && typeof side.team === "object" && side.team.id != null) return side.team.id;
  return side.teamId ?? null;
}

function slotOf(side, teams) {
  const id = teamIdOf(side);
  if (id != null && id === teams.awayId) return "away";
  if (id != null && id === teams.homeId) return "home";
  return null;
}

// {side: away|home|over|under, points} that buying `contractSide` of the
// market is, or null when the market does not read cleanly (the bets source's
// side_and_points_of, so a Place bet reads back as the line it was placed on).
function sideAndPointsOf(market, contractSide, betTypeId, teams) {
  const sides = sidesOf(market, contractSide);
  if (!sides) return null;
  const [held, other] = sides;
  if (betTypeId === BET_TYPE_TOTAL) {
    if (sides.some((side) => side.team || side.teamId)) return null;
    const word = String(held.description || "").trim().toLowerCase();
    const line = signedNumber(market.line);
    if ((word !== "over" && word !== "under") || line == null || line <= 0) return null;
    return { side: word, points: line };
  }
  const heldSlot = slotOf(held, teams);
  const otherSlot = slotOf(other, teams);
  if (!heldSlot || !otherSlot || heldSlot === otherSlot) return null;
  for (const side of sides) {
    const ordering = side.team && typeof side.team === "object" ? side.team.ordering : null;
    if ((ordering === "away" || ordering === "home") && ordering !== slotOf(side, teams)) return null;
  }
  if (betTypeId === BET_TYPE_MONEYLINE) return { side: heldSlot, points: null };
  const points = signedNumber(held.description);
  const otherPoints = signedNumber(other.description);
  const line = signedNumber(market.line);
  if (points == null || otherPoints == null || line == null || points !== -otherPoints || Math.abs(points) !== Math.abs(line)) return null;
  return { side: heldSlot, points };
}

// The side a row is on in the bets' words.
function rowSideOf(row) {
  if (row.betTypeId === BET_TYPE_TOTAL) return row.sideIndex === SIDE_AWAY_OR_OVER ? "over" : "under";
  return row.sideIndex === SIDE_AWAY_OR_OVER ? "away" : "home";
}

// The one contract on the event that is the row: {slug, side: yes|no, market},
// or {reason, detail}.
//   type  marketTypeOf(row).type
function contractOf(row, event, teams, type) {
  const markets = (Array.isArray(event.markets) ? event.markets : []).filter((market) => market && market.sportsMarketType === type);
  const wanted = rowSideOf(row);
  const found = [];
  for (const market of markets) {
    for (const contractSide of ["yes", "no"]) {
      const taken = sideAndPointsOf(market, contractSide, row.betTypeId, teams);
      if (!taken || taken.side !== wanted) continue;
      if (row.betTypeId !== BET_TYPE_MONEYLINE && taken.points !== row.points) continue;
      found.push({ slug: market.slug, side: contractSide, market });
    }
  }
  const label = `${wanted}${row.points == null ? "" : ` ${row.points}`}`;
  if (found.length === 0) return { reason: "no_market", detail: `no ${type} market on ${event.slug} reads as ${label} (${markets.length} listed)` };
  if (found.length > 1) return { reason: "no_market", detail: `${found.length} ${type} markets on ${event.slug} read as ${label}: ${found.map((entry) => `${entry.slug} ${entry.side}`).join(", ")}` };
  if (typeof found[0].slug !== "string" || !found[0].slug) return { reason: "no_market", detail: `the ${label} market on ${event.slug} has no slug` };
  return found[0];
}

function amountOf(amount) {
  const value = amount && typeof amount === "object" ? amount.value : amount;
  const number = Number(value);
  return value != null && value !== "" && Number.isFinite(number) ? number : null;
}

// What buying `side` costs now off GET /v1/markets/<slug>/book: {ask, depth}
// (ask = that side's own price, depth = contracts at it), or {reason, detail}.
function quoteOf(bookBody, side) {
  const book = bookBody && typeof bookBody === "object" ? bookBody.marketData : null;
  if (!book || typeof book !== "object") return { reason: "no_book", detail: "Polymarket US's book has no marketData" };
  if (book.state !== MARKET_STATE_OPEN) return { reason: "market_closed", detail: `the market is ${book.state || "in no state"}, not open` };
  const levels = (list) => (Array.isArray(list) ? list : [])
    .map((level) => ({ price: amountOf(level && level.px), qty: Number(level && level.qty) }))
    .filter((level) => level.price != null && level.price > 0 && level.price < 1 && Number.isFinite(level.qty) && level.qty > 0);
  if (side === "yes") {
    const offers = levels(book.offers);
    if (!offers.length) return { reason: "no_book", detail: "nothing is offered on YES" };
    const best = Math.min(...offers.map((level) => level.price));
    return { ask: best, depth: offers.filter((level) => level.price === best).reduce((sum, level) => sum + level.qty, 0) };
  }
  const bids = levels(book.bids);
  if (!bids.length) return { reason: "no_book", detail: "nothing is bid on YES, so nothing sells NO" };
  const best = Math.max(...bids.map((level) => level.price));
  return { ask: Math.round((1 - best) * 1e6) / 1e6, depth: bids.filter((level) => level.price === best).reduce((sum, level) => sum + level.qty, 0) };
}

// One contract's price with the taker fee, as a probability.
function allInOf(ask) {
  return Math.round((ask + TAKER_FEE_THETA * ask * (1 - ask)) * 1e6) / 1e6;
}

module.exports = {
  BOOK_NAME, VENUE, START_TOLERANCE_MS, EVENTS_PAGE_LIMIT, TAKER_FEE_THETA, MARKET_STATE_OPEN,
  isPolymarketUsBook, marketTypeOf, eventTeamsOf, matchingEvents, sideAndPointsOf, rowSideOf, contractOf, quoteOf, allInOf,
};
