// Run: node --test unabated_ticket/tests/polymarket.test.js
// server/polymarket.js: which Polymarket US contract a board line buys, and its price.
const test = require("node:test");
const assert = require("node:assert/strict");
const polymarket = require("../server/polymarket.js");

const TEAMS = { awayTeam: "Duquesne", homeTeam: "Washington State", awayId: 501, homeId: 502 };
const DUQ = { id: 501, safeName: "Duquesne", league: "cfb", ordering: "away" };
const WST = { id: 502, safeName: "Washington State", league: "cfb", ordering: "home" };
const ROW = { league: "cfb", period: "FG", periodTypeId: 1, betTypeId: 2, sideIndex: 1, points: -35.5 };

test("market types: league, period and bet type; a half moneyline refuses (tie rule unchecked)", () => {
  assert.equal(polymarket.marketTypeOf(ROW).type, "football_team_full_game_spread");
  assert.equal(polymarket.marketTypeOf({ ...ROW, league: "nba", betTypeId: 3, period: "1H", periodTypeId: 2 }).type, "basketball_team_first_half_total");
  assert.equal(polymarket.marketTypeOf({ ...ROW, betTypeId: 1 }).type, "football_team_full_game_winner");
  assert.match(polymarket.marketTypeOf({ ...ROW, betTypeId: 1, period: "1H", periodTypeId: 2 }).detail, /1H moneyline/);
  assert.match(polymarket.marketTypeOf({ ...ROW, league: "epl" }).detail, /league epl/);
});

test("a total's Over is YES and its Under NO; a team total never matches", () => {
  const total = { slug: "t", sportsMarketType: "football_team_full_game_total", line: 52.5,
    marketSides: [{ long: true, description: "Over" }, { long: false, description: "Under" }] };
  assert.deepEqual(polymarket.sideAndPointsOf(total, "yes", 3, TEAMS), { side: "over", points: 52.5 });
  assert.deepEqual(polymarket.sideAndPointsOf(total, "no", 3, TEAMS), { side: "under", points: 52.5 });
  const teamTotal = { ...total, marketSides: total.marketSides.map((side) => ({ ...side, teamId: 501 })) };
  assert.equal(polymarket.sideAndPointsOf(teamTotal, "yes", 3, TEAMS), null);
  const event = { slug: "e", markets: [total] };
  const under = polymarket.contractOf({ ...ROW, betTypeId: 3, sideIndex: 1, points: 52.5 }, event, TEAMS, total.sportsMarketType);
  assert.deepEqual([under.slug, under.side], ["t", "no"]);
});

test("a moneyline is the side whose team the row backs; an away/home that disagrees with the event fails closed", () => {
  const winner = { slug: "w", sportsMarketType: "football_team_full_game_winner",
    marketSides: [{ long: true, teamId: 502, team: WST }, { long: false, teamId: 501, team: DUQ }] };
  const event = { slug: "e", markets: [winner] };
  const away = polymarket.contractOf({ ...ROW, betTypeId: 1, sideIndex: 0, points: null }, event, TEAMS, winner.sportsMarketType);
  assert.deepEqual([away.slug, away.side], ["w", "no"]);
  const flipped = { ...winner, marketSides: [{ long: true, teamId: 502, team: { ...WST, ordering: "away" } }, winner.marketSides[1]] };
  assert.equal(polymarket.sideAndPointsOf(flipped, "yes", 1, TEAMS), null);
});

test("YES costs the best offer, NO 1 minus the best bid; the all-in adds theta p (1 - p)", () => {
  const book = { marketData: { state: "MARKET_STATE_OPEN",
    bids: [{ px: { value: "0.4100" }, qty: "10" }, { px: { value: "0.4150" }, qty: "7" }],
    offers: [{ px: { value: "0.4300" }, qty: "3" }, { px: { value: "0.4300" }, qty: "2" }, { px: { value: "0.4500" }, qty: "9" }] } };
  assert.deepEqual(polymarket.quoteOf(book, "yes"), { ask: 0.43, depth: 5 });
  assert.deepEqual(polymarket.quoteOf(book, "no"), { ask: 0.585, depth: 7 });
  assert.equal(polymarket.quoteOf({ marketData: { ...book.marketData, offers: [] } }, "yes").reason, "no_book");
  assert.equal(polymarket.quoteOf({}, "yes").reason, "no_book");
  assert.equal(polymarket.allInOf(0.5), 0.515);
});
