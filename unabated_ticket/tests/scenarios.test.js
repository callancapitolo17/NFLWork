// Run: node --test unabated_ticket/tests
// server/scenarios.js: the Bet Tracker's Live tab cards — a game in progress
// cut into the results that change the money, each with its P&L and its
// chance off the kickoff fair ladder. Boards are built by hand in the feed's
// shape (as teaser.test.js does); no network.
const test = require("node:test");
const assert = require("node:assert/strict");
const feed = require("../extension/feed.js");
const ladderLib = require("../extension/ladder.js");
const scenarios = require("../server/scenarios.js");

const NFL = 1;
const MLB = 5;
const FULL_GAME = 1;
const MONEYLINE = 1;
const SPREAD = 2;
const TOTAL = 3;
const AWAY_OR_OVER = 0;
const HOME_OR_UNDER = 1;
const BOOK = 59;
const FAIR_BOOK = 4;
const HOUR = 3600 * 1000;
const NOW = Date.parse("2026-10-05T19:00:00Z");
const KICKOFF = NOW - HOUR;

const nearly = (actual, expected, tol = 1e-9) =>
  assert.ok(Math.abs(actual - expected) <= tol, `expected ${expected}, got ${actual}`);

function probOf(american) {
  return american > 0 ? 100 / (american + 100) : -american / (-american + 100);
}

// One game: a book's main spread/total and the fair rungs [cut, bacr of the
// side that wins above the cut] (margin = away minus home).
function addGame(state, game) {
  const { eventId, away, home } = game;
  const leagueId = game.leagueId ?? NFL;
  const awayTeamId = eventId * 10 + 1;
  const homeTeamId = eventId * 10 + 2;
  state.teams[awayTeamId] = away;
  state.teams[homeTeamId] = home;
  state.events[eventId] = {
    eventId, leagueId, eventName: `${away} @ ${home}`, eventStart: game.startMs ?? KICKOFF,
    awayTeamId, homeTeamId, awayRotation: game.awayRotation, homeRotation: game.homeRotation,
    venueIds: { kalshiEventSuffixes: [], kalshiContracts: {}, novigOutcomes: {} },
  };
  const put = (fields) => {
    const line = {
      isAlt: false, leagueId, periodTypeId: FULL_GAME, eventId, marketId: eventId * 100 + fields.betTypeId, price: -110,
      sourceFormat: 1, sourcePrice: null, bacr: null, ge: null, liquidity: null, statusId: 1, sequenceNumber: null,
      isBlurred: false, modifiedOn: null, ...fields,
    };
    line.sideKey = `si${line.sideIndex}`;
    line.key = `${line.bookId}:${eventId}:${line.betTypeId}:${line.sideIndex}:${line.points}`;
    state.lines[line.key] = line;
  };
  put({ bookId: BOOK, betTypeId: SPREAD, sideIndex: AWAY_OR_OVER, points: game.spread ?? 2.5 });
  put({ bookId: BOOK, betTypeId: SPREAD, sideIndex: HOME_OR_UNDER, points: -(game.spread ?? 2.5) || 0 });
  for (const [cut, bacr] of game.marginRungs || []) put({ bookId: FAIR_BOOK, isAlt: true, betTypeId: SPREAD, sideIndex: AWAY_OR_OVER, points: -cut || 0, bacr });
  for (const [cut, bacr] of game.totalRungs || []) put({ bookId: FAIR_BOOK, isAlt: true, betTypeId: TOTAL, sideIndex: AWAY_OR_OVER, points: cut, bacr });
  for (const [sideIndex, bacr] of game.moneyline || []) put({ bookId: FAIR_BOOK, betTypeId: MONEYLINE, sideIndex, points: null, bacr });
  return state;
}

function emptyBoard() {
  return { teams: {}, teamIndex: {}, books: { [BOOK]: { id: BOOK, name: "Buckeye", isLive: true }, [FAIR_BOOK]: { id: FAIR_BOOK, name: "BetMGM", isLive: true } }, events: {}, lines: {} };
}

// The runner's remembered games: one describeLine row per event and its ladder.
function gamesOf(state, { noOdds } = {}) {
  const linesByEvent = ladderLib.groupLinesByEvent(Object.values(state.lines));
  const seen = new Set();
  const games = [];
  for (const line of Object.values(state.lines)) {
    if (line.isAlt || seen.has(line.eventId)) continue;
    seen.add(line.eventId);
    const lines = linesByEvent.get(line.eventId);
    const ladderOf = noOdds ? null : (period, axis) => ladderLib.buildLadder(lines, { periodTypeId: ladderLib.periodTypeIdOf(period), axis });
    games.push({ row: feed.describeLine(line, state), ladderOf, oddsAt: noOdds ? null : KICKOFF - 60 * 1000 });
  }
  return games;
}

// An open straight bet as the bets service serves it, joined by rotation and start.
function straightBet(id, fields) {
  return {
    id: `novig:${id}`, source: "novig_api", venue: "novig", status: "open", league: "nfl", betType: "spread", period: "FG",
    side: "away", points: 3.5, price: -110, rotation: 401, awayTeam: "BILLS", homeTeam: "CHIEFS", awayKey: null, homeKey: null,
    eventStart: new Date(KICKOFF).toISOString(), stake: 220, toWin: 200, isParlayLeg: false, parlayId: null,
    placedAt: "2026-10-05T12:00:00Z", approx: [], unmatchable: null, raw: {}, ...fields,
  };
}

function teaserLeg(parlayId, legIndex, fields) {
  return straightBet(`${parlayId}:leg${legIndex}`, {
    venue: "bfa", source: "bfa_api", stake: 200, toWin: 600, isParlayLeg: true, parlayId: `bfa:${parlayId}`, legIndex, legCount: 4,
    raw: { headerDescription: "4 TEAM TEASERS" }, ...fields,
  });
}

// Bills @ Chiefs, kicked off an hour ago, with rungs at every cut the bets need.
function billsChiefs() {
  return addGame(emptyBoard(), {
    eventId: 101, away: "Buffalo Bills", home: "Kansas City Chiefs", awayRotation: 401, homeRotation: 402,
    marginRungs: [[-9.5, -560], [-3.5, -150], [-1.5, 120], [-0.5, 125], [2.5, 200]],
    totalRungs: [[47.5, 104]],
  });
}

function bandRows(group) {
  return group.bands.map((band) => [band.label, band.pnl]);
}

test("a middle: Bills +3.5 and Chiefs -1.5 cut the margin into three results, each scored and priced", () => {
  const records = [
    straightBet("bills"),
    straightBet("chiefs", { side: "home", points: -1.5, rotation: 402, stake: 200, toWin: 210 }),
  ];
  const { games, coveredBetIds } = scenarios.buildScenarios({ records, games: gamesOf(billsChiefs()), now: NOW });
  assert.equal(games.length, 1);
  const [card] = games;
  assert.deepEqual([card.awayTeam, card.homeTeam, card.startMs], ["Buffalo Bills", "Kansas City Chiefs", KICKOFF]);
  assert.deepEqual(coveredBetIds, ["novig:bills", "novig:chiefs"]);
  assert.equal(card.groups.length, 1);
  const [group] = card.groups;
  assert.equal(group.title, "Result");
  assert.deepEqual(bandRows(group), [
    ["Kansas City Chiefs by 4+", -10],
    ["Kansas City Chiefs by 2-3", 410],
    ["Kansas City Chiefs by 1, tie or Buffalo Bills win", 0],
  ]);
  const above = (american) => probOf(american);
  nearly(group.bands[0].prob, 1 - above(-150));
  nearly(group.bands[1].prob, above(-150) - above(120));
  nearly(group.bands[2].prob, above(120));
  nearly(group.ev, group.bands.reduce((sum, band) => sum + band.prob * band.pnl, 0));
  assert.deepEqual([group.best, group.worst], [410, -10]);
});

test("a teaser leg adds its own cut and says where it dies; the ticket's dollars stay out of the P&L", () => {
  const records = [straightBet("bills"), teaserLeg("t1", 0, { points: 9.5 })];
  const [card] = scenarios.buildScenarios({ records, games: gamesOf(billsChiefs()), now: NOW }).games;
  const [group] = card.groups;
  assert.deepEqual(bandRows(group), [
    ["Kansas City Chiefs by 10+", -220],
    ["Kansas City Chiefs by 4-9", -220],
    ["Kansas City Chiefs by 3 or less, tie or Buffalo Bills win", 200],
  ]);
  assert.deepEqual(group.bands.map((band) => band.legs.map((leg) => leg.result)), [["lost"], ["won"], ["won"]]);
  assert.deepEqual(card.legs.map((leg) => [leg.label, leg.kind, leg.ticketStake, leg.ticketToWin, leg.group]), [["BILLS +9.5", "teaser", 200, 600, "FG:margin"]]);
});

test("a whole number gets its own push band; a total is its own group", () => {
  const state = addGame(billsChiefs(), { eventId: 101, away: "Buffalo Bills", home: "Kansas City Chiefs", awayRotation: 401, homeRotation: 402, marginRungs: [[-2.5, -140]] });
  const records = [
    straightBet("bills3", { points: 3 }),
    straightBet("over", { betType: "total", side: "over", points: 47.5, stake: 162, toWin: 150 }),
  ];
  const [card] = scenarios.buildScenarios({ records, games: gamesOf(state), now: NOW }).games;
  assert.deepEqual(card.groups.map((group) => group.title), ["Result", "Total points"]);
  assert.deepEqual(bandRows(card.groups[0]), [
    ["Kansas City Chiefs by 4+", -220],
    ["Kansas City Chiefs by 3", 0],
    ["Kansas City Chiefs by 2 or less, tie or Buffalo Bills win", 200],
  ]);
  nearly(card.groups[0].bands[1].prob, probOf(-150) - probOf(-140));
  assert.deepEqual(bandRows(card.groups[1]), [["47 or fewer points", -162], ["48+ points", 150]]);
});

test("a game the runner never saw before kickoff, or a missing rung, gets the card without chances", () => {
  const records = [straightBet("bills"), straightBet("deep", { points: 13.5 })];
  const noOdds = scenarios.buildScenarios({ records, games: gamesOf(billsChiefs(), { noOdds: true }), now: NOW }).games[0];
  assert.ok(noOdds.groups[0].bands.every((band) => band.prob === null));
  assert.equal(noOdds.groups[0].ev, null);
  const noRung = scenarios.buildScenarios({ records, games: gamesOf(billsChiefs()), now: NOW }).games[0];
  assert.deepEqual(noRung.groups[0].bands.map((band) => band.prob === null), [true, true, false]);
});

test("only started games get cards, and only open bets count; a bet off the ladder is named on the card", () => {
  const state = billsChiefs();
  addGame(state, { eventId: 102, away: "Dallas Cowboys", home: "New York Giants", awayRotation: 403, homeRotation: 404, startMs: NOW + HOUR, marginRungs: [[-3.5, -120]] });
  const records = [
    straightBet("bills"),
    straightBet("settled", { status: "won" }),
    straightBet("first-half", { period: "1H", stake: 50, toWin: 45 }),
    straightBet("no-stake", { stake: null, toWin: null }),
    straightBet("cowboys", { rotation: 403, awayTeam: "COWBOYS", homeTeam: "GIANTS", eventStart: new Date(NOW + HOUR).toISOString() }),
  ];
  const { games, coveredBetIds } = scenarios.buildScenarios({ records, games: gamesOf(state), now: NOW });
  assert.deepEqual(games.map((game) => game.eventId), [101]);
  assert.deepEqual(coveredBetIds, ["novig:bills", "novig:first-half", "novig:no-stake"]);
  assert.deepEqual(games[0].groups.map((group) => group.title), ["Result", "1st half result"]);
  assert.equal(games[0].bets.find((bet) => bet.id === "novig:no-stake").reason, "no stake on the record");
});

test("labels: MLB moneylines have no tie, totals count runs, an open band reads as a win", () => {
  const mlb = { league: "mlb", awayTeam: "Mariners", homeTeam: "Yankees" };
  assert.equal(scenarios.marginLabel(null, 0.5, mlb), "Yankees win");
  assert.equal(scenarios.marginLabel(0.5, null, mlb), "Mariners win");
  assert.equal(scenarios.marginLabel(-1.5, 1.5, mlb), "Yankees by 1 or Mariners by 1");
  assert.equal(scenarios.totalLabel(null, 8.5, mlb), "8 or fewer runs");
  assert.equal(scenarios.totalLabel(8.5, 10.5, mlb), "9-10 runs");
  const nfl = { league: "nfl", awayTeam: "Bills", homeTeam: "Chiefs" };
  assert.equal(scenarios.marginLabel(-0.5, 0.5, nfl), "Tie");
  assert.equal(scenarios.marginLabel(2.5, 6.5, nfl), "Bills by 3-6");
  assert.equal(scenarios.marginLabel(null, null, nfl), "Any result");
});

test("an MLB moneyline card: two results off the moneyline rungs", () => {
  const state = addGame(emptyBoard(), {
    eventId: 201, leagueId: MLB, away: "Seattle Mariners", home: "New York Yankees", awayRotation: 951, homeRotation: 952,
    spread: 1.5, moneyline: [[AWAY_OR_OVER, 135], [HOME_OR_UNDER, -135]],
  });
  const records = [straightBet("sea", { league: "mlb", betType: "moneyline", points: null, rotation: 951, awayTeam: "MARINERS", homeTeam: "YANKEES", stake: 100, toWin: 135 })];
  const [card] = scenarios.buildScenarios({ records, games: gamesOf(state), now: NOW }).games;
  assert.deepEqual(bandRows(card.groups[0]), [["New York Yankees win", -100], ["Seattle Mariners win", 135]]);
  nearly(card.groups[0].bands[1].prob, probOf(135));
});
