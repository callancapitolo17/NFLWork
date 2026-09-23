// Run: node --test unabated_ticket/tests
//
// The Bets tab's Attach picker (attach.js): step 1's game list, step 2's
// name plan, and the POST /pins.json body. Board rows are describeLine-shaped
// with the real Unabated team ids of tests/fixtures/teams_index.json.
const test = require("node:test");
const assert = require("node:assert/strict");
const fs = require("node:fs");
const path = require("node:path");
const bets = require("../extension/bets.js");
const teams = require("../extension/teams.js");
const attach = require("../extension/attach.js");
teams.loadIndex(JSON.parse(fs.readFileSync(path.join(__dirname, "fixtures", "teams_index.json"), "utf8")).leagues);

function boardRow(fields) {
  return Object.assign({ betType: "Total", period: "FG", sideIndex: 0, points: 50.5 }, fields, {
    eventStartMs: Date.parse(fields.eventStart), leagueLabel: fields.league.toUpperCase(),
  });
}

// Sat Sep 26 8:00 PM ET is 00:00 UTC on the 27th: the Eastern date decides the window.
const ABILENE_AT_TARLETON = boardRow({ league: "cfb", eventId: 7001, awayTeam: "Abilene Christian", homeTeam: "Tarleton State",
  awayTeamId: 1123, homeTeamId: 1196, awayRotation: 371, homeRotation: 372, eventStart: "2026-09-27T00:00:00Z" });
const LAMAR_AT_NORTHWESTERN = boardRow({ league: "cfb", eventId: 7002, awayTeam: "Lamar", homeTeam: "Northwestern State",
  awayTeamId: 1160, homeTeamId: 1177, awayRotation: 367, homeRotation: 368, eventStart: "2026-09-26T23:00:00Z" });
const CHATTANOOGA_NEXT_WEEK = boardRow({ league: "cfb", eventId: 7003, awayTeam: "Chattanooga", homeTeam: "Eastern Kentucky",
  awayTeamId: 791, homeTeamId: 1062, awayRotation: 401, homeRotation: 402, eventStart: "2026-10-03T23:00:00Z" });
const SOX_AT_TIGERS_G1 = boardRow({ league: "mlb", eventId: 8001, awayTeam: "Chicago White Sox", homeTeam: "Detroit Tigers",
  awayTeamId: 38, homeTeamId: 42, awayRotation: 951, homeRotation: 952, eventStart: "2026-09-26T17:10:00Z" });
const SOX_AT_TIGERS_G2 = boardRow({ league: "mlb", eventId: 8002, awayTeam: "Chicago White Sox", homeTeam: "Detroit Tigers",
  awayTeamId: 38, homeTeamId: 42, awayRotation: 953, homeRotation: 954, eventStart: "2026-09-26T22:40:00Z" });
const BOARD = [ABILENE_AT_TARLETON, LAMAR_AT_NORTHWESTERN, CHATTANOOGA_NEXT_WEEK, SOX_AT_TIGERS_G1, SOX_AT_TIGERS_G2];

function keyed(fields) {
  const [record] = bets.resolveTeamKeys([Object.assign({ status: "open", period: "FG", stake: 100, toWin: 90, price: -110 }, fields)], []);
  return record;
}

// BFA writes "Abilene Chr", which no rule resolves; "Tarleton St" resolves ("st" reads as "state").
const BFA_TOTAL = keyed({ id: "bfa:1", venue: "bfa", league: "cfb", betType: "total", period: "1H", side: "under", points: 30.5,
  awayTeam: "Abilene Chr", homeTeam: "Tarleton St", eventStart: "2026-09-27T00:00:00Z" });

test("step 1: the bet's league around its Eastern date, a game with a team that already resolves first", () => {
  assert.equal(BFA_TOTAL.awayKey, null);
  assert.equal(BFA_TOTAL.homeKey, "cfb:1196");
  const { scope, events, more } = attach.attachCandidates(BFA_TOTAL, BOARD, {});
  assert.equal(scope, "CFB · Fri Sep 25 to Sun Sep 27");
  assert.deepEqual(events.map(({ event, why }) => [event.eventId, why]), [[7001, "Tarleton St matches"], [7002, null]]);
  assert.equal(more, 0);
  assert.equal(events[0].event.homeKey, "cfb:1196");
});

test("step 1: a query of two letters searches every date of the league; one letter does not", () => {
  const searched = attach.attachCandidates(BFA_TOTAL, BOARD, { query: " chatt " });
  assert.equal(searched.scope, "CFB · every date");
  assert.deepEqual(searched.events.map(({ event }) => event.eventId), [7003]);
  assert.deepEqual(attach.attachCandidates(BFA_TOTAL, BOARD, { query: "c" }).events.map(({ event }) => event.eventId), [7001, 7002]);
});

test("step 1: a rotation fits too; a bet with no league or date searches every league and date", () => {
  const byRotation = keyed({ id: "bol:1", venue: "betonline", league: "cfb", betType: "spread", side: "away", points: 7.5,
    awayTeam: "ACU", homeTeam: null, rotation: 367, eventDate: "2026-09-26" });
  assert.deepEqual(attach.attachCandidates(byRotation, BOARD, {}).events.map(({ event, why }) => [event.eventId, why]),
    [[7002, "rot 367 matches"], [7001, null]]);
  const bare = keyed({ id: "x:1", venue: "betonline", league: null, betType: "total", side: "over", points: 8.5, awayTeam: "Somebody" });
  const everything = attach.attachCandidates(bare, BOARD, {});
  assert.equal(everything.scope, "Every league · every date");
  assert.equal(everything.events.length, BOARD.length);
});

test("step 1: more than a screen of games is capped and counted", () => {
  const many = Array.from({ length: attach.MAX_CANDIDATES + 3 }, (_, index) => boardRow({ league: "cfb", eventId: 9000 + index,
    awayTeam: `Away ${index}`, homeTeam: `Home ${index}`, awayTeamId: 1, homeTeamId: 2, eventStart: "2026-09-26T20:00:00Z" }));
  const { events, more } = attach.attachCandidates(BFA_TOTAL, many, {});
  assert.equal(events.length, attach.MAX_CANDIDATES);
  assert.equal(more, 3);
});

test("step 2: each venue name lines up with the same side; unknown names are learned, known ones are not written", () => {
  const plan = attach.attachPlan(BFA_TOTAL, attach.boardGames(BOARD, "cfb")[0], { swapped: false, lines: BOARD });
  assert.deepEqual(plan.names.map((name) => [name.venueTeamKey, name.unabatedTeamName, name.eventSide, name.status, name.was]), [
    ["Abilene Chr", "Abilene Christian", "away", "learn", null],
    ["Tarleton St", "Tarleton State", "home", "known", null],
  ]);
  assert.deepEqual(plan.crosswalk, [{
    venue: "bfa", league: "cfb", venueTeamKey: "Abilene Chr", venueTeamName: "Abilene Chr",
    unabatedTeamId: "1123", unabatedTeamName: "Abilene Christian", learnedFrom: "attached by you: bfa:1 on board event 7001",
  }]);
  assert.equal(plan.betOnGame, bets.describeBet(BFA_TOTAL));
});

test("step 2: Swap flips the sides, and a name that resolved to the other team is a fix that says what it was", () => {
  const plan = attach.attachPlan(BFA_TOTAL, attach.boardGames(BOARD, "cfb")[0], { swapped: true, lines: BOARD });
  assert.deepEqual(plan.names.map((name) => [name.venueTeamKey, name.unabatedTeamName, name.status, name.was]), [
    ["Abilene Chr", "Tarleton State", "learn", null],
    ["Tarleton St", "Abilene Christian", "fix", "Tarleton State"],
  ]);
  assert.deepEqual(plan.crosswalk.map((row) => [row.venueTeamKey, row.unabatedTeamId]), [["Abilene Chr", "1196"], ["Tarleton St", "1123"]]);
});

test("step 2: a one-team spread shows one name and restates the bet on the side it lands on", () => {
  const spread = keyed({ id: "bol:2", venue: "betonline", league: "cfb", betType: "spread", side: "away", points: 7.5,
    awayTeam: "Abilene Chr", homeTeam: null, eventDate: "2026-09-26" });
  const game = attach.boardGames(BOARD, "cfb")[0];
  const plan = attach.attachPlan(spread, game, { swapped: false, lines: BOARD });
  assert.deepEqual(plan.names.map((name) => [name.venueTeamKey, name.unabatedTeamName]), [["Abilene Chr", "Abilene Christian"]]);
  assert.equal(plan.betOnGame, "Abilene Christian +7.5");
  assert.equal(attach.attachPlan(spread, game, { swapped: true, lines: BOARD }).betOnGame, "Tarleton State +7.5");
  const firstHalfMoneyline = keyed({ ...spread, id: "bol:3", betType: "moneyline", period: "1H", points: null });
  assert.equal(attach.attachPlan(firstHalfMoneyline, game, { lines: BOARD }).betOnGame, "1H Abilene Christian to win");
  // Named on neither side: the rotation says which, or the step says it cannot tell.
  const byRotation = keyed({ ...spread, id: "bol:4", awayTeam: null, rotation: 372 });
  assert.equal(attach.attachPlan(byRotation, game, { lines: BOARD }).betOnGame, "Tarleton State +7.5");
  const blind = keyed({ ...spread, id: "bol:5", awayTeam: null });
  assert.equal(attach.attachPlan(blind, game, { lines: BOARD }).betOnGame, "side not known: it will not size a line on this game");
});

test("step 2: a doubleheader bet whose names both resolve teaches nothing; the attach only pins the game", () => {
  const doubleheader = keyed({ id: "wz:1", venue: "wagerzon", league: "mlb", betType: "spread", side: "away", points: -1.5,
    awayTeam: "Chicago White Sox", homeTeam: "Detroit Tigers", eventDate: "2026-09-26" });
  const [gameOne, gameTwo] = attach.attachCandidates(doubleheader, BOARD, {}).events.map(({ event }) => event);
  assert.deepEqual([gameOne.eventId, gameTwo.eventId], [8001, 8002]);
  const plan = attach.attachPlan(doubleheader, gameTwo, { lines: BOARD });
  assert.deepEqual(plan.names.map((name) => name.status), ["known", "known"]);
  assert.deepEqual(plan.crosswalk, []);
  assert.equal(plan.betOnGame, "Chicago White Sox -1.5");
});

test("pinRequest: the pin names the event, its start and both teams as strings; the plan's rows ride along", () => {
  const game = attach.boardGames(BOARD, "cfb")[0];
  const plan = attach.attachPlan(BFA_TOTAL, game, { lines: BOARD });
  assert.deepEqual(attach.pinRequest(BFA_TOTAL, game, plan), {
    pin: { betId: "bfa:1", venue: "bfa", league: "cfb", eventId: "7001", eventStart: "2026-09-27T00:00:00.000Z",
      awayTeamId: "1123", homeTeamId: "1196", awayTeamName: "Abilene Christian", homeTeamName: "Tarleton State" },
    crosswalk: plan.crosswalk,
  });
});
