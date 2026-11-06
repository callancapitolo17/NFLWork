// Run: node --test unabated_ticket/tests
// server/gamescore.js: the Live card's score line off a game's feed state.
const test = require("node:test");
const assert = require("node:assert/strict");
const gamescore = require("../server/gamescore.js");

const NOW = Date.parse("2026-11-06T02:10:00Z");
const READ_AT = NOW - 40 * 1000;

function gameState(fields) {
  return {
    statusId: 2, awayScore: 10, homeScore: 3, gameClock: "6:43", periodId: 5,
    snapshotAtMs: READ_AT, awayAbbr: "JAX", homeAbbr: "BAL", ...fields,
  };
}

function line(fields) {
  return {
    status: "live", awayAbbr: "JAX", homeAbbr: "BAL", awayScore: 10, homeScore: 3, clock: "6:43", period: "Q2",
    readAt: new Date(READ_AT).toISOString(), staleAt: new Date(READ_AT + gamescore.STALE_AFTER_MS).toISOString(), ...fields,
  };
}

test("a football game in progress shows its score, quarter and clock", () => {
  assert.deepEqual(gamescore.scoreLineOf(gameState(), "nfl", NOW), line());
  assert.deepEqual(gamescore.scoreLineOf(gameState({ periodId: 3 }), "cbb", NOW), line({ period: "2nd half" }));
});

test("a hockey game names its period 8 / 9 / 10 as P1 / P2 / P3", () => {
  assert.deepEqual(gamescore.scoreLineOf(gameState({ periodId: 8 }), "nhl", NOW), line({ period: "P1" }));
  assert.deepEqual(gamescore.scoreLineOf(gameState({ periodId: 9 }), "nhl", NOW), line({ period: "P2" }));
  assert.deepEqual(gamescore.scoreLineOf(gameState({ periodId: 10 }), "nhl", NOW), line({ period: "P3" }));
  assert.deepEqual(gamescore.scoreLineOf(gameState({ periodId: 10 }), "nba", NOW), line({ period: null }));
});

test("a period number not checked for the league shows no name, the clock still does", () => {
  assert.deepEqual(gamescore.scoreLineOf(gameState({ periodId: 11 }), "nhl", NOW), line({ period: null }));
  assert.deepEqual(gamescore.scoreLineOf(gameState({ periodId: 5 }), "nhl", NOW), line({ period: null }));
  assert.deepEqual(gamescore.scoreLineOf(gameState({ periodId: 9 }), "nfl", NOW), line({ period: null }));
  assert.deepEqual(gamescore.scoreLineOf(gameState({ gameClock: null, periodId: null }), "mlb", NOW), line({ clock: null, period: null }));
});

test("a final game keeps its score and shows no clock (the feed's end marker is not a time)", () => {
  assert.deepEqual(gamescore.scoreLineOf(gameState({ statusId: 3, gameClock: "E" }), "nfl", NOW), line({ status: "final", clock: null, period: null }));
});

test("an upcoming or unseen status is other, with no clock", () => {
  assert.deepEqual(gamescore.scoreLineOf(gameState({ statusId: 1, awayScore: 0, homeScore: 0 }), "nfl", NOW),
    line({ status: "other", awayScore: 0, homeScore: 0, clock: null, period: null }));
  assert.equal(gamescore.scoreLineOf(gameState({ statusId: 7 }), "nfl", NOW).status, "other");
});

test("a score older than five minutes, or with no cut time, reads stale", () => {
  assert.deepEqual(gamescore.scoreLineOf(gameState({ snapshotAtMs: NOW - gamescore.STALE_AFTER_MS - 1 }), "nfl", NOW), { status: "stale" });
  assert.deepEqual(gamescore.scoreLineOf(gameState({ snapshotAtMs: null }), "nfl", NOW), { status: "stale" });
});

test("no feed state, no score line", () => {
  assert.equal(gamescore.scoreLineOf(null, "nfl", NOW), null);
});
