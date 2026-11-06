// Run: node --test unabated_ticket/tests
// server/gamescore.js: the Live card's score line off a game's feed state.
const test = require("node:test");
const assert = require("node:assert/strict");
const gamescore = require("../server/gamescore.js");

const NOW = Date.parse("2026-11-06T02:10:00Z");

function gameState(fields) {
  return {
    statusId: 2, awayScore: 10, homeScore: 3, gameClock: "6:43", periodId: 5, overtimeNumber: null,
    snapshotAtMs: NOW - 40 * 1000, awayAbbr: "JAX", homeAbbr: "BAL", ...fields,
  };
}

test("a game in progress shows its score and clock; an unverified period number shows no name", () => {
  assert.deepEqual(gamescore.scoreLineOf(gameState(), "nfl", NOW), {
    status: "live", awayAbbr: "JAX", homeAbbr: "BAL", awayScore: 10, homeScore: 3, clock: "6:43", period: null,
    readAt: new Date(NOW - 40 * 1000).toISOString(),
  });
});

test("a period number in the league's verified table is named", (t) => {
  gamescore.PERIOD_NAMES.testleague = { 5: "2nd quarter" };
  t.after(() => delete gamescore.PERIOD_NAMES.testleague);
  assert.equal(gamescore.scoreLineOf(gameState(), "testleague", NOW).period, "2nd quarter");
  assert.equal(gamescore.scoreLineOf(gameState({ periodId: 6 }), "testleague", NOW).period, null);
});

test("a final game shows no clock (the feed's end marker is not a time)", () => {
  const line = gamescore.scoreLineOf(gameState({ statusId: 3, gameClock: "E" }), "nfl", NOW);
  assert.equal(line.status, "final");
  assert.equal(line.clock, null);
  assert.equal(line.period, null);
});

test("a score older than five minutes reads stale, not as the current score", () => {
  const old = NOW - gamescore.STALE_AFTER_MS - 1;
  assert.deepEqual(gamescore.scoreLineOf(gameState({ snapshotAtMs: old }), "nfl", NOW), { status: "stale", readAt: new Date(old).toISOString() });
  assert.equal(gamescore.scoreLineOf(gameState({ snapshotAtMs: null }), "nfl", NOW).status, "stale");
});

test("no feed state, no score line; an unknown status is named unknown", () => {
  assert.equal(gamescore.scoreLineOf(null, "nfl", NOW), null);
  assert.equal(gamescore.scoreLineOf(gameState({ statusId: 7 }), "nfl", NOW).status, "unknown");
});
