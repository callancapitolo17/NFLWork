// Run: node --test unabated_ticket/tests
// server/closefair.js: which open bets get a closing fair, what is resent, and
// which leagues the runner must load for them. The board lookup itself runs
// end to end in runner.test.js against the NFL slice.
const test = require("node:test");
const assert = require("node:assert/strict");
const closefair = require("../server/closefair.js");

function bet(overrides) {
  return Object.assign({
    id: "b1", venue: "kalshi", league: "mlb", status: "open", betType: "total", period: "FG", side: "over", points: 8.5,
    placedAt: "2026-10-04T15:00:00Z", isParlayLeg: false, unmatchable: null, approx: [],
  }, overrides);
}

test("leagueIdsOfOpenBets: every Unabated league on the open straight bets' paths, never a leg's or a settled bet's", () => {
  assert.deepEqual(closefair.leagueIdsOfOpenBets([bet()]), [5, 12]);
  assert.deepEqual(closefair.leagueIdsOfOpenBets([
    bet({ league: "nhl", isParlayLeg: true }), bet({ league: "nba", status: "won" }),
    bet({ league: "cbb", betType: "other" }), bet({ league: "wnba", unmatchable: "no game" }),
  ]), []);
});

test("unsentRows resends a bet only when its line, fair or reading time changed", () => {
  const row = { betId: "b1", lineKey: "m1:ms105:si0:tid6", fairAmerican: -120, fairObservedAt: "2026-10-04T16:00:00.000Z" };
  const sent = new Map();
  assert.deepEqual(closefair.unsentRows([row], sent), [row]);
  closefair.markSent([row], sent);
  assert.deepEqual(closefair.unsentRows([row], sent), []);
  const newer = { ...row, fairObservedAt: "2026-10-04T16:02:00.000Z" };
  assert.deepEqual(closefair.unsentRows([newer, { ...row, betId: "b2" }], sent), [newer, { ...row, betId: "b2" }]);
});

test("forgetClosed keeps only what was sent for bets still open", () => {
  const sent = new Map([["b1", "x"], ["b2", "y"]]);
  closefair.forgetClosed(sent, [bet({ id: "b1" }), bet({ id: "b2", status: "won" })]);
  assert.deepEqual([...sent.keys()], ["b1"]);
});

test("closingFairRows: nothing without a board, or without an open straight bet", () => {
  assert.deepEqual(closefair.closingFairRows({ records: [bet()], state: null, boardLines: [], leagueLoadedAt: {}, now: 0 }), []);
  assert.deepEqual(closefair.closingFairRows({ records: [bet({ status: "won" })], state: { lines: {} }, boardLines: [], leagueLoadedAt: {}, now: 0 }), []);
});
