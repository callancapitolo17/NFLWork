// Run: node --test unabated_ticket/tests/ladder.test.js
const test = require("node:test");
const assert = require("node:assert/strict");
const ladder = require("../extension/ladder.js");

const NFL = 1;
const PREMIER_LEAGUE = 28;
const FULL_GAME = 1;
const FIRST_HALF = 2;
const MONEYLINE = 1;
const SPREAD = 2;
const TOTAL = 3;
const AWAY_OR_OVER = 0;
const HOME_OR_UNDER = 1;

const nearly = (actual, expected, tol = 1e-9) =>
  assert.ok(Math.abs(actual - expected) <= tol, `expected ${expected}, got ${actual}`);

function feedLine(fields) {
  return { eventId: 501, leagueId: NFL, periodTypeId: FULL_GAME, bookId: 4, fromSnapshot: true, ...fields };
}

function probAt(built, cut) {
  return ladder.probAbove(built, cut);
}

test("groupLinesByEvent files every line under its event and skips lines with none", () => {
  const lines = [feedLine({ eventId: 1 }), feedLine({ eventId: 2 }), feedLine({ eventId: 1 }), feedLine({ eventId: null }), null];
  const byEvent = ladder.groupLinesByEvent(lines);
  assert.deepEqual(Array.from(byEvent.keys()), [1, 2]);
  assert.equal(byEvent.get(1).length, 2);
});

test("totals: the cut is the number, an Over fair is P(above), an Under fair its complement", () => {
  const built = ladder.buildLadder([
    feedLine({ betTypeId: TOTAL, sideIndex: AWAY_OR_OVER, points: 47.5, bacr: -150 }),
    feedLine({ betTypeId: TOTAL, sideIndex: HOME_OR_UNDER, points: 52.5, bacr: 300 }),
  ], { periodTypeId: FULL_GAME, axis: ladder.AXIS_TOTAL });
  assert.deepEqual(built.rungs.map(([cut]) => cut), [47.5, 52.5]);
  nearly(probAt(built, 47.5).prob, 0.6);
  nearly(probAt(built, 52.5).prob, 0.75);
});

test("the rung's fair is the median across books, both sides counted", () => {
  const built = ladder.buildLadder([
    feedLine({ bookId: 1, betTypeId: TOTAL, sideIndex: AWAY_OR_OVER, points: 47.5, bacr: -150 }),
    feedLine({ bookId: 2, betTypeId: TOTAL, sideIndex: AWAY_OR_OVER, points: 47.5, bacr: 100 }),
    feedLine({ bookId: 3, betTypeId: TOTAL, sideIndex: HOME_OR_UNDER, points: 47.5, bacr: 300 }),
    feedLine({ bookId: 4, betTypeId: TOTAL, sideIndex: AWAY_OR_OVER, points: 47.5, bacr: null }),
  ], { periodTypeId: FULL_GAME, axis: ladder.AXIS_TOTAL });
  nearly(probAt(built, 47.5).prob, 0.6);
});

test("margin: an away spread at a wins above -a, a home spread at h wins below h, the moneyline is the +/-0.5 cut", () => {
  const built = ladder.buildLadder([
    feedLine({ betTypeId: SPREAD, sideIndex: AWAY_OR_OVER, points: 3.5, bacr: -150 }),
    feedLine({ betTypeId: SPREAD, sideIndex: HOME_OR_UNDER, points: -3.5, bacr: 150 }),
    feedLine({ betTypeId: SPREAD, sideIndex: HOME_OR_UNDER, points: 7.5, bacr: -400 }),
    feedLine({ betTypeId: MONEYLINE, sideIndex: AWAY_OR_OVER, points: null, bacr: 150 }),
    feedLine({ betTypeId: MONEYLINE, sideIndex: HOME_OR_UNDER, points: null, bacr: -150 }),
    feedLine({ betTypeId: TOTAL, sideIndex: AWAY_OR_OVER, points: 47.5, bacr: -150 }),
  ], { periodTypeId: FULL_GAME, axis: ladder.AXIS_MARGIN });
  assert.deepEqual(built.rungs.map(([cut]) => cut), [-3.5, -0.5, 0.5, 7.5]);
  nearly(probAt(built, -3.5).prob, 0.6);
  nearly(probAt(built, 7.5).prob, 0.2);
  nearly(probAt(built, 0.5).prob, 0.4);
});

test("the two moneyline cuts share one fair and are not flat; any other repeat is", () => {
  const built = ladder.buildLadder([
    feedLine({ betTypeId: MONEYLINE, sideIndex: AWAY_OR_OVER, points: null, bacr: 150 }),
    feedLine({ betTypeId: MONEYLINE, sideIndex: HOME_OR_UNDER, points: null, bacr: -150 }),
    feedLine({ betTypeId: SPREAD, sideIndex: AWAY_OR_OVER, points: -1.5, bacr: 150 }),
  ], { periodTypeId: FULL_GAME, axis: ladder.AXIS_MARGIN });
  nearly(probAt(built, -0.5).prob, 0.4);
  assert.deepEqual(probAt(built, 0.5), { reason: ladder.REASON_FLAT });
  assert.deepEqual(probAt(built, 1.5), { reason: ladder.REASON_FLAT });
});

test("a flat-lined tail is no fair: every rung that repeats a neighbour's", () => {
  const over = (points, bacr) => feedLine({ periodTypeId: FIRST_HALF, betTypeId: TOTAL, sideIndex: AWAY_OR_OVER, points, bacr });
  const built = ladder.buildLadder([over(33.5, 250), over(34.5, 300), over(35.5, 405), over(36.5, 405), over(37.5, 405)],
    { periodTypeId: FIRST_HALF, axis: ladder.AXIS_TOTAL });
  nearly(probAt(built, 34.5).prob, 0.25);
  for (const cut of [35.5, 36.5, 37.5]) assert.deepEqual(probAt(built, cut), { reason: ladder.REASON_FLAT });
});

test("no rung: a number the ladder lacks, a whole number, another period", () => {
  const built = ladder.buildLadder([
    feedLine({ betTypeId: TOTAL, sideIndex: AWAY_OR_OVER, points: 47.5, bacr: -150 }),
    feedLine({ betTypeId: TOTAL, sideIndex: AWAY_OR_OVER, points: 48, bacr: -120 }),
    feedLine({ periodTypeId: FIRST_HALF, betTypeId: TOTAL, sideIndex: AWAY_OR_OVER, points: 23.5, bacr: -110 }),
  ], { periodTypeId: FULL_GAME, axis: ladder.AXIS_TOTAL });
  assert.deepEqual(built.rungs.map(([cut]) => cut), [47.5]);
  for (const cut of [36.5, 48, 23.5]) assert.deepEqual(probAt(built, cut), { reason: ladder.REASON_NO_RUNG });
  assert.deepEqual(ladder.probAbove(null, 47.5), { reason: ladder.REASON_NO_RUNG });
});

test("a line the changes stream added is not read: it may be a team total filed under the game total", () => {
  const built = ladder.buildLadder([
    feedLine({ betTypeId: TOTAL, sideIndex: AWAY_OR_OVER, points: 47.5, bacr: -150 }),
    feedLine({ betTypeId: TOTAL, sideIndex: AWAY_OR_OVER, points: 23.5, bacr: -110, fromSnapshot: false }),
    feedLine({ betTypeId: TOTAL, sideIndex: AWAY_OR_OVER, points: 47.5, bacr: 400, fromSnapshot: undefined }),
  ], { periodTypeId: FULL_GAME, axis: ladder.AXIS_TOTAL });
  assert.deepEqual(built.rungs, [[47.5, 0.6]]);
});

test("a soccer moneyline is three-way and feeds no margin cut", () => {
  const built = ladder.buildLadder([
    feedLine({ leagueId: PREMIER_LEAGUE, betTypeId: MONEYLINE, sideIndex: AWAY_OR_OVER, points: null, bacr: 250 }),
    feedLine({ leagueId: PREMIER_LEAGUE, betTypeId: SPREAD, sideIndex: AWAY_OR_OVER, points: 0.5, bacr: -120 }),
  ], { periodTypeId: FULL_GAME, axis: ladder.AXIS_MARGIN });
  assert.deepEqual(built.rungs.map(([cut]) => cut), [-0.5]);
});

test("a fair that is not an American price is skipped, not thrown on", () => {
  const built = ladder.buildLadder([
    feedLine({ betTypeId: TOTAL, sideIndex: AWAY_OR_OVER, points: 47.5, bacr: 50 }),
    feedLine({ betTypeId: TOTAL, sideIndex: AWAY_OR_OVER, points: 47.5, bacr: -150 }),
  ], { periodTypeId: FULL_GAME, axis: ladder.AXIS_TOTAL });
  nearly(probAt(built, 47.5).prob, 0.6);
});

test("periodTypeIdOf reads feed.PERIODS backwards; a period Unabated has no id for is null", () => {
  assert.equal(ladder.periodTypeIdOf("FG"), 1);
  assert.equal(ladder.periodTypeIdOf("1H"), 2);
  assert.equal(ladder.periodTypeIdOf("F5"), null);
  assert.equal(ladder.periodTypeIdOf(null), null);
});
