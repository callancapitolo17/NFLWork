// Run: node --test unabated_ticket/tests/novigquote.test.js
const test = require("node:test");
const assert = require("node:assert/strict");
const novigQuote = require("../server/novigquote.js");

const MAIN_0 = { key: "m0", bookId: 89, eventId: 7, periodTypeId: 1, betTypeId: 2, sideIndex: 0, points: 3.5 };
const MAIN_1 = { key: "m1", bookId: 89, eventId: 7, periodTypeId: 1, betTypeId: 2, sideIndex: 1, points: -3.5 };
const ALT_1 = { key: "a1", bookId: 89, eventId: 7, periodTypeId: 1, betTypeId: 2, sideIndex: 1, points: -6.5 };
const ID = (n) => `00000000-0000-0000-0000-00000000000${n}`;

function boardWith(outcomes) {
  return { lines: { m0: MAIN_0, m1: MAIN_1, a1: ALT_1 }, events: { 7: { venueIds: { novigOutcomes: outcomes } } } };
}

test("a spread line finds its outcome and the other side's at the flipped number", () => {
  const board = boardWith({
    [ID(1)]: { lineKey: "m1", mainKey: "m1", points: -3.5, sideIndex: 1 },
    [ID(2)]: { lineKey: "m0", mainKey: "m0", points: 3.5, sideIndex: 0 },
    [ID(3)]: { lineKey: "a1", mainKey: "m1", points: -6.5, sideIndex: 1 },
    [ID(4)]: { lineKey: null, mainKey: "m0", points: 6.5, sideIndex: 0 },
  });
  assert.deepEqual(novigQuote.novigOutcomesOf(MAIN_1, board), { outcomeId: ID(1), otherOutcomeId: ID(2) });
  // The other side's alt rung counts even when Unabated lists no price for it.
  assert.deepEqual(novigQuote.novigOutcomesOf(ALT_1, board), { outcomeId: ID(3), otherOutcomeId: ID(4) });
});

test("a total keeps the number on the other side", () => {
  const over = { ...MAIN_0, key: "t0", betTypeId: 3, points: 47.5 };
  const under = { ...MAIN_1, key: "t1", betTypeId: 3, points: 47.5 };
  const board = { lines: { t0: over, t1: under }, events: { 7: { venueIds: { novigOutcomes: {
    [ID(5)]: { lineKey: "t0", mainKey: "t0", points: 47.5, sideIndex: 0 },
    [ID(6)]: { lineKey: "t1", mainKey: "t1", points: 47.5, sideIndex: 1 },
  } } } } };
  assert.deepEqual(novigQuote.novigOutcomesOf(under, board), { outcomeId: ID(6), otherOutcomeId: ID(5) });
});

test("a moneyline, a missing id, a missing other side and another market are refused", () => {
  const board = boardWith({ [ID(1)]: { lineKey: "m1", mainKey: "m1", points: -3.5, sideIndex: 1 } });
  assert.equal(novigQuote.novigOutcomesOf({ ...MAIN_1, betTypeId: 1 }, board).reason, "no_outcome");
  assert.equal(novigQuote.novigOutcomesOf({ ...MAIN_1, bookId: 105 }, board).reason, "not_novig");
  assert.match(novigQuote.novigOutcomesOf(MAIN_1, board).detail, /0 Novig outcome ids for the other side/);
  assert.match(novigQuote.novigOutcomesOf(MAIN_0, board).detail, /0 Novig outcome ids for this line/);
  // An id on a first-half line never prices the full game.
  const half = { ...MAIN_0, key: "h0", periodTypeId: 2 };
  const mixed = { lines: { m1: MAIN_1, h0: half }, events: { 7: { venueIds: { novigOutcomes: {
    [ID(1)]: { lineKey: "m1", mainKey: "m1", points: -3.5, sideIndex: 1 },
    [ID(2)]: { lineKey: "h0", mainKey: "h0", points: 3.5, sideIndex: 0 },
  } } } } };
  assert.equal(novigQuote.novigOutcomesOf(MAIN_1, mixed).reason, "no_outcome");
});
