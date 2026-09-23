// Run: node --test unabated_ticket/tests
// The edge-move tag (#132): a per-line history and the four-case rule the
// fair decides. Probabilities: -150 = 60.0%, -160 = 61.5%, -151 = 60.2%,
// +199 = 33.4%, +215 = 31.7%, +205 = 32.8%, -140 = 58.3%.
const test = require("node:test");
const assert = require("node:assert/strict");
const edgemove = require("../extension/edgemove.js");
const { observe, forget, edgeMove, EDGE_MOVE_WINDOW_MS, MOVE_LABELS } = edgemove;

const T0 = Date.UTC(2026, 8, 22, 17, 0, 0);
const MIN = 60 * 1000;
const NOW = T0 + 5 * MIN;

function line(overrides = {}) {
  return { key: "m1:ms4:si0:tid6", points: 47.5, price: 199, sourceFormat: 1, sourcePrice: null, bacr: -150, ...overrides };
}

// Observe a sequence of {at, source, ...line fields} on one key.
function historyOf(observations) {
  const history = {};
  for (const { at, source = "snapshot", ...fields } of observations) observe(history, line(fields), { at, source });
  return history;
}

const pts = (delta) => Math.round(delta * 1000) / 10; // probability points, one decimal

test("fair moved toward the side, price flat: fair moved to you", () => {
  const history = historyOf([{ at: T0, bacr: -150 }, { at: T0 + 2 * MIN, bacr: -160 }]);
  const move = edgeMove(history["m1:ms4:si0:tid6"], NOW);
  assert.equal(move.kind, "fair_to_you");
  assert.equal(pts(move.fairDelta), 1.5);
  assert.equal(move.priceDelta, 0);
  assert.equal(move.sinceMs, 3 * MIN);
  assert.equal(move.source, "snapshot");
  assert.equal(move.from.bacr, -150);
  assert.equal(move.to.bacr, -160);
  assert.equal(MOVE_LABELS[move.kind], "fair moved to you");
});

test("price got better, fair flat: book moved away", () => {
  const history = historyOf([{ at: T0, price: 199 }, { at: T0 + 2 * MIN, price: 215, source: "stream" }]);
  const move = edgeMove(history["m1:ms4:si0:tid6"], NOW);
  assert.equal(move.kind, "book_away");
  assert.equal(move.fairDelta, 0);
  assert.equal(pts(move.priceDelta), 1.7);
  assert.equal(move.source, "stream");
  assert.equal(MOVE_LABELS[move.kind], "book moved away");
});

test("the fair decides: fair up stays green whether the price got better or worse", () => {
  const better = historyOf([{ at: T0 }, { at: T0 + MIN, bacr: -160, price: 215 }]);
  assert.equal(edgeMove(better["m1:ms4:si0:tid6"], NOW).kind, "fair_to_you");
  const worse = historyOf([{ at: T0 }, { at: T0 + MIN, bacr: -160, price: 185 }]);
  const move = edgeMove(worse["m1:ms4:si0:tid6"], NOW);
  assert.equal(move.kind, "fair_to_you");
  assert.ok(move.priceDelta < 0);
});

test("fair moved against the side while the price got better: red, even with the edge bigger", () => {
  // Fair 60.0% -> 58.3% (-1.7 pts); price +199 -> +250 (33.4% -> 28.6%, +4.9 pts):
  // EV per $1 grows from 0.60*2.99-1 = 0.79 to 0.583*3.5-1 = 1.04 and the row
  // would read as an improvement. It is the adverse-selection case.
  const history = historyOf([{ at: T0 }, { at: T0 + MIN, bacr: -140, price: 250, source: "stream" }]);
  const move = edgeMove(history["m1:ms4:si0:tid6"], NOW);
  assert.equal(move.kind, "fair_against");
  assert.equal(pts(move.fairDelta), -1.7);
  assert.equal(pts(move.priceDelta), 4.9);
  assert.equal(MOVE_LABELS[move.kind], "fair moved against you");
});

test("fair moved against the side with the price flat is red too", () => {
  const history = historyOf([{ at: T0 }, { at: T0 + MIN, bacr: -140 }]);
  assert.equal(edgeMove(history["m1:ms4:si0:tid6"], NOW).kind, "fair_against");
});

test("the amber tag turns red once the fair follows the book", () => {
  const history = historyOf([{ at: T0 }, { at: T0 + MIN, price: 215 }]);
  assert.equal(edgeMove(history["m1:ms4:si0:tid6"], T0 + 2 * MIN).kind, "book_away");
  observe(history, line({ price: 215, bacr: -140 }), { at: T0 + 3 * MIN, source: "snapshot" });
  const move = edgeMove(history["m1:ms4:si0:tid6"], T0 + 4 * MIN);
  assert.equal(move.kind, "fair_against");
  assert.equal(move.sinceMs, 3 * MIN); // the book's move was first, and is what "ago" counts from
  assert.equal(move.from.price, 199);
  assert.equal(move.to.bacr, -140);
});

test("a book that only shortened earns no tag", () => {
  const history = historyOf([{ at: T0, price: 199 }, { at: T0 + MIN, price: 185 }]);
  const move = edgeMove(history["m1:ms4:si0:tid6"], NOW);
  assert.equal(move.kind, "none");
  assert.ok(move.priceDelta < 0);
});

test("a whole-point fair change counts; a one-point one near even money does not", () => {
  const whole = historyOf([{ at: T0, bacr: -150 }, { at: T0 + MIN, bacr: -160 }]);
  assert.equal(edgeMove(whole["m1:ms4:si0:tid6"], NOW).kind, "fair_to_you");
  const one = historyOf([{ at: T0, bacr: -150 }, { at: T0 + MIN, bacr: -151 }]);
  const move = edgeMove(one["m1:ms4:si0:tid6"], NOW);
  assert.equal(move.kind, "none");
  assert.equal(pts(move.fairDelta), 0.2);
  // The sub-threshold change is still recorded, so a second point later adds up.
  assert.equal(one["m1:ms4:si0:tid6"].length, 2);
});

test("a Novig half-cent is a move (float boundary)", () => {
  // Under side of a total on Novig: 0.565 -> 0.56 is 0.00499999... in floats.
  const history = historyOf([
    { at: T0, sourceFormat: 4, sourcePrice: 0.565, price: -130 },
    { at: T0 + MIN, sourceFormat: 4, sourcePrice: 0.56, price: -127 },
  ]);
  const move = edgeMove(history["m1:ms4:si0:tid6"], NOW);
  assert.equal(move.kind, "book_away");
  assert.ok(Math.abs(move.priceDelta - 0.005) < 1e-6);
});

test("first sighting and nothing moved read as no tag", () => {
  assert.deepEqual(edgeMove(undefined, NOW), { kind: "none", fairDelta: null, priceDelta: null, sinceMs: null, source: null, from: null, to: null });
  const history = historyOf([{ at: T0 }, { at: T0 + MIN }, { at: T0 + 2 * MIN }]);
  assert.equal(history["m1:ms4:si0:tid6"].length, 1);
  assert.equal(edgeMove(history["m1:ms4:si0:tid6"], NOW).kind, "none");
});

test("an alt rung's history is snapshot-only and the tag says so", () => {
  const altKey = "m1:ms4:si0:tid6:alt48.5";
  const history = {};
  observe(history, line({ key: altKey, points: 48.5, price: 150, isAlt: true }), { at: T0, source: "snapshot" });
  observe(history, line({ key: altKey, points: 48.5, price: 165, isAlt: true }), { at: T0 + 2 * MIN, source: "snapshot" });
  const move = edgeMove(history[altKey], NOW);
  assert.equal(move.kind, "book_away");
  assert.equal(move.source, "snapshot");
  assert.equal(move.sinceMs, 3 * MIN);
});

test("a number move resets the history: the line reads as first seen", () => {
  const history = historyOf([{ at: T0, points: 47.5, price: 199 }, { at: T0 + MIN, points: 48.5, price: 215 }]);
  assert.equal(history["m1:ms4:si0:tid6"].length, 1);
  assert.equal(history["m1:ms4:si0:tid6"][0].points, 48.5);
  assert.equal(edgeMove(history["m1:ms4:si0:tid6"], NOW).kind, "none");
  // Moneylines have no number and never reset on it.
  const ml = historyOf([{ at: T0, points: null }, { at: T0 + MIN, points: null, bacr: -160 }]);
  assert.equal(edgeMove(ml["m1:ms4:si0:tid6"], NOW).kind, "fair_to_you");
});

test("a move older than the window is no longer a move", () => {
  const history = historyOf([{ at: T0 }, { at: T0 + MIN, bacr: -160 }]);
  assert.equal(edgeMove(history["m1:ms4:si0:tid6"], T0 + MIN + EDGE_MOVE_WINDOW_MS - 1).kind, "fair_to_you");
  assert.equal(edgeMove(history["m1:ms4:si0:tid6"], T0 + MIN + EDGE_MOVE_WINDOW_MS).kind, "none");
});

test("a line first seen inside the window compares to its first sighting", () => {
  const history = historyOf([{ at: NOW - 2 * MIN }, { at: NOW - MIN, bacr: -160 }]);
  const move = edgeMove(history["m1:ms4:si0:tid6"], NOW);
  assert.equal(move.kind, "fair_to_you");
  assert.equal(move.sinceMs, MIN);
});

test("a missing fair decides on the price alone", () => {
  const history = historyOf([{ at: T0, bacr: null }, { at: T0 + MIN, bacr: null, price: 215 }]);
  const move = edgeMove(history["m1:ms4:si0:tid6"], NOW);
  assert.equal(move.kind, "book_away");
  assert.equal(move.fairDelta, null);
  const flat = historyOf([{ at: T0, bacr: null }, { at: T0 + MIN, bacr: -150 }]);
  assert.equal(edgeMove(flat["m1:ms4:si0:tid6"], NOW).kind, "none");
});

test("pruning keeps one baseline older than the window and drops the rest", () => {
  const history = historyOf([
    { at: T0, bacr: -150 },
    { at: T0 + MIN, bacr: -155 },
    { at: T0 + 2 * MIN, bacr: -160 },
    { at: T0 + 2 * MIN + EDGE_MOVE_WINDOW_MS, bacr: -165 },
  ]);
  const entries = history["m1:ms4:si0:tid6"];
  assert.deepEqual(entries.map((entry) => entry.bacr), [-160, -165]);
  const move = edgeMove(entries, T0 + 3 * MIN + EDGE_MOVE_WINDOW_MS);
  assert.equal(move.kind, "fair_to_you");
  assert.equal(move.from.bacr, -160);
});

test("forget drops a key; observe and edgeMove fail loudly on bad input", () => {
  const history = historyOf([{ at: T0 }]);
  forget(history, "m1:ms4:si0:tid6");
  assert.deepEqual(history, {});
  assert.throws(() => observe(history, line({ price: null }), { at: T0, source: "snapshot" }), /expected a numeric price/);
  assert.throws(() => observe(history, line(), { at: T0, source: "page" }), /expected source snapshot\|stream/);
  assert.throws(() => observe(history, line(), { at: "now", source: "stream" }), /expected a numeric time/);
  assert.throws(() => observe(null, line(), { at: T0, source: "stream" }), /expected a history object/);
  assert.throws(() => edgeMove({}, NOW), /expected an array/);
  assert.throws(() => edgeMove([], NaN), /expected a numeric time/);
});
