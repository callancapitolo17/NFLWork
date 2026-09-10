// Run: node --test unabated_ticket/tests/kelly.test.js
const test = require("node:test");
const assert = require("node:assert/strict");
const kelly = require("../extension/kelly.js");

const nearly = (actual, expected, tol = 1e-9) =>
  assert.ok(Math.abs(actual - expected) <= tol, `expected ${expected}, got ${actual}`);

test("american -> decimal", () => {
  nearly(kelly.americanToDecimal(-400), 1.25);
  nearly(kelly.americanToDecimal(+150), 2.5);
  nearly(kelly.americanToDecimal(-110), 1 + 100 / 110);
  nearly(kelly.americanToDecimal(100), 2);
});

test("american -> probability", () => {
  nearly(kelly.americanToProb(-900), 0.9);
  nearly(kelly.americanToProb(+300), 0.25);
  nearly(kelly.americanToProb(-100), 0.5);
});

test("rejects prices inside (-100, 100) and non-numbers", () => {
  assert.throws(() => kelly.americanToDecimal(50), /American price/);
  assert.throws(() => kelly.americanToProb(NaN), /finite American price/);
  assert.throws(() => kelly.americanToProb("-110"), /finite American price/);
});

test("Kelly sheet worked example via edge: -400 at +12.5% edge, 30000 bankroll, quarter Kelly = 3750", () => {
  // -400 vs fair -900 is a 12.5% edge in the sheet's Mode A; Mode B takes the edge directly.
  const result = kelly.kellyStakeFromEdge({ bookPrice: -400, edgePct: 12.5, bankroll: 30000, multiplier: 0.25 });
  nearly(result.fullKellyFraction, 0.5);
  nearly(result.stake, 3750);
});

test("Unabated's Seattle example: -133 at +1.89% edge", () => {
  const result = kelly.kellyStakeFromEdge({ bookPrice: -133, edgePct: 1.89, bankroll: 30000, multiplier: 0.25 });
  // b = 100/133 = 0.7519; f = 0.0189 / 0.7519 = 0.025137; stake = 30000 * f * 0.25
  nearly(result.stake, 30000 * 0.25 * (0.0189 / (100 / 133)), 1e-6);
});

test("zero or negative edge gives stake 0", () => {
  assert.equal(kelly.kellyStakeFromEdge({ bookPrice: -110, edgePct: 0, bankroll: 1000, multiplier: 1 }).stake, 0);
  assert.equal(kelly.kellyStakeFromEdge({ bookPrice: -110, edgePct: -5.43, bankroll: 1000, multiplier: 1 }).stake, 0);
});

test("stake is not rounded", () => {
  const result = kelly.kellyStakeFromEdge({ bookPrice: 101, edgePct: 3.3, bankroll: 30000, multiplier: 0.25 });
  assert.notEqual(result.stake, Math.round(result.stake));
});

test("display cents prefer the exchange source price (Novig 0.525 vs rounded -111)", () => {
  nearly(kelly.bookProbOf({ bookPrice: -111, sourceFormat: 4, sourcePrice: 0.525 }), 0.525);
  nearly(kelly.bookProbOf({ bookPrice: -111, sourceFormat: 1, sourcePrice: null }), 111 / 211);
  nearly(kelly.bookProbOf({ bookPrice: -110, sourceFormat: 2, sourcePrice: 1.90909 }), 1 / 1.90909);
});

test("bankroll, multiplier and edge must be valid", () => {
  assert.throws(() => kelly.kellyStakeFromEdge({ bookPrice: -400, edgePct: 1, bankroll: 0, multiplier: 0.25 }), /bankroll/);
  assert.throws(() => kelly.kellyStakeFromEdge({ bookPrice: -400, edgePct: 1, bankroll: 100, multiplier: -1 }), /multiplier/);
  assert.throws(() => kelly.kellyStakeFromEdge({ bookPrice: -400, edgePct: null, bankroll: 100, multiplier: 0.25 }), /edgePct/);
});
