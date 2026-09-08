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
  nearly(kelly.americanToProb(100), 0.5);
});

test("rejects prices inside (-100, 100) and non-numbers", () => {
  assert.throws(() => kelly.americanToDecimal(50), /American price/);
  assert.throws(() => kelly.americanToProb(NaN), /finite American price/);
  assert.throws(() => kelly.americanToProb("-110"), /finite American price/);
});

test("Kelly sheet worked example: -400 vs fair -900, 30000 bankroll, quarter Kelly = 3750", () => {
  const result = kelly.kellyStake({
    bookPrice: -400,
    fairPrice: -900,
    bankroll: 30000,
    multiplier: 0.25,
  });
  nearly(result.fullKellyFraction, 0.5);
  nearly(result.fullKellyStake, 15000);
  nearly(result.stake, 3750);
  nearly(result.edge, 0.9 * 1.25 - 1); // +12.5% per $1
});

test("negative edge gives stake 0 and a negative edge fraction", () => {
  const result = kelly.kellyStake({
    bookPrice: -49900,
    fairPrice: -5854,
    bankroll: 30000,
    multiplier: 0.25,
  });
  assert.equal(result.stake, 0);
  assert.equal(result.fullKellyStake, 0);
  assert.equal(result.fullKellyFraction, 0);
  assert.ok(result.edge < 0);
});

test("zero edge (book price equals fair) gives stake 0", () => {
  const result = kelly.kellyStake({ bookPrice: -110, fairPrice: -110, bankroll: 1000, multiplier: 1 });
  assert.equal(result.stake, 0);
  nearly(result.edge, 0);
});

test("stake is not rounded", () => {
  // +101 vs fair -106: p = 106/206, b = 1.01
  const result = kelly.kellyStake({ bookPrice: 101, fairPrice: -106, bankroll: 30000, multiplier: 0.25 });
  const p = 106 / 206;
  const b = 1.01;
  const expectedFull = (p * b - (1 - p)) / b;
  nearly(result.stake, 30000 * expectedFull * 0.25);
  assert.notEqual(result.stake, Math.round(result.stake));
});

test("exchange source price beats the rounded American (Novig 0.525 vs -111)", () => {
  // Unabated rounds 0.525 -> -110.5 -> -111; converting back gives 52.6c. The exchange said 52.5c.
  nearly(kelly.bookProbOf({ bookPrice: -111, sourceFormat: 4, sourcePrice: 0.525 }), 0.525);
  nearly(kelly.bookProbOf({ bookPrice: -111, sourceFormat: 1, sourcePrice: null }), 111 / 211);
  nearly(kelly.bookDecimalOf({ bookPrice: -110, sourceFormat: 2, sourcePrice: 1.90909 }), 1.90909);
  const exact = kelly.kellyStake({ bookPrice: 105, sourceFormat: 4, sourcePrice: 0.48, fairPrice: -101, bankroll: 30000, multiplier: 0.25 });
  const rounded = kelly.kellyStake({ bookPrice: 105, fairPrice: -101, bankroll: 30000, multiplier: 0.25 });
  assert.ok(exact.stake > rounded.stake, "48.0c pays better than +105 (48.8c), so the exact stake is larger");
});

test("bankroll and multiplier must be positive", () => {
  assert.throws(() => kelly.kellyStake({ bookPrice: -400, fairPrice: -900, bankroll: 0, multiplier: 0.25 }), /bankroll/);
  assert.throws(() => kelly.kellyStake({ bookPrice: -400, fairPrice: -900, bankroll: 100, multiplier: -1 }), /multiplier/);
});
