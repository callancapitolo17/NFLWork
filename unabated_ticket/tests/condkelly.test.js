// Run: node --test unabated_ticket/tests/condkelly.test.js
//
// Pins the conditional Kelly numbers agreed in issue #130. Every card was
// measured live 2026-09-17/18 on a Kelly bankroll (bankroll x multiplier) of
// $8,000; probabilities are Unabated's fairs at the time.
const test = require("node:test");
const assert = require("node:assert/strict");
const condkelly = require("../extension/condkelly.js");
const kelly = require("../extension/kelly.js");

const KELLY_BANKROLL = 8000;
const CENT = 0.005;

const toTheCent = (actual, expected) =>
  assert.ok(Math.abs(actual - expected) <= CENT, `expected ${expected}, got ${actual}`);

function netOddsOf(american) {
  return kelly.americanToDecimal(american) - 1;
}

function heldBet(group, cut, direction, stake, american) {
  return { group, cut, direction, stake, toWin: stake * netOddsOf(american) };
}

function candidateBet(group, cut, direction, american, prob) {
  return { group, cut, direction, netOdds: netOddsOf(american), prob };
}

function stakeOf(candidate, held, ladders) {
  const solved = condkelly.solveStake({ kellyBankroll: KELLY_BANKROLL, candidate, held, ladders: ladders || {} });
  assert.equal(solved.reason, null);
  return solved.stake;
}

// UConn @ Southern Miss totals.
const UCONN_PROB_UNDER_52_5 = 1.0465 / 2.25;
const UCONN_LADDER = { FG: [[47.5, 1 - 0.3413], [52.5, 1 - UCONN_PROB_UNDER_52_5]] };

test("no held bets: the stake is kellyStakeFromEdge, to the cent", () => {
  const cases = [
    { bookPrice: 125, edgePct: 4.65 },
    { bookPrice: -133, edgePct: 1.89 },
    { bookPrice: 944, edgePct: 6 },
    { bookPrice: -400, edgePct: 12.5 },
  ];
  for (const { bookPrice, edgePct } of cases) {
    const standalone = kelly.kellyStakeFromEdge({ bookPrice, edgePct, bankroll: 32000, multiplier: 0.25 }).stake;
    const prob = (1 + edgePct / 100) / kelly.americanToDecimal(bookPrice);
    toTheCent(stakeOf(candidateBet("FG", 52.5, "below", bookPrice, prob), []), standalone);
  }
});

test("no held bets: no edge sizes to zero", () => {
  const prob = 0.99 / kelly.americanToDecimal(-110);
  assert.equal(stakeOf(candidateBet("FG", 52.5, "above", -110, prob), []), 0);
});

test("same line held: full size minus what is held, and zero once over size", () => {
  const candidate = candidateBet("FG", 52.5, "below", 125, UCONN_PROB_UNDER_52_5);
  toTheCent(stakeOf(candidate, [heldBet("FG", 52.5, "below", 181.5, 125)]), 116.10);
  assert.equal(stakeOf(candidate, [heldBet("FG", 52.5, "below", 400, 125)]), 0);
});

test("a held bet on the candidate's own cut needs no ladder rung", () => {
  const candidate = candidateBet("FG", 52.5, "below", 125, UCONN_PROB_UNDER_52_5);
  const solved = condkelly.solveStake({ kellyBankroll: KELLY_BANKROLL, candidate, held: [heldBet("FG", 52.5, "above", 100, -110)], ladders: {} });
  assert.equal(solved.reason, null);
  assert.ok(solved.stake > 0);
});

test("UConn: an easier number on the same side held, add $121.28", () => {
  const candidate = candidateBet("FG", 52.5, "below", 125, UCONN_PROB_UNDER_52_5);
  toTheCent(stakeOf(candidate, [heldBet("FG", 47.5, "below", 181.5, 203)], UCONN_LADDER), 121.28);
});

test("UConn: a harder number on the same side, add $24.10", () => {
  const candidate = candidateBet("FG", 47.5, "below", 203, 0.3413);
  toTheCent(stakeOf(candidate, [heldBet("FG", 52.5, "below", 181.5, 125)], UCONN_LADDER), 24.10);
});

test("New Mexico @ Oklahoma: two unders held against the over, bet $295.99", () => {
  const candidate = candidateBet("FG", 53.5, "above", 217, 1.0497 / 3.17);
  const held = [heldBet("FG", 39.5, "below", 214, 239), heldBet("FG", 48.5, "below", 6.2, 108)];
  const ladders = { FG: [[39.5, 1 - 100 / 300], [48.5, 1 - 137 / 237]] };
  toTheCent(stakeOf(candidate, held, ladders), 295.99);
});

test("Lions @ Bills: the same over held plus two unders at another number, add $188.32", () => {
  const candidate = candidateBet("FG", 61.5, "above", 213, 1.0719 / 3.13);
  const held = [
    heldBet("FG", 61.5, "above", 270, 212),
    heldBet("FG", 51.5, "below", 352, 127),
    heldBet("FG", 51.5, "below", 61, 150),
  ];
  toTheCent(stakeOf(candidate, held, { FG: [[51.5, 0.5745]] }), 188.32);
});

test("Chargers: a moneyline with the +3.5 already held adds nothing", () => {
  // Margin axis = away minus home, Chargers away: +3.5 wins above -3.5, the moneyline above +0.5.
  const candidate = candidateBet("FG", 0.5, "above", 213, 1.0503 / 3.13);
  assert.equal(stakeOf(candidate, [heldBet("FG", -3.5, "above", 400, 122)], { FG: [[-3.5, 0.4785]] }), 0);
});

test("cross-period: a 1H over held pairs worst-case with the FG over, $408.39 not $666.67", () => {
  const candidate = candidateBet("FG", 48.5, "above", 120, 0.5);
  toTheCent(stakeOf(candidate, []), 666.67);
  toTheCent(stakeOf(candidate, [heldBet("1H", 24.5, "above", 300, 110)], { "1H": [[24.5, 0.55]] }), 408.39);
});

test("cross-period: a bet on the other side in another period earns no hedge credit", () => {
  const candidate = candidateBet("FG", 48.5, "above", 120, 0.5);
  const stake = stakeOf(candidate, [heldBet("1H", 24.5, "below", 300, 110)], { "1H": [[24.5, 0.55]] });
  assert.ok(stake < 666.67, `worst case must not size above the standalone stake, got ${stake}`);
});

test("whole number: the push row comes from the two neighbouring half-point rungs", () => {
  // Over 52 held: wins above 52.5, pushes on 52, loses below 51.5.
  const candidate = candidateBet("FG", 52.5, "above", 100, 0.52);
  const ladders = { FG: [[51.5, 0.56]] };
  const withPush = stakeOf(candidate, [heldBet("FG", 52, "above", 100, 100)], ladders);
  const noPush = stakeOf(candidate, [heldBet("FG", 52.5, "above", 100, 100)], ladders);
  toTheCent(noPush, KELLY_BANKROLL * 0.04 - 100);
  assert.ok(withPush > noPush, `a bet that pushes 4% of the time is less held than one that loses, got ${withPush} vs ${noPush}`);
  assert.deepEqual(condkelly.cutsNeeded({ cut: 52, direction: "above" }), [51.5, 52.5]);
  assert.deepEqual(condkelly.cutsNeeded({ cut: 52.5, direction: "below" }), [52.5]);
});

test("whole-number candidate with nothing held is still the standalone stake", () => {
  const prob = 1.04 / 2;
  const ladders = { FG: [[51.5, 0.55], [52.5, 0.49]] };
  toTheCent(stakeOf(candidateBet("FG", 52, "above", 100, prob), [], ladders), KELLY_BANKROLL * 0.04);
});

test("a ladder crossing by up to half a point of probability is clamped, past that the calc is declined", () => {
  const candidate = candidateBet("FG", 52.5, "above", 100, 0.52);
  const held = [heldBet("FG", 51.5, "above", 100, -110)];
  const clamped = condkelly.solveStake({ kellyBankroll: KELLY_BANKROLL, candidate, held, ladders: { FG: [[51.5, 0.516]] } });
  assert.equal(clamped.reason, null);
  assert.ok(clamped.stake >= 0);
  const declined = condkelly.solveStake({ kellyBankroll: KELLY_BANKROLL, candidate, held, ladders: { FG: [[51.5, 0.51]] } });
  assert.deepEqual(declined, { stake: null, reason: condkelly.REASON_NOT_MONOTONE });
});

test("held bets that can already lose the Kelly bankroll decline the calc", () => {
  const candidate = candidateBet("FG", 52.5, "above", 100, 0.52);
  const solved = condkelly.solveStake({
    kellyBankroll: KELLY_BANKROLL, candidate, held: [heldBet("FG", 60.5, "above", 9000, 300)], ladders: { FG: [[60.5, 0.25]] },
  });
  assert.deepEqual(solved, { stake: null, reason: condkelly.REASON_HELD_RISK });
});

test("fails loudly on a missing rung, a quarter line and a bad probability", () => {
  const candidate = candidateBet("FG", 52.5, "above", 100, 0.52);
  assert.throws(() => condkelly.solveStake({ kellyBankroll: KELLY_BANKROLL, candidate, held: [heldBet("FG", 47.5, "below", 100, 200)], ladders: {} }),
    /expected a ladder probability for FG at 47.5/);
  assert.throws(() => condkelly.solveStake({ kellyBankroll: KELLY_BANKROLL, candidate: { ...candidate, cut: 52.25 }, held: [] }),
    /whole or half-point number, got 52.25/);
  assert.throws(() => condkelly.solveStake({ kellyBankroll: KELLY_BANKROLL, candidate: { ...candidate, prob: 1.2 }, held: [] }),
    /inside \(0, 1\), got 1.2/);
  assert.throws(() => condkelly.solveStake({ kellyBankroll: 0, candidate, held: [] }), /kellyBankroll: expected a positive number/);
});
