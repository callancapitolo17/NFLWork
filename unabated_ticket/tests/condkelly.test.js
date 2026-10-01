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

test("resultAt: win, push and loss of a bet at a result between its cuts", () => {
  const over52 = { cut: 52, direction: "above" };
  assert.deepEqual([51.25, 51.75, 52.75].map((value) => condkelly.resultAt(over52, value)), [-1, 0, 1]);
  const homeMinus3 = { cut: -3, direction: "below" };
  assert.deepEqual([-3.75, -3.25, -2.25].map((value) => condkelly.resultAt(homeMinus3, value)), [1, 0, -1]);
  const under44Half = { cut: 44.5, direction: "below" };
  assert.deepEqual([44.25, 44.75].map((value) => condkelly.resultAt(under44Half, value)), [1, -1]);
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

// ---- open teasers (2026-09-30, teasers plan section 14) ----------------------
//
// A ticket's leg on this game cuts its rows like a held bet (a half-point: a
// teaser push loses); its other legs ride on other games, enumerated.

function ticketLeg(group, cut, direction, stake, toWin, others) {
  return { group, cut, direction, stake, toWin, others: others || [] };
}

// The log growth written out by hand, [chance, P&L at stake x] per outcome,
// topped by a search over every cent.
function bestStakeByGrid(kellyBankroll, outcomes, maxStake) {
  let best = { stake: 0, growth: -Infinity };
  for (let cents = 0; cents <= maxStake * 100; cents += 1) {
    const stake = cents / 100;
    let growth = 0;
    for (const [chance, pnlAt] of outcomes) growth += chance * Math.log(1 + pnlAt(stake) / kellyBankroll);
    if (growth > best.growth) best = { stake, growth };
  }
  return best.stake;
}

test("a ticket whose other legs are all in plays like a straight at its half-point", () => {
  const candidate = candidateBet("FG", 52.5, "below", 125, UCONN_PROB_UNDER_52_5);
  const asTicket = condkelly.solveStake({ kellyBankroll: KELLY_BANKROLL, candidate, held: [], tickets: [ticketLeg("FG", 47.5, "below", 200, 600)], ladders: UCONN_LADDER });
  const asStraight = condkelly.solveStake({ kellyBankroll: KELLY_BANKROLL, candidate, held: [{ group: "FG", cut: 47.5, direction: "below", stake: 200, toWin: 600 }], ladders: UCONN_LADDER });
  assert.equal(asTicket.reason, null);
  toTheCent(asTicket.stake, asStraight.stake);
});

test("one other leg: the stake tops the growth averaged over that leg, not its average P&L", () => {
  const prob = UCONN_PROB_UNDER_52_5;
  const otherWins = 0.7;
  const candidate = candidateBet("FG", 52.5, "below", 125, prob);
  const tickets = [ticketLeg("FG", 52.5, "below", 200, 600, [{ factor: "g2", cut: 0.5, direction: "above" }])];
  const solved = condkelly.solveStake({ kellyBankroll: KELLY_BANKROLL, candidate, held: [], tickets, factors: { g2: [[0.5, otherWins]] }, ladders: {} });
  assert.equal(solved.reason, null);
  const byHand = bestStakeByGrid(KELLY_BANKROLL, [
    [prob * otherWins, (x) => 600 + 1.25 * x],
    [prob * (1 - otherWins), (x) => -200 + 1.25 * x],
    [1 - prob, (x) => -200 - x],
  ], 1500);
  assert.ok(Math.abs(solved.stake - byHand) <= 0.01, `expected ${byHand}, got ${solved.stake}`);
  const averaged = bestStakeByGrid(KELLY_BANKROLL, [
    [prob, (x) => otherWins * 600 - (1 - otherWins) * 200 + 1.25 * x],
    [1 - prob, (x) => -200 - x],
  ], 1500);
  assert.ok(Math.abs(solved.stake - averaged) > 1, `the gamble must not size like its average, both ${solved.stake}`);
});

test("legs shared across tickets are joint: two tickets on one other game are not two independent ones", () => {
  // The Over against two Under tickets, so both ways size above zero.
  const underWins = UCONN_PROB_UNDER_52_5;
  const candidate = candidateBet("FG", 52.5, "above", 100, 1 - underWins);
  const leg = (factor) => [{ factor, cut: 0.5, direction: "above" }];
  const shared = condkelly.solveStake({
    kellyBankroll: KELLY_BANKROLL, candidate, held: [], ladders: {}, factors: { g2: [[0.5, 0.7]] },
    tickets: [ticketLeg("FG", 52.5, "below", 200, 600, leg("g2")), ticketLeg("FG", 52.5, "below", 200, 600, leg("g2"))],
  });
  const apart = condkelly.solveStake({
    kellyBankroll: KELLY_BANKROLL, candidate, held: [], ladders: {}, factors: { g2: [[0.5, 0.7]], g3: [[0.5, 0.7]] },
    tickets: [ticketLeg("FG", 52.5, "below", 200, 600, leg("g2")), ticketLeg("FG", 52.5, "below", 200, 600, leg("g3"))],
  });
  const sharedByHand = bestStakeByGrid(KELLY_BANKROLL, [
    [underWins * 0.7, (x) => 1200 - x],
    [underWins * 0.3, (x) => -400 - x],
    [1 - underWins, (x) => -400 + x],
  ], 1500);
  const apartByHand = bestStakeByGrid(KELLY_BANKROLL, [
    [underWins * 0.49, (x) => 1200 - x],
    [underWins * 0.42, (x) => 400 - x],
    [underWins * 0.09, (x) => -400 - x],
    [1 - underWins, (x) => -400 + x],
  ], 1500);
  assert.ok(Math.abs(shared.stake - sharedByHand) <= 0.01, `shared: expected ${sharedByHand}, got ${shared.stake}`);
  assert.ok(Math.abs(apart.stake - apartByHand) <= 0.01, `apart: expected ${apartByHand}, got ${apart.stake}`);
  assert.ok(Math.abs(shared.stake - apart.stake) > 1, `joint and independent legs must size apart, both ${shared.stake}`);
});

// Sunday 9/27, 16:46 UTC board, Cal's real BFA teasers, K = $20,000 x 0.25:
// the walk-through that settled the design. Each other leg is Unabated's
// fair at its teased number, a whole American price on the leg's side.
const TEASER_KELLY_BANKROLL = 5000;
const TEASER_LEGS = {
  Titans: { cut: -8.5, direction: "above", american: -295 },
  "49ers": { cut: -1.5, direction: "below", american: -336 },
  Seahawks: { cut: 2.5, direction: "above", american: -317 },
  Browns: { cut: 7.5, direction: "below", american: -283 },
  Lions: { cut: -1.5, direction: "below", american: -267 },
  Bills: { cut: -1.5, direction: "below", american: -285 },
  Broncos: { cut: 7.5, direction: "below", american: -283 },
  Colts: { cut: 7.5, direction: "below", american: -295 },
};

function teaserOn(cut, direction, otherTeams) {
  const others = otherTeams.map((team) => ({ factor: team, cut: TEASER_LEGS[team].cut, direction: TEASER_LEGS[team].direction }));
  return ticketLeg("FG", cut, direction, 200, 600, others);
}

function teaserFactors() {
  const factors = {};
  for (const [team, { cut, direction, american }] of Object.entries(TEASER_LEGS)) {
    const wins = kelly.americanToProb(american);
    factors[team] = [[cut, direction === "above" ? wins : 1 - wins]];
  }
  return factors;
}

test("9/27 Colts +7.5 -270 at Novig: $315.90 alone, $0 with the three Colts +7.5 teasers", () => {
  const candidate = candidateBet("FG", 7.5, "below", -270, 1.0234 / kelly.americanToDecimal(-270));
  const alone = condkelly.solveStake({ kellyBankroll: TEASER_KELLY_BANKROLL, candidate, held: [], ladders: {} });
  toTheCent(alone.stake, 315.90);
  const tickets = [
    teaserOn(7.5, "below", ["Titans", "49ers", "Seahawks"]),
    teaserOn(7.5, "below", ["Titans", "Browns", "Lions"]),
    teaserOn(7.5, "below", ["Lions", "Seahawks", "Bills"]),
  ];
  const solved = condkelly.solveStake({ kellyBankroll: TEASER_KELLY_BANKROLL, candidate, held: [], tickets, factors: teaserFactors(), ladders: {} });
  assert.deepEqual(solved, { stake: 0, reason: null });
});

test("9/27 Chargers +6.5 +120 at Novig: $137.48 on the straights, $425.65 with the four Bills -1 teasers against", () => {
  const candidate = candidateBet("FG", -6.5, "above", 120, 1.0233 / kelly.americanToDecimal(120));
  const held = [
    { group: "FG", cut: -3.5, direction: "above", stake: 400, toWin: 488.88 },
    { group: "FG", cut: -14.5, direction: "below", stake: 375, toWin: 988.63 },
    { group: "FG", cut: -14.5, direction: "below", stake: 70, toWin: 179.99 },
  ];
  // Chargers @ Bills: P(margin above the cut) = the Chargers at +3.5 (+173), +14.5 (-235), +1.5 (+285).
  const ladders = { FG: [[-3.5, kelly.americanToProb(173)], [-14.5, kelly.americanToProb(-235)], [-1.5, kelly.americanToProb(285)]] };
  const straightsOnly = condkelly.solveStake({ kellyBankroll: TEASER_KELLY_BANKROLL, candidate, held, ladders });
  toTheCent(straightsOnly.stake, 137.48);
  const tickets = [
    teaserOn(-1.5, "below", ["Browns", "49ers", "Seahawks"]),
    teaserOn(-1.5, "below", ["Lions", "49ers", "Broncos"]),
    teaserOn(-1.5, "below", ["Lions", "Seahawks", "Colts"]),
    teaserOn(-1.5, "below", ["Broncos", "Titans", "Browns"]),
  ];
  const solved = condkelly.solveStake({ kellyBankroll: TEASER_KELLY_BANKROLL, candidate, held, tickets, factors: teaserFactors(), ladders });
  assert.equal(solved.reason, null);
  toTheCent(solved.stake, 425.65);
});

test("too many other games under the tickets decline the calc", () => {
  const candidate = candidateBet("FG", 52.5, "below", 125, UCONN_PROB_UNDER_52_5);
  const factors = {};
  const tickets = [];
  for (let game = 0; game < 17; game += 1) {
    factors[`g${game}`] = [[0.5, 0.7]];
    tickets.push(ticketLeg("FG", 52.5, "below", 10, 30, [{ factor: `g${game}`, cut: 0.5, direction: "above" }]));
  }
  const solved = condkelly.solveStake({ kellyBankroll: KELLY_BANKROLL, candidate, held: [], tickets, factors, ladders: {} });
  assert.deepEqual(solved, { stake: null, reason: condkelly.REASON_TOO_MANY_GAMES });
});

test("an other game's ladder rising with the cut declines the calc", () => {
  const candidate = candidateBet("FG", 52.5, "below", 125, UCONN_PROB_UNDER_52_5);
  const tickets = [
    ticketLeg("FG", 52.5, "below", 200, 600, [{ factor: "g2", cut: 0.5, direction: "above" }]),
    ticketLeg("FG", 52.5, "below", 200, 600, [{ factor: "g2", cut: 1.5, direction: "above" }]),
  ];
  const solved = condkelly.solveStake({ kellyBankroll: KELLY_BANKROLL, candidate, held: [], tickets, factors: { g2: [[0.5, 0.5], [1.5, 0.6]] }, ladders: {} });
  assert.deepEqual(solved, { stake: null, reason: condkelly.REASON_NOT_MONOTONE });
});

test("fails loudly on a whole-number ticket leg and a missing other-game rung", () => {
  const candidate = candidateBet("FG", 52.5, "below", 125, UCONN_PROB_UNDER_52_5);
  assert.throws(() => condkelly.solveStake({ kellyBankroll: KELLY_BANKROLL, candidate, held: [], tickets: [ticketLeg("FG", 52, "below", 200, 600)], ladders: {} }),
    /tickets\[0\]\.cut: expected a half-point number \(a teaser push loses\), got 52/);
  assert.throws(() => condkelly.solveStake({
    kellyBankroll: KELLY_BANKROLL, candidate, held: [], ladders: {}, factors: {},
    tickets: [ticketLeg("FG", 52.5, "below", 200, 600, [{ factor: "g2", cut: -1.5, direction: "below" }])],
  }), /expected a ladder probability for g2 at -1.5/);
});
