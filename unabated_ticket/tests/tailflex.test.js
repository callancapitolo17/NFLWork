// Run: node --test unabated_ticket/tests
// tailflex.js: the probit devig, the live measurement of c off two-sided
// exchange rungs, the tau samples and their rolling pool, and the shrink
// that sizes the stake and ranks the lines (issue #208). The worked card is NFL Chargers @ Seahawks (2026-09-30
// design session): Seahawks main fair 48.1%, c = 10%, $10k bankroll, quarter Kelly.
const test = require("node:test");
const assert = require("node:assert/strict");
const tailflex = require("../extension/tailflex.js");
const feed = require("../extension/feed.js");
const kelly = require("../extension/kelly.js");

const NOW = Date.parse("2026-09-30T12:00:00Z");
const KICKOFF = NOW + 6 * 3600 * 1000;
const NFL = 1;
const FULL_GAME = 1;
const SPREAD = 2;
const EXCHANGE_BOOK = 89;
const BANKROLL = 10000;
const QUARTER_KELLY = 0.25;

function americanOf(prob) {
  return prob >= 0.5 ? -Math.round((100 * prob) / (1 - prob)) : Math.round((100 * (1 - prob)) / prob);
}

// ---- normal distribution and devig -----------------------------------------

test("qnorm inverts pnorm; dnorm(0) and the one-cent fuselage", () => {
  for (const p of [0.001, 0.02, 0.1, 0.308, 0.481, 0.5, 0.75, 0.99]) {
    assert.ok(Math.abs(tailflex.pnorm(tailflex.qnorm(p)) - p) < 1e-6, `p = ${p}`);
  }
  assert.ok(Math.abs(tailflex.qnorm(0.975) - 1.959964) < 1e-5);
  assert.ok(Math.abs(tailflex.dnorm(0) - 0.398942) < 1e-6);
  assert.ok(Math.abs(tailflex.DZ_MAIN - 0.025066) < 1e-6);
  assert.throws(() => tailflex.qnorm(0), /strictly between 0 and 1/);
});

test("devigProbitTwoWay: symmetric vig comes off exactly; the sides sum to one", () => {
  const fairZ = -0.6;
  const vig = 0.04;
  const implied0 = tailflex.pnorm(fairZ + vig);
  const implied1 = tailflex.pnorm(-fairZ + vig);
  // pnorm is good to 1.5e-7, so "exactly" is to 1e-6.
  assert.ok(Math.abs(tailflex.devigProbitTwoWay(implied0, implied1) - tailflex.pnorm(fairZ)) < 1e-6);
  assert.ok(Math.abs(tailflex.devigProbitTwoWay(0.55, 0.5) + tailflex.devigProbitTwoWay(0.5, 0.55) - 1) < 1e-6);
  assert.ok(Math.abs(tailflex.devigProbitTwoWay(0.52, 0.52) - 0.5) < 1e-6);
});

test("impliedProbOf reads the exchange's exact figure before the rounded American price", () => {
  assert.equal(tailflex.impliedProbOf({ price: 270, sourceFormat: 4, sourcePrice: 0.27 }), 0.27);
  assert.equal(tailflex.impliedProbOf({ price: 270, sourceFormat: 2, sourcePrice: 4 }), 0.25);
  assert.ok(Math.abs(tailflex.impliedProbOf({ price: -110, sourceFormat: 1, sourcePrice: null }) - 110 / 210) < 1e-12);
  assert.equal(tailflex.impliedProbOf({ price: 50, sourceFormat: 1, sourcePrice: null }), null);
});

// ---- measuring c ------------------------------------------------------------

// A state of `events` NFL games, each with Unabated's own main line at -3
// (fair 50%) and one exchange quoting both sides of the main number plus
// five alt rungs. Unabated's side-0 fair at points x is z = 0.08 * (x + 3);
// the exchange's devigged fair sits sqrt((cTrue * dist)^2 + baseline^2) away
// (baseline alone at the main number), under symmetric probit vig.
function syntheticState({ events, cTrue, baseline = 0, vig = 0.02, overrides = () => ({}) }) {
  const state = feed.emptyState();
  state.books = { [EXCHANGE_BOOK]: { id: EXCHANGE_BOOK, hasLiquidity: true }, [feed.UNABATED_LINE_BOOK_ID]: { id: feed.UNABATED_LINE_BOOK_ID } };
  const add = (line) => {
    state.lines[`${line.eventId}:ms${line.bookId}:si${line.sideIndex}:${line.points}`] = {
      leagueId: NFL, periodTypeId: FULL_GAME, betTypeId: SPREAD, statusId: 1, isAlt: false,
      sequenceNumber: NOW - 60000, modifiedOn: new Date(NOW - 60000).toISOString(), ...line,
    };
  };
  for (let eventId = 1; eventId <= events; eventId += 1) {
    state.events[eventId] = { eventId, leagueId: NFL, eventStart: KICKOFF };
    add({ eventId, bookId: feed.UNABATED_LINE_BOOK_ID, sideIndex: 0, points: -3, price: 100, bacr: 100 });
    add({ eventId, bookId: feed.UNABATED_LINE_BOOK_ID, sideIndex: 1, points: 3, price: -100, bacr: -100 });
    for (const points of [-3, -7, -10, -14, -17, -21]) {
      const bacr = americanOf(tailflex.pnorm(0.08 * (points + 3)));
      const unabatedZ = tailflex.qnorm(kelly.americanToProb(bacr));
      const distSd = Math.abs(unabatedZ);
      const gap = points === -3 ? baseline : Math.sqrt((cTrue * distSd) ** 2 + baseline ** 2);
      const exchangeZ = unabatedZ + gap;
      const extra = overrides({ eventId, points });
      add({ eventId, bookId: EXCHANGE_BOOK, sideIndex: 0, points, bacr, price: 100, sourceFormat: 4,
        sourcePrice: tailflex.pnorm(exchangeZ + vig), isAlt: points !== -3, ...(extra.side0 || {}) });
      if (extra.dropSide1) continue;
      add({ eventId, bookId: EXCHANGE_BOOK, sideIndex: 1, points: -points, bacr: -bacr, price: 100, sourceFormat: 4,
        sourcePrice: tailflex.pnorm(-exchangeZ + vig), isAlt: points !== -3, ...(extra.side1 || {}) });
    }
  }
  return state;
}

const NFL_SPREAD_FG = tailflex.marketKeyOf({ leagueId: NFL, periodTypeId: FULL_GAME, betTypeId: SPREAD });

test("measureTailFlex recovers c from 125 two-sided rungs, net of the main-number baseline", () => {
  const measurement = tailflex.measureTailFlex(syntheticState({ events: 25, cTrue: 0.062, baseline: 0.02 }), { now: NOW });
  const market = measurement.markets[NFL_SPREAD_FG];
  assert.equal(market.measured, true);
  assert.equal(market.rungCount, 125);
  assert.ok(Math.abs(market.c - 0.062) < 1e-5, `c = ${market.c}`);
  assert.ok(Math.abs(market.baseline - 0.02) < 1e-5);
  assert.equal(tailflex.cOf(measurement, { leagueId: NFL, periodTypeId: FULL_GAME, betTypeId: SPREAD }), market.c);
});

test("fewer than 100 rungs past 0.3 SD falls back to 10%; an unmeasured market falls back too", () => {
  const measurement = tailflex.measureTailFlex(syntheticState({ events: 19, cTrue: 0.06 }), { now: NOW });
  const market = measurement.markets[NFL_SPREAD_FG];
  assert.equal(market.rungCount, 95);
  assert.equal(market.measured, false);
  assert.equal(market.c, tailflex.C_FALLBACK);
  assert.equal(tailflex.cOf(measurement, { leagueId: 2, periodTypeId: FULL_GAME, betTypeId: SPREAD }), 0.10);
  assert.equal(tailflex.isMeasured(measurement, { leagueId: 2, periodTypeId: FULL_GAME, betTypeId: SPREAD }), false);
});

test("cOfRungs: an RMS, not a median; rungs inside 0.3 SD never divide; the top 1% is dropped; noise below baseline is zero", () => {
  const rungs = [{ onMain: true, gap: 0.03, distSd: 0 }];
  for (let i = 0; i < 99; i += 1) rungs.push({ onMain: false, gap: Math.hypot(0.05 * 0.5, 0.03), distSd: 0.5 });
  rungs.push({ onMain: false, gap: 0.5, distSd: 0.05 }); // ratio 10 if it counted
  rungs.push({ onMain: false, gap: 2, distSd: 0.6 }); // a stale quote: the top 1% of 100 ratios, dropped
  const result = tailflex.cOfRungs(rungs);
  assert.equal(result.rungCount, 100);
  assert.ok(Math.abs(result.c - 0.05) < 1e-12, `c = ${result.c}`);
  // Half the rungs at 0.02 and half at 0.10: the median would say about 0.06, the RMS sqrt((0.02^2 + 0.10^2) / 2).
  const mixed = Array.from({ length: 200 }, (_, i) => ({ onMain: false, gap: i % 2 ? 0.02 : 0.10, distSd: 1 }));
  assert.ok(Math.abs(tailflex.cOfRungs(mixed).c - Math.sqrt((0.02 ** 2 + 0.10 ** 2) / 2)) < 0.001);
  const quiet = tailflex.cOfRungs([{ onMain: true, gap: 0.05, distSd: 0 }, ...Array.from({ length: 100 }, () => ({ onMain: false, gap: 0.01, distSd: 1 }))]);
  assert.equal(quiet.c, 0);
});

test("only fresh, two-sided, sanely-vigged rungs of unstarted games measure", () => {
  const base = { events: 25, cTrue: 0.06 };
  const count = (state, opts) => (tailflex.measureTailFlex(state, { now: NOW, ...opts }).markets[NFL_SPREAD_FG] || { rungCount: 0 }).rungCount;
  // One-sided quotes on the first five games: not measured.
  assert.equal(count(syntheticState({ ...base, overrides: ({ eventId }) => (eventId <= 5 ? { dropSide1: true } : {}) })), 100);
  // A crossed pair (sum 0.90 < 0.995) and a 24% overround are out.
  const bothSides = (sourcePrice) => ({ eventId }) => (eventId === 1 ? { side0: { sourcePrice }, side1: { sourcePrice } } : {});
  assert.equal(count(syntheticState({ ...base, overrides: bothSides(0.45) })), 120);
  assert.equal(count(syntheticState({ ...base, overrides: bothSides(0.62) })), 120);
  // Off the board or older than the line-age gate.
  assert.equal(count(syntheticState({ ...base, overrides: ({ eventId }) => (eventId === 1 ? { side1: { statusId: 2 } } : {}) })), 120);
  const stale = syntheticState({ ...base, overrides: ({ eventId }) => (eventId === 1 ? { side0: { sequenceNumber: NOW - 3 * 3600 * 1000, modifiedOn: new Date(NOW - 3 * 3600 * 1000).toISOString() } } : {}) });
  assert.equal(count(stale, { maxLineAgeMs: 3600 * 1000 }), 120);
  assert.equal(count(stale), 125);
  // A game that has started measures nothing.
  assert.equal(count(syntheticState(base), { now: KICKOFF + 1 }), 0);
});

// ---- the rank score -----------------------------------------------------------

const MEASUREMENT_C10 = {
  markets: {},
  mainFairs: new Map([["777:1:2:1", { points: -8.5, fairProb: 0.481 }]]),
};

function seahawksRow({ key, points, price, edgePct, book }) {
  const row = {
    key, eventId: 777, leagueId: NFL, periodTypeId: FULL_GAME, betTypeId: SPREAD, sideIndex: 1,
    points, price, edgePct, book: { id: book, name: `book ${book}` }, eventStartMs: KICKOFF, sideLabel: `Seattle Seahawks ${points}`,
  };
  return { ...row, stake: kelly.kellyStakeFromEdge({ bookPrice: price, edgePct, bankroll: BANKROLL, multiplier: QUARTER_KELLY }).stake };
}

const MINUS_14_5 = seahawksRow({ key: "novig-14.5", points: -14.5, price: 270, edgePct: 13.9, book: 89 });
const MINUS_27_5 = seahawksRow({ key: "hardrock-27.5", points: -27.5, price: 1300, edgePct: 17.3, book: 87 });
const MINUS_2_5 = seahawksRow({ key: "novig-2.5", points: -2.5, price: -257, edgePct: 1.0, book: 89 });

test("shrinkFactor: the share of a claimed edge to believe does not depend on the claimed edge; it falls with distance and with c", () => {
  const decimal = kelly.americanToDecimal(270);
  const fairProb = 1.139 / decimal;
  assert.ok(Math.abs(fairProb - 0.308) < 0.001);
  const zMain = tailflex.qnorm(0.481);
  const { shrink, sigmaE } = tailflex.shrinkFactor({ decimal, fairProb, zMain, c: 0.10, tau: 0.0185 });
  // sigma_e = 3.7 x dnorm(z) x sqrt(DZ_MAIN^2 + (0.1 x 0.45)^2) ~ 6.6% EV against tau 1.85%: ~7% believed.
  assert.equal(sigmaE, tailflex.edgeUncertainty({ decimal, fairProb, zMain, c: 0.10 }));
  assert.ok(Math.abs(sigmaE - 0.066) < 0.002, `sigmaE = ${sigmaE}`);
  assert.ok(Math.abs(shrink - 0.0185 ** 2 / (0.0185 ** 2 + sigmaE ** 2)) < 1e-12);
  assert.ok(shrink > 0.06 && shrink < 0.08, `shrink = ${shrink}`);
  // A wider c flexes more; a bigger tau believes more; at the main number only the one-cent fuselage is left (c drops out).
  assert.ok(tailflex.shrinkFactor({ decimal, fairProb, zMain, c: 0.2, tau: 0.0185 }).shrink < shrink);
  assert.ok(tailflex.shrinkFactor({ decimal, fairProb, zMain, c: 0.1, tau: 0.03 }).shrink > shrink);
  const atMain = tailflex.shrinkFactor({ decimal, fairProb, zMain: tailflex.qnorm(fairProb), c: 0.10, tau: 0.0185 }).shrink;
  assert.equal(atMain, tailflex.shrinkFactor({ decimal, fairProb, zMain: tailflex.qnorm(fairProb), c: 0, tau: 0.0185 }).shrink);
  // At the main, sigma_e = decimal x dnorm(z) x one cent: 3.7 x 0.351 x 0.025 = 3.3% against 1.85%: ~24%.
  assert.ok(atMain > 0.23 && atMain < 0.25, `atMain = ${atMain}`);
});

// A row ranked as edgerows.withStakeAndRank ranks it, on the stake it carries.
function ranked(row, measurement) {
  const shrunk = tailflex.shrinkOfRow(row, measurement);
  return { ...row, shrunkEdgePct: shrunk ? shrunk.edgePct : null, rankScore: tailflex.rankScoreOf({ ...row, shrunkEdgePct: shrunk ? shrunk.edgePct : null }) };
}

test("shrinkOfRow and rankScoreOf: the stake reads the shrunk edge, the rank is shrunk edge x stake; the issue's card picks the main", () => {
  // Issue #208's card: a CFB favorite, main -14 at 51% fair, c 9.4%, $35k at quarter Kelly, tau 1.85%.
  const main = { points: -14, price: 105, edgePct: 4.5 };
  const deep = { points: -24.5, price: 340, edgePct: 12.2 };
  const measurement = {
    markets: { "2:1:2": { c: 0.094, measured: true, tau: 0.0185, tauMeasured: false } },
    mainFairs: new Map([["55:1:2:0", { points: -14, fairProb: 0.51 }]]),
  };
  const rowOf = (line) => ({ key: `k${line.points}`, eventId: 55, leagueId: 2, periodTypeId: FULL_GAME, betTypeId: SPREAD, sideIndex: 0, ...line });
  const shrunkMain = tailflex.shrinkOfRow(rowOf(main), measurement);
  const shrunkDeep = tailflex.shrinkOfRow(rowOf(deep), measurement);
  assert.ok(Math.abs(shrunkMain.edgePct - 2.0) < 0.1, `main ${shrunkMain.edgePct}`);
  assert.ok(Math.abs(shrunkDeep.edgePct - 0.4) < 0.1, `deep ${shrunkDeep.edgePct}`);
  assert.deepEqual([shrunkMain.tau, shrunkMain.c], [0.0185, 0.094]);
  const sized = (line) => {
    const row = rowOf(line);
    const shrunk = tailflex.shrinkOfRow(row, measurement);
    return { ...row, shrunkEdgePct: shrunk.edgePct, stake: kelly.kellyStakeFromEdge({ bookPrice: row.price, edgePct: shrunk.edgePct, bankroll: 35000, multiplier: 0.25 }).stake };
  };
  const mainRow = sized(main);
  const deepRow = sized(deep);
  assert.ok(Math.abs(mainRow.stake - 170) < 3, `main stake ${mainRow.stake}`);
  assert.ok(Math.abs(deepRow.stake - 11) < 2, `deep stake ${deepRow.stake}`);
  const rankMain = tailflex.rankScoreOf(mainRow);
  const rankDeep = tailflex.rankScoreOf(deepRow);
  assert.ok(Math.abs(rankMain - (mainRow.shrunkEdgePct / 100) * mainRow.stake) < 1e-9);
  assert.ok(rankMain > 20 * rankDeep, `${rankMain} vs ${rankDeep}`);
  // At every tau from 1% to 3% the main still wins the card.
  for (const tau of [0.01, 0.03]) {
    const at = { ...measurement, markets: { "2:1:2": { ...measurement.markets["2:1:2"], tau } } };
    const score = (line) => { const row = rowOf(line); const shrunk = tailflex.shrinkOfRow(row, at); return tailflex.rankScoreOf({ ...row, shrunkEdgePct: shrunk.edgePct, stake: kelly.kellyStakeFromEdge({ bookPrice: row.price, edgePct: shrunk.edgePct, bankroll: 35000, multiplier: 0.25 }).stake }); };
    assert.ok(score(main) > score(deep), `tau ${tau}`);
  }
  // No stake, no positive edge, no American price: nothing.
  assert.equal(tailflex.rankScoreOf({ ...mainRow, stake: null }), null);
  assert.equal(tailflex.rankScoreOf({ ...mainRow, shrunkEdgePct: null }), null);
  assert.equal(tailflex.shrinkOfRow({ ...mainRow, edgePct: 0 }, measurement), null);
  assert.equal(tailflex.shrinkOfRow({ ...mainRow, price: 50 }, measurement), null);
});

test("groupEdges by the rank score picks -14.5 over the bigger-edge -27.5 and the -257 short", () => {
  const rows = [MINUS_27_5, MINUS_2_5, MINUS_14_5].map((row) => ranked(row, MEASUREMENT_C10));
  const [card] = feed.groupEdges(rows, (row) => row.rankScore);
  // The +1300 longshot at 17.3% now ranks last: 1% at the -257 short keeps more of its edge than 17.3% four SD out.
  assert.deepEqual(card.rows.map((row) => row.key), ["novig-14.5", "novig-2.5", "hardrock-27.5"]);
  assert.equal(card.best.key, "novig-14.5");
  // By raw edge the -27.5 would have won the card.
  assert.equal(feed.groupEdges(rows)[0].best.key, "hardrock-27.5");
});

test("the video's check: -12 at 4% (+156) vs -24 at 10% (+650) off a -8.5 main — the near rung wins at every c (the deep one only kept the card under the old keep)", () => {
  const measurementAt = (c) => ({ markets: { [NFL_SPREAD_FG]: { c, measured: true, tau: 0.0185, tauMeasured: true } }, mainFairs: new Map([["777:1:2:1", { points: -8.5, fairProb: 0.5 }]]) });
  const twelve = seahawksRow({ key: "-12", points: -12, price: 156, edgePct: 4, book: 89 });
  const twentyFour = seahawksRow({ key: "-24", points: -24, price: 650, edgePct: 10, book: 89 });
  const winner = (c) => feed.groupEdges([twelve, twentyFour].map((row) => ranked(row, measurementAt(c))), (row) => row.rankScore)[0].best.key;
  assert.equal(winner(0.04), "-12");
  assert.equal(winner(0.08), "-12");
  assert.equal(winner(0.10), "-12");
});

test("a moneyline shrinks on the fuselage alone; a side with no Unabated main line measures from even money", () => {
  const moneyline = { ...MINUS_14_5, betTypeId: 1, points: null };
  const decimal = kelly.americanToDecimal(270);
  const shrinkMl = tailflex.shrinkOfRow(moneyline, MEASUREMENT_C10).shrink;
  assert.ok(Math.abs(shrinkMl - tailflex.shrinkFactor({ decimal, fairProb: 1.139 / decimal, zMain: null, c: 0.1, tau: tailflex.TAU_FALLBACK }).shrink) < 1e-12);
  const orphan = tailflex.shrinkOfRow({ ...MINUS_14_5, eventId: 999 }, MEASUREMENT_C10).shrink;
  assert.ok(Math.abs(orphan - tailflex.shrinkFactor({ decimal, fairProb: 1.139 / decimal, zMain: 0, c: 0.1, tau: tailflex.TAU_FALLBACK }).shrink) < 1e-12);
  // No measurement at all: the fallbacks for both c and tau.
  assert.deepEqual([tailflex.cOf(null, MINUS_14_5), tailflex.tauOf(null, MINUS_14_5), tailflex.isTauMeasured(null, MINUS_14_5)], [0.10, 0.0185, false]);
});

// ---- measuring tau ------------------------------------------------------------

// Two exchanges on every rung of syntheticState: book 89 as built, book 90 a
// copy whose devigged fair sits `secondGap` probit from Unabated's. With
// `realEdge(eventId, points)` (probit) both exchanges and Unabated share
// that much of their difference from the ask.
function twoExchangeState({ events, realEdgeZ = () => 0, noiseU = () => 0, noiseX = () => 0 }) {
  const state = syntheticState({ events, cTrue: 0 });
  state.books[90] = { id: 90, hasLiquidity: true };
  const added = [];
  for (const [key, line] of Object.entries(state.lines)) {
    if (line.bookId !== EXCHANGE_BOOK) continue;
    // The exchange's ask: the fair plus the real edge the book missed, plus its own noise, under 2% vig.
    const z = tailflex.qnorm(kelly.americanToProb(line.bacr));
    const sign = line.sideIndex === 0 ? 1 : -1;
    const trueZ = sign * (Math.abs(z) * sign) ;
    void trueZ;
    const askZ = (book) => z - sign * realEdgeZ(line.eventId, line.points) + sign * noiseX(line.eventId, line.points, book) + 0.02;
    state.lines[key] = { ...line, sourcePrice: tailflex.pnorm(askZ(EXCHANGE_BOOK)) };
    added.push([key.replace(`ms${EXCHANGE_BOOK}`, "ms90"), { ...line, bookId: 90, sourcePrice: tailflex.pnorm(askZ(90)) }]);
    // Unabated's fair off the truth by its own noise (the bacr every book's line carries).
    const unabatedZ = z + sign * noiseU(line.eventId, line.points);
    const bacr = americanOf(tailflex.pnorm(unabatedZ));
    state.lines[key].bacr = bacr;
    added[added.length - 1][1].bacr = bacr;
  }
  for (const [key, line] of added) state.lines[key] = line;
  return state;
}

test("tauSamplesOfRungs: a rung needs two exchanges; each exchange is priced against Unabated and against the OTHER exchanges, both sides", () => {
  const one = tailflex.measureTailFlex(syntheticState({ events: 3, cTrue: 0.06 }), { now: NOW });
  assert.equal(one.tauSamples.length, 0);
  const two = tailflex.measureTailFlex(twoExchangeState({ events: 3 }), { now: NOW });
  // 3 events x 6 rungs x 2 exchanges x 2 sides.
  assert.equal(two.tauSamples.length, 72);
  const sample = two.tauSamples.find((item) => item.key === `${EXCHANGE_BOOK}|1:1:2:0:-7|0`);
  assert.ok(sample && sample.marketKey === NFL_SPREAD_FG);
  // Identical exchanges at Unabated's fair: the other exchange's read is exactly zero (priced at the
  // mid, not the ask, so the 2% vig never shows), Unabated's within its rounding to a whole American price.
  assert.ok(two.tauSamples.every((item) => Math.abs(item.edgeExchanges) < 1e-5));
  assert.ok(Math.abs(sample.edgeUnabated) < 0.01, `${sample.edgeUnabated}`);
  // Fewer than 100 samples: the fallback, the count reported.
  assert.deepEqual(tailflex.tauOfSamples(two.tauSamples.filter((item) => item.key.endsWith("|0"))), { tau: tailflex.TAU_FALLBACK, measured: false, rungCount: 36, covariance: null });
});

test("tauOfSamples recovers the spread of real edges when Unabated and the exchanges err independently; books at the truth give no tau", () => {
  // Deterministic "noise": a different phase per source, so the two errors are uncorrelated over the rungs.
  const wave = (seed) => (eventId, points) => 0.06 * Math.sin(seed + eventId * 1.7 + points * 0.9);
  const realEdgeZ = (eventId, points) => 0.05 * Math.cos(eventId * 2.3 + points * 1.1);
  const state = twoExchangeState({ events: 40, realEdgeZ, noiseU: wave(1), noiseX: (eventId, points, book) => wave(book)(eventId, points) });
  const measured = tailflex.tauOfSamples(tailflex.measureTailFlex(state, { now: NOW }).tauSamples);
  assert.equal(measured.measured, true);
  // The real edge in EV units: ~0.05 probit x dnorm(z) x decimal ~ 2-4% across the rungs.
  assert.ok(measured.tau > 0.015 && measured.tau < 0.045, `tau = ${measured.tau}`);
  // 960 built; a rung whose noisy Unabated fair rounds to its neighbour's whole American price is a flat tail and never measures (feed.flatFairRungKeys).
  assert.ok(measured.rungCount >= 700, `rungCount = ${measured.rungCount}`);
  // An exchange's own mispricing IS a real edge (both reads see the cheap
  // ask), so "no real edge" is every book at the truth with only Unabated
  // noisy: the two reads then share nothing, the covariance sits at zero
  // (Unabated's rounding to a whole American price apart) and tau falls back.
  const flat = tailflex.tauOfSamples(tailflex.measureTailFlex(twoExchangeState({ events: 40, noiseU: wave(1) }), { now: NOW }).tauSamples);
  assert.ok(flat.covariance !== null && Math.abs(flat.covariance) < 0.005 ** 2, `covariance = ${flat.covariance}`);
  assert.ok(flat.tau < 0.005 || !flat.measured, `tau = ${flat.tau} measured ${flat.measured}`);
});

test("tauOfSamples measures from 101 samples: the 1% trim comes off first, and 100 then leaves 99", () => {
  const rows = (count) => Array.from({ length: count }, (_, index) => { const edge = -0.02 + (0.04 * index) / (count - 1); return { marketKey: NFL_SPREAD_FG, edgeUnabated: edge, edgeExchanges: edge }; });
  assert.deepEqual([tailflex.tauOfSamples(rows(100)).measured, tailflex.tauOfSamples(rows(100)).rungCount], [false, 99]);
  assert.deepEqual([tailflex.tauOfSamples(rows(101)).measured, tailflex.tauOfSamples(rows(101)).rungCount], [true, 100]);
  // Two reads that disagree on where the edges are (covariance <= 0): the fallback, with the covariance reported.
  const crossed = rows(150).map((row) => ({ ...row, edgeExchanges: -row.edgeUnabated }));
  assert.deepEqual([tailflex.tauOfSamples(crossed).measured, tailflex.tauOfSamples(crossed).tau, tailflex.tauOfSamples(crossed).covariance < 0], [false, tailflex.TAU_FALLBACK, true]);
});

test("tauOfSamples drops the top 1% by either read, so one stale quote cannot own the covariance", () => {
  const samples = [];
  for (let i = 0; i < 200; i += 1) samples.push({ key: `k${i}`, marketKey: "m", edgeUnabated: -0.05 + 0.02 * Math.sin(i), edgeExchanges: -0.05 + 0.02 * Math.sin(i) });
  const clean = tailflex.tauOfSamples(samples);
  const spiked = tailflex.tauOfSamples([...samples.slice(0, 198), { key: "stale", marketKey: "m", edgeUnabated: 0.5, edgeExchanges: 0.5 }, { key: "stale2", marketKey: "m", edgeUnabated: -0.6, edgeExchanges: -0.6 }]);
  assert.ok(Math.abs(clean.tau - Math.sqrt(0.02 ** 2 / 2)) < 0.001, `clean tau = ${clean.tau}`);
  assert.ok(Math.abs(spiked.tau - clean.tau) < 0.001, `spiked tau = ${spiked.tau}`);
  assert.equal(spiked.rungCount, 198);
});

test("createTauPool: the latest reading per key, evicted after the window, persisted when changed or hourly; withPooledTau replaces tau and keeps c", () => {
  const pool = tailflex.createTauPool({ windowMs: 1000, heartbeatMs: 100 });
  const sample = (key, edgeUnabated) => ({ key, marketKey: NFL_SPREAD_FG, edgeUnabated, edgeExchanges: edgeUnabated - 0.01 });
  assert.deepEqual(pool.absorb([sample("a", -0.05), sample("b", -0.04)], 0).map((row) => row.key), ["a", "b"]);
  // Unchanged within the heartbeat: nothing to persist; changed: persisted; past the heartbeat: persisted again.
  assert.deepEqual(pool.absorb([sample("a", -0.05), sample("b", -0.03)], 50).map((row) => [row.key, row.seenAt]), [["b", 50]]);
  assert.deepEqual(pool.absorb([sample("a", -0.05)], 150).map((row) => [row.key, row.seenAt]), [["a", 150]]);
  assert.equal(pool.size(), 2);
  // "b" was last seen at 50: gone after the window; "a" stays.
  pool.absorb([sample("a", -0.05)], 1100);
  assert.equal(pool.size(), 1);
  // A persisted copy restores, and the pool measures what is inside the window.
  pool.load([{ key: "c", marketKey: NFL_SPREAD_FG, edgeUnabated: -0.02, edgeExchanges: -0.03, seenAt: 1100 }]);
  assert.equal(pool.size(), 2);
  // A persisted row older than the reading held (the service came back after the pool filled) does not replace it; a newer one does.
  pool.load([{ ...sample("a", -0.08), seenAt: 100 }]);
  assert.deepEqual(pool.absorb([sample("a", -0.05)], 1110), []);
  pool.load([{ ...sample("a", -0.08), seenAt: 1150 }]);
  assert.deepEqual(pool.absorb([sample("a", -0.08)], 1160), []);
  assert.deepEqual(pool.absorb([sample("a", -0.05)], 1170).map((row) => row.key), ["a"]);
  const measured = pool.measure(1100);
  assert.deepEqual(measured, { [NFL_SPREAD_FG]: { tau: tailflex.TAU_FALLBACK, measured: false, rungCount: 2, covariance: null } });
  const base = tailflex.measureTailFlex(syntheticState({ events: 25, cTrue: 0.062, baseline: 0.02 }), { now: NOW });
  const pooled = tailflex.withPooledTau(base, { [NFL_SPREAD_FG]: { tau: 0.021, measured: true, rungCount: 400, covariance: 0.021 ** 2 }, "9:9:9": { tau: 0.015, measured: true, rungCount: 150, covariance: 0.015 ** 2 } });
  assert.equal(pooled.markets[NFL_SPREAD_FG].c, base.markets[NFL_SPREAD_FG].c);
  assert.deepEqual([pooled.markets[NFL_SPREAD_FG].tau, pooled.markets[NFL_SPREAD_FG].tauMeasured, pooled.markets[NFL_SPREAD_FG].tauRungCount], [0.021, true, 400]);
  // A market the pool knows and this refresh does not: the pool's tau, the fallback c.
  assert.deepEqual([pooled.markets["9:9:9"].c, pooled.markets["9:9:9"].measured, pooled.markets["9:9:9"].tau], [0.10, false, 0.015]);
  // A market this refresh has and the pool does not: the fallback tau, never the one-refresh read.
  const unpooled = tailflex.withPooledTau(base, {});
  assert.deepEqual([unpooled.markets[NFL_SPREAD_FG].tau, unpooled.markets[NFL_SPREAD_FG].tauMeasured], [tailflex.TAU_FALLBACK, false]);
  // measureTailFlex itself carries no tau: one refresh is never sized on.
  assert.equal(base.markets[NFL_SPREAD_FG].tau, undefined);
  assert.equal(tailflex.tauOf(base, { leagueId: NFL, periodTypeId: FULL_GAME, betTypeId: SPREAD }), tailflex.TAU_FALLBACK);
  assert.equal(base.tauSamples.length, 0);
});
