// Run: node --test unabated_ticket/tests
// tailflex.js: the probit devig, the live measurement of c off two-sided
// exchange rungs, and the rank score (EV dollars after flex) that picks a
// card's best line. The worked card is NFL Chargers @ Seahawks (2026-09-30
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

test("keepFactor on the worked card: -14.5 +270 at 13.9% (fair 30.8%) keeps ~0.81, ~11.2% after flex", () => {
  const decimal = kelly.americanToDecimal(270);
  const fairProb = 1.139 / decimal;
  assert.ok(Math.abs(fairProb - 0.308) < 0.001);
  const keep = tailflex.keepFactor({ edge: 0.139, decimal, fairProb, zMain: tailflex.qnorm(0.481), c: 0.10 });
  assert.ok(Math.abs(keep - 0.809) < 0.002, `keep = ${keep}`);
  assert.ok(Math.abs(keep * 0.139 - 0.1125) < 0.0005);
  // A wider c flexes more; at the main number only the one-cent fuselage is left (c drops out).
  assert.ok(tailflex.keepFactor({ edge: 0.139, decimal, fairProb, zMain: tailflex.qnorm(0.481), c: 0.2 }) < keep);
  const atMain = tailflex.keepFactor({ edge: 0.139, decimal, fairProb, zMain: tailflex.qnorm(fairProb), c: 0.10 });
  assert.equal(atMain, tailflex.keepFactor({ edge: 0.139, decimal, fairProb, zMain: tailflex.qnorm(fairProb), c: 0 }));
  assert.ok(atMain > 0.94 && atMain < 0.96, `atMain = ${atMain}`);
});

test("rankOfRow = keep x edge x stake, the stake being the row's own Kelly stake on Unabated's edge", () => {
  const near = tailflex.rankOfRow(MINUS_14_5, MEASUREMENT_C10);
  const deep = tailflex.rankOfRow(MINUS_27_5, MEASUREMENT_C10);
  const short = tailflex.rankOfRow(MINUS_2_5, MEASUREMENT_C10);
  assert.ok(Math.abs(near.score - near.keep * 0.139 * MINUS_14_5.stake) < 1e-9);
  assert.ok(Math.abs(near.score - 14.47) < 0.02, `near = ${near.score}`);
  assert.ok(Math.abs(deep.keep - 0.26) < 0.005, `deep keep = ${deep.keep}`);
  assert.ok(Math.abs(deep.score - 1.49) < 0.02, `deep = ${deep.score}`);
  assert.ok(short.score < 0.1, `short = ${short.score}`);
  // No stake or no positive edge: no score.
  assert.equal(tailflex.rankOfRow({ ...MINUS_14_5, stake: null }, MEASUREMENT_C10), null);
  assert.equal(tailflex.rankOfRow({ ...MINUS_14_5, edgePct: 0 }, MEASUREMENT_C10), null);
});

test("groupEdges by the rank score picks -14.5 over the bigger-edge -27.5 and the -257 short; stakes untouched", () => {
  const rows = [MINUS_27_5, MINUS_2_5, MINUS_14_5].map((row) => ({ ...row, rankScore: tailflex.rankOfRow(row, MEASUREMENT_C10).score }));
  const [card] = feed.groupEdges(rows, (row) => row.rankScore);
  assert.deepEqual(card.rows.map((row) => row.key), ["novig-14.5", "hardrock-27.5", "novig-2.5"]);
  assert.equal(card.best.key, "novig-14.5");
  // By raw stake the -14.5 also leads, but by raw edge the -27.5 would have won the card.
  assert.equal(feed.groupEdges(rows)[0].best.key, "hardrock-27.5");
  for (const row of card.rows) {
    const raw = kelly.kellyStakeFromEdge({ bookPrice: row.price, edgePct: row.edgePct, bankroll: BANKROLL, multiplier: QUARTER_KELLY }).stake;
    assert.equal(row.stake, raw);
  }
});

test("the video's check: -12 at 4% (+156) vs -24 at 10% (+650) off a -8.5 main — the deep rung wins only below c ~6%", () => {
  const measurementAt = (c) => ({ markets: { [NFL_SPREAD_FG]: { c, measured: true } }, mainFairs: new Map([["777:1:2:1", { points: -8.5, fairProb: 0.5 }]]) });
  const twelve = seahawksRow({ key: "-12", points: -12, price: 156, edgePct: 4, book: 89 });
  const twentyFour = seahawksRow({ key: "-24", points: -24, price: 650, edgePct: 10, book: 89 });
  const winner = (c) => feed.groupEdges([twelve, twentyFour].map((row) => ({ ...row, rankScore: tailflex.rankOfRow(row, measurementAt(c)).score })), (row) => row.rankScore)[0].best.key;
  assert.equal(winner(0.04), "-24");
  assert.equal(winner(0.08), "-12");
  assert.equal(winner(0.10), "-12");
});

test("a moneyline ranks on the fuselage alone; a side with no Unabated main line measures from even money", () => {
  const moneyline = { ...MINUS_14_5, betTypeId: 1, points: null };
  const keepMl = tailflex.rankOfRow(moneyline, MEASUREMENT_C10).keep;
  const decimal = kelly.americanToDecimal(270);
  assert.ok(Math.abs(keepMl - tailflex.keepFactor({ edge: 0.139, decimal, fairProb: 1.139 / decimal, zMain: null, c: 0.1 })) < 1e-12);
  const orphan = tailflex.rankOfRow({ ...MINUS_14_5, eventId: 999 }, MEASUREMENT_C10).keep;
  assert.ok(Math.abs(orphan - tailflex.keepFactor({ edge: 0.139, decimal, fairProb: 1.139 / decimal, zMain: 0, c: 0.1 })) < 1e-12);
});
