// server/taker.js: the automated taker's per-row decision (T1, places nothing).
// Board: fixtures/v2_venue_ids_slice.json (real CFB event 123742, Duquesne @
// Washington State, with Kalshi rung ids), edges set on its Kalshi lines here.
"use strict";

const test = require("node:test");
const assert = require("node:assert/strict");
const fs = require("node:fs");
const path = require("node:path");
const feed = require("../extension/feed.js");
const pricecheck = require("../extension/pricecheck.js");
const taker = require("../server/taker.js");
const { moveLine } = require("./helpers.js");

const CFB = 2;
const EVENT_START_MS = 1789846200000;
const NOW = EVENT_START_MS - 24 * 3600 * 1000;
// Kalshi WSU -35.5 alt: -143, all-in 0.587157 (a 57c ask plus the fee).
const WSU_ALT_KEY = "420942308:ms105:si1:tid717:alt-35.5";
const DUQ_ALT_KEY = "420942310:ms105:si0:tid927:alt35.5";
const ML_WSU_KEY = "411181220:ms105:si1:tid717";
const EDGE_SETTINGS = { leagues: [CFB], periods: [1] };

function boardWithEdge(key, ge, fields = {}) {
  const json = JSON.parse(fs.readFileSync(path.join(__dirname, "fixtures", "v2_venue_ids_slice.json"), "utf8"));
  const state = feed.parseSnapshot(json, { leagueId: CFB });
  moveLine(state, key, { ge, bacr: -200, sequenceNumber: NOW - 60 * 1000, modifiedOn: new Date(NOW - 60 * 1000).toISOString(), ...fields });
  return state;
}

function rowOf(state, key) {
  const row = taker.candidateRows(state, EDGE_SETTINGS, NOW).find((candidate) => candidate.key === key);
  assert.ok(row, `expected ${key} among the candidates`);
  return row;
}

// A book whose best opposite bid makes `ask` on the given side.
function bookWithAsk(side, ask, depth = 1000) {
  const level = [[(1 - ask).toFixed(4), String(depth)]];
  return taker.parseKalshiOrderbook({ orderbook_fp: side === "yes" ? { no_dollars: level, yes_dollars: [] } : { yes_dollars: level, no_dollars: [] } });
}

function decide(state, key, quote, extra = {}) {
  const row = rowOf(state, key);
  return taker.decideRow({
    row, contract: taker.resolveKalshiContract(row, state), quote,
    priceIndex: pricecheck.buildPriceIndex(state, { now: NOW, maxLineAgeMs: taker.MAX_LINE_AGE_MS }),
    sizing: { matches: [], ladderOf: () => null, teasers: [] }, ...extra,
  });
}

test("all-in price adds the 7% taker fee unrounded, as Unabated shows Kalshi", () => {
  assert.equal(Math.round(taker.kalshiAllInPrice(0.22) * 1e5) / 1e5, 0.23201);
  assert.equal(Math.round(taker.kalshiAllInPrice(0.57) * 1e6) / 1e6, 0.587157);
  assert.equal(taker.probToAmerican(0.25), 300);
  assert.equal(taker.probToAmerican(0.6), -150);
});

test("a spread alt resolves to its rung's contract: Y buys YES, N buys NO", () => {
  const state = boardWithEdge(WSU_ALT_KEY, 0.05);
  assert.deepEqual(taker.resolveKalshiContract(rowOf(state, WSU_ALT_KEY), state),
    { ticker: "KXNCAAFSPREAD-26SEP19DUQWSU-WSU36", side: "yes", sourceKey: "Y-KXNCAAFSPREAD-26SEP19DUQWSU-WSU36", derived: false });
  const other = boardWithEdge(DUQ_ALT_KEY, 0.05);
  assert.deepEqual(taker.resolveKalshiContract(rowOf(other, DUQ_ALT_KEY), other),
    { ticker: "KXNCAAFSPREAD-26SEP19DUQWSU-WSU36", side: "no", sourceKey: "N-KXNCAAFSPREAD-26SEP19DUQWSU-WSU36", derived: false });
});

test("a line no rung names, or whose rung disagrees on the number, has no ticker", () => {
  const state = boardWithEdge(WSU_ALT_KEY, 0.05);
  const row = rowOf(state, WSU_ALT_KEY);
  assert.equal(taker.contractFromRungs({ ...row, key: "nothing" }).reason, "no_ticker");
  assert.match(taker.contractFromRungs({ ...row, points: -34.5 }).detail, /row is -34.5/);
});

test("a full-game moneyline's ticker is derived from the spread rungs' team codes", () => {
  const state = boardWithEdge(ML_WSU_KEY, 0.05, { price: -400 });
  const row = rowOf(state, ML_WSU_KEY);
  assert.deepEqual(taker.resolveKalshiContract(row, state), { ticker: "KXNCAAFGAME-26SEP19DUQWSU-WSU", side: "yes", sourceKey: null, derived: true });
  // Only WSU's code is on the spread rungs (Kalshi lists WSU-covers markets): the Duquesne side has none.
  assert.match(taker.derivedMoneylineContract({ ...row, sideIndex: 0 }, state).detail, /side 0/);
  assert.equal(taker.derivedMoneylineContract({ ...row, periodTypeId: 2 }, state).reason, "no_ticker");
});

test("orderbook: the ask is 1 - the best opposite bid; legacy cents parse too", () => {
  const book = taker.parseKalshiOrderbook({ orderbook_fp: { yes_dollars: [["0.4000", "50"], ["0.4200", "7"]], no_dollars: [["0.5500", "30"]] } });
  assert.deepEqual(taker.quoteOf(book, "yes"), { ask: 0.45, depth: 30 });
  assert.deepEqual(taker.quoteOf(book, "no"), { ask: 0.58, depth: 7 });
  const legacy = taker.parseKalshiOrderbook({ orderbook: { yes: [[40, 10]], no: null } });
  assert.equal(taker.quoteOf(legacy, "yes"), null);
  assert.deepEqual(taker.quoteOf(legacy, "no"), { ask: 0.6, depth: 10 });
});

test("at Unabated's own price, best at the number, a 5% edge bets quarter Kelly of $500 in contracts", () => {
  // Prophet Exchange and Novig sit at the same number at worse or no price; none is better than -143.
  const state = boardWithEdge(WSU_ALT_KEY, 0.05);
  for (const key of Object.keys(state.lines)) {
    if (!key.includes("ms105") && key.includes("si1:tid717:alt-35.5")) moveLine(state, key, { price: -200 });
  }
  const decision = decide(state, WSU_ALT_KEY, taker.quoteOf(bookWithAsk("yes", 0.57), "yes"));
  assert.equal(decision.reason, "bet", decision.detail);
  assert.equal(decision.liveAsk, 0.57);
  // At Unabated's own all-in the live edge is Unabated's edge (the fair is read off its cheaper price).
  assert.ok(Math.abs(decision.liveEdge - 0.05) < 0.0005, String(decision.liveEdge));
  assert.ok(decision.stake > 0 && decision.stake <= taker.BOT_BANKROLL * taker.BOT_KELLY);
  assert.equal(decision.contracts, Math.floor(decision.stake / decision.liveAllIn));
});

test("depth cuts the contracts; none resting at the ask is no_book", () => {
  const state = boardWithEdge(WSU_ALT_KEY, 0.05);
  for (const key of Object.keys(state.lines)) {
    if (!key.includes("ms105") && key.includes("si1:tid717:alt-35.5")) moveLine(state, key, { price: -200 });
  }
  assert.equal(decide(state, WSU_ALT_KEY, taker.quoteOf(bookWithAsk("yes", 0.57, 3), "yes")).contracts, 3);
  assert.equal(decide(state, WSU_ALT_KEY, null).reason, "no_book");
  assert.equal(decide(state, WSU_ALT_KEY, null, { quoteError: "HTTP 404" }).detail, "KXNCAAFSPREAD-26SEP19DUQWSU-WSU36: HTTP 404");
});

test("price rule: 0.1c worse than Unabated still trades, a cent worse is price_moved", () => {
  const state = boardWithEdge(WSU_ALT_KEY, 0.05);
  for (const key of Object.keys(state.lines)) {
    if (!key.includes("ms105") && key.includes("si1:tid717:alt-35.5")) moveLine(state, key, { price: -200 });
  }
  // 57.08c ask: all-in 0.58799, 0.08c worse than Unabated's 0.587157.
  const quote = { ask: 0.5708, depth: 1000 };
  assert.equal(decide(state, WSU_ALT_KEY, quote).reason, "bet");
  assert.equal(decide(state, WSU_ALT_KEY, { ask: 0.58, depth: 1000 }).reason, "price_moved");
});

test("edge floor 3%, cap 20%, at the live price", () => {
  const short = boardWithEdge(WSU_ALT_KEY, 0.031);
  // 57.08c: inside the 0.1c tolerance, but 3.1% shrinks to 2.95% there.
  assert.equal(decide(short, WSU_ALT_KEY, { ask: 0.5708, depth: 1000 }).reason, "edge_short");
  const big = boardWithEdge(WSU_ALT_KEY, 0.25);
  for (const key of Object.keys(big.lines)) {
    if (!key.includes("ms105") && key.includes("si1:tid717:alt-35.5")) moveLine(big, key, { price: -200 });
  }
  assert.equal(decide(big, WSU_ALT_KEY, { ask: 0.57, depth: 1000 }).reason, "cap");
});

test("one live book better at the number is not_best", () => {
  const state = boardWithEdge(WSU_ALT_KEY, 0.05);
  const novigKey = Object.keys(state.lines).find((key) => key.includes("ms89") && key.endsWith("si1:tid717:alt-35.5")) || "420942308:ms89:si1:tid717";
  moveLine(state, novigKey, { points: -35.5, price: -120, statusId: 1, isBlurred: false, modifiedOn: new Date(NOW - 1000).toISOString(), sequenceNumber: NOW - 1000 });
  const decision = decide(state, WSU_ALT_KEY, { ask: 0.57, depth: 1000 });
  assert.equal(decision.reason, "not_best");
  assert.match(decision.detail, /Novig -120/);
});

test("the decision signature changes with the ask, not with the clock", () => {
  const base = { reason: "bet", ticker: "T", liveAsk: 0.57, unabatedPrice: -143, contracts: 3, liveEdge: 0.05 };
  assert.equal(taker.decisionSignature(base), taker.decisionSignature({ ...base, liveEdge: 0.051 }));
  assert.notEqual(taker.decisionSignature(base), taker.decisionSignature({ ...base, liveAsk: 0.58 }));
});
