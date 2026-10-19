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
const betsLib = require("../extension/bets.js");
const edgeRows = require("../extension/edgerows.js");
const teams = require("../extension/teams.js");
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

// Main lines in the fixture: Unabated's ladders skip the main-number rung (52.5 / -38.5), as live.
const UNDER_MAIN_KEY = "433644483:ms105:si1:tid717";
const OVER_MAIN_KEY = "433644484:ms105:si0:tid927";
const WSU_MAIN_KEY = "420942308:ms105:si1:tid717";
const DUQ_MAIN_KEY = "420942310:ms105:si0:tid927";

test("a main number's ticker is derived from its own ladder's rungs: Kalshi names X.5 as X+1", () => {
  const derived = (key) => {
    const state = boardWithEdge(key, 0.05);
    return taker.resolveKalshiContract(rowOf(state, key), state);
  };
  assert.deepEqual(derived(UNDER_MAIN_KEY), { ticker: "KXNCAAFTOTAL-26SEP19DUQWSU-53", side: "no", sourceKey: null, derived: true });
  assert.deepEqual(derived(OVER_MAIN_KEY), { ticker: "KXNCAAFTOTAL-26SEP19DUQWSU-53", side: "yes", sourceKey: null, derived: true });
  assert.deepEqual(derived(WSU_MAIN_KEY), { ticker: "KXNCAAFSPREAD-26SEP19DUQWSU-WSU39", side: "yes", sourceKey: null, derived: true });
  assert.deepEqual(derived(DUQ_MAIN_KEY), { ticker: "KXNCAAFSPREAD-26SEP19DUQWSU-WSU39", side: "no", sourceKey: null, derived: true });
});

test("a main number refuses a ticker when its ladder cannot vouch for one", () => {
  const state = boardWithEdge(UNDER_MAIN_KEY, 0.05);
  const row = rowOf(state, UNDER_MAIN_KEY);
  // No rungs hang off the line.
  assert.match(taker.contractFromRungs({ ...row, key: "nothing" }).detail, /has no rungs/);
  // A whole number has no X.5 strike.
  assert.match(taker.contractFromRungs({ ...row, points: 52 }).detail, /not a half point/);
  // An alt the map misses is never derived.
  assert.equal(taker.contractFromRungs({ ...row, isAlt: true }).detail, "no Kalshi rung names this line");
  // The Over side's ladder carrying NO rungs is a contradiction, not a guess.
  assert.match(taker.contractFromRungs({ ...row, sideIndex: 0 }).detail, /carries N rungs/);
  // Rungs that disagree on the event.
  const contracts = { ...row.venueIds.kalshiContracts, "N-KXNCAAFTOTAL-26SEP20OTHER-50": { lineKey: null, mainKey: row.key, points: 49.5, sideIndex: 1 } };
  assert.match(taker.contractFromRungs({ ...row, venueIds: { ...row.venueIds, kalshiContracts: contracts } }).detail, /disagree/);
  // A spread main with no rung on its side of zero.
  const spreadState = boardWithEdge(WSU_MAIN_KEY, 0.05);
  const spreadRow = rowOf(spreadState, WSU_MAIN_KEY);
  assert.match(taker.contractFromRungs({ ...spreadRow, points: 1.5 }).detail, /no rungs on its side/);
  // Same-side rungs that contradict the sign (an N rung at a negative number).
  const flipped = Object.fromEntries(Object.entries(spreadRow.venueIds.kalshiContracts).filter(([, target]) => target.mainKey === spreadRow.key).map(([key, target]) => [key.replace(/^Y-/, "N-"), target]));
  assert.match(taker.contractFromRungs({ ...spreadRow, venueIds: { ...spreadRow.venueIds, kalshiContracts: flipped } }).detail, /reads N rungs/);
});

test("a close spread whose ladder crosses zero derives from the rungs on its own side only", () => {
  const state = boardWithEdge(WSU_MAIN_KEY, 0.05);
  const row = rowOf(state, WSU_MAIN_KEY);
  // WSU -1.5 main: -2.5 / -3.5 are YES on WSU, +0.5 / +1.5 are NO on Duquesne-covers markets.
  const kalshiContracts = {
    "Y-KXNCAAFSPREAD-26SEP19DUQWSU-WSU3": { lineKey: "a", mainKey: row.key, points: -2.5, sideIndex: 1 },
    "Y-KXNCAAFSPREAD-26SEP19DUQWSU-WSU4": { lineKey: "b", mainKey: row.key, points: -3.5, sideIndex: 1 },
    "N-KXNCAAFSPREAD-26SEP19DUQWSU-DUQ1": { lineKey: "c", mainKey: row.key, points: 0.5, sideIndex: 1 },
    "N-KXNCAAFSPREAD-26SEP19DUQWSU-DUQ2": { lineKey: "d", mainKey: row.key, points: 1.5, sideIndex: 1 },
  };
  assert.deepEqual(taker.contractFromRungs({ ...row, points: -1.5, venueIds: { ...row.venueIds, kalshiContracts } }),
    { ticker: "KXNCAAFSPREAD-26SEP19DUQWSU-WSU2", side: "yes", sourceKey: null, derived: true });
});

test("a derived main-number ticker must still pass Kalshi's market check", () => {
  const state = boardWithEdge(UNDER_MAIN_KEY, 0.05);
  const row = rowOf(state, UNDER_MAIN_KEY);
  const contract = taker.resolveKalshiContract(row, state);
  const market = { status: "active", ticker: contract.ticker, floor_strike: 52.5 };
  assert.equal(taker.marketMismatch(market, row, contract), null);
  assert.match(taker.marketMismatch({ ...market, floor_strike: 53.5 }, row, contract), /floor_strike/);
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

test("at Unabated's own price, best at the number, a 5% edge bets quarter Kelly of $2,000 in contracts", () => {
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
  // Quarter Kelly of $2,000 at a 5% edge on -143: $35.56, bought as 60 contracts at the 0.5872 all-in.
  assert.equal(Math.round(decision.stake * 100) / 100, 35.56);
  assert.deepEqual({ contracts: decision.contracts, cost: decision.cost, liveAllIn: decision.liveAllIn }, { contracts: 60, cost: 35.2294, liveAllIn: 0.5872 });
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

test("price rule: up to 1c better than Unabated trades (with the bigger edge); more is a stale fair, price_moved", () => {
  const state = boardWithEdge(WSU_ALT_KEY, 0.05);
  for (const key of Object.keys(state.lines)) {
    if (!key.includes("ms105") && key.includes("si1:tid717:alt-35.5")) moveLine(state, key, { price: -200 });
  }
  // 56.1c ask: all-in 0.578288, 0.89c better than Unabated's 0.587157.
  const better = decide(state, WSU_ALT_KEY, { ask: 0.561, depth: 1000 });
  assert.equal(better.reason, "bet", better.detail);
  assert.ok(better.liveEdge > 0.05, String(better.liveEdge));
  // 55c ask: all-in 0.567325, 1.98c better: the market ran away from Unabated's snapshot.
  const ranAway = decide(state, WSU_ALT_KEY, { ask: 0.55, depth: 1000 });
  assert.equal(ranAway.reason, "price_moved");
  assert.match(ranAway.detail, /better than Unabated/);
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

test("a thin book better at the number still makes not_best", () => {
  const state = boardWithEdge(WSU_ALT_KEY, 0.05);
  const novigKey = Object.keys(state.lines).find((key) => key.includes("ms89") && key.endsWith("si1:tid717:alt-35.5")) || "420942308:ms89:si1:tid717";
  moveLine(state, novigKey, { points: -35.5, price: -120, liquidity: 12, statusId: 1, isBlurred: false, modifiedOn: new Date(NOW - 1000).toISOString(), sequenceNumber: NOW - 1000 });
  const decision = decide(state, WSU_ALT_KEY, { ask: 0.57, depth: 1000 });
  assert.equal(decision.reason, "not_best");
  assert.match(decision.detail, /Novig -120/);
});

test("the decision signature changes with the ask, not with the clock", () => {
  const base = { reason: "bet", ticker: "T", liveAsk: 0.57, unabatedPrice: -143, contracts: 3, liveEdge: 0.05 };
  assert.equal(taker.decisionSignature(base), taker.decisionSignature({ ...base, liveEdge: 0.051 }));
  assert.notEqual(taker.decisionSignature(base), taker.decisionSignature({ ...base, liveAsk: 0.58 }));
  // A skipped row does not re-log on every tick of the book: only its reason or Unabated's price matter.
  const skipped = { ...base, reason: "price_moved" };
  assert.equal(taker.decisionSignature(skipped), taker.decisionSignature({ ...skipped, liveAsk: 0.58, contracts: 4 }));
  assert.notEqual(taker.decisionSignature(skipped), taker.decisionSignature({ ...skipped, unabatedPrice: -150 }));
});

test("a derived moneyline trades only when Kalshi's market names the row's team and not the other", () => {
  const state = boardWithEdge(ML_WSU_KEY, 0.05, { price: -400 });
  const row = rowOf(state, ML_WSU_KEY);
  assert.equal(taker.moneylineTeamMatches({ yes_sub_title: "Washington St." }, row), true);
  assert.equal(taker.moneylineTeamMatches({ yes_sub_title: "Duquesne" }, row), false);
  assert.equal(taker.moneylineTeamMatches({ yes_sub_title: "" }, row), false);
  assert.equal(taker.moneylineTeamMatches(null, row), false);
  // A shared first word names neither side for sure.
  const sameCity = { ...row, awayTeam: "Los Angeles Chargers", homeTeam: "Los Angeles Rams" };
  assert.equal(taker.moneylineTeamMatches({ yes_sub_title: "Los Angeles R" }, sameCity), false);
});

test("the market check: active, the row's strike, the series' period, the side's sign, and the YES title's team", () => {
  const state = boardWithEdge(WSU_ALT_KEY, 0.05);
  const wsuRow = rowOf(state, WSU_ALT_KEY);
  const ticker = "KXNCAAFSPREAD-26SEP19DUQWSU-WSU36";
  const yes = { side: "yes", ticker };
  const no = { side: "no", ticker };
  const spread = { ticker, status: "active", floor_strike: 35.5, yes_sub_title: "Washington St. wins by over 35.5 points" };
  assert.equal(taker.marketMismatch(spread, wsuRow, yes), null);
  assert.match(taker.marketMismatch({ ...spread, status: "closed" }, wsuRow, yes), /status "closed"/);
  assert.match(taker.marketMismatch({ ...spread, floor_strike: 34.5 }, wsuRow, yes), /floor_strike 34.5, row is -35.5/);
  assert.match(taker.marketMismatch({ ...spread, floor_strike: undefined }, wsuRow, yes), /floor_strike undefined/);
  assert.match(taker.marketMismatch(null, wsuRow, yes), /no market object/);
  // Kalshi answering another market, or a first-half series for a full-game row.
  assert.match(taker.marketMismatch({ ...spread, ticker: "KXNCAAFSPREAD-26SEP19DUQWSU-WSU42" }, wsuRow, yes), /answered KXNCAAFSPREAD-26SEP19DUQWSU-WSU42 for KXNCAAFSPREAD-26SEP19DUQWSU-WSU36/);
  const halfTicker = "KXNCAAF1HSPREAD-26SEP19DUQWSU-WSU36";
  assert.match(taker.marketMismatch({ ...spread, ticker: halfTicker }, wsuRow, { side: "yes", ticker: halfTicker }), /series KXNCAAF1HSPREAD is not the row's period/);
  // NO of a "wins by over" market is the dog's side: never for the favorite's row, and YES never for the dog's.
  assert.match(taker.marketMismatch(spread, wsuRow, no), /NO of a wins-by-over market for a -35.5 row/);
  const duqState = boardWithEdge(DUQ_ALT_KEY, 0.05);
  const duqRow = rowOf(duqState, DUQ_ALT_KEY);
  assert.equal(taker.marketMismatch(spread, duqRow, no), null);
  assert.match(taker.marketMismatch(spread, duqRow, yes), /YES of a wins-by-over market for a 35.5 row/);
  // The Duquesne +35.5 row buys NO of the same market: YES must still name Washington State.
  assert.match(taker.marketMismatch({ ...spread, yes_sub_title: "Duquesne wins by over 35.5 points" }, duqRow, no), /does not name side 1/);
  // A total: YES is Over, NO is Under, at the row's strike.
  const totalTicker = "KXNCAAFTOTAL-26SEP19DUQWSU-60";
  const total = { ticker: totalTicker, status: "active", floor_strike: 59.5, yes_sub_title: "Over 59.5 points scored" };
  const underRow = { ...wsuRow, betTypeId: 3, points: 59.5, sideIndex: 1, sideLabel: "Under" };
  const overRow = { ...underRow, sideIndex: 0, sideLabel: "Over" };
  assert.equal(taker.marketMismatch(total, underRow, { side: "no", ticker: totalTicker }), null);
  assert.equal(taker.marketMismatch(total, overRow, { side: "yes", ticker: totalTicker }), null);
  assert.match(taker.marketMismatch(total, underRow, { side: "yes", ticker: totalTicker }), /YES of an Over market for the Under row/);
  assert.match(taker.marketMismatch(total, overRow, { side: "no", ticker: totalTicker }), /NO of an Over market for the Over row/);
  assert.match(taker.marketMismatch({ ...total, floor_strike: 60.5 }, underRow, { side: "no", ticker: totalTicker }), /floor_strike 60.5/);
  // A moneyline: the team only.
  const mlState = boardWithEdge(ML_WSU_KEY, 0.05, { price: -400 });
  const mlRow = rowOf(mlState, ML_WSU_KEY);
  const mlTicker = "KXNCAAFGAME-26SEP19DUQWSU-WSU";
  assert.equal(taker.marketMismatch({ ticker: mlTicker, status: "active", yes_sub_title: "Washington St." }, mlRow, { side: "yes", ticker: mlTicker }), null);
  assert.match(taker.marketMismatch({ ticker: mlTicker, status: "active", yes_sub_title: "Duquesne" }, mlRow, { side: "yes", ticker: mlTicker }), /is not/);
});

// An open Kalshi bet of Cal's on the slice's game, as /bets.json lists it.
function heldKalshi(ticker, side, betType, points, contractSide = "yes") {
  return {
    id: `kalshi:${ticker}:${contractSide}`, source: "kalshi_api", venue: "kalshi", league: "cfb", status: "open",
    betType, period: "FG", side, points, stake: 100, toWin: 25, contracts: 100, price: -400,
    awayTeam: "Duquesne", homeTeam: "Washington St.", awayKey: null, homeKey: null,
    eventStart: "2026-09-19T19:30:00.000Z", eventDate: "2026-09-19", placedAt: "2026-09-18T19:30:00.000Z",
    isParlayLeg: false, parlayId: null, approx: [], unmatchable: null, rotation: null,
    venueIds: { marketTicker: ticker, eventTicker: ticker.split("-").slice(0, 2).join("-") }, raw: {},
  };
}

// decideRow with the runner's real sizing inputs: the held bets annotated
// against the row (bets.annotateRows) and the game's fair ladder.
function decideHolding(key, records) {
  const state = boardWithEdge(key, 0.05);
  const row = rowOf(state, key);
  for (const [league, list] of Object.entries(edgeRows.feedTeamsByLeague(state))) teams.registerTeams(league, list);
  const flags = betsLib.annotateRows([row], betsLib.resolveTeamKeys(records, []), { lines: edgeRows.boardLines(state) });
  const ladderOf = edgeRows.createLadderReaders(state)(row.eventId);
  return decide(state, key, { ask: 0.57, depth: 100000 }, { sizing: { matches: flags[0].matches, ladderOf, teasers: [] } });
}

test("a held bet on the row's axis the sizing math could not price refuses the stake; a held total on a spread row does not", () => {
  // WSU -41.5 held, deciding WSU -35.5: the same side at a number the fair ladder lacks.
  const sameSide = decideHolding(WSU_ALT_KEY, [heldKalshi("KXNCAAFSPREAD-26SEP19DUQWSU-WSU42", "home", "spread", -41.5)]);
  assert.equal(sameSide.reason, "size_declined");
  assert.match(sameSide.detail, /held .*-41\.5.* is not in the sizing math: no fair at -41\.5/);
  // A held total is another market (#129): left out on purpose, the spread row still bets.
  const otherMarket = decideHolding(WSU_ALT_KEY, [heldKalshi("KXNCAAFTOTAL-26SEP19DUQWSU-60", "over", "total", 59.5)]);
  assert.equal(otherMarket.reason, "bet", otherMarket.detail);
  // The same line held: conditional Kelly sizes it (here to nothing more), no refusal.
  const sameLine = decideHolding(WSU_ALT_KEY, [heldKalshi("KXNCAAFSPREAD-26SEP19DUQWSU-WSU36", "home", "spread", -35.5)]);
  assert.equal(sameLine.reason, "size_zero", sameLine.detail);
  // A smaller same-line holding leaves room: $20 at -133 held, the bot adds $15.04. At the old $500
  // bankroll this sized to $0 too, which is why the bot never bet a game Cal already held.
  const topUp = decideHolding(WSU_ALT_KEY, [{ ...heldKalshi("KXNCAAFSPREAD-26SEP19DUQWSU-WSU36", "home", "spread", -35.5), stake: 20, toWin: 15, contracts: 35, price: -133 }]);
  assert.deepEqual({ reason: topUp.reason, stake: topUp.stake, contracts: topUp.contracts, cost: topUp.cost }, { reason: "bet", stake: 15.04, contracts: 25, cost: 14.6789 });
});

test("the other side of the market held is never bet against, in either account: size_declined", () => {
  const ticker = "KXNCAAFSPREAD-26SEP19DUQWSU-WSU36";
  const primary = decideHolding(WSU_ALT_KEY, [heldKalshi(ticker, "away", "spread", 35.5, "no")]);
  assert.equal(primary.reason, "size_declined");
  assert.match(primary.detail, /the NO side of KXNCAAFSPREAD-26SEP19DUQWSU-WSU36 is held in the primary account: the bot never takes the other side/);
  // The bot's own subaccount: the id carries :s<N> (bets_service/sources/kalshi.py record_id_of).
  const botHeld = { ...heldKalshi(ticker, "away", "spread", 35.5, "no"), id: `kalshi:${ticker}:no:s1` };
  const bot = decideHolding(WSU_ALT_KEY, [botHeld]);
  assert.equal(bot.reason, "size_declined");
  assert.match(bot.detail, /is held in the bot's subaccount/);
  // A longer ticker sharing the prefix is another market, not the other side.
  assert.equal(taker.heldOtherSide([{ bet: { id: `kalshi:${ticker}:nox`, status: "open" } }], { ticker, side: "yes" }), null);
  assert.equal(taker.heldOtherSide([{ bet: { id: `kalshi:${ticker}:no:s1`, status: "lost" } }], { ticker, side: "yes" }), null);
});

test("every bot order id starts ubot- and is fixed per ticker, side and pass", () => {
  const id = taker.clientOrderIdOf("KXNCAAFSPREAD-26SEP19DUQWSU-WSU36", "yes", 1000);
  assert.match(id, /^ubot-[0-9a-f]{32}$/);
  assert.equal(taker.clientOrderIdOf("KXNCAAFSPREAD-26SEP19DUQWSU-WSU36", "yes", 1000), id);
  assert.notEqual(taker.clientOrderIdOf("KXNCAAFSPREAD-26SEP19DUQWSU-WSU36", "yes", 31000), id);
  assert.notEqual(taker.clientOrderIdOf("KXNCAAFSPREAD-26SEP19DUQWSU-WSU36", "no", 1000), id);
});

// ---- Kalshi's own start for a game (second clock beside Unabated's eventStart) ----
const CLMB_EVENT = "KXNCAAFSPREAD-26OCT16CLMBPENN";

test("kalshiStartMsOf: the earliest zoned start_date among milestones naming the event", () => {
  const body = { milestones: [
    { start_date: "2026-10-17T00:30:00Z", related_event_tickers: [CLMB_EVENT] },
    { start_date: "2026-10-16T23:00:00Z", primary_event_tickers: [CLMB_EVENT] },
    { start_date: "2026-10-16T20:00:00Z", related_event_tickers: ["KXNCAAFSPREAD-26OCT16OTHER"] },
    { start_date: "2026-10-16T21:00:00", related_event_tickers: [CLMB_EVENT] },
    { start_date: "not a time", related_event_tickers: [CLMB_EVENT] },
    null,
  ] };
  assert.equal(taker.kalshiStartMsOf(body, CLMB_EVENT), Date.parse("2026-10-16T23:00:00Z"));
  assert.equal(taker.kalshiStartMsOf({ milestones: [] }, CLMB_EVENT), null);
  assert.equal(taker.kalshiStartMsOf({}, CLMB_EVENT), null);
  assert.equal(taker.kalshiStartMsOf(null, CLMB_EVENT), null);
});

test("earliestStartMs: the earlier clock wins; Unabated's alone when Kalshi has none", () => {
  const unabated = Date.parse("2026-10-17T00:30:00Z");
  const kalshi = Date.parse("2026-10-16T23:00:00Z");
  assert.equal(taker.earliestStartMs(unabated, kalshi), kalshi);
  assert.equal(taker.earliestStartMs(kalshi, unabated), kalshi);
  assert.equal(taker.earliestStartMs(unabated, null), unabated);
  assert.equal(taker.earliestStartMs(NaN, null), null);
});

test("eventTickerOf: Kalshi's event_ticker, else the ticker less its last part", () => {
  assert.equal(taker.eventTickerOf({ event_ticker: CLMB_EVENT }, "X-Y-Z"), CLMB_EVENT);
  assert.equal(taker.eventTickerOf({}, `${CLMB_EVENT}-CLMB6`), CLMB_EVENT);
  assert.equal(taker.eventTickerOf(null, "NODASH"), null);
});
