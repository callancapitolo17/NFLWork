// server/taker.js, Novig: the outcome off Unabated's rung, the price off
// Novig's order-by-order book, 1c contracts, the market check, the other side
// of a held Novig market. Board: fixtures/v2_venue_ids_slice.json (CFB event
// 123742, Duquesne @ Washington State, with Novig rung outcome ids).
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
// Novig WSU -35.5 (main line): -127, Novig's own 0.56.
const WSU_KEY = "420942308:ms89:si1:tid717";
const WSU_OUTCOME = "01a0a05e-e864-7b73-9b80-9e502e4e98d5";
const DUQ_OUTCOME = "01a0a05e-e864-7b73-9b80-0000000000aa";
const MARKET_ID = "0b0b0b0b-1111-4222-8333-444444444444";
const NOVIG_EVENT_ID = "0e0e0e0e-1111-4222-8333-555555555555";
const EDGE_SETTINGS = { leagues: [CFB], periods: [1] };

function board(fields = {}) {
  const json = JSON.parse(fs.readFileSync(path.join(__dirname, "fixtures", "v2_venue_ids_slice.json"), "utf8"));
  const state = feed.parseSnapshot(json, { leagueId: CFB });
  moveLine(state, WSU_KEY, { ge: 0.06, bacr: -200, statusId: 1, sequenceNumber: NOW - 60 * 1000, modifiedOn: new Date(NOW - 60 * 1000).toISOString(), ...fields });
  return state;
}

function novigRow(state) {
  const row = taker.candidateRows(state, EDGE_SETTINGS, NOW, taker.NOVIG_BOOK_ID).find((candidate) => candidate.key === WSU_KEY);
  assert.ok(row, "expected the Novig WSU -35.5 row among the candidates");
  return row;
}

function market(fields = {}) {
  return {
    marketId: MARKET_ID, eventId: NOVIG_EVENT_ID, marketType: "SPREAD", strike: "-35.5", status: "OPEN",
    fee: { coefficient: "0.06", makerCredit: "0", charged: "WHEN_LIVE" },
    outcomes: [{ outcomeId: WSU_OUTCOME, name: "Washington St -35.5", status: "TBD" }, { outcomeId: DUQ_OUTCOME, name: "Duquesne +35.5", status: "TBD" }],
    ...fields,
  };
}

const EVENT = { eventId: NOVIG_EVENT_ID, status: "OPEN_PREGAME", league: "NCAAF", startsTs: EVENT_START_MS };

function contractOf(row, held = market()) {
  return taker.novigContractOf(row, taker.novigOutcomeIndex([held]));
}

// A book whose best resting order on Duquesne makes WSU cost `ask`.
function bookWithAsk(ask, qty = 100000) {
  return { marketId: MARKET_ID, seq: 1, orders: { [DUQ_OUTCOME]: [{ orderId: "a", price: (1 - ask).toFixed(3), qty }, { orderId: "b", price: (1 - ask - 0.01).toFixed(3), qty: 5 }], [WSU_OUTCOME]: [] } };
}

function decide(state, quote, sizing = { matches: [], ladderOf: () => null, teasers: [] }) {
  const row = novigRow(state);
  return taker.decideRow({
    row, contract: contractOf(row), quote, priceIndex: pricecheck.buildPriceIndex(state, { now: NOW, maxLineAgeMs: taker.MAX_LINE_AGE_MS }),
    sizing, venue: taker.VENUE_NOVIG,
  });
}

test("the outcome is the one Unabated's rung puts on this line, placed in its two-outcome market", () => {
  const row = novigRow(board());
  assert.deepEqual(taker.resolveNovigOutcome(row), { outcomeId: WSU_OUTCOME });
  assert.deepEqual(contractOf(row), { outcomeId: WSU_OUTCOME, marketId: MARKET_ID, otherOutcomeId: DUQ_OUTCOME });
  assert.equal(taker.novigContractOf(row, new Map()).reason, "no_outcome");
  assert.equal(contractOf(row, market({ outcomes: [{ outcomeId: WSU_OUTCOME }] })).reason, "no_outcome");
});

test("the ask is 1 minus the best resting order on the other outcome; depth sums that price", () => {
  const contract = { outcomeId: WSU_OUTCOME, otherOutcomeId: DUQ_OUTCOME };
  const book = { orders: { [DUQ_OUTCOME]: [{ price: "0.445", qty: 300 }, { price: "0.445", qty: 200 }, { price: "0.440", qty: 9000 }], [WSU_OUTCOME]: [{ price: "0.5", qty: 1 }] } };
  assert.deepEqual(taker.novigQuoteOf(book, contract), { ask: 0.555, depth: 500 });
  assert.equal(taker.novigQuoteOf({ orders: {} }, contract), null);
  assert.throws(() => taker.novigQuoteOf({}, contract), /expected \{orders\}/);
});

test("a clean Novig row bets in 1c contracts with no fee pregame", () => {
  const decision = decide(board(), taker.novigQuoteOf(bookWithAsk(0.555), { otherOutcomeId: DUQ_OUTCOME }));
  assert.equal(decision.reason, "bet", decision.detail);
  assert.equal(decision.venue, "novig");
  assert.equal(decision.ticker, WSU_OUTCOME);
  assert.equal(decision.liveAllIn, 0.555);
  // contracts pay 1c: the stake buys stake x 100 / ask of them, never Kalshi's stake / ask.
  assert.equal(decision.contracts, Math.floor((decision.stake * 100) / 0.555));
  assert.ok(Math.abs(decision.cost - (decision.contracts * 0.555) / 100) < 1e-4);
});

test("Novig's live price must be Unabated's: worse by over 0.1c or better by over 1c skips", () => {
  assert.equal(decide(board(), { ask: 0.565, depth: 100000 }).reason, "price_moved");
  assert.equal(decide(board(), { ask: 0.545, depth: 100000 }).reason, "price_moved");
  assert.equal(decide(board(), null).reason, "no_book");
});

test("the other outcome of a held Novig market declines the stake", () => {
  const held = { bet: { id: "novig:o1", venue: "novig", status: "open", venueIds: { marketId: MARKET_ID, outcomeId: DUQ_OUTCOME } } };
  const reason = taker.novigHeldOtherSide([held], { marketId: MARKET_ID, outcomeId: WSU_OUTCOME });
  assert.match(reason, /other outcome/);
  assert.equal(taker.novigHeldOtherSide([{ ...held, bet: { ...held.bet, venueIds: { marketId: MARKET_ID, outcomeId: WSU_OUTCOME } } }],
    { marketId: MARKET_ID, outcomeId: WSU_OUTCOME }), null);
});

test("the market check: open, pregame, no pregame fee, the row's type, period and home strike", () => {
  const row = novigRow(board());
  const contract = contractOf(row);
  assert.equal(taker.novigMarketMismatch(market(), EVENT, row, contract), null);
  assert.match(taker.novigMarketMismatch(market({ status: "CLOSED" }), EVENT, row, contract), /CLOSED/);
  assert.match(taker.novigMarketMismatch(market(), { ...EVENT, status: "OPEN_INGAME" }, row, contract), /OPEN_INGAME/);
  assert.match(taker.novigMarketMismatch(market({ fee: { charged: "ALWAYS" } }), EVENT, row, contract), /fee/);
  assert.match(taker.novigMarketMismatch(market({ strike: "35.5" }), EVENT, row, contract), /strike/);
  assert.match(taker.novigMarketMismatch(market({ marketType: "SPREAD_1H" }), EVENT, row, contract), /period/);
  assert.match(taker.novigMarketMismatch(market({ marketType: "TOTAL" }), EVENT, row, contract), /type/);
});

test("the outcome's own name must be the row's side: a swapped rung fails closed", () => {
  const row = novigRow(board());
  const contract = contractOf(row);
  const named = (wsuName) => market({ outcomes: [{ outcomeId: WSU_OUTCOME, name: wsuName, status: "TBD" }, { outcomeId: DUQ_OUTCOME, name: "Duquesne +35.5", status: "TBD" }] });
  assert.match(taker.novigMarketMismatch(named("Washington St +35.5"), EVENT, row, contract), /is 35.5, the row is -35.5/);
  assert.match(taker.novigMarketMismatch(named("Washington St"), EVENT, row, contract), /names no handicap/);
  assert.equal(taker.novigMarketMismatch(named("WSU -35.5"), EVENT, row, contract), null);
  assert.match(taker.novigMarketMismatch(named("WSU PK"), EVENT, row, contract), /pick'em/);
});

test("the clientId is a UUID per outcome and pass, and the line is named the way a bet record is", () => {
  const id = taker.novigClientIdOf(WSU_OUTCOME, NOW);
  assert.match(id, /^[0-9a-f]{8}-[0-9a-f]{4}-4[0-9a-f]{3}-8[0-9a-f]{3}-[0-9a-f]{12}$/);
  assert.equal(id, taker.novigClientIdOf(WSU_OUTCOME, NOW));
  assert.notEqual(id, taker.novigClientIdOf(WSU_OUTCOME, NOW + 30000));
  assert.deepEqual(taker.novigLineOf(novigRow(board())), {
    league: "cfb", betType: "spread", period: "FG", side: "home", points: -35.5,
    awayTeam: "Duquesne", homeTeam: "Washington State", awayTeamId: "927", homeTeamId: "717",
  });
});
