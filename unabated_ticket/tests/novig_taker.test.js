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

function decide(state, quote, sizing = { matches: [], ladderOf: () => null, teasers: [], shrink: 1, minShrunkEdgePct: 1 }) {
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

test("a Novig row meets the Min shrunk edge floor too (issue #208): the shrunk live edge gates it, the stake stays on the raw live edge", () => {
  const quote = taker.novigQuoteOf(bookWithAsk(0.555), { otherOutcomeId: DUQ_OUTCOME });
  const full = decide(board(), quote);
  const half = decide(board(), quote, { matches: [], ladderOf: () => null, teasers: [], shrink: 0.5, minShrunkEdgePct: 1 });
  assert.deepEqual([full.reason, half.reason, full.shrink, half.shrink], ["bet", "bet", 1, 0.5]);
  assert.ok(Math.abs(half.shrunkEdge - full.liveEdge / 2) < 1e-4, `${half.shrunkEdge} vs ${full.liveEdge}`);
  assert.deepEqual([half.stake, half.contracts], [full.stake, full.contracts]);
  const deep = decide(board(), quote, { matches: [], ladderOf: () => null, teasers: [], shrink: 0.1, minShrunkEdgePct: 1 });
  assert.equal(deep.reason, "edge_short");
  assert.match(deep.detail, /shrunk live edge .* under 1%/);
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

// ---- main numbers: the feed's ladder leaves out the main rung ----------------
// The fixture's Novig 52.5 total hangs 50.5/51.5/53.5/54.5 and no 52.5 rung,
// so the main line's outcome is found in Novig's catalog off its siblings.
const UNDER_KEY = "433644483:ms89:si1:tid717";
const OVER_KEY = "433644484:ms89:si0:tid927";
const DUQ_KEY = "420942310:ms89:si0:tid927";
const UNDER_515 = "01a0a059-6836-7b83-8813-18d718a0283c";
const OVER_515 = "01a0a059-6836-7b83-8813-18c4e5a44a2f";
const UNDER_535 = "01a0a059-6836-7b83-8813-19060702da02";
const OVER_535 = "01a0a059-6836-7b83-8813-18f82e7d1e57";
const TOTAL_MAIN_MARKET = "0c0c0c0c-1111-4222-8333-525252525252";
const TOTAL_MAIN_OVER = "0d0d0d0d-1111-4222-8333-0000000000a1";
const TOTAL_MAIN_UNDER = "0d0d0d0d-1111-4222-8333-0000000000a2";

function lineRow(state, key) {
  return feed.describeLine(state.lines[key], state);
}

function totalMarket(marketId, strike, over, under, fields = {}) {
  return market({
    marketId, marketType: "TOTAL", strike,
    outcomes: [{ outcomeId: over, name: `Over ${strike}`, status: "TBD" }, { outcomeId: under, name: `Under ${strike}`, status: "TBD" }],
    ...fields,
  });
}

function totalCatalog(extra = []) {
  return taker.novigOutcomeIndex([
    totalMarket("0c0c0c0c-1111-4222-8333-515151515151", "51.5", OVER_515, UNDER_515),
    totalMarket("0c0c0c0c-1111-4222-8333-535353535353", "53.5", OVER_535, UNDER_535),
    totalMarket(TOTAL_MAIN_MARKET, "52.5", TOTAL_MAIN_OVER, TOTAL_MAIN_UNDER),
    ...extra,
  ]);
}

test("a main total with no rung of its own takes the one catalog market of its siblings' event and type at its number", () => {
  const state = board();
  const under = lineRow(state, UNDER_KEY);
  assert.equal(taker.resolveNovigOutcome(under).reason, "no_outcome");
  const contract = taker.novigContractOf(under, totalCatalog());
  assert.deepEqual(contract, { outcomeId: TOTAL_MAIN_UNDER, marketId: TOTAL_MAIN_MARKET, otherOutcomeId: TOTAL_MAIN_OVER, derived: true });
  assert.equal(taker.novigMarketMismatch(totalMarket(TOTAL_MAIN_MARKET, "52.5", TOTAL_MAIN_OVER, TOTAL_MAIN_UNDER), EVENT, under, contract), null);
  const over = taker.novigContractOf(lineRow(state, OVER_KEY), totalCatalog());
  assert.equal(over.outcomeId, TOTAL_MAIN_OVER);
  assert.equal(over.otherOutcomeId, TOTAL_MAIN_UNDER);
});

test("a main number refuses when the catalog has no market, two markets, or names that cannot tell the side", () => {
  const under = lineRow(board(), UNDER_KEY);
  const withoutMain = taker.novigOutcomeIndex([totalMarket("0c0c0c0c-1111-4222-8333-515151515151", "51.5", OVER_515, UNDER_515)]);
  assert.match(taker.novigContractOf(under, withoutMain).detail, /0 open Novig TOTAL markets at strike 52.5/);
  const twin = totalMarket("0c0c0c0c-1111-4222-8333-999999999999", "52.5", "0d0d0d0d-1111-4222-8333-0000000000b1", "0d0d0d0d-1111-4222-8333-0000000000b2");
  assert.match(taker.novigContractOf(under, totalCatalog([twin])).detail, /2 open Novig TOTAL markets/);
  // Another period's market at the same number is another type: not a twin.
  const firstHalf = totalMarket("0c0c0c0c-1111-4222-8333-888888888888", "52.5", "0d0d0d0d-1111-4222-8333-0000000000c1", "0d0d0d0d-1111-4222-8333-0000000000c2", { marketType: "TOTAL_1H" });
  assert.equal(taker.novigContractOf(under, totalCatalog([firstHalf])).outcomeId, TOTAL_MAIN_UNDER);
  const unnamed = taker.novigOutcomeIndex([
    totalMarket("0c0c0c0c-1111-4222-8333-515151515151", "51.5", OVER_515, UNDER_515),
    market({ marketId: TOTAL_MAIN_MARKET, marketType: "TOTAL", strike: "52.5", outcomes: [{ outcomeId: TOTAL_MAIN_OVER, name: "Yes" }, { outcomeId: TOTAL_MAIN_UNDER, name: "No" }] }),
  ]);
  assert.match(taker.novigContractOf(under, unnamed).detail, /0 outcomes of market .* name the row's side/);
  assert.match(taker.novigContractOf(under, new Map()).detail, /no rung on its ladder is on an open pregame Novig market/);
});

test("a main number refuses when its ladder's rungs sit on different Novig events or market types", () => {
  const under = lineRow(board(), UNDER_KEY);
  const otherEvent = taker.novigOutcomeIndex([
    totalMarket("0c0c0c0c-1111-4222-8333-515151515151", "51.5", OVER_515, UNDER_515),
    totalMarket("0c0c0c0c-1111-4222-8333-535353535353", "53.5", OVER_535, UNDER_535, { eventId: "0e0e0e0e-1111-4222-8333-000000000000" }),
    totalMarket(TOTAL_MAIN_MARKET, "52.5", TOTAL_MAIN_OVER, TOTAL_MAIN_UNDER),
  ]);
  assert.match(taker.novigContractOf(under, otherEvent).detail, /rungs disagree: 2 events, 1 market types/);
});

test("a main spread missing its rung finds the market at the home side's handicap and the side by its name", () => {
  const state = board();
  // Drop the main-number rungs the fixture does carry, as the feed leaves them out.
  const outcomes = state.events[123742].venueIds.novigOutcomes;
  delete outcomes[WSU_OUTCOME];
  delete outcomes["01a0a05e-e864-7b73-9b80-9e64b163cc27"];
  const siblings = market({
    marketId: "0b0b0b0b-1111-4222-8333-343434343434", strike: "-34.5",
    outcomes: [{ outcomeId: "01a0a05e-e864-7b73-9b80-9fde29b0f1dc", name: "Washington St -34.5" }, { outcomeId: "01a0a05e-e864-7b73-9b80-9fe2c39c2f8d", name: "Duquesne +34.5" }],
  });
  const index = taker.novigOutcomeIndex([siblings, market()]);
  const duquesne = lineRow(state, DUQ_KEY);
  const contract = taker.novigContractOf(duquesne, index);
  assert.deepEqual(contract, { outcomeId: DUQ_OUTCOME, marketId: MARKET_ID, otherOutcomeId: WSU_OUTCOME, derived: true });
  assert.equal(taker.novigMarketMismatch(market(), EVENT, duquesne, contract), null);
  assert.equal(taker.novigContractOf(lineRow(state, WSU_KEY), index).outcomeId, WSU_OUTCOME);
  // A pick'em's names cannot tell the side: refused.
  const pick = market({ outcomes: [{ outcomeId: WSU_OUTCOME, name: "Washington St PK" }, { outcomeId: DUQ_OUTCOME, name: "Duquesne PK" }] });
  assert.equal(taker.novigContractOf(duquesne, taker.novigOutcomeIndex([siblings, pick])).reason, "no_outcome");
});

test("an alt line or a line its own rung names is never derived", () => {
  const state = board();
  const wsu = lineRow(state, WSU_KEY);
  // The rung on the main number is used as is, even with a catalog that lists more.
  assert.equal(taker.novigContractOf(wsu, taker.novigOutcomeIndex([market()])).derived, undefined);
  // An alt with no rung stays refused: only a main number has a market to look up.
  const alt = { ...lineRow(state, UNDER_KEY), isAlt: true, key: `${UNDER_KEY}:alt52.5` };
  assert.equal(taker.novigContractOf(alt, totalCatalog()).detail, "no Novig rung names this line");
});

test("a main spread refuses when Novig lists the game with home and away swapped, or the rungs' names disagree", () => {
  const state = board();
  const outcomes = state.events[123742].venueIds.novigOutcomes;
  delete outcomes[WSU_OUTCOME];
  delete outcomes["01a0a05e-e864-7b73-9b80-9e64b163cc27"];
  const duquesne = lineRow(state, DUQ_KEY);
  // Neutral site: Novig's home is Duquesne, so every strike is Duquesne's
  // handicap. Picking by number and name alone would buy Washington St -35.5
  // for the Duquesne +35.5 row; the siblings' strikes give the swap away.
  const swappedSibling = market({
    marketId: "0b0b0b0b-1111-4222-8333-343434343434", strike: "34.5",
    outcomes: [{ outcomeId: "01a0a05e-e864-7b73-9b80-9fde29b0f1dc", name: "Washington St -34.5" }, { outcomeId: "01a0a05e-e864-7b73-9b80-9fe2c39c2f8d", name: "Duquesne +34.5" }],
  });
  const swappedMain = market({ strike: "35.5" });
  const swapped = taker.novigContractOf(duquesne, taker.novigOutcomeIndex([swappedSibling, swappedMain]));
  assert.equal(swapped.reason, "no_outcome");
  assert.match(swapped.detail, /at strike 34.5, Unabated as 34.5 side 0/);
  // A main market whose outcome at the row's number names another team than the rungs do.
  const sibling = market({
    marketId: "0b0b0b0b-1111-4222-8333-343434343434", strike: "-34.5",
    outcomes: [{ outcomeId: "01a0a05e-e864-7b73-9b80-9fde29b0f1dc", name: "Washington St -34.5" }, { outcomeId: "01a0a05e-e864-7b73-9b80-9fe2c39c2f8d", name: "Duquesne +34.5" }],
  });
  const otherTeam = market({ outcomes: [{ outcomeId: WSU_OUTCOME, name: "Washington St -35.5" }, { outcomeId: DUQ_OUTCOME, name: "Idaho +35.5" }] });
  assert.match(taker.novigContractOf(duquesne, taker.novigOutcomeIndex([sibling, otherTeam])).detail, /0 outcomes of market/);
  // A rung whose own name is the other side.
  const misnamed = market({
    marketId: "0b0b0b0b-1111-4222-8333-343434343434", strike: "-34.5",
    outcomes: [{ outcomeId: "01a0a05e-e864-7b73-9b80-9fde29b0f1dc", name: "Duquesne +34.5" }, { outcomeId: "01a0a05e-e864-7b73-9b80-9fe2c39c2f8d", name: "Washington St -34.5" }],
  });
  assert.match(taker.novigContractOf(duquesne, taker.novigOutcomeIndex([misnamed, market()])).detail, /Novig lists rung/);
});

// ---- moneylines (the Ticket's Place only) ------------------------------------
// Unabated carries no Novig id on a moneyline: the game's other rungs find
// Novig's event, period and team names, and the event's one MONEY market is it.
const ML_WSU_KEY = "411181220:ms105:si1:tid717";
const MONEY_MARKET = "0f0f0f0f-1111-4222-8333-111111111111";
const MONEY_WSU = "0f0f0f0f-1111-4222-8333-0000000000e1";
const MONEY_DUQ = "0f0f0f0f-1111-4222-8333-0000000000e2";

function moneyMarket(fields = {}) {
  return market({
    marketId: MONEY_MARKET, marketType: "MONEY", strike: undefined,
    outcomes: [{ outcomeId: MONEY_DUQ, name: "Duquesne", status: "TBD" }, { outcomeId: MONEY_WSU, name: "Washington St", status: "TBD" }],
    ...fields,
  });
}

function moneylineRow(state, sideIndex = 1) {
  return { ...lineRow(state, ML_WSU_KEY), sideIndex };
}

test("a moneyline takes its event's one MONEY market, each outcome told by the team's name on Novig's spread rungs", () => {
  const state = board();
  const index = taker.novigOutcomeIndex([market(), moneyMarket()]);
  const wsu = taker.novigMoneylineContract(moneylineRow(state, 1), index, state);
  assert.deepEqual(wsu, { outcomeId: MONEY_WSU, marketId: MONEY_MARKET, otherOutcomeId: MONEY_DUQ, teamText: "washington st", derived: true });
  assert.equal(taker.novigMarketMismatch(moneyMarket(), EVENT, moneylineRow(state, 1), wsu), null);
  assert.deepEqual(taker.novigMoneylineContract(moneylineRow(state, 0), index, state),
    { outcomeId: MONEY_DUQ, marketId: MONEY_MARKET, otherOutcomeId: MONEY_WSU, teamText: "duquesne", derived: true });
  // The market check reads the outcome's name again: renamed since the catalog read, refused.
  const renamed = moneyMarket({ outcomes: [{ outcomeId: MONEY_DUQ, name: "Duquesne" }, { outcomeId: MONEY_WSU, name: "Wazzu" }] });
  assert.match(taker.novigMarketMismatch(renamed, EVENT, moneylineRow(state, 1), wsu), /is not the row's team/);
  // The bot never resolves one: its path refuses as before.
  assert.equal(taker.novigContractOf(moneylineRow(state, 1), index).reason, "no_outcome");
});

test("with no spread rung on Novig, a moneyline's outcomes are told by Unabated's team names", () => {
  const state = board();
  const index = taker.novigOutcomeIndex([totalMarket("0c0c0c0c-1111-4222-8333-515151515151", "51.5", OVER_515, UNDER_515), moneyMarket()]);
  assert.deepEqual(taker.novigMoneylineContract(moneylineRow(state, 1), index, state),
    { outcomeId: MONEY_WSU, marketId: MONEY_MARKET, otherOutcomeId: MONEY_DUQ, teamText: "washington st", derived: true });
  const strangers = moneyMarket({ outcomes: [{ outcomeId: MONEY_DUQ, name: "Team A" }, { outcomeId: MONEY_WSU, name: "Team B" }] });
  const unnamed = taker.novigOutcomeIndex([totalMarket("0c0c0c0c-1111-4222-8333-515151515151", "51.5", OVER_515, UNDER_515), strangers]);
  assert.match(taker.novigMoneylineContract(moneylineRow(state, 1), unnamed, state).detail, /do not name Duquesne and Washington State one each/);
});

test("a moneyline refuses on no market, two markets, a draw, names that disagree, or a spread rung read swapped", () => {
  const state = board();
  const row = moneylineRow(state, 1);
  const refused = (markets) => taker.novigMoneylineContract(row, taker.novigOutcomeIndex(markets), state);
  assert.match(refused([market()]).detail, /0 open Novig MONEY markets/);
  assert.match(refused([market(), moneyMarket(), moneyMarket({ marketId: "0f0f0f0f-1111-4222-8333-222222222222", outcomes: [{ outcomeId: "x1", name: "Duquesne" }, { outcomeId: "x2", name: "Washington St" }] })]).detail, /2 open Novig MONEY markets/);
  assert.match(refused([market(), moneyMarket({ outcomes: [...moneyMarket().outcomes, { outcomeId: "x3", name: "Draw" }] })]).detail, /3 outcomes/);
  // Both outcomes read as one team, or neither as either: refused.
  assert.match(refused([market(), moneyMarket({ outcomes: [{ outcomeId: MONEY_DUQ, name: "Washington St" }, { outcomeId: MONEY_WSU, name: "Washington St" }] })]).detail, /do not name/);
  assert.match(refused([market(), moneyMarket({ outcomes: [{ outcomeId: MONEY_DUQ, name: "Team A" }, { outcomeId: MONEY_WSU, name: "Team B" }] })]).detail, /do not name/);
  // One outcome told by Novig's spread name ("Washington St"), the other is the other team.
  assert.equal(refused([market(), moneyMarket({ outcomes: [{ outcomeId: MONEY_DUQ, name: "DUQ" }, { outcomeId: MONEY_WSU, name: "Washington St" }] })]).outcomeId, MONEY_WSU);
  // A neutral site listed home/away swapped: the spread rung's strike gives it away.
  assert.match(refused([market({ strike: "35.5" }), moneyMarket()]).detail, /Novig lists rung/);
  assert.match(refused([]).detail, /no Novig rung on this game/);
  assert.equal(taker.novigMoneylineContract({ ...row, betTypeId: 2 }, new Map(), state).reason, "no_outcome");
});

test("a half's moneyline refuses: its winner market may pay a tie as a loss", () => {
  const state = board();
  const row = { ...moneylineRow(state, 1), periodTypeId: 2, period: "1H" };
  const half = taker.novigOutcomeIndex([market(), moneyMarket({ marketType: "MONEY_1H" })]);
  assert.deepEqual(taker.novigMoneylineContract(row, half, state),
    { reason: "no_outcome", detail: "only full-game moneylines place: a half's winner market may pay a tie as a loss" });
});

test("a pick'em spread rung is skipped; both spread names known tell a differently spelled outcome by elimination", () => {
  const state = board();
  // Put a pick'em rung on the WSU spread ladder: its name cannot tell the side, the -35.5 rung still can.
  const outcomes = state.events[123742].venueIds.novigOutcomes;
  outcomes["pk-wsu"] = { ...outcomes[WSU_OUTCOME], points: 0, lineKey: null };
  const pick = market({ marketId: "0b0b0b0b-1111-4222-8333-000000000000", strike: "0", outcomes: [{ outcomeId: "pk-wsu", name: "Washington St PK" }, { outcomeId: "pk-duq", name: "Duquesne PK" }] });
  const index = taker.novigOutcomeIndex([market(), pick, moneyMarket({ outcomes: [{ outcomeId: MONEY_DUQ, name: "DUQ" }, { outcomeId: MONEY_WSU, name: "Washington St" }] })]);
  assert.deepEqual(taker.novigMoneylineContract(moneylineRow(state, 0), index, state),
    { outcomeId: MONEY_DUQ, marketId: MONEY_MARKET, otherOutcomeId: MONEY_WSU, teamText: "duq", derived: true });
});

test("a moneyline whose outcome names only part of a team, or a team sharing a city, refuses rather than guess", () => {
  const state = board();
  const row = { ...moneylineRow(state, 1), awayTeam: "Miami (OH)", homeTeam: "Miami" };
  const totals = totalMarket("0c0c0c0c-1111-4222-8333-515151515151", "51.5", OVER_515, UNDER_515);
  const miami = moneyMarket({ outcomes: [{ outcomeId: MONEY_DUQ, name: "Miami" }, { outcomeId: MONEY_WSU, name: "Miami" }] });
  assert.match(taker.novigMoneylineContract(row, taker.novigOutcomeIndex([totals, miami]), state).detail, /do not name/);
});
