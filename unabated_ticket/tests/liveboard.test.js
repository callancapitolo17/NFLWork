// server/liveboard.js and taker.js's live Place helpers: a live line has no
// venue ids, so its contract comes from the game's remembered pregame ids.
// Board: fixtures/v2_venue_ids_slice.json (CFB event 123742, Duquesne @
// Washington State, Kalshi and Novig rung ids on its spread and total).
"use strict";

const test = require("node:test");
const assert = require("node:assert/strict");
const fs = require("node:fs");
const path = require("node:path");
const feed = require("../extension/feed.js");
const liveboard = require("../server/liveboard.js");
const taker = require("../server/taker.js");

const EVENT_ID = 123742;
const raw = JSON.parse(fs.readFileSync(path.join(__dirname, "fixtures", "v2_venue_ids_slice.json"), "utf8"));
const pregameBoard = feed.parseSnapshot(raw, { leagueId: 2 });
const pregameLines = Object.values(pregameBoard.lines).filter((line) => line.eventId === EVENT_ID);
const remembered = liveboard.rememberPregameVenueIds(pregameBoard.events[EVENT_ID], pregameLines);

// The same game in play: its rows moved to live with no alternates or ids, the spread main at `spreadPoints` (WSU's side).
function liveBoardOf(spreadPoints = -35.5, totalPoints = 48.5) {
  const json = JSON.parse(JSON.stringify(raw));
  const rows = json.odds["lg2:pt1:pregame"].filter((row) => row.eventId === EVENT_ID);
  json.odds["lg2:pt1:pregame"] = json.odds["lg2:pt1:pregame"].filter((row) => row.eventId !== EVENT_ID);
  json.odds["lg2:pt1:live"] = rows.map((row) => {
    const live = { ...row, key: row.key.replace(":pregame:", ":live:"), live: true };
    for (const [sideKey, books] of Object.entries(live.sides)) {
      for (const book of Object.values(books)) {
        Object.assign(book, { alternateLines: null, sourceKey: null, sourceData: null });
        if (row.key.includes(":bt2:")) book.points = sideKey === "si1:tid717" ? spreadPoints : -spreadPoints;
        if (row.key.includes(":bt3:")) book.points = totalPoints;
      }
    }
    return live;
  });
  return feed.parseSnapshot(json, { leagueId: 2, phase: "live" });
}

function liveRowOf(board, bookId, betTypeId, sideIndex) {
  const line = Object.values(board.lines).find((candidate) => candidate.eventId === EVENT_ID && candidate.bookId === bookId
    && candidate.betTypeId === betTypeId && candidate.sideIndex === sideIndex && candidate.periodTypeId === 1);
  return feed.describeLine(line, board);
}

test("a live board's lines carry no venue ids of their own", () => {
  const live = liveBoardOf();
  assert.deepEqual(Object.keys(live.events[EVENT_ID].venueIds.kalshiContracts), []);
  assert.ok(remembered, "the pregame game is remembered");
  assert.ok(Object.keys(remembered.venueIds.kalshiContracts).length > 0);
});

test("remembered Kalshi rungs re-pointed at the live main find its contract: on a rung's number by the rung, else derived", () => {
  const live = liveBoardOf(-35.5, 48.5);
  const board = liveboard.liveQuoteBoard(live, EVENT_ID, remembered);
  const spread = taker.resolveKalshiContract(liveRowOf(board, taker.KALSHI_BOOK_ID, 2, 1), board);
  assert.deepEqual({ ticker: spread.ticker, side: spread.side, derived: spread.derived }, { ticker: "KXNCAAFSPREAD-26SEP19DUQWSU-WSU36", side: "yes", derived: false });
  const over = taker.resolveKalshiContract(liveRowOf(board, taker.KALSHI_BOOK_ID, 3, 0), board);
  assert.deepEqual({ ticker: over.ticker, side: over.side, derived: over.derived }, { ticker: "KXNCAAFTOTAL-26SEP19DUQWSU-49", side: "yes", derived: true });
  const under = taker.resolveKalshiContract(liveRowOf(board, taker.KALSHI_BOOK_ID, 3, 1), board);
  assert.deepEqual({ ticker: under.ticker, side: under.side }, { ticker: "KXNCAAFTOTAL-26SEP19DUQWSU-49", side: "no" });
  // Without the memory the same live line names nothing.
  assert.ok(taker.resolveKalshiContract(liveRowOf(live, taker.KALSHI_BOOK_ID, 2, 1), live).reason);
});

test("the live board's own ids, if a feed ever sends them, win over remembered ones", () => {
  const live = liveBoardOf();
  const own = { lineKey: "own", mainKey: "own", points: -1.5, sideIndex: 1 };
  live.events[EVENT_ID].venueIds.kalshiContracts["Y-KXNCAAFSPREAD-26SEP19DUQWSU-WSU36"] = own;
  const board = liveboard.liveQuoteBoard(live, EVENT_ID, remembered);
  assert.deepEqual(board.events[EVENT_ID].venueIds.kalshiContracts["Y-KXNCAAFSPREAD-26SEP19DUQWSU-WSU36"], own);
  assert.equal(liveboard.liveQuoteBoard(live, EVENT_ID, null), live, "no memory: the live board as it is");
});

test("a live spread whose favorite flipped derives off the team codes: YES on the favorite's code, NO on the other's", () => {
  const board = liveboard.liveQuoteBoard(liveBoardOf(3.5), EVENT_ID, remembered);
  // The fixture's Kalshi spread rungs name only WSU's code: WSU +3.5 would be NO on "DUQ wins by over 3.5", which Kalshi does not list.
  const wsuDog = taker.derivedSpreadFromTeamCodes(liveRowOf(board, taker.KALSHI_BOOK_ID, 2, 1), board);
  assert.deepEqual({ reason: wsuDog.reason, detail: wsuDog.detail }, { reason: "no_ticker", detail: "no spread rung names side 0's team code" });
  const duq = liveRowOf(board, taker.KALSHI_BOOK_ID, 2, 0);
  const duqDog = taker.derivedSpreadFromTeamCodes({ ...duq, points: 3.5 }, board);
  assert.deepEqual({ ticker: duqDog.ticker, side: duqDog.side }, { ticker: "KXNCAAFSPREAD-26SEP19DUQWSU-WSU4", side: "no" });
  const wsuFavorite = taker.derivedSpreadFromTeamCodes({ ...liveRowOf(board, taker.KALSHI_BOOK_ID, 2, 1), points: -10.5 }, board);
  assert.deepEqual({ ticker: wsuFavorite.ticker, side: wsuFavorite.side }, { ticker: "KXNCAAFSPREAD-26SEP19DUQWSU-WSU11", side: "yes" });
  assert.equal(taker.derivedSpreadFromTeamCodes({ ...duq, points: 0 }, board).reason, "no_ticker", "a pick'em names no side");
});

test("Novig's fee rate: none pregame, the market's coefficient in game, unreadable refused", () => {
  const market = (fee) => ({ fee });
  assert.equal(taker.novigFeeRateOf(market({ charged: "WHEN_LIVE", coefficient: "0.06" }), false), 0);
  assert.equal(taker.novigFeeRateOf(market({ charged: "ALWAYS", coefficient: "0.06" }), false), null);
  assert.equal(taker.novigFeeRateOf(market({ charged: "WHEN_LIVE", coefficient: "0.06" }), true), 0.06);
  assert.equal(taker.novigFeeRateOf(market({ charged: "ALWAYS", coefficient: "0.03" }), true), 0.03);
  assert.equal(taker.novigFeeRateOf(market({ charged: "NEVER" }), true), 0);
  assert.equal(taker.novigFeeRateOf(market({ charged: "SOMETIMES", coefficient: "0.06" }), true), null);
  assert.equal(taker.novigFeeRateOf(market({ charged: "WHEN_LIVE", coefficient: "x" }), true), null);
  assert.equal(taker.novigFeeRateOf({}, true), null);
});

// Novig's catalog for the game: the WSU -35.5 spread rung (fixture outcome id) pregame,
// and in game a new SPREAD market at -10.5, a TOTAL at 48.5 and the MONEY market.
const NOVIG_EVENT = "ev-1";
const WSU_SPREAD_RUNG = "01a0a05e-e864-7b73-9b80-9e502e4e98d5";
function novigMarket(marketId, marketType, strike, names) {
  return { marketId, eventId: NOVIG_EVENT, marketType, strike, status: "OPEN", fee: { charged: "WHEN_LIVE", coefficient: "0.06" },
    outcomes: names.map((name, at) => ({ outcomeId: `${marketId}-o${at}`, name })) };
}

test("a pregame catalog read remembers the game's Novig event, its market types and how Novig names each side", () => {
  const pregameIndex = taker.novigOutcomeIndex([{ ...novigMarket("m-pre", "SPREAD", "-35.5", []), outcomes: [{ outcomeId: WSU_SPREAD_RUNG, name: "WSU -35.5" }, { outcomeId: "duq", name: "DUQ +35.5" }] }]);
  const memory = taker.novigEventMemoryOf(pregameBoard.events[EVENT_ID], pregameBoard, pregameIndex);
  assert.deepEqual(memory, { novigEventId: NOVIG_EVENT, marketTypes: { "1:2": "SPREAD" }, sideNames: [null, "wsu"] });
  assert.equal(taker.novigEventMemoryOf(pregameBoard.events[EVENT_ID], pregameBoard, new Map()), null, "no rung in the read: nothing to remember");
});

test("a live Novig line's outcome is found on the remembered event by type, strike and name; anything unsure refuses", () => {
  const memory = { novigEventId: NOVIG_EVENT, marketTypes: { "1:2": "SPREAD", "1:3": "TOTAL" }, sideNames: ["duq", "wsu"] };
  const index = taker.novigOutcomeIndex([
    novigMarket("m-spread", "SPREAD", "-10.5", ["WSU -10.5", "DUQ +10.5"]),
    novigMarket("m-total", "TOTAL", "48.5", ["Over 48.5", "Under 48.5"]),
    novigMarket("m-money", "MONEY", null, ["WSU", "DUQ"]),
    { ...novigMarket("m-other", "SPREAD", "-10.5", ["X -10.5", "Y +10.5"]), eventId: "ev-2" },
  ]);
  const row = (fields) => ({ awayTeam: "Duquesne", homeTeam: "Washington State", period: "FG", periodTypeId: 1, ...fields });
  assert.deepEqual(taker.novigLiveContract(row({ betTypeId: 2, betType: "spread", sideIndex: 1, points: -10.5 }), index, memory),
    { outcomeId: "m-spread-o0", marketId: "m-spread", otherOutcomeId: "m-spread-o1", derived: true });
  assert.deepEqual(taker.novigLiveContract(row({ betTypeId: 2, betType: "spread", sideIndex: 0, points: 10.5 }), index, memory).outcomeId, "m-spread-o1");
  assert.deepEqual(taker.novigLiveContract(row({ betTypeId: 3, betType: "total", sideIndex: 1, points: 48.5 }), index, memory).outcomeId, "m-total-o1");
  const money = taker.novigLiveContract(row({ betTypeId: 1, betType: "moneyline", sideIndex: 0, points: null }), index, memory);
  assert.deepEqual({ outcomeId: money.outcomeId, teamText: money.teamText }, { outcomeId: "m-money-o1", teamText: "duq" });
  // A number Novig has no market at, a half moneyline, an unremembered market type, no memory, an unnamed side.
  assert.match(taker.novigLiveContract(row({ betTypeId: 2, betType: "spread", sideIndex: 1, points: -9.5 }), index, memory).detail, /0 open in-game Novig SPREAD markets/);
  assert.equal(taker.novigLiveContract(row({ betTypeId: 1, betType: "moneyline", periodTypeId: 2, sideIndex: 0 }), index, memory).reason, "no_outcome");
  assert.match(taker.novigLiveContract(row({ betTypeId: 2, betType: "spread", periodTypeId: 2, period: "1H", sideIndex: 1, points: -10.5 }), index, memory).detail, /no market type for 1H spread/);
  assert.match(taker.novigLiveContract(row({ betTypeId: 3, sideIndex: 0, points: 48.5 }), index, null).detail, /no pregame Novig read/);
  assert.match(taker.novigLiveContract(row({ betTypeId: 2, sideIndex: 0, points: 10.5 }), index, { ...memory, sideNames: [null, "wsu"] }).detail, /named side 0's team/);
});

test("novigMarketMismatch live: the event must be in game and the fee priceable; pregame checks unchanged", () => {
  const market = novigMarket("m-total", "TOTAL", "48.5", ["Over 48.5", "Under 48.5"]);
  const row = { betTypeId: 3, periodTypeId: 1, sideIndex: 0, points: 48.5, betType: "total", period: "FG" };
  const contract = { marketId: "m-total", outcomeId: "m-total-o0" };
  assert.equal(taker.novigMarketMismatch(market, { eventId: NOVIG_EVENT, status: "OPEN_INGAME" }, row, contract, { live: true }), null);
  assert.match(taker.novigMarketMismatch(market, { eventId: NOVIG_EVENT, status: "OPEN_PREGAME" }, row, contract, { live: true }), /not in game/);
  assert.match(taker.novigMarketMismatch({ ...market, fee: { charged: "WHEN_LIVE", coefficient: "?" } }, { eventId: NOVIG_EVENT, status: "OPEN_INGAME" }, row, contract, { live: true }), /fee/);
  assert.match(taker.novigMarketMismatch(market, { eventId: NOVIG_EVENT, status: "OPEN_INGAME" }, row, contract), /event status "OPEN_INGAME"/);
});
