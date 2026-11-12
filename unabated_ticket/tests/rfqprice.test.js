// server/rfqprice.js: a Kalshi RFQ leg priced off Unabated's fair ladder (#282).
// Boards: fixtures/v2_venue_ids_slice.json (real CFB event 123742, Duquesne @
// Washington State, with Kalshi rung ids; fairs set on its lines here), and a
// small hand-built game per league for each league's ticker shapes.
"use strict";

const test = require("node:test");
const assert = require("node:assert/strict");
const fs = require("node:fs");
const path = require("node:path");
const feed = require("../extension/feed.js");
const tailflex = require("../extension/tailflex.js");
const edgeRows = require("../extension/edgerows.js");
const taker = require("../server/taker.js");
const rfqprice = require("../server/rfqprice.js");
const { moveLine } = require("./helpers.js");

const CFB_START_MS = 1789846200000;
const NOW = CFB_START_MS - 24 * 3600 * 1000;
const CLOSE = 1e-6;

function near(actual, expected, message) {
  assert.ok(Math.abs(actual - expected) < CLOSE, `${message || ""} expected ${expected}, got ${actual}`);
}

// A Kalshi book with this YES bid and NO bid (dollars).
function bookOf(yesBid, noBid) {
  return taker.parseKalshiOrderbook({ orderbook_fp: { yes_dollars: yesBid == null ? [] : [[String(yesBid), "100"]], no_dollars: noBid == null ? [] : [[String(noBid), "100"]] } });
}

function contextOf(board, extra = {}) {
  return {
    board, now: NOW, measurement: null, ladderOf: edgeRows.createLadderReaders(board),
    staleReasonOf: () => null, fairAsOf: () => NOW - 30 * 1000, ...extra,
  };
}

function kalshiOf(marketTicker, floorStrike, fields = {}) {
  return { market: { ticker: marketTicker, status: "active", floor_strike: floorStrike }, book: bookOf(0.5, 0.45), startMs: null, ...fields };
}

function price(board, marketTicker, side, { context = {}, kalshi = {}, floorStrike = null } = {}) {
  const leg = { marketTicker, side };
  return rfqprice.priceLeg(leg, rfqprice.resolveLeg(marketTicker, board), contextOf(board, context), { ...kalshiOf(marketTicker, floorStrike), ...kalshi });
}

// ---- the real CFB game ---------------------------------------------------------

const WSU_MAIN_KEY = "420942308:ms105:si1:tid717";
const WSU_ALT_355 = "420942308:ms105:si1:tid717:alt-35.5";
const WSU_ALT_415 = "420942308:ms105:si1:tid717:alt-41.5";
const ML_WSU_KEY = "411181220:ms105:si1:tid717";
const OVER_ALT_535 = "433644484:ms105:si0:tid927:alt53.5";

// Fairs on WSU's spread ladder: -41.5 at +150 (WSU covers 40%), the -38.5
// main at -110, -35.5 at -150 (60%); WSU's moneyline at -2000; Over 53.5 at
// +110. The fixture already holds Over/Under 52.5 (P(Over) 0.5025).
function cfbBoard() {
  const json = JSON.parse(fs.readFileSync(path.join(__dirname, "fixtures", "v2_venue_ids_slice.json"), "utf8"));
  const state = feed.parseSnapshot(json, { leagueId: 2 });
  moveLine(state, WSU_ALT_415, { bacr: 150 });
  moveLine(state, WSU_MAIN_KEY, { bacr: -110 });
  moveLine(state, WSU_ALT_355, { bacr: -150 });
  moveLine(state, ML_WSU_KEY, { bacr: -2000 });
  moveLine(state, OVER_ALT_535, { bacr: 110 });
  return state;
}

test("legs param: TICKER:side in order; anything malformed says why (the runner's 400), never a shorter answer", () => {
  assert.deepEqual(rfqprice.parseLegsParam("KXNFLSPREAD-26NOV09KCDEN-KC7:yes,KXNFLTOTAL-26NOV09KCDEN-45:NO"),
    [{ marketTicker: "KXNFLSPREAD-26NOV09KCDEN-KC7", side: "yes" }, { marketTicker: "KXNFLTOTAL-26NOV09KCDEN-45", side: "no" }]);
  for (const bad of [null, "", "KXNFLSPREAD-26NOV09KCDEN-KC7", "KXNFLSPREAD-26NOV09KCDEN-KC7:maybe", "KXNFLSPREAD-26NOV09KCDEN-KC7:yes,", "bad ticker:yes"]) {
    assert.equal(typeof rfqprice.parseLegsParam(bad), "string", `expected a refusal for ${bad}`);
  }
  const many = Array.from({ length: rfqprice.MAX_LEGS + 1 }, () => "KXNFLGAME-26NOV09KCDEN-KC:yes").join(",");
  assert.match(rfqprice.parseLegsParam(many), /at most/);
});

test("ticker shapes: each league's full-game moneyline, spread and total; anything else is unsupported_series", () => {
  assert.deepEqual(rfqprice.parseTicker("KXNFLSPREAD-26NOV09KCDEN-KC7"),
    { series: "KXNFLSPREAD", leagueId: 1, betTypeId: 2, eventTicker: "KXNFLSPREAD-26NOV09KCDEN", eventSuffix: "26NOV09KCDEN", strikeCode: "KC", number: 6.5 });
  assert.deepEqual(rfqprice.parseTicker("KXNBATOTAL-26NOV09BOSNYK-221"),
    { series: "KXNBATOTAL", leagueId: 3, betTypeId: 3, eventTicker: "KXNBATOTAL-26NOV09BOSNYK", eventSuffix: "26NOV09BOSNYK", strikeCode: null, number: 220.5 });
  assert.equal(rfqprice.parseTicker("KXNCAAMBGAME-26NOV10DUKEUNC-DUKE").strikeCode, "DUKE");
  assert.equal(rfqprice.parseTicker("KXNCAAMBGAME-26NOV10DUKEUNC-DUKE").leagueId, 4);
  assert.equal(rfqprice.parseTicker("KXMLBSPREAD-26APR091905NYYBOS-NYY2").number, 1.5);
  assert.equal(rfqprice.parseTicker("KXNCAAFGAME-26SEP19DUQWSU-WSU").leagueId, 2);
  for (const ticker of ["KXNFL1HSPREAD-26NOV09KCDEN-KC4", "KXMLBF5TOTAL-26APR091905NYYBOS-5", "KXMVECROSSCATEGORY-S2026ABC-XYZ", "KXNFLGAME-26NOV09KCDEN", "KXNFLTOTAL-26NOV09KCDEN-KC", "KXNFLTOTAL-26NOV09KCDEN-0", "KXNFLSPREAD-26NOV09KCDEN-45", "KXNFLSPREAD-26NOV09KCDEN-KC0", "KXMLBRFI-26APR091905NYYBOS-1"]) {
    assert.equal(rfqprice.parseTicker(ticker).reason, "unsupported_series", ticker);
  }
});

test("CFB alt rung: WSU wins by over 35.5 is Unabated's -35.5 fair; NO is its complement; Kalshi's book rides along", () => {
  const board = cfbBoard();
  const ticker = "KXNCAAFSPREAD-26SEP19DUQWSU-WSU36";
  const yes = price(board, ticker, "yes", { floorStrike: 35.5, kalshi: { book: bookOf(0.58, 0.40) } });
  // The whole leg, so a field that drifts from the quoter's contract (#283) fails here.
  assert.deepEqual(yes, {
    market_ticker: ticker, side: "yes", ok: true, reason: null, detail: null,
    fair: 0.6, fair_lo: 0.6, fair_hi: 0.6, bracketed: false,
    // 3 points off the -38.5 main, measured from even money (no measurement passed): tail flex at the fallback c.
    sigma: Math.round(tailflex.edgeUncertainty({ decimal: 1, fairProb: 0.6, zMain: 0, c: tailflex.C_FALLBACK }) * 1e6) / 1e6,
    c: tailflex.C_FALLBACK, tau: tailflex.TAU_FALLBACK,
    line: { league: "CFB", market: "spread", period: "FG", team: "Washington State", over_under: null, number: -35.5, unabated_event_id: 123742, event_name: "Duquesne Dukes - DUQ @ Washington State Cougars - WST" },
    kickoff_ms: CFB_START_MS, fair_age_ms: 30 * 1000,
    kalshi_status: "active", kalshi_yes_bid: 0.58, kalshi_yes_ask: 0.6, kalshi_no_ask: 0.42,
  });
  const no = price(board, ticker, "no", { floorStrike: 35.5 });
  near(no.fair, 0.4);
  assert.deepEqual([no.line.team, no.line.number], ["Duquesne", 35.5]);
});

test("CFB derived main number: WSU39 (no Kalshi rung) prices at Unabated's -38.5 main", () => {
  const board = cfbBoard();
  const leg = price(board, "KXNCAAFSPREAD-26SEP19DUQWSU-WSU39", "yes", { floorStrike: 38.5 });
  assert.equal(leg.ok, true, leg.detail);
  near(leg.fair, 110 / 210);
  assert.equal(leg.bracketed, false);
});

test("CFB between rungs: WSU38 is bracketed by -38.5 and -35.5, fair the midpoint; past the ladder it is declined", () => {
  const board = cfbBoard();
  const leg = price(board, "KXNCAAFSPREAD-26SEP19DUQWSU-WSU38", "yes", { floorStrike: 37.5 });
  assert.equal(leg.ok, true, leg.detail);
  assert.equal(leg.bracketed, true);
  near(leg.fair_lo, 110 / 210);
  near(leg.fair_hi, 0.6);
  near(leg.fair, (110 / 210 + 0.6) / 2);
  const past = price(board, "KXNCAAFSPREAD-26SEP19DUQWSU-WSU50", "yes", { floorStrike: 49.5 });
  assert.deepEqual([past.ok, past.reason], [false, "past_ladder"]);
  assert.match(past.detail, /outside the margin ladder/);
});

test("CFB moneyline: YES on WSU is Unabated's -2000; NO is the complement (overtime: a CFB game cannot end level)", () => {
  const board = cfbBoard();
  const yes = price(board, "KXNCAAFGAME-26SEP19DUQWSU-WSU", "yes");
  assert.equal(yes.ok, true, yes.detail);
  near(yes.fair, 2000 / 2100);
  assert.deepEqual([yes.sigma, yes.line.market, yes.line.team, yes.line.number], [0, "ml", "Washington State", null]);
  near(price(board, "KXNCAAFGAME-26SEP19DUQWSU-WSU", "no").fair, 100 / 2100);
});

test("CFB totals: Over 53.5 alt and the 52.5 main; sigma is zero at Unabated's main and grows off it", () => {
  const board = cfbBoard();
  const mainFair = 0.5024875621890548;
  const measurement = tailflex.withPooledTau({ markets: {}, mainFairs: new Map([["123742:1:3:0", { points: 52.5, fairProb: mainFair }]]), tauSamples: [] }, {});
  const main = price(board, "KXNCAAFTOTAL-26SEP19DUQWSU-53", "yes", { floorStrike: 52.5, context: { measurement } });
  assert.equal(main.ok, true, main.detail);
  near(main.fair, mainFair);
  assert.deepEqual([main.sigma, main.c, main.tau, main.line.over_under, main.line.number], [0, tailflex.C_FALLBACK, tailflex.TAU_FALLBACK, "over", 52.5]);
  const alt = price(board, "KXNCAAFTOTAL-26SEP19DUQWSU-54", "yes", { floorStrike: 53.5, context: { measurement } });
  near(alt.fair, 100 / 210);
  const z = tailflex.qnorm(100 / 210);
  near(alt.sigma, tailflex.dnorm(z) * tailflex.C_FALLBACK * Math.abs(z - tailflex.qnorm(mainFair)));
  const under = price(board, "KXNCAAFTOTAL-26SEP19DUQWSU-54", "no", { floorStrike: 53.5, context: { measurement } });
  near(under.fair, 110 / 210);
  assert.equal(under.line.over_under, "under");
});

test("declines: stale feed, no game, two games, started, Kalshi unread, market closed or at another strike", () => {
  const board = cfbBoard();
  const ticker = "KXNCAAFSPREAD-26SEP19DUQWSU-WSU36";
  assert.equal(price(board, ticker, "yes", { floorStrike: 35.5, context: { staleReasonOf: () => "league 2 snapshot landed 400 s ago" } }).reason, "stale_feed");
  const noEvent = price(board, "KXNCAAFSPREAD-26SEP20AAABBB-AAA4", "yes", { floorStrike: 3.5 });
  assert.deepEqual([noEvent.reason, noEvent.fair, noEvent.line, noEvent.fair_age_ms], ["no_event", null,
    { league: "CFB", market: "spread", period: "FG", team: null, over_under: null, number: null, unabated_event_id: null, event_name: null }, 30 * 1000]);
  assert.equal(price(board, "KXNFLSPREAD-26SEP19DUQWSU-WSU36", "yes", { floorStrike: 35.5 }).reason, "no_event", "a suffix matches only in its own league");
  const twice = cfbBoard();
  twice.events[999] = { ...twice.events[123742], eventId: 999 };
  assert.equal(price(twice, ticker, "yes", { floorStrike: 35.5 }).reason, "ambiguous_event");
  const started = price(board, ticker, "yes", { floorStrike: 35.5, context: { now: CFB_START_MS } });
  assert.equal(started.reason, "started");
  // Kalshi's own start, when earlier, is the kickoff.
  const early = price(board, ticker, "yes", { floorStrike: 35.5, kalshi: { startMs: NOW - 1 } });
  assert.deepEqual([early.reason, early.kickoff_ms], ["started", NOW - 1]);
  assert.equal(price(board, ticker, "yes", { kalshi: { error: "Kalshi 503" } }).reason, "kalshi_unread");
  // No start anywhere: the game may be under way, so it is never priced as pregame.
  const noStart = cfbBoard();
  noStart.events[123742].eventStart = null;
  assert.equal(price(noStart, ticker, "yes", { floorStrike: 35.5 }).reason, "no_start");
  assert.equal(price(noStart, ticker, "yes", { floorStrike: 35.5, kalshi: { startMs: CFB_START_MS } }).ok, true);
  const stale = price(board, ticker, "yes", { floorStrike: 35.5, context: { staleReasonOf: () => "stale build" } });
  assert.deepEqual([stale.reason, stale.line.league, stale.fair], ["stale_feed", "CFB", null]);
  assert.equal(rfqprice.needsKalshi(rfqprice.resolveLeg(ticker, board), contextOf(board, { staleReasonOf: () => "stale build" })), false, "a stale league reads no Kalshi");
  assert.equal(rfqprice.needsKalshi(rfqprice.resolveLeg(ticker, board), contextOf(board)), true);
  assert.equal(price(board, ticker, "yes", { floorStrike: 36.5 }).reason, "market_mismatch");
  const closed = price(board, ticker, "yes", { kalshi: { market: { ticker, status: "closed", floor_strike: 35.5 } } });
  assert.deepEqual([closed.reason, closed.kalshi_status], ["market_mismatch", "closed"]);
  assert.equal(rfqprice.priceLeg({ marketTicker: ticker, side: "yes" }, rfqprice.resolveLeg(ticker, null), contextOf(board), {}).reason, "no_board");
});

test("no_team: Unabated's rung and the ticker disagree on the YES team, or the rung sits at another number", () => {
  const board = cfbBoard();
  const contracts = board.events[123742].venueIds.kalshiContracts;
  const ticker = "KXNCAAFSPREAD-26SEP19DUQWSU-WSU36";
  contracts[`Y-${ticker}`] = { ...contracts[`Y-${ticker}`], sideIndex: 0 };
  contracts[`N-${ticker}`] = { ...contracts[`N-${ticker}`], sideIndex: 1 };
  assert.equal(price(board, ticker, "yes", { floorStrike: 35.5 }).reason, "no_team");
  const moved = cfbBoard();
  const movedContracts = moved.events[123742].venueIds.kalshiContracts;
  movedContracts[`Y-${ticker}`] = { ...movedContracts[`Y-${ticker}`], points: -34.5 };
  assert.match(price(moved, ticker, "yes", { floorStrike: 35.5 }).detail, /sits at 34.5/);
});

test("a flat-lined rung at the cut is not a fair: flat_tail", () => {
  const board = cfbBoard();
  moveLine(board, WSU_ALT_415, { bacr: -110 }); // the -38.5 main's fair again
  assert.equal(price(board, "KXNCAAFSPREAD-26SEP19DUQWSU-WSU42", "yes", { floorStrike: 41.5 }).reason, "flat_tail");
});

// ---- one hand-built game per league ---------------------------------------------

// AWY @ HOM. Home spread fairs: -6.5 at +150 (40%), -3.5 at -110; moneylines
// AWY +150 / HOM -150; Over 44.5 at -120, Over 47.5 at +130. Kalshi rungs:
// HOM7 (the -6.5 alt) and the total's 48 (Over 47.5).
function leagueBoard(leagueId, series, suffix) {
  const eventId = 500 + leagueId;
  const line = (key, betTypeId, sideIndex, points, bacr, isAlt = false) => [key, { key, eventId, leagueId, periodTypeId: 1, betTypeId, sideIndex, points, bacr, bookId: 105, isAlt }];
  const lines = Object.fromEntries([
    line("hs", 2, 1, -3.5, -110), line("hs-alt", 2, 1, -6.5, 150, true),
    line("ml-a", 1, 0, null, 150), line("ml-h", 1, 1, null, -150),
    line("ov", 3, 0, 44.5, -120), line("ov-alt", 3, 0, 47.5, 130, true),
  ]);
  const contracts = {
    [`Y-${series}SPREAD-${suffix}-HOM7`]: { lineKey: "hs-alt", mainKey: "hs", points: -6.5, sideIndex: 1 },
    [`Y-${series}TOTAL-${suffix}-48`]: { lineKey: "ov-alt", mainKey: "ov", points: 47.5, sideIndex: 0 },
  };
  return {
    teams: { 1: "Away Team", 2: "Home Team" }, lines,
    events: { [eventId]: { eventId, leagueId, eventName: "AWY @ HOM", eventStart: CFB_START_MS, awayTeamId: 1, homeTeamId: 2, venueIds: { kalshiEventSuffixes: [suffix], kalshiContracts: contracts, novigOutcomes: {} } } },
  };
}

const LEAGUE_SHAPES = [
  { leagueId: 1, label: "NFL", series: "KXNFL", suffix: "26NOV09AWYHOM" },
  { leagueId: 2, label: "CFB", series: "KXNCAAF", suffix: "26NOV14AWYHOM" },
  { leagueId: 3, label: "NBA", series: "KXNBA", suffix: "26NOV09AWYHOM" },
  { leagueId: 4, label: "CBB", series: "KXNCAAMB", suffix: "26NOV10AWYHOM" },
  { leagueId: 5, label: "MLB", series: "KXMLB", suffix: "26APR091905AWYHOM" },
];

for (const shape of LEAGUE_SHAPES) {
  test(`${shape.label}: main, alt and derived main number spreads, both moneylines, totals`, () => {
    const board = leagueBoard(shape.leagueId, shape.series, shape.suffix);
    const spread = (strike) => `${shape.series}SPREAD-${shape.suffix}-${strike}`;
    const cases = [
      // [ticker, side, floor_strike, fair, number, team / over_under]
      [spread("HOM7"), "yes", 6.5, 0.4, -6.5, "Home Team"], // alt rung
      [spread("HOM4"), "yes", 3.5, 110 / 210, -3.5, "Home Team"], // main number, no Kalshi rung: derived
      [spread("HOM4"), "no", 3.5, 100 / 210, 3.5, "Away Team"],
      [spread("HOM5"), "yes", 4.5, (0.4 + 110 / 210) / 2, -4.5, "Home Team"], // between rungs
      [`${shape.series}GAME-${shape.suffix}-HOM`, "yes", null, 0.6, null, "Home Team"],
      [`${shape.series}GAME-${shape.suffix}-AWY`, "yes", null, 0.4, null, "Away Team"],
      [`${shape.series}TOTAL-${shape.suffix}-45`, "yes", 44.5, 120 / 220, 44.5, "over"],
      [`${shape.series}TOTAL-${shape.suffix}-48`, "no", 47.5, 130 / 230, 47.5, "under"],
    ];
    for (const [ticker, side, floorStrike, fair, number, who] of cases) {
      const leg = price(board, ticker, side, { floorStrike });
      assert.equal(leg.ok, true, `${ticker} ${side}: ${leg.reason} ${leg.detail}`);
      if (shape.leagueId === 1 && ticker.includes("GAME")) {
        // An NFL game can end level: YES is bracketed down by the tie ceiling, NO up.
        near(leg.fair_hi, fair, `${ticker} ${side}`);
        near(leg.fair_lo, fair * (1 - rfqprice.NFL_TIE_PROB_MAX), `${ticker} ${side}`);
        assert.equal(leg.bracketed, true);
      } else {
        near(leg.fair, fair, `${ticker} ${side}`);
      }
      assert.equal(leg.line.league, shape.label);
      assert.equal(leg.line.number, number, ticker);
      assert.equal(leg.line.team ?? leg.line.over_under, who, ticker);
    }
    assert.equal(price(board, spread("HOM21"), "yes", { floorStrike: 20.5 }).reason, "past_ladder");
    assert.equal(price(board, spread("XYZ4"), "yes", { floorStrike: 3.5 }).reason, "no_team");
    if (shape.leagueId === 1) {
      const no = price(board, `${shape.series}GAME-${shape.suffix}-HOM`, "no");
      near(no.fair_lo, 0.4);
      near(no.fair_hi, 1 - 0.6 * (1 - rfqprice.NFL_TIE_PROB_MAX));
    }
  });
}
