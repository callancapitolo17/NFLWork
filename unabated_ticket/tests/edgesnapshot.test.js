// Run: node --test unabated_ticket/tests
// server/edgesnapshot.js: the market at a bet's fill and close. Board:
// Seahawks @ 49ers, the 49ers -14.5 alt and the Seahawks +14.5 alt on the
// other side, Unabated's main 49ers -8.5; a Novig bet at +208 placed at
// T0+2min. The close runs end to end in runner.test.js against the NFL slice.
const test = require("node:test");
const assert = require("node:assert/strict");
const feed = require("../extension/feed.js");
const edgemove = require("../extension/edgemove.js");
const fillbooks = require("../server/fillbooks.js");
const edgesnapshot = require("../server/edgesnapshot.js");

const MIN = 60 * 1000;
const T0 = Date.UTC(2026, 9, 22, 17, 0, 0);
const KICKOFF = T0 + 3 * 60 * MIN;
const WATCHING_SINCE = T0 - 60 * MIN;
const iso = (ms) => new Date(ms).toISOString();

const EVENT_ID = 900;
const UNABATED = feed.UNABATED_LINE_BOOK_ID;
const NOVIG = 89;
const KALSHI = 105;
const DRAFTKINGS = 4;
const NINERS_ALT = (book) => `m1:ms${book}:si1:tid22:alt-14.5`;
const SEAHAWKS_ALT = (book) => `m2:ms${book}:si0:tid11:alt14.5`;
const UNABATED_MAIN = `m1:ms${UNABATED}:si1:tid22`;

function line(bookId, overrides = {}) {
  const mainKey = `m1:ms${bookId}:si1:tid22`;
  return {
    key: NINERS_ALT(bookId), isAlt: true, mainKey, mainPoints: -8.5,
    leagueId: 1, periodTypeId: 1, betTypeId: 2, eventId: EVENT_ID, marketId: "m1", bookId, sideKey: "si1:tid22", sideIndex: 1,
    points: -14.5, price: 223, sourceFormat: 1, sourcePrice: null, bacr: 188, ge: 0.12, liquidity: null, statusId: 1,
    sequenceNumber: T0, isBlurred: false, modifiedOn: null,
    ...overrides,
  };
}

function otherSide(bookId, overrides = {}) {
  return line(bookId, {
    key: SEAHAWKS_ALT(bookId), mainKey: `m2:ms${bookId}:si0:tid11`, mainPoints: 8.5, marketId: "m2",
    sideKey: "si0:tid11", sideIndex: 0, points: 14.5, price: -260, bacr: -188, ...overrides,
  });
}

function unabatedMain(overrides = {}) {
  return line(UNABATED, { key: UNABATED_MAIN, isAlt: false, mainKey: undefined, mainPoints: undefined, points: -8.5, price: -112, bacr: -112, ...overrides });
}

function boardState(lines) {
  const state = feed.emptyState();
  state.leagues = [1];
  state.teams = { 11: "Seattle Seahawks", 22: "San Francisco 49ers" };
  const book = (id, name, hasLiquidity = false) => ({ id, name, isLive: true, hasLiquidity });
  state.books = {
    [UNABATED]: book(UNABATED, "Unabated"), [NOVIG]: book(NOVIG, "Novig", true), [KALSHI]: book(KALSHI, "Kalshi", true),
    [DRAFTKINGS]: book(DRAFTKINGS, "DraftKings"),
  };
  state.events[EVENT_ID] = {
    eventId: EVENT_ID, leagueId: 1, eventName: "Seahawks @ 49ers", eventStart: KICKOFF,
    awayTeamId: 11, homeTeamId: 22, awayRotation: 101, homeRotation: 102, venueIds: null,
  };
  const main = line(DRAFTKINGS, { isAlt: false, key: `m1:ms${DRAFTKINGS}:si1:tid22`, points: -8.5, price: -110, bacr: -112 });
  for (const each of [main, ...lines]) state.lines[each.key] = each;
  return state;
}

function boardLinesOf(state) {
  const seen = new Set();
  const rows = [];
  for (const each of Object.values(state.lines)) {
    if (each.isAlt || seen.has(each.eventId) || each.bookId === UNABATED) continue;
    seen.add(each.eventId);
    rows.push(feed.describeLine(each, state));
  }
  return rows;
}

function bet(id, overrides = {}) {
  return {
    id, source: "novig_api", venue: "novig", league: "nfl", status: "open",
    eventStart: iso(KICKOFF), eventDate: null, awayTeam: "Seattle Seahawks", homeTeam: "San Francisco 49ers",
    awayKey: "nfl:11", homeKey: "nfl:22", betType: "spread", period: "FG", side: "home", points: -14.5,
    price: 208, stake: 150, toWin: 312, placedAt: iso(T0 + 2 * MIN), closedAt: null, unmatchable: null, approx: [],
    isParlayLeg: false, raw: {},
    ...overrides,
  };
}

// observations: [line, [{at, ...fields}]] pairs, through edgemove.observe.
function historyOf(observations) {
  const history = {};
  for (const [base, entries] of observations) {
    for (const { at, ...fields } of entries) edgemove.observe(history, { ...base, ...fields }, { at });
  }
  return history;
}

function capture({ records, lines, history, observingSince = WATCHING_SINCE, leagueObservingSince = { 1: WATCHING_SINCE }, now = T0 + 6 * MIN }) {
  const state = boardState(lines);
  return edgesnapshot.captureFillSnapshots({
    records, skipIds: new Set(), state, boardLines: boardLinesOf(state), history, observingSince, leagueObservingSince,
    linesOfEvent: (eventId) => Object.values(state.lines).filter((each) => each.eventId === eventId), now,
  });
}

const price = (value, sourcePrice = null, sourceFormat = 1) => ({ price: value, sourcePrice, sourceFormat, liquidity: null });

test("the fill: every book's two sides as they stood at the fill, Unabated's main number and fair, the venue's edge-move tag", () => {
  const lines = [line(NOVIG), otherSide(NOVIG), line(KALSHI), otherSide(KALSHI), line(DRAFTKINGS), unabatedMain()];
  const history = historyOf([
    // Novig lengthened its price with the fair unchanged inside the ten minutes before the fill: book moved away.
    [line(NOVIG), [{ at: T0 - 5 * MIN, price: 190, sourcePrice: 0.345, sourceFormat: 4 }, { at: T0 + MIN, price: 208, sourcePrice: 0.325, sourceFormat: 4 }]],
    [otherSide(NOVIG), [{ at: T0 - 5 * MIN, price: -240, sourcePrice: 0.7, sourceFormat: 4 }]],
    // Kalshi moved after the fill: the fill-time price is kept.
    [line(KALSHI), [{ at: T0, price: 215 }, { at: T0 + 4 * MIN, price: 190 }]],
    [otherSide(KALSHI), [{ at: T0, price: -250 }]],
    // DraftKings was first seen after the fill: left out.
    [line(DRAFTKINGS), [{ at: T0 + 3 * MIN, price: 200 }]],
    [unabatedMain(), [{ at: T0 - 30 * MIN, price: -112, bacr: -112 }, { at: T0 + 5 * MIN, price: -120, bacr: -120 }]],
  ]);
  const result = capture({ records: [bet("novig:1")], lines, history });
  assert.deepEqual(result.refusals, []);
  const [row] = result.saves;
  assert.deepEqual({ ...row, edgeMoveFairDelta: undefined, edgeMovePriceDelta: undefined }, {
    betId: "novig:1", placedAt: iso(T0 + 2 * MIN), lineKey: NINERS_ALT(NOVIG), points: -14.5,
    mainPoints: -8.5, mainFairAmerican: -112, pointsFromMain: -6,
    edgeMoveLineKey: NINERS_ALT(NOVIG), edgeMove: "book_away", edgeMoveFairDelta: undefined, edgeMovePriceDelta: undefined,
    books: [
      { bookId: NOVIG, book: "Novig", exchange: true, venue: true, side: price(208, 0.325, 4), other: price(-240, 0.7, 4) },
      { bookId: KALSHI, book: "Kalshi", exchange: true, venue: false, side: price(215), other: price(-250) },
    ],
  });
  assert.equal(row.edgeMoveFairDelta, 0);
  assert.ok(Math.abs(row.edgeMovePriceDelta - 0.02) < 1e-9);
});

test("a book with no line on the other side keeps a null other; Unabated's main moved since the fill leaves the main fields null", () => {
  const lines = [line(KALSHI), unabatedMain({ points: -9.5 })];
  const history = historyOf([
    [line(KALSHI), [{ at: T0, price: 215 }]],
    // The history restarts on a new number, so the fill-time main is gone.
    [unabatedMain(), [{ at: T0, price: -112 }]],
    [unabatedMain({ points: -9.5 }), [{ at: T0 + 4 * MIN, price: -105 }]],
  ]);
  const [row] = capture({ records: [bet("novig:1")], lines, history }).saves;
  assert.deepEqual(row.books, [{ bookId: KALSHI, book: "Kalshi", exchange: true, venue: false, side: price(215), other: null }]);
  assert.deepEqual([row.mainPoints, row.mainFairAmerican, row.pointsFromMain, row.edgeMove], [null, null, null, "none"]);
});

test("same refusals as the fill books: never a backfill, never a parlay leg, never past the history window", () => {
  const history = historyOf([[line(KALSHI), [{ at: T0, price: 215 }]]]);
  const early = capture({ records: [bet("novig:1")], lines: [line(KALSHI)], history, observingSince: T0 + 3 * MIN });
  assert.deepEqual(early, { saves: [], refusals: [{ betId: "novig:1", reason: fillbooks.REFUSED.notWatching }] });
  const gap = capture({ records: [bet("novig:1")], lines: [line(KALSHI)], history, leagueObservingSince: { 1: T0 + 3 * MIN } });
  assert.deepEqual(gap, { saves: [], refusals: [{ betId: "novig:1", reason: fillbooks.REFUSED.leagueGap }] });
  const leg = capture({ records: [bet("bfa:1", { isParlayLeg: true })], lines: [line(KALSHI)], history });
  assert.deepEqual(leg.refusals, [{ betId: "bfa:1", reason: fillbooks.REFUSED.notStraight }]);
  const late = capture({ records: [bet("novig:1")], lines: [line(KALSHI)], history, now: T0 + 2 * MIN + 11 * MIN });
  assert.deepEqual(late.refusals, [{ betId: "novig:1", reason: fillbooks.REFUSED.tooLate }]);
});

test("the other side: a spread's number flips, a total's stays, a moneyline has none", () => {
  assert.equal(edgesnapshot.oppositePoints({ betTypeId: 2, points: -14.5 }), 14.5);
  assert.equal(edgesnapshot.oppositePoints({ betTypeId: 3, points: 47.5 }), 47.5);
  assert.equal(edgesnapshot.oppositePoints({ betTypeId: 1, points: null }), null);
  const over = { bookId: 4, periodTypeId: 1, betTypeId: 3, sideIndex: 0, points: 47.5 };
  const under = { ...over, sideIndex: 1 };
  assert.equal(edgesnapshot.otherSideOf(over, [over, { ...under, points: 48.5 }, under]), under);
  const away = { bookId: 4, periodTypeId: 1, betTypeId: 1, sideIndex: 0, points: null };
  const home = { ...away, sideIndex: 1 };
  const draw = { ...away, sideIndex: 2 };
  assert.equal(edgesnapshot.otherSideOf(away, [draw, { ...home, periodTypeId: 2 }, home]), home);
  assert.equal(edgesnapshot.otherSideOf(draw, [away, home]), null, "a three-way moneyline's draw has no other side");
});

test("close snapshots are resent only when their books change, not on every new league snapshot", () => {
  const row = { betId: "b1", observedAt: iso(T0), eventStart: iso(KICKOFF), books: [{ bookId: 105, side: price(215), other: null }] };
  const sent = new Map();
  assert.deepEqual(edgesnapshot.unsentCloseRows([row], sent), [row]);
  edgesnapshot.markCloseSent([row], sent);
  assert.deepEqual(edgesnapshot.unsentCloseRows([row], sent), []);
  assert.deepEqual(edgesnapshot.unsentCloseRows([{ ...row, observedAt: iso(T0 + MIN) }], sent), []);
  const moved = { ...row, books: [{ bookId: 105, side: price(220), other: null }] };
  assert.deepEqual(edgesnapshot.unsentCloseRows([moved], sent), [moved]);
});
