// Run: node --test unabated_ticket/tests
// server/fillbooks.js: the best other book's price at the fill and how many
// books beat the price taken. Board: Seahawks @ 49ers, the 49ers -14.5 alt,
// a Novig bet at +208 placed at T0+2min.
const test = require("node:test");
const assert = require("node:assert/strict");
const feed = require("../extension/feed.js");
const edgemove = require("../extension/edgemove.js");
const { captureFillBooks, REFUSED } = require("../server/fillbooks.js");

const MIN = 60 * 1000;
const T0 = Date.UTC(2026, 8, 23, 17, 0, 0);
const KICKOFF = T0 + 3 * 60 * MIN;
const WATCHING_SINCE = T0 - 60 * MIN;
const iso = (ms) => new Date(ms).toISOString();

const EVENT_ID = 900;
const NOVIG = 89;
const KALSHI = 105;
const DRAFTKINGS = 4;
const THESCORE_US = 36;
const THESCORE_CA = 60;
const NINERS_ALT = (book) => `m1:ms${book}:si1:tid22:alt-14.5`;

function line(bookId, overrides = {}) {
  const mainKey = `m1:ms${bookId}:si1:tid22`;
  return {
    key: `${mainKey}:alt-14.5`, isAlt: true, mainKey, mainPoints: -8.5,
    leagueId: 1, periodTypeId: 1, betTypeId: 2, eventId: EVENT_ID, marketId: "m1", bookId, sideKey: "si1:tid22", sideIndex: 1,
    points: -14.5, price: 223, sourceFormat: 1, sourcePrice: null, bacr: 188, ge: 0.12, liquidity: 500, statusId: 1,
    sequenceNumber: T0, isBlurred: false, modifiedOn: null,
    ...overrides,
  };
}

function boardState(lines) {
  const state = feed.emptyState();
  state.leagues = [1];
  state.teams = { 11: "Seattle Seahawks", 22: "San Francisco 49ers" };
  const book = (id, name, isLive = true) => ({ id, name, isLive, hasLiquidity: false });
  state.books = {
    [NOVIG]: book(NOVIG, "Novig"), [KALSHI]: book(KALSHI, "Kalshi"), [DRAFTKINGS]: book(DRAFTKINGS, "DraftKings"),
    [THESCORE_US]: book(THESCORE_US, "TheScore US"), [THESCORE_CA]: book(THESCORE_CA, "TheScore CA"),
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
    if (each.isAlt || seen.has(each.eventId)) continue;
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

// observations: line key -> [{at, price}], through edgemove.observe.
function historyOf(observations) {
  const history = {};
  for (const [key, entries] of Object.entries(observations)) {
    const [, bookId] = /:ms(\d+):/.exec(key);
    for (const { at, ...fields } of entries) edgemove.observe(history, { ...line(Number(bookId)), key, ...fields }, { at });
  }
  return history;
}

function capture({ records, lines, history, observingSince = WATCHING_SINCE, leagueObservingSince = { 1: WATCHING_SINCE },
  now = T0 + 6 * MIN }) {
  const state = boardState(lines);
  return captureFillBooks({ records, skipIds: new Set(), state, boardLines: boardLinesOf(state), history, observingSince, leagueObservingSince, now });
}

test("reads each other book's price at the fill, leaves the bet's own book out, counts the books that beat it", () => {
  const history = historyOf({
    [NINERS_ALT(NOVIG)]: [{ at: T0, price: 208 }],
    [NINERS_ALT(KALSHI)]: [{ at: T0, price: 215 }, { at: T0 + 5 * MIN, price: 190 }],
    [NINERS_ALT(DRAFTKINGS)]: [{ at: T0, price: 200 }, { at: T0 + 4 * MIN, price: 230 }],
  });
  const result = capture({ records: [bet("novig:1")], lines: [line(NOVIG), line(KALSHI), line(DRAFTKINGS)], history });
  assert.deepEqual(result.refusals, []);
  assert.deepEqual(result.saves, [{
    betId: "novig:1", placedAt: iso(T0 + 2 * MIN), fillPrice: 208, bookCount: 2, booksBetter: 1, bestBook: "Kalshi", bestPrice: 215,
  }]);
});

test("one book under two country names counts once at a shared price", () => {
  const history = historyOf({
    [NINERS_ALT(THESCORE_US)]: [{ at: T0, price: 220 }],
    [NINERS_ALT(THESCORE_CA)]: [{ at: T0, price: 220 }],
  });
  const [row] = capture({ records: [bet("novig:1")], lines: [line(THESCORE_US), line(THESCORE_CA)], history }).saves;
  assert.equal(row.bookCount, 1);
  assert.equal(row.booksBetter, 1);
  assert.equal(row.bestBook, "TheScore US / CA");
});

test("a book first seen after the fill, blurred, or at a dead book is left out", () => {
  const history = historyOf({
    [NINERS_ALT(KALSHI)]: [{ at: T0 + 3 * MIN, price: 300 }],
    [NINERS_ALT(DRAFTKINGS)]: [{ at: T0, price: 300 }],
    [NINERS_ALT(THESCORE_US)]: [{ at: T0, price: 300 }],
  });
  const lines = [line(KALSHI), line(DRAFTKINGS, { isBlurred: true }), line(THESCORE_US, { statusId: 2 })];
  const [row] = capture({ records: [bet("novig:1")], lines, history }).saves;
  assert.deepEqual({ bookCount: row.bookCount, booksBetter: row.booksBetter, bestBook: row.bestBook, bestPrice: row.bestPrice },
    { bookCount: 0, booksBetter: 0, bestBook: null, bestPrice: null });
});

test("never a backfill: before the scanner watched, or with the league not loading, is refused", () => {
  const history = historyOf({ [NINERS_ALT(KALSHI)]: [{ at: T0, price: 215 }] });
  const early = capture({ records: [bet("novig:1")], lines: [line(KALSHI)], history, observingSince: T0 + 3 * MIN });
  assert.deepEqual(early, { saves: [], refusals: [{ betId: "novig:1", reason: REFUSED.notWatching }] });
  const gap = capture({ records: [bet("novig:1")], lines: [line(KALSHI)], history, leagueObservingSince: { 1: T0 + 3 * MIN } });
  assert.deepEqual(gap, { saves: [], refusals: [{ betId: "novig:1", reason: REFUSED.leagueGap }] });
});

test("parlay legs and closed bets are skipped; a game the board has not resolved waits", () => {
  const history = historyOf({ [NINERS_ALT(KALSHI)]: [{ at: T0, price: 215 }] });
  const leg = capture({ records: [bet("bfa:1", { isParlayLeg: true })], lines: [line(KALSHI)], history });
  assert.deepEqual(leg.refusals, [{ betId: "bfa:1", reason: REFUSED.notStraight }]);
  const closed = capture({ records: [bet("novig:2", { status: "settled" })], lines: [line(KALSHI)], history });
  assert.deepEqual(closed, { saves: [], refusals: [] });
  const offBoard = capture({ records: [bet("novig:3", { points: -20.5 })], lines: [line(KALSHI)], history });
  assert.deepEqual(offBoard, { saves: [], refusals: [] });
});

test("a venue Unabated lists by name only is its own book too", () => {
  const history = historyOf({ [NINERS_ALT(DRAFTKINGS)]: [{ at: T0, price: 230 }], [NINERS_ALT(KALSHI)]: [{ at: T0, price: 200 }] });
  const [row] = capture({ records: [bet("dk:1", { venue: "draftkings" })], lines: [line(DRAFTKINGS), line(KALSHI)], history }).saves;
  assert.deepEqual({ bookCount: row.bookCount, bestBook: row.bestBook }, { bookCount: 1, bestBook: "Kalshi" });
});

test("a Kalshi fill is compared all-in of its taker fee, as Unabated shows Kalshi", () => {
  // 30c before the fee is 31.47c all-in (0.30 + 0.07 x 0.30 x 0.70), +217.76.
  // DraftKings +220 (31.25%) beats the all-in price but not the pre-fee 30%;
  // Novig +215 (31.75%) beats neither.
  const history = historyOf({
    [NINERS_ALT(DRAFTKINGS)]: [{ at: T0, price: 220 }],
    [NINERS_ALT(NOVIG)]: [{ at: T0, price: 215 }],
  });
  const kalshiBet = bet("kalshi:1", { venue: "kalshi", price: 233, raw: { vwapCents: 30 } });
  const [row] = capture({ records: [kalshiBet], lines: [line(DRAFTKINGS), line(NOVIG)], history }).saves;
  assert.equal(row.fillPrice, 217.76);
  assert.equal(row.booksBetter, 1);
});

test("a bet first read after the history window is refused, never counted from a partial history", () => {
  const history = historyOf({ [NINERS_ALT(KALSHI)]: [{ at: T0, price: 215 }] });
  const late = capture({ records: [bet("novig:1")], lines: [line(KALSHI)], history, now: T0 + 2 * MIN + 11 * MIN });
  assert.deepEqual(late, { saves: [], refusals: [{ betId: "novig:1", reason: REFUSED.tooLate }] });
});
