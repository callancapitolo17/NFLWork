// Run: node --test unabated_ticket/tests
// pricecheck.js: the other books' prices at the same number — the "N books
// better · would skip" and "outlier, check" tags (Cal's picks, 2026-10-09).
const test = require("node:test");
const assert = require("node:assert/strict");
const fs = require("node:fs");
const path = require("node:path");
const feed = require("../extension/feed.js");
const pricecheck = require("../extension/pricecheck.js");
const edgeRows = require("../extension/edgerows.js");

const NOW = Date.parse("2026-10-09T12:00:00Z");
const HOUR = 3600 * 1000;

// A hand-built feed state: one game, a total at 47.5 (side 0 = Over), a
// moneyline, and books keyed by id. Each line: [bookId, sideIndex, points, price, extra].
function stateWith(lines, books) {
  const state = feed.emptyState();
  state.books = books || {
    1: { id: 1, name: "BetMGM", isLive: true }, 2: { id: 2, name: "Circa", isLive: true },
    3: { id: 3, name: "Pinnacle", isLive: true }, 4: { id: 4, name: "DraftKings", isLive: true },
    5: { id: 5, name: "FanDuel", isLive: true }, 6: { id: 6, name: "Dead Book", isLive: false },
    49: { id: 49, name: "Unabated", isLive: true },
  };
  state.events[100] = { eventId: 100, leagueId: 1, eventStart: NOW + 10 * HOUR };
  for (const [bookId, sideIndex, points, price, extra] of lines) {
    const betTypeId = points == null ? 1 : 3;
    const key = `m${betTypeId}:ms${bookId}:si${sideIndex}:${points}`;
    state.lines[key] = {
      key, isAlt: false, leagueId: 1, periodTypeId: 1, betTypeId, eventId: 100, bookId, sideIndex, points, price,
      statusId: 1, isBlurred: false, modifiedOn: new Date(NOW - HOUR).toISOString(), ...(extra || {}),
    };
  }
  return state;
}

function overSpec(bookId, price, points) {
  return { eventId: 100, periodTypeId: 1, betTypeId: 3, sideIndex: 0, points: points ?? 47.5, bookId, price };
}

test("counts books that beat the price at the same number, best first; 2+ reads amber", () => {
  const state = stateWith([
    [1, 0, 47.5, -112], [2, 0, 47.5, -105], [3, 0, 47.5, -107], [4, 0, 47.5, -115], [5, 0, 47.5, -112],
  ]);
  const index = pricecheck.buildPriceIndex(state, { now: NOW, maxLineAgeMs: 168 * HOUR });
  const check = pricecheck.comparePrice(index, overSpec(1, -112));
  assert.deepEqual(check.better, [{ bookName: "Circa", price: -105, liquidity: null }, { bookName: "Pinnacle", price: -107, liquidity: null }]);
  assert.deepEqual(check.worse, [{ bookName: "DraftKings", price: -115, liquidity: null }]);
  assert.equal(check.same, 1);
  assert.equal(check.bookCount, 5);
  assert.equal(check.wouldSkip, true);
  assert.equal(check.outlier, null);
  assert.deepEqual(pricecheck.priceCheckTag(check), { kind: "skip", label: "2 books better" });
});

test("one book better tags but does not skip; best price says so", () => {
  const state = stateWith([[1, 0, 47.5, -110], [2, 0, 47.5, -105], [3, 0, 47.5, -115]]);
  const index = pricecheck.buildPriceIndex(state, { now: NOW, maxLineAgeMs: null });
  const one = pricecheck.comparePrice(index, overSpec(1, -110));
  assert.equal(one.wouldSkip, false);
  assert.deepEqual(pricecheck.priceCheckTag(one), { kind: "better", label: "1 book better" });
  assert.deepEqual(pricecheck.priceCheckTag(pricecheck.comparePrice(index, overSpec(2, -105))), { kind: "best", label: "best of 3 books" });
});

test("only the same number and side count; another number or the other side is not compared", () => {
  const state = stateWith([[1, 0, 47.5, -110], [2, 0, 48.5, +120], [3, 1, 47.5, +150]]);
  const index = pricecheck.buildPriceIndex(state, { now: NOW, maxLineAgeMs: null });
  const check = pricecheck.comparePrice(index, overSpec(1, -110));
  assert.equal(check.bookCount, 1);
  assert.equal(pricecheck.priceCheckTag(check), null);
});

test("stale, dead, blurred, off-board and Unabated's own lines never count as a better price", () => {
  const state = stateWith([
    [1, 0, 47.5, -110],
    [2, 0, 47.5, +200, { modifiedOn: new Date(NOW - 200 * HOUR).toISOString() }],
    [6, 0, 47.5, +200],
    [3, 0, 47.5, +200, { isBlurred: true }],
    [4, 0, 47.5, +200, { statusId: 2 }],
    [49, 0, 47.5, +200],
  ]);
  const index = pricecheck.buildPriceIndex(state, { now: NOW, maxLineAgeMs: 168 * HOUR });
  assert.equal(pricecheck.comparePrice(index, overSpec(1, -110)).bookCount, 1);
});

test("a book on the number twice (main and alt) counts once at its better price", () => {
  const state = stateWith([[1, 0, 47.5, -110], [2, 0, 47.5, -115], [3, 0, 47.5, -120]]);
  state.lines["alt"] = { ...state.lines["m3:ms2:si0:47.5"], key: "alt", isAlt: true, price: -105, sequenceNumber: NOW - HOUR, modifiedOn: "0001-01-01T00:00:00" };
  const index = pricecheck.buildPriceIndex(state, { now: NOW, maxLineAgeMs: 168 * HOUR });
  const check = pricecheck.comparePrice(index, overSpec(1, -110));
  assert.deepEqual(check.better, [{ bookName: "Circa", price: -105, liquidity: null }]);
  assert.equal(check.bookCount, 3);
});

test("one book under two country names counts once at a shared price; a different price still counts", () => {
  const books = {
    1: { id: 1, name: "Circa", isLive: true }, 36: { id: 36, name: "TheScore US", isLive: true },
    60: { id: 60, name: "TheScore CA", isLive: true }, 107: { id: 107, name: "Polymarket", isLive: true },
    111: { id: 111, name: "Polymarket US", isLive: true },
  };
  const state = stateWith([[1, 0, 47.5, -112], [36, 0, 47.5, -105], [60, 0, 47.5, -105], [107, 0, 47.5, -104], [111, 0, 47.5, -103]], books);
  const index = pricecheck.buildPriceIndex(state, { now: NOW, maxLineAgeMs: null });
  const check = pricecheck.comparePrice(index, { ...overSpec(1, -112), bookName: "Circa" });
  assert.deepEqual(check.better, [
    { bookName: "Polymarket US", price: -103, liquidity: null }, { bookName: "Polymarket", price: -104, liquidity: null }, { bookName: "TheScore US / CA", price: -105, liquidity: null },
  ]);
  assert.equal(check.bookCount, 4);
  // Betting theScore US: its CA copy is you, not another book.
  const own = pricecheck.comparePrice(index, { ...overSpec(36, -105), bookName: "TheScore US" });
  assert.equal(own.bookCount, 4);
  assert.ok(!own.better.concat(own.worse).some((entry) => entry.bookName.startsWith("TheScore")));
});

test("a book under $30 resting at its price is thin: listed apart, never better, counted or an outlier", () => {
  const state = stateWith([
    [1, 0, 47.5, -110], [2, 0, 47.5, -105, { liquidity: 12 }], [3, 0, 47.5, -107, { liquidity: 30 }], [4, 0, 47.5, -115],
  ]);
  const index = pricecheck.buildPriceIndex(state, { now: NOW, maxLineAgeMs: null });
  const check = pricecheck.comparePrice(index, overSpec(1, -110));
  // $30 exactly counts; $12 does not, so one book better and no Skip.
  assert.deepEqual(check.better, [{ bookName: "Pinnacle", price: -107, liquidity: 30 }]);
  assert.deepEqual(check.thin, [{ bookName: "Circa", price: -105, liquidity: 12 }]);
  assert.equal(check.bookCount, 3);
  assert.equal(pricecheck.canSkip(check), true);
  assert.equal(check.wouldSkip, false);
  // A thin book alone beating you is no reason to skip.
  const onlyThin = pricecheck.comparePrice(pricecheck.buildPriceIndex(
    stateWith([[1, 0, 47.5, -110], [2, 0, 47.5, -105, { liquidity: 4 }], [4, 0, 47.5, -115]]), { now: NOW, maxLineAgeMs: null }), overSpec(1, -110));
  assert.equal(pricecheck.canSkip(onlyThin), false);
  assert.deepEqual(pricecheck.priceCheckTag(onlyThin), { kind: "best", label: "best of 2 books" });
});

test("one book under two names is thin only when both copies are", () => {
  const books = {
    1: { id: 1, name: "Novig", isLive: true }, 2: { id: 2, name: "FanDuel", isLive: true },
    107: { id: 107, name: "Polymarket", isLive: true }, 111: { id: 111, name: "Polymarket US", isLive: true },
  };
  const state = stateWith([[1, 0, 47.5, -110], [2, 0, 47.5, -103], [107, 0, 47.5, -105, { liquidity: 4 }], [111, 0, 47.5, -105, { liquidity: 800 }]], books);
  const check = pricecheck.comparePrice(pricecheck.buildPriceIndex(state, { now: NOW, maxLineAgeMs: null }), { ...overSpec(1, -110), bookName: "Novig" });
  assert.deepEqual(check.better, [{ bookName: "FanDuel", price: -103, liquidity: null }, { bookName: "Polymarket / US", price: -105, liquidity: 800 }]);
  assert.deepEqual(check.thin, []);
  assert.equal(check.wouldSkip, true);
});

test("a thin book is not the next best for the outlier check", () => {
  // +300 against a -290 market is an outlier even when a $5 offer sits at +280.
  const state = stateWith([[1, 0, null, 300], [2, 0, null, 280, { liquidity: 5 }], [3, 0, null, -290], [4, 0, null, -300]]);
  const index = pricecheck.buildPriceIndex(state, { now: NOW, maxLineAgeMs: null });
  const check = pricecheck.comparePrice(index, { eventId: 100, periodTypeId: 1, betTypeId: 1, sideIndex: 0, points: null, bookId: 1, price: 300 });
  assert.equal(check.outlier.nextBest.price, -290);
});

test("the Jackson St case: +300 against a -290 market is an outlier, tagged and never a skip", () => {
  const state = stateWith([[1, 0, null, +300], [2, 0, null, -290], [3, 0, null, -300], [4, 0, null, -310]]);
  const index = pricecheck.buildPriceIndex(state, { now: NOW, maxLineAgeMs: null });
  const check = pricecheck.comparePrice(index, { eventId: 100, periodTypeId: 1, betTypeId: 1, sideIndex: 0, points: null, bookId: 1, price: 300 });
  assert.deepEqual(check.outlier, { nextBest: { bookName: "Circa", price: -290 } });
  assert.equal(check.wouldSkip, false);
  assert.deepEqual(pricecheck.priceCheckTag(check), { kind: "outlier", label: "outlier, check · next best -290" });
});

test("a best price a few points clear of the market is not an outlier", () => {
  const state = stateWith([[1, 0, 47.5, +105], [2, 0, 47.5, -110], [3, 0, 47.5, -112]]);
  const index = pricecheck.buildPriceIndex(state, { now: NOW, maxLineAgeMs: null });
  assert.equal(pricecheck.comparePrice(index, overSpec(1, 105)).outlier, null);
});

test("a moneyline matches across books whatever points field they send", () => {
  const state = stateWith([[1, 0, null, +120], [2, 0, null, +130]]);
  state.lines["m1:ms2:si0:null"].points = 0;
  const index = pricecheck.buildPriceIndex(state, { now: NOW, maxLineAgeMs: null });
  const check = pricecheck.comparePrice(index, { eventId: 100, periodTypeId: 1, betTypeId: 1, sideIndex: 0, points: null, bookId: 1, price: 120 });
  assert.equal(check.better.length, 1);
});

test("the header counts outliers; empty when there are none", () => {
  const rows = [
    { priceCheck: { wouldSkip: true, outlier: null } },
    { priceCheck: { wouldSkip: false, outlier: { nextBest: {} } } }, { priceCheck: null },
  ];
  assert.equal(pricecheck.describePriceChecks(rows), "1 outlier to check");
  assert.equal(pricecheck.describePriceChecks([{ priceCheck: { wouldSkip: true, outlier: null } }]), "");
});

test("a skipped row folds only while a book beats it and its game is still to start", () => {
  const beaten = { better: [{ bookName: "Circa", price: -105 }] };
  const row = (eventId, priceCheck, eventStartMs) => ({
    eventId, periodTypeId: 1, betTypeId: 3, sideIndex: 0, points: 47.5, book: { id: 4 }, priceCheck, eventStartMs,
  });
  const skipped = row(1, beaten, NOW + HOUR);
  const nowBest = row(2, { better: [] }, NOW + HOUR);
  const started = row(3, beaten, NOW - HOUR);
  const notSkipped = row(4, beaten, NOW + HOUR);
  const skips = {};
  for (const r of [skipped, nowBest, started]) skips[pricecheck.skipKeyOf(r)] = r.eventStartMs;
  assert.equal(pricecheck.skipKeyOf(skipped), "1:1:3:0|47.5");
  const { listed, skipped: folded } = pricecheck.splitSkipped([skipped, nowBest, started, notSkipped], skips, NOW);
  assert.deepEqual(folded, [skipped]);
  assert.deepEqual(listed, [nowBest, started, notSkipped]);
  assert.equal(pricecheck.canSkip({ better: [] }), false);
  assert.equal(pricecheck.canSkip(null), false);
});

test("a skip covers the line at every book a book beats; the best price there still lists", () => {
  // Cal, 2026-10-16: skipped Pittsburgh -16.5 at FanDuel, then DraftKings'
  // worse -16.5 showed. Now it folds too; the other side and number do not.
  const beaten = { better: [{ bookName: "Circa", price: -105 }] };
  const row = (bookId, priceCheck, overrides = {}) => ({
    eventId: 7, periodTypeId: 1, betTypeId: 2, sideIndex: 0, points: -16.5,
    book: { id: bookId }, priceCheck, eventStartMs: NOW + HOUR, ...overrides,
  });
  const fanduel = row(3, beaten);
  const draftkings = row(5, beaten);
  const bestAtNumber = row(9, { better: [] });
  const otherSide = row(5, beaten, { sideIndex: 1, points: 16.5 });
  const otherNumber = row(5, beaten, { points: -17 });
  const skips = { [pricecheck.skipKeyOf(fanduel)]: fanduel.eventStartMs };
  const { listed, skipped } = pricecheck.splitSkipped([fanduel, draftkings, bestAtNumber, otherSide, otherNumber], skips, NOW);
  assert.deepEqual(skipped, [fanduel, draftkings]);
  assert.deepEqual(listed, [bestAtNumber, otherSide, otherNumber]);
  // A skip saved per book before this change covers the line too.
  const legacy = { "7:1:2:0|-16.5@3": fanduel.eventStartMs };
  assert.deepEqual(pricecheck.splitSkipped([draftkings], legacy, NOW).skipped, [draftkings]);
  assert.equal(pricecheck.skipLineOf("7:1:2:0|-16.5@3"), "7:1:2:0|-16.5");
  assert.equal(pricecheck.skipLineOf("7:1:1:0|ml"), "7:1:1:0|ml");
});

test("every Edges row carries its price check, on the real NFL slice", () => {
  const snapshot = JSON.parse(fs.readFileSync(path.join(__dirname, "fixtures", "v2_slice.json"), "utf8"));
  const state = feed.mergeStates([feed.parseSnapshot(snapshot, { leagueId: 1 })]);
  const now = Date.parse("2026-09-13T16:00:00Z");
  const edgeSettings = { ...edgeRows.DEFAULT_EDGE_SETTINGS, bookIds: null, maxLineAgeHours: 1e6, includeAlts: true, minEdgePct: 1 };
  const rows = edgeRows.sizedEdgeRows(state, {
    edgeSettings, stakeSettings: edgeRows.DEFAULT_STAKE_SETTINGS, records: [], boardLines: edgeRows.boardLines(state),
    ladderReaderOf: edgeRows.createLadderReaders(state), teasers: [], measurement: null,
    effective: edgeRows.effectiveFilter(edgeSettings, state, null, null), now, minEdgePct: 1,
  });
  assert.ok(rows.length > 0);
  for (const row of rows) {
    assert.ok(row.priceCheck, `${row.key} has a price check`);
    assert.ok(row.priceCheck.bookCount >= 1);
    // Your own book is never one of the others.
    assert.ok(!row.priceCheck.better.concat(row.priceCheck.worse).some((entry) => entry.bookName === row.book.name));
  }
});

test("the phone's /edges.json carries each row's check, its tag and the header count", () => {
  const edgesPayload = require("../server/edges_payload.js");
  const tailflex = require("../extension/tailflex.js");
  const snapshot = JSON.parse(fs.readFileSync(path.join(__dirname, "fixtures", "v2_slice.json"), "utf8"));
  const state = feed.mergeStates([feed.parseSnapshot(snapshot, { leagueId: 1 })]);
  const now = Date.parse("2026-09-13T16:00:00Z");
  const edgeSettings = { ...edgeRows.DEFAULT_EDGE_SETTINGS, bookIds: null, maxLineAgeHours: 1e6, includeAlts: true, minEdgePct: 1, groupByMarket: false };
  const payload = edgesPayload.buildEdgesPayload({
    feedState: state, scannerStatus: null, history: {}, betRecords: [], fillFairIndex: new Map(), knownStarts: {},
    stakeSettings: edgeRows.DEFAULT_STAKE_SETTINGS, edgeSettings, settingsStatus: {}, betsStatus: {},
    boardLines: edgeRows.boardLines(state), ladderReaderOf: edgeRows.createLadderReaders(state), teasers: [],
    measurement: tailflex.measureTailFlex(state, { now, maxLineAgeMs: 1e12 }), now,
  });
  // The slice's Bears ML: BetMGM and Sports Interaction -145 lead SouthPoint's -150.
  const southPoint = payload.items.find((row) => row.betTypeId === 1 && row.book.name === "SouthPoint" && row.sideIndex === 0);
  assert.deepEqual(southPoint.priceCheckTag, { kind: "skip", label: "2 books better" });
  assert.deepEqual(southPoint.priceCheck.better.map((entry) => entry.price), [-145, -145]);
  assert.equal(payload.priceCheck, "");
  assert.equal(southPoint.skipKey, pricecheck.skipKeyOf(southPoint, southPoint.book.id));
  assert.deepEqual(payload.skipped, []);

  // Skipped on either page: the row leaves the list for payload.skipped while its game is to start.
  const skipped = edgesPayload.buildEdgesPayload({
    feedState: state, scannerStatus: null, history: {}, betRecords: [], fillFairIndex: new Map(), knownStarts: {},
    stakeSettings: edgeRows.DEFAULT_STAKE_SETTINGS, edgeSettings, settingsStatus: {}, betsStatus: {},
    boardLines: edgeRows.boardLines(state), ladderReaderOf: edgeRows.createLadderReaders(state), teasers: [],
    measurement: tailflex.measureTailFlex(state, { now, maxLineAgeMs: 1e12 }), now,
    edgeSkips: [{ skipKey: southPoint.skipKey, eventStartMs: southPoint.eventStartMs }],
  });
  assert.deepEqual(skipped.skipped.map((row) => row.key), [southPoint.key]);
  assert.ok(!skipped.items.some((row) => row.key === southPoint.key));
  assert.equal(skipped.total, payload.total - 1);
});
