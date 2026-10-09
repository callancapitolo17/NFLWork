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

test("counts books that beat the price at the same number, best first; 2+ would skip", () => {
  const state = stateWith([
    [1, 0, 47.5, -112], [2, 0, 47.5, -105], [3, 0, 47.5, -107], [4, 0, 47.5, -115], [5, 0, 47.5, -112],
  ]);
  const index = pricecheck.buildPriceIndex(state, { now: NOW, maxLineAgeMs: 168 * HOUR });
  const check = pricecheck.comparePrice(index, overSpec(1, -112));
  assert.deepEqual(check.better, [{ bookName: "Circa", price: -105 }, { bookName: "Pinnacle", price: -107 }]);
  assert.deepEqual(check.worse, [{ bookName: "DraftKings", price: -115 }]);
  assert.equal(check.same, 1);
  assert.equal(check.bookCount, 5);
  assert.equal(check.wouldSkip, true);
  assert.equal(check.outlier, null);
  assert.deepEqual(pricecheck.priceCheckTag(check), { kind: "skip", label: "2 books better · would skip" });
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
  assert.deepEqual(check.better, [{ bookName: "Circa", price: -105 }]);
  assert.equal(check.bookCount, 3);
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

test("the header counts would-skip rows and outliers; empty when there are none", () => {
  const rows = [
    { priceCheck: { wouldSkip: true, outlier: null } }, { priceCheck: { wouldSkip: true, outlier: null } },
    { priceCheck: { wouldSkip: false, outlier: { nextBest: {} } } }, { priceCheck: null },
  ];
  assert.equal(pricecheck.describePriceChecks(rows), "2 would skip: 2+ books better · 1 outlier (test mode, nothing hidden)");
  assert.equal(pricecheck.describePriceChecks([{ priceCheck: { wouldSkip: false, outlier: null } }]), "");
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
  assert.deepEqual(southPoint.priceCheckTag, { kind: "skip", label: "2 books better · would skip" });
  assert.deepEqual(southPoint.priceCheck.better.map((entry) => entry.price), [-145, -145]);
  assert.equal(payload.priceCheck, "1 would skip: 2+ books better (test mode, nothing hidden)");
});
