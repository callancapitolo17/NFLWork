// Run: node --test unabated_ticket/tests
// server/edgelog.js: which board lines are logged as edges, and how the
// segments open, stay, close and are re-sent. The POST itself runs end to
// end in runner.test.js.
const test = require("node:test");
const assert = require("node:assert/strict");
const fs = require("node:fs");
const path = require("node:path");
const feed = require("../extension/feed.js");
const edgelog = require("../server/edgelog.js");

const NFL = 1;
const KICKOFF_MS = Date.parse("2026-09-13T17:00:00Z");
const OBSERVED_MS = KICKOFF_MS - 3600 * 1000;
const MINUTE_MS = 60 * 1000;

function sliceState() {
  const snapshot = JSON.parse(fs.readFileSync(path.join(__dirname, "fixtures", "v2_slice.json"), "utf8"));
  return feed.mergeStates([feed.parseSnapshot(snapshot, { leagueId: NFL })]);
}

function sliceEdgeLines(state, observedMs = OBSERVED_MS) {
  return edgelog.leagueEdgeLines({ lines: Object.values(state.lines), state, observedMs });
}

test("leagueEdgeLines: the Edges tab's own lines at a 1% floor, every live book, alts included, full game", () => {
  const state = sliceState();
  const logged = sliceEdgeLines(state);
  const edgesTab = feed.selectEdges(state, {
    minEdge: edgelog.MIN_LOGGED_EDGE, includeAlts: true, periods: new Set([1]), betTypes: new Set([1, 2, 3]), now: OBSERVED_MS,
  }).filter((row) => row.isBlurred !== true);
  assert.ok(logged.length > 5, "the slice has edges to log");
  assert.deepEqual(logged.map((row) => row.lineKey).sort(), edgesTab.map((row) => row.key).sort());
  assert.ok(logged.some((row) => row.isAlt), "alt rungs are logged");
  assert.ok(logged.every((row) => row.edge >= edgelog.MIN_LOGGED_EDGE && row.observedMs === OBSERVED_MS));
  const moneyline = logged.find((row) => row.betTypeId === 1);
  assert.equal(moneyline.points, null);
});

test("leagueEdgeLines: nothing once the snapshot landed at or after the start, or with no snapshot time", () => {
  const state = sliceState();
  assert.deepEqual(sliceEdgeLines(state, KICKOFF_MS), []);
  assert.deepEqual(sliceEdgeLines(state, Number.NaN), []);
});

test("leagueEdgeLines: blurred, clamped, Unabated's own, off-board and dead-book lines are not edges", () => {
  const state = sliceState();
  const [first] = sliceEdgeLines(state);
  const line = state.lines[first.lineKey];
  const loggedWith = (fields, bookFields = {}) => {
    const changed = { ...state, lines: { ...state.lines, [line.key]: { ...line, ...fields } },
      books: { ...state.books, [line.bookId]: { ...state.books[line.bookId], ...bookFields } } };
    return sliceEdgeLines(changed).some((row) => row.lineKey === line.key);
  };
  assert.equal(loggedWith({}), true);
  assert.equal(loggedWith({ isBlurred: true }), false);
  assert.equal(loggedWith({ bacr: 999900 }), false);
  assert.equal(loggedWith({ statusId: 2 }), false);
  assert.equal(loggedWith({ ge: 0.0099 }), false);
  assert.equal(loggedWith({ periodTypeId: 2 }), false);
  assert.equal(loggedWith({}, { isLive: false }), false);
});

function edgeLine(fields = {}) {
  return {
    lineKey: "m1:ms4:si0", eventId: 7, leagueId: NFL, bookId: 4, betTypeId: 2, sideIndex: 0, points: -3, isAlt: false,
    eventStartMs: KICKOFF_MS, price: 110, sourcePrice: 110, sourceFormat: 1, fairAmerican: -105, edge: 0.07, liquidity: null,
    observedMs: OBSERVED_MS, signature: "-3|110|110|-105", ...fields,
  };
}

function advance(segments, edgeLines, now = OBSERVED_MS, readLeagueIds = [NFL]) {
  return edgelog.advanceSegments({ edgeLines, readLeagueIds, segments, now });
}

test("advanceSegments: a new edge opens a segment, sent once, then again only every 5 minutes while it lasts", () => {
  const segments = new Map();
  const [opened] = advance(segments, [edgeLine()]);
  assert.deepEqual(opened, {
    lineKey: "m1:ms4:si0", segmentStart: new Date(OBSERVED_MS).toISOString(), lastSeen: new Date(OBSERVED_MS).toISOString(),
    eventId: 7, leagueId: NFL, bookId: 4, betTypeId: 2, sideIndex: 0, points: -3, isAlt: false,
    eventStart: new Date(KICKOFF_MS).toISOString(), price: 110, sourcePrice: 110, sourceFormat: 1, fairAmerican: -105, edge: 0.07, liquidity: null,
  });
  assert.deepEqual(advance(segments, [edgeLine()]), [], "the same snapshot again sends nothing");
  assert.deepEqual(advance(segments, [edgeLine({ observedMs: OBSERVED_MS + 4 * MINUTE_MS })]), []);
  const [confirmed] = advance(segments, [edgeLine({ observedMs: OBSERVED_MS + 5 * MINUTE_MS })]);
  assert.deepEqual([confirmed.segmentStart, confirmed.lastSeen],
    [new Date(OBSERVED_MS).toISOString(), new Date(OBSERVED_MS + 5 * MINUTE_MS).toISOString()]);
});

test("advanceSegments: a new number, price or fair closes the segment where it was last seen and opens the next", () => {
  const segments = new Map();
  advance(segments, [edgeLine()]);
  advance(segments, [edgeLine({ observedMs: OBSERVED_MS + MINUTE_MS })]);
  const rows = advance(segments, [edgeLine({ observedMs: OBSERVED_MS + 2 * MINUTE_MS, price: 105, signature: "-3|105|105|-105" })]);
  assert.deepEqual(rows.map((row) => [row.segmentStart, row.lastSeen, row.price]), [
    [new Date(OBSERVED_MS).toISOString(), new Date(OBSERVED_MS + MINUTE_MS).toISOString(), 110],
    [new Date(OBSERVED_MS + 2 * MINUTE_MS).toISOString(), new Date(OBSERVED_MS + 2 * MINUTE_MS).toISOString(), 105],
  ]);
  assert.equal(segments.size, 1);
});

test("advanceSegments: an edge gone from a league read this pass ends; one of a league not read stays open", () => {
  const segments = new Map();
  advance(segments, [edgeLine()]);
  assert.deepEqual(advance(segments, [], OBSERVED_MS + MINUTE_MS, []), [], "its league was not read: left open");
  assert.equal(segments.size, 1);
  advance(segments, [edgeLine({ observedMs: OBSERVED_MS + MINUTE_MS })], OBSERVED_MS + MINUTE_MS);
  const [ended] = advance(segments, [], OBSERVED_MS + 2 * MINUTE_MS);
  assert.equal(ended.lastSeen, new Date(OBSERVED_MS + MINUTE_MS).toISOString(), "ends at its last sighting, unsent until now");
  assert.equal(segments.size, 0);
  advance(segments, [edgeLine()]);
  assert.deepEqual(advance(segments, [], OBSERVED_MS + MINUTE_MS), [], "an end already sent is not sent again");
});

test("advanceSegments: kickoff ends every open segment, and a cached line of a started game opens none", () => {
  const segments = new Map();
  advance(segments, [edgeLine()]);
  advance(segments, [edgeLine({ observedMs: OBSERVED_MS + MINUTE_MS })]);
  const rows = advance(segments, [edgeLine({ observedMs: OBSERVED_MS + MINUTE_MS })], KICKOFF_MS, []);
  assert.deepEqual(rows.map((row) => row.lastSeen), [new Date(OBSERVED_MS + MINUTE_MS).toISOString()]);
  assert.equal(segments.size, 0);
  assert.deepEqual(advance(segments, [edgeLine()], KICKOFF_MS), []);
  assert.equal(segments.size, 0);
});

test("segmentKeyOf: one key per line and segment start", () => {
  const segments = new Map();
  const [row] = advance(segments, [edgeLine()]);
  assert.equal(edgelog.segmentKeyOf(row), `m1:ms4:si0|${new Date(OBSERVED_MS).toISOString()}`);
});
