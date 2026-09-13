// Row resolution in page.js, run against the real file: which grid row a
// capture classifies and watches against, and what the watcher reads back
// from it. page.js is a MAIN-world IIFE, so it is loaded into a vm sandbox
// with a fake window/document and exposes its internals only when the
// sandbox sets `__unabatedTicketExposeInternals` (never true on Unabated).
//
// Fixtures replay the two grids that misbehaved live on 2026-09-12: NFL
// CHI@CAR, where the market's rungs are sibling top-level rows sharing one
// grid key, the -2.5 row's entry carrying an alternateLines array with no
// rung in it, next to a main row whose 25-rung ladder holds -2.5; and CFB
// UAPB@ALCN, the same per-rung layout with no ladder anywhere.
const test = require("node:test");
const assert = require("node:assert/strict");
const fs = require("node:fs");
const path = require("node:path");
const vm = require("node:vm");

const PAGE_JS = process.env.PAGE_JS || path.join(__dirname, "..", "extension", "page.js");

function loadPage(pathname = "/nfl/odds") {
  const sandbox = {
    __unabatedTicketExposeInternals: true,
    console: { info() {}, warn() {}, error() {} },
    setInterval: () => 0,
    clearInterval() {},
    Element: class {},
    location: { origin: "https://tools.unabated.com", pathname, href: `https://tools.unabated.com${pathname}` },
    postMessage() {},
    addEventListener() {},
    document: { addEventListener() {}, removeEventListener() {}, querySelectorAll: () => [], querySelector: () => null },
  };
  sandbox.window = sandbox;
  vm.createContext(sandbox);
  vm.runInContext(fs.readFileSync(PAGE_JS, "utf8"), sandbox, { filename: "page.js" });
  assert.ok(sandbox.__unabatedTicketInternals, "page.js exposed nothing to the harness");
  return sandbox.__unabatedTicketInternals;
}

const BOOK = 99;
const BOOK_KEY = `ms${BOOK}`;
const CONTEXT = {
  fullOddsData: { teams: { 1: { name: "Chicago Bears" }, 2: { name: "Carolina Panthers" } } },
  marketSources: [{ id: BOOK, name: "Novig" }],
};

// One grid row of a spread market; `entry` is the book's entry on the HOME side.
function spreadRow(entry, { gridKey = "eid:125807:pid:tid5:eid:tid5:bt2:pt1:bst:livefalse:proj0", bestLines = { a: 1, b: 1 } } = {}) {
  return {
    gridKey, eventId: 125807, betTypeId: 2, periodTypeId: 1, bestLines,
    eventName: "Chicago Bears @ Carolina Panthers", eventStart: "2026-09-13T17:00:00Z",
    eventTeams: [{ id: 1, rotationNumber: 465 }, { id: 2, rotationNumber: 466 }],
    sides: { "si0:tid1": {}, "si1:tid2": { [BOOK_KEY]: entry } },
  };
}

function topNode(data, id) {
  return { id, data, level: 0, detail: false, parent: { data: null } };
}

function childNode(data, id, parent) {
  return { id, data, level: 1, detail: true, parent };
}

// A fake AG Grid API. `keyed` decides what getRowNode answers when several
// rows share a key (AG Grid's answer for duplicate ids is not specified).
function gridApi(nodes, keyed = "first") {
  return {
    nodes,
    forEachNode(fn) { nodes.forEach(fn); },
    getRowNode(key) {
      const hits = nodes.filter((node) => node.data && node.data.gridKey === key);
      return keyed === "last" ? hits[hits.length - 1] || null : hits[0] || null;
    },
    getDisplayedRowCount() { return nodes.length; },
    isDestroyed() { return false; },
  };
}

function capture(page, api, rowData, marketLine, sideIndex = 1) {
  const { gridApi: pickedApi, ...ticket } = page.buildTicket({ marketLine, sideIndex, rowData, context: CONTEXT, cellProps: null, gridApi: api });
  page.startWatching(ticket, pickedApi);
  return ticket;
}

// ---- NFL CHI@CAR: per-rung rows next to a main row carrying the ladder ----

function nflGrid() {
  const alt25 = { points: -2.5, price: 182, ge: 0.0483, bacr: 169 };
  const ladder = [-14.5, -13.5, -10.5, -9.5, -7.5, -6.5, -2.5, 3.5, 6.5].map((points) =>
    (points === -2.5 ? alt25 : { points, price: points < 0 ? 400 : -127 }));
  const main = spreadRow({ points: 3.5, price: -127, alternateLines: ladder });
  const clicked = { points: -2.5, price: 182, ge: 0.0483, bacr: 169, alternateLines: [null] };
  const rungs = [
    spreadRow(clicked),
    spreadRow({ points: -14.5, price: 1076, alternateLines: [] }),
    spreadRow({ points: -13.5, price: 809, alternateLines: [] }),
    spreadRow({ points: -6.5, price: 388, alternateLines: [] }),
  ];
  const nodes = [topNode(rungs[0], "r0"), topNode(main, "main"), ...rungs.slice(1).map((row, i) => topNode(row, `r${i + 1}`))];
  return { nodes, main, clicked, alt25, rungRow: rungs[0] };
}

test("NFL per-rung grid: the click on -2.5 resolves without a tie and prices the clicked rung", () => {
  const page = loadPage();
  const { nodes, clicked, rungRow } = nflGrid();
  const ticket = capture(page, gridApi(nodes), rungRow, clicked);
  assert.equal(ticket.rowResolution.ambiguous, false, ticket.rowResolution.trace);
  assert.equal(ticket.points, -2.5);
  assert.equal(ticket.price, 182);
  assert.equal(ticket.sideLabel, "Carolina Panthers -2.5");
  assert.equal(ticket.book.name, "Novig");
  assert.equal(ticket.isAlt, true, "a rung next to a main row carrying the ladder is an alt line");
});

for (const keyed of ["first", "last"]) {
  test(`NFL per-rung grid: the watcher reads the -2.5 rung whichever sibling the shared key answers with (${keyed})`, () => {
    const page = loadPage();
    const { nodes, clicked, rungRow } = nflGrid();
    capture(page, gridApi(nodes, keyed), rungRow, clicked);
    const read = page.readWatchedLine();
    assert.equal(read.offBoard, false, read.row);
    assert.equal(read.points, -2.5);
    assert.equal(read.price, 182);
  });
}

test("NFL per-rung grid: the rung gone from both the ladder and the rows reads off the board at the captured price", () => {
  const page = loadPage();
  const { nodes, clicked, rungRow, main, alt25 } = nflGrid();
  capture(page, gridApi(nodes), rungRow, clicked);
  main.sides["si1:tid2"][BOOK_KEY].alternateLines = main.sides["si1:tid2"][BOOK_KEY].alternateLines.filter((alt) => alt !== alt25);
  nodes.splice(0, 1);
  const read = page.readWatchedLine();
  assert.equal(read.offBoard, true);
  assert.equal(read.points, -2.5);
  assert.equal(read.price, 182);
});

// ---- CFB UAPB@ALCN: per-rung rows, no ladder anywhere ----

function cfbGrid() {
  const rows = [33.5, 37.5, 41.5, 42.5, 43.5, 44.5, 45.5, 46.5, 56.5].map((points) => {
    const price = points === 56.5 ? 199 : -1718;
    const row = spreadRow({ points, price, alternateLines: [] });
    row.betTypeId = 3;
    row.gridKey = "eid:123733:pid:tid:eid:tid:bt3:pt1:bst:livefalse:proj0";
    row.sides = { "si0:tid1": { [BOOK_KEY]: row.sides["si1:tid2"][BOOK_KEY] }, "si1:tid2": {} };
    return row;
  });
  return { nodes: rows.map((row, i) => topNode(row, `t${i}`)), rows };
}

test("CFB per-rung grid: the click on Over 56.5 picks its own row, no tie, not an alt", () => {
  const page = loadPage("/cfb/odds");
  const { nodes, rows } = cfbGrid();
  const row = rows[rows.length - 1];
  const ticket = capture(page, gridApi(nodes), row, row.sides["si0:tid1"][BOOK_KEY], 0);
  assert.equal(ticket.rowResolution.ambiguous, false, ticket.rowResolution.trace);
  assert.equal(ticket.points, 56.5);
  assert.equal(ticket.price, 199);
  assert.equal(ticket.sideLabel, "Over 56.5");
  assert.equal(ticket.isAlt, false);
  const read = page.readWatchedLine();
  assert.equal(read.offBoard, false, read.row);
  assert.equal(read.points, 56.5);
  assert.equal(read.price, 199);
});

// ---- an ordinary grid with an expanded Alts section ----

function altsGrid() {
  const rung = { points: -3.5, price: 150 };
  const parentData = spreadRow({ points: -5.5, price: -110, alternateLines: [{ points: -7.5, price: 210 }, rung, { points: -1.5, price: 110 }] });
  const parent = topNode(parentData, "p");
  // The child's entry is the rung itself; whether it is the SAME object as
  // the parent ladder's element is unknown, so the fixture uses a copy.
  const childRung = { points: -3.5, price: 150 };
  const childData = spreadRow(childRung, { gridKey: "p-alt-1", bestLines: null });
  const child = childNode(childData, "c1", parent);
  return { nodes: [parent, child], parentData, childData, childRung, rung };
}

test("Alts child click: the top-level parent carrying the ladder is picked over the child, as an alt line", () => {
  const page = loadPage();
  const { nodes, childData, childRung } = altsGrid();
  const ticket = capture(page, gridApi(nodes), childData, childRung);
  assert.equal(ticket.rowResolution.ambiguous, false, ticket.rowResolution.trace);
  assert.match(ticket.rowResolution.trace, /^picked \[top L0/);
  assert.equal(ticket.isAlt, true);
  assert.equal(ticket.watch.rowTop, true);
  assert.equal(ticket.points, -3.5);
});

test("Alts child click: after the section closes the watcher still reads the rung by number, never the parent's main line", () => {
  const page = loadPage();
  const { nodes, childData, childRung, rung, parentData } = altsGrid();
  capture(page, gridApi(nodes), childData, childRung);
  nodes.splice(1, 1);
  rung.price = 160;
  let read = page.readWatchedLine();
  assert.equal(read.points, -3.5);
  assert.equal(read.price, 160);
  assert.equal(read.rowTop, true);
  parentData.sides["si1:tid2"][BOOK_KEY].alternateLines = [];
  read = page.readWatchedLine();
  assert.equal(read.offBoard, true);
});

// ---- a plain main-line capture ----

test("main-line capture: a price move at the number is read; the number moving off the row reads off the board", () => {
  const page = loadPage();
  const entry = { points: -5.5, price: -110, alternateLines: [] };
  const row = spreadRow(entry);
  const ticket = capture(page, gridApi([topNode(row, "p")]), row, entry);
  assert.equal(ticket.isAlt, false);
  assert.equal(ticket.rowResolution.ambiguous, false);
  entry.price = -118;
  let read = page.readWatchedLine();
  assert.equal(read.price, -118);
  assert.equal(read.offBoard, false);
  entry.points = -6;
  entry.alternateLines = [{ points: -5.5, price: -100 }];
  read = page.readWatchedLine();
  assert.equal(read.points, -5.5, "the captured number is followed into the ladder");
  assert.equal(read.price, -100);
  entry.alternateLines = [];
  read = page.readWatchedLine();
  assert.equal(read.offBoard, true);
});
