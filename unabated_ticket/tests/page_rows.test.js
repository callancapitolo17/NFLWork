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

// A zone west of UTC, so a naive UTC timestamp read as local time lands
// hours late and the started-game gate below fails loudly (it passed in UTC).
process.env.TZ = "America/Los_Angeles";

const PAGE_JS = process.env.PAGE_JS || path.join(__dirname, "..", "extension", "page.js");
// selfcheck.js is injected immediately before page.js (manifest, MAIN world)
// and page.js takes its reference off the page global. The sandbox loads it
// the same way; `selfcheck` here is the same module under require().
const SELF_CHECK_JS = path.join(__dirname, "..", "extension", "selfcheck.js");
const selfcheck = require("../extension/selfcheck.js");

// `page.window` is the sandbox: `posted` collects every window.postMessage
// (what content.js would receive), `deliver(type, payload)` plays a message
// from content.js into page.js's listener, `watchTimers()` counts watch intervals started.
function loadPage(pathname = "/nfl/odds", { gridRoots = [], shells = [], rowElements = [], rendered = (selector) => selector.startsWith(".ag-row") } = {}) {
  const listeners = [];
  const posted = [];
  let timers = 0;
  const sandbox = {
    __unabatedTicketExposeInternals: true,
    console: { info() {}, warn() {}, error() {} },
    setInterval: () => { timers += 1; return timers; },
    clearInterval() {},
    Element: class {},
    location: { origin: "https://tools.unabated.com", pathname, href: `https://tools.unabated.com${pathname}` },
    postMessage(data) { posted.push(data); },
    addEventListener(type, fn) { if (type === "message") listeners.push(fn); },
    document: {
      addEventListener() {}, removeEventListener() {},
      // Only the locate flow and the self-check's locate probe use
      // querySelector; both ask for `.ag-row[row-id=...]` and a cell inside it.
      querySelector: (selector) => (rendered(selector) ? { selector } : null),
      querySelectorAll: (selector) => {
        if (selector === ".ag-root") return gridRoots;
        if (selector === ".ag-row") return rowElements;
        if (selector === ".odds-cell-action-shell") return shells;
        return [];
      },
    },
  };
  sandbox.window = sandbox;
  vm.createContext(sandbox);
  vm.runInContext(fs.readFileSync(SELF_CHECK_JS, "utf8"), sandbox, { filename: "selfcheck.js" });
  vm.runInContext(fs.readFileSync(PAGE_JS, "utf8"), sandbox, { filename: "page.js" });
  assert.ok(sandbox.__unabatedTicketInternals, "page.js exposed nothing to the harness");
  const timersAtLoad = timers; // the heartbeat
  // page.js checks event.source against its own `window`, which inside the
  // context is the contextified global, not the outer sandbox object.
  const innerWindow = vm.runInContext("window", sandbox);
  const deliver = (type, payload) => listeners.forEach((fn) => fn({ source: innerWindow, data: { source: "unabated-ticket", type, payload } }));
  return { ...sandbox.__unabatedTicketInternals, posted, deliver, watchTimers: () => timers - timersAtLoad };
}

// A rendered grid root whose one cell reaches `api` through fake React fiber
// props, the way allGridApis re-acquires a grid a resumed watcher never held.
function gridRoot(api) {
  const cell = { "__reactFiber$test": { memoizedProps: { api, node: { data: {} } }, return: null } };
  return { querySelectorAll: (selector) => (selector === ".ag-cell" ? [cell] : []) };
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
  // AG Grid row nodes carry setExpanded; the Alts probe and locate both need it.
  return { id, data, level: 0, detail: false, parent: { data: null }, setExpanded() {} };
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

// ---- resume after a navigation: the stored ticket rebuilds the watcher ----
//
// page.js dies with every navigation (one-click betting leaving the tab,
// Back, a reload or discard, an extension reload's takeover) while the ticket
// stays in storage; before 0.6.6 nothing re-attached, and the panel read
// "Not watching the line" with a live heartbeat for every such ticket.

function storedTicketAfterCapture() {
  const page = loadPage();
  const { nodes, clicked, rungRow } = nflGrid();
  const ticket = capture(page, gridApi(nodes), rungRow, clicked);
  const stored = JSON.parse(JSON.stringify(ticket)); // structured clone, as storage hands it back
  // Kickoff a day out in the grid's naive-UTC format, so the started-game
  // gate never depends on when the suite runs (the fixture's own date passed).
  stored.eventStart = new Date(Date.now() + 86400000).toISOString().slice(0, 19);
  return stored;
}

test("a freshly loaded page.js asks content.js for the stored ticket", () => {
  const page = loadPage();
  assert.ok(page.posted.some((message) => message.type === "resume_request"), "no resume_request posted on load");
});

test("resume: the watcher is rebuilt from the stored ticket with no grid API and reads the rung off the DOM's grid", () => {
  const ticket = storedTicketAfterCapture();
  const { nodes } = nflGrid();
  const page = loadPage("/nfl/odds", { gridRoots: [gridRoot(gridApi(nodes))] });
  page.deliver("resume", ticket);
  assert.equal(page.watchTimers(), 1, "one watch interval started");
  const read = page.readWatchedLine();
  assert.equal(read.offBoard, false, read.row);
  assert.equal(read.points, -2.5);
  assert.equal(read.price, 182);
});

test("resume: a second offer of the same ticket does not restart the watcher; a newer ticket replaces it, an older one does not", () => {
  const ticket = storedTicketAfterCapture();
  const { nodes } = nflGrid();
  const page = loadPage("/nfl/odds", { gridRoots: [gridRoot(gridApi(nodes))] });
  page.deliver("resume", ticket);
  page.deliver("resume", ticket);
  assert.equal(page.watchTimers(), 1);
  page.deliver("resume", { ...ticket, capturedAt: ticket.capturedAt + 1 });
  assert.equal(page.watchTimers(), 2);
  page.deliver("resume", ticket);
  assert.equal(page.watchTimers(), 2, "an older ticket read from storage must not replace a newer capture");
});

test("resume: a tab whose grid has not mounted yet says so, then reads once it has", () => {
  const ticket = storedTicketAfterCapture();
  const roots = [];
  const page = loadPage("/nfl/odds", { gridRoots: roots });
  page.deliver("resume", ticket);
  assert.throws(() => page.readWatchedLine(), /no odds grid reachable on the page/);
  roots.push(gridRoot(gridApi(nflGrid().nodes)));
  assert.equal(page.readWatchedLine().price, 182);
});

test("resume: a malformed or watch-less ticket, or one arriving before any capture, is ignored", () => {
  const page = loadPage();
  page.deliver("resume", null);
  page.deliver("resume", { capturedAt: "x", watch: {} });
  page.deliver("resume", { capturedAt: Date.now() });
  assert.equal(page.watchTimers(), 0);
});

test("resume: a tab on another league, or a ticket whose game has started, is left alone", () => {
  const ticket = storedTicketAfterCapture();
  const otherLeague = loadPage("/cfb/odds", { gridRoots: [gridRoot(gridApi(nflGrid().nodes))] });
  otherLeague.deliver("resume", ticket);
  assert.equal(otherLeague.watchTimers(), 0, "a CFB tab must not watch an NFL ticket");
  const started = loadPage("/nfl/odds", { gridRoots: [gridRoot(gridApi(nflGrid().nodes))] });
  started.deliver("resume", { ...ticket, eventStart: new Date(Date.now() - 60000).toISOString() });
  assert.equal(started.watchTimers(), 0, "a started game is not resumed");
  // The grid's own format: naive UTC, one hour ago.
  const naiveHourAgo = new Date(Date.now() - 3600000).toISOString().slice(0, 19);
  started.deliver("resume", { ...ticket, eventStart: naiveHourAgo });
  assert.equal(started.watchTimers(), 0, `naive UTC ${naiveHourAgo} read as local time`);
  started.deliver("resume", { ...ticket, eventStart: null });
  assert.equal(started.watchTimers(), 1, "no start time known: resumed");
});

// ---- the bundle self-check (#123) ----------------------------------------
//
// page.js dies with every navigation and Unabated deploys mid-session, so the
// shapes it reads are re-probed on each load and on a timer. These run the
// real page.js against a fiber-shaped board and assert what the panel would
// be told: every probe passing on today's shape, and the probe NAMED when one
// field is renamed out from under it.

const USER_SETTINGS = {
  gameOdds: { books: [{ marketSourceId: BOOK, isUnavailable: false }, { marketSourceId: 4, isUnavailable: true }] },
};

// A book entry carrying every field page.js reads off a line.
function probeLine(overrides = {}) {
  return {
    points: -2.5, price: 182, americanPrice: 182, ge: 0.0483, edge: { edge: 4.83 },
    bacr: 169, sourcePrice: 0.352, sourceFormat: 4, statusId: 1, alternateLines: [],
    ...overrides,
  };
}

// A rendered odds cell: the price shell's fiber walks up to the line
// renderer's props and then to the AG Grid cell renderer's (api, node,
// context) — the same walk a click makes.
function oddsCell({ api, node, marketLine, sideIndex = 1 }) {
  const context = { ...CONTEXT, userSettings: USER_SETTINGS };
  const gridFiber = { memoizedProps: { api, node, context }, return: null };
  const lineFiber = { memoizedProps: { marketLine, sideIndex, context }, return: gridFiber };
  return {
    "__reactFiber$test": lineFiber,
    getAttribute: (name) => (name === "data-side-index" ? String(sideIndex) : null),
  };
}

// One market's rows plus the cell rendered over the clicked one.
function probeBoard(line = probeLine()) {
  const row = spreadRow(line);
  const node = topNode(row, "r0");
  const api = gridApi([node]);
  return {
    line, row, node, api,
    shells: [oddsCell({ api, node, marketLine: line })],
    rowElements: [{ id: "r0" }],
    gridRoots: [gridRoot(api)],
  };
}

function probePage(board = probeBoard(), options = {}) {
  return loadPage("/nfl/odds", { gridRoots: board.gridRoots, shells: board.shells, rowElements: board.rowElements, ...options });
}

function lastReport(page) {
  const reports = page.posted.filter((message) => message.type === "selfcheck");
  assert.ok(reports.length, "page.js published no self-check");
  return reports[reports.length - 1].payload;
}

function stateOf(report, id) {
  const probe = report.probes.find((entry) => entry.id === id);
  assert.ok(probe, `no probe ${id} in the report`);
  return probe.state;
}

test("self-check: today's board passes every probe and says nothing", () => {
  const report = lastReport(probePage());
  const notPassing = report.probes.filter((probe) => probe.state !== "pass").map((probe) => `${probe.id}: ${probe.detail}`);
  assert.equal(notPassing.join(" | "), "", "probes did not pass on a board shaped like the live one");
  assert.equal(report.status, "ok");
  assert.equal(selfcheck.summaryOf(report), null);
  assert.equal(selfcheck.noteFor(report, "ticket"), null);
  assert.equal(selfcheck.noteFor(report, "edges"), null);
});

// The issue's acceptance test: rename one probed prop, and the panel must name
// that probe rather than showing an empty pane.
for (const [renamed, probeId, label] of [
  ["bacr", "fair", "Unabated's fair price (bacr)"],
  ["sourcePrice", "sourcePrice", "the exchange source price"],
]) {
  test(`self-check: renaming ${renamed} on the board's lines names that probe`, () => {
    const line = probeLine();
    line[`${renamed}Renamed`] = line[renamed];
    delete line[renamed];
    const report = lastReport(probePage(probeBoard(line)));
    assert.equal(report.status, "changed");
    assert.equal(report.failed.join(","), probeId);
    assert.equal(selfcheck.summaryOf(report), `Unabated changed: ${label}`);
    assert.match(selfcheck.noteFor(report, "ticket"), new RegExp(label.replace(/[()]/g, "\\$&")));
    // A line-shape failure is the ticket's, not the Edges list's.
    assert.equal(selfcheck.noteFor(report, "edges"), null);
  });
}

test("self-check: the edge field moving names the edge probe, whichever of edge/ge it was", () => {
  const line = probeLine();
  delete line.edge;
  delete line.ge;
  const report = lastReport(probePage(probeBoard(line)));
  assert.equal(report.failed.join(","), "edge");
  assert.equal(stateOf(report, "fair"), "pass");
});

test("self-check: a renamed cell class is a change, not a slow load", () => {
  const board = probeBoard();
  // Rows rendered; nothing in them matched .odds-cell-action-shell.
  const report = lastReport(probePage({ ...board, shells: [] }));
  assert.equal(report.status, "changed");
  assert.equal(report.failed.join(","), "cells");
  assert.match(selfcheck.summaryOf(report), /odds cells/);
  // Nothing else can be judged without a cell, and nothing else is blamed.
  assert.equal(stateOf(report, "fiber"), "unknown");
  assert.equal(stateOf(report, "books"), "unknown");
});

test("self-check: a tab with no grid yet is waiting, and the panel stays quiet", () => {
  const report = lastReport(loadPage("/nfl/odds"));
  assert.equal(report.status, "waiting");
  assert.equal(report.failed.join(","), "");
  assert.equal(selfcheck.summaryOf(report), null);
  assert.equal(report.probes.every((probe) => probe.state === "unknown"), true);
});

test("self-check: a grid up with no rows is an empty board, not a change", () => {
  const board = probeBoard();
  // The user filtered every game away, or the slate is over: the grid is
  // mounted and empty. Shouting "Unabated changed" here would train the
  // reader to ignore the banner.
  const report = lastReport(probePage({ ...board, shells: [], rowElements: [] }));
  assert.equal(report.status, "waiting");
  assert.equal(selfcheck.summaryOf(report), null);
});

test("self-check: every book unticked is the user's selection, not a moved field", () => {
  const board = probeBoard();
  const settings = { gameOdds: { books: [{ marketSourceId: BOOK, isUnavailable: true }] } };
  for (const fiber of [board.shells[0].__reactFiber$test, board.shells[0].__reactFiber$test.return]) {
    fiber.memoizedProps.context = { ...CONTEXT, userSettings: settings };
  }
  const report = lastReport(probePage(board));
  assert.equal(report.status, "ok");
  assert.equal(stateOf(report, "books"), "pass");
});

test("self-check: the book selection failing is the Edges tab's problem, not the ticket's", () => {
  const board = probeBoard();
  const cell = board.shells[0];
  // userSettings kept, gameOdds gone: the shape page.js reads the Unabated
  // book selection from (live 2026-09-10 it sat one level down under it).
  for (const fiber of [cell.__reactFiber$test, cell.__reactFiber$test.return]) {
    fiber.memoizedProps.context = { ...CONTEXT, userSettings: { theme: "dark" } };
  }
  const report = lastReport(probePage(board));
  assert.equal(report.failed.join(","), "books");
  assert.match(selfcheck.noteFor(report, "edges"), /book selection/);
  assert.equal(selfcheck.noteFor(report, "ticket"), null);
});

test("self-check: a row node without setExpanded fails the Alts probe", () => {
  const board = probeBoard();
  delete board.node.setExpanded;
  const report = lastReport(probePage(board));
  assert.equal(report.failed.join(","), "alts");
  assert.match(selfcheck.summaryOf(report), /Alts expander/);
  // Alts is the locate path too, so both views are told.
  assert.ok(selfcheck.noteFor(report, "ticket"));
  assert.ok(selfcheck.noteFor(report, "edges"));
});

test("self-check: the row element the locate path keys on going missing names the lookup", () => {
  const report = lastReport(probePage(probeBoard(), { rendered: () => false }));
  assert.equal(report.failed.join(","), "cellLookup");
  assert.match(report.probes.find((probe) => probe.id === "cellLookup").detail, /\.ag-row/);
});

test("self-check: the published report carries verdicts only, never the sampled lines", () => {
  const report = lastReport(probePage());
  const json = JSON.stringify(report);
  assert.equal(json.includes("alternateLines"), false, "a raw board line reached storage");
  assert.ok(json.length < 4096, `report is ${json.length} bytes; it is written to storage on every probe`);
});

test("self-check: it is re-run on a timer, so a mid-session deploy surfaces without a reload", () => {
  const board = probeBoard();
  const page = probePage(board);
  assert.equal(lastReport(page).status, "ok");
  // page.js registers the self-check interval at load; the callback is what a
  // later tick runs, so re-running it against a moved board is the same path.
  delete board.line.bacr;
  page.runSelfCheck();
  assert.equal(lastReport(page).failed.join(","), "fair");
});

test("self-check: a sample that throws is reported as a change, not swallowed", () => {
  const report = selfcheck.errorReport({ at: Date.now(), url: "https://tools.unabated.com/nfl/odds", build: "x", message: "boom" });
  assert.equal(report.status, "changed");
  assert.match(selfcheck.summaryOf(report), /reading the page/);
  assert.match(selfcheck.detailOf(report), /boom/);
});

test("self-check: an odds screen that has rendered nothing a minute after load is a change; a fresh one is waiting", () => {
  const bare = { path: "/nfl/odds", shells: 0, rowElements: 0, gridRoots: 0, lines: [] };
  assert.equal(selfcheck.run({ ...bare, sinceLoadMs: 5000 }).status, "waiting");
  const late = selfcheck.run({ ...bare, sinceLoadMs: 61000 });
  assert.equal(late.status, "changed");
  assert.equal(late.failed.join(","), "cells");
  // Not the odds screen: nothing is expected to render, however long ago.
  assert.equal(selfcheck.run({ ...bare, path: "/account", sinceLoadMs: 600000 }).status, "waiting");
});

test("self-check: a malformed sample never throws out of the module", () => {
  for (const sample of [null, undefined, 42, {}, { shells: 3, lines: "no", cell: "no" }]) {
    const report = selfcheck.run(sample);
    assert.ok(["ok", "changed", "waiting"].includes(report.status), `status ${report.status} for ${JSON.stringify(sample)}`);
  }
});

test("self-check: a best-line cell first in DOM order does not get the bundle blamed for its missing book", () => {
  const board = probeBoard();
  // A cell whose marketLine is a copy the row does not hold and which sits in
  // a column with no book id: bookIdOf falls all the way through on it.
  const orphan = oddsCell({ api: board.api, node: board.node, marketLine: { ...board.line } });
  const report = lastReport(probePage({ ...board, shells: [orphan, ...board.shells] }));
  assert.equal(report.status, "ok", report.probes.filter((probe) => probe.state !== "pass").map((probe) => `${probe.id}: ${probe.detail}`).join(" | "));
});
