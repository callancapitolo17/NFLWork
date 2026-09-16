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

// `page.window` is the sandbox: `posted` collects every window.postMessage
// (what content.js would receive), `deliver(type, payload)` plays a message
// from content.js into page.js's listener, `watchTimers()` counts watch intervals started.
function loadPage(pathname = "/nfl/odds", { gridRoots = [], rows = [], shells = [] } = {}) {
  const listeners = [];
  const posted = [];
  let timers = 0;
  const intervalCallbacks = [];
  const bySelector = { ".ag-root": gridRoots, ".ag-row": rows, ".odds-cell-action-shell": shells };
  const sandbox = {
    __unabatedTicketExposeInternals: true,
    console: { info() {}, warn() {}, error() {} },
    setInterval: (fn) => { timers += 1; intervalCallbacks.push(fn); return timers; },
    clearInterval() {},
    // The shape check retries while the grid has no rows; the harness never
    // runs those, it calls shapeCheckResult directly.
    setTimeout() { return 0; },
    Element: class {},
    location: { origin: "https://tools.unabated.com", pathname, href: `https://tools.unabated.com${pathname}` },
    postMessage(data) { posted.push(data); },
    addEventListener(type, fn) { if (type === "message") listeners.push(fn); },
    document: {
      addEventListener() {}, removeEventListener() {}, querySelector: () => null,
      querySelectorAll: (selector) => bySelector[selector] || [],
    },
  };
  sandbox.window = sandbox;
  vm.createContext(sandbox);
  vm.runInContext(fs.readFileSync(PAGE_JS, "utf8"), sandbox, { filename: "page.js" });
  assert.ok(sandbox.__unabatedTicketInternals, "page.js exposed nothing to the harness");
  const timersAtLoad = timers; // the heartbeat
  // page.js checks event.source against its own `window`, which inside the
  // context is the contextified global, not the outer sandbox object.
  const innerWindow = vm.runInContext("window", sandbox);
  const deliver = (type, payload) => listeners.forEach((fn) => fn({ source: innerWindow, data: { source: "unabated-ticket", type, payload } }));
  // The heartbeat is the one interval started at load.
  const heartbeat = () => intervalCallbacks[0]();
  return { ...sandbox.__unabatedTicketInternals, posted, deliver, heartbeat, watchTimers: () => timers - timersAtLoad };
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

// ---- load-time shape check ------------------------------------------------
//
// The one check page.js runs before you click anything: on a grid that has
// rendered rows, at least one price cell must still carry the React props a
// ticket is read from. Everything else in page.js reports itself at the point
// of use, so this stays a single count, not a probe matrix.

// A rendered price cell. `props` is what its fiber carries; the real cell's
// props are {marketLine, sideIndex, context} and only that shape is readable.
function cellShell(props) {
  return { "__reactFiber$test": { memoizedProps: props, return: null } };
}

const READABLE_CELL_PROPS = { marketLine: { points: -2.5, price: 182 }, sideIndex: 1, context: CONTEXT };

test("shape check: a grid with rows whose cells carry the ticket props passes", () => {
  const page = loadPage("/nfl/odds", { rows: [{}, {}, {}], shells: [cellShell(READABLE_CELL_PROPS)] });
  const result = page.shapeCheckResult();
  assert.equal(result.status, "ok");
  assert.equal(result.rows, 3);
  assert.equal(result.readable, 1);
});

test("shape check: a board with no rows is not checked at all, so nothing is claimed", () => {
  const page = loadPage("/nfl/odds", { rows: [], shells: [] });
  assert.equal(page.shapeCheckResult(), null);
});

test("shape check: rows but no price cells names the selector that moved", () => {
  const page = loadPage("/nfl/odds", { rows: [{}, {}], shells: [] });
  const result = page.shapeCheckResult();
  assert.equal(result.status, "changed");
  assert.match(result.message, /no price cells \(\.odds-cell-action-shell\) on a grid showing 2 rows/);
});

test("shape check: renaming the prop a ticket is read from fails the check and prints the shape found", () => {
  // Unabated's next bundle renames marketLine -> line: the cell still renders
  // and still has a fiber, it just no longer says what the ticket needs.
  const renamed = { line: { points: -2.5, price: 182 }, sideIndex: 1, context: CONTEXT };
  const page = loadPage("/nfl/odds", { rows: [{}, {}], shells: [cellShell(renamed), cellShell(renamed)] });
  const result = page.shapeCheckResult();
  assert.equal(result.status, "changed");
  assert.equal(result.shells, 2);
  assert.equal(result.readable, 0);
  assert.match(result.message, /2 price cells on 2 rows/);
  assert.match(result.detail, /line,sideIndex,context/, result.detail);
});

test("shape check: a price cell React no longer attaches a fiber to says so", () => {
  const page = loadPage("/nfl/odds", { rows: [{}], shells: [{}] });
  const result = page.shapeCheckResult();
  assert.equal(result.status, "changed");
  assert.match(result.detail, /no __reactFiber\$ key/);
});

test("shape check: page.js publishes 'checking' on load, so a previous load's verdict never stands", () => {
  const page = loadPage("/nfl/odds", { rows: [], shells: [] });
  const checks = page.posted.filter((message) => message.type === "pagecheck");
  assert.equal(checks.length, 1);
  assert.equal(checks[0].payload.status, "checking");
});

test("shape check: the verdict rides on every heartbeat, so a post content.js missed is not a silent miss", () => {
  // Extension reload re-injects page.js and content.js in two round trips; the
  // load-time post can land before content.js is listening.
  const renamed = { line: { points: -2.5 }, sideIndex: 1 };
  const page = loadPage("/nfl/odds", { rows: [{}], shells: [cellShell(renamed)] });
  page.posted.length = 0;
  page.heartbeat();
  const resent = page.posted.filter((message) => message.type === "pagecheck");
  assert.equal(resent.length, 1);
  assert.equal(resent[0].payload.status, "changed");
});

test("shape check: 'checking' is never re-sent, so a tab with no rows cannot keep clearing another tab's verdict", () => {
  const page = loadPage("/nfl/odds", { rows: [], shells: [] });
  page.posted.length = 0;
  page.heartbeat();
  assert.equal(page.posted.filter((message) => message.type === "pagecheck").length, 0);
});
