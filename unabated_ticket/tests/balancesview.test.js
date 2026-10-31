// Run: node --test unabated_ticket/tests/balancesview.test.js
//
// The Edges tab's Balances strip (balancesview.js), shared by the panel and
// the phone: which books show, in what order, and how a failed or missing
// read reads.
const test = require("node:test");
const assert = require("node:assert/strict");
const balancesView = require("../extension/balancesview.js");

const NOW = Date.parse("2026-10-31T20:00:00Z");
const FRESH = "2026-10-31T19:59:30Z";

function readingOf(cash, overrides) {
  return { cash, inBets: 10, available: cash, creditLimit: null, fetchedAt: FRESH, error: null, ...overrides };
}

test("rows follow Kalshi, Novig, Polymarket US with their sub-accounts; other books are left out", () => {
  const payload = { balances: {
    polymarket_us: readingOf(1402.55), novig_bot: readingOf(128.64, { inBets: null }), kalshi: readingOf(4210.18),
    wagerzon: readingOf(500), kalshi_bot: readingOf(765.25), novig: readingOf(11893.02), novig_place: readingOf(2720.76),
  } };
  const view = balancesView.viewOf(payload, null, NOW);
  assert.deepEqual(view.rows.map((row) => row.label),
    ["Kalshi", "Kalshi · Bot", "Novig", "Novig · Place", "Novig · Bot", "Polymarket US"]);
  assert.deepEqual(view.summary.map((item) => `${item.label} ${item.text}`), ["Kalshi $4,210", "Novig $11,893", "Poly $1,403"]);
  assert.equal(view.rows[4].inBetsText, "—");
  assert.equal(view.totalAvailableText, "$21,120");
  assert.deepEqual(view.notes, []);
});

test("a never-read book shows unavailable and stays out of the total", () => {
  const payload = { balances: { kalshi: readingOf(100), novig: readingOf(null, { error: "RuntimeError: login refused", fetchedAt: null }) } };
  const view = balancesView.viewOf(payload, null, NOW);
  assert.equal(view.rows[1].availableText, "unavailable");
  assert.equal(view.summary[1].text, "—");
  assert.equal(view.totalAvailableText, "$100");
  assert.match(view.notes.join("\n"), /Novig unavailable: RuntimeError: login refused/);
  assert.match(view.notes.join("\n"), /total leaves out/);
});

test("a failed or old read keeps its numbers and is flagged with its age", () => {
  const payload = { balances: {
    kalshi: readingOf(100, { error: "RuntimeError: 429", fetchedAt: "2026-10-31T19:54:00Z" }),
    novig: readingOf(50, { fetchedAt: "2026-10-31T19:40:00Z" }),
  } };
  const view = balancesView.viewOf(payload, null, NOW);
  assert.equal(view.rows[0].availableText, "$100.00");
  assert.ok(view.summary.every((item) => item.stale));
  assert.deepEqual(view.notes, ["Kalshi read 6 min ago (last read failed)", "Novig read 20 min ago"]);
});

test("available falls back to cash and a negative balance keeps its sign", () => {
  const payload = { balances: { kalshi: readingOf(-12.5, { available: null }) } };
  assert.equal(balancesView.viewOf(payload, null, NOW).rows[0].availableText, "-$12.50");
});

test("a failed fetch is named; before any load the strip has nothing to show", () => {
  assert.deepEqual(balancesView.viewOf(null, null, NOW), { loaded: false, summary: [], rows: [], totalAvailableText: "$0", notes: [] });
  assert.deepEqual(balancesView.viewOf(null, "HTTP 502", NOW).notes, ["Could not load balances: HTTP 502"]);
});
