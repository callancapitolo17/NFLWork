// Run: node --test unabated_ticket/tests/balancesview.test.js
//
// The Edges tab's Balances line (balancesview.js), shared by the panel and the
// phone: which books show, in what order, and how a failed or missing read reads.
const test = require("node:test");
const assert = require("node:assert/strict");
const balancesView = require("../extension/balancesview.js");

const NOW = Date.parse("2026-10-31T20:00:00Z");
const FRESH = "2026-10-31T19:59:30Z";

function readingOf(cash, overrides) {
  return { cash, inBets: 10, available: cash, creditLimit: null, fetchedAt: FRESH, error: null, ...overrides };
}

function textOf(view) {
  return view.items.map((item) => `${item.label} ${item.text}`);
}

test("Kalshi, Novig, Place, Poly in that order; bots and other books left out", () => {
  const payload = { balances: {
    polymarket_us: readingOf(980.3), novig_bot: readingOf(157.36), kalshi: readingOf(4872.23),
    wagerzon: readingOf(500), kalshi_bot: readingOf(955.58), novig: readingOf(24.29), novig_place: readingOf(1018.23),
  } };
  const view = balancesView.viewOf(payload, null, NOW);
  assert.deepEqual(textOf(view), ["Kalshi $4,872", "Novig $24", "Place $1,018", "Poly $980"]);
  assert.ok(view.items.every((item) => !item.stale));
});

test("a never-read book shows a dash in amber with the reason", () => {
  const payload = { balances: { kalshi: readingOf(null, { error: "RuntimeError: login refused", fetchedAt: null }) } };
  const [item] = balancesView.viewOf(payload, null, NOW).items;
  assert.equal(item.text, "—");
  assert.equal(item.stale, true);
  assert.equal(item.why, "RuntimeError: login refused");
});

test("a failed or old read keeps its number in amber and says how old it is", () => {
  const payload = { balances: {
    kalshi: readingOf(100, { error: "RuntimeError: 429", fetchedAt: "2026-10-31T19:54:00Z" }),
    novig: readingOf(50, { fetchedAt: "2026-10-31T19:40:00Z" }),
  } };
  const items = balancesView.viewOf(payload, null, NOW).items;
  assert.deepEqual(items.map((item) => [item.text, item.stale, item.why]), [
    ["$100", true, "last read failed (RuntimeError: 429); this is from 6 min ago"],
    ["$50", true, "read 20 min ago"],
  ]);
});

test("available falls back to cash and a negative balance keeps its sign", () => {
  const payload = { balances: { kalshi: readingOf(-12.5, { available: null }) } };
  assert.equal(balancesView.viewOf(payload, null, NOW).items[0].text, "-$13");
});

test("before any load the line has nothing; a failed fetch is carried", () => {
  assert.deepEqual(balancesView.viewOf(null, null, NOW), { loaded: false, items: [], error: null });
  assert.equal(balancesView.viewOf(null, "HTTP 502", NOW).error, "HTTP 502");
});
