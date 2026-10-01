// Run: node --test unabated_ticket/tests
// bet105.js: the requests the panel makes to app.bet105.ag and the push body
// it sends the bets service (2026-09-29). The wire shapes are pinned by the
// Python side (bets_service/tests/test_bet105.py on fixtures/bets/bet105_history.json).
const test = require("node:test");
const assert = require("node:assert/strict");
const bet105 = require("../extension/bet105.js");

test("historyUrl and the two requests: credentials go along, the CSRF header rides the POST", () => {
  assert.equal(bet105.historyUrl("prematch"), "https://app.bet105.ag/__bff/__partner-prematch/betLobbyV2/logic/");
  assert.equal(bet105.historyUrl("live"), "https://app.bet105.ag/__bff/__partner-live/betLobbyV2/logic/");
  assert.equal(bet105.CUSTOMERS_URL, "https://app.bet105.ag/__bff/api/customers");
  assert.deepEqual(bet105.FEEDS, ["prematch", "live"]);
  const customers = bet105.customersRequest();
  assert.equal(customers.method, "GET");
  assert.equal(customers.credentials, "include");
  const history = bet105.historyRequest("tok-1");
  assert.equal(history.method, "POST");
  assert.equal(history.credentials, "include");
  assert.equal(history.headers["X-Broker-CSRF"], "tok-1");
  assert.deepEqual(JSON.parse(history.body), { a: "getHistory", state: "0" });
});

test("csrfTokenOf: the token on 200; the fix on 401; Cloudflare on 403; a tokenless body is an error", () => {
  assert.deepEqual(bet105.csrfTokenOf(200, { csrfToken: "tok-1", customerData: {} }), { csrfToken: "tok-1" });
  assert.deepEqual(bet105.csrfTokenOf(401, null), { error: bet105.NOT_LOGGED_IN });
  assert.match(bet105.csrfTokenOf(403, null).error, /^Bet105 session check refused \(HTTP 403\) — Cloudflare/);
  assert.equal(bet105.csrfTokenOf(500, null).error, "Bet105 session check failed (HTTP 500)");
  assert.equal(bet105.csrfTokenOf(200, { customerData: {} }).error, "Bet105 session check answered without a csrfToken");
  assert.equal(bet105.csrfTokenOf(200, null).error, "Bet105 session check answered without a csrfToken");
});

test("betGroupsOf: a list, a dict of groups, an empty list; site errors and bad replies are errors", () => {
  assert.deepEqual(bet105.betGroupsOf("prematch", 200, { betGroups: [{ betGroupId: 1 }] }), { betGroups: [{ betGroupId: 1 }] });
  assert.deepEqual(bet105.betGroupsOf("prematch", 200, { betGroups: { 1: { betGroupId: 1 } } }), { betGroups: [{ betGroupId: 1 }] });
  assert.deepEqual(bet105.betGroupsOf("live", 200, { betGroups: [] }), { betGroups: [] });
  assert.equal(bet105.betGroupsOf("live", 200, { e: "Session expired" }).error, "Bet105 live history: Session expired");
  assert.equal(bet105.betGroupsOf("live", 200, {}).error, "Bet105 live history has no betGroups");
  assert.equal(bet105.betGroupsOf("live", 200, null).error, "Bet105 live history is not JSON");
  assert.equal(bet105.betGroupsOf("prematch", 401, null).error, bet105.NOT_LOGGED_IN);
  assert.equal(bet105.betGroupsOf("prematch", 502, null).error, "Bet105 prematch history failed (HTTP 502)");
});

test("pushBody carries both feeds or throws; errorBody carries the text", () => {
  const body = bet105.pushBody("2026-09-29T04:00:00Z", { prematch: [{ betGroupId: 1 }], live: [] });
  assert.deepEqual(body, { fetchedAt: "2026-09-29T04:00:00Z", feeds: { prematch: [{ betGroupId: 1 }], live: [] } });
  assert.throws(() => bet105.pushBody("2026-09-29T04:00:00Z", { prematch: [] }), /live missing/);
  assert.deepEqual(bet105.errorBody("not logged in"), { error: "not logged in" });
  assert.deepEqual(bet105.errorBody(""), { error: "unknown error" });
});

test("isDue: at once when never run, then every POLL_MS", () => {
  assert.equal(bet105.POLL_MS, 5 * 60 * 1000);
  assert.equal(bet105.isDue(0, 1000), true);
  assert.equal(bet105.isDue(1000, 1000 + bet105.POLL_MS - 1), false);
  assert.equal(bet105.isDue(1000, 1000 + bet105.POLL_MS), true);
});
