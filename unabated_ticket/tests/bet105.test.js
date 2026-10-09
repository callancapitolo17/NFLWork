// Run: node --test unabated_ticket/tests
// bet105.js: the requests the panel makes to app.bet105.ag and the push body
// it sends the bets service (open 2026-09-29, settled 2026-10-07). The wire shapes are pinned by the
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
  assert.ok(customers.signal instanceof AbortSignal, "the session check gives up instead of hanging the poll");
  assert.ok(history.signal instanceof AbortSignal, "getHistory gives up instead of hanging the poll");
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

test("settledRequest: the My Bets page's own wagers/search POST, CSRF header and an empty body", () => {
  assert.equal(bet105.SETTLED_URL, "https://app.bet105.ag/__bff/api/wagers/search");
  const settled = bet105.settledRequest("tok-1");
  assert.equal(settled.method, "POST");
  assert.equal(settled.credentials, "include");
  assert.equal(settled.headers["X-Broker-CSRF"], "tok-1");
  assert.equal(settled.headers["Content-Type"], "application/json");
  assert.deepEqual(JSON.parse(settled.body), {});
  assert.ok(settled.signal instanceof AbortSignal, "wagers/search gives up instead of hanging the poll");
});

// One wager as wagers/search sends it (the 2026-10-07 capture's shape; ids and
// amounts synthetic, the account fields made up).
function capturedWager() {
  return {
    risk: 117, isPoS: false, toWin: 100, xRate: 1, result: 100, agentId: 1, details: {}, gradeId: 1, usdRisk: 117,
    wagerId: 81000002, category: "SPORTS", userName: "someone", agentName: "an agent", gradeTime: "2026-10-06T01:56:19Z",
    isCashout: false, maxPayout: 217, placeTime: "2026-09-28T02:01:09Z", taxAmount: 0, usdResult: 100, customerId: 1,
    isFreePlay: false, updateTime: "2026-10-06T01:56:44Z", isAnonymous: 0, isCashBonus: false, productCode: "PreMatch",
    wagerStatus: "Win", currencyCode: "USD", ticketNumber: "91000002", adminUserName: "bet105-wager-bot",
    wagerDetails: "Atlanta Falcons vs New Orleans Saints / 1st Half / Total / Over 22.5 -117",
    dailyFigureTime: "2026-10-06T01:56:19Z", preGradeBalance: 1, preWagerBalance: 2, postGradeBalance: 3, postWagerBalance: 4,
    properties: {
      odds: 1.8547, fmtOdds: "-117", grades: ["W"], rrSizeMap: 1, teaserName: null, featuredBetOdds: null,
      fixedParlayName: null, fixedParlayOdds: null, balances: { preGrade: 1, preWager: 2, postGrade: 3, postWager: 4 },
      legs: [{
        odds: 1.8547, side: 1, legId: 71000002, sport: "Football", team1: "Atlanta Falcons", team2: "New Orleans Saints",
        figure: 22.5, league: "NFL", market: "Total", period: "1st Half", eventId: 194318121, fmtOdds: "-117",
        isIfLeg: false, premium: 0, sportId: 3, leagueId: 4, marketId: 5, periodId: "h1", pitcher1: null, pitcher2: null,
        startTime: "2026-10-06 00:15:00+00:00", contestant: "", disclaimer: "", team1Short: "ATL Falcons",
        team2Short: "NO Saints", wagerSelection: "Over 22.5", figureDescription: "22.5",
      }],
    },
  };
}

test("settledWagerOf keeps what the service reads, in the venue's nesting, and drops balances and account names", () => {
  const sent = bet105.settledWagerOf(capturedWager());
  assert.deepEqual(sent, {
    wagerId: 81000002, ticketNumber: "91000002", productCode: "PreMatch", category: "SPORTS", wagerStatus: "Win",
    placeTime: "2026-09-28T02:01:09Z", gradeTime: "2026-10-06T01:56:19Z", risk: 117, toWin: 100, result: 100,
    isFreePlay: false, isCashout: false,
    wagerDetails: "Atlanta Falcons vs New Orleans Saints / 1st Half / Total / Over 22.5 -117",
    properties: {
      odds: 1.8547, fmtOdds: "-117", grades: ["W"], teaserName: null, fixedParlayName: null,
      legs: [{
        legId: 71000002, eventId: 194318121, sportId: 3, leagueId: 4, league: "NFL", team1: "Atlanta Falcons",
        team2: "New Orleans Saints", periodId: "h1", period: "1st Half", marketId: 5, market: "Total", side: 1,
        figure: 22.5, fmtOdds: "-117", startTime: "2026-10-06 00:15:00+00:00",
      }],
    },
  });
  // A Pending wager has no gradeTime: left out, not sent as null.
  const pending = capturedWager();
  delete pending.gradeTime;
  assert.equal("gradeTime" in bet105.settledWagerOf(pending), false);
  assert.deepEqual(bet105.settledWagerOf({ ticketNumber: "1" }), { ticketNumber: "1", properties: { legs: [] } });
});

test("wagersOf: the list cut to the read fields; the site's refusal code and bad replies are errors", () => {
  const result = bet105.wagersOf(200, [capturedWager()]);
  assert.equal(result.wagers.length, 1);
  assert.equal(result.wagers[0].userName, undefined);
  assert.deepEqual(bet105.wagersOf(200, []), { wagers: [] });
  assert.equal(bet105.wagersOf(403, { code: "CSRF_FAILED" }).error, "Bet105 settled history refused (HTTP 403, CSRF_FAILED)");
  assert.match(bet105.wagersOf(403, null).error, /^Bet105 settled history refused \(HTTP 403\) — Cloudflare/);
  assert.equal(bet105.wagersOf(401, { code: "UNAUTHORIZED" }).error, bet105.NOT_LOGGED_IN);
  assert.equal(bet105.wagersOf(500, null).error, "Bet105 settled history failed (HTTP 500)");
  assert.equal(bet105.wagersOf(200, { betGroups: [] }).error, "Bet105 settled history is not a list");
  assert.equal(bet105.wagersOf(200, null).error, "Bet105 settled history is not a list");
});

test("pushBody carries both feeds and the settled list or throws; errorBody carries the text", () => {
  const body = bet105.pushBody("2026-09-29T04:00:00Z", { prematch: [{ betGroupId: 1 }], live: [] }, [{ ticketNumber: "1" }]);
  assert.deepEqual(body, {
    fetchedAt: "2026-09-29T04:00:00Z", feeds: { prematch: [{ betGroupId: 1 }], live: [] }, settled: [{ ticketNumber: "1" }],
  });
  assert.throws(() => bet105.pushBody("2026-09-29T04:00:00Z", { prematch: [] }, []), /live missing/);
  assert.throws(() => bet105.pushBody("2026-09-29T04:00:00Z", { prematch: [], live: [] }), /settled list/);
  assert.deepEqual(bet105.errorBody("not logged in"), { error: "not logged in" });
  assert.deepEqual(bet105.errorBody(""), { error: "unknown error" });
});

test("isDue: at once when never run, then every POLL_MS", () => {
  assert.equal(bet105.POLL_MS, 5 * 60 * 1000);
  assert.equal(bet105.isDue(0, 1000), true);
  assert.equal(bet105.isDue(1000, 1000 + bet105.POLL_MS - 1), false);
  assert.equal(bet105.isDue(1000, 1000 + bet105.POLL_MS), true);
});
