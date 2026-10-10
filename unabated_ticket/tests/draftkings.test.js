// Run: node --test unabated_ticket/tests
// draftkings.js: the token request, the My Bets socket messages, the paging rule,
// the trimmed push body and readAccount over a fake socket (2026-10-10). The bet
// shapes are pinned by the Python side (bets_service/tests/test_draftkings.py on
// fixtures/bets/draftkings_bets.json).
const test = require("node:test");
const assert = require("node:assert/strict");
const dk = require("../extension/draftkings.js");

const NOW = Date.parse("2026-10-10T19:00:00Z");
const DAY = 24 * 60 * 60 * 1000;

function bet(id, placedMs, extra = {}) {
  return {
    betId: id, receiptId: `r${id}`, type: "Single", status: "Unsettled", settlementStatus: "Open",
    numberOfBets: 1, numberOfSelections: 1, displayOdds: "+566", currency: "USD",
    placementDate: new Date(placedMs).toISOString(), stake: 75, potentialReturns: 499.5,
    liveActivity: { status: "Inactive" }, earlyExitStatus: "None",
    selections: [{
      selectionId: "s1", eventId: "34749110", marketId: "2_1", status: "Unsettled", settlementStatus: "Open",
      displayOdds: "+566", selectionDisplayName: "North Dakota State -14.5", marketDisplayName: "Spread Alternate",
      isInBetBoost: false, players: [], participants: [{ id: "139685", name: "North Dakota State", metadata: { rosettaTeamId: "9381" } }],
    }],
    combinations: [],
    ...extra,
  };
}

const EVENT = {
  sportId: "3", eventId: "34749110", eventStartDate: "2026-10-10T23:00:00.000Z", eventDisplayName: "North Dakota State @ UNLV",
  homeTeamName: "UNLV", awayTeamName: "North Dakota State", status: "NotStarted", leagueId: "87637", tags: ["SGP"],
  participants: [
    { id: "6310", name: "UNLV", venueRole: "Home", metadata: { retailRotNumber: "881" } },
    { id: "139685", name: "North Dakota State", venueRole: "Away", metadata: { retailRotNumber: "880" } },
  ],
  media: [{ providerName: "BetRadarV3" }],
};

test("the token request carries the session cookies and gives up instead of hanging", () => {
  const request = dk.jwtRequest();
  assert.equal(request.method, "GET");
  assert.equal(request.credentials, "include");
  assert.ok(request.signal instanceof AbortSignal);
  assert.equal(dk.JWT_URL, "https://gaming-us-ma.draftkings.com/api/wager/v1/generateEnterpriseJWT");
});

test("tokenOf: the token on 200; the fix on 401 / 403; a tokenless reply is an error", () => {
  assert.deepEqual(dk.tokenOf(200, { token: "abc", expiresIn: 600 }), { token: "abc" });
  assert.deepEqual(dk.tokenOf(401, null), { error: dk.NOT_LOGGED_IN });
  assert.deepEqual(dk.tokenOf(403, null), { error: dk.NOT_LOGGED_IN });
  assert.equal(dk.tokenOf(500, null).error, "DraftKings token failed (HTTP 500)");
  assert.equal(dk.tokenOf(200, {}).error, "DraftKings token reply carries no token");
});

test("socketUrl puts the token in the query, encoded", () => {
  assert.equal(dk.socketUrl("a.b+c"),
    "wss://gateway.northamerica-northeast2.prod.dkapis.com/dkusma/shelby/api/v1/websocket?format=json&jwt=a.b%2Bc");
});

test("betsMessage is the page's own BetsRequest, newest first, 25 at a time", () => {
  assert.deepEqual(dk.betsMessage("id-1", "Settled", 50), {
    jsonrpc: "2.0", method: "BetsRequest", id: "id-1",
    params: {
      filter: { status: "Settled" }, orderCriteria: { orderBy: "placementDate", direction: "DESC" },
      pagination: { count: 25, skip: 50 }, locale: "en", ScoreboardType: "EventScore",
    },
  });
  assert.equal(dk.initializeMessage("id-0").method, "InitializeBetsPageRequest");
});

test("pageOf: a BetsRequest reply, the initial page's nesting, a refusal and a bad reply", () => {
  assert.deepEqual(dk.pageOf({ result: { bets: [], events: {} }, id: "x" }), { bets: [], events: {} });
  assert.deepEqual(dk.pageOf({ result: { initial: { bets: [1], events: { e: {} } } } }), { bets: [1], events: { e: {} } });
  assert.match(dk.pageOf({ error: { code: 401, message: "Unauthorized" }, id: "x" }).error,
    /^DraftKings bets request refused: \{"code":401/);
  assert.equal(dk.pageOf({ result: {} }).error, "DraftKings bets reply carries no bets list");
  assert.equal(dk.pageOf(null).error, "DraftKings bets reply is not JSON");
});

test("isLastPage: a short page ends a list; a settled page past the window ends it too", () => {
  const full = (placedMs) => Array.from({ length: 25 }, (_, index) => bet(`b${index}`, placedMs));
  assert.equal(dk.isLastPage("Open", [bet("1", NOW)], NOW), true);
  assert.equal(dk.isLastPage("Open", full(NOW - 90 * DAY), NOW), false);
  assert.equal(dk.isLastPage("Settled", full(NOW - DAY), NOW), false);
  assert.equal(dk.isLastPage("Settled", full(NOW - 32 * DAY), NOW), true);
});

test("pushBody keeps what the service reads and drops the rest", () => {
  const body = dk.pushBody("2026-10-10T19:00:00.000Z", { Open: [bet("1", NOW, { bonus: { freeBetAmount: "25" } })], Settled: [] },
    { 34749110: EVENT });
  assert.deepEqual(Object.keys(body), ["fetchedAt", "open", "settled", "events"]);
  const [open] = body.open;
  assert.equal(open.betId, "1");
  assert.equal(open.freeBetAmount, 25);
  assert.equal(open.combinationCount, 0);
  assert.equal(open.liveActivity, undefined);
  assert.deepEqual(open.selections[0].participants, [{ id: "139685", name: "North Dakota State", venueRole: undefined }]);
  assert.equal(open.selections[0].players, undefined);
  const event = body.events["34749110"];
  assert.equal(event.homeTeamName, "UNLV");
  assert.equal(event.media, undefined);
  assert.deepEqual(event.participants[1], { id: "139685", name: "North Dakota State", venueRole: "Away" });
  assert.throws(() => dk.pushBody("t", { Open: [] }, {}), /Settled missing/);
});

test("an SGP group inside a parlay is flagged by its nested count", () => {
  const nested = bet("1", NOW);
  nested.selections[0].nestedSGPSelections = [{}, {}];
  assert.equal(dk.trimBet(nested).selections[0].nestedSelectionCount, 2);
});

// A socket that answers each request from `answer(message)`, after the open.
function fakeSocketClass(answer, log) {
  return class FakeSocket {
    constructor(url) {
      this.url = url;
      log.push({ url });
      setTimeout(() => this.onopen(), 0);
    }
    send(text) {
      const message = JSON.parse(text);
      log.push(message);
      setTimeout(() => {
        this.onmessage({ data: JSON.stringify({ result: { cashOutUpdate: { betId: "x" } }, id: message.id }) });
        this.onmessage({ data: JSON.stringify(answer(message)) });
      }, 0);
    }
    close() {
      log.push("closed");
    }
  };
}

function fakeFetch(status, body) {
  return async () => ({ status, json: async () => body });
}

test("readAccount: the token, the first message, every page of both lists, one push; the socket closes", async () => {
  const log = [];
  const openBets = Array.from({ length: 30 }, (_, index) => bet(`o${index}`, NOW - index * 1000));
  const answer = (message) => {
    if (message.method === "InitializeBetsPageRequest") return { result: { initial: { bets: [], events: {} } }, id: message.id };
    const { status } = message.params.filter;
    const { skip, count } = message.params.pagination;
    const list = status === "Open" ? openBets : [bet("s1", NOW - DAY, { settlementStatus: "Won" })];
    return { result: { bets: list.slice(skip, skip + count), events: { 34749110: EVENT } }, id: message.id };
  };
  const push = await dk.readAccount({ fetchImpl: fakeFetch(200, { token: "tok" }), WebSocketImpl: fakeSocketClass(answer, log), now: () => NOW });
  assert.equal(push.open.length, 30);
  assert.deepEqual(push.settled.map((item) => item.betId), ["s1"]);
  assert.deepEqual(Object.keys(push.events), ["34749110"]);
  assert.equal(push.fetchedAt, "2026-10-10T19:00:00.000Z");
  assert.match(log[0].url, /jwt=tok$/);
  assert.deepEqual(log.filter((entry) => entry.method === "BetsRequest").map((entry) => [entry.params.filter.status, entry.params.pagination.skip]),
    [["Open", 0], ["Open", 25], ["Settled", 0]]);
  assert.equal(log[log.length - 1], "closed");
});

test("readAccount: a logged-out token stops before the socket; a refused page throws and still closes it", async () => {
  const log = [];
  await assert.rejects(dk.readAccount({ fetchImpl: fakeFetch(401, null), WebSocketImpl: fakeSocketClass(() => ({}), log) }),
    new RegExp(dk.NOT_LOGGED_IN));
  assert.equal(log.length, 0);
  const refuse = (message) => (message.method === "InitializeBetsPageRequest"
    ? { result: { initial: { bets: [], events: {} } }, id: message.id }
    : { error: { message: "expired" }, id: message.id });
  await assert.rejects(dk.readAccount({ fetchImpl: fakeFetch(200, { token: "tok" }), WebSocketImpl: fakeSocketClass(refuse, log), now: () => NOW }),
    /DraftKings bets request refused/);
  assert.equal(log[log.length - 1], "closed");
});

test("isDue: first run at once, then every five minutes", () => {
  assert.equal(dk.isDue(0, NOW), true);
  assert.equal(dk.isDue(NOW, NOW + 60e3), false);
  assert.equal(dk.isDue(NOW, NOW + dk.POLL_MS), true);
});
