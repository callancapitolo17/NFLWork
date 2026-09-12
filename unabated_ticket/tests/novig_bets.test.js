// Run: node --test unabated_ticket/tests
// novig_bets.js: the Novig content-script source's normaliser (#116).
//
// Fixture: fixtures/bets/novig_bets.json — the three portfolio responses the
// Novig web app fetches for its Portfolio screen, in the shape novig_page.js
// hands novig_content.js. Shapes come from the app bundle (see the fixture's
// "provenance"), not a live capture; the NFL rows sit on the NFL slice
// matchup (Chicago Bears @ Carolina Panthers, 2026-09-13T17:00Z) so
// bets.test.js can match them against real board lines.
const test = require("node:test");
const assert = require("node:assert/strict");
const fs = require("node:fs");
const path = require("node:path");
const novig = require("../extension/novig_bets.js");

const READ_AT = "2026-09-11T21:00:01.200Z";

function fixture() {
  return JSON.parse(fs.readFileSync(path.join(__dirname, "fixtures", "bets", "novig_bets.json"), "utf8"));
}

// Every response applied in order, as the content script does.
function collected() {
  const pages = new Map();
  let out = null;
  for (const response of fixture().responses) out = novig.applyResponse(pages, response);
  return out;
}

function records() {
  const out = collected();
  return novig.normalizeNovig({ orders: out.orders, parlays: out.parlays, readAt: READ_AT });
}

function byId(list, id) {
  const record = list.find((r) => r.id === id);
  assert.ok(record, `record ${id} missing`);
  return record;
}

function pick(record, keys) {
  const out = {};
  for (const key of keys) out[key] = record[key];
  return out;
}

const GAME_KEYS = ["league", "eventStart", "eventDate", "awayTeam", "homeTeam", "awayKey", "homeKey", "betType", "period", "side", "points", "price", "stake", "toWin", "contracts", "status", "approx", "unmatchable"];

test("applyResponse: three lists, complete when every last page is short of its limit", () => {
  const out = collected();
  assert.equal(out.orders.length, 15);
  assert.equal(out.parlays.length, 2);
  assert.equal(out.complete, true);
  assert.equal(novig.applyResponse(new Map(), { operationName: "WalletBalance_Query", variables: {}, data: { wallet: [] } }), null);
});

test("applyResponse: a full page leaves the read incomplete; offset 0 restarts its list; parlay lists key on their where-clause", () => {
  const pages = new Map();
  const rows = (n, prefix) => Array.from({ length: n }, (_, i) => ({ id: `${prefix}${i}`, market: {}, outcome: {}, fills: [] }));
  let out = novig.applyResponse(pages, { operationName: "ActivePortfolioOrders_Query", variables: { offset: 0, limit: 2 }, data: { ActivePortfolioOrders_Query: rows(2, "a") } });
  assert.equal(out.complete, false);
  out = novig.applyResponse(pages, { operationName: "ActivePortfolioOrders_Query", variables: { offset: 2, limit: 2 }, data: { ActivePortfolioOrders_Query: rows(1, "b") } });
  assert.deepEqual(out.orders.map((r) => r.id), ["a0", "a1", "b0"]);
  assert.equal(out.complete, false); // the settled and parlay lists have not been seen in this tab
  novig.applyResponse(pages, { operationName: "SettledPortfolioOrders_Query", variables: { offset: 0, limit: 2 }, data: { SettledPortfolioOrders_Query: [] } });
  out = novig.applyResponse(pages, { operationName: "ParlayPortfolioQuery", variables: { offset: 0, limit: 2, where: { status: { _eq: "FILLED" } } }, data: { parlay: [] } });
  assert.equal(out.complete, true);
  out = novig.applyResponse(pages, { operationName: "ActivePortfolioOrders_Query", variables: { offset: 0, limit: 2 }, data: { ActivePortfolioOrders_Query: rows(1, "c") } });
  assert.deepEqual(out.orders.map((r) => r.id), ["c0"]);
  novig.applyResponse(pages, { operationName: "ParlayPortfolioQuery", variables: { offset: 0, limit: 15, where: { status: { _eq: "FILLED" } } }, data: { parlay: [{ id: "p1", legs: [] }] } });
  out = novig.applyResponse(pages, { operationName: "ParlayPortfolioQuery", variables: { offset: 0, limit: 15, where: { status: { _in: ["WIN"] } } }, data: { parlay: [{ id: "p2", legs: [] }] } });
  assert.deepEqual(out.parlays.map((r) => r.id), ["p1", "p2"]);
  assert.equal(out.complete, true);
});

test("normalize: a matched moneyline bid — index 0 is the HOME team, price is a probability, stake = contracts x price", () => {
  const record = byId(records(), "novig:o-ml-car");
  assert.deepEqual(pick(record, GAME_KEYS), {
    league: "nfl", eventStart: "2026-09-13T17:00:00.000Z", eventDate: "2026-09-13",
    awayTeam: "Chicago Bears", homeTeam: "Carolina Panthers", awayKey: null, homeKey: null,
    betType: "moneyline", period: "FG", side: "home", points: null, price: -138, stake: 58, toWin: 42, contracts: 100,
    status: "open", approx: [], unmatchable: null,
  });
  assert.equal(record.source, "novig_page");
  assert.equal(record.venue, "novig");
  assert.equal(record.placedAt, "2026-09-11T15:05:12.481Z");
  assert.equal(record.closedAt, null);
  assert.equal(record.sourceFetchedAt, READ_AT);
  assert.equal(record.raw.marketTitle, "Chicago Bears @ Carolina Panthers · MONEY · CAR");
});

test("normalize: spread bid = the named team at the description's number; the LAY of the same outcome is the other team at the negated number and 1 - price", () => {
  const list = records();
  const bid = byId(list, "novig:o-sp-chi");
  assert.deepEqual(pick(bid, ["side", "points", "price", "stake", "contracts", "status"]), { side: "away", points: -13.5, price: 150, stake: 80, contracts: 200, status: "open" });
  const lay = byId(list, "novig:o-sp-lay");
  assert.deepEqual(pick(lay, ["side", "points", "price", "stake", "contracts", "status"]), { side: "home", points: 13.5, price: -150, stake: 30, contracts: 50, status: "open" });
  assert.equal(lay.raw.isBid, false);
});

test("normalize: totals — index 0 / 'Over' is over; a partial fill is sized on the matched part; a resting order is open on its full size with the caveat", () => {
  const list = records();
  const partial = byId(list, "novig:o-tot-over");
  assert.deepEqual(pick(partial, ["side", "points", "price", "stake", "contracts", "status", "approx"]), { side: "over", points: 47.5, price: -108, stake: 31.2, contracts: 60, status: "open", approx: [] });
  const resting = byId(list, "novig:o-tot-under-rest");
  assert.deepEqual(pick(resting, ["side", "points", "price", "stake", "contracts", "status", "approx"]), { side: "under", points: 47.5, price: 122, stake: 9, contracts: 20, status: "open", approx: [novig.APPROX_UNMATCHED] });
});

test("normalize: MLB *_1H markets are the F5 period; the team resolves by competitor symbol, the number by description", () => {
  const record = byId(records(), "novig:o-f5");
  assert.deepEqual(pick(record, ["league", "betType", "period", "side", "points", "awayTeam", "homeTeam"]), { league: "mlb", betType: "spread", period: "F5", side: "away", points: 0.5, awayTeam: "Miami Marlins", homeTeam: "Minnesota Twins" });
});

test("normalize: NCAAF is cfb; the away side by symbol even though its index says home would be 0", () => {
  const record = byId(records(), "novig:o-cfb");
  assert.deepEqual(pick(record, ["league", "side", "awayTeam", "homeTeam", "price", "stake"]), { league: "cfb", side: "away", awayTeam: "Chattanooga", homeTeam: "Eastern Kentucky", price: 138, stake: 10.5 });
});

test("normalize: a spread with no number in the description falls back to the market strike in the home perspective", () => {
  const fx = fixture();
  const order = fx.responses[0].data.ActivePortfolioOrders_Query.find((o) => o.id === "o-sp-chi");
  const awayNoNumber = JSON.parse(JSON.stringify(order));
  awayNoNumber.outcome.description = "CHI";
  assert.equal(novig.normalizeOrder(awayNoNumber, READ_AT).points, -13.5);
  const homeNoNumber = JSON.parse(JSON.stringify(order));
  homeNoNumber.outcome = { id: "x", index: 0, description: "CAR", status: "TBD", competitor: { symbol: "CAR" } };
  assert.equal(novig.normalizeOrder(homeNoNumber, READ_AT).points, 13.5);
  assert.equal(novig.normalizeOrder(homeNoNumber, READ_AT).side, "home");
});

test("normalize: props, unsupported leagues and blobs without teams fail closed with a reason, never a guessed game", () => {
  const list = records();
  assert.equal(byId(list, "novig:o-prop").unmatchable, "not a game market");
  assert.equal(byId(list, "novig:o-prop").betType, "other");
  assert.equal(byId(list, "novig:o-atp").unmatchable, "league not supported (ATP)");
  assert.equal(byId(list, "novig:o-noteams").unmatchable, "unreadable Novig order (no teams on event ev-x)");
  for (const id of ["novig:o-prop", "novig:o-atp", "novig:o-noteams"]) {
    const record = byId(list, id);
    assert.equal(record.league, null);
    assert.equal(record.status, "open");
  }
});

test("normalize: settled orders — a bid on WIN won, a lay on LOSS won, PUSH pushed, an unfilled cancel is void, a wash is closed", () => {
  const list = records();
  assert.deepEqual(pick(byId(list, "novig:o-won"), ["status", "closedAt", "side", "price", "stake"]), { status: "won", closedAt: "2026-09-07T21:00:00.000Z", side: "home", price: -150, stake: 30 });
  assert.deepEqual(pick(byId(list, "novig:o-lay-won"), ["status", "side", "points", "price", "stake"]), { status: "won", side: "home", points: -3.5, price: 122, stake: 18 });
  assert.equal(byId(list, "novig:o-push").status, "push");
  assert.deepEqual(pick(byId(list, "novig:o-cancel"), ["status", "stake", "contracts", "closedAt"]), { status: "void", stake: 0, contracts: 0, closedAt: "2026-09-06T12:30:00.000Z" });
  assert.equal(byId(list, "novig:o-wash").status, "closed");
});

test("normalize: a cancelled order with fills before settlement stays open on the matched part; a cash-out is closed; REJECTED is void", () => {
  const order = fixture().responses[0].data.ActivePortfolioOrders_Query.find((o) => o.id === "o-tot-over");
  const cancelled = { ...order, status: "CANCELED" };
  assert.deepEqual(pick(novig.normalizeOrder(cancelled, READ_AT), ["status", "contracts", "stake"]), { status: "open", contracts: 60, stake: 31.2 });
  const cashedOut = JSON.parse(JSON.stringify(order));
  cashedOut.status = "FILLED";
  cashedOut.qty = 0;
  cashedOut.market.cash_out_requests = [{ id: "c1", created_at: "2026-09-11T18:00:00+00:00", status: "APPROVED" }];
  assert.equal(novig.normalizeOrder(cashedOut, READ_AT).status, "closed");
  assert.equal(novig.normalizeOrder({ ...order, status: "REJECTED" }, READ_AT).status, "void");
  assert.deepEqual(novig.normalizeOrder({ ...order, status: "PENDING", qty: order.originalQty, fills: [] }, READ_AT).approx, [novig.APPROX_PENDING]);
});

test("normalize: parlays — one record per leg with the parlay's wager as the stake, status from the parlay", () => {
  const list = records();
  const leg0 = byId(list, "novig:pl-1:0");
  assert.deepEqual(pick(leg0, ["isParlayLeg", "parlayId", "legIndex", "legCount", "status", "side", "betType", "price", "stake", "toWin", "contracts"]),
    { isParlayLeg: true, parlayId: "novig:pl-1", legIndex: 0, legCount: 2, status: "open", side: "home", betType: "moneyline", price: -122, stake: 25, toWin: 75, contracts: null });
  const leg1 = byId(list, "novig:pl-1:1");
  assert.deepEqual(pick(leg1, ["betType", "side", "points", "status"]), { betType: "total", side: "under", points: 47.5, status: "open" });
  assert.equal(byId(list, "novig:pl-2:0").status, "lost");
  assert.equal(byId(list, "novig:pl-2:0").closedAt, "2026-09-07T21:00:00.000Z");
});

test("normalize: a row seen twice keeps the later page's copy; rows without an id are skipped", () => {
  const order = fixture().responses[0].data.ActivePortfolioOrders_Query.find((o) => o.id === "o-ml-car");
  const list = novig.normalizeNovig({ orders: [order, { ...order, qty: 100, status: "OPEN", fills: [] }, { market: {}, outcome: {} }], parlays: [], readAt: READ_AT });
  assert.equal(list.length, 1);
  assert.deepEqual(list[0].approx, [novig.APPROX_UNMATCHED]);
});

test("probabilityToAmerican: both sides of even money; outside (0, 1) is null", () => {
  assert.equal(novig.probabilityToAmerican(0.5), -100);
  assert.equal(novig.probabilityToAmerican(0.58), -138);
  assert.equal(novig.probabilityToAmerican(0.4), 150);
  assert.equal(novig.probabilityToAmerican(0), null);
  assert.equal(novig.probabilityToAmerican(1), null);
});
