// Run: node --test unabated_ticket/tests
// server/phone/phoneview.js (the phone page's words and numbers) and the
// unmatched-bets list the runner adds to /edges.json for the phone's Bets
// tab. The rows come from the real runner over the NFL slice
// (fixtures/v2_slice.json, as in runner.test.js), so the ticket is checked on
// the exact shape /edges.json serves.
const test = require("node:test");
const assert = require("node:assert/strict");
const fs = require("node:fs");
const path = require("node:path");
const { SNAPSHOT_BASE_URL } = require("../extension/scanner.js");
const edgeRows = require("../extension/edgerows.js");
const runnerLib = require("../server/runner.js");
const phoneView = require("../server/phone/phoneview.js");

const KICKOFF_MS = Date.parse("2026-09-13T17:00:00Z");
const NOW = KICKOFF_MS - 3600 * 1000;
const SERVICE_URL = "http://127.0.0.1:8094";
const BEARS_SPREAD_SOUTHPOINT = "289357360:ms99:si0:tid6";
const NOVIG_PANTHERS_ALT = "289357357:ms89:si1:tid5:alt-2.5";
const sliceBody = fs.readFileSync(path.join(__dirname, "fixtures", "v2_slice.json"), "utf8");
const OPEN_SETTINGS = { leagues: [1], bookMode: "all", bookIds: null, maxLineAgeHours: 1e6, includeAlts: true, minEdgePct: 1, groupByMarket: false };

function response(status, body, headers = {}) {
  return { ok: status < 300, status, headers: { get: (name) => headers[name.toLowerCase()] ?? null }, json: async () => JSON.parse(body) };
}

// The Bears -2.5 held at BetOnline (runner.test.js BEARS_HELD) and a bet on a team the board does not know.
const BEARS_HELD = {
  id: "bol-1", source: "betonline", venue: "betonline", league: "nfl", status: "open", sourceFetchedAt: "2026-09-13T15:00:00Z",
  awayTeam: "Chicago Bears", homeTeam: "Carolina Panthers", awayKey: null, homeKey: null,
  betType: "spread", period: "FG", side: "away", points: -2.5, price: -110, stake: 200, toWin: 181.82,
  placedAt: "2026-09-13T14:00:00Z", closedAt: null, eventStart: "2026-09-13T17:00:00Z", isParlayLeg: false, unmatchable: null, approx: [],
};
const UNKNOWN_TEAM = { ...BEARS_HELD, id: "bol-2", awayTeam: "Nowhere Ducks", homeTeam: "Carolina Panthers", placedAt: "2026-09-13T14:30:00Z", stake: 50 };
const A_FUTURE = { ...BEARS_HELD, id: "kalshi-9", venue: "kalshi", source: "kalshi", betType: "other", unmatchable: "future: season wins", stake: null, placedAt: "2026-09-12T10:00:00Z" };

async function payloadFor(settings, bets) {
  const fetchImpl = async (url) => (url.split("?")[0] === SNAPSHOT_BASE_URL(1)
    ? response(200, sliceBody, { "content-length": "1000", "last-modified": new Date(NOW - 20000).toUTCString() })
    : response(404, "{}"));
  const betsBody = { generatedAt: "2026-09-13T15:59:30Z", sources: { betonline: { ok: true, fetchedAt: "2026-09-13T15:59:00Z", error: null, count: bets.length } }, bets, crosswalk: [], pins: [], fillFairs: [] };
  const serviceFetch = async (url) => response(200, JSON.stringify(url.endsWith("/settings.json") ? { settings, updatedAt: null } : betsBody));
  const runner = runnerLib.createRunner({
    fetchImpl, serviceFetch, betsServiceUrl: SERVICE_URL, now: () => NOW,
    timers: { setInterval: () => 1, clearInterval: () => {} }, logInfo: () => {}, logWarning: () => {},
  });
  await runner.start();
  return { payload: runner.edgesPayload(), betsBody };
}

test("formatting: edge, price with cents off the exchange's own price, points, kickoff and line age", () => {
  assert.equal(phoneView.fmtEdgePct(5.3), "+5.30%");
  assert.equal(phoneView.fmtEdgePct(-0.5), "-0.50%");
  // Novig at 34.5¢ shows the exchange's number, not a rounded American round trip.
  assert.equal(phoneView.fmtPriceBoth(190, 4, 0.345), "+190 · 34.5¢");
  assert.equal(phoneView.fmtPriceBoth(-110, 1, -110), "-110 · 52.4¢");
  assert.deepEqual([phoneView.fmtPoints(-2.5), phoneView.fmtPoints(3), phoneView.fmtPoints(null)], ["-2.5", "+3", ""]);
  assert.equal(phoneView.fmtUntil(NOW + 45 * 60000, NOW), "in 45m");
  assert.equal(phoneView.fmtUntil(NOW + 125 * 60000, NOW), "in 2h 5m");
  assert.equal(phoneView.fmtUntil(NOW + 3 * 86400000, NOW), "in 3d");
  assert.equal(phoneView.fmtUntil(NOW - 60000, NOW), "started");
  assert.deepEqual([phoneView.untilLevel(NOW + 3600000, NOW), phoneView.untilLevel(NOW + 6 * 3600000, NOW), phoneView.untilLevel(NOW + 20 * 3600000, NOW)],
    ["soon urgent", "soon", ""]);
  assert.deepEqual([phoneView.fmtLineAge(null, NOW), phoneView.fmtLineAge(NOW - 30000, NOW), phoneView.fmtLineAge(NOW - 4 * 60000, NOW), phoneView.fmtLineAge(NOW - 3 * 86400000, NOW)],
    ["line age unknown", "line just changed", "line 4m old", "line 3d old"]);
  assert.equal(phoneView.fmtLiquidity(1234.5), "liq $1,235");
  assert.equal(phoneView.marketText({ betType: "Spread", period: "1H", isAlt: true, mainPoints: -3, rotation: 466 }), "Spread · 1H · alt of -3 · rot 466");
  assert.equal(phoneView.describeMatchup({ awayTeam: "Chicago Bears", homeTeam: "Carolina Panthers", leagueLabel: "NFL" }), "Chicago Bears @ Carolina Panthers · NFL");
});

test("the ticket for a held sportsbook line: add, its standalone size, no contracts, the payout of the add", async () => {
  const { payload } = await payloadFor(OPEN_SETTINGS, [BEARS_HELD]);
  const row = payload.items.find((item) => item.key === BEARS_SPREAD_SOUTHPOINT);
  const ticket = phoneView.ticketView(row, NOW);
  assert.equal(ticket.book, "SouthPoint");
  assert.equal(ticket.price, "-110 · 52.4¢");
  assert.equal(ticket.edge, "+5.30%");
  assert.equal(ticket.edgeTier, "hot");
  assert.match(ticket.start, / · in 1h 0m$/);
  const block = ticket.stakeBlock;
  assert.deepEqual([block.label, block.stake, block.noEdge, block.position, block.positionAgainst, block.contracts],
    ["Add to your position", "$237.25", false, "held $200 · $437.25 alone", false, null]);
  // -110: $237.25 wins $215.68 and pays $452.93.
  assert.deepEqual([block.toWin, block.payout], ["$215.68", "$452.93"]);
});

test("the ticket for an exchange line: contracts at the exchange's price, cost under the stake", async () => {
  const { payload } = await payloadFor(OPEN_SETTINGS, []);
  const row = payload.items.find((item) => item.key === NOVIG_PANTHERS_ALT);
  const block = phoneView.stakeBlockView(row);
  assert.equal(block.label, "Bet");
  assert.equal(block.stake, edgeRows.fmtDollars(row.advice.bet));
  // floor($215.13 / 0.345) = 623 contracts; 623 x 0.345 = 214.935 rounds (in float) to $214.93.
  assert.deepEqual(block.contracts, { text: "623 contracts @ 34.5¢", cost: "$214.93", under: false });
  assert.equal(phoneView.contractsView(0.2, row).text, "under 1 contract @ 34.5¢");
  assert.equal(phoneView.contractsView(0, row), null);
});

test("the stake block's label at size, with nothing resting, and an against-only position", () => {
  const row = { price: 150, sourceFormat: 1, sourcePrice: 150 };
  const atSize = phoneView.stakeBlockView({ ...row, advice: { kind: "sized", bet: 0, alone: 120, verb: "add", held: 300, against: 0, teasers: { held: 0, against: 0 }, cappedAt: null } });
  assert.deepEqual([atSize.label, atSize.stake, atSize.noEdge, atSize.toWin], ["Already at full size", "$0.00", true, null]);
  const empty = phoneView.stakeBlockView({ ...row, advice: { kind: "none", bet: 0, alone: 80, verb: "bet", held: 0, against: 0, teasers: { held: 0, against: 0 }, cappedAt: 0 } });
  assert.equal(empty.label, "Nothing resting at this price");
  const against = phoneView.stakeBlockView({ ...row, advice: { kind: "sized", bet: 140, alone: 100, verb: "bet", held: 0, against: 200, teasers: { held: 0, against: 0 }, cappedAt: null } });
  assert.deepEqual([against.position, against.positionAgainst, against.payout], ["against $200 · $100 alone", true, "$350.00"]);
});

test("the stake block shows the uncapped stake beside a liquidity cap", () => {
  const row = { price: 150, sourceFormat: 1, sourcePrice: 150 };
  const capped = phoneView.stakeBlockView({ ...row, advice: { kind: "none", bet: 17, alone: 240, verb: "bet", held: 0, against: 0, teasers: { held: 0, against: 0 }, cappedAt: 17, uncapped: 240 } });
  assert.equal(capped.position, "all $17 liq");
  assert.equal(capped.uncapped, "uncapped $240.00");
  const novig = phoneView.stakeBlockView({ ...row, sourceFormat: 4, sourcePrice: 0.319, book: { id: 89, name: "Novig" }, advice: { kind: "none", bet: 17, alone: 240, verb: "bet", held: 0, against: 0, teasers: { held: 0, against: 0 }, cappedAt: 17, uncapped: 240 } });
  assert.equal(novig.uncapped, "uncapped $240.00 · 752 contracts @ 31.9¢");
});

test("a Polymarket ticket shows contracts at Unabated's price", () => {
  const row = { price: 150, sourceFormat: 1, sourcePrice: 150, book: { id: 0, name: "Polymarket US" } };
  assert.deepEqual(phoneView.contractsView(100, row), { text: "250 contracts @ 40.0¢", cost: "$100.00", under: false });
});

test("the runner names the open bets no board game matches; the Bets tab groups them as the panel does", async () => {
  const { payload, betsBody } = await payloadFor(OPEN_SETTINGS, [BEARS_HELD, UNKNOWN_TEAM, A_FUTURE]);
  const unmatched = payload.betsService.unmatched;
  assert.ok(payload.betsService.boardLineCount > 0);
  assert.deepEqual(unmatched.map((entry) => entry.betId).sort(), ["bol-2", "kalshi-9"]);
  const ducks = unmatched.find((entry) => entry.betId === "bol-2");
  assert.deepEqual([ducks.needsGame, ducks.needsFix, ducks.attachable], [true, false, true]);
  assert.match(ducks.reason, /Nowhere Ducks/);
  assert.deepEqual(payload.books.live.map((book) => book.name), ["BetMGM", "Kalshi", "Novig", "SouthPoint", "Sports Interaction"]);

  const held = edgeRows.applyBetsPayload({ records: [], crosswalk: [], pins: [], fillFairs: [] }, betsBody, NOW);
  const tab = phoneView.betsTabView(held.records, unmatched, payload.betsService.boardLineCount, { sources: betsBody.sources }, NOW);
  assert.equal(tab.atRisk, "$250.00");
  assert.equal(tab.caption, "at risk · 3 open bets · 1 venue · 1 with no stake reported");
  assert.deepEqual(tab.needsGame.map((entry) => entry.betId), ["bol-2"]);
  assert.deepEqual(tab.needsFix, []);
  assert.deepEqual(tab.offBoard.map((entry) => entry.betId), ["kalshi-9"]);
  // Newest first; the unmatched one marked.
  assert.deepEqual(tab.open.map((entry) => [entry.item.id, entry.unmatched]), [["bol-2", true], ["bol-1", false], ["kalshi-9", true]]);
  assert.equal(tab.open[1].item.meta, "Betonline · NFL · Chicago Bears @ Carolina Panthers");
  assert.equal(tab.matchNote, "1 of 2 open game bets match a game on the board.");
  assert.equal(tab.venueRows.find((row) => row.venue === "betonline").level, "green");
  // Before the runner is read nothing is grouped as unmatched, and the note says it was not checked.
  const unread = phoneView.betsTabView(held.records, null, 0, null, NOW);
  assert.deepEqual([unread.needsGame.length, unread.offBoard.length], [0, 0]);
  assert.match(unread.matchNote, /not checked/);
});

test("banners: the page's failed reads, then what the runner reports, worst first", () => {
  const since = NOW - 6 * 60000;
  const edgesPoll = {
    okAt: NOW - 7 * 60000, error: "server runner at http://127.0.0.1:8095 unreachable (Connection refused)", errorAt: NOW, failingSince: since,
    payload: { scanner: { error: "league 2: HTTP 404" }, betsService: { error: "ECONNREFUSED", unreachableSince: since }, settings: { source: "defaults", error: "timeout" } },
  };
  const banners = phoneView.banners(edgesPoll, { okAt: NOW, error: null }, NOW);
  assert.deepEqual(banners.map((banner) => banner.level), ["bad", "bad", "bad", "warn"]);
  assert.match(banners[0].text, /^Edges unreachable since .* \(6 min\): server runner at http:\/\/127\.0\.0\.1:8095 unreachable/);
  assert.equal(banners[1].text, "Scanner: league 2: HTTP 404");
  assert.match(banners[2].text, /^Runner cannot read the bets service since .*: stakes ignore your open bets\. ECONNREFUSED$/);
  assert.equal(banners[3].text, "Runner is using its defaults settings: timeout");
  assert.deepEqual(phoneView.banners({ payload: null, okAt: null, error: null }, null, NOW), []);
  assert.equal(phoneView.freshnessText({ okAt: NOW - 4000 }, { okAt: null }, NOW), "edges 4s ago · bets —");
});

test("settings: defaults fill nulls, a save sends only what changed, books travel with their mode, reset sends null", () => {
  const nothing = { bankroll: null, multiplier: null, leagues: null, periods: null, betTypes: null, bookMode: null, bookIds: null, minEdgePct: null,
    minStake: null, maxLineAgeHours: null, minLiquidityToWin: null, includeAlts: null, sortBy: null, groupByMarket: null };
  const effective = phoneView.effectiveSettings(nothing);
  assert.deepEqual([effective.bankroll, effective.multiplier, effective.bookMode, effective.minEdgePct, effective.groupByMarket],
    [edgeRows.DEFAULT_STAKE_SETTINGS.bankroll, edgeRows.DEFAULT_STAKE_SETTINGS.multiplier, "default", edgeRows.DEFAULT_EDGE_SETTINGS.minEdgePct, true]);
  assert.deepEqual(Object.keys(phoneView.settingsDefaults()).sort(), Object.keys(nothing).sort());

  // The form as it loads: nothing to send.
  const untouched = { ...effective, leagues: phoneView.leaguesForSports(phoneView.sportsOfLeagues(effective.leagues), effective.leagues) };
  assert.deepEqual(phoneView.settingsUpdate(nothing, untouched), {});
  assert.deepEqual(phoneView.settingsUpdate(nothing, { ...untouched, bankroll: 8000, periods: [2, 1] }), { bankroll: 8000, periods: [2, 1] });
  // Periods in another order are the same list; undefined is "not given".
  assert.deepEqual(phoneView.settingsUpdate(nothing, { periods: [1], bankroll: undefined }), {});
  assert.deepEqual(phoneView.settingsUpdate(nothing, { bookMode: "custom", bookIds: [89, 105] }), { bookMode: "custom", bookIds: [89, 105] });
  const custom = { ...nothing, bookMode: "custom", bookIds: [89] };
  assert.deepEqual(phoneView.settingsUpdate(custom, { bookMode: "all" }), { bookMode: "all", bookIds: null });
  assert.deepEqual(phoneView.settingsUpdate(custom, { bookMode: "custom", bookIds: [89] }), {});
  assert.deepEqual(phoneView.resetUpdate("bankroll"), { bankroll: null });
  assert.deepEqual(phoneView.resetUpdate("bookMode"), { bookMode: null, bookIds: null });
  assert.equal(phoneView.defaultText("bankroll"), `default: ${edgeRows.DEFAULT_STAKE_SETTINGS.bankroll}`);
  assert.equal(phoneView.defaultText("includeAlts"), "default: off");
  assert.equal(phoneView.defaultText("periods"), "default: FG");
});

test("leagues follow the sports ticked, and a stored NFL-only list survives a save that left the sports alone", () => {
  assert.deepEqual(phoneView.sportsOfLeagues([1, 5]), ["football", "baseball"]);
  assert.deepEqual(phoneView.leaguesForSports(["football"], [1]), [1]);
  assert.deepEqual(phoneView.leaguesForSports(["football", "hockey"], [1]), [1, 2, 6, 11]);
  assert.deepEqual(phoneView.leaguesForSports([], [1]), []);
});
