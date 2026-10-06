// Run: node --test unabated_ticket/tests
// server/runner.js and server/edges_payload.js: the headless Edges scan for
// the phone page (plan step 1). The scanner is fed the real NFL slice
// (fixtures/v2_slice.json) through an injected fetch, the bets service is a
// scripted fetch, and the clock sits an hour before the slice's kickoff — no
// live network anywhere.
const test = require("node:test");
const assert = require("node:assert/strict");
const fs = require("node:fs");
const path = require("node:path");
const { SNAPSHOT_BASE_URL } = require("../extension/scanner.js");
const edgeRows = require("../extension/edgerows.js");
const tailflex = require("../extension/tailflex.js");
const teaser = require("../extension/teaser.js");
const betsLib = require("../extension/bets.js");
const runnerLib = require("../server/runner.js");
const edgesPayload = require("../server/edges_payload.js");

const KICKOFF_MS = Date.parse("2026-09-13T17:00:00Z");
const NOW = KICKOFF_MS - 3600 * 1000;
const SERVICE_URL = "http://127.0.0.1:8094";
const BEARS_SPREAD_SOUTHPOINT = "289357360:ms99:si0:tid6";
const noTimers = { setInterval: () => 1, clearInterval: () => {} };
const quiet = { logInfo: () => {}, logWarning: () => {} };
const sliceBody = fs.readFileSync(path.join(__dirname, "fixtures", "v2_slice.json"), "utf8");

function response({ status = 200, body = "", headers = {} }) {
  return {
    ok: status >= 200 && status < 300,
    status,
    headers: { get: (name) => headers[name.toLowerCase()] ?? null },
    json: async () => JSON.parse(body),
  };
}

// The scanner always loads NFL and CFB (the panel's Teasers leagues) on top of the settings' leagues.
const NFL_AND_CFB = [SNAPSHOT_BASE_URL(1), SNAPSHOT_BASE_URL(2)];

// Unabated's CDN: the NFL slice for league 1, a 404 for anything else.
function unabatedFetch() {
  const calls = [];
  const impl = async (url) => {
    const base = url.split("?")[0];
    calls.push(base);
    if (base === SNAPSHOT_BASE_URL(1)) return response({ body: sliceBody, headers: { "content-length": "1000", "last-modified": new Date(NOW - 20000).toUTCString() } });
    return response({ status: 404, body: "not found" });
  };
  impl.calls = calls;
  return impl;
}

// An open BetOnline bet on the Bears -2.5 at -110 for $200, as /bets.json serves it.
const BEARS_HELD = {
  id: "bol-1", source: "betonline", venue: "betonline", league: "nfl", status: "open", sourceFetchedAt: "2026-09-13T15:00:00Z",
  awayTeam: "Chicago Bears", homeTeam: "Carolina Panthers", awayKey: null, homeKey: null,
  betType: "spread", period: "FG", side: "away", points: -2.5, price: -110, stake: 200, toWin: 181.82,
  placedAt: "2026-09-13T14:00:00Z", closedAt: null, eventStart: "2026-09-13T17:00:00Z", isParlayLeg: false, unmatchable: null, approx: [],
};

// An open $400 BFA 2-team teaser: Panthers +3.5 (rotation 466) on the slice's
// game, and a leg on a game that started elsewhere, which counts as won.
function teaserLeg(legIndex, fields) {
  return {
    id: `bfa:T1:leg${legIndex}`, source: "bfa_api", venue: "bfa", status: "open", league: "nfl", betType: "spread",
    period: "FG", side: "home", points: 3.5, price: -110, rotation: 466, awayTeam: null, homeTeam: "Carolina Panthers", awayKey: null, homeKey: null,
    eventStart: new Date(KICKOFF_MS).toISOString(), stake: 400, toWin: 1200, isParlayLeg: true, parlayId: "bfa:T1", legIndex,
    legCount: 2, placedAt: "2026-09-13T12:00:00Z", approx: ["side_from_rotation_parity"], unmatchable: null,
    raw: { headerDescription: "2 TEAM TEASERS" }, ...fields,
  };
}
const PANTHERS_TEASER = [teaserLeg(0), teaserLeg(1, { rotation: 901, homeTeam: "Elsewhere", eventStart: new Date(NOW - 2 * 3600 * 1000).toISOString() })];

// The slice's lines are days old by NOW, so the line-age gate is opened; every
// live book, alts on, a 1% floor so the slice's thinner lines list.
const OPEN_SETTINGS = { leagues: [1], bookMode: "all", bookIds: null, maxLineAgeHours: 1e6, includeAlts: true, minEdgePct: 1 };

// The bets service: scripted /bets.json and /settings.json bodies, or a thrown error.
function serviceFetch(script) {
  const calls = [];
  const impl = async (url) => {
    calls.push(url);
    const route = url.slice(SERVICE_URL.length);
    const answer = script[route];
    if (answer instanceof Error) throw answer;
    if (answer === undefined) return response({ status: 404, body: "{}" });
    return response({ body: JSON.stringify(typeof answer === "function" ? answer() : answer) });
  };
  impl.calls = calls;
  return impl;
}

function betsBody(bets) {
  return { generatedAt: "2026-09-13T15:59:30Z", sources: { betonline: { ok: true, fetchedAt: "2026-09-13T15:59:00Z", error: null, count: bets.length } }, bets, crosswalk: [], pins: [], fillFairs: [] };
}

async function startedRunner(script) {
  const fetchImpl = unabatedFetch();
  const runner = runnerLib.createRunner({ fetchImpl, serviceFetch: serviceFetch(script), betsServiceUrl: SERVICE_URL, now: () => NOW, timers: noTimers, ...quiet });
  await runner.start();
  return { runner, fetchImpl };
}

test("configFromEnv: loopback 8095 and the local bets service by default; a bad port or URL fails naming the variable", () => {
  assert.deepEqual(runnerLib.configFromEnv({}), { host: "127.0.0.1", port: 8095, betsServiceUrl: "http://127.0.0.1:8094" });
  assert.deepEqual(runnerLib.configFromEnv({ UNABATED_RUNNER_HOST: "100.64.0.7", UNABATED_RUNNER_PORT: "9000", BETS_SERVICE_URL: "http://127.0.0.1:8094/" }),
    { host: "100.64.0.7", port: 9000, betsServiceUrl: "http://127.0.0.1:8094" });
  assert.throws(() => runnerLib.configFromEnv({ UNABATED_RUNNER_PORT: "80a" }), /UNABATED_RUNNER_PORT must be a whole number from 1 to 65535, got "80a"/);
  assert.throws(() => runnerLib.configFromEnv({ BETS_SERVICE_URL: "ftp://x" }), /BETS_SERVICE_URL must be an http\(s\) URL/);
});

test("allowedHosts: the loopback names and the bound address, with the port; bare names only on 80", () => {
  assert.deepEqual(runnerLib.allowedHosts("127.0.0.1", 8095), ["127.0.0.1:8095", "localhost:8095"]);
  assert.deepEqual(runnerLib.allowedHosts("100.64.0.7", 8095), ["127.0.0.1:8095", "localhost:8095", "100.64.0.7:8095"]);
  assert.deepEqual(runnerLib.allowedHosts("0.0.0.0", 80), ["127.0.0.1:80", "localhost:80", "127.0.0.1", "localhost"]);
  const allowed = runnerLib.allowedHosts("127.0.0.1", 8095);
  assert.equal(runnerLib.hostAllowed("LOCALHOST:8095", allowed), true);
  assert.equal(runnerLib.hostAllowed("evil.example:8095", allowed), false);
  assert.equal(runnerLib.hostAllowed("127.0.0.1:8094", allowed), false);
  assert.equal(runnerLib.hostAllowed(undefined, allowed), false);
});

test("settingsFromService: nulls are the panel's defaults; bookMode picks default books, every live book or the ticks", () => {
  const nothing = edgesPayload.settingsFromService({ bankroll: null, multiplier: null, bookMode: null, bookIds: null, minEdgePct: null });
  assert.deepEqual(nothing.stakeSettings, edgeRows.DEFAULT_STAKE_SETTINGS);
  assert.deepEqual(nothing.edgeSettings, edgeRows.DEFAULT_EDGE_SETTINGS);
  assert.deepEqual(edgesPayload.settingsFromService(null), nothing);
  const set = edgesPayload.settingsFromService({ bankroll: 8000, multiplier: 0.5, bookMode: "all", minEdgePct: 2, includeAlts: true });
  assert.deepEqual(set.stakeSettings, { bankroll: 8000, multiplier: 0.5 });
  assert.deepEqual([set.edgeSettings.bookIds, set.edgeSettings.minEdgePct, set.edgeSettings.includeAlts], [null, 2, true]);
  assert.deepEqual(edgesPayload.settingsFromService({ bookMode: "custom", bookIds: [89] }).edgeSettings.bookIds, [89]);
  assert.equal(edgesPayload.settingsFromService({ bookMode: "default", bookIds: null }).edgeSettings.bookIds, undefined);
});

// The panel's own path (edgerows.js, tailflex.js, teaser.js) over the runner's state, for parity checks.
function panelCards(state, serviceSettings, bets) {
  const { stakeSettings, edgeSettings } = edgesPayload.settingsFromService(serviceSettings);
  const records = betsLib.resolveTeamKeys(bets, []);
  const boardLines = edgeRows.boardLines(state);
  const rows = edgeRows.listedEdgeRows(state, {
    edgeSettings, stakeSettings, effective: edgeRows.effectiveFilter(edgeSettings, state, null, null), records, boardLines,
    ladderReaderOf: edgeRows.createLadderReaders(state),
    teasers: teaser.openTeasers(records, boardLines, { now: NOW, ladderOf: teaser.teaserBoardOf(state).ladderOf }),
    measurement: tailflex.measureTailFlex(state, { now: NOW, maxLineAgeMs: edgeSettings.maxLineAgeHours * 3600 * 1000 }), now: NOW,
  });
  return { rows, cards: edgeRows.groupEdgeRows(rows, edgeSettings.sortBy) };
}

test("edges.json lists exactly what the panel's modules list for the same settings, sized against the held bet", async () => {
  const { runner, fetchImpl } = await startedRunner({ "/bets.json": betsBody([BEARS_HELD]), "/settings.json": { settings: OPEN_SETTINGS, updatedAt: "2026-10-01T12:00:00Z" } });
  assert.deepEqual(fetchImpl.calls, NFL_AND_CFB);
  const payload = runner.edgesPayload();
  assert.equal(payload.generatedAt, new Date(NOW).toISOString());
  assert.deepEqual([payload.settings.source, payload.settings.error, payload.settings.updatedAt], ["service", null, "2026-10-01T12:00:00Z"]);
  assert.equal(payload.settings.edges.bookIds, null);
  assert.deepEqual([payload.scanner.phase, payload.scanner.leagues, payload.scanner.leaguesLoaded], ["live", [1, 2], [1]]);
  assert.deepEqual([payload.betsService.openBets, payload.betsService.error], [1, null]);
  assert.deepEqual([payload.books.mode, payload.books.ids, payload.books.liveCount], ["all", null, 5]);
  assert.deepEqual([payload.grouped, payload.unit, payload.total, payload.maxShown], [true, "cards", 4, 200]);
  assert.equal(payload.tailFlex, "tail flex: NFL spr 10% · NFL tot 10%");

  // Parity: the panel's own path over the same state and settings.
  const { rows, cards } = panelCards(runner.scanner.getState(), OPEN_SETTINGS, [BEARS_HELD]);
  assert.deepEqual(payload.items.map((card) => [card.key, card.best.key, card.others.map((row) => row.key)]),
    cards.map((card) => [card.key, card.best.key, card.rows.slice(1).map((row) => row.key)]));
  assert.deepEqual(payload.items.map((card) => card.best.rail), cards.map((card) => edgeRows.stakeRail(card.best)));
  assert.deepEqual(payload.items.flatMap((card) => [card.best, ...card.others].map((row) => row.rankScore)),
    cards.flatMap((card) => card.rows.map((row) => row.rankScore)));
  assert.equal(rows.length, 14);

  const held = payload.items.find((card) => card.best.key === BEARS_SPREAD_SOUTHPOINT).best;
  assert.deepEqual(held.rail, { text: "add $237.25", note: "$437.25 alone", atSize: false });
  assert.deepEqual(held.badges, [{ kind: "held", text: "held $200" }]);
  assert.deepEqual([held.advice.kind, held.advice.verb, held.advice.held, held.advice.bet], ["sized", "add", 200, 237.25]);
  assert.deepEqual(held.advice.teasers, { held: 0, against: 0 });
  assert.equal(held.related.length, 1);
  assert.equal(held.related[0].tier, "same_line");
  assert.equal(held.move, null); // one snapshot seen: nothing has moved yet
  assert.equal(held.edgeTier, "hot");
  assert.ok(held.rankScore > 0);
  assert.ok(!("venueIds" in held) && !("bet" in held));
});

test("an open BFA teaser sizes the rows as on the panel: the Panthers lines at size, the Bears side bigger", async () => {
  const bets = [BEARS_HELD, ...PANTHERS_TEASER];
  const { runner } = await startedRunner({ "/bets.json": betsBody(bets), "/settings.json": { settings: { ...OPEN_SETTINGS, groupByMarket: false }, updatedAt: null } });
  const payload = runner.edgesPayload();
  const { rows } = panelCards(runner.scanner.getState(), { ...OPEN_SETTINGS, groupByMarket: false }, bets);
  assert.deepEqual(payload.items.map((row) => [row.key, row.advice.bet, row.advice.teasers, row.rail]),
    rows.map((row) => [row.key, row.bet.advice.bet, row.bet.advice.teasers, edgeRows.stakeRail(row)]));
  const panthers = payload.items.find((row) => row.key === "289357357:ms89:si1:tid5:alt-2.5");
  assert.deepEqual([panthers.advice.teasers, panthers.rail.atSize], [{ held: 400, against: 0 }, true]);
  assert.ok(panthers.badges.some((badge) => badge.text === "teasers $400"));
  const bears = payload.items.find((row) => row.key === BEARS_SPREAD_SOUTHPOINT);
  assert.deepEqual(bears.advice.teasers, { held: 0, against: 400 });
  assert.ok(bears.advice.bet > 237.25);
});

test("a flat list when cards are off, and the minimum edge and sort come from the settings row", async () => {
  const { runner } = await startedRunner({ "/bets.json": betsBody([]), "/settings.json": { settings: { ...OPEN_SETTINGS, groupByMarket: false, minEdgePct: 5, sortBy: "stake" }, updatedAt: null } });
  const payload = runner.edgesPayload();
  assert.deepEqual([payload.grouped, payload.unit, payload.total], [false, "lines", 3]);
  assert.deepEqual(payload.items.map((row) => row.key), [BEARS_SPREAD_SOUTHPOINT, "289357357:ms89:si1:tid5:alt-2.5", "289357357:ms89:si1:tid5:alt-6.5"]);
  assert.deepEqual(payload.items[0].rail, { text: "bet $437.25", note: null, atSize: false });
});

test("the bets service down: the panel's defaults, both errors named, nothing sized against bets", async () => {
  const down = new Error("connect ECONNREFUSED 127.0.0.1:8094");
  const { runner, fetchImpl } = await startedRunner({ "/bets.json": down, "/settings.json": down });
  const payload = runner.edgesPayload();
  const named = `${SERVICE_URL}/settings.json: ${down.message}`;
  assert.deepEqual([payload.settings.source, payload.settings.error], ["defaults", named]);
  assert.deepEqual(payload.settings.stake, edgeRows.DEFAULT_STAKE_SETTINGS);
  assert.equal(payload.settings.edges.bookIds, "default");
  assert.deepEqual([payload.betsService.error, payload.betsService.unreachableSince, payload.betsService.openBets], [`${SERVICE_URL}/bets.json: ${down.message}`, NOW, 0]);
  // The default leagues are every league; only NFL answers here.
  assert.equal(fetchImpl.calls.length, edgeRows.ALL_LEAGUE_IDS.length);
  assert.equal(payload.books.mode, "default");
  assert.equal(payload.total, 0); // the default books post no main line in the slice, and its lines are older than a week
  assert.equal(runner.health().settings.error, named);
  runner.stop();
});

test("a settings change of leagues restarts the scan; a failed settings poll keeps the last good settings", async () => {
  let settingsBody = { settings: { ...OPEN_SETTINGS, leagues: [1] }, updatedAt: null };
  const fetchImpl = unabatedFetch();
  const service = serviceFetch({ "/bets.json": betsBody([]), "/settings.json": () => settingsBody });
  const runner = runnerLib.createRunner({ fetchImpl, serviceFetch: service, betsServiceUrl: SERVICE_URL, now: () => NOW, timers: noTimers, ...quiet });
  await runner.start();
  assert.deepEqual(fetchImpl.calls, NFL_AND_CFB);
  await runner.pollSettings(); // unchanged leagues: no reload
  assert.deepEqual(fetchImpl.calls, NFL_AND_CFB);
  // Ticking CFB changes the list, not the scan: the scanner holds NFL and CFB anyway.
  settingsBody = { settings: { ...OPEN_SETTINGS, leagues: [1, 2] }, updatedAt: null };
  await runner.pollSettings();
  assert.deepEqual(fetchImpl.calls, NFL_AND_CFB);
  settingsBody = { settings: { ...OPEN_SETTINGS, leagues: [1, 3] }, updatedAt: null };
  await runner.pollSettings();
  await runner.scanLoaded();
  assert.deepEqual(runner.edgesPayload().scanner.leagues, [1, 2, 3]);
  assert.deepEqual(fetchImpl.calls.slice(2).sort(), [SNAPSHOT_BASE_URL(1), SNAPSHOT_BASE_URL(2), SNAPSHOT_BASE_URL(3)].sort());
  settingsBody = { nope: true };
  await runner.pollSettings();
  const payload = runner.edgesPayload();
  assert.deepEqual([payload.settings.source, payload.settings.error], ["service", "settings.json has no settings object"]);
  assert.deepEqual(payload.settings.edges.leagues, [1, 3]);
});

test("scannerLeaguesOf: the settings' leagues plus NFL and CFB, sorted, once each", () => {
  assert.deepEqual(runnerLib.scannerLeaguesOf([]), [1, 2]);
  assert.deepEqual(runnerLib.scannerLeaguesOf([5, 1]), [1, 2, 5]);
  // The open bets' leagues load too, for their closing fairs.
  assert.deepEqual(runnerLib.scannerLeaguesOf([1], [5, 12]), [1, 2, 5, 12]);
});

test("Football unticked: the scan still holds NFL for the teasers, the list shows none of it", async () => {
  const { runner, fetchImpl } = await startedRunner({ "/bets.json": betsBody([]), "/settings.json": { settings: { ...OPEN_SETTINGS, leagues: [3] }, updatedAt: null } });
  assert.deepEqual(fetchImpl.calls.sort(), [...NFL_AND_CFB, SNAPSHOT_BASE_URL(3)].sort());
  const payload = runner.edgesPayload();
  assert.deepEqual([payload.scanner.leaguesLoaded, payload.total, payload.tailFlex], [[1], 0, ""]);
});

test("HTTP: /edges.json and /health on loopback; a foreign Host is 403, another verb 405, another path 404", async () => {
  const { runner } = await startedRunner({ "/bets.json": betsBody([BEARS_HELD]), "/settings.json": { settings: OPEN_SETTINGS, updatedAt: null } });
  const server = runnerLib.createHttpServer(runner, { host: "127.0.0.1" });
  await new Promise((resolve) => server.listen(0, "127.0.0.1", resolve));
  const port = server.address().port;
  const get = (route, options) => new Promise((resolve, reject) => {
    const request = require("node:http").request({ host: "127.0.0.1", port, path: route, method: "GET", ...options }, (reply) => {
      let body = "";
      reply.on("data", (chunk) => { body += chunk; });
      reply.on("end", () => resolve({ status: reply.statusCode, headers: reply.headers, body: JSON.parse(body) }));
    });
    request.on("error", reject);
    request.end();
  });
  try {
    const edges = await get("/edges.json");
    assert.equal(edges.status, 200);
    assert.equal(edges.headers["cache-control"], "no-store");
    assert.deepEqual(edges.body, JSON.parse(JSON.stringify(runner.edgesPayload())));
    const health = await get("/health");
    assert.deepEqual([health.status, health.body.ok, health.body.scanner.phase], [200, true, "live"]);
    const foreign = await get("/edges.json", { headers: { Host: "evil.example" } });
    assert.equal(foreign.status, 403);
    assert.match(foreign.body.error, /Host must be one of/);
    assert.equal((await get("/edges.json", { headers: { Host: `localhost:${port}` } })).status, 200);
    assert.equal((await get("/edges.json", { method: "POST" })).status, 405);
    assert.equal((await get("/bets.json")).status, 404);
  } finally {
    await new Promise((resolve) => server.close(resolve));
    runner.stop();
  }
});

test("closing fairs: the open bet's fair on its own line, POSTed once, again only when it changes, never after kickoff", async () => {
  let clock = NOW;
  const posts = [];
  const service = async (url, init = {}) => {
    const route = url.slice(SERVICE_URL.length);
    if (route === "/closing_fairs.json") {
      posts.push(JSON.parse(init.body).rows);
      assert.deepEqual([init.method, init.headers["Content-Type"]], ["POST", "application/json"]);
      return response({ body: JSON.stringify({ ok: true, saved: 1 }) });
    }
    if (route === "/bets.json") return response({ body: JSON.stringify(betsBody([BEARS_HELD])) });
    if (route === "/settings.json") return response({ body: JSON.stringify({ settings: OPEN_SETTINGS, updatedAt: null }) });
    return response({ status: 404, body: "{}" });
  };
  const runner = runnerLib.createRunner({ fetchImpl: unabatedFetch(), serviceFetch: service, betsServiceUrl: SERVICE_URL, now: () => clock, timers: noTimers, ...quiet });
  await runner.start();
  await runner.postClosingFairs();
  assert.equal(posts.length, 1);
  // BetOnline is not in the slice, so the most recently changed book at -2.5 carries the fair.
  assert.deepEqual(posts[0], [{
    betId: "bol-1", lineKey: "289357360:ms105:si0:tid6", points: -2.5, fairAmerican: -123,
    fairObservedAt: new Date(NOW).toISOString(), eventStart: new Date(KICKOFF_MS).toISOString(),
  }]);
  await runner.postClosingFairs();
  assert.equal(posts.length, 1, "nothing changed, nothing resent");
  clock = KICKOFF_MS + 1000;
  await runner.postClosingFairs();
  assert.equal(posts.length, 1, "after kickoff the stored row is the close");
  assert.equal(runner.health().closingFairs.rowsSent, 1);
  runner.stop();
});

test("closing fairs: a league whose snapshot build is stale gives no close", async () => {
  const posts = [];
  const service = async (url, init = {}) => {
    const route = url.slice(SERVICE_URL.length);
    if (route === "/closing_fairs.json") { posts.push(JSON.parse(init.body).rows); return response({ body: "{}" }); }
    if (route === "/bets.json") return response({ body: JSON.stringify(betsBody([BEARS_HELD])) });
    if (route === "/settings.json") return response({ body: JSON.stringify({ settings: OPEN_SETTINGS, updatedAt: null }) });
    return response({ status: 404, body: "{}" });
  };
  // The CDN serves an NFL snapshot built an hour ago: the scanner marks the league stale.
  const staleFetch = async (url) => {
    if (url.split("?")[0] !== SNAPSHOT_BASE_URL(1)) return response({ status: 404, body: "not found" });
    return response({ body: sliceBody, headers: { "content-length": "1000", "last-modified": new Date(NOW - 3600 * 1000).toUTCString() } });
  };
  const runner = runnerLib.createRunner({ fetchImpl: staleFetch, serviceFetch: service, betsServiceUrl: SERVICE_URL, now: () => NOW, timers: noTimers, ...quiet });
  await runner.start();
  assert.deepEqual(runner.edgesPayload().scanner.staleLeagues, [1]);
  await runner.postClosingFairs();
  assert.deepEqual(posts, []);
  runner.stop();
});

test("an open bet's league is loaded for its close, and stays loaded after the bet settles", async () => {
  const mlbBet = { ...BEARS_HELD, id: "k-mlb", venue: "kalshi", league: "mlb", betType: "total", side: "over", points: 8.5 };
  let bets = [mlbBet];
  const fetchImpl = unabatedFetch();
  const service = serviceFetch({ "/bets.json": () => betsBody(bets), "/settings.json": { settings: { ...OPEN_SETTINGS, leagues: [1] }, updatedAt: null } });
  const runner = runnerLib.createRunner({ fetchImpl, serviceFetch: service, betsServiceUrl: SERVICE_URL, now: () => NOW, timers: noTimers, ...quiet });
  await runner.start();
  assert.deepEqual(runner.edgesPayload().scanner.leagues, [1, 2, 5, 12]);
  const loads = fetchImpl.calls.length;
  bets = [{ ...mlbBet, status: "won" }];
  await runner.pollBets();
  await runner.scanLoaded();
  assert.equal(fetchImpl.calls.length, loads, "no restart when the last MLB bet settles");
  assert.deepEqual(runner.edgesPayload().scanner.leagues, [1, 2, 5, 12]);
  runner.stop();
});
