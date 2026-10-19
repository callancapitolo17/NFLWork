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
  assert.deepEqual(runnerLib.configFromEnv({}), { host: "127.0.0.1", port: 8095, betsServiceUrl: "http://127.0.0.1:8094", kalshiBaseUrl: "https://api.elections.kalshi.com/trade-api/v2", novigBaseUrl: "https://api.novig.com" });
  assert.deepEqual(runnerLib.configFromEnv({ UNABATED_RUNNER_HOST: "100.64.0.7", UNABATED_RUNNER_PORT: "9000", BETS_SERVICE_URL: "http://127.0.0.1:8094/" }),
    { host: "100.64.0.7", port: 9000, betsServiceUrl: "http://127.0.0.1:8094", kalshiBaseUrl: "https://api.elections.kalshi.com/trade-api/v2", novigBaseUrl: "https://api.novig.com" });
  assert.throws(() => runnerLib.configFromEnv({ UNABATED_RUNNER_PORT: "80a" }), /UNABATED_RUNNER_PORT must be a whole number from 1 to 65535, got "80a"/);
  assert.throws(() => runnerLib.configFromEnv({ BETS_SERVICE_URL: "ftp://x" }), /BETS_SERVICE_URL must be an http\(s\) URL/);
  assert.throws(() => runnerLib.configFromEnv({ KALSHI_BASE_URL: "ftp://x" }), /KALSHI_BASE_URL must be an http\(s\) URL/);
  assert.throws(() => runnerLib.configFromEnv({ NOVIG_API_URL: "ftp://x" }), /NOVIG_API_URL must be an http\(s\) URL/);
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
    assert.deepEqual((await get("/teasers.json")).body, JSON.parse(JSON.stringify(runner.teasersPayload())));
    const missing = await get("/attach.json?betId=nope");
    assert.deepEqual([missing.status, missing.body.error], [404, 'no open bet "nope"']);
    assert.equal((await get("/attach.json?betId=bol-1&query=bears")).body.events[0].eventId, 125807);
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

// ---- board rows kept per league ------------------------------------------------------------

// The NFL slice re-labelled as CFB, with its own market and event ids so nothing collides with NFL's;
// `dropPregame` leaves out its pregame game, so CFB's board loses its one event.
function cfbSliceBody({ dropPregame = false } = {}) {
  const json = JSON.parse(sliceBody);
  if (dropPregame) {
    delete json.odds["lg1:pt1:pregame"];
    delete json.odds["lg1:pt2:pregame"];
  }
  return JSON.stringify(json).replace(/lg1:/g, "lg2:").replace(/"marketId":(\d+)/g, '"marketId":9$1').replace(/"eventId":(\d+)/g, '"eventId":9$1');
}

// NFL and CFB from the slice; `feed.failNfl` makes NFL's next downloads fail
// (it keeps its lines), `feed.dropCfbPregame` changes CFB's next snapshot.
function twoLeagueFetch(feed) {
  return async (url) => {
    const base = url.split("?")[0];
    const headers = { "content-length": "1000", "last-modified": new Date(NOW - 20000).toUTCString() };
    if (base === SNAPSHOT_BASE_URL(1)) return feed.failNfl ? response({ status: 503, body: "" }) : response({ body: sliceBody, headers });
    if (base === SNAPSHOT_BASE_URL(2)) return response({ body: cfbSliceBody({ dropPregame: feed.dropCfbPregame }), headers });
    return response({ status: 404, body: "not found" });
  };
}

test("boardLines: built league by league, the same rows a walk of the whole board gives; a league that did not land keeps its rows untouched", async () => {
  const feed = { failNfl: false, dropCfbPregame: false };
  const runner = runnerLib.createRunner({
    fetchImpl: twoLeagueFetch(feed), betsServiceUrl: SERVICE_URL, now: () => NOW, timers: noTimers, ...quiet,
    serviceFetch: serviceFetch({ "/bets.json": betsBody([BEARS_HELD]), "/settings.json": { settings: OPEN_SETTINGS } }),
  });
  await runner.start();
  const walked = () => edgeRows.boardLines(runner.scanner.getState());
  const rowsOf = (leagueId) => runner.boardLines().filter((row) => row.leagueId === leagueId);
  assert.deepEqual(runner.boardLines(), walked());
  const nflRows = rowsOf(1);
  assert.equal(rowsOf(2).length, 1);

  // Only CFB lands, without its game.
  Object.assign(feed, { failNfl: true, dropCfbPregame: true });
  await runner.scanner.resync();
  assert.deepEqual(runner.boardLines(), walked(), "the landed league's rows are rebuilt");
  assert.equal(rowsOf(2).length, 0);
  assert.ok(rowsOf(1).every((row, index) => row === nflRows[index]) && rowsOf(1).length === nflRows.length, "NFL did not land: its rows are the same objects");
  runner.stop();
});

test("kickoff odds: a game's lines are re-read when its league lands, kept when the league's download fails", async () => {
  const feed = { failNfl: false, dropCfbPregame: false };
  const clock = { ms: NOW };
  const runner = runnerLib.createRunner({
    fetchImpl: twoLeagueFetch(feed), betsServiceUrl: SERVICE_URL, now: () => clock.ms, timers: noTimers, ...quiet,
    serviceFetch: serviceFetch({ "/bets.json": betsBody([BEARS_HELD]), "/settings.json": { settings: OPEN_SETTINGS } }),
  });
  await runner.start();
  clock.ms = NOW + 60 * 1000;
  await runner.scanner.resync();
  clock.ms = NOW + 120 * 1000;
  feed.failNfl = true;
  await runner.scanner.resync();
  clock.ms = KICKOFF_MS + 3600 * 1000;
  const card = runner.scenariosPayload().games.find((game) => game.eventId === 125807);
  assert.equal(card.oddsAt, NOW + 60 * 1000, "the last NFL landing before kickoff, not the failed one after it");
  runner.stop();
});

// ---- /scenarios.json (the Bet Tracker's Live tab) ----------------------------------------

async function runnerOnClock(clock, bets, extra) {
  const runner = runnerLib.createRunner({
    fetchImpl: unabatedFetch(), betsServiceUrl: SERVICE_URL, now: () => clock.ms, timers: noTimers, ...quiet,
    serviceFetch: serviceFetch({ "/bets.json": { ...betsBody(bets), ...extra }, "/settings.json": { settings: OPEN_SETTINGS } }),
  });
  await runner.start();
  return runner;
}

test("scenarios: a game seen before kickoff keeps its kickoff chances once it is under way", async () => {
  const clock = { ms: NOW };
  const runner = await runnerOnClock(clock, [BEARS_HELD]);
  assert.deepEqual(runner.scenariosPayload().games, [], "nothing has started yet");
  clock.ms = KICKOFF_MS + 3600 * 1000;
  const payload = runner.scenariosPayload();
  assert.equal(payload.games.length, 1);
  const [card] = payload.games;
  assert.deepEqual([card.awayTeam, card.homeTeam, card.oddsAt], ["Chicago Bears", "Carolina Panthers", NOW]);
  assert.deepEqual(payload.coveredBetIds, ["bol-1"]);
  const [group] = card.groups;
  assert.deepEqual(group.bands.map((band) => band.pnl), [-200, 181.82]);
  assert.ok(group.bands.every((band) => band.prob > 0 && band.prob < 1), "both results priced off the kickoff ladder");
  runner.stop();
});

test("scenarios: started before the runner saw it, the card has no chances; a removed bet is never on a card", async () => {
  const clock = { ms: KICKOFF_MS + 3600 * 1000 };
  const runner = await runnerOnClock(clock, [BEARS_HELD]);
  const [card] = runner.scenariosPayload().games;
  assert.equal(card.oddsAt, null);
  assert.ok(card.groups[0].bands.every((band) => band.prob === null));
  runner.stop();
  const removed = await runnerOnClock(clock, [BEARS_HELD], { exclusions: [{ betId: "bol-1", venue: "betonline", excludedAt: "2026-09-13T16:00:00Z" }] });
  assert.deepEqual(removed.scenariosPayload().games, []);
  removed.stop();
});

// ---- /teasers.json, Attach and the shared Dismiss marks ---------------------------------

test("teasers.json: the panel's Teasers tab in its words, the open BFA teaser held fixed", async () => {
  const { runner } = await startedRunner({ "/bets.json": betsBody([BEARS_HELD, ...PANTHERS_TEASER]), "/settings.json": { settings: OPEN_SETTINGS, updatedAt: null } });
  const payload = runner.teasersPayload();
  runner.stop();
  assert.equal(payload.ready, true);
  assert.equal(payload.error, null);
  assert.equal(payload.status, "NFL · CFB unavailable · 0 games on Buckeye's board · 0 of 0 legs priced");
  assert.match(payload.empty, /^Fewer than 4 games/);
  assert.deepEqual(payload.summary.cells, [{ label: "Placed", value: "$400", small: null }]);
  assert.equal(payload.open.length, 1);
  assert.deepEqual([payload.open[0].size, payload.open[0].stake, payload.open[0].inPlay], ["2-team · $400 to win $1,200", 400, true]);
  assert.deepEqual(payload.open[0].legs[1], { label: "Elsewhere +3.5", note: "no board game; counted as won", win: "—", counted: true });
  assert.equal(payload.breakEven, "70.7% a leg breaks even at +300");
});

test("attach.json: step 1 lists the bet's games best first, step 2 the plan and the pins body; a settled or unknown bet is a 404", async () => {
  const { runner } = await startedRunner({ "/bets.json": betsBody([BEARS_HELD]), "/settings.json": { settings: OPEN_SETTINGS, updatedAt: null } });
  const pick = runner.attachPayload(new URLSearchParams({ betId: "bol-1", query: "" }));
  assert.equal(pick.scope, "NFL · Sat Sep 12 to Mon Sep 14");
  assert.deepEqual(pick.events[0], { eventId: 125807, label: "Chicago Bears @ Carolina Panthers", meta: "Sun Sep 13 1:00 PM ET · rot 465 / 466", why: "Chicago Bears matches" });
  const confirm = runner.attachPayload(new URLSearchParams({ betId: "bol-1", eventId: "125807", swapped: "0" }));
  assert.equal(confirm.label, "Chicago Bears @ Carolina Panthers");
  assert.deepEqual([confirm.pinRequest.pin.betId, confirm.pinRequest.pin.eventId], ["bol-1", "125807"]);
  assert.equal(confirm.venueLabel, betsLib.venueLabel("betonline"));
  assert.throws(() => runner.attachPayload(new URLSearchParams({ betId: "nope" })), (error) => error.status === 404);
  assert.throws(() => runner.attachPayload(new URLSearchParams({ betId: "bol-1", eventId: "1" })), /has left the board/);
  runner.stop();
});

test("a bet Dismissed on either page stops needing a game in edges.json, and lists as dismissed", async () => {
  const stray = { ...BEARS_HELD, id: "bol-stray", awayTeam: "Nowhere Owls", homeTeam: "Lost Llamas" };
  const flagged = await startedRunner({ "/bets.json": betsBody([stray]), "/settings.json": { settings: OPEN_SETTINGS, updatedAt: null } });
  const before = flagged.runner.edgesPayload().betsService.unmatched.find((entry) => entry.betId === "bol-stray");
  flagged.runner.stop();
  assert.deepEqual([before.needsGame, before.dismissed], [true, false]);
  const body = { ...betsBody([stray]), dismissals: [{ betId: "bol-stray", dismissedAt: "2026-09-13T15:00:00Z" }] };
  const { runner } = await startedRunner({ "/bets.json": body, "/settings.json": { settings: OPEN_SETTINGS, updatedAt: null } });
  const after = runner.edgesPayload().betsService.unmatched.find((entry) => entry.betId === "bol-stray");
  runner.stop();
  assert.deepEqual([after.needsGame, after.dismissed], [false, true]);
});

// ---- the automated taker (server/taker.js, T1: decides, places nothing) ----
// The CFB venue-id slice (Duquesne @ Washington State, Kalshi rung ids) served
// as league 2, the clock four days before its kickoff, a 5% Unabated edge set
// on Kalshi's WSU -35.5 alt (-143, 0.587157 all-in) and every other book at
// that number made worse, so the row is the best price.
const TAKER_NOW = Date.parse("2026-09-15T19:00:00Z");
const TAKER_TICKER = "KXNCAAFSPREAD-26SEP19DUQWSU-WSU36";

function takerVenueSlice() {
  const json = JSON.parse(fs.readFileSync(path.join(__dirname, "fixtures", "v2_venue_ids_slice.json"), "utf8"));
  const spreadRow = json.odds["lg2:pt1:pregame"].find((row) => row.key === "pt1:pregame:bt2:e123742");
  for (const [sideKey, books] of Object.entries(spreadRow.sides)) {
    for (const [bookKey, book] of Object.entries(books)) {
      for (const rung of book.alternateLines || []) {
        if (rung && rung.points === -41.5 && sideKey === "si1:tid717" && bookKey === "ms105") Object.assign(rung, { ge: 0.06, bacr: 100, sequenceNumber: TAKER_NOW - 60000 });
        if (!rung || rung.points !== -35.5 || sideKey !== "si1:tid717") continue;
        if (bookKey === "ms105") Object.assign(rung, { ge: 0.05, bacr: -160, sequenceNumber: TAKER_NOW - 60000 });
        else if (rung.americanPrice != null) Object.assign(rung, { americanPrice: -200, price: -200 });
      }
    }
  }
  return JSON.stringify(json);
}

function takerUnabatedFetch() {
  const body = takerVenueSlice();
  return async (url) => {
    if (url.split("?")[0] === SNAPSHOT_BASE_URL(2)) return response({ body, headers: { "content-length": "1000", "last-modified": new Date(TAKER_NOW - 20000).toUTCString() } });
    return response({ status: 404, body: "not found" });
  };
}

// Kalshi: a 57c YES ask with 40 contracts behind it on the WSU -35.5 market,
// 46c on WSU -41.5, and both markets active at their strikes. Its milestone
// for the game starts when Unabated's does, unless `kalshiStart` says another
// time (null: no milestone; an Error: the read fails).
const WSU_EVENT_TICKER = "KXNCAAFSPREAD-26SEP19DUQWSU";
const WSU_UNABATED_START = "2026-09-19T19:30:00.000Z";
function kalshiFetch(markets = {}, kalshiStart = WSU_UNABATED_START) {
  const calls = [];
  const impl = async (url) => {
    calls.push(url);
    if (url.includes("/milestones?")) {
      if (kalshiStart instanceof Error) return response({ status: 500, body: JSON.stringify({ error: kalshiStart.message }) });
      const named = new URL(url).searchParams.get("related_event_ticker");
      const milestones = kalshiStart === null || named !== WSU_EVENT_TICKER ? []
        : [{ id: "m-1", type: "football_game", start_date: kalshiStart, related_event_tickers: [WSU_EVENT_TICKER, "KXNCAAFGAME-26SEP19DUQWSU"] }];
      return response({ body: JSON.stringify({ milestones, cursor: "" }) });
    }
    if (url.endsWith("/markets/KXNCAAFSPREAD-26SEP19DUQWSU-WSU42/orderbook")) {
      return response({ body: JSON.stringify({ orderbook_fp: { yes_dollars: [["0.4400", "100"]], no_dollars: [["0.5400", "40"]] } }) });
    }
    if (url.endsWith(`/markets/${TAKER_TICKER}/orderbook`)) {
      return response({ body: JSON.stringify({ orderbook_fp: { yes_dollars: [["0.5500", "100"]], no_dollars: [["0.4300", "40"]] } }) });
    }
    if (url.endsWith("/markets/KXNCAAFSPREAD-26SEP19DUQWSU-WSU42")) {
      return response({ body: JSON.stringify({ market: { status: "active", floor_strike: 41.5, yes_sub_title: "Washington St. wins by over 41.5 points", ...markets.WSU42 } }) });
    }
    if (url.endsWith(`/markets/${TAKER_TICKER}`)) {
      return response({ body: JSON.stringify({ market: { status: "active", floor_strike: 35.5, yes_sub_title: "Washington St. wins by over 35.5 points", ...markets.WSU36 } }) });
    }
    return response({ status: 404, body: JSON.stringify({ error: "not found" }) });
  };
  impl.calls = calls;
  return impl;
}

// The bets service with the desk's unseen-fill rule (taker_orders.py): an
// order on a game whose earlier filled bot order is not among the request's
// heldIds is refused, else it fills in full.
function takerService(script) {
  const posts = [];
  const filledByEvent = new Map();
  const impl = async (url, init = {}) => {
    const route = url.slice(SERVICE_URL.length);
    if (init.method === "POST") {
      const body = JSON.parse(init.body);
      posts.push({ route, body });
      if (route === "/taker/order.json") {
        if (Date.parse(body.eventStart) <= TAKER_NOW) return response({ body: JSON.stringify({ status: "started", detail: `the game started at ${body.eventStart}` }) });
        const unseen = (filledByEvent.get(body.eventId) || []).find((heldId) => !body.heldIds.includes(heldId));
        if (unseen) return response({ body: JSON.stringify({ status: "unseen_fill", detail: `${unseen} is not yet among the held bets` }) });
        filledByEvent.set(body.eventId, [...(filledByEvent.get(body.eventId) || []), `kalshi:${body.ticker}:${body.side}`]);
        return response({ body: JSON.stringify({ status: "filled", count: body.count, filled: body.count, cost: body.count * body.allIn, orderId: "k-1" }) });
      }
      return response({ body: JSON.stringify({ ok: true, logged: body.rows.length }) });
    }
    const answer = script[route];
    if (answer instanceof Error) throw answer;
    if (answer === undefined) return response({ status: 404, body: "{}" });
    return response({ body: JSON.stringify(answer) });
  };
  impl.posts = posts;
  impl.script = script;
  return impl;
}

const TAKER_SETTINGS = { settings: { leagues: [2], periods: [1] }, updatedAt: null };

// The bets service as the taker needs it: its Kalshi poll succeeded 30 s ago.
function takerBetsBody(bets = [], kalshiSource = {}) {
  const body = betsBody(bets);
  body.sources.kalshi = { ok: true, fetchedAt: new Date(TAKER_NOW - 30 * 1000).toISOString(), error: null, count: bets.length, ...kalshiSource };
  return body;
}

test("taker: a best-price Kalshi edge at the live ask is sent as one ubot- limit order; the game's next add waits for the fill to be held", async () => {
  const kalshi = kalshiFetch();
  const service = takerService({ "/bets.json": takerBetsBody(), "/settings.json": TAKER_SETTINGS });
  const runner = runnerLib.createRunner({
    fetchImpl: takerUnabatedFetch(), serviceFetch: service, kalshiFetch: kalshi, kalshiBaseUrl: "https://kalshi.test/v2",
    betsServiceUrl: SERVICE_URL, now: () => TAKER_NOW, timers: noTimers, ...quiet,
  });
  await runner.start();
  await runner.takerPass();
  const payload = runner.takerPayload();
  assert.equal(payload.pass.skipped, null);
  assert.equal(payload.desk, "armed");
  // WSU -35.5 and WSU -41.5 are both bets: the better live edge (-41.5) goes
  // first, and the desk refuses -35.5 until the -41.5 fill is among its held bets.
  const decision = payload.decisions.find((row) => row.ticker === "KXNCAAFSPREAD-26SEP19DUQWSU-WSU42");
  const waiting = payload.decisions.find((row) => row.ticker === TAKER_TICKER);
  assert.equal(decision.reason, "bet", decision.detail);
  assert.equal(waiting.reason, "bet", waiting.detail);
  assert.equal(waiting.liveAsk, 0.57);
  assert.ok(decision.contracts >= 1 && decision.contracts <= 40, String(decision.contracts));
  assert.equal(payload.decisions[0], decision, "bets list first, best live edge first");
  // Both reached the desk, best live edge first; only the first was sent.
  const orders = service.posts.filter((post) => post.route === "/taker/order.json").map((post) => post.body);
  assert.deepEqual(orders.map((order) => order.ticker), [decision.ticker, TAKER_TICKER]);
  assert.match(orders[0].clientOrderId, /^ubot-[0-9a-f]{32}$/);
  assert.deepEqual({ ticker: orders[0].ticker, side: orders[0].side, ask: orders[0].ask, count: orders[0].count },
    { ticker: decision.ticker, side: "yes", ask: decision.liveAsk, count: decision.contracts });
  assert.deepEqual({ eventId: orders[0].eventId, heldIds: orders[0].heldIds }, { eventId: 123742, heldIds: [] });
  // The desk refuses a decision over a minute old, or priced on another Kalshi.
  assert.deepEqual({ decidedAt: orders[0].decidedAt, kalshiBaseUrl: orders[0].kalshiBaseUrl },
    { decidedAt: new Date(TAKER_NOW).toISOString(), kalshiBaseUrl: "https://kalshi.test/v2" });
  assert.equal(decision.order.status, "filled");
  assert.equal(waiting.order.status, "unseen_fill");
  // Each order was preceded by a read of its market; everything else was orderbooks.
  assert.ok(kalshi.calls.includes("https://kalshi.test/v2/markets/KXNCAAFSPREAD-26SEP19DUQWSU-WSU42"));
  assert.ok(kalshi.calls.every((url) => url.startsWith("https://kalshi.test/v2/markets/") || url.startsWith("https://kalshi.test/v2/milestones?")),
    "only market, orderbook and milestone reads");
  assert.equal(orders[0].eventStart, WSU_UNABATED_START, "Kalshi's start agrees with Unabated's");
  assert.equal(kalshi.calls.filter((url) => url.includes("/milestones?")).length, 1, "one milestone read per game, cached");
  const logged = service.posts.filter((post) => post.route === "/taker_decisions.json").flatMap((post) => post.body.rows);
  assert.equal(logged.filter((row) => row.rowKey === decision.rowKey).length, 1);
  // The next pass: the -41.5 fill is still not among the held bets, so the desk refuses every add on the game.
  await runner.takerPass();
  assert.equal(service.posts.filter((post) => post.route === "/taker/order.json").length, 4);
  assert.ok(runner.takerPayload().decisions.filter((row) => row.reason === "bet").every((row) => row.order.status === "unseen_fill"));
  assert.equal(runner.health().taker.byReason.bet >= 1, true);
  runner.stop();
});

test("taker: Kalshi's start, when earlier than Unabated's, is the one the order carries, so the desk refuses a started game", async () => {
  // Unabated still lists the game four days out; Kalshi says it kicked off a minute ago.
  const kickedOff = new Date(TAKER_NOW - 60 * 1000).toISOString();
  const service = takerService({ "/bets.json": takerBetsBody(), "/settings.json": TAKER_SETTINGS });
  const runner = runnerLib.createRunner({
    fetchImpl: takerUnabatedFetch(), serviceFetch: service, kalshiFetch: kalshiFetch({}, kickedOff), kalshiBaseUrl: "https://kalshi.test/v2",
    betsServiceUrl: SERVICE_URL, now: () => TAKER_NOW, timers: noTimers, ...quiet,
  });
  await runner.start();
  await runner.takerPass();
  runner.stop();
  const orders = service.posts.filter((post) => post.route === "/taker/order.json").map((post) => post.body);
  assert.ok(orders.length >= 1);
  assert.ok(orders.every((order) => order.eventStart === kickedOff), JSON.stringify(orders.map((order) => order.eventStart)));
  const decision = runner.takerPayload().decisions.find((row) => row.ticker === TAKER_TICKER);
  assert.equal(decision.order.status, "started");
  assert.deepEqual({ unabatedStart: decision.unabatedStart, kalshiStart: decision.kalshiStart }, { unabatedStart: WSU_UNABATED_START, kalshiStart: kickedOff });
});

test("taker: no Kalshi start (none listed, or the read fails) trades on Unabated's start alone (option B)", async () => {
  for (const kalshiStart of [null, new Error("boom")]) {
    const service = takerService({ "/bets.json": takerBetsBody(), "/settings.json": TAKER_SETTINGS });
    const runner = runnerLib.createRunner({
      fetchImpl: takerUnabatedFetch(), serviceFetch: service, kalshiFetch: kalshiFetch({}, kalshiStart), kalshiBaseUrl: "https://kalshi.test/v2",
      betsServiceUrl: SERVICE_URL, now: () => TAKER_NOW, timers: noTimers, ...quiet,
    });
    await runner.start();
    await runner.takerPass();
    runner.stop();
    const first = service.posts.find((post) => post.route === "/taker/order.json").body;
    assert.equal(first.eventStart, WSU_UNABATED_START, String(kalshiStart));
    const decision = runner.takerPayload().decisions.find((row) => row.order);
    assert.equal(decision.order.status, "filled");
    assert.equal(decision.kalshiStart, null);
  }
});

test("taker: a market that is closed, or at another strike, is market_mismatch and nothing is sent", async () => {
  const service = takerService({ "/bets.json": takerBetsBody(), "/settings.json": TAKER_SETTINGS });
  const runner = runnerLib.createRunner({
    fetchImpl: takerUnabatedFetch(), serviceFetch: service, kalshiBaseUrl: "https://kalshi.test/v2",
    kalshiFetch: kalshiFetch({ WSU42: { status: "closed" }, WSU36: { floor_strike: 34.5 } }),
    betsServiceUrl: SERVICE_URL, now: () => TAKER_NOW, timers: noTimers, ...quiet,
  });
  await runner.start();
  await runner.takerPass();
  const byTicker = new Map(runner.takerPayload().decisions.map((row) => [row.ticker, row]));
  assert.equal(byTicker.get("KXNCAAFSPREAD-26SEP19DUQWSU-WSU42").reason, "market_mismatch");
  assert.match(byTicker.get("KXNCAAFSPREAD-26SEP19DUQWSU-WSU42").detail, /status "closed"/);
  assert.match(byTicker.get(TAKER_TICKER).detail, /floor_strike 34.5, row is -35.5/);
  assert.equal(service.posts.filter((post) => post.route === "/taker/order.json").length, 0);
  runner.stop();
});

// ---- the Ticket's Place button: GET /kalshi_quote.json (places nothing) ----
async function quoteRunner(markets = {}, extraDeps = {}, kalshiStart = WSU_UNABATED_START) {
  const service = takerService({ "/bets.json": takerBetsBody(), "/settings.json": TAKER_SETTINGS });
  const kalshi = kalshiFetch(markets, kalshiStart);
  const runner = runnerLib.createRunner({
    fetchImpl: takerUnabatedFetch(), serviceFetch: service, kalshiFetch: kalshi, kalshiBaseUrl: "https://kalshi.test/v2",
    betsServiceUrl: SERVICE_URL, now: () => TAKER_NOW, timers: noTimers, ...quiet, ...extraDeps,
  });
  await runner.start();
  return { runner, kalshi, service };
}

const WSU_36_LINE = { eventId: "123742", periodTypeId: "1", betTypeId: "2", sideIndex: "1", points: "-35.5" };

test("kalshi_quote: the board's Kalshi line gives its contract, live ask, all-in and depth; nothing is sent", async () => {
  const { runner, kalshi, service } = await quoteRunner();
  const quote = await runner.kalshiQuote(new URLSearchParams(WSU_36_LINE));
  runner.stop();
  assert.deepEqual(
    { ok: quote.ok, ticker: quote.ticker, side: quote.side, ask: quote.ask, depth: quote.depth, eventId: quote.eventId },
    { ok: true, ticker: TAKER_TICKER, side: "yes", ask: 0.57, depth: 40, eventId: 123742 });
  assert.equal(quote.allIn, Math.round((0.57 + 0.07 * 0.57 * 0.43) * 1e6) / 1e6);
  assert.ok(quote.eventStart);
  assert.equal(quote.kalshiBaseUrl, "https://kalshi.test/v2", "the desk checks it trades the Kalshi this read");
  assert.ok(kalshi.calls.includes(`https://kalshi.test/v2/markets/${TAKER_TICKER}`), "the market is checked");
  assert.equal(service.posts.filter((post) => post.route === "/taker/order.json").length, 0);
});

test("kalshi_quote: while a scan restart is still downloading, the quote reads the board the list still shows", async () => {
  const unabated = takerUnabatedFetch();
  let hold = false;
  const held = [];
  const fetchImpl = (url) => (hold ? new Promise((resolve) => held.push(() => resolve(unabated(url)))) : unabated(url));
  const { runner } = await quoteRunner({}, { fetchImpl });
  hold = true;
  const restart = runner.scanner.start([2]);
  assert.equal(Object.keys(runner.scanner.getState().lines).length, 0, "the scanner's new board is still empty");
  const quote = await runner.kalshiQuote(new URLSearchParams(WSU_36_LINE));
  for (const release of held) release();
  await restart;
  runner.stop();
  assert.deepEqual({ ok: quote.ok, ticker: quote.ticker }, { ok: true, ticker: TAKER_TICKER });
});

test("kalshi_quote: the quote's eventStart is Kalshi's when earlier, so Place refuses a game Unabated still lists as pregame", async () => {
  const kickedOff = new Date(TAKER_NOW - 60 * 1000).toISOString();
  const { runner } = await quoteRunner({}, {}, kickedOff);
  const quote = await runner.kalshiQuote(new URLSearchParams(WSU_36_LINE));
  runner.stop();
  assert.deepEqual({ ok: quote.ok, eventStart: quote.eventStart, unabatedStart: quote.unabatedStart, kalshiStart: quote.kalshiStart },
    { ok: true, eventStart: kickedOff, unabatedStart: WSU_UNABATED_START, kalshiStart: kickedOff });
  const { runner: unlisted } = await quoteRunner({}, {}, null);
  const fallback = await unlisted.kalshiQuote(new URLSearchParams(WSU_36_LINE));
  unlisted.stop();
  assert.deepEqual({ eventStart: fallback.eventStart, kalshiStart: fallback.kalshiStart }, { eventStart: WSU_UNABATED_START, kalshiStart: null });
});

test("kalshi_quote: a closed market, a number the board lacks and bad params are refused", async () => {
  const { runner } = await quoteRunner({ WSU36: { status: "closed" } });
  const closed = await runner.kalshiQuote(new URLSearchParams(WSU_36_LINE));
  const missing = await runner.kalshiQuote(new URLSearchParams({ ...WSU_36_LINE, points: "-99.5" }));
  await assert.rejects(runner.kalshiQuote(new URLSearchParams({ ...WSU_36_LINE, sideIndex: "2" })), (error) => error.status === 400);
  runner.stop();
  assert.deepEqual([closed.ok, closed.reason], [false, "market_mismatch"]);
  assert.deepEqual([missing.ok, missing.reason], [false, "no_line"]);
});

// Novig Place reads Novig's public catalog the taker's way (novigFetch below):
// the WSU -35.5 outcome on a two-outcome SPREAD market, DUQ resting at 0.445.
function novigQuoteRunner(marketFields = {}) {
  const novig = novigFetch(marketFields);
  return quoteRunner({}, { novigFetch: novig, novigBaseUrl: "https://novig.test" }).then((made) => ({ ...made, novig }));
}

test("novig_quote: the board's Novig line gives its outcome, market, the other side's, the live ask and the line; nothing is sent", async () => {
  const { runner, service, novig } = await novigQuoteRunner();
  const quote = await runner.novigQuote(new URLSearchParams(WSU_36_LINE));
  runner.stop();
  assert.deepEqual(
    { ok: quote.ok, outcomeId: quote.outcomeId, otherOutcomeId: quote.otherOutcomeId, marketId: quote.marketId, price: quote.price, depth: quote.depth, novigBaseUrl: quote.novigBaseUrl },
    { ok: true, outcomeId: NOVIG_WSU_OUTCOME, otherOutcomeId: NOVIG_DUQ_OUTCOME, marketId: NOVIG_MARKET, price: 0.555, depth: 500000, novigBaseUrl: "https://novig.test" });
  assert.ok(quote.eventStart);
  assert.deepEqual(quote.line, {
    league: "cfb", betType: "spread", period: "FG", side: "home", points: -35.5,
    awayTeam: quote.line.awayTeam, homeTeam: quote.line.homeTeam, awayTeamId: quote.line.awayTeamId, homeTeamId: "717",
  });
  assert.ok(novig.calls.some((call) => call.url === `https://novig.test/v3/public/catalog/events/${NOVIG_EVENT}`), "the event is checked");
  assert.equal(service.posts.filter((post) => post.route.startsWith("/taker/")).length, 0);
});

test("novig_quote: a main total its ladder has no rung for is found in the catalog off its siblings and checked like any other", async () => {
  // The fixture's Novig 52.5 total hangs 50.5/51.5/53.5/54.5 rungs and no 52.5 one.
  const sibling = {
    marketId: "0c0c0c0c-1111-4222-8333-515151515151", eventId: NOVIG_EVENT, marketType: "TOTAL", strike: "51.5", status: "OPEN",
    fee: { coefficient: "0.06", makerCredit: "0", charged: "WHEN_LIVE" },
    outcomes: [{ outcomeId: "01a0a059-6836-7b83-8813-18c4e5a44a2f", name: "Over 51.5" }, { outcomeId: "01a0a059-6836-7b83-8813-18d718a0283c", name: "Under 51.5" }],
  };
  const mainId = "0c0c0c0c-1111-4222-8333-525252525252";
  const main = { ...sibling, marketId: mainId, strike: "52.5", outcomes: [{ outcomeId: "0d0d-over", name: "Over 52.5" }, { outcomeId: "0d0d-under", name: "Under 52.5" }] };
  const base = novigFetch();
  const impl = async (url, init) => {
    const bare = url.split("?")[0];
    if (bare === "https://novig.test/v3/public/catalog/markets") return response({ body: JSON.stringify({ items: [sibling, main] }) });
    if (bare === `https://novig.test/v3/public/catalog/markets/${mainId}/book`) {
      return response({ body: JSON.stringify({ marketId: mainId, seq: 3, orders: { "0d0d-over": [{ orderId: "o1", price: "0.48", qty: 2000 }], "0d0d-under": [] } }) });
    }
    if (bare === `https://novig.test/v3/public/catalog/markets/${mainId}`) return response({ body: JSON.stringify(main) });
    return base(url, init);
  };
  const { runner } = await quoteRunner({}, { novigFetch: impl, novigBaseUrl: "https://novig.test" });
  const quote = await runner.novigQuote(new URLSearchParams({ eventId: "123742", periodTypeId: "1", betTypeId: "3", sideIndex: "1", points: "52.5" }));
  runner.stop();
  assert.deepEqual(
    { ok: quote.ok, outcomeId: quote.outcomeId, otherOutcomeId: quote.otherOutcomeId, marketId: quote.marketId, price: quote.price, depth: quote.depth },
    { ok: true, outcomeId: "0d0d-under", otherOutcomeId: "0d0d-over", marketId: mainId, price: 0.52, depth: 2000 });
});

test("novig_quote: a market charging a pregame fee, a closed market and a number the board lacks are refused", async () => {
  const fee = await novigQuoteRunner({ fee: { coefficient: "0.06", makerCredit: "0", charged: "ALWAYS" } });
  const feeQuote = await fee.runner.novigQuote(new URLSearchParams(WSU_36_LINE));
  const missing = await fee.runner.novigQuote(new URLSearchParams({ ...WSU_36_LINE, points: "-99.5" }));
  fee.runner.stop();
  assert.deepEqual([feeQuote.ok, feeQuote.reason], [false, "market_mismatch"]);
  assert.deepEqual([missing.ok, missing.reason], [false, "no_line"]);
  const closed = await novigQuoteRunner({ status: "CLOSED" });
  const closedQuote = await closed.runner.novigQuote(new URLSearchParams(WSU_36_LINE));
  closed.runner.stop();
  assert.deepEqual([closedQuote.ok, closedQuote.reason], [false, "market_mismatch"]);
  assert.match(closedQuote.detail, /market status "CLOSED"/);
  // The catalog down, nothing resting on the other side, a market with one outcome: each refused, nothing thrown.
  const down = novigFetch();
  const downImpl = async (url, init) => (url.includes("/catalog/markets?") ? response({ status: 503, body: "down" }) : down(url, init));
  const noCatalog = await quoteRunner({}, { novigFetch: downImpl, novigBaseUrl: "https://novig.test" });
  const noCatalogQuote = await noCatalog.runner.novigQuote(new URLSearchParams(WSU_36_LINE));
  noCatalog.runner.stop();
  assert.deepEqual([noCatalogQuote.ok, noCatalogQuote.reason], [false, "no_outcome"]);
  assert.match(noCatalogQuote.detail, /catalog unread/);
  const empty = novigFetch();
  const emptyImpl = async (url, init) => (url.endsWith("/book")
    ? response({ body: JSON.stringify({ marketId: NOVIG_MARKET, seq: 1, orders: { [NOVIG_DUQ_OUTCOME]: [], [NOVIG_WSU_OUTCOME]: [] } }) }) : empty(url, init));
  const noBook = await quoteRunner({}, { novigFetch: emptyImpl, novigBaseUrl: "https://novig.test" });
  const noBookQuote = await noBook.runner.novigQuote(new URLSearchParams(WSU_36_LINE));
  noBook.runner.stop();
  assert.deepEqual([noBookQuote.ok, noBookQuote.reason, noBookQuote.outcomeId], [false, "no_book", NOVIG_WSU_OUTCOME]);
  assert.match(noBookQuote.detail, /nothing rests/);
});

test("novig_quote: a catalog under a minute old is reused, the market, event and book still read live; a missing outcome reads it fresh", async () => {
  let clock = TAKER_NOW;
  const logs = [];
  const novig = novigFetch();
  const { runner } = await quoteRunner({}, { novigFetch: novig, novigBaseUrl: "https://novig.test", now: () => clock, logInfo: (line) => logs.push(line) });
  const catalogReads = () => novig.calls.filter((call) => call.url.includes("/catalog/markets?")).length;
  const first = await runner.novigQuote(new URLSearchParams(WSU_36_LINE));
  clock += 59 * 1000;
  const reused = await runner.novigQuote(new URLSearchParams(WSU_36_LINE));
  assert.deepEqual([first.ok, reused.ok, catalogReads()], [true, true, 1], "the second quote reuses the first catalog");
  assert.equal(novig.calls.filter((call) => call.url.endsWith(`/markets/${NOVIG_MARKET}/book`)).length, 2, "each quote reads the book");
  assert.equal(novig.calls.filter((call) => call.url.endsWith(`/events/${NOVIG_EVENT}`)).length, 2, "each quote reads the event");
  assert.match(logs.at(-1), /ok @ 0\.555 .*catalog reused, 59s old/);
  clock += 2 * 1000;
  await runner.novigQuote(new URLSearchParams(WSU_36_LINE));
  assert.equal(catalogReads(), 2, "a catalog over a minute old is read again");
  runner.stop();
  // A reused catalog without the outcome (listed since): read fresh once, then found.
  let listed = false;
  const late = novigFetch();
  const lateImpl = async (url, init) => (url.includes("/catalog/markets?") && !listed
    ? response({ body: JSON.stringify({ items: [] }) }) : late(url, init));
  const { runner: lateRunner } = await quoteRunner({}, { novigFetch: lateImpl, novigBaseUrl: "https://novig.test" });
  const missing = await lateRunner.novigQuote(new URLSearchParams(WSU_36_LINE));
  listed = true;
  const found = await lateRunner.novigQuote(new URLSearchParams(WSU_36_LINE));
  lateRunner.stop();
  assert.deepEqual([missing.ok, missing.reason, found.ok, found.detail], [false, "no_outcome", true, undefined]);
});

test("novig_quote: an event the catalog lists that is not the market's is refused", async () => {
  const novig = novigFetch({ eventId: "another-event" });
  const { runner } = await quoteRunner({}, { novigFetch: novig, novigBaseUrl: "https://novig.test" });
  const quote = await runner.novigQuote(new URLSearchParams(WSU_36_LINE));
  runner.stop();
  assert.deepEqual([quote.ok, quote.reason], [false, "market_mismatch"]);
});

test("taker: a service without the bot key answers 503, the decision reads off", async () => {
  const service = takerService({ "/bets.json": takerBetsBody(), "/settings.json": TAKER_SETTINGS });
  const offService = async (url, init = {}) => (url.endsWith("/taker/order.json")
    ? response({ status: 503, body: JSON.stringify({ status: "off", error: "no bot Kalshi key is configured" }) })
    : service(url, init));
  const runner = runnerLib.createRunner({
    fetchImpl: takerUnabatedFetch(), serviceFetch: offService, kalshiFetch: kalshiFetch(), kalshiBaseUrl: "https://kalshi.test/v2",
    betsServiceUrl: SERVICE_URL, now: () => TAKER_NOW, timers: noTimers, ...quiet,
  });
  await runner.start();
  await runner.takerPass();
  await runner.takerPass();
  assert.equal(runner.takerPayload().decisions[0].order.status, "off");
  assert.match(runner.takerPayload().desk, /no bot Kalshi key/);
  runner.stop();
});

// Cal's open Kalshi bet on the slice's WSU -35.5 market, as /bets.json lists it.
const WSU_HELD = {
  id: `kalshi:${TAKER_TICKER}:yes`, source: "kalshi_api", venue: "kalshi", league: "cfb", status: "open", sourceFetchedAt: "2026-09-15T18:59:00Z",
  awayTeam: "Duquesne", homeTeam: "Washington St.", awayKey: null, homeKey: null, betType: "spread", period: "FG", side: "home", points: -35.5,
  price: -143, stake: 100, toWin: 70, contracts: 100, placedAt: "2026-09-15T18:00:00Z", closedAt: null, eventStart: "2026-09-19T19:30:00.000Z",
  eventDate: "2026-09-19", isParlayLeg: false, parlayId: null, unmatchable: null, approx: [], rotation: null,
  venueIds: { marketTicker: TAKER_TICKER, eventTicker: "KXNCAAFSPREAD-26SEP19DUQWSU" }, raw: {},
};

test("taker: a Removed bet (bet_exclusions) is not Cal's: it neither sizes the stake nor counts as held", async () => {
  async function decisionWith(betsBodyForTaker) {
    const runner = runnerLib.createRunner({
      fetchImpl: takerUnabatedFetch(), serviceFetch: takerService({ "/bets.json": betsBodyForTaker, "/settings.json": TAKER_SETTINGS }),
      kalshiFetch: kalshiFetch(), kalshiBaseUrl: "https://kalshi.test/v2", betsServiceUrl: SERVICE_URL, now: () => TAKER_NOW, timers: noTimers, ...quiet,
    });
    await runner.start();
    await runner.takerPass();
    const decision = runner.takerPayload().decisions.find((row) => row.ticker === TAKER_TICKER);
    runner.stop();
    return decision;
  }
  const held = await decisionWith(takerBetsBody([WSU_HELD]));
  assert.equal(held.reason, "size_zero", held.detail);
  assert.deepEqual(held.heldIds, [WSU_HELD.id]);
  const removed = await decisionWith({ ...takerBetsBody([WSU_HELD]), exclusions: [{ betId: WSU_HELD.id, venue: "kalshi", excludedAt: "2026-09-15T18:30:00Z" }] });
  assert.equal(removed.reason, "bet", removed.detail);
  assert.deepEqual(removed.heldIds, []);
});

test("taker: a decision chunk the bets service refuses (400) is dropped and counted, not retried forever", async () => {
  const service = takerService({ "/bets.json": takerBetsBody(), "/settings.json": TAKER_SETTINGS });
  let refused = 0;
  const refusingService = async (url, init = {}) => {
    if (url.endsWith("/taker_decisions.json") && refused === 0) {
      refused += 1;
      return response({ status: 400, body: JSON.stringify({ error: "rows[0].reason must be one of ..." }) });
    }
    return service(url, init);
  };
  const runner = runnerLib.createRunner({
    fetchImpl: takerUnabatedFetch(), serviceFetch: refusingService, kalshiFetch: kalshiFetch(), kalshiBaseUrl: "https://kalshi.test/v2",
    betsServiceUrl: SERVICE_URL, now: () => TAKER_NOW, timers: noTimers, ...quiet,
  });
  await runner.start();
  await runner.takerPass();
  const log = runner.health().taker.log;
  assert.equal(refused, 1);
  assert.ok(log.dropped >= 1, JSON.stringify(log));
  assert.deepEqual({ unsent: log.unsent, error: log.error }, { unsent: 0, error: null });
  runner.stop();
});

test("taker: holdings unknown (bets service down) decides nothing and reads no book", async () => {
  const kalshi = kalshiFetch();
  const service = takerService({ "/bets.json": new Error("connection refused"), "/settings.json": TAKER_SETTINGS });
  const runner = runnerLib.createRunner({
    fetchImpl: takerUnabatedFetch(), serviceFetch: service, kalshiFetch: kalshi, kalshiBaseUrl: "https://kalshi.test/v2",
    betsServiceUrl: SERVICE_URL, now: () => TAKER_NOW, timers: noTimers, ...quiet,
  });
  await runner.start();
  await runner.takerPass();
  assert.match(runner.takerPayload().pass.skipped, /bets service/);
  assert.deepEqual(runner.takerPayload().decisions, []);
  assert.equal(kalshi.calls.length, 0);
  runner.stop();
});

test("taker: a Kalshi poll that failed, or is older than 5 min, means unknown holdings: the pass decides nothing", async () => {
  for (const [source, expected] of [
    [{ ok: false, error: "HTTP 401" }, /Kalshi poll failed: HTTP 401/],
    [{ fetchedAt: new Date(TAKER_NOW - 6 * 60 * 1000).toISOString() }, /older than 5 min/],
  ]) {
    const kalshi = kalshiFetch();
    const service = takerService({ "/bets.json": takerBetsBody([], source), "/settings.json": TAKER_SETTINGS });
    const runner = runnerLib.createRunner({
      fetchImpl: takerUnabatedFetch(), serviceFetch: service, kalshiFetch: kalshi, kalshiBaseUrl: "https://kalshi.test/v2",
      betsServiceUrl: SERVICE_URL, now: () => TAKER_NOW, timers: noTimers, ...quiet,
    });
    await runner.start();
    await runner.takerPass();
    assert.match(runner.takerPayload().pass.skipped, expected);
    assert.equal(kalshi.calls.length, 0);
    runner.stop();
  }
  // No Kalshi source at all (the service has no Kalshi key) is the same.
  const service = takerService({ "/bets.json": betsBody([]), "/settings.json": TAKER_SETTINGS });
  const runner = runnerLib.createRunner({
    fetchImpl: takerUnabatedFetch(), serviceFetch: service, kalshiFetch: kalshiFetch(), kalshiBaseUrl: "https://kalshi.test/v2",
    betsServiceUrl: SERVICE_URL, now: () => TAKER_NOW, timers: noTimers, ...quiet,
  });
  await runner.start();
  await runner.takerPass();
  assert.match(runner.takerPayload().pass.skipped, /does not poll Kalshi/);
  runner.stop();
});

test("taker: a league snapshot that landed over 3 min ago is stale_feed for its rows, and the runner re-downloads every minute", async () => {
  let clock = TAKER_NOW;
  const kalshi = kalshiFetch();
  const service = takerService({ "/bets.json": takerBetsBody(), "/settings.json": TAKER_SETTINGS });
  service.posts.length = 0;
  const runner = runnerLib.createRunner({
    fetchImpl: takerUnabatedFetch(), serviceFetch: service, kalshiFetch: kalshi, kalshiBaseUrl: "https://kalshi.test/v2",
    betsServiceUrl: SERVICE_URL, now: () => clock, timers: noTimers, ...quiet,
  });
  await runner.start();
  assert.equal(runner.scanner.refreshEveryMs(9.7 * 1024 * 1024), 60 * 1000);
  clock += 4 * 60 * 1000;
  // The Kalshi poll must still be fresh for the pass to run at all.
  service.script["/bets.json"] = takerBetsBody([], { fetchedAt: new Date(clock - 30 * 1000).toISOString() });
  await runner.pollBets();
  await runner.takerPass();
  const decisions = runner.takerPayload().decisions;
  assert.ok(decisions.length > 0);
  assert.ok(decisions.every((row) => row.reason === "stale_feed"), JSON.stringify(decisions.map((row) => row.reason)));
  assert.match(decisions[0].detail, /landed 240 s ago/);
  assert.equal(kalshi.calls.length, 0, "no orderbook is read for a stale row");
  runner.stop();
});

test("taker: a failed orderbook read is no_book for that row, the pass goes on", async () => {
  const service = takerService({ "/bets.json": takerBetsBody(), "/settings.json": TAKER_SETTINGS });
  const runner = runnerLib.createRunner({
    fetchImpl: takerUnabatedFetch(), serviceFetch: service, kalshiFetch: async () => { throw new Error("socket hang up"); },
    kalshiBaseUrl: "https://kalshi.test/v2", betsServiceUrl: SERVICE_URL, now: () => TAKER_NOW, timers: noTimers, ...quiet,
  });
  await runner.start();
  await runner.takerPass();
  const decision = runner.takerPayload().decisions.find((row) => row.ticker === TAKER_TICKER);
  assert.equal(decision.reason, "no_book");
  assert.match(decision.detail, /socket hang up/);
  runner.stop();
});

// ---- the automated taker on Novig (taker.js VENUE_NOVIG) ----
// The same slice with a 6% edge on Novig's WSU -35.5 main line (-127, Novig's
// own 0.56) and Kalshi's edges left as they are; Novig's catalog lists the
// outcome Unabated's rung names on a two-outcome SPREAD market.
const NOVIG_WSU_OUTCOME = "01a0a05e-e864-7b73-9b80-9e502e4e98d5";
const NOVIG_DUQ_OUTCOME = "01a0a05e-e864-7b73-9b80-0000000000aa";
const NOVIG_MARKET = "0b0b0b0b-1111-4222-8333-444444444444";
const NOVIG_EVENT = "0e0e0e0e-1111-4222-8333-555555555555";

function novigTakerUnabatedFetch() {
  const json = JSON.parse(takerVenueSlice());
  const spreadRow = json.odds["lg2:pt1:pregame"].find((row) => row.key === "pt1:pregame:bt2:e123742");
  Object.assign(spreadRow.sides["si1:tid717"].ms89, { statusId: 1, ge: 0.06, bacr: -200, modifiedOn: new Date(TAKER_NOW - 60000).toISOString() });
  const body = JSON.stringify(json);
  return async (url) => {
    if (url.split("?")[0] === SNAPSHOT_BASE_URL(2)) return response({ body, headers: { "content-length": "1000", "last-modified": new Date(TAKER_NOW - 20000).toUTCString() } });
    return response({ status: 404, body: "not found" });
  };
}

function novigFetch(marketFields = {}) {
  const calls = [];
  const market = {
    marketId: NOVIG_MARKET, eventId: NOVIG_EVENT, marketType: "SPREAD", strike: "-35.5", status: "OPEN",
    fee: { coefficient: "0.06", makerCredit: "0", charged: "WHEN_LIVE" },
    outcomes: [{ outcomeId: NOVIG_WSU_OUTCOME, name: "WSU -35.5", status: "TBD" }, { outcomeId: NOVIG_DUQ_OUTCOME, name: "DUQ +35.5", status: "TBD" }],
    ...marketFields,
  };
  const impl = async (url, init = {}) => {
    calls.push({ url, userAgent: init.headers && init.headers["User-Agent"] });
    const bare = url.split("?")[0];
    if (bare === "https://novig.test/v3/public/catalog/markets") return response({ body: JSON.stringify({ items: [market] }) });
    if (bare === `https://novig.test/v3/public/catalog/markets/${NOVIG_MARKET}/book`) {
      return response({ body: JSON.stringify({ marketId: NOVIG_MARKET, seq: 9, orders: { [NOVIG_DUQ_OUTCOME]: [{ orderId: "r1", price: "0.445", qty: 500000 }], [NOVIG_WSU_OUTCOME]: [] } }) });
    }
    if (bare === `https://novig.test/v3/public/catalog/markets/${NOVIG_MARKET}`) return response({ body: JSON.stringify(market) });
    if (bare === `https://novig.test/v3/public/catalog/events/${NOVIG_EVENT}`) return response({ body: JSON.stringify({ eventId: NOVIG_EVENT, status: "OPEN_PREGAME" }) });
    return response({ status: 404, body: "{}" });
  };
  impl.calls = calls;
  return impl;
}

function novigTakerBetsBody(novigSources = {}) {
  const body = takerBetsBody();
  const fresh = { ok: true, fetchedAt: new Date(TAKER_NOW - 30 * 1000).toISOString(), error: null, count: 0 };
  body.sources.novig = { ...fresh, ...novigSources.novig };
  body.sources.novig_bot = { ...fresh, ...novigSources.novig_bot };
  return body;
}

function novigTakerService(betsBodyForTaker, deskReply = null) {
  const service = takerService({ "/bets.json": betsBodyForTaker, "/settings.json": TAKER_SETTINGS });
  const novigOrders = [];
  const impl = async (url, init = {}) => {
    if (url.endsWith("/taker/novig_order.json")) {
      const body = JSON.parse(init.body);
      novigOrders.push(body);
      if (deskReply) return response(deskReply);
      return response({ body: JSON.stringify({ status: "filled", qty: body.qty, filled: body.qty, cost: (body.qty * body.price) / 100, orderId: "n-1" }) });
    }
    return service(url, init);
  };
  impl.novigOrders = novigOrders;
  impl.posts = service.posts;
  return impl;
}

test("taker: a Novig edge at the live ask is sent to the Novig desk in 1c contracts, with the board line, after its market is read", async () => {
  const novig = novigFetch();
  const service = novigTakerService(novigTakerBetsBody());
  const runner = runnerLib.createRunner({
    fetchImpl: novigTakerUnabatedFetch(), serviceFetch: service, kalshiFetch: kalshiFetch(), kalshiBaseUrl: "https://kalshi.test/v2",
    novigFetch: novig, novigBaseUrl: "https://novig.test", betsServiceUrl: SERVICE_URL, now: () => TAKER_NOW, timers: noTimers, ...quiet,
  });
  await runner.start();
  await runner.takerPass();
  const payload = runner.takerPayload();
  assert.equal(payload.pass.novigSkipped, null);
  const decision = payload.decisions.find((row) => row.venue === "novig");
  assert.equal(decision.reason, "bet", decision.detail);
  assert.equal(decision.liveAsk, 0.555);
  assert.equal(decision.order.status, "filled");
  assert.equal(payload.novigDesk, "armed");
  const [order] = service.novigOrders;
  assert.match(order.clientId, /^[0-9a-f]{8}-[0-9a-f]{4}-4[0-9a-f]{3}-8[0-9a-f]{3}-[0-9a-f]{12}$/);
  assert.deepEqual({ outcomeId: order.outcomeId, marketId: order.marketId, price: order.price, novigBaseUrl: order.novigBaseUrl },
    { outcomeId: NOVIG_WSU_OUTCOME, marketId: NOVIG_MARKET, price: 0.555, novigBaseUrl: "https://novig.test" });
  assert.equal(order.qty, decision.contracts);
  assert.equal(order.qty, Math.floor((decision.stake * 100) / 0.555));
  assert.deepEqual({ league: order.line.league, side: order.line.side, points: order.line.points, homeTeamId: order.line.homeTeamId },
    { league: "cfb", side: "home", points: -35.5, homeTeamId: "717" });
  // Kalshi's rows still decided and sent in the same pass, on the same list.
  assert.ok(payload.decisions.some((row) => row.venue === "kalshi" && row.reason === "bet"));
  assert.ok(novig.calls.some((call) => call.url === `https://novig.test/v3/public/catalog/events/${NOVIG_EVENT}`), "the event is read before the order");
  assert.ok(novig.calls.every((call) => call.userAgent === "unabated-ticket-runner/1.0"));
  runner.stop();
});

test("taker: Novig decides nothing while the bot subaccount or Cal's Novig read is stale, and Kalshi carries on", async () => {
  for (const [sources, expected] of [
    [{ novig_bot: { ok: false, error: "HTTP 401" } }, /novig_bot read failed: HTTP 401/],
    [{ novig: { fetchedAt: new Date(TAKER_NOW - 6 * 60 * 1000).toISOString() } }, /novig read is older than 5 min/],
  ]) {
    const novig = novigFetch();
    const service = novigTakerService(novigTakerBetsBody(sources));
    const runner = runnerLib.createRunner({
      fetchImpl: novigTakerUnabatedFetch(), serviceFetch: service, kalshiFetch: kalshiFetch(), kalshiBaseUrl: "https://kalshi.test/v2",
      novigFetch: novig, novigBaseUrl: "https://novig.test", betsServiceUrl: SERVICE_URL, now: () => TAKER_NOW, timers: noTimers, ...quiet,
    });
    await runner.start();
    await runner.takerPass();
    const payload = runner.takerPayload();
    assert.match(payload.pass.novigSkipped, expected);
    assert.equal(payload.pass.skipped, null);
    assert.equal(novig.calls.length, 0);
    assert.equal(service.novigOrders.length, 0);
    assert.ok(payload.decisions.every((row) => row.venue === "kalshi"));
    runner.stop();
  }
});

test("taker: a Novig market charging a pregame fee is market_mismatch and nothing is sent", async () => {
  const service = novigTakerService(novigTakerBetsBody());
  const runner = runnerLib.createRunner({
    fetchImpl: novigTakerUnabatedFetch(), serviceFetch: service, kalshiFetch: kalshiFetch(), kalshiBaseUrl: "https://kalshi.test/v2",
    novigFetch: novigFetch({ fee: { charged: "ALWAYS" } }), novigBaseUrl: "https://novig.test", betsServiceUrl: SERVICE_URL,
    now: () => TAKER_NOW, timers: noTimers, ...quiet,
  });
  await runner.start();
  await runner.takerPass();
  const decision = runner.takerPayload().decisions.find((row) => row.venue === "novig");
  assert.equal(decision.reason, "market_mismatch");
  assert.match(decision.detail, /fee/);
  assert.equal(service.novigOrders.length, 0);
  runner.stop();
});

test("taker: /taker.json says what the Novig desk answered: off, stopped or out of location", async () => {
  for (const [deskReply, expected] of [
    [{ status: 503, body: JSON.stringify({ error: "no NOVIG_TAKER_KEY_ID" }) }, /^off: /],
    [{ body: JSON.stringify({ status: "stopped", detail: "stop file /data/novig_taker.stop is present" }) }, /^stopped: stop file/],
    [{ body: JSON.stringify({ status: "geo_blocked", detail: "Novig HTTP 451 GEOLOCATION_EXPIRED" }) }, /^geo_blocked: Novig HTTP 451/],
  ]) {
    const service = novigTakerService(novigTakerBetsBody(), deskReply);
    const runner = runnerLib.createRunner({
      fetchImpl: novigTakerUnabatedFetch(), serviceFetch: service, kalshiFetch: kalshiFetch(), kalshiBaseUrl: "https://kalshi.test/v2",
      novigFetch: novigFetch(), novigBaseUrl: "https://novig.test", betsServiceUrl: SERVICE_URL, now: () => TAKER_NOW, timers: noTimers, ...quiet,
    });
    await runner.start();
    await runner.takerPass();
    const payload = runner.takerPayload();
    assert.equal(service.novigOrders.length, 1);
    assert.match(payload.novigDesk, expected);
    runner.stop();
  }
});

test("taker: a stale Kalshi poll stops the Novig side too (its stake is sized on Cal's Kalshi bets)", async () => {
  const body = novigTakerBetsBody();
  body.sources.kalshi = { ...body.sources.kalshi, fetchedAt: new Date(TAKER_NOW - 6 * 60 * 1000).toISOString() };
  const novig = novigFetch();
  const service = novigTakerService(body);
  const runner = runnerLib.createRunner({
    fetchImpl: novigTakerUnabatedFetch(), serviceFetch: service, kalshiFetch: kalshiFetch(), kalshiBaseUrl: "https://kalshi.test/v2",
    novigFetch: novig, novigBaseUrl: "https://novig.test", betsServiceUrl: SERVICE_URL, now: () => TAKER_NOW, timers: noTimers, ...quiet,
  });
  await runner.start();
  await runner.takerPass();
  const payload = runner.takerPayload();
  assert.match(payload.pass.novigSkipped, /Kalshi poll is older than 5 min/);
  assert.equal(novig.calls.length, 0);
  assert.equal(service.novigOrders.length, 0);
  runner.stop();
});

// ---- The board thread (server/board_thread.js) ----
// A worker stand-in: the board side served over a real MessageChannel, so the
// messages are cloned as between threads; exit() plays the thread dying.
const { MessageChannel } = require("node:worker_threads");
const { EventEmitter } = require("node:events");
const boardThread = require("../server/board_thread.js");

function channelWorker(runner) {
  const { port1, port2 } = new MessageChannel();
  boardThread.serveBoardCalls(port2, runner, () => {});
  const worker = new EventEmitter();
  port1.on("message", (message) => worker.emit("message", message));
  worker.postMessage = (message) => port1.postMessage(message);
  worker.terminate = () => worker.exit(1);
  worker.exit = (code) => {
    port1.close();
    port2.close();
    worker.emit("exit", code);
  };
  return worker;
}

function boardClient(runners, clock, extra = {}) {
  const spawned = [];
  const client = boardThread.createBoardClient({
    now: () => clock.now, ...quiet, setTimeout: extra.setTimeout || (() => {}),
    createWorker: () => {
      const worker = channelWorker(runners[Math.min(spawned.length, runners.length - 1)]);
      spawned.push(worker);
      return worker;
    },
  });
  client.start();
  return { client, spawned };
}

// A runner whose payloads count their builds and say which build they are.
function countingRunner() {
  const builds = { edgesPayload: 0, teasersPayload: 0 };
  const runner = {
    builds, betsPolls: 0,
    edgesPayload: () => ({ build: (builds.edgesPayload += 1) }),
    teasersPayload: () => ({ build: (builds.teasersPayload += 1) }),
    scenariosPayload: () => ({ games: [] }),
    takerPayload: () => ({ decisions: [] }),
    pollBets: async () => { runner.betsPolls += 1; },
    attachPayload: () => { const error = new Error("no open bet \"x\""); error.status = 404; throw error; },
    kalshiQuote: async (params) => ({ ok: false, eventId: Number(params.get("eventId")) }),
    novigQuote: async () => { throw new Error("boom"); },
  };
  return runner;
}

test("board thread: the HTTP routes answer exactly what the runner builds, through the thread's messages", async () => {
  const { runner } = await startedRunner({ "/bets.json": betsBody([BEARS_HELD]), "/settings.json": { settings: OPEN_SETTINGS, updatedAt: null } });
  const clock = { now: NOW };
  const { client } = boardClient([runner], clock);
  const server = runnerLib.createHttpServer(client, { host: "127.0.0.1" });
  await new Promise((resolve) => server.listen(0, "127.0.0.1", resolve));
  const port = server.address().port;
  const get = (route) => new Promise((resolve, reject) => {
    const request = require("node:http").request({ host: "127.0.0.1", port, path: route, method: "GET" }, (reply) => {
      let body = "";
      reply.on("data", (chunk) => { body += chunk; });
      reply.on("end", () => resolve({ status: reply.statusCode, body: JSON.parse(body) }));
    });
    request.on("error", reject);
    request.end();
  });
  const same = (value) => JSON.parse(JSON.stringify(value));
  try {
    assert.deepEqual((await get("/edges.json")).body, same(runner.edgesPayload()));
    assert.deepEqual((await get("/edges.json?fresh=1")).body, same(runner.edgesPayload()));
    assert.deepEqual((await get("/teasers.json")).body, same(runner.teasersPayload()));
    assert.deepEqual((await get("/scenarios.json")).body, same(runner.scenariosPayload()));
    assert.equal((await get("/attach.json?betId=bol-1&query=bears")).body.events[0].eventId, 125807);
    const missing = await get("/attach.json?betId=nope");
    assert.deepEqual([missing.status, missing.body.error], [404, 'no open bet "nope"']);
    const health = await get("/health");
    assert.deepEqual([health.status, health.body.boardThread.state, health.body.scanner.phase], [200, "running", "starting"]);
  } finally {
    await new Promise((resolve) => server.close(resolve));
    client.stop();
    runner.stop();
  }
});

test("board thread: a copy under a minute old is served at once and rebuilt for the next ask; older, or fresh, waits for a build", async () => {
  const runner = countingRunner();
  const clock = { now: NOW };
  const { client } = boardClient([runner], clock);
  const read = async (opts) => JSON.parse((await client.edgesPayload(opts)).text).build;
  assert.equal(await read(), 1);
  assert.equal(await read(), 1);
  await new Promise((resolve) => setImmediate(resolve));
  await new Promise((resolve) => setImmediate(resolve));
  assert.equal(runner.builds.edgesPayload, 2);
  assert.equal(await read(), 2);
  await new Promise((resolve) => setImmediate(resolve));
  await new Promise((resolve) => setImmediate(resolve));
  clock.now += boardThread.PAYLOAD_SERVE_MAX_AGE_MS + 1;
  assert.equal(await read(), 4);
  assert.equal(await read({ fresh: true }), 5);
  assert.equal(runner.betsPolls, 1);
  client.stop();
});

test("board thread: asks at once share one build; errors keep their status, a bug is a 500", async () => {
  const runner = countingRunner();
  const { client } = boardClient([runner], { now: NOW });
  const both = await Promise.all([client.teasersPayload(), client.teasersPayload()]);
  assert.deepEqual(both.map((reply) => JSON.parse(reply.text).build), [1, 1]);
  assert.equal(runner.builds.teasersPayload, 1);
  await assert.rejects(client.attachPayload(new URLSearchParams("betId=x")), (error) => error.status === 404 && /no open bet/.test(error.message));
  await assert.rejects(client.novigQuote(new URLSearchParams("eventId=1")), (error) => error.status === 500 && error.message === "boom");
  assert.deepEqual(JSON.parse((await client.kalshiQuote(new URLSearchParams("eventId=7"))).text), { ok: false, eventId: 7 });
  client.stop();
});

test("board thread: a thread that dies fails what waits with 503, drops its copies, and a new thread starts after the pause", async () => {
  const first = countingRunner();
  let release;
  first.kalshiQuote = () => new Promise((resolve) => { release = resolve; });
  const second = countingRunner();
  const restarts = [];
  const { client, spawned } = boardClient([first, second], { now: NOW }, { setTimeout: (fn, ms) => restarts.push({ fn, ms }) });
  await client.edgesPayload();
  const waiting = client.kalshiQuote(new URLSearchParams("eventId=1"));
  await new Promise((resolve) => setImmediate(resolve));
  spawned[0].exit(1);
  await assert.rejects(waiting, (error) => error.status === 503);
  await assert.rejects(client.edgesPayload(), (error) => error.status === 503 && /restarting/.test(error.message));
  assert.equal(client.health().boardThread.state, "restarting");
  assert.deepEqual(restarts.map((entry) => entry.ms), [boardThread.BOARD_RESTART_DELAY_MS]);
  restarts[0].fn();
  assert.equal(JSON.parse((await client.edgesPayload()).text).build, 1);
  assert.deepEqual([second.builds.edgesPayload, client.health().boardThread.restarts, client.health().boardThread.lastExit.code], [1, 1, 1]);
  release({});
  client.stop();
  spawned[1].exit(1);
  assert.equal(restarts.length, 1);
});

test("board thread: the real worker boots board_thread.js with runner.js and posts its health", async () => {
  // Nothing listens on port 9: every read fails fast, which is all a boot needs.
  const dead = "http://127.0.0.1:9";
  const { Worker } = require("node:worker_threads");
  const worker = new Worker(path.join(__dirname, "..", "server", "board_thread.js"), {
    workerData: { role: "board", config: { betsServiceUrl: dead, kalshiBaseUrl: dead, novigBaseUrl: dead } }, stdout: true, stderr: true,
  });
  try {
    const health = await new Promise((resolve, reject) => {
      worker.once("error", reject);
      worker.on("message", (message) => { if (message.type === "health") resolve(message.health); });
    });
    assert.deepEqual([health.ok, health.scanner.phase], [true, "starting"]);
  } finally {
    await worker.terminate();
  }
});
