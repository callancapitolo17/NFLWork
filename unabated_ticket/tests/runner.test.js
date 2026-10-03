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

test("HTTP: the phone page's four files under the CSP; the Host rule and 405 cover them; nothing else on disk is reachable", async () => {
  const { runner } = await startedRunner({ "/bets.json": betsBody([]), "/settings.json": { settings: OPEN_SETTINGS, updatedAt: null } });
  const server = runnerLib.createHttpServer(runner, { host: "127.0.0.1" });
  await new Promise((resolve) => server.listen(0, "127.0.0.1", resolve));
  const port = server.address().port;
  const get = (route, options) => new Promise((resolve, reject) => {
    const request = require("node:http").request({ host: "127.0.0.1", port, path: route, method: "GET", ...options }, (reply) => {
      let body = "";
      reply.on("data", (chunk) => { body += chunk; });
      reply.on("end", () => resolve({ status: reply.statusCode, headers: reply.headers, body }));
    });
    request.on("error", reject);
    request.end();
  });
  try {
    const page = await get("/");
    assert.equal(page.status, 200);
    assert.match(page.headers["content-type"], /^text\/html/);
    assert.equal(page.headers["content-security-policy"], runnerLib.PHONE_CSP);
    assert.equal(page.headers["cache-control"], "no-store");
    assert.match(page.body, /<script src="phone_view.js"><\/script>/);
    for (const [route, type] of [["/phone.css", /^text\/css/], ["/phone_view.js", /^text\/javascript/], ["/phone.js", /^text\/javascript/]]) {
      const asset = await get(route);
      assert.deepEqual([asset.status, asset.headers["x-content-type-options"]], [200, "nosniff"], route);
      assert.match(asset.headers["content-type"], type, route);
    }
    assert.equal((await get("/", { headers: { Host: "evil.example" } })).status, 403);
    assert.equal((await get("/", { method: "POST" })).status, 405);
    for (const route of ["/index.html", "/runner.js", "/../runner.js", "/phone/phone.js", "/%2e%2e/runner.js"]) {
      assert.equal((await get(route)).status, 404, route);
    }
  } finally {
    await new Promise((resolve) => server.close(resolve));
    runner.stop();
  }
});

test("the phone page's files load only from this origin: no URL in them names another host", () => {
  const dir = path.join(__dirname, "..", "server", "phone");
  for (const file of ["index.html", "phone.css", "phone.js", "phone_view.js"]) {
    const text = fs.readFileSync(path.join(dir, file), "utf8");
    assert.doesNotMatch(text, /https?:\/\//, file);
    assert.doesNotMatch(text, /innerHTML|outerHTML|insertAdjacentHTML|document\.write/, file);
  }
  const loaded = runnerLib.loadPhoneFiles();
  assert.deepEqual(Object.keys(loaded).sort(), ["/", "/phone.css", "/phone.js", "/phone_view.js"]);
});

// ---- the phone page's view (server/phone/phone_view.js) ---------------------

const phoneView = require("../server/phone/phone_view.js");

test("phone view: the runner's real payload becomes cards carrying its own stake, edge and words — nothing re-derived", async () => {
  const { runner } = await startedRunner({ "/bets.json": betsBody([BEARS_HELD]), "/settings.json": { settings: OPEN_SETTINGS, updatedAt: null } });
  try {
    const payload = JSON.parse(JSON.stringify(runner.edgesPayload()));
    assert.ok(payload.items.length > 0, "the slice lists edges");
    const page = phoneView.pageView(payload, Date.parse(payload.generatedAt) + 12 * 1000, "America/Los_Angeles");
    assert.equal(page.updated, "updated 12s ago");
    assert.equal(page.stale, false);
    assert.equal(page.cards.length, payload.items.length);
    payload.items.forEach((card, index) => {
      const shown = page.cards[index];
      assert.equal(shown.stake, card.best.rail.text);
      assert.equal(shown.stakeNote, card.best.rail.note);
      assert.equal(shown.edge, `${card.best.edgePct.toFixed(2)}%`);
      assert.equal(shown.matchup, `${card.awayTeam} @ ${card.homeTeam}`);
      assert.equal(shown.badges.length, card.best.badges.length);
    });
    const bears = page.cards.find((card) => card.matchup === "Chicago Bears @ Carolina Panthers" && card.header.includes("Spread"));
    assert.ok(bears, "the held Bears spread is listed");
    assert.match(bears.when, /^Sun Sep 13 · 10:00 AM$/);
    // The slice serves NFL only, so CFB (loaded for the teasers) fails: a warning, not an alarm.
    assert.deepEqual(page.chips.map((chip) => chip.tone), ["warn", "ok", "plain"]);
    assert.match(page.chips[0].text, /^Feed ok · 1 leagues · feed unavailable for /);
  } finally {
    runner.stop();
  }
});

const MINIMAL_PAYLOAD = {
  generatedAt: "2026-10-03T17:26:10Z", grouped: true, unit: "cards", total: 0, items: [], tailFlex: "",
  scanner: { phase: "live", error: null, leaguesLoaded: [1, 3], leagueErrors: {}, snapshotBuiltAt: Date.parse("2026-10-03T17:25:00Z"), loading: { done: 2, total: 2 } },
  betsService: { okAt: Date.parse("2026-10-03T17:26:00Z"), error: null, unreachableSince: null, openBets: 210, sources: { kalshi: { ok: true }, novig: { ok: true } } },
  settings: { source: "service", error: null, stake: { bankroll: 30000, multiplier: 0.25 }, edges: { minEdgePct: 2.5 } },
};
const AT = Date.parse("2026-10-03T17:26:22Z");

test("phone view: the status chips name what is wrong — bets down is the loudest, a failing venue and a slow feed warn", () => {
  const healthy = phoneView.pageView(MINIMAL_PAYLOAD, AT);
  assert.deepEqual(healthy.chips.map((chip) => chip.text), ["Feed ok · 2 leagues", "Bets ok · 210 open", "$30,000 · ¼ Kelly · 2.5% min"]);
  assert.equal(healthy.empty, "No edges pass your filters right now.");
  assert.equal(healthy.summary, "0 cards · pregame");

  const betsDown = phoneView.pageView({ ...MINIMAL_PAYLOAD, betsService: { ...MINIMAL_PAYLOAD.betsService, error: "HTTP 500", unreachableSince: AT - 5 * 60 * 1000 } }, AT);
  assert.deepEqual(betsDown.chips[1], { text: "Bets service down 5 min ago — stakes use bets as of then", tone: "bad" });
  const neverRead = phoneView.pageView({ ...MINIMAL_PAYLOAD, betsService: { okAt: null, error: "connect ECONNREFUSED" } }, AT);
  assert.deepEqual(neverRead.chips[1], { text: "Bets service not read — stakes ignore open bets", tone: "bad" });
  const venueFailing = phoneView.pageView({ ...MINIMAL_PAYLOAD, betsService: { ...MINIMAL_PAYLOAD.betsService, sources: { kalshi: { ok: true }, novig: { ok: false } } } }, AT);
  assert.deepEqual(venueFailing.chips[1], { text: "Bets ok · novig failing", tone: "warn" });

  const loading = phoneView.pageView({ ...MINIMAL_PAYLOAD, scanner: { ...MINIMAL_PAYLOAD.scanner, loading: { done: 1, total: 2 } } }, AT);
  assert.deepEqual(loading.chips[0], { text: "Feed loading 1/2", tone: "warn" });
  const oldFeed = phoneView.pageView({ ...MINIMAL_PAYLOAD, scanner: { ...MINIMAL_PAYLOAD.scanner, snapshotBuiltAt: AT - 15 * 60 * 1000 } }, AT);
  assert.deepEqual(oldFeed.chips[0], { text: "Feed old · built 15 min ago", tone: "warn" });
  const feedError = phoneView.pageView({ ...MINIMAL_PAYLOAD, scanner: { ...MINIMAL_PAYLOAD.scanner, leaguesLoaded: [], error: "HTTP 403" } }, AT);
  assert.deepEqual(feedError.chips[0], { text: "Feed error: HTTP 403", tone: "bad" });
  const partly = phoneView.pageView({ ...MINIMAL_PAYLOAD, scanner: { ...MINIMAL_PAYLOAD.scanner, error: "feed unavailable for CFB" } }, AT);
  assert.deepEqual(partly.chips[0], { text: "Feed ok · 2 leagues · feed unavailable for CFB", tone: "warn" });
});

test("phone view: a list older than two minutes says so; an ungrouped list shows each line as a card of one", () => {
  const stale = phoneView.pageView(MINIMAL_PAYLOAD, Date.parse(MINIMAL_PAYLOAD.generatedAt) + 3 * 60 * 1000);
  assert.deepEqual([stale.updated, stale.stale], ["data from 3 min ago", true]);

  const line = {
    key: "k1", sideName: "Under", leagueLabel: "NFL", betType: "Total", period: "1H", eventStartMs: Date.parse("2026-10-04T17:00:00Z"),
    awayTeam: "New England Patriots", homeTeam: "Buffalo Bills", sideLabel: "Under 24.5", book: { id: 89, name: "Novig" }, price: 122, fair: -105,
    liquidity: 997.4, edgePct: 3.26, edgeTier: "warm", move: { kind: "book_away", label: "book moved away" },
    rail: { text: "add $88.10", note: "$200.41 alone", atSize: false },
    badges: [{ kind: "held", text: "held $112" }, { kind: "against", text: "against" }],
    related: [{ text: "Buffalo Bills -14.5 +292 · $390 · Novig", tag: "game · not sized" }],
  };
  const page = phoneView.pageView({ ...MINIMAL_PAYLOAD, grouped: false, unit: "lines", total: 5, items: [line] }, AT, "America/Los_Angeles");
  assert.equal(page.summary, "5 lines (top 1) · pregame");
  assert.deepEqual(page.cards[0], {
    key: "k1", header: "NFL · Total · 1st half", when: "Sun Oct 4 · 10:00 AM", matchup: "New England Patriots @ Buffalo Bills",
    badges: [{ text: "held $112", tone: "warn" }, { text: "against", tone: "bad" }],
    side: "Under 24.5", priceLine: "Novig +122 · fair −105", liquidity: "liq to win $997", edge: "3.26%", edgeTier: "warm",
    moveLabel: "book moved away", related: ["Buffalo Bills -14.5 +292 · $390 · Novig — game · not sized"],
    stake: "add $88.10", stakeNote: "$200.41 alone", atSize: false, othersText: "no other lines",
  });
});
