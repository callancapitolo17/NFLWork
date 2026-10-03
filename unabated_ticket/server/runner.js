#!/usr/bin/env node
// Unabated Ticket server runner — the panel's Edges scan, headless, for the
// phone page (plan: phone_page_plan.md, step 1).
//
//   node unabated_ticket/server/runner.js        (Node 18+, no build, no npm deps)
//
// Loads the extension's own pure modules unchanged (scanner.js, feed.js,
// edgerows.js, tailflex.js, teaser.js, betsview.js, bets.js, teams.js,
// ladder.js, condkelly.js, kelly.js, edgemove.js, fillfair.js) and runs the
// scan loop the panel runs while it is open and visible: a snapshot per
// league on start — the settings' leagues plus NFL and CFB, which the panel
// always loads for its Teasers tab and whose open BFA teasers size the Edges
// rows — each league re-downloaded on its own size-based cadence, a full
// resync every 10 min (scanner.js). It never pauses — the panel pauses when
// hidden; this is a panel that is always in view.
//
// Inputs (env):
//   UNABATED_RUNNER_HOST   bind address, default 127.0.0.1 (loopback only;
//                          the Tailscale address comes with the deploy step)
//   UNABATED_RUNNER_PORT   default 8095
//   BETS_SERVICE_URL       the bets service, default http://127.0.0.1:8094
// Network:
//   GET content.unabated.com league snapshots (public, no login), as the panel does.
//   GET <bets service>/bets.json every 30 s (the panel's cadence) — open bets
//       for conditional Kelly, the team crosswalk, pins and saved fill fairs.
//   GET <bets service>/settings.json every 10 s — bankroll, Kelly multiplier,
//       books, minimum edge, line age, min liq to win, alts, sort, cards
//       (bets.duckdb::edge_settings; null fields are the panel's defaults).
// Serves (HTTP on UNABATED_RUNNER_HOST:UNABATED_RUNNER_PORT, no auth):
//   GET /            the phone page (server/phone/, step 2): the Edges list
//                    read-only, fetched from /edges.json every 30 s; its
//                    three files are read once at start and served with a
//                    CSP that allows only this origin
//   GET /edges.json  the Edges list the panel would show for those settings
//                    (server/edges_payload.js documents the shape)
//   GET /health      {ok, generatedAt, uptimeSec, scanner, betsService, settings}
//   Every request whose Host header is not a name this runner serves on is
//   refused with 403 (DNS rebinding, #125); any verb but GET is 405.
// Side effects: none on disk — no DuckDB, no files, no writes to the bets
// service. State (feed, line history, bets) is in memory and rebuilt on
// restart. Logs state changes (not every poll) to stdout/stderr.

"use strict";

const fs = require("node:fs");
const http = require("node:http");
const nodePath = require("node:path");
const scannerLib = require("../extension/scanner.js");
const teams = require("../extension/teams.js");
const betsLib = require("../extension/bets.js");
const fillfair = require("../extension/fillfair.js");
const tailflex = require("../extension/tailflex.js");
const teaser = require("../extension/teaser.js");
const edgeRows = require("../extension/edgerows.js");
const edgesPayload = require("./edges_payload.js");

const DEFAULT_HOST = "127.0.0.1";
const DEFAULT_PORT = 8095;
const DEFAULT_BETS_SERVICE_URL = "http://127.0.0.1:8094";
// The panel polls the bets service every 30 s (panel.js BETS_POLL_MS).
const BETS_POLL_MS = 30 * 1000;
// A phone setting change shows within one scanner tick.
const SETTINGS_POLL_MS = 10 * 1000;
const SERVICE_TIMEOUT_MS = 10 * 1000;
const LOOPBACK_HOST_NAMES = ["127.0.0.1", "localhost"];
const WILDCARD_HOSTS = ["0.0.0.0", "::", ""];
const HTTP_DEFAULT_PORT = 80;
const MAX_PORT = 65535;
const HOUR_MS = 3600 * 1000;

function log(message) {
  console.log(`${new Date().toISOString()} [runner] ${message}`);
}

function logError(message) {
  console.error(`${new Date().toISOString()} [runner] ${message}`);
}

// The runner's config from the environment; throws naming the bad variable.
//   {host, port, betsServiceUrl}
function configFromEnv(env) {
  const host = env.UNABATED_RUNNER_HOST || DEFAULT_HOST;
  const rawPort = env.UNABATED_RUNNER_PORT || String(DEFAULT_PORT);
  const port = Number(rawPort);
  if (!Number.isInteger(port) || port < 1 || port > MAX_PORT) {
    throw new Error(`UNABATED_RUNNER_PORT must be a whole number from 1 to ${MAX_PORT}, got ${JSON.stringify(rawPort)}`);
  }
  const betsServiceUrl = (env.BETS_SERVICE_URL || DEFAULT_BETS_SERVICE_URL).replace(/\/+$/, "");
  if (!/^https?:\/\/\S+$/.test(betsServiceUrl)) {
    throw new Error(`BETS_SERVICE_URL must be an http(s) URL, got ${JSON.stringify(env.BETS_SERVICE_URL)}`);
  }
  return { host, port, betsServiceUrl };
}

// The Host values a request to this runner can carry: the loopback names and
// the address it is bound to, each with the port (a browser omits the
// scheme's default port, so on 80 the bare names too). Same rule as the bets
// service (service.py allowed_hosts): a DNS-rebound page sends its own name.
function allowedHosts(bindHost, port) {
  const names = LOOPBACK_HOST_NAMES.slice();
  const bound = String(bindHost).toLowerCase();
  if (!WILDCARD_HOSTS.includes(bound) && !names.includes(bound)) names.push(bound.includes(":") ? `[${bound}]` : bound);
  const withPort = names.map((name) => `${name}:${port}`);
  return port === HTTP_DEFAULT_PORT ? withPort.concat(names) : withPort;
}

function hostAllowed(hostHeader, allowed) {
  return allowed.includes(String(hostHeader || "").trim().toLowerCase());
}

// The leagues the scanner loads: the Edges list's, plus NFL and CFB, as the
// panel's scannerLeaguesOf does — the board the open teasers and the bets are
// matched against must be the panel's, or the stakes would differ.
function scannerLeaguesOf(edgeLeagues) {
  return Array.from(new Set([...edgeLeagues, ...teaser.TEASER_LEAGUE_IDS])).sort((a, b) => a - b);
}

// scanner.js passes the browser's fetch options ({cache: "no-cache",
// credentials: "include"}); Node has no HTTP cache and no cookie jar, so
// only the abort signal means anything here.
function nodeFetch(url, options) {
  return fetch(url, { signal: options && options.signal });
}

// One JSON GET of the bets service. Node's fetch says only "fetch failed" for
// a refused or reset connection; the reason is in its cause, so it is named.
async function getServiceJson(fetchImpl, url) {
  let response;
  try {
    response = await fetchImpl(url, { signal: AbortSignal.timeout(SERVICE_TIMEOUT_MS) });
  } catch (error) {
    const reason = error.cause ? ` (${error.cause.code || error.cause.message})` : "";
    throw new Error(`${url}: ${error.message}${reason}`, { cause: error });
  }
  if (!response.ok) throw new Error(`HTTP ${response.status} from ${url}`);
  return response.json();
}

// The scan loop, the bets poll and the settings poll, held in memory.
//   deps  {fetchImpl (Unabated snapshots), serviceFetch (the bets service),
//          betsServiceUrl, now, timers {setInterval, clearInterval},
//          logInfo / logWarning (message) -> void, default stdout / stderr}
// Returns {start, stop, pollBets, pollSettings, scanLoaded, edgesPayload, health, scanner}.
function createRunner(deps) {
  const now = deps.now || (() => Date.now());
  const timers = deps.timers || { setInterval: (fn, ms) => setInterval(fn, ms), clearInterval: (id) => clearInterval(id) };
  const serviceFetch = deps.serviceFetch || deps.fetchImpl;
  const logInfo = deps.logInfo || log;
  const logWarning = deps.logWarning || logError;
  const startedAt = now();

  let feedState = null;
  let scannerStatus = null;
  let history = {};
  let boardLinesCache = null;
  let ladderReaders = null;
  // The panel's other per-snapshot caches: the tail-flex measurement (dropped
  // on every scanner update, re-measured when the max line age changes) and
  // teaser.teaserBoardOf's pass, redone only when an NFL or CFB snapshot lands.
  let tailFlexCache = null;
  let teaserBoard = null;
  let teaserBoardLoadedAt = null;
  let teamsSpellingCount = teams.spellingCount();
  let held = { records: [], crosswalk: [], pins: [], fillFairs: [] };
  let fillFairIndex = new Map();
  const betsStatus = { okAt: null, error: null, unreachableSince: null, generatedAt: null, sources: {} };
  let settings = { ...edgesPayload.settingsFromService(null), source: "defaults", error: null, okAt: null, updatedAt: null };
  let scannedLeagues = null;
  // The scanner's current full load (scanner.start), resolved when it lands.
  let scanLoad = Promise.resolve();
  let pollTimers = [];
  let betsBusy = false;
  let settingsBusy = false;

  // Every snapshot carries Unabated's team list; register it so bet records
  // resolve to the board's team ids, and re-resolve when it grew (as the panel does).
  function registerFeedTeams() {
    if (!feedState || !feedState.teamIndex) return;
    for (const [league, list] of Object.entries(edgeRows.feedTeamsByLeague(feedState))) teams.registerTeams(league, list);
    const spellings = teams.spellingCount();
    if (spellings === teamsSpellingCount) return;
    teamsSpellingCount = spellings;
    held = { ...held, records: betsLib.resolveTeamKeys(held.records, held.crosswalk) };
  }

  const scanner = scannerLib.createScanner({
    fetchImpl: deps.fetchImpl, now, timers,
    onChange: (status, state, lineHistory) => {
      if (status.error && (!scannerStatus || scannerStatus.error !== status.error)) logWarning(`scanner: ${status.error}`);
      if (status.phase === "live" && (!scannerStatus || scannerStatus.phase !== "live")) {
        logInfo(`scanner live: ${status.leaguesLoaded.length} leagues, ${status.lineCount} lines (+${status.altLineCount} alts)`);
      }
      scannerStatus = status;
      feedState = state;
      history = lineHistory || {};
      boardLinesCache = null;
      ladderReaders = null;
      tailFlexCache = null;
      registerFeedTeams();
    },
  });

  function boardLines() {
    if (!feedState) return [];
    if (!boardLinesCache) boardLinesCache = edgeRows.boardLines(feedState);
    return boardLinesCache;
  }

  function ladderReaderOf(eventId) {
    if (!ladderReaders) ladderReaders = edgeRows.createLadderReaders(feedState);
    return ladderReaders(eventId);
  }

  // Measured once per scanner update, on the same line-age gate as the list (panel.js tailFlexMeasurement).
  function tailFlexMeasurement() {
    const maxLineAgeMs = settings.edgeSettings.maxLineAgeHours * HOUR_MS;
    if (!tailFlexCache || tailFlexCache.maxLineAgeMs !== maxLineAgeMs) {
      tailFlexCache = { maxLineAgeMs, measurement: tailflex.measureTailFlex(feedState, { now: now(), maxLineAgeMs }) };
    }
    return tailFlexCache.measurement;
  }

  // The board pass, redone when an NFL or CFB snapshot has landed since (panel.js currentTeaserBoard).
  function currentTeaserBoard() {
    const loadedAt = teaser.TEASER_LEAGUE_IDS.map((id) => (scannerStatus && scannerStatus.leagueLoadedAt[id]) || 0).join(",");
    if (!teaserBoard || loadedAt !== teaserBoardLoadedAt) {
      teaserBoard = teaser.teaserBoardOf(feedState);
      teaserBoardLoadedAt = loadedAt;
    }
    return teaserBoard;
  }

  // The open BFA teasers, each leg joined to its board game and priced
  // (panel.js openTeasersNow): the Edges rows are sized against them too.
  function openTeasersNow() {
    if (!feedState) return [];
    return teaser.openTeasers(held.records, boardLines(), { now: now(), ladderOf: currentTeaserBoard().ladderOf });
  }

  // One GET of /bets.json, applied exactly as the panel applies it. A failure
  // keeps the held bets and records since when the service has been unreachable.
  async function pollBets() {
    if (betsBusy) return;
    betsBusy = true;
    const at = now();
    try {
      const applied = edgeRows.applyBetsPayload(held, await getServiceJson(serviceFetch, `${deps.betsServiceUrl}/bets.json`), at);
      held = { records: applied.records, crosswalk: applied.crosswalk, pins: applied.pins, fillFairs: applied.fillFairs };
      fillFairIndex = fillfair.fairsByBetId(held.fillFairs);
      if (betsStatus.error) logInfo("bets service reachable again");
      Object.assign(betsStatus, { okAt: at, error: null, unreachableSince: null, generatedAt: applied.generatedAt, sources: applied.sources });
    } catch (error) {
      if (betsStatus.error !== error.message) logWarning(`bets poll failed: ${error.message}`);
      Object.assign(betsStatus, { error: error.message, unreachableSince: betsStatus.unreachableSince ?? at });
    } finally {
      betsBusy = false;
    }
  }

  // The scanner follows the settings' leagues (plus NFL and CFB); a change
  // restarts it, as the panel's league checkboxes do. Not awaited here: a
  // full load takes a while and must not hold the settings poll
  // (scanLoaded() waits for it).
  function followLeagues() {
    const leagues = scannerLeaguesOf(settings.edgeSettings.leagues);
    const signature = leagues.join(",");
    if (signature === scannedLeagues) return;
    scannedLeagues = signature;
    logInfo(`scanning ${leagues.length} league(s)`);
    scanLoad = scanner.start(leagues).catch((error) => logWarning(`scanner start failed: ${error.message}`));
  }

  // One GET of /settings.json. A failure keeps the last settings read (the
  // panel's defaults until one succeeds) and says so in /edges.json.
  async function pollSettings() {
    if (settingsBusy) return;
    settingsBusy = true;
    try {
      const body = await getServiceJson(serviceFetch, `${deps.betsServiceUrl}/settings.json`);
      if (!body || typeof body.settings !== "object" || body.settings === null) throw new Error("settings.json has no settings object");
      if (settings.error) logInfo("settings readable again");
      settings = { ...edgesPayload.settingsFromService(body.settings), source: "service", error: null, okAt: now(), updatedAt: body.updatedAt ?? null };
    } catch (error) {
      if (settings.error !== error.message) logWarning(`settings poll failed (keeping the ${settings.source} settings): ${error.message}`);
      settings = { ...settings, error: error.message };
    } finally {
      settingsBusy = false;
    }
    followLeagues();
  }

  // Bets and settings first, so the first scan uses the right leagues and the
  // first list is sized against what is held; then the polls, and resolve
  // once the first full load has landed.
  async function start() {
    await pollBets();
    await pollSettings();
    pollTimers = [
      timers.setInterval(() => pollBets(), BETS_POLL_MS),
      timers.setInterval(() => pollSettings(), SETTINGS_POLL_MS),
    ];
    await scanLoad;
  }

  function stop() {
    for (const id of pollTimers) timers.clearInterval(id);
    pollTimers = [];
    scanner.stop();
    scannedLeagues = null;
  }

  function edgesPayloadNow() {
    return edgesPayload.buildEdgesPayload({
      feedState, scannerStatus, history, betRecords: held.records, fillFairIndex,
      stakeSettings: settings.stakeSettings, edgeSettings: settings.edgeSettings,
      settingsStatus: { source: settings.source, error: settings.error, okAt: settings.okAt, updatedAt: settings.updatedAt },
      betsStatus: { ...betsStatus }, boardLines: boardLines(), ladderReaderOf,
      teasers: openTeasersNow(), measurement: feedState ? tailFlexMeasurement() : null, now: now(),
    });
  }

  function health() {
    return {
      ok: true, generatedAt: new Date(now()).toISOString(), uptimeSec: Math.round((now() - startedAt) / 1000),
      scanner: scannerStatus ? { phase: scannerStatus.phase, error: scannerStatus.error, lastSnapshotAt: scannerStatus.lastSnapshotAt } : { phase: "starting", error: null },
      betsService: { okAt: betsStatus.okAt, error: betsStatus.error },
      settings: { source: settings.source, okAt: settings.okAt, error: settings.error },
    };
  }

  return { start, stop, pollBets, pollSettings, scanLoaded: () => scanLoad, edgesPayload: edgesPayloadNow, health, scanner };
}

// The phone page's files, by URL path. A fixed map, never a path built from
// the request, so no URL reaches any other file on disk.
const PHONE_DIR = nodePath.join(__dirname, "phone");
const PHONE_FILES = Object.freeze({
  "/": { file: "index.html", contentType: "text/html; charset=utf-8" },
  "/phone.css": { file: "phone.css", contentType: "text/css; charset=utf-8" },
  "/phone_view.js": { file: "phone_view.js", contentType: "text/javascript; charset=utf-8" },
  "/phone.js": { file: "phone.js", contentType: "text/javascript; charset=utf-8" },
});
// Only this origin's own files and its own /edges.json: the page loads no
// fonts, scripts or images from anywhere else, and no other site may frame it.
const PHONE_CSP = "default-src 'none'; script-src 'self'; style-src 'self'; connect-src 'self'; img-src 'self'; base-uri 'none'; form-action 'none'; frame-ancestors 'none'";

//   {"/": {body: Buffer, contentType}, ...} — read once, so a missing file
//   fails the runner at start instead of on the phone.
function loadPhoneFiles(dir = PHONE_DIR) {
  const loaded = {};
  for (const [urlPath, entry] of Object.entries(PHONE_FILES)) {
    loaded[urlPath] = { body: fs.readFileSync(nodePath.join(dir, entry.file)), contentType: entry.contentType };
  }
  return loaded;
}

function sendStatic(response, asset) {
  response.writeHead(200, {
    "Content-Type": asset.contentType,
    "Content-Length": asset.body.length,
    "Cache-Control": "no-store",
    "Content-Security-Policy": PHONE_CSP,
    "X-Content-Type-Options": "nosniff",
    "Referrer-Policy": "no-referrer",
  });
  response.end(asset.body);
}

function sendJson(response, status, body) {
  const text = JSON.stringify(body);
  response.writeHead(status, { "Content-Type": "application/json", "Content-Length": Buffer.byteLength(text), "Cache-Control": "no-store" });
  response.end(text);
}

// The HTTP side: the phone page, GET /edges.json and /health, the Host
// allowlist first on every request. The port comes from the socket, not
// config, so a runner on an ephemeral or non-default port guards itself (as
// the bets service does).
function createHttpServer(runner, { host }, phoneFiles = loadPhoneFiles()) {
  const server = http.createServer((request, response) => {
    const allowed = allowedHosts(host, server.address().port);
    if (!hostAllowed(request.headers.host, allowed)) {
      sendJson(response, 403, { error: `Host must be one of ${JSON.stringify(allowed)}, got ${JSON.stringify(request.headers.host ?? null)}` });
      return;
    }
    if (request.method !== "GET") {
      sendJson(response, 405, { error: `only GET is served, got ${request.method}` });
      return;
    }
    const path = new URL(request.url, "http://runner.invalid").pathname;
    if (Object.hasOwn(phoneFiles, path)) return sendStatic(response, phoneFiles[path]);
    try {
      if (path === "/edges.json") return sendJson(response, 200, runner.edgesPayload());
      if (path === "/health") return sendJson(response, 200, runner.health());
    } catch (error) {
      logError(`${path} failed: ${error.stack || error.message}`);
      return sendJson(response, 500, { error: `${path} failed: ${error.message}` });
    }
    return sendJson(response, 404, { error: `no route for ${path}` });
  });
  return server;
}

async function main() {
  const config = configFromEnv(process.env);
  const runner = createRunner({ fetchImpl: nodeFetch, serviceFetch: fetch, betsServiceUrl: config.betsServiceUrl });
  const server = createHttpServer(runner, config);
  await new Promise((resolve, reject) => {
    server.once("error", reject);
    server.listen(config.port, config.host, resolve);
  });
  log(`serving http://${config.host}:${config.port}/ (phone page) and /edges.json (bets service ${config.betsServiceUrl})`);
  const shutdown = (signal) => {
    log(`stopping (${signal})`);
    runner.stop();
    server.close(() => process.exit(0));
  };
  process.on("SIGINT", () => shutdown("SIGINT"));
  process.on("SIGTERM", () => shutdown("SIGTERM"));
  await runner.start();
}

if (require.main === module) {
  main().catch((error) => {
    logError(`fatal: ${error.stack || error.message}`);
    process.exit(1);
  });
}

module.exports = {
  DEFAULT_HOST, DEFAULT_PORT, DEFAULT_BETS_SERVICE_URL, BETS_POLL_MS, SETTINGS_POLL_MS,
  PHONE_FILES, PHONE_CSP, configFromEnv, allowedHosts, hostAllowed, scannerLeaguesOf, nodeFetch, createRunner,
  createHttpServer, loadPhoneFiles,
};
