"use strict";
// The runner's board on its own thread (Cal, 2026-10-19: a slow scanner must
// not freeze the phone page or Place).
//
// The board thread (a worker_threads Worker running this file) holds the whole
// runner as createRunner builds it: the scanner and the board, the bets and
// settings polls, learning, closing fairs, the taker's passes and its order
// sends, and every page and quote built off the board. The main thread
// (runner.js main) holds only the HTTP server and this file's board client,
// which answers:
//   /edges.json, /teasers.json, /scenarios.json, /taker.json, /tau.json, /live_prices.json  from the last
//       copy the board thread built, at once, when that copy is under
//       PAYLOAD_SERVE_MAX_AGE_MS old, and asks for a new copy for the next
//       request (so a page polling every 15 s reads a list at most one poll
//       old); an older or missing copy, or ?fresh=1, waits for a build
//       (fresh asks that land while one fresh build runs share the next one).
//   /health  from the health the board thread posts every HEALTH_POST_MS, plus
//       the board thread's own state (running, restarting, restarts); ok is
//       false while it restarts or has posted nothing for HEALTH_STALE_MS.
//   /attach.json, /ticket_line.json (the open Ticket's 2 s refresh) and the Place quotes (/kalshi_quote.json, /novig_quote.json,
//   /polymarket_quote.json, /prophetx_quote.json, /bfa_quote.json)  one message to the
//       board thread, answered when it is free (after the work in progress).
// If the board thread dies, every waiting request answers 503, the copies are
// dropped and a new board thread starts after BOARD_RESTART_DELAY_MS (doubling
// while threads keep dying): a fresh
// board, as a process restart gives (the taker fails closed until it is live).
// Side effects: none of its own; the board thread's are runner.js's.

const { Worker, isMainThread, parentPort, workerData } = require("node:worker_threads");

// A copy this young is served without waiting; the phone polls every 15 s, the tracker every 60 s.
const PAYLOAD_SERVE_MAX_AGE_MS = 60 * 1000;
const HEALTH_POST_MS = 5 * 1000;
// No health from the board thread for this long: it is stuck, and /health says ok: false.
const HEALTH_STALE_MS = 3 * HEALTH_POST_MS;
// A dead board thread is replaced after BOARD_RESTART_DELAY_MS, doubling for
// each death in a row up to BOARD_RESTART_MAX_DELAY_MS (a thread that dies on
// every boot must not re-download every league every few seconds); a thread
// that lived BOARD_HEALTHY_RUN_MS starts the count again.
const BOARD_RESTART_DELAY_MS = 5 * 1000;
const BOARD_RESTART_MAX_DELAY_MS = 5 * 60 * 1000;
const BOARD_HEALTHY_RUN_MS = 10 * 60 * 1000;
// The pages the main thread serves from a copy, by the runner method that builds each.
const PAYLOAD_METHODS = { edges: "edgesPayload", teasers: "teasersPayload", scenarios: "scenariosPayload", taker: "takerPayload", tau: "tauPayload", livePrices: "livePricesPayload" };
const HTTP_SERVICE_UNAVAILABLE = 503;

// A reply body already serialized: runner.js sendJson writes its text as is.
class RawJson {
  constructor(text) {
    this.text = text;
  }
}

function httpError(status, message) {
  const error = new Error(message);
  error.status = status;
  return error;
}

// The main thread's view of the board thread: the runner methods the HTTP
// server calls, each answered off a copy or by one message.
//   options  {workerPath (default this file), workerData (runner.js config),
//             now, logInfo, logWarning, createWorker (path, options) -> Worker,
//             setTimeout}
function createBoardClient(options) {
  const now = options.now || (() => Date.now());
  const logInfo = options.logInfo || console.log;
  const logWarning = options.logWarning || console.error;
  const createWorker = options.createWorker || ((path, workerOptions) => new Worker(path, workerOptions));
  const later = options.setTimeout || setTimeout;
  const workerPath = options.workerPath || __filename;

  let worker = null;
  let stopped = false;
  let restarts = 0;
  // Deaths in a row (a thread that ran BOARD_HEALTHY_RUN_MS resets it).
  let deathsInARow = 0;
  let startedAt = null;
  let lastExit = null;
  let nextCallId = 1;
  // callId -> {resolve, reject}
  const pending = new Map();
  // payload name -> {text, builtAt}
  const copies = new Map();
  // payload name -> the build in flight (a Promise of its RawJson)
  const builds = new Map();
  // payload name -> {running: the fresh build in flight, queued: {rereadKeys, rereading, promise} | null}
  const freshBuilds = new Map();
  let workerHealth = null;
  let healthPostedAt = null;

  function start() {
    stopped = false;
    spawn();
  }

  function spawn() {
    startedAt = now();
    const spawned = createWorker(workerPath, { workerData: options.workerData });
    worker = spawned;
    spawned.on("message", (message) => onMessage(message));
    spawned.on("error", (error) => logWarning(`board thread failed: ${error.stack || error.message}`));
    spawned.on("exit", (code) => onExit(spawned, code));
  }

  function onMessage(message) {
    if (message.type === "health") {
      workerHealth = message.health;
      healthPostedAt = now();
      return;
    }
    if (message.type !== "reply") return;
    const waiting = pending.get(message.id);
    if (!waiting) return;
    pending.delete(message.id);
    if (message.ok) waiting.resolve(message.json);
    else waiting.reject(httpError(message.status || 500, message.error));
  }

  // Every waiting request fails, the copies built on the dead board go, and
  // (unless stopping) a new board thread starts after a pause that grows while
  // threads keep dying.
  function onExit(exited, code) {
    if (exited !== worker) return;
    worker = null;
    lastExit = { at: new Date(now()).toISOString(), code };
    for (const waiting of pending.values()) waiting.reject(httpError(HTTP_SERVICE_UNAVAILABLE, `the board thread stopped (exit code ${code})`));
    pending.clear();
    copies.clear();
    builds.clear();
    freshBuilds.clear();
    workerHealth = null;
    healthPostedAt = null;
    if (stopped) return;
    deathsInARow = now() - startedAt >= BOARD_HEALTHY_RUN_MS ? 1 : deathsInARow + 1;
    const delayMs = Math.min(BOARD_RESTART_DELAY_MS * 2 ** (deathsInARow - 1), BOARD_RESTART_MAX_DELAY_MS);
    logWarning(`board thread exited with code ${code}; starting a new one in ${delayMs / 1000} s`);
    later(() => {
      if (stopped || worker) return;
      restarts += 1;
      logInfo(`board thread restart ${restarts}`);
      spawn();
    }, delayMs);
  }

  // One runner method run on the board thread; resolves to its reply as JSON text.
  function call(method, args = []) {
    if (!worker) return Promise.reject(httpError(HTTP_SERVICE_UNAVAILABLE, "the board thread is restarting"));
    const id = nextCallId;
    nextCallId += 1;
    return new Promise((resolve, reject) => {
      pending.set(id, { resolve, reject });
      worker.postMessage({ type: "call", id, method, args });
    });
  }

  // One build of a page on the board thread, kept as the page's copy.
  function runBuild(name, fresh, reread) {
    const builtAt = now();
    return call("buildPayload", [name, fresh, reread]).then((text) => {
      copies.set(name, { text, builtAt });
      return new RawJson(text);
    });
  }

  // A build that is not fresh; concurrent asks for the same page share one.
  function build(name) {
    if (builds.has(name)) return builds.get(name);
    const shared = runBuild(name, false, null).finally(() => { if (builds.get(name) === shared) builds.delete(name); });
    builds.set(name, shared);
    return shared;
  }

  // A fresh build: one runs per page and one waits behind it. An ask while one
  // runs, or while the waiting build waits on its rereads, joins the waiting
  // build, which starts when the running one ends, so its bets read is sent
  // after the ask. A joining ask's reread starts at once and the waiting build
  // starts only once every such reread has landed.
  function freshBuild(name, reread) {
    let state = freshBuilds.get(name);
    if (!state) {
      state = { running: null, queued: null };
      freshBuilds.set(name, state);
    }
    if (!state.running && !state.queued) {
      state.running = startFreshBuild(name, state, reread);
      return state.running;
    }
    if (!state.queued) {
      const queued = { rereadKeys: new Set(), rereading: [], promise: null };
      queued.promise = state.running.catch(() => {}).then(() => allRereadsLanded(queued.rereading)).then(() => {
        state.queued = null;
        state.running = startFreshBuild(name, state, null);
        return state.running;
      });
      state.queued = queued;
    }
    if (reread && !state.queued.rereadKeys.has(reread)) {
      state.queued.rereadKeys.add(reread);
      state.queued.rereading.push(call("rereadLiveLine", [reread]).catch(() => {}));
    }
    return state.queued.promise;
  }

  // Waits for every reread in the list, including any pushed while it waits.
  async function allRereadsLanded(rereading) {
    let waitedFor = -1;
    while (waitedFor !== rereading.length) {
      waitedFor = rereading.length;
      await Promise.all(rereading);
    }
  }

  function startFreshBuild(name, state, reread) {
    const running = runBuild(name, true, reread).finally(() => { if (state.running === running) state.running = null; });
    return running;
  }

  // The page's copy when young enough (and a rebuild started for the next
  // ask), else the build it waits for.
  function payload(name, { fresh = false, reread = null } = {}) {
    if (fresh) return freshBuild(name, reread);
    const copy = copies.get(name);
    if (copy && now() - copy.builtAt <= PAYLOAD_SERVE_MAX_AGE_MS) {
      build(name).catch((error) => logWarning(`board thread: ${name} rebuild failed: ${error.message}`));
      return Promise.resolve(new RawJson(copy.text));
    }
    return build(name);
  }

  function health() {
    const silentMs = healthPostedAt === null ? null : now() - healthPostedAt;
    const boardThread = {
      state: worker ? "running" : (stopped ? "stopped" : "restarting"),
      startedAt: startedAt === null ? null : new Date(startedAt).toISOString(), restarts, lastExit,
      healthPostedAt: healthPostedAt === null ? null : new Date(healthPostedAt).toISOString(),
    };
    // Up while starting (no health posted yet); down while restarting or silent past HEALTH_STALE_MS.
    const ok = Boolean(worker) && (silentMs === null || silentMs <= HEALTH_STALE_MS);
    if (!workerHealth) return { ok, scanner: { phase: "starting", error: null }, boardThread };
    return { ...workerHealth, ok, boardThread };
  }

  function stop() {
    stopped = true;
    if (worker) worker.terminate();
  }

  return {
    start, stop, health,
    // The HTTP server reads the bets (and a ?reread= line) first on ?fresh=1; the fresh build does that on the board thread.
    pollBets: () => Promise.resolve(),
    rereadLiveLine: () => Promise.resolve(),
    edgesPayload: (opts) => payload("edges", opts),
    teasersPayload: (opts) => payload("teasers", opts),
    // The Live tab's filtered cards are built on demand; the unfiltered page is served off its copy.
    scenariosPayload: (params) => (params && [...params.keys()].length
      ? call("scenariosFiltered", [params.toString()]).then((text) => new RawJson(text))
      : payload("scenarios")),
    takerPayload: () => payload("taker"),
    tauPayload: () => payload("tau"),
    livePricesPayload: () => payload("livePrices"),
    attachPayload: (params) => call("attachPayload", [params.toString()]).then((text) => new RawJson(text)),
    ticketLinePayload: (params) => call("ticketLinePayload", [params.toString()]).then((text) => new RawJson(text)),
    kalshiQuote: (params) => call("kalshiQuote", [params.toString()]).then((text) => new RawJson(text)),
    rfqPrice: (params) => call("rfqPrice", [params.toString()]).then((text) => new RawJson(text)),
    novigQuote: (params) => call("novigQuote", [params.toString()]).then((text) => new RawJson(text)),
    polymarketQuote: (params) => call("polymarketQuote", [params.toString()]).then((text) => new RawJson(text)),
    prophetxQuote: (params) => call("prophetxQuote", [params.toString()]).then((text) => new RawJson(text)),
    bfaQuote: (params) => call("bfaQuote", [params.toString()]).then((text) => new RawJson(text)),
    teamKeys: (params) => call("teamKeys", [params.toString()]).then((text) => new RawJson(text)),
  };
}

// The board thread: the whole runner, answering the main thread's calls.
//   port    the parent port
//   runner  createRunner's result (built by runBoardThread)
function serveBoardCalls(port, runner, logWarning) {
  const methods = {
    async buildPayload(name, fresh, reread = null) {
      const method = PAYLOAD_METHODS[name];
      if (!method) throw httpError(404, `no page ${JSON.stringify(name)}`);
      // Side by side, so the reread's wait does not stack on the bets read; it never throws.
      if (fresh) await Promise.all([runner.pollBets(), reread && name === "edges" ? runner.rereadLiveLine(reread) : null]);
      return runner[method]();
    },
    rereadLiveLine: (key) => runner.rereadLiveLine(key),
    attachPayload: (search) => runner.attachPayload(new URLSearchParams(search)),
    ticketLinePayload: (search) => runner.ticketLinePayload(new URLSearchParams(search)),
    scenariosFiltered: (search) => runner.scenariosPayload(new URLSearchParams(search)),
    kalshiQuote: (search) => runner.kalshiQuote(new URLSearchParams(search)),
    rfqPrice: (search) => runner.rfqPrice(new URLSearchParams(search)),
    novigQuote: (search) => runner.novigQuote(new URLSearchParams(search)),
    polymarketQuote: (search) => runner.polymarketQuote(new URLSearchParams(search)),
    prophetxQuote: (search) => runner.prophetxQuote(new URLSearchParams(search)),
    bfaQuote: (search) => runner.bfaQuote(new URLSearchParams(search)),
    teamKeys: (search) => runner.teamKeys(new URLSearchParams(search)),
  };
  port.on("message", async (message) => {
    if (message.type !== "call") return;
    const method = methods[message.method];
    try {
      if (!method) throw httpError(404, `no board call ${JSON.stringify(message.method)}`);
      const result = await method(...message.args);
      port.postMessage({ type: "reply", id: message.id, ok: true, json: JSON.stringify(result) });
    } catch (error) {
      if (!error.status) logWarning(`board call ${message.method} failed: ${error.stack || error.message}`);
      port.postMessage({ type: "reply", id: message.id, ok: false, status: error.status || 500, error: error.message });
    }
  });
}

// The board thread's entry: runner.js's createRunner with its real network
// (as main() built it before the split), its calls served, its health posted.
function runBoardThread(port, config) {
  // Required here, not at the top: runner.js requires this file for the main thread's side.
  const runnerLib = require("./runner.js");
  const runner = runnerLib.createRunner({
    fetchImpl: runnerLib.nodeFetch, serviceFetch: fetch, kalshiFetch: fetch, kalshiBaseUrl: config.kalshiBaseUrl, betsServiceUrl: config.betsServiceUrl,
    novigBaseUrl: config.novigBaseUrl,
  });
  serveBoardCalls(port, runner, runnerLib.logError);
  const postHealth = () => port.postMessage({ type: "health", health: runner.health() });
  postHealth();
  setInterval(postHealth, HEALTH_POST_MS);
  // A start that throws ends the board thread; the main thread starts another.
  runner.start().catch((error) => {
    runnerLib.logError(`board thread: start failed: ${error.stack || error.message}`);
    process.exit(1);
  });
}

module.exports = {
  PAYLOAD_SERVE_MAX_AGE_MS, HEALTH_POST_MS, HEALTH_STALE_MS, BOARD_RESTART_DELAY_MS, BOARD_RESTART_MAX_DELAY_MS, BOARD_HEALTHY_RUN_MS, PAYLOAD_METHODS,
  RawJson, createBoardClient, serveBoardCalls,
};

// Last, so runner.js (required by runBoardThread) finds this file's exports complete.
if (!isMainThread && workerData && workerData.role === "board") runBoardThread(parentPort, workerData.config);
