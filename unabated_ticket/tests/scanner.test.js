// Run: node --test unabated_ticket/tests
// Drives the scanner loop with an injected fetch that serves the fixtures.
const test = require("node:test");
const assert = require("node:assert/strict");
const fs = require("node:fs");
const path = require("node:path");
const { createScanner, SNAPSHOT_URL, SNAPSHOT_BASE_URL } = require("../extension/scanner.js");
const edgemove = require("../extension/edgemove.js");

const fixture = (name) => fs.readFileSync(path.join(__dirname, "fixtures", name), "utf8");
const KICKOFF_MS = Date.parse("2026-09-13T17:00:00Z");
const NOW = KICKOFF_MS - 3600 * 1000;

function response({ status = 200, body = "", headers = {} }) {
  return {
    ok: status >= 200 && status < 300,
    status,
    headers: { get: (name) => headers[name.toLowerCase()] ?? null },
    json: async () => JSON.parse(body),
    text: async () => body,
  };
}

// A fetch that serves the NFL snapshot (built 20s before NOW), records every
// URL, and lets a test script the next responses.
function fakeFetch(overrides = {}) {
  const calls = [];
  const snapshotBuilt = new Date(NOW - 20 * 1000).toUTCString();
  // Snapshot URLs carry a cache-busting query; match on the base.
  const base = (url) => url.split("?")[0];
  const impl = async (url) => {
    calls.push(base(url));
    if (overrides[base(url)]) return overrides[base(url)]();
    if (base(url) === SNAPSHOT_BASE_URL(1)) return response({ body: fixture("v2_slice.json"), headers: { "last-modified": snapshotBuilt, "content-length": String(overrides.nflBytes || 1000) } });
    return response({ status: 404, body: "not found" });
  };
  impl.calls = calls;
  return impl;
}

const noTimers = { setInterval: () => 1, clearInterval: () => {} };

test("start loads the snapshot; a tick before the league is due fetches nothing, and nothing ever goes to the retired stream", async () => {
  const fetchImpl = fakeFetch();
  const seen = [];
  const scanner = createScanner({ fetchImpl, now: () => NOW, timers: noTimers, onChange: (status) => seen.push(status.phase) });
  await scanner.start([1]);
  const status = scanner.getStatus();
  assert.equal(status.phase, "live");
  assert.equal(status.error, null);
  assert.deepEqual(status.leaguesLoaded, [1]);
  assert.equal(status.lineCount, 62);
  await scanner.tick();
  assert.deepEqual(fetchImpl.calls, [SNAPSHOT_BASE_URL(1)]);
  assert.ok(!fetchImpl.calls.some((url) => url.includes("api-k.unabated.com")));
  assert.deepEqual(seen, ["loading", "live"]); // NFL landed, load complete
});

test("a league that fails to load is reported by name and the rest keep working", async () => {
  const fetchImpl = fakeFetch({ [SNAPSHOT_BASE_URL(2)]: () => response({ status: 503, body: "" }) });
  const scanner = createScanner({ fetchImpl, now: () => NOW, timers: noTimers });
  await scanner.start([1, 2]);
  const status = scanner.getStatus();
  assert.equal(status.phase, "live");
  assert.deepEqual(status.leaguesLoaded, [1]);
  assert.match(status.error, /feed unavailable for CFB \(HTTP 503\)/);
});

test("every league failing is a loud error, not an empty list", async () => {
  const fetchImpl = fakeFetch({ [SNAPSHOT_BASE_URL(1)]: () => { throw new Error("network down"); } });
  const scanner = createScanner({ fetchImpl, now: () => NOW, timers: noTimers });
  await scanner.start([1]);
  const status = scanner.getStatus();
  assert.equal(status.phase, "error");
  assert.match(status.error, /feed unavailable: NFL \(network down\)/);
  assert.equal(status.lineCount, 0);
});

test("a failed re-download keeps the last state and surfaces the message", async () => {
  let clock = NOW;
  let offline = false;
  const serve = fakeFetch();
  const fetchImpl = async (url, options) => {
    if (offline) throw new Error("offline");
    return serve(url, options);
  };
  const scanner = createScanner({ fetchImpl, now: () => clock, timers: noTimers });
  await scanner.start([1]);
  offline = true;
  clock += 61 * 1000; // past the 60 s tier of a small file
  await scanner.tick();
  const status = scanner.getStatus();
  assert.equal(status.phase, "live");
  assert.match(status.error, /feed unavailable: NFL \(offline\)/);
  assert.equal(status.lineCount, 62);
});

test("resume after a long pause resyncs; after a short one it only tops up leagues that are due", async () => {
  let clock = NOW;
  const fetchImpl = fakeFetch();
  const scanner = createScanner({ fetchImpl, now: () => clock, timers: noTimers });
  await scanner.start([1]);
  scanner.pause();
  clock += 30 * 1000;
  await scanner.resume();
  assert.equal(fetchImpl.calls.filter((u) => u === SNAPSHOT_BASE_URL(1)).length, 1);
  scanner.pause();
  clock += 5 * 60 * 1000;
  await scanner.resume();
  assert.equal(fetchImpl.calls.filter((u) => u === SNAPSHOT_BASE_URL(1)).length, 2);
});

test("switching leagues while a snapshot is still downloading discards the old load", async () => {
  let releaseCfb;
  const cfbGate = new Promise((resolve) => { releaseCfb = resolve; });
  const fetchImpl = fakeFetch({
    [SNAPSHOT_BASE_URL(2)]: async () => { await cfbGate; return response({ body: fixture("v2_slice.json").replace(/lg1:/g, "lg2:") }); },
  });
  const scanner = createScanner({ fetchImpl, now: () => NOW, timers: noTimers });
  const first = scanner.start([1, 2]); // CFB hangs
  const second = scanner.start([1]);   // user unticks CFB meanwhile
  await second;
  releaseCfb();
  await first;
  const status = scanner.getStatus();
  assert.deepEqual(status.leagues, [1]);
  assert.deepEqual(status.leaguesLoaded, [1]);
  assert.deepEqual(scanner.getState().leagues, [1]);
});

test("while every league is failing, ticks retry on the 30s throttle instead of every poll", async () => {
  let clock = NOW;
  const fetchImpl = fakeFetch({ [SNAPSHOT_BASE_URL(1)]: () => response({ status: 503, body: "" }) });
  const scanner = createScanner({ fetchImpl, now: () => clock, timers: noTimers });
  await scanner.start([1]);
  await scanner.tick();
  clock += 10 * 1000;
  await scanner.tick();
  assert.equal(fetchImpl.calls.filter((u) => u === SNAPSHOT_BASE_URL(1)).length, 2); // start + first throttled retry
  clock += 30 * 1000;
  await scanner.tick();
  assert.equal(fetchImpl.calls.filter((u) => u === SNAPSHOT_BASE_URL(1)).length, 3);
  assert.equal(scanner.getStatus().phase, "error");
});

test("a full load publishes each league as it lands, CFB last", async () => {
  let releaseCfb;
  const cfbGate = new Promise((resolve) => { releaseCfb = resolve; });
  const fetchImpl = fakeFetch({
    [SNAPSHOT_BASE_URL(2)]: async () => { await cfbGate; return response({ body: fixture("v2_slice.json").replace(/lg1:/g, "lg2:") }); },
  });
  const seen = [];
  const scanner = createScanner({ fetchImpl, now: () => NOW, timers: noTimers, onChange: (status) => seen.push([status.phase, status.leaguesLoaded.slice(), status.loading && status.loading.done]) });
  const done = scanner.start([2, 1]);
  await new Promise((resolve) => setImmediate(resolve));
  await new Promise((resolve) => setImmediate(resolve));
  // NFL is in before CFB has finished downloading.
  assert.deepEqual(scanner.getStatus().leaguesLoaded, [1]);
  assert.deepEqual(scanner.getStatus().staleLeagues, []);
  assert.equal(scanner.getStatus().lineCount, 62);
  assert.deepEqual(scanner.getStatus().loading, { done: 1, total: 2 });
  releaseCfb();
  await done;
  assert.deepEqual(scanner.getStatus().leaguesLoaded, [1, 2]);
  assert.equal(scanner.getStatus().loading, null);
  assert.equal(scanner.getStatus().phase, "live");
  assert.ok(seen.some(([phase, loaded]) => phase === "loading" && loaded.length === 1));
});

test("each league re-downloads its snapshot on a cadence set by its compressed size", async () => {
  let clock = NOW;
  const fetchImpl = fakeFetch();
  const scanner = createScanner({ fetchImpl, now: () => clock, timers: noTimers });
  await scanner.start([1]);
  assert.equal(scanner.refreshEveryMs(1000), 60 * 1000);
  assert.equal(scanner.refreshEveryMs(3 * 1024 * 1024), 120 * 1000);
  assert.equal(scanner.refreshEveryMs(9.7 * 1024 * 1024), 300 * 1000);
  const loads = () => fetchImpl.calls.filter((u) => u === SNAPSHOT_BASE_URL(1)).length;
  clock += 30 * 1000;
  await scanner.tick();
  assert.equal(loads(), 1); // 1 KB file: not due until 60s
  clock += 31 * 1000;
  await scanner.tick();
  assert.equal(loads(), 2);
  assert.equal(scanner.getStatus().phase, "live");
  assert.equal(scanner.getStatus().loading, null); // a refresh shows no progress counter
});

test("a large snapshot refreshes on the slow tier", async () => {
  let clock = NOW;
  const fetchImpl = fakeFetch({ nflBytes: 9.7 * 1024 * 1024 });
  const scanner = createScanner({ fetchImpl, now: () => clock, timers: noTimers });
  await scanner.start([1]);
  clock += 200 * 1000;
  await scanner.tick();
  assert.equal(fetchImpl.calls.filter((u) => u === SNAPSHOT_BASE_URL(1)).length, 1);
  clock += 101 * 1000;
  await scanner.tick();
  assert.equal(fetchImpl.calls.filter((u) => u === SNAPSHOT_BASE_URL(1)).length, 2);
});

test("snapshot URLs carry a query that changes every 30s so CloudFront cannot serve a stale edge copy", () => {
  const t = Date.UTC(2026, 8, 10, 22, 42, 0);
  assert.equal(SNAPSHOT_URL(5, t), `${SNAPSHOT_BASE_URL(5)}?t=${Math.floor(t / 30000)}`);
  assert.equal(SNAPSHOT_URL(5, t + 29 * 1000), SNAPSHOT_URL(5, t));
  assert.notEqual(SNAPSHOT_URL(5, t + 30 * 1000), SNAPSHOT_URL(5, t));
});

test("a full resync keeps every league listed while the others re-download", async () => {
  let releaseCfb = null;
  const fetchImpl = fakeFetch({
    [SNAPSHOT_BASE_URL(2)]: async () => {
      if (releaseCfb) await new Promise((resolve) => { releaseCfb = resolve; });
      // Same slice re-labelled as CFB with distinct market ids so its keys do not collide with NFL's.
      return response({ body: fixture("v2_slice.json").replace(/lg1:/g, "lg2:").replace(/"marketId":(\d+)/g, '"marketId":9$1') });
    },
  });
  const scanner = createScanner({ fetchImpl, now: () => NOW, timers: noTimers });
  await scanner.start([1, 2]);
  assert.deepEqual(scanner.getStatus().leaguesLoaded, [1, 2]);
  assert.equal(scanner.getStatus().lineCount, 124);
  releaseCfb = () => {};
  const resync = scanner.resync();
  await new Promise((resolve) => setImmediate(resolve));
  await new Promise((resolve) => setImmediate(resolve));
  // NFL has re-landed, CFB is still downloading: both are still in the state.
  assert.deepEqual(scanner.getStatus().leaguesLoaded, [1, 2]);
  assert.equal(scanner.getStatus().lineCount, 124);
  releaseCfb();
  await resync;
  assert.deepEqual(scanner.getStatus().leaguesLoaded, [1, 2]);
});

test("a league whose snapshot build is older than 15 min is reported stale", async () => {
  const fetchImpl = fakeFetch();
  const scanner = createScanner({ fetchImpl, now: () => NOW + 60 * 60 * 1000, timers: noTimers }); // snapshot built 1h before "now"
  await scanner.start([1]);
  assert.deepEqual(scanner.getStatus().staleLeagues, [1]);
});

test("status counts alt lines apart from main lines", async () => {
  const scanner = createScanner({ fetchImpl: fakeFetch(), now: () => NOW, timers: noTimers });
  await scanner.start([1]);
  const status = scanner.getStatus();
  assert.equal(status.lineCount, 62);
  assert.equal(status.altLineCount, 27);
});

// ---- per-line history for the edge-move tag (#132) ------------------------

// The v2 fixture with one edit applied to its parsed form.
function snapshotWith(edit) {
  const json = JSON.parse(fixture("v2_slice.json"));
  edit(json);
  return JSON.stringify(json);
}

const NOVIG_ML_KEY = "289357353:ms89:si0:tid6";
const KALSHI_ALT_KEY = "289357360:ms105:si0:tid6:alt-3.5";

test("history: every snapshot line gets one observation on start, and start clears it", async () => {
  let clock = NOW;
  const scanner = createScanner({ fetchImpl: fakeFetch(), now: () => clock, timers: noTimers });
  await scanner.start([1]);
  const history = scanner.getHistory();
  assert.equal(Object.keys(history).length, 62 + 27);
  assert.ok(Object.values(history).every((entries) => entries.length === 1 && entries[0].at === NOW));
  assert.equal(edgemove.edgeMove(history[NOVIG_ML_KEY], NOW).kind, "none");
  assert.equal(history[NOVIG_ML_KEY][0].bacr, -156);
  clock += 60 * 1000;
  await scanner.start([1]);
  assert.ok(Object.values(scanner.getHistory()).every((entries) => entries.length === 1 && entries[0].at === clock));
});

test("history: a re-downloaded snapshot records fair and alt-rung moves as snapshot observations and forgets pulled lines", async () => {
  let clock = NOW;
  let snapshots = 0;
  const fetchImpl = fakeFetch({
    [SNAPSHOT_BASE_URL(1)]: () => {
      snapshots += 1;
      const built = new Date(clock - 20 * 1000).toUTCString();
      if (snapshots === 1) return response({ body: fixture("v2_slice.json"), headers: { "last-modified": built, "content-length": "1000" } });
      return response({
        body: snapshotWith((json) => {
          const rows = json.odds["lg1:pt1:pregame"];
          // Novig moneyline: fair -156 -> -166 (60.9% -> 62.4%), price unchanged.
          rows.find((row) => row.key === "pt1:pregame:bt1:e125807").sides["si0:tid6"].ms89.bacr = -166;
          // Kalshi -3.5 rung: 46.7c -> 44.4c (+113 -> +125), fair unchanged. An
          // exchange line is measured on its exact source price, so both move.
          const kalshiSpread = rows.find((row) => row.key === "pt1:pregame:bt2:e125807").sides["si0:tid6"].ms105;
          const rung = kalshiSpread.alternateLines.find((alt) => alt && alt.points === -3.5);
          rung.americanPrice = 125;
          rung.price = 125;
          rung.sourcePrice = 0.444;
          // The first-half rows are gone.
          delete json.odds["lg1:pt2:pregame"];
        }),
        headers: { "last-modified": built, "content-length": "1000" },
      });
    },
  });
  const scanner = createScanner({ fetchImpl, now: () => clock, timers: noTimers });
  await scanner.start([1]);
  const firstHalfKeys = Object.values(scanner.getState().lines).filter((line) => line.periodTypeId === 2).map((line) => line.key);
  assert.equal(firstHalfKeys.length, 12);
  clock += 61 * 1000; // past the 60s tier
  await scanner.tick();
  assert.equal(snapshots, 2);
  const history = scanner.getHistory();
  const fair = edgemove.edgeMove(history[NOVIG_ML_KEY], clock);
  assert.equal(fair.kind, "fair_to_you");
  assert.equal(fair.sinceMs, 0);
  assert.equal(fair.from.bacr, -156);
  assert.equal(fair.to.bacr, -166);
  const alt = edgemove.edgeMove(history[KALSHI_ALT_KEY], clock);
  assert.equal(alt.kind, "book_away");
  assert.equal(alt.from.price, 113);
  assert.equal(alt.to.price, 125);
  // Every 1H line is gone from the state and from the history.
  const lines = scanner.getState().lines;
  assert.ok(firstHalfKeys.every((key) => lines[key] === undefined && history[key] === undefined));
});

// ---- observingSince: since when the history has no gap (fillfair.js) ------------------

test("observingSince: set when the first load lands, cleared by pause, set again only once the catch-up after resume lands", async () => {
  let clock = NOW;
  const serve = fakeFetch();
  // Every request takes 5 s on the clock, so "when resume() ran" and "when its catch-up landed" differ.
  const fetchImpl = async (url, options) => {
    clock += 5000;
    return serve(url, options);
  };
  const scanner = createScanner({ fetchImpl, now: () => clock, timers: noTimers });
  const loading = scanner.start([1]);
  assert.equal(scanner.getStatus().observingSince, null);
  await loading;
  assert.equal(scanner.getStatus().observingSince, NOW + 5000);
  clock += 61 * 1000;
  await scanner.tick(); // a due re-download inside the run keeps it
  assert.equal(scanner.getStatus().observingSince, NOW + 5000);
  scanner.pause();
  assert.equal(scanner.getStatus().observingSince, null);
  clock += 30 * 1000;
  await scanner.resume(); // short pause, nothing due: no request, so no run yet
  assert.equal(scanner.getStatus().observingSince, null);
  clock += 61 * 1000;
  const dueAt = clock;
  await scanner.tick(); // the next due re-download starts the run
  assert.equal(scanner.getStatus().observingSince, dueAt + 5000);
  scanner.pause();
  clock += 5 * 60 * 1000;
  const longResumeAt = clock;
  await scanner.resume(); // long pause: a full resync
  assert.equal(scanner.getStatus().observingSince, longResumeAt + 5000);
  await scanner.start([1]);
  assert.equal(scanner.getStatus().observingSince, clock); // start() began a new history, and a new run
});

test("observingSince: a pass that lands while the panel is hidden does not start a run", async () => {
  let clock = NOW;
  let releaseSnapshot = null;
  let hold = false;
  const serve = fakeFetch();
  const fetchImpl = async (url, options) => {
    if (hold) await new Promise((resolve) => { releaseSnapshot = resolve; });
    return serve(url, options);
  };
  const scanner = createScanner({ fetchImpl, now: () => clock, timers: noTimers });
  await scanner.start([1]);
  hold = true;
  clock += 61 * 1000;
  const inFlight = scanner.tick();
  await new Promise(setImmediate); // the re-download is now waiting on the network
  scanner.pause();
  clock += 10 * 1000;
  releaseSnapshot();
  await inFlight;
  assert.equal(scanner.getStatus().lastSnapshotAt, clock);
  assert.equal(scanner.getStatus().observingSince, null);
});

test("leagueObservingSince: a league whose loads keep failing while another lands starts its own run over", async () => {
  let clock = NOW;
  let cfbDown = false;
  const fetchImpl = fakeFetch({
    [SNAPSHOT_BASE_URL(2)]: () => {
      if (cfbDown) throw new Error("timed out");
      return response({ body: fixture("v2_venue_ids_slice.json"), headers: { "last-modified": new Date(clock - 20 * 1000).toUTCString(), "content-length": "1000" } });
    },
  });
  const scanner = createScanner({ fetchImpl, now: () => clock, timers: noTimers });
  await scanner.start([1, 2]);
  assert.deepEqual(scanner.getStatus().leagueObservingSince, { 1: NOW, 2: NOW });
  // One missed refresh and a quick recovery: the same run (60 s tier + the 120 s allowance).
  cfbDown = true;
  clock += 61 * 1000;
  await scanner.tick();
  cfbDown = false;
  clock += 31 * 1000;
  await scanner.tick();
  assert.deepEqual(scanner.getStatus().leagueObservingSince, { 1: NOW, 2: NOW });
  // Down for five minutes while NFL keeps landing: CFB's run starts over, NFL's does not.
  cfbDown = true;
  for (let step = 0; step < 10; step += 1) {
    clock += 31 * 1000;
    await scanner.tick();
  }
  cfbDown = false;
  clock += 31 * 1000;
  await scanner.tick();
  const status = scanner.getStatus();
  assert.equal(status.leagueObservingSince[1], NOW);
  assert.equal(status.leagueObservingSince[2], clock);
  assert.equal(status.observingSince, NOW);
  await scanner.start([1]);
  assert.deepEqual(scanner.getStatus().leagueObservingSince, { 1: clock });
});

test("leagueLoadedAt: when each league's last snapshot landed; a failed load keeps it; start() clears it", async () => {
  let clock = NOW;
  let cfbDown = false;
  const fetchImpl = fakeFetch({
    [SNAPSHOT_BASE_URL(2)]: () => {
      if (cfbDown) throw new Error("timed out");
      return response({ body: fixture("v2_venue_ids_slice.json"), headers: { "content-length": "1000" } });
    },
  });
  const scanner = createScanner({ fetchImpl, now: () => clock, timers: noTimers });
  await scanner.start([1, 2]);
  assert.deepEqual(scanner.getStatus().leagueLoadedAt, { 1: NOW, 2: NOW });
  cfbDown = true;
  clock += 61 * 1000;
  await scanner.tick();
  assert.deepEqual(scanner.getStatus().leagueLoadedAt, { 1: clock, 2: NOW });
  cfbDown = false;
  await scanner.start([1]);
  assert.deepEqual(scanner.getStatus().leagueLoadedAt, { 1: clock });
});

test("observingSince: more than two minutes without a successful observation starts a new run", async () => {
  let clock = NOW;
  let offline = false;
  const serve = fakeFetch();
  const fetchImpl = async (url, options) => {
    if (offline) throw new Error("offline");
    return serve(url, options);
  };
  const scanner = createScanner({ fetchImpl, now: () => clock, timers: noTimers });
  await scanner.start([1]);
  assert.equal(scanner.getStatus().observingSince, NOW);
  // 110 s without a success, then back: the same run.
  offline = true;
  clock += 50 * 1000;
  await scanner.tick();
  clock += 50 * 1000;
  await scanner.tick();
  offline = false;
  clock += 10 * 1000;
  await scanner.tick();
  assert.equal(scanner.getStatus().observingSince, NOW);
  const lastSuccessAt = clock;
  // 140 s without one: the history may have missed a move, so the run starts over.
  offline = true;
  clock += 70 * 1000;
  await scanner.tick();
  clock += 60 * 1000;
  await scanner.tick();
  offline = false;
  clock += 10 * 1000;
  await scanner.tick();
  assert.equal(clock - lastSuccessAt, 140 * 1000);
  assert.equal(scanner.getStatus().observingSince, clock);
});
