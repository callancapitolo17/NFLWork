// Run: node --test unabated_ticket/tests
// Drives the scanner loop with an injected fetch that serves the fixtures.
const test = require("node:test");
const assert = require("node:assert/strict");
const fs = require("node:fs");
const path = require("node:path");
const { createScanner, SNAPSHOT_URL, CHANGES_URL } = require("../extension/scanner.js");

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

// A fetch that serves the NFL snapshot (built 20s before NOW) and the changes
// fixture, records every URL, and lets a test script the next responses.
function fakeFetch(overrides = {}) {
  const calls = [];
  const snapshotBuilt = new Date(NOW - 20 * 1000).toUTCString();
  const impl = async (url) => {
    calls.push(url);
    if (overrides[url]) return overrides[url]();
    if (url === SNAPSHOT_URL(1)) return response({ body: fixture("v2_slice.json"), headers: { "last-modified": snapshotBuilt } });
    if (url.startsWith(CHANGES_URL)) return response({ body: fixture("changes_slice.json") });
    return response({ status: 404, body: "not found" });
  };
  impl.calls = calls;
  return impl;
}

const noTimers = { setInterval: () => 1, clearInterval: () => {} };

test("start loads the snapshot, seeds the cursor from Last-Modified, then polls", async () => {
  const fetchImpl = fakeFetch();
  const seen = [];
  const scanner = createScanner({ fetchImpl, now: () => NOW, timers: noTimers, onChange: (status) => seen.push(status.phase) });
  await scanner.start([1]);
  let status = scanner.getStatus();
  assert.equal(status.phase, "live");
  assert.equal(status.error, null);
  assert.deepEqual(status.leaguesLoaded, [1]);
  assert.equal(status.lineCount, 62);
  // Cursor = whole seconds of (NOW - 20s) since 2021-01-06, in ns.
  const expectedCursor = `${Math.floor((NOW - 20 * 1000 - Date.UTC(2021, 0, 6)) / 1000)}000000000`;
  assert.equal(status.cursor, expectedCursor);

  await scanner.tick();
  status = scanner.getStatus();
  assert.equal(fetchImpl.calls[1], `${CHANGES_URL}/${expectedCursor}`);
  assert.equal(status.pollCount, 1);
  assert.equal(status.lastPollLines, 14);
  assert.equal(status.cursor, "179164041314243100");
  assert.equal(status.lineCount, 66);
  assert.equal(scanner.getState().lines["289357360:ms4:si0:tid6"].points, -3.5);
  assert.deepEqual(seen, ["live", "live"]);
});

test("a stale snapshot leaves the cursor null so the first poll uses the server default", async () => {
  const fetchImpl = fakeFetch();
  const scanner = createScanner({ fetchImpl, now: () => NOW + 10 * 60 * 1000, timers: noTimers });
  await scanner.start([1]);
  assert.equal(scanner.getStatus().cursor, null);
  await scanner.tick();
  assert.equal(fetchImpl.calls[1], CHANGES_URL);
});

test("a rejected cursor triggers a resync from snapshots", async () => {
  let failNext = true;
  const fetchImpl = fakeFetch({
    [`${CHANGES_URL}/179164041314243100`]: () => {
      if (!failNext) return response({ body: fixture("changes_slice.json") });
      failNext = false;
      return response({ body: JSON.stringify({ latestTimestamp: 1, resultCode: "Failed", results: [] }) });
    },
  });
  const scanner = createScanner({ fetchImpl, now: () => NOW, timers: noTimers });
  await scanner.start([1]);
  await scanner.tick(); // moves the cursor to the fixture's latestTimestamp
  await scanner.tick(); // Failed -> resync
  const snapshotLoads = fetchImpl.calls.filter((url) => url === SNAPSHOT_URL(1)).length;
  assert.equal(snapshotLoads, 2);
  assert.equal(scanner.getStatus().phase, "live");
});

test("a league that fails to load is reported by name and the rest keep working", async () => {
  const fetchImpl = fakeFetch({ [SNAPSHOT_URL(2)]: () => response({ status: 503, body: "" }) });
  const scanner = createScanner({ fetchImpl, now: () => NOW, timers: noTimers });
  await scanner.start([1, 2]);
  const status = scanner.getStatus();
  assert.equal(status.phase, "live");
  assert.deepEqual(status.leaguesLoaded, [1]);
  assert.match(status.error, /feed unavailable for CFB \(HTTP 503\)/);
  await scanner.tick();
  assert.equal(scanner.getStatus().pollCount, 1);
});

test("every league failing is a loud error, not an empty list", async () => {
  const fetchImpl = fakeFetch({ [SNAPSHOT_URL(1)]: () => { throw new Error("network down"); } });
  const scanner = createScanner({ fetchImpl, now: () => NOW, timers: noTimers });
  await scanner.start([1]);
  const status = scanner.getStatus();
  assert.equal(status.phase, "error");
  assert.match(status.error, /feed unavailable: NFL \(network down\)/);
  assert.equal(status.lineCount, 0);
});

test("a network error mid-stream keeps the last state and surfaces the message", async () => {
  const fetchImpl = fakeFetch({
    [CHANGES_URL]: () => { throw new Error("offline"); },
  });
  const scanner = createScanner({ fetchImpl, now: () => NOW + 10 * 60 * 1000, timers: noTimers });
  await scanner.start([1]);
  await scanner.tick();
  const status = scanner.getStatus();
  assert.equal(status.phase, "live");
  assert.match(status.error, /feed unavailable: offline/);
  assert.equal(status.lineCount, 62);
});

test("resume after a long pause resyncs; after a short one it just polls", async () => {
  let clock = NOW;
  const fetchImpl = fakeFetch();
  const scanner = createScanner({ fetchImpl, now: () => clock, timers: noTimers });
  await scanner.start([1]);
  scanner.pause();
  clock += 30 * 1000;
  await scanner.resume();
  assert.equal(fetchImpl.calls.filter((u) => u === SNAPSHOT_URL(1)).length, 1);
  scanner.pause();
  clock += 5 * 60 * 1000;
  await scanner.resume();
  assert.equal(fetchImpl.calls.filter((u) => u === SNAPSHOT_URL(1)).length, 2);
});
