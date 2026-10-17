// One-shot CPU profile of the live runner, to name what stalls its event loop.
//
//   docker compose exec -T runner node - < unabated_ticket/server/profile_stalls.js
//   (or: node unabated_ticket/server/profile_stalls.js [seconds])
//
// Finds the runner's node process (cmdline ends in server/runner.js), sends it
// SIGUSR1 (Node opens its inspector on 127.0.0.1:9229), samples its CPU for
// PROFILE_SECONDS, then closes the inspector again. Changes nothing the runner
// does; reads no secrets; writes nothing to disk. Prints every stretch the
// loop was busy longer than LONG_RUN_MS, with the call path that held it.

const fs = require("fs");

const PROFILE_SECONDS = Number(process.argv[2] || process.env.PROFILE_SECONDS || 180);
const SAMPLING_INTERVAL_US = 5000;
const LONG_RUN_MS = 300;
const INSPECTOR_URL = "http://127.0.0.1:9229/json/list";
const IDLE_FRAMES = new Set(["(idle)"]);

// The node process itself, never a wrapper: with compose's `init: true` the
// container's PID 1 is docker-init (tini) whose cmdline also names runner.js,
// and a process without a SIGUSR1 handler dies of it. So the target must be
// exactly one process whose executable is node, running runner.js as its script.
function runnerPid() {
  const pids = fs.readdirSync("/proc").filter((name) => /^\d+$/.test(name) && Number(name) !== process.pid);
  const found = [];
  for (const pid of pids) {
    let args;
    let comm;
    try {
      args = fs.readFileSync(`/proc/${pid}/cmdline`, "utf8").split("\0").filter(Boolean);
      comm = fs.readFileSync(`/proc/${pid}/comm`, "utf8").trim();
    } catch { continue; }
    const isNode = comm === "node" && args.length >= 2 && /(^|\/)node$/.test(args[0]);
    if (isNode && args[1].endsWith("server/runner.js")) found.push({ pid: Number(pid), args });
  }
  if (found.length !== 1) {
    throw new Error(`expected exactly 1 node process running server/runner.js, found ${found.length}: `
      + JSON.stringify(found));
  }
  console.log(`target: pid ${found[0].pid}, command: ${found[0].args.join(" ")}`);
  return found[0].pid;
}

async function inspectorWsUrl() {
  for (let attempt = 0; attempt < 20; attempt += 1) {
    try {
      const targets = await (await fetch(INSPECTOR_URL)).json();
      if (targets.length) return targets[0].webSocketDebuggerUrl;
    } catch { /* not open yet */ }
    await new Promise((resolve) => setTimeout(resolve, 250));
  }
  throw new Error(`the inspector never answered at ${INSPECTOR_URL}`);
}

function connect(url) {
  const socket = new WebSocket(url);
  let nextId = 1;
  const pending = new Map();
  socket.onmessage = (event) => {
    const message = JSON.parse(event.data);
    if (message.id && pending.has(message.id)) {
      const { resolve, reject } = pending.get(message.id);
      pending.delete(message.id);
      if (message.error) reject(new Error(`${message.error.message}`));
      else resolve(message.result);
    }
  };
  const opened = new Promise((resolve, reject) => { socket.onopen = resolve; socket.onerror = reject; });
  const send = (method, params = {}) => new Promise((resolve, reject) => {
    const id = nextId++;
    pending.set(id, { resolve, reject });
    socket.send(JSON.stringify({ id, method, params }));
  });
  return { opened, send, close: () => socket.close() };
}

function frameLabel(node) {
  const frame = node.callFrame;
  const file = frame.url ? frame.url.replace(/^.*\/unabated_ticket\//, "") : "";
  return `${frame.functionName || "(anonymous)"}${file ? ` ${file}:${frame.lineNumber + 1}` : ""}`;
}

function summarize(profile) {
  const byId = new Map(profile.nodes.map((node) => [node.id, node]));
  const parentOf = new Map();
  for (const node of profile.nodes) for (const child of node.children || []) parentOf.set(child, node.id);
  const pathOf = (id) => {
    const path = [];
    for (let at = id; at != null; at = parentOf.get(at)) path.unshift(at);
    return path.slice(1); // drop (root)
  };

  // Consecutive busy samples = one stretch the loop never yielded (close enough at this sampling rate).
  const runs = [];
  let current = null;
  let busyMs = 0;
  profile.samples.forEach((id, index) => {
    const ms = (profile.timeDeltas[index + 1] ?? SAMPLING_INTERVAL_US) / 1000;
    const idle = IDLE_FRAMES.has(byId.get(id).callFrame.functionName);
    if (idle) { current = null; return; }
    busyMs += ms;
    if (!current) { current = { startMs: 0, ms: 0, samples: [] }; runs.push(current); }
    current.ms += ms;
    current.samples.push({ id, ms });
  });

  const totalMs = profile.timeDeltas.reduce((a, b) => a + b, 0) / 1000;
  console.log(`profiled ${(totalMs / 1000).toFixed(0)} s: busy ${(busyMs / 1000).toFixed(1)} s (${(100 * busyMs / totalMs).toFixed(0)}%)`);

  const inclusive = (samples) => {
    const ms = new Map();
    for (const { id, ms: sampleMs } of samples) {
      for (const at of new Set(pathOf(id).map((nodeId) => frameLabel(byId.get(nodeId))))) ms.set(at, (ms.get(at) || 0) + sampleMs);
    }
    return Array.from(ms.entries()).sort((a, b) => b[1] - a[1]);
  };
  const self = (samples) => {
    const ms = new Map();
    for (const { id, ms: sampleMs } of samples) {
      const label = frameLabel(byId.get(id));
      ms.set(label, (ms.get(label) || 0) + sampleMs);
    }
    return Array.from(ms.entries()).sort((a, b) => b[1] - a[1]);
  };

  const long = runs.filter((run) => run.ms >= LONG_RUN_MS).sort((a, b) => b.ms - a.ms);
  console.log(`${long.length} busy stretch(es) of ${LONG_RUN_MS} ms or more; longest ${long[0] ? Math.round(long[0].ms) : 0} ms`);
  for (const run of long.slice(0, 12)) {
    console.log(`\n=== ${Math.round(run.ms)} ms busy`);
    console.log("  inclusive (ms):");
    for (const [label, ms] of inclusive(run.samples).slice(0, 18)) console.log(`    ${String(Math.round(ms)).padStart(6)}  ${label}`);
    console.log("  self (ms):");
    for (const [label, ms] of self(run.samples).slice(0, 6)) console.log(`    ${String(Math.round(ms)).padStart(6)}  ${label}`);
  }

  const allBusy = runs.flatMap((run) => run.samples);
  console.log("\n=== whole window, inclusive (ms):");
  for (const [label, ms] of inclusive(allBusy).slice(0, 30)) console.log(`    ${String(Math.round(ms)).padStart(6)}  ${label}`);
}

async function main() {
  const pid = runnerPid();
  process.kill(pid, "SIGUSR1");
  const client = connect(await inspectorWsUrl());
  await client.opened;
  try {
    await client.send("Profiler.enable");
    await client.send("Profiler.setSamplingInterval", { interval: SAMPLING_INTERVAL_US });
    await client.send("Profiler.start");
    console.log(`runner pid ${pid}: profiling ${PROFILE_SECONDS} s ...`);
    await new Promise((resolve) => setTimeout(resolve, PROFILE_SECONDS * 1000));
    const { profile } = await client.send("Profiler.stop");
    summarize(profile);
  } finally {
    // inspector.close() waits for every session to leave, so schedule it and leave.
    await client.send("Runtime.evaluate", {
      expression: "setTimeout(() => process.mainModule.require('inspector').close(), 1000); 'closing'",
    }).catch((error) => console.log(`could not schedule the inspector close: ${error.message}`));
    client.close();
  }
  await new Promise((resolve) => setTimeout(resolve, 2500));
  try {
    await fetch(INSPECTOR_URL);
    console.log("WARNING: the inspector is still open on 127.0.0.1:9229 (closes on the next runner restart)");
  } catch {
    console.log("inspector closed");
  }
}

main().catch((error) => { console.error(`profile failed: ${error.stack || error.message}`); process.exit(1); });
