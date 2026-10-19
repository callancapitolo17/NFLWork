// One-shot CPU profile of the live runner, to name what stalls its event loop.
//
//   docker compose exec -T runner node - < unabated_ticket/server/profile_stalls.js
//   (or: node unabated_ticket/server/profile_stalls.js [seconds])
//
// Finds the runner's node process (cmdline ends in server/runner.js), sends it
// SIGUSR1 (Node opens its inspector on 127.0.0.1:9229), samples its CPU for
// PROFILE_SECONDS (its main thread, and its board thread when it has one),
// then closes the inspector again. Changes nothing the runner
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

async function inspectorTargets() {
  for (let attempt = 0; attempt < 20; attempt += 1) {
    try {
      const targets = await (await fetch(INSPECTOR_URL)).json();
      if (targets.length) return targets;
    } catch { /* not open yet */ }
    await new Promise((resolve) => setTimeout(resolve, 250));
  }
  throw new Error(`the inspector never answered at ${INSPECTOR_URL}`);
}

// A protocol session: send(method, params) resolves to the result; `receive`
// takes each raw message the transport delivers; events go to onEvent.
function protocolSession(transmit) {
  let nextId = 1;
  const pending = new Map();
  const listeners = [];
  const receive = (message) => {
    if (message.id && pending.has(message.id)) {
      const { resolve, reject } = pending.get(message.id);
      pending.delete(message.id);
      if (message.error) reject(new Error(`${message.error.message}`));
      else resolve(message.result);
      return;
    }
    if (message.method) for (const listener of listeners) listener(message);
  };
  const send = (method, params = {}) => new Promise((resolve, reject) => {
    const id = nextId++;
    pending.set(id, { resolve, reject });
    transmit(JSON.stringify({ id, method, params }));
  });
  return { send, receive, onEvent: (listener) => listeners.push(listener) };
}

function connect(url) {
  const socket = new WebSocket(url);
  const session = protocolSession((text) => socket.send(text));
  socket.onmessage = (event) => session.receive(JSON.parse(event.data));
  const opened = new Promise((resolve, reject) => { socket.onopen = resolve; socket.onerror = reject; });
  return { opened, send: session.send, onEvent: session.onEvent, close: () => socket.close() };
}

// The runner's board thread (server/board_thread.js) as its own session,
// reached through the main thread's NodeWorker domain. Null when the runner
// has no worker (a runner from before the board thread).
async function boardThreadSession(client) {
  const attached = new Promise((resolve) => {
    client.onEvent((message) => { if (message.method === "NodeWorker.attachedToWorker") resolve(message.params.sessionId); });
  });
  await client.send("NodeWorker.enable", { waitForDebuggerOnStart: false });
  const sessionId = await Promise.race([attached, new Promise((resolve) => setTimeout(() => resolve(null), 3000))]);
  if (!sessionId) return null;
  const session = protocolSession((message) => client.send("NodeWorker.sendMessageToWorker", { sessionId, message }));
  client.onEvent((message) => {
    if (message.method === "NodeWorker.receivedMessageFromWorker" && message.params.sessionId === sessionId) {
      session.receive(JSON.parse(message.params.message));
    }
  });
  return session;
}

async function startProfiler(session) {
  await session.send("Profiler.enable");
  await session.send("Profiler.setSamplingInterval", { interval: SAMPLING_INTERVAL_US });
  await session.send("Profiler.start");
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

async function inspectorAnswers() {
  try {
    await fetch(INSPECTOR_URL);
    return true;
  } catch {
    return false;
  }
}

// Whatever happens after SIGUSR1, the runner's inspector is closed (or a
// warning says it is still open): an inspector left open is code execution
// for anything on the VM's loopback until the runner restarts.
async function main() {
  if (await inspectorAnswers()) throw new Error(`something already answers at ${INSPECTOR_URL}: refusing to signal the runner`);
  const pid = runnerPid();
  process.kill(pid, "SIGUSR1");
  let client = null;
  try {
    const targets = await inspectorTargets();
    if (!targets[0].url.endsWith("server/runner.js")) throw new Error(`the inspector at ${INSPECTOR_URL} is not the runner's (${targets[0].url})`);
    client = connect(targets[0].webSocketDebuggerUrl);
    await client.opened;
    const board = await boardThreadSession(client);
    await startProfiler(client);
    if (board) await startProfiler(board);
    console.log(`runner pid ${pid}: profiling ${PROFILE_SECONDS} s (main thread${board ? " and board thread" : " only: no board thread"}) ...`);
    await new Promise((resolve) => setTimeout(resolve, PROFILE_SECONDS * 1000));
    const { profile } = await client.send("Profiler.stop");
    console.log("\n##### main thread (HTTP: every reply waits for its stalls)");
    summarize(profile);
    if (board) {
      const boardProfile = (await board.send("Profiler.stop")).profile;
      console.log("\n##### board thread (scanner, pages, quotes, taker)");
      summarize(boardProfile);
    }
  } finally {
    if (client) {
      // inspector.close() waits for every session to leave, so schedule it and leave.
      await client.send("Runtime.evaluate", {
        expression: "setTimeout(() => { try { process.mainModule.require('inspector').close(); } catch {} }, 1000); 'closing'",
      }).catch((error) => console.log(`could not schedule the inspector close: ${error.message}`));
      client.close();
    }
    await new Promise((resolve) => setTimeout(resolve, 2500));
    if (await inspectorAnswers()) console.log("WARNING: the inspector is still open on 127.0.0.1:9229 (closes on the next runner restart)");
    else console.log("inspector closed");
  }
}

main().catch((error) => { console.error(`profile failed: ${error.stack || error.message}`); process.exit(1); });
