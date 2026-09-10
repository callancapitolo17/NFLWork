// Unabated Ticket — cross-league edge scanner loop.
//
// Runs inside the side panel page (never the service worker): the panel is
// open whenever you are betting, setInterval is reliable there, and closing
// the panel stops every request. Snapshot per enabled league on start and
// every RESYNC_MS, then the changes stream every POLL_MS from the last cursor.
//
// Side effects: network requests to content.unabated.com and
// api-k.unabated.com only. No storage, no DOM. State lives in memory and is
// handed to the panel through onChange(status, state).
//
// Loaded as a plain <script> in panel.html (globalThis.UnabatedScanner) and via
// require() in tests, where fetch and timers are injected.

(function (root) {
  "use strict";

  const feed = typeof module !== "undefined" && module.exports ? require("./feed.js") : root.UnabatedFeed;

  const SNAPSHOT_URL = (leagueId) => `https://content.unabated.com/markets/v2/league/${leagueId}/odds.json`;
  const CHANGES_URL = "https://api-k.unabated.com/api/markets/changes/query";
  const POLL_MS = 10000;
  const RESYNC_MS = 10 * 60 * 1000;
  // A changes response carries at most this many ~1.3s batches; a full page
  // means there is more to read right away.
  const FULL_PAGE_BATCHES = 7;
  const MAX_PAGES_PER_POLL = 8;
  // The server rejects cursors older than a few minutes (Failed at 300s,
  // fine at 180s on 2026-09-10); a snapshot is regenerated every ~27s.
  const CURSOR_MAX_AGE_MS = 150 * 1000;
  // After a pause longer than this the cursor may be dead: resync instead.
  const PAUSE_RESYNC_MS = 120 * 1000;
  const FAILED_LEAGUE_RETRY_MS = 30 * 1000;
  // A hung fetch would otherwise hold `busy` forever and stall the loop silently.
  const SNAPSHOT_TIMEOUT_MS = 60 * 1000;
  const CHANGES_TIMEOUT_MS = 20 * 1000;

  function createScanner(deps) {
    const fetchImpl = (deps && deps.fetchImpl) || ((...args) => root.fetch(...args));
    const now = (deps && deps.now) || (() => Date.now());
    const timers = (deps && deps.timers) || { setInterval: (f, ms) => root.setInterval(f, ms), clearInterval: (id) => root.clearInterval(id) };
    const onChange = (deps && deps.onChange) || (() => {});

    let state = feed.emptyState();
    let leagues = [];
    let cursor = null;
    let pollTimer = null;
    let resyncTimer = null;
    let busy = false;
    let paused = false;
    const status = {
      phase: "idle", // idle | loading | live | error
      error: null,
      leagues: [],
      leaguesLoaded: [],
      leagueErrors: {},
      lastSnapshotAt: null,
      lastUpdateAt: null,
      lastPollAt: null,
      lastPollLines: 0,
      pollCount: 0,
      lineCount: 0,
      eventCount: 0,
      cursor: null,
    };
    let failedLeagueRetryAt = 0;
    // Bumped by start(); a load that began under an older generation is discarded.
    let generation = 0;

    async function fetchWithTimeout(url, options, timeoutMs) {
      const controller = typeof AbortController === "function" ? new AbortController() : null;
      const timer = controller ? setTimeout(() => controller.abort(), timeoutMs) : null;
      try {
        return await fetchImpl(url, controller ? { ...options, signal: controller.signal } : options);
      } catch (error) {
        throw new Error(controller && controller.signal.aborted ? `timed out after ${timeoutMs / 1000}s` : error.message);
      } finally {
        if (timer) clearTimeout(timer);
      }
    }

    function notify() {
      status.lineCount = feed.countLines(state);
      status.eventCount = Object.keys(state.events).length;
      status.cursor = cursor;
      onChange({ ...status }, state);
    }

    function setError(message) {
      status.error = message;
      status.phase = status.leaguesLoaded.length ? "live" : "error";
    }

    async function fetchSnapshot(leagueId) {
      // no-cache = revalidate with If-None-Match; a 304 serves the cached body.
      const response = await fetchWithTimeout(SNAPSHOT_URL(leagueId), { cache: "no-cache", credentials: "include" }, SNAPSHOT_TIMEOUT_MS);
      if (!response.ok) throw new Error(`HTTP ${response.status}`);
      const json = await response.json();
      const lastModified = response.headers && response.headers.get ? response.headers.get("last-modified") : null;
      const builtAt = lastModified ? Date.parse(lastModified) : NaN;
      return { state: feed.parseSnapshot(json, { leagueId }), builtAt: Number.isFinite(builtAt) ? builtAt : null };
    }

    // Load every requested league; keep whatever succeeds and report the rest
    // by name. The changes cursor starts at the OLDEST snapshot build time when
    // that is recent enough, so nothing between build and first poll is missed.
    async function loadSnapshots(leagueIds) {
      const startedUnder = generation;
      status.phase = status.leaguesLoaded.length ? status.phase : "loading";
      const results = await Promise.all(leagueIds.map((leagueId) => fetchSnapshot(leagueId).then(
        (loaded) => ({ leagueId, loaded }),
        (error) => ({ leagueId, error }),
      )));
      // start() ran meanwhile: these leagues are no longer what the panel wants.
      if (startedUnder !== generation) return false;
      const loadedStates = [];
      let oldestBuild = null;
      const errors = { ...status.leagueErrors };
      for (const result of results) {
        if (result.error) {
          errors[result.leagueId] = result.error.message;
          continue;
        }
        delete errors[result.leagueId];
        loadedStates.push(result.loaded.state);
        if (result.loaded.builtAt != null) oldestBuild = oldestBuild == null ? result.loaded.builtAt : Math.min(oldestBuild, result.loaded.builtAt);
      }
      status.leagueErrors = errors;
      if (loadedStates.length === 0) {
        setError(`feed unavailable: ${describeLeagueErrors(errors)}`);
        notify();
        return false;
      }
      // Keep already-loaded leagues that were not part of this load (partial retry).
      const keep = status.leaguesLoaded.length && leagueIds.length < leagues.length
        ? [stateWithout(state, leagueIds)]
        : [];
      state = feed.mergeStates([...keep, ...loadedStates]);
      status.leaguesLoaded = Array.from(new Set(state.leagues)).sort((a, b) => a - b);
      status.lastSnapshotAt = now();
      status.phase = "live";
      status.error = Object.keys(errors).length ? `feed unavailable for ${describeLeagueErrors(errors)}` : null;
      const recent = oldestBuild != null && now() - oldestBuild <= CURSOR_MAX_AGE_MS;
      cursor = recent ? feed.cursorFromDate(oldestBuild) : null;
      notify();
      return true;
    }

    function stateWithout(current, leagueIds) {
      const drop = new Set(leagueIds);
      const kept = feed.emptyState();
      kept.leagues = current.leagues.filter((id) => !drop.has(id));
      kept.teams = current.teams;
      kept.books = current.books;
      for (const [id, event] of Object.entries(current.events)) if (!drop.has(event.leagueId)) kept.events[id] = event;
      for (const [key, line] of Object.entries(current.lines)) if (!drop.has(line.leagueId)) kept.lines[key] = line;
      return kept;
    }

    function describeLeagueErrors(errors) {
      return Object.entries(errors).map(([id, message]) => `${(feed.LEAGUES[id] || { label: `league ${id}` }).label} (${message})`).join(", ");
    }

    async function fetchChangesPage() {
      const url = cursor ? `${CHANGES_URL}/${cursor}` : CHANGES_URL;
      const response = await fetchWithTimeout(url, { credentials: "include" }, CHANGES_TIMEOUT_MS);
      if (!response.ok) throw new Error(`changes HTTP ${response.status}`);
      return feed.parseChanges(await response.text());
    }

    // Read the stream until a page comes back short. A Failed result means the
    // cursor expired: resync from snapshots and restart from the server default.
    async function pollChanges() {
      let lines = 0;
      const startedUnder = generation;
      for (let page = 0; page < MAX_PAGES_PER_POLL; page += 1) {
        const parsed = await fetchChangesPage();
        if (startedUnder !== generation) return;
        if (!parsed.ok) {
          cursor = null;
          status.error = `changes cursor rejected (${parsed.resultCode}); resyncing`;
          await loadSnapshots(leagues);
          return;
        }
        const counts = feed.applyChanges(state, parsed);
        lines += counts.applied;
        cursor = parsed.cursor || cursor;
        if (parsed.batches < FULL_PAGE_BATCHES) break;
      }
      status.pollCount += 1;
      status.lastPollAt = now();
      status.lastPollLines = lines;
      if (lines > 0) status.lastUpdateAt = status.lastPollAt;
      status.error = Object.keys(status.leagueErrors).length ? `feed unavailable for ${describeLeagueErrors(status.leagueErrors)}` : null;
      notify();
    }

    async function retryFailedLeagues() {
      const failed = Object.keys(status.leagueErrors).map(Number).filter((id) => leagues.includes(id));
      if (!failed.length || now() - failedLeagueRetryAt < FAILED_LEAGUE_RETRY_MS) return;
      failedLeagueRetryAt = now();
      await loadSnapshots(failed);
    }

    // Failed leagues retry on the 30s throttle only (never a full reload every
    // tick while Unabated is down); loaded leagues poll the stream.
    async function tick() {
      if (busy || paused || !leagues.length) return;
      busy = true;
      try {
        await retryFailedLeagues();
        if (status.leaguesLoaded.length) await pollChanges();
      } catch (error) {
        setError(`feed unavailable: ${error.message}`);
        notify();
      } finally {
        busy = false;
      }
    }

    async function resync() {
      if (busy || paused || !leagues.length) return;
      busy = true;
      try {
        await loadSnapshots(leagues);
      } catch (error) {
        setError(`feed unavailable: ${error.message}`);
        notify();
      } finally {
        busy = false;
      }
    }

    function clearTimers() {
      if (pollTimer) timers.clearInterval(pollTimer);
      if (resyncTimer) timers.clearInterval(resyncTimer);
      pollTimer = null;
      resyncTimer = null;
    }

    async function start(leagueIds) {
      generation += 1;
      leagues = Array.from(new Set(leagueIds)).filter((id) => Number.isInteger(id));
      status.leagues = leagues.slice();
      clearTimers();
      state = feed.emptyState();
      status.leaguesLoaded = [];
      status.leagueErrors = {};
      status.error = null;
      cursor = null;
      paused = false;
      failedLeagueRetryAt = 0;
      if (!leagues.length) {
        status.phase = "idle";
        status.error = "no leagues enabled";
        notify();
        return;
      }
      pollTimer = timers.setInterval(tick, POLL_MS);
      resyncTimer = timers.setInterval(resync, RESYNC_MS);
      // Not resync(): an in-flight load from the old generation holds `busy`,
      // and its result is discarded anyway, so this generation loads now.
      try {
        await loadSnapshots(leagues);
      } catch (error) {
        setError(`feed unavailable: ${error.message}`);
        notify();
      }
    }

    function stop() {
      clearTimers();
      leagues = [];
      paused = false;
      status.phase = "idle";
      notify();
    }

    function pause() {
      paused = true;
    }

    // Back from hidden: continue the stream if the cursor is still fresh, else resync.
    async function resume() {
      if (!paused) return;
      paused = false;
      const lastSeen = Math.max(status.lastPollAt || 0, status.lastSnapshotAt || 0);
      if (!lastSeen || now() - lastSeen > PAUSE_RESYNC_MS) await resync();
      else await tick();
    }

    return {
      start, stop, pause, resume, tick, resync,
      getState: () => state,
      getStatus: () => ({ ...status }),
    };
  }

  const api = { createScanner, SNAPSHOT_URL, CHANGES_URL };
  if (typeof module !== "undefined" && module.exports) {
    module.exports = api;
  } else {
    root.UnabatedScanner = api;
  }
})(typeof globalThis !== "undefined" ? globalThis : this);
