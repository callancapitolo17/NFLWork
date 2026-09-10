// Unabated Ticket — cross-league edge scanner loop.
//
// Runs inside the side panel page (never the service worker): the panel is
// open whenever you are betting, setInterval is reliable there, and closing
// the panel stops every request. Snapshot per enabled league on start, then
// the changes stream every POLL_MS from the last cursor — and, because the
// anonymous stream is INCOMPLETE (measured 2026-09-10 over 3 min: 69 of 191
// NFL line changes delivered; Kalshi, Caesars, ProphetX, Polymarket and
// Underdog moves mostly missing), each league's snapshot is re-downloaded on
// its own cadence by file size so exchange prices never sit stale for long.
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
  // Full resync (books, teams, cursor) — the per-league refresh below is what
  // keeps prices current.
  const RESYNC_MS = 10 * 60 * 1000;
  // Per-league snapshot refresh cadence by compressed size: the snapshot is
  // regenerated every ~27s and is the only complete source, so small files
  // refresh often and CFB (9.7 MB) least. Unknown size -> the middle tier.
  const REFRESH_TIERS = [
    { maxBytes: 2 * 1024 * 1024, everyMs: 60 * 1000 },
    { maxBytes: 5 * 1024 * 1024, everyMs: 120 * 1000 },
    { maxBytes: Infinity, everyMs: 300 * 1000 },
  ];
  const UNKNOWN_SIZE_BYTES = 3 * 1024 * 1024;
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
  // ~27 leagues, 18 MB gzip per resync; a few at a time keeps peak memory
  // (each body is parsed in full) and the JSON.parse stalls bounded.
  const SNAPSHOT_CONCURRENCY = 4;
  // Loaded last: CFB is half the bytes of a full load (9.7 of ~18 MB gzip);
  // everything else shows up while it downloads.
  const LOAD_LAST_LEAGUE_IDS = [2];
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
      loading: null, // {done, total} while snapshots are downloading
    };
    let failedLeagueRetryAt = 0;
    // Bumped by start(); a load that began under an older generation is discarded.
    let generation = 0;
    // leagueId -> {loadedAt, bytes}: drives the per-league refresh cadence.
    let leagueMeta = {};

    function refreshEveryMs(bytes) {
      return REFRESH_TIERS.find((tier) => bytes <= tier.maxBytes).everyMs;
    }

    function leaguesDueForRefresh() {
      return leagues.filter((leagueId) => {
        const meta = leagueMeta[leagueId];
        return meta && now() - meta.loadedAt >= refreshEveryMs(meta.bytes);
      });
    }

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
      const contentLength = response.headers && response.headers.get ? Number(response.headers.get("content-length")) : NaN;
      const bytes = Number.isFinite(contentLength) && contentLength > 0 ? contentLength : UNKNOWN_SIZE_BYTES;
      return { state: feed.parseSnapshot(json, { leagueId }), builtAt: Number.isFinite(builtAt) ? builtAt : null, bytes };
    }

    function loadOrder(leagueIds) {
      const last = leagueIds.filter((id) => LOAD_LAST_LEAGUE_IDS.includes(id));
      return leagueIds.filter((id) => !LOAD_LAST_LEAGUE_IDS.includes(id)).concat(last);
    }

    // Merge one league's snapshot into the live state in place (a full
    // mergeStates per league would copy every line 29 times).
    function mergeInto(target, loaded) {
      target.leagues = target.leagues.filter((id) => !loaded.leagues.includes(id)).concat(loaded.leagues);
      Object.assign(target.teams, loaded.teams);
      Object.assign(target.books, loaded.books);
      for (const [id, event] of Object.entries(target.events)) if (loaded.leagues.includes(event.leagueId)) delete target.events[id];
      for (const [key, line] of Object.entries(target.lines)) if (loaded.leagues.includes(line.leagueId)) delete target.lines[key];
      Object.assign(target.events, loaded.events);
      Object.assign(target.lines, loaded.lines);
    }

    // Load every requested league, publishing each one the moment it lands
    // (the panel fills in league by league); keep whatever succeeds and
    // report the rest by name. The changes cursor starts at the OLDEST
    // snapshot build time when that is recent enough, so nothing between
    // build and first poll is missed.
    async function loadSnapshots(leagueIds) {
      const startedUnder = generation;
      const fullLoad = leagueIds.length >= leagues.length;
      const target = fullLoad ? feed.emptyState() : state;
      const ordered = loadOrder(leagueIds);
      // Progress counter only for a full load; a background refresh is silent.
      if (fullLoad) status.loading = { done: 0, total: ordered.length };
      if (!status.leaguesLoaded.length || fullLoad) status.phase = status.leaguesLoaded.length ? status.phase : "loading";
      const errors = { ...status.leagueErrors };
      let oldestBuild = null;
      let loadedCount = 0;
      await mapWithConcurrency(ordered, SNAPSHOT_CONCURRENCY, async (leagueId) => {
        let loaded = null;
        try {
          loaded = await fetchSnapshot(leagueId);
        } catch (error) {
          if (startedUnder !== generation) return;
          errors[leagueId] = error.message;
          if (status.loading) status.loading = { done: status.loading.done + 1, total: ordered.length };
          return;
        }
        // start() ran meanwhile: this league is no longer what the panel wants.
        if (startedUnder !== generation) return;
        delete errors[leagueId];
        loadedCount += 1;
        leagueMeta[leagueId] = { loadedAt: now(), bytes: loaded.bytes };
        if (loaded.builtAt != null) oldestBuild = oldestBuild == null ? loaded.builtAt : Math.min(oldestBuild, loaded.builtAt);
        mergeInto(target, loaded.state);
        if (fullLoad && target !== state) {
          // First league of a full (re)load: switch to the fresh state now so
          // the panel shows it, and the rest merge into it as they land.
          state = target;
        }
        status.leaguesLoaded = Array.from(new Set(state.leagues)).sort((a, b) => a - b);
        status.leagueErrors = { ...errors };
        if (status.loading) status.loading = { done: status.loading.done + 1, total: ordered.length };
        notify();
      });
      if (startedUnder !== generation) return false;
      status.loading = null;
      status.leagueErrors = errors;
      if (loadedCount === 0) {
        setError(`feed unavailable: ${describeLeagueErrors(errors)}`);
        notify();
        return false;
      }
      status.lastSnapshotAt = now();
      status.phase = "live";
      status.error = Object.keys(errors).length ? `feed unavailable for ${describeLeagueErrors(errors)}` : null;
      // Only a full load restarts the stream at the snapshot build time; a
      // per-league refresh leaves the cursor where the stream is.
      if (fullLoad) {
        const recent = oldestBuild != null && now() - oldestBuild <= CURSOR_MAX_AGE_MS;
        cursor = recent ? feed.cursorFromDate(oldestBuild) : null;
      }
      notify();
      return true;
    }

    async function mapWithConcurrency(items, width, worker) {
      const results = new Array(items.length);
      let next = 0;
      async function lane() {
        while (next < items.length) {
          const index = next;
          next += 1;
          results[index] = await worker(items[index]);
        }
      }
      await Promise.all(Array.from({ length: Math.min(width, items.length) }, lane));
      return results;
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
    // tick while Unabated is down); leagues past their refresh cadence are
    // re-downloaded; loaded leagues poll the stream.
    async function tick() {
      if (busy || paused || !leagues.length || status.loading) return;
      busy = true;
      try {
        await retryFailedLeagues();
        const due = leaguesDueForRefresh();
        if (due.length && due.length < leagues.length) await loadSnapshots(due);
        else if (due.length) await loadSnapshots(leagues);
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
      leagueMeta = {};
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
      refreshEveryMs,
    };
  }

  const api = { createScanner, SNAPSHOT_URL, CHANGES_URL };
  if (typeof module !== "undefined" && module.exports) {
    module.exports = api;
  } else {
    root.UnabatedScanner = api;
  }
})(typeof globalThis !== "undefined" ? globalThis : this);
