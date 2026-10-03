// Unabated Ticket — cross-league edge scanner loop.
//
// Runs inside the side panel page (never the service worker): the panel is
// open whenever you are betting, setInterval is reliable there, and closing
// the panel stops every request. Snapshot per enabled league on start, then
// each league's snapshot is re-downloaded on its own cadence by file size.
// Snapshots are the only source: the polling changes stream
// (api-k.unabated.com/api/markets/changes/query) has answered HTTP 410
// "Polling change subscriptions have been retired" since 2026-09-27, and its
// SSE replacement needs the logged-in session.
//
// Side effects: network requests to content.unabated.com only. No storage, no
// DOM. State lives in memory and is handed to the panel through
// onChange(status, state, history) — `history` is the per-line record of what
// each line was worth on every snapshot (edgemove.js, #132), also in memory only. `status.observingSince`
// says since when that record has no gap, and `status.leagueObservingSince`
// the same per league (fillfair.js reads a bet's fill off it only when the
// bet was placed after both). `status.leagueLoadedAt` is when each league's
// last snapshot landed.
//
// Loaded as a plain <script> in panel.html (globalThis.UnabatedScanner) and via
// require() in tests, where fetch and timers are injected.

(function (root) {
  "use strict";

  const feed = typeof module !== "undefined" && module.exports ? require("./feed.js") : root.UnabatedFeed;
  const edgemove = typeof module !== "undefined" && module.exports ? require("./edgemove.js") : root.UnabatedEdgeMove;

  const SNAPSHOT_BASE_URL = (leagueId) => `https://content.unabated.com/markets/v2/league/${leagueId}/odds.json`;
  // CloudFront serves the bare URL to gzip clients (every browser) from an
  // edge cache with no TTL from the origin: measured 2026-09-10 a copy 6.5 h
  // old ("age: 23253", last-modified 16:14 at 22:42 UTC) while curl without
  // gzip and any query string got the fresh file. A query that changes
  // every CACHE_BUST_SEC forces a miss, so the panel never prices off a
  // stale edge copy; within that window the browser cache still serves it.
  const CACHE_BUST_SEC = 30;
  const SNAPSHOT_URL = (leagueId, atMs) => `${SNAPSHOT_BASE_URL(leagueId)}?t=${Math.floor((atMs ?? Date.now()) / (CACHE_BUST_SEC * 1000))}`;
  // How often the loop checks which leagues are due a re-download.
  const TICK_MS = 10000;
  // Full resync (books, teams) — the per-league refresh below is what keeps
  // prices current.
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
  // A snapshot regenerates every ~27s; a build older than this is a stale
  // CDN copy (or an off-season league that stopped regenerating).
  const STALE_BUILD_MS = 15 * 60 * 1000;
  // After a pause longer than this the board is too old to top up league by
  // league: resync everything instead.
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

  function createScanner(deps) {
    const fetchImpl = (deps && deps.fetchImpl) || ((...args) => root.fetch(...args));
    const now = (deps && deps.now) || (() => Date.now());
    const timers = (deps && deps.timers) || { setInterval: (f, ms) => root.setInterval(f, ms), clearInterval: (id) => root.clearInterval(id) };
    const onChange = (deps && deps.onChange) || (() => {});

    let state = feed.emptyState();
    // line key -> observations (edgemove.observe); reset with the state.
    let history = {};
    let leagues = [];
    let tickTimer = null;
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
      lineCount: 0,
      altLineCount: 0,
      eventCount: 0,
      loading: null, // {done, total} while snapshots are downloading
      snapshotBuiltAt: null, // newest Last-Modified among loaded leagues
      staleLeagues: [], // leagues whose Last-Modified is older than STALE_BUILD_MS (a stale edge copy, or an off-season file)
      // Since when the history is a record of the board with no gap: the end
      // of the first successful observation after start(), after resume(),
      // or after more than PAUSE_RESYNC_MS without one (Unabated down, no
      // network). Not the time resume() is called: until the catch-up lands
      // the history still holds the pre-pause board. Null while paused.
      observingSince: null,
      // leagueId -> since when that league's snapshots have loaded without a
      // gap. One league can keep failing while the others land, which the
      // clock above cannot see.
      leagueObservingSince: {},
      // leagueId -> when its last snapshot landed: a reader of one league's
      // lines (the Teasers tab reads NFL and CFB) redoes its work only when
      // that league changed, not on every league's refresh.
      leagueLoadedAt: {},
    };
    // The last successful snapshot load, paused or not.
    let lastObservedAt = null;
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
        throw new Error(controller && controller.signal.aborted ? `timed out after ${timeoutMs / 1000}s` : error.message, { cause: error });
      } finally {
        if (timer) clearTimeout(timer);
      }
    }

    function notify() {
      status.lineCount = feed.countLines(state);
      status.altLineCount = feed.countAltLines(state);
      status.eventCount = Object.keys(state.events).length;
      onChange({ ...status }, state, history);
    }

    function setError(message) {
      status.error = message;
      status.phase = status.leaguesLoaded.length ? "live" : "error";
    }

    // A snapshot load succeeded. A pass that lands while
    // paused (it was in flight when the panel hid) moves the clock but never
    // starts a run: the panel is not watching.
    function markObserved() {
      const at = now();
      const gap = lastObservedAt == null || at - lastObservedAt > PAUSE_RESYNC_MS;
      if (!paused && (status.observingSince == null || gap)) status.observingSince = at;
      lastObservedAt = at;
    }

    // One league's snapshot landed: its run goes on, or starts over when the
    // previous load was more than one refresh interval plus PAUSE_RESYNC_MS
    // ago (its loads were failing).
    function markLeagueObserved(leagueId, previousLoadAt, bytes) {
      const at = now();
      const held = status.leagueObservingSince[leagueId];
      const gap = held == null || previousLoadAt == null || at - previousLoadAt > refreshEveryMs(bytes) + PAUSE_RESYNC_MS;
      status.leagueObservingSince = { ...status.leagueObservingSince, [leagueId]: gap ? at : held };
    }

    async function fetchSnapshot(leagueId) {
      // no-cache = revalidate with If-None-Match; a 304 serves the cached body.
      const response = await fetchWithTimeout(SNAPSHOT_URL(leagueId, now()), { cache: "no-cache", credentials: "include" }, SNAPSHOT_TIMEOUT_MS);
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
    // mergeStates per league would copy every line 29 times). Returns the
    // keys of the league's lines it dropped, re-listed or not.
    function mergeInto(target, loaded) {
      target.leagues = target.leagues.filter((id) => !loaded.leagues.includes(id)).concat(loaded.leagues);
      Object.assign(target.teams, loaded.teams);
      Object.assign(target.teamIndex ||= {}, loaded.teamIndex || {});
      Object.assign(target.books, loaded.books);
      for (const [id, event] of Object.entries(target.events)) if (loaded.leagues.includes(event.leagueId)) delete target.events[id];
      const dropped = [];
      for (const [key, line] of Object.entries(target.lines)) {
        if (!loaded.leagues.includes(line.leagueId)) continue;
        delete target.lines[key];
        dropped.push(key);
      }
      Object.assign(target.events, loaded.events);
      Object.assign(target.lines, loaded.lines);
      return dropped;
    }

    // One observation per line of a freshly merged snapshot, and the history
    // of every line the snapshot no longer lists is forgotten: a pulled line
    // must not keep a stale record for the session.
    function recordSnapshot(loaded, droppedKeys) {
      const at = now();
      for (const line of Object.values(loaded.lines)) edgemove.observe(history, line, { at });
      for (const key of droppedKeys) if (!state.lines[key]) edgemove.forget(history, key);
    }

    // Load every requested league, publishing each one the moment it lands
    // (the panel fills in league by league); keep whatever succeeds and
    // report the rest by name.
    async function loadSnapshots(leagueIds) {
      const startedUnder = generation;
      const fullLoad = leagueIds.length >= leagues.length;
      const ordered = loadOrder(leagueIds);
      // Progress counter only for a full load; a background refresh is silent.
      if (fullLoad) status.loading = { done: 0, total: ordered.length };
      if (!status.leaguesLoaded.length || fullLoad) status.phase = status.leaguesLoaded.length ? status.phase : "loading";
      const errors = { ...status.leagueErrors };
      let loadedCount = 0;
      await mapWithConcurrency(ordered, SNAPSHOT_CONCURRENCY, async (leagueId) => {
        let loaded;
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
        const previousLoadAt = leagueMeta[leagueId] ? leagueMeta[leagueId].loadedAt : null;
        leagueMeta[leagueId] = { loadedAt: now(), bytes: loaded.bytes, builtAt: loaded.builtAt };
        status.leagueLoadedAt = { ...status.leagueLoadedAt, [leagueId]: leagueMeta[leagueId].loadedAt };
        if (loaded.builtAt != null) {
          status.snapshotBuiltAt = Math.max(status.snapshotBuiltAt || 0, loaded.builtAt);
        }
        status.staleLeagues = leagues.filter((id) => leagueMeta[id] && leagueMeta[id].builtAt != null && now() - leagueMeta[id].builtAt > STALE_BUILD_MS);
        // In place, full load or refresh alike: the list never collapses to
        // one league while the others are still downloading.
        const dropped = mergeInto(state, loaded.state);
        recordSnapshot(loaded.state, dropped);
        markLeagueObserved(leagueId, previousLoadAt, loaded.bytes);
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
      markObserved();
      status.phase = "live";
      status.error = Object.keys(errors).length ? `feed unavailable for ${describeLeagueErrors(errors)}` : null;
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

    async function retryFailedLeagues() {
      const failed = Object.keys(status.leagueErrors).map(Number).filter((id) => leagues.includes(id));
      if (!failed.length || now() - failedLeagueRetryAt < FAILED_LEAGUE_RETRY_MS) return;
      failedLeagueRetryAt = now();
      await loadSnapshots(failed);
    }

    // Failed leagues retry on the 30s throttle only (never a full reload every
    // tick while Unabated is down); leagues past their refresh cadence are
    // re-downloaded.
    async function tick() {
      if (busy || paused || !leagues.length || status.loading) return;
      busy = true;
      try {
        await retryFailedLeagues();
        const due = leaguesDueForRefresh();
        if (due.length && due.length < leagues.length) await loadSnapshots(due);
        else if (due.length) await loadSnapshots(leagues);
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
      if (tickTimer) timers.clearInterval(tickTimer);
      if (resyncTimer) timers.clearInterval(resyncTimer);
      tickTimer = null;
      resyncTimer = null;
    }

    async function start(leagueIds) {
      generation += 1;
      leagues = Array.from(new Set(leagueIds)).filter((id) => Number.isInteger(id));
      status.leagues = leagues.slice();
      clearTimers();
      state = feed.emptyState();
      history = {};
      status.observingSince = null;
      status.leagueObservingSince = {};
      status.leagueLoadedAt = {};
      lastObservedAt = null;
      status.leaguesLoaded = [];
      status.leagueErrors = {};
      status.error = null;
      paused = false;
      failedLeagueRetryAt = 0;
      leagueMeta = {};
      if (!leagues.length) {
        status.phase = "idle";
        status.error = "no leagues enabled";
        notify();
        return;
      }
      tickTimer = timers.setInterval(tick, TICK_MS);
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
      status.observingSince = null;
      notify();
    }

    // The history is kept but has a gap from here until the first
    // observation after resume().
    function pause() {
      paused = true;
      status.observingSince = null;
    }

    // Back from hidden: top up the leagues that are due if the board is recent, else resync.
    async function resume() {
      if (!paused) return;
      paused = false;
      const lastSeen = status.lastSnapshotAt || 0;
      if (!lastSeen || now() - lastSeen > PAUSE_RESYNC_MS) await resync();
      else await tick();
    }

    return {
      start, stop, pause, resume, tick, resync,
      getState: () => state,
      getHistory: () => history,
      getStatus: () => ({ ...status }),
      refreshEveryMs,
    };
  }

  const api = { createScanner, SNAPSHOT_URL, SNAPSHOT_BASE_URL };
  if (typeof module !== "undefined" && module.exports) {
    module.exports = api;
  } else {
    root.UnabatedScanner = api;
  }
})(typeof globalThis !== "undefined" ? globalThis : this);
