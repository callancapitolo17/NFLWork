// Underdog for the Unabated Ticket extension (2026-11-09): the read of Cal's
// Underdog prediction slips, run inside his own logged-in app.underdogsports.com
// tab, and the body POSTed to the bets service (sources/underdog.py parses).
//
// Why from Chrome: Underdog's login server (login.underdogsports.com, Auth0)
// answers a script on the VM with a Cloudflare bot wall, token refreshes
// included (2026-11-09), so the VM cannot hold a login. The web app's API
// (api.underdogfantasy.com) takes the app's own 10-minute bearer token, so the
// read reuses the headers the app itself just sent (underdogpage.js captures
// them) and never touches the login: the app refreshes its own token, so its
// rotating refresh token is never reused and nothing logs out.
//
// Loaded in the Underdog tab's MAIN world before underdogpage.js (manifest
// content_scripts; exposes globalThis.UnabatedUnderdog), in background.js via
// importScripts, and via require() in tests/underdog.test.js. No DOM and no
// chrome.*: readAccount takes fetch as an argument.
//
// Endpoints (Cal's Chrome, 2026-11-07 / 11-09):
//   GET api.underdogfantasy.com/v10/user/active_entry_slips            every open slip
//   GET api.underdogfantasy.com/v10/user/settled_entry_slips?page=N    15 a page, newest
//       made first; meta {count, items, next, page}
// Headers the read needs: Authorization (the raw JWT, no "Bearer"), Client-Type,
// Client-Version (an old one answers 426 upgrade_required), Client-Device-Id.
// The app's location headers are not sent: reads answer without them.
// A settled page is ~600 KB whole (pick'em players and appearances), so each
// page is cut to its prediction slips and the lookups they name.

(function (root) {
  "use strict";

  const API_ORIGIN = "https://api.underdogfantasy.com";
  const ACTIVE_URL = `${API_ORIGIN}/v10/user/active_entry_slips`;
  const SETTLED_URL = `${API_ORIGIN}/v10/user/settled_entry_slips`;
  const APP_TAB_PATTERN = "https://app.underdogsports.com/*";
  const SERVICE_PATH = "/underdog.json";
  const PRODUCT_PREDICTION = "prediction_markets";
  // The headers copied off the app's own request; everything else (location,
  // request ids, referral links) stays with the app.
  const READ_HEADERS = ["Accept", "Client-Type", "Client-Version", "Client-Device-Id", "Authorization"];
  // A token this close to expiry is not used: the read could outlive it.
  const TOKEN_MARGIN_MS = 60 * 1000;
  // Settled pages back to slips made this long ago (the service keeps 30 days),
  // and never more than MAX_SETTLED_PAGES a read (the service refuses more:
  // sources/underdog.py MAX_SETTLED_PAGES, pinned equal by its tests).
  const SETTLED_LOOKBACK_MS = 31 * 24 * 60 * 60 * 1000;
  const MAX_SETTLED_PAGES = 20;
  const POLL_MINUTES = 2;
  // The app refreshes its token only when it makes a request; an idle tab's
  // token dies after 10 minutes. Reloading the tab makes the app refresh its
  // own login. At most this often, and never while Cal is looking at the tab.
  const RELOAD_MIN_MS = 5 * 60 * 1000;
  const RELOAD_WAIT_MS = 30 * 1000;
  const FETCH_TIMEOUT_MS = 30 * 1000;
  const HTTP_UNAUTHORIZED = 401;
  const HTTP_UPGRADE_REQUIRED = 426;
  const NO_TAB = "no app.underdogsports.com tab — open it in this Chrome, log in and leave the tab open";
  const STALE = "the Underdog tab's login has expired — reload the tab";

  // The read headers of one request the app sent, or null without a token.
  function readHeadersOf(sentHeaders) {
    const byLowerName = {};
    for (const [name, value] of Object.entries(sentHeaders || {})) byLowerName[name.toLowerCase()] = String(value);
    if (!byLowerName.authorization) return null;
    const headers = {};
    for (const name of READ_HEADERS) {
      const value = byLowerName[name.toLowerCase()];
      if (value !== undefined) headers[name] = value;
    }
    return headers;
  }

  // Epoch ms the bearer token expires (its exp claim), or null when it does not read.
  function tokenExpiryMs(authorization, decodeBase64 = root.atob) {
    try {
      const jwt = String(authorization).replace(/^Bearer\s+/i, "");
      const payload = jwt.split(".")[1].replace(/-/g, "+").replace(/_/g, "/");
      const exp = JSON.parse(decodeBase64(payload + "=".repeat((4 - (payload.length % 4)) % 4))).exp;
      return Number.isFinite(exp) ? exp * 1000 : null;
    } catch {
      return null;
    }
  }

  function isFresh(capture, nowMs) {
    return Boolean(capture && capture.headers && capture.expiresAtMs && capture.expiresAtMs - nowMs > TOKEN_MARGIN_MS);
  }

  function byId(items) {
    const map = new Map();
    for (const item of Array.isArray(items) ? items : []) if (item && item.id !== undefined) map.set(item.id, item);
    return map;
  }

  // One answer cut to its prediction slips and the options, lines, markets,
  // games and teams they name (what sources/underdog.py reads), meta kept.
  function trimPage(page) {
    const data = (page && page.data) || {};
    const slips = (Array.isArray(data.entry_slips) ? data.entry_slips : []).filter((slip) => slip && slip.product === PRODUCT_PREDICTION);
    const options = byId(data.over_under_options);
    const lines = byId(data.over_under_lines);
    const overUnders = byId(data.over_unders);
    const games = byId(data.games);
    const teams = byId(data.teams);
    const keep = { options: new Set(), lines: new Set(), overUnders: new Set(), games: new Set(), teams: new Set() };
    for (const slip of slips) {
      for (const group of slip.selection_groups || []) {
        const game = games.get(group.match_id);
        if (game) {
          keep.games.add(game.id);
          keep.teams.add(game.away_team_id);
          keep.teams.add(game.home_team_id);
        }
        for (const selection of group.selections || []) {
          const option = options.get(selection.option_id);
          if (!option) continue;
          keep.options.add(option.id);
          const line = lines.get(option.over_under_line_id);
          if (!line) continue;
          keep.lines.add(line.id);
          if (overUnders.has(line.over_under_id)) keep.overUnders.add(line.over_under_id);
        }
      }
    }
    const pickKept = (map, ids) => [...ids].filter((id) => map.has(id)).map((id) => map.get(id));
    return {
      data: {
        entry_slips: slips,
        over_under_options: pickKept(options, keep.options),
        over_under_lines: pickKept(lines, keep.lines),
        over_unders: pickKept(overUnders, keep.overUnders),
        games: pickKept(games, keep.games),
        teams: pickKept(teams, keep.teams),
      },
      meta: (page && page.meta) || {},
    };
  }

  // Epoch ms of the oldest slip made on a page, any product (they come newest first).
  function oldestCreatedMs(page) {
    const slips = (page && page.data && Array.isArray(page.data.entry_slips)) ? page.data.entry_slips : [];
    const stamps = slips.map((slip) => Date.parse(slip && slip.created_at)).filter(Number.isFinite);
    return stamps.length ? Math.min(...stamps) : null;
  }

  function statusError(what, status) {
    if (status === HTTP_UNAUTHORIZED) return STALE;
    if (status === HTTP_UPGRADE_REQUIRED) return `Underdog ${what} wants a newer app version (HTTP 426) — reload the Underdog tab`;
    return `Underdog ${what} failed (HTTP ${status})`;
  }

  // GET one answer; {page} or {error, status}.
  async function getPage(fetchImpl, url, headers, what) {
    const response = await fetchImpl(url, { method: "GET", headers, cache: "no-store", credentials: "omit",
      signal: AbortSignal.timeout(FETCH_TIMEOUT_MS) });
    if (response.status !== 200) return { error: statusError(what, response.status), status: response.status };
    const page = await response.json().catch(() => null);
    if (!page || !page.data || !Array.isArray(page.data.entry_slips)) return { error: `Underdog ${what} answered without data.entry_slips` };
    return { page };
  }

  // The whole read: the active list, then settled pages from page 1 until the
  // list ends (complete), a page's oldest slip is past the lookback, or the cap.
  // {push} or {error, stale}. A push is all of it or nothing — a half push
  // would close open slips (the service closes a held slip the pages passed).
  async function readAccount({ fetchImpl, headers, nowMs }) {
    const active = await getPage(fetchImpl, ACTIVE_URL, headers, "open slips");
    if (active.error) return { error: active.error, stale: active.status === HTTP_UNAUTHORIZED };
    const settled = [];
    let complete = false;
    for (let pageNumber = 1; pageNumber <= MAX_SETTLED_PAGES; pageNumber += 1) {
      const result = await getPage(fetchImpl, `${SETTLED_URL}?page=${pageNumber}`, headers, `settled page ${pageNumber}`);
      if (result.error) return { error: result.error, stale: result.status === HTTP_UNAUTHORIZED };
      settled.push(trimPage(result.page));
      const oldest = oldestCreatedMs(result.page);
      if (!(result.page.meta && result.page.meta.next) || oldest === null) {
        complete = true;
        break;
      }
      if (oldest < nowMs - SETTLED_LOOKBACK_MS) break;
    }
    return { push: { fetchedAt: new Date(nowMs).toISOString(), active: trimPage(active.page), settled, complete } };
  }

  // The bets service's base URL off the stored betsSettings — the panel's own
  // rule (betsview.sanitizeBetsSettings), copied because the service worker
  // cannot load betsview.js without bets/kelly/ladder/condkelly;
  // tests/underdog.test.js keeps the two equal.
  const DEFAULT_SERVICE_URL = "http://127.0.0.1:8094";
  function serviceUrlOf(storedBetsSettings) {
    const url = storedBetsSettings && storedBetsSettings.serviceUrl;
    return typeof url === "string" && /^https?:\/\/\S+$/.test(url) ? url.replace(/\/+$/, "") : DEFAULT_SERVICE_URL;
  }

  function errorBody(message) {
    return { error: String(message || "unknown error") };
  }

  // What the poll does with one tab's read: push it, revive the tab (its token
  // lapsed, the reader is not in it, or Chrome put it to sleep), or post the error.
  function stepAfterRead(result) {
    if (result && result.push) return "push";
    if (result && (result.stale || result.missing || result.unreachable)) return "reload";
    return "error";
  }

  // The error once a reload brought no usable login back. Two causes look the
  // same from here: the session is gone (log in), or the app no longer sends
  // its API calls by XMLHttpRequest, which is all underdogpage.js watches.
  function afterReloadError(result) {
    if (result && result.missing) return "the Underdog tab did not load the extension's reader — reload it";
    if (result && result.unreachable) return result.error;
    return "reloaded the Underdog tab and no login came back — log in at app.underdogsports.com in this Chrome "
      + "(if you are logged in, the app changed how it calls its API and the reader needs a fix)";
  }

  // Whether the poll may reload the Underdog tab to make the app refresh its login.
  function mayReload(tab, windowFocused, lastReloadAtMs, nowMs) {
    if (tab.active && windowFocused) return false;
    return !(lastReloadAtMs > 0) || nowMs - lastReloadAtMs >= RELOAD_MIN_MS;
  }

  const api = {
    API_ORIGIN, ACTIVE_URL, SETTLED_URL, APP_TAB_PATTERN, SERVICE_PATH, MAX_SETTLED_PAGES, POLL_MINUTES,
    RELOAD_MIN_MS, RELOAD_WAIT_MS, FETCH_TIMEOUT_MS, NO_TAB, STALE, serviceUrlOf, readHeadersOf, tokenExpiryMs,
    isFresh, trimPage, readAccount, errorBody, stepAfterRead, afterReloadError, mayReload,
  };
  if (typeof module !== "undefined" && module.exports) module.exports = api;
  else root.UnabatedUnderdog = api;
})(typeof globalThis !== "undefined" ? globalThis : this);
