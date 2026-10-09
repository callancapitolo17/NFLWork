// Bet105 for the Unabated Ticket panel (open 2026-09-29, settled 2026-10-07):
// the requests the panel makes to app.bet105.ag from Cal's own logged-in
// Chrome, and the body it POSTs to the bets service. Bet105 is a LinePros white-label behind
// Cloudflare, which challenges any request that is not the browser session's
// own — so the panel reads the account and the service (sources/bet105.py)
// parses and stores it. Pure: no DOM, no fetch, no chrome.* — loaded as a
// plain <script> in panel.html (exposes globalThis.UnabatedBet105) and via
// require() in tests/bet105.test.js. panel.js does the fetches.
//
// Read from the site's bundle and the My Bets captures of 2026-09-29 and 2026-10-07:
//   GET  /__bff/api/customers                      200 {csrfToken, customerData, ...}
//                                                  401 when not logged in
//   POST /__bff/__partner-{prematch|live}/betLobbyV2/logic/
//        {a: "getHistory", state: "0"}             200 {betGroups: [...]} — the open
//        header X-Broker-CSRF: <csrfToken>         bets; {e: "..."} on a site error
//   POST /__bff/api/wagers/search  {}              200 [wager, ...] — every wager the
//        header X-Broker-CSRF: <csrfToken>         account has, graded or not (the My
//                                                  Bets page's own call); 403
//                                                  {code: "CSRF_FAILED"} without the header
// getHistory answers the same for every other state ("1", "2", "3", "-1", ""):
// no groups, so the settled bets come from wagers/search.
// A push is all three reads or nothing (the service marks an open bet closed
// when a complete push lists it in neither, so a half push would close real bets).

(function (root) {
  "use strict";

  const ORIGIN = "https://app.bet105.ag";
  const CUSTOMERS_URL = `${ORIGIN}/__bff/api/customers`;
  const SETTLED_URL = `${ORIGIN}/__bff/api/wagers/search`;
  const FEEDS = ["prematch", "live"];
  const HISTORY_STATE_OPEN = "0";
  // The site's own My Plays view is polled by hand; every 5 min while the
  // panel is open keeps the flags current without hammering the account.
  const POLL_MS = 5 * 60 * 1000;
  // A request Bet105 never answers would hold the poll's busy flag forever and
  // stop every later read (2026-10-01: 17 h with no push), so each one gives up.
  const FETCH_TIMEOUT_MS = 30 * 1000;
  const SERVICE_PATH = "/bet105.json";
  const HTTP_UNAUTHORIZED = 401;
  const HTTP_FORBIDDEN = 403;
  const NOT_LOGGED_IN = "not logged in at app.bet105.ag — open it in this Chrome, log in, keep the panel open";
  // What the service reads of a wager (sources/bet105.py). The rest — balances,
  // the account and agent names — stays here, and the push stays far under the
  // service's 1 MB body cap (2026-10-07: 230 wagers, 170 KB cut vs 322 KB whole;
  // ~1,350 wagers fit, years at the account's pace).
  const SETTLED_WAGER_FIELDS = ["wagerId", "ticketNumber", "productCode", "category", "wagerStatus",
    "placeTime", "gradeTime", "risk", "toWin", "result", "isFreePlay", "isCashout", "wagerDetails"];
  const SETTLED_PROPERTY_FIELDS = ["odds", "fmtOdds", "grades", "teaserName", "fixedParlayName"];
  const SETTLED_LEG_FIELDS = ["legId", "eventId", "sportId", "leagueId", "league", "team1", "team2",
    "periodId", "period", "marketId", "market", "side", "figure", "fmtOdds", "startTime"];

  function historyUrl(feed) {
    return `${ORIGIN}/__bff/__partner-${feed}/betLobbyV2/logic/`;
  }

  function customersRequest() {
    return {
      method: "GET", credentials: "include", cache: "no-store", headers: { Accept: "application/json" },
      signal: AbortSignal.timeout(FETCH_TIMEOUT_MS),
    };
  }

  function historyRequest(csrfToken) {
    return {
      method: "POST", credentials: "include", cache: "no-store",
      headers: { "Content-Type": "application/json;charset=UTF-8", "X-Broker-CSRF": csrfToken },
      body: JSON.stringify({ a: "getHistory", state: HISTORY_STATE_OPEN }),
      signal: AbortSignal.timeout(FETCH_TIMEOUT_MS),
    };
  }

  function settledRequest(csrfToken) {
    return {
      method: "POST", credentials: "include", cache: "no-store",
      headers: { "Content-Type": "application/json", "X-Broker-CSRF": csrfToken },
      body: JSON.stringify({}),
      signal: AbortSignal.timeout(FETCH_TIMEOUT_MS),
    };
  }

  // What a status says about the session: the fix for a 401, Cloudflare for a 403.
  function statusError(what, status) {
    if (status === HTTP_UNAUTHORIZED) return NOT_LOGGED_IN;
    if (status === HTTP_FORBIDDEN) return `${what} refused (HTTP 403) — Cloudflare challenged the request`;
    return `${what} failed (HTTP ${status})`;
  }

  // {csrfToken} from the customers reply, or {error}.
  function csrfTokenOf(status, body) {
    if (status !== 200) return { error: statusError("Bet105 session check", status) };
    if (!body || typeof body !== "object" || typeof body.csrfToken !== "string" || !body.csrfToken) {
      return { error: "Bet105 session check answered without a csrfToken" };
    }
    return { csrfToken: body.csrfToken };
  }

  // {betGroups: [...]} from one feed's getHistory reply, or {error}.
  function betGroupsOf(feed, status, body) {
    if (status !== 200) return { error: statusError(`Bet105 ${feed} history`, status) };
    if (!body || typeof body !== "object") return { error: `Bet105 ${feed} history is not JSON` };
    if (body.e) return { error: `Bet105 ${feed} history: ${String(body.e)}` };
    const groups = body.betGroups;
    if (Array.isArray(groups)) return { betGroups: groups };
    if (groups && typeof groups === "object") return { betGroups: Object.values(groups) };
    return { error: `Bet105 ${feed} history has no betGroups` };
  }

  function pick(source, fields) {
    const out = {};
    if (!source || typeof source !== "object") return out;
    for (const field of fields) if (source[field] !== undefined) out[field] = source[field];
    return out;
  }

  // One wagers/search wager cut to SETTLED_*_FIELDS, the venue's own nesting kept.
  function settledWagerOf(wager) {
    const properties = (wager && wager.properties) || {};
    const legs = Array.isArray(properties.legs) ? properties.legs : [];
    return {
      ...pick(wager, SETTLED_WAGER_FIELDS),
      properties: { ...pick(properties, SETTLED_PROPERTY_FIELDS), legs: legs.map((leg) => pick(leg, SETTLED_LEG_FIELDS)) },
    };
  }

  // {wagers: [...]} from the wagers/search reply, or {error}. A refusal names
  // the site's own code when it sends one (CSRF_FAILED), else reads like the rest.
  function wagersOf(status, body) {
    if (status !== 200) {
      if (status !== HTTP_UNAUTHORIZED && body && typeof body === "object" && typeof body.code === "string") {
        return { error: `Bet105 settled history refused (HTTP ${status}, ${body.code})` };
      }
      return { error: statusError("Bet105 settled history", status) };
    }
    if (!Array.isArray(body)) return { error: "Bet105 settled history is not a list" };
    return { wagers: body.map(settledWagerOf) };
  }

  function pushBody(fetchedAt, groupsByFeed, settled) {
    const feeds = {};
    for (const feed of FEEDS) {
      if (!Array.isArray(groupsByFeed[feed])) throw new Error(`push needs every feed, ${feed} missing`);
      feeds[feed] = groupsByFeed[feed];
    }
    if (!Array.isArray(settled)) throw new Error("push needs the settled list");
    return { fetchedAt, feeds, settled };
  }

  function errorBody(message) {
    return { error: String(message || "unknown error") };
  }

  function isDue(lastRunAt, now) {
    return !(lastRunAt > 0) || now - lastRunAt >= POLL_MS;
  }

  const api = {
    ORIGIN, CUSTOMERS_URL, SETTLED_URL, FEEDS, POLL_MS, FETCH_TIMEOUT_MS, SERVICE_PATH, NOT_LOGGED_IN,
    historyUrl, customersRequest, historyRequest, settledRequest, csrfTokenOf, betGroupsOf,
    settledWagerOf, wagersOf, pushBody, errorBody, isDue,
  };
  if (typeof module !== "undefined" && module.exports) module.exports = api;
  else root.UnabatedBet105 = api;
})(typeof globalThis !== "undefined" ? globalThis : this);
