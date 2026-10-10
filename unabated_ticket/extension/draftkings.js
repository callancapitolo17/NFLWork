// DraftKings for the Unabated Ticket panel (2026-10-10): the panel reads the
// account's open and settled bets from Cal's own logged-in Chrome and POSTs
// them to the bets service (sources/draftkings.py parses and stores them).
// DraftKings sits behind Akamai's bot checks and asks for a code on a new
// login, so — like Bet105 — no script logs in; the browser's own session does.
// Loaded as a plain <script> in panel.html (exposes globalThis.UnabatedDraftKings)
// and via require() in tests/draftkings.test.js. No DOM and no chrome.*:
// readAccount takes fetch and WebSocket as arguments, so the tests drive it
// with fakes and panel.js passes the real ones.
//
// Read off Cal's My Bets HAR of 2026-10-10 (site US-MA-SB) and the page's own
// bundle (@draftkings/dk-my-bets-list 3.11.0). My Bets makes no HTTP call for
// the bets: it mints a token, then asks over a JSON-RPC websocket.
//   GET  https://gaming-us-ma.draftkings.com/api/wager/v1/generateEnterpriseJWT
//        (the session cookies)                      200 {token, expiresIn}
//   WSS  gateway.northamerica-northeast2.prod.dkapis.com/dkusma/shelby/api/v1/websocket
//        ?format=json&jwt=<token>
//     -> {jsonrpc: "2.0", method: "InitializeBetsPageRequest", id, params: {...}}
//        (what the page sends first; its reply is the first "All" page)
//     -> {jsonrpc: "2.0", method: "BetsRequest", id, params: {filter: {status:
//        "Open" | "Settled"}, orderCriteria: {orderBy: "placementDate",
//        direction: "DESC"}, pagination: {count, skip}, locale, ScoreboardType}}
//     <- {result: {bets: [...], events: {eventId: {...}}, cashOut, metadata}, id}
//        and, unasked, {result: {cashOutUpdate: {...}}, id} — ignored.
// The hosts carry the state (us-ma / dkusma): the account's own, Massachusetts.
// A push is both lists or nothing (the service closes an open bet a complete
// push lists nowhere, so a half push would close real bets).

(function (root) {
  "use strict";

  const JWT_URL = "https://gaming-us-ma.draftkings.com/api/wager/v1/generateEnterpriseJWT";
  const SOCKET_URL = "wss://gateway.northamerica-northeast2.prod.dkapis.com/dkusma/shelby/api/v1/websocket";
  const LISTS = ["Open", "Settled"];
  // The page asks 25 at a time; so do we.
  const PAGE_SIZE = 25;
  // 40 pages = 1,000 bets per list: far past the account's month, and a stop
  // if a reply ever ignored `skip` and repeated the same page forever.
  const MAX_PAGES_PER_LIST = 40;
  // The service keeps 30 days (config.RETENTION_DAYS); one day more so a bet
  // placed at the window's edge still gets its settlement.
  const SETTLED_LOOKBACK_MS = 31 * 24 * 60 * 60 * 1000;
  const POLL_MS = 5 * 60 * 1000;
  // A token or a reply DraftKings never sends would hold the poll's busy flag
  // forever (the Bet105 lesson of 2026-10-01), so each step gives up.
  const FETCH_TIMEOUT_MS = 30 * 1000;
  const SOCKET_TIMEOUT_MS = 30 * 1000;
  const SERVICE_PATH = "/draftkings.json";
  const HTTP_UNAUTHORIZED = 401;
  const HTTP_FORBIDDEN = 403;
  const NOT_LOGGED_IN = "not logged in at sportsbook.draftkings.com — open it in this Chrome, log in, keep the panel open";
  // What the service reads (sources/draftkings.py). The rest — cash-out
  // offers, scores, team colours — stays here, and the push stays far under
  // the service's 1 MB body cap.
  const BET_FIELDS = ["betId", "receiptId", "type", "status", "settlementStatus", "numberOfBets",
    "numberOfSelections", "displayOdds", "currency", "placementDate", "settlementDate", "stake",
    "potentialReturns", "returns"];
  const SELECTION_FIELDS = ["selectionId", "eventId", "marketId", "status", "settlementStatus", "displayOdds",
    "selectionDisplayName", "marketDisplayName"];
  const EVENT_FIELDS = ["eventId", "sportId", "leagueId", "eventStartDate", "eventDisplayName",
    "homeTeamName", "awayTeamName", "status"];

  function jwtRequest() {
    return {
      method: "GET", credentials: "include", cache: "no-store", headers: { Accept: "application/json" },
      signal: AbortSignal.timeout(FETCH_TIMEOUT_MS),
    };
  }

  function statusError(what, status) {
    if (status === HTTP_UNAUTHORIZED || status === HTTP_FORBIDDEN) return NOT_LOGGED_IN;
    return `${what} failed (HTTP ${status})`;
  }

  // {token} from the token reply, or {error}.
  function tokenOf(status, body) {
    if (status !== 200) return { error: statusError("DraftKings token", status) };
    if (!body || typeof body !== "object" || typeof body.token !== "string" || !body.token) {
      return { error: "DraftKings token reply carries no token" };
    }
    return { token: body.token };
  }

  function socketUrl(token) {
    return `${SOCKET_URL}?format=json&jwt=${encodeURIComponent(token)}`;
  }

  // The page's own first message, as the page sends it.
  function initializeMessage(id) {
    return {
      jsonrpc: "2.0", method: "InitializeBetsPageRequest", id,
      params: {
        betsRequest: { filter: { status: "All" }, pagination: { count: PAGE_SIZE, skip: 0 }, ScoreboardType: "EventScore" },
        cashOut: { information: true, pullOperations: true }, locale: "en",
      },
    };
  }

  // Settled newest-settled first (the bundle's sort keys include settlementDate),
  // so a bet held for months still lands in the window the day it settles.
  function orderByOf(list) {
    return list === "Settled" ? "settlementDate" : "placementDate";
  }

  function betsMessage(id, list, skip) {
    return {
      jsonrpc: "2.0", method: "BetsRequest", id,
      params: {
        filter: { status: list }, orderCriteria: { orderBy: orderByOf(list), direction: "DESC" },
        pagination: { count: PAGE_SIZE, skip }, locale: "en", ScoreboardType: "EventScore",
      },
    };
  }

  // {bets, events} of one reply, or {error}. `result.initial` is the
  // InitializeBetsPageRequest reply's nesting.
  function pageOf(reply) {
    if (!reply || typeof reply !== "object") return { error: "DraftKings bets reply is not JSON" };
    if (reply.error !== undefined) return { error: `DraftKings bets request refused: ${JSON.stringify(reply.error).slice(0, 300)}` };
    const result = reply.result && reply.result.initial ? reply.result.initial : reply.result;
    if (!result || !Array.isArray(result.bets)) return { error: "DraftKings bets reply carries no bets list" };
    const events = result.events && typeof result.events === "object" ? result.events : {};
    return { bets: result.bets, events };
  }

  function settledMs(bet) {
    const ms = Date.parse(bet && bet.settlementDate);
    return Number.isFinite(ms) ? ms : null;
  }

  // Done with a list: a short page, or (settled, newest-settled first) a page
  // that reaches back past the service's window.
  function isLastPage(list, bets, nowMs) {
    if (bets.length < PAGE_SIZE) return true;
    if (list !== "Settled") return false;
    const oldest = settledMs(bets[bets.length - 1]);
    return oldest !== null && oldest < nowMs - SETTLED_LOOKBACK_MS;
  }

  function pick(source, fields) {
    const out = {};
    if (!source || typeof source !== "object") return out;
    for (const field of fields) if (source[field] !== undefined) out[field] = source[field];
    return out;
  }

  function participantsOf(list) {
    return (Array.isArray(list) ? list : []).map((participant) => ({
      id: participant && participant.id, name: participant && participant.name,
      venueRole: participant && participant.venueRole,
    }));
  }

  function trimSelection(selection) {
    const out = { ...pick(selection, SELECTION_FIELDS), participants: participantsOf(selection && selection.participants) };
    // An SGP group inside a parlay (SGPx): its legs are nested, never read as one line.
    if (selection && Array.isArray(selection.nestedSGPSelections) && selection.nestedSGPSelections.length) {
      out.nestedSelectionCount = selection.nestedSGPSelections.length;
    }
    return out;
  }

  function trimBet(bet) {
    const bonus = bet && bet.bonus && typeof bet.bonus === "object" ? bet.bonus : null;
    return {
      ...pick(bet, BET_FIELDS),
      bonusType: bonus && typeof bonus.bonusType === "string" ? bonus.bonusType : null,
      freeBetAmount: bonus && Number(bonus.freeBetAmount) > 0 ? Number(bonus.freeBetAmount) : 0,
      combinationCount: Array.isArray(bet && bet.combinations) ? bet.combinations.length : 0,
      selections: (Array.isArray(bet && bet.selections) ? bet.selections : []).map(trimSelection),
    };
  }

  function trimEvent(event) {
    return { ...pick(event, EVENT_FIELDS), participants: participantsOf(event && event.participants) };
  }

  function pushBody(fetchedAt, betsByList, events) {
    const lists = {};
    for (const list of LISTS) {
      if (!Array.isArray(betsByList[list])) throw new Error(`push needs every list, ${list} missing`);
      lists[list.toLowerCase()] = betsByList[list].map(trimBet);
    }
    const trimmedEvents = {};
    for (const [eventId, event] of Object.entries(events || {})) trimmedEvents[eventId] = trimEvent(event);
    return { fetchedAt, open: lists.open, settled: lists.settled, events: trimmedEvents };
  }

  function errorBody(message) {
    return { error: String(message || "unknown error") };
  }

  function isDue(lastRunAt, now) {
    return !(lastRunAt > 0) || now - lastRunAt >= POLL_MS;
  }

  // One socket, request/reply by JSON-RPC id. Replies to other ids (the
  // cash-out pushes the page's first message subscribes to) are ignored.
  function openSocket(WebSocketImpl, url) {
    return new Promise((resolve, reject) => {
      const socket = new WebSocketImpl(url);
      const waiting = new Map();
      const failAll = (message) => {
        for (const { reject: rejectReply, timer } of waiting.values()) {
          clearTimeout(timer);
          rejectReply(new Error(message));
        }
        waiting.clear();
      };
      const openTimer = setTimeout(() => {
        reject(new Error(`DraftKings bets socket did not open within ${SOCKET_TIMEOUT_MS / 1000}s`));
        socket.close();
      }, SOCKET_TIMEOUT_MS);
      socket.onopen = () => {
        clearTimeout(openTimer);
        resolve({
          request(message) {
            return new Promise((resolveReply, rejectReply) => {
              const timer = setTimeout(() => {
                waiting.delete(message.id);
                rejectReply(new Error(`DraftKings did not answer ${message.method} within ${SOCKET_TIMEOUT_MS / 1000}s`));
              }, SOCKET_TIMEOUT_MS);
              waiting.set(message.id, { resolve: resolveReply, reject: rejectReply, timer });
              socket.send(JSON.stringify(message));
            });
          },
          close() {
            failAll("DraftKings bets socket closed");
            socket.close();
          },
        });
      };
      socket.onmessage = (event) => {
        let reply;
        try {
          reply = JSON.parse(event.data);
        } catch {
          return;
        }
        if (!reply || (reply.result && reply.result.cashOutUpdate)) return;
        const pending = waiting.get(reply.id);
        if (!pending) return;
        waiting.delete(reply.id);
        clearTimeout(pending.timer);
        pending.resolve(reply);
      };
      socket.onerror = () => {
        clearTimeout(openTimer);
        reject(new Error("DraftKings bets socket refused the connection"));
        failAll("DraftKings bets socket failed");
      };
      socket.onclose = (event) => {
        clearTimeout(openTimer);
        const code = event && event.code !== undefined ? ` (code ${event.code})` : "";
        reject(new Error(`DraftKings bets socket closed before opening${code}`));
        failAll(`DraftKings bets socket closed${code}`);
      };
    });
  }

  async function readList(socket, list, nextId, nowMs, events) {
    const bets = [];
    for (let page = 0; page < MAX_PAGES_PER_LIST; page += 1) {
      const result = pageOf(await socket.request(betsMessage(nextId(), list, page * PAGE_SIZE)));
      if (result.error) throw new Error(result.error);
      bets.push(...result.bets);
      Object.assign(events, result.events);
      if (isLastPage(list, result.bets, nowMs)) return bets;
    }
    throw new Error(`DraftKings ${list} bets ran past ${MAX_PAGES_PER_LIST} pages of ${PAGE_SIZE} — not pushed`);
  }

  // The token, the socket, the page's first message, then every page of each
  // list. Throws with the reason on any step; the socket is always closed.
  async function readAccount({ fetchImpl, WebSocketImpl, now = () => Date.now() }) {
    const tokenResponse = await fetchImpl(JWT_URL, jwtRequest());
    const session = tokenOf(tokenResponse.status, await tokenResponse.json().catch(() => null));
    if (session.error) throw new Error(session.error);
    const socket = await openSocket(WebSocketImpl, socketUrl(session.token));
    let counter = 0;
    const nextId = () => `unabated-ticket-${now()}-${(counter += 1)}`;
    try {
      const first = pageOf(await socket.request(initializeMessage(nextId())));
      if (first.error) throw new Error(first.error);
      const nowMs = now();
      const events = {};
      const betsByList = {};
      for (const list of LISTS) betsByList[list] = await readList(socket, list, nextId, nowMs, events);
      return pushBody(new Date(nowMs).toISOString(), betsByList, events);
    } finally {
      socket.close();
    }
  }

  const api = {
    JWT_URL, SOCKET_URL, LISTS, PAGE_SIZE, MAX_PAGES_PER_LIST, SETTLED_LOOKBACK_MS, POLL_MS,
    FETCH_TIMEOUT_MS, SOCKET_TIMEOUT_MS, SERVICE_PATH, NOT_LOGGED_IN,
    jwtRequest, tokenOf, socketUrl, initializeMessage, betsMessage, pageOf, isLastPage, trimBet, trimEvent,
    pushBody, errorBody, isDue, openSocket, readAccount,
  };
  if (typeof module !== "undefined" && module.exports) module.exports = api;
  else root.UnabatedDraftKings = api;
})(typeof globalThis !== "undefined" ? globalThis : this);
