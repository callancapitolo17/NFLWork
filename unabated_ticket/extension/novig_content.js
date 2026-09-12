// Unabated Ticket — Novig isolated-world bridge (#116).
//
// novig_page.js (MAIN world) mirrors the Portfolio responses the Novig app
// fetches for itself and posts them here. This script (ISOLATED world, loaded
// after novig_bets.js) keeps the pages of each list seen in this tab,
// normalises them with UnabatedNovigBets and writes the result to
// chrome.storage.local, where the side panel reads it and listens for
// storage.onChanged — the same direct-to-storage pattern as content.js on
// Unabated (a service worker is never on the hot path). No request is made
// from here.
//
// Side effects: writes chrome.storage.local {betsNovig: {bets, readAt, url,
// error, complete, pageSeenAt}}. `bets` are the records of every page this
// tab has seen since load; `readAt` is the time of the last portfolio
// response; `complete` says every list was seen to its end (a short last
// page), which lets the panel treat the read as authoritative for the venue;
// `pageSeenAt` is the last time a Novig tab was open at all. None on the page.

(function () {
  "use strict";

  if (window.__unabatedTicketNovigContentActive) return;
  window.__unabatedTicketNovigContentActive = true;

  const MESSAGE_SOURCE = "unabated-ticket-novig";
  const novig = globalThis.UnabatedNovigBets;
  const pages = new Map();
  let readAt = null;

  function setLocal(obj) {
    try {
      chrome.storage.local.set(obj, () => {
        if (chrome.runtime.lastError) console.info("[unabated-ticket] novig storage write failed:", chrome.runtime.lastError.message);
      });
    } catch (_error) {
      // Extension context invalidated (reloaded while the page stayed open).
    }
  }

  // Everything this tab has seen, normalised, with `error` when the latest
  // response could not be used — the held pages are never blanked by it.
  function publish(error, url, at) {
    const collected = pages.size ? novig.collectPages(pages) : { orders: [], parlays: [], complete: false };
    let bets = [];
    try {
      bets = novig.normalizeNovig({ orders: collected.orders, parlays: collected.parlays, readAt });
    } catch (normaliseError) {
      console.warn("[unabated-ticket] novig normalise failed:", normaliseError);
      error = `normalise failed: ${normaliseError.message}`;
    }
    setLocal({ betsNovig: { bets, readAt, url, error, complete: collected.complete && !error, pageSeenAt: at } });
  }

  function handlePortfolio(payload) {
    if (payload.error) {
      publish(`${payload.operationName}: ${payload.error}`, payload.url, payload.at);
      return;
    }
    if (!novig.applyResponse(pages, payload)) return;
    readAt = payload.at;
    publish(null, payload.url, payload.at);
  }

  window.addEventListener("message", (event) => {
    if (event.source !== window) return;
    const data = event.data;
    if (!data || data.source !== MESSAGE_SOURCE || data.type !== "portfolio") return;
    handlePortfolio(data.payload || {});
  });

  // A Novig tab is open: the panel can tell "no tab" from "tab open, no
  // Portfolio screen fetched yet" without a message from the main world.
  try {
    chrome.storage.local.get("betsNovig", (stored) => {
      if (chrome.runtime.lastError) return;
      const previous = stored.betsNovig && typeof stored.betsNovig === "object" ? stored.betsNovig : {};
      setLocal({ betsNovig: { bets: Array.isArray(previous.bets) ? previous.bets : [], readAt: previous.readAt ?? null, url: window.location.href, error: previous.error ?? null, complete: previous.complete === true, pageSeenAt: new Date().toISOString() } });
    });
  } catch (_error) {
    // Extension context invalidated.
  }

  console.info("[unabated-ticket] novig_content.js active (direct-to-storage)");
})();
