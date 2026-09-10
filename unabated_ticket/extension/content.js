// Unabated Ticket — isolated-world bridge.
//
// page.js (MAIN world) reads React fiber but cannot use chrome.* APIs, so it
// hands data here via window.postMessage. This script (ISOLATED world) writes
// chrome.storage.local directly — the side panel reads it and listens for
// storage.onChanged. background.js only sets the storage access level and the
// panel-open-on-click behavior; it is deliberately NOT on the hot path, so a
// dead/asleep service worker cannot stop a ticket from reaching the panel.
//
// Side effects: writes chrome.storage.local {ticket, error, watchStatus,
// pageReady, booksFilter}. None on the page.

(function () {
  "use strict";

  const MESSAGE_SOURCE = "unabated-ticket";
  const HANDLED_TYPES = new Set(["ticket", "watch", "error", "ready", "filters"]);

  function setSession(obj) {
    try {
      chrome.storage.local.set(obj, () => {
        if (chrome.runtime.lastError) console.info("[unabated-ticket] storage write failed:", chrome.runtime.lastError.message);
      });
    } catch (_error) {
      // Extension context invalidated (reloaded while the page stayed open).
    }
  }

  function sameCurrent(a, b) {
    if (a == null || b == null) return a == null && b == null;
    return a.price === b.price && a.points === b.points && a.fair === b.fair && a.offBoard === b.offBoard;
  }

  function handleTicket(ticket) {
    setSession({ ticket, error: null, watchStatus: null });
  }

  function handleError(payload) {
    setSession({ ticket: null, error: payload, watchStatus: null });
  }

  function handleReady(payload) {
    setSession({ pageReady: { url: payload.url, at: payload.at } });
  }

  // Books/bet-type filter for the Edges tab. A failed read keeps the last
  // good filter (a tab mid-load must not blank it) but records why, so the
  // panel can say the filter is stale rather than pretend it is current.
  function handleFilters(payload) {
    if (payload.error) {
      chrome.storage.local.get("booksFilter", (stored) => {
        if (chrome.runtime.lastError) return;
        const previous = stored.booksFilter || null;
        setSession({ booksFilter: { ...(previous || {}), bookIds: previous ? previous.bookIds : null, betTypeIds: previous ? previous.betTypeIds : null, lastError: payload.error, lastErrorAt: payload.at, at: previous ? previous.at : null } });
      });
      return;
    }
    setSession({ booksFilter: { bookIds: payload.bookIds, betTypeIds: payload.betTypeIds, betTypeReason: payload.betTypeReason, url: payload.url, at: payload.at, lastError: null, lastErrorAt: null } });
  }

  function handleWatch(payload) {
    chrome.storage.local.get("ticket", (stored) => {
      if (chrome.runtime.lastError) return;
      const ticket = stored.ticket;
      if (!ticket || ticket.capturedAt !== payload.capturedAt) return; // stale watcher
      if (payload.error) {
        setSession({ watchStatus: { capturedAt: ticket.capturedAt, seenAt: Date.now(), error: payload.error } });
        return;
      }
      const line = payload.line;
      const moved = ticket.price !== line.price || ticket.points !== line.points || line.offBoard;
      const current = moved ? line : null;
      const updates = { watchStatus: { capturedAt: ticket.capturedAt, seenAt: line.seenAt, error: null } };
      // Only rewrite the ticket when `current` changes, so the panel does not re-render every 5s.
      if (!sameCurrent(ticket.current, current)) updates.ticket = { ...ticket, current };
      setSession(updates);
    });
  }

  const handlers = { ticket: handleTicket, error: handleError, watch: handleWatch, ready: handleReady, filters: handleFilters };

  console.info("[unabated-ticket] content.js active (direct-to-storage)");
  window.addEventListener("message", (event) => {
    if (event.source !== window) return;
    const data = event.data;
    if (!data || data.source !== MESSAGE_SOURCE || !HANDLED_TYPES.has(data.type)) return;
    const handler = handlers[data.type];
    if (handler) handler(data.payload);
  });
})();
