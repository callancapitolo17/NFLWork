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
// pageReady, booksFilter, locateResult}; forwards {locate} requests (row or
// notification click) and the stored {ticket} (for the watcher to resume
// after a navigation) to page.js via window.postMessage. None on the page.

(function () {
  "use strict";

  // A re-injected copy (extension reloaded with this tab open) replaces the
  // dead one; if a live one is somehow present, it must not double-write.
  if (window.__unabatedTicketContentActive) return;
  window.__unabatedTicketContentActive = true;

  const MESSAGE_SOURCE = "unabated-ticket";
  const HANDLED_TYPES = new Set(["ticket", "watch", "error", "ready", "filters", "located", "resume_request"]);
  // A locate request older than this is left alone (the tab it targeted may
  // have been reloaded long after the click).
  const LOCATE_MAX_AGE_MS = 90 * 1000;

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

  // Unabated book selection, the Edges tab's default book filter. A failed read keeps the last
  // good filter (a tab mid-load must not blank it) but records why, so the
  // panel can say the filter is stale rather than pretend it is current.
  function handleFilters(payload) {
    if (payload.error) {
      chrome.storage.local.get("booksFilter", (stored) => {
        if (chrome.runtime.lastError) return;
        const previous = stored.booksFilter || null;
        setSession({ booksFilter: { ...(previous || {}), bookIds: previous ? previous.bookIds : null, lastError: payload.error, lastErrorAt: payload.at, at: previous ? previous.at : null, debug: payload.debug || (previous ? previous.debug : null) } });
      });
      return;
    }
    setSession({ booksFilter: { bookIds: payload.bookIds, url: payload.url, at: payload.at, lastError: null, lastErrorAt: null, debug: payload.debug || null } });
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

  function handleLocated(payload) {
    setSession({ locateResult: payload });
  }

  // The stored ticket, handed to page.js so a freshly loaded copy resumes
  // watching it (page.js dies with every navigation; the ticket does not).
  // Offered once on load and again on request, since the two scripts load
  // in no guaranteed order; page.js ignores a ticket it already watches.
  function offerStoredTicket() {
    chrome.storage.local.get("ticket", (stored) => {
      if (chrome.runtime.lastError || !stored.ticket) return;
      window.postMessage({ source: MESSAGE_SOURCE, type: "resume", payload: stored.ticket }, window.location.origin);
    });
  }

  const handlers = { ticket: handleTicket, error: handleError, watch: handleWatch, ready: handleReady, filters: handleFilters, located: handleLocated, resume_request: offerStoredTicket };

  function forwardLocate(locate) {
    if (!locate || typeof locate.at !== "number" || Date.now() - locate.at > LOCATE_MAX_AGE_MS) return;
    window.postMessage({ source: MESSAGE_SOURCE, type: "locate", payload: locate }, window.location.origin);
  }

  // Requests arrive live while this tab is open, or are picked up on load when
  // the panel had to navigate/open the tab (page.js waits for the grid).
  try {
    chrome.storage.onChanged.addListener((changes, area) => {
      if (area === "local" && changes.locate && changes.locate.newValue) forwardLocate(changes.locate.newValue);
    });
    chrome.storage.local.get(["locate", "locateResult"], (stored) => {
      if (chrome.runtime.lastError) return;
      const locate = stored.locate;
      const result = stored.locateResult;
      // Already answered (a plain reload within the max age): do not flash again.
      if (locate && result && result.key === locate.key && result.at >= locate.at) return;
      forwardLocate(locate);
    });
    offerStoredTicket();
  } catch (_error) {
    // Extension context invalidated.
  }

  console.info("[unabated-ticket] content.js active (direct-to-storage)");
  // Reloading the extension leaves this copy running in an open tab with a
  // dead chrome.* handle ("Extension context invalidated" on every watch
  // tick). One such failure retires this copy; the reloaded extension injects
  // a live one when the tab is reloaded.
  function onPageMessage(event) {
    if (event.source !== window) return;
    const data = event.data;
    if (!data || data.source !== MESSAGE_SOURCE || !HANDLED_TYPES.has(data.type)) return;
    const handler = handlers[data.type];
    if (!handler) return;
    try {
      handler(data.payload);
    } catch (error) {
      if (!/context invalidated/i.test(String(error && error.message))) throw error;
      window.removeEventListener("message", onPageMessage);
      window.__unabatedTicketContentActive = false;
      console.info("[unabated-ticket] content.js retired: extension was reloaded; reload this tab");
    }
  }
  window.addEventListener("message", onPageMessage);
})();
