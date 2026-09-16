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
// pageReady, booksFilter, locateResult, pageCheck}; forwards {locate} requests
// (row or notification click) and the stored {ticket} (for the watcher to
// resume after a navigation) to page.js via window.postMessage. On the page it
// leaves one global, __unabatedTicketContentRetire, the handoff a re-injected
// copy calls to retire the copy before it.

(function () {
  "use strict";

  // A re-injected copy (extension reloaded with this tab open) replaces the
  // copy before it. Chrome keys the isolated world on the extension id, which
  // a reload does not change, so both copies share this `window`: the old copy
  // publishes how to retire it and the new one calls that, synchronously,
  // before installing its own listeners. A boolean guard here instead would
  // make the NEW copy return and leave the dead one — whose chrome.* handle
  // throws "Extension context invalidated" on every write — registered on
  // `message`, which is how the tab used to need a reload (README, Install).
  if (typeof window.__unabatedTicketContentRetire === "function") {
    try {
      window.__unabatedTicketContentRetire();
    } catch (error) {
      console.info("[unabated-ticket] previous content.js did not retire cleanly:", error && error.message);
    }
  }

  let retired = false;

  const MESSAGE_SOURCE = "unabated-ticket";
  const HANDLED_TYPES = new Set(["ticket", "watch", "error", "ready", "filters", "located", "resume_request", "pagecheck"]);
  // A locate request older than this is left alone (the tab it targeted may
  // have been reloaded long after the click).
  const LOCATE_MAX_AGE_MS = 90 * 1000;

  function setSession(obj) {
    if (retired) return;
    try {
      chrome.storage.local.set(obj, () => {
        if (chrome.runtime.lastError) console.info("[unabated-ticket] storage write failed:", chrome.runtime.lastError.message);
      });
    } catch (_error) {
      // Extension context invalidated (reloaded while the page stayed open).
    }
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

  // Two Unabated tabs can both watch the resumed ticket; a tab that reads
  // the line outranks one that cannot (other game date, filtered row), so a
  // failure is dropped while a good read from within the last interval and
  // a half stands. Otherwise the panel would flap between the two every 5 s.
  const GOOD_READ_OUTRANKS_MS = 7500;

  // The watcher writes ONLY `watchStatus`, never the ticket. It used to do a
  // read-modify-write of the whole ticket to hang `current` off it, and the
  // capturedAt guard below tests the ticket it READ, not the one in storage
  // when the write lands: a click inside that get->set window was written by
  // handleTicket and then reverted ~1 ms later by the old ticket coming back,
  // leaving a stale ticket whose watcher was already gone. The two writers now
  // touch different values, so a fresh capture always stands. A late tick for
  // the previous capture can still land here, which is why every reader
  // (panel.js) matches watchStatus.capturedAt against the ticket's.
  function handleWatch(payload) {
    chrome.storage.local.get(["ticket", "watchStatus"], (stored) => {
      if (chrome.runtime.lastError) return;
      const ticket = stored.ticket;
      if (!ticket || ticket.capturedAt !== payload.capturedAt) return; // stale watcher
      if (payload.error) {
        const last = stored.watchStatus;
        const goodReadStands = last && !last.error && last.capturedAt === ticket.capturedAt
          && Date.now() - last.seenAt < GOOD_READ_OUTRANKS_MS;
        if (goodReadStands) return;
        setSession({ watchStatus: { capturedAt: ticket.capturedAt, seenAt: Date.now(), error: payload.error, current: null } });
        return;
      }
      const line = payload.line;
      const moved = ticket.price !== line.price || ticket.points !== line.points || line.offBoard;
      setSession({ watchStatus: { capturedAt: ticket.capturedAt, seenAt: line.seenAt, error: null, current: moved ? line : null } });
    });
  }

  function handleLocated(payload) {
    setSession({ locateResult: payload });
  }

  // page.js's once-per-load check that Unabated's price cells still carry the
  // props a ticket is read from. Stored as posted, including the "checking"
  // state it publishes first, so the panel's banner never outlives the page
  // load it described.
  function handlePageCheck(payload) {
    setSession({ pageCheck: payload });
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

  const handlers = { ticket: handleTicket, error: handleError, watch: handleWatch, ready: handleReady, filters: handleFilters, located: handleLocated, resume_request: offerStoredTicket, pagecheck: handlePageCheck };

  function forwardLocate(locate) {
    if (!locate || typeof locate.at !== "number" || Date.now() - locate.at > LOCATE_MAX_AGE_MS) return;
    window.postMessage({ source: MESSAGE_SOURCE, type: "locate", payload: locate }, window.location.origin);
  }

  function onStorageChanged(changes, area) {
    if (area === "local" && changes.locate && changes.locate.newValue) forwardLocate(changes.locate.newValue);
  }

  function onPageMessage(event) {
    if (event.source !== window) return;
    const data = event.data;
    if (!data || data.source !== MESSAGE_SOURCE || !HANDLED_TYPES.has(data.type)) return;
    const handler = handlers[data.type];
    if (!handler) return;
    try {
      handler(data.payload);
    } catch (error) {
      // The handoff above covers an extension reload that re-injects into this
      // tab; this covers one that cannot (the injection failed, or the tab was
      // in a state Chrome would not script). Retiring stops the dead copy from
      // holding the `message` listener the live one would otherwise never get.
      if (!/context invalidated/i.test(String(error && error.message))) throw error;
      retire("the extension was reloaded; reload this tab");
    }
  }

  // Published for the next injected copy to call (see the handoff at the top).
  // Idempotent, and never throws out: a dead chrome.* handle must not stop the
  // live copy from installing itself.
  function retire(why) {
    if (retired) return;
    retired = true;
    window.removeEventListener("message", onPageMessage);
    try {
      chrome.storage.onChanged.removeListener(onStorageChanged);
    } catch (_error) {
      // Extension context invalidated: the listener is already dead with it.
    }
    // Only if it is still OURS: the newer copy overwrites it right after
    // calling us, and deleting it then would leave the world with no handoff.
    if (window.__unabatedTicketContentRetire === retireForHandoff) delete window.__unabatedTicketContentRetire;
    console.info(`[unabated-ticket] content.js retired (${why})`);
  }

  function retireForHandoff() {
    retire("a newer copy took over");
  }

  window.addEventListener("message", onPageMessage);
  window.__unabatedTicketContentRetire = retireForHandoff;

  // Requests arrive live while this tab is open, or are picked up on load when
  // the panel had to navigate/open the tab (page.js waits for the grid).
  try {
    chrome.storage.onChanged.addListener(onStorageChanged);
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
})();
