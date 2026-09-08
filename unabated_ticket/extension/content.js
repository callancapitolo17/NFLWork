// Unabated Ticket — isolated-world bridge.
//
// page.js (MAIN world) cannot use chrome.* APIs, so it posts window messages;
// this script forwards them to background.js, which owns storage writes.
// Side effects: none on the page.

(function () {
  "use strict";

  const MESSAGE_SOURCE = "unabated-ticket";
  const FORWARDED_TYPES = new Set(["ticket", "watch", "error", "ready"]);

  function forward(message) {
    try {
      chrome.runtime.sendMessage(message, () => {
        // Reading lastError marks it handled; the extension may have been reloaded.
        if (chrome.runtime.lastError) console.info("[unabated-ticket] forward failed:", chrome.runtime.lastError.message);
      });
    } catch (_error) {
      // Extension context invalidated (reloaded while the page stayed open). Nothing to do.
    }
  }

  console.info("[unabated-ticket] content.js active");
  window.addEventListener("message", (event) => {
    if (event.source !== window) return;
    const data = event.data;
    if (!data || data.source !== MESSAGE_SOURCE || !FORWARDED_TYPES.has(data.type)) return;
    forward({ type: data.type, payload: data.payload });
  });
})();
