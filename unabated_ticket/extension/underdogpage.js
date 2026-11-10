// Underdog tab, MAIN world, document_start (2026-11-09): keeps the read headers
// of the newest request the Underdog web app sends to its API, so the
// extension's read (underdog.js readAccount) rides the app's own login and
// never refreshes it (the refresh token rotates; reusing it would log the app
// out). The app sends its API calls by XMLHttpRequest; the hook only watches,
// it never changes a request.
//
// Exposes one function, globalThis.__unabatedUnderdogRead(), which background.js
// calls through chrome.scripting: {push} | {error, stale}. The headers stay in
// this closure; the page already holds the same token.
//
// Side effects: wraps XMLHttpRequest.prototype.open / setRequestHeader / send.

(function () {
  "use strict";

  const lib = globalThis.UnabatedUnderdog;
  if (!lib || globalThis.__unabatedUnderdogRead) return;
  let capture = null;

  const originalOpen = XMLHttpRequest.prototype.open;
  const originalSetRequestHeader = XMLHttpRequest.prototype.setRequestHeader;
  const originalSend = XMLHttpRequest.prototype.send;
  const sentHeaders = new WeakMap();

  XMLHttpRequest.prototype.open = function (method, url, ...rest) {
    sentHeaders.set(this, { url: String(url), headers: {} });
    return originalOpen.call(this, method, url, ...rest);
  };
  XMLHttpRequest.prototype.setRequestHeader = function (name, value) {
    const sent = sentHeaders.get(this);
    if (sent) sent.headers[name] = value;
    return originalSetRequestHeader.call(this, name, value);
  };
  XMLHttpRequest.prototype.send = function (...args) {
    try {
      const sent = sentHeaders.get(this);
      if (sent && sent.url.startsWith(lib.API_ORIGIN)) {
        const headers = lib.readHeadersOf(sent.headers);
        if (headers) capture = { headers, expiresAtMs: lib.tokenExpiryMs(headers.Authorization) };
      }
    } catch (error) {
      console.warn("[unabated-ticket] underdog: could not read the app's request headers", error);
    }
    return originalSend.apply(this, args);
  };

  Object.defineProperty(globalThis, "__unabatedUnderdogRead", {
    value: async function () {
      const nowMs = Date.now();
      if (!lib.isFresh(capture, nowMs)) return { error: lib.STALE, stale: true };
      try {
        return await lib.readAccount({ fetchImpl: fetch.bind(globalThis), headers: capture.headers, nowMs });
      } catch (error) {
        const message = error && error.name === "TimeoutError"
          ? `Underdog did not answer within ${lib.FETCH_TIMEOUT_MS / 1000}s` : String((error && error.message) || error);
        return { error: message, stale: false };
      }
    },
    writable: false, configurable: false, enumerable: false,
  });
})();
