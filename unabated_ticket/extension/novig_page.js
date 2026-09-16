// Unabated Ticket — Novig main-world mirror (#116).
//
// Runs in the Novig web app's page context (world: MAIN, document_start, so
// the wrapper is in place before the app bundle captures window.fetch) and
// mirrors the responses the app fetches FOR ITSELF for its Portfolio screen:
// the Apollo queries ActivePortfolioOrders_Query, SettledPortfolioOrders_Query
// and ParlayPortfolioQuery, POSTed to https://api.novig.us/v1/graphql. Nothing
// here sends a request, reads a token, or touches the DOM: a watched response
// is cloned, parsed and handed to novig_content.js via window.postMessage.
// Subscriptions run over a WebSocket link and are not mirrored — the app
// refetches its lists over fetch after a fill, and that refetch is caught.
//
// Side effects: replaces window.fetch with a wrapper that forwards every call
// unchanged; posts {source: "unabated-ticket-novig", type: "portfolio"} messages.

(function () {
  "use strict";

  if (window.__unabatedTicketNovigActive) return;
  window.__unabatedTicketNovigActive = true;

  const MESSAGE_SOURCE = "unabated-ticket-novig";
  const GRAPHQL_URL_PREFIX = "https://api.novig.us/v1/graphql";
  const WATCHED = new Set(["ActivePortfolioOrders_Query", "SettledPortfolioOrders_Query", "ParlayPortfolioQuery"]);

  function post(type, payload) {
    window.postMessage({ source: MESSAGE_SOURCE, type, payload }, window.location.origin);
  }

  function requestUrl(input) {
    if (typeof input === "string") return input;
    if (input && typeof input.url === "string") return input.url;
    return String(input);
  }

  // The operations of one GraphQL POST body (Apollo may batch into an array).
  function operationsOf(bodyText) {
    let body;
    try {
      body = JSON.parse(bodyText);
    } catch (_error) {
      return [];
    }
    const items = Array.isArray(body) ? body : [body];
    return items.filter((item) => item && WATCHED.has(item.operationName));
  }

  async function mirror(bodyText, requestCopy, responseCopy) {
    const text = bodyText != null ? bodyText : requestCopy ? await requestCopy.text() : null;
    if (!text) return;
    const operations = operationsOf(text);
    if (!operations.length) return;
    const json = await responseCopy.json();
    const results = Array.isArray(json) ? json : [json];
    const at = new Date().toISOString();
    operations.forEach((operation, index) => {
      const result = results[Math.min(index, results.length - 1)] || {};
      const error = result.errors && result.errors.length ? String(result.errors[0].message || "GraphQL error") : null;
      post("portfolio", { operationName: operation.operationName, variables: operation.variables || {}, data: error ? null : result.data || null, error, at, url: window.location.href });
    });
  }

  // Both bodies are copied BEFORE anyone else can read them: the request's
  // before the real fetch consumes it, the response's synchronously in the
  // first .then so the app's own continuation cannot drain it first.
  const originalFetch = window.fetch;
  window.fetch = function unabatedTicketFetch(input, init) {
    let bodyText = null;
    let requestCopy = null;
    let watched;
    try {
      watched = requestUrl(input).startsWith(GRAPHQL_URL_PREFIX);
      if (watched) {
        if (init && typeof init.body === "string") bodyText = init.body;
        else if (input && typeof input.clone === "function") requestCopy = input.clone();
      }
    } catch (_error) {
      watched = false;
    }
    const pending = originalFetch.apply(this, arguments);
    if (watched) {
      pending.then((response) => {
        let responseCopy;
        try {
          responseCopy = response.clone();
        } catch (error) {
          console.info("[unabated-ticket] novig mirror could not clone a response:", error && error.message);
          return;
        }
        mirror(bodyText, requestCopy, responseCopy).catch((error) => console.info("[unabated-ticket] novig mirror skipped a response:", error && error.message));
      }, () => {});
    }
    return pending;
  };

  console.info("[unabated-ticket] novig_page.js active (mirroring portfolio queries)");
})();
