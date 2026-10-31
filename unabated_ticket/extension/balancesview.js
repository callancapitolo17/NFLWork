// The Balances line on the Edges tab, shared by the panel (panel.js) and the
// phone page (server/phone/phone.js) so both read the same. No fetch, no
// chrome.*; viewOf is pure, fillLine touches only the element it is handed.
// Loaded as a plain <script> (exposes globalThis.UnabatedBalancesView) and via
// require() in tests/balancesview.test.js. Nothing here writes anywhere.
//
// Input: the bets service's GET /balances.json body (bets_service/balances.py):
//   {generatedAt, balances: {venue: {cash, inBets, available, creditLimit, fetchedAt, error, failedAt?}}}
// with dollars or null; a venue never read has cash null.
// Output: viewOf() -> {loaded, items: [{label, text, stale, why}], error}; fillLine
// lays it out as one line: "Available  Kalshi $4,872  Novig $24  Place $1,018  Poly $980".
// Cal 2026-10-31: money he can bet now only (no in-bets, no total, no bot
// accounts); Novig's Place sub-account as its own item, labelled just "Place".

(function (root) {
  "use strict";

  // The line's items, in order. A venue the service sends that is not listed
  // here is not shown (the bots' wallets live on the tracker's Balances tab).
  const BALANCE_ITEMS = [
    { venue: "kalshi", label: "Kalshi" },
    { venue: "novig", label: "Novig" },
    { venue: "novig_place", label: "Place" },
    { venue: "polymarket_us", label: "Poly" },
  ];
  // A reading older than this, or one whose last read failed, is amber.
  const STALE_AFTER_MS = 5 * 60 * 1000;

  function fmtWholeDollars(value) {
    const text = Math.abs(value).toLocaleString("en-US", {
      style: "currency", currency: "USD", minimumFractionDigits: 0, maximumFractionDigits: 0,
    });
    return value < 0 ? `-${text}` : text;
  }

  function fmtAge(fetchedAt, nowMs) {
    const minutes = Math.max(0, Math.round((nowMs - Date.parse(fetchedAt)) / 60000));
    return minutes < 60 ? `${minutes} min ago` : `${Math.round(minutes / 60)} h ago`;
  }

  /** One book's item; `why` (the tooltip) says what is wrong with an amber one. */
  function itemOf(spec, reading, nowMs) {
    if (!reading || reading.cash === null || reading.cash === undefined) {
      return { label: spec.label, text: "—", stale: true, why: reading && reading.error ? reading.error : "no reading yet" };
    }
    const available = reading.available ?? reading.cash;
    const ageMs = reading.fetchedAt ? nowMs - Date.parse(reading.fetchedAt) : Infinity;
    const age = reading.fetchedAt ? fmtAge(reading.fetchedAt, nowMs) : "never";
    let why = null;
    if (reading.error) why = `last read failed (${reading.error}); this is from ${age}`;
    else if (ageMs > STALE_AFTER_MS) why = `read ${age}`;
    return { label: spec.label, text: fmtWholeDollars(available), stale: why !== null, why };
  }

  /**
   * The line's text from a /balances.json body (null while none has loaded).
   * `error` is the last fetch's failure message, or null.
   */
  function viewOf(payload, error, nowMs) {
    const balances = (payload && payload.balances) || {};
    const items = BALANCE_ITEMS
      .filter((spec) => spec.venue in balances)
      .map((spec) => itemOf(spec, balances[spec.venue], nowMs));
    return { loaded: payload !== null, items, error: error || null };
  }

  function make(doc, tag, className, text) {
    const node = doc.createElement(tag);
    if (className) node.className = className;
    if (text != null) node.textContent = text;
    return node;
  }

  /**
   * Lays `view` out in the line element. Hidden until the first
   * /balances.json has loaded or failed; a failed fetch keeps the last
   * numbers and turns the label amber with the reason as its tooltip.
   */
  function fillLine(line, view) {
    const doc = line.ownerDocument;
    line.hidden = !view.loaded && view.error === null;
    const label = make(doc, "span", view.error ? "balances-label stale" : "balances-label", "Available");
    if (view.error) label.title = `Could not load balances: ${view.error}`;
    const items = view.items.map((item) => {
      const part = make(doc, "span", item.stale ? "balances-item stale" : "balances-item", `${item.label} `);
      part.append(make(doc, "b", null, item.text));
      if (item.why) part.title = item.why;
      return part;
    });
    line.replaceChildren(label, ...items);
  }

  const api = { BALANCE_ITEMS, STALE_AFTER_MS, viewOf, fillLine };
  if (typeof module !== "undefined" && module.exports) module.exports = api;
  else root.UnabatedBalancesView = api;
})(typeof globalThis !== "undefined" ? globalThis : this);
