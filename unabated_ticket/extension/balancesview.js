// The Balances strip on the Edges tab, shared by the panel (panel.js) and the
// phone page (server/phone/phone.js) so both read the same. No fetch, no
// chrome.*; viewOf is pure, fillStrip touches only the <details> it is handed.
// Loaded as a plain <script> (exposes globalThis.UnabatedBalancesView) and via
// require() in tests/balancesview.test.js. Nothing here writes anywhere.
//
// Input: the bets service's GET /balances.json body (bets_service/balances.py):
//   {generatedAt, balances: {venue: {cash, inBets, available, creditLimit, fetchedAt, error, failedAt?}}}
// with dollars or null; a venue never read has cash null.
// Output: viewOf() -> {loaded, summary, rows, totalAvailableText, notes}, the
// strip's text; fillStrip(details, view) lays it out (both pages' markup:
// <details class="balances-strip"><summary></summary><div class="balances-body"></div></details>). Cal 2026-10-31: Kalshi, Novig and Polymarket US with their
// sub-accounts only; collapsed it shows the three main accounts.

(function (root) {
  "use strict";

  // Row order and words. A venue the service sends that is not listed here is
  // not shown (Cal picked these books; a new one is added here on purpose).
  const BALANCE_ROWS = [
    { venue: "kalshi", label: "Kalshi", summaryLabel: "Kalshi", isSubaccount: false },
    { venue: "kalshi_bot", label: "Kalshi · Bot", summaryLabel: null, isSubaccount: true },
    { venue: "novig", label: "Novig", summaryLabel: "Novig", isSubaccount: false },
    { venue: "novig_place", label: "Novig · Place", summaryLabel: null, isSubaccount: true },
    { venue: "novig_bot", label: "Novig · Bot", summaryLabel: null, isSubaccount: true },
    { venue: "polymarket_us", label: "Polymarket US", summaryLabel: "Poly", isSubaccount: false },
  ];
  // A reading older than this, or one whose last read failed, is flagged.
  const STALE_AFTER_MS = 5 * 60 * 1000;

  function fmtDollars(value, wholeDollars) {
    const digits = wholeDollars ? 0 : 2;
    const text = Math.abs(value).toLocaleString("en-US", {
      style: "currency", currency: "USD", minimumFractionDigits: digits, maximumFractionDigits: digits,
    });
    return value < 0 ? `-${text}` : text;
  }

  function fmtAge(fetchedAt, nowMs) {
    const minutes = Math.max(0, Math.round((nowMs - Date.parse(fetchedAt)) / 60000));
    if (minutes < 60) return `${minutes} min ago`;
    return `${Math.round(minutes / 60)} h ago`;
  }

  /** What a book can bet now: `available` when the service reports it, else cash. */
  function availableOf(reading) {
    return reading.available ?? reading.cash;
  }

  function rowOf(spec, reading, nowMs) {
    if (!reading || reading.cash === null || reading.cash === undefined) {
      return { ...spec, read: false, available: null, inBets: null, availableText: "unavailable", inBetsText: "",
        why: reading && reading.error ? reading.error : "no reading yet" };
    }
    const available = availableOf(reading);
    const ageMs = reading.fetchedAt ? nowMs - Date.parse(reading.fetchedAt) : Infinity;
    const stale = Boolean(reading.error) || ageMs > STALE_AFTER_MS;
    return {
      ...spec, read: true, available, inBets: reading.inBets,
      availableText: fmtDollars(available, false),
      inBetsText: reading.inBets === null ? "—" : fmtDollars(reading.inBets, false),
      stale, why: reading.error || null,
      ageText: reading.fetchedAt ? fmtAge(reading.fetchedAt, nowMs) : "never",
    };
  }

  /**
   * The strip's text from a /balances.json body (null while none has loaded).
   * `error` is the last fetch's failure message, or null.
   */
  function viewOf(payload, error, nowMs) {
    const balances = (payload && payload.balances) || {};
    const rows = BALANCE_ROWS
      .filter((spec) => spec.venue in balances)
      .map((spec) => rowOf(spec, balances[spec.venue], nowMs));
    const read = rows.filter((row) => row.read);
    const summary = rows
      .filter((row) => row.summaryLabel !== null)
      .map((row) => ({ label: row.summaryLabel, text: row.read ? fmtDollars(row.available, true) : "—", stale: !row.read || row.stale }));
    const notes = [];
    for (const row of rows) {
      if (!row.read) notes.push(`${row.label} unavailable: ${row.why}`);
      else if (row.stale) notes.push(`${row.label} read ${row.ageText}${row.why ? " (last read failed)" : ""}`);
    }
    if (read.length < rows.length) notes.push("The total leaves out books shown as unavailable.");
    if (error) notes.unshift(`Could not load balances: ${error}`);
    return {
      loaded: payload !== null,
      summary,
      rows,
      totalAvailableText: fmtDollars(read.reduce((total, row) => total + row.available, 0), true),
      notes,
    };
  }

  function make(doc, tag, className, text) {
    const node = doc.createElement(tag);
    if (className) node.className = className;
    if (text != null) node.textContent = text;
    return node;
  }

  /**
   * Lays `view` out in the strip. Collapsed: the main accounts' available
   * money; open: every account with available and in bets, then the notes.
   * The <details> keeps its open state across refreshes. Hidden until the
   * first /balances.json has loaded or failed.
   */
  function fillStrip(details, view) {
    const doc = details.ownerDocument;
    const summary = details.querySelector("summary");
    const body = details.querySelector(".balances-body");
    details.hidden = !view.loaded && view.notes.length === 0;
    summary.replaceChildren(make(doc, "span", "balances-label", "Balances"));
    for (const item of view.summary) {
      const part = make(doc, "span", item.stale ? "balances-item stale" : "balances-item", `${item.label} `);
      part.append(make(doc, "b", null, item.text));
      summary.append(part);
    }
    if (view.notes.length) summary.append(make(doc, "span", "balances-warn", "!"));

    const rows = view.rows.map((row) => {
      const tr = make(doc, "tr", row.isSubaccount ? "sub" : null);
      tr.append(make(doc, "td", null, row.label));
      tr.append(make(doc, "td", row.read ? (row.stale ? "num stale" : "num") : "num stale", row.availableText));
      tr.append(make(doc, "td", "num muted", row.inBetsText));
      if (row.why) tr.title = row.why;
      return tr;
    });
    const head = make(doc, "tr");
    for (const text of ["Book", "Available", "In bets"]) head.append(make(doc, "th", null, text));
    const total = make(doc, "tr", "total");
    total.append(make(doc, "td", null, "Total"), make(doc, "td", "num", view.totalAvailableText), make(doc, "td"));
    const table = make(doc, "table", "balances-table");
    table.append(head, ...rows, total);
    const notes = view.notes.map((text) => make(doc, "div", "balances-note", text));
    body.replaceChildren(table, ...notes);
  }

  const api = { BALANCE_ROWS, STALE_AFTER_MS, viewOf, fillStrip, fmtDollars };
  if (typeof module !== "undefined" && module.exports) module.exports = api;
  else root.UnabatedBalancesView = api;
})(typeof globalThis !== "undefined" ? globalThis : this);
