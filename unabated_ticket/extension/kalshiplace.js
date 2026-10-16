// Unabated Ticket — the Ticket's Place button on Kalshi (manual_place_plan.md,
// Cal 2026-10-16): a Kalshi line's Ticket gets an amount box filled with the
// stake to act on, which Cal can change, and a Place button; a second click
// sends ONE order through the bets service (POST /place_kalshi.json), which
// signs it with Cal's own Kalshi key on the server, checks the live ask is at
// most 0.1c worse than the Ticket's price and sends an immediate-or-cancel
// limit (bets_service/manual_orders.py). Nothing is ever placed without that
// second click.
//
// Shared by the panel (panel.js) and the phone page (phone.js), so both send
// the same request and read the reply the same way.
//
// Pure: no DOM, no fetch, no chrome.*, no clock. Loaded as a plain <script>
// after kelly.js (exposes globalThis.UnabatedKalshiPlace) and via require() in node.

(function (root) {
  "use strict";

  const inNode = typeof module !== "undefined" && module.exports;
  const kelly = inNode ? require("./kelly.js") : root.UnabatedKelly;

  // Unabated's book id for Kalshi.
  const KALSHI_BOOK_ID = 105;
  // Statuses the bets service answers; the first two put contracts in the account.
  const FILLED_STATUSES = ["filled", "partial"];

  // Does this Ticket get a Place button: a Kalshi line before its game, on the
  // pregame board (a live capture is priced off the in-game screen, which the
  // server's board does not carry).
  function canPlace({ bookId, live, eventId, periodTypeId, betTypeId, sideIndex }) {
    return bookId === KALSHI_BOOK_ID && !live && Number.isInteger(eventId)
      && Number.isInteger(periodTypeId) && Number.isInteger(betTypeId) && (sideIndex === 0 || sideIndex === 1);
  }

  // "$1,234.50" or "250" typed in the box, as dollars, or null.
  function parseStake(text) {
    const cleaned = String(text ?? "").replace(/[$,\s]/g, "");
    if (!/^\d+(\.\d{0,2})?$/.test(cleaned)) return null;
    const value = Number(cleaned);
    return value > 0 ? value : null;
  }

  // The amount box's starting text: the stake to act on, to the cent.
  function stakeText(stake) {
    return typeof stake === "number" && stake > 0 ? stake.toFixed(2) : "";
  }

  // The POST /place_kalshi.json body for one click, or {error} in words.
  //   line         {eventId, periodTypeId, betTypeId, sideIndex, points, price, sourceFormat, sourcePrice}
  //                the Ticket's line at the price it shows
  //   stakeInput   the amount box's text
  //   kalshiReadAt the Kalshi read the Ticket's stake was sized on (the bets
  //                payload's sources.kalshi.fetchedAt)
  //   clickId      a fresh id per click (the server sends each id once)
  function placeRequestOf({ line, stakeInput, kalshiReadAt, clickId }) {
    const stake = parseStake(stakeInput);
    if (stake == null) return { error: `"${stakeInput}" is not an amount in dollars` };
    if (!kalshiReadAt) return { error: "your Kalshi bets have not been read yet, so the stake cannot count them" };
    const ticketAllIn = kelly.bookProbOf({ bookPrice: line.price, sourceFormat: line.sourceFormat, sourcePrice: line.sourcePrice });
    if (!(ticketAllIn > 0 && ticketAllIn < 1)) return { error: "the Ticket has no Kalshi price" };
    return {
      body: {
        clickId, eventId: line.eventId, periodTypeId: line.periodTypeId, betTypeId: line.betTypeId, sideIndex: line.sideIndex,
        points: line.points ?? null, stake, ticketAllIn: Math.round(ticketAllIn * 1e6) / 1e6, kalshiReadAt,
      },
    };
  }

  // How the reply reads: "ok" (contracts bought), "bad" (Kalshi refused or
  // nothing filled), "warn" (no definite answer: check Kalshi), "refused"
  // (nothing sent; the button comes back).
  function replyTone(status) {
    if (FILLED_STATUSES.includes(status)) return "ok";
    if (status === "unknown") return "warn";
    if (status === "refused") return "refused";
    return "bad";
  }

  const api = { KALSHI_BOOK_ID, canPlace, parseStake, stakeText, placeRequestOf, replyTone };

  if (inNode) {
    module.exports = api;
  } else {
    root.UnabatedKalshiPlace = api;
  }
})(typeof globalThis !== "undefined" ? globalThis : this);
