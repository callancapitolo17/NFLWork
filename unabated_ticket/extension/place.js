// Unabated Ticket — the Ticket's Place button on Kalshi (manual_place_plan.md,
// Cal 2026-10-16) and Novig (novig_place_plan.md): a Kalshi or Novig line's
// Ticket gets an amount box filled with the stake to act on, which Cal can
// change, and a Place button; a second click sends ONE order through the bets
// service (POST /place_kalshi.json or /place_novig.json), which signs it on
// the server with Place's own key at that venue, checks the live price is at
// most 0.1c worse than the Ticket's and sends an immediate-or-cancel limit
// (bets_service/manual_orders.py, novig_place.py). Nothing is ever placed
// without that second click.
//
// Shared by the panel (panel.js) and the phone page (phone.js), so both send
// the same request and read the reply the same way.
//
// Pure: no DOM, no fetch, no chrome.*, no clock. Loaded as a plain <script>
// after kelly.js (exposes globalThis.UnabatedPlace) and via require() in node.

(function (root) {
  "use strict";

  const inNode = typeof module !== "undefined" && module.exports;
  const kelly = inNode ? require("./kelly.js") : root.UnabatedKelly;

  // Unabated's book ids.
  const KALSHI_BOOK_ID = 105;
  const NOVIG_BOOK_ID = 89;
  // Each venue Place sends to, by Unabated book id: its name in words, the
  // bets service route, and the bet types it places (Unabated carries a
  // Novig outcome id on spreads and totals only, never on a moneyline).
  const VENUES = {
    [KALSHI_BOOK_ID]: { venue: "kalshi", name: "Kalshi", path: "/place_kalshi.json", betTypeIds: [1, 2, 3] },
    [NOVIG_BOOK_ID]: { venue: "novig", name: "Novig", path: "/place_novig.json", betTypeIds: [2, 3] },
  };
  // The bets service's Novig reads a stake counts: Cal's app account and
  // Place's own sub-account, both required, and the bot's sub-account when
  // the service reads it.
  const NOVIG_REQUIRED_SOURCES = ["novig", "novig_place"];
  const NOVIG_OPTIONAL_SOURCES = ["novig_bot"];
  // Statuses the bets service answers; the first two put contracts in the account.
  const FILLED_STATUSES = ["filled", "partial"];

  // The venue a book id places at ({venue, name, path, betTypeIds}), or null.
  function venueOf(bookId) {
    return VENUES[bookId] || null;
  }

  // Does this Ticket get a Place button: a Kalshi or Novig line before its
  // game, on the pregame board (a live capture is priced off the in-game
  // screen, which the server's board does not carry).
  function canPlace({ bookId, live, eventId, periodTypeId, betTypeId, sideIndex }) {
    const venue = venueOf(bookId);
    return Boolean(venue) && venue.betTypeIds.includes(betTypeId) && !live && Number.isInteger(eventId)
      && Number.isInteger(periodTypeId) && (sideIndex === 0 || sideIndex === 1);
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

  // {kalshi, novig}: the read times a stake counted, off the bets payload's
  // sources ({name: {fetchedAt}}). Novig's is the OLDEST of its reads (a
  // Place fill after any of them may be missing from the stake), null while
  // either required read has never landed.
  function readAtsOf(sources) {
    const fetchedAtOf = (name) => (sources && sources[name] ? sources[name].fetchedAt || null : null);
    const required = NOVIG_REQUIRED_SOURCES.map(fetchedAtOf);
    const optional = NOVIG_OPTIONAL_SOURCES.filter((name) => sources && sources[name]).map(fetchedAtOf);
    const novigReads = [...required, ...optional];
    const novig = novigReads.every(Boolean)
      ? novigReads.reduce((oldest, readAt) => (Date.parse(readAt) < Date.parse(oldest) ? readAt : oldest))
      : null;
    return { kalshi: fetchedAtOf("kalshi"), novig };
  }

  // The request for one click, {path, body}, or {error} in words.
  //   line      {bookId, eventId, periodTypeId, betTypeId, sideIndex, points, price, sourceFormat, sourcePrice}
  //             the Ticket's line at the price it shows
  //   stakeInput  the amount box's text
  //   readAts   {kalshi, novig} (readAtsOf): the reads the Ticket's stake
  //             was sized on (null: never read)
  //   clickId   a fresh id per click (the server sends each id once)
  // A Kalshi order needs the Kalshi read; a Novig order needs the Novig read
  // and sends the Kalshi one too (the server checks for a Place at either
  // venue the stake has not seen).
  function placeRequestOf({ line, stakeInput, readAts, clickId }) {
    const venue = venueOf(line.bookId);
    if (!venue) return { error: "this book has no Place button" };
    const stake = parseStake(stakeInput);
    if (stake == null) return { error: `"${stakeInput}" is not an amount in dollars` };
    const venueReadAt = readAts ? readAts[venue.venue] : null;
    if (!venueReadAt) return { error: `your ${venue.name} bets have not been read yet, so the stake cannot count them` };
    const ticketAllIn = kelly.bookProbOf({ bookPrice: line.price, sourceFormat: line.sourceFormat, sourcePrice: line.sourcePrice });
    if (!(ticketAllIn > 0 && ticketAllIn < 1)) return { error: `the Ticket has no ${venue.name} price` };
    const body = {
      clickId, eventId: line.eventId, periodTypeId: line.periodTypeId, betTypeId: line.betTypeId, sideIndex: line.sideIndex,
      points: line.points ?? null, stake, ticketAllIn: Math.round(ticketAllIn * 1e6) / 1e6,
    };
    if (venue.venue === "kalshi") return { path: venue.path, body: { ...body, kalshiReadAt: venueReadAt } };
    return { path: venue.path, body: { ...body, novigReadAt: venueReadAt, kalshiReadAt: readAts.kalshi || null } };
  }

  // How the reply reads: "ok" (contracts bought), "bad" (the venue refused or
  // nothing filled), "warn" (no definite answer: check the venue), "refused"
  // (nothing sent; the button comes back).
  function replyTone(status) {
    if (FILLED_STATUSES.includes(status)) return "ok";
    if (status === "unknown") return "warn";
    if (status === "refused") return "refused";
    return "bad";
  }

  const api = { KALSHI_BOOK_ID, NOVIG_BOOK_ID, venueOf, canPlace, parseStake, stakeText, readAtsOf, placeRequestOf, replyTone };

  if (inNode) {
    module.exports = api;
  } else {
    root.UnabatedPlace = api;
  }
})(typeof globalThis !== "undefined" ? globalThis : this);
