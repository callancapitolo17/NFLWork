// Unabated Ticket — the Ticket's Place button on Kalshi (manual_place_plan.md,
// Cal 2026-10-16), Novig (novig_place_plan.md), Polymarket US (Cal
// 2026-10-21), ProphetX (2026-10-31) and BFA (2026-11-02): a Kalshi, Novig,
// Polymarket US, ProphetX or BFA line's
// Ticket gets an amount box filled with the stake to act on, which Cal can
// change, and a Place button; a second click sends ONE order through the bets
// service (POST /place_kalshi.json, /place_novig.json, /place_polymarket.json or
// /place_prophetx.json or /place_bfa.json),
// which signs it on the server with its key at that venue, checks the
// live price is at most 0.1c worse than the Ticket's and sends an
// immediate-or-cancel limit (bets_service/manual_orders.py, novig_place.py,
// polymarket_place.py, prophetx_place.py) or, at BFA, one straight bet on
// its login (bfa_place.py; whole dollars, main full-game lines). Nothing is
// ever placed without that second click.
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
  // Polymarket US (the CFTC app), not polymarket.com (107): read off the feed's
  // marketSources on 2026-10-22.
  const POLYMARKET_US_BOOK_ID = 111;
  // Prophet Exchange (fillfair.js BOOK_ID_OF_VENUE).
  const PROPHETX_BOOK_ID = 66;
  // feed.PERIODS: 1 is the full game.
  const FULL_GAME_PERIOD_TYPE_ID = 1;
  // Buckeye, which is BFA (teaser.js BUCKEYE_BOOK_ID).
  const BFA_BOOK_ID = 59;
  // Each venue Place sends to, by Unabated book id: its name in words, the
  // bets service route, and the bet types it places (every one: a moneyline
  // the feed names no contract for is found in the venue's own market list
  // by the server runner, Cal 2026-10-19).
  const VENUES = {
    // live: places a Ticket priced off Unabated's in-game fair at a break
    // (Cal, 2026-11-04); the bets service finds the line on the runner's live
    // board and refuses a fair the panel has not seen ready lately.
    [KALSHI_BOOK_ID]: { venue: "kalshi", name: "Kalshi", path: "/place_kalshi.json", betTypeIds: [1, 2, 3], live: true },
    [NOVIG_BOOK_ID]: { venue: "novig", name: "Novig", path: "/place_novig.json", betTypeIds: [1, 2, 3] },
    [POLYMARKET_US_BOOK_ID]: { venue: "polymarket_us", name: "Polymarket US", path: "/place_polymarket.json", betTypeIds: [1, 2, 3] },
    [PROPHETX_BOOK_ID]: { venue: "prophetx", name: "ProphetX", path: "/place_prophetx.json", betTypeIds: [1, 2, 3] },
    // BFA takes whole dollars and books a straight whole or not at all.
    [BFA_BOOK_ID]: { venue: "bfa", name: "BFA", path: "/place_bfa.json", betTypeIds: [1, 2, 3], wholeDollars: true, fullGameOnly: true, mainLinesOnly: true },
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

  // Does this Ticket get a Place button: a line of a Place venue; a live
  // capture (priced off the in-game screen at a break) only at a venue that
  // places live; at BFA a main full-game line only (BFA's board lists no alt
  // or period lines).
  function canPlace({ bookId, live, eventId, periodTypeId, betTypeId, sideIndex, isAlt }) {
    const venue = venueOf(bookId);
    return Boolean(venue) && venue.betTypeIds.includes(betTypeId) && (!live || venue.live === true) && Number.isInteger(eventId)
      && Number.isInteger(periodTypeId) && (!venue.fullGameOnly || periodTypeId === FULL_GAME_PERIOD_TYPE_ID)
      && (!venue.mainLinesOnly || isAlt !== true)
      && (sideIndex === 0 || sideIndex === 1);
  }

  // "$1,234.50" or "250" typed in the box, as dollars, or null.
  function parseStake(text) {
    const cleaned = String(text ?? "").replace(/[$,\s]/g, "");
    if (!/^\d+(\.\d{0,2})?$/.test(cleaned)) return null;
    const value = Number(cleaned);
    return value > 0 ? value : null;
  }

  // The amount box's starting text: the stake to act on, to the cent, or in
  // whole dollars (rounded down) at a venue that takes no cents.
  function stakeText(stake, bookId) {
    if (!(typeof stake === "number" && stake > 0)) return "";
    const venue = venueOf(bookId);
    if (venue && venue.wholeDollars) return Math.floor(stake) >= 1 ? String(Math.floor(stake)) : "";
    return stake.toFixed(2);
  }

  // {kalshi, novig, polymarket_us, prophetx, bfa}: the read times a stake counted, off the bets payload's
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
    return {
      kalshi: fetchedAtOf("kalshi"), novig, polymarket_us: fetchedAtOf("polymarket_us"), prophetx: fetchedAtOf("prophetx"),
      bfa: fetchedAtOf("bfa"),
    };
  }

  // The request for one click, {path, body}, or {error} in words.
  //   line      {bookId, eventId, periodTypeId, betTypeId, sideIndex, points, price, sourceFormat, sourcePrice}
  //             the Ticket's line at the price it shows
  //   stakeInput  the amount box's text
  //   readAts   {kalshi, novig, polymarket_us, prophetx, bfa} (readAtsOf): the reads the Ticket's stake
  //             was sized on (null: never read)
  //   clickId   a fresh id per click (the server sends each id once)
  //   liveFairReadAtMs  a live Ticket only: when the Unabated tab last said
  //             its fair was ready (live.liveFairReadAtMs), null when it does not now
  // An order needs its own venue's read and sends the other venues' too, or
  // null when a venue was never read (the server checks for a Place at any
  // venue the stake has not seen).
  function placeRequestOf({ line, stakeInput, readAts, clickId, liveFairReadAtMs = null }) {
    const venue = venueOf(line.bookId);
    if (!venue) return { error: "this book has no Place button" };
    if (line.live && venue.live !== true) return { error: `${venue.name} does not place live yet` };
    if (line.live && typeof liveFairReadAtMs !== "number") {
      return { error: "the in-game fair is not ready on the Unabated tab now (the break ended, or the tab stopped reading)" };
    }
    const stake = parseStake(stakeInput);
    if (stake == null) return { error: `"${stakeInput}" is not an amount in dollars` };
    if (venue.wholeDollars && !Number.isInteger(stake)) return { error: `${venue.name} takes whole dollars, not ${stakeInput}` };
    const venueReadAt = readAts ? readAts[venue.venue] : null;
    if (!venueReadAt) return { error: `your ${venue.name} bets have not been read yet, so the stake cannot count them` };
    const ticketAllIn = kelly.bookProbOf({ bookPrice: line.price, sourceFormat: line.sourceFormat, sourcePrice: line.sourcePrice });
    if (!(ticketAllIn > 0 && ticketAllIn < 1)) return { error: `the Ticket has no ${venue.name} price` };
    const body = {
      clickId, eventId: line.eventId, periodTypeId: line.periodTypeId, betTypeId: line.betTypeId, sideIndex: line.sideIndex,
      points: line.points ?? null, stake, ticketAllIn: Math.round(ticketAllIn * 1e6) / 1e6,
    };
    return {
      path: venue.path,
      body: {
        ...body, kalshiReadAt: readAts.kalshi || null, novigReadAt: readAts.novig || null,
        polymarketReadAt: readAts.polymarket_us || null, prophetxReadAt: readAts.prophetx || null,
        bfaReadAt: readAts.bfa || null,
        ...(line.live ? { live: true, liveFairReadAt: new Date(liveFairReadAtMs).toISOString() } : {}),
      },
    };
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

  const api = { KALSHI_BOOK_ID, NOVIG_BOOK_ID, POLYMARKET_US_BOOK_ID, PROPHETX_BOOK_ID, BFA_BOOK_ID, venueOf, canPlace, parseStake, stakeText, readAtsOf, placeRequestOf, replyTone };

  if (inNode) {
    module.exports = api;
  } else {
    root.UnabatedPlace = api;
  }
})(typeof globalThis !== "undefined" ? globalThis : this);
