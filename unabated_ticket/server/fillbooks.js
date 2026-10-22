// Other books at the fill (2026-10-13, Cal's picks 1A 2B 4B): when a new
// straight bet shows up, what the rest of the market offered at the moment it
// was placed — the best other book's price at the bet's number and how many
// books beat the price taken. Stored only; nothing reads it back yet.
//
// Read off the scanner's line history (edgemove.js), the way fillfair.js
// reads the fill-time fair: every book's line at the bet's number (same
// game, period, bet type, side and number), its newest observation at or
// before placedAt. Only while the scanner watched through the fill — a bet
// placed before it was watching, or while its league's snapshot was not
// loading, gets nothing. Never a backfill.
//
// Counted the Edges tab's way (pricecheck.js): every book on Unabated, by
// price as Unabated shows it (an exchange's is already all-in of its fee),
// one book under two country names once at a shared price; a Kalshi fill is
// compared with its taker fee added, as Unabated shows Kalshi. The bet's own
// venue's line is left out of the other books. A book counts only while its
// line is still on the board, unblurred, at a live book, and had been seen by
// the fill; a book whose number moved since the fill is left out (the history
// restarts on a new number).
//
// Pure: no fetch, no clock, no timers. Required by server/runner.js and
// tests/fillbooks.test.js.

"use strict";

const feed = require("../extension/feed.js");
const fillfair = require("../extension/fillfair.js");
const pricecheck = require("../extension/pricecheck.js");
const bets = require("../extension/bets.js");
const edgemove = require("../extension/edgemove.js");

const STATUS_ON_BOARD = 1;
// A placement up to this far past our clock is a venue clock running ahead;
// it is looked at again next pass (fillfair.js uses the same grace).
const CLOCK_SKEW_GRACE_MS = 60 * 1000;
// Kalshi's taker fee per contract is 7% x P x (1 - P). Unabated shows Kalshi
// all-in (ask + fee, matched to six decimals 2026-09-21), so a Kalshi fill is
// put on that basis before it is compared with the board.
const KALSHI_TAKER_FEE_RATE = 0.07;

// Why a bet gets no row. Every one is final: the history never gains back
// the moment of a fill.
const REFUSED = Object.freeze({
  notStraight: "a parlay or teaser leg has no price of its own",
  noBoardMarket: "the board lists no market of its bet type and period",
  noPlacedAt: "no placed time on the record",
  noPrice: "no American price on the record",
  future: "placed after now (a clock or time-zone error)",
  notWatching: "placed before the scanner was watching",
  leagueGap: "placed while its league's snapshot was not loading",
  tooLate: "placed longer ago than the line history keeps every change",
});

// The probability the bet was filled at, on Unabated's basis: an exchange's
// exact fill when the record keeps it, else its American price; a Kalshi fill
// with its taker fee added.
function fillProbOf(bet) {
  const exact = fillfair.fillPriceOf(bet);
  const prob = exact && exact.sourcePrice != null ? exact.sourcePrice : pricecheck.impliedProb(bet.price);
  if (prob == null || bet.venue !== "kalshi") return prob;
  return prob + KALSHI_TAKER_FEE_RATE * prob * (1 - prob);
}

// A probability as an American price, to two decimals.
function probToAmerican(prob) {
  const american = prob >= 0.5 ? -100 * prob / (1 - prob) : 100 * (1 - prob) / prob;
  return Math.round(american * 100) / 100;
}

function refusalOf(bet, placedMs, observingSince, now) {
  if (bet.isParlayLeg === true) return REFUSED.notStraight;
  if (!fillfair.onBoardMarket(bet)) return REFUSED.noBoardMarket;
  if (!Number.isFinite(placedMs)) return REFUSED.noPlacedAt;
  if (fillProbOf(bet) == null) return REFUSED.noPrice;
  if (placedMs > now + CLOCK_SKEW_GRACE_MS) return REFUSED.future;
  if (placedMs < observingSince) return REFUSED.notWatching;
  // The history keeps a line's changes for EDGE_MOVE_WINDOW_MS (edgemove.js
  // prune). Past that, a line that moved since the fill no longer holds its
  // fill-time price and would silently drop out of the count; refuse the bet
  // rather than store a partial one.
  if (now - placedMs > edgemove.EDGE_MOVE_WINDOW_MS) return REFUSED.tooLate;
  return null;
}

// The newest observation at or before placedMs, or null when the line was
// first seen after it.
function entryAtFill(entries, placedMs) {
  if (!Array.isArray(entries) || entries.length === 0 || entries[0].at > placedMs) return null;
  let entry = entries[0];
  for (const candidate of entries) if (candidate.at <= placedMs) entry = candidate;
  return entry;
}

// The bet's own venue on the board: by Unabated's book id where fillfair.js
// knows it, else by name (DraftKings, Polymarket US), so the price taken
// never reads as another book's.
function isOwnBook(line, book, venue) {
  if (line.bookId === fillfair.BOOK_ID_OF_VENUE[venue]) return true;
  return String(book.name || "").toLowerCase() === bets.venueLabel(venue).toLowerCase();
}

// Can this book's line stand as another book's price right now.
function bookLineCounts(line, state, venue) {
  if (!line || line.bookId === feed.UNABATED_LINE_BOOK_ID) return false;
  if (line.statusId !== STATUS_ON_BOARD || line.isBlurred) return false;
  const book = state.books[line.bookId];
  return Boolean(book && book.isLive) && !isOwnBook(line, book, venue);
}

// The row for one bet, from its same_line rows: {row} or {reason}.
function decide(bet, placedMs, rows, state, history, leagueObservingSince) {
  const leagueId = rows[0].leagueId;
  const leagueSince = leagueObservingSince[leagueId];
  if (typeof leagueSince !== "number" || placedMs < leagueSince) return { reason: REFUSED.leagueGap };
  const prices = [];
  for (const row of rows) {
    if (!bookLineCounts(state.lines[row.key], state, bet.venue)) continue;
    const entry = entryAtFill(history[row.key], placedMs);
    if (!entry || pricecheck.impliedProb(entry.price) == null) continue;
    prices.push({ bookId: row.book.id, bookName: row.book.name, price: entry.price });
  }
  const others = pricecheck.mergeSameNameBooks(prices, null)
    .sort((a, b) => pricecheck.impliedProb(a.price) - pricecheck.impliedProb(b.price));
  const fillProb = fillProbOf(bet);
  const best = others[0] || null;
  return { row: {
    betId: bet.id,
    placedAt: bet.placedAt,
    fillPrice: probToAmerican(fillProb),
    bookCount: others.length,
    booksBetter: others.filter((entry) => pricecheck.impliedProb(entry.price) < fillProb).length,
    bestBook: best ? best.bookName : null,
    bestPrice: best ? best.price : null,
  } };
}

// One capture pass: {saves, refusals}. A save is the row POST
// /fill_books.json takes — {betId, placedAt, fillPrice, bookCount,
// booksBetter, bestBook, bestPrice}; a refusal {betId, reason} is final. A
// bet whose game the board has not resolved yet is in neither and is looked
// at again next pass.
//   records         bet records (bets.js contract); only open ones are read
//   skipIds         Set of bet ids already saved or decided
//   state           the scanner's feed state {lines, events, books, teams}
//   boardLines      one describeLine row per event (runner boardLines())
//   history         edgemove history: line key -> entries
//   observingSince  ms since when the history has no gap, or null
//   leagueObservingSince  leagueId -> ms (scanner status)
//   now             ms
function captureFillBooks({ records, skipIds, state, boardLines, history, observingSince, leagueObservingSince, now }) {
  const result = { saves: [], refusals: [] };
  if (!state || typeof observingSince !== "number") return result;
  const pending = [];
  for (const bet of records || []) {
    if (bet.status !== "open" || bet.unmatchable || skipIds.has(bet.id)) continue;
    const placedMs = Date.parse(bet.placedAt);
    const reason = refusalOf(bet, placedMs, observingSince, now);
    if (reason) result.refusals.push({ betId: bet.id, reason });
    else if (placedMs <= now) pending.push({ bet, placedMs });
  }
  if (!pending.length) return result;
  const rowsByBet = fillfair.sameLineRows(pending.map((item) => item.bet), state, boardLines);
  for (const { bet, placedMs } of pending) {
    const rows = rowsByBet.get(bet.id);
    if (!rows || !rows.length) continue;
    const decision = decide(bet, placedMs, rows, state, history || {}, leagueObservingSince || {});
    if (decision.row) result.saves.push(decision.row);
    else result.refusals.push({ betId: bet.id, reason: decision.reason });
  }
  return result;
}

module.exports = { REFUSED, captureFillBooks, refusalOf, entryAtFill, isOwnBook };
