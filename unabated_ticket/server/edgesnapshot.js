// Why each edge was right or wrong (Cal, 2026-10-22): per straight bet, what
// the market looked like at the fill and at the close, so a leak can later be
// told apart — a market-wide mispricing, a rich alt rung, or an early-week
// error. Stored only; nothing reads it back but /bets.json?snapshots=1.
//
// Fill snapshot (one per bet, insert-only, bets.duckdb::bet_fill_snapshots):
//   books        every book's price on BOTH sides at the bet's number (the
//                bet's side at its number, the other side at the opposite
//                number), read off the line history at the fill; raw prices,
//                devigged in analysis, never here.
//   main         Unabated's own main number for the bet's side and its fair
//                there, and the bet's points minus that main number.
//   edge move    the edgemove.js tag the Edges tab showed at the fill, on the
//                bet's own venue's line when the board carries it.
// Close snapshot (newest pre-start reading wins, bet_close_snapshots): the
// same books record off the board's current state, sent with the closing
// fairs (server/closefair.js) while the game is still to start.
//
// Same rules as fillbooks.js: only while the scanner watched through the fill,
// never a backfill. A book counts while its line is on the board, unblurred,
// at a live book; Unabated's own line is not a book (its fair is `main`).
//
// Pure: no fetch, no clock, no timers. Required by server/runner.js and
// tests/edgesnapshot.test.js.

"use strict";

const feed = require("../extension/feed.js");
const fillfair = require("../extension/fillfair.js");
const edgemove = require("../extension/edgemove.js");
const fillbooks = require("./fillbooks.js");

const STATUS_ON_BOARD = 1;
const BET_TYPE_MONEYLINE = 1;
const BET_TYPE_SPREAD = 2;

// eventId -> that game's lines, one walk of each league the events are in
// (the scanner's per-league index), never the whole board.
//   eventIds       the games wanted
//   state          the scanner's feed state (events carry their leagueId)
//   leagueLinesOf  leagueId -> that league's lines (scanner.getLeagueLines)
function eventLinesOf(eventIds, state, leagueLinesOf) {
  const wanted = new Set(eventIds);
  const byEvent = new Map();
  const leagueIds = new Set();
  for (const eventId of wanted) {
    const event = state.events[eventId];
    if (event) leagueIds.add(event.leagueId);
  }
  for (const leagueId of leagueIds) {
    for (const line of leagueLinesOf(leagueId)) {
      if (!wanted.has(line.eventId)) continue;
      if (!byEvent.has(line.eventId)) byEvent.set(line.eventId, []);
      byEvent.get(line.eventId).push(line);
    }
  }
  return byEvent;
}

// The other side's number: a spread's sign flips, a total's stays, a moneyline has none.
function oppositePoints(line) {
  if (line.betTypeId === BET_TYPE_MONEYLINE || line.points == null) return null;
  return line.betTypeId === BET_TYPE_SPREAD ? -line.points : line.points;
}

// The same book's line on the other side of the same market at the opposite number, or null.
function otherSideOf(ownLine, eventLines) {
  const points = oppositePoints(ownLine);
  for (const line of eventLines || []) {
    if (line.bookId !== ownLine.bookId || line.periodTypeId !== ownLine.periodTypeId || line.betTypeId !== ownLine.betTypeId) continue;
    if (line.sideIndex == null || line.sideIndex === ownLine.sideIndex) continue;
    if (ownLine.betTypeId === BET_TYPE_MONEYLINE || line.points === points) return line;
  }
  return null;
}

function lineOnBoard(line, state) {
  if (!line || line.bookId === feed.UNABATED_LINE_BOOK_ID) return false;
  if (line.statusId !== STATUS_ON_BOARD || line.isBlurred) return false;
  const book = state.books[line.bookId];
  return Boolean(book && book.isLive);
}

function priceFields(observed, liquidity) {
  return {
    price: observed.price, sourcePrice: observed.sourcePrice ?? null, sourceFormat: observed.sourceFormat ?? 1,
    liquidity: liquidity ?? null,
  };
}

// Every book's two sides at the bet's number, one entry per book, sorted by
// book id: {bookId, book, exchange, venue (the bet's own book), side, other}.
// `priceAt(line)` is the price to record for a board line — {price,
// sourcePrice, sourceFormat, liquidity} — or null when there is none (at the
// fill: the line was first seen after it). A book with neither side is left out.
function booksRecord(ownRows, bet, state, eventLines, priceAt) {
  const books = [];
  for (const row of ownRows) {
    const ownLine = state.lines[row.key];
    if (!lineOnBoard(ownLine, state)) continue;
    const book = state.books[ownLine.bookId];
    const otherLine = otherSideOf(ownLine, eventLines.get(ownLine.eventId));
    const side = priceAt(ownLine);
    const other = lineOnBoard(otherLine, state) ? priceAt(otherLine) : null;
    if (!side && !other) continue;
    books.push({
      bookId: ownLine.bookId, book: book.name, exchange: book.hasLiquidity === true,
      venue: fillbooks.isOwnBook(ownLine, book, bet.venue), side, other,
    });
  }
  return books.sort((a, b) => a.bookId - b.bookId);
}

// Unabated's own main line on the bet's side of the market, or null.
function unabatedMainLine(ownLine, eventLines) {
  for (const line of eventLines || []) {
    if (line.bookId !== feed.UNABATED_LINE_BOOK_ID || line.isAlt) continue;
    if (line.periodTypeId === ownLine.periodTypeId && line.betTypeId === ownLine.betTypeId && line.sideIndex === ownLine.sideIndex) return line;
  }
  return null;
}

// {mainPoints, mainFairAmerican, pointsFromMain} at the fill, nulls when
// Unabated's main line was not seen by then or its number moved since (the
// history restarts on a new number). A moneyline has no number.
function mainAtFill(bet, ownLine, eventLines, history, placedMs) {
  const none = { mainPoints: null, mainFairAmerican: null, pointsFromMain: null };
  if (ownLine.betTypeId === BET_TYPE_MONEYLINE) return none;
  const mainLine = unabatedMainLine(ownLine, eventLines);
  const entry = mainLine ? fillbooks.entryAtFill(history[mainLine.key], placedMs) : null;
  if (!entry || entry.points == null) return none;
  const fair = fillfair.isWholeAmerican(entry.bacr) ? entry.bacr : null;
  const pointsFromMain = typeof bet.points === "number" ? Math.round((bet.points - entry.points) * 100) / 100 : null;
  return { mainPoints: entry.points, mainFairAmerican: fair, pointsFromMain };
}

// The edgemove.js tag at the fill: the history up to the placement, read as
// the Edges tab read it then. On the bet's own venue's line when the board
// carries it, else the preferred line at the number.
function edgeMoveAtFill(bet, ownRows, state, history, placedMs) {
  const ordered = fillfair.byPreference(ownRows, bet);
  const venueRow = ordered.find((row) => {
    const line = state.lines[row.key];
    return line && fillbooks.isOwnBook(line, state.books[line.bookId] || {}, bet.venue);
  });
  const row = venueRow || ordered[0];
  const entries = (history[row.key] || []).filter((entry) => entry.at <= placedMs);
  const move = edgemove.edgeMove(entries, placedMs);
  return { edgeMoveLineKey: row.key, edgeMove: move.kind, edgeMoveFairDelta: move.fairDelta, edgeMovePriceDelta: move.priceDelta };
}

// One capture pass: {saves, refusals}, fillbooks.captureFillBooks' contract.
// A save is the row POST /fill_snapshots.json takes — {betId, placedAt,
// lineKey, points, mainPoints, mainFairAmerican, pointsFromMain,
// edgeMoveLineKey, edgeMove, edgeMoveFairDelta, edgeMovePriceDelta, books}.
//   leagueLinesOf  leagueId -> that league's lines (scanner.getLeagueLines)
//   the rest as captureFillBooks
function captureFillSnapshots({ records, skipIds, state, boardLines, history, observingSince, leagueObservingSince, leagueLinesOf, now }) {
  const result = { saves: [], refusals: [] };
  if (!state || typeof observingSince !== "number") return result;
  const pending = [];
  for (const bet of records || []) {
    if (bet.status !== "open" || bet.unmatchable || skipIds.has(bet.id)) continue;
    const placedMs = Date.parse(bet.placedAt);
    const reason = fillbooks.refusalOf(bet, placedMs, observingSince, now);
    if (reason) result.refusals.push({ betId: bet.id, reason });
    else if (placedMs <= now) pending.push({ bet, placedMs });
  }
  if (!pending.length) return result;
  const rowsByBet = fillfair.sameLineRows(pending.map((item) => item.bet), state, boardLines);
  const eventIds = [...rowsByBet.values()].map((rows) => rows[0].eventId);
  const eventLines = eventLinesOf(eventIds, state, leagueLinesOf);
  for (const { bet, placedMs } of pending) {
    const rows = rowsByBet.get(bet.id);
    if (!rows || !rows.length) continue;
    const leagueSince = (leagueObservingSince || {})[rows[0].leagueId];
    if (typeof leagueSince !== "number" || placedMs < leagueSince) {
      result.refusals.push({ betId: bet.id, reason: fillbooks.REFUSED.leagueGap });
      continue;
    }
    const lineHistory = history || {};
    const priceAt = (line) => {
      const entry = fillbooks.entryAtFill(lineHistory[line.key], placedMs);
      return entry ? priceFields(entry, null) : null;
    };
    const preferred = fillfair.byPreference(rows, bet)[0];
    result.saves.push({
      betId: bet.id, placedAt: bet.placedAt, lineKey: preferred.key, points: preferred.points ?? null,
      ...mainAtFill(bet, state.lines[preferred.key], eventLines.get(preferred.eventId), lineHistory, placedMs),
      ...edgeMoveAtFill(bet, rows, state, lineHistory, placedMs),
      books: booksRecord(rows, bet, state, eventLines, priceAt),
    });
  }
  return result;
}

// The close snapshot of each closing-fair row (closefair.closingFairRows):
// {betId, observedAt, eventStart, books} off the board as it is now, with
// each line's resting liquidity. Same observedAt and eventStart as the
// closing fair, so the service keeps them by the same newest-pre-start rule.
//   fairRows   the closing-fair rows of this pass
//   rowsByBet  fillfair.sameLineRows of the same bets
function closeSnapshotRows({ fairRows, records, rowsByBet, state, leagueLinesOf }) {
  if (!fairRows.length) return [];
  const betsById = new Map((records || []).map((bet) => [bet.id, bet]));
  const eventIds = fairRows.map((row) => (rowsByBet.get(row.betId) || [])[0]).filter(Boolean).map((row) => row.eventId);
  const eventLines = eventLinesOf(eventIds, state, leagueLinesOf);
  const priceAt = (line) => priceFields(line, line.liquidity);
  const out = [];
  for (const fair of fairRows) {
    const bet = betsById.get(fair.betId);
    const rows = rowsByBet.get(fair.betId);
    if (!bet || !rows || !rows.length) continue;
    out.push({
      betId: fair.betId, observedAt: fair.fairObservedAt, eventStart: fair.eventStart,
      books: booksRecord(rows, bet, state, eventLines, priceAt),
    });
  }
  return out;
}

// A close snapshot is resent only when its reading time or its books changed.
function closeSignatureOf(row) {
  return `${row.observedAt}|${JSON.stringify(row.books)}`;
}

function unsentCloseRows(rows, sent) {
  return rows.filter((row) => sent.get(row.betId) !== closeSignatureOf(row));
}

function markCloseSent(rows, sent) {
  for (const row of rows) sent.set(row.betId, closeSignatureOf(row));
}

module.exports = {
  captureFillSnapshots, closeSnapshotRows, unsentCloseRows, markCloseSent,
  eventLinesOf, otherSideOf, oppositePoints,
};
