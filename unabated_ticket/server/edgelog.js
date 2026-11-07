// Every edge's life to kickoff (Cal, 2026-11-07, "Decide when to place
// bets"): should a bet wait for kickoff? Betting early gives up about a third
// of the edge (the fair moves toward the book), but waiting loses the bet
// whenever the book moves first. His own bets are too few to measure the
// second, so this logs every edge on the board, not only the ones he bet:
// for each NFL / CFB / CBB full-game line at a live book where Unabated shows
// at least MIN_LOGGED_EDGE, when the edge was there and at what price and
// fair, until the edge is gone or the game starts. Scored in analysis against
// the kickoff fair ladders.duckdb already keeps (laddersnapshot.js).
//
// A segment is a run of league snapshots in which one line was an edge at the
// same number, price and fair: a change to any of them starts a new segment,
// and a line that stops being an edge (off the board, or its edge under the
// floor) ends one. The runner sends a segment when it starts, again
// every LAST_SEEN_CONFIRM_MS while it lasts, and once more when it ends, so
// the stored last_seen trails the truth by at most that long. A restart forgets
// the open segments and starts new ones at the next snapshot: analysis joins
// a line's back-to-back segments at the same price and fair.
//
// The Edges tab's own line rules (feed.selectEdges): on the board, a live
// book, not Unabated's own line, not blurred, a fair that is not the ±999900 clamp, a
// valid American price, the game still to start, an alt only inside the
// 15-85% depth cap and off a flat-lined tail. Unabated's price as shown (no
// exchange live reads), so every snapshot reads the same way. No liquidity
// floor: the resting liquidity at the segment's start is stored instead.
// Stored only: nothing touches a stake, an order or the Edges list.
//
// Pure: no fetch, no clock, no timers. Required by server/runner.js and
// tests/edgelog.test.js.

"use strict";

const feed = require("../extension/feed.js");
const laddersnapshot = require("./laddersnapshot.js");

const EDGE_LOG_LEAGUE_IDS = laddersnapshot.LADDER_LEAGUE_IDS;
// Unabated's `ge` is a fraction: 0.01 = a 1% edge, the Edges tab's default floor.
const MIN_LOGGED_EDGE = 0.01;
const FULL_GAME_PERIOD_ID = 1;
const BET_TYPE_MONEYLINE = 1;
const LOGGED_BET_TYPE_IDS = new Set([1, 2, 3]);
const STATUS_ON_BOARD = 1;
const MIN_AMERICAN_MAGNITUDE = 100;
// An open segment's last_seen is re-sent this often while it lasts.
const LAST_SEEN_CONFIRM_MS = 5 * 60 * 1000;

function isLoggedAlt(line, state, flatRungKeys) {
  const main = state.lines[line.mainKey];
  if (main && main.points === line.points) return false;
  if (!feed.altWithinDepthCap(line.bacr, line.price)) return false;
  return !flatRungKeys.has(feed.rungKeyOf(line));
}

function isLoggedLine(line, state, observedMs, flatRungKeys) {
  if (line.bookId === feed.UNABATED_LINE_BOOK_ID || line.statusId !== STATUS_ON_BOARD || line.isBlurred === true) return false;
  if (line.periodTypeId !== FULL_GAME_PERIOD_ID || !LOGGED_BET_TYPE_IDS.has(line.betTypeId)) return false;
  if (line.sideIndex !== 0 && line.sideIndex !== 1) return false;
  if (typeof line.ge !== "number" || line.ge < MIN_LOGGED_EDGE) return false;
  if (!Number.isInteger(line.bacr) || feed.isClampedFair(line.bacr)) return false;
  if (!Number.isInteger(line.price) || Math.abs(line.price) < MIN_AMERICAN_MAGNITUDE) return false;
  if (line.betTypeId !== BET_TYPE_MONEYLINE && typeof line.points !== "number") return false;
  const book = state.books[line.bookId];
  if (!book || !book.isLive) return false;
  const event = state.events[line.eventId];
  if (!event || !Number.isFinite(event.eventStart) || observedMs >= event.eventStart) return false;
  return !line.isAlt || isLoggedAlt(line, state, flatRungKeys);
}

// The edge lines of one league's snapshot: [{lineKey, eventId, leagueId,
// bookId, betTypeId, sideIndex, points, isAlt, eventStartMs, price,
// sourcePrice, sourceFormat, fairAmerican, edge, liquidity, observedMs,
// signature}]. `signature` is what starts a new segment when it changes: the
// number (a main line keeps its key when the book moves it), price or fair.
//   lines       the league's lines (scanner.getLeagueLines)
//   state       the scanner's feed state (books, events, and lines for an alt's main)
//   observedMs  when the league's snapshot landed (scanner leagueLoadedAt)
function leagueEdgeLines({ lines, state, observedMs }) {
  if (!Number.isFinite(observedMs)) return [];
  const flatRungKeys = feed.flatFairRungKeys(lines);
  const out = [];
  for (const line of lines) {
    if (!isLoggedLine(line, state, observedMs, flatRungKeys)) continue;
    out.push({
      lineKey: line.key, eventId: line.eventId, leagueId: line.leagueId, bookId: line.bookId, betTypeId: line.betTypeId,
      sideIndex: line.sideIndex, points: line.betTypeId === BET_TYPE_MONEYLINE ? null : line.points, isAlt: line.isAlt === true,
      eventStartMs: state.events[line.eventId].eventStart, price: line.price, sourcePrice: line.sourcePrice ?? null,
      sourceFormat: line.sourceFormat ?? null, fairAmerican: line.bacr, edge: line.ge, liquidity: line.liquidity ?? null,
      observedMs, signature: `${line.points ?? ""}|${line.price}|${line.sourcePrice ?? ""}|${line.bacr}`,
    });
  }
  return out;
}

function wireRow(segment) {
  const { edgeLine } = segment;
  return {
    lineKey: edgeLine.lineKey, segmentStart: new Date(segment.startMs).toISOString(), lastSeen: new Date(segment.lastSeenMs).toISOString(),
    eventId: edgeLine.eventId, leagueId: edgeLine.leagueId, bookId: edgeLine.bookId, betTypeId: edgeLine.betTypeId,
    sideIndex: edgeLine.sideIndex, points: edgeLine.points, isAlt: edgeLine.isAlt, eventStart: new Date(edgeLine.eventStartMs).toISOString(),
    price: edgeLine.price, sourcePrice: edgeLine.sourcePrice, sourceFormat: edgeLine.sourceFormat,
    fairAmerican: edgeLine.fairAmerican, edge: edgeLine.edge, liquidity: edgeLine.liquidity,
  };
}

// Advance the open segments by one read of some leagues and return the rows
// to send now, each what POST /edge_segments.json takes. A segment counts as
// sent once its row is returned: the runner queues the rows and retries a
// failed POST itself (the service keeps the newest last_seen per segment).
//   edgeLines      leagueEdgeLines of every league read this pass
//   readLeagueIds  those leagues: an open segment of a league not read this
//                  pass (stale, or no longer loaded) is left as it is
//   segments       lineKey -> {edgeLine, startMs, lastSeenMs, sentLastSeenMs} (mutated)
//   now            ms; a segment whose game has started ends where it was last seen
function advanceSegments({ edgeLines, readLeagueIds, segments, now }) {
  const rows = [];
  const send = (segment) => {
    segment.sentLastSeenMs = segment.lastSeenMs;
    rows.push(wireRow(segment));
  };
  const endSegment = (lineKey, segment) => {
    if (segment.sentLastSeenMs !== segment.lastSeenMs) send(segment);
    segments.delete(lineKey);
  };
  const seen = new Set();
  for (const edgeLine of edgeLines) {
    // A league's lines are re-read only when its snapshot lands again, so a game can start in between.
    if (now >= edgeLine.eventStartMs) continue;
    seen.add(edgeLine.lineKey);
    const open = segments.get(edgeLine.lineKey);
    if (open && open.edgeLine.signature === edgeLine.signature) {
      open.lastSeenMs = Math.max(open.lastSeenMs, edgeLine.observedMs);
      if (open.lastSeenMs - open.sentLastSeenMs >= LAST_SEEN_CONFIRM_MS) send(open);
      continue;
    }
    if (open) endSegment(edgeLine.lineKey, open);
    const segment = { edgeLine, startMs: edgeLine.observedMs, lastSeenMs: edgeLine.observedMs, sentLastSeenMs: null };
    segments.set(edgeLine.lineKey, segment);
    send(segment);
  }
  const read = new Set(readLeagueIds);
  for (const [lineKey, segment] of Array.from(segments)) {
    const started = now >= segment.edgeLine.eventStartMs;
    if (started || (read.has(segment.edgeLine.leagueId) && !seen.has(lineKey))) endSegment(lineKey, segment);
  }
  return rows;
}

// The queue key of a row: a later row of the same segment replaces an unsent earlier one.
function segmentKeyOf(row) {
  return `${row.lineKey}|${row.segmentStart}`;
}

module.exports = {
  EDGE_LOG_LEAGUE_IDS, MIN_LOGGED_EDGE, LAST_SEEN_CONFIRM_MS, leagueEdgeLines, advanceSegments, segmentKeyOf,
};
