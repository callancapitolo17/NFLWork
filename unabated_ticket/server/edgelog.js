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
// same number, price (Unabated's American price and the book's own, e.g. an
// exchange's cents) and fair, for the same scheduled start: a change to any of
// them starts a new segment, and a line that stops being an edge (off the
// board, or its edge under the floor) ends one. The runner sends a segment
// when it starts, again once LAST_SEEN_CONFIRM_MS has passed while it lasts,
// and once more when it ends, so the stored last_seen trails the truth by that
// long plus one league refresh. A restart forgets the open segments and starts
// new ones at the next snapshot: analysis joins a line's back-to-back segments
// that agree on every one of those fields.
//
// The lines are the Edges tab's own (feed.selectEdges at MIN_LOGGED_EDGE, alts
// included, lines changed within the taker's MAX_LINE_AGE_MS so a dead feed's
// stale price is not logged as an edge that lasted), less what this log cannot
// store: blurred lines (no real fair) and a fair or price that is not a whole
// American price. Unabated's price as shown (no exchange live reads), so every
// snapshot reads the same way. No liquidity floor: the resting liquidity at
// the segment's start is stored instead. Stored only: nothing touches a stake,
// an order or the Edges list.
//
// Pure: no fetch, no clock, no timers. Required by server/runner.js and
// tests/edgelog.test.js.

"use strict";

const feed = require("../extension/feed.js");
const fillfair = require("../extension/fillfair.js");
const laddersnapshot = require("./laddersnapshot.js");
const taker = require("./taker.js");

const EDGE_LOG_LEAGUE_IDS = laddersnapshot.LADDER_LEAGUE_IDS;
// Unabated's `ge` is a fraction: 0.01 = a 1% edge, the Edges tab's default floor.
const MIN_LOGGED_EDGE = 0.01;
const FULL_GAME_PERIOD_ID = 1;
const BET_TYPE_MONEYLINE = 1;
const LOGGED_BET_TYPE_IDS = [1, 2, 3];
// An open segment's last_seen is re-sent once this has passed while it lasts.
const LAST_SEEN_CONFIRM_MS = 5 * 60 * 1000;

// The edge lines of one league's snapshot: [{lineKey, eventId, leagueId,
// periodTypeId, bookId, betTypeId, sideIndex, points, isAlt, eventStartMs,
// price, sourcePrice, sourceFormat, fairAmerican, edge, liquidity,
// observedMs, signature}]. `signature` is what starts a new segment when it
// changes (a main line keeps its key when the book moves its number).
//   leagueId    the league to read
//   lines       that league's lines (scanner.getLeagueLines)
//   state       the scanner's feed state
//   observedMs  when the league's snapshot landed (scanner leagueLoadedAt)
function leagueEdgeLines({ leagueId, lines, state, observedMs }) {
  if (!state || !Number.isFinite(observedMs)) return [];
  // selectEdges walks every line it is given, and an alt's main is in its own
  // league, so it reads this league's lines alone, not the whole board.
  const leagueState = { ...state, lines: Object.fromEntries(lines.map((line) => [line.key, line])) };
  const rows = feed.selectEdges(leagueState, {
    leagueIds: new Set([leagueId]), periods: new Set([FULL_GAME_PERIOD_ID]), betTypes: new Set(LOGGED_BET_TYPE_IDS),
    minEdge: MIN_LOGGED_EDGE, includeAlts: true, now: observedMs, maxLineAgeMs: taker.MAX_LINE_AGE_MS, minLiquidityToWin: null,
  });
  const out = [];
  for (const row of rows) {
    const line = state.lines[row.key];
    if (line.isBlurred === true || (line.sideIndex !== 0 && line.sideIndex !== 1)) continue;
    if (!fillfair.isWholeAmerican(line.bacr) || !fillfair.isWholeAmerican(line.price)) continue;
    const points = line.betTypeId === BET_TYPE_MONEYLINE ? null : line.points;
    const eventStartMs = state.events[line.eventId].eventStart;
    out.push({
      lineKey: line.key, eventId: line.eventId, leagueId: line.leagueId, periodTypeId: line.periodTypeId, bookId: line.bookId,
      betTypeId: line.betTypeId, sideIndex: line.sideIndex, points, isAlt: line.isAlt === true, eventStartMs,
      price: line.price, sourcePrice: line.sourcePrice ?? null, sourceFormat: line.sourceFormat ?? null,
      fairAmerican: line.bacr, edge: line.ge, liquidity: line.liquidity ?? null, observedMs,
      signature: JSON.stringify([points, line.price, line.sourcePrice ?? null, line.bacr, eventStartMs]),
    });
  }
  return out;
}

function wireRow(segment) {
  const { edgeLine } = segment;
  return {
    lineKey: edgeLine.lineKey, segmentStart: new Date(segment.startMs).toISOString(), lastSeen: new Date(segment.lastSeenMs).toISOString(),
    eventId: edgeLine.eventId, leagueId: edgeLine.leagueId, periodTypeId: edgeLine.periodTypeId, bookId: edgeLine.bookId, betTypeId: edgeLine.betTypeId,
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
