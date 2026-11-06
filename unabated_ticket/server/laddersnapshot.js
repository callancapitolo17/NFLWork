// Unabated's full-game fair ladders, stored before kickoff (issue #264, part
// of the alt-spread model epic #263). The runner keeps only 10 minutes of
// line history in memory, wiped by every deploy; this keeps, per NFL / CFB /
// CBB game and market (spread, total, moneyline), Unabated's fair on both
// sides of every number any book posts, raw (never devigged here), in
// ladders.duckdb via the bets service's POST /ladders.json:
//   opener  the first reading the runner sends for that game and market,
//           never overwritten (the service keeps the first row ever stored,
//           so a restart's re-send is a no-op; a game already on the board
//           when storing began has its first sighting, not Unabated's opener)
//   close   the newest pregame reading: sent when the ladder changes, and
//           re-sent unchanged every CLOSE_CONFIRM_MS so the stored
//           observed_at stays near kickoff; the service keeps the newest
//           observation before the start. Whether it was near enough to the
//           start to count as the close is decided in analysis.
// Read-only: nothing here or downstream touches a stake or an order.
//
// A rung's fair: Unabated's `bacr` is its fair for THAT side at THAT number,
// and every book's line at the number carries it. The stored fair is the
// most common whole American value across those lines (ties: Unabated's own
// line's, else the lowest book id's), with how many lines carried a fair
// and how many agreed, so a disagreement stays visible. Blurred lines and
// Unabated's ±999900 placeholder are not fairs.
//
// Pure: no fetch, no clock, no timers. Required by server/runner.js and
// tests/laddersnapshot.test.js.

"use strict";

const feed = require("../extension/feed.js");
const fillfair = require("../extension/fillfair.js");

// NFL, CFB, CBB (Cal, 2026-11-06: CBB too, its season just started).
const LADDER_LEAGUE_IDS = [1, 2, 4];
const FULL_GAME_PERIOD_ID = 1;
const BET_TYPE_MONEYLINE = 1;
const LADDER_BET_TYPE_IDS = [1, 2, 3];
const READING_OPENER = "opener";
const READING_CLOSE = "close";
// An unchanged close is re-sent this often, so its observed_at trails the
// last pregame snapshot by at most this much.
const CLOSE_CONFIRM_MS = 10 * 60 * 1000;

function isUsableFair(line) {
  return line.isBlurred !== true && fillfair.isWholeAmerican(line.bacr) && !feed.isClampedFair(line.bacr);
}

// {fair, books, agreeing} of one (side, number)'s lines: the most common
// fair, how many lines carried a usable fair, how many of them had that one.
function rungFairOf(lines) {
  const countByFair = new Map();
  for (const line of lines) countByFair.set(line.bacr, (countByFair.get(line.bacr) || 0) + 1);
  const agreeing = Math.max(...countByFair.values());
  const tied = lines.filter((line) => countByFair.get(line.bacr) === agreeing);
  const unabatedLine = tied.find((line) => line.bookId === feed.UNABATED_LINE_BOOK_ID);
  const lowestBookLine = tied.reduce((best, line) => (line.bookId < best.bookId ? line : best));
  return { fair: (unabatedLine || lowestBookLine).bacr, books: lines.length, agreeing };
}

// Unabated's own main line on one side ({points, fair}), or nulls.
function mainOfSide(lines, sideIndex) {
  for (const line of lines) {
    if (line.bookId !== feed.UNABATED_LINE_BOOK_ID || line.isAlt || line.sideIndex !== sideIndex) continue;
    return { points: line.points ?? null, fair: isUsableFair(line) ? line.bacr : null };
  }
  return { points: null, fair: null };
}

// The ladder of one market of one game's full-game lines, or null when no
// line carries a usable fair. rungs: [[sideIndex, points, fair, books,
// agreeing], ...] sorted by side then number (a moneyline's points are null).
function marketLadderOf(eventLines, betTypeId) {
  const marketLines = eventLines.filter((line) => line.periodTypeId === FULL_GAME_PERIOD_ID && line.betTypeId === betTypeId
    && (line.sideIndex === 0 || line.sideIndex === 1));
  const linesByRung = new Map();
  for (const line of marketLines) {
    if (!isUsableFair(line)) continue;
    const points = betTypeId === BET_TYPE_MONEYLINE ? null : line.points;
    if (betTypeId !== BET_TYPE_MONEYLINE && (typeof points !== "number" || !Number.isFinite(points))) continue;
    const rungKey = `${line.sideIndex}|${points}`;
    if (!linesByRung.has(rungKey)) linesByRung.set(rungKey, { sideIndex: line.sideIndex, points, lines: [] });
    linesByRung.get(rungKey).lines.push(line);
  }
  if (!linesByRung.size) return null;
  const rungs = Array.from(linesByRung.values()).map(({ sideIndex, points, lines }) => {
    const { fair, books, agreeing } = rungFairOf(lines);
    return [sideIndex, points, fair, books, agreeing];
  }).sort((a, b) => a[0] - b[0] || (a[1] ?? 0) - (b[1] ?? 0));
  const main0 = mainOfSide(marketLines, 0);
  const main1 = mainOfSide(marketLines, 1);
  return {
    betTypeId, mainPoints0: main0.points, mainFair0: main0.fair, mainPoints1: main1.points, mainFair1: main1.fair, rungs,
  };
}

// Every ladder of one league's games still to start as of its snapshot:
// [{eventId, leagueId, periodTypeId, betTypeId, eventName, awayTeamId,
// homeTeamId, eventStartMs, observedMs, mainPoints0, mainFair0, mainPoints1,
// mainFair1, rungs, signature}].
//   linesByEvent  eventId -> that game's lines (the runner's per-league index)
//   events        the scanner state's events
//   observedMs    when the league's snapshot landed (scanner leagueLoadedAt)
function leagueLadders({ leagueId, linesByEvent, events, observedMs }) {
  const ladders = [];
  if (!Number.isFinite(observedMs)) return ladders;
  for (const [eventId, eventLines] of linesByEvent) {
    const event = events[eventId];
    if (!event || event.leagueId !== leagueId || !Number.isFinite(event.eventStart) || observedMs >= event.eventStart) continue;
    for (const betTypeId of LADDER_BET_TYPE_IDS) {
      const ladder = marketLadderOf(eventLines, betTypeId);
      if (!ladder) continue;
      ladders.push({
        eventId: event.eventId, leagueId, periodTypeId: FULL_GAME_PERIOD_ID, eventName: event.eventName ?? null,
        awayTeamId: event.awayTeamId ?? null, homeTeamId: event.homeTeamId ?? null,
        eventStartMs: event.eventStart, observedMs, ...ladder,
        signature: JSON.stringify([ladder.mainPoints0, ladder.mainFair0, ladder.mainPoints1, ladder.mainFair1, ladder.rungs]),
      });
    }
  }
  return ladders;
}

function ladderKeyOf(ladder) {
  return `${ladder.eventId}|${ladder.periodTypeId}|${ladder.betTypeId}`;
}

function rowOf(ladder, reading) {
  return {
    reading, eventId: ladder.eventId, leagueId: ladder.leagueId, periodTypeId: ladder.periodTypeId, betTypeId: ladder.betTypeId,
    eventName: ladder.eventName, awayTeamId: ladder.awayTeamId, homeTeamId: ladder.homeTeamId,
    eventStart: new Date(ladder.eventStartMs).toISOString(), observedAt: new Date(ladder.observedMs).toISOString(),
    mainPoints0: ladder.mainPoints0, mainFair0: ladder.mainFair0, mainPoints1: ladder.mainPoints1, mainFair1: ladder.mainFair1,
    rungs: ladder.rungs,
  };
}

// The rows to POST now: each ladder's opener until one was taken, and its
// close when the ladder changed since the last close taken or that one is
// CLOSE_CONFIRM_MS old. A game that started since its league's snapshot is
// left out (the service refuses a reading at or after the start anyway).
//   ladders  leagueLadders' output for the leagues read
//   sent     ladderKey -> {openerSent, closeSignature, closeObservedMs} (markSent)
function unsentRows(ladders, sent, now) {
  const rows = [];
  for (const ladder of ladders) {
    if (now >= ladder.eventStartMs) continue;
    const known = sent.get(ladderKeyOf(ladder)) || {};
    if (!known.openerSent) rows.push(rowOf(ladder, READING_OPENER));
    const changed = known.closeSignature !== ladder.signature;
    const confirmDue = !Number.isFinite(known.closeObservedMs) || ladder.observedMs - known.closeObservedMs >= CLOSE_CONFIRM_MS;
    if (changed || confirmDue) rows.push({ ...rowOf(ladder, READING_CLOSE), signature: ladder.signature });
  }
  return rows;
}

// Record the rows the service took, so they are not re-sent.
function markSent(rows, sent) {
  for (const row of rows) {
    const key = ladderKeyOf(row);
    const known = sent.get(key) || {};
    if (row.reading === READING_OPENER) known.openerSent = true;
    else Object.assign(known, { closeSignature: row.signature, closeObservedMs: Date.parse(row.observedAt) });
    known.eventStartMs = Date.parse(row.eventStart);
    sent.set(key, known);
  }
}

// Forget the ladders of games that started: nothing more is sent for them.
function forgetStarted(sent, now) {
  for (const [key, known] of sent) if (Number.isFinite(known.eventStartMs) && now >= known.eventStartMs) sent.delete(key);
}

// The row as POST /ladders.json takes it (the runner-only signature dropped).
function wireRow(row) {
  const { signature, ...wire } = row;
  return wire;
}

module.exports = {
  LADDER_LEAGUE_IDS, CLOSE_CONFIRM_MS, READING_OPENER, READING_CLOSE,
  marketLadderOf, leagueLadders, unsentRows, markSent, forgetStarted, wireRow,
};
