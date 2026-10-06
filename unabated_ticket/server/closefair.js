// Closing fairs for CLV (2026-10-05, Cal: "a flaw we're not using CLV").
//
// A bet's closing fair is Unabated's fair (`bacr`) for the bet's own side and
// number on the last snapshot before its game started. The runner keeps
// reading it while the game is still to start and POSTs each change to the
// bets service, which keeps the newest pre-start observation per bet
// (bets.duckdb::bet_closing_fairs). Once the game starts nothing more is
// sent, so the stored row is the close. The tracker decides whether the last
// observation was near enough to the start to count (trackerstats.js).
//
// The line is found the way fillfair.js finds a fill's: the bet's same_line
// rows on the board (any book at its number, its own venue first), so a bet
// on a number the market moved off is priced at that alt rung when the board
// carries it, and has no close when it does not.
//
// Pure: no fetch, no clock, no timers. Required by server/runner.js and
// tests/closefair.test.js.

"use strict";

const feed = require("../extension/feed.js");
const fillfair = require("../extension/fillfair.js");

// A row is resent only when what it says changed: the fair, or the time it
// was last confirmed (a new snapshot of its league).
function signatureOf(row) {
  return `${row.lineKey}|${row.fairAmerican}|${row.fairObservedAt}`;
}

// The open straight bets the board can price: not a parlay or teaser leg
// (a leg has no stake of its own to measure), not unmatchable.
function closingCandidates(records) {
  return (records || []).filter((bet) => bet.status === "open" && !bet.unmatchable && bet.isParlayLeg !== true
    && Number.isFinite(Date.parse(bet.placedAt)) && fillfair.onBoardMarket(bet));
}

// The first row, in preference order, with a whole American fair on a game
// still to start that the bet was placed before, and its league's snapshot
// time — {row, observedMs} or null. A league whose snapshot was built long
// ago (scanner staleLeagues: a stale CDN copy) is skipped: its load time would
// pass an old fair off as a fresh close.
function pregameRow(rows, bet, leagueLoadedAt, staleLeagues, now) {
  const placedMs = Date.parse(bet.placedAt);
  for (const row of fillfair.byPreference(rows, bet)) {
    const startMs = row.eventStartMs;
    if (!Number.isFinite(startMs) || now >= startMs || placedMs >= startMs) continue;
    if (!fillfair.isWholeAmerican(row.fair) || staleLeagues.has(row.leagueId)) continue;
    const observedMs = leagueLoadedAt[row.leagueId];
    if (!Number.isFinite(observedMs) || observedMs >= startMs) continue;
    return { row, observedMs };
  }
  return null;
}

// The rows POST /closing_fairs.json takes, one per open bet whose line the
// board prices right now before its game starts:
//   {betId, lineKey, points, fairAmerican, fairObservedAt, eventStart}, times ISO.
//   records        bet records (bets.js contract)
//   state          the scanner's feed state
//   boardLines     one describeLine row per event (runner boardLines())
//   leagueLoadedAt leagueId -> ms the league's last snapshot landed (scanner status)
//   staleLeagues   league ids whose snapshot build is stale (scanner status)
//   now            ms
function closingFairRows({ records, state, boardLines, leagueLoadedAt, staleLeagues, now }) {
  if (!state) return [];
  const candidates = closingCandidates(records);
  if (!candidates.length) return [];
  const rowsByBet = fillfair.sameLineRows(candidates, state, boardLines);
  const out = [];
  for (const bet of candidates) {
    const rows = rowsByBet.get(bet.id);
    if (!rows) continue;
    const found = pregameRow(rows, bet, leagueLoadedAt || {}, new Set(staleLeagues || []), now);
    if (!found) continue;
    out.push({
      betId: bet.id, lineKey: found.row.key, points: found.row.points ?? null, fairAmerican: found.row.fair,
      fairObservedAt: new Date(found.observedMs).toISOString(), eventStart: new Date(found.row.eventStartMs).toISOString(),
    });
  }
  return out;
}

// The rows that differ from what was last sent for their bet.
//   sent  Map betId -> signature, updated by markSent once a POST succeeds
function unsentRows(rows, sent) {
  return rows.filter((row) => sent.get(row.betId) !== signatureOf(row));
}

function markSent(rows, sent) {
  for (const row of rows) sent.set(row.betId, signatureOf(row));
}

// Drop what was sent for bets no longer open, so the map stays the size of the open book.
function forgetClosed(sent, records) {
  const openIds = new Set((records || []).filter((bet) => bet.status === "open").map((bet) => bet.id));
  for (const betId of sent.keys()) if (!openIds.has(betId)) sent.delete(betId);
}

// Unabated league ids of the open bets' leagues (bet.league is feed.LEAGUES'
// `path`; one path can name several ids, e.g. every soccer league), so the
// runner loads them even when the Edges settings leave them out.
function leagueIdsOfOpenBets(records) {
  const paths = new Set(closingCandidates(records).map((bet) => bet.league).filter(Boolean));
  return Object.entries(feed.LEAGUES).filter(([, league]) => paths.has(league.path)).map(([id]) => Number(id));
}

module.exports = { closingFairRows, unsentRows, markSent, forgetClosed, leagueIdsOfOpenBets };
