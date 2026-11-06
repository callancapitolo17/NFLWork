// The score line on a Bet Tracker Live card: a game's latest score, clock and
// period off Unabated's free league snapshot (feed.gameStateOf), as words.
// Pure: no DOM, no fetch, no disk.
//
// Input  gameState {statusId, awayScore, homeScore, gameClock, periodId,
//        snapshotAtMs, awayAbbr, homeAbbr} or null, the league path ("nfl",
//        ...), now ms
// Output {status, awayAbbr, homeAbbr, awayScore, homeScore, clock, period,
//        readAt, staleAt}: status "live" / "final" / "other" (upcoming or a
//        status not seen yet: no score is shown for it); staleAt is when the
//        reading stops counting, so the tracker can age a copy it kept.
//        {status: "stale"} once no snapshot in STALE_AFTER_MS carried the
//        game; null when the runner never read one.
// Side effects: none.
//
// Period names: eventPeriodTypeId is a per-sport number, named only for the
// leagues whose numbers were checked against live games; any other league or
// number shows no period name (the clock and score still show), never a guess.

"use strict";

const live = require("../extension/live.js");

// Read live 2026-11-06 (NFL, CFB, NCAAB, NHL in progress, NCAAB and NHL final); the
// extension's own between-quarters reader (page.js) uses 3 as final too.
const STATUS_IN_PROGRESS = 2;
const STATUS_FINAL = 3;
// A score older than this reads as missing, not as the current score.
const STALE_AFTER_MS = 5 * 60 * 1000;

// League path -> {eventPeriodTypeId: name} for the numbers checked against
// live games; other sports are added once their own live games are read.
//
// Football and NCAAB use live.js's names (Q1-Q4 = 4-7, halves 2 / 3): live.js
// reads these at every NFL break; on 2026-11-06 an NFL game went 5 (10-3,
// 6:43) to 6 (17-3, 13:59), CFB games sat at 2 with 00:00 (halftime), and
// NCAAB was at 3 with minutes left and at the final the same night.
//
// Hockey (league 6) read live 2026-11-06: games started 00:00Z were at 10
// with the clock counting down in the 3rd at 02:08Z and final at 10 by
// 03:46Z; a 01:00Z game sat at 8 with 00:00 at 02:08Z (end of the 1st); a
// 02:30Z game was at 9 with 14:10 left at 03:46Z. The league's own market
// periods are pt8 / pt9 / pt10. Overtime and the shootout were not seen, so
// they stay unnamed. The "nhl" path also carries Olympic hockey (league 11),
// not read live; no Olympic tournament before 2030.
function namesOf(periodIds) {
  return Object.fromEntries(periodIds.map((id) => [id, live.PERIOD_ENDED[id]]));
}
const PERIOD_NAMES = {
  nfl: namesOf([2, 4, 5, 6, 7]),
  cfb: namesOf([2, 4, 5, 6, 7]),
  cbb: namesOf([2, 3]),
  nhl: { 8: "P1", 9: "P2", 10: "P3" },
};

function periodNameOf(league, periodId) {
  return PERIOD_NAMES[league]?.[periodId] ?? null;
}

function statusOf(statusId) {
  if (statusId === STATUS_IN_PROGRESS) return "live";
  if (statusId === STATUS_FINAL) return "final";
  return "other";
}

function scoreLineOf(gameState, league, now) {
  if (!gameState) return null;
  if (gameState.snapshotAtMs === null || now - gameState.snapshotAtMs > STALE_AFTER_MS) return { status: "stale" };
  const status = statusOf(gameState.statusId);
  return {
    status,
    awayAbbr: gameState.awayAbbr ?? null,
    homeAbbr: gameState.homeAbbr ?? null,
    awayScore: gameState.awayScore,
    homeScore: gameState.homeScore,
    // A final game's clock is the feed's end marker ("E", "00:00"), not a time.
    clock: status === "live" ? gameState.gameClock : null,
    period: status === "live" ? periodNameOf(league, gameState.periodId) : null,
    readAt: new Date(gameState.snapshotAtMs).toISOString(),
    staleAt: new Date(gameState.snapshotAtMs + STALE_AFTER_MS).toISOString(),
  };
}

module.exports = { STALE_AFTER_MS, scoreLineOf };
