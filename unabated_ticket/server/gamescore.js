// The score line on a Bet Tracker Live card: a game's latest score, clock and
// period off Unabated's free league snapshot (feed.gameStateOf), as words.
// Pure: no DOM, no fetch, no disk.
//
// Input  gameState {statusId, awayScore, homeScore, gameClock, periodId,
//        overtimeNumber, snapshotAtMs, awayAbbr, homeAbbr} or null, the
//        league path ("nfl", ...), now ms
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
// Overtime the same way, off overtimeNumber (NAMED_OVERTIMES).

"use strict";

const live = require("../extension/live.js");

// Read live 2026-11-06 (NFL, CFB, NCAAB, NHL in progress, NCAAB final); the
// extension's own between-quarters reader (page.js) uses 3 as final too.
const STATUS_IN_PROGRESS = 2;
const STATUS_FINAL = 3;
// A score older than this reads as missing, not as the current score.
const STALE_AFTER_MS = 5 * 60 * 1000;

// League path -> the eventPeriodTypeIds named there, with live.js's names
// (Q1-Q4 = 4-7, halves 2 / 3). Football: live.js reads these at every NFL
// break; on 2026-11-06 an NFL game went 5 (10-3, 6:43) to 6 (17-3, 13:59)
// and CFB games sat at 2 with 00:00 (halftime). NCAAB: 3 with minutes left
// and at the final the same night. NBA, read 2026-11-06 every 30 s: one game
// went 5 (36-33, 8:53) to 2 with 00:00 at halftime (53-53) to 6 at 12:00,
// another 6 (65-91, 1:54) to 7 (75-96, 11:45). NBA overtime is unread, so it
// shows no name. Other sports are added once their own live games are read.
const NAMED_PERIOD_IDS = {
  nfl: [2, 4, 5, 6, 7],
  cfb: [2, 4, 5, 6, 7],
  cbb: [2, 3],
  nba: [2, 4, 5, 6, 7],
};

// League path -> overtimeNumber -> its name. Overtime is not a period id:
// while overtimeNumber is set, eventPeriodTypeId stays at the last regulation
// period (NCAAB read 2026-11-07: 3 into OT, overtimeNumber null -> 1 -> 2), so
// an unchecked league or count shows no name, never the regulation period's
// (football and NBA overtime unread).
const NAMED_OVERTIMES = {
  cbb: { 1: "OT", 2: "2OT" },
};

function periodNameOf(league, periodId, overtimeNumber) {
  if (overtimeNumber > 0) return NAMED_OVERTIMES[league]?.[overtimeNumber] ?? null;
  const ids = NAMED_PERIOD_IDS[league];
  if (!ids || !ids.includes(periodId)) return null;
  return live.PERIOD_ENDED[periodId] ?? null;
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
    period: status === "live" ? periodNameOf(league, gameState.periodId, gameState.overtimeNumber) : null,
    readAt: new Date(gameState.snapshotAtMs).toISOString(),
    staleAt: new Date(gameState.snapshotAtMs + STALE_AFTER_MS).toISOString(),
  };
}

module.exports = { STALE_AFTER_MS, scoreLineOf };
