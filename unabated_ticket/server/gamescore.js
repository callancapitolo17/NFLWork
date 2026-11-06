// The score line on a Bet Tracker Live card: a game's latest score, clock and
// period off Unabated's free league snapshot (feed.gameStateOf), as words.
// Pure: no DOM, no fetch, no disk.
//
// Input  gameState {statusId, awayScore, homeScore, gameClock, periodId,
//        overtimeNumber, snapshotAtMs, awayAbbr, homeAbbr} or null, the league path ("nfl", ...),
//        now ms
// Output {status, awayAbbr, homeAbbr, awayScore, homeScore, clock, period, readAt} or
//        {status: "stale", readAt} when no snapshot in STALE_AFTER_MS carried
//        the game, or null when the runner has never read one.
// Side effects: none.
//
// Period names: eventPeriodTypeId is a bare per-sport number. Only numbers
// checked against a live game in that sport go in PERIOD_NAMES; a number not
// in it shows no period name (the clock and score still show), never a guess.

"use strict";

// Read live 2026-11-06 (NFL, CFB, NCAAB, NHL in progress, NCAAB final).
const STATUS_UPCOMING = 1;
const STATUS_IN_PROGRESS = 2;
const STATUS_FINAL = 3;
// A score older than this reads as missing, not as the current score.
const STALE_AFTER_MS = 5 * 60 * 1000;

// league path -> {eventPeriodTypeId: name}. Each sport's entries are added by
// its own check of a live game; none are verified yet.
const PERIOD_NAMES = {};

function periodNameOf(league, periodId) {
  const names = PERIOD_NAMES[league];
  if (!names || periodId === null) return null;
  return names[periodId] ?? null;
}

function statusOf(statusId) {
  if (statusId === STATUS_IN_PROGRESS) return "live";
  if (statusId === STATUS_FINAL) return "final";
  if (statusId === STATUS_UPCOMING) return "upcoming";
  return "unknown";
}

function scoreLineOf(gameState, league, now) {
  if (!gameState) return null;
  const readAt = gameState.snapshotAtMs === null ? null : new Date(gameState.snapshotAtMs).toISOString();
  if (gameState.snapshotAtMs === null || now - gameState.snapshotAtMs > STALE_AFTER_MS) return { status: "stale", readAt };
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
    readAt,
  };
}

module.exports = { STALE_AFTER_MS, PERIOD_NAMES, periodNameOf, scoreLineOf };
