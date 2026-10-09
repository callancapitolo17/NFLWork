// Unabated Ticket server — the body of GET /edges.json, from what the runner
// holds: the scanner's feed state and line history, the open bets from the
// bets service, and the settings row. Pure: no I/O, no clock (the runner
// passes `now`), so the node tests drive it on fixtures.
//
// Every number and word comes from the panel's own modules (edgerows.js,
// tailflex.js, betsview.js, kelly.js): the list, the stakes (sized against
// the open straights and open BFA teasers), the tail-flex rank that picks a
// card's best line, the rail, the badges, the related bets and the
// edge-move tag are what the panel's Edges tab shows for the same settings.
// Pregame only — the Live block needs the logged-in Unabated screen
// (page.js) and has no server equivalent.
//
// Shape (one object; ms epochs unless named ...At ISO):
//   {generatedAt, settings: {source, error, okAt, updatedAt, stake: {bankroll,
//    multiplier}, edges: {...edge settings, bookIds: "default" | null (every
//    live book) | [ids]}}, scanner: {phase, error, leagues, leaguesLoaded,
//    leagueErrors, lineCount, altLineCount, eventCount, snapshotBuiltAt,
//    staleLeagues, loading, lastSnapshotAt}, betsService: {okAt, error,
//    unreachableSince, generatedAt, openBets, sources, boardLineCount,
//    unmatched: [{betId, reason, attachable, needsGame, needsFix, dismissed}]
//    (the open bets no board game matches, bets.unmatchedReasons — the
//    panel's Bets tab lists; a dismissed one never needs a game)},
//    books: {mode, ids, names, liveCount, live: [{id, name}] (every live
//    book, for the phone's book picker)}, tailFlex ("tail flex: NFL spr 7.1% · …", the panel's
//    header line, "" when nothing lists), grouped, unit: "cards" | "lines",
//    total, maxShown, items: [card] when grouped else [row]}
//   card  {key, sideName, eventId, league, leagueLabel, eventName, eventStart,
//          eventStartMs, awayTeam, homeTeam, betType, period, bookCount,
//          lineCount, best: row, others: [row]}  (best is the highest rank
//          score; the panel's expander shows each of `others`, in rank
//          order, with its standalone `stake`, not its rail)
//   row   the selectEdges fields the panel renders, plus edgeTier, stake (the
//         standalone, liquidity-capped Kelly stake), rankScore (tail-flex EV
//         dollars, null when unranked), bookProb, advice {kind, bet, alone,
//         verb, held, against, teasers {held, against}, reason, cappedAt, uncapped},
//         rail {text, note, atSize}, badges [{kind, text, title?}], related
//         [{tier, inMath, tag, text, fairThen, title?}], move {kind, label,
//         detail, sinceFill} | null

"use strict";

const kelly = require("../extension/kelly.js");
const betsLib = require("../extension/bets.js");
const betsView = require("../extension/betsview.js");
const edgeRows = require("../extension/edgerows.js");
const pricecheck = require("../extension/pricecheck.js");

// The settings row as the panel's settings objects (edgerows.js, shared with
// the panel's sync); re-exported here for the runner and its tests.
const { settingsFromService } = edgeRows;

function bookProbOrNull(row) {
  try {
    return kelly.bookProbOf({ bookPrice: row.price, sourceFormat: row.sourceFormat, sourcePrice: row.sourcePrice });
  } catch (_error) {
    return null;
  }
}

// The edge-move tag as data, or null (edgeRows.moveTag: only on a line held in this direction).
function moveView(row, moveContext) {
  const tag = edgeRows.moveTag(row, row.bet, moveContext);
  if (!tag) return null;
  return { kind: tag.kind, label: tag.label, detail: tag.detail, sinceFill: Boolean(tag.reading.baseline) };
}

// One sized Edges row as JSON: no venueIds map (every row of an event shares
// it by reference and the phone never reads it), no bet records beyond their words.
function rowView(row, moveContext) {
  const advice = row.bet.advice;
  return {
    key: row.key, league: row.league, leagueLabel: row.leagueLabel, eventId: row.eventId, eventName: row.eventName,
    eventStart: row.eventStart, eventStartMs: row.eventStartMs, awayTeam: row.awayTeam, homeTeam: row.homeTeam,
    rotation: row.rotation ?? null, marketId: row.marketId,
    betTypeId: row.betTypeId, betType: row.betType, periodTypeId: row.periodTypeId, period: row.period,
    sideIndex: row.sideIndex, sideKey: row.sideKey, sideLabel: row.sideLabel, points: row.points,
    isAlt: row.isAlt, mainPoints: row.mainPoints,
    book: { id: row.book.id, name: row.book.name }, price: row.price, sourceFormat: row.sourceFormat, sourcePrice: row.sourcePrice,
    bookProb: bookProbOrNull(row), fair: row.fair, edgePct: row.edgePct, edgeTier: edgeRows.edgeTier(row.edgePct),
    liquidity: row.liquidity, modifiedMs: row.modifiedMs, isBlurred: row.isBlurred,
    openerPrice: row.openerPrice, openerPoints: row.openerPoints,
    stake: row.stake, rankScore: row.rankScore,
    advice: {
      kind: advice.kind, bet: advice.bet, alone: advice.alone, verb: advice.verb,
      held: advice.held, against: advice.against, teasers: advice.teasers, reason: advice.reason, cappedAt: advice.cappedAt,
      uncapped: advice.uncapped ?? null,
    },
    rail: edgeRows.stakeRail(row),
    badges: betsView.badges(row.bet),
    related: betsView.relatedLines(row.bet, moveContext.fillFairIndex),
    move: moveView(row, moveContext),
    // Other books at the same number (pricecheck.js): the comparison and its tag.
    priceCheck: row.priceCheck ?? null,
    priceCheckTag: pricecheck.priceCheckTag(row.priceCheck),
  };
}

function cardView(group, moveContext) {
  return {
    key: group.key, sideName: group.sideName, eventId: group.eventId, league: group.league, leagueLabel: group.leagueLabel,
    eventName: group.eventName, eventStart: group.eventStart, eventStartMs: group.eventStartMs,
    awayTeam: group.awayTeam, homeTeam: group.homeTeam, betType: group.betType, period: group.period,
    bookCount: group.bookCount, lineCount: group.rows.length,
    best: rowView(group.best, moveContext),
    others: group.rows.slice(1).map((row) => rowView(row, moveContext)),
  };
}

function scannerView(status) {
  if (!status) return { phase: "starting", error: null };
  return {
    phase: status.phase, error: status.error, leagues: status.leagues, leaguesLoaded: status.leaguesLoaded,
    leagueErrors: status.leagueErrors, lineCount: status.lineCount, altLineCount: status.altLineCount,
    eventCount: status.eventCount, snapshotBuiltAt: status.snapshotBuiltAt, staleLeagues: status.staleLeagues,
    loading: status.loading, lastSnapshotAt: status.lastSnapshotAt,
  };
}

// The open bets no board game matches, as the panel's Bets tab lists them
// (bets.unmatchedReasons), by bet id: the phone joins them to its own
// /bets.json. `knownStarts` is the runner's {betId: startMs} memory of each
// bet's matched game (bets.matchedStarts), as the panel keeps it;
// `dismissedIds` the bets Dismissed on either page (bets.duckdb::bet_dismissals),
// which stop flagging red and list folded with Restore.
function unmatchedView(betRecords, boardLines, knownStarts, now, dismissedIds) {
  return betsLib.unmatchedReasons(betRecords, boardLines, now, { dismissedIds: dismissedIds || [], knownStarts: knownStarts || {} })
    .map(({ bet, reason, attachable, needsGame, needsFix, dismissed }) => ({ betId: bet.id, reason, attachable, needsGame, needsFix, dismissed: Boolean(dismissed) }));
}

// The /edges.json body.
//   input  {feedState, scannerStatus, history, betRecords, fillFairIndex, knownStarts, dismissedIds,
//           stakeSettings, edgeSettings, settingsStatus {source, error, okAt, updatedAt},
//           betsStatus {okAt, error, unreachableSince, generatedAt, sources},
//           boardLines, ladderReaderOf, teasers (teaser.openTeasers),
//           measurement (tailflex.measureTailFlex over feedState; null with no feed), now}
function buildEdgesPayload(input) {
  const { feedState, edgeSettings, stakeSettings, now } = input;
  // No Unabated tab on a server, so "follow Unabated" is every live book.
  const effective = edgeRows.effectiveFilter(edgeSettings, feedState, null, null);
  const rows = edgeRows.listedEdgeRows(feedState, {
    edgeSettings, stakeSettings, effective, records: input.betRecords,
    boardLines: input.boardLines, ladderReaderOf: input.ladderReaderOf, teasers: input.teasers,
    measurement: input.measurement, now,
  });
  const grouped = edgeSettings.groupByMarket;
  const items = grouped ? edgeRows.groupEdgeRows(rows, edgeSettings.sortBy) : rows;
  const moveContext = { history: input.history || {}, fillFairIndex: input.fillFairIndex, now };
  const bookName = (id) => (feedState && feedState.books[id] ? feedState.books[id].name : `book ${id}`);
  const ids = effective.bookIds ? Array.from(effective.bookIds).sort((a, b) => a - b) : null;
  const liveBooks = edgeRows.liveBooks(feedState);
  return {
    generatedAt: new Date(now).toISOString(),
    settings: {
      ...input.settingsStatus,
      stake: stakeSettings,
      edges: { ...edgeSettings, bookIds: edgeSettings.bookIds === undefined ? "default" : edgeSettings.bookIds },
    },
    scanner: scannerView(input.scannerStatus),
    betsService: {
      ...input.betsStatus, openBets: input.betRecords.filter((record) => record.status === "open").length,
      boardLineCount: input.boardLines.length, unmatched: unmatchedView(input.betRecords, input.boardLines, input.knownStarts, now, input.dismissedIds),
    },
    books: {
      mode: effective.mode, ids, names: ids ? ids.map(bookName) : null, liveCount: liveBooks.length,
      live: liveBooks.map((book) => ({ id: book.id, name: book.name })),
    },
    tailFlex: feedState ? edgeRows.describeTailFlex(rows, input.measurement) : "",
    priceCheck: pricecheck.describePriceChecks(rows),
    grouped,
    unit: grouped ? "cards" : "lines",
    total: items.length,
    maxShown: edgeRows.MAX_EDGE_ROWS,
    items: items.slice(0, edgeRows.MAX_EDGE_ROWS).map((item) => (grouped ? cardView(item, moveContext) : rowView(item, moveContext))),
  };
}

module.exports = { EDGE_SETTING_KEYS: edgeRows.EDGE_SETTING_KEYS, settingsFromService, rowView, cardView, unmatchedView, buildEdgesPayload };
