"use strict";
// Novig Place (novig_place_plan.md): which Novig outcome a board line buys,
// and the outcome on the other side of the same market (the one whose
// resting orders set the price, and the one the bets service checks Cal does
// not already hold).
//
// Unabated carries Novig outcome ids only on spread and total rungs
// (feed.js noteRungVenueIds -> event.venueIds.novigOutcomes: id -> {lineKey,
// mainKey, points, sideIndex}); a moneyline has none, so it cannot be placed.
//
// Pure: no I/O, no clock.

const NOVIG_BOOK_ID = 89;
const BET_TYPE_SPREAD = 2;
const BET_TYPE_TOTAL = 3;

// The other side's number: a spread flips sign (-3.5 / +3.5), a total keeps
// it (over 47.5 / under 47.5).
function otherSidePoints(betTypeId, points) {
  return betTypeId === BET_TYPE_SPREAD ? -points : points;
}

// {outcomeId, otherOutcomeId} for a board line, or {reason, detail}.
//   line   a board line (feedState.lines value): {key, bookId, eventId,
//          periodTypeId, betTypeId, sideIndex, points}
//   board  the feed state: {lines, events}
function novigOutcomesOf(line, board) {
  if (line.bookId !== NOVIG_BOOK_ID) return { reason: "not_novig", detail: `book ${line.bookId} is not Novig` };
  if (line.betTypeId !== BET_TYPE_SPREAD && line.betTypeId !== BET_TYPE_TOTAL) {
    return { reason: "no_outcome", detail: "Unabated carries Novig outcome ids on spreads and totals only" };
  }
  const event = board.events[line.eventId];
  const outcomes = event && event.venueIds ? event.venueIds.novigOutcomes || {} : {};
  const sameMarket = (target) => {
    const main = board.lines[target.mainKey];
    return Boolean(main) && main.periodTypeId === line.periodTypeId && main.betTypeId === line.betTypeId;
  };
  const entries = Object.entries(outcomes).filter(([, target]) => sameMarket(target));
  const mine = entries.filter(([, target]) => target.lineKey === line.key
    && target.sideIndex === line.sideIndex && target.points === line.points);
  const otherPoints = otherSidePoints(line.betTypeId, line.points);
  const other = entries.filter(([, target]) => target.sideIndex === 1 - line.sideIndex && target.points === otherPoints);
  if (mine.length !== 1) {
    return { reason: "no_outcome", detail: `the board holds ${mine.length} Novig outcome ids for this line at ${line.points}` };
  }
  if (other.length !== 1) {
    return { reason: "no_outcome", detail: `the board holds ${other.length} Novig outcome ids for the other side at ${otherPoints}` };
  }
  if (other[0][0] === mine[0][0]) return { reason: "no_outcome", detail: `both sides read as outcome ${mine[0][0]}` };
  return { outcomeId: mine[0][0], otherOutcomeId: other[0][0] };
}

module.exports = { NOVIG_BOOK_ID, novigOutcomesOf, otherSidePoints };
