// Unabated Ticket server — the board a live Place quote reads (Cal 2026-11-06): Unabated's live board
// for one game, with the venue ids its pregame board carried laid over it.
// Pure functions: no I/O, no clock.
//
// Why: Unabated's free feed sends its in-game lines (lg<id>:pt<n>:live) with
// no venue ids at all, no Kalshi sourceKey and no Novig sourceData, and no
// alternate lines (live recon of JAX @ BAL, 2026-11-06). The pregame board
// names each Kalshi contract and Novig outcome on the game's rungs, and those
// name the venue's game and market: Kalshi's tickers derive from any rung's
// series, event and team code (taker.derivedMainNumberContract,
// spreadTeamCodes), and Novig's event from any rung's outcome. So the runner
// remembers each game's pregame ids (rememberPregameVenueIds) and a live quote
// re-points them at the live line of the same market and side
// (liveQuoteBoard). The venue's own market read still checks every contract
// before an order (taker.marketMismatch / novigMarketMismatch).

// What a remembered rung target needs of its pregame main line to find the
// live line of the same market and side.
function mainShapeOf(line) {
  return { periodTypeId: line.periodTypeId, betTypeId: line.betTypeId, sideIndex: line.sideIndex, bookId: line.bookId };
}

// A game's pregame venue ids and the shape of every main line they hang off:
// {venueIds, mains: {lineKey: {periodTypeId, betTypeId, sideIndex, bookId}}},
// or null when the game carries none.
//   event  the pregame board's event (its venueIds); lines  the game's pregame lines
function rememberPregameVenueIds(event, lines) {
  const venueIds = event && event.venueIds;
  if (!venueIds) return null;
  const targets = [...Object.values(venueIds.kalshiContracts || {}), ...Object.values(venueIds.novigOutcomes || {})];
  if (targets.length === 0) return null;
  const byKey = new Map((lines || []).map((line) => [line.key, line]));
  const mains = {};
  for (const target of targets) {
    const main = byKey.get(target.mainKey);
    if (main) mains[target.mainKey] = mainShapeOf(main);
  }
  return { venueIds, mains };
}

// The one live line of each remembered main's market and side, by main key.
// Two live lines at one market and side (should not happen) leave it out.
function liveMainsByPregameKey(remembered, liveLines) {
  const byShape = new Map();
  for (const line of liveLines) {
    if (line.isAlt) continue;
    const shape = JSON.stringify(mainShapeOf(line));
    byShape.set(shape, byShape.has(shape) ? null : line);
  }
  const found = new Map();
  for (const [mainKey, shape] of Object.entries(remembered.mains)) {
    const line = byShape.get(JSON.stringify(shape));
    if (line) found.set(mainKey, line);
  }
  return found;
}

// One venue's remembered targets re-pointed at the live mains: a rung's
// target keeps its number and side, its mainKey becomes the live line's, and
// its lineKey the live line's only when the rung sits on the live number.
// A target whose market has no live line is left out.
function repointTargets(targets, liveMains) {
  const repointed = {};
  for (const [venueId, target] of Object.entries(targets || {})) {
    const live = liveMains.get(target.mainKey);
    if (!live) continue;
    const onLiveNumber = target.points != null && target.points === live.points;
    repointed[venueId] = { lineKey: onLiveNumber ? live.key : null, mainKey: live.key, points: target.points, sideIndex: target.sideIndex };
  }
  return repointed;
}

// The live board with one game's remembered pregame venue ids laid over it,
// re-pointed at its live lines. Without a memory of the game, the live board
// as it is (its event has no ids; contracts that need them refuse).
//   liveBoard  the scanner's live board of the game's league
//   remembered rememberPregameVenueIds' answer, or null
function liveQuoteBoard(liveBoard, eventId, remembered) {
  const event = liveBoard.events[eventId];
  if (!event || !remembered) return liveBoard;
  const liveLines = Object.values(liveBoard.lines).filter((line) => line.eventId === eventId);
  const liveMains = liveMainsByPregameKey(remembered, liveLines);
  // Ids the live lines carry themselves (none in Unabated's feed as of 2026-11-06) win over remembered ones.
  const own = event.venueIds || {};
  const suffixes = [...(remembered.venueIds.kalshiEventSuffixes || []), ...(own.kalshiEventSuffixes || [])];
  const venueIds = {
    ...remembered.venueIds,
    kalshiEventSuffixes: Array.from(new Set(suffixes)),
    kalshiContracts: { ...repointTargets(remembered.venueIds.kalshiContracts, liveMains), ...(own.kalshiContracts || {}) },
    novigOutcomes: { ...repointTargets(remembered.venueIds.novigOutcomes, liveMains), ...(own.novigOutcomes || {}) },
  };
  return { ...liveBoard, events: { ...liveBoard.events, [eventId]: { ...event, venueIds } } };
}

module.exports = { rememberPregameVenueIds, liveQuoteBoard };
