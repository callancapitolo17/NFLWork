// Unabated's fair ladder for one market of one game (#130): the chance the
// result lands above each half-point number, read off the scanner state's
// lines. Pure: no DOM, no fetch, no chrome.* — loaded as a plain <script> in
// panel.html after feed.js (exposes globalThis.UnabatedLadder) and via
// require() in tests/ladder.test.js.
//
// Inputs
//   lines  scannerState.lines values (feed.js): every book line and alt rung,
//          each {eventId, leagueId, periodTypeId, betTypeId, sideIndex,
//          points, bacr}. `bacr` is Unabated's fair for THAT side at THAT
//          number, as a whole American price; side 0 = away / Over, side 1 =
//          home / Under.
// Outputs a ladder {axis, periodTypeId, rungs: [[cut, probAbove], ...]}
//          sorted by cut, and probAbove(ladder, cut) -> {prob} or {reason}.
//          Nothing here writes anywhere.
//
// Two axes. "total": cut = the number, above = Over. "margin" = away score
// minus home score: an away spread at points a wins above -a, a home spread
// at points h wins below h, and the moneyline is the +0.5 (away) / -0.5
// (home) cut, so spreads and moneylines share one ladder.
//
// Only half-point rungs are kept: a whole-number rung's fair is conditional
// on no push, so a push is priced from the two half-point rungs around it.

(function (root) {
  "use strict";

  const kelly = typeof module !== "undefined" && module.exports ? require("./kelly.js") : root.UnabatedKelly;
  const feed = typeof module !== "undefined" && module.exports ? require("./feed.js") : root.UnabatedFeed;

  const AXIS_TOTAL = "total";
  const AXIS_MARGIN = "margin";
  const BET_TYPE_MONEYLINE = 1;
  const BET_TYPE_SPREAD = 2;
  const BET_TYPE_TOTAL = 3;
  const SIDE_AWAY_OR_OVER = 0;
  const MONEYLINE_CUT = 0.5;
  // Soccer's moneyline is three-way: neither side is the +/-0.5 cut.
  const SPORT_WITH_THREE_WAY_MONEYLINE = "soccer";
  const REASON_NO_RUNG = "no_rung";
  // Unabated flat-lines deep tails: Lions @ Bills 1H showed 19.8% on every
  // rung from 35.5 to 89 (2026-09-18). A rung whose fair repeats a
  // neighbour's is not a fair.
  const REASON_FLAT = "flat";
  const FLAT_TOLERANCE = 1e-9;

  function isHalfPoint(cut) {
    return Number.isInteger(cut * 2) && !Number.isInteger(cut);
  }

  // Map(eventId -> lines), one pass over the scanner state's lines, so a
  // ladder is built from its own event's lines only.
  function groupLinesByEvent(lines) {
    const byEvent = new Map();
    for (const line of lines) {
      if (line == null || line.eventId == null) continue;
      if (!byEvent.has(line.eventId)) byEvent.set(line.eventId, []);
      byEvent.get(line.eventId).push(line);
    }
    return byEvent;
  }

  function isThreeWayMoneyline(line) {
    const league = feed.LEAGUES[line.leagueId];
    return Boolean(league && league.sport === SPORT_WITH_THREE_WAY_MONEYLINE);
  }

  // Where one line sits: {axis, cut, winsAbove}, or null when it names no cut.
  function cutOfLine(line) {
    const winsAbove = line.sideIndex === SIDE_AWAY_OR_OVER;
    if (line.betTypeId === BET_TYPE_MONEYLINE) {
      if (isThreeWayMoneyline(line)) return null;
      return { axis: AXIS_MARGIN, cut: winsAbove ? MONEYLINE_CUT : -MONEYLINE_CUT, winsAbove };
    }
    if (typeof line.points !== "number" || !Number.isFinite(line.points)) return null;
    if (line.betTypeId === BET_TYPE_TOTAL) return { axis: AXIS_TOTAL, cut: line.points, winsAbove };
    if (line.betTypeId === BET_TYPE_SPREAD) return { axis: AXIS_MARGIN, cut: winsAbove ? -line.points : line.points, winsAbove };
    return null;
  }

  function median(values) {
    const sorted = values.slice().sort((a, b) => a - b);
    const middle = Math.floor(sorted.length / 2);
    return sorted.length % 2 === 1 ? sorted[middle] : (sorted[middle - 1] + sorted[middle]) / 2;
  }

  // The fair ladder of one (period, axis) from one event's lines. Per rung:
  // the median across books (and both sides) of the fair, as P(above the cut).
  function buildLadder(eventLines, { periodTypeId, axis }) {
    const samplesByCut = new Map();
    for (const line of eventLines || []) {
      if (line.periodTypeId !== periodTypeId) continue;
      if (typeof line.bacr !== "number" || Math.abs(line.bacr) < 100) continue;
      const position = cutOfLine(line);
      if (!position || position.axis !== axis || !isHalfPoint(position.cut)) continue;
      const probWin = kelly.americanToProb(line.bacr);
      if (!samplesByCut.has(position.cut)) samplesByCut.set(position.cut, []);
      samplesByCut.get(position.cut).push(position.winsAbove ? probWin : 1 - probWin);
    }
    const rungs = Array.from(samplesByCut.entries())
      .map(([cut, samples]) => [cut, median(samples)])
      .sort((a, b) => a[0] - b[0]);
    return { axis, periodTypeId, rungs };
  }

  // The two moneyline cuts legitimately share one fair when the game cannot
  // end level, so that pair is never "flat".
  function isMoneylinePair(ladder, cutA, cutB) {
    return ladder.axis === AXIS_MARGIN && Math.abs(cutA) === MONEYLINE_CUT && Math.abs(cutB) === MONEYLINE_CUT;
  }

  // P(result > cut) at one half-point rung: {prob}, or {reason} when the
  // ladder has no rung there or the rung repeats a neighbour's fair.
  function probAbove(ladder, cut) {
    const rungs = ladder && Array.isArray(ladder.rungs) ? ladder.rungs : [];
    const index = rungs.findIndex(([rungCut]) => rungCut === cut);
    if (index < 0) return { reason: REASON_NO_RUNG };
    const prob = rungs[index][1];
    for (const neighbour of [rungs[index - 1], rungs[index + 1]]) {
      if (!neighbour || isMoneylinePair(ladder, cut, neighbour[0])) continue;
      if (Math.abs(neighbour[1] - prob) <= FLAT_TOLERANCE) return { reason: REASON_FLAT };
    }
    return { prob };
  }

  // "1H" -> 2: a bet record names its period the way feed.PERIODS does. A
  // period Unabated has no id for (Kalshi's F5, I1) has no ladder.
  function periodTypeIdOf(period) {
    for (const [id, name] of Object.entries(feed.PERIODS)) if (name === period) return Number(id);
    return null;
  }

  const api = { AXIS_TOTAL, AXIS_MARGIN, REASON_NO_RUNG, REASON_FLAT, groupLinesByEvent, cutOfLine, buildLadder, probAbove, periodTypeIdOf };

  if (typeof module !== "undefined" && module.exports) {
    module.exports = api;
  } else {
    root.UnabatedLadder = api;
  }
})(typeof globalThis !== "undefined" ? globalThis : this);
