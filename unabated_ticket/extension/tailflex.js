// Tail flex for the Unabated Ticket panel: how far to trust Unabated's edge on
// a line far from the main number, used ONLY to rank lines. Pure: no DOM, no
// fetch, no chrome.* — loaded as a plain <script> in panel.html after feed.js
// (exposes globalThis.UnabatedTailFlex) and via require() in
// tests/tailflex.test.js.
//
// The "wing" model (Unabated Live 2026-09-11, 49:26: the fuselage is solid,
// the wings flex). Per line with Unabated fair p, decimal odds d, edge e:
//   z = qnorm(p), z_main = qnorm(fair of the same side at Unabated's main
//   number), dist_sd = |z - z_main|
//   sigma_z = sqrt(DZ_MAIN^2 + (c * dist_sd)^2)
//   sigma_e = d * dnorm(z) * sigma_z                  (edge uncertainty)
//   keep    = e^2 / (e^2 + sigma_e^2)                 (Baker & McHale 2013)
//   rank score = keep * e * stake                     (EV dollars after flex)
// The stake itself stays on Unabated's raw edge (user decision 2026-09-30):
// the keep factor only picks a card's best line and orders lines inside it.
//
// c is measured live per (league, period, bet type) from exchanges quoting
// BOTH sides of a rung: devig the pair (probit additive shift), take the
// probit gap to Unabated's fair, remove the gap the same exchanges show at
// the main number, divide by the rung's distance, and take the RMS (c is a
// standard deviation; a median of |miss| runs about a third low and floored
// a quarter of rungs at zero). Measured 2026-09-30, RMS with the top 1%
// dropped: NFL spreads 7.1%, totals 5.5%; CFB spreads 9.4%, totals 9.5%
// (the median had said 2.7 / 3.0 / 5.8 / 7.2). Exchanges are a measuring stick only;
// no fair is ever blended with an exchange price (user decision 2026-09-30).
//
// Inputs   scannerState (feed.js): lines {bookId, eventId, leagueId,
//          periodTypeId, betTypeId, sideIndex, points, price, sourceFormat,
//          sourcePrice, bacr, statusId, isAlt, ...}, books {id: {hasLiquidity}},
//          events {id: {eventStart}}.
// Outputs  measureTailFlex -> {markets, mainFairs}; rankOfRow -> {keep, score}.
//          Nothing here writes anywhere.

(function (root) {
  "use strict";

  const kelly = typeof module !== "undefined" && module.exports ? require("./kelly.js") : root.UnabatedKelly;
  const feed = typeof module !== "undefined" && module.exports ? require("./feed.js") : root.UnabatedFeed;

  // A one-cent error at the median, in probit units (the video's "fuselage").
  const DZ_MAIN = 0.01 / dnorm(0);
  // Used for a market with too few two-sided exchange rungs to measure c.
  const C_FALLBACK = 0.10;
  // An RMS of fewer ratios than this swings by points from slate to slate.
  const MIN_RUNGS = 100;
  // The largest ratios dropped before the RMS: one stale quote (live max 0.94
  // vs an RMS of 0.10 on CFB spreads) would otherwise move c by points.
  const TRIM_TOP_SHARE = 0.01;
  // ratio = excess gap / dist_sd; near the main number the division blows up.
  const MIN_DIST_SD = 0.3;
  // An exchange pair whose implied probabilities sum outside this is crossed
  // or one side is stale (the same envelope idea as unabated_edge's gate).
  const IMPLIED_SUM_MIN = 0.995;
  const IMPLIED_SUM_MAX = 1.20;
  const BET_TYPE_MONEYLINE = 1;
  const BET_TYPE_SPREAD = 2;
  const STATUS_ON_BOARD = 1;

  // ---- normal distribution -------------------------------------------------

  function dnorm(z) {
    return Math.exp(-0.5 * z * z) / Math.sqrt(2 * Math.PI);
  }

  // Abramowitz & Stegun 7.1.26 erf, |error| < 1.5e-7.
  function pnorm(z) {
    const x = Math.abs(z) / Math.SQRT2;
    const t = 1 / (1 + 0.3275911 * x);
    const poly = t * (0.254829592 + t * (-0.284496736 + t * (1.421413741 + t * (-1.453152027 + t * 1.061405429))));
    const erf = 1 - poly * Math.exp(-x * x);
    return z >= 0 ? 0.5 * (1 + erf) : 0.5 * (1 - erf);
  }

  // Acklam's inverse normal CDF, relative error < 1.2e-9.
  const QNORM_A = [-3.969683028665376e+01, 2.209460984245205e+02, -2.759285104469687e+02, 1.383577518672690e+02, -3.066479806614716e+01, 2.506628277459239e+00];
  const QNORM_B = [-5.447609879822406e+01, 1.615858368580409e+02, -1.556989798598866e+02, 6.680131188771972e+01, -1.328068155288572e+01];
  const QNORM_C = [-7.784894002430293e-03, -3.223964580411365e-01, -2.400758277161838e+00, -2.549732539343734e+00, 4.374664141464968e+00, 2.938163982698783e+00];
  const QNORM_D = [7.784695709041462e-03, 3.224671290700398e-01, 2.445134137142996e+00, 3.754408661907416e+00];
  const QNORM_P_LOW = 0.02425;

  function qnorm(p) {
    if (!(p > 0 && p < 1)) throw new Error(`qnorm: expected a probability strictly between 0 and 1, got ${p}`);
    if (p > 1 - QNORM_P_LOW) return -qnorm(1 - p);
    if (p < QNORM_P_LOW) {
      const [c, d] = [QNORM_C, QNORM_D];
      const q = Math.sqrt(-2 * Math.log(p));
      return (((((c[0] * q + c[1]) * q + c[2]) * q + c[3]) * q + c[4]) * q + c[5]) / ((((d[0] * q + d[1]) * q + d[2]) * q + d[3]) * q + 1);
    }
    const [a, b] = [QNORM_A, QNORM_B];
    const q = p - 0.5;
    const r = q * q;
    return (((((a[0] * r + a[1]) * r + a[2]) * r + a[3]) * r + a[4]) * r + a[5]) * q / (((((b[0] * r + b[1]) * r + b[2]) * r + b[3]) * r + b[4]) * r + 1);
  }

  // ---- devig ---------------------------------------------------------------

  // A book's implied probability for one side (kelly.bookProbOf: the exact
  // exchange figure when published, else the American price), or null.
  function impliedProbOf(line) {
    try {
      return kelly.bookProbOf({ bookPrice: line.price, sourceFormat: line.sourceFormat, sourcePrice: line.sourcePrice });
    } catch (_error) {
      return null;
    }
  }

  // Probit additive shift (as unabated_edge / Tools.R): find k with
  // pnorm(z0 - k) + pnorm(z1 - k) = 1. For two outcomes that is k = (z0 + z1)/2,
  // so side 0's fair in probit is (z0 - z1)/2.
  function devigProbitTwoWayZ(impliedProb0, impliedProb1) {
    return (qnorm(impliedProb0) - qnorm(impliedProb1)) / 2;
  }

  function devigProbitTwoWay(impliedProb0, impliedProb1) {
    return pnorm(devigProbitTwoWayZ(impliedProb0, impliedProb1));
  }

  // ---- keys ----------------------------------------------------------------

  function marketKeyOf({ leagueId, periodTypeId, betTypeId }) {
    return `${leagueId}:${periodTypeId}:${betTypeId}`;
  }

  function sideKeyOf({ eventId, periodTypeId, betTypeId, sideIndex }) {
    return `${eventId}:${periodTypeId}:${betTypeId}:${sideIndex}`;
  }

  // The other side of the same bet: a spread's -x is the other team's +x, a
  // total's other side sits on the same number.
  function mirroredPoints(line) {
    return line.betTypeId === BET_TYPE_SPREAD ? -line.points : line.points;
  }

  // ---- Unabated's main number ----------------------------------------------

  function isUnabatedMainLine(line) {
    return line.bookId === feed.UNABATED_LINE_BOOK_ID && !line.isAlt;
  }

  // Record Unabated's own main line (book 49) as {points, fairProb} under its
  // (event, period, bet type, side). Its price is its fair.
  function noteUnabatedMainFair(fairs, line) {
    const fairAmerican = line.bacr ?? line.price;
    if (!kelly.isAmericanPrice(fairAmerican) || feed.isClampedFair(fairAmerican)) return;
    fairs.set(sideKeyOf(line), { points: line.points, fairProb: kelly.americanToProb(fairAmerican) });
  }

  // ---- measuring c -----------------------------------------------------------

  function isFresh(line, now, maxLineAgeMs) {
    if (line.statusId !== STATUS_ON_BOARD) return false;
    if (maxLineAgeMs == null) return true;
    const changedMs = feed.lineChangedMs(line);
    return changedMs != null && now - changedMs <= maxLineAgeMs;
  }

  function isMeasurableExchangeLine(line, state, now, maxLineAgeMs) {
    const book = state.books[line.bookId];
    if (!book || !book.hasLiquidity || line.betTypeId === BET_TYPE_MONEYLINE || line.points == null) return false;
    const event = state.events[line.eventId];
    if (!event || event.eventStart == null || event.eventStart <= now) return false;
    return isFresh(line, now, maxLineAgeMs);
  }

  // Every exchange rung quoted fresh on both sides at the mirrored number of
  // a game not yet started, as {marketKey, gap, distSd, onMain}: gap = probit
  // distance between Unabated's side-0 fair and the exchange's devigged side-0
  // fair. The probit is symmetric, so side 0 measures the tail side as well.
  // One pass over the state's lines (~190k on an NFL + CFB Saturday) collects
  // both Unabated's main fairs and the exchange lines.
  function twoSidedExchangeRungs(state, now, maxLineAgeMs) {
    const mainFairs = new Map();
    const exchangeLines = new Map();
    for (const key in state.lines) {
      const line = state.lines[key];
      if (isUnabatedMainLine(line)) noteUnabatedMainFair(mainFairs, line);
      else if (isMeasurableExchangeLine(line, state, now, maxLineAgeMs)) exchangeLines.set(`${line.bookId}|${sideKeyOf(line)}|${line.points}`, line);
    }
    const rungs = [];
    for (const side0 of exchangeLines.values()) {
      if (side0.sideIndex !== 0) continue;
      const side1 = exchangeLines.get(`${side0.bookId}|${sideKeyOf({ ...side0, sideIndex: 1 })}|${mirroredPoints(side0)}`);
      const main = mainFairs.get(sideKeyOf(side0));
      if (!side1 || !main || !kelly.isAmericanPrice(side0.bacr) || feed.isClampedFair(side0.bacr)) continue;
      const implied0 = impliedProbOf(side0);
      const implied1 = impliedProbOf(side1);
      if (implied0 == null || implied1 == null) continue;
      const impliedSum = implied0 + implied1;
      if (impliedSum < IMPLIED_SUM_MIN || impliedSum > IMPLIED_SUM_MAX) continue;
      const unabatedZ = qnorm(kelly.americanToProb(side0.bacr));
      rungs.push({
        marketKey: marketKeyOf(side0),
        gap: Math.abs(unabatedZ - devigProbitTwoWayZ(implied0, implied1)),
        distSd: Math.abs(unabatedZ - qnorm(main.fairProb)),
        onMain: side0.points === main.points,
      });
    }
    return { rungs, mainFairs };
  }

  function rootMeanSquare(values) {
    return Math.sqrt(values.reduce((sum, value) => sum + value * value, 0) / values.length);
  }

  // RMS of the values after dropping the largest TRIM_TOP_SHARE of them.
  function trimmedRootMeanSquare(values) {
    const sorted = [...values].sort((a, b) => a - b);
    const keepCount = sorted.length - Math.floor(sorted.length * TRIM_TOP_SHARE);
    return rootMeanSquare(sorted.slice(0, keepCount));
  }

  // c for one market from its rungs: the gap at the main number is the
  // exchanges' own noise (baseline, RMS); what is left further out, per SD
  // of distance, is the flex. Fewer than MIN_RUNGS usable ratios -> fallback.
  function cOfRungs(rungs) {
    const mainGaps = rungs.filter((rung) => rung.onMain).map((rung) => rung.gap);
    const baseline = mainGaps.length ? rootMeanSquare(mainGaps) : 0;
    const ratios = rungs.filter((rung) => !rung.onMain && rung.distSd >= MIN_DIST_SD)
      .map((rung) => Math.sqrt(Math.max(0, rung.gap * rung.gap - baseline * baseline)) / rung.distSd);
    const measured = ratios.length >= MIN_RUNGS;
    return { c: measured ? trimmedRootMeanSquare(ratios) : C_FALLBACK, measured, rungCount: ratios.length, baseline };
  }

  // {markets: {marketKey: {c, measured, rungCount, baseline}}, mainFairs}.
  // A market absent from `markets` has no measurable rung at all: cOf falls
  // back for it. opts {now, maxLineAgeMs} — the panel's own line-age gate.
  function measureTailFlex(state, opts) {
    const options = opts || {};
    const now = typeof options.now === "number" ? options.now : Date.now();
    const maxLineAgeMs = typeof options.maxLineAgeMs === "number" && options.maxLineAgeMs > 0 ? options.maxLineAgeMs : null;
    const { rungs, mainFairs } = twoSidedExchangeRungs(state, now, maxLineAgeMs);
    const byMarket = new Map();
    for (const rung of rungs) {
      if (!byMarket.has(rung.marketKey)) byMarket.set(rung.marketKey, []);
      byMarket.get(rung.marketKey).push(rung);
    }
    const markets = {};
    for (const [marketKey, rungs] of byMarket) markets[marketKey] = cOfRungs(rungs);
    return { markets, mainFairs };
  }

  // The c in use for a row's market: measured, else C_FALLBACK.
  function cOf(measurement, row) {
    const market = measurement ? measurement.markets[marketKeyOf(row)] : null;
    return market ? market.c : C_FALLBACK;
  }

  function isMeasured(measurement, row) {
    const market = measurement ? measurement.markets[marketKeyOf(row)] : null;
    return Boolean(market && market.measured);
  }

  // ---- the rank score --------------------------------------------------------

  // keep in (0, 1]: the share of the edge that survives its own uncertainty.
  // zMain null = no ladder (a moneyline): only the one-cent fuselage applies.
  function keepFactor({ edge, decimal, fairProb, zMain, c }) {
    const z = qnorm(fairProb);
    const distSd = zMain == null ? 0 : Math.abs(z - zMain);
    const sigmaZ = Math.sqrt(DZ_MAIN * DZ_MAIN + (c * distSd) * (c * distSd));
    const sigmaE = decimal * dnorm(z) * sigmaZ;
    return (edge * edge) / (edge * edge + sigmaE * sigmaE);
  }

  // A listed row's {keep, score}, score = keep x edge x stake (EV dollars
  // after flex); null when the row has no stake or no positive edge. The fair
  // is read back from the edge (Unabated's ge = p x decimal - 1, exactly).
  // A side Unabated has no main line for is measured from even money, where
  // a main number sits.
  function rankOfRow(row, measurement) {
    if (typeof row.stake !== "number" || row.edgePct == null || row.edgePct <= 0) return null;
    if (!kelly.isAmericanPrice(row.price)) return null;
    const edge = row.edgePct / 100;
    const decimal = kelly.americanToDecimal(row.price);
    const fairProb = (1 + edge) / decimal;
    if (!(fairProb > 0 && fairProb < 1)) return null;
    let zMain = null;
    if (row.betTypeId !== BET_TYPE_MONEYLINE) {
      const main = measurement ? measurement.mainFairs.get(sideKeyOf(row)) : null;
      zMain = main ? qnorm(main.fairProb) : 0;
    }
    const keep = keepFactor({ edge, decimal, fairProb, zMain, c: cOf(measurement, row) });
    return { keep, score: keep * edge * row.stake };
  }

  const api = {
    DZ_MAIN, C_FALLBACK, MIN_RUNGS, TRIM_TOP_SHARE, MIN_DIST_SD, IMPLIED_SUM_MIN, IMPLIED_SUM_MAX,
    dnorm, pnorm, qnorm, impliedProbOf, devigProbitTwoWay, marketKeyOf,
    cOfRungs, measureTailFlex, cOf, isMeasured, keepFactor, rankOfRow,
  };

  if (typeof module !== "undefined" && module.exports) {
    module.exports = api;
  } else {
    root.UnabatedTailFlex = api;
  }
})(typeof globalThis !== "undefined" ? globalThis : this);
