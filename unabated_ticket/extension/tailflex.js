// Tail flex for the Unabated Ticket panel: how much of Unabated's claimed edge
// on a line to believe, less the farther the line sits from the main number.
// The shrunk edge picks WHICH rung of a ladder to bet (the rank score) and
// whether a rung is bet at all (the Min shrunk edge floor, a setting); the
// stake stays Kelly on Unabated's raw edge (issue #208, Cal 2026-10-23: the
// rung choice was the leak, and a chosen rung is bet on Unabated's own
// number). Pure: no DOM, no fetch, no chrome.* — loaded as a plain
// <script> in panel.html after feed.js (exposes globalThis.UnabatedTailFlex)
// and via require() in tests/tailflex.test.js.
//
// The "wing" model (Unabated Live 2026-09-11, 49:26: the fuselage is solid,
// the wings flex). Per line with Unabated fair p, decimal odds d, edge e:
//   z = qnorm(p), z_main = qnorm(fair of the same side at Unabated's main
//   number), dist_sd = |z - z_main|
//   sigma_z = c * dist_sd          (zero at the main number: Cal 2026-10-23,
//                                   Unabated's main fair is taken as right)
//   sigma_e = d * dnorm(z) * sigma_z                  (edge uncertainty, EV units)
//   shrink  = tau^2 / (tau^2 + sigma_e^2)             (empirical Bayes)
//   shrunk edge = e * shrink                          (the floor reads it)
//   rank score  = shrunk edge * stake                 (EV dollars believed on the
//                                                      raw-edge stake; picks a card's line)
// tau is how much REAL edges vary from line to line (a standard deviation in
// EV units); sigma_e is how wrong THIS line's claimed edge can be. The
// earlier keep = e^2 / (e^2 + sigma_e^2) (Baker & McHale 2013) let a big
// claimed edge vouch for itself, so the deepest rungs — the noisiest, and
// the ones a scan of thousands of rungs lists for their noise (the
// optimizer's curse, Smith & Winkler 2006) — won the card. The shrink does
// not depend on the claimed edge: the main number wins a card unless a deep
// rung's claimed edge is many times larger.
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
// tau is measured per (league, period, bet type) from the rungs where TWO OR
// MORE exchanges quote both sides. Each exchange's own devigged price (its
// mid) at such a rung is priced twice: against Unabated's fair (x = p_U * d -
// 1) and against the OTHER exchanges' devigged fair (y = p_O * d - 1; the
// exchange being priced never votes on itself). Both read the same real edge with their own
// errors, so their covariance across rungs is the variance of real edges:
// tau = sqrt(cov(x, y)). One refresh holds too few such rungs, so the caller
// pools them over a rolling window (createTauPool, 24 h; the server runner
// keeps it across restarts through the bets service) and the measurement
// takes the pooled tau (withPooledTau). Under MIN_RUNGS pooled rungs, or a
// covariance at or below zero, tau falls back to TAU_FALLBACK (1.85%, the fit
// on Cal's bets 9/25 to 10/4 against the exchange close; issue #208) and the
// header says so.
//
// Inputs   scannerState (feed.js): lines {bookId, eventId, leagueId,
//          periodTypeId, betTypeId, sideIndex, points, price, sourceFormat,
//          sourcePrice, bacr, statusId, isAlt, ...}, books {id: {hasLiquidity}},
//          events {id: {eventStart}}.
// Outputs  measureTailFlex -> {markets, mainFairs, tauSamples}; shrinkOfRow ->
//          {shrink, edgePct, tau, c}; rankScoreOf -> the believed growth of
//          a sized row's stake, in dollars.
//          Nothing here writes anywhere.

(function (root) {
  "use strict";

  const kelly = typeof module !== "undefined" && module.exports ? require("./kelly.js") : root.UnabatedKelly;
  const feed = typeof module !== "undefined" && module.exports ? require("./feed.js") : root.UnabatedFeed;

  // Used for a market with too few two-sided exchange rungs to measure c.
  const C_FALLBACK = 0.10;
  // An RMS of fewer ratios than this swings by points from slate to slate;
  // a covariance over fewer distinct rungs than this swings more (a rung's
  // two sides and its exchanges are readings of one real edge, so tau
  // counts rungs, never samples).
  const MIN_RUNGS = 100;
  // tau until a market's pool holds MIN_RUNGS distinct rungs with a positive
  // covariance: 1.85%, fitted on Cal's 255 bets of 9/25 to 10/4 that had both
  // a saved fill fair and a real exchange close (issue #208; 1% to 3% fit
  // about as well — a soft number, which is why tau is measured).
  const TAU_FALLBACK = 0.0185;
  // How long a two-sided rung stays in the tau pool after it was last seen.
  // One refresh holds about a hundred such rungs per market on a game day;
  // a day's distinct rungs hold thousands (issue #208, open question 4).
  const TAU_POOL_WINDOW_MS = 24 * 3600 * 1000;
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
  // a game not yet started, as {marketKey, rungKey, bookId, gap, distSd,
  // onMain, unabatedZ, exchangeZ, implied0, implied1}: gap = probit distance
  // between Unabated's side-0 fair and the exchange's devigged side-0 fair
  // (exchangeZ); rungKey names the rung without the book, so the rungs of two
  // exchanges at one number meet (tau). The probit is symmetric, so side 0
  // measures the tail side as well. One pass over the state's lines (~190k on
  // an NFL + CFB Saturday) collects both Unabated's main fairs and the
  // exchange lines.
  function twoSidedExchangeRungs(state, now, maxLineAgeMs) {
    const mainFairs = new Map();
    const exchangeLines = new Map();
    for (const key in state.lines) {
      const line = state.lines[key];
      if (isUnabatedMainLine(line)) noteUnabatedMainFair(mainFairs, line);
      else if (isMeasurableExchangeLine(line, state, now, maxLineAgeMs)) exchangeLines.set(`${line.bookId}|${sideKeyOf(line)}|${line.points}`, line);
    }
    // A flat-lined Unabated tail is no fair to measure against (feed.flatFairRungKeys).
    const flatRungKeys = feed.flatFairRungKeys(Object.values(state.lines));
    const rungs = [];
    for (const side0 of exchangeLines.values()) {
      if (side0.sideIndex !== 0 || flatRungKeys.has(feed.rungKeyOf(side0))) continue;
      const side1 = exchangeLines.get(`${side0.bookId}|${sideKeyOf({ ...side0, sideIndex: 1 })}|${mirroredPoints(side0)}`);
      const main = mainFairs.get(sideKeyOf(side0));
      if (!side1 || !main || !kelly.isAmericanPrice(side0.bacr) || feed.isClampedFair(side0.bacr)) continue;
      const implied0 = impliedProbOf(side0);
      const implied1 = impliedProbOf(side1);
      if (implied0 == null || implied1 == null) continue;
      const impliedSum = implied0 + implied1;
      if (impliedSum < IMPLIED_SUM_MIN || impliedSum > IMPLIED_SUM_MAX) continue;
      const unabatedZ = qnorm(kelly.americanToProb(side0.bacr));
      const exchangeZ = devigProbitTwoWayZ(implied0, implied1);
      rungs.push({
        marketKey: marketKeyOf(side0),
        rungKey: `${sideKeyOf(side0)}:${side0.points}`,
        bookId: side0.bookId,
        gap: Math.abs(unabatedZ - exchangeZ),
        distSd: Math.abs(unabatedZ - qnorm(main.fairProb)),
        onMain: side0.points === main.points,
        unabatedZ, exchangeZ, implied0, implied1,
      });
    }
    return { rungs, mainFairs };
  }

  // ---- measuring tau ---------------------------------------------------------

  // The tau samples of one refresh: for every rung two or more exchanges
  // quote both sides, each exchange's devigged price on each side (its mid,
  // not its ask: the ask's hold varies with the rung's probability and would
  // read as edge variance both reads share) priced against Unabated's fair
  // (edgeUnabated) and against the other exchanges' devigged fair
  // (edgeExchanges, the mean of their probit fairs; the exchange being
  // priced is left out). Both are edges in EV units at that mid; tau is
  // about how they vary together across rungs, not their level.
  // [{key, marketKey, edgeUnabated, edgeExchanges}], key = "<bookId>|<rungKey>|<sideIndex>".
  function tauSamplesOfRungs(rungs) {
    const byRung = new Map();
    for (const rung of rungs) {
      if (!byRung.has(rung.rungKey)) byRung.set(rung.rungKey, []);
      byRung.get(rung.rungKey).push(rung);
    }
    const samples = [];
    for (const atRung of byRung.values()) {
      if (new Set(atRung.map((rung) => rung.bookId)).size < 2) continue;
      for (const rung of atRung) {
        const others = atRung.filter((other) => other.bookId !== rung.bookId);
        const othersZ = others.reduce((sum, other) => sum + other.exchangeZ, 0) / others.length;
        const unabatedProb = pnorm(rung.unabatedZ);
        const othersProb = pnorm(othersZ);
        const decimal0 = 1 / pnorm(rung.exchangeZ);
        const decimal1 = 1 / pnorm(-rung.exchangeZ);
        samples.push({
          key: `${rung.bookId}|${rung.rungKey}|0`, marketKey: rung.marketKey,
          edgeUnabated: unabatedProb * decimal0 - 1, edgeExchanges: othersProb * decimal0 - 1,
        });
        samples.push({
          key: `${rung.bookId}|${rung.rungKey}|1`, marketKey: rung.marketKey,
          edgeUnabated: (1 - unabatedProb) * decimal1 - 1, edgeExchanges: (1 - othersProb) * decimal1 - 1,
        });
      }
    }
    return samples;
  }

  // The rung a tau sample was read at, off its key "<bookId>|<rungKey>|<sideIndex>"
  // (tauSamplesOfRungs). Any other shape is a bug upstream, never a rung.
  function rungOfSampleKey(key) {
    const parts = typeof key === "string" ? key.split("|") : [];
    if (parts.length !== 3) throw new Error(`tau sample key must be "<bookId>|<rungKey>|<sideIndex>", got ${JSON.stringify(key)}`);
    return parts[1];
  }

  // tau for one market from its samples: the covariance of the two edge
  // reads after dropping the TRIM_TOP_SHARE largest (by either read: a stale
  // quote at a tail reads as a 50% edge and would own the covariance).
  // Fewer than MIN_RUNGS distinct rungs among the kept samples (a rung's two
  // sides at two exchanges are four samples of one real edge, so samples
  // would reach 100 at 25 rungs), or a covariance at or below zero (the two
  // reads do not agree on where the edges are), fall back; `covariance` is
  // reported either way so the fallback can be seen for what it is, and
  // `sampleCount` beside `rungCount` so the pool's depth can be.
  function tauOfSamples(samples) {
    const sorted = [...samples].sort((a, b) => Math.max(Math.abs(a.edgeUnabated), Math.abs(a.edgeExchanges)) - Math.max(Math.abs(b.edgeUnabated), Math.abs(b.edgeExchanges)));
    const kept = sorted.slice(0, sorted.length - Math.floor(sorted.length * TRIM_TOP_SHARE));
    const sampleCount = kept.length;
    const rungCount = new Set(kept.map((sample) => rungOfSampleKey(sample.key))).size;
    if (rungCount < MIN_RUNGS) return { tau: TAU_FALLBACK, measured: false, rungCount, sampleCount, covariance: null };
    const meanX = kept.reduce((sum, sample) => sum + sample.edgeUnabated, 0) / sampleCount;
    const meanY = kept.reduce((sum, sample) => sum + sample.edgeExchanges, 0) / sampleCount;
    const covariance = kept.reduce((sum, sample) => sum + (sample.edgeUnabated - meanX) * (sample.edgeExchanges - meanY), 0) / sampleCount;
    if (!(covariance > 0)) return { tau: TAU_FALLBACK, measured: false, rungCount, sampleCount, covariance };
    return { tau: Math.sqrt(covariance), measured: true, rungCount, sampleCount, covariance };
  }

  // {marketKey: tauOfSamples(that market's samples)}.
  function tauByMarketOf(samples) {
    const byMarket = new Map();
    for (const sample of samples) {
      if (!byMarket.has(sample.marketKey)) byMarket.set(sample.marketKey, []);
      byMarket.get(sample.marketKey).push(sample);
    }
    const markets = {};
    for (const [marketKey, marketSamples] of byMarket) markets[marketKey] = tauOfSamples(marketSamples);
    return markets;
  }

  // A rolling pool of tau samples keyed by sample key: the latest reading of
  // each rung side, kept TAU_POOL_WINDOW_MS after it was last seen. Held by
  // the caller (one per runner), never by this module.
  //   absorb(samples, seenAt)  -> the samples whose reading changed or whose
  //                               stored seenAt is older than `heartbeatMs`
  //                               (what the runner has to persist)
  //   load(rows)               -> restore rows [{key, marketKey, edgeUnabated,
  //                               edgeExchanges, seenAt (ms)}] (a persisted
  //                               copy; a held newer reading wins)
  //   measure(now)             -> tauByMarketOf(the samples still inside the window)
  //   size()                   -> samples held
  function createTauPool({ windowMs = TAU_POOL_WINDOW_MS, heartbeatMs = 3600 * 1000 } = {}) {
    const samples = new Map();
    function evict(now) {
      for (const [key, sample] of samples) {
        if (now - sample.seenAt > windowMs) samples.delete(key);
      }
    }
    function absorb(fresh, seenAt) {
      const changed = [];
      for (const sample of fresh) {
        const held = samples.get(sample.key);
        const same = held && held.edgeUnabated === sample.edgeUnabated && held.edgeExchanges === sample.edgeExchanges;
        const stored = { ...sample, seenAt, persistedAt: same ? held.persistedAt : null };
        if (!same || seenAt - held.persistedAt >= heartbeatMs) {
          stored.persistedAt = seenAt;
          changed.push({ key: sample.key, marketKey: sample.marketKey, edgeUnabated: sample.edgeUnabated, edgeExchanges: sample.edgeExchanges, seenAt });
        }
        samples.set(sample.key, stored);
      }
      evict(seenAt);
      return changed;
    }
    // A persisted row never replaces a fresher reading this process already
    // absorbed (the service comes back after the pool filled).
    function load(rows) {
      for (const row of rows) {
        const held = samples.get(row.key);
        if (held && held.seenAt > row.seenAt) continue;
        samples.set(row.key, { key: row.key, marketKey: row.marketKey, edgeUnabated: row.edgeUnabated, edgeExchanges: row.edgeExchanges, seenAt: row.seenAt, persistedAt: row.seenAt });
      }
    }
    function measure(now) {
      evict(now);
      return tauByMarketOf(Array.from(samples.values()));
    }
    return { absorb, load, measure, size: () => samples.size };
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

  // {markets: {marketKey: {c, measured, rungCount, baseline}}, mainFairs,
  // tauSamples}. A market absent from `markets` has no measurable rung at
  // all: cOf falls back for it. tau is never read off one refresh: the caller
  // absorbs `tauSamples` into a pool (or reads a served pool) and adds the
  // tau fields through withPooledTau; without that step tauOf falls back.
  // opts {now, maxLineAgeMs} — the panel's own line-age gate.
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
    for (const [marketKey, marketRungs] of byMarket) markets[marketKey] = cOfRungs(marketRungs);
    return { markets, mainFairs, tauSamples: tauSamplesOfRungs(rungs) };
  }

  function tauFields(measured) {
    return { tau: measured.tau, tauMeasured: measured.measured, tauRungCount: measured.rungCount, tauCovariance: measured.covariance };
  }

  // The measurement with each market's tau fields (tau, tauMeasured,
  // tauRungCount, tauCovariance) from the pool (tauPool.measure(now), or a
  // served copy of it: {marketKey: {tau, measured, rungCount, covariance}}).
  // A market the pool has no samples for gets the fallback.
  function withPooledTau(measurement, pooled) {
    const markets = {};
    for (const [marketKey, market] of Object.entries(measurement.markets)) {
      markets[marketKey] = { ...market, ...tauFields(pooled[marketKey] || tauOfSamples([])) };
    }
    // A market the pool measured but this refresh has no rung for still sizes
    // rows (c falls back; tau is the pool's).
    for (const [marketKey, measured] of Object.entries(pooled)) {
      if (!markets[marketKey]) markets[marketKey] = { c: C_FALLBACK, measured: false, rungCount: 0, baseline: 0, ...tauFields(measured) };
    }
    return { ...measurement, markets };
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

  // The tau in use for a row's market: pooled and measured, else TAU_FALLBACK.
  function tauOf(measurement, row) {
    const market = measurement ? measurement.markets[marketKeyOf(row)] : null;
    return market && typeof market.tau === "number" ? market.tau : TAU_FALLBACK;
  }

  function isTauMeasured(measurement, row) {
    const market = measurement ? measurement.markets[marketKeyOf(row)] : null;
    return Boolean(market && market.tauMeasured);
  }

  // ---- the shrink --------------------------------------------------------------

  // sigma_e: how wrong the claimed edge can be, in EV units: nothing at the
  // main number (Unabated's main fair is taken as right, so a main line and a
  // moneyline keep their whole edge), growing with distance from it at c.
  function edgeUncertainty({ decimal, fairProb, zMain, c }) {
    const z = qnorm(fairProb);
    const distSd = zMain == null ? 0 : Math.abs(z - zMain);
    return decimal * dnorm(z) * c * distSd;
  }

  // {shrink, sigmaE}: shrink in (0, 1] is the share of a claimed edge to
  // believe, given how much real edges vary (tau) and how wrong this one can
  // be (sigma_e). It does not depend on the claimed edge itself.
  function shrinkFactor({ decimal, fairProb, zMain, c, tau }) {
    const sigmaE = edgeUncertainty({ decimal, fairProb, zMain, c });
    return { shrink: (tau * tau) / (tau * tau + sigmaE * sigmaE), sigmaE };
  }

  // A row's shrink and shrunk edge: {shrink, edgePct (shrunk, in percent),
  // tau, c}, or null when the row has no positive edge or no American price. The fair is read back from the edge (Unabated's ge = p x decimal -
  // 1, exactly). A side Unabated has no main line for is measured from even
  // money, where a main number sits. The raw edge stays on the row: the
  // listing floor (Min edge, the bot's 3%) reads it; the stake and the rank
  // read this.
  function shrinkOfRow(row, measurement) {
    if (row.edgePct == null || row.edgePct <= 0 || !kelly.isAmericanPrice(row.price)) return null;
    const edge = row.edgePct / 100;
    const decimal = kelly.americanToDecimal(row.price);
    const fairProb = (1 + edge) / decimal;
    if (!(fairProb > 0 && fairProb < 1)) return null;
    let zMain = null;
    if (row.betTypeId !== BET_TYPE_MONEYLINE) {
      const main = measurement ? measurement.mainFairs.get(sideKeyOf(row)) : null;
      zMain = main ? qnorm(main.fairProb) : 0;
    }
    const c = cOf(measurement, row);
    const tau = tauOf(measurement, row);
    const { shrink } = shrinkFactor({ decimal, fairProb, zMain, c, tau });
    return { shrink, edgePct: row.edgePct * shrink, tau, c };
  }

  // A sized row's rank score: the expected log growth, in dollars, of the
  // stake that would be placed (Kelly on the raw edge; edgerows.withStakeAndRank
  // sets it), believed at the shrunk edge — the second-order Kelly expansion
  //   stake x shrunkEdge - stake^2 x b / (2 x bankroll),   b = decimal - 1.
  // The first term is the EV believed; the second is what the stake's size
  // costs in variance, so a big stake on a little-trusted claim cannot win a
  // card through its size (shrunk edge x stake alone let a 10% claim at
  // shrink 0.15 outrank a 4% claim at 0.80; the 2026-10-23 review). It goes
  // negative once the stake is past twice the Kelly the shrunk edge supports
  // (shrink under multiplier / 2); nothing gates on that (Cal, 2026-10-23: no
  // minimum trust) — it only picks a card's best line and orders the lines
  // inside it. null when the row has no stake or no shrunk edge; a bankroll
  // that is not a positive number is a bug.
  function rankScoreOf(row, bankroll) {
    if (!(typeof bankroll === "number" && bankroll > 0)) throw new Error(`rankScoreOf: expected a positive bankroll, got ${bankroll}`);
    if (typeof row.stake !== "number" || typeof row.shrunkEdgePct !== "number") return null;
    const netOdds = kelly.americanToDecimal(row.price) - 1;
    const believedEdge = row.shrunkEdgePct / 100;
    return row.stake * believedEdge - (row.stake * row.stake * netOdds) / (2 * bankroll);
  }

  const api = {
    C_FALLBACK, TAU_FALLBACK, TAU_POOL_WINDOW_MS, MIN_RUNGS, TRIM_TOP_SHARE, MIN_DIST_SD, IMPLIED_SUM_MIN, IMPLIED_SUM_MAX,
    dnorm, pnorm, qnorm, impliedProbOf, devigProbitTwoWay, marketKeyOf,
    cOfRungs, tauSamplesOfRungs, tauOfSamples, tauByMarketOf, createTauPool, measureTailFlex, withPooledTau,
    cOf, isMeasured, tauOf, isTauMeasured, edgeUncertainty, shrinkFactor, shrinkOfRow, rankScoreOf,
  };

  if (typeof module !== "undefined" && module.exports) {
    module.exports = api;
  } else {
    root.UnabatedTailFlex = api;
  }
})(typeof globalThis !== "undefined" ? globalThis : this);
