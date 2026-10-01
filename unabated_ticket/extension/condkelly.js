// Conditional Kelly for the Unabated Ticket panel (#130): the stake for the
// next bet GIVEN the bets already held on the same market of the same game.
// Pure functions, no DOM, no feed or bet-record knowledge — it sees cuts,
// directions, dollars and probabilities.
//
// The rule: held bets are fixed; choose x >= 0 that maximizes
//   sum over outcome rows of  prob * ln(1 + pnl / K)
// where K = bankroll * multiplier (the scale the held bets were sized on).
// With no held bets the maximum is kelly.kellyStakeFromEdge, to the cent.
//
// Inputs (solveStake)
//   kellyBankroll  K, dollars
//   candidate      { group, cut, direction, netOdds, prob } — the new bet
//   held           [{ group, cut, direction, stake, toWin }] — all on the
//                  candidate's market axis (the caller filters)
//   tickets        [{ group, cut, direction, stake, toWin, others }] — open
//                  teasers with a leg on this axis (2026-09-30). The leg wins
//                  past its half-point `cut` (a teaser push loses); the
//                  ticket pays +toWin only when that leg AND every leg in
//                  `others` win, else -stake. An other leg { factor, cut,
//                  direction } sits on another game, named by `factor` and
//                  independent of this one. Optional.
//   factors        { [factor]: [[cut, probAbove], ...] } — P(result > cut) on
//                  each other game, at every cut its legs need. Optional.
//   ladders        { [group]: [[cut, probAbove], ...] } — P(result > cut) at
//                  half-point cuts, for every cut a held bet or ticket needs
// A `group` is one period of the axis ("FG", "1H"). `cut` + `direction`
// ("above" | "below") say when a bet wins: Over 52.5 = above 52.5; on the
// margin axis (away minus home) an away bet at points a = above -a, a home
// bet at points h = below h, a moneyline = above +0.5 / below -0.5.
// A whole-number cut k pushes between k-0.5 and k+0.5.
// A ticket's other legs make its P&L on a row of this game a gamble, and
// under the log objective a gamble cannot be replaced by its average: the
// other games are enumerated, grouped by which tickets they leave alive, and
// the score is the chance-weighted sum of the growth over those groups.
// Output { stake, reason }: `reason` is set (and stake null) when the calc is
// declined. Nothing here writes anywhere. teaser.js reads cutsNeeded and
// resultAt to put a straight bet held under a teaser leg on its game's rows.
//
// Loaded two ways: as a plain <script> in panel.html (exposes
// globalThis.UnabatedCondKelly), before teaser.js, and via require() in
// tests/condkelly.test.js.

(function (root) {
  "use strict";

  // Unabated's fair is a whole American price, so neighbouring rungs can
  // cross by a hair; a slice this negative is rounding, past it the ladder
  // is wrong and the calc is declined.
  const MAX_NEGATIVE_SLICE = 0.005;
  const REASON_NOT_MONOTONE = "ladder not monotone";
  const REASON_HELD_RISK = "held risk exceeds the Kelly bankroll";
  // A clamped ladder can leave the new bet no losing row; there is no
  // maximum to find, so the calc is declined rather than sized to infinity.
  const REASON_NO_LOSING_ROW = "ladder leaves the new bet no losing outcome";
  // The joint outcomes of the other games under the tickets one calc may
  // enumerate. A 4-team ticket rides on three other games; 9/27's eight
  // tickets put seven under one row (128 outcomes).
  const MAX_OTHER_GAME_OUTCOMES = 1 << 16;
  const REASON_TOO_MANY_GAMES = "the open teasers ride on too many other games to size";
  const GOLDEN_RATIO = (Math.sqrt(5) - 1) / 2;
  const SEARCH_MAX_ITERATIONS = 200;
  const SEARCH_TOLERANCE_DOLLARS = 1e-7;
  // The probe for "is the slope at 0 positive", as a fraction of K.
  const SLOPE_PROBE_FRACTION = 1e-7;
  // Stay strictly inside the stake that would take a row's wealth to zero.
  const UPPER_BOUND_SHRINK = 1 - 1e-9;
  const PROB_EPSILON = 1e-12;
  const HALF_POINT = 0.5;

  function assertFinite(value, label) {
    if (typeof value !== "number" || !Number.isFinite(value)) {
      throw new Error(`condkelly ${label}: expected a finite number, got ${value}`);
    }
  }

  function assertPositive(value, label) {
    assertFinite(value, label);
    if (value <= 0) throw new Error(`condkelly ${label}: expected a positive number, got ${value}`);
  }

  function assertBetShape(bet, label) {
    assertFinite(bet.cut, `${label}.cut`);
    if (!Number.isInteger(bet.cut * 2)) {
      throw new Error(`condkelly ${label}.cut: expected a whole or half-point number, got ${bet.cut}`);
    }
    if (bet.direction !== "above" && bet.direction !== "below") {
      throw new Error(`condkelly ${label}.direction: expected "above" or "below", got ${bet.direction}`);
    }
    if (typeof bet.group !== "string" || bet.group === "") {
      throw new Error(`condkelly ${label}.group: expected a period name, got ${bet.group}`);
    }
  }

  function isWholeNumber(cut) {
    return Number.isInteger(cut);
  }

  // A teaser leg loses on a push, so the caller passes the half-point one
  // step against it (Bills -1 -> -1.5): a whole number here is a caller bug.
  function assertHalfPointLeg(leg, label) {
    assertFinite(leg.cut, `${label}.cut`);
    if (!Number.isInteger(leg.cut * 2) || isWholeNumber(leg.cut)) {
      throw new Error(`condkelly ${label}.cut: expected a half-point number (a teaser push loses), got ${leg.cut}`);
    }
    if (leg.direction !== "above" && leg.direction !== "below") {
      throw new Error(`condkelly ${label}.direction: expected "above" or "below", got ${leg.direction}`);
    }
  }

  function assertTicketShape(ticket, label) {
    assertBetShape(ticket, label);
    assertHalfPointLeg(ticket, label);
    assertPositive(ticket.stake, `${label}.stake`);
    assertPositive(ticket.toWin, `${label}.toWin`);
    if (!Array.isArray(ticket.others)) throw new Error(`condkelly ${label}.others: expected an array, got ${ticket.others}`);
    ticket.others.forEach((leg, i) => {
      if (typeof leg.factor !== "string" || leg.factor === "") {
        throw new Error(`condkelly ${label}.others[${i}].factor: expected a game name, got ${leg.factor}`);
      }
      assertHalfPointLeg(leg, `${label}.others[${i}]`);
    });
  }

  // Where a bet stops losing and starts winning. A half-point number is one
  // cut; a whole number k loses past k-0.5 on its wrong side, wins past
  // k+0.5 on its right side, and pushes between.
  function winLoseCuts(bet) {
    if (!isWholeNumber(bet.cut)) return { winCut: bet.cut, loseCut: bet.cut };
    if (bet.direction === "above") return { winCut: bet.cut + HALF_POINT, loseCut: bet.cut - HALF_POINT };
    return { winCut: bet.cut - HALF_POINT, loseCut: bet.cut + HALF_POINT };
  }

  // The half-point cuts a bet splits its axis at: what the caller must
  // supply a ladder probability for.
  function cutsNeeded(bet) {
    const { winCut, loseCut } = winLoseCuts(bet);
    return winCut === loseCut ? [winCut] : [Math.min(winCut, loseCut), Math.max(winCut, loseCut)];
  }

  // +1 win, 0 push, -1 loss for a result sitting at `value`.
  function resultAt(bet, value) {
    const { winCut, loseCut } = winLoseCuts(bet);
    if (bet.direction === "above") return value > winCut ? 1 : value < loseCut ? -1 : 0;
    return value < winCut ? 1 : value > loseCut ? -1 : 0;
  }

  function ladderMap(ladders, group) {
    const pairs = ladders && Array.isArray(ladders[group]) ? ladders[group] : [];
    return new Map(pairs.map(([cut, probAbove]) => [cut, probAbove]));
  }

  function ladderProbAbove(ladder, group, cut) {
    const probAbove = ladder.get(cut);
    if (typeof probAbove !== "number" || !(probAbove >= 0 && probAbove <= 1)) {
      throw new Error(`condkelly: expected a ladder probability for ${group} at ${cut}, got ${probAbove}`);
    }
    return probAbove;
  }

  // The candidate's win chance comes from Unabated's edge, not the ladder, so
  // it overrides the ladder at its own cut — a held bet on that cut then needs
  // no rung. On a whole number the edge's chance is conditional on no push:
  // the push mass comes from the two neighbouring rungs.
  function candidateOverrides(candidate, ladder) {
    if (!isWholeNumber(candidate.cut)) {
      return [[candidate.cut, candidate.direction === "above" ? candidate.prob : 1 - candidate.prob]];
    }
    const lower = candidate.cut - HALF_POINT;
    const upper = candidate.cut + HALF_POINT;
    const pushProb = Math.max(0, ladderProbAbove(ladder, candidate.group, lower) - ladderProbAbove(ladder, candidate.group, upper));
    const winProb = candidate.prob * (1 - pushProb);
    if (candidate.direction === "above") return [[lower, winProb + pushProb], [upper, winProb]];
    return [[lower, 1 - winProb], [upper, 1 - winProb - pushProb]];
  }

  // The outcome rows of one axis cut at every cut in `probAboveByCut` (cut
  // -> P(result > cut)): the gaps between the sorted cuts, each with its
  // chance and a value strictly inside it. {rows: [{prob, value}]} or
  // {reason} when the chances do not fall as the cut rises.
  function rowsBetweenCuts(probAboveByCut) {
    const cuts = Array.from(probAboveByCut.keys()).sort((a, b) => a - b);
    const edges = [-Infinity, ...cuts, Infinity];
    const probsAbove = [1, ...cuts.map((cut) => probAboveByCut.get(cut)), 0];
    const rows = [];
    let totalProb = 0;
    for (let i = 0; i < edges.length - 1; i += 1) {
      const slice = probsAbove[i] - probsAbove[i + 1];
      if (slice < -MAX_NEGATIVE_SLICE) return { reason: REASON_NOT_MONOTONE };
      const prob = Math.max(0, slice);
      totalProb += prob;
      // Any value strictly inside the gap: every cut is a half-point, so a
      // quarter past the lower edge (or before the upper) is inside.
      const value = edges[i] === -Infinity ? edges[i + 1] - HALF_POINT / 2 : edges[i] + HALF_POINT / 2;
      rows.push({ prob, value });
    }
    // Clamped slices leave the total a hair over 1.
    const live = rows.filter((row) => row.prob > PROB_EPSILON);
    for (const row of live) row.prob /= totalProb;
    return { rows: live };
  }

  // The held bets' P&L at a result `value`: +toWin, 0 on a push, -stake.
  function heldPnlAt(heldBets, value) {
    let pnl = 0;
    for (const bet of heldBets) {
      const result = resultAt(bet, value);
      pnl += result > 0 ? bet.toWin : result < 0 ? -bet.stake : 0;
    }
    return pnl;
  }

  // One group's outcome rows: the gaps between its sorted cuts — the
  // candidate's, the held bets' and the tickets' legs'. Each row carries its
  // chance, its value (the tickets' P&L is added per survival class), the
  // held bets' P&L there, and whether the candidate wins (+1), pushes (0) or
  // loses (-1) there — null when the candidate is in another group. Returns
  // { rows } or { reason }.
  function groupRows(group, heldBets, groupTickets, candidate, ladders) {
    const ladder = ladderMap(ladders, group);
    const overrides = new Map(candidate ? candidateOverrides(candidate, ladder) : []);
    const probAboveByCut = new Map(overrides);
    const cutsHeld = [...heldBets.flatMap(cutsNeeded), ...groupTickets.map((ticket) => ticket.cut)];
    for (const cut of cutsHeld) {
      if (!overrides.has(cut)) probAboveByCut.set(cut, ladderProbAbove(ladder, group, cut));
    }
    const built = rowsBetweenCuts(probAboveByCut);
    if (built.reason) return built;
    return {
      rows: built.rows.map((row) => ({
        prob: row.prob, value: row.value, heldPnl: heldPnlAt(heldBets, row.value),
        candidateResult: candidate ? resultAt(candidate, row.value) : null,
      })),
    };
  }

  // The other games under the tickets, each cut at every number one of their
  // legs needs: {games: [{factor, rows}]} or {reason}.
  function otherGames(tickets, factors) {
    const cutsByFactor = new Map();
    for (const ticket of tickets) {
      for (const leg of ticket.others) {
        if (!cutsByFactor.has(leg.factor)) cutsByFactor.set(leg.factor, new Set());
        cutsByFactor.get(leg.factor).add(leg.cut);
      }
    }
    const games = [];
    let outcomeCount = 1;
    for (const [factor, cuts] of cutsByFactor) {
      const ladder = ladderMap(factors, factor);
      const built = rowsBetweenCuts(new Map(Array.from(cuts, (cut) => [cut, ladderProbAbove(ladder, factor, cut)])));
      if (built.reason) return built;
      outcomeCount *= built.rows.length;
      if (outcomeCount > MAX_OTHER_GAME_OUTCOMES) return { reason: REASON_TOO_MANY_GAMES };
      games.push({ factor, rows: built.rows });
    }
    return { games };
  }

  // Which tickets the other games leave alive, and how likely each way is.
  // The games are independent, so a joint outcome's chance is the product of
  // its rows'; outcomes leaving the same tickets alive are pooled, since that
  // is all the score reads. A ticket with no other leg is always alive.
  //   {classes: [{prob, alive: [bool per ticket]}]} or {reason}
  function survivalClasses(tickets, factors) {
    const found = otherGames(tickets, factors);
    if (found.reason) return found;
    const games = found.games;
    const byAlive = new Map();
    const valueByFactor = new Map();
    const visit = (gameIndex, prob) => {
      if (gameIndex === games.length) {
        const alive = tickets.map((ticket) => ticket.others.every((leg) => resultAt(leg, valueByFactor.get(leg.factor)) > 0));
        const key = alive.map(Number).join("");
        if (byAlive.has(key)) byAlive.get(key).prob += prob;
        else byAlive.set(key, { prob, alive });
        return;
      }
      for (const row of games[gameIndex].rows) {
        valueByFactor.set(games[gameIndex].factor, row.value);
        visit(gameIndex + 1, prob * row.prob);
      }
    };
    visit(0, 1);
    return { classes: Array.from(byAlive.values()) };
  }

  // The tickets' P&L at a result `value` of this game in `group`: +toWin
  // where the leg wins and the other games left the ticket alive, else -stake.
  function ticketsPnlAt(tickets, alive, group, value) {
    let pnl = 0;
    tickets.forEach((ticket, index) => {
      if (ticket.group !== group) return;
      pnl += resultAt(ticket, value) > 0 && alive[index] ? ticket.toWin : -ticket.stake;
    });
    return pnl;
  }

  function rowPnl(row, stake, netOdds) {
    if (row.candidateResult === null || row.candidateResult === 0) return row.heldPnl;
    return row.heldPnl + (row.candidateResult > 0 ? stake * netOdds : -stake);
  }

  // Expected log growth of K at `stake`. One group: a plain sum. Several
  // groups (periods): Unabated says nothing about how two periods move
  // together, so take the worst case the fairs allow — sort each group's rows
  // by P&L (the candidate's group depends on the stake, so inside here) and
  // pair them by cumulative probability: bad with bad, good with good.
  function expectedLogGrowth(kellyBankroll, groups, stake, netOdds) {
    const sorted = groups.map((rows) => rows
      .map((row) => ({ prob: row.prob, pnl: rowPnl(row, stake, netOdds) }))
      .sort((a, b) => a.pnl - b.pnl));
    const index = sorted.map(() => 0);
    const remaining = sorted.map((rows) => rows[0].prob);
    let total = 0;
    for (;;) {
      let step = Infinity;
      let pnl = 0;
      for (let g = 0; g < sorted.length; g += 1) {
        step = Math.min(step, remaining[g]);
        pnl += sorted[g][index[g]].pnl;
      }
      const wealth = kellyBankroll + pnl;
      if (wealth <= 0) return -Infinity;
      total += step * Math.log(wealth / kellyBankroll);
      for (let g = 0; g < sorted.length; g += 1) {
        remaining[g] -= step;
        if (remaining[g] > PROB_EPSILON) continue;
        index[g] += 1;
        if (index[g] >= sorted[g].length) return total;
        remaining[g] = sorted[g][index[g]].prob;
      }
    }
  }

  // The largest stake that keeps K + P&L positive on every row: { upper } or
  // { reason } when the held bets alone can already lose K. The worst joint
  // row pairs every group's worst row, and only a row the candidate loses
  // gets worse with x.
  function stakeUpperBound(kellyBankroll, groups, candidateGroupIndex) {
    let worstElsewhere = 0;
    groups.forEach((rows, g) => {
      if (g !== candidateGroupIndex) worstElsewhere += Math.min(...rows.map((row) => row.heldPnl));
    });
    const candidateRows = groups[candidateGroupIndex];
    const room = (row) => kellyBankroll + worstElsewhere + row.heldPnl;
    if (candidateRows.some((row) => room(row) <= 0)) return { reason: REASON_HELD_RISK };
    const losing = candidateRows.filter((row) => row.candidateResult < 0);
    if (losing.length === 0) return { reason: REASON_NO_LOSING_ROW };
    return { upper: Math.min(...losing.map(room)) * UPPER_BOUND_SHRINK };
  }

  // The score is a single hill (a minimum of concave functions), so a bounded
  // golden-section search finds its top.
  function goldenSectionMax(score, upper) {
    let lo = 0;
    let hi = upper;
    let left = hi - GOLDEN_RATIO * (hi - lo);
    let right = lo + GOLDEN_RATIO * (hi - lo);
    let scoreLeft = score(left);
    let scoreRight = score(right);
    for (let i = 0; i < SEARCH_MAX_ITERATIONS && hi - lo > SEARCH_TOLERANCE_DOLLARS; i += 1) {
      if (scoreLeft > scoreRight) {
        hi = right;
        right = left;
        scoreRight = scoreLeft;
        left = hi - GOLDEN_RATIO * (hi - lo);
        scoreLeft = score(left);
      } else {
        lo = left;
        left = right;
        scoreLeft = scoreRight;
        right = lo + GOLDEN_RATIO * (hi - lo);
        scoreRight = score(right);
      }
    }
    return (lo + hi) / 2;
  }

  function solveStake({ kellyBankroll, candidate, held, tickets, factors, ladders }) {
    assertPositive(kellyBankroll, "kellyBankroll");
    assertBetShape(candidate, "candidate");
    assertPositive(candidate.netOdds, "candidate.netOdds");
    assertFinite(candidate.prob, "candidate.prob");
    if (!(candidate.prob > 0 && candidate.prob < 1)) {
      throw new Error(`condkelly candidate.prob: expected a probability inside (0, 1), got ${candidate.prob}`);
    }
    const heldBets = Array.isArray(held) ? held : [];
    heldBets.forEach((bet, i) => {
      assertBetShape(bet, `held[${i}]`);
      assertPositive(bet.stake, `held[${i}].stake`);
      assertPositive(bet.toWin, `held[${i}].toWin`);
    });
    const heldTickets = Array.isArray(tickets) ? tickets : [];
    heldTickets.forEach((ticket, i) => assertTicketShape(ticket, `tickets[${i}]`));

    const groupNames = [candidate.group];
    for (const bet of [...heldBets, ...heldTickets]) if (!groupNames.includes(bet.group)) groupNames.push(bet.group);
    const groups = [];
    for (const name of groupNames) {
      const built = groupRows(name, heldBets.filter((bet) => bet.group === name), heldTickets.filter((ticket) => ticket.group === name),
        name === candidate.group ? candidate : null, ladders);
      if (built.reason) return { stake: null, reason: built.reason };
      groups.push({ name, rows: built.rows });
    }
    const survival = survivalClasses(heldTickets, factors);
    if (survival.reason) return { stake: null, reason: survival.reason };
    // One outcome table per survival class: every row with the tickets' P&L in.
    const scenarios = survival.classes.map((survivalClass) => ({
      prob: survivalClass.prob,
      groups: groups.map(({ name, rows }) => rows.map((row) => ({
        prob: row.prob, candidateResult: row.candidateResult,
        heldPnl: row.heldPnl + ticketsPnlAt(heldTickets, survivalClass.alive, name, row.value),
      }))),
    }));

    let upper = Infinity;
    for (const scenario of scenarios) {
      const bound = stakeUpperBound(kellyBankroll, scenario.groups, 0);
      if (bound.reason) return { stake: null, reason: bound.reason };
      upper = Math.min(upper, bound.upper);
    }
    // Each class's growth is concave in the stake (a minimum of concave
    // functions, as above) and a sum of concave functions is concave, so the
    // search below still finds the top.
    const score = (stake) => {
      let total = 0;
      for (const scenario of scenarios) {
        const growth = expectedLogGrowth(kellyBankroll, scenario.groups, stake, candidate.netOdds);
        if (growth === -Infinity) return -Infinity;
        total += scenario.prob * growth;
      }
      return total;
    };
    const probe = Math.min(kellyBankroll * SLOPE_PROBE_FRACTION, upper / 2);
    if (score(probe) <= score(0)) return { stake: 0, reason: null };
    return { stake: goldenSectionMax(score, upper), reason: null };
  }

  const api = { REASON_NOT_MONOTONE, REASON_HELD_RISK, REASON_TOO_MANY_GAMES, cutsNeeded, resultAt, solveStake };

  if (typeof module !== "undefined" && module.exports) {
    module.exports = api;
  } else {
    root.UnabatedCondKelly = api;
  }
})(typeof globalThis !== "undefined" ? globalThis : this);
