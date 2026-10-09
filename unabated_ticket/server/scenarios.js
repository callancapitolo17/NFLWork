// Live scenarios for the Bet Tracker's Live tab: for every game in progress
// that holds open bets, the results that change the money, each with its
// P&L and its chance at kickoff (Unabated's fair, frozen at the last snapshot
// before the start). Pure: no DOM, no fetch, no disk — the server runner
// (runner.js) calls it with its own board and bets.
//
// Inputs (buildScenarios)
//   records   the bets service's open and closed bet records (runner `held.records`,
//             team keys resolved); only open ones count
//   games     [{row, ladderOf, oddsAt}] — every board event the runner knows:
//             `row` one feed.describeLine row of the event (eventId, league,
//             awayTeam, homeTeam, eventStartMs, venueIds, rotations),
//             `ladderOf(period, axis)` the event's fair ladder as of `oddsAt`
//             (ladder.js shape) or null when the runner never saw the game
//             before it started, `oddsAt` ms or null
//   now       ms
// Output {games: [game], coveredBetIds} — `game` = {eventId, league,
//   awayTeam, homeTeam, startMs, oddsAt, groups, bets, legs}; each group is
//   one (period, axis) ladder: {key, period, axis, title, bands, ev, best,
//   worst} and each band {label, low, high, prob, pnl, legs}. `low` / `high`
//   are the half-point cuts the band sits between (null = open-ended); `prob`
//   is null when a rung is missing. `legs` per band are the parlay and teaser
//   legs on that market and how each fares there (won / lost / push); their
//   tickets' dollars are never in `pnl` (the other legs are other games).
//   `bets` lists every open straight bet on the game, each placed on a group
//   or with the `reason` it is not; `coveredBetIds` is every bet id shown on
//   a card, so the tracker can list the rest on their own.
// Side effects: none.

"use strict";

const betsLib = require("../extension/bets.js");
const condkelly = require("../extension/condkelly.js");
const ladderLib = require("../extension/ladder.js");

const AXIS_MARGIN = betsLib.AXIS_MARGIN;
const AXIS_TOTAL = betsLib.AXIS_TOTAL;
const FULL_GAME = "FG";
// Leagues whose full game can end level; elsewhere a margin of 0 is never a result.
const LEAGUES_WITH_TIES = new Set(["nfl", "soccer"]);
const PERIOD_WORDS = { FG: "", "1H": "1st half ", "2H": "2nd half ", "1Q": "1st quarter ", "2Q": "2nd quarter ", "3Q": "3rd quarter ", "4Q": "4th quarter " };
const RESULT_WON = "won";
const RESULT_LOST = "lost";
const RESULT_PUSH = "push";
// Two P&L values within a cent are the same band.
const SAME_PNL_TOLERANCE = 0.005;

function isOpen(record) {
  return Boolean(record) && record.status === "open" && !record.unmatchable && !record.mergedInto;
}

function resultWord(result) {
  return result > 0 ? RESULT_WON : result < 0 ? RESULT_LOST : RESULT_PUSH;
}

// ---- labels -------------------------------------------------------------------

// The whole-number results a band between two half-point cuts holds:
// [first, last], either end ±Infinity when the band is open.
function integerRange(low, high) {
  return [low === null ? -Infinity : low + 0.5, high === null ? Infinity : high - 0.5];
}

// "Bills by 3", "Bills by 2-5", "Bills by 4+", "Bills win" (first = 1, open top).
function teamBy(team, first, last) {
  if (last === Infinity) return first <= 1 ? `${team} win` : `${team} by ${first}+`;
  return first === last ? `${team} by ${first}` : `${team} by ${first}-${last}`;
}

// A margin band (away minus home) in words: "Chiefs by 2-3", "Bills win",
// "Chiefs by 3 or less, tie or Bills win".
function marginLabel(low, high, game) {
  const [first, last] = integerRange(low, high);
  const away = game.awayTeam || "Away";
  const home = game.homeTeam || "Home";
  if (first >= 1) return teamBy(away, first, last);
  if (last <= -1) return teamBy(home, -last, -first);
  if (first === -Infinity && last === Infinity) return "Any result";
  const upTo = (team, most) => (most === Infinity ? `${team} win` : most === 1 ? `${team} by 1` : `${team} by ${most} or less`);
  const parts = [];
  if (first <= -1) parts.push(upTo(home, -first));
  if (LEAGUES_WITH_TIES.has(game.league)) parts.push("tie");
  if (last >= 1) parts.push(upTo(away, last));
  if (parts.length === 0) return "Tie";
  const text = parts.length === 1 ? parts[0] : `${parts.slice(0, -1).join(", ")} or ${parts[parts.length - 1]}`;
  return text.charAt(0).toUpperCase() + text.slice(1);
}

// A total band in words: "48+ points", "47 or fewer", "44-47 points".
function totalLabel(low, high, game) {
  const unit = game.league === "mlb" ? "runs" : "points";
  const [rawFirst, last] = integerRange(low, high);
  const first = Math.max(0, rawFirst);
  if (last === Infinity) return `${first}+ ${unit}`;
  if (first === 0) return `${last} or fewer ${unit}`;
  return first === last ? `${first} ${unit}` : `${first}-${last} ${unit}`;
}

function groupTitle(period, axis, game) {
  const prefix = PERIOD_WORDS[period] ?? `${period} `;
  const text = axis === AXIS_TOTAL ? `${prefix}total ${game.league === "mlb" ? "runs" : "points"}` : `${prefix}result`;
  return text.charAt(0).toUpperCase() + text.slice(1);
}

function legLabel(record) {
  return betsLib.describeBet({ ...record, isParlayLeg: false, price: null });
}

// ---- bands --------------------------------------------------------------------

// Every half-point cut the positions split the axis at, ascending.
function cutsOf(positions) {
  const cuts = new Set();
  for (const position of positions) for (const cut of condkelly.cutsNeeded(position)) cuts.add(cut);
  return Array.from(cuts).sort((a, b) => a - b);
}

// A result inside the band (every result in it scores the same).
function sampleValue(low, high) {
  if (low === null && high === null) return 0;
  if (low === null) return high - 0.5;
  if (high === null) return low + 0.5;
  return (low + high) / 2;
}

// P(low < result < high) off the ladder, or null when a rung is missing.
function bandProb(ladder, low, high) {
  if (!ladder) return null;
  const above = (cut) => {
    if (cut === null) return null;
    const fair = ladderLib.probAbove(ladder, cut);
    return fair.reason ? undefined : fair.prob;
  };
  const lowAbove = low === null ? 1 : above(low);
  const highAbove = high === null ? 0 : above(high);
  if (lowAbove === undefined || highAbove === undefined) return null;
  return Math.max(0, lowAbove - highAbove);
}

function sameOutcome(a, b) {
  if (Math.abs(a.pnl - b.pnl) > SAME_PNL_TOLERANCE) return false;
  return a.legs.length === b.legs.length && a.legs.every((leg, index) => leg.id === b.legs[index].id && leg.result === b.legs[index].result);
}

// Neighbouring bands with the same P&L and the same leg results are one band.
function mergeBands(bands) {
  const merged = [];
  for (const band of bands) {
    const last = merged[merged.length - 1];
    if (last && sameOutcome(last, band)) {
      last.high = band.high;
      last.prob = last.prob === null || band.prob === null ? null : last.prob + band.prob;
    } else {
      merged.push({ ...band });
    }
  }
  return merged;
}

// One (period, axis) group: the bands its straight bets and legs cut it into,
// scored, merged and labelled (labels read the merged ends).
function groupOf(game, period, axis, straights, legs, ladder) {
  const cuts = cutsOf(straights.concat(legs));
  const edges = [null, ...cuts, null];
  const raw = [];
  for (let index = 0; index < edges.length - 1; index += 1) {
    const low = edges[index];
    const high = edges[index + 1];
    const value = sampleValue(low, high);
    let pnl = 0;
    for (const bet of straights) {
      const result = condkelly.resultAt(bet, value);
      pnl += result > 0 ? bet.toWin : result < 0 ? -bet.stake : 0;
    }
    const legResults = legs.map((leg) => ({ id: leg.id, label: leg.label, kind: leg.kind, result: resultWord(condkelly.resultAt(leg, value)) }));
    raw.push({ low, high, pnl, prob: bandProb(ladder, low, high), legs: legResults });
  }
  const label = axis === AXIS_TOTAL ? totalLabel : marginLabel;
  const bands = mergeBands(raw).map((band) => ({ ...band, label: label(band.low, band.high, game) }));
  const priced = bands.every((band) => band.prob !== null);
  return {
    key: `${period}:${axis}`, period, axis, title: groupTitle(period, axis, game), bands,
    ev: priced ? bands.reduce((sum, band) => sum + band.prob * band.pnl, 0) : null,
    best: Math.max(...bands.map((band) => band.pnl)),
    worst: Math.min(...bands.map((band) => band.pnl)),
  };
}

// FG first, margin before total, then the other periods in feed order.
function groupOrder(a, b) {
  const rank = (group) => (group.period === FULL_GAME ? 0 : 1);
  return rank(a) - rank(b) || a.period.localeCompare(b.period) || (a.axis === AXIS_MARGIN ? -1 : 1) - (b.axis === AXIS_MARGIN ? -1 : 1);
}

// ---- games --------------------------------------------------------------------

function betView(record, extra) {
  return {
    id: record.id, label: betsLib.describeBet(record), venue: betsLib.venueLabel(record.venue),
    stake: record.stake ?? null, toWin: record.toWin ?? null, ...extra,
  };
}

function legView(record, extra) {
  return {
    id: record.id, parlayId: record.parlayId, label: legLabel(record), venue: betsLib.venueLabel(record.venue),
    ticketStake: record.stake ?? null, ticketToWin: record.toWin ?? null, legCount: record.legCount ?? null,
    kind: /teaser/i.test(String((record.raw && (record.raw.headerDescription || record.raw.type)) || "")) ? "teaser" : "parlay",
    ...extra,
  };
}

// One game card from the matches on its row: straights placed on their
// groups, legs placed on theirs, the rest named with their reason.
function gameCard(game, matches) {
  const { row } = game;
  const info = { league: row.league, awayTeam: row.awayTeam, homeTeam: row.homeTeam };
  const byGroup = new Map();
  const groupFor = (period, axis) => {
    const key = `${period}:${axis}`;
    if (!byGroup.has(key)) byGroup.set(key, { period, axis, straights: [], legs: [] });
    return byGroup.get(key);
  };
  const bets = [];
  const legs = [];
  for (const { bet, position } of matches) {
    if (bet.isParlayLeg) {
      const legPosition = betsLib.teaserLegPositionOf(bet, row);
      const view = legView(bet, { group: legPosition.reason ? null : `${legPosition.period}:${legPosition.axis}`, reason: legPosition.reason || null });
      legs.push(view);
      if (!legPosition.reason) groupFor(legPosition.period, legPosition.axis).legs.push({ ...legPosition, id: bet.id, label: view.label, kind: view.kind });
      continue;
    }
    if (position.reason) {
      bets.push(betView(bet, { group: null, reason: position.reason }));
      continue;
    }
    bets.push(betView(bet, { group: `${position.period}:${position.axis}`, reason: null }));
    groupFor(position.period, position.axis).straights.push(position);
  }
  const groups = Array.from(byGroup.values())
    .map(({ period, axis, straights, legs: groupLegs }) => groupOf(info, period, axis, straights, groupLegs, game.ladderOf ? game.ladderOf(period, axis) : null))
    .sort(groupOrder);
  return {
    eventId: row.eventId, ...info, startMs: row.eventStartMs, oddsAt: game.oddsAt ?? null,
    groups, bets, legs,
  };
}

// Every game in progress with an open bet on it, earliest start first.
function buildScenarios({ records, games, now }) {
  const open = (records || []).filter(isOpen);
  const known = games || [];
  const live = known.filter((game) => typeof game.row.eventStartMs === "number" && game.row.eventStartMs <= now);
  if (open.length === 0 || live.length === 0) return { games: [], coveredBetIds: [] };
  // The whole board decides each bet's game (venue ids, two possible games);
  // only the started rows get cards.
  const boardRows = known.map((game) => game.row);
  const flags = betsLib.annotateRows(live.map((game) => game.row), open, { lines: boardRows });
  const cards = [];
  const covered = new Set();
  live.forEach((game, index) => {
    const matches = flags[index].matches;
    if (matches.length === 0) return;
    const card = gameCard(game, matches);
    for (const bet of card.bets) covered.add(bet.id);
    for (const leg of card.legs) covered.add(leg.id);
    cards.push(card);
  });
  cards.sort((a, b) => a.startMs - b.startMs || String(a.eventId).localeCompare(String(b.eventId)));
  return { games: cards, coveredBetIds: Array.from(covered).sort() };
}

module.exports = { buildScenarios, marginLabel, totalLabel, mergeBands, bandProb, RESULT_WON, RESULT_LOST, RESULT_PUSH };
