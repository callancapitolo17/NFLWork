// Unabated Ticket server — one Kalshi RFQ leg priced off Unabated (#282, epic
// #281): the Python parlay quoter (kalshi_parlay_quoter, #283) asks
// GET /rfq_price.json with an RFQ's leg tickers and multiplies the fairs.
// Pure pricing: no orders, no writes, no clock (now comes in), no fetch (the
// runner reads Kalshi's market, book and start and passes them in).
//
// Why it reads the way it does:
//   - Kalshi names a strike X.5 as X+1 ("-52" is Over 51.5, "WSU36" is WSU
//     wins by over 35.5), as the taker reads it (taker.js
//     derivedMainNumberContract); Kalshi's floor_strike must agree.
//   - The game is the one board event of the series' league whose Kalshi
//     rungs carry the ticker's event suffix (feed.noteRungVenueIds). A game
//     Unabated lists no Kalshi rung for is not matched by name: no_event.
//   - The YES team is read off Unabated's own rung for the ticker, the game's
//     spread rungs' team codes and the ticker's own codes; any two readings
//     that disagree decline (no_team).
//   - "margin" = away minus home (ladder.js), so YES on the away team at X is
//     P(margin > X), on the home team 1 - P(margin > -X); a moneyline is the
//     +/-0.5 cut. A cut with no usable rung is bracketed by the rungs around
//     it: fair is the midpoint of [fair_lo, fair_hi]. Past the ladder: declined.
//   - An NFL game can end level and Kalshi's NO on a team then wins, but
//     Unabated's two-way moneyline fair pushes the tie, so an NFL moneyline is
//     bracketed: YES in [p x (1 - NFL_TIE_PROB_MAX), p], NO the complement.
//     The other four leagues cannot tie.
//   - sigma is tail flex (tailflex.edgeUncertainty in probability units,
//     decimal 1): 0 at Unabated's main number and on moneylines. A bracketed
//     leg's own spread is fair_hi - fair_lo, not in sigma.
//   - `ok` alone says whether a leg may be quoted: a declined leg still
//     carries whatever was priced before the decline, for the dry-run log.

"use strict";

const feed = require("../extension/feed.js");
const betsLib = require("../extension/bets.js");
const ladderLib = require("../extension/ladder.js");
const tailflex = require("../extension/tailflex.js");
const taker = require("./taker.js");

const BET_TYPE_MONEYLINE = 1;
const BET_TYPE_SPREAD = 2;
const BET_TYPE_TOTAL = 3;
const PERIOD_FULL_GAME = 1;
const FULL_GAME = feed.PERIODS[PERIOD_FULL_GAME];
const SIDE_AWAY_OR_OVER = 0;
const MONEYLINE_CUT = 0.5;
const BET_TYPE_IDS = { moneyline: BET_TYPE_MONEYLINE, spread: BET_TYPE_SPREAD, total: BET_TYPE_TOTAL };
const MARKET_NAMES = { [BET_TYPE_MONEYLINE]: "ml", [BET_TYPE_SPREAD]: "spread", [BET_TYPE_TOTAL]: "total" };
// Cal 2026-11-07: every league Unabated has fairs for, full game only.
const RFQ_LEAGUE_IDS = [1, 2, 3, 4, 5];
const NFL_LEAGUE_ID = 1;
// The most an NFL game's tie can weigh: about 2 of 272 games a season end
// level under the 10-minute overtime, so 1% is a generous ceiling.
const NFL_TIE_PROB_MAX = 0.01;
// bets.js GAME_SERIES (the Tracker's series table), full game, the RFQ leagues.
const SERIES = Object.fromEntries(Object.entries(betsLib.GAME_SERIES).flatMap(([series, spec]) => {
  const leagueId = RFQ_LEAGUE_IDS.find((id) => feed.LEAGUES[id].path === spec.league);
  if (spec.period !== FULL_GAME || leagueId == null || spec.fixedPoints != null) return [];
  return [[series, { leagueId, betTypeId: BET_TYPE_IDS[spec.betType] }]];
}));
const SIDES = new Set(["yes", "no"]);
// More legs than this is not a parlay request, it is a bug.
const MAX_LEGS = 50;
// GET /rfq_fairs.json prices the quoter's whole leg table in chunks of this:
// 200 tickers stay under Node's 16 KB request-line limit.
const MAX_FAIR_LEGS = 200;
const TICKER_RE = /^[A-Z0-9.]+-[A-Z0-9]+-[A-Z0-9.]+$/;
const SPREAD_STRIKE_RE = /^([A-Z]+)(\d+)$/;
const MONEYLINE_STRIKE_RE = /^[A-Z]+$/;
const TOTAL_STRIKE_RE = /^[1-9]\d*$/;
const MARKET_STATUS_ACTIVE = "active";
const PROB_DECIMALS = 1e6;

function roundProb(value) {
  return value == null ? null : Math.round(value * PROB_DECIMALS) / PROB_DECIMALS;
}

// "T1:yes,T2:no" -> [{marketTicker, side}], request order kept, or a string
// saying why the list is malformed (the runner answers it with a 400).
function parseLegsParam(text, maxLegs = MAX_LEGS) {
  if (typeof text !== "string" || !text.trim()) return "legs is required: TICKER:yes|no, comma-separated";
  const legs = [];
  for (const part of text.split(",")) {
    const cut = part.lastIndexOf(":");
    const marketTicker = cut > 0 ? part.slice(0, cut).trim() : "";
    const side = cut > 0 ? part.slice(cut + 1).trim().toLowerCase() : "";
    if (!TICKER_RE.test(marketTicker) || !SIDES.has(side)) return `leg ${JSON.stringify(part)} is not TICKER:yes|no`;
    legs.push({ marketTicker, side });
  }
  if (legs.length > maxLegs) return `${legs.length} legs, at most ${maxLegs}`;
  return legs;
}

// The ticker's own facts: {series, leagueId, betTypeId, eventTicker,
// eventSuffix, strikeCode, number} (number: the strike as Kalshi's
// floor_strike, X.5; null on a moneyline), or {reason, detail}.
function parseTicker(marketTicker) {
  const parts = String(marketTicker || "").split("-");
  if (parts.length !== 3) return { reason: "unsupported_series", detail: `${marketTicker} is not series-event-market` };
  const [series, eventSuffix, strike] = parts;
  const known = SERIES[series];
  if (!known) return { reason: "unsupported_series", detail: `series ${series} is not a full-game moneyline, spread or total of NFL, CFB, NBA, CBB or MLB` };
  const base = { series, ...known, eventTicker: `${series}-${eventSuffix}`, eventSuffix };
  if (known.betTypeId === BET_TYPE_TOTAL) {
    if (!TOTAL_STRIKE_RE.test(strike)) return { reason: "unsupported_series", detail: `total strike ${strike} is not a positive whole number` };
    return { ...base, strikeCode: null, number: Number(strike) - 0.5 };
  }
  if (known.betTypeId === BET_TYPE_SPREAD) {
    const match = SPREAD_STRIKE_RE.exec(strike);
    if (!match || Number(match[2]) < 1) return { reason: "unsupported_series", detail: `spread strike ${strike} is not TEAM + number` };
    return { ...base, strikeCode: match[1], number: Number(match[2]) - 0.5 };
  }
  if (!MONEYLINE_STRIKE_RE.test(strike)) return { reason: "unsupported_series", detail: `moneyline strike ${strike} is not a team code` };
  return { ...base, strikeCode: strike, number: null };
}

// The board events of the ticker's league whose Kalshi rungs carry its event suffix.
function eventsOfSuffix(board, leagueId, eventSuffix) {
  return Object.values(board.events || {}).filter((event) => event && event.leagueId === leagueId
    && event.venueIds && Array.isArray(event.venueIds.kalshiEventSuffixes) && event.venueIds.kalshiEventSuffixes.includes(eventSuffix));
}

// Which Unabated side (0 away, 1 home) the YES team is: {side} or {reason, detail}.
function yesTeamSideOf(parsed, marketTicker, event, board) {
  const readings = [];
  if (parsed.betTypeId === BET_TYPE_SPREAD) {
    // A Y- rung sits on the YES team's side, an N- rung on the other team's.
    const contracts = event.venueIds.kalshiContracts || {};
    for (const prefix of ["Y", "N"]) {
      const target = contracts[`${prefix}-${marketTicker}`];
      if (!target || (target.sideIndex !== 0 && target.sideIndex !== 1)) continue;
      if (Math.abs(target.points) !== parsed.number) {
        return { reason: "no_team", detail: `Unabated's rung for ${marketTicker} sits at ${Math.abs(target.points)}, the ticker says ${parsed.number}` };
      }
      readings.push({ side: prefix === "Y" ? target.sideIndex : 1 - target.sideIndex, source: "rung" });
    }
  }
  const codes = taker.spreadTeamCodes({ venueIds: event.venueIds }, board);
  if (!codes.reason) {
    const index = codes.codes.indexOf(parsed.strikeCode);
    if (index >= 0 && codes.codes[1 - index] !== parsed.strikeCode) readings.push({ side: index, source: "team_codes" });
  }
  const fromTicker = taker.kalshiYesTeamSide(marketTicker);
  if (!fromTicker.reason) readings.push({ side: fromTicker.side, source: "ticker" });
  if (readings.length === 0) return { reason: "no_team", detail: `no reading of ${parsed.strikeCode}'s side: ${fromTicker.reason}` };
  if (new Set(readings.map((reading) => reading.side)).size > 1) {
    return { reason: "no_team", detail: `${parsed.strikeCode} reads as both sides: ${readings.map((r) => `${r.source} ${r.side}`).join(", ")}` };
  }
  return { side: readings[0].side };
}

// Map a ticker onto the board: {parsed, event, yesSide} (yesSide null on a
// total) or {reason, detail, parsed?, event?}. No Kalshi data needed.
function resolveLeg(marketTicker, board) {
  const parsed = parseTicker(marketTicker);
  if (parsed.reason) return parsed;
  if (!board) return { reason: "no_board", detail: "the board has not loaded yet", parsed };
  const events = eventsOfSuffix(board, parsed.leagueId, parsed.eventSuffix);
  if (events.length === 0) return { reason: "no_event", detail: `no ${feed.LEAGUES[parsed.leagueId].label} game on the board carries Kalshi event ${parsed.eventSuffix}`, parsed };
  if (events.length > 1) return { reason: "ambiguous_event", detail: `games ${events.map((event) => event.eventId).join(", ")} all carry Kalshi event ${parsed.eventSuffix}`, parsed };
  const event = events[0];
  if (parsed.betTypeId === BET_TYPE_TOTAL) return { parsed, event, yesSide: null };
  const team = yesTeamSideOf(parsed, marketTicker, event, board);
  if (team.reason) return { ...team, parsed, event };
  return { parsed, event, yesSide: team.side };
}

// Where YES sits on the ladder: {axis, cut, yesAbove}; YES is P(above cut)
// when yesAbove, else 1 - P(above cut).
function yesCutOf(parsed, yesSide) {
  if (parsed.betTypeId === BET_TYPE_TOTAL) return { axis: ladderLib.AXIS_TOTAL, cut: parsed.number, yesAbove: true };
  const number = parsed.betTypeId === BET_TYPE_MONEYLINE ? MONEYLINE_CUT : parsed.number;
  return yesSide === SIDE_AWAY_OR_OVER
    ? { axis: ladderLib.AXIS_MARGIN, cut: number, yesAbove: true }
    : { axis: ladderLib.AXIS_MARGIN, cut: -number, yesAbove: false };
}

// P(above cut) off the ladder: {lo, hi, bracketed} or {reason, detail}.
function probAboveOf(ladder, cut) {
  const exact = ladderLib.probAbove(ladder, cut);
  if (exact.prob != null) return { lo: exact.prob, hi: exact.prob, bracketed: false };
  if (exact.reason === ladderLib.REASON_FLAT) return { reason: "flat_tail", detail: `the rung at ${cut} repeats a neighbour's fair` };
  const usable = ladderLib.usableRungs(ladder);
  if (!usable.some(([rungCut]) => rungCut < cut) || !usable.some(([rungCut]) => rungCut > cut)) {
    const span = usable.length ? `${usable[0][0]} to ${usable[usable.length - 1][0]}` : "no usable rung";
    return { reason: "past_ladder", detail: `cut ${cut} is outside the ${ladder.axis} ladder (${span})` };
  }
  return { ...ladderLib.probBounds(ladder, cut), bracketed: true };
}

// The Unabated side and number the leg's own side is: {sideIndex, points}
// (points null on a moneyline).
function legLineOf(parsed, yesSide, side) {
  const yes = side === "yes";
  if (parsed.betTypeId === BET_TYPE_TOTAL) return { sideIndex: yes ? 0 : 1, points: parsed.number };
  const sideIndex = yes ? yesSide : 1 - yesSide;
  if (parsed.betTypeId === BET_TYPE_MONEYLINE) return { sideIndex, points: null };
  return { sideIndex, points: yes ? -parsed.number : parsed.number };
}

// Tail-flex sigma of the leg's fair in probability units, with its c and tau.
// z_main is Unabated's main fair on the leg's side, else even money (tailflex.shrinkOfRow's rule).
function sigmaOf({ parsed, event, legLine, fair, measurement }) {
  const market = { leagueId: parsed.leagueId, periodTypeId: PERIOD_FULL_GAME, betTypeId: parsed.betTypeId };
  const c = tailflex.cOf(measurement, market);
  const tau = tailflex.tauOf(measurement, market);
  if (parsed.betTypeId === BET_TYPE_MONEYLINE || !(fair > 0 && fair < 1)) return { sigma: 0, c, tau };
  const side = { eventId: event.eventId, periodTypeId: PERIOD_FULL_GAME, betTypeId: parsed.betTypeId, sideIndex: legLine.sideIndex };
  const main = measurement && measurement.mainFairs ? measurement.mainFairs.get(tailflex.sideKeyOf(side)) : null;
  const zMain = main ? tailflex.qnorm(main.fairProb) : 0;
  return { sigma: tailflex.edgeUncertainty({ decimal: 1, fairProb: fair, zMain, c }), c, tau };
}

function teamNameOf(board, event, sideIndex) {
  const teamId = sideIndex === 0 ? event.awayTeamId : event.homeTeamId;
  const teams = board.teams || {};
  return teams[teamId] || teams[String(teamId)] || null;
}

// The line the leg resolved to; before the game is found only its league, market and period.
function lineFieldsOf(board, parsed, event, legLine) {
  const isTotal = parsed.betTypeId === BET_TYPE_TOTAL;
  const resolved = Boolean(event && legLine);
  return {
    league: feed.LEAGUES[parsed.leagueId].label,
    market: MARKET_NAMES[parsed.betTypeId],
    period: FULL_GAME,
    team: !resolved || isTotal ? null : teamNameOf(board, event, legLine.sideIndex),
    over_under: resolved && isTotal ? (legLine.sideIndex === SIDE_AWAY_OR_OVER ? "over" : "under") : null,
    number: resolved ? legLine.points : null,
    unabated_event_id: event ? event.eventId : null,
    event_name: event ? event.eventName : null,
  };
}

// Kalshi's book as dollars: YES bid, YES ask, NO ask.
function kalshiPricesOf(book) {
  const yesAsk = taker.quoteOf(book, "yes");
  const noAsk = taker.quoteOf(book, "no");
  return {
    kalshi_yes_bid: book.yesBids.length ? book.yesBids[0].price : null,
    kalshi_yes_ask: yesAsk ? yesAsk.ask : null,
    kalshi_no_ask: noAsk ? noAsk.ask : null,
  };
}

function emptyLeg(marketTicker, side) {
  return {
    market_ticker: marketTicker, side, ok: false, reason: null, detail: null,
    fair: null, fair_lo: null, fair_hi: null, bracketed: null, sigma: null, c: null, tau: null,
    line: null, kickoff_ms: null, fair_age_ms: null,
    kalshi_status: null, kalshi_yes_bid: null, kalshi_yes_ask: null, kalshi_no_ask: null,
  };
}

// Whether a leg needs Kalshi read: one resolved to a game, in a fresh league.
function needsKalshi(resolved, context) {
  return !resolved.reason && !context.staleReasonOf(resolved.parsed.leagueId);
}

// One leg's answer (the endpoint's leg shape).
//   leg       {marketTicker, side}
//   resolved  resolveLeg's result
//   context   {board, ladderOf: eventId -> (period, axis) -> ladder, measurement (tail flex, pooled tau),
//              now, staleReasonOf: leagueId -> string|null, fairAsOf: leagueId -> ms|null}
//   kalshi    {market, book, startMs} or {error}; ignored unless needsKalshi.
//             {unchecked: true} (/rfq_fairs.json) prices off the board alone:
//             ok then says only that Unabated has a fair for a game it has
//             not seen start; kalshi_status stays null and the quoter must
//             confirm through /rfq_price.json before it relies on the market.
function priceLeg(leg, resolved, context, kalshi) {
  const out = emptyLeg(leg.marketTicker, leg.side);
  const decline = (reason, detail) => ({ ...out, ok: false, reason, detail });
  const { parsed } = resolved;
  if (!parsed) return decline(resolved.reason, resolved.detail);
  out.line = lineFieldsOf(context.board, parsed, resolved.event, null);
  if (resolved.reason === "no_board") return decline(resolved.reason, resolved.detail);
  const fairAsOf = context.fairAsOf(parsed.leagueId);
  out.fair_age_ms = typeof fairAsOf === "number" ? Math.max(0, context.now - fairAsOf) : null;
  const stale = context.staleReasonOf(parsed.leagueId);
  if (stale) return decline("stale_feed", stale);
  if (resolved.reason) return decline(resolved.reason, resolved.detail);
  const { event, yesSide } = resolved;
  const legLine = legLineOf(parsed, yesSide, leg.side);
  out.line = lineFieldsOf(context.board, parsed, event, legLine);

  const where = yesCutOf(parsed, yesSide);
  const ladder = context.ladderOf(event.eventId)(FULL_GAME, where.axis);
  const above = ladder ? probAboveOf(ladder, where.cut) : { reason: "past_ladder", detail: "no ladder for the game" };
  if (above.reason) return decline(above.reason, above.detail);
  let yesLo = where.yesAbove ? above.lo : 1 - above.hi;
  const yesHi = where.yesAbove ? above.hi : 1 - above.lo;
  let bracketed = above.bracketed;
  if (parsed.betTypeId === BET_TYPE_MONEYLINE && parsed.leagueId === NFL_LEAGUE_ID) {
    yesLo *= 1 - NFL_TIE_PROB_MAX;
    bracketed = true;
  }
  const lo = leg.side === "yes" ? yesLo : 1 - yesHi;
  const hi = leg.side === "yes" ? yesHi : 1 - yesLo;
  const fair = (lo + hi) / 2;
  Object.assign(out, { fair: roundProb(fair), fair_lo: roundProb(lo), fair_hi: roundProb(hi), bracketed });
  const flex = sigmaOf({ parsed, event, legLine, fair, measurement: context.measurement });
  Object.assign(out, { sigma: roundProb(flex.sigma), c: flex.c, tau: flex.tau });

  out.kickoff_ms = taker.earliestStartMs(event.eventStart, kalshi.startMs);
  if (out.kickoff_ms == null) return decline("no_start", "neither Unabated nor Kalshi lists the game's start");
  if (context.now >= out.kickoff_ms) return decline("started", `the game started at ${new Date(out.kickoff_ms).toISOString()}`);
  if (kalshi.unchecked) return { ...out, ok: true };
  if (kalshi.error) return decline("kalshi_unread", kalshi.error);
  const market = kalshi.market;
  out.kalshi_status = market.status ?? null;
  Object.assign(out, kalshiPricesOf(kalshi.book));
  if (market.status !== MARKET_STATUS_ACTIVE) return decline("market_mismatch", `${leg.marketTicker} is ${market.status}`);
  if (market.ticker && market.ticker !== leg.marketTicker) return decline("market_mismatch", `Kalshi answered ${market.ticker} for ${leg.marketTicker}`);
  if (parsed.number != null && Number(market.floor_strike) !== parsed.number) {
    return decline("market_mismatch", `floor_strike ${market.floor_strike}, the ticker says ${parsed.number}`);
  }
  return { ...out, ok: true };
}

module.exports = { RFQ_LEAGUE_IDS, MAX_LEGS, MAX_FAIR_LEGS, NFL_TIE_PROB_MAX, parseLegsParam, parseTicker, resolveLeg, needsKalshi, priceLeg };
