// Unabated Ticket server — one Kalshi RFQ leg priced off Unabated (#282, epic
// #281): the Python parlay quoter (kalshi_parlay_quoter, #283) asks
// GET /rfq_price.json with an RFQ's leg tickers and multiplies the fairs.
// Pure pricing: no orders, no writes, no clock (now comes in), no fetch (the
// runner reads Kalshi's market, book and start and passes them in).
//
// Per leg (a Kalshi market ticker and the side the RFQ buys, yes or no):
//   1. the ticker: a full-game moneyline, spread or total series of NFL, CFB,
//      NBA, CBB or MLB (SERIES); its event suffix ("26SEP19DUQWSU") and its
//      strike. Kalshi names a strike X.5 as X+1 ("-52" is Over 51.5, "WSU36"
//      is WSU wins by over 35.5), as the taker reads it (taker.js
//      derivedMainNumberContract), and Kalshi's floor_strike must agree.
//   2. the game: the one board event of the series' league whose Kalshi rungs
//      carry that event suffix (feed.noteRungVenueIds). A game Unabated lists
//      no Kalshi rung for is not matched by name: no_event.
//   3. the YES team (spreads, moneylines): off Unabated's own rung for this
//      ticker when it has one, else the game's spread rungs' team codes
//      (taker.spreadTeamCodes), else the ticker's team codes
//      (taker.kalshiYesTeamSide). Two readings that disagree: no_team.
//   4. the fair: P(YES) off Unabated's fair ladder (ladder.js) at the cut;
//      "margin" = away minus home, so YES on the away team at X is
//      P(margin > X), on the home team 1 - P(margin > -X), a moneyline the
//      +/-0.5 cut. A cut with no usable rung is bracketed by the rungs around
//      it (probBounds): fair is the midpoint, `bracketed` true; a cut past
//      the ladder's last usable rung is declined (past_ladder), and so is a
//      flat-lined rung (flat_tail). A NO leg is 1 - P(YES); on a football
//      moneyline that NO includes the tie, which the ladder's +/-0.5 cuts
//      price exactly.
//   5. sigma: tail flex (tailflex.js) in probability units, dnorm(z) x c x
//      |z - z_main|, z the leg's fair and z_main Unabated's main fair on the
//      same side (0 at the main number and on moneylines; even money when
//      Unabated has no main line for the side, as tailflex.shrinkOfRow does).
//      c and pooled tau ride along for the quoter's margin (#284).
//   6. Kalshi: the market must be active at the ticker's strike; its live bid
//      and ask on YES, and the NO ask, as dollars.
//
// Decline reasons, in the order they are checked:
//   unsupported_series  not a full-game ML/spread/total series of the five leagues, or an unreadable ticker
//   no_board            the runner has no board yet
//   stale_feed          the league's snapshot is stale, failed or too old (the taker's own gate)
//   no_event            no board event of the league carries the ticker's Kalshi event
//   ambiguous_event     more than one does
//   no_team             the YES team cannot be read, or two readings disagree
//   flat_tail           the rung at the cut repeats a neighbour's fair (Unabated's flat deep tail)
//   past_ladder         no usable rung on one side of the cut
//   started             the earlier of Unabated's and Kalshi's start has passed
//   kalshi_unread       Kalshi's market or book could not be read
//   market_mismatch     Kalshi's market is not active, or its floor_strike is not the ticker's

"use strict";

const feed = require("../extension/feed.js");
const ladderLib = require("../extension/ladder.js");
const tailflex = require("../extension/tailflex.js");
const taker = require("./taker.js");

const BET_TYPE_MONEYLINE = 1;
const BET_TYPE_SPREAD = 2;
const BET_TYPE_TOTAL = 3;
const PERIOD_FULL_GAME = 1;
const SIDE_AWAY_OR_OVER = 0;
const MONEYLINE_CUT = 0.5;
const MARKET_NAMES = { [BET_TYPE_MONEYLINE]: "ml", [BET_TYPE_SPREAD]: "spread", [BET_TYPE_TOTAL]: "total" };
// Kalshi series prefix -> Unabated league id (feed.LEAGUES); Cal 2026-11-07:
// every league Unabated has fairs for, full game only.
const LEAGUE_OF_SERIES_PREFIX = { KXNFL: 1, KXNCAAF: 2, KXNBA: 3, KXNCAAMB: 4, KXMLB: 5 };
const BET_TYPE_OF_SERIES_KIND = { GAME: BET_TYPE_MONEYLINE, SPREAD: BET_TYPE_SPREAD, TOTAL: BET_TYPE_TOTAL };
const SERIES = Object.fromEntries(Object.entries(LEAGUE_OF_SERIES_PREFIX).flatMap(([prefix, leagueId]) =>
  Object.entries(BET_TYPE_OF_SERIES_KIND).map(([kind, betTypeId]) => [`${prefix}${kind}`, { leagueId, betTypeId }])));
const RFQ_LEAGUE_IDS = Object.values(LEAGUE_OF_SERIES_PREFIX);
const SIDES = new Set(["yes", "no"]);
// A leg is "<market ticker>:<side>"; an RFQ of more legs than this is not a parlay request, it is a bug.
const MAX_LEGS = 50;
const TICKER_RE = /^[A-Z0-9.]+-[A-Z0-9]+-[A-Z0-9.]+$/;
const SPREAD_STRIKE_RE = /^([A-Z]+)(\d+)$/;
const MONEYLINE_STRIKE_RE = /^[A-Z]+$/;
const TOTAL_STRIKE_RE = /^\d+$/;
const MARKET_STATUS_ACTIVE = "active";
const PROB_DECIMALS = 1e6;

const REASONS = ["unsupported_series", "no_board", "stale_feed", "no_event", "ambiguous_event", "no_team", "flat_tail", "past_ladder", "started", "kalshi_unread", "market_mismatch"];

function badRequest(message) {
  const error = new Error(message);
  error.status = 400;
  return error;
}

function roundProb(value) {
  return value == null ? null : Math.round(value * PROB_DECIMALS) / PROB_DECIMALS;
}

// "T1:yes,T2:no" -> [{marketTicker, side}], request order kept. A malformed
// list is the caller's bug: a 400, never a silently shorter answer.
function parseLegsParam(text) {
  if (typeof text !== "string" || !text.trim()) throw badRequest("legs is required: TICKER:yes|no, comma-separated");
  const legs = text.split(",").map((part) => {
    const cut = part.lastIndexOf(":");
    const marketTicker = cut > 0 ? part.slice(0, cut).trim() : "";
    const side = cut > 0 ? part.slice(cut + 1).trim().toLowerCase() : "";
    if (!TICKER_RE.test(marketTicker) || !SIDES.has(side)) throw badRequest(`leg ${JSON.stringify(part)} is not TICKER:yes|no`);
    return { marketTicker, side };
  });
  if (legs.length > MAX_LEGS) throw badRequest(`${legs.length} legs, at most ${MAX_LEGS}`);
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
    if (!TOTAL_STRIKE_RE.test(strike)) return { reason: "unsupported_series", detail: `total strike ${strike} is not a whole number` };
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

// Unabated's own rung for this spread ticker: {side, points} of the YES team,
// or null when no rung names it. A Y- rung sits on the YES team's side, an
// N- rung on the other team's.
function yesTeamFromOwnRung(event, marketTicker) {
  const contracts = event.venueIds.kalshiContracts || {};
  const readings = [];
  for (const prefix of ["Y", "N"]) {
    const target = contracts[`${prefix}-${marketTicker}`];
    if (!target || (target.sideIndex !== 0 && target.sideIndex !== 1)) continue;
    readings.push({ side: prefix === "Y" ? target.sideIndex : 1 - target.sideIndex, absPoints: Math.abs(target.points) });
  }
  return readings;
}

// Which Unabated side (0 away, 1 home) the YES team is: {side, source} or {reason, detail}.
function yesTeamSideOf(parsed, marketTicker, event, board) {
  const readings = [];
  if (parsed.betTypeId === BET_TYPE_SPREAD) {
    for (const reading of yesTeamFromOwnRung(event, marketTicker)) {
      if (reading.absPoints !== parsed.number) {
        return { reason: "no_team", detail: `Unabated's rung for ${marketTicker} sits at ${reading.absPoints}, the ticker says ${parsed.number}` };
      }
      readings.push({ side: reading.side, source: "rung" });
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
  const sides = new Set(readings.map((reading) => reading.side));
  if (sides.size > 1) return { reason: "no_team", detail: `${parsed.strikeCode} reads as both sides: ${readings.map((r) => `${r.source} ${r.side}`).join(", ")}` };
  return readings[0];
}

// Map a ticker onto the board: {parsed, event, yesSide} (yesSide null on a
// total) or {reason, detail, parsed?}. No Kalshi data needed.
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

// Where YES sits on the ladder: {axis, cut, yesAbove}. YES is P(above cut)
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
  const usable = ladder.rungs.filter(([rungCut]) => ladderLib.probAbove(ladder, rungCut).prob != null);
  const below = usable.some(([rungCut]) => rungCut < cut);
  const above = usable.some(([rungCut]) => rungCut > cut);
  if (!below || !above) {
    const span = usable.length ? `${usable[0][0]} to ${usable[usable.length - 1][0]}` : "no usable rung";
    return { reason: "past_ladder", detail: `cut ${cut} is outside the ${ladder.axis} ladder (${span})` };
  }
  const bounds = ladderLib.probBounds(ladder, cut);
  return { lo: bounds.lo, hi: bounds.hi, bracketed: true };
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

// Tail-flex sigma of the leg's fair, in probability units, with its c and tau.
function sigmaOf({ parsed, event, legLine, fair, measurement }) {
  const market = { leagueId: parsed.leagueId, periodTypeId: PERIOD_FULL_GAME, betTypeId: parsed.betTypeId };
  const c = tailflex.cOf(measurement, market);
  const tau = tailflex.tauOf(measurement, market);
  if (parsed.betTypeId === BET_TYPE_MONEYLINE || !(fair > 0 && fair < 1)) return { sigma: 0, c, tau, zMainSource: null };
  const sideKey = `${event.eventId}:${PERIOD_FULL_GAME}:${parsed.betTypeId}:${legLine.sideIndex}`;
  const main = measurement && measurement.mainFairs ? measurement.mainFairs.get(sideKey) : null;
  const z = tailflex.qnorm(fair);
  const zMain = main ? tailflex.qnorm(main.fairProb) : 0;
  return { sigma: tailflex.dnorm(z) * c * Math.abs(z - zMain), c, tau, zMainSource: main ? "unabated_main" : "even_money" };
}

function teamNameOf(board, event, sideIndex) {
  const teamId = sideIndex === 0 ? event.awayTeamId : event.homeTeamId;
  const teams = board.teams || {};
  return teams[teamId] || teams[String(teamId)] || null;
}

function lineFieldsOf(board, parsed, event, legLine) {
  const isTotal = parsed.betTypeId === BET_TYPE_TOTAL;
  return {
    league: feed.LEAGUES[parsed.leagueId].label,
    market: MARKET_NAMES[parsed.betTypeId],
    period: "FG",
    team: isTotal || !event ? null : teamNameOf(board, event, legLine.sideIndex),
    over_under: isTotal ? (legLine.sideIndex === SIDE_AWAY_OR_OVER ? "over" : "under") : null,
    number: legLine.points,
    unabated_event_id: event ? event.eventId : null,
    event_name: event ? event.eventName : null,
  };
}

// Kalshi's book as dollars: YES bid, YES ask, NO ask (1 - the YES bid).
function kalshiPricesOf(book) {
  const yesBid = book && book.yesBids.length ? book.yesBids[0].price : null;
  const noBid = book && book.noBids.length ? book.noBids[0].price : null;
  return {
    kalshi_yes_bid: yesBid,
    kalshi_yes_ask: noBid == null ? null : roundProb(1 - noBid),
    kalshi_no_ask: yesBid == null ? null : roundProb(1 - yesBid),
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

// One leg's answer (the endpoint's leg shape).
//   leg          {marketTicker, side}
//   resolved     resolveLeg's result
//   context      {board, linesOfEvent: eventId -> lines, measurement (tail flex, pooled tau),
//                 now, staleReasonOf: leagueId -> string|null, leagueLoadedAt: leagueId -> ms|null}
//   kalshi       {market, book, startMs, error} (error: the read that failed)
function priceLeg(leg, resolved, context, kalshi) {
  const out = emptyLeg(leg.marketTicker, leg.side);
  const decline = (reason, detail) => ({ ...out, ok: false, reason, detail });
  if (resolved.reason && !resolved.parsed) return decline(resolved.reason, resolved.detail);
  const { parsed } = resolved;
  if (resolved.reason === "no_board") return decline(resolved.reason, resolved.detail);
  const stale = context.staleReasonOf(parsed.leagueId);
  const loadedAt = context.leagueLoadedAt(parsed.leagueId);
  out.fair_age_ms = typeof loadedAt === "number" ? Math.max(0, context.now - loadedAt) : null;
  if (stale) return decline("stale_feed", stale);
  if (resolved.reason) return decline(resolved.reason, resolved.detail);
  const { event, yesSide } = resolved;
  const legLine = legLineOf(parsed, yesSide, leg.side);
  out.line = lineFieldsOf(context.board, parsed, event, legLine);

  const where = yesCutOf(parsed, yesSide);
  const ladder = ladderLib.buildLadder(context.linesOfEvent(event.eventId), { periodTypeId: PERIOD_FULL_GAME, axis: where.axis });
  const above = probAboveOf(ladder, where.cut);
  if (above.reason) return decline(above.reason, above.detail);
  // P(YES) bounds, then the leg's own side.
  const yesLo = where.yesAbove ? above.lo : 1 - above.hi;
  const yesHi = where.yesAbove ? above.hi : 1 - above.lo;
  const lo = leg.side === "yes" ? yesLo : 1 - yesHi;
  const hi = leg.side === "yes" ? yesHi : 1 - yesLo;
  const fair = (lo + hi) / 2;
  Object.assign(out, { fair: roundProb(fair), fair_lo: roundProb(lo), fair_hi: roundProb(hi), bracketed: above.bracketed });
  const flex = sigmaOf({ parsed, event, legLine, fair, measurement: context.measurement });
  Object.assign(out, { sigma: roundProb(flex.sigma), c: flex.c, tau: flex.tau });

  out.kickoff_ms = taker.earliestStartMs(event.eventStart, kalshi.startMs);
  if (out.kickoff_ms != null && context.now >= out.kickoff_ms) return decline("started", `the game started at ${new Date(out.kickoff_ms).toISOString()}`);
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

module.exports = {
  SERIES, RFQ_LEAGUE_IDS, REASONS, MAX_LEGS,
  parseLegsParam, parseTicker, resolveLeg, yesCutOf, probAboveOf, legLineOf, priceLeg, kalshiPricesOf,
};
