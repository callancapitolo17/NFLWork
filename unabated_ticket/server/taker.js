// Unabated Ticket server — the automated taker's decision, per Edges row
// (plan: automated_taker_plan.md). Every decision ends in a reason; a "bet"
// decision is sent by the runner to the bets service's order desk
// (bets_service/taker_orders.py), which holds the bot's Kalshi key, the stop
// file and the $500 ceiling. Live from merge (Cal, 2026-10-14: no shadow run).
//
// Cal's rules (2026-10-01 and 2026-10-14), each a named constant below:
//   - Kalshi rows Unabated shows an edge on, every league and market the
//     settings scan (alts included, moneylines through a derived ticker);
//   - fair = Unabated's own (`ge` on the row), no cross-check;
//   - Kalshi's live all-in price (ask + the 7% taker fee, the way Unabated
//     shows Kalshi) at most 0.1c worse than Unabated's shown price, and at
//     most 1c better: a price that ran away from Unabated's since its
//     snapshot is news Unabated's fair has not caught up with, not an edge;
//   - live edge at least 3% and at most 20%;
//   - no live book beats the price at that number (pricecheck.js, strictly:
//     one better book is enough to skip; outliers and Cal's Skips trade);
//   - lines changed within 7 days, right up to the start;
//   - $500 at 0.25 Kelly (Cal, 2026-10-14; was $1,000), conditional on every open bet (betsview.js
//     stakeAdvice), contracts cut to the depth at the ask.
//
// Reasons, in the order they are checked:
//   stale_feed    the row's league snapshot is stale or its last load failed
//   no_ticker     no Kalshi contract for the row (or an ambiguous one)
//   no_book       the orderbook read failed or has no ask on the row's side
//   price_moved   live all-in more than PRICE_TOLERANCE worse than Unabated's,
//                 or more than BETTER_PRICE_TOLERANCE better
//   edge_short    live edge under MIN_LIVE_EDGE
//   not_best      another live book beats the live price at the number
//   cap           live edge over MAX_LIVE_EDGE
//   size_declined conditional Kelly refused the held bets (fail closed)
//   size_zero     the stake buys no contract at the ask, or no depth
//   market_mismatch  at send time, Kalshi's market for the ticker is not
//                 active, or its strike or its YES team is not the row's
//                 (marketMismatch; the runner reads the market before every order)
//   bet           buy `contracts` at `liveAsk` (the runner sends it)
//
// Pure: no fetch, no clock (now comes in), no DOM. node-tested.

"use strict";

const crypto = require("node:crypto");
const feed = require("../extension/feed.js");
const kelly = require("../extension/kelly.js");
const betsView = require("../extension/betsview.js");
const pricecheck = require("../extension/pricecheck.js");

const KALSHI_BOOK_ID = 105;
const BOT_BANKROLL = 500;
const BOT_KELLY = 0.25;
const MIN_LIVE_EDGE = 0.03;
const MAX_LIVE_EDGE = 0.20;
// Unabated's Kalshi price sometimes lags the book by a tick; 0.1c worse still trades (Cal, 2026-10-14).
const PRICE_TOLERANCE = 0.001;
// A live price more than this much BETTER than Unabated's shown price is a
// market that moved since Unabated's snapshot (news), so its fair is stale:
// the edge the bot would compute there is its own, not Unabated's.
const BETTER_PRICE_TOLERANCE = 0.01;
// Kalshi's market status that takes orders (GET /markets/<ticker>).
const MARKET_STATUS_ACTIVE = "active";
const MAX_LINE_AGE_MS = 7 * 24 * 3600 * 1000;
// Kalshi's taker fee per contract is 0.07 * P * (1 - P); Unabated adds it unrounded.
const KALSHI_TAKER_FEE_RATE = 0.07;
// Unabated's sourceFormat for a price sent as a probability (Kalshi's all-in).
const SOURCE_FORMAT_PROBABILITY = 4;
const BET_TYPE_MONEYLINE = 1;
const BET_TYPE_SPREAD = 2;
const BET_TYPE_TOTAL = 3;
const PERIOD_FULL_GAME = 1;
const ALL_BET_TYPES = [1, 2, 3];
// "WSU36" on a spread contract: Kalshi's team code, then the strike.
const SPREAD_STRIKE_RE = /^([A-Z]+)(\d+)$/;
const KALSHI_CONTRACT_RE = /^([YN])-([A-Z0-9]+)-([A-Z0-9]+)-([A-Z0-9.]+)$/;
const SPREAD_SERIES_RE = /^(KX[A-Z0-9]+?)SPREAD$/;
const PRICE_DECIMALS = 1e4;

const REASONS = ["stale_feed", "no_ticker", "no_book", "price_moved", "edge_short", "not_best", "cap", "size_declined", "size_zero", "market_mismatch", "bet"];

function roundPrice(value) {
  return Math.round(value * PRICE_DECIMALS) / PRICE_DECIMALS;
}

// The price a taker pays per contract, fee included: 22c -> 0.23201.
function kalshiAllInPrice(ask) {
  return ask + KALSHI_TAKER_FEE_RATE * ask * (1 - ask);
}

// A probability as an American price, unrounded (0.25 -> +300, 0.6 -> -150).
function probToAmerican(prob) {
  if (!(prob > 0 && prob < 1)) throw new Error(`probToAmerican: expected a probability in (0, 1), got ${prob}`);
  return prob < 0.5 ? (100 * (1 - prob)) / prob : (-100 * prob) / (1 - prob);
}

// Unabated's shown Kalshi price as a probability: its sourcePrice when the
// feed sent one (Kalshi's all-in, unrounded), else the American price.
function unabatedAllInOf(row) {
  const sent = row.sourceFormat === SOURCE_FORMAT_PROBABILITY && typeof row.sourcePrice === "number" && row.sourcePrice > 0 && row.sourcePrice < 1;
  return sent ? row.sourcePrice : kelly.americanToProb(row.price);
}

// The rows the taker looks at: Kalshi lines at MIN_LIVE_EDGE and up by
// Unabated's own edge, in the settings' leagues and periods, alts included,
// every market, lines changed within MAX_LINE_AGE_MS. No liquidity gate: the
// live depth sizes the order instead.
//   edgeSettings  the settings the runner holds (edgerows.js shape)
function candidateRows(feedState, edgeSettings, now) {
  if (!feedState) return [];
  return feed.selectEdges(feedState, {
    leagueIds: new Set(edgeSettings.leagues),
    periods: new Set(edgeSettings.periods),
    betTypes: new Set(ALL_BET_TYPES),
    bookIds: new Set([KALSHI_BOOK_ID]),
    now,
    maxLineAgeMs: MAX_LINE_AGE_MS,
    includeAlts: true,
    minLiquidityToWin: null,
    minEdge: MIN_LIVE_EDGE,
  });
}

// The Kalshi contract a spread or total row buys, off the rung ids
// (feed.noteRungVenueIds): the one sourceKey whose listed line IS this row
// and whose number and side agree. "Y-<ticker>" buys YES, "N-<ticker>" NO.
//   {ticker, side: "yes" | "no", sourceKey, derived: false} or {reason, detail}
function contractFromRungs(row) {
  const contracts = (row.venueIds && row.venueIds.kalshiContracts) || {};
  const found = Object.entries(contracts).filter(([, target]) => target.lineKey === row.key);
  if (found.length === 0) return { reason: "no_ticker", detail: "no Kalshi rung names this line" };
  if (found.length > 1) return { reason: "no_ticker", detail: `${found.length} Kalshi rungs name this line: ${found.map(([key]) => key).join(", ")}` };
  const [sourceKey, target] = found[0];
  if (target.points !== row.points || target.sideIndex !== row.sideIndex) {
    return { reason: "no_ticker", detail: `rung ${sourceKey} is ${target.points} side ${target.sideIndex}, row is ${row.points} side ${row.sideIndex}` };
  }
  const match = KALSHI_CONTRACT_RE.exec(sourceKey);
  if (!match) return { reason: "no_ticker", detail: `rung id ${sourceKey} is not a Kalshi contract` };
  return { ticker: `${match[2]}-${match[3]}-${match[4]}`, side: match[1] === "Y" ? "yes" : "no", sourceKey, derived: false };
}

// Kalshi's team code for each Unabated side, read off the game's full-game
// spread rungs: "Y-KXNCAAFSPREAD-26SEP19DUQWSU-WSU36" on side 1 says WSU is
// side 1 (YES = WSU covers); the "N-" twin on side 0 says the same. Every
// rung must agree. Returns {series, eventSuffix, codes: [side0, side1]} or a
// {reason, detail}.
function spreadTeamCodes(row, feedState) {
  const contracts = (row.venueIds && row.venueIds.kalshiContracts) || {};
  const codes = [null, null];
  let series = null;
  let eventSuffix = null;
  for (const [sourceKey, target] of Object.entries(contracts)) {
    const main = feedState.lines[target.mainKey];
    if (!main || main.betTypeId !== BET_TYPE_SPREAD || main.periodTypeId !== PERIOD_FULL_GAME) continue;
    const match = KALSHI_CONTRACT_RE.exec(sourceKey);
    const strike = match ? SPREAD_STRIKE_RE.exec(match[4]) : null;
    const spreadSeries = match ? SPREAD_SERIES_RE.exec(match[2]) : null;
    if (!strike || !spreadSeries) return { reason: "no_ticker", detail: `spread rung ${sourceKey} has no team code` };
    if (series !== null && (series !== spreadSeries[1] || eventSuffix !== match[3])) {
      return { reason: "no_ticker", detail: `spread rungs name two events: ${series}-${eventSuffix} and ${spreadSeries[1]}-${match[3]}` };
    }
    series = spreadSeries[1];
    eventSuffix = match[3];
    const teamSide = match[1] === "Y" ? target.sideIndex : 1 - target.sideIndex;
    if (codes[teamSide] !== null && codes[teamSide] !== strike[1]) {
      return { reason: "no_ticker", detail: `side ${teamSide} reads as both ${codes[teamSide]} and ${strike[1]}` };
    }
    codes[teamSide] = strike[1];
  }
  if (series === null) return { reason: "no_ticker", detail: "no full-game Kalshi spread rung to read the team codes from" };
  return { series, eventSuffix, codes };
}

// A full-game moneyline's contract, which the feed never carries: the
// game series beside the spread series ("KXNCAAFSPREAD" -> "KXNCAAFGAME"),
// the event suffix, and the row's team code; YES on that team. Only one
// side's code is needed (Kalshi's spread markets may list one team only).
// Inferred from Kalshi's naming, so every one is marked `derived` and the
// orderbook read is its existence check.
function derivedMoneylineContract(row, feedState) {
  if (row.periodTypeId !== PERIOD_FULL_GAME) return { reason: "no_ticker", detail: "only full-game moneylines have a derivable Kalshi ticker" };
  const read = spreadTeamCodes(row, feedState);
  if (read.reason) return read;
  const code = read.codes[row.sideIndex];
  if (!code) return { reason: "no_ticker", detail: `no spread rung names side ${row.sideIndex}'s team code` };
  if (read.codes[1 - row.sideIndex] === code) return { reason: "no_ticker", detail: `both sides read as ${code}` };
  return { ticker: `${read.series}GAME-${read.eventSuffix}-${code}`, side: "yes", sourceKey: null, derived: true };
}

function resolveKalshiContract(row, feedState) {
  if (row.betTypeId === BET_TYPE_MONEYLINE) return derivedMoneylineContract(row, feedState);
  return contractFromRungs(row);
}

// One side's resting bids as [{price, quantity}] best first, from Kalshi's
// orderbook body: orderbook_fp.*_dollars string dollars, or the legacy
// integer-cents orderbook.* (unabated_edge/venues/kalshi.py).
function parseLevels(dollarLevels, centLevels) {
  const levels = Array.isArray(dollarLevels) && dollarLevels.length
    ? dollarLevels.map(([price, quantity]) => ({ price: Number(price), quantity: Number(quantity) }))
    : Array.isArray(centLevels) ? centLevels.map(([cents, quantity]) => ({ price: cents / 100, quantity: Number(quantity) })) : [];
  return levels.filter((level) => level.price > 0 && level.price < 1 && level.quantity > 0).sort((a, b) => b.price - a.price);
}

function parseKalshiOrderbook(body) {
  if (!body || typeof body !== "object") throw new Error("Kalshi orderbook: expected an object");
  const dollars = body.orderbook_fp || {};
  const cents = body.orderbook || {};
  return { yesBids: parseLevels(dollars.yes_dollars, cents.yes), noBids: parseLevels(dollars.no_dollars, cents.no) };
}

// What buying `side` costs now: the ask is 1 - the best bid on the other
// side, and the depth is that bid's size. Null when nothing rests there.
//   {ask, depth}
function quoteOf(book, side) {
  const opposite = side === "yes" ? book.noBids : book.yesBids;
  if (!opposite.length) return null;
  return { ask: roundPrice(1 - opposite[0].price), depth: opposite[0].quantity };
}

// One row's decision.
//   row         a candidateRows row
//   contract    resolveKalshiContract's result
//   quote       quoteOf's result, or null (read failed or empty side)
//   quoteError  why the orderbook read failed, when it did
//   priceIndex  pricecheck.buildPriceIndex over the same feed
//   sizing      {matches (bets.annotateRows for the row), ladderOf, teasers}
// Returns {reason, detail, ...the numbers it got to}.
function decideRow({ row, contract, quote, quoteError, priceIndex, sizing }) {
  const unabatedAllIn = unabatedAllInOf(row);
  // Unabated's edge is on its price, which it shows twice: the American
  // (-143) and Kalshi's all-in (0.587157 = -142.2). Which one `ge` is measured
  // on is not documented, so the fair is read off the cheaper of the two:
  // the lower fair, never an edge the bot gave itself.
  const fairProb = (1 + row.edgePct / 100) * Math.min(unabatedAllIn, kelly.americanToProb(row.price));
  const decision = {
    reason: null, detail: null,
    ticker: contract.ticker || null, side: contract.side || null, derivedTicker: contract.derived === true,
    unabatedPrice: row.price, unabatedAllIn: roundPrice(unabatedAllIn), unabatedEdge: row.edgePct / 100, fairProb: roundPrice(fairProb),
    liveAsk: null, liveAllIn: null, liveEdge: null, depth: null, betterBooks: null,
    stake: null, contracts: null, cost: null,
  };
  const done = (reason, detail = null) => ({ ...decision, reason, detail });
  if (contract.reason) return done(contract.reason, contract.detail);
  if (quoteError) return done("no_book", `${contract.ticker}: ${quoteError}`);
  if (!quote) return done("no_book", `nothing rests on the ${contract.side} side of ${contract.ticker}`);

  const liveAllIn = kalshiAllInPrice(quote.ask);
  const liveEdge = fairProb / liveAllIn - 1;
  Object.assign(decision, { liveAsk: quote.ask, liveAllIn: roundPrice(liveAllIn), liveEdge: roundPrice(liveEdge), depth: quote.depth });
  if (liveAllIn > unabatedAllIn + PRICE_TOLERANCE) {
    return done("price_moved", `live ${(liveAllIn * 100).toFixed(2)}c vs Unabated ${(unabatedAllIn * 100).toFixed(2)}c`);
  }
  if (liveAllIn < unabatedAllIn - BETTER_PRICE_TOLERANCE) {
    return done("price_moved", `live ${(liveAllIn * 100).toFixed(2)}c is over ${BETTER_PRICE_TOLERANCE * 100}c better than Unabated's ${(unabatedAllIn * 100).toFixed(2)}c: its fair is stale`);
  }
  if (liveEdge < MIN_LIVE_EDGE) return done("edge_short", `live edge ${(liveEdge * 100).toFixed(2)}%`);

  const liveAmerican = probToAmerican(liveAllIn);
  const check = pricecheck.comparePrice(priceIndex, {
    eventId: row.eventId, periodTypeId: row.periodTypeId, betTypeId: row.betTypeId, sideIndex: row.sideIndex,
    points: row.points, bookId: row.book.id, bookName: row.book.name, price: liveAmerican,
  });
  if (check && check.better.length > 0) {
    decision.betterBooks = check.better;
    return done("not_best", check.better.map((book) => `${book.bookName} ${book.price}`).join(", "));
  }
  if (liveEdge > MAX_LIVE_EDGE) return done("cap", `live edge ${(liveEdge * 100).toFixed(2)}%`);

  const advice = betsView.stakeAdvice({
    line: row, price: liveAmerican, edgePct: liveEdge * 100, bankroll: BOT_BANKROLL, multiplier: BOT_KELLY,
    matches: sizing.matches, ladderOf: sizing.ladderOf, liquidity: null, teasers: sizing.teasers,
  });
  if (advice.kind === "declined") return done("size_declined", advice.reason);
  const stake = betsView.suggestedBetAmount(advice);
  const contracts = Math.floor(Math.min(stake / liveAllIn, quote.depth));
  Object.assign(decision, { stake, contracts, cost: roundPrice(contracts * liveAllIn) });
  if (contracts < 1) return done("size_zero", `stake $${stake.toFixed(2)} at ${(liveAllIn * 100).toFixed(2)}c, depth ${quote.depth}`);
  return done("bet");
}

// Words of a team name, lowercase, without the filler Kalshi and Unabated
// spell differently ("Washington St." / "Washington State Cougars").
const TEAM_FILLER_WORDS = new Set(["st", "state", "university", "univ", "the", "of", "u"]);
function teamWords(name) {
  return String(name || "").toLowerCase().split(/[^a-z0-9]+/).filter((word) => word && !TEAM_FILLER_WORDS.has(word));
}

// Does Kalshi's market's YES title name the row's game's team on `sideIndex`?
// Its first word must be a word of that team and not of the other team
// ("Los Angeles R" against "Los Angeles C" fails: neither side is sure).
//   market  Kalshi's GET /markets/<ticker> `market` object (yes_sub_title)
function titleNamesSide(market, row, sideIndex) {
  const title = market && (market.yes_sub_title || market.title);
  const first = teamWords(title)[0];
  if (!first) return false;
  const named = sideIndex === 0 ? row.awayTeam : row.homeTeam;
  const other = sideIndex === 0 ? row.homeTeam : row.awayTeam;
  return teamWords(named).includes(first) && !teamWords(other).includes(first);
}

// Does Kalshi's market for a derived moneyline name the row's team?
function moneylineTeamMatches(market, row) {
  return titleNamesSide(market, row, row.sideIndex);
}

// Why Kalshi's market is not the one the row prices, or null when it is:
// the market must be active; a spread's or total's floor_strike must be the
// row's number ("wins by over 35.5" for -35.5, "Over 59.5" for 59.5); a
// spread's or moneyline's YES title must name the team YES is on (the row's
// own on a Y- contract, the other on an N- one). Read before every order so
// a ticker that drifted from the row (a re-keyed rung, a closed market) is
// never paid for.
//   market    Kalshi's GET /markets/<ticker> `market` object
//   contract  resolveKalshiContract's {side}
function marketMismatch(market, row, contract) {
  if (!market || typeof market !== "object") return "no market object";
  if (market.status !== MARKET_STATUS_ACTIVE) return `status "${market.status}"`;
  if (row.betTypeId === BET_TYPE_MONEYLINE) {
    return titleNamesSide(market, row, row.sideIndex) ? null : `"${market.yes_sub_title}" is not ${row.sideLabel}`;
  }
  const strike = Number(market.floor_strike);
  if (!Number.isFinite(strike) || strike !== Math.abs(row.points)) return `floor_strike ${market.floor_strike}, row is ${row.points}`;
  if (row.betTypeId === BET_TYPE_SPREAD) {
    const yesSide = contract.side === "yes" ? row.sideIndex : 1 - row.sideIndex;
    if (!titleNamesSide(market, row, yesSide)) return `"${market.yes_sub_title}" does not name side ${yesSide}'s team`;
  } else if (row.betTypeId !== BET_TYPE_TOTAL) {
    return `bet type ${row.betTypeId} has no market check`;
  }
  return null;
}

// The order's client_order_id: "ubot-" (every bot order carries it, Cal's
// tag) and a hash of ticker, side and pass, so one pass can never send the
// same order twice (the desk refuses a repeat id).
function clientOrderIdOf(ticker, side, passStartedMs) {
  const digest = crypto.createHash("sha256").update(`${ticker}|${side}|${passStartedMs}`).digest("hex");
  return `ubot-${digest.slice(0, 32)}`;
}

// What makes a decision worth logging again: its outcome and the prices it
// was taken at. The same row deciding the same way every pass logs once.
function decisionSignature(decision) {
  const orderStatus = decision.order ? decision.order.status : "";
  return [decision.reason, decision.ticker, decision.liveAsk, decision.unabatedPrice, decision.contracts, orderStatus].join("|");
}

module.exports = {
  KALSHI_BOOK_ID, BOT_BANKROLL, BOT_KELLY, MIN_LIVE_EDGE, MAX_LIVE_EDGE, PRICE_TOLERANCE, BETTER_PRICE_TOLERANCE,
  MAX_LINE_AGE_MS, KALSHI_TAKER_FEE_RATE, REASONS,
  kalshiAllInPrice, probToAmerican, unabatedAllInOf, candidateRows,
  contractFromRungs, spreadTeamCodes, derivedMoneylineContract, resolveKalshiContract,
  parseKalshiOrderbook, quoteOf, decideRow, decisionSignature, moneylineTeamMatches, marketMismatch, clientOrderIdOf,
};
