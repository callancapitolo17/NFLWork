// Unabated Ticket server — the automated taker's decision, per Edges row
// (plan: automated_taker_plan.md). Every decision ends in a reason; a "bet"
// decision is sent by the runner to the bets service's order desk
// (bets_service/taker_orders.py), which holds the bot's Kalshi key, the stop
// file and the $500 ceiling. Live from merge (Cal, 2026-10-14: no shadow run).
//
// Cal's rules (2026-10-01 and 2026-10-14), each a named constant below:
//   - Kalshi rows Unabated shows an edge on, every league and market the
//     settings scan (alts included, moneylines and main numbers through a derived ticker);
//   - fair = Unabated's own (`ge` on the row), no cross-check;
//   - Kalshi's live all-in price (ask + the 7% taker fee, the way Unabated
//     shows Kalshi) at most 0.1c worse than Unabated's shown price, and at
//     most 1c better: a price that ran away from Unabated's since its
//     snapshot is news Unabated's fair has not caught up with, not an edge;
//   - live edge at least 3% and at most 20%;
//   - no live book beats the price at that number (pricecheck.js, strictly:
//     one better book is enough to skip, a thin one under pricecheck's $30
//     floor too; outliers and Cal's Skips trade);
//   - lines changed within 7 days, right up to the start;
//   - $2,000 at 0.25 Kelly (Cal, 2026-10-16), conditional on every open bet (betsview.js stakeAdvice),
//     contracts cut to the depth at the ask; the order desk caps the money at risk.
//
// Novig (2026-10-16, automated_taker_plan.md § Novig; Cal: its own $500 sub-account,
// $2,000 at 0.25 Kelly, the better all-in price first when both books trade a number):
// the same rules on Novig's rows (book 89), with Novig's own facts. The
// contract is the Novig outcome Unabated's rung names (feed.js novigOutcomes;
// a main number with no rung is looked up in Novig's catalog off its
// siblings; Novig moneylines carry no rung, so they never trade); buying it costs 1
// minus the best resting order on the OTHER outcome (every Novig order buys
// an outcome); pregame game markets charge takers nothing, so all-in = ask;
// a contract pays 1c, so the stake buys stake x 100 / ask contracts; the
// market must be OPEN, its event OPEN_PREGAME, its fee not charged pregame,
// and its strike the row's (novigMarketMismatch). The desk is
// bets_service/novig_taker_desk.py.
//
// Reasons, in the order they are checked:
//   stale_feed    the row's league snapshot is stale or its last load failed
//   no_ticker     no Kalshi contract for the row (or an ambiguous one)
//   no_outcome    no Novig outcome for the row (or an ambiguous one), or no
//                 open pregame Novig market lists it
//   no_book       the orderbook read failed or has no ask on the row's side
//   price_moved   live all-in more than PRICE_TOLERANCE worse than Unabated's,
//                 or more than BETTER_PRICE_TOLERANCE better
//   edge_short    live edge under MIN_LIVE_EDGE
//   not_best      another live book beats the live price at the number
//   cap           live edge over MAX_LIVE_EDGE
//   size_declined the held bets refuse the stake (fail closed): conditional
//                 Kelly declined; a bet on the game, on the row's axis, that
//                 the math could not price (unpricedHolding: the stake would
//                 be sized as if it were not held); or the other side of the
//                 same Kalshi market (or Novig market) is held
//   size_zero     the stake buys no contract at the ask, or no depth
//   market_mismatch  at send time, Kalshi's market for the ticker is not
//                 active, or its strike or its YES team is not the row's
//                 (marketMismatch; the runner reads the market before every order)
//   bet           buy `contracts` at `liveAsk` (the runner sends it)
// At send time the order carries the earlier of Unabated's and Kalshi's start
// for the game (kalshiStartMsOf, GET /milestones), and the desk refuses it
// once passed ("started"); Unabated's alone when Kalshi lists none.
//
// Pure: no fetch, no clock (now comes in), no DOM. node-tested.

"use strict";

const crypto = require("node:crypto");
const feed = require("../extension/feed.js");
const kelly = require("../extension/kelly.js");
const betsLib = require("../extension/bets.js");
const betsView = require("../extension/betsview.js");
const pricecheck = require("../extension/pricecheck.js");

const KALSHI_BOOK_ID = 105;
const NOVIG_BOOK_ID = 89;
// What Kalshi's Kelly divides by; the money at risk is the order desk's ceiling (taker_orders.BOT_BANKROLL).
const KELLY_BANKROLL = 2000;
// Novig's Kelly bankroll (Cal, 2026-10-16: $2,000 like Kalshi's), its own parameter, apart from the Novig desk's $500 cap.
const NOVIG_KELLY_BANKROLL = 2000;
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
// A full-game series carries no period token: KXNCAAFSPREAD, not KXNCAAF1HTOTAL.
const FULL_GAME_SERIES_RE = /^KX[A-Z]+(SPREAD|TOTAL|GAME)$/;
const SIDE_AWAY_OR_OVER = 0;
const PRICE_DECIMALS = 1e4;

// Novig: a contract pays 1c, so $1 of payout is 100 contracts.
const NOVIG_QTY_PER_DOLLAR = 100;
const NOVIG_MARKET_OPEN = "OPEN";
const NOVIG_EVENT_PREGAME = "OPEN_PREGAME";
// A market whose fee is charged ALWAYS charges a pregame taker; WHEN_LIVE charges only in-game.
const NOVIG_FEE_ALWAYS = "ALWAYS";
// Novig's market types: the full game's carry the bare name (docs.novig.com GET /v3/types/markets).
const NOVIG_MARKET_TYPE_OF_BET_TYPE = { 1: "MONEY", 2: "SPREAD", 3: "TOTAL" };
// A total outcome's side word: "Over" / "Over 47.5" / "O 47.5" (Novig's app labels a total "O 7.5").
const NOVIG_OVER_NAME_RE = /^(over|o)\b/i;
const NOVIG_UNDER_NAME_RE = /^(under|u)\b/i;
const VENUE_KALSHI = "kalshi";
const VENUE_NOVIG = "novig";

const REASONS = ["stale_feed", "no_ticker", "no_outcome", "no_book", "price_moved", "edge_short", "not_best", "cap", "size_declined", "size_zero", "market_mismatch", "bet"];

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
//   bookId        KALSHI_BOOK_ID (default) or NOVIG_BOOK_ID
function candidateRows(feedState, edgeSettings, now, bookId = KALSHI_BOOK_ID) {
  if (!feedState) return [];
  return feed.selectEdges(feedState, {
    leagueIds: new Set(edgeSettings.leagues),
    periods: new Set(edgeSettings.periods),
    betTypes: new Set(ALL_BET_TYPES),
    bookIds: new Set([bookId]),
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
  if (found.length === 0 && !row.isAlt) return derivedMainNumberContract(row);
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

// A Kalshi main line's contract, which the feed never carries: Unabated's
// ladder leaves out the rung at the main number (the fixture's 52.5 main
// hangs 50.5/51.5/53.5/54.5; Place refused a live Under 43.5 main on
// 2026-10-17), and main lines carry no id. Every rung on the line's OWN
// ladder, on its side of zero, shares its contract side, series, event and
// (spreads) team code, and Kalshi names a strike X.5 as X+1 (51.5 -> "-52",
// -35.5 -> "WSU36"), so the main number's ticker is the siblings' with the
// row's strike. Inferred from Kalshi's naming like the moneyline: `derived`,
// and marketMismatch's ticker, floor_strike, side and team checks are its
// existence check before any order.
//   {ticker, side, sourceKey: null, derived: true} or {reason, detail}
function derivedMainNumberContract(row) {
  const contracts = (row.venueIds && row.venueIds.kalshiContracts) || {};
  const seriesKind = row.betTypeId === BET_TYPE_SPREAD ? "SPREAD" : row.betTypeId === BET_TYPE_TOTAL ? "TOTAL" : null;
  if (!seriesKind) return { reason: "no_ticker", detail: `bet type ${row.betTypeId} has no main-number ticker` };
  // A spread ladder that crosses zero trades the other team's market past it
  // (-1.5 main: -2.5 is "A wins by over 2.5", +0.5 is NO on "B wins by over
  // 0.5"), so only the rungs on the main number's side of zero vouch for it.
  const onRowSide = (target) => seriesKind !== "SPREAD" || (target.points < 0) === (row.points < 0);
  const siblings = Object.entries(contracts).filter(([, target]) => target.mainKey === row.key && onRowSide(target));
  if (siblings.length === 0) return { reason: "no_ticker", detail: "no Kalshi rung names this line, and its ladder has no rungs on its side to derive the main number from" };
  const absPoints = Math.abs(row.points);
  const strikeNumber = absPoints + 0.5;
  if (!Number.isInteger(strikeNumber)) return { reason: "no_ticker", detail: `main number ${row.points} is not a half point` };
  let shared = null;
  for (const [sourceKey] of siblings) {
    const match = KALSHI_CONTRACT_RE.exec(sourceKey);
    if (!match || !match[2].endsWith(seriesKind)) return { reason: "no_ticker", detail: `ladder rung ${sourceKey} is not a Kalshi ${seriesKind} contract` };
    const strike = seriesKind === "SPREAD" ? SPREAD_STRIKE_RE.exec(match[4]) : /^\d+$/.test(match[4]) ? [match[4], ""] : null;
    if (!strike) return { reason: "no_ticker", detail: `ladder rung ${sourceKey} has no readable strike` };
    const shape = { prefix: match[1], series: match[2], eventSuffix: match[3], teamCode: strike[1] };
    if (shared && (shared.prefix !== shape.prefix || shared.series !== shape.series || shared.eventSuffix !== shape.eventSuffix || shared.teamCode !== shape.teamCode)) {
      return { reason: "no_ticker", detail: `ladder rungs disagree: ${shared.prefix}-${shared.series}-${shared.eventSuffix}-${shared.teamCode} and ${sourceKey}` };
    }
    shared = shape;
  }
  const yes = shared.prefix === "Y";
  // YES is "<team> wins by over X": the favorite's side, a negative number.
  if (seriesKind === "SPREAD" && yes !== (row.points < 0)) {
    return { reason: "no_ticker", detail: `main number ${row.points} reads ${shared.prefix} rungs on its own side of zero` };
  }
  if (seriesKind === "TOTAL" && yes !== (row.sideIndex === SIDE_AWAY_OR_OVER)) {
    return { reason: "no_ticker", detail: `${row.sideLabel || `side ${row.sideIndex}`} ladder carries ${shared.prefix} rungs` };
  }
  return { ticker: `${shared.series}-${shared.eventSuffix}-${shared.teamCode}${strikeNumber}`, side: yes ? "yes" : "no", sourceKey: null, derived: true };
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

// Every book that beats the live price, thin ones included: pricecheck's
// $30 floor keeps a thin book off the panel's tags, but the bot still skips
// on it (the floor was added after the taker's rules were set, and Cal has
// not moved it onto the bot).
function booksBeating(check, price) {
  if (!check) return [];
  const yourProb = pricecheck.impliedProb(price);
  const thinBeating = (check.thin || []).filter((book) => pricecheck.impliedProb(book.price) < yourProb);
  return check.better.concat(thinBeating);
}

// What differs per venue in a decision: the all-in price of an ask, how many
// contracts a stake buys, what they cost, and which held bet is the other
// side of the market. A Kalshi contract pays $1, a Novig one 1c.
const VENUE_RULES = {
  [VENUE_KALSHI]: {
    kellyBankroll: KELLY_BANKROLL,
    allInOf: kalshiAllInPrice,
    contractsFor: (stake, allIn, depth) => Math.floor(Math.min(stake / allIn, depth)),
    costOf: (contracts, allIn) => contracts * allIn,
    heldOtherSide,
  },
  [VENUE_NOVIG]: {
    kellyBankroll: NOVIG_KELLY_BANKROLL,
    allInOf: (ask) => ask,
    contractsFor: (stake, allIn, depth) => Math.floor(Math.min((stake * NOVIG_QTY_PER_DOLLAR) / allIn, depth)),
    costOf: (contracts, allIn) => (contracts * allIn) / NOVIG_QTY_PER_DOLLAR,
    heldOtherSide: novigHeldOtherSide,
  },
};

// One row's decision.
//   row         a candidateRows row
//   contract    resolveKalshiContract's (or resolveNovigOutcome's) result
//   quote       quoteOf's (or novigQuoteOf's) result, or null (read failed or empty side)
//   quoteError  why the orderbook read failed, when it did
//   priceIndex  pricecheck.buildPriceIndex over the same feed
//   sizing      {matches (bets.annotateRows for the row), ladderOf, teasers}
//   venue       VENUE_KALSHI (default) or VENUE_NOVIG
// Returns {reason, detail, ...the numbers it got to}.
function decideRow({ row, contract, quote, quoteError, priceIndex, sizing, venue = VENUE_KALSHI }) {
  const rules = VENUE_RULES[venue];
  if (!rules) throw new Error(`decideRow: unknown venue ${venue}`);
  const unabatedAllIn = unabatedAllInOf(row);
  // Unabated's edge is on its price, which it shows twice: the American
  // (-143) and Kalshi's all-in (0.587157 = -142.2). Which one `ge` is measured
  // on is not documented, so the fair is read off the cheaper of the two:
  // the lower fair, never an edge the bot gave itself.
  const fairProb = (1 + row.edgePct / 100) * Math.min(unabatedAllIn, kelly.americanToProb(row.price));
  const decision = {
    reason: null, detail: null, venue,
    // A Novig decision's "ticker" is the outcome it buys (the log's ticker column), its side null.
    ticker: contract.ticker || contract.outcomeId || null, side: contract.side || null, derivedTicker: contract.derived === true,
    marketId: contract.marketId || null,
    unabatedPrice: row.price, unabatedAllIn: roundPrice(unabatedAllIn), unabatedEdge: row.edgePct / 100, fairProb: roundPrice(fairProb),
    liveAsk: null, liveAllIn: null, liveEdge: null, depth: null, betterBooks: null,
    stake: null, contracts: null, cost: null,
  };
  const done = (reason, detail = null) => ({ ...decision, reason, detail });
  if (contract.reason) return done(contract.reason, contract.detail);
  if (quoteError) return done("no_book", `${decision.ticker}: ${quoteError}`);
  if (!quote) return done("no_book", `nothing rests to sell ${venue === VENUE_NOVIG ? "outcome" : `the ${contract.side} side of`} ${decision.ticker}`);

  const liveAllIn = rules.allInOf(quote.ask);
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
  const beating = booksBeating(check, liveAmerican);
  if (beating.length > 0) {
    decision.betterBooks = beating;
    return done("not_best", beating.map((book) => `${book.bookName} ${book.price}`).join(", "));
  }
  if (liveEdge > MAX_LIVE_EDGE) return done("cap", `live edge ${(liveEdge * 100).toFixed(2)}%`);

  const otherSideHeld = rules.heldOtherSide(sizing.matches, contract);
  if (otherSideHeld) return done("size_declined", otherSideHeld);
  const advice = betsView.stakeAdvice({
    line: row, price: liveAmerican, edgePct: liveEdge * 100, bankroll: rules.kellyBankroll, multiplier: BOT_KELLY,
    matches: sizing.matches, ladderOf: sizing.ladderOf, liquidity: null, teasers: sizing.teasers,
  });
  if (advice.kind === "declined") return done("size_declined", advice.reason);
  const blind = unpricedHolding(advice, row);
  if (blind) return done("size_declined", blind);
  const stake = betsView.suggestedBetAmount(advice);
  const contracts = rules.contractsFor(stake, liveAllIn, quote.depth);
  const cost = rules.costOf(contracts, liveAllIn);
  Object.assign(decision, { stake, contracts, cost: roundPrice(cost) });
  if (contracts < 1) return done("size_zero", `stake $${stake.toFixed(2)} at ${(liveAllIn * 100).toFixed(2)}c, depth ${quote.depth}`);
  return done("bet");
}

// Kalshi nets one subaccount per market: a NO buy on a market where YES is
// held there sells that YES instead. The bot trades its own subaccount, so it
// would only net against its own bet there; it still never takes the other
// side of a held bet in any subaccount, Cal's primary included (unchanged
// since 2026-10-14). Kalshi record ids are kalshi:<ticker>:<side>, plus
// :s<N> in numbered subaccount N (bets_service/sources/kalshi.py record_id_of).
//   matches  bets.annotateRows matches for the row ({bet})
function heldOtherSide(matches, contract) {
  const otherSide = contract.side === "yes" ? "no" : "yes";
  const otherId = `kalshi:${contract.ticker}:${otherSide}`;
  const isOtherSide = (id) => id === otherId || (typeof id === "string" && id.startsWith(`${otherId}:s`));
  const held = (matches || []).find((match) => match.bet && isOtherSide(match.bet.id) && match.bet.status === "open");
  if (!held) return null;
  const where = held.bet.id === otherId ? "in the primary account" : "in the bot's subaccount";
  return `the ${otherSide.toUpperCase()} side of ${contract.ticker} is held ${where}: the bot never takes the other side of a held bet`;
}

// A held bet on the game that conditional Kelly left out of the math, when
// that is not by design: another market (#129) and an other-period hedge
// (#130) are left out on purpose; a bet on the row's axis the math could not
// price (a number the fair ladder lacks, a side it could not place, a
// parlay leg) would size the stake as if the bet were not held. Null when
// the math saw everything it should.
//   advice  betsView.stakeAdvice's result (matches annotated with inMath, note, position)
function unpricedHolding(advice, row) {
  const rowPosition = betsLib.linePosition(row);
  const rowAxis = rowPosition && !rowPosition.reason ? rowPosition.axis : null;
  for (const match of advice.matches || []) {
    if (match.inMath) continue;
    const position = match.position && !match.position.reason ? match.position : null;
    const axis = position ? position.axis : axisOfBetType(match.bet && match.bet.betType);
    if (rowAxis !== null && axis != null && axis !== rowAxis) continue;
    if (position && rowAxis !== null && position.period !== rowPosition.period && position.direction !== rowPosition.direction) continue;
    return `held ${match.label || match.bet.id} is not in the sizing math: ${match.note}`;
  }
  const teaserOut = (advice.teaserGroups || []).find((group) => !group.inMath);
  return teaserOut ? `open teaser not in the sizing math: ${teaserOut.note}` : null;
}

function axisOfBetType(betType) {
  if (betType === "total") return betsLib.AXIS_TOTAL;
  if (betType === "spread" || betType === "moneyline") return betsLib.AXIS_MARGIN;
  return null;
}

// Words of a team name, lowercase, without the filler Kalshi and Unabated
// spell differently ("Washington St." / "Washington State Cougars").
// Apostrophes and okinas are dropped, not split on, and accents are folded:
// Kalshi writes "Hawai'i" and "San José St" where Unabated writes "Hawaii"
// and "San Jose State".
const TEAM_FILLER_WORDS = new Set(["st", "state", "university", "univ", "the", "of", "u"]);
const ACCENT_MARKS_RE = /[\u0300-\u036f]/g;
const APOSTROPHES_RE = /['\u2018\u2019\u02bb\u02bc`]/g;
function teamWords(name) {
  const folded = String(name || "").normalize("NFD").replace(ACCENT_MARKS_RE, "").replace(APOSTROPHES_RE, "").toLowerCase();
  return folded.split(/[^a-z0-9]+/).filter((word) => word && !TEAM_FILLER_WORDS.has(word));
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
//   contract  resolveKalshiContract's {side, ticker}
function marketMismatch(market, row, contract) {
  if (!market || typeof market !== "object") return "no market object";
  if (market.status !== MARKET_STATUS_ACTIVE) return `status "${market.status}"`;
  if (market.ticker && market.ticker !== contract.ticker) return `Kalshi answered ${market.ticker} for ${contract.ticker}`;
  const series = String(contract.ticker || "").split("-")[0];
  if (FULL_GAME_SERIES_RE.test(series) !== (row.periodTypeId === PERIOD_FULL_GAME)) {
    return `series ${series} is not the row's period (periodTypeId ${row.periodTypeId})`;
  }
  if (row.betTypeId === BET_TYPE_MONEYLINE) {
    return titleNamesSide(market, row, row.sideIndex) ? null : `"${market.yes_sub_title}" is not ${row.sideLabel}`;
  }
  const strike = Number(market.floor_strike);
  if (!Number.isFinite(strike) || strike !== Math.abs(row.points)) return `floor_strike ${market.floor_strike}, row is ${row.points}`;
  if (row.betTypeId === BET_TYPE_SPREAD) {
    // YES is "<team> wins by over X": the favorite's side, a negative number on the row.
    if ((contract.side === "yes") !== (row.points < 0)) return `${contract.side.toUpperCase()} of a wins-by-over market for a ${row.points} row`;
    const yesSide = contract.side === "yes" ? row.sideIndex : 1 - row.sideIndex;
    if (!titleNamesSide(market, row, yesSide)) return `"${market.yes_sub_title}" does not name side ${yesSide}'s team`;
  } else if (row.betTypeId === BET_TYPE_TOTAL) {
    // YES is "Over X": the row's Over side.
    if ((contract.side === "yes") !== (row.sideIndex === SIDE_AWAY_OR_OVER)) return `${contract.side.toUpperCase()} of an Over market for the ${row.sideLabel} row`;
  } else {
    return `bet type ${row.betTypeId} has no market check`;
  }
  return null;
}

// ---- Novig ------------------------------------------------------------------

// The Novig outcome a row buys, off the rung ids (feed.noteRungVenueIds): the
// one outcome whose listed line IS this row, at the row's number and side.
//   {outcomeId} or {reason: "no_outcome", detail}
function resolveNovigOutcome(row) {
  const outcomes = (row.venueIds && row.venueIds.novigOutcomes) || {};
  const found = Object.entries(outcomes).filter(([, target]) => target.lineKey === row.key);
  if (found.length === 0) return { reason: "no_outcome", detail: "no Novig rung names this line" };
  if (found.length > 1) return { reason: "no_outcome", detail: `${found.length} Novig rungs name this line` };
  const [outcomeId, target] = found[0];
  if (target.points !== row.points || target.sideIndex !== row.sideIndex) {
    return { reason: "no_outcome", detail: `outcome ${outcomeId} is ${target.points} side ${target.sideIndex}, row is ${row.points} side ${row.sideIndex}` };
  }
  return { outcomeId };
}

// outcome id -> {market} over Novig's public catalog listing of open
// pregame markets (GET /v3/public/catalog/markets?eventStatus=OPEN_PREGAME
// pages, each market with its outcomes).
//   markets  [Novig market]
function novigOutcomeIndex(markets) {
  const index = new Map();
  for (const market of markets || []) {
    for (const outcome of (market && market.outcomes) || []) {
      if (outcome && typeof outcome.outcomeId === "string") index.set(outcome.outcomeId, { market });
    }
  }
  return index;
}

// The outcome resolved and placed in its market: {outcomeId, marketId,
// otherOutcomeId} or {reason: "no_outcome", detail}. A market with other than
// two outcomes has no "other side" to read the price off.
function novigContractOf(row, index) {
  const resolved = resolveNovigOutcome(row);
  if (resolved.reason && !row.isAlt && !hasOwnRung(row)) return derivedNovigMainContract(row, index);
  if (resolved.reason) return resolved;
  const listed = index.get(resolved.outcomeId);
  if (!listed) return { reason: "no_outcome", detail: `outcome ${resolved.outcomeId} is on no open pregame Novig market` };
  const outcomes = listed.market.outcomes || [];
  if (outcomes.length !== 2) return { reason: "no_outcome", detail: `market ${listed.market.marketId} has ${outcomes.length} outcomes, expected 2` };
  const other = outcomes.find((outcome) => outcome.outcomeId !== resolved.outcomeId);
  return { outcomeId: resolved.outcomeId, marketId: listed.market.marketId, otherOutcomeId: other.outcomeId };
}

// Does any Novig rung claim this line? Two rungs on it is ambiguous, not missing.
function hasOwnRung(row) {
  const outcomes = (row.venueIds && row.venueIds.novigOutcomes) || {};
  return Object.values(outcomes).some((target) => target.lineKey === row.key);
}

// A Novig main line's outcome, which the feed often lacks: Unabated's ladder
// can leave out the rung at the main number (a 52.5 total hanging only
// 50.5/51.5/53.5/54.5), and Novig's outcome ids are random, so nothing
// derives them from the siblings the way Kalshi's tickers are (#190).
// Instead the siblings' outcomes find Novig's event and market type in the
// catalog listing already read, and the one open market of that event and
// type at the row's strike is the main number's (a spread's strike is the
// home side's handicap). The outcome is picked by its name, the same reading
// novigMarketMismatch checks before any order.
//   {outcomeId, marketId, otherOutcomeId, derived: true} or {reason: "no_outcome", detail}
function derivedNovigMainContract(row, index) {
  if (row.betTypeId !== BET_TYPE_SPREAD && row.betTypeId !== BET_TYPE_TOTAL) {
    return { reason: "no_outcome", detail: "no Novig rung names this line" };
  }
  const outcomes = (row.venueIds && row.venueIds.novigOutcomes) || {};
  const siblings = Object.entries(outcomes).filter(([, target]) => target.mainKey === row.key)
    .map(([outcomeId, target]) => ({ outcomeId, target, listed: index.get(outcomeId) })).filter((sibling) => sibling.listed);
  if (siblings.length === 0) {
    return { reason: "no_outcome", detail: "no Novig rung names this line, and no rung on its ladder is on an open pregame Novig market" };
  }
  // Each sibling must read on Novig the way Unabated lists it: its market's
  // strike at the rung's number from Unabated's home side, its outcome's name
  // the rung's side. A game Novig lists with home and away swapped (a neutral
  // site) fails here; past this point the side is picked by name and number
  // alone, which a swap would pass.
  const siblingNames = new Set();
  for (const { outcomeId, target, listed } of siblings) {
    const outcome = (listed.market.outcomes || []).find((candidate) => candidate.outcomeId === outcomeId);
    const siblingRow = { betTypeId: row.betTypeId, sideIndex: target.sideIndex, points: target.points };
    if (target.sideIndex !== row.sideIndex || Number(listed.market.strike) !== novigStrikeOf(siblingRow) || !novigOutcomeNamesRow(outcome, siblingRow)) {
      return { reason: "no_outcome", detail: `Novig lists rung ${outcomeId} as "${outcome && outcome.name}" at strike ${listed.market.strike}, Unabated as ${target.points} side ${target.sideIndex}` };
    }
    siblingNames.add(novigOutcomeTeamText(outcome));
  }
  if (siblingNames.size !== 1) return { reason: "no_outcome", detail: `the ladder's Novig rungs name ${siblingNames.size} different sides` };
  const [siblingName] = siblingNames;
  const siblingMarkets = siblings.map((sibling) => sibling.listed.market);
  const eventIds = new Set(siblingMarkets.map((market) => market.eventId));
  const marketTypes = new Set(siblingMarkets.map((market) => market.marketType));
  if (eventIds.size !== 1 || marketTypes.size !== 1) {
    return { reason: "no_outcome", detail: `the ladder's Novig rungs disagree: ${eventIds.size} events, ${marketTypes.size} market types` };
  }
  const [eventId] = eventIds;
  const [marketType] = marketTypes;
  const strike = novigStrikeOf(row);
  const candidates = new Map();
  for (const listed of index.values()) {
    const market = listed.market;
    if (market.eventId === eventId && market.marketType === marketType && Number(market.strike) === strike) candidates.set(market.marketId, market);
  }
  if (candidates.size !== 1) {
    return { reason: "no_outcome", detail: `${candidates.size} open Novig ${marketType} markets at strike ${strike} on event ${eventId}` };
  }
  const [market] = candidates.values();
  const marketOutcomes = market.outcomes || [];
  if (marketOutcomes.length !== 2) return { reason: "no_outcome", detail: `market ${market.marketId} has ${marketOutcomes.length} outcomes, expected 2` };
  // The side by its name and number, and by the same team (or Over/Under) words its rungs carry.
  const named = marketOutcomes.filter((outcome) => novigOutcomeNamesRow(outcome, row) && novigOutcomeTeamText(outcome) === siblingName);
  if (named.length !== 1) {
    return { reason: "no_outcome", detail: `${named.length} outcomes of market ${market.marketId} name the row's side (${marketOutcomes.map((outcome) => `"${outcome.name}"`).join(", ")})` };
  }
  const other = marketOutcomes.find((outcome) => outcome.outcomeId !== named[0].outcomeId);
  return { outcomeId: named[0].outcomeId, marketId: market.marketId, otherOutcomeId: other.outcomeId, derived: true };
}

// A Novig market's strike for a row: a total's is the total, a spread's the
// home side's handicap (an away row's -points).
function novigStrikeOf(row) {
  if (row.betTypeId !== BET_TYPE_SPREAD) return row.points;
  return row.sideIndex === SIDE_AWAY_OR_OVER ? -row.points : row.points;
}

// An outcome's name less its number: the team ("Washington St -35.5" ->
// "washington st") or the total's side word ("Over 52.5" -> "over").
function novigOutcomeTeamText(outcome) {
  return String((outcome && outcome.name) || "").trim().toLowerCase().replace(/\s*[+-]?\d+(\.\d+)?$/, "");
}

// Does the outcome's name say it is the row's side? A total's "Over"/"O 47.5"
// or "Under"/"U 47.5"; a spread's last word its handicap, a pick'em never.
function novigOutcomeNamesRow(outcome, row) {
  const name = String((outcome && outcome.name) || "").trim();
  if (row.betTypeId === BET_TYPE_TOTAL) {
    return (row.sideIndex === SIDE_AWAY_OR_OVER ? NOVIG_OVER_NAME_RE : NOVIG_UNDER_NAME_RE).test(name);
  }
  const handicap = novigNameHandicap(name);
  return handicap !== null && handicap !== 0 && handicap === row.points;
}

// ---- Place's moneylines (never the bot's) ------------------------------------

// A half's moneyline can tie, and whether Kalshi's or Novig's two-way half
// winner voids a tie or pays it as a loss is unchecked: bought at Unabated's
// tie-no-bet fair, a contract that loses the tie overpays by about P(tie).
// Refused until checked.
const HALF_MONEYLINE_REFUSAL = "only full-game moneylines place: a half's winner market may pay a tie as a loss";

// A Novig full-game moneyline's outcome, for the Ticket's Place button only
// (Cal, 2026-10-19). Unabated gives no Novig id on a moneyline, so the game's
// other Novig rungs name its event, whose one two-outcome "MONEY" market it is.
// Each outcome is told by the team's name on Novig's spread rungs, else by
// Unabated's team name; one told makes the other the other team.
//   feedState  the board, for each rung's own line (its market and period)
//   {outcomeId, marketId, otherOutcomeId, teamText, derived: true} or {reason: "no_outcome", detail}
function novigMoneylineContract(row, index, feedState) {
  const refuse = (detail) => ({ reason: "no_outcome", detail });
  if (row.betTypeId !== BET_TYPE_MONEYLINE) return refuse(`bet type ${row.betTypeId} is not a moneyline`);
  if (row.periodTypeId !== PERIOD_FULL_GAME) return refuse(HALF_MONEYLINE_REFUSAL);
  const outcomes = (row.venueIds && row.venueIds.novigOutcomes) || {};
  const rungs = Object.entries(outcomes)
    .map(([outcomeId, target]) => ({ outcomeId, target, listed: index.get(outcomeId), line: feedState.lines[target.mainKey] }))
    .filter((rung) => rung.listed && rung.line);
  if (rungs.length === 0) return refuse("no Novig rung on this game is on an open pregame Novig market");
  const eventIds = new Set(rungs.map((rung) => rung.listed.market.eventId));
  if (eventIds.size !== 1) return refuse(`the game's Novig rungs are on ${eventIds.size} Novig events`);
  const [eventId] = eventIds;
  const moneyType = NOVIG_MARKET_TYPE_OF_BET_TYPE[BET_TYPE_MONEYLINE];
  const sideNames = novigSideNamesOf(rungs);
  if (sideNames.reason) return sideNames;
  const candidates = new Map();
  for (const listed of index.values()) {
    const market = listed.market;
    if (market.eventId === eventId && market.marketType === moneyType) candidates.set(market.marketId, market);
  }
  if (candidates.size !== 1) return refuse(`${candidates.size} open Novig ${moneyType} markets on event ${eventId}`);
  const [market] = candidates.values();
  const marketOutcomes = market.outcomes || [];
  if (marketOutcomes.length !== 2) return refuse(`market ${market.marketId} has ${marketOutcomes.length} outcomes, expected 2 (a draw makes it 3-way)`);
  const sides = marketOutcomes.map((outcome) => novigMoneylineSideOf(outcome, row, sideNames.names));
  // Two outcomes on the game's own event: one told, the other is the other team.
  const unknown = sides.indexOf(null);
  if (unknown !== -1 && sides[1 - unknown] !== null) sides[unknown] = 1 - sides[1 - unknown];
  if (!sides.includes(0) || !sides.includes(1)) {
    return refuse(`market ${market.marketId}'s outcomes (${marketOutcomes.map((outcome) => `"${outcome.name}"`).join(", ")}) do not name ${row.awayTeam} and ${row.homeTeam} one each`);
  }
  const chosen = marketOutcomes[sides.indexOf(row.sideIndex)];
  const other = marketOutcomes[sides.indexOf(1 - row.sideIndex)];
  return { outcomeId: chosen.outcomeId, marketId: market.marketId, otherOutcomeId: other.outcomeId, teamText: novigOutcomeTeamText(chosen), derived: true };
}

// How Novig names each Unabated side's team, off the game's spread rungs
// ("Virginia Tech +3.5" on side 0 -> names[0] "virginia tech"). A rung that
// does not read on Novig the way Unabated lists it (a neutral site listed
// home/away swapped) refuses, as do two names for one side or one name for
// both. A side no spread rung names is null. {names: [side0, side1]} or {reason, detail}.
function novigSideNamesOf(rungs) {
  const names = [null, null];
  for (const { outcomeId, target, listed, line } of rungs) {
    // A pick'em rung's name cannot tell the side; the others can.
    if (line.betTypeId !== BET_TYPE_SPREAD || target.points === 0) continue;
    const outcome = (listed.market.outcomes || []).find((candidate) => candidate.outcomeId === outcomeId);
    const rungRow = { betTypeId: BET_TYPE_SPREAD, sideIndex: target.sideIndex, points: target.points };
    if (Number(listed.market.strike) !== novigStrikeOf(rungRow) || !novigOutcomeNamesRow(outcome, rungRow)) {
      return { reason: "no_outcome", detail: `Novig lists rung ${outcomeId} as "${outcome && outcome.name}" at strike ${listed.market.strike}, Unabated as ${target.points} side ${target.sideIndex}` };
    }
    const name = novigOutcomeTeamText(outcome);
    if (!name) return { reason: "no_outcome", detail: `rung ${outcomeId}'s name "${outcome.name}" carries no team` };
    if (names[target.sideIndex] !== null && names[target.sideIndex] !== name) {
      return { reason: "no_outcome", detail: `side ${target.sideIndex} reads on Novig as both "${names[target.sideIndex]}" and "${name}"` };
    }
    names[target.sideIndex] = name;
  }
  if (names[0] !== null && names[0] === names[1]) return { reason: "no_outcome", detail: `both sides read on Novig as "${names[0]}"` };
  return { names };
}

// Which Unabated side a moneyline outcome is (0 away, 1 home), or null: the
// side whose Novig spread name it carries; with no spread name to go on, the
// one team whose name holds every word of the outcome's and the other does
// not (teamWords drops "St"/"State" and the like).
function novigMoneylineSideOf(outcome, row, names) {
  const text = novigOutcomeTeamText(outcome);
  if (text && text === names[0]) return 0;
  if (text && text === names[1]) return 1;
  // Both sides named and this is neither: a third name.
  if (names[0] !== null && names[1] !== null) return null;
  const words = teamWords(text);
  if (!words.length) return null;
  const holds = (team) => { const teamSet = new Set(teamWords(team)); return words.every((word) => teamSet.has(word)); };
  const away = holds(row.awayTeam);
  const home = holds(row.homeTeam);
  if (away === home) return null;
  const side = away ? 0 : 1;
  // A name the spreads give the other side never stands for this one.
  return names[side] !== null && names[side] !== text ? null : side;
}

// The Kalshi event a moneyline's market is listed under, for the Ticket's
// Place button when its team code cannot be read off a full-game spread rung
// (derivedMoneylineContract): the league prefix and event suffix every Kalshi
// rung on the game shares ("KXNCAAF1HTOTAL-26SEP19DUQWSU" -> KXNCAAF,
// 26SEP19DUQWSU) under the game series ("KXNCAAFGAME"). {eventTicker} or {reason: "no_ticker", detail}.
function kalshiMoneylineEventOf(row) {
  const refuse = (detail) => ({ reason: "no_ticker", detail });
  if (row.betTypeId !== BET_TYPE_MONEYLINE) return refuse(`bet type ${row.betTypeId} is not a moneyline`);
  if (row.periodTypeId !== PERIOD_FULL_GAME) return refuse(HALF_MONEYLINE_REFUSAL);
  const contracts = (row.venueIds && row.venueIds.kalshiContracts) || {};
  const prefixes = new Set();
  const suffixes = new Set();
  for (const sourceKey of Object.keys(contracts)) {
    const match = KALSHI_CONTRACT_RE.exec(sourceKey);
    const prefix = match ? /^(KX[A-Z]+?)(1H|2H)?(SPREAD|TOTAL)$/.exec(match[2]) : null;
    // Another series (a quarter, a team total) names no league prefix we can read: skipped.
    if (!prefix) continue;
    prefixes.add(prefix[1]);
    suffixes.add(match[3]);
  }
  if (prefixes.size === 0) return refuse("no Kalshi spread or total rung on this game names its Kalshi event");
  if (prefixes.size !== 1 || suffixes.size !== 1) {
    return refuse(`the game's Kalshi rungs name ${prefixes.size} leagues and ${suffixes.size} events: ${[...prefixes].join(", ")} / ${[...suffixes].join(", ")}`);
  }
  const [prefix] = prefixes;
  const [suffix] = suffixes;
  return { eventTicker: `${prefix}GAME-${suffix}` };
}

// The moneyline's market among Kalshi's markets of that event (GET
// /markets?event_ticker=...): exactly two markets (a tie market makes it
// 3-way, refused), one whose YES title names each team (titleNamesSide);
// YES on the row's team's. marketMismatch re-checks the title before any order.
//   {ticker, side: "yes", sourceKey: null, derived: true} or {reason: "no_ticker", detail}
function kalshiMoneylineFromMarkets(row, markets, eventTicker) {
  const listed = (Array.isArray(markets) ? markets : []).filter((market) => market && market.event_ticker === eventTicker);
  if (listed.length !== 2) {
    return { reason: "no_ticker", detail: `Kalshi lists ${listed.length} markets on ${eventTicker}, expected 2 (one per team)` };
  }
  const named = listed.filter((market) => titleNamesSide(market, row, row.sideIndex));
  const other = listed.filter((market) => titleNamesSide(market, row, 1 - row.sideIndex));
  if (named.length !== 1 || other.length !== 1 || named[0] === other[0]) {
    return { reason: "no_ticker", detail: `${eventTicker}'s markets (${listed.map((market) => `"${market.yes_sub_title}"`).join(", ")}) do not name ${row.awayTeam} and ${row.homeTeam} one each` };
  }
  return { ticker: named[0].ticker, side: "yes", sourceKey: null, derived: true };
}

// Novig's book (GET /v3/public/catalog/markets/<id>/book): {orders: {outcomeId:
// [{orderId, price, qty}] best first}}. Every resting order BUYS its outcome,
// so buying one outcome crosses the best resting order on the other, at 1 -
// its price; the depth is every order resting at that best price, in 1c contracts.
//   {ask, depth} or null when nothing rests on the other outcome
function novigQuoteOf(book, contract) {
  if (!book || typeof book !== "object" || !book.orders || typeof book.orders !== "object") throw new Error("Novig book: expected {orders}");
  const resting = (book.orders[contract.otherOutcomeId] || [])
    .map((order) => ({ price: Number(order.price), qty: Number(order.qty) }))
    .filter((order) => order.price > 0 && order.price < 1 && Number.isInteger(order.qty) && order.qty > 0);
  if (!resting.length) return null;
  const best = Math.max(...resting.map((order) => order.price));
  const depth = resting.filter((order) => order.price === best).reduce((sum, order) => sum + order.qty, 0);
  return { ask: Math.round((1 - best) * 1000) / 1000, depth };
}

// Novig records the outcome held; another open Novig bet on the same market,
// on the other outcome, is the other side of this one (Cal's app account,
// the bot's or Place's subaccount alike).
//   matches  bets.annotateRows matches for the row ({bet})
function novigHeldOtherSide(matches, contract) {
  const held = (matches || []).find((match) => match.bet && match.bet.venue === VENUE_NOVIG && match.bet.status === "open"
    && match.bet.venueIds && match.bet.venueIds.marketId === contract.marketId && match.bet.venueIds.outcomeId !== contract.outcomeId);
  if (!held) return null;
  return `the other outcome of Novig market ${contract.marketId} is held (${held.bet.id}): the bot never takes the other side of a held bet`;
}

// Why Novig's market is not the one the row prices, or null when it is: the
// market OPEN, its event OPEN_PREGAME (an in-game take pays a fee and the
// event's resting orders are voided), its fee not charged pregame, the
// outcome on it, its type the row's market (the bare type is the full
// game's), and its strike the row's number (a spread's strike is the home
// side's handicap; a total's is the total). Read before every order.
//   market  GET /v3/public/catalog/markets/<id>; event  GET .../events/<id>
function novigMarketMismatch(market, event, row, contract) {
  if (!market || typeof market !== "object") return "no market object";
  if (market.marketId !== contract.marketId) return `Novig answered market ${market.marketId} for ${contract.marketId}`;
  if (market.status !== NOVIG_MARKET_OPEN) return `market status "${market.status}"`;
  if (!event || event.eventId !== market.eventId) return `no event ${market.eventId} read`;
  if (event.status !== NOVIG_EVENT_PREGAME) return `event status "${event.status}"`;
  if (!market.fee || market.fee.charged === NOVIG_FEE_ALWAYS) return `the market charges a taker fee pregame (${market.fee && market.fee.charged})`;
  if (!(market.outcomes || []).some((outcome) => outcome.outcomeId === contract.outcomeId)) return `outcome ${contract.outcomeId} is not on the market`;
  const bareType = NOVIG_MARKET_TYPE_OF_BET_TYPE[row.betTypeId];
  if (!bareType || typeof market.marketType !== "string" || !market.marketType.startsWith(bareType)) return `market type ${market.marketType} is not the row's ${row.betType}`;
  if ((market.marketType === bareType) !== (row.periodTypeId === PERIOD_FULL_GAME)) return `market type ${market.marketType} is not the row's period (${row.period})`;
  const strike = Number(market.strike);
  // The outcome's own name decides its side, apart from the rung map that
  // picked it: a swapped rung would otherwise buy the other side at this
  // side's fair. A name that says nothing readable fails closed.
  const name = String((market.outcomes.find((outcome) => outcome.outcomeId === contract.outcomeId) || {}).name || "").trim();
  if (row.betTypeId === BET_TYPE_SPREAD) {
    const homePoints = row.sideIndex === SIDE_AWAY_OR_OVER ? -row.points : row.points;
    if (!Number.isFinite(strike) || strike !== homePoints) return `strike ${market.strike} (home), row is ${row.points} side ${row.sideIndex}`;
    const handicap = novigNameHandicap(name);
    if (handicap === null) return `outcome "${name}" names no handicap`;
    // At a pick'em both outcomes read 0, so the name cannot tell a swapped rung: refused.
    if (handicap === 0) return `outcome "${name}" is a pick'em: its name cannot confirm the side`;
    if (handicap !== row.points) return `outcome "${name}" is ${handicap}, the row is ${row.points}`;
  } else if (row.betTypeId === BET_TYPE_MONEYLINE) {
    // The outcome must still carry the team name it was picked by (novigMoneylineContract).
    if (!contract.teamText || novigOutcomeTeamText({ name }) !== contract.teamText) return `outcome "${name}" is not the row's team ("${contract.teamText}")`;
  } else if (row.betTypeId === BET_TYPE_TOTAL) {
    if (!Number.isFinite(strike) || strike !== row.points) return `strike ${market.strike}, row is ${row.points}`;
    const sideWord = row.sideIndex === SIDE_AWAY_OR_OVER ? "over" : "under";
    const pattern = row.sideIndex === SIDE_AWAY_OR_OVER ? NOVIG_OVER_NAME_RE : NOVIG_UNDER_NAME_RE;
    if (!pattern.test(name)) return `outcome "${name}" is not the row's ${sideWord}`;
  } else {
    return `bet type ${row.betTypeId} has no Novig market check`;
  }
  return null;
}

// A spread outcome name's handicap ("Washington St -35.5" -> -35.5, "PK" or
// a bare 0 -> 0), or null when its last word is no number.
function novigNameHandicap(name) {
  const last = String(name).trim().split(/\s+/).pop() || "";
  if (/^(pk|pick|pick'em)$/i.test(last)) return 0;
  if (!/^[+-]?\d+(\.\d+)?$/.test(last)) return null;
  return Number(last);
}

// The order's clientId: a UUID (Novig requires one) hashed from outcome and
// pass, so one pass can never send the same order twice (the desk refuses a
// repeat). Novig does not dedupe it, so it is a guard of ours only.
function novigClientIdOf(outcomeId, passStartedMs) {
  const hex = crypto.createHash("sha256").update(`nbot|${outcomeId}|${passStartedMs}`).digest("hex");
  return `${hex.slice(0, 8)}-${hex.slice(8, 12)}-4${hex.slice(13, 16)}-8${hex.slice(17, 20)}-${hex.slice(20, 32)}`;
}

// The board line as a bet record names it (bets_service sources), for the
// desk to keep: the bot subaccount's read names its bets with it.
function novigLineOf(row) {
  const betType = String(row.betType || "").toLowerCase();
  const isTotal = row.betTypeId === BET_TYPE_TOTAL;
  const side = isTotal ? (row.sideIndex === SIDE_AWAY_OR_OVER ? "over" : "under") : (row.sideIndex === SIDE_AWAY_OR_OVER ? "away" : "home");
  const idOrNull = (id) => (id == null ? null : String(id));
  return {
    league: row.league, betType, period: row.period, side, points: row.points ?? null,
    awayTeam: row.awayTeam ?? null, homeTeam: row.homeTeam ?? null, awayTeamId: idOrNull(row.awayTeamId), homeTeamId: idOrNull(row.homeTeamId),
  };
}

// Kalshi's own start for a game, read off GET /milestones?related_event_ticker=<event>:
// the earliest start_date among the milestones that name the event (in
// related_event_tickers or primary_event_tickers), as epoch ms, or null when
// none does or none carries a zoned time. A backup clock beside Unabated's
// eventStart (Cal, 2026-10-16), read off a different feed.
//   body         Kalshi's {milestones: [{start_date, related_event_tickers, primary_event_tickers}]}
//   eventTicker  the market's event_ticker ("KXNCAAFSPREAD-26OCT16CLMBPENN")
function kalshiStartMsOf(body, eventTicker) {
  const milestones = body && Array.isArray(body.milestones) ? body.milestones : [];
  let earliest = null;
  for (const milestone of milestones) {
    if (!milestone || typeof milestone !== "object") continue;
    const named = [milestone.related_event_tickers, milestone.primary_event_tickers]
      .some((tickers) => Array.isArray(tickers) && tickers.includes(eventTicker));
    if (!named) continue;
    const start = milestone.start_date;
    // A time without a zone is read nowhere: Unabated's naive times are what failed here.
    if (typeof start !== "string" || !/(?:Z|[+-]\d\d:\d\d)$/.test(start)) continue;
    const startMs = Date.parse(start);
    if (Number.isFinite(startMs) && (earliest === null || startMs < earliest)) earliest = startMs;
  }
  return earliest;
}

// The start the desk refuses on: the earlier of Unabated's and Kalshi's, or
// Unabated's alone when Kalshi has none (Cal, 2026-10-16: option B).
function earliestStartMs(unabatedStartMs, kalshiStartMs) {
  const starts = [unabatedStartMs, kalshiStartMs].filter((ms) => Number.isFinite(ms));
  return starts.length ? Math.min(...starts) : null;
}

// The event a market belongs to: Kalshi's event_ticker, else the ticker less its last part.
function eventTickerOf(market, ticker) {
  if (market && typeof market.event_ticker === "string" && market.event_ticker) return market.event_ticker;
  const cut = String(ticker || "").lastIndexOf("-");
  return cut > 0 ? ticker.slice(0, cut) : null;
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
  // A skipped row logs again when its reason or Unabated's price changes, not on every tick of the live book.
  const isBet = decision.reason === "bet";
  return [decision.reason, decision.ticker, isBet ? decision.liveAsk : "", decision.unabatedPrice, isBet ? decision.contracts : "", orderStatus].join("|");
}

module.exports = {
  KALSHI_BOOK_ID, NOVIG_BOOK_ID, VENUE_KALSHI, VENUE_NOVIG, NOVIG_EVENT_PREGAME, KELLY_BANKROLL, NOVIG_KELLY_BANKROLL, BOT_KELLY, MIN_LIVE_EDGE, MAX_LIVE_EDGE, PRICE_TOLERANCE, BETTER_PRICE_TOLERANCE,
  MAX_LINE_AGE_MS, KALSHI_TAKER_FEE_RATE, REASONS,
  kalshiAllInPrice, probToAmerican, unabatedAllInOf, candidateRows,
  contractFromRungs, derivedMainNumberContract, spreadTeamCodes, derivedMoneylineContract, resolveKalshiContract,
  parseKalshiOrderbook, quoteOf, decideRow, decisionSignature, moneylineTeamMatches, marketMismatch, clientOrderIdOf,
  heldOtherSide, unpricedHolding, kalshiStartMsOf, earliestStartMs, eventTickerOf,
  novigMoneylineContract, kalshiMoneylineEventOf, kalshiMoneylineFromMarkets,
  resolveNovigOutcome, novigOutcomeIndex, novigContractOf, novigQuoteOf, novigHeldOtherSide, novigMarketMismatch,
  novigClientIdOf, novigLineOf,
};
