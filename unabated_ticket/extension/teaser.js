// Buckeye 6-point teasers for the Unabated Ticket panel's Teasers tab
// (docs/2026-09-27-unabated-ticket-teasers-plan.md). Pure: no DOM, no fetch,
// no chrome.* — loaded as a plain <script> in panel.html after bets.js and
// ladder.js (exposes globalThis.UnabatedTeaser) and via require() in
// tests/teaser.test.js. Nothing here writes anywhere.
//
// Inputs
//   feedState   the scanner's state (feed.js): events and every book's lines.
//   records     bet records (bets.js contract). BFA is Buckeye, so an open BFA
//               teaser is a placed ticket. It arrives as one record per leg:
//               isParlayLeg, parlayId, legCount, the ticket's stake and toWin
//               on every leg, the leg's teased points, rotation and
//               eventStart, raw.headerDescription "4 TEAM TEASERS".
//   boardLines  one feed.describeLine row per board event (the panel's
//               boardLines()): what bets.js joins a leg to its game on.
// Outputs
//   teaserBoardOf one pass over the board: Buckeye's lines and each game's
//                 fair ladder (the panel redoes it only when NFL or CFB reloads)
//   teaserLegs    every Buckeye side on the board teased 6 points, with its
//                 win chance or the reason it has none
//   openTeasers   the open BFA teasers as tickets, each leg joined to its game
//   planTeasers   the tickets to bet, kept while nothing real changes
//   describePlan  those tickets' EVs and the summary at the current fairs
//   describeLegs  the Legs list: the best leg of each game and its standing
//
// The rules (user decisions 2026-09-27/28):
//   - A leg is Buckeye's main full-game spread or total, NFL or CFB, teased
//     6 points: spread +6, Over -6, Under +6. Buckeye's price is ignored
//     (-120 and +100 tease the same).
//   - A push loses, so a leg wins exactly when the half-point one step
//     against it wins: Bills -1 is priced at -1.5, Browns +8 at +7.5, Over
//     38 at 38.5. The fair is Unabated's ladder at that half-point, built
//     from the event's lines of that bet type alone: the moneyline stays out
//     of the spread ladder (Unabated's moneyline fair and its alt spreads
//     disagree by ~2.5 points near zero).
//   - Tickets are 4 teams at +300: all four legs win -> +3 per dollar,
//     anything else -> -1.
//   - One leg per game, so legs are independent. The pool is each game's
//     best leg, the top POOL_SIZE by win chance; every 4-leg combination of
//     the pool is a candidate ticket.
//   - The set maximizes E[ln(1 + pnl / K)], K = bankroll x multiplier: add
//     the $200 ticket that most raises it, repeat until none does; the last
//     may be partial (whole dollars); no ticket twice. Open teasers count as
//     part of the set: an open ticket's own four legs are never offered
//     again, and once an open 4-team ticket under $200 is in play (the
//     list's partial, placed) no second partial is.
//   - An open BFA teaser is fixed P&L in that objective (+toWin when every
//     leg wins, -stake otherwise) until its last game starts. Its legs join
//     their games through bets.js's matcher, else by rotation within 12 h of
//     BFA's start. A started leg,
//     a leg no board game matches and a leg with no fair count as won: the
//     ticket keeps its full weight on the legs still to play, and no fair
//     has to be saved. A game with an open leg offers new legs only on that
//     leg's market, and its legs at different numbers share rows — the game
//     is one outcome cut at every number, as in condkelly.js.
//   - The list holds still: it is built on reference fairs, each leg's fair
//     as of the build, replaced only when the leg moves REBUILD_FAIR_MOVE.
//     The list is rebuilt only when what it is built on changes: the pool,
//     a pool leg's number, an open teaser, K, or a reference fair. Smaller
//     ticks move the EVs in place (describePlan). Built on the same fairs,
//     the list after placing its first ticket is the same list less that
//     ticket.

(function (root) {
  "use strict";

  const inNode = typeof module !== "undefined" && module.exports;
  const feed = inNode ? require("./feed.js") : root.UnabatedFeed;
  const ladderLib = inNode ? require("./ladder.js") : root.UnabatedLadder;
  const bets = inNode ? require("./bets.js") : root.UnabatedBets;

  // Unabated's market source id for Buckeye, the book the BFA account bets at.
  const BUCKEYE_BOOK_ID = 59;
  // NFL and CFB: the scanner loads them whatever the Edges tab shows.
  const TEASER_LEAGUE_IDS = [1, 2];
  const FULL_GAME_PERIOD_ID = 1;
  const FULL_GAME = "FG";
  const BET_TYPE_SPREAD = 2;
  const BET_TYPE_TOTAL = 3;
  const STATUS_ON_BOARD = 1;
  const SIDE_AWAY_OR_OVER = 0;
  const TEASER_POINTS = 6;
  const LEGS_PER_TICKET = 4;
  // +300: a winning ticket returns $3 profit per $1 staked.
  const TICKET_NET_ODDS = 3;
  // Buckeye's limit per ticket (user, 2026-09-27).
  const TICKET_MAX_STAKE = 200;
  // 2^10 = 1,024 joint outcomes (26 ms on the 9/27 board); 12 legs took 221 ms and added nothing.
  const POOL_SIZE = 10;
  // One percentage point of a leg's win chance.
  const REBUILD_FAIR_MOVE = 0.01;
  // The win chance four equal legs need for +300 to break even: 0.25^(1/4) = 70.7%.
  const BREAK_EVEN_WIN = (1 / (1 + TICKET_NET_ODDS)) ** (1 / LEGS_PER_TICKET);
  const HALF_POINT = 0.5;
  // A value strictly inside the gap between two half-point cuts.
  const QUARTER_POINT = 0.25;
  // condkelly.js's tolerance: Unabated's fairs are whole American prices, so
  // neighbouring rungs can cross by a hair; past this the ladder is wrong.
  const MAX_NEGATIVE_SLICE = 0.005;
  // The joint outcomes one build may enumerate: the pool's 1,024 and room for
  // four more games open teasers still ride on (past it the pool sheds legs).
  // A build at 2^16 took 0.8-1.5 s on the panel's thread; 2^14 stays under ~0.3 s.
  const MAX_OUTCOMES = 1 << 14;
  const GOLDEN_RATIO = (Math.sqrt(5) - 1) / 2;
  const SEARCH_MAX_ITERATIONS = 200;
  const SEARCH_TOLERANCE_DOLLARS = 1e-7;
  // Stay strictly inside the stake that would take the worst outcome's wealth to zero.
  const UPPER_BOUND_SHRINK = 1 - 1e-9;
  const BFA_VENUE = "bfa";
  const TEASER_HEADER_RE = /TEASER/i;
  // How far BFA's start may sit from the board's for a rotation-only join:
  // wide enough for a shifted clock (7 h once), narrow enough that next
  // week's game with the same rotation never fits (bets.js's own
  // "start time differs" window).
  const ROTATION_JOIN_WINDOW_MS = 12 * 3600 * 1000;
  // An open teaser's leg: still to play and priced (in the math), or counted as won.
  const LEG_LIVE = "live";
  const LEG_STARTED = "started";
  const LEG_OFF_BOARD = "off_board";
  const LEG_UNPRICED = "unpriced";
  // Where a game's leg stands in the Legs list.
  const STANDING_POOL = "pool";
  const STANDING_OUT = "out";
  const STANDING_BELOW = "below_break_even";
  const STANDING_OTHER_MARKET = "other_market";
  const STANDING_UNPRICED = "unpriced";
  const REASON_FEW_LEGS = `fewer than ${LEGS_PER_TICKET} games with a priced Buckeye leg`;
  const REASON_HELD_RISK = "open teasers can already lose the Kelly bankroll";
  const REASON_NO_STAKE = "no stake on the record";

  // ---- small helpers -------------------------------------------------------

  function signedPoints(points) {
    if (points === 0) return "pk";
    return points > 0 ? `+${points}` : `${points}`;
  }

  function americanText(price) {
    return price > 0 ? `+${price}` : `${price}`;
  }

  function isPositiveNumber(value) {
    return typeof value === "number" && Number.isFinite(value) && value > 0;
  }

  function compareText(a, b) {
    return a < b ? -1 : a > b ? 1 : 0;
  }

  function cutKey(eventId, axis, cut) {
    return `${eventId}|${axis}|${cut}`;
  }

  function factorKey(eventId, axis) {
    return `${eventId}|${axis}`;
  }

  function winOf(probAbove, direction) {
    return direction === "above" ? probAbove : 1 - probAbove;
  }

  // ---- legs ------------------------------------------------------------------

  // Spread +6 on the side's own number; a total moves 6 toward the side.
  function teasedPointsOf(betTypeId, sideIndex, points) {
    if (betTypeId === BET_TYPE_SPREAD) return points + TEASER_POINTS;
    return sideIndex === SIDE_AWAY_OR_OVER ? points - TEASER_POINTS : points + TEASER_POINTS;
  }

  // A whole number pushes on itself and a push loses, so the leg wins
  // exactly when the half-point one step against it wins.
  function halfPointAgainst(cut, direction) {
    if (!Number.isInteger(cut)) return cut;
    return direction === "above" ? cut + HALF_POINT : cut - HALF_POINT;
  }

  // Where a teased leg sits on its market's axis and when it wins:
  // {axis, direction, winCut}. condkelly.js's convention: an away spread at
  // points a wins above -a, a home spread at h below h, an Over above its
  // number, an Under below. winCut is always a half-point.
  function legPosition(betTypeId, sideIndex, teasedPoints) {
    const direction = sideIndex === SIDE_AWAY_OR_OVER ? "above" : "below";
    if (betTypeId === BET_TYPE_TOTAL) {
      return { axis: ladderLib.AXIS_TOTAL, direction, winCut: halfPointAgainst(teasedPoints, direction) };
    }
    // `|| 0` folds a teased pick'em's -0 into 0.
    const cut = (sideIndex === SIDE_AWAY_OR_OVER ? -teasedPoints : teasedPoints) || 0;
    return { axis: ladderLib.AXIS_MARGIN, direction, winCut: halfPointAgainst(cut, direction) };
  }

  // The number the leg is priced at, in its own terms: Bills -1 -> -1.5,
  // Browns +8 -> +7.5, Over 38 -> 38.5, Under 44 -> 43.5.
  function pricedPointsOf(betTypeId, sideIndex, teasedPoints) {
    if (!Number.isInteger(teasedPoints)) return teasedPoints;
    if (betTypeId === BET_TYPE_TOTAL && sideIndex === SIDE_AWAY_OR_OVER) return teasedPoints + HALF_POINT;
    return teasedPoints - HALF_POINT;
  }

  // Buckeye's main full-game spread or total on an NFL or CFB board, with a number.
  function isBuckeyeTeaserLine(line) {
    return Boolean(line) && line.bookId === BUCKEYE_BOOK_ID && !line.isAlt
      && line.periodTypeId === FULL_GAME_PERIOD_ID
      && (line.betTypeId === BET_TYPE_SPREAD || line.betTypeId === BET_TYPE_TOTAL)
      && TEASER_LEAGUE_IDS.includes(line.leagueId) && line.statusId === STATUS_ON_BOARD
      && typeof line.points === "number" && Number.isFinite(line.points);
  }

  // One pass over the board, done once per NFL or CFB snapshot (the panel
  // keeps it while other leagues refresh; a pass over NFL + CFB is ~100 ms):
  // Buckeye's candidate lines, and ladderOf(eventId, axis) — Unabated's fair
  // ladder for that game and market, built on first use from the game's
  // full-game lines of that bet type alone (no moneyline in the spread
  // ladder). A game outside NFL and CFB has no ladder: its legs read "no fair".
  //   {lines, ladderOf}
  function teaserBoardOf(feedState) {
    const lines = [];
    const spreadAndTotalLinesByEvent = new Map();
    for (const line of Object.values(feedState.lines)) {
      if (!TEASER_LEAGUE_IDS.includes(line.leagueId) || line.periodTypeId !== FULL_GAME_PERIOD_ID) continue;
      if (line.betTypeId !== BET_TYPE_SPREAD && line.betTypeId !== BET_TYPE_TOTAL) continue;
      if (!spreadAndTotalLinesByEvent.has(line.eventId)) spreadAndTotalLinesByEvent.set(line.eventId, []);
      spreadAndTotalLinesByEvent.get(line.eventId).push(line);
      if (isBuckeyeTeaserLine(line)) lines.push(line);
    }
    const ladders = new Map();
    const ladderOf = (eventId, axis) => {
      const key = factorKey(eventId, axis);
      if (!ladders.has(key)) {
        const betTypeId = axis === ladderLib.AXIS_TOTAL ? BET_TYPE_TOTAL : BET_TYPE_SPREAD;
        const marketLines = (spreadAndTotalLinesByEvent.get(eventId) || []).filter((line) => line.betTypeId === betTypeId);
        ladders.set(key, ladderLib.buildLadder(marketLines, { periodTypeId: FULL_GAME_PERIOD_ID, axis }));
      }
      return ladders.get(key);
    };
    return { lines, ladderOf };
  }

  // ladder.probAbove's reason code in words, at the number the leg is priced at.
  function noFairText(reason, pricedPoints, isSpread) {
    const at = isSpread ? signedPoints(pricedPoints) : `${pricedPoints}`;
    return reason === ladderLib.REASON_FLAT ? `Unabated's fair is flat at ${at}` : `no Unabated fair at ${at}`;
  }

  function legLabel(betTypeId, sideIndex, points, awayTeam, homeTeam) {
    if (betTypeId === BET_TYPE_TOTAL) return `${sideIndex === SIDE_AWAY_OR_OVER ? "Over" : "Under"} ${points}`;
    const team = sideIndex === SIDE_AWAY_OR_OVER ? awayTeam : homeTeam;
    return `${team || "?"} ${signedPoints(points)}`;
  }

  function legOf(line, feedState, readLadder) {
    const described = feed.describeLine(line, feedState);
    const isSpread = line.betTypeId === BET_TYPE_SPREAD;
    const teasedPoints = teasedPointsOf(line.betTypeId, line.sideIndex, line.points);
    const position = legPosition(line.betTypeId, line.sideIndex, teasedPoints);
    const pricedPoints = pricedPointsOf(line.betTypeId, line.sideIndex, teasedPoints);
    const fair = ladderLib.probAbove(readLadder(line.eventId, position.axis), position.winCut);
    return {
      key: `${line.eventId}:bt${line.betTypeId}:si${line.sideIndex}`,
      eventId: line.eventId,
      leagueId: line.leagueId,
      leagueLabel: described.leagueLabel,
      betTypeId: line.betTypeId,
      sideIndex: line.sideIndex,
      axis: position.axis,
      direction: position.direction,
      winCut: position.winCut,
      teasedPoints,
      pricedPoints,
      bookPoints: line.points,
      bookPrice: line.price,
      label: legLabel(line.betTypeId, line.sideIndex, teasedPoints, described.awayTeam, described.homeTeam),
      fromLabel: `from ${isSpread ? signedPoints(line.points) : line.points}`,
      bookLabel: `${isSpread ? signedPoints(line.points) : line.points} ${americanText(line.price)}`,
      matchup: `${described.awayTeam || "?"} @ ${described.homeTeam || "?"}`,
      rotation: described.rotation,
      eventStartMs: described.eventStartMs,
      modifiedMs: described.modifiedMs,
      probAbove: fair.reason ? null : fair.prob,
      win: fair.reason ? null : winOf(fair.prob, position.direction),
      reason: fair.reason ? noFairText(fair.reason, pricedPoints, isSpread) : null,
    };
  }

  function byWinThenKey(a, b) {
    return (b.win ?? -1) - (a.win ?? -1) || compareText(a.key, b.key);
  }

  // Every Buckeye side on the board that can be teased right now: on the
  // board, game not started, and — with maxLineAgeMs — changed by Buckeye
  // within that window (the Edges "Max line age h": a dead feed keeps old
  // numbers). Priced legs carry `win`; the others `reason`. Best first.
  //   options  {now, maxLineAgeMs, board?} — `board` is teaserBoardOf(feedState),
  //            made here when not passed
  function teaserLegs(feedState, options) {
    const { now, maxLineAgeMs } = options;
    const board = options.board || teaserBoardOf(feedState);
    const legs = [];
    for (const line of board.lines) {
      const event = feedState.events[line.eventId];
      if (!event || event.eventStart == null || event.eventStart <= now) continue;
      if (maxLineAgeMs != null) {
        const changedMs = feed.lineChangedMs(line);
        if (changedMs == null || now - changedMs > maxLineAgeMs) continue;
      }
      legs.push(legOf(line, feedState, board.ladderOf));
    }
    return legs.sort(byWinThenKey);
  }

  // ---- open BFA teasers --------------------------------------------------------

  function isOpenBfaTeaserLeg(record) {
    if (!record || record.venue !== BFA_VENUE || record.status !== "open") return false;
    if (!record.isParlayLeg || !record.parlayId) return false;
    // The open list names the ticket in headerDescription, the history in type.
    const raw = record.raw || {};
    return TEASER_HEADER_RE.test(String(raw.headerDescription || raw.type || ""));
  }

  // bet id -> the board row of its game, through bets.js's matcher (rotation
  // + start: BFA rotations are Unabated's own), else by rotation alone within
  // ROTATION_JOIN_WINDOW_MS (rowByRotation). A leg on two board games, or
  // none, has no row.
  function boardRowsByBetId(legRecords, boardLines) {
    const rows = boardLines || [];
    const rowByBetId = new Map();
    bets.annotateRows(rows, legRecords, { lines: rows }).forEach((flag, index) => {
      for (const match of flag.matches) rowByBetId.set(match.bet.id, rows[index]);
    });
    for (const record of legRecords) {
      if (rowByBetId.has(record.id)) continue;
      const row = rowByRotation(record, rows);
      if (row) rowByBetId.set(record.id, row);
    }
    return rowByBetId;
  }

  // The one board game of the leg's league whose away or home rotation is
  // the leg's and whose start is within ROTATION_JOIN_WINDOW_MS of BFA's,
  // or null. A leg the matcher could not place — BFA's clock more than 30
  // min off the board's (it moved once already), a team name resolving to
  // another team — would otherwise count as won, and the ticket it is on
  // would come back into the list as if never placed.
  function rowByRotation(record, rows) {
    const startMs = record.eventStart ? Date.parse(record.eventStart) : NaN;
    if (record.rotation == null || !Number.isFinite(startMs)) return null;
    const found = new Map();
    for (const row of rows) {
      if (row.league !== record.league || row.eventId == null || typeof row.eventStartMs !== "number") continue;
      if (row.awayRotation !== record.rotation && row.homeRotation !== record.rotation) continue;
      if (Math.abs(row.eventStartMs - startMs) > ROTATION_JOIN_WINDOW_MS) continue;
      found.set(row.eventId, row);
    }
    return found.size === 1 ? found.values().next().value : null;
  }

  // "San Francisco 49ers -1.5" on the board's spelling when the leg joined a
  // game and its side is known, else BFA's own ("SAN FRANCISCO 49ERS -1.5").
  function placedLegLabel(record, row, sideIndex) {
    if (record.betType === "total") return `${record.side === "over" ? "Over" : "Under"} ${record.points}`;
    const boardTeam = row && sideIndex != null ? (sideIndex === SIDE_AWAY_OR_OVER ? row.awayTeam : row.homeTeam) : null;
    const team = boardTeam || record.awayTeam || record.homeTeam || "?";
    return typeof record.points === "number" ? `${team} ${signedPoints(record.points)}` : team;
  }

  function counted(base, state, note) {
    return { ...base, state, note: `${note}; counted as won` };
  }

  // One leg of an open teaser: live (priced, still to play) or counted as won.
  function placedLegOf(record, row, now, ladderOf) {
    const base = { record, row, eventId: row ? row.eventId : null, label: placedLegLabel(record, row, null), win: null };
    if (!row) return counted(base, LEG_OFF_BOARD, record.unmatchable ? `not matched (${record.unmatchable})` : "no board game");
    if (!(row.eventStartMs > now)) return counted(base, LEG_STARTED, "started");
    const position = bets.teaserLegPositionOf(record, row);
    if (position.reason) return counted(base, LEG_UNPRICED, position.reason);
    const sideIndex = position.direction === "above" ? SIDE_AWAY_OR_OVER : 1;
    const withSide = { ...base, label: placedLegLabel(record, row, sideIndex), matchup: `${row.awayTeam || "?"} @ ${row.homeTeam || "?"}` };
    if (position.period !== FULL_GAME) return counted(withSide, LEG_UNPRICED, `${position.period} leg`);
    const winCut = halfPointAgainst(position.cut, position.direction);
    const fair = ladderLib.probAbove(ladderOf(row.eventId, position.axis), winCut);
    if (fair.reason) {
      const pricedPoints = record.betType === "total" ? winCut : (sideIndex === SIDE_AWAY_OR_OVER ? -winCut : winCut);
      return counted(withSide, LEG_UNPRICED, noFairText(fair.reason, pricedPoints, record.betType !== "total"));
    }
    return {
      ...withSide, state: LEG_LIVE, note: null, axis: position.axis, direction: position.direction, winCut,
      probAbove: fair.prob, win: winOf(fair.prob, position.direction),
    };
  }

  function liveLegsOf(ticket) {
    return ticket.legs.filter((leg) => leg.state === LEG_LIVE);
  }

  function placedTicketOf(id, legRecords, rowByBetId, now, ladderOf) {
    const ordered = legRecords.slice().sort((a, b) => (a.legIndex ?? 0) - (b.legIndex ?? 0));
    const first = ordered[0];
    const legs = ordered.map((record) => placedLegOf(record, rowByBetId.get(record.id) || null, now, ladderOf));
    const hasDollars = isPositiveNumber(first.stake) && isPositiveNumber(first.toWin);
    const ticket = {
      id, stake: first.stake, toWin: first.toWin, legCount: first.legCount ?? ordered.length, placedAt: first.placedAt ?? null, legs,
      // In the objective while it has its dollars and a leg still to play.
      inPlay: hasDollars && legs.some((leg) => leg.state === LEG_LIVE),
      reason: hasDollars ? null : REASON_NO_STAKE,
    };
    // Legs counted as won leave the ticket riding on the live ones alone.
    ticket.winAll = ticket.inPlay ? liveLegsOf(ticket).reduce((product, leg) => product * leg.win, 1) : null;
    return ticket;
  }

  // Every open BFA teaser as a ticket {id, stake, toWin, legCount, placedAt,
  // legs, inPlay, winAll, reason}, oldest first. `winAll` is the chance
  // every live leg covers. options {now, ladderOf}.
  function openTeasers(records, boardLines, options) {
    const { now, ladderOf } = options;
    const legsByTicket = new Map();
    for (const record of records || []) {
      if (!isOpenBfaTeaserLeg(record)) continue;
      if (!legsByTicket.has(record.parlayId)) legsByTicket.set(record.parlayId, []);
      legsByTicket.get(record.parlayId).push(record);
    }
    const rowByBetId = boardRowsByBetId(Array.from(legsByTicket.values()).flat(), boardLines);
    return Array.from(legsByTicket, ([id, legRecords]) => placedTicketOf(id, legRecords, rowByBetId, now, ladderOf))
      .sort((a, b) => compareText(a.placedAt || "", b.placedAt || "") || compareText(a.id, b.id));
  }

  // ---- what a build is made of ---------------------------------------------------

  // cut key -> P(result > cut) at every number a leg is priced at now.
  function liveValuation(legs, placed) {
    const valuation = new Map();
    for (const leg of legs) if (leg.probAbove != null) valuation.set(cutKey(leg.eventId, leg.axis, leg.winCut), leg.probAbove);
    for (const ticket of placed) {
      for (const leg of liveLegsOf(ticket)) valuation.set(cutKey(leg.eventId, leg.axis, leg.winCut), leg.probAbove);
    }
    return valuation;
  }

  // The reference fairs: each number's fair as the last build had it, unless
  // it has since moved REBUILD_FAIR_MOVE or more (then the live one); a new
  // number takes its live fair.
  function referenceValuation(refs, live) {
    const valuation = new Map();
    for (const [key, liveProb] of live) {
      const ref = refs ? refs.get(key) : undefined;
      valuation.set(key, ref !== undefined && Math.abs(liveProb - ref) < REBUILD_FAIR_MOVE ? ref : liveProb);
    }
    return valuation;
  }

  // eventId -> the markets (axes) open teasers still ride on in that game.
  function openMarketsByEvent(tickets) {
    const markets = new Map();
    for (const ticket of tickets) {
      for (const leg of liveLegsOf(ticket)) {
        if (!markets.has(leg.eventId)) markets.set(leg.eventId, new Set());
        markets.get(leg.eventId).add(leg.axis);
      }
    }
    return markets;
  }

  // A game with an open leg offers new legs on that leg's market only (none
  // when open legs sit on both).
  function isOffered(leg, openMarkets) {
    const markets = openMarkets.get(leg.eventId);
    return !markets || (markets.size === 1 && markets.has(leg.axis));
  }

  function byPoolOrder(a, b) {
    return b.win - a.win || compareText(a.leg.key, b.leg.key);
  }

  // What a build is made of, on one valuation: the pool ({leg, win} per
  // game, the top POOL_SIZE), the open teasers in play, and the valuation.
  function selectInputs(legs, placed, valuation) {
    const inPlay = placed.filter((ticket) => ticket.inPlay);
    const openMarkets = openMarketsByEvent(inPlay);
    const bestByEvent = new Map();
    for (const leg of legs) {
      if (leg.probAbove == null || !isOffered(leg, openMarkets)) continue;
      const candidate = { leg, win: winOf(valuation.get(cutKey(leg.eventId, leg.axis, leg.winCut)), leg.direction) };
      const held = bestByEvent.get(leg.eventId);
      if (!held || byPoolOrder(candidate, held) < 0) bestByEvent.set(leg.eventId, candidate);
    }
    const pool = Array.from(bestByEvent.values()).sort(byPoolOrder).slice(0, POOL_SIZE);
    return { pool, placed: inPlay, valuation };
  }

  // Everything a build is made of, as text: equal keys build equal lists.
  function inputKey(inputs, kellyBankroll) {
    return JSON.stringify({
      kellyBankroll,
      pool: inputs.pool.map(({ leg, win }) => [leg.key, leg.winCut, win]),
      placed: inputs.placed.map((ticket) => [ticket.id, ticket.stake, ticket.toWin,
        liveLegsOf(ticket).map((leg) => [leg.eventId, leg.axis, leg.winCut, leg.direction, inputs.valuation.get(cutKey(leg.eventId, leg.axis, leg.winCut))])]),
    });
  }

  // ---- joint outcomes ------------------------------------------------------------
  //
  // Every game a pool leg or a live open leg is on is one factor, cut at every
  // number a leg there is priced at; a factor's rows are the gaps between its
  // cuts, and games are independent. An outcome picks one row per factor
  // (mixed radix: row = floor(outcome / stride) % rows).

  // The rows of a factor a leg wins: row r holds the results between
  // cuts[r - 1] and cuts[r].
  function rowsWon(cuts, leg) {
    const won = new Uint8Array(cuts.length + 1);
    for (let row = 0; row <= cuts.length; row += 1) {
      const value = row === 0 ? cuts[0] - QUARTER_POINT : cuts[row - 1] + QUARTER_POINT;
      won[row] = (leg.direction === "above" ? value > leg.winCut : value < leg.winCut) ? 1 : 0;
    }
    return won;
  }

  function outcomesWon(factor, leg, outcomeCount) {
    const won = rowsWon(factor.cuts, leg);
    const wins = new Uint8Array(outcomeCount);
    for (let outcome = 0; outcome < outcomeCount; outcome += 1) wins[outcome] = won[Math.floor(outcome / factor.stride) % factor.rows];
    return wins;
  }

  // {factors, outcomeCount, poolWins, placedWins}: which outcomes each pool
  // leg wins, and which each open ticket wins (all of its live legs).
  function outcomeStructure(poolLegs, placedTickets) {
    const byKey = new Map();
    const addCut = (leg) => {
      const key = factorKey(leg.eventId, leg.axis);
      if (!byKey.has(key)) byKey.set(key, { eventId: leg.eventId, axis: leg.axis, matchup: leg.matchup, cuts: new Set() });
      byKey.get(key).cuts.add(leg.winCut);
    };
    poolLegs.forEach(addCut);
    for (const ticket of placedTickets) liveLegsOf(ticket).forEach(addCut);
    let outcomeCount = 1;
    const factors = Array.from(byKey.values()).map((factor) => {
      const cuts = Array.from(factor.cuts).sort((a, b) => a - b);
      const shaped = { eventId: factor.eventId, axis: factor.axis, matchup: factor.matchup, cuts, rows: cuts.length + 1, stride: outcomeCount };
      outcomeCount *= shaped.rows;
      return shaped;
    });
    const structure = { factors, outcomeCount, poolWins: [], placedWins: [] };
    if (outcomeCount > MAX_OUTCOMES) return structure;
    const factorOf = new Map(factors.map((factor) => [factorKey(factor.eventId, factor.axis), factor]));
    const winsOf = (leg) => outcomesWon(factorOf.get(factorKey(leg.eventId, leg.axis)), leg, outcomeCount);
    structure.poolWins = poolLegs.map(winsOf);
    structure.placedWins = placedTickets.map((ticket) => {
      const wins = new Uint8Array(outcomeCount).fill(1);
      for (const leg of liveLegsOf(ticket)) {
        const legWins = winsOf(leg);
        for (let outcome = 0; outcome < outcomeCount; outcome += 1) wins[outcome] &= legWins[outcome];
      }
      return wins;
    });
    return structure;
  }

  // One factor's row chances off the valuation: {probs} or {reason}.
  function factorRowProbs(factor, valuation) {
    const above = factor.cuts.map((cut) => valuation.get(cutKey(factor.eventId, factor.axis, cut)));
    const missing = factor.cuts.find((cut, index) => typeof above[index] !== "number");
    if (missing !== undefined) {
      throw new Error(`teaser: expected a fair at ${missing} on event ${factor.eventId} (${factor.axis}), found none`);
    }
    const edges = [1, ...above, 0];
    const probs = [];
    for (let row = 0; row < factor.rows; row += 1) {
      const slice = edges[row] - edges[row + 1];
      if (slice < -MAX_NEGATIVE_SLICE) {
        return { reason: `Unabated's fair on ${factor.matchup} is higher at ${factor.cuts[row]} than at ${factor.cuts[row - 1]} (${factor.axis})` };
      }
      probs.push(Math.max(0, slice));
    }
    // Clamped slices leave the total a hair over 1.
    const total = probs.reduce((sum, prob) => sum + prob, 0);
    return { probs: probs.map((prob) => prob / total) };
  }

  // {prob}: each outcome's chance, or {reason} when a ladder is not monotone.
  function outcomeProbs(structure, valuation) {
    const rowProbs = [];
    for (const factor of structure.factors) {
      const rows = factorRowProbs(factor, valuation);
      if (rows.reason) return rows;
      rowProbs.push(rows.probs);
    }
    const prob = new Float64Array(structure.outcomeCount);
    for (let outcome = 0; outcome < structure.outcomeCount; outcome += 1) {
      let chance = 1;
      structure.factors.forEach((factor, index) => {
        chance *= rowProbs[index][Math.floor(outcome / factor.stride) % factor.rows];
      });
      prob[outcome] = chance;
    }
    return { prob };
  }

  // The open tickets' P&L in every outcome: +toWin when every live leg wins, -stake otherwise.
  function placedPnl(structure, placedTickets) {
    const pnl = new Float64Array(structure.outcomeCount);
    placedTickets.forEach((ticket, index) => {
      const wins = structure.placedWins[index];
      for (let outcome = 0; outcome < pnl.length; outcome += 1) pnl[outcome] += wins[outcome] ? ticket.toWin : -ticket.stake;
    });
    return pnl;
  }

  // ---- the greedy ------------------------------------------------------------------

  // Every `size`-leg combination of `count` pool legs, in lexicographic order.
  function legCombinations(count, size) {
    const combos = [];
    const picked = [];
    function choose(start) {
      if (picked.length === size) {
        combos.push(picked.slice());
        return;
      }
      for (let index = start; index < count; index += 1) {
        picked.push(index);
        choose(index + 1);
        picked.pop();
      }
    }
    choose(0);
    return combos;
  }

  // The outcomes a ticket wins: all of its legs win.
  function winningOutcomes(combo, poolWins, outcomeCount) {
    const indexes = [];
    for (let outcome = 0; outcome < outcomeCount; outcome += 1) {
      let allWin = 1;
      for (const leg of combo) allWin &= poolWins[leg][outcome];
      if (allWin) indexes.push(outcome);
    }
    return Int32Array.from(indexes);
  }

  function membership(indexes, outcomeCount) {
    const member = new Uint8Array(outcomeCount);
    for (const outcome of indexes) member[outcome] = 1;
    return member;
  }

  // K plus the P&L of the worst outcome that can happen.
  function worstWealth(pnl, prob, kellyBankroll) {
    let worst = Infinity;
    for (let outcome = 0; outcome < pnl.length; outcome += 1) if (prob[outcome] > 0) worst = Math.min(worst, pnl[outcome]);
    return kellyBankroll + worst;
  }

  // E[ln(1 + pnl / K)]; every outcome's wealth must be positive.
  function expectedLogGrowth(pnl, prob, kellyBankroll) {
    let total = 0;
    for (let outcome = 0; outcome < pnl.length; outcome += 1) {
      if (prob[outcome] > 0) total += prob[outcome] * Math.log(1 + pnl[outcome] / kellyBankroll);
    }
    return total;
  }

  // The unused ticket that raises the growth most at `stake`: {index,
  // growth}, or null. Each outcome's lose term is shared by every ticket, so
  // a ticket's growth is their sum plus its own winning outcomes' gains.
  function bestTicketAt(stake, winning, used, pnl, prob, kellyBankroll) {
    const gain = new Float64Array(pnl.length);
    let loseTotal = 0;
    for (let outcome = 0; outcome < pnl.length; outcome += 1) {
      if (prob[outcome] === 0) continue;
      const lose = prob[outcome] * Math.log(1 + (pnl[outcome] - stake) / kellyBankroll);
      const win = prob[outcome] * Math.log(1 + (pnl[outcome] + stake * TICKET_NET_ODDS) / kellyBankroll);
      loseTotal += lose;
      gain[outcome] = win - lose;
    }
    let best = null;
    for (let index = 0; index < winning.length; index += 1) {
      if (used[index]) continue;
      let growth = loseTotal;
      for (const outcome of winning[index]) growth += gain[outcome];
      if (!best || growth > best.growth) best = { index, growth };
    }
    return best;
  }

  // The growth is concave in the stake, so a bounded golden-section search finds its top.
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

  // The partial last ticket: the whole-dollar stake up to `upper` that
  // maximizes the growth, or 0 when a dollar does not raise it above `base`.
  function partialStake(winningIndexes, pnl, prob, kellyBankroll, base, upper) {
    const wins = membership(winningIndexes, pnl.length);
    const growthAt = (stake) => {
      let total = 0;
      for (let outcome = 0; outcome < pnl.length; outcome += 1) {
        if (prob[outcome] === 0) continue;
        const ticketPnl = wins[outcome] ? stake * TICKET_NET_ODDS : -stake;
        total += prob[outcome] * Math.log(1 + (pnl[outcome] + ticketPnl) / kellyBankroll);
      }
      return total;
    };
    const stake = Math.floor(goldenSectionMax(growthAt, upper));
    return stake >= 1 && growthAt(stake) > base ? stake : 0;
  }

  // Greedy Kelly over the set: {tickets: [{legIndexes, stake}], reason,
  // partialHeldBack} — the last the partial ticket's stake when one is due
  // but not allowed, else 0.
  //   takenCombos     "i,j,k,l" pool-index keys already open at BFA: no ticket twice
  //   partialAllowed  false once this list's partial is open (one partial per set)
  function greedyTickets(structure, prob, heldPnl, kellyBankroll, { takenCombos, partialAllowed }) {
    const combos = legCombinations(structure.poolWins.length, LEGS_PER_TICKET);
    const winning = combos.map((combo) => winningOutcomes(combo, structure.poolWins, structure.outcomeCount));
    const used = new Uint8Array(combos.length);
    combos.forEach((combo, index) => {
      if (takenCombos.has(combo.join(","))) used[index] = 1;
    });
    const pnl = Float64Array.from(heldPnl);
    const tickets = [];
    const add = (index, stake) => {
      used[index] = 1;
      tickets.push({ legIndexes: combos[index], stake });
      const wins = membership(winning[index], pnl.length);
      for (let outcome = 0; outcome < pnl.length; outcome += 1) pnl[outcome] += wins[outcome] ? stake * TICKET_NET_ODDS : -stake;
    };
    for (;;) {
      const room = worstWealth(pnl, prob, kellyBankroll);
      if (room <= 0) return { tickets, reason: tickets.length ? null : REASON_HELD_RISK, partialHeldBack: 0 };
      const base = expectedLogGrowth(pnl, prob, kellyBankroll);
      // A full ticket only while losing it everywhere keeps every wealth positive.
      const cap = Math.min(TICKET_MAX_STAKE, room * UPPER_BOUND_SHRINK);
      const full = cap === TICKET_MAX_STAKE;
      const best = bestTicketAt(full ? TICKET_MAX_STAKE : cap / 2, winning, used, pnl, prob, kellyBankroll);
      if (!best) break;
      if (full && best.growth > base) {
        add(best.index, TICKET_MAX_STAKE);
        continue;
      }
      const stake = partialStake(winning[best.index], pnl, prob, kellyBankroll, base, cap);
      if (!partialAllowed) return { tickets, reason: null, partialHeldBack: stake };
      if (stake > 0) add(best.index, stake);
      break;
    }
    return { tickets, reason: null, partialHeldBack: 0 };
  }

  function legIdentity(leg) {
    return `${leg.eventId}|${leg.axis}|${leg.winCut}|${leg.direction}`;
  }

  // The pool combinations open teasers already are — a 4-leg ticket whose
  // legs are all live and exactly four pool legs (same game, market, number
  // and side) — as "i,j,k,l" keys. The greedy never offers one again:
  // without this, a ticket placed at BFA came back at the top of the list on
  // a board with standout legs, and following the list stacked one
  // combination 4 times. A 5-team ticket with a leg started is not one.
  function openComboKeys(poolLegs, placedTickets) {
    const poolIndexByLeg = new Map(poolLegs.map((leg, index) => [legIdentity(leg), index]));
    const keys = new Set();
    for (const ticket of placedTickets) {
      if (ticket.legs.length !== LEGS_PER_TICKET || liveLegsOf(ticket).length !== LEGS_PER_TICKET) continue;
      const indexes = liveLegsOf(ticket).map((leg) => poolIndexByLeg.get(legIdentity(leg)));
      if (indexes.includes(undefined)) continue;
      const distinct = Array.from(new Set(indexes)).sort((a, b) => a - b);
      if (distinct.length === LEGS_PER_TICKET) keys.add(distinct.join(","));
    }
    return keys;
  }

  // The joint outcomes within MAX_OUTCOMES: while over it, the pool gives up
  // its weakest leg on a game no open teaser rides on (dropping one on an open
  // game saves nothing, its rows stay). {legs, structure}, or {overBudget}
  // when even four pool legs do not fit.
  function structureWithinBudget(poolLegs, placedTickets) {
    const openGames = new Set();
    for (const ticket of placedTickets) for (const leg of liveLegsOf(ticket)) openGames.add(factorKey(leg.eventId, leg.axis));
    let legs = poolLegs;
    for (;;) {
      const structure = outcomeStructure(legs, placedTickets);
      if (structure.outcomeCount <= MAX_OUTCOMES) return { legs, structure };
      let drop = -1;
      for (let index = legs.length - 1; index >= 0 && drop < 0; index -= 1) {
        if (!openGames.has(factorKey(legs[index].eventId, legs[index].axis))) drop = index;
      }
      if (drop < 0 || legs.length <= LEGS_PER_TICKET) return { overBudget: structure };
      legs = legs.filter((_, index) => index !== drop);
    }
  }

  // One build: {key, refs, kellyBankroll, pool, placed, structure, tickets,
  // reason, partialHeldBack}. `tickets` index into `pool`; `refs` is the
  // valuation it was built on; `partialHeldBack` the partial ticket's stake
  // when one is due but this list's partial is already open, else 0.
  function buildTeaserSet(inputs, kellyBankroll, key) {
    const build = {
      key, refs: inputs.valuation, kellyBankroll, pool: inputs.pool.map(({ leg }) => leg), placed: inputs.placed,
      structure: null, tickets: [], reason: null, partialHeldBack: 0,
    };
    if (build.pool.length < LEGS_PER_TICKET) return { ...build, reason: REASON_FEW_LEGS };
    const fitted = structureWithinBudget(build.pool, inputs.placed);
    if (fitted.overBudget) {
      return { ...build, reason: `open teasers ride on too many games to size with ${LEGS_PER_TICKET} pool legs (${fitted.overBudget.factors.length} games, ${fitted.overBudget.outcomeCount} outcomes)` };
    }
    const pool = fitted.legs;
    const probs = outcomeProbs(fitted.structure, inputs.valuation);
    if (probs.reason) return { ...build, pool, reason: probs.reason };
    // One partial ticket per set: an open 4-team ticket under $200 still in
    // play is taken to be this list's partial, placed, and no second one is
    // offered — else placing it brought a smaller one, and so on. Another
    // kind of teaser under $200 (a 2-team) holds nothing back.
    const partialAllowed = !inputs.placed.some((ticket) => ticket.legCount === LEGS_PER_TICKET && ticket.stake < TICKET_MAX_STAKE);
    const picked = greedyTickets(fitted.structure, probs.prob, placedPnl(fitted.structure, inputs.placed), kellyBankroll,
      { takenCombos: openComboKeys(pool, inputs.placed), partialAllowed });
    return {
      ...build, pool, structure: fitted.structure, tickets: picked.tickets, reason: picked.reason, partialHeldBack: picked.partialHeldBack,
    };
  }

  // The list to show: `previous` kept while what it is built on is unchanged
  // (its reference fairs carried forward), else a new build on the
  // reference fairs. Returns {build, rebuilt}.
  //   legs          teaserLegs(...) now
  //   placed        openTeasers(...) now
  //   kellyBankroll K = bankroll x multiplier
  //   previous      the last build, or null
  function planTeasers({ legs, placed, kellyBankroll, previous }) {
    const valuation = referenceValuation(previous ? previous.refs : null, liveValuation(legs, placed));
    const inputs = selectInputs(legs, placed, valuation);
    const key = inputKey(inputs, kellyBankroll);
    // Carrying the valuation forward freezes a number that is not in the
    // list too, so a leg ticking around the pool's edge cannot flip it in and out.
    if (previous && previous.key === key) return { build: { ...previous, refs: valuation }, rebuilt: false };
    return { build: buildTeaserSet(inputs, kellyBankroll, key), rebuilt: true };
  }

  // ---- at the current fairs ------------------------------------------------------------

  // The build's tickets and summary at the current fairs: the list holds
  // still while the win chances, EVs and summary move with every tick.
  //   {tickets: [{number, legs, stake, winAll, ev}], summary, reason}
  //   summary {count, stake, expected, placedCount, placedStake, expectedAll,
  //            makesMoney, allLose} — the last three over the tickets to bet
  //            and every open teaser in play; null without a build to size.
  function describePlan(build, legs, placed) {
    const liveLegByKey = new Map(legs.map((leg) => [leg.key, leg]));
    const poolNow = build.pool.map((leg) => liveLegByKey.get(leg.key) || leg);
    const tickets = build.tickets.map((ticket, index) => {
      const ticketLegs = ticket.legIndexes.map((legIndex) => poolNow[legIndex]);
      const winAll = ticketLegs.reduce((product, leg) => product * leg.win, 1);
      return { number: index + 1, legs: ticketLegs, stake: ticket.stake, winAll, ev: winAll * (1 + TICKET_NET_ODDS) - 1 };
    });
    return { tickets, summary: summaryOf(build, tickets, legs, placed), reason: build.reason };
  }

  function summaryOf(build, tickets, legs, placed) {
    const stake = tickets.reduce((sum, ticket) => sum + ticket.stake, 0);
    const expected = tickets.reduce((sum, ticket) => sum + ticket.stake * ticket.ev, 0);
    const placedStake = build.placed.reduce((sum, ticket) => sum + ticket.stake, 0);
    const summary = { count: tickets.length, stake, expected, placedCount: build.placed.length, placedStake, expectedAll: null, makesMoney: null, allLose: null };
    if (!build.structure || tickets.length + build.placed.length === 0) return summary;
    // Live fairs where the numbers are still priced, the build's own elsewhere.
    const valuation = new Map([...build.refs, ...liveValuation(legs, placed)]);
    const probs = outcomeProbs(build.structure, valuation);
    const prob = probs.prob || outcomeProbs(build.structure, build.refs).prob;
    const pnl = placedPnl(build.structure, build.placed);
    const anyWin = new Uint8Array(build.structure.outcomeCount);
    build.structure.placedWins.forEach((wins) => wins.forEach((won, outcome) => { anyWin[outcome] |= won; }));
    for (const ticket of build.tickets) {
      const wins = membership(winningOutcomes(ticket.legIndexes, build.structure.poolWins, build.structure.outcomeCount), pnl.length);
      for (let outcome = 0; outcome < pnl.length; outcome += 1) {
        pnl[outcome] += wins[outcome] ? ticket.stake * TICKET_NET_ODDS : -ticket.stake;
        anyWin[outcome] |= wins[outcome];
      }
    }
    let expectedAll = 0;
    let makesMoney = 0;
    let allLose = 0;
    for (let outcome = 0; outcome < pnl.length; outcome += 1) {
      expectedAll += prob[outcome] * pnl[outcome];
      if (pnl[outcome] > 0) makesMoney += prob[outcome];
      if (!anyWin[outcome]) allLose += prob[outcome];
    }
    return { ...summary, expectedAll, makesMoney, allLose };
  }

  // The Legs list: one row per game — its pool leg, else its best priced leg
  // (on the open teasers' market when the game has one), else a leg with
  // the reason it has no fair. Priced rows best first, then the unpriced.
  //   rows [{leg, standing, inTickets, openTickets, note}]
  function describeLegs(legs, build, placed) {
    const inPlay = placed.filter((ticket) => ticket.inPlay);
    const openMarkets = openMarketsByEvent(inPlay);
    const poolIndexByKey = new Map(build.pool.map((leg, index) => [leg.key, index]));
    const ticketsByPoolIndex = new Map();
    for (const ticket of build.tickets) for (const index of ticket.legIndexes) ticketsByPoolIndex.set(index, (ticketsByPoolIndex.get(index) || 0) + 1);
    const openTicketsByEvent = new Map();
    for (const ticket of inPlay) {
      for (const eventId of new Set(liveLegsOf(ticket).map((leg) => leg.eventId))) openTicketsByEvent.set(eventId, (openTicketsByEvent.get(eventId) || 0) + 1);
    }
    const legsByEvent = new Map();
    for (const leg of legs) {
      if (!legsByEvent.has(leg.eventId)) legsByEvent.set(leg.eventId, []);
      legsByEvent.get(leg.eventId).push(leg);
    }
    const rows = [];
    for (const [eventId, eventLegs] of legsByEvent) {
      const openTickets = openTicketsByEvent.get(eventId) || 0;
      const poolLeg = eventLegs.find((leg) => poolIndexByKey.has(leg.key));
      const priced = eventLegs.filter((leg) => leg.win != null).sort(byWinThenKey);
      const offered = priced.filter((leg) => isOffered(leg, openMarkets));
      if (poolLeg) {
        const inTickets = ticketsByPoolIndex.get(poolIndexByKey.get(poolLeg.key)) || 0;
        rows.push({ leg: poolLeg, standing: STANDING_POOL, inTickets, openTickets, note: null });
      } else if (offered.length) {
        const best = offered[0];
        const standing = best.win >= BREAK_EVEN_WIN ? STANDING_OUT : STANDING_BELOW;
        rows.push({ leg: best, standing, inTickets: 0, openTickets, note: standing === STANDING_OUT ? `not in the top ${build.pool.length}` : "below break-even" });
      } else if (priced.length) {
        rows.push({ leg: priced[0], standing: STANDING_OTHER_MARKET, inTickets: 0, openTickets, note: "open teasers use the game's other market" });
      } else {
        rows.push({ leg: eventLegs[0], standing: STANDING_UNPRICED, inTickets: 0, openTickets, note: eventLegs[0].reason });
      }
    }
    return rows.sort((a, b) => {
      const pricedA = a.leg.win != null;
      const pricedB = b.leg.win != null;
      if (pricedA !== pricedB) return pricedA ? -1 : 1;
      return pricedA ? byWinThenKey(a.leg, b.leg) : (a.leg.eventStartMs - b.leg.eventStartMs) || compareText(a.leg.key, b.leg.key);
    });
  }

  const api = {
    TEASER_LEAGUE_IDS, TEASER_POINTS, LEGS_PER_TICKET, TICKET_NET_ODDS, TICKET_MAX_STAKE, POOL_SIZE, BREAK_EVEN_WIN,
    LEG_LIVE, LEG_STARTED,
    STANDING_POOL, STANDING_OUT, STANDING_BELOW, STANDING_OTHER_MARKET, STANDING_UNPRICED,
    REASON_FEW_LEGS,
    teaserBoardOf, teaserLegs, openTeasers, planTeasers, describePlan, describeLegs,
  };

  if (inNode) {
    module.exports = api;
  } else {
    root.UnabatedTeaser = api;
  }
})(typeof globalThis !== "undefined" ? globalThis : this);
