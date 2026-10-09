// Unabated Ticket — the Edges list's rows, without the DOM: the default and
// stored settings, which books the list is restricted to, the standalone
// stake and its tail-flex rank score, conditional-Kelly sizing against the
// open bets and open BFA teasers (#130), the Min suggested bet gate, the
// sort, the market cards, the tail-flex header, the edge tier, the stake
// rail's words, the edge-move tag's reading and words (#132), and how one
// /bets.json body is applied to the held bet records.
//
// Two callers, one set of numbers: the side panel (panel.js, which keeps only
// the DOM writes, chrome.storage, the timers and the per-scanner-update
// caches) and the headless server runner (server/runner.js, phone page plan
// step 1), which must list exactly what the panel lists for the same
// settings. Moved here verbatim from panel.js on 2026-10-02; behaviour unchanged.
//
// Pure: no DOM, no fetch, no chrome.*, no clock — every function that needs
// "now" takes it. Inputs are the scanner's feed state (feed.js shape: lines,
// events, books, teams), its per-line history (edgemove.js), normalised bet
// records (bets.js contract), the open BFA teasers (teaser.openTeasers, built
// by the caller), the tail-flex measurement (tailflex.measureTailFlex, built by
// the caller) and saved fill fairs (fillfair.fairsByBetId). Outputs plain
// objects / strings; nothing here writes anywhere.
//
// Loaded as a plain <script> in panel.html after tailflex.js, betsview.js and
// fillfair.js (exposes globalThis.UnabatedEdgeRows) and via require() in node.

(function (root) {
  "use strict";

  const inNode = typeof module !== "undefined" && module.exports;
  const feed = inNode ? require("./feed.js") : root.UnabatedFeed;
  const kelly = inNode ? require("./kelly.js") : root.UnabatedKelly;
  const betsLib = inNode ? require("./bets.js") : root.UnabatedBets;
  const betsView = inNode ? require("./betsview.js") : root.UnabatedBetsView;
  const ladderLib = inNode ? require("./ladder.js") : root.UnabatedLadder;
  const edgemove = inNode ? require("./edgemove.js") : root.UnabatedEdgeMove;
  const fillfair = inNode ? require("./fillfair.js") : root.UnabatedFillFair;
  const tailflex = inNode ? require("./tailflex.js") : root.UnabatedTailFlex;

  const DEFAULT_STAKE_SETTINGS = { bankroll: 30000, multiplier: 0.25 };
  // maxLineAgeHours: a "live" book's line unchanged for a week is a dead feed
  // (live 2026-09-10: Buckeye -110 on a 44.5 total, 96 days old, "+36.67%").
  const ALL_LEAGUE_IDS = Object.keys(feed.LEAGUES).map(Number);
  // bookIds undefined (never ticked) = DEFAULT_BOOK_NAMES; null = follow the
  // Unabated selection page.js publishes (all live books until one exists),
  // set by its button; an array = the user's own ticks in the panel.
  // Alt lines (#113) are off until asked for; an alt lists only while
  // Unabated's fair and the book's price are 15-85% (feed.altWithinDepthCap; the 7-point cap
  // went 2026-09-30), and inside that the tail flex (tailflex.js) ranks deep
  // rungs down. minLiquidityToWin
  // $100: an exchange line is listed when its resting money can win $100
  // (feed.liquidityCanWin), so a thin longshot stays and a thin favorite goes. 0 = off.
  // minEdgePct 2.5 since 0.16.1 (user, 2026-10-02; was 1.0).
  const DEFAULT_EDGE_SETTINGS = {
    leagues: ALL_LEAGUE_IDS, periods: [1], betTypes: [1, 2, 3], bookIds: undefined, minEdgePct: 2.5, maxLineAgeHours: 168, sortBy: "edge",
    minStake: 0,
    minLiquidityToWin: 100,
    includeAlts: false,
    // One card per (game, market, side) with its best line; the flat list is the toggle off.
    groupByMarket: true,
  };
  // The books the Edges list starts on until you tick your own (the user's
  // list, 2026-09-15). By NAME, not id: BetOnline Direct, Bookmaker-Internal,
  // Poly US Ing and Polymarket US are listed in the panel but absent from the
  // anonymous feed their ids could be read from. A name the feed does not
  // carry (a book not listed today) simply ticks nothing.
  const DEFAULT_BOOK_NAMES = [
    "Bet105", "BetOnline", "BetOnline Direct", "Bookmaker", "Bookmaker-Internal", "Buckeye", "Kalshi",
    "Novig", "NoVig-Internal", "Poly US Ing", "Polymarket", "Polymarket US", "Prophet Exchange",
    "Underdog Prediction Market",
  ];
  const SORT_KEYS = ["edge", "stake", "start", "exposure"];
  // The list shows the top this many rows or cards.
  const MAX_EDGE_ROWS = 200;
  const BET_TYPE_MONEYLINE = 1;
  const BET_TYPE_SPREAD = 2;

  // ---- settings --------------------------------------------------------------

  // Bankroll and Kelly multiplier as stored: a missing or non-numeric value
  // (or 0) falls back to the default.
  function sanitizeStakeSettings(stored) {
    const source = stored && typeof stored === "object" ? stored : {};
    return {
      bankroll: Number(source.bankroll) || DEFAULT_STAKE_SETTINGS.bankroll,
      multiplier: Number(source.multiplier) || DEFAULT_STAKE_SETTINGS.multiplier,
    };
  }

  // The Edges settings as stored, every field checked; anything missing or
  // invalid is the default.
  function sanitizeEdgeSettings(stored) {
    const base = { ...DEFAULT_EDGE_SETTINGS };
    if (!stored || typeof stored !== "object") return base;
    if (Array.isArray(stored.leagues)) base.leagues = stored.leagues.filter((id) => feed.LEAGUES[id]);
    if (Array.isArray(stored.periods) && stored.periods.length) base.periods = stored.periods.filter((id) => feed.PERIODS[id]);
    if (Array.isArray(stored.betTypes) && stored.betTypes.length) base.betTypes = stored.betTypes.filter((id) => feed.BET_TYPES[id]);
    // A stored null is the Unabated-selection choice; no key at all is the default books.
    if (Array.isArray(stored.bookIds)) base.bookIds = stored.bookIds.filter((id) => Number.isInteger(id));
    else if (stored.bookIds === null) base.bookIds = null;
    if (typeof stored.minEdgePct === "number" && stored.minEdgePct >= 0) base.minEdgePct = stored.minEdgePct;
    if (typeof stored.minStake === "number" && stored.minStake >= 0) base.minStake = stored.minStake;
    if (typeof stored.maxLineAgeHours === "number" && stored.maxLineAgeHours > 0) base.maxLineAgeHours = stored.maxLineAgeHours;
    // The old altMinLiquidity was a stake floor, not a to-win one, so it does not carry over.
    if (typeof stored.minLiquidityToWin === "number" && stored.minLiquidityToWin >= 0) base.minLiquidityToWin = stored.minLiquidityToWin;
    if (SORT_KEYS.includes(stored.sortBy)) base.sortBy = stored.sortBy;
    if (typeof stored.includeAlts === "boolean") base.includeAlts = stored.includeAlts;
    if (typeof stored.groupByMarket === "boolean") base.groupByMarket = stored.groupByMarket;
    return base;
  }

  // ---- feed teams and the bets service ---------------------------------------

  // The snapshot's team list by league path ("nfl" -> [team]), for
  // teams.registerTeams: every snapshot carries Unabated's teams and, per game
  // row, a second spelling of each (#118), so bet records resolve to the
  // board's team ids with no hand-written table (#116).
  function feedTeamsByLeague(feedState) {
    const byLeague = {};
    if (!feedState || !feedState.teamIndex) return byLeague;
    for (const team of Object.values(feedState.teamIndex)) {
      const league = feed.LEAGUES[team.leagueId];
      if (!league) continue;
      (byLeague[league.path] ||= []).push(team);
    }
    return byLeague;
  }

  // One /bets.json body applied to what the caller holds. Throws when the
  // body has no bets array (an error page, a proxy) and then nothing changes.
  //   held     {records, crosswalk, pins, fillFairs} as last applied
  //   payload  the parsed /bets.json body; now  epoch ms (the retention prune)
  // Returns {records, crosswalk, pins, fillFairs, generatedAt, sources}:
  // records merged on native id (betsview.mergeServicePayload: pins applied,
  // team keys filled, pruned); the service's crosswalk and pins are the truth,
  // and a payload without them (an older service) keeps the held ones; saved
  // fill fairs never change, so the served rows are merged into the held
  // ones, not swapped in, and a payload without them keeps the held array
  // (the same object, so a caller can tell nothing changed).
  function applyBetsPayload(held, payload, now) {
    if (!payload || !Array.isArray(payload.bets)) throw new Error("bets.json has no bets array");
    const crosswalk = Array.isArray(payload.crosswalk) ? payload.crosswalk : held.crosswalk;
    const pins = Array.isArray(payload.pins) ? payload.pins : held.pins;
    const records = betsView.mergeServicePayload(held.records, { ...payload, crosswalk, pins }, now);
    const fillFairs = Array.isArray(payload.fillFairs) ? fillfair.mergeFillFairs(held.fillFairs, payload.fillFairs, records) : held.fillFairs;
    return {
      records, crosswalk, pins, fillFairs,
      generatedAt: payload.generatedAt ?? null,
      exclusions: Array.isArray(payload.exclusions) ? payload.exclusions : null,
      sources: payload.sources && typeof payload.sources === "object" ? payload.sources : {},
    };
  }

  // ---- books -----------------------------------------------------------------

  // Every live book in the feed but Unabated's own line, by name.
  function liveBooks(feedState) {
    if (!feedState) return [];
    return Object.values(feedState.books).filter((book) => book.isLive && book.id !== feed.UNABATED_LINE_BOOK_ID)
      .sort((a, b) => a.name.localeCompare(b.name));
  }

  function defaultBookIds(feedState) {
    if (!feedState) return [];
    const wanted = new Set(DEFAULT_BOOK_NAMES);
    return Object.values(feedState.books).filter((book) => wanted.has(book.name)).map((book) => book.id);
  }

  // Which books the list is restricted to: the user's own ticks when they
  // have made any, else the default books until "My Unabated selection" is
  // chosen, which follows `unabatedSelection` (the book ids page.js read off
  // the open Unabated tab, or null when there is none), else every live
  // book. Bet types are the settings' own.
  //   {mode: "custom"|"default"|"unabated"|"all", bookIds: Set|null, betTypeIds: Set, filter: booksFilter}
  function effectiveFilter(edgeSettings, feedState, unabatedSelection, booksFilter) {
    let mode;
    let bookIds;
    if (Array.isArray(edgeSettings.bookIds)) {
      mode = "custom";
      bookIds = new Set(edgeSettings.bookIds);
    } else if (edgeSettings.bookIds === undefined) {
      mode = "default";
      bookIds = new Set(defaultBookIds(feedState));
    } else if (unabatedSelection) {
      mode = "unabated";
      bookIds = new Set(unabatedSelection);
    } else {
      mode = "all";
      bookIds = null;
    }
    return { mode, bookIds, betTypeIds: new Set(edgeSettings.betTypes), filter: booksFilter };
  }

  // Everything selectEdges needs except the edge threshold (list and alerts
  // differ there). leagueIds: the scanner also holds NFL and CFB for the
  // Teasers tab, which the Edges tab lists only when Football is ticked.
  function edgeSelectionOptions(edgeSettings, effective, now) {
    return {
      leagueIds: new Set(edgeSettings.leagues),
      periods: new Set(edgeSettings.periods),
      betTypes: effective.betTypeIds,
      bookIds: effective.bookIds,
      now,
      maxLineAgeMs: edgeSettings.maxLineAgeHours * 3600 * 1000,
      includeAlts: edgeSettings.includeAlts,
      minLiquidityToWin: edgeSettings.minLiquidityToWin,
    };
  }

  // ---- sizing ----------------------------------------------------------------

  // Quarter-Kelly with nothing held, never more than the line's resting
  // liquidity: it sorts "by stake" and feeds the tail-flex rank score.
  function stakeFor(row, stakeSettings) {
    if (row.edgePct == null) return null;
    try {
      const kellyStake = kelly.kellyStakeFromEdge({ bookPrice: row.price, edgePct: row.edgePct, bankroll: stakeSettings.bankroll, multiplier: stakeSettings.multiplier }).stake;
      return betsView.capAtLiquidity(kellyStake, row.liquidity).stake;
    } catch (_error) {
      return null;
    }
  }

  // The standalone stake (sized on Unabated's raw edge, never on the flexed
  // one) and the tail-flex rank score: EV dollars after flex, which picks a
  // card's best line and orders the lines inside it. `measurement` is
  // tailflex.measureTailFlex over the same feed state.
  function withStakeAndRank(row, measurement, stakeSettings) {
    const sized = { ...row, stake: stakeFor(row, stakeSettings) };
    const rank = tailflex.rankOfRow(sized, measurement);
    return { ...sized, rankScore: rank ? rank.score : null };
  }

  // One describeLine-shaped row per event on the board (main lines only):
  // the matcher's ambiguity check, its venue id join (each row carries its
  // event's venueIds) and the unmatched list only need to know which games
  // exist, not every book's price. Callers cache it per feed update (73 ms
  // over the 140,755 NFL + CFB lines of 2026-09-15).
  function boardLines(feedState) {
    if (!feedState) return [];
    const seen = new Set();
    const rows = [];
    for (const line of Object.values(feedState.lines)) {
      if (line.isAlt || seen.has(line.eventId)) continue;
      seen.add(line.eventId);
      rows.push(feed.describeLine(line, feedState));
    }
    return rows;
  }

  // Unabated's fair ladders for sizing against held bets (#130): eventId ->
  // (period, axis) -> ladder.buildLadder result, built from the feed's lines
  // on first use and cached inside the returned reader. Make a new one per
  // feed update; a reader over a null state reads no ladder.
  function createLadderReaders(feedState) {
    let linesByEvent = null;
    const ladders = new Map();
    return (eventId) => (period, axis) => {
      const periodTypeId = ladderLib.periodTypeIdOf(period);
      if (!feedState || eventId == null || periodTypeId == null) return null;
      const cacheKey = `${eventId}|${periodTypeId}|${axis}`;
      if (!ladders.has(cacheKey)) {
        if (!linesByEvent) linesByEvent = ladderLib.groupLinesByEvent(Object.values(feedState.lines));
        ladders.set(cacheKey, ladderLib.buildLadder(linesByEvent.get(eventId), { periodTypeId, axis }));
      }
      return ladders.get(cacheKey);
    };
  }

  // Each row gets `bet` = {tier, matches, advice} from the open bet records:
  // what you hold on that market and the stake sized against it (conditional
  // Kelly, #130), open BFA teasers on the game included. No row is ever
  // hidden for being bet — the edge still being there after you bet it is
  // information, and the stake column carries the top-up.
  //   context  {records, boardLines, stakeSettings, ladderReaderOf: eventId -> (period, axis) -> ladder,
  //             teasers (teaser.openTeasers over the same records and board)}
  function withBetFlags(rows, context) {
    const flags = betsLib.annotateRows(rows, context.records, { lines: context.boardLines });
    const teasers = context.teasers;
    return rows.map((row, index) => {
      const flag = flags[index];
      const advice = betsView.stakeAdvice({
        line: row, price: row.price, edgePct: row.edgePct,
        bankroll: context.stakeSettings.bankroll, multiplier: context.stakeSettings.multiplier,
        matches: flag.matches, ladderOf: context.ladderReaderOf(row.eventId), liquidity: row.liquidity, teasers,
      });
      return { ...row, bet: { tier: flag.tier, matches: advice.matches, advice } };
    });
  }

  // Sort key for "by my exposure": dollars in the math on the market, held or
  // against, straight bets and teaser stakes alike.
  function exposureDollars(row) {
    if (!row.bet) return 0;
    const { held, against, teasers } = row.bet.advice;
    return held + against + (teasers ? teasers.held + teasers.against : 0);
  }

  // Min suggested bet: gates on what the rail says to bet now (the stake
  // sized against what is held), not the standalone size. The list and alerts share it.
  function meetsMinStake(row, minStake) {
    return minStake === 0 || betsView.suggestedBetAmount(row.bet.advice) >= minStake;
  }

  // The selected lines at `minEdgePct` and up, each with its standalone
  // `stake`, its `rankScore` and its `bet` flag, in selectEdges' order (edge,
  // then start).
  //   context  {edgeSettings, stakeSettings, effective (effectiveFilter), minEdgePct, measurement,
  //             records, boardLines, ladderReaderOf, teasers, now}
  function sizedEdgeRows(feedState, context) {
    if (!feedState) return [];
    const selected = feed.selectEdges(feedState, { ...edgeSelectionOptions(context.edgeSettings, context.effective, context.now), minEdge: context.minEdgePct / 100 })
      .map((row) => withStakeAndRank(row, context.measurement, context.stakeSettings));
    return withBetFlags(selected, context);
  }

  // Rows in the list's order, in place: selectEdges already sorts by edge.
  function sortEdgeRows(rows, sortBy) {
    if (sortBy === "stake") rows.sort((a, b) => (b.stake ?? -1) - (a.stake ?? -1) || b.edgePct - a.edgePct);
    if (sortBy === "start") rows.sort((a, b) => a.eventStartMs - b.eventStartMs || b.edgePct - a.edgePct);
    if (sortBy === "exposure") rows.sort((a, b) => exposureDollars(b) - exposureDollars(a) || b.edgePct - a.edgePct);
    return rows;
  }

  // What the Edges list shows (flat): sized rows at the list's minimum edge
  // that meet the Min suggested bet, in the list's sort. Same context as
  // sizedEdgeRows minus minEdgePct, which is the settings' own.
  function listedEdgeRows(feedState, context) {
    const settings = context.edgeSettings;
    const rows = sizedEdgeRows(feedState, { ...context, minEdgePct: settings.minEdgePct })
      .filter((row) => meetsMinStake(row, settings.minStake));
    return sortEdgeRows(rows, settings.sortBy);
  }

  // Cards: the best line of each (game, market, side) is always the highest
  // tail-flex rank score (EV dollars after flex); the list's sort orders the
  // cards through that line.
  function groupEdgeRows(rows, sortBy) {
    const groups = feed.groupEdges(rows, (row) => row.rankScore);
    if (sortBy === "stake") groups.sort((a, b) => (b.best.stake ?? -1) - (a.best.stake ?? -1) || b.best.edgePct - a.best.edgePct);
    if (sortBy === "edge") groups.sort((a, b) => b.best.edgePct - a.best.edgePct || a.eventStartMs - b.eventStartMs);
    if (sortBy === "start") groups.sort((a, b) => a.eventStartMs - b.eventStartMs || b.best.edgePct - a.best.edgePct);
    if (sortBy === "exposure") groups.sort((a, b) => exposureDollars(b.best) - exposureDollars(a.best) || b.best.edgePct - a.best.edgePct);
    return groups;
  }

  // "tail flex: NFL spr 6.2% · CFB 1H tot 10%": the c in use for every
  // spread/total market on the list, in list order; a market too thin to
  // measure shows the fallback. "" for an empty list.
  function describeTailFlex(rows, measurement) {
    if (!rows.length) return "";
    const parts = new Map();
    for (const row of rows) {
      if (row.betTypeId === BET_TYPE_MONEYLINE) continue;
      const key = tailflex.marketKeyOf(row);
      if (parts.has(key)) continue;
      const c = tailflex.cOf(measurement, row);
      const cText = tailflex.isMeasured(measurement, row) ? `${(c * 100).toFixed(1)}%` : `${Math.round(c * 100)}%`;
      const period = row.period === "FG" ? "" : `${row.period} `;
      parts.set(key, `${row.leagueLabel} ${period}${row.betTypeId === BET_TYPE_SPREAD ? "spr" : "tot"} ${cText}`);
    }
    return parts.size ? `tail flex: ${Array.from(parts.values()).join(" · ")}` : "";
  }

  function fmtDollars(value) {
    return value.toLocaleString("en-US", { style: "currency", currency: "USD", minimumFractionDigits: 2, maximumFractionDigits: 2 });
  }

  // The rail under the edge: the number to act on, with the verb on it, then
  // one small line — that it is all the liquidity there is, and what the
  // stake would be with nothing held. "add $250" is not the same instruction
  // as "bet $250" and must not look like it. `row` carries `stake` and `bet`.
  //   {text, note, atSize}  text "add $188.32" / "bet $437.25" / "—"; note
  //                         "all $17 liq · $270.05 alone" or null; atSize when
  //                         the held bets leave nothing to add
  function stakeRail(row) {
    const advice = row.bet ? row.bet.advice : null;
    const words = betsView.stakeAdviceWords(advice);
    if (!words) return { text: row.stake == null ? "—" : `bet ${fmtDollars(row.stake)}`, note: null, atSize: false };
    const noteText = [words.cap, words.alone].filter(Boolean).join(" · ");
    return { text: `${words.verb} ${fmtDollars(advice.bet)}`, note: noteText || null, atSize: advice.bet === 0 };
  }

  // Edge magnitude in three steps, so a +6% and a +1.1% never read the same:
  // the row's left stripe and the figure both take their colour from here.
  function edgeTier(edgePct) {
    if (edgePct >= 4) return "hot";
    if (edgePct >= 2) return "warm";
    return "thin";
  }

  // ---- why the edge grew (#132) ----------------------------------------------

  function fmtAmerican(price) {
    return price > 0 ? `+${price}` : `${price}`;
  }

  // A history entry's fair as a probability, "33.7%"; the American price when
  // it is not one, "?" when there is none.
  function fmtFairEntry(entry) {
    if (entry.bacr == null) return "?";
    try {
      return `${(kelly.americanToProb(entry.bacr) * 100).toFixed(1)}%`;
    } catch (_error) {
      return fmtAmerican(entry.bacr);
    }
  }

  // " +208" for a fill's price, "" for a record with none.
  function fmtBetPrice(bet) {
    return typeof bet.price === "number" ? ` ${fmtAmerican(bet.price)}` : "";
  }

  // The tag shows only on a line the user already holds in the same direction
  // (user decision, 2026-09-23): it exists for the adverse selection of ADDING
  // to a position, and a first bet is not a top-up. stakeAdvice says "add"
  // exactly when held dollars are on the row's direction. The history behind
  // the tag is still recorded for every line (scanner.js), so a line bet later
  // is tagged at once.
  function heldInThisDirection(bet) {
    return Boolean(bet && bet.advice && bet.advice.verb === "add");
  }

  // What the tag reads on a line held in this direction: {move, baseline,
  // recent}. With a saved fill fair on the line (fillfair.baselineOf: the
  // earliest open bet on this very line that has one) the move is since that
  // fill; without one it is the last ten minutes, as before. `recent` is
  // always the ten-minute move. `line` needs key, points and the openers (an
  // Edges row or a raw feed line); `bet` is its {tier, matches, advice}.
  //   context  {history (scanner.getHistory()), fillFairIndex (fillfair.fairsByBetId), now}
  function tagReading(line, bet, context) {
    const recent = edgemove.edgeMove(context.history[line.key], context.now);
    const baseline = fillfair.baselineOf(bet.matches, context.fillFairIndex, line.marketId);
    const entries = context.history[line.key];
    if (!baseline || !entries || entries.length === 0) return { move: recent, baseline: null, recent };
    // An Edges row carries its book as `book`; the Ticket's raw feed line as `bookId`.
    const rowBookId = line.book ? line.book.id : line.bookId;
    return { move: fillfair.moveSinceFill(baseline, entries[entries.length - 1], rowBookId), baseline, recent };
  }

  // The mover, on the card under the tag: the fair then and now when the
  // fair decided (green / red), the price when it was the book (amber),
  // named from your fill when the tag reads since one. Both ends with cents,
  // the last ten minutes and the opener stay in the panel's tooltip.
  function moveDetail(reading) {
    const move = reading.move;
    const fromPrice = typeof move.from.price === "number" ? fmtAmerican(move.from.price) : "?";
    const what = move.kind === "book_away"
      ? `price ${fromPrice} → ${fmtAmerican(move.to.price)}`
      : `fair ${fmtFairEntry(move.from)} → ${fmtFairEntry(move.to)}`;
    return reading.baseline ? `since your${fmtBetPrice(reading.baseline.bet)} bet: ${what}` : what;
  }

  // The tag on a line, or null when the line is not held in this direction
  // or nothing moved: {kind, label, detail, reading} — kind fair_to_you /
  // book_away / fair_against, label its words (edgemove.MOVE_LABELS), detail
  // the mover (moveDetail), reading the tagReading behind it. Same context as tagReading.
  function moveTag(line, bet, context) {
    if (!heldInThisDirection(bet)) return null;
    const reading = tagReading(line, bet, context);
    if (reading.move.kind === "none") return null;
    return { kind: reading.move.kind, label: edgemove.MOVE_LABELS[reading.move.kind], detail: moveDetail(reading), reading };
  }

  // The tag's words for an alert body, or null; the row must carry `bet`
  // (withBetFlags). Same context as tagReading.
  function moveWords(row, context) {
    const tag = moveTag(row, row.bet, context);
    if (!tag) return null;
    return tag.reading.baseline ? `${tag.label} since your${fmtBetPrice(tag.reading.baseline.bet)} bet` : tag.label;
  }

  const api = {
    DEFAULT_STAKE_SETTINGS, ALL_LEAGUE_IDS, DEFAULT_EDGE_SETTINGS, DEFAULT_BOOK_NAMES, SORT_KEYS, MAX_EDGE_ROWS,
    sanitizeStakeSettings, sanitizeEdgeSettings, feedTeamsByLeague, applyBetsPayload,
    liveBooks, defaultBookIds, effectiveFilter, edgeSelectionOptions,
    stakeFor, withStakeAndRank, boardLines, createLadderReaders, withBetFlags, exposureDollars, meetsMinStake,
    sizedEdgeRows, sortEdgeRows, listedEdgeRows, groupEdgeRows, describeTailFlex, edgeTier, fmtDollars, stakeRail,
    fmtAmerican, fmtFairEntry, fmtBetPrice, heldInThisDirection, tagReading, moveDetail, moveTag, moveWords,
  };

  if (inNode) {
    module.exports = api;
  } else {
    root.UnabatedEdgeRows = api;
  }
})(typeof globalThis !== "undefined" ? globalThis : this);
