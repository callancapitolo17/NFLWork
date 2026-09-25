// Since your first fill (the #132 follow-up, user decisions 2026-09-23): the
// fair Unabated showed for a bet's own line at the moment the bet was
// placed, saved once per bet by the bets service, so the edge-move tag on a
// line you are adding to compares now against your first fill instead of
// only the last ten minutes.
//
// Capture (captureFillFairs): the scanner's in-memory history (edgemove.js)
// holds what every line was worth on every snapshot and stream update. A
// new open bet's fill-time fair is the newest observation of its line at or
// before placedAt, read only when the history is a gap-free record of that
// moment: the scanner has watched without a gap since before placedAt
// (scanner.js status.observingSince, and the line's own league
// status.leagueObservingSince) and the line had been seen by then.
// Never a backfill: a bet placed before the panel was watching (open when
// this shipped, placed with the panel hidden or closed) gets nothing and
// keeps the ten-minute tag.
//
// Display (mergeFillFairs, fairsByBetId, baselineOf, moveSinceFill): a line
// held in the same direction is measured from the EARLIEST open straight bet
// on that very line with a saved fair — the fair then against the fair now,
// and on the bet's own book its fill price against the price now — through
// edgemove.classifyMove, the same four cases and the same 0.5-point
// threshold as the ten-minute tag.
//
// Pure: no DOM, no fetch, no chrome.*. Loaded as a plain <script> in
// panel.html after kelly.js, edgemove.js, feed.js and bets.js (exposes
// globalThis.UnabatedFillFair) and via require() in tests/fillfair.test.js.

(function (root) {
  "use strict";

  const inNode = typeof module !== "undefined" && module.exports;
  const kelly = inNode ? require("./kelly.js") : root.UnabatedKelly;
  const feed = inNode ? require("./feed.js") : root.UnabatedFeed;
  const bets = inNode ? require("./bets.js") : root.UnabatedBets;
  const edgemove = inNode ? require("./edgemove.js") : root.UnabatedEdgeMove;

  // Unabated's book id for each bets-service venue it prices (the snapshot's
  // marketSources, 2026-09-23). BFA and Wagerzon are not on the board, so
  // their bets read another book's line at the number.
  const BOOK_ID_OF_VENUE = { kalshi: 105, novig: 89, betonline: 9, prophetx: 66 };
  // An exchange line's price is a probability (feed sourceFormat 4); the bet
  // record keeps its exact fill the same way.
  const EXCHANGE_PROBABILITY_FORMAT = 4;
  const AMERICAN_FORMAT = 1;
  // A placement up to this far past our clock is a venue clock running a
  // little ahead (a fill can reach the panel seconds after it happens): it is
  // looked at again next pass. Further out is a time-zone error, refused —
  // waiting it out would read the history at the wrong moment.
  const CLOCK_SKEW_GRACE_MS = 60 * 1000;
  // The markets the board carries: a bet on anything else (a prop, a team
  // total, Kalshi's F5) can never find a line of its own.
  const BOARD_BET_TYPES = new Set(Object.values(feed.BET_TYPES).map((name) => name.toLowerCase()));
  const BOARD_PERIODS = new Set(Object.values(feed.PERIODS));
  // Why a bet gets no saved fair. Every one is final: the history only ever
  // loses the moment of a fill, it never gains it back.
  const REFUSED = Object.freeze({
    noBoardMarket: "the board lists no market of its bet type and period",
    noPlacedAt: "no placed time on the record",
    future: "placed after now (a clock or time-zone error)",
    notWatching: "placed before the panel was watching",
    tieCaveat: "a Kalshi NO that also wins on a tie has no line of its own",
    leagueGap: "placed while its league's snapshot was not loading",
    notSeen: "its line was first seen after the bet",
    noFair: "no Unabated fair on its line at the fill",
  });

  // bacr is a whole American price (172,172 of 172,172 live NFL/CFB/MLB
  // lines, 2026-09-23), and the service stores it as one.
  function isWholeAmerican(value) {
    return Number.isInteger(value) && kelly.isAmericanPrice(value);
  }

  // The reasons a bet is refused before its line is looked up: they need
  // only the record and the clocks. A NO that also wins on a tie pays on
  // "the other team or a tie", which no line's fair prices.
  function refusalBeforeLookup(bet, placedMs, observingSince, now) {
    if (!BOARD_BET_TYPES.has(bet.betType) || !BOARD_PERIODS.has(bet.period)) return REFUSED.noBoardMarket;
    if (!Number.isFinite(placedMs)) return REFUSED.noPlacedAt;
    if (placedMs > now + CLOCK_SKEW_GRACE_MS) return REFUSED.future;
    if (placedMs < observingSince) return REFUSED.notWatching;
    if (Array.isArray(bet.approx) && bet.approx.includes(bets.TIE_CAVEAT)) return REFUSED.tieCaveat;
    return null;
  }

  // A cheap filter before the matcher: the line is on the bet's period, bet
  // type and number. The matcher (bets.tierOf) then decides the side.
  function couldBeSameLine(bet, line) {
    const betType = feed.BET_TYPES[line.betTypeId];
    return feed.PERIODS[line.periodTypeId] === bet.period
      && typeof betType === "string" && betType.toLowerCase() === bet.betType
      && bets.samePoints(line.points, bet.points);
  }

  // bet id -> the board rows (feed.describeLine) each bet is `same_line` on:
  // its game, period, bet type, side and number, at any book. A bet whose
  // game is not on the board, or is ambiguous there, has none.
  function sameLineRows(betList, state, boardLines) {
    const eventFlags = bets.annotateRows(boardLines, betList, { lines: boardLines });
    const betsByEvent = new Map();
    boardLines.forEach((row, index) => {
      for (const match of eventFlags[index].matches) {
        if (!betsByEvent.has(row.eventId)) betsByEvent.set(row.eventId, []);
        betsByEvent.get(row.eventId).push(match.bet);
      }
    });
    const out = new Map();
    if (!betsByEvent.size) return out;
    const rows = [];
    for (const line of Object.values(state.lines)) {
      // Snapshot lines only, as the fair ladder (ladder.js): the changes
      // stream files an event's team totals under its game total's bet type
      // (feed.js), so a stream-only line at the bet's number can be another
      // market — and a saved fair is permanent.
      if (line.bookId === feed.UNABATED_LINE_BOOK_ID || line.fromSnapshot !== true) continue;
      const onEvent = betsByEvent.get(line.eventId);
      if (onEvent && onEvent.some((bet) => couldBeSameLine(bet, line))) rows.push(feed.describeLine(line, state));
    }
    const lineFlags = bets.annotateRows(rows, betList, { lines: boardLines });
    rows.forEach((row, index) => {
      for (const match of lineFlags[index].matches) {
        if (match.tier !== "same_line") continue;
        if (!out.has(match.bet.id)) out.set(match.bet.id, []);
        out.get(match.bet.id).push(row);
      }
    });
    return out;
  }

  // The bet's own venue first; then the book that changed its line most
  // recently. The fair is the same at every book on a number (NFL 5,796 of
  // 5,796 multi-book rungs, MLB 1,394 of 1,394, CFB 95.3% — the rest mostly
  // alt rungs untouched for hours still carrying an older fair; 2026-09-23).
  function byPreference(rows, bet) {
    const ownBookId = BOOK_ID_OF_VENUE[bet.venue];
    const rank = (row) => (row.book.id === ownBookId ? 0 : 1);
    const changedMs = (row) => row.modifiedMs ?? 0;
    return rows.slice().sort((a, b) => rank(a) - rank(b) || changedMs(b) - changedMs(a));
  }

  // What a line's history says about the moment `placedMs`: the newest
  // observation at or before it, when the line had been seen by then and the
  // fair was a price. {entry} or {reason}.
  function fillTimeEntry(entries, placedMs) {
    if (!Array.isArray(entries) || entries.length === 0 || entries[0].at > placedMs) return { reason: REFUSED.notSeen };
    let entry = entries[0];
    for (const candidate of entries) if (candidate.at <= placedMs) entry = candidate;
    if (!isWholeAmerican(entry.bacr)) return { reason: REFUSED.noFair };
    return { entry };
  }

  // The first line, in preference order, whose league was loading through
  // the fill and whose history holds it: {save} or {reason} (the most
  // preferred line's).
  function decide(bet, placedMs, rows, history, leagueObservingSince) {
    let reason = null;
    for (const row of rows) {
      const leagueSince = leagueObservingSince[row.leagueId];
      if (typeof leagueSince !== "number" || placedMs < leagueSince) {
        reason ||= REFUSED.leagueGap;
        continue;
      }
      const found = fillTimeEntry(history[row.key], placedMs);
      if (found.entry) {
        return { save: {
          betId: bet.id, lineKey: row.key, points: row.points ?? null, fairAmerican: found.entry.bacr,
          fairObservedAt: new Date(found.entry.at).toISOString(), placedAt: bet.placedAt,
        } };
      }
      reason ||= found.reason;
    }
    return { reason };
  }

  // One capture pass: {saves, refusals}. A save is the row POST
  // /fill_fairs.json takes — {betId, lineKey, points, fairAmerican,
  // fairObservedAt, placedAt}, times ISO; a refusal {betId, reason} is final
  // (REFUSED). A bet with no
  // line on the board yet — its game not resolved, no book at its number —
  // is in neither and is looked at again next pass.
  //   records         bet records (bets.js contract); only open ones are read
  //   skipIds         Set of bet ids already saved or decided
  //   state           the scanner's feed state {lines, events, books, teams}
  //   boardLines      one describeLine row per event (the panel's boardLines())
  //   history         edgemove history: line key -> entries
  //   observingSince  ms since when the history has no gap, or null
  //   leagueObservingSince  leagueId -> ms since when that league's snapshots
  //                   have loaded without a gap (scanner status)
  //   now             ms
  function captureFillFairs({ records, skipIds, state, boardLines, history, observingSince, leagueObservingSince, now }) {
    const result = { saves: [], refusals: [] };
    if (!state || typeof observingSince !== "number") return result;
    const pending = [];
    for (const bet of records) {
      if (bet.status !== "open" || bet.unmatchable || skipIds.has(bet.id)) continue;
      const placedMs = Date.parse(bet.placedAt);
      const reason = refusalBeforeLookup(bet, placedMs, observingSince, now);
      if (reason) result.refusals.push({ betId: bet.id, reason });
      else if (placedMs <= now) pending.push({ bet, placedMs });
    }
    if (!pending.length) return result;
    const rowsByBet = sameLineRows(pending.map((item) => item.bet), state, boardLines);
    for (const { bet, placedMs } of pending) {
      const rows = rowsByBet.get(bet.id);
      if (!rows) continue;
      const decision = decide(bet, placedMs, byPreference(rows, bet), history || {}, leagueObservingSince || {});
      if (decision.save) result.saves.push(decision.save);
      else result.refusals.push({ betId: bet.id, reason: decision.reason });
    }
    return result;
  }

  // ---- display ---------------------------------------------------------------

  // The saved fairs to hold after a service reply: the served rows, plus any
  // held row the reply lacks — a /bets.json that left before a POST landed
  // answers without the new row — for the bets still in the records. A saved
  // row never changes, so the union is exact, and the records bound it to
  // the retention window.
  function mergeFillFairs(held, served, records) {
    const byBet = new Map();
    for (const row of [...held, ...served]) if (row && typeof row.betId === "string") byBet.set(row.betId, row);
    const recordIds = new Set(records.map((record) => record.id));
    return Array.from(byBet.values()).filter((row) => recordIds.has(row.betId));
  }

  // bet id -> saved row, over the rows the bets service serves
  // ({betId, fairAmerican, ...}); a malformed row is skipped.
  function fairsByBetId(rows) {
    const index = new Map();
    for (const row of Array.isArray(rows) ? rows : []) {
      if (!row || typeof row.betId !== "string" || !isWholeAmerican(row.fairAmerican)) continue;
      index.set(row.betId, row);
    }
    return index;
  }

  // The exact probability an exchange fill was made at, or null: Kalshi
  // keeps its VWAP in cents, Novig its 0-1 price.
  function exactFillProb(bet) {
    const raw = bet.raw && typeof bet.raw === "object" ? bet.raw : {};
    let prob = null;
    if (bet.venue === "kalshi" && typeof raw.vwapCents === "number") prob = raw.vwapCents / 100;
    if (bet.venue === "novig" && typeof raw.probability === "number") prob = raw.probability;
    return prob != null && prob > 0 && prob < 1 ? prob : null;
  }

  // A bet's fill price as an observation's price fields, on the basis
  // edgemove compares prices: an exchange fill as its exact probability, a
  // sportsbook fill as its American price. Null when the record carries no
  // price (a parlay leg).
  function fillPriceOf(bet) {
    if (!kelly.isAmericanPrice(bet.price)) return null;
    const exact = exactFillProb(bet);
    if (exact == null) return { price: bet.price, sourceFormat: AMERICAN_FORMAT, sourcePrice: null };
    return { price: bet.price, sourceFormat: EXCHANGE_PROBABILITY_FORMAT, sourcePrice: exact };
  }

  // "289357353" from "289357353:ms89:si0:tid6:alt-3.5": Unabated's market for
  // one side of one bet type in one period of one game, shared by every book
  // (1,599 of 1,608 live NFL sides, 2026-09-23).
  function marketIdOfLineKey(lineKey) {
    return typeof lineKey === "string" ? lineKey.split(":")[0] : null;
  }

  // The bet a line held in this direction is measured from: the EARLIEST
  // open straight bet on this very line (same_line — period, bet type, side
  // and number, any venue) with a saved fair, as {bet, fair, placedMs}; null
  // when none has one (the ten-minute tag stays). A parlay leg is not a
  // position on the line and has no fill price of its own. A fair read off
  // another market than the row's is never used: it is permanent, and a bet
  // the matcher later puts on another game must not carry it there.
  //   matches      the row's matchBets matches (open bets only)
  //   fairs        fairsByBetId(...)
  //   rowMarketId  the row's marketId
  function baselineOf(matches, fairs, rowMarketId) {
    let baseline = null;
    for (const match of matches || []) {
      if (match.tier !== "same_line" || !match.bet || match.bet.isParlayLeg === true) continue;
      const fair = fairs.get(match.bet.id);
      const placedMs = Date.parse(match.bet.placedAt);
      if (!fair || !Number.isFinite(placedMs)) continue;
      if (rowMarketId == null || marketIdOfLineKey(fair.lineKey) !== String(rowMarketId)) continue;
      if (!baseline || placedMs < baseline.placedMs) baseline = { bet: match.bet, fair, placedMs };
    }
    return baseline;
  }

  // The four cases from the fill to `current`, the line's newest observation
  // ({price, sourceFormat, sourcePrice, bacr}): the saved fair against the
  // fair now, and the bet's own fill price against the price now — only on
  // the bet's own book. A Novig fill against a Kalshi row is the spread
  // between two books, not a book moving away, so on another book (and for
  // a record with no price) the fair decides alone and `from.price` is null.
  // {kind, fairDelta, priceDelta, from, to}.
  function moveSinceFill(baseline, current, rowBookId) {
    const onOwnBook = rowBookId != null && BOOK_ID_OF_VENUE[baseline.bet.venue] === rowBookId;
    const price = (onOwnBook && fillPriceOf(baseline.bet)) || { price: null, sourceFormat: AMERICAN_FORMAT, sourcePrice: null };
    const from = { ...price, bacr: baseline.fair.fairAmerican };
    return { ...edgemove.classifyMove(from, current), from, to: current };
  }

  const api = { BOOK_ID_OF_VENUE, REFUSED, captureFillFairs, fillTimeEntry, mergeFillFairs, fairsByBetId, fillPriceOf, baselineOf, moveSinceFill };

  if (inNode) {
    module.exports = api;
  } else {
    root.UnabatedFillFair = api;
  }
})(typeof globalThis !== "undefined" ? globalThis : this);
