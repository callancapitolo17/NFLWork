// Parsers for Unabated's two public market feeds. Pure functions: no DOM, no
// fetch, no chrome.* — loaded as a plain <script> in panel.html (exposes
// globalThis.UnabatedFeed) and via require() in tests/feed.test.js.
//
// Feeds (verified 2026-09-10, no login needed):
//   Snapshot  https://content.unabated.com/markets/v2/league/{id}/odds.json
//             {odds: {"lg1:pt1:pregame": [row, ...]}, teams: {id: {name}}, marketSources: [{id, name, ...}]}
//             row.sides["si0:tid6"]["ms4"] = line {marketId, points, americanPrice, price, sourcePrice,
//             sourceFormat, bacr, ge, liquidity, statusId, sequenceNumber}
//   Changes   https://api-k.unabated.com/api/markets/changes/query[/{cursor}]
//             {latestTimestamp, resultCode, results: [{latestTimestamp, marketLineChanges: [{gameOdds:
//             {gameOddsEvents: {"lg1:pt1:pregame": [{eventId, eventStart, gameOddsMarketSourcesLines:
//             {"si0:ms4:an0": {"bt2": line}}}]}}}]}]}
//
// One line is identified by (marketId, book, sideKey): a marketId is one
// (event, period, bet type, side) market and the same key can carry several
// rows' worth of updates (team totals reuse bt3), so event+betType is NOT a
// key. Updates apply only when their sequenceNumber is newer than the line we
// hold — the stream replays old lines and the snapshot may already be ahead.

(function (root) {
  "use strict";

  const LEAGUES = {
    1: { label: "NFL", path: "nfl" },
    2: { label: "CFB", path: "cfb" },
    5: { label: "MLB", path: "mlb" },
  };
  const BET_TYPES = { 1: "Moneyline", 2: "Spread", 3: "Total" };
  const PERIODS = { 1: "FG", 2: "1H", 3: "2H", 4: "1Q", 5: "2Q", 6: "3Q", 7: "4Q" };
  // ms49 is Unabated's own line, not a book anyone can bet.
  const UNABATED_LINE_BOOK_ID = 49;
  const STATUS_ON_BOARD = 1;
  // Cursor = nanoseconds since 2021-01-06T00:00:00Z (decoded from the feed's own
  // latestTimestamp vs modifiedOn pairs, 2026-09-10). Kept as a STRING: it is
  // above 2^53 and JSON.parse would round it.
  const CURSOR_EPOCH_MS = Date.UTC(2021, 0, 6);
  const GAME_ROW_KEY = /^pt(\d+):pregame:bt([123]):e(\d+)$/;
  const LEAGUE_KEY = /^lg(\d+):pt(\d+):(pregame|live)$/;
  const LINE_KEY_RE = /^si(\d):ms(\d+):an(\d+)$/;

  // ---- small helpers -------------------------------------------------------

  function numberOrNull(value) {
    return typeof value === "number" && Number.isFinite(value) ? value : null;
  }

  function parseLeagueKey(key) {
    const match = LEAGUE_KEY.exec(key);
    if (!match) return null;
    return { leagueId: Number(match[1]), periodTypeId: Number(match[2]), phase: match[3] };
  }

  // "2026-09-13T17:00:00" (snapshot, naive UTC) or "2026-09-13T17:00:00+00:00" (changes).
  function parseEventStart(value) {
    if (typeof value !== "string" || !value) return null;
    const hasZone = /(?:Z|[+-]\d\d:\d\d)$/.test(value);
    const ms = Date.parse(hasZone ? value : `${value}Z`);
    return Number.isFinite(ms) ? ms : null;
  }

  function sideIndexOf(sideKey) {
    const match = /^si(\d)/.exec(sideKey || "");
    return match ? Number(match[1]) : null;
  }

  function lineKeyOf({ marketId, bookId, sideKey }) {
    return `${marketId}:ms${bookId}:${sideKey}`;
  }

  // Sequence numbers are monotonic per line (the same feed cannot go backwards).
  function isNewer(candidateSeq, heldSeq) {
    if (typeof candidateSeq !== "number") return false;
    if (typeof heldSeq !== "number") return true;
    return candidateSeq > heldSeq;
  }

  function emptyState() {
    return { leagues: [], teams: {}, books: {}, events: {}, lines: {} };
  }

  // ---- snapshot ------------------------------------------------------------

  function normalizeSnapshotLine(raw, context) {
    const price = numberOrNull(raw.americanPrice) ?? numberOrNull(raw.price);
    if (price == null || raw.marketId == null) return null;
    return {
      key: lineKeyOf({ marketId: raw.marketId, bookId: context.bookId, sideKey: context.sideKey }),
      leagueId: context.leagueId,
      periodTypeId: context.periodTypeId,
      betTypeId: context.betTypeId,
      eventId: context.eventId,
      marketId: raw.marketId,
      bookId: context.bookId,
      sideKey: context.sideKey,
      sideIndex: sideIndexOf(context.sideKey),
      points: numberOrNull(raw.points),
      price,
      sourceFormat: numberOrNull(raw.sourceFormat) ?? 1,
      sourcePrice: numberOrNull(raw.sourcePrice),
      bacr: numberOrNull(raw.bacr),
      ge: numberOrNull(raw.ge),
      liquidity: numberOrNull(raw.liquidity),
      statusId: numberOrNull(raw.statusId),
      sequenceNumber: numberOrNull(raw.sequenceNumber),
      isBlurred: raw.isBlurred === true,
      modifiedOn: raw.modifiedOn ?? null,
    };
  }

  function ingestSnapshotRow(state, row, leagueId, periodTypeId, counts) {
    const rowKey = typeof row.key === "string" ? GAME_ROW_KEY.exec(row.key) : null;
    // Only pregame moneyline/spread/total game rows; team totals (bt4), props and
    // live rows share the structure but are not what the scanner lists.
    if (!rowKey || row.live === true) {
      counts.skippedRows += 1;
      return;
    }
    const betTypeId = Number(rowKey[2]);
    const eventId = row.eventId;
    if (eventId == null) {
      counts.skippedRows += 1;
      return;
    }
    if (!state.events[eventId]) {
      const teams = row.eventTeams || {};
      state.events[eventId] = {
        eventId,
        leagueId,
        eventName: row.eventName ?? null,
        eventStart: parseEventStart(row.eventStart),
        awayTeamId: teams[0] ? teams[0].id ?? null : null,
        homeTeamId: teams[1] ? teams[1].id ?? null : null,
        awayRotation: teams[0] ? teams[0].rotationNumber ?? null : null,
        homeRotation: teams[1] ? teams[1].rotationNumber ?? null : null,
      };
    }
    counts.rows += 1;
    for (const [sideKey, books] of Object.entries(row.sides || {})) {
      if (!books) continue;
      for (const [bookKey, raw] of Object.entries(books)) {
        const bookId = Number(bookKey.replace(/^ms/, ""));
        if (!Number.isInteger(bookId) || !raw) continue;
        const line = normalizeSnapshotLine(raw, { leagueId, periodTypeId, betTypeId, eventId, bookId, sideKey });
        if (!line) {
          counts.skippedLines += 1;
          continue;
        }
        state.lines[line.key] = line;
        counts.lines += 1;
      }
    }
  }

  // Parse one league's odds.json. Returns a state for that league only; merge
  // several with mergeStates. Books carry the flags the scanner filters on.
  function parseSnapshot(json, { leagueId }) {
    if (!json || typeof json !== "object" || !json.odds || typeof json.odds !== "object") {
      throw new Error("snapshot: expected an object with an `odds` map");
    }
    const state = emptyState();
    state.leagues = [leagueId];
    const counts = { rows: 0, skippedRows: 0, lines: 0, skippedLines: 0 };
    for (const [teamId, team] of Object.entries(json.teams || {})) {
      if (team && team.name) state.teams[teamId] = team.name;
    }
    for (const source of Array.isArray(json.marketSources) ? json.marketSources : []) {
      if (!source || source.id == null) continue;
      state.books[source.id] = {
        id: source.id,
        name: source.name || `book ${source.id}`,
        // What the odds screen itself shows: an active source whose game-odds status is live.
        isLive: source.isActive === true && source.statusId === STATUS_ON_BOARD,
        hasLiquidity: source.hasLiquidity === true,
      };
    }
    for (const [leagueKey, rows] of Object.entries(json.odds)) {
      const parsed = parseLeagueKey(leagueKey);
      if (!parsed || parsed.leagueId !== leagueId || parsed.phase !== "pregame" || !Array.isArray(rows)) continue;
      for (const row of rows) ingestSnapshotRow(state, row, leagueId, parsed.periodTypeId, counts);
    }
    if (counts.lines === 0) {
      throw new Error(`snapshot: no game lines found for league ${leagueId} (rows seen: ${counts.rows + counts.skippedRows})`);
    }
    state.counts = counts;
    return state;
  }

  function mergeStates(states) {
    const merged = emptyState();
    for (const state of states) {
      merged.leagues.push(...state.leagues);
      Object.assign(merged.teams, state.teams);
      Object.assign(merged.books, state.books);
      Object.assign(merged.events, state.events);
      Object.assign(merged.lines, state.lines);
    }
    return merged;
  }

  // ---- changes -------------------------------------------------------------

  // The top-level latestTimestamp is the next cursor. Read it from the raw text
  // so it stays exact (it does not fit in a double).
  function extractCursor(text) {
    const match = /"latestTimestamp"\s*:\s*(\d+)/.exec(text);
    return match ? match[1] : null;
  }

  function cursorFromDate(date) {
    const ms = date instanceof Date ? date.getTime() : Number(date);
    if (!Number.isFinite(ms) || ms < CURSOR_EPOCH_MS) return null;
    // Whole seconds only: the ms part would need BigInt to stay exact and the
    // server is happy with a cursor up to ~3 minutes old.
    return `${Math.floor((ms - CURSOR_EPOCH_MS) / 1000)}000000000`;
  }

  function normalizeChangeLine(raw, context) {
    const price = numberOrNull(raw.price);
    if (price == null || raw.marketId == null || typeof raw.sideKey !== "string") return null;
    return {
      key: lineKeyOf({ marketId: raw.marketId, bookId: context.bookId, sideKey: raw.sideKey }),
      leagueId: context.leagueId,
      periodTypeId: context.periodTypeId,
      betTypeId: context.betTypeId,
      eventId: context.eventId,
      eventStart: context.eventStart,
      marketId: raw.marketId,
      bookId: context.bookId,
      sideKey: raw.sideKey,
      sideIndex: sideIndexOf(raw.sideKey),
      points: numberOrNull(raw.points),
      price,
      sourceFormat: numberOrNull(raw.sourceFormat) ?? 1,
      sourcePrice: numberOrNull(raw.sourcePrice),
      bacr: numberOrNull(raw.bacr),
      ge: numberOrNull(raw.ge),
      statusId: numberOrNull(raw.statusId),
      sequenceNumber: numberOrNull(raw.sequenceNumber),
      isBlurred: raw.isBlurred === true,
      modifiedOn: raw.modifiedOn ?? null,
    };
  }

  // Parse one changes response (object or raw text). `ok` is false when the
  // server rejected the cursor (resultCode "Failed": too old, > ~3 min) — the
  // caller must resync from a snapshot.
  function parseChanges(input) {
    const text = typeof input === "string" ? input : null;
    const json = text ? JSON.parse(text) : input;
    if (!json || typeof json !== "object" || !Array.isArray(json.results)) {
      throw new Error("changes: expected an object with a `results` array");
    }
    const cursor = text ? extractCursor(text) : (json.latestTimestamp != null ? String(json.latestTimestamp) : null);
    const lines = [];
    let batches = 0;
    for (const result of json.results) {
      if (!result || !Array.isArray(result.marketLineChanges)) continue;
      batches += 1;
      for (const change of result.marketLineChanges) {
        const events = change && change.gameOdds && change.gameOdds.gameOddsEvents;
        if (!events || typeof events !== "object") continue;
        for (const [leagueKey, eventList] of Object.entries(events)) {
          const parsed = parseLeagueKey(leagueKey);
          if (!parsed || parsed.phase !== "pregame" || !Array.isArray(eventList)) continue;
          for (const event of eventList) {
            if (!event || event.eventId == null) continue;
            const eventStart = parseEventStart(event.eventStart);
            for (const [lineKey, byBetType] of Object.entries(event.gameOddsMarketSourcesLines || {})) {
              const keyMatch = LINE_KEY_RE.exec(lineKey);
              if (!keyMatch || !byBetType) continue;
              const bookId = Number(keyMatch[2]);
              for (const [betTypeKey, raw] of Object.entries(byBetType)) {
                const betTypeId = Number(betTypeKey.replace(/^bt/, ""));
                if (!BET_TYPES[betTypeId] || !raw) continue;
                const line = normalizeChangeLine(raw, {
                  leagueId: parsed.leagueId, periodTypeId: parsed.periodTypeId, betTypeId,
                  eventId: event.eventId, eventStart, bookId,
                });
                if (line) lines.push(line);
              }
            }
          }
        }
      }
    }
    return { ok: json.resultCode === "Success", resultCode: json.resultCode ?? null, cursor, batches, lines };
  }

  // Apply parsed changes to a merged state in place. Lines for events the
  // snapshot never listed are skipped (we have no teams/name for them); lines
  // for known events are added when new and replaced when their sequence
  // number is newer. Snapshot-only fields (liquidity) carry over on replace.
  function applyChanges(state, changes) {
    const counts = { applied: 0, added: 0, stale: 0, unknownEvent: 0, otherLeague: 0 };
    const leagues = new Set(state.leagues);
    for (const line of changes.lines) {
      if (!leagues.has(line.leagueId)) {
        counts.otherLeague += 1;
        continue;
      }
      const event = state.events[line.eventId];
      if (!event) {
        counts.unknownEvent += 1;
        continue;
      }
      if (line.eventStart != null && event.eventStart !== line.eventStart) event.eventStart = line.eventStart;
      const held = state.lines[line.key];
      if (held && !isNewer(line.sequenceNumber, held.sequenceNumber)) {
        counts.stale += 1;
        continue;
      }
      const { eventStart, ...fields } = line;
      state.lines[line.key] = { liquidity: held ? held.liquidity : null, ...fields };
      counts.applied += 1;
      if (!held) counts.added += 1;
    }
    return counts;
  }

  // ---- edge selection ------------------------------------------------------

  function signedPoints(points) {
    return points > 0 ? `+${points}` : `${points}`;
  }

  function sideLabelOf(line, event, teams) {
    if (line.betTypeId === 3) {
      return `${line.sideIndex === 0 ? "Over" : "Under"} ${line.points}`;
    }
    const teamId = line.sideIndex === 0 ? event.awayTeamId : event.homeTeamId;
    const team = teams[teamId] || teams[String(teamId)] || `team ${teamId}`;
    if (line.betTypeId === 2) return `${team} ${signedPoints(line.points)}`;
    return team;
  }

  // Ticket-shaped so the panel can describe an edge row with the same wording
  // as a captured ticket (sideLabel, betType, sideIndex, homeTeam, awayTeam...).
  function describeLine(line, state) {
    const event = state.events[line.eventId];
    const league = LEAGUES[line.leagueId] || { label: `league ${line.leagueId}`, path: String(line.leagueId) };
    const book = state.books[line.bookId] || { id: line.bookId, name: `book ${line.bookId}`, isLive: false, hasLiquidity: false };
    return {
      key: line.key,
      leagueId: line.leagueId,
      league: league.path,
      leagueLabel: league.label,
      eventId: line.eventId,
      eventName: event.eventName,
      eventStart: event.eventStart == null ? null : new Date(event.eventStart).toISOString(),
      eventStartMs: event.eventStart,
      betTypeId: line.betTypeId,
      betType: BET_TYPES[line.betTypeId],
      periodTypeId: line.periodTypeId,
      period: PERIODS[line.periodTypeId] || `pt${line.periodTypeId}`,
      sideIndex: line.sideIndex,
      sideKey: line.sideKey,
      sideLabel: sideLabelOf(line, event, state.teams),
      awayTeam: state.teams[event.awayTeamId] ?? null,
      homeTeam: state.teams[event.homeTeamId] ?? null,
      homeAway: line.sideIndex === 0 ? "Away" : "Home",
      rotation: line.sideIndex === 0 ? event.awayRotation : event.homeRotation,
      points: line.points,
      book: { id: book.id, name: book.name, hasLiquidity: book.hasLiquidity },
      price: line.price,
      sourceFormat: line.sourceFormat,
      sourcePrice: line.sourcePrice,
      fair: line.bacr,
      // ge has 4 decimals; round so 0.0156 * 100 prints as 1.56, not 1.5599999.
      edgePct: line.ge == null ? null : Math.round(line.ge * 1e6) / 1e4,
      liquidity: line.liquidity,
      marketId: line.marketId,
      isBlurred: line.isBlurred,
    };
  }

  // Lines worth listing: on the board, edge known and >= minEdge (a fraction),
  // period/bet type enabled, book allowed, game not started. Sorted by edge.
  function selectEdges(state, options) {
    const opts = options || {};
    const minEdge = typeof opts.minEdge === "number" ? opts.minEdge : 0.01;
    const periods = opts.periods instanceof Set ? opts.periods : new Set([1]);
    const betTypes = opts.betTypes instanceof Set ? opts.betTypes : new Set([1, 2, 3]);
    const bookIds = opts.bookIds instanceof Set ? opts.bookIds : null;
    const now = typeof opts.now === "number" ? opts.now : Date.now();
    const rows = [];
    for (const line of Object.values(state.lines)) {
      if (line.bookId === UNABATED_LINE_BOOK_ID) continue;
      if (line.statusId !== STATUS_ON_BOARD) continue;
      if (line.ge == null || line.ge < minEdge) continue;
      if (!periods.has(line.periodTypeId) || !betTypes.has(line.betTypeId)) continue;
      const book = state.books[line.bookId];
      if (bookIds ? !bookIds.has(line.bookId) : !(book && book.isLive)) continue;
      const event = state.events[line.eventId];
      if (!event || event.eventStart == null || event.eventStart <= now) continue;
      if (line.betTypeId !== 1 && line.points == null) continue;
      rows.push(describeLine(line, state));
    }
    rows.sort((a, b) => b.edgePct - a.edgePct || a.eventStartMs - b.eventStartMs);
    return rows;
  }

  function countLines(state) {
    return Object.keys(state.lines).length;
  }

  const api = {
    LEAGUES, BET_TYPES, PERIODS, UNABATED_LINE_BOOK_ID, CURSOR_EPOCH_MS,
    parseLeagueKey, parseEventStart, lineKeyOf, emptyState,
    parseSnapshot, mergeStates, extractCursor, cursorFromDate, parseChanges, applyChanges,
    describeLine, selectEdges, countLines,
  };

  if (typeof module !== "undefined" && module.exports) {
    module.exports = api;
  } else {
    root.UnabatedFeed = api;
  }
})(typeof globalThis !== "undefined" ? globalThis : this);
