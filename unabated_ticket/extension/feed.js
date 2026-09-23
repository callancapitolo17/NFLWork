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
// One line is identified by (marketId, book, sideKey). Snapshot game rows are
// unique per (event, period, bet type), but the changes stream tags OTHER
// markets of the same event with the same bt key (e.g. team totals under
// bt3) — only the marketId tells them apart, so event+betType is NOT a key
// for updates. Updates apply only when their sequenceNumber is newer than the
// line we hold — the stream replays old lines and the snapshot may already be
// ahead.
//
// Alternate lines (issue #113): each snapshot line carries alternateLines[]
// with the same price/edge fields at other points. They are expanded into
// lines keyed (marketId, book, sideKey, points) — `<mainKey>:alt<points>` —
// flagged isAlt with mainPoints = the book's own main-line points. Facts
// measured on the live NFL file 2026-09-11 that shape the parsing:
//   - the changes stream carries NO alt updates (2,603 keys, all an0, no
//     alternateLines), so alts refresh only with the per-league snapshot;
//   - every alt's modifiedOn is the sentinel "0001-01-01T00:00:00", but its
//     sequenceNumber is the change time in epoch ms (on 10,182 main lines it
//     trailed modifiedOn by a median 1.2 s), so alt freshness reads from it;
//   - an alt's marketId can be null (Fanatics) and its marketSourceId can
//     name another book (Sports Interaction mirrors BetMGM's id 4), so the
//     key uses the parent line's marketId and the ms<id> the alt sits under;
//   - `stn` is the market's standard number, not this book's main points
//     (Hard Rock: stn 47.5 on a 48.0 main), so mainPoints comes from the
//     parent line; and alternateLines can hold null entries.

(function (root) {
  "use strict";

  // Every team-sport league the v2 feed served on 2026-09-10 (ids probed
  // 1-70; labels read off the fixtures' team names). `path` is the
  // tools.unabated.com odds screen for the row click: nfl/cfb/mlb are
  // verified, the rest follow the site's nav (nba, cbb, nhl, wnba, soccer)
  // and are unverified without a login. Tennis (9, 10) and combat (22) key
  // sides on people and use other bet types, so they are not listed.
  const LEAGUES = {
    1: { label: "NFL", path: "nfl", sport: "football" },
    2: { label: "CFB", path: "cfb", sport: "football" },
    3: { label: "NBA", path: "nba", sport: "basketball" },
    4: { label: "CBB", path: "cbb", sport: "basketball" },
    7: { label: "WNBA", path: "wnba", sport: "basketball" },
    5: { label: "MLB", path: "mlb", sport: "baseball" },
    12: { label: "WBC", path: "mlb", sport: "baseball" },
    6: { label: "NHL", path: "nhl", sport: "hockey" },
    11: { label: "Olympic hockey", path: "nhl", sport: "hockey" },
    21: { label: "Intl soccer", path: "soccer", sport: "soccer" },
    25: { label: "MLS", path: "soccer", sport: "soccer" },
    26: { label: "La Liga", path: "soccer", sport: "soccer" },
    27: { label: "Serie A", path: "soccer", sport: "soccer" },
    28: { label: "Premier League", path: "soccer", sport: "soccer" },
    29: { label: "Europa League", path: "soccer", sport: "soccer" },
    30: { label: "Bundesliga", path: "soccer", sport: "soccer" },
    31: { label: "Ligue 1", path: "soccer", sport: "soccer" },
    32: { label: "Liga MX", path: "soccer", sport: "soccer" },
    33: { label: "Mexico (lg33)", path: "soccer", sport: "soccer" },
    34: { label: "Primeira Liga", path: "soccer", sport: "soccer" },
    35: { label: "Belgian Pro League", path: "soccer", sport: "soccer" },
    36: { label: "Eredivisie", path: "soccer", sport: "soccer" },
    37: { label: "EFL Championship", path: "soccer", sport: "soccer" },
    38: { label: "Serie B", path: "soccer", sport: "soccer" },
    39: { label: "Scottish Premiership", path: "soccer", sport: "soccer" },
    41: { label: "English cups", path: "soccer", sport: "soccer" },
    42: { label: "German cup", path: "soccer", sport: "soccer" },
    43: { label: "Danish Superliga", path: "soccer", sport: "soccer" },
    44: { label: "Swiss Super League", path: "soccer", sport: "soccer" },
  };
  const SPORTS = {
    football: "Football", basketball: "Basketball", baseball: "Baseball", hockey: "Hockey", soccer: "Soccer",
  };

  function leagueIdsOfSport(sport) {
    return Object.entries(LEAGUES).filter(([, league]) => league.sport === sport).map(([id]) => Number(id));
  }
  const BET_TYPES = { 1: "Moneyline", 2: "Spread", 3: "Total" };
  const PERIODS = { 1: "FG", 2: "1H", 3: "2H", 4: "1Q", 5: "2Q", 6: "3Q", 7: "4Q" };
  // ms49 is Unabated's own line, not a book anyone can bet.
  const UNABATED_LINE_BOOK_ID = 49;
  const STATUS_ON_BOARD = 1;
  // The feed writes this in place of an unknown modifiedOn (every alt line).
  const MODIFIED_ON_UNKNOWN_PREFIX = "0001-";
  // A sequenceNumber below this cannot be an epoch-ms change time (2021-01-06,
  // the changes cursor epoch); anything older is a counter, not a clock.
  const SEQUENCE_AS_EPOCH_MS_MIN = Date.UTC(2021, 0, 6);
  // Cursor = nanoseconds since 2021-01-06T00:00:00Z (decoded from the feed's own
  // latestTimestamp vs modifiedOn pairs, 2026-09-10). Kept as a STRING: it is
  // above 2^53 and JSON.parse would round it.
  const CURSOR_EPOCH_MS = Date.UTC(2021, 0, 6);
  const GAME_ROW_KEY = /^pt(\d+):pregame:bt([123]):e(\d+)$/;
  const LEAGUE_KEY = /^lg(\d+):pt(\d+):(pregame|live)$/;
  const LINE_KEY_RE = /^si(\d):ms(\d+):an(\d+)$/;
  // Venue ids on alternate-line rungs (#118 step 2; measured on the live NFL
  // and CFB files 2026-09-15). Every Kalshi rung's sourceKey is the contract
  // side + full market ticker, "Y-KXNCAAFSPREAD-26SEP19DUQWSU-WSU36" (the
  // other side of the same contract reads "N-…"; 12,318 of 12,318 matched the
  // pattern below); every Novig rung's sourceData is the Novig outcome id (a
  // UUID). Main lines never carry either field and the changes stream has
  // neither, so venue ids refresh with the snapshot only. The fields are
  // undocumented: any other shape is "no id", never an error.
  const KALSHI_BOOK_ID = 105;
  const NOVIG_BOOK_ID = 89;
  const KALSHI_CONTRACT_RE = /^([YN])-([A-Z0-9]+)-([A-Z0-9]+)-([A-Z0-9.]+)$/;
  const NOVIG_OUTCOME_ID_RE = /^[0-9a-f]{8}-[0-9a-f]{4}-[0-9a-f]{4}-[0-9a-f]{4}-[0-9a-f]{12}$/;

  // ---- small helpers -------------------------------------------------------

  function numberOrNull(value) {
    return typeof value === "number" && Number.isFinite(value) ? value : null;
  }

  function stringOrNull(value) {
    return typeof value === "string" && value !== "" ? value : null;
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

  function altLineKeyOf({ marketId, bookId, sideKey, points }) {
    return `${lineKeyOf({ marketId, bookId, sideKey })}:alt${points}`;
  }

  // The feed's modifiedOn, or null when it is missing or the sentinel.
  function parseModifiedOn(value) {
    if (typeof value !== "string" || value.startsWith(MODIFIED_ON_UNKNOWN_PREFIX)) return null;
    return parseEventStart(value);
  }

  // When the book last changed the line, in epoch ms. Main lines read their
  // modifiedOn; an alt has only the sentinel, so it reads its sequenceNumber
  // (epoch ms of the change, see the header). Null when unknowable.
  function lineChangedMs(line) {
    const fromModifiedOn = parseModifiedOn(line.modifiedOn);
    if (fromModifiedOn != null) return fromModifiedOn;
    if (!line.isAlt) return null;
    const sequence = line.sequenceNumber;
    return typeof sequence === "number" && sequence >= SEQUENCE_AS_EPOCH_MS_MIN ? sequence : null;
  }

  // Sequence numbers are monotonic per line (the same feed cannot go backwards).
  function isNewer(candidateSeq, heldSeq) {
    if (typeof candidateSeq !== "number") return false;
    if (typeof heldSeq !== "number") return true;
    return candidateSeq > heldSeq;
  }

  function emptyState() {
    return { leagues: [], teams: {}, teamIndex: {}, books: {}, events: {}, lines: {} };
  }

  // A row's eventName spells both teams in a form the `teams` map does not
  // carry (#118, measured 2026-09-12): CFB "Prairie View A&M Panthers - PRV @
  // Baylor Bears - BAY" (the map says "Prairie View"), CBB the same with the
  // abbreviation empty ("UConn Huskies - @ Michigan Wolverines -"), NFL / NBA
  // / NHL "Ravens Baltimore BAL @ Cowboys Dallas DAL", MLB "Dodgers Los Angeles
  // @ Marlins Miami", MLS the plain club name. Returns [away, home] spellings
  // (side 0 of eventTeams is away), or null when the name is not "a @ b".
  // Per side: drop a trailing " - <abbr>" (or a bare " -"), else a trailing
  // token equal to that team's own abbreviation, else keep it whole.
  function teamSpellingsFromEventName(eventName, abbreviations) {
    if (typeof eventName !== "string") return null;
    const sides = eventName.split(" @ ");
    if (sides.length !== 2) return null;
    return sides.map((side, index) => {
      const withoutSuffix = side.replace(/ -(?: \S+)?$/, "").trim();
      if (withoutSuffix !== side.trim()) return withoutSuffix;
      const abbreviation = abbreviations[index];
      const trailing = abbreviation ? ` ${abbreviation}` : null;
      if (trailing && withoutSuffix.endsWith(trailing) && withoutSuffix.length > trailing.length) {
        return withoutSuffix.slice(0, -trailing.length).trim();
      }
      return withoutSuffix;
    }).map((spelling) => (spelling === "" ? null : spelling));
  }

  // Register the eventName spellings on the event's teams in the league's
  // team index, as a second name for the same team id.
  function noteEventNameSpellings(state, row) {
    const teams = row.eventTeams || {};
    const ids = [teams[0] ? teams[0].id ?? null : null, teams[1] ? teams[1].id ?? null : null];
    const entries = ids.map((id) => (id == null ? null : state.teamIndex[String(id)] ?? null));
    const spellings = teamSpellingsFromEventName(row.eventName, entries.map((entry) => (entry ? entry.abbreviation : null)));
    if (!spellings) return;
    entries.forEach((entry, index) => {
      if (entry && spellings[index] && entry.eventName == null) entry.eventName = spellings[index];
    });
  }

  // ---- venue ids -----------------------------------------------------------

  // The Kalshi event suffix of a rung's sourceKey ("26SEP19DUQWSU"), or null.
  // Kept whole: Kalshi's team codes are not Unabated's abbreviations (CFB
  // `ELONURI` vs `ELO`+`RIL`), so the suffix is only ever compared as a string.
  function kalshiEventSuffixOf(sourceKey) {
    const match = typeof sourceKey === "string" ? KALSHI_CONTRACT_RE.exec(sourceKey) : null;
    return match ? match[3] : null;
  }

  function isNovigOutcomeId(sourceData) {
    return typeof sourceData === "string" && NOVIG_OUTCOME_ID_RE.test(sourceData);
  }

  function emptyVenueIds() {
    return { kalshiEventSuffixes: [], kalshiContracts: {}, novigOutcomes: {} };
  }

  // Record one rung's venue id on its event: id -> {lineKey, mainKey, points,
  // sideIndex}. `points` and `sideIndex` are the contract's own strike and
  // Unabated side — fixed for the id, so the bet matcher reads them to tell
  // which side a bet whose team names do not resolve is on (#118 step 3).
  // `lineKey` is the listed line at the rung's number — the alt, or the main
  // line when the rung sits on the main number (normalizeAltLine drops that
  // rung, but its id is exactly what a main-line bet joins on) — or null when
  // the rung is not listed (unpriced); the id still names the event and
  // market, so it is kept. `alt` is normalizeAltLine's result for the rung.
  // Applies to Kalshi and Novig alike. The map is rebuilt only by the next
  // snapshot, while the changes stream can move a main line off `points`, so
  // a joiner must check `points` against the line's current points.
  function noteRungVenueIds(venueIds, rung, mainLine, alt) {
    if (!rung || typeof rung !== "object") return;
    const points = numberOrNull(rung.points);
    const onMainNumber = points != null && points === mainLine.points;
    const lineKey = alt ? alt.key : onMainNumber ? mainLine.key : null;
    const target = { lineKey, mainKey: mainLine.key, points, sideIndex: mainLine.sideIndex };
    if (mainLine.bookId === KALSHI_BOOK_ID) {
      const suffix = kalshiEventSuffixOf(rung.sourceKey);
      if (!suffix) return;
      if (!venueIds.kalshiEventSuffixes.includes(suffix)) venueIds.kalshiEventSuffixes.push(suffix);
      venueIds.kalshiContracts[rung.sourceKey] = target;
    } else if (mainLine.bookId === NOVIG_BOOK_ID && isNovigOutcomeId(rung.sourceData)) {
      venueIds.novigOutcomes[rung.sourceData] = target;
    }
  }

  // ---- snapshot ------------------------------------------------------------

  function normalizeSnapshotLine(raw, context) {
    const price = numberOrNull(raw.americanPrice) ?? numberOrNull(raw.price);
    if (price == null || raw.marketId == null) return null;
    return {
      key: lineKeyOf({ marketId: raw.marketId, bookId: context.bookId, sideKey: context.sideKey }),
      isAlt: false,
      // A snapshot game row is the game's own market. The changes stream is
      // not: it files team totals under the game total's bt3 (see the header),
      // so only snapshot lines may feed a fair ladder (ladder.js, #130).
      fromSnapshot: true,
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

  // One alternateLines[] entry, keyed under its parent main line. Null when
  // the entry is null, unpriced, has no points, hangs off a main line with no
  // points (nothing to measure distance from), or sits at the main line's own
  // points (that would be the same bet listed twice). Moneylines have no alts.
  function normalizeAltLine(raw, mainLine) {
    if (!raw || typeof raw !== "object" || mainLine.points == null) return null;
    const points = numberOrNull(raw.points);
    const price = numberOrNull(raw.americanPrice) ?? numberOrNull(raw.price);
    if (points == null || price == null || points === mainLine.points) return null;
    return {
      key: altLineKeyOf({ marketId: mainLine.marketId, bookId: mainLine.bookId, sideKey: mainLine.sideKey, points }),
      isAlt: true,
      fromSnapshot: true,
      mainKey: mainLine.key,
      mainPoints: mainLine.points,
      leagueId: mainLine.leagueId,
      periodTypeId: mainLine.periodTypeId,
      betTypeId: mainLine.betTypeId,
      eventId: mainLine.eventId,
      marketId: mainLine.marketId,
      bookId: mainLine.bookId,
      sideKey: mainLine.sideKey,
      sideIndex: mainLine.sideIndex,
      points,
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
      // The book's own id for the rung (see KALSHI_CONTRACT_RE): Kalshi fills
      // sourceKey, Novig and ProphetX fill sourceData, many books both.
      sourceData: stringOrNull(raw.sourceData),
      sourceKey: stringOrNull(raw.sourceKey),
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
        venueIds: emptyVenueIds(),
      };
      noteEventNameSpellings(state, row);
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
        // Spread/total alts only: a moneyline has no points to be an alt of.
        if (betTypeId === 1 || !Array.isArray(raw.alternateLines)) continue;
        for (const rawAlt of raw.alternateLines) {
          const alt = normalizeAltLine(rawAlt, line);
          noteRungVenueIds(state.events[eventId].venueIds, rawAlt, line, alt);
          if (!alt) {
            counts.skippedAltLines += 1;
            continue;
          }
          state.lines[alt.key] = alt;
          counts.altLines += 1;
        }
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
    const counts = { rows: 0, skippedRows: 0, lines: 0, skippedLines: 0, altLines: 0, skippedAltLines: 0 };
    for (const [teamId, team] of Object.entries(json.teams || {})) {
      if (!team || !team.name) continue;
      state.teams[teamId] = team.name;
      // The full entry, for the team index teams.js builds at runtime (#116).
      // eventName is filled from the team's game row (noteEventNameSpellings); null until one is seen.
      state.teamIndex[teamId] = { id: team.id ?? Number(teamId), name: team.name, abbreviation: team.abbreviation ?? null, leagueId: team.leagueId ?? leagueId, eventName: null };
    }
    for (const source of Array.isArray(json.marketSources) ? json.marketSources : []) {
      if (!source || source.id == null) continue;
      state.books[source.id] = {
        id: source.id,
        name: source.name || `book ${source.id}`,
        // An active source enabled for game odds. marketSources.statusId is
        // NOT a liveness flag: on 2026-09-11 Caesars, Bet365, Fliff, Bet105,
        // BetOnline and Underdog Prediction Market carried statusId 2 or 3
        // with lines changed minutes earlier, and the old `statusId == 1`
        // rule hid them from the Books filter. Dead feeds (Matchbook, the
        // pool books) are isActive false or caught by the max-line-age gate.
        isLive: source.isActive === true && source.isEnabledForGameOdds !== false,
        hasLiquidity: source.hasLiquidity === true,
      };
    }
    let leagueKeysSeen = 0;
    for (const [leagueKey, rows] of Object.entries(json.odds)) {
      const parsed = parseLeagueKey(leagueKey);
      if (!parsed || parsed.leagueId !== leagueId) continue;
      leagueKeysSeen += 1;
      if (parsed.phase !== "pregame" || !Array.isArray(rows)) continue;
      for (const row of rows) ingestSnapshotRow(state, row, leagueId, parsed.periodTypeId, counts);
    }
    // No lg<id> key while other leagues' keys exist is the wrong file or a
    // schema change; an empty `odds` map (Serie B off-season, 2026-09-10) or
    // keys with no game rows is an empty slate and must not read as a failure.
    if (leagueKeysSeen === 0 && Object.keys(json.odds).length > 0) {
      throw new Error(`snapshot: no lg${leagueId} odds keys in the file (keys: ${Object.keys(json.odds).slice(0, 5).join(", ") || "none"})`);
    }
    state.counts = counts;
    return state;
  }

  function mergeStates(states) {
    const merged = emptyState();
    for (const state of states) {
      merged.leagues.push(...state.leagues);
      Object.assign(merged.teams, state.teams);
      Object.assign(merged.teamIndex, state.teamIndex || {});
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
      isAlt: false,
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
  // number is newer. Snapshot-only fields (liquidity, fromSnapshot) carry over
  // on replace: an update shares its snapshot line's key, so it is the same
  // market; a line the stream ADDS may be another market of the event.
  // Alt lines are never in the stream, so they stay as the last snapshot left
  // them; selectEdges reads the main line's CURRENT points for the distance
  // gate and the same-points dedupe, so a main line that moves onto an alt's
  // number hides that alt until the next snapshot refresh replaces it.
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
      state.lines[line.key] = { liquidity: held ? held.liquidity : null, fromSnapshot: held ? held.fromSnapshot === true : false, ...fields };
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
      awayTeamId: event.awayTeamId ?? null,
      homeTeamId: event.homeTeamId ?? null,
      homeAway: line.sideIndex === 0 ? "Away" : "Home",
      rotation: line.sideIndex === 0 ? event.awayRotation : event.homeRotation,
      awayRotation: event.awayRotation,
      homeRotation: event.homeRotation,
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
      // When the book last changed this line (modifiedOn; an alt's sequenceNumber); null if unknown.
      modifiedMs: lineChangedMs(line),
      isAlt: line.isAlt === true,
      // The book's main-line points this alt hangs off (current main line when held); null on a main line.
      mainPoints: line.isAlt ? currentMainPoints(line, state) : null,
      // The EVENT's venue id map (noteRungVenueIds), by reference — every row
      // of the event shares it, so the bet matcher can join a bet on its
      // Kalshi / Novig id (#118 step 3). Null for an event with no snapshot.
      venueIds: event.venueIds ?? null,
    };
  }

  function currentMainPoints(altLine, state) {
    const main = state.lines[altLine.mainKey];
    return main && main.points != null ? main.points : altLine.mainPoints;
  }

  function positiveNumberOrNull(value) {
    return typeof value === "number" && Number.isFinite(value) && value > 0 ? value : null;
  }

  // Alt-only gates. An alt is listed only when includeAlts is on, the main
  // line is not currently sitting on the same number (same bet twice), and it
  // is within altMaxDistance points of the book's current main number (deep
  // ladders are extrapolated fairs and a few-dollar stake).
  function altPassesGates(line, state, opts) {
    if (!opts.includeAlts) return false;
    const main = state.lines[line.mainKey];
    if (main && main.points === line.points) return false;
    if (opts.altMaxDistance != null) {
      const mainPoints = currentMainPoints(line, state);
      if (mainPoints == null || Math.abs(line.points - mainPoints) > opts.altMaxDistance) return false;
    }
    return true;
  }

  // Dollars won per dollar staked at an American price: +2000 -> 20, -110 -> 0.909.
  function winPerDollarStaked(americanPrice) {
    return americanPrice > 0 ? americanPrice / 100 : 100 / Math.abs(americanPrice);
  }

  // A line that reports liquidity, i.e. an exchange, is listed only when the
  // money resting at its price can win at least minLiquidityToWin: $20 at
  // +2000 wins $400 and is worth a look, $20 at +100 wins $20 and is not
  // (Cal, 2026-09-23). A flat stake floor hid longshots whose whole Kelly bet
  // is small. Main lines and alts alike; books with no liquidity figure pass.
  function liquidityCanWin(line, minLiquidityToWin) {
    if (minLiquidityToWin == null || line.liquidity == null) return true;
    return line.liquidity * winPerDollarStaked(line.price) >= minLiquidityToWin;
  }

  // Lines worth listing: on the board, edge known and >= minEdge (a fraction),
  // period/bet type enabled, book allowed, game not started, and — when
  // maxLineAgeMs is set — changed by the book within that window (a 96-day-old
  // line at a "live" book is a dead feed, and its 36% "edge" is not bettable;
  // a line whose change time is unknowable is excluded too). Alt lines
  // additionally pass altPassesGates. Every line passes liquidityCanWin.
  // Sorted by edge.
  function selectEdges(state, options) {
    const opts = options || {};
    const minEdge = typeof opts.minEdge === "number" ? opts.minEdge : 0.01;
    const periods = opts.periods instanceof Set ? opts.periods : new Set([1]);
    const betTypes = opts.betTypes instanceof Set ? opts.betTypes : new Set([1, 2, 3]);
    const bookIds = opts.bookIds instanceof Set ? opts.bookIds : null;
    const now = typeof opts.now === "number" ? opts.now : Date.now();
    const maxLineAgeMs = positiveNumberOrNull(opts.maxLineAgeMs);
    const minLiquidityToWin = positiveNumberOrNull(opts.minLiquidityToWin);
    const altOpts = {
      includeAlts: opts.includeAlts === true,
      altMaxDistance: positiveNumberOrNull(opts.altMaxDistance),
    };
    const rows = [];
    for (const line of Object.values(state.lines)) {
      if (line.bookId === UNABATED_LINE_BOOK_ID) continue;
      if (line.isAlt && !altPassesGates(line, state, altOpts)) continue;
      if (line.statusId !== STATUS_ON_BOARD) continue;
      if (line.ge == null || line.ge < minEdge) continue;
      if (!periods.has(line.periodTypeId) || !betTypes.has(line.betTypeId)) continue;
      const book = state.books[line.bookId];
      if (bookIds ? !bookIds.has(line.bookId) : !(book && book.isLive)) continue;
      const event = state.events[line.eventId];
      if (!event || event.eventStart == null || event.eventStart <= now) continue;
      if (line.betTypeId !== 1 && line.points == null) continue;
      // Not a valid American price: nothing downstream (cents, Kelly) can use it.
      if (Math.abs(line.price) < 100) continue;
      if (!liquidityCanWin(line, minLiquidityToWin)) continue;
      if (maxLineAgeMs != null) {
        const changedMs = lineChangedMs(line);
        if (changedMs == null || now - changedMs > maxLineAgeMs) continue;
      }
      rows.push(describeLine(line, state));
    }
    rows.sort((a, b) => b.edgePct - a.edgePct || a.eventStartMs - b.eventStartMs);
    return rows;
  }

  // ---- grouping ------------------------------------------------------------

  function groupKeyOf(row) {
    return `${row.eventId}:pt${row.periodTypeId}:bt${row.betTypeId}:si${row.sideIndex}`;
  }

  // The side without its number: "Idaho Vandals", "Over".
  function sideNameOf(row) {
    if (row.betTypeId === 3) return row.sideIndex === 0 ? "Over" : "Under";
    const team = row.sideIndex === 0 ? row.awayTeam : row.homeTeam;
    return team || row.sideLabel;
  }

  // One card per (game, period, bet type, side): a +EV opinion is
  // directional, so the two sides of a market are two cards. Rows inside a
  // card sort by rankOf(row) descending (the panel passes the Kelly stake,
  // which already taxes longshots), edge as the tie-break; `best` is the
  // first. Cards come back in the same order by their best line. Rows are
  // the selectEdges output (any extra fields, e.g. stake, ride along).
  function groupEdges(rows, rankOf) {
    const rank = typeof rankOf === "function" ? rankOf : (row) => row.edgePct;
    const byKey = new Map();
    for (const row of rows) {
      const key = groupKeyOf(row);
      let group = byKey.get(key);
      if (!group) {
        group = {
          key,
          eventId: row.eventId,
          leagueId: row.leagueId,
          league: row.league,
          leagueLabel: row.leagueLabel,
          eventName: row.eventName,
          eventStart: row.eventStart,
          eventStartMs: row.eventStartMs,
          awayTeam: row.awayTeam,
          homeTeam: row.homeTeam,
          betTypeId: row.betTypeId,
          betType: row.betType,
          periodTypeId: row.periodTypeId,
          period: row.period,
          sideIndex: row.sideIndex,
          sideName: sideNameOf(row),
          rows: [],
          bookCount: 0,
          best: null,
        };
        byKey.set(key, group);
      }
      group.rows.push(row);
    }
    const compare = (a, b) => (rank(b) ?? -Infinity) - (rank(a) ?? -Infinity) || b.edgePct - a.edgePct;
    const groups = Array.from(byKey.values());
    for (const group of groups) {
      group.rows.sort(compare);
      group.best = group.rows[0];
      group.bookCount = new Set(group.rows.map((row) => row.book.id)).size;
    }
    groups.sort((a, b) => compare(a.best, b.best) || a.eventStartMs - b.eventStartMs);
    return groups;
  }

  // Main lines only; alts are counted apart so the header can say both.
  function countLines(state) {
    let count = 0;
    for (const line of Object.values(state.lines)) if (!line.isAlt) count += 1;
    return count;
  }

  function countAltLines(state) {
    let count = 0;
    for (const line of Object.values(state.lines)) if (line.isAlt) count += 1;
    return count;
  }

  const api = {
    LEAGUES, SPORTS, leagueIdsOfSport, BET_TYPES, PERIODS, UNABATED_LINE_BOOK_ID, CURSOR_EPOCH_MS,
    kalshiEventSuffixOf,
    parseLeagueKey, parseEventStart, parseModifiedOn, lineChangedMs, lineKeyOf, altLineKeyOf, emptyState, teamSpellingsFromEventName,
    parseSnapshot, mergeStates, extractCursor, cursorFromDate, parseChanges, applyChanges,
    describeLine, selectEdges, groupEdges, groupKeyOf, countLines, countAltLines,
  };

  if (typeof module !== "undefined" && module.exports) {
    module.exports = api;
  } else {
    root.UnabatedFeed = api;
  }
})(typeof globalThis !== "undefined" ? globalThis : this);
