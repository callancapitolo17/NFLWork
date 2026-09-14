// Unabated Ticket — page-world script (runs in the page's own JS world).
//
// Why MAIN world: the ticket is read from React fiber props and the AG Grid
// row data hanging off the odds-cell DOM nodes. Those expando properties are
// invisible from an isolated-world content script, so this file runs in the
// page world and hands results to content.js via window.postMessage.
//
// Side effects: one capture-phase click listener on document (Unabated's
// own handler still runs), one 5s interval while a ticket is being watched
// (resumed from the stored ticket on load, so a navigation does not end it),
// and a 10s heartbeat that also publishes the user's Unabated book selection
// (read from the grid's React context) as the Edges tab's default book filter. The only DOM touch is the locate flash:
// a 2.5s outline on the cell an Edges row or notification pointed at.

(function () {
  "use strict";

  const MESSAGE_SOURCE = "unabated-ticket";
  // Each injected copy has an id; a newer copy (re-injected after an extension
  // reload) posts a takeover and every older copy retires itself.
  const INSTANCE_ID = `${Date.now()}-${Math.random().toString(36).slice(2, 8)}`;
  let retired = false;
  const intervals = [];
  const CELL_SHELL_SELECTOR = ".odds-cell-action-shell";
  // The "..." menu button is a child of the shell; opening a menu is not picking a bet.
  const MORE_BUTTON_SELECTOR = ".odds-cell-more-button";
  const MAX_FIBER_HOPS = 40;
  const WATCH_INTERVAL_MS = 5000;
  const BET_TYPE_NAMES = { 1: "Moneyline", 2: "Spread", 3: "Total" };
  // Same names feed.PERIODS uses; the bet matcher compares a ticket's period
  // to a bet's, so an unknown id is named, never defaulted to full game.
  const PERIOD_NAMES = { 1: "FG", 2: "1H", 3: "2H", 4: "1Q", 5: "2Q", 6: "3Q", 7: "4Q" };

  // ---- messaging -----------------------------------------------------------

  function post(type, payload) {
    if (retired) return;
    window.postMessage({ source: MESSAGE_SOURCE, type, payload, instanceId: INSTANCE_ID }, window.location.origin);
  }

  // ---- React fiber helpers -------------------------------------------------

  function fiberOf(element) {
    if (!element) return null;
    const key = Object.keys(element).find((k) => k.startsWith("__reactFiber$"));
    return key ? element[key] : null;
  }

  // Walk up the fiber tree and return the first memoizedProps object that
  // satisfies `predicate`. Stops after MAX_FIBER_HOPS.
  function findProps(startFiber, predicate) {
    let fiber = startFiber;
    for (let hop = 0; fiber && hop < MAX_FIBER_HOPS; hop += 1) {
      const props = fiber.memoizedProps;
      if (props && typeof props === "object" && predicate(props)) return props;
      fiber = fiber.return;
    }
    return null;
  }

  function isLineProps(props) {
    return props.marketLine && typeof props.marketLine === "object" && "sideIndex" in props;
  }

  function isGridCellProps(props) {
    return props.api && props.node && props.node.data && typeof props.api.forEachNode === "function";
  }

  // ---- ticket building -----------------------------------------------------

  function leagueFromUrl() {
    // "/cfb/odds" -> "cfb"
    const segments = window.location.pathname.split("/").filter(Boolean);
    return segments[0] || null;
  }

  function requireNumber(value, label) {
    if (typeof value !== "number" || !Number.isFinite(value)) {
      throw new Error(`${label} is missing or not a number (got ${JSON.stringify(value)})`);
    }
    return value;
  }

  function bookPriceOf(marketLine) {
    // Exchanges (Novig, Kalshi, ProphetX) carry no americanPrice; `price` is the American price there.
    const raw = marketLine.americanPrice ?? marketLine.price;
    return requireNumber(raw, "book price");
  }

  // bacr = Unabated fair American price at this book's points; shown for
  // information only, the stake is sized from Unabated's edge.
  function fairPriceOrNull(marketLine) {
    const fair = marketLine.bacr;
    return typeof fair === "number" && Number.isFinite(fair) ? fair : null;
  }

  // Script build, so a stale copy of page.js in an old tab shows itself in the panel.
  const PAGE_SCRIPT_BUILD = "0.6.6";

  // What the clicked object actually carried, for the panel's no-edge detail:
  // decides between "Unabated never priced it" and "the field moved". Space
  // separated so the panel can wrap it.
  function describeLineFields(marketLine) {
    const show = (value) => {
      try {
        return JSON.stringify(value) ?? "undefined";
      } catch (_error) {
        return "(unserialisable)";
      }
    };
    return `keys=[${Object.keys(marketLine).slice(0, 40).join(" ")}] edge=${show(marketLine.edge)} ge=${show(marketLine.ge)} `
      + `bacr=${show(marketLine.bacr)} price=${show(marketLine.price)} americanPrice=${show(marketLine.americanPrice)} statusId=${show(marketLine.statusId)}`;
  }

  // The row's own sides entry for this book at the clicked points: the main
  // line, or the ladder rung at those points. The cell's prop can be a copy
  // the screen made without the feed's ge (live 2026-09-12: a Novig main
  // line listed on the Edges tab carried neither edge nor ge on the cell).
  function rowLineFor(marketLine, rowData, sideKey, bookKey) {
    const main = rowData.sides && rowData.sides[sideKey] && rowData.sides[sideKey][bookKey];
    if (!main) return null;
    if (main.points === marketLine.points) return main;
    return altLineAt(main, marketLine.points);
  }

  // Edge and fair for the clicked cell: the cell's object first, then the
  // row's entry for the same book and points. Both null when neither has
  // one — the panel then tries the Edges feed before calling it unpriced.
  function edgeForCell(marketLine, rowData, sideKey, bookKey) {
    const own = edgePctOf(marketLine);
    if (own != null) return { edgePct: own, fair: fairPriceOrNull(marketLine) };
    const entry = rowLineFor(marketLine, rowData, sideKey, bookKey);
    const fromRow = entry && entry !== marketLine ? edgePctOf(entry) : null;
    if (fromRow != null) return { edgePct: fromRow, fair: fairPriceOrNull(entry) };
    return { edgePct: null, fair: fairPriceOrNull(marketLine) };
  }

  // Unabated's sourceFormat: 1 = American, 2 = decimal (1.909), 4 = probability (0.525).
  // Exchanges publish probability; Unabated's `price` is that rounded to a whole
  // American number, so the panel prefers the exact source when it is 2 or 4.
  function sourcePriceOf(marketLine) {
    const format = marketLine.sourceFormat;
    const value = marketLine.sourcePrice;
    if ((format === 2 || format === 4) && typeof value === "number" && Number.isFinite(value) && value > 0) {
      return { sourceFormat: format, sourcePrice: value };
    }
    return { sourceFormat: 1, sourcePrice: null };
  }

  // The screen's computed edge when the cell carries one, else the feed's own
  // fraction on the same object (0.0296 = +2.96%) — the number the Edges tab
  // and Unabated's own % come from. Alternate-line objects never carry
  // `edge`, and main-line cells sometimes lack it too (live 2026-09-11: a
  // line listed with an edge on the Edges tab captured as "no fair"). A line
  // Unabated has not priced has neither, so this returns null.
  function edgePctOf(marketLine) {
    const edge = marketLine.edge && marketLine.edge.edge;
    if (typeof edge === "number" && Number.isFinite(edge)) return edge;
    const ge = marketLine.ge;
    return typeof ge === "number" && Number.isFinite(ge) ? Math.round(ge * 1e6) / 1e4 : null;
  }

  function bookNameOf(bookId, context, cellProps) {
    const fromCell = cellProps && cellProps.marketSource && cellProps.marketSource.id === bookId
      ? cellProps.marketSource.name
      : null;
    if (fromCell) return fromCell;
    const sources = context && context.marketSources;
    if (Array.isArray(sources)) {
      const hit = sources.find((s) => s && s.id === bookId);
      if (hit && hit.name) return hit.name;
    } else if (sources && typeof sources === "object") {
      const hit = sources[bookId] ?? sources[String(bookId)];
      if (typeof hit === "string") return hit;
      if (hit && hit.name) return hit.name;
    }
    return `book ${bookId}`;
  }

  function lineIdOf(line) {
    return line.marketLineId ?? line.id ?? null;
  }

  function bookIdOfKey(bookKey) {
    const parsed = Number(bookKey.replace(/^ms/, ""));
    return Number.isInteger(parsed) ? parsed : null;
  }

  // Which "ms<id>" entry under this side IS `marketLine` — the entry itself
  // or one of its alternateLines (same object). Exact when it hits.
  function bookIdByIdentity(marketLine, rowData, sideKey) {
    const books = rowData.sides && rowData.sides[sideKey];
    if (!books) return null;
    for (const [bookKey, line] of Object.entries(books)) {
      if (!line) continue;
      const ladder = Array.isArray(line.alternateLines) ? line.alternateLines : [];
      if (line === marketLine || ladder.includes(marketLine)) return bookIdOfKey(bookKey);
    }
    return null;
  }

  // Which "ms<id>" entry under this side carries the same line id.
  function bookIdByLineId(marketLine, rowData, sideKey) {
    const books = rowData.sides && rowData.sides[sideKey];
    const wantedId = lineIdOf(marketLine);
    if (!books || wantedId == null) return null;
    for (const [bookKey, line] of Object.entries(books)) {
      if (line && String(lineIdOf(line)) === String(wantedId)) return bookIdOfKey(bookKey);
    }
    return null;
  }

  // Identity first, then the column the cell sits in, then the line's own
  // marketSourceId: an alternate-line object can name ANOTHER book there
  // (Sports Interaction's alts carry BetMGM's id 4, feed 2026-09-11), and
  // best-line cells sit in a column with no book id, so no single field is
  // enough on its own.
  function bookIdOf(marketLine, cellProps, rowData, sideKey) {
    const byIdentity = bookIdByIdentity(marketLine, rowData, sideKey);
    if (byIdentity != null) return byIdentity;
    if (cellProps && cellProps.marketSource && typeof cellProps.marketSource.id === "number") {
      return cellProps.marketSource.id;
    }
    if (typeof marketLine.marketSourceId === "number") return marketLine.marketSourceId;
    const byLineId = bookIdByLineId(marketLine, rowData, sideKey);
    if (byLineId != null) return byLineId;
    const colId = cellProps && cellProps.colDef && cellProps.colDef.colId;
    const parsed = Number(colId);
    if (Number.isInteger(parsed)) return parsed;
    throw new Error("could not identify the book for this line");
  }

  function teamNameOf(teamId, context) {
    const teams = context && context.fullOddsData && context.fullOddsData.teams;
    const team = teams && (teams[teamId] ?? teams[String(teamId)]);
    if (team && team.name) return team.name;
    throw new Error(`team ${teamId} not found in fullOddsData.teams`);
  }

  function signedPoints(points) {
    return points > 0 ? `+${points}` : `${points}`;
  }

  function sideLabelOf(betType, sideIndex, points, rowData, context) {
    if (betType === "Total") {
      const overUnder = sideIndex === 0 ? "Over" : "Under";
      return `${overUnder} ${requireNumber(points, "total points")}`;
    }
    const eventTeam = rowData.eventTeams && rowData.eventTeams[sideIndex];
    if (!eventTeam || eventTeam.id == null) {
      throw new Error(`eventTeams[${sideIndex}] missing on the row`);
    }
    const team = teamNameOf(eventTeam.id, context);
    if (betType === "Spread") {
      return `${team} ${signedPoints(requireNumber(points, "spread points"))}`;
    }
    return team;
  }

  // The key inside rowData.sides for this side, e.g. "si0:tid771".
  function sideKeyOf(rowData, sideIndex) {
    const sides = rowData.sides || {};
    const prefix = `si${sideIndex}:`;
    const key = Object.keys(sides).find((k) => k.startsWith(prefix));
    if (!key) throw new Error(`no sides key starting with ${prefix} on the row`);
    return key;
  }

  // Team name or null; the panel falls back to eventName when a lookup fails.
  function teamNameOrNull(sideIdx, rowData, context) {
    const eventTeam = rowData.eventTeams && rowData.eventTeams[sideIdx];
    if (!eventTeam || eventTeam.id == null) return null;
    try { return teamNameOf(eventTeam.id, context); } catch (_error) { return null; }
  }

  // An alternate-line cell's marketLine is one of the main line's
  // alternateLines, not the sides entry itself; the watcher must then
  // re-find it by points inside that ladder. Null for a main-line cell (the
  // same object, or the same points, as the sides entry).
  function altPointsOf(marketLine, rowData, sideKey, bookKey) {
    const main = bookEntryOf(rowData, sideKey, bookKey);
    if (!main || main === marketLine || main.points === marketLine.points) return null;
    return typeof marketLine.points === "number" ? marketLine.points : null;
  }

  // The row's entry for this book on this side. The side key is matched by
  // its "si<index>:" prefix when the exact key is absent, so a key read off
  // an Alts child row still finds the side on the top-level row (and the
  // other way round) should the two spell the team part differently.
  function bookEntryOf(rowData, sideKey, bookKey) {
    const sides = rowData.sides;
    if (!sides || !sideKey) return null;
    let side = sides[sideKey];
    if (!side) {
      const prefix = sideKey.slice(0, sideKey.indexOf(":") + 1);
      const key = prefix ? Object.keys(sides).find((k) => k.startsWith(prefix)) : null;
      side = key ? sides[key] : null;
    }
    return (side && side[bookKey]) || null;
  }

  // Counted by rung, not array length: a per-rung row's entry can carry an
  // alternateLines array holding no rung at all (live 2026-09-12, NFL CHI@CAR:
  // the -2.5 row traced "ladder 0" yet ranked as carrying one, tying the main
  // row whose 25-rung ladder held -2.5 — the trace fired on every capture).
  function carriesLadder(entry) {
    return !!entry && Array.isArray(entry.alternateLines) && entry.alternateLines.some(Boolean);
  }

  function sameMarketRow(a, b) {
    return a.eventId === b.eventId && a.betTypeId === b.betTypeId && (a.periodTypeId ?? 1) === (b.periodTypeId ?? 1);
  }

  // A market's MAIN row is a top-level grid row. The rows of an expanded Alts
  // section are children (AG Grid detail/tree rows: `detail`, `level` > 0, a
  // parent that carries data) whose entry for the book is one rung, so any
  // lookup that takes "a row of this market" can land on a rung — live
  // 2026-09-12 the watcher followed Kalshi's 1H Under 2.5 (+2242) after a
  // capture of Under 19.5 (+265). Rank rows so a top-level row always beats a
  // child, a row carrying the market's bestLines beats one without, and a row
  // carrying this book's ladder beats one whose entry is a lone rung.
  function nodeIsTopLevel(node) {
    if (!node || node.detail === true) return false;
    if (typeof node.level === "number" && node.level > 0) return false;
    return !(node.parent && node.parent.data);
  }

  // `fitsLine(entry)` says whether the row's entry for the book IS the line
  // being resolved (capture: the clicked object; the watcher: the captured
  // number). It outranks every shape signal: a grid that lists a market's
  // rungs as sibling TOP-LEVEL rows sharing one grid key (Unabated's CFB
  // alt-lines view, live 2026-09-12: eight "Over" rows 33.5 .. 56.5, every
  // one top-level with bestLines and no ladder) ties every row on shape, and
  // grid order then picked the 33.5 row for a click on 56.5.
  function rowRank(node, sideKey, bookKey, fitsLine) {
    const data = node.data || {};
    const entry = bookEntryOf(data, sideKey, bookKey);
    const hasBestLines = !!data.bestLines && Object.keys(data.bestLines).length > 0;
    return (fitsLine && entry && fitsLine(entry) ? 8 : 0)
      + (nodeIsTopLevel(node) ? 4 : 0) + (hasBestLines ? 2 : 0) + (carriesLadder(entry) ? 1 : 0);
  }

  // The clicked line, on the row's entry itself or inside its ladder: the
  // same object, or the same number. By number as well because an Alts
  // child's entry is the rung itself and is not known to be the SAME object
  // as the parent ladder's element — on identity alone the child could
  // outrank its top-level parent, and the watcher would follow the parent's
  // main number once the Alts section closed (the 2026-09-12 bug again).
  function fitsClickedLine(marketLine) {
    const points = typeof marketLine.points === "number" ? marketLine.points : null;
    const same = (line) => line === marketLine || (points != null && line.points === points);
    return (entry) => same(entry)
      || (Array.isArray(entry.alternateLines) && entry.alternateLines.some((alt) => alt && same(alt)));
  }

  // The captured number, wherever the row carries it: any entry when the
  // market has no number (a moneyline), else the entry at that number
  // itself or its ladder's rung at it. One rule for capture's pick and for
  // every watch tick, whichever sibling row a keyed lookup answers with.
  function fitsWatchedLine(watch) {
    if (typeof watch.points !== "number") return () => true;
    return (entry) => !!lineOnEntry(entry, watch.points);
  }

  // The line at `points` on a book's entry: the entry itself when it sits
  // at that number (a main line, or a per-rung row's own rung), else the
  // rung inside its ladder, else null (the number is off the board).
  function lineOnEntry(entry, points) {
    if (typeof points !== "number") return entry;
    if (entry.points === points) return entry;
    return altLineAt(entry, points);
  }

  // One line per candidate row, for the panel's trace when the watched
  // number is not the captured one: which rows the market had, what each
  // one's entry for the book carried, and which one was chosen.
  function describeMarketNode(node, sideKey, bookKey, marketLine) {
    const data = node.data || {};
    const entry = bookEntryOf(data, sideKey, bookKey);
    const ladder = Array.isArray(entry && entry.alternateLines) ? entry.alternateLines.filter(Boolean) : [];
    const shape = `${nodeIsTopLevel(node) ? "top" : "child"}${typeof node.level === "number" ? ` L${node.level}` : ""}${node.detail ? " detail" : ""}${node.parent && node.parent.data ? " parented" : ""}`;
    const key = data.gridKey != null ? String(data.gridKey) : node.id != null ? `#${node.id}` : "?";
    const entryText = entry
      ? `entry ${entry.points ?? "-"} ${entry.americanPrice ?? entry.price ?? "?"}${entry === marketLine ? " (=clicked)" : ""}, ladder ${ladder.length}${ladder.some((alt) => alt === marketLine) ? " (has clicked)" : ""}`
      : "no entry";
    return `[${shape}; key ${key}; ${entryText}; bestLines ${data.bestLines ? Object.keys(data.bestLines).length : 0}]`;
  }

  // Every distinct grid API on the page: the master grid and, when an Alts
  // section is open, its detail grid, which mounts as its own AG Grid. One
  // fiber walk per grid root rather than per rendered cell — a busy slate
  // renders hundreds of cells and this runs on a watch tick that missed.
  const GRID_ROOT_SELECTOR = ".ag-root";
  const GRID_ROOT_LIMIT = 20;

  function apiInside(container) {
    for (const selector of GRID_PROBE_SELECTORS) {
      for (const element of container.querySelectorAll(selector)) {
        const fiber = fiberOf(element);
        const props = fiber && findProps(fiber, isGridCellProps);
        if (props && props.api) return props.api;
      }
    }
    return null;
  }

  function allGridApis(seed) {
    const apis = [];
    const seen = new Set();
    const add = (api) => {
      if (api && !seen.has(api) && !apiIsDead(api)) {
        seen.add(api);
        apis.push(api);
      }
    };
    add(seed);
    const roots = Array.from(document.querySelectorAll(GRID_ROOT_SELECTOR)).slice(0, GRID_ROOT_LIMIT);
    // No .ag-root (a grid that mounts differently): fall back to the cell scan.
    for (const root of roots.length ? roots : [document]) add(apiInside(root));
    return apis;
  }

  // Every row of this market across the reachable grids, best-ranked first.
  // `apis` lets a caller that already probed the page reuse that list.
  function rankedMarketNodes(seedApi, identity, sideKey, bookKey, apis, fitsLine) {
    const found = [];
    for (const api of apis || allGridApis(seedApi)) {
      api.forEachNode((node) => {
        if (!node.data || !sameMarketRow(node.data, identity)) return;
        found.push({ node, api, rank: rowRank(node, sideKey, bookKey, fitsLine), order: found.length });
      });
    }
    found.sort((a, b) => b.rank - a.rank || a.order - b.order);
    return found;
  }

  // The row to classify, price and watch against: the market's best-ranked
  // row that carries an entry for the book, else the clicked row itself.
  // `trace` lists every candidate for the panel. `ambiguous` means two rows
  // TIED at the best rank, so the pick came down to grid order — deliberately
  // not "the pick was not top-level", which would be true of every row on a
  // grid that groups its rows and would show the trace on every ticket.
  // `top` is the pick's shape, which the watcher compares against later.
  function ladderRowFor(gridApi, rowData, sideKey, bookKey, marketLine) {
    const own = { rowData, gridApi, trace: "clicked row only", ambiguous: false, top: true };
    // Without an event id the identity match would accept any row; keep the cell's own.
    if (rowData.eventId == null) return own;
    const ranked = rankedMarketNodes(gridApi, rowData, sideKey, bookKey, undefined, fitsClickedLine(marketLine));
    const trace = ranked.slice(0, 8).map(({ node }) => describeMarketNode(node, sideKey, bookKey, marketLine)).join(" ");
    const withEntry = ranked.filter(({ node }) => bookEntryOf(node.data, sideKey, bookKey));
    if (!withEntry.length) return { ...own, trace: `no row carries ${bookKey}: ${trace}` };
    const pick = withEntry[0];
    const tied = withEntry.filter((candidate) => candidate.rank === pick.rank).length > 1;
    return {
      rowData: pick.node.data,
      gridApi: pick.api,
      trace: `picked ${describeMarketNode(pick.node, sideKey, bookKey, marketLine)} of ${trace}`,
      ambiguous: tied,
      top: nodeIsTopLevel(pick.node),
    };
  }

  function buildTicket({ marketLine, sideIndex, rowData: cellRowData, context, cellProps, gridApi: cellGridApi }) {
    const betTypeId = cellRowData.betTypeId;
    const betType = BET_TYPE_NAMES[betTypeId];
    if (!betType) throw new Error(`unsupported betTypeId ${betTypeId} (only moneyline, spread, total)`);

    const points = marketLine.points ?? null;
    const sideKey = sideKeyOf(cellRowData, sideIndex);
    const bookId = bookIdOf(marketLine, cellProps, cellRowData, sideKey);
    const bookKey = `ms${bookId}`;
    const { rowData, gridApi, trace: rowTrace, ambiguous: rowAmbiguous, top: rowTop } = ladderRowFor(cellGridApi, cellRowData, sideKey, bookKey, marketLine);
    const altPoints = altPointsOf(marketLine, rowData, sideKey, bookKey);
    const edge = edgeForCell(marketLine, rowData, sideKey, bookKey);
    const rotation = rowData.eventTeams && rowData.eventTeams[sideIndex]
      ? rowData.eventTeams[sideIndex].rotationNumber ?? null
      : null;

    return {
      capturedAt: Date.now(),
      league: leagueFromUrl(),
      eventId: rowData.eventId ?? null,
      eventStart: rowData.eventStart ?? null,
      eventName: rowData.eventName ?? null,
      betType,
      period: PERIOD_NAMES[rowData.periodTypeId ?? 1] || `pt${rowData.periodTypeId}`,
      sideIndex,
      sideLabel: sideLabelOf(betType, sideIndex, points, rowData, context),
      // Side 0 is the away team / Over, side 1 the home team / Under.
      awayTeam: teamNameOrNull(0, rowData, context),
      homeTeam: teamNameOrNull(1, rowData, context),
      homeAway: sideIndex === 0 ? "Away" : "Home",
      rotation,
      points,
      book: { id: bookId, name: bookNameOf(bookId, context, cellProps) },
      price: bookPriceOf(marketLine),
      ...sourcePriceOf(marketLine),
      fair: edge.fair,
      // Null when neither the cell nor the row carried an edge: the panel
      // sizes from the Edges feed's copy of this line when it has one at the
      // same price, else shows noEdgeDetail as "No Unabated fair".
      edgePct: edge.edgePct,
      noEdgeDetail: edge.edgePct == null ? `script ${PAGE_SCRIPT_BUILD}; cell fields ${describeLineFields(marketLine)}` : null,
      isAlt: altPoints != null,
      // Which grid row the ticket was classified and priced against, and the
      // rows it was chosen from; the panel shows it when the watched number
      // is not the captured one, or when the pick was not a lone top-level row.
      rowResolution: { trace: rowTrace, ambiguous: rowAmbiguous, build: PAGE_SCRIPT_BUILD },
      // Watcher handle: how to find this same line again through the grid API.
      watch: {
        gridKey: rowData.gridKey ?? null, sideKey, bookKey,
        // The picked row's shape; a watch tick reading a row of the OTHER
        // shape is reading a different rung, whatever the grid's layout.
        rowTop,
        // The captured number: every watch tick re-finds the line by it
        // (the row's entry at it, or its ladder's rung), so a grid that
        // lists rungs as sibling rows under one key cannot hand the watcher
        // another rung, and a number gone from the row reads off the board.
        points,
        // Row identity for when the grid key stops resolving (rebuilt grid, re-keyed row).
        eventId: rowData.eventId ?? null, betTypeId, periodTypeId: rowData.periodTypeId ?? 1,
      },
      current: null,
      // The grid API the ladder row lives on: what the watcher must poll.
      gridApi,
    };
  }

  // ---- reading the clicked cell ---------------------------------------------

  function sideIndexOf(shell, lineProps) {
    const attr = shell.getAttribute("data-side-index");
    const fromAttr = attr == null ? NaN : Number(attr);
    if (Number.isInteger(fromAttr)) return fromAttr;
    if (lineProps && Number.isInteger(lineProps.sideIndex)) return lineProps.sideIndex;
    throw new Error("side index missing on the clicked cell");
  }

  // Fast path: the cell's own fiber props carry marketLine + sideIndex + context.
  function readViaFiber(shell) {
    const fiber = fiberOf(shell);
    if (!fiber) throw new Error("no React fiber on the clicked cell");
    const lineProps = findProps(fiber, isLineProps);
    if (!lineProps) throw new Error("no marketLine props found above the clicked cell");
    const cellProps = findProps(fiber, isGridCellProps);
    if (!cellProps) throw new Error("no AG Grid cell props (api/node) found above the clicked cell");
    return {
      marketLine: lineProps.marketLine,
      sideIndex: sideIndexOf(shell, lineProps),
      rowData: cellProps.node.data,
      context: lineProps.context ?? cellProps.context ?? null,
      cellProps,
      gridApi: cellProps.api,
    };
  }

  // AG Grid's own .ag-cell wrappers are not React-rendered on Unabated (no
  // fiber on them); the price shells inside are, and their fiber walks up to
  // the cell renderer props (api, node, context) — the same path capture uses.
  // .ag-cell stays as a fallback for a grid that mounts differently.
  const GRID_PROBE_SELECTORS = [CELL_SHELL_SELECTOR, ".ag-cell"];

  function anyGridApi() {
    for (const selector of GRID_PROBE_SELECTORS) {
      for (const element of document.querySelectorAll(selector)) {
        const fiber = fiberOf(element);
        if (!fiber) continue;
        const props = findProps(fiber, isGridCellProps);
        if (!props) continue;
        // The line renderer's context carries userSettings; the grid context is the fallback.
        const lineProps = findProps(fiber, isLineProps);
        const context = (lineProps && lineProps.context) || props.context || null;
        return { api: props.api, context };
      }
    }
    throw new Error(`could not reach the AG Grid API from any rendered cell (tried ${GRID_PROBE_SELECTORS.join(", ")})`);
  }

  // Fallback: data-marketline-id on the shell + a scan of every row's sides.
  function readViaMarketLineId(shell) {
    const marketLineId = shell.getAttribute("data-marketline-id");
    if (!marketLineId) throw new Error("no data-marketline-id on the clicked cell");
    const { api, context } = anyGridApi();
    const sideIndex = sideIndexOf(shell, null);
    let hit = null;
    api.forEachNode((node) => {
      if (hit || !node.data || !node.data.sides) return;
      for (const [sideKey, books] of Object.entries(node.data.sides)) {
        if (!sideKey.startsWith(`si${sideIndex}:`) || !books) continue;
        for (const line of Object.values(books)) {
          if (line && String(lineIdOf(line)) === String(marketLineId)) {
            hit = { marketLine: line, rowData: node.data };
            return;
          }
        }
      }
    });
    if (!hit) throw new Error(`marketLineId ${marketLineId} not found in any row`);
    return { marketLine: hit.marketLine, sideIndex, rowData: hit.rowData, context, cellProps: null, gridApi: api };
  }

  function readCell(shell) {
    try {
      return readViaFiber(shell);
    } catch (fiberError) {
      try {
        return readViaMarketLineId(shell);
      } catch (fallbackError) {
        throw new Error(`fiber: ${fiberError.message}; fallback: ${fallbackError.message}`);
      }
    }
  }

  // ---- line watcher --------------------------------------------------------

  let watcher = null; // { ticket, gridApi, timer }

  function stopWatching() {
    if (watcher && watcher.timer) clearInterval(watcher.timer);
    watcher = null;
  }

  function apiIsDead(api) {
    return !api || (typeof api.isDestroyed === "function" && api.isDestroyed());
  }

  function rowCountOf(api) {
    try {
      return typeof api.getDisplayedRowCount === "function" ? api.getDisplayedRowCount() : "?";
    } catch (_error) {
      return "?";
    }
  }

  // The watched row, surviving what the plain grid-key lookup does not: a
  // grid Unabated rebuilt (the held API answers with no rows — re-acquire one
  // from the DOM), a row it re-keyed (its key embeds flags like livefalse —
  // find it by event, bet type and period), and the tab having been steered
  // to another league by a locate (say which, instead of "row gone").
  function watchedRowNode() {
    const { ticket } = watcher;
    const { gridKey, eventId, betTypeId, periodTypeId } = ticket.watch;
    const shownLeague = leagueFromUrl();
    if (ticket.league && shownLeague && shownLeague !== ticket.league) {
      throw new Error(`the Unabated tab is showing ${shownLeague.toUpperCase()}; this line is on ${ticket.league.toUpperCase()}`);
    }
    if (!gridKey && eventId == null) throw new Error("no gridKey on the ticket");
    const { sideKey, bookKey } = ticket.watch;
    const identity = { eventId, betTypeId, periodTypeId: periodTypeId ?? 1 };
    // The key lookup is trusted only when it answers with a row of the SAME
    // shape capture picked: a re-keyed grid can answer with an Alts child row,
    // whose entry for the book is one rung. Compared rather than forced to
    // top-level, because capture legitimately picks a child when the book
    // prices no main line — forcing it would send every tick down the rescan.
    const wantTop = ticket.watch.rowTop ?? true;
    const fitsLine = fitsWatchedLine(ticket.watch);
    // ... and, when the keyed row carries the book, with the captured line
    // on it: rungs listed as sibling rows share one key, so the key alone
    // can answer with another rung.
    const byKey = (api) => {
      const node = gridKey ? api.getRowNode(gridKey) : null;
      if (!node || !node.data || nodeIsTopLevel(node) !== wantTop) return null;
      const entry = bookEntryOf(node.data, sideKey, bookKey);
      return entry && !fitsLine(entry) ? null : node;
    };
    let api = watcher.gridApi;
    let node = apiIsDead(api) ? null : byKey(api);
    if (!node) {
      // The held API may belong to a grid that no longer exists; a rendered
      // cell always reaches the live ones. Probed once and reused below.
      const apis = allGridApis(api);
      if (apis.length === 0) throw new Error("no odds grid reachable on the page (still loading, or Unabated changed its grid)");
      for (const other of apis) {
        node = byKey(other);
        if (node) {
          api = other;
          watcher.gridApi = api;
          break;
        }
      }
      if (!node && eventId != null) {
        // Best-ranked row that actually carries this book, as capture did:
        // the top-ranked row overall may not price this book at all, and
        // taking it would report "book line no longer on the row" while a
        // row that carries it sits further down.
        const ranked = rankedMarketNodes(api, identity, sideKey, bookKey, apis, fitsLine);
        const best = ranked.find((candidate) => bookEntryOf(candidate.node.data, sideKey, bookKey));
        if (best) {
          node = best.node;
          api = best.api;
          watcher.gridApi = api;
        }
      }
    }
    if (!node || !node.data) {
      throw new Error(`row no longer in the grid (${rowCountOf(api)} rows shown, event ${eventId} bt${betTypeId} pt${periodTypeId} not among them)`);
    }
    if (node.data.gridKey && node.data.gridKey !== gridKey) ticket.watch.gridKey = node.data.gridKey;
    return node;
  }

  function readWatchedLine() {
    const { ticket } = watcher;
    const { sideKey, bookKey } = ticket.watch;
    const node = watchedRowNode();
    const bookLine = bookEntryOf(node.data, sideKey, bookKey);
    if (!bookLine) throw new Error("book line no longer on the row");
    // Which row this read came from, and whether its shape matches the row
    // capture picked: the panel shows the trace when the two disagree.
    const row = describeMarketNode(node, sideKey, bookKey, null);
    const rowTop = nodeIsTopLevel(node);
    const line = lineOnEntry(bookLine, ticket.watch.points);
    // Neither the entry nor its ladder offers that number: off the board at the captured price.
    if (!line) {
      return {
        price: ticket.price, sourceFormat: ticket.sourceFormat, sourcePrice: ticket.sourcePrice,
        points: ticket.points, fair: ticket.fair, edgePct: ticket.edgePct, offBoard: true, seenAt: Date.now(), row, rowTop,
      };
    }
    return {
      price: bookPriceOf(line),
      ...sourcePriceOf(line),
      points: line.points ?? null,
      fair: fairPriceOrNull(line),
      edgePct: edgePctOf(line),
      offBoard: line.statusId === 2,
      seenAt: Date.now(),
      row,
      rowTop,
    };
  }

  function altLineAt(bookLine, points) {
    const ladder = Array.isArray(bookLine.alternateLines) ? bookLine.alternateLines : [];
    return ladder.find((alt) => alt && alt.points === points) || null;
  }

  function watchTick() {
    if (!watcher) return;
    try {
      post("watch", { capturedAt: watcher.ticket.capturedAt, line: readWatchedLine() });
    } catch (error) {
      post("watch", { capturedAt: watcher.ticket.capturedAt, error: error.message });
    }
  }

  function startWatching(ticket, gridApi) {
    stopWatching();
    watcher = { ticket, gridApi, timer: setInterval(watchTick, WATCH_INTERVAL_MS) };
    intervals.push(watcher.timer);
  }


  // ---- Unabated book selection publish --------------------------------------

  // The odds screen's book selection lives under context.userSettings.gameOdds
  // as entries carrying isUnavailable (false = the book is shown). Live on
  // 2026-09-10 gameOdds had 6 top-level entries with no isUnavailable on them,
  // so the book entries sit one level down; search up to 3 levels for the
  // first array/object whose members carry the flag, and when nothing does,
  // report the keys seen so the panel header shows the real shape.
  const BOOK_ENTRY_SEARCH_DEPTH = 3;

  function bookEntriesOf(node, depth) {
    if (!node || typeof node !== "object" || depth > BOOK_ENTRY_SEARCH_DEPTH) return null;
    const members = Array.isArray(node) ? node.map((entry) => ({ entry, key: null })) : Object.entries(node).map(([key, entry]) => ({ entry, key }));
    const flagged = members.filter(({ entry }) => entry && typeof entry === "object" && "isUnavailable" in entry);
    if (flagged.length) return flagged;
    for (const { entry } of members) {
      const found = bookEntriesOf(entry, depth + 1);
      if (found) return found;
    }
    return null;
  }

  function describeShape(node) {
    if (Array.isArray(node)) return `array[${node.length}]${node.length ? ` of {${Object.keys(node[0] || {}).slice(0, 8).join(",")}}` : ""}`;
    if (node && typeof node === "object") return `{${Object.keys(node).slice(0, 12).join(",")}}`;
    return typeof node;
  }

  function enabledBookIdsOf(userSettings) {
    const gameOdds = userSettings && userSettings.gameOdds;
    if (!gameOdds || typeof gameOdds !== "object") {
      throw new Error(`userSettings.gameOdds missing (userSettings keys: ${Object.keys(userSettings || {}).slice(0, 12).join(",") || "none"})`);
    }
    const entries = bookEntriesOf(gameOdds, 0);
    if (!entries) throw new Error(`no isUnavailable entries under userSettings.gameOdds; shape ${describeShape(gameOdds)}`);
    const ids = [];
    for (const { entry, key } of entries) {
      if (entry.isUnavailable !== false) continue;
      const id = Number(entry.marketSourceId ?? entry.id ?? key);
      if (Number.isInteger(id)) ids.push(id);
    }
    if (!ids.length) throw new Error(`no enabled books among ${entries.length} isUnavailable entries under userSettings.gameOdds (first: ${describeShape(entries[0].entry)})`);
    return ids;
  }

  let lastFiltersSignature = null;

  // What the selection was read from, for the panel's click-to-expand line:
  // the gameOdds entry fields with true/false counts per boolean field (so a
  // wrong flag shows up as "33 of 33 false") and the first entry.
  function filterDiagnostic(userSettings) {
    const out = { userSettingsKeys: Object.keys(userSettings || {}).slice(0, 20) };
    const gameOdds = userSettings && userSettings.gameOdds;
    out.gameOddsShape = describeShape(gameOdds);
    const entries = bookEntriesOf(gameOdds, 0) || [];
    out.entryCount = entries.length;
    const counts = {};
    for (const { entry } of entries) {
      for (const [key, value] of Object.entries(entry)) {
        if (typeof value !== "boolean") continue;
        counts[key] = counts[key] || { true: 0, false: 0 };
        counts[key][value ? "true" : "false"] += 1;
      }
    }
    out.booleanFields = counts;
    out.firstEntry = entries.length ? JSON.stringify(entries[0].entry).slice(0, 400) : null;
    return out;
  }

  function publishFilters() {
    let payload;
    let context = null;
    try {
      context = anyGridApi().context;
      const bookIds = enabledBookIdsOf(context && context.userSettings);
      payload = { bookIds, error: null, url: window.location.href, at: Date.now() };
    } catch (error) {
      payload = { bookIds: null, error: error.message, url: window.location.href, at: Date.now() };
    }
    try {
      payload.debug = context ? filterDiagnostic(context.userSettings) : { error: "no grid context" };
    } catch (error) {
      payload.debug = { error: error.message };
    }
    const signature = JSON.stringify([payload.bookIds, payload.error]);
    if (signature !== lastFiltersSignature) {
      lastFiltersSignature = signature;
      console.info("[unabated-ticket] books filter", payload);
    }
    post("filters", payload);
  }


  // ---- locate: scroll the grid to a line and flash its cell -----------------

  const LOCATE_ATTEMPTS = 20;
  const LOCATE_RETRY_MS = 1000;
  // After expanding a row's Alts, its cells mount over a few frames.
  const ALT_CELL_ATTEMPTS = 8;
  const ALT_CELL_RETRY_MS = 250;
  const FLASH_MS = 2500;
  const FLASH_ATTR = "data-unabated-ticket-flash";
  let lastLocateAt = 0;

  // The market's main row (top-level, never an Alts child) on any reachable grid.
  function findRowNode(api, request) {
    const best = rankedMarketNodes(api, request, null, null)[0];
    return best ? best.node : null;
  }

  function cellShellFor(node, request) {
    const rowId = node.data.gridKey ?? node.id;
    if (rowId == null) return null;
    const rowSelector = `.ag-row[row-id="${String(rowId).replace(/"/g, '\\"')}"]`;
    const inBookColumn = document.querySelector(`${rowSelector} .ag-cell[col-id="${request.bookId}"] ${CELL_SHELL_SELECTOR}[data-side-index="${request.sideIndex}"]`);
    if (inBookColumn) return inBookColumn;
    // Book column may be scrolled out / hidden: fall back to any shell on the row for that side.
    return document.querySelector(`${rowSelector} ${CELL_SHELL_SELECTOR}[data-side-index="${request.sideIndex}"]`);
  }

  // Alt cells live in the row's expanded Alts section (AG Grid master/detail:
  // node.setExpanded). Rendered shells are matched on their fiber props —
  // points, side, book and, when the cell's row data carries them, event and
  // bet type — never on DOM position, which differs per layout.
  function expandAlts(node) {
    if (typeof node.setExpanded !== "function") throw new Error("this grid row cannot be expanded (no setExpanded on the row node)");
    if (node.expanded !== true) node.setExpanded(true);
  }

  function shellBookId(marketLine, cellProps, rowData, sideKey) {
    try {
      return bookIdOf(marketLine, cellProps, rowData, sideKey);
    } catch (_error) {
      return null;
    }
  }

  function altCellShellFor(request) {
    for (const shell of document.querySelectorAll(CELL_SHELL_SELECTOR)) {
      const fiber = fiberOf(shell);
      const lineProps = fiber && findProps(fiber, isLineProps);
      if (!lineProps || !lineProps.marketLine || lineProps.marketLine.points !== request.points) continue;
      const attr = shell.getAttribute("data-side-index");
      const sideIndex = attr == null ? lineProps.sideIndex : Number(attr);
      if (sideIndex !== request.sideIndex) continue;
      const cellProps = findProps(fiber, isGridCellProps);
      const rowData = (cellProps && cellProps.node && cellProps.node.data) || {};
      if (rowData.eventId != null && rowData.eventId !== request.eventId) continue;
      if (rowData.betTypeId != null && rowData.betTypeId !== request.betTypeId) continue;
      if (shellBookId(lineProps.marketLine, cellProps, rowData, request.sideKey) !== request.bookId) continue;
      return shell;
    }
    return null;
  }

  function locateAltCell(node, request, attempt) {
    if (request.at !== lastLocateAt) return;
    const shell = altCellShellFor(request);
    if (shell) {
      flash(shell);
      reportLocate(request, true, null);
      return;
    }
    if (attempt < ALT_CELL_ATTEMPTS) {
      setTimeout(() => locateAltCell(node, request, attempt + 1), ALT_CELL_RETRY_MS);
      return;
    }
    // Loud, and still useful: outline the main-line cell so the row is found.
    const mainShell = cellShellFor(node, request);
    if (mainShell) flash(mainShell);
    reportLocate(request, false, `row expanded but no ${request.bookName} cell at ${request.points} rendered${mainShell ? " (its main-line cell is outlined; open the row's Alts and look for that number)" : ""}`);
  }

  function flash(element) {
    element.setAttribute(FLASH_ATTR, "1");
    const previousOutline = element.style.outline;
    const previousOffset = element.style.outlineOffset;
    element.style.outline = "3px solid #f59e0b";
    element.style.outlineOffset = "1px";
    element.scrollIntoView({ block: "center", inline: "center" });
    setTimeout(() => {
      element.style.outline = previousOutline;
      element.style.outlineOffset = previousOffset;
      element.removeAttribute(FLASH_ATTR);
    }, FLASH_MS);
  }

  function reportLocate(request, ok, message) {
    post("located", { key: request.key, sideLabel: request.sideLabel, bookName: request.bookName, leagueLabel: request.leagueLabel, ok, message, at: Date.now() });
  }

  // Retries while the grid loads (a navigated tab has no rows for a few seconds).
  function locateLine(request, attempt) {
    if (request.at !== lastLocateAt) return; // superseded by a newer request
    let api = null;
    let node = null;
    try {
      api = anyGridApi().api;
      node = findRowNode(api, request);
    } catch (_error) {
      api = null;
    }
    if (!node) {
      if (attempt < LOCATE_ATTEMPTS) {
        setTimeout(() => locateLine(request, attempt + 1), LOCATE_RETRY_MS);
        return;
      }
      reportLocate(request, false, api
        ? "row is not on the grid (hidden by your bet-type or period filter, or the game left the board)"
        : "the odds grid never appeared on this tab");
      return;
    }
    if (request.isAlt) {
      try {
        expandAlts(node);
      } catch (error) {
        reportLocate(request, false, `could not open the row's Alts: ${error.message}`);
        return;
      }
    }
    try {
      if (typeof api.ensureNodeVisible === "function") api.ensureNodeVisible(node, "middle");
      if (typeof api.ensureColumnVisible === "function") api.ensureColumnVisible(String(request.bookId));
    } catch (error) {
      reportLocate(request, false, `grid scroll failed: ${error.message}`);
      return;
    }
    if (request.isAlt) {
      locateAltCell(node, request, 0);
      return;
    }
    // The row renders on the next frame after ensureNodeVisible.
    setTimeout(() => {
      const shell = cellShellFor(node, request);
      if (!shell) {
        reportLocate(request, false, "row found but its cell did not render (book column hidden?)");
        return;
      }
      flash(shell);
      reportLocate(request, true, null);
    }, 50);
  }

  function onLocateMessage(payload) {
    if (!payload || typeof payload.at !== "number" || payload.at <= lastLocateAt) return;
    lastLocateAt = payload.at;
    const wantedPath = `/${payload.league}/`;
    if (!window.location.pathname.startsWith(wantedPath)) return; // another tab (or this one, mid-navigation) will handle it
    locateLine(payload, 0);
  }

  // The ticket outlives this script: content.js keeps it in storage while
  // this copy dies with any navigation (one-click betting leaving the tab,
  // Back, a tab reload or discard, an extension reload's takeover). Without
  // this the panel showed "Not watching the line" for every ticket after
  // such an event, with a live heartbeat and nothing to re-click for. The
  // watcher is rebuilt from the stored ticket's identity; watchedRowNode
  // finds the grid from the DOM and the row by event/bet type/period/number,
  // so no grid API from the capture is needed.
  // Not resumed: a ticket for a league this tab is not showing (a second
  // tab on another league would otherwise post "wrong league" every 5 s
  // against the tab that can read it), and one whose game has started (its
  // pre-game line is gone; a scan every 5 s for it would never end).
  // The grid row's eventStart is naive UTC ("2026-09-13T17:00:00"), and
  // Date.parse reads a naive date-time as LOCAL time: 7 h late in PDT.
  // Same rule as feed.parseEventStart and betsview.ticketAsLine.
  function eventStartMs(value) {
    if (typeof value !== "string" || !value) return null;
    const hasZone = /(?:Z|[+-]\d\d:\d\d)$/.test(value);
    const ms = Date.parse(hasZone ? value : `${value}Z`);
    return Number.isFinite(ms) ? ms : null;
  }

  function onResumeMessage(ticket) {
    if (retired || !ticket || typeof ticket.capturedAt !== "number" || !ticket.watch) return;
    // Already watching it, or watching a NEWER capture: an offer read from
    // storage just before a fresh click on this tab must not replace it.
    if (watcher && watcher.ticket.capturedAt >= ticket.capturedAt) return;
    const shownLeague = leagueFromUrl();
    if (ticket.league && shownLeague && shownLeague !== ticket.league) return;
    const startMs = eventStartMs(ticket.eventStart);
    if (startMs != null && startMs <= Date.now()) return;
    startWatching(ticket, null);
    console.info("[unabated-ticket] resumed watching", ticket.sideLabel, "captured", new Date(ticket.capturedAt).toLocaleTimeString());
  }

  function retire() {
    retired = true;
    stopWatching();
    for (const timer of intervals) clearInterval(timer);
    document.removeEventListener("pointerdown", onClickCapture, true);
    document.removeEventListener("click", onClickCapture, true);
    console.info("[unabated-ticket] page.js retired (a newer copy took over)");
  }

  window.addEventListener("message", (event) => {
    if (event.source !== window) return;
    const data = event.data;
    if (!data || data.source !== MESSAGE_SOURCE) return;
    if (data.type === "takeover" && data.instanceId !== INSTANCE_ID && !retired) retire();
    if (retired) return;
    if (data.type === "locate") onLocateMessage(data.payload);
    if (data.type === "resume") onResumeMessage(data.payload);
  });

  // ---- click capture -------------------------------------------------------

  // One-click betting: Unabated's own handler opens the book's deeplink, and
  // it can do so on mouse-down and can navigate THIS tab away before a `click`
  // event ever fires. So capture on pointerdown (capture phase on document runs
  // before any page handler) and keep `click` only as a fallback for keyboard
  // or synthetic activation. The two are deduped per shell.
  const CAPTURE_DEDUPE_MS = 1500;
  const PRIMARY_BUTTON = 0;
  const MIDDLE_BUTTON = 1;
  let lastCapture = { shell: null, at: 0 };

  function onClickCapture(event) {
    if (retired) return;
    const target = event.target instanceof Element ? event.target : null;
    const shell = target && target.closest(CELL_SHELL_SELECTOR);
    if (!shell) return;
    if (target.closest(MORE_BUTTON_SELECTOR)) return;
    if (event.type === "pointerdown" && event.button !== PRIMARY_BUTTON && event.button !== MIDDLE_BUTTON) return;
    if (lastCapture.shell === shell && Date.now() - lastCapture.at < CAPTURE_DEDUPE_MS) return;
    lastCapture = { shell, at: Date.now() };
    try {
      const cell = readCell(shell);
      const { gridApi, ...ticket } = buildTicket(cell);
      startWatching(ticket, gridApi);
      console.info("[unabated-ticket] captured", ticket);
      post("ticket", ticket);
    } catch (error) {
      stopWatching();
      console.info("[unabated-ticket] capture failed:", error.message);
      post("error", { message: error.message, kind: error.kind || "read_failed", at: Date.now() });
    }
  }

  document.addEventListener("pointerdown", onClickCapture, true);
  document.addEventListener("click", onClickCapture, true);
  // Row resolution under test (tests/page_rows.test.js): the harness sets
  // this flag on its fake window before loading the file. Never set on
  // tools.unabated.com, so production exposes nothing.
  if (window.__unabatedTicketExposeInternals === true) {
    window.__unabatedTicketInternals = { buildTicket, startWatching, readWatchedLine };
  }
  // Heartbeat so the panel can show whether this script is alive on the tab,
  // plus the books/bet-type filter for the Edges tab (grid may not be up yet
  // on the first tick; the error is published and the next tick retries).
  const HEARTBEAT_MS = 10000;
  function heartbeat() {
    if (retired) return;
    post("ready", { url: window.location.href, at: Date.now() });
    publishFilters();
  }
  post("takeover", { at: Date.now() });
  post("resume_request", { at: Date.now() });
  heartbeat();
  intervals.push(setInterval(heartbeat, HEARTBEAT_MS));
  console.info("[unabated-ticket] page.js active on", window.location.href);
})();
