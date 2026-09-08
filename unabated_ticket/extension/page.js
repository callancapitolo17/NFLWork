// Unabated Ticket — page-world script (runs in the page's own JS world).
//
// Why MAIN world: the ticket is read from React fiber props and the AG Grid
// row data hanging off the odds-cell DOM nodes. Those expando properties are
// invisible from an isolated-world content script, so this file runs in the
// page world and hands results to content.js via window.postMessage.
//
// Side effects: none on the page. One capture-phase click listener on
// document (Unabated's own handler still runs), one 5s interval while a
// ticket is being watched. Never touches the DOM.

(function () {
  "use strict";

  const MESSAGE_SOURCE = "unabated-ticket";
  const CELL_SHELL_SELECTOR = ".odds-cell-action-shell";
  // The "..." menu button is a child of the shell; opening a menu is not picking a bet.
  const MORE_BUTTON_SELECTOR = ".odds-cell-more-button";
  const MAX_FIBER_HOPS = 40;
  const WATCH_INTERVAL_MS = 5000;
  const BET_TYPE_NAMES = { 1: "Moneyline", 2: "Spread", 3: "Total" };

  // ---- messaging -----------------------------------------------------------

  function post(type, payload) {
    window.postMessage({ source: MESSAGE_SOURCE, type, payload }, window.location.origin);
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

  function fairPriceOf(marketLine) {
    // bacr = Unabated fair American price at this book's points. Null is a
    // normal condition (lopsided moneylines, exchange-only lines), not a parse
    // failure, so it gets its own error kind for the panel.
    const fair = marketLine.bacr;
    if (typeof fair !== "number" || !Number.isFinite(fair)) {
      const error = new Error("Unabated has no fair price for this line");
      error.kind = "no_fair";
      throw error;
    }
    return fair;
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

  function edgePctOf(marketLine) {
    const edge = marketLine.edge && marketLine.edge.edge;
    return typeof edge === "number" && Number.isFinite(edge) ? edge : null;
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

  // Which "ms<id>" entry under this side holds `marketLine` (same object, or same line id).
  function bookIdFromSides(marketLine, rowData, sideKey) {
    const books = rowData.sides && rowData.sides[sideKey];
    if (!books) return null;
    const wantedId = lineIdOf(marketLine);
    for (const [bookKey, line] of Object.entries(books)) {
      const sameObject = line === marketLine;
      const sameId = wantedId != null && line && String(lineIdOf(line)) === String(wantedId);
      if (sameObject || sameId) {
        const parsed = Number(bookKey.replace(/^ms/, ""));
        if (Number.isInteger(parsed)) return parsed;
      }
    }
    return null;
  }

  function bookIdOf(marketLine, cellProps, rowData, sideKey) {
    if (typeof marketLine.marketSourceId === "number") return marketLine.marketSourceId;
    if (cellProps && cellProps.marketSource && typeof cellProps.marketSource.id === "number") {
      return cellProps.marketSource.id;
    }
    // Best-line cells sit in a column with no book id; find the line inside the row's sides instead.
    const fromSides = bookIdFromSides(marketLine, rowData, sideKey);
    if (fromSides != null) return fromSides;
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

  function buildTicket({ marketLine, sideIndex, rowData, context, cellProps }) {
    const betTypeId = rowData.betTypeId;
    const betType = BET_TYPE_NAMES[betTypeId];
    if (!betType) throw new Error(`unsupported betTypeId ${betTypeId} (only moneyline, spread, total)`);

    const points = marketLine.points ?? null;
    const sideKey = sideKeyOf(rowData, sideIndex);
    const bookId = bookIdOf(marketLine, cellProps, rowData, sideKey);
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
      fair: fairPriceOf(marketLine),
      edgePct: edgePctOf(marketLine),
      // Watcher handle: how to find this same line again through the grid API.
      watch: { gridKey: rowData.gridKey ?? null, sideKey, bookKey: `ms${bookId}` },
      current: null,
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

  function anyGridApi() {
    const cells = document.querySelectorAll(".ag-cell");
    for (const cell of cells) {
      const props = findProps(fiberOf(cell), isGridCellProps);
      if (props) return { api: props.api, context: props.context ?? null };
    }
    throw new Error("could not reach the AG Grid API from any rendered cell");
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

  function readWatchedLine() {
    const { ticket, gridApi } = watcher;
    const { gridKey, sideKey, bookKey } = ticket.watch;
    if (!gridKey) throw new Error("no gridKey on the ticket");
    if (typeof gridApi.isDestroyed === "function" && gridApi.isDestroyed()) {
      throw new Error("grid was destroyed");
    }
    const node = gridApi.getRowNode(gridKey);
    if (!node || !node.data) throw new Error("row no longer in the grid");
    const books = node.data.sides && node.data.sides[sideKey];
    const line = books && books[bookKey];
    if (!line) throw new Error("book line no longer on the row");
    return {
      price: bookPriceOf(line),
      ...sourcePriceOf(line),
      points: line.points ?? null,
      fair: typeof line.bacr === "number" ? line.bacr : null,
      offBoard: line.statusId === 2,
      seenAt: Date.now(),
    };
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
  }

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
    const target = event.target instanceof Element ? event.target : null;
    const shell = target && target.closest(CELL_SHELL_SELECTOR);
    if (!shell) return;
    if (target.closest(MORE_BUTTON_SELECTOR)) return;
    if (event.type === "pointerdown" && event.button !== PRIMARY_BUTTON && event.button !== MIDDLE_BUTTON) return;
    if (lastCapture.shell === shell && Date.now() - lastCapture.at < CAPTURE_DEDUPE_MS) return;
    lastCapture = { shell, at: Date.now() };
    try {
      const cell = readCell(shell);
      const ticket = buildTicket(cell);
      startWatching(ticket, cell.gridApi);
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
  // Heartbeat so the panel can show whether this script is alive on the tab.
  post("ready", { url: window.location.href, at: Date.now() });
  setInterval(() => post("ready", { url: window.location.href, at: Date.now() }), 10000);
  console.info("[unabated-ticket] page.js active on", window.location.href);
})();
