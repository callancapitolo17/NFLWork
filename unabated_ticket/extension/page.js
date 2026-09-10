// Unabated Ticket — page-world script (runs in the page's own JS world).
//
// Why MAIN world: the ticket is read from React fiber props and the AG Grid
// row data hanging off the odds-cell DOM nodes. Those expando properties are
// invisible from an isolated-world content script, so this file runs in the
// page world and hands results to content.js via window.postMessage.
//
// Side effects: one capture-phase click listener on document (Unabated's
// own handler still runs), one 5s interval while a ticket is being watched,
// and a 10s heartbeat that also publishes the user's Unabated book selection
// + bet-type filter (read from the grid's React context and localStorage) so
// the Edges tab can filter on them. The only DOM touch is the locate flash:
// a 2.5s outline on the cell an Edges row or notification pointed at.

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

  // bacr = Unabated fair American price at this book's points; shown for
  // information only, the stake is sized from Unabated's edge.
  function fairPriceOrNull(marketLine) {
    const fair = marketLine.bacr;
    return typeof fair === "number" && Number.isFinite(fair) ? fair : null;
  }

  // Unabated's own edge % (EV per $1 staked) for this line. Null is a normal
  // condition (lopsided moneylines, exchange-only lines Unabated has not
  // priced), not a parse failure, so it gets its own error kind for the panel.
  function requireEdgePct(marketLine) {
    const edge = edgePctOf(marketLine);
    if (edge == null) {
      const error = new Error("Unabated has no edge for this line");
      error.kind = "no_fair";
      throw error;
    }
    return edge;
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
      fair: fairPriceOrNull(marketLine),
      edgePct: requireEdgePct(marketLine),
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
      fair: fairPriceOrNull(line),
      edgePct: edgePctOf(line),
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


  // ---- books / bet-type filter publish ------------------------------------

  // The odds screen's book selection lives in context.userSettings.gameOdds:
  // one entry per market source, isUnavailable === false when the book is
  // shown. Accepts an array or an id-keyed object. Throws when unreadable so
  // the panel says "no books filter yet" rather than silently showing all.
  function enabledBookIdsOf(userSettings) {
    const gameOdds = userSettings && userSettings.gameOdds;
    if (!gameOdds || typeof gameOdds !== "object") throw new Error("userSettings.gameOdds missing");
    const entries = Array.isArray(gameOdds)
      ? gameOdds.map((entry) => ({ entry, key: null }))
      : Object.entries(gameOdds).map(([key, entry]) => ({ entry, key }));
    const ids = [];
    for (const { entry, key } of entries) {
      if (!entry || typeof entry !== "object" || entry.isUnavailable !== false) continue;
      const id = Number(entry.marketSourceId ?? entry.id ?? key);
      if (Number.isInteger(id)) ids.push(id);
    }
    if (!ids.length) throw new Error(`no enabled books in userSettings.gameOdds (${entries.length} entries)`);
    return ids;
  }

  const ODDS_FILTER_STORAGE_KEY = "oddsFilterContext:preferences";

  // Bet types selected in the odds-screen filter. Shape unverified (needs a
  // login to see), so this looks for any betType* array of ids or {id}
  // objects and reports null + reason when nothing matches.
  function selectedBetTypeIdsOf() {
    let raw = null;
    try { raw = window.localStorage.getItem(ODDS_FILTER_STORAGE_KEY); } catch (_error) { /* storage blocked */ }
    if (!raw) return { betTypeIds: null, reason: `localStorage ${ODDS_FILTER_STORAGE_KEY} missing` };
    let parsed;
    try { parsed = JSON.parse(raw); } catch (_error) { return { betTypeIds: null, reason: `${ODDS_FILTER_STORAGE_KEY} is not JSON` }; }
    const found = findBetTypeIds(parsed, 0);
    if (!found) return { betTypeIds: null, reason: `no betType ids found under ${ODDS_FILTER_STORAGE_KEY} (keys: ${Object.keys(parsed || {}).join(",")})` };
    return { betTypeIds: found, reason: null };
  }

  function idsFromArray(values) {
    const ids = [];
    for (const value of values) {
      const id = typeof value === "number" ? value : Number(value && (value.betTypeId ?? value.id ?? value.value));
      if (Number.isInteger(id)) ids.push(id);
    }
    return ids.length ? ids : null;
  }

  function findBetTypeIds(node, depth) {
    if (!node || typeof node !== "object" || depth > 4) return null;
    for (const [key, value] of Object.entries(node)) {
      if (/bettype/i.test(key) && Array.isArray(value)) {
        const ids = idsFromArray(value);
        if (ids) return ids;
      }
    }
    for (const value of Object.values(node)) {
      const ids = findBetTypeIds(value, depth + 1);
      if (ids) return ids;
    }
    return null;
  }

  let lastFiltersSignature = null;

  function publishFilters() {
    let payload;
    try {
      const { context } = anyGridApi();
      const bookIds = enabledBookIdsOf(context && context.userSettings);
      const betTypes = selectedBetTypeIdsOf();
      payload = { bookIds, betTypeIds: betTypes.betTypeIds, betTypeReason: betTypes.reason, error: null, url: window.location.href, at: Date.now() };
    } catch (error) {
      payload = { bookIds: null, betTypeIds: null, betTypeReason: null, error: error.message, url: window.location.href, at: Date.now() };
    }
    const signature = JSON.stringify([payload.bookIds, payload.betTypeIds, payload.error]);
    if (signature !== lastFiltersSignature) {
      lastFiltersSignature = signature;
      console.info("[unabated-ticket] books filter", payload);
    }
    post("filters", payload);
  }


  // ---- locate: scroll the grid to a line and flash its cell -----------------

  const LOCATE_ATTEMPTS = 20;
  const LOCATE_RETRY_MS = 1000;
  const FLASH_MS = 2500;
  const FLASH_ATTR = "data-unabated-ticket-flash";
  let lastLocateAt = 0;

  function findRowNode(api, request) {
    let hit = null;
    api.forEachNode((node) => {
      if (hit || !node.data) return;
      const data = node.data;
      if (data.eventId !== request.eventId || data.betTypeId !== request.betTypeId) return;
      if ((data.periodTypeId ?? 1) !== request.periodTypeId) return;
      hit = node;
    });
    return hit;
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
    try {
      if (typeof api.ensureNodeVisible === "function") api.ensureNodeVisible(node, "middle");
      if (typeof api.ensureColumnVisible === "function") api.ensureColumnVisible(String(request.bookId));
    } catch (error) {
      reportLocate(request, false, `grid scroll failed: ${error.message}`);
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

  window.addEventListener("message", (event) => {
    if (event.source !== window) return;
    const data = event.data;
    if (!data || data.source !== MESSAGE_SOURCE || data.type !== "locate") return;
    onLocateMessage(data.payload);
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
  // Heartbeat so the panel can show whether this script is alive on the tab,
  // plus the books/bet-type filter for the Edges tab (grid may not be up yet
  // on the first tick; the error is published and the next tick retries).
  const HEARTBEAT_MS = 10000;
  function heartbeat() {
    post("ready", { url: window.location.href, at: Date.now() });
    publishFilters();
  }
  heartbeat();
  setInterval(heartbeat, HEARTBEAT_MS);
  console.info("[unabated-ticket] page.js active on", window.location.href);
})();
