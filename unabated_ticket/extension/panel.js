// Unabated Ticket — side panel: the Ticket tab (one captured bet) and the
// Edges tab (every positive-edge line across the enabled leagues).
//
// Reads: chrome.storage.local {ticket, error, watchStatus, pageReady,
// booksFilter, locateResult} (written by content.js) and {bankroll,
// multiplier, edges, activeTab} (settings, written here).
// Writes: chrome.storage.local settings and {locate} (row click, via
// locate.js, which also focuses the Unabated tab). Re-renders on storage.onChanged.
// Network: the scanner (scanner.js) fetches Unabated's public feeds while this
// page is open; it pauses when the panel is hidden and stops when it closes.

(function () {
  "use strict";

  const DEFAULT_SETTINGS = { bankroll: 30000, multiplier: 0.25 };
  const DEFAULT_EDGE_SETTINGS = { leagues: [1, 2, 5], periods: [1], minEdgePct: 1.0, sortBy: "edge" };
  // Watcher heartbeats every 5s; past this with no heartbeat, the Unabated tab is gone.
  const WATCH_STALE_MS = 15000;
  // page.js republishes the books filter every 10s while an Unabated tab is open.
  const BOOKS_FILTER_STALE_MS = 6 * 60 * 60 * 1000;
  const MAX_EDGE_ROWS = 200;
  const kelly = globalThis.UnabatedKelly;
  const feed = globalThis.UnabatedFeed;

  const el = (id) => document.getElementById(id);
  const view = {
    ticket: el("ticket"), error: el("error"), empty: el("empty"),
    warning: el("warning"), sideLabel: el("side-label"), betLine: el("bet-line"),
    eventLine: el("event-line"), startLine: el("start-line"),
    book: el("book"), price: el("price"), fair: el("fair"), edge: el("edge"),
    stake: el("stake"), fullKelly: el("full-kelly"),
    copy: el("copy"), copyStatus: el("copy-status"),
    errorTitle: el("error-title"), errorDetail: el("error-detail"), errorHint: el("error-hint"),
    bankroll: el("bankroll"), multiplier: el("multiplier"), settingsError: el("settings-error"),
    pageStatus: el("page-status"),
    tabs: el("tabs"), tabTicket: el("tab-ticket"), tabEdges: el("tab-edges"), edgesCount: el("edges-count"),
    edgesError: el("edges-error"), edgesStatus: el("edges-status"), edgesFilter: el("edges-filter"), edgesLocate: el("edges-locate"),
    edgesLeagues: el("edges-leagues"), edgesPeriods: el("edges-periods"), edgesMin: el("edges-min"), edgesSort: el("edges-sort"),
    edgesSettingsError: el("edges-settings-error"), edgesList: el("edges-list"), edgesEmpty: el("edges-empty"),
  };
  // page.js heartbeats every 10s; past this it is not running on any Unabated tab.
  const PAGE_READY_STALE_MS = 25000;

  let state = {
    ticket: null, error: null, watchStatus: null, pageReady: null, settings: { ...DEFAULT_SETTINGS },
    edgeSettings: { ...DEFAULT_EDGE_SETTINGS }, booksFilter: null, activeTab: "ticket",
    locateResult: null, locating: null,
  };
  let lastCopyText = "";
  let scannerStatus = null;
  let scannerState = null;
  const scanner = globalThis.UnabatedScanner.createScanner({
    onChange: (status, feedState) => {
      scannerStatus = status;
      scannerState = feedState;
      renderEdges();
    },
  });

  // ---- formatting ----------------------------------------------------------

  function fmtAmerican(price) {
    return price > 0 ? `+${price}` : `${price}`;
  }

  // Prediction-market style: implied probability in cents (what Kalshi/Novig
  // show). Uses the exchange's exact source price when the line carries one, so
  // it matches Unabated's screen instead of a rounded American round-trip.
  function fmtCents(line) {
    return `${(kelly.bookProbOf(line) * 100).toFixed(1)}\u00a2`;
  }

  function fmtPriceBoth(line) {
    return `${fmtAmerican(line.bookPrice)} \u00b7 ${fmtCents(line)}`;
  }

  function asBookLine(price, sourceFormat, sourcePrice) {
    return { bookPrice: price, sourceFormat, sourcePrice };
  }

  function fmtDollars(value) {
    return value.toLocaleString("en-US", { style: "currency", currency: "USD", minimumFractionDigits: 2, maximumFractionDigits: 2 });
  }

  function fmtPct(fraction) {
    const pct = fraction * 100;
    return `${pct > 0 ? "+" : ""}${pct.toFixed(2)}%`;
  }

  function fmtStart(eventStart) {
    if (!eventStart) return "";
    // Unabated's eventStart is naive UTC ("2026-09-12T23:30:00").
    const date = new Date(eventStart.endsWith("Z") ? eventStart : `${eventStart}Z`);
    if (Number.isNaN(date.getTime())) return eventStart;
    return date.toLocaleString([], { weekday: "short", month: "short", day: "numeric", hour: "numeric", minute: "2-digit" });
  }

  function fmtPoints(points) {
    return points == null ? "" : `${points > 0 ? "+" : ""}${points}`;
  }

  // ---- pricing -------------------------------------------------------------

  // The line the stake is computed from: the current line if it moved, else the captured one.
  function pricedLine(ticket) {
    const current = ticket.current;
    if (!current) return { price: ticket.price, sourceFormat: ticket.sourceFormat, sourcePrice: ticket.sourcePrice, fair: ticket.fair, edgePct: ticket.edgePct, points: ticket.points, moved: false };
    return { price: current.price, sourceFormat: current.sourceFormat, sourcePrice: current.sourcePrice, fair: current.fair, edgePct: current.edgePct, points: current.points, moved: true };
  }

  // Stake from Unabated's own edge for the line being priced (captured, or current if it moved).
  function computeStake(ticket, settings) {
    const line = pricedLine(ticket);
    if (line.edgePct == null) return { line, result: null, reason: "Unabated has no edge at the new line" };
    try {
      const result = kelly.kellyStakeFromEdge({ bookPrice: line.price, edgePct: line.edgePct, bankroll: settings.bankroll, multiplier: settings.multiplier });
      return { line, result, reason: null };
    } catch (error) {
      return { line, result: null, reason: error.message };
    }
  }

  // ---- side wording --------------------------------------------------------

  // "Total · Over 55.5 combined points" / "Spread · Oregon (away) vs Oklahoma State"
  function describeSide(ticket) {
    const rotation = ticket.rotation != null ? ` \u00b7 rot ${ticket.rotation}` : "";
    if (ticket.betType === "Total") {
      const overUnder = ticket.sideIndex === 0 ? "Over" : "Under";
      return `Total \u00b7 ${overUnder} ${ticket.points} combined points${rotation}`;
    }
    const opponent = ticket.sideIndex === 0 ? ticket.homeTeam : ticket.awayTeam;
    const vs = opponent ? ` vs ${opponent}` : "";
    const where = ticket.homeAway ? ` (${ticket.homeAway.toLowerCase()})` : "";
    return `${ticket.betType} \u00b7 ${ticket.homeAway === "Away" ? ticket.awayTeam || "" : ticket.homeTeam || ""}${where}${vs}${rotation}`;
  }

  // "Villanova Wildcats @ Louisville Cardinals · CFB", falling back to Unabated's event name.
  function describeMatchup(ticket) {
    const league = (ticket.league || "").toUpperCase();
    if (ticket.awayTeam && ticket.homeTeam) return `${ticket.awayTeam} @ ${ticket.homeTeam}${league ? ` \u00b7 ${league}` : ""}`;
    return ticket.eventName || "";
  }

  // ---- rendering -----------------------------------------------------------

  function show(which) {
    view.ticket.hidden = which !== "ticket";
    view.error.hidden = which !== "error";
    view.empty.hidden = which !== "empty";
  }

  function watcherIsLive(ticket, watchStatus) {
    if (!watchStatus || watchStatus.capturedAt !== ticket.capturedAt) {
      // No heartbeat yet: live only during the first interval after capture.
      return Date.now() - ticket.capturedAt < WATCH_STALE_MS;
    }
    return !watchStatus.error && Date.now() - watchStatus.seenAt < WATCH_STALE_MS;
  }

  function renderWarning(ticket, watchStatus) {
    const messages = [];
    let bad = false;
    if (ticket.current && ticket.current.offBoard) {
      messages.push("Off the board at this book.");
      bad = true;
    } else if (ticket.current) {
      const pts = ticket.current.points != null ? ` at ${fmtPoints(ticket.current.points)}` : "";
      messages.push(`Line moved: now ${fmtAmerican(ticket.current.price)}${pts} (captured ${fmtAmerican(ticket.price)}${ticket.points != null ? ` at ${fmtPoints(ticket.points)}` : ""}). Stake re-sized.`);
    }
    if (!watcherIsLive(ticket, watchStatus)) {
      const why = watchStatus && watchStatus.error ? `: ${watchStatus.error}` : "";
      messages.push(`Not watching the line${why}. Showing the captured price.`);
    }
    view.warning.hidden = messages.length === 0;
    view.warning.textContent = messages.join(" ");
    view.warning.classList.toggle("bad", bad);
  }

  function renderTicket() {
    const { ticket, settings, watchStatus } = state;
    renderWarning(ticket, watchStatus);

    view.sideLabel.textContent = ticket.sideLabel;
    view.betLine.textContent = describeSide(ticket);
    view.eventLine.textContent = describeMatchup(ticket);
    view.startLine.textContent = fmtStart(ticket.eventStart);

    const { line, result, reason } = computeStake(ticket, settings);
    view.book.textContent = ticket.book.name;
    view.price.textContent = fmtPriceBoth(asBookLine(line.price, line.sourceFormat, line.sourcePrice));
    // The fair is Unabated's own American number; there is no more exact source for it.
    view.fair.textContent = line.fair == null ? "unknown" : fmtPriceBoth(asBookLine(line.fair, 1, null));
    view.edge.textContent = line.edgePct == null ? "—" : fmtPct(line.edgePct / 100);

    view.stake.classList.remove("no-edge");
    if (!result) {
      view.stake.textContent = "—";
      view.stake.classList.add("no-edge");
      view.fullKelly.textContent = `Cannot size: ${reason}`;
    } else if (result.stake <= 0) {
      view.stake.textContent = fmtDollars(0);
      view.stake.classList.add("no-edge");
      view.fullKelly.textContent = "No edge at this price.";
    } else {
      view.stake.textContent = fmtDollars(result.stake);
      view.fullKelly.textContent = "";
    }

    const stakeText = result ? result.stake.toFixed(2) : "n/a";
    lastCopyText = `${ticket.sideLabel} ${fmtPriceBoth(asBookLine(line.price, line.sourceFormat, line.sourcePrice))} @ ${ticket.book.name} | fair ${line.fair == null ? "?" : fmtPriceBoth(asBookLine(line.fair, 1, null))} | edge ${line.edgePct == null ? "?" : fmtPct(line.edgePct / 100)} | stake $${stakeText} | ${describeMatchup(ticket)}`;
    view.copyStatus.textContent = "";
    show("ticket");
  }

  // Two kinds of capture error need opposite advice: no_fair is Unabated
  // having no number for that line (normal); read_failed means the page changed.
  const ERROR_COPY = {
    no_fair: {
      title: "No Unabated fair for this line",
      hint: "Unabated has not priced this line, so there is nothing to size against. This is normal for lopsided moneylines and exchange-only lines. Pick a line that shows an edge %.",
    },
    read_failed: {
      title: "Could not read this cell",
      hint: "Click the price again. If it keeps failing, Unabated's page changed; see README troubleshooting.",
    },
  };

  function render() {
    if (state.error) {
      const copy = ERROR_COPY[state.error.kind] || ERROR_COPY.read_failed;
      view.errorTitle.textContent = copy.title;
      view.errorHint.textContent = copy.hint;
      view.errorDetail.textContent = `${state.error.message} (${new Date(state.error.at).toLocaleTimeString()})`;
      show("error");
      return;
    }
    if (!state.ticket) {
      const ready = state.pageReady;
      const alive = ready && Date.now() - ready.at < PAGE_READY_STALE_MS;
      view.pageStatus.textContent = alive
        ? `Capture script active on ${ready.url}`
        : "Capture script not detected. Reload the Unabated tab; if this persists, see README troubleshooting.";
      show("empty");
      return;
    }
    renderTicket();
  }


  // ---- Edges tab -----------------------------------------------------------

  function fmtAge(ms) {
    if (ms < 1000) return "just now";
    if (ms < 60 * 1000) return `${Math.round(ms / 1000)}s ago`;
    if (ms < 60 * 60 * 1000) return `${Math.round(ms / 60000)}m ago`;
    return `${Math.round(ms / 3600000)}h ago`;
  }

  function fmtUntil(startMs) {
    const mins = Math.round((startMs - Date.now()) / 60000);
    if (mins < 60) return `in ${mins}m`;
    if (mins < 48 * 60) return `in ${Math.floor(mins / 60)}h ${mins % 60}m`;
    return `in ${Math.round(mins / 1440)}d`;
  }

  function fmtLiquidity(value) {
    return value == null ? "" : `liq ${value.toLocaleString("en-US", { style: "currency", currency: "USD", maximumFractionDigits: 0 })}`;
  }

  // Which books and bet types the list is restricted to. The user's own
  // Unabated selection wins when page.js has published one; otherwise every
  // live book and moneyline/spread/total, and the header says so.
  function effectiveFilter() {
    const filter = state.booksFilter;
    const fresh = filter && typeof filter.at === "number" && Date.now() - filter.at < BOOKS_FILTER_STALE_MS;
    const bookIds = fresh && Array.isArray(filter.bookIds) && filter.bookIds.length ? new Set(filter.bookIds) : null;
    const betTypeIds = fresh && Array.isArray(filter.betTypeIds) && filter.betTypeIds.length
      ? new Set(filter.betTypeIds.filter((id) => feed.BET_TYPES[id]))
      : null;
    return { bookIds, betTypeIds: betTypeIds && betTypeIds.size ? betTypeIds : null, fresh: Boolean(fresh), filter };
  }

  function describeFilter(effective) {
    const filter = effective.filter;
    const parts = [];
    if (effective.bookIds) {
      parts.push(`your ${effective.bookIds.size} Unabated books (read ${fmtAge(Date.now() - filter.at)})`);
    } else if (filter && filter.lastError) {
      parts.push(`no books filter yet: showing all live books (page read failed: ${filter.lastError})`);
    } else {
      parts.push("no books filter yet: showing all live books (open an Unabated odds tab to publish your selection)");
    }
    if (effective.betTypeIds) {
      parts.push(`bet types: ${Array.from(effective.betTypeIds).map((id) => feed.BET_TYPES[id]).join("/")}`);
    } else if (effective.fresh && filter.betTypeReason) {
      parts.push(`bet-type filter unreadable (${filter.betTypeReason}); showing ML/spread/total`);
    } else {
      parts.push("ML/spread/total");
    }
    if (effective.fresh && filter.lastError && effective.bookIds) parts.push(`latest page read failed: ${filter.lastError}`);
    return parts.join(" · ");
  }

  function stakeFor(row) {
    if (row.edgePct == null) return null;
    try {
      return kelly.kellyStakeFromEdge({ bookPrice: row.price, edgePct: row.edgePct, bankroll: state.settings.bankroll, multiplier: state.settings.multiplier }).stake;
    } catch (_error) {
      return null;
    }
  }

  function currentEdgeRows() {
    if (!scannerState) return [];
    const effective = effectiveFilter();
    const settings = state.edgeSettings;
    const rows = feed.selectEdges(scannerState, {
      minEdge: settings.minEdgePct / 100,
      periods: new Set(settings.periods),
      betTypes: effective.betTypeIds || new Set([1, 2, 3]),
      bookIds: effective.bookIds,
      now: Date.now(),
    }).map((row) => ({ ...row, stake: stakeFor(row) }));
    if (settings.sortBy === "stake") rows.sort((a, b) => (b.stake ?? -1) - (a.stake ?? -1) || b.edgePct - a.edgePct);
    if (settings.sortBy === "start") rows.sort((a, b) => a.eventStartMs - b.eventStartMs || b.edgePct - a.edgePct);
    return rows;
  }

  function renderEdgeRow(row) {
    const li = document.createElement("li");
    li.className = `edge-row${row.isBlurred ? " blurred" : ""}`;
    li.dataset.key = row.key;
    const top = document.createElement("div");
    top.className = "edge-top";
    const side = document.createElement("span");
    side.className = "edge-side";
    side.textContent = row.sideLabel;
    const pct = document.createElement("span");
    pct.className = "edge-pct";
    pct.textContent = fmtPct(row.edgePct / 100);
    top.append(side, pct);

    const bet = document.createElement("div");
    bet.className = "muted";
    bet.textContent = `${describeSide(row)}${row.period === "FG" ? "" : ` · ${row.period}`}`;
    const matchup = document.createElement("div");
    matchup.className = "muted";
    matchup.textContent = `${describeMatchup(row)} · ${fmtStart(row.eventStart)} · ${fmtUntil(row.eventStartMs)}`;

    const bottom = document.createElement("div");
    bottom.className = "edge-bottom";
    const book = document.createElement("span");
    book.className = "edge-book";
    book.textContent = `${row.book.name} ${fmtPriceBoth(asBookLine(row.price, row.sourceFormat, row.sourcePrice))}`;
    const liquidity = document.createElement("span");
    liquidity.className = "muted";
    liquidity.textContent = fmtLiquidity(row.liquidity);
    const stake = document.createElement("span");
    stake.className = "edge-stake";
    stake.textContent = row.stake == null ? "—" : fmtDollars(row.stake);
    bottom.append(book, liquidity, stake);

    li.append(top, bet, matchup, bottom);
    return li;
  }

  function renderEdgesStatus(rows) {
    const status = scannerStatus;
    if (!status) {
      view.edgesStatus.textContent = "Starting the scanner…";
      return;
    }
    const leagues = status.leaguesLoaded.map((id) => (feed.LEAGUES[id] || { label: `league ${id}` }).label);
    const updated = status.lastUpdateAt ? `updated ${fmtAge(Date.now() - status.lastUpdateAt)}` : (status.lastSnapshotAt ? `snapshot ${fmtAge(Date.now() - status.lastSnapshotAt)}` : "no data yet");
    const polled = status.lastPollAt ? ` · polled ${fmtAge(Date.now() - status.lastPollAt)}` : "";
    view.edgesStatus.textContent = status.phase === "loading"
      ? "Loading snapshots…"
      : `${leagues.join(" · ") || "no leagues"} · ${status.lineCount.toLocaleString()} lines · ${updated}${polled}`;
    view.edgesError.hidden = !status.error;
    view.edgesError.textContent = status.error || "";
    view.edgesCount.hidden = rows.length === 0;
    view.edgesCount.textContent = String(rows.length);
  }

  function renderEdges() {
    const rows = currentEdgeRows();
    renderEdgesStatus(rows);
    view.edgesFilter.textContent = describeFilter(effectiveFilter());
    view.edgesList.replaceChildren(...rows.slice(0, MAX_EDGE_ROWS).map(renderEdgeRow));
    const status = scannerStatus;
    if (rows.length === 0) {
      view.edgesEmpty.hidden = false;
      view.edgesEmpty.textContent = !status || status.phase !== "live"
        ? (status && status.phase === "error" ? "Nothing to list: the feed is unavailable (see above)." : "Waiting for the first snapshot…")
        : `No line at or above ${state.edgeSettings.minEdgePct}% edge right now.`;
    } else {
      view.edgesEmpty.hidden = rows.length > MAX_EDGE_ROWS ? false : true;
      view.edgesEmpty.textContent = rows.length > MAX_EDGE_ROWS ? `Showing the top ${MAX_EDGE_ROWS} of ${rows.length}; raise the minimum edge to see fewer.` : "";
    }
  }


  // ---- row click -> locate on the Unabated tab -----------------------------

  function locateRequestOf(row) {
    return {
      key: row.key, league: row.league, leagueLabel: row.leagueLabel, eventId: row.eventId,
      betTypeId: row.betTypeId, periodTypeId: row.periodTypeId, sideKey: row.sideKey, sideIndex: row.sideIndex,
      bookId: row.book.id, bookName: row.book.name, marketId: row.marketId, points: row.points, price: row.price,
      sideLabel: row.sideLabel, matchup: describeMatchup(row),
    };
  }

  function renderLocate() {
    const result = state.locateResult;
    const pending = state.locating;
    if (pending && (!result || result.at < pending.at)) {
      view.edgesLocate.hidden = false;
      view.edgesLocate.textContent = `Locating ${pending.sideLabel} @ ${pending.bookName} on the ${pending.leagueLabel} tab…`;
      return;
    }
    if (result && Date.now() - result.at < 60000) {
      view.edgesLocate.hidden = false;
      view.edgesLocate.textContent = result.ok
        ? `On the ${result.leagueLabel} tab: ${result.sideLabel} @ ${result.bookName} is highlighted. Click the price there to bet.`
        : `Could not show ${result.sideLabel} @ ${result.bookName}: ${result.message}`;
      return;
    }
    view.edgesLocate.hidden = true;
  }

  view.edgesList.addEventListener("click", async (event) => {
    const li = event.target.closest("li.edge-row");
    if (!li) return;
    const row = currentEdgeRows().find((r) => r.key === li.dataset.key);
    if (!row) return;
    const request = locateRequestOf(row);
    state.locating = { ...request, at: Date.now() };
    state.locateResult = null;
    renderLocate();
    try {
      await globalThis.UnabatedLocate.locateLine(request);
    } catch (error) {
      state.locating = null;
      state.locateResult = { ...request, at: Date.now(), ok: false, message: `could not focus an Unabated tab (${error.message})` };
      renderLocate();
    }
  });

  // ---- tabs ----------------------------------------------------------------

  function showTab(name) {
    state.activeTab = name === "edges" ? "edges" : "ticket";
    view.tabTicket.hidden = state.activeTab !== "ticket";
    view.tabEdges.hidden = state.activeTab !== "edges";
    for (const button of view.tabs.querySelectorAll("button[data-tab]")) {
      button.classList.toggle("active", button.dataset.tab === state.activeTab);
    }
  }

  view.tabs.addEventListener("click", (event) => {
    const button = event.target.closest("button[data-tab]");
    if (!button) return;
    showTab(button.dataset.tab);
    chrome.storage.local.set({ activeTab: state.activeTab });
    if (state.activeTab === "edges") renderEdges();
  });

  // ---- edge settings -------------------------------------------------------

  function readEdgeSettingInputs() {
    const leagues = Array.from(view.edgesLeagues.querySelectorAll("input:checked")).map((input) => Number(input.dataset.league));
    const periods = Array.from(view.edgesPeriods.querySelectorAll("input:checked")).map((input) => Number(input.dataset.period));
    const minEdgePct = Number(view.edgesMin.value);
    if (!Number.isFinite(minEdgePct) || minEdgePct < 0) return { error: "Minimum edge must be zero or more." };
    if (!periods.length) return { error: "Pick at least one period." };
    return { settings: { leagues, periods, minEdgePct, sortBy: view.edgesSort.value } };
  }

  function fillEdgeSettingInputs() {
    const settings = state.edgeSettings;
    for (const input of view.edgesLeagues.querySelectorAll("input")) input.checked = settings.leagues.includes(Number(input.dataset.league));
    for (const input of view.edgesPeriods.querySelectorAll("input")) input.checked = settings.periods.includes(Number(input.dataset.period));
    view.edgesMin.value = settings.minEdgePct;
    view.edgesSort.value = settings.sortBy;
  }

  function onEdgeSettingsInput() {
    const parsed = readEdgeSettingInputs();
    view.edgesSettingsError.textContent = parsed.error || "";
    if (parsed.error) return;
    const leaguesChanged = parsed.settings.leagues.join(",") !== state.edgeSettings.leagues.join(",");
    state.edgeSettings = parsed.settings;
    chrome.storage.local.set({ edges: parsed.settings });
    if (leaguesChanged) scanner.start(parsed.settings.leagues).catch((error) => console.error("[unabated-ticket] scanner restart failed", error));
    renderEdges();
  }

  function sanitizeEdgeSettings(stored) {
    const base = { ...DEFAULT_EDGE_SETTINGS };
    if (!stored || typeof stored !== "object") return base;
    if (Array.isArray(stored.leagues)) base.leagues = stored.leagues.filter((id) => feed.LEAGUES[id]);
    if (Array.isArray(stored.periods) && stored.periods.length) base.periods = stored.periods.filter((id) => feed.PERIODS[id]);
    if (typeof stored.minEdgePct === "number" && stored.minEdgePct >= 0) base.minEdgePct = stored.minEdgePct;
    if (["edge", "stake", "start"].includes(stored.sortBy)) base.sortBy = stored.sortBy;
    return base;
  }

  // ---- settings ------------------------------------------------------------

  function readSettingInputs() {
    const bankroll = Number(view.bankroll.value);
    const multiplier = Number(view.multiplier.value);
    if (!Number.isFinite(bankroll) || bankroll <= 0) return { error: "Bankroll must be a positive number." };
    if (!Number.isFinite(multiplier) || multiplier <= 0 || multiplier > 1) return { error: "Multiplier must be between 0 and 1." };
    return { settings: { bankroll, multiplier } };
  }

  function onSettingsInput() {
    const parsed = readSettingInputs();
    view.settingsError.textContent = parsed.error || "";
    if (parsed.error) return;
    state.settings = parsed.settings;
    chrome.storage.local.set(parsed.settings);
    render();
    renderEdges();
  }

  function fillSettingInputs() {
    view.bankroll.value = state.settings.bankroll;
    view.multiplier.value = state.settings.multiplier;
  }

  // ---- wiring --------------------------------------------------------------

  async function load() {
    const local = await chrome.storage.local.get(DEFAULT_SETTINGS);
    state.settings = { bankroll: Number(local.bankroll) || DEFAULT_SETTINGS.bankroll, multiplier: Number(local.multiplier) || DEFAULT_SETTINGS.multiplier };
    fillSettingInputs();
    const relay = await chrome.storage.local.get(["ticket", "error", "watchStatus", "pageReady", "booksFilter", "edges", "activeTab", "locateResult"]);
    state.ticket = relay.ticket || null;
    state.error = relay.error || null;
    state.watchStatus = relay.watchStatus || null;
    state.pageReady = relay.pageReady || null;
    state.booksFilter = relay.booksFilter || null;
    state.locateResult = relay.locateResult || null;
    state.edgeSettings = sanitizeEdgeSettings(relay.edges);
    fillEdgeSettingInputs();
    showTab(relay.activeTab === "edges" ? "edges" : "ticket");
    render();
    renderEdges();
    renderLocate();
    await scanner.start(state.edgeSettings.leagues);
  }

  chrome.storage.onChanged.addListener((changes, area) => {
    if (area !== "local") return;
    if ("ticket" in changes) {
      const previous = state.ticket;
      state.ticket = changes.ticket.newValue || null;
      // A fresh capture (not the watcher rewriting `current`) brings the Ticket tab forward.
      if (state.ticket && (!previous || previous.capturedAt !== state.ticket.capturedAt)) showTab("ticket");
    }
    if ("error" in changes) {
      state.error = changes.error.newValue || null;
      if (state.error) showTab("ticket");
    }
    if ("watchStatus" in changes) state.watchStatus = changes.watchStatus.newValue || null;
    if ("pageReady" in changes) state.pageReady = changes.pageReady.newValue || null;
    if ("ticket" in changes || "error" in changes || "watchStatus" in changes || "pageReady" in changes) render();
    if ("booksFilter" in changes) {
      state.booksFilter = changes.booksFilter.newValue || null;
      renderEdges();
    }
    if ("locateResult" in changes) {
      state.locateResult = changes.locateResult.newValue || null;
      if (state.locateResult && state.locating && state.locateResult.at >= state.locating.at) state.locating = null;
      renderLocate();
    }
  });

  view.edgesLeagues.addEventListener("change", onEdgeSettingsInput);
  view.edgesPeriods.addEventListener("change", onEdgeSettingsInput);
  view.edgesMin.addEventListener("input", onEdgeSettingsInput);
  view.edgesSort.addEventListener("change", onEdgeSettingsInput);

  // Nothing polls while the panel is hidden; back in view, the scanner catches up or resyncs.
  document.addEventListener("visibilitychange", () => {
    if (document.hidden) scanner.pause();
    else scanner.resume().catch((error) => console.error("[unabated-ticket] scanner resume failed", error));
  });

  view.bankroll.addEventListener("input", onSettingsInput);
  view.multiplier.addEventListener("input", onSettingsInput);

  view.copy.addEventListener("click", async () => {
    try {
      await navigator.clipboard.writeText(lastCopyText);
      view.copyStatus.textContent = "Copied";
    } catch (error) {
      view.copyStatus.textContent = `Copy failed: ${error.message}`;
    }
  });

  // Re-evaluate the "not watching" state and the edge ages even when no event arrives.
  setInterval(() => {
    if (!state.error) render();
    if (state.activeTab === "edges") {
      renderEdges();
      renderLocate();
    }
  }, 5000);

  load().catch((error) => {
    view.errorDetail.textContent = error.message;
    show("error");
  });
})();
