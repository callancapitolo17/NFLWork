// Unabated Ticket — side panel: the Ticket tab (one captured bet) and the
// Edges tab (every positive-edge line across the enabled leagues).
//
// Reads: chrome.storage.local {ticket, error, watchStatus, pageCheck, pageReady,
// booksFilter, locateResult} (written by content.js) and {bankroll,
// multiplier, edges, alerts, alertLog, activeTab, betsService, betsSettings}
// (written here) and {betsNovig} (written by novig_content.js on
// app.novig.us, #116).
// Writes: chrome.storage.local settings, {locate} (row click, via locate.js,
// which also focuses the Unabated tab), {alertLog, alertTargets} and Chrome
// notifications for new edges; {betsService} after every bets-service poll
// (records + the service's team crosswalk); POST /crosswalk.json to the
// bets service with the team rows an id join taught, DELETE it on Clear
// (#118 step 4). Re-renders on storage.onChanged.
// Network: the scanner (scanner.js) fetches Unabated's public feeds while this
// page is open; it pauses when the panel is hidden and stops when it closes.
// The bets tab (#114) polls the local bets service (betsSettings.serviceUrl,
// default http://127.0.0.1:8094) every 30 s on the same visibility rule —
// never from the service worker. Matching is bets.js; presentation helpers
// are betsview.js.

(function () {
  "use strict";

  const kelly = globalThis.UnabatedKelly;
  const feed = globalThis.UnabatedFeed;
  const betsLib = globalThis.UnabatedBets;
  const betsView = globalThis.UnabatedBetsView;
  const ladderLib = globalThis.UnabatedLadder;
  // Bets service poll cadence while the panel is visible (plan § Storage).
  const BETS_POLL_MS = 30 * 1000;
  const DEFAULT_SETTINGS = { bankroll: 30000, multiplier: 0.25 };
  // maxLineAgeHours: a "live" book's line unchanged for a week is a dead feed
  // (live 2026-09-10: Buckeye -110 on a 44.5 total, 96 days old, "+36.67%").
  const ALL_LEAGUE_IDS = Object.keys(feed.LEAGUES).map(Number);
  // bookIds undefined (never ticked) = DEFAULT_BOOK_NAMES; null = follow the
  // Unabated selection page.js publishes (all live books until one exists),
  // set by its button; an array = the user's own ticks in the panel.
  // Alt lines (#113) are off until asked for; altMaxDistance 7 points keeps
  // NFL/CFB spreads to about a touchdown off the number (live 2026-09-11 the
  // median NFL alt "edge" sat 13 points out, +400 and up). minLiquidity $100
  // (every exchange line, main or alt) is Kalshi's median alt depth ($129)
  // with its thin tail cut. 0 = off.
  const DEFAULT_EDGE_SETTINGS = {
    leagues: ALL_LEAGUE_IDS, periods: [1], betTypes: [1, 2, 3], bookIds: undefined, minEdgePct: 1.0, maxLineAgeHours: 168, sortBy: "edge",
    minStake: 0,
    minLiquidity: 100,
    includeAlts: false, altMaxDistance: 7,
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
  // Off until the list has been watched for a session (plan, 2026-09-10).
  const DEFAULT_ALERT_SETTINGS = { enabled: false, minEdgePct: 2.0 };
  const ALERT_EVENT_COOLDOWN_MS = 5 * 60 * 1000;
  const ALERT_LOG_TTL_MS = 24 * 60 * 60 * 1000;
  const ALERT_TARGETS_KEPT = 50;
  // Watcher heartbeats every 5s; past this with no heartbeat, the Unabated tab is gone.
  const WATCH_STALE_MS = 15000;
  // page.js republishes the books filter every 10s while an Unabated tab is open.
  const BOOKS_FILTER_STALE_MS = 6 * 60 * 60 * 1000;
  const MAX_EDGE_ROWS = 200;

  const el = (id) => document.getElementById(id);
  const view = {
    ticket: el("ticket"), error: el("error"), empty: el("empty"),
    warning: el("warning"), rowTrace: el("row-trace"), sideLabel: el("side-label"), betLine: el("bet-line"),
    eventLine: el("event-line"), startLine: el("start-line"),
    book: el("book"), price: el("price"), fair: el("fair"), edge: el("edge"),
    stake: el("stake"), contracts: el("contracts"), fullKelly: el("full-kelly"), stakeExposure: el("stake-exposure"), payoutRow: el("payout-row"), profit: el("profit"), payout: el("payout"),
    copy: el("copy"), copyStatus: el("copy-status"),
    errorTitle: el("error-title"), errorDetail: el("error-detail"), errorHint: el("error-hint"),
    bankroll: el("bankroll"), multiplier: el("multiplier"), settingsError: el("settings-error"),
    pageStatus: el("page-status"),
    rowTraceReason: el("row-trace-reason"), rowTraceDetail: el("row-trace-detail"),
    tabs: el("tabs"), tabTicket: el("tab-ticket"), tabEdges: el("tab-edges"), edgesCount: el("edges-count"),
    edgesToolbar: el("edges-toolbar"), edgesControls: el("edges-controls"),
    filtersToggle: el("filters-toggle"), filtersSummary: el("filters-summary"),
    settingsToggle: el("settings-toggle"), settings: el("settings"),
    backToEdges: el("back-to-edges"), stakeLabel: el("stake-label"), betsBannerHead: el("bets-banner-head"),
    betsCount: el("bets-count"), betsRisk: el("bets-risk"), betsRiskCaption: el("bets-risk-caption"),
    edgesError: el("edges-error"), edgesStatus: el("edges-status"), edgesFilter: el("edges-filter"), edgesFilterDebug: el("edges-filter-debug"), edgesLocate: el("edges-locate"),
    edgesSports: el("edges-sports"), edgesBetTypes: el("edges-bettypes"), edgesBooks: el("edges-books"), edgesBooksMode: el("edges-books-mode"),
    booksDefault: el("books-default"), booksUnabated: el("books-unabated"), booksAll: el("books-all"), booksNone: el("books-none"), edgesPeriods: el("edges-periods"), edgesMin: el("edges-min"), edgesMinStake: el("edges-min-stake"), edgesMaxAge: el("edges-max-age"), edgesSort: el("edges-sort"),
    edgesIncludeAlts: el("edges-include-alts"), edgesAltDistance: el("edges-alt-distance"), edgesMinLiquidity: el("edges-min-liquidity"), edgesGroup: el("edges-group"),
    edgesSettingsError: el("edges-settings-error"), edgesList: el("edges-list"), edgesEmpty: el("edges-empty"),
    alertsEnabled: el("alerts-enabled"), alertsMin: el("alerts-min"),
    betsHeader: el("bets-header"), betsBanner: el("bets-banner"),
    tabBets: el("tab-bets"), betsService: el("bets-service"), betsSources: el("bets-sources"),
    betsUrl: el("bets-url"), betsSettingsError: el("bets-settings-error"),
    betsOpen: el("bets-open"), betsOpenCount: el("bets-open-count"), betsOpenEmpty: el("bets-open-empty"),
    betsUnmatched: el("bets-unmatched"), betsUnmatchedCount: el("bets-unmatched-count"), betsUnmatchedEmpty: el("bets-unmatched-empty"),
    betsCrosswalk: el("bets-crosswalk"), betsCrosswalkCount: el("bets-crosswalk-count"), betsCrosswalkEmpty: el("bets-crosswalk-empty"),
    betsCrosswalkClear: el("bets-crosswalk-clear"),
    shapeBanner: el("shape-banner"),
  };
  // page.js heartbeats every 10s; past this it is not running on any Unabated tab.
  const PAGE_READY_STALE_MS = 25000;

  let state = {
    ticket: null, error: null, watchStatus: null, pageReady: null, settings: { ...DEFAULT_SETTINGS },
    // page.js's once-per-load shape check: {status: "checking"|"ok"|"changed", ...}.
    pageCheck: null,
    edgeSettings: { ...DEFAULT_EDGE_SETTINGS }, booksFilter: null, activeTab: "ticket",
    locateResult: null, locating: null,
    alertSettings: { ...DEFAULT_ALERT_SETTINGS },
    betsSettings: { ...betsView.DEFAULT_BETS_SETTINGS },
    // What the last bets-service poll left: {payload: {generatedAt, sources},
    // okAt, error, errorAt, unreachableSince}; null before the first poll.
    betsService: null,
    // What novig_content.js last wrote: {bets, readAt, url, error, complete, pageSeenAt} or null.
    betsNovig: null,
    // Normalised bet records (bets.js contract), team keys resolved, pruned to the retention window.
    betRecords: [],
    // The team crosswalk the bets service holds (#118 step 4), as last served
    // or stored: [{venue, league, venueTeamKey, venueTeamName, unabatedTeamId,
    // unabatedTeamName, learnedFrom, learnedAt}]. Keys resolve through it first.
    crosswalk: [],
  };
  // key -> {price, at}: what has been alerted (or seen at baseline); persisted.
  let alertLog = {};
  let eventAlertAt = {};
  // The first pass after a scanner (re)start records what is already on the
  // board without notifying, so opening the panel is not twenty pings.
  let alertsBaselined = false;
  // processAlerts awaits storage + notifications; a poll landing mid-run must
  // not start a second pass that notifies the same line twice.
  let alertsBusy = false;
  let lastCopyText = "";
  // The row a locate (or a capture) last came from, so returning to the Edges
  // list shows where you were rather than only restoring the scroll offset.
  let lastClickedKey = null;
  let scannerStatus = null;
  let scannerState = null;
  let boardLinesCache = null;
  // Unabated's fair ladders for sizing against held bets (#130): the feed's
  // lines grouped by event, and each (event, period, axis) ladder built from
  // them on first use. Both are dropped on every scanner update.
  let linesByEventCache = null;
  let ladderCache = new Map();
  const teamsLib = globalThis.UnabatedTeams;
  let teamsSpellingCount = 0;

  // Every snapshot carries Unabated's team list and, per game row, a second
  // spelling of each team (feed.teamSpellingsFromEventName, #118): register
  // both as the team index (teams.js), persist it, and fill keys on bet
  // records that were waiting for it (#116 — no hand-written team tables).
  // The trigger counts spellings, not teams: a new eventName spelling for a
  // team already indexed must persist and re-resolve too.
  function registerFeedTeams(feedState) {
    if (!feedState || !feedState.teamIndex) return;
    const byLeague = {};
    for (const team of Object.values(feedState.teamIndex)) {
      const league = feed.LEAGUES[team.leagueId];
      if (!league) continue;
      (byLeague[league.path] ||= []).push(team);
    }
    for (const [league, list] of Object.entries(byLeague)) teamsLib.registerTeams(league, list);
    const spellings = teamsLib.spellingCount();
    if (spellings === teamsSpellingCount) return;
    teamsSpellingCount = spellings;
    chrome.storage.local.set({ teamsIndex: teamsLib.exportIndex() });
    state.betRecords = betsLib.resolveTeamKeys(state.betRecords, state.crosswalk);
  }

  const scanner = globalThis.UnabatedScanner.createScanner({
    onChange: (status, feedState) => {
      scannerStatus = status;
      scannerState = feedState;
      boardLinesCache = null;
      linesByEventCache = null;
      ladderCache = new Map();
      registerFeedTeams(feedState);
      learnCrosswalk().catch((error) => console.error("[unabated-ticket] crosswalk learn failed", error));
      renderEdges();
      // A ticket sized from the feed (or waiting for it) follows the feed's
      // updates; one the screen priced is left alone (a re-render clears the copy status).
      if (state.ticket && !state.error && pricedLine(state.ticket).edgeFrom !== "screen") render();
      processAlerts().catch((error) => console.error("[unabated-ticket] alerts failed", error));
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

  // Does this feed line describe the ticket's line: same game, bet type,
  // period, side and book. Points are the caller's (the captured or current number).
  function feedLineMatchesTicket(line, ticket) {
    return line.eventId === ticket.eventId
      && line.betTypeId === ticket.watch.betTypeId
      && line.periodTypeId === ticket.watch.periodTypeId
      && line.sideIndex === ticket.sideIndex
      && line.bookId === ticket.book.id;
  }

  // Every feed line with the ticket's game, bet type, period, side, book and
  // points. Matched on the fields both sides carry, never on the feed key: an
  // alt rung's cell object can lack marketId, and whether page.js classed the
  // cell as an alt does not decide which key the feed filed the line under
  // (live 2026-09-12: an Under 46.5 +213 Novig rung the Edges tab listed came
  // back "no copy" through the key). Normally one line: the feed drops an alt
  // sitting on its main line's points. Two means two MARKETS — the changes
  // stream tags an event's team totals bt3 like its game total, and only the
  // marketId the ticket may lack tells them apart — so the caller refuses.
  function feedLinesFor(ticket, points) {
    if (!scannerState || !ticket.watch || ticket.eventId == null) return [];
    const matches = [];
    for (const line of Object.values(scannerState.lines)) {
      if (line.points === points && feedLineMatchesTicket(line, ticket)) matches.push(line);
    }
    return matches;
  }

  // The Edges feed's one copy of the ticket's line, or null. The screen cell
  // can carry no edge at all, and this is the same `ge` the Edges tab sizes from.
  function feedLineFor(ticket, points) {
    const matches = feedLinesFor(ticket, points);
    return matches.length === 1 ? matches[0] : null;
  }

  // Every number the feed holds for the ticket's game, side and book, for the
  // no-edge view: tells a game the feed lacks apart from a missing rung.
  function feedPointsHeld(ticket) {
    if (!scannerState || !ticket.watch || ticket.eventId == null) return [];
    const points = [];
    for (const line of Object.values(scannerState.lines)) {
      if (feedLineMatchesTicket(line, ticket)) points.push(line.points);
    }
    return points.sort((a, b) => a - b);
  }

  // The line the stake is computed from: the current line if it moved, else
  // the captured one. edgeFrom says where the edge came from: "screen" (the
  // cell or its row), "feed" (the Edges feed at the same price), or null
  // (none — feedLine then carries the feed's copy when there is one, so the
  // panel can say what price the feed has instead).
  function pricedLine(ticket) {
    const current = currentOf(ticket);
    const line = current
      ? { price: current.price, sourceFormat: current.sourceFormat, sourcePrice: current.sourcePrice, fair: current.fair, edgePct: current.edgePct, points: current.points, moved: true }
      : { price: ticket.price, sourceFormat: ticket.sourceFormat, sourcePrice: ticket.sourcePrice, fair: ticket.fair, edgePct: ticket.edgePct, points: ticket.points, moved: false };
    line.edgeFrom = line.edgePct == null ? null : "screen";
    line.feedLine = null;
    if (line.edgePct != null) return line;
    const held = feedLineFor(ticket, line.points);
    if (!held) return line;
    line.feedLine = held;
    // An edge is for one price: the feed's number only applies at the price the cell shows.
    if (held.price !== line.price || held.ge == null) return line;
    line.edgePct = Math.round(held.ge * 1e6) / 1e4;
    if (line.fair == null) line.fair = held.bacr;
    line.edgeFrom = "feed";
    return line;
  }

  function noEdgeReason(line) {
    const held = line.feedLine;
    if (!held) return "Unabated has no edge at this line";
    if (held.ge == null) return `Unabated has no edge at this line (the Edges feed has it at ${fmtAmerican(held.price)} with no edge either)`;
    return `Unabated has no edge at ${fmtAmerican(line.price)} (the Edges feed has this line at ${fmtAmerican(held.price)} with ${fmtPct(held.ge)}; the cell shows a different price)`;
  }

  // Stake from Unabated's own edge for the line being priced (captured, or current if it moved).
  function computeStake(ticket, settings) {
    const line = pricedLine(ticket);
    if (line.edgePct == null) return { line, result: null, reason: noEdgeReason(line) };
    try {
      const result = kelly.kellyStakeFromEdge({ bookPrice: line.price, edgePct: line.edgePct, bankroll: settings.bankroll, multiplier: settings.multiplier });
      return { line, result, reason: null };
    } catch (error) {
      return { line, result: null, reason: error.message };
    }
  }

  // ---- side wording --------------------------------------------------------

  // "Total · Over 55.5 combined points" / "Spread · Oregon (away) vs Oklahoma State"
  // " · 1H" on anything but the full game: a first-half ticket read as a
  // full-game one until 2026-09-12 (the bet banner said 1H, the heading did not).
  function periodSuffix(ticket) {
    return ticket.period && ticket.period !== "FG" ? ` \u00b7 ${ticket.period}` : "";
  }

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
    const league = ticket.leagueLabel || (ticket.league || "").toUpperCase();
    if (ticket.awayTeam && ticket.homeTeam) return `${ticket.awayTeam} @ ${ticket.homeTeam}${league ? ` \u00b7 ${league}` : ""}`;
    return ticket.eventName || "";
  }

  // ---- rendering -----------------------------------------------------------

  function show(which) {
    view.ticket.hidden = which !== "ticket";
    view.error.hidden = which !== "error";
    view.empty.hidden = which !== "empty";
  }

  // The watched line when it has moved off the captured one, else null.
  // It rides on `watchStatus`, not on the ticket: content.js's watcher writes
  // only that key, so a tick in flight can never put a stale ticket back over
  // a fresh capture. A late tick for the PREVIOUS capture is what capturedAt
  // guards here.
  function currentOf(ticket) {
    const status = state.watchStatus;
    if (!ticket || !status || status.capturedAt !== ticket.capturedAt) return null;
    return status.current || null;
  }

  function watcherIsLive(ticket, watchStatus) {
    if (!watchStatus || watchStatus.capturedAt !== ticket.capturedAt) {
      // No heartbeat yet: live only during the first interval after capture.
      return Date.now() - ticket.capturedAt < WATCH_STALE_MS;
    }
    return !watchStatus.error && Date.now() - watchStatus.seenAt < WATCH_STALE_MS;
  }

  function renderWarning(ticket, watchStatus, line) {
    const messages = [];
    let bad = false;
    const current = currentOf(ticket);
    if (current && current.offBoard) {
      messages.push("Off the board at this book.");
      bad = true;
    } else if (current) {
      const pts = current.points != null ? ` at ${fmtPoints(current.points)}` : "";
      messages.push(`Line moved: now ${fmtAmerican(current.price)}${pts} (captured ${fmtAmerican(ticket.price)}${ticket.points != null ? ` at ${fmtPoints(ticket.points)}` : ""}). Stake re-sized.`);
    }
    if (line.edgeFrom === "feed") {
      const age = fmtLineAge(feed.lineChangedMs(line.feedLine)).replace(/^line /, "");
      messages.push(`Edge from the Edges feed (the screen cell carried none): same line at the same price, feed copy ${age}.`);
    }
    if (betsView.sourcesUnavailable(state.betsService && state.betsService.payload, Date.now(), pageSources())) {
      messages.push("Bet sources unavailable (no venue has reported in the last hour), so bet flags may be missing; see the Bets tab.");
    }
    if (!pageScriptAlive()) {
      messages.push("No Unabated tab is running the capture script, so new clicks will not reach this panel. Open an odds tab, or reload the one you have.");
      bad = true;
    } else if (!watcherIsLive(ticket, watchStatus)) {
      const why = watchStatus && watchStatus.error ? `: ${watchStatus.error}` : "";
      messages.push(`Not watching the line${why}. Showing the captured price.`);
    }
    view.warning.hidden = messages.length === 0;
    view.warning.textContent = messages.join(" ");
    view.warning.classList.toggle("bad", bad);
    renderRowTrace(ticket);
  }

  // Why the feed gave no line, for the no-edge view: nothing for the game,
  // no rung at this number, or two markets at it (a team total beside the total).
  function describeFeedMiss(ticket, points) {
    const atNumber = feedLinesFor(ticket, points).length;
    if (atNumber > 1) return `${atNumber} markets at this number for this game, side and book; refusing to guess which is the ticket's`;
    const held = feedPointsHeld(ticket);
    const where = points == null ? "no line" : `no line at ${fmtPoints(points)}`;
    if (held.length === 0) return `${where} for this game, side and book; it holds nothing for them`;
    return `${where} for this game, side and book; it holds ${held.map(fmtPoints).join(", ")}`;
  }

  // A captured line with no edge anywhere — not the cell, not its row, not the
  // Edges feed at that price — is unpriced: the same view a failed read uses,
  // with the cell's fields from page.js so a moved field can be spotted.
  function renderNoEdge(ticket, line) {
    const copy = ERROR_COPY.no_fair;
    view.errorTitle.textContent = copy.title;
    view.errorHint.textContent = copy.hint;
    const feedNote = line.feedLine
      ? ` Edges feed: ${fmtAmerican(line.feedLine.price)}${line.feedLine.ge == null ? ", no edge" : ` with ${fmtPct(line.feedLine.ge)}`}.`
      : scannerState ? ` Edges feed: ${describeFeedMiss(ticket, line.points)}.` : " Edges feed: not loaded yet.";
    view.errorDetail.textContent = `${ticket.sideLabel} ${fmtAmerican(ticket.price)} @ ${ticket.book.name}: ${ticket.noEdgeDetail || "no edge on the cell"} (${new Date(ticket.capturedAt).toLocaleTimeString()}).${feedNote}`;
    show("error");
  }

  // Which grid rows the capture chose between and which one the watcher is
  // reading. Shown only when the panel has a concrete reason to doubt it is
  // following the clicked rung (live 2026-09-12: Under 19.5 captured, Under
  // 2.5 watched), never on an ordinary line move:
  //   - the watcher is reading a row of a different SHAPE than capture picked
  //     (top-level vs an Alts child) — the shape of that bug, and a comparison
  //     rather than an absolute, so a grid where every row is a child is quiet;
  //   - an alt ticket's watched points changed, which cannot happen while the
  //     watcher is re-finding the rung by its number;
  //   - capture found two rows TIED at the best rank, so grid order decided.
  function rowTraceReason(ticket) {
    const current = currentOf(ticket);
    if (current && current.rowTop != null && ticket.watch && ticket.watch.rowTop != null
      && current.rowTop !== ticket.watch.rowTop) {
      return "the watcher is reading a different row than the capture";
    }
    if (ticket.isAlt && current && !current.offBoard && current.points !== ticket.points) {
      return "the watched rung is not the captured number";
    }
    if (ticket.rowResolution && ticket.rowResolution.ambiguous) {
      return "two grid rows tied for this market";
    }
    return null;
  }

  function renderRowTrace(ticket) {
    const resolution = ticket.rowResolution;
    const reason = resolution ? rowTraceReason(ticket) : null;
    view.rowTrace.hidden = !reason;
    if (!reason) return;
    const watched = currentOf(ticket);
    const watching = watched && watched.row ? ` Watching ${watched.row}.` : "";
    view.rowTraceReason.textContent = `Row trace: ${reason}`;
    view.rowTraceDetail.textContent = `Script ${resolution.build}: ${resolution.trace}.${watching}`;
  }

  function renderTicket() {
    const { ticket, settings, watchStatus } = state;
    const { line, result, reason } = computeStake(ticket, settings);
    if (!line.moved && line.edgePct == null) {
      renderNoEdge(ticket, line);
      return;
    }
    renderWarning(ticket, watchStatus, line);

    view.sideLabel.textContent = ticket.sideLabel;
    view.betLine.textContent = `${describeSide(ticket)}${periodSuffix(ticket)}${ticket.isAlt ? " \u00b7 alt line" : ""}`;
    view.eventLine.textContent = describeMatchup(ticket);
    view.startLine.textContent = fmtStart(ticket.eventStart);
    const betFlag = ticketBetFlag(ticket, line);
    renderBetBanner(betFlag);

    view.book.textContent = ticket.book.name;
    view.price.textContent = fmtPriceBoth(asBookLine(line.price, line.sourceFormat, line.sourcePrice));
    // The fair is Unabated's own American number; there is no more exact source for it.
    view.fair.textContent = line.fair == null ? "unknown" : fmtPriceBoth(asBookLine(line.fair, 1, null));
    view.edge.textContent = line.edgePct == null ? "—" : fmtPct(line.edgePct / 100);

    view.stake.classList.remove("no-edge");
    view.payoutRow.hidden = true;
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
    // Sets view.stake to the number to act on when held bets changed it.
    const advice = betFlag.advice;
    renderStakeExposure(advice);

    // Payout = stake x decimal odds at the book's American price; "to win" is
    // the profit on top of it. Both describe the number shown above them, so a
    // top-up prices the top-up and an at-size line shows no payout at all.
    let payoutText = "";
    const acted = result ? betsView.suggestedBetAmount(advice) : null;
    if (acted != null && acted > 0) {
      const payout = acted * kelly.americanToDecimal(line.price);
      view.profit.textContent = fmtDollars(payout - acted);
      view.payout.textContent = fmtDollars(payout);
      view.payoutRow.hidden = false;
      payoutText = ` | to win $${(payout - acted).toFixed(2)} | payout $${payout.toFixed(2)}`;
    }
    const contractsText = renderContracts(acted, line);

    const stakeText = acted != null ? acted.toFixed(2) : "n/a";
    lastCopyText = `${ticket.sideLabel}${periodSuffix(ticket)} ${fmtPriceBoth(asBookLine(line.price, line.sourceFormat, line.sourcePrice))} @ ${ticket.book.name} | fair ${line.fair == null ? "?" : fmtPriceBoth(asBookLine(line.fair, 1, null))} | edge ${line.edgePct == null ? "?" : fmtPct(line.edgePct / 100)} | stake $${stakeText}${copyExposureText(advice)}${contractsText}${payoutText} | ${describeMatchup(ticket)}`;
    view.copyStatus.textContent = "";
    show("ticket");
  }

  // Under the dollar figure, the limit order it means on an exchange: "466
  // contracts @ 53¢ · $246.98", with the count floored so the cost never
  // passes the stake (kelly.contractOrder). Only a line priced in contracts
  // (Kalshi, Novig) gets the row; a sportsbook line keeps just the dollars.
  // `acted` is the number to act on, so a top-up shows the top-up's contracts.
  // Returns what Copy appends, "" when there is no row.
  function renderContracts(acted, line) {
    view.contracts.hidden = true;
    view.contracts.classList.remove("under");
    view.contracts.replaceChildren();
    if (acted == null || acted <= 0) return "";
    const order = kelly.contractOrder({ stake: acted, bookPrice: line.price, sourceFormat: line.sourceFormat, sourcePrice: line.sourcePrice });
    if (!order) return "";
    view.contracts.hidden = false;
    if (order.contracts === 0) {
      view.contracts.classList.add("under");
      view.contracts.textContent = `under 1 contract @ ${order.priceCents}\u00a2`;
      return ` | under 1 contract @ ${order.priceCents}\u00a2`;
    }
    const count = document.createElement("span");
    count.textContent = `${order.contracts.toLocaleString("en-US")} contract${order.contracts === 1 ? "" : "s"} @ ${order.priceCents}\u00a2`;
    const cost = document.createElement("span");
    cost.className = "cost";
    cost.textContent = ` \u00b7 ${fmtDollars(order.costDollars)}`;
    view.contracts.append(count, cost);
    return ` | ${order.contracts} contract${order.contracts === 1 ? "" : "s"} @ ${order.priceCents}\u00a2`;
  }

  // The ticket's open bets and what they do to its stake: the row-style flag
  // {tier, matches, advice}. `line` is the priced line (the current number
  // when the line moved), so the bets are sized against what is on offer now.
  function ticketBetFlag(ticket, line) {
    // The ticket's event's venue id map lets an id-joined bet whose team
    // names do not resolve be placed on a side (bets.sideIndexByVenueId).
    const event = scannerState && ticket.eventId != null ? scannerState.events[ticket.eventId] : null;
    const matchLine = { ...betsView.ticketAsLine(ticket), points: line.points ?? ticket.points ?? null, venueIds: event ? event.venueIds ?? null : null };
    const { matches } = betsLib.matchBets(matchLine, state.betRecords, { lines: boardLines() });
    const advice = betsView.stakeAdvice({
      line: matchLine, price: line.price, edgePct: line.edgePct,
      bankroll: state.settings.bankroll, multiplier: state.settings.multiplier,
      matches, ladderOf: ladderReader(ticket.eventId),
    });
    return { tier: matches.length ? matches[0].tier : null, matches: advice.matches, advice };
  }

  // Every open bet on this line's game, bets in the math first: this line,
  // same side, the other side (red); then the ones that are not sized, grey,
  // with why. Nothing when none.
  function renderBetBanner(flag) {
    const { shown, more } = betsView.bannerLines(betsView.relatedLines(flag));
    const items = shown.map((related) => {
      const div = document.createElement("div");
      const against = related.tier === "opposite" || related.tier === "related_opposite";
      div.className = `bet-match tier-${related.tier}${!related.inMath ? " not-sized" : against ? " bad" : ""}`;
      const kind = document.createElement("span");
      kind.className = "k";
      kind.textContent = related.tag;
      const text = document.createElement("span");
      text.textContent = related.text;
      div.append(kind, text);
      return div;
    });
    if (more > 0) {
      const div = document.createElement("div");
      div.className = "bet-match more";
      div.textContent = `+${more} more on this game (Bets tab)`;
      items.push(div);
    }
    view.betsBanner.replaceChildren(...items);
    view.betsBanner.hidden = items.length === 0;
    view.betsBannerHead.hidden = items.length === 0;
  }

  // What the Copy button adds after "stake $X": the verb and the standalone
  // size, so the clipboard says the number was sized against held bets.
  function copyExposureText(advice) {
    const line = betsView.stakeAdviceLine(advice);
    return line ? ` (${line})` : "";
  }

  // Under the stake, the same three pieces as the Edges rail: the position in
  // dollars, the number to act on (the big figure, its label says "Bet" or
  // "Add to your position"), and what the stake would be alone.
  function renderStakeExposure(advice) {
    const words = betsView.stakeAdviceWords(advice);
    const position = [
      advice.held > 0 ? `held ${betsLib.formatStake(advice.held)}` : null,
      advice.against > 0 ? `against ${betsLib.formatStake(advice.against)}` : null,
      words ? words.alone : null,
    ].filter(Boolean);
    view.stakeExposure.classList.toggle("against", advice.against > 0 && advice.held === 0);
    view.stakeExposure.hidden = position.length === 0;
    view.stakeExposure.textContent = position.join(" \u00b7 ");
    view.stakeLabel.textContent = !words ? "Bet"
      : advice.bet === 0 ? "Already at full size"
        : words.verb === "add" ? "Add to your position" : "Bet";
    // The stake shown is the number to act on, not the standalone size.
    if (words) view.stake.textContent = fmtDollars(advice.bet);
  }

  // Two kinds of capture error need opposite advice: no_fair is Unabated
  // having no number for that line (normal); read_failed means the page changed.
  const ERROR_COPY = {
    no_fair: {
      title: "No Unabated fair for this line",
      hint: "Neither the clicked cell, its grid row, nor the Edges feed (at this price) carries an edge for this line, so there is nothing to size against. This is normal for lopsided moneylines and exchange-only lines. If the Edges tab lists this line at this price, the detail above says which fields the cell carried — send it along.",
    },
    read_failed: {
      title: "Could not read this cell",
      hint: "Click the price again. If it keeps failing, Unabated's page changed; see README troubleshooting.",
    },
  };

  // One loud state above the tabs when page.js's load check found Unabated's
  // price cells no longer carrying the props a ticket is read from. Gated on
  // the capture script being alive, so a verdict from a tab that has since
  // been closed does not outlive it — and on a board with rows, since the
  // check publishes "checking" and stops when there are none.
  function renderShapeBanner() {
    const check = state.pageCheck;
    const changed = Boolean(check && check.status === "changed" && pageScriptAlive());
    view.shapeBanner.hidden = !changed;
    if (!changed) return;
    const detail = check.detail ? ` ${check.detail}.` : "";
    view.shapeBanner.textContent = `Unabated changed: ${check.message}.${detail} Clicking a price will not produce a ticket until page.js is updated for the new bundle.`;
  }

  function render() {
    renderShapeBanner();
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

  // How long ago the book last changed this line; the reader's stale-line tell.
  function fmtLineAge(modifiedMs) {
    if (modifiedMs == null) return "line age unknown";
    const ms = Date.now() - modifiedMs;
    if (ms < 60 * 1000) return "line just changed";
    if (ms < 60 * 60 * 1000) return `line ${Math.round(ms / 60000)}m old`;
    if (ms < 48 * 60 * 60 * 1000) return `line ${Math.round(ms / 3600000)}h old`;
    return `line ${Math.round(ms / 86400000)}d old`;
  }

  function fmtLiquidity(value) {
    return value == null ? "" : `liq ${value.toLocaleString("en-US", { style: "currency", currency: "USD", maximumFractionDigits: 0 })}`;
  }

  function pageScriptAlive() {
    const ready = state.pageReady;
    return Boolean(ready && Date.now() - ready.at < PAGE_READY_STALE_MS);
  }

  function liveBooks() {
    if (!scannerState) return [];
    return Object.values(scannerState.books).filter((book) => book.isLive && book.id !== feed.UNABATED_LINE_BOOK_ID)
      .sort((a, b) => a.name.localeCompare(b.name));
  }

  function unabatedSelection() {
    const filter = state.booksFilter;
    const fresh = filter && typeof filter.at === "number" && Date.now() - filter.at < BOOKS_FILTER_STALE_MS;
    return fresh && Array.isArray(filter.bookIds) && filter.bookIds.length ? filter.bookIds : null;
  }

  function defaultBookIds() {
    if (!scannerState) return [];
    const wanted = new Set(DEFAULT_BOOK_NAMES);
    return Object.values(scannerState.books).filter((book) => wanted.has(book.name)).map((book) => book.id);
  }

  // Which books the list is restricted to: the user's own ticks when they
  // have made any, else the default books until "My Unabated selection" is
  // chosen, which follows the selection page.js published, else every live
  // book. Bet types are the panel's own checkboxes.
  function effectiveFilter() {
    const settings = state.edgeSettings;
    const selection = unabatedSelection();
    let mode;
    let bookIds;
    if (Array.isArray(settings.bookIds)) {
      mode = "custom";
      bookIds = new Set(settings.bookIds);
    } else if (settings.bookIds === undefined) {
      mode = "default";
      bookIds = new Set(defaultBookIds());
    } else if (selection) {
      mode = "unabated";
      bookIds = new Set(selection);
    } else {
      mode = "all";
      bookIds = null;
    }
    return { mode, bookIds, betTypeIds: new Set(settings.betTypes), filter: state.booksFilter };
  }

  function describeFilter(effective) {
    const filter = effective.filter;
    const live = liveBooks().length;
    const parts = [];
    if (effective.mode === "custom") {
      parts.push(`books: your ${effective.bookIds.size} ticks below (of ${live} live)`);
    } else if (effective.mode === "default") {
      parts.push(`books: the ${effective.bookIds.size} default books (of ${live} live; tick below to change)`);
    } else if (effective.mode === "unabated") {
      parts.push(`books: your Unabated selection, ${effective.bookIds.size} books (read ${fmtAge(Date.now() - filter.at)})`);
    } else if (!pageScriptAlive()) {
      parts.push(`books: all ${live} live (no Unabated odds tab is running the capture script; open or reload one to default to your selection, or tick books below)`);
    } else if (filter && filter.lastError) {
      parts.push(`books: all ${live} live (Unabated selection unreadable ${fmtAge(Date.now() - (filter.lastErrorAt || 0))}: ${filter.lastError})`);
    } else {
      parts.push(`books: all ${live} live (waiting for the Unabated tab's first read)`);
    }
    parts.push(`bets: ${Array.from(effective.betTypeIds).map((id) => feed.BET_TYPES[id]).join("/") || "none"}`);
    if (state.edgeSettings.minLiquidity > 0) parts.push(`exchange liq ≥ ${fmtDollars(state.edgeSettings.minLiquidity)}`);
    parts.push(describeAltFilter(state.edgeSettings));
    return parts.join(" · ");
  }

  function describeAltFilter(settings) {
    if (!settings.includeAlts) return "alts: off";
    const gates = [];
    if (settings.altMaxDistance > 0) gates.push(`within ${settings.altMaxDistance} pts of main`);
    return `alts: on${gates.length ? ` (${gates.join(", ")})` : " (no gates)"}`;
  }

  // ---- books checkboxes ----------------------------------------------------

  let booksListSignature = null;

  // Rebuild the checkbox list only when the set of live books changes (a
  // resync), otherwise just sync the ticks, so a click never loses its target.
  function renderBooksList(effective) {
    const books = liveBooks();
    const signature = books.map((book) => book.id).join(",");
    if (signature !== booksListSignature) {
      booksListSignature = signature;
      view.edgesBooks.replaceChildren(...books.map((book) => {
        const label = document.createElement("label");
        const input = document.createElement("input");
        input.type = "checkbox";
        input.dataset.book = String(book.id);
        label.append(input, ` ${book.name}`);
        label.title = book.name;
        return label;
      }));
    }
    for (const input of view.edgesBooks.querySelectorAll("input")) {
      const id = Number(input.dataset.book);
      input.checked = effective.bookIds ? effective.bookIds.has(id) : true;
    }
    const count = effective.bookIds ? effective.bookIds.size : books.length;
    const source = { custom: "your ticks", default: "default books", unabated: "Unabated selection", all: "all live" }[effective.mode];
    view.edgesBooksMode.textContent = `${count} of ${books.length} (${source})`;
    view.booksUnabated.disabled = !unabatedSelection();
  }

  function setBookIds(bookIds) {
    state.edgeSettings = { ...state.edgeSettings, bookIds };
    chrome.storage.local.set({ edges: state.edgeSettings });
    // A newly ticked book brings lines the alert log has never seen: baseline them, don't ping.
    alertsBaselined = false;
    renderEdges();
    processAlerts().catch((error) => console.error("[unabated-ticket] alerts failed", error));
  }

  // The dropdown stays open while you tick; a click anywhere else closes it.
  document.addEventListener("click", (event) => {
    const dropdown = document.getElementById("edges-books-dropdown");
    if (dropdown && dropdown.open && !dropdown.contains(event.target)) dropdown.open = false;
  });

  view.edgesBooks.addEventListener("change", () => {
    const ticked = Array.from(view.edgesBooks.querySelectorAll("input:checked")).map((input) => Number(input.dataset.book));
    setBookIds(ticked);
  });
  view.booksDefault.addEventListener("click", () => setBookIds(undefined));
  view.booksUnabated.addEventListener("click", () => setBookIds(null));
  view.booksAll.addEventListener("click", () => setBookIds(liveBooks().map((book) => book.id)));
  view.booksNone.addEventListener("click", () => setBookIds([]));

  // Quarter-Kelly with nothing held, never more than the line's resting
  // liquidity: it picks a card's best line and sorts "by stake".
  function stakeFor(row) {
    if (row.edgePct == null) return null;
    try {
      const kellyStake = kelly.kellyStakeFromEdge({ bookPrice: row.price, edgePct: row.edgePct, bankroll: state.settings.bankroll, multiplier: state.settings.multiplier }).stake;
      return betsView.capAtLiquidity(kellyStake, row.liquidity).stake;
    } catch (_error) {
      return null;
    }
  }

  // Everything selectEdges needs except the edge threshold (list and alerts differ there).
  function edgeSelectionOptions(effective) {
    const settings = state.edgeSettings;
    return {
      periods: new Set(settings.periods),
      betTypes: effective.betTypeIds,
      bookIds: effective.bookIds,
      now: Date.now(),
      maxLineAgeMs: settings.maxLineAgeHours * 3600 * 1000,
      includeAlts: settings.includeAlts,
      altMaxDistance: settings.altMaxDistance,
      minLiquidity: settings.minLiquidity,
    };
  }

  // (period, axis) -> Unabated's fair ladder for one event, built from the
  // feed's lines on first use and kept until the next scanner update.
  function ladderReader(eventId) {
    return (period, axis) => {
      const periodTypeId = ladderLib.periodTypeIdOf(period);
      if (!scannerState || eventId == null || periodTypeId == null) return null;
      const cacheKey = `${eventId}|${periodTypeId}|${axis}`;
      if (!ladderCache.has(cacheKey)) {
        if (!linesByEventCache) linesByEventCache = ladderLib.groupLinesByEvent(Object.values(scannerState.lines));
        ladderCache.set(cacheKey, ladderLib.buildLadder(linesByEventCache.get(eventId), { periodTypeId, axis }));
      }
      return ladderCache.get(cacheKey);
    };
  }

  // Each row gets `bet` = {tier, matches, advice} from the open bet records:
  // what you hold on that market and the stake sized against it (conditional
  // Kelly, #130). No row is ever hidden for being bet — the edge still being
  // there after you bet it is information, and the stake column carries the top-up.
  function withBetFlags(rows) {
    const flags = betsLib.annotateRows(rows, state.betRecords, { lines: boardLines() });
    return rows.map((row, index) => {
      const flag = flags[index];
      const advice = betsView.stakeAdvice({
        line: row, price: row.price, edgePct: row.edgePct,
        bankroll: state.settings.bankroll, multiplier: state.settings.multiplier,
        matches: flag.matches, ladderOf: ladderReader(row.eventId), liquidity: row.liquidity,
      });
      return { ...row, bet: { tier: flag.tier, matches: advice.matches, advice } };
    });
  }

  // Sort key for "by my exposure": dollars in the math on the market, held or against.
  function exposureDollars(row) {
    return row.bet ? row.bet.advice.held + row.bet.advice.against : 0;
  }

  // Min suggested bet: gates on what the rail says to bet now (the stake
  // sized against what is held), not the standalone size. The list and alerts share it.
  function meetsMinStake(row) {
    const minStake = state.edgeSettings.minStake;
    return minStake === 0 || betsView.suggestedBetAmount(row.bet.advice) >= minStake;
  }

  function currentEdgeRows() {
    if (!scannerState) return [];
    const effective = effectiveFilter();
    const settings = state.edgeSettings;
    const selected = feed.selectEdges(scannerState, { ...edgeSelectionOptions(effective), minEdge: settings.minEdgePct / 100 })
      .map((row) => ({ ...row, stake: stakeFor(row) }));
    const rows = withBetFlags(selected).filter(meetsMinStake);
    if (settings.sortBy === "stake") rows.sort((a, b) => (b.stake ?? -1) - (a.stake ?? -1) || b.edgePct - a.edgePct);
    if (settings.sortBy === "start") rows.sort((a, b) => a.eventStartMs - b.eventStartMs || b.edgePct - a.edgePct);
    if (settings.sortBy === "exposure") rows.sort((a, b) => exposureDollars(b) - exposureDollars(a) || b.edgePct - a.edgePct);
    return rows;
  }

  // Cards: the best line of each (game, market, side) is always the highest
  // stake; the panel's sort orders the cards through that line.
  function groupsOf(rows) {
    const groups = feed.groupEdges(rows, (row) => row.stake);
    const sortBy = state.edgeSettings.sortBy;
    if (sortBy === "edge") groups.sort((a, b) => b.best.edgePct - a.best.edgePct || a.eventStartMs - b.eventStartMs);
    if (sortBy === "start") groups.sort((a, b) => a.eventStartMs - b.eventStartMs || b.best.edgePct - a.best.edgePct);
    if (sortBy === "exposure") groups.sort((a, b) => exposureDollars(b.best) - exposureDollars(a.best) || b.best.edgePct - a.best.edgePct);
    return groups;
  }

  // Cards the user has opened; survives the 5s re-render, not a panel reload.
  const expandedGroups = new Set();

  // "held $300" and/or "against $200", or "game", with every match's label as the tooltip.
  function betBadges(flag) {
    return betsView.badges(flag).map(({ kind, text }) => {
      const badge = document.createElement("span");
      badge.className = `tag ${kind}`;
      badge.textContent = text;
      badge.title = flag.matches.map((match) => match.label).join("\n");
      return badge;
    });
  }

  // Edge magnitude in three steps, so a +6% and a +1.1% never read the same:
  // the row's left stripe and the figure both take their colour from here.
  function edgeTier(edgePct) {
    if (edgePct >= 4) return "hot";
    if (edgePct >= 2) return "warm";
    return "thin";
  }

  // The rail under the edge: the number to act on, with the verb on it, then
  // one small line — that it is all the liquidity there is, and what the
  // stake would be with nothing held. "add $250" is not the same instruction
  // as "bet $250" and must not look like it.
  function fillStakeCell(cell, row) {
    const advice = row.bet ? row.bet.advice : null;
    const words = betsView.stakeAdviceWords(advice);
    cell.classList.toggle("at-size", Boolean(words) && advice.bet === 0);
    if (!words) {
      cell.textContent = row.stake == null ? "—" : `bet ${fmtDollars(row.stake)}`;
      return;
    }
    cell.append(`${words.verb} ${fmtDollars(advice.bet)}`);
    const noteText = [words.cap, words.alone].filter(Boolean).join(" · ");
    if (!noteText) return;
    const note = document.createElement("small");
    note.textContent = noteText;
    cell.append(" ", note);
  }

  // The bets already on this game, as a labelled section of the row rather
  // than a loose line. Every line names the bet and how it relates to this
  // line; past three the rest are a count, as the Ticket banner does it.
  const RELATED_LINES_ON_A_ROW = 3;

  function relatedBlock(flag) {
    const all = betsView.relatedLines(flag);
    if (!all.length) return null;
    const lines = all.slice(0, RELATED_LINES_ON_A_ROW);
    const block = document.createElement("div");
    block.className = `related-block${all.some((line) => line.inMath && (line.tier === "opposite" || line.tier === "related_opposite")) ? " against" : ""}`;
    const head = document.createElement("div");
    head.className = "related-head";
    head.textContent = `Related bets · ${all.length}`;
    block.append(head, ...lines.map((line) => {
      const div = document.createElement("div");
      div.className = `related-line tier-${line.tier}${line.inMath ? "" : " not-sized"}`;
      const tag = document.createElement("span");
      tag.className = "related-tag";
      tag.textContent = line.tag;
      const text = document.createElement("span");
      text.textContent = line.text;
      div.append(tag, text);
      return div;
    }));
    if (all.length > lines.length) {
      const more = document.createElement("div");
      more.className = "related-more";
      more.textContent = `+${all.length - lines.length} more on this game (Bets tab)`;
      block.append(more);
    }
    return block;
  }

  // Time to first pitch is a decision, not a footnote: inside 12 hours it
  // warms, inside 2 it goes red.
  function untilEl(startMs) {
    const span = document.createElement("span");
    const hours = Number.isFinite(startMs) ? (startMs - Date.now()) / 3600000 : null;
    span.className = hours == null ? "" : hours <= 2 ? "soon urgent" : hours <= 12 ? "soon" : "";
    span.textContent = fmtUntil(startMs);
    return span;
  }

  // One book's line, as a full row.
  function renderEdgeRow(row) {
    const li = document.createElement("li");
    const tier = edgeTier(row.edgePct);
    li.className = `edge-row tier-${tier}${row.isBlurred ? " blurred" : ""}${row.key === lastClickedKey ? " last-clicked" : ""}`;
    li.dataset.key = row.key;
    li.append(...rowParts(row, tier));
    return li;
  }

  // The two columns every row and card share: content, then the rail of
  // numbers that line up down the list.
  function rowParts(row, tier) {
    const main = document.createElement("div");

    const side = document.createElement("div");
    side.className = "edge-side";
    side.append(...betBadges(row.bet));
    if (row.isAlt) {
      const badge = document.createElement("span");
      badge.className = "tag";
      badge.textContent = "alt";
      badge.title = `Alternate line; this book's main number is ${fmtPoints(row.mainPoints)}`;
      side.append(badge);
    }
    side.append(row.sideLabel);

    const meta = document.createElement("div");
    meta.className = "edge-meta";
    meta.append(
      `${row.betType}${row.period === "FG" ? "" : ` · ${row.period}`}${row.isAlt ? ` · alt of ${fmtPoints(row.mainPoints)}` : ""}`
        + `${row.rotation != null ? ` · rot ${row.rotation}` : ""} · `,
      `${describeMatchup(row)} · ${fmtStart(row.eventStart)} · `,
      untilEl(row.eventStartMs),
    );

    const book = document.createElement("div");
    book.className = "edge-book";
    const price = document.createElement("span");
    price.className = "price";
    price.textContent = `${row.book.name} ${fmtPriceBoth(asBookLine(row.price, row.sourceFormat, row.sourcePrice))}`;
    const age = document.createElement("span");
    age.className = "age";
    age.textContent = ` · ${[fmtLineAge(row.modifiedMs), fmtLiquidity(row.liquidity)].filter(Boolean).join(" · ")}`;
    book.append(price, age);
    main.append(side, meta, book);

    const rail = document.createElement("div");
    rail.className = "edge-rail";
    const pct = document.createElement("span");
    pct.className = `edge-pct tier-${tier}`;
    pct.textContent = fmtPct(row.edgePct / 100);
    const stake = document.createElement("span");
    stake.className = "edge-stake";
    fillStakeCell(stake, row);
    rail.append(pct, stake);

    return [main, rail, ...[relatedBlock(row.bet)].filter(Boolean)];
  }

  // A line inside a card's expander: price and edge only.
  function renderGroupLine(row, cardSideLabel) {
    const li = document.createElement("li");
    li.className = `group-line${row.isBlurred ? " blurred" : ""}${row.key === lastClickedKey ? " last-clicked" : ""}`;
    li.dataset.key = row.key;

    const main = document.createElement("div");
    const price = document.createElement("span");
    price.className = "gl-price";
    const rung = cardSideLabel && row.sideLabel !== cardSideLabel ? `${row.sideLabel} · ` : "";
    price.textContent = `${rung}${row.book.name} ${fmtPriceBoth(asBookLine(row.price, row.sourceFormat, row.sourcePrice))}`
      + `${row.isAlt ? ` · alt of ${fmtPoints(row.mainPoints)}` : ""}`;
    const age = document.createElement("span");
    age.className = "gl-age";
    age.textContent = [fmtLineAge(row.modifiedMs), fmtLiquidity(row.liquidity)].filter(Boolean).join(" · ");
    main.append(price, age);

    const rail = document.createElement("div");
    rail.className = "gl-rail";
    const edge = document.createElement("span");
    edge.className = "gl-edge";
    edge.textContent = fmtPct(row.edgePct / 100);
    const stake = document.createElement("span");
    stake.className = "gl-stake";
    stake.textContent = row.stake == null ? "—" : fmtDollars(row.stake);
    rail.append(edge, stake);

    li.append(main, rail);
    return li;
  }

  // A card is a row built from the market's best line, with the other books
  // and rungs behind its expander. Same two columns as an ungrouped row, so
  // the edge and the stake stay in one column down the whole list.
  function renderGroupCard(group) {
    const best = group.best;
    const tier = edgeTier(best.edgePct);
    const li = document.createElement("li");
    li.className = `edge-row tier-${tier}${best.isBlurred ? " blurred" : ""}`
      + `${group.rows.some((row) => row.key === lastClickedKey) ? " last-clicked" : ""}`;
    li.dataset.group = group.key;
    li.dataset.key = best.key;
    li.append(...rowParts(best, tier));

    if (group.rows.length > 1) {
      const expanded = expandedGroups.has(group.key);
      const more = document.createElement("button");
      more.type = "button";
      more.className = "group-more";
      more.dataset.group = group.key;
      more.textContent = expanded
        ? "▾ hide the other lines"
        : `▸ ${group.bookCount} book${group.bookCount === 1 ? "" : "s"} · ${group.rows.length} lines (+${group.rows.length - 1})`;
      li.append(more);
      if (expanded) {
        const lines = document.createElement("ol");
        lines.className = "group-lines";
        lines.append(...group.rows.slice(1).map((row) => renderGroupLine(row, best.sideLabel)));
        li.append(lines);
      }
    }
    return li;
  }

  function renderEdgesStatus(rows) {
    const status = scannerStatus;
    if (!status) {
      view.edgesStatus.textContent = "Starting the scanner…";
      return;
    }
    const loadedSports = Array.from(new Set(status.leaguesLoaded.map((id) => (feed.LEAGUES[id] || {}).sport).filter(Boolean)));
    const leagues = status.leaguesLoaded.length > 4
      ? [`${status.leaguesLoaded.length} leagues (${loadedSports.map((sport) => feed.SPORTS[sport]).join(", ")})`]
      : status.leaguesLoaded.map((id) => (feed.LEAGUES[id] || { label: `league ${id}` }).label);
    // "built" is the newest snapshot's Last-Modified: minutes or hours here means a stale edge copy, not a slow poll.
    const stale = status.staleLeagues && status.staleLeagues.length
      ? ` (${status.staleLeagues.length} stale: ${status.staleLeagues.map((id) => (feed.LEAGUES[id] || { label: id }).label).join(", ")})`
      : "";
    const built = status.snapshotBuiltAt ? `snapshot built ${fmtAge(Date.now() - status.snapshotBuiltAt)}${stale}` : "no data yet";
    const updated = status.lastUpdateAt ? `${built} · stream ${fmtAge(Date.now() - status.lastUpdateAt)}` : built;
    const polled = status.lastPollAt ? ` · polled ${fmtAge(Date.now() - status.lastPollAt)}` : "";
    const loading = status.loading ? ` · loading ${status.loading.done}/${status.loading.total} leagues` : "";
    view.edgesStatus.textContent = status.phase === "loading" && !status.leaguesLoaded.length
      ? `Loading snapshots…${loading}`
      : `${leagues.join(" · ") || "no leagues"} · ${status.lineCount.toLocaleString()} lines${status.altLineCount ? ` (+${status.altLineCount.toLocaleString()} alts)` : ""} · ${updated}${polled}${loading}`;
    view.edgesError.hidden = !status.error;
    view.edgesError.textContent = status.error || "";
    view.edgesCount.hidden = rows.length === 0;
    view.edgesCount.textContent = String(rows.length);
  }

  function bookNameOf(id) {
    const book = scannerState && scannerState.books[id];
    return book ? `${book.name} (${id})` : `book ${id}`;
  }

  // Books in the filter by name, then what page.js read them from.
  function renderFilterDebug() {
    if (view.edgesFilterDebug.hidden) return;
    const filter = state.booksFilter;
    if (!filter) {
      view.edgesFilterDebug.textContent = "no booksFilter published yet";
      return;
    }
    const names = Array.isArray(filter.bookIds) ? filter.bookIds.map(bookNameOf).join(", ") : "none";
    view.edgesFilterDebug.textContent = `Unabated selection as read from the tab: ${names}\n\n` +
      JSON.stringify({ lastError: filter.lastError, debug: filter.debug }, null, 1);
  }

  view.edgesFilter.addEventListener("click", () => {
    view.edgesFilterDebug.hidden = !view.edgesFilterDebug.hidden;
    renderFilterDebug();
  });

  // Books, markets, minimum edge — the filter in words, for the header chip
  // that stands in for the whole control block.
  function summariseFilter(effective) {
    const sports = Array.from(new Set(state.edgeSettings.leagues.map((id) => (feed.LEAGUES[id] || {}).sport).filter(Boolean)))
      .map((sport) => feed.SPORTS[sport]);
    const books = effective.bookIds ? `${effective.bookIds.size} books` : `all ${liveBooks().length} books`;
    return [
      sports.length > 3 ? `${sports.length} sports` : sports.join("/") || "no sport",
      state.edgeSettings.periods.map((id) => feed.PERIODS[id] || `pt${id}`).join("/"),
      Array.from(effective.betTypeIds).map((id) => feed.BET_TYPES[id]).join("/") || "no market",
      books,
      `≥${state.edgeSettings.minEdgePct}%`,
      state.edgeSettings.minStake > 0 ? `bet ≥${fmtDollars(state.edgeSettings.minStake)}` : null,
      state.edgeSettings.includeAlts ? "+alts" : null,
    ].filter(Boolean).join(" · ");
  }

  function renderEdges() {
    const rows = currentEdgeRows();
    const grouped = state.edgeSettings.groupByMarket;
    const items = grouped ? groupsOf(rows) : rows;
    renderEdgesStatus(items);
    const effective = effectiveFilter();
    view.edgesFilter.textContent = describeFilter(effective);
    view.filtersSummary.textContent = summariseFilter(effective);
    renderBooksList(effective);
    renderFilterDebug();
    // The scanner rebuilds this list every few seconds; without this the reader
    // is thrown back to the top of it mid-scroll.
    const scrollTop = view.tabEdges.scrollTop;
    view.edgesList.replaceChildren(...items.slice(0, MAX_EDGE_ROWS).map((item) => (grouped ? renderGroupCard(item) : renderEdgeRow(item))));
    view.tabEdges.scrollTop = scrollTop;
    const status = scannerStatus;
    const unit = grouped ? "cards" : "lines";
    if (rows.length === 0) {
      view.edgesEmpty.hidden = false;
      view.edgesEmpty.textContent = !status || (status.phase !== "live" && !status.leaguesLoaded.length)
        ? (status && status.phase === "error" ? "Nothing to list: the feed is unavailable (see above)." : "Waiting for the first snapshot…")
        : `No line at or above ${state.edgeSettings.minEdgePct}% edge${state.edgeSettings.minStake > 0 ? ` with a suggested bet of ${fmtDollars(state.edgeSettings.minStake)} or more` : ""} right now.`;
    } else {
      view.edgesEmpty.hidden = items.length > MAX_EDGE_ROWS ? false : true;
      view.edgesEmpty.textContent = items.length > MAX_EDGE_ROWS ? `Showing the top ${MAX_EDGE_ROWS} of ${items.length} ${unit}; raise the minimum edge to see fewer.` : "";
    }
  }


  // ---- row click -> locate on the Unabated tab -----------------------------

  function locateRequestOf(row) {
    return {
      key: row.key, league: row.league, leagueLabel: row.leagueLabel, eventId: row.eventId,
      betTypeId: row.betTypeId, periodTypeId: row.periodTypeId, sideKey: row.sideKey, sideIndex: row.sideIndex,
      bookId: row.book.id, bookName: row.book.name, marketId: row.marketId, points: row.points, price: row.price,
      isAlt: row.isAlt === true, mainPoints: row.mainPoints,
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
    const more = event.target.closest("button.group-more");
    if (more) {
      if (expandedGroups.has(more.dataset.group)) expandedGroups.delete(more.dataset.group);
      else expandedGroups.add(more.dataset.group);
      renderEdges();
      return;
    }
    const li = event.target.closest("li.edge-row, li.group-line");
    if (!li) return;
    const row = currentEdgeRows().find((r) => r.key === li.dataset.key);
    if (!row) return;
    lastClickedKey = row.key;
    for (const marked of view.edgesList.querySelectorAll(".last-clicked")) marked.classList.remove("last-clicked");
    const card = li.closest("li.edge-row");
    if (card) card.classList.add("last-clicked");
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


  // ---- alerts --------------------------------------------------------------

  // Same line = same market, book, side and points; a price change on it is
  // an update to the same alert key and only notifies again if it improved.
  function alertKeyOf(row) {
    return `${row.marketId}:${row.book.id}:${row.sideKey}:${row.points}`;
  }

  // What one notification is about: a line (flat list) or a card's best line
  // (grouped), with the rule for notifying the same key again.
  function alertItems() {
    const rows = alertRows();
    if (!state.edgeSettings.groupByMarket) {
      return rows.map((row) => ({
        key: alertKeyOf(row), row, summary: null,
        improvedOn: (previous) => priceImproved(row.price, previous.price),
      }));
    }
    return groupsOf(rows).map((group) => ({
      key: `group:${group.key}`, row: group.best,
      summary: `${group.bookCount} book${group.bookCount === 1 ? "" : "s"} \u00b7 ${group.rows.length} line${group.rows.length === 1 ? "" : "s"}`,
      // A card pings again only when its best line got better by the card's
      // own ranking, the stake: the best rung pulled and a +944 longshot
      // taking over is a worse card, not news, whatever its edge %.
      improvedOn: (previous) => typeof previous.stake === "number" && typeof group.best.stake === "number" && group.best.stake > previous.stake,
    }));
  }

  function priceImproved(newPrice, oldPrice) {
    try {
      return kelly.americanToDecimal(newPrice) > kelly.americanToDecimal(oldPrice);
    } catch (_error) {
      return false;
    }
  }

  let iconDataUrl = null;
  // chrome.notifications needs an iconUrl; drawn here so the repo carries no binary.
  function notificationIcon() {
    if (iconDataUrl) return iconDataUrl;
    const canvas = document.createElement("canvas");
    canvas.width = 64;
    canvas.height = 64;
    const ctx = canvas.getContext("2d");
    ctx.fillStyle = "#f59e0b";
    ctx.fillRect(0, 0, 64, 64);
    ctx.fillStyle = "#111827";
    ctx.font = "bold 40px system-ui, sans-serif";
    ctx.textAlign = "center";
    ctx.textBaseline = "middle";
    ctx.fillText("U", 32, 34);
    iconDataUrl = canvas.toDataURL("image/png");
    return iconDataUrl;
  }

  async function rememberAlertTarget(notificationId, row) {
    const stored = await chrome.storage.local.get("alertTargets");
    const targets = stored.alertTargets && typeof stored.alertTargets === "object" ? stored.alertTargets : {};
    targets[notificationId] = locateRequestOf(row);
    const ids = Object.keys(targets);
    for (const id of ids.slice(0, Math.max(0, ids.length - ALERT_TARGETS_KEPT))) delete targets[id];
    await chrome.storage.local.set({ alertTargets: targets });
  }

  async function notifyEdge(row, summary) {
    const notificationId = `edge:${row.key}:${Date.now()}`;
    await rememberAlertTarget(notificationId, row);
    const stake = row.stake ?? stakeFor(row);
    const message = [
      `${fmtPct(row.edgePct / 100)} edge`,
      stake == null ? null : `stake ${fmtDollars(stake)}`,
      summary,
      describeMatchup(row),
      fmtUntil(row.eventStartMs),
    ].filter(Boolean).join(" · ");
    await new Promise((resolve) => {
      chrome.notifications.create(notificationId, {
        type: "basic",
        iconUrl: notificationIcon(),
        title: `${row.sideLabel} ${fmtAmerican(row.price)} @ ${row.book.name}${row.isAlt ? ` (alt of ${fmtPoints(row.mainPoints)})` : ""}`,
        message,
        priority: 1,
      }, () => {
        if (chrome.runtime.lastError) console.error("[unabated-ticket] notification failed:", chrome.runtime.lastError.message);
        resolve();
      });
    });
  }

  function pruneAlertLog(now) {
    for (const [key, entry] of Object.entries(alertLog)) {
      if (!entry || typeof entry.at !== "number" || now - entry.at > ALERT_LOG_TTL_MS) delete alertLog[key];
    }
  }

  function alertRows() {
    const selected = feed.selectEdges(scannerState, { ...edgeSelectionOptions(effectiveFilter()), minEdge: state.alertSettings.minEdgePct / 100 })
      .map((row) => ({ ...row, stake: stakeFor(row) }));
    // A line whose stake against what you hold is $0 has nothing to act on, and
    // one below the Min suggested bet is hidden from the list, so neither alerts.
    return withBetFlags(selected).filter((row) => betsView.suggestedBetAmount(row.bet.advice) > 0 && meetsMinStake(row));
  }

  // Runs after every scanner update. Baseline first, then one notification
  // per line that first crosses the alert threshold (or improves its price),
  // at most one per event per ALERT_EVENT_COOLDOWN_MS.
  async function processAlerts() {
    if (!scannerState || !scannerStatus || scannerStatus.phase !== "live") return;
    if (!state.alertSettings.enabled) {
      alertsBaselined = false;
      return;
    }
    if (alertsBusy) return;
    alertsBusy = true;
    try {
      await processAlertsOnce();
    } finally {
      alertsBusy = false;
    }
  }

  async function processAlertsOnce() {
    const now = Date.now();
    const items = alertItems();
    pruneAlertLog(now);
    if (!alertsBaselined) {
      for (const item of items) alertLog[item.key] = { price: item.row.price, stake: item.row.stake, at: now, baseline: true };
      alertsBaselined = true;
      await chrome.storage.local.set({ alertLog });
      return;
    }
    let fired = 0;
    for (const item of items) {
      const previous = alertLog[item.key];
      if (previous && !item.improvedOn(previous)) continue;
      const lastForEvent = eventAlertAt[item.row.eventId] || 0;
      if (now - lastForEvent < ALERT_EVENT_COOLDOWN_MS) continue;
      await notifyEdge(item.row, item.summary);
      alertLog[item.key] = { price: item.row.price, stake: item.row.stake, at: now };
      eventAlertAt[item.row.eventId] = now;
      fired += 1;
    }
    if (fired) await chrome.storage.local.set({ alertLog });
  }

  function readAlertSettingInputs() {
    const minEdgePct = Number(view.alertsMin.value);
    if (!Number.isFinite(minEdgePct) || minEdgePct < 0) return { error: "Alert edge must be zero or more." };
    return { settings: { enabled: view.alertsEnabled.checked, minEdgePct } };
  }

  function fillAlertSettingInputs() {
    view.alertsEnabled.checked = state.alertSettings.enabled;
    view.alertsMin.value = state.alertSettings.minEdgePct;
  }

  function onAlertSettingsInput() {
    const parsed = readAlertSettingInputs();
    view.edgesSettingsError.textContent = parsed.error || "";
    if (parsed.error) return;
    const thresholdChanged = parsed.settings.minEdgePct !== state.alertSettings.minEdgePct;
    state.alertSettings = parsed.settings;
    chrome.storage.local.set({ alerts: parsed.settings });
    // A new threshold re-baselines so lowering it does not fire for everything already listed.
    if (thresholdChanged) alertsBaselined = false;
    processAlerts().catch((error) => console.error("[unabated-ticket] alerts failed", error));
  }

  function sanitizeAlertSettings(stored) {
    const base = { ...DEFAULT_ALERT_SETTINGS };
    if (!stored || typeof stored !== "object") return base;
    if (typeof stored.enabled === "boolean") base.enabled = stored.enabled;
    if (typeof stored.minEdgePct === "number" && stored.minEdgePct >= 0) base.minEdgePct = stored.minEdgePct;
    return base;
  }

  // ---- tabs ----------------------------------------------------------------

  const TABS = ["ticket", "edges", "bets"];
  const paneOf = { ticket: view.tabTicket, edges: view.tabEdges, bets: view.tabBets };
  // Each pane scrolls on its own, but a hidden element is not guaranteed to
  // keep its scrollTop, so the position is remembered explicitly. Without this
  // the Edges list went back to the top every time a capture brought the
  // Ticket tab forward, which is the whole reason for the split panes.
  const paneScroll = { ticket: 0, edges: 0, bets: 0 };

  function showTab(name) {
    const next = TABS.includes(name) ? name : "ticket";
    const previous = state.activeTab;
    if (paneOf[previous] && !paneOf[previous].hidden) paneScroll[previous] = paneOf[previous].scrollTop;
    state.activeTab = next;
    for (const tab of TABS) paneOf[tab].hidden = tab !== next;
    for (const button of view.tabs.querySelectorAll("button[data-tab]")) {
      button.classList.toggle("active", button.dataset.tab === next);
    }
    // The filter toolbar belongs to the Edges tab; it is chrome, not content.
    view.edgesToolbar.hidden = next !== "edges";
    if (next !== "edges") setFiltersOpen(false);
    view.backToEdges.hidden = next !== "ticket" || !lastClickedKey;
    paneOf[next].scrollTop = paneScroll[next];
  }

  view.settingsToggle.addEventListener("click", () => {
    showTab("ticket");
    chrome.storage.local.set({ activeTab: state.activeTab });
    view.settings.scrollIntoView({ block: "end", behavior: "smooth" });
    view.settings.classList.remove("flash");
    // Restart the animation: the class has to leave the element and come back.
    void view.settings.offsetWidth;
    view.settings.classList.add("flash");
  });

  // The drawer, its chip's caret and its aria state are one thing; leaving the
  // tab closed the drawer and left the chip claiming it was open.
  function setFiltersOpen(open) {
    view.edgesControls.hidden = !open;
    view.filtersToggle.setAttribute("aria-expanded", String(open));
    view.filtersToggle.querySelector(".caret").textContent = open ? "\u25B4" : "\u25BE";
  }

  view.filtersToggle.addEventListener("click", () => setFiltersOpen(view.edgesControls.hidden));

  view.backToEdges.addEventListener("click", () => {
    showTab("edges");
    chrome.storage.local.set({ activeTab: state.activeTab });
    renderEdges();
  });

  view.tabs.addEventListener("click", (event) => {
    const button = event.target.closest("button[data-tab]");
    if (!button) return;
    showTab(button.dataset.tab);
    chrome.storage.local.set({ activeTab: state.activeTab });
    if (state.activeTab === "edges") renderEdges();
    if (state.activeTab === "bets") renderBets();
  });

  // ---- edge settings -------------------------------------------------------

  function readEdgeSettingInputs() {
    const leagues = Array.from(view.edgesSports.querySelectorAll("input:checked")).flatMap((input) => feed.leagueIdsOfSport(input.dataset.sport));
    const periods = Array.from(view.edgesPeriods.querySelectorAll("input:checked")).map((input) => Number(input.dataset.period));
    const betTypes = Array.from(view.edgesBetTypes.querySelectorAll("input:checked")).map((input) => Number(input.dataset.bettype));
    const minEdgePct = Number(view.edgesMin.value);
    const minStake = Number(view.edgesMinStake.value);
    const maxLineAgeHours = Number(view.edgesMaxAge.value);
    const altMaxDistance = Number(view.edgesAltDistance.value);
    const minLiquidity = Number(view.edgesMinLiquidity.value);
    if (!Number.isFinite(minEdgePct) || minEdgePct < 0) return { error: "Minimum edge must be zero or more." };
    if (!Number.isFinite(minStake) || minStake < 0) return { error: "Minimum suggested bet must be zero (off) or more." };
    if (!Number.isFinite(maxLineAgeHours) || maxLineAgeHours <= 0) return { error: "Max line age must be above zero hours." };
    if (!Number.isFinite(altMaxDistance) || altMaxDistance < 0) return { error: "Max points from main must be zero (off) or more." };
    if (!Number.isFinite(minLiquidity) || minLiquidity < 0) return { error: "Min liquidity must be zero (off) or more." };
    if (!periods.length) return { error: "Pick at least one period." };
    if (!betTypes.length) return { error: "Pick at least one bet type." };
    return {
      settings: {
        ...state.edgeSettings, leagues, periods, betTypes, minEdgePct, minStake, maxLineAgeHours, minLiquidity, sortBy: view.edgesSort.value,
        includeAlts: view.edgesIncludeAlts.checked, altMaxDistance,
        groupByMarket: view.edgesGroup.checked,
      },
    };
  }

  function fillEdgeSettingInputs() {
    const settings = state.edgeSettings;
    for (const input of view.edgesSports.querySelectorAll("input")) {
      input.checked = feed.leagueIdsOfSport(input.dataset.sport).some((id) => settings.leagues.includes(id));
    }
    for (const input of view.edgesPeriods.querySelectorAll("input")) input.checked = settings.periods.includes(Number(input.dataset.period));
    for (const input of view.edgesBetTypes.querySelectorAll("input")) input.checked = settings.betTypes.includes(Number(input.dataset.bettype));
    view.edgesMin.value = settings.minEdgePct;
    view.edgesMinStake.value = settings.minStake;
    view.edgesMaxAge.value = settings.maxLineAgeHours;
    view.edgesSort.value = settings.sortBy;
    view.edgesIncludeAlts.checked = settings.includeAlts;
    view.edgesAltDistance.value = settings.altMaxDistance;
    view.edgesMinLiquidity.value = settings.minLiquidity;
    view.edgesGroup.checked = settings.groupByMarket;
  }

  function onEdgeSettingsInput() {
    const parsed = readEdgeSettingInputs();
    view.edgesSettingsError.textContent = parsed.error || "";
    if (parsed.error) return;
    const before = state.edgeSettings;
    const leaguesChanged = parsed.settings.leagues.join(",") !== before.leagues.join(",");
    // Widening periods, bet types, the liquidity floor or the alt gates exposes lines the alert log has never seen.
    const scopeChanged = leaguesChanged
      || parsed.settings.periods.join(",") !== before.periods.join(",")
      || parsed.settings.betTypes.join(",") !== before.betTypes.join(",")
      || parsed.settings.includeAlts !== before.includeAlts
      || parsed.settings.altMaxDistance !== before.altMaxDistance
      || parsed.settings.minLiquidity !== before.minLiquidity
      // Alert keys differ between the flat list and cards.
      || parsed.settings.groupByMarket !== before.groupByMarket;
    state.edgeSettings = parsed.settings;
    chrome.storage.local.set({ edges: parsed.settings });
    if (scopeChanged) alertsBaselined = false;
    if (leaguesChanged) {
      scanner.start(parsed.settings.leagues).catch((error) => console.error("[unabated-ticket] scanner restart failed", error));
    }
    renderEdges();
    if (scopeChanged && !leaguesChanged) processAlerts().catch((error) => console.error("[unabated-ticket] alerts failed", error));
  }

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
    if (["edge", "stake", "start", "exposure"].includes(stored.sortBy)) base.sortBy = stored.sortBy;
    if (typeof stored.includeAlts === "boolean") base.includeAlts = stored.includeAlts;
    if (typeof stored.altMaxDistance === "number" && stored.altMaxDistance >= 0) base.altMaxDistance = stored.altMaxDistance;
    // altMinLiquidity is the pre-rename key, from when the floor gated alts only.
    const storedMinLiquidity = typeof stored.minLiquidity === "number" ? stored.minLiquidity : stored.altMinLiquidity;
    if (typeof storedMinLiquidity === "number" && storedMinLiquidity >= 0) base.minLiquidity = storedMinLiquidity;
    if (typeof stored.groupByMarket === "boolean") base.groupByMarket = stored.groupByMarket;
    return base;
  }

  // ---- bets (#114) ---------------------------------------------------------

  // One describeLine-shaped row per event on the board (main lines only):
  // the matcher's ambiguity check, its venue id join (each row carries its
  // event's venueIds) and the unmatched list only need to know which games
  // exist, not every book's price.
  // Built once per scanner update (73 ms over the 140,755 NFL + CFB lines of
  // 2026-09-15) and shared by the Edges rows, the alert pass, the Ticket
  // banner and the Bets tab. The rows only name games — teams, start, venue
  // ids — so a copy from the last update is current; every tier reads the
  // row or ticket being matched, never these.
  function boardLines() {
    if (!scannerState) return [];
    if (boardLinesCache) return boardLinesCache;
    const seen = new Set();
    const rows = [];
    for (const line of Object.values(scannerState.lines)) {
      if (line.isAlt || seen.has(line.eventId)) continue;
      seen.add(line.eventId);
      rows.push(feed.describeLine(line, scannerState));
    }
    boardLinesCache = rows;
    return rows;
  }

  function betsPayload() {
    return state.betsService ? state.betsService.payload : null;
  }

  // Venues read by a content script rather than the service (#116).
  function pageSources() {
    return state.betsNovig ? { novig: state.betsNovig } : {};
  }

  // A new Novig read from storage: merge its records (complete reads are
  // authoritative for the venue) and refresh every view that shows a flag.
  function applyNovigRead(betsNovig) {
    state.betsNovig = betsNovig && typeof betsNovig === "object" ? betsNovig : null;
    if (state.betsNovig) state.betRecords = betsView.mergePageSource(state.betRecords, "novig", state.betsNovig, Date.now(), state.crosswalk);
  }

  let betsPollBusy = false;
  let betsPollTimer = null;

  // One GET of /bets.json. Success replaces the source status and merges the
  // records (bets.js dedupe on native id, team keys, retention prune); failure
  // keeps everything and records since when the service has been unreachable.
  async function pollBets() {
    if (betsPollBusy || document.hidden) return;
    betsPollBusy = true;
    const now = Date.now();
    const previous = state.betsService || {};
    try {
      const response = await fetch(`${state.betsSettings.serviceUrl}/bets.json`, { cache: "no-store" });
      if (!response.ok) throw new Error(`HTTP ${response.status}`);
      const payload = await response.json();
      if (!payload || !Array.isArray(payload.bets)) throw new Error("bets.json has no bets array");
      // The service's crosswalk is the truth; a payload without one (an older service) keeps the stored rows.
      if (Array.isArray(payload.crosswalk)) state.crosswalk = payload.crosswalk;
      state.betRecords = betsView.mergeServicePayload(state.betRecords, { ...payload, crosswalk: state.crosswalk }, now);
      state.betsService = {
        payload: { generatedAt: payload.generatedAt ?? null, sources: payload.sources && typeof payload.sources === "object" ? payload.sources : {} },
        okAt: now, error: null, errorAt: null, unreachableSince: null,
      };
    } catch (error) {
      if (previous.error !== error.message) console.warn("[unabated-ticket] bets service poll failed:", error.message);
      state.betsService = {
        payload: previous.payload || null, okAt: previous.okAt ?? null,
        error: error.message, errorAt: now, unreachableSince: previous.unreachableSince ?? now,
      };
    } finally {
      betsPollBusy = false;
    }
    await persistBets();
    renderBetsFlags();
    learnCrosswalk().catch((error) => console.error("[unabated-ticket] crosswalk learn failed", error));
  }

  // The records, the service state and the crosswalk, as one stored object.
  function persistBets() {
    return chrome.storage.local.set({ betsService: { ...state.betsService, bets: state.betRecords, crosswalk: state.crosswalk } });
  }

  // Every surface that shows a bet flag, after the records or the crosswalk changed.
  function renderBetsFlags() {
    renderBetsHeader();
    if (!state.error) render();
    renderEdges();
    if (state.activeTab === "bets") renderBets();
  }

  // ---- team crosswalk (#118 step 4) -----------------------------------------
  //
  // The panel is the only side that sees both a bet and the board, so it
  // learns; the service keeps the table. After every poll and every scanner
  // update, the open bets the board joined by venue id teach what their
  // venue calls both teams (bets.learnCrosswalk); new rows go to the service
  // in one POST, its reply is the whole table, and the records are re-keyed
  // through it. A refused bet (its name resolves to another team) is logged
  // once. A failed POST waits a minute before the next try.
  const CROSSWALK_RETRY_MS = 60 * 1000;
  let crosswalkBusy = false;
  let crosswalkRetryAt = 0;
  let crosswalkLastError = null;
  const crosswalkConflictsLogged = new Set();

  async function learnCrosswalk() {
    if (crosswalkBusy || Date.now() < crosswalkRetryAt || !scannerState) return;
    const { learned, conflicts } = betsLib.learnCrosswalk(state.betRecords, boardLines(), state.crosswalk);
    for (const conflict of conflicts) {
      const tag = `${conflict.betId}:${conflict.side}`;
      if (crosswalkConflictsLogged.has(tag)) continue;
      crosswalkConflictsLogged.add(tag);
      console.warn(`[unabated-ticket] crosswalk: not learning ${conflict.venue} ${conflict.league} "${conflict.venueTeamKey}" from ${conflict.betId}: ${conflict.reason}`);
    }
    if (!learned.length) return;
    crosswalkBusy = true;
    try {
      const reply = await postCrosswalk("POST", { rows: learned });
      console.info(`[unabated-ticket] crosswalk: learned ${reply.learned} row(s), ${reply.conflicts.length} refused by the service, ${reply.crosswalk.length} held`);
      await applyCrosswalk(reply.crosswalk);
      crosswalkLastError = null;
    } catch (error) {
      crosswalkRetryAt = Date.now() + CROSSWALK_RETRY_MS;
      if (crosswalkLastError !== error.message) console.warn("[unabated-ticket] crosswalk: service write failed:", error.message);
      crosswalkLastError = error.message;
    } finally {
      crosswalkBusy = false;
    }
  }

  // POST (learn) or DELETE (clear) /crosswalk.json; the reply carries the table.
  async function postCrosswalk(method, body) {
    const response = await fetch(`${state.betsSettings.serviceUrl}/crosswalk.json`, {
      method, cache: "no-store",
      headers: body ? { "Content-Type": "application/json" } : {},
      body: body ? JSON.stringify(body) : undefined,
    });
    if (!response.ok) throw new Error(`HTTP ${response.status}`);
    const reply = await response.json();
    if (!reply || !Array.isArray(reply.crosswalk)) throw new Error("crosswalk.json has no crosswalk array");
    return reply;
  }

  // The served table replaces the held one; keys are rebuilt from scratch so
  // a cleared row takes its key back and a new one applies everywhere.
  async function applyCrosswalk(crosswalk) {
    state.crosswalk = crosswalk;
    state.betRecords = betsLib.rekeyRecords(state.betRecords, state.crosswalk);
    await persistBets();
    renderBetsFlags();
  }

  // Clear is two clicks: the first arms the button ("Clear 68 rows?") for a
  // few seconds, the second deletes. No dialog — a side panel is not a place
  // to rely on window.confirm — and a re-render while armed leaves it armed.
  const CROSSWALK_CLEAR_ARM_MS = 6000;
  let crosswalkClearArmedUntil = 0;
  let crosswalkClearTimer = null;

  function crosswalkClearArmed() {
    return Date.now() < crosswalkClearArmedUntil;
  }

  function renderCrosswalkClearButton(count) {
    const button = view.betsCrosswalkClear;
    button.disabled = count === 0;
    button.classList.toggle("armed", crosswalkClearArmed());
    button.textContent = crosswalkClearArmed() ? `Clear ${count} row${count === 1 ? "" : "s"}?` : "Clear";
  }

  async function clearCrosswalk() {
    const count = state.crosswalk.length;
    if (!count) return;
    if (!crosswalkClearArmed()) {
      crosswalkClearArmedUntil = Date.now() + CROSSWALK_CLEAR_ARM_MS;
      clearTimeout(crosswalkClearTimer);
      crosswalkClearTimer = setTimeout(() => renderCrosswalkClearButton(state.crosswalk.length), CROSSWALK_CLEAR_ARM_MS);
      renderCrosswalkClearButton(count);
      return;
    }
    crosswalkClearArmedUntil = 0;
    view.betsCrosswalkClear.disabled = true;
    try {
      const reply = await postCrosswalk("DELETE", null);
      console.info(`[unabated-ticket] crosswalk: cleared ${reply.cleared} row(s)`);
      await applyCrosswalk([]);
      crosswalkConflictsLogged.clear();
      view.betsSettingsError.textContent = "";
    } catch (error) {
      console.warn("[unabated-ticket] crosswalk: clear failed:", error.message);
      view.betsSettingsError.textContent = `Could not clear the crosswalk: ${error.message}`;
    } finally {
      renderCrosswalkClearButton(state.crosswalk.length);
    }
  }

  function renderCrosswalk() {
    const rows = betsView.crosswalkRows(state.crosswalk);
    view.betsCrosswalkCount.textContent = rows.length ? String(rows.length) : "";
    renderCrosswalkClearButton(rows.length);
    view.betsCrosswalk.replaceChildren(...rows.map((row) => {
      const li = document.createElement("li");
      li.title = row.title;
      const main = document.createElement("div");
      const what = document.createElement("div");
      what.className = "bet-what";
      what.textContent = row.what;
      const meta = document.createElement("div");
      meta.className = "bet-meta";
      meta.textContent = row.meta;
      main.append(what, meta);
      li.append(main);
      return li;
    }));
    view.betsCrosswalkEmpty.hidden = rows.length > 0;
    view.betsCrosswalkEmpty.textContent = "Nothing learned yet. Rows appear when an open bet joins the board by its Kalshi event or Novig outcome id.";
  }

  function startBetsPolling() {
    if (betsPollTimer != null) return;
    betsPollTimer = setInterval(() => pollBets().catch((error) => console.error("[unabated-ticket] bets poll failed", error)), BETS_POLL_MS);
    pollBets().catch((error) => console.error("[unabated-ticket] bets poll failed", error));
  }

  function renderBetsHeader() {
    const now = Date.now();
    const open = state.betRecords.filter((bet) => bet.status === "open").length;
    view.betsCount.hidden = open === 0;
    view.betsCount.textContent = String(open);
    view.betsHeader.textContent = betsView.headerLine(state.betRecords, betsPayload(), now, pageSources());
    view.betsHeader.classList.toggle("bad", betsView.serviceStatus(state.betsService, now).unreachable);
  }

  // One line per venue: a dot for freshness, what it holds, how old the last
  // successful pull is. A venue that failed says why, in place of its count.
  function renderBetsSources(now) {
    const service = betsView.serviceStatus(state.betsService, now);
    view.betsService.hidden = !service.unreachable;
    view.betsService.textContent = service.unreachable ? `${service.text}. Start it with unabated_ticket/bets_service/run.sh; the last records it served are still shown.` : "";
    view.betsSources.replaceChildren(...betsView.sourceRows(betsPayload(), now, pageSources()).map((row) => {
      const div = document.createElement("div");
      div.className = `venue fresh-${row.level}`;
      const bets = row.count == null ? null : `${row.count} bets`;
      const trouble = row.error || row.note;
      const note = trouble ? [bets, trouble].filter(Boolean).join(" · ") : bets || "no bets";
      div.title = row.fetchedAt ? `last pull ${new Date(row.fetchedAt).toLocaleTimeString([], { hour: "numeric", minute: "2-digit" })}` : note;
      const dot = document.createElement("span");
      dot.className = "vdot";
      const name = document.createElement("span");
      const venue = document.createElement("span");
      venue.className = "vname";
      venue.textContent = row.venue;
      const detail = document.createElement("span");
      detail.className = "vnote";
      detail.textContent = note;
      name.append(venue, detail);
      const age = document.createElement("span");
      age.className = "vage";
      age.textContent = row.configured ? row.ageText : "—";
      div.append(dot, name, age);
      return div;
    }));
  }

  // `trailing` is the unmatched reason; its presence is also what colours the
  // bet's left edge, so an unmatched bet is visible in the open list too.
  function betItem(bet, trailing, unmatched) {
    const li = document.createElement("li");
    li.className = unmatched ? "unmatched" : "matched";
    const main = document.createElement("div");
    const what = document.createElement("div");
    what.className = "bet-what";
    what.textContent = `${betsLib.describeBet(bet)}`;
    const meta = document.createElement("div");
    meta.className = "bet-meta";
    const venue = bet.venue ? bet.venue.charAt(0).toUpperCase() + bet.venue.slice(1) : "unknown venue";
    const game = bet.awayTeam && bet.homeTeam ? `${bet.awayTeam} @ ${bet.homeTeam}` : null;
    meta.textContent = [venue, game].filter(Boolean).join(" · ");
    main.append(what, meta);

    const rail = document.createElement("div");
    const stake = document.createElement("span");
    stake.className = "bet-stake";
    stake.textContent = bet.stake == null ? "—" : fmtDollars(bet.stake);
    rail.append(stake);
    if (bet.placedAt) {
      const when = document.createElement("small");
      when.className = "bet-when";
      when.textContent = betsLib.formatPlacedAt(bet.placedAt);
      rail.append(when);
    }
    li.append(main, rail);
    if (trailing) {
      const extra = document.createElement("div");
      extra.className = "bet-reason";
      extra.textContent = trailing;
      li.append(extra);
    }
    return li;
  }

  function renderBets() {
    const now = Date.now();
    const scrollTop = view.tabBets.scrollTop;
    renderBetsSources(now);
    const open = state.betRecords.filter((bet) => bet.status === "open")
      .sort((a, b) => Date.parse(b.placedAt || 0) - Date.parse(a.placedAt || 0));
    const unmatched = betsLib.unmatchedReasons(state.betRecords, boardLines());
    const unmatchedIds = new Set(unmatched.map(({ bet }) => bet.id));

    // What the tab opens with: the money, before the plumbing. A bet whose
    // venue reported no stake is counted separately rather than as zero.
    const priced = open.filter((bet) => typeof bet.stake === "number");
    const atRisk = priced.reduce((total, bet) => total + bet.stake, 0);
    const venues = betsView.sourceRows(betsPayload(), now, pageSources()).filter((row) => row.configured).length;
    view.betsRisk.textContent = fmtDollars(atRisk);
    view.betsRiskCaption.textContent = [
      `at risk · ${open.length} open bet${open.length === 1 ? "" : "s"}`,
      `${venues} venue${venues === 1 ? "" : "s"}`,
      priced.length === open.length ? null : `${open.length - priced.length} with no stake reported`,
    ].filter(Boolean).join(" · ");

    view.betsOpenCount.textContent = open.length ? String(open.length) : "";
    view.betsOpen.replaceChildren(...open.map((bet) => betItem(bet, null, unmatchedIds.has(bet.id))));
    view.betsOpenEmpty.hidden = open.length > 0;
    view.betsOpenEmpty.textContent = state.betsService && state.betsService.okAt != null ? "No open bets." : "No bets loaded yet.";
    view.betsUnmatchedCount.textContent = unmatched.length ? String(unmatched.length) : "";
    view.betsUnmatched.replaceChildren(...unmatched.map(({ bet, reason }) => betItem(bet, reason, true)));
    view.betsUnmatchedEmpty.hidden = unmatched.length > 0;
    view.betsUnmatchedEmpty.textContent = open.length ? "Every open bet matches a game on the board." : "";
    renderCrosswalk();
    view.tabBets.scrollTop = scrollTop;
  }

  function readBetsSettingInputs() {
    const serviceUrl = view.betsUrl.value.trim().replace(/\/+$/, "");
    if (!/^https?:\/\/\S+$/.test(serviceUrl)) return { error: "Service URL must start with http:// or https://." };
    return { settings: { serviceUrl } };
  }

  function fillBetsSettingInputs() {
    view.betsUrl.value = state.betsSettings.serviceUrl;
  }

  function onBetsSettingsInput() {
    const parsed = readBetsSettingInputs();
    view.betsSettingsError.textContent = parsed.error || "";
    if (parsed.error) return;
    const urlChanged = parsed.settings.serviceUrl !== state.betsSettings.serviceUrl;
    state.betsSettings = parsed.settings;
    chrome.storage.local.set({ betsSettings: parsed.settings });
    if (urlChanged) pollBets().catch((error) => console.error("[unabated-ticket] bets poll failed", error));
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
    const relay = await chrome.storage.local.get(["ticket", "error", "watchStatus", "pageReady", "pageCheck", "booksFilter", "edges", "alerts", "alertLog", "activeTab", "locateResult", "betsService", "betsSettings", "betsNovig", "teamsIndex"]);
    // The team index from the last session, so bet records resolve before the first snapshot lands.
    teamsLib.loadIndex(relay.teamsIndex);
    teamsSpellingCount = teamsLib.spellingCount();
    state.ticket = relay.ticket || null;
    state.error = relay.error || null;
    state.watchStatus = relay.watchStatus || null;
    state.pageReady = relay.pageReady || null;
    state.pageCheck = relay.pageCheck || null;
    state.booksFilter = relay.booksFilter || null;
    state.locateResult = relay.locateResult || null;
    state.edgeSettings = sanitizeEdgeSettings(relay.edges);
    state.alertSettings = sanitizeAlertSettings(relay.alerts);
    alertLog = relay.alertLog && typeof relay.alertLog === "object" ? relay.alertLog : {};
    state.betsSettings = betsView.sanitizeBetsSettings(relay.betsSettings);
    const storedBets = relay.betsService && typeof relay.betsService === "object" ? relay.betsService : null;
    if (storedBets) {
      // Team keys resolved and the retention window applied on every load, so
      // a grown teams.js table, a grown crosswalk and a passed month all take effect.
      state.crosswalk = Array.isArray(storedBets.crosswalk) ? storedBets.crosswalk : [];
      state.betRecords = betsLib.pruneForRetention(betsLib.rekeyRecords(Array.isArray(storedBets.bets) ? storedBets.bets : [], state.crosswalk), Date.now());
      state.betsService = {
        payload: storedBets.payload || null, okAt: storedBets.okAt ?? null,
        error: storedBets.error ?? null, errorAt: storedBets.errorAt ?? null, unreachableSince: storedBets.unreachableSince ?? null,
      };
    }
    applyNovigRead(relay.betsNovig);
    fillEdgeSettingInputs();
    fillAlertSettingInputs();
    fillBetsSettingInputs();
    showTab(relay.activeTab);
    renderBetsHeader();
    render();
    renderEdges();
    renderLocate();
    if (state.activeTab === "bets") renderBets();
    startBetsPolling();
    await scanner.start(state.edgeSettings.leagues);
  }

  chrome.storage.onChanged.addListener((changes, area) => {
    if (area !== "local") return;
    if ("ticket" in changes) {
      const previous = state.ticket;
      state.ticket = changes.ticket.newValue || null;
      // The watcher no longer rewrites the ticket, so every write here is a
      // capture; the capturedAt check stays as the cheap guard against one.
      if (state.ticket && (!previous || previous.capturedAt !== state.ticket.capturedAt)) showTab("ticket");
    }
    if ("error" in changes) {
      state.error = changes.error.newValue || null;
      if (state.error) showTab("ticket");
    }
    if ("watchStatus" in changes) state.watchStatus = changes.watchStatus.newValue || null;
    if ("pageReady" in changes) state.pageReady = changes.pageReady.newValue || null;
    if ("pageCheck" in changes) state.pageCheck = changes.pageCheck.newValue || null;
    if ("ticket" in changes || "error" in changes || "watchStatus" in changes || "pageReady" in changes || "pageCheck" in changes) render();
    if ("pageReady" in changes) renderEdges();
    if ("booksFilter" in changes) {
      state.booksFilter = changes.booksFilter.newValue || null;
      renderEdges();
    }
    if ("locateResult" in changes) {
      state.locateResult = changes.locateResult.newValue || null;
      if (state.locateResult && state.locating && state.locateResult.at >= state.locating.at) state.locating = null;
      renderLocate();
    }
    // novig_content.js wrote a read of the Novig Portfolio screen (#116).
    if ("betsNovig" in changes) {
      applyNovigRead(changes.betsNovig.newValue);
      renderBetsHeader();
      if (!state.error) render();
      renderEdges();
      if (state.activeTab === "bets") renderBets();
    }
  });

  view.edgesSports.addEventListener("change", onEdgeSettingsInput);
  view.edgesPeriods.addEventListener("change", onEdgeSettingsInput);
  view.edgesBetTypes.addEventListener("change", onEdgeSettingsInput);
  view.edgesMin.addEventListener("input", onEdgeSettingsInput);
  view.edgesMinStake.addEventListener("input", onEdgeSettingsInput);
  view.edgesMaxAge.addEventListener("input", onEdgeSettingsInput);
  view.edgesSort.addEventListener("change", onEdgeSettingsInput);
  view.edgesIncludeAlts.addEventListener("change", onEdgeSettingsInput);
  view.edgesAltDistance.addEventListener("input", onEdgeSettingsInput);
  view.edgesMinLiquidity.addEventListener("input", onEdgeSettingsInput);
  view.edgesGroup.addEventListener("change", onEdgeSettingsInput);
  view.alertsEnabled.addEventListener("change", onAlertSettingsInput);
  view.alertsMin.addEventListener("input", onAlertSettingsInput);
  view.betsUrl.addEventListener("change", onBetsSettingsInput);
  view.betsCrosswalkClear.addEventListener("click", () => clearCrosswalk().catch((error) => console.error("[unabated-ticket] crosswalk clear failed", error)));

  // Nothing polls while the panel is hidden; back in view, the scanner catches
  // up or resyncs and the bets service is polled at once (its tick skips hidden).
  document.addEventListener("visibilitychange", () => {
    if (document.hidden) {
      scanner.pause();
      return;
    }
    scanner.resume().catch((error) => console.error("[unabated-ticket] scanner resume failed", error));
    pollBets().catch((error) => console.error("[unabated-ticket] bets poll failed", error));
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
    // Unconditionally: the banner is above the tabs and its liveness gate is a
    // clock, so it must keep re-evaluating while a failed capture is on screen
    // — which is exactly the state a new Unabated bundle puts the panel in.
    renderShapeBanner();
    if (!state.error) render();
    renderBetsHeader();
    if (state.activeTab === "edges") {
      renderEdges();
      renderLocate();
    }
    if (state.activeTab === "bets") renderBets();
  }, 5000);

  load().catch((error) => {
    view.errorDetail.textContent = error.message;
    show("error");
  });
})();
