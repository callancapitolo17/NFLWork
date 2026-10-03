// Unabated Ticket — side panel: the Ticket tab (one captured bet), the Edges
// tab (every positive-edge line across the enabled leagues), the Bets tab
// (open bets from the bets service) and the Teasers tab (Buckeye 6-point
// teasers, teaser.js).
//
// Reads: chrome.storage.local {ticket, error, watchStatus, pageCheck, pageReady,
// booksFilter, locateResult, "liveEdges:<league>"} (written by content.js) and {bankroll,
// multiplier, edges, alerts, alertLog, activeTab, betsService, betsSettings,
// teaserRefs, teaserBlocked} (written here; teaserRefs = the Teasers list's
// reference fairs, rewritten on each rebuild; teaserBlocked = the CFB markets
// marked can't tease, each dropped once its game starts).
// Writes: chrome.storage.local settings, {locate} (row click, via locate.js,
// which also focuses the Unabated tab), {alertLog, alertTargets} and Chrome
// notifications for new edges; {betsService} after every bets-service poll
// (records + the service's team crosswalk + the saved fill fairs); POST
// /crosswalk.json to the bets service with the team rows an id join taught,
// DELETE it on Clear (#118 step 4); POST /fill_fairs.json with the fair a
// new bet's line had when it was placed (fillfair.js, the service stores it
// once per bet, insert-only). Re-renders on storage.onChanged.
// Network: the scanner (scanner.js) fetches Unabated's public snapshots — the
// Edges tab's leagues plus NFL and CFB for the Teasers tab — while this page
// is open; it pauses when the panel is hidden and stops when it closes.
// The bets tab (#114) polls the local bets service (betsSettings.serviceUrl,
// default http://127.0.0.1:8094) every 30 s on the same visibility rule —
// never from the service worker. Matching is bets.js; presentation helpers
// are betsview.js; the Edges rows' selection, sizing, sort and words are
// edgerows.js (shared with the server runner). Bet105 (2026-09-29) is the
// one venue this page reads itself: every 5 min while visible it fetches the
// account's open bets from app.bet105.ag on Cal's own login in this Chrome
// (bet105.js) and POSTs them to the service's /bet105.json, which parses and
// stores them like any other.

(function () {
  "use strict";

  const kelly = globalThis.UnabatedKelly;
  const feed = globalThis.UnabatedFeed;
  const betsLib = globalThis.UnabatedBets;
  const betsView = globalThis.UnabatedBetsView;
  const attachLib = globalThis.UnabatedAttach;
  const live = globalThis.UnabatedLive;
  const teaserLib = globalThis.UnabatedTeaser;
  const bet105 = globalThis.UnabatedBet105;
  // Bets service poll cadence while the panel is visible (plan § Storage).
  const BETS_POLL_MS = 30 * 1000;
  // The Edges list's defaults, book list and row logic (selection, sizing,
  // tail-flex rank, sort, cards, edge-move reading) live in edgerows.js,
  // shared with the headless server runner so both list the same lines at
  // the same stakes.
  const edgeRows = globalThis.UnabatedEdgeRows;
  const DEFAULT_SETTINGS = edgeRows.DEFAULT_STAKE_SETTINGS;
  const DEFAULT_EDGE_SETTINGS = edgeRows.DEFAULT_EDGE_SETTINGS;
  // Off until the list has been watched for a session (plan, 2026-09-10).
  const DEFAULT_ALERT_SETTINGS = { enabled: false, minEdgePct: 2.0 };
  const ALERT_EVENT_COOLDOWN_MS = 5 * 60 * 1000;
  const ALERT_LOG_TTL_MS = 24 * 60 * 60 * 1000;
  const ALERT_TARGETS_KEPT = 50;
  // Watcher heartbeats every 5s; past this with no heartbeat, the Unabated tab is gone.
  const WATCH_STALE_MS = 15000;
  // page.js republishes the books filter every 10s while an Unabated tab is open.
  const BOOKS_FILTER_STALE_MS = 6 * 60 * 60 * 1000;
  const MAX_EDGE_ROWS = edgeRows.MAX_EDGE_ROWS;
  // Open bets never net against a live line: conditional Kelly needs a live
  // fair ladder, which the screen does not expose (plan 2026-09-27, § Sizing).
  const LIVE_NOT_SIZED_NOTE = "live · not sized";

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
    betsCount: el("bets-count"), betsAlert: el("bets-alert"), betsTabButton: el("bets-tab-button"),
    betsRisk: el("bets-risk"), betsRiskCaption: el("bets-risk-caption"),
    edgesError: el("edges-error"), edgesStatus: el("edges-status"), edgesFilter: el("edges-filter"), edgesTailFlex: el("edges-tailflex"), edgesFilterDebug: el("edges-filter-debug"), edgesLocate: el("edges-locate"),
    edgesSports: el("edges-sports"), edgesBetTypes: el("edges-bettypes"), edgesBooks: el("edges-books"), edgesBooksMode: el("edges-books-mode"),
    booksDefault: el("books-default"), booksUnabated: el("books-unabated"), booksAll: el("books-all"), booksNone: el("books-none"), edgesPeriods: el("edges-periods"), edgesMin: el("edges-min"), edgesMinStake: el("edges-min-stake"), edgesMaxAge: el("edges-max-age"), edgesSort: el("edges-sort"),
    edgesIncludeAlts: el("edges-include-alts"), edgesMinToWin: el("edges-min-to-win"), edgesGroup: el("edges-group"),
    edgesSettingsError: el("edges-settings-error"), edgesList: el("edges-list"), edgesEmpty: el("edges-empty"),
    edgesLive: el("edges-live"), edgesLiveHead: el("edges-live-head"), edgesLiveList: el("edges-live-list"),
    edgesLiveCount: el("edges-live-count"), edgesPregameLabel: el("edges-pregame-label"),
    alertsEnabled: el("alerts-enabled"), alertsMin: el("alerts-min"),
    betsHeader: el("bets-header"), betsBanner: el("bets-banner"),
    tabBets: el("tab-bets"), betsService: el("bets-service"), betsSources: el("bets-sources"),
    betsUrl: el("bets-url"), betsSettingsError: el("bets-settings-error"),
    betsOpen: el("bets-open"), betsOpenCount: el("bets-open-count"), betsOpenEmpty: el("bets-open-empty"), betsMatchNote: el("bets-match-note"),
    betsNeedsBanner: el("bets-needs-banner"), betsNeedsBlock: el("bets-needs-block"), betsNeeds: el("bets-needs"), betsNeedsCount: el("bets-needs-count"),
    betsFixBanner: el("bets-fix-banner"), betsFixBlock: el("bets-fix-block"), betsFix: el("bets-fix"), betsFixCount: el("bets-fix-count"),
    betsOffboard: el("bets-offboard"), betsOffboardList: el("bets-offboard-list"), betsOffboardCount: el("bets-offboard-count"),
    betsCrosswalk: el("bets-crosswalk"), betsCrosswalkCount: el("bets-crosswalk-count"), betsCrosswalkEmpty: el("bets-crosswalk-empty"),
    betsCrosswalkClear: el("bets-crosswalk-clear"),
    shapeBanner: el("shape-banner"),
    tabTeasers: el("tab-teasers"), teasersCount: el("teasers-count"),
    teasersService: el("teasers-service"), teasersWarning: el("teasers-warning"), teasersStatus: el("teasers-status"), teasersEmpty: el("teasers-empty"),
    teasersSummary: el("teasers-summary"), teasersSummaryLabel: el("teasers-summary-label"), teasersSummaryStake: el("teasers-summary-stake"),
    teasersSummaryCells: el("teasers-summary-cells"), teasersSummaryNote: el("teasers-summary-note"),
    teasersOpen: el("teasers-open"), teasersOpenCount: el("teasers-open-count"), teasersOpenNote: el("teasers-open-note"), teasersOpenList: el("teasers-open-list"),
    teasersListLabel: el("teasers-list-label"), teasersListCount: el("teasers-list-count"), teasersList: el("teasers-list"),
    teasersMore: el("teasers-more"), teasersMoreLabel: el("teasers-more-label"), teasersMoreList: el("teasers-more-list"),
    teasersLegsLabel: el("teasers-legs-label"), teasersLegsCount: el("teasers-legs-count"), teasersLegs: el("teasers-legs"),
    teasersLegsMore: el("teasers-legs-more"), teasersLegsMoreLabel: el("teasers-legs-more-label"), teasersLegsMoreList: el("teasers-legs-more-list"),
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
    // Normalised bet records (bets.js contract), team keys resolved, pruned to the retention window.
    betRecords: [],
    // The team crosswalk the bets service holds (#118 step 4), as last served
    // or stored: [{venue, league, venueTeamKey, venueTeamName, unabatedTeamId,
    // unabatedTeamName, learnedFrom, learnedAt, pinnedBetId}]. Keys resolve through it first.
    crosswalk: [],
    // Cal's manual attaches as the bets service serves them: [{betId, venue,
    // league, eventId, eventStart, awayTeamId, homeTeamId, awayTeamName,
    // homeTeamName, pinnedAt}]. A pin decides its bet's game (bets.applyPins).
    pins: [],
    // Bet ids Cal dismissed from the red flag (the Dismiss chip): panel view
    // state, kept with the records in chrome.storage.local, pruned to open bets.
    dismissedBetIds: [],
    // {betId: startMs}: the start of the board event each open bet last
    // matched (bets.matchedStarts), so a bet with no start of its own reads as
    // started once its finished game leaves the board. Kept like dismissals.
    knownStarts: {},
    // The fair each bet's line had when it was placed, as the bets service
    // serves it (bets.duckdb::bet_fill_fairs): [{betId, lineKey, points,
    // fairAmerican, fairObservedAt, placedAt, capturedAt}]. Set through
    // setFillFairs, which keeps fillFairIndex in step.
    fillFairs: [],
  };
  // key -> alertLogEntry {lineKey, price, edgePct, stake, rankScore, at}: what has been alerted (or seen at baseline); persisted.
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
  // The scanner's per-line history (edgemove.js, #132): what each line was
  // worth on every snapshot, read for the edge-move tag.
  let scannerHistory = {};
  let boardLinesCache = null;
  // Unabated's fair ladders for sizing against held bets (#130):
  // edgeRows.createLadderReaders over the current feed state, built on first
  // use and dropped on every scanner update.
  let ladderReaders = null;
  // The Teasers tab: the ticket list as last built (teaser.planTeasers keeps
  // it while nothing real changes), and teaser.teaserBoardOf's one pass over
  // the board with the NFL/CFB load times it was made at — redone only when
  // one of those leagues' snapshots lands, not on every league's refresh.
  let teaserBuild = null;
  let teaserBoard = null;
  let teaserBoardLoadedAt = null;
  // The CFB markets Cal marked can't tease (teaser.marketKeyOf -> the game's
  // start in ms), kept in chrome.storage.local ("teaserBlocked"): Buckeye
  // keeps some college games off its teaser menu. A mark ends when its game starts.
  let teaserBlocked = {};
  // The tail-flex c per (league, period, bet type), measured off the
  // exchanges' two-sided rungs (tailflex.js); dropped on every scanner update
  // and when the max line age changes.
  const tailflex = globalThis.UnabatedTailFlex;
  let tailFlexCache = null;
  const teamsLib = globalThis.UnabatedTeams;
  let teamsSpellingCount = 0;
  const fillfair = globalThis.UnabatedFillFair;
  // bet id -> saved fill-fair row, over state.fillFairs.
  let fillFairIndex = new Map();

  function setFillFairs(rows) {
    state.fillFairs = rows;
    fillFairIndex = fillfair.fairsByBetId(rows);
  }

  // Every snapshot carries Unabated's team list and, per game row, a second
  // spelling of each team (feed.teamSpellingsFromEventName, #118): register
  // both as the team index (teams.js), persist it, and fill keys on bet
  // records that were waiting for it (#116 — no hand-written team tables).
  // The trigger counts spellings, not teams: a new eventName spelling for a
  // team already indexed must persist and re-resolve too.
  function registerFeedTeams(feedState) {
    if (!feedState || !feedState.teamIndex) return;
    for (const [league, list] of Object.entries(edgeRows.feedTeamsByLeague(feedState))) teamsLib.registerTeams(league, list);
    const spellings = teamsLib.spellingCount();
    if (spellings === teamsSpellingCount) return;
    teamsSpellingCount = spellings;
    chrome.storage.local.set({ teamsIndex: teamsLib.exportIndex() });
    state.betRecords = betsLib.resolveTeamKeys(state.betRecords, state.crosswalk);
  }

  const scanner = globalThis.UnabatedScanner.createScanner({
    onChange: (status, feedState, history) => {
      scannerStatus = status;
      scannerState = feedState;
      scannerHistory = history || {};
      boardLinesCache = null;
      ladderReaders = null;
      tailFlexCache = null;
      registerFeedTeams(feedState);
      if (noteMatchedStarts()) persistBets().catch((error) => console.error("[unabated-ticket] bets persist failed", error));
      learnCrosswalk().catch((error) => console.error("[unabated-ticket] crosswalk learn failed", error));
      captureFillFairs().catch((error) => console.error("[unabated-ticket] fill fair capture failed", error));
      renderEdges();
      renderTeasers();
      // A ticket sized from the feed (or waiting for it) follows the feed's
      // updates; one the screen priced is left alone (a re-render clears the copy status).
      if (state.ticket && !state.error && pricedLine(state.ticket).edgeFrom !== "screen") render();
      processAlerts().catch((error) => console.error("[unabated-ticket] alerts failed", error));
    },
  });

  // ---- formatting ----------------------------------------------------------

  const fmtAmerican = edgeRows.fmtAmerican;

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

  const fmtDollars = edgeRows.fmtDollars;

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
  // sitting on its main line's points. Two means two MARKETS that only the
  // marketId the ticket may lack tells apart, so the caller refuses.
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
    if (betsView.sourcesUnavailable(state.betsService && state.betsService.payload, Date.now())) {
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
    // A live capture names the break its fair was set at instead of a kickoff long past.
    view.startLine.textContent = ticket.live ? `live · ${live.checkpointLabel(ticket.live)}` : fmtStart(ticket.eventStart);
    const betFlag = ticketBetFlag(ticket, line);
    renderBetBanner(betFlag);

    view.book.textContent = ticket.book.name;
    view.price.textContent = fmtPriceBoth(asBookLine(line.price, line.sourceFormat, line.sourcePrice));
    // The fair is Unabated's own American number; there is no more exact source for it.
    view.fair.textContent = line.fair == null ? "unknown" : fmtPriceBoth(asBookLine(line.fair, 1, null));
    if (ticket.live && line.fair != null) view.fair.append(` · Unabated live fair, ${live.checkpointLabel(ticket.live)}`);
    // The edge-move tag reads the feed's copy of this line, and only at the
    // price being sized: the feed's history says nothing about another price.
    const feedCopy = feedLineFor(ticket, line.points);
    const ticketMoveParts = feedCopy && feedCopy.price === line.price ? moveParts(feedCopy, betFlag) : [];
    view.edge.replaceChildren(line.edgePct == null ? "\u2014" : fmtPct(line.edgePct / 100), ...ticketMoveParts.flatMap((part) => [" ", part]));

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

  // Under the dollar figure, the order it means on an exchange: "1,127
  // contracts @ 23.2¢ · $261.48", sized straight off Unabated's price for the
  // book with the count floored so the cost never passes the stake
  // (kelly.contractOrder). Only a line priced in contracts (Kalshi, Novig)
  // gets the row; a sportsbook line keeps just the dollars. `acted` is the
  // number to act on, so a top-up shows the top-up's contracts. Returns what
  // Copy appends, "" when there is no row.
  function renderContracts(acted, line) {
    view.contracts.hidden = true;
    view.contracts.classList.remove("under");
    view.contracts.replaceChildren();
    if (acted == null || acted <= 0) return "";
    const order = kelly.contractOrder({ stake: acted, bookPrice: line.price, sourceFormat: line.sourceFormat, sourcePrice: line.sourcePrice });
    if (!order) return "";
    view.contracts.hidden = false;
    const priceText = `${order.priceCents.toFixed(1)}\u00a2`;
    if (order.contracts === 0) {
      view.contracts.classList.add("under");
      view.contracts.textContent = `under 1 contract @ ${priceText}`;
      return ` | under 1 contract @ ${priceText}`;
    }
    const count = document.createElement("span");
    count.textContent = `${order.contracts.toLocaleString("en-US")} contract${order.contracts === 1 ? "" : "s"} @ ${priceText}`;
    const cost = document.createElement("span");
    cost.className = "cost";
    cost.textContent = ` \u00b7 ${fmtDollars(order.costDollars)}`;
    view.contracts.append(count, cost);
    return ` | ${order.contracts} contract${order.contracts === 1 ? "" : "s"} @ ${priceText}`;
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
    // Liquidity is for one price, like the edge: only the feed's copy of the
    // line at the price being sized caps the stake, as it does on the Edges row.
    const feedCopy = feedLineFor(ticket, line.points);
    const liquidity = feedCopy && feedCopy.price === line.price ? feedCopy.liquidity : null;
    const advice = betsView.stakeAdvice({
      line: matchLine, price: line.price, edgePct: line.edgePct,
      bankroll: state.settings.bankroll, multiplier: state.settings.multiplier,
      matches, ladderOf: ladderReader(ticket.eventId), liquidity, teasers: openTeasersNow(),
    });
    return { tier: matches.length ? matches[0].tier : null, matches: advice.matches, advice };
  }

  // The open BFA teasers, each leg joined to its board game and priced
  // (teaser.openTeasers): the Edges rows, alerts and the Ticket size the next
  // bet against the ones on its game (teasers plan section 14). None before
  // the scanner holds a board; a leg whose game is not on it yet (CFB still
  // loading, a league that failed) leaves its ticket out (betsview.js) until
  // it is, unless BFA says that game has started.
  function openTeasersNow() {
    if (!scannerState) return [];
    return teaserLib.openTeasers(state.betRecords, boardLines(), { now: Date.now(), ladderOf: currentTeaserBoard().ladderOf });
  }

  // Every open bet on this line's game, bets in the math first: this line,
  // same side, the other side (red); then the ones that are not sized, grey,
  // with why. Nothing when none.
  function renderBetBanner(flag) {
    const { shown, more } = betsView.bannerLines(betsView.relatedLines(flag, fillFairIndex));
    const items = shown.map((related) => {
      const div = document.createElement("div");
      const against = related.tier === "opposite" || related.tier === "related_opposite";
      div.className = `bet-match tier-${related.tier}${!related.inMath ? " not-sized" : against ? " bad" : ""}`;
      if (related.title) div.title = related.title;
      const kind = document.createElement("span");
      kind.className = "k";
      kind.textContent = related.tag;
      div.append(kind, relatedText(related));
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

  // Under the stake, the same pieces as the Edges rail: the position in
  // dollars, the number to act on (the big figure, its label says "Bet" or
  // "Add to your position"), that liquidity capped it, and what the stake
  // would be alone.
  function renderStakeExposure(advice) {
    const words = betsView.stakeAdviceWords(advice);
    const teasers = advice.teasers || { held: 0, against: 0 };
    const position = [
      advice.held > 0 ? `held ${betsLib.formatStake(advice.held)}` : null,
      advice.against > 0 ? `against ${betsLib.formatStake(advice.against)}` : null,
      teasers.held > 0 ? `teasers ${betsLib.formatStake(teasers.held)}` : null,
      teasers.against > 0 ? `teasers against ${betsLib.formatStake(teasers.against)}` : null,
      words ? words.cap : null,
      words ? words.alone : null,
    ].filter(Boolean);
    const onlyAgainst = advice.held === 0 && teasers.held === 0 && (advice.against > 0 || teasers.against > 0);
    view.stakeExposure.classList.toggle("against", onlyAgainst);
    view.stakeExposure.hidden = position.length === 0;
    view.stakeExposure.textContent = position.join(" \u00b7 ");
    view.stakeLabel.textContent = !words ? "Bet"
      : advice.bet === 0 && advice.cappedAt === 0 ? "Nothing resting at this price"
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

  // A live line's last price change: seconds matter between quarters.
  function fmtChangedAgo(modifiedMs) {
    if (modifiedMs == null) return null;
    const seconds = Math.max(0, Math.round((Date.now() - modifiedMs) / 1000));
    if (seconds < 90) return `changed ${seconds}s ago`;
    return `changed ${Math.round(seconds / 60)}m ago`;
  }

  function fmtLiquidity(value) {
    return value == null ? "" : `liq ${value.toLocaleString("en-US", { style: "currency", currency: "USD", maximumFractionDigits: 0 })}`;
  }

  // ---- why the edge grew (#132) --------------------------------------------

  const edgemove = globalThis.UnabatedEdgeMove;
  const MOVE_TAG_CLASS = { fair_to_you: "move-fair", book_away: "move-book", fair_against: "move-against" };

  const fmtFairEntry = edgeRows.fmtFairEntry;
  const fmtBetPrice = edgeRows.fmtBetPrice;

  // What the edge-move reading needs from the panel: the scanner's history,
  // the saved fill fairs and the clock.
  function moveContext() {
    return { history: scannerHistory, fillFairIndex, now: Date.now() };
  }

  // "fair 33.7% → 35.6% · price +199 · 33.4¢ → +215 · 31.7¢": a move's two ends.
  function fmtMoveNumbers(move) {
    const price = (entry) => (typeof entry.price === "number" ? fmtPriceBoth(asBookLine(entry.price, entry.sourceFormat, entry.sourcePrice)) : "?");
    return `fair ${fmtFairEntry(move.from)} \u2192 ${fmtFairEntry(move.to)} \u00b7 price ${price(move.from)} \u2192 ${price(move.to)}`;
  }

  // "opened +185", or "opened -120 at -3" when the book has moved its number
  // since (#126: a price on another number is not comparable); null when the
  // line carries no opener (every alt rung).
  function openerText(line) {
    const openerPrice = line.openerPrice ?? null;
    const openerPoints = line.openerPoints ?? null;
    if (openerPrice == null) return null;
    const sameNumber = openerPoints == null || openerPoints === line.points;
    return `opened ${fmtAmerican(openerPrice)}${sameNumber ? "" : ` at ${fmtPoints(openerPoints)}`}`;
  }

  // The last ten minutes, for the since-fill tooltip: the tag it would show
  // on its own and its numbers. An amber there is the fair lag itself: the
  // book moved and Unabated's fair (~1-2 min behind, #126) has not answered.
  function recentMoveText(move) {
    if (move.kind === "none") return "last 10 min: nothing moved";
    const lag = move.kind === "book_away" ? " \u2014 Unabated's fair runs ~1-2 min behind the book; wait a snapshot" : "";
    return `last 10 min: ${edgemove.MOVE_LABELS[move.kind]}, ${fmtMoveNumbers(move)}, moved ${fmtAge(move.sinceMs)}${lag}`;
  }

  // The numbers behind the tag. Since a fill: "since your Novig +208 bet
  // (Sep 23, 12:03 PM): fair 36.0% → 34.7% · price … | last 10 min: … |
  // opened +185". Ten minutes only: "fair 33.7% → 35.6% · price +199 → +215
  // · moved 2m ago · opened +185" — "ago" is when the panel first SAW the
  // move: a snapshot observation can be up to one refresh interval after the
  // book moved.
  function moveTooltip(reading, line) {
    const opener = openerText(line);
    if (!reading.baseline) {
      const move = reading.move;
      return [fmtMoveNumbers(move), `moved ${fmtAge(move.sinceMs)}`, opener].filter(Boolean).join(" \u00b7 ");
    }
    const bet = reading.baseline.bet;
    const placed = new Date(reading.baseline.placedMs).toLocaleString([], { month: "short", day: "numeric", hour: "numeric", minute: "2-digit" });
    const numbers = reading.move.from.price == null
      ? `fair ${fmtFairEntry(reading.move.from)} \u2192 ${fmtFairEntry(reading.move.to)} \u00b7 price not compared (another book, or no fill price)`
      : fmtMoveNumbers(reading.move);
    const since = `since your ${betsLib.venueLabel(bet.venue)}${fmtBetPrice(bet)} bet (${placed}): ${numbers}`;
    return [since, recentMoveText(reading.recent), opener].filter(Boolean).join(" | ");
  }

  // One small tag naming why the edge on this line is what it is (the fair
  // decides — see edgemove.js), its detail line under it and the numbers in
  // the tooltip, for a rail or the Ticket's Edge fact; [] when the line is
  // not held in this direction or nothing moved. `bet` is the row's {tier,
  // matches, advice} (withBetFlags / ticketBetFlag).
  function moveParts(line, bet) {
    const moveTag = edgeRows.moveTag(line, bet, moveContext());
    if (!moveTag) return [];
    const tag = document.createElement("span");
    tag.className = `tag ${MOVE_TAG_CLASS[moveTag.kind]}`;
    tag.textContent = moveTag.label;
    tag.title = moveTooltip(moveTag.reading, line);
    const detail = document.createElement("small");
    detail.className = "move-detail";
    detail.textContent = moveTag.detail;
    detail.title = tag.title;
    return [tag, detail];
  }

  // The tag's words for an alert body, or null; alert rows come through
  // withBetFlags, so row.bet is set.
  function moveWords(row) {
    return edgeRows.moveWords(row, moveContext());
  }

  function pageScriptAlive() {
    const ready = state.pageReady;
    return Boolean(ready && Date.now() - ready.at < PAGE_READY_STALE_MS);
  }

  function liveBooks() {
    return edgeRows.liveBooks(scannerState);
  }

  function unabatedSelection() {
    const filter = state.booksFilter;
    const fresh = filter && typeof filter.at === "number" && Date.now() - filter.at < BOOKS_FILTER_STALE_MS;
    return fresh && Array.isArray(filter.bookIds) && filter.bookIds.length ? filter.bookIds : null;
  }

  // Which books the list is restricted to (edgeRows.effectiveFilter): the
  // user's ticks, the default books, the Unabated selection page.js read off
  // the open tab, or every live book.
  function effectiveFilter() {
    return edgeRows.effectiveFilter(state.edgeSettings, scannerState, unabatedSelection(), state.booksFilter);
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
    if (state.edgeSettings.minLiquidityToWin > 0) parts.push(`exchange liq wins ≥ ${fmtDollars(state.edgeSettings.minLiquidityToWin)}`);
    parts.push(describeAltFilter(state.edgeSettings));
    return parts.join(" · ");
  }

  // "tail flex: NFL spr 6.2% · CFB 1H tot 10%": the c in use for every
  // spread/total market on the list, in list order; a market too thin to
  // measure shows the fallback.
  function describeTailFlex(rows) {
    if (!scannerState || !rows.length) return "";
    return edgeRows.describeTailFlex(rows, tailFlexMeasurement());
  }

  function describeAltFilter(settings) {
    if (!settings.includeAlts) return "alts: off";
    const minPct = Math.round(feed.ALT_MIN_PROB * 100);
    return `alts: on (fair and price ${minPct}-${100 - minPct}%)`;
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

  // Quarter-Kelly with nothing held, never more than the line's resting liquidity.
  function stakeFor(row) {
    return edgeRows.stakeFor(row, state.settings);
  }

  // Everything selectEdges needs except the edge threshold (list and alerts differ there).
  function edgeSelectionOptions(effective) {
    return edgeRows.edgeSelectionOptions(state.edgeSettings, effective, Date.now());
  }

  // (period, axis) -> Unabated's fair ladder for one event, kept until the next scanner update.
  function ladderReader(eventId) {
    if (!ladderReaders) ladderReaders = edgeRows.createLadderReaders(scannerState);
    return ladderReaders(eventId);
  }

  // What conditional Kelly sizes a row against: the open bets, the board, the
  // bankroll, the fair ladders and the open BFA teasers.
  function betFlagContext() {
    return { records: state.betRecords, boardLines: boardLines(), stakeSettings: state.settings, ladderReaderOf: ladderReader, teasers: openTeasersNow() };
  }

  // Each row gets `bet` = {tier, matches, advice} (edgeRows.withBetFlags).
  function withBetFlags(rows) {
    return edgeRows.withBetFlags(rows, betFlagContext());
  }

  function meetsMinStake(row) {
    return edgeRows.meetsMinStake(row, state.edgeSettings.minStake);
  }

  // Measured once per scanner update, on the same line-age gate as the list.
  function tailFlexMeasurement() {
    const maxLineAgeMs = state.edgeSettings.maxLineAgeHours * 3600 * 1000;
    if (!tailFlexCache || tailFlexCache.maxLineAgeMs !== maxLineAgeMs) {
      tailFlexCache = { maxLineAgeMs, measurement: tailflex.measureTailFlex(scannerState, { now: Date.now(), maxLineAgeMs }) };
    }
    return tailFlexCache.measurement;
  }

  // The standalone stake and the tail-flex rank score (edgeRows.withStakeAndRank).
  function withStakeAndRank(row, measurement) {
    return edgeRows.withStakeAndRank(row, measurement, state.settings);
  }

  function currentEdgeRows() {
    if (!scannerState) return [];
    return edgeRows.listedEdgeRows(scannerState, {
      ...betFlagContext(), edgeSettings: state.edgeSettings, effective: effectiveFilter(), measurement: tailFlexMeasurement(), now: Date.now(),
    });
  }

  // ---- live edges (the Unabated live screen, between quarters) -------------

  // league path -> the latest live_edges payload page.js posted for it.
  const liveStore = {};

  function livePayloads() {
    return Object.values(liveStore);
  }

  // Held bets on the game show on a live row, but none is netted: the stake
  // is the line's own quarter-Kelly (LIVE_NOT_SIZED_NOTE says why).
  function withLiveBetFlags(rows) {
    let flags = null;
    try {
      flags = betsLib.annotateRows(rows, state.betRecords, { lines: boardLines() });
    } catch (error) {
      console.info("[unabated-ticket] live bet flags skipped:", error.message);
    }
    return rows.map((row, index) => {
      const flag = flags ? flags[index] : { tier: null, matches: [] };
      const advice = betsView.stakeAdvice({
        line: row, price: row.price, edgePct: row.edgePct,
        bankroll: state.settings.bankroll, multiplier: state.settings.multiplier,
        matches: [], ladderOf: null, liquidity: row.liquidity,
      });
      const matches = (flag.matches || []).map((match) => ({ ...match, inMath: false, note: LIVE_NOT_SIZED_NOTE }));
      return { ...row, bet: { tier: flag.tier, matches, advice } };
    });
  }

  function liveRowsAt(minEdgePct) {
    const options = { ...edgeSelectionOptions(effectiveFilter()), minEdgePct, now: Date.now() };
    const selected = live.selectLiveEdges(livePayloads(), options).map((row) => ({ ...row, stake: stakeFor(row) }));
    return withLiveBetFlags(selected);
  }

  function currentLiveRows() {
    return liveRowsAt(state.edgeSettings.minEdgePct).filter(meetsMinStake);
  }

  function fmtAgoShort(ms) {
    const seconds = Math.max(0, Math.round(ms / 1000));
    return seconds < 90 ? `0:${String(seconds).padStart(2, "0")}` : `${Math.round(seconds / 60)}m`;
  }

  // The Live block's one-line header: what the screen is doing right now.
  function liveHeadText(liveState, rows) {
    if (liveState.stale) return `last read ${fmtAgoShort(liveState.readAgoMs)} ago · bring the Unabated live tab to the front`;
    const ready = liveState.games.filter((game) => game.fairReady);
    if (!ready.length) {
      const names = liveState.games.length === 1 ? liveState.games[0].matchup : `${liveState.games.length} games`;
      return `${names} · in play · edges return at the next break`;
    }
    const newestFairMs = ready.reduce((max, game) => {
      const produced = Date.parse(/[zZ]$/.test(game.producedUtc || "") ? game.producedUtc : `${game.producedUtc}Z`);
      return Number.isFinite(produced) ? Math.max(max, produced) : max;
    }, 0);
    const games = ready.map((game) => `${game.matchup} · ${game.checkpoint}`).join("; ");
    const fairAge = newestFairMs ? ` · fair set ${fmtAgoShort(Date.now() - newestFairMs)} ago` : "";
    const none = rows.length ? "" : ` · nothing at ${state.edgeSettings.minEdgePct}% or more`;
    return `${games}${fairAge}${none}`;
  }

  // Rendered on every Edges render and once a second while a live game is on
  // a screen, so the staleness clock moves on its own.
  function renderLive() {
    const liveState = live.liveView(livePayloads(), Date.now());
    const shown = liveState.games.length > 0;
    view.edgesLive.hidden = !shown;
    view.edgesPregameLabel.hidden = !shown;
    const rows = shown ? currentLiveRows() : [];
    view.edgesLiveCount.hidden = rows.length === 0 || liveState.stale;
    view.edgesLiveCount.textContent = `${rows.length} live`;
    if (!shown) {
      view.edgesLiveList.replaceChildren();
      return;
    }
    const anyReady = liveState.games.some((game) => game.fairReady);
    view.edgesLiveHead.className = `live-head${liveState.stale ? " stale" : anyReady ? "" : " inplay"}`;
    const dot = document.createElement("span");
    dot.className = "dot";
    const label = document.createElement("b");
    label.textContent = "Live";
    const text = document.createElement("span");
    text.append(label, ` · ${liveHeadText(liveState, rows)}`);
    view.edgesLiveHead.replaceChildren(dot, text);
    view.edgesLiveList.replaceChildren(...rows.slice(0, MAX_EDGE_ROWS).map((row) => renderEdgeRow({ ...row, liveStale: liveState.stale })));
  }

  // Cards: the best line of each (game, market, side) is always the highest
  // tail-flex rank score (EV dollars after flex); the panel's sort orders the
  // cards through that line.
  function groupsOf(rows) {
    return edgeRows.groupEdgeRows(rows, state.edgeSettings.sortBy);
  }

  // Cards the user has opened; survives the 5s re-render, not a panel reload.
  const expandedGroups = new Set();

  // "held $300" and/or "against $200", "teasers $600", or "game"; the tooltip
  // is every match's label, a teasers chip's its tickets.
  function betBadges(flag) {
    return betsView.badges(flag).map(({ kind, text, title }) => {
      const badge = document.createElement("span");
      badge.className = `tag ${kind}`;
      badge.textContent = text;
      badge.title = title || flag.matches.map((match) => match.label).join("\n");
      return badge;
    });
  }

  // Edge magnitude in three steps (hot / warm / thin): the row's stripe and figure colour.
  const edgeTier = edgeRows.edgeTier;

  // The rail under the edge (edgeRows.stakeRail): the number to act on with
  // its verb, then one small line — all the liquidity there is, and what the
  // stake would be with nothing held.
  function fillStakeCell(cell, row) {
    // A live row from a tab that stopped reading: the number may be gone.
    if (row.liveStale) {
      cell.textContent = "—";
      return;
    }
    const rail = edgeRows.stakeRail(row);
    cell.classList.toggle("at-size", rail.atSize);
    cell.append(rail.text);
    if (!rail.note) return;
    const note = document.createElement("small");
    note.textContent = rail.note;
    cell.append(" ", note);
  }

  // The bets already on this game, as a labelled section of the row rather
  // than a loose line. Every line names the bet and how it relates to this
  // line; past three the rest are a count, as the Ticket banner does it.
  const RELATED_LINES_ON_A_ROW = 3;

  // The bet's own words, then "fair then 36.0%" when its fill fair was saved.
  function relatedText(related) {
    const text = document.createElement("span");
    text.append(related.text);
    if (related.fairThen) {
      const fairThen = document.createElement("span");
      fairThen.className = "fair-then";
      fairThen.textContent = related.fairThen;
      text.append(" \u00b7 ", fairThen);
    }
    return text;
  }

  function relatedBlock(flag) {
    const all = betsView.relatedLines(flag, fillFairIndex);
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
      // A teasers line lists its tickets on hover.
      if (line.title) div.title = line.title;
      const tag = document.createElement("span");
      tag.className = "related-tag";
      tag.textContent = line.tag;
      div.append(tag, relatedText(line));
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
    li.className = `edge-row tier-${tier}${row.isBlurred ? " blurred" : ""}${row.liveStale ? " stale" : ""}${row.key === lastClickedKey ? " last-clicked" : ""}`;
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
    if (row.live) {
      const liveTag = document.createElement("span");
      liveTag.className = "tag live";
      liveTag.textContent = "live";
      liveTag.title = `Priced off Unabated's in-game fair at the ${row.live.checkpoint}`;
      side.append(liveTag);
    }
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
    const marketText = `${row.betType}${row.period === "FG" ? "" : ` · ${row.period}`}${row.isAlt ? ` · alt of ${fmtPoints(row.mainPoints)}` : ""}`
      + `${row.rotation != null ? ` · rot ${row.rotation}` : ""} · `;
    if (row.live) meta.append(marketText, `${row.live.matchup} · ${row.live.checkpoint}`);
    else meta.append(marketText, `${describeMatchup(row)} · ${fmtStart(row.eventStart)} · `, untilEl(row.eventStartMs));

    const book = document.createElement("div");
    book.className = "edge-book";
    const price = document.createElement("span");
    price.className = "price";
    price.textContent = `${row.book.name} ${fmtPriceBoth(asBookLine(row.price, row.sourceFormat, row.sourcePrice))}`;
    const age = document.createElement("span");
    age.className = "age";
    const ageParts = row.live
      ? [`fair ${fmtAmerican(row.fair)}`, fmtChangedAgo(row.modifiedMs), fmtLiquidity(row.liquidity)]
      : [fmtLineAge(row.modifiedMs), fmtLiquidity(row.liquidity)];
    age.textContent = ` · ${ageParts.filter(Boolean).join(" · ")}`;
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
    // The edge-move tag reads the pregame history; a live line has none.
    rail.append(pct, ...(row.live ? [] : moveParts(row, row.bet)), stake);

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
    rail.append(edge, ...moveParts(row, row.bet), stake);

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
    const loading = status.loading ? ` · loading ${status.loading.done}/${status.loading.total} leagues` : "";
    view.edgesStatus.textContent = status.phase === "loading" && !status.leaguesLoaded.length
      ? `Loading snapshots…${loading}`
      : `${leagues.join(" · ") || "no leagues"} · ${status.lineCount.toLocaleString()} lines${status.altLineCount ? ` (+${status.altLineCount.toLocaleString()} alts)` : ""} · ${built}${loading}`;
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
    renderLive();
    const rows = currentEdgeRows();
    const grouped = state.edgeSettings.groupByMarket;
    const items = grouped ? groupsOf(rows) : rows;
    renderEdgesStatus(items);
    const effective = effectiveFilter();
    view.edgesFilter.textContent = describeFilter(effective);
    view.edgesTailFlex.textContent = describeTailFlex(rows);
    view.edgesTailFlex.hidden = view.edgesTailFlex.textContent === "";
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
      // Only a live row says so: page.js then looks among the live grid rows.
      ...(row.live ? { live: true } : {}),
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

  view.edgesLiveList.addEventListener("click", async (event) => {
    const li = event.target.closest("li.edge-row");
    if (!li) return;
    const row = currentLiveRows().find((r) => r.key === li.dataset.key);
    if (!row) return;
    lastClickedKey = row.key;
    const request = locateRequestOf(row);
    state.locating = { ...request, at: Date.now() };
    state.locateResult = null;
    renderLocate();
    renderLive();
    try {
      await globalThis.UnabatedLocate.locateLine(request);
    } catch (error) {
      state.locating = null;
      state.locateResult = { ...request, at: Date.now(), ok: false, message: `could not focus an Unabated tab (${error.message})` };
      renderLocate();
    }
  });

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
      // own ranking, the rank score: the best rung pulled and a +944 longshot
      // taking over is a worse card, not news, whatever its edge %. The line
      // itself must have changed (another line, or its price or edge moved):
      // the score also moves with every re-measure of c, which is not news.
      improvedOn: (previous) => cardLineChanged(group.best, previous)
        && typeof previous.rankScore === "number" && typeof group.best.rankScore === "number" && group.best.rankScore > previous.rankScore,
    }));
  }

  function cardLineChanged(best, previous) {
    return best.key !== previous.lineKey || best.price !== previous.price || best.edgePct !== previous.edgePct;
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
      // Why it grew (#132): a re-fire on a better edge says which improvement it was.
      moveWords(row),
      stake == null ? null : `stake ${fmtDollars(stake)}`,
      summary,
      describeMatchup(row),
      row.live ? null : fmtUntil(row.eventStartMs),
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
    const measurement = tailFlexMeasurement();
    const selected = feed.selectEdges(scannerState, { ...edgeSelectionOptions(effectiveFilter()), minEdge: state.alertSettings.minEdgePct / 100 })
      .map((row) => withStakeAndRank(row, measurement));
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

  // Live edges alert once per line per break (live.liveAlertKey carries the
  // fair's production time); a tab that stopped reading alerts nothing.
  function liveAlertItems() {
    if (live.liveView(livePayloads(), Date.now()).stale) return [];
    return liveRowsAt(state.alertSettings.minEdgePct)
      .filter((row) => betsView.suggestedBetAmount(row.bet.advice) > 0 && meetsMinStake(row))
      .map((row) => ({ key: `live:${live.liveAlertKey(row)}`, row, summary: `live · ${row.live.checkpoint}`, improvedOn: () => false }));
  }

  // What the next pass compares against: the line, its price and edge, and the card's rank score.
  function alertLogEntry(row, now, baseline) {
    const entry = { lineKey: row.key, price: row.price, edgePct: row.edgePct, stake: row.stake, rankScore: row.rankScore ?? null, at: now };
    return baseline ? { ...entry, baseline: true } : entry;
  }

  async function processAlertsOnce() {
    const now = Date.now();
    const items = alertItems().concat(liveAlertItems());
    pruneAlertLog(now);
    if (!alertsBaselined) {
      for (const item of items) alertLog[item.key] = alertLogEntry(item.row, now, true);
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
      alertLog[item.key] = alertLogEntry(item.row, now, false);
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

  const TABS = ["ticket", "edges", "bets", "teasers"];
  const paneOf = { ticket: view.tabTicket, edges: view.tabEdges, bets: view.tabBets, teasers: view.tabTeasers };
  // Each pane scrolls on its own, but a hidden element is not guaranteed to
  // keep its scrollTop, so the position is remembered explicitly. Without this
  // the Edges list went back to the top every time a capture brought the
  // Ticket tab forward, which is the whole reason for the split panes.
  const paneScroll = { ticket: 0, edges: 0, bets: 0, teasers: 0 };

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
    if (state.activeTab === "teasers") renderTeasers();
  });

  // ---- edge settings -------------------------------------------------------

  function readEdgeSettingInputs() {
    const leagues = Array.from(view.edgesSports.querySelectorAll("input:checked")).flatMap((input) => feed.leagueIdsOfSport(input.dataset.sport));
    const periods = Array.from(view.edgesPeriods.querySelectorAll("input:checked")).map((input) => Number(input.dataset.period));
    const betTypes = Array.from(view.edgesBetTypes.querySelectorAll("input:checked")).map((input) => Number(input.dataset.bettype));
    const minEdgePct = Number(view.edgesMin.value);
    const minStake = Number(view.edgesMinStake.value);
    const maxLineAgeHours = Number(view.edgesMaxAge.value);
    const minLiquidityToWin = Number(view.edgesMinToWin.value);
    if (!Number.isFinite(minEdgePct) || minEdgePct < 0) return { error: "Minimum edge must be zero or more." };
    if (!Number.isFinite(minStake) || minStake < 0) return { error: "Minimum suggested bet must be zero (off) or more." };
    if (!Number.isFinite(maxLineAgeHours) || maxLineAgeHours <= 0) return { error: "Max line age must be above zero hours." };
    if (!Number.isFinite(minLiquidityToWin) || minLiquidityToWin < 0) return { error: "Min liq to win must be zero (off) or more." };
    if (!periods.length) return { error: "Pick at least one period." };
    if (!betTypes.length) return { error: "Pick at least one bet type." };
    return {
      settings: {
        ...state.edgeSettings, leagues, periods, betTypes, minEdgePct, minStake, maxLineAgeHours, minLiquidityToWin, sortBy: view.edgesSort.value,
        includeAlts: view.edgesIncludeAlts.checked,
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
    view.edgesMinToWin.value = settings.minLiquidityToWin;
    view.edgesSort.value = settings.sortBy;
    view.edgesIncludeAlts.checked = settings.includeAlts;
    view.edgesGroup.checked = settings.groupByMarket;
  }

  function onEdgeSettingsInput() {
    const parsed = readEdgeSettingInputs();
    view.edgesSettingsError.textContent = parsed.error || "";
    if (parsed.error) return;
    const before = state.edgeSettings;
    const leaguesChanged = parsed.settings.leagues.join(",") !== before.leagues.join(",");
    // Widening periods, bet types, the liquidity floor or the alts toggle exposes lines the alert log has never seen.
    const scopeChanged = leaguesChanged
      || parsed.settings.periods.join(",") !== before.periods.join(",")
      || parsed.settings.betTypes.join(",") !== before.betTypes.join(",")
      || parsed.settings.includeAlts !== before.includeAlts
      || parsed.settings.minLiquidityToWin !== before.minLiquidityToWin
      // Alert keys differ between the flat list and cards.
      || parsed.settings.groupByMarket !== before.groupByMarket;
    state.edgeSettings = parsed.settings;
    chrome.storage.local.set({ edges: parsed.settings });
    if (scopeChanged) alertsBaselined = false;
    // Ticking Football on or off leaves the scanner's leagues as they are: it loads football for the Teasers tab anyway.
    const scannerLeaguesChanged = scannerLeaguesOf(parsed.settings.leagues).join(",") !== scannerLeaguesOf(before.leagues).join(",");
    if (scannerLeaguesChanged) {
      scanner.start(scannerLeaguesOf(parsed.settings.leagues)).catch((error) => console.error("[unabated-ticket] scanner restart failed", error));
    }
    renderEdges();
    // Max line age decides which Buckeye legs can be teased.
    renderTeasers();
    // A restarted scanner runs the alert pass on its first snapshot.
    if (scopeChanged && !scannerLeaguesChanged) processAlerts().catch((error) => console.error("[unabated-ticket] alerts failed", error));
  }

  // The leagues the scanner loads: the Edges tab's, plus NFL and CFB for the Teasers tab.
  function scannerLeaguesOf(edgeLeagues) {
    return Array.from(new Set([...edgeLeagues, ...teaserLib.TEASER_LEAGUE_IDS])).sort((a, b) => a - b);
  }

  const sanitizeEdgeSettings = edgeRows.sanitizeEdgeSettings;

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
    boardLinesCache = edgeRows.boardLines(scannerState);
    return boardLinesCache;
  }

  function betsPayload() {
    return state.betsService ? state.betsService.payload : null;
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
      // Records merged, crosswalk / pins / fill fairs taken (edgeRows.applyBetsPayload); throws on a body with no bets.
      const held = { records: state.betRecords, crosswalk: state.crosswalk, pins: state.pins, fillFairs: state.fillFairs };
      const applied = edgeRows.applyBetsPayload(held, await response.json(), now);
      state.crosswalk = applied.crosswalk;
      state.pins = applied.pins;
      state.betRecords = applied.records;
      state.dismissedBetIds = betsView.keepDismissedOpen(state.dismissedBetIds, state.betRecords);
      noteMatchedStarts();
      if (applied.fillFairs !== state.fillFairs) setFillFairs(applied.fillFairs);
      state.betsService = {
        payload: { generatedAt: applied.generatedAt, sources: applied.sources },
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
    captureFillFairs().catch((error) => console.error("[unabated-ticket] fill fair capture failed", error));
    pollBet105().catch((error) => console.error("[unabated-ticket] bet105 poll failed", error));
  }

  // ---- Bet105 (read here, stored by the service) -----------------------------
  //
  // On the bets tick, at most every bet105.POLL_MS: the account's open bets
  // from both LinePros feeds, then one POST to the service. A read that fails
  // (not logged in, Cloudflare, the site down) is POSTed as an error so the
  // Bets tab's Bet105 row turns red with the fix; nothing half-read is ever
  // pushed (the service closes an open bet a complete push no longer lists).
  let bet105Busy = false;
  let bet105LastRunAt = 0;
  let bet105LastError = null;

  async function pollBet105() {
    if (bet105Busy || document.hidden || !bet105.isDue(bet105LastRunAt, Date.now())) return;
    bet105Busy = true;
    bet105LastRunAt = Date.now();
    try {
      const push = await readBet105();
      const reply = await serviceRequest("POST", bet105.SERVICE_PATH, push);
      if (bet105LastError !== null) console.info("[unabated-ticket] bet105: reading again", reply);
      bet105LastError = null;
    } catch (error) {
      const message = error.name === "TimeoutError" ? `Bet105 did not answer within ${bet105.FETCH_TIMEOUT_MS / 1000}s` : error.message;
      if (bet105LastError !== message) console.warn("[unabated-ticket] bet105:", message);
      bet105LastError = message;
      await serviceRequest("POST", bet105.SERVICE_PATH, bet105.errorBody(message)).catch(() => {});
    } finally {
      bet105Busy = false;
    }
  }

  // The session check (its reply carries the CSRF token the history POST
  // needs), then getHistory on each feed. Throws with the reason on any step.
  async function readBet105() {
    const customers = await fetch(bet105.CUSTOMERS_URL, bet105.customersRequest());
    const session = bet105.csrfTokenOf(customers.status, await customers.json().catch(() => null));
    if (session.error) throw new Error(session.error);
    const groupsByFeed = {};
    for (const feedName of bet105.FEEDS) {
      const response = await fetch(bet105.historyUrl(feedName), bet105.historyRequest(session.csrfToken));
      const result = bet105.betGroupsOf(feedName, response.status, await response.json().catch(() => null));
      if (result.error) throw new Error(result.error);
      groupsByFeed[feedName] = result.betGroups;
    }
    return bet105.pushBody(new Date().toISOString(), groupsByFeed);
  }

  // The records, the service state, the crosswalk, the pins and the saved fill fairs, as one stored object.
  function persistBets() {
    return chrome.storage.local.set({ betsService: { ...state.betsService, bets: state.betRecords, crosswalk: state.crosswalk, pins: state.pins, fillFairs: state.fillFairs, dismissed: state.dismissedBetIds, knownStarts: state.knownStarts } });
  }

  // Every surface that shows a bet flag, after the records or the crosswalk
  // changed — the Teasers tab too: open BFA teasers are its placed tickets.
  function renderBetsFlags() {
    renderBetsHeader();
    if (!state.error) render();
    renderEdges();
    if (state.activeTab === "bets") renderBets();
    renderTeasers();
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
      await applyServiceTables({ crosswalk: reply.crosswalk });
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
    const reply = await serviceRequest(method, "/crosswalk.json", body);
    if (!reply || !Array.isArray(reply.crosswalk)) throw new Error("crosswalk.json has no crosswalk array");
    return reply;
  }

  // One request to a bets-service write route (crosswalk, pins, fill fairs):
  // the parsed reply, or an Error carrying the HTTP status and the service's
  // own error text — none when the body is not JSON (an older service's 404).
  async function serviceRequest(method, path, body) {
    const response = await fetch(`${state.betsSettings.serviceUrl}${path}`, {
      method, cache: "no-store",
      headers: body ? { "Content-Type": "application/json" } : {},
      body: body ? JSON.stringify(body) : undefined,
    });
    const reply = await response.json().catch(() => null);
    if (!response.ok) {
      const error = new Error(`HTTP ${response.status}${reply && reply.error ? `: ${reply.error}` : ""}`);
      error.status = response.status;
      throw error;
    }
    return reply;
  }

  // The served crosswalk and/or pins replace the held ones; pins are applied
  // first (a pin can give a record its league) and keys rebuilt from scratch,
  // so a cleared row takes its key back and a new one applies everywhere.
  async function applyServiceTables(tables) {
    if (Array.isArray(tables.crosswalk)) state.crosswalk = tables.crosswalk;
    if (Array.isArray(tables.pins)) state.pins = tables.pins;
    state.betRecords = betsLib.rekeyRecords(betsLib.applyPins(state.betRecords, state.pins), state.crosswalk);
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
      await applyServiceTables({ crosswalk: [] });
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
    view.betsCrosswalkEmpty.textContent = "Nothing learned yet. Rows appear when an open bet joins the board by its Kalshi event or Novig outcome id, or when you attach one.";
  }

  // ---- fill fairs (since your first fill, 2026-09-23) -------------------------
  //
  // After every poll and every scanner update, each new open bet gets the
  // fair its line had when it was placed, read off the scanner's history
  // (fillfair.captureFillFairs) and POSTed to the bets service, which keeps
  // the first one per bet for good (bets.duckdb::bet_fill_fairs). A bet is
  // looked at until it is decided — saved, or refused for a reason that
  // cannot change (placed before the panel was watching, its line first seen
  // after it) — once per panel session. Rows go one per request: the service
  // refuses a whole request on its first bad row, and a good row must not go
  // down with it. A 400 drops that row; any other failure keeps every unsent
  // row and retries a minute later — they are the fill's numbers and do not age.
  const FILL_FAIRS_RETRY_MS = 60 * 1000;
  const HTTP_BAD_REQUEST = 400;
  let fillFairsBusy = false;
  let fillFairsRetryAt = 0;
  let fillFairsLastError = null;
  // Bet ids saved or refused this session; captured rows the service has not stored yet.
  const fillFairsDecided = new Set();
  const fillFairsUnsent = new Map();

  async function captureFillFairs() {
    if (!scannerState) return;
    const liveStatus = scanner.getStatus();
    const { saves, refusals } = fillfair.captureFillFairs({
      records: state.betRecords, skipIds: new Set([...fillFairIndex.keys(), ...fillFairsDecided]),
      state: scannerState, boardLines: boardLines(), history: scannerHistory,
      // Live, not the last onChange copy: pause() clears it before any notify.
      observingSince: liveStatus.observingSince, leagueObservingSince: liveStatus.leagueObservingSince, now: Date.now(),
    });
    for (const save of saves) {
      fillFairsDecided.add(save.betId);
      fillFairsUnsent.set(save.betId, save);
    }
    for (const refusal of refusals) fillFairsDecided.add(refusal.betId);
    logFillFairDecisions(saves, refusals);
    await sendFillFairs();
  }

  // One line per pass that decided anything: "2 captured; not saved: 14
  // placed before the panel was watching".
  function logFillFairDecisions(saves, refusals) {
    if (!saves.length && !refusals.length) return;
    const counts = new Map();
    for (const { reason } of refusals) counts.set(reason, (counts.get(reason) || 0) + 1);
    const refused = Array.from(counts, ([reason, count]) => `${count} ${reason}`).join(", ");
    console.info(`[unabated-ticket] fill fairs: ${saves.length} captured${refused ? `; not saved: ${refused}` : ""}`);
  }

  async function sendFillFairs() {
    if (fillFairsBusy || fillFairsUnsent.size === 0 || Date.now() < fillFairsRetryAt) return;
    fillFairsBusy = true;
    let served = null;
    try {
      for (const row of Array.from(fillFairsUnsent.values())) {
        const outcome = await sendFillFairRow(row);
        if (outcome.retryLater) break;
        if (outcome.served) served = outcome.served;
      }
    } finally {
      fillFairsBusy = false;
    }
    if (!served) return;
    setFillFairs(fillfair.mergeFillFairs(state.fillFairs, served, state.betRecords));
    await persistBets();
    renderBetsFlags();
  }

  // One row: {served} (the reply's rows) once the service holds it; {} when
  // the service refused it (a 400: dropped and logged); {retryLater: true}
  // when the service failed or could not be reached.
  async function sendFillFairRow(row) {
    try {
      const reply = await serviceRequest("POST", "/fill_fairs.json", { rows: [row] });
      if (!reply || !Array.isArray(reply.fillFairs)) throw new Error("fill_fairs.json has no fillFairs array");
      fillFairsUnsent.delete(row.betId);
      fillFairsLastError = null;
      console.info(`[unabated-ticket] fill fairs: ${row.betId} ${reply.saved ? "saved" : "already held by the service"}`);
      return { served: reply.fillFairs };
    } catch (error) {
      if (error.status === HTTP_BAD_REQUEST) {
        fillFairsUnsent.delete(row.betId);
        console.error(`[unabated-ticket] fill fairs: the service refused ${row.betId}, dropped:`, error.message, row);
        return {};
      }
      fillFairsRetryAt = Date.now() + FILL_FAIRS_RETRY_MS;
      if (fillFairsLastError !== error.message) console.warn("[unabated-ticket] fill fairs: service write failed:", error.message);
      fillFairsLastError = error.message;
      return { retryLater: true };
    }
  }

  function startBetsPolling() {
    if (betsPollTimer != null) return;
    betsPollTimer = setInterval(() => pollBets().catch((error) => console.error("[unabated-ticket] bets poll failed", error)), BETS_POLL_MS);
    pollBets().catch((error) => console.error("[unabated-ticket] bets poll failed", error));
  }

  // Every open bet no board line matches, with why and whether it needs a
  // game or a code fix (the red flag). Before the first snapshot the board
  // is empty and nothing is attachable, so until the board is known only a
  // bet its source could not read flags.
  function unmatchedNow() {
    return betsLib.unmatchedReasons(state.betRecords, boardLines(), Date.now(),
      { dismissedIds: state.dismissedBetIds, knownStarts: state.knownStarts });
  }

  // Remember the start of the board event each open bet matches now; true
  // when that changed what is held (the caller persists it).
  function noteMatchedStarts() {
    const before = JSON.stringify(state.knownStarts);
    state.knownStarts = betsView.keepKnownStartsOpen(state.knownStarts, betsLib.matchedStarts(state.betRecords, boardLines()), state.betRecords);
    return JSON.stringify(state.knownStarts) !== before;
  }

  // Dismiss stops a bet flagging red; Restore flags it again. Both are panel
  // view state (persisted with the records), so they take effect at once.
  async function setDismissed(betId, dismissed) {
    const others = state.dismissedBetIds.filter((id) => id !== betId);
    state.dismissedBetIds = dismissed ? [...others, betId] : others;
    await persistBets();
    renderBetsFlags();
  }

  function renderBetsHeader() {
    const now = Date.now();
    const open = state.betRecords.filter((bet) => bet.status === "open").length;
    const unmatched = unmatchedNow();
    const needsGame = unmatched.filter((entry) => entry.needsGame).length;
    const needsFix = unmatched.filter((entry) => entry.needsFix).length;
    const flagged = needsGame + needsFix;
    view.betsCount.hidden = open === 0;
    view.betsCount.textContent = String(open);
    view.betsAlert.hidden = flagged === 0;
    view.betsAlert.textContent = String(flagged);
    view.betsAlert.title = [
      needsGame ? `${needsGame} open bet${needsGame === 1 ? "" : "s"} not matched to a game` : null,
      needsFix ? `${needsFix} open bet${needsFix === 1 ? "" : "s"} needing a code fix` : null,
    ].filter(Boolean).join(" · ");
    view.betsTabButton.classList.toggle("flagged", flagged > 0);
    view.betsHeader.textContent = betsView.headerLine(state.betRecords, betsPayload(), now, needsGame, needsFix);
    view.betsHeader.classList.toggle("bad", flagged > 0 || betsView.serviceStatus(state.betsService, now).unreachable);
  }

  // One line per venue: a dot for freshness, what it holds, how old the last
  // successful pull is. A venue that failed says why, in place of its count.
  function renderBetsSources(now) {
    const service = betsView.serviceStatus(state.betsService, now);
    view.betsService.hidden = !service.unreachable;
    view.betsService.textContent = service.unreachable ? `${service.text}. Start it with unabated_ticket/bets_service/run.sh; the last records it served are still shown.` : "";
    view.betsSources.replaceChildren(...betsView.sourceRows(betsPayload(), now).map((row) => {
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

  // ---- Bets tab rows and the Attach flow ------------------------------------
  //
  // An open bet no board line matches is listed twice: in the open list with
  // a red edge, and in "Needs a game" or "Needs a code fix" (the red flag,
  // bets.unmatchedReasons needsGame / needsFix) or the folded "Not on the
  // board" with its reason. A code fix has no Attach. Attach opens
  // a two-step panel under the listed row (attach.js): pick the game, then
  // confirm what the venue's names mean. The attach is POSTed to the bets
  // service (/pins.json), whose reply — every pin and the whole crosswalk —
  // replaces the held tables. The open panel lives in `attachState`, so the
  // 30 s poll's re-render redraws it, search focus included.
  let attachState = null;

  function openAttach(betId) {
    attachState = { betId, step: "pick", query: "", eventId: null, swapped: false, busy: false, error: null, focusSearch: true };
    renderBets();
  }

  function closeAttach() {
    attachState = null;
    renderBets();
  }

  function updateAttach(changes) {
    attachState = { ...attachState, ...changes };
    renderBets();
  }

  function makeEl(tag, className, text) {
    const node = document.createElement(tag);
    if (className) node.className = className;
    if (text != null) node.textContent = text;
    return node;
  }

  function makeButton(className, text, onClick) {
    const button = makeEl("button", className, text);
    button.type = "button";
    button.addEventListener("click", (event) => {
      event.stopPropagation();
      onClick();
    });
    return button;
  }

  // POST (attach) or DELETE (undo) /pins.json; the reply carries every pin and the whole crosswalk.
  async function sendPins(method, body, betId) {
    const query = betId ? `?betId=${encodeURIComponent(betId)}` : "";
    const reply = await serviceRequest(method, `/pins.json${query}`, body);
    if (!reply || !Array.isArray(reply.pins) || !Array.isArray(reply.crosswalk)) throw new Error("pins.json reply has no pins or crosswalk array");
    return reply;
  }

  async function submitAttach(bet, game, plan) {
    updateAttach({ busy: true, error: null });
    try {
      const reply = await sendPins("POST", attachLib.pinRequest(bet, game, plan));
      console.info(`[unabated-ticket] attach: ${bet.id} -> board event ${game.eventId}, ${plan.crosswalk.length} name(s) taught`);
      attachState = null;
      await applyServiceTables({ crosswalk: reply.crosswalk, pins: reply.pins });
    } catch (error) {
      console.warn("[unabated-ticket] attach failed:", error.message);
      updateAttach({ busy: false, error: `Attach failed: ${error.message}. Is the bets service running and up to date?` });
    }
  }

  async function undoAttach(bet) {
    try {
      const reply = await sendPins("DELETE", null, bet.id);
      console.info(`[unabated-ticket] attach undone: ${bet.id}, ${reply.removedRows} name(s) removed`);
      view.betsSettingsError.textContent = "";
      await applyServiceTables({ crosswalk: reply.crosswalk, pins: reply.pins });
    } catch (error) {
      console.warn("[unabated-ticket] undo failed:", error.message);
      view.betsSettingsError.textContent = `Undo failed: ${error.message}`;
    }
  }

  // Step 1's scope line and game list for the current query, into their two
  // holders. Typing redraws only these, never the whole tab.
  function fillCandidates(bet, scopeHolder, listHolder) {
    const { scope, events, more } = attachLib.attachCandidates(bet, boardLines(), { query: attachState.query });
    scopeHolder.textContent = scope;
    const list = makeEl("ol", "candidates");
    for (const { event, why } of events) {
      const item = makeEl("li", why ? "candidate best" : "candidate");
      item.tabIndex = 0;
      item.append(makeEl("div", "c-game", attachLib.gameLabel(event)), makeEl("div", "c-meta", attachLib.gameMeta(event)));
      if (why) item.append(makeEl("span", "tag held c-why", why));
      const choose = () => updateAttach({ step: "confirm", eventId: event.eventId, swapped: false, error: null });
      item.addEventListener("click", choose);
      item.addEventListener("keydown", (keyEvent) => {
        if (keyEvent.key === "Enter") choose();
      });
      list.append(item);
    }
    const notes = [];
    if (!events.length) notes.push(makeEl("div", "muted", "No game on the board fits. Search a team by name."));
    if (more) notes.push(makeEl("div", "muted", `${more} more: type to narrow the list.`));
    listHolder.replaceChildren(list, ...notes);
  }

  // Step 1: the games to pick from, a search box, best fit first.
  function attachPickStep(bet, panel) {
    const head = makeEl("div", "attach-head", "Which game?");
    const scopeHolder = makeEl("span", "scope");
    head.append(scopeHolder);
    const search = makeEl("input", "attach-search");
    search.type = "search";
    search.placeholder = "Search a team";
    search.value = attachState.query;
    const listHolder = makeEl("div", "candidates-holder");
    search.addEventListener("input", () => {
      attachState.query = search.value;
      fillCandidates(bet, scopeHolder, listHolder);
    });
    fillCandidates(bet, scopeHolder, listHolder);
    panel.append(head, search, listHolder,
      makeEl("div", "muted", "Not listed? The game may not be on the board yet. The bet stays flagged until you attach it."));
  }

  // Step 2: the picked game, what the venue's names mean on it, the bet restated.
  function attachConfirmStep(bet, panel) {
    const game = attachLib.boardGames(boardLines(), null).find((candidate) => candidate.eventId === attachState.eventId);
    if (!game) {
      panel.append(makeEl("div", "attach-error", "That game has left the board."),
        makeButton("linkbtn", "Pick another game", () => updateAttach({ step: "pick", eventId: null })));
      return;
    }
    const plan = attachLib.attachPlan(bet, game, { swapped: attachState.swapped, lines: boardLines() });
    const gameHead = makeEl("div", "attach-head", "Game");
    gameHead.append(makeButton("linkbtn", "Change", () => updateAttach({ step: "pick", eventId: null, error: null })));
    const chosen = makeEl("div", "chosen", attachLib.gameLabel(game));
    chosen.append(makeEl("span", "muted", ` · ${attachLib.gameMeta(game)}`));
    panel.append(gameHead, chosen);
    if (plan.names.length) {
      const namesHead = makeEl("div", "attach-head", `${betsLib.venueLabel(bet.venue)} calls them`);
      namesHead.append(makeButton("linkbtn", "Swap", () => updateAttach({ swapped: !attachState.swapped })));
      const rows = makeEl("div", "map-rows");
      for (const name of plan.names) {
        const row = makeEl("div", "map-row");
        const target = makeEl("span", null, `${name.unabatedTeamName} `);
        target.append(makeEl("span", "muted", name.was ? `${name.eventSide}, was ${name.was}` : name.eventSide));
        row.append(makeEl("span", "vn", name.venueTeamName), makeEl("span", "arrow", "→"), target,
          makeEl("span", name.status === attachLib.STATUS_KNOWN ? "tag" : `tag ${name.status}`, name.status));
        rows.append(row);
      }
      panel.append(namesHead, rows);
    } else {
      panel.append(makeEl("div", "muted", "The bet names no team, so nothing is learned. The attach pins the game."));
    }
    const yourBet = makeEl("div", "your-bet");
    yourBet.append(makeEl("span", "k", "Your bet"), document.createTextNode(plan.betOnGame));
    const learned = plan.crosswalk.map((row) => `"${row.venueTeamName}"`).join(" and ");
    const submit = makeButton("btn primary sm", plan.crosswalk.length ? "Attach and learn" : "Attach", () => submitAttach(bet, game, plan));
    submit.disabled = attachState.busy;
    const actions = makeEl("div", "attach-actions");
    actions.append(submit, makeEl("span", "muted", plan.crosswalk.length
      ? `Next time ${betsLib.venueLabel(bet.venue)} writes ${learned}, it matches on its own.`
      : "Nothing to learn: every name already matches. This pins the bet to this game."));
    panel.append(yourBet, actions);
  }

  function attachPanel(bet) {
    const panel = makeEl("div", "attach-panel");
    if (attachState.step === "confirm") attachConfirmStep(bet, panel);
    else attachPickStep(bet, panel);
    if (attachState.error) panel.append(makeEl("div", "attach-error", attachState.error));
    return panel;
  }

  // One bet row. options: {reason, unmatched, attachable, quiet, flagged, dismissed}. `reason`
  // (an unmatched list row) adds the reason chip and, when attachable, the
  // Attach button and panel; `unmatched` colours the left edge, so an
  // unmatched bet is visible in the open list too; `quiet` greys a row that
  // does not flag; `flagged` (a red row) adds the Dismiss chip; `dismissed`
  // (a row Cal dismissed, in the folded list) adds the tag and Restore. A
  // pinned, matched open bet says "attached" with an Undo.
  function betItem(bet, options) {
    const { reason = null, unmatched = false, attachable = false, quiet = false, flagged = false, dismissed = false } = options || {};
    const li = makeEl("li", unmatched ? "unmatched" : "matched");
    const main = makeEl("div");
    const what = makeEl("div", "bet-what", betsLib.describeBet(bet));
    const pinned = Boolean(bet.pin) && !unmatched;
    if (pinned) what.append(makeEl("span", "tag pinned", "attached"));
    if (dismissed && reason) what.append(makeEl("span", "tag dismissed", "dismissed"));
    const meta = makeEl("div", "bet-meta");
    const venue = bet.venue ? betsLib.venueLabel(bet.venue) : "unknown venue";
    const game = pinned && bet.pin.awayTeamName && bet.pin.homeTeamName ? `${bet.pin.awayTeamName} @ ${bet.pin.homeTeamName}`
      : bet.awayTeam && bet.homeTeam ? `${bet.awayTeam} @ ${bet.homeTeam}` : bet.awayTeam || bet.homeTeam || null;
    const league = reason && bet.league ? String(bet.league).toUpperCase() : null;
    meta.textContent = [venue, league, game].filter(Boolean).join(" · ");
    if (pinned) {
      meta.append(document.createTextNode(" · "));
      meta.append(makeButton("linkbtn", "Undo", () => undoAttach(bet)));
    }
    if (dismissed && reason) {
      meta.append(document.createTextNode(" · "));
      const restore = makeButton("linkbtn", "Restore", () => setDismissed(bet.id, false));
      restore.title = "Flag this bet again.";
      meta.append(restore);
    }
    main.append(what, meta);

    const rail = makeEl("div");
    rail.append(makeEl("span", "bet-stake", bet.stake == null ? "—" : fmtDollars(bet.stake)));
    if (bet.placedAt) rail.append(makeEl("small", "bet-when", betsLib.formatPlacedAt(bet.placedAt)));
    li.append(main, rail);
    if (!reason) return li;

    const reasonRow = makeEl("div", quiet ? "reason-row quiet" : "reason-row");
    reasonRow.append(makeEl("span", "bet-reason", reason));
    if (flagged) {
      const dismiss = makeButton("dismiss-chip", "Dismiss", () => setDismissed(bet.id, true));
      dismiss.title = "Stop flagging this bet. It moves to Not on the board and still does not size your next bet.";
      reasonRow.append(dismiss);
    }
    const panelOpen = attachState != null && attachState.betId === bet.id;
    if (attachable && panelOpen) {
      // Not while the attach is being saved: its reply, or its error, lands in this panel.
      const cancel = makeButton("attach-btn open", "Cancel", closeAttach);
      cancel.disabled = attachState.busy;
      reasonRow.append(cancel);
    } else if (attachable) {
      reasonRow.append(makeButton("attach-btn", "Attach", () => openAttach(bet.id)));
    }
    li.append(reasonRow);
    if (attachable && panelOpen) li.append(attachPanel(bet));
    return li;
  }

  function renderBets() {
    const now = Date.now();
    const scrollTop = view.tabBets.scrollTop;
    // The search box is rebuilt on every render; keep typing where it was.
    const active = document.activeElement;
    const searchCaret = active && active.classList && active.classList.contains("attach-search") ? active.selectionStart : null;
    renderBetsSources(now);
    const open = state.betRecords.filter((bet) => bet.status === "open")
      .sort((a, b) => Date.parse(b.placedAt || 0) - Date.parse(a.placedAt || 0));
    const unmatched = unmatchedNow();
    const unmatchedIds = new Set(unmatched.map(({ bet }) => bet.id));
    const needsGame = unmatched.filter((entry) => entry.needsGame);
    const needsFix = unmatched.filter((entry) => entry.needsFix);
    const offBoard = unmatched.filter((entry) => !entry.needsGame && !entry.needsFix);
    // A bet that matched, settled or stopped being attachable closes its panel (never mid-POST).
    if (attachState && !attachState.busy && !unmatched.some((entry) => entry.attachable && entry.bet.id === attachState.betId)) attachState = null;

    // What the tab opens with: the money, before the plumbing. A bet whose
    // venue reported no stake is counted separately rather than as zero.
    const priced = open.filter((bet) => typeof bet.stake === "number");
    const atRisk = priced.reduce((total, bet) => total + bet.stake, 0);
    const venues = betsView.sourceRows(betsPayload(), now).filter((row) => row.configured).length;
    view.betsRisk.textContent = fmtDollars(atRisk);
    view.betsRiskCaption.textContent = [
      `at risk · ${open.length} open bet${open.length === 1 ? "" : "s"}`,
      `${venues} venue${venues === 1 ? "" : "s"}`,
      priced.length === open.length ? null : `${open.length - priced.length} with no stake reported`,
    ].filter(Boolean).join(" · ");

    view.betsNeedsBanner.hidden = needsGame.length === 0;
    view.betsNeedsBanner.textContent = betsView.needsGameBanner(needsGame.length);
    view.betsNeedsBlock.hidden = needsGame.length === 0;
    view.betsNeedsCount.textContent = needsGame.length ? String(needsGame.length) : "";
    view.betsNeeds.replaceChildren(...needsGame.map(({ bet, reason, attachable }) => betItem(bet, { reason, unmatched: true, attachable, flagged: true })));
    view.betsFixBanner.hidden = needsFix.length === 0;
    view.betsFixBanner.textContent = betsView.needsFixBanner(needsFix.length);
    view.betsFixBlock.hidden = needsFix.length === 0;
    view.betsFixCount.textContent = needsFix.length ? String(needsFix.length) : "";
    view.betsFix.replaceChildren(...needsFix.map(({ bet, reason }) => betItem(bet, { reason, unmatched: true, flagged: true })));

    view.betsOpenCount.textContent = open.length ? String(open.length) : "";
    view.betsOpen.replaceChildren(...open.map((bet) => betItem(bet, { unmatched: unmatchedIds.has(bet.id) })));
    view.betsOpenEmpty.hidden = open.length > 0;
    view.betsOpenEmpty.textContent = state.betsService && state.betsService.okAt != null ? "No open bets." : "No bets loaded yet.";
    // Say so when nothing is wrong, so "all matched" never looks like "not checked".
    // A bet its source could not read is still a game bet; a prop or a future is not.
    const isGameBet = (bet) => !bet.unmatchable || betsLib.isParseFailure(bet);
    const gameBets = open.filter(isGameBet).length;
    const unmatchedGameBets = unmatched.filter((entry) => isGameBet(entry.bet)).length;
    view.betsMatchNote.hidden = gameBets === 0;
    view.betsMatchNote.textContent = boardLines().length === 0 ? "Waiting for the board to load before checking which game each bet is on."
      : unmatchedGameBets === 0 ? "Every open game bet matches a game on the board."
        : `${gameBets - unmatchedGameBets} of ${gameBets} open game bets match a game on the board.`;

    view.betsOffboard.hidden = offBoard.length === 0;
    view.betsOffboardCount.textContent = offBoard.length ? String(offBoard.length) : "";
    view.betsOffboardList.replaceChildren(...offBoard.map(({ bet, reason, attachable, dismissed }) => betItem(bet, { reason, unmatched: true, attachable, quiet: true, dismissed })));
    // An Attach opened from the folded list keeps the fold open.
    if (attachState && offBoard.some((entry) => entry.bet.id === attachState.betId)) view.betsOffboard.open = true;

    renderCrosswalk();
    view.tabBets.scrollTop = scrollTop;
    const search = view.tabBets.querySelector(".attach-search");
    if (search && (searchCaret != null || (attachState && attachState.focusSearch))) {
      search.focus();
      const caret = searchCaret ?? search.value.length;
      search.setSelectionRange(caret, caret);
      attachState.focusSearch = false;
    }
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

  // ---- Teasers tab (Buckeye 6-point teasers, teaser.js) ----------------------
  //
  // Every render recomputes the legs and the open BFA teasers; the ticket list
  // is rebuilt only when teaser.planTeasers says what it is built on changed,
  // so it holds still while the EVs move. It is computed on every scanner
  // update, bets poll and 5 s tick even with the tab hidden, so the tab's
  // count stays current. Placed tickets come from the bets service's BFA
  // open bets (BFA is Buckeye): nothing to mark here.

  const HOUR_MS = 3600 * 1000;
  // Ticket cards shown before the fold; the rest sit behind "N more tickets".
  const TEASER_TICKETS_SHOWN = 3;
  // Ticket number -> the text its Copy button puts on the clipboard.
  const teaserCopyTexts = new Map();
  // The last failure logged, so a persistent one is logged once, not every 5 s.
  let teaserLastError = null;

  // The board pass, redone when an NFL or CFB snapshot has landed since (a
  // scanner restart clears the load times, so it counts too).
  function currentTeaserBoard() {
    const loadedAt = teaserLib.TEASER_LEAGUE_IDS.map((id) => (scannerStatus && scannerStatus.leagueLoadedAt[id]) || 0).join(",");
    if (!teaserBoard || loadedAt !== teaserBoardLoadedAt) {
      teaserBoard = teaserLib.teaserBoardOf(scannerState);
      teaserBoardLoadedAt = loadedAt;
    }
    return teaserBoard;
  }

  // Both football boards in, or failed. Before that a list would be NFL's
  // alone while CFB — loaded last, the biggest file — is still coming, and a
  // scanner restart's empty board would throw the list's reference fairs away.
  function teaserBoardReady() {
    if (!scannerStatus) return false;
    return teaserLib.TEASER_LEAGUE_IDS.every((id) => scannerStatus.leaguesLoaded.includes(id)
      || Boolean(scannerStatus.leagueErrors && scannerStatus.leagueErrors[id]));
  }

  // The reference fairs of the last build in chrome.storage.local
  // ("teaserRefs": [cut key, P(above)] pairs), so a reopened panel builds on
  // the same fairs and the list holds still across a reopen too.
  function persistTeaserRefs(refs) {
    chrome.storage.local.set({ teaserRefs: Array.from(refs) });
  }

  // The can't-tease marks whose game is still to start; the rest leave storage too.
  function liveTeaserBlocks(now) {
    const live = Object.fromEntries(Object.entries(teaserBlocked).filter(([, startMs]) => startMs > now));
    if (Object.keys(live).length !== Object.keys(teaserBlocked).length) {
      teaserBlocked = live;
      chrome.storage.local.set({ teaserBlocked });
    }
    return new Set(Object.keys(live));
  }

  // Can't tease marks the leg's market (its spread or total, both sides);
  // Restore removes the mark. Either way the list rebuilds on the next render.
  function setTeaserBlocked(leg, blocked) {
    const next = { ...teaserBlocked };
    if (blocked) next[teaserLib.marketKeyOf(leg)] = leg.eventStartMs;
    else delete next[teaserLib.marketKeyOf(leg)];
    teaserBlocked = next;
    chrome.storage.local.set({ teaserBlocked });
    renderTeasers();
  }

  function teaserModel() {
    if (!scannerState || !teaserBoardReady()) return null;
    const now = Date.now();
    const board = currentTeaserBoard();
    const blocked = liveTeaserBlocks(now);
    const legs = teaserLib.teaserLegs(scannerState, { now, maxLineAgeMs: state.edgeSettings.maxLineAgeHours * HOUR_MS, board });
    const placed = teaserLib.openTeasers(state.betRecords, boardLines(), { now, ladderOf: board.ladderOf });
    const straights = teaserLib.heldStraights(state.betRecords, boardLines(), { now, ladderOf: board.ladderOf });
    const kellyBankroll = state.settings.bankroll * state.settings.multiplier;
    const planned = teaserLib.planTeasers({ legs, placed, straights, kellyBankroll, previous: teaserBuild, blocked });
    teaserBuild = planned.build;
    if (planned.rebuilt) persistTeaserRefs(teaserBuild.refs);
    return {
      legs, placed, build: teaserBuild,
      plan: teaserLib.describePlan(teaserBuild, legs, placed),
      legRows: teaserLib.describeLegs(legs, teaserBuild, placed, straights, blocked),
    };
  }

  function fmtWholeDollars(value) {
    return value.toLocaleString("en-US", { style: "currency", currency: "USD", maximumFractionDigits: 0 });
  }

  function fmtSignedDollars(value) {
    return `${value >= 0 ? "+" : "-"}${fmtWholeDollars(Math.abs(value))}`;
  }

  function fmtWin(win) {
    return `${(win * 100).toFixed(1)}%`;
  }

  function fmtEv(ev) {
    return `EV ${ev >= 0 ? "+" : ""}${(ev * 100).toFixed(1)}%`;
  }

  function fmtWholePct(fraction) {
    return `${Math.round(fraction * 100)}%`;
  }

  function kellyWords(multiplier) {
    const names = { 1: "Full", 0.5: "Half", 0.25: "Quarter" };
    return `${names[multiplier] || `${multiplier}x`} Kelly`;
  }

  function plural(count, word) {
    return `${count} ${word}${count === 1 ? "" : "s"}`;
  }

  function teaserTicketSize() {
    return `${teaserLib.LEGS_PER_TICKET}-team · pays +${teaserLib.TICKET_NET_ODDS * 100}`;
  }

  // The bets service's BFA row (betsview.sourceRows), or null before the service has been reached.
  function bfaSourceRow(now) {
    if (!state.betsService || state.betsService.okAt == null) return null;
    return betsView.sourceRows(betsPayload(), now).find((row) => row.venue === "bfa") || null;
  }

  function renderTeasers() {
    let model = null;
    let failure = null;
    try {
      model = teaserModel();
      teaserLastError = null;
    } catch (error) {
      failure = error;
      if (teaserLastError !== error.message) console.error("[unabated-ticket] teasers failed", error);
      teaserLastError = error.message;
    }
    const count = model ? model.plan.tickets.length : 0;
    view.teasersCount.hidden = count === 0;
    view.teasersCount.textContent = String(count);
    if (state.activeTab !== "teasers") return;
    const scrollTop = view.tabTeasers.scrollTop;
    renderTeasersBanners(failure);
    renderTeasersStatus(model);
    renderTeasersSummary(model);
    renderTeasersOpen(model);
    renderTeasersTickets(model);
    renderTeasersLegs(model);
    view.tabTeasers.scrollTop = scrollTop;
  }

  // Red: the tab failed, or the bets service is down (the list stays up, but
  // an open teaser may be missing from it). Amber: no BFA account is read, or
  // BFA's last pull failed.
  function renderTeasersBanners(failure) {
    const now = Date.now();
    const red = [];
    if (failure) red.push(`Teasers failed: ${failure.message}`);
    const service = betsView.serviceStatus(state.betsService, now);
    if (service.unreachable) {
      red.push(`${service.text.charAt(0).toUpperCase()}${service.text.slice(1)}. Open BFA teasers may be missing, so a ticket already placed may be suggested again.`);
    }
    view.teasersService.hidden = red.length === 0;
    view.teasersService.textContent = red.join(" ");
    const bfa = bfaSourceRow(now);
    const amber = !bfa || service.unreachable ? null
      : !bfa.configured ? "The bets service reads no BFA account, so tickets already placed at Buckeye are not known here."
        : bfa.error ? `BFA's last pull failed (${bfa.error}); open teasers are as of ${bfa.ageText} ago.` : null;
    view.teasersWarning.hidden = !amber;
    view.teasersWarning.textContent = amber || "";
  }

  // "NFL · CFB loading", then "NFL · CFB · 15 games on Buckeye's board · 56 of 56 legs priced".
  function renderTeasersStatus(model) {
    if (!scannerStatus) {
      view.teasersStatus.textContent = "Starting the scanner…";
      return;
    }
    const leagues = teaserLib.TEASER_LEAGUE_IDS.map((id) => {
      const label = feed.LEAGUES[id].label;
      if (scannerStatus.leaguesLoaded.includes(id)) return label;
      return scannerStatus.leagueErrors && scannerStatus.leagueErrors[id] ? `${label} unavailable` : `${label} loading`;
    });
    if (!model) {
      view.teasersStatus.textContent = leagues.join(" · ");
      return;
    }
    const games = new Set(model.legs.map((leg) => leg.eventId)).size;
    const priced = model.legs.filter((leg) => leg.win != null).length;
    view.teasersStatus.textContent = `${leagues.join(" · ")} · ${plural(games, "game")} on Buckeye's board · ${priced} of ${plural(model.legs.length, "leg")} priced`;
  }

  function teasersEmptyText(model) {
    if (!model) return "Waiting for Buckeye's NFL and CFB board…";
    const { plan, build } = model;
    if (plan.tickets.length) return null;
    if (build.reason === teaserLib.REASON_FEW_LEGS) return `Fewer than ${teaserLib.LEGS_PER_TICKET} games with a priced Buckeye leg right now: nothing to tease.`;
    if (build.reason) return `No tickets: ${build.reason}.`;
    if (build.partialHeldBack > 0) {
      return `Nothing more to bet: the next ticket would be a partial ${fmtWholeDollars(build.partialHeldBack)}, and a 4-team ticket under $200 is already open at BFA (one partial ticket per set).`;
    }
    return build.placed.length
      ? "Nothing more to bet: no other ticket raises the Kelly growth with the open teasers held."
      : "No ticket worth betting right now: no 4-team ticket raises the Kelly growth at these fairs.";
  }

  function summaryCell(label, value, small) {
    const cell = makeEl("div", "payout-cell");
    const valueEl = makeEl("div", "payout-value", value);
    if (small) valueEl.append(" ", makeEl("small", null, small));
    cell.append(makeEl("div", "payout-label", label), valueEl);
    return cell;
  }

  // Tickets, dollars, expected profit, the chance to make money and the
  // chance every ticket loses — over the open teasers too when there are any.
  function renderTeasersSummary(model) {
    const empty = teasersEmptyText(model);
    view.teasersEmpty.hidden = !empty;
    view.teasersEmpty.textContent = empty || "";
    const summary = model ? model.plan.summary : null;
    const shown = Boolean(summary) && (summary.count > 0 || summary.placedCount > 0);
    view.teasersSummary.hidden = !shown;
    if (!shown) return;
    view.teasersSummaryLabel.textContent = summary.count === 0 ? "Nothing more to bet"
      : summary.placedCount ? `Bet ${summary.count} more` : `Bet ${plural(summary.count, "ticket")}`;
    view.teasersSummaryStake.textContent = fmtWholeDollars(summary.stake);
    const cells = [];
    if (summary.placedCount) {
      // The open tickets still in the math; the Open at BFA header below counts every open one.
      const allOpen = model.placed.length === summary.placedCount;
      cells.push(summaryCell(allOpen ? "Placed" : "Placed, in play", fmtWholeDollars(summary.placedStake)));
      if (summary.expectedAll != null) cells.push(summaryCell(`Expected, all ${summary.count + summary.placedCount}`, fmtSignedDollars(summary.expectedAll)));
    } else {
      cells.push(summaryCell("Expected", fmtSignedDollars(summary.expected), summary.stake > 0 ? fmtWholePct(summary.expected / summary.stake) : null));
    }
    if (summary.makesMoney != null) cells.push(summaryCell("Makes money", fmtWholePct(summary.makesMoney)));
    if (summary.allLose != null) cells.push(summaryCell("All lose", fmtWholePct(summary.allLose)));
    view.teasersSummaryCells.replaceChildren(...cells);
    const straights = summary.straights;
    view.teasersSummaryNote.textContent = `${kellyWords(state.settings.multiplier)} on ${fmtWholeDollars(state.settings.bankroll)} · `
      + (summary.placedCount ? "open BFA teasers held fixed" : "tickets that share a leg are sized together")
      + (straights.count ? ` · counts ${fmtWholeDollars(straights.stake)} of straight bets on ${plural(straights.games, "game")}` : "");
  }

  // What Copy puts on the clipboard: BFA's own "[rotation] side" for each leg.
  function teaserCopyText(ticket) {
    const legs = ticket.legs.map((leg) => `${leg.rotation != null ? `[${leg.rotation}] ` : ""}${leg.label}`).join(" / ");
    return `Buckeye ${teaserLib.LEGS_PER_TICKET}-team ${teaserLib.TEASER_POINTS}-pt teaser ${fmtWholeDollars(ticket.stake)}: ${legs} | ${fmtEv(ticket.ev)}`;
  }

  function teaserTicketCard(ticket) {
    const card = makeEl("li", "ticket-card");
    const head = makeEl("div", "tk-head");
    const rail = makeEl("span", "tk-rail");
    rail.append(makeEl("span", "tk-stake", fmtWholeDollars(ticket.stake)), makeEl("span", "tk-ev", fmtEv(ticket.ev)));
    head.append(makeEl("span", "tk-num", `#${ticket.number}`), makeEl("span", "tk-size", teaserTicketSize()), rail);
    const legs = makeEl("ol", "tk-legs");
    for (const leg of ticket.legs) {
      const row = makeEl("li");
      row.append(makeEl("span", "tl-leg", leg.label), makeEl("span", "tl-from", leg.fromLabel), makeEl("span", "tl-win", fmtWin(leg.win)));
      legs.append(row);
    }
    const actions = makeEl("div", "tk-actions");
    const copy = makeEl("button", "btn sm", "Copy");
    copy.type = "button";
    copy.dataset.copyTicket = String(ticket.number);
    actions.append(makeEl("span", "tk-copied"), copy);
    card.append(head, legs, actions);
    teaserCopyTexts.set(ticket.number, teaserCopyText(ticket));
    return card;
  }

  function renderTeasersTickets(model) {
    teaserCopyTexts.clear();
    const tickets = model ? model.plan.tickets : [];
    view.teasersListLabel.hidden = tickets.length === 0;
    view.teasersListCount.textContent = tickets.length ? String(tickets.length) : "";
    view.teasersList.replaceChildren(...tickets.slice(0, TEASER_TICKETS_SHOWN).map(teaserTicketCard));
    const rest = tickets.slice(TEASER_TICKETS_SHOWN);
    view.teasersMore.hidden = rest.length === 0;
    view.teasersMoreLabel.textContent = `${plural(rest.length, "more ticket")} · ${fmtWholeDollars(rest.reduce((sum, ticket) => sum + ticket.stake, 0))}`;
    view.teasersMoreList.replaceChildren(...rest.map(teaserTicketCard));
  }

  // One open BFA teaser, read-only: its legs at BFA's numbers, each at its
  // current fair or why it counts as won.
  function openTeaserCard(ticket) {
    const card = makeEl("li", `ticket-card placed${ticket.inPlay ? "" : " done"}`);
    const head = makeEl("div", "tk-head");
    const rail = makeEl("span", "tk-rail");
    const hasDollars = ticket.reason == null;
    if (ticket.inPlay) rail.append(makeEl("span", "tk-ev", fmtEv(ticket.winAll * (1 + ticket.toWin / ticket.stake) - 1)));
    const size = hasDollars ? `${ticket.legCount}-team · ${fmtWholeDollars(ticket.stake)} to win ${fmtWholeDollars(ticket.toWin)}` : `${ticket.legCount}-team`;
    head.append(makeEl("span", "tk-size", size), rail);
    const legs = makeEl("ol", "tk-legs");
    for (const leg of ticket.legs) {
      const row = makeEl("li", leg.state === teaserLib.LEG_LIVE ? null : "counted");
      row.append(makeEl("span", "tl-leg", leg.label), makeEl("span", "tl-from", leg.note || ""), makeEl("span", "tl-win", leg.win == null ? "—" : fmtWin(leg.win)));
      legs.append(row);
    }
    const actions = makeEl("div", "tk-actions");
    if (!ticket.inPlay) actions.append(makeEl("span", "muted", outOfMathText(ticket)));
    actions.append(makeEl("span", "tag held", ticket.placedAt ? `placed ${betsLib.formatPlacedAt(ticket.placedAt)}` : "placed"));
    card.append(head, legs, actions);
    return card;
  }

  // Why an open ticket is not in the math, by its legs: every game started,
  // or no leg joined a game Buckeye's board prices (a basketball teaser, CFB
  // not loaded, a start the board disagrees with by over 12 h).
  function outOfMathText(ticket) {
    if (ticket.reason) return `${ticket.reason}: out of the math`;
    if (ticket.legs.every((leg) => leg.state === teaserLib.LEG_STARTED)) return "every game has started: out of the math";
    return "no leg priced on the board: out of the math";
  }

  // Always shown once the board is up, so "none open · BFA pulled 41 s ago"
  // says how fresh the answer is right after a ticket is placed.
  function renderTeasersOpen(model) {
    view.teasersOpen.hidden = !model;
    if (!model) return;
    const placed = model.placed;
    const bfa = bfaSourceRow(Date.now());
    const pull = !bfa ? "bets service not reached yet"
      : !bfa.configured ? "no BFA account read"
        : bfa.fetchedAt ? `BFA pulled ${bfa.ageText} ago` : "no BFA pull yet";
    const dollars = (tickets) => tickets.reduce((sum, ticket) => sum + (typeof ticket.stake === "number" ? ticket.stake : 0), 0);
    const inPlay = placed.filter((ticket) => ticket.inPlay);
    const amount = inPlay.length === placed.length ? fmtWholeDollars(dollars(placed))
      : `${fmtWholeDollars(dollars(placed))}, ${fmtWholeDollars(dollars(inPlay))} in play`;
    view.teasersOpenCount.textContent = placed.length ? String(placed.length) : "";
    view.teasersOpenNote.textContent = [placed.length ? `${amount} · each until its last game starts` : "none open", pull].join(" · ");
    view.teasersOpenList.replaceChildren(...placed.map(openTeaserCard));
  }

  function legStandingText(row) {
    const open = row.openTickets ? ` · ${row.openTickets} open` : "";
    if (row.standing !== teaserLib.STANDING_POOL) return `${row.note}${open}`;
    return `${row.inTickets ? `in ${plural(row.inTickets, "ticket")}` : "in no ticket"}${open}`;
  }

  // The straight bets on a pool leg's game, the Edges tab's tags: `held $X`
  // on the leg's side, `against $Y` on the other, or a bare `game` when the
  // game has straights and none of them counts. Every bet in the tooltip.
  function teaserStraightTags(straights) {
    if (!straights) return [];
    const tooltip = straights.counted.map((bet) => bet.label)
      .concat(straights.leftOut.map((bet) => `${bet.label} — not counted: ${bet.reason}`)).join("\n");
    const tag = (kind, text) => {
      const el = makeEl("span", `tag ${kind}`, text);
      el.title = tooltip;
      return el;
    };
    const tags = [];
    if (straights.held > 0) tags.push(tag("held", `held ${betsLib.formatStake(straights.held)}`));
    if (straights.against > 0) tags.push(tag("against", `against ${betsLib.formatStake(straights.against)}`));
    if (tags.length === 0 && straights.leftOut.length > 0) tags.push(tag("game", "game"));
    return tags;
  }

  // A college leg on show (in the pool or above break-even) offers Can't
  // tease; a marked one shows the tag and Restore. NFL legs never: Buckeye
  // teases every NFL game.
  function teaserBlockControl(row) {
    if (row.standing === teaserLib.STANDING_BLOCKED) return makeButton("linkbtn", "Restore", () => setTeaserBlocked(row.leg, false));
    const onShow = row.standing === teaserLib.STANDING_POOL || row.standing === teaserLib.STANDING_OUT;
    if (!onShow || !teaserLib.canBlock(row.leg)) return null;
    const chip = makeButton("dismiss-chip", "Can't tease", () => setTeaserBlocked(row.leg, true));
    chip.title = `Buckeye won't tease this. Takes the game's ${teaserLib.marketNameOf(row.leg)} (both sides) out of the tickets until the game starts.`;
    return chip;
  }

  function teaserLegRow(row) {
    const { leg } = row;
    const pool = row.standing === teaserLib.STANDING_POOL;
    const item = makeEl("li", `edge-row leg-row ${pool ? "tier-hot" : "tier-thin out"}`);
    const main = makeEl("div");
    const meta = makeEl("div", "edge-meta");
    meta.append(`${leg.matchup} · ${leg.leagueLabel} · ${fmtStart(new Date(leg.eventStartMs).toISOString())} · `, untilEl(leg.eventStartMs));
    const book = makeEl("div", "edge-book");
    book.append(makeEl("span", "price", `Buckeye ${leg.bookLabel}`),
      makeEl("span", "age", ` · teased ${teaserLib.TEASER_POINTS} · ${fmtLineAge(leg.modifiedMs)}`));
    const side = makeEl("div", "edge-side", leg.label);
    if (row.standing === teaserLib.STANDING_BLOCKED) side.append(makeEl("span", "tag dismissed", "can't tease"));
    main.append(side, meta, book);
    const rail = makeEl("div", "edge-rail");
    rail.append(makeEl("span", `edge-pct tier-${pool ? "hot" : "thin"}`, leg.win == null ? "—" : fmtWin(leg.win)),
      makeEl("span", "leg-in", legStandingText(row)), ...teaserStraightTags(row.straights));
    const blockControl = teaserBlockControl(row);
    if (blockControl) rail.append(blockControl);
    item.append(main, rail);
    return item;
  }

  function breakEvenDivider() {
    const divider = makeEl("li", "be-divider");
    divider.setAttribute("role", "separator");
    divider.append(makeEl("span", null, `${fmtWin(teaserLib.BREAK_EVEN_WIN)} a leg breaks even at +${teaserLib.TICKET_NET_ODDS * 100}`));
    return divider;
  }

  // The pool's legs and the priced legs above break-even, best first, with
  // the break-even line where it falls; the rest of the games behind a fold,
  // the markets marked can't tease first.
  function renderTeasersLegs(model) {
    const rows = model ? model.legRows : [];
    const gameCount = new Set(rows.map((row) => row.leg.eventId)).size;
    view.teasersLegsLabel.hidden = rows.length === 0;
    view.teasersLegsCount.textContent = gameCount ? String(gameCount) : "";
    const shownStandings = [teaserLib.STANDING_POOL, teaserLib.STANDING_OUT];
    const shown = rows.filter((row) => shownStandings.includes(row.standing));
    const blockedRows = rows.filter((row) => row.standing === teaserLib.STANDING_BLOCKED);
    const otherFolded = rows.filter((row) => !shownStandings.includes(row.standing) && row.standing !== teaserLib.STANDING_BLOCKED);
    const folded = [...blockedRows, ...otherFolded];
    const items = [];
    let dividerPlaced = false;
    for (const row of shown) {
      if (!dividerPlaced && row.leg.win < teaserLib.BREAK_EVEN_WIN) {
        items.push(breakEvenDivider());
        dividerPlaced = true;
      }
      items.push(teaserLegRow(row));
    }
    if (!dividerPlaced && folded.length) items.push(breakEvenDivider());
    view.teasersLegs.replaceChildren(...items);
    view.teasersLegsMore.hidden = folded.length === 0;
    const counts = [
      [teaserLib.STANDING_BELOW, "below break-even"],
      [teaserLib.STANDING_OTHER_MARKET, "on an open teaser's other market"],
      [teaserLib.STANDING_UNPRICED, "with no fair"],
    ].map(([standing, words]) => [otherFolded.filter((row) => row.standing === standing).length, words])
      .filter(([count]) => count > 0).map(([count, words]) => `${count} ${words}`);
    const labelParts = [];
    if (otherFolded.length) labelParts.push(`${plural(otherFolded.length, "more game")}: ${counts.join(", ")}`);
    if (blockedRows.length) labelParts.push(`${blockedRows.length} can't tease`);
    view.teasersLegsMoreLabel.textContent = labelParts.join(" · ");
    view.teasersLegsMoreList.replaceChildren(...folded.map(teaserLegRow));
  }

  view.tabTeasers.addEventListener("click", async (event) => {
    const button = event.target.closest("button[data-copy-ticket]");
    if (!button) return;
    const status = button.parentElement.querySelector(".tk-copied");
    try {
      await navigator.clipboard.writeText(teaserCopyTexts.get(Number(button.dataset.copyTicket)) || "");
      status.textContent = "Copied";
    } catch (error) {
      status.textContent = `Copy failed: ${error.message}`;
    }
  });

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
    renderTeasers();
  }

  function fillSettingInputs() {
    view.bankroll.value = state.settings.bankroll;
    view.multiplier.value = state.settings.multiplier;
  }

  // ---- wiring --------------------------------------------------------------

  async function load() {
    const local = await chrome.storage.local.get(DEFAULT_SETTINGS);
    state.settings = edgeRows.sanitizeStakeSettings(local);
    fillSettingInputs();
    const relay = await chrome.storage.local.get(["ticket", "error", "watchStatus", "pageReady", "pageCheck", "booksFilter", "edges", "alerts", "alertLog", "activeTab", "locateResult", "betsService", "betsSettings", "teamsIndex", "teaserRefs", "teaserBlocked"]);
    // The reference fairs the last Teasers list was built on: a seed with no
    // key, so the first plan rebuilds — on these fairs, where still within a
    // point of the live ones — and the list is the one this panel last showed.
    if (Array.isArray(relay.teaserRefs)) {
      const pairs = relay.teaserRefs.filter((pair) => Array.isArray(pair) && typeof pair[0] === "string" && typeof pair[1] === "number");
      teaserBuild = { key: null, refs: new Map(pairs) };
    }
    if (relay.teaserBlocked && typeof relay.teaserBlocked === "object") {
      teaserBlocked = Object.fromEntries(Object.entries(relay.teaserBlocked).filter(([, startMs]) => typeof startMs === "number"));
    }
    const leaguePaths = Array.from(new Set(Object.values(feed.LEAGUES).map((league) => league.path)));
    const storedLive = await chrome.storage.local.get(leaguePaths.map(live.storageKeyOf));
    for (const [key, payload] of Object.entries(storedLive)) {
      if (payload) liveStore[key.slice(live.STORAGE_PREFIX.length)] = payload;
    }
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
      state.pins = Array.isArray(storedBets.pins) ? storedBets.pins : [];
      const storedRecords = Array.isArray(storedBets.bets) ? storedBets.bets : [];
      state.betRecords = betsLib.pruneForRetention(betsLib.rekeyRecords(betsLib.applyPins(storedRecords, state.pins), state.crosswalk), Date.now());
      state.dismissedBetIds = betsView.keepDismissedOpen(storedBets.dismissed, state.betRecords);
      state.knownStarts = betsView.keepKnownStartsOpen(storedBets.knownStarts, {}, state.betRecords);
      setFillFairs(fillfair.mergeFillFairs([], Array.isArray(storedBets.fillFairs) ? storedBets.fillFairs : [], state.betRecords));
      state.betsService = {
        payload: storedBets.payload || null, okAt: storedBets.okAt ?? null,
        error: storedBets.error ?? null, errorAt: storedBets.errorAt ?? null, unreachableSince: storedBets.unreachableSince ?? null,
      };
    }
    fillEdgeSettingInputs();
    fillAlertSettingInputs();
    fillBetsSettingInputs();
    showTab(relay.activeTab);
    renderBetsHeader();
    render();
    renderEdges();
    renderLocate();
    if (state.activeTab === "bets") renderBets();
    renderTeasers();
    startBetsPolling();
    await scanner.start(scannerLeaguesOf(state.edgeSettings.leagues));
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
    const liveKeys = Object.keys(changes).filter((key) => key.startsWith(live.STORAGE_PREFIX));
    if (liveKeys.length) {
      for (const key of liveKeys) {
        const payload = changes[key].newValue || null;
        const league = key.slice(live.STORAGE_PREFIX.length);
        if (payload) liveStore[league] = payload;
        else delete liveStore[league];
      }
      renderLive();
      processAlerts().catch((error) => console.error("[unabated-ticket] alerts failed", error));
    }
    if ("locateResult" in changes) {
      state.locateResult = changes.locateResult.newValue || null;
      if (state.locateResult && state.locating && state.locateResult.at >= state.locating.at) state.locating = null;
      renderLocate();
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
  view.edgesMinToWin.addEventListener("input", onEdgeSettingsInput);
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
    // Even hidden: a game that starts leaves the list, and the tab's count follows.
    renderTeasers();
  }, 5000);

  // The Live block's clock: staleness and "fair set 0:38 ago" move without a
  // new payload. Idle (no DOM work) unless a live game is on a screen.
  setInterval(() => {
    if (state.activeTab === "edges" && livePayloads().length) renderLive();
  }, 1000);

  load().catch((error) => {
    view.errorDetail.textContent = error.message;
    show("error");
  });
})();
