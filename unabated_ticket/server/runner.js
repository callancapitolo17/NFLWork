#!/usr/bin/env node
// Unabated Ticket server runner — the panel's Edges scan, headless, for the
// phone page (plan: phone_page_plan.md, step 1).
//
//   node unabated_ticket/server/runner.js        (Node 18+, no build, no npm deps)
//
// Loads the extension's own pure modules unchanged (scanner.js, feed.js,
// edgerows.js, tailflex.js, teaser.js, betsview.js, bets.js, teams.js,
// ladder.js, condkelly.js, kelly.js, edgemove.js, fillfair.js) and runs the
// scan loop the panel runs while it is open and visible: a snapshot per
// league on start — the settings' leagues plus NFL and CFB, which the panel
// always loads for its Teasers tab and whose open BFA teasers size the Edges
// rows, plus the leagues of open straight bets, whose closing fairs it saves —
// each league re-downloaded on its own size-based cadence, a full
// resync every 10 min (scanner.js). It never pauses — the panel pauses when
// hidden; this is a panel that is always in view.
//
// Inputs (env):
//   UNABATED_RUNNER_HOST   bind address, default 127.0.0.1 (loopback only;
//                          the Tailscale address comes with the deploy step)
//   UNABATED_RUNNER_PORT   default 8095
//   BETS_SERVICE_URL       the bets service, default http://127.0.0.1:8094
//   KALSHI_BASE_URL        Kalshi's API, default https://api.elections.kalshi.com/trade-api/v2
// Network:
//   GET content.unabated.com league snapshots (public, no login), as the panel does.
//   GET <bets service>/bets.json every 30 s (the panel's cadence) — open bets
//       for conditional Kelly, the team crosswalk, pins and saved fill fairs.
//   GET <bets service>/settings.json every 10 s — bankroll, Kelly multiplier,
//       books, minimum edge, line age, min liq to win, alts, sort, cards
//       (bets.duckdb::edge_settings; null fields are the panel's defaults).
//   POST <bets service>/closing_fairs.json every 30 s, only the rows that
//       changed — each open straight bet's Unabated fair on its own line
//       while its game is still to start (server/closefair.js); the service
//       keeps the newest pre-start one, the close CLV is measured against.
//   POST <bets service>/crosswalk.json and /fill_fairs.json — what the panel
//       learns while it is open, learned here all the time: after every bets
//       poll and scanner update, the team names of open bets the board joined
//       by venue id (bets.learnCrosswalk) and each new open bet's fill-time
//       fair off the line history (fillfair.captureFillFairs). Both tables are
//       insert-only on the service, so the panel learning the same rows is harmless.
//   POST <bets service>/fill_books.json — each new open bet's other books at
//       its number when it was placed (server/fillbooks.js): how many posted
//       it, how many beat the price taken, the best of them. Runner only.
//   Every league's snapshot every 60 s (scanner maxRefreshMs: the taker
//       trades against the fair, so the panel's 2-5 min tiers are too slow).
//   GET <Kalshi>/markets/<ticker>/orderbook every 30 s for each Kalshi Edges
//       row the automated taker looks at, and GET <Kalshi>/markets/<ticker>
//       before every order, so the market is active and its strike and team
//       are the row's (public, no key; server/taker.js marketMismatch).
//   POST <bets service>/taker/order.json — each "bet" decision, with the ids
//       of the held bets its stake was sized on: the bets service's order
//       desk places it on Kalshi on the bot's own key (MONEY; stop file, $500
//       ceiling and the unseen-fill rule there; bets_service/taker_orders.py).
//   GET <Novig>/v3/public/catalog/markets?league=..&eventStatus=OPEN_PREGAME
//       and /v3/public/catalog/markets/<id>/book every 30 s while Novig Edges
//       rows are up and the bets service reads the bot's Novig subaccount;
//       /v3/public/catalog/markets/<id> and /events/<id> before every Novig
//       order (public, no key; taker.novigMarketMismatch).
//   POST <bets service>/taker/novig_order.json — each Novig "bet" decision
//       (MONEY, on the bot's own Novig subaccount; bets_service/novig_taker_desk.py).
//   POST <bets service>/taker_decisions.json — each taker decision that
//       changed since it was last logged (taker.duckdb::taker_decisions).
//   GET <Kalshi>/milestones?related_event_ticker=<event> before an order or a
//       Place quote, cached per event for 10 min: Kalshi's own start for the
//       game. The order and the quote carry the earlier of it and Unabated's
//       eventStart, which the desks refuse on once passed (Unabated's alone
//       when Kalshi lists none; public, no key; taker.kalshiStartMsOf).
//   GET <Kalshi>/markets/<ticker>/orderbook and /markets/<ticker> once per
//       GET /kalshi_quote.json (Cal's Place button, never on a timer).
// Serves (HTTP on UNABATED_RUNNER_HOST:UNABATED_RUNNER_PORT, no auth):
//   GET /edges.json  the Edges list the panel would show for those settings
//                    (server/edges_payload.js documents the shape)
//   GET /teasers.json  the panel's Teasers tab, as words (server/teasers_payload.js)
//   ?fresh=1 on /edges.json or /teasers.json reads /bets.json first, so a
//                    Dismiss, Can't tease or Attach the page just wrote shows at once
//   GET /attach.json?betId=..&query=..  Attach step 1: the board games to pick
//                    from for an open bet, best fit first (attach.js)
//   GET /attach.json?betId=..&eventId=..&swapped=0|1  Attach step 2: the plan
//                    and the POST /pins.json body the page sends the bets service
//   GET /scenarios.json  the Bet Tracker's Live tab: every game in progress
//                    with open bets, its results and the P&L and kickoff
//                    chance of each (server/scenarios.js documents the shape)
//   GET /taker.json  the automated taker's last pass: every Kalshi row it
//                    looked at, what it decided and what each order did
//   GET /kalshi_quote.json?eventId=..&periodTypeId=..&betTypeId=..&sideIndex=..&points=..
//                    the Ticket's Place button (bets_service/manual_orders.py):
//                    the board's Kalshi line with those fields, the contract it
//                    buys, its live ask and Kalshi's market check (the taker's
//                    own resolveKalshiContract / quoteOf / marketMismatch) and
//                    the Kalshi API it read; places nothing. points omitted
//                    (and ignored) on a moneyline.
//   GET /novig_quote.json?eventId=..&periodTypeId=..&betTypeId=..&sideIndex=..&points=..
//                    the same for Novig Place (bets_service/novig_place.py):
//                    the board's Novig line, its outcome, market and the other
//                    side's outcome, its live price, Novig's market check
//                    (the taker's own novigContractOf / novigQuoteOf /
//                    novigMarketMismatch), the Novig API read and the line
//                    as a bet record names it; places nothing. Spreads and
//                    totals only.
//   GET /health      {ok, generatedAt, uptimeSec, scanner, betsService, settings, closingFairs, taker}
//   Every request whose Host header is not a name this runner serves on is
//   refused with 403 (DNS rebinding, #125); any verb but GET is 405.
// Threads: main() runs createRunner on a worker thread (server/board_thread.js)
// and keeps only the HTTP server, which answers the pages from that thread's
// last copy, so no board walk holds up a phone reply.
// Side effects: none on disk — no DuckDB, no files. Its writes go through the
// bets service: the closing fairs (bets.duckdb::bet_closing_fairs), the
// crosswalk rows, fill fairs and fill books above (insert-only tables) and
// the taker's decisions (taker.duckdb::taker_decisions, append). State (feed, line
// history, bets, each game's last pregame odds for
// the Live tab) is in memory and rebuilt on restart: a game already under way
// when the runner starts has its card but no kickoff chances. Logs state
// changes (not every poll) to stdout/stderr.

"use strict";

const http = require("node:http");
const { monitorEventLoopDelay } = require("node:perf_hooks");
const scannerLib = require("../extension/scanner.js");
const teams = require("../extension/teams.js");
const betsLib = require("../extension/bets.js");
const betsView = require("../extension/betsview.js");
const fillfair = require("../extension/fillfair.js");
const tailflex = require("../extension/tailflex.js");
const teaser = require("../extension/teaser.js");
const edgeRows = require("../extension/edgerows.js");
const edgesPayload = require("./edges_payload.js");
const closefair = require("./closefair.js");
const fillbooks = require("./fillbooks.js");
const ladderLib = require("../extension/ladder.js");
const scenarios = require("./scenarios.js");
const teasersPayload = require("./teasers_payload.js");
const attachLib = require("../extension/attach.js");
const pricecheck = require("../extension/pricecheck.js");
const taker = require("./taker.js");
const feedLib = require("../extension/feed.js");
const boardThread = require("./board_thread.js");

const DEFAULT_HOST = "127.0.0.1";
const DEFAULT_PORT = 8095;
// A moneyline has no number: the Place button's line matches on side alone.
const MONEYLINE_BET_TYPE_ID = 1;
const DEFAULT_BETS_SERVICE_URL = "http://127.0.0.1:8094";
const DEFAULT_KALSHI_BASE_URL = "https://api.elections.kalshi.com/trade-api/v2";
const DEFAULT_NOVIG_BASE_URL = "https://api.novig.com";
// Novig's CloudFront refused an unnamed client once (2026-10-10): every Novig read names itself.
const NOVIG_USER_AGENT = "unabated-ticket-runner/1.0";
// Novig's public catalog: one page holds up to 5000 markets; a slate is a few hundred.
const NOVIG_MARKETS_PAGE_LIMIT = 5000;
const NOVIG_MARKETS_MAX_PAGES = 10;
// A Place quote reuses Novig's catalog when it is at most this old: an outcome
// never moves to another market, and the quote still reads the market, its
// event and its book live (2026-10-17: a fresh catalog plus three reads one
// after another ran past the bets service's deadline, so the first click failed).
const NOVIG_CATALOG_REUSE_MS = 60 * 1000;
// feed.js league path -> Novig's canonical league name (GET /v3/types/leagues).
const NOVIG_LEAGUE_OF = { nfl: "NFL", cfb: "NCAAF", nba: "NBA", cbb: "NCAAB", wnba: "WNBA", mlb: "MLB", nhl: "NHL" };
// The bets service's sources the Novig stake is sized on: Cal's app account and the bot's subaccount.
const NOVIG_HOLDINGS_SOURCES = ["novig", "novig_bot"];
// The taker decides on every Kalshi row this often; Kalshi's public reads allow far more.
const TAKER_PASS_MS = 30 * 1000;
// Two reads in flight stay well under Kalshi's 20 reads/s basic limit at ~100 ms a read.
const KALSHI_BOOK_READS_AT_ONCE = 2;
// The taker trades against Unabated's fair, which lags the market by a minute
// or two already: every league re-downloads at least this often (the panel's
// size tiers go up to 5 min), and a row whose league snapshot landed longer
// ago than TAKER_MAX_SNAPSHOT_AGE_MS is stale_feed, not a bet.
const TAKER_SNAPSHOT_REFRESH_MS = 60 * 1000;
const TAKER_MAX_SNAPSHOT_AGE_MS = 3 * 60 * 1000;
// The stake is conditional on the held bets, read off the bets service's
// Kalshi poll (60 s): a poll that failed, or none this recent, and the
// holdings are unknown, so the pass decides nothing.
const TAKER_MAX_HELD_AGE_MS = 5 * 60 * 1000;
// Kalshi's start for a game (GET /milestones) is read once per event this
// often; a read that failed is retried sooner. Starts move by days, not minutes.
const KALSHI_START_TTL_MS = 10 * 60 * 1000;
const KALSHI_START_RETRY_MS = 60 * 1000;
// The desk sends, cancels and reads fills, each up to 10 s at Kalshi.
const TAKER_ORDER_TIMEOUT_MS = 45 * 1000;
// Decisions held for the bets service while it is down; past this the oldest go (logged).
const TAKER_LOG_MAX_UNSENT = 5000;
const TAKER_LOG_ROWS_PER_POST = 200;
// The runner re-reads the bets every 10 s (the panel and phone every 30 s): the taker
// sizes against these, and the venues' own polls are 10-20 s, so a fill is held fast.
const BETS_POLL_MS = 10 * 1000;
// A phone setting change shows within one scanner tick.
const SETTINGS_POLL_MS = 10 * 1000;
// Snapshots land every 1-5 min per league, so a 30 s pass misses no change.
const CLOSING_FAIRS_POST_MS = 30 * 1000;
const SERVICE_TIMEOUT_MS = 10 * 1000;
const LOOPBACK_HOST_NAMES = ["127.0.0.1", "localhost"];
const WILDCARD_HOSTS = ["0.0.0.0", "::", ""];
const HTTP_DEFAULT_PORT = 80;
const MAX_PORT = 65535;
const HOUR_MS = 3600 * 1000;
// A started game keeps its card this long after its start, whether or not
// it is still on the board: past any game's length plus a slow grader's lag.
const KICKOFF_MEMO_HOURS = 12;
// A failed crosswalk or fill-fair write waits this long before the next try (panel.js).
const SERVICE_WRITE_RETRY_MS = 60 * 1000;
const HTTP_BAD_REQUEST = 400;

function log(message) {
  console.log(`${new Date().toISOString()} [runner] ${message}`);
}

function logError(message) {
  console.error(`${new Date().toISOString()} [runner] ${message}`);
}

// The runner's config from the environment; throws naming the bad variable.
//   {host, port, betsServiceUrl, kalshiBaseUrl, novigBaseUrl}
function configFromEnv(env) {
  const host = env.UNABATED_RUNNER_HOST || DEFAULT_HOST;
  const rawPort = env.UNABATED_RUNNER_PORT || String(DEFAULT_PORT);
  const port = Number(rawPort);
  if (!Number.isInteger(port) || port < 1 || port > MAX_PORT) {
    throw new Error(`UNABATED_RUNNER_PORT must be a whole number from 1 to ${MAX_PORT}, got ${JSON.stringify(rawPort)}`);
  }
  const betsServiceUrl = (env.BETS_SERVICE_URL || DEFAULT_BETS_SERVICE_URL).replace(/\/+$/, "");
  if (!/^https?:\/\/\S+$/.test(betsServiceUrl)) {
    throw new Error(`BETS_SERVICE_URL must be an http(s) URL, got ${JSON.stringify(env.BETS_SERVICE_URL)}`);
  }
  const kalshiBaseUrl = (env.KALSHI_BASE_URL || DEFAULT_KALSHI_BASE_URL).replace(/\/+$/, "");
  if (!/^https?:\/\/\S+$/.test(kalshiBaseUrl)) {
    throw new Error(`KALSHI_BASE_URL must be an http(s) URL, got ${JSON.stringify(env.KALSHI_BASE_URL)}`);
  }
  const novigBaseUrl = (env.NOVIG_API_URL || DEFAULT_NOVIG_BASE_URL).replace(/\/+$/, "");
  if (!/^https?:\/\/\S+$/.test(novigBaseUrl)) {
    throw new Error(`NOVIG_API_URL must be an http(s) URL, got ${JSON.stringify(env.NOVIG_API_URL)}`);
  }
  return { host, port, betsServiceUrl, kalshiBaseUrl, novigBaseUrl };
}

// The Host values a request to this runner can carry: the loopback names and
// the address it is bound to, each with the port (a browser omits the
// scheme's default port, so on 80 the bare names too). Same rule as the bets
// service (service.py allowed_hosts): a DNS-rebound page sends its own name.
function allowedHosts(bindHost, port) {
  const names = LOOPBACK_HOST_NAMES.slice();
  const bound = String(bindHost).toLowerCase();
  if (!WILDCARD_HOSTS.includes(bound) && !names.includes(bound)) names.push(bound.includes(":") ? `[${bound}]` : bound);
  const withPort = names.map((name) => `${name}:${port}`);
  return port === HTTP_DEFAULT_PORT ? withPort.concat(names) : withPort;
}

function hostAllowed(hostHeader, allowed) {
  return allowed.includes(String(hostHeader || "").trim().toLowerCase());
}

// The leagues the scanner loads: the Edges list's, plus NFL and CFB, as the
// panel's scannerLeaguesOf does — the board the open teasers and the bets are
// matched against must be the panel's, or the stakes would differ — plus the
// open bets' leagues, so each has a closing fair. A league outside the
// settings never reaches the Edges list (edgerows.js leagueIds).
function scannerLeaguesOf(edgeLeagues, betLeagueIds = []) {
  return Array.from(new Set([...edgeLeagues, ...teaser.TEASER_LEAGUE_IDS, ...betLeagueIds])).sort((a, b) => a - b);
}

// scanner.js passes the browser's fetch options ({cache: "no-cache",
// credentials: "include"}); Node has no HTTP cache and no cookie jar, so
// only the abort signal means anything here.
function nodeFetch(url, options) {
  return fetch(url, { signal: options && options.signal });
}

// One JSON GET of the bets service. Node's fetch says only "fetch failed" for
// a refused or reset connection; the reason is in its cause, so it is named.
async function getServiceJson(fetchImpl, url, init = {}) {
  let response;
  try {
    response = await fetchImpl(url, { ...init, signal: AbortSignal.timeout(SERVICE_TIMEOUT_MS) });
  } catch (error) {
    const reason = error.cause ? ` (${error.cause.code || error.cause.message})` : "";
    throw new Error(`${url}: ${error.message}${reason}`, { cause: error });
  }
  if (!response.ok) throw new Error(`HTTP ${response.status} from ${url}`);
  return response.json();
}

// One JSON request to a bets-service write route: the parsed reply, or an
// Error carrying the HTTP status and the service's own error text.
async function sendServiceJson(fetchImpl, url, method, body, timeoutMs = SERVICE_TIMEOUT_MS) {
  let response;
  try {
    response = await fetchImpl(url, {
      method, headers: { "Content-Type": "application/json" }, body: JSON.stringify(body), signal: AbortSignal.timeout(timeoutMs),
    });
  } catch (error) {
    const reason = error.cause ? ` (${error.cause.code || error.cause.message})` : "";
    throw new Error(`${method} ${url}: ${error.message}${reason}`, { cause: error });
  }
  const reply = await response.json().catch(() => null);
  if (!response.ok) {
    const error = new Error(`HTTP ${response.status} from ${method} ${url}${reply && reply.error ? `: ${reply.error}` : ""}`);
    error.status = response.status;
    throw error;
  }
  return reply;
}

// The scan loop, the bets poll and the settings poll, held in memory.
//   deps  {fetchImpl (Unabated snapshots), serviceFetch (the bets service),
//          kalshiFetch (Kalshi's orderbooks, default serviceFetch), kalshiBaseUrl,
//          novigFetch (Novig's public catalog, default serviceFetch), novigBaseUrl,
//          betsServiceUrl, now, timers {setInterval, clearInterval},
//          logInfo / logWarning (message) -> void, default stdout / stderr}
// Returns {start, stop, pollBets, pollSettings, scanLoaded, learn, takerPass,
//          postTakerDecisions, edgesPayload, teasersPayload, attachPayload,
//          scenariosPayload, takerPayload, health, scanner, boardLines (the
//          board's rows, read by the tests)}.
function createRunner(deps) {
  const now = deps.now || (() => Date.now());
  const timers = deps.timers || { setInterval: (fn, ms) => setInterval(fn, ms), clearInterval: (id) => clearInterval(id) };
  const serviceFetch = deps.serviceFetch || deps.fetchImpl;
  const kalshiFetch = deps.kalshiFetch || serviceFetch;
  const kalshiBaseUrl = deps.kalshiBaseUrl || DEFAULT_KALSHI_BASE_URL;
  const novigFetch = deps.novigFetch || serviceFetch;
  const novigBaseUrl = deps.novigBaseUrl || DEFAULT_NOVIG_BASE_URL;
  const logInfo = deps.logInfo || log;
  const logWarning = deps.logWarning || logError;
  const startedAt = now();

  let feedState = null;
  let scannerStatus = null;
  let history = {};
  let boardLinesCache = null;
  // leagueId -> {version, rows}: each league's board rows as of its last
  // snapshot (scanner.getLeagueVersion; boardLines rebuilds only the league that landed).
  const boardRowsByLeague = new Map();
  // leagueId -> the league version rememberKickoffOdds last read.
  const kickoffLeagueVersions = new Map();
  let ladderReaders = null;
  // The panel's other per-snapshot caches: the tail-flex measurement (dropped
  // on every scanner update, re-measured when the max line age changes) and
  // teaser.teaserBoardOf's pass, redone only when an NFL or CFB snapshot lands.
  let tailFlexCache = null;
  let teaserBoard = null;
  let teaserBoardLoadedAt = null;
  let teamsSpellingCount = teams.spellingCount();
  let held = { records: [], crosswalk: [], pins: [], fillFairs: [] };
  // Bet ids the tracker's Bets page removed (not Cal's): never on a card.
  let excludedBetIds = new Set();
  // eventId -> {row, lines, oddsAt}: each game's board row and its lines as
  // of the last snapshot before it started (lines null when the runner never
  // saw it pregame). The Live tab's chances are these kickoff odds.
  const kickoffMemo = new Map();
  let fillFairIndex = new Map();
  // {betId: startMs} of the board game each open bet last matched (the
  // panel's knownStarts): a bet with no start of its own whose game has
  // started stops flagging as needing a game.
  let knownStarts = {};
  // The shared Dismiss and Can't tease marks (bets.duckdb::bet_dismissals,
  // teaser_blocks), as /bets.json last served them.
  let dismissedIds = [];
  let teaserBlocks = [];
  // The Edges lines Skipped on either page (bets.duckdb::edge_skips).
  let edgeSkips = [];
  // eventTicker -> {startMs, expiresAt, error?}: Kalshi's start per game (kalshiStartOf).
  const kalshiStarts = new Map();
  // The Teasers list's last build (teaser.planTeasers `previous`): the list
  // is rebuilt only when what it is built on changed, so it holds still.
  let teaserBuild = null;
  let teaserLastError = null;
  // Learning (panel.js learnCrosswalk / captureFillFairs): one pass at a time,
  // a minute's wait after a failed write, each conflict logged once, each bet
  // looked at until its fill fair is decided, unsent rows kept for the retry.
  let learnBusy = false;
  let crosswalkRetryAt = 0;
  let crosswalkLastError = null;
  const crosswalkConflictsLogged = new Set();
  let fillFairsRetryAt = 0;
  let fillFairsLastError = null;
  const fillFairsDecided = new Set();
  const fillFairsUnsent = new Map();
  let fillBooksRetryAt = 0;
  let fillBooksLastError = null;
  const fillBooksDecided = new Set();
  const fillBooksUnsent = new Map();
  const betsStatus = { okAt: null, error: null, unreachableSince: null, generatedAt: null, sources: {} };
  let settings = { ...edgesPayload.settingsFromService(null), source: "defaults", error: null, okAt: null, updatedAt: null };
  let scannedLeagues = null;
  // The scanner's current full load (scanner.start), resolved when it lands.
  let scanLoad = Promise.resolve();
  let pollTimers = [];
  let betsBusy = false;
  let settingsBusy = false;
  let closingBusy = false;
  // betId -> what was last POSTed for it (closefair.unsentRows).
  const sentClosingFairs = new Map();
  const closingStatus = { okAt: null, error: null, rowsSent: 0 };
  // Leagues added for open bets' closing fairs. It only grows (until a
  // restart): every change restarts the scan, which wipes the board and the
  // line history, so a league is not dropped when its last bet settles.
  const betLeagueIds = new Set();
  // The automated taker (server/taker.js): the last pass's decisions, each
  // row's last logged signature (a row logs again only when its decision
  // changes), and the decisions the bets service has not taken yet.
  let takerBusy = false;
  let takerPostBusy = false;
  // `skipped` is the Kalshi side's (as before Novig), `novigSkipped` the Novig side's.
  let takerLastPass = { startedAt: null, finishedAt: null, skipped: "not run yet", novigSkipped: "not run yet", rows: 0, byReason: {}, coverage: {} };
  let takerDecisions = [];
  const takerLogged = new Map();
  let takerUnsent = [];
  let takerDeskStatus = "not asked yet";
  let novigDeskStatus = "not asked yet";
  // Novig's catalog as last read: {leagues: Set, readAt, index (novigOutcomeIndex)}.
  let novigCatalog = null;
  // How long the event loop stalled since the last Novig quote (logged with it):
  // a quote slow on a stalled loop is the runner's own work, not Novig.
  const loopDelay = monitorEventLoopDelay({ resolution: 20 });
  loopDelay.enable();
  const takerLogStatus = { okAt: null, error: null, rowsSent: 0, dropped: 0 };

  // Every snapshot carries Unabated's team list; register it so bet records
  // resolve to the board's team ids, and re-resolve when it grew (as the panel does).
  function registerFeedTeams() {
    if (!feedState || !feedState.teamIndex) return;
    for (const [league, list] of Object.entries(edgeRows.feedTeamsByLeague(feedState))) teams.registerTeams(league, list);
    const spellings = teams.spellingCount();
    if (spellings === teamsSpellingCount) return;
    teamsSpellingCount = spellings;
    held = { ...held, records: betsLib.resolveTeamKeys(held.records, held.crosswalk) };
  }

  const scanner = scannerLib.createScanner({
    fetchImpl: deps.fetchImpl, now, timers, maxRefreshMs: TAKER_SNAPSHOT_REFRESH_MS,
    onChange: (status, state, lineHistory) => {
      if (status.error && (!scannerStatus || scannerStatus.error !== status.error)) logWarning(`scanner: ${status.error}`);
      if (status.phase === "live" && (!scannerStatus || scannerStatus.phase !== "live")) {
        logInfo(`scanner live: ${status.leaguesLoaded.length} leagues, ${status.lineCount} lines (+${status.altLineCount} alts)`);
      }
      scannerStatus = status;
      feedState = state;
      history = lineHistory || {};
      boardLinesCache = null;
      ladderReaders = null;
      tailFlexCache = null;
      registerFeedTeams();
      noteMatchedStarts();
      rememberKickoffOdds();
      learn();
    },
  });

  // Keep every board game's row, and its lines while it has not started, so
  // a game in progress keeps its card and its kickoff chances whether or not
  // the board still lists it. Games started more than KICKOFF_MEMO_HOURS ago go.
  // Only leagues whose snapshot landed since the last call are looked at: an
  // unchanged league's rows and lines are the ones already remembered.
  function rememberKickoffOdds() {
    const at = now();
    const leagueIds = boardLeagueIds();
    for (const leagueId of leagueIds) {
      const version = scanner.getLeagueVersion(leagueId);
      if (kickoffLeagueVersions.get(leagueId) === version) continue;
      kickoffLeagueVersions.set(leagueId, version);
      const linesByEvent = ladderLib.groupLinesByEvent(scanner.getLeagueLines(leagueId));
      for (const row of leagueBoardRows(leagueId)) {
        if (typeof row.eventStartMs !== "number") continue;
        const pregame = row.eventStartMs > at;
        const known = kickoffMemo.get(row.eventId);
        if (pregame) kickoffMemo.set(row.eventId, { row, lines: linesByEvent.get(row.eventId) || [], oddsAt: at });
        else if (!known) kickoffMemo.set(row.eventId, { row, lines: null, oddsAt: null });
      }
    }
    for (const leagueId of kickoffLeagueVersions.keys()) if (!leagueIds.includes(leagueId)) kickoffLeagueVersions.delete(leagueId);
    for (const [eventId, entry] of kickoffMemo) {
      if (at - entry.row.eventStartMs > KICKOFF_MEMO_HOURS * HOUR_MS) kickoffMemo.delete(eventId);
    }
  }

  // Remember each open bet's matched game start, as the panel's noteMatchedStarts does.
  function noteMatchedStarts() {
    knownStarts = betsView.keepKnownStartsOpen(knownStarts, betsLib.matchedStarts(held.records, boardLines()), held.records);
  }

  // The board's leagues in board order (state.leagues: a refreshed league moves to the end, as its lines do).
  function boardLeagueIds() {
    return feedState ? Array.from(new Set(feedState.leagues)) : [];
  }

  // One league's board rows (edgeRows.boardLines over its lines), rebuilt only when its snapshot landed again.
  // A row keeps the book and team names of its own league's last snapshot (accepted: rows only name games).
  function leagueBoardRows(leagueId) {
    const version = scanner.getLeagueVersion(leagueId);
    const cached = boardRowsByLeague.get(leagueId);
    if (cached && cached.version === version) return cached.rows;
    const rows = edgeRows.boardLines(feedState, scanner.getLeagueLines(leagueId));
    boardRowsByLeague.set(leagueId, { version, rows });
    return rows;
  }

  // Every league's rows, in the order a pass over state.lines gives them
  // (lines sit in league merge order), one row per event.
  function boardLines() {
    if (!feedState) return [];
    if (boardLinesCache) return boardLinesCache;
    const leagueIds = boardLeagueIds();
    for (const leagueId of boardRowsByLeague.keys()) if (!leagueIds.includes(leagueId)) boardRowsByLeague.delete(leagueId);
    const seen = new Set();
    const rows = [];
    for (const leagueId of leagueIds) {
      for (const row of leagueBoardRows(leagueId)) {
        if (seen.has(row.eventId)) continue;
        seen.add(row.eventId);
        rows.push(row);
      }
    }
    boardLinesCache = rows;
    return rows;
  }

  function ladderReaderOf(eventId) {
    if (!ladderReaders) ladderReaders = edgeRows.createLadderReaders(feedState);
    return ladderReaders(eventId);
  }

  // Measured once per scanner update, on the same line-age gate as the list (panel.js tailFlexMeasurement).
  function tailFlexMeasurement() {
    const maxLineAgeMs = settings.edgeSettings.maxLineAgeHours * HOUR_MS;
    if (!tailFlexCache || tailFlexCache.maxLineAgeMs !== maxLineAgeMs) {
      tailFlexCache = { maxLineAgeMs, measurement: tailflex.measureTailFlex(feedState, { now: now(), maxLineAgeMs }) };
    }
    return tailFlexCache.measurement;
  }

  // The board pass, redone when an NFL or CFB snapshot has landed since (panel.js currentTeaserBoard).
  function currentTeaserBoard() {
    const loadedAt = teaser.TEASER_LEAGUE_IDS.map((id) => (scannerStatus && scannerStatus.leagueLoadedAt[id]) || 0).join(",");
    if (!teaserBoard || loadedAt !== teaserBoardLoadedAt) {
      teaserBoard = teaser.teaserBoardOf(feedState);
      teaserBoardLoadedAt = loadedAt;
    }
    return teaserBoard;
  }

  // The open BFA teasers, each leg joined to its board game and priced
  // (panel.js openTeasersNow): the Edges rows are sized against them too.
  function openTeasersNow() {
    if (!feedState) return [];
    return teaser.openTeasers(held.records, boardLines(), { now: now(), ladderOf: currentTeaserBoard().ladderOf });
  }

  // One GET of /bets.json, applied exactly as the panel applies it. A failure
  // keeps the held bets and records since when the service has been unreachable.
  async function pollBets() {
    if (betsBusy) return;
    betsBusy = true;
    const at = now();
    try {
      const applied = edgeRows.applyBetsPayload(held, await getServiceJson(serviceFetch, `${deps.betsServiceUrl}/bets.json`), at);
      held = { records: applied.records, crosswalk: applied.crosswalk, pins: applied.pins, fillFairs: applied.fillFairs };
      if (Array.isArray(applied.exclusions)) excludedBetIds = new Set(applied.exclusions.map((row) => row.betId));
      if (Array.isArray(applied.dismissals)) dismissedIds = applied.dismissals.map((row) => row.betId);
      if (Array.isArray(applied.teaserBlocks)) teaserBlocks = applied.teaserBlocks;
      if (Array.isArray(applied.edgeSkips)) edgeSkips = applied.edgeSkips;
      fillFairIndex = fillfair.fairsByBetId(held.fillFairs);
      noteMatchedStarts();
      // Before the first settings read the scan has not started; start() starts it.
      if (scannedLeagues !== null) followLeagues();
      if (betsStatus.error) logInfo("bets service reachable again");
      Object.assign(betsStatus, { okAt: at, error: null, unreachableSince: null, generatedAt: applied.generatedAt, sources: applied.sources });
    } catch (error) {
      if (betsStatus.error !== error.message) logWarning(`bets poll failed: ${error.message}`);
      Object.assign(betsStatus, { error: error.message, unreachableSince: betsStatus.unreachableSince ?? at });
    } finally {
      betsBusy = false;
    }
    learn();
  }

  // Crosswalk rows, then fill fairs, then other books at the fill; never
  // throws (a failure is logged and retried later).
  function learn() {
    if (learnBusy || !feedState) return;
    learnBusy = true;
    learnCrosswalk()
      .then(() => captureFillFairs())
      .then(() => captureFillBooks())
      .catch((error) => logWarning(`learning failed: ${error.stack || error.message}`))
      .finally(() => { learnBusy = false; });
  }

  // The open bets the board joined by venue id teach what their venue calls
  // both teams; the service's reply is the whole table, and the records are
  // re-keyed through it (panel.js learnCrosswalk + applyServiceTables).
  async function learnCrosswalk() {
    if (now() < crosswalkRetryAt) return;
    const { learned, conflicts } = betsLib.learnCrosswalk(held.records, boardLines(), held.crosswalk);
    for (const conflict of conflicts) {
      const tag = `${conflict.betId}:${conflict.side}`;
      if (crosswalkConflictsLogged.has(tag)) continue;
      crosswalkConflictsLogged.add(tag);
      logInfo(`crosswalk: not learning ${conflict.venue} ${conflict.league} "${conflict.venueTeamKey}" from ${conflict.betId}: ${conflict.reason}`);
    }
    if (!learned.length) return;
    try {
      const reply = await sendServiceJson(serviceFetch, `${deps.betsServiceUrl}/crosswalk.json`, "POST", { rows: learned });
      if (!reply || !Array.isArray(reply.crosswalk)) throw new Error("crosswalk.json reply has no crosswalk array");
      logInfo(`crosswalk: learned ${reply.learned} row(s), ${(reply.conflicts || []).length} refused by the service, ${reply.crosswalk.length} held`);
      held = { ...held, crosswalk: reply.crosswalk, records: betsLib.rekeyRecords(betsLib.applyPins(held.records, held.pins), reply.crosswalk) };
      crosswalkLastError = null;
    } catch (error) {
      crosswalkRetryAt = now() + SERVICE_WRITE_RETRY_MS;
      if (crosswalkLastError !== error.message) logWarning(`crosswalk: service write failed: ${error.message}`);
      crosswalkLastError = error.message;
    }
  }

  // Each new open bet's fair when it was placed, read off the line history
  // while the scanner watched through the fill (fillfair.js), one row per
  // POST: the service refuses a whole request on its first bad row. A 400
  // drops the row; any other failure keeps every unsent row for a minute.
  async function captureFillFairs() {
    const liveStatus = scanner.getStatus();
    const { saves, refusals } = fillfair.captureFillFairs({
      records: held.records, skipIds: new Set([...fillFairIndex.keys(), ...fillFairsDecided]),
      state: feedState, boardLines: boardLines(), history,
      observingSince: liveStatus.observingSince, leagueObservingSince: liveStatus.leagueObservingSince, now: now(),
    });
    for (const save of saves) {
      fillFairsDecided.add(save.betId);
      fillFairsUnsent.set(save.betId, save);
    }
    for (const refusal of refusals) fillFairsDecided.add(refusal.betId);
    if (saves.length) logInfo(`fill fairs: ${saves.length} captured`);
    if (fillFairsUnsent.size === 0 || now() < fillFairsRetryAt) return;
    for (const row of Array.from(fillFairsUnsent.values())) {
      try {
        const reply = await sendServiceJson(serviceFetch, `${deps.betsServiceUrl}/fill_fairs.json`, "POST", { rows: [row] });
        if (!reply || !Array.isArray(reply.fillFairs)) throw new Error("fill_fairs.json reply has no fillFairs array");
        fillFairsUnsent.delete(row.betId);
        fillFairsLastError = null;
        held = { ...held, fillFairs: fillfair.mergeFillFairs(held.fillFairs, reply.fillFairs, held.records) };
        fillFairIndex = fillfair.fairsByBetId(held.fillFairs);
      } catch (error) {
        if (error.status === HTTP_BAD_REQUEST) {
          fillFairsUnsent.delete(row.betId);
          logWarning(`fill fairs: the service refused ${row.betId}, dropped: ${error.message}`);
          continue;
        }
        fillFairsRetryAt = now() + SERVICE_WRITE_RETRY_MS;
        if (fillFairsLastError !== error.message) logWarning(`fill fairs: service write failed: ${error.message}`);
        fillFairsLastError = error.message;
        return;
      }
    }
  }

  // Each new open bet's other books at its number when it was placed
  // (fillbooks.js), one row per POST, the same way as the fill fairs. Decided
  // ids live in memory only: after a restart every open bet was placed before
  // the scanner watched, so nothing is re-captured, and the service keeps the
  // first row regardless.
  async function captureFillBooks() {
    const liveStatus = scanner.getStatus();
    const { saves, refusals } = fillbooks.captureFillBooks({
      records: held.records, skipIds: fillBooksDecided,
      state: feedState, boardLines: boardLines(), history,
      observingSince: liveStatus.observingSince, leagueObservingSince: liveStatus.leagueObservingSince, now: now(),
    });
    for (const save of saves) {
      fillBooksDecided.add(save.betId);
      fillBooksUnsent.set(save.betId, save);
    }
    for (const refusal of refusals) fillBooksDecided.add(refusal.betId);
    if (saves.length) logInfo(`fill books: ${saves.length} captured`);
    if (fillBooksUnsent.size === 0 || now() < fillBooksRetryAt) return;
    for (const row of Array.from(fillBooksUnsent.values())) {
      try {
        const reply = await sendServiceJson(serviceFetch, `${deps.betsServiceUrl}/fill_books.json`, "POST", { rows: [row] });
        if (!reply || reply.ok !== true) throw new Error("fill_books.json reply is not ok");
        fillBooksUnsent.delete(row.betId);
        fillBooksLastError = null;
      } catch (error) {
        if (error.status === HTTP_BAD_REQUEST) {
          fillBooksUnsent.delete(row.betId);
          logWarning(`fill books: the service refused ${row.betId}, dropped: ${error.message}`);
          continue;
        }
        fillBooksRetryAt = now() + SERVICE_WRITE_RETRY_MS;
        if (fillBooksLastError !== error.message) logWarning(`fill books: service write failed: ${error.message}`);
        fillBooksLastError = error.message;
        return;
      }
    }
  }

  // The scanner follows the settings' leagues (plus NFL and CFB); a change
  // restarts it, as the panel's league checkboxes do. Not awaited here: a
  // full load takes a while and must not hold the settings poll
  // (scanLoaded() waits for it).
  function followLeagues() {
    for (const id of closefair.leagueIdsOfOpenBets(held.records)) betLeagueIds.add(id);
    const leagues = scannerLeaguesOf(settings.edgeSettings.leagues, betLeagueIds);
    const signature = leagues.join(",");
    if (signature === scannedLeagues) return;
    scannedLeagues = signature;
    logInfo(`scanning ${leagues.length} league(s)`);
    scanLoad = scanner.start(leagues).catch((error) => logWarning(`scanner start failed: ${error.message}`));
  }

  // One GET of /settings.json. A failure keeps the last settings read (the
  // panel's defaults until one succeeds) and says so in /edges.json.
  async function pollSettings() {
    if (settingsBusy) return;
    settingsBusy = true;
    try {
      const body = await getServiceJson(serviceFetch, `${deps.betsServiceUrl}/settings.json`);
      if (!body || typeof body.settings !== "object" || body.settings === null) throw new Error("settings.json has no settings object");
      if (settings.error) logInfo("settings readable again");
      settings = { ...edgesPayload.settingsFromService(body.settings), source: "service", error: null, okAt: now(), updatedAt: body.updatedAt ?? null };
    } catch (error) {
      if (settings.error !== error.message) logWarning(`settings poll failed (keeping the ${settings.source} settings): ${error.message}`);
      settings = { ...settings, error: error.message };
    } finally {
      settingsBusy = false;
    }
    followLeagues();
  }

  // POST the closing-fair rows that changed since the last successful POST.
  // A failure keeps them unsent, so the next pass retries them.
  async function postClosingFairs() {
    if (closingBusy || !feedState || !scannerStatus) return;
    closingBusy = true;
    try {
      const rows = closefair.unsentRows(closefair.closingFairRows({
        records: held.records, state: feedState, boardLines: boardLines(),
        leagueLoadedAt: scannerStatus.leagueLoadedAt || {}, staleLeagues: scannerStatus.staleLeagues, now: now(),
      }), sentClosingFairs);
      closefair.forgetClosed(sentClosingFairs, held.records);
      if (!rows.length) return;
      await getServiceJson(serviceFetch, `${deps.betsServiceUrl}/closing_fairs.json`, {
        method: "POST", headers: { "Content-Type": "application/json" }, body: JSON.stringify({ rows }),
      });
      closefair.markSent(rows, sentClosingFairs);
      if (closingStatus.error) logInfo("closing fairs saving again");
      Object.assign(closingStatus, { okAt: now(), error: null, rowsSent: closingStatus.rowsSent + rows.length });
    } catch (error) {
      if (closingStatus.error !== error.message) logWarning(`closing fairs POST failed: ${error.message}`);
      closingStatus.error = error.message;
    } finally {
      closingBusy = false;
    }
  }

  // One Kalshi orderbook, parsed (taker.parseKalshiOrderbook). Throws naming
  // the ticker and what went wrong; a 404 is a ticker Kalshi does not list.
  async function readKalshiBook(ticker) {
    const body = await getServiceJson(kalshiFetch, `${kalshiBaseUrl}/markets/${encodeURIComponent(ticker)}/orderbook`);
    return taker.parseKalshiOrderbook(body);
  }

  // id -> {book} | {error} for every id, KALSHI_BOOK_READS_AT_ONCE reads of
  // `readOne(id)` in flight (Kalshi tickers and Novig market ids alike).
  async function readBooks(ids, readOne) {
    const books = new Map();
    const queue = Array.from(ids);
    const worker = async () => {
      while (queue.length) {
        const id = queue.shift();
        try {
          books.set(id, { book: await readOne(id) });
        } catch (error) {
          books.set(id, { error: error.message });
        }
      }
    };
    await Promise.all(Array.from({ length: KALSHI_BOOK_READS_AT_ONCE }, worker));
    return books;
  }

  // Why the taker must not decide on any venue right now, or null: no live
  // board, or the bets service unread or unreachable. Fail closed.
  function boardBlocker() {
    if (!feedState || !scannerStatus || scannerStatus.phase !== "live") return "the board is not live";
    if (betsStatus.okAt === null) return "the bets service has not been read yet";
    if (betsStatus.error) return `the bets service is unreachable: ${betsStatus.error}`;
    return null;
  }

  // Why the Novig side must not decide, or null: the held Novig bets are
  // unknown (Cal's app account or the bot's subaccount unread, failed or too
  // old), or the bets service reads no bot subaccount (no bot key: no desk).
  function novigBlocker() {
    for (const name of NOVIG_HOLDINGS_SOURCES) {
      const source = betsStatus.sources && betsStatus.sources[name];
      if (!source) return `the bets service does not read ${name}`;
      if (source.ok !== true) return `the bets service's ${name} read failed: ${source.error || "no detail"}`;
      const fetchedAt = Date.parse(source.fetchedAt || "");
      if (!Number.isFinite(fetchedAt) || now() - fetchedAt > TAKER_MAX_HELD_AGE_MS) {
        return `the bets service's ${name} read is older than ${TAKER_MAX_HELD_AGE_MS / 60000} min (${source.fetchedAt || "never"})`;
      }
    }
    return null;
  }

  // Why the Kalshi side must not decide right now, or null: boardBlocker, or
  // the held Kalshi bets are unknown (its poll failed or is too old).
  function takerBlocker() {
    const board = boardBlocker();
    if (board) return board;
    const kalshiSource = betsStatus.sources && betsStatus.sources.kalshi;
    if (!kalshiSource) return "the bets service does not poll Kalshi, so the held Kalshi bets are unknown";
    if (kalshiSource.ok !== true) return `the bets service's Kalshi poll failed: ${kalshiSource.error || "no detail"}`;
    const fetchedAt = Date.parse(kalshiSource.fetchedAt || "");
    if (!Number.isFinite(fetchedAt) || now() - fetchedAt > TAKER_MAX_HELD_AGE_MS) {
      return `the bets service's Kalshi poll is older than ${TAKER_MAX_HELD_AGE_MS / 60000} min (${kalshiSource.fetchedAt || "never"})`;
    }
    return null;
  }

  // Why a row's league snapshot is not one to trade on, or null: the
  // scanner marked it stale or failed, or it landed too long ago.
  function staleSnapshotReason(leagueId, passStartedAt) {
    if ((scannerStatus.staleLeagues || []).includes(leagueId)) return `league ${leagueId} snapshot is a stale build`;
    if (scannerStatus.leagueErrors && scannerStatus.leagueErrors[leagueId]) return `league ${leagueId} failed to load: ${scannerStatus.leagueErrors[leagueId]}`;
    const loadedAt = scannerStatus.leagueLoadedAt && scannerStatus.leagueLoadedAt[leagueId];
    if (typeof loadedAt !== "number") return `league ${leagueId} has not loaded`;
    const ageMs = passStartedAt - loadedAt;
    if (ageMs > TAKER_MAX_SNAPSHOT_AGE_MS) return `league ${leagueId} snapshot landed ${Math.round(ageMs / 1000)} s ago`;
    return null;
  }

  // One taker pass (server/taker.js): every Kalshi row at the bot's minimum
  // edge, its contract, its live book, its decision; each "bet" goes to the
  // desk. The decisions that changed since they were last logged are queued
  // for the bets service.
  async function takerPass() {
    if (takerBusy) return;
    takerBusy = true;
    const startedAt = now();
    try {
      const blocker = takerBlocker();
      // A Novig stake is sized on every held bet, Cal's Kalshi bets included: it needs the Kalshi poll fresh too.
      const novigSkip = blocker || novigBlocker();
      if (blocker && novigSkip) {
        takerLastPass = { ...takerLastPass, startedAt, finishedAt: now(), skipped: blocker, novigSkipped: novigSkip };
        return;
      }
      // The board as of the pass start: a snapshot landing during the book reads must not mix in.
      const board = feedState;
      const pass = {
        board, startedAt, priceIndex: pricecheck.buildPriceIndex(board, { now: startedAt, maxLineAgeMs: taker.MAX_LINE_AGE_MS }),
        // The same held bets the Edges list sizes against: Cal's Removed bets (bet_exclusions) are not his.
        heldRecords: held.records.filter((record) => !excludedBetIds.has(record.id)), teasers: openTeasersNow(),
      };
      const kalshiPart = blocker ? { rows: [], decisions: [] } : await kalshiDecisions(pass);
      const novigPart = novigSkip ? { rows: [], decisions: [] } : await novigDecisions(pass);
      takerDecisions = kalshiPart.decisions.concat(novigPart.decisions);
      const rows = kalshiPart.rows.concat(novigPart.rows);
      await sendTakerOrders(takerDecisions, rows, startedAt);
      queueChangedDecisions(takerDecisions);
      takerLastPass = {
        startedAt, finishedAt: now(), skipped: blocker, novigSkipped: novigSkip, rows: rows.length,
        byReason: countBy(takerDecisions, (decision) => decision.reason),
        coverage: takerCoverage(takerDecisions),
      };
    } catch (error) {
      logWarning(`taker pass failed: ${error.stack || error.message}`);
      takerLastPass = { ...takerLastPass, startedAt, finishedAt: now(), skipped: `the pass failed: ${error.message}`, novigSkipped: `the pass failed: ${error.message}` };
    } finally {
      takerBusy = false;
    }
    await postTakerDecisions();
  }

  // One venue's rows decided: the contract each buys, its live book, the
  // stake conditional on the held bets. `resolve(row)` gives the contract,
  // `quoteFor(contract)` the live {ask, depth} or {error}.
  function decideRows(pass, rows, venue, resolve, quoteFor) {
    const contracts = rows.map((row) => {
      const stale = staleSnapshotReason(row.leagueId, pass.startedAt);
      return stale ? { reason: "stale_feed", detail: stale } : resolve(row);
    });
    const flags = betsLib.annotateRows(rows, pass.heldRecords, { lines: boardLines() });
    const decidedAt = now();
    const decisions = rows.map((row, index) => {
      const contract = contracts[index];
      const read = contract.reason ? null : quoteFor(contract);
      const decision = taker.decideRow({
        row, contract, quote: read && !read.error ? read.quote : null, quoteError: read ? read.error || null : null,
        priceIndex: pass.priceIndex, venue,
        sizing: { matches: flags[index].matches, ladderOf: ladderReaderOf(row.eventId), teasers: pass.teasers },
      });
      // The held bets the stake was sized on, by id: the desk refuses an add
      // on a game whose earlier bot fill (either venue) is not among them.
      const heldIds = flags[index].matches.map((match) => match.bet.id);
      const extra = venue === taker.VENUE_NOVIG ? { line: taker.novigLineOf(row) } : {};
      return { ...takerRowFields(row), decidedAt: new Date(decidedAt).toISOString(), ...decision, ...extra, heldIds };
    });
    return { rows, decisions };
  }

  async function kalshiDecisions(pass) {
    const rows = taker.candidateRows(pass.board, settings.edgeSettings, pass.startedAt);
    const tickers = new Set();
    for (const row of rows) {
      if (staleSnapshotReason(row.leagueId, pass.startedAt)) continue;
      const contract = taker.resolveKalshiContract(row, pass.board);
      if (contract.ticker) tickers.add(contract.ticker);
    }
    const books = await readBooks(tickers, readKalshiBook);
    return decideRows(pass, rows, taker.VENUE_KALSHI, (row) => taker.resolveKalshiContract(row, pass.board), (contract) => {
      const read = books.get(contract.ticker);
      if (!read) return null;
      return read.book ? { quote: taker.quoteOf(read.book, contract.side) } : { error: read.error };
    });
  }

  // The Novig rows: Novig's open pregame markets for their leagues (one
  // catalog read), each row's outcome on them, the books of those markets.
  // A failed catalog read leaves every row no_outcome, so nothing is sent.
  async function novigDecisions(pass) {
    const rows = taker.candidateRows(pass.board, settings.edgeSettings, pass.startedAt, taker.NOVIG_BOOK_ID);
    if (!rows.length) return { rows, decisions: [] };
    let index = new Map();
    let catalogError = null;
    try {
      index = (await novigCatalogFor(rows, 0)).index;
    } catch (error) {
      catalogError = error.message;
    }
    const resolve = (row) => {
      if (catalogError) return { reason: "no_outcome", detail: `Novig's catalog unread: ${catalogError}` };
      return taker.novigContractOf(row, index);
    };
    const marketIds = new Set();
    for (const row of rows) {
      if (staleSnapshotReason(row.leagueId, pass.startedAt)) continue;
      const contract = resolve(row);
      if (contract.marketId) marketIds.add(contract.marketId);
    }
    const books = await readBooks(marketIds, (marketId) => getNovigJson(`${novigMarketUrl(marketId)}/book`));
    return decideRows(pass, rows, taker.VENUE_NOVIG, resolve, (contract) => {
      const read = books.get(contract.marketId);
      if (!read) return null;
      if (read.error) return { error: read.error };
      try {
        return { quote: taker.novigQuoteOf(read.book, contract) };
      } catch (error) {
        return { error: error.message };
      }
    });
  }

  // One public Novig GET as JSON; throws naming the URL and status.
  function getNovigJson(url) {
    return getServiceJson(novigFetch, url, { headers: { Accept: "application/json", "User-Agent": NOVIG_USER_AGENT } });
  }

  // Novig's outcome index for the rows' leagues: the last catalog read when it
  // covers them and is at most maxAgeMs old, else a fresh read (kept for the
  // next quote); maxAgeMs 0 always reads. {index, readAt, reused}.
  async function novigCatalogFor(rows, maxAgeMs) {
    const leagues = rows.map((row) => NOVIG_LEAGUE_OF[row.league]).filter(Boolean);
    const cached = novigCatalog;
    if (maxAgeMs > 0 && leagues.length && cached && now() - cached.readAt <= maxAgeMs && leagues.every((league) => cached.leagues.has(league))) {
      return { index: cached.index, readAt: cached.readAt, reused: true };
    }
    const readAt = now();
    const index = taker.novigOutcomeIndex(await readNovigMarkets(rows));
    novigCatalog = { leagues: new Set(leagues), readAt, index };
    return { index, readAt, reused: false };
  }

  // Every open pregame Novig market in the rows' leagues, every page (never a partial list).
  async function readNovigMarkets(rows) {
    const leagues = Array.from(new Set(rows.map((row) => NOVIG_LEAGUE_OF[row.league]).filter(Boolean))).sort();
    if (!leagues.length) return [];
    const markets = [];
    let after = null;
    for (let page = 0; page < NOVIG_MARKETS_MAX_PAGES; page += 1) {
      const query = new URLSearchParams({ league: leagues.join(","), eventStatus: taker.NOVIG_EVENT_PREGAME, limit: String(NOVIG_MARKETS_PAGE_LIMIT) });
      if (after) query.set("after", after);
      const body = await getNovigJson(`${novigBaseUrl}/v3/public/catalog/markets?${query}`);
      if (!body || !Array.isArray(body.items)) throw new Error("Novig markets: no items list");
      markets.push(...body.items);
      after = body.next || null;
      if (!after) return markets;
    }
    throw new Error(`Novig markets: more than ${NOVIG_MARKETS_MAX_PAGES} pages`);
  }

  // The Ticket's Place button (bets_service/manual_orders.py asks, Cal
  // clicked): the board's one Kalshi line with these fields, the contract it
  // buys, what buying it costs now and whether Kalshi's market is the line's.
  // The board is used only to find the contract, so its age does not matter
  // here: the price Cal acts on is the Ticket's, checked against the live ask
  // by the bets service. A line the board lacks or holds twice, or no
  // contract, answers {ok: false, reason, detail}; bad params are a 400.
  //   params  URLSearchParams {eventId, periodTypeId, betTypeId, sideIndex, points?}
  async function kalshiQuoteNow(params) {
    const spec = kalshiLineSpecOf(params);
    const board = feedState;
    if (!board) return { ok: false, reason: "no_board", detail: "the board has not loaded yet" };
    const found = linesOfEvent(board, spec.eventId).filter((line) => line.bookId === taker.KALSHI_BOOK_ID
      && line.eventId === spec.eventId && line.periodTypeId === spec.periodTypeId && line.betTypeId === spec.betTypeId
      && line.sideIndex === spec.sideIndex && (spec.betTypeId === MONEYLINE_BET_TYPE_ID || (line.points ?? null) === spec.points));
    if (found.length !== 1) {
      return { ok: false, reason: "no_line", detail: `the board holds ${found.length} Kalshi lines for game ${spec.eventId} at that market, side and number` };
    }
    const row = feedLib.describeLine(found[0], board);
    const fields = { eventId: row.eventId, eventStart: row.eventStart, label: `${row.sideLabel} ${row.period}`.trim(), eventName: row.eventName };
    const contract = await kalshiPlaceContractOf(row, board);
    if (contract.reason) return { ok: false, ...fields, reason: contract.reason, detail: contract.detail };
    const position = { ...fields, ticker: contract.ticker, side: contract.side, derivedTicker: contract.derived === true };
    let book;
    let market;
    try {
      book = await readKalshiBook(contract.ticker);
      market = await readKalshiMarket(contract.ticker);
    } catch (error) {
      return { ok: false, ...position, reason: "no_book", detail: error.message };
    }
    const mismatch = taker.marketMismatch(market, row, contract);
    if (mismatch) return { ok: false, ...position, reason: "market_mismatch", detail: `${contract.ticker}: ${mismatch}` };
    const kalshiStartMs = await kalshiStartOf(taker.eventTickerOf(market, contract.ticker));
    Object.assign(position, startFields(row.eventStart, kalshiStartMs));
    const quote = taker.quoteOf(book, contract.side);
    if (!quote) return { ok: false, ...position, reason: "no_book", detail: `nothing rests on the ${contract.side} side of ${contract.ticker}` };
    return {
      ok: true, ...position, ask: quote.ask, allIn: Math.round(taker.kalshiAllInPrice(quote.ask) * 1e6) / 1e6, depth: quote.depth,
      kalshiBaseUrl,
    };
  }

  // The Kalshi contract Place buys: the bot's own reading of the rungs, and
  // for a moneyline whose team code no full-game spread rung gives (a game
  // without one, a half), Kalshi's own market list for the game's winner
  // event (taker.kalshiMoneylineEventOf / kalshiMoneylineFromMarkets).
  async function kalshiPlaceContractOf(row, board) {
    const contract = taker.resolveKalshiContract(row, board);
    if (!contract.reason || row.betTypeId !== MONEYLINE_BET_TYPE_ID) return contract;
    const event = taker.kalshiMoneylineEventOf(row);
    if (event.reason) return { reason: event.reason, detail: `${contract.detail}; ${event.detail}` };
    let body;
    try {
      body = await getServiceJson(kalshiFetch, `${kalshiBaseUrl}/markets?limit=20&event_ticker=${encodeURIComponent(event.eventTicker)}`);
    } catch (error) {
      return { reason: "no_ticker", detail: `${contract.detail}; Kalshi's markets on ${event.eventTicker} unread: ${error.message}` };
    }
    return taker.kalshiMoneylineFromMarkets(row, body && body.markets, event.eventTicker);
  }

  // The Ticket's Place button on Novig (bets_service/novig_place.py asks):
  // the board's one Novig line with these fields, the outcome it buys and the
  // other side's (the taker's own novigContractOf over Novig's open pregame
  // catalog), what buying it costs now (novigQuoteOf) and whether Novig's
  // market and event are the line's and take a free pregame order
  // (novigMarketMismatch). As on Kalshi, the board only finds the outcome;
  // the price Cal acts on is checked by the bets service. Places nothing.
  //   params  URLSearchParams {eventId, periodTypeId, betTypeId, sideIndex, points}
  async function novigQuoteNow(params) {
    const spec = kalshiLineSpecOf(params);
    const startedAt = now();
    const { answer, catalog } = await novigQuoteOfSpec(spec);
    const stalledMs = Math.round(loopDelay.max / 1e6);
    loopDelay.reset();
    logInfo(`novig quote ${spec.eventId} ${answer.label || `bet type ${spec.betTypeId} side ${spec.sideIndex} ${spec.points}`}: `
      + `${answer.ok ? `ok @ ${answer.price}` : answer.reason} in ${now() - startedAt} ms`
      + (catalog ? ` (catalog ${catalog.reused ? `reused, ${Math.round((startedAt - catalog.readAt) / 1000)}s old` : "read"}, `
        + `market reads ${catalog.readsMs} ms)` : "")
      + `; event loop stalled up to ${stalledMs} ms since the last quote`);
    return answer;
  }

  // One game's lines, read from its league's lines only: a Place quote walking
  // the whole ~1.2M-line board took ~1.3 s per click (2026-10-19 review).
  // A game the board does not hold has none. Right after a scan restart the
  // scanner holds a new board while `board` is still the last published one:
  // that board is walked whole, as before, so the quote matches the list Cal sees.
  function linesOfEvent(board, eventId) {
    const event = board.events[eventId];
    if (!event) return [];
    const lines = board === scanner.getState() ? scanner.getLeagueLines(event.leagueId) : Object.values(board.lines);
    return lines.filter((line) => line.eventId === eventId);
  }

  // {answer: novigQuoteNow's answer, catalog: {reused, readAt, readsMs} for its log line once the outcome was found, else null}.
  async function novigQuoteOfSpec(spec) {
    const board = feedState;
    if (!board) return { answer: { ok: false, reason: "no_board", detail: "the board has not loaded yet" }, catalog: null };
    const found = linesOfEvent(board, spec.eventId).filter((line) => line.bookId === taker.NOVIG_BOOK_ID
      && line.eventId === spec.eventId && line.periodTypeId === spec.periodTypeId && line.betTypeId === spec.betTypeId
      && line.sideIndex === spec.sideIndex && (spec.betTypeId === MONEYLINE_BET_TYPE_ID || (line.points ?? null) === spec.points));
    if (found.length !== 1) {
      return { answer: { ok: false, reason: "no_line", detail: `the board holds ${found.length} Novig lines for game ${spec.eventId} at that market, side and number` }, catalog: null };
    }
    const row = feedLib.describeLine(found[0], board);
    const fields = { eventId: row.eventId, eventStart: row.eventStart, label: `${row.sideLabel} ${row.period}`.trim(), eventName: row.eventName };
    let catalog;
    let contract;
    try {
      // A moneyline has no Novig rung of its own: Place finds it off the game's other rungs (never the bot).
      const contractOf = (index) => (row.betTypeId === MONEYLINE_BET_TYPE_ID
        ? taker.novigMoneylineContract(row, index, board) : taker.novigContractOf(row, index));
      catalog = await novigCatalogFor([row], NOVIG_CATALOG_REUSE_MS);
      contract = contractOf(catalog.index);
      // A market listed since the reused read: one fresh read before giving up.
      if (contract.reason && catalog.reused) {
        catalog = await novigCatalogFor([row], 0);
        contract = contractOf(catalog.index);
      }
    } catch (error) {
      return { answer: { ok: false, ...fields, reason: "no_outcome", detail: `Novig's catalog unread: ${error.message}` }, catalog: null };
    }
    if (contract.reason) return { answer: { ok: false, ...fields, reason: contract.reason, detail: contract.detail }, catalog: null };
    const position = { ...fields, outcomeId: contract.outcomeId, otherOutcomeId: contract.otherOutcomeId, marketId: contract.marketId };
    // The book, the market and its event at once (the catalog names the event; the check confirms it is the market's).
    const listedEventId = catalog.index.get(contract.outcomeId).market.eventId;
    const readsStartedAt = now();
    const [bookRead, mismatch] = await Promise.all([
      getNovigJson(`${novigMarketUrl(contract.marketId)}/book`).then((book) => ({ book }), (error) => ({ error })),
      readNovigMarketMismatch(row, contract, listedEventId),
    ]);
    const timing = { reused: catalog.reused, readAt: catalog.readAt, readsMs: now() - readsStartedAt };
    let quote;
    try {
      if (bookRead.error) throw bookRead.error;
      quote = taker.novigQuoteOf(bookRead.book, contract);
    } catch (error) {
      return { answer: { ok: false, ...position, reason: "no_book", detail: error.message }, catalog: timing };
    }
    if (mismatch) return { answer: { ok: false, ...position, reason: "market_mismatch", detail: mismatch }, catalog: timing };
    if (!quote) return { answer: { ok: false, ...position, reason: "no_book", detail: `nothing rests on the other side of Novig market ${contract.marketId}` }, catalog: timing };
    return { answer: { ok: true, ...position, price: quote.ask, depth: quote.depth, novigBaseUrl, line: taker.novigLineOf(row) }, catalog: timing };
  }

  // Kalshi's start for an event as epoch ms, or null when it lists none or
  // the read failed (never throws: option B trades on Unabated's start then).
  async function kalshiStartOf(eventTicker) {
    if (!eventTicker) return null;
    const cached = kalshiStarts.get(eventTicker);
    if (cached && now() < cached.expiresAt) return cached.startMs;
    for (const [key, entry] of kalshiStarts) if (now() >= entry.expiresAt) kalshiStarts.delete(key);
    const url = `${kalshiBaseUrl}/milestones?limit=20&related_event_ticker=${encodeURIComponent(eventTicker)}`;
    try {
      const startMs = taker.kalshiStartMsOf(await getServiceJson(kalshiFetch, url), eventTicker);
      kalshiStarts.set(eventTicker, { startMs, expiresAt: now() + KALSHI_START_TTL_MS });
      return startMs;
    } catch (error) {
      if (!cached || cached.error !== error.message) logWarning(`kalshi start ${eventTicker}: ${error.message}`);
      kalshiStarts.set(eventTicker, { startMs: null, error: error.message, expiresAt: now() + KALSHI_START_RETRY_MS });
      return null;
    }
  }

  // The start an order or a quote carries: the earlier of Unabated's and
  // Kalshi's, with both kept for the log.
  //   unabatedStart  the row's eventStart, ISO UTC
  function startFields(unabatedStart, kalshiStartMs) {
    const startMs = taker.earliestStartMs(Date.parse(unabatedStart), kalshiStartMs);
    return {
      eventStart: startMs === null ? unabatedStart : new Date(startMs).toISOString(),
      unabatedStart, kalshiStart: kalshiStartMs === null ? null : new Date(kalshiStartMs).toISOString(),
    };
  }

  // Kalshi's market object for a ticker (public GET /markets/<ticker>).
  async function readKalshiMarket(ticker) {
    const body = await getServiceJson(kalshiFetch, `${kalshiBaseUrl}/markets/${encodeURIComponent(ticker)}`);
    if (!body || typeof body.market !== "object" || body.market === null) throw new Error(`Kalshi market ${ticker}: no market object`);
    return body.market;
  }

  // Send every "bet" decision to the bets service's order desk, best live
  // edge first, one at a time (each changes the open exposure the next is
  // capped by), each after Kalshi's market is read and found to be the
  // row's (taker.marketMismatch). The outcome rides on the decision as
  // `order`. Mutates `decisions`.
  // Kalshi and Novig bets are one list (Cal, 2026-10-16: the better all-in
  // price first): a number both books trade at a tie sends the better live
  // edge first, and the other waits on the desk's unseen_fill. A venue whose
  // desk is off stops only that venue's sends.
  async function sendTakerOrders(decisions, rows, passStartedAt) {
    const rowByKey = new Map(rows.map((row) => [row.key, row]));
    const bets = decisions.filter((decision) => decision.reason === "bet").sort((a, b) => b.liveEdge - a.liveEdge);
    const venuesOff = new Set();
    for (const decision of bets) {
      if (venuesOff.has(decision.venue)) continue;
      const row = rowByKey.get(decision.rowKey);
      if (decision.venue === taker.VENUE_NOVIG) {
        if (!(await novigMarketConfirmed(decision, row))) continue;
        decision.order = await sendNovigTakerOrder(decision, passStartedAt);
      } else {
        if (!(await marketConfirmed(decision, row))) continue;
        decision.order = await sendTakerOrder(decision, passStartedAt);
      }
      if (decision.order.status === "off") venuesOff.add(decision.venue);
    }
  }

  function novigMarketUrl(marketId) {
    return `${novigBaseUrl}/v3/public/catalog/markets/${encodeURIComponent(marketId)}`;
  }

  // Why Novig's market and its event, read now, are not the row's or do not
  // take a free pregame order (taker.novigMarketMismatch), an unread market
  // included (fail closed); null when they are. Shared by the bot's orders
  // and the Place button's quotes, so neither can skip a check the other has.
  //   listedEventId  the event the catalog listed the market under, when known:
  //                  then the market and the event are read at once
  async function readNovigMarketMismatch(row, contract, listedEventId = null) {
    const eventUrl = (eventId) => `${novigBaseUrl}/v3/public/catalog/events/${encodeURIComponent(eventId)}`;
    let market;
    let event;
    try {
      if (typeof listedEventId === "string" && listedEventId) {
        [market, event] = await Promise.all([getNovigJson(novigMarketUrl(contract.marketId)), getNovigJson(eventUrl(listedEventId))]);
      } else {
        market = await getNovigJson(novigMarketUrl(contract.marketId));
        event = await getNovigJson(eventUrl(market.eventId));
      }
    } catch (error) {
      return `Novig market ${contract.marketId} unread: ${error.message}`;
    }
    const mismatch = taker.novigMarketMismatch(market, event, row, contract);
    return mismatch === null ? null : `Novig market ${contract.marketId}: ${mismatch}`;
  }

  // The bot's order goes only when readNovigMarketMismatch finds nothing; else the decision reads market_mismatch.
  async function novigMarketConfirmed(decision, row) {
    const mismatch = await readNovigMarketMismatch(row, { outcomeId: decision.ticker, marketId: decision.marketId });
    if (mismatch === null) return true;
    Object.assign(decision, { reason: "market_mismatch", detail: mismatch });
    return false;
  }

  // One POST /taker/novig_order.json; never throws. {status, ...} as the desk
  // answers it, "off" when the service has no bot Novig key, "error" on no answer.
  async function sendNovigTakerOrder(decision, passStartedAt) {
    const body = {
      clientId: taker.novigClientIdOf(decision.ticker, passStartedAt), outcomeId: decision.ticker, marketId: decision.marketId,
      price: decision.liveAsk, qty: decision.contracts, rowKey: decision.rowKey, eventId: decision.eventId, eventStart: decision.eventStart,
      heldIds: decision.heldIds, decidedAt: decision.decidedAt, novigBaseUrl, line: decision.line,
    };
    try {
      const reply = await sendServiceJson(serviceFetch, `${deps.betsServiceUrl}/taker/novig_order.json`, "POST", body, TAKER_ORDER_TIMEOUT_MS);
      novigDeskStatus = reply.status === "geo_blocked" ? `geo_blocked: ${reply.detail}`
        : reply.status === "stopped" ? `stopped: ${reply.detail}` : "armed";
      if (!["ceiling", "stopped", "unseen_fill", "stale_holdings"].includes(reply.status)) {
        logInfo(`taker novig order ${body.clientId} ${decision.ticker} x${decision.contracts} @ ${decision.liveAsk}: ${reply.status}${reply.detail ? ` (${reply.detail})` : ""}`);
      }
      return { clientId: body.clientId, ...reply };
    } catch (error) {
      if (error.status === 503) {
        novigDeskStatus = `off: the bets service places no Novig orders (${error.message})`;
        return { status: "off", detail: error.message };
      }
      logWarning(`taker novig order ${body.clientId} ${decision.ticker}: ${error.message}`);
      return { status: "error", detail: error.message };
    }
  }

  // Kalshi's market for the decision's ticker is active and prices the row
  // (its strike, its YES team), else the decision becomes market_mismatch
  // (fail closed, read failures included).
  async function marketConfirmed(decision, row) {
    let market;
    try {
      market = await readKalshiMarket(decision.ticker);
    } catch (error) {
      Object.assign(decision, { reason: "market_mismatch", detail: `${decision.ticker} unread: ${error.message}` });
      return false;
    }
    const mismatch = taker.marketMismatch(market, row, { side: decision.side, ticker: decision.ticker });
    if (mismatch === null) {
      // The order carries the earlier start; the desk refuses it once passed ("started").
      Object.assign(decision, startFields(decision.eventStart, await kalshiStartOf(taker.eventTickerOf(market, decision.ticker))));
      return true;
    }
    Object.assign(decision, { reason: "market_mismatch", detail: `${decision.ticker}: ${mismatch}` });
    return false;
  }

  // One POST /taker/order.json; never throws. {status, ...} as the desk
  // answers it, "off" when the service has no bot key, "error" on no answer.
  async function sendTakerOrder(decision, passStartedAt) {
    const body = {
      clientOrderId: taker.clientOrderIdOf(decision.ticker, decision.side, passStartedAt),
      ticker: decision.ticker, side: decision.side, ask: decision.liveAsk, allIn: decision.liveAllIn,
      count: decision.contracts, rowKey: decision.rowKey, eventId: decision.eventId, eventStart: decision.eventStart,
      heldIds: decision.heldIds, decidedAt: decision.decidedAt, kalshiBaseUrl,
    };
    try {
      const reply = await sendServiceJson(serviceFetch, `${deps.betsServiceUrl}/taker/order.json`, "POST", body, TAKER_ORDER_TIMEOUT_MS);
      takerDeskStatus = "armed";
      if (!["ceiling", "stopped", "unseen_fill", "started"].includes(reply.status)) logInfo(`taker order ${body.clientOrderId} ${decision.ticker} ${decision.side} x${decision.contracts} @ ${decision.liveAsk}: ${reply.status}`);
      return { clientOrderId: body.clientOrderId, ...reply };
    } catch (error) {
      if (error.status === 503) {
        takerDeskStatus = `off: the bets service places no orders (${error.message})`;
        return { status: "off", detail: error.message };
      }
      logWarning(`taker order ${body.clientOrderId} ${decision.ticker}: ${error.message}`);
      return { status: "error", detail: error.message };
    }
  }

  // What a decision row says about the line, for the log and /taker.json.
  function takerRowFields(row) {
    return {
      rowKey: row.key, leagueId: row.leagueId, league: row.leagueLabel, eventId: row.eventId, eventName: row.eventName,
      eventStart: row.eventStart, period: row.period, betType: row.betType, sideLabel: row.sideLabel, points: row.points, isAlt: row.isAlt,
    };
  }

  // Per league: rows looked at, rows with a ticker (derived moneylines and main numbers apart),
  // rows without one. T1's coverage check.
  function takerCoverage(decisions) {
    const byLeague = {};
    for (const decision of decisions) {
      const league = byLeague[decision.league] || (byLeague[decision.league] = { rows: 0, withTicker: 0, derivedTicker: 0, noTicker: 0 });
      league.rows += 1;
      if (decision.reason === "no_ticker") league.noTicker += 1;
      else if (decision.ticker) league.withTicker += 1;
      if (decision.derivedTicker) league.derivedTicker += 1;
    }
    return byLeague;
  }

  // Queue each decision whose signature changed; forget rows that left the
  // list, so one coming back logs again.
  function queueChangedDecisions(decisions) {
    const present = new Set();
    for (const decision of decisions) {
      present.add(decision.rowKey);
      const signature = taker.decisionSignature(decision);
      if (takerLogged.get(decision.rowKey) === signature) continue;
      takerLogged.set(decision.rowKey, signature);
      takerUnsent.push(decision);
    }
    for (const rowKey of Array.from(takerLogged.keys())) if (!present.has(rowKey)) takerLogged.delete(rowKey);
    if (takerUnsent.length > TAKER_LOG_MAX_UNSENT) {
      const dropped = takerUnsent.length - TAKER_LOG_MAX_UNSENT;
      takerUnsent = takerUnsent.slice(dropped);
      takerLogStatus.dropped += dropped;
      logWarning(`taker log: the bets service is behind, dropped the ${dropped} oldest unsent decision(s)`);
    }
  }

  // POST the queued decisions, oldest first; a failure keeps the rest for the next pass.
  async function postTakerDecisions() {
    if (takerPostBusy || takerUnsent.length === 0) return;
    takerPostBusy = true;
    try {
      while (takerUnsent.length) {
        const chunk = takerUnsent.slice(0, TAKER_LOG_ROWS_PER_POST);
        try {
          await sendServiceJson(serviceFetch, `${deps.betsServiceUrl}/taker_decisions.json`, "POST", { rows: chunk });
        } catch (error) {
          // A chunk the service refuses would be refused forever and hold every later row behind it.
          if (error.status !== HTTP_BAD_REQUEST) throw error;
          logWarning(`taker log: the bets service refused ${chunk.length} decision row(s), dropped: ${error.message}`);
          takerLogStatus.dropped += chunk.length;
        }
        takerUnsent = takerUnsent.slice(chunk.length);
        if (takerLogStatus.error) logInfo("taker log: saving again");
        Object.assign(takerLogStatus, { okAt: now(), error: null, rowsSent: takerLogStatus.rowsSent + chunk.length });
      }
    } catch (error) {
      if (takerLogStatus.error !== error.message) logWarning(`taker log: POST failed: ${error.message}`);
      takerLogStatus.error = error.message;
    } finally {
      takerPostBusy = false;
    }
  }

  // /taker.json: the rules, the last pass, its decisions (would-bets first,
  // then by live edge) and the log's state.
  function takerPayloadNow() {
    const order = (decision) => (decision.reason === "bet" ? 0 : 1);
    const decisions = takerDecisions.slice().sort((a, b) => order(a) - order(b) || (b.liveEdge ?? -Infinity) - (a.liveEdge ?? -Infinity));
    return {
      generatedAt: new Date(now()).toISOString(),
      desk: takerDeskStatus,
      novigDesk: novigDeskStatus,
      rules: {
        kellyBankroll: taker.KELLY_BANKROLL, novigKellyBankroll: taker.NOVIG_KELLY_BANKROLL, kelly: taker.BOT_KELLY, minLiveEdge: taker.MIN_LIVE_EDGE, maxLiveEdge: taker.MAX_LIVE_EDGE,
        priceTolerance: taker.PRICE_TOLERANCE, maxLineAgeHours: taker.MAX_LINE_AGE_MS / HOUR_MS,
      },
      pass: { ...takerLastPass, startedAt: isoOrNull(takerLastPass.startedAt), finishedAt: isoOrNull(takerLastPass.finishedAt) },
      decisions,
      log: { ...takerLogStatus, unsent: takerUnsent.length },
    };
  }

  // Bets and settings first, so the first scan uses the right leagues and the
  // first list is sized against what is held; then the polls, and resolve
  // once the first full load has landed.
  async function start() {
    await pollBets();
    await pollSettings();
    pollTimers = [
      timers.setInterval(() => pollBets(), BETS_POLL_MS),
      timers.setInterval(() => pollSettings(), SETTINGS_POLL_MS),
      timers.setInterval(() => postClosingFairs(), CLOSING_FAIRS_POST_MS),
      timers.setInterval(() => takerPass(), TAKER_PASS_MS),
    ];
    await scanLoad;
  }

  function stop() {
    for (const id of pollTimers) timers.clearInterval(id);
    pollTimers = [];
    scanner.stop();
    scannedLeagues = null;
  }

  function edgesPayloadNow() {
    return edgesPayload.buildEdgesPayload({
      feedState, scannerStatus, history, betRecords: held.records, fillFairIndex, knownStarts,
      stakeSettings: settings.stakeSettings, edgeSettings: settings.edgeSettings,
      settingsStatus: { source: settings.source, error: settings.error, okAt: settings.okAt, updatedAt: settings.updatedAt },
      betsStatus: { ...betsStatus }, boardLines: boardLines(), ladderReaderOf,
      teasers: openTeasersNow(), measurement: feedState ? tailFlexMeasurement() : null, now: now(),
      dismissedIds, edgeSkips,
    });
  }

  // Both football boards in, or failed (panel.js teaserBoardReady): before
  // that a list would be NFL's alone while CFB is still coming.
  function teaserBoardReady() {
    if (!scannerStatus || !feedState) return false;
    return teaser.TEASER_LEAGUE_IDS.every((id) => scannerStatus.leaguesLoaded.includes(id)
      || Boolean(scannerStatus.leagueErrors && scannerStatus.leagueErrors[id]));
  }

  // The Teasers tab's model, built as panel.js teaserModel builds it, with
  // the Can't tease marks whose game is still to start.
  function teaserModel(at) {
    if (!teaserBoardReady()) return null;
    const board = currentTeaserBoard();
    const blocked = new Set(teaserBlocks.filter((block) => block.eventStartMs > at).map((block) => block.marketKey));
    const legs = teaser.teaserLegs(feedState, { now: at, maxLineAgeMs: settings.edgeSettings.maxLineAgeHours * HOUR_MS, board });
    const placed = teaser.openTeasers(held.records, boardLines(), { now: at, ladderOf: board.ladderOf });
    const straights = teaser.heldStraights(held.records, boardLines(), { now: at, ladderOf: board.ladderOf });
    const kellyBankroll = settings.stakeSettings.bankroll * settings.stakeSettings.multiplier;
    const planned = teaser.planTeasers({ legs, placed, straights, kellyBankroll, previous: teaserBuild, blocked });
    teaserBuild = planned.build;
    return {
      legs, placed, build: teaserBuild,
      plan: teaser.describePlan(teaserBuild, legs, placed),
      legRows: teaser.describeLegs(legs, teaserBuild, placed, straights, blocked),
    };
  }

  function teasersPayloadNow() {
    const at = now();
    let model = null;
    let error = null;
    try {
      model = teaserModel(at);
      teaserLastError = null;
    } catch (thrown) {
      error = `Teasers failed: ${thrown.message}`;
      if (teaserLastError !== thrown.message) logWarning(`teasers failed: ${thrown.stack || thrown.message}`);
      teaserLastError = thrown.message;
    }
    return teasersPayload.buildTeasersPayload({ model, error, scannerStatus, stakeSettings: settings.stakeSettings, now: at });
  }

  // The Bets tab's Attach, off the board this runner holds (attach.js):
  // step 1 lists the games for `query`; with `eventId`, step 2 is the plan
  // for that game and the POST /pins.json body. Throws a 400/404-tagged
  // Error naming what was wrong.
  //   params  URLSearchParams {betId, query?, eventId?, swapped?}
  function attachPayloadNow(params) {
    const betId = params.get("betId");
    const bet = held.records.find((record) => record.id === betId && record.status === "open");
    if (!bet) throw httpError(404, `no open bet ${JSON.stringify(betId)}`);
    const lines = boardLines();
    const eventIdText = params.get("eventId");
    if (eventIdText == null) {
      const { scope, events, more } = attachLib.attachCandidates(bet, lines, { query: params.get("query") || "" });
      return {
        betId, scope, more,
        events: events.map(({ event, why }) => ({ eventId: event.eventId, label: attachLib.gameLabel(event), meta: attachLib.gameMeta(event), why: why || null })),
      };
    }
    const game = attachLib.boardGames(lines, null).find((candidate) => String(candidate.eventId) === eventIdText);
    if (!game) throw httpError(404, `game ${eventIdText} has left the board`);
    const swapped = params.get("swapped") === "1";
    const plan = attachLib.attachPlan(bet, game, { swapped, lines });
    return {
      betId, eventId: game.eventId, swapped, label: attachLib.gameLabel(game), meta: attachLib.gameMeta(game),
      venueLabel: betsLib.venueLabel(bet.venue), names: plan.names, betOnGame: plan.betOnGame,
      learnCount: plan.crosswalk.length, learnedNames: plan.crosswalk.map((row) => row.venueTeamName),
      pinRequest: attachLib.pinRequest(bet, game, plan),
    };
  }

  // The Live tab's cards, off the remembered games (scenarios.js).
  function scenariosPayloadNow() {
    const at = now();
    const games = Array.from(kickoffMemo.values(), (entry) => {
      const ladders = new Map();
      const ladderOf = entry.lines === null ? null : (period, axis) => {
        const periodTypeId = ladderLib.periodTypeIdOf(period);
        if (periodTypeId == null) return null;
        const key = `${periodTypeId}|${axis}`;
        if (!ladders.has(key)) ladders.set(key, ladderLib.buildLadder(entry.lines, { periodTypeId, axis }));
        return ladders.get(key);
      };
      return { row: entry.row, ladderOf, oddsAt: entry.oddsAt };
    });
    const records = held.records.filter((record) => !excludedBetIds.has(record.id));
    return {
      generatedAt: new Date(at).toISOString(),
      scanner: scannerStatus ? { phase: scannerStatus.phase, error: scannerStatus.error } : { phase: "starting", error: null },
      betsService: { okAt: betsStatus.okAt, error: betsStatus.error },
      ...scenarios.buildScenarios({ records, games, now: at }),
    };
  }

  function health() {
    return {
      ok: true, generatedAt: new Date(now()).toISOString(), uptimeSec: Math.round((now() - startedAt) / 1000),
      scanner: scannerStatus ? { phase: scannerStatus.phase, error: scannerStatus.error, lastSnapshotAt: scannerStatus.lastSnapshotAt } : { phase: "starting", error: null },
      betsService: { okAt: betsStatus.okAt, error: betsStatus.error },
      settings: { source: settings.source, okAt: settings.okAt, error: settings.error },
      closingFairs: { ...closingStatus },
      taker: {
        desk: takerDeskStatus, novigDesk: novigDeskStatus, lastPassAt: isoOrNull(takerLastPass.finishedAt),
        skipped: takerLastPass.skipped, novigSkipped: takerLastPass.novigSkipped,
        rows: takerLastPass.rows, byReason: takerLastPass.byReason, log: { okAt: takerLogStatus.okAt, error: takerLogStatus.error, dropped: takerLogStatus.dropped, unsent: takerUnsent.length },
      },
    };
  }

  return {
    start, stop, pollBets, pollSettings, postClosingFairs, scanLoaded: () => scanLoad, learn, takerPass, postTakerDecisions,
    edgesPayload: edgesPayloadNow, teasersPayload: teasersPayloadNow, attachPayload: attachPayloadNow,
    scenariosPayload: scenariosPayloadNow, takerPayload: takerPayloadNow, kalshiQuote: kalshiQuoteNow, novigQuote: novigQuoteNow, health, scanner,
    boardLines,
  };
}

// The line the Place button names, from the query: ids are whole numbers,
// points a number or absent (a moneyline). Throws a 400-tagged Error.
function kalshiLineSpecOf(params) {
  const whole = (name, allowed) => {
    const text = params.get(name);
    const value = Number(text);
    if (text === null || text === "" || !Number.isInteger(value) || (allowed && !allowed.includes(value))) {
      throw httpError(400, `${name} must be ${allowed ? `one of ${allowed.join(", ")}` : "a whole number"}, got ${JSON.stringify(text)}`);
    }
    return value;
  };
  const pointsText = params.get("points");
  const points = pointsText === null || pointsText === "" ? null : Number(pointsText);
  if (points !== null && !Number.isFinite(points)) throw httpError(400, `points must be a number, got ${JSON.stringify(pointsText)}`);
  return {
    eventId: whole("eventId"), periodTypeId: whole("periodTypeId"), betTypeId: whole("betTypeId", [1, 2, 3]),
    sideIndex: whole("sideIndex", [0, 1]), points,
  };
}

function isoOrNull(ms) {
  return typeof ms === "number" ? new Date(ms).toISOString() : null;
}

function countBy(items, keyOf) {
  const counts = {};
  for (const item of items) counts[keyOf(item)] = (counts[keyOf(item)] || 0) + 1;
  return counts;
}

function httpError(status, message) {
  const error = new Error(message);
  error.status = status;
  return error;
}

// body: a value to serialize, or a RawJson the board thread already serialized.
function sendJson(response, status, body) {
  const text = body instanceof boardThread.RawJson ? body.text : JSON.stringify(body);
  response.writeHead(status, { "Content-Type": "application/json", "Content-Length": Buffer.byteLength(text), "Cache-Control": "no-store" });
  response.end(text);
}

// The HTTP side: GET /edges.json, /teasers.json, /attach.json, /scenarios.json, /taker.json, /kalshi_quote.json, /novig_quote.json and /health, the Host allowlist first on
// every request. The port comes from the socket, not config, so a runner on
// an ephemeral or non-default port guards itself (as the bets service does).
function createHttpServer(runner, { host }) {
  const server = http.createServer(async (request, response) => {
    const allowed = allowedHosts(host, server.address().port);
    if (!hostAllowed(request.headers.host, allowed)) {
      sendJson(response, 403, { error: `Host must be one of ${JSON.stringify(allowed)}, got ${JSON.stringify(request.headers.host ?? null)}` });
      return;
    }
    if (request.method !== "GET") {
      sendJson(response, 405, { error: `only GET is served, got ${request.method}` });
      return;
    }
    const url = new URL(request.url, "http://runner.invalid");
    const path = url.pathname;
    try {
      // A page that just wrote a mark or a pin asks for the list on bets read now, not up to 30 s ago.
      const fresh = url.searchParams.get("fresh") === "1" && (path === "/edges.json" || path === "/teasers.json");
      // A runner on this thread reads here; the board client's fresh build reads on the board thread.
      if (fresh) await runner.pollBets();
      if (path === "/edges.json") return sendJson(response, 200, await runner.edgesPayload({ fresh }));
      if (path === "/teasers.json") return sendJson(response, 200, await runner.teasersPayload({ fresh }));
      if (path === "/attach.json") return sendJson(response, 200, await runner.attachPayload(url.searchParams));
      if (path === "/scenarios.json") return sendJson(response, 200, await runner.scenariosPayload());
      if (path === "/taker.json") return sendJson(response, 200, await runner.takerPayload());
      if (path === "/kalshi_quote.json") return sendJson(response, 200, await runner.kalshiQuote(url.searchParams));
      if (path === "/novig_quote.json") return sendJson(response, 200, await runner.novigQuote(url.searchParams));
      if (path === "/health") return sendJson(response, 200, await runner.health());
    } catch (error) {
      if (error.status) return sendJson(response, error.status, { error: error.message });
      logError(`${path} failed: ${error.stack || error.message}`);
      return sendJson(response, 500, { error: `${path} failed: ${error.message}` });
    }
    return sendJson(response, 404, { error: `no route for ${path}` });
  });
  return server;
}

// The main thread: the HTTP server and the board client; the runner itself
// (createRunner) runs on the board thread (board_thread.js), so no board work
// holds up a reply the main thread can give from a copy.
async function main() {
  const config = configFromEnv(process.env);
  const runner = boardThread.createBoardClient({ workerData: { role: "board", config }, logInfo: log, logWarning: logError });
  const server = createHttpServer(runner, config);
  await new Promise((resolve, reject) => {
    server.once("error", reject);
    server.listen(config.port, config.host, resolve);
  });
  log(`serving http://${config.host}:${config.port}/edges.json (bets service ${config.betsServiceUrl})`);
  const shutdown = (signal) => {
    log(`stopping (${signal})`);
    runner.stop();
    server.close(() => process.exit(0));
  };
  process.on("SIGINT", () => shutdown("SIGINT"));
  process.on("SIGTERM", () => shutdown("SIGTERM"));
  runner.start();
}

if (require.main === module) {
  main().catch((error) => {
    logError(`fatal: ${error.stack || error.message}`);
    process.exit(1);
  });
}

module.exports = {
  DEFAULT_HOST, DEFAULT_PORT, DEFAULT_BETS_SERVICE_URL, DEFAULT_KALSHI_BASE_URL, DEFAULT_NOVIG_BASE_URL, BETS_POLL_MS, SETTINGS_POLL_MS, CLOSING_FAIRS_POST_MS, TAKER_PASS_MS,
  configFromEnv, allowedHosts, hostAllowed, scannerLeaguesOf, nodeFetch, createRunner, createHttpServer, logError,
};
