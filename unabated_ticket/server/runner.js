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
// rows, plus CBB, whose fair ladders it stores with NFL's and CFB's
// (server/laddersnapshot.js), plus the leagues of open straight bets, whose
// closing fairs it saves — each league re-downloaded on its own size-based cadence, a full
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
//   POST <bets service>/fill_snapshots.json — each new open bet's market at
//       the fill (server/edgesnapshot.js): every book's two sides at its
//       number, Unabated's main number and fair, the edge-move tag. Runner only.
//   POST <bets service>/close_snapshots.json every 30 s with the closing
//       fairs — every book's two sides at each open bet's number before its
//       game starts; the service keeps the newest pre-start one.
//   POST <bets service>/ladders.json every 30 s — Unabated's full-game fair
//       ladder (spread, total, moneyline) of every NFL / CFB / CBB game still
//       to start (server/laddersnapshot.js): an opener once per game and
//       market, then the close whenever the ladder changes or every 10 min
//       unchanged; ladders.duckdb keeps the first opener and the newest
//       pregame close (issue #264). Stored only, nothing reads it back.
//   POST <bets service>/edge_segments.json every 30 s — every NFL / CFB / CBB
//       full-game line showing a 1%+ edge, as segments at one number, price,
//       fair and start from first sighting to the edge's end or kickoff
//       (server/edgelog.js), into ladders.duckdb: whether an edge seen early
//       was still there later. Stored only, nothing reads it back.
//   GET <bets service>/tau_samples.json once on start and POST /tau_samples.json
//       every 30 s with the readings that changed (an unchanged one hourly;
//       each pass absorbs the board's current refresh first): the tail
//       flex's tau pool (tailflex.createTauPool, issue #208), the day's
//       two-sided exchange rungs, so a restart keeps the measured tau.
//       GET /tau.json here serves the pooled tau per market (the panel ranks
//       and gates on it through the bets service's pass-through).
//   Every league's snapshot every 60 s (scanner maxRefreshMs: the taker
//       trades against the fair, so the panel's 2-5 min tiers are too slow).
//   GET <Kalshi>/markets/<ticker>/orderbook every 30 s for each Kalshi Edges
//       row the automated taker looks at, and GET <Kalshi>/markets/<ticker>
//       before every order, so the market is active and its strike and team
//       are the row's (public, no key; server/taker.js marketMismatch).
//   POST <bets service>/taker/order.json — each "bet" decision, with the ids
//       of the held bets its stake was sized on: the bets service's order
//       desk places it on Kalshi on the bot's own key (MONEY; stop file, open-cost
//       ceiling and the unseen-fill rule there; bets_service/taker_orders.py).
//   POST <bets service>/novig_read.json — Novig's catalog, signed with the
//       bot's Novig key on Novig's keyed routes (bets_service/
//       novig_catalog_reads.py): /v3/catalog/markets?league=..&eventStatus=
//       OPEN_PREGAME every 30 s while Novig Edges rows are up, and
//       /v3/catalog/markets/<id> and /events/<id> before every Novig order and
//       Place quote (taker.novigMarketMismatch). The unkeyed /v3/public routes
//       share one small per-IP throttle, and Place quotes 429'd on them
//       (2026-10-28, 10-30); the key's own read bucket is 64 at once, 16 a second.
//   POST <bets service>/novig_books.json — the books of those rows' Novig
//       markets (and of a Place quote's), off the bets service's websocket on
//       the bot's Novig key (bets_service/novig_book_stream.py). Polling the
//       public /book per market shared Novig's per-IP public throttle with
//       Place and got Place HTTP 429 (2026-10-27). A Place quote reads the
//       keyed /book (novig_read.json) only when the stream has no book for it.
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
//   The live price pass (livePricePass), right after each taker pass: for every
//       Edges line at Kalshi, Novig and Polymarket US not read in the last
//       LIVE_REREAD_MS (the taker's own reads count), GET <Kalshi>/markets/<ticker>
//       once per ticker per 10 min and /markets/<ticker>/orderbook (at most
//       LIVE_KALSHI_READS_PER_PASS a pass, two at a time), POST <bets
//       service>/novig_books.json with background: true (subscribed behind
//       the bot's and Place's markets, no wait) and, for a market the stream
//       has no book of yet, POST <bets service>/novig_read.json for its keyed
//       /v3/catalog/markets/<id>/book (at most LIVE_NOVIG_REST_READS_PER_PASS
//       a pass, one at a time, stopping at the first failure), Polymarket US's GET
//       /v1/events (reused 5 min per game, at most
//       LIVE_POLYMARKET_EVENT_READS_PER_PASS fresh a pass) and POST <bets
//       service>/polymarket_books.json with background: true (the books off
//       its websocket, never the gateway's /book; at most
//       LIVE_POLYMARKET_READS_PER_PASS rows a pass), and POST <bets service>/prophetx_live.json
//       naming every game with a ProphetX Edges line (answered at once from the
//       service's own reads on its PROPHETX_LIVE_* key pair). Reads only; served at
//       /live_prices.json and laid over /edges.json (extension/liveprices.js).
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
//                    chance of each (server/scenarios.js documents the shape);
//                    ?not_origin=&not_venue=&not_league=&not_kind= (each may
//                    repeat) rebuilds the cards without the bets the Live
//                    tab's filter hides
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
//   GET /polymarket_quote.json?eventId=..&periodTypeId=..&betTypeId=..&sideIndex=..&points=..
//                    the same for Polymarket US Place (bets_service/polymarket_place.py):
//                    the board's Polymarket US line, the market slug and contract
//                    side found in Polymarket US's own event list (polymarket.js),
//                    its live ask, the all-in price with the taker fee and the
//                    line as a bet record names it; places nothing.
//   GET /prophetx_quote.json?eventId=..&periodTypeId=..&betTypeId=..&sideIndex=..&points=..
//                    the same for ProphetX Place (bets_service/prophetx_place.py):
//                    the board's ProphetX line, the line as a bet record names it
//                    and the ProphetX strike ids the feed carries for its game
//                    (prophetx.js); ProphetX itself is read by the bets service
//                    on its key. Places nothing.
//   GET /bfa_quote.json?eventId=..&periodTypeId=..&betTypeId=..&sideIndex=..&points=..
//                    the same for BFA Place (bets_service/bfa_place.py): the
//                    board's one Buckeye (book 59) line, its side's rotation
//                    (BFA's rotations are Unabated's own), whether it is an
//                    alt, and the line as a bet record names it. BFA's board
//                    is read by the bets service on its login. Places nothing.
//   GET /live_prices.json  the live price overlay's reads (extension/liveprices.js
//                    createLineView's input): {generatedAt, maxAgeMs, bookIds (the
//                    board's Kalshi, Novig, Polymarket US and, once the bets service
//                    reads it, ProphetX book ids), reads:
//                    {lineKey: {venue, points, allIn, depthDollars, readAtMs} |
//                    {venue, points, empty: true, readAtMs}}, lastPass}. /edges.json
//                    already lists through it; the panel reads it via the bets service.
//   GET /health      {ok, generatedAt, uptimeSec, scanner, betsService, settings, closingFairs, closeSnapshots, tauPool, livePrices, taker}
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
const edgesnapshot = require("./edgesnapshot.js");
const laddersnapshot = require("./laddersnapshot.js");
const edgelog = require("./edgelog.js");
const ladderLib = require("../extension/ladder.js");
const scenarios = require("./scenarios.js");
const teasersPayload = require("./teasers_payload.js");
const attachLib = require("../extension/attach.js");
const pricecheck = require("../extension/pricecheck.js");
const taker = require("./taker.js");
const polymarket = require("./polymarket.js");
const prophetx = require("./prophetx.js");
const liveboard = require("./liveboard.js");
const feedLib = require("../extension/feed.js");
const boardThread = require("./board_thread.js");
const livePrices = require("../extension/liveprices.js");

const DEFAULT_HOST = "127.0.0.1";
const DEFAULT_PORT = 8095;
// A moneyline has no number: the Place button's line matches on side alone.
const MONEYLINE_BET_TYPE_ID = 1;
const SPREAD_BET_TYPE_ID = 2;
// Unabated's book id for Buckeye, which is BFA (teaser.js BUCKEYE_BOOK_ID).
const BUCKEYE_BOOK_ID = 59;
const DEFAULT_BETS_SERVICE_URL = "http://127.0.0.1:8094";
const DEFAULT_KALSHI_BASE_URL = "https://api.elections.kalshi.com/trade-api/v2";
const DEFAULT_NOVIG_BASE_URL = "https://api.novig.com";
// Polymarket US's public gateway (no key): events, markets and order books.
const DEFAULT_POLYMARKET_GATEWAY_URL = "https://gateway.polymarket.us";
// Novig's catalog: one page holds up to 5000 markets; a slate is a few hundred.
const NOVIG_MARKETS_PAGE_LIMIT = 5000;
const NOVIG_CATALOG_MARKETS_PATH = "/v3/catalog/markets";
// One keyed catalog read through the bets service: Novig's own 10 s timeout there, plus the hop.
const NOVIG_READ_TIMEOUT_MS = 12 * 1000;
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
// How long a Novig book ask may wait for the stream's first snapshot of a game it was not watching.
const NOVIG_STREAM_WAIT_MS = 2000;
// The live price pass (livePricePass, run right after each taker pass): how long a read
// (its own or the taker's) stands before the line is read again, and how
// long a Kalshi market's check (taker.marketMismatch) and a Polymarket US
// event list are reused (neither changes inside it; the book read every pass
// still shows a closed market as empty).
const LIVE_REREAD_MS = 20 * 1000;
const KALSHI_MARKET_CHECK_TTL_MS = 10 * 60 * 1000;
const POLYMARKET_EVENTS_REUSE_MS = 5 * 60 * 1000;
// The most Novig markets one POST /novig_books.json may name (service.py MAX_NOVIG_BOOK_MARKETS).
const NOVIG_BOOK_MARKETS_PER_ASK = 500;
// Rows whose Polymarket US contract the live pass looks up at once (the
// event-list reads behind them are the only gateway reads it sends).
const POLYMARKET_CONTRACT_LOOKUPS_AT_ONCE = 2;
// Fresh Polymarket US event lists one live pass may read off the public
// gateway (the rest wait a pass; a list is reused POLYMARKET_EVENTS_REUSE_MS).
// The gateway allows about 5 reads then a 10 s wait per address (2026-11-05),
// and Place's own event-list read shares it, so the pass keeps well under.
const LIVE_POLYMARKET_EVENT_READS_PER_PASS = 3;
// The most slugs one POST /polymarket_books.json may name (service.py MAX_POLYMARKET_BOOK_SLUGS).
const POLYMARKET_BOOK_SLUGS_PER_ASK = 500;
// How long Place's quote waits for a market the stream had not watched yet
// (polymarket_book_stream.MAX_WAIT_SEC caps it), and the service's answer on top.
const POLYMARKET_STREAM_WAIT_MS = 2000;
const POLYMARKET_STREAM_READ_TIMEOUT_MS = POLYMARKET_STREAM_WAIT_MS + 3000;
// Per-market reads a live price pass sends at most, best-need first (the rest
// wait a pass): Kalshi's and Polymarket US's read budgets are shared with the
// taker and Place, and a read stands for LIVE_READ_MAX_AGE_MS (2 min, four
// passes), so up to ~4x this many lines stay read at each venue.
const LIVE_KALSHI_READS_PER_PASS = 120;
const LIVE_POLYMARKET_READS_PER_PASS = 60;
// Novig's stream adds only ~15 markets a minute (16 tokens a subscribe,
// refilling 4 a second), so on a big slate, and after every restart, most
// Novig lines would wait many minutes for a stream book (2026-11-06: 3 of 218
// read 7 min after a deploy). A market without one is read off the keyed REST
// /book instead, at most this many a pass (~2 reads/s at a pass every 30 s,
// of the key's 16/s read budget the bot and Place share), oldest read first.
const LIVE_NOVIG_REST_READS_PER_PASS = 60;
// A live all-in price is kept to 6 decimals (Polymarket US's fee math rounds there).
const LIVE_PRICE_DECIMALS = 1e6;
// The venue name the overlay's ProphetX reads carry, and the most games one
// POST /prophetx_live.json names (prophetx_live.py MAX_GAMES_PER_ASK).
const VENUE_PROPHETX = "prophetx";
const PROPHETX_LIVE_GAMES_PER_ASK = 300;
// Time for the bets service to answer, on top of that wait.
const NOVIG_STREAM_ANSWER_MS = 3000;
const NOVIG_STREAM_READ_TIMEOUT_MS = NOVIG_STREAM_WAIT_MS + NOVIG_STREAM_ANSWER_MS;
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
// Tau samples (tailflex createTauPool, issue #208) are POSTed in chunks of
// this: one multi-row UPSERT under the bets service's lock, kept short so
// /bets.json never waits behind it (the service refuses over
// MAX_TAU_SAMPLE_ROWS_PER_POST per request).
const TAU_SAMPLES_PER_POST = 500;
// How often the unsent tau samples are POSTed and the pool's GET retried.
const TAU_SAMPLES_POST_MS = 30 * 1000;
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
// Close snapshots per POST: each carries every book's two sides (a few KB),
// so a full slate of open bets stays far under the service's 1 MB body cap.
const CLOSE_SNAPSHOT_ROWS_PER_POST = 25;
// Ladder rows per POST /ladders.json: a CBB spread ladder is a few KB, so 40
// stay far under the service's 1 MB body cap.
const LADDER_ROWS_PER_POST = 40;
// Edge segment rows are ~400 bytes: 500 to a POST is ~200 KB, under the
// service's 1 MB body cap and its 1000-row limit.
const EDGE_SEGMENT_ROWS_PER_POST = 500;
// Edge segment rows held while the service is down; past this the oldest go
// (counted in /health), so an outage cannot grow the runner without bound.
const EDGE_SEGMENT_MAX_UNSENT = 50000;
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
// open bets' leagues, so each has a closing fair, plus the leagues whose fair
// ladders are stored (laddersnapshot.js: NFL, CFB, CBB). A league outside the
// settings never reaches the Edges list (edgerows.js leagueIds) or the taker
// (taker.js leagueIds), and tail flex measures per league, so loading CBB for
// its ladders changes no other league's rows or stakes.
function scannerLeaguesOf(edgeLeagues, betLeagueIds = []) {
  return Array.from(new Set([...edgeLeagues, ...teaser.TEASER_LEAGUE_IDS, ...laddersnapshot.LADDER_LEAGUE_IDS, ...betLeagueIds]))
    .sort((a, b) => a - b);
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
//          novigRead (path, query) -> {status, body} (Novig's keyed catalog,
//          default POST <bets service>/novig_read.json), novigBaseUrl,
//          polymarketFetch (Polymarket US's public gateway, default serviceFetch), polymarketGatewayUrl,
//          betsServiceUrl, now, timers {setInterval, clearInterval},
//          logInfo / logWarning (message) -> void, default stdout / stderr,
//          edgeSegmentMaxUnsent (default EDGE_SEGMENT_MAX_UNSENT; tests)}
// Returns {start, stop, pollBets, pollSettings, postLadders, postEdgeSegments, scanLoaded, learn, takerPass,
//          postTakerDecisions, edgesPayload, teasersPayload, attachPayload,
//          scenariosPayload, takerPayload, health, scanner, boardLines (the
//          board's rows, read by the tests)}.
function createRunner(deps) {
  const now = deps.now || (() => Date.now());
  const timers = deps.timers || { setInterval: (fn, ms) => setInterval(fn, ms), clearInterval: (id) => clearInterval(id) };
  const serviceFetch = deps.serviceFetch || deps.fetchImpl;
  const kalshiFetch = deps.kalshiFetch || serviceFetch;
  const kalshiBaseUrl = deps.kalshiBaseUrl || DEFAULT_KALSHI_BASE_URL;
  const novigRead = deps.novigRead || ((path, query) => sendServiceJson(serviceFetch, `${deps.betsServiceUrl}/novig_read.json`, "POST",
    { path, query }, NOVIG_READ_TIMEOUT_MS));
  const novigBaseUrl = deps.novigBaseUrl || DEFAULT_NOVIG_BASE_URL;
  const polymarketFetch = deps.polymarketFetch || serviceFetch;
  const polymarketGatewayUrl = deps.polymarketGatewayUrl || DEFAULT_POLYMARKET_GATEWAY_URL;
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
  // leagueId -> {version, byEvent}: each league's lines grouped by game as of
  // its last snapshot (the kickoff odds and the edge snapshots read it).
  const linesByEventOfLeague = new Map();
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
  // eventId -> {row, lines, oddsAt, gameState}: each game's board row and its
  // lines as of the last snapshot before it started (lines null when the
  // runner never saw it pregame), and its score and clock as of the newest
  // snapshot that carried it. The Live tab's chances are these kickoff odds.
  const kickoffMemo = new Map();
  // eventId -> the game's venue ids as its last pregame snapshot carried them
  // (liveboard.rememberPregameVenueIds), for a live Place quote: Unabated's
  // live lines carry none. Kept and pruned with kickoffMemo.
  const pregameVenueMemo = new Map();
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
  // One insert-only per-bet capture (fill books, fill snapshots): ids decided
  // this run, rows not yet taken by the service, the retry clock, the last error.
  const newCapture = () => ({ decided: new Set(), unsent: new Map(), retryAt: 0, lastError: null });
  const fillBooksCapture = newCapture();
  const fillSnapshotsCapture = newCapture();
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
  // betId -> what was last POSTed to /close_snapshots.json (edgesnapshot.unsentCloseRows).
  const sentCloseSnapshots = new Map();
  const closeSnapshotStatus = { okAt: null, error: null, rowsSent: 0 };
  // The stored fair ladders (laddersnapshot.js): leagueId -> {version,
  // ladders} as of the league's last snapshot, ladderKey -> what the service
  // took (laddersnapshot.markSent), and the last pass's outcome.
  const ladderCache = new Map();
  const sentLadders = new Map();
  let laddersBusy = false;
  const ladderStatus = { okAt: null, error: null, rowsSent: 0 };
  // Every edge's life to kickoff (edgelog.js): leagueId -> {version,
  // edgeLines} as of the league's last snapshot, lineKey -> its open segment,
  // segment key -> the newest row the service has not taken yet.
  const edgeLineCache = new Map();
  const edgeSegments = new Map();
  const unsentEdgeRows = new Map();
  let edgeLogBusy = false;
  const edgeLogStatus = { okAt: null, error: null, rowsSent: 0, dropped: 0 };
  const edgeSegmentMaxUnsent = deps.edgeSegmentMaxUnsent ?? EDGE_SEGMENT_MAX_UNSENT;
  // The tail flex's rolling tau pool (tailflex.createTauPool): every refresh's
  // two-sided exchange rungs where two or more exchanges quote both sides,
  // kept 24 h, loaded from the bets service on start (GET /tau_samples.json)
  // and persisted there as readings change (POST /tau_samples.json), so a
  // restart starts with the day's rungs, not an empty pool and the fallback.
  const tauPool = tailflex.createTauPool();
  // Unsent samples by key, the newest reading of each (the service keeps the
  // newest seen_at anyway), a reading older than the pool's window dropped
  // before each POST, so a long service outage holds at most one row per
  // pooled sample, never a row per refresh and never more than the pool.
  const tauUnsent = new Map();
  let tauBusy = false;
  const tauStatus = { loaded: false, loadError: null, loadSkipped: 0, okAt: null, error: null, rowsSent: 0 };
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
  // The same for Novig's in-game markets (a live Place reads it, nothing else).
  let novigLiveCatalog = null;
  // Unabated eventId -> its Novig event as the last pregame catalog read
  // showed it (taker.novigEventMemoryOf), for a live Place: in game the
  // board's lines carry no Novig id. Kept until kickoffMemo lets the game go.
  const novigEventMemo = new Map();
  // The live price overlay (extension/liveprices.js): lineKey -> the last read
  // of that line's exchange book, {venue, points, allIn, depthDollars,
  // readAtMs} or {venue, points, empty: true, readAtMs}; ticker -> Kalshi's
  // market check {mismatch, checkedAt}; Polymarket US event-list query ->
  // {events, readAt} | {error, readAt}; Polymarket US line (row key, points,
  // start) -> its contract {slug, side, startMs}, found once and kept until
  // the game starts, so the event-list budget goes only to lines not yet
  // found; and the last pass's counts.
  const liveReads = new Map();
  const kalshiMarketChecks = new Map();
  const polymarketEventLists = new Map();
  const polymarketLineContracts = new Map();
  let liveBusy = false;
  // Whether the bets service reads ProphetX's live prices (false while it answers 503: no live key pair).
  let prophetxLiveOn = false;
  let liveLastPass = { startedAt: null, finishedAt: null, rows: 0, byVenue: {}, error: null };
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
    fetchImpl: deps.fetchImpl, now, timers, maxRefreshMs: TAKER_SNAPSHOT_REFRESH_MS, keepLive: true,
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
      const linesByEvent = leagueLinesByEvent(leagueId);
      for (const row of leagueBoardRows(leagueId)) {
        if (typeof row.eventStartMs !== "number") continue;
        const pregame = row.eventStartMs > at;
        const known = kickoffMemo.get(row.eventId);
        const gameState = gameStateWithAbbreviations(row);
        if (pregame) {
          const lines = linesByEvent.get(row.eventId) || [];
          kickoffMemo.set(row.eventId, { row, lines, oddsAt: at, gameState });
          const venueIds = liveboard.rememberPregameVenueIds(feedState.events[row.eventId], lines);
          if (venueIds) pregameVenueMemo.set(row.eventId, venueIds);
        } else if (!known) kickoffMemo.set(row.eventId, { row, lines: null, oddsAt: null, gameState });
        else known.gameState = gameState;
      }
    }
    for (const leagueId of kickoffLeagueVersions.keys()) if (!leagueIds.includes(leagueId)) kickoffLeagueVersions.delete(leagueId);
    for (const [eventId, entry] of kickoffMemo) {
      if (at - entry.row.eventStartMs > KICKOFF_MEMO_HOURS * HOUR_MS) kickoffMemo.delete(eventId);
    }
    for (const eventId of pregameVenueMemo.keys()) if (!kickoffMemo.has(eventId)) pregameVenueMemo.delete(eventId);
  }

  // The game's latest score and clock (feed.gameStateOf) plus its teams' short
  // names for the Live card's score line, or null when the feed has no event.
  function gameStateWithAbbreviations(row) {
    const event = feedState && feedState.events[row.eventId];
    if (!event || !event.gameState) return null;
    const abbreviationOf = (teamId) => {
      const team = teamId == null ? null : feedState.teamIndex[teamId];
      return team ? team.abbreviation ?? null : null;
    };
    return { ...event.gameState, awayAbbr: abbreviationOf(row.awayTeamId), homeAbbr: abbreviationOf(row.homeTeamId) };
  }

  // One league's lines by game, regrouped only when its snapshot landed again.
  function leagueLinesByEvent(leagueId) {
    const version = scanner.getLeagueVersion(leagueId);
    const cached = linesByEventOfLeague.get(leagueId);
    if (cached && cached.version === version) return cached.byEvent;
    const byEvent = ladderLib.groupLinesByEvent(scanner.getLeagueLines(leagueId));
    linesByEventOfLeague.set(leagueId, { version, byEvent });
    for (const known of linesByEventOfLeague.keys()) if (!boardLeagueIds().includes(known)) linesByEventOfLeague.delete(known);
    return byEvent;
  }

  // One game's lines as of its league's last snapshot (edgesnapshot.js linesOfEvent).
  function linesOfEventNow(eventId) {
    const event = feedState && feedState.events[eventId];
    return event ? leagueLinesByEvent(event.leagueId).get(eventId) || [] : [];
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

  // Measured once per scanner update, on the same line-age gate as the list
  // (panel.js tailFlexMeasurement): c from this refresh, tau from the pool
  // after this refresh's samples join it (the changed ones queued for the
  // bets service).
  function tailFlexMeasurement() {
    const maxLineAgeMs = settings.edgeSettings.maxLineAgeHours * HOUR_MS;
    if (!tailFlexCache || tailFlexCache.maxLineAgeMs !== maxLineAgeMs) {
      const at = now();
      const measured = tailflex.measureTailFlex(feedState, { now: at, maxLineAgeMs });
      for (const sample of tauPool.absorb(measured.tauSamples, at)) tauUnsent.set(sample.key, sample);
      tailFlexCache = { maxLineAgeMs, measurement: tailflex.withPooledTau(measured, tauPool.measure(at)) };
    }
    return tailFlexCache.measurement;
  }

  // The row's tail-flex shrink for the bot's Min shrunk edge floor (taker.decideRow sizing.shrink).
  function shrinkOfRow(row) {
    const shrunk = tailflex.shrinkOfRow(row, tailFlexMeasurement());
    return shrunk ? shrunk.shrink : null;
  }

  // One GET of /tau_samples.json into the pool, once: the day's rungs as the
  // service kept them. Until it succeeds the pool holds only what this
  // process has seen, and /health says so.
  async function loadTauPool() {
    if (tauStatus.loaded) return;
    try {
      const body = await getServiceJson(serviceFetch, `${deps.betsServiceUrl}/tau_samples.json`);
      if (!body || !Array.isArray(body.rows)) throw new Error("tau_samples.json has no rows array");
      const { kept, skipped } = tauPool.load(body.rows);
      tailFlexCache = null;
      Object.assign(tauStatus, { loaded: true, loadError: null, loadSkipped: skipped });
      logInfo(`tau pool: loaded ${kept} sample(s) from the bets service`);
      if (skipped) logWarning(`tau pool: skipped ${skipped} persisted row(s) of another shape than a tau sample (a bug upstream; they are not pooled)`);
    } catch (error) {
      if (tauStatus.loadError !== error.message) logWarning(`tau pool load failed (retrying): ${error.message}`);
      tauStatus.loadError = error.message;
    }
  }

  // Absorb the board's current refresh into the pool (the measurement is
  // cached per refresh, so this costs nothing when /edges.json or the taker
  // already built it), then POST the queued tau samples, oldest first; a
  // failure keeps the rest for the next pass. A chunk the service refuses
  // (400) is dropped, as the taker log does. A reading that changed while
  // its chunk was in flight stays queued; one older than the pool's window
  // is dropped unsent (the service would prune it in the same write).
  async function postTauSamples() {
    if (tauBusy) return;
    tauBusy = true;
    try {
      await loadTauPool();
      if (feedState) tailFlexMeasurement();
      const at = now();
      for (const [key, sample] of tauUnsent) {
        if (at - sample.seenAt > tailflex.TAU_POOL_WINDOW_MS) tauUnsent.delete(key);
      }
      while (tauUnsent.size) {
        const sent = Array.from(tauUnsent.entries()).slice(0, TAU_SAMPLES_PER_POST);
        const chunk = sent.map(([, sample]) => ({ ...sample, seenAt: new Date(sample.seenAt).toISOString() }));
        try {
          await sendServiceJson(serviceFetch, `${deps.betsServiceUrl}/tau_samples.json`, "POST", { rows: chunk });
        } catch (error) {
          if (error.status !== HTTP_BAD_REQUEST) throw error;
          logWarning(`tau samples: the bets service refused ${chunk.length} row(s), dropped: ${error.message}`);
        }
        for (const [key, sample] of sent) {
          if (tauUnsent.get(key) === sample) tauUnsent.delete(key);
        }
        if (tauStatus.error) logInfo("tau samples: saving again");
        Object.assign(tauStatus, { okAt: now(), error: null, rowsSent: tauStatus.rowsSent + chunk.length });
      }
    } catch (error) {
      if (tauStatus.error !== error.message) logWarning(`tau samples POST failed: ${error.message}`);
      tauStatus.error = error.message;
    } finally {
      tauBusy = false;
    }
  }

  // /tau.json: the pooled tau per market, what the panel ranks and gates on
  // through the bets service's pass-through (the pool's state is on /health).
  function tauPayloadNow() {
    const at = now();
    return { generatedAt: new Date(at).toISOString(), markets: tauPool.measure(at) };
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
      .then(() => captureFillBooksAndSnapshots())
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
      await postCrosswalkRows(learned);
      crosswalkLastError = null;
    } catch (error) {
      crosswalkRetryAt = now() + SERVICE_WRITE_RETRY_MS;
      if (crosswalkLastError !== error.message) logWarning(`crosswalk: service write failed: ${error.message}`);
      crosswalkLastError = error.message;
    }
  }

  // POST crosswalk rows (insert-only at the service: a held key is never
  // rewritten) and re-key the held bets through the whole table it replies
  // with. Throws on a failed write.
  async function postCrosswalkRows(rows) {
    const reply = await sendServiceJson(serviceFetch, `${deps.betsServiceUrl}/crosswalk.json`, "POST", { rows });
    if (!reply || !Array.isArray(reply.crosswalk)) throw new Error("crosswalk.json reply has no crosswalk array");
    logInfo(`crosswalk: learned ${reply.learned} row(s), ${(reply.conflicts || []).length} refused by the service, ${reply.crosswalk.length} held`);
    held = { ...held, crosswalk: reply.crosswalk, records: betsLib.rekeyRecords(betsLib.applyPins(held.records, held.pins), reply.crosswalk) };
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
  // (fillbooks.js) and its market at the fill (edgesnapshot.js: every book's
  // two sides at its number, Unabated's main number and fair, the edge-move
  // tag), one row per POST, the same way as the fill fairs. Both read one
  // lookup of the board for the bets either still has to decide. Decided ids
  // live in memory only: after a restart every open bet was placed before
  // the scanner watched, so nothing is re-captured, and the service keeps the
  // first row regardless.
  async function captureFillBooksAndSnapshots() {
    const liveStatus = scanner.getStatus();
    const shared = {
      records: held.records, state: feedState, boardLines: boardLines(), history,
      observingSince: liveStatus.observingSince, leagueObservingSince: liveStatus.leagueObservingSince, now: now(),
    };
    const undecided = held.records.filter((bet) => bet.status === "open" && !bet.unmatchable
      && !(fillBooksCapture.decided.has(bet.id) && fillSnapshotsCapture.decided.has(bet.id)));
    if (typeof liveStatus.observingSince !== "number" || !undecided.length) return;
    const rowsByBet = fillfair.sameLineRows(undecided, feedState, shared.boardLines);
    const books = fillbooks.captureFillBooks({ ...shared, skipIds: fillBooksCapture.decided, rowsByBet });
    const snapshots = edgesnapshot.captureFillSnapshots({
      ...shared, skipIds: fillSnapshotsCapture.decided, rowsByBet, linesOfEvent: linesOfEventNow,
    });
    await sendCaptured(fillBooksCapture, books, "fill books", "fill_books.json");
    await sendCaptured(fillSnapshotsCapture, snapshots, "fill snapshots", "fill_snapshots.json");
  }

  // Record a capture pass's saves and refusals, then POST the unsent rows one
  // per request (the service refuses a whole request on its first bad row): a
  // 400 drops the row, any other failure keeps every unsent row for a minute.
  async function sendCaptured(capture, { saves, refusals }, label, path) {
    for (const save of saves) {
      capture.decided.add(save.betId);
      capture.unsent.set(save.betId, save);
    }
    for (const refusal of refusals) capture.decided.add(refusal.betId);
    if (saves.length) logInfo(`${label}: ${saves.length} captured`);
    if (capture.unsent.size === 0 || now() < capture.retryAt) return;
    for (const row of Array.from(capture.unsent.values())) {
      try {
        const reply = await sendServiceJson(serviceFetch, `${deps.betsServiceUrl}/${path}`, "POST", { rows: [row] });
        if (!reply || reply.ok !== true) throw new Error(`${path} reply is not ok`);
        capture.unsent.delete(row.betId);
        capture.lastError = null;
      } catch (error) {
        if (error.status === HTTP_BAD_REQUEST) {
          capture.unsent.delete(row.betId);
          logWarning(`${label}: the service refused ${row.betId}, dropped: ${error.message}`);
          continue;
        }
        capture.retryAt = now() + SERVICE_WRITE_RETRY_MS;
        if (capture.lastError !== error.message) logWarning(`${label}: service write failed: ${error.message}`);
        capture.lastError = error.message;
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

  // POST the closing-fair rows that changed since the last successful POST,
  // then the close snapshots of the same bets off the same lookup. A failure
  // keeps them unsent, so the next pass retries them.
  async function postClosingFairs() {
    if (closingBusy || !feedState || !scannerStatus) return;
    closingBusy = true;
    try {
      const pass = closefair.closingPass({
        records: held.records, state: feedState, boardLines: boardLines(),
        leagueLoadedAt: scannerStatus.leagueLoadedAt || {}, staleLeagues: scannerStatus.staleLeagues, now: now(),
      });
      await postClosingFairRows(pass.rows);
      await postCloseSnapshots(pass);
    } catch (error) {
      if (closingStatus.error !== error.message) logWarning(`closing pass failed: ${error.stack || error.message}`);
      closingStatus.error = error.message;
    } finally {
      closingBusy = false;
    }
  }

  async function postClosingFairRows(passRows) {
    try {
      const rows = closefair.unsentRows(passRows, sentClosingFairs);
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
    }
  }

  // The close snapshots that changed, CLOSE_SNAPSHOT_ROWS_PER_POST to a
  // request so the body stays far under the service's size cap however many
  // bets are open; a request that fails keeps its rows (and the rest) for the
  // next pass, and a 400 (a row the service refuses) drops only its own rows.
  async function postCloseSnapshots(pass) {
    const rows = edgesnapshot.unsentCloseRows(edgesnapshot.closeSnapshotRows({
      fairRows: pass.rows, records: held.records, rowsByBet: pass.rowsByBet, state: feedState, linesOfEvent: linesOfEventNow,
    }), sentCloseSnapshots);
    closefair.forgetClosed(sentCloseSnapshots, held.records);
    for (let start = 0; start < rows.length; start += CLOSE_SNAPSHOT_ROWS_PER_POST) {
      const chunk = rows.slice(start, start + CLOSE_SNAPSHOT_ROWS_PER_POST);
      try {
        await sendServiceJson(serviceFetch, `${deps.betsServiceUrl}/close_snapshots.json`, "POST", { rows: chunk });
        edgesnapshot.markCloseSent(chunk, sentCloseSnapshots);
        if (closeSnapshotStatus.error) logInfo("close snapshots saving again");
        Object.assign(closeSnapshotStatus, { okAt: now(), error: null, rowsSent: closeSnapshotStatus.rowsSent + chunk.length });
      } catch (error) {
        if (closeSnapshotStatus.error !== error.message) logWarning(`close snapshots POST failed: ${error.message}`);
        closeSnapshotStatus.error = error.message;
        if (error.status === HTTP_BAD_REQUEST) {
          edgesnapshot.markCloseSent(chunk, sentCloseSnapshots);
          continue;
        }
        return;
      }
    }
  }

  // Each stored league's ladders, rebuilt only when its snapshot landed again;
  // a stale league (a CDN copy built long ago) is not read, as closefair does.
  function currentLadders() {
    const staleLeagues = new Set(scannerStatus.staleLeagues || []);
    const loadedAt = scannerStatus.leagueLoadedAt || {};
    const ladders = [];
    for (const leagueId of laddersnapshot.LADDER_LEAGUE_IDS) {
      if (!boardLeagueIds().includes(leagueId) || staleLeagues.has(leagueId)) continue;
      const version = scanner.getLeagueVersion(leagueId);
      let cached = ladderCache.get(leagueId);
      if (!cached || cached.version !== version) {
        cached = { version, ladders: laddersnapshot.leagueLadders({
          leagueId, linesByEvent: leagueLinesByEvent(leagueId), events: feedState.events, observedMs: loadedAt[leagueId],
        }) };
        ladderCache.set(leagueId, cached);
      }
      ladders.push(...cached.ladders);
    }
    return ladders;
  }

  // POST the ladder rows the service has not taken, LADDER_ROWS_PER_POST to
  // a request. A failed request keeps its rows (and the rest) for the next
  // pass; a 400 (a row the service refuses) drops only its own rows. Stored
  // only: nothing reads the ladders back.
  async function postLadders() {
    if (laddersBusy || !feedState || !scannerStatus) return;
    laddersBusy = true;
    try {
      const at = now();
      laddersnapshot.forgetStarted(sentLadders, at);
      const rows = laddersnapshot.unsentRows(currentLadders(), sentLadders, at);
      for (let start = 0; start < rows.length; start += LADDER_ROWS_PER_POST) {
        const chunk = rows.slice(start, start + LADDER_ROWS_PER_POST);
        try {
          await sendServiceJson(serviceFetch, `${deps.betsServiceUrl}/ladders.json`, "POST", { rows: chunk.map(laddersnapshot.wireRow) });
          laddersnapshot.markSent(chunk, sentLadders);
          if (ladderStatus.error) logInfo("ladders saving again");
          Object.assign(ladderStatus, { okAt: now(), error: null, rowsSent: ladderStatus.rowsSent + chunk.length });
        } catch (error) {
          if (ladderStatus.error !== error.message) logWarning(`ladders POST failed: ${error.message}`);
          ladderStatus.error = error.message;
          if (error.status !== HTTP_BAD_REQUEST) return;
          laddersnapshot.markSent(chunk, sentLadders);
        }
      }
    } catch (error) {
      if (ladderStatus.error !== error.message) logWarning(`ladder pass failed: ${error.stack || error.message}`);
      ladderStatus.error = error.message;
    } finally {
      laddersBusy = false;
    }
  }

  // Each logged league's edge lines, rebuilt only when its snapshot landed
  // again; a stale league is not read (its segments stay open, unchanged).
  function currentEdgeLines() {
    const staleLeagues = new Set(scannerStatus.staleLeagues || []);
    const loadedAt = scannerStatus.leagueLoadedAt || {};
    const edgeLines = [];
    const readLeagueIds = [];
    for (const leagueId of edgelog.EDGE_LOG_LEAGUE_IDS) {
      if (!boardLeagueIds().includes(leagueId) || staleLeagues.has(leagueId)) continue;
      const version = scanner.getLeagueVersion(leagueId);
      let cached = edgeLineCache.get(leagueId);
      if (!cached || cached.version !== version) {
        cached = { version, edgeLines: edgelog.leagueEdgeLines({
          leagueId, lines: scanner.getLeagueLines(leagueId), state: feedState, observedMs: loadedAt[leagueId],
        }) };
        edgeLineCache.set(leagueId, cached);
      }
      edgeLines.push(...cached.edgeLines);
      readLeagueIds.push(leagueId);
    }
    return { edgeLines, readLeagueIds };
  }

  // Queue the rows; a newer row of a segment replaces its unsent one and
  // moves to the back. Past edgeSegmentMaxUnsent the oldest go, logged.
  function queueEdgeRows(rows) {
    for (const row of rows) {
      const key = edgelog.segmentKeyOf(row);
      unsentEdgeRows.delete(key);
      unsentEdgeRows.set(key, row);
    }
    let dropped = 0;
    for (const key of unsentEdgeRows.keys()) {
      if (unsentEdgeRows.size <= edgeSegmentMaxUnsent) break;
      unsentEdgeRows.delete(key);
      dropped += 1;
    }
    if (!dropped) return;
    edgeLogStatus.dropped += dropped;
    logWarning(`edge segments: dropped the ${dropped} oldest unsent rows (more than ${edgeSegmentMaxUnsent} held)`);
  }

  function firstUnsentEdgeRows() {
    const chunk = [];
    for (const entry of unsentEdgeRows) {
      if (chunk.length >= EDGE_SEGMENT_ROWS_PER_POST) break;
      chunk.push(entry);
    }
    return chunk;
  }

  // The chunk index of the row a 400 names ("rows[N] ..."), or null.
  function refusedRowIndex(error, chunkLength) {
    const match = /rows\[(\d+)\]/.exec(error.message);
    const index = match ? Number(match[1]) : null;
    return index !== null && index < chunkLength ? index : null;
  }

  // Advance the edge segments and POST the rows the service has not taken,
  // EDGE_SEGMENT_ROWS_PER_POST to a request. A failed request keeps its rows
  // (and the rest) for the next pass. A 400 drops only the row it names
  // (the service refuses the whole request over one bad row, which would be
  // refused again), logged and counted, and the rest are sent again; one that
  // names no row drops its chunk. Stored only: nothing reads the segments back.
  async function postEdgeSegments() {
    if (edgeLogBusy || !feedState || !scannerStatus) return;
    edgeLogBusy = true;
    try {
      const { edgeLines, readLeagueIds } = currentEdgeLines();
      queueEdgeRows(edgelog.advanceSegments({ edgeLines, readLeagueIds, segments: edgeSegments, now: now() }));
      while (unsentEdgeRows.size) {
        const chunk = firstUnsentEdgeRows();
        // Taken, or refused: a row the service refuses would be refused again.
        let finished = chunk;
        try {
          await sendServiceJson(serviceFetch, `${deps.betsServiceUrl}/edge_segments.json`, "POST", { rows: chunk.map(([, row]) => row) });
          if (edgeLogStatus.error) logInfo("edge segments saving again");
          Object.assign(edgeLogStatus, { okAt: now(), error: null, rowsSent: edgeLogStatus.rowsSent + chunk.length });
        } catch (error) {
          if (edgeLogStatus.error !== error.message) logWarning(`edge segments POST failed: ${error.message}`);
          edgeLogStatus.error = error.message;
          if (error.status !== HTTP_BAD_REQUEST) return;
          const index = refusedRowIndex(error, chunk.length);
          if (index !== null) finished = [chunk[index]];
          edgeLogStatus.dropped += finished.length;
          logWarning(`edge segments: dropped ${finished.length} row(s) the service refused`);
        }
        for (const [key, row] of finished) if (unsentEdgeRows.get(key) === row) unsentEdgeRows.delete(key);
      }
    } catch (error) {
      if (edgeLogStatus.error !== error.message) logWarning(`edge segment pass failed: ${error.stack || error.message}`);
      edgeLogStatus.error = error.message;
    } finally {
      edgeLogBusy = false;
    }
  }

  // One Kalshi orderbook, parsed (taker.parseKalshiOrderbook). Throws naming
  // the ticker and what went wrong; a 404 is a ticker Kalshi does not list.
  async function readKalshiBook(ticker) {
    const body = await getServiceJson(kalshiFetch, `${kalshiBaseUrl}/markets/${encodeURIComponent(ticker)}/orderbook`);
    return taker.parseKalshiOrderbook(body);
  }

  // id -> {book} | {error} for every id, KALSHI_BOOK_READS_AT_ONCE reads of
  // `readOne(id)` in flight (Kalshi tickers; Novig's books come off the stream).
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
        sizing: {
          matches: flags[index].matches, ladderOf: ladderReaderOf(row.eventId), teasers: pass.teasers,
          shrink: shrinkOfRow(row), minShrunkEdgePct: settings.edgeSettings.minShrunkEdgePct,
        },
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
    keepTakerKalshiReads(rows, pass.board, books, now());
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
    const books = await readStreamBooks(Array.from(marketIds), false);
    if (!catalogError) keepTakerNovigReads(rows, resolve, index, books, now());
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

  // ---- the live price overlay (extension/liveprices.js) ---------------------
  //
  // Unabated's copies of the exchanges' prices lag, sometimes for hours, so
  // every Edges line at Kalshi, Novig, Polymarket US and ProphetX is read off the
  // exchange's own book here and the reads are laid over the board for the
  // Edges list (/edges.json here, /live_prices.json for the panel). The
  // taker decides on its own reads as before (Cal, 2026-11-02: the bot stays
  // as it is); the reads it already made count here, so the 3%+ lines are not
  // read twice. Each venue's line is shown at its book only when the book's
  // market is the line's (Kalshi's market check, Novig's listed market and
  // outcome name, Polymarket US's own sides): a rung Unabated filed under the
  // wrong contract stays unchecked rather than painting that contract's price.

  // The venue the pass reads a board book at, or null for a book it does not read.
  function liveVenueOf(board, bookId) {
    if (bookId === taker.KALSHI_BOOK_ID) return taker.VENUE_KALSHI;
    if (bookId === taker.NOVIG_BOOK_ID) return taker.VENUE_NOVIG;
    if (polymarket.isPolymarketUsBook(board.books[bookId])) return polymarket.VENUE;
    if (bookId === prophetx.PROPHETX_BOOK_ID) return VENUE_PROPHETX;
    return null;
  }

  // The board books the pass reads. `served`: only those whose reads are laid
  // over the board, so ProphetX is left out (its lines untagged, at
  // Unabated's price) until the bets service answers on the live key pair.
  function liveBookIdsOf(board, served = false) {
    return Object.values(board.books).map((book) => book.id)
      .filter((id) => liveVenueOf(board, id) !== null && !(served && id === prophetx.PROPHETX_BOOK_ID && !prophetxLiveOn))
      .sort((a, b) => a - b);
  }

  // The lines the Edges list could show at the read venues: the list's own
  // selection (its settings, leagues, markets, alts, line age, minimum edge)
  // by Unabated's price, with no liquidity gate (the live depth decides it).
  function liveCandidateRows(board, at) {
    const edgeSettings = settings.edgeSettings;
    const effective = edgeRows.effectiveFilter(edgeSettings, board, null, null);
    const bookIds = liveBookIdsOf(board).filter((id) => (effective.bookIds ? effective.bookIds.has(id) : board.books[id].isLive));
    if (!bookIds.length) return [];
    return feedLib.selectEdges(board, {
      ...edgeRows.edgeSelectionOptions(edgeSettings, effective, at), bookIds: new Set(bookIds),
      minLiquidityToWin: null, minEdge: edgeSettings.minEdgePct / 100,
    });
  }

  // One line's read: quote {allIn, depthDollars}, or null when nothing rests to buy.
  function recordLiveRead(row, venue, quote, at) {
    const points = row.points ?? null;
    liveReads.set(row.key, quote
      ? { venue, points, allIn: Math.round(quote.allIn * LIVE_PRICE_DECIMALS) / LIVE_PRICE_DECIMALS, depthDollars: Math.round(quote.depthDollars * 100) / 100, readAtMs: at }
      : { venue, points, empty: true, readAtMs: at });
  }

  // What buying costs off each venue's quote, as {allIn, depthDollars} (or null: nothing rests).
  function kalshiLiveQuote(book, side) {
    const quote = taker.quoteOf(book, side);
    if (!quote) return null;
    const allIn = taker.kalshiAllInPrice(quote.ask);
    return { allIn, depthDollars: quote.depth * allIn };
  }

  function novigLiveQuote(book, contract) {
    const quote = taker.novigQuoteOf(book, contract);
    if (!quote) return null;
    return { allIn: quote.ask, depthDollars: (quote.depth * quote.ask) / taker.NOVIG_QTY_PER_DOLLAR };
  }

  // Kalshi's market check for a row off the cached market read: the mismatch
  // (null when the market is the row's), or undefined when no fresh read is held.
  function kalshiCheckOf(row, contract, at) {
    const check = kalshiMarketChecks.get(contract.ticker);
    if (!check || at - check.checkedAt > KALSHI_MARKET_CHECK_TTL_MS) return undefined;
    return taker.marketMismatch(check.market, row, contract);
  }

  // The row's Novig market off the catalog listing, or why it is not the row's
  // line (a reason the overlay leaves the row unchecked for).
  function novigListedMismatch(row, contract, index) {
    const listed = index.get(contract.outcomeId);
    if (!listed) return `outcome ${contract.outcomeId} is on no listed market`;
    if (listed.market.status !== taker.NOVIG_MARKET_OPEN) return `market status "${listed.market.status}"`;
    return taker.novigSideMismatch(listed.market, row, contract);
  }

  // The taker's own Kalshi reads, kept for the overlay where Kalshi's market check already passed for the row.
  function keepTakerKalshiReads(rows, board, books, at) {
    for (const row of rows) {
      const contract = taker.resolveKalshiContract(row, board);
      const read = contract.ticker ? books.get(contract.ticker) : null;
      if (!read || !read.book || kalshiCheckOf(row, contract, at) !== null) continue;
      recordLiveRead(row, taker.VENUE_KALSHI, kalshiLiveQuote(read.book, contract.side), at);
    }
  }

  // The taker's own Novig reads, kept for the overlay where the listed market is the row's.
  function keepTakerNovigReads(rows, resolve, index, books, at) {
    for (const row of rows) {
      const contract = resolve(row);
      const read = contract.marketId ? books.get(contract.marketId) : null;
      if (!read || read.error || novigListedMismatch(row, contract, index) !== null) continue;
      try {
        recordLiveRead(row, taker.VENUE_NOVIG, novigLiveQuote(read.book, contract), at);
      } catch (_error) {
        // A malformed book is the taker's no_book; the overlay leaves the row unchecked.
      }
    }
  }

  // One live price pass: every Edges line at the read venues that has no read
  // newer than LIVE_REREAD_MS, oldest read first (a line never read first,
  // best edge first), each venue's book read and kept. Kalshi and Polymarket
  // US are read per market at a few at a time, at most LIVE_*_READS_PER_PASS
  // a pass, so the taker's and Place's reads keep their room; Novig's books
  // come off the stream, subscribed behind the bot's and Place's markets, and
  // off the keyed REST /book for a market the stream has not reached yet.
  async function livePricePass() {
    const board = feedState;
    if (liveBusy || !board || !scannerStatus || scannerStatus.phase !== "live") return;
    liveBusy = true;
    const startedAt = now();
    try {
      for (const [key, read] of liveReads) {
        if (startedAt - read.readAtMs > livePrices.LIVE_READ_MAX_AGE_MS) liveReads.delete(key);
      }
      for (const [ticker, check] of kalshiMarketChecks) {
        if (startedAt - check.checkedAt > KALSHI_MARKET_CHECK_TTL_MS) kalshiMarketChecks.delete(ticker);
      }
      for (const [query, list] of polymarketEventLists) {
        if (startedAt - list.readAt > POLYMARKET_EVENTS_REUSE_MS) polymarketEventLists.delete(query);
      }
      for (const [line, contract] of polymarketLineContracts) {
        if (startedAt >= contract.startMs) polymarketLineContracts.delete(line);
      }
      const readAtOf = (row) => {
        const read = liveReads.get(row.key);
        return read && (read.points ?? null) === (row.points ?? null) ? read.readAtMs : -Infinity;
      };
      const rows = liveCandidateRows(board, startedAt);
      const due = rows.filter((row) => startedAt - readAtOf(row) > LIVE_REREAD_MS)
        .sort((a, b) => readAtOf(a) - readAtOf(b));
      const ofVenue = (venue) => due.filter((row) => liveVenueOf(board, row.book.id) === venue);
      const [kalshiCounts, novigCounts, polymarketCounts, prophetxCounts] = await Promise.all([
        liveKalshiReads(board, ofVenue(taker.VENUE_KALSHI).slice(0, LIVE_KALSHI_READS_PER_PASS)),
        liveNovigReads(board, ofVenue(taker.VENUE_NOVIG)),
        livePolymarketReads(ofVenue(polymarket.VENUE).slice(0, LIVE_POLYMARKET_READS_PER_PASS)),
        // Every ProphetX row, due or not: the bets service answers its last reads at once and reads on its own pace.
        liveProphetxReads(board, rows.filter((row) => liveVenueOf(board, row.book.id) === VENUE_PROPHETX)),
      ]);
      liveLastPass = {
        startedAt, finishedAt: now(), rows: rows.length, due: due.length, error: null,
        byVenue: { kalshi: kalshiCounts, novig: novigCounts, polymarket_us: polymarketCounts, prophetx: prophetxCounts },
      };
    } catch (error) {
      logWarning(`live price pass failed: ${error.stack || error.message}`);
      liveLastPass = { ...liveLastPass, startedAt, finishedAt: now(), error: error.message };
    } finally {
      liveBusy = false;
    }
  }

  function newLiveCounts(rows) {
    return { rows, live: 0, empty: 0, unchecked: 0, lastError: null };
  }

  function noteUnchecked(counts, why) {
    counts.unchecked += 1;
    counts.lastError = why;
  }

  function noteRead(counts, row, venue, quote, at) {
    recordLiveRead(row, venue, quote, at);
    if (quote) counts.live += 1;
    else counts.empty += 1;
  }

  async function liveKalshiReads(board, rows) {
    const counts = newLiveCounts(rows.length);
    const at = now();
    const wanted = [];
    for (const row of rows) {
      const contract = taker.resolveKalshiContract(row, board);
      if (contract.reason) noteUnchecked(counts, `${row.key}: ${contract.detail}`);
      else wanted.push({ row, contract });
    }
    const unchecked = new Set(wanted.filter(({ row, contract }) => kalshiCheckOf(row, contract, at) === undefined).map(({ contract }) => contract.ticker));
    const markets = await readBooks(unchecked, readKalshiMarket);
    for (const [ticker, read] of markets) {
      if (read.book) kalshiMarketChecks.set(ticker, { market: read.book, checkedAt: now() });
    }
    const checked = wanted.filter(({ row, contract }) => {
      const mismatch = kalshiCheckOf(row, contract, now());
      if (mismatch === null) return true;
      const read = markets.get(contract.ticker);
      noteUnchecked(counts, `${contract.ticker}: ${mismatch === undefined ? `market unread: ${read ? read.error : "no read"}` : mismatch}`);
      return false;
    });
    const books = await readBooks(new Set(checked.map(({ contract }) => contract.ticker)), readKalshiBook);
    const readAt = now();
    for (const { row, contract } of checked) {
      const read = books.get(contract.ticker);
      if (!read || !read.book) noteUnchecked(counts, `${contract.ticker}: ${read ? read.error : "no read"}`);
      else noteRead(counts, row, taker.VENUE_KALSHI, kalshiLiveQuote(read.book, contract.side), readAt);
    }
    return counts;
  }

  async function liveNovigReads(board, rows) {
    const counts = newLiveCounts(rows.length);
    if (!rows.length) return counts;
    let index;
    try {
      index = (await novigCatalogFor(rows, NOVIG_CATALOG_REUSE_MS)).index;
    } catch (error) {
      counts.unchecked = rows.length;
      counts.lastError = `Novig's catalog unread: ${error.message}`;
      return counts;
    }
    const wanted = [];
    for (const row of rows) {
      const contract = row.betTypeId === MONEYLINE_BET_TYPE_ID ? taker.novigMoneylineContract(row, index, board) : taker.novigContractOf(row, index);
      const mismatch = contract.reason ? contract.detail : novigListedMismatch(row, contract, index);
      if (mismatch) noteUnchecked(counts, `${row.key}: ${mismatch}`);
      else wanted.push({ row, contract });
    }
    const marketIds = Array.from(new Set(wanted.map(({ contract }) => contract.marketId)));
    const books = new Map();
    for (let start = 0; start < marketIds.length; start += NOVIG_BOOK_MARKETS_PER_ASK) {
      const read = await readStreamBooks(marketIds.slice(start, start + NOVIG_BOOK_MARKETS_PER_ASK), false, true);
      for (const [marketId, view] of read) books.set(marketId, view);
    }
    const restBooks = await readNovigRestBooks(marketIds.filter((marketId) => {
      const read = books.get(marketId);
      return !read || (read.error && !read.final);
    }));
    const readAt = now();
    for (const { row, contract } of wanted) {
      const read = restBooks.get(contract.marketId) || books.get(contract.marketId);
      if (!read) {
        noteUnchecked(counts, `${contract.marketId}: no stream answer`);
      } else if (read.error) {
        // The stream saw the market go live or leave OPEN: nothing there to buy pregame.
        if (read.final) noteRead(counts, row, taker.VENUE_NOVIG, null, readAt);
        else noteUnchecked(counts, read.error);
      } else {
        try {
          noteRead(counts, row, taker.VENUE_NOVIG, novigLiveQuote(read.book, contract), readAt);
        } catch (error) {
          noteUnchecked(counts, `${contract.marketId}: ${error.message}`);
        }
      }
    }
    return counts;
  }

  // The keyed REST books (novig_read.json) of the first
  // LIVE_NOVIG_REST_READS_PER_PASS markets the stream has no book of:
  // marketId -> {book} | {error}. The rest wait a pass (the stream may have
  // them by then). The catalog the rows were matched on lists only open
  // pregame markets, read at most NOVIG_CATALOG_REUSE_MS ago.
  //
  // One read at a time, so of novig_read.json's four sessions three stay free
  // for the taker's catalog read and a Place quote's market and event; the
  // first failed read (a 429, a slow service) stops the rest of the pass.
  async function readNovigRestBooks(marketIds) {
    const books = new Map();
    for (const marketId of marketIds.slice(0, LIVE_NOVIG_REST_READS_PER_PASS)) {
      try {
        books.set(marketId, { book: await getNovigJson(`${novigMarketPath(marketId)}/book`) });
      } catch (error) {
        books.set(marketId, { error: `Novig's keyed book of market ${marketId} unread: ${error.message}` });
        break;
      }
    }
    return books;
  }

  // Polymarket US's live prices: each row's contract found once in its event
  // list (reused POLYMARKET_EVENTS_REUSE_MS, at most LIVE_POLYMARKET_EVENT_READS_PER_PASS
  // fresh lists a pass) and then remembered until its game starts (a line's
  // market does not change; a closed one reads as an empty book), then every
  // contract's book in ONE ask of the bets service's book stream
  // (background: behind Place's markets, no wait). No
  // book is read off the public gateway here: its 429s once refused Place.
  async function livePolymarketReads(rows) {
    const counts = newLiveCounts(rows.length);
    const found = [];
    const queue = [];
    for (const row of rows) {
      const known = polymarketLineContracts.get(polymarketLineKey(row));
      if (known) found.push({ row, slug: known.slug, side: known.side });
      else queue.push(row);
    }
    // One read per list a pass, shared by the rows that need it, so two lookups of one game spend one read.
    const freshEventReads = new Map();
    const readEventsWithinBudget = (marketType, startMs) => {
      const key = `${marketType}|${startMs}`;
      if (freshEventReads.has(key)) return freshEventReads.get(key);
      if (polymarketEventListHeld(marketType, startMs)) return readPolymarketEventsReused(marketType, startMs);
      if (freshEventReads.size >= LIVE_POLYMARKET_EVENT_READS_PER_PASS) {
        return Promise.reject(new Error(`this pass already read ${LIVE_POLYMARKET_EVENT_READS_PER_PASS} Polymarket US event lists; next pass`));
      }
      const read = readPolymarketEventsReused(marketType, startMs);
      freshEventReads.set(key, read);
      return read;
    };
    const lookUp = async (row) => {
      const contract = await polymarketContractOfRow(row, readEventsWithinBudget);
      if (contract.reason) return noteUnchecked(counts, `${row.key}: ${contract.detail}`);
      const { slug, side } = contract.contract;
      polymarketLineContracts.set(polymarketLineKey(row), { slug, side, startMs: contract.startMs });
      found.push({ row, slug, side });
    };
    const worker = async () => {
      while (queue.length) await lookUp(queue.shift());
    };

    await Promise.all(Array.from({ length: POLYMARKET_CONTRACT_LOOKUPS_AT_ONCE }, worker));
    const slugs = Array.from(new Set(found.map((entry) => entry.slug))).slice(0, POLYMARKET_BOOK_SLUGS_PER_ASK);
    const books = await readPolymarketStreamBooks(slugs, { background: true });
    for (const { row, slug, side } of found) {
      const read = books.get(slug) || { error: `${slug}: more than ${POLYMARKET_BOOK_SLUGS_PER_ASK} markets in one pass` };
      if (read.error) {
        noteUnchecked(counts, `${slug}: ${read.error}`);
        continue;
      }
      const quote = polymarket.quoteOf(read.book, side);
      if (quote.reason === "no_book" || quote.reason === "market_closed") {
        noteRead(counts, row, polymarket.VENUE, null, now());
        continue;
      }
      if (quote.reason) {
        noteUnchecked(counts, `${slug}: ${quote.detail}`);
        continue;
      }
      const allIn = polymarket.allInOf(quote.ask);
      noteRead(counts, row, polymarket.VENUE, { allIn, depthDollars: quote.depth * allIn }, now());
    }
    return counts;
  }

  // One Polymarket US line: a moved number or start is another line, found afresh.
  function polymarketLineKey(row) {
    return `${row.key}|${row.points ?? ""}|${row.eventStart}`;
  }

  // Polymarket US's books off the bets service's stream (POST
  // /polymarket_books.json): slug -> {book: {marketData}} (the REST /book's
  // body, so polymarket.quoteOf reads either) | {error}. The stream
  // unreachable, off or down is an error on every slug, never an old book.
  //   priority    the Place quote: subscribed first, waiting POLYMARKET_STREAM_WAIT_MS for a new market
  //   background  the live price pass: behind the rest, no wait
  async function readPolymarketStreamBooks(slugs, { priority = false, background = false } = {}) {
    const books = new Map();
    if (!slugs.length) return books;
    let reply;
    try {
      const body = { slugs, waitMs: background ? 0 : POLYMARKET_STREAM_WAIT_MS, priority, background };
      reply = await sendServiceJson(serviceFetch, `${deps.betsServiceUrl}/polymarket_books.json`, "POST", body,
        POLYMARKET_STREAM_READ_TIMEOUT_MS);
    } catch (error) {
      for (const slug of slugs) books.set(slug, { error: `Polymarket US stream books unread: ${error.message}` });
      return books;
    }
    for (const slug of slugs) {
      const view = reply && reply.books && typeof reply.books === "object" ? reply.books[slug] : null;
      if (view && view.marketData && typeof view.marketData === "object") books.set(slug, { book: { marketData: view.marketData }, ageMs: view.ageMs });
      else books.set(slug, { error: view && typeof view.error === "string" ? view.error : "the stream answered no book for it" });
    }
    return books;
  }

  // ProphetX's reads off the bets service (POST /prophetx_live.json,
  // prophetx_live.py, its own key pair): each game the rows are on is named
  // with the ProphetX strike ids Unabated carries for it, and each row is
  // matched to the one ProphetX side that reads as its line (the line's own
  // strike id, when Unabated names one, must be that side's). The service
  // never waits: a game it has not read yet leaves its rows unchecked.
  async function liveProphetxReads(board, rows) {
    const counts = newLiveCounts(rows.length);
    if (!rows.length) return counts;
    const rowsByGame = new Map();
    for (const row of rows) {
      if (!rowsByGame.has(row.eventId)) rowsByGame.set(row.eventId, []);
      rowsByGame.get(row.eventId).push(row);
    }
    const games = Array.from(rowsByGame, ([eventId, gameRows]) => ({
      gameKey: String(eventId), eventStart: gameRows[0].eventStart, line: taker.novigLineOf(gameRows[0]),
      strikeIds: prophetx.strikeIdsOf(board.lines[gameRows[0].key] || {}, linesOfEvent(board, eventId)).strikeIds,
    })).slice(0, PROPHETX_LIVE_GAMES_PER_ASK);
    let reply;
    try {
      reply = await sendServiceJson(serviceFetch, `${deps.betsServiceUrl}/prophetx_live.json`, "POST", { games }, SERVICE_TIMEOUT_MS);
      prophetxLiveOn = true;
    } catch (error) {
      if (error.status === 503) prophetxLiveOn = false;
      counts.unchecked = rows.length;
      counts.lastError = `ProphetX's live prices unread: ${error.message}`;
      return counts;
    }
    const gameReads = reply && reply.games && typeof reply.games === "object" ? reply.games : {};
    for (const row of rows) {
      const read = gameReads[String(row.eventId)];
      const matched = prophetxLiveSideOf(board, row, read);
      if (matched.reason) {
        noteUnchecked(counts, `${row.key}: ${matched.reason}`);
        continue;
      }
      const { side, readAtMs } = matched;
      const quote = side.price == null ? null : { allIn: pricecheck.impliedProb(side.price), depthDollars: side.depth };
      noteRead(counts, row, VENUE_PROPHETX, quote, readAtMs);
    }
    return counts;
  }

  // The one side of a ProphetX game read that is the row's line, with the
  // read's time, or {reason} why the row stays unchecked.
  function prophetxLiveSideOf(board, row, read) {
    if (!read) return { reason: "ProphetX game not read yet" };
    if (read.error) return { reason: read.error };
    const readAtMs = Date.parse(read.readAt);
    if (!Number.isFinite(readAtMs) || !Array.isArray(read.lines)) return { reason: "ProphetX read has no time or lines" };
    const want = taker.novigLineOf(row);
    if (read.totalsOnly && want.betType !== "total") return { reason: "ProphetX lists the teams the other way round: totals only" };
    const points = want.betType === "moneyline" ? null : want.points;
    const found = read.lines.filter((side) => side.betType === want.betType && side.period === want.period
      && side.side === want.side && (side.points ?? null) === points);
    if (found.length !== 1) return { reason: `${found.length} ProphetX sides read as ${want.side} ${points ?? ""} ${want.period}` };
    const ownStrike = (board.lines[row.key] || {}).sourceData;
    if (prophetx.isProphetxStrikeId(ownStrike) && ownStrike !== found[0].strikeId) {
      return { reason: `Unabated's ProphetX strike ${ownStrike} is not the one that reads as this line (${found[0].strikeId})` };
    }
    if (found[0].price != null && pricecheck.impliedProb(found[0].price) == null) return { reason: `ProphetX price ${found[0].price} is no American price` };
    return { side: found[0], readAtMs };
  }

  // The /live_prices.json body (liveprices.createLineView's input): the read
  // venues' book ids on the board and every line's last read, for the panel.
  function livePricesPayloadNow() {
    const board = feedState;
    const reads = {};
    for (const [key, read] of liveReads) {
      reads[key] = read.empty
        ? { venue: read.venue, points: read.points, empty: true, readAtMs: read.readAtMs }
        : { venue: read.venue, points: read.points, allIn: read.allIn, depthDollars: read.depthDollars, readAtMs: read.readAtMs };
    }
    return {
      generatedAt: new Date(now()).toISOString(), maxAgeMs: livePrices.LIVE_READ_MAX_AGE_MS,
      bookIds: board ? liveBookIdsOf(board, true) : [], reads, lastPass: liveStatus(),
    };
  }

  function liveStatus() {
    return {
      startedAt: isoOrNull(liveLastPass.startedAt), finishedAt: isoOrNull(liveLastPass.finishedAt),
      rows: liveLastPass.rows, due: liveLastPass.due ?? null, reads: liveReads.size, error: liveLastPass.error, byVenue: liveLastPass.byVenue,
    };
  }

  // The view /edges.json reads the board through: the reads held now, or none before the board loads.
  function liveLineView() {
    return feedState ? livePrices.createLineView(livePricesPayloadNow(), now()) : null;
  }

  // Novig's live books off the bets service's stream: marketId -> {book} |
  // {error, final?} (taker.novigStreamBookRead) for every market id. The
  // stream unreachable, off or down is an error on every market (never a
  // stale book). `priority` (the Place quote) subscribes ahead of the bot's.
  //   background  the live price pass: subscribed after the bot's and Place's markets, no wait for a snapshot
  async function readStreamBooks(marketIds, priority, background = false) {
    const books = new Map();
    if (!marketIds.length) return books;
    let reply;
    try {
      const body = background ? { marketIds, waitMs: 0, priority: false, background: true } : { marketIds, waitMs: NOVIG_STREAM_WAIT_MS, priority };
      reply = await sendServiceJson(serviceFetch, `${deps.betsServiceUrl}/novig_books.json`, "POST", body, NOVIG_STREAM_READ_TIMEOUT_MS);
    } catch (error) {
      for (const marketId of marketIds) books.set(marketId, { error: `Novig's stream books unread: ${error.message}` });
      return books;
    }
    for (const marketId of marketIds) books.set(marketId, taker.novigStreamBookRead(reply, marketId));
    return books;
  }

  // One keyed Novig catalog GET (through the bets service) as JSON; throws
  // naming the path, Novig's status and its error code (a 429 reads as one).
  async function getNovigJson(path, query = null) {
    const { status, body } = await novigRead(path, query);
    if (status === 200) return body;
    const code = body && typeof body === "object" && body.code ? ` ${body.code}` : "";
    throw new Error(`HTTP ${status}${code} from Novig's keyed catalog GET ${path}`);
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
    rememberNovigEvents(index);
    return { index, readAt, reused: false };
  }

  // The in-game Novig catalog of one live row's league, reused for maxAgeMs (a live Place only).
  async function novigLiveCatalogFor(row, maxAgeMs) {
    const league = NOVIG_LEAGUE_OF[row.league];
    const cached = novigLiveCatalog;
    if (maxAgeMs > 0 && league && cached && now() - cached.readAt <= maxAgeMs && cached.leagues.has(league)) {
      return { index: cached.index, readAt: cached.readAt, reused: true };
    }
    const readAt = now();
    const index = taker.novigOutcomeIndex(await readNovigMarkets([row], taker.NOVIG_EVENT_INGAME));
    novigLiveCatalog = { leagues: new Set(league ? [league] : []), readAt, index };
    return { index, readAt, reused: false };
  }

  // Each board game's Novig event off a fresh pregame catalog read, kept for a
  // live Place once the game starts (its last pregame reading wins).
  function rememberNovigEvents(index) {
    if (!feedState) return;
    for (const [eventId, event] of Object.entries(feedState.events)) {
      const outcomes = event.venueIds && event.venueIds.novigOutcomes;
      if (!outcomes || Object.keys(outcomes).length === 0) continue;
      const memory = taker.novigEventMemoryOf(event, feedState, index);
      if (memory) novigEventMemo.set(Number(eventId), memory);
    }
    for (const eventId of novigEventMemo.keys()) if (!kickoffMemo.has(eventId) && !(feedState.events[eventId])) novigEventMemo.delete(eventId);
  }

  // Every open Novig market of `eventStatus` (pregame by default) in the rows' leagues, every page (never a partial list).
  async function readNovigMarkets(rows, eventStatus = taker.NOVIG_EVENT_PREGAME) {
    const leagues = Array.from(new Set(rows.map((row) => NOVIG_LEAGUE_OF[row.league]).filter(Boolean))).sort();
    if (!leagues.length) return [];
    const markets = [];
    let after = null;
    for (let page = 0; page < NOVIG_MARKETS_MAX_PAGES; page += 1) {
      const query = { league: leagues.join(","), eventStatus, limit: String(NOVIG_MARKETS_PAGE_LIMIT) };
      if (after) query.after = after;
      const body = await getNovigJson(NOVIG_CATALOG_MARKETS_PATH, query);
      if (!body || !Array.isArray(body.items)) throw new Error("Novig markets: no items list");
      markets.push(...body.items);
      after = body.next || null;
      if (!after) return markets;
    }
    throw new Error(`Novig markets: more than ${NOVIG_MARKETS_MAX_PAGES} pages`);
  }

  // The board a Place quote finds its line on: the pregame board, or for a
  // live Place (spec.live) the live board of the league whose last snapshot
  // has the game in play (scanner keepLive; its rows never reach the Edges
  // list or the bot). {board} or {board: null, reason, detail}.
  function quoteBoardOf(spec) {
    if (!spec.live) {
      return feedState ? { board: feedState } : { board: null, reason: "no_board", detail: "the board has not loaded yet" };
    }
    const leagueId = scanner.liveLeagueOfEvent(spec.eventId);
    const live = leagueId == null ? null : scanner.getLiveBoard(leagueId);
    if (!live) return { board: null, reason: "not_live", detail: `game ${spec.eventId} is not on Unabated's live board` };
    return { board: liveboard.liveQuoteBoard(live.state, spec.eventId, pregameVenueMemo.get(spec.eventId) || null) };
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
    const { board, reason, detail } = quoteBoardOf(spec);
    if (!board) return { ok: false, reason, detail };
    const found = linesOfEvent(board, spec.eventId).filter((line) => line.bookId === taker.KALSHI_BOOK_ID
      && line.eventId === spec.eventId && line.periodTypeId === spec.periodTypeId && line.betTypeId === spec.betTypeId
      && line.sideIndex === spec.sideIndex && (spec.betTypeId === MONEYLINE_BET_TYPE_ID || (line.points ?? null) === spec.points));
    if (found.length !== 1) {
      return { ok: false, reason: "no_line", detail: `the board holds ${found.length} Kalshi lines for game ${spec.eventId} at that market, side and number` };
    }
    const row = feedLib.describeLine(found[0], board);
    const fields = { eventId: row.eventId, eventStart: row.eventStart, label: `${row.sideLabel} ${row.period}`.trim(), eventName: row.eventName, live: spec.live };
    const contract = await kalshiPlaceContractOf(row, board, spec.live);
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
    // exchangeIndex: the Kalshi cash pot (shard) the market draws on, so the
    // bets service can top that pot up before ordering (null when Kalshi omits it).
    return {
      ok: true, ...position, ask: quote.ask, allIn: Math.round(taker.kalshiAllInPrice(quote.ask) * 1e6) / 1e6, depth: quote.depth,
      exchangeIndex: Number.isInteger(market.exchange_index) ? market.exchange_index : null, kalshiBaseUrl,
    };
  }

  // The Kalshi contract Place buys: the bot's own reading of the rungs, and
  // for a moneyline whose team code no full-game spread rung gives (a game
  // without one, a half), Kalshi's own market list for the game's winner
  // event (taker.kalshiMoneylineEventOf / kalshiMoneylineFromMarkets).
  async function kalshiPlaceContractOf(row, board, live = false) {
    const contract = taker.resolveKalshiContract(row, board);
    // Live, the favorite may have changed sides of zero since the pregame ladder the ids came from.
    if (contract.reason && live && row.betTypeId === SPREAD_BET_TYPE_ID) {
      const derived = taker.derivedSpreadFromTeamCodes(row, board);
      return derived.reason ? { reason: derived.reason, detail: `${contract.detail}; ${derived.detail}` } : derived;
    }
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
      + `${answer.ok ? `ok @ ${answer.price}` : `${answer.reason} (${answer.detail})`} in ${now() - startedAt} ms`
      + (catalog ? ` (catalog ${catalog.reused ? `reused, ${Math.round((startedAt - catalog.readAt) / 1000)}s old` : "read"}, `
        + `market reads ${catalog.readsMs} ms, book from ${catalog.bookSource})` : "")
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
    const { board, reason: boardReason, detail: boardDetail } = quoteBoardOf(spec);
    if (!board) return { answer: { ok: false, reason: boardReason, detail: boardDetail }, catalog: null };
    const found = linesOfEvent(board, spec.eventId).filter((line) => line.bookId === taker.NOVIG_BOOK_ID
      && line.eventId === spec.eventId && line.periodTypeId === spec.periodTypeId && line.betTypeId === spec.betTypeId
      && line.sideIndex === spec.sideIndex && (spec.betTypeId === MONEYLINE_BET_TYPE_ID || (line.points ?? null) === spec.points));
    if (found.length !== 1) {
      return { answer: { ok: false, reason: "no_line", detail: `the board holds ${found.length} Novig lines for game ${spec.eventId} at that market, side and number` }, catalog: null };
    }
    const row = feedLib.describeLine(found[0], board);
    const fields = { eventId: row.eventId, eventStart: row.eventStart, label: `${row.sideLabel} ${row.period}`.trim(), eventName: row.eventName, live: spec.live };
    let catalog;
    let contract;
    try {
      // A moneyline has no Novig rung of its own: Place finds it off the game's other rungs (never the bot).
      // Live, no line carries a Novig id: the game's remembered Novig event names the market (novigLiveContract).
      const contractOf = (index) => (spec.live ? taker.novigLiveContract(row, index, novigEventMemo.get(row.eventId) || null)
        : row.betTypeId === MONEYLINE_BET_TYPE_ID ? taker.novigMoneylineContract(row, index, board) : taker.novigContractOf(row, index));
      const catalogFor = (maxAgeMs) => (spec.live ? novigLiveCatalogFor(row, maxAgeMs) : novigCatalogFor([row], maxAgeMs));
      catalog = await catalogFor(NOVIG_CATALOG_REUSE_MS);
      contract = contractOf(catalog.index);
      // A market listed since the reused read: one fresh read before giving up.
      if (contract.reason && catalog.reused) {
        catalog = await catalogFor(0);
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
    // The book stream drops a market once its game starts, so a live Place reads the book over keyed REST.
    const [streamRead, { mismatch, market: checkedMarket }] = await Promise.all([
      spec.live ? Promise.resolve({ error: "a live Place reads the book over keyed REST" })
        : readStreamBooks([contract.marketId], true).then((books) => books.get(contract.marketId)),
      readNovigMarketCheck(row, contract, listedEventId, spec.live),
    ]);
    // Without a stream book (the stream down or still subscribing) Place reads the book once over keyed REST;
    // the stream seeing the market live or not OPEN is final.
    const bookRead = streamRead.book || streamRead.final ? streamRead
      : await getNovigJson(`${novigMarketPath(contract.marketId)}/book`).then((book) => ({ book }), (error) => ({ error }));
    const timing = {
      reused: catalog.reused, readAt: catalog.readAt, readsMs: now() - readsStartedAt,
      bookSource: streamRead.book || streamRead.final ? "stream" : `keyed REST (stream: ${streamRead.error})`,
    };
    let quote;
    try {
      if (bookRead.error) throw typeof bookRead.error === "string" ? new Error(bookRead.error) : bookRead.error;
      quote = taker.novigQuoteOf(bookRead.book, contract);
    } catch (error) {
      return { answer: { ok: false, ...position, reason: "no_book", detail: error.message }, catalog: timing };
    }
    if (mismatch) return { answer: { ok: false, ...position, reason: "market_mismatch", detail: mismatch }, catalog: timing };
    if (!quote) return { answer: { ok: false, ...position, reason: "no_book", detail: `nothing rests on the other side of Novig market ${contract.marketId}` }, catalog: timing };
    // The taker fee per contract the bets service adds to the price (0 pregame; novigFeeRateOf).
    const feeRate = taker.novigFeeRateOf(checkedMarket, spec.live);
    return { answer: { ok: true, ...position, price: quote.ask, feeRate, depth: quote.depth, novigBaseUrl, line: taker.novigLineOf(row) }, catalog: timing };
  }

  // The Ticket's Place button at Polymarket US (bets_service/polymarket_place.py
  // asks): the board's one Polymarket US line with these fields, the contract
  // it buys, found in Polymarket US's own list (polymarket.js: the game by its
  // two teams through the crosswalk, the market by its own sides), what buying
  // it costs now and the line as a bet record names it. Places nothing.
  //   params  URLSearchParams {eventId, periodTypeId, betTypeId, sideIndex, points?}
  async function polymarketQuoteNow(params) {
    const spec = kalshiLineSpecOf(params);
    const startedAt = now();
    const answer = await polymarketQuoteOfSpec(spec);
    logInfo(`polymarket quote ${spec.eventId} ${answer.label || `bet type ${spec.betTypeId} side ${spec.sideIndex} ${spec.points}`}: `
      + `${answer.ok ? `ok ${answer.slug} ${answer.side} @ ${answer.ask} (book from ${answer.bookFrom})` : `${answer.reason}: ${answer.detail}`} in ${now() - startedAt} ms`);
    return answer;
  }

  async function polymarketQuoteOfSpec(spec) {
    const { board, reason: boardReason, detail: boardDetail } = quoteBoardOf(spec);
    if (!board) return { ok: false, reason: boardReason, detail: boardDetail };
    const found = linesOfEvent(board, spec.eventId).filter((line) => polymarket.isPolymarketUsBook(board.books[line.bookId])
      && line.eventId === spec.eventId && line.periodTypeId === spec.periodTypeId && line.betTypeId === spec.betTypeId
      && line.sideIndex === spec.sideIndex && (spec.betTypeId === MONEYLINE_BET_TYPE_ID || (line.points ?? null) === spec.points));
    if (found.length !== 1) {
      return { ok: false, reason: "no_line", detail: `the board holds ${found.length} Polymarket US lines for game ${spec.eventId} at that market, side and number` };
    }
    const row = feedLib.describeLine(found[0], board);
    const fields = { eventId: row.eventId, eventStart: row.eventStart, label: `${row.sideLabel} ${row.period}`.trim(), eventName: row.eventName };
    const foundContract = await polymarketContractOfRow(row, readPolymarketEvents);
    if (foundContract.reason) return { ok: false, ...fields, reason: foundContract.reason, detail: foundContract.detail };
    const { contract, event, startMs } = foundContract;
    // Live, Polymarket US must itself show the game in play (its event's `live`, recon 2026-11-06).
    if (spec.live && event.live !== true) {
      return { ok: false, ...fields, reason: "not_live", detail: `Polymarket US shows ${event.slug} live=${JSON.stringify(event.live)}, ended=${JSON.stringify(event.ended)}` };
    }
    if (foundContract.learnedRow) {
      const { learnedRow } = foundContract;
      logInfo(`polymarket quote ${row.eventId}: matched ${event.slug} by one team; "${learnedRow.venueTeamName}" is Unabated team ${learnedRow.unabatedTeamId}`);
      postCrosswalkRows([learnedRow]).catch((error) => logWarning(`crosswalk: polymarket name not saved: ${error.message}`));
    }
    const eventStartMs = Math.min(startMs, Date.parse(event.startTime));
    const position = { ...fields, eventStart: new Date(eventStartMs).toISOString(), unabatedStart: row.eventStart,
      polymarketStart: event.startTime, slug: contract.slug, side: contract.side, eventSlug: event.slug, live: spec.live };
    const read = await polymarketPlaceBook(contract.slug);
    if (read.error) return { ok: false, ...position, reason: "no_book", detail: read.error };
    const { book } = read;
    const quote = polymarket.quoteOf(book, contract.side);
    if (quote.reason) return { ok: false, ...position, reason: quote.reason, detail: `${contract.slug}: ${quote.detail}` };
    return {
      ok: true, ...position, ask: quote.ask, allIn: polymarket.allInOf(quote.ask), depth: quote.depth, bookFrom: read.from,
      feeTheta: polymarket.TAKER_FEE_THETA, polymarketGatewayUrl, line: taker.novigLineOf(row),
    };
  }

  // The Ticket's Place button at ProphetX (bets_service/prophetx_place.py
  // asks): the board's one ProphetX line with these fields, the line as a bet
  // record names it, and the ProphetX strike ids Unabated carries for its game
  // (prophetx.js). ProphetX's API needs the bets service's key for every read,
  // so the strike and its price are found there. Places nothing.
  //   params  URLSearchParams {eventId, periodTypeId, betTypeId, sideIndex, points?}
  function prophetxQuoteNow(params) {
    const spec = kalshiLineSpecOf(params);
    const answer = prophetxQuoteOfSpec(spec);
    logInfo(`prophetx quote ${spec.eventId} ${answer.label || `bet type ${spec.betTypeId} side ${spec.sideIndex} ${spec.points}`}: `
      + `${answer.ok ? `ok, ${answer.strikeIds.length} strike ids${answer.lineStrikeId ? " (the line's own)" : ""}` : `${answer.reason}: ${answer.detail}`}`);
    return answer;
  }

  function prophetxQuoteOfSpec(spec) {
    const { board, reason: boardReason, detail: boardDetail } = quoteBoardOf(spec);
    if (!board) return { ok: false, reason: boardReason, detail: boardDetail };
    const gameLines = linesOfEvent(board, spec.eventId);
    const found = gameLines.filter((line) => line.bookId === prophetx.PROPHETX_BOOK_ID
      && line.eventId === spec.eventId && line.periodTypeId === spec.periodTypeId && line.betTypeId === spec.betTypeId
      && line.sideIndex === spec.sideIndex && (spec.betTypeId === MONEYLINE_BET_TYPE_ID || (line.points ?? null) === spec.points));
    if (found.length !== 1) {
      return { ok: false, reason: "no_line", detail: `the board holds ${found.length} ProphetX lines for game ${spec.eventId} at that market, side and number` };
    }
    const row = feedLib.describeLine(found[0], board);
    return {
      ok: true, eventId: row.eventId, eventStart: row.eventStart, label: `${row.sideLabel} ${row.period}`.trim(),
      eventName: row.eventName, price: found[0].price, line: taker.novigLineOf(row), live: spec.live, ...prophetx.strikeIdsOf(found[0], gameLines),
    };
  }

  // A venue's two team names as Unabated team keys ("<league>:<team id>"),
  // through the same resolver every bet record goes through (teams.js + the
  // held crosswalk, bets.resolveTeamKeys), so a Place check reads a team the
  // way bet matching does. ProphetX Place asks (bets_service/prophetx_place.py):
  // its game's names, to compare with the board line's team ids. A key is
  // null when nothing resolves the name. Reads nothing.
  //   params  URLSearchParams {venue, league, away, home}
  function teamKeysNow(params) {
    const text = (name) => {
      const value = params.get(name);
      if (typeof value !== "string" || value.trim() === "") throw httpError(400, `${name} must be a non-empty string, got ${JSON.stringify(value)}`);
      return value;
    };
    const record = { venue: text("venue"), league: text("league"), awayTeam: text("away"), homeTeam: text("home") };
    const resolved = betsLib.resolveTeamKeys([record], held.crosswalk)[0];
    return { ok: true, ...record, awayKey: resolved.awayKey ?? null, homeKey: resolved.homeKey ?? null };
  }

  // The Ticket's Place button at BFA (bets_service/bfa_place.py asks): the
  // board's one Buckeye line with these fields, its side's rotation and the
  // line as a bet record names it. BFA is read by the bets service on its own
  // login, so the game and the live price are found there. Places nothing.
  //   params  URLSearchParams {eventId, periodTypeId, betTypeId, sideIndex, points?}
  function bfaQuoteNow(params) {
    const spec = kalshiLineSpecOf(params);
    const answer = bfaQuoteOfSpec(spec);
    logInfo(`bfa quote ${spec.eventId} ${answer.label || `bet type ${spec.betTypeId} side ${spec.sideIndex} ${spec.points}`}: `
      + `${answer.ok ? `ok, rotation ${answer.rotation} at ${answer.price}` : `${answer.reason}: ${answer.detail}`}`);
    return answer;
  }

  function bfaQuoteOfSpec(spec) {
    const board = feedState;
    if (!board) return { ok: false, reason: "no_board", detail: "the board has not loaded yet" };
    const found = linesOfEvent(board, spec.eventId).filter((line) => line.bookId === BUCKEYE_BOOK_ID
      && line.eventId === spec.eventId && line.periodTypeId === spec.periodTypeId && line.betTypeId === spec.betTypeId
      && line.sideIndex === spec.sideIndex && (spec.betTypeId === MONEYLINE_BET_TYPE_ID || (line.points ?? null) === spec.points));
    if (found.length !== 1) {
      return { ok: false, reason: "no_line", detail: `the board holds ${found.length} BFA lines for game ${spec.eventId} at that market, side and number` };
    }
    const row = feedLib.describeLine(found[0], board);
    return {
      ok: true, eventId: row.eventId, eventStart: row.eventStart, label: `${row.sideLabel} ${row.period}`.trim(),
      eventName: row.eventName, price: found[0].price, isAlt: row.isAlt, rotation: row.rotation ?? null,
      awayRotation: row.awayRotation ?? null, homeRotation: row.homeRotation ?? null, line: taker.novigLineOf(row),
    };
  }

  // The Polymarket US contract a board row buys, found in Polymarket US's own
  // list (polymarket.js: the game by its two teams through the crosswalk, the
  // market by its own sides): {contract, event, startMs, learnedRow} or
  // {reason, detail}. Shared by Place's quote (which saves a learnedRow, a
  // name matched by one team) and the live price pass (which saves nothing).
  //   readEvents  (marketType, startMs) -> Promise of Polymarket US's events around the start (throws naming why)
  async function polymarketContractOfRow(row, readEvents) {
    const type = polymarket.marketTypeOf(row);
    if (type.reason) return { reason: type.reason, detail: type.detail };
    const startMs = Date.parse(row.eventStart);
    if (!Number.isFinite(startMs)) return { reason: "no_start", detail: "the board has no start for this game" };
    if (row.awayTeamId == null || row.homeTeamId == null) return { reason: "no_teams", detail: "the board names no team ids for this game" };
    let events;
    try {
      events = await readEvents(type.type, startMs);
    } catch (error) {
      return { reason: "no_market", detail: `Polymarket US's events unread: ${error.message}` };
    }
    const keysOf = (record) => betsLib.resolveTeamKeys([record], held.crosswalk)[0];
    const boardKeys = { away: teams.keyOf(row.league, row.awayTeamId), home: teams.keyOf(row.league, row.homeTeamId) };
    const exact = polymarket.matchingEvents(events, row, keysOf, boardKeys);
    // No event names both teams: a college name nothing resolves yet may still
    // be the game by its other team (polymarket.oneTeamMatches).
    const games = exact.length ? exact : polymarket.oneTeamMatches(events, row, keysOf, boardKeys);
    if (games.length !== 1) {
      const byOneTeam = exact.length ? "" : ` (${games.length} by one team and start)`;
      return { reason: "no_market", detail: `${exact.length} of Polymarket US's ${events.length} ${type.type} events are ${row.eventName} at its start${byOneTeam}` };
    }
    const { event, teams: eventTeams, learnedRow } = games[0];
    const contract = polymarket.contractOf(row, event, eventTeams, type.type);
    if (contract.reason) return { reason: contract.reason, detail: contract.detail };
    return { contract, event, startMs, learnedRow: learnedRow || null };
  }

  // Polymarket US's active events of one market type starting within
  // START_TOLERANCE_MS of startMs (public gateway). Throws naming why; a full
  // page may be a cut list, so it throws rather than miss the game.
  async function readPolymarketEvents(marketType, startMs) {
    const query = new URLSearchParams({
      limit: String(polymarket.EVENTS_PAGE_LIMIT), active: "true", sportsMarketTypes: marketType,
      startTimeMin: new Date(startMs - polymarket.START_TOLERANCE_MS).toISOString(),
      startTimeMax: new Date(startMs + polymarket.START_TOLERANCE_MS).toISOString(),
    });
    const body = await getServiceJson(polymarketFetch, `${polymarketGatewayUrl}/v1/events?${query}`, { headers: { Accept: "application/json" } });
    if (!body || !Array.isArray(body.events)) throw new Error("Polymarket US events: no events list");
    if (body.events.length >= polymarket.EVENTS_PAGE_LIMIT) throw new Error(`${body.events.length} events in one page may be a cut list`);
    return body.events;
  }

  // readPolymarketEvents for the live price pass, each query's answer (or
  // failure) reused for POLYMARKET_EVENTS_REUSE_MS: one game's ladder is one
  // read, not one per rung. Place reads its own, fresh.
  function polymarketEventListHeld(marketType, startMs) {
    const held_ = polymarketEventLists.get(`${marketType}|${startMs}`);
    return Boolean(held_ && now() - held_.readAt <= POLYMARKET_EVENTS_REUSE_MS);
  }

  async function readPolymarketEventsReused(marketType, startMs) {
    const key = `${marketType}|${startMs}`;
    const held_ = polymarketEventLists.get(key);
    if (held_ && now() - held_.readAt <= POLYMARKET_EVENTS_REUSE_MS) {
      if (held_.error) throw new Error(held_.error);
      return held_.events;
    }
    try {
      const events = await readPolymarketEvents(marketType, startMs);
      polymarketEventLists.set(key, { events, readAt: now() });
      return events;
    } catch (error) {
      polymarketEventLists.set(key, { error: error.message, readAt: now() });
      throw error;
    }
  }

  // Place's book: the stream's (subscribed first; a market it had not watched
  // gets POLYMARKET_STREAM_WAIT_MS for its first book), else ONE read of the
  // public gateway, which the live pass no longer reads books off.
  //   -> {book, from: "stream" | "gateway"} | {error}
  async function polymarketPlaceBook(slug) {
    const streamed = (await readPolymarketStreamBooks([slug], { priority: true })).get(slug);
    if (streamed.book) return { book: streamed.book, from: "stream" };
    logInfo(`polymarket quote ${slug}: no stream book (${streamed.error}), reading the gateway`);
    try {
      return { book: await readPolymarketBook(slug), from: "gateway" };
    } catch (error) {
      return { error: `${error.message} (stream: ${streamed.error})` };
    }
  }

  async function readPolymarketBook(slug) {
    return getServiceJson(polymarketFetch, `${polymarketGatewayUrl}/v1/markets/${encodeURIComponent(slug)}/book`, { headers: { Accept: "application/json" } });
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

  function novigMarketPath(marketId) {
    return `${NOVIG_CATALOG_MARKETS_PATH}/${encodeURIComponent(marketId)}`;
  }

  // Why Novig's market and its event, read now, are not the row's or do not
  // take a free pregame order (taker.novigMarketMismatch), an unread market
  // included (fail closed); null when they are. Shared by the bot's orders
  // and the Place button's quotes, so neither can skip a check the other has.
  //   listedEventId  the event the catalog listed the market under, when known:
  //                  then the market and the event are read at once
  async function readNovigMarketMismatch(row, contract, listedEventId = null) {
    return (await readNovigMarketCheck(row, contract, listedEventId, false)).mismatch;
  }

  // {mismatch: why Novig's market is not the row's (null when it is), market: the market read, or null}.
  //   live  a live Place: the event must be in game and the fee is priced (taker.novigMarketMismatch)
  async function readNovigMarketCheck(row, contract, listedEventId, live) {
    const eventPath = (eventId) => `/v3/catalog/events/${encodeURIComponent(eventId)}`;
    let market;
    let event;
    try {
      if (typeof listedEventId === "string" && listedEventId) {
        [market, event] = await Promise.all([getNovigJson(novigMarketPath(contract.marketId)), getNovigJson(eventPath(listedEventId))]);
      } else {
        market = await getNovigJson(novigMarketPath(contract.marketId));
        event = await getNovigJson(eventPath(market.eventId));
      }
    } catch (error) {
      return { mismatch: `Novig market ${contract.marketId} unread: ${error.message}`, market: null };
    }
    const mismatch = taker.novigMarketMismatch(market, event, row, contract, { live });
    return { mismatch: mismatch === null ? null : `Novig market ${contract.marketId}: ${mismatch}`, market };
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
      // The market's Kalshi cash pot: the desk tops it up from the bot's other pot before ordering (null = Kalshi omits it).
      decision.exchangeIndex = Number.isInteger(market.exchange_index) ? market.exchange_index : null;
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
      heldIds: decision.heldIds, decidedAt: decision.decidedAt, kalshiBaseUrl, exchangeIndex: decision.exchangeIndex ?? null,
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
      timers.setInterval(() => postLadders(), CLOSING_FAIRS_POST_MS),
      timers.setInterval(() => postEdgeSegments(), CLOSING_FAIRS_POST_MS),
      // The live price pass right after each taker pass, so the lines the taker just read are not read again.
      timers.setInterval(async () => {
        await takerPass();
        await livePricePass();
      }, TAKER_PASS_MS),
      timers.setInterval(() => postTauSamples(), TAU_SAMPLES_POST_MS),
    ];
    await loadTauPool();
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
      dismissedIds, edgeSkips, lineView: liveLineView(),
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

  // The Live tab's cards, off the remembered games (scenarios.js). With the
  // Live tab's filter in `params`, `games` holds the cards rebuilt from only
  // the bets it keeps; coveredBetIds / finishedBetIds stay every bet's, since
  // the tracker splits its Open tab and badges off them.
  function scenariosPayloadNow(params) {
    const filter = params ? scenarios.liveFilterOf(params) : null;
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
      return { row: entry.row, ladderOf, oddsAt: entry.oddsAt, gameState: entry.gameState };
    });
    const records = held.records.filter((record) => !excludedBetIds.has(record.id));
    const everyBet = scenarios.buildScenarios({ records, games, now: at });
    return {
      generatedAt: new Date(at).toISOString(),
      scanner: scannerStatus ? { phase: scannerStatus.phase, error: scannerStatus.error } : { phase: "starting", error: null },
      betsService: { okAt: betsStatus.okAt, error: betsStatus.error },
      ...everyBet,
      games: filter ? scenarios.buildScenarios({ records: scenarios.recordsInFilter(records, filter), games, now: at }).games : everyBet.games,
      filtered: filter !== null,
    };
  }

  function health() {
    return {
      ok: true, generatedAt: new Date(now()).toISOString(), uptimeSec: Math.round((now() - startedAt) / 1000),
      scanner: scannerStatus ? { phase: scannerStatus.phase, error: scannerStatus.error, lastSnapshotAt: scannerStatus.lastSnapshotAt } : { phase: "starting", error: null },
      betsService: { okAt: betsStatus.okAt, error: betsStatus.error },
      settings: { source: settings.source, okAt: settings.okAt, error: settings.error },
      closingFairs: { ...closingStatus },
      closeSnapshots: { ...closeSnapshotStatus },
      ladders: { ...ladderStatus, tracked: sentLadders.size },
      edgeLog: { ...edgeLogStatus, openSegments: edgeSegments.size, unsent: unsentEdgeRows.size },
      tauPool: { ...tauStatus, size: tauPool.size(), unsent: tauUnsent.size },
      livePrices: liveStatus(),
      taker: {
        desk: takerDeskStatus, novigDesk: novigDeskStatus, lastPassAt: isoOrNull(takerLastPass.finishedAt),
        skipped: takerLastPass.skipped, novigSkipped: takerLastPass.novigSkipped,
        rows: takerLastPass.rows, byReason: takerLastPass.byReason, log: { okAt: takerLogStatus.okAt, error: takerLogStatus.error, dropped: takerLogStatus.dropped, unsent: takerUnsent.length },
      },
    };
  }

  return {
    start, stop, pollBets, pollSettings, postClosingFairs, postLadders, postEdgeSegments, postTauSamples, scanLoaded: () => scanLoad, learn, takerPass, postTakerDecisions, livePricePass,
    edgesPayload: edgesPayloadNow, livePricesPayload: livePricesPayloadNow, teasersPayload: teasersPayloadNow, attachPayload: attachPayloadNow, tauPayload: tauPayloadNow,
    scenariosPayload: scenariosPayloadNow, takerPayload: takerPayloadNow, kalshiQuote: kalshiQuoteNow, novigQuote: novigQuoteNow, polymarketQuote: polymarketQuoteNow, prophetxQuote: prophetxQuoteNow, bfaQuote: bfaQuoteNow, teamKeys: teamKeysNow, health, scanner,
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
  const liveText = params.get("live");
  if (liveText !== null && liveText !== "0" && liveText !== "1") throw httpError(400, `live must be 0 or 1, got ${JSON.stringify(liveText)}`);
  return {
    eventId: whole("eventId"), periodTypeId: whole("periodTypeId"), betTypeId: whole("betTypeId", [1, 2, 3]),
    sideIndex: whole("sideIndex", [0, 1]), points, live: liveText === "1",
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

// The HTTP side: GET /edges.json, /teasers.json, /attach.json, /scenarios.json, /taker.json, /tau.json, /live_prices.json, /kalshi_quote.json, /novig_quote.json, /polymarket_quote.json, /prophetx_quote.json, /bfa_quote.json, /team_keys.json and /health, the Host allowlist first on
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
      if (path === "/scenarios.json") return sendJson(response, 200, await runner.scenariosPayload(url.searchParams));
      if (path === "/taker.json") return sendJson(response, 200, await runner.takerPayload());
      if (path === "/tau.json") return sendJson(response, 200, await runner.tauPayload());
      if (path === "/live_prices.json") return sendJson(response, 200, await runner.livePricesPayload());
      if (path === "/kalshi_quote.json") return sendJson(response, 200, await runner.kalshiQuote(url.searchParams));
      if (path === "/novig_quote.json") return sendJson(response, 200, await runner.novigQuote(url.searchParams));
      if (path === "/polymarket_quote.json") return sendJson(response, 200, await runner.polymarketQuote(url.searchParams));
      if (path === "/prophetx_quote.json") return sendJson(response, 200, await runner.prophetxQuote(url.searchParams));
      if (path === "/bfa_quote.json") return sendJson(response, 200, await runner.bfaQuote(url.searchParams));
      if (path === "/team_keys.json") return sendJson(response, 200, await runner.teamKeys(url.searchParams));
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
