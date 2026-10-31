"""Entry point: the local bets service the Unabated Ticket panel polls (#114).

    python -m unabated_ticket.bets_service.service      (or bets_service/run.sh)

Inputs:  each registered Source (sources/kalshi.py; sources/betonline.py when
         its cookie file exists; sources/novig.py when its token file exists;
         sources/bfa.py and sources/wagerzon.py when their logins are configured;
         sources/polymarket_us.py when its API key is configured)
         on its own poll_sec; plus the extension's pushes for Bet105 and DraftKings
         (below).
Outputs: HTTP on 127.0.0.1:8094 (loopback only, no auth):
           GET /bets.json[?days=N]  {generatedAt, sources: {name: {fetchedAt, ok,
                                    error, count}}, bets: [records open + settled
                                    within N days (default RETENTION_DAYS=30)],
                                    crosswalk: [team_crosswalk rows, newest first],
                                    pins: [bet_pins rows, newest first],
                                    fillFairs: [bet_fill_fairs rows of those bets],
                                    closingFairs: [bet_closing_fairs rows of those bets],
                                    exclusions: [bet_exclusions rows, newest first],
                                    dismissals: [{betId, dismissedAt}] (open bets only),
                                    teaserBlocks: [{marketKey, eventStartMs, blockedAt}]
                                    (games still to start),
                                    edgeSkips: [{skipKey, eventStartMs, skippedAt}]
                                    (games still to start)}
           GET /bets.json?snapshots=1  the same plus fillSnapshots and closeSnapshots:
                                    [bet_fill_snapshots / bet_close_snapshots rows of those
                                    bets] (every book's two sides at each bet's number at the
                                    fill and at the close; opt-in, they are large)
           GET /health              {ok, generatedAt, uptimeSec, sources}
           POST /crosswalk.json     body {rows: [{venue, league, venueTeamKey,
                                    unabatedTeamId, venueTeamName?, unabatedTeamName?,
                                    learnedFrom?}]} -> {ok, learned, conflicts, crosswalk}
                                    (#118 step 4: the panel's lessons from id joins;
                                    Content-Type must be application/json — a web page
                                    cannot send that cross-origin without a preflight
                                    this server never answers, so no site can write here)
           DELETE /crosswalk.json   -> {ok, cleared, crosswalk: []}
           POST /pins.json          body {pin: {betId, venue, league, eventId, eventStart?,
                                    awayTeamId?, homeTeamId?, awayTeamName?, homeTeamName?},
                                    crosswalk: [at most 2 rows, the crosswalk shape]}
                                    -> {ok, pins, crosswalk}; 404 when no bet has that id,
                                    400 when the pin's venue is not the bet's
                                    (Cal's manual attach: the pin plus the team names it
                                    teaches, which REPLACE a held key)
           DELETE /pins.json?betId= -> {ok, removedPin, removedRows, pins, crosswalk}
                                    (Undo: the pin and the rows it taught)
           POST /fill_fairs.json    body {rows: [{betId, lineKey, points, fairAmerican,
                                    fairObservedAt, placedAt}]} -> {ok, saved, fillFairs}
                                    (the fair a new bet's line had when it was placed,
                                    read by the panel off its line history; same
                                    Content-Type guard as the crosswalk)
           POST /fill_books.json    body {rows: [{betId, placedAt, fillPrice, bookCount,
                                    booksBetter, bestBook, bestPrice}]} -> {ok, saved}
                                    (the other books at a new bet's number when it was
                                    placed, server/fillbooks.js; insert-only; same
                                    Content-Type guard)
           POST /fill_snapshots.json body {rows: [{betId, placedAt, lineKey, points, mainPoints,
                                    mainFairAmerican, pointsFromMain, edgeMoveLineKey, edgeMove,
                                    edgeMoveFairDelta, edgeMovePriceDelta, books}]} -> {ok, saved}
                                    (a new bet's market at the fill, server/edgesnapshot.js;
                                    insert-only; same Content-Type guard)
           POST /close_snapshots.json body {rows: [{betId, observedAt, eventStart, books}]}
                                    -> {ok, saved} (every book's two sides at each open bet's
                                    number before its game starts, server/edgesnapshot.js; a
                                    row older than the stored one is ignored; same guard)
           POST /closing_fairs.json body {rows: [{betId, lineKey, points, fairAmerican,
                                    fairObservedAt, eventStart}]} -> {ok, saved}
                                    (the server runner's pre-start read of each open bet's
                                    line, server/closefair.js; a row older than the stored
                                    one is ignored; same Content-Type guard)
           POST /tau_samples.json   body {rows: [{key, marketKey, edgeUnabated, edgeExchanges,
                                    seenAt}]} -> {ok, saved}
                                    (the server runner's tail-flex tau samples, issue #208:
                                    newest seenAt wins, rows past 25 h pruned; same
                                    Content-Type guard)
           GET /tau_samples.json    {generatedAt, rows} — the samples of the last 25 h,
                                    the runner's pool on a restart
           POST /exclusions.json    body {betIds: [...], excluded: bool} -> {ok, changed,
                                    exclusions} (the Bet Tracker's Bets page removing or
                                    restoring bets that are not Cal's; BFA and Wagerzon
                                    bets only — anything else is a 400, an unknown id a
                                    404; same Content-Type guard)
           POST /dismissals.json    body {betIds: [...], dismissed: bool} -> {ok, changed,
                                    dismissals} (the panel's and the phone's Dismiss /
                                    Restore on an open bet no board game matches: stop or
                                    resume flagging it; an unknown id is a 404; same
                                    Content-Type guard)
           POST /teaser_blocks.json body {marketKey, eventStartMs, blocked: bool} -> {ok,
                                    changed, teaserBlocks} (the Teasers tab's Can't tease /
                                    Restore on a spread or total; same Content-Type guard)
           POST /edge_skips.json    body {skipKey, eventStartMs, skipped: bool} -> {ok,
                                    changed, edgeSkips} (the Edges Skip / Restore on a
                                    line a book beats; same Content-Type guard)
           POST /bet105.json        body {fetchedAt, feeds: {prematch: [betGroup], live:
                                    [betGroup]}, settled: [wager]} -> {ok, count, settled,
                                    closed, skipped}, or {error} -> {ok, recorded: "error"}
                                    (a PUSHED source: the extension reads Bet105 from
                                    Cal's own Chrome because Cloudflare challenges anything
                                    else — sources/bet105.py parses; a graded wager settles
                                    its bet with the venue's result and settle time, an open
                                    bet neither read lists is closed with no result, and the
                                    push is logged as that source's run)
           POST /draftkings.json    body {fetchedAt, open: [bet], settled: [bet], events:
                                    {eventId: event}} -> {ok, count, closed, skipped}, or
                                    {error} -> {ok, recorded: "error"} (pushed like Bet105:
                                    the extension reads DraftKings' My Bets socket from Cal's
                                    own Chrome because Akamai and a login code stop anything
                                    else — sources/draftkings.py parses; a settled bet
                                    carries the venue's result, an open bet neither list
                                    carries is closed with no result)
           GET /settings.json       {settings: {bankroll, multiplier, leagues, periods, betTypes,
                                    bookMode, bookIds, minEdgePct, minStake, maxLineAgeHours,
                                    minLiquidityToWin, includeAlts, sortBy, groupByMarket,
                                    minShrunkEdgePct},
                                    updatedAt} — the Edges settings the server
                                    runner (unabated_ticket/server/runner.js) reads every cycle;
                                    null = the panel's default (extension/edgerows.js)
           POST /place_teaser.json  body {stake, legs: [4 x {league, betType, side, rotation, points,
                                    eventStart, label}]} -> {ok, status: placed | refused |
                                    unconfirmed, message, ticket?} (the Teasers tab's Place
                                    button: ONE 4-team 6-point teaser at BFA, real money —
                                    bfa_teaser.py checks every number on BFA's board first and
                                    confirms off the open bets; 503 when no BFA account is read)
           POST /place_kalshi.json  body {clickId, eventId, periodTypeId, betTypeId, sideIndex, points,
                                    stake, ticketAllIn, kalshiReadAt} -> {ok, status: filled |
                                    partial | unfilled | refused | rejected | unknown, message, ...}
                                    (the Ticket's Place button: ONE immediate-or-cancel Kalshi
                                    order on Cal's own key, real money — manual_orders.py checks
                                    the live ask against the Ticket's price first; 503 when no
                                    Kalshi key is configured)
           POST /place_novig.json   body {clickId, eventId, periodTypeId, betTypeId (2|3),
                                    sideIndex, points, stake, ticketAllIn, novigReadAt,
                                    kalshiReadAt} -> {ok, status: filled | partial | unfilled |
                                    refused | rejected | unknown, message, ...} (the Ticket's Place
                                    button on Novig: ONE immediate-or-cancel order from Place's own
                                    Novig sub-account, real money — novig_place.py checks the live
                                    price against the Ticket's first; 503 while Place has no Novig
                                    order client)
           POST /place_polymarket.json  body {clickId, eventId, periodTypeId, betTypeId, sideIndex,
                                    points, stake, ticketAllIn, polymarketReadAt, kalshiReadAt,
                                    novigReadAt} -> {ok, status, message, ...} (the Ticket's Place
                                    button at Polymarket US: ONE immediate-or-cancel order on Cal's
                                    main account with the POLYMARKET_US_* key, real money —
                                    polymarket_place.py checks the live all-in price against the
                                    Ticket's first; 503 while that key is unset or refused)
           POST /novig_books.json   body {marketIds: [..], waitMs, priority?} -> {stream:
                                    {connected, lastFrameAgeMs, error, ...}, books: {marketId:
                                    {status, started, seq, orders: {outcomeId: [{order, price,
                                    qty}]}} | {error}}}: Novig's live books off the bot key's
                                    websocket (novig_book_stream.py; reads only), each market
                                    subscribed on first ask (priority: ahead of the rest, the
                                    Place quote), waiting up to waitMs (capped at
                                    novig_book_stream.MAX_WAIT_SEC) for a new market's snapshot.
                                    Loopback only (the runner); 503 {error} with no bot Novig key
           POST /novig_read.json    body {path, query?} -> {status, body}: one GET of Novig's
                                    keyed catalog (/v3/catalog/markets[?league,eventStatus,
                                    limit,after], /markets/<id>, /markets/<id>/book,
                                    /events/<id>) signed with the bot's Novig key, Novig's
                                    status and body passed back as they came (a 429 too), so
                                    the runner's checks never touch the per-IP public throttle
                                    (novig_catalog_reads.py). Loopback only (the runner); 400
                                    on any other path or query, 502 when Novig did not answer,
                                    503 {error} with no bot Novig key
           PUT /settings.json       body {settings: {<any of those fields>: value or null}} ->
                                    {ok, settings, updatedAt}; a field left out keeps its value,
                                    null resets it to the default, an unknown field or a bad
                                    value is a 400 naming it (same Content-Type guard as the
                                    POSTs: a web page cannot send JSON cross-origin)
           GET /balances.json       {generatedAt, balances: {venue: {cash, inBets, available,
                                    creditLimit, fetchedAt, error}}}, the live per-venue
                                    account balance read after each poll (balances.py;
                                    memory only, a failed read keeps the last numbers)
           GET /edges.json          the server runner's Edges list (config.RUNNER_URL,
                                    server/edges_payload.js documents it), passed through
                                    as it came; 502 {error, runnerUrl} when the runner
                                    cannot be reached within config.RUNNER_TIMEOUT_SEC or
                                    answers anything but 200 (phone page plan step 2:
                                    one origin for the page, its reads and its PUT)
           GET /scenarios.json      the server runner's Live-tab scenarios (server/scenarios.js
                                    documents them), passed through the same way, with the
                                    same 502 (the Bet Tracker's Live tab)
           GET /teasers.json        the server runner's Teasers list (server/teasers_payload.js
                                    documents it), passed through the same way (the phone's
                                    Teasers tab)
           GET /attach.json?...     the server runner's Attach picker (step 1: games for a
                                    bet and a search; step 2: the plan and the /pins.json
                                    body), query string passed on; the runner's own 400/404
                                    (no such open bet, a game gone) passes through as it came
           GET / and the phone page's files   STATIC_FILES, a fixed map of URL path ->
                                    file (server/phone/ and the extension's pure modules
                                    the page loads under /ext/; the Bet Tracker,
                                    server/tracker/, under /tracker/); any other path is a 404,
                                    so no request can name a file outside it
         Every verb refuses a request whose Host header is not the loopback
         name the service is serving on, or a name in BETS_EXTRA_ALLOWED_HOSTS
         (the tailnet name `tailscale serve` forwards; bare or :443), with
         403: a page at evil.example whose DNS flips to 127.0.0.1 is
         SAME-ORIGIN with this server, so neither CORS nor the JSON
         Content-Type guard applies to it.
Side effects: UPSERTs records into bets.duckdb::bets and APPENDs a row to
bets.duckdb::source_runs per poll (see store.py); INSERTs into / DELETEs from
bets.duckdb::team_crosswalk on the two crosswalk routes; UPSERTs / DELETEs
bets.duckdb::bet_pins and the crosswalk rows a pin taught on the two pin
routes; INSERTs into bets.duckdb::bet_fill_fairs on POST /fill_fairs.json —
insert-only, a bet that has a saved fair keeps it; INSERTs into
bets.duckdb::bet_fill_books on POST /fill_books.json (insert-only); INSERTs into
bets.duckdb::bet_fill_snapshots on POST /fill_snapshots.json (insert-only); UPSERTs
bets.duckdb::bet_closing_fairs on POST /closing_fairs.json and
bets.duckdb::bet_close_snapshots on POST /close_snapshots.json (newest observation
wins); UPSERTs / prunes bets.duckdb::tau_samples on POST /tau_samples.json; INSERTs into / DELETEs from bets.duckdb::bet_exclusions on POST
/exclusions.json; INSERTs into / DELETEs from bets.duckdb::bet_dismissals on
POST /dismissals.json; UPSERTs / DELETEs bets.duckdb::teaser_blocks on POST
/teaser_blocks.json and bets.duckdb::edge_skips on POST /edge_skips.json;
rotating log at
bets_service.log. A poll that raises writes a failed source_runs row and
leaves `bets` untouched — a dark source never blanks the list. PUT
/settings.json UPSERTs the one row of bets.duckdb::edge_settings. POST
/place_teaser.json places a wager at BFA (money, not disk) and UPSERTs the BFA open
bets it read while confirming into bets.duckdb::bets, so the panel's next poll has the
ticket. POST /taker_decisions.json APPENDs the server runner's taker decisions to
taker.duckdb::taker_decisions (taker_store.py; a separate file from bets.duckdb).
POST /place_kalshi.json places Cal's own Kalshi order (money) and INSERTs /
UPDATEs bets.duckdb::manual_orders (manual_orders.py).
POST /place_novig.json places a Novig order from Place's own sub-account (money)
and INSERTs / UPDATEs bets.duckdb::manual_orders, venue 'novig' (novig_place.py).
POST /place_polymarket.json places a Polymarket US order on Cal's main account
(money) and INSERTs / UPDATEs bets.duckdb::manual_orders, venue 'polymarket_us'
(polymarket_place.py).
POST /taker/order.json places a bot order on Kalshi (money) on the bot's own
key, locked to the bot's own Kalshi subaccount, when one is configured, and
INSERTs / UPDATEs taker.duckdb::taker_orders (taker_orders.py: stop file, open-cost
ceiling at 20% of the Kelly bankroll, limit at the ask, rest cancelled); the desk's own thread cancels bot
orders left resting on Kalshi, reads final fills and frees settled markets
every 30 s, runner or no runner.
POST /taker/novig_order.json places a bot order on Novig (money) on the bot's
own trading key, which trades only in the bot's Novig subaccount, and INSERTs /
UPDATEs taker.duckdb::novig_taker_orders (novig_taker_desk.py: stop files,
holdings freshness, its own $500 ceiling, one IOC buy); its thread reads final
fills and frees settled markets every 30 s. The bot's (and the Novig Place
button's) subaccount is read as its own source, novig_bot / novig_place
(sources/novig_subaccount.py).
"""
import json
import logging
import math
import re
import signal
import threading
import time
import urllib.error
import urllib.request
from collections.abc import Callable
from datetime import datetime, timezone
from http.server import BaseHTTPRequestHandler, ThreadingHTTPServer
from pathlib import Path
from urllib.parse import parse_qs, urlparse

from unabated_ticket.bets_service import bfa_teaser, config
from unabated_ticket.bets_service.balances import BalanceBook, balances_payload, read_balance_once
from unabated_ticket.bets_service.log_setup import setup_logging
from unabated_ticket.bets_service.manual_orders import ManualKalshiDesk, read_runner_quote
from unabated_ticket.bets_service.manual_orders import validate_place_request as validate_kalshi_place_request
from unabated_ticket.bets_service.polymarket_orders import PolymarketError, PolymarketOrderClient
from unabated_ticket.bets_service.polymarket_place import (RUNNER_POLYMARKET_QUOTE_PATH,
                                                           RUNNER_POLYMARKET_QUOTE_TIMEOUT_SEC,
                                                           ManualPolymarketDesk, validate_polymarket_place_request)
from unabated_ticket.bets_service.novig_place import (NOVIG_PLACE_READ_SOURCES, RUNNER_NOVIG_QUOTE_PATH,
                                                      RUNNER_NOVIG_QUOTE_TIMEOUT_SEC, ManualNovigDesk,
                                                      validate_novig_place_request)
from unabated_ticket.bets_service.normalize import parse_iso_ms, utc_now_iso
from unabated_ticket.bets_service.novig_book_stream import NovigBookStream, stream_url
from unabated_ticket.bets_service.novig_catalog_reads import NovigCatalogReads, pooled_transport, validate_read_request
from unabated_ticket.bets_service.novig_orders import NovigError, NovigOrderClient
from unabated_ticket.bets_service.novig_taker_desk import NovigTakerDesk, account_problem, novig_unseen_fill
from unabated_ticket.bets_service.novig_taker_desk import validate_order_request as validate_novig_order_request
from unabated_ticket.bets_service.sources import Source, bet105, draftkings
from unabated_ticket.bets_service.sources.betonline import source_if_configured as betonline_source_if_configured
from unabated_ticket.bets_service.sources.bfa import BFASource, normalize_open_bets
from unabated_ticket.bets_service.sources.bfa import source_if_configured as bfa_source_if_configured
from unabated_ticket.bets_service.sources.kalshi import PRIMARY_SUBACCOUNT, KalshiSource
from unabated_ticket.bets_service.sources.novig import source_if_connected as novig_source_if_connected
from unabated_ticket.bets_service.sources.novig_subaccount import NovigSubaccountSource
from unabated_ticket.bets_service.sources.polymarket_us import source_if_configured as polymarket_us_source_if_configured
from unabated_ticket.bets_service.sources.wagerzon import source_if_configured as wagerzon_source_if_configured
from unabated_ticket.bets_service.store import EDGE_SETTINGS_FIELDS, BetsStore
from unabated_ticket.bets_service.taker_orders import (KalshiOrderClient, TakerDesk, key_lock_problem,
                                                      validate_order_request)
from unabated_ticket.bets_service.taker_store import DECISION_REASONS, TakerStore

log = logging.getLogger(__name__)

POLL_TICK_SEC = 1.0
# After a failed poll a source waits at least this long before the next try: the
# fast cadences (Kalshi, Novig and Polymarket US every 10 s) must not hammer a venue
# that is down or refusing us. One success puts it back on its own cadence.
FAILED_POLL_RETRY_SEC = 60.0
# Sources polled on a poll thread of their own, one lane per tuple: Kalshi's first
# poll after a restart is the slowest read in the service (throttled GETs per live
# market); Novig's and Polymarket US's 10 s cadences only hold when nothing else
# runs between their polls (2026-10-29: in the shared loop BFA, BetOnline and
# Wagerzon ran between them and every venue waited ~35 s). Novig's three sources
# share a lane: each read is ~1 s or less, and one queue keeps their requests apart.
OWN_LANES = (("kalshi",), ("novig", "novig_bot", "novig_place"), ("polymarket_us",))
# A Place that may hold contracts reads its venue's sources before it answers
# (Cal, 2026-10-29), so the panel's and phone's /bets.json read right after the
# reply already counts it and the next Place on the game is not refused as
# unseen (manual_orders.unseen_manual_fill compares the Place time with the
# venue's read). Novig's three accounts are read out of turn, never behind
# BFA, BetOnline or Wagerzon in the shared lane.
PLACE_READ_SOURCES = {
    "kalshi": ("kalshi",),
    "novig": ("novig", "novig_place", "novig_bot"),
    "polymarket_us": ("polymarket_us",),
}
# Statuses after which an order may hold contracts (store._SELECT_MANUAL_MAYBE_FILLED;
# Novig's "unfilled" too, which the store still counts until Novig's own read).
PLACE_READ_STATUSES = ("filled", "partial", "unknown")
NOVIG_PLACE_READ_STATUSES = PLACE_READ_STATUSES + ("unfilled",)
# Reads already running are waited for at most this long in all; past it the
# reply goes and the next scheduled poll counts the fill, as before.
# The phone waits 90 s for a Place, and a Novig click can take ~72 s itself.
PLACE_READ_LOCK_WAIT_SEC = 10.0
MAX_DAYS = 3650
# A crosswalk POST is a few rows per poll; anything near this is not the panel.
MAX_BODY_BYTES = 1024 * 1024
MAX_CROSSWALK_ROWS_PER_POST = 1000
# The only Host headers a request can carry (DNS rebinding sends its own name).
LOOPBACK_HOST_NAMES = ("127.0.0.1", "localhost")
HTTP_DEFAULT_PORT = 80
# config.EXTRA_ALLOWED_HOSTS names arrive through `tailscale serve`'s HTTPS
# proxy, so the browser sends them bare (443 is the scheme's default) or, from
# a client that spells it out, with :443 — never with this service's port.
HTTPS_DEFAULT_PORT = 443
CROSSWALK_REQUIRED_FIELDS = ("venue", "league", "venueTeamKey", "unabatedTeamId")
CROSSWALK_OPTIONAL_FIELDS = ("venueTeamName", "unabatedTeamName", "learnedFrom")
# A bet names at most two teams, so one attach teaches at most two rows.
MAX_CROSSWALK_ROWS_PER_PIN = 2
PIN_REQUIRED_FIELDS = ("betId", "venue", "league", "eventId")
PIN_OPTIONAL_ID_FIELDS = ("awayTeamId", "homeTeamId")
PIN_OPTIONAL_TEXT_FIELDS = ("eventStart", "awayTeamName", "homeTeamName")
# Every `rows` POST that _validate_rows reads (fill snapshots, close snapshots,
# tau samples) is refused past this; the runner posts tau samples 500 at a
# time (runner.js TAU_SAMPLES_PER_POST), the snapshot POSTs carry a handful.
MAX_ROWS_PER_POST = 1000
# The runner posts its taker decisions 200 at a time (runner.js TAKER_LOG_ROWS_PER_POST).
MAX_TAKER_DECISION_ROWS_PER_POST = 1000
# A tau sample's seenAt is the runner's own clock; one further ahead of this
# server's clock than this would never be pruned or replaced (seen_at only
# moves forward), so it is refused.
MAX_TAU_SAMPLE_FUTURE_MS = 60 * 1000
# Cal's ask (2026-10-05): only these books carry bets that are not his.
EXCLUDABLE_VENUES = ("bfa", "wagerzon")
# A parlay's legs go in one request; a dozen is already a long ticket.
MAX_EXCLUSION_IDS_PER_POST = 50
# The panel's first sync sends every bet it had dismissed locally: a handful.
MAX_DISMISSAL_IDS_PER_POST = 200
# teaser.js marketKeyOf: "<Unabated eventId>:bt2" (spread) or ":bt3" (total).
TEASER_MARKET_KEY_PATTERN = re.compile(r"^[A-Za-z0-9_-]{1,64}:bt[23]$")
# pricecheck.js skipKeyOf: one line (game, period, market, side, number) at
# every book. A panel not yet reloaded still sends the old '@<bookId>' suffix;
# the store keeps only the line (store.edge_skip_line).
EDGE_SKIP_KEY_PATTERN = re.compile(r"^[A-Za-z0-9_-]{1,64}:\d{1,4}:\d{1,4}:\d{1,2}\|(ml|-?\d{1,4}(\.\d{1,3})?)(@\d{1,6})?$")
# American odds run from -100 down and +100 up; the gap between is no price.
MIN_AMERICAN_MAGNITUDE = 100
# The Edges settings' allowed values, as the panel's inputs enforce them
# (panel.js readSettingInputs / readEdgeSettingInputs). feed.PERIODS ids run
# 1 (full game) to 7 (4Q); feed.BET_TYPES are 1 moneyline, 2 spread, 3 total.
SETTINGS_PERIOD_IDS = range(1, 8)
SETTINGS_BET_TYPE_IDS = (1, 2, 3)
SETTINGS_BOOK_MODES = ("default", "all", "custom")
SETTINGS_SORT_KEYS = ("edge", "stake", "start", "exposure")
# Why POST /place_polymarket.json answers 503 when no Polymarket US key is configured.
POLYMARKET_PLACE_OFF_REASON = "no Polymarket US key is configured (POLYMARKET_US_KEY_ID / POLYMARKET_US_SECRET_KEY)"
# POST /novig_books.json: the most markets one ask may name (a taker pass names its candidates).
MAX_NOVIG_BOOK_MARKETS = 500
# Why POST /novig_books.json answers 503 when no bot Novig key is configured.
NOVIG_BOOK_STREAM_OFF_REASON = "no bot Novig key is configured (NOVIG_TAKER_KEY_ID / NOVIG_TAKER_PRIVATE_KEY_PATH)"
# Why POST /novig_read.json answers 503: its reads are signed with the same key.
NOVIG_CATALOG_READS_OFF_REASON = NOVIG_BOOK_STREAM_OFF_REASON
# Why POST /place_novig.json answers 503 when no Place key is configured.
NOVIG_PLACE_OFF_REASON = "no Place Novig key is configured (NOVIG_PLACE_KEY_ID / NOVIG_PLACE_PRIVATE_KEY_PATH)"

# The phone page (phone page plan step 2): every file the service serves, by
# exact URL path. Nothing is resolved from the request, so `..`, encoded
# slashes or a symlink cannot reach any other file. The extension modules are
# the panel's own pure ones, loaded by the page as plain <script>s in the
# order their globals need (index.html).
UNABATED_TICKET_DIR = Path(__file__).resolve().parent.parent
PHONE_DIR = UNABATED_TICKET_DIR / "server" / "phone"
TRACKER_DIR = UNABATED_TICKET_DIR / "server" / "tracker"
EXTENSION_DIR = UNABATED_TICKET_DIR / "extension"
PHONE_EXTENSION_MODULES = (
    "kelly.js", "place.js", "feed.js", "teams.js", "bets.js", "ladder.js", "condkelly.js", "betsview.js",
    "edgemove.js", "fillfair.js", "tailflex.js", "edgerows.js", "teaser.js", "teaserview.js",
    "balancesview.js",
)
STATIC_FILES: dict[str, Path] = {
    "/": PHONE_DIR / "index.html",
    "/index.html": PHONE_DIR / "index.html",
    "/phone.css": PHONE_DIR / "phone.css",
    "/phone.js": PHONE_DIR / "phone.js",
    "/phoneview.js": PHONE_DIR / "phoneview.js",
    **{f"/ext/{name}": EXTENSION_DIR / name for name in PHONE_EXTENSION_MODULES},
    # The Bet Tracker (server/tracker/): daily P&L and analysis off /bets.json.
    "/tracker": TRACKER_DIR / "index.html",
    "/tracker/": TRACKER_DIR / "index.html",
    "/tracker/tracker.css": TRACKER_DIR / "tracker.css",
    "/tracker/tracker.js": TRACKER_DIR / "tracker.js",
    "/tracker/trackerstats.js": TRACKER_DIR / "trackerstats.js",
}
STATIC_CONTENT_TYPES = {
    ".html": "text/html; charset=utf-8",
    ".css": "text/css; charset=utf-8",
    ".js": "text/javascript; charset=utf-8",
}
# The page loads only its own files and talks only to this origin.
PAGE_SECURITY_HEADERS = {
    "Content-Security-Policy": "default-src 'self'; frame-ancestors 'none'; base-uri 'none'; form-action 'none'",
    "X-Content-Type-Options": "nosniff",
    "Referrer-Policy": "no-referrer",
}
# The runner is local (or on the tailnet), never behind the sandbox's HTTP proxy.
_RUNNER_OPENER = urllib.request.build_opener(urllib.request.ProxyHandler({}))

# Sources with no poll here: the extension POSTs their records (/bet105.json,
# /draftkings.json).
# Listed so the panel reads "no completed poll yet" before the first push, not
# "no source configured".
PUSHED_SOURCES = (bet105.VENUE, draftkings.VENUE)


def _now() -> datetime:
    return datetime.now(timezone.utc)


def _iso(value: datetime) -> str:
    return value.strftime("%Y-%m-%dT%H:%M:%SZ")


def allowed_hosts(port: int, extra_hosts: tuple[str, ...] = ()) -> tuple[str, ...]:
    """The Host values a request to this port can carry. A browser omits the
    scheme's default port, so on 80 the bare names are valid too.
    `extra_hosts` (config.EXTRA_ALLOWED_HOSTS, validated there) are the names
    a reverse proxy on this machine forwards — the tailnet name behind
    `tailscale serve` — each bare and with :443."""
    with_port = tuple(f"{name}:{port}" for name in LOOPBACK_HOST_NAMES)
    loopback = with_port + LOOPBACK_HOST_NAMES if port == HTTP_DEFAULT_PORT else with_port
    proxied = tuple(value for name in extra_hosts for value in (name, f"{name}:{HTTPS_DEFAULT_PORT}"))
    return loopback + proxied


def host_allowed(host_header: str | None, port: int, extra_hosts: tuple[str, ...] = ()) -> bool:
    """Whether a request's Host names this service. A browser sends the name
    the page was loaded from, so a rebound evil.example does not match — the
    one check that survives the attacker being same-origin."""
    return (host_header or "").strip().lower() in allowed_hosts(port, extra_hosts)


def run_source_once(source: Source, store: BetsStore) -> bool:
    """One poll of one source: upsert its records and log the run. Never raises."""
    started_at = _now()
    try:
        records = source.fetch()
        written = store.upsert_bets(records, started_at)
        if written:  # unchanged records are not rewritten, so this is what moved
            log.info("source %s: %d of %d record(s) changed", source.name, written, len(records))
        store.log_source_run(source.name, started_at, _now(), True, None, len(records))
        return True
    except Exception as error:  # a source failure is a failed run, not a dead service
        log.warning("source %s failed: %s: %s", source.name, type(error).__name__, error)
        store.log_source_run(source.name, started_at, _now(), False,
                             f"{type(error).__name__}: {error}", 0)
        return False


def next_poll_delay(poll_sec: float, succeeded: bool) -> float:
    """Seconds until a source's next poll: its own cadence, or at least
    FAILED_POLL_RETRY_SEC after a failed one."""
    if succeeded:
        return poll_sec
    return max(poll_sec, FAILED_POLL_RETRY_SEC)


class SourceReader:
    """Runs each source's poll under that source's own lock, so no source is
    ever read by two threads at once: its poll lane on the cadence, and a Place
    that reads its venue out of turn (read_now)."""

    def __init__(self, sources: list[Source], store: BetsStore):
        self._sources = {source.name: source for source in sources}
        self._locks = {source.name: threading.Lock() for source in sources}
        self._store = store

    def run(self, source: Source) -> bool:
        """One poll of `source` (run_source_once), waiting for a read already running."""
        with self._locks[source.name]:
            return run_source_once(source, self._store)

    def read_now(self, names: tuple[str, ...], lock_wait_sec: float) -> list[str]:
        """Poll each registered source in `names` now, in order, waiting at
        most `lock_wait_sec` IN ALL for reads of them already running (a first
        Kalshi poll after a restart takes minutes; one budget so Novig's three
        cannot stack three waits onto a click). Returns the names NOT read, so
        the caller can say why; a name no source has is skipped silently (a
        venue whose sub-account is not read)."""
        missed = []
        deadline = time.monotonic() + lock_wait_sec
        for name in names:
            source = self._sources.get(name)
            if source is None:
                continue
            if not self._locks[name].acquire(timeout=max(0.0, deadline - time.monotonic())):
                missed.append(name)
                continue
            try:
                if not run_source_once(source, self._store):
                    missed.append(name)
            finally:
                self._locks[name].release()
        return missed


def poll_loop(sources: list[Source], store: BetsStore, stop: threading.Event,
              balance_book: BalanceBook | None = None, reader: SourceReader | None = None) -> None:
    """Run each source on its own cadence until `stop` (backed off to
    FAILED_POLL_RETRY_SEC after a failure). A store write that
    raises (disk full, a locked file) is logged and retried next poll — the
    thread must outlive it, or the HTTP side would serve ageing data with no
    failed run to show for it."""
    reader = reader or SourceReader(sources, store)
    next_due = {source.name: 0.0 for source in sources}
    # Balances first: a first Kalshi poll takes minutes and the loop is serial, so
    # reading them only after each fetch left the Balances tab empty after a restart.
    if balance_book is not None:
        for source in sources:
            read_balance_once(source, balance_book)
    while not stop.is_set():
        now = time.monotonic()
        for source in sources:
            if now >= next_due[source.name]:
                try:
                    succeeded = reader.run(source)
                except Exception:  # noqa: BLE001 — the poll thread is the service's heartbeat
                    log.exception("store write failed for source %s; retrying next poll", source.name)
                    succeeded = False
                if balance_book is not None:
                    read_balance_once(source, balance_book)
                next_due[source.name] = time.monotonic() + next_poll_delay(source.poll_sec, succeeded)
        stop.wait(POLL_TICK_SEC)


def poll_lanes(sources: list[Source]) -> list[list[Source]]:
    """The sources grouped one list per poll thread: each OWN_LANES lane with
    any of its sources registered, then every other source together in one
    shared loop."""
    lanes = [[source for source in sources if source.name in lane_names] for lane_names in OWN_LANES]
    own_names = {name for lane_names in OWN_LANES for name in lane_names}
    shared = [source for source in sources if source.name not in own_names]
    return [lane for lane in lanes if lane] + ([shared] if shared else [])


# A registered source that has not finished a poll yet (the first Kalshi poll
# takes ~1-2 min: one throttled GET per market and per event). Without this
# entry the panel would read "no source configured" — the text for venues
# with no source at all (#115-#117).
NO_POLL_YET = {"fetchedAt": None, "ok": False, "error": "no completed poll yet", "count": 0}


def source_status(store: BetsStore, source_names: list[str]) -> dict[str, dict]:
    status = store.source_status()
    for name in source_names:
        status.setdefault(name, dict(NO_POLL_YET))
    return status


def bets_payload(store: BetsStore, days: int, source_names: list[str] = (), snapshots: bool = False) -> dict:
    now = _now()
    payload = {"generatedAt": _iso(now), "sources": source_status(store, list(source_names)),
            "bets": store.load_bets(days, now), "crosswalk": store.load_crosswalk(), "pins": store.load_pins(),
            "fillFairs": store.load_fill_fairs(days, now), "closingFairs": store.load_closing_fairs(days, now),
            "exclusions": store.load_exclusions(), "dismissals": store.load_dismissals(),
            "teaserBlocks": store.load_teaser_blocks(now), "edgeSkips": store.load_edge_skips(now)}
    if snapshots:
        payload["fillSnapshots"] = store.load_fill_snapshots(days, now)
        payload["closeSnapshots"] = store.load_close_snapshots(days, now)
    return payload


def _non_empty_string(value: object) -> bool:
    return isinstance(value, str) and value != ""


def _id_as_string(value: object) -> object:
    """Unabated's ids are numeric; the store keeps them as strings. Anything
    else is returned as it came, for the caller's type check to name."""
    if isinstance(value, int) and not isinstance(value, bool):
        return str(value)
    return value


def _is_iso_timestamp(value: str) -> bool:
    try:
        datetime.fromisoformat(value.replace("Z", "+00:00"))
    except ValueError:
        return False
    return True


def validate_pin_request(body: object) -> tuple[dict, list[dict]] | str:
    """(pin, crosswalk rows) of a POST /pins.json body, or an error message.
    Every crosswalk row must be the pin's own venue and league: an attach
    teaches only what that bet's venue calls that game's teams."""
    if not isinstance(body, dict) or not isinstance(body.get("pin"), dict):
        return "body must be an object with a `pin` object"
    raw = body["pin"]
    pin: dict = {}
    for field in PIN_REQUIRED_FIELDS:
        value = _id_as_string(raw.get(field)) if field == "eventId" else raw.get(field)
        if not _non_empty_string(value):
            return f"pin.{field} must be a non-empty string"
        pin[field] = value
    for field in PIN_OPTIONAL_ID_FIELDS:
        value = _id_as_string(raw.get(field))
        if value is not None and not _non_empty_string(value):
            return f"pin.{field} must be a non-empty string, a number or null"
        pin[field] = value
    for field in PIN_OPTIONAL_TEXT_FIELDS:
        value = raw.get(field)
        if value is not None and not isinstance(value, str):
            return f"pin.{field} must be a string or null"
        pin[field] = value
    if pin["eventStart"] is not None and not _is_iso_timestamp(pin["eventStart"]):
        return f"pin.eventStart must be an ISO timestamp or null, got {pin['eventStart']!r}"
    raw_rows = body.get("crosswalk", [])
    if not isinstance(raw_rows, list):
        return "crosswalk must be an array"
    if len(raw_rows) > MAX_CROSSWALK_ROWS_PER_PIN:
        return f"at most {MAX_CROSSWALK_ROWS_PER_PIN} crosswalk rows per pin, got {len(raw_rows)}"
    rows = validate_crosswalk_rows({"rows": raw_rows})
    if isinstance(rows, str):
        return f"crosswalk: {rows}"
    for index, row in enumerate(rows):
        if (row["venue"], row["league"]) != (pin["venue"], pin["league"]):
            return (f"crosswalk[{index}] is {row['venue']}/{row['league']}, "
                    f"the pin is {pin['venue']}/{pin['league']}")
    return pin, rows


def validate_crosswalk_rows(body: object) -> list[dict] | str:
    """The rows of a POST /crosswalk.json body in the store's shape, or an
    error message naming the first bad row. `unabatedTeamId` may arrive as an
    int (Unabated's ids are numeric) and is stored as its string."""
    if not isinstance(body, dict) or not isinstance(body.get("rows"), list):
        return "body must be an object with a `rows` array"
    raw_rows = body["rows"]
    if len(raw_rows) > MAX_CROSSWALK_ROWS_PER_POST:
        return f"at most {MAX_CROSSWALK_ROWS_PER_POST} rows per request, got {len(raw_rows)}"
    rows: list[dict] = []
    for index, raw in enumerate(raw_rows):
        if not isinstance(raw, dict):
            return f"rows[{index}] must be an object"
        row = {"unabatedTeamId": _id_as_string(raw.get("unabatedTeamId"))}
        for field in CROSSWALK_REQUIRED_FIELDS:
            value = row.get(field, raw.get(field))
            if not _non_empty_string(value):
                return f"rows[{index}].{field} must be a non-empty string"
            row[field] = value
        for field in CROSSWALK_OPTIONAL_FIELDS:
            value = raw.get(field)
            if value is not None and not isinstance(value, str):
                return f"rows[{index}].{field} must be a string or null"
            row[field] = value
        rows.append(row)
    return rows


def _is_number(value: object) -> bool:
    return isinstance(value, (int, float)) and not isinstance(value, bool) and math.isfinite(value)


def _is_whole_american(value: object) -> bool:
    return isinstance(value, int) and not isinstance(value, bool) and abs(value) >= MIN_AMERICAN_MAGNITUDE


def _iso_ms(value: object) -> int | None:
    """Epoch ms of an ISO time string, None for anything else (parse_iso_ms takes strings only)."""
    return parse_iso_ms(value) if isinstance(value, str) else None


def _validate_fill_fair_row(index: int, raw: object) -> dict | str:
    """One row of a POST /fill_fairs.json body in the store's shape, or what
    was expected and what was found."""
    if not isinstance(raw, dict):
        return f"rows[{index}] must be an object, got {type(raw).__name__}"
    for field in ("betId", "lineKey"):
        if not _non_empty_string(raw.get(field)):
            return f"rows[{index}].{field} must be a non-empty string, got {raw.get(field)!r}"
    points = raw.get("points")
    if points is not None and not _is_number(points):
        return f"rows[{index}].points must be a number or null, got {points!r}"
    fair = raw.get("fairAmerican")
    if not _is_whole_american(fair):
        return f"rows[{index}].fairAmerican must be a whole American price (an integer <= -100 or >= 100), got {fair!r}"
    observed_ms = _iso_ms(raw.get("fairObservedAt"))
    placed_ms = _iso_ms(raw.get("placedAt"))
    if observed_ms is None:
        return f"rows[{index}].fairObservedAt must be an ISO time, got {raw.get('fairObservedAt')!r}"
    if placed_ms is None:
        return f"rows[{index}].placedAt must be an ISO time, got {raw.get('placedAt')!r}"
    if observed_ms > placed_ms:
        return (f"rows[{index}].fairObservedAt must be at or before placedAt (the fair is read at the fill), "
                f"got {raw['fairObservedAt']} after {raw['placedAt']}")
    return {"betId": raw["betId"], "lineKey": raw["lineKey"], "points": points, "fairAmerican": fair,
            "fairObservedAt": raw["fairObservedAt"], "placedAt": raw["placedAt"]}


def _iso_datetime(value: object) -> datetime | None:
    """An aware datetime from an ISO time string, None for anything else."""
    ms = _iso_ms(value)
    return None if ms is None else datetime.fromtimestamp(ms / 1000, tz=timezone.utc)


def _validate_taker_decision_row(index: int, raw: object) -> dict | str:
    """One row of a POST /taker_decisions.json body as {decision (as sent),
    decided_at, event_start (parsed)}, or what was expected and what was found. Only the columns taker_store.py
    queries on are checked; the rest of the decision is kept as sent."""
    if not isinstance(raw, dict):
        return f"rows[{index}] must be an object, got {type(raw).__name__}"
    if not _non_empty_string(raw.get("rowKey")):
        return f"rows[{index}].rowKey must be a non-empty string, got {raw.get('rowKey')!r}"
    if raw.get("reason") not in DECISION_REASONS:
        return f"rows[{index}].reason must be one of {list(DECISION_REASONS)}, got {raw.get('reason')!r}"
    decided_at = _iso_datetime(raw.get("decidedAt"))
    if decided_at is None:
        return f"rows[{index}].decidedAt must be an ISO time, got {raw.get('decidedAt')!r}"
    event_start = _iso_datetime(raw.get("eventStart"))
    if raw.get("eventStart") is not None and event_start is None:
        return f"rows[{index}].eventStart must be an ISO time or null, got {raw.get('eventStart')!r}"
    for field in ("eventId", "unabatedPrice", "contracts"):
        value = raw.get(field)
        if value is not None and (not isinstance(value, int) or isinstance(value, bool)):
            return f"rows[{index}].{field} must be an integer or null, got {value!r}"
    for field in ("liveAllIn", "liveEdge"):
        if raw.get(field) is not None and not _is_number(raw.get(field)):
            return f"rows[{index}].{field} must be a number or null, got {raw.get(field)!r}"
    for field in ("league", "ticker", "side"):
        if raw.get(field) is not None and not isinstance(raw.get(field), str):
            return f"rows[{index}].{field} must be a string or null, got {raw.get(field)!r}"
    return {"decision": raw, "decided_at": decided_at, "event_start": event_start}


def validate_taker_decisions(body: object) -> list[dict] | str:
    """The rows of a POST /taker_decisions.json body, or an error naming the
    first bad row (the whole request is refused; the runner keeps its queue)."""
    if not isinstance(body, dict) or not isinstance(body.get("rows"), list):
        return f"body must be an object with a `rows` array, got {type(body).__name__}"
    raw_rows = body["rows"]
    if len(raw_rows) > MAX_TAKER_DECISION_ROWS_PER_POST:
        return f"at most {MAX_TAKER_DECISION_ROWS_PER_POST} rows per request, got {len(raw_rows)}"
    rows: list[dict] = []
    for index, raw in enumerate(raw_rows):
        row = _validate_taker_decision_row(index, raw)
        if isinstance(row, str):
            return row
        rows.append(row)
    return rows


def validate_fill_fair_rows(body: object) -> list[dict] | str:
    """The rows of a POST /fill_fairs.json body, or an error naming the first
    bad row. Each row: betId and lineKey non-empty strings; points a number or
    null (a moneyline); fairAmerican a whole American price (Unabated's bacr
    always is); fairObservedAt and placedAt ISO times, the fair observed at or
    before the placement — a saved fair is permanent, so a row that cannot be
    the fill's is refused rather than stored."""
    return _validate_rows(body, _validate_fill_fair_row)


def _is_american_price(value: object) -> bool:
    return _is_number(value) and abs(value) >= MIN_AMERICAN_MAGNITUDE


def _is_count(value: object) -> bool:
    return isinstance(value, int) and not isinstance(value, bool) and value >= 0


def _validate_fill_book_row(index: int, raw: object) -> dict | str:
    """One row of a POST /fill_books.json body in the store's shape, or what
    was expected and what was found."""
    if not isinstance(raw, dict):
        return f"rows[{index}] must be an object, got {type(raw).__name__}"
    if not _non_empty_string(raw.get("betId")):
        return f"rows[{index}].betId must be a non-empty string, got {raw.get('betId')!r}"
    if _iso_ms(raw.get("placedAt")) is None:
        return f"rows[{index}].placedAt must be an ISO time, got {raw.get('placedAt')!r}"
    if not _is_american_price(raw.get("fillPrice")):
        return f"rows[{index}].fillPrice must be an American price (<= -100 or >= 100), got {raw.get('fillPrice')!r}"
    for field in ("bookCount", "booksBetter"):
        if not _is_count(raw.get(field)):
            return f"rows[{index}].{field} must be a whole number >= 0, got {raw.get(field)!r}"
    if raw["booksBetter"] > raw["bookCount"]:
        return f"rows[{index}].booksBetter must be at most bookCount, got {raw['booksBetter']} > {raw['bookCount']}"
    best_book, best_price = raw.get("bestBook"), raw.get("bestPrice")
    if raw["bookCount"] == 0:
        if best_book is not None or best_price is not None:
            return f"rows[{index}] with bookCount 0 must have null bestBook and bestPrice, got {best_book!r}, {best_price!r}"
    else:
        if not _non_empty_string(best_book):
            return f"rows[{index}].bestBook must be a non-empty string, got {best_book!r}"
        if not _is_american_price(best_price):
            return f"rows[{index}].bestPrice must be an American price (<= -100 or >= 100), got {best_price!r}"
    return {"betId": raw["betId"], "placedAt": raw["placedAt"], "fillPrice": raw["fillPrice"],
            "bookCount": raw["bookCount"], "booksBetter": raw["booksBetter"], "bestBook": best_book, "bestPrice": best_price}


def validate_fill_book_rows(body: object) -> list[dict] | str:
    """The rows of a POST /fill_books.json body, or an error naming the first
    bad row: a saved row is permanent, so a malformed one is refused."""
    return _validate_rows(body, _validate_fill_book_row)


EDGE_MOVE_KINDS = ("none", "fair_to_you", "book_away", "fair_against")
MAX_SNAPSHOT_BOOKS = 200


def _number_or_none(value: object) -> bool:
    return value is None or _is_number(value)


def _validate_snapshot_price(where: str, raw: object) -> str | None:
    """What is wrong with one side's price in a books record, or None. Null is
    a side the book did not list."""
    if raw is None:
        return None
    if not isinstance(raw, dict):
        return f"{where} must be an object or null, got {type(raw).__name__}"
    if not _is_number(raw.get("price")):
        return f"{where}.price must be a number, got {raw.get('price')!r}"
    for field in ("sourcePrice", "liquidity"):
        if not _number_or_none(raw.get(field)):
            return f"{where}.{field} must be a number or null, got {raw.get(field)!r}"
    if not _is_count(raw.get("sourceFormat")):
        return f"{where}.sourceFormat must be a whole number, got {raw.get('sourceFormat')!r}"
    return None


def _validate_snapshot_books(where: str, raw: object) -> list[dict] | str:
    """A books record ([{bookId, book, exchange, venue, side, other}]) as
    sent, or what was expected and what was found."""
    if not isinstance(raw, list):
        return f"{where} must be a list, got {type(raw).__name__}"
    if len(raw) > MAX_SNAPSHOT_BOOKS:
        return f"{where} holds at most {MAX_SNAPSHOT_BOOKS} books, got {len(raw)}"
    books: list[dict] = []
    for index, entry in enumerate(raw):
        at = f"{where}[{index}]"
        if not isinstance(entry, dict):
            return f"{at} must be an object, got {type(entry).__name__}"
        if not _is_count(entry.get("bookId")):
            return f"{at}.bookId must be a whole number, got {entry.get('bookId')!r}"
        if not _non_empty_string(entry.get("book")):
            return f"{at}.book must be a non-empty string, got {entry.get('book')!r}"
        for field in ("exchange", "venue"):
            if not isinstance(entry.get(field), bool):
                return f"{at}.{field} must be true or false, got {entry.get(field)!r}"
        for field in ("side", "other"):
            error = _validate_snapshot_price(f"{at}.{field}", entry.get(field))
            if error:
                return error
        books.append({key: entry.get(key) for key in ("bookId", "book", "exchange", "venue", "side", "other")})
    return books


def _validate_fill_snapshot_row(index: int, raw: object) -> dict | str:
    """One row of a POST /fill_snapshots.json body in the store's shape, or
    what was expected and what was found."""
    if not isinstance(raw, dict):
        return f"rows[{index}] must be an object, got {type(raw).__name__}"
    for field in ("betId", "lineKey", "edgeMoveLineKey"):
        if not _non_empty_string(raw.get(field)):
            return f"rows[{index}].{field} must be a non-empty string, got {raw.get(field)!r}"
    if _iso_ms(raw.get("placedAt")) is None:
        return f"rows[{index}].placedAt must be an ISO time, got {raw.get('placedAt')!r}"
    for field in ("points", "mainPoints", "pointsFromMain", "edgeMoveFairDelta", "edgeMovePriceDelta"):
        if not _number_or_none(raw.get(field)):
            return f"rows[{index}].{field} must be a number or null, got {raw.get(field)!r}"
    main_fair = raw.get("mainFairAmerican")
    if main_fair is not None and not _is_whole_american(main_fair):
        return f"rows[{index}].mainFairAmerican must be a whole American price or null, got {main_fair!r}"
    if raw.get("edgeMove") not in EDGE_MOVE_KINDS:
        return f"rows[{index}].edgeMove must be one of {list(EDGE_MOVE_KINDS)}, got {raw.get('edgeMove')!r}"
    books = _validate_snapshot_books(f"rows[{index}].books", raw.get("books"))
    if isinstance(books, str):
        return books
    fields = ("betId", "placedAt", "lineKey", "points", "mainPoints", "mainFairAmerican", "pointsFromMain",
              "edgeMoveLineKey", "edgeMove", "edgeMoveFairDelta", "edgeMovePriceDelta")
    return {**{field: raw.get(field) for field in fields}, "books": books}


def _validate_close_snapshot_row(index: int, raw: object) -> dict | str:
    """One row of a POST /close_snapshots.json body in the store's shape, or
    what was expected and what was found: read before its game's start."""
    if not isinstance(raw, dict):
        return f"rows[{index}] must be an object, got {type(raw).__name__}"
    if not _non_empty_string(raw.get("betId")):
        return f"rows[{index}].betId must be a non-empty string, got {raw.get('betId')!r}"
    observed_ms = _iso_ms(raw.get("observedAt"))
    start_ms = _iso_ms(raw.get("eventStart"))
    if observed_ms is None:
        return f"rows[{index}].observedAt must be an ISO time, got {raw.get('observedAt')!r}"
    if start_ms is None:
        return f"rows[{index}].eventStart must be an ISO time, got {raw.get('eventStart')!r}"
    if observed_ms >= start_ms:
        return (f"rows[{index}].observedAt must be before eventStart (a close is read pregame), "
                f"got {raw['observedAt']} at or after {raw['eventStart']}")
    books = _validate_snapshot_books(f"rows[{index}].books", raw.get("books"))
    if isinstance(books, str):
        return books
    return {"betId": raw["betId"], "observedAt": raw["observedAt"], "eventStart": raw["eventStart"], "books": books}


def _validate_rows(body: object, validate_row) -> list[dict] | str:
    """The rows of a POST body, or an error naming the first bad row (the
    whole request is refused)."""
    if not isinstance(body, dict) or not isinstance(body.get("rows"), list):
        return f"body must be an object with a `rows` array, got {type(body).__name__}"
    raw_rows = body["rows"]
    if len(raw_rows) > MAX_ROWS_PER_POST:
        return f"at most {MAX_ROWS_PER_POST} rows per request, got {len(raw_rows)}"
    rows: list[dict] = []
    for index, raw in enumerate(raw_rows):
        row = validate_row(index, raw)
        if isinstance(row, str):
            return row
        rows.append(row)
    return rows


def validate_fill_snapshot_rows(body: object) -> list[dict] | str:
    """The rows of a POST /fill_snapshots.json body, or an error naming the
    first bad row: a saved snapshot is permanent, so a malformed one is refused."""
    return _validate_rows(body, _validate_fill_snapshot_row)


def validate_close_snapshot_rows(body: object) -> list[dict] | str:
    """The rows of a POST /close_snapshots.json body, or an error naming the first bad row."""
    return _validate_rows(body, _validate_close_snapshot_row)


def _is_tau_sample_key(value: object) -> bool:
    """The shape the runner's pool reads a rung off: "<bookId>|<rungKey>|<sideIndex>"
    (tailflex.js rungOfSampleKey). A stored key of another shape would make
    every pool measurement throw until the row is pruned."""
    if not _non_empty_string(value):
        return False
    parts = value.split("|")
    return len(parts) == 3 and all(parts) and parts[2] in ("0", "1")


def _validate_tau_sample_row(index: int, raw: object) -> dict | str:
    """One row of a POST /tau_samples.json body in the store's shape, or what
    was expected and what was found. seenAt is refused further ahead of this
    server's clock than MAX_TAU_SAMPLE_FUTURE_MS."""
    if not isinstance(raw, dict):
        return f"rows[{index}] must be an object, got {type(raw).__name__}"
    if not _is_tau_sample_key(raw.get("key")):
        return f"rows[{index}].key must be \"<bookId>|<rungKey>|<sideIndex>\", got {raw.get('key')!r}"
    if not _non_empty_string(raw.get("marketKey")):
        return f"rows[{index}].marketKey must be a non-empty string, got {raw.get('marketKey')!r}"
    for field in ("edgeUnabated", "edgeExchanges"):
        if not _is_number(raw.get(field)):
            return f"rows[{index}].{field} must be a finite number, got {raw.get(field)!r}"
    seen_at_ms = _iso_ms(raw.get("seenAt"))
    if seen_at_ms is None:
        return f"rows[{index}].seenAt must be an ISO time, got {raw.get('seenAt')!r}"
    latest_allowed_ms = int(_now().timestamp() * 1000) + MAX_TAU_SAMPLE_FUTURE_MS
    if seen_at_ms > latest_allowed_ms:
        return (f"rows[{index}].seenAt must be within {MAX_TAU_SAMPLE_FUTURE_MS // 1000}s of this server's clock "
                f"({utc_now_iso()}), got {raw.get('seenAt')!r}")
    return {"key": raw["key"], "marketKey": raw["marketKey"], "edgeUnabated": raw["edgeUnabated"],
            "edgeExchanges": raw["edgeExchanges"], "seenAt": raw["seenAt"]}


def validate_tau_sample_rows(body: object) -> list[dict] | str:
    """The rows of a POST /tau_samples.json body, or an error naming the first
    bad row (the whole request is refused; the runner drops that chunk)."""
    return _validate_rows(body, _validate_tau_sample_row)


def _validate_closing_fair_row(index: int, raw: object) -> dict | str:
    """One row of a POST /closing_fairs.json body in the store's shape, or what
    was expected and what was found."""
    if not isinstance(raw, dict):
        return f"rows[{index}] must be an object, got {type(raw).__name__}"
    for field in ("betId", "lineKey"):
        if not _non_empty_string(raw.get(field)):
            return f"rows[{index}].{field} must be a non-empty string, got {raw.get(field)!r}"
    points = raw.get("points")
    if points is not None and not _is_number(points):
        return f"rows[{index}].points must be a number or null, got {points!r}"
    fair = raw.get("fairAmerican")
    if not _is_whole_american(fair):
        return f"rows[{index}].fairAmerican must be a whole American price (an integer <= -100 or >= 100), got {fair!r}"
    observed_ms = _iso_ms(raw.get("fairObservedAt"))
    start_ms = _iso_ms(raw.get("eventStart"))
    if observed_ms is None:
        return f"rows[{index}].fairObservedAt must be an ISO time, got {raw.get('fairObservedAt')!r}"
    if start_ms is None:
        return f"rows[{index}].eventStart must be an ISO time, got {raw.get('eventStart')!r}"
    if observed_ms >= start_ms:
        return (f"rows[{index}].fairObservedAt must be before eventStart (a close is read pregame), "
                f"got {raw['fairObservedAt']} at or after {raw['eventStart']}")
    return {"betId": raw["betId"], "lineKey": raw["lineKey"], "points": points, "fairAmerican": fair,
            "fairObservedAt": raw["fairObservedAt"], "eventStart": raw["eventStart"]}


def validate_closing_fair_rows(body: object) -> list[dict] | str:
    """The rows of a POST /closing_fairs.json body, or an error naming the
    first bad row: the fill-fair shape with eventStart in place of placedAt,
    each observed before its game's start."""
    return _validate_rows(body, _validate_closing_fair_row)


def validate_exclusion_request(body: object) -> tuple[list[str], bool] | str:
    """(bet ids, excluded) of a POST /exclusions.json body, or what was
    expected and what was found."""
    if not isinstance(body, dict):
        return f"body must be an object, got {type(body).__name__}"
    bet_ids = body.get("betIds")
    if not isinstance(bet_ids, list) or not bet_ids or not all(_non_empty_string(bet_id) for bet_id in bet_ids):
        return f"betIds must be a non-empty list of bet ids, got {bet_ids!r}"
    if len(bet_ids) > MAX_EXCLUSION_IDS_PER_POST:
        return f"at most {MAX_EXCLUSION_IDS_PER_POST} betIds per request, got {len(bet_ids)}"
    if not isinstance(body.get("excluded"), bool):
        return f"excluded must be true or false, got {body.get('excluded')!r}"
    return list(dict.fromkeys(bet_ids)), body["excluded"]


def validate_dismissal_request(body: object) -> tuple[list[str], bool] | str:
    """(bet ids, dismissed) of a POST /dismissals.json body, or what was
    expected and what was found."""
    if not isinstance(body, dict):
        return f"body must be an object, got {type(body).__name__}"
    bet_ids = body.get("betIds")
    if not isinstance(bet_ids, list) or not bet_ids or not all(_non_empty_string(bet_id) for bet_id in bet_ids):
        return f"betIds must be a non-empty list of bet ids, got {bet_ids!r}"
    if len(bet_ids) > MAX_DISMISSAL_IDS_PER_POST:
        return f"at most {MAX_DISMISSAL_IDS_PER_POST} betIds per request, got {len(bet_ids)}"
    if not isinstance(body.get("dismissed"), bool):
        return f"dismissed must be true or false, got {body.get('dismissed')!r}"
    return list(dict.fromkeys(bet_ids)), body["dismissed"]


def validate_teaser_block_request(body: object) -> tuple[str, datetime, bool] | str:
    """(market key, game start, blocked) of a POST /teaser_blocks.json body,
    or what was expected and what was found."""
    if not isinstance(body, dict):
        return f"body must be an object, got {type(body).__name__}"
    market_key = body.get("marketKey")
    if not isinstance(market_key, str) or not TEASER_MARKET_KEY_PATTERN.match(market_key):
        return f"marketKey must look like '<eventId>:bt2' or '<eventId>:bt3', got {market_key!r}"
    event_start_ms = body.get("eventStartMs")
    if not _is_whole_number(event_start_ms) or event_start_ms <= 0:
        return f"eventStartMs must be the game's start in whole epoch milliseconds, got {event_start_ms!r}"
    if not isinstance(body.get("blocked"), bool):
        return f"blocked must be true or false, got {body.get('blocked')!r}"
    return market_key, datetime.fromtimestamp(event_start_ms / 1000, timezone.utc), body["blocked"]


def validate_edge_skip_request(body: object) -> tuple[str, datetime, bool] | str:
    """(skip key, game start, skipped) of a POST /edge_skips.json body, or
    what was expected and what was found."""
    if not isinstance(body, dict):
        return f"body must be an object, got {type(body).__name__}"
    skip_key = body.get("skipKey")
    if not isinstance(skip_key, str) or not EDGE_SKIP_KEY_PATTERN.match(skip_key):
        return f"skipKey must look like '<eventId>:<period>:<betType>:<side>|<points|ml>', got {skip_key!r}"
    event_start_ms = body.get("eventStartMs")
    if not _is_whole_number(event_start_ms) or event_start_ms <= 0:
        return f"eventStartMs must be the game's start in whole epoch milliseconds, got {event_start_ms!r}"
    if not isinstance(body.get("skipped"), bool):
        return f"skipped must be true or false, got {body.get('skipped')!r}"
    return skip_key, datetime.fromtimestamp(event_start_ms / 1000, timezone.utc), body["skipped"]


def _is_whole_number(value: object) -> bool:
    return isinstance(value, int) and not isinstance(value, bool)


def _id_list_error(field: str, value: object, allowed: object = None, non_empty: bool = False) -> str | None:
    """Why `value` is not a list of whole ids (within `allowed` when given), or None."""
    if not isinstance(value, list) or not all(_is_whole_number(item) and item >= 0 for item in value):
        return f"settings.{field} must be a list of whole ids, got {value!r}"
    if non_empty and not value:
        return f"settings.{field} must name at least one id, got []"
    if allowed is not None and any(item not in allowed for item in value):
        return f"settings.{field} must hold only {list(allowed)}, got {value!r}"
    return None


def _settings_value_error(field: str, value: object) -> str | None:
    """Why `value` is not a valid non-null value of settings `field`, or None."""
    if field == "bankroll" or field == "maxLineAgeHours":
        return None if _is_number(value) and value > 0 else f"settings.{field} must be a number above 0, got {value!r}"
    if field == "multiplier":
        return None if _is_number(value) and 0 < value <= 1 else f"settings.multiplier must be above 0 and at most 1, got {value!r}"
    if field in ("minEdgePct", "minShrunkEdgePct", "minStake", "minLiquidityToWin"):
        return None if _is_number(value) and value >= 0 else f"settings.{field} must be a number, 0 or more, got {value!r}"
    if field in ("includeAlts", "groupByMarket"):
        return None if isinstance(value, bool) else f"settings.{field} must be true or false, got {value!r}"
    if field == "leagues":
        return _id_list_error(field, value)
    if field == "periods":
        return _id_list_error(field, value, SETTINGS_PERIOD_IDS, non_empty=True)
    if field == "betTypes":
        return _id_list_error(field, value, SETTINGS_BET_TYPE_IDS, non_empty=True)
    if field == "bookIds":
        return _id_list_error(field, value)
    if field == "bookMode":
        return None if value in SETTINGS_BOOK_MODES else f"settings.bookMode must be one of {list(SETTINGS_BOOK_MODES)}, got {value!r}"
    if field == "sortBy":
        return None if value in SETTINGS_SORT_KEYS else f"settings.sortBy must be one of {list(SETTINGS_SORT_KEYS)}, got {value!r}"
    raise AssertionError(f"no validator for settings field {field!r}")


def validate_settings_update(body: object, held: dict) -> dict | str:
    """The full settings row after applying a PUT /settings.json body to the
    `held` one, or an error naming the first problem. Fields left out keep
    their held value; null resets one to the default. bookIds is set exactly
    when bookMode is 'custom' (the panel's own ticks), checked on the result."""
    if not isinstance(body, dict) or not isinstance(body.get("settings"), dict):
        return f"body must be an object with a `settings` object, got {type(body).__name__}"
    update = body["settings"]
    unknown = sorted(set(update) - set(EDGE_SETTINGS_FIELDS))
    if unknown:
        return f"unknown settings field(s) {unknown}; expected some of {list(EDGE_SETTINGS_FIELDS)}"
    for field, value in update.items():
        error = None if value is None else _settings_value_error(field, value)
        if error:
            return error
    merged = {field: held.get(field) for field in EDGE_SETTINGS_FIELDS}
    merged.update(update)
    if (merged["bookMode"] == "custom") != (merged["bookIds"] is not None):
        return (f"settings.bookIds must be a list exactly when bookMode is 'custom', "
                f"got bookMode {merged['bookMode']!r} with bookIds {merged['bookIds']!r}")
    return merged


# The runner's JSON routes the bets service passes through (one origin for the pages).
RUNNER_ROUTES = ("/edges.json", "/scenarios.json", "/teasers.json", "/attach.json", "/tau.json")
# The runner's own answer about the request (no such open bet, a game gone
# from the board) passes through as it came; any other failure is a 502.
RUNNER_PASSTHROUGH_ERRORS = (400, 404)


def fetch_runner_json(runner_url: str, path: str, timeout_sec: float, query: str = "") -> tuple[int, bytes]:
    """(status, body) for GET <path>?<query> (path one of RUNNER_ROUTES): the
    runner's body as it came on a 200 or a 400/404, else 502 with a JSON error
    naming the runner URL and what failed. Network: one GET of
    <runner_url><path>, no proxy, `timeout_sec`."""
    url = f"{runner_url}{path}{'?' + query if query else ''}"
    try:
        with _RUNNER_OPENER.open(url, timeout=timeout_sec) as response:
            return 200, response.read()
    except urllib.error.HTTPError as error:
        if error.code in RUNNER_PASSTHROUGH_ERRORS:
            return error.code, error.read()
        problem = f"answered HTTP {error.code}"
    except (urllib.error.URLError, OSError) as error:  # refused, reset, DNS, timeout
        problem = f"unreachable ({getattr(error, 'reason', None) or error})"
    error_body = {"error": f"server runner at {runner_url} {problem} for {path}; start it with "
                           f"node unabated_ticket/server/runner.js", "runnerUrl": runner_url}
    return 502, json.dumps(error_body).encode()


def health_payload(store: BetsStore, started_at: float, source_names: list[str] = (),
                   novig_book_stream: NovigBookStream | None = None) -> dict:
    payload = {"ok": True, "generatedAt": _iso(_now()), "uptimeSec": int(time.monotonic() - started_at),
               "sources": source_status(store, list(source_names))}
    if novig_book_stream is not None:
        payload["novigStream"] = novig_book_stream.status()
    return payload


def validate_novig_books_request(body: object) -> tuple[list[str], float, bool] | str:
    """(market ids, wait seconds, priority) or why the body is refused."""
    if not isinstance(body, dict):
        return "expected a JSON object body"
    market_ids = body.get("marketIds")
    if not isinstance(market_ids, list) or len(market_ids) > MAX_NOVIG_BOOK_MARKETS:
        return f"marketIds must be a list of at most {MAX_NOVIG_BOOK_MARKETS} market ids"
    for index, market_id in enumerate(market_ids):
        if not _non_empty_string(market_id):
            return f"marketIds[{index}] must be a non-empty string, got {market_id!r}"
    wait_ms = body.get("waitMs", 0)
    if not _is_number(wait_ms) or wait_ms < 0:
        return f"waitMs must be a number >= 0, got {wait_ms!r}"
    priority = body.get("priority", False)
    if not isinstance(priority, bool):
        return f"priority must be true or false, got {priority!r}"
    return list(dict.fromkeys(market_ids)), wait_ms / 1000, priority


def parse_days(query: str) -> int | str:
    """`days` from a query string, or an error message."""
    raw = parse_qs(query).get("days")
    if not raw:
        return config.RETENTION_DAYS
    try:
        days = int(raw[0])
    except ValueError:
        return f"days must be an integer, got {raw[0]!r}"
    if not 0 <= days <= MAX_DAYS:
        return f"days must be between 0 and {MAX_DAYS}, got {days}"
    return days


def make_handler(store: BetsStore, started_at: float, source_names: list[str] = (),
                 runner_url: str | None = None, runner_timeout_sec: float | None = None,
                 extra_hosts: tuple[str, ...] | None = None,
                 teaser_placer: bfa_teaser.BFATeaserPlacer | None = None,
                 balance_book: BalanceBook | None = None,
                 taker_store: TakerStore | None = None,
                 taker_desk: TakerDesk | None = None,
                 taker_desk_off_reason: str = "no bot Kalshi key is configured (KALSHI_TAKER_API_KEY_ID / "
                                             "KALSHI_TAKER_PRIVATE_KEY_PATH)",
                 kalshi_place_desk: ManualKalshiDesk | None = None,
                 kalshi_place_off_reason: str = "no Kalshi trade key for Place is configured "
                                                "(KALSHI_PLACE_API_KEY_ID / KALSHI_PLACE_PRIVATE_KEY_PATH)",
                 novig_taker_desk: NovigTakerDesk | None = None,
                 novig_taker_desk_off_reason: str = "no bot Novig key is configured (NOVIG_TAKER_KEY_ID / "
                                                    "NOVIG_TAKER_PRIVATE_KEY_PATH)",
                 novig_place_desk: ManualNovigDesk | None = None,
                 novig_place_off_reason: str = NOVIG_PLACE_OFF_REASON,
                 polymarket_place_desk: ManualPolymarketDesk | None = None,
                 polymarket_place_off_reason: str = POLYMARKET_PLACE_OFF_REASON,
                 novig_book_stream: NovigBookStream | None = None,
                 novig_catalog_reads: NovigCatalogReads | None = None,
                 source_reader: SourceReader | None = None) -> type[BaseHTTPRequestHandler]:
    """The request handler class. `runner_url` / `runner_timeout_sec` /
    `extra_hosts` default to config.RUNNER_URL / RUNNER_TIMEOUT_SEC /
    EXTRA_ALLOWED_HOSTS (the tests pass their own). `teaser_placer` is the
    Place button's (None: no BFA account is read, POST /place_teaser.json 503s).
    `taker_store` is the taker's decision log (None: POST /taker_decisions.json 503s).
    `taker_desk` places the bot's orders (None: POST /taker/order.json 503s with `taker_desk_off_reason`).
    `kalshi_place_desk` places Cal's own Place-button orders (None: POST
    /place_kalshi.json 503s with `kalshi_place_off_reason`).
    `novig_taker_desk` places the bot's Novig orders (None: POST
    /taker/novig_order.json 503s with `novig_taker_desk_off_reason`). `novig_place_desk`
    places Cal's Novig Place-button orders (None: POST /place_novig.json 503s
    with `novig_place_off_reason`). `polymarket_place_desk` places Cal's
    Polymarket US Place-button orders (None: POST /place_polymarket.json 503s
    with `polymarket_place_off_reason`). `novig_book_stream` serves Novig's live
    books (None: POST /novig_books.json 503s), `novig_catalog_reads` the keyed
    catalog reads on the same key (None: POST /novig_read.json 503s). `source_reader` reads a
    Place's venue before its reply (None: the reply goes at once and the next
    scheduled poll counts the fill)."""
    names = list(source_names)
    proxied_host_names = tuple(extra_hosts) if extra_hosts is not None else config.EXTRA_ALLOWED_HOSTS
    edges_runner_url = (runner_url or config.RUNNER_URL).rstrip("/")
    edges_timeout_sec = runner_timeout_sec if runner_timeout_sec is not None else config.RUNNER_TIMEOUT_SEC

    class BetsHandler(BaseHTTPRequestHandler):
        def do_GET(self) -> None:  # noqa: N802 — http.server's name
            if self._refused_foreign_host():
                return
            url = urlparse(self.path)
            if url.path == "/health":
                self._send_json(200, health_payload(store, started_at, names, novig_book_stream))
                return
            if url.path == "/settings.json":
                self._send_json(200, store.load_edge_settings())
                return
            if url.path == "/balances.json":
                self._send_json(200, balances_payload(balance_book or BalanceBook()))
                return
            if url.path == "/tau_samples.json":
                now = _now()
                self._send_json(200, {"generatedAt": _iso(now),
                                      "rows": store.load_tau_samples(config.TAU_SAMPLES_RETENTION_HOURS, now)})
                return
            if url.path == "/bets.json":
                days = parse_days(url.query)
                if isinstance(days, str):
                    self._send_json(400, {"error": days})
                    return
                snapshots = parse_qs(url.query).get("snapshots") == ["1"]
                self._send_json(200, bets_payload(store, days, names, snapshots))
                return
            if url.path in RUNNER_ROUTES:
                status, body = fetch_runner_json(edges_runner_url, url.path, edges_timeout_sec, url.query)
                self._send_bytes(status, body, "application/json")
                return
            if url.path in STATIC_FILES:
                self._send_static(STATIC_FILES[url.path])
                return
            self._send_json(404, {"error": f"no route for {url.path}"})

        def _send_static(self, path: Path) -> None:
            try:
                body = path.read_bytes()
            except OSError as error:
                log.error("phone page file %s unreadable: %s", path, error)
                self._send_json(404, {"error": f"phone page file {path.name} is missing on this server"})
                return
            self._send_bytes(200, body, STATIC_CONTENT_TYPES[path.suffix], PAGE_SECURITY_HEADERS)

        def do_POST(self) -> None:  # noqa: N802 — http.server's name
            if self._refused_foreign_host():
                return
            url = urlparse(self.path)
            if url.path not in ("/crosswalk.json", "/pins.json", "/fill_fairs.json", "/fill_books.json", "/closing_fairs.json",
                                "/fill_snapshots.json", "/close_snapshots.json", "/tau_samples.json",
                                "/exclusions.json", "/dismissals.json", "/teaser_blocks.json", "/edge_skips.json", "/bet105.json",
                                "/draftkings.json", "/place_teaser.json", "/place_kalshi.json", "/place_novig.json",
                                "/place_polymarket.json", "/taker_decisions.json",
                                "/taker/order.json", "/taker/novig_order.json", "/novig_books.json", "/novig_read.json"):
                self._send_json(404, {"error": f"no route for POST {url.path}"})
                return
            body = self._read_json_body()
            if isinstance(body, tuple):
                self._send_json(*body)
                return
            if url.path == "/pins.json":
                self._post_pin(body)
                return
            if url.path == "/fill_fairs.json":
                self._save_fill_fairs(body)
                return
            if url.path == "/fill_books.json":
                self._save_fill_books(body)
                return
            if url.path == "/closing_fairs.json":
                self._save_closing_fairs(body)
                return
            if url.path == "/fill_snapshots.json":
                self._save_rows(validate_fill_snapshot_rows(body), store.save_fill_snapshots)
                return
            if url.path == "/close_snapshots.json":
                self._save_rows(validate_close_snapshot_rows(body), store.save_close_snapshots)
                return
            if url.path == "/tau_samples.json":
                self._save_rows(validate_tau_sample_rows(body),
                                lambda rows, now: store.save_tau_samples(rows, config.TAU_SAMPLES_RETENTION_HOURS, now))
                return
            if url.path == "/exclusions.json":
                self._set_exclusions(body)
                return
            if url.path == "/dismissals.json":
                self._set_dismissals(body)
                return
            if url.path == "/teaser_blocks.json":
                self._set_teaser_block(body)
                return
            if url.path == "/edge_skips.json":
                self._set_edge_skip(body)
                return
            if url.path == "/bet105.json":
                self._push_bet105(body)
                return
            if url.path == "/draftkings.json":
                self._push_draftkings(body)
                return
            if url.path == "/place_teaser.json":
                self._place_teaser(body)
                return
            if url.path == "/place_kalshi.json":
                self._place_kalshi(body)
                return
            if url.path == "/place_novig.json":
                self._place_novig(body)
                return
            if url.path == "/place_polymarket.json":
                self._place_polymarket(body)
                return
            if (url.path in ("/taker_decisions.json", "/taker/order.json", "/taker/novig_order.json", "/novig_books.json",
                             "/novig_read.json")
                    and self._refused_non_loopback()):
                return
            if url.path == "/novig_books.json":
                self._novig_books(body)
                return
            if url.path == "/novig_read.json":
                self._novig_read(body)
                return
            if url.path == "/taker_decisions.json":
                self._log_taker_decisions(body)
                return
            if url.path == "/taker/order.json":
                self._place_taker_order(body)
                return
            if url.path == "/taker/novig_order.json":
                self._place_novig_taker_order(body)
                return
            rows = validate_crosswalk_rows(body)
            if isinstance(rows, str):
                self._send_json(400, {"error": rows})
                return
            result = store.learn_crosswalk(rows, _now())
            self._send_json(200, {"ok": True, "learned": result["learned"], "conflicts": result["conflicts"],
                                  "crosswalk": store.load_crosswalk()})

        # PUT /settings.json: merge the body into the held row, UPSERT it, and
        # reply with what is now stored.
        def do_PUT(self) -> None:  # noqa: N802 — http.server's name
            if self._refused_foreign_host():
                return
            url = urlparse(self.path)
            if url.path != "/settings.json":
                self._send_json(404, {"error": f"no route for PUT {url.path}"})
                return
            body = self._read_json_body()
            if isinstance(body, tuple):
                self._send_json(*body)
                return
            merged = validate_settings_update(body, store.load_edge_settings()["settings"])
            if isinstance(merged, str):
                self._send_json(400, {"error": merged})
                return
            store.save_edge_settings(merged, _now())
            self._send_json(200, {"ok": True, **store.load_edge_settings()})

        def _post_pin(self, body: object) -> None:
            request = validate_pin_request(body)
            if isinstance(request, str):
                self._send_json(400, {"error": request})
                return
            pin, rows = request
            venue = store.bet_venue(pin["betId"])
            if venue is None:
                self._send_json(404, {"error": f"no bet with id {pin['betId']!r}"})
                return
            # The rows carry the pin's venue, so the pin must carry the bet's.
            if venue != pin["venue"]:
                self._send_json(400, {"error": f"bet {pin['betId']!r} is {venue}, the pin says {pin['venue']}"})
                return
            store.pin_bet(pin, rows, _now())
            self._send_json(200, {"ok": True, "pins": store.load_pins(), "crosswalk": store.load_crosswalk()})

        # POST /fill_fairs.json: INSERT the rows (first capture wins) and reply
        # with every saved fair of the bets /bets.json serves.
        def _save_fill_fairs(self, body: object) -> None:
            rows = validate_fill_fair_rows(body)
            if isinstance(rows, str):
                self._send_json(400, {"error": rows})
                return
            now = _now()
            saved = store.save_fill_fairs(rows, now)
            self._send_json(200, {"ok": True, "saved": saved,
                                  "fillFairs": store.load_fill_fairs(config.RETENTION_DAYS, now)})

        # POST /fill_books.json: INSERT the rows (first capture wins).
        def _save_fill_books(self, body: object) -> None:
            rows = validate_fill_book_rows(body)
            if isinstance(rows, str):
                self._send_json(400, {"error": rows})
                return
            self._send_json(200, {"ok": True, "saved": store.save_fill_books(rows, _now())})

        # POST /fill_snapshots.json (insert-only), /close_snapshots.json
        # (newest pre-start wins) and /tau_samples.json (newest seenAt wins,
        # rows past the retention window pruned in the same write): validated
        # rows, or a 400 naming the first bad one.
        def _save_rows(self, rows: list[dict] | str, save) -> None:
            if isinstance(rows, str):
                self._send_json(400, {"error": rows})
                return
            self._send_json(200, {"ok": True, "saved": save(rows, _now())})

        # POST /closing_fairs.json: UPSERT the rows, each kept only when newer.
        def _save_closing_fairs(self, body: object) -> None:
            rows = validate_closing_fair_rows(body)
            if isinstance(rows, str):
                self._send_json(400, {"error": rows})
                return
            self._send_json(200, {"ok": True, "saved": store.save_closing_fairs(rows, _now())})

        # POST /exclusions.json: remove (or restore) BFA / Wagerzon bets from
        # the tracker's numbers. Every id must be a stored bet at one of
        # EXCLUDABLE_VENUES, or nothing changes.
        def _set_exclusions(self, body: object) -> None:
            request = validate_exclusion_request(body)
            if isinstance(request, str):
                self._send_json(400, {"error": request})
                return
            bet_ids, excluded = request
            venue_by_bet = {bet_id: store.bet_venue(bet_id) for bet_id in bet_ids}
            unknown = [bet_id for bet_id, venue in venue_by_bet.items() if venue is None]
            if unknown:
                self._send_json(404, {"error": f"no bet with id(s) {unknown}"})
                return
            other = {bet_id: venue for bet_id, venue in venue_by_bet.items() if venue not in EXCLUDABLE_VENUES}
            if other:
                self._send_json(400, {"error": f"only {list(EXCLUDABLE_VENUES)} bets can be removed, got {other}"})
                return
            changed = store.set_exclusions(bet_ids, venue_by_bet, excluded, _now())
            self._send_json(200, {"ok": True, "changed": changed, "exclusions": store.load_exclusions()})

        # POST /dismissals.json: stop (or resume) flagging open bets no board
        # game matches. Every id must be a stored bet, or nothing changes.
        def _set_dismissals(self, body: object) -> None:
            request = validate_dismissal_request(body)
            if isinstance(request, str):
                self._send_json(400, {"error": request})
                return
            bet_ids, dismissed = request
            unknown = [bet_id for bet_id in bet_ids if store.bet_venue(bet_id) is None]
            if unknown:
                self._send_json(404, {"error": f"no bet with id(s) {unknown}"})
                return
            changed = store.set_dismissals(bet_ids, dismissed, _now())
            self._send_json(200, {"ok": True, "changed": changed, "dismissals": store.load_dismissals()})

        # POST /teaser_blocks.json: mark (or clear) one market Can't tease.
        def _set_teaser_block(self, body: object) -> None:
            request = validate_teaser_block_request(body)
            if isinstance(request, str):
                self._send_json(400, {"error": request})
                return
            market_key, event_start, blocked = request
            now = _now()
            changed = store.set_teaser_block(market_key, event_start, blocked, now)
            self._send_json(200, {"ok": True, "changed": changed, "teaserBlocks": store.load_teaser_blocks(now)})

        # POST /edge_skips.json: skip (or restore) one Edges line.
        def _set_edge_skip(self, body: object) -> None:
            request = validate_edge_skip_request(body)
            if isinstance(request, str):
                self._send_json(400, {"error": request})
                return
            skip_key, event_start, skipped = request
            now = _now()
            changed = store.set_edge_skip(skip_key, event_start, skipped, now)
            self._send_json(200, {"ok": True, "changed": changed, "edgeSkips": store.load_edge_skips(now)})

        # POST /bet105.json: the extension's read of the account, as one source
        # run — a complete push UPSERTs its records (a graded wager over its open
        # record) and closes the open ones neither read lists; an error push is a
        # failed run and the records stand.
        def _log_taker_decisions(self, body: object) -> None:
            if taker_store is None:
                self._send_json(503, {"error": "the taker log is not open on this service"})
                return
            rows = validate_taker_decisions(body)
            if isinstance(rows, str):
                self._send_json(400, {"error": rows})
                return
            self._send_json(200, {"ok": True, "logged": taker_store.append_decisions(rows, _now())})

        # POST /taker/order.json: one bot order (taker_orders.py). Its outcome
        # is always a 200 with {status}; 503 only when the bot has no key.
        def _place_taker_order(self, body: object) -> None:
            if taker_desk is None:
                self._send_json(503, {"status": "off", "error": taker_desk_off_reason})
                return
            request = validate_order_request(body)
            if isinstance(request, str):
                self._send_json(400, {"error": request})
                return
            self._send_json(200, taker_desk.place(request))

        # POST /taker/novig_order.json: one bot order on Novig (novig_taker_desk.py),
        # a 200 with {status} always; 503 only when the bot has no Novig key.
        def _place_novig_taker_order(self, body: object) -> None:
            if novig_taker_desk is None:
                self._send_json(503, {"status": "off", "error": novig_taker_desk_off_reason})
                return
            request = validate_novig_order_request(body)
            if isinstance(request, str):
                self._send_json(400, {"error": request})
                return
            self._send_json(200, novig_taker_desk.place(request))

        # POST /novig_books.json: Novig's live books for the runner (novig_book_stream.py).
        def _novig_books(self, body: object) -> None:
            if novig_book_stream is None:
                self._send_json(503, {"error": NOVIG_BOOK_STREAM_OFF_REASON})
                return
            request = validate_novig_books_request(body)
            if isinstance(request, str):
                self._send_json(400, {"error": request})
                return
            market_ids, wait_sec, priority = request
            self._send_json(200, novig_book_stream.books_for(market_ids, wait_sec, priority))

        # POST /novig_read.json: one keyed Novig catalog GET for the runner (novig_catalog_reads.py).
        def _novig_read(self, body: object) -> None:
            if novig_catalog_reads is None:
                self._send_json(503, {"error": NOVIG_CATALOG_READS_OFF_REASON})
                return
            request = validate_read_request(body)
            if isinstance(request, str):
                self._send_json(400, {"error": request})
                return
            path, query = request
            try:
                reply = novig_catalog_reads.read(path, query)
            except OSError as error:
                self._send_json(502, {"error": f"Novig did not answer GET {path}: {error}"})
                return
            self._send_json(200, reply)

        def _push_bet105(self, body: object) -> None:
            push = bet105.validate_push(body)
            if isinstance(push, str):
                self._send_json(400, {"error": push})
                return
            started_at = _now()
            if "error" in push:
                store.log_source_run(bet105.VENUE, started_at, _now(), False, push["error"], 0)
                self._send_json(200, {"ok": True, "recorded": "error"})
                return
            stored = store.load_bets(config.RETENTION_DAYS, started_at)
            open_records = bet105.normalize_bet105(push["feeds"], push["fetchedAt"])
            known_feeds = bet105.feeds_by_native_id(stored + open_records)
            settled, skipped = bet105.normalize_settled(push["settled"], push["fetchedAt"], known_feeds)
            if skipped:
                log.warning("bet105: %d settled wager(s) not read: %s", len(skipped), "; ".join(skipped))
            records = bet105.merge_settled(open_records, settled)
            closed = bet105.closed_by_absence(stored, {record["id"] for record in records}, _iso(started_at))
            store.upsert_bets(records + closed, started_at)
            store.log_source_run(bet105.VENUE, started_at, _now(), True, None, len(records))
            self._send_json(200, {"ok": True, "count": len(records), "settled": len(settled),
                                  "closed": len(closed), "skipped": len(skipped)})

        # POST /draftkings.json: the extension's read of the account, as one source run —
        # a complete push UPSERTs its records (a settled bet over its open copy) and
        # closes the open ones neither list carries; an error push is a failed run and
        # the records stand.
        def _push_draftkings(self, body: object) -> None:
            push = draftkings.validate_push(body)
            if isinstance(push, str):
                self._send_json(400, {"error": push})
                return
            started_at = _now()
            if "error" in push:
                store.log_source_run(draftkings.VENUE, started_at, _now(), False, push["error"], 0)
                self._send_json(200, {"ok": True, "recorded": "error"})
                return
            records, skipped = draftkings.normalize_draftkings(push)
            if skipped:
                log.warning("draftkings: %d settled bet(s) not read: %s", len(skipped), "; ".join(skipped))
            stored = store.load_bets(config.RETENTION_DAYS, started_at)
            closed = draftkings.closed_by_absence(stored, {record["id"] for record in records}, _iso(started_at))
            store.upsert_bets(records + closed, started_at)
            store.log_source_run(draftkings.VENUE, started_at, _now(), True, None, len(records))
            self._send_json(200, {"ok": True, "count": len(records), "closed": len(closed), "skipped": len(skipped)})

        # POST /place_teaser.json: one teaser at BFA (bfa_teaser.py). An exception out of
        # place() is raised before the wager goes out, so it reads "Not placed". The open
        # bets read while confirming are stored at once, the BFA source's own way.
        # The phone page and a panel pointed at the server place through the
        # tailnet name `tailscale serve` forwards (Cal, 2026-10-09); do_POST's
        # Host allowlist has already refused every other name.
        def _place_teaser(self, body: object) -> None:
            if teaser_placer is None:
                self._send_json(503, {"error": "no BFA account is read: set BFA_USERNAME and BFA_PASSWORD "
                                               "(bet_logger/.env in the main checkout)"})
                return
            request = bfa_teaser.validate_place_request(body)
            if isinstance(request, str):
                self._send_json(400, {"error": request})
                return
            try:
                result = teaser_placer.place(request)
            except Exception as error:  # noqa: BLE001 — place() answers for itself once the POST is out
                log.exception("bfa teaser: placement failed before the wager went out")
                result = {"status": bfa_teaser.STATUS_REFUSED, "message": f"Not placed: {type(error).__name__}: {error}"}
            open_wagers = result.pop("openWagers", None)
            if open_wagers:
                self._store_bfa_open_bets(open_wagers)
            self._send_json(200, {"ok": result["status"] == bfa_teaser.STATUS_PLACED, **result})

        # POST /place_kalshi.json: one order on Cal's own Kalshi key
        # (manual_orders.py). Like Place on a teaser, the phone page and a
        # panel pointed at the server send it through the tailnet name.
        def _place_kalshi(self, body: object) -> None:
            if kalshi_place_desk is None:
                self._send_json(503, {"error": kalshi_place_off_reason})
                return
            request = validate_kalshi_place_request(body)
            if isinstance(request, str):
                self._send_json(400, {"error": request})
                return
            try:
                result = kalshi_place_desk.place(request)
            except Exception as error:  # noqa: BLE001 — logged; the row (if any) says how far it got
                log.exception("kalshi place: failed")
                result = {"status": "unknown", "message": f"The order failed on this server ({type(error).__name__}: "
                                                          f"{error}). Check Kalshi before placing again."}
            if result["status"] in PLACE_READ_STATUSES:
                self._read_place_venue("kalshi")
            self._send_json(200, {"ok": result["status"] in ("filled", "partial"), **result})

        # POST /place_novig.json: one order from Place's Novig sub-account
        # (novig_place.py), through the tailnet name like Place on Kalshi.
        def _place_novig(self, body: object) -> None:
            if novig_place_desk is None:
                self._send_json(503, {"error": novig_place_off_reason})
                return
            request = validate_novig_place_request(body)
            if isinstance(request, str):
                self._send_json(400, {"error": request})
                return
            try:
                result = novig_place_desk.place(request)
            except Exception as error:  # noqa: BLE001 — logged; the row (if any) says how far it got
                log.exception("novig place: failed")
                result = {"status": "unknown", "message": f"The order failed on this server ({type(error).__name__}: "
                                                          f"{error}). Check Novig before placing again."}
            if result["status"] in NOVIG_PLACE_READ_STATUSES:
                self._read_place_venue("novig")
            self._send_json(200, {"ok": result["status"] in ("filled", "partial"), **result})

        # POST /place_polymarket.json: one order on Cal's main Polymarket US
        # account on the POLYMARKET_US_* key (polymarket_place.py), through the
        # tailnet name like Place on Kalshi.
        def _place_polymarket(self, body: object) -> None:
            if polymarket_place_desk is None:
                self._send_json(503, {"error": polymarket_place_off_reason})
                return
            request = validate_polymarket_place_request(body)
            if isinstance(request, str):
                self._send_json(400, {"error": request})
                return
            try:
                result = polymarket_place_desk.place(request)
            except Exception as error:  # noqa: BLE001 — logged; the row (if any) says how far it got
                log.exception("polymarket place: failed")
                result = {"status": "unknown", "message": f"The order failed on this server ({type(error).__name__}: "
                                                          f"{error}). Check Polymarket US before placing again."}
            if result["status"] in PLACE_READ_STATUSES:
                self._read_place_venue("polymarket_us")
            self._send_json(200, {"ok": result["status"] in ("filled", "partial"), **result})

        def _read_place_venue(self, venue: str) -> None:
            """Read the venue a Place just traded at (PLACE_READ_SOURCES) before
            the reply goes. Never raises: the order already stands."""
            if source_reader is None:
                return
            try:
                missed = source_reader.read_now(PLACE_READ_SOURCES[venue], PLACE_READ_LOCK_WAIT_SEC)
            except Exception:  # noqa: BLE001 — the next scheduled poll counts the fill
                log.exception("place read: reading %s after a Place failed", venue)
                return
            if missed:
                log.warning("place read: %s not read after a Place (busy or failed); "
                            "the next scheduled poll counts it", ", ".join(missed))

        def _store_bfa_open_bets(self, open_wagers: list[dict]) -> None:
            try:
                store.upsert_bets(normalize_open_bets(open_wagers, utc_now_iso()), _now())
            except Exception:  # noqa: BLE001 — the ticket stands; the next BFA poll stores it
                log.exception("bfa teaser: storing the open bets read at confirmation failed")

        def do_DELETE(self) -> None:  # noqa: N802 — http.server's name
            if self._refused_foreign_host():
                return
            url = urlparse(self.path)
            if url.path == "/pins.json":
                self._delete_pin(url.query)
                return
            if url.path != "/crosswalk.json":
                self._send_json(404, {"error": f"no route for DELETE {url.path}"})
                return
            cleared = store.clear_crosswalk()
            self._send_json(200, {"ok": True, "cleared": cleared, "crosswalk": []})

        def _delete_pin(self, query: str) -> None:
            bet_ids = parse_qs(query).get("betId")
            if not bet_ids or not _non_empty_string(bet_ids[0]):
                self._send_json(400, {"error": "betId query parameter required"})
                return
            result = store.unpin_bet(bet_ids[0])
            self._send_json(200, {"ok": True, **result, "pins": store.load_pins(), "crosswalk": store.load_crosswalk()})

        # DNS-rebinding guard, first thing on every verb: a page whose name
        # resolves to 127.0.0.1 is same-origin with this server, so nothing
        # else here stops it reading every bet or clearing the crosswalk. Its
        # Host is still its own name. The port comes from the socket, not
        # config, so a service on a non-default port guards itself.
        # The taker's routes are the runner's alone, which calls on loopback:
        # the tailnet name (phone, any tailnet device) can never send a bot order.
        def _refused_non_loopback(self) -> bool:
            port = self.server.server_address[1]
            host = self.headers.get("Host")
            if host_allowed(host, port, ()):
                return False
            self._send_json(403, {"error": f"the taker's routes are loopback-only, got Host {host!r}"})
            return True

        def _refused_foreign_host(self) -> bool:
            port = self.server.server_address[1]
            host = self.headers.get("Host")
            if host_allowed(host, port, proxied_host_names):
                return False
            allowed = list(allowed_hosts(port, proxied_host_names))
            self._send_json(403, {"error": f"Host must be one of {allowed}, got {host!r}"})
            return True

        # The parsed JSON body, or a (status, error payload) tuple to send. The
        # Content-Type check is the write guard: without CORS headers here a
        # browser page can only reach this server with a "simple" request
        # (form or text/plain), never application/json, and the extension
        # page is exempt from CORS for its host permission.
        def _read_json_body(self) -> object | tuple[int, dict]:
            content_type = self.headers.get("Content-Type", "")
            if not content_type.split(";")[0].strip().lower() == "application/json":
                return 415, {"error": "Content-Type must be application/json"}
            try:
                length = int(self.headers.get("Content-Length", ""))
            except ValueError:
                return 411, {"error": "Content-Length required"}
            if length < 0 or length > MAX_BODY_BYTES:
                return 413, {"error": f"body must be at most {MAX_BODY_BYTES} bytes, got {length}"}
            try:
                return json.loads(self.rfile.read(length))
            except ValueError as error:
                return 400, {"error": f"body is not JSON: {error}"}

        def _send_json(self, status: int, payload: dict) -> None:
            self._send_bytes(status, json.dumps(payload).encode(), "application/json")

        def _send_bytes(self, status: int, body: bytes, content_type: str, extra_headers: dict | None = None) -> None:
            self.send_response(status)
            self.send_header("Content-Type", content_type)
            self.send_header("Content-Length", str(len(body)))
            self.send_header("Cache-Control", "no-store")
            for name, value in (extra_headers or {}).items():
                self.send_header(name, value)
            self.end_headers()
            self.wfile.write(body)

        def log_message(self, format: str, *args: object) -> None:  # noqa: A002
            log.debug("http %s", format % args)

    return BetsHandler


def teaser_placer_for(sources: list[Source]) -> bfa_teaser.BFATeaserPlacer | None:
    """The Place button's placer on the BFA source's session, or None when BFA is not read."""
    bfa_source = next((source for source in sources if isinstance(source, BFASource)), None)
    return bfa_teaser.BFATeaserPlacer(bfa_source) if bfa_source is not None else None


def open_taker_desk(store: TakerStore | None,
                    lock_check: Callable[[KalshiOrderClient, int], tuple[str | None, float | None]] = key_lock_problem,
                    ) -> tuple[TakerDesk | None, str]:
    """(the bot's order desk, "") or (None, why): without its own key, its
    subaccount, Kalshi's word that the key is locked to that subaccount, or its
    log there are no orders, and POST /taker/order.json answers 503 with that
    why. The lock check is two signed GETs (taker_orders.key_lock_problem)."""
    if not config.KALSHI_TAKER_API_KEY_ID or not config.KALSHI_TAKER_PRIVATE_KEY_PATH:
        reason = "no bot Kalshi key is configured (KALSHI_TAKER_API_KEY_ID / KALSHI_TAKER_PRIVATE_KEY_PATH)"
        log.info("taker: %s, the bot places no orders", reason)
        return None, reason
    if store is None:
        reason = f"the taker log {config.TAKER_DB_PATH} is not open (its ceiling lives there)"
        log.error("taker: %s, the bot places no orders", reason)
        return None, reason
    subaccount = config.KALSHI_TAKER_SUBACCOUNT
    if subaccount is None:
        reason = "no bot subaccount is configured (KALSHI_TAKER_SUBACCOUNT): the bot trades only in its own"
        log.error("taker: %s, the bot places no orders", reason)
        return None, reason
    if not Path(config.KALSHI_TAKER_PRIVATE_KEY_PATH).exists():
        reason = f"no bot key file at {config.KALSHI_TAKER_PRIVATE_KEY_PATH}"
        log.info("taker: %s yet, the bot places no orders", reason)
        return None, reason
    try:
        client = KalshiOrderClient(config.KALSHI_TAKER_API_KEY_ID, config.KALSHI_TAKER_PRIVATE_KEY_PATH,
                                   config.KALSHI_BASE_URL)
    except (OSError, ValueError) as error:
        reason = f"the bot key {config.KALSHI_TAKER_PRIVATE_KEY_PATH} could not be read: {error}"
        log.error("taker: %s, the bot places no orders", reason)
        return None, reason
    problem, balance = lock_check(client, subaccount)
    if problem is not None:
        reason = f"the bot key is not proven locked to subaccount {subaccount}: {problem}"
        log.error("taker: %s, the bot places no orders", reason)
        return None, reason
    shown_balance = "unknown" if balance is None else f"${balance:,.2f}"
    log.info("taker: ARMED on its own key, locked to Kalshi subaccount %d (available %s); stop file %s",
             subaccount, shown_balance, config.TAKER_STOP_PATH)
    return TakerDesk(client, store, config.TAKER_STOP_PATH, subaccount,
                     other_venue_unseen=lambda event_id, held_ids: novig_unseen_fill(store, event_id, held_ids)), ""


def novig_trading_client(key_id: str | None, key_path: str | None, variable: str) -> tuple[NovigOrderClient | None, str]:
    """(a client signing with one Novig trading key, "") or (None, why)."""
    if not key_id or not key_path:
        return None, f"no {variable} Novig key is configured"
    if not Path(key_path).exists():
        return None, f"no {variable} Novig key file at {key_path}"
    try:
        return NovigOrderClient(key_id, key_path, config.NOVIG_API_URL), ""
    except (OSError, ValueError, NovigError) as error:
        return None, f"the {variable} Novig key {key_path} could not be read: {type(error).__name__}"


def open_novig_taker_desk(store: TakerStore | None, holdings_status: Callable[[], dict[str, dict]],
                          check: Callable[[NovigOrderClient], tuple[str | None, float | None]] = account_problem,
                          ) -> tuple[NovigTakerDesk | None, str]:
    """(the bot's Novig desk, "") or (None, why): without its own trading key,
    Novig's word that it is a trading key on its own subaccount (it reads
    that balance and is refused the key listing), or its log, nothing is sent."""
    client, reason = novig_trading_client(config.NOVIG_TAKER_KEY_ID, config.NOVIG_TAKER_PRIVATE_KEY_PATH, "bot")
    if client is None:
        log.info("novig taker: %s, the bot places no Novig orders", reason)
        return None, reason
    if store is None:
        reason = f"the taker log {config.TAKER_DB_PATH} is not open (its ceiling lives there)"
        log.error("novig taker: %s, the bot places no Novig orders", reason)
        return None, reason
    problem, balance = check(client)
    if problem is not None:
        reason = f"the bot Novig key is not proven a trading key on its own subaccount: {problem}"
        log.error("novig taker: %s, the bot places no Novig orders", reason)
        return None, reason
    log.info("novig taker: ARMED on its own trading key %s (subaccount balance $%s); stop files %s, %s",
             config.NOVIG_TAKER_KEY_ID[:8], f"{balance:,.2f}", config.TAKER_STOP_PATH, config.NOVIG_TAKER_STOP_PATH)
    return NovigTakerDesk(client, store, (config.TAKER_STOP_PATH, config.NOVIG_TAKER_STOP_PATH), holdings_status), ""


def novig_subaccount_sources(taker_store: TakerStore | None, store: BetsStore) -> list[Source]:
    """The bot's and the Place button's Novig subaccounts, each read on its
    own trading key when one is configured. The bot's bets are named by the
    board lines its desk kept (taker.duckdb), Place's by the lines its orders
    bought (bets.duckdb::manual_orders), which also keep Place's settle times."""
    sources: list[Source] = []
    bot_client, reason = novig_trading_client(config.NOVIG_TAKER_KEY_ID, config.NOVIG_TAKER_PRIVATE_KEY_PATH, "bot")
    if bot_client is not None and taker_store is not None:
        sources.append(NovigSubaccountSource("bot", bot_client, taker_store.novig_order_contexts))
    elif bot_client is not None:
        log.error("novig_bot: the taker log is not open, so the bot subaccount's bets cannot be named; not read")
    else:
        log.info("novig_bot: %s, not read", reason)
    place_client, reason = novig_trading_client(config.NOVIG_PLACE_KEY_ID, config.NOVIG_PLACE_PRIVATE_KEY_PATH, "Place")
    if place_client is not None:
        sources.append(NovigSubaccountSource(
            "place", place_client, store.novig_place_order_contexts, note_settled=store.note_novig_place_settled,
            lost_orders=store.lost_novig_place_orders,
            adopt_order=lambda client_order_id, order_id: store.adopt_novig_place_order(client_order_id, order_id, _now()),
            unconfirmed_orders=store.unconfirmed_novig_place_orders,
            confirm_order=lambda client_order_id, outcome, filled, paid: store.confirm_novig_place_order(
                client_order_id, outcome, filled, paid, _now())))
    else:
        log.info("novig_place: %s, not read", reason)
    return sources


def open_novig_book_stream() -> NovigBookStream | None:
    """The live Novig book stream on the bot's Novig key, or None without one.
    Reads only: it runs whether or not the bot desk is armed or stopped."""
    client, reason = novig_trading_client(config.NOVIG_TAKER_KEY_ID, config.NOVIG_TAKER_PRIVATE_KEY_PATH, "bot")
    if client is None:
        log.info("novig stream: %s, no live Novig books and no keyed catalog reads (the runner's Novig checks fail closed)", reason)
        return None
    return NovigBookStream(client, stream_url(config.NOVIG_API_URL))


def open_novig_catalog_reads() -> NovigCatalogReads | None:
    """The runner's keyed Novig catalog reads on the bot's Novig key (its own
    client and sessions, apart from the desk's), or None without that key."""
    if not config.NOVIG_TAKER_KEY_ID or not config.NOVIG_TAKER_PRIVATE_KEY_PATH:
        log.info("novig reads: no bot Novig key is configured, the runner reads no Novig catalog")
        return None
    try:
        client = NovigOrderClient(config.NOVIG_TAKER_KEY_ID, config.NOVIG_TAKER_PRIVATE_KEY_PATH, config.NOVIG_API_URL,
                                  transport=pooled_transport())
    except (OSError, ValueError, NovigError) as error:
        log.error("novig reads: the bot Novig key %s could not be read (%s), the runner reads no Novig catalog",
                  config.NOVIG_TAKER_PRIVATE_KEY_PATH, type(error).__name__)
        return None
    log.info("novig reads: keyed catalog reads on the bot's key %s", config.NOVIG_TAKER_KEY_ID[:8])
    return NovigCatalogReads(client)


def kalshi_subaccounts_read(bot_subaccount: int | None) -> tuple[int, ...]:
    """The Kalshi subaccounts the Bets tab reads: the primary, plus the bot's
    when one is configured, so the bot's bets size Cal's and its own."""
    return (PRIMARY_SUBACCOUNT,) if bot_subaccount is None else (PRIMARY_SUBACCOUNT, bot_subaccount)


def open_kalshi_place_desk(store: BetsStore) -> tuple[ManualKalshiDesk | None, str]:
    """(the Place button's Kalshi desk on Cal's MAIN account, "") or (None, why).
    It signs only with its own trade key (KALSHI_PLACE_*): the read key the
    Kalshi source polls with is read-only, so it is never a fallback."""
    if not config.KALSHI_PLACE_API_KEY_ID or not config.KALSHI_PLACE_PRIVATE_KEY_PATH:
        reason = ("no Kalshi trade key for Place is configured (KALSHI_PLACE_API_KEY_ID / "
                  "KALSHI_PLACE_PRIVATE_KEY_PATH)")
        log.info("kalshi place: %s, the Place button is off", reason)
        return None, reason
    try:
        client = KalshiOrderClient(config.KALSHI_PLACE_API_KEY_ID, config.KALSHI_PLACE_PRIVATE_KEY_PATH,
                                   config.KALSHI_BASE_URL)
    except (OSError, ValueError) as error:
        reason = f"the Place key {config.KALSHI_PLACE_PRIVATE_KEY_PATH} could not be read: {error}"
        log.error("kalshi place: %s, the Place button is off", reason)
        return None, reason
    log.info("kalshi place: on, signing with Place key %s on the main account", config.KALSHI_PLACE_API_KEY_ID[:8])
    return ManualKalshiDesk(client, store, lambda line: read_runner_quote(config.RUNNER_URL, line)), ""


def open_novig_place_desk(store: BetsStore, source_names: list[str],
                          check: Callable[[NovigOrderClient], tuple[str | None, float | None]] = account_problem,
                          ) -> tuple[ManualNovigDesk | None, str]:
    """(the Place button's Novig desk on Place's own sub-account, "") or (None, why).
    It signs only with Place's own trading key (NOVIG_PLACE_*), proven by Novig
    a trading key on its own sub-account (it reads that balance and is refused
    the key listing), so the money in that sub-account is all it can spend.
    The stake counts Cal's app account and Place's sub-account, and the bot's
    when it is read, so each of those reads must be fresh before an order."""
    client, reason = novig_trading_client(config.NOVIG_PLACE_KEY_ID, config.NOVIG_PLACE_PRIVATE_KEY_PATH, "Place")
    if client is None:
        log.info("novig place: %s, the Novig Place button is off", reason)
        return None, reason
    missing = [name for name in NOVIG_PLACE_READ_SOURCES if name not in source_names]
    if missing:
        reason = f"the bets service does not read {', '.join(missing)}, so a stake could miss a Novig bet"
        log.error("novig place: %s, the Novig Place button is off", reason)
        return None, reason
    problem, balance = check(client)
    if problem is not None:
        reason = f"the Place Novig key is not proven a trading key on its own sub-account: {problem}"
        log.error("novig place: %s, the Novig Place button is off", reason)
        return None, reason
    read_sources = NOVIG_PLACE_READ_SOURCES + (("novig_bot",) if "novig_bot" in source_names else ())
    log.info("novig place: on, signing with Place key %s (sub-account balance $%s); stakes count %s",
             config.NOVIG_PLACE_KEY_ID[:8], f"{balance:,.2f}", ", ".join(read_sources))
    return ManualNovigDesk(client, store,
                           lambda line: read_runner_quote(config.RUNNER_URL, line, RUNNER_NOVIG_QUOTE_PATH,
                                                          RUNNER_NOVIG_QUOTE_TIMEOUT_SEC),
                           read_sources=read_sources), ""


def open_polymarket_place_desk(store: BetsStore, source_names: list[str],
                               client_factory: Callable[[str, str], PolymarketOrderClient] = PolymarketOrderClient,
                               ) -> tuple[ManualPolymarketDesk | None, str]:
    """(the Place button's Polymarket US desk on Cal's main account, "") or (None, why).
    It signs with the key the polymarket_us source reads with (POLYMARKET_US_*,
    Cal 2026-10-22), proven to sign by one balance read. The stake counts that
    source's read, so it must be read here."""
    if not config.POLYMARKET_US_KEY_ID or not config.POLYMARKET_US_SECRET_KEY:
        log.info("polymarket place: %s, the Polymarket US Place button is off", POLYMARKET_PLACE_OFF_REASON)
        return None, POLYMARKET_PLACE_OFF_REASON
    if "polymarket_us" not in source_names:
        reason = "the bets service does not read polymarket_us, so a stake could miss a Polymarket US bet"
        log.error("polymarket place: %s, the Polymarket US Place button is off", reason)
        return None, reason
    try:
        client = client_factory(config.POLYMARKET_US_KEY_ID, config.POLYMARKET_US_SECRET_KEY)
        status, reply = client.balance()
    except (OSError, PolymarketError) as error:
        reason = f"the POLYMARKET_US_* key could not be used: {error}"
        log.error("polymarket place: %s, the Polymarket US Place button is off", reason)
        return None, reason
    if status != 200:
        reason = f"Polymarket US refused the key's balance read: HTTP {status} {str(reply)[:200]}"
        log.error("polymarket place: %s, the Polymarket US Place button is off", reason)
        return None, reason
    log.info("polymarket place: on, signing with key %s on the main account", config.POLYMARKET_US_KEY_ID[:8])
    return ManualPolymarketDesk(client, store,
                                lambda line: read_runner_quote(config.RUNNER_URL, line, RUNNER_POLYMARKET_QUOTE_PATH,
                                                               RUNNER_POLYMARKET_QUOTE_TIMEOUT_SEC)), ""


def open_taker_store(db_path: Path) -> TakerStore | None:
    """The taker's decision log, or None (logged) when its file cannot be
    opened: the log is the bot's, so a bad path must not take Cal's bets down
    with it; the runner keeps its decisions queued and /taker_decisions.json 503s."""
    try:
        return TakerStore(db_path)
    except Exception as error:  # duckdb raises its own IOException family, not OSError
        log.error("taker log %s could not be opened, POST /taker_decisions.json will 503: %s", db_path, error)
        return None


def serve(sources: list[Source], store: BetsStore, host: str, port: int,
          taker_store: TakerStore | None = None) -> None:
    """Run the poll thread and the HTTP server until SIGINT/SIGTERM. The
    taker log is opened here unless main() opened it for the Novig sources."""
    stop = threading.Event()
    started_at = time.monotonic()
    source_names = [source.name for source in sources] + list(PUSHED_SOURCES)
    balance_book = BalanceBook()
    if taker_store is None:
        taker_store = open_taker_store(config.TAKER_DB_PATH)
    taker_desk, taker_desk_off_reason = open_taker_desk(taker_store)
    novig_taker_desk, novig_taker_desk_off_reason = open_novig_taker_desk(
        taker_store, lambda: source_status(store, source_names))
    kalshi_place_desk, kalshi_place_off_reason = open_kalshi_place_desk(store)
    novig_place_desk, novig_place_off_reason = open_novig_place_desk(store, source_names)
    polymarket_place_desk, polymarket_place_off_reason = open_polymarket_place_desk(store, source_names)
    novig_book_stream = open_novig_book_stream()
    novig_catalog_reads = open_novig_catalog_reads()
    source_reader = SourceReader(sources, store)
    server = ThreadingHTTPServer((host, port), make_handler(store, started_at, source_names,
                                                            teaser_placer=teaser_placer_for(sources),
                                                            taker_store=taker_store, taker_desk=taker_desk,
                                                            taker_desk_off_reason=taker_desk_off_reason,
                                                            kalshi_place_desk=kalshi_place_desk,
                                                            kalshi_place_off_reason=kalshi_place_off_reason,
                                                            novig_taker_desk=novig_taker_desk,
                                                            novig_taker_desk_off_reason=novig_taker_desk_off_reason,
                                                            novig_place_desk=novig_place_desk,
                                                            novig_place_off_reason=novig_place_off_reason,
                                                            polymarket_place_desk=polymarket_place_desk,
                                                            polymarket_place_off_reason=polymarket_place_off_reason,
                                                            novig_book_stream=novig_book_stream,
                                                            novig_catalog_reads=novig_catalog_reads,
                                                            balance_book=balance_book,
                                                            source_reader=source_reader))
    server.daemon_threads = True
    pollers = [threading.Thread(target=poll_loop, args=(lane, store, stop, balance_book, source_reader),
                                name=f"poll-{lane[0].name}" if len(lane) == 1 else "poll", daemon=True)
               for lane in poll_lanes(sources)]
    # The desk sweeps its own orders (resting cancels, final fills, settlement) whether or not the runner sends any.
    background_threads = [threading.Thread(target=desk.run_refresh_loop, args=(stop,), name=name, daemon=True)
                   for desk, name in ((taker_desk, "taker-desk"), (novig_taker_desk, "novig-taker-desk")) if desk is not None]
    if novig_book_stream is not None:
        background_threads.append(threading.Thread(target=novig_book_stream.run, args=(stop,), name="novig-stream", daemon=True))

    def request_stop(_signum, _frame) -> None:
        log.info("stopping")
        stop.set()
        threading.Thread(target=server.shutdown, daemon=True).start()

    signal.signal(signal.SIGINT, request_stop)
    signal.signal(signal.SIGTERM, request_stop)
    for poller in pollers:
        poller.start()
    for background_thread in background_threads:
        background_thread.start()
    log.info("bets service on http://%s:%d (sources: %s); phone page at / with /edges.json from %s; "
             "extra allowed hosts: %s", host, port, ", ".join(source_names), config.RUNNER_URL,
             ", ".join(config.EXTRA_ALLOWED_HOSTS) or "none")
    try:
        server.serve_forever()
    finally:
        server.server_close()
        stop.set()
        for poller in pollers:
            poller.join(timeout=5)
        for background_thread in background_threads:
            background_thread.join(timeout=5)
        store.close()
        if taker_store is not None:
            taker_store.close()


# Venues polled only when their credentials are present; each factory returns
# None (logging the fix) when they are not. check_sources.py reads this same
# list, so the login check can never test a different set than the service runs.
OPTIONAL_SOURCE_FACTORIES: tuple[tuple[str, Callable[[], Source | None]], ...] = (
    ("betonline", betonline_source_if_configured),
    ("novig", novig_source_if_connected),
    ("bfa", bfa_source_if_configured),
    ("wagerzon", wagerzon_source_if_configured),
    ("polymarket_us", polymarket_us_source_if_configured),
)


def main() -> None:
    setup_logging()
    store = BetsStore(config.DB_PATH, config.SOURCE_RUNS_RETENTION_DAYS)
    sources: list[Source] = [KalshiSource(known_tickers=store.load_kalshi_tickers,
                                          final_tickers=store.load_kalshi_final_tickers,
                                          subaccounts=kalshi_subaccounts_read(config.KALSHI_TAKER_SUBACCOUNT))]
    for _venue, factory in OPTIONAL_SOURCE_FACTORIES:
        optional = factory()
        if optional is not None:
            sources.append(optional)
    taker_store = open_taker_store(config.TAKER_DB_PATH)
    sources.extend(novig_subaccount_sources(taker_store, store))
    serve(sources, store, config.BIND_HOST, config.PORT, taker_store=taker_store)


if __name__ == "__main__":
    main()
