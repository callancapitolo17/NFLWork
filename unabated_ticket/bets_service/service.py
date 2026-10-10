"""Entry point: the local bets service the Unabated Ticket panel polls (#114).

    python -m unabated_ticket.bets_service.service      (or bets_service/run.sh)

Inputs:  each registered Source (sources/kalshi.py; sources/betonline.py when
         its cookie file exists; sources/novig.py when its token file exists;
         sources/bfa.py and sources/wagerzon.py when their logins are configured;
         sources/polymarket_us.py when its API key is configured)
         on its own poll_sec; plus the extension's pushes for Bet105 (below).
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
                                    (games still to start)}
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
           POST /closing_fairs.json body {rows: [{betId, lineKey, points, fairAmerican,
                                    fairObservedAt, eventStart}]} -> {ok, saved}
                                    (the server runner's pre-start read of each open bet's
                                    line, server/closefair.js; a row older than the stored
                                    one is ignored; same Content-Type guard)
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
           POST /bet105.json        body {fetchedAt, feeds: {prematch: [betGroup], live:
                                    [betGroup]}, settled: [wager]} -> {ok, count, settled,
                                    closed, skipped}, or {error} -> {ok, recorded: "error"}
                                    (the one PUSHED source: the extension reads Bet105 from
                                    Cal's own Chrome because Cloudflare challenges anything
                                    else — sources/bet105.py parses; a graded wager settles
                                    its bet with the venue's result and settle time, an open
                                    bet neither read lists is closed with no result, and the
                                    push is logged as that source's run)
           GET /settings.json       {settings: {bankroll, multiplier, leagues, periods, betTypes,
                                    bookMode, bookIds, minEdgePct, minStake, maxLineAgeHours,
                                    minLiquidityToWin, includeAlts, sortBy, groupByMarket},
                                    updatedAt} — the Edges settings the server
                                    runner (unabated_ticket/server/runner.js) reads every cycle;
                                    null = the panel's default (extension/edgerows.js)
           POST /place_teaser.json  body {stake, legs: [4 x {league, betType, side, rotation, points,
                                    eventStart, label}]} -> {ok, status: placed | refused |
                                    unconfirmed, message, ticket?} (the Teasers tab's Place
                                    button: ONE 4-team 6-point teaser at BFA, real money —
                                    bfa_teaser.py checks every number on BFA's board first and
                                    confirms off the open bets; 503 when no BFA account is read)
           PUT /settings.json       body {settings: {<any of those fields>: value or null}} ->
                                    {ok, settings, updatedAt}; a field left out keeps its value,
                                    null resets it to the default, an unknown field or a bad
                                    value is a 400 naming it (same Content-Type guard as the
                                    POSTs: a web page cannot send JSON cross-origin)
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
insert-only, a bet that has a saved fair keeps it; UPSERTs
bets.duckdb::bet_closing_fairs on POST /closing_fairs.json (newest observation
wins); INSERTs into / DELETEs from bets.duckdb::bet_exclusions on POST
/exclusions.json; INSERTs into / DELETEs from bets.duckdb::bet_dismissals on
POST /dismissals.json; UPSERTs / DELETEs bets.duckdb::teaser_blocks on POST
/teaser_blocks.json; rotating log at
bets_service.log. A poll that raises writes a failed source_runs row and
leaves `bets` untouched — a dark source never blanks the list. PUT
/settings.json UPSERTs the one row of bets.duckdb::edge_settings. POST
/place_teaser.json places a wager at BFA (money, not disk) and UPSERTs the BFA open
bets it read while confirming into bets.duckdb::bets, so the panel's next poll has the
ticket.
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
from unabated_ticket.bets_service.log_setup import setup_logging
from unabated_ticket.bets_service.normalize import parse_iso_ms, utc_now_iso
from unabated_ticket.bets_service.sources import Source, bet105
from unabated_ticket.bets_service.sources.betonline import source_if_configured as betonline_source_if_configured
from unabated_ticket.bets_service.sources.bfa import BFASource, normalize_open_bets
from unabated_ticket.bets_service.sources.bfa import source_if_configured as bfa_source_if_configured
from unabated_ticket.bets_service.sources.kalshi import KalshiSource
from unabated_ticket.bets_service.sources.novig import source_if_connected as novig_source_if_connected
from unabated_ticket.bets_service.sources.polymarket_us import source_if_configured as polymarket_us_source_if_configured
from unabated_ticket.bets_service.sources.wagerzon import source_if_configured as wagerzon_source_if_configured
from unabated_ticket.bets_service.store import EDGE_SETTINGS_FIELDS, BetsStore

log = logging.getLogger(__name__)

POLL_TICK_SEC = 1.0
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
# A fill-fair POST carries the bets a capture pass decided: a handful.
MAX_FILL_FAIR_ROWS_PER_POST = 1000
# Cal's ask (2026-10-05): only these books carry bets that are not his.
EXCLUDABLE_VENUES = ("bfa", "wagerzon")
# A parlay's legs go in one request; a dozen is already a long ticket.
MAX_EXCLUSION_IDS_PER_POST = 50
# The panel's first sync sends every bet it had dismissed locally: a handful.
MAX_DISMISSAL_IDS_PER_POST = 200
# teaser.js marketKeyOf: "<Unabated eventId>:bt2" (spread) or ":bt3" (total).
TEASER_MARKET_KEY_PATTERN = re.compile(r"^[A-Za-z0-9_-]{1,64}:bt[23]$")
# American odds run from -100 down and +100 up; the gap between is no price.
MIN_AMERICAN_MAGNITUDE = 100
# The Edges settings' allowed values, as the panel's inputs enforce them
# (panel.js readSettingInputs / readEdgeSettingInputs). feed.PERIODS ids run
# 1 (full game) to 7 (4Q); feed.BET_TYPES are 1 moneyline, 2 spread, 3 total.
SETTINGS_PERIOD_IDS = range(1, 8)
SETTINGS_BET_TYPE_IDS = (1, 2, 3)
SETTINGS_BOOK_MODES = ("default", "all", "custom")
SETTINGS_SORT_KEYS = ("edge", "stake", "start", "exposure")
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
    "kelly.js", "feed.js", "teams.js", "bets.js", "ladder.js", "condkelly.js", "betsview.js",
    "edgemove.js", "fillfair.js", "tailflex.js", "edgerows.js", "teaser.js", "teaserview.js",
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

# Sources with no poll here: the extension POSTs their records (/bet105.json).
# Listed so the panel reads "no completed poll yet" before the first push, not
# "no source configured".
PUSHED_SOURCES = (bet105.VENUE,)


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


def poll_loop(sources: list[Source], store: BetsStore, stop: threading.Event) -> None:
    """Run each source on its own cadence until `stop`. A store write that
    raises (disk full, a locked file) is logged and retried next poll — the
    thread must outlive it, or the HTTP side would serve ageing data with no
    failed run to show for it."""
    next_due = {source.name: 0.0 for source in sources}
    while not stop.is_set():
        now = time.monotonic()
        for source in sources:
            if now >= next_due[source.name]:
                try:
                    run_source_once(source, store)
                except Exception:  # noqa: BLE001 — the poll thread is the service's heartbeat
                    log.exception("store write failed for source %s; retrying next poll", source.name)
                next_due[source.name] = time.monotonic() + source.poll_sec
        stop.wait(POLL_TICK_SEC)


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


def bets_payload(store: BetsStore, days: int, source_names: list[str] = ()) -> dict:
    now = _now()
    return {"generatedAt": _iso(now), "sources": source_status(store, list(source_names)),
            "bets": store.load_bets(days, now), "crosswalk": store.load_crosswalk(), "pins": store.load_pins(),
            "fillFairs": store.load_fill_fairs(days, now), "closingFairs": store.load_closing_fairs(days, now),
            "exclusions": store.load_exclusions(), "dismissals": store.load_dismissals(),
            "teaserBlocks": store.load_teaser_blocks(now)}


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


def validate_fill_fair_rows(body: object) -> list[dict] | str:
    """The rows of a POST /fill_fairs.json body, or an error naming the first
    bad row. Each row: betId and lineKey non-empty strings; points a number or
    null (a moneyline); fairAmerican a whole American price (Unabated's bacr
    always is); fairObservedAt and placedAt ISO times, the fair observed at or
    before the placement — a saved fair is permanent, so a row that cannot be
    the fill's is refused rather than stored."""
    if not isinstance(body, dict) or not isinstance(body.get("rows"), list):
        return f"body must be an object with a `rows` array, got {type(body).__name__}"
    raw_rows = body["rows"]
    if len(raw_rows) > MAX_FILL_FAIR_ROWS_PER_POST:
        return f"at most {MAX_FILL_FAIR_ROWS_PER_POST} rows per request, got {len(raw_rows)}"
    rows: list[dict] = []
    for index, raw in enumerate(raw_rows):
        row = _validate_fill_fair_row(index, raw)
        if isinstance(row, str):
            return row
        rows.append(row)
    return rows


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
    if not isinstance(body, dict) or not isinstance(body.get("rows"), list):
        return f"body must be an object with a `rows` array, got {type(body).__name__}"
    raw_rows = body["rows"]
    if len(raw_rows) > MAX_FILL_FAIR_ROWS_PER_POST:
        return f"at most {MAX_FILL_FAIR_ROWS_PER_POST} rows per request, got {len(raw_rows)}"
    rows: list[dict] = []
    for index, raw in enumerate(raw_rows):
        row = _validate_closing_fair_row(index, raw)
        if isinstance(row, str):
            return row
        rows.append(row)
    return rows


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
    if field in ("minEdgePct", "minStake", "minLiquidityToWin"):
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
RUNNER_ROUTES = ("/edges.json", "/scenarios.json", "/teasers.json", "/attach.json")
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


def health_payload(store: BetsStore, started_at: float, source_names: list[str] = ()) -> dict:
    return {"ok": True, "generatedAt": _iso(_now()), "uptimeSec": int(time.monotonic() - started_at),
            "sources": source_status(store, list(source_names))}


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
                 teaser_placer: bfa_teaser.BFATeaserPlacer | None = None) -> type[BaseHTTPRequestHandler]:
    """The request handler class. `runner_url` / `runner_timeout_sec` /
    `extra_hosts` default to config.RUNNER_URL / RUNNER_TIMEOUT_SEC /
    EXTRA_ALLOWED_HOSTS (the tests pass their own). `teaser_placer` is the
    Place button's (None: no BFA account is read, POST /place_teaser.json 503s)."""
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
                self._send_json(200, health_payload(store, started_at, names))
                return
            if url.path == "/settings.json":
                self._send_json(200, store.load_edge_settings())
                return
            if url.path == "/bets.json":
                days = parse_days(url.query)
                if isinstance(days, str):
                    self._send_json(400, {"error": days})
                    return
                self._send_json(200, bets_payload(store, days, names))
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
            if url.path not in ("/crosswalk.json", "/pins.json", "/fill_fairs.json", "/closing_fairs.json",
                                "/exclusions.json", "/dismissals.json", "/teaser_blocks.json", "/bet105.json",
                                "/place_teaser.json"):
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
            if url.path == "/closing_fairs.json":
                self._save_closing_fairs(body)
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
            if url.path == "/bet105.json":
                self._push_bet105(body)
                return
            if url.path == "/place_teaser.json":
                self._place_teaser(body)
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

        # POST /bet105.json: the extension's read of the account, as one source
        # run — a complete push UPSERTs its records (a graded wager over its open
        # record) and closes the open ones neither read lists; an error push is a
        # failed run and the records stand.
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


def serve(sources: list[Source], store: BetsStore, host: str, port: int) -> None:
    """Run the poll thread and the HTTP server until SIGINT/SIGTERM."""
    stop = threading.Event()
    started_at = time.monotonic()
    source_names = [source.name for source in sources] + list(PUSHED_SOURCES)
    server = ThreadingHTTPServer((host, port), make_handler(store, started_at, source_names,
                                                            teaser_placer=teaser_placer_for(sources)))
    server.daemon_threads = True
    poller = threading.Thread(target=poll_loop, args=(sources, store, stop), name="poll", daemon=True)

    def request_stop(_signum, _frame) -> None:
        log.info("stopping")
        stop.set()
        threading.Thread(target=server.shutdown, daemon=True).start()

    signal.signal(signal.SIGINT, request_stop)
    signal.signal(signal.SIGTERM, request_stop)
    poller.start()
    log.info("bets service on http://%s:%d (sources: %s); phone page at / with /edges.json from %s; "
             "extra allowed hosts: %s", host, port, ", ".join(source_names), config.RUNNER_URL,
             ", ".join(config.EXTRA_ALLOWED_HOSTS) or "none")
    try:
        server.serve_forever()
    finally:
        server.server_close()
        stop.set()
        poller.join(timeout=5)
        store.close()


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
    sources: list[Source] = [KalshiSource(known_tickers=store.load_kalshi_tickers)]
    for _venue, factory in OPTIONAL_SOURCE_FACTORIES:
        optional = factory()
        if optional is not None:
            sources.append(optional)
    serve(sources, store, config.BIND_HOST, config.PORT)


if __name__ == "__main__":
    main()
