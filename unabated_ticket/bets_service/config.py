"""Config for the bets service. Loaded from .env files or the environment.

Reads, in this order (the first place a key is found wins):
  1. process environment
  2. unabated_ticket/bets_service/.env            (gitignored)
  3. <main checkout>/kalshi_draft/.env             (the bots' shared credentials)
  4. <main checkout>/bet_logger/.env               (the sheet scrapers' book logins: BFA, Wagerzon;
                                                     the Polymarket US API key)
Side effects: none — pure constants; auth_client.configure() happens in the
Kalshi source.
"""
import os
import re
from pathlib import Path

PKG_DIR = Path(__file__).parent
_RAW_ROOT = PKG_DIR.parent.parent
# Worktree-aware root (kalshi_rfi/config.py pattern): credentials, bet_logger's
# cookie file and kalshi_draft/auth.py live in the main checkout, not the
# worktree. Claude desktop puts worktrees under .claude/worktrees/<name>, which
# the ".worktrees" marker alone missed (no Kalshi credentials from a worktree).
WORKTREE_MARKERS = ("/.claude/worktrees/", "/.worktrees/")


def main_checkout_root(raw_root: Path) -> Path:
    path = f"{raw_root}/"
    for marker in WORKTREE_MARKERS:
        if marker in path:
            return Path(path.split(marker)[0])
    return raw_root


PROJECT_ROOT = main_checkout_root(_RAW_ROOT)


def _load_env(path: Path) -> dict[str, str]:
    env: dict[str, str] = {}
    if not path.exists():
        return env
    for raw in path.read_text().splitlines():
        line = raw.strip()
        if not line or line.startswith("#") or "=" not in line:
            continue
        key, value = line.split("=", 1)
        env[key.strip()] = value.strip().strip('"').strip("'")
    return env


_SERVICE_FILE_ENV = _load_env(PKG_DIR / ".env")
_DRAFT_FILE_ENV = _load_env(PROJECT_ROOT / "kalshi_draft" / ".env")
_BET_LOGGER_FILE_ENV = _load_env(PROJECT_ROOT / "bet_logger" / ".env")
_ENV_LOOKUP_ORDER = (os.environ, _SERVICE_FILE_ENV, _DRAFT_FILE_ENV, _BET_LOGGER_FILE_ENV)


def _get(key: str, default: str | None = None) -> str | None:
    for env in _ENV_LOOKUP_ORDER:
        if key in env:
            return env[key]
    return default


# Kalshi credentials (same names as every bot so one .env works everywhere).
KALSHI_API_KEY_ID = _get("KALSHI_API_KEY_ID")
KALSHI_PRIVATE_KEY_PATH = _get("KALSHI_PRIVATE_KEY_PATH")
# The Ticket's Place button's own trade key (manual_orders.py, Cal 2026-10-16):
# a key on his MAIN account with trade scope, made for Place alone. The read key
# above is read-only, so with this unset the Place button is off.
KALSHI_PLACE_API_KEY_ID = _get("KALSHI_PLACE_API_KEY_ID")
KALSHI_PLACE_PRIVATE_KEY_PATH = _get("KALSHI_PLACE_PRIVATE_KEY_PATH")
KALSHI_BASE_URL = _get("KALSHI_BASE_URL",
                       "https://api.elections.kalshi.com/trade-api/v2")

# Where the store lives. The server deploy (unabated_ticket/deploy/) points it
# at a writable data directory, since the repo is mounted read-only there.
DB_PATH = Path(_get("BETS_DB_PATH", str(PKG_DIR / "bets.duckdb")))
# The automated taker's decision log (taker_store.py), its own file so the
# bot's log never shares a write path with Cal's bet history.
TAKER_DB_PATH = Path(_get("TAKER_DB_PATH", str(PKG_DIR / "taker.duckdb")))
# The bot's own Kalshi key (taker_orders.py): made for the bot alone so it can
# be revoked alone, and locked to the bot's own subaccount. All three set = the
# bot places orders; any missing = it does not.
KALSHI_TAKER_API_KEY_ID = _get("KALSHI_TAKER_API_KEY_ID")
KALSHI_TAKER_PRIVATE_KEY_PATH = _get("KALSHI_TAKER_PRIVATE_KEY_PATH")
FIRST_NUMBERED_SUBACCOUNT = 1
LAST_NUMBERED_SUBACCOUNT = 63


def _subaccount_number(raw: str | None) -> int | None:
    """A numbered Kalshi subaccount (1-63), or None when unset. Raises on
    anything else: 0 is the primary account, which the bot must never trade."""
    if raw is None or not raw.strip():
        return None
    if raw.strip().isdigit() and FIRST_NUMBERED_SUBACCOUNT <= int(raw) <= LAST_NUMBERED_SUBACCOUNT:
        return int(raw)
    raise ValueError(f"KALSHI_TAKER_SUBACCOUNT: expected a subaccount number "
                     f"{FIRST_NUMBERED_SUBACCOUNT}-{LAST_NUMBERED_SUBACCOUNT}, got {raw!r}")


# The bot's Kalshi subaccount (Cal, 2026-10-16): the money in it is all the bot
# can lose, since its key is locked there. The Kalshi source reads it next to
# the primary account (sources/kalshi.py); the desk arms only with it set.
KALSHI_TAKER_SUBACCOUNT = _subaccount_number(_get("KALSHI_TAKER_SUBACCOUNT"))
# The stop switch: while this file exists the bot sends no order (`touch` it to stop, `rm` to resume).
TAKER_STOP_PATH = Path(_get("TAKER_STOP_PATH", str(PKG_DIR / "taker.stop")))

# Novig's v3 trading API (novig_orders.py; plan: automated_taker_plan.md § Novig).
NOVIG_API_URL = _get("NOVIG_API_URL", "https://api.novig.com")
# The bot's Novig subaccount, addressed by its one trading key (Cal, 2026-10-16:
# $500 in a subaccount of its own). A trading key can trade only there and can
# never move money. Both set = the Novig desk arms (novig_taker_desk.py).
NOVIG_TAKER_KEY_ID = _get("NOVIG_TAKER_KEY_ID")
NOVIG_TAKER_PRIVATE_KEY_PATH = _get("NOVIG_TAKER_PRIVATE_KEY_PATH")
# The Novig desk's own stop switch, next to the shared TAKER_STOP_PATH (which stops
# both bots): the Novig setup puts it in place so arming Novig never touches Kalshi's.
NOVIG_TAKER_STOP_PATH = Path(_get("NOVIG_TAKER_STOP_PATH", str(PKG_DIR / "novig_taker.stop")))
# The Ticket's Novig Place button's own subaccount (Cal, 2026-10-16), apart
# from the bot's. Read here so its bets count in sizing; Place itself is built on top.
NOVIG_PLACE_KEY_ID = _get("NOVIG_PLACE_KEY_ID")
NOVIG_PLACE_PRIVATE_KEY_PATH = _get("NOVIG_PLACE_PRIVATE_KEY_PATH")

# HTTP. Loopback only: the service has no auth and serves the user's own bets.
# Not configurable on purpose: on the server, `tailscale serve` on the host
# proxies the tailnet name to this loopback port, so nothing listens publicly.
BIND_HOST = "127.0.0.1"
PORT = int(_get("BETS_SERVICE_PORT", "8094"))

# A host name: dot-separated labels of lowercase letters, digits and inner
# hyphens. No scheme, port, path, wildcard or IPv6 literal.
_HOST_NAME_PATTERN = re.compile(r"[a-z0-9]([a-z0-9-]{0,61}[a-z0-9])?(\.[a-z0-9]([a-z0-9-]{0,61}[a-z0-9])?)*")


def parse_extra_allowed_hosts(raw: str | None) -> tuple[str, ...]:
    """BETS_EXTRA_ALLOWED_HOSTS ("a.ts.net,b.ts.net") as a tuple of host names.

    Empty or unset is (). Raises ValueError naming the first bad entry: the
    service must refuse to start rather than run with an allowlist it misread.
    """
    names: list[str] = []
    for entry in (raw or "").split(","):
        name = entry.strip()
        if not name:
            continue
        if name != name.lower():
            raise ValueError(f"BETS_EXTRA_ALLOWED_HOSTS: expected a lowercase host name, got {name!r}")
        if not _HOST_NAME_PATTERN.fullmatch(name) or len(name) > 253:
            raise ValueError(
                "BETS_EXTRA_ALLOWED_HOSTS: expected a bare host name like mlb-stack.tail1234.ts.net "
                f"(no scheme, port, path or wildcard), got {name!r}")
        if name not in names:
            names.append(name)
    return tuple(names)


# Extra names a request's Host may carry besides the loopback ones (#125
# allowlist). On the server, `tailscale serve` proxies
# https://<vm>.<tailnet>.ts.net to 127.0.0.1:PORT and forwards the browser's
# Host header unchanged, so the tailnet name must be listed here or the phone
# gets a 403. Each name is accepted bare (browsers omit :443) and with :443.
EXTRA_ALLOWED_HOSTS = parse_extra_allowed_hosts(_get("BETS_EXTRA_ALLOWED_HOSTS", ""))

# The server runner (unabated_ticket/server/runner.js) whose /edges.json the
# service proxies for the phone page, so the page has one origin. Short
# timeout: the runner answers from memory, so a slow reply means it is stuck.
RUNNER_URL = _get("UNABATED_RUNNER_URL", "http://127.0.0.1:8095").rstrip("/")
RUNNER_TIMEOUT_SEC = float(_get("UNABATED_RUNNER_TIMEOUT_SEC", "5"))

# /bets.json returns open bets plus settled/closed ones within this window.
RETENTION_DAYS = int(_get("BETS_RETENTION_DAYS", "30"))
# How long a source_runs row is kept. Only the latest run per source is ever
# read (the freshness the panel shows); the table appends ~2,000 rows/day and
# is scanned on every /bets.json and /health, so it is pruned to this window
# (0 = never prune).
SOURCE_RUNS_RETENTION_DAYS = int(_get("BETS_SOURCE_RUNS_RETENTION_DAYS", "7"))

# Kalshi cadence: fills + positions every poll; a full fills re-pull (the
# reconcile) once per RECONCILE_SEC; the fills poll overlaps the previous
# window by FILLS_OVERLAP_SEC (trade_id dedupe absorbs re-delivery), which must
# stay >= POLL_SEC. Every 10 s (2026-10-17) so a bot or Place fill is held within
# seconds: a poll is ~5 GETs (fills + positions per subaccount, a lookup per new
# ticker), under 1/s against Kalshi's published 20/s (Basic) read limit.
KALSHI_POLL_SEC = float(_get("BETS_KALSHI_POLL_SEC", "10"))
KALSHI_RECONCILE_SEC = float(_get("BETS_KALSHI_RECONCILE_SEC", "3600"))
KALSHI_FILLS_OVERLAP_SEC = 60
# Minimum gap between market/event lookups (kalshi_draft.auth precedent).
KALSHI_LOOKUP_GAP_SEC = 0.6

# Novig (issue #116): the service's own Auth0 refresh token, minted once by
# `python -m unabated_ticket.bets_service.sources.novig_auth connect` and
# rewritten on rotation. The source registers only when the file exists.
NOVIG_TOKEN_PATH = Path(_get("NOVIG_TOKEN_PATH", str(PKG_DIR / "novig_token.json")))
# Every 20 s (2026-10-17): Novig publishes no limit and 403s bursts, and each poll
# pages the whole active and settled portfolio, so it stays slower than Kalshi.
NOVIG_POLL_SEC = float(_get("BETS_NOVIG_POLL_SEC", "20"))
# The app's own REST base (novig.com bundle, 2026-09-22); the Portfolio feed
# lives under it. Its GraphQL is behind a query allowlist and is not used.
NOVIG_REST_URL = _get("NOVIG_REST_URL", "https://api.novig.us/nbx/v1")
# BetOnline (#115): the Keycloak refresh token lives in bet_logger's recon cookie
# file in the MAIN checkout (shared with scraper_betonline.py + its LaunchAgent);
# the paged report is polled every POLL_SEC over the last HISTORY_DAYS.
BETONLINE_COOKIES_PATH = Path(_get("BETS_BETONLINE_COOKIES_PATH",
                                   str(PROJECT_ROOT / "bet_logger" / "recon_betonline_cookies.json")))
BETONLINE_POLL_SEC = float(_get("BETS_BETONLINE_POLL_SEC", "300"))
BETONLINE_HISTORY_DAYS = int(_get("BETS_BETONLINE_HISTORY_DAYS", str(RETENTION_DAYS + 1)))

# BFA (2026-09-23): the account's own Keycloak password login — the sheet scraper's
# credentials in bet_logger/.env, no token file (see sources/bfa.py). Cadence and
# window are constants. The open bets every BFA_POLL_SEC: BFA is Buckeye, and a teaser
# placed there must reach the panel's Teasers tab before it is suggested again (this
# plus the panel's 30 s poll: ~90 s). The history — what settles a bet once it leaves
# the open list, one GET per 100 wagers — every BFA_HISTORY_POLL_SEC.
BFA_USERNAME = _get("BFA_USERNAME")
BFA_PASSWORD = _get("BFA_PASSWORD")
BFA_POLL_SEC = 60.0
BFA_HISTORY_POLL_SEC = 300.0
BFA_HISTORY_DAYS = RETENTION_DAYS + 1

# Wagerzon (2026-09-23): the C account. Its login has sat in the primary WAGERZON_*
# slot since 2026-06-26 (bet_logger/scraper_wagerzon.py); WAGERZONC_* wins when set.
# Six Mon-Sun weeks of HistoryHelper cover the 30-day window with margin.
WAGERZON_USERNAME = _get("WAGERZONC_USERNAME") or _get("WAGERZON_USERNAME")
WAGERZON_PASSWORD = _get("WAGERZONC_PASSWORD") or _get("WAGERZON_PASSWORD")
WAGERZON_POLL_SEC = 300.0
WAGERZON_HISTORY_WEEKS = 6

# Polymarket US (2026-09-23): the CFTC app's own API key, created by Cal at
# polymarket.us/developer and kept in bet_logger/.env (see sources/polymarket_us.py).
# Every 20 s (2026-10-17; no published limit, so slower than Kalshi); activities are read back one day past
# the retention window, and further while an open or just-settled position's fills are older.
POLYMARKET_US_KEY_ID = _get("POLYMARKET_US_KEY_ID")
POLYMARKET_US_SECRET_KEY = _get("POLYMARKET_US_SECRET_KEY")
POLYMARKET_US_POLL_SEC = 20.0
POLYMARKET_US_HISTORY_DAYS = RETENTION_DAYS + 1
# The Ticket's Place button at Polymarket US (polymarket_place.py, Cal 2026-10-21):
# its OWN key on Cal's main account, made for Place alone at polymarket.us/developer;
# the read key above is never a fallback. The secret (the base64 "Secret Key" shown
# once) lives in a FILE (scp'd to the VM's data dir), never in .env. Unset = Place off.
POLYMARKET_PLACE_KEY_ID = _get("POLYMARKET_PLACE_KEY_ID")
POLYMARKET_PLACE_SECRET_PATH = _get("POLYMARKET_PLACE_SECRET_PATH")

# Bet105 (2026-09-29): no credentials and no poll here. Cloudflare challenges any
# request that is not the browser session's own, so the extension reads the account
# from Cal's logged-in Chrome (extension/bet105.js, every 5 min while the panel is
# open) and POSTs it to /bet105.json; sources/bet105.py parses it.

# DraftKings (2026-10-10): the same — no credentials and no poll here. Akamai's bot
# checks and a login code stop anything but Cal's own Chrome, so the extension reads
# My Bets' socket there (extension/draftkings.js) and POSTs to /draftkings.json;
# sources/draftkings.py parses it.

# Logging
LOG_PATH = Path(_get("BETS_SERVICE_LOG_PATH", str(PKG_DIR / "bets_service.log")))
LOG_LEVEL = _get("BETS_SERVICE_LOG_LEVEL", "INFO")
LOG_ROTATE_MAX_BYTES = int(_get("BETS_SERVICE_LOG_ROTATE_MAX_BYTES", str(10 * 1024 * 1024)))
LOG_ROTATE_BACKUPS = int(_get("BETS_SERVICE_LOG_ROTATE_BACKUPS", "3"))
# Also log to stderr: "1" yes, "0" no, unset = only when stderr is a terminal.
# The container sets 1 so `docker compose logs` shows the service.
_LOG_CONSOLE_RAW = _get("BETS_SERVICE_LOG_CONSOLE", "")
if _LOG_CONSOLE_RAW not in ("", "0", "1"):
    raise ValueError(f"BETS_SERVICE_LOG_CONSOLE: expected 0, 1 or unset, got {_LOG_CONSOLE_RAW!r}")
LOG_CONSOLE = None if _LOG_CONSOLE_RAW == "" else _LOG_CONSOLE_RAW == "1"
