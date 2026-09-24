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
DB_PATH = PKG_DIR / "bets.duckdb"


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
KALSHI_BASE_URL = _get("KALSHI_BASE_URL",
                       "https://api.elections.kalshi.com/trade-api/v2")

# HTTP. Loopback only: the service has no auth and serves the user's own bets.
BIND_HOST = "127.0.0.1"
PORT = int(_get("BETS_SERVICE_PORT", "8094"))

# /bets.json returns open bets plus settled/closed ones within this window.
RETENTION_DAYS = int(_get("BETS_RETENTION_DAYS", "30"))
# How long a source_runs row is kept. Only the latest run per source is ever
# read (the freshness the panel shows); the table appends ~2,000 rows/day and
# is scanned on every /bets.json and /health, so it is pruned to this window
# (0 = never prune).
SOURCE_RUNS_RETENTION_DAYS = int(_get("BETS_SOURCE_RUNS_RETENTION_DAYS", "7"))

# Kalshi cadence: fills + positions every poll; a full fills re-pull (the
# reconcile) once per RECONCILE_SEC; the fills poll overlaps the previous
# window by FILLS_OVERLAP_SEC (trade_id dedupe absorbs re-delivery).
KALSHI_POLL_SEC = float(_get("BETS_KALSHI_POLL_SEC", "60"))
KALSHI_RECONCILE_SEC = float(_get("BETS_KALSHI_RECONCILE_SEC", "3600"))
KALSHI_FILLS_OVERLAP_SEC = 60
# Minimum gap between market/event lookups (kalshi_draft.auth precedent).
KALSHI_LOOKUP_GAP_SEC = 0.6

# Novig (issue #116): the service's own Auth0 refresh token, minted once by
# `python -m unabated_ticket.bets_service.sources.novig_auth connect` and
# rewritten on rotation. The source registers only when the file exists.
NOVIG_TOKEN_PATH = Path(_get("NOVIG_TOKEN_PATH", str(PKG_DIR / "novig_token.json")))
NOVIG_POLL_SEC = float(_get("BETS_NOVIG_POLL_SEC", "60"))
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
# window are constants: the history is one GET per 100 wagers.
BFA_USERNAME = _get("BFA_USERNAME")
BFA_PASSWORD = _get("BFA_PASSWORD")
BFA_POLL_SEC = 300.0
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
# An exchange like Kalshi, so the same cadence; activities are read back one day past
# the retention window, and further while an open or just-settled position's fills are older.
POLYMARKET_US_KEY_ID = _get("POLYMARKET_US_KEY_ID")
POLYMARKET_US_SECRET_KEY = _get("POLYMARKET_US_SECRET_KEY")
POLYMARKET_US_POLL_SEC = 60.0
POLYMARKET_US_HISTORY_DAYS = RETENTION_DAYS + 1

# Logging
LOG_PATH = Path(_get("BETS_SERVICE_LOG_PATH", str(PKG_DIR / "bets_service.log")))
LOG_LEVEL = _get("BETS_SERVICE_LOG_LEVEL", "INFO")
LOG_ROTATE_MAX_BYTES = int(_get("BETS_SERVICE_LOG_ROTATE_MAX_BYTES", str(10 * 1024 * 1024)))
LOG_ROTATE_BACKUPS = int(_get("BETS_SERVICE_LOG_ROTATE_BACKUPS", "3"))
