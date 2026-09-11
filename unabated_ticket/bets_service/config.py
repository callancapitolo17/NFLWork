"""Config for the bets service. Loaded from .env files or the environment.

Reads, in this order (the first place a key is found wins):
  1. process environment
  2. unabated_ticket/bets_service/.env            (gitignored)
  3. <main checkout>/kalshi_draft/.env             (the bots' shared credentials)
Side effects: none — pure constants; auth_client.configure() happens in the
Kalshi source.
"""
import os
from pathlib import Path

PKG_DIR = Path(__file__).parent
_RAW_ROOT = PKG_DIR.parent.parent
# Worktree-aware root (kalshi_rfi/config.py pattern): credentials and
# kalshi_draft/auth.py live in the main checkout, not the worktree.
PROJECT_ROOT = (Path(str(_RAW_ROOT).split(".worktrees")[0].rstrip("/"))
                if ".worktrees" in str(_RAW_ROOT) else _RAW_ROOT)
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


def _get(key: str, default: str | None = None) -> str | None:
    return os.environ.get(
        key, _SERVICE_FILE_ENV.get(key, _DRAFT_FILE_ENV.get(key, default)))


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

# Kalshi cadence: fills + positions every poll; a full fills re-pull (the
# reconcile) once per RECONCILE_SEC; the fills poll overlaps the previous
# window by FILLS_OVERLAP_SEC (trade_id dedupe absorbs re-delivery).
KALSHI_POLL_SEC = float(_get("BETS_KALSHI_POLL_SEC", "60"))
KALSHI_RECONCILE_SEC = float(_get("BETS_KALSHI_RECONCILE_SEC", "3600"))
KALSHI_FILLS_OVERLAP_SEC = 60
# Minimum gap between market/event lookups (kalshi_draft.auth precedent).
KALSHI_LOOKUP_GAP_SEC = 0.6

# Logging
LOG_PATH = Path(_get("BETS_SERVICE_LOG_PATH", str(PKG_DIR / "bets_service.log")))
LOG_LEVEL = _get("BETS_SERVICE_LOG_LEVEL", "INFO")
LOG_ROTATE_MAX_BYTES = int(_get("BETS_SERVICE_LOG_ROTATE_MAX_BYTES", str(10 * 1024 * 1024)))
LOG_ROTATE_BACKUPS = int(_get("BETS_SERVICE_LOG_ROTATE_BACKUPS", "3"))
