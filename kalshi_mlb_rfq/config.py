"""Config knobs for the Kalshi MLB RFQ bot. Loaded from .env (or environment)."""

import os
from pathlib import Path

PKG_DIR = Path(__file__).parent
_RAW_PROJECT_ROOT = PKG_DIR.parent
# Worktree-aware: when running from .claude/worktrees/<name>/kalshi_mlb_rfq/,
# resolve PROJECT_ROOT to the main repo so we find Answer Keys/, mlb_sgp/, etc.
# Matches mlb_sgp/db.py:19-20 pattern.
if ".worktrees" in str(_RAW_PROJECT_ROOT):
    PROJECT_ROOT = Path(str(_RAW_PROJECT_ROOT).split(".worktrees")[0].rstrip("/"))
else:
    PROJECT_ROOT = _RAW_PROJECT_ROOT
DB_PATH = PKG_DIR / "kalshi_mlb_rfq.duckdb"
LOG_PATH = PKG_DIR / "bot.log"
KILL_FILE = PKG_DIR / ".kill"


def _load_env(path: Path) -> dict[str, str]:
    """Parse a .env file; returns a dict of key→value."""
    env: dict[str, str] = {}
    if not path.exists():
        return env
    for raw in path.read_text().splitlines():
        line = raw.strip()
        if not line or line.startswith("#") or "=" not in line:
            continue
        k, v = line.split("=", 1)
        env[k.strip()] = v.strip().strip('"').strip("'")
    return env


_FILE_ENV = _load_env(PKG_DIR / ".env")


def _get(key: str, default: str | None = None) -> str | None:
    return os.environ.get(key, _FILE_ENV.get(key, default))


# Credentials
KALSHI_API_KEY_ID = _get("KALSHI_API_KEY_ID")
KALSHI_PRIVATE_KEY_PATH = _get("KALSHI_PRIVATE_KEY_PATH")
KALSHI_USER_ID = _get("KALSHI_USER_ID")
KALSHI_BASE_URL = _get(
    "KALSHI_BASE_URL", "https://api.elections.kalshi.com/trade-api/v2"
)

# MVE
MVE_COLLECTION_TICKER = _get("MVE_COLLECTION_TICKER", "KXMVECROSSCATEGORY-R")

# Sizing
BANKROLL = float(_get("BANKROLL", "1000.0"))
KELLY_FRACTION = float(_get("KELLY_FRACTION", "0.25"))
# Create-time Kelly sizing: we estimate the maker's yes_ask before any quote
# arrives by solving for the worst price that would still pass MIN_EV_PCT
# after fees (binary search in main._worst_acceptable_yes_ask). Kelly at
# that worst price = the smallest count we'd ever defensibly take, locked
# into the RFQ's `contracts` field. If the maker quotes better, we still get
# that count at a better price (under-sized vs ideal but safe). If worse,
# the accept-time MIN_EV_PCT gate declines.
#
# KELLY_CREATE_EV_FLOOR_PCT lets you target a different EV cliff for sizing
# than for acceptance (e.g., MIN_EV_PCT=5% accept threshold, but size as if
# the cliff were 8% → ~30% smaller Kelly count, more conservative).
# Default = MIN_EV_PCT so sizing and acceptance share one cliff.
KELLY_CREATE_EV_FLOOR_PCT = float(
    _get("KELLY_CREATE_EV_FLOOR_PCT", _get("MIN_EV_PCT", "0.05"))
)

# Per-accept gates
MIN_EV_PCT = float(_get("MIN_EV_PCT", "0.05"))
MIN_FAIR_PROB = float(_get("MIN_FAIR_PROB", "0.05"))
MAX_FAIR_PROB = float(_get("MAX_FAIR_PROB", "0.95"))
MAX_GAME_EXPOSURE_PCT = float(_get("MAX_GAME_EXPOSURE_PCT", "0.10"))
DAILY_EXPOSURE_CAP_USD = float(_get("DAILY_EXPOSURE_CAP_USD", "200.0"))
MAX_PREDICTION_STALENESS_SEC = int(_get("MAX_PREDICTION_STALENESS_SEC", "600"))
MAX_BOOK_STALENESS_SEC = int(_get("MAX_BOOK_STALENESS_SEC", "60"))
COMBO_COOLDOWN_SEC = int(_get("COMBO_COOLDOWN_SEC", "30"))
POSITIONS_HEALTH_RETRIES = int(_get("POSITIONS_HEALTH_RETRIES", "2"))
MIN_FILL_RATIO = float(_get("MIN_FILL_RATIO", "0.50"))
FILL_RATIO_WINDOW = int(_get("FILL_RATIO_WINDOW", "50"))
TIPOFF_CANCEL_MIN = int(_get("TIPOFF_CANCEL_MIN", "5"))

# Loops
RFQ_REFRESH_SEC = int(_get("RFQ_REFRESH_SEC", "30"))
QUOTE_POLL_SEC = int(_get("QUOTE_POLL_SEC", "2"))
RISK_SWEEP_SEC = int(_get("RISK_SWEEP_SEC", "10"))
PIPELINE_REFRESH_SEC = int(_get("PIPELINE_REFRESH_SEC", "600"))
MAX_LIVE_RFQS = int(_get("MAX_LIVE_RFQS", "80"))

# Notifications
NOTIFY_WEBHOOK_URL = _get("NOTIFY_WEBHOOK_URL")

# Vig fallbacks
DK_VIG_FALLBACK = float(_get("DK_VIG_FALLBACK", "0.125"))
FD_VIG_FALLBACK = float(_get("FD_VIG_FALLBACK", "0.18"))
PX_VIG_FALLBACK = float(_get("PX_VIG_FALLBACK", "0.05"))
NOVIG_VIG_FALLBACK = float(_get("NOVIG_VIG_FALLBACK", "0.05"))

# Source path for the answer-key DB
ANSWER_KEY_DB = PROJECT_ROOT / "Answer Keys" / "mlb_mm.duckdb"

# SGP scraper cadence (bot-driven independent of dashboard refresh)
SGP_REFRESH_SEC = int(_get("SGP_REFRESH_SEC", "60"))
SGP_SCRAPER_TIMEOUT_SEC = int(_get("SGP_SCRAPER_TIMEOUT_SEC", "90"))
SGP_MIN_INTERVAL_SEC = int(_get("SGP_MIN_INTERVAL_SEC", "30"))

# Bot market DB (sibling to bot state DB, holds mlb_target_lines + mlb_sgp_odds)
BOT_MARKET_DB = Path(_get("BOT_MARKET_DB",
                          str(PKG_DIR / "kalshi_mlb_rfq_market.duckdb")))

# Path to mlb_sgp directory (for subprocess scraper invocation).
# Resolved through PROJECT_ROOT so worktree contexts use the main repo's venv.
MLB_SGP_DIR = Path(_get("MLB_SGP_DIR", str(PROJECT_ROOT / "mlb_sgp")))

# Cross-book book-count gate (drop candidate if fewer than N books priced).
MIN_BOOK_COUNT_FOR_BLEND = int(_get("MIN_BOOK_COUNT_FOR_BLEND", "2"))

# Logging (operational log rotation)
LOG_LEVEL = _get("LOG_LEVEL", "INFO")
LOG_MAX_BYTES = int(_get("LOG_MAX_BYTES", str(50 * 1024 * 1024)))  # 50 MB
LOG_BACKUP_COUNT = int(_get("LOG_BACKUP_COUNT", "5"))
