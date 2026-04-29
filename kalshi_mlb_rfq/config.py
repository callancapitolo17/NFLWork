"""Config knobs for the Kalshi MLB RFQ bot. Loaded from .env (or environment)."""

import os
from pathlib import Path

PKG_DIR = Path(__file__).parent
PROJECT_ROOT = PKG_DIR.parent
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

# Per-accept gates
MIN_EV_PCT = float(_get("MIN_EV_PCT", "0.05"))
MAX_QUOTE_DEVIATION = float(_get("MAX_QUOTE_DEVIATION", "0.15"))
MIN_FAIR_PROB = float(_get("MIN_FAIR_PROB", "0.05"))
MAX_FAIR_PROB = float(_get("MAX_FAIR_PROB", "0.95"))
MAX_GAME_EXPOSURE_PCT = float(_get("MAX_GAME_EXPOSURE_PCT", "0.10"))
DAILY_EXPOSURE_CAP_USD = float(_get("DAILY_EXPOSURE_CAP_USD", "200.0"))
LINE_MOVE_THRESHOLD = float(_get("LINE_MOVE_THRESHOLD", "0.5"))
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
ANSWER_KEY_DB = PROJECT_ROOT / "Answer Keys" / "mlb.duckdb"
