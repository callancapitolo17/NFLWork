"""Config for the Kalshi MLB MM (maker) bot. Loaded from .env or environment."""
import os
from pathlib import Path

PKG_DIR = Path(__file__).parent
_RAW_ROOT = PKG_DIR.parent
PROJECT_ROOT = (Path(str(_RAW_ROOT).split(".worktrees")[0].rstrip("/"))
                if ".worktrees" in str(_RAW_ROOT) else _RAW_ROOT)
DB_PATH = PKG_DIR / "kalshi_mlb_mm.duckdb"
MARKET_DB = PKG_DIR / "kalshi_mlb_mm_market.duckdb"
KILL_FILE = PKG_DIR / ".kill"


def _load_env(path: Path) -> dict[str, str]:
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


def _get(key, default=None):
    return os.environ.get(key, _FILE_ENV.get(key, default))


# Credentials
KALSHI_API_KEY_ID = _get("KALSHI_API_KEY_ID")
KALSHI_PRIVATE_KEY_PATH = _get("KALSHI_PRIVATE_KEY_PATH")
KALSHI_USER_ID = _get("KALSHI_USER_ID")
KALSHI_BASE_URL = _get("KALSHI_BASE_URL", "https://api.elections.kalshi.com/trade-api/v2")
MVE_COLLECTION_TICKER = _get("MVE_COLLECTION_TICKER", "KXMVECROSSCATEGORY-R")

# Pricing
TARGET_ROI = float(_get("TARGET_ROI", "0.05"))
QUOTE_HYSTERESIS = float(_get("QUOTE_HYSTERESIS", "0.005"))

# Risk (master dial = BANKROLL)
BANKROLL = float(_get("BANKROLL", "500.0"))
DAILY_EXPOSURE_CAP_PCT = float(_get("DAILY_EXPOSURE_CAP_PCT", "0.75"))
MAX_GAME_EXPOSURE_PCT = float(_get("MAX_GAME_EXPOSURE_PCT", "0.10"))
MAX_RFQ_CONTRACTS = int(_get("MAX_RFQ_CONTRACTS", "5"))
MAX_OPEN_QUOTES = int(_get("MAX_OPEN_QUOTES", "25"))
FAIR_DRIFT_TOLERANCE = float(_get("FAIR_DRIFT_TOLERANCE", "0.02"))
MIN_FAIR_PROB = float(_get("MIN_FAIR_PROB", "0.05"))
MAX_FAIR_PROB = float(_get("MAX_FAIR_PROB", "0.95"))

# Freshness / circuit breaker
MAX_BOOK_STALENESS_SEC = int(_get("MAX_BOOK_STALENESS_SEC", "60"))
BOOK_MOVE_CB_THRESHOLD = float(_get("BOOK_MOVE_CB_THRESHOLD", "0.03"))
TIPOFF_CANCEL_MIN = int(_get("TIPOFF_CANCEL_MIN", "5"))

# Book consensus gate (v1 correlation defense — mirrors MLB answer-key dashboard
# pattern). v1.1: see docs/superpowers/specs/2026-05-26-kalshi-mlb-mm-design.md §13
# for the explicit correlation-premium gate (deferred enhancement).
BOOK_CONSENSUS_BAND = float(_get("BOOK_CONSENSUS_BAND", "0.02"))
MIN_AGREEING_BOOKS = int(_get("MIN_AGREEING_BOOKS", "3"))

# Loops (seconds)
DISCOVERY_SEC = int(_get("DISCOVERY_SEC", "2"))
CONFIRM_SEC = int(_get("CONFIRM_SEC", "2"))
RISK_SWEEP_SEC = int(_get("RISK_SWEEP_SEC", "10"))
RECONCILE_SWEEP_SEC = int(_get("RECONCILE_SWEEP_SEC", "30"))
SGP_REFRESH_SEC = int(_get("SGP_REFRESH_SEC", "60"))
SGP_SCRAPER_TIMEOUT_SEC = int(_get("SGP_SCRAPER_TIMEOUT_SEC", "90"))

# Adverse-selection halts (H4)
VOID_RATE_HALT_THRESHOLD = float(_get("VOID_RATE_HALT_THRESHOLD", "0.25"))
VOID_RATE_WINDOW_HOURS = int(_get("VOID_RATE_WINDOW_HOURS", "1"))
PER_CREATOR_FILL_HALT = int(_get("PER_CREATOR_FILL_HALT", "10"))
PER_CREATOR_WINDOW_HOURS = int(_get("PER_CREATOR_WINDOW_HOURS", "24"))

# Per-combo concentration controls (H8 / H9)
MAX_COMBO_EXPOSURE_USD = float(_get("MAX_COMBO_EXPOSURE_USD", "10.0"))
COMBO_COOLDOWN_SEC = int(_get("COMBO_COOLDOWN_SEC", "60"))

# Vig fallbacks (only used if a combo lacks full 4-side devig)
DK_VIG_FALLBACK = float(_get("DK_VIG_FALLBACK", "0.125"))
FD_VIG_FALLBACK = float(_get("FD_VIG_FALLBACK", "0.18"))
PX_VIG_FALLBACK = float(_get("PX_VIG_FALLBACK", "0.05"))
NOVIG_VIG_FALLBACK = float(_get("NOVIG_VIG_FALLBACK", "0.05"))

# Paths
MLB_SGP_DIR = Path(_get("MLB_SGP_DIR", str(PROJECT_ROOT / "mlb_sgp")))
NOTIFY_WEBHOOK_URL = _get("NOTIFY_WEBHOOK_URL")


def daily_exposure_cap_usd() -> float:
    return BANKROLL * DAILY_EXPOSURE_CAP_PCT
