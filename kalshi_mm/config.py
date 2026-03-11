"""
Configuration for Kalshi CBB 1H Market Maker.
All settings loaded from environment or defaults.
"""

import os
from pathlib import Path

# --- Paths ---
PROJECT_ROOT = Path(__file__).parent.parent
MM_DIR = Path(__file__).parent
KALSHI_DRAFT_DIR = PROJECT_ROOT / "kalshi_draft"
ANSWER_KEYS_DIR = PROJECT_ROOT / "Answer Keys"
CBB_DB_PATH = ANSWER_KEYS_DIR / "CBB Answer Key" / "cbb.duckdb"
MM_DB_PATH = MM_DIR / "kalshi_mm.duckdb"
BOOKMAKER_SCRAPER = PROJECT_ROOT / "bookmaker_odds" / "scraper.py"
BET105_SCRAPER = PROJECT_ROOT / "bet105_odds" / "scraper.py"

# --- Kalshi API ---
KALSHI_BASE_URL = os.getenv("KALSHI_BASE_URL",
                            "https://api.elections.kalshi.com/trade-api/v2")
# Set to "https://demo-api.kalshi.co/trade-api/v2" for demo environment

# --- Credentials (loaded from .env) ---
def load_env():
    """Load credentials from .env file."""
    env_path = MM_DIR / ".env"
    env = {}
    if env_path.exists():
        with open(env_path) as f:
            for line in f:
                line = line.strip()
                if "=" in line and not line.startswith("#"):
                    key, val = line.split("=", 1)
                    env[key.strip()] = val.strip()
    return env

_env = load_env()
KALSHI_API_KEY_ID = _env.get("KALSHI_API_KEY_ID", "")
KALSHI_PRIVATE_KEY_PATH = _env.get("KALSHI_PRIVATE_KEY_PATH", "")

# --- Quoting Parameters ---
HALF_SPREAD_CENTS = int(_env.get("HALF_SPREAD_CENTS", "5"))
CONTRACT_SIZE = int(_env.get("CONTRACT_SIZE", "2"))
SKEW_PER_CONTRACT = int(_env.get("SKEW_PER_CONTRACT", "1"))

# --- Risk Limits ---
MAX_POSITION_PER_MARKET = int(_env.get("MAX_POSITION_PER_MARKET", "5"))
MAX_POSITION_PER_EVENT = int(_env.get("MAX_POSITION_PER_EVENT", "8"))  # Across all strikes for one game
MAX_TOTAL_EXPOSURE_DOLLARS = float(_env.get("MAX_TOTAL_EXPOSURE", "50.0"))
MAX_MARKETS = int(_env.get("MAX_MARKETS", "10"))
MAX_STALENESS_SEC = int(_env.get("MAX_STALENESS_SEC", "600"))
MIN_QUOTE_SPREAD_CENTS = int(_env.get("MIN_QUOTE_SPREAD", "4"))
MIN_FAIR_VALUE = 10  # Don't quote if fair_yes < 10 cents
MAX_FAIR_VALUE = 90  # Don't quote if fair_yes > 90 cents
TIPOFF_PULLBACK_MIN = 30  # Pull quotes this many min before tipoff

# --- Line Move Monitoring ---
LINE_MOVE_THRESHOLD = float(_env.get("LINE_MOVE_THRESHOLD", "1.0"))

# --- Timing (seconds) ---
QUOTE_CYCLE_SEC = 10
FILL_POLL_SEC = 30
MONITOR_CYCLE_SEC = 60

# --- Kalshi Series ---
SPREAD_SERIES = "KXNCAAMB1HSPREAD"
