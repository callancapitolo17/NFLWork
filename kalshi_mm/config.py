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
CBB_DB_PATH = ANSWER_KEYS_DIR / "cbb.duckdb"
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
MIN_EV_PCT = float(_env.get("MIN_EV_PCT", "0.05"))  # 5% minimum EV to quote
CONTRACT_SIZE = int(_env.get("CONTRACT_SIZE", "5"))
SKEW_PER_CONTRACT = int(_env.get("SKEW_PER_CONTRACT", "0"))  # Disabled — Kelly handles inventory

# --- Operational Limits ---
MAX_STALENESS_SEC = int(_env.get("MAX_STALENESS_SEC", "600"))
MIN_QUOTE_SPREAD_CENTS = int(_env.get("MIN_QUOTE_SPREAD", "4"))
MIN_FAIR_VALUE = 10  # Don't quote if fair_yes < 10 cents
MAX_FAIR_VALUE = 90  # Don't quote if fair_yes > 90 cents
TIPOFF_PULLBACK_MIN = 5  # Cancel all orders & stop quoting/taking this many min before tipoff

# --- Line Move Monitoring ---
LINE_MOVE_THRESHOLD = float(_env.get("LINE_MOVE_THRESHOLD", "0.5"))

# --- Timing (seconds) ---
QUOTE_CYCLE_SEC = 10
FILL_POLL_SEC = 30
MONITOR_CYCLE_SEC = 60
PIPELINE_REFRESH_SEC = int(_env.get("PIPELINE_REFRESH_SEC", "600"))  # Auto-refresh predictions

# --- Taker Parameters ---
MIN_TAKE_EV_PCT = float(_env.get("MIN_TAKE_EV_PCT", "0.05"))  # 5% minimum EV after fees
TAKE_CONTRACT_SIZE = int(_env.get("TAKE_CONTRACT_SIZE", "5"))
TAKE_POLL_SEC = float(_env.get("TAKE_POLL_SEC", "2.0"))
TAKER_FEE_RATE = 0.07  # Kalshi taker fee: 7% * P * (1-P)
MAKER_FEE_RATE = float(_env.get("MAKER_FEE_RATE", "0.0"))  # Currently 0% for most accounts

# --- Kalshi Series ---
MARKET_SERIES = {
    "spreads": "KXNCAAMB1HSPREAD",
    "totals": "KXNCAAMB1HTOTAL",
    "moneyline": "KXNCAAMB1HWINNER",
    "race_to_10": "KXNCAAMBFIRST10",
}
ENABLED_MARKET_TYPES = set(_env.get("ENABLED_MARKETS", "spreads,totals,moneyline,race_to_10").split(","))
SPREAD_SERIES = MARKET_SERIES["spreads"]  # backward compat
MAX_BOOK_STALENESS_SEC = int(_env.get("MAX_BOOK_STALENESS_SEC", "10"))

# --- Kelly Sizing ---
BANKROLL = float(_env.get("BANKROLL", "1000.0"))
KELLY_FRACTION = float(_env.get("KELLY_FRACTION", "0.25"))  # Quarter Kelly
USE_KELLY_SIZING = _env.get("USE_KELLY_SIZING", "true").lower() == "true"

# --- Hard Position Limits (circuit breakers, independent of Kelly) ---
# Max filled exposure per game per market type as fraction of bankroll.
# e.g., 0.30 = 30% of $2K bankroll = $600 max on spreads for one game.
MAX_GAME_TYPE_EXPOSURE_PCT = float(_env.get("MAX_GAME_TYPE_EXPOSURE_PCT", "0.30"))
