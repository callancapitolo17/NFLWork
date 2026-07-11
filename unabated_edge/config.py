import os
from pathlib import Path

PKG_DIR = Path(__file__).resolve().parent


def _load_env(path: Path) -> dict[str, str]:
    env = {}
    if path.exists():
        for line in path.read_text().splitlines():
            line = line.strip()
            if line and not line.startswith("#") and "=" in line:
                k, v = line.split("=", 1)
                env[k.strip()] = v.strip().strip('"').strip("'")
    return env


_FILE_ENV = _load_env(PKG_DIR / ".env")


def _get(key, default=None):
    return os.environ.get(key, _FILE_ENV.get(key, default))


UNABATED_SNAPSHOT_URL = "https://content.unabated.com/markets/game-odds/b_gameodds.json"
UNABATED_CHANGES_URL = "https://api-k.unabated.com/api/markets/changes/query"
# v2 per-league odds file (discovered 2026-07-08 via app-bundle recon). Anchors are
# UNBLURRED here even anonymously; includes alternateLines ladders. The legacy
# changes/query delta feed does NOT carry soccer (lg21) at all — v2 polling replaces it.
UNABATED_V2_LEAGUE_URL = "https://content.unabated.com/markets/v2/league/{league_id}/odds.json"
V2_POLL_SEC = float(_get("UNABATED_V2_POLL_SEC", "5"))
# Kalshi executed-trades tape poll cadence (rides every Nth main tick)
TRADES_POLL_SEC = float(_get("KALSHI_TRADES_POLL_SEC", "30"))
UNABATED_TOKEN = _get("UNABATED_AT_PROD")
SHARP_BOOK_PRICE_ID = 7
ANCHOR_SOURCE_IDS = [7, 6, 68]

BANKROLL = float(_get("BANKROLL", "1000.0"))
KELLY_FRACTION = float(_get("KELLY_FRACTION", "0.25"))
MIN_EV_PCT = float(_get("MIN_EV_PCT", "0.03"))
# Absolute per-contract EV floor (dollars), gates alongside MIN_EV_PCT so a tiny
# absolute edge on a cheap longshot can't clear on percentage alone (ev_pct=ev/ask
# inflates as ask->0, surfacing the noisiest devig estimates).
MIN_EV_DOLLARS = float(_get("MIN_EV_DOLLARS", "0.02"))
MAX_STALENESS_SEC = int(_get("MAX_STALENESS_SEC", "20"))  # RESERVED — Plan 2 live-execution staleness gate (not yet enforced)
KICKOFF_CUTOFF_MIN = int(_get("KICKOFF_CUTOFF_MIN", "3"))
PER_MATCH_CAP_PCT = float(_get("PER_MATCH_CAP_PCT", "0.03"))

KALSHI_API_KEY_ID = _get("KALSHI_API_KEY_ID")
KALSHI_PRIVATE_KEY_PATH = _get("KALSHI_PRIVATE_KEY_PATH")
KALSHI_BASE_URL = _get("KALSHI_BASE_URL", "https://api.elections.kalshi.com/trade-api/v2")

MARKET_DB_PATH = PKG_DIR / "unabated_edge_market.duckdb"
RESEARCH_DB_PATH = PKG_DIR / "unabated_edge_research.duckdb"
LOG_PATH = PKG_DIR / "bot.log"
KILL_FILE = PKG_DIR / ".kill"
