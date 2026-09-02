"""Config for the one-sided Kalshi YRFI maker. Loaded from .env or environment.

Reads: kalshi_rfi/.env (gitignored) then process env (env wins).
Side effects: none — pure constants; auth_client.configure() happens in main.
"""
import os
from pathlib import Path

PKG_DIR = Path(__file__).parent
_RAW_ROOT = PKG_DIR.parent
# Worktree-aware root (kalshi_mlb_mm pattern): creds + kalshi_draft/auth.py
# live in the main checkout, not the worktree.
PROJECT_ROOT = (Path(str(_RAW_ROOT).split(".worktrees")[0].rstrip("/"))
                if ".worktrees" in str(_RAW_ROOT) else _RAW_ROOT)
DB_PATH = PKG_DIR / "kalshi_rfi.duckdb"
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


# Credentials (same names as the other bots so one .env works everywhere)
KALSHI_API_KEY_ID = _get("KALSHI_API_KEY_ID")
KALSHI_PRIVATE_KEY_PATH = _get("KALSHI_PRIVATE_KEY_PATH")
KALSHI_BASE_URL = _get("KALSHI_BASE_URL",
                       "https://api.elections.kalshi.com/trade-api/v2")

# Mode: off | shadow | live. Live additionally requires RFI_LIVE_ACK=1
# (dead-man switch, same contract as the unabated_edge maker).
RFI_MODE = _get("RFI_MODE", "shadow")
RFI_LIVE_ACK = _get("RFI_LIVE_ACK")

# Strategy: we ONLY buy YES (YRFI) as a resting maker bid at fair − margin.
# Thesis (user decision 2026-08-25): retail flow buys NRFI (NO), which
# matches against resting YES bids — we are the other side at a discount.
MARGIN_CENTS = int(_get("RFI_MARGIN_CENTS", "3"))
# Post-fee floor: quote only if fair − price − maker_fee ≥ this many cents.
MIN_EDGE_CENTS = float(_get("RFI_MIN_EDGE_CENTS", "2"))

# Risk caps, in worst-case dollars (a YES bid's worst case = its cost).
PER_GAME_CAP_USD = float(_get("RFI_PER_GAME_CAP_USD", "10"))
DAILY_CAP_USD = float(_get("RFI_DAILY_CAP_USD", "100"))
DAILY_ROLL_HOUR_ET = float(_get("RFI_DAILY_ROLL_HOUR_ET", "6"))

# Cadence
CYCLE_SEC = float(_get("RFI_CYCLE_SEC", "30"))
# Book consensus TTL: refresh faster near first pitch (lineup-news window).
# Near-window 60s (was 120s): fetches are genuinely live since the
# structure_ttl_sec=0 fix, so the extra cadence buys real staleness cuts
# in the exact window where news pick-offs happen.
FAIR_REFRESH_FAR_SEC = float(_get("RFI_FAIR_REFRESH_FAR_SEC", "300"))
FAIR_REFRESH_NEAR_SEC = float(_get("RFI_FAIR_REFRESH_NEAR_SEC", "60"))
NEAR_WINDOW_MIN = float(_get("RFI_NEAR_WINDOW_MIN", "60"))
SETTLEMENT_POLL_SEC = float(_get("RFI_SETTLEMENT_POLL_SEC", "600"))
# Stand down on a game for this long after each fill. At a 1c margin the
# Kalshi jump guard re-quotes every ~30s, which took one game to 5 fills /
# 45 contracts in 60 seconds against a $5 cap (bug 1, observed 2026-08-27).
POST_FILL_COOLDOWN_SEC = float(_get("RFI_POST_FILL_COOLDOWN_SEC", "60"))
# One live 4-book fair fetch's wall budget; books still running past it
# count as not-priced this refresh (finding 8: one hung book must never
# freeze quote management for every other game).
FAIR_FETCH_WALL_SEC = float(_get("RFI_FAIR_FETCH_WALL_SEC", "15"))

# Quote window: only games starting within this horizon; pull quotes this
# many seconds before first pitch. 600s (was 120s): late scratches land
# T-10min to T-2min and books' derivative markets lag them — the sharpest
# pick-off window is not worth the marginal fill volume.
QUOTE_HORIZON_HOURS = float(_get("RFI_QUOTE_HORIZON_HOURS", "12"))
PULL_BEFORE_START_SEC = float(_get("RFI_PULL_BEFORE_START_SEC", "600"))

# Between fair refreshes, Kalshi's own market is the cheap movement guard:
# a mid move ≥ this many cents since our fair was computed cancels the quote
# and forces a refresh (constituent-jump idea, single-market version).
KALSHI_JUMP_CENTS = float(_get("RFI_KALSHI_JUMP_CENTS", "3"))

# Book consensus gate (same semantics as kalshi_mlb_mm issue #20).
MIN_BOOKS = int(_get("RFI_MIN_BOOKS", "2"))
SIGMA_Z_MAX = float(_get("RFI_SIGMA_Z_MAX", "0.07"))
# The 5 books with working period-aware I1 hooks (issue #87; Caesars mapped
# 2026-08-25 off its "Any Run In 1st Inning?" Yes/No market). ProphetX is
# name-mapped but has never carried a priced I1 line, and 403s at the event
# stage as of 2026-08-25 — it declines cleanly but costs wire calls, so it
# stays deliberately absent.
_DEFAULT_BOOKS = "draftkings,fanduel,betmgm,novig,caesars"
BOOKS = tuple(b.strip() for b in _get("RFI_BOOKS", _DEFAULT_BOOKS).split(",")
              if b.strip())

# Logging
LOG_PATH = Path(_get("RFI_LOG_PATH", str(PKG_DIR / "bot.log")))
LOG_LEVEL = _get("RFI_LOG_LEVEL", "INFO")
LOG_ROTATE_MAX_BYTES = int(_get("RFI_LOG_ROTATE_MAX_BYTES", str(10 * 1024 * 1024)))
LOG_ROTATE_BACKUPS = int(_get("RFI_LOG_ROTATE_BACKUPS", "3"))
