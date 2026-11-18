"""Config for the dry-run parlay quoter. Environment only (the VM mounts the
repo read-only, so there is no .env inside the package to write)."""
import os
from pathlib import Path

PKG_DIR = Path(__file__).parent
PROJECT_ROOT = PKG_DIR.parent


def _get(key: str, default: str | None = None) -> str | None:
    return os.environ.get(key, default)


def _key_from_env_file(path: str | None) -> str | None:
    """KALSHI_API_KEY_ID only, from a KEY=value file (on the VM:
    bets_service/.env), so no other secret in that file enters this process."""
    if not path or not Path(path).exists():
        return None
    for raw in Path(path).read_text().splitlines():
        k, _, v = raw.strip().partition("=")
        if k.strip() == "KALSHI_API_KEY_ID":
            return v.strip().strip('"').strip("'") or None
    return None


# Kalshi read access: the websocket handshake is signed. A read-only key is
# enough; this process never places or quotes, and makes no Kalshi REST call
# (KALSHI_BASE_URL is only what auth_client.configure() requires).
KALSHI_API_KEY_ID = (_get("KALSHI_API_KEY_ID")
                     or _key_from_env_file(_get("PARLAY_KALSHI_ENV_FILE")))
KALSHI_PRIVATE_KEY_PATH = _get("KALSHI_PRIVATE_KEY_PATH")
KALSHI_BASE_URL = _get("KALSHI_BASE_URL",
                       "https://api.elections.kalshi.com/trade-api/v2")

# Decision log. On the VM this lives under /data (the only writable mount).
DB_PATH = Path(_get("PARLAY_QUOTER_DB_PATH",
                    str(PKG_DIR / "kalshi_parlay_quoter.duckdb")))

# The RFQ feed (feed.py). A socket with no frame for this long is dead:
# the feed is ~9,900 frames a second and Kalshi pings every ~10 s.
KALSHI_WS_URL = _get("KALSHI_WS_URL", "wss://api.elections.kalshi.com/trade-api/ws/v2")
FEED_IDLE_TIMEOUT_SEC = float(_get("PARLAY_FEED_IDLE_TIMEOUT_SEC", "30"))
# Kept RFQs waiting for the decision loop; beyond it they are counted
# queue_dropped (about 40 s of the ~475 a second measured 2026-11-13).
FEED_QUEUE_MAX = int(_get("PARLAY_FEED_QUEUE_MAX", "20000"))

# Door: RFQs below this dollar size (or, when sized in contracts, with fewer
# contracts than this many dollars) are counted per minute, never decided
# (the firehose is mostly ~$10 algorithmic spam). 0 keeps every size.
MIN_RFQ_TARGET_COST_USD = float(_get("PARLAY_MIN_RFQ_TARGET_COST_USD", "50"))

# The runner's leg pricers (#282), loopback only.
RUNNER_URL = _get("UNABATED_RUNNER_URL", "http://127.0.0.1:8095").rstrip("/")

# Leg table (legtable.py): every REFRESH_SEC the runner's league snapshot
# times are read off /health, and the legs seen in the last LEG_TTL whose
# league changed are re-priced by board-only runner calls.
TABLE_REFRESH_SEC = float(_get("PARLAY_TABLE_REFRESH_SEC", "2"))
TABLE_HEALTH_TIMEOUT_SEC = float(_get("PARLAY_TABLE_HEALTH_TIMEOUT_SEC", "3"))
# A leg is re-priced when its league's snapshot changes, and at the latest
# this long after it was priced.
TABLE_MAX_REUSE_SEC = float(_get("PARLAY_TABLE_MAX_REUSE_SEC", "90"))
TABLE_LEG_TTL_SEC = float(_get("PARLAY_TABLE_LEG_TTL_SEC", "600"))
# Fairs are served while the snapshot times were read under this long ago,
# and for this long after a league's snapshot changed while its legs are
# re-priced. 15 s, not 2x the refresh: the runner's board thread is busy for
# 5-14 s at a time building Edges, and /health and a re-price wait behind
# that. Every row logs the fair's true age.
TABLE_MAX_AGE_SEC = float(_get("PARLAY_TABLE_MAX_AGE_SEC", "15"))
TABLE_FETCH_TIMEOUT_SEC = float(_get("PARLAY_TABLE_FETCH_TIMEOUT_SEC", "15"))
# An RFQ with a leg the table has not priced yet waits this long for it
# (new legs are priced within ~0.25 s when the runner answers; RFQs live
# ~10 s), then declines pricer_error.
PENDING_MAX_WAIT_SEC = float(_get("PARLAY_PENDING_MAX_WAIT_SEC", "5"))
# An RFQ this old when the decision loop reaches it is counted
# expired_in_queue and not decided: its requester is gone.
MAX_RFQ_AGE_SEC = float(_get("PARLAY_MAX_RFQ_AGE_SEC", "15"))
# RFQs waiting that way beyond this are logged pricer_backlog.
MAX_BACKLOG = int(_get("PARLAY_MAX_BACKLOG", "2000"))

# Confirm: a parlay the table would quote is re-priced through
# /rfq_price.json, which reads Kalshi's market, book and start per leg. Those
# reads share Kalshi's 20 reads/s limit with the taker bot, so only this many
# would-quote parlays a second are confirmed; the rest are counted
# would_quote_unconfirmed. Off (0) by default since 2026-11-13: on the VM
# most confirms timed out (16 confirmed, 288 timed out in 80 minutes) and
# their reads run on the taker's thread. To come back reading the runner's
# Kalshi book stream instead.
CONFIRMS_PER_SEC = float(_get("PARLAY_CONFIRMS_PER_SEC", "0"))
CONFIRM_WORKERS = int(_get("PARLAY_CONFIRM_WORKERS", "2"))
CONFIRM_TIMEOUT_SEC = float(_get("PARLAY_CONFIRM_TIMEOUT_SEC", "3"))
# After a confirm the runner failed or timed out on, no confirm is sent for
# this long: the runner keeps working a call the quoter gave up on, and its
# Kalshi reads queue behind 4 slots, so pressing on would pile reads up.
CONFIRM_BACKOFF_SEC = float(_get("PARLAY_CONFIRM_BACKOFF_SEC", "30"))

# Decision log: confirmed parlays are stored in full; of every other RFQ
# this percent is stored (chosen by a hash of its id) and all are counted
# per minute. 0.5% of ~475 RFQs a second plus 1 confirm a second is ~3.4
# stored decisions a second, measured at ~190 MB a day (2026-11-13).
STORED_SAMPLE_PCT = float(_get("PARLAY_STORED_SAMPLE_PCT", "0.5"))

FLUSH_SEC = float(_get("PARLAY_FLUSH_SEC", "1.0"))
STOP_PATH = Path(_get("PARLAY_STOP_PATH", str(PKG_DIR / ".stop")))
