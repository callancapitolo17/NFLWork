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


# Kalshi read access: the WS handshake and the REST gap-fill are signed. A
# read-only key is enough; this process never places or quotes.
KALSHI_API_KEY_ID = (_get("KALSHI_API_KEY_ID")
                     or _key_from_env_file(_get("PARLAY_KALSHI_ENV_FILE")))
KALSHI_PRIVATE_KEY_PATH = _get("KALSHI_PRIVATE_KEY_PATH")
KALSHI_BASE_URL = _get("KALSHI_BASE_URL",
                       "https://api.elections.kalshi.com/trade-api/v2")

# Decision log. On the VM this lives under /data (the only writable mount).
DB_PATH = Path(_get("PARLAY_QUOTER_DB_PATH",
                    str(PKG_DIR / "kalshi_parlay_quoter.duckdb")))

# The runner's leg pricer (#282), loopback only.
RUNNER_URL = _get("UNABATED_RUNNER_URL", "http://127.0.0.1:8095").rstrip("/")
PRICER_TIMEOUT_SEC = float(_get("PARLAY_PRICER_TIMEOUT_SEC", "3"))
PRICER_WORKERS = int(_get("PARLAY_PRICER_WORKERS", "8"))
# Each (ticker, side) is priced at most once per this many seconds; RFQs
# re-ask the same legs constantly (~100 in-league RFQs/s, 2026-11-07). 10 s,
# not 3: the runner reads Kalshi per uncached leg (4 in flight) and could
# not keep ~2k live legs fresh every 3 s (2026-11-12 VM run).
PRICE_CACHE_SEC = float(_get("PARLAY_PRICE_CACHE_SEC", "10"))
# After a runner failure, skip calls for this long (logged pricer_error).
PRICER_BACKOFF_SEC = float(_get("PARLAY_PRICER_BACKOFF_SEC", "5"))
# In-scope RFQs waiting for a pricer worker beyond this are logged as
# pricer_backlog instead of queueing without bound.
MAX_BACKLOG = int(_get("PARLAY_MAX_BACKLOG", "200"))

# Door: RFQs below this dollar size (or, when sized in contracts, with fewer
# contracts than this many dollars) are counted per minute, not logged as
# rows (the firehose is mostly ~$10 algorithmic spam). 0 logs every size.
MIN_RFQ_TARGET_COST_USD = float(_get("PARLAY_MIN_RFQ_TARGET_COST_USD", "50"))

# How long an RFQ id stays in the "already decided" set. Kalshi RFQs live
# ~10s p50; a decision is logged once per RFQ id.
SEEN_TTL_SEC = float(_get("PARLAY_SEEN_TTL_SEC", "900"))
MIRROR_TTL_SEC = float(_get("PARLAY_MIRROR_TTL_SEC", "60"))

# Legless WS frames need a GET /markets/{combo}; cap those per second.
MAX_MARKET_FETCHES_PER_SEC = float(_get("PARLAY_MAX_MARKET_FETCHES_PER_SEC", "5"))

TICK_SEC = float(_get("PARLAY_TICK_SEC", "0.25"))
FLUSH_SEC = float(_get("PARLAY_FLUSH_SEC", "1.0"))
STOP_PATH = Path(_get("PARLAY_STOP_PATH", str(PKG_DIR / ".stop")))
