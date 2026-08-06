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


def _get_bool(key: str, default: str) -> bool:
    """Env/file flag as a bool. Anything other than a recognised true-ish
    word is False, so a typo disables a gate rather than silently enabling
    one (fail toward the pre-ticket behaviour, never toward surprise)."""
    return str(_get(key, default)).strip().lower() in ("1", "true", "yes", "on")


# Credentials
KALSHI_API_KEY_ID = _get("KALSHI_API_KEY_ID")
KALSHI_PRIVATE_KEY_PATH = _get("KALSHI_PRIVATE_KEY_PATH")
KALSHI_USER_ID = _get("KALSHI_USER_ID")
KALSHI_BASE_URL = _get("KALSHI_BASE_URL", "https://api.elections.kalshi.com/trade-api/v2")
MVE_COLLECTION_TICKER = _get("MVE_COLLECTION_TICKER", "KXMVECROSSCATEGORY-R")

# Pricing
# Tightened 5%→3% (2026-06-18) to test competitiveness: 13 quotes floated over
# ~8 days got 0 accepts. Competitor quotes are invisible on Kalshi (403 to
# non-creators), so the only way to learn whether we're being outbid is to
# quote tighter and watch for fills. 3% per side = ~6% gross spread.
TARGET_ROI = float(_get("TARGET_ROI", "0.03"))
QUOTE_HYSTERESIS = float(_get("QUOTE_HYSTERESIS", "0.005"))
# Uncertainty-scaled margin (issue #19): per side,
#   margin_pts = max(p·(1−(1+TARGET_ROI)^−n_games), MIN_MARGIN_PTS + K_SIGMA·σ)
# where σ = sample stddev of the consensus books' fairs for the combo (prob
# points). The constant-ROI cushion ≈ 0.029·p collapses at longshot fairs
# (~0.3¢ at p=0.10) while absolute fair error does NOT shrink with p — every
# book loads margin onto longshots (favorite-longshot vig distribution) and
# options MMs widen deep-OTM quotes for the same reason. σ prices visible
# disagreement; MIN_MARGIN_PTS covers error σ cannot see (shared devig/
# correlation bias, mirrored books — a 2-book set can agree by coincidence).
# First-principles defaults (#13 markouts descoped): calibrate later from
# settlement P&L (#12) and the daily report's demand curve (#14).
MIN_MARGIN_PTS = float(_get("MIN_MARGIN_PTS", "0.01"))
K_SIGMA = float(_get("K_SIGMA", "1.0"))

# Risk (master dial = BANKROLL)
BANKROLL = float(_get("BANKROLL", "500.0"))
DAILY_EXPOSURE_CAP_PCT = float(_get("DAILY_EXPOSURE_CAP_PCT", "0.75"))
MAX_GAME_EXPOSURE_PCT = float(_get("MAX_GAME_EXPOSURE_PCT", "0.10"))
# Per-fill exposure cap = max DOLLARS at risk on a single fill (our cost basis
# on the side we'd hold = our max loss, since a binary settles to 0/1). The
# maker can't choose fill size — only quote-or-skip — so this is the one
# per-fill size lever. Replaces the old contract-count proxy MAX_RFQ_CONTRACTS.
MAX_FILL_EXPOSURE_PCT = float(_get("MAX_FILL_EXPOSURE_PCT", "0.10"))
MAX_OPEN_QUOTES = int(_get("MAX_OPEN_QUOTES", "25"))
FAIR_DRIFT_TOLERANCE = float(_get("FAIR_DRIFT_TOLERANCE", "0.02"))
MIN_FAIR_PROB = float(_get("MIN_FAIR_PROB", "0.05"))
MAX_FAIR_PROB = float(_get("MAX_FAIR_PROB", "0.95"))

# Freshness / circuit breaker
# Observed live 2026-06-08: a full SGP cycle is ~150-165s end-to-end (60s tick
# gap + ~90-105s scrape runtime), so a 60s staleness window left the bot
# books-stale (idle) for most of every cycle. 180s covers the worst observed
# inter-fetch gap; pickoff defense in the stale tail remains margin + circuit
# breaker + last-look.
MAX_BOOK_STALENESS_SEC = int(_get("MAX_BOOK_STALENESS_SEC", "180"))
BOOK_MOVE_CB_THRESHOLD = float(_get("BOOK_MOVE_CB_THRESHOLD", "0.03"))
TIPOFF_CANCEL_MIN = int(_get("TIPOFF_CANCEL_MIN", "5"))

# Book consensus gate (issue #20, rescoped 2026-07-25): z-space dispersion
# threshold, replacing the old absolute ±2¢ outlier band (BOOK_CONSENSUS_BAND,
# removed). Books' devigged combo fairs go through norm.ppf; we quote only if
# the sample stddev (ddof=1) of those z-values is <= SIGMA_Z_MAX. Constant
# width in z-space = the same amount of DISAGREEMENT at every price level —
# the tolerated absolute gap tightens automatically at the tails (~2¢ at
# p=0.50 → ~0.6¢ at p=0.08), where the ±2¢ band tolerated 25% relative
# disagreement. No outlier removal: a dissenting book is as likely the
# informed one (news mid-propagation) as a broken scrape, so large dispersion
# DECLINES the quote instead of outvoting the dissenter (books suspend on
# news; a false decline is ~free, a false quote is not). Default 0.07 keeps
# continuity with the old gate at p=0.50 (2-book set 4¢ apart: z-gap ≈ 0.100
# → sample stddev ≈ 0.071); calibrate later from #12/#14 data. v1.1: see
# docs/superpowers/specs/2026-05-26-kalshi-mlb-mm-design.md §13 for the
# explicit correlation-premium gate (deferred enhancement).
SIGMA_Z_MAX = float(_get("SIGMA_Z_MAX", "0.07"))
# Lowered 3→2 (2026-06-18, user-approved): at 3 books we quoted ~13 times in
# 8 days (0 fills) — too few to test competitiveness or gather data. Measured
# ~5 quotable tuples @3 vs ~40-47 @2 books (~8×). Trade-off: 2-book consensus
# is weaker (more model risk); caveat: a DK+Novig pair is effectively ONE
# independent source since Novig mirrors DK (see [[novig_sgp_scraping]]) — a
# DK/Novig-independence guard is a noted follow-up if those pairs dominate fills.
MIN_AGREEING_BOOKS = int(_get("MIN_AGREEING_BOOKS", "2"))

# Correlation sanity vs Kalshi's own single-leg markets (issue #23, spec §13).
# Every leg of a combo IS its own 2-way Kalshi market trading in REAL TIME,
# which makes its devigged marginal the one pricing input we have that is both
# independent of the books and fast. Two tests on the combo fair:
#
#   Frechet:  max(0, Σp − (n−1)) <= combo_fair <= min(p)      [always true]
#   premium:  combo_fair / Πp  ∈ [CORR_PREMIUM_MIN, MAX]      [heuristic]
#
# Frechet gates by default: it is parameter-free, and a violation means our
# fair is arithmetically not a joint probability of these legs. The premium
# band ships LOG-ONLY (enabled=False) because same-game legs legitimately run
# well above 1× (run line + moneyline ~2×) and a guessed band would decline
# real business; switch it on once the `corr_sanity_check` firehose shows the
# true premium distribution. This is a LEVEL check against an outside anchor
# and is deliberately independent of SIGMA_Z_MAX, which is a DISPERSION check
# among the books themselves — tightly-agreeing books can still be jointly
# wrong, and only this gate can see that.
CORR_SANITY_FRECHET_ENABLED = _get_bool("CORR_SANITY_FRECHET_ENABLED", "true")
CORR_SANITY_PREMIUM_ENABLED = _get_bool("CORR_SANITY_PREMIUM_ENABLED", "false")
CORR_PREMIUM_MIN = float(_get("CORR_PREMIUM_MIN", "0.5"))
CORR_PREMIUM_MAX = float(_get("CORR_PREMIUM_MAX", "2.0"))

# Constituent-jump circuit breaker (issue #23 items 1-2). Our books refresh
# every ~150-165s; the combo's constituent Kalshi singles trade in real time.
# If a constituent's devigged mid moves past this threshold AFTER we quoted,
# the market has moved and our resting quote is the stale side.
CONSTITUENT_JUMP_THRESHOLD = float(_get("CONSTITUENT_JUMP_THRESHOLD", "0.03"))
# Mode flag, flipped by the live-pricing pivot (#54) as designed:
#   "unconditional" (default, post-#54) — every quote is priced from a live
#                   fetch, so the fair was fresh at placement and ANY
#                   constituent jump during the resting window means the
#                   market moved after us.
#   "book_quiet"    (pre-#54 behavior) — also required our book consensus to
#                   have stayed within CONSTITUENT_BOOK_QUIET_MAX, separating
#                   "we are the stale ones" from "the whole market moved
#                   together". Only made sense while quotes were priced off a
#                   <=150s cache; kept as an env override.
CONSTITUENT_JUMP_MODE = _get("CONSTITUENT_JUMP_MODE", "unconditional")
CONSTITUENT_BOOK_QUIET_MAX = float(_get("CONSTITUENT_BOOK_QUIET_MAX", "0.01"))
# Wall-clock ceiling on the constituent poll. ALL ticks share one thread, so a
# slow poll inside the risk sweep delays the confirm tick — and Kalshi allows
# only 2s to confirm in High Volatility Markets. Blowing a confirm window is a
# far worse failure than polling fewer tickers this pass, so the poll yields
# first and the iteration order rotates (main._rotated_poll_order) so no ticker
# is starved across sweeps. At ~100ms/GET this covers ~10 tickers per 10s
# sweep; a large resting book is therefore swept over several passes. The real
# fix for wide coverage is the WebSocket feed the transport interface allows.
CONSTITUENT_POLL_BUDGET_SEC = float(_get("CONSTITUENT_POLL_BUDGET_SEC", "1.0"))

# Loops (seconds)
DISCOVERY_SEC = int(_get("DISCOVERY_SEC", "2"))
CONFIRM_SEC = int(_get("CONFIRM_SEC", "2"))
RISK_SWEEP_SEC = int(_get("RISK_SWEEP_SEC", "10"))
RECONCILE_SWEEP_SEC = int(_get("RECONCILE_SWEEP_SEC", "30"))
SGP_REFRESH_SEC = int(_get("SGP_REFRESH_SEC", "60"))
# Settlement sweep (issue #12): populate fills.realized_pnl once markets
# settle. Only matters hours post-game, so a slow cadence is plenty.
SETTLEMENT_SWEEP_SEC = int(_get("SETTLEMENT_SWEEP_SEC", "600"))
SGP_SCRAPER_TIMEOUT_SEC = int(_get("SGP_SCRAPER_TIMEOUT_SEC", "90"))
# Issue #50: structure-only warming cadence. Keeps every book's
# events/structure TTL caches (and Caesars' 240s WAF token) warm so an RFQ
# never pays a cold-structure penalty. Must stay under STRUCTURE_TTL_SEC
# (kalshi_common/sgp_service.py) and under the CZR token TTL.
STRUCTURE_WARM_SEC = int(_get("STRUCTURE_WARM_SEC", "120"))
# Issue #50: per-book wall budget for LIVE (on-demand) pricing fetches.
# 10s keeps warm Novig (p95 ~9s at #42's baseline) barely inside while a
# hung book is dropped instead of stalling the combo to the sweep budget.
ON_DEMAND_DEADLINE_SEC = float(_get("ON_DEMAND_DEADLINE_SEC", "10.0"))

# Run-time book-health alerting (issue #37). Keys on consecutive FAILED
# fetches, never on data age — an age rule would false-fire by design once
# #57 slows the sweep to background structure-warming.
BOOK_ALERT_ENABLED = _get_bool("BOOK_ALERT_ENABLED", "true")
# One 403 is noise (a book hiccups); three consecutive is a dead book.
BOOK_ALERT_STREAK = int(_get("BOOK_ALERT_STREAK", "3"))
# Which fetch paths count toward a book's health. Post-#53 every quote is
# priced by an on-demand fetch and the sweep is legacy warming, so #57 flips
# this to "on_demand" in CONFIG — no code change.
BOOK_ALERT_PATHS = tuple(
    p.strip() for p in _get("BOOK_ALERT_PATHS", "sweep,on_demand").split(",")
    if p.strip())
# Rule B's floor deliberately REUSES MIN_AGREEING_BOOKS rather than adding a
# second knob: an alert that fires at a different count than the gate it is
# warning about is worse than no alert.

# Adverse-selection halts (H4)
VOID_RATE_HALT_THRESHOLD = float(_get("VOID_RATE_HALT_THRESHOLD", "0.25"))
VOID_RATE_WINDOW_HOURS = int(_get("VOID_RATE_WINDOW_HOURS", "1"))
# NOTE (verified live 2026-06-08): Kalshi anonymizes creator_id to "" in the
# market-wide RFQ poll, so the per-creator halt is currently INERT (empty id
# short-circuits to no-halt). Kept wired in case Kalshi populates the field;
# the research firehose captures rfq_raw so we'll see it if that changes.
PER_CREATOR_FILL_HALT = int(_get("PER_CREATOR_FILL_HALT", "10"))
PER_CREATOR_WINDOW_HOURS = int(_get("PER_CREATOR_WINDOW_HOURS", "24"))

# Per-combo concentration controls (H8 / H9)
# Per-combo exposure cap. Must be >= max_fill_exposure_usd() (else a single max
# fill self-blocks). At BANKROLL=$500: per-fill $50 = per-combo $50 = per-game
# $50; diversification is across distinct games under the $375 daily cap.
MAX_COMBO_EXPOSURE_USD = float(_get("MAX_COMBO_EXPOSURE_USD", "50.0"))
COMBO_COOLDOWN_SEC = int(_get("COMBO_COOLDOWN_SEC", "60"))

# Reconcile max-age fallback (N11): fills older than this with positions API
# persistently down get marked reconciled=TRUE with recorded values.
MAX_RECONCILE_AGE_SEC = int(_get("MAX_RECONCILE_AGE_SEC", "300"))

# Vig fallbacks (only used if a combo lacks full 4-side devig)
DK_VIG_FALLBACK = float(_get("DK_VIG_FALLBACK", "0.125"))
FD_VIG_FALLBACK = float(_get("FD_VIG_FALLBACK", "0.18"))
PX_VIG_FALLBACK = float(_get("PX_VIG_FALLBACK", "0.05"))
NOVIG_VIG_FALLBACK = float(_get("NOVIG_VIG_FALLBACK", "0.05"))

# Paths
MLB_SGP_DIR = Path(_get("MLB_SGP_DIR", str(PROJECT_ROOT / "mlb_sgp")))
NOTIFY_WEBHOOK_URL = _get("NOTIFY_WEBHOOK_URL")

# Logging
LOG_PATH = PKG_DIR / "bot.log"
LOG_LEVEL = _get("LOG_LEVEL", "INFO")
LOG_ROTATE_MAX_BYTES = int(_get("LOG_ROTATE_MAX_BYTES", str(50 * 1024 * 1024)))
LOG_ROTATE_BACKUPS = int(_get("LOG_ROTATE_BACKUPS", "5"))

# Research firehose
RESEARCH_DB_PATH = PKG_DIR / "kalshi_mlb_mm_research.duckdb"
RESEARCH_BUFFER_MAX = int(_get("RESEARCH_BUFFER_MAX", "5000"))
RESEARCH_FLUSH_WARN_RATE_LIMIT_SEC = int(_get("RESEARCH_FLUSH_WARN_RATE_LIMIT_SEC", "60"))


def daily_exposure_cap_usd() -> float:
    return BANKROLL * DAILY_EXPOSURE_CAP_PCT


def max_fill_exposure_usd() -> float:
    return BANKROLL * MAX_FILL_EXPOSURE_PCT
