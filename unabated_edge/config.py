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
# Measured full-slate MLB tick (176+ KXMLBTOTAL markets, one book fetch per
# market) runs ~37s end to end. This threshold must exceed the worst-case
# tick duration or the maker watchdog self-triggers every iteration once
# quotes rest (each adapter is judged on its own last-successful-tick clock,
# not a shared iteration-start one — see runner.py's per-adapter now()).
MAX_STALENESS_SEC = int(_get("MAX_STALENESS_SEC", "120"))
KICKOFF_CUTOFF_MIN = int(_get("KICKOFF_CUTOFF_MIN", "3"))
PER_MATCH_CAP_PCT = float(_get("PER_MATCH_CAP_PCT", "0.03"))
# Beyond this many hours out from first pitch, skip the per-tick Kalshi book
# + trades fetch for an event (line_snapshots still capture it). Shortens a
# full-slate tick and caps capture volume for games nobody can trade yet.
BOOK_CAPTURE_HORIZON_HOURS = float(_get("BOOK_CAPTURE_HORIZON_HOURS", "12"))
# Go-forward retention for the high-volume capture tables (book_snapshots,
# kalshi_trades), scoped to CAPTURE_PRUNE_SPORTS below — it prunes old rows
# for those sports on every restart/day, including rows that predate this
# feature. It does NOT touch line_snapshots, and it does not touch any sport
# absent from CAPTURE_PRUNE_SPORTS.
CAPTURE_RETENTION_DAYS = int(_get("CAPTURE_RETENTION_DAYS", "14"))
# Sports the automated prune is allowed to delete from, comma-separated.
# Soccer is deliberately excluded: the WC-era book_snapshots/kalshi_trades
# rows (2026-07-10 -> 19) in the shared market DB are a preserved backtest
# archive (GitHub issue #9), not disposable capture — an unscoped prune
# would delete a large chunk of it on first launch from main. Add a sport
# here only when its accumulated history is genuinely disposable.
CAPTURE_PRUNE_SPORTS = tuple(
    s.strip() for s in _get("CAPTURE_PRUNE_SPORTS", "mlb").split(",") if s.strip()
)

# ---- maker (unabated_edge/maker/) — spec docs/superpowers/specs/2026-07-10-wc-totals-maker-design.md ----
MAKER_MODE = _get("MAKER_MODE", "off")            # off | shadow | live
MAKER_LIVE_ACK = _get("MAKER_LIVE_ACK")           # dead-man switch: must be "1" for live
ROI_MARGIN = float(_get("ROI_MARGIN", "0.03"))
PICKOFF_BUFFER_CENTS = int(_get("PICKOFF_BUFFER_CENTS", "1"))
MAX_MARGIN_CENTS = int(_get("MAX_MARGIN_CENTS", "5"))
ALT_MARGIN_MULT = float(_get("ALT_MARGIN_MULT", "1.5"))
ALT_SIZE_MULT = float(_get("ALT_SIZE_MULT", "0.5"))
ALT_OVERROUND_MIN = float(_get("ALT_OVERROUND_MIN", "1.01"))
ALT_OVERROUND_MAX = float(_get("ALT_OVERROUND_MAX", "1.15"))
# Feed-integrity envelope on EVERY anchor rung (main and alt) at ladder build
# (issue #73): raw implied sum < 1 -> anchor_crossed; outside [MIN, MAX] ->
# anchor_overround; rejected rungs are never devigged. The envelope strictly
# contains the alt band above, so the maker's tighter alt_overround gate keeps
# its original behavior for everything that passes here.
ANCHOR_OVERROUND_MIN = float(_get("ANCHOR_OVERROUND_MIN", "1.005"))
ANCHOR_OVERROUND_MAX = float(_get("ANCHOR_OVERROUND_MAX", "1.20"))
QUOTE_PULL_MIN = float(_get("QUOTE_PULL_MIN", "3"))
MAX_QUOTE_PCT = float(_get("MAX_QUOTE_PCT", "0.30"))
# Conservative tuition-run defaults (finding #74/F2/F8): at the documented
# BANKROLL=1000, MATCH_CAP_PCT=0.03 -> $30/game and GLOBAL_CAP_PCT=0.075 ->
# $75 global. The OLD defaults (0.40/0.75 -> $400/$750) were far above the
# real ~$300 account and above what HARD_STOP_DOLLARS ($50) implied was the
# leash -- the interval ledger, not the hard stop, is what actually bounds
# a single game's worst case, so it must be the honest small number.
# Arithmetic (BANKROLL * pct): 1000*0.03=30.0, 1000*0.075=75.0.
# THE OPERATOR MUST CONFIRM these before going live -- this is risk
# appetite, not something to accept on a default.
MATCH_CAP_PCT = float(_get("MATCH_CAP_PCT", "0.03"))
GLOBAL_CAP_PCT = float(_get("GLOBAL_CAP_PCT", "0.075"))
DAILY_LOSS_HALT_PCT = float(_get("DAILY_LOSS_HALT_PCT", "0.40"))
# The trading day rolls in US-Eastern at this hour, not UTC midnight (finding
# #75/F6): UTC midnight is 8pm ET -- squarely mid-slate -- so a naive
# now.date() reset hands a fresh loss budget in the middle of a losing
# night's games. All MLB/soccer games are settled by 6am ET, so rolling
# then keeps a full evening slate inside one trading day.
DAILY_ROLL_HOUR_ET = float(_get("DAILY_ROLL_HOUR_ET", "6"))
FILL_BURST_N = int(_get("FILL_BURST_N", "3"))
COOLOFF_MIN = float(_get("COOLOFF_MIN", "10"))
MAKER_MAX_CONTRACTS = int(_get("MAKER_MAX_CONTRACTS", "2"))   # hard per-quote contract ceiling (tuition run)
HARD_STOP_DOLLARS = float(_get("HARD_STOP_DOLLARS", "50"))    # cumulative realized+unrealized loss halt (mark-to-anchor)
ANCHOR_STALE_SEC = float(_get("ANCHOR_STALE_SEC", "180"))     # WITHIN ANCHOR_STALE_FARK_SEC of kickoff: pull match if even the freshest anchor rung's modifiedOn is older than this (frozen-feed guard)
ANCHOR_STALE_FARK_SEC = float(_get("ANCHOR_STALE_FARK_SEC", "7200"))  # far-from-kickoff cutoff (default 2h): beyond this the sharp total legitimately sits unchanged for hours, so tolerate an old-but-latest number and rely on the poll-success watchdog as the dead-feed guard
# Frozen-CDN detector (Finding #77): an HTTP 200 serving a byte-identical
# feed_signature (feed.py) for longer than this is treated as a dead feed —
# note_success is withheld so MAX_STALENESS_SEC's watchdog trips and pulls
# quotes, even though every individual poll "succeeded". Deliberately much
# longer than a normal tick so a genuinely calm pre-game market (numbers not
# moving, but the connection is alive) never false-trips this.
FEED_FROZEN_SEC = float(_get("FEED_FROZEN_SEC", "600"))
TOUCH_JOIN_MIN_EDGE_CENTS = float(_get("TOUCH_JOIN_MIN_EDGE_CENTS", "1.0"))       # net edge (fair - touch - maker fee) to join the crowd's best bid
TOUCH_JOIN_ALT_MIN_EDGE_CENTS = float(_get("TOUCH_JOIN_ALT_MIN_EDGE_CENTS", "1.5"))  # alt rungs: less-trusted fair demands more edge
TOUCH_JOIN_EXIT_EDGE_CENTS = float(_get("TOUCH_JOIN_EXIT_EDGE_CENTS", "0.25"))    # hold a resting touch-join (queue position) until edge decays below this
MAKER_DB_PATH = PKG_DIR / "unabated_edge_maker.duckdb"

KALSHI_API_KEY_ID = _get("KALSHI_API_KEY_ID")
KALSHI_PRIVATE_KEY_PATH = _get("KALSHI_PRIVATE_KEY_PATH")
KALSHI_BASE_URL = _get("KALSHI_BASE_URL", "https://api.elections.kalshi.com/trade-api/v2")

MARKET_DB_PATH = PKG_DIR / "unabated_edge_market.duckdb"
RESEARCH_DB_PATH = PKG_DIR / "unabated_edge_research.duckdb"
LOG_PATH = PKG_DIR / "bot.log"
KILL_FILE = PKG_DIR / ".kill"
