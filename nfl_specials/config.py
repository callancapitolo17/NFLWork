"""Fixed constants for the NFL trifecta/superfecta pricer.

Deployment-specific values come from env vars read ONCE here; every pricing
and sizing constant is a plain named number so the whole model fits on one
screen.
"""
from __future__ import annotations

import os
from pathlib import Path

PACKAGE_DIR = Path(__file__).resolve().parent
REPO_ROOT = PACKAGE_DIR.parent

# App state (latest board, quote history, placed bets, settings).
STATE_DB_PATH = Path(os.environ.get("NFL_SPECIALS_DB", str(PACKAGE_DIR / "nfl_specials.duckdb")))

# Local web page. Loopback only: it can place real Wagerzon bets.
APP_HOST = "127.0.0.1"
APP_PORT = 8096

# DraftKings prices SGPs only to a real, non-headless Chrome (issue #102), so
# calculateBets goes through dk_price_sidecar/ (loopback service owning that
# Chrome). Start it with dk_price_sidecar/run.sh before refreshing.
DK_SIDECAR_URL = os.environ.get("DK_SIDECAR_URL", "http://127.0.0.1:8095")

# --- Wagerzon ---------------------------------------------------------------
WZ_BASE_URL = "https://backend.wagerzon.com"
# "NFL WEEK 4 - SPECIALS": the week number changes, so the league id is
# resolved from the live catalog by this pattern on every scrape.
WZ_SPECIALS_DESCRIPTION_RE = r"^NFL WEEK \d+ - SPECIALS$"
# A special is a one-sided prop on the "home" slot; Wagerzon's preflight
# echoes Play=5 (home moneyline) with the posted odds.
WZ_SPECIAL_PLAY = 5

# --- Sizing -------------------------------------------------------------------
DEFAULT_BANKROLL = 1000.0
DEFAULT_KELLY_FRACTION = 0.25
# Measured 2026-10-02 with preview-only ConfirmWagerHelper calls on a special:
# $15-$19 -> MINWAGERONLINE, $20 accepted.
WZ_MIN_STAKE = 20.0
# Wagerzon's maximum risk per special (Cal, 2026-10-02).
WZ_MAX_STAKE = 250.0
