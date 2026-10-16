"""Fixed constants for the NFL trifecta/superfecta pricer.

Deployment-specific values come from env vars read ONCE here; every pricing
and sizing constant is a plain named number so the whole model fits on one
screen.
"""
from __future__ import annotations

import os
import re
from pathlib import Path

PACKAGE_DIR = Path(__file__).resolve().parent
REPO_ROOT = PACKAGE_DIR.parent

# App state (latest board, quote history, placed bets, settings).
STATE_DB_PATH = Path(os.environ.get("NFL_SPECIALS_DB", str(PACKAGE_DIR / "nfl_specials.duckdb")))

# Local web page. Loopback only: it can place real Wagerzon bets.
APP_HOST = "127.0.0.1"
APP_PORT = 8096

# --- Server (Oracle VM) ------------------------------------------------------
# The HTTPS port `tailscale serve` publishes the page on, so a browser on the
# tailnet sends Host "<tailnet name>:8443".
TAILNET_HTTPS_PORT = 8443
# dk_relay's port on the Mac's Tailscale address.
DK_RELAY_PORT = 8097

def parse_extra_allowed_hosts(raw: str | None) -> tuple[str, ...]:
    """NFL_SPECIALS_EXTRA_ALLOWED_HOSTS ("a.ts.net,b.ts.net") as host names; empty or unset is ()."""
    names = []
    for entry in (raw or "").split(","):
        name = entry.strip()
        if name and name not in names:
            names.append(name)
    return tuple(names)


def parse_dk_proxy(raw: str | None) -> str | None:
    """NFL_SPECIALS_DK_PROXY ("http://<Mac tailnet IP>:8097") or None when unset."""
    value = (raw or "").strip()
    if not value:
        return None
    if not re.fullmatch(r"http://[^/\s]+:\d+", value):
        raise ValueError(f"NFL_SPECIALS_DK_PROXY: expected http://<host>:<port>, got {value!r}")
    return value


# The VM's tailnet name, which `tailscale serve` forwards in the Host header.
EXTRA_ALLOWED_HOSTS = parse_extra_allowed_hosts(os.environ.get("NFL_SPECIALS_EXTRA_ALLOWED_HOSTS"))
# DraftKings answers "Access Denied" to data-center IPs, so on a server its
# requests go through a proxy at home. None: straight to DraftKings.
DK_PROXY_URL = parse_dk_proxy(os.environ.get("NFL_SPECIALS_DK_PROXY"))
# On a server, the env file holding the Wagerzon login (the bets service's .env,
# which also holds other books' keys). None: bet_logger/.env.
WAGERZON_ENV_PATH = os.environ.get("NFL_SPECIALS_WAGERZON_ENV") or None

# --- Wagerzon ---------------------------------------------------------------
WZ_BASE_URL = "https://backend.wagerzon.com"
# "NFL WEEK 4 - SPECIALS": the week number changes, so the league id is
# resolved from the live catalog by this pattern on every scrape.
WZ_SPECIALS_DESCRIPTION_RE = r"^NFL WEEK \d+ - SPECIALS$"

# --- Sizing -------------------------------------------------------------------
DEFAULT_BANKROLL = 1000.0
DEFAULT_KELLY_FRACTION = 0.25
# Measured 2026-10-02 with preview-only ConfirmWagerHelper calls on a special:
# $15-$19 -> MINWAGERONLINE, $20 accepted.
WZ_MIN_STAKE = 20.0
# Wagerzon's maximum risk per special (Cal, 2026-10-02).
WZ_MAX_STAKE = 250.0
