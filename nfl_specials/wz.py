"""Wagerzon side: read the week's NFL trifectas/superfectas and place them.

Reads:  ActiveLeaguesHelper (find "NFL WEEK <n> - SPECIALS"), then
        NewScheduleHelper for that league (JSON; needs the XHR header).
Places: wagerzon_odds/single_placer.place_single — ConfirmWagerHelper
        preflight with a $0.01 win-drift check, then PostWagerMultipleHelper.
        A special is a one-sided prop on the "home" slot, so it goes in as
        Play=5 with no points (verified against the preflight echo).
Credentials: WAGERZON[_X]_USERNAME / _PASSWORD from bet_logger/.env (the main
checkout's copy is found from a worktree too). No DB writes here.
"""
from __future__ import annotations

import logging
import re
import sys
import time
from dataclasses import dataclass
from typing import Callable, TypeVar

import requests
from dotenv import load_dotenv

from nfl_specials import config


def _load_wagerzon_env() -> None:
    """bet_logger/.env lives in the main checkout; a worktree walks up to it."""
    for candidate in (config.REPO_ROOT, *config.REPO_ROOT.parents):
        env_file = candidate / "bet_logger" / ".env"
        if env_file.exists():
            load_dotenv(env_file)
            return


_load_wagerzon_env()
# The Wagerzon modules import their siblings by bare name.
sys.path.insert(0, str(config.REPO_ROOT / "wagerzon_odds"))
import single_placer  # noqa: E402
import wagerzon_auth  # noqa: E402
from scraper_v2 import fetch_active_leagues  # noqa: E402
from wagerzon_accounts import list_accounts  # noqa: E402

log = logging.getLogger("nfl_specials.wz")

SCHEDULE_URL = config.WZ_BASE_URL + "/wager/NewScheduleHelper.aspx?WT=0&lg={league_id}"
FECTA_WORDS = re.compile(r"\b(TRIFECTA|SUPERFECTA)\b")
# Reads only. Wagerzon's login page occasionally stalls past the 15 s read
# timeout and answers in under a second on the next try (2026-10-02).
READ_ATTEMPTS = 3
READ_RETRY_PAUSE_SECONDS = 5
T = TypeVar("T")


@dataclass(frozen=True)
class WzSpecial:
    wz_game_id: int        # Wagerzon idgm — what a bet is placed against
    rotation: int          # Wagerzon rotation number (hnum)
    description: str       # e.g. "SEAHAWKS TRIFECTA (1Q, 1H & GM)"
    wz_american: int       # posted price


def account_labels() -> list[str]:
    return [account.label for account in list_accounts()]


def _session(account_label: str | None = None):
    accounts = list_accounts()
    if not accounts:
        raise RuntimeError("no Wagerzon accounts configured (WAGERZON_USERNAME / _PASSWORD)")
    account = next((a for a in accounts if a.label == account_label), accounts[0])
    return wagerzon_auth.get_session(account)


def specials_league_id(session) -> int:
    catalog = fetch_active_leagues(session) or []
    for row in catalog:
        if re.match(config.WZ_SPECIALS_DESCRIPTION_RE, (row.get("Description") or "").strip()):
            return int(row["IdLeague"])
    raise RuntimeError("no 'NFL WEEK <n> - SPECIALS' league in Wagerzon's catalog")


def _retry_reads(read: Callable[[], T]) -> T:
    """Retry a READ on timeouts / dropped connections. Never wrap a placement:
    a retried submission can place the bet twice."""
    for attempt in range(1, READ_ATTEMPTS + 1):
        try:
            return read()
        except (requests.Timeout, requests.ConnectionError) as exc:
            if attempt == READ_ATTEMPTS:
                raise
            log.warning("Wagerzon read failed (%s), attempt %d/%d; retrying in %ds",
                        type(exc).__name__, attempt, READ_ATTEMPTS, READ_RETRY_PAUSE_SECONDS)
            time.sleep(READ_RETRY_PAUSE_SECONDS)
    raise AssertionError("unreachable")


def fetch_fecta_specials() -> list[WzSpecial]:
    """Every trifecta/superfecta currently posted in the week's specials."""
    return _retry_reads(_fetch_fecta_specials_once)


def _fetch_fecta_specials_once() -> list[WzSpecial]:
    session = _session()
    url = SCHEDULE_URL.format(league_id=specials_league_id(session))
    resp = session.get(url, timeout=30, headers={"Accept": "application/json",
                                                 "X-Requested-With": "XMLHttpRequest"})
    resp.raise_for_status()
    leagues = (resp.json().get("result") or {}).get("listLeagues") or [[]]
    specials = []
    for league in leagues[0]:
        for game in league.get("Games", []):
            description = (game.get("htm") or "").strip()
            if not FECTA_WORDS.search(description.upper()):
                continue
            line = (game.get("GameLines") or [{}])[0] or {}
            posted = line.get("oddsh") or line.get("odds")
            if not posted:
                continue
            specials.append(WzSpecial(wz_game_id=int(game["idgm"]), rotation=int(game["hnum"]),
                                      description=description, wz_american=int(posted)))
    return specials


def place_fecta(account_label: str, special: WzSpecial, risk: float) -> dict:
    """Place `risk` dollars on `special` at its posted price.

    Returns single_placer's result: status 'placed' (with ticket_number) or
    'price_moved' / 'rejected' / 'auth_error' / 'network_error' / 'orphaned'.
    A changed price is refused by the placer's preflight, never filled.
    """
    bet = {
        "idgm": special.wz_game_id,
        "play": config.WZ_SPECIAL_PLAY,
        "line": None,
        "american_odds": special.wz_american,
        "wz_odds_at_place": special.wz_american,
        "actual_size": risk,
        "bet_hash": f"nfl-fecta-{special.wz_game_id}",
        "market": "nfl_fecta",
        "bet_on": special.description,
    }
    return single_placer.place_single(account_label, bet)
