# wagerzon_odds/parlay_placer.py
"""Wagerzon parlay placer — pure REST, no browser.

See docs/superpowers/specs/2026-04-26-mlb-parlay-auto-placement-design.md
for the full design. This module provides:
    - ParlaySpec / Leg / PlacementResult dataclasses
    - encode_sel / encode_detail_data leg-encoding helpers
    - _get_session / _clear_session session management
    - place_parlays(specs, dry_run=False) — top-level entry point
"""
from __future__ import annotations
from dataclasses import dataclass
from typing import Optional
import os
import re
import requests

WAGERZON_BASE_URL = "https://backend.wagerzon.com"

_CACHED_SESSION: Optional[requests.Session] = None


@dataclass
class Leg:
    """One leg of a parlay.

    play codes: 0=away spread, 1=home spread, 2=over, 3=under,
                4=away ML, 5=home ML
    points: signed line value (negative for fav spread / over)
    odds: integer American odds (e.g. 117 = +117, -130 = -130)
    pitcher: 0 = action / 3 = listed (mirrors what ConfirmWagerHelper expects)
    """
    idgm: int
    play: int
    points: float
    odds: int
    pitcher: int = 0


@dataclass
class ParlaySpec:
    parlay_hash: str
    legs: list[Leg]
    amount: float
    expected_win: float
    expected_risk: float


@dataclass
class PlacementResult:
    parlay_hash: str
    status: str
    ticket_number: Optional[str] = None
    idwt: Optional[int] = None
    actual_win: Optional[float] = None
    actual_risk: Optional[float] = None
    error_msg: str = ""
    error_msg_key: str = ""
    raw_response: Optional[str] = None


def encode_sel(legs: list[Leg]) -> str:
    """Encode legs as Wagerzon's `sel` parameter: play_idgm_points_odds, ..."""
    parts = []
    for leg in legs:
        # points printed without trailing .0 if integer
        pts = int(leg.points) if leg.points == int(leg.points) else leg.points
        parts.append(f"{leg.play}_{leg.idgm}_{pts}_{leg.odds}")
    return ",".join(parts)


def encode_detail_data(legs: list[Leg], amount: float) -> list[dict]:
    """Build the `detailData` JSON array for ConfirmWagerHelper / PostWager."""
    return [
        {
            "Amount": str(int(amount)) if amount == int(amount) else str(amount),
            "RiskWin": 0,
            "TeaserPointsPurchased": 0,
            "IdGame": leg.idgm,
            "Play": leg.play,
            "Pitcher": leg.pitcher,
            "Points": {
                "BuyPoints": 0,
                "BuyPointsDesc": "",
                "LineDesc": "",
                "selected": True,
            },
        }
        for leg in legs
    ]


def _get_session() -> requests.Session:
    """Return a logged-in Wagerzon session, reusing a cached one if present.

    Mirrors the login pattern in wagerzon_odds/parlay_pricer.py: GET base URL,
    if not already authenticated, parse __VIEWSTATE etc. and form-POST
    Account/Password.

    Resets via _clear_session() (called from auth_error retry path).
    """
    global _CACHED_SESSION
    if _CACHED_SESSION is not None:
        return _CACHED_SESSION

    session = requests.Session()
    resp = session.get(WAGERZON_BASE_URL, timeout=15)
    resp.raise_for_status()
    if "NewSchedule" in resp.url or "Welcome" in resp.url:
        _CACHED_SESSION = session
        return session

    html = resp.text
    fields = {}
    for name in ("__VIEWSTATE", "__VIEWSTATEGENERATOR", "__EVENTVALIDATION",
                 "__EVENTTARGET", "__EVENTARGUMENT"):
        m = re.search(rf'(?:name|id)="{name}"[^>]*value="([^"]*)"', html)
        if m:
            fields[name] = m.group(1)
    fields["Account"] = os.environ["WAGERZON_USERNAME"]
    fields["Password"] = os.environ["WAGERZON_PASSWORD"]
    fields["BtnSubmit"] = ""

    resp = session.post(WAGERZON_BASE_URL, data=fields, timeout=15)
    resp.raise_for_status()
    _CACHED_SESSION = session
    return session


def _clear_session() -> None:
    """Force the next _get_session() call to re-login."""
    global _CACHED_SESSION
    _CACHED_SESSION = None
