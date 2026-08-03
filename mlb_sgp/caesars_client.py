"""Caesars (Liberty / American Wagering, ZeroFlucs engine) HTTP transport
for the MLB SGP library.

Caesars prices arbitrary same-game combos via a plain REST call:
    POST /sb/v2/bets/details   body {channel, channelDetail, includeMultiLine,
                                     combinationSelections: [], legs: [...]}
    -> resp.parlays[0].price.decimal   == correlated SGP (ZeroFlucs)
`combinationSelections` MUST be empty — the server infers the SGP from the
same-event legs. Works LOGGED OUT. (The be-push socket.io is a login-only
bet-referral channel, NOT the quote path.)

Every request to api.americanwagering.com needs an AWS-WAF token. The token is
minted BROWSER-FREE by running AWS WAF's real challenge.js under `node` (see
caesars_waf.py — no Chromium/Playwright, ~1.5s) and CACHED to disk (short TTL,
~5 min). The token-broker validates the token (a tabs GET must return JSON)
before handing it out — an unverified WAF token 403s — so a cold cycle simply
yields no rows rather than garbage.

Headers that matter (all required): x-aws-waf-token, x-platform:cordova-desktop,
x-app-version, x-unique-device-id, a NON-headless Chrome UA, and
Origin: https://sportsbook.caesars.com.

See recon_caesars_betdetails_live.py for the end-to-end proof (3 MLB games).
"""
from __future__ import annotations
import logging

import json
import time
from dataclasses import dataclass
from pathlib import Path

from curl_cffi import requests

from mlb_sgp._shared import (RETRY_BACKGROUND, RETRY_LIVE, BookTransportError,
                             PriceCallTallyMixin, RetryProfile,
                             request_with_retry)

logger = logging.getLogger(__name__)


def _looks_like_json(resp) -> bool:
    """True when the body opens like JSON.

    Caesars' AWS-WAF challenge is served as HTTP 200 with an HTML
    interstitial, so status alone cannot tell a good response from a blocked
    one. Passed to ``request_with_retry(body_ok=...)`` so a challenge is
    RETRIED (asking again after a token mint is exactly how this book
    recovers) rather than treated as a usable response.
    """
    return (getattr(resp, "text", "") or "").strip().startswith(("{", "["))


BOOK = "caesars"
DEFAULT_STATE = "nj"
API = "https://api.americanwagering.com"
# A real (non-Headless) Chrome UA — CloudFront blocks HeadlessChrome.
REAL_UA = ("Mozilla/5.0 (Macintosh; Intel Mac OS X 10_15_7) AppleWebKit/537.36 "
           "(KHTML, like Gecko) Chrome/149.0.0.0 Safari/537.36")
APP_VERSION = "7.49.0"
# Tabs feed (URL-encoded "SCHEDULE|Games ⚾") used both to list events and to
# validate a freshly-minted token.
TABS_PATH = "v4/sports/baseball/tabs/SCHEDULE%7CGames%20%E2%9A%BE"

_THIS_DIR = Path(__file__).resolve().parent
TOKEN_CACHE = _THIS_DIR / ".caesars_token.json"
TOKEN_TTL_SEC = 240          # refresh comfortably inside the ~5 min WAF TTL


@dataclass(frozen=True)
class Event:
    event_id: str
    competition_id: str
    home_team: str
    away_team: str
    start_time: str          # ISO-8601 UTC (event.startTime)


class CaesarsClient(PriceCallTallyMixin):
    BOOK = BOOK

    def __init__(self, state: str = DEFAULT_STATE, verbose: bool = False):
        self.state = state
        self.verbose = verbose
        self.session = requests.Session(impersonate="chrome")
        self._token = ""
        self._device = ""
        self._cookies: dict = {}
        self._minted_at = 0.0
        # Wall-clock cost of the most recent token mint (0.0 until one runs).
        self.last_mint_sec = 0.0

    # --- hosts -------------------------------------------------------------- #
    @property
    def _base(self) -> str:
        return f"{API}/regions/us/locations/{self.state}/brands/czr"

    @property
    def _sb(self) -> str:
        return f"{self._base}/sb"

    def _headers(self) -> dict:
        h = {
            "User-Agent": REAL_UA,
            "Accept": "application/json, text/plain, */*",
            "Content-Type": "application/json",
            "Origin": "https://sportsbook.caesars.com",
            "Referer": "https://sportsbook.caesars.com/",
            "X-Aws-Waf-Token": self._token,
            "x-platform": "cordova-desktop",
            "x-app-version": APP_VERSION,
        }
        if self._device:
            h["x-unique-device-id"] = self._device
        return h

    # --- token broker ------------------------------------------------------- #
    def ensure_token(self) -> bool:
        """Guarantee a *validated* WAF token is loaded. Returns False if one
        can't be obtained this cycle (caller then produces no rows)."""
        if self._token and (time.time() - self._minted_at) < TOKEN_TTL_SEC:
            return True
        if self._load_cache():
            return True
        return self._mint()

    def _load_cache(self) -> bool:
        try:
            c = json.loads(TOKEN_CACHE.read_text())
        except Exception:
            return False
        if time.time() - c.get("minted_at", 0) >= TOKEN_TTL_SEC:
            return False
        self._token = c.get("token", "")
        self._device = c.get("device", "")
        self._cookies = c.get("cookies", {})
        self._minted_at = c.get("minted_at", 0)
        if not self._token or not self._validate():
            return False
        return True

    def _save_cache(self):
        try:
            TOKEN_CACHE.write_text(json.dumps({
                "token": self._token, "device": self._device,
                "cookies": self._cookies, "minted_at": self._minted_at,
            }))
        except Exception:
            pass

    def _validate(self) -> bool:
        """A verified token returns JSON from the tabs feed; an unverified one
        gets a CloudFront/WAF challenge (non-JSON)."""
        try:
            r = self.session.get(f"{self._sb}/{TABS_PATH}", headers=self._headers(),
                                  cookies=self._cookies, impersonate="chrome", timeout=20)
            return r.status_code == 200 and r.text.strip().startswith("{")
        except Exception:
            return False

    def _mint(self) -> bool:
        """Browser-free mint of the aws-waf-token + device id.

        Runs AWS WAF's real challenge.js under `node` (no Chromium/Playwright)
        via caesars_waf.mint_token_browser_free, then validates the token.
        Sub-second in practice and works in the bot venvs (node on PATH +
        curl_cffi only). `last_mint_sec` carries the most recent mint's
        wall-clock cost — the latency an RFQ pays on a cold token, which #50's
        pre-warm work needs to measure.
        """
        try:
            from mlb_sgp.caesars_waf import mint_token_browser_free
        except ImportError:
            from caesars_waf import mint_token_browser_free  # cwd=mlb_sgp path
        def _record_duration(seconds: float) -> None:
            self.last_mint_sec = seconds

        for attempt in range(3):
            try:
                tok, dev = mint_token_browser_free(state=self.state,
                                                   on_duration=_record_duration)
            except Exception as e:
                # Auth-stage surprise: WARNING (rare, and a dead WAF token is
                # how this book rots silently — see issue #32's meta-bug).
                logger.warning("caesars: WAF mint error (attempt %d/3): %r",
                               attempt + 1, e)
                continue
            if not tok:
                logger.warning("caesars: WAF mint returned no token "
                               "(attempt %d/3, %.2fs)",
                               attempt + 1, self.last_mint_sec)
                continue
            self._token = tok
            self._device = dev
            self._cookies = {"aws-waf-token": tok}
            self._minted_at = time.time()
            if self._validate():
                self._save_cache()
                logger.info("caesars: WAF token minted+validated browser-free "
                            "(attempt %d/3, mint %.2fs)",
                            attempt + 1, self.last_mint_sec)
                return True
            logger.warning("caesars: WAF token minted but NOT validated "
                           "(attempt %d/3) — WAF reject or IP throttle",
                           attempt + 1)
        return False

    # --- REST --------------------------------------------------------------- #
    # NOTE: the module-level _looks_like_json (top of this file) is what makes
    # a 200-with-HTML AWS-WAF challenge retryable instead of an instant failure.
    def _get_json(self, url: str, stage: str,
                  profile: RetryProfile = RETRY_BACKGROUND):
        """GET + decode, or raise ``BookTransportError``.

        Returns None for a 404 at the ``structure`` stage only — one delisted
        event is a skip, but a 404 on the tabs FEED means the endpoint moved
        and the book is down. A WAF challenge comes back 200 with HTML, so
        "not JSON" is a transport failure too: before issue #33 that quietly
        became an empty slate.

        Issue #34 replaced this method's bespoke 4x1.5s loop with the shared
        helper. The behavior change worth knowing: that loop retried EVERY
        non-JSON outcome including a 403 auth gate, burning 6s before
        admitting the book was dead. 4xx now fails on the first attempt; only
        5xx/429/connection errors and the 200-with-HTML WAF challenge (via
        ``body_ok``) are retried.
        """
        r = request_with_retry(
            lambda: self.session.get(url, headers=self._headers(),
                                     cookies=self._cookies,
                                     impersonate="chrome", timeout=25),
            profile=profile, book=BOOK, stage=stage,
            body_ok=_looks_like_json, detail="no usable JSON")
        status = getattr(r, "status_code", 200)
        if status == 404 and stage == "structure":
            return None
        if status == 200 and _looks_like_json(r):
            return r.json()
        raise BookTransportError(
            BOOK, stage, status_code=status,
            detail="response was not usable JSON")

    def list_events(self, profile: RetryProfile = RETRY_BACKGROUND) -> list[Event]:
        if not self.ensure_token():
            raise BookTransportError(
                BOOK, "auth",
                detail="AWS-WAF token could not be minted/validated")
        tj = self._get_json(f"{self._sb}/{TABS_PATH}", "events", profile)
        if not tj:
            return []
        out: list[Event] = []
        for comp in tj.get("competitions", []):
            cid = comp.get("id")
            cname = (comp.get("name") or "")
            # baseball tabs also carry KBO/NPB/CPBL/college/MiLB — keep MLB only.
            if cname.strip().upper() != "MLB":
                continue
            for ev in comp.get("events", []):
                # Skip in-play games: their markets are renamed "... Live" and
                # the lines move mid-game. The bots/dashboard price pregame only.
                # NOTE: `started` is the in-play signal; `tradedInPlay` is a
                # capability flag (True for ALL games that support live betting,
                # including pregame) — do NOT filter on it or you drop everything.
                if ev.get("started"):
                    continue
                md = ev.get("metadata", {}) or {}
                home = md.get("homeTeamName")
                away = md.get("awayTeamName")
                if not (home and away):
                    continue
                out.append(Event(
                    event_id=ev.get("id"),
                    competition_id=ev.get("competitionId") or cid,
                    home_team=home, away_team=away,
                    start_time=ev.get("startTime", "") or "",
                ))
        return out

    def fetch_event(self, event_id: str,
                    profile: RetryProfile = RETRY_BACKGROUND) -> dict | None:
        """Full event detail (event.keyMarketGroups[].markets[] incl. alts).

        None means Caesars no longer lists this event (404); transport
        failures raise instead — see issue #33.
        """
        if not self.ensure_token():
            raise BookTransportError(
                BOOK, "auth",
                detail="AWS-WAF token could not be minted/validated")
        j = self._get_json(f"{self._sb}/v4/events/{event_id}", "structure",
                           profile)
        if not j:
            return None
        return j.get("event") or j

    def price_combo(self, legs: list[dict]) -> dict | None:
        """POST /sb/v2/bets/details for a same-game combo; read the correlated
        decimal from parlays[0]. Returns {"decimal","american"} or None."""
        body = {
            "channel": "desktop",
            "channelDetail": "cordova-desktop",
            "includeMultiLine": False,
            "combinationSelections": [],
            "legs": legs,
        }
        try:
            # RETRY_LIVE on every path — price calls are the fan-out surface.
            r = request_with_retry(
                lambda: self.session.post(
                    f"{self._sb}/v2/bets/details", headers=self._headers(),
                    cookies=self._cookies, json=body, impersonate="chrome",
                    timeout=30),
                profile=RETRY_LIVE, book=BOOK, stage="price",
                body_ok=_looks_like_json)
        except BookTransportError:
            # Per-combo price failures stay row-drops (they are partial), but
            # they must be TALLIED so an all-fail cycle still gets a verdict.
            self.price_calls.record(False)
            return None
        if (getattr(r, "status_code", 200) != 200
                or not (getattr(r, "text", "") or "").strip().startswith("{")):
            self.price_calls.record(False)
            return None
        try:
            parlays = r.json().get("parlays") or []
        except Exception:
            self.price_calls.record(False)
            return None
        if not parlays:
            self.price_calls.record(False)
            return None
        price = (parlays[0].get("price") or {})
        dec = price.get("decimal") if price.get("decimal") is not None else price.get("d")
        if not dec or dec <= 1.0:
            self.price_calls.record(False)
            return None
        am = price.get("american")
        try:
            am = int(str(am).replace("+", "")) if am is not None else None
        except (TypeError, ValueError):
            am = None
        self.price_calls.record(True)
        return {"decimal": float(dec), "american": am}
