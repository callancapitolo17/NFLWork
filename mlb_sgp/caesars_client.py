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

import json
import time
from dataclasses import dataclass
from pathlib import Path

from curl_cffi import requests

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


class CaesarsClient:
    def __init__(self, state: str = DEFAULT_STATE, verbose: bool = False):
        self.state = state
        self.verbose = verbose
        self.session = requests.Session(impersonate="chrome")
        self._token = ""
        self._device = ""
        self._cookies: dict = {}
        self._minted_at = 0.0

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
        ~1.5s and works in the bot venvs (node on PATH + curl_cffi only).
        """
        try:
            from mlb_sgp.caesars_waf import mint_token_browser_free
        except ImportError:
            from caesars_waf import mint_token_browser_free  # cwd=mlb_sgp path
        for attempt in range(3):
            try:
                tok, dev = mint_token_browser_free(state=self.state, verbose=self.verbose)
            except Exception as e:
                if self.verbose:
                    print(f"  [czr] mint error (attempt {attempt+1}): {e!r}")
                continue
            if not tok:
                if self.verbose:
                    print(f"  [czr] node mint returned no token (attempt {attempt+1})")
                continue
            self._token = tok
            self._device = dev
            self._cookies = {"aws-waf-token": tok}
            self._minted_at = time.time()
            if self._validate():
                self._save_cache()
                if self.verbose:
                    print(f"  [czr] token minted+validated browser-free (attempt {attempt+1})")
                return True
            if self.verbose:
                print(f"  [czr] token minted but NOT validated (attempt {attempt+1}) "
                      f"— WAF reject or IP throttle")
        return False

    # --- REST --------------------------------------------------------------- #
    def _get_json(self, url: str, tries: int = 4):
        for _ in range(tries):
            try:
                r = self.session.get(url, headers=self._headers(), cookies=self._cookies,
                                     impersonate="chrome", timeout=25)
                if r.status_code == 200 and r.text.strip().startswith(("{", "[")):
                    return r.json()
            except Exception:
                pass
            time.sleep(1.5)
        return None

    def list_events(self) -> list[Event]:
        if not self.ensure_token():
            return []
        tj = self._get_json(f"{self._sb}/{TABS_PATH}")
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

    def fetch_event(self, event_id: str) -> dict | None:
        """Full event detail (event.keyMarketGroups[].markets[] incl. alts)."""
        if not self.ensure_token():
            return None
        j = self._get_json(f"{self._sb}/v4/events/{event_id}")
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
            r = self.session.post(f"{self._sb}/v2/bets/details", headers=self._headers(),
                                  cookies=self._cookies, json=body, impersonate="chrome",
                                  timeout=30)
        except Exception:
            return None
        if r.status_code != 200 or not r.text.strip().startswith("{"):
            return None
        try:
            parlays = r.json().get("parlays") or []
        except Exception:
            return None
        if not parlays:
            return None
        price = (parlays[0].get("price") or {})
        dec = price.get("decimal") if price.get("decimal") is not None else price.get("d")
        if not dec or dec <= 1.0:
            return None
        am = price.get("american")
        try:
            am = int(str(am).replace("+", "")) if am is not None else None
        except (TypeError, ValueError):
            am = None
        return {"decimal": float(dec), "american": am}
