"""BetMGM (Entain "CDS" platform) HTTP transport for the MLB SGP library.

Thin client over a curl_cffi Chrome-TLS session. Three responsibilities:

  1. ``accessid()``  — harvest the static, per-state ``publicAccessId`` that
     gates every CDS read. The unlock is the ``x-bwin-sports-api: prod``
     header on the ``clientconfig`` endpoint; without it the endpoint
     returns the host-app shell (no accessid). The id is long-lived and
     identical across re-harvests, so we cache it and only refresh on a
     ``400 "Access id ..."`` from the API.

  2. ``list_events()`` — the cheap "Gridable" fixtures call: enough to map
     BetMGM fixtures to canonical games (teams + start time). One request.

  3. ``fetch_markets(fixture_id)`` / ``price_picks(...)`` — the per-game
     full market tree (``offerMapping=All``, all alt run lines + alt totals)
     and the Angstrom bet-builder price POST.

Geo note: confirmed reachable from a non-legal-state IP (CA/ID) with plain
curl_cffi. Only the SPA bootstrap is IP-geo-gated, which we bypass by going
straight to the JSON endpoints. NJ omits the pubid in its config, so the
default state is ``pa``.

See ``recon_betmgm_sgp.py`` (``mgm`` mode) for the end-to-end proof.
"""
from __future__ import annotations

import re
import uuid
from dataclasses import dataclass

from curl_cffi import requests

from mlb_sgp._shared import (RETRY_BACKGROUND, RETRY_LIVE, BookTransportError,
                             PriceCallTallyMixin, RetryProfile, check_response,
                             json_or_raise, request_with_retry)

BOOK = "betmgm"
MLB_SPORT_ID = "23"          # Entain CDS sport id for baseball
COUNTRY = "US"
DEFAULT_STATE = "pa"         # NJ omits the pubid; PA returns a clean static id


@dataclass(frozen=True)
class Event:
    """One BetMGM MLB fixture, reduced to what match_events needs."""
    event_id: str
    home_team: str
    away_team: str
    start_time: str          # ISO-8601 UTC (fixture.startDate)


class BetMGMClient(PriceCallTallyMixin):
    BOOK = BOOK

    def __init__(self, state: str = DEFAULT_STATE, verbose: bool = False):
        self.state = state
        self.verbose = verbose
        self.session = requests.Session(impersonate="chrome")
        self.session.headers.update({
            "Accept": "application/json, text/plain, */*",
            "Accept-Language": "en-US,en;q=0.9",
            "Referer": f"https://www.{state}.betmgm.com/",
        })
        self._accessid = ""

    # --- hosts / common params ------------------------------------------- #
    @property
    def _host(self) -> str:
        return f"www.{self.state}.betmgm.com"

    @property
    def _base(self) -> str:
        return f"https://{self._host}/cds-api"

    def _common_q(self) -> str:
        return (f"x-bwin-accessid={self.accessid()}"
                f"&lang=en&country={COUNTRY}&userCountry={COUNTRY}")

    # --- accessid harvest ------------------------------------------------ #
    def accessid(self, refresh: bool = False,
                 profile: RetryProfile = RETRY_BACKGROUND) -> str:
        """Return the cached static accessid, harvesting it on first use.

        The ``x-bwin-sports-api: prod`` header is the unlock — it flips the
        clientconfig response from the 11 KB host-app shell to the ~97 KB
        sports config carrying ``msApp.publicAccessId``.

        Every CDS read is gated on this id, so a failed harvest is a
        session-level AUTH failure, not an empty slate: it raises
        ``BookTransportError`` rather than returning "" and letting every
        downstream call quietly produce no rows (issue #33).
        """
        if self._accessid and not refresh:
            return self._accessid
        burl = f"http%3A%2F%2F{self._host}%2Fen%2Fsports%2Fbaseball-23"
        url = (f"https://{self._host}/en/api/clientconfig"
               f"?browserUrl={burl}&x-from-product=host-app")
        def _get():
            return self.session.get(url, headers={
                "Referer": f"https://{self._host}/en/sports/baseball-23",
                "x-bwin-browser-url": burl,
                "x-from-product": "host-app",
                "x-bwin-sports-api": "prod",
            }, timeout=25)

        # Every CDS read is gated on this call, so retrying it is the single
        # highest-leverage retry at this book: one blip here loses the whole
        # cycle, not one event.
        r = request_with_retry(_get, profile=profile, book=BOOK, stage="auth",
                               detail="accessid harvest failed")
        check_response(BOOK, "auth", r)
        ids = re.findall(r'"publicAccessId":"([^"]+)"', getattr(r, "text", "") or "")
        if not ids:
            raise BookTransportError(
                BOOK, "auth",
                detail=(f"no publicAccessId in clientconfig for state="
                        f"{self.state} (NJ omits it — try pa/mi/co/va/tn)"))
        self._accessid = ids[0]
        return self._accessid

    # --- event listing --------------------------------------------------- #
    def list_events(self, take: int = 120,
                    profile: RetryProfile = RETRY_BACKGROUND) -> list[Event]:
        """List upcoming MLB fixtures (cheap Gridable call) as Event records.

        ``sportIds=23`` is *all* baseball (KBO, NPB, CPBL, college, Mexican
        leagues, ...), so we filter to the MLB competition (id 75 / name
        "MLB") client-side. ``take`` is generous because foreign overnight
        games pad the unfiltered list.
        """
        self.accessid(profile=profile)   # raises BookTransportError if it can't harvest
        url = (f"{self._base}/bettingoffer/fixtures?{self._common_q()}"
               f"&fixtureTypes=Standard&state=Latest&offerMapping=Filtered"
               f"&offerCategories=Gridable&fixtureCategories=Gridable"
               f"&sportIds={MLB_SPORT_ID}&skip=0&take={take}&sortBy=Tags")
        r = request_with_retry(lambda: self.session.get(url, timeout=30),
                               profile=profile, book=BOOK, stage="events")
        check_response(BOOK, "events", r)
        out: list[Event] = []
        for fx in json_or_raise(BOOK, "events", r).get("fixtures", []):
            if not _is_mlb(fx):
                continue
            home, away = _split_home_away(fx)
            if not (home and away):
                continue
            out.append(Event(
                event_id=str(fx.get("id")),
                home_team=home,
                away_team=away,
                start_time=fx.get("startDate", "") or "",
            ))
        return out

    # --- full market tree (alts) ----------------------------------------- #
    def fetch_markets(self, fixture_id: str,
                      profile: RetryProfile = RETRY_BACKGROUND) -> list[dict]:
        """Return the full ``optionMarkets`` list for one fixture.

        ``offerMapping=All`` returns every alt run line + alt total + F5
        market (the Gridable listing only carries the main line). Refreshes
        the accessid once on a 400 gate error; a 404 means the fixture is
        gone (skip it), anything else non-200 raises.
        """
        self.accessid(profile=profile)   # raises BookTransportError if it can't harvest
        r = self._get_fixture_markets(fixture_id, profile)
        if _is_accessid_gate(r):
            self.accessid(refresh=True, profile=profile)
            r = self._get_fixture_markets(fixture_id, profile)
        if not check_response(BOOK, "structure", r, allow_404=True):
            return []
        fxs = json_or_raise(BOOK, "structure", r).get("fixtures", [])
        return fxs[0].get("optionMarkets", []) if fxs else []

    def _get_fixture_markets(self, fixture_id: str, profile: RetryProfile):
        """One raw GET of the full market tree (retried per ``profile``).

        Split out so fetch_markets can re-issue it once after refreshing a
        stale accessid without a loop whose exit branch is unreachable. That
        accessid re-issue is an APPLICATION-level retry on a 400 gate — the
        transport retry inside this call is a separate, orthogonal concern.
        """
        url = (f"{self._base}/bettingoffer/fixtures?{self._common_q()}"
               f"&fixtureIds={fixture_id}&offerMapping=All&state=Latest")
        return request_with_retry(lambda: self.session.get(url, timeout=30),
                                  profile=profile, book=BOOK, stage="structure")

    # --- bet-builder price ----------------------------------------------- #
    def price_picks(self, fixture_id: str,
                    legs: list[tuple[int, int]]) -> dict | None:
        """POST a same-game combo (legs share one pickGroupId) and read the
        correlated Angstrom price.

        ``legs`` is a list of ``(market_id, option_id)`` tuples. Returns
        ``{"decimal", "american", "provider"}`` or ``None`` when BetMGM
        declines to price the combo (combo-prevention, suspended, etc.).
        """
        pg = str(uuid.uuid4())
        body = {
            "tv1Picks": [
                {"fixtureId": str(fixture_id), "gameId": mid, "resultId": oid,
                 "useLiveFallback": False, "pickGroupId": pg}
                for (mid, oid) in legs
            ],
            "tv2Picks": [],
        }
        url = f"{self._base}/bettingoffer/picks?{self._common_q()}"
        try:
            # RETRY_LIVE on every path — price calls are the fan-out surface.
            r = request_with_retry(
                lambda: self.session.post(url, json=body, timeout=30),
                profile=RETRY_LIVE, book=BOOK, stage="price")
        except BookTransportError as e:
            # Per-combo price failures stay row-drops (they are partial), but
            # they must be TALLIED so an all-fail cycle still gets a verdict.
            self.price_calls.record(False)
            if self.verbose:
                print(f"  [mgm] price_picks error: {e!r}")
            return None
        if getattr(r, "status_code", 200) != 200:
            self.price_calls.record(False)
            return None
        try:
            groups = r.json().get("betBuilderPricingGroups", {}) or {}
        except Exception:
            self.price_calls.record(False)
            return None
        if not groups:
            self.price_calls.record(False)
            return None
        grp = next(iter(groups.values()))
        odds = grp.get("odds") or {}
        dec = odds.get("odds")
        if not dec or dec <= 1.0:
            self.price_calls.record(False)
            return None
        self.price_calls.record(True)
        return {
            "decimal": float(dec),
            "american": odds.get("americanOdds"),
            "provider": grp.get("providerId"),
        }


# --------------------------------------------------------------------------- #
# Module helpers
# --------------------------------------------------------------------------- #
def _is_accessid_gate(resp) -> bool:
    """True for BetMGM's ``400 "Access id ..."`` stale-credential response —
    the one non-200 worth retrying rather than raising on."""
    if getattr(resp, "status_code", 200) != 400:
        return False
    return "ccess id" in (getattr(resp, "text", "") or "").lower()


def _is_mlb(fixture: dict) -> bool:
    """True only for Major League Baseball fixtures (competition name 'MLB').

    Filters out KBO / NPB / CPBL / college / Mexican-league baseball that
    also live under baseball's sportId.
    """
    comp = fixture.get("competition") or {}
    cname = comp.get("name")
    if isinstance(cname, dict):
        cname = cname.get("value", "")
    return (cname or "").strip().upper() == "MLB"


def _split_home_away(fixture: dict) -> tuple[str, str]:
    """Return (home_team, away_team) full names for a fixture.

    BetMGM fixture names are ``"<Away> at <Home>"``; we parse that first and
    fall back to the participant order (away is listed first). Returns
    ("", "") when neither yields two teams.
    """
    name = ""
    nm = fixture.get("name")
    if isinstance(nm, dict):
        name = nm.get("value", "") or ""
    if " at " in name:
        away, home = name.split(" at ", 1)
        return home.strip(), away.strip()

    teams = [
        p for p in fixture.get("participants", [])
        if (p.get("properties") or {}).get("type") != "Player"
    ]
    if len(teams) >= 2:
        def _pname(p):
            n = p.get("name") or {}
            return (n.get("value") or "").strip() if isinstance(n, dict) else ""
        away, home = _pname(teams[0]), _pname(teams[1])
        if home and away:
            return home, away
    return "", ""
