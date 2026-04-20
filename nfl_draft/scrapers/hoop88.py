"""Hoop88 NFL Draft markets adapter.

Hoop88 posts NFL Draft futures under SportType=FOOTBALL,
SportSubType='NFL Draft 2026', Period='Prop' (PeriodNumber=0). Each
sub-market (10th Overall Pick, 1st Cornerback, etc.) is a separate
`propDescription` value that has to be queried individually -- the API does
NOT return all props in one call.

Discovery path (Phase 1 recon):
  1. /cloud/api/League/Get_SportsLeagues enumerates every sport/league.
     NFL Draft markets appear as entries with SportSubType='NFL Draft 2026'
     and SportSubType2 = the per-prop name ('10th Overall Pick',
     '1st Cornerback', 'Mr Irrelevant Position', ...).
  2. /cloud/api/Lines/Get_LeagueLines2 accepts that SportSubType2 as
     `propDescription` and returns Lines[*] = [metadata_string,
     [Contestant, ...]].

As of 2026-04-20 the account-visible board contains 13 markets:
  - "2nd Overall Pick" through "10th Overall Pick" (9 pick_outright markets)
  - "1st Cornerback" / "1st Offensive Line" / "1st Wide Receiver"
    (3 first_at_position markets)
  - "Mr Irrelevant Position" (prop - last-pick position; unmapped)

Notably Hoop88 does NOT post "1st Overall Pick" for this account tier, nor
any Top-N / Round-N range markets. The parser handles whatever propDescriptions
come back so new markets surface automatically once the fixture is refreshed.

Response shape (saved in tests/fixtures/hoop88/draft_markets.json):

    {
      "book": "hoop88",
      "meta": {"sport_type": "FOOTBALL", "sport_sub_type": "NFL Draft 2026",
               "period": "Prop", "markets": [...]},
      "data": {
        "markets": {
          "10th Overall Pick": {
            "Lines": [
              [
                "Odds to Win 10th Overall Pick|...|...|10th Overall Pick|...",
                [
                  {"ContestantName": "Caleb Downs       ",
                   "MoneyLine": 290, "ContestantNum": 235790619, ...},
                  ...
                ]
              ]
            ],
            "HourLimit": [], ...
          },
          "1st Cornerback": {...},
          ...
        }
      }
    }

Each Line is a 2-element list: [metadata_pipe_string, contestant_list].
The metadata field is noise for our purposes -- the propDescription key in
`data.markets` is the canonical market name we use as book_label.
"""

from __future__ import annotations

import os
import re
from datetime import datetime
from typing import List, Optional

from nfl_draft.scrapers._base import OddsRow


# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------

def _american(raw) -> Optional[int]:
    """Hoop88 ships MoneyLine as a signed int (e.g. -290, 2200). 0 means no
    line / suspended. Returns None for any unparseable / zero value."""
    if raw is None:
        return None
    try:
        iv = int(raw)
    except (TypeError, ValueError):
        return None
    if iv == 0:
        return None
    return iv


def _clean_name(name: str) -> str:
    """Hoop88 pads ContestantName with trailing spaces to a fixed width."""
    return (name or "").strip()


# Classify book_label (= propDescription sent to the API) into a market_group.
# Mirrors bookmaker.py/_classify_market shape so the MARKET_MAP builder in
# config/markets.py can re-key off the same strings.

# "10th Overall Pick" / "1st Overall Pick" / "2nd Overall Pick"
_PICK_LABEL_RE = re.compile(
    r"^(\d+)(?:st|nd|rd|th)\s+Overall\s+Pick$", re.IGNORECASE,
)
# "1st Wide Receiver" / "1st Cornerback" / "1st Offensive Line"
# (Hoop88 drops the 'Selected' suffix vs. BM's '1st Cornerback Selected')
_FIRST_POS_LABEL_RE = re.compile(
    r"^1st\s+(.+?)\s*$", re.IGNORECASE,
)

# Hoop88 position words -> canonical abbrev used in build_market_id.
POSITION_MAP = {
    "quarterback": "QB", "qb": "QB",
    "running back": "RB", "rb": "RB",
    "wide receiver": "WR", "wr": "WR",
    "tight end": "TE", "te": "TE",
    "cornerback": "CB", "cb": "CB",
    "safety": "S",
    "linebacker": "LB", "lb": "LB",
    "offensive line": "OL", "offensive lineman": "OL", "ol": "OL",
    "defensive back": "DB",
    "defensive tackle": "DT",
    "defensive end": "DE",
    "edge": "EDGE",
}


def _classify_label(label: str) -> str:
    """Return market_group for a Hoop88 propDescription / book_label.

    Returns:
      'pick_outright' for 'Nth Overall Pick'
      'first_at_position' for '1st <position>'
      'prop_<slug>' for anything else (e.g. 'Mr Irrelevant Position')
    """
    s = (label or "").strip()
    if _PICK_LABEL_RE.match(s):
        return "pick_outright"
    m = _FIRST_POS_LABEL_RE.match(s)
    if m:
        pos_raw = m.group(1).strip().lower()
        # Only return 'first_at_position' if the position is known --
        # otherwise the '1st <something>' string is actually a different
        # kind of prop we don't want to miscategorize.
        if pos_raw in POSITION_MAP:
            return "first_at_position"
    # Fallback: treat as generic prop, keyed by slugified label so quarantine
    # rows are grouped sensibly.
    slug = re.sub(r"[^a-z0-9]+", "_", s.lower()).strip("_")
    return f"prop_{slug or 'hoop88'}"


# ---------------------------------------------------------------------------
# Public parser
# ---------------------------------------------------------------------------

def parse_response(raw: dict) -> List[OddsRow]:
    """Walk Hoop88's per-market envelope and emit one OddsRow per contestant.

    Accepts either the full recon envelope (with 'book'/'meta'/'data') or the
    inner `data` dict directly. Silently drops markets with empty Lines,
    contestants with MoneyLine=0, and empty names (API quirk).
    """
    data = raw.get("data") if isinstance(raw, dict) and "data" in raw else raw
    if not isinstance(data, dict):
        return []

    markets = data.get("markets") or {}
    if not isinstance(markets, dict):
        return []

    now = datetime.now()
    rows: List[OddsRow] = []

    for market_name, resp in markets.items():
        if not isinstance(resp, dict):
            continue
        lines = resp.get("Lines") or []
        # Each entry is [metadata_str, contestant_list]. Some entries might
        # have just one element if the market was suspended.
        for entry in lines:
            if not isinstance(entry, list) or len(entry) < 2:
                continue
            contestants = entry[1]
            if not isinstance(contestants, list):
                continue
            group = _classify_label(market_name)
            for c in contestants:
                if not isinstance(c, dict):
                    continue
                name = _clean_name(c.get("ContestantName"))
                if not name:
                    continue
                american = _american(c.get("MoneyLine"))
                if american is None:
                    continue
                rows.append(OddsRow(
                    book="hoop88",
                    book_label=market_name,
                    book_subject=name,
                    american_odds=american,
                    fetched_at=now,
                    market_group=group,
                ))
    return rows


# ---------------------------------------------------------------------------
# Live fetch
# ---------------------------------------------------------------------------

def _login() -> tuple["requests.Session", str]:
    """Login to Hoop88 and return (authenticated Session, JWT token).

    The login flow mirrors hoop88_odds/scraper.py::login() exactly (same
    endpoint + form fields). We inline it here instead of importing from
    the production scraper because that module imports `dotenv` at
    top-level, which isn't installed in this venv. Keeping the flow in
    sync is a known maintenance cost -- if Hoop88 changes the auth path,
    both files need updating.
    """
    import requests

    # Ensure env is loaded so HOOP88_USERNAME / PASSWORD are visible.
    from nfl_draft.scrapers._recon_util import load_env
    load_env()

    hoop88_url = os.getenv("HOOP88_URL", "https://hoop88.com")
    username = os.getenv("HOOP88_USERNAME")
    password = os.getenv("HOOP88_PASSWORD")
    if not username or not password:
        raise RuntimeError("HOOP88_USERNAME / HOOP88_PASSWORD missing from bet_logger/.env")
    api_base = f"{hoop88_url}/cloud/api"

    sess = requests.Session()
    sess.headers.update({
        "X-Requested-With": "XMLHttpRequest",
        "Content-Type": "application/x-www-form-urlencoded; charset=UTF-8",
    })
    # Prime the Cloudflare cookie by hitting the site root.
    sess.get(hoop88_url, timeout=15)

    resp = sess.post(f"{api_base}/System/authenticateCustomer", data={
        "customerID": username,
        "password": password,
        "state": "true",
        "multiaccount": "1",
        "response_type": "code",
        "client_id": username,
        "domain": hoop88_url.replace("https://", "").replace("http://", ""),
        "redirect_uri": hoop88_url.replace("https://", "").replace("http://", ""),
        "operation": "authenticateCustomer",
        "RRO": "1",
    }, timeout=15)

    if resp.status_code == 204:
        raise RuntimeError("Hoop88 login failed -- bad credentials (HTTP 204)")
    resp.raise_for_status()

    token = resp.json().get("code")
    if not token:
        raise RuntimeError("Hoop88 login succeeded but no token in response")
    sess.headers["Authorization"] = f"Bearer {token}"
    return sess, token


def _enumerate_draft_markets(sess, customer_id: str, api_base: str) -> list[str]:
    """Query Get_SportsLeagues and return every SportSubType2 under NFL Draft 2026.

    Returning a live list (rather than hardcoding the 13 markets we saw at
    recon time) means the scraper auto-picks up new picks/positions as Hoop88
    posts them -- the 1st Overall Pick, Top N, Round 1, etc. markets will
    surface automatically if the book starts offering them without any code
    change here.
    """
    resp = sess.post(f"{api_base}/League/Get_SportsLeagues", data={
        "customerID": customer_id,
        "wagerType": "Straight",
        "office": "COINDEVIL",
        "placeLateFlag": "false",
        "operation": "Get_SportsLeagues",
        "RRO": "1",
        "agentSite": "0",
    }, timeout=20)
    resp.raise_for_status()
    leagues = resp.json().get("Leagues") or []
    markets: list[str] = []
    for e in leagues:
        if (e.get("SportType") == "FOOTBALL"
                and e.get("SportSubType") == "NFL Draft 2026"
                and e.get("SportSubType2")):
            markets.append(e["SportSubType2"])
    # dedupe preserving order
    return list(dict.fromkeys(markets))


def _fetch_market(sess, customer_id: str, api_base: str, prop_desc: str) -> dict:
    """POST Get_LeagueLines2 for one propDescription and return parsed JSON."""
    body = {
        "customerID": customer_id,
        "operation": "Get_LeagueLines2",
        "sportType": "FOOTBALL",
        "sportSubType": "NFL Draft 2026",
        "period": "Prop",
        "hourFilter": "0",
        "propDescription": prop_desc,
        "wagerType": "Straight",
        "keyword": "",
        "office": "COINDEVIL",
        "correlationID": "",
        "periodNumber": "0",
        "grouping": "",
        "periods": "0",
        "rotOrder": "0",
        "placeLateFlag": "false",
        "RRO": "1",
        "agentSite": "0",
    }
    resp = sess.post(f"{api_base}/Lines/Get_LeagueLines2", data=body, timeout=20)
    resp.raise_for_status()
    return resp.json()


def fetch_raw() -> dict:
    """Live fetch: login, enumerate draft props, fetch each, return envelope.

    Returns the same shape as the fixture: the parser walks
    `data.markets[propDescription]` -> Lines list.
    """
    from datetime import timezone

    sess, _token = _login()
    # After _login(), env is loaded so these are populated.
    customer_id = os.getenv("HOOP88_USERNAME")
    api_base = f"{os.getenv('HOOP88_URL', 'https://hoop88.com')}/cloud/api"

    market_names = _enumerate_draft_markets(sess, customer_id, api_base)

    market_responses: dict[str, dict] = {}
    for name in market_names:
        try:
            market_responses[name] = _fetch_market(sess, customer_id, api_base, name)
        except Exception as exc:
            # Skip individual market failures; the rest still come through.
            # Don't surface the exception message in case it contains creds.
            print(f"[hoop88] {name!r} fetch failed: {type(exc).__name__}")

    return {
        "book": "hoop88",
        "captured_at": datetime.now(timezone.utc).strftime("%Y-%m-%dT%H:%M:%SZ"),
        "meta": {
            "method": "rest",
            "sport_type": "FOOTBALL",
            "sport_sub_type": "NFL Draft 2026",
            "period": "Prop",
            "period_number": 0,
            "markets": market_names,
        },
        "data": {
            "markets": market_responses,
        },
    }


def fetch_draft_odds() -> List[OddsRow]:
    """Thin wrapper: live fetch -> parse."""
    return parse_response(fetch_raw())


# ---------------------------------------------------------------------------
# Exports for MARKET_MAP builder
# ---------------------------------------------------------------------------

PICK_LABEL_RE = _PICK_LABEL_RE
FIRST_POS_LABEL_RE = _FIRST_POS_LABEL_RE
