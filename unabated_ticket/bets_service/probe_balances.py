"""One-shot, read-only probe: where does each venue report the account balance?

    cd ~/NFLWork && kalshi_draft/venv/bin/python3 /tmp/probe_balances.py

Throwaway discovery for the tracker's "balance by book" panel; removed once the
fields are confirmed. Uses each bets-service source's own login (config.py).

Outputs: per venue, each tried GET with its HTTP status, and for JSON answers only
         the fields whose NAME looks like money (balance, available, cash, risk,
         credit, funds...), with numbers rounded to whole units. Tokens, passwords,
         ids and any text values are never printed.
Side effects: GETs only, plus each source's normal login. BetOnline/Novig may
         rotate their refresh token file exactly as a service poll does. No DuckDB.
"""
import os
import re
import sys
from datetime import datetime, timezone
from urllib.parse import urlsplit

sys.path.insert(0, os.getcwd())

MONEY_KEY = re.compile(r"bal|avail|cash|risk|credit|fund|wallet|pending|exposure|equity|buying|withdraw|amount|value",
                       re.I)
SECRET_KEY = re.compile(r"pass|token|secret|key|email|phone|name|address|ssn|birth|^id$", re.I)
PATH_IN_BUNDLE = re.compile(r"""["'`]((?:https?://[a-z0-9.-]+)?/[A-Za-z0-9_/.{}$-]*"""
                            r"""(?:[Bb]alance|[Ww]allet|[Ff]unds|[Bb]ankroll)[A-Za-z0-9_/.-]*)["'`]""")
SCRIPT_SRC = re.compile(r"""<script[^>]+src=["']([^"']+)["']""", re.I)
MAX_BUNDLES = 25
MAX_DISCOVERED = 12
TIMEOUT_SEC = 20


def http_session():
    try:
        from curl_cffi import requests as cffi_requests
        return cffi_requests.Session(impersonate="chrome")
    except ImportError:
        import requests
        return requests.Session()


def as_number(value):
    if isinstance(value, bool):
        return None
    if isinstance(value, (int, float)):
        return value
    if isinstance(value, str):
        cleaned = value.replace(",", "").replace("$", "").strip()
        try:
            return float(cleaned)
        except ValueError:
            return None
    return None


def money_fields(body, prefix=""):
    """(dotted key, rounded number or type name) for money-named keys only."""
    found = []
    if isinstance(body, dict):
        for key, value in body.items():
            path = f"{prefix}.{key}" if prefix else str(key)
            if isinstance(value, (dict, list)):
                found.extend(money_fields(value, path))
            elif MONEY_KEY.search(str(key)) and not SECRET_KEY.search(str(key)):
                number = as_number(value)
                shown = round(number) if number is not None else f"<{type(value).__name__}>"
                found.append((path, shown, type(value).__name__))
    elif isinstance(body, list):
        for index, item in enumerate(body[:3]):
            found.extend(money_fields(item, f"{prefix}[{index}]"))
    return found


def report(venue, label, status, body):
    print(f"  {venue} GET {label} -> HTTP {status}")
    if body is None:
        return
    fields = money_fields(body)
    if not fields:
        top = list(body.keys())[:15] if isinstance(body, dict) else type(body).__name__
        print(f"      no money-named fields; top-level keys: {top}")
    for path, shown, kind in fields:
        print(f"      {path} = {shown}   ({kind})")


def json_or_none(response):
    try:
        return response.json()
    except Exception:  # noqa: BLE001 — a non-JSON answer is reported by status alone
        return None


def discover_paths(site_url):
    """Public site JS (no login): quoted URL paths mentioning balance/wallet/funds."""
    session = http_session()
    try:
        home = session.get(site_url, timeout=TIMEOUT_SEC)
    except Exception as exc:  # noqa: BLE001
        print(f"  bundle scan {site_url}: {type(exc).__name__}")
        return []
    sources = SCRIPT_SRC.findall(home.text)[:MAX_BUNDLES]
    found = set(PATH_IN_BUNDLE.findall(home.text))
    for src in sources:
        url = src if src.startswith("http") else site_url.rstrip("/") + "/" + src.lstrip("/")
        try:
            found.update(PATH_IN_BUNDLE.findall(session.get(url, timeout=TIMEOUT_SEC).text))
        except Exception:  # noqa: BLE001
            continue
    paths = sorted(path for path in found if "{" not in path and "$" not in path)
    print(f"  bundle scan {site_url}: HTTP {home.status_code}, {len(sources)} scripts, "
          f"{len(paths)} candidate paths: {paths[:30]}")
    return paths[:MAX_DISCOVERED]


def full_urls(paths, api_bases):
    urls = []
    for path in paths:
        if path.startswith("http"):
            if urlsplit(path).netloc in {urlsplit(base).netloc for base in api_bases}:
                urls.append(path)
        else:
            urls.extend(base.rstrip("/") + path for base in api_bases)
    return urls


def probe(venue, fn):
    print(f"\n== {venue}")
    try:
        fn()
    except Exception as exc:  # noqa: BLE001 — one venue failing never stops the rest
        print(f"  FAILED: {type(exc).__name__}: {str(exc)[:200]}")


def probe_kalshi():
    from kalshi_common import auth_client
    from unabated_ticket.bets_service.sources.kalshi import KalshiSource
    KalshiSource()
    status, body, _headers = auth_client.api("GET", "/portfolio/balance")
    report("kalshi", "/portfolio/balance", status, body if isinstance(body, dict) else None)


def probe_wagerzon():
    from unabated_ticket.bets_service.sources.wagerzon import BASE_URL, WagerzonSource
    source = WagerzonSource()
    body = source._get_json(f"{BASE_URL}/wager/PlayerInfoHelper.aspx", None)
    report("wagerzon", "/wager/PlayerInfoHelper.aspx", 200, body)


def probe_bfa():
    from unabated_ticket.bets_service.sources.bfa import BFASource
    headers, player_id = BFASource().authorized()
    session = http_session()
    guesses = ["/player/api/GetPlayerBalance", "/player/api/GetPlayerInfo", "/history/api/GetPlayerBalance",
               "/account/api/GetPlayerBalance", "/wagering/api/v1/balance"]
    paths = guesses + discover_paths("https://bfagaming.com/")
    for url in dict.fromkeys(full_urls(paths, ["https://api.bfagaming.com"])):
        response = session.get(url, headers=headers, params={"playerId": player_id}, timeout=TIMEOUT_SEC)
        report("bfa", url, response.status_code, json_or_none(response))


def probe_betonline():
    from unabated_ticket.bets_service.sources.betonline import BetOnlineSource, api_headers
    access_token, session = BetOnlineSource()._authenticated_session()
    headers = api_headers(access_token, datetime.now(timezone.utc))
    guesses = ["/customer/api/customer/balance", "/report/api/report/get-balance", "/customer/api/balance",
               "/wallet/api/balance"]
    paths = guesses + discover_paths("https://www.betonline.ag/")
    for url in dict.fromkeys(full_urls(paths, ["https://api.betonline.ag"])):
        response = session.get(url, headers=headers, timeout=TIMEOUT_SEC)
        report("betonline", url, response.status_code, json_or_none(response))


def probe_novig():
    from unabated_ticket.bets_service import config
    from unabated_ticket.bets_service.sources.novig import CLIENT_CAPABILITIES_HEADER
    from unabated_ticket.bets_service.sources.novig_auth import NovigAuth
    bearer = NovigAuth(config.NOVIG_TOKEN_PATH).token().token
    headers = {"Authorization": f"Bearer {bearer}", "Accept": "application/json",
               CLIENT_CAPABILITIES_HEADER[0]: CLIENT_CAPABILITIES_HEADER[1]}
    base = config.NOVIG_REST_URL.rstrip("/")
    guesses = ["/portfolio/balance", "/portfolio/balances", "/portfolio/summary", "/wallet", "/wallet/balance",
               "/user/balance", "/account/balance"]
    paths = guesses + discover_paths("https://novig.com/")
    session = http_session()
    for url in dict.fromkeys(full_urls(paths, [base, "https://api.novig.us"])):
        response = session.get(url, headers=headers, timeout=TIMEOUT_SEC)
        report("novig", url, response.status_code, json_or_none(response))


def probe_polymarket_us():
    from unabated_ticket.bets_service.sources.polymarket_us import API_BASE_URL, PolymarketUSSource
    source = PolymarketUSSource()
    guesses = ["/v1/account/balances", "/v1/account/balance", "/v1/portfolio/balance", "/v1/portfolio/balances",
               "/v1/portfolio/summary", "/v1/account"]
    paths = guesses + [path for path in discover_paths("https://polymarket.us/") if path.startswith("/v1/")]
    for path in dict.fromkeys(paths):
        response = source._session.get(f"{API_BASE_URL}{path}", headers=source._signed_headers(path),
                                       timeout=TIMEOUT_SEC)
        report("polymarket_us", path, response.status_code, json_or_none(response))


if __name__ == "__main__":
    for venue, fn in [("kalshi", probe_kalshi), ("wagerzon", probe_wagerzon), ("bfa", probe_bfa),
                      ("betonline", probe_betonline), ("novig", probe_novig),
                      ("polymarket_us", probe_polymarket_us)]:
        probe(venue, fn)
    print("\nDone. Bet105 and DraftKings are read from Chrome; not probed here.")
