#!/usr/bin/env python3
"""Capture raw account data from BFA and Wagerzon (the C account) so their bet
sources for the bets service can be built on real rows. The cloud sessions
that write the parsers cannot reach either book, so this runs on the Mac and
the two output files are attached to the project thread.

    cd /Users/callancapitolo/NFLWork
    bet_logger/venv/bin/python3 unabated_ticket/bets_service/recon_capture.py [--out DIR]

Standalone on purpose (imports only `requests`): a downloaded copy runs too,
as long as the working directory is inside the NFLWork checkout.

Inputs:  credentials from the environment, else bet_logger/.env found by
         walking up from the working directory (then from this file's own
         directory):
           BFA_USERNAME / BFA_PASSWORD
           WAGERZONC_USERNAME / WAGERZONC_PASSWORD, else WAGERZON_USERNAME /
           WAGERZON_PASSWORD (the C account's login has lived in the primary
           slot since 2026-06-26 — bet_logger/scraper_wagerzon.py).
Outputs: <out>/bfa_capture.json      every GetPlayerHistory wager of the last
                                     CAPTURE_DAYS (through tomorrow, so a bet
                                     placed today is inside the window), raw,
                                     plus a count per result and per type.
         <out>/wagerzon_capture.json HistoryHelper weeks 0 and 1 raw, the
                                     OpenBets.aspx page, and every open /
                                     pending helper endpoint the page or its
                                     scripts name, fetched and stored raw.
         DIR defaults to ~/Downloads/bets_recon. A venue that fails is recorded
         under "error" in its own file and the other venue still runs.
Side effects: two logins. BFA mints a fresh Keycloak session (bet_logger's
         recon_bfa_auth.json is neither read nor written). No token, cookie or
         password is written anywhere.
"""
import argparse
import base64
import hashlib
import json
import os
import re
import secrets
import sys
from datetime import datetime, timedelta, timezone
from pathlib import Path
from urllib.parse import parse_qs, urljoin, urlparse

import requests

CAPTURE_DAYS = 30
HTTP_TIMEOUT_SEC = 20
DEFAULT_OUT_DIR = Path.home() / "Downloads" / "bets_recon"
ENV_FILE_RELATIVE = Path("bet_logger") / ".env"
# Page bodies are stored for reading, not replaying; this keeps a file small.
MAX_STORED_TEXT_BYTES = 300_000
USER_AGENT = ("Mozilla/5.0 (Macintosh; Intel Mac OS X 10_15_7) "
              "AppleWebKit/537.36 (KHTML, like Gecko) Chrome/145.0.0.0 Safari/537.36")

# BFA (bet_logger/recon_bfa.py + scraper_bfa.py, kept in step by hand).
BFA_KEYCLOAK_BASE = "https://auth.bfagaming.com/realms/players_realm/protocol/openid-connect"
BFA_CLIENT_ID = "bfagaming"
BFA_REDIRECT_URI = "https://bfagaming.com/"
BFA_HISTORY_URL = "https://api.bfagaming.com/history/api/GetPlayerHistory"
BFA_RECORDS_PER_PAGE = 100
BFA_MAX_PAGES = 50
BFA_HEADERS = {"Accept": "application/json", "Origin": "https://bfagaming.com",
               "Referer": "https://bfagaming.com/", "User-Agent": USER_AGENT}

# Wagerzon (bet_logger/scraper_wagerzon.py; OpenBets.aspx from Cal, 2026-09-22).
WAGERZON_BASE_URL = "https://backend.wagerzon.com"
WAGERZON_HISTORY_URL = f"{WAGERZON_BASE_URL}/wager/HistoryHelper.aspx"
WAGERZON_OPEN_BETS_URL = f"{WAGERZON_BASE_URL}/wager/OpenBets.aspx"
# The helper the open-bets widget most likely calls, by analogy with
# History.aspx -> HistoryHelper.aspx; tried first, then whatever the page names.
WAGERZON_OPEN_BETS_HELPER_GUESS = f"{WAGERZON_BASE_URL}/wager/OpenBetsHelper.aspx"
WAGERZON_HISTORY_WEEKS = (0, 1)
WAGERZON_XHR_HEADERS = {"X-Requested-With": "XMLHttpRequest", "Accept": "application/json"}
WAGERZON_ASPNET_HIDDEN_FIELDS = ("__VIEWSTATE", "__VIEWSTATEGENERATOR", "__EVENTVALIDATION",
                                 "__EVENTTARGET", "__EVENTARGUMENT")
OPEN_BETS_HELPER_NAME_RE = re.compile(r"open|pending", re.IGNORECASE)
ASPX_REFERENCE_RE = re.compile(r"[\w./-]*\w\.aspx(?:\?[^\s\"'<>)]*)?")
SCRIPT_SRC_RE = re.compile(r"<script[^>]+src=[\"']([^\"']+)[\"']", re.IGNORECASE)
MAX_HELPER_FETCHES = 15


# ---- credentials --------------------------------------------------------------------

def read_env_file(path: Path) -> dict[str, str]:
    values: dict[str, str] = {}
    for raw_line in path.read_text().splitlines():
        line = raw_line.strip()
        if not line or line.startswith("#") or "=" not in line:
            continue
        key, value = line.split("=", 1)
        values[key.strip()] = value.strip().strip('"').strip("'")
    return values


def find_env_file() -> Path | None:
    """bet_logger/.env, walking up from the working directory, then from here."""
    for start in (Path.cwd(), Path(__file__).resolve().parent):
        for directory in (start, *start.parents):
            candidate = directory / ENV_FILE_RELATIVE
            if candidate.exists():
                return candidate
    return None


def load_credentials() -> dict[str, str]:
    env_file = find_env_file()
    file_values = read_env_file(env_file) if env_file else {}
    merged = {**file_values, **os.environ}
    print(f"credentials: {env_file or 'environment only (no bet_logger/.env found)'}")
    return merged


def required(values: dict[str, str], *keys: str) -> tuple[str, str]:
    """(username, password) from the first key pair that is set."""
    for username_key, password_key in zip(keys[::2], keys[1::2]):
        username, password = values.get(username_key), values.get(password_key)
        if username and password:
            return username, password
    raise RuntimeError(f"expected {' or '.join(keys[::2])} (+ password) in the environment "
                       "or bet_logger/.env, found none set")


# ---- BFA ------------------------------------------------------------------------------

def pkce_pair() -> tuple[str, str]:
    verifier = secrets.token_urlsafe(32)
    challenge = base64.urlsafe_b64encode(hashlib.sha256(verifier.encode()).digest()).rstrip(b"=").decode()
    return verifier, challenge


def jwt_payload(token: str) -> dict:
    parts = token.split(".")
    if len(parts) != 3:
        return {}
    padded = parts[1] + "=" * (-len(parts[1]) % 4)
    return json.loads(base64.urlsafe_b64decode(padded))


def bfa_login(username: str, password: str) -> tuple[str, str]:
    """Keycloak password login (PKCE) -> (access_token, player_id)."""
    session = requests.Session()
    session.headers["User-Agent"] = USER_AGENT
    verifier, challenge = pkce_pair()
    login_page = session.get(
        f"{BFA_KEYCLOAK_BASE}/auth",
        params={"client_id": BFA_CLIENT_ID, "redirect_uri": BFA_REDIRECT_URI, "response_type": "code",
                "scope": "openid", "code_challenge": challenge, "code_challenge_method": "S256"},
        timeout=HTTP_TIMEOUT_SEC)
    if login_page.status_code != 200:
        raise RuntimeError(f"BFA login page: HTTP {login_page.status_code}")
    action = re.search(r'action="([^"]+)"', login_page.text)
    if not action:
        raise RuntimeError("BFA login page carries no form action")
    submitted = session.post(action.group(1).replace("&amp;", "&"),
                             data={"username": username, "password": password},
                             allow_redirects=False, timeout=HTTP_TIMEOUT_SEC)
    if submitted.status_code not in (302, 303):
        detail = "invalid credentials" if "Invalid username or password" in submitted.text else "no redirect"
        raise RuntimeError(f"BFA login failed: HTTP {submitted.status_code} ({detail})")
    code = parse_qs(urlparse(submitted.headers.get("Location", "")).query).get("code", [None])[0]
    if not code:
        raise RuntimeError("BFA login redirect carries no auth code")
    token_form = {"grant_type": "authorization_code", "client_id": BFA_CLIENT_ID, "code": code,
                  "redirect_uri": BFA_REDIRECT_URI, "code_verifier": verifier}
    exchanged = session.post(f"{BFA_KEYCLOAK_BASE}/token", data=token_form, timeout=HTTP_TIMEOUT_SEC)
    if exchanged.status_code != 200:
        raise RuntimeError(f"BFA token exchange failed: HTTP {exchanged.status_code} {exchanged.text[:200]}")
    access_token = exchanged.json()["access_token"]
    player_id = jwt_payload(access_token).get("player_id")
    if not player_id:
        raise RuntimeError("BFA access token carries no player_id")
    return access_token, str(player_id)


def bfa_history_pages(access_token: str, player_id: str, start_date: str, end_date: str) -> list[dict]:
    """Every page body of GetPlayerHistory over the window, raw."""
    headers = {"Authorization": f"Bearer {access_token}", **BFA_HEADERS}
    pages: list[dict] = []
    wagers_seen = 0
    for page in range(BFA_MAX_PAGES):
        response = requests.get(BFA_HISTORY_URL, headers=headers, timeout=HTTP_TIMEOUT_SEC,
                                params={"playerId": player_id, "startDate": start_date, "endDate": end_date,
                                        "page": page, "recordsByPage": BFA_RECORDS_PER_PAGE})
        if response.status_code != 200:
            raise RuntimeError(f"BFA history page {page}: HTTP {response.status_code} {response.text[:200]}")
        body = response.json()
        pages.append(body)
        page_wagers = body.get("wagers") or []
        wagers_seen += len(page_wagers)
        if not page_wagers or wagers_seen >= int(body.get("totalRecords") or 0):
            return pages
    raise RuntimeError(f"BFA history did not end within {BFA_MAX_PAGES} pages")


def count_by(rows: list[dict], key: str) -> dict[str, int]:
    counts: dict[str, int] = {}
    for row in rows:
        value = str(row.get(key))
        counts[value] = counts.get(value, 0) + 1
    return dict(sorted(counts.items()))


def capture_bfa(credentials: dict[str, str]) -> dict:
    username, password = required(credentials, "BFA_USERNAME", "BFA_PASSWORD")
    access_token, player_id = bfa_login(username, password)
    today = datetime.now()
    start_date = (today - timedelta(days=CAPTURE_DAYS)).strftime("%Y-%m-%d")
    end_date = (today + timedelta(days=1)).strftime("%Y-%m-%d")
    pages = bfa_history_pages(access_token, player_id, start_date, end_date)
    wagers = [wager for page in pages for wager in page.get("wagers") or []]
    return {
        "window": {"startDate": start_date, "endDate": end_date},
        "pageBodyKeys": sorted(pages[0].keys()) if pages else [],
        "totalRecords": pages[0].get("totalRecords") if pages else None,
        "wagerCount": len(wagers),
        "countByResult": count_by(wagers, "result"),
        "countByType": count_by(wagers, "type"),
        "wagers": wagers,
    }


# ---- Wagerzon -------------------------------------------------------------------------

def wagerzon_login(session: requests.Session, username: str, password: str) -> None:
    """ASP.NET form login: the hidden fields are posted back exactly as served."""
    login_page = session.get(WAGERZON_BASE_URL, timeout=HTTP_TIMEOUT_SEC)
    login_page.raise_for_status()
    fields = {}
    for name in WAGERZON_ASPNET_HIDDEN_FIELDS:
        match = re.search(rf'(?:name|id)="{name}"[^>]*value="([^"]*)"', login_page.text)
        if match:
            fields[name] = match.group(1)
    if "__VIEWSTATE" not in fields:
        raise RuntimeError("Wagerzon login page carries no __VIEWSTATE")
    fields.update({"Account": username, "Password": password, "BtnSubmit": ""})
    session.post(WAGERZON_BASE_URL, data=fields, timeout=HTTP_TIMEOUT_SEC).raise_for_status()


def stored_body(response: requests.Response) -> dict:
    """A response as JSON when it parses, else its text (truncated), for the file."""
    content_type = response.headers.get("Content-Type", "")
    entry = {"url": response.url, "status": response.status_code, "contentType": content_type}
    try:
        entry["json"] = response.json()
    except ValueError:
        text = response.text
        entry["text"] = text[:MAX_STORED_TEXT_BYTES]
        entry["textTruncated"] = len(text) > MAX_STORED_TEXT_BYTES
    return entry


def wagerzon_history(session: requests.Session, week: int) -> dict:
    response = session.get(WAGERZON_HISTORY_URL, params={"week": week}, headers=WAGERZON_XHR_HEADERS,
                           timeout=HTTP_TIMEOUT_SEC)
    entry = stored_body(response)
    if "json" not in entry:
        raise RuntimeError(f"Wagerzon HistoryHelper week {week} returned no JSON (HTTP {response.status_code}, "
                           f"{response.url}) — the login did not take")
    return entry


def same_host(url: str) -> bool:
    return urlparse(url).netloc == urlparse(WAGERZON_BASE_URL).netloc


def aspx_references(text: str, base_url: str) -> list[str]:
    return sorted({urljoin(base_url, match) for match in ASPX_REFERENCE_RE.findall(text)})


def open_bets_helper_candidates(session: requests.Session, page: requests.Response) -> tuple[list[str], dict]:
    """Every .aspx the OpenBets page or its same-host scripts name whose name reads
    open/pending — the endpoints its widget could be loading rows from."""
    references: dict[str, list[str]] = {page.url: aspx_references(page.text, page.url)}
    for script_src in SCRIPT_SRC_RE.findall(page.text):
        script_url = urljoin(page.url, script_src)
        if not same_host(script_url):
            continue
        script = session.get(script_url, timeout=HTTP_TIMEOUT_SEC)
        references[script_url] = aspx_references(script.text, script_url) if script.status_code == 200 else []
    candidates = [WAGERZON_OPEN_BETS_HELPER_GUESS]
    for found in references.values():
        for url in found:
            name = Path(urlparse(url).path).name
            if OPEN_BETS_HELPER_NAME_RE.search(name) and url not in candidates and url != page.url:
                candidates.append(url)
    return candidates[:MAX_HELPER_FETCHES], references


def capture_wagerzon(credentials: dict[str, str]) -> dict:
    username, password = required(credentials, "WAGERZONC_USERNAME", "WAGERZONC_PASSWORD",
                                  "WAGERZON_USERNAME", "WAGERZON_PASSWORD")
    session = requests.Session()
    session.headers["User-Agent"] = USER_AGENT
    wagerzon_login(session, username, password)
    history = {f"week{week}": wagerzon_history(session, week) for week in WAGERZON_HISTORY_WEEKS}
    open_bets_page = session.get(WAGERZON_OPEN_BETS_URL, timeout=HTTP_TIMEOUT_SEC)
    candidates, references = open_bets_helper_candidates(session, open_bets_page)
    helpers = [stored_body(session.get(url, headers=WAGERZON_XHR_HEADERS, timeout=HTTP_TIMEOUT_SEC))
               for url in candidates]
    return {
        "history": history,
        "openBetsPage": stored_body(open_bets_page),
        "aspxReferencesByFile": references,
        "openBetsHelpers": helpers,
        "helpersReturningJson": [helper["url"] for helper in helpers if "json" in helper],
    }


# ---- entry point ----------------------------------------------------------------------

def run_capture(name: str, capture, credentials: dict[str, str], out_dir: Path) -> None:
    print(f"\n== {name}")
    result = {"capturedAt": datetime.now(timezone.utc).strftime("%Y-%m-%dT%H:%M:%SZ"), "venue": name}
    try:
        result.update(capture(credentials))
    except Exception as error:  # noqa: BLE001 — recorded in the file; the other venue still runs
        result["error"] = f"{type(error).__name__}: {error}"
        print(f"FAILED: {result['error']}")
    path = out_dir / f"{name}_capture.json"
    path.write_text(json.dumps(result, indent=1, default=str))
    print(f"wrote {path} ({path.stat().st_size:,} bytes)")


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    parser.add_argument("--out", type=Path, default=DEFAULT_OUT_DIR, help=f"output directory (default {DEFAULT_OUT_DIR})")
    args = parser.parse_args()
    args.out.mkdir(parents=True, exist_ok=True)
    credentials = load_credentials()
    run_capture("bfa", capture_bfa, credentials, args.out)
    run_capture("wagerzon", capture_wagerzon, credentials, args.out)
    print(f"\nAttach both files in {args.out} to the project thread.")
    return 0


if __name__ == "__main__":
    sys.exit(main())
