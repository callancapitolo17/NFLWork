"""Novig account auth for the bets service: a one-time connect that mints the
service's OWN Auth0 refresh token, and the refresh that turns it into 30-min
access tokens for api.novig.us.

    python -m unabated_ticket.bets_service.sources.novig_auth connect [--browser]

Why its own token: the web app keeps a rotating refresh token in localStorage;
Auth0 revokes a whole chain when a rotated token is reused, so copying the
app's token out of the browser would log the app out and kill the poller the
next time either side refreshed. A separate login = a separate chain.

Inputs:  the Novig web app's public Auth0 client (domain, client id, audience —
         read off app.novig.us's bundle, 2026-09-11; there is no secret) and
         one interactive login by the user (PKCE authorization-code flow).
Outputs: NOVIG_TOKEN_PATH (gitignored JSON: refresh_token, obtained_at,
         auth_id) — rewritten in place whenever Auth0 rotates the token.
Side effects: writes that file (mode 0600); no DuckDB, no other disk writes.
"""
import argparse
import base64
import hashlib
import json
import logging
import os
import secrets
import sys
import time
import urllib.error
import urllib.request
from dataclasses import dataclass
from pathlib import Path
from urllib.parse import parse_qs, urlencode, urlparse

from unabated_ticket.bets_service import config

log = logging.getLogger(__name__)

AUTH0_DOMAIN = "auth.novig.us"
AUTH0_CLIENT_ID = "gwTWbh0EewYND7LUCX0fKyv5xTzVpTCi"
AUTH0_AUDIENCE = "https://api.novig.us"
AUTH0_SCOPE = "openid profile email offline_access"
# The only callback the app registers (its redirect_uri is window.location.origin).
REDIRECT_URI = "https://app.novig.us"
TOKEN_URL = f"https://{AUTH0_DOMAIN}/oauth/token"
AUTHORIZE_URL = f"https://{AUTH0_DOMAIN}/authorize"
HTTP_TIMEOUT_SEC = 20
# Refresh when this close to expiry, so a poll never starts on a dying token.
REFRESH_MARGIN_SEC = 120


class NovigAuthError(RuntimeError):
    """Token exchange or refresh failed; the message says which and why."""


@dataclass
class AccessToken:
    token: str
    expires_at: float
    auth_id: str


def jwt_subject(access_token: str) -> str:
    """The `sub` claim (Novig's user.auth_id). Decoded without verification —
    Hasura verifies the token; this only names the account for the where-clause."""
    try:
        payload = access_token.split(".")[1]
        payload += "=" * (-len(payload) % 4)
        claims = json.loads(base64.urlsafe_b64decode(payload))
        return str(claims["sub"])
    except (IndexError, KeyError, ValueError) as error:
        raise NovigAuthError(f"access token has no readable sub claim: {error}") from error


def _token_request(body: dict) -> dict:
    """POST to Auth0's token endpoint with the standard library — `connect` must
    run from any venv (mlb_sgp's has Playwright but no requests)."""
    request = urllib.request.Request(TOKEN_URL, data=json.dumps(body).encode(),
                                     headers={"Content-Type": "application/json", "Accept": "application/json"})
    try:
        with urllib.request.urlopen(request, timeout=HTTP_TIMEOUT_SEC) as response:
            status, text = response.status, response.read().decode()
    except urllib.error.HTTPError as error:
        status, text = error.code, error.read().decode(errors="replace")
    except urllib.error.URLError as error:
        raise NovigAuthError(f"Auth0 {body.get('grant_type')} unreachable: {error.reason}") from error
    try:
        payload = json.loads(text)
    except ValueError:
        payload = {}
    if status != 200 or "access_token" not in payload:
        detail = payload.get("error_description") or payload.get("error") or text[:200]
        raise NovigAuthError(f"Auth0 {body.get('grant_type')} failed ({status}): {detail}")
    return payload


def load_token_file(path: Path) -> dict | None:
    if not path.exists():
        return None
    try:
        stored = json.loads(path.read_text())
    except ValueError as error:
        raise NovigAuthError(f"{path} is not valid JSON: {error}") from error
    if not isinstance(stored, dict) or not stored.get("refresh_token"):
        raise NovigAuthError(f"{path} has no refresh_token; run `novig_auth connect`")
    return stored


def save_token_file(path: Path, refresh_token: str, auth_id: str) -> None:
    """Atomic rewrite, owner-only permissions."""
    tmp = path.with_suffix(path.suffix + ".tmp")
    tmp.write_text(json.dumps({"refresh_token": refresh_token, "auth_id": auth_id,
                               "obtained_at": time.strftime("%Y-%m-%dT%H:%M:%SZ", time.gmtime())}, indent=1))
    os.chmod(tmp, 0o600)
    os.replace(tmp, path)


class NovigAuth:
    """Holds the refresh token and hands out a live access token.

    `token()` refreshes when the cached access token is within
    REFRESH_MARGIN_SEC of expiry; a rotated refresh token is persisted at once
    (losing it would orphan the chain). `request` is injectable for tests.
    """

    def __init__(self, token_path: Path, request=_token_request, clock=time.time):
        self._path = token_path
        self._request = request
        self._clock = clock
        self._access: AccessToken | None = None

    def token(self) -> AccessToken:
        if self._access and self._access.expires_at - self._clock() > REFRESH_MARGIN_SEC:
            return self._access
        stored = load_token_file(self._path)
        if stored is None:
            raise NovigAuthError(f"no Novig token at {self._path}; run "
                                 "`python -m unabated_ticket.bets_service.sources.novig_auth connect`")
        payload = self._request({"grant_type": "refresh_token", "client_id": AUTH0_CLIENT_ID,
                                 "refresh_token": stored["refresh_token"]})
        auth_id = jwt_subject(payload["access_token"])
        rotated = payload.get("refresh_token")
        if rotated and rotated != stored["refresh_token"]:
            save_token_file(self._path, rotated, auth_id)
            log.info("novig: refresh token rotated and saved")
        self._access = AccessToken(payload["access_token"], self._clock() + float(payload.get("expires_in", 1800)), auth_id)
        return self._access


# ---- connect (one-time, interactive) --------------------------------------------

def pkce_pair() -> tuple[str, str]:
    verifier = base64.urlsafe_b64encode(secrets.token_bytes(48)).rstrip(b"=").decode()
    challenge = base64.urlsafe_b64encode(hashlib.sha256(verifier.encode()).digest()).rstrip(b"=").decode()
    return verifier, challenge


def authorize_url(challenge: str, state: str) -> str:
    return AUTHORIZE_URL + "?" + urlencode({
        "response_type": "code", "client_id": AUTH0_CLIENT_ID, "redirect_uri": REDIRECT_URI,
        "scope": AUTH0_SCOPE, "audience": AUTH0_AUDIENCE, "code_challenge": challenge,
        "code_challenge_method": "S256", "state": state, "prompt": "login",
    })


def code_from_redirect(redirected_url: str, expected_state: str) -> str:
    query = parse_qs(urlparse(redirected_url.strip()).query)
    if query.get("error"):
        raise NovigAuthError(f"login refused: {query['error'][0]}: {query.get('error_description', [''])[0]}")
    if query.get("state", [None])[0] != expected_state:
        raise NovigAuthError("state mismatch — paste the URL from THIS run's login, not an older one")
    code = query.get("code", [None])[0]
    if not code:
        raise NovigAuthError("no ?code= in that URL (the app may have stripped it — try --browser)")
    return code


def _redirect_via_paste(url: str) -> str:
    print("\n1. Open this URL in a browser and log in to Novig:\n\n" + url + "\n")
    print("2. After login you land on https://app.novig.us/?code=...&state=... — copy that")
    print("   address IMMEDIATELY (the app removes it once it loads).\n")
    return input("Paste the redirected URL here: ")


def _redirect_via_browser(url: str) -> str:
    """Headed Playwright window; the app never loads because the callback
    request is intercepted, so the code cannot be stripped. The user types
    their own credentials; nothing here reads them."""
    try:
        from playwright.sync_api import sync_playwright
    except ImportError as error:  # pragma: no cover
        raise NovigAuthError("--browser needs Playwright (mlb_sgp/venv has it); use the paste flow instead") from error
    landed: list[str] = []
    with sync_playwright() as p:
        browser = p.chromium.launch(headless=False)
        page = browser.new_page()

        def intercept(route, request):
            landed.append(request.url)
            route.fulfill(status=200, content_type="text/html",
                          body="<h2>Novig connected — you can close this window.</h2>")

        page.route(f"{REDIRECT_URI}/**", intercept)
        page.route(REDIRECT_URI, intercept)
        page.goto(url)
        deadline = time.time() + 300
        while not landed and time.time() < deadline:
            page.wait_for_timeout(250)
        browser.close()
    if not landed:
        raise NovigAuthError("login did not complete within 5 minutes")
    return landed[0]


def connect(token_path: Path, use_browser: bool) -> None:
    verifier, challenge = pkce_pair()
    state = secrets.token_urlsafe(16)
    url = authorize_url(challenge, state)
    redirected = _redirect_via_browser(url) if use_browser else _redirect_via_paste(url)
    code = code_from_redirect(redirected, state)
    payload = _token_request({"grant_type": "authorization_code", "client_id": AUTH0_CLIENT_ID,
                              "code_verifier": verifier, "code": code, "redirect_uri": REDIRECT_URI})
    if not payload.get("refresh_token"):
        raise NovigAuthError("Auth0 returned no refresh_token (offline_access not granted)")
    auth_id = jwt_subject(payload["access_token"])
    save_token_file(token_path, payload["refresh_token"], auth_id)
    print(f"\nSaved {token_path} for account {auth_id}. Start the bets service; the novig source is now registered.")


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description="Connect the bets service to a Novig account (one-time).")
    sub = parser.add_subparsers(dest="command", required=True)
    connect_parser = sub.add_parser("connect", help="log in once and save a refresh token")
    connect_parser.add_argument("--browser", action="store_true",
                                help="open a Playwright window that captures the callback (no copy/paste race)")
    args = parser.parse_args(argv)
    try:
        connect(config.NOVIG_TOKEN_PATH, args.browser)
    except NovigAuthError as error:
        print(f"connect failed: {error}", file=sys.stderr)
        return 1
    return 0


if __name__ == "__main__":
    sys.exit(main())
