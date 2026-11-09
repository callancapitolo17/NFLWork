"""Underdog account auth for the bets service: a one-time connect that mints the
service's OWN Auth0 refresh token, and the refresh that turns it into 10-min
access tokens for api.underdogfantasy.com.

    python -m unabated_ticket.bets_service.sources.underdog_auth connect

Why its own token (the novig_auth.py rule): the web app keeps a rotating
refresh token in its `session_refresh` cookie; Auth0 revokes a whole chain when
a rotated token is reused, so copying the app's token would log the app out
and kill the poller the next time either side refreshed. A separate login is a
separate chain.

Why a password grant: Underdog's own web client logs in with Auth0's password
grant (its access tokens carry gty ["refresh_token", "password"], read
2026-11-07), so the same client takes `grant_type=password` here. The password
is typed at a getpass prompt, sent once to Auth0 and never stored or logged;
only the refresh token is kept.

Inputs:  the web app's public Auth0 client (domain, client id, audience, read
         off the app's own access token 2026-11-07; there is no secret) and
         the account's email + password, typed once by the user.
Outputs: UNDERDOG_TOKEN_PATH (gitignored JSON: refresh_token, obtained_at,
         user_id) — rewritten in place whenever Auth0 rotates the token.
Side effects: writes that file (mode 0600); no DuckDB, no other disk writes.
"""
import argparse
import base64
import getpass
import json
import logging
import os
import sys
import time
import urllib.error
import urllib.request
from dataclasses import dataclass
from pathlib import Path

from unabated_ticket.bets_service import config

log = logging.getLogger(__name__)

AUTH0_DOMAIN = "login.underdogsports.com"
AUTH0_CLIENT_ID = "cQvYz1T2BAFbix4dYR37dyD9O0Thf1s6"
AUTH0_AUDIENCE = "https://api.underdogfantasy.com"
AUTH0_SCOPE = "offline_access"
TOKEN_URL = f"https://{AUTH0_DOMAIN}/oauth/token"
HTTP_TIMEOUT_SEC = 20
# The app's tokens live 10 min; refresh when this close to expiry so a poll never starts on a dying one.
REFRESH_MARGIN_SEC = 60
DEFAULT_TOKEN_LIFE_SEC = 600


class UnderdogAuthError(RuntimeError):
    """Login or refresh failed; the message says which and why."""


@dataclass
class AccessToken:
    token: str
    expires_at: float
    user_id: str


def jwt_claim(access_token: str, claim: str) -> str:
    """One claim of the access token, decoded without verification — Underdog's
    API verifies it; this only names the account in the token file and logs."""
    try:
        payload = access_token.split(".")[1]
        payload += "=" * (-len(payload) % 4)
        claims = json.loads(base64.urlsafe_b64decode(payload))
        return str(claims[claim])
    except (IndexError, KeyError, ValueError) as error:
        raise UnderdogAuthError(f"access token has no readable {claim} claim: {error}") from error


def _token_request(body: dict) -> dict:
    """POST to Auth0's token endpoint with the standard library (`connect` runs
    from any Python on the VM). The error text never carries the request body."""
    request = urllib.request.Request(TOKEN_URL, data=json.dumps(body).encode(),
                                     headers={"Content-Type": "application/json", "Accept": "application/json"})
    try:
        with urllib.request.urlopen(request, timeout=HTTP_TIMEOUT_SEC) as response:
            status, text = response.status, response.read().decode()
    except urllib.error.HTTPError as error:
        status, text = error.code, error.read().decode(errors="replace")
    except urllib.error.URLError as error:
        raise UnderdogAuthError(f"Auth0 {body.get('grant_type')} unreachable: {error.reason}") from error
    try:
        payload = json.loads(text)
    except ValueError:
        payload = {}
    if status != 200 or "access_token" not in payload:
        detail = payload.get("error_description") or payload.get("error") or text[:200]
        raise UnderdogAuthError(f"Auth0 {body.get('grant_type')} failed ({status}): {detail}")
    return payload


def load_token_file(path: Path) -> dict | None:
    if not path.exists():
        return None
    try:
        stored = json.loads(path.read_text())
    except ValueError as error:
        raise UnderdogAuthError(f"{path} is not valid JSON: {error}") from error
    if not isinstance(stored, dict) or not stored.get("refresh_token"):
        raise UnderdogAuthError(f"{path} has no refresh_token; run `underdog_auth connect`")
    return stored


def save_token_file(path: Path, refresh_token: str, user_id: str) -> None:
    """Atomic rewrite, owner-only permissions."""
    tmp = path.with_suffix(path.suffix + ".tmp")
    tmp.write_text(json.dumps({"refresh_token": refresh_token, "user_id": user_id,
                               "obtained_at": time.strftime("%Y-%m-%dT%H:%M:%SZ", time.gmtime())}, indent=1))
    os.chmod(tmp, 0o600)
    os.replace(tmp, path)


class UnderdogAuth:
    """Holds the refresh token and hands out a live access token.

    `token()` refreshes when the cached access token is within
    REFRESH_MARGIN_SEC of expiry; a rotated refresh token is saved at once
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
            raise UnderdogAuthError(f"no Underdog token at {self._path}; run "
                                    "`python -m unabated_ticket.bets_service.sources.underdog_auth connect`")
        payload = self._request({"grant_type": "refresh_token", "client_id": AUTH0_CLIENT_ID,
                                 "refresh_token": stored["refresh_token"]})
        user_id = jwt_claim(payload["access_token"], "sub")
        rotated = payload.get("refresh_token")
        if rotated and rotated != stored["refresh_token"]:
            save_token_file(self._path, rotated, user_id)
            log.info("underdog: refresh token rotated and saved")
        life = float(payload.get("expires_in", DEFAULT_TOKEN_LIFE_SEC))
        self._access = AccessToken(payload["access_token"], self._clock() + life, user_id)
        return self._access

    def forget(self) -> None:
        """Drop the cached access token (the API refused it): the next token() refreshes."""
        self._access = None


def connect(token_path: Path, email: str, password: str, request=_token_request) -> str:
    """Log in once with the password grant and save the refresh token. Returns the account's user id."""
    payload = request({"grant_type": "password", "client_id": AUTH0_CLIENT_ID, "audience": AUTH0_AUDIENCE,
                       "scope": AUTH0_SCOPE, "username": email, "password": password})
    if not payload.get("refresh_token"):
        raise UnderdogAuthError("Auth0 returned no refresh_token (offline_access not granted)")
    user_id = jwt_claim(payload["access_token"], "sub")
    save_token_file(token_path, payload["refresh_token"], user_id)
    return user_id


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description="Connect the bets service to an Underdog account (one-time).")
    sub = parser.add_subparsers(dest="command", required=True)
    sub.add_parser("connect", help="log in once (email + password at a prompt) and save a refresh token")
    parser.parse_args(argv)
    email = input("Underdog email: ").strip()
    password = getpass.getpass("Underdog password (not stored): ")
    try:
        user_id = connect(config.UNDERDOG_TOKEN_PATH, email, password)
    except UnderdogAuthError as error:
        print(f"connect failed: {error}", file=sys.stderr)
        return 1
    print(f"\nSaved {config.UNDERDOG_TOKEN_PATH} for account {user_id}. "
          "Restart the bets service; the underdog source is now registered.")
    return 0


if __name__ == "__main__":
    sys.exit(main())
