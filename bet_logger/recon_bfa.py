"""
BFA Gaming Auth Capture Script
Performs Keycloak OIDC login via pure HTTP (no browser needed).

Usage:
    cd bet_logger
    ./venv/bin/python3 recon_bfa.py              # primary account
    ./venv/bin/python3 recon_bfa.py --account j   # second account (BFAJ)

Saves auth JSON with refresh_token + player_id for scraper_bfa.py.
"""

import os
import sys
import json
import re
import base64
import hashlib
import secrets
import argparse
from urllib.parse import urlparse, parse_qs
import requests
from dotenv import load_dotenv

load_dotenv()

SCRIPT_DIR = os.path.dirname(__file__)

ACCOUNTS = {
    'default': {
        'username_env': 'BFA_USERNAME',
        'password_env': 'BFA_PASSWORD',
        'auth_file': os.path.join(SCRIPT_DIR, 'recon_bfa_auth.json'),
    },
    'j': {
        'username_env': 'BFAJ_USERNAME',
        'password_env': 'BFAJ_PASSWORD',
        'auth_file': os.path.join(SCRIPT_DIR, 'recon_bfaj_auth.json'),
    },
}

KEYCLOAK_BASE = "https://auth.bfagaming.com/realms/players_realm/protocol/openid-connect"
CLIENT_ID = "bfagaming"
REDIRECT_URI = "https://bfagaming.com/"

USER_AGENT = (
    "Mozilla/5.0 (Macintosh; Intel Mac OS X 10_15_7) "
    "AppleWebKit/537.36 (KHTML, like Gecko) Chrome/145.0.0.0 Safari/537.36"
)


def decode_jwt_payload(token: str) -> dict:
    """Decode JWT payload without verification."""
    parts = token.split('.')
    if len(parts) != 3:
        return {}
    payload = parts[1]
    payload += '=' * (4 - len(payload) % 4)
    try:
        return json.loads(base64.urlsafe_b64decode(payload))
    except Exception:
        return {}


def generate_pkce():
    """Generate PKCE code_verifier and code_challenge."""
    verifier = secrets.token_urlsafe(32)
    challenge = base64.urlsafe_b64encode(
        hashlib.sha256(verifier.encode()).digest()
    ).rstrip(b'=').decode()
    return verifier, challenge


def run_recon(account: dict):
    username = os.getenv(account['username_env'])
    password = os.getenv(account['password_env'])
    auth_file = account['auth_file']

    if not username or not password:
        print(f"ERROR: {account['username_env']}/{account['password_env']} not set in .env")
        sys.exit(1)

    session = requests.Session()
    session.headers['User-Agent'] = USER_AGENT

    # Step 1: Generate PKCE and start auth flow
    code_verifier, code_challenge = generate_pkce()

    auth_url = (
        f"{KEYCLOAK_BASE}/auth"
        f"?client_id={CLIENT_ID}"
        f"&redirect_uri={REDIRECT_URI}"
        f"&response_type=code"
        f"&scope=openid"
        f"&code_challenge={code_challenge}"
        f"&code_challenge_method=S256"
    )

    print("Starting Keycloak auth flow...")
    resp = session.get(auth_url, allow_redirects=True, timeout=15)
    if resp.status_code != 200:
        print(f"ERROR: Failed to load Keycloak login page (HTTP {resp.status_code})")
        sys.exit(1)

    # Step 2: Extract the login form action URL
    action_match = re.search(r'action="([^"]+)"', resp.text)
    if not action_match:
        print("ERROR: Could not find login form action URL")
        sys.exit(1)

    login_url = action_match.group(1).replace('&amp;', '&')

    # Step 3: Submit credentials
    print("Submitting credentials...")
    resp = session.post(
        login_url,
        data={"username": username, "password": password},
        allow_redirects=False,
        timeout=15,
    )

    # Keycloak returns 302 redirect to BFA with ?code=xxx
    if resp.status_code not in (302, 303):
        print(f"ERROR: Login failed (HTTP {resp.status_code})")
        if "Invalid username or password" in resp.text:
            print("  Invalid credentials")
        sys.exit(1)

    redirect_url = resp.headers.get('Location', '')
    parsed = urlparse(redirect_url)
    code = parse_qs(parsed.query).get('code', [None])[0]

    if not code:
        print(f"ERROR: No auth code in redirect URL: {redirect_url[:100]}")
        sys.exit(1)

    print("Login successful, got auth code")

    # Step 4: Exchange auth code for tokens
    print("Exchanging code for tokens...")
    resp = session.post(
        f"{KEYCLOAK_BASE}/token",
        data={
            "grant_type": "authorization_code",
            "client_id": CLIENT_ID,
            "code": code,
            "redirect_uri": REDIRECT_URI,
            "code_verifier": code_verifier,
        },
        timeout=15,
    )

    if resp.status_code != 200:
        # Retry without PKCE in case Keycloak doesn't require it
        resp = session.post(
            f"{KEYCLOAK_BASE}/token",
            data={
                "grant_type": "authorization_code",
                "client_id": CLIENT_ID,
                "code": code,
                "redirect_uri": REDIRECT_URI,
            },
            timeout=15,
        )

    if resp.status_code != 200:
        print(f"ERROR: Token exchange failed (HTTP {resp.status_code})")
        print(f"  {resp.text[:200]}")
        sys.exit(1)

    tokens = resp.json()
    refresh_token = tokens.get("refresh_token")
    access_token = tokens.get("access_token")

    if not refresh_token:
        print("ERROR: No refresh_token in token response")
        sys.exit(1)

    # Extract player_id from JWT
    payload = decode_jwt_payload(access_token)
    player_id = payload.get("player_id")

    if not player_id:
        print("WARNING: No player_id in JWT, will need fallback")

    # Save auth
    auth_data = {
        "refresh_token": refresh_token,
        "player_id": str(player_id) if player_id else "",
    }
    with open(auth_file, "w") as f:
        json.dump(auth_data, f, indent=2)

    print(f"Saved auth to {auth_file}")
    print(f"  Player ID: {player_id}")
    print("Done.")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='Capture BFA auth tokens')
    parser.add_argument('--account', choices=list(ACCOUNTS.keys()), default='default',
                        help='Which BFA account to authenticate (default or j)')
    args = parser.parse_args()

    account = ACCOUNTS[args.account]
    print(f"Account: {args.account}")
    run_recon(account)
