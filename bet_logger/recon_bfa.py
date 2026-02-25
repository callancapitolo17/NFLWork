"""
BFA Gaming Auth Capture Script
Launches Chrome to capture Keycloak auth tokens for the API scraper.
Uses persistent profile + auto-login from .env credentials.
Fully automated — no manual interaction needed.

Usage:
    cd bet_logger
    ./venv/bin/python3 recon_bfa.py

Saves recon_bfa_auth.json with refresh_token + player_id for scraper_bfa.py.
"""

from playwright.sync_api import sync_playwright
import os
import sys
import json
import re
import base64
from dotenv import load_dotenv

load_dotenv()

BFA_URL = "https://bfagaming.com"
BFA_MY_BETS = "https://bfagaming.com/my-bets"
SCRIPT_DIR = os.path.dirname(__file__)
PROFILE_DIR = os.path.join(SCRIPT_DIR, ".bfa_bet_profile")
AUTH_FILE = os.path.join(SCRIPT_DIR, "recon_bfa_auth.json")

BFA_USERNAME = os.getenv("BFA_USERNAME")
BFA_PASSWORD = os.getenv("BFA_PASSWORD")


def decode_jwt_payload(token: str) -> dict:
    """Decode JWT payload without verification (just base64)."""
    parts = token.split('.')
    if len(parts) != 3:
        return {}
    payload = parts[1]
    payload += '=' * (4 - len(payload) % 4)
    try:
        return json.loads(base64.urlsafe_b64decode(payload))
    except Exception:
        return {}


def run_recon():
    captured = {"refresh_token": None, "access_token": None, "player_id": None}

    def handle_response(response):
        """Capture token endpoint responses."""
        if response.request.resource_type not in ("xhr", "fetch"):
            return

        url = response.url

        # Capture Keycloak token refresh response
        if "openid-connect/token" in url and response.status == 200:
            try:
                body = response.text()
                data = json.loads(body)
                if "refresh_token" in data:
                    captured["refresh_token"] = data["refresh_token"]
                    captured["access_token"] = data.get("access_token", "")
                    print(f"  [AUTH] Captured refresh token")

                    # Extract player_id from JWT payload
                    payload = decode_jwt_payload(captured["access_token"])
                    pid = payload.get("player_id")
                    if pid:
                        captured["player_id"] = pid
                        print(f"  [AUTH] Player ID: {pid}")
            except Exception:
                pass

        # Fallback: get player_id from API URL params
        if "playerId=" in url and not captured["player_id"]:
            pid_match = re.search(r'playerId=(\d+)', url)
            if pid_match:
                captured["player_id"] = pid_match.group(1)
                print(f"  [API]  Player ID: {captured['player_id']}")

    with sync_playwright() as p:
        context = p.chromium.launch_persistent_context(
            user_data_dir=PROFILE_DIR,
            channel="chrome",
            headless=True,
            viewport={"width": 1400, "height": 900},
            args=["--disable-blink-features=AutomationControlled"],
        )
        page = context.pages[0] if context.pages else context.new_page()
        page.on("response", handle_response)

        # Navigate to BFA
        print("Navigating to BFA...")
        page.goto(BFA_URL, wait_until="domcontentloaded", timeout=60000)
        page.wait_for_timeout(3000)

        # Check if login is needed
        needs_login = page.locator('text=You must log in').count() > 0
        login_btn = page.locator(
            '.header-auth-container button:has-text("Log In"), '
            'button:has-text("Log in")'
        )

        if needs_login or login_btn.count() > 0:
            if not BFA_USERNAME or not BFA_PASSWORD:
                print("ERROR: Login required but BFA_USERNAME/BFA_PASSWORD not set in .env")
                context.close()
                sys.exit(1)

            print("Login required. Auto-logging in via Keycloak...")

            # Click login button → redirects to Keycloak
            header_login = page.locator('.header-auth-container button:has-text("Log In")')
            if header_login.count() > 0:
                header_login.click()
            else:
                login_btn.first.click()

            page.wait_for_url('**/realms/**', timeout=15000)
            page.wait_for_timeout(2000)

            # Fill Keycloak login form
            username_field = page.locator('#username')
            if username_field.count() > 0:
                username_field.fill(BFA_USERNAME)
            else:
                page.fill('input[name="username"]', BFA_USERNAME)

            page.wait_for_timeout(500)

            password_field = page.locator('#password')
            if password_field.count() > 0:
                password_field.fill(BFA_PASSWORD)
            else:
                page.fill('input[name="password"]', BFA_PASSWORD)

            page.wait_for_timeout(500)

            # Submit
            login_submit = page.locator('#kc-login')
            if login_submit.count() > 0:
                login_submit.click()
            else:
                page.click('input[type="submit"]')

            # Wait for redirect back to BFA
            page.wait_for_url('**/bfagaming.com/**', timeout=30000)
            page.wait_for_load_state('networkidle', timeout=30000)
            page.wait_for_timeout(3000)
            print("Login successful")
        else:
            print("Already logged in (persistent session)")

        # Navigate to my-bets — triggers token refresh automatically
        print("Navigating to my-bets (captures token refresh)...")
        page.goto(BFA_MY_BETS, wait_until="domcontentloaded", timeout=60000)
        page.wait_for_timeout(8000)

        # Save auth if captured
        if captured["refresh_token"] and captured["player_id"]:
            auth_data = {
                "refresh_token": captured["refresh_token"],
                "player_id": captured["player_id"],
            }
            with open(AUTH_FILE, "w") as f:
                json.dump(auth_data, f, indent=2)
            print(f"Saved auth to {AUTH_FILE}")
        else:
            print("WARNING: Could not capture auth data.")
            if not captured["refresh_token"]:
                print("  Missing refresh token — session may have expired")
            if not captured["player_id"]:
                print("  Missing player_id")
            context.close()
            sys.exit(1)

        context.close()
        print("Done.")


if __name__ == "__main__":
    run_recon()
