"""
Bet105.ag Recon Script
- Uses real Chrome (not Playwright Chromium) to bypass Cloudflare
- Persistent profile so login survives between runs
- Intercepts WebSocket to pandora.ganchrow.com to extract auth params
- Saves session to .bet105_session.json for the scraper

Usage:
    python recon_bet105.py

First run: browser opens, you log in manually (handle MFA).
Subsequent runs: reuses saved browser profile, extracts params automatically.
"""

import json
import os
import re
import sys
import time

from playwright.sync_api import sync_playwright

SITE_URL = "https://sportsbook.bet105.ag"
PREMATCH_URL = "https://sportsbook.bet105.ag/#prematch"
PROFILE_DIR = os.path.join(os.path.dirname(__file__), ".bet105_profile")
SESSION_PATH = os.path.join(os.path.dirname(__file__), ".bet105_session.json")


def run_recon():
    captured = {
        "prematch_key": None,
        "user_id": None,
        "group_id": None,
        "partner_id": None,
    }

    def on_ws_frame_sent(payload):
        """Parse outgoing WebSocket frames to extract auth params."""
        text = payload if isinstance(payload, str) else payload.decode("utf-8", errors="ignore")

        # subscribeSystemEvents contains userId, groupId, partnerId
        if "subscribeSystemEvents" in text:
            try:
                # Extract JSON array from Socket.IO frame: 42["subscribeSystemEvents",{...}]
                match = re.search(r'\["subscribeSystemEvents"\s*,\s*(\{[^}]+\})\]', text)
                if match:
                    data = json.loads(match.group(1))
                    captured["user_id"] = data.get("userId")
                    captured["group_id"] = data.get("groupId")
                    captured["partner_id"] = str(data.get("partnerId", "111"))
                    print(f"  Captured: userId={captured['user_id']}, "
                          f"groupId={captured['group_id']}, partnerId={captured['partner_id']}")
            except (json.JSONDecodeError, AttributeError):
                pass

        # subscribe to eventData room contains prematchKey
        if "prematch.main." in text and ".eventData" in text:
            try:
                match = re.search(r'prematch\.main\.([A-Za-z0-9+/=]+)\.eventData', text)
                if match:
                    captured["prematch_key"] = match.group(1)
                    print(f"  Captured: prematchKey={captured['prematch_key']}")
            except AttributeError:
                pass

    with sync_playwright() as p:
        context = p.chromium.launch_persistent_context(
            user_data_dir=PROFILE_DIR,
            channel="chrome",
            headless=False,
            viewport={"width": 1400, "height": 900},
            args=["--disable-blink-features=AutomationControlled"],
        )
        page = context.pages[0] if context.pages else context.new_page()

        # Attach WebSocket listener
        def on_websocket(ws):
            if "pandora.ganchrow.com" in ws.url:
                print(f"  WebSocket opened: {ws.url[:80]}")
                ws.on("framesent", lambda payload: on_ws_frame_sent(payload))

        page.on("websocket", on_websocket)

        # Navigate to bet105
        print("Navigating to Bet105.ag...")
        page.goto(SITE_URL, wait_until="domcontentloaded", timeout=60000)

        # Check if already logged in by looking for common logged-in indicators
        print("\n" + "=" * 60)
        print("1. Wait for Cloudflare check to pass")
        print("2. Log in if needed (handle MFA if prompted)")
        print("3. Navigate to the SPORTS / PREMATCH section")
        print("4. Wait for odds to load, then press ENTER here")
        print("=" * 60)
        input()

        # If we haven't captured params yet, try navigating to prematch
        if not captured["prematch_key"]:
            print("\nNavigating to prematch section...")
            page.goto(PREMATCH_URL, wait_until="domcontentloaded", timeout=60000)
            # Wait for WebSocket to connect and send auth
            time.sleep(10)

        # Check if we got everything
        if all(captured.values()):
            print(f"\nAll params captured successfully!")
        else:
            missing = [k for k, v in captured.items() if not v]
            print(f"\nMissing params: {missing}")
            print("Try clicking on different sports pages to trigger WebSocket.")
            print("Press ENTER when ready.")
            input()
            time.sleep(5)

        context.close()

    if all(captured.values()):
        session = {
            "prematch_key": captured["prematch_key"],
            "user_id": int(captured["user_id"]),
            "group_id": int(captured["group_id"]),
            "partner_id": captured["partner_id"],
            "captured_at": time.strftime("%Y-%m-%dT%H:%M:%SZ", time.gmtime()),
        }
        with open(SESSION_PATH, "w") as f:
            json.dump(session, f, indent=2)
        print(f"\nSaved session to {SESSION_PATH}")
        print(json.dumps(session, indent=2))
        return session
    else:
        print("\nFailed to capture all params. Session not saved.")
        return None


if __name__ == "__main__":
    run_recon()
