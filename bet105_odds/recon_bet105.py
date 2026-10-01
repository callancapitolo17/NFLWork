"""
Bet105 recon: one manual login in plain Chrome, two captures.

Why plain Chrome first (2026-09-28): under Playwright's launch the site crawled,
because Playwright relays every frame of the odds socket (a nonstop stream) to
Python. So this script starts an ordinary Chrome on the persistent profile with a
debugging port and attaches only after you have logged in — the
mlb_sgp/quick_recon.py pattern. Sockets opened before the attach are not relayed.

You pass Cloudflare's check and the login's Turnstile yourself in that window; this
script never types credentials and saves no cookies (the bets service reads Bet105
through the Unabated Ticket extension in your own Chrome — unabated_ticket/README.md
§ Bets). Captures:
- The bets API calls the site makes while you click through My Plays
  -> .bet105_recon_api.json (request + response bodies: the shapes the bets parser
  is pinned to; cookie and CSRF header values redacted; local only). This is how
  unabated_ticket/tests/fixtures/bets/bet105_history.json was made (2026-09-29).
- WebSocket frames to pandora.ganchrow.com after one visit to the board -> the odds
  scraper's params (.bet105_session.json: prematch key, user id, group id, partner id).

The profile and both outputs live in bet105_odds/ of the MAIN checkout, even when this
script runs from a worktree. Both outputs are gitignored.

Usage:
    python recon_bet105.py
"""

import json
import re
import subprocess
import time
from pathlib import Path

from playwright.sync_api import sync_playwright

# My Plays is the site's My Bets view. The home page renders the live board and the
# casino carousels, which keep a CPU core busy the whole time you are logging in.
SITE_URL = "https://app.bet105.ag/my-plays"
# The odds scraper's params ride the prematch socket, which only the board subscribes
# to; My Plays does not. Visited once, at the end, after the bets calls are captured.
ODDS_PAGE_URL = "https://app.bet105.ag/sports/home"
SITE_HOST_MARKER = "bet105.ag"
CHROME_BINARY = "/Applications/Google Chrome.app/Contents/MacOS/Google Chrome"
# Not 9222: mlb_sgp/quick_recon.py attaches there, and two Chromes cannot share it.
DEBUG_PORT = 9235
DEBUG_PORT_WAIT_SEC = 30
# The profile was created by Playwright, whose Chrome encrypts cookies with a mock
# keychain; the same flags keep it readable (and avoid a macOS keychain prompt).
CHROME_ARGS = ["--no-first-run", "--no-default-browser-check", "--use-mock-keychain", "--password-store=basic"]
# Claude desktop puts worktrees under .claude/worktrees/<name>; the profile and the
# captures belong to the main checkout.
WORKTREE_MARKERS = ("/.claude/worktrees/", "/.worktrees/")
# The bets endpoints the site calls (read from its bundle 2026-09-27): the account
# session check (its body carries the CSRF token), and the LinePros "logic" POST
# whose `a` names the action (getHistory is the open bets).
CUSTOMERS_PATH = "/__bff/api/customers"
BETS_API_MARKERS = (CUSTOMERS_PATH, "/betLobbyV2/logic/", "/betLobbyV3/")
REDACTED_HEADERS = ("cookie", "x-broker-csrf")
HISTORY_WAIT_SEC = 20
ODDS_PARAMS_WAIT_SEC = 20


def main_checkout_dir(here: Path) -> Path:
    """This directory's twin in the main checkout when `here` is inside a worktree."""
    path = f"{here}/"
    for marker in WORKTREE_MARKERS:
        if marker in path:
            main_root, rest = path.split(marker, 1)
            inside_worktree = rest.split("/", 1)[1]  # drop the worktree's own name
            return Path(main_root) / inside_worktree.rstrip("/")
    return here


OUT_DIR = main_checkout_dir(Path(__file__).resolve().parent)
PROFILE_DIR = OUT_DIR / ".bet105_profile"
SESSION_PATH = OUT_DIR / ".bet105_session.json"
API_CAPTURE_PATH = OUT_DIR / ".bet105_recon_api.json"


def utc_now_iso() -> str:
    return time.strftime("%Y-%m-%dT%H:%M:%SZ", time.gmtime())


def launch_chrome() -> subprocess.Popen:
    """An ordinary Chrome on the profile with a debugging port; nothing attached yet."""
    return subprocess.Popen(
        [CHROME_BINARY, f"--user-data-dir={PROFILE_DIR}", f"--remote-debugging-port={DEBUG_PORT}",
         *CHROME_ARGS, SITE_URL],
        stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)


def connect(playwright):
    """Attach over CDP, retrying while Chrome opens its debugging port."""
    deadline = time.monotonic() + DEBUG_PORT_WAIT_SEC
    while True:
        try:
            return playwright.chromium.connect_over_cdp(f"http://127.0.0.1:{DEBUG_PORT}")
        except Exception:  # noqa: BLE001 — the port is not up yet
            if time.monotonic() > deadline:
                raise RuntimeError(f"Chrome did not open debugging port {DEBUG_PORT} within "
                                   f"{DEBUG_PORT_WAIT_SEC}s (is another Chrome using {PROFILE_DIR}?)")
            time.sleep(0.5)


def site_page(context):
    pages = [page for page in context.pages if SITE_HOST_MARKER in page.url]
    if not pages:
        raise RuntimeError(f"no {SITE_HOST_MARKER} tab open; open {SITE_URL} in that Chrome window")
    return pages[0]


def capture_ws_params(text: str, captured: dict) -> None:
    """Outgoing Socket.IO frames carry the odds scraper's params."""
    # subscribeSystemEvents contains userId, groupId, partnerId
    if "subscribeSystemEvents" in text:
        match = re.search(r'\["subscribeSystemEvents"\s*,\s*(\{[^}]+\})\]', text)
        if match:
            try:
                data = json.loads(match.group(1))
            except json.JSONDecodeError:
                return
            captured["user_id"] = data.get("userId")
            captured["group_id"] = data.get("groupId")
            captured["partner_id"] = str(data.get("partnerId", "111"))
    # the eventData room subscription carries the prematch key
    if "prematch.main." in text and ".eventData" in text:
        match = re.search(r'prematch\.main\.([A-Za-z0-9+/=]+)\.eventData', text)
        if match:
            captured["prematch_key"] = match.group(1)


def is_bets_api(url: str) -> bool:
    return any(marker in url for marker in BETS_API_MARKERS)


def action_of(post_data: str | None) -> str:
    if not post_data:
        return ""
    try:
        return str(json.loads(post_data).get("a", ""))
    except (json.JSONDecodeError, AttributeError):
        return ""


def redacted_headers(headers: dict) -> dict:
    return {name: ("<redacted>" if name.lower() in REDACTED_HEADERS else value) for name, value in headers.items()}


def call_entry(response) -> dict:
    """One bets API request + its response body. Read after the clicks, not inside the
    event callback (a body read there deadlocks the sync API — quick_recon.py)."""
    request = response.request
    entry = {"method": request.method, "url": response.url, "requestHeaders": redacted_headers(request.headers),
             "postData": request.post_data, "status": response.status, "body": None}
    try:
        text = response.text()
    except Exception as error:  # noqa: BLE001 — a body that cannot be read is still a captured call
        entry["bodyError"] = f"{type(error).__name__}: {error}"
        return entry
    if CUSTOMERS_PATH in response.url:
        # The account's own profile: keep the shape the extension reads, not the person.
        try:
            body = json.loads(text)
            entry["bodyKeys"] = sorted(body.keys()) if isinstance(body, dict) else type(body).__name__
            entry["hasCsrfToken"] = isinstance(body, dict) and bool(body.get("csrfToken"))
        except json.JSONDecodeError:
            entry["bodyKeys"] = "not JSON"
        return entry
    entry["body"] = text
    return entry


def save_api_capture(calls: list) -> None:
    API_CAPTURE_PATH.write_text(json.dumps({"capturedAt": utc_now_iso(), "calls": calls}, indent=2))
    API_CAPTURE_PATH.chmod(0o600)
    actions = sorted({action_of(call["postData"]) for call in calls} - {""})
    print(f"Saved {len(calls)} bets API calls to {API_CAPTURE_PATH} (actions: {', '.join(actions) or 'none'})")
    if "getHistory" not in actions:
        print("WARNING: no getHistory call captured — was My Plays open when you pressed ENTER?")


def save_session_params(captured: dict) -> None:
    if not all(captured.values()):
        missing = [key for key, value in captured.items() if not value]
        print(f"Odds params NOT saved; missing {missing}.")
        return
    session = {
        "prematch_key": captured["prematch_key"],
        "user_id": int(captured["user_id"]),
        "group_id": int(captured["group_id"]),
        "partner_id": captured["partner_id"],
        "captured_at": utc_now_iso(),
    }
    SESSION_PATH.write_text(json.dumps(session, indent=2))
    print(f"Saved odds params to {SESSION_PATH}")


def capture_history_calls(page, responses: list) -> None:
    """Reload My Plays so getHistory fires WHILE we listen (its first call fired
    during login, before we attached), then wait for the open-bets call to land."""
    responses.clear()
    page.goto(SITE_URL, wait_until="domcontentloaded", timeout=60000)
    deadline = time.monotonic() + HISTORY_WAIT_SEC
    while time.monotonic() < deadline:
        if any(action_of(response.request.post_data) == "getHistory" for response in responses):
            break
        page.wait_for_timeout(500)


def capture_odds_params(page) -> dict:
    """Go to the board with a socket listener; stop once the four params are in."""
    captured = {"prematch_key": None, "user_id": None, "group_id": None, "partner_id": None}

    def on_websocket(ws):
        if "ganchrow.com" in ws.url:
            ws.on("framesent", lambda payload: capture_ws_params(
                payload if isinstance(payload, str) else payload.decode("utf-8", errors="ignore"), captured))

    page.on("websocket", on_websocket)
    page.goto(ODDS_PAGE_URL, wait_until="domcontentloaded", timeout=60000)
    deadline = time.monotonic() + ODDS_PARAMS_WAIT_SEC
    while not all(captured.values()) and time.monotonic() < deadline:
        page.wait_for_timeout(500)
    page.remove_listener("websocket", on_websocket)
    return captured


def run_recon() -> None:
    chrome = launch_chrome()
    try:
        print("\n" + "=" * 60)
        print("STEP 1 — in the Chrome window that just opened:")
        print("  pass the Cloudflare check and log in (MFA if asked).")
        print("  Press ENTER here when you see My Plays logged in.")
        print("=" * 60)
        input()

        with sync_playwright() as p:
            browser = connect(p)
            context = browser.contexts[0]
            page = site_page(context)
            responses = []
            page.on("response", lambda response: responses.append(response) if is_bets_api(response.url) else None)

            print("Reloading My Plays so the bets calls fire while we listen...")
            capture_history_calls(page, responses)

            print("=" * 60)
            print("STEP 2 — click any other My Plays view you want captured, then press ENTER here.")
            print("=" * 60)
            input()

            calls = [call_entry(response) for response in responses]
            print("Visiting the board once for the odds params...")
            odds_params = capture_odds_params(page)
            browser.close()
    finally:
        chrome.terminate()

    print()
    save_api_capture(calls)
    save_session_params(odds_params)


if __name__ == "__main__":
    run_recon()
