#!/usr/bin/env python3
"""
Caesars Sportsbook SGP Network Recon Tool
==========================================

GOAL of this recon: determine whether we can programmatically
  (a) list today's MLB events,
  (b) discover selection IDs for a run-line (spread) + total (over/under),
  (c) obtain a priced SGP slip with correlated odds.

Caesars runs the Liberty / American Wagering platform (Caesars acquired
William Hill US). It uses the ZeroFlucs correlation engine in-house for SGP.

---------------------------------------------------------------------------
WHAT THE 2026-06 RECON FOUND (see REPORT in the conversation for full detail)
---------------------------------------------------------------------------
API gateway host:   api.americanwagering.com   (fronted by AWS CloudFront)
API path family:    /regions/{country}/locations/{region}/brands/czr/sb/v3/...
                    (also a newer /v4/ feed: /v4/sports/{sportId}/competitions/...)
Event detail:       GET .../sb/v3/events/{event_id}
                      -> JSON {name, markets:[{id,name,active,display,line,
                               selections:[{id,name,active,display,price:{d}}]}]}
Betslip deeplink:   https://sportsbook.caesars.com/us/{state}/bet/betslip?selectionIds=<uuid>[,<uuid>...]
Realtime feed:      Push Technology Diffusion over STOMP-over-WebSocket (wss://),
                    not plain REST. Prices + SGP combined odds arrive as topic
                    pushes, not a single REST response.

TWO HARD WALLS (both confirmed live, 2026-06-11):
  1. AWS WAF token gate. EVERY request to api.americanwagering.com returns
     HTTP 403 (CloudFront "request could not be satisfied") unless it carries
     an `X-Aws-Waf-Token` header. The SPA mints that token by running the AWS
     WAF captcha SDK in-browser (4ad3fec456d9.edge.sdk.awswaf.com). Confirmed
     in the bundle: every fetch is wrapped as
        fetch(url, {...headers, "X-Aws-Waf-Token": <token>})
     curl_cffi impersonate="chrome" alone does NOT pass (still 403) because
     the token requires executing the challenge JS.
  2. GeoComply PLC (Player Location Check). The bundle pulls GeoComply PLC
     heavily (logger.geocomply.net, myip.geocomply.com, solution "PLC").
     Pricing/placement is region-gated behind a geolocation token.

CONCLUSION: a pure curl_cffi/requests REST scraper is NOT viable for Caesars.
The only workable path is a real browser (Playwright) that lets the AWS WAF
SDK mint the token, then scrapes the Diffusion WS pushes / betslip DOM, or
re-uses the short-lived WAF token for direct REST event-detail reads.

This file ships TWO entry points:
  * probe_rest()         (--rest) demonstrates the 403 WAF wall, no browser.
  * run_browser_recon()  (default) Playwright capture mirroring
                         recon_draftkings_sgp.py / recon_novig_sgp.py.

Usage:
    cd mlb_sgp
    python recon_caesars_sgp.py --rest      # re-confirm the WAF 403 wall (fast)
    python recon_caesars_sgp.py             # browser capture (default)
    python recon_caesars_sgp.py --verbose

Browser capture flow (you drive it):
    1. Browser opens sportsbook.caesars.com (set your legal state in the URL).
    2. Allow GeoComply / location if prompted; log in if needed.
    3. Navigate to an MLB game.
    4. Open the Same Game Parlay / Parlay Builder.
    5. Add a run-line (spread) leg, then a total (over/under) leg.
    6. Watch the [MATCH]/[WS] lines for the SGP combined-odds push.
    7. Everything is saved to recon_caesars_sgp.json.
"""

import argparse
import json
import os

STATE = os.environ.get("CAESARS_STATE", "nj")
CAESARS_URL = f"https://sportsbook.caesars.com/us/{STATE}/bet/baseball"

API_HOST = "api.americanwagering.com"
EVENT_DETAIL_TMPL = (
    "https://api.americanwagering.com/regions/us/locations/{jurisdiction}"
    "/brands/czr/sb/v3/events/{event_id}"
)

PROFILE_DIR = os.environ.get(
    "CAESARS_PROFILE", os.path.join(os.path.dirname(__file__), ".caesars_profile")
)
OUT_PATH = os.path.join(os.path.dirname(__file__), "recon_caesars_sgp.json")

TRACKED_DOMAINS = (
    "americanwagering.com",
    "caesars.com",
    "l5y.app",
    "awswaf.com",
    "geocomply",
)

SGP_KEYWORDS = [
    "sgp", "parlay", "combined", "combinedodds", "combined_odds",
    "correlat", "selection", "selectionids", "betslip", "slip", "wager",
    "quote", "payout", "price", "builder", "sgm", "multi", "accumulator",
    "potentialpayout", "shortcode", "diffusion",
]

JSON_PREVIEW_LIMIT = 12_000
TAIL_GRACE_MS = 400


def probe_rest():
    from curl_cffi import requests as cr

    hdrs = {
        "User-Agent": "Mozilla/5.0 (Macintosh; Intel Mac OS X 10_15_7) "
        "AppleWebKit/537.36 (KHTML, like Gecko) Chrome/120.0 Safari/537.36",
        "Accept": "application/json",
        "Origin": "https://sportsbook.caesars.com",
        "Referer": "https://sportsbook.caesars.com/",
    }

    print("=" * 70)
    print("REST PROBE -- api.americanwagering.com (expect 403 AWS-WAF wall)")
    print("=" * 70)

    candidates = [
        f"https://{API_HOST}/regions/us/locations/{STATE}/brands/czr/sb/v3/sports",
        f"https://{API_HOST}/regions/us/locations/{STATE}/brands/czr/sb/v3/sports/baseball",
        f"https://{API_HOST}/regions/us/locations/{STATE}/brands/czr/sb/v3/eventgroups",
        EVENT_DETAIL_TMPL.format(
            jurisdiction=STATE, event_id="00000000-0000-0000-0000-000000000000"
        ),
    ]

    for u in candidates:
        try:
            r = cr.get(u, headers=hdrs, impersonate="chrome", timeout=20)
            waf = r.headers.get("x-amzn-waf-action")
            body = r.text[:120].replace("\n", " ")
            blocked = "Request blocked" in r.text or "could not be satisfied" in r.text
            tag = "WAF-BLOCKED" if (r.status_code == 403 and blocked) else "?"
            print(f"\n[{r.status_code}] {tag} waf-action={waf}")
            print(f"    {u}")
            print(f"    {body}")
        except Exception as e:
            print(f"\nERR {u}\n    {e}")

    print("\n--- AWS WAF token issuer reachability ---")
    for h in (
        "https://4ad3fec456d9.edge.sdk.awswaf.com/4ad3fec456d9/",
        "https://4ad3fec456d9.edge.captcha-sdk.awswaf.com/4ad3fec456d9/",
    ):
        try:
            rr = cr.get(h, impersonate="chrome", timeout=15)
            print(f"[{rr.status_code}] {h}  (live = needs challenge params)")
        except Exception as e:
            print(f"ERR {h} {e}")

    print(
        "\nVERDICT: every gateway path 403s without X-Aws-Waf-Token. "
        "A headless REST scraper cannot pass. Use the browser capture (default)."
    )


def run_auto_recon(verbose_default=False):
    """
    Headless / automated recon. No human drives the betslip.

    Flow:
      1. Launch Chrome (headed=False unless CAESARS_HEADED=1), navigate to the
         Caesars sportsbook for STATE. Let the AWS WAF SDK mint the
         `aws-waf-token` cookie in-browser.
      2. Capture EVERY request to api.americanwagering.com and record which
         headers the SPA attaches (esp. X-Aws-Waf-Token) + which paths return
         200 vs 403. This reveals the live API shape and whether the line
         surface is reachable from THIS IP at all.
      3. Pull the aws-waf-token cookie and the live request headers, then
         REPLAY a direct REST event-detail read from inside the page context
         (page.request) and from a fresh curl_cffi session, to prove whether
         token-replay works headlessly.

    Everything is dumped to recon_caesars_sgp.json.
    """
    from playwright.sync_api import sync_playwright

    captured = {
        "api_requests": [],   # every api.americanwagering.com request seen
        "waf_token": None,
        "cookies": [],
        "events": [],
        "replay": {},
        "ws": [],
        "geo": {},
    }

    headed = os.environ.get("CAESARS_HEADED", "0") == "1"

    def _short(url, n=140):
        for p in ("https://", "http://"):
            if url.startswith(p):
                url = url[len(p):]
                break
        return url if len(url) <= n else url[: n - 3] + "..."

    def on_request(req):
        u = req.url
        if "americanwagering.com" not in u:
            return
        hdrs = req.headers
        captured["api_requests"].append({
            "method": req.method,
            "url": u,
            "has_waf": "x-aws-waf-token" in {k.lower() for k in hdrs},
            "waf_token": hdrs.get("x-aws-waf-token") or hdrs.get("X-Aws-Waf-Token"),
            "req_headers": {k: (v[:80] if k.lower() in ("x-aws-waf-token", "authorization", "cookie") else v)
                            for k, v in hdrs.items()},
            "post_data": (req.post_data or "")[:600] if req.method != "GET" else None,
        })

    def on_response(resp):
        u = resp.url
        if "americanwagering.com" not in u:
            return
        for r in captured["api_requests"]:
            if r["url"] == u and "status" not in r:
                r["status"] = resp.status
                ct = resp.headers.get("content-type", "")
                r["content_type"] = ct
                # capture a JSON body sample for any 200 that looks like odds
                if resp.status == 200 and "json" in ct:
                    try:
                        body = resp.text()
                        r["body_len"] = len(body)
                        low = body.lower()
                        if any(k in low for k in ("selection", "market", "price", "event", "competition")):
                            r["body_sample"] = body[:6000]
                            captured["events"].append({"url": u, "body": body})
                    except Exception:
                        pass
                elif resp.status in (403, 405):
                    try:
                        b = resp.text()
                        r["block_body"] = b[:300]
                        r["waf_action"] = resp.headers.get("x-amzn-waf-action")
                    except Exception:
                        pass
                break

    def on_websocket(ws):
        rec = {"url": ws.url, "frames": []}
        captured["ws"].append(rec)
        print(f"  [WS] {_short(ws.url)}")

        def on_frame(payload, direction):
            try:
                text = payload if isinstance(payload, str) else payload.decode("utf-8", "ignore")
            except Exception:
                return
            low = text.lower()
            if any(kw in low for kw in SGP_KEYWORDS):
                rec["frames"].append({"dir": direction, "data": text[:3000]})

        ws.on("framereceived", lambda p: on_frame(p, "recv"))
        ws.on("framesent", lambda p: on_frame(p, "sent"))

    with sync_playwright() as pw:
        REAL_UA = ("Mozilla/5.0 (Macintosh; Intel Mac OS X 10_15_7) "
                   "AppleWebKit/537.36 (KHTML, like Gecko) "
                   "Chrome/149.0.0.0 Safari/537.36")
        ctx = pw.chromium.launch_persistent_context(
            user_data_dir=PROFILE_DIR,
            channel="chrome",
            headless=not headed,
            viewport={"width": 1400, "height": 900},
            user_agent=REAL_UA,
            args=[
                "--disable-blink-features=AutomationControlled",
                f"--user-agent={REAL_UA}",
            ],
        )
        # Stealth: strip navigator.webdriver before any page script runs.
        ctx.add_init_script(
            "Object.defineProperty(navigator,'webdriver',{get:()=>undefined});"
        )
        page = ctx.pages[0] if ctx.pages else ctx.new_page()
        page.on("request", on_request)
        page.on("response", on_response)
        page.on("websocket", on_websocket)

        print("=" * 70)
        print(f"AUTO RECON  (headed={headed}, STATE={STATE})")
        print(f"  -> {CAESARS_URL}")
        print("=" * 70)
        try:
            page.goto(CAESARS_URL, wait_until="commit", timeout=90000)
        except Exception as e:
            print(f"  nav: {e}")

        # Give the WAF SDK time to run its challenge + the SPA to fire API calls.
        for i in range(8):
            page.wait_for_timeout(2500)
            cookies = ctx.cookies()
            waf = next((c for c in cookies if c["name"] == "aws-waf-token"), None)
            n_api = len(captured["api_requests"])
            n_200 = sum(1 for r in captured["api_requests"] if r.get("status") == 200)
            n_403 = sum(1 for r in captured["api_requests"] if r.get("status") == 403)
            print(f"  t+{(i+1)*2.5:.0f}s  waf_cookie={'YES' if waf else 'no '}  "
                  f"api_reqs={n_api} (200={n_200} 403={n_403})")
            if waf:
                captured["waf_token"] = waf["value"]
            if waf and n_200 >= 2:
                break

        # Re-navigate now that the WAF token + cookies are warm, to force the
        # SPA to fire the real MLB odds feed (first load often sits on splash).
        for path in ("/bet/baseball", "/bet/baseball/competition/mlb-04f8d"):
            try:
                page.goto(f"https://sportsbook.caesars.com/us/{STATE}{path}",
                          wait_until="commit", timeout=45000)
                page.wait_for_timeout(6000)
            except Exception as e:
                print(f"  renav {path}: {e}")
            n_200 = sum(1 for r in captured["api_requests"] if r.get("status") == 200)
            print(f"  after renav {path}: 200s={n_200} events={len(captured['events'])}")
            if len(captured["events"]) > 2:
                break

        # Try to navigate to the baseball / MLB listing to force odds API calls.
        try:
            page.wait_for_timeout(1500)
            # Record final URL + page title to detect geo redirect / "not available"
            captured["geo"]["final_url"] = page.url
            captured["geo"]["title"] = page.title()
            body_text = page.evaluate("() => document.body ? document.body.innerText.slice(0,1200) : ''")
            captured["geo"]["body_text"] = body_text
            print(f"\n  final_url: {page.url}")
            print(f"  title: {captured['geo']['title']}")
            print(f"  body[:400]: {body_text[:400]!r}")
        except Exception as e:
            print(f"  geo probe: {e}")

        captured["cookies"] = [
            {"name": c["name"], "domain": c["domain"], "len": len(c["value"])}
            for c in ctx.cookies()
        ]

        # ---- REPLAY: use the in-page request context (carries WAF cookie) ----
        # Find any event_id seen in captured bodies/urls to test a real read.
        import re as _re
        uuid_re = _re.compile(r"[a-f0-9]{8}-[a-f0-9]{4}-[a-f0-9]{4}-[a-f0-9]{4}-[a-f0-9]{12}", _re.I)
        seen_ids = []
        for ev in captured["events"]:
            seen_ids += uuid_re.findall(ev["body"])
        seen_ids = list(dict.fromkeys(seen_ids))[:5]
        captured["replay"]["candidate_event_ids"] = seen_ids

        base = f"https://{API_HOST}/regions/us/locations/{STATE}/brands/czr"
        test_urls = [
            f"{base}/sb/v3/sports",
            f"{base}/sb/v4/sports/baseball/tabs",
            f"{base}/sb/v4/sports/baseball/quick-picks",
            f"{base}/sb/v4/navigation-items",
            f"{base}/gw/growth/v4/sports/baseball/banners",
        ]
        for eid in seen_ids[:3]:
            test_urls.append(f"{base}/sb/v3/events/{eid}")
            test_urls.append(f"{base}/sb/v4/events/{eid}")

        # Reconstruct the exact header set the SPA attaches to successful odds
        # calls (discovered: x-platform / x-app-version / x-unique-device-id are
        # required by the CloudFront rule, alongside the WAF token).
        spa_hdrs = {}
        for r in captured["api_requests"]:
            if r.get("status") == 200 and r.get("has_waf"):
                for hk, hv in r.get("req_headers", {}).items():
                    if hk.lower() in ("x-platform", "x-app-version",
                                      "x-unique-device-id", "x-aws-waf-token",
                                      "content-type"):
                        spa_hdrs[hk] = hv
        # x-aws-waf-token in req_headers is truncated to 80 chars; use live cookie.
        if captured.get("waf_token"):
            spa_hdrs["x-aws-waf-token"] = captured["waf_token"]
        spa_hdrs.setdefault("x-platform", "cordova-desktop")
        spa_hdrs.setdefault("x-app-version", "7.49.0")
        spa_hdrs["accept"] = "application/json"
        captured["replay"]["spa_headers_used"] = {
            k: (v[:20] + "..." if k == "x-aws-waf-token" else v)
            for k, v in spa_hdrs.items()
        }

        captured["replay"]["in_page"] = []
        for u in test_urls:
            try:
                resp = page.request.get(u, headers=spa_hdrs, timeout=20000)
                txt = resp.text()
                captured["replay"]["in_page"].append({
                    "url": u, "status": resp.status, "len": len(txt),
                    "sample": txt[:1500],
                })
                print(f"  [REPLAY in-page] {resp.status}  {_short(u)}  ({len(txt)}B)")
            except Exception as e:
                captured["replay"]["in_page"].append({"url": u, "error": str(e)[:200]})
                print(f"  [REPLAY in-page] ERR {_short(u)}  {e}")

        # ---- REPLAY: fresh curl_cffi with the WAF cookie + observed headers ----
        if captured["waf_token"]:
            try:
                from curl_cffi import requests as cr
                hdrs = {
                    "User-Agent": page.evaluate("() => navigator.userAgent"),
                    "Accept": "application/json",
                    "Origin": "https://sportsbook.caesars.com",
                    "Referer": "https://sportsbook.caesars.com/",
                    "X-Aws-Waf-Token": captured["waf_token"],
                }
                hdrs.update(spa_hdrs)
                cookie_str = "; ".join(
                    f"{c['name']}={c['value']}" for c in ctx.cookies()
                    if "americanwagering" in c["domain"] or c["name"] == "aws-waf-token"
                )
                captured["replay"]["curl_cffi"] = []
                for u in test_urls:
                    r = cr.get(u, headers=hdrs, cookies={"aws-waf-token": captured["waf_token"]},
                               impersonate="chrome", timeout=20)
                    captured["replay"]["curl_cffi"].append({
                        "url": u, "status": r.status_code, "len": len(r.text),
                        "sample": r.text[:600],
                    })
                    print(f"  [REPLAY curl] {r.status_code}  {_short(u)}  ({len(r.text)}B)")
            except Exception as e:
                print(f"  [REPLAY curl] ERR {e}")

        with open(OUT_PATH, "w") as f:
            json.dump(captured, f, indent=2, default=str)
        print(f"\nSaved -> {OUT_PATH}")
        ctx.close()


def run_browser_recon(verbose_default=False):
    from playwright.sync_api import sync_playwright

    state = {"verbose": verbose_default}
    current_phase = {"name": "startup", "requests": []}
    all_phases = []
    pending_responses = []

    def _short(url, n=120):
        for p in ("https://", "http://"):
            if url.startswith(p):
                url = url[len(p):]
                break
        return url if len(url) <= n else url[: n - 3] + "..."

    def _is_tracked(url):
        return any(d in url for d in TRACKED_DOMAINS)

    def _fmt(nb):
        return f"{nb}B" if nb < 1024 else f"{nb/1024:.1f}KB"

    def handle_request(request):
        if not _is_tracked(request.url):
            return
        if request.resource_type in ("image", "stylesheet", "font", "media"):
            return
        current_phase["requests"].append({
            "url": request.url,
            "method": request.method,
            "resource_type": request.resource_type,
            "post_data": request.post_data,
        })

    def handle_response(response):
        if not _is_tracked(response.url):
            return
        for entry in current_phase["requests"]:
            if entry["url"] == response.url and "status" not in entry:
                entry["status"] = response.status
                ct = response.headers.get("content-type", "")
                entry["content_type"] = ct
                if "json" in ct or "html" in ct:
                    pending_responses.append((response, entry))
                break

    def handle_websocket(ws):
        rec = {"url": ws.url, "frames": []}
        current_phase.setdefault("websockets", []).append(rec)
        print(f"    [WS OPEN] {_short(ws.url)}")

        def on_frame(payload, direction):
            try:
                text = payload if isinstance(payload, str) else payload.decode("utf-8", "ignore")
            except Exception:
                return
            low = text.lower()
            if any(kw in low for kw in SGP_KEYWORDS):
                rec["frames"].append({"dir": direction, "data": text[:4000]})
                print(f"    [WS {direction}] {_short(ws.url, 60)}  {text[:160]}")

        ws.on("framereceived", lambda p: on_frame(p, "recv"))
        ws.on("framesent", lambda p: on_frame(p, "sent"))

    def flush_pending_responses():
        matches = []
        for response, entry in pending_responses:
            try:
                ct = entry.get("content_type", "")
                body = response.text()
                entry["response_size"] = len(body)
                low = body.lower()
                if "json" in ct:
                    entry["response_preview"] = body[:JSON_PREVIEW_LIMIT]
                    found = [kw for kw in SGP_KEYWORDS if kw in low]
                    if found:
                        entry["keywords_found"] = found
                        matches.append((entry, found))
                elif "html" in ct:
                    snips = []
                    for kw in SGP_KEYWORDS:
                        idx = low.find(kw)
                        if idx >= 0:
                            snips.append({"keyword": kw, "context": body[max(0, idx-80):idx+160]})
                    if snips:
                        entry["snippets"] = snips[:5]
            except Exception:
                pass
        pending_responses.clear()
        for entry, kws in matches:
            print(f"    [MATCH] {entry['method']} {_short(entry['url'])}"
                  f"  ({_fmt(entry.get('response_size', 0))})  kw=[{','.join(sorted(set(kws)))}]")

    def start_phase(name):
        nonlocal current_phase
        if current_phase["requests"] or current_phase.get("websockets"):
            all_phases.append(current_phase)
            save_snapshot()
        current_phase = {"name": name, "requests": []}
        print(f"\n--- Capturing: {name} ---")

    def save_snapshot():
        try:
            with open(OUT_PATH, "w") as f:
                json.dump(all_phases + [current_phase], f, indent=2, default=str)
        except Exception as e:
            print(f"  (save failed: {e})")

    def show_phase_results():
        reqs = current_phase["requests"]
        flagged = [r for r in reqs if "keywords_found" in r or "snippets" in r]
        wsr = current_phase.get("websockets", [])
        ws_hits = sum(len(w["frames"]) for w in wsr)
        print(f"=== Phase '{current_phase['name']}': {len(reqs)} reqs, "
              f"{len(flagged)} flagged, {len(wsr)} WS ({ws_hits} sgp frames) ===")
        for req in flagged:
            kws = req.get("keywords_found") or [s["keyword"] for s in req.get("snippets", [])]
            print(f"  *** [{req['method']}] {_short(req['url'])}  "
                  f"({req.get('status','?')}, {_fmt(req.get('response_size',0))})  "
                  f"kw=[{','.join(sorted(set(kws)))}]")
            if req.get("post_data"):
                p = req["post_data"]
                print(f"        POST: {p[:400]}{'...' if len(p) > 400 else ''}")
            if "response_preview" in req:
                print(f"        body: {req['response_preview'][:400]}")
        if not flagged and not ws_hits:
            print("  (nothing flagged this phase)")

    def prompt(label):
        try:
            raw = input(f"\n[ENTER=next  s=skip  v=verbose  q=quit]  >>> {label}  ").strip().lower()
        except EOFError:
            return False
        if raw in ("q", "quit"):
            return False
        if raw in ("v", "verbose"):
            state["verbose"] = not state["verbose"]
        return True

    with sync_playwright() as pw:
        context = pw.chromium.launch_persistent_context(
            user_data_dir=PROFILE_DIR,
            channel="chrome",
            headless=False,
            viewport={"width": 1400, "height": 900},
            args=["--disable-blink-features=AutomationControlled"],
        )
        page = context.pages[0] if context.pages else context.new_page()
        page.on("request", handle_request)
        page.on("response", handle_response)
        page.on("websocket", handle_websocket)

        start_phase("page_load")
        print(f"Navigating to {CAESARS_URL} ...")
        try:
            page.goto(CAESARS_URL, wait_until="commit", timeout=120000)
        except Exception as e:
            print(f"  Navigation slow ({e}), continuing...")
        page.wait_for_timeout(TAIL_GRACE_MS)

        print("\n" + "=" * 60)
        print("STEP 1: Caesars runs an AWS WAF challenge + GeoComply PLC.")
        print("        Allow location if prompted, log in if needed, then")
        print("        navigate to MLB and let the games load.")
        print("=" * 60)
        cont = prompt("MLB games visible")
        page.wait_for_timeout(TAIL_GRACE_MS)
        flush_pending_responses()
        show_phase_results()

        if cont:
            start_phase("open_game")
            print("\nSTEP 2: Click into one MLB game.")
            cont = prompt("game page loaded")
            page.wait_for_timeout(TAIL_GRACE_MS)
            flush_pending_responses()
            show_phase_results()

        if cont:
            start_phase("open_sgp")
            print("\nSTEP 3: Open the Same Game Parlay / Parlay Builder.")
            cont = prompt("SGP builder open")
            page.wait_for_timeout(TAIL_GRACE_MS)
            flush_pending_responses()
            show_phase_results()

        if cont:
            start_phase("add_runline")
            print("\nSTEP 4: Add a RUN LINE (spread) leg, e.g. Team -1.5.")
            cont = prompt("run-line leg on slip")
            page.wait_for_timeout(TAIL_GRACE_MS)
            flush_pending_responses()
            show_phase_results()

        if cont:
            start_phase("add_total")
            print("\nSTEP 5: Add a TOTAL (Over/Under) leg, same game.")
            print("        >>> The SGP combined-odds push fires here.")
            cont = prompt("total leg on slip")
            page.wait_for_timeout(TAIL_GRACE_MS)
            flush_pending_responses()
            show_phase_results()

        if cont:
            start_phase("inspect_slip")
            print("\nInspecting betslip DOM for combined SGP odds + selection IDs...")
            try:
                info = page.evaluate(r"""() => {
                    const out = {found:false, text:'', american:'', decimal:'', selectionIds:[]};
                    const sels = ['[class*="betslip" i]','[id*="betslip" i]',
                                  '[class*="parlay" i]','[class*="slip" i]',
                                  '[class*="sgp" i]','[class*="quickbet" i]'];
                    for (const s of sels) {
                        const el = document.querySelector(s);
                        if (el && el.innerText.trim().length > 5) {
                            out.found = true; out.selector = s;
                            out.text = el.innerText.substring(0, 3000);
                            const a = el.innerText.match(/[+-]\d{3,}/g); if (a) out.american = a.join(', ');
                            const d = el.innerText.match(/\b\d+\.\d{2}\b/g); if (d) out.decimal = d.join(', ');
                            const html = el.innerHTML;
                            const re = /([a-f0-9]{8}-[a-f0-9]{4}-[a-f0-9]{4}-[a-f0-9]{4}-[a-f0-9]{12})/gi;
                            let m; while ((m = re.exec(html)) !== null) out.selectionIds.push(m[1]);
                            break;
                        }
                    }
                    const u = new URL(window.location.href);
                    const sid = u.searchParams.get('selectionIds');
                    if (sid) out.urlSelectionIds = sid;
                    return out;
                }""")
                if info.get("found"):
                    print(f"  slip selector: {info.get('selector')}")
                    if info.get("american"): print(f"  american odds: {info['american']}")
                    if info.get("decimal"):  print(f"  decimal odds:  {info['decimal']}")
                    if info.get("selectionIds"): print(f"  selectionIds (dom): {info['selectionIds']}")
                    if info.get("urlSelectionIds"): print(f"  selectionIds (url): {info['urlSelectionIds']}")
                    print(f"\n  slip text:\n{info['text'][:1500]}")
                    current_phase["dom_slip"] = info
                else:
                    print("  No slip element found in DOM.")
            except Exception as e:
                print(f"  DOM inspect failed: {e}")
            page.wait_for_timeout(TAIL_GRACE_MS)
            flush_pending_responses()
            show_phase_results()

        all_phases.append(current_phase)
        save_snapshot()
        print(f"\nSaved capture to {OUT_PATH}")
        print("Press ENTER to close the browser.")
        try:
            input()
        except EOFError:
            pass
        context.close()


if __name__ == "__main__":
    ap = argparse.ArgumentParser()
    ap.add_argument("--rest", action="store_true",
                    help="just re-confirm the AWS WAF 403 wall (no browser)")
    ap.add_argument("--auto", action="store_true",
                    help="headless auto recon: mint WAF token + replay REST")
    ap.add_argument("--verbose", action="store_true")
    args = ap.parse_args()
    if args.rest:
        probe_rest()
    elif args.auto:
        run_auto_recon(verbose_default=args.verbose)
    else:
        run_browser_recon(verbose_default=args.verbose)
