#!/usr/bin/env python3
"""
DraftKings SGP Network Recon Tool

Captures network traffic while you manually build a Same Game Parlay on DraftKings.
Goal: discover the API endpoint that calculates combined SGP odds.

Usage:
    cd mlb_sgp
    python recon_draftkings_sgp.py

Steps:
    1. Browser opens to DraftKings MLB page
    2. Press ENTER -> navigate to a game
    3. Press ENTER -> click "Same Game Parlay" tab
    4. Press ENTER -> add a spread/run line leg
    5. Press ENTER -> add a total (over/under) leg
    6. Press ENTER -> inspect the betslip for combined SGP odds
    7. Results saved to recon_dk_sgp.json
"""

from playwright.sync_api import sync_playwright
import os
import json

# ---------------------------------------------------------------------------
# Config
# ---------------------------------------------------------------------------
DK_MLB_URL = "https://sportsbook.draftkings.com/leagues/baseball/mlb"
PROFILE_DIR = os.path.join(os.path.dirname(__file__), ".dk_profile")

# Keywords that signal we found the SGP pricing endpoint
SGP_KEYWORDS = [
    "sgp", "parlay", "payout", "towin", "combinedOdds",
    "correlat", "selection", "betslip", "wager",
]

# Max bytes of a JSON response body to store in the capture log
JSON_PREVIEW_LIMIT = 10_000


def _format_size(n_bytes):
    """Human-readable byte size (e.g. '1.2KB')."""
    if n_bytes < 1024:
        return f"{n_bytes}B"
    return f"{n_bytes / 1024:.1f}KB"


def run_recon():
    # ------------------------------------------------------------------
    # Phase tracking — same pattern as recon_parlay_slip.py.
    # `current_phase` holds the active phase dict; `all_phases` collects
    # completed phases so we can dump everything to JSON at the end.
    # ------------------------------------------------------------------
    current_phase = {"name": "startup", "requests": []}
    all_phases = []

    # --- Playwright event handlers ---

    def handle_request(request):
        """Called for every outgoing network request in the browser."""
        # Only care about DraftKings traffic
        if "draftkings.com" not in request.url:
            return
        # Skip static assets — images, CSS, fonts add noise
        if request.resource_type in ("image", "stylesheet", "font"):
            return

        entry = {
            "url": request.url,
            "method": request.method,
            "resource_type": request.resource_type,
            "post_data": request.post_data,
        }
        current_phase["requests"].append(entry)

    # We can't call response.text() inside the event handler — it deadlocks
    # because Playwright's sync API can't nest sync calls inside callbacks.
    # Instead, store (response, entry) pairs and process bodies after each phase.
    pending_responses = []

    def handle_response(response):
        """Called for every incoming response. Store status and queue body reads."""
        if "draftkings.com" not in response.url:
            return

        for entry in current_phase["requests"]:
            if entry["url"] == response.url and "status" not in entry:
                entry["status"] = response.status
                content_type = response.headers.get("content-type", "")
                entry["content_type"] = content_type
                # Queue for body reading after the phase completes
                if "json" in content_type or "html" in content_type:
                    pending_responses.append((response, entry))
                break

    def flush_pending_responses():
        """Read response bodies outside the event handler (safe to call sync)."""
        for response, entry in pending_responses:
            try:
                content_type = entry.get("content_type", "")
                body = response.text()

                if "json" in content_type:
                    entry["response_preview"] = body[:JSON_PREVIEW_LIMIT]
                    entry["response_size"] = len(body)
                    body_lower = body.lower()
                    found = [kw for kw in SGP_KEYWORDS if kw.lower() in body_lower]
                    if found:
                        entry["keywords_found"] = found

                elif "html" in content_type:
                    entry["response_size"] = len(body)
                    body_lower = body.lower()
                    for keyword in SGP_KEYWORDS:
                        idx = body_lower.find(keyword.lower())
                        if idx >= 0:
                            snippet = body[max(0, idx - 100):idx + 200]
                            if "snippets" not in entry:
                                entry["snippets"] = []
                            entry["snippets"].append({
                                "keyword": keyword,
                                "context": snippet,
                            })
            except Exception:
                pass
        pending_responses.clear()

    # --- Phase helpers ---

    def start_phase(name):
        """Finish the current phase and begin a new one."""
        nonlocal current_phase
        if current_phase["requests"]:
            all_phases.append(current_phase)
        current_phase = {"name": name, "requests": []}
        print(f"\n--- Capturing: {name} ---")

    def show_phase_results():
        """Print a categorised summary of every request captured in the
        current phase.  Format mirrors recon_parlay_slip.py but adds the
        [API] / [XHR] / [PAGE] / [JS] tags from the spec."""
        reqs = current_phase["requests"]
        print(f"\n=== Phase: {current_phase['name']} ({len(reqs)} requests captured) ===")
        if not reqs:
            print("  (no requests captured)")
            return

        for req in reqs:
            rtype = req.get("resource_type", "?")
            status = req.get("status", "?")
            method = req["method"]
            url = req["url"]
            size_str = _format_size(req["response_size"]) if "response_size" in req else ""

            # Categorise for the label
            if rtype in ("xhr", "fetch"):
                if "response_preview" in req:
                    label = "API"
                else:
                    label = "XHR"
            elif rtype == "script":
                label = "JS"
            elif rtype == "document":
                label = "PAGE"
            else:
                label = rtype.upper()

            has_json = "response_preview" in req
            has_keywords = "keywords_found" in req
            has_snippets = "snippets" in req

            # Highlight requests that look SGP-related
            star = " ***" if (has_keywords or has_snippets) else ""
            content_note = f" JSON" if has_json else ""
            print(f"  [{label}]  {method} {url[:140]}  ({status}, {size_str}{content_note}){star}")

            # POST body
            if req.get("post_data"):
                post = req["post_data"]
                preview = post[:500] + ("..." if len(post) > 500 else "")
                print(f"         POST data: {preview}")

            # Keywords found in JSON
            if has_keywords:
                print(f"         Keywords found: {', '.join(req['keywords_found'])}")

            # JSON response preview
            if has_json:
                print(f"         Response preview: {req['response_preview'][:500]}")

            # HTML snippets
            if has_snippets:
                print(f"         SGP-RELATED KEYWORDS IN HTML:")
                for s in req["snippets"][:3]:
                    print(f"           '{s['keyword']}': ...{s['context'][:150]}...")

    # ------------------------------------------------------------------
    # Main browser session
    # ------------------------------------------------------------------
    with sync_playwright() as p:
        # launch_persistent_context keeps cookies/local-storage between runs
        # so DraftKings doesn't hit you with repeated captchas.
        context = p.chromium.launch_persistent_context(
            user_data_dir=PROFILE_DIR,
            channel="chrome",
            headless=False,
            viewport={"width": 1400, "height": 900},
            args=["--disable-blink-features=AutomationControlled"],
        )
        page = context.pages[0] if context.pages else context.new_page()
        page.on("request", handle_request)
        page.on("response", handle_response)

        # ---- Phase 1: page_load ----
        start_phase("page_load")
        print(f"Navigating to {DK_MLB_URL}...")
        # DK's SPA is heavy — don't crash on slow load, just wait
        try:
            page.goto(DK_MLB_URL, wait_until="commit", timeout=120000)
        except Exception as e:
            print(f"  Navigation slow ({e}), continuing anyway...")
        page.wait_for_timeout(5000)

        print("\n" + "=" * 60)
        print("MLB schedule should be visible (or loading).")
        print("Wait until it loads, pick a game, then press ENTER.")
        print("=" * 60)
        input()
        flush_pending_responses()
        show_phase_results()

        # ---- Phase 2: open_game ----
        start_phase("open_game")
        print("=" * 60)
        print("STEP 1: Click into a specific MLB game.")
        print("        (Click the game row / 'More Wagers' link)")
        print("Press ENTER after the game page loads.")
        print("=" * 60)
        input()
        page.wait_for_timeout(2000)
        flush_pending_responses()
        show_phase_results()

        # ---- Phase 3: open_sgp ----
        start_phase("open_sgp")
        print("\n" + "=" * 60)
        print("STEP 2: Click the 'Same Game Parlay' (SGP) tab.")
        print("        It should be near the top of the game page.")
        print("Press ENTER after clicking.")
        print("=" * 60)
        input()
        page.wait_for_timeout(2000)
        flush_pending_responses()
        show_phase_results()

        # ---- Phase 4: add_spread ----
        start_phase("add_spread")
        print("\n" + "=" * 60)
        print("STEP 3: Add a SPREAD / RUN LINE leg.")
        print("        (e.g., click Team -1.5)")
        print("Press ENTER after clicking.")
        print("=" * 60)
        input()
        page.wait_for_timeout(2000)
        flush_pending_responses()
        show_phase_results()

        # ---- Phase 5: add_total ----
        start_phase("add_total")
        print("\n" + "=" * 60)
        print("STEP 4: Add a TOTAL (Over/Under) leg.")
        print("        (e.g., click Over 8.5)")
        print("        >>> This is where the SGP pricing API should fire.")
        print("Press ENTER after clicking.")
        print("=" * 60)
        input()
        page.wait_for_timeout(2000)
        flush_pending_responses()
        show_phase_results()

        # ---- Phase 6: inspect_slip ----
        start_phase("inspect_slip")
        print("\n" + "=" * 60)
        print("Inspecting betslip DOM for SGP odds...")
        print("=" * 60)

        # --- DOM inspection: betslip text, odds, data attributes ---
        try:
            slip_info = page.evaluate(r"""() => {
                const result = {found: false, selector: null, text: '', odds: '', dataAttrs: [], selectionIds: []};

                // Common betslip selectors on DraftKings
                const selectors = [
                    '#betslip', '.betslip', '.bet-slip',
                    '[class*="betslip"]', '[id*="betslip"]',
                    '[class*="BetSlip"]', '[id*="BetSlip"]',
                    '[class*="bet-slip"]',
                    '[class*="parlay"]', '[id*="parlay"]',
                    '[class*="Parlay"]',
                    '[class*="sgp"]', '[id*="sgp"]',
                    '[class*="SGP"]',
                    '[class*="slip"]', '[id*="slip"]',
                    '[class*="Slip"]',
                    '[class*="wager"]',
                ];

                for (const sel of selectors) {
                    const el = document.querySelector(sel);
                    if (el && el.innerText.trim().length > 5) {
                        result.found = true;
                        result.selector = sel;
                        result.text = el.innerText.substring(0, 3000);

                        // Look for American odds pattern (+XXX or -XXX)
                        const oddsMatch = el.innerText.match(/[+-]\\d{3,}/g);
                        if (oddsMatch) result.odds = oddsMatch.join(', ');

                        // Look for decimal odds pattern (e.g. 3.45)
                        const decimalMatch = el.innerText.match(/\\b\\d+\\.\\d{2}\\b/g);
                        if (decimalMatch) result.decimalOdds = decimalMatch.join(', ');

                        // Grab data- attributes from the slip and its children
                        const allEls = el.querySelectorAll('*');
                        for (const child of allEls) {
                            for (const attr of child.attributes) {
                                if (attr.name.startsWith('data-')) {
                                    result.dataAttrs.push({
                                        tag: child.tagName,
                                        attr: attr.name,
                                        value: attr.value.substring(0, 200),
                                    });
                                }
                            }
                        }
                        // Keep data attrs manageable
                        result.dataAttrs = result.dataAttrs.slice(0, 50);

                        // Look for selection IDs
                        const idPattern = /selection[_-]?id["\s:=]+["']?(\w+)/gi;
                        let m;
                        const html = el.innerHTML;
                        while ((m = idPattern.exec(html)) !== null) {
                            result.selectionIds.push(m[1]);
                        }

                        break;
                    }
                }

                // Fallback: search entire body for SGP odds
                if (!result.found) {
                    const body = document.body.innerText;
                    // Look for "SGP" or "Same Game Parlay" nearby odds
                    const sgpIdx = body.toLowerCase().indexOf('same game parlay');
                    if (sgpIdx >= 0) {
                        result.selector = 'body (Same Game Parlay text)';
                        result.text = body.substring(Math.max(0, sgpIdx - 200), sgpIdx + 500);
                        result.found = true;
                    }
                }

                return result;
            }""")

            if slip_info["found"]:
                print(f"  Found betslip element: {slip_info['selector']}")
                if slip_info.get("odds"):
                    print(f"  American odds in slip: {slip_info['odds']}")
                if slip_info.get("decimalOdds"):
                    print(f"  Decimal odds in slip: {slip_info['decimalOdds']}")
                if slip_info.get("selectionIds"):
                    print(f"  Selection IDs: {slip_info['selectionIds']}")
                if slip_info.get("dataAttrs"):
                    print(f"  Data attributes ({len(slip_info['dataAttrs'])}):")
                    for da in slip_info["dataAttrs"][:10]:
                        print(f"    <{da['tag']}> {da['attr']}=\"{da['value']}\"")
                print(f"\n  Betslip text:\n{slip_info['text'][:2000]}")
            else:
                print("  No betslip element found in DOM.")
        except Exception as e:
            print(f"  Error inspecting DOM: {e}")

        # --- Scan inline scripts for SGP-related logic ---
        try:
            js_hits = page.evaluate("""() => {
                const scripts = document.querySelectorAll('script');
                const found = [];
                for (const s of scripts) {
                    const text = s.textContent || '';
                    if (text.match(/sgp|parlay|combinedOdds|correlat|payout|towin/i)) {
                        found.push(text.substring(0, 2000));
                    }
                }
                return found;
            }""")
            if js_hits:
                print(f"\n  Found {len(js_hits)} inline scripts with SGP-related code:")
                for i, js in enumerate(js_hits):
                    print(f"\n  --- Script {i + 1} ---")
                    print(f"  {js[:1000]}")
        except Exception as e:
            print(f"  Error scanning scripts: {e}")

        flush_pending_responses()
        show_phase_results()

        # ------------------------------------------------------------------
        # Final summary across ALL phases
        # ------------------------------------------------------------------
        all_phases.append(current_phase)

        print("\n" + "=" * 60)
        print("FINAL SUMMARY — All unique API endpoints")
        print("=" * 60)

        # Collect unique endpoints and note which phases they appeared in
        endpoint_phases = {}  # url -> set of phase names
        endpoint_details = {}  # url -> latest entry dict
        for phase in all_phases:
            for req in phase["requests"]:
                url = req["url"]
                if url not in endpoint_phases:
                    endpoint_phases[url] = set()
                    endpoint_details[url] = req
                endpoint_phases[url].add(phase["name"])
                # Keep the richest entry (one with response info)
                if "response_preview" in req:
                    endpoint_details[url] = req

        for url, phases in sorted(endpoint_phases.items(), key=lambda x: x[0]):
            detail = endpoint_details[url]
            method = detail.get("method", "?")
            status = detail.get("status", "?")
            size_str = _format_size(detail["response_size"]) if "response_size" in detail else ""
            has_json = "response_preview" in detail
            has_keywords = "keywords_found" in detail

            # Highlight endpoints from the leg-adding phases
            is_leg_phase = bool(phases & {"add_spread", "add_total"})
            marker = " <<<< LEG PHASE" if is_leg_phase else ""
            star = " ***" if has_keywords else ""

            content_note = " JSON" if has_json else ""
            phase_list = ", ".join(sorted(phases))
            print(f"  {method} {url[:140]}  ({status}, {size_str}{content_note}) [{phase_list}]{marker}{star}")

            if has_keywords:
                print(f"       Keywords: {', '.join(detail['keywords_found'])}")

        # ------------------------------------------------------------------
        # Save to JSON
        # ------------------------------------------------------------------
        out_path = os.path.join(os.path.dirname(__file__), "recon_dk_sgp.json")
        with open(out_path, "w") as f:
            # Convert sets to lists for JSON serialization
            json.dump(all_phases, f, indent=2, default=str)
        print(f"\nSaved all captured data to {out_path}")

        print("\n" + "=" * 60)
        print("Done! Press ENTER to close browser.")
        print("=" * 60)
        input()
        context.close()


if __name__ == "__main__":
    run_recon()
