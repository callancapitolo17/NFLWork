#!/usr/bin/env python3
"""
ProphetX SGP Network Recon Tool

Captures network traffic while you manually explore ProphetX's MLB markets
and (if available) the Parlay Builder. Mirrors the DraftKings recon pattern,
but with faster interaction: live keyword alerts, compact phase summaries,
incremental JSON saves, and skippable phases.

Public intel going in (2026-04):
    - ProphetX is a P2P exchange; singles are P2P matched.
    - Parlay Builder is RFQ-based (third-party MMs quote) and was last
      reported as mobile-only / beta. We still run a web capture because
      product status changes; if no SGP surface appears, capture singles.
    - Nuxt bundle hints at endpoints: /api/v3/events/,
      /api/v2/public/get_multiple_markets, /api/v1/market-lines/,
      /api/v1/markets/moneyline/odds-ladder, /v1/auth/login.
    - GeoComply geofence on API. User must be in a legal state (FL ok).

Usage:
    cd mlb_sgp
    python recon_prophetx_sgp.py            # compact output
    python recon_prophetx_sgp.py --verbose  # dump every request

At any prompt, type:
    [ENTER]  advance to next phase
    s        skip current phase
    v        toggle verbose for the next summary
    q        quit and save what we have so far
"""

from playwright.sync_api import sync_playwright
import argparse
import os
import json
import time

# ---------------------------------------------------------------------------
# Config
# ---------------------------------------------------------------------------
PROPHETX_URL = "https://www.prophetx.co"
PROFILE_DIR = os.path.join(os.path.dirname(__file__), ".prophetx_profile")
OUT_PATH = os.path.join(os.path.dirname(__file__), "recon_prophetx_sgp.json")

# Domains we care about (everything else gets filtered out as noise)
TRACKED_DOMAINS = (
    "prophetx.co",
    "betprophet.co",
    "prophet.co",
)

# Keywords that signal we found the SGP / parlay / RFQ pricing endpoint
SGP_KEYWORDS = [
    "sgp", "parlay", "rfq", "quote", "combined", "combinedodds",
    "correlat", "selection", "betslip", "wager", "payout",
    "market_line", "marketline", "sgm", "builder",
]

# How long after ENTER to wait for tail-end responses to arrive
TAIL_GRACE_MS = 300

# Max bytes of a JSON response body to store in the capture log
JSON_PREVIEW_LIMIT = 10_000


def _format_size(n_bytes):
    if n_bytes < 1024:
        return f"{n_bytes}B"
    return f"{n_bytes / 1024:.1f}KB"


def _short_url(url, max_len=120):
    """Drop scheme+host so the path is readable."""
    for prefix in ("https://", "http://"):
        if url.startswith(prefix):
            url = url[len(prefix):]
            break
    if len(url) > max_len:
        url = url[: max_len - 3] + "..."
    return url


def _is_tracked(url):
    return any(d in url for d in TRACKED_DOMAINS)


def hotkey_input(prompt):
    """Returns one of: 'next', 'skip', 'verbose', 'quit'."""
    raw = input(prompt).strip().lower()
    if raw in ("q", "quit", "exit"):
        return "quit"
    if raw in ("s", "skip"):
        return "skip"
    if raw in ("v", "verbose"):
        return "verbose"
    return "next"


def run_recon(verbose_default=False):
    state = {"verbose": verbose_default}
    current_phase = {"name": "startup", "requests": []}
    all_phases = []
    pending_responses = []

    # --- Playwright event handlers ---

    def handle_request(request):
        if not _is_tracked(request.url):
            return
        if request.resource_type in ("image", "stylesheet", "font", "media"):
            return
        entry = {
            "url": request.url,
            "method": request.method,
            "resource_type": request.resource_type,
            "post_data": request.post_data,
        }
        current_phase["requests"].append(entry)

    def handle_response(response):
        if not _is_tracked(response.url):
            return
        # Attach status/content-type to the in-flight request entry
        for entry in current_phase["requests"]:
            if entry["url"] == response.url and "status" not in entry:
                entry["status"] = response.status
                content_type = response.headers.get("content-type", "")
                entry["content_type"] = content_type
                if "json" in content_type or "html" in content_type:
                    pending_responses.append((response, entry))
                break

    def flush_pending_responses():
        """Drain queued responses, annotate keywords, emit live alerts."""
        matches = []
        for response, entry in pending_responses:
            try:
                content_type = entry.get("content_type", "")
                body = response.text()

                if "json" in content_type:
                    entry["response_preview"] = body[:JSON_PREVIEW_LIMIT]
                    entry["response_size"] = len(body)
                    body_lower = body.lower()
                    found = [kw for kw in SGP_KEYWORDS if kw in body_lower]
                    if found:
                        entry["keywords_found"] = found
                        matches.append((entry, found))

                elif "html" in content_type:
                    entry["response_size"] = len(body)
                    body_lower = body.lower()
                    snippets = []
                    for keyword in SGP_KEYWORDS:
                        idx = body_lower.find(keyword)
                        if idx >= 0:
                            snippets.append({
                                "keyword": keyword,
                                "context": body[max(0, idx - 100):idx + 200],
                            })
                    if snippets:
                        entry["snippets"] = snippets[:5]
            except Exception:
                pass
        pending_responses.clear()

        # Live alerts — make SGP-flagged requests impossible to miss
        for entry, kws in matches:
            print(f"    [MATCH] {entry['method']} {_short_url(entry['url'])}"
                  f"  ({_format_size(entry.get('response_size', 0))})"
                  f"  kw=[{','.join(sorted(set(kws)))}]")

    # --- Phase helpers ---

    def start_phase(name):
        nonlocal current_phase
        if current_phase["requests"]:
            all_phases.append(current_phase)
            save_snapshot()
        current_phase = {"name": name, "requests": []}
        print(f"\n--- Capturing: {name} ---")

    def save_snapshot():
        """Write incremental JSON so we never lose data mid-session."""
        try:
            with open(OUT_PATH, "w") as f:
                json.dump(all_phases + [current_phase], f, indent=2, default=str)
        except Exception as e:
            print(f"  (save snapshot failed: {e})")

    def show_phase_results():
        reqs = current_phase["requests"]
        flagged = [r for r in reqs if "keywords_found" in r or "snippets" in r]
        print(f"=== Phase '{current_phase['name']}': "
              f"{len(reqs)} requests, {len(flagged)} flagged ===")
        if not reqs:
            print("  (no requests captured)")
            return

        # Always show flagged first (these are the signal)
        for req in flagged:
            method = req["method"]
            url = _short_url(req["url"])
            status = req.get("status", "?")
            size = _format_size(req.get("response_size", 0))
            kws = req.get("keywords_found") or [s["keyword"] for s in req.get("snippets", [])]
            print(f"  *** [{method}] {url}  ({status}, {size})  kw=[{','.join(sorted(set(kws)))}]")
            if req.get("post_data"):
                post = req["post_data"]
                preview = post[:400] + ("..." if len(post) > 400 else "")
                print(f"        POST: {preview}")
            if "response_preview" in req:
                print(f"        body: {req['response_preview'][:400]}")

        if state["verbose"]:
            print("  (verbose: all unflagged requests)")
            for req in reqs:
                if req in flagged:
                    continue
                print(f"    [{req['method']}] {_short_url(req['url'])}  "
                      f"({req.get('status', '?')}, "
                      f"{_format_size(req.get('response_size', 0))})")
        else:
            unflagged = len(reqs) - len(flagged)
            if unflagged:
                print(f"  ({unflagged} unflagged requests hidden — run with --verbose or type 'v')")

    def prompt_phase(label):
        """Prompt the user; return True to continue, False to quit."""
        while True:
            action = hotkey_input(f"\n[ENTER=next  s=skip  v=verbose  q=quit]  >>> {label}  ")
            if action == "quit":
                return False
            if action == "verbose":
                state["verbose"] = not state["verbose"]
                print(f"  verbose -> {state['verbose']}")
                continue
            if action == "skip":
                print("  (skipping phase — partial capture kept)")
            return True

    # ------------------------------------------------------------------
    # Main browser session
    # ------------------------------------------------------------------
    with sync_playwright() as p:
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

        # ---- Phase 1: page_load + login ----
        start_phase("page_load")
        print(f"Navigating to {PROPHETX_URL} ...")
        try:
            page.goto(PROPHETX_URL, wait_until="commit", timeout=120000)
        except Exception as e:
            print(f"  Navigation slow ({e}), continuing anyway...")
        page.wait_for_timeout(TAIL_GRACE_MS)

        print("\n" + "=" * 60)
        print("STEP 1: If you're not logged in, log in now.")
        print("        GeoComply will check geolocation — allow if prompted.")
        print("        Then navigate to the MLB section.")
        print("Press ENTER when you can see MLB games.")
        print("=" * 60)
        if not prompt_phase("logged in + MLB lobby visible"):
            goto_save_and_exit = True
        else:
            goto_save_and_exit = False
        page.wait_for_timeout(TAIL_GRACE_MS)
        flush_pending_responses()
        show_phase_results()

        # ---- Phase 2: open_game ----
        if not goto_save_and_exit:
            start_phase("open_game")
            print("\n" + "=" * 60)
            print("STEP 2: Click into a specific MLB game (one happening soon).")
            print("        We want to see what endpoints fire for game markets.")
            print("=" * 60)
            if not prompt_phase("game page is loaded"):
                goto_save_and_exit = True
            page.wait_for_timeout(TAIL_GRACE_MS)
            flush_pending_responses()
            show_phase_results()

        # ---- Phase 3: try_parlay_builder ----
        if not goto_save_and_exit:
            start_phase("try_parlay_builder")
            print("\n" + "=" * 60)
            print("STEP 3: Look for a 'Parlay Builder', 'SGP', 'Same Game")
            print("        Parlay', or similar. It may be under a tab, a")
            print("        toggle, or a button on the game page.")
            print("        >>> If you CAN'T find one on the web: type 's' to")
            print("            skip. That's a crucial data point on its own.")
            print("=" * 60)
            if not prompt_phase("parlay builder opened (or confirmed missing)"):
                goto_save_and_exit = True
            page.wait_for_timeout(TAIL_GRACE_MS)
            flush_pending_responses()
            show_phase_results()

        # ---- Phase 4: add_leg1 ----
        if not goto_save_and_exit:
            start_phase("add_leg1")
            print("\n" + "=" * 60)
            print("STEP 4: Add ONE leg to the slip.")
            print("        Prefer a spread/runline or total so we can see")
            print("        whether the RFQ quote fires for correlated markets.")
            print("        If no Parlay Builder exists, add ANY single bet")
            print("        to the slip — we'll capture the pricing endpoint.")
            print("=" * 60)
            if not prompt_phase("first leg is on the slip"):
                goto_save_and_exit = True
            page.wait_for_timeout(TAIL_GRACE_MS)
            flush_pending_responses()
            show_phase_results()

        # ---- Phase 5: add_leg2 ----
        if not goto_save_and_exit:
            start_phase("add_leg2")
            print("\n" + "=" * 60)
            print("STEP 5: Add a SECOND leg from the SAME game.")
            print("        >>> If ProphetX supports SGP, adding the 2nd leg")
            print("            should fire an RFQ or quote endpoint.")
            print("        If no SGP: add any second leg anyway and note")
            print("        whether the slip treats it as a parlay or refuses.")
            print("=" * 60)
            if not prompt_phase("second leg is on the slip"):
                goto_save_and_exit = True
            page.wait_for_timeout(TAIL_GRACE_MS)
            flush_pending_responses()
            show_phase_results()

        # ---- Phase 6: inspect_slip ----
        if not goto_save_and_exit:
            start_phase("inspect_slip")
            print("\n" + "=" * 60)
            print("Inspecting bet slip DOM ...")
            print("=" * 60)

            try:
                slip_info = page.evaluate(r"""() => {
                    const result = {found: false, selector: null, text: '',
                                    odds: '', decimalOdds: '', dataAttrs: [], selectionIds: []};
                    const selectors = [
                        '[class*="betslip"]', '[id*="betslip"]',
                        '[class*="BetSlip"]', '[id*="BetSlip"]',
                        '[class*="bet-slip"]',
                        '[class*="parlay"]', '[id*="parlay"]',
                        '[class*="Parlay"]',
                        '[class*="builder"]', '[class*="Builder"]',
                        '[class*="slip"]', '[id*="slip"]',
                        '[class*="wager"]', '[class*="Wager"]',
                        '[class*="quote"]', '[class*="rfq"]',
                    ];
                    for (const sel of selectors) {
                        const el = document.querySelector(sel);
                        if (el && el.innerText.trim().length > 5) {
                            result.found = true;
                            result.selector = sel;
                            result.text = el.innerText.substring(0, 3000);
                            const oddsMatch = el.innerText.match(/[+-]\d{3,}/g);
                            if (oddsMatch) result.odds = oddsMatch.join(', ');
                            const decimalMatch = el.innerText.match(/\b\d+\.\d{2}\b/g);
                            if (decimalMatch) result.decimalOdds = decimalMatch.join(', ');
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
                            result.dataAttrs = result.dataAttrs.slice(0, 50);
                            const idPattern = /selection[_-]?id["\s:=]+["']?([\w-]+)/gi;
                            let m;
                            const html = el.innerHTML;
                            while ((m = idPattern.exec(html)) !== null) {
                                result.selectionIds.push(m[1]);
                            }
                            break;
                        }
                    }
                    return result;
                }""")
                if slip_info["found"]:
                    print(f"  Slip element: {slip_info['selector']}")
                    if slip_info.get("odds"):
                        print(f"  American odds: {slip_info['odds']}")
                    if slip_info.get("decimalOdds"):
                        print(f"  Decimal odds: {slip_info['decimalOdds']}")
                    if slip_info.get("selectionIds"):
                        print(f"  Selection IDs: {slip_info['selectionIds']}")
                    if slip_info.get("dataAttrs"):
                        print(f"  {len(slip_info['dataAttrs'])} data-attrs:")
                        for da in slip_info["dataAttrs"][:10]:
                            print(f"    <{da['tag']}> {da['attr']}=\"{da['value']}\"")
                    print(f"\n  Slip text:\n{slip_info['text'][:1500]}")
                    current_phase["dom_slip"] = slip_info
                else:
                    print("  No slip element found in DOM.")
            except Exception as e:
                print(f"  DOM inspection failed: {e}")

            page.wait_for_timeout(TAIL_GRACE_MS)
            flush_pending_responses()
            show_phase_results()

        # ------------------------------------------------------------------
        # Final cross-phase summary
        # ------------------------------------------------------------------
        all_phases.append(current_phase)

        print("\n" + "=" * 60)
        print("FINAL SUMMARY — unique URLs across all phases")
        print("=" * 60)

        endpoint_phases = {}
        endpoint_details = {}
        for phase in all_phases:
            for req in phase["requests"]:
                url = req["url"]
                if url not in endpoint_phases:
                    endpoint_phases[url] = set()
                    endpoint_details[url] = req
                endpoint_phases[url].add(phase["name"])
                if "response_preview" in req:
                    endpoint_details[url] = req

        # Flagged endpoints first
        flagged_urls = [u for u, d in endpoint_details.items() if "keywords_found" in d or "snippets" in d]
        other_urls = [u for u in endpoint_details if u not in flagged_urls]

        print(f"\n  Flagged endpoints ({len(flagged_urls)}):")
        for url in sorted(flagged_urls):
            d = endpoint_details[url]
            kws = d.get("keywords_found") or [s["keyword"] for s in d.get("snippets", [])]
            phase_list = ", ".join(sorted(endpoint_phases[url]))
            print(f"    *** {d['method']} {_short_url(url)}  [{phase_list}]  kw=[{','.join(sorted(set(kws)))}]")

        if state["verbose"]:
            print(f"\n  Other endpoints ({len(other_urls)}):")
            for url in sorted(other_urls):
                d = endpoint_details[url]
                phase_list = ", ".join(sorted(endpoint_phases[url]))
                print(f"    {d['method']} {_short_url(url)}  ({d.get('status', '?')})  [{phase_list}]")
        else:
            print(f"\n  + {len(other_urls)} unflagged endpoints (run with --verbose to list)")

        save_snapshot()
        print(f"\nSaved all captured data to {OUT_PATH}")
        print("\nDone. Press ENTER to close browser (or 'q').")
        try:
            input()
        except EOFError:
            pass
        context.close()


if __name__ == "__main__":
    ap = argparse.ArgumentParser()
    ap.add_argument("--verbose", action="store_true", help="dump every captured request in phase summaries")
    args = ap.parse_args()
    run_recon(verbose_default=args.verbose)
