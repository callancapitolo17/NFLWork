#!/usr/bin/env python3
"""
Novig SGP Network Recon Tool

Captures network traffic while you manually explore Novig's MLB markets
and (if available) an SGP product. Mirrors the DraftKings recon pattern
with speed improvements: live keyword alerts, compact phase summaries,
incremental JSON saves, and skippable phases.

Public intel going in (2026-04):
    - Novig is a sweepstakes P2P exchange (Novig Coins / Novig Cash).
    - Per Sportico (Dec 2025), Novig SGP odds are set by a third-party
      provider — NOT peer-matched. So any SGP capture reflects that
      book's price, not an organic exchange quote.
    - MLB SGP availability in 2026 is ambiguous across reviews; confirming
      this one way or the other is the PRIMARY goal of this recon.
    - Web app at app.novig.us; marketing at www.novig.us. Neither exposes
      markets without login. FL-eligible.
    - No community scraping projects or documented API exist publicly,
      so this recon is also the first real look at the network surface.

Usage:
    cd mlb_sgp
    python recon_novig_sgp.py            # compact output
    python recon_novig_sgp.py --verbose  # dump every request

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

# ---------------------------------------------------------------------------
# Config
# ---------------------------------------------------------------------------
NOVIG_URL = "https://novig.com/events"
# Land directly on the markets lobby. .us variants were dead ends (app.novig.us
# is a 2KB SPA shell; www.novig.us is the marketing site). novig.com is the
# real sportsbook; /events is the markets list page.
PROFILE_DIR = os.path.join(os.path.dirname(__file__), ".novig_profile")
OUT_PATH = os.path.join(os.path.dirname(__file__), "recon_novig_sgp.json")

TRACKED_DOMAINS = (
    "novig.com",
    "novig.us",
    "novig.io",
)

SGP_KEYWORDS = [
    "sgp", "parlay", "combined", "combinedodds", "combined_odds",
    "correlat", "selection", "betslip", "slip", "wager",
    "quote", "payout", "builder", "third_party", "third-party",
    "provider", "sgm", "multi", "accumulator",
]

TAIL_GRACE_MS = 300
JSON_PREVIEW_LIMIT = 10_000


def _format_size(n_bytes):
    if n_bytes < 1024:
        return f"{n_bytes}B"
    return f"{n_bytes / 1024:.1f}KB"


def _short_url(url, max_len=120):
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
        for entry in current_phase["requests"]:
            if entry["url"] == response.url and "status" not in entry:
                entry["status"] = response.status
                content_type = response.headers.get("content-type", "")
                entry["content_type"] = content_type
                if "json" in content_type or "html" in content_type:
                    pending_responses.append((response, entry))
                break

    def flush_pending_responses():
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

        for entry, kws in matches:
            print(f"    [MATCH] {entry['method']} {_short_url(entry['url'])}"
                  f"  ({_format_size(entry.get('response_size', 0))})"
                  f"  kw=[{','.join(sorted(set(kws)))}]")

    def start_phase(name):
        nonlocal current_phase
        if current_phase["requests"]:
            all_phases.append(current_phase)
            save_snapshot()
        current_phase = {"name": name, "requests": []}
        print(f"\n--- Capturing: {name} ---")

    def save_snapshot():
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
        print(f"Navigating to {NOVIG_URL} ...")
        try:
            page.goto(NOVIG_URL, wait_until="commit", timeout=120000)
        except Exception as e:
            print(f"  Navigation slow ({e}), continuing anyway...")
        page.wait_for_timeout(TAIL_GRACE_MS)

        print("\n" + "=" * 60)
        print("STEP 1: Log in to Novig if needed. Allow geolocation if")
        print("        prompted. Then navigate to MLB.")
        print("        (If app.novig.us redirects you, try novig.us/app")
        print("         or whatever your mobile app uses.)")
        print("Press ENTER when MLB games are visible.")
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
            print("STEP 2: Click into a specific MLB game.")
            print("        We want to see what endpoints fire for game markets.")
            print("=" * 60)
            if not prompt_phase("game page is loaded"):
                goto_save_and_exit = True
            page.wait_for_timeout(TAIL_GRACE_MS)
            flush_pending_responses()
            show_phase_results()

        # ---- Phase 3: try_sgp ----
        if not goto_save_and_exit:
            start_phase("try_sgp")
            print("\n" + "=" * 60)
            print("STEP 3: Look for Same Game Parlay / SGP / Parlay Builder.")
            print("        Could be a tab, a toggle, or a button.")
            print("        >>> If NOT present on web — type 's' to skip.")
            print("            That confirms MLB SGP is mobile-only.")
            print("=" * 60)
            if not prompt_phase("SGP interface opened (or confirmed missing)"):
                goto_save_and_exit = True
            page.wait_for_timeout(TAIL_GRACE_MS)
            flush_pending_responses()
            show_phase_results()

        # ---- Phase 4: add_leg1 ----
        if not goto_save_and_exit:
            start_phase("add_leg1")
            print("\n" + "=" * 60)
            print("STEP 4: Add ONE leg (spread/runline/total preferred).")
            print("        If no SGP exists, add any single to the slip so")
            print("        we capture the pricing/quote endpoint for singles.")
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
            print("        >>> Adding the 2nd correlated leg is where an SGP")
            print("            quote endpoint should fire (if one exists).")
            print("        Watch the [MATCH] lines below for keyword hits.")
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
                        '[class*="slip"]', '[id*="slip"]',
                        '[class*="Slip"]',
                        '[class*="wager"]', '[class*="Wager"]',
                        '[class*="builder"]', '[class*="Builder"]',
                        '[class*="sgp"]', '[id*="sgp"]',
                        '[class*="SGP"]',
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
