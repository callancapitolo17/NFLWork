#!/usr/bin/env python3
"""Quick recon — captures ALL DK network requests as you click."""
import json
from playwright.sync_api import sync_playwright

captures = []

with sync_playwright() as p:
    browser = p.chromium.connect_over_cdp("http://localhost:9222")
    ctx = browser.contexts[0]
    pages = ctx.pages
    print(f"Found {len(pages)} pages")

    # Can't read response body in callback (deadlocks in sync API)
    # Instead capture URL + status, then we'll analyze
    def on_req(request):
        url = request.url
        if "draftkings" in url and request.resource_type in ("xhr", "fetch"):
            method = request.method
            post = request.post_data
            print(f"  [{method}] {url[:120]}")
            if post:
                print(f"    POST: {post[:300]}")
            captures.append({"url": url, "method": method, "post_data": post})

    for pg in pages:
        pg.on("request", on_req)

    print("\nClick SGP legs in Chrome. Press ENTER when done.\n")
    input()
    browser.close()

with open("quick_recon.json", "w") as f:
    json.dump(captures, f, indent=2)
print(f"\nSaved {len(captures)} captures to quick_recon.json")
