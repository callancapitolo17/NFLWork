#!/usr/bin/env python3
"""Fetch Torvik ratings via Playwright (bypasses Cloudflare).
Outputs JSON to stdout: [[team, adj_o, adj_d, barthag, record], ...]
Called by shared.R when cbbdata doesn't have current season data.
"""
import sys
import json

try:
    from playwright.sync_api import sync_playwright
except ImportError:
    print("[]")
    sys.exit(0)

year = int(sys.argv[1]) if len(sys.argv) > 1 else 2026

try:
    with sync_playwright() as p:
        browser = p.chromium.launch(headless=True)
        context = browser.new_context(
            user_agent="Mozilla/5.0 (Macintosh; Intel Mac OS X 10_15_7) AppleWebKit/537.36"
        )
        page = context.new_page()
        page.add_init_script('Object.defineProperty(navigator, "webdriver", {get: () => undefined})')

        page.goto(
            f"https://barttorvik.com/trank.php?year={year}&t=0&json=1",
            timeout=30000,
            wait_until="networkidle",
        )
        page.wait_for_timeout(5000)

        text = page.inner_text("body").strip()
        browser.close()

        if text.startswith("["):
            data = json.loads(text)
            # Output simplified: [team, adj_o, adj_d, barthag, record]
            out = []
            for row in data:
                if len(row) >= 5:
                    out.append([row[0], row[1], row[2], row[3], row[4]])
            print(json.dumps(out))
        else:
            print("[]")
except Exception:
    print("[]")
