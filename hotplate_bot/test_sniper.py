#!/usr/bin/env python3
"""
Regression tests for the in-page sniper, run against a mock drop page.

The real Hotplate page cannot be used as a test fixture (one drop, one shot), so
these reproduce the only behaviour that matters: buttons that sit disabled through
a waiting room and then enable all at once. Every test asserts on latency measured
INSIDE the page (sniper clickedAt minus the page's own openedAt), which is the
number the race is actually decided by.

Usage:  python test_sniper.py
Side effects: none. Launches a headless browser against in-memory HTML, touches
no network and no files.
"""

from __future__ import annotations

import sys

from playwright.sync_api import Page, sync_playwright

from hotplate_snipe import BOT_DIR, IN_PAGE_SNIPER_JS, SNIPE_POLL_BACKSTOP_MS, launch_browser

# Generous enough not to flake on a loaded CI box, tight enough that a regression
# to Python-side polling (tens of ms, plus CDP round-trips) would fail the test.
MAX_ACCEPTABLE_CLICK_LATENCY_MS = 60

MOCK_PAGE_TEMPLATE = """
<html><body>
  <h1>Butter and Crumble</h1>
  <div id="menu">
    <div class="card"><h3>Pistachio Croissant</h3><p>$8</p><button disabled>Add to Cart</button></div>
    <div class="card"><h3>Cinnamon Roll</h3><p>$6</p><button disabled>Add to Cart</button></div>
    <div class="card"><h3>Sourdough Loaf</h3><p>$12</p><button disabled>Add to Cart</button></div>
  </div>
  <script>
    window.__clicks = [];
    window.__openedAt = null;

    function wireClickRecorders() {
      document.querySelectorAll('#menu button').forEach(function (button) {
        button.addEventListener('click', function () {
          window.__clicks.push(button.closest('.card').querySelector('h3').textContent);
        });
      });
    }
    wireClickRecorders();

    setTimeout(function () {
      __OPEN_BEHAVIOUR__
      window.__openedAt = Date.now();
    }, __OPEN_DELAY_MS__);
  </script>
</body></html>
"""

# The ordinary case: the same button elements simply stop being disabled.
ENABLE_IN_PLACE = "document.querySelectorAll('#menu button').forEach(function (b) { b.disabled = false; });"

# The nasty case: a framework re-render swaps in brand new nodes at open time,
# detaching anything the sniper cached during the waiting room.
REPLACE_WHOLE_MENU = """
    document.getElementById('menu').innerHTML =
      '<div class="card"><h3>Pistachio Croissant</h3><p>$8</p><button>Add to Cart</button></div>' +
      '<div class="card"><h3>Cinnamon Roll</h3><p>$6</p><button>Add to Cart</button></div>';
    wireClickRecorders();
"""


def build_mock_page(open_delay_ms: int, open_behaviour: str) -> str:
    return MOCK_PAGE_TEMPLATE.replace("__OPEN_DELAY_MS__", str(open_delay_ms)).replace(
        "__OPEN_BEHAVIOUR__", open_behaviour
    )


def run_sniper(page: Page, target_items: list[str], max_wait_ms: int = 15000) -> dict:
    return page.evaluate(
        IN_PAGE_SNIPER_JS,
        {
            "itemNames": target_items,
            "buttonPattern": "add to cart|add|order now|reserve",
            "maxWaitMs": max_wait_ms,
            "backstopMs": SNIPE_POLL_BACKSTOP_MS,
        },
    )


def check(condition: bool, description: str) -> bool:
    print(f"  {'PASS' if condition else 'FAIL'}  {description}")
    return condition


def test_clicks_promptly_when_enabled_in_place(page: Page) -> bool:
    page.set_content(build_mock_page(1200, ENABLE_IN_PLACE))
    result = run_sniper(page, ["Pistachio Croissant", "Cinnamon Roll"])
    opened_at = page.evaluate("window.__openedAt")
    clicks = page.evaluate("window.__clicks")
    latency_ms = result["clickedAt"] - opened_at

    print(f"\n[enable-in-place]  latency={latency_ms}ms  clicks={clicks}")
    return all(
        [
            check(result["status"] == "clicked", "sniper reports a click"),
            check(clicks == ["Pistachio Croissant"], "clicked exactly the top-priority item, once"),
            check(
                0 <= latency_ms <= MAX_ACCEPTABLE_CLICK_LATENCY_MS,
                f"clicked within {MAX_ACCEPTABLE_CLICK_LATENCY_MS}ms of open (got {latency_ms}ms)",
            ),
        ]
    )


def test_survives_full_rerender_at_open(page: Page) -> bool:
    page.set_content(build_mock_page(1200, REPLACE_WHOLE_MENU))
    result = run_sniper(page, ["Pistachio Croissant", "Cinnamon Roll"])
    opened_at = page.evaluate("window.__openedAt")
    clicks = page.evaluate("window.__clicks")
    latency_ms = result["clickedAt"] - opened_at

    print(f"\n[rerender-at-open]  latency={latency_ms}ms  clicks={clicks}")
    return all(
        [
            check(result["status"] == "clicked", "sniper recovers from a detached cache and clicks"),
            check(clicks == ["Pistachio Croissant"], "clicked exactly the top-priority item, once"),
            check(
                0 <= latency_ms <= MAX_ACCEPTABLE_CLICK_LATENCY_MS,
                f"clicked within {MAX_ACCEPTABLE_CLICK_LATENCY_MS}ms of open (got {latency_ms}ms)",
            ),
        ]
    )


def test_falls_through_to_second_choice(page: Page) -> bool:
    """Top choice sold out (button never appears) — the sniper should take number two."""
    page.set_content(build_mock_page(1200, ENABLE_IN_PLACE))
    page.evaluate(
        "document.querySelectorAll('.card')[0].querySelector('button').remove();"
    )
    result = run_sniper(page, ["Pistachio Croissant", "Cinnamon Roll"])
    clicks = page.evaluate("window.__clicks")

    print(f"\n[top-choice-sold-out]  clicks={clicks}")
    return all(
        [
            check(result["status"] == "clicked", "sniper still clicks"),
            check(clicks == ["Cinnamon Roll"], "fell through to the second-priority item"),
        ]
    )


def test_does_not_click_while_disabled(page: Page) -> bool:
    """The waiting room must be respected — an early click is a wasted entry."""
    page.set_content(build_mock_page(999999, ENABLE_IN_PLACE))
    result = run_sniper(page, ["Pistachio Croissant"], max_wait_ms=1500)
    clicks = page.evaluate("window.__clicks")

    print(f"\n[never-opens]  status={result['status']}  clicks={clicks}")
    return all(
        [
            check(result["status"] == "timeout", "times out rather than clicking a disabled button"),
            check(clicks == [], "no click was fired while buttons were disabled"),
        ]
    )


def test_ignores_items_not_in_config(page: Page) -> bool:
    page.set_content(build_mock_page(800, ENABLE_IN_PLACE))
    run_sniper(page, ["Sourdough Loaf"])
    clicks = page.evaluate("window.__clicks")

    print(f"\n[item-targeting]  clicks={clicks}")
    return check(clicks == ["Sourdough Loaf"], "clicked only the configured item, not its neighbours")


def test_console_recon_pairs_items_to_buttons(page: Page) -> bool:
    """
    console_recon.js walks button -> item name; the sniper walks item name -> button.

    They must agree, because recon's whole job is telling you what the sniper will do.
    """
    snippet = (BOT_DIR / "console_recon.js").read_text()

    page.set_content(build_mock_page(999999, ENABLE_IN_PLACE))
    locked = page.evaluate(snippet)
    pairs_locked = {entry["itemName"]: entry["clickableNow"] for entry in locked["items"]}

    page.evaluate("document.querySelectorAll('#menu button').forEach(function (b) { b.disabled = false; });")
    opened = page.evaluate(snippet)
    pairs_open = {entry["itemName"]: entry["clickableNow"] for entry in opened["items"]}

    print(f"\n[console-recon]  locked={pairs_locked}  open={pairs_open}")
    return all(
        [
            check(
                set(pairs_locked) == {"Pistachio Croissant", "Cinnamon Roll", "Sourdough Loaf"},
                "paired every button back to its own item name",
            ),
            check(not any(pairs_locked.values()), "reports nothing clickable during the waiting room"),
            check(all(pairs_open.values()), "reports buttons clickable once the drop opens"),
        ]
    )


def main() -> int:
    tests = [
        test_clicks_promptly_when_enabled_in_place,
        test_survives_full_rerender_at_open,
        test_falls_through_to_second_choice,
        test_does_not_click_while_disabled,
        test_ignores_items_not_in_config,
        test_console_recon_pairs_items_to_buttons,
    ]

    with sync_playwright() as playwright:
        browser = launch_browser(playwright, headless=True)
        page = browser.new_page()
        results = [test(page) for test in tests]
        browser.close()

    passed = sum(1 for result in results if result)
    print(f"\n{'=' * 60}\n{passed}/{len(results)} tests passed\n{'=' * 60}")
    return 0 if passed == len(results) else 1


if __name__ == "__main__":
    sys.exit(main())
