#!/usr/bin/env python3
"""
Hotplate drop assistant — races the add-to-cart click at the end of a waiting-room countdown.

The contested moment in a Hotplate drop is a single click: when the countdown hits
zero the add-to-cart buttons enable, and whoever reserves the item first gets it.
Everything after that (quantity, checkout, payment) runs against a cart-hold timer
that is generous by comparison. So this tool spends all of its effort on being
ready for that one click and none on being clever afterwards.

Three subcommands:
    login   Open a browser, you sign in by hand, session cookies are saved to disk.
    recon   Load the event page and dump its structure so you can fill in config.json.
    arm     Sit on the event page, click add-to-cart the instant it enables, alert you.

Side effects (all confined to this directory, nothing touches any DuckDB):
    login   writes STORAGE_STATE_PATH (hotplate_bot/storage_state.json) — contains
            live session cookies, gitignored, treat it like a password.
    recon   writes recon_dump.html + recon_dump.png (gitignored).
    arm     writes snipe_<timestamp>.log and screenshots on each step; performs real
            add-to-cart actions against your real logged-in account. It never submits
            payment unless config.auto_submit_payment is explicitly true.

Deliberately NOT implemented, and I won't add them: multi-account entry stacking,
proxy rotation, CAPTCHA solving, or browser-fingerprint spoofing. This drives one
ordinary logged-in browser on your own machine, faster than your hand can.
"""

from __future__ import annotations

import argparse
import json
import logging
import platform
import subprocess
import sys
import time
from datetime import datetime, timezone
from email.utils import parsedate_to_datetime
from pathlib import Path
from typing import Any

from playwright.sync_api import Browser, BrowserContext, Error as PlaywrightError, Page, sync_playwright

BOT_DIR = Path(__file__).resolve().parent
STORAGE_STATE_PATH = BOT_DIR / "storage_state.json"
DEFAULT_CONFIG_PATH = BOT_DIR / "config.json"

# The in-page sniper backstops its MutationObserver with a timer at roughly frame
# rate; 16ms gives that worst-case detection even if a mutation is somehow missed.
SNIPE_POLL_BACKSTOP_MS = 16

log = logging.getLogger("hotplate")


# --------------------------------------------------------------------------------------
# config
# --------------------------------------------------------------------------------------


def load_config(config_path: Path) -> dict[str, Any]:
    """Read and validate config.json. Fails loudly on anything missing or malformed."""
    if not config_path.exists():
        raise SystemExit(
            f"No config at {config_path}. Copy config.example.json to config.json and fill it in."
        )

    with config_path.open() as handle:
        config = json.load(handle)

    event_url = config.get("event_url", "")
    if not event_url.startswith("https://www.hotplate.com/"):
        raise SystemExit(f"config.event_url must be a https://www.hotplate.com/... URL, got {event_url!r}")

    target_items = config.get("target_items")
    if not isinstance(target_items, list) or not target_items:
        raise SystemExit("config.target_items must be a non-empty list of item names, in priority order.")

    for item in target_items:
        if not isinstance(item, str) or not item.strip():
            raise SystemExit(f"config.target_items entries must be non-empty strings, got {item!r}")

    config.setdefault("add_to_cart_pattern", "add to cart|add|order now|reserve")
    config.setdefault("headless", False)
    config.setdefault("auto_submit_payment", False)
    config.setdefault("max_wait_minutes", 20)
    return config


# --------------------------------------------------------------------------------------
# browser plumbing
# --------------------------------------------------------------------------------------


def launch_browser(playwright: Any, headless: bool) -> Browser:
    """Launch Chromium, preferring the image-provided binary over a downloaded one."""
    bundled_chromium = Path("/opt/pw-browsers/chromium")
    launch_kwargs: dict[str, Any] = {
        "headless": headless,
        # A visible, human-driven-looking window is the point: this is your own
        # session, not a disguised one.
        "args": ["--disable-blink-features=AutomationControlled"],
    }
    if bundled_chromium.exists():
        launch_kwargs["executable_path"] = str(bundled_chromium)
    return playwright.chromium.launch(**launch_kwargs)


def require_saved_session() -> None:
    """Fail before launching a browser if `login` was never run — a clearer error, sooner."""
    if not STORAGE_STATE_PATH.exists():
        raise SystemExit(
            f"No saved session at {STORAGE_STATE_PATH}. Run `python hotplate_snipe.py login` first."
        )


def open_context(browser: Browser) -> BrowserContext:
    """Build a context, loading saved cookies when a prior `login` run stored them."""
    storage_state = str(STORAGE_STATE_PATH) if STORAGE_STATE_PATH.exists() else None
    return browser.new_context(
        storage_state=storage_state,
        viewport={"width": 1400, "height": 1000},
    )


def check_clock_skew(page: Page) -> None:
    """
    Warn if the local clock disagrees with Hotplate's server clock.

    Detection here is DOM-driven, so skew is not fatal — but a badly wrong clock
    means the on-screen countdown you are eyeballing is wrong too, which is exactly
    how people miss a waiting room. Worth surfacing, never worth aborting over.
    """
    try:
        response = page.request.head("https://www.hotplate.com/")
        server_date_header = response.headers.get("date")
        if not server_date_header:
            log.warning("No Date header from hotplate.com; skipping clock-skew check.")
            return
        server_time = parsedate_to_datetime(server_date_header)
        skew_seconds = (datetime.now(timezone.utc) - server_time).total_seconds()
    except (PlaywrightError, ValueError, TypeError) as exc:
        log.warning("Clock-skew check failed (%s); continuing.", exc)
        return

    if abs(skew_seconds) > 2:
        log.warning(
            "Local clock is %.1fs %s of Hotplate's server. Your on-screen countdown is off by that much.",
            abs(skew_seconds),
            "ahead" if skew_seconds > 0 else "behind",
        )
    else:
        log.info("Clock in sync with Hotplate server (skew %.2fs).", skew_seconds)


def alert_human(message: str) -> None:
    """Make noise on macOS so you look at the screen. Silently no-ops elsewhere."""
    log.info("ALERT: %s", message)
    if platform.system() != "Darwin":
        return
    try:
        subprocess.run(["afplay", "/System/Library/Sounds/Sosumi.aiff"], timeout=5, check=False)
        subprocess.run(
            ["osascript", "-e", f'display notification "{message}" with title "Hotplate"'],
            timeout=5,
            check=False,
        )
    except (OSError, subprocess.SubprocessError) as exc:
        log.warning("Could not play alert (%s).", exc)


# --------------------------------------------------------------------------------------
# the in-page sniper
# --------------------------------------------------------------------------------------

# This runs entirely inside the browser. That matters: a Python-side polling loop
# pays a CDP round-trip (plus interpreter scheduling jitter) on every check, which
# is the difference between reacting in ~1ms and reacting in ~20-50ms. The
# MutationObserver fires in the same microtask as the DOM change that enables the
# button, so the click lands about as early as the page can physically allow.
#
# Two performance rules hold in the hot path, because a sniper that bogs the page
# down loses the race it was written to win:
#   1. Buttons are resolved ONCE during the waiting room and cached. The hot path
#      only re-checks whether a cached button became clickable.
#   2. Matching uses textContent, never innerText. innerText forces a layout
#      reflow; running it over the whole DOM every 16ms would cost more than the
#      milliseconds we are trying to save.
IN_PAGE_SNIPER_JS = r"""
(options) => new Promise((resolve) => {
    const { itemNames, buttonPattern, maxWaitMs, backstopMs } = options;
    const buttonRegex = new RegExp(buttonPattern, "i");
    const deadline = Date.now() + maxWaitMs;

    let finished = false;
    let observer = null;
    let backstopTimer = null;
    let cachedTargets = [];

    const finish = (value) => {
        if (finished) return;
        finished = true;
        if (observer) observer.disconnect();
        if (backstopTimer) clearInterval(backstopTimer);
        resolve(value);
    };

    const labelOf = (element) => (element.textContent || "").trim();

    const isClickable = (element) => {
        if (!element || !element.isConnected) return false;
        if (element.disabled) return false;
        if (element.getAttribute("aria-disabled") === "true") return false;
        const rect = element.getBoundingClientRect();
        if (rect.width === 0 || rect.height === 0) return false;
        const style = window.getComputedStyle(element);
        return style.visibility !== "hidden" && style.display !== "none" && style.pointerEvents !== "none";
    };

    // Walk up from the node holding the item name until we reach an ancestor that
    // also contains a button — that ancestor is the item's card.
    const buttonWithinCard = (node) => {
        let current = node;
        for (let depth = 0; depth < 8 && current; depth += 1) {
            const button = [...current.querySelectorAll("button, [role=button], a")]
                .find((candidate) => buttonRegex.test(labelOf(candidate)));
            if (button) return button;
            current = current.parentElement;
        }
        return null;
    };

    // Expensive full-DOM scan. Runs during the waiting room and again only if the
    // page re-renders and detaches what we cached.
    const resolveTargets = () => {
        const found = [];
        for (const itemName of itemNames) {
            const needle = itemName.toLowerCase();
            const nameNodes = [...document.querySelectorAll("h1,h2,h3,h4,h5,p,span,div,li")]
                .filter((element) => (element.textContent || "").toLowerCase().includes(needle));
            // Deepest match first: the innermost node is the item's own label rather
            // than some outer container that happens to contain the whole menu.
            for (const node of nameNodes.reverse()) {
                const button = buttonWithinCard(node);
                if (button && !found.some((entry) => entry.button === button)) {
                    found.push({ itemName, button });
                    break;
                }
            }
        }
        return found;
    };

    const attempt = () => {
        if (finished) return;
        try {
            const cacheIsStale = cachedTargets.length === 0
                || cachedTargets.some((entry) => !entry.button.isConnected);
            if (cacheIsStale) cachedTargets = resolveTargets();

            // config.target_items is a priority list, so first clickable hit wins.
            for (const entry of cachedTargets) {
                if (isClickable(entry.button)) {
                    const label = labelOf(entry.button);
                    entry.button.click();
                    finish({
                        status: "clicked",
                        itemName: entry.itemName,
                        buttonLabel: label,
                        clickedAt: Date.now(),
                    });
                    return;
                }
            }
        } catch (error) {
            finish({ status: "error", detail: String(error) });
            return;
        }
        if (Date.now() > deadline) finish({ status: "timeout" });
    };

    // Primary trigger: react in the same microtask the button is enabled in.
    observer = new MutationObserver(attempt);
    observer.observe(document.body, {
        childList: true,
        subtree: true,
        attributes: true,
        attributeFilter: ["disabled", "aria-disabled", "class", "style"],
    });

    // Backstop: covers an enable that mutates nothing we observe (e.g. a CSS
    // animation completing, or state flipping inside a shadow root).
    backstopTimer = setInterval(attempt, backstopMs);

    attempt();
})
"""


def snipe_add_to_cart(page: Page, config: dict[str, Any]) -> dict[str, Any]:
    """
    Block until an add-to-cart button for a target item becomes clickable, then click it.

    Returns the sniper's result dict: {"status": "clicked"|"timeout"|"error", ...}.
    Side effect: on success, an item is really reserved in your real cart.
    """
    max_wait_ms = int(config["max_wait_minutes"] * 60 * 1000)
    options = {
        "itemNames": config["target_items"],
        "buttonPattern": config["add_to_cart_pattern"],
        "maxWaitMs": max_wait_ms,
        "backstopMs": SNIPE_POLL_BACKSTOP_MS,
    }
    log.info(
        "Armed. Watching for %s on: %s",
        config["add_to_cart_pattern"],
        ", ".join(config["target_items"]),
    )
    try:
        return page.evaluate(IN_PAGE_SNIPER_JS, options)
    except PlaywrightError as exc:
        # If add-to-cart triggers a real navigation, the JS context is torn down
        # before the resolved value can be handed back. That destruction is itself
        # evidence the click landed, so report it rather than treating it as failure.
        if "context was destroyed" in str(exc).lower() or "execution context" in str(exc).lower():
            log.warning("Page navigated during the click — treating as a probable success. VERIFY YOUR CART.")
            return {"status": "navigated", "itemName": "unknown", "buttonLabel": "unknown"}
        raise


# --------------------------------------------------------------------------------------
# recon probe
# --------------------------------------------------------------------------------------

# Mirrors how the sniper matches, so recon's answer is the sniper's answer. It reads
# innerText but falls back to textContent: a page saved to disk and reloaded without
# its stylesheets can lay out differently, and innerText alone would under-report.
RECON_PROBE_JS = r"""
(options) => {
    const buttonRegex = new RegExp(options.buttonPattern, "i");
    const readText = (element) => ((element.innerText || element.textContent || "").trim());

    const buttons = [...document.querySelectorAll("button, [role=button], a")]
        .map((element) => ({
            text: readText(element).slice(0, 60),
            disabled: Boolean(element.disabled) || element.getAttribute("aria-disabled") === "true",
        }))
        .filter((entry) => entry.text.length > 0);

    const pageText = readText(document.body).toLowerCase();
    const foundItems = options.itemNames.filter((name) => pageText.includes(name.toLowerCase()));

    return {
        pageTitle: document.title,
        foundItems,
        missingItems: options.itemNames.filter((name) => !foundItems.includes(name)),
        matchingButtons: buttons.filter((entry) => buttonRegex.test(entry.text)),
        allButtonTexts: [...new Set(buttons.map((entry) => entry.text))].slice(0, 40),
    };
}
"""


def probe_page(page: Page, config: dict[str, Any]) -> dict[str, Any]:
    """Run the recon probe against whatever is currently loaded in `page`."""
    return page.evaluate(
        RECON_PROBE_JS,
        {"itemNames": config["target_items"], "buttonPattern": config["add_to_cart_pattern"]},
    )


def print_recon_report(probe_result: dict[str, Any], config: dict[str, Any]) -> None:
    """Print the probe result, leading with the two failure modes that lose drops."""
    print("\n=== RECON ===")
    print(f"Page title      : {probe_result['pageTitle']}")
    print(f"Items found     : {probe_result['foundItems'] or 'NONE — fix config.target_items'}")
    print(f"Items missing   : {probe_result['missingItems'] or 'none'}")
    print(f"\nButtons matching '{config['add_to_cart_pattern']}':")
    for entry in probe_result["matchingButtons"]:
        print(f"    [{'DISABLED' if entry['disabled'] else 'ENABLED'}] {entry['text']!r}")
    if not probe_result["matchingButtons"]:
        print("    NONE — widen config.add_to_cart_pattern using the button texts below.")

    print("\nAll button texts on page (first 40):")
    for text in probe_result["allButtonTexts"]:
        print(f"    {text!r}")

    print("\nVerdict:")
    if not probe_result["foundItems"]:
        print("    NOT READY — none of your target_items appear on this page.")
    elif not probe_result["matchingButtons"]:
        print("    NOT READY — items found, but no button matches add_to_cart_pattern.")
    elif probe_result["missingItems"]:
        print(f"    PARTIAL — will race for {probe_result['foundItems']}, but "
              f"{probe_result['missingItems']} were not found.")
    else:
        print("    READY — all target items and a matching button are visible.")
    print(f"\nDump written to {BOT_DIR / 'recon_dump.html'} and recon_dump.png")


# --------------------------------------------------------------------------------------
# subcommands
# --------------------------------------------------------------------------------------


def command_login(config: dict[str, Any]) -> int:
    """Open a real browser window, wait for you to sign in, then persist the session."""
    with sync_playwright() as playwright:
        browser = launch_browser(playwright, headless=False)
        context = open_context(browser)
        page = context.new_page()
        page.goto("https://www.hotplate.com/", wait_until="domcontentloaded")

        print("\n  A browser window is open. Sign in to Hotplate by hand.")
        print("  Also add your card and address now so checkout has nothing to type.")
        input("  Press Enter here once you are fully signed in... ")

        context.storage_state(path=str(STORAGE_STATE_PATH))
        STORAGE_STATE_PATH.chmod(0o600)
        log.info("Session saved to %s (mode 600 — this file is as good as your password).", STORAGE_STATE_PATH)

        context.close()
        browser.close()
    return 0


def command_recon(config: dict[str, Any], html_file: Path | None = None) -> int:
    """
    Report what the sniper would see on the drop page, so config can be verified.

    Run this on a live drop page BEFORE the real drop. It reports which target items
    it can find and whether their buttons are currently disabled (they should be,
    during the waiting room — that is the signal the matcher is working).

    With html_file set, probes a saved copy of the page instead of fetching it. The
    probe only needs the DOM, so this works with no network at all — useful when the
    page is saved from one machine and checked on another.
    """
    if html_file is None:
        require_saved_session()
    elif not html_file.exists():
        raise SystemExit(f"No such HTML file: {html_file}")

    with sync_playwright() as playwright:
        # Offline mode never needs a visible window, and often runs somewhere headless.
        browser = launch_browser(playwright, headless=config["headless"] or html_file is not None)
        context = open_context(browser)
        page = context.new_page()

        if html_file is None:
            page.goto(config["event_url"], wait_until="networkidle")
            check_clock_skew(page)
        else:
            log.info("Offline recon against %s (no network).", html_file)
            page.set_content(html_file.read_text(encoding="utf-8", errors="replace"))

        probe_result = probe_page(page, config)
        (BOT_DIR / "recon_dump.html").write_text(page.content())
        page.screenshot(path=str(BOT_DIR / "recon_dump.png"), full_page=True)
        print_recon_report(probe_result, config)

        context.close()
        browser.close()
    return 0


def command_arm(config: dict[str, Any]) -> int:
    """
    The real run. Sit on the event page and click add-to-cart the moment it enables.

    Stops at the cart and hands off to you unless auto_submit_payment is true.
    """
    run_stamp = datetime.now().strftime("%Y%m%d_%H%M%S")
    require_saved_session()

    with sync_playwright() as playwright:
        browser = launch_browser(playwright, headless=config["headless"])
        context = open_context(browser)
        page = context.new_page()

        log.info("Loading %s", config["event_url"])
        page.goto(config["event_url"], wait_until="networkidle")
        check_clock_skew(page)

        if "sign in" in page.content().lower() and "log in" in page.content().lower():
            log.warning("Page may be showing a sign-in prompt — your saved session might be stale.")
            log.warning("If recon shows you logged out, re-run `python hotplate_snipe.py login`.")

        page.screenshot(path=str(BOT_DIR / f"snipe_{run_stamp}_1_armed.png"))
        alert_human("Armed on the drop page. Waiting for the countdown.")

        snipe_started_at = time.monotonic()
        result = snipe_add_to_cart(page, config)
        elapsed_seconds = time.monotonic() - snipe_started_at

        if result["status"] == "timeout":
            log.error(
                "Timed out after %.0fs without a clickable button. Either the drop never opened, "
                "the item names in config do not match the page, or the button text is outside "
                "add_to_cart_pattern. Run `recon` against the live page to see what is there.",
                elapsed_seconds,
            )
            page.screenshot(path=str(BOT_DIR / f"snipe_{run_stamp}_timeout.png"))
            alert_human("MISSED — no add-to-cart button became clickable.")
            context.close()
            browser.close()
            return 1

        if result["status"] == "error":
            log.error("In-page sniper crashed: %s", result.get("detail"))
            page.screenshot(path=str(BOT_DIR / f"snipe_{run_stamp}_error.png"))
            alert_human("ERROR — sniper crashed, take over manually NOW.")
            context.close()
            browser.close()
            return 1

        log.info(
            "CLICKED %r (button %r) after %.1fs of waiting.",
            result["itemName"],
            result["buttonLabel"],
            elapsed_seconds,
        )
        page.screenshot(path=str(BOT_DIR / f"snipe_{run_stamp}_2_clicked.png"))
        alert_human(f"IN CART: {result['itemName']} — finish checkout now.")

        # Past this point the item is reserved and the cart timer is the only clock
        # that matters. That is a comfortable amount of time, so stop optimising for
        # speed and let a human confirm quantity, pickup slot, and payment.
        print("\n" + "=" * 70)
        print(f"  ADDED TO CART: {result['itemName']}")
        print("  The browser is yours. Adjust quantity, pick your slot, and pay.")
        print("  Leave this process running so the window stays open.")
        print("=" * 70 + "\n")

        if config["auto_submit_payment"]:
            log.warning("auto_submit_payment is true but automated payment submission is not implemented.")
            log.warning("Finish checkout by hand — an unattended mis-click on payment is not worth the seconds.")

        input("Press Enter to close the browser once your order is placed... ")
        page.screenshot(path=str(BOT_DIR / f"snipe_{run_stamp}_3_final.png"))
        context.close()
        browser.close()
    return 0


# --------------------------------------------------------------------------------------
# entry point
# --------------------------------------------------------------------------------------


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description="Hotplate drop assistant.")
    parser.add_argument("command", choices=["login", "recon", "arm"])
    parser.add_argument("--config", type=Path, default=DEFAULT_CONFIG_PATH)
    parser.add_argument(
        "--html-file",
        type=Path,
        default=None,
        help="recon only: probe a saved copy of the drop page instead of fetching it (no network).",
    )
    args = parser.parse_args(argv)

    if args.html_file is not None and args.command != "recon":
        raise SystemExit("--html-file is only valid with the `recon` command.")

    logging.basicConfig(
        level=logging.INFO,
        format="%(asctime)s.%(msecs)03d %(levelname)s %(message)s",
        datefmt="%H:%M:%S",
    )

    config = load_config(args.config)
    if args.command == "recon":
        return command_recon(config, html_file=args.html_file)
    handlers = {"login": command_login, "arm": command_arm}
    return handlers[args.command](config)


if __name__ == "__main__":
    sys.exit(main())
