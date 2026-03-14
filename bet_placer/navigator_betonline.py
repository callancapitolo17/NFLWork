#!/usr/bin/env python3
"""
BetOnline bet navigator.
Opens real Chrome (bypasses Cloudflare), navigates to CBB schedule,
clicks odds buttons to build betslip, fills amounts.

BetOnline DOM patterns (from sportsbook recon):
  - CBB URL: /sportsbook/basketball/ncaa
  - Game rows: a.event-row-container (href=/sportsbook/basketball/ncaa/game/{id})
  - Team names: p[data-testid="participant-name"] (short names: "Wake Forest", "Clemson")
  - Odds buttons: div.odd.medium (data-testid="odd-component-{odds}")
    - Spread: two <p> tags — line ("+1.5") and odds ("-120")
    - ML: single <p> — just odds ("+100")
    - Total: "O 237.5" / "U 237.5" prefix
  - Betslip: div.betslip-container (right sidebar)
  - Amount input: Risk input in betslip after selection is added

Auth: Keycloak OAuth2 via persistent Chrome profile. Cloudflare requires
real Chrome channel (not Playwright Chromium).
"""

import os
import time

from playwright.sync_api import sync_playwright

from base_navigator import (BaseNavigator, parse_market, update_bet_status,
                           cleanup_singleton_lock, wait_for_confirmation, _REPO_ROOT)

BETONLINE_URL = "https://www.betonline.ag"
CBB_URL = f"{BETONLINE_URL}/sportsbook/basketball/ncaa"

BETONLINE_USERNAME = os.getenv("BETONLINE_USERNAME")
BETONLINE_PASSWORD = os.getenv("BETONLINE_PASSWORD")

# Persistent browser profile — reuse bet_logger's profile if available
# (already has Cloudflare cookies + Keycloak session)
_BET_LOGGER_PROFILE = str(_REPO_ROOT / "bet_logger" / ".betonline_profile")
_NAV_PROFILE = str(_REPO_ROOT / "bet_placer" / ".betonline_nav_profile")
BETONLINE_PROFILE_DIR = _BET_LOGGER_PROFILE if os.path.isdir(_BET_LOGGER_PROFILE) else _NAV_PROFILE


class BetOnlineNavigator(BaseNavigator):

    BOOK_NAME = "betonline"

    def place_bet(self, bet_data: dict):
        """Place a single bet. Delegates to place_bets()."""
        self.place_bets([bet_data])

    def place_bets(self, bets: list[dict], timeout: int = 300):
        """Place multiple bets in a single BetOnline browser session.

        Strategy:
        1. Open Chrome with persistent profile (Cloudflare + Keycloak session)
        2. Navigate to CBB schedule page
        3. For each bet: find game row, click correct odds button
        4. Fill amounts in betslip
        5. Wait for user confirmation via dashboard (or timeout)
        """
        print(f"\n  BetOnline Batch: {len(bets)} bet(s)", flush=True)
        for i, bet in enumerate(bets):
            print(f"  [{i+1}] {self._format_bet_summary(bet)}", flush=True)

        # Clean up stale browser lock before launching
        cleanup_singleton_lock(BETONLINE_PROFILE_DIR)

        with sync_playwright() as p:
            # Must use real Chrome to bypass Cloudflare
            context = p.chromium.launch_persistent_context(
                user_data_dir=BETONLINE_PROFILE_DIR,
                channel="chrome",
                headless=False,
                no_viewport=True,
                args=[
                    "--disable-blink-features=AutomationControlled",
                    "--window-size=1400,900",
                    "--window-position=50,50",
                ],
            )
            page = context.pages[0] if context.pages else context.new_page()

            # Navigate to CBB
            print(f"  Navigating to {CBB_URL}...", flush=True)
            page.goto(CBB_URL, wait_until="domcontentloaded", timeout=60000)

            # Wait for Cloudflare challenge to pass
            self._wait_for_cloudflare(page)

            # Check auth
            if not self._is_logged_in(page):
                self._login(page)
                page.goto(CBB_URL, wait_until="domcontentloaded", timeout=30000)

            # Wait for game rows to load
            self._wait_for_games(page)

            # Clear any leftover bets from previous sessions
            self._clear_betslip(page)

            # Click odds for each bet
            clicked_bets = []
            for i, bet in enumerate(bets):
                parsed = parse_market(bet["market"])
                print(f"\n  --- Bet {i+1}/{len(bets)} ---", flush=True)
                print(f"  {self._format_bet_summary(bet)}", flush=True)

                try:
                    clicked = self._click_odds_for_bet(page, bet, parsed)
                    if clicked:
                        clicked_bets.append(bet)
                        time.sleep(2)  # Wait for betslip to update
                    else:
                        self._print_manual_instructions(bet, parsed)
                        if bet.get("bet_hash"):
                            update_bet_status(bet["bet_hash"], "nav_error", only_if="navigating")
                except Exception as e:
                    print(f"  Failed to add bet: {e}", flush=True)
                    if bet.get("bet_hash"):
                        update_bet_status(bet["bet_hash"], "nav_error", only_if="navigating")

            # Fill amounts after all bets are in betslip
            if clicked_bets:
                print(f"\n  Filling amounts for {len(clicked_bets)} bets...", flush=True)
                time.sleep(2)

                # BetOnline auto-switches to Parlay tab with 2+ bets.
                # Click "Straight" tab to fill individual amounts.
                if len(clicked_bets) > 1:
                    self._switch_to_straight_tab(page)

                self._fill_all_amounts(page, clicked_bets)

            print(f"\n  {len(clicked_bets)}/{len(bets)} bets added to betslip.", flush=True)
            print(f"  Confirm via dashboard. Browser open for {timeout}s.", flush=True)

            # Wait for user to confirm bets via dashboard, or timeout
            hashes = [b["bet_hash"] for b in clicked_bets if b.get("bet_hash")]
            if hashes:
                wait_for_confirmation(hashes, timeout=timeout)

            context.close()

    # =========================================================================
    # AUTH
    # =========================================================================

    def _wait_for_cloudflare(self, page):
        """Wait for Cloudflare challenge to pass (usually automatic with persistent profile)."""
        for i in range(30):
            url = page.url
            title = page.title()
            # Cloudflare challenge pages have specific titles
            if "just a moment" in title.lower() or "challenge" in title.lower():
                if i == 0:
                    print("  Waiting for Cloudflare...", flush=True)
                time.sleep(2)
                continue
            # Check if sportsbook content loaded
            if "sportsbook" in url.lower() or "betonline" in title.lower():
                if i > 0:
                    print(f"  Cloudflare passed ({i*2}s)", flush=True)
                return
            time.sleep(1)
        print("  Cloudflare may not have passed — continuing anyway", flush=True)

    def _is_logged_in(self, page) -> bool:
        """Check if user is authenticated."""
        try:
            # Check for balance display (logged in indicator)
            balance = page.locator('[data-testid="balance"], .balance, [class*="balance"]')
            if balance.count() > 0 and balance.first.is_visible(timeout=5000):
                print("  Already authenticated (session active)", flush=True)
                return True
        except Exception:
            pass

        try:
            # Check for login/signup buttons (not logged in)
            login_btn = page.locator('button:has-text("Log In"), a:has-text("Log In"), [data-testid="login"]')
            if login_btn.count() > 0 and login_btn.first.is_visible(timeout=3000):
                print("  Not logged in (login button visible)", flush=True)
                return False
        except Exception:
            pass

        # Ambiguous — check for deposit button (only visible when logged in)
        try:
            deposit = page.locator('button:has-text("Deposit"), a:has-text("Deposit")')
            if deposit.count() > 0 and deposit.first.is_visible(timeout=3000):
                print("  Already authenticated (deposit button visible)", flush=True)
                return True
        except Exception:
            pass

        # Default to logged in (persistent profile usually works)
        print("  Auth status unclear — assuming logged in", flush=True)
        return True

    def _login(self, page):
        """Log in to BetOnline via browser (Keycloak)."""
        print("  Logging in to BetOnline...", flush=True)

        try:
            # Click login button
            login_btn = page.locator('button:has-text("Log In"), a:has-text("Log In"), [data-testid="login"]')
            if login_btn.count() > 0:
                login_btn.first.click()
                time.sleep(3)

            # Fill Keycloak login form
            url = page.url
            if "auth" in url.lower() or "login" in url.lower():
                # Keycloak form
                user_field = page.locator('input[name="username"], input[type="text"]').first
                pass_field = page.locator('input[name="password"], input[type="password"]').first

                if BETONLINE_USERNAME and BETONLINE_PASSWORD:
                    user_field.click()
                    page.keyboard.type(BETONLINE_USERNAME, delay=50)
                    pass_field.click()
                    page.keyboard.type(BETONLINE_PASSWORD, delay=50)

                    # Submit
                    submit = page.locator('input[type="submit"], button[type="submit"], #kc-login')
                    if submit.count() > 0:
                        submit.first.click()
                    else:
                        page.keyboard.press("Enter")

                    time.sleep(5)
                    print(f"  Login submitted. URL: {page.url}", flush=True)
                    return

            print("  Auto-login failed. Please log in manually (60 seconds).", flush=True)
            time.sleep(60)

        except Exception as e:
            print(f"  Login error: {e}", flush=True)
            print("  Please log in manually (60 seconds).", flush=True)
            time.sleep(60)

    # =========================================================================
    # BETSLIP MANAGEMENT
    # =========================================================================

    def _switch_to_straight_tab(self, page):
        """Switch betslip from Parlay to Straight tab.

        Uses Playwright locator (not evaluate) — works inside simplebar wrapper.
        """
        try:
            straight = page.locator('text="Straight"').first
            if straight.is_visible(timeout=2000):
                straight.click()
                time.sleep(1)
                print("  Switched to Straight tab", flush=True)
                return
        except Exception:
            pass

    def _clear_betslip(self, page):
        """Remove all existing bets from the betslip.

        Uses Playwright locators (not evaluate) — works inside simplebar wrapper.
        1. Click the trash icon (bet-summary__icon) to remove all bets
        2. Handle confirmation dialog
        3. Clear bulk "All Selections" inputs
        """
        try:
            trash = page.locator('[class*="bet-summary__icon"]').first
            if trash.is_visible(timeout=2000):
                trash.click()
                time.sleep(1)
                # Handle confirmation dialog if one appears
                for label in ("Remove", "Yes", "OK", "Clear", "Delete"):
                    try:
                        btn = page.locator(f'button:has-text("{label}")')
                        if btn.count() > 0 and btn.first.is_visible(timeout=1000):
                            btn.first.click()
                            time.sleep(1)
                            break
                    except Exception:
                        continue
                print("  Cleared leftover bets from betslip", flush=True)
        except Exception:
            pass

        self._clear_bulk_inputs(page)

    def _clear_bulk_inputs(self, page):
        """Clear the Risk All / Win All Selections bulk inputs.

        BetOnline betslip structure (Svelte):
        - Bulk area has <span> "Risk All Selections" / "Win All Selections"
          followed by div.risk__amount > input.risk__input
        - Per-bet area has <span> "Risk" / "Win"
          followed by div.risk__amount > input.risk__input
        """
        # Bulk inputs are the first two risk__input elements (Risk All, Win All)
        # Per-bet inputs follow after. We locate bulk by their parent container
        # having the "All Selections" text.
        try:
            # Find containers with "All Selections" text, then find inputs inside
            bulk_containers = page.locator(':text("All Selections")').locator('xpath=..')
            bc_count = bulk_containers.count()
            for i in range(bc_count):
                container = bulk_containers.nth(i)
                inp = container.locator('input.risk__input')
                if inp.count() > 0 and inp.first.is_visible(timeout=500):
                    val = inp.first.input_value(timeout=500)
                    if val and val != "0.00" and val != "" and val != "0":
                        inp.first.click()
                        time.sleep(0.2)
                        page.keyboard.press("Meta+a")
                        page.keyboard.press("Backspace")
                        page.keyboard.press("Tab")
                        time.sleep(0.3)
                        print(f"  Cleared bulk input (was ${val})", flush=True)
        except Exception:
            pass

    def _count_betslip_inputs(self, page) -> int:
        """Count per-bet Risk inputs in the betslip.

        BetOnline Svelte betslip:
        - All inputs have class risk__input
        - 4 visible = 2 bulk (Risk All, Win All) + 2 per-bet (Risk, Win) for 1 bet
        - 6 visible = 2 bulk + 4 per-bet (Risk, Win) × 2 for 2 bets
        - Per-bet count = (total visible risk__input - 2 bulk) / 2
        """
        try:
            total = page.locator('input.risk__input').count()
            visible = 0
            for i in range(total):
                if page.locator('input.risk__input').nth(i).is_visible(timeout=300):
                    visible += 1
            # 2 bulk inputs always present, remaining are pairs (Risk + Win) per bet
            per_bet = max(0, (visible - 2) // 2)
            return per_bet
        except Exception:
            return 0

    # =========================================================================
    # NAVIGATION
    # =========================================================================

    def _wait_for_games(self, page):
        """Wait for game rows to render on the schedule page."""
        for i in range(20):
            time.sleep(1)
            count = page.locator("a.event-row-container").count()
            if count > 0:
                print(f"  Schedule loaded ({count} games, {i+1}s)", flush=True)
                return
        print("  Warning: no game rows found — page may not have loaded", flush=True)

    # =========================================================================
    # ODDS CLICKING
    # =========================================================================

    def _click_odds_for_bet(self, page, bet_data: dict, parsed) -> bool:
        """Find and click the correct odds button for a bet.

        Strategy:
        1. Find all game rows (a.event-row-container)
        2. Match by team names in participant-name elements
        3. Within matched game row, find the correct odds button by market type + value
        4. Click with React-compatible event dispatch
        """
        away = bet_data["away_team"]
        home = bet_data["home_team"]
        bet_on = bet_data["bet_on"]
        line = bet_data.get("line")
        odds = bet_data.get("odds")

        # Find matching game row
        game_row = self._find_game_row(page, away, home)
        if not game_row:
            return False

        # Find and click the correct odds button within this game row
        return self._click_odds_in_row(page, game_row, bet_data, parsed)

    def _find_game_row(self, page, away: str, home: str):
        """Find the game row matching the given teams."""
        rows = page.locator("a.event-row-container")
        count = rows.count()

        away_terms = self._search_terms(away)
        home_terms = self._search_terms(home)

        for i in range(count):
            row = rows.nth(i)
            try:
                # Get all participant names in this row
                names = row.locator('[data-testid="participant-name"]')
                name_count = names.count()
                if name_count < 2:
                    continue

                name_texts = []
                for j in range(name_count):
                    name_texts.append(names.nth(j).inner_text(timeout=1000).strip())

                # Match: first participant = away, second = home
                row_text = " ".join(name_texts).upper()
                has_away = any(t.upper() in row_text for t in away_terms)
                has_home = any(t.upper() in row_text for t in home_terms)

                if has_away and has_home:
                    row.scroll_into_view_if_needed()
                    print(f"  Found game: {' @ '.join(name_texts)}", flush=True)
                    return row

            except Exception:
                continue

        print(f"  Could not find game: {away} @ {home}", flush=True)
        # Diagnostic
        try:
            for i in range(min(count, 10)):
                row = rows.nth(i)
                names = row.locator('[data-testid="participant-name"]')
                texts = [names.nth(j).inner_text(timeout=500).strip()
                         for j in range(names.count())]
                print(f"    [{i}] {' @ '.join(texts)}", flush=True)
        except Exception:
            pass
        return None

    def _click_odds_in_row(self, page, game_row, bet_data: dict, parsed) -> bool:
        """Click the correct odds button within a game row.

        BetOnline layout:
        - Each game row has 12 odds buttons (FG + 1H):
          [away_spread, home_spread, away_ml, home_ml, over, under] × 2 (FG then 1H)
        - Spread: div.odd with two <p> tags (line + odds)
        - ML: div.odd with one <p> tag (just odds)
        - Total: div.odd with "O"/"U" prefix on line
        """
        bet_on = bet_data["bet_on"]
        line = bet_data.get("line")
        odds = bet_data.get("odds")

        # Get all odds buttons in this row
        odds_btns = game_row.locator("div.odd")
        btn_count = odds_btns.count()

        if btn_count == 0:
            print(f"  No odds buttons in game row", flush=True)
            return False

        # Determine which half of buttons to search based on period
        # First 6 = FG, next 6 = 1H
        is_first_half = parsed.period == "1h"
        if btn_count == 12:
            search_start = 6 if is_first_half else 0
            search_end = 12 if is_first_half else 6
        else:
            # Fallback: search all buttons
            search_start = 0
            search_end = btn_count

        # Score each button in the relevant half
        best_btn = None
        best_score = 0

        is_away = (bet_on.lower() in bet_data["away_team"].lower()
                   or bet_data["away_team"].lower() in bet_on.lower())

        for j in range(search_start, search_end):
            btn = odds_btns.nth(j)
            try:
                text = btn.evaluate("el => el.textContent").strip().replace("\n", "")
            except Exception:
                continue

            if not text:
                continue

            score = self._score_odds_match(text, parsed, bet_on, line, odds,
                                           bet_data["away_team"], bet_data["home_team"])

            # Position bonus: within each 6-button group
            # Order: [away_spread, home_spread, away_ml, home_ml, over, under]
            local_idx = j - search_start
            expected_idx = self._expected_button_index(parsed, bet_on, is_away)
            if expected_idx is not None and local_idx == expected_idx:
                score += 2

            if score > best_score:
                best_score = score
                best_btn = btn

        if best_btn and best_score >= 3:
            try:
                text = best_btn.evaluate("el => el.textContent").strip().replace("\n", "")
                # Count betslip entries before click
                bets_before = self._count_betslip_entries(page)

                best_btn.scroll_into_view_if_needed()
                time.sleep(0.3)
                best_btn.click()
                time.sleep(1.5)

                # Check if a new bet appeared in the betslip
                bets_after = self._count_betslip_entries(page)
                if bets_after > bets_before:
                    print(f"  Clicked odds: '{text}' (score={best_score}, betslip: {bets_before}→{bets_after})", flush=True)
                    return True

                # Click registered but betslip count didn't change — may still work
                print(f"  Clicked odds: '{text}' (score={best_score}, betslip unchanged at {bets_after})", flush=True)
                return True

            except Exception as e:
                print(f"  Click failed: {e}", flush=True)

        print(f"  No matching odds button (best_score={best_score})", flush=True)
        # Debug: show buttons in the searched range
        for j in range(search_start, min(search_end, search_start + 10)):
            try:
                text = odds_btns.nth(j).evaluate("el => el.textContent").strip().replace("\n", "")
                if text:
                    print(f"    [{j}] '{text}'", flush=True)
            except Exception:
                pass
        return False

    def _count_betslip_entries(self, page) -> int:
        """Count how many bet entries are in the betslip."""
        return self._count_betslip_inputs(page)

    def _expected_button_index(self, parsed, bet_on: str, is_away: bool) -> int | None:
        """Expected button position in BetOnline's 6-button layout.

        Order: [away_spread, home_spread, away_ml, home_ml, over, under]
        """
        if parsed.market_type in ("spreads", "alternate_spreads"):
            return 0 if is_away else 1
        elif parsed.market_type == "h2h":
            return 2 if is_away else 3
        elif parsed.market_type in ("totals", "alternate_totals"):
            return 4 if bet_on == "Over" else 5
        return None

    def _score_odds_match(self, text: str, parsed, bet_on, line, odds, away, home) -> int:
        """Score how well button text matches target bet."""
        score = 0
        text = text.replace("\n", "")
        text_lower = text.lower()

        odds_str = str(odds) if odds else ""
        try:
            if odds and int(odds) > 0:
                odds_str = f"+{odds}"
        except (ValueError, TypeError):
            pass

        if line is not None:
            abs_line = abs(line)
            line_str = str(abs_line) if abs_line == int(abs_line) else str(abs_line)

            if parsed.market_type in ("spreads", "alternate_spreads"):
                sign = "+" if line >= 0 else "-"
                target = f"{sign}{line_str}"
                if target in text:
                    score += 3

            elif parsed.market_type in ("totals", "alternate_totals"):
                prefix = "O" if bet_on == "Over" else "U"
                if prefix in text and line_str in text:
                    score += 3

            if line_str in text:
                score += 1

        if odds_str and odds_str in text:
            score += 2

        # Total direction
        if parsed.market_type in ("totals", "alternate_totals"):
            if bet_on == "Over" and ("O " in text or text.startswith("O")):
                score += 1
            elif bet_on == "Under" and ("U " in text or text.startswith("U")):
                score += 1

        return score

    @staticmethod
    def _search_terms(team_name: str) -> list[str]:
        """Generate search terms from canonical team name.

        BetOnline uses short names without mascots (e.g., "Wake Forest", "Oklahoma State").
        """
        words = team_name.split()
        mascots = {"wildcats", "cougars", "hokies", "bulldogs", "eagles", "bears",
                   "tigers", "hawks", "huskies", "cardinals", "knights", "rams",
                   "wolves", "panthers", "broncos", "bobcats", "terriers", "gaels",
                   "warriors", "demon", "deacons", "blue", "devils", "wave",
                   "green", "orange", "red", "golden", "fighting", "irish"}

        while words and words[-1].lower() in mascots:
            words = words[:-1]

        terms = []
        for n in range(len(words), 0, -1):
            prefix = " ".join(words[:n])
            if len(prefix) >= 3:
                terms.append(prefix)

        return terms if terms else [team_name]

    # =========================================================================
    # BETSLIP + AMOUNT FILLING
    # =========================================================================

    def _fill_all_amounts(self, page, bets: list[dict]):
        """Fill risk amounts in the betslip for all added bets.

        BetOnline betslip DOM (Svelte, NOT MUI):
        - All amount inputs share class: input.risk__input
        - Layout order (top to bottom):
            [0] Risk All Selections  (bulk)
            [1] Win All Selections   (bulk)
            [2] Risk  (bet 1)
            [3] Win   (bet 1)
            [4] Risk  (bet 2)
            [5] Win   (bet 2)
            ...
        - Per-bet Risk inputs are at indices 2, 4, 6, ... (every other, starting at 2)
        - Parent structure: <span>label</span> + <div class="risk__amount"> <input>
        """
        sizes = [int(b.get("recommended_size", 0)) for b in bets]

        # Take screenshot before filling
        try:
            page.screenshot(path="/tmp/betonline_before_fill.png")
            print("  Screenshot: /tmp/betonline_before_fill.png", flush=True)
        except Exception:
            pass

        # Clear bulk inputs first
        self._clear_bulk_inputs(page)

        # Get all visible risk__input elements
        all_risk = page.locator('input.risk__input')
        total = all_risk.count()

        # Filter to visible only
        visible_indices = []
        for i in range(total):
            try:
                if all_risk.nth(i).is_visible(timeout=300):
                    visible_indices.append(i)
            except Exception:
                pass

        print(f"  {len(visible_indices)} visible risk__input elements (2 bulk + {len(visible_indices)-2} per-bet)", flush=True)

        # Per-bet Risk inputs start at index 2 (skip bulk), every other (Risk, Win, Risk, Win...)
        per_bet_risk_indices = visible_indices[2::2]  # indices 2, 4, 6, ...

        print(f"  {len(per_bet_risk_indices)} per-bet Risk inputs for {len(bets)} bets", flush=True)

        if not per_bet_risk_indices:
            print("  No per-bet Risk inputs found!", flush=True)
            # Debug
            for i in visible_indices:
                try:
                    parent_text = all_risk.nth(i).evaluate(
                        "el => el.parentElement.parentElement.querySelector('span')?.textContent || 'no-span'"
                    )
                    print(f"    input[{i}] label='{parent_text}'", flush=True)
                except Exception:
                    pass
            return

        # Fill each per-bet Risk input
        for bet_num, inp_idx in enumerate(per_bet_risk_indices):
            if bet_num >= len(sizes) or not sizes[bet_num]:
                continue
            val_str = str(sizes[bet_num])
            field = all_risk.nth(inp_idx)
            try:
                field.click()
                time.sleep(0.3)
                page.keyboard.press("Meta+a")
                page.keyboard.press("Backspace")
                time.sleep(0.1)
                page.keyboard.type(val_str, delay=80)
                page.keyboard.press("Tab")
                time.sleep(0.5)
                print(f"  Bet {bet_num+1}: ${sizes[bet_num]} filled", flush=True)
                if bet_num < len(bets) and bets[bet_num].get("bet_hash"):
                    update_bet_status(bets[bet_num]["bet_hash"], "ready_to_confirm")
                # Clear bulk inputs after each fill to prevent auto-calc override
                self._clear_bulk_inputs(page)
            except Exception as e:
                print(f"  Bet {bet_num+1}: failed to fill ${sizes[bet_num]}: {e}", flush=True)

        # Final clear of bulk inputs
        self._clear_bulk_inputs(page)

        # Verify: re-read per-bet Risk values
        for bet_num, inp_idx in enumerate(per_bet_risk_indices):
            if bet_num >= len(sizes) or not sizes[bet_num]:
                continue
            try:
                actual = all_risk.nth(inp_idx).input_value(timeout=1000)
                expected = f"{sizes[bet_num]}.00"
                if actual != expected and actual != str(sizes[bet_num]):
                    print(f"  WARNING: Bet {bet_num+1} Risk=${actual}, expected ${sizes[bet_num]} — refilling", flush=True)
                    field = all_risk.nth(inp_idx)
                    field.click()
                    time.sleep(0.2)
                    page.keyboard.press("Meta+a")
                    page.keyboard.press("Backspace")
                    page.keyboard.type(str(sizes[bet_num]), delay=80)
                    page.keyboard.press("Tab")
                    time.sleep(0.3)
                    self._clear_bulk_inputs(page)
                else:
                    print(f"  Bet {bet_num+1}: verified ${actual}", flush=True)
            except Exception:
                pass

        # Take screenshot after filling
        try:
            page.screenshot(path="/tmp/betonline_after_fill.png")
            print("  Screenshot: /tmp/betonline_after_fill.png", flush=True)
        except Exception:
            pass

    def _print_manual_instructions(self, bet_data: dict, parsed):
        """Print instructions for manual placement."""
        bet_on = bet_data["bet_on"]
        line = bet_data.get("line")
        odds = bet_data.get("odds")
        size = bet_data.get("recommended_size", 0)

        print(f"\n  --- MANUAL PLACEMENT ---")
        print(f"  Game:   {bet_data['away_team']} @ {bet_data['home_team']}")
        print(f"  Market: {parsed.market_type.upper()} | Period: {parsed.period.upper()}")
        print(f"  Pick:   {bet_on}" + (f" {line}" if line else ""))
        print(f"  Odds:   {'+' if odds and odds > 0 else ''}{odds}")
        print(f"  Amount: ${size:.0f}")
        print(f"  ---")
