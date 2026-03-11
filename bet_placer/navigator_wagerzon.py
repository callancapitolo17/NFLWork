#!/usr/bin/env python3
"""
Wagerzon bet navigator.
Uses the NewScheduleHelper API to find internal game/line IDs,
then constructs a direct CreateWager URL to pre-fill the bet slip.
User confirms manually.

URL pattern discovered via recon:
  https://backend.wagerzon.com/wager/CreateWager.aspx?sel={sel_type}_{idgm}_{line}_{odds}&WT=0&lg={league_id}

API response structure (key fields):
  game.idgm     = internal game ID (used in sel parameter)
  game.idgp     = parent game group ID
  game.idlg     = league ID
  game.vnum     = away rotation number
  game.hnum     = home rotation number
  game.vtm      = away team name (Wagerzon format)
  game.htm      = home team name (Wagerzon format)
  child.idgm    = derivative line ID (1H, team totals, etc.)
  child.idgmtyp = type: 10=FG, 15=1H, 25=Alt, 35=TT(FG), 66=TT(1H)
"""

import os
import re
import sys
import time
import requests
from pathlib import Path

from playwright.sync_api import sync_playwright

from base_navigator import BaseNavigator, parse_market, update_bet_status, _REPO_ROOT

# Add wagerzon_odds to path for config (use main repo, not worktree)
sys.path.insert(0, str(_REPO_ROOT / "wagerzon_odds"))
from config import WAGERZON_BASE_URL, WAGERZON_HELPER_URL, SPORTS

WAGERZON_USERNAME = os.getenv("WAGERZON_USERNAME")
WAGERZON_PASSWORD = os.getenv("WAGERZON_PASSWORD")

# Map market periods to Wagerzon child game types
PERIOD_TO_GMTYP = {
    "fg": 10,   # Full game (parent)
    "h1": 15,   # First half
}

# Map market types to which line fields to use for the sel parameter
MARKET_SEL_FIELDS = {
    "spreads": {"away": ("vsprdt", "vsprdoddst"), "home": ("hsprdt", "hsprdoddst")},
    "totals": {"Over": ("ovt", "ovoddst"), "Under": ("unt", "unoddst")},
    "h2h": {"away": ("voddst", None), "home": ("hoddst", None)},
    "team_totals": {"Over": ("ovt", "ovoddst"), "Under": ("unt", "unoddst")},
    "alternate_spreads": {"away": ("vsprdt", "vsprdoddst"), "home": ("hsprdt", "hsprdoddst")},
    "alternate_totals": {"Over": ("ovt", "ovoddst"), "Under": ("unt", "unoddst")},
}


class WagerzonNavigator(BaseNavigator):

    BOOK_NAME = "wagerzon"

    def place_bet(self, bet_data: dict):
        """Place a single bet. Delegates to place_bets()."""
        self.place_bets([bet_data])

    def place_bets(self, bets: list[dict]):
        """Place multiple bets in a single Wagerzon browser session.

        Strategy: Navigate to the schedule page, click odds cells in the DOM
        to add each bet to the betslip, then fill amounts. All bets end up
        in one betslip for the user to confirm.
        """
        print(f"\n  Wagerzon Batch: {len(bets)} bet(s)", flush=True)
        for i, bet in enumerate(bets):
            print(f"  [{i+1}] {self._format_bet_summary(bet)}", flush=True)

        # Pre-fetch game data via API to get rotation numbers and team names
        api_cache = {}  # (away, home) -> api_game_data
        try:
            session = self._api_login()
            url = f"{WAGERZON_HELPER_URL}?WT=0&{SPORTS['cbb']['url_params']}"
            resp = session.get(url, timeout=30, headers={
                "Accept": "application/json, text/plain, */*",
                "X-Requested-With": "XMLHttpRequest",
            })
            api_data = resp.json()
            leagues_wrapper = api_data.get("result", {}).get("listLeagues", [[]])[0]

            for league in (leagues_wrapper or []):
                for game in league.get("Games", []):
                    api_cache[(game.get("vtm", "").upper(), game.get("htm", "").upper())] = game
        except Exception as e:
            print(f"  API pre-fetch failed: {e}")

        with sync_playwright() as p:
            browser = p.chromium.launch(
                headless=False,
                channel="chrome",
            )
            context = browser.new_context(
                viewport={"width": 1920, "height": 1080},
            )
            page = context.new_page()

            self._login(page)

            # Navigate to CBB schedule page
            schedule_url = f"{WAGERZON_BASE_URL}/wager/NewSchedule.aspx?WT=0&{SPORTS['cbb']['url_params']}"
            print(f"  Navigating to schedule page...", flush=True)
            page.goto(schedule_url)
            page.wait_for_load_state("networkidle")

            # Wait for game data to render — the page fetches via XHR
            # Look for elements that indicate games have loaded
            for wait_attempt in range(10):
                time.sleep(2)
                # Check if any game-like content appeared
                links = page.locator("a")
                link_count = links.count()
                if link_count > 20:  # Schedule page should have many links
                    print(f"  Schedule loaded ({link_count} links found)", flush=True)
                    break
                if wait_attempt == 9:
                    print(f"  Warning: schedule may not have loaded fully ({link_count} links)", flush=True)
                    # Dump page title and URL for debugging
                    print(f"  Page URL: {page.url}", flush=True)
                    try:
                        print(f"  Page title: {page.title()}", flush=True)
                    except Exception:
                        pass

            # Click odds cells on the schedule page to add to betslip
            success_count = 0
            for i, bet in enumerate(bets):
                parsed = parse_market(bet["market"])
                print(f"\n  --- Bet {i+1}/{len(bets)} ---", flush=True)
                print(f"  {self._format_bet_summary(bet)}", flush=True)

                try:
                    clicked = self._click_odds_on_schedule(page, bet, parsed, api_cache)
                    if clicked:
                        time.sleep(1)
                        success_count += 1
                        if bet.get("bet_hash"):
                            update_bet_status(bet["bet_hash"], "ready_to_confirm")
                    else:
                        self._print_manual_instructions(bet, parsed)
                        if bet.get("bet_hash"):
                            update_bet_status(bet["bet_hash"], "nav_error")
                except Exception as e:
                    print(f"  Failed to add bet: {e}", flush=True)
                    if bet.get("bet_hash"):
                        update_bet_status(bet["bet_hash"], "nav_error")

            # Click Continue to navigate to CreateWager page with all selected bets
            if success_count > 0:
                time.sleep(1)
                print(f"\n  Clicking Continue to open betslip...", flush=True)
                try:
                    continue_btn = page.locator("input.btnContinue")
                    if continue_btn.is_visible(timeout=5000):
                        continue_btn.click()
                        page.wait_for_load_state("networkidle")
                        time.sleep(3)
                        print(f"  CreateWager page: {page.url}", flush=True)

                        # Fill amounts on the CreateWager page
                        self._fill_betslip_amounts(page, bets)
                    else:
                        print(f"  Continue button not found", flush=True)
                except Exception as e:
                    print(f"  Error clicking Continue: {e}", flush=True)

            print(f"\n  {success_count}/{len(bets)} bets added to betslip.", flush=True)
            print("  Review amounts and confirm on the book's site.", flush=True)
            print("  Browser will stay open for 10 minutes.", flush=True)
            try:
                input()
            except EOFError:
                time.sleep(600)

            browser.close()

    def _click_odds_on_schedule(self, page, bet_data: dict, parsed, api_cache: dict) -> bool:
        """Find and click the correct odds button on the Wagerzon schedule page.

        DOM structure (discovered via inspection):
        - Team names: <span class="Team">KANSAS STATE</span>
        - Odds buttons: <div class="btn btn-odds ...">+11-110</div>
          - onclick handler, cursor:pointer
          - Spread format: +11-110 (signedLine + price)
          - Total format: o169-110 / u169-110
          - ML format: +480 / -700
          - Half-points: ½ character (e.g. +1½-110, o151½-110)
        """
        away = bet_data["away_team"]
        home = bet_data["home_team"]
        bet_on = bet_data["bet_on"]
        line = bet_data.get("line")
        odds = bet_data.get("odds")

        # Build target text patterns
        target_texts = self._build_odds_text_patterns(parsed, bet_on, line, odds, away, home)
        print(f"  Looking for odds: {target_texts}", flush=True)

        # Strategy: search all btn-odds buttons on the page, score each one,
        # and verify the match is in the correct game by checking nearby team names.
        odds_buttons = page.locator("div.btn-odds")
        btn_count = odds_buttons.count()
        print(f"  Found {btn_count} odds buttons on page", flush=True)

        best = None
        best_score = 0

        for j in range(btn_count):
            btn = odds_buttons.nth(j)
            try:
                # Use evaluate to get text regardless of viewport visibility.
                # is_visible() fails for off-screen buttons on scrollable pages.
                text = btn.evaluate("el => el.textContent").strip()
            except Exception:
                continue

            if not text:
                continue

            score = self._score_odds_match(text, parsed, bet_on, line, odds, away, home)
            if score < 3:
                continue

            # Verify this button belongs to the right game by checking nearby team names
            # Walk up to find a container that also has the team names
            is_correct_game = page.evaluate("""(el) => {
                // Walk up to find a game-level container
                let container = el.parentElement;
                for (let i = 0; i < 10 && container; i++) {
                    const teamSpans = container.querySelectorAll('span.Team');
                    if (teamSpans.length >= 2) {
                        const teams = Array.from(teamSpans).map(s => s.textContent.trim().toUpperCase());
                        return teams;
                    }
                    container = container.parentElement;
                }
                return [];
            }""", btn.element_handle())

            # Check if our teams are in the container
            if is_correct_game:
                away_terms = self._search_terms(away)
                home_terms = self._search_terms(home)
                joined = " ".join(is_correct_game)
                # Use ALL search terms — Wagerzon uses short names (e.g., "TULANE")
                # while canonical names include mascots ("Tulane Green Wave").
                # The single-word term may be the only one that matches.
                has_away = any(t.upper() in joined for t in away_terms)
                has_home = any(t.upper() in joined for t in home_terms)
                if has_away and has_home:
                    score += 3  # Bonus for correct game
                else:
                    continue  # Wrong game, skip

            if score > best_score:
                best_score = score
                best = (btn, text)

        if best and best_score >= 3:
            btn, text = best
            btn.scroll_into_view_if_needed()
            time.sleep(0.5)
            btn.click()
            print(f"  Clicked odds: '{text}' (score={best_score})", flush=True)
            return True

        print(f"  No matching odds button found (best_score={best_score})", flush=True)
        # Debug: show all odds buttons
        for j in range(min(btn_count, 20)):
            try:
                text = odds_buttons.nth(j).inner_text(timeout=200).strip()
                if text:
                    print(f"    [{j}] '{text}'", flush=True)
            except Exception:
                pass
        return False

    def _build_odds_text_patterns(self, parsed, bet_on, line, odds, away, home) -> list[str]:
        """Build expected Wagerzon odds text patterns.

        Wagerzon format: +11-110, o169-110, u151½-110, +480
        """
        patterns = []
        odds_str = str(odds) if odds else ""
        if odds and int(odds) > 0:
            odds_str = f"+{odds}"

        if line is not None:
            abs_line = abs(line)
            # Integer lines
            if abs_line == int(abs_line):
                line_int = str(int(abs_line))
                line_half = ""
            else:
                line_int = str(int(abs_line))
                line_half = f"{line_int}½"

            if parsed.market_type in ("spreads", "alternate_spreads"):
                sign = "+" if line >= 0 else "-"
                patterns.append(f"{sign}{line_int}{odds_str}")
                if line_half:
                    patterns.append(f"{sign}{line_half}{odds_str}")

            elif parsed.market_type in ("totals", "alternate_totals"):
                prefix = "o" if bet_on == "Over" else "u"
                patterns.append(f"{prefix}{line_int}{odds_str}")
                if line_half:
                    patterns.append(f"{prefix}{line_half}{odds_str}")

        if parsed.market_type == "h2h":
            patterns.append(odds_str)

        return patterns

    def _score_odds_match(self, text: str, parsed, bet_on, line, odds, away, home) -> int:
        """Score how well a link text matches our target bet."""
        score = 0
        # inner_text() returns newlines between line and odds (e.g., "+4\n-110")
        text = text.replace("\n", "")
        text_lower = text.lower()

        odds_str = str(odds) if odds else ""
        if odds and int(odds) > 0:
            odds_str = f"+{odds}"

        if line is not None:
            abs_line = abs(line)
            line_int = str(int(abs_line))
            line_half = f"{line_int}½" if abs_line != int(abs_line) else ""

            if parsed.market_type in ("spreads", "alternate_spreads"):
                sign = "+" if line >= 0 else "-"
                target = f"{sign}{line_half or line_int}{odds_str}"
                if target.lower() in text_lower:
                    score += 5

            elif parsed.market_type in ("totals", "alternate_totals"):
                prefix = "o" if bet_on == "Over" else "u"
                target = f"{prefix}{line_half or line_int}{odds_str}"
                if target.lower() in text_lower:
                    score += 5

            # Partial matches
            if line_int in text:
                score += 1
            if line_half and line_half in text:
                score += 1

        # Check odds presence
        if odds_str and odds_str in text:
            score += 2

        # For totals: check o/u prefix
        if parsed.market_type in ("totals", "alternate_totals"):
            if bet_on == "Over" and text_lower.startswith("o"):
                score += 1
            elif bet_on == "Under" and text_lower.startswith("u"):
                score += 1

        return score

    @staticmethod
    def _search_terms(team_name: str) -> list[str]:
        """Generate search terms from team name.

        Canonical names include mascots (e.g. "Kansas State Wildcats").
        Wagerzon uses short names (e.g. "KANSAS STATE").
        Generate progressively shorter prefixes, stripping mascots first.
        """
        words = team_name.split()
        terms = []

        # Common mascot words to strip
        mascots = {"wildcats", "cougars", "hokies", "deacons", "demon",
                   "bulldogs", "eagles", "bears", "tigers", "hawks", "huskies",
                   "cardinals", "knights", "rams", "wolves", "panthers",
                   "terriers", "gaels", "warriors", "broncos", "bobcats"}

        # Strip trailing mascot words
        while words and words[-1].lower() in mascots:
            words = words[:-1]

        # Generate prefixes from longest to shortest
        for n in range(len(words), 0, -1):
            prefix = " ".join(words[:n])
            if len(prefix) >= 3:
                terms.append(prefix)

        return terms if terms else [team_name]

    def _fill_betslip_amounts(self, page, bets: list[dict]):
        """Fill risk amounts on Wagerzon's CreateWager page.

        The CreateWager page has a checkbox to toggle between same amount
        for all bets vs individual amounts per selection. When bets have
        different sizes, we uncheck "same amount" to get per-bet inputs.
        """
        print(f"\n  Filling betslip amounts...", flush=True)
        time.sleep(1)

        # Check if bets have different sizes
        sizes = [int(b.get("recommended_size", 0)) for b in bets if b.get("recommended_size")]
        all_same = len(set(sizes)) <= 1

        # If different sizes, look for and uncheck "same amount" checkbox
        if not all_same:
            try:
                checkbox = page.locator("input[type='checkbox']")
                if checkbox.count() > 0 and checkbox.first.is_visible(timeout=2000):
                    # If checked, uncheck it to enable per-bet amounts
                    if checkbox.first.is_checked():
                        checkbox.first.uncheck()
                        print(f"  Unchecked 'same amount' for individual sizing", flush=True)
                        time.sleep(1)
            except Exception as e:
                print(f"  Could not toggle checkbox: {e}", flush=True)

        # Now find amount inputs — try multiple selectors since unchecking
        # the "same amount" checkbox can change input types (number → text)
        amount_inputs = None
        count = 0
        for sel in [
            "input[type='number'][placeholder='Amount']",
            "input[type='number']",
            "input[placeholder='Amount']",
            "input[placeholder='Amount ']",  # trailing space variant
            "input[name*='Amount' i]",
            "input[id*='Amount' i]",
        ]:
            candidate = page.locator(sel)
            c = candidate.count()
            if c > 0:
                amount_inputs = candidate
                count = c
                print(f"  Found {c} amount input(s) via: {sel}", flush=True)
                break

        if count == 0:
            # Wagerzon's CreateWager page uses type=text inputs with no
            # placeholder/name attributes after unchecking "same amount".
            # Find all text inputs (excluding checkbox/submit) — those are
            # the per-bet amount fields.
            text_inputs = page.locator(
                "input:visible[type='text']:not([type='checkbox']):not([type='submit'])"
            )
            tc = text_inputs.count()
            print(f"  Found {tc} text inputs as amount candidates", flush=True)

            if tc >= len(bets):
                # Fill per bet
                for i, bet in enumerate(bets):
                    size = int(bet.get("recommended_size", 0))
                    if not size or i >= tc:
                        continue
                    try:
                        field = text_inputs.nth(i)
                        field.evaluate(f"""el => {{
                            el.click();
                            el.focus();
                            el.select();
                            el.value = '';
                            document.execCommand('insertText', false, '{size}');
                            el.dispatchEvent(new Event('input', {{bubbles: true}}));
                            el.dispatchEvent(new Event('change', {{bubbles: true}}));
                        }}""")
                        print(f"  Filled bet {i+1}: ${size}", flush=True)
                    except Exception as e:
                        print(f"  Could not fill bet {i+1}: {e}", flush=True)
                return
            elif tc == 1:
                # Single input — fill with first bet's size
                size = sizes[0] if sizes else 0
                if size:
                    field = text_inputs.first
                    field.evaluate(f"""el => {{
                        el.click(); el.focus(); el.value = '';
                        document.execCommand('insertText', false, '{size}');
                        el.dispatchEvent(new Event('input', {{bubbles: true}}));
                        el.dispatchEvent(new Event('change', {{bubbles: true}}));
                    }}""")
                    print(f"  Filled single input: ${size}", flush=True)
                return

            print(f"  No amount inputs found on CreateWager page.", flush=True)
            # Diagnostic
            try:
                all_inputs = page.locator("input:visible")
                ic = all_inputs.count()
                print(f"  Visible inputs ({ic}):", flush=True)
                for k in range(min(ic, 10)):
                    attrs = page.evaluate("""el => ({
                        type: el.type, placeholder: el.placeholder, id: el.id, name: el.name
                    })""", all_inputs.nth(k).element_handle())
                    print(f"    [{k}] type={attrs.get('type')} placeholder='{attrs.get('placeholder')}'", flush=True)
            except Exception:
                pass
            return

        print(f"  Found {count} amount input(s)", flush=True)

        if count == 1 and all_same:
            # Single amount, same for all
            if sizes:
                field = amount_inputs.first
                field.click()
                time.sleep(0.3)
                page.keyboard.press("Meta+a")
                page.keyboard.press("Backspace")
                time.sleep(0.2)
                page.keyboard.type(str(sizes[0]), delay=80)
                page.keyboard.press("Tab")
                print(f"  Filled amount: ${sizes[0]}", flush=True)
        elif count == 1 and not all_same:
            # Still only one input even after unchecking — fill with largest
            max_size = max(sizes) if sizes else 0
            if max_size:
                field = amount_inputs.first
                field.click()
                time.sleep(0.3)
                page.keyboard.press("Meta+a")
                page.keyboard.press("Backspace")
                time.sleep(0.2)
                page.keyboard.type(str(max_size), delay=80)
                page.keyboard.press("Tab")
                print(f"  Filled amount: ${max_size} (single input, different sizes requested)", flush=True)
                print(f"  Adjust individual amounts manually: {['$'+str(s) for s in sizes]}", flush=True)
        else:
            # Multiple amount inputs — fill per bet
            for i, bet in enumerate(bets):
                size = int(bet.get("recommended_size", 0))
                if not size or i >= count:
                    continue
                try:
                    field = amount_inputs.nth(i)
                    if field.is_visible(timeout=1000):
                        field.click()
                        time.sleep(0.3)
                        page.keyboard.press("Meta+a")
                        page.keyboard.press("Backspace")
                        time.sleep(0.2)
                        page.keyboard.type(str(size), delay=80)
                        page.keyboard.press("Tab")
                        print(f"  Filled bet {i+1}: ${size}", flush=True)
                except Exception as e:
                    print(f"  Could not fill bet {i+1}: {e}", flush=True)

    def _api_login(self) -> requests.Session:
        """Login to Wagerzon via REST API (for game lookup)."""
        session = requests.Session()
        resp = session.get(WAGERZON_BASE_URL, timeout=15)
        resp.raise_for_status()

        if "NewSchedule" in resp.url:
            return session

        html = resp.text
        fields = {}
        for name in ["__VIEWSTATE", "__VIEWSTATEGENERATOR", "__EVENTVALIDATION",
                      "__EVENTTARGET", "__EVENTARGUMENT"]:
            match = re.search(rf'(?:name|id)="{name}"[^>]*value="([^"]*)"', html)
            if match:
                fields[name] = match.group(1)

        fields["Account"] = WAGERZON_USERNAME
        fields["Password"] = WAGERZON_PASSWORD
        fields["BtnSubmit"] = ""

        resp = session.post(WAGERZON_BASE_URL, data=fields, timeout=15)
        resp.raise_for_status()
        return session

    def _build_wager_url(self, bet_data: dict, parsed) -> tuple[str | None, int | None]:
        """Find the internal game ID via API and construct CreateWager URL.

        Returns (url, league_id) or (None, None) if game not found.
        """
        print("  Fetching game data from Wagerzon API...")

        try:
            session = self._api_login()
        except Exception as e:
            print(f"  API login failed: {e}")
            return None, None

        # Fetch CBB odds
        url = f"{WAGERZON_HELPER_URL}?WT=0&{SPORTS['cbb']['url_params']}"
        try:
            resp = session.get(url, timeout=30, headers={
                "Accept": "application/json, text/plain, */*",
                "X-Requested-With": "XMLHttpRequest",
            })
            data = resp.json()
        except Exception as e:
            print(f"  API fetch failed: {e}")
            return None, None

        leagues_wrapper = data.get("result", {}).get("listLeagues", [[]])[0]
        if not leagues_wrapper:
            print("  No leagues in API response")
            return None, None

        # Search through all leagues for matching game
        away_target = bet_data["away_team"].upper()
        home_target = bet_data["home_team"].upper()
        bet_on = bet_data["bet_on"]
        line = bet_data.get("line")
        odds = bet_data.get("odds")

        for league in leagues_wrapper:
            for game in league.get("Games", []):
                away_raw = game.get("vtm", "").upper()
                home_raw = game.get("htm", "").upper()

                # Match by team name (partial match — API uses abbreviated names)
                if not (away_target in away_raw or away_raw in away_target or
                        home_target in home_raw or home_raw in home_target):
                    continue

                # Check both directions (sometimes away/home are flipped)
                if not (home_target in home_raw or home_raw in home_target):
                    continue

                print(f"  Found game: {away_raw} @ {home_raw} (idgm={game.get('idgm')})")

                # Find the right line (parent or derivative)
                target_gmtyp = PERIOD_TO_GMTYP.get(parsed.period, 10)

                if target_gmtyp == 10:
                    # Full game — use parent
                    return self._construct_url(
                        game, game.get("GameLines", [{}])[0],
                        parsed, bet_data, game.get("idlg")
                    )
                else:
                    # Derivative — search GameChilds
                    for child in game.get("GameChilds", []):
                        if child.get("idgmtyp") == target_gmtyp:
                            # For team totals, match the team name
                            if parsed.market_type in ("team_totals", "alternate_team_totals"):
                                child_team = child.get("vtm", "").upper()
                                target_team = bet_on.upper() if bet_on not in ("Over", "Under") else ""
                                if parsed.side == "home":
                                    target_team = home_target
                                elif parsed.side == "away":
                                    target_team = away_target

                                if target_team and target_team not in child_team:
                                    continue

                            print(f"  Found derivative: idgm={child.get('idgm')} (type={target_gmtyp})")
                            return self._construct_url(
                                child, child.get("GameLines", [{}])[0],
                                parsed, bet_data, child.get("idlg")
                            )

                    print(f"  Game found but no matching derivative (type={target_gmtyp})")
                    # Fall back to parent game
                    return self._construct_url(
                        game, game.get("GameLines", [{}])[0],
                        parsed, bet_data, game.get("idlg")
                    )

        print(f"  Game not found in API response: {bet_data['away_team']} @ {bet_data['home_team']}")
        return None, None

    def _construct_url(self, game: dict, game_line: dict, parsed, bet_data: dict,
                       league_id: int | None) -> tuple[str | None, int | None]:
        """Build the CreateWager.aspx URL from game data.

        URL format: CreateWager.aspx?sel={sel_type}_{idgm}_{line}_{odds}&WT=0&lg={league_id}
        """
        idgm = game.get("idgm")
        if not idgm:
            return None, None

        bet_on = bet_data["bet_on"]
        line = bet_data.get("line")
        odds = bet_data.get("odds")

        # Determine sel_type (1 seems to be the standard for straight bets)
        sel_type = 1

        # Build the sel value: {sel_type}_{idgm}_{line}_{odds}
        line_str = str(line) if line is not None else "0"
        odds_str = str(odds) if odds else "-110"

        sel = f"{sel_type}_{idgm}_{line_str}_{odds_str}"
        lg = league_id or 43  # Default to CBB league

        url = f"{WAGERZON_BASE_URL}/wager/CreateWager.aspx?sel={sel}&WT=0&lg={lg}"
        print(f"  Constructed URL: {url}")
        return url, lg

    def _login(self, page):
        """Log in to Wagerzon via browser."""
        print("  Logging in to Wagerzon...")
        page.goto(WAGERZON_BASE_URL)
        page.wait_for_load_state("networkidle")

        if "NewSchedule" in page.url or "History" in page.url or "Create" in page.url:
            print("  Already authenticated (session restored)")
            return

        try:
            username_field = page.locator(
                'input[type="text"], input[name="username"], input[name="txtUsername"]'
            ).first
            username_field.fill(WAGERZON_USERNAME, timeout=10000)

            password_field = page.locator(
                'input[type="password"], input[name="password"], input[name="txtPassword"]'
            ).first
            password_field.fill(WAGERZON_PASSWORD)

            login_button = page.locator(
                'input[type="submit"], button[type="submit"], '
                'input[value="Login"], button:has-text("Login")'
            ).first
            login_button.click()

            page.wait_for_load_state("networkidle")
            time.sleep(2)
            print(f"  Login complete. URL: {page.url}")

        except Exception as e:
            print(f"  Auto-login failed: {e}")
            print("  Please log in manually. You have 60 seconds.")
            try:
                page.wait_for_url("**/wager/**", timeout=60000)
            except Exception:
                print("  Continuing anyway...")

    def _fill_amount(self, page, bet_data: dict):
        """Try to fill the bet amount on the wager page."""
        size = bet_data.get("recommended_size", 0)
        if not size:
            return

        print(f"  Looking for amount input field...")
        time.sleep(2)

        # Try common selectors for bet amount input
        for selector in [
            "input[name*='Risk']", "input[name*='risk']",
            "input[name*='Amount']", "input[name*='amount']",
            "input[name*='Stake']", "input[name*='stake']",
            "input[name*='Wager']", "input[name*='wager']",
            "input[type='number']",
            "input[type='text'][id*='risk']",
            "input[type='text'][id*='amount']",
        ]:
            try:
                field = page.locator(selector).first
                if field.is_visible(timeout=2000):
                    field.fill(str(int(size)))
                    print(f"  Filled amount: ${int(size)} (selector: {selector})")
                    return
            except Exception:
                continue

        print(f"  Could not find amount input. Enter ${int(size)} manually.")

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
