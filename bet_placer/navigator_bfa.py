#!/usr/bin/env python3
"""
BFA Gaming bet navigator.
Opens visible browser, navigates to CBB game, selects the correct odds cell,
and pre-fills the bet slip. User confirms manually.

Approach (inspired by Wagerzon navigator):
  1. API lookup: find BFA game ID + verify line via public REST API (no auth)
  2. Persistent browser profile: Keycloak session cookies survive across runs
  3. Direct URL navigation: go to /league?slug=ncaa_b (skip sidebar clicking)
  4. DOM odds selection: click the correct odds cell using column index
  5. Bet slip filling: enter amount in the betslip input

BFA is a Blazor Server SPA (MudBlazor components, SignalR WebSocket).
Key DOM patterns:
  - Game rows: article.fixture-view-component > article.threeopt-fixture-component
  - Team name: span.team-name-value
  - Odds cells: article.odd-button > button.mud-button-root
  - 3 cells per team row: spread(0), moneyline(1), total(2)
  - "More wagers" link in div.more-wagers-link-content

Public REST API (no auth) for game lookup + line verification.
"""

import os
import shutil
import time
import requests

from playwright.sync_api import sync_playwright

from base_navigator import (
    BaseNavigator, parse_market, update_bet_status, cleanup_singleton_lock,
    wait_for_confirmation, _REPO_ROOT,
)

BFA_USERNAME = os.getenv("BFA_USERNAME")
BFA_PASSWORD = os.getenv("BFA_PASSWORD")

BFA_BASE_URL = "https://bfagaming.com"
BFA_API_BASE = "https://api.bfagaming.com"
BFA_API_PARAMS = {"playerId": 0, "agentId": 0, "fixtureType": 0}
BFA_API_HEADERS = {
    "Authorization": "Bearer",
    "User-Agent": "bfa-client/1.0.0",
    "Referer": "https://bfagaming.com/",
}
BFA_CBB_LEAGUE_ID = 493

# Persistent browser profile directory (shared with recon_bfa.py pattern)
BFA_PROFILE_DIR = str(_REPO_ROOT / "bet_placer" / ".bfa_nav_profile")

# Map our period codes to BFA API period names
PERIOD_DISPLAY = {
    "fg": "Game",
    "h1": "1st Half",
    "h2": "2nd Half",
}

# Map our market types to BFA market type IDs
MARKET_TYPE_ID = {
    "spreads": 2,
    "alternate_spreads": 2,
    "totals": 3,
    "alternate_totals": 3,
    "h2h": 1,
    "team_totals": None,  # 4=home TT, 5=away TT
    "alternate_team_totals": None,
}

# Column index for each market type in the 3-cell odds row
MARKET_COLUMN = {
    "spreads": 0,
    "alternate_spreads": 0,
    "h2h": 1,
    "totals": 2,
    "alternate_totals": 2,
    "team_totals": 2,
    "alternate_team_totals": 2,
}


class BFANavigator(BaseNavigator):

    BOOK_NAME = "bfa"

    def place_bet(self, bet_data: dict):
        """Place a single bet. Delegates to place_bets()."""
        self.place_bets([bet_data])

    def place_bets(self, bets: list[dict], timeout: int = 300):
        """Place multiple bets in a single BFA browser session."""
        print(f"\n  BFA Batch: {len(bets)} bet(s)", flush=True)
        for i, bet in enumerate(bets):
            print(f"  [{i+1}] {self._format_bet_summary(bet)}", flush=True)

        # Pre-fetch all game data via API (deduplicated by game)
        api_cache = {}  # (home, away) -> (game_id, detail, bfa_names)
        for bet in bets:
            key = (bet["home_team"], bet["away_team"])
            if key not in api_cache:
                parsed = parse_market(bet["market"])
                game_id, detail, bfa_names = self._lookup_game_via_api(*key)
                if detail:
                    self._verify_line_from_api(detail, bet, parsed)
                api_cache[key] = (game_id, detail, bfa_names)

        with sync_playwright() as p:
            cleanup_singleton_lock(BFA_PROFILE_DIR)

            context = p.chromium.launch_persistent_context(
                user_data_dir=BFA_PROFILE_DIR,
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

            self._navigate_to_league(page)

            if not self._is_logged_in(page):
                self._login(page)
                self._navigate_to_league(page)

            # Verify page actually loaded content (stale sessions pass
            # _is_logged_in but don't render game data)
            time.sleep(3)  # Give Blazor time to render
            team_spans = page.locator("span.team-name-value")
            if team_spans.count() == 0:
                print("  Page content didn't load — clearing cookies and re-login...", flush=True)
                context.clear_cookies()
                self._login(page)
                self._navigate_to_league(page)
                time.sleep(3)

            # Clear any stale bets from previous runs
            self._clear_betslip(page)

            # Track whether we've left the league page (e.g., after "More wagers")
            on_league_page = True

            # Phase 1: Click all odds cells to add bets to betslip
            clicked_bets = []
            for i, bet in enumerate(bets):
                parsed = parse_market(bet["market"])
                key = (bet["home_team"], bet["away_team"])
                _, detail, bfa_names = api_cache.get(key, (None, None, []))

                print(f"\n  --- Bet {i+1}/{len(bets)} ---", flush=True)
                print(f"  {self._format_bet_summary(bet)}", flush=True)

                # If we navigated away (expanded view), go back to league page
                if not on_league_page:
                    self._navigate_to_league(page)
                    on_league_page = True

                try:
                    needs_expand = parsed.period != "fg"
                    clicked = self._select_odds_cell(page, bet, parsed, bfa_names)
                    if needs_expand:
                        on_league_page = False
                    if clicked:
                        clicked_bets.append(bet)
                        time.sleep(2)
                    else:
                        self._print_manual_instructions(bet, parsed)
                        if bet.get("bet_hash"):
                            update_bet_status(bet["bet_hash"], "nav_error")
                except Exception as e:
                    print(f"  Failed to add bet: {e}", flush=True)
                    if bet.get("bet_hash"):
                        update_bet_status(bet["bet_hash"], "nav_error")

            # Phase 2: Fill all amounts AFTER all odds are clicked.
            # SPA betslips re-render when new bets are added, wiping values.
            if clicked_bets:
                print(f"\n  Filling amounts for {len(clicked_bets)} bets...", flush=True)
                time.sleep(2)
                self._fill_all_amounts(page, clicked_bets)

            print(f"\n  {len(clicked_bets)}/{len(bets)} bets added to betslip.", flush=True)
            print("  Confirm manually on the book's site.", flush=True)
            print(f"  Browser will close after {timeout}s or when bets are confirmed.", flush=True)

            confirmed_hashes = [b["bet_hash"] for b in clicked_bets if b.get("bet_hash")]
            if confirmed_hashes:
                wait_for_confirmation(confirmed_hashes, timeout=timeout)

            context.close()

    # =========================================================================
    # API GAME LOOKUP
    # =========================================================================

    def _lookup_game_via_api(self, home_team: str, away_team: str) -> tuple[int | None, dict | None, list[str]]:
        """Look up BFA internal game ID and fetch full game detail.

        Returns (game_id, game_detail, bfa_names) where bfa_names is
        [away_name, home_name] from the API (for DOM matching).
        """
        print("  Fetching game data from BFA API...")

        try:
            session = requests.Session()
            session.headers.update(BFA_API_HEADERS)

            games = self._fetch_game_list(session)
            if not games:
                print("  No games found via BFA API")
                return None, None, []

            home_lower = home_team.lower()
            away_lower = away_team.lower()

            for game in games:
                if game.get("type") != 1:
                    continue
                fixtures = game.get("fixtures", [])
                if not fixtures:
                    continue
                contestants = fixtures[0].get("contestants", [])
                if len(contestants) < 2:
                    continue

                names = [c.get("name", "").lower() for c in contestants]
                home_match = any(home_lower in n or n in home_lower for n in names)
                away_match = any(away_lower in n or n in away_lower for n in names)

                if home_match and away_match:
                    game_id = game.get("id")
                    bfa_names = [c.get("name", "") for c in contestants]
                    print(f"  Found game: {bfa_names[0]} @ {bfa_names[1]} (id={game_id})")
                    detail = self._fetch_game_detail(session, game_id)
                    return game_id, detail, bfa_names

            print(f"  Could not find game via BFA API: {away_team} @ {home_team}")
            return None, None, []

        except Exception as e:
            print(f"  BFA API lookup failed: {e}")
            return None, None, []

    def _fetch_game_list(self, session: requests.Session) -> list:
        """Fetch game list from BFA API (popular + league fallback)."""
        resp = session.get(
            f"{BFA_API_BASE}/oddsservice/events/popular/ncaa_b",
            params={**BFA_API_PARAMS, "set": "Auto"},
            timeout=15,
        )
        if resp.status_code == 200:
            data = resp.json()
            games = [g for g in data.get("games", data if isinstance(data, list) else []) if isinstance(g, dict)]
            if games:
                return games

        resp = session.get(
            f"{BFA_API_BASE}/oddsservice/events/leagues",
            params={**BFA_API_PARAMS, "id": BFA_CBB_LEAGUE_ID, "set": "Auto"},
            timeout=15,
        )
        if resp.status_code == 200:
            data = resp.json()
            return [g for g in data.get("games", data if isinstance(data, list) else []) if isinstance(g, dict)]

        return []

    def _fetch_game_detail(self, session: requests.Session, game_id: int) -> dict | None:
        """Fetch full game detail (all markets, periods, alt lines)."""
        try:
            resp = session.get(
                f"{BFA_API_BASE}/oddsservice/event/{game_id}",
                params=BFA_API_PARAMS,
                timeout=15,
            )
            resp.raise_for_status()
            detail = resp.json()
            print(f"  Fetched game detail: {len(detail.get('markets', []))} markets")
            return detail
        except Exception as e:
            print(f"  Game detail fetch failed: {e}")
            return None

    def _verify_line_from_api(self, game_detail: dict, bet_data: dict, parsed):
        """Check if the line/odds we want still match the API data."""
        target_line = bet_data.get("line")
        target_odds = bet_data.get("odds")
        bet_on = bet_data["bet_on"]

        if target_line is None and target_odds is None:
            return

        fixture = game_detail.get("fixtures", [{}])[0]
        contestants = fixture.get("contestants", [])
        if len(contestants) < 2:
            return
        away_id = contestants[0]["id"]
        home_id = contestants[1]["id"]

        period_map = {p["name"]: p["number"] for p in game_detail.get("periods", [])}
        target_period = PERIOD_DISPLAY.get(parsed.period, "Game")
        period_num = period_map.get(target_period)

        market_type_id = MARKET_TYPE_ID.get(parsed.market_type)
        if parsed.market_type in ("team_totals", "alternate_team_totals"):
            market_type_id = 4 if parsed.side == "home" else 5

        if market_type_id is None or period_num is None:
            return

        target_market = next(
            (m for m in game_detail.get("markets", [])
             if m.get("type") == market_type_id and m.get("periodNumber") == period_num),
            None,
        )
        if not target_market:
            return

        is_alt = "alternate" in parsed.market_type
        for odd in target_market.get("odds", []):
            if odd.get("status", 0) != 0:
                continue
            if not is_alt and odd.get("index", 0) != 0:
                continue

            matched = False
            if parsed.market_type in ("spreads", "alternate_spreads", "h2h"):
                is_home = bet_on.lower() in bet_data["home_team"].lower() or bet_data["home_team"].lower() in bet_on.lower()
                if odd.get("contestantId") == (home_id if is_home else away_id):
                    matched = True
            elif parsed.market_type in ("totals", "alternate_totals", "team_totals", "alternate_team_totals"):
                if bet_on == "Over" and odd.get("side") == 4:
                    matched = True
                elif bet_on == "Under" and odd.get("side") == 5:
                    matched = True

            if not matched:
                continue
            if is_alt and target_line is not None:
                if abs(odd.get("line", 0) - target_line) > 0.01:
                    continue

            api_line = odd.get("line")
            api_odds = odd.get("price")
            line_ok = True
            if target_line is not None and api_line is not None:
                if abs(api_line - target_line) > 0.01:
                    print(f"  WARNING: Line moved! Expected {target_line}, API has {api_line}")
                    line_ok = False
            if target_odds is not None and api_odds is not None:
                if api_odds != target_odds:
                    print(f"  WARNING: Odds moved! Expected {target_odds}, API has {api_odds}")
                    line_ok = False
            if line_ok:
                print(f"  Line verified: {api_line} ({api_odds})")
            return

    # =========================================================================
    # NAVIGATION (direct URL — no sidebar clicking)
    # =========================================================================

    def _navigate_to_league(self, page):
        """Navigate directly to the NCAA(B) league page and select 'Games' tab.

        The league page has a scroller with tabs like 'Second half', 'Games', etc.
        Default is often 'Second half' which only shows live/in-progress games.
        Clicking 'Games' shows all upcoming games.
        """
        league_url = f"{BFA_BASE_URL}/league?slug=ncaa_b"
        print(f"  Navigating to {league_url}...")
        page.goto(league_url, wait_until="domcontentloaded")

        # Wait for Blazor to render game fixtures (poll for team name spans)
        for i in range(30):
            time.sleep(1)
            try:
                if page.locator("span.team-name-value").count() > 2:
                    print(f"  League page rendered ({i + 1}s)")
                    break
            except Exception:
                pass
        else:
            print("  League page may not have fully rendered — continuing")

        # Click 'Games' tab to show all upcoming games (default may be 'Second half')
        games_tab = page.locator("div.scroller-item:has-text('Games')")
        try:
            if games_tab.count() > 0 and games_tab.first.is_visible(timeout=3000):
                games_tab.first.click()
                time.sleep(3)
                print(f"  Selected 'Games' tab ({page.locator('span.team-name-value').count()} team spans)")
        except Exception:
            print("  'Games' tab not found — using default view")

    # =========================================================================
    # LOGIN (persistent profile handles most cases)
    # =========================================================================

    def _is_logged_in(self, page) -> bool:
        """Check if the user is already authenticated.

        IMPORTANT: Check negative signals FIRST. BFA shows "My Account" text
        in the header even when NOT logged in, so we must check for "Log In"
        button before checking positive indicators.
        """
        try:
            # Negative signals — check these FIRST (they're definitive)
            login_btn = page.locator('button:has-text("Log In")')
            if login_btn.count() > 0 and login_btn.first.is_visible(timeout=3000):
                print("  Not logged in (Log In button visible)")
                return False

            for text in ["Sign Up", "Register", "Join"]:
                if page.locator(f"text='{text}'").count() > 0:
                    print(f"  Not logged in ('{text}' visible)")
                    return False

            # Positive signals — only trust these after ruling out negatives
            for selector in [
                ".balance",
                ".user-balance",
                ".header-balance",
                "text='Deposit'",
            ]:
                if page.locator(selector).count() > 0:
                    print("  Already authenticated (session active)")
                    return True

        except Exception:
            pass

        # Conservative: if we can't confirm auth, assume NOT logged in
        print("  Auth status unclear — will attempt login")
        return False

    def _login(self, page):
        """Log in to BFA.

        BFA has two possible login flows:
          A. Keycloak SSO: redirects to auth.bfagaming.com
          B. Internal Blazor login: redirects to bfagaming.com/authentication/login

        With persistent profile, sessions may survive across runs.
        """
        print("  Logging in...")

        # Find and click the login button
        login_btn = None
        for sel in [
            'button:has-text("Log In")',
            'a:has-text("Log In")',
            'text="Log In"',
        ]:
            loc = page.locator(sel)
            if loc.count() > 0:
                login_btn = loc.first
                break

        if not login_btn:
            print("  Could not find Log In button")
            return

        try:
            login_btn.click(no_wait_after=True)

            # Wait for login page to load (either Keycloak or internal)
            for _ in range(10):
                time.sleep(1)
                url = page.url
                if "auth.bfagaming.com" in url or "/authentication/login" in url:
                    break

            url = page.url
            print(f"  Login page: {url}")

            if "auth.bfagaming.com" in url:
                # Keycloak SSO flow
                page.fill("#username", BFA_USERNAME, timeout=10000)
                page.fill("#password", BFA_PASSWORD)
                page.click("#kc-login")
                page.wait_for_url("**/bfagaming.com/**", timeout=30000)
                time.sleep(3)

            elif "/authentication/login" in url:
                # Internal Blazor login page — find username/password inputs
                time.sleep(2)  # Let Blazor render

                # Try common input patterns for Blazor/MudBlazor login forms
                filled = False
                for user_sel, pass_sel in [
                    ("input[type='text']", "input[type='password']"),
                    ("input[autocomplete='username']", "input[autocomplete='current-password']"),
                    ("#username", "#password"),
                    ("input.mud-input-slot[type='text']", "input.mud-input-slot[type='password']"),
                    ("input[placeholder*='user' i]", "input[placeholder*='pass' i]"),
                    ("input[placeholder*='email' i]", "input[placeholder*='pass' i]"),
                ]:
                    user_input = page.locator(user_sel)
                    pass_input = page.locator(pass_sel)
                    if user_input.count() > 0 and pass_input.count() > 0:
                        # Use mouse click + keyboard type for Blazor compatibility
                        user_input.first.click()
                        time.sleep(0.3)
                        page.keyboard.press("Meta+a")
                        page.keyboard.type(BFA_USERNAME, delay=50)
                        time.sleep(0.3)

                        pass_input.first.click()
                        time.sleep(0.3)
                        page.keyboard.press("Meta+a")
                        page.keyboard.type(BFA_PASSWORD, delay=50)
                        time.sleep(0.3)

                        # Submit: look for a login/submit button, or press Enter
                        submit = None
                        for btn_sel in [
                            'button:has-text("Log In")',
                            'button:has-text("Login")',
                            'button:has-text("Sign In")',
                            'button[type="submit"]',
                        ]:
                            loc = page.locator(btn_sel)
                            if loc.count() > 0:
                                submit = loc.first
                                break

                        if submit:
                            submit.click()
                        else:
                            page.keyboard.press("Enter")

                        filled = True
                        print(f"  Credentials entered (selectors: {user_sel})")
                        break

                if not filled:
                    print("  Could not find username/password inputs on login page")

                # Wait for redirect back to main site
                for _ in range(15):
                    time.sleep(1)
                    if "/authentication/login" not in page.url:
                        break

                time.sleep(3)

            # Verify login succeeded
            login_gone = page.locator('button:has-text("Log In")').count() == 0
            if login_gone:
                print("  Login complete")
            else:
                print("  Login may not have succeeded — please log in manually (60s)")
                for _ in range(60):
                    time.sleep(1)
                    if page.locator('button:has-text("Log In")').count() == 0:
                        print("  Manual login detected")
                        break

        except Exception as e:
            print(f"  Auto-login failed: {e}")
            print("  Please log in manually (60s).")
            try:
                for _ in range(60):
                    time.sleep(1)
                    if page.locator('button:has-text("Log In")').count() == 0:
                        print("  Manual login detected")
                        break
            except Exception:
                print("  Continuing anyway...")

    # =========================================================================
    # MARKET SELECTION + BET SLIP
    # =========================================================================

    def _select_odds_cell(self, page, bet_data: dict, parsed, bfa_names: list[str] = None) -> bool:
        """Find the game on the page and click the correct odds cell (without filling amount).

        Returns True if odds cell was successfully clicked, False otherwise.
        """
        game_fixture = self._find_game_fixture(page, bet_data, bfa_names)
        if not game_fixture:
            return False

        if parsed.period != "fg":
            self._click_more_wagers(page, game_fixture, bet_data)
            self._select_period_tab(page, parsed.period)
            return self._click_odds_in_expanded(page, bet_data, parsed)
        else:
            return self._click_odds_in_row(page, game_fixture, bet_data, parsed)

    def _find_game_fixture(self, page, bet_data: dict, bfa_names: list[str] = None):
        """Find the fixture-view-component containing the target game.

        Uses BFA API names (most reliable) then falls back to dashboard names.
        Scrolls down to trigger lazy loading if game isn't found initially.
        """
        away = bet_data["away_team"]
        home = bet_data["home_team"]

        # Build search pairs: BFA API names first (most reliable), then dashboard names
        name_pairs = []
        if bfa_names and len(bfa_names) >= 2:
            name_pairs.append((bfa_names[0], bfa_names[1]))
        name_pairs.append((away, home))

        # Try finding the game, scrolling down if needed
        for attempt in range(3):
            if attempt > 0:
                # Scroll down to trigger lazy loading
                print(f"  Scrolling down (attempt {attempt + 1})...")
                page.evaluate("window.scrollTo(0, document.body.scrollHeight)")
                time.sleep(2)

            result = self._search_fixtures(page, name_pairs)
            if result:
                return result

        print(f"  Could not find game on page: {away} @ {home}")
        return None

    def _search_fixtures(self, page, name_pairs: list[tuple[str, str]]):
        """Search visible fixtures for a matching game."""
        for away_n, home_n in name_pairs:
            for team in [away_n, home_n]:
                other = home_n if team == away_n else away_n
                for term in self._search_terms(team):
                    matches = page.locator(f"span.team-name-value:has-text('{term}')")
                    try:
                        count = matches.count()
                        for i in range(count):
                            el = matches.nth(i)
                            if not el.is_visible(timeout=2000):
                                continue
                            fixture = el.locator(
                                "xpath=ancestor::article[contains(@class, 'fixture-view-component')]"
                            )
                            if fixture.count() == 0:
                                continue
                            # Validate: other team must also be in this fixture
                            fixture_text = fixture.first.inner_text()
                            if self._team_in_text(other, fixture_text):
                                el.scroll_into_view_if_needed()
                                print(f"  Found game on page (matched '{term}')")
                                return fixture.first
                    except Exception:
                        continue
        return None

    @staticmethod
    def _search_terms(team_name: str) -> list[str]:
        """Generate progressively shorter search terms from front of name.
        "Oregon Ducks" -> ["Oregon Ducks", "Oregon"]
        "South Dakota State" -> ["South Dakota State", "South Dakota"]
        """
        words = team_name.split()
        return [" ".join(words[:n]) for n in range(len(words), 0, -1)
                if len(" ".join(words[:n])) >= 3]

    @staticmethod
    def _team_in_text(team_name: str, text: str) -> bool:
        """Check if any prefix of the team name appears in text."""
        text_lower = text.lower()
        words = team_name.split()
        for n in range(len(words), 0, -1):
            t = " ".join(words[:n]).lower()
            if len(t) >= 3 and t in text_lower:
                return True
        return False

    def _click_odds_in_row(self, page, fixture, bet_data: dict, parsed) -> bool:
        """Click the correct odds cell within a game fixture row.

        The main game list shows 3 odds cells per team: spread(0), ML(1), total(2).
        """
        bet_on = bet_data["bet_on"]
        column_idx = MARKET_COLUMN.get(parsed.market_type, 0)

        # Determine which team row (first fixture-entry = away, second = home)
        if parsed.market_type in ("totals", "alternate_totals"):
            team_row_idx = 0 if bet_on == "Over" else 1
        elif parsed.market_type in ("team_totals", "alternate_team_totals"):
            team_row_idx = 0 if parsed.side == "away" else 1
        else:
            is_away = (bet_on.lower() in bet_data["away_team"].lower()
                       or bet_data["away_team"].lower() in bet_on.lower())
            team_row_idx = 0 if is_away else 1

        entries = fixture.locator("div.fixture-entry")
        try:
            if entries.count() <= team_row_idx:
                print(f"  Not enough fixture entries (found {entries.count()})")
                return False

            entry = entries.nth(team_row_idx)
            odds_buttons = entry.locator("article.odd-button button")

            if odds_buttons.count() <= column_idx:
                print(f"  Not enough odds cells (found {odds_buttons.count()}, need index {column_idx})")
                return False

            target_btn = odds_buttons.nth(column_idx)
            if target_btn.is_visible(timeout=3000):
                btn_text = target_btn.inner_text()
                print(f"  Odds cell text: '{btn_text.strip()}'")
                target_btn.click()
                print(f"  Clicked odds cell (row={team_row_idx}, col={column_idx})")
                return True

        except Exception as e:
            print(f"  Error clicking odds cell: {e}")

        return False

    def _click_more_wagers(self, page, fixture, bet_data: dict):
        """Click 'More wagers' to expand the game to show derivative markets."""
        more = fixture.locator("text='More wagers'")
        try:
            if more.is_visible(timeout=3000):
                more.click()
                time.sleep(3)
                print(f"  Expanded game via 'More wagers'")
                return True
        except Exception:
            pass
        print(f"  'More wagers' not found — derivative markets may not be available")
        return False

    def _select_period_tab(self, page, period: str):
        """Try to click a period tab on the expanded game view."""
        period_labels = {
            "h1": ["1st Half", "1H", "First Half"],
            "h2": ["2nd Half", "2H", "Second Half"],
        }
        for label in period_labels.get(period, []):
            tab = page.locator(f"text=/{label}/i").first
            try:
                if tab.is_visible(timeout=3000):
                    tab.click()
                    time.sleep(2)
                    print(f"  Selected period: {label}")
                    return True
            except Exception:
                continue

        print(f"  Could not find period tab for '{period}'")
        return False

    def _expand_market_panel(self, page, parsed) -> bool:
        """Expand the correct MudBlazor accordion panel for the target market.

        BFA uses collapsible panels like "1st Half Alternate Total Points",
        "1st Half Spread", etc. Must click the panel header to reveal odds.
        """
        # Build search terms for the panel header text
        market_keywords = {
            "spreads": ["Spread"],
            "alternate_spreads": ["Alternate Spread", "Alt Spread"],
            "totals": ["Total Points", "Total"],
            "alternate_totals": ["Alternate Total", "Alt Total"],
            "h2h": ["Moneyline", "Money Line", "ML"],
            "team_totals": ["Team Total"],
            "alternate_team_totals": ["Alternate Team Total", "Alt Team Total"],
        }

        keywords = market_keywords.get(parsed.market_type, [])
        if not keywords:
            return False

        # Find and click the matching expand panel
        panels = page.locator("div.mud-expand-panel")
        try:
            panel_count = panels.count()
            print(f"  Found {panel_count} market panels")

            for i in range(panel_count):
                panel = panels.nth(i)
                try:
                    header_text = panel.locator("div.mud-expand-panel-text").inner_text(timeout=1000)
                except Exception:
                    continue

                for kw in keywords:
                    if kw.lower() in header_text.lower():
                        # Check if already expanded
                        panel_classes = panel.get_attribute("class") or ""
                        if "mud-panel-expanded" not in panel_classes:
                            panel.locator("div.mud-expand-panel-header").click()
                            time.sleep(1)
                            print(f"  Expanded panel: '{header_text}'")
                        else:
                            print(f"  Panel already expanded: '{header_text}'")
                        return True
        except Exception as e:
            print(f"  Error expanding market panel: {e}")

        print(f"  Could not find market panel for {parsed.market_type}")
        return False

    def _click_odds_in_expanded(self, page, bet_data: dict, parsed) -> bool:
        """Click odds cell in the expanded game view (after More wagers).

        BFA expanded view shows market groups in collapsible MudBlazor panels.
        Strategy:
          1. Expand the correct market panel (e.g., "1st Half Alternate Total Points")
          2. Search for buttons containing the line number
          3. Among matches, prefer one that also contains the odds
          4. Handle half-point display (65.5 vs 65½)
        """
        # First, expand the correct market panel
        self._expand_market_panel(page, parsed)
        time.sleep(1)

        line = bet_data.get("line")
        odds = bet_data.get("odds")
        bet_on = bet_data["bet_on"]

        # Build line search variants (65.5, +65.5, -65.5, O 65.5, U 65.5)
        line_variants = []
        if line is not None:
            abs_line = abs(line)
            line_plain = str(abs_line)
            # Handle integer lines (e.g., 4.0 -> "4")
            if abs_line == int(abs_line):
                line_plain = str(int(abs_line))
            line_variants.append(line_plain)
            # Half-point unicode variant (65.5 -> 65½)
            if abs_line != int(abs_line):
                line_variants.append(str(int(abs_line)) + "½")

        odds_str = ""
        if odds:
            odds_str = f"+{odds}" if int(odds) > 0 else str(odds)

        # Determine if we need Over/Under context
        is_over = bet_on == "Over" or (parsed.side == "away" and parsed.market_type in ("team_totals", "alternate_team_totals"))
        is_under = bet_on == "Under" or (parsed.side == "home" and parsed.market_type in ("team_totals", "alternate_team_totals"))

        # Search all visible odd-button buttons in the expanded panel
        all_buttons = page.locator("article.odd-button button")
        try:
            total = all_buttons.count()
            print(f"  Scanning {total} odds buttons in expanded view...")

            best_match = None
            best_score = 0

            for i in range(min(total, 50)):
                btn = all_buttons.nth(i)
                try:
                    if not btn.is_visible(timeout=500):
                        continue
                    btn_text = btn.inner_text().strip()
                except Exception:
                    continue

                # Check if any line variant appears in button text
                line_match = False
                for variant in line_variants:
                    if variant in btn_text:
                        line_match = True
                        break

                if not line_match:
                    continue

                score = 1  # matched line

                # Bonus: odds match
                if odds_str and odds_str in btn_text:
                    score += 2

                # Bonus/penalty: Over/Under direction
                btn_upper = btn_text.upper()
                if is_over and ("O " in btn_upper or "OVER" in btn_upper):
                    score += 1
                elif is_under and ("U " in btn_upper or "UNDER" in btn_upper):
                    score += 1
                elif is_over and ("U " in btn_upper or "UNDER" in btn_upper):
                    score -= 2  # wrong direction
                elif is_under and ("O " in btn_upper or "OVER" in btn_upper):
                    score -= 2

                if score > best_score:
                    best_score = score
                    best_match = (btn, btn_text)

            if best_match and best_score > 0:
                btn, btn_text = best_match
                btn.click()
                print(f"  Clicked expanded odds: '{btn_text}' (score={best_score})")
                return True

            print(f"  No matching odds button found (line_variants={line_variants}, odds={odds_str})")
        except Exception as e:
            print(f"  Error scanning expanded odds: {e}")

        return False

    def _clear_betslip(self, page):
        """Remove all existing selections from the betslip.

        BFA has a "Remove all selections" link in the betslip when bets exist.
        Must clear stale bets from previous runs before adding new ones.
        """
        try:
            remove_link = page.locator("text='Remove all selections'")
            if remove_link.count() > 0 and remove_link.first.is_visible(timeout=2000):
                remove_link.first.click()
                time.sleep(1)
                print("  Cleared stale bets from betslip")
            else:
                # Also try close/X buttons on individual bet cards
                close_btns = page.locator(".bet-slip-card .mud-icon-button, .remove-selection")
                if close_btns.count() > 0:
                    for i in range(close_btns.count()):
                        try:
                            close_btns.nth(0).click()  # Always click first (list shrinks)
                            time.sleep(0.3)
                        except Exception:
                            break
                    print(f"  Cleared {close_btns.count()} stale bet cards")
        except Exception as e:
            print(f"  Betslip clear attempt: {e}")

    def _fill_all_amounts(self, page, bets: list[dict]):
        """Fill amounts for all bets at once, after all odds are clicked.

        BFA betslip has "Straight" tab with individual bet cards. Each card has
        a "Risk" input (left) and "To Win" input (right) at the same Y position.
        We take the leftmost input per Y-row to get only Risk inputs.

        Must be called AFTER all bets are added — SPA re-renders clear values.
        BFA uses Blazor Server (SignalR) — must use real mouse/keyboard, not
        execCommand or manual event dispatch.
        """
        sizes = [int(b.get("recommended_size", 0)) for b in bets]

        # Make sure "Straight" tab is selected (not parlay)
        try:
            straight_tab = page.locator("text='Straight'")
            if straight_tab.count() > 0 and straight_tab.first.is_visible(timeout=2000):
                straight_tab.first.click()
                time.sleep(0.5)
        except Exception:
            pass

        # Find all visible number inputs, collect with position info
        all_number_inputs = page.locator("input[type='number']")
        count = all_number_inputs.count()

        visible_inputs = []
        for idx in range(count):
            candidate = all_number_inputs.nth(idx)
            try:
                info = candidate.evaluate("""el => {
                    const r = el.getBoundingClientRect();
                    let label = '';
                    let parent = el.parentElement;
                    for (let j = 0; j < 5 && parent; j++) {
                        const lbl = parent.querySelector('label');
                        if (lbl) { label = lbl.textContent.trim().toLowerCase(); break; }
                        parent = parent.parentElement;
                    }
                    return { t: r.top, l: r.left, w: r.width, label: label };
                }""")
                if info["w"] > 0:
                    visible_inputs.append((info["t"], info["l"], info["label"], idx, candidate))
            except Exception:
                continue

        visible_inputs.sort(key=lambda x: (x[0], x[1]))

        # Take leftmost input per Y-row (Risk is left, To Win is right)
        risk_inputs = []
        seen_rows = set()
        for y, x, label, idx, field in visible_inputs:
            row_key = round(y / 10)
            if label and any(kw in label for kw in ["risk", "wager", "stake", "amount"]):
                risk_inputs.append((y, idx, field))
                seen_rows.add(row_key)
            elif row_key not in seen_rows:
                risk_inputs.append((y, idx, field))
                seen_rows.add(row_key)

        print(f"  Found {len(risk_inputs)} risk inputs for {len(bets)} bets", flush=True)

        for i, (y_pos, idx, field) in enumerate(risk_inputs):
            if i >= len(sizes) or not sizes[i]:
                continue
            val_str = str(sizes[i])
            try:
                field.evaluate("el => el.scrollIntoView({block: 'center'})")
                time.sleep(0.8)

                rect = field.evaluate("""el => {
                    const r = el.getBoundingClientRect();
                    return { x: r.left + r.width / 2, y: r.top + r.height / 2 };
                }""")

                page.mouse.click(rect["x"], rect["y"], click_count=3)
                time.sleep(0.5)

                focused_ok = field.evaluate("el => document.activeElement === el")
                if not focused_ok:
                    page.mouse.click(rect["x"], rect["y"], click_count=3)
                    time.sleep(0.5)

                page.keyboard.type(val_str, delay=80)
                time.sleep(0.3)
                page.keyboard.press("Tab")
                print(f"  Bet {i+1}: ${sizes[i]} filled", flush=True)
                if bets[i].get("bet_hash"):
                    update_bet_status(bets[i]["bet_hash"], "ready_to_confirm")
                time.sleep(0.5)
            except Exception as e:
                print(f"  Bet {i+1}: failed to fill ${sizes[i]}: {e}", flush=True)

        if not risk_inputs:
            print(f"  No risk inputs found — enter amounts manually.", flush=True)

    def _print_manual_instructions(self, bet_data: dict, parsed):
        """Print instructions for manual placement."""
        bet_on = bet_data["bet_on"]
        line = bet_data.get("line")
        odds = bet_data.get("odds")
        size = bet_data.get("recommended_size", 0)

        print(f"\n  --- MANUAL STEP REQUIRED ---")
        print(f"  Game:   {bet_data['away_team']} @ {bet_data['home_team']}")
        print(f"  Select: {parsed.market_type.upper()} | Period: {parsed.period.upper()}")
        print(f"  Pick:   {bet_on}" + (f" {line}" if line else ""))
        print(f"  Odds:   {'+' if odds and odds > 0 else ''}{odds}")
        print(f"  Amount: ${size:.0f}")
        print(f"  ---")
