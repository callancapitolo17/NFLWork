#!/usr/bin/env python3
"""
Hoop88 bet navigator.
Opens visible browser, navigates to CBB game, selects the correct odds cell,
and pre-fills the bet slip. User confirms manually.

Approach (same hybrid as BFA/Wagerzon navigators):
  1. API lookup: find game via REST API + verify line/odds
  2. Persistent browser profile: session cookies survive across runs
  3. Sidebar navigation: BASKETBALL → NCAA via [data-sport-sub-type] attributes
  4. DOM odds selection: .group-2/.group-3/.group-4 .line-play by market + row
  5. Bet slip filling: input[data-field="risk"] with event dispatching

Hoop88 DOM patterns (from hoop88_correlation/scraper_hoop88_odds.py):
  - Game panels: [data-panel="line"].GAME
  - Team names: [data-field="first-team"] (away), [data-field="second-team"] (home)
  - Line rows: .lines .line (first=away, second=home)
  - Odds groups: .group-2 (spread), .group-3 (ML), .group-4 (total)
  - Clickable cell: .line-play within each group
  - Betslip: .slide-up[data-panel="bets"], input[data-field="risk"]

REST API (JWT auth) for game lookup + line verification.
"""

import os
import time
import requests

from playwright.sync_api import sync_playwright

from base_navigator import (
    BaseNavigator, parse_market, update_bet_status, cleanup_singleton_lock,
    wait_for_confirmation, _REPO_ROOT,
)

HOOP88_URL = os.getenv("HOOP88_URL", "https://hoop88.com")
HOOP88_USERNAME = os.getenv("HOOP88_USERNAME")
HOOP88_PASSWORD = os.getenv("HOOP88_PASSWORD")

API_BASE = f"{HOOP88_URL}/cloud/api"

# Persistent browser profile directory
HOOP88_PROFILE_DIR = str(_REPO_ROOT / "bet_placer" / ".hoop88_nav_profile")

# Sidebar data-sport-sub-type values for CBB navigation
# NOTE: API uses sportSubType="NCAA" but the sidebar DOM uses "NCAA Basketball"
CBB_SPORT_TYPE = "BASKETBALL"
CBB_SPORT_SUB_TYPE = "NCAA Basketball"
CBB_API_SUB_TYPE = "NCAA"

# Period mapping: our codes → Hoop88 sidebar labels + API params
PERIOD_CONFIG = {
    "fg": {"label": CBB_SPORT_SUB_TYPE, "api_name": "Game", "api_num": 0},
    "h1": {"label": "1st Half", "api_name": "1st Half", "api_num": 1},
}

# Map market types to Hoop88 DOM group selectors
MARKET_GROUP = {
    "spreads": ".group-2",
    "alternate_spreads": ".group-2",
    "h2h": ".group-3",
    "totals": ".group-4",
    "alternate_totals": ".group-4",
    "team_totals": ".group-4",
    "alternate_team_totals": ".group-4",
}

# Major sidebar sections (used to scope period search within correct league)
MAJOR_SECTIONS = ["NFL", "NCAA Football", "NBA", "NHL", "NCAA Basketball"]


class Hoop88Navigator(BaseNavigator):

    BOOK_NAME = "hoop88"

    def place_bet(self, bet_data: dict):
        """Place a single bet. Delegates to place_bets()."""
        self.place_bets([bet_data])

    def place_bets(self, bets: list[dict], timeout: int = 300):
        """Place multiple bets in a single Hoop88 browser session.

        Groups bets by period to minimize sidebar switching.
        """
        from collections import defaultdict

        print(f"\n  Hoop88 Batch: {len(bets)} bet(s)", flush=True)
        for i, bet in enumerate(bets):
            print(f"  [{i+1}] {self._format_bet_summary(bet)}", flush=True)

        # Pre-fetch all game data via API (deduplicated by game)
        api_cache = {}  # (home, away) -> (api_game, hoop88_names)
        for bet in bets:
            key = (bet["home_team"], bet["away_team"])
            if key not in api_cache:
                parsed = parse_market(bet["market"])
                api_game, hoop88_names = self._lookup_game_via_api(*key)
                if api_game:
                    self._verify_line_from_api(api_game, bet, parsed)
                api_cache[key] = (api_game, hoop88_names)

        # Group bets by period for efficient sidebar navigation
        period_groups = defaultdict(list)
        for bet in bets:
            parsed = parse_market(bet["market"])
            period_groups[parsed.period].append(bet)

        with sync_playwright() as p:
            cleanup_singleton_lock(HOOP88_PROFILE_DIR)

            context = p.chromium.launch_persistent_context(
                user_data_dir=HOOP88_PROFILE_DIR,
                channel="chrome",
                headless=False,
                viewport={"width": 1920, "height": 1080},
                args=["--disable-blink-features=AutomationControlled"],
            )
            page = context.pages[0] if context.pages else context.new_page()

            page.goto(HOOP88_URL, wait_until="domcontentloaded")
            time.sleep(3)

            if not self._is_logged_in(page):
                self._login(page)

            self._navigate_to_sports(page)

            # Phase 1: Click all odds cells to add bets to betslip
            clicked_bets = []  # Track which bets were successfully clicked
            bet_idx = 0
            for period, period_bets in period_groups.items():
                self._navigate_to_cbb(page, period)

                for bet in period_bets:
                    bet_idx += 1
                    parsed = parse_market(bet["market"])
                    key = (bet["home_team"], bet["away_team"])
                    _, hoop88_names = api_cache.get(key, (None, []))

                    print(f"\n  --- Bet {bet_idx}/{len(bets)} ---", flush=True)
                    print(f"  {self._format_bet_summary(bet)}", flush=True)

                    try:
                        page.evaluate("window.scrollTo(0, 0)")
                        time.sleep(1)
                        game_panel = self._find_game_panel(page, bet, hoop88_names)
                        if game_panel:
                            clicked = self._click_odds_cell(page, game_panel, bet, parsed)
                            if clicked:
                                clicked_bets.append(bet)
                                time.sleep(2)
                            else:
                                self._print_manual_instructions(bet, parsed)
                                if bet.get("bet_hash"):
                                    update_bet_status(bet["bet_hash"], "nav_error")
                        else:
                            self._print_manual_instructions(bet, parsed)
                            if bet.get("bet_hash"):
                                update_bet_status(bet["bet_hash"], "nav_error")
                    except Exception as e:
                        print(f"  Failed to add bet: {e}", flush=True)
                        if bet.get("bet_hash"):
                            update_bet_status(bet["bet_hash"], "nav_error")

            # Phase 2: Fill all amounts AFTER all bets are in the betslip.
            # The SPA re-renders the betslip when a new bet is added, which
            # wipes previously set values. Fill all at once at the end.
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

    def _lookup_game_via_api(self, home_team: str, away_team: str) -> tuple[dict | None, list[str]]:
        """Look up game via Hoop88 REST API and return (game_dict, [away_name, home_name]).

        Uses same API patterns as hoop88_odds/scraper.py.
        """
        print("  Fetching game data from Hoop88 API...")

        try:
            session = requests.Session()
            session.headers.update({
                "X-Requested-With": "XMLHttpRequest",
                "Content-Type": "application/x-www-form-urlencoded; charset=UTF-8",
            })

            # Step 1: Get landing page cookies (Cloudflare)
            session.get(HOOP88_URL, timeout=15)

            # Step 2: Authenticate
            token = self._api_login(session)
            if not token:
                return None, []
            session.headers["Authorization"] = f"Bearer {token}"

            # Step 3: Fetch CBB lines
            lines = self._api_fetch_lines(session)
            if not lines:
                print("  No CBB lines from Hoop88 API")
                return None, []

            # Step 4: Match game by team name
            home_lower = home_team.lower()
            away_lower = away_team.lower()

            for game in lines:
                if game.get("Status") != "O":
                    continue

                away_raw = (game.get("Team1ID") or "").strip()
                home_raw = (game.get("Team2ID") or "").strip()
                if not away_raw or not home_raw:
                    continue

                away_match = away_lower in away_raw.lower() or away_raw.lower() in away_lower
                home_match = home_lower in home_raw.lower() or home_raw.lower() in home_lower

                if away_match and home_match:
                    print(f"  Found game: {away_raw} @ {home_raw}")
                    return game, [away_raw, home_raw]

            print(f"  Could not find game via API: {away_team} @ {home_team}")
            return None, []

        except Exception as e:
            print(f"  Hoop88 API lookup failed: {e}")
            return None, []

    def _api_login(self, session: requests.Session) -> str | None:
        """Authenticate via REST API and return JWT token."""
        try:
            domain = HOOP88_URL.replace("https://", "").replace("http://", "")
            resp = session.post(f"{API_BASE}/System/authenticateCustomer", data={
                "customerID": HOOP88_USERNAME,
                "password": HOOP88_PASSWORD,
                "state": "true",
                "multiaccount": "1",
                "response_type": "code",
                "client_id": HOOP88_USERNAME,
                "domain": domain,
                "redirect_uri": domain,
                "operation": "authenticateCustomer",
                "RRO": "1",
            }, timeout=15)

            if resp.status_code == 204:
                print("  API login failed — bad credentials")
                return None
            resp.raise_for_status()
            token = resp.json().get("code")
            if token:
                print("  API authenticated")
            return token
        except Exception as e:
            print(f"  API login failed: {e}")
            return None

    def _api_fetch_lines(self, session: requests.Session) -> list:
        """Fetch CBB game lines from Get_LeagueLines2."""
        try:
            resp = session.post(f"{API_BASE}/Lines/Get_LeagueLines2", data={
                "sportType": CBB_SPORT_TYPE,
                "sportSubType": CBB_API_SUB_TYPE,
                "period": "Game",
                "periodNumber": "0",
                "propDescription": "Game",
                "wagerType": "Straight",
                "office": "COINDEVIL",
                "customerID": HOOP88_USERNAME,
                "hourFilter": "0",
                "keyword": "",
                "correlationID": "",
                "grouping": "",
                "periods": "0",
                "rotOrder": "0",
                "placeLateFlag": "false",
                "RRO": "1",
                "agentSite": "0",
                "operation": "Get_LeagueLines2",
            }, timeout=30)
            resp.raise_for_status()
            lines = resp.json().get("Lines", [])
            print(f"  API returned {len(lines)} CBB lines")
            return lines
        except Exception as e:
            print(f"  API fetch failed: {e}")
            return []

    def _verify_line_from_api(self, game: dict, bet_data: dict, parsed):
        """Check if the line/odds from the API still match our target."""
        target_line = bet_data.get("line")
        target_odds = bet_data.get("odds")
        bet_on = bet_data["bet_on"]

        if target_line is None and target_odds is None:
            return

        # Determine spread direction
        spread_val = game.get("Spread")
        favored = (game.get("FavoredTeamID") or "").strip()
        away_raw = (game.get("Team1ID") or "").strip()

        if spread_val is not None and parsed.market_type in ("spreads", "alternate_spreads"):
            is_betting_away = (bet_on.lower() in bet_data["away_team"].lower()
                               or bet_data["away_team"].lower() in bet_on.lower())
            if favored.lower() == away_raw.lower():
                api_away_spread = -abs(spread_val)
            else:
                api_away_spread = abs(spread_val)
            api_line = api_away_spread if is_betting_away else -api_away_spread
            api_odds = game.get("SpreadAdj1") if is_betting_away else game.get("SpreadAdj2")

            if target_line is not None and api_line is not None:
                if abs(api_line - target_line) > 0.01:
                    print(f"  WARNING: Line moved! Expected {target_line}, API has {api_line}")
                else:
                    print(f"  Line verified: {api_line} ({api_odds})")

        elif parsed.market_type in ("totals", "alternate_totals"):
            api_total = game.get("TotalPoints")
            is_over = bet_on == "Over"
            api_odds = game.get("TtlPtsAdj1") if is_over else game.get("TtlPtsAdj2")
            if target_line is not None and api_total is not None:
                if abs(api_total - target_line) > 0.01:
                    print(f"  WARNING: Total moved! Expected {target_line}, API has {api_total}")
                else:
                    print(f"  Line verified: {api_total} ({api_odds})")

        elif parsed.market_type == "h2h":
            is_betting_away = (bet_on.lower() in bet_data["away_team"].lower()
                               or bet_data["away_team"].lower() in bet_on.lower())
            api_odds = game.get("MoneyLine1") if is_betting_away else game.get("MoneyLine2")
            if target_odds is not None and api_odds is not None:
                if api_odds != target_odds:
                    print(f"  WARNING: ML moved! Expected {target_odds}, API has {api_odds}")
                else:
                    print(f"  ML verified: {api_odds}")

    # =========================================================================
    # NAVIGATION
    # =========================================================================

    def _navigate_to_sports(self, page):
        """Navigate to the Hoop88 sports page."""
        sports_url = f"{HOOP88_URL}/sports.html"
        print(f"  Navigating to {sports_url}...")
        page.goto(sports_url, wait_until="domcontentloaded")

        # Wait for sidebar to render
        for i in range(20):
            time.sleep(1)
            try:
                if (page.locator("[data-sport-sub-type]").count() > 0
                        or page.locator('[data-panel="line"]').count() > 0):
                    print(f"  Sports page rendered ({i + 1}s)")
                    return
            except Exception:
                pass

        print("  Sports page may not have fully rendered — continuing")

    def _is_logged_in(self, page) -> bool:
        """Check if user is already authenticated."""
        try:
            # Login form visible = NOT logged in
            login_form = page.locator('input[name="customerID"]')
            if login_form.count() > 0 and login_form.first.is_visible(timeout=3000):
                return False
        except Exception:
            pass
        # No login form = logged in (persistent profile)
        print("  Already authenticated (session active)")
        return True

    def _login(self, page):
        """Log in to Hoop88 via browser form.

        With persistent profile, this only runs on first use or session expiry.
        """
        print("  Logging in to Hoop88...")

        try:
            page.fill('input[name="customerID"]', HOOP88_USERNAME)
            page.fill('input[name="Password"]', HOOP88_PASSWORD)
            page.click('button[data-action="login"]')

            # Wait for login form to disappear
            page.locator('input[name="customerID"]').wait_for(state="hidden", timeout=15000)
            time.sleep(2)
            print("  Login complete")

        except Exception as e:
            print(f"  Auto-login failed: {e}")
            print("  Please log in manually. You have 60 seconds.")
            try:
                page.locator('input[name="customerID"]').wait_for(state="hidden", timeout=60000)
            except Exception:
                print("  Continuing anyway...")

    def _navigate_to_cbb(self, page, period: str):
        """Navigate to CBB via sidebar: BASKETBALL → NCAA [→ period].

        Clicking BASKETBALL expands the section and auto-loads the default
        league's games. Then we click the NCAA sub-type to switch to CBB.
        If NCAA isn't in the sidebar (no CBB games today), we stay on the
        default view.

        For derivative periods (1H), we walk the sidebar DOM to find the
        correct "1st Half" that belongs to the NCAA section (same approach
        as hoop88_correlation/scraper_hoop88_odds.py navigate_to_period).
        """
        print("  Navigating to CBB section via sidebar...")

        # Click BASKETBALL to expand (auto-loads default league games)
        try:
            basketball = page.locator("text=BASKETBALL").first
            if basketball.is_visible(timeout=5000):
                basketball.click()
                time.sleep(3)
                print("  Expanded BASKETBALL section")
        except Exception:
            print("  BASKETBALL section may already be expanded")

        # Get the period config
        period_cfg = PERIOD_CONFIG.get(period, PERIOD_CONFIG["fg"])
        target_label = period_cfg["label"]

        if target_label == CBB_SPORT_SUB_TYPE:
            # Full game: click NCAA Basketball sub-type via JS (matches reference scraper)
            clicked = page.evaluate(f'''() => {{
                const elements = document.querySelectorAll('[data-sport-sub-type="{CBB_SPORT_SUB_TYPE}"]');
                for (const elem of elements) {{
                    elem.click();
                    return true;
                }}
                return false;
            }}''')
            if clicked:
                time.sleep(3)
                print(f"  Selected {CBB_SPORT_SUB_TYPE} (full game)")
                return
            print(f"  {CBB_SPORT_SUB_TYPE} not in sidebar — CBB may not have games today")
        else:
            # Derivative period: find the correct "1st Half" after NCAA in sidebar
            clicked = page.evaluate(f'''() => {{
                const allElements = Array.from(document.querySelectorAll('[data-sport-sub-type]'));
                const leagueSelector = "{CBB_SPORT_SUB_TYPE}";
                const periodSelector = "{target_label}";
                const majorSections = {MAJOR_SECTIONS};

                let leagueIdx = -1;
                for (let i = 0; i < allElements.length; i++) {{
                    if (allElements[i].getAttribute('data-sport-sub-type') === leagueSelector) {{
                        leagueIdx = i;
                        break;
                    }}
                }}

                if (leagueIdx === -1) return false;

                for (let i = leagueIdx + 1; i < allElements.length; i++) {{
                    const subType = allElements[i].getAttribute('data-sport-sub-type');
                    if (majorSections.includes(subType) && subType !== leagueSelector) break;
                    if (subType === periodSelector) {{
                        allElements[i].click();
                        return true;
                    }}
                }}
                return false;
            }}''')

            if clicked:
                time.sleep(3)
                print(f"  Selected period: {target_label}")
                return

            print(f"  Period '{target_label}' not found under NCAA")

    # =========================================================================
    # MARKET SELECTION + BET SLIP
    # =========================================================================

    def _find_game_panel(self, page, bet_data: dict, hoop88_names: list[str] = None):
        """Find the [data-panel="line"].GAME panel containing the target game.

        Searches by team name in [data-field="first-team"] (away) and
        [data-field="second-team"] (home).
        """
        away = bet_data["away_team"]
        home = bet_data["home_team"]

        # Build search pairs: API names first (most reliable), then dashboard names
        name_pairs = []
        if hoop88_names and len(hoop88_names) >= 2:
            name_pairs.append((hoop88_names[0], hoop88_names[1]))
        name_pairs.append((away, home))

        # Wait for game panels to load
        try:
            page.locator('[data-panel="line"].GAME').first.wait_for(timeout=10000)
        except Exception:
            print("  No game panels found on page")
            return None

        panels = page.locator('[data-panel="line"].GAME')
        count = panels.count()
        print(f"  Found {count} game panels on page")

        for away_n, home_n in name_pairs:
            for i in range(count):
                panel = panels.nth(i)
                try:
                    first_team = panel.locator('[data-field="first-team"]').inner_text().strip()
                    second_team = panel.locator('[data-field="second-team"]').inner_text().strip()

                    away_match = self._fuzzy_team_match(away_n, first_team)
                    home_match = self._fuzzy_team_match(home_n, second_team)

                    if away_match and home_match:
                        panel.scroll_into_view_if_needed()
                        print(f"  Found game: {first_team} @ {second_team}")
                        return panel
                except Exception:
                    continue

        # Also try matching by just the first word of each team name
        for i in range(count):
            panel = panels.nth(i)
            try:
                first_team = panel.locator('[data-field="first-team"]').inner_text().strip()
                second_team = panel.locator('[data-field="second-team"]').inner_text().strip()

                # Last resort: first-word match for both teams
                away_first = away.split()[0].lower()
                home_first = home.split()[0].lower()
                if (away_first in first_team.lower() and home_first in second_team.lower()):
                    panel.scroll_into_view_if_needed()
                    print(f"  Found game (first-word match): {first_team} @ {second_team}")
                    return panel
            except Exception:
                continue

        print(f"  Could not find game on page: {away} @ {home}", flush=True)
        # Diagnostic: show what teams are on the page
        try:
            for i in range(min(count, 15)):
                panel = panels.nth(i)
                ft = panel.locator('[data-field="first-team"]').inner_text(timeout=500).strip()
                st = panel.locator('[data-field="second-team"]').inner_text(timeout=500).strip()
                print(f"    Panel {i}: {ft} @ {st}", flush=True)
        except Exception:
            pass
        return None

    def _click_odds_cell(self, page, panel, bet_data: dict, parsed) -> bool:
        """Click the correct odds cell within a game panel.

        Layout: .lines .line (first=away, second=home)
        Groups: .group-2 (spread), .group-3 (ML), .group-4 (total)
        Clickable: .line-play within the group
        """
        bet_on = bet_data["bet_on"]
        group_sel = MARKET_GROUP.get(parsed.market_type, ".group-2")

        # Determine which row: first .line = away (Team1), second .line = home (Team2)
        if parsed.market_type in ("totals", "alternate_totals"):
            row_idx = 0 if bet_on == "Over" else 1
        elif parsed.market_type in ("team_totals", "alternate_team_totals"):
            row_idx = 0 if parsed.side == "away" else 1
        else:
            is_away = (bet_on.lower() in bet_data["away_team"].lower()
                       or bet_data["away_team"].lower() in bet_on.lower())
            row_idx = 0 if is_away else 1

        rows = panel.locator(".lines .line")
        try:
            if rows.count() <= row_idx:
                print(f"  Not enough line rows (found {rows.count()})")
                return False

            row = rows.nth(row_idx)
            cell = row.locator(f"{group_sel} .line-play")

            if cell.count() > 0 and cell.first.is_visible(timeout=3000):
                cell_text = cell.first.inner_text().strip()
                print(f"  Odds cell text: '{cell_text}'")
                cell.first.click()
                print(f"  Clicked odds cell (row={row_idx}, group={group_sel})")
                return True

        except Exception as e:
            print(f"  Error clicking odds cell: {e}")

        return False

    def _fill_all_amounts(self, page, bets: list[dict]):
        """Fill amounts for all bets at once, after all odds are clicked.

        Finds all visible risk inputs in the betslip and fills them in order.
        Must be called AFTER all bets are added — SPA re-renders clear values.
        """
        import time as _time
        sizes = [int(b.get("recommended_size", 0)) for b in bets]

        # Find visible risk inputs in betslip
        selector = '.slide-up[data-panel="bets"] input[data-field="risk"]:not([disabled])'
        fields = page.locator(selector)
        count = fields.count()

        # Collect visible inputs in order (top to bottom by Y position)
        visible_inputs = []
        for idx in range(count):
            candidate = fields.nth(idx)
            try:
                rect = candidate.evaluate("el => { const r = el.getBoundingClientRect(); return {t: r.top, w: r.width}; }")
                if rect["w"] > 0:
                    visible_inputs.append((rect["t"], idx, candidate))
            except Exception:
                continue

        # Sort by Y position (top to bottom = bet 1, bet 2, ...)
        visible_inputs.sort(key=lambda x: x[0])

        print(f"  Found {len(visible_inputs)} visible risk inputs for {len(bets)} bets", flush=True)

        for i, (y_pos, idx, field) in enumerate(visible_inputs):
            if i >= len(sizes) or not sizes[i]:
                continue
            val_str = str(sizes[i])
            try:
                field.evaluate(f"""el => {{
                    el.scrollIntoView();
                    el.click();
                    el.focus();
                    el.select();
                    el.value = '';
                    document.execCommand('insertText', false, '{val_str}');
                    el.dispatchEvent(new Event('input', {{bubbles: true}}));
                    el.dispatchEvent(new Event('change', {{bubbles: true}}));
                }}""")
                print(f"  Bet {i+1}: ${sizes[i]} filled (y={y_pos:.0f})", flush=True)
                if bets[i].get("bet_hash"):
                    update_bet_status(bets[i]["bet_hash"], "ready_to_confirm")
                _time.sleep(0.5)  # Brief pause between fills
            except Exception as e:
                print(f"  Bet {i+1}: failed to fill ${sizes[i]}: {e}", flush=True)

    # Common state/word abbreviations used by sportsbooks
    _ABBREVS = {
        "va": "virginia", "n": "north", "s": "south", "e": "east", "w": "west",
        "st": "state", "jax": "jacksonville", "fla": "florida", "la": "louisiana",
        "ark": "arkansas", "miss": "mississippi", "tenn": "tennessee",
        "conn": "connecticut", "mass": "massachusetts", "mich": "michigan",
        "minn": "minnesota", "wis": "wisconsin", "ill": "illinois",
        "ind": "indiana", "penn": "pennsylvania", "md": "maryland",
        "ala": "alabama", "ga": "georgia", "colo": "colorado",
        "ore": "oregon", "wash": "washington", "okla": "oklahoma",
        "neb": "nebraska", "tex": "texas",
        "nc": "north carolina", "sc": "south carolina",
    }

    @staticmethod
    def _fuzzy_team_match(canonical: str, page_name: str) -> bool:
        """Fuzzy match team names. Handles abbreviations like 'VA Tech' vs 'Virginia Tech'."""
        c = canonical.lower()
        p = page_name.lower()

        # Direct substring match (either direction)
        if c in p or p in c:
            return True

        # First-word match (e.g. "Wake" in "Wake Forest")
        c_first = c.split()[0]
        if c_first in p:
            return True

        # Expand abbreviations in page name, then re-check
        p_words = p.split()
        expanded = []
        for w in p_words:
            expanded.append(Hoop88Navigator._ABBREVS.get(w, w))
        expanded_name = " ".join(expanded)
        if expanded_name in c or c in expanded_name:
            return True

        # Check if page words (expanded) are prefixes of canonical words
        c_words = c.split()
        if len(p_words) <= len(c_words):
            match_count = 0
            for pw in expanded:
                for cw in c_words:
                    if cw.startswith(pw) or pw.startswith(cw):
                        match_count += 1
                        break
            if match_count == len(p_words):
                return True

        return False

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
