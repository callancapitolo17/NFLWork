#!/usr/bin/env python3
"""
Base navigator class for bet auto-placement.
Shared logic: market parsing, game matching via DuckDB, login helpers.
"""

import os
import re
import duckdb
from pathlib import Path
from dataclasses import dataclass
from typing import Optional

from dotenv import load_dotenv

# Resolve the main repo root (handles both direct runs and worktree runs).
# Worktrees live under .claude/worktrees/<name>/ — detect and resolve to actual repo.
_THIS_DIR = Path(__file__).parent.resolve()
_REPO_ROOT = _THIS_DIR.parent
if ".claude/worktrees" in str(_REPO_ROOT):
    # Running inside a worktree — resolve to the real repo root
    _REPO_ROOT = Path(str(_REPO_ROOT).split(".claude/worktrees")[0])

# Load credentials from bet_logger/.env
ENV_PATH = _REPO_ROOT / "bet_logger" / ".env"
load_dotenv(ENV_PATH)

# DuckDB paths for game matching (always in the main repo, not worktree)
WAGERZON_DB = _REPO_ROOT / "wagerzon_odds" / "wagerzon.duckdb"
HOOP88_DB = _REPO_ROOT / "hoop88_odds" / "hoop88.duckdb"
BFA_DB = _REPO_ROOT / "bfa_odds" / "bfa.duckdb"
DASHBOARD_DB = _REPO_ROOT / "Answer Keys" / "CBB Dashboard" / "cbb_dashboard.duckdb"


@dataclass
class ParsedMarket:
    """Decomposed market identifier."""
    market_type: str   # spreads, totals, h2h, team_totals, alternate_spreads, alternate_totals
    period: str        # h1, h2, fg, q1, q2, q3, q4
    side: Optional[str] = None  # away, home (for team_totals only)


def parse_market(market: str) -> ParsedMarket:
    """Parse market string into components.

    Examples:
        "spreads_h1"              -> type=spreads, period=h1
        "totals_fg"               -> type=totals, period=fg
        "h2h_h1"                  -> type=h2h, period=h1
        "team_totals_away_h1"     -> type=team_totals, period=h1, side=away
        "team_totals_home_h2"     -> type=team_totals, period=h2, side=home
        "alternate_spreads_h1"    -> type=alternate_spreads, period=h1
        "alternate_totals_fg"     -> type=alternate_totals, period=fg
        "alternate_team_totals_away_h1" -> type=alternate_team_totals, period=h1, side=away
    """
    # Team totals with side (home/away) embedded
    m = re.match(r'^(alternate_team_totals|team_totals)_(away|home)_(\w+)$', market)
    if m:
        return ParsedMarket(market_type=m.group(1), period=m.group(3), side=m.group(2))

    # Alternate markets (alternate_spreads, alternate_totals)
    m = re.match(r'^(alternate_spreads|alternate_totals)_(\w+)$', market)
    if m:
        return ParsedMarket(market_type=m.group(1), period=m.group(2))

    # Standard markets (spreads, totals, h2h)
    m = re.match(r'^(\w+?)_(\w+)$', market)
    if m:
        return ParsedMarket(market_type=m.group(1), period=m.group(2))

    return ParsedMarket(market_type=market, period="fg")


def lookup_game(bookmaker: str, home_team: str, away_team: str) -> Optional[dict]:
    """Look up a book-specific game record by canonical team names.

    Returns dict with game_id and other book-specific fields, or None.
    """
    db_map = {
        "wagerzon": WAGERZON_DB,
        "hoop88": HOOP88_DB,
        "bfa": BFA_DB,
    }
    db_path = db_map.get(bookmaker)
    if not db_path or not db_path.exists():
        return None

    con = None
    try:
        con = duckdb.connect(str(db_path), read_only=True)
        row = con.execute("""
            SELECT DISTINCT game_id, game_date, game_time, away_team, home_team
            FROM cbb_odds
            WHERE home_team = ? AND away_team = ?
            LIMIT 1
        """, [home_team, away_team]).fetchone()

        if row:
            return {
                "game_id": row[0],
                "game_date": row[1],
                "game_time": row[2],
                "away_team": row[3],
                "home_team": row[4],
            }
    except Exception as e:
        print(f"  Game lookup failed for {bookmaker}: {e}")
    finally:
        if con:
            con.close()

    return None


def update_bet_status(bet_hash: str, status: str):
    """Update placed_bets.status in the dashboard DuckDB (with retry for lock contention)."""
    import time as _time
    for attempt in range(5):
        con = None
        try:
            con = duckdb.connect(str(DASHBOARD_DB))
            con.execute(
                "UPDATE placed_bets SET status = ? WHERE bet_hash = ?",
                [status, bet_hash],
            )
            return
        except Exception as e:
            if attempt < 4 and "lock" in str(e).lower():
                _time.sleep(0.5 * (attempt + 1))
            else:
                print(f"  Status update failed: {e}")
        finally:
            if con:
                con.close()


class BaseNavigator:
    """Base class for book-specific navigators."""

    BOOK_NAME = "unknown"

    def place_bets(self, bets: list[dict]):
        """Place multiple bets in one browser session. Override for batch betslip.

        Default: falls back to place_bet() per bet (N sessions).
        Subclasses should override for true batch support (1 session, all bets in betslip).
        """
        for bet in bets:
            self.place_bet(bet)

    def place_bet(self, bet_data: dict):
        """Launch browser, navigate to book, pre-fill bet. Override in subclass."""
        raise NotImplementedError

    def _format_bet_summary(self, bet_data: dict) -> str:
        """Human-readable summary of the bet being placed."""
        parsed = parse_market(bet_data["market"])
        period_label = {
            "h1": "1H", "h2": "2H", "fg": "FG",
            "q1": "1Q", "q2": "2Q", "q3": "3Q", "q4": "4Q",
        }.get(parsed.period, parsed.period.upper())

        type_label = {
            "spreads": "Spread",
            "totals": "Total",
            "h2h": "ML",
            "team_totals": "Team Total",
            "alternate_spreads": "Alt Spread",
            "alternate_totals": "Alt Total",
            "alternate_team_totals": "Alt Team Total",
        }.get(parsed.market_type, parsed.market_type)

        line_str = ""
        if bet_data.get("line") is not None:
            line_str = f" {bet_data['line']}"

        odds = bet_data.get("odds", "")
        try:
            if odds and int(odds) > 0:
                odds = f"+{odds}"
        except (ValueError, TypeError):
            pass

        size = bet_data.get("recommended_size", "?")

        return (
            f"{bet_data['away_team']} @ {bet_data['home_team']} | "
            f"{period_label} {type_label} | "
            f"{bet_data['bet_on']}{line_str} ({odds}) | "
            f"${size}"
        )
