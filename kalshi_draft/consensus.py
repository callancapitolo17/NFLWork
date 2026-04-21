"""
Mock Draft Consensus Scraper + Comparison to Kalshi Odds.
Scrapes Tankathon NFL Draft big board (static HTML, reliable).
"""

import re
import sys
from datetime import datetime
from pathlib import Path

import requests
from bs4 import BeautifulSoup

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))
from nfl_draft.lib import db as nfl_db

import db


CONSENSUS_URL = "https://www.tankathon.com/nfl/big_board"


def scrape_consensus_board():
    """Scrape the big board from Tankathon.

    HTML structure: .mock-row.nfl containers, each with:
      - .mock-row-pick-number (rank)
      - .mock-row-player (name + "POS | School" concatenated)
    """
    print(f"Scraping consensus board from {CONSENSUS_URL}...")

    headers = {
        "User-Agent": "Mozilla/5.0 (Macintosh; Intel Mac OS X 10_15_7) AppleWebKit/537.36"
    }

    try:
        resp = requests.get(CONSENSUS_URL, headers=headers, timeout=30)
        resp.raise_for_status()
    except Exception as e:
        print(f"  Failed to fetch consensus board: {e}")
        return []

    soup = BeautifulSoup(resp.text, "html.parser")
    players = []

    # Each prospect is in a .mock-row.nfl div
    rows = soup.select(".mock-row.nfl")

    for row in rows:
        try:
            # Rank from .mock-row-pick-number
            rank_el = row.select_one(".mock-row-pick-number")
            if not rank_el:
                continue
            rank_text = rank_el.get_text(strip=True)
            rank_match = re.search(r"\d+", rank_text)
            if not rank_match:
                continue
            rank = int(rank_match.group())

            # Player info from .mock-row-player
            player_el = row.select_one(".mock-row-player")
            if not player_el:
                continue

            raw = player_el.get_text(strip=True)
            # Format: "Arvell ReeseLB | Ohio State"
            # Split on known position abbreviations
            match = re.match(
                r"^(.+?)(QB|RB|WR|TE|OT|IOL|DL|EDGE|LB|CB|S|K|P|LS|EDGE/LB)\s*\|\s*(.+)$",
                raw
            )
            if match:
                name = match.group(1).strip()
                position = match.group(2).strip()
                school = match.group(3).strip()
            else:
                # Fallback: try splitting on pipe
                parts = raw.split("|")
                if len(parts) >= 2:
                    name_pos = parts[0].strip()
                    school = parts[1].strip()
                    # Try to separate name from position
                    pos_match = re.search(r"(QB|RB|WR|TE|OT|IOL|DL|EDGE|LB|CB|S)$", name_pos)
                    if pos_match:
                        position = pos_match.group(1)
                        name = name_pos[:pos_match.start()].strip()
                    else:
                        name = name_pos
                        position = ""
                else:
                    name = raw
                    position = ""
                    school = ""

            if rank and name and len(name) > 2:
                players.append({
                    "rank": rank,
                    "player_name": name,
                    "position": position,
                    "school": school,
                })

        except (ValueError, AttributeError):
            continue

    # Deduplicate by rank
    seen_ranks = set()
    unique_players = []
    for p in players:
        if p["rank"] not in seen_ranks:
            seen_ranks.add(p["rank"])
            unique_players.append(p)

    unique_players.sort(key=lambda x: x["rank"])
    print(f"  Scraped {len(unique_players)} players from consensus board")
    return unique_players


def fuzzy_match_name(name1, name2):
    """Check if two player names match (handles first/last name variations)."""
    def normalize(n):
        n = n.lower().strip()
        n = re.sub(r"[^a-z\s]", "", n)
        n = re.sub(r"\s+(jr|sr|ii|iii|iv)\.?$", "", n)
        return n.split()

    parts1 = normalize(name1)
    parts2 = normalize(name2)

    if not parts1 or not parts2:
        return False

    # Exact match after normalization
    if parts1 == parts2:
        return True

    # Last name match + first initial
    if parts1[-1] == parts2[-1]:
        if parts1[0][0] == parts2[0][0]:
            return True

    # First and last name match (ignoring middle)
    if parts1[0] == parts2[0] and parts1[-1] == parts2[-1]:
        return True

    return False


def save_consensus(players):
    """Save consensus board to DuckDB."""
    if not players:
        return

    con = db.get_connection()
    fetch_time = datetime.now().strftime("%Y-%m-%d %H:%M:%S")

    con.executemany("""
        INSERT INTO consensus_board VALUES (?, ?, ?, ?, ?, ?)
    """, [
        (fetch_time, p["rank"], p["player_name"], p["position"], p["school"], "tankathon")
        for p in players
    ])

    con.close()
    print(f"  Saved {len(players)} consensus entries")


def compare_to_market():
    """Compare consensus rankings to Kalshi odds and print disagreements."""
    consensus = db.get_latest_consensus()
    latest = db.get_latest_odds()

    if consensus is None or consensus.empty or latest is None or latest.empty:
        print("Need both consensus and market data for comparison")
        return

    draft1 = latest[latest["series_ticker"] == "KXNFLDRAFT1"]
    if draft1.empty:
        print("No #1 pick series data")
        return

    print("\nConsensus vs Kalshi #1 Pick:")
    print("-" * 60)

    for _, player in consensus.head(20).iterrows():
        name = player["player_name"]
        rank = player["rank"]

        # Find in Kalshi
        kalshi_prob = None
        for _, market in draft1.iterrows():
            if fuzzy_match_name(name, market["candidate"]):
                kalshi_prob = market["last_price"]
                break

        if kalshi_prob is not None:
            marker = " ***" if (rank <= 5 and kalshi_prob < 5) or (rank > 10 and kalshi_prob > 10) else ""
            print(f"  #{rank:2d} {name:<25s} Kalshi: {kalshi_prob:3d}%{marker}")
        else:
            print(f"  #{rank:2d} {name:<25s} Kalshi: not found")


def run():
    """Full consensus pipeline."""
    nfl_db.init_schema()
    players = scrape_consensus_board()
    if players:
        save_consensus(players)
        compare_to_market()
    else:
        print("No consensus data scraped")


if __name__ == "__main__":
    run()
