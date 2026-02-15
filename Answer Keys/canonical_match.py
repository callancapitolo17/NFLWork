"""
Shared team name resolution for offshore scrapers.
Two-layer approach: ESPN dictionary lookup + game-level matching fallback.
"""

import json
import re
import duckdb
from pathlib import Path

ANSWER_KEYS_DIR = Path(__file__).parent


def strip_punct(name):
    """Remove punctuation and normalize whitespace for matching."""
    return re.sub(r"['\.\-]", "", name).strip().lower()


def load_team_dict(sport="cbb"):
    """Load ESPN team dictionary from DuckDB.
    Returns {stripped_name: odds_api_name} with multiple keys per team.
    """
    if sport == "nfl":
        return {}

    db_path = ANSWER_KEYS_DIR / f"{sport}.duckdb"
    if not db_path.exists():
        return {}

    conn = duckdb.connect(str(db_path), read_only=True)
    try:
        tables = [r[0] for r in conn.execute("SHOW TABLES").fetchall()]
        if "cbb_team_dict" not in tables:
            return {}
        df = conn.execute(
            "SELECT short_name, nickname, abbreviation, odds_api_name "
            "FROM cbb_team_dict WHERE odds_api_name IS NOT NULL"
        ).fetchdf()
    except Exception:
        return {}
    finally:
        conn.close()

    lookup = {}
    for _, row in df.iterrows():
        oa = row["odds_api_name"]
        for col in ["short_name", "nickname", "abbreviation"]:
            val = row[col]
            if val and str(val).strip():
                lookup[str(val).strip().lower()] = oa
                lookup[strip_punct(str(val))] = oa
    return lookup


def load_canonical_games(sport):
    """Load canonical game list from JSON (created by run.py)."""
    path = ANSWER_KEYS_DIR / f".canonical_games_{sport}.json"
    if not path.exists():
        return []
    with open(path) as f:
        return json.load(f)


def _dict_lookup(name_raw, name_stripped, team_dict):
    """Try dict lookup with raw, stripped, and without trailing state suffix."""
    result = team_dict.get(name_raw) or team_dict.get(name_stripped)
    if result:
        return result
    # Try without trailing 2-letter state suffix (e.g., "st marys ca" -> "st marys")
    parts = name_stripped.rsplit(" ", 1)
    if len(parts) == 2 and len(parts[1]) == 2:
        result = team_dict.get(parts[0])
    return result


def resolve_team_names(scraped_away, scraped_home, team_dict, canonical_games):
    """Resolve scraped team names to Odds API format.

    Layer 1: Dictionary lookup with punctuation stripping.
    Layer 2: Game-level matching (both teams must match same canonical game).
    Returns (resolved_away, resolved_home).
    """
    # Layer 1: dict lookup (try raw lowercase, then stripped, then without state suffix)
    sa_raw = scraped_away.strip().lower()
    sh_raw = scraped_home.strip().lower()
    sa_stripped = strip_punct(scraped_away)
    sh_stripped = strip_punct(scraped_home)

    away = _dict_lookup(sa_raw, sa_stripped, team_dict)
    home = _dict_lookup(sh_raw, sh_stripped, team_dict)

    if away and home:
        return (away, home)

    # Layer 2: game-level matching against canonical game list
    for cg in canonical_games:
        ca = cg["away_team"]
        ch = cg["home_team"]
        ca_s = strip_punct(ca)
        ch_s = strip_punct(ch)

        # Canonical name without mascot (last word)
        ca_no_mascot = ca_s.rsplit(" ", 1)[0] if " " in ca_s else ca_s
        ch_no_mascot = ch_s.rsplit(" ", 1)[0] if " " in ch_s else ch_s

        # For each side: use dict result if available, else try substring match
        # Check both directions: scraped in canonical, canonical in scraped
        away_ok = away or (sa_stripped in ca_s) or (ca_no_mascot == sa_stripped) or (ca_no_mascot in sa_stripped)
        home_ok = home or (sh_stripped in ch_s) or (ch_no_mascot == sh_stripped) or (ch_no_mascot in sh_stripped)

        if away_ok and home_ok:
            return (away or ca, home or ch)

    # Fallback: return dict results or originals
    if not away:
        print(f"  Warning: Could not resolve away team '{scraped_away}'")
    if not home:
        print(f"  Warning: Could not resolve home team '{scraped_home}'")
    return (away or scraped_away, home or scraped_home)
