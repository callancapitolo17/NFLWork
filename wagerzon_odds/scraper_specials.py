"""
Wagerzon MLB Specials Scraper

Fetches all MLB specials (lg=4899) from the NewScheduleHelper JSON endpoint
and writes one row per prop to wagerzon.duckdb / wagerzon_specials table.

The pricer (Answer Keys/mlb_triple_play.R) reads from this table.

Usage:
    python3 scraper_specials.py            # one-shot scrape, writes snapshot
    python3 scraper_specials.py --dry-run  # fetch + parse, don't write
"""
from __future__ import annotations

import argparse
import json
import logging
import re
import sys
from datetime import datetime, timezone
from pathlib import Path
from typing import Optional

import duckdb
import requests
from dotenv import load_dotenv

# Load .env from the canonical location (bet_logger/.env in the repo root).
# scraper_v2.login() reads WAGERZON_USERNAME/PASSWORD from env; we ensure they
# are loaded before importing scraper_v2 so worktree runs find the creds too.
_repo_root = Path(__file__).resolve().parent.parent
# Worktrees live at <repo>/.worktrees/<branch>/; real repo lives at <repo>/
# Walk up until we find bet_logger/.env or exhaust parents.
for _candidate in [_repo_root, _repo_root.parent, _repo_root.parent.parent]:
    _env = _candidate / "bet_logger" / ".env"
    if _env.exists():
        load_dotenv(_env)
        break

# Reuse the existing login + base-URL logic
sys.path.insert(0, str(Path(__file__).parent))
from scraper_v2 import login  # noqa: E402
from config import WAGERZON_BASE_URL  # noqa: E402

DB_PATH = Path(__file__).parent / "wagerzon.duckdb"
SPECIALS_URL_TPL = WAGERZON_BASE_URL + "/wager/NewSchedule" + "Helper.aspx?WT=0&lg={lg}"

# Regex used to extract the subject team. Anchors on known prop-type tokens
# so multi-word team names (WHITE SOX, RED SOX, BLUE JAYS) work. Extend the
# alternation when new prop types are added (e.g. DOUBLE-PLAY, MEGA).
PROP_TYPE_RE = re.compile(r"^(.+?)\s+(TRIPLE-PLAY|GRAND-SLAM)\b")
SECTION_RE = re.compile(r"\b(TRIPLE-PLAY|GRAND-SLAM)\b")

logging.basicConfig(level=logging.INFO, format="%(asctime)s [%(levelname)s] %(message)s")
log = logging.getLogger("scraper_specials")


def extract_team_from_htm(htm: str) -> Optional[str]:
    """Return the subject team (e.g. 'DODGERS', 'WHITE SOX') from a prop
    description, or None if the prop is cross-game / multi-team / unsupported.
    """
    if not htm:
        return None
    m = PROP_TYPE_RE.match(htm)
    if not m:
        return None
    team = m.group(1).strip()
    # Cross-game props ("CHC-SDG G-S", "MIA-LAD TRIPLE-PLAY") have a hyphenated
    # abbreviation in the team slot — those are game-level, not single-team.
    if "-" in team:
        return None
    return team


def extract_prop_type(htm: str) -> Optional[str]:
    m = SECTION_RE.search(htm or "")
    return m.group(1) if m else None


def parse_specials_json(payload: dict, sport: str, league_id: int) -> list[dict]:
    """Parse the NewScheduleHelper response into row dicts.

    Returns a list of dicts ready for DuckDB insert. Drops rows that aren't
    priceable single-team triple-plays or grand-slams (cross-game props,
    multi-team parlays, HR/hits combos, etc.). Filtering happens here so the
    DB stays focused; raw cross-game props can be added later by relaxing the
    filter without changing the schema.
    """
    scraped_at = datetime.now(timezone.utc).replace(microsecond=0)
    rows: list[dict] = []

    leagues_outer = payload.get("result", {}).get("listLeagues", [])
    if not leagues_outer or not leagues_outer[0]:
        return rows
    for league in leagues_outer[0]:
        # We only want MLB - SPECIALS leagues; other section headers (e.g.
        # MLB - REGULAR) are scraped by other tools.
        desc = (league.get("Description") or "").upper()
        if "SPECIALS" not in desc:
            continue
        for g in league.get("Games", []):
            htm = g.get("htm") or ""
            prop_type = extract_prop_type(htm)
            if prop_type not in ("TRIPLE-PLAY", "GRAND-SLAM"):
                continue
            team = extract_team_from_htm(htm)
            if team is None:
                continue  # cross-game props skipped
            game_lines = g.get("GameLines") or [{}]
            odds_str = (game_lines[0] or {}).get("odds") or ""
            oddsh    = (game_lines[0] or {}).get("oddsh") or ""
            try:
                # oddsh has explicit sign, e.g. "+190" or "-110"; odds is unsigned
                odds = int(oddsh) if oddsh else int(odds_str)
            except (TypeError, ValueError):
                log.warning("Could not parse odds for hnum=%s: oddsh=%r odds=%r",
                            g.get("hnum"), oddsh, odds_str)
                continue
            gmdt = g.get("gmdt") or ""
            gmtm = g.get("gmtm") or ""
            game_date = datetime.strptime(gmdt, "%Y%m%d").date() if gmdt else None
            try:
                # Wagerzon's gmtm appears to be local Eastern time — the scraper
                # stores it verbatim (no tz conversion). Pricer's 12-hour filter
                # uses commence_time from mlb_consensus_temp instead, which is
                # authoritative UTC.
                game_time = datetime.strptime(f"{gmdt} {gmtm}", "%Y%m%d %H:%M:%S")
            except ValueError:
                game_time = None
            rows.append({
                "scraped_at":      scraped_at,
                "sport":           sport,
                "league_id":       league_id,
                "game_id":         g.get("idgm"),
                "rotation_number": g.get("hnum"),
                "prop_type":       prop_type,
                "description":     htm,
                "team":            team,
                "line":            None,  # not applicable for SCR/period legs
                "odds":            odds,
                "game_date":       game_date,
                "game_time":       game_time,
            })
    return rows


def fetch_specials_json(lg: int) -> dict:
    s = requests.Session()
    login(s)
    url = SPECIALS_URL_TPL.format(lg=lg)
    # XHR headers required: without them the endpoint returns the HTML login page
    r = s.get(url, timeout=30, headers={
        "Accept": "application/json, text/plain, */*",
        "X-Requested-With": "XMLHttpRequest",
    })
    r.raise_for_status()
    return r.json()


def write_rows(con: duckdb.DuckDBPyConnection, rows: list[dict]) -> None:
    """Idempotent: ensures table exists, then bulk inserts."""
    con.execute("""
        CREATE TABLE IF NOT EXISTS wagerzon_specials (
            scraped_at      TIMESTAMP,
            sport           VARCHAR,
            league_id       INTEGER,
            game_id         INTEGER,
            rotation_number INTEGER,
            prop_type       VARCHAR,
            description     VARCHAR,
            team            VARCHAR,
            line            DOUBLE,
            odds            INTEGER,
            game_date       DATE,
            game_time       TIMESTAMP
        )
    """)
    if not rows:
        log.info("No rows to write.")
        return
    con.executemany(
        "INSERT INTO wagerzon_specials VALUES "
        "(?,?,?,?,?,?,?,?,?,?,?,?)",
        [
            (r["scraped_at"], r["sport"], r["league_id"], r["game_id"],
             r["rotation_number"], r["prop_type"], r["description"],
             r["team"], r["line"], r["odds"], r["game_date"], r["game_time"])
            for r in rows
        ],
    )


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--lg", type=int, default=4899,
                    help="Wagerzon league id (4899=MLB-SPECIALS)")
    ap.add_argument("--sport", default="mlb")
    # Positional `sport` is what `Answer Keys/run.py` passes to every scraper
    # (e.g. `python scraper_specials.py mlb`). When present, it overrides
    # `--sport`. Manual callers can still use `--sport <name>` or omit both.
    ap.add_argument("sport_positional", nargs="?", default=None,
                    help="Sport name (positional). Overrides --sport when provided.")
    ap.add_argument("--dry-run", action="store_true",
                    help="Fetch + parse, don't write to DB")
    ap.add_argument("--cron", action="store_true",
                    help="Quiet mode: log only WARN+, suitable for scheduled runs")
    args = ap.parse_args()

    sport = args.sport_positional or args.sport

    if args.cron:
        logging.getLogger().setLevel(logging.WARNING)

    try:
        log.info("Fetching specials lg=%d", args.lg)
        payload = fetch_specials_json(args.lg)
        rows = parse_specials_json(payload, sport=sport, league_id=args.lg)
        log.info("Parsed %d priceable rows", len(rows))
    except Exception:
        log.exception("Scrape failed")
        return 2

    if args.dry_run:
        for r in rows[:5]:
            print(r)
        log.info("Dry-run: not writing.")
        return 0

    if not rows:
        # Empty result is NOT a failure — it's a successful determination that
        # Wagerzon has nothing posted. Returning 0 lets the orchestrator (run.py)
        # treat empty days as a clean run rather than a scraper error.
        log.warning("No priceable specials found — not writing snapshot")
        return 0

    con = duckdb.connect(str(DB_PATH))
    try:
        # `with con:` makes the snapshot atomic — commits on success,
        # rolls back on exception so a partial-failure never leaves the
        # table with a half-written snapshot for the latest scraped_at.
        with con:
            write_rows(con, rows)
        log.info("Wrote %d rows to %s/wagerzon_specials", len(rows), DB_PATH.name)
    finally:
        con.close()
    return 0


if __name__ == "__main__":
    sys.exit(main())
