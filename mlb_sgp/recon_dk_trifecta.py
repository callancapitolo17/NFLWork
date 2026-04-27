"""
DraftKings Trifecta Recon Script

Captures DK's `parlays/v1/sgp/events/{event_id}` response for one or two
known MLB games, plus an attempted `calculateBets` POST with candidate
selection IDs. Output guides the Plan #2 leg-resolver implementation:
- "First Team to Score" market name + selection ID format
- "1st 5 Innings Winner" market — 2-way ML or 3-way?
- "Team Total Over/Under" market for the GRAND-SLAM 4th leg
- Whether DK accepts these legs as a single SGP

Usage:
    1. Find one or two DK MLB game URLs (e.g. https://sportsbook.draftkings.com/event/<event_id>)
    2. Update DK_EVENT_IDS below
    3. Run: python3 mlb_sgp/recon_dk_trifecta.py
    4. Inspect the resulting recon_dk_trifecta_<event_id>.json files
    5. Plan #2 reads these to populate the leg resolvers

Manual one-shot artifact. Not part of any pipeline.
"""
from __future__ import annotations

import json
import logging
import sys
from pathlib import Path
from typing import Optional

try:
    from curl_cffi import requests
except ImportError:
    print("ERROR: curl_cffi not installed. The DK SGP scraper venv has it.", file=sys.stderr)
    print("       cd mlb_sgp && python -m venv venv && venv/bin/pip install curl_cffi", file=sys.stderr)
    sys.exit(1)

# ===== CONFIG =====
# DK event IDs to characterize. If empty, the script auto-discovers up to
# MAX_AUTO_EVENTS upcoming MLB games via the same eventgroup API the
# production DK SGP scraper uses (mlb_sgp/scraper_draftkings_sgp.py).
# Override this list to pin specific games, e.g. for re-running on a known fixture.
DK_EVENT_IDS: list[str] = [
    # "32109876",  # populate to target specific games; otherwise leave empty
]
MAX_AUTO_EVENTS = 2

# DK API endpoints (public, no auth required). The full SGP path lives under
# /sites/US-SB/api/sportscontent/ — same prefix the production scraper uses.
SGP_EVENTS_URL_TPL = (
    "https://sportsbook-nash.draftkings.com/sites/US-SB/api/sportscontent/"
    "parlays/v1/sgp/events/{event_id}"
)

OUT_DIR = Path(__file__).parent
logging.basicConfig(level=logging.INFO, format="%(asctime)s [%(levelname)s] %(message)s")
log = logging.getLogger("recon_dk_trifecta")


def init_session():
    """Reuse the production DK SGP scraper's session init: curl_cffi + Chrome
    impersonation + a warmup GET to seed cookies. Importing here so that
    smoke-tests of this module don't require curl_cffi at import time."""
    sys.path.insert(0, str(Path(__file__).parent))
    from scraper_draftkings_sgp import init_session as _init  # type: ignore
    return _init()


def discover_dk_events(session, max_events: int) -> list[dict]:
    """Discover today's MLB events via the same eventgroup API the production
    scraper uses. Returns a list of {dk_event_id, name, start_time} dicts.
    Limited to max_events to keep recon fast."""
    sys.path.insert(0, str(Path(__file__).parent))
    from scraper_draftkings_sgp import fetch_dk_events  # type: ignore
    events = fetch_dk_events(session)
    log.info("DK eventgroup returned %d MLB events", len(events))
    return events[:max_events]


def fetch_event_selections(session, event_id: str) -> Optional[dict]:
    """Fetch the full SGP event payload (~2MB JSON of all selection IDs)."""
    url = SGP_EVENTS_URL_TPL.format(event_id=event_id)
    log.info("GET %s", url)
    try:
        r = session.get(url, timeout=30)
        r.raise_for_status()
        return r.json()
    except Exception as e:
        log.exception("Failed to fetch event %s: %s", event_id, e)
        return None


def summarize_markets(payload: dict) -> dict:
    """Walk the payload and pull out market-name / selection-id summaries
    relevant to triple-play / grand-slam legs:
      - "First Team to Score" or similar
      - "1st 5 Innings Winner"
      - "Game Moneyline" (h2h)
      - Team Total Over/Under
    Returns a dict suitable for Plan #2's resolver design."""
    summary = {
        "first_to_score": [],
        "first_5_innings_winner": [],
        "moneyline": [],
        "team_totals": [],
        "other_markets_found": [],
    }
    # The DK SGP payload structure varies; this is a permissive walker
    # that records anything that looks like a market.
    def walk(node, path=""):
        if isinstance(node, dict):
            name = (node.get("marketName") or node.get("name") or "").strip()
            if name:
                lname = name.lower()
                bucket = None
                if "first" in lname and "score" in lname:
                    bucket = "first_to_score"
                elif "1st 5" in lname or "first 5" in lname:
                    bucket = "first_5_innings_winner"
                elif "moneyline" in lname or lname in ("h2h", "game"):
                    bucket = "moneyline"
                elif "team total" in lname:
                    bucket = "team_totals"
                else:
                    summary["other_markets_found"].append({"path": path, "name": name})
                if bucket:
                    selections = node.get("selections") or node.get("outcomes") or []
                    summary[bucket].append({
                        "path": path,
                        "name": name,
                        "selections": [
                            {"name": s.get("label") or s.get("name"),
                             "id": s.get("id") or s.get("selectionId")}
                            for s in selections
                            if isinstance(s, dict)
                        ],
                    })
            for k, v in node.items():
                walk(v, f"{path}.{k}" if path else k)
        elif isinstance(node, list):
            for i, v in enumerate(node):
                walk(v, f"{path}[{i}]")
    walk(payload)
    return summary


def main() -> int:
    session = init_session()

    # Resolve which event IDs to target. Manual override wins; otherwise auto.
    if DK_EVENT_IDS:
        targets = [{"dk_event_id": eid, "name": "(manual override)", "start_time": ""}
                   for eid in DK_EVENT_IDS]
    else:
        targets = discover_dk_events(session, MAX_AUTO_EVENTS)
        if not targets:
            log.error("No upcoming MLB events found. Specials may not be posted, "
                      "or DK eventgroup API changed.")
            return 1

    for evt in targets:
        event_id = str(evt["dk_event_id"])
        log.info("=== Reconning event %s — %s @ %s ===",
                 event_id, evt.get("name", ""), evt.get("start_time", ""))
        payload = fetch_event_selections(session, event_id)
        if payload is None:
            log.warning("Skipping %s (fetch failed)", event_id)
            continue
        out_full = OUT_DIR / f"recon_dk_trifecta_{event_id}.json"
        out_full.write_text(json.dumps(payload, indent=2, default=str))
        log.info("Saved full payload to %s (%d bytes)", out_full.name, out_full.stat().st_size)

        summary = summarize_markets(payload)
        out_summary = OUT_DIR / f"recon_dk_trifecta_{event_id}_summary.json"
        out_summary.write_text(json.dumps(summary, indent=2, default=str))
        n_first = len(summary["first_to_score"])
        n_5inn  = len(summary["first_5_innings_winner"])
        n_ml    = len(summary["moneyline"])
        n_tt    = len(summary["team_totals"])
        log.info("Summary: first_to_score=%d  first_5_innings=%d  moneyline=%d  team_totals=%d",
                 n_first, n_5inn, n_ml, n_tt)
        log.info("Wrote market summary to %s", out_summary.name)
    log.info("Done.")
    return 0


if __name__ == "__main__":
    sys.exit(main())
