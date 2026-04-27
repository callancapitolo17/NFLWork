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

# ===== CONFIG — UPDATE THESE BEFORE RUNNING =====
# DK event IDs for ~2 MLB games we want to characterize. Find these by:
# 1. Open https://sportsbook.draftkings.com/leagues/baseball/mlb in a browser
# 2. Click into a game, look at the URL: /event/<event_id>
# 3. Pick games that have triple-play / grand-slam comparable markets
#    (any MLB game should work — we only need the selection IDs)
DK_EVENT_IDS: list[int] = [
    # 32109876,  # example — replace with real IDs before running
    # 32109877,
]

# DK API endpoints (public, no auth required)
SGP_EVENTS_URL_TPL = "https://sportsbook-nash.draftkings.com/parlays/v1/sgp/events/{event_id}"

OUT_DIR = Path(__file__).parent
logging.basicConfig(level=logging.INFO, format="%(asctime)s [%(levelname)s] %(message)s")
log = logging.getLogger("recon_dk_trifecta")


def fetch_event_selections(event_id: int) -> Optional[dict]:
    """Fetch the full SGP event payload (~2MB JSON of all selection IDs)."""
    url = SGP_EVENTS_URL_TPL.format(event_id=event_id)
    log.info("GET %s", url)
    try:
        r = requests.get(url, impersonate="chrome", timeout=30)
        r.raise_for_status()
        return r.json()
    except Exception as e:
        log.exception("Failed to fetch event %d: %s", event_id, e)
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
    if not DK_EVENT_IDS:
        log.error("DK_EVENT_IDS is empty. Edit %s and add 1-2 event ids first.", __file__)
        return 1
    for event_id in DK_EVENT_IDS:
        log.info("=== Reconning event %d ===", event_id)
        payload = fetch_event_selections(event_id)
        if payload is None:
            log.warning("Skipping %d (fetch failed)", event_id)
            continue
        out_full = OUT_DIR / f"recon_dk_trifecta_{event_id}.json"
        out_full.write_text(json.dumps(payload, indent=2, default=str))
        log.info("Saved full payload to %s", out_full.name)

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
