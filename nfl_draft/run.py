"""NFL Draft Portal orchestrator.

Modes:
  --mode scrape: pull all venues' odds -> draft_odds + invoke legacy
                 edge_detector and consensus.
  --mode trades: poll Kalshi trade tape with cursor + dedup -> kalshi_trades
                 (stubbed in Task 20; implemented in Task 21).

Per-book error isolation: a failure in one scraper does NOT stop the others.
Seed runs at the start of every invocation, so dict edits in config/*.py
take effect on the next cron tick without manual re-seeding.
"""

import argparse
import sys
import traceback
from pathlib import Path

# Make the repo root importable so `import kalshi_draft.edge_detector` works
# when run.py is invoked as a module (`python -m nfl_draft.run`) or directly.
sys.path.insert(0, str(Path(__file__).resolve().parent.parent))

from nfl_draft.lib import seed
from nfl_draft.lib.quarantine import write_or_quarantine

SCRAPERS = {
    "kalshi": "nfl_draft.scrapers.kalshi",
    "draftkings": "nfl_draft.scrapers.draftkings",
    "fanduel": "nfl_draft.scrapers.fanduel",
    "bookmaker": "nfl_draft.scrapers.bookmaker",
    "wagerzon": "nfl_draft.scrapers.wagerzon",
    "hoop88": "nfl_draft.scrapers.hoop88",
}


def run_scrape(book: str) -> None:
    """Run --mode scrape for one or all books.

    After scrapers finish, invoke legacy edge_detector + consensus so the
    existing Kalshi dashboard tabs (detected_edges, consensus_board) stay
    fresh on the same cron tick.
    """
    seed.run()  # always reseed in case config dicts edited
    targets = list(SCRAPERS.keys()) if book == "all" else [book]
    total_mapped = 0
    total_unmapped = 0
    for name in targets:
        try:
            module = __import__(SCRAPERS[name], fromlist=["fetch_draft_odds"])
            print(f"[scrape] {name}: fetching...")
            rows = module.fetch_draft_odds()
            mapped, unmapped = write_or_quarantine(rows)
            total_mapped += mapped
            total_unmapped += unmapped
            print(f"[scrape] {name}: mapped={mapped} unmapped={unmapped}")
        except NotImplementedError:
            print(f"[scrape] {name}: NOT IMPLEMENTED, skipping")
        except Exception as e:
            print(f"[scrape] {name}: ERROR - {e}")
            traceback.print_exc()
    print(f"[scrape] TOTAL: mapped={total_mapped} unmapped={total_unmapped}")

    # Trigger legacy edge detector (writes to detected_edges).
    # detect_all_edges() is the library entry point; __main__ also calls
    # init_schema() first, which we replicate for safety.
    try:
        from nfl_draft.lib import db as nfl_db
        import kalshi_draft.edge_detector as edge_detector
        nfl_db.init_schema()
        edge_detector.detect_all_edges()
        print("[scrape] edge_detector: done")
    except Exception as e:
        print(f"[scrape] edge_detector: ERROR - {e}")
        traceback.print_exc()

    # Trigger legacy consensus scraper (writes to consensus_board).
    # consensus.run() is the library entry point.
    try:
        import kalshi_draft.consensus as consensus
        consensus.run()
        print("[scrape] consensus: done")
    except Exception as e:
        print(f"[scrape] consensus: ERROR - {e}")
        traceback.print_exc()


def run_trades() -> None:
    """Run --mode trades. Fully implemented in Task 21."""
    seed.run()
    from nfl_draft.scrapers import kalshi
    print("[trades] polling...")
    try:
        trades = kalshi.fetch_trades()
        print(f"[trades] ingested {len(trades)} trades")
    except NotImplementedError:
        print("[trades] NOT IMPLEMENTED (Task 21)")


def main() -> None:
    p = argparse.ArgumentParser(description="NFL Draft Portal orchestrator")
    p.add_argument("--mode", choices=["scrape", "trades"], required=True)
    p.add_argument(
        "--book",
        default="all",
        help="Scraper name ('kalshi', 'draftkings', ...) or 'all' (default).",
    )
    args = p.parse_args()
    if args.mode == "scrape":
        run_scrape(args.book)
    else:
        run_trades()


if __name__ == "__main__":
    main()
