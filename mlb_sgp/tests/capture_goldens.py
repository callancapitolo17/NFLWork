"""Capture frozen golden fixtures from a LIVE slate (issue #42).

Run this by hand, never from a test. It drives each book's real
``price_sgps`` once against a small target grid, records every payload the
orchestrator fetched, and writes a gzipped
``tests/golden/<book>_golden.json.gz`` holding both the recorded payloads and
the resulting ``PricedRow`` set.

``test_golden_baselines.py`` then replays those payloads offline forever, so
the pinned rows are reproducible without a network call and without touching
any DuckDB.

Only recapture when a book's call pattern changes on purpose — and only from
a slate where that book is healthy. Recapturing mid-repair bakes a broken
book's output in as "correct".

    python3 mlb_sgp/tests/capture_goldens.py                 # all books
    python3 mlb_sgp/tests/capture_goldens.py --books novig

Side effects: writes ONLY to mlb_sgp/tests/golden/. No DB is opened.
"""
from __future__ import annotations

import argparse
import sys
from datetime import datetime, timezone
from pathlib import Path

_REPO_ROOT = Path(__file__).resolve().parents[2]
for _p in (str(_REPO_ROOT), str(_REPO_ROOT / "mlb_sgp"),
           str(_REPO_ROOT / "Answer Keys")):
    if _p not in sys.path:
        sys.path.insert(0, _p)

from mlb_sgp._shared import TargetLine                       # noqa: E402
from mlb_sgp.tests import golden_replay as gr                # noqa: E402
from mlb_sgp.verify_books import ALL_BOOKS, pick_game        # noqa: E402

# One game, a couple of rungs: enough to exercise the grid, the ML x total
# phase and the sanity filter without a wire-call burst that would trip
# Novig's rate limit mid-capture.
CAPTURE_SPREADS = (-1.5,)
CAPTURE_TOTALS = (8.5, 9.0)


class _Patcher:
    """monkeypatch-alike for use outside pytest."""

    def __init__(self):
        self._undo = []

    def setattr(self, target, name, value):
        self._undo.append((target, name, getattr(target, name)))
        setattr(target, name, value)

    def undo(self):
        for target, name, original in reversed(self._undo):
            setattr(target, name, original)
        self._undo.clear()


def build_client(book: str):
    from kalshi_common.sgp_service import SGPService
    service = SGPService(books=(book,), health_db_path=None)
    return service._ensure_client(book).client


def capture_book(book: str, targets) -> dict | None:
    sink: dict = {}
    patcher = _Patcher()
    try:
        client = build_client(book)
        fetchers = gr.record(book, patcher, client, sink)
        rows = gr.price_with(book, targets, client, fetchers)
    except Exception as exc:
        print(f"  {book}: capture FAILED — {exc!r}")
        return None
    finally:
        patcher.undo()

    if not rows:
        print(f"  {book}: produced 0 rows — NOT written "
              f"(a golden pinning an empty result is worse than none)")
        return None
    print(f"  {book}: {len(rows)} rows, "
          f"{sum(len(v) for v in sink.values())} recorded payloads")
    return {"book": book,
            "captured_at": datetime.now(timezone.utc).isoformat(),
            "targets": [{"game_id": t.game_id, "home_team": t.home_team,
                         "away_team": t.away_team,
                         "commence_time": t.commence_time.isoformat(),
                         "period": t.period, "spread": t.spread,
                         "total": t.total} for t in targets],
            "payloads": sink,
            "rows": gr.sort_rows(gr.rows_to_json(rows))}


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--books", default=",".join(ALL_BOOKS))
    args = parser.parse_args()

    game = pick_game()
    print(f"game: {game.away_team} @ {game.home_team}  {game.commence_time}")
    targets = [TargetLine(game_id=game.game_id, home_team=game.home_team,
                          away_team=game.away_team,
                          commence_time=game.commence_time,
                          period="FG", spread=s, total=t)
               for s in CAPTURE_SPREADS for t in CAPTURE_TOTALS]

    gr.GOLDEN_DIR.mkdir(parents=True, exist_ok=True)
    written = 0
    for book in (b.strip() for b in args.books.split(",") if b.strip()):
        blob = capture_book(book, targets)
        if blob is None:
            continue
        gr.save_golden(book, blob)
        written += 1
    print(f"\nwrote {written} golden fixture(s) to {gr.GOLDEN_DIR}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
