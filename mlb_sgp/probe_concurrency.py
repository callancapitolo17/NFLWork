"""Concurrency-ceiling probe for the MLB SGP scrapers.

Times real ``price_sgps`` runs on a small fixed target set at
escalating ``parallelism`` levels and reports rows/sec + miss-rate per
level. The knee (throughput stops scaling, or misses spike) is the
book's empirical ceiling; ship the last clean level.

MANUAL-ONLY tool — never invoked by the bots. Budgets per the speedup
spec (docs/superpowers/specs/2026-06-08-sgp-scraper-speedup-design.md):
DK <= 150 pricing calls, PX <= 40 RFQs, cooldown between levels, and a
post-run health check at parallelism 1. Run in the morning window
(~7-8am PT, no MLB games in progress).

Usage (from repo root):
  mlb_sgp/venv/bin/python -m mlb_sgp.probe_concurrency --book dk
  mlb_sgp/venv/bin/python -m mlb_sgp.probe_concurrency --book px \
      --db kalshi_mlb_rfq/kalshi_mlb_rfq_market.duckdb
"""
from __future__ import annotations
import argparse
import sys
import time
from pathlib import Path

_REPO_ROOT = Path(__file__).resolve().parents[1]
_MLB_SGP_DIR = Path(__file__).resolve().parent
for p in (str(_REPO_ROOT), str(_MLB_SGP_DIR)):
    if p not in sys.path:
        sys.path.insert(0, p)

from mlb_sgp._shared import TargetLine, load_target_lines  # noqa: E402

DEFAULT_DB = str(_REPO_ROOT / "kalshi_mlb_rfq" / "kalshi_mlb_rfq_market.duckdb")
COOLDOWN_SEC = 20

# Per-book ramp plans. calls-per-level = targets_per_level * 4 combos.
#   DK: 5 levels x 6 targets x 4 = 120 calls + 4 health = 124 <= 150
#   PX: 4 levels x 2 targets x 4 =  32 calls + 4 health =  36 <= 40
BOOK_PLANS = {
    "dk": {"levels": [2, 4, 8, 12, 16], "targets_per_level": 6, "budget": 150},
    "px": {"levels": [2, 3, 4, 6],      "targets_per_level": 2, "budget": 40},
}


def pick_probe_targets(targets: list[TargetLine], n_games: int) -> list[TargetLine]:
    """One target per game (first n_games games, input order), preferring
    the main-ish spread (smallest |spread|) so the book always offers it."""
    by_game: dict[str, list[TargetLine]] = {}
    order: list[str] = []
    for t in targets:
        if t.game_id not in by_game:
            order.append(t.game_id)
        by_game.setdefault(t.game_id, []).append(t)
    picked = []
    for gid in order[:n_games]:
        picked.append(min(by_game[gid], key=lambda t: abs(t.spread)))
    return picked


def summarize_level(level: int, n_targets: int, n_rows: int,
                    elapsed_sec: float) -> dict:
    expected = n_targets * 4
    return {
        "level": level,
        "n_targets": n_targets,
        "n_rows": n_rows,
        "expected_rows": expected,
        "elapsed_sec": round(elapsed_sec, 2),
        "rows_per_sec": round(n_rows / elapsed_sec, 2) if elapsed_sec > 0 else 0.0,
        "miss_pct": (expected - n_rows) / expected * 100 if expected else 0.0,
    }


def _book_module(book: str):
    if book == "dk":
        from mlb_sgp import draftkings
        return draftkings
    if book == "px":
        from mlb_sgp import prophetx
        return prophetx
    raise SystemExit(f"unknown book {book!r} (use dk|px)")


def _build_client(book: str):
    if book == "dk":
        from mlb_sgp.dk_client import DraftKingsClient
        return DraftKingsClient()
    from mlb_sgp.prophetx_client import ProphetXClient
    return ProphetXClient()


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--book", required=True, choices=("dk", "px"))
    ap.add_argument("--db", default=DEFAULT_DB,
                    help="DuckDB holding mlb_target_lines (read-only)")
    args = ap.parse_args()

    plan = BOOK_PLANS[args.book]
    all_targets = [t for t in load_target_lines(args.db) if t.period == "FG"]
    if not all_targets:
        print(f"no FG target lines in {args.db} — run while the bot has "
              f"populated mlb_target_lines (or pass --db).")
        return 2
    probe_targets = pick_probe_targets(all_targets, plan["targets_per_level"])
    print(f"book={args.book} probe_targets={len(probe_targets)} "
          f"games={[t.game_id[:8] for t in probe_targets]}")

    mod = _book_module(args.book)
    client = _build_client(args.book)  # ONE persistent session for the whole ramp

    results = []
    for level in plan["levels"]:
        t0 = time.monotonic()
        rows = mod.price_sgps(probe_targets, periods=("FG",),
                              client=client, parallelism=level)
        s = summarize_level(level, len(probe_targets), len(rows),
                            time.monotonic() - t0)
        results.append(s)
        print(f"  level={s['level']:<3} rows={s['n_rows']}/{s['expected_rows']} "
              f"elapsed={s['elapsed_sec']}s rows/sec={s['rows_per_sec']} "
              f"miss={s['miss_pct']:.1f}%")
        # Stop the ramp on degradation: any misses at this level when the
        # previous level was clean, or throughput regressing >20%.
        if len(results) >= 2:
            prev = results[-2]
            deg_miss = s["miss_pct"] > prev["miss_pct"]
            deg_tput = s["rows_per_sec"] < prev["rows_per_sec"] * 0.8
            if deg_miss or deg_tput:
                print(f"  DEGRADATION at level {level} "
                      f"(miss {prev['miss_pct']:.1f}->{s['miss_pct']:.1f}%, "
                      f"tput {prev['rows_per_sec']}->{s['rows_per_sec']}). "
                      f"Stopping ramp.")
                break
        time.sleep(COOLDOWN_SEC)

    # Post-probe health check: one target at parallelism 1 must price.
    hc = mod.price_sgps(probe_targets[:1], periods=("FG",),
                        client=client, parallelism=1)
    print(f"health check: {len(hc)}/4 rows "
          f"{'OK' if hc else '** FAILED — book may have tripped a defense **'}")

    clean = [r for r in results if r["miss_pct"] == results[0]["miss_pct"]]
    if clean:
        print(f"recommended parallelism: {clean[-1]['level']} "
              f"(last level matching baseline miss-rate)")
    return 0


if __name__ == "__main__":
    sys.exit(main())
