"""One-shot go-live gate: per book, what fraction of the Kalshi spread×total
ladder (BOTH signed grids, all N) comes back as a full 4-cell grid?

Run AFTER the bot has done at least one both_teams=True SGP scrape today.
Reads the bot market DB read-only; prints a table; persists nothing.
"""
import duckdb
from collections import defaultdict
from pathlib import Path

MARKET_DB = Path(__file__).resolve().parents[1] / "kalshi_mlb_rfq_market.duckdb"
LABELS = ["Home Spread + Over", "Home Spread + Under",
          "Away Spread + Over", "Away Spread + Under"]

def main():
    con = duckdb.connect(str(MARKET_DB), read_only=True)
    try:
        df = con.execute(
            "SELECT game_id, bookmaker, spread_line, total_line, combo "
            "FROM mlb_sgp_odds"
        ).fetchdf()
    finally:
        con.close()
    if df.empty:
        print("mlb_sgp_odds empty — run a both_teams=True SGP scrape first.")
        return

    # A (game, book, spread_line, total_line) cell-group is 'full' iff all 4
    # combo labels are present.
    full = defaultdict(int); seen = defaultdict(int)
    neg_full = defaultdict(int); pos_full = defaultdict(int)
    neg_seen = defaultdict(int); pos_seen = defaultdict(int)
    grp = df.groupby(["game_id", "bookmaker", "spread_line", "total_line"])
    for (game, book, sl, tl), sub in grp:
        labels = set(sub["combo"].unique())
        is_full = all(l in labels for l in LABELS)
        seen[book] += 1; full[book] += int(is_full)
        if sl < 0:
            neg_seen[book] += 1; neg_full[book] += int(is_full)
        else:
            pos_seen[book] += 1; pos_full[book] += int(is_full)

    print(f"{'book':12} {'all grids':>14} {'home(-) grids':>16} {'away(+) grids':>16}")
    for book in sorted(seen):
        def pct(f, s): return f"{f}/{s} ({100*f/s:.0f}%)" if s else "0/0"
        print(f"{book:12} {pct(full[book],seen[book]):>14} "
              f"{pct(neg_full[book],neg_seen[book]):>16} "
              f"{pct(pos_full[book],pos_seen[book]):>16}")
    print("\nThe away(+) column is the new capability — if it is ~0% for a "
          "book, that book has no dog-side alt SGP coverage and its away-margin "
          "combos will drop (by design).")

if __name__ == "__main__":
    main()
