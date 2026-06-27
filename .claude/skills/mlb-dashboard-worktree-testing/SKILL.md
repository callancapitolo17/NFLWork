---
name: mlb-dashboard-worktree-testing
description: >-
  How to test the MLB dashboard / pricing pipeline from inside a git worktree
  without touching main's live data. Load this when testing, rendering, or
  re-pricing the MLB dashboard, running the MLB pipeline from a worktree, or
  using seed_test_data.sh / test_dashboard.sh.
---

# Testing the MLB Dashboard / Pipeline from a Worktree

The worktree is self-contained — every MLB DB path derives from the running
code's location (`Tools.R::nflwork_root()` / `derive_repo_root()`), so it
reads/writes its own seeded DBs and never touches main's live data. On `main`,
paths fall back to `~/NFLWork` byte-identically.

## Steps

1. **Seed the worktree with live data (CoW clone — instant on APFS):**
   ```
   "Answer Keys/MLB Dashboard/seed_test_data.sh"
   ```
   This copy-on-write clones the live DuckDBs into the worktree. Seeded
   `*.duckdb` files are gitignored and removed with the worktree.

2. **Render + serve the dashboard on port :8093:**
   ```
   "Answer Keys/MLB Dashboard/test_dashboard.sh"
   ```
   Flags:
   - `--mlb` — re-price (re-run the model) before rendering
   - `--pipeline` — full re-scrape + re-price

## Why this is safe

Because the running code derives its repo root from its own location, a worktree
checkout reads and writes only the DBs seeded inside that worktree. Main's live
DBs are never read or mutated. This is why testing from a worktree is preferred
over testing against `~/NFLWork` directly.

## Reminder

- **Never symlink the `.duckdb` files** into the worktree — DuckDB writes WAL
  files next to the database path, so removing the worktree would lose
  uncommitted WAL data. `seed_test_data.sh` *copies* (CoW), it does not symlink.
- For final verification, prefer testing from `main` after merge.
