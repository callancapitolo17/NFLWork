# FG `-1.5` Favorite Correlated Parlay Leak — Investigation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Identify the specific mechanism causing FG `-1.5` fav correlated parlays to bleed `-8.3%` ROI on 204 bets (ex-bad-week, `-11.5%` on the full 222-bet sample), with statistical rigor sufficient to recommend a concrete sizing or selection rule.

**Architecture:** Hypothesis-driven analytical investigation. Build one Python analysis script + one R helper to extract sample-level correlations, test six concrete mechanisms against the bet history, rank hypotheses by evidence strength, and write a findings memo. No production code changes — the deliverable is a memo and an updated memory note.

**Tech Stack:** Python 3 (pandas, scipy, numpy, openpyxl, duckdb), R (for `mlb_game_samples` correlation extraction), DuckDB (data join layer).

---

## File Structure

- Create: `docs/superpowers/analysis/2026-05-11-fg-fav-leak/`
  - `01_extract_bets.py` — load 627 MLB parlays from xlsx + join to placed_parlays edge_pct
  - `02_baseline.py` — bootstrap CI on FG-fav ROI; binomial tests vs model
  - `03_hypothesis_tests.py` — runs H1-H6 against the dataset, emits a `findings_data.json`
  - `04_sample_correlation.R` — computes empirical fav+over correlation from `mlb_game_samples`, compared to model's predicted correlation_factor
  - `05_write_memo.py` — reads `findings_data.json` + writes `findings.md`
  - `findings.md` — output memo (generated)
- Read-only: `Answer Keys/mlb_correlated_parlay.R`, `Answer Keys/Tools.R`, `MLB Dashboard/mlb_dashboard.duckdb`, `Answer Keys/mlb_mm.duckdb`, `/tmp/bet_tracking.xlsx` (or re-fetched)

---

## Worktree / Version Control

- Branch: `analysis/fg-fav-leak-2026-05-11` (already created)
- Worktree: `/Users/callancapitolo/NFLWork/.worktrees/fg-fav-leak` (already created)
- Commits: one per Task. Commit messages prefixed `analysis(fg-fav-leak): <task>`.
- This is **analysis-only** — no production code changes. Merge strategy: keep the branch as a frozen record, or merge to `main` to preserve the memo under `docs/superpowers/analysis/`. User picks at the end.
- Cleanup if branch is merged: `git worktree remove .worktrees/fg-fav-leak && git branch -d analysis/fg-fav-leak-2026-05-11`
- Cleanup if branch is kept: leave worktree in place, commit any follow-up directly there.

---

## Documentation

- The findings memo (`docs/superpowers/analysis/2026-05-11-fg-fav-leak/findings.md`) IS the documentation deliverable.
- If the investigation produces an actionable rule (e.g., "skip fav+over when total ≥ 9.5"), add it to `Answer Keys/CLAUDE.md` under a new "Known model biases" subsection.
- Update memory file `/Users/callancapitolo/.claude-personal/projects/-Users-callancapitolo-NFLWork/memory/mlb_parlay_edge_overestimation.md` with the final verdict.
- If any hypothesis tests conclusively, also add a one-line entry under "Reference" in `MEMORY.md`.

---

## Pre-merge Review Checklist

- All six hypothesis tests ran without errors and produced numeric output (not NaN-everywhere).
- Bootstrap CIs reported with seed for reproducibility.
- Findings memo cites every claim with the underlying number + sample size.
- No production `.R` / `.py` files modified — `git diff main..HEAD --stat` should show only `docs/` and `.gitignore` additions.

---

## Pre-Flight: Data Sources Sanity Check

Before any task, verify these files/tables exist and have the expected row counts. If any of these fails, STOP and report — the data foundation is broken.

- [ ] **Step 0.1:** Confirm `/tmp/bet_tracking.xlsx` exists and has `MLB Summary` tab with 625 bets. If missing, re-fetch via Google Drive MCP `download_file_content` with `exportMimeType: application/vnd.openxmlformats-officedocument.spreadsheetml.sheet`.

Run: `ls -la /tmp/bet_tracking.xlsx && /Users/callancapitolo/NFLWork/bet_logger/venv/bin/python3 -c "import openpyxl; w=openpyxl.load_workbook('/tmp/bet_tracking.xlsx', data_only=True); print(w.sheetnames, w['MLB Summary']['B5'].value)"`
Expected: file present, prints `['Sheet1', 'Sheet8', 'Shared', 'Stats', 'Sheet6', 'CBB Summary', '_cbb_ak', 'MLB Summary'] 625`

- [ ] **Step 0.2:** Confirm dashboard DB has 126 `placed` parlays with `edge_pct`.

Run: `/Users/callancapitolo/NFLWork/bet_logger/venv/bin/python3 -c "import duckdb; c=duckdb.connect('/Users/callancapitolo/NFLWork/Answer Keys/MLB Dashboard/mlb_dashboard.duckdb', read_only=True); print(c.execute(\"SELECT COUNT(*) FROM placed_parlays WHERE status='placed' AND edge_pct IS NOT NULL\").fetchone())"`
Expected: prints `(126,)` (or close to it — exact count may drift if dashboard logs new bets between sessions)

- [ ] **Step 0.3:** Confirm `mlb_mm.duckdb` has `mlb_game_samples` with at least 3000 rows.

Run: `/Users/callancapitolo/NFLWork/bet_logger/venv/bin/python3 -c "import duckdb; c=duckdb.connect('/Users/callancapitolo/NFLWork/Answer Keys/mlb_mm.duckdb', read_only=True); print(c.execute('SELECT COUNT(DISTINCT game_id) FROM mlb_game_samples').fetchone())"`
Expected: prints a number ≥ 20 distinct game_ids (one set of resamples per game in the latest slate)

---

## Tasks

### Task 1: Project setup + extract canonical FG-fav-`1.5` dataset

**Files:**
- Create: `docs/superpowers/analysis/2026-05-11-fg-fav-leak/01_extract_bets.py`
- Create: `docs/superpowers/analysis/2026-05-11-fg-fav-leak/.gitignore`
- Output: `/tmp/fg_fav_leak/bets.parquet` (gitignored — data not source)

**What this does:** Pulls all 627 MLB correlated parlays from xlsx, joins to dashboard `placed_parlays.edge_pct` where available, parses slice/line_type/total/ou, saves a clean dataset for downstream analysis.

- [ ] **Step 1.1: Create directory + .gitignore**

```bash
cd /Users/callancapitolo/NFLWork/.worktrees/fg-fav-leak
mkdir -p docs/superpowers/analysis/2026-05-11-fg-fav-leak
mkdir -p /tmp/fg_fav_leak
cat > docs/superpowers/analysis/2026-05-11-fg-fav-leak/.gitignore <<'EOF'
# Intermediate outputs not committed (memo + scripts only)
*.parquet
*.csv
findings_data.json
EOF
```

- [ ] **Step 1.2: Write `01_extract_bets.py`**

```python
#!/usr/bin/env python3
"""Extract 627 MLB correlated parlays + join to placed_parlays edge_pct."""
import openpyxl, duckdb, re, sys
import pandas as pd
import numpy as np
from pathlib import Path

XLSX = Path('/tmp/bet_tracking.xlsx')
DASH_DB = Path('/Users/callancapitolo/NFLWork/Answer Keys/MLB Dashboard/mlb_dashboard.duckdb')
OUT = Path('/tmp/fg_fav_leak/bets.parquet')

# --- Load xlsx ---
wb = openpyxl.load_workbook(XLSX, data_only=True)
ws = wb['Sheet1']
rows = []
for r in range(2, ws.max_row + 1):
    sport = ws.cell(row=r, column=3).value
    bet_type = ws.cell(row=r, column=5).value
    desc = ws.cell(row=r, column=4).value
    if sport != 'MLB' or bet_type != 'Parlay' or not desc or 'PARLAY (2 TEAMS)' not in str(desc):
        continue
    rows.append({
        'date': ws.cell(row=r, column=1).value,
        'platform': ws.cell(row=r, column=2).value,
        'desc': str(desc),
        'odds_amer': ws.cell(row=r, column=7).value,
        'stake': ws.cell(row=r, column=8).value,
        'dec_odds': ws.cell(row=r, column=9).value,
        'result': ws.cell(row=r, column=10).value,
        'pnl': ws.cell(row=r, column=13).value,
    })
df = pd.DataFrame(rows)
df['date'] = pd.to_datetime(df['date'])
df['stake'] = pd.to_numeric(df['stake'], errors='coerce').fillna(0)
df['pnl'] = pd.to_numeric(df['pnl'], errors='coerce').fillna(0)
df['dec_odds'] = pd.to_numeric(df['dec_odds'], errors='coerce')

# --- Parse description ---
df['is_f5'] = df['desc'].apply(lambda d: bool(re.search(r'\b(1H|F5|1ST 5|FIRST 5)\b', d, re.IGNORECASE)))
df['slice'] = df['is_f5'].map({True: 'F5', False: 'FG'})

def line_type(d):
    n = d.replace('½', '.5').replace('1/2', '.5')
    if re.search(r'\-1\.5', n): return 'fav-1.5'
    if re.search(r'\+1\.5', n): return 'dog+1.5'
    if re.search(r'[+\-]\.?5', n): return 'half(F5)'
    return 'other'
df['line_type'] = df['desc'].apply(line_type)
df['ou'] = df['desc'].apply(lambda d: 'over' if re.search(r'TOTAL O', d, re.IGNORECASE) else ('under' if re.search(r'TOTAL U', d, re.IGNORECASE) else 'unk'))

def total_num(d):
    m = re.search(r'TOTAL [ou](\d+)([½.]?\d*)', d, re.IGNORECASE)
    if m:
        return int(m.group(1)) + (0.5 if ('½' in m.group(2) or '.5' in m.group(2)) else 0)
    return None
df['total_line'] = df['desc'].apply(total_num)

# Extract spread leg American odds. Pattern: `TEAM ±1½±NNN` or similar
def parse_spread_amer(d):
    n = d.replace('½', '.5')
    m = re.search(r'[+\-]1\.5\s*([+\-]\d+|EV)', n)
    if m:
        v = m.group(1).replace('EV', '+100')
        try: return int(v)
        except: return None
    return None
df['spread_amer'] = df['desc'].apply(parse_spread_amer)

# Extract total leg American odds. Pattern: `TOTAL [ou]N±NNN`
def parse_total_amer(d):
    n = d.replace('½', '.5')
    m = re.search(r'TOTAL [ou]\d+\.?5?\s*([+\-]\d+|EV)', n, re.IGNORECASE)
    if m:
        v = m.group(1).replace('EV', '+100')
        try: return int(v)
        except: return None
    return None
df['total_amer'] = df['desc'].apply(parse_total_amer)

# --- Teams ---
def teams(d):
    m = re.search(r'\(([^()]+?)\s+vrs\s+([^()]+?)\)', d, re.IGNORECASE)
    return (m.group(1).strip().upper(), m.group(2).strip().upper()) if m else ('', '')
df[['away_full', 'home_full']] = df['desc'].apply(lambda d: pd.Series(teams(d)))
df['week'] = df['date'].dt.to_period('W-SUN').apply(lambda p: p.start_time.date())

# --- Join to placed_parlays for edge_pct (only ~125 rows have this) ---
con = duckdb.connect(str(DASH_DB), read_only=True)
pp = con.execute("""
  SELECT parlay_hash, home_team, away_team, combo,
         spread_line, total_line as pp_total_line, fair_odds, wz_odds,
         edge_pct, actual_size, placed_at, ticket_number
  FROM placed_parlays
  WHERE status='placed' AND edge_pct IS NOT NULL
""").fetchdf()
con.close()

def short(s): return str(s).upper().split()[-1] if pd.notna(s) else ''
pp['home_short'] = pp['home_team'].apply(short)
pp['away_short'] = pp['away_team'].apply(short)
df['home_short'] = df['home_full'].apply(lambda s: s.split()[-1] if s else '')
df['away_short'] = df['away_full'].apply(lambda s: s.split()[-1] if s else '')
pp['date_key'] = pd.to_datetime(pp['placed_at']).dt.date

# Match: home_short + away_short + slice + ou + line_type + date ±2 days
def combo_to_meta(combo):
    is_f5 = combo.startswith('F5 ')
    ou = 'over' if 'Over' in combo else 'under'
    return is_f5, ou
pp[['is_f5', 'ou']] = pp['combo'].apply(lambda c: pd.Series(combo_to_meta(c)))
pp['slice'] = pp['is_f5'].map({True: 'F5', False: 'FG'})
def pp_lt(s):
    if abs(s) < 1: return 'half(F5)'
    return 'fav-1.5' if s < 0 else 'dog+1.5'
pp['line_type'] = pp['spread_line'].apply(pp_lt)

# Build matches
df['edge_pct'] = np.nan
df['fair_odds'] = np.nan
df['wz_odds_pp'] = np.nan
df['parlay_hash'] = np.nan
for i, p in pp.iterrows():
    cand = df[(df['home_short'] == p['home_short']) & (df['away_short'] == p['away_short']) &
              (df['slice'] == p['slice']) & (df['ou'] == p['ou']) & (df['line_type'] == p['line_type'])]
    if len(cand) == 0: continue
    dd = (cand['date'].dt.date - p['date_key']).apply(lambda d: abs(d.days))
    cand = cand[dd <= 2]
    if len(cand) == 0: continue
    j = cand.index[0]
    df.at[j, 'edge_pct'] = p['edge_pct']
    df.at[j, 'fair_odds'] = p['fair_odds']
    df.at[j, 'wz_odds_pp'] = p['wz_odds']
    df.at[j, 'parlay_hash'] = p['parlay_hash']

print(f"Loaded {len(df)} bets. Matched to placed_parlays: {df['edge_pct'].notna().sum()}")
print(f"Slice counts: {df['slice'].value_counts().to_dict()}")
print(f"FG-fav-1.5: n={((df['slice']=='FG')&(df['line_type']=='fav-1.5')).sum()}")
print(f"FG-fav-1.5 + over: n={((df['slice']=='FG')&(df['line_type']=='fav-1.5')&(df['ou']=='over')).sum()}")
print(f"FG-fav-1.5 + under: n={((df['slice']=='FG')&(df['line_type']=='fav-1.5')&(df['ou']=='under')).sum()}")

df.to_parquet(OUT)
print(f"Saved {OUT}")
```

- [ ] **Step 1.3: Run extraction**

Run: `/Users/callancapitolo/NFLWork/bet_logger/venv/bin/python3 docs/superpowers/analysis/2026-05-11-fg-fav-leak/01_extract_bets.py`
Expected: prints `Loaded 627 bets. Matched to placed_parlays: ~110-125`, `FG-fav-1.5: n=222`, `FG-fav-1.5 + over: n=≥150` (most of the FG-fav are +over; under is sparse).

- [ ] **Step 1.4: Verify spread/total American odds parsed for ≥95% of FG-fav rows**

Run: `/Users/callancapitolo/NFLWork/bet_logger/venv/bin/python3 -c "import pandas as pd; d=pd.read_parquet('/tmp/fg_fav_leak/bets.parquet'); fg=d[(d.slice=='FG')&(d.line_type=='fav-1.5')]; print('total parsed:', fg.total_amer.notna().mean(), 'spread parsed:', fg.spread_amer.notna().mean())"`
Expected: both ≥ 0.95. If lower, fix the regex in `parse_spread_amer` / `parse_total_amer` before continuing.

- [ ] **Step 1.5: Commit**

```bash
git add docs/superpowers/analysis/2026-05-11-fg-fav-leak/01_extract_bets.py docs/superpowers/analysis/2026-05-11-fg-fav-leak/.gitignore docs/superpowers/plans/2026-05-11-fg-fav-leak-investigation.md
git commit -m "analysis(fg-fav-leak): extract canonical bet dataset"
```

---

### Task 2: Statistical baseline — bootstrap CI on FG-fav ROI

**Files:**
- Create: `docs/superpowers/analysis/2026-05-11-fg-fav-leak/02_baseline.py`

**What this does:** Computes bootstrapped 95% CIs on ROI for FG-fav (and reference slices) on (a) the full sample and (b) ex-bad-week. Also: binomial test on hit rate vs break-even. Outputs raw numbers we'll cite in the memo.

- [ ] **Step 2.1: Write `02_baseline.py`**

```python
#!/usr/bin/env python3
"""Bootstrap CIs + binomial tests for FG-fav and reference slices."""
import pandas as pd, numpy as np
from scipy import stats
from pathlib import Path
import json

BETS = Path('/tmp/fg_fav_leak/bets.parquet')
OUT_JSON = Path('/tmp/fg_fav_leak/baseline.json')
SEED = 42

df = pd.read_parquet(BETS)
df = df[df['result'].isin(['win', 'loss', 'push'])].copy()  # settled only

bad_mask = (df['date'] >= '2026-04-27') & (df['date'] <= '2026-05-03')

def bootstrap_roi(sub, n_boot=10000, seed=SEED):
    rng = np.random.default_rng(seed)
    stakes = sub['stake'].values
    pnls = sub['pnl'].values
    rois = np.empty(n_boot)
    for i in range(n_boot):
        idx = rng.integers(0, len(sub), len(sub))
        s = stakes[idx].sum()
        rois[i] = pnls[idx].sum() / s if s > 0 else 0
    return float(np.percentile(rois, 2.5)), float(np.percentile(rois, 97.5)), float((rois >= 0).mean())

def slice_stats(sub, label):
    n = len(sub)
    settled = sub[sub['result'].isin(['win', 'loss'])]
    wins = (settled['result'] == 'win').sum()
    losses = (settled['result'] == 'loss').sum()
    pushes = (sub['result'] == 'push').sum()
    wagered = sub['stake'].sum()
    pnl = sub['pnl'].sum()
    roi = pnl / wagered if wagered > 0 else 0
    hit_rate = wins / (wins + losses) if (wins + losses) > 0 else 0
    avg_dec = settled['dec_odds'].mean()
    breakeven = 1 / avg_dec if avg_dec and avg_dec > 0 else None
    # binomial test: is observed hit_rate worse than break-even?
    if breakeven and wins + losses > 0:
        p_lt_be = stats.binomtest(int(wins), int(wins + losses), breakeven, alternative='less').pvalue
    else:
        p_lt_be = None
    ci_lo, ci_hi, p_roi_pos = bootstrap_roi(sub)
    return {
        'label': label, 'n': n, 'wins': int(wins), 'losses': int(losses), 'pushes': int(pushes),
        'wagered': round(wagered, 2), 'pnl': round(pnl, 2), 'roi': round(roi, 4),
        'hit_rate': round(hit_rate, 4), 'break_even': round(breakeven, 4) if breakeven else None,
        'p_value_hit_lt_breakeven': round(p_lt_be, 4) if p_lt_be is not None else None,
        'ci_lo_roi': round(ci_lo, 4), 'ci_hi_roi': round(ci_hi, 4), 'p_roi_pos': round(p_roi_pos, 4),
    }

results = {}
for label, sub in [
    ('FG-fav-1.5 (all)', df[(df['slice'] == 'FG') & (df['line_type'] == 'fav-1.5')]),
    ('FG-fav-1.5 (ex-bad-week)', df[(df['slice'] == 'FG') & (df['line_type'] == 'fav-1.5') & ~bad_mask]),
    ('FG-fav-1.5 + over', df[(df['slice'] == 'FG') & (df['line_type'] == 'fav-1.5') & (df['ou'] == 'over')]),
    ('FG-fav-1.5 + under', df[(df['slice'] == 'FG') & (df['line_type'] == 'fav-1.5') & (df['ou'] == 'under')]),
    ('FG-dog+1.5', df[(df['slice'] == 'FG') & (df['line_type'] == 'dog+1.5')]),
    ('F5 half', df[(df['slice'] == 'F5') & (df['line_type'] == 'half(F5)')]),
]:
    results[label] = slice_stats(sub, label)
    print(f"\n{label}: n={results[label]['n']}, ROI={results[label]['roi']*100:.2f}%, CI=[{results[label]['ci_lo_roi']*100:.1f}, {results[label]['ci_hi_roi']*100:.1f}], P(hit < BE)={results[label]['p_value_hit_lt_breakeven']}")

with open(OUT_JSON, 'w') as f:
    json.dump(results, f, indent=2)
print(f"\nSaved {OUT_JSON}")
```

- [ ] **Step 2.2: Run baseline**

Run: `/Users/callancapitolo/NFLWork/bet_logger/venv/bin/python3 docs/superpowers/analysis/2026-05-11-fg-fav-leak/02_baseline.py`
Expected: Prints each slice's n, ROI, CI, p-value. The FG-fav-1.5 (all) row should show n≈222, ROI≈-11.5%, CI roughly [-25%, +5%], p-value for hit-rate-below-break-even likely > 0.05 (so we can't say it's statistically worse than fair on hit rate alone). Confirms what we already saw.

- [ ] **Step 2.3: Commit**

```bash
git add docs/superpowers/analysis/2026-05-11-fg-fav-leak/02_baseline.py
git commit -m "analysis(fg-fav-leak): baseline bootstrap CIs and binomial tests"
```

---

### Task 3: Hypothesis 1 — over vs under within FG-fav-`1.5`

**Mechanism tested:** If the bleed is entirely on `+over` (the model claims fav+over is positively correlated), but `+under` is fine or profitable, that's strong evidence the model overestimates fav+over correlation specifically. If both bleed equally, the issue is the spread leg or independent leg-level pricing, not correlation.

**Files:**
- Modify: `docs/superpowers/analysis/2026-05-11-fg-fav-leak/03_hypothesis_tests.py` (create as a single file, append each hypothesis function as we go)

- [ ] **Step 3.1: Create `03_hypothesis_tests.py` with H1**

```python
#!/usr/bin/env python3
"""Six hypothesis tests for the FG-fav-1.5 leak. Run as a single script.
Each H{n} function writes to results dict; final json saved at end."""
import pandas as pd, numpy as np
from scipy import stats
from pathlib import Path
import json

BETS = Path('/tmp/fg_fav_leak/bets.parquet')
OUT_JSON = Path('/tmp/fg_fav_leak/hypotheses.json')
SEED = 42

df = pd.read_parquet(BETS)
df = df[df['result'].isin(['win', 'loss', 'push'])].copy()
fg_fav = df[(df['slice'] == 'FG') & (df['line_type'] == 'fav-1.5')].copy()

results = {}

def bootstrap_roi(sub, n_boot=10000, seed=SEED):
    rng = np.random.default_rng(seed)
    s, p = sub['stake'].values, sub['pnl'].values
    out = np.empty(n_boot)
    for i in range(n_boot):
        idx = rng.integers(0, len(sub), len(sub))
        denom = s[idx].sum()
        out[i] = p[idx].sum() / denom if denom > 0 else 0
    return float(np.percentile(out, 2.5)), float(np.percentile(out, 97.5))

# ---------- H1: over vs under within FG-fav-1.5 ----------
def h1():
    out = {}
    for ou, sub in fg_fav.groupby('ou'):
        if ou == 'unk': continue
        settled = sub[sub['result'].isin(['win', 'loss'])]
        wins = int((settled['result'] == 'win').sum())
        losses = int((settled['result'] == 'loss').sum())
        if wins + losses == 0: continue
        breakeven = 1 / settled['dec_odds'].mean()
        p = stats.binomtest(wins, wins + losses, breakeven, alternative='less').pvalue
        ci_lo, ci_hi = bootstrap_roi(sub) if len(sub) >= 5 else (None, None)
        out[ou] = {
            'n': len(sub), 'wins': wins, 'losses': losses,
            'wagered': round(sub['stake'].sum(), 2), 'pnl': round(sub['pnl'].sum(), 2),
            'roi': round(sub['pnl'].sum() / sub['stake'].sum(), 4),
            'hit_rate': round(wins / (wins + losses), 4), 'break_even': round(breakeven, 4),
            'p_hit_lt_breakeven': round(p, 4),
            'ci_lo_roi': round(ci_lo, 4) if ci_lo is not None else None,
            'ci_hi_roi': round(ci_hi, 4) if ci_hi is not None else None,
        }
    # 2-proportion z-test: is over hit-rate worse than under hit-rate?
    if 'over' in out and 'under' in out and out['over']['n'] >= 5 and out['under']['n'] >= 5:
        n1, w1 = out['over']['wins'] + out['over']['losses'], out['over']['wins']
        n2, w2 = out['under']['wins'] + out['under']['losses'], out['under']['wins']
        pooled = (w1 + w2) / (n1 + n2)
        se = np.sqrt(pooled * (1 - pooled) * (1/n1 + 1/n2))
        z = (w1/n1 - w2/n2) / se if se > 0 else 0
        p = stats.norm.cdf(z)  # one-sided: over < under
        out['two_proportion_test'] = {'z': round(z, 3), 'p_value': round(p, 4)}
    return out

results['H1_over_vs_under'] = h1()
print("H1:", json.dumps(results['H1_over_vs_under'], indent=2))
```

- [ ] **Step 3.2: Run H1 standalone**

Run: `/Users/callancapitolo/NFLWork/bet_logger/venv/bin/python3 docs/superpowers/analysis/2026-05-11-fg-fav-leak/03_hypothesis_tests.py`
Expected: Two sub-rows printed (`over`, `under`). Likely outcome: over has n≈200 with ROI in the -10 to -15% range, under has n≈20 with ROI we'd need to see. If over is significantly worse than under (p < 0.10), H1 has support.

- [ ] **Step 3.3: Commit**

```bash
git add docs/superpowers/analysis/2026-05-11-fg-fav-leak/03_hypothesis_tests.py
git commit -m "analysis(fg-fav-leak): H1 over-vs-under test within FG-fav"
```

---

### Task 4: Hypothesis 2 — Calibration on the 125 placed_parlays subset

**Mechanism tested:** Direct test of "is the model's predicted hit rate accurate." We have ~30 FG-fav rows with stored `edge_pct` + `fair_odds`. Bucket by model-predicted probability and compare to realized hit rate. If predicted vs realized line lies along the 45° line, model is calibrated. If they diverge systematically (especially in high-prob bins), the model has a bias.

- [ ] **Step 4.1: Append H2 function to `03_hypothesis_tests.py`**

```python
# ---------- H2: calibration on placed_parlays subset ----------
def h2():
    sub = fg_fav[fg_fav['edge_pct'].notna() & fg_fav['fair_odds'].notna()].copy()
    if len(sub) < 10:
        return {'status': 'insufficient_data', 'n': len(sub)}
    def amer_to_dec(a):
        a = float(a)
        return 1 + a/100 if a > 0 else 1 + 100/abs(a)
    sub['fair_dec'] = sub['fair_odds'].apply(amer_to_dec)
    sub['model_prob'] = 1 / sub['fair_dec']
    sub['hit'] = sub['result'].apply(lambda r: 1 if r == 'win' else (0 if r == 'loss' else np.nan))
    settled = sub.dropna(subset=['hit']).copy()
    # bin by model_prob into quartiles
    settled['mp_q'] = pd.qcut(settled['model_prob'], q=min(4, max(2, len(settled)//5)), duplicates='drop')
    cal = settled.groupby('mp_q', observed=True).agg(
        n=('hit', 'count'),
        model_prob=('model_prob', 'mean'),
        actual_hit=('hit', 'mean'),
    ).reset_index(drop=True)
    cal['gap_pp'] = (cal['model_prob'] - cal['actual_hit']) * 100
    # Sum-of-model-prob vs observed wins (overall calibration test)
    expected = settled['model_prob'].sum()
    observed = settled['hit'].sum()
    return {
        'n_settled': len(settled),
        'expected_wins': round(float(expected), 2),
        'observed_wins': int(observed),
        'gap_wins': round(float(observed - expected), 2),
        'bins': cal.round(4).to_dict(orient='records'),
    }

results['H2_calibration'] = h2()
print("H2:", json.dumps(results['H2_calibration'], indent=2))
```

- [ ] **Step 4.2: Run + check**

Run: `/Users/callancapitolo/NFLWork/bet_logger/venv/bin/python3 docs/superpowers/analysis/2026-05-11-fg-fav-leak/03_hypothesis_tests.py`
Expected: Prints bins. If model is calibrated, gap_pp should be small (within ±5pp) across bins. If gap_pp is large positive (model > actual) especially in higher-prob bins, that's evidence of overestimation.

- [ ] **Step 4.3: Commit**

```bash
git add docs/superpowers/analysis/2026-05-11-fg-fav-leak/03_hypothesis_tests.py
git commit -m "analysis(fg-fav-leak): H2 calibration test on placed subset"
```

---

### Task 5: Hypothesis 3 — Total-line stratification

**Mechanism tested:** If the bleed is concentrated at high total lines (e.g., total ≥ 9), that's actionable: skip fav+over when total exceeds a threshold. If the bleed is uniform across totals, it's a deeper modeling issue.

- [ ] **Step 5.1: Append H3 function**

```python
# ---------- H3: total-line stratification ----------
def h3():
    sub = fg_fav[fg_fav['total_line'].notna()].copy()
    sub['tot_bucket'] = pd.cut(sub['total_line'], bins=[0, 7, 8, 9, 20], labels=['≤7', '7-8', '8-9', '9+'])
    out = []
    for bucket, g in sub.groupby('tot_bucket', observed=True):
        if len(g) == 0: continue
        settled = g[g['result'].isin(['win', 'loss'])]
        wins = int((settled['result'] == 'win').sum())
        losses = int((settled['result'] == 'loss').sum())
        if wins + losses == 0: continue
        wagered = g['stake'].sum()
        pnl = g['pnl'].sum()
        be = 1 / settled['dec_odds'].mean()
        ci_lo, ci_hi = bootstrap_roi(g) if len(g) >= 10 else (None, None)
        out.append({
            'bucket': str(bucket), 'n': len(g), 'wins': wins, 'losses': losses,
            'wagered': round(wagered, 2), 'pnl': round(pnl, 2),
            'roi': round(pnl / wagered, 4) if wagered > 0 else None,
            'hit_rate': round(wins / (wins + losses), 4), 'break_even': round(be, 4),
            'ci_lo_roi': round(ci_lo, 4) if ci_lo is not None else None,
            'ci_hi_roi': round(ci_hi, 4) if ci_hi is not None else None,
        })
    return out

results['H3_total_line_stratification'] = h3()
print("H3:", json.dumps(results['H3_total_line_stratification'], indent=2))
```

- [ ] **Step 5.2: Run + check**

Run: `/Users/callancapitolo/NFLWork/bet_logger/venv/bin/python3 docs/superpowers/analysis/2026-05-11-fg-fav-leak/03_hypothesis_tests.py`
Expected: 3-4 buckets printed. Flag any bucket where CI excludes zero (i.e., ROI confidently positive or negative).

- [ ] **Step 5.3: Commit**

```bash
git add docs/superpowers/analysis/2026-05-11-fg-fav-leak/03_hypothesis_tests.py
git commit -m "analysis(fg-fav-leak): H3 total-line stratification"
```

---

### Task 6: Hypothesis 4 — Effective WZ shave reality check

**Mechanism tested:** The model assumes `wz_dec = american_to_dec(spread) × american_to_dec(total) × 0.989` when the exact-API price is stale (>15 min). If WZ's effective shave is actually different — say, 0.95 — every fallback-path bet has its edge overstated by ~3.9%. Back-derive the shave from realized WZ tickets and check.

- [ ] **Step 6.1: Append H4 function**

```python
# ---------- H4: effective WZ shave check ----------
def h4():
    def amer_to_dec(a):
        if pd.isna(a): return None
        a = float(a)
        return 1 + a/100 if a > 0 else 1 + 100/abs(a)

    sub = fg_fav.copy()
    sub['spread_dec'] = sub['spread_amer'].apply(amer_to_dec)
    sub['total_dec'] = sub['total_amer'].apply(amer_to_dec)
    sub['independent_dec'] = sub['spread_dec'] * sub['total_dec']
    sub['effective_shave'] = sub['dec_odds'] / sub['independent_dec']
    valid = sub.dropna(subset=['effective_shave'])
    valid = valid[(valid['effective_shave'] > 0.5) & (valid['effective_shave'] < 1.1)]  # exclude pushes (dec_odds=0) and parsing errors
    return {
        'n': len(valid),
        'mean_shave': round(valid['effective_shave'].mean(), 4),
        'median_shave': round(valid['effective_shave'].median(), 4),
        'p10_shave': round(valid['effective_shave'].quantile(0.10), 4),
        'p90_shave': round(valid['effective_shave'].quantile(0.90), 4),
        'model_assumed_shave': 0.989,
        'gap': round(0.989 - valid['effective_shave'].mean(), 4),
    }

results['H4_wz_shave_check'] = h4()
print("H4:", json.dumps(results['H4_wz_shave_check'], indent=2))
```

- [ ] **Step 6.2: Run + check**

Run: `/Users/callancapitolo/NFLWork/bet_logger/venv/bin/python3 docs/superpowers/analysis/2026-05-11-fg-fav-leak/03_hypothesis_tests.py`
Expected: Mean effective shave printed. If it's close to 0.989, model assumption is fine. If it's systematically lower (e.g., 0.97), every bet's edge is overstated by ~2pp on the fallback path — that alone explains 2pp of the bleed.

- [ ] **Step 6.3: Commit**

```bash
git add docs/superpowers/analysis/2026-05-11-fg-fav-leak/03_hypothesis_tests.py
git commit -m "analysis(fg-fav-leak): H4 effective WZ shave check"
```

---

### Task 7: Hypothesis 5 — Push EV accounting

**Mechanism tested:** When the `-1.5` fav wins by exactly 1, the spread leg pushes and the parlay collapses to a single-leg total bet at a much worse implied price. The model treats `joint_prob = n_hits / n_resolved` (conditional on no push), which ignores this collapse. Quantify how much this matters.

- [ ] **Step 7.1: Append H5 function**

```python
# ---------- H5: push EV accounting ----------
def h5():
    sub = fg_fav.copy()
    pushes = sub[sub['result'] == 'push']
    wins = sub[sub['result'] == 'win']
    losses = sub[sub['result'] == 'loss']
    push_rate = len(pushes) / len(sub) if len(sub) > 0 else 0
    # Counterfactual: if pushes had been "win the single total leg" (since spread leg pushes, parlay collapses to total)
    # we'd have earned (total_dec - 1) * stake instead of 0. But that requires the total leg to have ALSO won.
    # Without per-bet leg-level outcome data we can only put an upper bound.
    # Compute hypothetical EV: assume total_leg hits at its break-even rate (1/total_dec) inside pushes.
    def amer_to_dec(a):
        if pd.isna(a): return None
        a = float(a)
        return 1 + a/100 if a > 0 else 1 + 100/abs(a)
    pushes = pushes.copy()
    pushes['total_dec'] = pushes['total_amer'].apply(amer_to_dec)
    pushes['total_breakeven_prob'] = 1 / pushes['total_dec']
    pushes['expected_payout_if_collapsed'] = (pushes['total_dec'] - 1) * pushes['stake'] * pushes['total_breakeven_prob']
    pushes['expected_loss_if_collapsed'] = -pushes['stake'] * (1 - pushes['total_breakeven_prob'])
    pushes['expected_net_if_collapsed'] = pushes['expected_payout_if_collapsed'] + pushes['expected_loss_if_collapsed']
    return {
        'n_total': len(sub),
        'n_pushes': len(pushes),
        'push_rate': round(push_rate, 4),
        'realized_push_pnl': 0.0,  # by definition
        'counterfactual_net_if_pushed_to_single': round(pushes['expected_net_if_collapsed'].sum(), 2) if len(pushes) > 0 else 0,
        'note': 'WZ pushes refund the full stake; this is the SAME as treating push as 0. The counterfactual quantifies opportunity cost if WZ instead paid the surviving leg.',
    }

results['H5_push_ev'] = h5()
print("H5:", json.dumps(results['H5_push_ev'], indent=2))
```

- [ ] **Step 7.2: Run + check**

Run: `/Users/callancapitolo/NFLWork/bet_logger/venv/bin/python3 docs/superpowers/analysis/2026-05-11-fg-fav-leak/03_hypothesis_tests.py`
Expected: Push rate around 0.022-0.025. Counterfactual likely small (a few dollars) since pushes are rare. Likely conclusion: H5 is NOT a major driver.

- [ ] **Step 7.3: Commit**

```bash
git add docs/superpowers/analysis/2026-05-11-fg-fav-leak/03_hypothesis_tests.py
git commit -m "analysis(fg-fav-leak): H5 push EV accounting"
```

---

### Task 8: Hypothesis 6 — Sample-based correlation vs realized correlation

**Mechanism tested:** The model computes `correlation_factor = joint_prob / (leg_prob_1 * leg_prob_2)` from `mlb_game_samples`. Test whether this factor is consistent with what realized 2026 MLB games actually produce, for fav+over combos specifically. This is the deepest test — directly verifies the correlation model.

**Why R:** The `mlb_game_samples` table is structured for R consumption; reading it directly with Python is doable but messy. The R helper extracts the sample distributions for the games we bet on.

- [ ] **Step 8.1: Write `04_sample_correlation.R`**

```r
#!/usr/bin/env Rscript
# Extract sample-based correlation factor for FG fav+over combos,
# for the games we placed bets on, then compare to realized outcomes.

suppressPackageStartupMessages({
  library(duckdb)
  library(dplyr)
  library(jsonlite)
})

ARGV <- commandArgs(trailingOnly = TRUE)
OUT <- if (length(ARGV) >= 1) ARGV[1] else "/tmp/fg_fav_leak/sample_corr.json"

con <- dbConnect(duckdb(), dbdir = "/Users/callancapitolo/NFLWork/Answer Keys/mlb_mm.duckdb", read_only = TRUE)
samples <- tbl(con, "mlb_game_samples") %>% collect()
dbDisconnect(con)

# For each distinct game in samples, compute:
#   P(home_margin_full > 1.5)  -- away-cover-flip, but we'll do home-cover here
#   P(total_full > 8.5)        -- canonical mid total
#   joint
#   correlation_factor

if (!"home_margin_period_Full" %in% names(samples)) {
  stop("mlb_game_samples missing home_margin_period_Full column")
}
if (!"game_total_period_Full" %in% names(samples)) {
  stop("mlb_game_samples missing game_total_period_Full column")
}

# Sample-based per-game correlation factor for fav+over (using home -1.5 as proxy for "fav")
per_game <- samples %>%
  group_by(game_id) %>%
  summarise(
    n = n(),
    p_home_cover_1_5 = mean(home_margin_period_Full > 1.5, na.rm = TRUE),
    p_away_cover_1_5 = mean(home_margin_period_Full < -1.5, na.rm = TRUE),
    p_over_8_5      = mean(game_total_period_Full > 8.5, na.rm = TRUE),
    p_home_cov_AND_over = mean(home_margin_period_Full > 1.5 & game_total_period_Full > 8.5, na.rm = TRUE),
    p_away_cov_AND_over = mean(home_margin_period_Full < -1.5 & game_total_period_Full > 8.5, na.rm = TRUE),
  ) %>%
  mutate(
    corr_factor_home_fav_over = p_home_cov_AND_over / (p_home_cover_1_5 * p_over_8_5),
    corr_factor_away_fav_over = p_away_cov_AND_over / (p_away_cover_1_5 * p_over_8_5),
  )

summary_out <- list(
  n_games = nrow(per_game),
  corr_home_fav_over = list(
    mean = mean(per_game$corr_factor_home_fav_over, na.rm = TRUE),
    median = median(per_game$corr_factor_home_fav_over, na.rm = TRUE),
    p10 = quantile(per_game$corr_factor_home_fav_over, 0.10, na.rm = TRUE),
    p90 = quantile(per_game$corr_factor_home_fav_over, 0.90, na.rm = TRUE)
  ),
  corr_away_fav_over = list(
    mean = mean(per_game$corr_factor_away_fav_over, na.rm = TRUE),
    median = median(per_game$corr_factor_away_fav_over, na.rm = TRUE),
    p10 = quantile(per_game$corr_factor_away_fav_over, 0.10, na.rm = TRUE),
    p90 = quantile(per_game$corr_factor_away_fav_over, 0.90, na.rm = TRUE)
  )
)

writeLines(toJSON(summary_out, auto_unbox = TRUE, pretty = TRUE), OUT)
cat("Wrote", OUT, "\n")
```

- [ ] **Step 8.2: Run R script**

Run: `Rscript docs/superpowers/analysis/2026-05-11-fg-fav-leak/04_sample_correlation.R`
Expected: prints `Wrote /tmp/fg_fav_leak/sample_corr.json`. Open the json — mean `corr_factor_home_fav_over` should be in the 1.05-1.20 range. If it's > 1.25, the model is systematically claiming a stronger correlation than typical.

- [ ] **Step 8.3: Append H6 function in Python to consume R output + compare to realized**

```python
# ---------- H6: sample-correlation vs realized ----------
def h6():
    import json
    try:
        sample = json.load(open('/tmp/fg_fav_leak/sample_corr.json'))
    except FileNotFoundError:
        return {'status': 'sample_corr.json missing — run 04_sample_correlation.R first'}
    # Realized "correlation factor" from FG-fav-1.5 over bets:
    # We can approximate by treating each bet's outcome as a sample of (covers, over_hits).
    # But we don't have leg-level outcomes per bet — only joint outcome (win = both hit).
    # So we estimate realized P(joint hit) and compare to model's expected joint.
    # Realized P(joint hit) ≈ wins / settled for FG-fav-over.
    fg_fav_over = fg_fav[fg_fav['ou'] == 'over']
    settled = fg_fav_over[fg_fav_over['result'].isin(['win', 'loss'])]
    realized_joint = (settled['result'] == 'win').mean() if len(settled) > 0 else None
    # Model's implied joint = 1 / fair_dec, averaged across the matched subset
    matched = fg_fav_over[fg_fav_over['fair_odds'].notna()].copy()
    def amer_to_dec(a):
        a = float(a)
        return 1 + a/100 if a > 0 else 1 + 100/abs(a)
    if len(matched) > 0:
        matched['model_joint'] = matched['fair_odds'].apply(lambda a: 1 / amer_to_dec(a))
        model_joint_mean = float(matched['model_joint'].mean())
    else:
        model_joint_mean = None
    return {
        'sample_corr_factor_home_fav_over_mean': sample.get('corr_home_fav_over', {}).get('mean'),
        'sample_corr_factor_away_fav_over_mean': sample.get('corr_away_fav_over', {}).get('mean'),
        'realized_joint_hit_prob_fav_over': round(realized_joint, 4) if realized_joint is not None else None,
        'model_predicted_joint_hit_prob_fav_over (matched subset)': round(model_joint_mean, 4) if model_joint_mean else None,
        'gap_pp_model_minus_realized': round((model_joint_mean - realized_joint) * 100, 2) if (model_joint_mean and realized_joint) else None,
    }

results['H6_sample_corr_vs_realized'] = h6()
print("H6:", json.dumps(results['H6_sample_corr_vs_realized'], indent=2))

# ---------- WRITE ALL RESULTS ----------
import json as _json
with open(OUT_JSON, 'w') as f:
    _json.dump(results, f, indent=2)
print(f"\nAll hypothesis results saved to {OUT_JSON}")
```

- [ ] **Step 8.4: Run full hypothesis suite**

Run: `/Users/callancapitolo/NFLWork/bet_logger/venv/bin/python3 docs/superpowers/analysis/2026-05-11-fg-fav-leak/03_hypothesis_tests.py`
Expected: All six hypotheses print without errors. `/tmp/fg_fav_leak/hypotheses.json` is written.

- [ ] **Step 8.5: Commit**

```bash
git add docs/superpowers/analysis/2026-05-11-fg-fav-leak/04_sample_correlation.R docs/superpowers/analysis/2026-05-11-fg-fav-leak/03_hypothesis_tests.py
git commit -m "analysis(fg-fav-leak): H6 sample-corr vs realized"
```

---

### Task 9: Synthesize findings into memo

**Files:**
- Create: `docs/superpowers/analysis/2026-05-11-fg-fav-leak/05_write_memo.py`
- Create: `docs/superpowers/analysis/2026-05-11-fg-fav-leak/findings.md` (output)

- [ ] **Step 9.1: Write `05_write_memo.py`**

```python
#!/usr/bin/env python3
"""Read baseline + hypotheses JSON, write a memo with ranked findings."""
import json
from pathlib import Path

ROOT = Path('/tmp/fg_fav_leak')
baseline = json.load(open(ROOT / 'baseline.json'))
hyp = json.load(open(ROOT / 'hypotheses.json'))
OUT = Path(__file__).parent / 'findings.md'

lines = []
lines.append("# FG `-1.5` Favorite Correlated Parlay Leak — Findings\n")
lines.append(f"**Date:** 2026-05-11\n")
lines.append("**Source data:** 627 MLB correlated parlays (2026-03-28 → 2026-05-10) from `bet_tracking.xlsx` joined to dashboard `placed_parlays` for stored model edge predictions.\n")

# ---- Baseline ----
lines.append("## Baseline\n")
bf = baseline['FG-fav-1.5 (all)']
bx = baseline['FG-fav-1.5 (ex-bad-week)']
lines.append(f"- **FG-fav-1.5 (all):** n={bf['n']}, ROI = {bf['roi']*100:.2f}%, 95% CI [{bf['ci_lo_roi']*100:.1f}%, {bf['ci_hi_roi']*100:.1f}%], P(ROI ≥ 0) = {bf['p_roi_pos']}")
lines.append(f"- **FG-fav-1.5 (ex-bad-week):** n={bx['n']}, ROI = {bx['roi']*100:.2f}%, 95% CI [{bx['ci_lo_roi']*100:.1f}%, {bx['ci_hi_roi']*100:.1f}%]")
for ref in ['FG-dog+1.5', 'F5 half']:
    r = baseline[ref]
    lines.append(f"- **{ref}:** n={r['n']}, ROI = {r['roi']*100:.2f}%, 95% CI [{r['ci_lo_roi']*100:.1f}%, {r['ci_hi_roi']*100:.1f}%]")
lines.append("")

# ---- Hypothesis verdicts ----
lines.append("## Hypothesis Verdicts\n")

# H1
h1 = hyp['H1_over_vs_under']
o = h1.get('over', {}); u = h1.get('under', {})
twop = h1.get('two_proportion_test', {})
lines.append("### H1 — Over vs Under within FG-fav-1.5")
lines.append(f"- `+over`: n={o.get('n')}, ROI = {o.get('roi', 0)*100:.2f}%, CI [{o.get('ci_lo_roi', 0)*100:.1f}%, {o.get('ci_hi_roi', 0)*100:.1f}%]")
lines.append(f"- `+under`: n={u.get('n')}, ROI = {u.get('roi', 0)*100:.2f}%, CI [{u.get('ci_lo_roi', 0)*100:.1f}%, {u.get('ci_hi_roi', 0)*100:.1f}%]" if u else "- `+under`: insufficient data")
if twop:
    lines.append(f"- 2-proportion test (over hit < under hit): z={twop['z']}, p={twop['p_value']}")
lines.append("")

# H2
h2 = hyp['H2_calibration']
if h2.get('status') == 'insufficient_data':
    lines.append(f"### H2 — Calibration: INSUFFICIENT DATA (n={h2['n']})\n")
else:
    lines.append("### H2 — Calibration on placed_parlays subset")
    lines.append(f"- n_settled = {h2['n_settled']}, expected wins = {h2['expected_wins']}, observed wins = {h2['observed_wins']} (gap = {h2['gap_wins']:+.2f})")
    lines.append("- By model_prob quartile:")
    for b in h2['bins']:
        lines.append(f"  - n={b['n']}, mean model_prob={b['model_prob']:.3f}, actual hit={b['actual_hit']:.3f}, gap={b['gap_pp']:+.1f} pp")
    lines.append("")

# H3
h3 = hyp['H3_total_line_stratification']
lines.append("### H3 — Total-line stratification")
for b in h3:
    ci = f"[{b['ci_lo_roi']*100:.1f}%, {b['ci_hi_roi']*100:.1f}%]" if b.get('ci_lo_roi') is not None else "(n too small for CI)"
    lines.append(f"- total {b['bucket']}: n={b['n']}, ROI={b['roi']*100:.2f}%, CI {ci}")
lines.append("")

# H4
h4 = hyp['H4_wz_shave_check']
lines.append("### H4 — Effective WZ shave reality check")
lines.append(f"- Sample n = {h4['n']}, model assumes shave = {h4['model_assumed_shave']}")
lines.append(f"- Mean realized shave = {h4['mean_shave']}, median = {h4['median_shave']}, p10 = {h4['p10_shave']}, p90 = {h4['p90_shave']}")
lines.append(f"- Gap (model assumption minus realized mean) = {h4['gap']:+.4f} (positive → model overestimates parlay value)")
lines.append("")

# H5
h5 = hyp['H5_push_ev']
lines.append("### H5 — Push EV accounting")
lines.append(f"- n_pushes = {h5['n_pushes']}, push_rate = {h5['push_rate']:.4f}")
lines.append(f"- Counterfactual net if WZ paid the surviving leg = ${h5.get('counterfactual_net_if_pushed_to_single', 0):.2f}")
lines.append("- Note: WZ actually refunds the full stake on a leg push (no payout), so this is OPPORTUNITY cost, not realized loss.")
lines.append("")

# H6
h6 = hyp['H6_sample_corr_vs_realized']
lines.append("### H6 — Sample correlation vs realized")
lines.append(f"- Sample corr_factor (home fav + over, mean across games): {h6.get('sample_corr_factor_home_fav_over_mean')}")
lines.append(f"- Sample corr_factor (away fav + over, mean): {h6.get('sample_corr_factor_away_fav_over_mean')}")
lines.append(f"- Model-predicted joint prob (matched subset): {h6.get('model_predicted_joint_hit_prob_fav_over (matched subset)')}")
lines.append(f"- Realized joint hit rate: {h6.get('realized_joint_hit_prob_fav_over')}")
lines.append(f"- Gap (model - realized, pp): {h6.get('gap_pp_model_minus_realized')}")
lines.append("")

lines.append("## Ranking (highest evidence first)\n")
lines.append("_To be filled in after running. Use this template:_\n")
lines.append("1. **<Hypothesis>** — <evidence summary> — <recommended action>")
lines.append("2. ...")
lines.append("")

lines.append("## Recommended Action\n")
lines.append("_To be filled in after ranking._\n")
lines.append("Candidates:")
lines.append("- Skip FG-fav-1.5 + over when total ≥ X (if H3 finds a confident high-total leak)")
lines.append("- Cut Kelly multiplier on FG-fav-1.5 + over by 0.5x (broad shrinkage)")
lines.append("- Recalibrate WZ shave constant in `mlb_correlated_parlay.R:35` (if H4 finds gap)")
lines.append("- Recalibrate correlation_factor with shrinkage prior (if H6 finds gap)")
lines.append("- No action; continue collecting data (if all tests fail to reject the null)")
lines.append("")

OUT.write_text('\n'.join(lines))
print(f"Wrote {OUT}")
```

- [ ] **Step 9.2: Run memo generator**

Run: `/Users/callancapitolo/NFLWork/bet_logger/venv/bin/python3 docs/superpowers/analysis/2026-05-11-fg-fav-leak/05_write_memo.py`
Expected: prints `Wrote .../findings.md`. The memo will have all baseline + hypothesis numbers but the Ranking and Recommended Action sections still as templates — those need human judgment.

- [ ] **Step 9.3: Read the memo + fill in Ranking and Recommended Action sections**

Open `findings.md` and:
1. Rank the six hypotheses by evidence strength (which had the most lopsided result?). Use criteria: confidence interval excludes zero, p-value < 0.10, magnitude of gap.
2. Pick one or two recommended actions backed by the data.
3. Quote specific numbers in support of each recommendation.

- [ ] **Step 9.4: Commit final memo**

```bash
git add docs/superpowers/analysis/2026-05-11-fg-fav-leak/05_write_memo.py docs/superpowers/analysis/2026-05-11-fg-fav-leak/findings.md
git commit -m "analysis(fg-fav-leak): findings memo + recommendation"
```

---

### Task 10: Update memory + (conditionally) CLAUDE.md

**Files:**
- Modify: `/Users/callancapitolo/.claude-personal/projects/-Users-callancapitolo-NFLWork/memory/mlb_parlay_edge_overestimation.md`
- Conditionally modify: `/Users/callancapitolo/NFLWork/Answer Keys/CLAUDE.md` (only if a hypothesis tested conclusive)

- [ ] **Step 10.1: Rewrite memory file with final verdict**

Replace the existing content with a section structured as: "What we tested → What we found → What changed (if anything)". The current memory file says "model calibration is perfect on small sample" — update to whatever the deeper analysis actually shows.

- [ ] **Step 10.2: If H3 found a confident high-total leak (CI excludes zero), add to Answer Keys/CLAUDE.md**

Under a new heading `## Known model biases`, add a bullet:
> - **FG fav `-1.5` + over at total ≥ X is structurally -EV** (n=Y bets, ROI=Z%, 95% CI [a,b]). The model's correlation_factor for this combo is uncalibrated to current MLB run-scoring. Skip via filter in `mlb_correlated_parlay.R` or apply Kelly multiplier 0.5x for this sub-slice.

If no hypothesis was conclusive, skip this step.

- [ ] **Step 10.3: Commit**

```bash
git add /Users/callancapitolo/.claude-personal/projects/-Users-callancapitolo-NFLWork/memory/mlb_parlay_edge_overestimation.md
# Add Answer Keys/CLAUDE.md only if modified
if git diff --quiet "Answer Keys/CLAUDE.md"; then
  echo "CLAUDE.md not modified"
else
  git add "Answer Keys/CLAUDE.md"
fi
git commit -m "analysis(fg-fav-leak): update memory and known-bias docs"
```

---

### Task 11: Present findings + decide on merge

- [ ] **Step 11.1: Render `findings.md` inline in conversation** (per user preference — never use `open` or external viewer).

- [ ] **Step 11.2: Ask user:** "Merge `analysis/fg-fav-leak-2026-05-11` to `main` so the memo is preserved in repo history, or leave as a branch?"

- [ ] **Step 11.3: If merge approved:**

```bash
cd /Users/callancapitolo/NFLWork
git checkout main
git merge --no-ff analysis/fg-fav-leak-2026-05-11
git worktree remove .worktrees/fg-fav-leak
git branch -d analysis/fg-fav-leak-2026-05-11
```

NEVER push to remote without explicit user approval beyond this.

---

## Stop Conditions / Honest Outs

This investigation should **stop and report inconclusive** if:
- Step 1.3 yields < 200 FG-fav-1.5 bets (data corruption — investigate before proceeding).
- All six hypotheses produce CIs that include zero and p-values > 0.10. In that case: the bleed is statistically indistinguishable from variance at this sample size. The honest output is "need ~2x more bets to resolve; here's the threshold N for statistical power 0.8."

This investigation should **flag a high-priority code fix** if:
- H4 finds effective WZ shave < 0.95 (i.e., gap > 0.04) — that's a direct sizing bug.
- H6 finds model-predicted joint > realized joint by > 5pp with confident CI — that's correlation overstatement, needs a shrinkage prior.
- H3 finds one total-line bucket with CI fully below zero — actionable filter.

Anything in between (gap exists but CI overlaps zero) → recommend Kelly multiplier shrinkage on the sub-slice, no code change beyond config.

---

## Expected Total Time

- Task 1: 30 min (data extraction has tricky regex parsing)
- Tasks 2-8: 15-20 min each (~2 hours)
- Task 9: 30 min (memo synthesis, requires human judgment on ranking)
- Task 10: 15 min
- Task 11: 10 min

**Total: ~3.5 hours.**
