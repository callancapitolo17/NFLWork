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
