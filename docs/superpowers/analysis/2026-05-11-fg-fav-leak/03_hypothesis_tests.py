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
    else:
        out['two_proportion_test'] = {'status': 'skipped — under arm has insufficient sample (n_under < 5)'}
    return out

results['H1_over_vs_under'] = h1()
print("H1:", json.dumps(results['H1_over_vs_under'], indent=2))

# Save partial results (will be overwritten as more H{n} are added in later tasks)
with open(OUT_JSON, 'w') as f:
    json.dump(results, f, indent=2)
print(f"\nSaved partial results to {OUT_JSON}")
