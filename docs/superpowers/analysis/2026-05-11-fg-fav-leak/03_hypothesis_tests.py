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
    # bin by model_prob into quartiles (or fewer if N small)
    n_bins = min(4, max(2, len(settled) // 5))
    settled['mp_q'] = pd.qcut(settled['model_prob'], q=n_bins, duplicates='drop')
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
print("H2:", json.dumps(results['H2_calibration'], indent=2, default=str))

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

# Save partial results (will be overwritten as more H{n} are added in later tasks)
with open(OUT_JSON, 'w') as f:
    json.dump(results, f, indent=2)
print(f"\nSaved partial results to {OUT_JSON}")
