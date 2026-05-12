#!/usr/bin/env python3
"""Read baseline + hypotheses JSON, emit complete findings.md."""
import json
from pathlib import Path

ROOT = Path('/tmp/fg_fav_leak')
baseline = json.load(open(ROOT / 'baseline.json'))
hyp = json.load(open(ROOT / 'hypotheses.json'))
OUT = Path(__file__).parent / 'findings.md'

def f_roi(x):
    if x is None:
        return "n/a"
    return f"{x*100:+.2f}%"

def f_pct(x):
    if x is None:
        return "n/a"
    return f"{x*100:.2f}%"

def f_pp(x):
    if x is None:
        return "n/a"
    return f"{x:+.1f} pp"

lines = []
lines.append("# FG `-1.5` Favorite Correlated Parlay Leak — Findings\n")
lines.append("**Date:** 2026-05-11  ")
lines.append("**Branch:** `analysis/fg-fav-leak-2026-05-11`  ")
lines.append("**Data:** 627 MLB correlated parlays (2026-03-28 → 2026-05-10) from `bet_tracking.xlsx`, joined to dashboard `placed_parlays` for stored model `edge_pct` where available (~119 matched).\n")
lines.append("## TL;DR\n")
lines.append("The FG `-1.5` favorite + over leak (-$1,575, -11.6% ROI on 222 bets) is concentrated almost entirely in one sub-slice: **bets where the total line ≤ 7.0** (n=33, ROI -73%, bootstrap 95% CI [-97%, -38%] fully below zero). The other three total-line buckets (7-8, 8-9, 9+) are statistically indistinguishable from break-even.\n")
lines.append("The mechanism is **NOT** correlation factor overstatement (H6: model corr is mild, ~1.04 at low totals). It is **marginal-probability overestimation** of the over leg at low totals: the simulator predicts ~28% joint hit rate at total_line=6.5, but realized rate at total ≤7 is 9.4% — a 19pp gap. Most likely cause: the resampled distribution overstates how often games with low-total lines go over, or WZ's over line at low totals is sharper than the model believes.\n")
lines.append("**Recommended action:** Skip FG `-1.5` favorite + over parlays when the total line is ≤ 7.0. On the observed sample this would have eliminated -$1,721 of losses while leaving the other +EV combos untouched. Plus follow-up: investigate why the simulator overstates over probability at low totals.\n")

# ---- Baseline ----
lines.append("## Baseline (Task 2)\n")
lines.append("| Slice | n | ROI | 95% CI ROI | Hit Rate | Break-even |")
lines.append("|---|---:|---:|---|---:|---:|")
for k in ['FG-fav-1.5 (all)', 'FG-fav-1.5 (ex-bad-week)', 'FG-fav-1.5 + over', 'FG-fav-1.5 + under', 'FG-dog+1.5', 'F5 half']:
    r = baseline[k]
    lines.append(f"| {k} | {r['n']} | {f_roi(r['roi'])} | [{f_roi(r['ci_lo_roi'])}, {f_roi(r['ci_hi_roi'])}] | {f_pct(r['hit_rate'])} | {f_pct(r['break_even'])} |")
lines.append("\nNote: hit rate is slightly **above** break-even for FG-fav-1.5 — the loss is in the odds distribution, not raw win frequency.\n")

# ---- H1 ----
h1 = hyp['H1_over_vs_under']
o = h1['over']
u = h1['under']
lines.append("## H1 — Over vs Under within FG-fav-1.5\n")
lines.append(f"- `+over`: n={o['n']}, ROI {f_roi(o['roi'])}, hit_rate {f_pct(o['hit_rate'])}, CI [{f_roi(o['ci_lo_roi'])}, {f_roi(o['ci_hi_roi'])}]")
lines.append(f"- `+under`: n={u['n']} — **insufficient sample** ({u['n']} settled).")
two_prop = h1.get('two_proportion_test', {})
two_prop_status = two_prop.get('status', 'see JSON') if isinstance(two_prop, dict) else 'see JSON'
lines.append(f"- 2-proportion test: {two_prop_status}\n")
lines.append("**Verdict:** Moot. The under arm of FG-fav-1.5 is essentially never bet (the model rarely flags it as +EV). The entire FG-fav-1.5 leak universe is the over arm.\n")

# ---- H2 ----
h2 = hyp['H2_calibration']
lines.append("## H2 — Calibration on placed_parlays subset\n")
if h2.get('status') == 'insufficient_data':
    lines.append(f"INSUFFICIENT DATA (n={h2.get('n', 'n/a')}).\n")
else:
    lines.append(f"- n_settled = {h2['n_settled']}, expected wins = {h2['expected_wins']}, observed wins = {h2['observed_wins']}, gap = {h2['gap_wins']:+.2f}")
    lines.append("- By model_prob quartile:")
    for b in h2['bins']:
        lines.append(f"  - n={b['n']}: model_prob {f_pct(b['model_prob'])}, actual hit {f_pct(b['actual_hit'])}, gap {f_pp(b['gap_pp'])}")
    lines.append("")
    lines.append("**Verdict:** On the small matched subset (n=30), the model is well-calibrated overall (7.06 expected vs 7 observed wins). Per-bin gaps are noisy but no monotonic over- or under-prediction pattern. The leak is NOT a calibration miss on probability.\n")

# ---- H3 ----
lines.append("## H3 — Total-line stratification (SMOKING GUN)\n")
lines.append("| total bucket | n | wins | losses | ROI | hit rate | break-even | 95% CI ROI |")
lines.append("|---|---:|---:|---:|---:|---:|---:|---|")
for b in hyp['H3_total_line_stratification']:
    ci = f"[{f_roi(b['ci_lo_roi'])}, {f_roi(b['ci_hi_roi'])}]" if b.get('ci_lo_roi') is not None else "(n<10 — no CI)"
    lines.append(f"| {b['bucket']} | {b['n']} | {b['wins']} | {b['losses']} | {f_roi(b['roi'])} | {f_pct(b['hit_rate'])} | {f_pct(b['break_even'])} | {ci} |")
lines.append("\n**Verdict:** The `≤7` bucket is catastrophic and its CI does not include zero — the leak is real, not variance. The other buckets are noisy break-even. This single bucket accounts for **more than 100%** of the FG-fav-1.5 loss.\n")

# ---- H4 ----
h4 = hyp['H4_wz_shave_check']
lines.append("## H4 — Effective WZ shave reality check\n")
lines.append(f"- n = {h4['n']}, model assumes shave = {h4['model_assumed_shave']}")
lines.append(f"- Mean realized shave = {h4['mean_shave']}, median = {h4['median_shave']}, p10 = {h4['p10_shave']}, p90 = {h4['p90_shave']}")
lines.append(f"- Gap (model - realized mean) = {h4['gap']:+.4f}")
lines.append("\n**Verdict:** Shave matches model assumption almost exactly (median 0.9909 vs assumed 0.989). NOT the leak.\n")

# ---- H5 ----
h5 = hyp['H5_push_ev']
lines.append("## H5 — Push EV accounting\n")
lines.append(f"- n_pushes = {h5['n_pushes']}, push_rate = {h5['push_rate']:.4f}")
lines.append(f"- Counterfactual net if WZ collapsed pushes to single leg: ${h5['counterfactual_net_if_pushed_to_single']:.2f}")
lines.append("\n**Verdict:** Push rate is only 2.3%. WZ refunds full stake on push. Opportunity cost from push-collapse is essentially zero. NOT the leak.\n")

# ---- H6 ----
h6 = hyp['H6_sample_corr_vs_realized']
lines.append("## H6 — Sample correlation vs realized\n")
lines.append("Sample-derived correlation factor and predicted joint probability by total line (from `mlb_game_samples`):")
lines.append("| total_line | mean corr (home fav) | mean corr (away fav) | predicted joint (home) | predicted joint (away) |")
lines.append("|---:|---:|---:|---:|---:|")
for r in h6['sample_corr_by_total_line']:
    lines.append(f"| {r['total_line']} | {r['mean_corr_home']:.3f} | {r['mean_corr_away']:.3f} | {r['mean_joint_home']:.3f} | {r['mean_joint_away']:.3f} |")
lines.append("")
lines.append("Realized joint hit rate by total bucket (from our actual bets):")
lines.append("| total bucket | n | realized joint hit rate |")
lines.append("|---|---:|---:|")
for r in h6['realized_joint_hit_prob_fav_over_by_total']:
    lines.append(f"| {r['bucket']} | {r['n']} | {f_pct(r['realized_joint_hit_prob'])} |")
lines.append("")
lines.append("**Critical comparison:** at total_line=6.5 the model's predicted joint hit prob (samples) is ~28%; realized rate at total ≤7 (our bets) is only 9.4%. Gap: ~19 pp.\n")
lines.append("**Verdict:** Correlation factor (1.04-1.15 range) is mild and consistent — NOT overstated. The leak is in the **underlying marginal probability** the simulator assigns to OVER at low totals. Either the resampled distribution overstates how often low-scoring games go over their (low) totals, or WZ's over line at total ≤7 is sharper than the simulator expects.\n")

# ---- Ranking ----
lines.append("## Ranking of Hypotheses by Evidence Strength\n")
lines.append("1. **H3 (total-line stratification): CONFIRMED.** CI [-97%, -38%] for the ≤7 bucket excludes zero. n=33 is modest but the effect size is large enough that this is signal, not noise. This is the only hypothesis with a CI fully on one side of zero.")
lines.append("2. **H6 (marginal probability mismatch at low totals): SUPPORTING.** Model-predicted joint at total=6.5 (~28%) is ~3× higher than realized (~9.4%) at total ≤7. Cannot quantify per-bet gap rigorously without re-running the pricer, but the directional evidence is consistent with H3 and identifies the mechanism (over-leg marginal, not correlation).")
lines.append("3. **H4, H5: ELIMINATED.** Shave matches model; push EV is negligible. Both ruled out as drivers.")
lines.append("4. **H1: MOOT.** Under arm has 1 bet; cannot test.")
lines.append("5. **H2: NEUTRAL.** Model calibration on n=30 is fine overall, but this subset over-represents games where the bet was placed (selection effect). Doesn't contradict H3/H6.\n")

# ---- Recommended Action ----
lines.append("## Recommended Action\n")
lines.append("**Primary recommendation:** Skip FG `-1.5` favorite + over parlays when total_line ≤ 7.0.\n")
lines.append("Implementation: add a filter in `Answer Keys/mlb_correlated_parlay.R`'s combo-construction loop that drops Home Spread + Over / Away Spread + Over combos where `row$total_line <= 7.0` for FG. (Spread + Under and dog combos remain eligible.)\n")
lines.append("**Expected impact on the observed sample:** -$1,721 of losses avoided (the ≤7 bucket's P&L), while leaving the other three buckets — net +$659 — untouched. Net P&L improvement on this slice: +$1,721 if perfectly executed.\n")
lines.append("**Caveats / known limits:**")
lines.append("- n=33 in the ≤7 bucket is small. The point estimate -73% ROI is noisy; CI [-97%, -38%] is wide. The directional conclusion (sub-slice is unprofitable) is robust, but the precise magnitude is uncertain. Re-evaluate after another month of data.")
lines.append("- The mechanism (H6 marginal-probability mismatch) implies that the simulator's run-distribution at low-total games is wrong. A more thorough fix would investigate why and recalibrate the underlying samples — but the filter recommendation is the immediate +EV move.")
lines.append("- Generalization: it is possible other slices (F5, FG-dog) also have a low-total bleed but we didn't subset-test them. Worth a follow-up stratification.")
lines.append("- This analysis is conditional on the WZ shave assumption holding (H4 confirmed it does, within 0.4pp).\n")
lines.append("**Secondary recommendation:** Investigate the simulator's over-distribution at low-total games. Specifically:")
lines.append("- For games where MLB Odds API posts a total ≤ 7.0, does `mlb_game_samples` predict over% above the over leg's implied break-even prob from WZ's line?")
lines.append("- If yes, the upstream sample generation is biased toward higher-scoring outcomes than current MLB delivers. Possible fixes: (a) restrict historical samples to 2024+ for low-total games, (b) introduce a total-line conditional adjustment in the resampler, (c) skip parlay generation when the simulator's over-prob disagrees with the line by >X%.\n")

# ---- Stop conditions reached ----
lines.append("## Stop Conditions Reached\n")
lines.append("Per the plan's stop conditions:")
lines.append("- H3 found one total-line bucket with CI fully below zero → **actionable filter triggered** (skip total ≤ 7).")
lines.append("- H4 + H5 ruled out as drivers (shave OK, push EV ~$0).")
lines.append("- H6 explains the mechanism (over-marginal mispriced at low totals) without contradicting H3.")
lines.append("- No hypothesis tested inconclusive across the board, so the \"need 2x more bets\" exit condition is not met.")

OUT.write_text('\n'.join(lines))
print(f"Wrote {OUT}")
