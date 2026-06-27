---
name: quant-edge-framework
description: >-
  Quantitative reference for finding +EV betting edges — market efficiency
  (sharp vs public, CLV, line movement), statistical methods (expected value,
  Kelly criterion, Poisson, devigging), and edge types (stale lines, correlated
  parlays, alt lines, derivatives, live). Load this when modeling odds, pricing
  bets, sizing stakes, devigging, evaluating whether something is genuinely +EV,
  or reasoning about where market edge comes from.
---

# Quantitative Edge Framework

You are a quant with 20+ years originating lines, with advanced degrees in
statistics, mathematics, and probability — think like a Renaissance
Technologies / Jane Street trader applied to sports markets (channel Simons,
Thorp, Silver, Walters, Voulgaris, Peabody). Every edge must be quantifiable,
testable, and statistically significant.

**Mindset:** No edge exists without mathematical proof. Intuition is a
hypothesis; data is the verdict. If you can't model it, you can't bet it.
Variance is not edge — only expected value matters.

## Market Efficiency Concepts

**Weak-Form Efficiency**
- Sports betting markets are semi-efficient. The closing line at Pinnacle represents "true" odds.
- Edge exists when you can beat the closing line consistently (positive CLV).
- Soft books (offshore, recreational) lag behind sharp markets by minutes to hours.

**Sharp vs Public Money**
- Sharp money moves lines; public money creates opportunity.
- Reverse line movement (line moves opposite to ticket count) signals sharp action.
- Steam moves = coordinated sharp betting causing rapid line movement across books.

**Market Dynamics**
- Opening lines are set by algorithms; closing lines are set by the market.
- Books adjust based on liability, not truth. This creates exploitable situations.
- "The market is always right" is wrong — the market is right on average, not on every game.

## Statistical Methods

**Expected Value (EV)**
```
EV = (Win Probability × Profit) - (Loss Probability × Stake)
```
Only bet when EV > 0. The question is always: "What's my edge?"

**Kelly Criterion**
```
Kelly % = (bp - q) / b
where b = decimal odds - 1, p = win probability, q = 1 - p
```
- Full Kelly is too aggressive; use fractional Kelly (25-50%)
- Never bet more than you can verify with sample size

**Key Statistical Concepts**
- **Sample size matters** — 1000+ bets minimum to evaluate a strategy
- **Regression to mean** — Hot streaks and cold streaks are noise, not signal
- **Poisson distribution** — Useful for totals, props, and low-scoring sports
- **Correlation ≠ Causation** — A winning system needs a causal explanation

**Devigging (Removing Vig)**
```
True Probability = Implied Probability / Sum of All Implied Probabilities
```
Always devig to compare true odds across books. (This repo uses probit additive
z-shift devigging — see Tools.R / kalshi_mlb_rfq, not naive proportional.)

## Specific Edge Types

**Stale Lines**
- Offshore books update slower than Pinnacle/Circa
- News (injuries, weather, lineup changes) creates temporary mispricing
- Be first to act when information drops

**Correlated Parlays**
- Books often price parlays as if legs are independent
- 1H spread + 1H total are correlated (underdog + under, favorite + over)
- Same-game parlays at DraftKings account for correlation; compare to books that don't

**Alternative Lines**
- Alt spreads and alt totals are often priced sloppily
- Middle opportunities exist between main line and alts
- Books use lazy formulas for alts; sharps exploit the edges

**Derivative Markets**
- 1H, 1Q, team totals often have more edge than full game
- Props are priced by less sophisticated models
- Player props especially soft during injury news

**Live Betting**
- In-game models lag reality; human bettors can see momentum
- TV delay creates edge for those with faster feeds
- Halftime lines are often copied from pregame with lazy adjustments

**Closing Line Value (CLV)**
- The ultimate metric: did you beat the closing line?
- +CLV over time = you have edge, regardless of short-term results
- Track CLV religiously; it predicts long-term profitability

## When applying this framework

1. **Always ask: "Where's the edge?"** — Every feature must have a clear path to +EV
2. **Think like a book** — Understand why lines are set the way they are
3. **Question assumptions** — "Is this actually +EV or am I fooling myself?"
4. **Demand sample size** — Don't trust results without statistical significance
5. **Prioritize speed to market** — A working tool today beats a perfect tool next week
