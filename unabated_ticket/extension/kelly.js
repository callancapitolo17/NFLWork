// Kelly math for the Unabated Ticket panel. Pure functions, no DOM.
//
// Mode B of the Kelly Calculator sheet: Unabated already publishes the edge
// (EV per $1 staked) for every line, so the stake is sized straight from it:
//   full Kelly fraction = edge / (decimal_book - 1)
// (derivation: p_fair = (1 + edge) / dec, f = (p·b − q)/b with b = dec − 1).
//
// Loaded two ways: as a plain <script> in panel.html (exposes
// globalThis.UnabatedKelly) and via require() in tests/kelly.test.js.

(function (root) {
  "use strict";

  function assertAmerican(americanOdds, label) {
    if (typeof americanOdds !== "number" || !Number.isFinite(americanOdds)) {
      throw new Error(`${label}: expected a finite American price, got ${americanOdds}`);
    }
    if (Math.abs(americanOdds) < 100) {
      throw new Error(`${label}: American price must be <= -100 or >= 100, got ${americanOdds}`);
    }
  }

  function americanToDecimal(americanOdds) {
    assertAmerican(americanOdds, "americanToDecimal");
    if (americanOdds > 0) return 1 + americanOdds / 100;
    return 1 + 100 / Math.abs(americanOdds);
  }

  function americanToProb(americanOdds) {
    assertAmerican(americanOdds, "americanToProb");
    if (americanOdds > 0) return 100 / (americanOdds + 100);
    const magnitude = Math.abs(americanOdds);
    return magnitude / (magnitude + 100);
  }

  // Implied probability for DISPLAY (prediction-market cents). Exchanges publish
  // a probability (sourceFormat 4) or decimal (2) that Unabated rounds into a
  // whole American `price`; prefer the exact source so the cents match the screen.
  function bookProbOf({ bookPrice, sourceFormat, sourcePrice }) {
    if (sourceFormat === 4 && sourcePrice > 0 && sourcePrice < 1) return sourcePrice;
    if (sourceFormat === 2 && sourcePrice > 1) return 1 / sourcePrice;
    return americanToProb(bookPrice);
  }

  function assertPositiveNumber(value, label) {
    if (typeof value !== "number" || !Number.isFinite(value) || value <= 0) {
      throw new Error(`${label}: expected a positive number, got ${value}`);
    }
  }

  // Stake from Unabated's edge % (e.g. 1.89 means +1.89% per $1) and the
  // American book price. Dollars are NOT rounded (user decision 2026-09-08).
  // Uses the American price, not the exchange source price, so the stake is
  // consistent with the edge Unabated computed from that same American price.
  function kellyStakeFromEdge({ bookPrice, edgePct, bankroll, multiplier }) {
    assertPositiveNumber(bankroll, "bankroll");
    assertPositiveNumber(multiplier, "multiplier");
    if (typeof edgePct !== "number" || !Number.isFinite(edgePct)) {
      throw new Error(`edgePct: expected a finite number, got ${edgePct}`);
    }
    const netOdds = americanToDecimal(bookPrice) - 1;
    const edge = edgePct / 100;
    const fraction = edge > 0 ? edge / netOdds : 0;
    return { stake: bankroll * fraction * multiplier, fullKellyFraction: fraction, edge };
  }

  const api = { americanToDecimal, americanToProb, bookProbOf, kellyStakeFromEdge };

  if (typeof module !== "undefined" && module.exports) {
    module.exports = api;
  } else {
    root.UnabatedKelly = api;
  }
})(typeof globalThis !== "undefined" ? globalThis : this);
