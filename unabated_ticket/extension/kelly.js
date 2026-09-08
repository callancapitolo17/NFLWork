// Kelly math for the Unabated Ticket panel. Pure functions, no DOM.
//
// Mode A of the Kelly Calculator sheet: the fair price is already no-vig
// (Unabated's `bacr` at the book's points), so p_fair comes straight from it.
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

  // Below this the "edge" is floating-point noise (e.g. -110 vs fair -110
  // comes out at ~1e-13), not a bet.
  const ZERO_EDGE_EPSILON = 1e-9;

  // Full-Kelly fraction of bankroll. 0 when the bet has no edge.
  function fullKellyFraction(bookPrice, fairPrice) {
    const decimalBook = americanToDecimal(bookPrice);
    const fairProb = americanToProb(fairPrice);
    const netOdds = decimalBook - 1;
    const fraction = (fairProb * netOdds - (1 - fairProb)) / netOdds;
    return fraction > ZERO_EDGE_EPSILON ? fraction : 0;
  }

  // Edge per $1 staked: p_fair * dec_book - 1. Negative means -EV.
  function edgeFraction(bookPrice, fairPrice) {
    return americanToProb(fairPrice) * americanToDecimal(bookPrice) - 1;
  }

  function assertPositiveNumber(value, label) {
    if (typeof value !== "number" || !Number.isFinite(value) || value <= 0) {
      throw new Error(`${label}: expected a positive number, got ${value}`);
    }
  }

  // Returns { stake, fullKellyStake, fullKellyFraction, edge }. Dollars are
  // NOT rounded (user decision 2026-09-08); the panel formats them.
  function kellyStake({ bookPrice, fairPrice, bankroll, multiplier }) {
    assertPositiveNumber(bankroll, "bankroll");
    assertPositiveNumber(multiplier, "multiplier");
    const fraction = fullKellyFraction(bookPrice, fairPrice);
    const fullKellyStake = bankroll * fraction;
    return {
      stake: fullKellyStake * multiplier,
      fullKellyStake,
      fullKellyFraction: fraction,
      edge: edgeFraction(bookPrice, fairPrice),
    };
  }

  const api = {
    americanToDecimal,
    americanToProb,
    fullKellyFraction,
    edgeFraction,
    kellyStake,
  };

  if (typeof module !== "undefined" && module.exports) {
    module.exports = api;
  } else {
    root.UnabatedKelly = api;
  }
})(typeof globalThis !== "undefined" ? globalThis : this);
