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

  // Exchanges (Kalshi, Novig) publish a probability, and a contract there costs
  // that probability in whole cents and pays $1. Unabated marks such lines with
  // sourceFormat 4; every other book takes dollars, so this is the one rule
  // that says a line trades in contracts.
  const EXCHANGE_PROBABILITY_FORMAT = 4;

  function isContractMarket({ sourceFormat, sourcePrice }) {
    return sourceFormat === EXCHANGE_PROBABILITY_FORMAT && sourcePrice > 0 && sourcePrice < 1;
  }

  // Exchange order books run 1¢..99¢; a probability outside that has no
  // contract to buy, so it is "not priced in contracts", not an error.
  const MIN_CONTRACT_PROB = 0.01;
  const MAX_CONTRACT_PROB = 0.99;
  // Float noise (261.69 / 0.232012 = 1127.999...) must not drop a contract
  // that the stake covers exactly.
  const FLOOR_EPSILON = 1e-9;

  // The order that spends `stake` on an exchange line: how many contracts at
  // Unabated's price, and what they cost. Null for a line that is not priced
  // in contracts. The price is Unabated's exact number for the book, taken as
  // the all-in cost of one contract (for Kalshi that number already carries
  // Kalshi's fee: a 22¢ ask shows as 23.2¢ = +331; user decision 2026-09-22),
  // so contracts = floor(stake / price), never rounded up past Kelly, and the
  // cost reconciles to the exchange's own Cost line. The leftover is under
  // one contract.
  function contractOrder({ stake, bookPrice, sourceFormat, sourcePrice }) {
    if (!isContractMarket({ sourceFormat, sourcePrice })) return null;
    if (typeof stake !== "number" || !Number.isFinite(stake) || stake < 0) {
      throw new Error(`contractOrder: expected a non-negative stake, got ${stake}`);
    }
    const priceProb = bookProbOf({ bookPrice, sourceFormat, sourcePrice });
    if (priceProb < MIN_CONTRACT_PROB || priceProb > MAX_CONTRACT_PROB) return null;
    const contracts = Math.floor(stake / priceProb + FLOOR_EPSILON);
    const costDollars = Math.round(contracts * priceProb * 100) / 100;
    return { contracts, priceCents: priceProb * 100, costDollars };
  }

  const api = { americanToDecimal, americanToProb, bookProbOf, kellyStakeFromEdge, isContractMarket, contractOrder };

  if (typeof module !== "undefined" && module.exports) {
    module.exports = api;
  } else {
    root.UnabatedKelly = api;
  }
})(typeof globalThis !== "undefined" ? globalThis : this);
