// Unabated Ticket — the exchanges' own live prices over Unabated's copy of
// them (Cal, 2026-11-04: "build it"). Unabated's Novig, Kalshi and
// Polymarket US prices lag the real books, sometimes for hours: on
// 2026-10-31 the Novig bot found 973 lines priced worse than Unabated showed
// and 575 with nothing resting at all, against 368 it bet. The server runner
// reads each listed line's live book (runner.js livePricePass) and serves the
// reads at GET /live_prices.json; this module lays them over the board, so
// the Edges list, its stakes and "best of N" on the phone and the panel price
// the exchange's real offer:
//   - a fresh read with an offer: the line's price is the live all-in price
//     (Kalshi's fee and Polymarket US's included; Novig charges none
//     pregame), its liquidity the dollars resting at that price, and its
//     edge re-measured against Unabated's fair (fairProbOf);
//   - a fresh read with nothing resting: the line is gone (null);
//   - no fresh read (the runner has not reached it, the read failed, or the
//     number moved since): Unabated's line as it is, tagged "unchecked".
// Lines at every other book pass untouched (liveCheck null).
//
// The fair is the automated taker's (taker.js decideRow): Unabated's edge on
// the cheaper of its two shown prices, so the fair is never one the overlay
// gave itself. The taker itself does not use this module (Cal, 2026-11-02:
// the bot stays as it is).
//
// Pure: no DOM, no fetch, no chrome.*, no clock (now comes in). Loaded as a
// plain <script> before edgerows.js (exposes globalThis.UnabatedLivePrices)
// and via require() in node.

(function (root) {
  "use strict";

  // A read older than this is not laid over the board: the runner reads every
  // listed line about every 30 s, so two minutes without one means the line
  // left its list, the read keeps failing, or the runner is down.
  const LIVE_READ_MAX_AGE_MS = 2 * 60 * 1000;
  const SOURCE_FORMAT_PROBABILITY = 4;
  const LIVE = "live";
  const UNCHECKED = "unchecked";

  // An American price as a probability, or null when it is no price.
  function probOfAmerican(americanOdds) {
    if (typeof americanOdds !== "number" || !Number.isFinite(americanOdds) || Math.abs(americanOdds) < 100) return null;
    return americanOdds > 0 ? 100 / (americanOdds + 100) : -americanOdds / (-americanOdds + 100);
  }

  // A probability as a whole American price (0.25 -> +300, 0.43 -> +133), as
  // Unabated shows every price; the exact price stays in sourcePrice, which is
  // what the stake and the cents are sized on (kelly.bookProbOf).
  function americanOfProb(prob) {
    return Math.round(prob < 0.5 ? (100 * (1 - prob)) / prob : (-100 * prob) / (1 - prob));
  }

  // Unabated's fair win probability for the line: its edge (ge) on the cheaper
  // of the two prices it shows, the sourcePrice (an exchange's unrounded
  // all-in, when sent) and the American. Null when the line has no edge or price.
  function fairProbOf(line) {
    if (typeof line.ge !== "number" || !Number.isFinite(line.ge)) return null;
    const american = probOfAmerican(line.price);
    if (american == null) return null;
    const sent = line.sourceFormat === SOURCE_FORMAT_PROBABILITY && typeof line.sourcePrice === "number"
      && line.sourcePrice > 0 && line.sourcePrice < 1 ? line.sourcePrice : american;
    return (1 + line.ge) * Math.min(sent, american);
  }

  function isUsableRead(read) {
    if (!read || typeof read !== "object" || typeof read.readAtMs !== "number") return false;
    if (read.empty === true) return true;
    return typeof read.allIn === "number" && read.allIn > 0 && read.allIn < 1
      && typeof read.depthDollars === "number" && read.depthDollars >= 0;
  }

  // The line at the live price, or null when nothing rests. A line with no
  // fair (no edge) keeps Unabated's edge field as null.
  //   read  {allIn, depthDollars, readAtMs} | {empty: true, readAtMs}
  function liveLineOf(line, read) {
    if (read.empty === true) return null;
    const fair = fairProbOf(line);
    const ge = fair == null ? null : Math.round((fair / read.allIn - 1) * 1e4) / 1e4;
    return {
      ...line,
      price: americanOfProb(read.allIn),
      sourceFormat: SOURCE_FORMAT_PROBABILITY,
      sourcePrice: read.allIn,
      ge,
      liquidity: Math.round(read.depthDollars * 100) / 100,
      liveReadAtMs: read.readAtMs,
    };
  }

  // The view the Edges list reads each board line through, or null for none
  // (no live prices served yet: every line as Unabated has it, untagged).
  //   livePrices  the /live_prices.json body: {bookIds: [ids the runner reads],
  //               reads: {lineKey: {points, allIn, depthDollars, readAtMs} | {points, empty: true, readAtMs}}}
  //   now         ms epoch on the reads' clock
  // Returns (line) -> {line, liveCheck: "live" | "unchecked" | null} | null (the line is gone),
  // with `readKeys`: the keys of the only lines whose price or presence it can change.
  // A read whose number is not the line's (a main line moved) is no read.
  function createLineView(livePrices, now) {
    if (!livePrices || !Array.isArray(livePrices.bookIds) || !livePrices.reads || typeof livePrices.reads !== "object") return null;
    const covered = new Set(livePrices.bookIds);
    const reads = livePrices.reads;
    const view = (line) => {
      if (!covered.has(line.bookId)) return { line, liveCheck: null };
      const read = Object.prototype.hasOwnProperty.call(reads, line.key) ? reads[line.key] : null;
      if (!isUsableRead(read) || now - read.readAtMs > LIVE_READ_MAX_AGE_MS || (read.points ?? null) !== (line.points ?? null)) {
        return { line, liveCheck: UNCHECKED };
      }
      const live = liveLineOf(line, read);
      return live ? { line: live, liveCheck: LIVE } : null;
    };
    view.readKeys = Object.keys(reads);
    return view;
  }

  const api = { LIVE_READ_MAX_AGE_MS, LIVE, UNCHECKED, fairProbOf, liveLineOf, createLineView };

  if (typeof module !== "undefined" && module.exports) {
    module.exports = api;
  } else {
    root.UnabatedLivePrices = api;
  }
})(typeof globalThis !== "undefined" ? globalThis : this);
