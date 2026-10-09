// Unabated Ticket — the other books' prices at the same number (2026-10-09,
// Cal's picks): how many books beat the price you are about to take on the
// same game, period, market, side and number, and whether the price is an
// outlier that beats every other book by a mile (the 2026-10-03 Jackson St
// ML at Novig +300 against a -309 fair: almost certainly a bad price).
//
// Test mode for now: the panel and the phone only TAG rows ("3 books better ·
// would skip"); nothing is hidden until Cal has watched it and says so.
//
// Rules (Cal's picks 1a 2a 5a 6a):
//   - same number only, compared by price (Unabated's price as-is, so an
//     exchange's is already all-in of its fee);
//   - every book on Unabated counts, not only the ticked ones;
//   - would skip = 2 or more books beat the price;
//   - outlier = the price beats the next best book by more than
//     OUTLIER_PROB_GAP in implied probability; it is tagged, never skipped.
//
// Pure: no DOM, no fetch, no chrome.*, no clock (now comes in). Input is the
// scanner's feed state (feed.js shape). Loaded as a plain <script> before
// edgerows.js (exposes globalThis.UnabatedPriceCheck) and via require() in node.

(function (root) {
  "use strict";

  const inNode = typeof module !== "undefined" && module.exports;
  const feed = inNode ? require("./feed.js") : root.UnabatedFeed;

  const STATUS_ON_BOARD = 1;
  // "Multiple books" beating the price is Cal's skip rule (2026-10-09).
  const BOOKS_BETTER_TO_SKIP = 2;
  // A normal best price leads the next book by a few probability points; a
  // stale-but-real line by maybe 5. Beyond 10 points the price is more likely
  // wrong than a gift, so it reads "outlier, check".
  const OUTLIER_PROB_GAP = 0.10;

  // A book's price as a probability, or null when it is no American price.
  function impliedProb(americanOdds) {
    if (typeof americanOdds !== "number" || !Number.isFinite(americanOdds) || Math.abs(americanOdds) < 100) return null;
    return americanOdds > 0 ? 100 / (americanOdds + 100) : -americanOdds / (-americanOdds + 100);
  }

  const BET_TYPE_MONEYLINE = 1;

  // One number on one side of one market of one game, whoever posts it. A
  // moneyline has no number, whatever points field a book sends with it.
  function rungKeyOf({ eventId, periodTypeId, betTypeId, sideIndex, points }) {
    return `${eventId}:${periodTypeId}:${betTypeId}:${sideIndex}|${betTypeId === BET_TYPE_MONEYLINE ? "ml" : points}`;
  }

  // Can this line stand as another book's price: a live book's line on the
  // board, priced, not blurred (a blurred price is Unabated's placeholder),
  // not Unabated's own line, and changed within maxLineAgeMs when set (a dead
  // feed's week-old price must not read as "better").
  function countsAsPrice(line, state, now, maxLineAgeMs) {
    if (line.bookId === feed.UNABATED_LINE_BOOK_ID || line.statusId !== STATUS_ON_BOARD || line.isBlurred) return false;
    if (impliedProb(line.price) == null) return false;
    const book = state.books[line.bookId];
    if (!book || !book.isLive) return false;
    if (maxLineAgeMs == null) return true;
    const changedMs = feed.lineChangedMs(line);
    return changedMs != null && now - changedMs <= maxLineAgeMs;
  }

  // rung key -> Map(bookId -> {bookId, bookName, price}), every book's best
  // price at each number (a book can sit on a number twice: main and alt).
  // Build once per feed update.
  //   options  {now, maxLineAgeMs (null = no age gate)}
  function buildPriceIndex(state, options) {
    const index = new Map();
    if (!state) return index;
    const now = options.now;
    const maxLineAgeMs = typeof options.maxLineAgeMs === "number" && options.maxLineAgeMs > 0 ? options.maxLineAgeMs : null;
    for (const line of Object.values(state.lines)) {
      if (!countsAsPrice(line, state, now, maxLineAgeMs)) continue;
      const key = rungKeyOf(line);
      if (!index.has(key)) index.set(key, new Map());
      const byBook = index.get(key);
      const held = byBook.get(line.bookId);
      if (held && impliedProb(held.price) <= impliedProb(line.price)) continue;
      byBook.set(line.bookId, { bookId: line.bookId, bookName: state.books[line.bookId].name, price: line.price });
    }
    return index;
  }

  // How the price compares with every other book at the same number.
  //   spec  {eventId, periodTypeId, betTypeId, sideIndex, points, bookId, price}
  // Returns null when the price is unusable, else
  //   {better: [{bookName, price}] best first, worse: [...] best first, same: n,
  //    bookCount (books at the number, yours included), wouldSkip, outlier: {nextBest}|null}
  // "better" means a strictly lower implied probability for the same bet.
  function comparePrice(index, spec) {
    const yourProb = impliedProb(spec.price);
    if (yourProb == null) return null;
    const others = Array.from((index.get(rungKeyOf(spec)) || new Map()).values())
      .filter((entry) => entry.bookId !== spec.bookId)
      .sort((a, b) => impliedProb(a.price) - impliedProb(b.price));
    const better = [];
    const worse = [];
    let same = 0;
    for (const entry of others) {
      const prob = impliedProb(entry.price);
      const shown = { bookName: entry.bookName, price: entry.price };
      if (prob < yourProb) better.push(shown);
      else if (prob > yourProb) worse.push(shown);
      else same += 1;
    }
    const nextBest = others.length ? others[0] : null;
    const isOutlier = better.length === 0 && nextBest != null && impliedProb(nextBest.price) - yourProb > OUTLIER_PROB_GAP;
    return {
      better,
      worse,
      same,
      bookCount: others.length + 1,
      wouldSkip: better.length >= BOOKS_BETTER_TO_SKIP,
      outlier: isOutlier ? { nextBest: { bookName: nextBest.bookName, price: nextBest.price } } : null,
    };
  }

  // The tag a row or ticket shows: {kind, label}, or null when no other book
  // posts the number (nothing to compare). kind: best | better | skip | outlier.
  function priceCheckTag(check) {
    if (!check || check.bookCount < 2) return null;
    const fmt = (price) => (price > 0 ? `+${price}` : `${price}`);
    if (check.outlier) return { kind: "outlier", label: `outlier, check · next best ${fmt(check.outlier.nextBest.price)}` };
    const n = check.better.length;
    if (n === 0) return { kind: "best", label: check.same > 0 ? `best, tied · ${check.bookCount} books` : `best of ${check.bookCount} books` };
    if (check.wouldSkip) return { kind: "skip", label: `${n} books better · would skip` };
    return { kind: "better", label: "1 book better" };
  }

  // The Edges header's count over the listed rows: "3 would skip: 2+ books
  // better · 1 outlier (test mode)", or "" when there is neither.
  function describePriceChecks(rows) {
    let skip = 0;
    let outliers = 0;
    for (const row of rows) {
      if (!row.priceCheck) continue;
      if (row.priceCheck.wouldSkip) skip += 1;
      if (row.priceCheck.outlier) outliers += 1;
    }
    if (skip === 0 && outliers === 0) return "";
    const parts = [];
    if (skip) parts.push(`${skip} would skip: 2+ books better`);
    if (outliers) parts.push(`${outliers} outlier${outliers === 1 ? "" : "s"}`);
    return `${parts.join(" · ")} (test mode, nothing hidden)`;
  }

  const api = {
    BOOKS_BETTER_TO_SKIP, OUTLIER_PROB_GAP,
    impliedProb, rungKeyOf, buildPriceIndex, comparePrice, priceCheckTag, describePriceChecks,
  };

  if (inNode) {
    module.exports = api;
  } else {
    root.UnabatedPriceCheck = api;
  }
})(typeof globalThis !== "undefined" ? globalThis : this);
