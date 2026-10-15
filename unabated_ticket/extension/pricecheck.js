// Unabated Ticket — the other books' prices at the same number (2026-10-09,
// Cal's picks): how many books beat the price you are about to take on the
// same game, period, market, side and number, and whether the price is an
// outlier that beats every other book by a mile (the 2026-10-03 Jackson St
// ML at Novig +300 against a -309 fair: almost certainly a bad price).
//
// A row a book beats gets a Skip chip (Cal, 2026-10-13: he decides, nothing
// is skipped for him). A skipped row folds into the Edges tab's "N skipped"
// line until no book beats its price any more or its game starts, whichever
// comes first; Restore brings it back at once (splitSkipped).
//
// Rules (Cal's picks 1a 2a 5a 6a):
//   - same number only, compared by price (Unabated's price as-is, so an
//     exchange's is already all-in of its fee);
//   - every book on Unabated counts, not only the ticked ones;
//   - one book under two country names (TheScore US and TheScore CA post one
//     feed) counts once when both post the same price at the number (Cal,
//     2026-10-13: "books with the same name should be matched");
//   - 2 or more books beating the price reads amber ("3 books better");
//   - outlier = the price beats the next best book by more than
//     OUTLIER_PROB_GAP in implied probability; it is tagged, never skipped;
//   - a book reporting under MIN_OTHER_BOOK_LIQUIDITY dollars resting at its
//     price (exchanges; sportsbooks report none and always count) is listed
//     greyed as thin and counts toward nothing: not better, not the outlier
//     check, not the book count (Cal, 2026-10-15: a $4 Kalshi offer must not
//     bring up Skip).
//
// Pure: no DOM, no fetch, no chrome.*, no clock (now comes in). Input is the
// scanner's feed state (feed.js shape). Loaded as a plain <script> before
// edgerows.js (exposes globalThis.UnabatedPriceCheck) and via require() in node.

(function (root) {
  "use strict";

  const inNode = typeof module !== "undefined" && module.exports;
  const feed = inNode ? require("./feed.js") : root.UnabatedFeed;

  const STATUS_ON_BOARD = 1;
  // "Multiple books" beating the price reads amber (Cal, 2026-10-09).
  const BOOKS_BETTER_TO_SKIP = 2;
  // A normal best price leads the next book by a few probability points; a
  // stale-but-real line by maybe 5. Beyond 10 points the price is more likely
  // wrong than a gift, so it reads "outlier, check".
  const OUTLIER_PROB_GAP = 0.10;
  // Dollars resting at the price (the feed's liquidity, the "liq" figure)
  // below which another book's price is too thin to act on.
  const MIN_OTHER_BOOK_LIQUIDITY = 30;

  // Does this book report too little money at its price to count.
  function isThin(entry) {
    return entry.liquidity != null && entry.liquidity < MIN_OTHER_BOOK_LIQUIDITY;
  }

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

  // A book's name without a trailing country tag: "TheScore US" and
  // "TheScore CA" are one book. Two books sharing it count once only where
  // they post the same price (Polymarket and Polymarket US price differently).
  const COUNTRY_SUFFIX = /\s*\(?\b(us|usa|ca|canada)\)?\s*$/i;
  function sameNameKey(bookName) {
    return String(bookName || "").replace(COUNTRY_SUFFIX, "").trim().toLowerCase();
  }

  // Drops each book that repeats a same-name book's price at this number,
  // naming the survivor after both ("TheScore US / CA"). yours, when given,
  // drops the same-name copy of your own price.
  function mergeSameNameBooks(entries, yours) {
    const kept = new Map();
    for (const entry of entries) {
      const key = `${sameNameKey(entry.bookName)}|${entry.price}`;
      if (yours && key === `${sameNameKey(yours.bookName)}|${yours.price}`) continue;
      const held = kept.get(key);
      if (!held) {
        kept.set(key, { ...entry });
        continue;
      }
      const suffix = entry.bookName.slice(sameNameKey(entry.bookName).length).trim();
      held.bookName = suffix ? `${held.bookName} / ${suffix}` : held.bookName;
    }
    return Array.from(kept.values());
  }

  // rung key -> Map(bookId -> {bookId, bookName, price, liquidity}), every book's best
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
      byBook.set(line.bookId, { bookId: line.bookId, bookName: state.books[line.bookId].name, price: line.price, liquidity: line.liquidity ?? null });
    }
    return index;
  }

  // How the price compares with every other book at the same number.
  //   spec  {eventId, periodTypeId, betTypeId, sideIndex, points, bookId, bookName, price}
  // Returns null when the price is unusable, else
  //   {better: [{bookName, price, liquidity}] best first, worse: [...] best first, same: n,
  //    sameAs: [...] (the n books at your price), thin: [...] best first (isThin,
  //    in none of the others), bookCount (counted books at the number, yours
  //    included), wouldSkip, outlier: {nextBest}|null}
  // "better" means a strictly lower implied probability for the same bet;
  // liquidity is null for a book that reports none.
  function comparePrice(index, spec) {
    const yourProb = impliedProb(spec.price);
    if (yourProb == null) return null;
    const atNumber = Array.from((index.get(rungKeyOf(spec)) || new Map()).values())
      .filter((entry) => entry.bookId !== spec.bookId);
    const yours = spec.bookName ? { bookName: spec.bookName, price: spec.price } : null;
    const merged = mergeSameNameBooks(atNumber, yours)
      .sort((a, b) => impliedProb(a.price) - impliedProb(b.price));
    const others = merged.filter((entry) => !isThin(entry));
    const thin = merged.filter(isThin).map((entry) => ({ bookName: entry.bookName, price: entry.price, liquidity: entry.liquidity }));
    const better = [];
    const worse = [];
    const sameAs = [];
    for (const entry of others) {
      const prob = impliedProb(entry.price);
      const shown = { bookName: entry.bookName, price: entry.price, liquidity: entry.liquidity ?? null };
      if (prob < yourProb) better.push(shown);
      else if (prob > yourProb) worse.push(shown);
      else sameAs.push(shown);
    }
    const nextBest = others.length ? others[0] : null;
    const isOutlier = better.length === 0 && nextBest != null && impliedProb(nextBest.price) - yourProb > OUTLIER_PROB_GAP;
    return {
      better,
      worse,
      same: sameAs.length,
      sameAs,
      thin,
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
    if (check.wouldSkip) return { kind: "skip", label: `${n} books better` };
    return { kind: "better", label: "1 book better" };
  }

  // The Edges header's count of outliers over the listed rows: "1 outlier to
  // check", or "" when there is none.
  function describePriceChecks(rows) {
    const outliers = rows.filter((row) => row.priceCheck && row.priceCheck.outlier).length;
    if (outliers === 0) return "";
    return `${outliers} outlier${outliers === 1 ? "" : "s"} to check`;
  }

  // One line at one book: what a Skip is keyed on (bets.duckdb::edge_skips).
  // A different number or a different book is a different line.
  function skipKeyOf(line, bookId) {
    return `${rungKeyOf(line)}@${bookId}`;
  }

  // A Skip is offered only where a book beats the price.
  function canSkip(check) {
    return Boolean(check && check.better.length > 0);
  }

  // The rows Cal skipped, apart from the rest. A skip holds only while a book
  // still beats the row's price and its game is still to start; after that the
  // row lists again by itself (the restore rule, Cal 2026-10-09/13).
  //   rows   Edges rows carrying priceCheck, book.id and eventStartMs
  //   skips  {skipKey: eventStartMs}; now  epoch ms
  // Returns {listed, skipped}, each in the rows' order.
  function splitSkipped(rows, skips, now) {
    const listed = [];
    const skipped = [];
    for (const row of rows) {
      const eventStartMs = skips[skipKeyOf(row, row.book.id)];
      const held = typeof eventStartMs === "number" && eventStartMs > now && canSkip(row.priceCheck);
      (held ? skipped : listed).push(row);
    }
    return { listed, skipped };
  }

  const api = {
    BOOKS_BETTER_TO_SKIP, OUTLIER_PROB_GAP, MIN_OTHER_BOOK_LIQUIDITY,
    impliedProb, isThin, rungKeyOf, sameNameKey, mergeSameNameBooks, buildPriceIndex, comparePrice, priceCheckTag, describePriceChecks,
    skipKeyOf, canSkip, splitSkipped,
  };

  if (inNode) {
    module.exports = api;
  } else {
    root.UnabatedPriceCheck = api;
  }
})(typeof globalThis !== "undefined" ? globalThis : this);
