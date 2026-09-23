// Why an edge grew (issue #132): since the previous observation of a line,
// what did Unabated's fair (`bacr`) do for the row's side, and what did the
// book's price do? The fair decides the tag (spec addition, 2026-09-22):
//   fair moved toward the side, whatever the price did  -> fair moved to you (green): the sharps agree, the price is a bonus
//   fair unchanged, price got better                    -> book moved away (amber): the book is shading against the side
//                                                          and the fair has not answered yet — Unabated's fair is ~1-2 min
//                                                          behind (v2 rebuild ~1/min + ~40 s ingest lag, #126); wait a snapshot
//   fair moved against the side, whatever the price did -> fair moved against you (red): the market is moving against the
//                                                          side and the book is ahead of the fair; the edge can still read
//                                                          bigger when the price improved by more than the fair fell —
//                                                          adding on it is adverse selection
//   nothing moved                                       -> no tag
// The tag informs only; the stake stays what condkelly.js computes.
//
// Pure functions, no DOM: a per-line history the scanner keeps in memory
// (`observe`), and `edgeMove` reading a tag off it. Loaded as a plain
// <script> in panel.html (globalThis.UnabatedEdgeMove) and via require() in
// tests/edgemove.test.js.

(function (root) {
  "use strict";

  const kelly = typeof module !== "undefined" && module.exports ? require("./kelly.js") : root.UnabatedKelly;

  // How far back a move still earns a tag: the fair lag plus two snapshots
  // on the slowest refresh tier (5 min, scanner.js REFRESH_TIERS), so "wait a
  // snapshot before adding" can be followed while the tag is still up.
  const EDGE_MOVE_WINDOW_MS = 10 * 60 * 1000;
  // Half a point of probability. `bacr` is a whole American price, so a
  // one-point fair change near even money (-110 -> -111, 0.23 pts) is
  // rounding, not a move.
  const MOVE_THRESHOLD = 0.005;
  // Novig prices in half cents: 0.565 - 0.56 is 0.00499999... in floats.
  const FLOAT_EPSILON = 1e-9;
  const SOURCES = ["snapshot", "stream"];
  const MOVE_LABELS = Object.freeze({ fair_to_you: "fair moved to you", book_away: "book moved away", fair_against: "fair moved against you" });
  const NO_MOVE = Object.freeze({ kind: "none", fairDelta: null, priceDelta: null, sinceMs: null, source: null, from: null, to: null });

  function entryOf(line, at, source) {
    return {
      at,
      source,
      points: line.points ?? null,
      price: line.price,
      sourceFormat: line.sourceFormat ?? 1,
      sourcePrice: line.sourcePrice ?? null,
      bacr: line.bacr ?? null,
    };
  }

  function sameValues(a, b) {
    return a.price === b.price && a.sourcePrice === b.sourcePrice && a.bacr === b.bacr;
  }

  // Index of the newest entry at or before `boundary`; 0 when none is (the
  // line was first seen inside the window and compares to its first sighting).
  function referenceIndex(entries, boundary) {
    let index = 0;
    for (let i = 0; i < entries.length; i += 1) if (entries[i].at <= boundary) index = i;
    return index;
  }

  // Entries older than the window go, except the newest of them: it is the
  // baseline the window compares against.
  function prune(entries, at) {
    const keepFrom = referenceIndex(entries, at - EDGE_MOVE_WINDOW_MS);
    return keepFrom === 0 ? entries : entries.slice(keepFrom);
  }

  // Record one observation of `line` in `history` (key -> entries, newest
  // last). An entry is added only when the price, the exchange source price
  // or the fair differ from the last one. A line whose NUMBER moved starts
  // afresh: a price at 48.5 is not comparable to one at 47.5, so it reads as
  // first seen rather than as a move. `source` says what observed it — the
  // anonymous changes stream misses most exchange moves and never carries an
  // alt rung, so for those every entry is a snapshot, up to one refresh
  // interval after the book moved.
  function observe(history, line, { at, source }) {
    if (!history || typeof history !== "object") throw new Error(`observe: expected a history object, got ${history}`);
    if (!line || typeof line.key !== "string") throw new Error(`observe: expected a line with a key, got ${line && line.key}`);
    if (typeof line.price !== "number" || !Number.isFinite(line.price)) throw new Error(`observe: expected a numeric price on ${line.key}, got ${line.price}`);
    if (typeof at !== "number" || !Number.isFinite(at)) throw new Error(`observe: expected a numeric time, got ${at}`);
    if (!SOURCES.includes(source)) throw new Error(`observe: expected source ${SOURCES.join("|")}, got ${source}`);
    const entry = entryOf(line, at, source);
    const entries = history[line.key];
    if (!entries || entries[entries.length - 1].points !== entry.points) {
      history[line.key] = [entry];
      return;
    }
    if (!sameValues(entries[entries.length - 1], entry)) entries.push(entry);
    history[line.key] = prune(entries, at);
  }

  function forget(history, key) {
    delete history[key];
  }

  // A fair the panel cannot price (missing, or a price no American odds
  // describe) is "unknown", never a crash in the render path.
  function fairProbOf(entry) {
    if (entry.bacr == null) return null;
    try {
      return kelly.americanToProb(entry.bacr);
    } catch (_error) {
      return null;
    }
  }

  // The exchange's exact source price when there is one, so a Kalshi cent or
  // a Novig half cent is measured as itself, not as its whole-American rounding.
  function bookProbOf(entry) {
    try {
      return kelly.bookProbOf({ bookPrice: entry.price, sourceFormat: entry.sourceFormat, sourcePrice: entry.sourcePrice });
    } catch (_error) {
      return null;
    }
  }

  // The tag for one line's entries at `now`. Deltas are in probability and
  // signed so that positive is edge-increasing for the row's side: fairDelta
  // = the fair's chance of the side now minus then; priceDelta = the book's
  // implied chance then minus now (the book lengthened the side's odds).
  // The reference is the newest entry at or before now - window (else the
  // first sighting); the reference being the current entry, or one entry in
  // all, is `none`. The fair decides (see the header): at least the
  // threshold up is `fair_to_you`, at least the threshold down is
  // `fair_against`, both whatever the price did; a fair inside the threshold
  // (or unknown at either end) is `book_away` when the price got better by
  // the threshold, else `none`. A book that only shortened never earns a tag.
  // `sinceMs` and `source` are the first entry after the reference: when the
  // move was first observed, and by what.
  function edgeMove(entries, now) {
    if (entries == null) return NO_MOVE;
    if (!Array.isArray(entries)) throw new Error(`edgeMove: expected an array of entries, got ${typeof entries}`);
    if (typeof now !== "number" || !Number.isFinite(now)) throw new Error(`edgeMove: expected a numeric time, got ${now}`);
    if (entries.length < 2) return NO_MOVE;
    const reference = referenceIndex(entries, now - EDGE_MOVE_WINDOW_MS);
    if (reference === entries.length - 1) return NO_MOVE;
    const from = entries[reference];
    const to = entries[entries.length - 1];
    const firstChange = entries[reference + 1];
    const fairFrom = fairProbOf(from);
    const fairTo = fairProbOf(to);
    const fairDelta = fairFrom == null || fairTo == null ? null : fairTo - fairFrom;
    const bookFrom = bookProbOf(from);
    const bookTo = bookProbOf(to);
    const priceDelta = bookFrom == null || bookTo == null ? null : bookFrom - bookTo;
    const fairUp = fairDelta != null && fairDelta >= MOVE_THRESHOLD - FLOAT_EPSILON;
    const fairDown = fairDelta != null && fairDelta <= -(MOVE_THRESHOLD - FLOAT_EPSILON);
    const priceUp = priceDelta != null && priceDelta >= MOVE_THRESHOLD - FLOAT_EPSILON;
    let kind = "none";
    if (fairUp) kind = "fair_to_you";
    else if (fairDown) kind = "fair_against";
    else if (priceUp) kind = "book_away";
    return { kind, fairDelta, priceDelta, sinceMs: now - firstChange.at, source: firstChange.source, from, to };
  }

  const api = { EDGE_MOVE_WINDOW_MS, MOVE_THRESHOLD, MOVE_LABELS, observe, forget, edgeMove };

  if (typeof module !== "undefined" && module.exports) {
    module.exports = api;
  } else {
    root.UnabatedEdgeMove = api;
  }
})(typeof globalThis !== "undefined" ? globalThis : this);
