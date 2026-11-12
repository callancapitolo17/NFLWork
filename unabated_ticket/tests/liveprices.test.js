// Run: node --test unabated_ticket/tests
// liveprices.js: the exchanges' own live prices laid over Unabated's copy of
// them on the Edges list (Cal, 2026-11-04).
const test = require("node:test");
const assert = require("node:assert/strict");
const feed = require("../extension/feed.js");
const pricecheck = require("../extension/pricecheck.js");
const livePrices = require("../extension/liveprices.js");

const NOW = Date.parse("2026-10-31T15:00:00Z");
const HOUR = 3600 * 1000;
const NOVIG = 89;
const KALSHI = 105;
const BET105 = 7;

// One game, a total at 47.5; each line [bookId, sideIndex, price, extra].
function stateWith(lines) {
  const state = feed.emptyState();
  state.books = {
    [NOVIG]: { id: NOVIG, name: "Novig", isLive: true, hasLiquidity: true },
    [KALSHI]: { id: KALSHI, name: "Kalshi", isLive: true, hasLiquidity: true },
    [BET105]: { id: BET105, name: "Bet105", isLive: true },
  };
  state.events[100] = { eventId: 100, leagueId: 2, eventName: "Wagner @ Merrimack", eventStart: NOW + HOUR };
  for (const [bookId, sideIndex, price, extra] of lines) {
    const key = `m3:ms${bookId}:si${sideIndex}`;
    state.lines[key] = {
      key, isAlt: false, leagueId: 2, periodTypeId: 1, betTypeId: 3, eventId: 100, bookId, sideIndex, sideKey: `si${sideIndex}`,
      points: 47.5, price, statusId: 1, isBlurred: false, modifiedOn: new Date(NOW - HOUR).toISOString(), ...(extra || {}),
    };
  }
  return state;
}

// Cal's 2026-10-31 refusal: Unabated showed Novig Under 47.5 at -115 (53.5c) with
// $209 and a 11.87% edge; Novig's real book had nothing, then 60.5c.
const MERRIMACK_UNDER = { ge: 0.1187, sourceFormat: 4, sourcePrice: 0.535, liquidity: 209.4 };

function servedReads(reads) {
  return { bookIds: [NOVIG, KALSHI], reads };
}

test("the fair is Unabated's edge on the cheaper of its two prices, as the taker reads it", () => {
  const line = { price: -115, ge: 0.1187, sourceFormat: 4, sourcePrice: 0.535 };
  // -115 is 53.488%, under the 53.5c sourcePrice: the fair is 1.1187 x 0.53488 = 59.84% (the Ticket's -149).
  assert.ok(Math.abs(livePrices.fairProbOf(line) - 1.1187 * (115 / 215)) < 1e-12);
  assert.equal(livePrices.fairProbOf({ price: -115, ge: null }), null);
  assert.equal(livePrices.fairProbOf({ price: 50, ge: 0.02 }), null);
});

test("a live offer reprices the line, its edge and its liquidity; nothing resting drops it", () => {
  const state = stateWith([[NOVIG, 1, -115, MERRIMACK_UNDER]]);
  const line = state.lines[`m3:ms${NOVIG}:si1`];
  const live = livePrices.liveLineOf(line, { allIn: 0.605, depthDollars: 111.85, readAtMs: NOW - 5000 });
  // -153.16 shown whole, as Unabated shows every price; the exact price is sourcePrice.
  assert.equal(live.price, -153);
  assert.equal(live.sourceFormat, 4);
  assert.equal(live.sourcePrice, 0.605);
  // 59.84% fair against a 60.5c offer: a 1.1% loss, not an 11.87% edge.
  assert.equal(live.ge, Math.round((1.1187 * (115 / 215) / 0.605 - 1) * 1e4) / 1e4);
  assert.ok(live.ge < 0);
  assert.equal(live.liquidity, 111.85);
  assert.equal(live.liveReadAtMs, NOW - 5000);
  assert.equal(livePrices.liveLineOf(line, { empty: true, readAtMs: NOW }), null);
});

test("the line view: other books untouched, a fresh read applied, a stale, missing or moved read unchecked", () => {
  const state = stateWith([[NOVIG, 1, -115, MERRIMACK_UNDER], [KALSHI, 1, -110], [BET105, 1, -157]]);
  const view = livePrices.createLineView(servedReads({
    [`m3:ms${NOVIG}:si1`]: { venue: "novig", points: 47.5, allIn: 0.52, depthDollars: 300, readAtMs: NOW - 30 * 1000 },
  }), NOW);
  assert.deepEqual(view.readKeys, [`m3:ms${NOVIG}:si1`]);
  const novig = view(state.lines[`m3:ms${NOVIG}:si1`]);
  assert.equal(novig.liveCheck, "live");
  assert.equal(novig.line.sourcePrice, 0.52);
  assert.deepEqual(view(state.lines[`m3:ms${BET105}:si1`]), { line: state.lines[`m3:ms${BET105}:si1`], liveCheck: null });
  assert.equal(view(state.lines[`m3:ms${KALSHI}:si1`]).liveCheck, "unchecked");

  const stale = livePrices.createLineView(servedReads({
    [`m3:ms${NOVIG}:si1`]: { points: 47.5, allIn: 0.52, depthDollars: 300, readAtMs: NOW - livePrices.LIVE_READ_MAX_AGE_MS - 1 },
  }), NOW);
  assert.deepEqual(stale(state.lines[`m3:ms${NOVIG}:si1`]), { line: state.lines[`m3:ms${NOVIG}:si1`], liveCheck: "unchecked" });
  // A main line that moved to another number since the read: the read is of another bet.
  const moved = livePrices.createLineView(servedReads({
    [`m3:ms${NOVIG}:si1`]: { points: 48.5, allIn: 0.52, depthDollars: 300, readAtMs: NOW },
  }), NOW);
  assert.equal(moved(state.lines[`m3:ms${NOVIG}:si1`]).liveCheck, "unchecked");
  const empty = livePrices.createLineView(servedReads({ [`m3:ms${NOVIG}:si1`]: { points: 47.5, empty: true, readAtMs: NOW } }), NOW);
  assert.equal(empty(state.lines[`m3:ms${NOVIG}:si1`]), null);
  assert.equal(livePrices.createLineView(null, NOW), null);
});

test("selectEdges through the view: the ghost is gone, a real offer lists at its own price and depth, unread lines are tagged", () => {
  const state = stateWith([
    [NOVIG, 1, -115, MERRIMACK_UNDER],
    [KALSHI, 1, -110, { ge: 0.05, liquidity: 5000 }],
    [NOVIG, 0, -105, { ge: 0.04, liquidity: 1000 }],
    [BET105, 1, -157, { ge: 0.02 }],
  ]);
  const options = { now: NOW, minEdge: 0.01, bookIds: new Set([NOVIG, KALSHI, BET105]), betTypes: new Set([3]), minLiquidityToWin: 100 };
  const before = feed.selectEdges(state, options).map((row) => [row.key, row.liveCheck]);
  assert.deepEqual(before, [[`m3:ms${NOVIG}:si1`, null], [`m3:ms${KALSHI}:si1`, null], [`m3:ms${NOVIG}:si0`, null], [`m3:ms${BET105}:si1`, null]]);

  const lineView = livePrices.createLineView(servedReads({
    [`m3:ms${NOVIG}:si1`]: { points: 47.5, empty: true, readAtMs: NOW - 1000 },
    // Kalshi -110 (52.38%) at a 5% edge is a 55.0% fair; the live 52c offer is a 5.77% edge, $40 resting.
    [`m3:ms${KALSHI}:si1`]: { points: 47.5, allIn: 0.52, depthDollars: 40, readAtMs: NOW - 1000 },
  }), NOW);
  const rows = feed.selectEdges(state, { ...options, lineView });
  // Kalshi's $40 wins $36.9 at 52c: under the $100 Min liq to win, so it goes too.
  assert.deepEqual(rows.map((row) => [row.key, row.liveCheck]), [[`m3:ms${NOVIG}:si0`, "unchecked"], [`m3:ms${BET105}:si1`, null]]);

  const roomier = feed.selectEdges(state, { ...options, minLiquidityToWin: 20, lineView });
  const kalshi = roomier.find((row) => row.key === `m3:ms${KALSHI}:si1`);
  assert.equal(kalshi.liveCheck, "live");
  assert.equal(kalshi.edgePct, Math.round(((1.05 * (110 / 210)) / 0.52 - 1) * 1e4) / 1e2);
  assert.equal(kalshi.liquidity, 40);
  assert.equal(kalshi.liveReadAtMs, NOW - 1000);
});

test("a ghost exchange price no longer beats another book in best of N", () => {
  // Novig's copy at +100 beats Kalshi's -110; Novig's live book has nothing.
  const state = stateWith([[NOVIG, 1, 100, { liquidity: 500 }], [KALSHI, 1, -110, { liquidity: 500 }]]);
  const spec = { eventId: 100, periodTypeId: 1, betTypeId: 3, sideIndex: 1, points: 47.5, bookId: KALSHI, bookName: "Kalshi", price: -110 };
  const raw = pricecheck.comparePrice(pricecheck.buildPriceIndex(state, { now: NOW, maxLineAgeMs: null }), spec);
  assert.equal(raw.better.length, 1);
  const lineView = livePrices.createLineView(servedReads({ [`m3:ms${NOVIG}:si1`]: { points: 47.5, empty: true, readAtMs: NOW } }), NOW);
  const live = pricecheck.comparePrice(pricecheck.buildPriceIndex(state, { now: NOW, maxLineAgeMs: null, lineView }), spec);
  assert.equal(live.better.length, 0);
});

test("a live price is a whole American number (Cal, 2026-11-05: +132.558... on the phone Ticket)", () => {
  const state = stateWith([[NOVIG, 1, -115, MERRIMACK_UNDER]]);
  const line = state.lines[`m3:ms${NOVIG}:si1`];
  const live = livePrices.liveLineOf(line, { allIn: 0.43, depthDollars: 98, readAtMs: NOW });
  assert.equal(live.price, 133);
  assert.equal(live.sourcePrice, 0.43);
});
