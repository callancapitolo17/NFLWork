// Run: node --test unabated_ticket/tests
// server/laddersnapshot.js: Unabated's full-game fair ladders read off a
// league snapshot, and which opener / close rows the runner sends. The POST
// itself runs end to end in runner.test.js.
const test = require("node:test");
const assert = require("node:assert/strict");
const fs = require("node:fs");
const path = require("node:path");
const feed = require("../extension/feed.js");
const ladder = require("../extension/ladder.js");
const laddersnapshot = require("../server/laddersnapshot.js");

const NFL = 1;
const KICKOFF_MS = Date.parse("2026-09-13T17:00:00Z");
const OBSERVED_MS = KICKOFF_MS - 3600 * 1000;
const BEARS_PANTHERS = 125807;

function sliceState() {
  const snapshot = JSON.parse(fs.readFileSync(path.join(__dirname, "fixtures", "v2_slice.json"), "utf8"));
  return feed.mergeStates([feed.parseSnapshot(snapshot, { leagueId: NFL })]);
}

function sliceLadders(observedMs = OBSERVED_MS) {
  const state = sliceState();
  return laddersnapshot.leagueLadders({
    leagueId: NFL, linesByEvent: ladder.groupLinesByEvent(Object.values(state.lines)), events: state.events, observedMs,
  });
}

function line(fields) {
  return { periodTypeId: 1, betTypeId: 2, sideIndex: 0, points: -3.5, bacr: 110, bookId: 4, isAlt: false, isBlurred: false, ...fields };
}

test("leagueLadders: the slice's game has a moneyline, spread and total ladder with Unabated's main on both sides", () => {
  const ladders = sliceLadders();
  assert.deepEqual(ladders.map((l) => [l.eventId, l.betTypeId]), [[BEARS_PANTHERS, 1], [BEARS_PANTHERS, 2], [BEARS_PANTHERS, 3]]);
  const [moneyline, spread, total] = ladders;
  assert.deepEqual(moneyline.rungs, [[0, null, -156, 9, 9], [1, null, 156, 8, 8]]);
  assert.deepEqual([spread.mainPoints0, spread.mainFair0, spread.mainPoints1, spread.mainFair1], [-3, -103, 3, 103]);
  assert.deepEqual([total.mainPoints0, total.mainFair0, total.mainPoints1, total.mainFair1], [47.5, 108, 47.5, -108]);
  assert.ok(spread.rungs.length > 10, "the alt rungs are in");
  assert.deepEqual([spread.leagueId, spread.periodTypeId, spread.eventStartMs, spread.observedMs], [NFL, 1, KICKOFF_MS, OBSERVED_MS]);
});

test("leagueLadders: a game whose snapshot landed at or after its start has no ladder", () => {
  assert.deepEqual(sliceLadders(KICKOFF_MS), []);
  assert.deepEqual(sliceLadders(Number.NaN), []);
});

test("marketLadderOf: the most common fair wins, ties go to Unabated's own line, then the lowest book id", () => {
  const ladderOf = (lines) => laddersnapshot.marketLadderOf(lines, 2).rungs;
  assert.deepEqual(ladderOf([line({ bookId: 4, bacr: 110 }), line({ bookId: 5, bacr: 110 }), line({ bookId: 6, bacr: 112 })]),
    [[0, -3.5, 110, 3, 2]]);
  assert.deepEqual(ladderOf([line({ bookId: 4, bacr: 110 }), line({ bookId: feed.UNABATED_LINE_BOOK_ID, bacr: 112 })]),
    [[0, -3.5, 112, 2, 1]]);
  assert.deepEqual(ladderOf([line({ bookId: 9, bacr: 110 }), line({ bookId: 4, bacr: 112 })]), [[0, -3.5, 112, 2, 1]]);
});

test("marketLadderOf: blurred lines, the ±999900 clamp, other periods and a draw side are not fairs", () => {
  const lines = [
    line({ isBlurred: true }), line({ bacr: 999900 }), line({ bacr: null }), line({ periodTypeId: 2 }), line({ sideIndex: 2 }),
    line({ points: null }),
  ];
  assert.equal(laddersnapshot.marketLadderOf(lines, 2), null);
  const kept = laddersnapshot.marketLadderOf([...lines, line({ sideIndex: 1, points: 3.5, bacr: -110 }), line({ points: -7.5, bacr: 250, isAlt: true })], 2);
  assert.deepEqual(kept.rungs, [[0, -7.5, 250, 1, 1], [1, 3.5, -110, 1, 1]]);
  // No Unabated line in the market: no main.
  assert.deepEqual([kept.mainPoints0, kept.mainFair0, kept.mainPoints1, kept.mainFair1], [null, null, null, null]);
});

test("unsentRows: an opener once, a close on a change or every CLOSE_CONFIRM_MS, nothing for a started game", () => {
  const sent = new Map();
  const [spread] = sliceLadders().filter((l) => l.betTypeId === 2);
  const first = laddersnapshot.unsentRows([spread], sent, OBSERVED_MS);
  assert.deepEqual(first.map((row) => row.reading), ["opener", "close"]);
  assert.equal(first[0].observedAt, new Date(OBSERVED_MS).toISOString());
  assert.equal(first[0].eventStart, new Date(KICKOFF_MS).toISOString());
  assert.equal("signature" in laddersnapshot.wireRow(first[1]), false);
  laddersnapshot.markSent(first, sent);
  assert.deepEqual(laddersnapshot.unsentRows([spread], sent, OBSERVED_MS), []);
  // A newer snapshot with the same ladder: not re-sent until the confirmation is due.
  const later = { ...spread, observedMs: OBSERVED_MS + laddersnapshot.CLOSE_CONFIRM_MS - 1 };
  assert.deepEqual(laddersnapshot.unsentRows([later], sent, later.observedMs), []);
  const due = { ...spread, observedMs: OBSERVED_MS + laddersnapshot.CLOSE_CONFIRM_MS };
  assert.deepEqual(laddersnapshot.unsentRows([due], sent, due.observedMs).map((row) => row.reading), ["close"]);
  // A changed ladder: the close at once.
  const moved = { ...spread, observedMs: OBSERVED_MS + 60 * 1000, signature: `${spread.signature}x` };
  assert.deepEqual(laddersnapshot.unsentRows([moved], sent, moved.observedMs).map((row) => row.reading), ["close"]);
  // Kicked off since the snapshot: nothing, and its record goes.
  assert.deepEqual(laddersnapshot.unsentRows([moved], sent, KICKOFF_MS), []);
  laddersnapshot.forgetStarted(sent, KICKOFF_MS);
  assert.equal(sent.size, 0);
});
