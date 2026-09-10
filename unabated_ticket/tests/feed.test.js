// Run: node --test unabated_ticket/tests
// Fixtures are real slices of both feeds captured 2026-09-10 (NFL event
// 125807, Chicago Bears @ Carolina Panthers, kickoff 2026-09-13T17:00:00Z).
const test = require("node:test");
const assert = require("node:assert/strict");
const fs = require("node:fs");
const path = require("node:path");
const feed = require("../extension/feed.js");

const fixture = (name) => fs.readFileSync(path.join(__dirname, "fixtures", name), "utf8");
const snapshotJson = () => JSON.parse(fixture("v2_slice.json"));
const changesText = () => fixture("changes_slice.json");

const KICKOFF_MS = Date.parse("2026-09-13T17:00:00Z");
const BEFORE_KICKOFF = KICKOFF_MS - 3600 * 1000;
const NFL = 1;

function loadedState() {
  return feed.mergeStates([feed.parseSnapshot(snapshotJson(), { leagueId: NFL })]);
}

test("league key and event start parsing", () => {
  assert.deepEqual(feed.parseLeagueKey("lg1:pt2:pregame"), { leagueId: 1, periodTypeId: 2, phase: "pregame" });
  assert.equal(feed.parseLeagueKey("props:lg1"), null);
  assert.equal(feed.parseEventStart("2026-09-13T17:00:00"), KICKOFF_MS);
  assert.equal(feed.parseEventStart("2026-09-13T17:00:00+00:00"), KICKOFF_MS);
  assert.equal(feed.parseEventStart(null), null);
});

test("snapshot: game rows only, books keyed by id with the live flag", () => {
  const state = feed.parseSnapshot(snapshotJson(), { leagueId: NFL });
  // 3 full-game rows + 1 first-half row; the bt4 team-total row is skipped and
  // the whole lg1:pt1:live group is never visited.
  assert.equal(state.counts.rows, 4);
  assert.equal(state.counts.skippedRows, 1);
  assert.ok(!Object.values(state.lines).some((l) => l.marketId === 297848222));
  assert.equal(state.teams["6"], "Chicago Bears");
  assert.equal(state.teams["5"], "Carolina Panthers");
  assert.equal(state.books[4].name, "BetMGM");
  assert.equal(state.books[4].isLive, true);
  assert.equal(state.books[89].isLive, true); // Novig
  assert.equal(state.books[52].isLive, false); // Matchbook, isActive false
  assert.equal(state.books[69].isLive, false); // Sports Interaction, statusId 2
  assert.equal(state.books[49].isLive, false); // Unabated line
  const event = state.events[125807];
  assert.equal(event.eventStart, KICKOFF_MS);
  assert.equal(event.awayTeamId, 6);
  assert.equal(event.homeTeamId, 5);
  assert.equal(event.awayRotation, 465);
});

test("snapshot line carries ge, bacr, sourcePrice, liquidity and the side key", () => {
  const state = feed.parseSnapshot(snapshotJson(), { leagueId: NFL });
  const mgmMoneyline = state.lines["289357353:ms4:si0:tid6"];
  assert.equal(mgmMoneyline.price, -145);
  assert.equal(mgmMoneyline.ge, 0.0296);
  assert.equal(mgmMoneyline.bacr, -156);
  assert.equal(mgmMoneyline.betTypeId, 1);
  assert.equal(mgmMoneyline.periodTypeId, 1);
  assert.equal(mgmMoneyline.sideIndex, 0);
  assert.equal(mgmMoneyline.points, null);
  const kalshiMoneyline = state.lines["289357353:ms105:si0:tid6"];
  assert.equal(kalshiMoneyline.sourceFormat, 4);
  assert.equal(kalshiMoneyline.sourcePrice, 0.6168);
  assert.equal(kalshiMoneyline.liquidity, 83613.6);
  const offBoard = state.lines["289357350:ms52:si1:tid5"];
  assert.equal(offBoard.statusId, 2);
  assert.equal(offBoard.ge, null);
});

test("snapshot for a league the file does not carry fails loudly; an empty slate does not", () => {
  assert.throws(() => feed.parseSnapshot(snapshotJson(), { leagueId: 2 }), /no lg2 odds keys in the file/);
  assert.throws(() => feed.parseSnapshot({ nope: true }, { leagueId: 1 }), /expected an object with an `odds` map/);
  const offSeason = feed.parseSnapshot({ odds: { "lg5:pt1:pregame": [] }, teams: {}, marketSources: [] }, { leagueId: 5 });
  assert.equal(feed.countLines(offSeason), 0);
  assert.deepEqual(offSeason.leagues, [5]);
});

test("selectEdges: maxLineAgeMs drops lines the book has not touched, and lines with no modifiedOn", () => {
  const state = loadedState();
  const DAY = 86400 * 1000;
  // The three edges were last changed Aug 29-30; kickoff-1h is Sep 13.
  assert.equal(feed.selectEdges(state, { now: BEFORE_KICKOFF, maxLineAgeMs: 7 * DAY }).length, 0);
  assert.equal(feed.selectEdges(state, { now: BEFORE_KICKOFF, maxLineAgeMs: 30 * DAY }).length, 3);
  const row = feed.selectEdges(state, { now: BEFORE_KICKOFF })[0];
  assert.equal(row.modifiedMs, Date.parse("2026-08-29T14:40:50.899Z"));
  state.lines["289357360:ms99:si0:tid6"].modifiedOn = null;
  assert.equal(feed.selectEdges(state, { now: BEFORE_KICKOFF, maxLineAgeMs: 30 * DAY }).length, 2);
});

test("selectEdges skips a line whose price is not a valid American number", () => {
  const state = loadedState();
  state.lines["289357360:ms99:si0:tid6"].price = -95;
  assert.ok(!feed.selectEdges(state, { now: BEFORE_KICKOFF }).some((r) => r.key === "289357360:ms99:si0:tid6"));
});

test("selectEdges: live books, full game, >= 1%, sorted by edge, ticket-shaped rows", () => {
  const rows = feed.selectEdges(loadedState(), { now: BEFORE_KICKOFF });
  assert.deepEqual(rows.map((r) => [r.book.name, r.betType, r.sideLabel, r.edgePct]), [
    ["SouthPoint", "Spread", "Chicago Bears -2.5", 5.3],
    ["BetMGM", "Moneyline", "Chicago Bears", 2.96],
    ["SouthPoint", "Moneyline", "Chicago Bears", 1.56],
  ]);
  const top = rows[0];
  assert.equal(top.league, "nfl");
  assert.equal(top.leagueLabel, "NFL");
  assert.equal(top.awayTeam, "Chicago Bears");
  assert.equal(top.homeTeam, "Carolina Panthers");
  assert.equal(top.homeAway, "Away");
  assert.equal(top.rotation, 465);
  assert.equal(top.price, -110);
  assert.equal(top.fair, -123);
  assert.equal(top.points, -2.5);
  assert.equal(top.period, "FG");
  assert.equal(top.eventStart, "2026-09-13T17:00:00.000Z");
  // Matchbook's +5900 at -1.5 (ge 33.36) is a dead feed and must not appear: the book is not live.
  assert.ok(!rows.some((r) => r.book.id === 52));
});

test("selectEdges: started games drop out; user book filter replaces the live-book default", () => {
  const state = loadedState();
  assert.equal(feed.selectEdges(state, { now: KICKOFF_MS }).length, 0);
  const onlyMgm = feed.selectEdges(state, { now: BEFORE_KICKOFF, bookIds: new Set([4]) });
  assert.deepEqual(onlyMgm.map((r) => r.key), ["289357353:ms4:si0:tid6"]);
});

test("selectEdges: periods, bet types and the threshold are honoured", () => {
  const state = loadedState();
  const withFirstHalf = feed.selectEdges(state, { now: BEFORE_KICKOFF, minEdge: 0.005, periods: new Set([1, 2]) });
  const firstHalf = withFirstHalf.filter((r) => r.period === "1H");
  assert.deepEqual(firstHalf.map((r) => [r.book.name, r.sideLabel, r.price, r.edgePct]), [["Novig", "Over 23.5", 100, 0.5]]);
  const totalsOnly = feed.selectEdges(state, { now: BEFORE_KICKOFF, betTypes: new Set([3]) });
  assert.equal(totalsOnly.length, 0);
  const strict = feed.selectEdges(state, { now: BEFORE_KICKOFF, minEdge: 0.03 });
  assert.deepEqual(strict.map((r) => r.edgePct), [5.3]);
});

test("changes: cursor is read exactly from the text, lines are flattened per bet type", () => {
  const text = changesText();
  const parsed = feed.parseChanges(text);
  assert.equal(parsed.ok, true);
  assert.equal(parsed.cursor, "179164041314243100");
  assert.equal(feed.extractCursor(text), "179164041314243100");
  assert.equal(parsed.batches, 2);
  assert.equal(parsed.lines.length, 15);
  const mgmSpread = parsed.lines.find((l) => l.key === "289357360:ms4:si0:tid6");
  assert.equal(mgmSpread.points, -3.5);
  assert.equal(mgmSpread.price, -105);
  assert.equal(mgmSpread.ge, -0.0961);
  assert.equal(mgmSpread.sideIndex, 0);
  assert.equal(mgmSpread.eventStart, KICKOFF_MS);
  assert.equal(mgmSpread.betTypeId, 2);
  assert.equal(mgmSpread.periodTypeId, 1);
  const failed = feed.parseChanges(JSON.stringify({ latestTimestamp: 1, resultCode: "Failed", results: [] }));
  assert.equal(failed.ok, false);
  assert.throws(() => feed.parseChanges("{}"), /expected an object with a `results` array/);
});

test("cursorFromDate encodes whole seconds since 2021-01-06 in nanoseconds", () => {
  assert.equal(feed.cursorFromDate(new Date("2026-09-10T15:41:22.637Z")), "179163682000000000");
  assert.equal(feed.cursorFromDate(new Date("2020-01-01T00:00:00Z")), null);
});

test("applyChanges overwrites newer lines, adds new ones, skips other leagues and replays", () => {
  const state = loadedState();
  const before = state.lines["289357360:ms4:si0:tid6"];
  assert.equal(before.points, -14);
  assert.equal(before.price, 400);
  const changes = feed.parseChanges(changesText());
  const counts = feed.applyChanges(state, changes);
  assert.deepEqual(counts, { applied: 14, added: 4, stale: 0, unknownEvent: 0, otherLeague: 1 });
  const after = state.lines["289357360:ms4:si0:tid6"];
  assert.equal(after.points, -3.5);
  assert.equal(after.price, -105);
  assert.equal(after.ge, -0.0961);
  assert.equal(after.sequenceNumber, 1789055238590);
  assert.equal(after.liquidity, null);
  // A first-quarter line the snapshot slice never carried is now known.
  assert.equal(state.lines["366866367:ms4:si0:tid6"].periodTypeId, 4);
  // Replaying the same batch changes nothing.
  assert.deepEqual(feed.applyChanges(state, changes), { applied: 0, added: 0, stale: 14, unknownEvent: 0, otherLeague: 1 });
  // The edge list is unaffected: the moved lines were all negative edge.
  assert.equal(feed.selectEdges(state, { now: BEFORE_KICKOFF }).length, 3);
});

test("applyChanges ignores an older sequence number and lines for unknown events", () => {
  const state = loadedState();
  const changes = feed.parseChanges(changesText());
  const stale = { ...changes, lines: changes.lines.map((l) => ({ ...l, sequenceNumber: 1 })) };
  const counts = feed.applyChanges(state, stale);
  assert.equal(counts.stale, 10); // every line that already existed
  assert.equal(counts.added, 4);
  assert.equal(state.lines["289357360:ms4:si0:tid6"].points, -14);
  const unknown = { ...changes, lines: changes.lines.map((l) => ({ ...l, eventId: 1 })) };
  assert.equal(feed.applyChanges(loadedState(), unknown).unknownEvent, 14);
});

test("a changes line whose snapshot edge existed but is now null lists nothing", () => {
  const state = loadedState();
  const changes = feed.parseChanges(changesText());
  const southPointKey = "289357360:ms99:si0:tid6";
  const overwrite = { ...changes.lines[1], key: southPointKey, bookId: 99, ge: null, bacr: null, sequenceNumber: 9e12 };
  feed.applyChanges(state, { ...changes, lines: [overwrite] });
  assert.equal(state.lines[southPointKey].ge, null);
  assert.ok(!feed.selectEdges(state, { now: BEFORE_KICKOFF }).some((r) => r.key === southPointKey));
});
