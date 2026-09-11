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
  // statusId 2 is not "dead": Caesars/Underdog PM/Bet365 carry 2 or 3 while live (2026-09-11).
  assert.equal(state.books[69].isLive, true); // Sports Interaction, statusId 2
  assert.equal(state.books[49].isLive, true); // Unabated line: active, but selectEdges skips UNABATED_LINE_BOOK_ID
  const disabled = feed.parseSnapshot({ odds: {}, teams: {}, marketSources: [{ id: 7, name: "x", isActive: true, statusId: 1, isEnabledForGameOdds: false }] }, { leagueId: 1 });
  assert.equal(disabled.books[7].isLive, false);
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
  // Serie B on 2026-09-10: teams listed, `odds` an empty object.
  const noOdds = feed.parseSnapshot({ odds: {}, teams: { 1: { name: "Genoa CFC" } }, marketSources: [] }, { leagueId: 38 });
  assert.equal(feed.countLines(noOdds), 0);
});

test("selectEdges: maxLineAgeMs drops lines the book has not touched, and lines with no modifiedOn", () => {
  const state = loadedState();
  const DAY = 86400 * 1000;
  // The four edges were last changed Aug 29-30; kickoff-1h is Sep 13.
  assert.equal(feed.selectEdges(state, { now: BEFORE_KICKOFF, maxLineAgeMs: 7 * DAY }).length, 0);
  assert.equal(feed.selectEdges(state, { now: BEFORE_KICKOFF, maxLineAgeMs: 30 * DAY }).length, 4);
  const row = feed.selectEdges(state, { now: BEFORE_KICKOFF })[0];
  assert.equal(row.modifiedMs, Date.parse("2026-08-29T14:40:50.899Z"));
  state.lines["289357360:ms99:si0:tid6"].modifiedOn = null;
  assert.equal(feed.selectEdges(state, { now: BEFORE_KICKOFF, maxLineAgeMs: 30 * DAY }).length, 3);
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
    ["Sports Interaction", "Moneyline", "Chicago Bears", 2.96], // statusId 2, live (mirrors BetMGM)
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
  assert.equal(feed.selectEdges(state, { now: BEFORE_KICKOFF }).length, 4);
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

test("league table: every entry has a label, an odds-screen path and a sport; sports group them", () => {
  for (const [id, league] of Object.entries(feed.LEAGUES)) {
    assert.ok(Number.isInteger(Number(id)), `league id ${id}`);
    assert.ok(league.label && league.path && feed.SPORTS[league.sport], `league ${id} incomplete`);
  }
  assert.deepEqual(feed.leagueIdsOfSport("football"), [1, 2]);
  assert.ok(feed.leagueIdsOfSport("soccer").includes(28));
  assert.equal(feed.leagueIdsOfSport("tennis").length, 0);
});

// ---- alternate lines (issue #113) --------------------------------------------
// The fixture's alternateLines are real entries from the live NFL file of
// 2026-09-11 for the same event, slimmed to the fields the parser reads, with
// two deliberate edits: a null entry in Kalshi's away-spread ladder (nulls do
// appear live) and nothing else. Every alt's modifiedOn is the feed's
// "0001-01-01T00:00:00" sentinel, exactly as served.

const ALT_CHANGED_MS = 1789147708947; // Kalshi Bears -20.5 sequenceNumber = 2026-09-11T17:28:28.947Z

test("snapshot: alternateLines expand into alt-keyed lines under their main line", () => {
  const state = feed.parseSnapshot(snapshotJson(), { leagueId: NFL });
  assert.equal(state.counts.lines, 62);
  assert.equal(state.counts.altLines, 27);
  // 1 null entry + 3 alts sitting on their main line's own points (Novig
  // -2.5 / 2.5, Kalshi 46.5) + 5 alts on rows the parser never lists.
  assert.equal(state.counts.skippedAltLines, 9);
  assert.equal(feed.countLines(state), 62);
  assert.equal(feed.countAltLines(state), 27);
  const alt = state.lines["289357360:ms105:si0:tid6:alt-20.5"];
  assert.equal(alt.isAlt, true);
  assert.equal(alt.mainKey, "289357360:ms105:si0:tid6");
  assert.equal(alt.mainPoints, -2.5);
  assert.equal(alt.points, -20.5);
  assert.equal(alt.price, 944);
  assert.equal(alt.ge, 0.0642);
  assert.equal(alt.bacr, 881);
  assert.equal(alt.sourceFormat, 4);
  assert.equal(alt.sourcePrice, 0.095733);
  assert.equal(alt.liquidity, 264.51);
  assert.equal(alt.betTypeId, 2);
  assert.equal(alt.sideIndex, 0);
  assert.equal(state.lines["289357360:ms105:si0:tid6"].isAlt, false);
  // No alt sits on the main line's own points, and moneylines carry none.
  assert.equal(state.lines["289357360:ms89:si0:tid6:alt-2.5"], undefined);
  assert.ok(!Object.values(state.lines).some((l) => l.isAlt && l.betTypeId === 1));
});

test("an alt's change time is its sequenceNumber; the sentinel modifiedOn is not a date", () => {
  const state = feed.parseSnapshot(snapshotJson(), { leagueId: NFL });
  const alt = state.lines["289357360:ms105:si0:tid6:alt-20.5"];
  assert.equal(alt.modifiedOn, "0001-01-01T00:00:00");
  assert.equal(feed.parseModifiedOn(alt.modifiedOn), null);
  assert.equal(feed.lineChangedMs(alt), ALT_CHANGED_MS);
  // A main line still reads modifiedOn only: a counter-sized sequence is never a clock.
  const main = state.lines["289357360:ms105:si0:tid6"];
  assert.equal(feed.lineChangedMs(main), Date.parse("2026-09-10T15:39:34.824Z"));
  assert.equal(feed.lineChangedMs({ ...main, modifiedOn: null }), null);
  assert.equal(feed.lineChangedMs({ ...alt, sequenceNumber: 12345 }), null);
  assert.equal(feed.lineChangedMs({ ...alt, sequenceNumber: null }), null);
});

test("selectEdges lists no alt unless includeAlts is on", () => {
  const state = loadedState();
  assert.equal(feed.selectEdges(state, { now: BEFORE_KICKOFF }).some((r) => r.isAlt), false);
  const rows = feed.selectEdges(state, { now: BEFORE_KICKOFF, includeAlts: true });
  const alts = rows.filter((r) => r.isAlt);
  assert.equal(rows.length - alts.length, 4); // the four main-line edges are still there
  assert.deepEqual(alts.slice(0, 4).map((r) => [r.book.name, r.sideLabel, r.price, r.edgePct, r.mainPoints]), [
    ["Kalshi", "Carolina Panthers -9.5", 625, 12.58, 2.5],
    ["Kalshi", "Over 64.5", 840, 9.56, 46.5],
    ["Kalshi", "Over 61.5", 573, 6.83, 46.5],
    ["Kalshi", "Chicago Bears -20.5", 944, 6.42, -2.5],
  ]);
  assert.equal(alts.length, 16);
  // Matchbook's Over 8.5 at +112 (ge 1.1193) is a dead feed: the book is not live.
  assert.ok(!alts.some((r) => r.book.id === 52));
  const top = alts[0];
  assert.equal(top.isAlt, true);
  assert.equal(top.key, "289357357:ms105:si1:tid5:alt-9.5");
  assert.equal(top.marketId, 289357357);
  assert.equal(top.modifiedMs, 1789147708929);
  assert.equal(rows.find((r) => !r.isAlt).mainPoints, null);
});

test("altMaxDistance keeps alts within N points of the book's main number", () => {
  const state = loadedState();
  const within7 = feed.selectEdges(state, { now: BEFORE_KICKOFF, includeAlts: true, altMaxDistance: 7 }).filter((r) => r.isAlt);
  assert.deepEqual(within7.map((r) => [r.sideLabel, r.mainPoints]), [
    ["Carolina Panthers -2.5", 2.5],
    ["Carolina Panthers -4.5", 2.5],
    ["Over 54.5", 47.5], // exactly 7 away is kept
    ["Over 50.5", 47.5],
    ["Chicago Bears -4.5", -2.5],
    ["Chicago Bears -9.5", -2.5],
  ]);
  const within2 = feed.selectEdges(state, { now: BEFORE_KICKOFF, includeAlts: true, altMaxDistance: 2 }).filter((r) => r.isAlt);
  assert.deepEqual(within2.map((r) => r.sideLabel), ["Chicago Bears -4.5"]);
  // 0 or a non-number means no distance gate.
  assert.equal(feed.selectEdges(state, { now: BEFORE_KICKOFF, includeAlts: true, altMaxDistance: 0 }).filter((r) => r.isAlt).length, 16);
});

test("altMinLiquidity drops thin exchange alts and leaves books with no liquidity figure alone", () => {
  const state = loadedState();
  const rows = feed.selectEdges(state, { now: BEFORE_KICKOFF, includeAlts: true, altMinLiquidity: 100 }).filter((r) => r.isAlt);
  assert.equal(rows.length, 14);
  // Kalshi Panthers -9.5 ($35 resting) and Over 64.5 ($86) are gone; Novig's alts carry no liquidity and stay.
  assert.ok(!rows.some((r) => r.key === "289357357:ms105:si1:tid5:alt-9.5"));
  assert.ok(!rows.some((r) => r.key === "289357345:ms105:si0:tid6:alt64.5"));
  assert.ok(rows.some((r) => r.key === "289357345:ms105:si0:tid6:alt61.5"));
  assert.ok(rows.some((r) => r.key === "289357357:ms89:si1:tid5:alt-2.5"));
});

test("an alt is hidden while the main line sits on its number, and distance follows the moved main line", () => {
  const state = loadedState();
  const opts = { now: BEFORE_KICKOFF, includeAlts: true, altMaxDistance: 7 };
  assert.ok(feed.selectEdges(state, opts).some((r) => r.key === "289357360:ms89:si0:tid6:alt-4.5"));
  // Novig moves its Bears main line from -2.5 to -4.5 (a changes-stream update).
  const main = state.lines["289357360:ms89:si0:tid6"];
  feed.applyChanges(state, { lines: [{ ...main, points: -4.5, price: 130, sequenceNumber: main.sequenceNumber + 1, eventStart: null }] });
  assert.equal(state.lines["289357360:ms89:si0:tid6"].points, -4.5);
  const rows = feed.selectEdges(state, opts);
  assert.ok(!rows.some((r) => r.key === "289357360:ms89:si0:tid6:alt-4.5"));
  // -9.5 is now 5 from the main number; -13.5 (9 away) is still out; mainPoints reports the current main.
  const nineHalf = rows.find((r) => r.key === "289357360:ms89:si0:tid6:alt-9.5");
  assert.equal(nineHalf.mainPoints, -4.5);
  assert.ok(!rows.some((r) => r.key === "289357360:ms89:si0:tid6:alt-13.5"));
  // The stream never carries alts, so the alt itself is untouched.
  assert.equal(state.lines["289357360:ms89:si0:tid6:alt-9.5"].price, 245);
});

test("maxLineAgeMs applies to alts through their sequenceNumber", () => {
  const state = loadedState();
  const DAY = 86400 * 1000;
  // Alts changed 2026-09-11; kickoff-1h is 2026-09-13T16:00Z, so they are ~46h old.
  const fresh = feed.selectEdges(state, { now: BEFORE_KICKOFF, includeAlts: true, maxLineAgeMs: 3 * DAY }).filter((r) => r.isAlt);
  assert.equal(fresh.length, 16);
  const strict = feed.selectEdges(state, { now: BEFORE_KICKOFF, includeAlts: true, maxLineAgeMs: DAY }).filter((r) => r.isAlt);
  assert.equal(strict.length, 0);
  state.lines["289357357:ms105:si1:tid5:alt-9.5"].sequenceNumber = null;
  const unknowable = feed.selectEdges(state, { now: BEFORE_KICKOFF, includeAlts: true, maxLineAgeMs: 3 * DAY }).filter((r) => r.isAlt);
  assert.equal(unknowable.length, 15);
  assert.ok(!unknowable.some((r) => r.key === "289357357:ms105:si1:tid5:alt-9.5"));
});

test("alt lines honour the same board, book, bet-type, period and price gates as main lines", () => {
  const state = loadedState();
  state.lines["289357357:ms105:si1:tid5:alt-9.5"].statusId = 2;
  const rows = feed.selectEdges(state, { now: BEFORE_KICKOFF, includeAlts: true });
  assert.ok(!rows.some((r) => r.key === "289357357:ms105:si1:tid5:alt-9.5"));
  const totalsOnly = feed.selectEdges(state, { now: BEFORE_KICKOFF, includeAlts: true, betTypes: new Set([3]) });
  assert.ok(totalsOnly.every((r) => r.betType === "Total"));
  assert.equal(totalsOnly.filter((r) => r.isAlt).length, 4);
  const novigOnly = feed.selectEdges(state, { now: BEFORE_KICKOFF, includeAlts: true, bookIds: new Set([89]) });
  assert.ok(novigOnly.every((r) => r.book.id === 89));
  assert.equal(feed.selectEdges(state, { now: KICKOFF_MS, includeAlts: true }).length, 0);
});

test("a ladder the book pulls is gone from the next snapshot parse", () => {
  const first = feed.parseSnapshot(snapshotJson(), { leagueId: NFL });
  const again = snapshotJson();
  const spreadRow = again.odds["lg1:pt1:pregame"].find((row) => row.key === "pt1:pregame:bt2:e125807");
  spreadRow.sides["si0:tid6"].ms105.alternateLines = [];
  const second = feed.parseSnapshot(again, { leagueId: NFL });
  assert.equal(feed.countAltLines(second), 24);
  assert.equal(second.lines["289357360:ms105:si0:tid6:alt-20.5"], undefined);
});

// ---- grouping by market ---------------------------------------------------------

const kellyForTests = require("../extension/kelly.js");
const stakeOf = (row) => kellyForTests.kellyStakeFromEdge({ bookPrice: row.price, edgePct: row.edgePct, bankroll: 30000, multiplier: 0.25 }).stake;

test("groupEdges: one card per (game, period, bet type, side) with books and lines counted", () => {
  const rows = feed.selectEdges(loadedState(), { now: BEFORE_KICKOFF, includeAlts: true });
  const groups = feed.groupEdges(rows);
  assert.deepEqual(groups.map((g) => [g.key, g.sideName, g.betType, g.bookCount, g.rows.length]), [
    ["125807:pt1:bt2:si1", "Carolina Panthers", "Spread", 2, 6],
    ["125807:pt1:bt3:si0", "Over", "Total", 2, 4],
    ["125807:pt1:bt2:si0", "Chicago Bears", "Spread", 3, 7],
    ["125807:pt1:bt1:si0", "Chicago Bears", "Moneyline", 3, 3],
  ]);
  assert.equal(rows.length, groups.reduce((n, g) => n + g.rows.length, 0));
  const bears = groups[2];
  assert.equal(bears.league, "nfl");
  assert.equal(bears.awayTeam, "Chicago Bears");
  assert.equal(bears.homeTeam, "Carolina Panthers");
  assert.equal(bears.eventStartMs, KICKOFF_MS);
  assert.equal(bears.period, "FG");
  assert.equal(feed.groupKeyOf(bears.best), bears.key);
});

test("groupEdges: with no rank function the best line is the highest edge; cards follow their best", () => {
  const rows = feed.selectEdges(loadedState(), { now: BEFORE_KICKOFF, includeAlts: true });
  const groups = feed.groupEdges(rows);
  assert.deepEqual(groups.map((g) => [g.best.sideLabel, g.best.book.name, g.best.edgePct]), [
    ["Carolina Panthers -9.5", "Kalshi", 12.58],
    ["Over 64.5", "Kalshi", 9.56],
    ["Chicago Bears -20.5", "Kalshi", 6.42],
    ["Chicago Bears", "BetMGM", 2.96],
  ]);
  assert.ok(groups[2].rows.every((row, i, all) => i === 0 || all[i - 1].edgePct >= row.edgePct));
});

test("groupEdges ranked by Kelly stake picks the bettable main line over the +944 rung", () => {
  const rows = feed.selectEdges(loadedState(), { now: BEFORE_KICKOFF, includeAlts: true }).map((row) => ({ ...row, stake: stakeOf(row) }));
  const groups = feed.groupEdges(rows, (row) => row.stake);
  assert.deepEqual(groups.map((g) => [g.sideName, g.betType, g.best.sideLabel, g.best.book.name, g.best.isAlt]), [
    ["Chicago Bears", "Spread", "Chicago Bears -2.5", "SouthPoint", false],
    ["Chicago Bears", "Moneyline", "Chicago Bears", "BetMGM", false],
    ["Carolina Panthers", "Spread", "Carolina Panthers -2.5", "Novig", true],
    ["Over", "Total", "Over 61.5", "Kalshi", true],
  ]);
  // Inside a card the lines fall by stake, and the rows keep their extra fields.
  const bears = groups[0];
  assert.ok(bears.rows.every((row, i, all) => i === 0 || all[i - 1].stake >= row.stake));
  assert.equal(typeof bears.rows[3].stake, "number");
});

test("groupEdges: alts off collapses the four main-line edges to two cards; empty in, empty out; a null rank sorts last", () => {
  const rows = feed.selectEdges(loadedState(), { now: BEFORE_KICKOFF });
  const groups = feed.groupEdges(rows);
  // SouthPoint -2.5 alone; BetMGM, Sports Interaction and SouthPoint moneylines share a card.
  assert.deepEqual(groups.map((g) => [g.betType, g.rows.length, g.bookCount]), [["Spread", 1, 1], ["Moneyline", 3, 3]]);
  assert.deepEqual(feed.groupEdges([]), []);
  const ranked = feed.groupEdges(rows, (row) => (row.book.id === 99 ? null : row.edgePct));
  assert.equal(ranked[ranked.length - 1].best.book.id, 99);
});
