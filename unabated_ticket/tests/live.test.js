// Run: node --test unabated_ticket/tests
// live.js: page.js's live_edges payloads -> the Edges tab's Live block.
const test = require("node:test");
const assert = require("node:assert/strict");
const live = require("../extension/live.js");

const NOW = Date.parse("2026-09-28T00:53:00Z");

function game(overrides = {}) {
  return {
    eventId: 125807, eventName: "Eagles Philadelphia PHI @ Bears Chicago CHI", eventStart: "2026-09-28T00:15:00",
    awayTeam: "Philadelphia Eagles", homeTeam: "Chicago Bears", awayTeamId: 21, homeTeamId: 6,
    fairStatus: "ready", checkpointType: "EndOfPeriod", producedUtc: "2026-09-28T00:52:10.000",
    eventPeriodTypeId: 4, gameClock: "00:00",
    ...overrides,
  };
}

function row(overrides = {}) {
  return {
    key: "live:125807:pt1:bt2:si1:ms1", eventId: 125807, eventName: null, periodTypeId: 1, betTypeId: 2,
    sideIndex: 1, sideKey: "si1:tid6", sideLabel: "Chicago Bears +7.5", bookId: 1, bookName: "DraftKings",
    isAlt: false, mainPoints: null, points: 7.5, price: -110, sourceFormat: 1, sourcePrice: null,
    fair: -124, edgePct: 5.8, liquidity: null, modifiedOn: "2026-09-28T00:52:48.000", marketId: 1,
    ...overrides,
  };
}

function payload(overrides = {}) {
  return { league: "nfl", at: NOW - 1000, error: null, games: [game()], rows: [row()], ...overrides };
}

test("checkpointLabel: the period that just ended, halftime for the 2nd quarter", () => {
  assert.equal(live.checkpointLabel({ eventPeriodTypeId: 4 }), "end of Q1");
  assert.equal(live.checkpointLabel({ eventPeriodTypeId: 5 }), "halftime");
  assert.equal(live.checkpointLabel({ eventPeriodTypeId: 6 }), "end of Q3");
  assert.equal(live.checkpointLabel({ eventPeriodTypeId: null, checkpointType: "EndOfPeriod" }), "end of period");
  assert.equal(live.checkpointLabel({}), "break");
});

test("selectLiveEdges: an Edges row with the screen's edge and fair, a live block, and the league's id", () => {
  const [edge] = live.selectLiveEdges([payload()], { now: NOW });
  assert.equal(edge.league, "nfl");
  assert.equal(edge.leagueId, 1);
  assert.equal(edge.betType, "Spread");
  assert.equal(edge.period, "FG");
  assert.equal(edge.edgePct, 5.8);
  assert.equal(edge.fair, -124);
  assert.equal(edge.book.name, "DraftKings");
  assert.equal(edge.live.checkpoint, "end of Q1");
  assert.equal(edge.live.matchup, "Philadelphia Eagles @ Chicago Bears");
  // Naive UTC read as UTC: 00:52:48 is 12 s before NOW.
  assert.equal(NOW - edge.modifiedMs, 12000);
  assert.equal(edge.eventStartMs, Date.parse("2026-09-28T00:15:00Z"));
});

test("selectLiveEdges: only a ready fair lists; an old payload lists nothing", () => {
  assert.equal(live.selectLiveEdges([payload({ games: [game({ fairStatus: "expired" })] })], { now: NOW }).length, 0);
  assert.equal(live.selectLiveEdges([payload({ at: NOW - live.LIVE_DROP_MS - 1 })], { now: NOW }).length, 0);
});

test("selectLiveEdges: the Edges tab's filters apply", () => {
  const rows = [
    row(),
    row({ key: "k-thin", edgePct: 0.8 }),
    row({ key: "k-book", bookId: 2, bookName: "FanDuel" }),
    row({ key: "k-total", betTypeId: 3 }),
    row({ key: "k-1q", periodTypeId: 4 }),
    row({ key: "k-alt-near", isAlt: true, mainPoints: 7.5, points: 10.5 }),
    row({ key: "k-alt-far", isAlt: true, mainPoints: 7.5, points: 17.5 }),
    row({ key: "k-thin-liq", liquidity: 50 }),
  ];
  const keys = (opts) => live.selectLiveEdges([payload({ rows })], { now: NOW, ...opts }).map((r) => r.key).sort();
  assert.deepEqual(keys({ minEdgePct: 1, bookIds: new Set([1]), betTypes: new Set([2]), periods: new Set([1]), minLiquidityToWin: 100 }), ["live:125807:pt1:bt2:si1:ms1"]);
  // No distance cap on alts (removed 2026-09-30): the 10-point rung lists too.
  assert.deepEqual(keys({ includeAlts: true, bookIds: new Set([1]), betTypes: new Set([2]), periods: new Set([1]) }),
    ["k-alt-far", "k-alt-near", "k-thin", "k-thin-liq", "live:125807:pt1:bt2:si1:ms1"]);
});

test("selectLiveEdges: best edge first across leagues", () => {
  const nfl = payload({ rows: [row({ key: "a", edgePct: 2 }), row({ key: "b", edgePct: 6 })] });
  const cfb = payload({ league: "cfb", games: [game({ eventId: 7 })], rows: [row({ key: "c", eventId: 7, edgePct: 4 })] });
  assert.deepEqual(live.selectLiveEdges([nfl, cfb], { now: NOW }).map((r) => r.key), ["b", "c", "a"]);
});

test("liveView: games with their checkpoint, stale past two heartbeats, hidden once dropped", () => {
  const fresh = live.liveView([payload()], NOW);
  assert.equal(fresh.games.length, 1);
  assert.equal(fresh.games[0].fairReady, true);
  assert.equal(fresh.games[0].checkpoint, "end of Q1");
  assert.equal(fresh.stale, false);
  const quiet = live.liveView([payload({ at: NOW - live.LIVE_STALE_MS - 1 })], NOW);
  assert.equal(quiet.stale, true);
  assert.equal(live.liveView([payload({ at: NOW - live.LIVE_DROP_MS - 1 })], NOW).games.length, 0);
  assert.equal(live.liveView([], NOW).games.length, 0);
});

test("liveAlertKey: one alert per line per break", () => {
  const [first] = live.selectLiveEdges([payload()], { now: NOW });
  const [again] = live.selectLiveEdges([payload()], { now: NOW });
  const [nextBreak] = live.selectLiveEdges([payload({ games: [game({ producedUtc: "2026-09-28T01:30:00.000" })] })], { now: NOW });
  assert.equal(live.liveAlertKey(first), live.liveAlertKey(again));
  assert.notEqual(live.liveAlertKey(first), live.liveAlertKey(nextBreak));
});

test("storageKeyOf matches what content.js writes", () => {
  assert.equal(live.storageKeyOf("nfl"), "liveEdges:nfl");
  assert.ok(live.storageKeyOf("nfl").startsWith(live.STORAGE_PREFIX));
});

test("liveView: a hidden tab gets a minute and a margin before it reads stale (Chrome throttles its timers)", () => {
  const hidden = (ageMs) => live.liveView([payload({ visible: false, at: NOW - ageMs })], NOW).stale;
  assert.equal(hidden(30 * 1000), false);
  assert.equal(hidden(live.LIVE_STALE_HIDDEN_MS + 1), true);
  assert.equal(live.liveView([payload({ visible: true, at: NOW - 30 * 1000 })], NOW).stale, true);
});
