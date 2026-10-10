// Run: node --test unabated_ticket/tests
// live.js: page.js's live_edges payloads -> the Edges tab's Live block.
const test = require("node:test");
const assert = require("node:assert/strict");
const live = require("../extension/live.js");
const fs = require("node:fs");
const path = require("node:path");
const teams = require("../extension/teams.js");
const betsLib = require("../extension/bets.js");
const edgeRows = require("../extension/edgerows.js");
// Real Unabated team ids (teams_index.json), so a bet's team names resolve.
teams.loadIndex(JSON.parse(fs.readFileSync(path.join(__dirname, "fixtures", "teams_index.json"), "utf8")).leagues);

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
    row({ key: "k-alt-deep", isAlt: true, mainPoints: 7.5, points: 24.5, fair: 1500 }),
    row({ key: "k-thin-liq", liquidity: 50 }),
  ];
  const keys = (opts) => live.selectLiveEdges([payload({ rows })], { now: NOW, ...opts }).map((r) => r.key).sort();
  assert.deepEqual(keys({ minEdgePct: 1, bookIds: new Set([1]), betTypes: new Set([2]), periods: new Set([1]), minLiquidityToWin: 100 }), ["live:125807:pt1:bt2:si1:ms1"]);
  // Alts are capped by the fair (10-90%), not points: the 10-point rung at -124 lists, the +1500 one does not.
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

test("liveFairState: every line of a ready game as a ladder line, the in-game fair as its bacr; expired games and old payloads give none", () => {
  const ready = payload({ rows: [row(), row({ key: "k-under", sideIndex: 0, points: -7.5, fair: 124, edgePct: -3.1 })] });
  const expired = payload({ league: "cfb", games: [game({ eventId: 7, fairStatus: "expired" })], rows: [row({ key: "k-expired", eventId: 7 })] });
  const old = payload({ at: NOW - live.LIVE_DROP_MS - 1, rows: [row({ key: "k-old" })] });
  assert.deepEqual(live.liveFairState([ready, expired, old], NOW), {
    lines: {
      "live:125807:pt1:bt2:si1:ms1": { eventId: 125807, leagueId: 1, periodTypeId: 1, betTypeId: 2, sideIndex: 1, points: 7.5, bacr: -124 },
      "k-under": { eventId: 125807, leagueId: 1, periodTypeId: 1, betTypeId: 2, sideIndex: 0, points: -7.5, bacr: 124 },
    },
  });
});

// ---- live stakes sized against the open bets (conditional Kelly on the in-game ladder)

const EAGLES_ID = 24;
const BEARS_ID = 6;

// Novig's full-game total at halftime: main 44.5, alt rungs 41.5 and 47.5,
// both sides, fair P(over) 62% / 50% / 38%.
function totalRow(sideIndex, points, fair, price, edgePct) {
  const isAlt = points !== 44.5;
  return row({
    key: `live:125807:pt1:bt3:si${sideIndex}:ms89${isAlt ? `:alt${points}` : ""}`, betTypeId: 3, sideIndex, sideKey: `si${sideIndex}`,
    sideLabel: `${sideIndex === 0 ? "Over" : "Under"} ${points}`, bookId: 89, bookName: "Novig",
    isAlt, mainPoints: isAlt ? 44.5 : null, points, price, fair, edgePct, marketId: 3,
  });
}

function halftimeTotals() {
  return payload({
    games: [game({ awayTeamId: EAGLES_ID, homeTeamId: BEARS_ID, eventPeriodTypeId: 5 })],
    rows: [
      totalRow(0, 41.5, -163, -170, -1.5), totalRow(1, 41.5, 163, 150, -4.9),
      totalRow(0, 44.5, 100, -105, -2.4), totalRow(1, 44.5, 100, -105, -2.4),
      totalRow(0, 47.5, 163, 185, 8.4), totalRow(1, 47.5, -163, -210, -4.2),
    ],
  });
}

// A $300 Kalshi total bet placed at the break. Kalshi names no start time,
// only the Eastern date, and the game has left the pregame board.
function kalshiTotalHeld(side, points) {
  return betsLib.resolveTeamKeys([{
    id: `kalshi:${side}${points}`, source: "kalshi_api", venue: "kalshi", league: "nfl", status: "open", sourceFetchedAt: "2026-09-28T00:50:00Z",
    awayTeam: "PHI Eagles", homeTeam: "CHI Bears", awayKey: null, homeKey: null,
    betType: "total", period: "FG", side, points, price: 100, stake: 300, toWin: 300,
    placedAt: "2026-09-28T00:40:00Z", closedAt: null, eventStart: null, eventDate: "2026-09-27", isParlayLeg: false, unmatchable: null, approx: [],
    venueIds: { eventTicker: "KXNFLTOTAL-26SEP27PHICHI", marketTicker: `KXNFLTOTAL-26SEP27PHICHI-${Math.floor(points)}` },
  }]);
}

// The Over 47.5 (+185 at a +163 fair) sized against `records` on the in-game ladder.
function sizedOver475(records) {
  const listed = live.selectLiveEdges([halftimeTotals()], { now: NOW, includeAlts: true, minEdgePct: 1 });
  assert.deepEqual(listed.map((r) => r.key), ["live:125807:pt1:bt3:si0:ms89:alt47.5"]);
  const [sized] = edgeRows.withBetFlags(listed, {
    records, boardLines: [], stakeSettings: { bankroll: 30000, multiplier: 0.25 },
    ladderReaderOf: edgeRows.createLadderReaders(live.liveFairState([halftimeTotals()], NOW)), teasers: [],
  });
  const { kind, bet, alone, verb, held, against, reason } = sized.bet.advice;
  return {
    tier: sized.bet.tier, advice: { kind, bet, alone, verb, held, against, reason },
    matches: sized.bet.matches.map((match) => ({ id: match.bet.id, tier: match.tier, inMath: match.inMath, note: match.note })),
  };
}

const OVER_475_ALONE = 340.5405405405406;

test("live sizing: a held Over 44.5 cuts the Over 47.5 from $341 to $167", () => {
  assert.deepEqual(sizedOver475(kalshiTotalHeld("over", 44.5)), {
    tier: "same_side",
    advice: { kind: "sized", bet: 167.18, alone: OVER_475_ALONE, verb: "add", held: 300, against: 0, reason: null },
    matches: [{ id: "kalshi:over44.5", tier: "same_side", inMath: true, note: null }],
  });
});

test("live sizing: a held Under 44.5 lifts the Over 47.5 to $508", () => {
  assert.deepEqual(sizedOver475(kalshiTotalHeld("under", 44.5)), {
    tier: "opposite",
    advice: { kind: "sized", bet: 507.93, alone: OVER_475_ALONE, verb: "bet", held: 0, against: 300, reason: null },
    matches: [{ id: "kalshi:under44.5", tier: "opposite", inMath: true, note: null }],
  });
});

test("live sizing: a held bet at a number the screen has no fair for is named and left out", () => {
  assert.deepEqual(sizedOver475(kalshiTotalHeld("over", 51.5)), {
    tier: "same_side",
    advice: { kind: "none", bet: OVER_475_ALONE, alone: OVER_475_ALONE, verb: "bet", held: 0, against: 0, reason: null },
    matches: [{ id: "kalshi:over51.5", tier: "same_side", inMath: false, note: "no fair at 51.5" }],
  });
});
