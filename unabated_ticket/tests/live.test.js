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
const kelly = require("../extension/kelly.js");
const betsView = require("../extension/betsview.js");
// Real Unabated team ids (teams_index.json), so a bet's team names resolve.
teams.loadIndex(JSON.parse(fs.readFileSync(path.join(__dirname, "fixtures", "teams_index.json"), "utf8")).leagues);

const NOW = Date.parse("2026-09-28T00:53:00Z");

function game(overrides = {}) {
  return {
    eventId: 125807, eventName: "Eagles Philadelphia PHI @ Bears Chicago CHI", eventStart: "2026-09-28T00:15:00",
    awayTeam: "Philadelphia Eagles", homeTeam: "Chicago Bears", awayTeamId: 24, homeTeamId: 6,
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
  // End of Q1: the first quarter is decided, nothing else.
  assert.deepEqual([edge.live.eventPeriodTypeId, edge.live.settledPeriods], [4, ["1Q"]]);
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

// One Novig line at halftime: the main number, or an alt rung off it.
function novigRow(betTypeId, sideIndex, points, mainPoints, fair, price, edgePct) {
  const isAlt = points !== mainPoints;
  return row({
    key: `live:125807:pt1:bt${betTypeId}:si${sideIndex}:ms89${isAlt ? `:alt${points}` : ""}`, betTypeId, sideIndex, sideKey: `si${sideIndex}`,
    sideLabel: null, bookId: 89, bookName: "Novig", isAlt, mainPoints: isAlt ? mainPoints : null, points, price, fair, edgePct, marketId: betTypeId,
  });
}

// Both sides of three numbers, fair 62% / 50% / 38% past them: the total
// around 44.5 and the Eagles spread around -7.5. The candidate is the far
// number's first side at +185 against a +163 fair.
function halftimeLadders() {
  const total = (sideIndex, points, fair, price, edgePct) => novigRow(3, sideIndex, points, 44.5, fair, price, edgePct);
  const spread = (sideIndex, points, fair, price, edgePct) => novigRow(2, sideIndex, points, sideIndex === 0 ? -7.5 : 7.5, fair, price, edgePct);
  return payload({
    games: [game({ eventPeriodTypeId: 5 })],
    rows: [
      total(0, 41.5, -163, -170, -1.5), total(1, 41.5, 163, 150, -4.9),
      total(0, 44.5, 100, -105, -2.4), total(1, 44.5, 100, -105, -2.4),
      total(0, 47.5, 163, 185, 8.4), total(1, 47.5, -163, -210, -4.2),
      spread(0, -4.5, -163, -170, -1.5), spread(1, 4.5, 163, 150, -4.9),
      spread(0, -7.5, 100, -105, -2.4), spread(1, 7.5, 100, -105, -2.4),
      spread(0, -10.5, 163, 185, 8.4), spread(1, 10.5, -163, -210, -4.2),
    ],
  });
}

// A $300 open bet on the game at even money, team keys resolved. Kalshi
// names no start time, only the Eastern date; the game has left the
// pregame board (boardLines is empty below).
function heldBet(fields) {
  return betsLib.resolveTeamKeys([{
    source: "kalshi_api", venue: "kalshi", league: "nfl", status: "open", sourceFetchedAt: "2026-09-28T00:50:00Z",
    awayTeam: "PHI Eagles", homeTeam: "CHI Bears", awayKey: null, homeKey: null,
    betType: "total", period: "FG", price: 100, stake: 300, toWin: 300,
    placedAt: "2026-09-28T00:40:00Z", closedAt: null, eventStart: null, eventDate: "2026-09-27", isParlayLeg: false, unmatchable: null, approx: [],
    venueIds: { eventTicker: "KXNFLTOTAL-26SEP27PHICHI", marketTicker: "KXNFLTOTAL-26SEP27PHICHI-44" },
    ...fields,
  }]);
}

// The listed row with `key`, sized against `records` on the in-game ladder
// the way panel.js does it (a bounded ladder reader over the live payload).
function sizedLiveRow(key, records, livePayload = halftimeLadders()) {
  const listed = live.selectLiveEdges([livePayload], { now: NOW, includeAlts: true, minEdgePct: 1 });
  const [sized] = edgeRows.withBetFlags(listed.filter((r) => r.key === key), {
    records, boardLines: [], stakeSettings: edgeRows.DEFAULT_STAKE_SETTINGS,
    ladderReaderOf: edgeRows.createLadderReaders(live.liveFairState([livePayload], NOW), { bounded: true }), teasers: [],
  });
  return sized;
}

// That row's flag, the matches cut down to what the tests read.
function sizedLive(key, records, livePayload = halftimeLadders()) {
  const sized = sizedLiveRow(key, records, livePayload);
  return {
    tier: sized.bet.tier, advice: { ...sized.bet.advice, matches: undefined },
    matches: sized.bet.matches.map((match) => ({ id: match.bet.id, tier: match.tier, inMath: match.inMath, note: match.note, ...(match.bound ? { bound: match.bound } : {}) })),
  };
}

const OVER_475 = "live:125807:pt1:bt3:si0:ms89:alt47.5";
const EAGLES_MINUS_105 = "live:125807:pt1:bt2:si0:ms89:alt-10.5";
const ALONE = kelly.kellyStakeFromEdge({ bookPrice: 185, edgePct: 8.4, ...edgeRows.DEFAULT_STAKE_SETTINGS }).stake;

function sizedAdvice(bet, verb, held, against) {
  const kind = bet === ALONE ? "none" : "sized";
  return {
    kind, bet, alone: ALONE, verb, held, against, teasers: { held: 0, against: 0 }, reason: null,
    matches: undefined, teaserGroups: [], cappedAt: null, uncapped: bet,
  };
}

test("live sizing: the listed rows are the two far numbers", () => {
  assert.deepEqual(live.selectLiveEdges([halftimeLadders()], { now: NOW, includeAlts: true, minEdgePct: 1 }).map((r) => r.key).sort(), [EAGLES_MINUS_105, OVER_475]);
  assert.equal(ALONE, 340.5405405405406);
});

test("live sizing: a held Over 44.5 (Kalshi, placed at the break) cuts the Over 47.5 from $341 to $167", () => {
  assert.deepEqual(sizedLive(OVER_475, heldBet({ id: "k-over", side: "over", points: 44.5 })), {
    tier: "same_side", advice: sizedAdvice(167.18, "add", 300, 0),
    matches: [{ id: "k-over", tier: "same_side", inMath: true, note: null }],
  });
});

test("live sizing: a held Under 44.5 lifts the Over 47.5 to $508", () => {
  assert.deepEqual(sizedLive(OVER_475, heldBet({ id: "k-under", side: "under", points: 44.5 })), {
    tier: "opposite", advice: sizedAdvice(507.93, "bet", 0, 300),
    matches: [{ id: "k-under", tier: "opposite", inMath: true, note: null }],
  });
});

test("live sizing: a pregame Over 44.5 (Novig, with a start time) counts the same as a live one", () => {
  const pregame = heldBet({
    id: "n-over", source: "novig_rest", venue: "novig", awayTeam: "Philadelphia Eagles", homeTeam: "Chicago Bears",
    side: "over", points: 44.5, placedAt: "2026-09-27T18:00:00Z", eventStart: "2026-09-28T00:15:00.000Z", venueIds: null,
  });
  assert.deepEqual(sizedLive(OVER_475, pregame), {
    tier: "same_side", advice: sizedAdvice(167.18, "add", 300, 0),
    matches: [{ id: "n-over", tier: "same_side", inMath: true, note: null }],
  });
});

test("live sizing: spreads read the same way: a held Eagles -7.5 cuts the Eagles -10.5 to $167, a held Bears +7.5 lifts it to $508", () => {
  const spreadBet = (id, side, points) => heldBet({
    id, betType: "spread", side, points, venueIds: { eventTicker: "KXNFLSPREAD-26SEP27PHICHI", marketTicker: "KXNFLSPREAD-26SEP27PHICHI-PHI7" },
  });
  assert.deepEqual(sizedLive(EAGLES_MINUS_105, spreadBet("k-eagles", "away", -7.5)), {
    tier: "same_side", advice: sizedAdvice(167.18, "add", 300, 0),
    matches: [{ id: "k-eagles", tier: "same_side", inMath: true, note: null }],
  });
  assert.deepEqual(sizedLive(EAGLES_MINUS_105, spreadBet("k-bears", "home", 7.5)), {
    tier: "opposite", advice: sizedAdvice(507.93, "bet", 0, 300),
    matches: [{ id: "k-bears", tier: "opposite", inMath: true, note: null }],
  });
});

// ---- a held bet at a number the screen shows no fair for: bounded, not left out
//
// P(total > 51.5) is somewhere between 0 (past the ladder's end) and the
// row's own 38.0% at 47.5 (Unabated's edge); the stake is the smaller of
// the two ends. Checked by hand with condkelly at each end: $326.92 with the
// bet certainly lost, $125.95 with it winning exactly when the Over 47.5 does.

test("live sizing: a held Over 51.5 past the ladder's end is sized at the worst case its bounds allow, $125.95, and says so", () => {
  const held = heldBet({ id: "k-over-515", side: "over", points: 51.5 });
  assert.deepEqual(sizedLive(OVER_475, held), {
    tier: "same_side", advice: sizedAdvice(125.95, "add", 300, 0),
    matches: [{ id: "k-over-515", tier: "same_side", inMath: true, note: null, bound: "no fair at 51.5 \u00b7 bounded 0%\u201338%" }],
  });
  const row = sizedLiveRow(OVER_475, held);
  assert.deepEqual(betsView.relatedLines({ matches: row.bet.matches }).map((line) => [line.tag, line.inMath]), [["same side \u00b7 no fair at 51.5 \u00b7 bounded 0%\u201338%", true]]);
  assert.deepEqual(edgeRows.stakeRail(row), { text: "add $125.95", note: "$340.54 alone", atSize: false });
});

test("live sizing: a held Under 51.5 past the ladder's end reads its bound from its own side, 62–100%, and takes the smaller end, $354.16", () => {
  assert.deepEqual(sizedLive(OVER_475, heldBet({ id: "k-under-515", side: "under", points: 51.5 })), {
    tier: "opposite", advice: sizedAdvice(354.16, "bet", 0, 300),
    matches: [{ id: "k-under-515", tier: "opposite", inMath: true, note: null, bound: "no fair at 51.5 \u00b7 bounded 62%\u2013100%" }],
  });
});

test("live sizing: a held Over 42.5 between two rungs is bracketed by them (50–62%); the worst end is the 44.5 rung's, $167.18", () => {
  const sized = sizedLive(OVER_475, heldBet({ id: "k-over-425", side: "over", points: 42.5 }));
  assert.deepEqual([sized.advice.bet, sized.matches[0].bound], [167.18, "no fair at 42.5 \u00b7 bounded 50%\u201362%"]);
});

test("live sizing: two held Overs past the ladder share one bracket; the end that puts both on the row's own cut leaves nothing to add", () => {
  // Corners by hand: $313.30 (both lost), $112.33 (48.5 at 38%, 49.5 lost), $0 (both at 38%: $600 already on the Over 47.5's outcome);
  // the fourth corner has the fairs rising with the number and is not a ladder.
  const held = [...heldBet({ id: "k-485", side: "over", points: 48.5 }), ...heldBet({ id: "k-495", side: "over", points: 49.5 })];
  const sized = sizedLive(OVER_475, held);
  assert.deepEqual([sized.advice.kind, sized.advice.bet, sized.advice.held], ["sized", 0, 600]);
  assert.deepEqual(edgeRows.stakeRail(sizedLiveRow(OVER_475, held)), { text: "add $0.00", note: "$340.54 alone", atSize: true });
});

test("live sizing: past six bounded numbers the calc is declined, and a live row then shows no stake at all", () => {
  const records = [48.5, 49.5, 50.5, 51.5, 52.5, 53.5, 54.5].flatMap((points) => heldBet({ id: `k-${points}`, side: "over", points }));
  const sized = sizedLive(OVER_475, records);
  assert.deepEqual([sized.advice.kind, sized.advice.reason], ["declined", betsView.REASON_TOO_MANY_BOUNDS]);
  assert.ok(sized.matches.every((match) => match.inMath === false && match.note === betsView.REASON_TOO_MANY_BOUNDS));
  assert.deepEqual(edgeRows.stakeRail(sizedLiveRow(OVER_475, records)), { text: "\u2014", note: betsView.REASON_TOO_MANY_BOUNDS, atSize: false });
});

// A whole-number live line: the push between 45.5 and 46.5 has no rung of
// its own, so it is bracketed by the 44.5 and 47.5 rungs (0 to 12 points of
// push). Checked by hand: $88.88 with no push, $60.85 with the most.
test("live sizing: a whole-number live Over 46 prices its push off the rungs around it, worst case the most push, $60.85", () => {
  const withOver46 = halftimeLadders();
  withOver46.rows.push(novigRow(3, 0, 46, 44.5, 136, 150, 6));
  const OVER_46 = "live:125807:pt1:bt3:si0:ms89:alt46";
  const sized = sizedLive(OVER_46, heldBet({ id: "k-over", side: "over", points: 44.5 }), withOver46);
  assert.deepEqual([sized.advice.kind, sized.advice.bet, sized.advice.verb, sized.advice.held, sized.advice.alone], ["sized", 60.85, "add", 300, 300]);
  assert.deepEqual(sized.matches, [{ id: "k-over", tier: "same_side", inMath: true, note: null }]);
});

// ---- periods the game has finished

test("settledPeriods: the quarters and halves decided at each break; the full game never", () => {
  assert.deepEqual(live.settledPeriods(4), ["1Q"]);
  assert.deepEqual(live.settledPeriods(5), ["1Q", "2Q", "1H"]);
  assert.deepEqual(live.settledPeriods(2), ["1Q", "2Q", "1H"]);
  assert.deepEqual(live.settledPeriods(6), ["1Q", "2Q", "1H", "3Q"]);
  assert.deepEqual(live.settledPeriods(7), ["1Q", "2Q", "1H", "3Q", "4Q"]);
  assert.deepEqual(live.settledPeriods(3), ["1Q", "2Q", "1H", "3Q", "4Q", "2H"]);
  assert.deepEqual(live.settledPeriods(null), []);
});

test("live sizing: a held 1H Over at halftime is settled, not a position: named grey, the row sized alone", () => {
  assert.deepEqual(sizedLive(OVER_475, heldBet({ id: "k-1h-over", side: "over", points: 24.5, period: "1H" })), {
    tier: "related_same", advice: sizedAdvice(ALONE, "bet", 0, 0),
    matches: [{ id: "k-1h-over", tier: "related_same", inMath: false, note: "1H \u00b7 settled" }],
  });
});

// Across periods the corners of a bound are not the worst case: the
// worst-case pairing puts the smallest stake inside it (a 1H Over with no
// live fair, FG Over 47.5: $326.92 and $354.16 at the ends, $126.15 at a
// 38% chance, scanned by hand), so these decline instead.

test("live sizing: a held 1H Over at the end of Q1 with no 1H line on the screen declines the calc: no number, the reason named", () => {
  const endOfQ1 = { ...halftimeLadders(), games: [game({ eventPeriodTypeId: 4 })] };
  const held = heldBet({ id: "k-1h-over", side: "over", points: 24.5, period: "1H" });
  const reason = "no live fair at 24.5 for the 1H bet";
  assert.deepEqual(sizedLive(OVER_475, held, endOfQ1), {
    tier: "related_same", advice: { ...sizedAdvice(ALONE, "bet", 0, 0), kind: "declined", reason },
    matches: [{ id: "k-1h-over", tier: "related_same", inMath: false, note: reason }],
  });
  assert.deepEqual(edgeRows.stakeRail(sizedLiveRow(OVER_475, held, endOfQ1)), { text: "\u2014", note: reason, atSize: false });
});

// The end of Q1 with Novig's 1H total on the screen too: both sides of 24.5 at even money.
function endOfQ1WithFirstHalf() {
  const base = halftimeLadders();
  const firstHalf = (sideIndex, fair, price, edgePct) => row({
    key: `live:125807:pt2:bt3:si${sideIndex}:ms89`, periodTypeId: 2, betTypeId: 3, sideIndex, sideKey: `si${sideIndex}`, sideLabel: null,
    bookId: 89, bookName: "Novig", isAlt: false, mainPoints: null, points: 24.5, price, fair, edgePct, marketId: 3,
  });
  return { ...base, games: [game({ eventPeriodTypeId: 4 })], rows: [...base.rows, firstHalf(0, 100, -105, -2.4), firstHalf(1, 100, -105, -2.4)] };
}

// The worst-case pairing puts the 1H win on the FG rows the Over 47.5 wins
// on, so a 50% 1H Over cuts it exactly as a same-period Over 44.5 does.
test("live sizing: a held 1H Over at a number the screen does price is in the math across periods as pregame (worst-case pairing), $167.18", () => {
  const sized = sizedLive(OVER_475, heldBet({ id: "k-1h-over", side: "over", points: 24.5, period: "1H" }), endOfQ1WithFirstHalf());
  assert.deepEqual(sized, {
    tier: "related_same", advice: sizedAdvice(167.18, "add", 300, 0),
    matches: [{ id: "k-1h-over", tier: "related_same", inMath: true, note: null }],
  });
});

test("live sizing: a bounded fair in the row's period while a 1H bet is in the math declines the calc: the corners are not the worst case then", () => {
  const held = [...heldBet({ id: "k-1h-over", side: "over", points: 24.5, period: "1H" }), ...heldBet({ id: "k-over-515", side: "over", points: 51.5 })];
  const sized = sizedLive(OVER_475, held, endOfQ1WithFirstHalf());
  const reason = "no live fair at 51.5 with a 1H bet held";
  assert.deepEqual([sized.advice.kind, sized.advice.reason], ["declined", reason]);
  assert.deepEqual(sized.matches.map((match) => [match.id, match.inMath, match.note]), [["k-over-515", false, reason], ["k-1h-over", false, reason]]);
});

test("hasReadyFair: the game's fair is ready in a payload not yet dropped; expired, another game or an old payload is not", () => {
  const ready = payload();
  assert.equal(live.hasReadyFair([ready], 125807, NOW), true);
  assert.equal(live.hasReadyFair([ready], 7, NOW), false);
  // Given the production time a live ticket was captured at: only that fair counts.
  assert.equal(live.hasReadyFair([ready], 125807, NOW, "2026-09-28T00:52:10.000"), true);
  assert.equal(live.hasReadyFair([ready], 125807, NOW, "2026-09-28T00:40:00.000"), false);
  assert.equal(live.hasReadyFair([ready], 125807, NOW, null), true);
  assert.equal(live.hasReadyFair([payload({ games: [game({ fairStatus: "expired" })] })], 125807, NOW), false);
  assert.equal(live.hasReadyFair([payload({ at: NOW - live.LIVE_DROP_MS - 1 })], 125807, NOW), false);
  assert.equal(live.hasReadyFair([], 125807, NOW), false);
});
