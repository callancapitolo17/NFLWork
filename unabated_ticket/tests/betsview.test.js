// Run: node --test unabated_ticket/tests
// betsview.js: the panel's bet-history presentation helpers (#114).
const test = require("node:test");
const assert = require("node:assert/strict");
const view = require("../extension/betsview.js");
const fs = require("node:fs");
const path = require("node:path");
const teams = require("../extension/teams.js");
const betsLib = require("../extension/bets.js");
const kelly = require("../extension/kelly.js");
const ladderLib = require("../extension/ladder.js");
// The runtime team index the panel builds from Unabated's snapshots, from a
// captured copy (fixtures/teams_index.json) — keys are "<league>:<Unabated id>".
teams.loadIndex(JSON.parse(fs.readFileSync(path.join(__dirname, "fixtures", "teams_index.json"), "utf8")).leagues);
const key = (league, name) => teams.teamKey(league, name);

const NOW = Date.parse("2026-09-11T21:00:00Z");
const iso = (msAgo) => new Date(NOW - msAgo).toISOString();

function payloadWith(sources) {
  return { generatedAt: iso(0), sources, bets: [] };
}

function record(id, status, overrides) {
  return Object.assign({
    id, source: "kalshi_api", venue: "kalshi", league: "cfb", status, sourceFetchedAt: iso(0),
    awayTeam: "Chattanooga", homeTeam: "Eastern Kentucky", awayKey: null, homeKey: null,
    betType: "spread", period: "FG", side: "away", points: -5.5, price: 138, stake: 168, placedAt: iso(3600e3),
    closedAt: null, unmatchable: null, approx: [],
  }, overrides);
}

test("fmtAgeShort: seconds, minutes, hours, days; never negative", () => {
  assert.equal(view.fmtAgeShort(20e3), "20 s");
  assert.equal(view.fmtAgeShort(3 * 60e3), "3 min");
  assert.equal(view.fmtAgeShort(2 * 3600e3), "2 h");
  assert.equal(view.fmtAgeShort(3 * 86400e3), "3 d");
  assert.equal(view.fmtAgeShort(-5000), "0 s");
});

test("freshnessLevel: green under 5 min, amber under 60 min, red past that or never", () => {
  assert.equal(view.freshnessLevel(0), "green");
  assert.equal(view.freshnessLevel(5 * 60e3 - 1), "green");
  assert.equal(view.freshnessLevel(5 * 60e3), "amber");
  assert.equal(view.freshnessLevel(60 * 60e3 - 1), "amber");
  assert.equal(view.freshnessLevel(60 * 60e3), "red");
  assert.equal(view.freshnessLevel(null), "red");
});

test("sourceRows: one row per venue; unconfigured venues say so; a failed poll keeps its fetchedAt and shows the error", () => {
  const rows = view.sourceRows(payloadWith({
    kalshi: { fetchedAt: iso(20e3), ok: true, error: null, count: 14 },
    novig: { fetchedAt: iso(90 * 60e3), ok: false, error: "HTTP 401", count: 3 },
  }), NOW);
  assert.deepEqual(rows.map((row) => row.venue), ["kalshi", "betonline", "novig", "prophetx", "bfa", "wagerzon"]);
  const [kalshi, betonline, novig, prophetx, bfa, wagerzon] = rows;
  assert.equal(kalshi.configured, true);
  assert.equal(kalshi.level, "green");
  assert.equal(kalshi.ageText, "20 s");
  assert.equal(kalshi.count, 14);
  assert.equal(kalshi.error, null);
  assert.equal(betonline.configured, false);
  assert.equal(betonline.note, "no source configured");
  assert.equal(betonline.level, "none");
  assert.equal(novig.level, "red");
  assert.equal(novig.error, "HTTP 401");
  assert.equal(novig.count, 3);
  assert.equal(prophetx.configured, false);
  assert.equal(bfa.configured, false);
  assert.equal(wagerzon.configured, false);
});

test("sourceRows: a configured source that never succeeded is red with 'never'", () => {
  const [kalshi] = view.sourceRows(payloadWith({ kalshi: { fetchedAt: null, ok: false, error: "boom", count: 0 } }), NOW);
  assert.equal(kalshi.level, "red");
  assert.equal(kalshi.ageText, "never");
  assert.equal(kalshi.error, "boom");
});

test("sourceRows: no payload at all is six unconfigured rows", () => {
  assert.equal(view.sourceRows(null, NOW).filter((row) => row.configured).length, 0);
});

test("serviceStatus: not reached, unreachable since, reached", () => {
  assert.deepEqual(view.serviceStatus(null, NOW), { unreachable: true, text: "bets service not reached yet" });
  const down = view.serviceStatus({ okAt: NOW - 600e3, error: "Failed to fetch", errorAt: NOW - 10e3, unreachableSince: NOW - 180e3 }, NOW);
  assert.equal(down.unreachable, true);
  assert.match(down.text, /^bets service unreachable since \d{1,2}:\d{2}(?: [AP]M)? \(3 min\): Failed to fetch$/);
  const up = view.serviceStatus({ okAt: NOW - 20e3, error: null, errorAt: null, unreachableSince: null }, NOW);
  assert.deepEqual(up, { unreachable: false, text: "bets service reached 20 s ago" });
});

test("sourcesUnavailable: true with no sources or every source red; false while one is amber", () => {
  assert.equal(view.sourcesUnavailable(null, NOW), true);
  assert.equal(view.sourcesUnavailable(payloadWith({}), NOW), true);
  assert.equal(view.sourcesUnavailable(payloadWith({ kalshi: { fetchedAt: iso(2 * 3600e3), ok: true, error: null, count: 1 } }), NOW), true);
  assert.equal(view.sourcesUnavailable(payloadWith({ kalshi: { fetchedAt: iso(30 * 60e3), ok: true, error: null, count: 1 } }), NOW), false);
  assert.equal(view.sourcesUnavailable(payloadWith({
    kalshi: { fetchedAt: iso(2 * 3600e3), ok: true, error: null, count: 1 },
    novig: { fetchedAt: iso(10e3), ok: true, error: null, count: 1 },
  }), NOW), false);
});

test("headerLine: open count, then every venue with its age or a dash", () => {
  const records = [record("a", "open"), record("b", "open"), record("c", "won"), record("d", "closed")];
  const payload = payloadWith({ kalshi: { fetchedAt: iso(20e3), ok: true, error: null, count: 4 } });
  assert.equal(view.headerLine(records, payload, NOW), "bets: 2 open · kalshi 20 s · betonline — · novig — · prophetx — · bfa — · wagerzon —");
  assert.equal(view.headerLine([], null, NOW), "bets: 0 open · kalshi — · betonline — · novig — · prophetx — · bfa — · wagerzon —");
});

test("bannerLines: at most five, strongest first as given, and the count of the rest", () => {
  const matches = Array.from({ length: 7 }, (_, i) => ({ tier: "same_game", label: `m${i}` }));
  const { shown, more } = view.bannerLines(matches);
  assert.deepEqual(shown.map((m) => m.label), ["m0", "m1", "m2", "m3", "m4"]);
  assert.equal(more, 2);
  assert.deepEqual(view.bannerLines(matches.slice(0, 3)), { shown: matches.slice(0, 3), more: 0 });
  assert.deepEqual(view.bannerLines([]), { shown: [], more: 0 });
});

// ---- #130: stake advice (conditional Kelly against held bets) ---------------------
//
// The cards below were measured live 2026-09-17/18 on a Kelly bankroll of
// $8,000 (32000 x 0.25). Bets go through the real matcher, so tiers, positions
// and the math are exercised together.

const SIZING = { bankroll: 32000, multiplier: 0.25 };
const KICKOFF = "2026-09-20T17:00:00.000Z";

function nflLine(overrides) {
  return Object.assign({
    league: "nfl", eventId: 700001, awayTeam: "Detroit Lions", homeTeam: "Buffalo Bills",
    eventStart: KICKOFF, eventStartMs: Date.parse(KICKOFF), betType: "Total", period: "FG", sideIndex: 0, points: 61.5, rotation: null,
  }, overrides);
}

function heldRecord(id, overrides) {
  const base = record(id, "open", Object.assign({
    league: "nfl", awayTeam: "Detroit Lions", homeTeam: "Buffalo Bills", eventStart: KICKOFF, isParlayLeg: false,
  }, overrides));
  const netOdds = base.price > 0 ? base.price / 100 : 100 / Math.abs(base.price);
  if (!("toWin" in overrides)) base.toWin = base.stake == null ? null : base.stake * netOdds;
  return betsLib.resolveTeamKeys([base])[0];
}

// A ladder as ladder.buildLadder returns it; neighbours differ so no rung is flat.
function ladderStub(byPeriodAxis) {
  return (period, axis) => {
    const rungs = byPeriodAxis[`${period} ${axis}`];
    return rungs ? { axis, rungs } : null;
  };
}

function adviceFor(line, price, edgePct, records, ladderOf) {
  const { matches } = betsLib.matchBets(line, records);
  return view.stakeAdvice({ line, price, edgePct, ...SIZING, matches, ladderOf });
}

const LIONS_BILLS_HELD = [
  heldRecord("over-61.5", { betType: "total", side: "over", points: 61.5, price: 212, stake: 270 }),
  heldRecord("under-51.5-a", { betType: "total", side: "under", points: 51.5, price: 127, stake: 352 }),
  heldRecord("under-51.5-b", { betType: "total", side: "under", points: 51.5, price: 150, stake: 61 }),
];
const LIONS_BILLS_LADDER = ladderStub({ "FG total": [[50.5, 0.61], [51.5, 0.5745], [52.5, 0.54]] });

test("bets.js and ladder.js name the two market axes the same: a position's axis is the ladder's key", () => {
  assert.equal(betsLib.AXIS_TOTAL, ladderLib.AXIS_TOTAL);
  assert.equal(betsLib.AXIS_MARGIN, ladderLib.AXIS_MARGIN);
});

test("stakeAdvice: the real ladder feeds it — fairs read off feed lines, a rung the feed lacks is no fair", () => {
  const over = (points, bacr) => ({ eventId: 700001, leagueId: 1, periodTypeId: 1, betTypeId: 3, sideIndex: 0, points, bacr, fromSnapshot: true });
  const eventLines = [over(50.5, -156), over(51.5, -135), over(52.5, -117)];
  const ladderOf = (period, axis) => ladderLib.buildLadder(eventLines, { periodTypeId: ladderLib.periodTypeIdOf(period), axis });
  const sized = adviceFor(nflLine(), 213, 7.19, LIONS_BILLS_HELD, ladderOf);
  assert.equal(sized.kind, "sized");
  // -135 is 57.447%, the whole American price behind the card's 57.45%: two cents off its $188.32.
  assert.equal(sized.bet, 188.34);
  const noRung = adviceFor(nflLine(), 213, 7.19, LIONS_BILLS_HELD, (period, axis) => ladderLib.buildLadder(eventLines.slice(0, 1), { periodTypeId: 1, axis }));
  assert.deepEqual(noRung.matches.filter((match) => !match.inMath).map((match) => match.note), ["no fair at 51.5", "no fair at 51.5"]);
  assert.deepEqual([noRung.held, noRung.against], [270, 0]);
});

test("stakeAdvice: nothing held is the standalone Kelly stake, exactly", () => {
  const advice = adviceFor(nflLine(), 213, 7.19, [], LIONS_BILLS_LADDER);
  const standalone = kelly.kellyStakeFromEdge({ bookPrice: 213, edgePct: 7.19, ...SIZING }).stake;
  assert.deepEqual(advice, { kind: "none", bet: standalone, alone: standalone, verb: "bet", held: 0, against: 0, reason: null, matches: [], cappedAt: null });
  assert.equal(view.stakeAdviceWords(advice), null);
  assert.equal(view.suggestedBetAmount(advice), standalone);
});

test("stakeAdvice: Lions @ Bills — the same over held and two unders at another number, add $188.32", () => {
  const advice = adviceFor(nflLine(), 213, 7.19, LIONS_BILLS_HELD, LIONS_BILLS_LADDER);
  assert.equal(advice.kind, "sized");
  assert.equal(advice.bet, 188.32);
  assert.equal(advice.verb, "add");
  assert.equal(advice.held, 270);
  assert.equal(advice.against, 413);
  assert.ok(advice.matches.every((match) => match.inMath && match.note === null));
  assert.deepEqual(view.badges({ tier: "same_line", matches: advice.matches, advice }),
    [{ kind: "held", text: "held $270" }, { kind: "against", text: "against $413" }]);
  assert.deepEqual(view.stakeAdviceWords(advice), { verb: "add", bet: "$188.32", alone: "$270.05 alone", cap: null });
  assert.equal(view.stakeAdviceLine(advice), "add $188.32, $270.05 alone");
  assert.equal(view.suggestedBetAmount(advice), 188.32);
});

test("stakeAdvice: a held spread never changes the size of a total — shown grey as game · not sized (#129)", () => {
  const billsSpread = heldRecord("bills-8.5", { betType: "spread", side: "home", points: -8.5, price: 212, stake: 472 });
  const advice = adviceFor(nflLine(), 213, 7.19, LIONS_BILLS_HELD.concat([billsSpread]), LIONS_BILLS_LADDER);
  assert.equal(advice.bet, 188.32);
  assert.equal(advice.held, 270);
  assert.equal(advice.against, 413);
  const last = advice.matches[advice.matches.length - 1];
  assert.equal(last.bet.id, "bills-8.5");
  assert.equal(last.inMath, false);
  const lines = view.relatedLines({ matches: advice.matches });
  assert.deepEqual(lines.map((line) => [line.inMath, line.tag]),
    [[true, "this line"], [true, "other side"], [true, "other side"], [false, "game · not sized"]]);
});

test("stakeAdvice: Chargers moneyline with the Chargers +3.5 held — in the math as the same side, add $0", () => {
  const line = nflLine({ awayTeam: "Los Angeles Chargers", betType: "Moneyline", sideIndex: 0, points: null });
  const held = heldRecord("chargers+3.5", { awayTeam: "Los Angeles Chargers", betType: "spread", side: "away", points: 3.5, price: 122, stake: 400 });
  const advice = adviceFor(line, 213, 5.03, [held], ladderStub({ "FG margin": [[-6.5, 0.56], [-3.5, 0.4785], [-0.5, 0.34]] }));
  assert.equal(advice.kind, "sized");
  assert.equal(advice.bet, 0);
  assert.equal(advice.verb, "add");
  assert.equal(advice.matches[0].tier, "same_side");
  assert.deepEqual(view.badges({ tier: "same_side", matches: advice.matches, advice }), [{ kind: "held", text: "held $400" }]);
  assert.deepEqual(view.stakeAdviceWords(advice), { verb: "add", bet: "$0", alone: "$188.92 alone", cap: null });
  assert.equal(view.suggestedBetAmount(advice), 0);
});

test("stakeAdvice: a Novig line with $17 resting says add $17, never the $32.58 Kelly wants", () => {
  // Live 2026-09-22: $35 held on Portland Fire, the rail said "add $32.58" on a line with $17 behind it.
  const advice = view.stakeAdvice({ line: nflLine(), price: 213, edgePct: 7.19, ...SIZING, matches: betsLib.matchBets(nflLine(), LIONS_BILLS_HELD).matches, ladderOf: LIONS_BILLS_LADDER, liquidity: 17 });
  assert.equal(advice.kind, "sized");
  assert.equal(advice.bet, 17);
  assert.equal(advice.cappedAt, 17);
  assert.equal(view.suggestedBetAmount(advice), 17);
  assert.deepEqual(view.stakeAdviceWords(advice), { verb: "add", bet: "$17", alone: "$270.05 alone", cap: "all $17 liq" });
  assert.equal(view.stakeAdviceLine(advice), "add $17, all $17 liq, $270.05 alone");
});

test("stakeAdvice: nothing held and thin liquidity still caps, and says so", () => {
  const advice = view.stakeAdvice({ line: nflLine(), price: 213, edgePct: 7.19, ...SIZING, matches: [], ladderOf: LIONS_BILLS_LADDER, liquidity: 50 });
  assert.equal(advice.kind, "none");
  assert.equal(advice.bet, 50);
  assert.deepEqual(view.stakeAdviceWords(advice), { verb: "bet", bet: "$50", alone: null, cap: "all $50 liq" });
});

test("capAtLiquidity: deep or unreported liquidity leaves the stake alone", () => {
  assert.deepEqual(view.capAtLiquidity(32.58, 17), { stake: 17, cappedAt: 17 });
  assert.deepEqual(view.capAtLiquidity(32.58, 5000), { stake: 32.58, cappedAt: null });
  assert.deepEqual(view.capAtLiquidity(32.58, null), { stake: 32.58, cappedAt: null });
  assert.deepEqual(view.capAtLiquidity(32.58, undefined), { stake: 32.58, cappedAt: null });
  assert.deepEqual(view.capAtLiquidity(null, 17), { stake: null, cappedAt: null });
  // Nothing resting: the stake is $0 and says why (the Ticket labels it "Nothing resting at this price").
  assert.deepEqual(view.capAtLiquidity(32.58, 0), { stake: 0, cappedAt: 0 });
});

test("stakeAdvice: a 1H over held sizes the FG over on the worst case, $408.39 not $666.67", () => {
  const line = nflLine({ points: 48.5 });
  const held = heldRecord("1h-over", { betType: "total", period: "1H", side: "over", points: 24.5, price: 110, stake: 300 });
  const advice = adviceFor(line, 120, 10, [held], ladderStub({ "1H total": [[23.5, 0.6], [24.5, 0.55], [25.5, 0.5]] }));
  assert.equal(advice.matches[0].tier, "related_same");
  assert.equal(advice.bet, 408.39);
  assert.equal(advice.held, 300);
  assert.equal(view.stakeAdviceLine(advice), "add $408.39, $666.67 alone");
});

test("stakeAdvice: the OTHER direction in another period is left out — no hedge credit, no worst-case penalty", () => {
  const line = nflLine({ points: 48.5 });
  const held = heldRecord("1h-under", { betType: "total", period: "1H", side: "under", points: 24.5, price: 110, stake: 300 });
  const advice = adviceFor(line, 120, 10, [held], ladderStub({ "1H total": [[23.5, 0.6], [24.5, 0.55], [25.5, 0.5]] }));
  const standalone = kelly.kellyStakeFromEdge({ bookPrice: 120, edgePct: 10, ...SIZING }).stake;
  assert.equal(advice.matches[0].tier, "related_opposite");
  // The worst-case pairing would have said $407.22 here; the measured 1H/FG link wants about $800.
  assert.deepEqual([advice.kind, advice.bet, advice.against, advice.verb], ["none", standalone, 0, "bet"]);
  assert.deepEqual([advice.matches[0].inMath, advice.matches[0].note], [false, "other period \u00b7 not sized"]);
  assert.deepEqual(view.badges({ tier: "related_opposite", matches: advice.matches, advice }), [{ kind: "against", text: "against" }]);
  // Beside a same-direction 1H bet, only that one sizes the row.
  const over = heldRecord("1h-over", { betType: "total", period: "1H", side: "over", points: 24.5, price: 110, stake: 300 });
  const both = adviceFor(line, 120, 10, [over, held], ladderStub({ "1H total": [[23.5, 0.6], [24.5, 0.55], [25.5, 0.5]] }));
  assert.deepEqual([both.bet, both.held, both.against], [408.39, 300, 0]);
});

test("stakeAdvice: the same line held needs no ladder at all", () => {
  const line = nflLine({ sideIndex: 1, points: 52.5 });
  const held = heldRecord("under-52.5", { betType: "total", side: "under", points: 52.5, price: 125, stake: 181.5 });
  const advice = adviceFor(line, 125, 4.65, [held], () => null);
  assert.equal(advice.bet, 116.1);
  const atSize = adviceFor(line, 125, 4.65, [Object.assign({}, held, { stake: 400, toWin: 500 })], () => null);
  assert.equal(atSize.bet, 0);
});

test("stakeAdvice guards: a bet with no fair, a flat rung, no stake or a parlay leg is left out and named", () => {
  const line = nflLine();
  const under = (id, overrides) => heldRecord(id, Object.assign({ betType: "total", side: "under", points: 36.5, price: 300, stake: 50 }, overrides));
  const standalone = kelly.kellyStakeFromEdge({ bookPrice: 213, edgePct: 7.19, ...SIZING }).stake;
  const cases = [
    [under("no-rung", {}), LIONS_BILLS_LADDER, "no fair at 36.5"],
    [under("flat", {}), ladderStub({ "FG total": [[35.5, 0.802], [36.5, 0.802], [37.5, 0.802]] }), "no fair at 36.5"],
    [under("leg", { isParlayLeg: true }), LIONS_BILLS_LADDER, "parlay leg"],
    [under("no-stake", { stake: null }), LIONS_BILLS_LADDER, "no stake on the record"],
  ];
  for (const [held, ladderOf, note] of cases) {
    const advice = adviceFor(line, 213, 7.19, [held], ladderOf);
    assert.equal(advice.kind, "none", held.id);
    assert.equal(advice.bet, standalone, held.id);
    assert.equal(advice.against, 0, held.id);
    assert.deepEqual([advice.matches[0].inMath, advice.matches[0].note], [false, note], held.id);
    // Still on the other side of this line: flagged, bare, with no dollars.
    assert.deepEqual(view.badges({ tier: advice.matches[0].tier, matches: advice.matches, advice }), [{ kind: "against", text: "against" }], held.id);
    assert.deepEqual(view.relatedLines({ matches: advice.matches })[0], { tier: advice.matches[0].tier, inMath: false, tag: note, text: advice.matches[0].label }, held.id);
  }
});

test("stakeAdvice guards: a period Unabated has no ladder for (Kalshi's F5) is no fair", () => {
  const held = heldRecord("f5-over", { betType: "total", period: "F5", side: "over", points: 36.5, price: -300, stake: 50 });
  const advice = adviceFor(nflLine(), 213, 7.19, [held], LIONS_BILLS_LADDER);
  assert.deepEqual([advice.kind, advice.held, advice.matches[0].inMath, advice.matches[0].note], ["none", 0, false, "no fair at 36.5"]);
  assert.deepEqual(view.badges({ tier: advice.matches[0].tier, matches: advice.matches, advice }), [{ kind: "held", text: "held" }]);
});

test("stakeAdvice: a ladder that is not monotone declines the whole calc and shows the standalone stake", () => {
  // The over's own chance at 61.5 is 34.2%; a 25% fair for Over 51.5 cannot sit below it.
  const advice = adviceFor(nflLine(), 213, 7.19, LIONS_BILLS_HELD, ladderStub({ "FG total": [[50.5, 0.3], [51.5, 0.25], [52.5, 0.2]] }));
  const standalone = kelly.kellyStakeFromEdge({ bookPrice: 213, edgePct: 7.19, ...SIZING }).stake;
  assert.equal(advice.kind, "declined");
  assert.equal(advice.reason, "ladder not monotone");
  assert.equal(advice.bet, standalone);
  assert.deepEqual([advice.held, advice.against, advice.verb], [0, 0, "bet"]);
  assert.ok(advice.matches.every((match) => !match.inMath && match.note === "ladder not monotone"));
  assert.equal(view.stakeAdviceWords(advice), null);
});

test("stakeAdvice: an edge that implies a certain win declines instead of throwing in the render loop", () => {
  const line = nflLine({ betType: "Moneyline", sideIndex: 1, points: null });
  const held = heldRecord("bills-ml", { betType: "moneyline", side: "home", points: null, price: -4000, stake: 500 });
  const advice = adviceFor(line, -5000, 2.5, [held], () => null);
  assert.deepEqual([advice.kind, advice.reason], ["declined", "edge implies a certain win"]);
  assert.equal(advice.bet, advice.alone);
});

test("stakeAdvice: a whole-number row needs the two rungs around it to price its push", () => {
  const line = nflLine({ points: 61 });
  const declined = adviceFor(line, 213, 7.19, LIONS_BILLS_HELD.slice(1), LIONS_BILLS_LADDER);
  assert.equal(declined.kind, "declined");
  assert.equal(declined.reason, "no fair at 60.5 to price the push");
  const ladderOf = ladderStub({ "FG total": [[50.5, 0.61], [51.5, 0.5745], [52.5, 0.54], [60.5, 0.37], [61.5, 0.33], [62.5, 0.3]] });
  assert.equal(adviceFor(line, 213, 7.19, LIONS_BILLS_HELD.slice(1), ladderOf).kind, "sized");
});

test("stakeAdvice: a price with no edge of its own is $0 even against a held bet (hedge sizing is deferred)", () => {
  const line = nflLine({ points: 52.5 });
  const under = heldRecord("under-52.5", { betType: "total", side: "under", points: 52.5, price: 125, stake: 600 });
  const advice = adviceFor(line, -115, 0, [under], () => null);
  assert.deepEqual([advice.kind, advice.bet, advice.against], ["none", 0, 600]);
  const unknownEdge = adviceFor(line, -115, null, [under], () => null);
  assert.equal(unknownEdge.bet, null);
  assert.equal(view.suggestedBetAmount(unknownEdge), 0);
  assert.equal(view.suggestedBetAmount(null), 0);
});

test("badges: a plain game marker for another market, nothing with no match", () => {
  const spread = heldRecord("bills-8.5", { betType: "spread", side: "home", points: -8.5, price: 212, stake: 472 });
  const advice = adviceFor(nflLine(), 213, 7.19, [spread], LIONS_BILLS_LADDER);
  assert.deepEqual(view.badges({ tier: "same_game", matches: advice.matches, advice }), [{ kind: "game", text: "game" }]);
  assert.deepEqual(view.badges({ tier: null, matches: [], advice }), []);
  assert.deepEqual(view.badges(null), []);
});

test("stakeAdviceWords: no `alone` line when the held bets left the number where it was", () => {
  assert.deepEqual(view.stakeAdviceWords({ kind: "sized", verb: "bet", bet: 183.04, alone: 183.0412 }), { verb: "bet", bet: "$183.04", alone: null, cap: null });
  assert.equal(view.stakeAdviceLine({ kind: "sized", verb: "bet", bet: 183.04, alone: 183.0412 }), "bet $183.04");
  assert.equal(view.stakeAdviceLine({ kind: "none", bet: 183.04, alone: 183.04 }), null);
});

test("relatedLines: bets in the math carry their tier's tag, the rest carry why they are not", () => {
  const flag = {
    tier: "same_line",
    matches: [
      { tier: "same_line", label: "Chattanooga -5.5 +138 · 42.0¢ · $300 · Kalshi", inMath: true, note: null },
      { tier: "related_opposite", label: "1H Eastern Kentucky +2.5 -110 · 52.4¢ · $200 · Kalshi", inMath: true, note: null },
      { tier: "same_game", label: "1H Under 22.5 -104 · 51.0¢ · $255 · Kalshi", inMath: false, note: "game · not sized" },
    ],
  };
  assert.deepEqual(view.relatedLines(flag), [
    { tier: "same_line", inMath: true, tag: "this line", text: "Chattanooga -5.5 +138 · 42.0¢ · $300 · Kalshi" },
    { tier: "related_opposite", inMath: true, tag: "other side", text: "1H Eastern Kentucky +2.5 -110 · 52.4¢ · $200 · Kalshi" },
    { tier: "same_game", inMath: false, tag: "game · not sized", text: "1H Under 22.5 -104 · 51.0¢ · $255 · Kalshi" },
  ]);
  assert.deepEqual(view.relatedLines(null), []);
});

test("mergeServicePayload: fresh wins on the same id, team keys are filled, old settled bets are pruned", () => {
  const stored = [
    record("kalshi:x:yes", "open", { sourceFetchedAt: iso(60e3), stake: 100 }),
    record("kalshi:old:yes", "won", { sourceFetchedAt: iso(60e3), closedAt: iso(40 * 86400e3) }),
    record("kalshi:gone:yes", "open", { sourceFetchedAt: iso(60e3) }),
  ];
  const payload = { generatedAt: iso(0), sources: {}, bets: [record("kalshi:x:yes", "open", { stake: 168 })] };
  const merged = view.mergeServicePayload(stored, payload, NOW);
  const ids = merged.map((r) => r.id).sort();
  assert.deepEqual(ids, ["kalshi:gone:yes", "kalshi:x:yes"]);
  const x = merged.find((r) => r.id === "kalshi:x:yes");
  assert.equal(x.stake, 168);
  assert.equal(x.awayKey, key("cfb", "Chattanooga"));
  assert.equal(x.homeKey, key("cfb", "Eastern Kentucky"));
});

test("mergeServicePayload: the payload's team crosswalk keys a record before its names do (#118 step 4)", () => {
  const stored = [record("kalshi:x:yes", "open", { sourceFetchedAt: iso(60e3), awayTeam: "Chatt (venue spelling)", awayKey: null, homeKey: null })];
  const payload = {
    generatedAt: iso(0), sources: {}, bets: [record("kalshi:x:yes", "open", { awayTeam: "Chatt (venue spelling)", awayKey: null, homeKey: null })],
    crosswalk: [{ venue: "kalshi", league: "cfb", venueTeamKey: "Chatt (venue spelling)", unabatedTeamId: "41", learnedAt: iso(0) }],
  };
  const [x] = view.mergeServicePayload(stored, payload, NOW);
  assert.equal(x.awayKey, "cfb:41");
  assert.equal(x.homeKey, key("cfb", "Eastern Kentucky"));
  assert.deepEqual(view.crosswalkOf(payload), payload.crosswalk);
  assert.deepEqual(view.crosswalkOf({ bets: [] }), []);
  const [plain] = view.mergeServicePayload(stored, { ...payload, crosswalk: undefined }, NOW);
  assert.equal(plain.awayKey, null);
});

test("crosswalkRows: one Bets-tab line per served row — venue spelling, Unabated name, venue, league, when", () => {
  const rows = view.crosswalkRows([
    { venue: "novig", league: "cfb", venueTeamKey: "nv-1", venueTeamName: "Wazzu", unabatedTeamId: "717", unabatedTeamName: "Washington State", learnedFrom: "novig:o on board event 1", learnedAt: "2026-09-15T19:00:00Z" },
    { venue: "kalshi", league: "nfl", venueTeamKey: "PIT Steelers", venueTeamName: null, unabatedTeamId: "3", unabatedTeamName: null, learnedFrom: null, learnedAt: null },
    null, "junk",
  ]);
  assert.deepEqual(rows, [
    { what: "Wazzu → Washington State", meta: "Novig · CFB · learned Sep 15 3:00 PM", title: "from novig:o on board event 1" },
    { what: "PIT Steelers → team 3", meta: "Kalshi · NFL", title: "" },
  ]);
  assert.deepEqual(view.crosswalkRows(undefined), []);
});

test("mergeServicePayload: a venue's successful pull drops its stored records the payload no longer lists; a failed or absent pull keeps them", () => {
  const stored = [
    record("kalshi:gone:yes", "open", { sourceFetchedAt: iso(60e3) }),
    record("kalshi:kept:yes", "open", { sourceFetchedAt: iso(60e3) }),
    record("novig:gone:yes", "open", { venue: "novig", source: "novig_api", sourceFetchedAt: iso(60e3) }),
  ];
  const kalshiOk = { fetchedAt: iso(0), ok: true, error: null, count: 1 };
  const novigFailed = { fetchedAt: iso(3600e3), ok: false, error: "HTTP 503", count: 0 };
  const payload = { generatedAt: iso(0), sources: { kalshi: kalshiOk, novig: novigFailed }, bets: [record("kalshi:kept:yes", "open")] };
  assert.deepEqual(view.mergeServicePayload(stored, payload, NOW).map((r) => r.id).sort(), ["kalshi:kept:yes", "novig:gone:yes"]);
  // The same payload with kalshi's pull failed keeps the missing record.
  const failed = { ...payload, sources: { kalshi: { ...kalshiOk, ok: false, error: "HTTP 503" } } };
  assert.deepEqual(view.mergeServicePayload(stored, failed, NOW).map((r) => r.id).sort(), ["kalshi:gone:yes", "kalshi:kept:yes", "novig:gone:yes"]);
  assert.deepEqual([...view.venuesWithFreshPull(payload)], ["kalshi"]);
  assert.deepEqual([...view.venuesWithFreshPull({ sources: null })], []);
});

test("ticketAsLine: naive-UTC eventStart parsed, period defaults to FG, fields the matcher reads", () => {
  const line = view.ticketAsLine({
    league: "cfb", eventId: 900001, eventStart: "2026-09-12T20:00:00", awayTeam: "Chattanooga", homeTeam: "Eastern Kentucky",
    betType: "Spread", sideIndex: 0, points: -5.5, rotation: 301,
  });
  assert.equal(line.eventStartMs, Date.parse("2026-09-12T20:00:00Z"));
  assert.equal(line.period, "FG");
  assert.equal(line.betType, "Spread");
  assert.equal(line.rotation, 301);
  assert.equal(view.ticketAsLine({ league: "nfl", eventStart: null, betType: "Total", sideIndex: 1, points: 44.5, period: "1H" }).period, "1H");
  assert.equal(view.ticketAsLine({ league: "nfl", eventStart: null, betType: "Total", sideIndex: 1, points: 44.5 }).eventStartMs, null);
});

test("sanitizeBetsSettings: defaults, a trailing slash trimmed, junk ignored", () => {
  assert.deepEqual(view.sanitizeBetsSettings(null), { serviceUrl: "http://127.0.0.1:8094" });
  assert.deepEqual(view.sanitizeBetsSettings({ serviceUrl: "http://localhost:9000/", hideBet: true }), { serviceUrl: "http://localhost:9000" });
  assert.deepEqual(view.sanitizeBetsSettings({ serviceUrl: "not a url" }), { serviceUrl: "http://127.0.0.1:8094" });
});
