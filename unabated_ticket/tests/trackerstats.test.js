// Run: node --test unabated_ticket/tests
// server/tracker/trackerstats.js: records -> tickets (P&L, Pacific settle
// day, parlays collapsed to one ticket, fair-at-fill edge) and the summaries
// the Bet Tracker's Overview and Analysis pages show.
const test = require("node:test");
const assert = require("node:assert/strict");
const stats = require("../server/tracker/trackerstats.js");

function straight(overrides) {
  return Object.assign({
    id: "novig:1", venue: "novig", league: "nfl", betType: "spread", period: "FG", side: "home", points: -3.5,
    awayTeam: "Jets", homeTeam: "Bills", price: -110, stake: 110, toWin: 100, status: "won",
    placedAt: "2026-10-04T15:00:00Z", closedAt: "2026-10-04T20:30:00Z", eventStart: "2026-10-04T17:00:00Z",
    isParlayLeg: false, parlayId: null, raw: {},
  }, overrides);
}

test("P&L by status: won pays toWin, lost costs the stake, push and void are flat", () => {
  const tickets = stats.buildTickets([
    straight({ id: "a", status: "won" }), straight({ id: "b", status: "lost" }),
    straight({ id: "c", status: "push" }), straight({ id: "d", status: "void" }),
  ], []);
  const byId = Object.fromEntries(tickets.map((t) => [t.id, t.pnl]));
  assert.deepEqual(byId, { a: 100, b: -110, c: 0, d: 0 });
});

test("open bets split into live (started), upcoming, and no start time", () => {
  const now = Date.parse("2026-10-04T18:00:00Z");
  const tickets = stats.buildTickets([
    straight({ id: "later", status: "open", closedAt: null, eventStart: "2026-10-04T20:25:00Z" }),
    straight({ id: "early", status: "open", closedAt: null, eventStart: "2026-10-04T17:00:00Z" }),
    straight({ id: "kickoff", status: "open", closedAt: null, eventStart: "2026-10-04T18:00:00Z" }),
    straight({ id: "future", status: "open", closedAt: null, eventStart: null }),
  ], []);
  const { live, upcoming, noStart } = stats.splitOpenByStart(tickets, now);
  assert.deepEqual(live.map((t) => t.id), ["early", "kickoff"]);
  assert.deepEqual(upcoming.map((t) => t.id), ["later"]);
  assert.deepEqual(noStart.map((t) => t.id), ["future"]);
});

test("open, closed-early and unknown bets carry no P&L and stay out of the summary", () => {
  const tickets = stats.buildTickets([
    straight({ id: "open", status: "open", closedAt: null }), straight({ id: "sold", status: "closed" }),
    straight({ id: "gone", status: "unknown" }), straight({ id: "won" }),
  ], []);
  assert.equal(tickets.filter((t) => t.pnl !== null).length, 1);
  const total = stats.summarize(tickets);
  assert.equal(total.bets, 1);
  assert.equal(total.pnl, 100);
  assert.deepEqual(stats.exclusions(tickets), { open: 1, noResult: 2 });
});

test("the venue's own pnl wins over stake/toWin, a sold position with one counts, its merged side is skipped", () => {
  const tickets = stats.buildTickets([
    straight({ id: "kalshi:R:yes", venue: "kalshi", stake: 236, toWin: 10964, status: "won", pnl: 5449.46 }),
    straight({ id: "kalshi:R:no", venue: "kalshi", status: "closed", stake: 0, toWin: 0, mergedInto: "kalshi:R:yes" }),
    straight({ id: "kalshi:P:no", venue: "kalshi", status: "closed", stake: 194.01, pnl: 24.84 }),
    straight({ id: "kalshi:M:no", venue: "kalshi", status: "closed", stake: 732.54, pnl: -632.67 }),
  ], []);
  assert.deepEqual(tickets.map((t) => t.id).sort(), ["kalshi:M:no", "kalshi:P:no", "kalshi:R:yes"]);
  const total = stats.summarize(tickets);
  assert.equal(total.bets, 3);
  assert.equal(total.wins, 2);
  assert.equal(total.losses, 1);
  assert.equal(Math.round(total.pnl * 100) / 100, 4841.63);
  assert.deepEqual(stats.exclusions(tickets), { open: 0, noResult: 0 });
});

test("a won bet with no toWin is paid off its American price", () => {
  const [ticket] = stats.buildTickets([straight({ price: 150, stake: 100, toWin: null })], []);
  assert.equal(ticket.pnl, 150);
});

test("the settle day is the Pacific calendar day, not UTC", () => {
  // 02:30 UTC on Oct 5 is 19:30 PDT on Oct 4.
  const [ticket] = stats.buildTickets([straight({ closedAt: "2026-10-05T02:30:00Z" })], []);
  assert.equal(ticket.settledDay, "2026-10-04");
  assert.equal(ticket.weekday, "Sun");
  // Standard time: 07:59 UTC on Jan 15 is 23:59 PST on Jan 14.
  assert.equal(stats.pacificDay(Date.parse("2026-01-15T07:59:00Z")), "2026-01-14");
});

test("a parlay's legs collapse to one ticket on the ticket's stake; a teaser is named by its text", () => {
  const leg = (index, raw) => straight({
    id: "bfa:9:leg" + index, venue: "bfa", isParlayLeg: true, parlayId: "bfa:9", legIndex: index, legCount: 2,
    stake: 200, toWin: 600, status: "lost", raw,
  });
  const [parlay] = stats.buildTickets([leg(0, { headerDescription: "PARLAY (2 TEAMS)", parlayPrice: 264 }), leg(1, {})], []);
  assert.equal(parlay.kind, "Parlay");
  assert.equal(parlay.pnl, -200);
  assert.equal(parlay.displayPrice, 264);
  assert.equal(parlay.event, "2-leg parlay");
  const [teaser] = stats.buildTickets([leg(0, { headerDescription: "2 TEAM TEASERS" }), leg(1, {})], []);
  assert.equal(teaser.kind, "Teaser");
  assert.equal(teaser.fairProb, null);
});

test("a parlay settles with its last leg: Sunday and Monday legs land on Monday", () => {
  const leg = (index, closedAt) => straight({
    id: "wagerzon:5:leg" + index, venue: "wagerzon", isParlayLeg: true, parlayId: "wagerzon:5", legIndex: index,
    legCount: 2, stake: 100, toWin: 260, status: "won", closedAt,
  });
  // Leg 0 kicks off Sunday 10:00 PDT, leg 1 Monday 17:15 PDT.
  const [parlay] = stats.buildTickets([leg(0, "2026-10-04T17:00:00Z"), leg(1, "2026-10-06T00:15:00Z")], []);
  assert.equal(parlay.settledDay, "2026-10-05");
});

test("a lost parlay is decided at its earliest losing leg, not its last game", () => {
  const leg = (index, closedAt, legResult) => straight({
    id: "wagerzon:6:leg" + index, venue: "wagerzon", isParlayLeg: true, parlayId: "wagerzon:6", legIndex: index,
    legCount: 3, stake: 100, toWin: 500, status: "lost", closedAt, raw: { legResult },
  });
  // Sunday leg loses; Monday's leg (still to play) has no result yet.
  const [parlay] = stats.buildTickets([
    leg(0, "2026-10-04T17:00:00Z", "WIN"), leg(1, "2026-10-04T20:25:00Z", "LOSE"), leg(2, "2026-10-06T00:15:00Z", ""),
  ], []);
  assert.equal(parlay.settledDay, "2026-10-04");
  assert.equal(parlay.closedAt, "2026-10-04T20:25:00.000Z");
});

test("a Kalshi multivariate combo (the bots' RFQ fills) is its own type, not a straight", () => {
  const [combo] = stats.buildTickets([straight({
    id: "kalshi:KXMVECROSSCATEGORY-X:yes", venue: "kalshi", league: null, betType: "other", side: null,
    raw: { series: "KXMVECROSSCATEGORY", marketTitle: "yes Dodgers, yes Over 8.5" },
  })], []);
  assert.equal(combo.kind, "Kalshi combo");
  assert.equal(combo.market, "Kalshi combo");
  assert.equal(combo.selection, "yes Dodgers, yes Over 8.5");
  assert.deepEqual(stats.groupBy([combo], "kind").map((r) => r.label), ["Kalshi combo"]);
});

test("edge and expected P&L come from the saved fair and the ticket's actual payout", () => {
  // -110 pays 100 on 110 (decimal 1.909); a fair of -130 is p = 0.5652.
  const [ticket] = stats.buildTickets([straight()], [{ betId: "novig:1", fairAmerican: -130 }]);
  const fairProb = 130 / 230;
  assert.ok(Math.abs(ticket.edge - (fairProb * (210 / 110) - 1)) < 1e-12);
  assert.ok(Math.abs(ticket.expected - 110 * ticket.edge) < 1e-9);
  assert.equal(ticket.edgeBucket, "6% and up");
  const [noFair] = stats.buildTickets([straight()], []);
  assert.equal(noFair.edge, null);
  assert.equal(noFair.edgeBucket, stats.NO_FAIR);
});

test("summary: expected and z cover only bets with a fair; ROI interval widens with fewer bets", () => {
  const fairs = [{ betId: "a", fairAmerican: -120 }];
  const tickets = stats.buildTickets([straight({ id: "a" }), straight({ id: "b", status: "lost" })], fairs);
  const total = stats.summarize(tickets);
  assert.equal(total.bets, 2);
  assert.equal(total.withFair, 1);
  assert.equal(total.fairPnl, 100);
  assert.ok(total.z > 0);
  assert.ok(total.ciHalf > 0);
  const many = stats.summarize(stats.buildTickets(
    Array.from({ length: 100 }, (_, i) => straight({ id: "m" + i, status: i % 2 ? "won" : "lost" })), []));
  assert.ok(many.ciHalf < total.ciHalf);
});

test("daily series fills empty days and keeps P&L on the settle day", () => {
  const tickets = stats.buildTickets([
    straight({ id: "a", closedAt: "2026-10-02T20:00:00Z" }),
    straight({ id: "b", status: "lost", closedAt: "2026-10-04T20:00:00Z" }),
  ], []);
  const days = stats.dailySeries(tickets, "2026-10-02", "2026-10-04");
  assert.deepEqual(days.map((d) => [d.day, d.bets, d.pnl]), [["2026-10-02", 1, 100], ["2026-10-03", 0, 0], ["2026-10-04", 1, -110]]);
  assert.equal(stats.firstSettledDay(tickets), "2026-10-02");
  assert.equal(stats.inDayRange(tickets, "2026-10-03", "2026-10-04").length, 1);
});

test("groupBy orders ordinal groups by their scale and others by handle; an unknown key fails loudly", () => {
  const tickets = stats.buildTickets([
    straight({ id: "a", price: 250, stake: 100, toWin: 250 }),
    straight({ id: "b", price: -250, stake: 250, toWin: 100 }),
    straight({ id: "c", venue: "kalshi", stake: 500, toWin: 450 }),
  ], []);
  assert.deepEqual(stats.groupBy(tickets, "oddsBucket").map((r) => r.label), ["−200 or shorter", "−120 to +120", "+200 and up"]);
  assert.deepEqual(stats.groupBy(tickets, "venue").map((r) => r.label), ["Kalshi", "Novig"]);
  assert.throws(() => stats.groupBy(tickets, "nope"), /unknown group nope/);
});

test("calibration bins straight bets with a fair by fair probability, pushes left out", () => {
  const records = [
    straight({ id: "a", status: "won" }), straight({ id: "b", status: "lost" }), straight({ id: "c", status: "push" }),
  ];
  const fairs = records.map((r) => ({ betId: r.id, fairAmerican: -122 }));  // p = 0.5495
  const [bin] = stats.calibration(stats.buildTickets(records, fairs));
  assert.equal(bin.bets, 2);
  assert.equal(bin.low, 0.5);
  assert.equal(bin.actual, 0.5);
});

test("selection labels: spread, total with period, moneyline, and the venue's text when unparsed", () => {
  const label = (overrides) => stats.buildTickets([straight(overrides)], [])[0].selection;
  assert.equal(label({}), "Bills −3.5");
  assert.equal(label({ betType: "total", side: "over", points: 47.5, period: "1H" }), "1H Over 47.5");
  assert.equal(label({ betType: "moneyline", side: "away", points: null }), "Jets ML");
  assert.equal(label({ betType: "other", side: null, raw: { description: "Bills to score first" } }), "Bills to score first");
});

test("rangeBounds: Today and Yesterday are single Pacific days; presets end today", () => {
  const today = "2026-10-06";
  assert.deepEqual(stats.rangeBounds("Today", today, null, "2026-01-02"), { first: today, last: today });
  assert.deepEqual(stats.rangeBounds("Yesterday", "2026-03-01", null, null), { first: "2026-02-28", last: "2026-02-28" });
  assert.deepEqual(stats.rangeBounds("7D", today, null, null), { first: "2026-09-30", last: today });
  assert.deepEqual(stats.rangeBounds("YTD", today, null, null), { first: "2026-01-01", last: today });
  assert.deepEqual(stats.rangeBounds("All", today, null, "2026-01-02"), { first: "2026-01-02", last: today });
  assert.throws(() => stats.rangeBounds("2W", today, null, null), /unknown date range 2W/);
});

test("rangeBounds: Custom takes the picked days, swaps a backwards pair, fills a blank end", () => {
  const today = "2026-10-06";
  assert.deepEqual(stats.rangeBounds("Custom", today, { first: "2026-09-01", last: "2026-09-15" }, "2026-01-02"), { first: "2026-09-01", last: "2026-09-15" });
  assert.deepEqual(stats.rangeBounds("Custom", today, { first: "2026-09-15", last: "2026-09-01" }, "2026-01-02"), { first: "2026-09-01", last: "2026-09-15" });
  assert.deepEqual(stats.rangeBounds("Custom", today, { first: null, last: "2026-02-30" }, "2026-01-02"), { first: "2026-01-02", last: today });
});

test("rangeBounds: Custom clamps to [first settled day, today] so a half-typed year can't span millennia", () => {
  const today = "2026-10-06";
  assert.deepEqual(stats.rangeBounds("Custom", today, { first: "0201-10-01", last: "2026-09-15" }, "2026-01-02"), { first: "2026-01-02", last: "2026-09-15" });
  assert.deepEqual(stats.rangeBounds("Custom", today, { first: "2026-09-01", last: "9999-01-01" }, "2026-01-02"), { first: "2026-09-01", last: today });
  assert.deepEqual(stats.rangeBounds("Custom", today, { first: "2026-09-15", last: "0201-01-01" }, "2026-01-02"), { first: "2026-01-02", last: "2026-09-15" });
  assert.deepEqual(stats.rangeBounds("Custom", today, { first: "2025-01-01", last: "2025-02-01" }, null), { first: today, last: today });
});

test("a removed record marks its ticket, and one removed leg marks the whole parlay", () => {
  const leg = (index) => straight({ id: "bfa:9:leg" + index, venue: "bfa", isParlayLeg: true, parlayId: "bfa:9", legIndex: index, legCount: 2 });
  const tickets = stats.buildTickets([straight({ id: "wz:1", venue: "wagerzon" }), straight({ id: "wz:2", venue: "wagerzon" }), leg(0), leg(1)],
    [], [{ betId: "wz:1" }, { betId: "bfa:9:leg1" }]);
  const byId = Object.fromEntries(tickets.map((t) => [t.id, t]));
  assert.deepEqual([byId["wz:1"].excluded, byId["wz:2"].excluded, byId["bfa:9"].excluded], [true, false, true]);
  assert.deepEqual(byId["bfa:9"].betIds, ["bfa:9:leg0", "bfa:9:leg1"]);
});

test("sortRows: numbers and text both ways, missing keys last, ties stable", () => {
  const rows = [{ id: "a", v: 2 }, { id: "b", v: null }, { id: "c", v: 10 }, { id: "d", v: 2 }, { id: "e", v: NaN }];
  const ids = (list) => list.map((r) => r.id).join("");
  assert.equal(ids(stats.sortRows(rows, (r) => r.v, "asc")), "adcbe");
  assert.equal(ids(stats.sortRows(rows, (r) => r.v, "desc")), "cadbe");
  const names = [{ n: "novig" }, { n: "BFA" }, { n: "Kalshi" }];
  assert.deepEqual(stats.sortRows(names, (r) => r.n, "asc").map((r) => r.n), ["BFA", "Kalshi", "novig"]);
  assert.equal(rows[0].id, "a", "the input is not reordered");
});
