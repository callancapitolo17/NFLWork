// Run: node --test unabated_ticket/tests
// betsview.js: the panel's bet-history presentation helpers (#114).
const test = require("node:test");
const assert = require("node:assert/strict");
const view = require("../extension/betsview.js");

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
  assert.deepEqual(rows.map((row) => row.venue), ["kalshi", "betonline", "novig", "prophetx"]);
  const [kalshi, betonline, novig, prophetx] = rows;
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
});

test("sourceRows: a configured source that never succeeded is red with 'never'", () => {
  const [kalshi] = view.sourceRows(payloadWith({ kalshi: { fetchedAt: null, ok: false, error: "boom", count: 0 } }), NOW);
  assert.equal(kalshi.level, "red");
  assert.equal(kalshi.ageText, "never");
  assert.equal(kalshi.error, "boom");
});

test("sourceRows: no payload at all is four unconfigured rows", () => {
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
  assert.equal(view.headerLine(records, payload, NOW), "bets: 2 open · kalshi 20 s · betonline — · novig — · prophetx —");
  assert.equal(view.headerLine([], null, NOW), "bets: 0 open · kalshi — · betonline — · novig — · prophetx —");
});

test("bannerLines: at most five, strongest first as given, and the count of the rest", () => {
  const matches = Array.from({ length: 7 }, (_, i) => ({ tier: "same_game", label: `m${i}` }));
  const { shown, more } = view.bannerLines(matches);
  assert.deepEqual(shown.map((m) => m.label), ["m0", "m1", "m2", "m3", "m4"]);
  assert.equal(more, 2);
  assert.deepEqual(view.bannerLines(matches.slice(0, 3)), { shown: matches.slice(0, 3), more: 0 });
  assert.deepEqual(view.bannerLines([]), { shown: [], more: 0 });
});

test("badgeText: BET for same_line and same_side, OTHER SIDE, GAME, null otherwise", () => {
  assert.equal(view.badgeText("same_line"), "BET");
  assert.equal(view.badgeText("same_side"), "BET");
  assert.equal(view.badgeText("opposite"), "OTHER SIDE");
  assert.equal(view.badgeText("same_game"), "GAME");
  assert.equal(view.badgeText(null), null);
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
  assert.equal(x.awayKey, "cfb:chattanooga");
  assert.equal(x.homeKey, "cfb:eastern-kentucky");
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
  assert.deepEqual(view.sanitizeBetsSettings(null), { serviceUrl: "http://127.0.0.1:8094", hideBet: false });
  assert.deepEqual(view.sanitizeBetsSettings({ serviceUrl: "http://localhost:9000/", hideBet: true }), { serviceUrl: "http://localhost:9000", hideBet: true });
  assert.deepEqual(view.sanitizeBetsSettings({ serviceUrl: "not a url", hideBet: "yes" }), { serviceUrl: "http://127.0.0.1:8094", hideBet: false });
});
