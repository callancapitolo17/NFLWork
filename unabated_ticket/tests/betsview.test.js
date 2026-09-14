// Run: node --test unabated_ticket/tests
// betsview.js: the panel's bet-history presentation helpers (#114).
const test = require("node:test");
const assert = require("node:assert/strict");
const view = require("../extension/betsview.js");
const fs = require("node:fs");
const path = require("node:path");
const teams = require("../extension/teams.js");
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

test("badgeText / badgeKind: dollars held win, then dollars against, then a plain game marker", () => {
  const flag = (tier, held, against) => ({ tier, matches: [], exposure: { held, against, heldBets: [], againstBets: [] } });
  assert.equal(view.badgeText(flag("same_line", 300, 0)), "held $300");
  assert.equal(view.badgeText(flag("same_side", 12.5, 20)), "held $12.50");
  assert.equal(view.badgeText(flag("opposite", 0, 200)), "against $200");
  assert.equal(view.badgeText(flag("same_game", 0, 0)), "game");
  assert.equal(view.badgeText({ tier: null, matches: [], exposure: { held: 0, against: 0 } }), null);
  assert.equal(view.badgeText(null), null);
  // A venue that gave no stake: the tier still names the kind, without dollars.
  assert.equal(view.badgeText(flag("opposite", 0, 0)), "against");
  assert.equal(view.badgeText(flag("same_side", 0, 0)), "held");
  assert.equal(view.badgeKind(flag("opposite", 0, 0)), "against");
  assert.deepEqual(["same_line", "opposite", "same_game"].map((tier) => view.badgeKind(flag(tier, tier === "same_line" ? 1 : 0, tier === "opposite" ? 1 : 0))), ["held", "against", "game"]);
});

test("stakeAdviceLine: the verb says what the number is, then what is held and full size", () => {
  const exposure = (held, against) => ({ held, against, heldBets: [], againstBets: [] });
  assert.equal(view.stakeAdviceLine(view.stakeAdvice(500, exposure(0, 0))), null);
  assert.equal(view.stakeAdviceLine(view.stakeAdvice(600, exposure(350, 0))), "add $250, $350 held · full size $600");
  assert.equal(view.stakeAdviceLine(view.stakeAdvice(520, exposure(600, 0))), "bet $0, $600 held · full size $520");
  assert.equal(view.stakeAdviceLine(view.stakeAdvice(null, exposure(600, 0))), "bet $0, $600 held · full size none here");
  assert.equal(view.stakeAdviceLine(view.stakeAdvice(500, exposure(0, 200))), "bet $500, $200 on the other side (net $300 on this side)");
  assert.equal(view.stakeAdviceLine(view.stakeAdvice(100, exposure(0, 200))), "bet $100, $200 on the other side (still $100 against)");
  assert.equal(view.stakeAdviceLine(view.stakeAdvice(null, exposure(0, 200))), "bet $0, $200 on the other side");
  assert.deepEqual(view.stakeAdviceWords(view.stakeAdvice(600, exposure(350.5, 0))),
    { verb: "add", have: "$350.50", target: "$600", bet: "$249.50", note: null, against: false });
});

test("stakeAdvice: none, add the difference, at size when held covers the stake, reverse with the net", () => {
  const exposure = (held, against) => ({ held, against, heldBets: [], againstBets: [] });
  assert.deepEqual(view.stakeAdvice(500, exposure(0, 0)), { kind: "none" });
  assert.deepEqual(view.stakeAdvice(500, exposure(300, 0)), { kind: "add", add: 200, held: 300, stake: 500 });
  assert.deepEqual(view.stakeAdvice(520, exposure(600, 0)), { kind: "at_size", held: 600, stake: 520 });
  assert.deepEqual(view.stakeAdvice(null, exposure(600, 0)), { kind: "at_size", held: 600, stake: null });
  assert.deepEqual(view.stakeAdvice(0, exposure(100, 0)), { kind: "at_size", held: 100, stake: null });
  assert.deepEqual(view.stakeAdvice(500, exposure(0, 200)), { kind: "reverse", against: 200, stake: 500, net: 300 });
  assert.deepEqual(view.stakeAdvice(100, exposure(0, 200)), { kind: "reverse", against: 200, stake: 100, net: -100 });
  assert.deepEqual(view.stakeAdvice(500, exposure(300, 200)).kind, "add");
});

test("relatedLines: one line per match, each naming the bet and its tier", () => {
  const flag = {
    tier: "same_line",
    matches: [
      { tier: "same_line", label: "You bet this: Chattanooga -5.5 +138 · 42.0¢ · $300 @ Kalshi · Sep 10 2:15 PM" },
      { tier: "opposite", label: "You are on the OTHER side: Eastern Kentucky +5.5 -150 · 60.0¢ · $200 @ Kalshi" },
    ],
  };
  assert.deepEqual(view.relatedLines(flag), [
    { tier: "same_line", text: "You bet this: Chattanooga -5.5 +138 · 42.0¢ · $300 @ Kalshi · Sep 10 2:15 PM" },
    { tier: "opposite", text: "You are on the OTHER side: Eastern Kentucky +5.5 -150 · 60.0¢ · $200 @ Kalshi" },
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

// ---- page-sourced venues (#116: Novig content script) ----

function novigRead(overrides) {
  return Object.assign({ bets: [record("novig:1", "open", { venue: "novig", source: "novig_page" })], readAt: iso(20e3), url: "https://app.novig.us/portfolio", error: null, complete: true, pageSeenAt: iso(0) }, overrides);
}

test("sourceRows: a content-script venue is configured with its read age; a stale or missing read carries the refresh hint", () => {
  const fresh = view.sourceRows(null, NOW, { novig: novigRead() }).find((row) => row.venue === "novig");
  assert.equal(fresh.configured, true);
  assert.equal(fresh.level, "green");
  assert.equal(fresh.ageText, "20 s");
  assert.equal(fresh.count, 1);
  assert.equal(fresh.note, null);
  const stale = view.sourceRows(null, NOW, { novig: novigRead({ readAt: iso(2 * 3600e3), pageSeenAt: iso(2 * 3600e3) }) }).find((row) => row.venue === "novig");
  assert.equal(stale.level, "red");
  assert.equal(stale.note, "open app.novig.us and its Portfolio screen in a tab to refresh");
  const tabOpen = view.sourceRows(null, NOW, { novig: novigRead({ readAt: iso(2 * 3600e3), pageSeenAt: iso(30e3) }) }).find((row) => row.venue === "novig");
  assert.equal(tabOpen.note, "Novig tab is open — open its Portfolio screen to refresh");
  const never = view.sourceRows(null, NOW, { novig: novigRead({ readAt: null, bets: [] }) }).find((row) => row.venue === "novig");
  assert.equal(never.ageText, "never");
  assert.equal(never.count, 0);
  const errored = view.sourceRows(null, NOW, { novig: novigRead({ error: "ActivePortfolioOrders_Query: boom" }) }).find((row) => row.venue === "novig");
  assert.equal(errored.error, "ActivePortfolioOrders_Query: boom");
  assert.equal(view.sourceRows(null, NOW, {}).find((row) => row.venue === "novig").configured, false);
});

test("headerLine and sourcesUnavailable: a fresh Novig read counts as a live source", () => {
  assert.equal(view.headerLine([], null, NOW, { novig: novigRead() }), "bets: 0 open · kalshi — · betonline — · novig 20 s · prophetx —");
  assert.equal(view.sourcesUnavailable(null, NOW, { novig: novigRead() }), false);
  assert.equal(view.sourcesUnavailable(null, NOW, { novig: novigRead({ readAt: iso(2 * 3600e3) }) }), true);
});

test("mergePageSource: a complete read is authoritative for its venue; an incomplete one only adds; other venues untouched", () => {
  const stored = [record("novig:old", "open", { venue: "novig" }), record("kalshi:k", "open")];
  const complete = view.mergePageSource(stored, "novig", novigRead(), NOW);
  assert.deepEqual(complete.map((r) => r.id).sort(), ["kalshi:k", "novig:1"]);
  const partial = view.mergePageSource(stored, "novig", novigRead({ complete: false }), NOW);
  assert.deepEqual(partial.map((r) => r.id).sort(), ["kalshi:k", "novig:1", "novig:old"]);
  assert.ok(complete.every((r) => r.awayKey === key("cfb", "Chattanooga")));
  assert.deepEqual(view.mergePageSource(stored, "novig", null, NOW).map((r) => r.id).sort(), ["kalshi:k", "novig:old"]);
});
