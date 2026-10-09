// Run: node --test unabated_ticket/tests
// edgerows.js: the Edges list's rows without the DOM — what the side panel
// and the server runner (server/runner.js) both list for the same settings.
// Driven by the real NFL slice (fixtures/v2_slice.json: Chicago Bears @
// Carolina Panthers, kickoff 2026-09-13T17:00Z), an hour before kickoff.
const test = require("node:test");
const assert = require("node:assert/strict");
const fs = require("node:fs");
const path = require("node:path");
const feed = require("../extension/feed.js");
const kelly = require("../extension/kelly.js");
const teams = require("../extension/teams.js");
const betsLib = require("../extension/bets.js");
const edgemove = require("../extension/edgemove.js");
const tailflex = require("../extension/tailflex.js");
const teaser = require("../extension/teaser.js");
const edgeRows = require("../extension/edgerows.js");

const KICKOFF_MS = Date.parse("2026-09-13T17:00:00Z");
const NOW = KICKOFF_MS - 3600 * 1000;
const NFL = 1;
const NOVIG = 89;
const KALSHI = 105;
const SOUTHPOINT = 99;
const BEARS_SPREAD_SOUTHPOINT = "289357360:ms99:si0:tid6";
const PANTHERS_ALT_NOVIG = "289357357:ms89:si1:tid5:alt-2.5";
const PANTHERS_ALT_65_NOVIG = "289357357:ms89:si1:tid5:alt-6.5";
const BEARS_ML_BETMGM = "289357353:ms4:si0:tid6";

function sliceState() {
  const snapshot = JSON.parse(fs.readFileSync(path.join(__dirname, "fixtures", "v2_slice.json"), "utf8"));
  return feed.mergeStates([feed.parseSnapshot(snapshot, { leagueId: NFL })]);
}

// The slice's lines are days old by NOW (SouthPoint's main spread changed
// 2026-08-29), so the age gate is opened and every live book allowed; the
// list's floor is 1% so the slice's thinner lines are in it.
function openSettings(overrides) {
  return { ...edgeRows.DEFAULT_EDGE_SETTINGS, bookIds: null, maxLineAgeHours: 1e6, includeAlts: true, minEdgePct: 1, ...overrides };
}

function contextFor(state, edgeSettings, records, teasers) {
  return {
    edgeSettings, stakeSettings: edgeRows.DEFAULT_STAKE_SETTINGS, records: records || [],
    boardLines: edgeRows.boardLines(state), ladderReaderOf: edgeRows.createLadderReaders(state), teasers: teasers || [],
    measurement: tailflex.measureTailFlex(state, { now: NOW, maxLineAgeMs: edgeSettings.maxLineAgeHours * 3600 * 1000 }),
    effective: edgeRows.effectiveFilter(edgeSettings, state, null, null), now: NOW,
  };
}

function registerSliceTeams(state) {
  for (const [league, list] of Object.entries(edgeRows.feedTeamsByLeague(state))) teams.registerTeams(league, list);
}

// An open BetOnline bet on the Bears -2.5, team keys resolved off the slice's own team list.
function bearsSpreadHeld(state) {
  registerSliceTeams(state);
  return betsLib.resolveTeamKeys([{
    id: "bol-1", source: "betonline", venue: "betonline", league: "nfl", status: "open", sourceFetchedAt: "2026-09-13T15:00:00Z",
    awayTeam: "Chicago Bears", homeTeam: "Carolina Panthers", awayKey: null, homeKey: null,
    betType: "spread", period: "FG", side: "away", points: -2.5, price: -110, stake: 200, toWin: 181.82,
    placedAt: "2026-09-13T14:00:00Z", closedAt: null, eventStart: "2026-09-13T17:00:00Z", isParlayLeg: false, unmatchable: null, approx: [],
  }]);
}

// An open $400 BFA 2-team teaser: Panthers +3.5 (rotation 466) on the slice's
// game and a leg on a game already started elsewhere, which counts as won.
function panthersTeaserOpen(state) {
  const leg = (legIndex, fields) => ({
    id: `bfa:T1:leg${legIndex}`, source: "bfa_api", venue: "bfa", status: "open", league: "nfl", betType: "spread",
    period: "FG", side: "home", points: 3.5, price: -110, rotation: 466, awayTeam: null, homeTeam: "Carolina Panthers", awayKey: null, homeKey: null,
    eventStart: new Date(KICKOFF_MS).toISOString(), stake: 400, toWin: 1200, isParlayLeg: true, parlayId: "bfa:T1", legIndex,
    legCount: 2, placedAt: "2026-09-13T12:00:00Z", approx: ["side_from_rotation_parity"], unmatchable: null,
    raw: { headerDescription: "2 TEAM TEASERS" }, ...fields,
  });
  const records = [leg(0), leg(1, { rotation: 901, homeTeam: "Elsewhere", eventStart: new Date(NOW - 2 * 3600 * 1000).toISOString() })];
  const tickets = teaser.openTeasers(records, edgeRows.boardLines(state), { now: NOW, ladderOf: teaser.teaserBoardOf(state).ladderOf });
  return { records, tickets };
}

test("defaults are the panel's: $30,000 at quarter Kelly, 2.5% edge, 168 h, $100 to win, no alts, cards, every league, the default books", () => {
  assert.deepEqual(edgeRows.DEFAULT_STAKE_SETTINGS, { bankroll: 30000, multiplier: 0.25 });
  const defaults = edgeRows.DEFAULT_EDGE_SETTINGS;
  assert.deepEqual(defaults.leagues, Object.keys(feed.LEAGUES).map(Number));
  assert.deepEqual({ ...defaults, leagues: null }, {
    leagues: null, periods: [1], betTypes: [1, 2, 3], bookIds: undefined, minEdgePct: 2.5, maxLineAgeHours: 168, sortBy: "edge",
    minStake: 0, minLiquidityToWin: 100, includeAlts: false, groupByMarket: true,
  });
  assert.equal(edgeRows.MAX_EDGE_ROWS, 200);
  assert.ok(edgeRows.DEFAULT_BOOK_NAMES.includes("Novig") && edgeRows.DEFAULT_BOOK_NAMES.includes("Kalshi"));
});

test("sanitizeStakeSettings: numbers kept, numeric strings read, missing / zero / junk is the default", () => {
  assert.deepEqual(edgeRows.sanitizeStakeSettings({ bankroll: 8000, multiplier: 0.5 }), { bankroll: 8000, multiplier: 0.5 });
  assert.deepEqual(edgeRows.sanitizeStakeSettings({ bankroll: "8000" }), { bankroll: 8000, multiplier: 0.25 });
  assert.deepEqual(edgeRows.sanitizeStakeSettings({ bankroll: 0, multiplier: "x" }), { bankroll: 30000, multiplier: 0.25 });
  assert.deepEqual(edgeRows.sanitizeStakeSettings(null), { bankroll: 30000, multiplier: 0.25 });
});

test("sanitizeEdgeSettings: unknown leagues and periods dropped, a stored null book choice kept, junk is the default, the retired alt distance ignored", () => {
  const settings = edgeRows.sanitizeEdgeSettings({ leagues: [1, 999], periods: [1, 2, 42], bookIds: null, minEdgePct: -1, sortBy: "nope", includeAlts: true, altMaxDistance: 3 });
  assert.deepEqual(settings.leagues, [1]);
  assert.deepEqual(settings.periods, [1, 2]);
  assert.equal(settings.bookIds, null);
  assert.equal(settings.minEdgePct, 2.5);
  assert.equal(settings.sortBy, "edge");
  assert.equal(settings.includeAlts, true);
  assert.equal("altMaxDistance" in settings, false);
  assert.equal(edgeRows.sanitizeEdgeSettings({}).bookIds, undefined);
  assert.deepEqual(edgeRows.sanitizeEdgeSettings({ bookIds: [89, "x", 105] }).bookIds, [89, 105]);
});

test("effectiveFilter: default books by name, your ticks, the Unabated selection, else every live book", () => {
  const state = sliceState();
  const base = edgeRows.DEFAULT_EDGE_SETTINGS;
  const byDefault = edgeRows.effectiveFilter(base, state, null, null);
  assert.equal(byDefault.mode, "default");
  assert.deepEqual(Array.from(byDefault.bookIds).sort((a, b) => a - b), [NOVIG, KALSHI]);
  assert.deepEqual(Array.from(byDefault.betTypeIds), [1, 2, 3]);
  const custom = edgeRows.effectiveFilter({ ...base, bookIds: [SOUTHPOINT] }, state, [NOVIG], null);
  assert.deepEqual([custom.mode, Array.from(custom.bookIds)], ["custom", [SOUTHPOINT]]);
  const unabated = edgeRows.effectiveFilter({ ...base, bookIds: null }, state, [NOVIG], { at: NOW });
  assert.deepEqual([unabated.mode, Array.from(unabated.bookIds)], ["unabated", [NOVIG]]);
  const all = edgeRows.effectiveFilter({ ...base, bookIds: null }, state, null, null);
  assert.deepEqual([all.mode, all.bookIds, all.filter], ["all", null, null]);
  assert.deepEqual(edgeRows.liveBooks(state).map((book) => book.name), ["BetMGM", "Kalshi", "Novig", "SouthPoint", "Sports Interaction"]);
  assert.deepEqual(Array.from(edgeRows.edgeSelectionOptions(base, byDefault, NOW).leagueIds), base.leagues);
});

test("listedEdgeRows: nothing held, every row's stake is quarter Kelly off its edge, its rank score keep x edge x stake, edge order", () => {
  const state = sliceState();
  const rows = edgeRows.listedEdgeRows(state, contextFor(state, openSettings()));
  assert.equal(rows.length, 14);
  assert.deepEqual(rows.slice(0, 3).map((row) => row.key), [PANTHERS_ALT_NOVIG, PANTHERS_ALT_65_NOVIG, BEARS_SPREAD_SOUTHPOINT]);
  const measurement = tailflex.measureTailFlex(state, { now: NOW, maxLineAgeMs: 1e6 * 3600 * 1000 });
  for (const row of rows) {
    const standalone = kelly.kellyStakeFromEdge({ bookPrice: row.price, edgePct: row.edgePct, bankroll: 30000, multiplier: 0.25 }).stake;
    assert.equal(row.stake, standalone);
    assert.equal(row.rankScore, tailflex.rankOfRow(row, measurement).score);
    assert.ok(row.rankScore > 0 && row.rankScore < row.edgePct / 100 * row.stake);
    assert.deepEqual([row.bet.tier, row.bet.advice.kind, row.bet.advice.bet], [null, "none", standalone]);
  }
  for (let i = 1; i < rows.length; i += 1) assert.ok(rows[i - 1].edgePct >= rows[i].edgePct);
});

test("listedEdgeRows: the settings' gates — min edge, alts off, the default books, max line age, leagues", () => {
  const state = sliceState();
  const list = (overrides) => edgeRows.listedEdgeRows(state, contextFor(state, openSettings(overrides))).map((row) => row.key);
  assert.deepEqual(list({ minEdgePct: 5 }), [PANTHERS_ALT_NOVIG, PANTHERS_ALT_65_NOVIG, BEARS_SPREAD_SOUTHPOINT]);
  assert.equal(list({ minEdgePct: 2.5 }).length, 8);
  assert.deepEqual(list({ includeAlts: false }), [BEARS_SPREAD_SOUTHPOINT, BEARS_ML_BETMGM, "289357353:ms69:si0:tid6", "289357353:ms99:si0:tid6"]);
  // The default books in the slice are Novig and Kalshi, which post only alts here.
  assert.deepEqual(list({ bookIds: undefined, includeAlts: false }), []);
  // The panel's own default, a week: Novig's and Kalshi's alts (changed 2026-09-11) stay, the late-August main lines go.
  const withinAWeek = list({ maxLineAgeHours: 168 });
  assert.equal(withinAWeek.length, 10);
  assert.ok(withinAWeek.every((key) => key.includes(":alt")));
  // Football unticked: the scanner still holds NFL for the Teasers tab, the Edges list does not show it.
  assert.deepEqual(list({ leagues: [] }), []);
});

test("listedEdgeRows: sort by stake and by start", () => {
  const state = sliceState();
  const byStake = edgeRows.listedEdgeRows(state, contextFor(state, openSettings({ sortBy: "stake" })));
  assert.equal(byStake[0].key, BEARS_SPREAD_SOUTHPOINT);
  for (let i = 1; i < byStake.length; i += 1) assert.ok(byStake[i - 1].stake >= byStake[i].stake);
  assert.deepEqual(edgeRows.sortEdgeRows([{ eventStartMs: 2, edgePct: 1 }, { eventStartMs: 1, edgePct: 1 }], "start").map((row) => row.eventStartMs), [1, 2]);
});

test("conditional Kelly: a held Bears -2.5 makes the same line an add, the other side bigger, and a minimum suggested bet hides the $0 adds", () => {
  const state = sliceState();
  const records = bearsSpreadHeld(state);
  const rows = edgeRows.listedEdgeRows(state, contextFor(state, openSettings(), records));
  const byKey = new Map(rows.map((row) => [row.key, row]));
  const same = byKey.get(BEARS_SPREAD_SOUTHPOINT);
  assert.deepEqual([same.bet.tier, same.bet.advice.kind, same.bet.advice.verb, same.bet.advice.held, same.bet.advice.bet], ["same_line", "sized", "add", 200, 237.25]);
  assert.deepEqual(edgeRows.stakeRail(same), { text: "add $237.25", note: "$437.25 alone", atSize: false });
  const other = byKey.get(PANTHERS_ALT_NOVIG);
  assert.deepEqual([other.bet.tier, other.bet.advice.against, other.bet.advice.bet], ["opposite", 200, 328.86]);
  assert.deepEqual(edgeRows.stakeRail(other), { text: "bet $328.86", note: "$215.13 alone", atSize: false });
  const atSize = rows.filter((row) => row.bet.advice.bet === 0);
  assert.equal(atSize.length, 5);
  assert.ok(atSize.every((row) => edgeRows.stakeRail(row).atSize));
  const gated = edgeRows.listedEdgeRows(state, contextFor(state, openSettings({ minStake: 1 }), records));
  assert.equal(gated.length, rows.length - 5);
  assert.equal(edgeRows.exposureDollars(same), 200);
  const byExposure = edgeRows.listedEdgeRows(state, contextFor(state, openSettings({ sortBy: "exposure" }), records));
  assert.equal(edgeRows.exposureDollars(byExposure[0]), 200);
  assert.equal(edgeRows.exposureDollars(byExposure[byExposure.length - 1]), 0);
});

test("open BFA teasers size the rows too: a Panthers +3.5 teaser puts the Panthers lines at size and grows the Bears side", () => {
  const state = sliceState();
  registerSliceTeams(state);
  const { records, tickets } = panthersTeaserOpen(state);
  assert.equal(tickets.length, 1);
  const rows = edgeRows.listedEdgeRows(state, contextFor(state, openSettings(), records, tickets));
  const byKey = new Map(rows.map((row) => [row.key, row]));
  const panthers = byKey.get(PANTHERS_ALT_NOVIG);
  assert.deepEqual([panthers.bet.advice.verb, panthers.bet.advice.bet, panthers.bet.advice.teasers], ["add", 0, { held: 400, against: 0 }]);
  assert.deepEqual(edgeRows.stakeRail(panthers), { text: "add $0.00", note: "$215.13 alone", atSize: true });
  assert.equal(edgeRows.exposureDollars(panthers), 400);
  const bears = byKey.get(BEARS_SPREAD_SOUTHPOINT);
  assert.deepEqual([bears.bet.advice.verb, bears.bet.advice.bet, bears.bet.advice.teasers], ["bet", 1187.54, { held: 0, against: 400 }]);
  // Without the tickets the same records size nothing: the teaser legs are the tickets' own.
  const unsized = edgeRows.listedEdgeRows(state, contextFor(state, openSettings(), records, []));
  const bearsUnsized = unsized.find((row) => row.key === BEARS_SPREAD_SOUTHPOINT);
  assert.equal(bearsUnsized.bet.advice.bet, bearsUnsized.stake);
});

test("groupEdgeRows: one card per (game, period, bet type, side), its best line the highest tail-flex rank score", () => {
  const state = sliceState();
  const rows = edgeRows.listedEdgeRows(state, contextFor(state, openSettings()));
  const cards = edgeRows.groupEdgeRows(rows, "edge");
  assert.deepEqual(cards.map((card) => [card.key, card.rows.length, card.best.key]), [
    ["125807:pt1:bt2:si1", 4, PANTHERS_ALT_NOVIG],
    ["125807:pt1:bt2:si0", 5, BEARS_SPREAD_SOUTHPOINT],
    ["125807:pt1:bt1:si0", 3, BEARS_ML_BETMGM],
    ["125807:pt1:bt3:si0", 2, "289357345:ms89:si0:tid6:alt54.5"],
  ]);
  for (const card of cards) {
    for (let i = 1; i < card.rows.length; i += 1) assert.ok(card.rows[i - 1].rankScore >= card.rows[i].rankScore);
  }
  // The Bears card: SouthPoint's main -2.5 at $437 ranks far above Novig's
  // deep alts, though the -13.5 alt is within a point of its edge.
  const bearsCard = cards[1];
  assert.ok(bearsCard.rows[0].rankScore > 100 * bearsCard.rows[bearsCard.rows.length - 1].rankScore);
  assert.deepEqual(edgeRows.groupEdgeRows(rows, "stake").map((card) => card.best.key), [BEARS_SPREAD_SOUTHPOINT, BEARS_ML_BETMGM, PANTHERS_ALT_NOVIG, "289357345:ms89:si0:tid6:alt54.5"]);
  assert.equal(edgeRows.groupEdgeRows(rows, "start").length, 4);
});

test("describeTailFlex: the c per spread / total market on the list, the fallback when too thin to measure; empty list, empty text", () => {
  const state = sliceState();
  const context = contextFor(state, openSettings());
  const rows = edgeRows.listedEdgeRows(state, context);
  assert.equal(edgeRows.describeTailFlex(rows, context.measurement), "tail flex: NFL spr 10% · NFL tot 10%");
  assert.equal(edgeRows.describeTailFlex(rows.filter((row) => row.betTypeId === 1), context.measurement), "");
  assert.equal(edgeRows.describeTailFlex([], context.measurement), "");
});

test("stakeRail: no flag is the standalone bet, no stake is a dash, a liquidity cap says so", () => {
  assert.deepEqual(edgeRows.stakeRail({ stake: 437.25 }), { text: "bet $437.25", note: null, atSize: false });
  assert.deepEqual(edgeRows.stakeRail({ stake: null }), { text: "—", note: null, atSize: false });
  const capped = { kind: "none", bet: 17, alone: 32.58, verb: "bet", held: 0, against: 0, teasers: { held: 0, against: 0 }, reason: null, matches: [], cappedAt: 17 };
  assert.deepEqual(edgeRows.stakeRail({ stake: 17, bet: { tier: null, matches: [], advice: capped } }), { text: "bet $17.00", note: "all $17 liq", atSize: false });
  assert.equal(edgeRows.edgeTier(4), "hot");
  assert.equal(edgeRows.edgeTier(2.5), "warm");
  assert.equal(edgeRows.edgeTier(1.2), "thin");
});

test("moveTag: only on a line held in this direction, and it names the fair's move", () => {
  const state = sliceState();
  const records = bearsSpreadHeld(state);
  const line = state.lines[BEARS_SPREAD_SOUTHPOINT];
  const history = {};
  edgemove.observe(history, line, { at: NOW - 5 * 60 * 1000 });
  // Unabated's fair on the Bears moves from its slice value to -150: well over half a point of probability.
  edgemove.observe(history, { ...line, bacr: -150 }, { at: NOW - 60 * 1000 });
  const rows = edgeRows.listedEdgeRows(state, contextFor(state, openSettings(), records));
  const held = rows.find((row) => row.key === BEARS_SPREAD_SOUTHPOINT);
  const context = { history, fillFairIndex: new Map(), now: NOW };
  const tag = edgeRows.moveTag(held, held.bet, context);
  assert.equal(tag.kind, "fair_to_you");
  assert.equal(tag.label, "fair moved to you");
  assert.match(tag.detail, /^fair \d+\.\d% → 60\.0%$/);
  assert.equal(edgeRows.moveWords(held, context), "fair moved to you");
  const notHeld = rows.find((row) => row.key === PANTHERS_ALT_NOVIG);
  assert.equal(edgeRows.moveTag(notHeld, notHeld.bet, context), null);
});

test("applyBetsPayload: a body with no bets throws; an older service's missing crosswalk, pins and fill fairs keep the held ones", () => {
  const held = { records: [], crosswalk: [{ venue: "novig" }], pins: [{ betId: "x" }], fillFairs: [] };
  assert.throws(() => edgeRows.applyBetsPayload(held, { sources: {} }, NOW), /bets\.json has no bets array/);
  const older = edgeRows.applyBetsPayload(held, { bets: [], generatedAt: "2026-09-13T15:59:00Z" }, NOW);
  assert.equal(older.crosswalk, held.crosswalk);
  assert.equal(older.pins, held.pins);
  assert.equal(older.fillFairs, held.fillFairs);
  assert.deepEqual([older.generatedAt, older.sources], ["2026-09-13T15:59:00Z", {}]);
  const current = edgeRows.applyBetsPayload(held, { bets: [], crosswalk: [], pins: [], fillFairs: [], sources: { kalshi: { ok: true } } }, NOW);
  assert.deepEqual([current.crosswalk, current.pins, current.sources], [[], [], { kalshi: { ok: true } }]);
});

test("serviceSettingsOf: every field for PUT /settings.json; settingsFromService gives the same settings back", () => {
  const stake = { bankroll: 8000, multiplier: 0.5 };
  for (const bookIds of [undefined, null, [89, 59]]) {
    const edges = { ...edgeRows.DEFAULT_EDGE_SETTINGS, bookIds, minEdgePct: 3, includeAlts: true, leagues: [1, 2] };
    const body = edgeRows.serviceSettingsOf(stake, edges);
    assert.deepEqual(Object.keys(body).sort(), ["bankroll", "multiplier", ...edgeRows.EDGE_SETTING_KEYS, "bookMode", "bookIds"].sort());
    assert.equal(body.bookMode, bookIds === undefined ? "default" : bookIds === null ? "all" : "custom");
    assert.deepEqual(edgeRows.settingsFromService(body), { stakeSettings: stake, edgeSettings: edges });
  }
});

test("applyBetsPayload: the shared Dismiss and Can't tease marks, null from a service that has none", () => {
  const held = { records: [], crosswalk: [], pins: [], fillFairs: [] };
  const old = edgeRows.applyBetsPayload(held, { bets: [] }, 0);
  assert.deepEqual([old.dismissals, old.teaserBlocks], [null, null]);
  const marks = { dismissals: [{ betId: "a", dismissedAt: "x" }], teaserBlocks: [{ marketKey: "1:bt2", eventStartMs: 5, blockedAt: "y" }] };
  const now = edgeRows.applyBetsPayload(held, { bets: [], ...marks }, 0);
  assert.deepEqual([now.dismissals, now.teaserBlocks], [marks.dismissals, marks.teaserBlocks]);
});
