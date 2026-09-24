// Run: node --test unabated_ticket/tests
// fillfair.js: the fair a bet's own line had when it was placed (capture) and
// the since-your-first-fill reading of the edge-move tag (display). The board
// is the 49ers card that motivated it (2026-09-23): Seahawks @ 49ers, a 49ers
// -14.5 alt of the -8.5 main, held on Novig at +208. Probabilities: +178 =
// 36.0%, +188 = 34.7%, +170 = 37.0%, +208 = 32.5%, +223 = 31.0%.
const test = require("node:test");
const assert = require("node:assert/strict");
const feed = require("../extension/feed.js");
const edgemove = require("../extension/edgemove.js");
const fillfair = require("../extension/fillfair.js");
const { captureFillFairs, fillTimeEntry, fairsByBetId, fillPriceOf, baselineOf, moveSinceFill, REFUSED } = fillfair;

const MIN = 60 * 1000;
const T0 = Date.UTC(2026, 8, 23, 17, 0, 0);
const KICKOFF = T0 + 3 * 60 * MIN;
const WATCHING_SINCE = T0 - 60 * MIN;
const iso = (ms) => new Date(ms).toISOString();

const EVENT_ID = 900;
const NOVIG = 89;
const KALSHI = 105;
const DRAFTKINGS = 4;
const NINERS_ALT = (book) => `m1:ms${book}:si1:tid22:alt-14.5`;
const SEAHAWKS_ALT = (book) => `m1:ms${book}:si0:tid11:alt14.5`;

// One feed line of the Seahawks @ 49ers spread market (feed.normalizeAltLine's shape).
function line(bookId, overrides = {}) {
  const sideIndex = overrides.sideIndex ?? 1;
  const points = overrides.points ?? (sideIndex === 1 ? -14.5 : 14.5);
  const sideKey = sideIndex === 1 ? "si1:tid22" : "si0:tid11";
  const isAlt = overrides.isAlt ?? true;
  const mainKey = `m1:ms${bookId}:${sideKey}`;
  return {
    key: isAlt ? `${mainKey}:alt${points}` : mainKey, isAlt, mainKey, mainPoints: sideIndex === 1 ? -8.5 : 8.5,
    leagueId: 1, periodTypeId: 1, betTypeId: 2, eventId: EVENT_ID, marketId: "m1", bookId, sideKey, sideIndex, points,
    price: 223, sourceFormat: 1, sourcePrice: null, bacr: 188, ge: 0.12, liquidity: 500, statusId: 1,
    sequenceNumber: T0, isBlurred: false, modifiedOn: null,
    ...overrides,
  };
}

function boardState(lines) {
  const state = feed.emptyState();
  state.leagues = [1];
  state.teams = { 11: "Seattle Seahawks", 22: "San Francisco 49ers" };
  state.books = {
    [NOVIG]: { id: NOVIG, name: "Novig", isLive: true, hasLiquidity: true },
    [KALSHI]: { id: KALSHI, name: "Kalshi", isLive: true, hasLiquidity: true },
    [DRAFTKINGS]: { id: DRAFTKINGS, name: "DraftKings", isLive: true, hasLiquidity: false },
  };
  state.events[EVENT_ID] = {
    eventId: EVENT_ID, leagueId: 1, eventName: "Seahawks @ 49ers", eventStart: KICKOFF,
    awayTeamId: 11, homeTeamId: 22, awayRotation: 101, homeRotation: 102, venueIds: null,
  };
  const main = line(DRAFTKINGS, { isAlt: false, points: -8.5, price: -110, bacr: -112 });
  for (const each of [main, ...lines]) state.lines[each.key] = each;
  return state;
}

// The panel's boardLines(): one describeLine row per event, from a main line.
function boardLinesOf(state) {
  const seen = new Set();
  const rows = [];
  for (const each of Object.values(state.lines)) {
    if (each.isAlt || seen.has(each.eventId)) continue;
    seen.add(each.eventId);
    rows.push(feed.describeLine(each, state));
  }
  return rows;
}

// A normalised bet record (bets.js contract) on the 49ers -14.5.
function bet(id, overrides = {}) {
  return {
    id, source: "novig_api", venue: "novig", league: "nfl", status: "open",
    eventStart: iso(KICKOFF), eventDate: null, awayTeam: "Seattle Seahawks", homeTeam: "San Francisco 49ers",
    awayKey: "nfl:11", homeKey: "nfl:22", betType: "spread", period: "FG", side: "home", points: -14.5,
    price: 208, stake: 150, toWin: 312, placedAt: iso(T0 + 2 * MIN), closedAt: null, unmatchable: null, approx: [],
    isParlayLeg: false, raw: { probability: 0.3247 },
    ...overrides,
  };
}

// history: line key -> [{at, ...line fields}] observations, through edgemove.observe.
function historyOf(observations) {
  const history = {};
  for (const [key, entries] of Object.entries(observations)) {
    for (const { at, ...fields } of entries) {
      const [, bookId] = /:ms(\d+):/.exec(key);
      const sideIndex = key.includes(":si1:") ? 1 : 0;
      edgemove.observe(history, { ...line(Number(bookId), { sideIndex }), key, ...fields }, { at, source: "snapshot" });
    }
  }
  return history;
}

function capture({ records, lines, history, observingSince = WATCHING_SINCE, now = T0 + 6 * MIN, skipIds = new Set() }) {
  const state = boardState(lines);
  return captureFillFairs({ records, skipIds, state, boardLines: boardLinesOf(state), history, observingSince, now });
}

test("saves the fair the bet's own line showed at the fill, not the one it shows now", () => {
  const history = historyOf({
    [NINERS_ALT(NOVIG)]: [
      { at: T0, price: 208, sourceFormat: 4, sourcePrice: 0.3247, bacr: 178 },
      { at: T0 + 5 * MIN, price: 223, sourceFormat: 4, sourcePrice: 0.3096, bacr: 188 },
    ],
  });
  const result = capture({ records: [bet("novig:1")], lines: [line(NOVIG)], history });
  assert.deepEqual(result.refusals, []);
  assert.deepEqual(result.saves, [{
    betId: "novig:1", lineKey: NINERS_ALT(NOVIG), points: -14.5, fairAmerican: 178, fairObservedAt: T0, placedAt: iso(T0 + 2 * MIN),
  }]);
});

test("prefers the bet's own venue; a venue off the board reads the most recently changed book at the number", () => {
  const history = historyOf({
    [NINERS_ALT(NOVIG)]: [{ at: T0, bacr: 178 }],
    [NINERS_ALT(KALSHI)]: [{ at: T0, bacr: 179 }],
    [NINERS_ALT(DRAFTKINGS)]: [{ at: T0, bacr: 180 }],
  });
  const lines = [
    line(NOVIG, { sequenceNumber: T0 - 60 * MIN }),
    line(KALSHI, { sequenceNumber: T0 - MIN }),
    line(DRAFTKINGS, { sequenceNumber: T0 - 30 * MIN }),
  ];
  const [own] = capture({ records: [bet("novig:1")], lines, history }).saves;
  assert.equal(own.lineKey, NINERS_ALT(NOVIG));
  assert.equal(own.fairAmerican, 178);
  // BFA is not on Unabated's board: the freshest line at -14.5 (Kalshi, changed a minute ago).
  const [other] = capture({ records: [bet("bfa:1", { venue: "bfa", raw: {} })], lines, history }).saves;
  assert.equal(other.lineKey, NINERS_ALT(KALSHI));
  assert.equal(other.fairAmerican, 179);
});

test("the own venue's line first seen after the fill falls back to a book that was watched through it", () => {
  const history = historyOf({
    [NINERS_ALT(NOVIG)]: [{ at: T0 + 4 * MIN, bacr: 188 }],
    [NINERS_ALT(KALSHI)]: [{ at: T0 - 20 * MIN, bacr: 180 }, { at: T0 + MIN, bacr: 178 }, { at: T0 + 5 * MIN, bacr: 188 }],
  });
  const [save] = capture({ records: [bet("novig:1")], lines: [line(NOVIG), line(KALSHI)], history }).saves;
  assert.equal(save.lineKey, NINERS_ALT(KALSHI));
  assert.equal(save.fairAmerican, 178); // the newest observation at or before the fill
  assert.equal(save.fairObservedAt, T0 + MIN);
});

test("refusals are final and say why: before the panel was watching, first seen after, no fair, no or future time", () => {
  const history = historyOf({
    [NINERS_ALT(NOVIG)]: [{ at: T0 + 3 * MIN, bacr: 188 }],
    [SEAHAWKS_ALT(NOVIG)]: [{ at: T0, bacr: null }],
  });
  const lines = [line(NOVIG), line(NOVIG, { sideIndex: 0, bacr: null })];
  const records = [
    bet("novig:before", { placedAt: iso(WATCHING_SINCE - MIN) }),
    bet("novig:unseen"),
    bet("novig:nofair", { side: "away", points: 14.5 }),
    bet("novig:future", { placedAt: iso(T0 + 7 * MIN) }),
    bet("novig:notime", { placedAt: null }),
  ];
  const result = capture({ records, lines, history });
  assert.deepEqual(result.saves, []);
  assert.deepEqual(result.refusals.sort((a, b) => a.betId.localeCompare(b.betId)), [
    { betId: "novig:before", reason: REFUSED.notWatching },
    { betId: "novig:future", reason: REFUSED.future },
    { betId: "novig:nofair", reason: REFUSED.noFair },
    { betId: "novig:notime", reason: REFUSED.noPlacedAt },
    { betId: "novig:unseen", reason: REFUSED.notSeen },
  ]);
  assert.equal(REFUSED.notWatching, "placed before the panel was watching");
});

test("a bet with no line on the board yet is neither saved nor refused, so the next pass looks again", () => {
  const history = historyOf({ [NINERS_ALT(NOVIG)]: [{ at: T0, bacr: 178 }] });
  const records = [
    bet("novig:other-game", { awayTeam: "Dallas Cowboys", homeTeam: "New York Giants", awayKey: "nfl:33", homeKey: "nfl:44" }),
    bet("novig:other-number", { points: -15.5 }),
  ];
  assert.deepEqual(capture({ records, lines: [line(NOVIG)], history }), { saves: [], refusals: [] });
});

test("the other side of the number is not the bet's line", () => {
  const history = historyOf({
    [NINERS_ALT(NOVIG)]: [{ at: T0, bacr: 178 }],
    [SEAHAWKS_ALT(NOVIG)]: [{ at: T0, bacr: -210 }],
  });
  const lines = [line(NOVIG), line(NOVIG, { sideIndex: 0, bacr: -210 })];
  const [save] = capture({ records: [bet("novig:seahawks", { side: "away", points: 14.5 })], lines, history }).saves;
  assert.equal(save.lineKey, SEAHAWKS_ALT(NOVIG));
  assert.equal(save.fairAmerican, -210);
});

test("nothing is read without a gap-free history, and saved, decided, settled or unmatchable bets are skipped", () => {
  const history = historyOf({ [NINERS_ALT(NOVIG)]: [{ at: T0, bacr: 178 }] });
  const lines = [line(NOVIG)];
  assert.deepEqual(capture({ records: [bet("novig:1")], lines, history, observingSince: null }), { saves: [], refusals: [] });
  const records = [
    bet("novig:saved"),
    bet("novig:won", { status: "won" }),
    bet("novig:odd", { unmatchable: "unreadable Novig order" }),
  ];
  assert.deepEqual(capture({ records, lines, history, skipIds: new Set(["novig:saved"]) }), { saves: [], refusals: [] });
});

test("fillTimeEntry: the newest observation at or before the fill", () => {
  const entries = [{ at: T0, bacr: 170 }, { at: T0 + MIN, bacr: 178 }, { at: T0 + 2 * MIN, bacr: 188 }];
  assert.equal(fillTimeEntry(entries, T0 + MIN).entry.bacr, 178);
  assert.equal(fillTimeEntry(entries, T0 + 90 * 1000).entry.bacr, 178);
  assert.deepEqual(fillTimeEntry(entries, T0 - 1), { reason: REFUSED.notSeen });
  assert.deepEqual(fillTimeEntry(undefined, T0), { reason: REFUSED.notSeen });
  assert.deepEqual(fillTimeEntry([{ at: T0, bacr: 50 }], T0), { reason: REFUSED.noFair });
});

// ---- display ----------------------------------------------------------------

test("fairsByBetId keeps well-formed rows only", () => {
  const index = fairsByBetId([
    { betId: "novig:1", fairAmerican: 178 }, { betId: "novig:2", fairAmerican: 50 }, { fairAmerican: 178 }, null,
  ]);
  assert.deepEqual(Array.from(index.keys()), ["novig:1"]);
  assert.equal(fairsByBetId(undefined).size, 0);
});

test("baselineOf: the earliest open bet on this very line that has a saved fair", () => {
  const first = bet("novig:1", { placedAt: iso(T0 + 2 * MIN) });
  const second = bet("novig:2", { placedAt: iso(T0 + 40 * MIN), price: 217 });
  const unsaved = bet("novig:0", { placedAt: iso(T0) });
  const otherNumber = bet("novig:3", { placedAt: iso(T0 - 60 * MIN), points: -13.5 });
  const fairs = fairsByBetId([
    { betId: "novig:1", fairAmerican: 178 }, { betId: "novig:2", fairAmerican: 183 }, { betId: "novig:3", fairAmerican: 165 },
  ]);
  const matches = [
    { tier: "same_line", bet: second }, { tier: "same_line", bet: unsaved }, { tier: "same_line", bet: first },
    { tier: "same_side", bet: otherNumber },
  ];
  const baseline = baselineOf(matches, fairs);
  assert.equal(baseline.bet.id, "novig:1");
  assert.equal(baseline.fair.fairAmerican, 178);
  assert.equal(baseline.placedMs, T0 + 2 * MIN);
  assert.equal(baselineOf([{ tier: "same_side", bet: otherNumber }], fairs), null);
  assert.equal(baselineOf(undefined, fairs), null);
});

test("the 49ers card: the price improved to +223 while the fair fell from 36.0% to 34.7% since the +208 fill: red", () => {
  const baseline = { bet: bet("novig:1"), fair: { fairAmerican: 178 }, placedMs: T0 + 2 * MIN };
  const now = { price: 223, sourceFormat: 4, sourcePrice: 0.3096, bacr: 188 };
  const move = moveSinceFill(baseline, now);
  assert.equal(move.kind, "fair_against");
  assert.equal(Math.round(move.fairDelta * 1000) / 10, -1.2); // 34.72% - 35.97%: the displays round to 34.7 / 36.0
  assert.equal(Math.round(move.priceDelta * 1000) / 10, 1.5);
  assert.deepEqual(move.from, { price: 208, sourceFormat: 4, sourcePrice: 0.3247, bacr: 178 });
  assert.equal(move.to, now);
  // The same price with the fair holding is the book moving away; the fair rising is green.
  assert.equal(moveSinceFill(baseline, { ...now, bacr: 178 }).kind, "book_away");
  assert.equal(moveSinceFill(baseline, { ...now, bacr: 170 }).kind, "fair_to_you");
  // The same rule as the ten-minute tag, applied to the fill and now.
  assert.deepEqual(edgemove.classifyMove(move.from, now), { kind: move.kind, fairDelta: move.fairDelta, priceDelta: move.priceDelta });
});

test("a fill with no price on its record is decided on the fair alone", () => {
  const baseline = { bet: bet("bfa:1", { venue: "bfa", price: null, raw: {} }), fair: { fairAmerican: 178 }, placedMs: T0 };
  const flat = moveSinceFill(baseline, { price: 223, sourceFormat: 1, sourcePrice: null, bacr: 178 });
  assert.equal(flat.kind, "none");
  assert.equal(flat.priceDelta, null);
  assert.equal(moveSinceFill(baseline, { price: 223, sourceFormat: 1, sourcePrice: null, bacr: 188 }).kind, "fair_against");
});

test("fillPriceOf: an exchange fill is its exact probability, a sportsbook fill its American price", () => {
  assert.deepEqual(fillPriceOf(bet("kalshi:x:yes", { venue: "kalshi", price: 355, raw: { vwapCents: 22 } })),
    { price: 355, sourceFormat: 4, sourcePrice: 0.22 });
  assert.deepEqual(fillPriceOf(bet("novig:1")), { price: 208, sourceFormat: 4, sourcePrice: 0.3247 });
  assert.deepEqual(fillPriceOf(bet("betonline:1", { venue: "betonline", price: -110, raw: {} })),
    { price: -110, sourceFormat: 1, sourcePrice: null });
  assert.deepEqual(fillPriceOf(bet("kalshi:y:yes", { venue: "kalshi", price: 355, raw: null })),
    { price: 355, sourceFormat: 1, sourcePrice: null });
  assert.equal(fillPriceOf(bet("bfa:parlay", { venue: "bfa", price: null })), null);
});
