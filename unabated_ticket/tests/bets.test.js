// Run: node --test unabated_ticket/tests
//
// Fixtures:
//   fixtures/bets/kalshi_fixture.json — a read-only pull of the account's own
//   Kalshi fills, unsettled positions and the public market/event payloads
//   (2026-09-11), one or two fills per in-scope series: KXNCAAF1HTOTAL (a NO
//   at position_fp -500 and a YES), KXNCAAFTOTAL (a YES open, a NO settled),
//   KXNCAAFSPREAD (a YES open, the plan's TXAM39 example settled), KXNFLGAME
//   and KXNCAAFGAME moneylines, plus KXNFLOROTY (future), KXMLBRFI (bot
//   series) and KXNEXTTEAMNFL (sell fills, position sold to 0). Ids stripped.
//   The account has no NO moneyline and no NFL spread/total fills, so those
//   are synthetic fills on the same grammar (event/market payloads mirror the
//   real ones; KXNFLSPREAD/KXNFL1HTOTAL grammar read off the public markets
//   endpoint on 2026-09-11).
//   fixtures/v2_slice.json — the NFL feed slice (event 125807, Chicago Bears @
//   Carolina Panthers, 2026-09-13T17:00Z) whose rows, described by
//   feed.describeLine, are the real lines the matcher sees.
const test = require("node:test");
const assert = require("node:assert/strict");
const fs = require("node:fs");
const path = require("node:path");
const feed = require("../extension/feed.js");
const bets = require("../extension/bets.js");
const teams = require("../extension/teams.js");
// The runtime team index the panel builds from Unabated's snapshots, from a
// captured copy (fixtures/teams_index.json) — keys are "<league>:<Unabated id>".
teams.loadIndex(JSON.parse(fs.readFileSync(path.join(__dirname, "fixtures", "teams_index.json"), "utf8")).leagues);
const key = (league, name) => teams.teamKey(league, name);

const FETCHED_AT = "2026-09-11T21:00:00Z";
const NOW = Date.parse(FETCHED_AT);
const NFL = 1;

function kalshiFixture() {
  return JSON.parse(fs.readFileSync(path.join(__dirname, "fixtures", "bets", "kalshi_fixture.json"), "utf8"));
}

function normalizedFixture(extra) {
  const fx = kalshiFixture();
  return bets.normalizeKalshi({
    fills: fx.fills.concat(extra && extra.fills ? extra.fills : []),
    positions: fx.positions.concat(extra && extra.positions ? extra.positions : []),
    markets: Object.assign({}, fx.markets, extra && extra.markets),
    events: Object.assign({}, fx.events, extra && extra.events),
    fetchedAt: FETCHED_AT,
  });
}

function byId(records, id) {
  const record = records.find((r) => r.id === id);
  assert.ok(record, `record ${id} missing`);
  return record;
}

// One NFL slice row per (bet type, period, side, points), described like the
// panel does — books post different spread and total numbers for the same game.
function nflRows() {
  const snapshot = JSON.parse(fs.readFileSync(path.join(__dirname, "fixtures", "v2_slice.json"), "utf8"));
  const state = feed.mergeStates([feed.parseSnapshot(snapshot, { leagueId: NFL })]);
  const rows = {};
  for (const line of Object.values(state.lines)) {
    if (line.isAlt) continue;
    const row = feed.describeLine(line, state);
    rows[`${row.betType} ${row.period} ${row.sideIndex} ${row.points}`] = row;
  }
  return rows;
}

// A describeLine-shaped row for a game not in the NFL slice.
function cfbRow(overrides) {
  return Object.assign({
    league: "cfb", eventId: 900001, awayTeam: "Chattanooga", homeTeam: "Eastern Kentucky",
    eventStart: "2026-09-12T20:00:00.000Z", eventStartMs: Date.parse("2026-09-12T20:00:00Z"),
    betType: "Spread", period: "FG", sideIndex: 0, points: -5.5, rotation: null,
  }, overrides);
}

function fill(ticker, side, countFp, yesPrice, createdTime, action) {
  return {
    ticker, market_ticker: ticker, side, action: action || "buy", count_fp: String(countFp),
    yes_price_dollars: yesPrice.toFixed(4), no_price_dollars: (1 - yesPrice).toFixed(4),
    created_time: createdTime, ts: Math.floor(Date.parse(createdTime) / 1000), is_taker: true,
    fee_cost: "0.000000", exchange_index: 0, book_side: "bid", outcome_side: side,
  };
}

function position(ticker, positionFp, totalTraded) {
  return {
    ticker, position_fp: positionFp.toFixed(2), total_traded_dollars: totalTraded.toFixed(6),
    market_exposure_dollars: "0.000000", realized_pnl_dollars: "0.000000", fees_paid_dollars: "0.000000",
    exchange_index: 0, last_updated_ts: "2026-09-11T15:05:04.400227Z",
  };
}

// Synthetic Kalshi payloads for the slice's game on the verified grammar.
const CHICAR_EVENT_ML = { event_ticker: "KXNFLGAME-26SEP13CHICAR", series_ticker: "KXNFLGAME", title: "CHI Bears vs CAR Panthers", sub_title: "CHI vs CAR (Sep 13)", mutually_exclusive: true };
const CHICAR_EVENT_SPREAD = { event_ticker: "KXNFLSPREAD-26SEP13CHICAR", series_ticker: "KXNFLSPREAD", title: "Chicago vs Carolina: Spread", sub_title: "CHI vs CAR (Sep 13)", mutually_exclusive: false };
const CHICAR_EVENT_1H = { event_ticker: "KXNFL1HTOTAL-26SEP13CHICAR", series_ticker: "KXNFL1HTOTAL", title: "CHI Bears vs CAR Panthers: 1st Half Total", sub_title: "CHI vs CAR (Sep 13)", mutually_exclusive: false };
function market(ticker, eventTicker, strikeType, floorStrike, title) {
  return { ticker, event_ticker: eventTicker, title, yes_sub_title: title, strike_type: strikeType, floor_strike: floorStrike,
    custom_strike: null, status: "active", result: "", expected_expiration_time: "2026-09-13T21:00:00Z", close_time: "2026-09-15T17:00:00Z", rules_primary: "" };
}
const CHICAR = {
  fills: [
    fill("KXNFLGAME-26SEP13CHICAR-CAR", "yes", 100, 0.58, "2026-09-11T15:05:04Z"),
    fill("KXNFLGAME-26SEP13CHICAR-CAR", "no", 100, 0.55, "2026-09-11T15:06:00Z"),
    fill("KXNFLSPREAD-26SEP13CHICAR-CHI14", "yes", 200, 0.40, "2026-09-11T15:07:00Z"),
    fill("KXNFL1HTOTAL-26SEP13CHICAR-21", "yes", 50, 0.50, "2026-09-11T15:08:00Z"),
    fill("KXNFL1HTOTAL-26SEP13CHICAR-21", "no", 60, 0.50, "2026-09-11T15:09:00Z"),
  ],
  positions: [
    position("KXNFLSPREAD-26SEP13CHICAR-CHI14", 200, 80),
    position("KXNFL1HTOTAL-26SEP13CHICAR-21", 50, 55),
  ],
  markets: {
    "KXNFLGAME-26SEP13CHICAR-CAR": market("KXNFLGAME-26SEP13CHICAR-CAR", "KXNFLGAME-26SEP13CHICAR", "structured", null, "Carolina wins"),
    "KXNFLSPREAD-26SEP13CHICAR-CHI14": market("KXNFLSPREAD-26SEP13CHICAR-CHI14", "KXNFLSPREAD-26SEP13CHICAR", "greater", 13.5, "Chicago wins by over 13.5 points?"),
    "KXNFL1HTOTAL-26SEP13CHICAR-21": market("KXNFL1HTOTAL-26SEP13CHICAR-21", "KXNFL1HTOTAL-26SEP13CHICAR", "greater", 20.5, "Will there be over 20.5 1H points scored?"),
  },
  events: {
    "KXNFLGAME-26SEP13CHICAR": CHICAR_EVENT_ML,
    "KXNFLSPREAD-26SEP13CHICAR": CHICAR_EVENT_SPREAD,
    "KXNFL1HTOTAL-26SEP13CHICAR": CHICAR_EVENT_1H,
  },
};

// ---- normalizeKalshi ---------------------------------------------------------

test("normalize: YES spread = the named team at -floor_strike, stake = position x VWAP", () => {
  const record = byId(normalizedFixture(), "kalshi:KXNCAAFSPREAD-26SEP12CHATEKY-CHAT6:yes");
  assert.equal(record.source, "kalshi_api");
  assert.equal(record.venue, "kalshi");
  assert.equal(record.league, "cfb");
  assert.equal(record.betType, "spread");
  assert.equal(record.period, "FG");
  assert.equal(record.awayTeam, "Chattanooga");
  assert.equal(record.homeTeam, "Eastern Kentucky");
  assert.equal(record.awayKey, key("cfb", "Chattanooga"));
  assert.equal(record.homeKey, key("cfb", "Eastern Kentucky"));
  assert.equal(record.side, "away");
  assert.equal(record.points, -5.5);
  assert.equal(record.price, 138); // two buys at 42c
  assert.equal(record.contracts, 400); // position_fp, not the fills
  assert.equal(record.stake, 168);
  assert.equal(record.toWin, 232);
  assert.equal(record.placedAt, "2026-09-11T17:14:55.910367Z");
  assert.equal(record.status, "open");
  assert.equal(record.eventStart, null); // football suffix carries the date only
  assert.equal(record.eventDate, "2026-09-12");
  assert.deepEqual(record.approx, []);
  assert.equal(record.unmatchable, null);
  assert.equal(record.sourceFetchedAt, FETCHED_AT);
  assert.equal(record.raw.ticker, "KXNCAAFSPREAD-26SEP12CHATEKY-CHAT6");
  // #118: the event ticker's suffix is what Unabated's Kalshi rungs carry in sourceKey.
  assert.deepEqual(record.venueIds, { marketTicker: "KXNCAAFSPREAD-26SEP12CHATEKY-CHAT6", eventTicker: "KXNCAAFSPREAD-26SEP12CHATEKY" });
  // No market payload: the ticker is known, the event ticker is not guessed from it.
  const [noMarket] = bets.normalizeKalshi({ fills: [fill("KXNFLGAME-26SEP20PITNE-NE", "yes", 10, 0.5, "2026-09-11T12:00:00Z")], positions: [], markets: {}, events: {}, fetchedAt: FETCHED_AT });
  assert.deepEqual(noMarket.venueIds, { marketTicker: "KXNFLGAME-26SEP20PITNE-NE", eventTicker: null });
});

test("normalize: the plan's TXAM39 example — YES on the home code is home -38.5, settled won", () => {
  const record = byId(normalizedFixture(), "kalshi:KXNCAAFSPREAD-26SEP05MOSUTXAM-TXAM39:yes");
  assert.equal(record.awayTeam, "Missouri St.");
  assert.equal(record.homeTeam, "Texas A&M");
  assert.equal(record.side, "home");
  assert.equal(record.points, -38.5);
  assert.equal(record.price, -108);
  assert.equal(record.stake, 196.56);
  assert.equal(record.status, "won");
  assert.equal(record.closedAt, "2026-09-06T02:00:00Z");
});

test("normalize: NO spread = the other team at +floor_strike (both signs covered)", () => {
  const fx = kalshiFixture();
  const noFills = [
    fill("KXNCAAFSPREAD-26SEP12CHATEKY-CHAT6", "no", 10, 0.42, "2026-09-11T18:00:00Z"),
    fill("KXNCAAFSPREAD-26SEP05MOSUTXAM-TXAM39", "no", 10, 0.52, "2026-09-02T23:00:00Z"),
  ];
  const records = bets.normalizeKalshi({ fills: noFills, positions: [position("KXNCAAFSPREAD-26SEP12CHATEKY-CHAT6", -10, 5.8)], markets: fx.markets, events: fx.events, fetchedAt: FETCHED_AT });
  const chat = byId(records, "kalshi:KXNCAAFSPREAD-26SEP12CHATEKY-CHAT6:no");
  assert.equal(chat.side, "home");
  assert.equal(chat.points, 5.5);
  assert.equal(chat.price, -138); // NO bought at 58c
  assert.equal(chat.status, "open");
  const txam = byId(records, "kalshi:KXNCAAFSPREAD-26SEP05MOSUTXAM-TXAM39:no");
  assert.equal(txam.side, "away");
  assert.equal(txam.points, 38.5);
  assert.equal(txam.status, "lost");
});

test("normalize: totals — YES = Over floor_strike, NO = Under; 1H from the series", () => {
  const records = normalizedFixture();
  const over = byId(records, "kalshi:KXNCAAFTOTAL-26SEP12RICEND-60:yes");
  assert.equal(over.betType, "total");
  assert.equal(over.period, "FG");
  assert.equal(over.side, "over");
  assert.equal(over.points, 59.5);
  assert.equal(over.price, 213); // 32c
  assert.equal(over.stake, 160);
  assert.equal(over.status, "open");
  const under1h = byId(records, "kalshi:KXNCAAF1HTOTAL-26SEP12OKLAMICH-23:no");
  assert.equal(under1h.period, "1H");
  assert.equal(under1h.side, "under");
  assert.equal(under1h.points, 22.5);
  assert.equal(under1h.price, -104); // NO at 51c
  assert.equal(under1h.contracts, 500); // position_fp -500
  assert.equal(under1h.stake, 255);
  assert.equal(under1h.awayTeam, "Oklahoma");
  assert.equal(under1h.homeTeam, "Michigan");
  const underSettled = byId(records, "kalshi:KXNCAAFTOTAL-26SEP05DRKEMONT-62:no");
  assert.equal(underSettled.side, "under");
  assert.equal(underSettled.points, 61.5);
  assert.equal(underSettled.price, -103); // VWAP of 45 @ 50c and 200 @ 51c
  assert.equal(underSettled.contracts, 245);
  assert.equal(underSettled.status, "won");
});

test("normalize: moneylines — YES is the named team; NFL event title 'PIT Steelers vs NE Patriots'", () => {
  const records = normalizedFixture();
  const ne = byId(records, "kalshi:KXNFLGAME-26SEP20PITNE-NE:yes");
  assert.equal(ne.league, "nfl");
  assert.equal(ne.betType, "moneyline");
  assert.equal(ne.awayTeam, "PIT Steelers");
  assert.equal(ne.homeTeam, "NE Patriots");
  assert.equal(ne.awayKey, key("nfl", "Pittsburgh Steelers"));
  assert.equal(ne.homeKey, key("nfl", "New England Patriots"));
  assert.equal(ne.side, "home");
  assert.equal(ne.points, null);
  assert.equal(ne.price, -178); // VWAP of two fills at 64c
  assert.equal(ne.contracts, 550);
  assert.equal(ne.stake, 352);
  assert.equal(ne.eventDate, "2026-09-20");
  assert.deepEqual(ne.approx, []);
  const ecu = byId(records, "kalshi:KXNCAAFGAME-26SEP05ECUALA-ECU:yes");
  assert.equal(ecu.side, "away");
  assert.equal(ecu.status, "lost");
  assert.equal(ecu.price, 3233);
});

test("normalize: NO moneyline is the other team with the tie caveat", () => {
  const record = byId(normalizedFixture(CHICAR), "kalshi:KXNFLGAME-26SEP13CHICAR-CAR:no");
  assert.equal(record.side, "away");
  assert.equal(record.awayTeam, "CHI Bears");
  assert.equal(record.awayKey, key("nfl", "Chicago Bears"));
  assert.equal(record.points, null);
  assert.deepEqual(record.approx, [bets.TIE_CAVEAT]);
  assert.equal(record.price, 122); // NO at 45c
  assert.equal(bets.describeBet(record), "NO CAR Panthers ≈ CHI Bears or tie +122 · 45.0¢");
});

test("normalize: position_fp 0 with trades is closed; sell-only fills price from the sells", () => {
  const record = byId(normalizedFixture(), "kalshi:KXNEXTTEAMNFL-26MCROSBY-LV:yes");
  assert.equal(record.status, "closed");
  assert.equal(record.contracts, 0);
  assert.equal(record.closedAt, "2026-08-24T11:50:26.544087Z");
  assert.equal(record.unmatchable, "not a game market");
  assert.equal(record.price, -285); // sold at 74c
});

test("normalize: futures are 'other' / not a game market; unknown series fails closed", () => {
  const unknown = [fill("KXNBAGAME-26OCT20LALBOS-BOS", "yes", 10, 0.50, "2026-09-11T12:00:00Z")];
  const records = normalizedFixture({ fills: unknown });
  const future = byId(records, "kalshi:KXNFLOROTY-27-MWAS:yes");
  assert.equal(future.betType, "other");
  assert.equal(future.league, null);
  assert.equal(future.side, null);
  assert.equal(future.unmatchable, "not a game market");
  assert.equal(future.status, "open");
  assert.equal(future.stake, 100); // 5000 x 2c
  assert.equal(bets.describeBet(future), "Will Mike Washington Jr. win the Offensive Rookie of the Year? +4900 · 2.0¢");
  const nba = byId(records, "kalshi:KXNBAGAME-26OCT20LALBOS-BOS:yes");
  assert.equal(nba.betType, "other");
  assert.equal(nba.league, null);
  assert.equal(nba.unmatchable, "unknown Kalshi series");
});

test("normalize: a known series with a strike its event cannot name fails closed", () => {
  const fx = kalshiFixture();
  const fills = [fill("KXNFLGAME-26SEP20PITNE-XYZ", "yes", 10, 0.50, "2026-09-11T12:00:00Z")];
  const markets = { "KXNFLGAME-26SEP20PITNE-XYZ": Object.assign({}, fx.markets["KXNFLGAME-26SEP20PITNE-NE"], { ticker: "KXNFLGAME-26SEP20PITNE-XYZ" }) };
  const [record] = bets.normalizeKalshi({ fills, positions: [], markets, events: fx.events, fetchedAt: FETCHED_AT });
  assert.equal(record.unmatchable, "unreadable Kalshi market (strike XYZ not in PIT/NE)");
  assert.equal(record.league, null);
});

test("normalize: MLB suffix carries HHMM Eastern -> eventStart UTC; RFI is an I1 total at 0.5", () => {
  const record = byId(normalizedFixture(), "kalshi:KXMLBRFI-26SEP062210WSHLAD:yes");
  assert.equal(record.league, "mlb");
  assert.equal(record.period, "I1");
  assert.equal(record.side, "over");
  assert.equal(record.points, 0.5);
  assert.equal(record.eventDate, "2026-09-06");
  assert.equal(record.eventStart, "2026-09-07T02:10:00.000Z"); // 10:10 PM EDT
  assert.equal(record.status, "lost");
  assert.equal(record.awayKey, key("mlb", "Washington Nationals")); // "WSH Nationals" via the Kalshi event title
  assert.equal(record.homeKey, key("mlb", "Los Angeles Dodgers"));
});

test("parseEventSuffix: football date only, MLB with time, doubleheader marker kept apart", () => {
  assert.deepEqual(bets.parseEventSuffix("KXNCAAFSPREAD-26SEP12CHATEKY"), { eventDate: "2026-09-12", eventStart: null, gameNumber: null });
  assert.deepEqual(bets.parseEventSuffix("KXMLBRFI-26SEP062210WSHLAD"), { eventDate: "2026-09-06", eventStart: "2026-09-07T02:10:00.000Z", gameNumber: null });
  assert.equal(bets.parseEventSuffix("KXMLBGAME-26SEP041410DETCLEG2").gameNumber, 2);
  assert.equal(bets.parseEventSuffix("KXMLBGAME-26JAN151300DETCLE").eventStart, "2026-01-15T18:00:00.000Z"); // EST
  assert.equal(bets.parseEventSuffix("KXNFLOROTY-27"), null);
});

test("centsToAmerican: both sides of even money", () => {
  assert.equal(bets.centsToAmerican(32), 213);
  assert.equal(bets.centsToAmerican(64), -178);
  assert.equal(bets.centsToAmerican(50), -100);
  assert.equal(bets.centsToAmerican(0), null);
  assert.equal(bets.centsToAmerican(100), null);
});

// ---- tiers on the real NFL slice rows -------------------------------------------

test("NFL slice: moneyline YES on Carolina is same_line on the Carolina row, opposite on Chicago", () => {
  const rows = nflRows();
  const records = normalizedFixture(CHICAR).filter((r) => r.id === "kalshi:KXNFLGAME-26SEP13CHICAR-CAR:yes");
  const home = bets.matchBets(rows["Moneyline FG 1 null"], records);
  assert.equal(home.matches.length, 1);
  assert.equal(home.matches[0].tier, "same_line");
  assert.equal(home.matches[0].label, "CAR Panthers -138 · 58.0¢ · $58 · Kalshi");
  const away = bets.matchBets(rows["Moneyline FG 0 null"], records);
  assert.equal(away.matches[0].tier, "opposite");
  assert.equal(away.matches[0].label, "CAR Panthers -138 · 58.0¢ · $58 · Kalshi");
  assert.deepEqual(home.unmatched, []);
});

test("NFL slice: NO on the Carolina market matches the Chicago row as same_line with the tie caveat", () => {
  const rows = nflRows();
  const records = normalizedFixture(CHICAR).filter((r) => r.id === "kalshi:KXNFLGAME-26SEP13CHICAR-CAR:no");
  const { matches } = bets.matchBets(rows["Moneyline FG 0 null"], records);
  assert.equal(matches[0].tier, "same_line");
  assert.equal(matches[0].label, "NO CAR Panthers ≈ CHI Bears or tie +122 · 45.0¢ · $45 · Kalshi");
  assert.equal(bets.matchBets(rows["Moneyline FG 1 null"], records).matches[0].tier, "opposite");
});

test("NFL slice: spread Chicago -13.5 is same_side on Chicago -14, opposite at a different number on Carolina +14", () => {
  const rows = nflRows();
  const records = normalizedFixture(CHICAR).filter((r) => r.id === "kalshi:KXNFLSPREAD-26SEP13CHICAR-CHI14:yes");
  const sameSide = bets.matchBets(rows["Spread FG 0 -14"], records).matches[0];
  assert.equal(sameSide.tier, "same_side");
  assert.equal(sameSide.label, "Chicago -13.5 +150 · 40.0¢ · $80 · Kalshi · now -14");
  assert.equal(bets.matchBets(rows["Spread FG 0 -3"], records).matches[0].label, "Chicago -13.5 +150 · 40.0¢ · $80 · Kalshi · now -3");
  const opposite = bets.matchBets(rows["Spread FG 1 14"], records).matches[0];
  assert.equal(opposite.tier, "opposite");
  assert.equal(opposite.label, "Chicago -13.5 +150 · 40.0¢ · $80 · Kalshi · now -14");
  const sameGame = bets.matchBets(rows["Total FG 0 47"], records).matches[0];
  assert.equal(sameGame.tier, "same_game");
  assert.equal(sameGame.label, "Chicago -13.5 +150 · 40.0¢ · $80 · Kalshi");
});

test("NFL slice: 1H total YES Over 20.5 is same_line on the 1H Over row; FG total is only same_game", () => {
  const rows = nflRows();
  const records = normalizedFixture(CHICAR).filter((r) => r.id === "kalshi:KXNFL1HTOTAL-26SEP13CHICAR-21:yes");
  const over = bets.matchBets(rows["Total 1H 0 20.5"], records).matches[0];
  assert.equal(over.tier, "same_line");
  assert.equal(over.label, "1H Over 20.5 -100 · 50.0¢ · $25 · Kalshi");
  assert.equal(bets.matchBets(rows["Total 1H 1 20.5"], records).matches[0].tier, "opposite");
  assert.equal(bets.matchBets(rows["Total FG 0 47"], records).matches[0].tier, "same_game");
});

test("NFL slice: a NO total (sold to the other side by the position) is excluded as closed", () => {
  const rows = nflRows();
  // position_fp +50 means only the YES side is open; the NO fills are a closed leg.
  const records = normalizedFixture(CHICAR).filter((r) => r.id === "kalshi:KXNFL1HTOTAL-26SEP13CHICAR-21:no");
  assert.equal(records[0].status, "closed");
  assert.equal(bets.matchBets(rows["Total 1H 1 20.5"], records).matches.length, 0);
});

test("NFL slice: the real NE bet is a different game and never matches", () => {
  const rows = nflRows();
  const records = normalizedFixture().filter((r) => r.league === "nfl");
  for (const row of Object.values(rows)) assert.equal(bets.matchBets(row, records).matches.length, 0);
});

// ---- tiers on CFB rows for the fixture games ------------------------------------

test("CFB: every tier for the Chattanooga -5.5 bet, exact labels", () => {
  const records = normalizedFixture();
  const sameLine = bets.matchBets(cfbRow({ sideIndex: 0, points: -5.5 }), records).matches;
  assert.equal(sameLine.length, 1);
  assert.equal(sameLine[0].tier, "same_line");
  assert.equal(sameLine[0].label, "Chattanooga -5.5 +138 · 42.0¢ · $168 · Kalshi");
  const sameSide = bets.matchBets(cfbRow({ sideIndex: 0, points: -6.5 }), records).matches[0];
  assert.equal(sameSide.tier, "same_side");
  assert.equal(sameSide.label, "Chattanooga -5.5 +138 · 42.0¢ · $168 · Kalshi · now -6.5");
  const opposite = bets.matchBets(cfbRow({ sideIndex: 1, points: 5.5 }), records).matches[0];
  assert.equal(opposite.tier, "opposite");
  assert.equal(opposite.label, "Chattanooga -5.5 +138 · 42.0¢ · $168 · Kalshi");
  const oppositeOther = bets.matchBets(cfbRow({ sideIndex: 1, points: 6.5 }), records).matches[0];
  assert.equal(oppositeOther.label, "Chattanooga -5.5 +138 · 42.0¢ · $168 · Kalshi · now -6.5");
  const sameGame = bets.matchBets(cfbRow({ betType: "Total", sideIndex: 0, points: 48.5 }), records).matches[0];
  assert.equal(sameGame.tier, "same_game");
  assert.equal(sameGame.label, "Chattanooga -5.5 +138 · 42.0¢ · $168 · Kalshi");
});

test("CFB: NO 1H total is same_line on the 1H Under row and opposite on the 1H Over row", () => {
  const records = normalizedFixture();
  const game = { eventId: 900002, awayTeam: "Oklahoma", homeTeam: "Michigan", betType: "Total", period: "1H", points: 22.5 };
  const under = bets.matchBets(cfbRow(Object.assign({ sideIndex: 1 }, game)), records).matches[0];
  assert.equal(under.tier, "same_line");
  assert.equal(under.label, "1H Under 22.5 -104 · 51.0¢ · $255 · Kalshi");
  const over = bets.matchBets(cfbRow(Object.assign({ sideIndex: 0 }, game)), records).matches[0];
  assert.equal(over.tier, "opposite");
  const fullGame = bets.matchBets(cfbRow(Object.assign({ sideIndex: 1 }, game, { period: "FG", points: 48.5 })), records).matches[0];
  assert.equal(fullGame.tier, "same_game");
  assert.equal(fullGame.label, "1H Under 22.5 -104 · 51.0¢ · $255 · Kalshi");
});

test("CFB: team order on the venue side does not matter", () => {
  const records = normalizedFixture();
  const swapped = cfbRow({ awayTeam: "Eastern Kentucky", homeTeam: "Chattanooga", sideIndex: 1, points: -5.5 });
  assert.equal(bets.matchBets(swapped, records).matches[0].tier, "same_line");
});

test("game match: football bets match by Eastern date within one day; settled bets never match", () => {
  const records = normalizedFixture();
  // 2026-09-13T02:00Z is still Sep 12 in Eastern time.
  assert.equal(bets.matchBets(cfbRow({ eventStartMs: Date.parse("2026-09-13T02:00:00Z") }), records).matches[0].tier, "same_line");
  assert.equal(bets.matchBets(cfbRow({ eventStartMs: Date.parse("2026-09-14T00:30:00Z") }), records).matches[0].tier, "same_line"); // Sep 13 ET
  assert.equal(bets.matchBets(cfbRow({ eventStartMs: Date.parse("2026-09-15T00:30:00Z") }), records).matches.length, 0); // Sep 14 ET
  assert.equal(bets.matchBets(cfbRow({ eventStartMs: null }), records).matches.length, 0);
  const settledGame = cfbRow({ awayTeam: "Missouri St.", homeTeam: "Texas A&M", eventStartMs: Date.parse("2026-09-05T20:00:00Z"), sideIndex: 1, points: -38.5 });
  assert.equal(bets.matchBets(settledGame, records).matches.length, 0);
});

test("game match: a bet with eventStart uses the 30-minute tolerance", () => {
  const withStart = Object.assign(byId(normalizedFixture(), "kalshi:KXNCAAFSPREAD-26SEP12CHATEKY-CHAT6:yes"), { eventStart: "2026-09-12T20:00:00Z" });
  assert.equal(bets.matchBets(cfbRow({ eventStartMs: Date.parse("2026-09-12T20:29:00Z") }), [withStart]).matches.length, 1);
  assert.equal(bets.matchBets(cfbRow({ eventStartMs: Date.parse("2026-09-12T20:31:00Z") }), [withStart]).matches.length, 0);
});

test("game match: two board events for one bet is ambiguous — unmatched, never a guess", () => {
  const records = normalizedFixture();
  const saturday = cfbRow({ eventId: 1 });
  const sunday = cfbRow({ eventId: 2, eventStartMs: Date.parse("2026-09-13T20:00:00Z") });
  const result = bets.matchBets(saturday, records, { lines: [saturday, sunday] });
  assert.equal(result.matches.length, 0);
  assert.equal(result.unmatched.length, 1);
  assert.equal(result.unmatched[0].reason, "ambiguous game");
  assert.equal(result.unmatched[0].bet.id, "kalshi:KXNCAAFSPREAD-26SEP12CHATEKY-CHAT6:yes");
});

test("game match: a closed position on the same market never matches", () => {
  const fx = kalshiFixture();
  const closed = bets.normalizeKalshi({
    fills: fx.fills.filter((f) => f.ticker === "KXNCAAFSPREAD-26SEP12CHATEKY-CHAT6"),
    positions: [position("KXNCAAFSPREAD-26SEP12CHATEKY-CHAT6", 0, 336)],
    markets: fx.markets, events: fx.events, fetchedAt: FETCHED_AT,
  });
  assert.equal(closed[0].status, "closed");
  assert.equal(bets.matchBets(cfbRow(), closed).matches.length, 0);
});

// ---- annotateRows / unmatchedReasons ---------------------------------------------

test("annotateRows: strongest tier per row, with the dollars held and against on the market", () => {
  const records = normalizedFixture();
  const rows = [
    cfbRow({ sideIndex: 0, points: -5.5 }),
    cfbRow({ sideIndex: 0, points: -6.5 }),
    cfbRow({ sideIndex: 1, points: 5.5 }),
    cfbRow({ betType: "Total", sideIndex: 0, points: 48.5 }),
    cfbRow({ eventId: 900003, awayTeam: "Lehigh", homeTeam: "Georgetown" }),
  ];
  const annotated = bets.annotateRows(rows, records);
  assert.deepEqual(annotated.map((a) => a.tier), ["same_line", "same_side", "opposite", "same_game", null]);
  assert.ok(annotated.every((a) => a.exposure && typeof a.exposure.held === "number" && typeof a.exposure.against === "number"));
  assert.ok(annotated[0].exposure.held > 0 && annotated[0].exposure.against === 0);
  assert.ok(annotated[2].exposure.against > 0 && annotated[2].exposure.held === 0);
  assert.deepEqual(annotated[4].exposure, { held: 0, against: 0, heldBets: [], againstBets: [] });
  assert.equal(annotated[0].matches[0].bet.id, "kalshi:KXNCAAFSPREAD-26SEP12CHATEKY-CHAT6:yes");
});

test("exposureOf: held sums same_line + same_side stakes, against sums opposite, same_game adds nothing", () => {
  const bet = (stake) => ({ stake });
  const exposure = bets.exposureOf([
    { tier: "same_line", bet: bet(100) }, { tier: "same_side", bet: bet(50.5) },
    { tier: "opposite", bet: bet(20) }, { tier: "same_game", bet: bet(999) }, { tier: "opposite", bet: bet(null) },
  ]);
  assert.equal(exposure.held, 150.5);
  assert.equal(exposure.against, 20);
  assert.equal(exposure.heldBets.length, 2);
  assert.equal(exposure.againstBets.length, 2);
  assert.deepEqual(bets.exposureOf([]), { held: 0, against: 0, heldBets: [], againstBets: [] });
});

test("unmatchedReasons: every open bet with no match, with why", () => {
  const nbaFill = [fill("KXNBAGAME-26OCT20LALBOS-BOS", "yes", 10, 0.50, "2026-09-11T12:00:00Z")];
  const records = normalizedFixture({ fills: nbaFill }).concat([
    Object.assign({}, byId(normalizedFixture(), "kalshi:KXNCAAFGAME-26SEP12LINWMOSU-LINW:yes"), { id: "kalshi:synthetic:unknown-team", awayTeam: "Springfield", awayKey: null }),
  ]);
  const board = [cfbRow(), cfbRow({ eventId: 1 }), cfbRow({ eventId: 2, eventStartMs: Date.parse("2026-09-13T20:00:00Z") })];
  const reasons = new Map(bets.unmatchedReasons(records, board).map((u) => [u.bet.id, u.reason]));
  assert.equal(reasons.get("kalshi:KXNCAAFSPREAD-26SEP12CHATEKY-CHAT6:yes"), "ambiguous game");
  // A Kalshi bet carries its event suffix, so its reason names both tiers (#118 step 3).
  assert.equal(reasons.get("kalshi:KXNCAAFTOTAL-26SEP12RICEND-60:yes"), "by id: Kalshi event 26SEP12RICEND not on any board ladder; by name: no event on the board yet");
  assert.equal(reasons.get("kalshi:KXNFLGAME-26SEP20PITNE-NE:yes"), "league not on the scanner");
  assert.equal(reasons.get("kalshi:KXNFLOROTY-27-MWAS:yes"), "not a game market");
  assert.equal(reasons.get("kalshi:KXNBAGAME-26OCT20LALBOS-BOS:yes"), "unknown Kalshi series");
  assert.equal(reasons.get("kalshi:synthetic:unknown-team"), "by id: Kalshi event 26SEP12LINWMOSU not on any board ladder; by name: team not recognised (Springfield)");
  assert.equal(reasons.has("kalshi:KXNCAAFSPREAD-26SEP05MOSUTXAM-TXAM39:yes"), false); // settled: not a problem
  assert.equal(reasons.has("kalshi:KXNEXTTEAMNFL-26MCROSBY-LV:yes"), false); // closed
  const unambiguous = new Map(bets.unmatchedReasons(records, [cfbRow()]).map((u) => [u.bet.id, u.reason]));
  assert.equal(unambiguous.has("kalshi:KXNCAAFSPREAD-26SEP12CHATEKY-CHAT6:yes"), false);
});

// ---- retention / dedupe -----------------------------------------------------------

test("pruneForRetention: open bets always, settled and closed within the window", () => {
  const records = normalizedFixture();
  const ids = (list) => list.map((r) => r.id).sort();
  const kept = bets.pruneForRetention(records, NOW, 30);
  assert.deepEqual(ids(kept), ids(records)); // everything settled in the last 30 days
  const tight = bets.pruneForRetention(records, NOW, 3);
  assert.equal(tight.some((r) => r.id === "kalshi:KXNCAAFTOTAL-26SEP05DRKEMONT-62:no"), false); // settled Sep 5
  assert.equal(tight.some((r) => r.id === "kalshi:KXNEXTTEAMNFL-26MCROSBY-LV:yes"), false); // closed Aug 24
  assert.equal(tight.some((r) => r.id === "kalshi:KXNCAAFSPREAD-26SEP12CHATEKY-CHAT6:yes"), true); // open
  const later = bets.pruneForRetention(records, NOW + 40 * 24 * 3600 * 1000);
  assert.equal(later.every((r) => r.status === "open"), true);
  const undated = [{ id: "x", status: "lost", closedAt: null }];
  assert.equal(bets.pruneForRetention(undated, NOW, 1).length, 1);
});

test("resolveTeamKeys: fills null keys from names the way the service leaves them", () => {
  const records = normalizedFixture();
  const fromService = records.map((r) => Object.assign({}, r, { awayKey: null, homeKey: null }));
  const resolved = bets.resolveTeamKeys(fromService);
  assert.deepEqual(resolved, records);
  assert.equal(fromService[0].awayKey, null); // input untouched
  const ne = resolved.find((r) => r.id === "kalshi:KXNFLGAME-26SEP20PITNE-NE:yes");
  assert.equal(ne.awayKey, key("nfl", "Pittsburgh Steelers"));
  assert.equal(ne.homeKey, key("nfl", "New England Patriots"));
  const unknown = bets.resolveTeamKeys([{ league: "cfb", awayTeam: "Springfield", homeTeam: "Lehigh", awayKey: null, homeKey: null }])[0];
  assert.equal(unknown.awayKey, null);
  assert.equal(unknown.homeKey, key("cfb", "Lehigh"));
  const kept = { league: "cfb", awayTeam: "Lehigh", homeTeam: "Drake", awayKey: "cfb:custom", homeKey: null };
  assert.equal(bets.resolveTeamKeys([kept])[0].awayKey, "cfb:custom");
  const future = fromService.find((r) => r.league === null);
  assert.equal(bets.resolveTeamKeys([future])[0], future);
});

test("dedupeByNativeId: newest record per id across consecutive payloads", () => {
  const first = { generatedAt: "2026-09-11T20:00:00Z", bets: [
    { id: "kalshi:a:yes", sourceFetchedAt: "2026-09-11T20:00:00Z", status: "open", contracts: 100 },
    { id: "kalshi:b:yes", sourceFetchedAt: "2026-09-11T20:00:00Z", status: "open", contracts: 5 },
  ] };
  const second = { generatedAt: "2026-09-11T20:01:00Z", bets: [
    { id: "kalshi:a:yes", sourceFetchedAt: "2026-09-11T20:01:00Z", status: "closed", contracts: 0 },
  ] };
  const merged = bets.dedupeByNativeId([first, second]);
  assert.equal(merged.length, 2);
  assert.equal(merged.find((r) => r.id === "kalshi:a:yes").status, "closed");
  assert.equal(merged.find((r) => r.id === "kalshi:b:yes").contracts, 5);
  // An older payload replayed after a newer one does not roll the record back.
  assert.equal(bets.dedupeByNativeId([second, first]).find((r) => r.id === "kalshi:a:yes").status, "closed");
  assert.equal(bets.dedupeByNativeId([first.bets]).length, 2);
});

// ---- Novig (#116): records from the content-script source against the NFL slice ----

const novigSource = require("../extension/novig_bets.js");

// The fixture's NFL rows sit on the slice's game (Chicago Bears @ Carolina
// Panthers, 2026-09-13T17:00Z); team keys resolved the way the panel does.
function novigRecords() {
  const fx = JSON.parse(fs.readFileSync(path.join(__dirname, "fixtures", "bets", "novig_bets.json"), "utf8"));
  const pages = new Map();
  let out = null;
  for (const response of fx.responses) out = novigSource.applyResponse(pages, response);
  return bets.resolveTeamKeys(novigSource.normalizeNovig({ orders: out.orders, parlays: out.parlays, readAt: FETCHED_AT }));
}

test("Novig: moneyline bid on Carolina is same_line on the Carolina row and opposite on Chicago", () => {
  const rows = nflRows();
  const records = novigRecords().filter((r) => r.id === "novig:o-ml-car");
  const home = bets.matchBets(rows["Moneyline FG 1 null"], records);
  assert.equal(home.matches.length, 1);
  assert.equal(home.matches[0].tier, "same_line");
  assert.equal(home.matches[0].label, "Carolina Panthers -138 · 58.0¢ · $58 · Novig");
  const away = bets.matchBets(rows["Moneyline FG 0 null"], records);
  assert.equal(away.matches[0].tier, "opposite");
  assert.equal(away.matches[0].label, "Carolina Panthers -138 · 58.0¢ · $58 · Novig");
});

test("Novig: spread bid Chicago -13.5 is same_side on Chicago -14 and opposite at a different number on Carolina +14", () => {
  const rows = nflRows();
  const records = novigRecords().filter((r) => r.id === "novig:o-sp-chi");
  const sameSide = bets.matchBets(rows["Spread FG 0 -14"], records).matches[0];
  assert.equal(sameSide.tier, "same_side");
  assert.equal(sameSide.label, "Chicago Bears -13.5 +150 · 40.0¢ · $80 · Novig · now -14");
  const opposite = bets.matchBets(rows["Spread FG 1 14"], records).matches[0];
  assert.equal(opposite.tier, "opposite");
  assert.equal(opposite.label, "Chicago Bears -13.5 +150 · 40.0¢ · $80 · Novig · now -14");
});

test("Novig: the LAY of Chicago -13.5 is Carolina +13.5 — same_side on Carolina +14, opposite on Chicago -14", () => {
  const rows = nflRows();
  const records = novigRecords().filter((r) => r.id === "novig:o-sp-lay");
  const sameSide = bets.matchBets(rows["Spread FG 1 14"], records).matches[0];
  assert.equal(sameSide.tier, "same_side");
  assert.equal(sameSide.label, "Carolina Panthers +13.5 -150 · 60.0¢ · $30 · Novig · now +14");
  assert.equal(bets.matchBets(rows["Spread FG 0 -14"], records).matches[0].tier, "opposite");
});

test("formatPlacedAt: Eastern month, day and 12-hour clock", () => {
  assert.equal(bets.formatPlacedAt("2026-09-10T18:15:00Z"), "Sep 10 2:15 PM");
  assert.equal(bets.formatPlacedAt("2026-01-15T18:00:00Z"), "Jan 15 1:00 PM");
  assert.equal(bets.formatPlacedAt("not a date"), "unknown time");
});

test("Novig: a resting Under and a parlay leg flag the game; settled, void and closed orders never match", () => {
  const rows = nflRows();
  const records = novigRecords();
  const resting = bets.matchBets(rows["Total FG 1 47"], records.filter((r) => r.id === "novig:o-tot-under-rest")).matches[0];
  assert.equal(resting.tier, "same_side");
  const leg = bets.matchBets(rows["Moneyline FG 1 null"], records.filter((r) => r.id === "novig:pl-1:0")).matches[0];
  assert.equal(leg.tier, "same_line");
  assert.equal(leg.label, "Carolina Panthers -122 · 55.0¢ (parlay leg) · $25 · Novig");
  for (const id of ["novig:o-won", "novig:o-lay-won", "novig:o-push", "novig:o-cancel", "novig:o-wash", "novig:pl-2:0"]) {
    assert.equal(bets.matchBets(rows["Moneyline FG 1 null"], records.filter((r) => r.id === id)).matches.length, 0, id);
  }
});

test("Novig: unmatched reasons — props, unsupported leagues and unreadable blobs carry the source's reason", () => {
  const reasons = Object.fromEntries(bets.unmatchedReasons(novigRecords(), Object.values(nflRows())).map(({ bet, reason }) => [bet.id, reason]));
  assert.equal(reasons["novig:o-prop"], "not a game market");
  assert.equal(reasons["novig:o-atp"], "league not supported (ATP)");
  assert.equal(reasons["novig:o-noteams"], "unreadable Novig order (no teams on event ev-x)");
  assert.equal(reasons["novig:o-f5"], "league not on the scanner");
  assert.equal(reasons["novig:o-ml-car"], undefined);
});

// ---- BetOnline records (no game date; rotation + placed-time window) --------------

// A record as bets_service/sources/betonline.py emits it (fixture
// betonline_history.json grammar): eventStart / eventDate null, the rotation
// number, approx game_date_unknown, keys null until resolveTeamKeys.
function betonlineRecord(overrides) {
  return Object.assign({
    id: "betonline:996000009-1", source: "betonline_api", venue: "betonline", league: "nfl",
    eventStart: null, eventDate: null, awayTeam: "Chicago Bears", homeTeam: null, awayKey: null, homeKey: null,
    rotation: 465, betType: "spread", period: "FG", side: "away", points: 3, price: -110,
    stake: 110, toWin: 100, contracts: null, placedAt: "2026-09-11T19:05:44Z", status: "open", closedAt: null,
    isParlayLeg: false, parlayId: null, legIndex: null, legCount: null,
    approx: ["game_date_unknown", "side_from_rotation_parity"], unmatchable: null, sourceFetchedAt: FETCHED_AT,
    raw: { description: "Desktop - FOOTBALL - 465 Chicago Bears +3 -110 for GAME " },
  }, overrides);
}

test("BetOnline: a first-half total naming both teams is same_line on the slice's 1H Under row by team pair + placed window", () => {
  const rows = nflRows();
  const bet = bets.resolveTeamKeys([betonlineRecord({
    id: "betonline:996000010-1", homeTeam: "Carolina Panthers", betType: "total", period: "1H", side: "under",
    points: 20.5, price: -110, rotation: 466, approx: ["game_date_unknown"],
  })])[0];
  const under = bets.matchBets(rows["Total 1H 1 20.5"], [bet], { lines: Object.values(rows) });
  assert.equal(under.matches[0].tier, "same_line");
  assert.equal(under.matches[0].label, "1H Under 20.5 -110 · 52.4¢ · $110 · Betonline");
  const over = bets.matchBets(rows["Total 1H 0 20.5"], [bet], { lines: Object.values(rows) });
  assert.equal(over.matches[0].tier, "opposite");
  assert.equal(bets.matchBets(rows["Total FG 1 46"], [bet], { lines: Object.values(rows) }).matches[0].tier, "same_game");
});

test("BetOnline: a one-team spread with unresolved keys matches event 125807 by rotation 465 alone", () => {
  const rows = nflRows();
  const bet = betonlineRecord();  // keys left null: only the rotation can carry it
  const chicago = bets.matchBets(rows["Spread FG 0 -3"], [bet], { lines: Object.values(rows) });
  assert.equal(chicago.matches[0].tier, "same_side");  // Chicago +3 held, row is Chicago -3
  const carolina = bets.matchBets(rows["Spread FG 1 3"], [bet], { lines: Object.values(rows) });
  assert.equal(carolina.matches[0].tier, "opposite");
  assert.deepEqual(bets.unmatchedReasons([bet], Object.values(rows)), []);
});

test("BetOnline: the rotation match is refused when the one resolved team is not in the row's game", () => {
  const rows = nflRows();
  const bet = betonlineRecord({ awayTeam: "Philadelphia Eagles", awayKey: "nfl:phi" });
  assert.equal(bets.matchBets(rows["Spread FG 0 -3"], [bet], { lines: Object.values(rows) }).matches.length, 0);
  assert.deepEqual(bets.unmatchedReasons([bet], Object.values(rows)).map((u) => u.reason), ["no event on the board yet"]);
  // Undecidable when the row's own names are unknown to teams.js: the rotation stands.
  const unknownRow = Object.assign({}, rows["Spread FG 0 -3"], {
    awayTeamId: null, homeTeamId: null, awayTeam: "Chi. (Unabated spelling)", homeTeam: "Car. (Unabated spelling)",
  });
  assert.equal(bets.matchBets(unknownRow, [bet], { lines: [unknownRow] }).matches.length, 1);
});

test("BetOnline: the placed-time window is 12 h back and 14 days ahead of placedAt", () => {
  const rows = nflRows();
  const lines = Object.values(rows);
  const kickoff = Date.parse("2026-09-13T17:00:00Z");
  const at = (ms) => betonlineRecord({ placedAt: new Date(ms).toISOString() });
  assert.equal(bets.matchBets(rows["Spread FG 0 -3"], [at(kickoff - 14 * 24 * 3600 * 1000)], { lines }).matches.length, 1);
  assert.equal(bets.matchBets(rows["Spread FG 0 -3"], [at(kickoff - 15 * 24 * 3600 * 1000)], { lines }).matches.length, 0);
  assert.equal(bets.matchBets(rows["Spread FG 0 -3"], [at(kickoff + 11 * 3600 * 1000)], { lines }).matches.length, 1);
  assert.equal(bets.matchBets(rows["Spread FG 0 -3"], [at(kickoff + 13 * 3600 * 1000)], { lines }).matches.length, 0);
  // Without the flag a dateless record never matches, whatever placedAt says.
  assert.equal(bets.matchBets(rows["Spread FG 0 -3"], [betonlineRecord({ approx: [] })], { lines }).matches.length, 0);
});

test("BetOnline: a matched bet whose team name teams.js does not know is not listed as unmatched", () => {
  const rows = nflRows();
  const bet = bets.resolveTeamKeys([betonlineRecord({ awayTeam: "Chi Bears (BOL spelling)" })])[0];
  assert.equal(bet.awayKey, null);
  assert.deepEqual(bets.unmatchedReasons([bet], Object.values(rows)), []);
  const off = bets.resolveTeamKeys([betonlineRecord({ awayTeam: "Chi Bears (BOL spelling)", rotation: 999 })])[0];
  assert.equal(bets.unmatchedReasons([off], Object.values(rows))[0].reason, "team not recognised (Chi Bears (BOL spelling))");
});

// ---- venue id joins (#118 step 3) ------------------------------------------------
// fixtures/v2_venue_ids_slice.json: real CFB event 123742, Duquesne (927) @
// Washington State (717), 2026-09-19T19:30Z. Its Kalshi rungs carry the
// suffix 26SEP19DUQWSU; Novig's rung on its -35.5 main number carries outcome
// 01a0a05e-…-9e502e4e98d5. Records below spell the teams the way no team
// table knows and leave both keys null: only the id can carry them.

const CFB_LEAGUE = 2;
const DUQ_WSU_SUFFIX = "26SEP19DUQWSU";
const NOVIG_WSU_MINUS_35_5 = "01a0a05e-e864-7b73-9b80-9e502e4e98d5";

function venueState() {
  const snapshot = JSON.parse(fs.readFileSync(path.join(__dirname, "fixtures", "v2_venue_ids_slice.json"), "utf8"));
  return feed.mergeStates([feed.parseSnapshot(snapshot, { leagueId: CFB_LEAGUE })]);
}

// Main-line rows keyed "<betType> <sideIndex> <book id>".
function venueRows(state) {
  const board = state || venueState();
  const rows = {};
  for (const line of Object.values(board.lines)) {
    if (line.isAlt) continue;
    const row = feed.describeLine(line, board);
    rows[`${row.betType} ${row.sideIndex} ${row.book.id}`] = row;
  }
  return rows;
}

// The Kalshi bet as normalizeKalshi builds it from the real grammar, with the
// team names then made unrecognisable.
function kalshiDuqWsu(ticker, side, floorStrike, options) {
  const opts = options || {};
  const seriesEvent = ticker.split("-").slice(0, 2).join("-");
  const [record] = bets.normalizeKalshi({
    fills: [fill(ticker, side, 100, 0.5, "2026-09-15T15:00:00Z")],
    positions: [position(ticker, side === "yes" ? 100 : -100, 50)],
    markets: { [ticker]: market(ticker, seriesEvent, floorStrike == null ? "structured" : "greater", floorStrike, ticker) },
    events: { [seriesEvent]: { event_ticker: seriesEvent, title: "Duquesne vs Washington St.: Spread", sub_title: opts.subTitle || "DUQ vs WSU (Sep 19)" } },
    fetchedAt: FETCHED_AT,
  });
  return Object.assign(record, { awayTeam: "Duquesne Dukes (venue spelling)", homeTeam: "Wazzu (venue spelling)", awayKey: null, homeKey: null },
    "venueIds" in opts ? { venueIds: opts.venueIds } : {});
}

function novigDuqWsu(overrides) {
  return Object.assign({
    id: "novig:o-wsu", source: "novig_api", venue: "novig", rotation: null, price: -110, stake: 55, toWin: 50, contracts: 100,
    placedAt: "2026-09-15T15:00:00Z", status: "open", closedAt: null, isParlayLeg: false, parlayId: null, legIndex: null, legCount: null,
    approx: [], sourceFetchedAt: FETCHED_AT,
    venueIds: { marketId: "m-1", outcomeId: NOVIG_WSU_MINUS_35_5, eventId: "e-1", gameId: "g-1" },
    league: "cfb", eventStart: "2026-09-19T19:30:00.000Z", eventDate: "2026-09-19",
    awayTeam: "Duquesne Dukes (venue spelling)", homeTeam: "Wazzu (venue spelling)", awayKey: null, homeKey: null,
    betType: "spread", period: "FG", side: "home", points: -35.5, unmatchable: null,
  }, overrides);
}

test("id join: a Kalshi spread whose team names resolve nowhere matches its game on the event suffix alone", () => {
  const rows = venueRows();
  const lines = Object.values(rows);
  const bet = kalshiDuqWsu("KXNCAAFSPREAD-26SEP19DUQWSU-WSU36", "yes", 35.5);
  assert.equal(bet.venueIds.eventTicker, "KXNCAAFSPREAD-26SEP19DUQWSU");
  const rowWithoutIds = Object.assign({}, rows["Spread 1 105"], { venueIds: null });
  assert.equal(bets.matchBets(rowWithoutIds, [bet]).matches.length, 0); // no ids on the board: the names alone cannot place it
  // Tier from the bet's own number: WSU -35.5 is the Novig main number, Kalshi's main sits at -38.5.
  assert.equal(bets.matchBets(rows["Spread 1 89"], [bet], { lines }).matches[0].tier, "same_line");
  const kalshiRow = bets.matchBets(rows["Spread 1 105"], [bet], { lines }).matches[0];
  assert.equal(kalshiRow.tier, "same_side");
  assert.ok(kalshiRow.label.endsWith("· now -38.5"));
  assert.equal(bets.matchBets(rows["Spread 0 89"], [bet], { lines }).matches[0].tier, "opposite");
  assert.equal(bets.matchBets(rows["Total 0 105"], [bet], { lines }).matches[0].tier, "same_game");
  assert.deepEqual(bets.unmatchedReasons([bet], lines), []);
  assert.equal(bets.annotateRows(lines, [bet]).filter((flag) => flag.tier).length, lines.length);
  // The NO side of the same contract is Duquesne +35.5: placed off the contract's own strike, not a name.
  const no = kalshiDuqWsu("KXNCAAFSPREAD-26SEP19DUQWSU-WSU36", "no", 35.5);
  assert.deepEqual([no.side, no.points], ["away", 35.5]);
  assert.equal(bets.matchBets(rows["Spread 0 89"], [no], { lines }).matches[0].tier, "same_line");
  assert.equal(bets.matchBets(rows["Spread 1 89"], [no], { lines }).matches[0].tier, "opposite");
});

test("id join: a Kalshi moneyline joins through the suffix its event's spread and total rungs carry", () => {
  const rows = venueRows();
  const lines = Object.values(rows);
  const bet = kalshiDuqWsu("KXNCAAFGAME-26SEP19DUQWSU-WSU", "yes", null);
  assert.equal(bet.betType, "moneyline");
  // Joined, but no rung carries a moneyline contract to say which side an unkeyed bet is on.
  assert.equal(bets.matchBets(rows["Moneyline 1 105"], [bet], { lines }).matches[0].tier, "same_game");
  assert.deepEqual(bets.unmatchedReasons([bet], lines), []);
  // With its teams keyed (one of them is enough to place the side) the tier is exact.
  const keyed = Object.assign({}, bet, { homeKey: "cfb:717" });
  assert.equal(bets.matchBets(rows["Moneyline 1 105"], [keyed], { lines }).matches[0].tier, "same_line");
});

test("id join: a Novig spread joins on its outcome id; a LAY of that outcome is tiered as the other side", () => {
  const rows = venueRows();
  const lines = Object.values(rows);
  const back = novigDuqWsu();
  assert.equal(bets.matchBets(rows["Spread 1 89"], [back], { lines }).matches[0].tier, "same_line");
  assert.equal(bets.matchBets(rows["Spread 0 89"], [back], { lines }).matches[0].tier, "opposite");
  // The lay's outcome id names WSU -35.5, the side it is against: the tier reads the bet, not the id's target.
  const lay = novigDuqWsu({ id: "novig:o-lay", side: "away", points: 35.5 });
  assert.equal(bets.matchBets(rows["Spread 0 89"], [lay], { lines }).matches[0].tier, "same_line");
  assert.equal(bets.matchBets(rows["Spread 1 89"], [lay], { lines }).matches[0].tier, "opposite");
});

test("id join: a main line the changes stream moved off the id's number is same_side, not same_line", () => {
  const state = venueState();
  const main = state.lines["420942308:ms89:si1:tid717"];
  feed.applyChanges(state, { lines: [{ ...main, points: -36.5, price: -105, sequenceNumber: main.sequenceNumber + 1, eventStart: null }] });
  // The id map still says -35.5 (it rebuilds only with the snapshot); the row is at -36.5 now.
  assert.equal(state.events[123742].venueIds.novigOutcomes[NOVIG_WSU_MINUS_35_5].points, -35.5);
  const rows = venueRows(state);
  const match = bets.matchBets(rows["Spread 1 89"], [novigDuqWsu()], { lines: Object.values(rows) }).matches[0];
  assert.equal(match.tier, "same_side");
  assert.ok(match.label.endsWith("· now -36.5"));
});

test("id join: a suffix on two board events is ambiguous — named, never a guess", () => {
  const rows = venueRows();
  const twin = Object.assign({}, rows["Spread 1 105"], { eventId: 999999 });
  const lines = Object.values(rows).concat([twin]);
  const bet = kalshiDuqWsu("KXNCAAFSPREAD-26SEP19DUQWSU-WSU36", "yes", 35.5);
  const result = bets.matchBets(rows["Spread 1 105"], [bet], { lines });
  assert.equal(result.matches.length, 0);
  assert.equal(result.unmatched[0].reason, "ambiguous game (Kalshi event 26SEP19DUQWSU on 2 board events)");
  assert.deepEqual(bets.unmatchedReasons([bet], lines).map((u) => u.reason), ["ambiguous game (Kalshi event 26SEP19DUQWSU on 2 board events)"]);
});

test("id join: when the id and the team names point at different events, the id wins", () => {
  const rows = venueRows();
  // Another event on the same day whose teams ARE the bet's resolved names, and no venue ids.
  const decoy = Object.assign({}, rows["Spread 1 105"], { eventId: 555555, awayTeamId: 11, homeTeamId: 12, venueIds: null });
  const lines = Object.values(rows).concat([decoy]);
  const bet = Object.assign(kalshiDuqWsu("KXNCAAFSPREAD-26SEP19DUQWSU-WSU36", "yes", 35.5), { awayKey: "cfb:11", homeKey: "cfb:12" });
  assert.equal(bets.matchBets(decoy, [bet], { lines: [decoy] }).matches.length, 1); // the name rule alone would take the decoy
  assert.equal(bets.matchBets(decoy, [bet], { lines }).matches.length, 0);
  assert.equal(bets.matchBets(rows["Spread 1 89"], [bet], { lines }).matches[0].tier, "same_line");
  const flags = bets.annotateRows(lines, [bet]);
  assert.equal(flags[lines.length - 1].tier, null);
  assert.deepEqual(bets.unmatchedReasons([bet], lines), []);
  // The Edges list may show only the decoy's row: the whole board still decides the game.
  assert.equal(bets.annotateRows([decoy], [bet])[0].tier, "same_side"); // names alone would size the decoy
  assert.equal(bets.annotateRows([decoy], [bet], { lines })[0].tier, null);
});

test("id join: malformed, missing or other-league ids are no id — never an error, never a join", () => {
  const rows = venueRows();
  const lines = Object.values(rows);
  const malformedKalshi = [undefined, null, "", "KXNCAAFSPREAD_26SEP19DUQWSU", "-26SEP19DUQWSU", "KXNCAAFSPREAD-26sep19duqwsu", 42]
    .map((eventTicker) => kalshiDuqWsu("KXNCAAFSPREAD-26SEP19DUQWSU-WSU36", "yes", 35.5, { venueIds: { marketTicker: "KXNCAAFSPREAD-26SEP19DUQWSU-WSU36", eventTicker } }))
    .map((record, index) => Object.assign(record, { id: `kalshi:bad${index}` }));
  const malformedNovig = ["01A0A05E-E864-7B73-9B80-9E502E4E98D5", "oc-sp-chi", 7, null]
    .map((outcomeId, index) => novigDuqWsu({ id: `novig:bad${index}`, venueIds: { outcomeId } }));
  const noMap = [novigDuqWsu({ id: "novig:nomap", venueIds: null }), novigDuqWsu({ id: "novig:string", venueIds: "01a0a05e" })];
  const otherLeague = novigDuqWsu({ id: "novig:nfl", league: "nfl" });
  const all = malformedKalshi.concat(malformedNovig, noMap, [otherLeague]);
  assert.equal(bets.matchBets(rows["Spread 1 89"], all, { lines }).matches.length, 0);
  const reasons = bets.unmatchedReasons(all, lines).map((u) => u.reason);
  assert.equal(reasons.length, all.length);
  // No usable id: the plain name reason, no "by id" tier.
  assert.equal(reasons[0], "team not recognised (Duquesne Dukes (venue spelling), Wazzu (venue spelling))");
  assert.ok(reasons.slice(0, all.length - 1).every((reason) => !reason.startsWith("by id")));
  // An id on an event of another league is no join: the NFL bet is not placed on the CFB game.
  assert.ok(reasons[all.length - 1].startsWith("by id: Novig outcome only on another league's board ladder; by name:"));
  // A row with a malformed map is skipped, not thrown on.
  const brokenRow = Object.assign({}, rows["Spread 1 89"], { venueIds: { kalshiEventSuffixes: "26SEP19DUQWSU", novigOutcomes: 3 } });
  assert.equal(bets.matchBets(brokenRow, [novigDuqWsu()], { lines: [brokenRow] }).matches.length, 0);
});

test("id join: unmatched reasons name the id tier and the name tier that both missed", () => {
  const rows = venueRows();
  const lines = Object.values(rows);
  const kalshiElsewhere = kalshiDuqWsu("KXNCAAFSPREAD-26SEP26ELONURI-URI7", "yes", 6.5, { subTitle: "ELON vs URI (Sep 26)" });
  const novigOffLadder = novigDuqWsu({ id: "novig:off", venueIds: { outcomeId: "01a0a05e-e864-7b73-9b80-000000000000" }, eventDate: "2026-09-26", eventStart: null });
  const novigMoneyline = novigDuqWsu({ id: "novig:ml", betType: "moneyline", points: null, eventDate: "2026-09-26", eventStart: null });
  const reasons = bets.unmatchedReasons([kalshiElsewhere, novigOffLadder, novigMoneyline], lines).map((u) => u.reason);
  assert.deepEqual(reasons, [
    "by id: Kalshi event 26SEP26ELONURI not on any board ladder; by name: team not recognised (Duquesne Dukes (venue spelling), Wazzu (venue spelling))",
    "by id: Novig outcome not on any board ladder; by name: team not recognised (Duquesne Dukes (venue spelling), Wazzu (venue spelling))",
    // Moneylines have no rungs, so a Novig moneyline has no id tier to report.
    "team not recognised (Duquesne Dukes (venue spelling), Wazzu (venue spelling))",
  ]);
});

test("id join: BetOnline records have no venue ids and stay on the name / rotation path on an id-carrying board", () => {
  const rows = venueRows();
  const lines = Object.values(rows);
  const bet = betonlineRecord({
    league: "cfb", awayTeam: "Duquesne", homeTeam: "Washington State", awayKey: "cfb:927", homeKey: "cfb:717", rotation: null,
    side: "home", points: -35.5, placedAt: "2026-09-15T15:00:00Z", approx: ["game_date_unknown"],
    // A stray Kalshi-shaped map is ignored: only Kalshi and Novig records join by id.
    venueIds: { eventTicker: "KXNCAAFSPREAD-26SEP19DUQWSU", outcomeId: NOVIG_WSU_MINUS_35_5 },
  });
  assert.equal(bets.matchBets(rows["Spread 1 89"], [bet], { lines }).matches[0].tier, "same_line");
  const lost = Object.assign({}, bet, { awayKey: null, awayTeam: "Duq (BOL spelling)" });
  assert.deepEqual(bets.unmatchedReasons([lost], lines).map((u) => u.reason), ["team not recognised (Duq (BOL spelling))"]);
});

test("id join: a bet with ONE team keyed still gets its side — its own key, or the other team's on the opposite side", () => {
  const rows = venueRows();
  const lines = Object.values(rows);
  // Duquesne moneyline with only the OTHER team (WSU) keyed: before, every row was same_game and held $0.
  const duqMoneyline = Object.assign(kalshiDuqWsu("KXNCAAFGAME-26SEP19DUQWSU-DUQ", "yes", null, { subTitle: "DUQ vs WSU (Sep 19)" }), { homeKey: "cfb:717" });
  assert.equal(duqMoneyline.side, "away");
  assert.equal(bets.matchBets(rows["Moneyline 0 105"], [duqMoneyline], { lines }).matches[0].tier, "same_line");
  assert.equal(bets.matchBets(rows["Moneyline 1 105"], [duqMoneyline], { lines }).matches[0].tier, "opposite");
  // NO on WSU36 = Duquesne +35.5 with only Duquesne keyed: the WSU row carries the "against" dollars.
  const duqPlus = Object.assign(kalshiDuqWsu("KXNCAAFSPREAD-26SEP19DUQWSU-WSU36", "no", 35.5), { awayKey: "cfb:927" });
  assert.equal(bets.matchBets(rows["Spread 0 89"], [duqPlus], { lines }).matches[0].tier, "same_line");
  assert.equal(bets.exposureOf(bets.matchBets(rows["Spread 1 89"], [duqPlus], { lines }).matches).against, duqPlus.stake);
  // A key for a team not in the game is ignored; the keyed WSU side still places the bet.
  const wrongOther = Object.assign(kalshiDuqWsu("KXNCAAFSPREAD-26SEP19DUQWSU-WSU36", "yes", 35.5), { awayKey: "cfb:1", homeKey: "cfb:717" });
  assert.equal(bets.matchBets(rows["Spread 0 89"], [wrongOther], { lines }).matches[0].tier, "opposite");
});

test("id join: a row without an eventId falls back to the name rule", () => {
  const rows = venueRows();
  const lines = Object.values(rows);
  const keyed = Object.assign(kalshiDuqWsu("KXNCAAFSPREAD-26SEP19DUQWSU-WSU36", "yes", 35.5), { awayKey: "cfb:927", homeKey: "cfb:717" });
  const ticketRow = Object.assign({}, rows["Spread 1 89"], { eventId: null, venueIds: null });
  assert.equal(bets.matchBets(ticketRow, [keyed], { lines }).matches[0].tier, "same_line");
});

// ---- team crosswalk (#118 step 4) -------------------------------------------
//
// Duquesne is Unabated id 927 and Washington State 717 on the venue-ids
// slice; the records above spell them so no name resolves.

const DUQ_KEY = "cfb:927";
const WSU_KEY = "cfb:717";

function crosswalkKeys(rows) {
  return rows.map((row) => [row.venue, row.venueTeamKey, row.unabatedTeamId, row.unabatedTeamName]);
}

test("crosswalk: a Kalshi id join teaches both event-title names; the rows carry where they came from", () => {
  const lines = Object.values(venueRows());
  const bet = kalshiDuqWsu("KXNCAAFSPREAD-26SEP19DUQWSU-WSU36", "yes", 35.5);
  const { learned, conflicts } = bets.learnCrosswalk([bet], lines, []);
  assert.deepEqual(conflicts, []);
  assert.deepEqual(crosswalkKeys(learned), [
    ["kalshi", "Duquesne Dukes (venue spelling)", "927", "Duquesne"],
    ["kalshi", "Wazzu (venue spelling)", "717", "Washington State"],
  ]);
  assert.deepEqual(learned.map((row) => [row.league, row.venueTeamName, row.learnedFrom]), [
    ["cfb", "Duquesne Dukes (venue spelling)", `${bet.id} on board event 123742`],
    ["cfb", "Wazzu (venue spelling)", `${bet.id} on board event 123742`],
  ]);
  // Two bets on the game teach the same two rows once; rows already held teach nothing.
  const twice = bets.learnCrosswalk([bet, kalshiDuqWsu("KXNCAAFGAME-26SEP19DUQWSU-WSU", "yes", null)], lines, []);
  assert.equal(twice.learned.length, 2);
  assert.deepEqual(bets.learnCrosswalk([bet], lines, learned).learned, []);
});

test("crosswalk: a Novig id join keys on Novig's team ids, naming the team as Novig does", () => {
  const lines = Object.values(venueRows());
  const bet = novigDuqWsu({
    awayTeamVenue: { id: "nv-duq", name: "Duquesne Dukes (venue spelling)", shortName: "Duquesne", symbol: "DUQ" },
    homeTeamVenue: { id: "nv-wsu", name: "Wazzu (venue spelling)", shortName: "Wazzu", symbol: "WSU" },
  });
  const { learned, conflicts } = bets.learnCrosswalk([bet], lines, []);
  assert.deepEqual(conflicts, []);
  assert.deepEqual(crosswalkKeys(learned), [["novig", "nv-duq", "927", "Duquesne"], ["novig", "nv-wsu", "717", "Washington State"]]);
  assert.deepEqual(learned.map((row) => row.venueTeamName), ["Duquesne Dukes (venue spelling)", "Wazzu (venue spelling)"]);
  // A Novig record without team ids (an older service) keys on the name, like Kalshi; without names on one side it teaches nothing.
  const noIds = novigDuqWsu({ awayTeamVenue: { id: null, name: "Duquesne Dukes (venue spelling)" }, homeTeamVenue: null });
  assert.deepEqual(crosswalkKeys(bets.learnCrosswalk([noIds], lines, []).learned), [
    ["novig", "Duquesne Dukes (venue spelling)", "927", "Duquesne"], ["novig", "Wazzu (venue spelling)", "717", "Washington State"]]);
  assert.deepEqual(bets.learnCrosswalk([novigDuqWsu({ homeTeam: null })], lines, []).learned, []);
});

test("crosswalk: nothing is learned from a name match, an ambiguous id, a closed bet, or a row without both Unabated ids", () => {
  const rows = venueRows();
  const lines = Object.values(rows);
  const bet = kalshiDuqWsu("KXNCAAFSPREAD-26SEP19DUQWSU-WSU36", "yes", 35.5);
  // The names resolve and the id is not on the board: a name match teaches nothing.
  const byName = Object.assign({}, bet, { awayKey: DUQ_KEY, homeKey: WSU_KEY, venueIds: { eventTicker: "KXNCAAFSPREAD-26SEP26ELONURI" } });
  assert.equal(bets.matchBets(rows["Spread 1 89"], [byName], { lines }).matches.length, 1);
  assert.deepEqual(bets.learnCrosswalk([byName], lines, []).learned, []);
  // The suffix on two board events.
  const twin = Object.assign({}, rows["Spread 1 105"], { eventId: 999999 });
  assert.deepEqual(bets.learnCrosswalk([bet], lines.concat([twin]), []), { learned: [], conflicts: [] });
  // A settled bet, and a board row missing a team id.
  assert.deepEqual(bets.learnCrosswalk([Object.assign({}, bet, { status: "won" })], lines, []).learned, []);
  const idless = lines.map((line) => Object.assign({}, line, { homeTeamId: null }));
  assert.deepEqual(bets.learnCrosswalk([bet], idless, []).learned, []);
});

test("crosswalk: a venue name that already resolves to a DIFFERENT team is a conflict — neither side is learned", () => {
  const lines = Object.values(venueRows());
  // The venue's away team reads as the Bears (nfl slice key) — a spelling teams.js knows for another team.
  const swapped = Object.assign(kalshiDuqWsu("KXNCAAFSPREAD-26SEP19DUQWSU-WSU36", "yes", 35.5), { awayTeam: "Washington St.", homeTeam: "Duquesne" });
  assert.equal(key("cfb", "Washington St."), WSU_KEY);
  const { learned, conflicts } = bets.learnCrosswalk([swapped], lines, []);
  assert.deepEqual(learned, []);
  assert.equal(conflicts.length, 1);
  assert.equal(conflicts[0].side, "away");
  assert.equal(conflicts[0].venueTeamKey, "Washington St.");
  assert.match(conflicts[0].reason, /resolves to cfb:717, the joined event says cfb:927/);
  // A held row with another id for the same venue team is a conflict too, never rewritten.
  const held = [{ venue: "kalshi", league: "cfb", venueTeamKey: "Duquesne Dukes (venue spelling)", unabatedTeamId: "1" }];
  const result = bets.learnCrosswalk([kalshiDuqWsu("KXNCAAFSPREAD-26SEP19DUQWSU-WSU36", "yes", 35.5)], lines, held);
  assert.deepEqual(result.learned, []);
  assert.match(result.conflicts[0].reason, /crosswalk holds cfb:1, the joined event says cfb:927/);
  // A name that resolves to the SAME team the board names is fine (it confirms the row).
  const agreeing = Object.assign(kalshiDuqWsu("KXNCAAFSPREAD-26SEP19DUQWSU-WSU36", "yes", 35.5), { homeTeam: "Washington St." });
  assert.deepEqual(crosswalkKeys(bets.learnCrosswalk([agreeing], lines, []).learned), [
    ["kalshi", "Duquesne Dukes (venue spelling)", "927", "Duquesne"], ["kalshi", "Washington St.", "717", "Washington State"]]);
});

test("crosswalk: a later bet on a game with NO ladder matches through the crosswalk alone, and a cleared table takes the keys back", () => {
  const rows = venueRows();
  const lines = Object.values(rows);
  const taught = bets.learnCrosswalk([kalshiDuqWsu("KXNCAAFSPREAD-26SEP19DUQWSU-WSU36", "yes", 35.5)], lines, []).learned;
  // The venue lists no ladder for the game: no venue ids on the board, the id has nothing to join.
  const bare = lines.map((line) => Object.assign({}, line, { venueIds: null }));
  const later = kalshiDuqWsu("KXNCAAFGAME-26SEP19DUQWSU-WSU", "yes", null, { venueIds: { marketTicker: "KXNCAAFGAME-26SEP19DUQWSU-WSU", eventTicker: "KXNCAAFGAME-26SEP19DUQWSU" } });
  assert.deepEqual(bets.unmatchedReasons([later], bare).map((u) => u.reason),
    ["by id: Kalshi event 26SEP19DUQWSU not on any board ladder; by name: team not recognised (Duquesne Dukes (venue spelling), Wazzu (venue spelling))"]);
  const [keyed] = bets.resolveTeamKeys([later], taught);
  assert.deepEqual([keyed.awayKey, keyed.homeKey], [DUQ_KEY, WSU_KEY]);
  assert.equal(later.awayKey, null); // input untouched
  const bareRows = venueRows();
  const bareMl = Object.assign({}, bareRows["Moneyline 1 105"], { venueIds: null });
  assert.equal(bets.matchBets(bareMl, [keyed], { lines: bare }).matches[0].tier, "same_line");
  assert.deepEqual(bets.unmatchedReasons([keyed], bare), []);
  // The crosswalk is applied whatever the name resolves to, and only for its own venue.
  const novigSame = Object.assign(novigDuqWsu(), { venueIds: null });
  assert.deepEqual(bets.resolveTeamKeys([novigSame], taught).map((r) => [r.awayKey, r.homeKey]), [[null, null]]);
  // Cleared: rekeyRecords drops the learned keys again (resolveTeamKeys alone keeps a key it finds).
  assert.deepEqual(bets.resolveTeamKeys([keyed], []).map((r) => [r.awayKey, r.homeKey]), [[DUQ_KEY, WSU_KEY]]);
  assert.deepEqual(bets.rekeyRecords([keyed], []).map((r) => [r.awayKey, r.homeKey]), [[null, null]]);
  // Malformed served rows are skipped, never thrown on.
  const junk = [null, 7, { venue: "kalshi" }, { venue: "kalshi", league: "cfb", venueTeamKey: "Wazzu (venue spelling)", unabatedTeamId: "" }];
  assert.deepEqual(bets.resolveTeamKeys([later], junk).map((r) => [r.awayKey, r.homeKey]), [[null, null]]);
});

test("crosswalk: venueTeamOf reads Novig's team id, then its name, then the record's own name", () => {
  const novig = novigDuqWsu({ awayTeamVenue: { id: "nv-duq", name: "Duquesne", shortName: "DUQ" }, homeTeamVenue: { id: "", name: null, shortName: "Wazzu" } });
  assert.deepEqual(bets.venueTeamOf(novig, "away"), { key: "nv-duq", name: "Duquesne" });
  assert.deepEqual(bets.venueTeamOf(novig, "home"), { key: "Wazzu", name: "Wazzu" });
  assert.deepEqual(bets.venueTeamOf(kalshiDuqWsu("KXNCAAFGAME-26SEP19DUQWSU-WSU", "yes", null), "home"), { key: "Wazzu (venue spelling)", name: "Wazzu (venue spelling)" });
  assert.equal(bets.venueTeamOf({ awayTeam: null, awayTeamVenue: null }, "away"), null);
  assert.equal(bets.venueTeamOf({ awayTeam: "", awayTeamVenue: { id: null, name: "" } }, "away"), null);
});

test("BetOnline: a reused rotation whose named teams are not in next week's game is not ambiguous on the Edges board", () => {
  const rows = nflRows();
  const thisWeek = rows["Total FG 1 46"];
  const nextWeek = Object.assign({}, thisWeek, { eventId: 999001, awayTeamId: 999, homeTeamId: 998,
    eventStartMs: thisWeek.eventStartMs + 7 * 24 * 3600 * 1000 });
  const bet = bets.resolveTeamKeys([betonlineRecord({
    id: "betonline:reused-rotation", homeTeam: "Carolina Panthers", betType: "total", side: "under", points: 46, rotation: 466,
  })])[0];
  assert.ok(bet.awayKey && bet.homeKey);
  const flags = bets.annotateRows([thisWeek], [bet], { lines: [thisWeek, nextWeek] });
  assert.equal(flags[0].tier, "same_line");
  assert.equal(flags[0].exposure.held, bet.stake);
  assert.equal(bets.matchBets(nextWeek, [bet], { lines: [thisWeek, nextWeek] }).matches.length, 0);
});
