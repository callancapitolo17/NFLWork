// Run: node --test unabated_ticket/tests/teaser.test.js
//
// Buckeye teasers (docs/2026-09-27-unabated-ticket-teasers-plan.md) on
// synthetic boards: every game is built here, its Unabated fair ladder from
// rungs whose `bacr` is chosen to land on round win chances (-400 = 80%).
const test = require("node:test");
const assert = require("node:assert/strict");
const feed = require("../extension/feed.js");
const teaser = require("../extension/teaser.js");

const NFL = 1;
const FULL_GAME = 1;
const MONEYLINE = 1;
const SPREAD = 2;
const TOTAL = 3;
const AWAY_OR_OVER = 0;
const HOME_OR_UNDER = 1;
const BUCKEYE = 59;
// Any other book: its rungs are the fair ladder.
const FAIR_BOOK = 4;
const HOUR = 3600 * 1000;
const NOW = Date.parse("2026-09-27T16:00:00Z");
const KICKOFF = Date.parse("2026-09-27T17:00:00Z");
const LINE_AGE = { maxLineAgeMs: 168 * HOUR };

const nearly = (actual, expected, tol = 1e-9) =>
  assert.ok(Math.abs(actual - expected) <= tol, `expected ${expected}, got ${actual}`);

function probOf(american) {
  return american > 0 ? 100 / (american + 100) : -american / (-american + 100);
}

// One game: Buckeye's spread (the away side's points) and total, and the fair
// rungs [cut, bacr of the side that wins above the cut]. A margin rung at cut
// c is the away spread at -c; a total rung is the Over at c.
function addGame(state, game) {
  const { eventId, away, home, awayRotation, homeRotation } = game;
  const awayTeamId = eventId * 10 + 1;
  const homeTeamId = eventId * 10 + 2;
  state.teams[awayTeamId] = away;
  state.teams[homeTeamId] = home;
  state.events[eventId] = {
    eventId, leagueId: NFL, eventName: `${away} @ ${home}`, eventStart: game.startMs ?? KICKOFF,
    awayTeamId, homeTeamId, awayRotation, homeRotation, venueIds: { kalshiEventSuffixes: [], kalshiContracts: {}, novigOutcomes: {} },
  };
  const put = (fields) => {
    const line = {
      isAlt: false, leagueId: NFL, periodTypeId: FULL_GAME, eventId, marketId: eventId * 100 + fields.betTypeId,
      price: -110, sourceFormat: 1, sourcePrice: null, bacr: null, ge: null, liquidity: null, statusId: 1,
      sequenceNumber: null, isBlurred: false, modifiedOn: new Date(NOW - HOUR).toISOString().slice(0, 19), ...fields,
    };
    line.sideKey = `si${line.sideIndex}`;
    line.key = `${line.bookId}:${eventId}:${line.betTypeId}:${line.sideIndex}:${line.points}`;
    state.lines[line.key] = line;
  };
  if (game.spread != null) {
    put({ bookId: BUCKEYE, betTypeId: SPREAD, sideIndex: AWAY_OR_OVER, points: game.spread });
    put({ bookId: BUCKEYE, betTypeId: SPREAD, sideIndex: HOME_OR_UNDER, points: -game.spread || 0 });
  }
  if (game.total != null) {
    put({ bookId: BUCKEYE, betTypeId: TOTAL, sideIndex: AWAY_OR_OVER, points: game.total });
    put({ bookId: BUCKEYE, betTypeId: TOTAL, sideIndex: HOME_OR_UNDER, points: game.total });
  }
  for (const [cut, bacr] of game.marginRungs || []) put({ bookId: FAIR_BOOK, isAlt: true, betTypeId: SPREAD, sideIndex: AWAY_OR_OVER, points: -cut || 0, bacr });
  for (const [cut, bacr] of game.totalRungs || []) put({ bookId: FAIR_BOOK, isAlt: true, betTypeId: TOTAL, sideIndex: AWAY_OR_OVER, points: cut, bacr });
  for (const [sideIndex, bacr] of game.moneyline || []) put({ bookId: FAIR_BOOK, betTypeId: MONEYLINE, sideIndex, points: null, bacr });
  return state;
}

function emptyBoard() {
  return { leagues: [NFL], teams: {}, teamIndex: {}, books: { [BUCKEYE]: { id: BUCKEYE, name: "Buckeye", isLive: true }, [FAIR_BOOK]: { id: FAIR_BOOK, name: "BetMGM", isLive: true } }, events: {}, lines: {} };
}

function legsOf(state, now) {
  return teaser.teaserLegs(state, { now: now ?? NOW, ...LINE_AGE });
}

function legByLabel(legs, label) {
  const leg = legs.find((candidate) => candidate.label === label);
  assert.ok(leg, `leg ${label} missing; have ${legs.map((candidate) => candidate.label).join(", ")}`);
  return leg;
}

// Games whose away side, teased from +2.5 to +8.5, wins with the given chance.
// Each game's other three legs sit far below break-even, so the pool is these legs.
function awayLegBoard(wins) {
  const state = emptyBoard();
  wins.forEach((bacr, index) => {
    addGame(state, {
      eventId: 101 + index, away: `Away ${index + 1}`, home: `Home ${index + 1}`, awayRotation: 401 + 2 * index, homeRotation: 402 + 2 * index,
      spread: 2.5, marginRungs: [[-8.5, bacr], [3.5, 150]],
    });
  });
  return state;
}

// One describeLine row per event, the way the panel's boardLines() builds them.
function boardLinesOf(state) {
  const seen = new Set();
  const rows = [];
  for (const line of Object.values(state.lines)) {
    if (line.isAlt || seen.has(line.eventId)) continue;
    seen.add(line.eventId);
    rows.push(feed.describeLine(line, state));
  }
  return rows;
}

// One leg of an open BFA teaser as the bets service serves it: the teased
// number, the rotation, the ticket's dollars on every leg.
function bfaLeg(parlayId, legIndex, fields) {
  return {
    id: `bfa:${parlayId}:leg${legIndex}`, source: "bfa_api", venue: "bfa", status: "open", league: "nfl", betType: "spread",
    period: "FG", side: "away", points: 8.5, price: -110, rotation: 401, awayTeam: "AWAY", homeTeam: null, awayKey: null, homeKey: null,
    eventStart: new Date(KICKOFF).toISOString(), stake: 200, toWin: 600, isParlayLeg: true, parlayId: `bfa:${parlayId}`, legIndex,
    legCount: 4, placedAt: "2026-09-27T15:30:00Z", approx: ["side_from_rotation_parity"], unmatchable: null,
    raw: { headerDescription: "4 TEAM TEASERS" }, ...fields,
  };
}

function openTeasersOf(state, records, now) {
  return teaser.openTeasers(records, boardLinesOf(state), { now: now ?? NOW, ladderOf: teaser.teaserBoardOf(state).ladderOf });
}

function plan(state, options) {
  const legs = legsOf(state);
  const placed = options && options.placed ? options.placed : [];
  return teaser.planTeasers({ legs, placed, kellyBankroll: (options && options.kellyBankroll) || 5000, previous: (options && options.previous) || null });
}

test("break-even: four equal legs at +300 need 70.7% each", () => {
  nearly(teaser.BREAK_EVEN_WIN, 0.25 ** 0.25);
  nearly(teaser.BREAK_EVEN_WIN, 0.7071, 1e-4);
});

test("teased numbers: spreads +6 for both sides, Over -6, Under +6, each read at its own cut", () => {
  const state = addGame(emptyBoard(), {
    eventId: 1, away: "Seattle Seahawks", home: "Washington Commanders", awayRotation: 461, homeRotation: 462,
    spread: -8.5, total: 44.5,
    // away -2.5 wins above 2.5; home +14.5 wins below 14.5; Over 38.5 above 38.5; Under 50.5 below 50.5.
    marginRungs: [[2.5, -300], [14.5, 400]],
    totalRungs: [[38.5, -400], [50.5, 150]],
  });
  const legs = legsOf(state);
  const away = legByLabel(legs, "Seattle Seahawks -2.5");
  assert.deepEqual([away.axis, away.direction, away.winCut, away.fromLabel, away.bookLabel], ["margin", "above", 2.5, "from -8.5", "-8.5 -110"]);
  nearly(away.win, 0.75);
  const home = legByLabel(legs, "Washington Commanders +14.5");
  assert.deepEqual([home.direction, home.winCut], ["below", 14.5]);
  nearly(home.win, 1 - probOf(400));
  const over = legByLabel(legs, "Over 38.5");
  assert.deepEqual([over.axis, over.direction, over.winCut, over.fromLabel], ["total", "above", 38.5, "from 44.5"]);
  nearly(over.win, 0.8);
  const under = legByLabel(legs, "Under 50.5");
  assert.deepEqual([under.direction, under.winCut], ["below", 50.5]);
  nearly(under.win, 1 - probOf(150));
});

test("a push loses: a whole number reads the half-point one step against the bet", () => {
  // Bills -1 (from -7) is priced at -1.5, Chargers +13 at +12.5, Over 38 at 38.5, Under 50 at 49.5.
  // The rungs one step the OTHER way carry different fairs, so reading them would show.
  const state = addGame(emptyBoard(), {
    eventId: 2, away: "Los Angeles Chargers", home: "Buffalo Bills", awayRotation: 477, homeRotation: 478,
    spread: 7, total: 44,
    marginRungs: [[-12.5, -900], [-13.5, -1900], [-1.5, 300], [-0.5, 250]],
    totalRungs: [[37.5, -500], [38.5, -400], [49.5, 200], [50.5, 300]],
  });
  const legs = legsOf(state);
  const bills = legByLabel(legs, "Buffalo Bills -1");
  assert.deepEqual([bills.winCut, bills.pricedPoints], [-1.5, -1.5]);
  nearly(bills.win, 1 - probOf(300));
  const chargers = legByLabel(legs, "Los Angeles Chargers +13");
  assert.deepEqual([chargers.winCut, chargers.pricedPoints], [-12.5, 12.5]);
  nearly(chargers.win, probOf(-900));
  const over = legByLabel(legs, "Over 38");
  assert.deepEqual([over.winCut, over.pricedPoints], [38.5, 38.5]);
  nearly(over.win, 0.8);
  const under = legByLabel(legs, "Under 50");
  assert.deepEqual([under.winCut, under.pricedPoints], [49.5, 49.5]);
  nearly(under.win, 1 - probOf(200));
});

test("the spread ladder leaves the moneyline out", () => {
  // Home -6.5 teases to -0.5: the -0.5 cut, where the home moneyline would also sit.
  const state = addGame(emptyBoard(), {
    eventId: 3, away: "New York Jets", home: "Detroit Lions", awayRotation: 465, homeRotation: 466,
    spread: 6.5, marginRungs: [[-0.5, 300]], moneyline: [[AWAY_OR_OVER, 150], [HOME_OR_UNDER, -150]],
  });
  const lions = legByLabel(legsOf(state), "Detroit Lions -0.5");
  nearly(lions.win, 1 - probOf(300));
});

test("no rung and a flat rung are grey with the reason, never priced", () => {
  const state = addGame(emptyBoard(), {
    eventId: 4, away: "Arizona Cardinals", home: "San Francisco 49ers", awayRotation: 481, homeRotation: 482,
    spread: 7.5, total: 48,
    // 49ers -1.5 has no rung at -1.5; Over 42 sits at 42.5 whose neighbour repeats its fair.
    marginRungs: [[-13.5, -800]],
    totalRungs: [[41.5, 110], [42.5, 110]],
  });
  const legs = legsOf(state);
  const niners = legByLabel(legs, "San Francisco 49ers -1.5");
  assert.equal(niners.win, null);
  assert.equal(niners.reason, "no Unabated fair at -1.5");
  const over = legByLabel(legs, "Over 42");
  assert.equal(over.win, null);
  assert.equal(over.reason, "Unabated's fair is flat at 42.5");
});

test("legs: started games, lines off the board, stale lines and other books are not teased", () => {
  const state = addGame(emptyBoard(), { eventId: 5, away: "A", home: "B", awayRotation: 1, homeRotation: 2, spread: 2.5, marginRungs: [[-8.5, -300]] });
  addGame(state, { eventId: 6, away: "C", home: "D", awayRotation: 3, homeRotation: 4, spread: 2.5, marginRungs: [[-8.5, -300]], startMs: NOW - HOUR });
  assert.deepEqual(legsOf(state).map((leg) => leg.eventId).sort(), [5, 5]);
  for (const line of Object.values(state.lines)) if (line.eventId === 5 && line.bookId === BUCKEYE && line.sideIndex === 0) line.statusId = 2;
  assert.deepEqual(legsOf(state).map((leg) => leg.label), ["B +3.5"]);
  for (const line of Object.values(state.lines)) if (line.bookId === BUCKEYE) line.modifiedOn = new Date(NOW - 200 * HOUR).toISOString().slice(0, 19);
  assert.equal(legsOf(state).length, 0);
});

test("grading: all four legs win pays +300, anything else loses the stake", () => {
  const { build } = plan(awayLegBoard([-400, -400, -400, -400]), { kellyBankroll: 100000 });
  assert.equal(build.tickets.length, 1);
  const view = teaser.describePlan(build, legsOf(awayLegBoard([-400, -400, -400, -400])), []);
  const [ticket] = view.tickets;
  nearly(ticket.winAll, 0.8 ** 4);
  nearly(ticket.ev, 0.8 ** 4 * 4 - 1);
  nearly(view.summary.expected, ticket.stake * (0.8 ** 4 * 4 - 1), 1e-6);
  nearly(view.summary.allLose, 1 - 0.8 ** 4);
  nearly(view.summary.makesMoney, 0.8 ** 4);
});

// The greedy adds a full $200 ticket while that still raises the growth; the
// first one that does not is the partial last ticket at its own optimum.
test("one ticket under the cap is the single-bet Kelly stake, in whole dollars", () => {
  const kellyBankroll = 300;
  const win = 0.8 ** 4;
  const kellyStake = kellyBankroll * (win * 3 - (1 - win)) / 3;
  const { build } = plan(awayLegBoard([-400, -400, -400, -400]), { kellyBankroll });
  assert.deepEqual(build.tickets.map((ticket) => ticket.stake), [Math.floor(kellyStake)]);
  assert.equal(Math.floor(kellyStake), 63);
});

test("the greedy never puts more than $200 on a ticket nor takes a ticket twice", () => {
  const state = awayLegBoard([-400, -400, -400, -400, -400, -400, -400, -400, -400, -400, -400, -400]);
  const { build } = plan(state, { kellyBankroll: 200000 });
  assert.equal(build.pool.length, teaser.POOL_SIZE);
  assert.ok(build.tickets.length > 20, `expected a big set, got ${build.tickets.length}`);
  const seen = new Set();
  build.tickets.forEach((ticket, index) => {
    assert.ok(Number.isInteger(ticket.stake) && ticket.stake >= 1 && ticket.stake <= teaser.TICKET_MAX_STAKE, `stake ${ticket.stake}`);
    if (index < build.tickets.length - 1) assert.equal(ticket.stake, teaser.TICKET_MAX_STAKE);
    assert.equal(ticket.legIndexes.length, 4);
    const combo = ticket.legIndexes.join(",");
    assert.ok(!seen.has(combo), `ticket ${combo} twice`);
    seen.add(combo);
  });
});

test("the pool is each game's best leg, the top 10 by win chance", () => {
  const state = awayLegBoard([-400, -300, -400, -300, -400, -300, -400, -300, -400, -300, -250, -250]);
  const { build } = plan(state);
  assert.equal(build.pool.length, 10);
  assert.equal(new Set(build.pool.map((leg) => leg.eventId)).size, 10);
  assert.ok(!build.pool.some((leg) => leg.eventId === 111 || leg.eventId === 112), "the two 71.4% legs are left out");
});

test("fewer than four priced games: no tickets, and the reason", () => {
  const { build } = plan(awayLegBoard([-400, -400, -400]));
  assert.deepEqual(build.tickets, []);
  assert.equal(build.reason, teaser.REASON_FEW_LEGS);
});

test("open BFA teasers: legs group into tickets and join their games by rotation", () => {
  const state = awayLegBoard([-400, -400, -400, -400]);
  const records = [
    ...[0, 1, 2, 3].map((index) => bfaLeg(1, index, { rotation: 401 + 2 * index })),
    // The home side of game 2, joined by its even rotation, at its own teased number.
    ...[0, 1].map((index) => bfaLeg(2, index, { rotation: 402 + 2 * index, side: "home", points: -2.5 + 6, awayTeam: null, homeTeam: "HOME" })),
    // Not open teasers: a BFA parlay, a settled teaser, another venue's leg.
    bfaLeg(3, 0, { raw: { headerDescription: "PARLAY (2 TEAMS)" } }),
    bfaLeg(4, 0, { status: "lost" }),
    bfaLeg(5, 0, { venue: "betonline" }),
  ];
  const placed = openTeasersOf(state, records);
  assert.deepEqual(placed.map((ticket) => ticket.id), ["bfa:1", "bfa:2"]);
  const [four, two] = placed;
  assert.deepEqual(four.legs.map((leg) => [leg.state, leg.eventId, leg.label]), [
    ["live", 101, "Away 1 +8.5"], ["live", 102, "Away 2 +8.5"], ["live", 103, "Away 3 +8.5"], ["live", 104, "Away 4 +8.5"],
  ]);
  assert.deepEqual([four.stake, four.toWin, four.inPlay], [200, 600, true]);
  nearly(four.winAll, 0.8 ** 4);
  assert.deepEqual(two.legs.map((leg) => [leg.eventId, leg.label, leg.winCut, leg.direction]), [
    [101, "Home 1 +3.5", 3.5, "below"], [102, "Home 2 +3.5", 3.5, "below"],
  ]);
  nearly(two.legs[0].win, 1 - probOf(150));
});

test("an open ticket lowers the next suggestion on the legs it shares", () => {
  // Five 80% legs, K = $1,000: alone, $200 on legs 1-2-3-4 then $84 on 1-2-3-5.
  const state = awayLegBoard([-400, -400, -400, -400, -400]);
  const alone = plan(state, { kellyBankroll: 1000 });
  const stakesOf = (build) => build.tickets.map((ticket) => [ticket.legIndexes.join(""), ticket.stake]);
  assert.deepEqual(stakesOf(alone.build), [["0123", 200], ["0124", 84]]);
  // With $200 open on 1-2-3-5, legs 1-2-3 already ride: 1-2-3-4 drops to $84.
  const records = [401, 403, 405, 409].map((rotation, index) => bfaLeg(9, index, { rotation }));
  const withOpen = plan(state, { kellyBankroll: 1000, placed: openTeasersOf(state, records), previous: alone.build });
  assert.equal(withOpen.rebuilt, true);
  assert.deepEqual(stakesOf(withOpen.build), [["0123", 84]]);
});

// Four standout legs (86-83%) and six ordinary ones (76-71%): the board where
// an open ticket's own combination came back at the top of the list.
const STANDOUT_BOARD = [-614, -567, -525, -488, -317, -300, -285, -270, -257, -245];

// A build's ticket as BFA's open list serves it once placed: awayLegBoard's
// away +8.5 legs, joined back to their games by rotation.
function openTicketOf(build, ticket, parlayId) {
  return ticket.legIndexes.map((poolIndex, legIndex) => {
    const leg = build.pool[poolIndex];
    return bfaLeg(parlayId, legIndex, { rotation: 401 + 2 * (leg.eventId - 101), stake: ticket.stake, toWin: ticket.stake * 3 });
  });
}

function ticketList(build) {
  return build.tickets.map((ticket) => `${ticket.legIndexes.map((index) => build.pool[index].key).join("/")} $${ticket.stake}`);
}

test("an open ticket's own four legs are never offered again: placing the top ticket leaves the list less that ticket", () => {
  const state = awayLegBoard(STANDOUT_BOARD);
  const first = plan(state);
  const before = ticketList(first.build);
  assert.ok(before.length >= 3, `expected a list, got ${before.length}`);
  const placed = openTeasersOf(state, openTicketOf(first.build, first.build.tickets[0], 1));
  const second = plan(state, { placed, previous: first.build });
  assert.equal(second.rebuilt, true);
  assert.deepEqual(ticketList(second.build), before.slice(1));
});

test("once the whole list is open at BFA, nothing more is offered (one partial ticket per set)", () => {
  const state = awayLegBoard(STANDOUT_BOARD);
  const first = plan(state);
  assert.ok(first.build.tickets.at(-1).stake < teaser.TICKET_MAX_STAKE, "the list ends on a partial ticket");
  const records = first.build.tickets.flatMap((ticket, index) => openTicketOf(first.build, ticket, 100 + index));
  const all = plan(state, { placed: openTeasersOf(state, records), previous: first.build });
  assert.deepEqual(all.build.tickets, []);
  assert.equal(all.build.reason, null);
});

test("past the outcome budget the pool sheds its weakest legs on games no open teaser rides on, not the list", () => {
  // Ten 80% games, then eight 60% ones that two open teasers ride on: 18 games
  // would be 2^18 outcomes, so four of the ten free pool legs go.
  const state = awayLegBoard([...Array(10).fill(-400), ...Array(8).fill(-150)]);
  const records = [0, 1, 2, 3].map((index) => bfaLeg(1, index, { rotation: 421 + 2 * index }))
    .concat([0, 1, 2, 3].map((index) => bfaLeg(2, index, { rotation: 429 + 2 * index })));
  const { build } = plan(state, { placed: openTeasersOf(state, records) });
  assert.equal(build.reason, null);
  assert.equal(build.pool.length, 6);
  assert.equal(build.structure.outcomeCount, 2 ** 14);
  assert.ok(build.tickets.length > 0);
});

test("a started leg and a leg no board game matches count as won; a ticket with none still to play leaves the math", () => {
  const state = awayLegBoard([-400, -400, -400, -400]);
  state.events[103].eventStart = NOW - HOUR;
  const records = [
    bfaLeg(7, 0, { rotation: 401 }), bfaLeg(7, 1, { rotation: 403 }),
    bfaLeg(7, 2, { rotation: 405, eventStart: new Date(NOW - HOUR).toISOString() }),
    bfaLeg(7, 3, { rotation: 999 }),
    bfaLeg(8, 0, { rotation: 405, eventStart: new Date(NOW - HOUR).toISOString() }), bfaLeg(8, 1, { rotation: 998 }),
  ];
  const [seven, eight] = openTeasersOf(state, records);
  assert.deepEqual(seven.legs.map((leg) => leg.state), ["live", "live", "started", "off_board"]);
  assert.equal(seven.legs[3].note, "no board game; counted as won");
  assert.equal(seven.inPlay, true);
  nearly(seven.winAll, 0.8 * 0.8);
  assert.deepEqual(eight.legs.map((leg) => leg.state), ["started", "off_board"]);
  assert.equal(eight.inPlay, false);
});

test("a game with an open leg offers new legs on that market only", () => {
  const state = awayLegBoard([-400, -400, -400, -400]);
  // Game 101 also has a total whose Over is its best leg by far.
  addGame(state, { eventId: 101, away: "Away 1", home: "Home 1", awayRotation: 401, homeRotation: 402, total: 44.5, totalRungs: [[38.5, -900], [50.5, 900]] });
  assert.equal(plan(state).build.pool.find((leg) => leg.eventId === 101).label, "Over 38.5");
  const records = [0, 1, 2, 3].map((index) => bfaLeg(6, index, { rotation: 401 + 2 * index }));
  const { build } = plan(state, { placed: openTeasersOf(state, records) });
  assert.equal(build.pool.find((leg) => leg.eventId === 101).label, "Away 1 +8.5");
});

test("an open leg at another number shares its game's rows with the new leg", () => {
  // Buckeye moved game 104 from +2 (open ticket: +8, priced +7.5) to +1 (new leg: +7, priced +6.5).
  const state = awayLegBoard([-400, -400, -400]);
  addGame(state, { eventId: 104, away: "Away 4", home: "Home 4", awayRotation: 407, homeRotation: 408, spread: 1, marginRungs: [[-6.5, -300], [-7.5, -400], [2.5, 150]] });
  const records = [0, 1, 2].map((index) => bfaLeg(5, index, { rotation: 401 + 2 * index }))
    .concat(bfaLeg(5, 3, { rotation: 407, points: 8 }));
  const placed = openTeasersOf(state, records);
  assert.deepEqual(placed[0].legs.map((leg) => leg.winCut), [-8.5, -8.5, -8.5, -7.5]);
  const { build } = plan(state, { placed, kellyBankroll: 100000 });
  const game = build.structure.factors.find((factor) => factor.eventId === 104);
  assert.deepEqual(game.cuts, [-7.5, -6.5]);
  assert.equal(build.structure.outcomeCount, 2 * 2 * 2 * 3);
  const view = teaser.describePlan(build, legsOf(state), placed);
  // The new ticket (+7) wins only where the open one (+8) does, so every
  // ticket loses exactly when the open one does: independent legs would say less.
  const firstThree = 0.8 ** 3;
  nearly(view.summary.allLose, 1 - firstThree * 0.8);
  assert.equal(view.tickets.length, 1);
  nearly(view.tickets[0].winAll, firstThree * 0.75);
});

test("the list holds still while fairs tick, and is rebuilt on a 1-point move", () => {
  const state = awayLegBoard([-400, -400, -400, -400, -300, -300]);
  const first = plan(state, { kellyBankroll: 50000 });
  const rung = Object.values(state.lines).find((line) => line.eventId === 101 && line.bookId === FAIR_BOOK && line.points === 8.5);
  rung.bacr = -412; // 80.0% -> 80.5%: under a point
  const ticked = plan(state, { kellyBankroll: 50000, previous: first.build });
  assert.equal(ticked.rebuilt, false);
  assert.deepEqual(ticked.build.tickets, first.build.tickets);
  const liveEv = teaser.describePlan(ticked.build, legsOf(state), []).tickets.find((ticket) => ticket.legs.some((leg) => leg.eventId === 101));
  nearly(liveEv.legs.find((leg) => leg.eventId === 101).win, probOf(-412));
  rung.bacr = -429; // 81.1%: a point or more since the build
  assert.equal(plan(state, { kellyBankroll: 50000, previous: ticked.build }).rebuilt, true);
});

test("the list is rebuilt when K changes or a pool leg's number moves", () => {
  const state = awayLegBoard([-400, -400, -400, -400]);
  const first = plan(state, { kellyBankroll: 300 });
  assert.equal(plan(state, { kellyBankroll: 300, previous: first.build }).rebuilt, false);
  assert.equal(plan(state, { kellyBankroll: 400, previous: first.build }).rebuilt, true);
  // +2.5 -> +3 changes nothing a push-loses leg wins on (+9 wins when +8.5 does)...
  const buckeyeAway = Object.values(state.lines).find((line) => line.eventId === 101 && line.bookId === BUCKEYE && line.sideIndex === 0);
  buckeyeAway.points = 3;
  assert.equal(plan(state, { kellyBankroll: 300, previous: first.build }).rebuilt, false);
  // ...+3.5 does: the leg is now +9.5.
  buckeyeAway.points = 3.5;
  addGame(state, { eventId: 101, away: "Away 1", home: "Home 1", awayRotation: 401, homeRotation: 402, marginRungs: [[-9.5, -500]] });
  const moved = plan(state, { kellyBankroll: 300, previous: first.build });
  assert.equal(moved.rebuilt, true);
  assert.ok(moved.build.pool.some((leg) => leg.label === "Away 1 +9.5"));
});

test("describeLegs: each game's pool leg with its ticket count, then the rest, the unpriced last", () => {
  const state = awayLegBoard([...Array(10).fill(-400), -300, -233]);
  addGame(state, { eventId: 120, away: "No Fair", home: "Anywhere", awayRotation: 499, homeRotation: 500, spread: 2.5 });
  const { build } = plan(state, { kellyBankroll: 300 });
  const rows = teaser.describeLegs(legsOf(state), build, []);
  assert.deepEqual(rows.map((row) => row.standing), [...Array(10).fill("pool"), "out", "below_break_even", "unpriced"]);
  assert.deepEqual(rows.slice(10).map((row) => [row.leg.label, row.note]), [
    ["Away 11 +8.5", "not in the top 10"], ["Away 12 +8.5", "below break-even"], ["No Fair +8.5", "no Unabated fair at +8.5"],
  ]);
  assert.deepEqual(rows.slice(0, 5).map((row) => row.inTickets), [1, 1, 1, 1, 0]);
});
