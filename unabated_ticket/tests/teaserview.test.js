// Run: node --test unabated_ticket/tests/teaserview.test.js
//
// The Teasers tab's words (teaserview.js), shared by the panel and the phone,
// on hand-built teaser.js outputs: the same inputs must read the same way on
// both pages.
const test = require("node:test");
const assert = require("node:assert/strict");
const teaser = require("../extension/teaser.js");
const teaserView = require("../extension/teaserview.js");

const NFL = 1;
const CFB = 2;
const BET_TYPE_SPREAD = 2;
const BET_TYPE_TOTAL = 3;

function legOf(overrides) {
  return {
    key: "e1:bt2:0", eventId: "e1", leagueId: NFL, betTypeId: BET_TYPE_SPREAD, leagueLabel: "NFL",
    label: "Bills +7.5", fromLabel: "from +1.5", bookLabel: "+1.5 -110", matchup: "Bills @ Jets",
    eventStartMs: Date.parse("2026-10-11T17:00:00Z"), modifiedMs: Date.parse("2026-10-11T12:00:00Z"), win: 0.75,
    ...overrides,
  };
}

function rowOf(standing, legOverrides, rowOverrides) {
  return { leg: legOf(legOverrides), standing, inTickets: 0, openTickets: 0, note: null, straights: null, ...rowOverrides };
}

test("the empty list says why: no board yet, too few legs, a partial held back, nothing worth betting", () => {
  assert.equal(teaserView.emptyText(null), "Waiting for Buckeye's NFL and CFB board…");
  assert.equal(teaserView.emptyText({ ticketCount: 2, reason: null, partialHeldBack: 0, placedCount: 0 }), null);
  assert.match(teaserView.emptyText({ ticketCount: 0, reason: teaser.REASON_FEW_LEGS, partialHeldBack: 0, placedCount: 0 }), /^Fewer than 4 games/);
  assert.match(teaserView.emptyText({ ticketCount: 0, reason: null, partialHeldBack: 120, placedCount: 1 }), /partial \$120/);
  assert.match(teaserView.emptyText({ ticketCount: 0, reason: null, partialHeldBack: 0, placedCount: 2 }), /with the open teasers held/);
  assert.match(teaserView.emptyText({ ticketCount: 0, reason: null, partialHeldBack: 0, placedCount: 0 }), /^No ticket worth betting/);
});

test("planStateOf reads the list's state off a plan and its build", () => {
  const state = teaserView.planStateOf({ tickets: [{}, {}] }, { reason: null, partialHeldBack: 0, placed: [{}] });
  assert.deepEqual(state, { ticketCount: 2, reason: null, partialHeldBack: 0, placedCount: 1 });
});

test("summary: a fresh list shows Expected with its return, an open set shows Placed and Expected over all", () => {
  const straights = { count: 0, stake: 0, games: 0 };
  const fresh = teaserView.summaryView({ count: 3, placedCount: 0, stake: 600, expected: 60, makesMoney: 0.41, allLose: 0.2, straights }, 0, { bankroll: 5000, multiplier: 0.25 });
  assert.equal(fresh.label, "Bet 3 tickets");
  assert.equal(fresh.stake, "$600");
  assert.deepEqual(fresh.cells[0], { label: "Expected", value: "+$60", small: "10%" });
  assert.deepEqual(fresh.cells.map((cell) => cell.label), ["Expected", "Makes money", "All lose"]);
  assert.equal(fresh.note, "Quarter Kelly on $5,000 · tickets that share a leg are sized together");

  const open = teaserView.summaryView({
    count: 1, placedCount: 2, stake: 200, placedStake: 400, expectedAll: -12, makesMoney: null, allLose: null,
    straights: { count: 2, stake: 150, games: 1 },
  }, 3, { bankroll: 5000, multiplier: 1 });
  assert.equal(open.label, "Bet 1 more");
  assert.deepEqual(open.cells, [
    { label: "Placed, in play", value: "$400", small: null },
    { label: "Expected, all 3", value: "-$12", small: null },
  ]);
  assert.equal(open.note, "Full Kelly on $5,000 · open BFA teasers held fixed · counts $150 of straight bets on 1 game");
  assert.equal(teaserView.summaryView({ count: 0, placedCount: 0, straights }, 0, { bankroll: 1, multiplier: 1 }), null);
});

test("ticket cards: a ticket to bet and an open one read in the panel's words", () => {
  const pending = teaserView.pendingTicketView({ number: 2, stake: 200, ev: 0.0512, legs: [legOf({})] });
  assert.deepEqual(pending, {
    number: 2, size: "4-team · pays +300", stake: "$200", ev: "EV +5.1%",
    legs: [{ label: "Bills +7.5", from: "from +1.5", win: "75.0%" }],
  });
  assert.equal(teaserView.moreTicketsText([{ stake: 200 }, { stake: 150 }]), "2 more tickets · $350");

  const done = teaserView.openTicketView({
    inPlay: false, reason: null, legCount: 4, stake: 100, toWin: 300, placedAt: null,
    legs: [{ label: "A", note: "started", win: null, state: teaser.LEG_STARTED }],
  });
  assert.equal(done.size, "4-team · $100 to win $300");
  assert.equal(done.ev, null);
  assert.equal(done.outOfMath, "every game has started: out of the math");
  assert.deepEqual(done.legs, [{ label: "A", note: "started", win: "—", counted: true }]);
  assert.equal(done.placedTag, "placed");
});

test("Open at BFA header: count, dollars in play and how fresh BFA's pull is", () => {
  assert.deepEqual(teaserView.openBlockView([], null), { count: "", note: "none open · bets service not reached yet" });
  const placed = [{ inPlay: true, stake: 200 }, { inPlay: false, stake: 100 }];
  assert.deepEqual(teaserView.openBlockView(placed, { configured: true, fetchedAt: 1, ageText: "41 s" }), {
    count: "2", note: "$300, $200 in play · each until its last game starts · BFA pulled 41 s ago",
  });
});

test("Can't tease shows on college legs on show only; a marked market offers Restore", () => {
  assert.equal(teaserView.blockControl(rowOf(teaser.STANDING_POOL, { leagueId: NFL })), null);
  assert.equal(teaserView.blockControl(rowOf(teaser.STANDING_BELOW, { leagueId: CFB })), null);
  const block = teaserView.blockControl(rowOf(teaser.STANDING_OUT, { leagueId: CFB, betTypeId: BET_TYPE_TOTAL }));
  assert.equal(block.action, "block");
  assert.match(block.title, /game's total \(both sides\)/);
  assert.equal(teaserView.blockControl(rowOf(teaser.STANDING_BLOCKED, { leagueId: CFB })).action, "restore");
});

test("legs layout: break-even divider where the shown legs cross it, the rest folded with counts", () => {
  const rows = [
    rowOf(teaser.STANDING_POOL, { eventId: "a", win: 0.76 }, { inTickets: 2 }),
    rowOf(teaser.STANDING_OUT, { eventId: "b", win: 0.70 }, { note: "outside the top 10" }),
    rowOf(teaser.STANDING_BELOW, { eventId: "c", win: 0.6 }, { note: "below break-even" }),
    rowOf(teaser.STANDING_BLOCKED, { eventId: "d", leagueId: CFB }, { note: "can't tease" }),
  ];
  const layout = teaserView.legsLayout(rows);
  assert.equal(layout.gameCount, 4);
  assert.deepEqual(layout.items.map((item) => (item.divider ? "divider" : item.row.leg.eventId)), ["a", "divider", "b"]);
  assert.deepEqual(layout.folded.map((row) => row.leg.eventId), ["d", "c"]);
  assert.equal(layout.foldLabel, "1 more game: 1 below break-even · 1 can't tease");
});

test("a leg row: its words, its straight-bet tags and its standing", () => {
  const straights = { held: 50, against: 0, counted: [{ label: "Bills +3" }], leftOut: [] };
  const rowView = teaserView.legRowView(rowOf(teaser.STANDING_POOL, {}, { inTickets: 1, openTickets: 2, straights }));
  assert.equal(rowView.pool, true);
  assert.equal(rowView.book, "Buckeye +1.5 -110");
  assert.equal(rowView.teased, "teased 6");
  assert.equal(rowView.standing, "in 1 ticket · 2 open");
  assert.deepEqual(rowView.tags, [{ kind: "held", text: "held $50", title: "Bills +3" }]);
  assert.equal(rowView.control, null);
});
