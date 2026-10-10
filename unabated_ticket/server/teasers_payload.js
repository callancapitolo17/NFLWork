// Unabated Ticket server — the body of GET /teasers.json: the panel's
// Teasers tab (Buckeye 6-point 4-team teasers) as words, for the phone page.
// Pure: no I/O, no clock (the runner passes `now`), so the node tests drive
// it on fixtures.
//
// The runner builds the model exactly as panel.js teaserModel does (teaser.js
// legs, open BFA teasers, held straights, planTeasers on the previous build so
// the list holds still, describePlan / describeLegs) with the Can't tease
// marks the bets service holds; every word here is teaserview.js's, the
// module the panel renders from, so both pages read the same.
//
// Shape (ms epochs):
//   {generatedAt, error (the tab failed: its message, else null), ready
//    (both football boards in or failed), status, empty (why the list is
//    empty, or null), ticketCount, summary (teaserview.summaryView, or null),
//    open: [ticket] (open BFA teasers; [] before the board is ready),
//    tickets: [{...pendingTicketView, signature, place: {body} | {error}}],
//    moreLabel (the fold under the first teaserview.TICKETS_SHOWN, "" when none),
//    legs: {gameCount, items: [{divider: true} | {row}], folded: [row],
//    foldLabel}, breakEven}
//   open ticket  teaserview.openTicketView plus {stake, inPlay}, so the page
//                can call openBlockView with its own clock and BFA row
//   row          teaserview.legRowView plus {marketKey} (the Can't tease key)
// The bets service's state and BFA's freshness are not here: the page reads
// /bets.json itself and words them with betsview.js, as the panel does.

"use strict";

const teaser = require("../extension/teaser.js");
const teaserView = require("../extension/teaserview.js");

function openTicketOf(ticket) {
  return { ...teaserView.openTicketView(ticket), stake: typeof ticket.stake === "number" ? ticket.stake : null, inPlay: Boolean(ticket.inPlay) };
}

// A ticket to bet, with what its Place button would send (teaser.placeRequestOf).
function pendingTicketOf(ticket) {
  const request = teaser.placeRequestOf(ticket);
  return {
    ...teaserView.pendingTicketView(ticket),
    signature: teaser.ticketSignatureOf(ticket),
    place: request.error ? { error: request.error } : { body: request.body },
  };
}

function legRowOf(row) {
  return { ...teaserView.legRowView(row), marketKey: teaser.marketKeyOf(row.leg) };
}

function legsOf(legRows) {
  const layout = teaserView.legsLayout(legRows);
  return {
    gameCount: layout.gameCount,
    items: layout.items.map((item) => (item.divider ? { divider: true } : { row: legRowOf(item.row) })),
    folded: layout.folded.map(legRowOf),
    foldLabel: layout.foldLabel,
  };
}

// The /teasers.json body.
//   input  {model: null (board not ready) | {legs, placed, build, plan, legRows},
//           error (a thrown model's message, or null), scannerStatus,
//           stakeSettings {bankroll, multiplier}, now}
function buildTeasersPayload(input) {
  const { model, stakeSettings } = input;
  const tickets = model ? model.plan.tickets : [];
  return {
    generatedAt: new Date(input.now).toISOString(),
    error: input.error || null,
    ready: Boolean(model),
    status: teaserView.statusText(input.scannerStatus, model ? model.legs : null),
    empty: teaserView.emptyText(model ? teaserView.planStateOf(model.plan, model.build) : null),
    ticketCount: tickets.length,
    summary: model ? teaserView.summaryView(model.plan.summary, model.placed.length, stakeSettings) : null,
    open: model ? model.placed.map(openTicketOf) : [],
    tickets: tickets.map(pendingTicketOf),
    moreLabel: tickets.length > teaserView.TICKETS_SHOWN ? teaserView.moreTicketsText(tickets.slice(teaserView.TICKETS_SHOWN)) : "",
    legs: legsOf(model ? model.legRows : []),
    breakEven: teaserView.breakEvenText(),
  };
}

module.exports = { buildTeasersPayload };
