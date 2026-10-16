// Run: node --test unabated_ticket/tests/kalshiplace.test.js
const test = require("node:test");
const assert = require("node:assert/strict");
const kalshiPlace = require("../extension/kalshiplace.js");

const LINE = { eventId: 123742, periodTypeId: 1, betTypeId: 2, sideIndex: 1, points: -35.5, price: -143, sourceFormat: 4, sourcePrice: 0.587157 };

test("only a pregame Kalshi line gets the button", () => {
  const ids = { eventId: 1, periodTypeId: 1, betTypeId: 3, sideIndex: 0 };
  assert.equal(kalshiPlace.canPlace({ ...ids, bookId: 105, live: null }), true);
  assert.equal(kalshiPlace.canPlace({ ...ids, bookId: 105, live: { period: 2 } }), false);
  assert.equal(kalshiPlace.canPlace({ ...ids, bookId: 100, live: null }), false);
  assert.equal(kalshiPlace.canPlace({ ...ids, eventId: null, bookId: 105 }), false);
});

test("the amount box reads dollars, commas and a $ allowed, cents at most", () => {
  assert.equal(kalshiPlace.parseStake("$1,234.50"), 1234.5);
  assert.equal(kalshiPlace.parseStake(" 250 "), 250);
  assert.equal(kalshiPlace.parseStake("0"), null);
  assert.equal(kalshiPlace.parseStake("12.345"), null);
  assert.equal(kalshiPlace.parseStake("-5"), null);
  assert.equal(kalshiPlace.stakeText(281.734), "281.73");
  assert.equal(kalshiPlace.stakeText(0), "");
});

test("the request carries the line, the typed stake, the Ticket's all-in price and the Kalshi read", () => {
  const { body } = kalshiPlace.placeRequestOf({ line: LINE, stakeInput: "150", kalshiReadAt: "2026-10-16T21:00:00Z", clickId: "abc12345" });
  assert.deepEqual(body, {
    clickId: "abc12345", eventId: 123742, periodTypeId: 1, betTypeId: 2, sideIndex: 1, points: -35.5,
    stake: 150, ticketAllIn: 0.587157, kalshiReadAt: "2026-10-16T21:00:00Z",
  });
  // An American-only price is its implied probability; a moneyline sends points null.
  const moneyline = kalshiPlace.placeRequestOf({ line: { ...LINE, betTypeId: 1, points: undefined, sourceFormat: 1, sourcePrice: null, price: 150 }, stakeInput: "10", kalshiReadAt: "x", clickId: "abc12345" });
  assert.deepEqual([moneyline.body.points, moneyline.body.ticketAllIn], [null, 0.4]);
});

test("no amount or no Kalshi read is an error in words, nothing to send", () => {
  assert.match(kalshiPlace.placeRequestOf({ line: LINE, stakeInput: "abc", kalshiReadAt: "x", clickId: "abc12345" }).error, /not an amount/);
  assert.match(kalshiPlace.placeRequestOf({ line: LINE, stakeInput: "10", kalshiReadAt: null, clickId: "abc12345" }).error, /not been read/);
});

test("reply tones", () => {
  assert.deepEqual(["filled", "partial", "unfilled", "rejected", "unknown", "refused"].map(kalshiPlace.replyTone),
    ["ok", "ok", "bad", "bad", "warn", "refused"]);
});
