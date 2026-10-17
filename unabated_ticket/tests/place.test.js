// Run: node --test unabated_ticket/tests/place.test.js
const test = require("node:test");
const assert = require("node:assert/strict");
const place = require("../extension/place.js");

const LINE = { bookId: 105, eventId: 123742, periodTypeId: 1, betTypeId: 2, sideIndex: 1, points: -35.5, price: -143, sourceFormat: 4, sourcePrice: 0.587157 };
const NOVIG_LINE = { ...LINE, bookId: 89, price: -208, sourcePrice: 0.675 };
const READS = { kalshi: "2026-10-16T21:00:00Z", novig: "2026-10-16T21:00:10Z" };

test("only a pregame Kalshi line, or a pregame Novig spread or total, gets the button", () => {
  const ids = { eventId: 1, periodTypeId: 1, betTypeId: 3, sideIndex: 0 };
  assert.equal(place.canPlace({ ...ids, bookId: 105, live: null }), true);
  assert.equal(place.canPlace({ ...ids, bookId: 105, live: { period: 2 } }), false);
  assert.equal(place.canPlace({ ...ids, bookId: 100, live: null }), false);
  assert.equal(place.canPlace({ ...ids, eventId: null, bookId: 105 }), false);
  assert.equal(place.canPlace({ ...ids, betTypeId: 1, bookId: 105, live: null }), true);
  assert.equal(place.canPlace({ ...ids, bookId: 89, live: null }), true);
  // Unabated carries no Novig outcome id on a moneyline.
  assert.equal(place.canPlace({ ...ids, betTypeId: 1, bookId: 89, live: null }), false);
  assert.equal(place.canPlace({ ...ids, bookId: 89, live: { period: 2 } }), false);
});

test("the amount box reads dollars, commas and a $ allowed, cents at most", () => {
  assert.equal(place.parseStake("$1,234.50"), 1234.5);
  assert.equal(place.parseStake(" 250 "), 250);
  assert.equal(place.parseStake("0"), null);
  assert.equal(place.parseStake("12.345"), null);
  assert.equal(place.parseStake("-5"), null);
  assert.equal(place.stakeText(281.734), "281.73");
  assert.equal(place.stakeText(0), "");
});

test("a Kalshi request carries the line, the typed stake, the Ticket's all-in price and the Kalshi read", () => {
  const request = place.placeRequestOf({ line: LINE, stakeInput: "150", readAts: READS, clickId: "abc12345" });
  assert.equal(request.path, "/place_kalshi.json");
  assert.deepEqual(request.body, {
    clickId: "abc12345", eventId: 123742, periodTypeId: 1, betTypeId: 2, sideIndex: 1, points: -35.5,
    stake: 150, ticketAllIn: 0.587157, kalshiReadAt: "2026-10-16T21:00:00Z",
  });
  // An American-only price is its implied probability; a moneyline sends points null.
  const moneyline = place.placeRequestOf({ line: { ...LINE, betTypeId: 1, points: undefined, sourceFormat: 1, sourcePrice: null, price: 150 }, stakeInput: "10", readAts: READS, clickId: "abc12345" });
  assert.deepEqual([moneyline.body.points, moneyline.body.ticketAllIn], [null, 0.4]);
});

test("a Novig request goes to its own route with the Novig read and the Kalshi read", () => {
  const request = place.placeRequestOf({ line: NOVIG_LINE, stakeInput: "50", readAts: READS, clickId: "abc12345" });
  assert.equal(request.path, "/place_novig.json");
  assert.deepEqual(request.body, {
    clickId: "abc12345", eventId: 123742, periodTypeId: 1, betTypeId: 2, sideIndex: 1, points: -35.5,
    stake: 50, ticketAllIn: 0.675, novigReadAt: "2026-10-16T21:00:10Z", kalshiReadAt: "2026-10-16T21:00:00Z",
  });
  const noKalshi = place.placeRequestOf({ line: NOVIG_LINE, stakeInput: "50", readAts: { novig: READS.novig, kalshi: null }, clickId: "abc12345" });
  assert.equal(noKalshi.body.kalshiReadAt, null);
});

test("no amount or no read at the venue is an error in words, nothing to send", () => {
  assert.match(place.placeRequestOf({ line: LINE, stakeInput: "abc", readAts: READS, clickId: "abc12345" }).error, /not an amount/);
  assert.match(place.placeRequestOf({ line: LINE, stakeInput: "10", readAts: { ...READS, kalshi: null }, clickId: "abc12345" }).error, /Kalshi bets have not been read/);
  assert.match(place.placeRequestOf({ line: NOVIG_LINE, stakeInput: "10", readAts: { ...READS, novig: null }, clickId: "abc12345" }).error, /Novig bets have not been read/);
  assert.match(place.placeRequestOf({ line: { ...LINE, bookId: 100 }, stakeInput: "10", readAts: READS, clickId: "abc12345" }).error, /no Place button/);
});

test("reply tones", () => {
  assert.deepEqual(["filled", "partial", "unfilled", "rejected", "unknown", "refused"].map(place.replyTone),
    ["ok", "ok", "bad", "bad", "warn", "refused"]);
});

test("Novig's read is the oldest of the app, Place and (when read) bot reads, null until both required land", () => {
  const at = (iso) => ({ fetchedAt: iso });
  const sources = { kalshi: at("2026-10-16T21:00:00Z"), novig: at("2026-10-16T21:00:10Z"), novig_place: at("2026-10-16T20:59:50Z") };
  assert.deepEqual(place.readAtsOf(sources), { kalshi: "2026-10-16T21:00:00Z", novig: "2026-10-16T20:59:50Z" });
  assert.equal(place.readAtsOf({ ...sources, novig_bot: at("2026-10-16T20:58:00Z") }).novig, "2026-10-16T20:58:00Z");
  assert.equal(place.readAtsOf({ ...sources, novig_bot: at(null) }).novig, null, "a bot read that never landed counts as missing");
  assert.equal(place.readAtsOf({ kalshi: sources.kalshi, novig: sources.novig }).novig, null);
  assert.deepEqual(place.readAtsOf({}), { kalshi: null, novig: null });
});
