// Run: node --test unabated_ticket/tests/prophetx.test.js
const test = require("node:test");
const assert = require("node:assert/strict");
const prophetx = require("../server/prophetx.js");

const STRIKE = "0123456789abcdef0123456789abcdef";
const OTHER = "fedcba9876543210fedcba9876543210";

test("the quote sends the line's own strike id first, then the game's other ProphetX strike ids", () => {
  const line = { bookId: 66, sourceData: STRIKE };
  const gameLines = [
    line, { bookId: 66, sourceData: OTHER }, { bookId: 66, sourceData: null },
    { bookId: 89, sourceData: "9b2f3c3e-0000-4000-8000-000000000000" }, { bookId: 66, sourceData: "not-hex" },
  ];
  assert.deepEqual(prophetx.strikeIdsOf(line, gameLines), { lineStrikeId: STRIKE, strikeIds: [STRIKE, OTHER] });
});

test("a main line carries no strike id of its own; a game with no ProphetX rung sends none", () => {
  const main = { bookId: 66, sourceData: null };
  assert.deepEqual(prophetx.strikeIdsOf(main, [main, { bookId: 66, sourceData: OTHER }]), { lineStrikeId: null, strikeIds: [OTHER] });
  assert.deepEqual(prophetx.strikeIdsOf(main, [main]), { lineStrikeId: null, strikeIds: [] });
});

test("at most MAX_STRIKE_IDS ids go out", () => {
  const many = Array.from({ length: 30 }, (_, index) => ({ bookId: 66, sourceData: index.toString(16).padStart(32, "0") }));
  assert.equal(prophetx.strikeIdsOf({ bookId: 66 }, many).strikeIds.length, prophetx.MAX_STRIKE_IDS);
});
