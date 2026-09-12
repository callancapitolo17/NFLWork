// Run: node --test unabated_ticket/tests
const test = require("node:test");
const assert = require("node:assert/strict");
const teams = require("../extension/teams.js");

test("normalizeName: case, periods, spaces, trailing St. -> State", () => {
  assert.equal(teams.normalizeName("  Missouri  St. "), "missouri state");
  assert.equal(teams.normalizeName("Missouri State"), "missouri state");
  assert.equal(teams.normalizeName("St. John's"), "st john's");
  assert.equal(teams.normalizeName("Texas A&M"), "texas a&m");
  assert.equal(teams.normalizeName(""), null);
  assert.equal(teams.normalizeName(null), null);
});

test("NFL: Unabated full names and every Kalshi spelling seen on the wire", () => {
  assert.equal(teams.teamKey("nfl", "Philadelphia Eagles"), "nfl:phi");
  assert.equal(teams.teamKey("nfl", "Carolina Panthers"), "nfl:car");
  assert.equal(teams.teamKey("nfl", "Chicago Bears"), "nfl:chi");
  assert.equal(teams.teamKey("nfl", "PIT Steelers"), "nfl:pit"); // KXNFLGAME event title
  assert.equal(teams.teamKey("nfl", "NE Patriots"), "nfl:ne");
  assert.equal(teams.teamKey("nfl", "New England"), "nfl:ne"); // market title "New England wins"
  assert.equal(teams.teamKey("nfl", "DEN Broncos"), "nfl:den"); // KXNFL1HTOTAL event title
  assert.equal(teams.teamKey("nfl", "KC Chiefs"), "nfl:kc");
  assert.equal(teams.teamKey("nfl", "Detroit"), "nfl:det"); // KXNFLSPREAD event title
  assert.equal(teams.teamKey("nfl", "Buffalo"), "nfl:buf");
});

test("NFL: nickname fallback lands an unexpected abbreviation; ambiguous cities stay null", () => {
  assert.equal(teams.teamKey("nfl", "LA Rams"), "nfl:lar");
  assert.equal(teams.teamKey("nfl", "SF 49ers"), "nfl:sf");
  assert.equal(teams.teamKey("nfl", "New York"), null);
  assert.equal(teams.teamKey("nfl", "Los Angeles"), null);
});

test("CFB: Kalshi short names, State/St. agree, no nickname fallback", () => {
  assert.equal(teams.teamKey("cfb", "Texas A&M"), "cfb:texas-a&m");
  assert.equal(teams.teamKey("cfb", "Missouri St."), "cfb:missouri-state");
  assert.equal(teams.teamKey("cfb", "Missouri State"), "cfb:missouri-state");
  assert.equal(teams.teamKey("cfb", "Penn St."), "cfb:penn-state");
  assert.equal(teams.teamKey("cfb", "Eastern Kentucky"), "cfb:eastern-kentucky");
  assert.equal(teams.teamKey("cfb", "Kentucky"), null);
  assert.equal(teams.teamKey("cfb", "Notre Dame"), "cfb:notre-dame");
  assert.equal(teams.teamKey("cfb", "Michigan State"), null);
});

test("unknown league or name is null, never a guess", () => {
  assert.equal(teams.teamKey("mlb", "Washington"), null);
  assert.equal(teams.teamKey("nfl", "Springfield Isotopes"), null);
  assert.equal(teams.teamKey(undefined, "Chicago Bears"), null);
  assert.deepEqual(teams.knownLeagues(), ["nfl", "cfb"]);
});
