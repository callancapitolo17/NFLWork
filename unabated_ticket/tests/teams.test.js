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
  assert.equal(teams.teamKey("cfb", "Kentucky"), "cfb:kentucky"); // seeded from Novig; "Eastern Kentucky" stays its own key
  assert.equal(teams.teamKey("cfb", "Notre Dame"), "cfb:notre-dame");
  assert.equal(teams.teamKey("cfb", "Michigan State"), null); // not seeded, and no nickname fallback
  assert.equal(teams.teamKey("cfb", "Wildcats"), null);
});

test("MLB: full names (Kalshi / Novig), codes, cities that name one team, the Athletics without a city", () => {
  assert.equal(teams.teamKey("mlb", "Miami Marlins"), "mlb:mia");
  assert.equal(teams.teamKey("mlb", "Minnesota Twins"), "mlb:min");
  assert.equal(teams.teamKey("mlb", "MIA Marlins"), "mlb:mia");
  assert.equal(teams.teamKey("mlb", "St. Louis"), "mlb:stl");
  assert.equal(teams.teamKey("mlb", "Athletics"), "mlb:ath");
  assert.equal(teams.teamKey("mlb", "Oakland Athletics"), "mlb:ath");
  assert.equal(teams.teamKey("mlb", "OAK Athletics"), "mlb:ath");
  assert.equal(teams.teamKey("mlb", "Los Angeles D"), "mlb:lad"); // KXMLBRFI event title truncation
  assert.equal(teams.teamKey("mlb", "Los Angeles"), null);
  assert.equal(teams.teamKey("mlb", "Chicago"), null); // Cubs or White Sox
  assert.equal(teams.teamKey("mlb", "New York"), null);
  assert.equal(teams.teamKey("mlb", "Yankees"), "mlb:nyy"); // nickname fallback
  assert.equal(teams.teamKey("mlb", "Chicago Bears"), null); // an NFL nickname is not an MLB team
});

test("CFB: Novig's spellings resolve and agree with Kalshi's St. forms; WNBA table", () => {
  assert.equal(teams.teamKey("cfb", "Portland State"), "cfb:portland-state");
  assert.equal(teams.teamKey("cfb", "Portland St."), "cfb:portland-state");
  assert.equal(teams.teamKey("cfb", "North Carolina A&T"), "cfb:north-carolina-a&t");
  assert.equal(teams.teamKey("cfb", "Texas A&M"), "cfb:texas-a&m");
  assert.equal(teams.teamKey("wnba", "Golden State Valkyries"), "wnba:gsv");
  assert.equal(teams.teamKey("wnba", "Valkyries"), "wnba:gsv");
  assert.equal(teams.teamKey("wnba", "Connecticut"), "wnba:conn");
});

test("unknown league or name is null, never a guess", () => {
  assert.equal(teams.teamKey("nba", "Washington"), null);
  assert.equal(teams.teamKey("nfl", "Springfield Isotopes"), null);
  assert.equal(teams.teamKey(undefined, "Chicago Bears"), null);
  assert.deepEqual(teams.knownLeagues(), ["nfl", "mlb", "wnba", "cfb"]);
});
