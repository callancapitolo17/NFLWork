// Run: node --test unabated_ticket/tests
// teams.js: the runtime team index (from Unabated's snapshot team lists) and
// the venue-spelling resolution rules, on fixtures/teams_index.json.
const test = require("node:test");
const assert = require("node:assert/strict");
const fs = require("node:fs");
const path = require("node:path");
const teams = require("../extension/teams.js");

const INDEX = JSON.parse(fs.readFileSync(path.join(__dirname, "fixtures", "teams_index.json"), "utf8")).leagues;
teams.loadIndex(INDEX);
const idOf = (league, name) => INDEX[league].find((t) => t.name === name).id;
const keyOf = (league, name) => `${league}:${idOf(league, name)}`;

test("normalizeName: case, periods, spaces, trailing St. -> State", () => {
  assert.equal(teams.normalizeName("  Missouri  St. "), "missouri state");
  assert.equal(teams.normalizeName("Missouri State"), "missouri state");
  assert.equal(teams.normalizeName("St. John's"), "st john's");
  assert.equal(teams.normalizeName("Texas A&M"), "texas a&m");
  assert.equal(teams.normalizeName(""), null);
  assert.equal(teams.normalizeName(null), null);
});

test("index: every league of the fixture registers; keys are league:id; exportIndex round-trips", () => {
  assert.deepEqual(teams.knownLeagues().sort(), ["cfb", "mlb", "nba", "nfl", "nhl", "wnba"]);
  assert.equal(teams.teamCount("cfb"), 265);
  assert.equal(teams.teamKey("nfl", "Carolina Panthers"), keyOf("nfl", "Carolina Panthers"));
  assert.equal(teams.keyOf("nfl", 5), "nfl:5");
  const exported = teams.exportIndex();
  assert.equal(exported.nfl.length, 32);
  assert.deepEqual(Object.keys(exported.nfl[0]).sort(), ["abbreviation", "id", "name"]);
});

test("NFL: Unabated full names and every Kalshi spelling seen on the wire", () => {
  assert.equal(teams.teamKey("nfl", "Philadelphia Eagles"), keyOf("nfl", "Philadelphia Eagles"));
  assert.equal(teams.teamKey("nfl", "PIT Steelers"), keyOf("nfl", "Pittsburgh Steelers")); // KXNFLGAME event title: code + nickname
  assert.equal(teams.teamKey("nfl", "NE Patriots"), keyOf("nfl", "New England Patriots"));
  assert.equal(teams.teamKey("nfl", "New England"), keyOf("nfl", "New England Patriots")); // market title "New England wins"
  assert.equal(teams.teamKey("nfl", "DEN Broncos"), keyOf("nfl", "Denver Broncos"));
  assert.equal(teams.teamKey("nfl", "Detroit"), keyOf("nfl", "Detroit Lions")); // KXNFLSPREAD event title, city only
  assert.equal(teams.teamKey("nfl", "Rams"), keyOf("nfl", "Los Angeles Rams")); // nickname alone is unique
  assert.equal(teams.teamKey("nfl", "Los Angeles"), null); // Rams or Chargers
  assert.equal(teams.teamKey("nfl", "New York"), null);
});

test("MLB: Kalshi codes and truncations, Novig full names, the Athletics", () => {
  assert.equal(teams.teamKey("mlb", "Miami Marlins"), keyOf("mlb", "Miami Marlins"));
  assert.equal(teams.teamKey("mlb", "WSH Nationals"), keyOf("mlb", "Washington Nationals")); // Kalshi code + nickname
  assert.equal(teams.teamKey("mlb", "Los Angeles D"), keyOf("mlb", "Los Angeles Dodgers")); // KXMLBRFI title truncation (alias)
  assert.equal(teams.teamKey("mlb", "Athletics"), keyOf("mlb", "Oakland Athletics"));
  assert.equal(teams.teamKey("mlb", "OAK Athletics"), keyOf("mlb", "Oakland Athletics"));
  assert.equal(teams.teamKey("mlb", "St. Louis"), keyOf("mlb", "St. Louis Cardinals"));
  assert.equal(teams.teamKey("mlb", "Chicago"), null); // Cubs or White Sox
});

test("CFB: Kalshi St. forms and Novig spellings both land on Unabated's team; ambiguity is null", () => {
  assert.equal(teams.teamKey("cfb", "Missouri St."), keyOf("cfb", "Missouri State"));
  assert.equal(teams.teamKey("cfb", "Penn St."), keyOf("cfb", "Penn State"));
  assert.equal(teams.teamKey("cfb", "Texas A&M"), keyOf("cfb", "Texas A&M"));
  assert.equal(teams.teamKey("cfb", "South Alabama"), keyOf("cfb", "South Alabama")); // the Tulane opponent the hand table lacked
  assert.equal(teams.teamKey("cfb", "Grambling St."), keyOf("cfb", "Grambling")); // query starts with the team
  assert.equal(teams.teamKey("cfb", "Middle Tennessee"), keyOf("cfb", "Middle Tennessee State")); // team starts with the query
  assert.equal(teams.teamKey("cfb", "Southern Mississippi"), keyOf("cfb", "Southern Miss")); // alias: "Southern" is its own team
  assert.equal(teams.teamKey("cfb", "Southern Mississippi State"), null); // no rule extends "Southern Miss" by "issippi State"
  assert.equal(teams.teamKey("cfb", "Tarleton State"), keyOf("cfb", "Tarleton")); // institutional suffix
  assert.equal(teams.teamKey("cfb", "Lindenwood"), keyOf("cfb", "Lindenwood University"));
  assert.equal(teams.teamKey("cfb", "Long Island University"), keyOf("cfb", "Long Island"));
  assert.equal(teams.teamKey("cfb", "UAlbany"), keyOf("cfb", "Albany")); // alias
  assert.equal(teams.teamKey("cfb", "North Carolina State"), keyOf("cfb", "NC State")); // alias
  assert.equal(teams.teamKey("cfb", "Kentucky"), keyOf("cfb", "Kentucky")); // exact beats "Eastern Kentucky"
  assert.equal(teams.teamKey("cfb", "Miami"), null); // Miami Florida or Miami Ohio
  assert.equal(teams.teamKey("cfb", "Southern"), keyOf("cfb", "Southern")); // Southern University, an exact name
  assert.equal(teams.teamKey("cfb", "Springfield Isotopes"), null);
});

test("WNBA / NBA / NHL come from the same snapshots", () => {
  assert.equal(teams.teamKey("wnba", "Golden State Valkyries"), keyOf("wnba", "Golden State Valkyries"));
  assert.equal(teams.teamKey("wnba", "Valkyries"), keyOf("wnba", "Golden State Valkyries"));
  assert.equal(teams.teamKey("nba", "LA Clippers"), keyOf("nba", "LA Clippers"));
  assert.equal(teams.teamKey("nhl", "Vegas Golden Knights"), keyOf("nhl", "Vegas Golden Knights"));
});

test("unknown league or name is null, never a guess; registering again refreshes, never duplicates", () => {
  assert.equal(teams.teamKey("soccer", "Arsenal"), null);
  assert.equal(teams.teamKey(undefined, "Chicago Bears"), null);
  teams.registerTeams("nfl", [{ id: 5, name: "Carolina Panthers", abbreviation: "CAR" }]);
  assert.equal(teams.teamCount("nfl"), 32);
  teams.registerTeams("nfl", [{ id: 9999, name: "Test Team" }, { id: null, name: "skipped" }, { id: 7, name: 12 }]);
  assert.equal(teams.teamCount("nfl"), 33);
  assert.equal(teams.teamKey("nfl", "Test Team"), "nfl:9999");
});
