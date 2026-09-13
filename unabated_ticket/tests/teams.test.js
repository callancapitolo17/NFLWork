// Run: node --test unabated_ticket/tests
// teams.js: the runtime team index (from Unabated's snapshot team lists plus
// the eventName spelling of each team, #118) and the venue-spelling
// resolution rules, on fixtures/teams_index.json.
const test = require("node:test");
const assert = require("node:assert/strict");
const fs = require("node:fs");
const path = require("node:path");
const teams = require("../extension/teams.js");

const INDEX = JSON.parse(fs.readFileSync(path.join(__dirname, "fixtures", "teams_index.json"), "utf8")).leagues;
teams.loadIndex(INDEX);
const idOf = (league, name) => INDEX[league].find((t) => t.name === name).id;
const keyOf = (league, name) => `${league}:${idOf(league, name)}`;
const eventNameOf = (league, name) => INDEX[league].find((t) => t.name === name).eventName;
// The CFB index under a league path with NO aliases: what the rules alone
// resolve, so a test can show which spellings still need an ALIASES row.
const NO_ALIAS_LEAGUE = "cfb_rules_only";
teams.registerTeams(NO_ALIAS_LEAGUE, INDEX.cfb);
const rulesOnlyKey = (name) => teams.teamKey(NO_ALIAS_LEAGUE, name);
const rulesOnlyKeyOf = (name) => `${NO_ALIAS_LEAGUE}:${idOf("cfb", name)}`;

test("normalizeName: case, periods, spaces, trailing St. -> State", () => {
  assert.equal(teams.normalizeName("  Missouri  St. "), "missouri state");
  assert.equal(teams.normalizeName("Missouri State"), "missouri state");
  assert.equal(teams.normalizeName("St. John's"), "st john's");
  assert.equal(teams.normalizeName("Texas A&M"), "texas a&m");
  assert.equal(teams.normalizeName(""), null);
  assert.equal(teams.normalizeName(null), null);
});

test("index: every league of the fixture registers; keys are league:id; exportIndex round-trips", () => {
  assert.deepEqual(teams.knownLeagues().filter((l) => l !== NO_ALIAS_LEAGUE && l !== "probe").sort(), ["cfb", "mlb", "nba", "nfl", "nhl", "wnba"]);
  assert.equal(teams.teamCount("cfb"), 265);
  assert.equal(teams.teamKey("nfl", "Carolina Panthers"), keyOf("nfl", "Carolina Panthers"));
  assert.equal(teams.keyOf("nfl", 5), "nfl:5");
  const exported = teams.exportIndex();
  assert.equal(exported.nfl.length, 32);
  // The persisted entry carries the eventName spelling, or loadIndex would lose it next session.
  assert.deepEqual(Object.keys(exported.nfl[0]).sort(), ["abbreviation", "eventName", "id", "name"]);
  assert.equal(exported.nfl.find((t) => t.id === 5).eventName, "Panthers Carolina");
});

test("eventName spellings: both names of a team resolve to the same key, and the count tells a new spelling from a refresh", () => {
  assert.equal(eventNameOf("nfl", "Baltimore Ravens"), "Ravens Baltimore");
  assert.equal(teams.teamKey("nfl", "Ravens Baltimore"), keyOf("nfl", "Baltimore Ravens"));
  assert.equal(eventNameOf("cfb", "Prairie View"), "Prairie View A&M Panthers");
  assert.equal(teams.teamKey("cfb", "Prairie View A&M Panthers"), keyOf("cfb", "Prairie View"));
  // Teams with no game row this session have no second spelling and still resolve by name.
  assert.equal(eventNameOf("nhl", "Vegas Golden Knights"), null);
  assert.equal(teams.teamKey("nhl", "Vegas Golden Knights"), keyOf("nhl", "Vegas Golden Knights"));
  const before = teams.spellingCount();
  teams.registerTeams("nhl", [{ id: idOf("nhl", "Vegas Golden Knights"), name: "Vegas Golden Knights", abbreviation: "VGK" }]);
  assert.equal(teams.spellingCount(), before); // same team, same spellings: a no-op refresh
  teams.registerTeams("nhl", [{ id: idOf("nhl", "Vegas Golden Knights"), name: "Vegas Golden Knights", abbreviation: "VGK", eventName: "Golden Knights Vegas" }]);
  assert.equal(teams.spellingCount(), before + 1);
  assert.equal(teams.teamKey("nhl", "Golden Knights Vegas"), keyOf("nhl", "Vegas Golden Knights"));
  // A refresh WITHOUT the spelling (the team's row went live) keeps the one held.
  teams.registerTeams("nhl", [{ id: idOf("nhl", "Vegas Golden Knights"), name: "Vegas Golden Knights", abbreviation: "VGK", eventName: null }]);
  assert.equal(teams.exportIndex().nhl.find((t) => t.name === "Vegas Golden Knights").eventName, "Golden Knights Vegas");
  assert.equal(teams.spellingCount(), before + 1);
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
  assert.equal(teams.teamKey("cfb", "Lafayette"), keyOf("cfb", "Lafayette")); // exact still beats the UL Lafayette alias target
  // Texas A&M is not Texas: the A&M alias must stay one spelling, never a rule.
  assert.equal(teams.teamKey("cfb", "Texas A&M"), keyOf("cfb", "Texas A&M"));
  assert.equal(teams.teamKey("cfb", "Alabama A&M"), keyOf("cfb", "Alabama A&M"));
  assert.equal(teams.teamKey("cfb", "Kentucky"), keyOf("cfb", "Kentucky")); // exact beats "Eastern Kentucky"
  assert.equal(teams.teamKey("cfb", "Miami"), null); // Miami Florida or Miami Ohio
  assert.equal(teams.teamKey("cfb", "Southern"), keyOf("cfb", "Southern")); // Southern University, an exact name
  assert.equal(teams.teamKey("cfb", "Springfield Isotopes"), null);
});

test("CFB eventName spellings (#118): the retired aliases resolve by rule alone; the kept aliases still cannot", () => {
  // Retired 2026-09-12 — each resolves with no ALIASES row (NO_ALIAS_LEAGUE has none).
  assert.equal(eventNameOf("cfb", "Prairie View"), "Prairie View A&M Panthers");
  assert.equal(rulesOnlyKey("Prairie View A&M"), rulesOnlyKeyOf("Prairie View"));
  assert.equal(eventNameOf("cfb", "SE Louisiana"), "Southeastern Louisiana Lions");
  assert.equal(rulesOnlyKey("Southeastern Louisiana"), rulesOnlyKeyOf("SE Louisiana"));
  // Kept — the rules alone are null and the alias is still what resolves the venue spelling.
  assert.equal(eventNameOf("cfb", "Albany"), "Albany Great Danes"); // nothing derives "UAlbany"
  assert.equal(rulesOnlyKey("UAlbany"), null);
  assert.equal(teams.teamKey("cfb", "UAlbany"), keyOf("cfb", "Albany"));
  assert.equal(eventNameOf("cfb", "Southern Miss"), "Southern Miss Golden Eagles"); // nothing turns Mississippi into Miss
  assert.equal(rulesOnlyKey("Southern Mississippi"), null);
  assert.equal(teams.teamKey("cfb", "Southern Mississippi"), keyOf("cfb", "Southern Miss"));
  assert.equal(eventNameOf("cfb", "NC State"), "North Carolina State Wolfpack"); // starts with the query, but so does "North Carolina" + State
  assert.equal(rulesOnlyKey("North Carolina State"), null);
  assert.equal(teams.teamKey("cfb", "North Carolina State"), keyOf("cfb", "NC State"));
  assert.equal(eventNameOf("cfb", "UL Lafayette"), "Louisiana-Lafayette Ragin' Cajuns"); // "Louisiana" hits Louisiana Tech (both forms) and SE Louisiana
  assert.equal(eventNameOf("cfb", "Louisiana Tech"), "Louisiana Tech Bulldogs");
  assert.equal(eventNameOf("cfb", "UL Monroe"), "Louisiana-Monroe Warhawks");
  assert.equal(rulesOnlyKey("Louisiana"), null);
  assert.equal(teams.teamKey("cfb", "Louisiana"), keyOf("cfb", "UL Lafayette"));
  // Venue spellings the long forms newly resolve (each verified by hand 2026-09-12).
  assert.equal(teams.teamKey("cfb", "Tennessee-Martin"), keyOf("cfb", "UT Martin")); // "Tennessee-Martin Skyhawks"
  assert.equal(teams.teamKey("cfb", "Miami (OH)"), keyOf("cfb", "Miami Ohio")); // "Miami (OH) RedHawks"
  assert.equal(teams.teamKey("cfb", "Southeast Missouri State"), keyOf("cfb", "SE Missouri State"));
  assert.equal(teams.teamKey("cfb", "Eastern Washington"), keyOf("cfb", "East Washington"));
  assert.equal(teams.teamKey("cfb", "UCF"), keyOf("cfb", "Central Florida")); // "UCF Knights"
  // Safety probes: the extra spellings add candidates, never a guess.
  assert.equal(teams.teamKey("cfb", "Miami"), null); // Miami Florida / Miami Ohio, plus both long forms
  assert.equal(teams.teamKey("cfb", "Carolina"), null); // North / South / East / Western / Coastal Carolina
  assert.equal(teams.teamKey("cfb", "Texas"), keyOf("cfb", "Texas")); // exact beats "Texas Longhorns" and Texas A&M / Tech / State
  assert.equal(teams.teamKey("cfb", "Southern"), keyOf("cfb", "Southern")); // exact beats "Southern University Jaguars"
  assert.equal(teams.teamKey("cfb", "Texas A&M"), keyOf("cfb", "Texas A&M"));
  assert.equal(teams.teamKey("cfb", "Alabama A&M"), keyOf("cfb", "Alabama A&M"));
  for (const nickname of ["Tigers", "Bulldogs", "Panthers"]) assert.equal(teams.teamKey("cfb", nickname), null, nickname);
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

test("a list name always beats a colliding eventName spelling, whichever registers first", () => {
  teams.registerTeams("probe", [{ id: 1, name: "Alpha", abbreviation: null, eventName: "Beta" }, { id: 2, name: "Beta", abbreviation: null, eventName: null }]);
  assert.equal(teams.teamKey("probe", "Beta"), "probe:2"); // Beta's own name overwrote Alpha's spelling
  teams.registerTeams("probe", [{ id: 3, name: "Gamma", abbreviation: null, eventName: "Beta" }]);
  assert.equal(teams.teamKey("probe", "Beta"), "probe:2"); // and a later colliding spelling is not registered
  assert.equal(teams.teamKey("probe", "Gamma"), "probe:3");
});
