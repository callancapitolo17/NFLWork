// Team-name -> stable key tables for matching a bet from one venue against a
// line on Unabated. Pure: no DOM, no fetch, no chrome.* — loaded as a plain
// <script> in panel.html (exposes globalThis.UnabatedTeams, before bets.js)
// and via require() in tests/teams.test.js.
//
// Input:  a league path as feed.LEAGUES uses it ("nfl", "cfb", "mlb", ...) and
//         a team name as some venue wrote it.
// Output: a key like "nfl:ne" / "cfb:texas-a&m", or null when the name is not
//         in the table. Null is the ONLY answer for an unknown name — a wrong
//         key would match a bet to the wrong game, and a null shows up in the
//         panel's unmatched list with the raw name so the table can grow.
//
// Spellings seeded (2026-09-11):
//   NFL  Unabated full names ("Philadelphia Eagles"); Kalshi event titles use
//        "PIT Steelers vs NE Patriots" (KXNFLGAME), "DEN Broncos vs KC Chiefs"
//        (KXNFL1H*), and city only "Detroit vs Buffalo" (KXNFLSPREAD/TOTAL);
//        Kalshi market titles say "New England wins". Only PIT/NE/DEN/KC/DET/
//        BUF abbreviations are seen on the wire; the rest are the standard
//        ones, and every NFL name also resolves by its nickname alone (unique
//        within the league), so an unexpected abbreviation still lands.
//   CFB  the Kalshi short names in tests/fixtures/bets/kalshi_fixture.json plus
//        a few from the account's other fills ("Penn St.", "Lehigh", ...).
//        Unabated's CFB snapshot was not captured; "Missouri State" and
//        "Missouri St." normalise to the same key. No nickname fallback: CFB
//        nicknames repeat ("Eastern Kentucky" must never become "Kentucky").
//   MLB  full names as kalshi_common/leg_types._MLB_CODE_TO_TEAM spells them
//        (Kalshi codes ported, ATH/OAK both = "Athletics"); Novig's
//        game.awayTeam.name / homeTeam.name are full names too (#116), and
//        nicknames are unique in MLB so the nickname fallback applies.

(function (root) {
  "use strict";

  // [key, Unabated full name, city, Kalshi abbreviation]. City is used as an
  // alias only when it names one team in the league (not "New York",
  // "Los Angeles").
  const NFL_TEAMS = [
    ["ari", "Arizona Cardinals", "Arizona", "ARI"],
    ["atl", "Atlanta Falcons", "Atlanta", "ATL"],
    ["bal", "Baltimore Ravens", "Baltimore", "BAL"],
    ["buf", "Buffalo Bills", "Buffalo", "BUF"],
    ["car", "Carolina Panthers", "Carolina", "CAR"],
    ["chi", "Chicago Bears", "Chicago", "CHI"],
    ["cin", "Cincinnati Bengals", "Cincinnati", "CIN"],
    ["cle", "Cleveland Browns", "Cleveland", "CLE"],
    ["dal", "Dallas Cowboys", "Dallas", "DAL"],
    ["den", "Denver Broncos", "Denver", "DEN"],
    ["det", "Detroit Lions", "Detroit", "DET"],
    ["gb", "Green Bay Packers", "Green Bay", "GB"],
    ["hou", "Houston Texans", "Houston", "HOU"],
    ["ind", "Indianapolis Colts", "Indianapolis", "IND"],
    ["jax", "Jacksonville Jaguars", "Jacksonville", "JAX"],
    ["kc", "Kansas City Chiefs", "Kansas City", "KC"],
    ["lv", "Las Vegas Raiders", "Las Vegas", "LV"],
    ["lac", "Los Angeles Chargers", null, "LAC"],
    ["lar", "Los Angeles Rams", null, "LAR"],
    ["mia", "Miami Dolphins", "Miami", "MIA"],
    ["min", "Minnesota Vikings", "Minnesota", "MIN"],
    ["ne", "New England Patriots", "New England", "NE"],
    ["no", "New Orleans Saints", "New Orleans", "NO"],
    ["nyg", "New York Giants", null, "NYG"],
    ["nyj", "New York Jets", null, "NYJ"],
    ["phi", "Philadelphia Eagles", "Philadelphia", "PHI"],
    ["pit", "Pittsburgh Steelers", "Pittsburgh", "PIT"],
    ["sf", "San Francisco 49ers", "San Francisco", "SF"],
    ["sea", "Seattle Seahawks", "Seattle", "SEA"],
    ["tb", "Tampa Bay Buccaneers", "Tampa Bay", "TB"],
    ["ten", "Tennessee Titans", "Tennessee", "TEN"],
    ["was", "Washington Commanders", "Washington", "WAS"],
  ];

  // [key, full name, city, Kalshi code]; "Athletics" carries no city (Oakland /
  // Sacramento / Las Vegas depending on the year) and answers to OAK and ATH.
  const MLB_TEAMS = [
    ["ari", "Arizona Diamondbacks", "Arizona", "ARI"],
    ["atl", "Atlanta Braves", "Atlanta", "ATL"],
    ["bal", "Baltimore Orioles", "Baltimore", "BAL"],
    ["bos", "Boston Red Sox", "Boston", "BOS"],
    ["chc", "Chicago Cubs", null, "CHC"],
    ["cws", "Chicago White Sox", null, "CWS"],
    ["cin", "Cincinnati Reds", "Cincinnati", "CIN"],
    ["cle", "Cleveland Guardians", "Cleveland", "CLE"],
    ["col", "Colorado Rockies", "Colorado", "COL"],
    ["det", "Detroit Tigers", "Detroit", "DET"],
    ["hou", "Houston Astros", "Houston", "HOU"],
    ["kc", "Kansas City Royals", "Kansas City", "KC"],
    ["laa", "Los Angeles Angels", null, "LAA"],
    ["lad", "Los Angeles Dodgers", null, "LAD"],
    ["mia", "Miami Marlins", "Miami", "MIA"],
    ["mil", "Milwaukee Brewers", "Milwaukee", "MIL"],
    ["min", "Minnesota Twins", "Minnesota", "MIN"],
    ["nym", "New York Mets", null, "NYM"],
    ["nyy", "New York Yankees", null, "NYY"],
    ["ath", "Athletics", null, "ATH"],
    ["phi", "Philadelphia Phillies", "Philadelphia", "PHI"],
    ["pit", "Pittsburgh Pirates", "Pittsburgh", "PIT"],
    ["sd", "San Diego Padres", "San Diego", "SD"],
    ["sf", "San Francisco Giants", "San Francisco", "SF"],
    ["sea", "Seattle Mariners", "Seattle", "SEA"],
    ["stl", "St. Louis Cardinals", "St. Louis", "STL"],
    ["tb", "Tampa Bay Rays", "Tampa Bay", "TB"],
    ["tex", "Texas Rangers", "Texas", "TEX"],
    ["tor", "Toronto Blue Jays", "Toronto", "TOR"],
    ["was", "Washington Nationals", "Washington", "WAS"],
  ];
  // Other spellings of the same MLB team: [alias, key]. "Los Angeles D" is
  // how the Kalshi KXMLBRFI event title truncates the Dodgers (fixture);
  // the Angels' truncation has not been seen on the wire, so it is not guessed.
  const MLB_ALIASES = [["Oakland Athletics", "ath"], ["Sacramento Athletics", "ath"], ["OAK Athletics", "ath"], ["AZ Diamondbacks", "ari"], ["WSH Nationals", "was"], ["Los Angeles D", "lad"]];

  // Kalshi short names; the key is the slug of the normalised name.
  const CFB_TEAMS = [
    "Alabama", "Chattanooga", "Drake", "East Carolina", "Eastern Kentucky", "Georgetown",
    "Grambling St.", "Lehigh", "Lindenwood", "Louisville", "Marshall", "Michigan",
    "Missouri St.", "Montana", "Notre Dame", "Oklahoma", "Penn St.", "Rice", "TCU",
    "Texas A&M", "Villanova",
  ];

  // Lowercase, single spaces, no periods; a trailing "st" becomes "state" so
  // Kalshi's "Missouri St." and Unabated's "Missouri State" agree. A leading
  // "St." (Saint) is left alone.
  function normalizeName(name) {
    if (typeof name !== "string") return null;
    const flat = name.toLowerCase().replace(/\./g, "").replace(/\s+/g, " ").trim();
    if (flat === "") return null;
    return flat.replace(/ st$/, " state");
  }

  function slugOf(normalized) {
    return normalized.replace(/ /g, "-");
  }

  // A pro league whose nicknames are unique: full name, "<CODE> <Nickname>",
  // the city when it names one team, and the nickname alone as the fallback.
  function buildProTable(league, teamRows, extraAliases) {
    const aliases = new Map();
    const nicknames = new Map();
    for (const [key, fullName, city, abbreviation] of teamRows) {
      const leagueKey = `${league}:${key}`;
      const nickname = fullName.split(" ").pop();
      aliases.set(normalizeName(fullName), leagueKey);
      aliases.set(normalizeName(`${abbreviation} ${nickname}`), leagueKey);
      if (city) aliases.set(normalizeName(city), leagueKey);
      nicknames.set(normalizeName(nickname), leagueKey);
    }
    for (const [alias, key] of extraAliases || []) aliases.set(normalizeName(alias), `${league}:${key}`);
    return { aliases, nicknames };
  }

  function buildCfbTable() {
    const aliases = new Map();
    for (const name of CFB_TEAMS) {
      const normalized = normalizeName(name);
      aliases.set(normalized, `cfb:${slugOf(normalized)}`);
    }
    return { aliases, nicknames: null };
  }

  const TABLES = { nfl: buildProTable("nfl", NFL_TEAMS), mlb: buildProTable("mlb", MLB_TEAMS, MLB_ALIASES), cfb: buildCfbTable() };

  // Key for (league, name) or null. The nickname fallback (last word) exists
  // only for leagues whose nicknames are unique.
  function teamKey(league, name) {
    const table = TABLES[league];
    const normalized = normalizeName(name);
    if (!table || normalized == null) return null;
    const exact = table.aliases.get(normalized);
    if (exact) return exact;
    if (!table.nicknames) return null;
    return table.nicknames.get(normalized.split(" ").pop()) ?? null;
  }

  function knownLeagues() {
    return Object.keys(TABLES);
  }

  const api = { teamKey, normalizeName, knownLeagues };

  if (typeof module !== "undefined" && module.exports) {
    module.exports = api;
  } else {
    root.UnabatedTeams = api;
  }
})(typeof globalThis !== "undefined" ? globalThis : this);
