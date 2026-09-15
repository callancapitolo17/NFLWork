// Team-name -> stable key resolution for matching a bet from one venue against
// a line on Unabated. Pure: no DOM, no fetch, no chrome.* — loaded as a plain
// <script> in panel.html (exposes globalThis.UnabatedTeams, before bets.js)
// and via require() in tests/teams.test.js.
//
// The index is NOT hand-written (user decision 2026-09-11): every league
// snapshot the scanner parses carries Unabated's own team list (id, name,
// abbreviation), and the panel registers it here (registerTeams) and
// persists it (exportIndex / loadIndex, chrome.storage.local `teamsIndex`)
// so it is there before the first snapshot of the next session. A key is
// "<league>:<Unabated team id>" — the same id the board's own lines carry, so
// the board side never needs a name match at all (feed.describeLine's
// awayTeamId / homeTeamId).
//
// Input:  a league path as feed.LEAGUES uses it ("nfl", "cfb", "mlb", ...) and
//         a team name as some venue wrote it.
// Output: the key, or null when the name does not resolve. Null is the ONLY
//         answer for an unknown or ambiguous name — a wrong key could match a
//         bet to the wrong game, and a null shows up in the panel's unmatched
//         list with the raw name. Resolution, in order:
//   1. exact normalised name ("Missouri St." and "Missouri State" agree);
//   2. a hand alias for a venue spelling that cannot be derived (ALIASES);
//   3. the name minus a leading code token ("PIT Steelers" -> "Steelers",
//      "OAK Athletics" -> "Athletics", Kalshi's event-title style), again by 1-2;
//   4. a UNIQUE word-boundary containment against the league's names:
//      "Steelers" ends "Pittsburgh Steelers", "New England" starts "New
//      England Patriots", "Middle Tennessee" starts "Middle Tennessee State";
//      the other direction ("Grambling St." starts with "Grambling") only
//      when what is left is an institutional suffix (State, University,
//      College) — "Southern Mississippi" must not become "Southern", the
//      university. Two candidates ("Los Angeles", "Miami") is null. A false
//      positive here still cannot match the wrong game: bets.js also needs
//      the opponent and the start time to agree.

(function (root) {
  "use strict";

  // Venue spellings that no rule derives: [league, venue spelling, Unabated name].
  const ALIASES = [
    ["cfb", "UAlbany", "Albany"],                       // Novig
    ["cfb", "North Carolina State", "NC State"],        // Novig
    ["cfb", "Southern Mississippi", "Southern Miss"],   // Novig
    // Measured against the live CFB snapshot 2026-09-12: these three blocked
    // 8 of 111 open bets. None can be a rule. "A&M" must never join
    // INSTITUTION_SUFFIXES — Texas A&M is not Texas, and Unabated lists only
    // the teams playing this week, so the school a query extends may simply
    // be absent. "Louisiana" is genuinely ambiguous by containment (Louisiana
    // Tech and SE Louisiana both hit) and resolves to null without this line.
    ["cfb", "Prairie View A&M", "Prairie View"],        // Novig
    ["cfb", "Louisiana", "UL Lafayette"],               // Novig (the Ragin' Cajuns)
    ["cfb", "Southeastern Louisiana", "SE Louisiana"],  // Novig
    ["mlb", "Los Angeles D", "Los Angeles Dodgers"],    // Kalshi KXMLBRFI event title truncation
  ];
  const CODE_TOKEN_RE = /^[A-Z][A-Z0-9&]{1,4}$/;
  // What a venue may append to a team's name without naming a different team.
  const INSTITUTION_SUFFIXES = new Set(["state", "university", "college"]);

  // league -> { byName: Map(normalised -> key), names: [[normalised, key]] }
  const index = new Map();
  const aliasByLeague = new Map();
  for (const [league, spelling, target] of ALIASES) {
    if (!aliasByLeague.has(league)) aliasByLeague.set(league, new Map());
    aliasByLeague.get(league).set(normalizeName(spelling), target);
  }

  // Lowercase, single spaces, no periods; a trailing "st" becomes "state" so
  // Kalshi's "Missouri St." and Unabated's "Missouri State" agree. A leading
  // "St." (Saint) is left alone.
  function normalizeName(name) {
    if (typeof name !== "string") return null;
    const flat = name.toLowerCase().replace(/\./g, "").replace(/\s+/g, " ").trim();
    if (flat === "") return null;
    return flat.replace(/ st$/, " state");
  }

  function keyOf(league, teamId) {
    return `${league}:${teamId}`;
  }

  // Add (or refresh) a league's teams: [{id, name, abbreviation}].
  function registerTeams(league, teams) {
    if (!league || !Array.isArray(teams)) return;
    if (!index.has(league)) index.set(league, { byName: new Map(), byId: new Map() });
    const table = index.get(league);
    for (const team of teams) {
      if (!team || team.id == null || typeof team.name !== "string") continue;
      const normalized = normalizeName(team.name);
      if (!normalized) continue;
      table.byId.set(String(team.id), { id: team.id, name: team.name, abbreviation: team.abbreviation ?? null });
      table.byName.set(normalized, keyOf(league, team.id));
    }
  }

  // {league: [{id, name, abbreviation}]} — what the panel persists.
  function exportIndex() {
    const out = {};
    for (const [league, table] of index) out[league] = Array.from(table.byId.values());
    return out;
  }

  function loadIndex(stored) {
    if (!stored || typeof stored !== "object") return;
    for (const [league, teams] of Object.entries(stored)) registerTeams(league, teams);
  }

  function teamCount(league) {
    const table = index.get(league);
    return table ? table.byId.size : 0;
  }

  function exactKey(table, league, normalized) {
    const direct = table.byName.get(normalized);
    if (direct) return direct;
    const aliases = aliasByLeague.get(league);
    const target = aliases ? aliases.get(normalized) : null;
    return target ? table.byName.get(normalizeName(target)) ?? null : null;
  }

  // The one league name the query starts or ends on a word boundary, or that
  // the query extends by an institutional suffix; null when none or more than one.
  function containmentKey(table, normalized) {
    let found = null;
    for (const [candidate, key] of table.byName) {
      const extendsCandidate = normalized.startsWith(`${candidate} `) && INSTITUTION_SUFFIXES.has(normalized.slice(candidate.length + 1));
      const hit = candidate.startsWith(`${normalized} `) || candidate.endsWith(` ${normalized}`) || extendsCandidate;
      if (!hit) continue;
      if (found && found !== key) return null;
      found = key;
    }
    return found;
  }

  // Key for (league, name) or null.
  function teamKey(league, name) {
    const table = index.get(league);
    const normalized = normalizeName(name);
    if (!table || normalized == null) return null;
    const exact = exactKey(table, league, normalized);
    if (exact) return exact;
    const [first, ...rest] = String(name).trim().split(/\s+/);
    if (rest.length && CODE_TOKEN_RE.test(first)) {
      const remainder = normalizeName(rest.join(" "));
      const stripped = exactKey(table, league, remainder) || containmentKey(table, remainder);
      if (stripped) return stripped;
    }
    return containmentKey(table, normalized);
  }

  function knownLeagues() {
    return Array.from(index.keys());
  }

  const api = { teamKey, keyOf, normalizeName, registerTeams, exportIndex, loadIndex, teamCount, knownLeagues };

  if (typeof module !== "undefined" && module.exports) {
    module.exports = api;
  } else {
    root.UnabatedTeams = api;
  }
})(typeof globalThis !== "undefined" ? globalThis : this);
