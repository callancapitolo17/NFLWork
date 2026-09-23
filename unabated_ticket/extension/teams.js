// Team-name -> stable key resolution for matching a bet from one venue against
// a line on Unabated. Pure: no DOM, no fetch, no chrome.* — loaded as a plain
// <script> in panel.html (exposes globalThis.UnabatedTeams, before bets.js)
// and via require() in tests/teams.test.js.
//
// The index is NOT hand-written (user decision 2026-09-11): every league
// snapshot the scanner parses carries Unabated's own team list (id, name,
// abbreviation) and, on each game row's eventName, a SECOND spelling of both
// teams ("Prairie View A&M Panthers" where the list says "Prairie View" —
// #118, feed.teamSpellingsFromEventName); the panel registers both here
// (registerTeams) and persists them (exportIndex / loadIndex,
// chrome.storage.local `teamsIndex`) so they are there before the first
// snapshot of the next session. A key is "<league>:<Unabated team id>" — the
// same id the board's own lines carry, so the board side never needs a name
// match at all (feed.describeLine's awayTeamId / homeTeamId).
//
// Input:  a league path as feed.LEAGUES uses it ("nfl", "cfb", "mlb", ...) and
//         a team name as some venue wrote it.
// Output: the key, or null when the name does not resolve. Null is the ONLY
//         answer for an unknown or ambiguous name — a wrong key could match a
//         bet to the wrong game, and a null shows up in the panel's unmatched
//         list with the raw name. Resolution, in order:
//   1. exact normalised name ("Missouri St." and "Missouri State" agree),
//      against either registered spelling of a team;
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
  // The eventName spellings (#118) retired two college aliases on 2026-09-12
  // ("Prairie View A&M" and "Southeastern Louisiana" now start Unabated's own
  // long forms "Prairie View A&M Panthers" / "Southeastern Louisiana Lions");
  // teams.test.js keeps a case per retired alias. The four college rows that
  // stay were each re-measured against that day's CFB snapshot: "UAlbany" (the
  // long form is "Albany Great Danes"), "Southern Mississippi" ("Southern Miss
  // Golden Eagles" — nothing turns Mississippi into Miss), "North Carolina
  // State" (containment hits BOTH "North Carolina State Wolfpack" and, by the
  // institutional-suffix rule, "North Carolina"), and "Louisiana" (Louisiana
  // Tech, Louisiana Tech Bulldogs and SE Louisiana all hit). "A&M" must never
  // join INSTITUTION_SUFFIXES — Texas A&M is not Texas, and Unabated lists
  // only the teams playing this week, so the school a query extends may
  // simply be absent.
  const ALIASES = [
    ["cfb", "UAlbany", "Albany"],                       // Novig
    ["cfb", "North Carolina State", "NC State"],        // Novig
    ["cfb", "Southern Mississippi", "Southern Miss"],   // Novig
    ["cfb", "Louisiana", "UL Lafayette"],               // Novig (the Ragin' Cajuns)
    ["mlb", "Los Angeles D", "Los Angeles Dodgers"],    // Kalshi KXMLBRFI event title truncation
    ["mlb", "ARI DBACKS", "Arizona Diamondbacks"],      // Wagerzon (the only MLB spelling of its 30 that no rule keys, 2026-09-23)
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

  // Add (or refresh) a league's teams: [{id, name, abbreviation, eventName}].
  // Both spellings resolve to the same key. A refresh without an eventName
  // (the team's row has gone live, so this snapshot carries no eventName for
  // it) keeps the spelling already held, or the next session would lose it.
  function registerTeams(league, teams) {
    if (!league || !Array.isArray(teams)) return;
    if (!index.has(league)) index.set(league, { byName: new Map(), byId: new Map() });
    const table = index.get(league);
    for (const team of teams) {
      if (!team || team.id == null || typeof team.name !== "string") continue;
      const normalized = normalizeName(team.name);
      if (!normalized) continue;
      const held = table.byId.get(String(team.id));
      const eventName = typeof team.eventName === "string" && team.eventName !== "" ? team.eventName : held ? held.eventName : null;
      table.byId.set(String(team.id), { id: team.id, name: team.name, abbreviation: team.abbreviation ?? null, eventName });
      table.byName.set(normalized, keyOf(league, team.id));
      const normalizedEventName = normalizeName(eventName);
      if (!normalizedEventName) continue;
      // A list name always wins over a second spelling: if some other team's
      // name already reads the same, the eventName spelling is not registered
      // (a collision must not turn an exact hit into the wrong team).
      const heldKey = table.byName.get(normalizedEventName);
      if (heldKey && heldKey !== keyOf(league, team.id)) continue;
      table.byName.set(normalizedEventName, keyOf(league, team.id));
    }
  }

  // {league: [{id, name, abbreviation, eventName}]} — what the panel persists.
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

  // Distinct registered spellings across every league — grows when a team OR
  // a second spelling arrives, so the panel can tell either apart from a no-op refresh.
  function spellingCount() {
    let total = 0;
    for (const table of index.values()) total += table.byName.size;
    return total;
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

  const api = { teamKey, keyOf, normalizeName, registerTeams, exportIndex, loadIndex, teamCount, spellingCount, knownLeagues };

  if (typeof module !== "undefined" && module.exports) {
    module.exports = api;
  } else {
    root.UnabatedTeams = api;
  }
})(typeof globalThis !== "undefined" ? globalThis : this);
