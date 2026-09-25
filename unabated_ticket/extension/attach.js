// The Bets tab's Attach picker: Cal attaches an open bet the matcher could not
// place to its game on the board, and the attach teaches the venue's team
// names. Pure: no DOM, no fetch, no chrome.* — loaded as a plain <script> in
// panel.html after teams.js and bets.js (exposes globalThis.UnabatedAttach)
// and via require() in tests/attach.test.js.
//
// Inputs
//   bet    a bets.js record with team keys resolved (resolveTeamKeys): id,
//          venue, league (may be null), awayTeam / homeTeam and their keys,
//          side, betType, period, points, rotation, eventStart / eventDate /
//          placedAt.
//   lines  the board as the panel's boardLines() gives it, one
//          feed.describeLine row per event: eventId, league, leagueLabel,
//          awayTeam / homeTeam, awayTeamId / homeTeamId, awayRotation /
//          homeRotation, eventStartMs.
// Outputs
//   attachCandidates(bet, lines, {query}) -> {scope, events: [{event, why}], more}
//          step 1: the games to pick from, best fit first.
//   attachPlan(bet, event, {swapped, lines}) -> {names, crosswalk, betOnGame}
//          step 2: each venue name lined up with an event team and tagged
//          known / learn / fix, the crosswalk rows the attach writes, and the
//          bet restated on the game.
//   pinRequest(bet, event, plan) -> the POST /pins.json body.
//   gameLabel(event) / gameMeta(event) -> "Away @ Home" / "Sat Sep 26 8:00 PM ET · rot 371 / 372".
// Nothing here writes anywhere; the panel POSTs pinRequest to the bets
// service, which stores the pin and the rows (bets.duckdb::bet_pins,
// team_crosswalk).

(function (root) {
  "use strict";

  const isNode = typeof module !== "undefined" && module.exports;
  const teams = isNode ? require("./teams.js") : root.UnabatedTeams;
  const bets = isNode ? require("./bets.js") : root.UnabatedBets;

  // More than a screen of games is noise: the search box narrows the rest.
  const MAX_CANDIDATES = 25;
  // A one-letter query matches half the board.
  const MIN_QUERY_LENGTH = 2;
  const STATUS_KNOWN = "known";
  const STATUS_LEARN = "learn";
  const STATUS_FIX = "fix";

  // ---- dates --------------------------------------------------------------

  // "Fri Sep 25" for a "YYYY-MM-DD" calendar date (read at UTC noon, so no
  // zone can move it to another day).
  function dateLabel(dateString) {
    const [year, month, day] = dateString.split("-").map(Number);
    const formatter = new Intl.DateTimeFormat("en-US", { timeZone: "UTC", weekday: "short", month: "short", day: "numeric" });
    const parts = {};
    for (const part of formatter.formatToParts(new Date(Date.UTC(year, month - 1, day, 12)))) parts[part.type] = part.value;
    return `${parts.weekday} ${parts.month} ${parts.day}`;
  }

  // ---- board games --------------------------------------------------------

  // One entry per board event of the league (every league when null), each
  // with its two team keys.
  function boardGames(lines, league) {
    const games = new Map();
    for (const line of lines || []) {
      if (line.eventId == null || games.has(line.eventId)) continue;
      if (league != null && line.league !== league) continue;
      games.set(line.eventId, {
        eventId: line.eventId, league: line.league, leagueLabel: line.leagueLabel || String(line.league).toUpperCase(),
        awayTeam: line.awayTeam ?? null, homeTeam: line.homeTeam ?? null,
        awayTeamId: line.awayTeamId ?? null, homeTeamId: line.homeTeamId ?? null,
        awayKey: line.awayTeamId != null ? teams.keyOf(line.league, line.awayTeamId) : null,
        homeKey: line.homeTeamId != null ? teams.keyOf(line.league, line.homeTeamId) : null,
        awayRotation: line.awayRotation ?? null, homeRotation: line.homeRotation ?? null,
        eventStartMs: typeof line.eventStartMs === "number" ? line.eventStartMs : null,
      });
    }
    return Array.from(games.values());
  }

  // Why a game is a good fit for the bet, or null: one of the bet's teams
  // already resolves to one of its teams, or the bet's rotation is one of its.
  function fitReason(bet, game) {
    for (const side of ["away", "home"]) {
      const key = side === "away" ? bet.awayKey : bet.homeKey;
      if (key && (key === game.awayKey || key === game.homeKey)) {
        return `${side === "away" ? bet.awayTeam : bet.homeTeam} matches`;
      }
    }
    if (bet.rotation != null && (bet.rotation === game.awayRotation || bet.rotation === game.homeRotation)) {
      return `rot ${bet.rotation} matches`;
    }
    return null;
  }

  function nameContains(game, query) {
    const wanted = teams.normalizeName(query);
    return [game.awayTeam, game.homeTeam].some((name) => typeof name === "string" && teams.normalizeName(name).includes(wanted));
  }

  function inWindow(game, window) {
    if (game.eventStartMs == null) return false;
    const date = bets.easternDateOf(game.eventStartMs);
    return date >= window.fromDate && date <= window.toDate;
  }

  // Best fit first, then the start nearest the bet's (soonest first when the
  // bet has no start).
  function compareCandidates(betStartMs) {
    return (a, b) => {
      if ((a.why == null) !== (b.why == null)) return a.why == null ? 1 : -1;
      const aStart = a.event.eventStartMs ?? Infinity;
      const bStart = b.event.eventStartMs ?? Infinity;
      if (Number.isFinite(betStartMs)) return Math.abs(aStart - betStartMs) - Math.abs(bStart - betStartMs);
      return aStart - bStart;
    };
  }

  // Step 1. With no query: the bet's league around its date. With a query of
  // two letters or more: every game of the league whose team name contains
  // it, at any date. Either way best fit first.
  function attachCandidates(bet, lines, options) {
    const query = options && typeof options.query === "string" ? options.query.trim() : "";
    const searching = query.length >= MIN_QUERY_LENGTH;
    const games = boardGames(lines, bet.league ?? null);
    const leagueText = bet.league ? String(bet.league).toUpperCase() : "Every league";
    // The same dates the red flag checks for a game (bets.betDateWindow).
    const window = searching ? null : bets.betDateWindow(bet);
    const listed = searching ? games.filter((game) => nameContains(game, query))
      : window ? games.filter((game) => inWindow(game, window)) : games;
    const scope = window ? `${leagueText} · ${dateLabel(window.fromDate)} to ${dateLabel(window.toDate)}` : `${leagueText} · every date`;
    const betStartMs = bet.eventStart ? Date.parse(bet.eventStart) : NaN;
    const ranked = listed.map((game) => ({ event: game, why: fitReason(bet, game) })).sort(compareCandidates(betStartMs));
    return { scope, events: ranked.slice(0, MAX_CANDIDATES), more: Math.max(0, ranked.length - MAX_CANDIDATES) };
  }

  // ---- step 2 -------------------------------------------------------------

  function otherSide(side) {
    return side === "away" ? "home" : "away";
  }

  // The board's name for a team key, or the key itself.
  function teamNameOfKey(key, lines) {
    for (const line of lines || []) {
      if (line.awayTeamId != null && teams.keyOf(line.league, line.awayTeamId) === key) return line.awayTeam;
      if (line.homeTeamId != null && teams.keyOf(line.league, line.homeTeamId) === key) return line.homeTeam;
    }
    return key;
  }

  function signedPoints(points) {
    return points > 0 ? `+${points}` : String(points);
  }

  // The event side the bet sits on: its venue side through the mapping, or
  // its rotation when the bet names no team on that side; null when neither says.
  function eventSideOfBet(bet, event, swapped) {
    if (bet.side !== "away" && bet.side !== "home") return null;
    if (bets.venueTeamOf(bet, bet.side)) return swapped ? otherSide(bet.side) : bet.side;
    if (bet.rotation != null && bet.rotation === event.awayRotation) return "away";
    if (bet.rotation != null && bet.rotation === event.homeRotation) return "home";
    return null;
  }

  // "Abilene Christian +7.5", "1H Under 30.5 -110": the bet as it reads on the picked game.
  function betOnGame(bet, event, swapped) {
    if (bet.betType !== "spread" && bet.betType !== "moneyline") return bets.describeBet(bet);
    const side = eventSideOfBet(bet, event, swapped);
    if (!side) return "side not known: it will not size a line on this game";
    const period = bet.period && bet.period !== "FG" ? `${bet.period} ` : "";
    const team = side === "away" ? event.awayTeam : event.homeTeam;
    if (bet.betType === "moneyline") return `${period}${team} to win`;
    return typeof bet.points === "number" ? `${period}${team} ${signedPoints(bet.points)}` : `${period}${team}`;
  }

  // Step 2. Each team the venue names, lined up with the event team on the
  // same side (the other side when swapped): `known` when it already resolves
  // to that team, `learn` when it resolves to nothing, `fix` when it resolves
  // to another team (`was` names it). The learn and fix rows are what the
  // attach writes to the crosswalk, in the event's league.
  function attachPlan(bet, event, options) {
    const swapped = Boolean(options && options.swapped);
    const lines = options && options.lines;
    const names = [];
    for (const venueSide of ["away", "home"]) {
      const venueTeam = bets.venueTeamOf(bet, venueSide);
      if (!venueTeam) continue;
      const eventSide = swapped ? otherSide(venueSide) : venueSide;
      const teamId = eventSide === "away" ? event.awayTeamId : event.homeTeamId;
      if (teamId == null) continue;
      const targetKey = teams.keyOf(event.league, teamId);
      // Keys resolved in another league (a bet with no league) say nothing about this one.
      const heldKey = bet.league === event.league ? (venueSide === "away" ? bet.awayKey : bet.homeKey) ?? null : null;
      const status = heldKey === targetKey ? STATUS_KNOWN : heldKey == null ? STATUS_LEARN : STATUS_FIX;
      names.push({
        venueSide, eventSide, venueTeamKey: venueTeam.key, venueTeamName: venueTeam.name,
        unabatedTeamId: String(teamId), unabatedTeamName: eventSide === "away" ? event.awayTeam : event.homeTeam,
        status, was: status === STATUS_FIX ? teamNameOfKey(heldKey, lines) : null,
      });
    }
    const crosswalk = names.filter((name) => name.status !== STATUS_KNOWN).map((name) => ({
      venue: bet.venue, league: event.league, venueTeamKey: name.venueTeamKey, venueTeamName: name.venueTeamName,
      unabatedTeamId: name.unabatedTeamId, unabatedTeamName: name.unabatedTeamName,
      learnedFrom: `attached by you: ${bet.id} on board event ${event.eventId}`,
    }));
    return { names, crosswalk, betOnGame: betOnGame(bet, event, swapped) };
  }

  // ---- labels -------------------------------------------------------------

  function gameLabel(event) {
    return `${event.awayTeam || "away team"} @ ${event.homeTeam || "home team"}`;
  }

  // "Sat Sep 26 8:00 PM ET"
  function gameStartLabel(ms) {
    if (typeof ms !== "number" || !Number.isFinite(ms)) return "start unknown";
    const formatter = new Intl.DateTimeFormat("en-US", {
      timeZone: "America/New_York", weekday: "short", month: "short", day: "numeric", hour: "numeric", minute: "2-digit", hour12: true,
    });
    const parts = {};
    for (const part of formatter.formatToParts(new Date(ms))) parts[part.type] = part.value;
    return `${parts.weekday} ${parts.month} ${parts.day} ${parts.hour}:${parts.minute} ${parts.dayPeriod} ET`;
  }

  // "Sat Sep 26 8:00 PM ET · rot 371 / 372": what tells two games of one pair apart.
  function gameMeta(event) {
    const rotations = event.awayRotation != null && event.homeRotation != null ? `rot ${event.awayRotation} / ${event.homeRotation}` : null;
    return [gameStartLabel(event.eventStartMs), rotations].filter(Boolean).join(" · ");
  }

  // The POST /pins.json body for an attach.
  function pinRequest(bet, event, plan) {
    return {
      pin: {
        betId: bet.id, venue: bet.venue, league: event.league, eventId: String(event.eventId),
        eventStart: event.eventStartMs == null ? null : new Date(event.eventStartMs).toISOString(),
        awayTeamId: event.awayTeamId == null ? null : String(event.awayTeamId),
        homeTeamId: event.homeTeamId == null ? null : String(event.homeTeamId),
        awayTeamName: event.awayTeam, homeTeamName: event.homeTeam,
      },
      crosswalk: plan.crosswalk,
    };
  }

  const api = {
    MAX_CANDIDATES, STATUS_KNOWN,
    attachCandidates, attachPlan, pinRequest, boardGames, gameLabel, gameMeta,
  };

  if (isNode) {
    module.exports = api;
  } else {
    root.UnabatedAttach = api;
  }
})(typeof globalThis !== "undefined" ? globalThis : this);
