// Presentation helpers for the bet-history flags in the Unabated Ticket panel
// (#114): source freshness, the header line under the tabs, the Ticket-tab
// banner truncation, the service-payload merge, and the ticket -> line shape.
// Pure: no DOM, no fetch, no chrome.* — loaded as a plain <script> in
// panel.html after bets.js (exposes globalThis.UnabatedBetsView) and via
// require() in tests/betsview.test.js. panel.js keeps only the DOM writes.
//
// Inputs
//   payload      the last /bets.json body the panel fetched:
//                {generatedAt, sources: {kalshi: {fetchedAt, ok, error, count}}, bets}
//   serviceState what panel.js remembers about the service itself:
//                {okAt, error, errorAt, unreachableSince} (ms epochs, error text)
//   records      normalised bet records (bets.js contract)
// Outputs plain objects / strings; nothing here writes anywhere.

(function (root) {
  "use strict";

  const bets = typeof module !== "undefined" && module.exports ? require("./bets.js") : root.UnabatedBets;

  // Every venue the plan registers a source for; #115/#116/#117 fill the
  // missing ones. Listing them keeps the Bets tab honest about coverage.
  const VENUES = ["kalshi", "betonline", "novig", "prophetx"];
  const FRESH_MS = 5 * 60 * 1000;
  const STALE_MS = 60 * 60 * 1000;
  const BANNER_MAX_LINES = 5;
  const DEFAULT_BETS_SETTINGS = { serviceUrl: "http://127.0.0.1:8094", hideBet: false };
  const BADGE_TEXT = { same_line: "BET", same_side: "BET", opposite: "OTHER SIDE", same_game: "GAME" };

  // "20 s" / "3 min" / "2 h" / "3 d" — the header line's short form.
  function fmtAgeShort(ms) {
    const age = Math.max(0, ms);
    if (age < 60 * 1000) return `${Math.round(age / 1000)} s`;
    if (age < 60 * 60 * 1000) return `${Math.round(age / 60000)} min`;
    if (age < 48 * 60 * 60 * 1000) return `${Math.round(age / 3600000)} h`;
    return `${Math.round(age / 86400000)} d`;
  }

  // green under 5 min, amber under 60 min, red past that or never fetched.
  function freshnessLevel(ageMs) {
    if (ageMs == null) return "red";
    if (ageMs < FRESH_MS) return "green";
    if (ageMs < STALE_MS) return "amber";
    return "red";
  }

  // One row per venue: what the service reported for it, or "no source configured".
  function sourceRows(payload, now) {
    const sources = payload && payload.sources && typeof payload.sources === "object" ? payload.sources : {};
    return VENUES.map((venue) => {
      const source = sources[venue];
      if (!source) {
        return { venue, configured: false, level: "none", ageMs: null, ageText: "—", fetchedAt: null, count: null, error: null, note: "no source configured" };
      }
      const fetchedMs = source.fetchedAt ? Date.parse(source.fetchedAt) : NaN;
      const ageMs = Number.isFinite(fetchedMs) ? now - fetchedMs : null;
      return {
        venue, configured: true, level: freshnessLevel(ageMs), ageMs,
        ageText: ageMs == null ? "never" : fmtAgeShort(ageMs),
        fetchedAt: Number.isFinite(fetchedMs) ? source.fetchedAt : null,
        count: typeof source.count === "number" ? source.count : 0,
        error: source.ok === false ? (source.error || "last poll failed") : null,
        note: null,
      };
    });
  }

  // The service itself: reachable, or unreachable since when and why.
  function serviceStatus(serviceState, now) {
    if (!serviceState || (serviceState.okAt == null && serviceState.errorAt == null)) {
      return { unreachable: true, text: "bets service not reached yet" };
    }
    if (serviceState.error) {
      const since = serviceState.unreachableSince ?? serviceState.errorAt;
      const sinceText = new Date(since).toLocaleTimeString([], { hour: "numeric", minute: "2-digit" });
      return { unreachable: true, text: `bets service unreachable since ${sinceText} (${fmtAgeShort(now - since)}): ${serviceState.error}` };
    }
    return { unreachable: false, text: `bets service reached ${fmtAgeShort(now - serviceState.okAt)} ago` };
  }

  // The Ticket tab's warning fires when nothing can vouch for the flags:
  // no source has ever reported, or every one that has is past the stale bound.
  function sourcesUnavailable(payload, now) {
    const configured = sourceRows(payload, now).filter((row) => row.configured);
    return configured.length === 0 || configured.every((row) => row.level === "red");
  }

  function openCount(records) {
    return records.filter((record) => record.status === "open").length;
  }

  // "bets: 14 open · kalshi 20 s · betonline — · novig — · prophetx —"
  function headerLine(records, payload, now) {
    const venues = sourceRows(payload, now).map((row) => `${row.venue} ${row.configured ? row.ageText : "—"}`);
    return [`bets: ${openCount(records)} open`, ...venues].join(" · ");
  }

  // At most `max` matches for the banner, strongest first (matchBets already
  // sorts), and how many were cut.
  function bannerLines(matches, max) {
    const limit = typeof max === "number" ? max : BANNER_MAX_LINES;
    return { shown: matches.slice(0, limit), more: Math.max(0, matches.length - limit) };
  }

  function badgeText(tier) {
    return BADGE_TEXT[tier] || null;
  }

  // Venues whose latest service poll succeeded: the payload is then the whole
  // store window for that venue (the service keeps records across its own
  // failed polls), so a stored record it no longer lists is gone for good — a
  // reset service DB, a purged fill — and must not flag lines forever.
  function venuesWithFreshPull(payload) {
    const sources = payload && payload.sources && typeof payload.sources === "object" ? payload.sources : {};
    return new Set(Object.keys(sources).filter((venue) => sources[venue] && sources[venue].ok === true));
  }

  // Records to keep after a poll: the stored ones and the fresh payload deduped
  // on native id (newest wins), minus stored records of a venue whose pull
  // succeeded without them; team keys filled, then the retention prune.
  function mergeServicePayload(storedRecords, payload, now) {
    const fresh = venuesWithFreshPull(payload);
    const listed = new Set((payload.bets || []).map((record) => record.id));
    const kept = (storedRecords || []).filter((record) => !fresh.has(record.venue) || listed.has(record.id));
    const merged = bets.dedupeByNativeId([kept, payload]);
    return bets.pruneForRetention(bets.resolveTeamKeys(merged), now);
  }

  // The captured ticket as the describeLine-shaped row the matcher reads.
  // Unabated's eventStart is naive UTC ("2026-09-12T23:30:00"); tickets
  // captured before page.js carried `period` are full game.
  function ticketAsLine(ticket) {
    const eventStart = ticket.eventStart == null ? null : String(ticket.eventStart);
    const eventStartMs = eventStart == null ? null : Date.parse(eventStart.endsWith("Z") ? eventStart : `${eventStart}Z`);
    return {
      league: ticket.league, eventId: ticket.eventId ?? null,
      awayTeam: ticket.awayTeam ?? null, homeTeam: ticket.homeTeam ?? null,
      eventStart, eventStartMs: Number.isFinite(eventStartMs) ? eventStartMs : null,
      betType: ticket.betType, period: ticket.period || "FG",
      sideIndex: ticket.sideIndex, points: ticket.points ?? null, rotation: ticket.rotation ?? null,
    };
  }

  function sanitizeBetsSettings(stored) {
    const base = { ...DEFAULT_BETS_SETTINGS };
    if (!stored || typeof stored !== "object") return base;
    if (typeof stored.serviceUrl === "string" && /^https?:\/\/\S+$/.test(stored.serviceUrl)) base.serviceUrl = stored.serviceUrl.replace(/\/+$/, "");
    if (typeof stored.hideBet === "boolean") base.hideBet = stored.hideBet;
    return base;
  }

  const api = {
    VENUES, FRESH_MS, STALE_MS, BANNER_MAX_LINES, DEFAULT_BETS_SETTINGS,
    fmtAgeShort, freshnessLevel, sourceRows, serviceStatus, sourcesUnavailable, openCount, headerLine,
    bannerLines, badgeText, venuesWithFreshPull, mergeServicePayload, ticketAsLine, sanitizeBetsSettings,
  };

  if (typeof module !== "undefined" && module.exports) {
    module.exports = api;
  } else {
    root.UnabatedBetsView = api;
  }
})(typeof globalThis !== "undefined" ? globalThis : this);
