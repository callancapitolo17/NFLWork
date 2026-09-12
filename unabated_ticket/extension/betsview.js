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
//   pageSources  venues read by a content script instead of the service, keyed
//                by venue: {novig: {bets, readAt, url, error, complete,
//                pageSeenAt}} — what novig_content.js writes to storage (#116)
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
  const DEFAULT_BETS_SETTINGS = { serviceUrl: "http://127.0.0.1:8094" };
  // How a page-sourced venue is refreshed, for the Bets tab when its read is
  // old or missing: the second form when its tab is open but has not shown
  // the screen the content script mirrors.
  const PAGE_SOURCE_HINT = {
    novig: { closed: "open app.novig.us and its Portfolio screen in a tab to refresh", open: "Novig tab is open — open its Portfolio screen to refresh" },
  };

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

  // A venue read by a content script: fresh as of its last portfolio
  // response; the hint says how to refresh once that is old or absent.
  function pageSourceRow(venue, source, now) {
    const readMs = source && source.readAt ? Date.parse(source.readAt) : NaN;
    const ageMs = Number.isFinite(readMs) ? now - readMs : null;
    const level = freshnessLevel(ageMs);
    const seenMs = source && source.pageSeenAt ? Date.parse(source.pageSeenAt) : NaN;
    const tabOpen = Number.isFinite(seenMs) && now - seenMs < FRESH_MS;
    const hint = PAGE_SOURCE_HINT[venue] || { closed: "open the venue's site in a tab to refresh", open: "the venue's tab is open — open its bets screen to refresh" };
    const note = level === "red" ? (tabOpen ? hint.open : hint.closed) : null;
    return {
      venue, configured: true, level, ageMs,
      ageText: ageMs == null ? "never" : fmtAgeShort(ageMs),
      fetchedAt: Number.isFinite(readMs) ? source.readAt : null,
      count: source && Array.isArray(source.bets) ? source.bets.length : 0,
      error: source && source.error ? source.error : null,
      note,
    };
  }

  // One row per venue: what the service reported for it, what a content
  // script wrote for it, or "no source configured".
  function sourceRows(payload, now, pageSources) {
    const sources = payload && payload.sources && typeof payload.sources === "object" ? payload.sources : {};
    const pages = pageSources && typeof pageSources === "object" ? pageSources : {};
    return VENUES.map((venue) => {
      const source = sources[venue];
      if (!source && pages[venue]) return pageSourceRow(venue, pages[venue], now);
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
  function sourcesUnavailable(payload, now, pageSources) {
    const configured = sourceRows(payload, now, pageSources).filter((row) => row.configured);
    return configured.length === 0 || configured.every((row) => row.level === "red");
  }

  function openCount(records) {
    return records.filter((record) => record.status === "open").length;
  }

  // "bets: 14 open · kalshi 20 s · betonline — · novig — · prophetx —"
  function headerLine(records, payload, now, pageSources) {
    const venues = sourceRows(payload, now, pageSources).map((row) => `${row.venue} ${row.configured ? row.ageText : "—"}`);
    return [`bets: ${openCount(records)} open`, ...venues].join(" · ");
  }

  // At most `max` matches for the banner, strongest first (matchBets already
  // sorts), and how many were cut.
  function bannerLines(matches, max) {
    const limit = typeof max === "number" ? max : BANNER_MAX_LINES;
    return { shown: matches.slice(0, limit), more: Math.max(0, matches.length - limit) };
  }

  // The badge's kind: held (warning tint), against (red), game (outline).
  // Dollars decide when a venue supplied them; a bet without a stake still
  // gets its kind from the tier, so an other-side position is never unflagged.
  function badgeKind(flag) {
    if (!flag || !flag.tier) return null;
    const exposure = flag.exposure || { held: 0, against: 0 };
    if (exposure.held > 0) return "held";
    if (exposure.against > 0) return "against";
    if (flag.tier === "same_line" || flag.tier === "same_side") return "held";
    if (flag.tier === "opposite") return "against";
    return "game";
  }

  // The row badge: what you already have on this market, in dollars when the
  // venue gave a stake ("held $300", "against $200"), bare otherwise. `held`
  // wins over `against` when both exist (the against bets stay in the tooltip);
  // a same_game match is a plain "game" marker — it does not change the size.
  function badgeText(flag) {
    const kind = badgeKind(flag);
    if (!kind) return null;
    const exposure = flag.exposure || { held: 0, against: 0 };
    const dollars = kind === "held" ? exposure.held : kind === "against" ? exposure.against : 0;
    return dollars > 0 ? `${kind} ${bets.formatStake(dollars)}` : kind;
  }

  // What to do with a Kelly stake given what you already hold on the market.
  //   none     nothing held either way: the stake stands
  //   add      held less than the stake: top up by `add` (the number to act on)
  //   at_size  held the stake or more (or no stake could be computed): nothing to add
  //   reverse  on the other side only: the stake stands and `net` is what is
  //            left after it cancels the against position (negative = still net against)
  function stakeAdvice(stake, exposure) {
    const held = exposure && exposure.held > 0 ? exposure.held : 0;
    const against = exposure && exposure.against > 0 ? exposure.against : 0;
    const sized = typeof stake === "number" && stake > 0;
    if (held > 0) {
      if (sized && stake > held) return { kind: "add", add: Math.round((stake - held) * 100) / 100, held, stake };
      return { kind: "at_size", held, stake: sized ? stake : null };
    }
    if (against > 0) return { kind: "reverse", against, stake: sized ? stake : null, net: sized ? Math.round((stake - against) * 100) / 100 : null };
    return { kind: "none" };
  }

  // "you hold Texas A&M -38.5 -110 · Kalshi · Sep 10 2:15 PM" — one line per
  // held or against bet, for the row's third line and the ticket's facts.
  function positionLines(flag) {
    const exposure = flag && flag.exposure;
    if (!exposure) return [];
    const describe = (bet) => `${bets.describeBet(bet)} · ${bets.formatStake(bet.stake)} · ${bet.venue ? bet.venue.charAt(0).toUpperCase() + bet.venue.slice(1) : "unknown venue"}${bet.placedAt ? ` · ${bets.formatPlacedAt(bet.placedAt)}` : ""}`;
    return [
      ...exposure.heldBets.map((bet) => `you hold ${describe(bet)}`),
      ...exposure.againstBets.map((bet) => `other side: ${describe(bet)}`),
    ];
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

  // Records to keep after a content-script venue wrote its read: the stored
  // ones and the read deduped on native id (newest wins). A COMPLETE read
  // (every list seen to its end) is authoritative for that venue, so a stored
  // record it no longer lists is dropped; an incomplete read only adds.
  function mergePageSource(storedRecords, venue, pageSource, now) {
    const read = pageSource && Array.isArray(pageSource.bets) ? pageSource.bets : [];
    const listed = new Set(read.map((record) => record.id));
    const authoritative = !!(pageSource && pageSource.complete === true);
    const kept = (storedRecords || []).filter((record) => record.venue !== venue || !authoritative || listed.has(record.id));
    const merged = bets.dedupeByNativeId([kept, read]);
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
    return base;
  }

  const api = {
    VENUES, FRESH_MS, STALE_MS, BANNER_MAX_LINES, DEFAULT_BETS_SETTINGS, PAGE_SOURCE_HINT,
    fmtAgeShort, freshnessLevel, sourceRows, serviceStatus, sourcesUnavailable, openCount, headerLine,
    bannerLines, badgeText, badgeKind, stakeAdvice, positionLines, venuesWithFreshPull, mergeServicePayload, mergePageSource, ticketAsLine, sanitizeBetsSettings,
  };

  if (typeof module !== "undefined" && module.exports) {
    module.exports = api;
  } else {
    root.UnabatedBetsView = api;
  }
})(typeof globalThis !== "undefined" ? globalThis : this);
