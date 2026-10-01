// Unabated Ticket — live edges (pure; runs in the panel and under node).
//
// page.js posts what Unabated's live screen priced off its between-quarters
// in-game fair (`live_edges`, stored per league by content.js as
// "liveEdges:<league>"). This module turns those payloads into Edges-tab rows
// shaped like feed.describeLine's, applies the tab's filters, and says what
// the Live block's header shows. It prices nothing: edge and fair are the
// screen's own numbers.
//
// Side effects: none.

(function (root) {
  "use strict";

  const feed = typeof module !== "undefined" && module.exports ? require("./feed.js") : root.UnabatedFeed;

  // page.js posts at least every 5 s from a visible tab; past two of those
  // the tab has stopped reading and the numbers may be gone. A hidden tab's
  // timers run about once a minute (Chrome's intensive throttling), so it
  // gets a minute and a margin; its grid updates still post at once.
  const LIVE_STALE_MS = 10 * 1000;
  const LIVE_STALE_HIDDEN_MS = 75 * 1000;
  // A payload this old is from a tab closed long ago: show nothing for it.
  const LIVE_DROP_MS = 2 * 60 * 1000;
  const FAIR_READY = "ready";
  const STORAGE_PREFIX = "liveEdges:";

  // eventPeriodTypeId on a live row -> the part of the game that just ended.
  const PERIOD_ENDED = { 2: "1st half", 3: "2nd half", 4: "Q1", 5: "Q2", 6: "Q3", 7: "Q4" };

  function storageKeyOf(league) {
    return `${STORAGE_PREFIX}${league}`;
  }

  function leagueIdOfPath(path) {
    for (const [id, league] of Object.entries(feed.LEAGUES)) {
      if (league.path === path) return Number(id);
    }
    return null;
  }

  function isFresh(payload, now) {
    return Boolean(payload) && typeof payload.at === "number" && now - payload.at <= LIVE_DROP_MS;
  }

  // "end of Q1", "halftime"; the checkpoint type when the period is unknown.
  function checkpointLabel(game) {
    const ended = PERIOD_ENDED[game.eventPeriodTypeId];
    if (ended === "Q2" || ended === "1st half") return "halftime";
    if (ended) return `end of ${ended}`;
    if (typeof game.checkpointType === "string" && game.checkpointType) {
      return game.checkpointType.replace(/([a-z])([A-Z])/g, "$1 $2").toLowerCase();
    }
    return "break";
  }

  function matchupOf(game) {
    if (game.awayTeam && game.homeTeam) return `${game.awayTeam} @ ${game.homeTeam}`;
    return game.eventName || `event ${game.eventId}`;
  }

  function parseUtcMs(value) {
    if (typeof value !== "string" || !value) return null;
    // Unabated's timestamps are naive UTC ("2026-09-28T01:02:03.123"): read as UTC, never local.
    const ms = Date.parse(/[zZ]|[+-]\d\d:?\d\d$/.test(value) ? value : `${value}Z`);
    return Number.isFinite(ms) ? ms : null;
  }

  // What the Live block shows, from the stored payloads (any leagues).
  //   {games: [{...game, matchup, checkpoint, fairReady}], stale, readAgoMs, error}
  // Empty `games` = no live game on any screen: the block hides.
  function liveView(payloads, now) {
    const fresh = (payloads || []).filter((payload) => isFresh(payload, now));
    const games = [];
    for (const payload of fresh) {
      for (const game of payload.games || []) {
        games.push({ ...game, league: payload.league, matchup: matchupOf(game), checkpoint: checkpointLabel(game), fairReady: game.fairStatus === FAIR_READY });
      }
    }
    const newest = fresh.reduce((best, payload) => (!best || payload.at > best.at ? payload : best), null);
    const readAgoMs = newest ? now - newest.at : null;
    const staleAfterMs = newest && newest.visible === false ? LIVE_STALE_HIDDEN_MS : LIVE_STALE_MS;
    const errors = fresh.map((payload) => payload.error).filter(Boolean);
    return { games, stale: readAgoMs != null && readAgoMs > staleAfterMs, readAgoMs, error: errors.length ? errors.join("; ") : null };
  }

  function passesFilters(row, opts) {
    if (opts.minEdgePct != null && row.edgePct < opts.minEdgePct) return false;
    if (opts.periods && !opts.periods.has(row.periodTypeId)) return false;
    if (opts.betTypes && !opts.betTypes.has(row.betTypeId)) return false;
    if (opts.bookIds && !opts.bookIds.has(row.bookId)) return false;
    if (row.isAlt && !(opts.includeAlts && feed.altWithinDepthCap(row.fair, row.price))) return false;
    // Not a valid American price: nothing downstream (Kelly, cents) can use it.
    if (Math.abs(row.price) < 100) return false;
    if (opts.minLiquidityToWin != null && opts.minLiquidityToWin > 0 && row.liquidity != null) {
      const winPerDollar = row.price > 0 ? row.price / 100 : 100 / -row.price;
      if (row.liquidity * winPerDollar < opts.minLiquidityToWin) return false;
    }
    return true;
  }

  // One page.js row as an Edges row (feed.describeLine's shape plus `live`).
  function describeLiveRow(raw, game, payload) {
    const leagueId = leagueIdOfPath(payload.league);
    const league = feed.LEAGUES[leagueId] || { label: payload.league, path: payload.league };
    const eventStartMs = parseUtcMs(game.eventStart);
    return {
      key: raw.key,
      leagueId,
      league: league.path,
      leagueLabel: league.label,
      eventId: raw.eventId,
      eventName: raw.eventName ?? game.eventName ?? null,
      eventStart: eventStartMs == null ? null : new Date(eventStartMs).toISOString(),
      eventStartMs,
      betTypeId: raw.betTypeId,
      betType: feed.BET_TYPES[raw.betTypeId],
      periodTypeId: raw.periodTypeId,
      period: feed.PERIODS[raw.periodTypeId] || `pt${raw.periodTypeId}`,
      sideIndex: raw.sideIndex,
      sideKey: raw.sideKey,
      sideLabel: raw.sideLabel || `side ${raw.sideIndex}`,
      awayTeam: game.awayTeam ?? null,
      homeTeam: game.homeTeam ?? null,
      awayTeamId: game.awayTeamId ?? null,
      homeTeamId: game.homeTeamId ?? null,
      homeAway: raw.sideIndex === 0 ? "Away" : "Home",
      rotation: null,
      awayRotation: null,
      homeRotation: null,
      points: raw.points,
      book: { id: raw.bookId, name: raw.bookName || `book ${raw.bookId}`, hasLiquidity: raw.liquidity != null },
      bookId: raw.bookId,
      price: raw.price,
      sourceFormat: raw.sourceFormat,
      sourcePrice: raw.sourcePrice,
      fair: raw.fair,
      edgePct: Math.round(raw.edgePct * 100) / 100,
      liquidity: raw.liquidity,
      marketId: raw.marketId,
      isBlurred: false,
      modifiedMs: parseUtcMs(raw.modifiedOn),
      openerPrice: null,
      openerPoints: null,
      isAlt: raw.isAlt === true,
      mainPoints: raw.mainPoints ?? null,
      venueIds: null,
      live: {
        checkpoint: checkpointLabel(game),
        producedMs: parseUtcMs(game.producedUtc),
        producedUtc: game.producedUtc ?? null,
        matchup: matchupOf(game),
      },
    };
  }

  // The listable live rows across every stored payload whose fair is ready,
  // best edge first. Stale payloads still list (greyed by the panel); a
  // payload past LIVE_DROP_MS lists nothing.
  //   opts  {minEdgePct, periods, betTypes, bookIds, includeAlts,
  //          minLiquidityToWin, now}
  function selectLiveEdges(payloads, opts) {
    const options = opts || {};
    const now = typeof options.now === "number" ? options.now : Date.now();
    const rows = [];
    for (const payload of payloads || []) {
      if (!isFresh(payload, now)) continue;
      const games = new Map((payload.games || []).map((game) => [game.eventId, game]));
      for (const raw of payload.rows || []) {
        const game = games.get(raw.eventId);
        if (!game || game.fairStatus !== FAIR_READY) continue;
        const row = describeLiveRow(raw, game, payload);
        if (!passesFilters({ ...row, bookId: raw.bookId }, options)) continue;
        rows.push(row);
      }
    }
    rows.sort((a, b) => b.edgePct - a.edgePct);
    return rows;
  }

  // One alert per line per break: the fair's production time is in the key,
  // so the next break alerts again on the same line.
  function liveAlertKey(row) {
    return `${row.key}@${row.live.producedUtc || "?"}`;
  }

  const api = {
    LIVE_STALE_MS, LIVE_STALE_HIDDEN_MS, LIVE_DROP_MS, STORAGE_PREFIX,
    storageKeyOf, checkpointLabel, liveView, selectLiveEdges, liveAlertKey,
  };

  if (typeof module !== "undefined" && module.exports) {
    module.exports = api;
  } else {
    root.UnabatedLive = api;
  }
})(typeof globalThis !== "undefined" ? globalThis : this);
