// Unabated Ticket phone page — what the page shows, as plain data, from one
// GET /edges.json body (server/edges_payload.js documents its shape). Pure:
// no DOM, no fetch, the clock passed in, so the node tests drive it on a
// fixture and phone.js only turns the result into elements.
//
// Every figure comes from the payload as the runner computed it (the panel's
// own modules); this file only picks fields and words, it never re-derives an
// edge or a stake.
//
//   pageView(payload, nowMs, timeZone?) ->
//     {updated, chips: [{text, tone}], tailFlex, summary, cards: [card], footer}
//   card {key, header, when, matchup, badges: [{text, tone}], side, priceLine,
//         liquidity, edge, edgeTier, moveLabel, related: [text], stake,
//         stakeNote, atSize, othersText}
//   tone: "ok" | "warn" | "bad" | "plain"

(function (root) {
  "use strict";

  // The runner holds the last good read; past this the page says the data is old.
  const STALE_AFTER_MS = 2 * 60 * 1000;

  function fmtAmerican(price) {
    if (typeof price !== "number") return "—";
    return price > 0 ? `+${price}` : `−${Math.abs(price)}`;
  }

  function fmtPoints(points) {
    if (typeof points !== "number") return "";
    return points > 0 ? `+${points}` : `${points}`;
  }

  function fmtWholeDollars(value) {
    return `$${Math.round(value).toLocaleString("en-US")}`;
  }

  // "12s ago" / "4 min ago" / "2 h ago" for an epoch in ms; "never" for none.
  function fmtAgo(thenMs, nowMs) {
    if (typeof thenMs !== "number") return "never";
    const seconds = Math.max(0, Math.round((nowMs - thenMs) / 1000));
    if (seconds < 60) return `${seconds}s ago`;
    const minutes = Math.round(seconds / 60);
    if (minutes < 60) return `${minutes} min ago`;
    return `${Math.round(minutes / 60)} h ago`;
  }

  // "Sun Oct 11 · 6:30 AM" in the phone's own time zone (or the one passed).
  function fmtStart(startMs, timeZone) {
    if (typeof startMs !== "number") return "";
    const date = new Date(startMs);
    const day = date.toLocaleDateString("en-US", { weekday: "short", month: "short", day: "numeric", timeZone });
    const time = date.toLocaleTimeString("en-US", { hour: "numeric", minute: "2-digit", timeZone });
    return `${day.replace(",", "")} · ${time}`;
  }

  const PERIOD_WORDS = { FG: "Full game", "1H": "1st half", "2H": "2nd half", "1Q": "1st quarter" };

  function periodWords(period) {
    return PERIOD_WORDS[period] || period || "";
  }

  // The side as the panel names it ("Under 48.5", "Philadelphia Eagles +3").
  function sideText(card, row) {
    if (row.sideLabel) return row.sideLabel;
    const points = row.betType === "Moneyline" ? "" : ` ${fmtPoints(row.points)}`;
    return `${card.sideName}${points}`.trim();
  }

  function scannerChip(scanner, nowMs) {
    if (!scanner) return { text: "Feed: no data", tone: "bad" };
    const leagues = Array.isArray(scanner.leaguesLoaded) ? scanner.leaguesLoaded.length : 0;
    const errors = scanner.leagueErrors ? Object.keys(scanner.leagueErrors).length : 0;
    // The scanner sets `error` both when nothing loaded and when only some
    // leagues failed; only the first leaves the list empty.
    if (scanner.error && leagues === 0) return { text: `Feed error: ${scanner.error}`, tone: "bad" };
    if (scanner.error) return { text: `Feed ok · ${leagues} leagues · ${scanner.error}`, tone: "warn" };
    if (scanner.loading && scanner.loading.done < scanner.loading.total) {
      return { text: `Feed loading ${scanner.loading.done}/${scanner.loading.total}`, tone: "warn" };
    }
    if (typeof scanner.snapshotBuiltAt === "number" && nowMs - scanner.snapshotBuiltAt > 10 * 60 * 1000) {
      return { text: `Feed old · built ${fmtAgo(scanner.snapshotBuiltAt, nowMs)}`, tone: "warn" };
    }
    if (errors) return { text: `Feed ok · ${leagues} leagues · ${errors} failing`, tone: "warn" };
    return { text: `Feed ok · ${leagues} leagues`, tone: "ok" };
  }

  // The stakes are sized against open bets, so an unreadable bets service is
  // the loudest warning the page has: stakes then ignore what is already held.
  function betsChip(betsService, nowMs) {
    if (!betsService || betsService.okAt == null) {
      return { text: "Bets service not read — stakes ignore open bets", tone: "bad" };
    }
    if (betsService.error) {
      return { text: `Bets service down ${fmtAgo(betsService.unreachableSince ?? betsService.okAt, nowMs)} — stakes use bets as of then`, tone: "bad" };
    }
    const failing = Object.entries(betsService.sources || {})
      .filter(([, source]) => source && source.ok === false)
      .map(([name]) => name);
    if (failing.length) return { text: `Bets ok · ${failing.join(", ")} failing`, tone: "warn" };
    return { text: `Bets ok · ${betsService.openBets ?? 0} open`, tone: "ok" };
  }

  function settingsChip(settings) {
    const stake = settings && settings.stake ? settings.stake : null;
    const minEdge = settings && settings.edges ? settings.edges.minEdgePct : null;
    if (!stake) return { text: "Settings: none", tone: "warn" };
    const kellyText = stake.multiplier === 0.25 ? "¼ Kelly" : `${stake.multiplier}× Kelly`;
    const parts = [fmtWholeDollars(stake.bankroll), kellyText];
    if (typeof minEdge === "number") parts.push(`${minEdge}% min`);
    return { text: parts.join(" · "), tone: settings.error ? "warn" : "plain" };
  }

  function badgeView(badge) {
    const tone = badge.kind === "against" ? "bad" : "warn";
    return { text: badge.text, tone };
  }

  function cardView(card, timeZone) {
    const row = card.best;
    const others = Array.isArray(card.others) ? card.others.length : 0;
    return {
      key: card.key,
      header: [card.leagueLabel, card.betType, periodWords(card.period)].filter(Boolean).join(" · "),
      when: fmtStart(card.eventStartMs, timeZone),
      matchup: `${card.awayTeam} @ ${card.homeTeam}`,
      badges: (row.badges || []).map(badgeView),
      side: sideText(card, row),
      priceLine: `${row.book ? row.book.name : "?"} ${fmtAmerican(row.price)} · fair ${fmtAmerican(row.fair)}`,
      liquidity: typeof row.liquidity === "number" ? `liq to win ${fmtWholeDollars(row.liquidity)}` : null,
      edge: `${row.edgePct.toFixed(2)}%`,
      edgeTier: row.edgeTier || "thin",
      moveLabel: row.move && row.move.label ? row.move.label : null,
      related: (row.related || []).map((item) => `${item.text}${item.tag ? ` — ${item.tag}` : ""}`),
      stake: row.rail ? row.rail.text : "—",
      stakeNote: row.rail ? row.rail.note : null,
      atSize: Boolean(row.rail && row.rail.atSize),
      othersText: others === 0 ? "no other lines" : `${others} more line${others === 1 ? "" : "s"}`,
    };
  }

  // Grouped payloads list cards; an ungrouped one lists lines, each shown as a
  // card of one so the page has a single layout.
  function cardsOf(payload) {
    const items = Array.isArray(payload.items) ? payload.items : [];
    if (payload.grouped) return items;
    return items.map((row) => ({ ...row, best: row, others: [] }));
  }

  function pageView(payload, nowMs, timeZone) {
    const generatedMs = payload && payload.generatedAt ? Date.parse(payload.generatedAt) : null;
    const stale = typeof generatedMs === "number" && nowMs - generatedMs > STALE_AFTER_MS;
    const cards = cardsOf(payload || {}).map((card) => cardView(card, timeZone));
    const total = payload && typeof payload.total === "number" ? payload.total : cards.length;
    const shown = cards.length < total ? ` (top ${cards.length})` : "";
    return {
      updated: stale ? `data from ${fmtAgo(generatedMs, nowMs)}` : `updated ${fmtAgo(generatedMs, nowMs)}`,
      stale,
      chips: [scannerChip(payload.scanner, nowMs), betsChip(payload.betsService, nowMs), settingsChip(payload.settings)],
      tailFlex: payload.tailFlex || "",
      summary: `${total} ${payload.unit || "cards"}${shown} · pregame`,
      cards,
      empty: cards.length === 0 ? "No edges pass your filters right now." : null,
    };
  }

  const api = { STALE_AFTER_MS, pageView, cardView, fmtAmerican, fmtAgo, fmtStart };
  if (typeof module !== "undefined" && module.exports) {
    module.exports = api;
  } else {
    root.UnabatedPhoneView = api;
  }
})(typeof globalThis !== "undefined" ? globalThis : this);
