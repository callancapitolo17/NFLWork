// Bet Tracker — the numbers, without the DOM. Pure: no fetch, no clock, no
// storage. Turns the bets service's /bets.json body into tickets (one row per
// straight bet, one per parlay or teaser) and summarises them for the
// Overview and Analysis pages (server/tracker/tracker.js renders them).
//
// Loaded as a plain <script> by server/tracker/index.html (exposes
// globalThis.UnabatedTrackerStats) and via require() in
// tests/trackerstats.test.js.
//
// Inputs
//   records    /bets.json `bets` (the normalised record, docs/2026-09-11-issue-114-bet-history-plan.md)
//   fillFairs  /bets.json `fillFairs` ({betId, fairAmerican, ...}): Unabated's
//              fair for the bet's side when it was placed
//   removedBets /bets.json `exclusions` ({betId}): bets Cal removed on the Bets
//              page; a ticket with any removed record is marked `excluded`
// Conventions
//   P&L lands on the Pacific calendar day the bet SETTLED (closedAt).
//   A bet counts toward P&L only when it is won, lost, push or void; an open
//   bet is exposure, and a position sold or cashed out before settlement, or
//   a Bet105 bet that left the open list with no graded wager to settle it
//   ("closed"), or a bet whose result the venue no longer shows ("unknown")
//   has no known P&L and is counted as excluded.
//   Expected P&L, edge and calibration use only bets with a saved fair.
//   A record carrying the venue's own `pnl` (Kalshi: priced off the market's
//   net position, after fees) counts that number, and a "closed" one with it
//   (sold before settlement) counts too: a gain as a win, a loss as a loss.
//   A record with `mergedInto` is the other side of such a position and is
//   skipped; its trades are already in the record it names.

(function (root) {
  "use strict";

  const PACIFIC_TZ = "America/Los_Angeles";
  const DAY_MS = 24 * 60 * 60 * 1000;
  const HOUR_MS = 60 * 60 * 1000;
  const Z_95 = 1.96;
  const SETTLED_STATUSES = new Set(["won", "lost", "push", "void"]);
  const WEEKDAYS = ["Sun", "Mon", "Tue", "Wed", "Thu", "Fri", "Sat"];

  const VENUE_NAMES = {
    kalshi: "Kalshi", novig: "Novig", betonline: "BetOnline", bfa: "BFA",
    wagerzon: "Wagerzon", polymarket_us: "Polymarket US", bet105: "Bet105",
  };
  const BET_TYPE_NAMES = { moneyline: "Moneyline", spread: "Spread", total: "Total", other: "Other" };
  const KIND_NAMES = { straight: "Straight", parlay: "Parlay", teaser: "Teaser", kalshiCombo: "Kalshi combo" };
  // Kalshi's multivariate combo series: the MLB bots' RFQ fills (kalshi_mlb_mm,
  // kalshi_mlb_rfq) on the same account land in the bets service as these.
  const KALSHI_COMBO_SERIES = "KXMVECROSSCATEGORY";

  const ODDS_BUCKETS = ["−200 or shorter", "−199 to −121", "−120 to +120", "+121 to +199", "+200 and up"];
  const EDGE_BUCKETS = ["Under 0%", "0 to 2%", "2 to 4%", "4 to 6%", "6% and up"];
  const TIMING_BUCKETS = ["Under 1h before", "1 to 6h before", "6 to 24h before", "1 day+ before", "Live or unknown"];
  const STAKE_BUCKETS = ["Under $150", "$150 to $299", "$300 to $449", "$450 and up"];
  const NO_FAIR = "No fair saved";

  // Group-by dimensions of the Analysis page; `order` fixes the row order of
  // ordinal ones, the rest sort by handle.
  const GROUPS = [
    { key: "venue", label: "Venue" },
    { key: "league", label: "League" },
    { key: "market", label: "Market" },
    { key: "period", label: "Period" },
    { key: "kind", label: "Type", order: ["Straight", "Parlay", "Teaser", "Kalshi combo"] },
    { key: "oddsBucket", label: "Odds", order: ODDS_BUCKETS },
    { key: "edgeBucket", label: "Edge", order: EDGE_BUCKETS.concat(NO_FAIR) },
    { key: "timing", label: "Timing", order: TIMING_BUCKETS },
    { key: "weekday", label: "Weekday", order: WEEKDAYS },
    { key: "stakeBucket", label: "Stake", order: STAKE_BUCKETS },
  ];

  // ---- dates ----------------------------------------------------------------

  const pacificDateFormat = new Intl.DateTimeFormat("en-CA", {
    timeZone: PACIFIC_TZ, year: "numeric", month: "2-digit", day: "2-digit",
  });

  /** Epoch ms -> "YYYY-MM-DD" on the Pacific calendar. */
  function pacificDay(epochMs) {
    return pacificDateFormat.format(new Date(epochMs));
  }

  /** "YYYY-MM-DD" -> epoch ms of that date at 00:00 UTC (day arithmetic only, never a real instant). */
  function dayKeyToUtc(dayKey) {
    const [year, month, day] = dayKey.split("-").map(Number);
    return Date.UTC(year, month - 1, day);
  }

  function utcToDayKey(utcMs) {
    return new Date(utcMs).toISOString().slice(0, 10);
  }

  function addDays(dayKey, days) {
    return utcToDayKey(dayKeyToUtc(dayKey) + days * DAY_MS);
  }

  function weekdayOf(dayKey) {
    return WEEKDAYS[new Date(dayKeyToUtc(dayKey)).getUTCDay()];
  }

  function parseMs(iso) {
    if (!iso) return null;
    const ms = Date.parse(iso);
    return Number.isFinite(ms) ? ms : null;
  }

  // ---- prices ---------------------------------------------------------------

  function americanToDecimal(american) {
    if (!Number.isFinite(american) || Math.abs(american) < 100) return null;
    return american > 0 ? 1 + american / 100 : 1 + 100 / -american;
  }

  function decimalToAmerican(decimal) {
    if (!Number.isFinite(decimal) || decimal <= 1) return null;
    return decimal >= 2 ? Math.round((decimal - 1) * 100) : -Math.round(100 / (decimal - 1));
  }

  /** The decimal payout actually on the ticket: stake + toWin over stake, else the quoted price. */
  function payoutDecimal(stake, toWin, american) {
    if (stake > 0 && Number.isFinite(toWin) && toWin > 0) return 1 + toWin / stake;
    return americanToDecimal(american);
  }

  // ---- buckets --------------------------------------------------------------

  function oddsBucket(american) {
    if (!Number.isFinite(american)) return null;
    if (american <= -200) return ODDS_BUCKETS[0];
    if (american <= -121) return ODDS_BUCKETS[1];
    if (american <= 120) return ODDS_BUCKETS[2];
    if (american <= 199) return ODDS_BUCKETS[3];
    return ODDS_BUCKETS[4];
  }

  function edgeBucket(edge) {
    if (edge === null) return NO_FAIR;
    if (edge < 0) return EDGE_BUCKETS[0];
    if (edge < 0.02) return EDGE_BUCKETS[1];
    if (edge < 0.04) return EDGE_BUCKETS[2];
    if (edge < 0.06) return EDGE_BUCKETS[3];
    return EDGE_BUCKETS[4];
  }

  function timingBucket(placedMs, eventStartMs) {
    if (placedMs === null || eventStartMs === null || eventStartMs <= placedMs) return TIMING_BUCKETS[4];
    const hours = (eventStartMs - placedMs) / HOUR_MS;
    if (hours < 1) return TIMING_BUCKETS[0];
    if (hours < 6) return TIMING_BUCKETS[1];
    if (hours < 24) return TIMING_BUCKETS[2];
    return TIMING_BUCKETS[3];
  }

  function stakeBucket(stake) {
    if (stake < 150) return STAKE_BUCKETS[0];
    if (stake < 300) return STAKE_BUCKETS[1];
    if (stake < 450) return STAKE_BUCKETS[2];
    return STAKE_BUCKETS[3];
  }

  // ---- records -> tickets ---------------------------------------------------

  function venueName(venue) {
    return VENUE_NAMES[venue] || String(venue || "Unknown");
  }

  function leagueName(league) {
    return league ? String(league).toUpperCase() : "Unknown";
  }

  function signedNumber(points) {
    return (points > 0 ? "+" : points < 0 ? "−" : "") + Math.abs(points);
  }

  function eventLabel(record) {
    if (record.awayTeam && record.homeTeam) return record.awayTeam + " @ " + record.homeTeam;
    return "";
  }

  /** "BUF −6.5", "1H Over 47.5", "PHI ML"; the venue's own text when the record has no parsed side. */
  function selectionLabel(record) {
    const raw = record.raw || {};
    const fallback = String(raw.description || raw.marketTitle || raw.title || record.id || "");
    if (!record.side || record.betType === "other") return fallback;
    const prefix = record.period && record.period !== "FG" ? record.period + " " : "";
    if (record.betType === "total") {
      const side = record.side === "over" ? "Over" : "Under";
      return prefix + side + (Number.isFinite(record.points) ? " " + record.points : "");
    }
    const team = record.side === "away" ? record.awayTeam : record.homeTeam;
    if (!team) return fallback;
    if (record.betType === "moneyline") return prefix + team + " ML";
    return prefix + team + (Number.isFinite(record.points) ? " " + signedNumber(record.points) : "");
  }

  function pnlOf(status, stake, toWin, decimal) {
    if (status === "lost") return -stake;
    if (status === "push" || status === "void") return 0;
    if (status !== "won") return null;
    if (Number.isFinite(toWin)) return toWin;
    return decimal ? stake * (decimal - 1) : null;
  }

  function isTeaserTicket(legs) {
    return legs.some((leg) => /teaser/i.test(JSON.stringify(leg.raw || {})));
  }

  /** The ticket's fields that do not depend on whether it is a straight or a parlay. */
  function finishTicket(ticket, fairAmerican) {
    const decimal = payoutDecimal(ticket.stake, ticket.toWin, ticket.price);
    const fairDecimal = americanToDecimal(fairAmerican);
    const fairProb = fairDecimal ? 1 / fairDecimal : null;
    const edge = fairProb !== null && decimal ? fairProb * decimal - 1 : null;
    // Outcome variance of a win/lose bet: stake² · decimal² · p(1−p). With no
    // fair, the break-even probability stands in so the CI still has a width.
    const winProb = fairProb !== null ? fairProb : decimal ? 1 / decimal : 0.5;
    const placedMs = parseMs(ticket.placedAt);
    const closedMs = parseMs(ticket.closedAt);
    const hasVenuePnl = Number.isFinite(ticket.venuePnl)
      && (SETTLED_STATUSES.has(ticket.status) || ticket.status === "closed");
    const settled = hasVenuePnl || SETTLED_STATUSES.has(ticket.status);
    const settledDay = settled && closedMs !== null ? pacificDay(closedMs) : null;
    const pnl = !settled ? null : hasVenuePnl ? ticket.venuePnl : pnlOf(ticket.status, ticket.stake, ticket.toWin, decimal);
    return Object.assign(ticket, {
      decimal, fairAmerican: fairDecimal ? fairAmerican : null, fairProb, edge,
      expected: edge !== null ? ticket.stake * edge : null,
      variance: decimal ? ticket.stake * ticket.stake * decimal * decimal * winProb * (1 - winProb) : 0,
      placedMs, closedMs, settled, settledDay, pnl,
      displayPrice: Number.isFinite(ticket.price) ? ticket.price : decimalToAmerican(decimal),
      oddsBucket: oddsBucket(Number.isFinite(ticket.price) ? ticket.price : decimalToAmerican(decimal)),
      edgeBucket: edgeBucket(edge),
      timing: timingBucket(placedMs, parseMs(ticket.eventStart)),
      stakeBucket: stakeBucket(ticket.stake),
      weekday: settledDay ? weekdayOf(settledDay) : null,
    });
  }

  function isKalshiCombo(record) {
    return record.venue === "kalshi" && (record.raw || {}).series === KALSHI_COMBO_SERIES;
  }

  function straightTicket(record, fairAmerican) {
    const combo = isKalshiCombo(record);
    return finishTicket({
      id: record.id, venueKey: record.venue, venue: venueName(record.venue), league: leagueName(record.league),
      kind: combo ? KIND_NAMES.kalshiCombo : KIND_NAMES.straight,
      market: combo ? KIND_NAMES.kalshiCombo : BET_TYPE_NAMES[record.betType] || "Other",
      period: record.period || "FG", event: eventLabel(record), selection: selectionLabel(record),
      price: record.price, stake: Number(record.stake) || 0, toWin: record.toWin,
      venuePnl: Number.isFinite(record.pnl) ? record.pnl : null,
      status: record.status, placedAt: record.placedAt, closedAt: record.closedAt, eventStart: record.eventStart || null,
      legCount: 1, betIds: [record.id],
    }, fairAmerican);
  }

  /** When a multi-leg ticket was decided: a lost one at its earliest losing leg's
   * close (the venue's per-leg result, raw.legResult "LOSE"), any other at its
   * latest leg's close; the first leg's closedAt when no leg says. */
  function multiLegClosedAt(legs, status) {
    const closeMs = (leg) => parseMs(leg.closedAt);
    const losing = status === "lost"
      ? legs.filter((leg) => String((leg.raw || {}).legResult || "").toUpperCase() === "LOSE")
      : [];
    const decidingCloses = (losing.length ? losing : legs).map(closeMs).filter((ms) => ms !== null);
    if (!decidingCloses.length) return legs[0].closedAt;
    const decidedMs = losing.length ? Math.min(...decidingCloses) : Math.max(...decidingCloses);
    return new Date(decidedMs).toISOString();
  }

  /** One ticket from a parlay's or teaser's legs: every leg carries the ticket's stake, toWin and status. */
  function multiLegTicket(parlayId, legs) {
    const sorted = legs.slice().sort((a, b) => (a.legIndex || 0) - (b.legIndex || 0));
    const first = sorted[0];
    const kind = isTeaserTicket(sorted) ? KIND_NAMES.teaser : KIND_NAMES.parlay;
    const leagues = new Set(sorted.map((leg) => leagueName(leg.league)));
    const starts = sorted.map((leg) => parseMs(leg.eventStart)).filter((ms) => ms !== null);
    const closedAt = multiLegClosedAt(sorted, first.status);
    const legCount = first.legCount || sorted.length;
    const raw = first.raw || {};
    return finishTicket({
      id: parlayId, venueKey: first.venue, venue: venueName(first.venue),
      league: leagues.size === 1 ? leagues.values().next().value : "Mixed",
      kind, market: kind, period: "FG", event: legCount + "-leg " + kind.toLowerCase(),
      selection: sorted.map(selectionLabel).filter(Boolean).join(" · "),
      price: Number.isFinite(raw.parlayPrice) ? raw.parlayPrice : null,
      stake: Number(first.stake) || 0, toWin: first.toWin,
      status: first.status, placedAt: first.placedAt,
      closedAt,
      eventStart: starts.length ? new Date(Math.min(...starts)).toISOString() : null,
      legCount, betIds: sorted.map((leg) => leg.id),
    }, null);
  }

  /** /bets.json records -> tickets, newest placement first. */
  function buildTickets(records, fillFairs, removedBets) {
    const fairByBet = new Map((fillFairs || []).map((row) => [row.betId, row.fairAmerican]));
    const removedIds = new Set((removedBets || []).map((row) => row.betId));
    const legsByParlay = new Map();
    const tickets = [];
    for (const record of records || []) {
      if (record.mergedInto) continue;
      if (record.isParlayLeg && record.parlayId) {
        if (!legsByParlay.has(record.parlayId)) legsByParlay.set(record.parlayId, []);
        legsByParlay.get(record.parlayId).push(record);
        continue;
      }
      tickets.push(straightTicket(record, fairByBet.has(record.id) ? fairByBet.get(record.id) : null));
    }
    for (const [parlayId, legs] of legsByParlay) tickets.push(multiLegTicket(parlayId, legs));
    for (const ticket of tickets) ticket.excluded = ticket.betIds.some((id) => removedIds.has(id));
    return tickets.sort((a, b) => (b.closedMs || b.placedMs || 0) - (a.closedMs || a.placedMs || 0));
  }

  // ---- summaries ------------------------------------------------------------

  /** Totals of settled tickets with a known P&L. Expected/z cover only those with a saved fair. */
  function summarize(tickets) {
    const total = {
      bets: 0, wins: 0, losses: 0, pushes: 0, handle: 0, pnl: 0, variance: 0,
      withFair: 0, fairHandle: 0, fairPnl: 0, expected: 0, fairVariance: 0,
    };
    for (const ticket of tickets) {
      if (ticket.pnl === null) continue;
      total.bets += 1;
      const result = ticket.status === "closed" ? (ticket.pnl > 0 ? "won" : ticket.pnl < 0 ? "lost" : "push") : ticket.status;
      if (result === "won") total.wins += 1;
      else if (result === "lost") total.losses += 1;
      else total.pushes += 1;
      total.handle += ticket.stake;
      total.pnl += ticket.pnl;
      total.variance += ticket.variance;
      if (ticket.expected !== null) {
        total.withFair += 1;
        total.fairHandle += ticket.stake;
        total.fairPnl += ticket.pnl;
        total.expected += ticket.expected;
        total.fairVariance += ticket.variance;
      }
    }
    total.roi = total.handle ? total.pnl / total.handle : 0;
    total.ciHalf = total.handle ? (Z_95 * Math.sqrt(total.variance)) / total.handle : 0;
    total.expRoi = total.fairHandle ? total.expected / total.fairHandle : null;
    total.z = total.fairVariance > 0 ? (total.fairPnl - total.expected) / Math.sqrt(total.fairVariance) : null;
    total.winRate = total.wins + total.losses ? total.wins / (total.wins + total.losses) : null;
    return total;
  }

  /** Counts of tickets that are not in P&L: open, and closed/unknown with no result. */
  function exclusions(tickets) {
    let open = 0; let noResult = 0;
    for (const ticket of tickets) {
      if (ticket.status === "open") open += 1;
      else if (ticket.pnl === null) noResult += 1;
    }
    return { open, noResult };
  }

  /** Settled tickets whose Pacific settle day is in [firstDay, lastDay]. */
  function inDayRange(tickets, firstDay, lastDay) {
    return tickets.filter((t) => t.settledDay !== null && t.pnl !== null
      && (firstDay === null || t.settledDay >= firstDay) && t.settledDay <= lastDay);
  }

  /** Every Pacific day from firstDay to lastDay with its totals (days with no bets included). */
  function dailySeries(tickets, firstDay, lastDay) {
    const byDay = new Map();
    for (const ticket of tickets) {
      if (ticket.settledDay === null || ticket.pnl === null) continue;
      if (!byDay.has(ticket.settledDay)) byDay.set(ticket.settledDay, []);
      byDay.get(ticket.settledDay).push(ticket);
    }
    const days = [];
    for (let day = firstDay; day <= lastDay; day = addDays(day, 1)) {
      days.push(Object.assign({ day }, summarize(byDay.get(day) || [])));
    }
    return days;
  }

  /** The earliest settle day among the tickets, or null. */
  function firstSettledDay(tickets) {
    let first = null;
    for (const ticket of tickets) {
      if (ticket.settledDay !== null && (first === null || ticket.settledDay < first)) first = ticket.settledDay;
    }
    return first;
  }

  const RANGE_DAYS = { "7D": 7, "30D": 30, "90D": 90 };
  const DAY_KEY_RE = /^\d{4}-\d{2}-\d{2}$/;

  function isDayKey(value) {
    return typeof value === "string" && DAY_KEY_RE.test(value) && utcToDayKey(dayKeyToUtc(value)) === value;
  }

  /**
   * The Pacific days [first, last] a date-range choice covers. `today` is a
   * Pacific day key; `custom` is { first, last } (either may be missing, and
   * they are swapped if entered backwards); `firstSettled` is the earliest
   * settle day, which "All" starts from.
   *
   * Custom days are clamped to [firstSettled, today]: a date input reports
   * half-typed years such as 0201-10-01, and an unclamped range would make
   * the daily series (one bar per day) millions of days long and hang the page.
   */
  function rangeBounds(range, today, custom, firstSettled) {
    if (range === "Today") return { first: today, last: today };
    if (range === "Yesterday") { const day = addDays(today, -1); return { first: day, last: day }; }
    if (RANGE_DAYS[range]) return { first: addDays(today, 1 - RANGE_DAYS[range]), last: today };
    if (range === "YTD") return { first: today.slice(0, 4) + "-01-01", last: today };
    if (range === "Custom") {
      const floor = firstSettled && firstSettled < today ? firstSettled : today;
      const clamp = (day) => (day < floor ? floor : day > today ? today : day);
      const first = clamp(custom && isDayKey(custom.first) ? custom.first : floor);
      const last = clamp(custom && isDayKey(custom.last) ? custom.last : today);
      return first <= last ? { first, last } : { first: last, last: first };
    }
    if (range === "All") return { first: firstSettled || today, last: today };
    throw new Error("unknown date range " + range + "; expected Today, Yesterday, 7D, 30D, 90D, YTD, All or Custom");
  }

  /** Rows of the Analysis breakdown: one per value of `groupKey`, with totals. */
  function groupBy(tickets, groupKey) {
    const spec = GROUPS.find((g) => g.key === groupKey);
    if (!spec) throw new Error("unknown group " + groupKey + "; expected one of " + GROUPS.map((g) => g.key).join(", "));
    const buckets = new Map();
    for (const ticket of tickets) {
      if (ticket.pnl === null) continue;
      const label = ticket[groupKey] || "Unknown";
      if (!buckets.has(label)) buckets.set(label, []);
      buckets.get(label).push(ticket);
    }
    const rows = [...buckets].map(([label, group]) => Object.assign({ label }, summarize(group)));
    if (spec.order) {
      const rank = (label) => (spec.order.includes(label) ? spec.order.indexOf(label) : spec.order.length);
      return rows.sort((a, b) => rank(a.label) - rank(b.label));
    }
    return rows.sort((a, b) => b.handle - a.handle);
  }

  /** Straight bets with a fair, push/void left out, in 10-point fair-probability bins. */
  function calibration(tickets) {
    const bins = [];
    for (let low = 0; low < 1; low = Math.round((low + 0.1) * 10) / 10) {
      const high = Math.round((low + 0.1) * 10) / 10;
      const inBin = tickets.filter((t) => t.kind === KIND_NAMES.straight && t.fairProb !== null
        && (t.status === "won" || t.status === "lost") && t.fairProb >= low && t.fairProb < high);
      if (!inBin.length) continue;
      const expected = inBin.reduce((sum, t) => sum + t.fairProb, 0) / inBin.length;
      const actual = inBin.filter((t) => t.status === "won").length / inBin.length;
      const ciHalf = Z_95 * Math.sqrt((expected * (1 - expected)) / inBin.length);
      bins.push({ low, high, bets: inBin.length, expected, actual, ciHalf });
    }
    return bins;
  }

  /**
   * Open tickets by whether their game has started: live = event start at or
   * before nowMs (a parlay's is its earliest leg), upcoming = start after
   * nowMs, noStart = no start time (BetOnline's report has none, Kalshi NFL
   * and CFB tickers carry only a date, futures and Kalshi combos none), so
   * they can't be placed in either. Live and upcoming by start, earliest first.
   */
  function splitOpenByStart(openTickets, nowMs) {
    const withStart = [];
    const noStart = [];
    for (const ticket of openTickets) {
      const startMs = parseMs(ticket.eventStart);
      if (startMs === null) noStart.push(ticket);
      else withStart.push({ ticket, startMs });
    }
    withStart.sort((a, b) => a.startMs - b.startMs);
    return {
      live: withStart.filter((row) => row.startMs <= nowMs).map((row) => row.ticket),
      upcoming: withStart.filter((row) => row.startMs > nowMs).map((row) => row.ticket),
      noStart,
    };
  }

  /**
   * A copy of `rows` ordered by `keyOf(row)`: numbers numerically, strings
   * case-insensitively, "asc" or "desc". Rows whose key is null, undefined or
   * NaN go last in either direction, so a "—" never tops a sorted column.
   * Stable: ties keep their incoming order.
   */
  function sortRows(rows, keyOf, direction) {
    const sign = direction === "desc" ? -1 : 1;
    const missing = (key) => key === null || key === undefined || (typeof key === "number" && Number.isNaN(key));
    const keyed = rows.map((row, index) => ({ row, index, key: keyOf(row) }));
    keyed.sort((a, b) => {
      const aMissing = missing(a.key);
      const bMissing = missing(b.key);
      if (aMissing || bMissing) return aMissing === bMissing ? a.index - b.index : (aMissing ? 1 : -1);
      const order = typeof a.key === "number" && typeof b.key === "number"
        ? a.key - b.key
        : String(a.key).localeCompare(String(b.key), undefined, { sensitivity: "base", numeric: true });
      return order !== 0 ? sign * order : a.index - b.index;
    });
    return keyed.map((entry) => entry.row);
  }

  const api = {
    PACIFIC_TZ, GROUPS, WEEKDAYS, NO_FAIR,
    pacificDay, addDays, dayKeyToUtc, weekdayOf, isDayKey, rangeBounds,
    americanToDecimal, decimalToAmerican,
    buildTickets, summarize, exclusions, inDayRange, dailySeries, firstSettledDay, groupBy, calibration,
    venueName, splitOpenByStart, sortRows,
  };
  if (typeof module !== "undefined" && module.exports) module.exports = api;
  else root.UnabatedTrackerStats = api;
})(typeof globalThis !== "undefined" ? globalThis : this);
