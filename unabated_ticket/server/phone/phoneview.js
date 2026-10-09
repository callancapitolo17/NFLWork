// Unabated Ticket phone page — its words and numbers, without the DOM (phone
// page plan steps 2-3). Pure: nothing here places, sends or stores
// anything. Every stake, edge tier, badge and rail word is already in the
// /edges.json row (server/edges_payload.js, from the panel's own
// edgerows.js / betsview.js); this module only formats what the panel's
// panel.js formats inline (prices with cents, start times, line age, the
// Ticket's stake block and contracts), builds the Bets tab's groups from
// /bets.json plus the runner's unmatched list, and turns the Settings form
// into a PUT /settings.json body.
//
// Loaded as a plain <script> in server/phone/index.html after the extension
// modules it reads (exposes globalThis.UnabatedPhoneView) and via require()
// in tests/phoneview.test.js.
//
// Inputs
//   row        one /edges.json line (edges_payload.js rowView)
//   payload    the last /edges.json body; records = /bets.json bets after
//              edgeRows.applyBetsPayload
//   now        epoch ms (the caller's clock; never read here)
// Outputs plain objects and strings.

(function (root) {
  "use strict";

  const inNode = typeof module !== "undefined" && module.exports;
  const kelly = inNode ? require("../../extension/kelly.js") : root.UnabatedKelly;
  const feed = inNode ? require("../../extension/feed.js") : root.UnabatedFeed;
  const betsLib = inNode ? require("../../extension/bets.js") : root.UnabatedBets;
  const betsView = inNode ? require("../../extension/betsview.js") : root.UnabatedBetsView;
  const edgeRows = inNode ? require("../../extension/edgerows.js") : root.UnabatedEdgeRows;

  const MINUTE_MS = 60 * 1000;
  const HOUR_MS = 60 * MINUTE_MS;
  const DAY_MS = 24 * HOUR_MS;
  // Inside these, time to kickoff turns amber, then red (panel.js untilEl).
  const SOON_HOURS = 12;
  const URGENT_HOURS = 2;
  // A price the panel reads as American with no exchange source price.
  const SOURCE_FORMAT_AMERICAN = 1;
  // The polls (plan step 2): edges every 15 s, bets every 30 s while visible.
  const EDGES_POLL_MS = 15 * 1000;
  const BETS_POLL_MS = 30 * 1000;
  // Past this with no successful read, a source is called out as stale.
  const EDGES_STALE_MS = 2 * EDGES_POLL_MS + 15 * 1000;

  const fmtAmerican = edgeRows.fmtAmerican;
  const fmtDollars = edgeRows.fmtDollars;

  // ---- formatting (panel.js's own helpers, with the clock passed in) ----------

  // "+3.21%" from Unabated's edge percent (3.21).
  function fmtEdgePct(edgePct) {
    return `${edgePct > 0 ? "+" : ""}${edgePct.toFixed(2)}%`;
  }

  // "+215 · 31.7¢": the American price and the prediction-market cents, from
  // the exchange's exact source price when the line carries one (kelly.bookProbOf).
  function fmtPriceBoth(price, sourceFormat, sourcePrice) {
    const cents = kelly.bookProbOf({ bookPrice: price, sourceFormat, sourcePrice }) * 100;
    return `${fmtAmerican(price)} · ${cents.toFixed(1)}¢`;
  }

  function fmtPoints(points) {
    return points == null ? "" : `${points > 0 ? "+" : ""}${points}`;
  }

  // Unabated's eventStart is naive UTC ("2026-09-12T23:30:00"); shown in the phone's zone.
  function fmtStart(eventStart) {
    if (!eventStart) return "";
    const date = new Date(eventStart.endsWith("Z") ? eventStart : `${eventStart}Z`);
    if (Number.isNaN(date.getTime())) return eventStart;
    return date.toLocaleString([], { weekday: "short", month: "short", day: "numeric", hour: "numeric", minute: "2-digit" });
  }

  function fmtUntil(startMs, now) {
    const minutes = Math.round((startMs - now) / MINUTE_MS);
    if (minutes < 0) return "started";
    if (minutes < 60) return `in ${minutes}m`;
    if (minutes < 48 * 60) return `in ${Math.floor(minutes / 60)}h ${minutes % 60}m`;
    return `in ${Math.round(minutes / 1440)}d`;
  }

  // "soon urgent" inside 2 h, "soon" inside 12 h, "" otherwise (the CSS classes).
  function untilLevel(startMs, now) {
    if (!Number.isFinite(startMs)) return "";
    const hours = (startMs - now) / HOUR_MS;
    if (hours <= URGENT_HOURS) return "soon urgent";
    if (hours <= SOON_HOURS) return "soon";
    return "";
  }

  // "just now" / "20s ago" / "4m ago" / "3h ago".
  function fmtAge(ms) {
    if (ms < 1000) return "just now";
    if (ms < MINUTE_MS) return `${Math.round(ms / 1000)}s ago`;
    if (ms < HOUR_MS) return `${Math.round(ms / MINUTE_MS)}m ago`;
    return `${Math.round(ms / HOUR_MS)}h ago`;
  }

  // How long ago the book last changed this line: the stale-line tell.
  function fmtLineAge(modifiedMs, now) {
    if (modifiedMs == null) return "line age unknown";
    const ms = now - modifiedMs;
    if (ms < MINUTE_MS) return "line just changed";
    if (ms < HOUR_MS) return `line ${Math.round(ms / MINUTE_MS)}m old`;
    if (ms < 2 * DAY_MS) return `line ${Math.round(ms / HOUR_MS)}h old`;
    return `line ${Math.round(ms / DAY_MS)}d old`;
  }

  function fmtLiquidity(value) {
    return value == null ? "" : `liq ${value.toLocaleString("en-US", { style: "currency", currency: "USD", maximumFractionDigits: 0 })}`;
  }

  function fmtClock(ms) {
    return new Date(ms).toLocaleTimeString([], { hour: "numeric", minute: "2-digit" });
  }

  // "Spread · 1H · alt of -3 · rot 466"
  function marketText(row) {
    return [
      row.betType,
      row.period && row.period !== "FG" ? row.period : null,
      row.isAlt ? `alt of ${fmtPoints(row.mainPoints)}` : null,
      row.rotation != null ? `rot ${row.rotation}` : null,
    ].filter(Boolean).join(" · ");
  }

  // "Chicago Bears @ Carolina Panthers · NFL", else Unabated's event name.
  function describeMatchup(row) {
    const league = row.leagueLabel || String(row.league || "").toUpperCase();
    if (row.awayTeam && row.homeTeam) return `${row.awayTeam} @ ${row.homeTeam}${league ? ` · ${league}` : ""}`;
    return row.eventName || "";
  }

  // "Novig +215 · 31.7¢"
  function bookPriceText(row) {
    return `${row.book.name} ${fmtPriceBoth(row.price, row.sourceFormat, row.sourcePrice)}`;
  }

  // "line 4m old · liq $1,200"
  function lineAgeText(row, now) {
    return [fmtLineAge(row.modifiedMs, now), fmtLiquidity(row.liquidity)].filter(Boolean).join(" · ");
  }

  // ---- the ticket sheet --------------------------------------------------------

  // "1,127 contracts @ 23.2¢" + its cost, for a line priced in contracts
  // (Kalshi, Novig, Polymarket); null for a sportsbook (panel.js renderContracts).
  function contractsView(acted, row) {
    if (acted == null || acted <= 0) return null;
    const order = kelly.contractOrder({ stake: acted, bookPrice: row.price, sourceFormat: row.sourceFormat, sourcePrice: row.sourcePrice, bookName: row.book ? row.book.name : null });
    if (!order) return null;
    const priceText = `${order.priceCents.toFixed(1)}¢`;
    if (order.contracts === 0) return { text: `under 1 contract @ ${priceText}`, cost: null, under: true };
    const plural = order.contracts === 1 ? "" : "s";
    return { text: `${order.contracts.toLocaleString("en-US")} contract${plural} @ ${priceText}`, cost: fmtDollars(order.costDollars), under: false };
  }

  // The Ticket's stake block for a row (panel.js renderTicket +
  // renderStakeExposure): the label, the number to act on, the position line
  // under it, the uncapped stake when liquidity cut it, the contracts and the
  // payout of that number.
  function stakeBlockView(row) {
    const advice = row.advice;
    const words = betsView.stakeAdviceWords(advice);
    const acted = betsView.suggestedBetAmount(advice);
    const teasers = advice.teasers || { held: 0, against: 0 };
    const position = [
      advice.held > 0 ? `held ${betsLib.formatStake(advice.held)}` : null,
      advice.against > 0 ? `against ${betsLib.formatStake(advice.against)}` : null,
      teasers.held > 0 ? `teasers ${betsLib.formatStake(teasers.held)}` : null,
      teasers.against > 0 ? `teasers against ${betsLib.formatStake(teasers.against)}` : null,
      words ? words.cap : null,
      words ? words.alone : null,
    ].filter(Boolean);
    const label = !words ? "Bet"
      : advice.bet === 0 && advice.cappedAt === 0 ? "Nothing resting at this price"
        : advice.bet === 0 ? "Already at full size"
          : words.verb === "add" ? "Add to your position" : "Bet";
    const payout = acted > 0 ? acted * kelly.americanToDecimal(row.price) : null;
    return {
      label, stake: fmtDollars(acted), noEdge: acted <= 0,
      position: position.join(" · "),
      positionAgainst: advice.held === 0 && teasers.held === 0 && (advice.against > 0 || teasers.against > 0),
      contracts: contractsView(acted, row),
      uncapped: betsView.uncappedLine(advice, { price: row.price, sourceFormat: row.sourceFormat, sourcePrice: row.sourcePrice, bookName: row.book ? row.book.name : null }),
      toWin: payout == null ? null : fmtDollars(payout - acted),
      payout: payout == null ? null : fmtDollars(payout),
    };
  }

  // Everything the ticket sheet shows for one row.
  function ticketView(row, now) {
    return {
      sideLabel: row.sideLabel,
      betLine: `${marketText(row)}${row.isAlt ? " · alt line" : ""}`,
      matchup: describeMatchup(row),
      start: `${fmtStart(row.eventStart)} · ${fmtUntil(row.eventStartMs, now)}`,
      book: row.book.name,
      price: fmtPriceBoth(row.price, row.sourceFormat, row.sourcePrice),
      lineAge: lineAgeText(row, now),
      fair: row.fair == null ? "unknown" : fmtPriceBoth(row.fair, SOURCE_FORMAT_AMERICAN, null),
      edge: fmtEdgePct(row.edgePct),
      edgeTier: row.edgeTier,
      stakeBlock: stakeBlockView(row),
    };
  }

  // ---- status at the top ---------------------------------------------------------

  // "NFL · CFB · 2,140 lines (+9,812 alts) · snapshot built 3m ago"
  function scannerText(scanner, now) {
    if (!scanner || scanner.phase === "starting") return "Runner starting the scanner…";
    const loaded = scanner.leaguesLoaded || [];
    const label = (id) => (feed.LEAGUES[id] || { label: `league ${id}` }).label;
    const leagues = loaded.length > 4 ? `${loaded.length} leagues` : loaded.map(label).join(" · ");
    if (scanner.phase === "loading" && !loaded.length) {
      return `Loading snapshots…${scanner.loading ? ` ${scanner.loading.done}/${scanner.loading.total} leagues` : ""}`;
    }
    const stale = scanner.staleLeagues && scanner.staleLeagues.length ? ` (${scanner.staleLeagues.length} stale: ${scanner.staleLeagues.map(label).join(", ")})` : "";
    const built = scanner.snapshotBuiltAt ? `snapshot built ${fmtAge(now - scanner.snapshotBuiltAt)}${stale}` : "no data yet";
    const alts = scanner.altLineCount ? ` (+${scanner.altLineCount.toLocaleString("en-US")} alts)` : "";
    return `${leagues || "no leagues"} · ${(scanner.lineCount || 0).toLocaleString("en-US")} lines${alts} · ${built}`;
  }

  // "4 books: Kalshi, Novig, …" / "all 12 live books" / "default books (9)"
  function booksText(books) {
    if (!books) return "";
    if (!books.ids) return `all ${books.liveCount} live books`;
    const names = (books.names || []).join(", ");
    return `${books.mode === "default" ? "default books" : "your books"} (${books.ids.length}): ${names}`;
  }

  // What one poll left: {payload, okAt, error, errorAt, failingSince}.
  //   "unreachable since 2:14 PM (6 min): <error>"
  function failingText(what, poll, now) {
    const since = poll.failingSince ?? poll.errorAt;
    return `${what} unreachable since ${fmtClock(since)} (${betsView.fmtAgeShort(now - since)}): ${poll.error}`;
  }

  // The red / amber strips above the views, worst first: the page's own two
  // reads, then what the runner says about its feed, the bets service and
  // its settings. [{level: "bad" | "warn", text}]
  function banners(edgesPoll, betsPoll, now) {
    const out = [];
    if (edgesPoll && edgesPoll.error) out.push({ level: "bad", text: failingText("Edges", edgesPoll, now) });
    if (betsPoll && betsPoll.error) out.push({ level: "bad", text: failingText("Bets service", betsPoll, now) });
    const payload = edgesPoll && edgesPoll.payload;
    if (!payload) return out;
    if (payload.scanner && payload.scanner.error) out.push({ level: "bad", text: `Scanner: ${payload.scanner.error}` });
    const service = payload.betsService || {};
    if (service.error) {
      const since = service.unreachableSince != null ? ` since ${fmtClock(service.unreachableSince)}` : "";
      out.push({ level: "bad", text: `Runner cannot read the bets service${since}: stakes ignore your open bets. ${service.error}` });
    }
    const settings = payload.settings || {};
    if (settings.error) out.push({ level: "warn", text: `Runner is using its ${settings.source} settings: ${settings.error}` });
    if (edgesPoll.okAt != null && now - edgesPoll.okAt > EDGES_STALE_MS && !edgesPoll.error) {
      out.push({ level: "warn", text: `Edges last read ${fmtAge(now - edgesPoll.okAt)}` });
    }
    return out;
  }

  // "edges 4s ago · bets 12s ago"
  function freshnessText(edgesPoll, betsPoll, now) {
    const part = (name, poll) => {
      if (!poll || poll.okAt == null) return `${name} —`;
      return `${name} ${fmtAge(now - poll.okAt)}`;
    };
    return [part("edges", edgesPoll), part("bets", betsPoll)].join(" · ");
  }

  // ---- Bets tab --------------------------------------------------------------------

  // One bet as a list item: the pick, venue · game, stake and when.
  function betItemView(bet) {
    const pinned = Boolean(bet.pin);
    const game = pinned && bet.pin.awayTeamName && bet.pin.homeTeamName ? `${bet.pin.awayTeamName} @ ${bet.pin.homeTeamName}`
      : bet.awayTeam && bet.homeTeam ? `${bet.awayTeam} @ ${bet.homeTeam}` : bet.awayTeam || bet.homeTeam || null;
    const venue = bet.venue ? betsLib.venueLabel(bet.venue) : "unknown venue";
    return {
      id: bet.id,
      what: betsLib.describeBet(bet),
      meta: [venue, bet.league ? String(bet.league).toUpperCase() : null, game].filter(Boolean).join(" · "),
      stake: bet.stake == null ? "—" : fmtDollars(bet.stake),
      when: bet.placedAt ? betsLib.formatPlacedAt(bet.placedAt) : "",
      pinned,
    };
  }

  // "Every open game bet matches a game on the board." / "3 of 5 …", or null
  // when there is no open game bet; says so when it could not be checked.
  function matchNoteText(unmatched, boardLineCount, gameBets, unmatchedGameBets) {
    if (gameBets === 0) return null;
    if (unmatched == null) return "Edges runner not read yet, so which game each bet is on is not checked.";
    if (!boardLineCount) return "Waiting for the runner's board to load before checking which game each bet is on.";
    if (unmatchedGameBets === 0) return "Every open game bet matches a game on the board.";
    return `${gameBets - unmatchedGameBets} of ${gameBets} open game bets match a game on the board.`;
  }

  // The Bets tab as the panel groups it (panel.js renderBets): money at
  // risk, Needs a game, Needs a code fix, the open list (unmatched ones
  // marked) and the folded Not on the board. `unmatched` is the runner's
  // list ([{betId, reason, attachable, needsGame, needsFix}]) or null when
  // the runner was not read, in which case nothing is grouped as unmatched.
  function betsTabView(records, unmatched, boardLineCount, sourcesPayload, now) {
    const open = records.filter((bet) => bet.status === "open")
      .sort((a, b) => Date.parse(b.placedAt || 0) - Date.parse(a.placedAt || 0));
    const byId = new Map(open.map((bet) => [bet.id, bet]));
    const entries = (unmatched || []).filter((entry) => byId.has(entry.betId))
      .map((entry) => ({ ...entry, item: betItemView(byId.get(entry.betId)) }));
    const unmatchedIds = new Set(entries.map((entry) => entry.betId));
    const priced = open.filter((bet) => typeof bet.stake === "number");
    const venues = betsView.sourceRows(sourcesPayload, now).filter((row) => row.configured).length;
    const isGameBet = (bet) => !bet.unmatchable || betsLib.isParseFailure(bet);
    const gameBets = open.filter(isGameBet).length;
    const unmatchedGameBets = entries.filter((entry) => isGameBet(byId.get(entry.betId))).length;
    return {
      atRisk: fmtDollars(priced.reduce((total, bet) => total + bet.stake, 0)),
      caption: [
        `at risk · ${open.length} open bet${open.length === 1 ? "" : "s"}`,
        `${venues} venue${venues === 1 ? "" : "s"}`,
        priced.length === open.length ? null : `${open.length - priced.length} with no stake reported`,
      ].filter(Boolean).join(" · "),
      needsGame: entries.filter((entry) => entry.needsGame),
      needsFix: entries.filter((entry) => entry.needsFix),
      offBoard: entries.filter((entry) => !entry.needsGame && !entry.needsFix),
      open: open.map((bet) => ({ item: betItemView(bet), unmatched: unmatchedIds.has(bet.id) })),
      matchNote: matchNoteText(unmatched, boardLineCount, gameBets, unmatchedGameBets),
      venueRows: betsView.sourceRows(sourcesPayload, now),
    };
  }

  // ---- Settings ----------------------------------------------------------------------

  // Every field of bets.duckdb::edge_settings with the panel's default for it
  // (edgerows.js): what a null shows as and what "reset" returns to.
  function settingsDefaults() {
    const edges = edgeRows.DEFAULT_EDGE_SETTINGS;
    return {
      ...edgeRows.DEFAULT_STAKE_SETTINGS,
      leagues: edges.leagues, periods: edges.periods, betTypes: edges.betTypes,
      bookMode: "default", bookIds: null,
      minEdgePct: edges.minEdgePct, minStake: edges.minStake, maxLineAgeHours: edges.maxLineAgeHours,
      minLiquidityToWin: edges.minLiquidityToWin, includeAlts: edges.includeAlts, sortBy: edges.sortBy,
      groupByMarket: edges.groupByMarket,
    };
  }

  // The value each field is in effect: the stored one, else the default.
  function effectiveSettings(held) {
    const defaults = settingsDefaults();
    const out = {};
    for (const field of Object.keys(defaults)) {
      const stored = held ? held[field] : null;
      out[field] = stored == null ? defaults[field] : stored;
    }
    return out;
  }

  function sameValue(a, b) {
    const norm = (value) => (Array.isArray(value) ? JSON.stringify([...value].sort((x, y) => x - y)) : JSON.stringify(value));
    return norm(a) === norm(b);
  }

  // The PUT /settings.json `settings` object for a submitted form: only the
  // fields whose value differs from what is in effect, so an untouched field
  // keeps following the default. A field the form left undefined is not
  // sent. bookMode and bookIds travel together (the service requires bookIds
  // exactly when bookMode is "custom").
  //   form  {field: value | undefined}
  function settingsUpdate(held, form) {
    const current = effectiveSettings(held);
    const update = {};
    for (const [field, value] of Object.entries(form)) {
      if (value === undefined || field === "bookMode" || field === "bookIds") continue;
      if (!sameValue(value, current[field])) update[field] = value;
    }
    const mode = form.bookMode === undefined ? current.bookMode : form.bookMode;
    const ids = mode === "custom" ? (form.bookIds === undefined ? current.bookIds || [] : form.bookIds) : null;
    if (mode !== current.bookMode || !sameValue(ids, current.bookIds)) {
      update.bookMode = mode;
      update.bookIds = ids;
    }
    return update;
  }

  // The body that resets one field to the panel's default; a book choice resets both halves.
  function resetUpdate(field) {
    return field === "bookMode" || field === "bookIds" ? { bookMode: null, bookIds: null } : { [field]: null };
  }

  // A field's default as shown beside it ("default 30000", "default all sports").
  function defaultText(field) {
    const value = settingsDefaults()[field];
    if (field === "leagues") return "default: every sport";
    if (field === "periods") return `default: ${value.map((id) => feed.PERIODS[id]).join(", ")}`;
    if (field === "betTypes") return `default: ${value.map((id) => feed.BET_TYPES[id]).join(", ")}`;
    if (field === "bookMode") return "default: the default books";
    if (typeof value === "boolean") return `default: ${value ? "on" : "off"}`;
    return `default: ${value}`;
  }

  // Sports ticked for a league list (a sport is on when any of its leagues is, as the panel's checkboxes read).
  function sportsOfLeagues(leagues) {
    return Object.keys(feed.SPORTS).filter((sport) => feed.leagueIdsOfSport(sport).some((id) => leagues.includes(id)));
  }

  // The league list for the sports ticked: `currentLeagues` itself while the
  // ticks are the sports it already covers (so a stored NFL-only list is not
  // widened to all football by a save that never touched it), else every
  // league of each ticked sport, as the panel's checkboxes write it.
  function leaguesForSports(sports, currentLeagues) {
    const unchanged = sports.slice().sort().join(",") === sportsOfLeagues(currentLeagues).sort().join(",");
    if (unchanged) return currentLeagues;
    return sports.flatMap((sport) => feed.leagueIdsOfSport(sport));
  }

  // The runner's unmatched list with the page's own Dismiss / Restore marks
  // laid over it, so a mark shows at once rather than when the runner next
  // reads the bets service. `dismissedIds` null (an older service) keeps the
  // runner's flags as they are.
  function applyDismissals(unmatched, dismissedIds) {
    if (!Array.isArray(unmatched) || !Array.isArray(dismissedIds)) return unmatched;
    const dismissed = new Set(dismissedIds);
    return unmatched.map((entry) => {
      const isDismissed = dismissed.has(entry.betId);
      const flag = entry.flagUnlessDismissed;
      return { ...entry, dismissed: isDismissed, needsGame: flag === "game" && !isDismissed, needsFix: flag === "fix" && !isDismissed };
    });
  }

  // ---- ticket extras ------------------------------------------------------------

  // The panel's "Line moved" note: the ticket's line now at another price or
  // number than when it was opened, or null. `opened` {price, points}.
  function lineMovedText(opened, row) {
    if (!opened || !row) return null;
    if (opened.price === row.price && opened.points === row.points) return null;
    const at = (points) => (points != null ? ` at ${fmtPoints(points)}` : "");
    return `Line moved: now ${fmtAmerican(row.price)}${at(row.points)} (opened ${fmtAmerican(opened.price)}${at(opened.points)}). Stake re-sized.`;
  }

  // ---- Teasers tab ---------------------------------------------------------------

  // "Bills @ Jets · NFL · Sun, Oct 11, 1:00 PM · " (the time to kickoff follows).
  function teaserLegMeta(rowView) {
    const start = Number.isFinite(rowView.eventStartMs) ? fmtStart(new Date(rowView.eventStartMs).toISOString()) : "";
    return `${rowView.matchup} · ${rowView.leagueLabel} · ${start} · `;
  }

  // The bets service's BFA row (betsview.sourceRows), or null before /bets.json was read.
  function bfaRowOf(betsPoll, now) {
    if (!betsPoll || betsPoll.okAt == null) return null;
    return betsView.sourceRows(betsPoll.payload, now).find((row) => row.venue === "bfa") || null;
  }

  // The Teasers tab's amber note (panel.js renderTeasersBanners), or null.
  function teasersWarningText(bfaRow) {
    if (!bfaRow) return null;
    if (!bfaRow.configured) return "The bets service reads no BFA account, so tickets already placed at Buckeye are not known here.";
    if (bfaRow.error) return `BFA's last pull failed (${bfaRow.error}); open teasers are as of ${bfaRow.ageText} ago.`;
    return null;
  }

  const api = {
    EDGES_POLL_MS, BETS_POLL_MS, EDGES_STALE_MS,
    fmtEdgePct, fmtPriceBoth, fmtPoints, fmtStart, fmtUntil, untilLevel, fmtAge, fmtLineAge, fmtLiquidity, fmtClock,
    marketText, describeMatchup, bookPriceText, lineAgeText, contractsView, stakeBlockView, ticketView,
    scannerText, booksText, banners, freshnessText, betItemView, betsTabView, applyDismissals,
    lineMovedText, teaserLegMeta, bfaRowOf, teasersWarningText,
    settingsDefaults, effectiveSettings, settingsUpdate, resetUpdate, defaultText, sportsOfLeagues, leaguesForSports,
    fmtDollars, fmtAmerican,
  };

  if (inNode) {
    module.exports = api;
  } else {
    root.UnabatedPhoneView = api;
  }
})(typeof globalThis !== "undefined" ? globalThis : this);
