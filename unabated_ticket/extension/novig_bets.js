// Novig bet source for the Unabated Ticket panel (#116): the pure normaliser
// that turns the portfolio responses the Novig web app fetches for itself
// (mirrored by novig_page.js, delivered by novig_content.js) into the
// normalised bet records bets.js matches. Pure: no DOM, no fetch, no chrome.*
// — loaded as a plain <script> in the Novig tab's isolated world (before
// novig_content.js; exposes globalThis.UnabatedNovigBets) and via require()
// in tests/novig_bets.test.js. Team keys are left null: the panel fills them
// with bets.resolveTeamKeys() so the team table lives only in teams.js.
//
// Input   {orders, parlays, readAt} — the rows of ActivePortfolioOrders_Query
//         / SettledPortfolioOrders_Query and ParlayPortfolioQuery (Apollo
//         queries the app POSTs to api.novig.us/v1/graphql).
// Output  records in the contract of docs/2026-09-11-issue-114-bet-history-plan.md
//         (id "novig:<order id>" or "novig:<parlay id>:<leg index>", source
//         "novig_page", venue "novig").
//
// Novig facts these rules rest on (read off the app.novig.us bundle, 2026-09-11;
// see tests/fixtures/bets/novig_bets.json "provenance"):
//   order         qty = contracts still resting, originalQty = placed, so the
//                 matched size is originalQty - qty; price is a 0-1
//                 probability; one contract pays $1; isBid true BACKS the
//                 outcome, false LAYS it (the card grades an ask on a LOSS
//                 outcome as a win) — in a two-way market a lay at p is the
//                 other side at 1 - p. status PENDING/OPEN/FILLED/CANCELED/
//                 REJECTED; fills [{cost, qty, isWash, isTaker, created_at}].
//   outcome       index 0 = HOME or OVER, 1 = AWAY or UNDER (HomeAwayIndex /
//                 OverUnderIndex); description carries the side's own number
//                 ("CHI -13.5", "Over 47.5"); status TBD/WIN/LOSS/PUSH.
//   market        type MONEY/SPREAD/TOTAL (+ _1H = first half, which on MLB is
//                 the first five innings — the SGP scrapers' finding), status
//                 OPEN/CLOSED/SETTLED, strike is the market's number in the
//                 home perspective; cash_out_requests[0] APPROVED = cashed out.
//   event/game    event.type GAME|FUTURE, game.league "NFL"/"NCAAF"/"MLB"/...,
//                 game.awayTeam / homeTeam {name, symbol, short_name},
//                 scheduled_start ISO with a +00:00 offset.
//   parlay        status FILLED (open) / WIN / LOSS / PUSH / UNFILLED; wager in
//                 dollars; each leg {price, outcome{... market{... event}}}.

(function (root) {
  "use strict";

  const SOURCE = "novig_page";
  const VENUE = "novig";
  const EASTERN = "America/New_York";

  // Novig game.league -> the feed.LEAGUES path the scanner uses.
  const LEAGUES = {
    NFL: "nfl", NCAAF: "cfb", NBA: "nba", NCAAB: "cbb", WNBA: "wnba", MLB: "mlb", NHL: "nhl",
    MLS: "soccer", EPL: "soccer", Bundesliga: "soccer", "Serie A": "soccer", "La Liga": "soccer",
    "Ligue 1": "soccer", "Champions League": "soccer", "Europa League": "soccer", "FIFA Club World Cup": "soccer",
  };
  // Novig market.type -> {betType, half}. Everything else is "other".
  const MARKET_TYPES = {
    MONEY: { betType: "moneyline", half: false },
    SPREAD: { betType: "spread", half: false },
    TOTAL: { betType: "total", half: false },
    MONEY_1H: { betType: "moneyline", half: true },
    SPREAD_1H: { betType: "spread", half: true },
    TOTAL_1H: { betType: "total", half: true },
  };
  const REASON_NOT_GAME = "not a game market";
  const APPROX_UNMATCHED = "novig_order_unmatched";
  const APPROX_PENDING = "novig_order_pending";
  const DESCRIPTION_NUMBER_RE = /([+-]?\d+(?:\.\d+)?)\s*$/;

  // ---- helpers -----------------------------------------------------------------

  function toNumber(value) {
    const n = typeof value === "number" ? value : parseFloat(value);
    return Number.isFinite(n) ? n : null;
  }

  // Novig sends "+00:00" offsets; the contract wants "...Z".
  function isoUtc(value) {
    const ms = typeof value === "string" ? Date.parse(value) : NaN;
    return Number.isFinite(ms) ? new Date(ms).toISOString() : null;
  }

  function easternDateOf(iso) {
    const ms = Date.parse(iso);
    if (!Number.isFinite(ms)) return null;
    const parts = {};
    const formatter = new Intl.DateTimeFormat("en-US", { timeZone: EASTERN, year: "numeric", month: "2-digit", day: "2-digit" });
    for (const part of formatter.formatToParts(new Date(ms))) parts[part.type] = part.value;
    return `${parts.year}-${parts.month}-${parts.day}`;
  }

  function probabilityToAmerican(probability) {
    if (!(probability > 0 && probability < 1)) return null;
    if (probability >= 0.5) return -Math.round((probability / (1 - probability)) * 100);
    return Math.round(((1 - probability) / probability) * 100);
  }

  function roundCents(dollars) {
    return Math.round(dollars * 100) / 100;
  }

  function teamName(team) {
    if (!team) return null;
    return team.name || team.short_name || team.symbol || null;
  }

  function matchupOf(market) {
    const game = market.event && market.event.game;
    return game ? `${teamName(game.awayTeam)} @ ${teamName(game.homeTeam)}` : null;
  }

  // ---- side / points -----------------------------------------------------------

  // Which team an outcome names, as 0 = away / 1 = home: the competitor symbol
  // against the game's teams first, Novig's index convention (0 = home)
  // second. Null when neither resolves — a wrong side would flag the wrong line.
  function teamIndexOf(outcome, game) {
    const symbol = outcome.competitor ? outcome.competitor.symbol : null;
    if (symbol && game.homeTeam && game.awayTeam) {
      if (symbol === game.homeTeam.symbol) return 1;
      if (symbol === game.awayTeam.symbol) return 0;
    }
    if (outcome.index === 0) return 1;
    if (outcome.index === 1) return 0;
    return null;
  }

  function descriptionNumber(description) {
    const m = typeof description === "string" ? DESCRIPTION_NUMBER_RE.exec(description.trim()) : null;
    return m ? toNumber(m[1]) : null;
  }

  // {side, points} of the OUTCOME (what a bid backs), or an error string.
  function outcomeSide(betType, outcome, market, game) {
    if (betType === "total") {
      const description = typeof outcome.description === "string" ? outcome.description.trim().toLowerCase() : "";
      const side = description.startsWith("over") ? "over" : description.startsWith("under") ? "under"
        : outcome.index === 0 ? "over" : outcome.index === 1 ? "under" : null;
      const points = descriptionNumber(outcome.description) ?? toNumber(market.strike);
      if (!side) return `total outcome ${outcome.id} names neither over nor under`;
      if (points == null) return `total outcome ${outcome.id} has no number`;
      return { side, points };
    }
    const teamIndex = teamIndexOf(outcome, game);
    if (teamIndex == null) return `outcome ${outcome.id} names no team`;
    const side = teamIndex === 0 ? "away" : "home";
    if (betType === "moneyline") return { side, points: null };
    const fromDescription = descriptionNumber(outcome.description);
    const strike = toNumber(market.strike);
    const points = fromDescription ?? (strike == null ? null : teamIndex === 1 ? strike : -strike);
    if (points == null) return `spread outcome ${outcome.id} has no number`;
    return { side, points };
  }

  // A lay is the other side: away<->home / over<->under, a spread's number negated.
  function laySide(betType, backed) {
    const flip = { away: "home", home: "away", over: "under", under: "over" };
    return { side: flip[backed.side], points: betType === "spread" ? -backed.points : backed.points };
  }

  // ---- status ------------------------------------------------------------------

  function gradeSettled(outcomeStatus, isBid) {
    if (outcomeStatus === "PUSH") return "push";
    if (outcomeStatus === "WIN") return isBid ? "won" : "lost";
    if (outcomeStatus === "LOSS") return isBid ? "lost" : "won";
    return "unknown";
  }

  // Port of the app's getOrderStateFromCardData onto the contract's vocabulary.
  function orderStatus(order, market, outcome, matchedQty) {
    const fills = Array.isArray(order.fills) ? order.fills : [];
    if (order.status === "REJECTED") return "void";
    if (order.qty === 0 && fills.length > 0 && fills.every((fill) => fill.isWash)) return "closed";
    const cashOut = Array.isArray(market.cash_out_requests) && market.cash_out_requests.length ? market.cash_out_requests[0] : null;
    if (cashOut && cashOut.status === "APPROVED" && Date.parse(cashOut.created_at) > Date.parse(order.created_at)) return "closed";
    if (market.status === "SETTLED") {
      if (order.status === "CANCELED" && matchedQty <= 0) return "void";
      return gradeSettled(outcome.status, order.isBid === true);
    }
    if (order.status === "CANCELED") return matchedQty > 0 ? "open" : "void";
    if (order.status === "OPEN" || order.status === "FILLED" || order.status === "PENDING") return "open";
    return "unknown";
  }

  function parlayStatus(status) {
    if (status === "FILLED") return "open";
    if (status === "WIN") return "won";
    if (status === "LOSS") return "lost";
    if (status === "PUSH") return "push";
    if (status === "UNFILLED") return "void";
    return "unknown";
  }

  // ---- records -----------------------------------------------------------------

  function unmatchable(base, reason) {
    return Object.assign(base, {
      league: null, eventStart: null, eventDate: null, awayTeam: null, homeTeam: null, awayKey: null, homeKey: null,
      betType: "other", period: null, side: null, points: null, unmatchable: reason,
    });
  }

  // league / event / teams / market shape shared by an order and a parlay leg.
  function gameFields(market, outcome, isBid) {
    const event = market.event || {};
    const game = event.game || null;
    if (event.type === "FUTURE" || market.player) return { unmatchable: REASON_NOT_GAME };
    const spec = MARKET_TYPES[market.type];
    if (!spec) return { unmatchable: REASON_NOT_GAME };
    if (!game || !game.awayTeam || !game.homeTeam) return { unmatchable: `unreadable Novig order (no teams on event ${event.id || "?"})` };
    const novigLeague = game.league || event.league || null;
    const league = LEAGUES[novigLeague] || null;
    if (!league) return { unmatchable: `league not supported (${novigLeague || "unknown"})` };
    const backed = outcomeSide(spec.betType, outcome, market, game);
    if (typeof backed === "string") return { unmatchable: `unreadable Novig order (${backed})` };
    const taken = isBid ? backed : laySide(spec.betType, backed);
    const eventStart = isoUtc(game.scheduled_start || event.scheduled_start);
    return {
      league, eventStart, eventDate: eventStart ? easternDateOf(eventStart) : null,
      awayTeam: teamName(game.awayTeam), homeTeam: teamName(game.homeTeam), awayKey: null, homeKey: null,
      betType: spec.betType, period: spec.half ? (league === "mlb" ? "F5" : "1H") : "FG",
      side: taken.side, points: taken.points, unmatchable: null,
    };
  }

  function normalizeOrder(order, readAt) {
    const market = order.market || {};
    const outcome = order.outcome || {};
    const originalQty = toNumber(order.originalQty) ?? 0;
    const remainingQty = toNumber(order.qty) ?? 0;
    const matchedQty = Math.max(originalQty - remainingQty, 0);
    const isBid = order.isBid === true;
    const outcomePrice = toNumber(order.price);
    // The side taken pays price per contract on a bid, 1 - price on a lay.
    const takenPrice = outcomePrice == null ? null : isBid ? outcomePrice : 1 - outcomePrice;
    const status = orderStatus(order, market, outcome, matchedQty);
    const approx = [];
    // A resting or pending order is a bet you are trying to place: flagged as
    // open on its full size, with the caveat that nothing has matched yet.
    if (status === "open" && matchedQty === 0) approx.push(order.status === "PENDING" ? APPROX_PENDING : APPROX_UNMATCHED);
    const sizedQty = status === "open" && matchedQty === 0 ? originalQty : matchedQty;
    const fills = Array.isArray(order.fills) ? order.fills : [];
    const base = {
      id: `novig:${order.id}`,
      source: SOURCE,
      venue: VENUE,
      rotation: null,
      price: takenPrice == null ? null : probabilityToAmerican(takenPrice),
      stake: takenPrice == null ? null : roundCents(sizedQty * takenPrice),
      toWin: takenPrice == null ? null : roundCents(sizedQty * (1 - takenPrice)),
      contracts: sizedQty,
      placedAt: isoUtc(order.created_at),
      status,
      closedAt: status === "open" ? null : isoUtc(order.updated_at) || isoUtc(order.created_at),
      isParlayLeg: false, parlayId: null, legIndex: null, legCount: null,
      approx,
      sourceFetchedAt: readAt,
      raw: {
        orderId: order.id, orderStatus: order.status, isBid, outcomePrice, originalQty, remainingQty,
        fillCount: fills.length, marketType: market.type || null, marketStatus: market.status || null,
        strike: toNumber(market.strike), outcomeIndex: outcome.index ?? null,
        outcomeDescription: outcome.description || null, outcomeStatus: outcome.status || null,
        novigLeague: market.event && market.event.game ? market.event.game.league : null,
        marketTitle: [matchupOf(market), market.type, outcome.description].filter(Boolean).join(" · "),
      },
    };
    const game = gameFields(market, outcome, isBid);
    if (game.unmatchable) return unmatchable(base, game.unmatchable);
    return Object.assign(base, game);
  }

  function normalizeParlay(parlay, readAt) {
    const legs = Array.isArray(parlay.legs) ? parlay.legs : [];
    const status = parlayStatus(parlay.status);
    const wager = toNumber(parlay.wager);
    const parlayPrice = toNumber(parlay.price);
    return legs.map((leg, legIndex) => {
      const outcome = leg.outcome || {};
      const market = outcome.market || {};
      const legPrice = toNumber(leg.price);
      const base = {
        id: `novig:${parlay.id}:${legIndex}`,
        source: SOURCE,
        venue: VENUE,
        rotation: null,
        price: legPrice == null ? null : probabilityToAmerican(legPrice),
        stake: wager,
        toWin: wager == null || !(parlayPrice > 0) ? null : roundCents(wager / parlayPrice - wager),
        contracts: null,
        placedAt: isoUtc(parlay.created_at),
        status,
        closedAt: status === "open" ? null : isoUtc(parlay.updated_at) || isoUtc(parlay.created_at),
        isParlayLeg: true, parlayId: `novig:${parlay.id}`, legIndex, legCount: legs.length,
        approx: [],
        sourceFetchedAt: readAt,
        raw: {
          parlayId: parlay.id, parlayStatus: parlay.status, parlayPrice, legPrice,
          marketType: market.type || null, strike: toNumber(market.strike), outcomeIndex: outcome.index ?? null,
          outcomeDescription: outcome.description || null, outcomeStatus: outcome.status || null,
          novigLeague: market.event && market.event.game ? market.event.game.league : null,
          marketTitle: [matchupOf(market), market.type, outcome.description].filter(Boolean).join(" · "),
        },
      };
      const game = gameFields(market, outcome, true);
      if (game.unmatchable) return unmatchable(base, game.unmatchable);
      return Object.assign(base, game);
    });
  }

  // Every record the mirrored responses describe. Orders and parlay legs
  // dedupe on their native id (a row seen in two pages keeps the later one).
  function normalizeNovig(input) {
    const readAt = input.readAt || null;
    const byId = new Map();
    for (const order of input.orders || []) {
      if (!order || !order.id) continue;
      const record = normalizeOrder(order, readAt);
      byId.set(record.id, record);
    }
    for (const parlay of input.parlays || []) {
      if (!parlay || !parlay.id) continue;
      for (const record of normalizeParlay(parlay, readAt)) byId.set(record.id, record);
    }
    return Array.from(byId.values());
  }

  // ---- response bookkeeping (what novig_content.js holds between responses) -----

  const WATCHED_OPERATIONS = {
    ActivePortfolioOrders_Query: "orders",
    SettledPortfolioOrders_Query: "orders",
    ParlayPortfolioQuery: "parlays",
  };

  // Rows of one mirrored response, or null when it is not a portfolio page.
  function rowsOf(operationName, data) {
    const kind = WATCHED_OPERATIONS[operationName];
    if (!kind || !data || typeof data !== "object") return null;
    const rows = kind === "parlays" ? data.parlay : data[operationName];
    return Array.isArray(rows) ? { kind, rows } : null;
  }

  // Pages are keyed on (operation, where-clause): the app reuses
  // ParlayPortfolioQuery for its active and settled lists with different
  // filters, and a page at offset 0 restarts that list.
  function pageKey(operationName, variables) {
    const where = variables && (variables.where ?? variables.item_where);
    return `${operationName}|${where === undefined ? "" : JSON.stringify(where)}`;
  }

  // Apply one response to the held pages (a Map the caller owns). Returns
  // {orders, parlays, complete} or null when the response is not a portfolio
  // page. `complete` is true only when all three portfolio operations have
  // been seen in this tab AND every held list's last page came back short of
  // its limit — nothing the app lists is beyond what was seen. (An
  // Active-only refetch elsewhere in the app must not count as complete, or
  // the panel would drop the parlay legs it holds.)
  function applyResponse(pages, response) {
    const found = rowsOf(response.operationName, response.data);
    if (!found) return null;
    const variables = response.variables || {};
    const offset = typeof variables.offset === "number" ? variables.offset : 0;
    const limit = typeof variables.limit === "number" ? variables.limit : null;
    const key = pageKey(response.operationName, variables);
    if (offset === 0 || !pages.has(key)) pages.set(key, { kind: found.kind, pages: new Map() });
    pages.get(key).pages.set(offset, { rows: found.rows, limit });
    return collectPages(pages);
  }

  function collectPages(pages) {
    const seenOperations = new Set(Array.from(pages.keys()).map((key) => key.split("|")[0]));
    const out = { orders: [], parlays: [], complete: Object.keys(WATCHED_OPERATIONS).every((name) => seenOperations.has(name)) };
    for (const list of pages.values()) {
      const offsets = Array.from(list.pages.keys()).sort((a, b) => a - b);
      for (const offset of offsets) out[list.kind].push(...list.pages.get(offset).rows);
      const last = list.pages.get(offsets[offsets.length - 1]);
      if (last.limit != null && last.rows.length >= last.limit) out.complete = false;
    }
    return out;
  }

  const api = {
    LEAGUES, MARKET_TYPES, WATCHED_OPERATIONS, APPROX_UNMATCHED, APPROX_PENDING, REASON_NOT_GAME,
    normalizeNovig, normalizeOrder, normalizeParlay, probabilityToAmerican,
    rowsOf, pageKey, applyResponse, collectPages,
  };

  if (typeof module !== "undefined" && module.exports) {
    module.exports = api;
  } else {
    root.UnabatedNovigBets = api;
  }
})(typeof globalThis !== "undefined" ? globalThis : this);
