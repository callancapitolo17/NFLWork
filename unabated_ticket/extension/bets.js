// Bet history for the Unabated Ticket panel: the normalised bet record, the
// Kalshi normaliser, and the matcher that flags a line as already bet, bet on
// the other side, bet in another period of the same market, or on a game you
// already have a position in (#114), and where each matched bet sits on its
// market's axis so the next stake can be sized against it (#130). Pure:
// no DOM, no fetch, no chrome.* — loaded as a plain <script> in panel.html
// after teams.js (exposes globalThis.UnabatedBets) and via require() in
// tests/bets.test.js.
//
// Inputs
//   Bets   normalised records (contract in docs/2026-09-11-issue-114-bet-history-plan.md),
//          produced here by normalizeKalshi() or by a venue's own source.
//   Lines  the ticket / edge-row shape feed.describeLine builds: league,
//          awayTeam, homeTeam, eventStartMs, eventId, betType ("Moneyline" /
//          "Spread" / "Total"), period ("FG", "1H", ...), sideIndex (0 = away
//          or Over, 1 = home or Under), points, rotation, and venueIds — the
//          line's EVENT's venue id map (feed.js noteRungVenueIds), which
//          describeLine hands every row of the event by reference. A row
//          shaped by hand (a captured ticket) has none and still matches
//          through the board rows passed as options.lines.
// Outputs matches per line ({tier, bet, label, position}), per-row annotations for the
//          Edges list, the unmatched list with a reason per bet, the retention
//          prune, the native-id dedupe, and the team-crosswalk rows an id join
//          teaches (#118 step 4; the bets service stores them). Nothing here
//          writes anywhere.
//
// normalizeKalshi is a PORT TARGET: the plan runs the Kalshi normaliser in the
// Python bets service (phase 2). It is written here first so the node tests on
// tests/fixtures/bets/kalshi_fixture.json pin the semantics — side, points,
// price, stake, status — that the Python port must reproduce. Keep it small.
// The service leaves awayKey / homeKey null and the panel fills them with
// resolveTeamKeys() on load, so the team table lives only in teams.js.
//
// Kalshi facts the normaliser relies on (recon 2026-09-11, plan § Recon):
//   event_ticker  <SERIES>-<YYMMMDD>[HHMM]<AWAY><HOME>[G1|G2]; football suffixes
//                 carry the date only (ET), MLB suffixes carry HHMM ET too.
//   event.title   "Missouri St. vs Texas A&M: Spread" / "PIT Steelers vs NE
//                 Patriots" — first team away, second home; the ": Market"
//                 tail is absent on moneyline events.
//   event.sub_title "MOSU vs TXAM (Sep 5)" — the codes in away/home order; a
//                 spread or moneyline strike names its YES team by code.
//   market        strike_type/floor_strike: spread "-TXAM39" -> floor 38.5 ->
//                 YES = TXAM -38.5, NO = other team +38.5; total "-52" ->
//                 floor 51.5 -> YES = Over, NO = Under; moneyline YES = the
//                 team, NO = the other team OR A TIE (approx flag).
//   fills         count_fp, side yes/no, action buy/sell, yes/no_price_dollars,
//                 created_time. positions: signed position_fp (negative = NO),
//                 position_fp 0 with total_traded_dollars > 0 = closed.

(function (root) {
  "use strict";

  const teams = typeof module !== "undefined" && module.exports ? require("./teams.js") : root.UnabatedTeams;

  // related_*: the same market in ANOTHER period of the game (a 1H total on
  // an FG total row). same_game: a different market — never sized off (#129).
  const TIER_RANK = { same_line: 0, same_side: 1, opposite: 2, related_same: 3, related_opposite: 4, same_game: 5 };
  // Spreads and moneylines are cuts on one axis, the margin (away minus
  // home); totals are cuts on the other. Bets on one axis size each other.
  const AXIS_TOTAL = "total";
  const AXIS_MARGIN = "margin";
  const AXIS_OF_BET_TYPE = { total: AXIS_TOTAL, spread: AXIS_MARGIN, moneyline: AXIS_MARGIN };
  const SIDE_AWAY_OR_OVER = 0;
  const SIDE_HOME_OR_UNDER = 1;
  const MONEYLINE_CUT = 0.5;
  const LEAGUE_WITH_THREE_WAY_MONEYLINE = "soccer";
  // Why a matched bet cannot be placed on its axis (#130 guards): named on
  // the card, never guessed around.
  const REASON_PARLAY_LEG = "parlay leg";
  const REASON_NO_STAKE = "no stake on the record";
  const REASON_BET_TYPE = "bet type not sized";
  const REASON_TIE_CAVEAT = "Kalshi NO also wins on a tie";
  const REASON_THREE_WAY = "three-way moneyline";
  const REASON_QUARTER_LINE = "quarter line";
  const REASON_NO_SIDE = "side not resolved";
  const REASON_NO_NUMBER = "no number on the line";
  const REASON_NO_PERIOD = "no period on the record";
  const START_TOLERANCE_MS = 30 * 60 * 1000;
  const DATE_TOLERANCE_DAYS = 1;
  const DAY_MS = 24 * 3600 * 1000;
  const RETENTION_DAYS_DEFAULT = 30;
  const EASTERN = "America/New_York";
  const TIE_CAVEAT = "kalshi_no_side_includes_tie";
  // A venue that gives no game date at all (BetOnline's report, #115): the bet
  // matches an event that starts within this window around its placed time —
  // 12 h back for live bets and a report clock read in the wrong offset, 14
  // days ahead because the account bets NFL totals up to two weeks out.
  const DATE_UNKNOWN = "game_date_unknown";
  const PLACED_WINDOW_BEFORE_MS = 12 * 3600 * 1000;
  const PLACED_WINDOW_AFTER_MS = 14 * DAY_MS;
  // Leagues whose games can end level, so a NO on a team market also wins on a tie.
  const LEAGUES_WITH_TIES = new Set(["nfl", "cfb", "soccer"]);

  // Series whose tickers this module reads. Verified: the five the account
  // traded (fixture) plus KXNFLSPREAD / KXNFLTOTAL / KXNFL1HTOTAL /
  // KXNFL1HSPREAD / KXNCAAF1HSPREAD read off the public markets endpoint on
  // 2026-09-11. MLB series follow kalshi_common/leg_types (bots' tickers);
  // only KXMLBRFI is in the fixture. `fixedPoints` overrides floor_strike.
  const GAME_SERIES = {
    KXNFLGAME: { league: "nfl", betType: "moneyline", period: "FG" },
    KXNFLSPREAD: { league: "nfl", betType: "spread", period: "FG" },
    KXNFLTOTAL: { league: "nfl", betType: "total", period: "FG" },
    KXNFL1HSPREAD: { league: "nfl", betType: "spread", period: "1H" },
    KXNFL1HTOTAL: { league: "nfl", betType: "total", period: "1H" },
    KXNCAAFGAME: { league: "cfb", betType: "moneyline", period: "FG" },
    KXNCAAFSPREAD: { league: "cfb", betType: "spread", period: "FG" },
    KXNCAAFTOTAL: { league: "cfb", betType: "total", period: "FG" },
    KXNCAAF1HSPREAD: { league: "cfb", betType: "spread", period: "1H" },
    KXNCAAF1HTOTAL: { league: "cfb", betType: "total", period: "1H" },
    KXMLBGAME: { league: "mlb", betType: "moneyline", period: "FG" },
    KXMLBSPREAD: { league: "mlb", betType: "spread", period: "FG" },
    KXMLBTOTAL: { league: "mlb", betType: "total", period: "FG" },
    KXMLBF5: { league: "mlb", betType: "moneyline", period: "F5" },
    KXMLBF5SPREAD: { league: "mlb", betType: "spread", period: "F5" },
    KXMLBF5TOTAL: { league: "mlb", betType: "total", period: "F5" },
    KXMLBRFI: { league: "mlb", betType: "total", period: "I1", fixedPoints: 0.5 },
  };
  // Futures, props and the bots' combos seen in the account: shown, never matched.
  const NON_GAME_SERIES = new Set([
    "KXJOINCLUB", "KXWCAWARD", "KXNEXTTEAMNFL", "KXNFLOROTY", "KXSTARTINGQBWEEK1",
    "KXWCGOALLEADER", "KXMVECROSSCATEGORY",
  ]);
  const REASON_NOT_GAME = "not a game market";
  const REASON_UNKNOWN_SERIES = "unknown Kalshi series";
  const REASON_AMBIGUOUS = "ambiguous game";
  const REASON_LEAGUE_OFF = "league not on the scanner";

  // Venue ids the board's rungs carry, in the shapes feed.js accepts (its
  // KALSHI_CONTRACT_RE event group and NOVIG_OUTCOME_ID_RE): a bet id in any
  // other shape can never be on the board, so it is "no id", not a miss.
  const KALSHI_EVENT_SUFFIX_RE = /^[A-Z0-9]+$/;
  const NOVIG_OUTCOME_ID_RE = /^[0-9a-f]{8}-[0-9a-f]{4}-[0-9a-f]{4}-[0-9a-f]{4}-[0-9a-f]{12}$/;
  const JOIN_KALSHI_EVENT = "kalshi_event";
  const JOIN_NOVIG_OUTCOME = "novig_outcome";
  const JOIN_NAME = "name";

  const MONTHS = { JAN: 0, FEB: 1, MAR: 2, APR: 3, MAY: 4, JUN: 5, JUL: 6, AUG: 7, SEP: 8, OCT: 9, NOV: 10, DEC: 11 };
  const EVENT_SUFFIX_RE = /^(\d{2})([A-Z]{3})(\d{2})(\d{4})?([A-Z0-9]+?)(G[12])?$/;
  const EVENT_TITLE_RE = /^(.+?) vs (.+?)(?:: .+)?$/;
  const SUB_TITLE_CODES_RE = /^([A-Z0-9]+) vs ([A-Z0-9]+)/;
  const SPREAD_STRIKE_RE = /^([A-Z]+[A-Z0-9]*?)(\d+)$/;

  // ---- time helpers (Eastern) ------------------------------------------------

  function easternParts(ms) {
    const formatter = new Intl.DateTimeFormat("en-US", {
      timeZone: EASTERN, year: "numeric", month: "2-digit", day: "2-digit",
      hour: "2-digit", minute: "2-digit", hourCycle: "h23",
    });
    const parts = {};
    for (const part of formatter.formatToParts(new Date(ms))) parts[part.type] = part.value;
    return { year: Number(parts.year), month: Number(parts.month) - 1, day: Number(parts.day),
      hour: Number(parts.hour), minute: Number(parts.minute) };
  }

  // "YYYY-MM-DD" of an instant in Eastern time.
  function easternDateOf(ms) {
    const p = easternParts(ms);
    return `${p.year}-${String(p.month + 1).padStart(2, "0")}-${String(p.day).padStart(2, "0")}`;
  }

  // Eastern wall-clock -> UTC ms. The guess is read back in Eastern to find the
  // offset in force at that moment (EDT or EST).
  function easternToUtcMs(year, month, day, hour, minute) {
    const guess = Date.UTC(year, month, day, hour, minute);
    const p = easternParts(guess);
    const asIfUtc = Date.UTC(p.year, p.month, p.day, p.hour, p.minute);
    return guess - (asIfUtc - guess);
  }

  function dateStringToUtcMs(dateString) {
    const [y, m, d] = dateString.split("-").map(Number);
    return Date.UTC(y, m - 1, d);
  }

  // "Sep 10 2:15 PM" in Eastern time, the way the plan's labels read.
  function formatPlacedAt(iso) {
    const ms = Date.parse(iso);
    if (!Number.isFinite(ms)) return "unknown time";
    const formatter = new Intl.DateTimeFormat("en-US", {
      timeZone: EASTERN, month: "short", day: "numeric", hour: "numeric", minute: "2-digit", hour12: true,
    });
    const parts = {};
    for (const part of formatter.formatToParts(new Date(ms))) parts[part.type] = part.value;
    return `${parts.month} ${parts.day} ${parts.hour}:${parts.minute} ${parts.dayPeriod}`;
  }

  // ---- price / money helpers -------------------------------------------------

  function centsToAmerican(cents) {
    if (!(cents > 0 && cents < 100)) return null;
    if (cents >= 50) return -Math.round((cents / (100 - cents)) * 100);
    return Math.round(((100 - cents) / cents) * 100);
  }

  function signedNumber(value) {
    return value > 0 ? `+${value}` : `${value}`;
  }

  // Every price in the panel reads in both worlds: the American number and the
  // prediction-market cents the exchanges actually quote (panel.js fmtPriceBoth,
  // kelly.bookProbOf). A bet record carries only the American price it was
  // placed at, so the cents are derived from it.
  function priceBoth(american) {
    const decimal = american > 0 ? 1 + american / 100 : 1 + 100 / Math.abs(american);
    return `${signedNumber(american)} \u00b7 ${(100 / decimal).toFixed(1)}\u00a2`;
  }

  function roundCents(dollars) {
    return Math.round(dollars * 100) / 100;
  }

  function formatStake(stake) {
    if (stake == null) return "$?";
    return Number.isInteger(stake) ? `$${stake}` : `$${stake.toFixed(2)}`;
  }

  // Venues whose display name is not the capitalised key.
  const VENUE_LABELS = { bfa: "BFA", polymarket_us: "Polymarket US" };

  function venueLabel(venue) {
    if (!venue) return "unknown venue";
    return VENUE_LABELS[venue] || venue.charAt(0).toUpperCase() + venue.slice(1);
  }

  // ---- Kalshi ticker parsing -------------------------------------------------

  function seriesOf(ticker) {
    return ticker.split("-")[0];
  }

  // {eventDate, eventStart, gameNumber} from the event-ticker suffix, or null.
  function parseEventSuffix(eventTicker) {
    const suffix = eventTicker.slice(eventTicker.indexOf("-") + 1);
    const m = EVENT_SUFFIX_RE.exec(suffix);
    if (!m || !(m[2] in MONTHS)) return null;
    const year = 2000 + Number(m[1]);
    const month = MONTHS[m[2]];
    const day = Number(m[3]);
    const eventDate = `${year}-${String(month + 1).padStart(2, "0")}-${String(day).padStart(2, "0")}`;
    let eventStart = null;
    if (m[4]) {
      const hour = Number(m[4].slice(0, 2));
      const minute = Number(m[4].slice(2));
      eventStart = new Date(easternToUtcMs(year, month, day, hour, minute)).toISOString();
    }
    return { eventDate, eventStart, gameNumber: m[6] ? Number(m[6][1]) : null };
  }

  // {awayTeam, homeTeam, codes: [awayCode, homeCode]} from the public event payload.
  function parseEventTeams(event) {
    const title = EVENT_TITLE_RE.exec(event.title || "");
    const codes = SUB_TITLE_CODES_RE.exec(event.sub_title || "");
    if (!title || !codes) return null;
    return { awayTeam: title[1], homeTeam: title[2], codes: [codes[1], codes[2]] };
  }

  // Which side (0 away / 1 home) a spread or moneyline strike names as YES, and
  // the strike's own number. Returns {yesIndex, points} or an error string.
  function parseStrike(marketTicker, market, spec) {
    const strike = marketTicker.slice(market.event_ticker.length + 1);
    const codes = spec.codes;
    if (spec.betType === "total") {
      const points = spec.fixedPoints ?? market.floor_strike;
      return typeof points === "number" ? { yesIndex: null, points } : `no floor_strike on ${marketTicker}`;
    }
    if (spec.betType === "moneyline") {
      const yesIndex = codes.indexOf(strike);
      return yesIndex >= 0 ? { yesIndex, points: null } : `strike ${strike} not in ${codes.join("/")}`;
    }
    const m = SPREAD_STRIKE_RE.exec(strike);
    const yesIndex = m ? codes.indexOf(m[1]) : -1;
    if (yesIndex < 0) return `strike ${strike} not in ${codes.join("/")}`;
    if (typeof market.floor_strike !== "number") return `no floor_strike on ${marketTicker}`;
    return { yesIndex, points: market.floor_strike };
  }

  // Side + points for one Kalshi (market, yes|no) in the record's convention:
  // points is the side's own number.
  function sideOfContract(spec, strike, contractSide, league) {
    const approx = [];
    if (spec.betType === "total") {
      return { side: contractSide === "yes" ? "over" : "under", points: strike.points, approx };
    }
    const teamIndex = contractSide === "yes" ? strike.yesIndex : 1 - strike.yesIndex;
    const side = teamIndex === 0 ? "away" : "home";
    if (spec.betType === "moneyline") {
      if (contractSide === "no" && LEAGUES_WITH_TIES.has(league)) approx.push(TIE_CAVEAT);
      return { side, points: null, approx };
    }
    const points = contractSide === "yes" ? -strike.points : strike.points;
    return { side, points, approx };
  }

  // ---- Kalshi fills -> records -------------------------------------------------

  function toNumber(value) {
    const n = typeof value === "number" ? value : parseFloat(value);
    return Number.isFinite(n) ? n : 0;
  }

  function fillPriceCents(fill) {
    return toNumber(fill.side === "yes" ? fill.yes_price_dollars : fill.no_price_dollars) * 100;
  }

  // VWAP in cents over buys (sells only when there were no buys), net contracts,
  // and the first fill time, for one (ticker, side).
  function aggregateFills(fills) {
    const sorted = fills.slice().sort((a, b) => Date.parse(a.created_time) - Date.parse(b.created_time));
    let bought = 0, boughtCost = 0, sold = 0, soldCost = 0;
    for (const fill of sorted) {
      const count = toNumber(fill.count_fp);
      const cents = fillPriceCents(fill);
      if (fill.action === "sell") { sold += count; soldCost += count * cents; }
      else { bought += count; boughtCost += count * cents; }
    }
    const vwapCents = bought > 0 ? boughtCost / bought : sold > 0 ? soldCost / sold : null;
    return { vwapCents, netContracts: bought - sold, firstFillAt: sorted[0].created_time, fillCount: sorted.length };
  }

  function statusOf(contractSide, market, position, netContracts) {
    const result = market ? market.result : null;
    if (position) {
      const positionFp = toNumber(position.position_fp);
      if (positionFp === 0 && toNumber(position.total_traded_dollars) > 0) return "closed";
      const openSide = positionFp > 0 ? "yes" : positionFp < 0 ? "no" : null;
      if (openSide !== contractSide) return "closed";
    } else if (netContracts <= 0) {
      return "closed";
    }
    if (result === "yes" || result === "no") return result === contractSide ? "won" : "lost";
    if (result === "void") return "void";
    if (result) return "unknown";
    return position || netContracts > 0 ? "open" : "unknown";
  }

  function closedAtOf(status, position, market, lastFillAt) {
    if (status === "open") return null;
    if (status === "closed") return position ? position.last_updated_ts ?? lastFillAt : lastFillAt;
    return (market && market.expected_expiration_time) || lastFillAt;
  }

  function unmatchableRecord(base, reason) {
    return Object.assign(base, {
      league: null, eventStart: null, eventDate: null, awayTeam: null, homeTeam: null, awayKey: null, homeKey: null,
      betType: "other", period: null, side: null, points: null, approx: [], unmatchable: reason,
    });
  }

  // One record per (ticker, side) seen in fills. `markets`/`events` are keyed by
  // ticker / event_ticker (the service's cached public GETs); a game-series
  // market missing either fails closed with a specific reason.
  function normalizeKalshi(input) {
    const fills = input.fills || [];
    const positions = new Map((input.positions || []).map((p) => [p.ticker, p]));
    const markets = input.markets || {};
    const events = input.events || {};
    const fetchedAt = input.fetchedAt ?? null;
    const groups = new Map();
    for (const fill of fills) {
      const key = `${fill.ticker}:${fill.side}`;
      if (!groups.has(key)) groups.set(key, []);
      groups.get(key).push(fill);
    }
    const records = [];
    for (const [key, groupFills] of groups) {
      const ticker = groupFills[0].ticker;
      const contractSide = groupFills[0].side;
      const series = seriesOf(ticker);
      const market = markets[ticker] || null;
      const position = positions.get(ticker) || null;
      const agg = aggregateFills(groupFills);
      const status = statusOf(contractSide, market, position, agg.netContracts);
      const contracts = position && status === "open" ? Math.abs(toNumber(position.position_fp)) : Math.max(agg.netContracts, 0);
      const entryDollars = agg.vwapCents == null ? null : agg.vwapCents / 100;
      const lastFillAt = groupFills.map((f) => f.created_time).sort().pop();
      const base = {
        id: `kalshi:${key}`,
        source: "kalshi_api",
        venue: "kalshi",
        rotation: null,
        price: agg.vwapCents == null ? null : centsToAmerican(agg.vwapCents),
        stake: entryDollars == null ? null : roundCents(contracts * entryDollars),
        toWin: entryDollars == null ? null : roundCents(contracts * (1 - entryDollars)),
        contracts,
        placedAt: agg.firstFillAt,
        status,
        closedAt: closedAtOf(status, position, market, lastFillAt),
        isParlayLeg: false, parlayId: null, legIndex: null, legCount: null,
        sourceFetchedAt: fetchedAt,
        // Kalshi's own ids (#118): the event ticker's suffix is what
        // Unabated's Kalshi alt rungs carry inside sourceKey.
        venueIds: { marketTicker: ticker, eventTicker: market ? market.event_ticker ?? null : null },
        raw: {
          ticker, series, contractSide, fillCount: agg.fillCount, vwapCents: agg.vwapCents,
          positionFp: position ? toNumber(position.position_fp) : null,
          totalTradedDollars: position ? toNumber(position.total_traded_dollars) : null,
          marketTitle: market ? market.title : null, marketStatus: market ? market.status : null,
          marketResult: market ? market.result : null,
          eventTitle: market && events[market.event_ticker] ? events[market.event_ticker].title : null,
        },
      };
      const spec = GAME_SERIES[series];
      if (!spec) {
        records.push(unmatchableRecord(base, NON_GAME_SERIES.has(series) ? REASON_NOT_GAME : REASON_UNKNOWN_SERIES));
        continue;
      }
      const event = market ? events[market.event_ticker] : null;
      if (!market || !event) {
        records.push(unmatchableRecord(base, `unreadable Kalshi market (no ${market ? "event" : "market"} payload for ${ticker})`));
        continue;
      }
      const suffix = parseEventSuffix(market.event_ticker);
      const eventTeams = parseEventTeams(event);
      if (!suffix || !eventTeams) {
        records.push(unmatchableRecord(base, `unreadable Kalshi market (${suffix ? "event title" : "event suffix"} on ${market.event_ticker})`));
        continue;
      }
      const strike = parseStrike(ticker, market, { betType: spec.betType, fixedPoints: spec.fixedPoints, codes: eventTeams.codes });
      if (typeof strike === "string") {
        records.push(unmatchableRecord(base, `unreadable Kalshi market (${strike})`));
        continue;
      }
      const contract = sideOfContract(spec, strike, contractSide, spec.league);
      records.push(Object.assign(base, {
        league: spec.league,
        eventStart: suffix.eventStart,
        eventDate: suffix.eventDate,
        awayTeam: eventTeams.awayTeam,
        homeTeam: eventTeams.homeTeam,
        awayKey: teams.teamKey(spec.league, eventTeams.awayTeam),
        homeKey: teams.teamKey(spec.league, eventTeams.homeTeam),
        betType: spec.betType,
        period: spec.period,
        side: contract.side,
        points: contract.points,
        approx: contract.approx,
        unmatchable: null,
      }));
    }
    return records;
  }

  // ---- game match ------------------------------------------------------------

  // A board line carries Unabated's team ids, which ARE the keys; a line
  // shaped by hand (tests, a captured ticket) resolves by name.
  function lineTeamKeys(line) {
    return {
      away: line.awayTeamId != null ? teams.keyOf(line.league, line.awayTeamId) : teams.teamKey(line.league, line.awayTeam),
      home: line.homeTeamId != null ? teams.keyOf(line.league, line.homeTeamId) : teams.teamKey(line.league, line.homeTeam),
    };
  }

  function sameTeamPair(bet, keys) {
    if (!bet.awayKey || !bet.homeKey || !keys.away || !keys.home) return false;
    return (bet.awayKey === keys.away && bet.homeKey === keys.home) || (bet.awayKey === keys.home && bet.homeKey === keys.away);
  }

  function rotationMatches(bet, line) {
    if (bet.rotation == null) return false;
    return [line.rotation, line.awayRotation, line.homeRotation].some((r) => r != null && r === bet.rotation);
  }

  function timeMatches(bet, line) {
    if (typeof line.eventStartMs !== "number") return false;
    if (bet.eventStart) {
      const betMs = Date.parse(bet.eventStart);
      return Number.isFinite(betMs) && Math.abs(betMs - line.eventStartMs) <= START_TOLERANCE_MS;
    }
    if (bet.eventDate) {
      const delta = Math.abs(dateStringToUtcMs(easternDateOf(line.eventStartMs)) - dateStringToUtcMs(bet.eventDate));
      return delta <= DATE_TOLERANCE_DAYS * DAY_MS;
    }
    if (!(bet.approx && bet.approx.includes(DATE_UNKNOWN)) || !bet.placedAt) return false;
    const placedMs = Date.parse(bet.placedAt);
    return Number.isFinite(placedMs)
      && line.eventStartMs >= placedMs - PLACED_WINDOW_BEFORE_MS
      && line.eventStartMs <= placedMs + PLACED_WINDOW_AFTER_MS;
  }

  // A rotation match must put every team the bet names (BetOnline names only
  // its own team on a spread, both on a total) in the row's game: the same
  // rotation number comes round again the next week, and a dateless bet's
  // 14-day window would otherwise accept both weeks as "ambiguous game".
  // Decidable only when both of the row's teams resolve; otherwise the
  // rotation stands on its own.
  function knownTeamFits(bet, keys) {
    if (!keys.away || !keys.home) return true;
    const known = [bet.awayKey, bet.homeKey].filter(Boolean);
    return known.every((key) => key === keys.away || key === keys.home);
  }

  function gameMatches(bet, line) {
    if (bet.unmatchable || bet.league !== line.league) return false;
    const keys = lineTeamKeys(line);
    if (!sameTeamPair(bet, keys) && !(rotationMatches(bet, line) && knownTeamFits(bet, keys))) return false;
    return timeMatches(bet, line);
  }

  function isMatchable(bet) {
    return bet.status === "open" && !bet.unmatchable;
  }

  function eventIdentity(line) {
    if (line.eventId != null) return `id:${line.eventId}`;
    const keys = lineTeamKeys(line);
    return `${line.league}|${[keys.away, keys.home].sort().join("/")}|${line.eventStartMs}`;
  }

  // Distinct board events a bet's game rule accepts; two means ambiguous.
  function candidateEvents(bet, lines) {
    const found = new Set();
    for (const line of lines) if (gameMatches(bet, line)) found.add(eventIdentity(line));
    return found;
  }

  // ---- venue id join (#118 step 3) --------------------------------------------
  //
  // Before the name rule, a bet whose venue ids the board's rungs carry joins
  // its game EXACTLY: a Kalshi bet on its event-ticker suffix, a Novig bet on
  // its outcome id. No team name has to resolve. Exact or nothing: an id on
  // no board event falls through to the name rule, an id on two board events
  // is ambiguous (never a guess), and an id join is final — the name rule is
  // not consulted, so a team pair that points at another event cannot win.
  //
  // The join decides WHICH GAME; the tier still comes from the bet's own
  // betType / period / side / points against the row's current line
  // (tierOf). The id map's lineKey / mainKey are not read: the map is as old
  // as the last snapshot while the changes stream moves main lines, and a
  // Novig lay's outcome id names the side the bet is AGAINST. Only the
  // contract's fixed strike and side are read, to orient a spread bet whose
  // team names do not resolve (sideIndexByVenueId).

  // "26SEP19DUQWSU" from "KXNCAAFSPREAD-26SEP19DUQWSU": everything after the
  // first "-" (series names carry none), kept whole — Kalshi's CFB team codes
  // are not Unabated's abbreviations, so the suffix is never split.
  function kalshiEventSuffixOfTicker(eventTicker) {
    if (typeof eventTicker !== "string") return null;
    const dash = eventTicker.indexOf("-");
    if (dash <= 0) return null;
    const suffix = eventTicker.slice(dash + 1);
    return KALSHI_EVENT_SUFFIX_RE.test(suffix) ? suffix : null;
  }

  // The bet's venue id the board can carry: {join, id}, or null. Novig
  // moneylines never have one on the board (moneylines have no rungs).
  function betVenueId(bet) {
    const ids = bet.venueIds;
    if (!ids || typeof ids !== "object") return null;
    if (bet.venue === "kalshi") {
      const suffix = kalshiEventSuffixOfTicker(ids.eventTicker);
      return suffix ? { join: JOIN_KALSHI_EVENT, id: suffix } : null;
    }
    if (bet.venue === "novig" && bet.betType !== "moneyline"
      && typeof ids.outcomeId === "string" && NOVIG_OUTCOME_ID_RE.test(ids.outcomeId)) {
      return { join: JOIN_NOVIG_OUTCOME, id: ids.outcomeId };
    }
    return null;
  }

  function addIdToIndex(byId, id, identity, league) {
    if (typeof id !== "string") return;
    if (!byId.has(id)) byId.set(id, new Map());
    byId.get(id).set(identity, league);
  }

  // id -> Map(event identity -> league), over the board rows' venueIds. Rows
  // of one event share one map, so each event is read once.
  function venueIdIndex(lines) {
    const index = { [JOIN_KALSHI_EVENT]: new Map(), [JOIN_NOVIG_OUTCOME]: new Map() };
    const seen = new Set();
    for (const line of lines) {
      const ids = line.venueIds;
      if (!ids || typeof ids !== "object" || line.eventId == null) continue;
      const identity = eventIdentity(line);
      if (seen.has(identity)) continue;
      seen.add(identity);
      const suffixes = Array.isArray(ids.kalshiEventSuffixes) ? ids.kalshiEventSuffixes : [];
      for (const suffix of suffixes) addIdToIndex(index[JOIN_KALSHI_EVENT], suffix, identity, line.league);
      const outcomes = ids.novigOutcomes && typeof ids.novigOutcomes === "object" ? Object.keys(ids.novigOutcomes) : [];
      for (const outcomeId of outcomes) addIdToIndex(index[JOIN_NOVIG_OUTCOME], outcomeId, identity, line.league);
    }
    return index;
  }

  // A board: its rows, their venue id index, and each bet's game decided once.
  function boardOf(lines) {
    return { lines, venueIds: venueIdIndex(lines), games: new Map() };
  }

  // {join, venueId, idInOtherLeague, events}: the board event identities the bet's game is —
  // by venue id when its id is on an event of the bet's own league, else by
  // the name rule. One event = matched, several = ambiguous, none = a miss.
  function resolveGame(bet, board) {
    const venueId = betVenueId(bet);
    const carriers = venueId ? board.venueIds[venueId.join].get(venueId.id) : null;
    const idEvents = new Set();
    let idInOtherLeague = false;
    for (const [identity, league] of carriers || []) {
      if (league === bet.league) idEvents.add(identity);
      else idInOtherLeague = true;
    }
    if (idEvents.size > 0) return { join: venueId.join, venueId, events: idEvents };
    return { join: JOIN_NAME, venueId, idInOtherLeague, events: candidateEvents(bet, board.lines) };
  }

  function gameOf(bet, board) {
    if (!board.games.has(bet)) board.games.set(bet, resolveGame(bet, board));
    return board.games.get(bet);
  }

  // Is this line on the bet's game? An id join names the event outright; the
  // name rule is re-checked on the line itself, which a captured ticket's
  // hand-shaped row may not share with the board.
  // A line with no eventId (a hand-shaped row) has no identity the board's
  // id join could name, so it falls back to the name rule.
  function lineInGame(bet, line, game) {
    if (game.join === JOIN_NAME || line.eventId == null) return gameMatches(bet, line);
    return bet.league === line.league && game.events.has(eventIdentity(line));
  }

  function venueIdLabel(venueId) {
    return venueId.join === JOIN_KALSHI_EVENT ? `Kalshi event ${venueId.id}` : "Novig outcome";
  }

  function ambiguousReason(game) {
    if (game.join === JOIN_NAME) return REASON_AMBIGUOUS;
    return `${REASON_AMBIGUOUS} (${venueIdLabel(game.venueId)} on ${game.events.size} board events)`;
  }

  // ---- tiers -----------------------------------------------------------------

  function lineBetType(line) {
    return typeof line.betType === "string" ? line.betType.toLowerCase() : null;
  }

  // "over"/"under" or the team key the line's side names.
  function lineSideKey(line) {
    if (lineBetType(line) === "total") return line.sideIndex === 0 ? "over" : "under";
    const keys = lineTeamKeys(line);
    return line.sideIndex === 0 ? keys.away : keys.home;
  }

  // Team bets only: a total's side is read straight off the record.
  function betSideKey(bet) {
    return bet.side === "away" ? bet.awayKey : bet.homeKey;
  }

  function betOtherSideKey(bet) {
    return bet.side === "away" ? bet.homeKey : bet.awayKey;
  }

  function samePoints(a, b) {
    return (a == null && b == null) || (typeof a === "number" && a === b);
  }

  // The side (0 away or Over / 1 home or Under, Unabated's frame) a bet sits
  // on in the row's game, or null when nothing says which. A team bet reads
  // its team keys, then its venue contract, then its rotation — never the
  // record's own away/home, which is the VENUE's frame.
  function betSideIndexOn(bet, line) {
    if (bet.betType === "total") return bet.side === "over" ? SIDE_AWAY_OR_OVER : bet.side === "under" ? SIDE_HOME_OR_UNDER : null;
    return sideIndexByTeamKeys(bet, line) ?? sideIndexByVenueId(bet, line) ?? sideIndexByRotation(bet, line);
  }

  // How a bet relates to a row. Same axis, same period: same_line (type and
  // number too), same_side, opposite — a moneyline held against a spread on
  // the same team is the same side (#130). Same axis, another period:
  // related_same / related_opposite. Another axis: same_game.
  function tierOf(bet, line) {
    const lineType = lineBetType(line);
    const axis = AXIS_OF_BET_TYPE[bet.betType];
    if (axis == null || axis !== AXIS_OF_BET_TYPE[lineType]) return "same_game";
    if (lineSideKey(line) == null) return "same_game";
    const betSideIndex = betSideIndexOn(bet, line);
    if (betSideIndex == null) return "same_game";
    const sameDirection = betSideIndex === line.sideIndex;
    if (bet.period !== line.period) return sameDirection ? "related_same" : "related_opposite";
    if (!sameDirection) return "opposite";
    return bet.betType === lineType && samePoints(bet.points, line.points) ? "same_line" : "same_side";
  }

  // The side (0 away / 1 home) the bet's team keys put it on in the row's
  // game: its own team's key where the game has it, else the other team's
  // key on the opposite side. One key is enough — an id join reaches here
  // with names teams.js half-resolved (#118 step 3). A key the game does not
  // have names no side (a name an id join overruled).
  function sideIndexByTeamKeys(bet, line) {
    const keys = lineTeamKeys(line);
    const own = betSideKey(bet);
    const other = betOtherSideKey(bet);
    if (own != null && own === keys.away) return 0;
    if (own != null && own === keys.home) return 1;
    if (other != null && other === keys.away) return 1;
    if (other != null && other === keys.home) return 0;
    return null;
  }

  // BetOnline names its team as it likes: the bet's rotation IS its team, so
  // the row's away/home rotations say which side it sits on.
  function sideIndexByRotation(bet, line) {
    if (bet.rotation == null) return null;
    return bet.rotation === line.awayRotation ? 0 : bet.rotation === line.homeRotation ? 1 : null;
  }

  // The bet's own contract on the row's event id map (#118 step 3): a Kalshi
  // market ticker (either contract side, "Y-"/"N-") or a Novig outcome id.
  // Its target carries the contract's fixed strike and Unabated side, so the
  // bet sits on that side at the same number and on the other side at the
  // negated one — a Kalshi NO or a Novig lay. Spreads only (a total needs no
  // side, a moneyline has no rung); a 0 strike cannot tell the sides apart.
  // A Kalshi moneyline whose names do not resolve stays same_game.
  function sideIndexByVenueId(bet, line) {
    if (bet.betType !== "spread" || typeof bet.points !== "number" || bet.points === 0) return null;
    const ids = bet.venueIds;
    const map = line.venueIds;
    if (!ids || typeof ids !== "object" || !map || typeof map !== "object") return null;
    const targets = [];
    if (bet.venue === "kalshi" && typeof ids.marketTicker === "string" && map.kalshiContracts) {
      targets.push(map.kalshiContracts[`Y-${ids.marketTicker}`], map.kalshiContracts[`N-${ids.marketTicker}`]);
    } else if (bet.venue === "novig" && typeof ids.outcomeId === "string" && map.novigOutcomes) {
      targets.push(map.novigOutcomes[ids.outcomeId]);
    }
    for (const target of targets) {
      if (!target || (target.sideIndex !== 0 && target.sideIndex !== 1)) continue;
      if (target.points === bet.points) return target.sideIndex;
      if (target.points === -bet.points) return 1 - target.sideIndex;
    }
    return null;
  }

  // ---- labels ----------------------------------------------------------------

  function betTeamName(bet, side) {
    return side === "away" ? bet.awayTeam : bet.homeTeam;
  }

  // "Chattanooga -5.5 +138", "1H Under 22.5 -104", "NO PIT Steelers ≈ NE Patriots or tie -178".
  function describeBet(bet) {
    const period = bet.period && bet.period !== "FG" ? `${bet.period} ` : "";
    let pick;
    if (bet.betType === "total") {
      pick = `${bet.side === "over" ? "Over" : "Under"} ${bet.points}`;
    } else if (bet.betType === "spread") {
      pick = `${betTeamName(bet, bet.side)} ${signedNumber(bet.points)}`;
    } else if (bet.betType === "other") {
      pick = (bet.raw && bet.raw.marketTitle) || bet.id;
    } else if (bet.approx && bet.approx.includes(TIE_CAVEAT)) {
      const yesTeam = betTeamName(bet, bet.side === "away" ? "home" : "away");
      pick = `NO ${yesTeam} ≈ ${betTeamName(bet, bet.side)} or tie`;
    } else {
      pick = betTeamName(bet, bet.side);
    }
    const price = bet.price == null ? "" : ` ${priceBoth(bet.price)}`;
    const parlay = bet.isParlayLeg ? " (parlay leg)" : "";
    return `${period}${pick}${price}${parlay}`;
  }

  // The other side of a spread carries the negated number; a total or a
  // moneyline carries the same one.
  function oppositeSameNumber(bet, line) {
    if (bet.betType === "spread") return typeof bet.points === "number" && bet.points === -line.points;
    return samePoints(bet.points, line.points);
  }

  function linePointsLabel(line) {
    return lineBetType(line) === "spread" ? signedNumber(line.points) : `${line.points}`;
  }

  // How a bet relates to the line it matched. The words are a TAG the panel
  // renders beside the label, not a sentence in front of it: the block has a
  // heading and a colour, so "You are on the OTHER side:" was three quarters
  // of the line saying what a red chip says (user decision 2026-09-14).
  const TIER_LABELS = {
    same_line: "this line",
    same_side: "same side",
    opposite: "other side",
    // The bet's own text starts with its period ("1H Over 20.5"), so the tag
    // only has to name the direction.
    related_same: "same side",
    related_opposite: "other side",
    same_game: "game",
  };

  function tierLabel(tier) {
    return TIER_LABELS[tier] || tier;
  }

  // The line's number as the BET's side would write it. The row's number is in
  // the row's frame, and on a spread the two sides are negatives of each other
  // — printing "now +14" beside a bet on "Chicago -13.5" reads as a 27.5-point
  // move in the wrong direction. Totals do not flip.
  function linePointsAsBetSide(tier, line) {
    if (tier !== "opposite" || lineBetType(line) !== "spread") return linePointsLabel(line);
    return signedNumber(-line.points);
  }

  // The bet itself: what, how much, where — and where the line sits now when
  // it has moved off the bet's number, which is the whole point of those two
  // tiers. No placed-at: it never told you which bet was which (user decision
  // 2026-09-14).
  function labelOf(tier, bet, line) {
    const parts = [describeBet(bet), formatStake(bet.stake), venueLabel(bet.venue)];
    // Only a bet of the row's own type has a number the row can have moved
    // off: a moneyline held against a spread row does not.
    const sameType = bet.betType === lineBetType(line);
    const numberMoved = sameType && (tier === "same_side" || (tier === "opposite" && !oppositeSameNumber(bet, line)));
    if (numberMoved) parts.push(`now ${linePointsAsBetSide(tier, line)}`);
    return parts.join(" · ");
  }

  // ---- positions on the market axis (#130) -----------------------------------
  //
  // Conditional Kelly needs to know WHEN each bet wins. Every spread,
  // moneyline and alt number is a cut on the margin (away score minus home
  // score); every total is a cut on the total. A position is {axis, cut,
  // direction}: the bet wins when the result lands "above" or "below" the
  // cut. Over 52.5 = above 52.5; an away bet at points a = above -a; a home
  // bet at points h = below h; a moneyline = above +0.5 (away) / below -0.5
  // (home). A whole-number cut pushes on the number (condkelly.js).

  // {axis, cut, direction} for a bet type, side and number, or {reason}.
  function axisPosition(betType, sideIndex, points) {
    const direction = sideIndex === SIDE_AWAY_OR_OVER ? "above" : "below";
    if (betType === "moneyline") return { axis: AXIS_MARGIN, cut: sideIndex === SIDE_AWAY_OR_OVER ? MONEYLINE_CUT : -MONEYLINE_CUT, direction };
    if (typeof points !== "number" || !Number.isFinite(points)) return { reason: REASON_NO_NUMBER };
    if (!Number.isInteger(points * 2)) return { reason: REASON_QUARTER_LINE };
    if (betType === "total") return { axis: AXIS_TOTAL, cut: points, direction };
    // `|| 0` folds a pick'em's -0 into 0.
    return { axis: AXIS_MARGIN, cut: (sideIndex === SIDE_AWAY_OR_OVER ? -points : points) || 0, direction };
  }

  // The period names the group a bet is sized in; without one it is left out.
  function hasPeriodName(period) {
    return typeof period === "string" && period !== "";
  }

  // The row's own position: {axis, period, cut, direction} or {reason}.
  function linePosition(line) {
    const betType = lineBetType(line);
    if (AXIS_OF_BET_TYPE[betType] == null) return { reason: REASON_BET_TYPE };
    if (betType === "moneyline" && line.league === LEAGUE_WITH_THREE_WAY_MONEYLINE) return { reason: REASON_THREE_WAY };
    if (line.sideIndex !== SIDE_AWAY_OR_OVER && line.sideIndex !== SIDE_HOME_OR_UNDER) return { reason: REASON_NO_SIDE };
    if (!hasPeriodName(line.period)) return { reason: REASON_NO_PERIOD };
    const position = axisPosition(betType, line.sideIndex, line.points);
    return position.reason ? position : { ...position, period: line.period };
  }

  // A matched bet's position in the row's game, with its dollars: {axis,
  // period, cut, direction, stake, toWin}, or {reason} for a bet that is
  // left out of the sizing and named.
  function positionOf(bet, line) {
    if (AXIS_OF_BET_TYPE[bet.betType] == null) return { reason: REASON_BET_TYPE };
    if (bet.isParlayLeg) return { reason: REASON_PARLAY_LEG };
    if (!hasPeriodName(bet.period)) return { reason: REASON_NO_PERIOD };
    if (!(typeof bet.stake === "number" && bet.stake > 0) || !(typeof bet.toWin === "number" && bet.toWin > 0)) return { reason: REASON_NO_STAKE };
    if (bet.betType === "moneyline") {
      if (Array.isArray(bet.approx) && bet.approx.includes(TIE_CAVEAT)) return { reason: REASON_TIE_CAVEAT };
      if (bet.league === LEAGUE_WITH_THREE_WAY_MONEYLINE) return { reason: REASON_THREE_WAY };
    }
    const sideIndex = betSideIndexOn(bet, line);
    if (sideIndex == null) return { reason: REASON_NO_SIDE };
    const position = axisPosition(bet.betType, sideIndex, bet.points);
    return position.reason ? position : { ...position, period: bet.period, stake: bet.stake, toWin: bet.toWin };
  }

  // ---- public API ------------------------------------------------------------

  // Matches for one line, strongest tier first. options.lines is the whole
  // board (defaults to [line]) — the rows whose venueIds a bet's id joins on,
  // and where a bet whose game rule accepts two distinct board events is
  // ambiguous and lands in `unmatched`, never in `matches`.
  // Closed and settled bets never match and are not listed here.
  function matchBets(line, bets, options) {
    const lines = options && Array.isArray(options.lines) ? options.lines : [line];
    return matchOnBoard(line, bets, boardOf(lines));
  }

  function matchOnBoard(line, bets, board) {
    const matches = [];
    const unmatched = [];
    for (const bet of bets) {
      if (!isMatchable(bet)) continue;
      const game = gameOf(bet, board);
      if (!lineInGame(bet, line, game)) continue;
      if (game.events.size > 1) {
        unmatched.push({ bet, reason: ambiguousReason(game) });
        continue;
      }
      const tier = tierOf(bet, line);
      matches.push({ tier, bet, label: labelOf(tier, bet, line), position: positionOf(bet, line) });
    }
    matches.sort((a, b) => TIER_RANK[a.tier] - TIER_RANK[b.tier] || (b.bet.stake ?? 0) - (a.bet.stake ?? 0));
    return { matches, unmatched };
  }

  // Per Edges row: the strongest tier (or null) and its matches. A bet you
  // hold never hides a line — it changes the size of the next one (see
  // betsview.stakeAdvice, which also counts the dollars held and against).
  // options.lines is the whole board (defaults to rows), the same one the
  // Ticket banner and the unmatched list decide games on: a bet's id event
  // may have no listed edge row while its team names fit one that is listed.
  function annotateRows(rows, bets, options) {
    const board = boardOf(options && Array.isArray(options.lines) ? options.lines : rows);
    return rows.map((row) => {
      const { matches } = matchOnBoard(row, bets, board);
      return { tier: matches.length ? matches[0].tier : null, matches };
    });
  }

  function unresolvedTeamNames(bet) {
    const names = [];
    if (bet.awayTeam != null && !bet.awayKey) names.push(bet.awayTeam);
    if (bet.homeTeam != null && !bet.homeKey) names.push(bet.homeTeam);
    return names;
  }

  // Why the name rule found no event for a bet.
  function nameMissReason(bet, leaguesOnBoard) {
    const unresolved = unresolvedTeamNames(bet);
    if (unresolved.length) return `team not recognised (${unresolved.join(", ")})`;
    if (!leaguesOnBoard.has(bet.league)) return REASON_LEAGUE_OFF;
    return "no event on the board yet";
  }

  // Every OPEN bet that matches no line on the board, with why. Closed and
  // settled bets are not problems, so they are not listed. A bet that carries
  // a venue id says both tiers failed — "by id: Kalshi event 26SEP19DUQWSU
  // not on any board ladder; by name: team not recognised (…)" — since an id
  // misses when the venue lists no ladder for the game (Kalshi rungs sat on
  // 126 of 346 CFB events, 2026-09-12) while the name rule may still see it.
  function unmatchedReasons(bets, lines) {
    const leaguesOnBoard = new Set(lines.map((line) => line.league));
    const board = boardOf(lines);
    const out = [];
    for (const bet of bets) {
      if (bet.status !== "open") continue;
      if (bet.unmatchable) { out.push({ bet, reason: bet.unmatchable }); continue; }
      // Matched (by id, team pair or rotation) is not a problem, whatever the
      // team table makes of the names; the diagnoses below explain a miss.
      const game = gameOf(bet, board);
      if (game.events.size === 1) continue;
      if (game.events.size > 1) { out.push({ bet, reason: ambiguousReason(game) }); continue; }
      const nameReason = nameMissReason(bet, leaguesOnBoard);
      if (!game.venueId || nameReason === REASON_LEAGUE_OFF) { out.push({ bet, reason: nameReason }); continue; }
      const idMiss = game.idInOtherLeague ? "only on another league's board ladder" : "not on any board ladder";
      out.push({ bet, reason: `by id: ${venueIdLabel(game.venueId)} ${idMiss}; by name: ${nameReason}` });
    }
    return out;
  }

  // Open bets, plus closed/settled ones whose closedAt is within the window.
  // A non-open bet with no closedAt is kept: dropping it would be silent.
  function pruneForRetention(bets, now, days) {
    const windowMs = (typeof days === "number" ? days : RETENTION_DAYS_DEFAULT) * DAY_MS;
    return bets.filter((bet) => {
      if (bet.status === "open") return true;
      const closedMs = bet.closedAt ? Date.parse(bet.closedAt) : NaN;
      return !Number.isFinite(closedMs) || now - closedMs <= windowMs;
    });
  }

  // ---- team crosswalk (#118 step 4) -------------------------------------------
  //
  // Every unambiguous id join is a lesson: the bet's venue names both teams
  // (Novig by team id, Kalshi by its event-title name) and the joined board
  // event carries both Unabated team ids. Remembered as (venue, league, venue
  // team) -> Unabated team id, the lesson keys a LATER bet of that venue on
  // either team before teams.js sees the name — a bet on a game the venue
  // lists no ladder for, where the id has nothing to join, or a Kalshi
  // moneyline whose names resolve nowhere (no rung carries its contract, so
  // without a key it tiered same_game and held $0). The bets service owns the
  // table (bets.duckdb::team_crosswalk, served with /bets.json); this module
  // only decides what to learn (learnCrosswalk) and applies what was learned
  // (resolveTeamKeys). Fail-closed: learn only from a join on exactly ONE
  // board event whose row carries both Unabated ids, from a bet naming both
  // venue teams; never learn a venue team whose name already resolves to a
  // DIFFERENT team (a swapped away/home — a neutral site — would otherwise
  // write a row that then needs opponent and time to agree before it could
  // mismatch, so it mostly yields no match; better not written at all); never
  // relearn a held venue team as another id. Conflicts are reported, not
  // written.

  // The venue's own team on one side of the bet, {key, name}, or null. Novig
  // records carry the venue's team id (awayTeamVenue.id); Kalshi has no team
  // ids, so its event-title name is the key — stable per team ("PIT
  // Steelers"), unlike its per-event codes.
  function venueTeamOf(bet, side) {
    const venueTeam = side === "away" ? bet.awayTeamVenue : bet.homeTeamVenue;
    const name = side === "away" ? bet.awayTeam : bet.homeTeam;
    const plainName = typeof name === "string" && name !== "" ? name : null;
    if (venueTeam && typeof venueTeam === "object") {
      const id = typeof venueTeam.id === "string" && venueTeam.id !== "" ? venueTeam.id : null;
      const venueName = (typeof venueTeam.name === "string" && venueTeam.name) || (typeof venueTeam.shortName === "string" && venueTeam.shortName) || plainName;
      if (id) return { key: id, name: venueName };
      return venueName ? { key: venueName, name: venueName } : null;
    }
    return plainName ? { key: plainName, name: plainName } : null;
  }

  function crosswalkKeyOf(venue, league, venueTeamKey) {
    return `${venue}|${league}|${venueTeamKey}`;
  }

  // crosswalk key -> row over the served rows; a malformed row is skipped.
  function crosswalkIndex(rows) {
    const index = new Map();
    for (const row of Array.isArray(rows) ? rows : []) {
      if (!row || typeof row !== "object") continue;
      if (typeof row.venue !== "string" || typeof row.league !== "string" || typeof row.venueTeamKey !== "string") continue;
      if (row.unabatedTeamId == null || row.unabatedTeamId === "") continue;
      index.set(crosswalkKeyOf(row.venue, row.league, row.venueTeamKey), row);
    }
    return index;
  }

  // The team key from the Unabated team id the venue itself sent for one
  // side of a bet, or null. Nothing is guessed: the id is the venue's, and a
  // non-numeric or missing one leaves the side to the crosswalk and names.
  function venueUnabatedKeyOf(bet, side) {
    const venueTeam = side === "away" ? bet.awayTeamVenue : bet.homeTeamVenue;
    const id = venueTeam && typeof venueTeam === "object" ? venueTeam.unabatedId : null;
    return typeof id === "string" && /^\d+$/.test(id) ? teams.keyOf(bet.league, id) : null;
  }

  // The team key a crosswalk row gives one side of a bet, or null.
  function learnedKeyOf(bet, side, index) {
    if (index.size === 0 || !bet.venue) return null;
    const venueTeam = venueTeamOf(bet, side);
    const row = venueTeam ? index.get(crosswalkKeyOf(bet.venue, bet.league, venueTeam.key)) : null;
    return row ? teams.keyOf(bet.league, row.unabatedTeamId) : null;
  }

  // Fill awayKey / homeKey: a crosswalk row for the venue team first (it was
  // learned from an id join and is applied whatever the name resolves to),
  // then the venue's own copy of Unabated's team id where it sends one
  // (Novig's team objects carry `unabatedId`; awayTeamVenue.unabatedId), else
  // the raw team name through teams.js where the key is still null. The bets
  // service leaves both keys null (the team table lives here, not in
  // Python); records that already carry keys and have no crosswalk row, or
  // have no league / no name, are returned unchanged.
  function resolveTeamKeys(records, crosswalk) {
    const index = crosswalkIndex(crosswalk);
    return records.map((record) => {
      if (!record.league) return record;
      const resolved = Object.assign({}, record);
      const awayLearned = learnedKeyOf(record, "away", index) || venueUnabatedKeyOf(record, "away");
      const homeLearned = learnedKeyOf(record, "home", index) || venueUnabatedKeyOf(record, "home");
      if (awayLearned) resolved.awayKey = awayLearned;
      else if (resolved.awayKey == null && resolved.awayTeam != null) resolved.awayKey = teams.teamKey(record.league, record.awayTeam);
      if (homeLearned) resolved.homeKey = homeLearned;
      else if (resolved.homeKey == null && resolved.homeTeam != null) resolved.homeKey = teams.teamKey(record.league, record.homeTeam);
      return resolved;
    });
  }

  // Keys from scratch: what a record resolves to under THIS crosswalk, so a
  // cleared table takes its keys back and a grown one applies to records that
  // were keyed by name before it was learned.
  function rekeyRecords(records, crosswalk) {
    return resolveTeamKeys(records.map((record) => (record.league ? Object.assign({}, record, { awayKey: null, homeKey: null }) : record)), crosswalk);
  }

  // Both sides of one id-joined bet as crosswalk rows to write, or the
  // conflict that stops the whole bet (a name pointing at another team on
  // either side makes the venue's orientation suspect on both). The bet's
  // own contract on the joined row is the orientation check that needs no
  // name: it says which Unabated side the bet sits on, and the bet's own
  // `side` says which venue side — they must be the same side, or the
  // venue's away/home is not Unabated's for this game.
  function lessonOf(bet, row, known) {
    const rows = [];
    const contractSideIndex = sideIndexByVenueId(bet, row);
    const betSideIndex = bet.side === "away" ? SIDE_AWAY_OR_OVER : bet.side === "home" ? SIDE_HOME_OR_UNDER : null;
    if (contractSideIndex != null && betSideIndex != null && contractSideIndex !== betSideIndex) {
      return { rows: [], conflict: { betId: bet.id, venue: bet.venue, league: bet.league, side: bet.side, venueTeamKey: null, boardKey: null,
        reason: `the bet's own contract puts it on Unabated side ${contractSideIndex}, its venue side is ${bet.side}` } };
    }
    for (const side of ["away", "home"]) {
      const venueTeam = venueTeamOf(bet, side);
      if (!venueTeam) return { rows: [], conflict: null };
      const teamId = String(side === "away" ? row.awayTeamId : row.homeTeamId);
      const boardKey = teams.keyOf(bet.league, teamId);
      const name = side === "away" ? bet.awayTeam : bet.homeTeam;
      const nameKey = name == null ? null : teams.teamKey(bet.league, name);
      const conflict = { betId: bet.id, venue: bet.venue, league: bet.league, side, venueTeamKey: venueTeam.key, boardKey };
      if (nameKey && nameKey !== boardKey) {
        return { rows: [], conflict: Object.assign(conflict, { reason: `name "${name}" resolves to ${nameKey}, the joined event says ${boardKey}` }) };
      }
      const held = known.get(crosswalkKeyOf(bet.venue, bet.league, venueTeam.key));
      if (held && String(held.unabatedTeamId) !== teamId) {
        return { rows: [], conflict: Object.assign(conflict, { reason: `crosswalk holds ${teams.keyOf(bet.league, held.unabatedTeamId)}, the joined event says ${boardKey}` }) };
      }
      if (held) continue;
      const unabatedTeamName = side === "away" ? row.awayTeam : row.homeTeam;
      rows.push({
        venue: bet.venue, league: bet.league, venueTeamKey: venueTeam.key, venueTeamName: venueTeam.name,
        unabatedTeamId: teamId, unabatedTeamName: typeof unabatedTeamName === "string" ? unabatedTeamName : null,
        learnedFrom: `${bet.id} on board event ${row.eventId}`,
      });
    }
    return { rows, conflict: null };
  }

  // What the open bets teach against this board: {learned, conflicts}.
  // `learned` holds rows not yet in `crosswalk` (one per venue team, two
  // bets on one game teach the same rows once) for the service to write;
  // `conflicts` the bets that were refused and why. Nothing is written here.
  function learnCrosswalk(betRecords, lines, crosswalk) {
    const known = crosswalkIndex(crosswalk);
    const board = boardOf(lines);
    const rowByIdentity = new Map();
    for (const line of lines) {
      if (line.eventId == null) continue;
      const identity = eventIdentity(line);
      if (!rowByIdentity.has(identity)) rowByIdentity.set(identity, line);
    }
    const learned = new Map();
    const conflicts = [];
    for (const bet of betRecords) {
      if (!isMatchable(bet) || !bet.venue || !bet.league) continue;
      const game = gameOf(bet, board);
      if (game.join === JOIN_NAME || game.events.size !== 1) continue;
      const row = rowByIdentity.get(game.events.values().next().value);
      if (!row || row.awayTeamId == null || row.homeTeamId == null) continue;
      const lesson = lessonOf(bet, row, known);
      if (lesson.conflict) {
        conflicts.push(lesson.conflict);
        continue;
      }
      for (const entry of lesson.rows) learned.set(crosswalkKeyOf(entry.venue, entry.league, entry.venueTeamKey), entry);
    }
    return { learned: Array.from(learned.values()), conflicts };
  }

  // Newest record per id across consecutive service payloads (arrays of
  // records or {bets: [...]}); a later payload wins a tie on sourceFetchedAt.
  function dedupeByNativeId(payloads) {
    const newest = new Map();
    for (const payload of payloads) {
      const records = Array.isArray(payload) ? payload : payload && Array.isArray(payload.bets) ? payload.bets : [];
      for (const record of records) {
        const held = newest.get(record.id);
        if (!held || Date.parse(record.sourceFetchedAt || 0) >= Date.parse(held.sourceFetchedAt || 0)) newest.set(record.id, record);
      }
    }
    return Array.from(newest.values());
  }

  const api = {
    TIE_CAVEAT, GAME_SERIES, RETENTION_DAYS_DEFAULT,
    normalizeKalshi, parseEventSuffix, centsToAmerican,
    AXIS_TOTAL, AXIS_MARGIN,
    matchBets, annotateRows, linePosition, unmatchedReasons, pruneForRetention, dedupeByNativeId, resolveTeamKeys,
    rekeyRecords, learnCrosswalk, venueTeamOf,
    describeBet, formatPlacedAt, formatStake, tierLabel, venueLabel,
  };

  if (typeof module !== "undefined" && module.exports) {
    module.exports = api;
  } else {
    root.UnabatedBets = api;
  }
})(typeof globalThis !== "undefined" ? globalThis : this);
