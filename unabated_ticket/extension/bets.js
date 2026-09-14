// Bet history for the Unabated Ticket panel: the normalised bet record, the
// Kalshi normaliser, and the matcher that flags a line as already bet, bet on
// the other side, or on a game you already have a position in (#114). Pure:
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
//          or Over, 1 = home or Under), points, rotation.
// Outputs matches per line ({tier, bet, label}), per-row annotations for the
//          Edges list, the unmatched list with a reason per bet, the retention
//          prune, and the native-id dedupe. Nothing here writes anywhere.
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

  const TIER_RANK = { same_line: 0, same_side: 1, opposite: 2, same_game: 3 };
  const START_TOLERANCE_MS = 30 * 60 * 1000;
  const DATE_TOLERANCE_DAYS = 1;
  const DAY_MS = 24 * 3600 * 1000;
  const RETENTION_DAYS_DEFAULT = 30;
  const EASTERN = "America/New_York";
  const TIE_CAVEAT = "kalshi_no_side_includes_tie";
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

  function venueLabel(venue) {
    return venue ? venue.charAt(0).toUpperCase() + venue.slice(1) : "unknown venue";
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
    if (!bet.eventDate) return false;
    const delta = Math.abs(dateStringToUtcMs(easternDateOf(line.eventStartMs)) - dateStringToUtcMs(bet.eventDate));
    return delta <= DATE_TOLERANCE_DAYS * DAY_MS;
  }

  function gameMatches(bet, line) {
    if (bet.unmatchable || bet.league !== line.league) return false;
    if (!sameTeamPair(bet, lineTeamKeys(line)) && !rotationMatches(bet, line)) return false;
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

  function betSideKey(bet) {
    if (bet.betType === "total") return bet.side;
    return bet.side === "away" ? bet.awayKey : bet.homeKey;
  }

  function betOtherSideKey(bet) {
    if (bet.betType === "total") return bet.side === "over" ? "under" : "over";
    return bet.side === "away" ? bet.homeKey : bet.awayKey;
  }

  function samePoints(a, b) {
    return (a == null && b == null) || (typeof a === "number" && a === b);
  }

  function tierOf(bet, line) {
    if (bet.betType !== lineBetType(line) || bet.period !== line.period) return "same_game";
    const lineKey = lineSideKey(line);
    if (lineKey == null) return "same_game";
    if (betSideKey(bet) === lineKey) return samePoints(bet.points, line.points) ? "same_line" : "same_side";
    if (betOtherSideKey(bet) === lineKey) return "opposite";
    return "same_game";
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
    same_game: "game",
  };

  function tierLabel(tier) {
    return TIER_LABELS[tier] || tier;
  }

  // The bet itself: what, how much, where — and the line's own number when it
  // differs from the bet's, which is the whole point of those two tiers. No
  // placed-at: it never told you which bet was which (user decision 2026-09-14).
  function labelOf(tier, bet, line) {
    const parts = [describeBet(bet), formatStake(bet.stake), venueLabel(bet.venue)];
    const numberMoved = tier === "same_side" || (tier === "opposite" && !oppositeSameNumber(bet, line));
    if (numberMoved) parts.push(`now ${linePointsLabel(line)}`);
    return parts.join(" · ");
  }

  // ---- public API ------------------------------------------------------------

  // Matches for one line, strongest tier first. options.lines is the whole
  // board (defaults to [line]) — a bet whose game rule accepts two distinct
  // board events is ambiguous and lands in `unmatched`, never in `matches`.
  // Closed and settled bets never match and are not listed here.
  function matchBets(line, bets, options) {
    const board = options && Array.isArray(options.lines) ? options.lines : [line];
    const matches = [];
    const unmatched = [];
    for (const bet of bets) {
      if (!isMatchable(bet) || !gameMatches(bet, line)) continue;
      if (candidateEvents(bet, board).size > 1) {
        unmatched.push({ bet, reason: "ambiguous game" });
        continue;
      }
      const tier = tierOf(bet, line);
      matches.push({ tier, bet, label: labelOf(tier, bet, line) });
    }
    matches.sort((a, b) => TIER_RANK[a.tier] - TIER_RANK[b.tier] || (b.bet.stake ?? 0) - (a.bet.stake ?? 0));
    return { matches, unmatched };
  }

  // Dollars already risked on this line's market, from its matches: `held` is
  // the same direction (same_line + same_side — a different number is still
  // the same opinion), `against` the other side. Stake is dollars risked at
  // every venue, so both compare directly with a Kelly stake. same_game
  // matches carry no dollars here: they do not change how this line is sized.
  function exposureOf(matches) {
    const heldBets = [];
    const againstBets = [];
    for (const match of matches) {
      if (match.tier === "same_line" || match.tier === "same_side") heldBets.push(match.bet);
      else if (match.tier === "opposite") againstBets.push(match.bet);
    }
    const sum = (list) => roundCents(list.reduce((total, bet) => total + (typeof bet.stake === "number" ? bet.stake : 0), 0));
    return { held: sum(heldBets), against: sum(againstBets), heldBets, againstBets };
  }

  // Per Edges row: the strongest tier (or null), its matches, and the dollars
  // already on the market. A bet you hold never hides a line — it changes the
  // size of the next one (see betsview.stakeAdvice).
  function annotateRows(rows, bets) {
    return rows.map((row) => {
      const { matches } = matchBets(row, bets, { lines: rows });
      const tier = matches.length ? matches[0].tier : null;
      return { tier, matches, exposure: exposureOf(matches) };
    });
  }

  function unresolvedTeamNames(bet) {
    const names = [];
    if (bet.awayTeam != null && !bet.awayKey) names.push(bet.awayTeam);
    if (bet.homeTeam != null && !bet.homeKey) names.push(bet.homeTeam);
    return names;
  }

  // Every OPEN bet that matches no line on the board, with why. Closed and
  // settled bets are not problems, so they are not listed.
  function unmatchedReasons(bets, lines) {
    const leaguesOnBoard = new Set(lines.map((line) => line.league));
    const out = [];
    for (const bet of bets) {
      if (bet.status !== "open") continue;
      if (bet.unmatchable) { out.push({ bet, reason: bet.unmatchable }); continue; }
      const unresolved = unresolvedTeamNames(bet);
      if (unresolved.length) { out.push({ bet, reason: `team not recognised (${unresolved.join(", ")})` }); continue; }
      if (!leaguesOnBoard.has(bet.league)) { out.push({ bet, reason: "league not on the scanner" }); continue; }
      const candidates = candidateEvents(bet, lines).size;
      if (candidates > 1) out.push({ bet, reason: "ambiguous game" });
      else if (candidates === 0) out.push({ bet, reason: "no event on the board yet" });
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

  // Fill null awayKey / homeKey from the raw team names. The bets service
  // leaves both null (the team table lives here, not in Python); records that
  // already carry keys, or have no league / no name, are returned unchanged.
  function resolveTeamKeys(records) {
    return records.map((record) => {
      if (!record.league) return record;
      const resolved = Object.assign({}, record);
      if (resolved.awayKey == null && resolved.awayTeam != null) resolved.awayKey = teams.teamKey(record.league, record.awayTeam);
      if (resolved.homeKey == null && resolved.homeTeam != null) resolved.homeKey = teams.teamKey(record.league, record.homeTeam);
      return resolved;
    });
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
    matchBets, annotateRows, exposureOf, unmatchedReasons, pruneForRetention, dedupeByNativeId, resolveTeamKeys,
    describeBet, formatPlacedAt, formatStake, tierLabel,
  };

  if (typeof module !== "undefined" && module.exports) {
    module.exports = api;
  } else {
    root.UnabatedBets = api;
  }
})(typeof globalThis !== "undefined" ? globalThis : this);
