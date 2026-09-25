// Presentation helpers for the bet-history flags in the Unabated Ticket panel
// (#114): source freshness, the header line under the tabs, the Ticket-tab
// banner truncation, the service-payload merge, the ticket -> line shape, and
// the stake advice — the next bet sized against the bets already held on its
// market (#130, conditional Kelly). Pure: no DOM, no fetch, no chrome.* —
// loaded as a plain <script> in panel.html after bets.js, ladder.js and
// condkelly.js (exposes globalThis.UnabatedBetsView) and via require() in
// tests/betsview.test.js. panel.js keeps only the DOM writes.
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

  const inNode = typeof module !== "undefined" && module.exports;
  const bets = inNode ? require("./bets.js") : root.UnabatedBets;
  const kelly = inNode ? require("./kelly.js") : root.UnabatedKelly;
  const ladderLib = inNode ? require("./ladder.js") : root.UnabatedLadder;
  const condkelly = inNode ? require("./condkelly.js") : root.UnabatedCondKelly;

  // Every venue the plan registers a source for (#117 ProphetX still to come;
  // BFA, Wagerzon and Polymarket US added 2026-09-23). Listing them keeps the
  // Bets tab honest about coverage.
  const VENUES = ["kalshi", "betonline", "novig", "prophetx", "bfa", "wagerzon", "polymarket_us"];
  const FRESH_MS = 5 * 60 * 1000;
  const STALE_MS = 60 * 60 * 1000;
  const BANNER_MAX_LINES = 5;
  const DEFAULT_BETS_SETTINGS = { serviceUrl: "http://127.0.0.1:8094" };
  // Tiers on the row's own direction and on the other one; same_game is neither.
  const HELD_TIERS = new Set(["same_line", "same_side", "related_same"]);
  const AGAINST_TIERS = new Set(["opposite", "related_opposite"]);
  // A bet on another market of the game: shown, never sized off (#129).
  const NOTE_OTHER_MARKET = "game \u00b7 not sized";
  // The other direction in another period is a hedge in the real world (1H
  // and FG totals move together, measured link 0.69), but the worst-case
  // pairing would size it as though both bets won together: hold 1H Under
  // $300, new FG Over, K $8,000 -> $407 against $667 alone, where the measured
  // link wants $802. Left out, the stake is the standalone one: no credit for
  // the hedge and no penalty either (user decision 2026-09-19).
  const NOTE_OTHER_PERIOD_HEDGE = "other period \u00b7 not sized";
  // A big edge on a heavy favourite ((1 + edge) / decimal >= 1) leaves the
  // new bet no losing outcome to weigh the held bets against.
  const REASON_CERTAIN_WIN = "edge implies a certain win";
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

  // "bets: 14 open · 2 not matched to a game · kalshi 20 s · betonline —"
  // needsGame is how many open bets need a game (bets.unmatchedReasons).
  function headerLine(records, payload, now, needsGame) {
    const venues = sourceRows(payload, now).map((row) => `${row.venue} ${row.configured ? row.ageText : "—"}`);
    const flagged = needsGame > 0 ? `${needsGame} not matched to a game` : null;
    return [`bets: ${openCount(records)} open`, flagged, ...venues].filter(Boolean).join(" · ");
  }

  // The Bets tab's red banner, or "" when no open bet needs a game.
  function needsGameBanner(count) {
    if (!(count > 0)) return "";
    if (count === 1) return "1 open bet is not matched to a game. It is left out when the panel sizes your next bet on that game. Attach it below.";
    return `${count} open bets are not matched to a game. They are left out when the panel sizes your next bet on that game. Attach them below.`;
  }

  // At most `max` matches for the banner, strongest first (matchBets already
  // sorts), and how many were cut.
  function bannerLines(matches, max) {
    const limit = typeof max === "number" ? max : BANNER_MAX_LINES;
    return { shown: matches.slice(0, limit), more: Math.max(0, matches.length - limit) };
  }

  // ---- stake advice (#130) ---------------------------------------------------
  //
  // The next bet is sized GIVEN the bets already held on its market of the
  // game (conditional Kelly, condkelly.js), not sized alone and then adjusted
  // by subtracting dollars: dollars at different prices and numbers are not
  // comparable. A matched bet is "in the math" when it sits on the row's
  // axis (totals with totals, spreads and moneylines together), in the row's
  // period or on the row's direction in another one, and Unabated has a fair
  // at its number; every other match is left out and says why. Bets on
  // another market of the game never size it (#129).

  function roundCents(dollars) {
    return Math.round(dollars * 100) / 100;
  }

  // The stake with nothing held: kelly.kellyStakeFromEdge, or null when the
  // line cannot be sized at all (no edge known, not an American price).
  function standaloneStake({ price, edgePct, bankroll, multiplier }) {
    if (edgePct == null) return null;
    try {
      return kelly.kellyStakeFromEdge({ bookPrice: price, edgePct, bankroll, multiplier }).stake;
    } catch (_error) {
      return null;
    }
  }

  // The bet's number as its own side writes it: "36.5", "-8.5".
  function betNumberLabel(bet) {
    if (typeof bet.points !== "number") return "its number";
    return bet.betType === "spread" && bet.points > 0 ? `+${bet.points}` : `${bet.points}`;
  }

  // P(above) at every cut the bet needs that the row does not set itself:
  // {pairs: [[cut, prob], ...]} or {reason}. A cut on the row's own number
  // takes the row's chance (from Unabated's edge), so it needs no rung.
  function ladderPairsFor(position, rowPosition, ladderOf) {
    const rowCuts = position.period === rowPosition.period ? condkelly.cutsNeeded(rowPosition) : [];
    const pairs = [];
    for (const cut of condkelly.cutsNeeded(position)) {
      if (rowCuts.includes(cut)) continue;
      const found = ladderLib.probAbove(ladderOf(position.period, position.axis), cut);
      if (found.reason) return { reason: found.reason };
      pairs.push([cut, found.prob]);
    }
    return { pairs };
  }

  // Why one match is not in the math, or null when it is: {note} | {pairs}.
  function sizingOf(match, rowPosition, ladderOf) {
    if (match.tier === "same_game") return { note: NOTE_OTHER_MARKET };
    if (rowPosition.reason) return { note: rowPosition.reason };
    const position = match.position;
    if (!position || position.reason) return { note: position ? position.reason : "position unknown" };
    if (position.axis !== rowPosition.axis) return { note: NOTE_OTHER_MARKET };
    if (position.period !== rowPosition.period && position.direction !== rowPosition.direction) return { note: NOTE_OTHER_PERIOD_HEDGE };
    const found = ladderPairsFor(position, rowPosition, ladderOf);
    if (found.reason) return { note: `no fair at ${betNumberLabel(match.bet)}` };
    return { pairs: found.pairs };
  }

  function heldBetOf(position) {
    return { group: position.period, cut: position.cut, direction: position.direction, stake: position.stake, toWin: position.toWin };
  }

  // A whole-number row pushes on its number: the push mass comes from the two
  // rungs around it. {pairs} or {reason}.
  function rowPushPairs(rowPosition, ladderOf) {
    const cuts = condkelly.cutsNeeded(rowPosition);
    if (cuts.length === 1) return { pairs: [] };
    const pairs = [];
    for (const cut of cuts) {
      const found = ladderLib.probAbove(ladderOf(rowPosition.period, rowPosition.axis), cut);
      if (found.reason) return { reason: `no fair at ${cut} to price the push` };
      pairs.push([cut, found.prob]);
    }
    return { pairs };
  }

  // What to bet on a line given the open bets matched to it.
  //   line       the row / ticket as the matcher saw it (bets.linePosition reads it)
  //   price      the book's American price; edgePct Unabated's edge in percent
  //   matches    bets.matchBets(...).matches for the line
  //   ladderOf   (period, axis) -> ladder.buildLadder(...) result or null
  // Returns
  //   kind      "none"     nothing held is in the math: `bet` is the standalone stake
  //             "sized"    conditional Kelly ran: `bet` is the stake given the held bets
  //             "declined" held bets belong in the math but the calc was refused
  //                        (`reason`): `bet` is the standalone stake
  //   bet       the dollars to act on (null when the line cannot be sized)
  //   alone     the standalone stake; verb "add" when anything in the math is
  //             on this direction, else "bet"; held / against the dollars in
  //             the math by direction
  //   matches   the input matches, in-math first, each with inMath and note
  //   cappedAt  the line's liquidity when it cut `bet`, else null (capAtLiquidity)
  function stakeAdvice({ line, price, edgePct, bankroll, multiplier, matches, ladderOf, liquidity }) {
    const advice = sizeAgainstHeld({ line, price, edgePct, bankroll, multiplier, matches, ladderOf });
    const capped = capAtLiquidity(advice.bet, liquidity);
    return { ...advice, bet: capped.stake, cappedAt: capped.cappedAt };
  }

  // A stake can never be more than is resting at the price: an exchange line
  // with $17 behind it takes $17, whatever Kelly wants. `liquidity` null or
  // absent (a book that reports none) leaves the stake alone.
  //   {stake, cappedAt}  cappedAt is the liquidity when it cut the stake, else null
  function capAtLiquidity(stake, liquidity) {
    const reported = typeof liquidity === "number" && Number.isFinite(liquidity) && liquidity >= 0;
    if (!reported || typeof stake !== "number" || stake <= liquidity) return { stake, cappedAt: null };
    return { stake: roundCents(liquidity), cappedAt: liquidity };
  }

  function sizeAgainstHeld({ line, price, edgePct, bankroll, multiplier, matches, ladderOf }) {
    const alone = standaloneStake({ price, edgePct, bankroll, multiplier });
    const rowPosition = bets.linePosition(line);
    const readLadder = typeof ladderOf === "function" ? ladderOf : () => null;
    const annotated = (matches || []).map((match) => {
      const sizing = sizingOf(match, rowPosition, readLadder);
      return { ...match, inMath: !sizing.note, note: sizing.note || null, pairs: sizing.pairs || [] };
    });
    const inMath = annotated.filter((match) => match.inMath);
    const sumStakes = (list) => roundCents(list.reduce((total, match) => total + match.position.stake, 0));
    const held = sumStakes(inMath.filter((match) => match.position.direction === rowPosition.direction));
    const against = sumStakes(inMath.filter((match) => match.position.direction !== rowPosition.direction));
    const advice = { kind: "none", bet: alone, alone, verb: held > 0 ? "add" : "bet", held, against, reason: null };
    const finish = (list) => ({ ...advice, matches: orderedForDisplay(list) });
    // No edge of its own: $0 as today. Sizing a hedge is deferred (user, 2026-09-19).
    if (inMath.length === 0 || alone == null || !(alone > 0)) return finish(annotated);

    const decline = (reason) => {
      Object.assign(advice, { kind: "declined", reason, held: 0, against: 0, verb: "bet" });
      return finish(annotated.map((match) => (match.inMath ? { ...match, inMath: false, note: reason } : match)));
    };
    const winProb = (1 + edgePct / 100) / kelly.americanToDecimal(price);
    if (!(winProb < 1)) return decline(REASON_CERTAIN_WIN);
    const push = rowPushPairs(rowPosition, readLadder);
    if (push.reason) return decline(push.reason);
    const ladders = {};
    const addPairs = (period, pairs) => { ladders[period] = (ladders[period] || []).concat(pairs); };
    addPairs(rowPosition.period, push.pairs);
    for (const match of inMath) addPairs(match.position.period, match.pairs);
    const solved = condkelly.solveStake({
      kellyBankroll: bankroll * multiplier,
      candidate: {
        group: rowPosition.period, cut: rowPosition.cut, direction: rowPosition.direction,
        netOdds: kelly.americanToDecimal(price) - 1,
        prob: winProb,
      },
      held: inMath.map((match) => heldBetOf(match.position)),
      ladders,
    });
    if (solved.reason) return decline(solved.reason);
    Object.assign(advice, { kind: "sized", bet: roundCents(solved.stake) });
    return finish(annotated);
  }

  // In-math bets first, each group in the matcher's order (strongest tier,
  // then stake); the working `pairs` are dropped.
  function orderedForDisplay(annotated) {
    const strip = ({ pairs: _pairs, ...match }) => match;
    return annotated.filter((match) => match.inMath).concat(annotated.filter((match) => !match.inMath)).map(strip);
  }

  // The dollars the rail tells the user to bet now. The Min suggested bet
  // filter, alerts, the Ticket and the Copy line all read this one number.
  function suggestedBetAmount(advice) {
    return advice && typeof advice.bet === "number" ? advice.bet : 0;
  }

  // The row's badges: your position in dollars. `held $X` for the bets in the
  // math on this direction, `against $X` for the other — both when both
  // exist. A bet that is on a direction but not in the math (no stake, no
  // fair) still flags it, bare. Anything else on the game is a plain `game`.
  //   flag  {tier, matches, advice} as panel.js builds it
  function badges(flag) {
    if (!flag || !flag.tier) return [];
    const advice = flag.advice || { held: 0, against: 0 };
    const tiers = new Set((flag.matches || []).map((match) => match.tier));
    const onDirection = (tierSet) => Array.from(tiers).some((tier) => tierSet.has(tier));
    const out = [];
    if (advice.held > 0) out.push({ kind: "held", text: `held ${bets.formatStake(advice.held)}` });
    else if (onDirection(HELD_TIERS)) out.push({ kind: "held", text: "held" });
    if (advice.against > 0) out.push({ kind: "against", text: `against ${bets.formatStake(advice.against)}` });
    else if (onDirection(AGAINST_TIERS)) out.push({ kind: "against", text: "against" });
    return out.length ? out : [{ kind: "game", text: "game" }];
  }

  // The advice as words, the same on the row, the Ticket block and the Copy
  // text. `verb` says what the number IS — "add" when topping up a position,
  // "bet" otherwise (user choice 2026-09-13); `alone` is the one small line
  // under it, only when the held bets changed the number; `cap` says the
  // number is all the liquidity there is, whenever liquidity cut it.
  //   {verb, bet, alone, cap}  display strings, or null when nothing is held in
  //                            the math and liquidity did not cut the stake
  function stakeAdviceWords(advice) {
    if (!advice) return null;
    const cap = advice.cappedAt != null ? `all ${bets.formatStake(roundCents(advice.cappedAt))} liq` : null;
    if (advice.kind !== "sized") {
      return cap ? { verb: advice.verb, bet: bets.formatStake(advice.bet), alone: null, cap } : null;
    }
    const changed = roundCents(advice.alone) !== advice.bet;
    return { verb: advice.verb, bet: bets.formatStake(advice.bet), alone: changed ? `${bets.formatStake(roundCents(advice.alone))} alone` : null, cap };
  }

  // One line for the clipboard: "add $188.32, $183 alone", "add $17, all $17 liq, $71.06 alone".
  function stakeAdviceLine(advice) {
    const words = stakeAdviceWords(advice);
    if (!words) return null;
    return [`${words.verb} ${words.bet}`, words.cap, words.alone].filter(Boolean).join(", ");
  }

  // The related bets for one line: the bet itself and a tag. A bet in the
  // math carries how it relates (`this line`, `same side`, `other side`); one
  // that is not carries why (`game · not sized`, `no fair at 36.5`) and is
  // rendered grey.
  function relatedLines(flag) {
    const matches = flag && Array.isArray(flag.matches) ? flag.matches : [];
    return matches.map((match) => ({
      tier: match.tier, inMath: match.inMath === true,
      tag: match.inMath === true ? bets.tierLabel(match.tier) : match.note || bets.tierLabel(match.tier),
      text: match.label,
    }));
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
  // succeeded without them; Cal's pins applied (a pin can give a record its
  // league, so before the keys), team keys filled (through the payload's team
  // crosswalk first, #118 step 4), then the retention prune.
  function mergeServicePayload(storedRecords, payload, now) {
    const fresh = venuesWithFreshPull(payload);
    const listed = new Set((payload.bets || []).map((record) => record.id));
    const kept = (storedRecords || []).filter((record) => !fresh.has(record.venue) || listed.has(record.id));
    const merged = bets.dedupeByNativeId([kept, payload]);
    const pinned = bets.applyPins(merged, pinsOf(payload));
    return bets.pruneForRetention(bets.resolveTeamKeys(pinned, crosswalkOf(payload)), now);
  }

  // The team crosswalk a /bets.json payload carries, or none.
  function crosswalkOf(payload) {
    return payload && Array.isArray(payload.crosswalk) ? payload.crosswalk : [];
  }

  // The manual attaches (pins) a /bets.json payload carries, or none.
  function pinsOf(payload) {
    return payload && Array.isArray(payload.pins) ? payload.pins : [];
  }

  // One Bets-tab line per crosswalk row: what the venue calls the team, what
  // Unabated calls it, and where it was learned. Newest first as served.
  function crosswalkRows(crosswalk) {
    return (Array.isArray(crosswalk) ? crosswalk : []).filter((row) => row && typeof row === "object").map((row) => ({
      what: `${row.venueTeamName || row.venueTeamKey} \u2192 ${row.unabatedTeamName || `team ${row.unabatedTeamId}`}`,
      meta: [bets.venueLabel(row.venue), typeof row.league === "string" ? row.league.toUpperCase() : null,
        row.learnedAt ? `${row.pinnedBetId ? "attached by you" : "learned"} ${bets.formatPlacedAt(row.learnedAt)}` : null].filter(Boolean).join(" \u00b7 "),
      title: row.learnedFrom ? `from ${row.learnedFrom}` : "",
    }));
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
      awayTeamId: ticket.awayTeamId ?? null, homeTeamId: ticket.homeTeamId ?? null,
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
    VENUES, FRESH_MS, STALE_MS, BANNER_MAX_LINES, DEFAULT_BETS_SETTINGS,
    fmtAgeShort, freshnessLevel, sourceRows, serviceStatus, sourcesUnavailable, openCount, headerLine,
    bannerLines, badges, relatedLines, stakeAdvice, capAtLiquidity, suggestedBetAmount, stakeAdviceWords, stakeAdviceLine, venuesWithFreshPull, mergeServicePayload, crosswalkOf, pinsOf, crosswalkRows, needsGameBanner, ticketAsLine, sanitizeBetsSettings,
  };

  if (typeof module !== "undefined" && module.exports) {
    module.exports = api;
  } else {
    root.UnabatedBetsView = api;
  }
})(typeof globalThis !== "undefined" ? globalThis : this);
