// The Teasers tab's words and layout, shared by the panel (panel.js) and the
// phone page (server/phone/phone.js), so both show the same list in the same
// words. Pure: no DOM, no fetch, no chrome.* — loaded as a plain <script>
// after teaser.js and bets.js (exposes globalThis.UnabatedTeaserView) and via
// require() in tests/teaserview.test.js. Nothing here writes anywhere.
//
// Inputs are teaser.js's own outputs (or their JSON, as the server runner
// serves them in /teasers.json): planned tickets (describePlan), open BFA
// teasers (openTeasers), the Legs rows (describeLegs) and the plan summary.
// Outputs are plain objects of text the two pages lay out in their DOM.

(function (root) {
  "use strict";

  const inNode = typeof module !== "undefined" && module.exports;
  const teaserLib = inNode ? require("./teaser.js") : root.UnabatedTeaser;
  const betsLib = inNode ? require("./bets.js") : root.UnabatedBets;
  const feed = inNode ? require("./feed.js") : root.UnabatedFeed;

  // The list shows its first tickets open; the rest sit behind a fold.
  const TICKETS_SHOWN = 3;
  const KELLY_NAMES = { 1: "Full", 0.5: "Half", 0.25: "Quarter" };
  // The Legs list shows these standings; the others are folded.
  const SHOWN_STANDINGS = [teaserLib.STANDING_POOL, teaserLib.STANDING_OUT];
  const FOLD_COUNT_WORDS = [
    [teaserLib.STANDING_BELOW, "below break-even"],
    [teaserLib.STANDING_OTHER_MARKET, "on an open teaser's other market"],
    [teaserLib.STANDING_UNPRICED, "with no fair"],
  ];

  // ---- numbers -----------------------------------------------------------------

  function fmtWholeDollars(value) {
    return value.toLocaleString("en-US", { style: "currency", currency: "USD", maximumFractionDigits: 0 });
  }

  function fmtSignedDollars(value) {
    return `${value >= 0 ? "+" : "-"}${fmtWholeDollars(Math.abs(value))}`;
  }

  function fmtWin(win) {
    return `${(win * 100).toFixed(1)}%`;
  }

  function fmtEv(ev) {
    return `EV ${ev >= 0 ? "+" : ""}${(ev * 100).toFixed(1)}%`;
  }

  function fmtWholePct(fraction) {
    return `${Math.round(fraction * 100)}%`;
  }

  function kellyWords(multiplier) {
    return `${KELLY_NAMES[multiplier] || `${multiplier}x`} Kelly`;
  }

  function plural(count, word) {
    return `${count} ${word}${count === 1 ? "" : "s"}`;
  }

  function ticketSizeText() {
    return `${teaserLib.LEGS_PER_TICKET}-team · pays +${teaserLib.TICKET_NET_ODDS * 100}`;
  }

  function breakEvenText() {
    return `${fmtWin(teaserLib.BREAK_EVEN_WIN)} a leg breaks even at +${teaserLib.TICKET_NET_ODDS * 100}`;
  }

  // ---- status and summary ------------------------------------------------------

  // "NFL · CFB loading", then "NFL · CFB · 15 games on Buckeye's board · 56 of 56 legs priced".
  //   scannerStatus  the scanner's status (leaguesLoaded, leagueErrors), or null before it starts
  //   legs           teaserLegs(...) now, or null before the board is ready
  function statusText(scannerStatus, legs) {
    if (!scannerStatus) return "Starting the scanner…";
    const leagues = teaserLib.TEASER_LEAGUE_IDS.map((id) => {
      const label = feed.LEAGUES[id].label;
      if (scannerStatus.leaguesLoaded.includes(id)) return label;
      return scannerStatus.leagueErrors && scannerStatus.leagueErrors[id] ? `${label} unavailable` : `${label} loading`;
    });
    if (!legs) return leagues.join(" · ");
    const games = new Set(legs.map((leg) => leg.eventId)).size;
    const priced = legs.filter((leg) => leg.win != null).length;
    return `${leagues.join(" · ")} · ${plural(games, "game")} on Buckeye's board · ${priced} of ${plural(legs.length, "leg")} priced`;
  }

  // What emptyText reads off a build: {ticketCount, reason, partialHeldBack, placedCount}.
  //   plan   describePlan(...)
  //   build  planTeasers(...).build
  function planStateOf(plan, build) {
    return { ticketCount: plan.tickets.length, reason: build.reason, partialHeldBack: build.partialHeldBack, placedCount: build.placed.length };
  }

  // Why the list is empty, or null when it has tickets.
  //   plan  planStateOf(...), or null before the board is ready
  function emptyText(plan) {
    if (!plan) return "Waiting for Buckeye's NFL and CFB board…";
    if (plan.ticketCount) return null;
    if (plan.reason === teaserLib.REASON_FEW_LEGS) return `Fewer than ${teaserLib.LEGS_PER_TICKET} games with a priced Buckeye leg right now: nothing to tease.`;
    if (plan.reason) return `No tickets: ${plan.reason}.`;
    if (plan.partialHeldBack > 0) {
      return `Nothing more to bet: the next ticket would be a partial ${fmtWholeDollars(plan.partialHeldBack)}, and a 4-team ticket under $200 is already open at BFA (one partial ticket per set).`;
    }
    return plan.placedCount
      ? "Nothing more to bet: no other ticket raises the Kelly growth with the open teasers held."
      : "No ticket worth betting right now: no 4-team ticket raises the Kelly growth at these fairs.";
  }

  // Tickets, dollars, expected profit, the chance to make money and the
  // chance every ticket loses — over the open teasers too when there are any.
  // Null when there is nothing to show.
  //   summary      describePlan(...).summary
  //   openCount    every open BFA teaser (the in-play ones are summary.placedCount)
  //   stake        {bankroll, multiplier}
  // Returns {label, stake, cells: [{label, value, small}], note}.
  function summaryView(summary, openCount, stake) {
    if (!summary || (summary.count === 0 && summary.placedCount === 0)) return null;
    const label = summary.count === 0 ? "Nothing more to bet"
      : summary.placedCount ? `Bet ${summary.count} more` : `Bet ${plural(summary.count, "ticket")}`;
    const cells = [];
    if (summary.placedCount) {
      // The open tickets still in the math; the Open at BFA header counts every open one.
      const allOpen = openCount === summary.placedCount;
      cells.push({ label: allOpen ? "Placed" : "Placed, in play", value: fmtWholeDollars(summary.placedStake), small: null });
      if (summary.expectedAll != null) cells.push({ label: `Expected, all ${summary.count + summary.placedCount}`, value: fmtSignedDollars(summary.expectedAll), small: null });
    } else {
      cells.push({ label: "Expected", value: fmtSignedDollars(summary.expected), small: summary.stake > 0 ? fmtWholePct(summary.expected / summary.stake) : null });
    }
    if (summary.makesMoney != null) cells.push({ label: "Makes money", value: fmtWholePct(summary.makesMoney), small: null });
    if (summary.allLose != null) cells.push({ label: "All lose", value: fmtWholePct(summary.allLose), small: null });
    const straights = summary.straights;
    const note = `${kellyWords(stake.multiplier)} on ${fmtWholeDollars(stake.bankroll)} · `
      + (summary.placedCount ? "open BFA teasers held fixed" : "tickets that share a leg are sized together")
      + (straights.count ? ` · counts ${fmtWholeDollars(straights.stake)} of straight bets on ${plural(straights.games, "game")}` : "");
    return { label, stake: fmtWholeDollars(summary.stake), cells, note };
  }

  // ---- tickets -----------------------------------------------------------------

  // A ticket to bet: {number, size, stake, ev, legs: [{label, from, win}]}.
  function pendingTicketView(ticket) {
    return {
      number: ticket.number, size: ticketSizeText(), stake: fmtWholeDollars(ticket.stake), ev: fmtEv(ticket.ev),
      legs: ticket.legs.map((leg) => ({ label: leg.label, from: leg.fromLabel, win: fmtWin(leg.win) })),
    };
  }

  // "4 more tickets · $800": the fold under the first TICKETS_SHOWN.
  function moreTicketsText(rest) {
    return `${plural(rest.length, "more ticket")} · ${fmtWholeDollars(rest.reduce((sum, ticket) => sum + ticket.stake, 0))}`;
  }

  // Why an open ticket is not in the math, by its legs: every game started,
  // or no leg joined a game Buckeye's board prices (a basketball teaser, CFB
  // not loaded, a start the board disagrees with by over 12 h).
  function outOfMathText(ticket) {
    if (ticket.reason) return `${ticket.reason}: out of the math`;
    if (ticket.legs.every((leg) => leg.state === teaserLib.LEG_STARTED)) return "every game has started: out of the math";
    return "no leg priced on the board: out of the math";
  }

  // One open BFA teaser, read-only: its legs at BFA's numbers, each at its
  // current fair or why it counts as won.
  //   {inPlay, size, ev (null when not in play), legs: [{label, note, win, counted}], outOfMath, placedTag}
  function openTicketView(ticket) {
    const hasDollars = ticket.reason == null;
    return {
      inPlay: Boolean(ticket.inPlay),
      size: hasDollars ? `${ticket.legCount}-team · ${fmtWholeDollars(ticket.stake)} to win ${fmtWholeDollars(ticket.toWin)}` : `${ticket.legCount}-team`,
      ev: ticket.inPlay ? fmtEv(ticket.winAll * (1 + ticket.toWin / ticket.stake) - 1) : null,
      legs: ticket.legs.map((leg) => ({
        label: leg.label, note: leg.note || "", win: leg.win == null ? "—" : fmtWin(leg.win), counted: leg.state !== teaserLib.LEG_LIVE,
      })),
      outOfMath: ticket.inPlay ? null : outOfMathText(ticket),
      placedTag: ticket.placedAt ? `placed ${betsLib.formatPlacedAt(ticket.placedAt)}` : "placed",
    };
  }

  // The Open at BFA header: {count, note}. Always shown once the board is up,
  // so "none open · BFA pulled 41 s ago" says how fresh the answer is.
  //   bfaRow  betsview.sourceRows' BFA row, or null before the bets service is reached
  function openBlockView(placed, bfaRow) {
    const pull = !bfaRow ? "bets service not reached yet"
      : !bfaRow.configured ? "no BFA account read"
        : bfaRow.fetchedAt ? `BFA pulled ${bfaRow.ageText} ago` : "no BFA pull yet";
    const dollars = (tickets) => tickets.reduce((sum, ticket) => sum + (typeof ticket.stake === "number" ? ticket.stake : 0), 0);
    const inPlay = placed.filter((ticket) => ticket.inPlay);
    const amount = inPlay.length === placed.length ? fmtWholeDollars(dollars(placed))
      : `${fmtWholeDollars(dollars(placed))}, ${fmtWholeDollars(dollars(inPlay))} in play`;
    return {
      count: placed.length ? String(placed.length) : "",
      note: [placed.length ? `${amount} · each until its last game starts` : "none open", pull].join(" · "),
    };
  }

  // ---- legs --------------------------------------------------------------------

  function legStandingText(row) {
    const open = row.openTickets ? ` · ${row.openTickets} open` : "";
    if (row.standing !== teaserLib.STANDING_POOL) return `${row.note}${open}`;
    return `${row.inTickets ? `in ${plural(row.inTickets, "ticket")}` : "in no ticket"}${open}`;
  }

  // The straight bets on a pool leg's game, the Edges tab's tags: `held $X`
  // on the leg's side, `against $Y` on the other, or a bare `game` when the
  // game has straights and none of them counts. Every bet in the title.
  //   [{kind, text, title}]
  function straightTags(straights) {
    if (!straights) return [];
    const title = straights.counted.map((bet) => bet.label)
      .concat(straights.leftOut.map((bet) => `${bet.label} — not counted: ${bet.reason}`)).join("\n");
    const tags = [];
    if (straights.held > 0) tags.push({ kind: "held", text: `held ${betsLib.formatStake(straights.held)}`, title });
    if (straights.against > 0) tags.push({ kind: "against", text: `against ${betsLib.formatStake(straights.against)}`, title });
    if (tags.length === 0 && straights.leftOut.length > 0) tags.push({ kind: "game", text: "game", title });
    return tags;
  }

  // A college leg on show (in the pool or above break-even) offers Can't
  // tease; a marked one offers Restore. NFL legs never: Buckeye teases every
  // NFL game. {action: "block" | "restore", title} or null.
  function blockControl(row) {
    if (row.standing === teaserLib.STANDING_BLOCKED) return { action: "restore", title: "Tease this market again." };
    const onShow = row.standing === teaserLib.STANDING_POOL || row.standing === teaserLib.STANDING_OUT;
    if (!onShow || !teaserLib.canBlock(row.leg)) return null;
    return {
      action: "block",
      title: `Buckeye won't tease this. Takes the game's ${teaserLib.marketNameOf(row.leg)} (both sides) out of the tickets until the game starts.`,
    };
  }

  // The Legs list's layout: the pool's legs and the priced legs above
  // break-even, best first, with the break-even divider where it falls; the
  // rest of the games behind a fold, the markets marked can't tease first.
  //   {gameCount, items: [{divider: true} | {row}], folded: [row], foldLabel}
  function legsLayout(rows) {
    const gameCount = new Set(rows.map((row) => row.leg.eventId)).size;
    const shown = rows.filter((row) => SHOWN_STANDINGS.includes(row.standing));
    const blockedRows = rows.filter((row) => row.standing === teaserLib.STANDING_BLOCKED);
    const otherFolded = rows.filter((row) => !SHOWN_STANDINGS.includes(row.standing) && row.standing !== teaserLib.STANDING_BLOCKED);
    const folded = [...blockedRows, ...otherFolded];
    const items = [];
    let dividerPlaced = false;
    for (const row of shown) {
      if (!dividerPlaced && row.leg.win < teaserLib.BREAK_EVEN_WIN) {
        items.push({ divider: true });
        dividerPlaced = true;
      }
      items.push({ row });
    }
    if (!dividerPlaced && folded.length) items.push({ divider: true });
    const counts = FOLD_COUNT_WORDS
      .map(([standing, words]) => [otherFolded.filter((row) => row.standing === standing).length, words])
      .filter(([count]) => count > 0).map(([count, words]) => `${count} ${words}`);
    const labelParts = [];
    if (otherFolded.length) labelParts.push(`${plural(otherFolded.length, "more game")}: ${counts.join(", ")}`);
    if (blockedRows.length) labelParts.push(`${blockedRows.length} can't tease`);
    return { gameCount, items, folded, foldLabel: labelParts.join(" · ") };
  }

  // One leg row's words: {pool, blocked, side, matchup, leagueLabel,
  // eventStartMs, modifiedMs, book, teased, win, standing, tags, control}.
  // The start and line age are formatted by each page (its own clock).
  function legRowView(row) {
    const { leg } = row;
    const pool = row.standing === teaserLib.STANDING_POOL;
    return {
      pool, blocked: row.standing === teaserLib.STANDING_BLOCKED,
      side: leg.label, matchup: leg.matchup, leagueLabel: leg.leagueLabel, eventStartMs: leg.eventStartMs, modifiedMs: leg.modifiedMs,
      book: `Buckeye ${leg.bookLabel}`, teased: `teased ${teaserLib.TEASER_POINTS}`,
      win: leg.win == null ? "—" : fmtWin(leg.win), standing: legStandingText(row),
      tags: straightTags(row.straights), control: blockControl(row),
    };
  }

  const api = {
    TICKETS_SHOWN,
    fmtWholeDollars, fmtSignedDollars, fmtWin, fmtEv, fmtWholePct, kellyWords, plural, ticketSizeText, breakEvenText,
    statusText, planStateOf, emptyText, summaryView, pendingTicketView, moreTicketsText, outOfMathText, openTicketView, openBlockView,
    legStandingText, straightTags, blockControl, legsLayout, legRowView,
  };

  if (inNode) {
    module.exports = api;
  } else {
    root.UnabatedTeaserView = api;
  }
})(typeof globalThis !== "undefined" ? globalThis : this);
