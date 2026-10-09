// Unabated Ticket phone page — the DOM side (phone page plan steps 2-3).
// One order route: a Teasers ticket's Place (POST /place_teaser.json, a real
// teaser at BFA, the panel's same two-step confirm). Everything else you bet
// in the book's own app.
//
// Reads (same origin, the bets service that served this page):
//   GET /edges.json     every 15 s — the server runner's Edges list, proxied
//                       by the bets service (502 {error} when the runner is down)
//   GET /teasers.json   every 15 s — the runner's Teasers tab (teasers_payload.js)
//   GET /bets.json      every 30 s — open bets, each venue's last poll, the
//                       crosswalk and the shared Dismiss marks
//   GET /attach.json    while an Attach is open — the runner's games and plan
//   GET /settings.json  on load, on opening Settings and after a save
// The polls run only while the page is visible, and refresh at once when it
// becomes visible again. After a write the lists are read with ?fresh=1, so
// the runner re-reads the bets service first and the change shows at once.
// Writes (Content-Type application/json, all on bets.duckdb through the
// bets service, the same rows the desktop panel writes):
//   PUT /settings.json         the Settings form: changed fields, or null to reset one
//   POST /dismissals.json      Dismiss / Restore an unmatched bet's red flag
//   POST /teaser_blocks.json   Can't tease / Restore a college teaser market
//   POST /pins.json, DELETE /pins.json?betId=   Attach a bet to a game / Undo
//   DELETE /crosswalk.json     Clear the learned team names (two taps)
//   POST /place_teaser.json    Place a Teasers ticket at BFA (money, not just rows)
// Nothing is stored in the browser.
// Words and numbers come from phoneview.js and the extension's own pure
// modules (edgerows.js, teaserview.js, betsview.js, bets.js, kelly.js, feed.js).

(function () {
  "use strict";

  const phoneView = globalThis.UnabatedPhoneView;
  const edgeRows = globalThis.UnabatedEdgeRows;
  const feed = globalThis.UnabatedFeed;
  const betsLib = globalThis.UnabatedBets;
  const betsView = globalThis.UnabatedBetsView;
  const teaserView = globalThis.UnabatedTeaserView;

  const FETCH_TIMEOUT_MS = 10 * 1000;
  // BFA usually answers a placement in seconds but can take about a minute; a
  // timeout here reads "unconfirmed", never "not placed".
  const PLACE_TIMEOUT_MS = 90 * 1000;
  // The panel's Place timings: Bet $X stays armed 6 s, and a confirm within
  // 0.6 s of arming (the second tap of a double tap) is ignored.
  const TEASER_CONFIRM_MS = 6000;
  const TEASER_CONFIRM_MIN_MS = 600;
  // The runner reads /settings.json every 10 s; one more edges read after that shows a save's effect.
  const RUNNER_SETTINGS_LAG_MS = 11 * 1000;
  const RELATED_LINES_ON_A_ROW = 3;
  const MOVE_TAG_CLASS = { fair_to_you: "move-fair", book_away: "move-book", fair_against: "move-against" };
  const PRICE_CHECK_TAG_CLASS = { best: "check-best", better: "", skip: "check-skip", outlier: "check-outlier" };
  // How many worse books the ticket lists before folding the rest into "N more".
  const OTHERS_WORSE_SHOWN = 3;
  // The panel offers full game and first half; a stored other period is shown too.
  const OFFERED_PERIODS = [1, 2];
  const SORT_LABELS = { edge: "edge", stake: "stake", start: "start", exposure: "my exposure" };
  const BOOK_MODE_LABELS = { default: "Default books", all: "All live books", custom: "Pick books" };
  // Typing in Attach's search asks the runner once the typing pauses.
  const ATTACH_SEARCH_DELAY_MS = 250;
  // Clear is two taps: the first arms it for this long (panel.js).
  const CROSSWALK_CLEAR_ARM_MS = 6000;

  const el = (id) => document.getElementById(id);
  const view = {
    fresh: el("fresh"), banners: el("banners"),
    edgesStatus: el("edges-status"), edgesBooks: el("edges-books"), edgesTailFlex: el("edges-tailflex"), edgesPriceCheck: el("edges-pricecheck"),
    edgesList: el("edges-list"), edgesEmpty: el("edges-empty"), edgesCount: el("edges-count"),
    betsCount: el("bets-count"), betsAlert: el("bets-alert"), betsTab: el("tab-bets"),
    betsRisk: el("bets-risk"), betsRiskCaption: el("bets-risk-caption"),
    betsNeedsBlock: el("bets-needs-block"), betsNeeds: el("bets-needs"), betsNeedsCount: el("bets-needs-count"),
    betsFixBlock: el("bets-fix-block"), betsFix: el("bets-fix"), betsFixCount: el("bets-fix-count"),
    betsSources: el("bets-sources"), betsOpen: el("bets-open"), betsOpenCount: el("bets-open-count"),
    betsOpenEmpty: el("bets-open-empty"), betsMatchNote: el("bets-match-note"),
    betsOffboard: el("bets-offboard"), betsOffboardList: el("bets-offboard-list"), betsOffboardCount: el("bets-offboard-count"),
    settingsForm: el("settings-form"), settingsSave: el("settings-save"), settingsStatus: el("settings-status"),
    sheet: el("sheet"), sheetBackdrop: el("sheet-backdrop"), sheetClose: el("sheet-close"), ticketGone: el("ticket-gone"),
    ticketSide: el("ticket-side"), ticketBadges: el("ticket-badges"), ticketBetLine: el("ticket-bet-line"),
    ticketMatchup: el("ticket-matchup"), ticketStart: el("ticket-start"),
    ticketBooksLabel: el("ticket-books-label"), ticketBooks: el("ticket-books"),
    ticketBook: el("ticket-book"), ticketPrice: el("ticket-price"), ticketFair: el("ticket-fair"), ticketEdge: el("ticket-edge"),
    ticketStakeLabel: el("ticket-stake-label"), ticketStake: el("ticket-stake"), ticketContracts: el("ticket-contracts"), ticketUncapped: el("ticket-uncapped"),
    ticketExposure: el("ticket-exposure"), ticketPayoutRow: el("ticket-payout-row"), ticketToWin: el("ticket-to-win"), ticketPayout: el("ticket-payout"),
    ticketRelatedBlock: el("ticket-related-block"), ticketRelated: el("ticket-related"),
    ticketOthersBlock: el("ticket-others-block"), ticketOthers: el("ticket-others"), ticketOthersCount: el("ticket-others-count"),
    ticketMoved: el("ticket-moved"),
    teasersTab: el("tab-teasers"), teasersCount: el("teasers-count"), viewTeasers: el("view-teasers"),
    teasersError: el("teasers-error"), teasersWarning: el("teasers-warning"), teasersStatus: el("teasers-status"), teasersEmpty: el("teasers-empty"),
    teasersSummary: el("teasers-summary"), teasersSummaryLabel: el("teasers-summary-label"), teasersSummaryStake: el("teasers-summary-stake"),
    teasersSummaryCells: el("teasers-summary-cells"), teasersSummaryNote: el("teasers-summary-note"),
    teasersOpen: el("teasers-open"), teasersOpenCount: el("teasers-open-count"), teasersOpenNote: el("teasers-open-note"), teasersOpenList: el("teasers-open-list"),
    teasersListLabel: el("teasers-list-label"), teasersListCount: el("teasers-list-count"), teasersList: el("teasers-list"),
    teasersMore: el("teasers-more"), teasersMoreLabel: el("teasers-more-label"), teasersMoreList: el("teasers-more-list"),
    teasersLegsLabel: el("teasers-legs-label"), teasersLegsCount: el("teasers-legs-count"), teasersLegs: el("teasers-legs"),
    teasersLegsMore: el("teasers-legs-more"), teasersLegsMoreLabel: el("teasers-legs-more-label"), teasersLegsMoreList: el("teasers-legs-more-list"),
    viewBets: el("view-bets"),
    betsCrosswalkCount: el("bets-crosswalk-count"), betsCrosswalkClear: el("bets-crosswalk-clear"), betsCrosswalkError: el("bets-crosswalk-error"),
    betsCrosswalk: el("bets-crosswalk"), betsCrosswalkEmpty: el("bets-crosswalk-empty"),
  };

  const state = {
    activeView: "edges",
    // What each poll left: {payload, okAt, error, errorAt, failingSince}.
    edges: { payload: null, okAt: null, error: null, errorAt: null, failingSince: null },
    teasers: { payload: null, okAt: null, error: null, errorAt: null, failingSince: null },
    bets: { payload: null, okAt: null, error: null, errorAt: null, failingSince: null },
    held: { records: [], crosswalk: [], pins: [], fillFairs: [] },
    // The shared Dismiss marks (bet ids) /bets.json or a Dismiss reply last
    // carried; null from a service that has none.
    dismissedIds: null,
    // Market keys whose Can't tease / Restore is being saved.
    teaserMarksSaving: new Set(),
    // Ticket signature -> its Place button's state across re-renders:
    // {phase: "confirm" | "placing" | "placed" | "refused" | "unconfirmed",
    // message, legsMatch, armedAt}.
    teaserPlace: new Map(),
    // Bet ids whose Dismiss, Restore or Undo is being saved.
    betsSaving: new Set(),
    // The open Attach: {betId, step: "pick" | "confirm", query, picks (the
    // runner's step-1 reply), eventId, swapped, plan (its step-2 reply), busy, error}.
    attach: null,
    betsError: null,
    crosswalkArmedUntil: 0,
    // GET /settings.json: {settings: {field: value | null}, updatedAt}.
    settingsHeld: null,
    settingsDirty: false,
    settingsBusy: false,
    // The open ticket: the row's key and its card's key (null for a flat row), plus the row last shown.
    ticket: null,
  };
  let timers = [];

  // ---- small DOM helpers ----------------------------------------------------------

  function makeEl(tag, className, text) {
    const node = document.createElement(tag);
    if (className) node.className = className;
    if (text != null) node.textContent = text;
    return node;
  }

  function tagEl(kind, text, title) {
    const tag = makeEl("span", `tag ${kind}`, text);
    if (title) tag.title = title;
    return tag;
  }

  // ---- reads ------------------------------------------------------------------------

  // One same-origin JSON request; a non-2xx is an Error carrying the body's `error`.
  async function requestJson(method, path, body, timeoutMs = FETCH_TIMEOUT_MS) {
    const options = { method, cache: "no-store", signal: AbortSignal.timeout(timeoutMs) };
    if (body !== undefined) {
      options.headers = { "Content-Type": "application/json" };
      options.body = JSON.stringify(body);
    }
    let response;
    try {
      response = await fetch(path, options);
    } catch (error) {
      throw new Error(error.name === "TimeoutError" ? `no answer in ${timeoutMs / 1000} s` : error.message, { cause: error });
    }
    let parsed;
    try {
      parsed = await response.json();
    } catch (_error) {
      parsed = null;
    }
    if (!response.ok) {
      const error = new Error(parsed && parsed.error ? parsed.error : `HTTP ${response.status}`);
      error.status = response.status;
      throw error;
    }
    return parsed;
  }

  function pollSucceeded(poll, payload) {
    Object.assign(poll, { payload, okAt: Date.now(), error: null, errorAt: null, failingSince: null });
  }

  function pollFailed(poll, error) {
    const now = Date.now();
    Object.assign(poll, { error: error.message, errorAt: now, failingSince: poll.failingSince ?? now });
  }

  // `fresh` asks the runner to re-read the bets service first (after a write here).
  async function pollEdges(fresh) {
    try {
      pollSucceeded(state.edges, await requestJson("GET", fresh === true ? "/edges.json?fresh=1" : "/edges.json"));
    } catch (error) {
      pollFailed(state.edges, error);
    }
    renderAll();
  }

  async function pollTeasers(fresh) {
    try {
      pollSucceeded(state.teasers, await requestJson("GET", fresh === true ? "/teasers.json?fresh=1" : "/teasers.json"));
    } catch (error) {
      pollFailed(state.teasers, error);
    }
    renderAll();
  }

  // The records as the panel holds them (edgeRows.applyBetsPayload: merged, pins applied).
  async function pollBets() {
    try {
      const body = await requestJson("GET", "/bets.json");
      const applied = edgeRows.applyBetsPayload(state.held, body, Date.now());
      state.held = { records: applied.records, crosswalk: applied.crosswalk, pins: applied.pins, fillFairs: applied.fillFairs };
      if (state.betsSaving.size === 0) state.dismissedIds = applied.dismissals ? applied.dismissals.map((row) => row.betId) : null;
      pollSucceeded(state.bets, { generatedAt: applied.generatedAt, sources: applied.sources });
    } catch (error) {
      pollFailed(state.bets, error);
    }
    renderAll();
  }

  // After a write: the bets now, then both lists on bets the runner read now.
  async function refreshAfterWrite() {
    await pollBets();
    await Promise.all([pollEdges(true), pollTeasers(true)]);
  }

  async function loadSettings() {
    try {
      state.settingsHeld = await requestJson("GET", "/settings.json");
      if (!state.settingsDirty) renderSettings();
    } catch (error) {
      setSettingsStatus(`Could not read the settings: ${error.message}`, "bad");
    }
  }

  function refreshNow() {
    pollEdges();
    pollTeasers();
    pollBets();
  }

  // Polls only while the page is in view; coming back refreshes at once.
  function onVisibility() {
    for (const id of timers) clearInterval(id);
    timers = [];
    if (document.visibilityState !== "visible") return;
    refreshNow();
    timers = [
      setInterval(pollEdges, phoneView.EDGES_POLL_MS), setInterval(pollTeasers, phoneView.EDGES_POLL_MS),
      setInterval(pollBets, phoneView.BETS_POLL_MS),
    ];
  }

  // ---- header ------------------------------------------------------------------------

  function renderTop(now) {
    view.fresh.textContent = phoneView.freshnessText(state.edges, state.bets, now);
    view.banners.replaceChildren(...phoneView.banners(state.edges, state.bets, now)
      .map((banner) => makeEl("div", `warning${banner.level === "bad" ? " bad" : ""}`, banner.text)));
  }

  // ---- Edges -------------------------------------------------------------------------

  function railEl(row, withMove) {
    const rail = makeEl("div", "edge-rail");
    rail.append(makeEl("span", `edge-pct tier-${row.edgeTier}`, phoneView.fmtEdgePct(row.edgePct)));
    if (withMove && row.move) {
      rail.append(tagEl(MOVE_TAG_CLASS[row.move.kind] || "", row.move.label), makeEl("small", "move-detail", row.move.detail));
    }
    const stake = makeEl("span", `edge-stake${row.rail.atSize ? " at-size" : ""}`, row.rail.text);
    if (row.rail.note) stake.append(" ", makeEl("small", null, row.rail.note));
    rail.append(stake);
    return rail;
  }

  function relatedBlock(row) {
    const all = row.related || [];
    if (!all.length) return null;
    const against = all.some((line) => line.inMath && (line.tier === "opposite" || line.tier === "related_opposite"));
    const block = makeEl("div", `related-block${against ? " against" : ""}`);
    block.append(makeEl("div", "related-head", `Related bets · ${all.length}`));
    for (const line of all.slice(0, RELATED_LINES_ON_A_ROW)) block.append(relatedLineEl(line, "related-line"));
    if (all.length > RELATED_LINES_ON_A_ROW) block.append(makeEl("div", "related-more", `+${all.length - RELATED_LINES_ON_A_ROW} more (open the ticket)`));
    return block;
  }

  // One related bet: its tag, its words and "fair then" when its fill fair was saved.
  function relatedLineEl(line, baseClass) {
    const against = line.tier === "opposite" || line.tier === "related_opposite";
    const className = baseClass === "bet-match"
      ? `bet-match tier-${line.tier}${!line.inMath ? " not-sized" : against ? " bad" : ""}`
      : `related-line tier-${line.tier}${line.inMath ? "" : " not-sized"}`;
    const div = makeEl("div", className);
    if (line.title) div.title = line.title;
    const text = makeEl("span", null, line.text);
    if (line.fairThen) text.append(" · ", makeEl("span", "fair-then", line.fairThen));
    div.append(makeEl("span", baseClass === "bet-match" ? "k" : "related-tag", line.tag), text);
    return div;
  }

  // The content column every row and card share: badges and side, market and game, book and price.
  function mainEl(row, now) {
    const main = makeEl("div");
    const side = makeEl("div", "edge-side");
    for (const badge of row.badges || []) side.append(tagEl(badge.kind, badge.text, badge.title));
    if (row.isAlt) side.append(tagEl("", "alt", `Alternate line; this book's main number is ${phoneView.fmtPoints(row.mainPoints)}`));
    side.append(row.sideLabel);
    const meta = makeEl("div", "edge-meta", `${phoneView.marketText(row)} · ${phoneView.describeMatchup(row)} · ${phoneView.fmtStart(row.eventStart)} · `);
    meta.append(makeEl("span", phoneView.untilLevel(row.eventStartMs, now), phoneView.fmtUntil(row.eventStartMs, now)));
    const book = makeEl("div", "edge-book");
    book.append(makeEl("span", "price", phoneView.bookPriceText(row)), makeEl("span", "age", ` · ${phoneView.lineAgeText(row, now)}`));
    main.append(side, meta, book);
    const checkTag = priceCheckTagEl(row);
    if (checkTag) {
      const checkLine = makeEl("div", "edge-pricecheck");
      checkLine.append(checkTag);
      main.append(checkLine);
    }
    return main;
  }

  // The other-books tag the runner computed (pricecheck.priceCheckTag), or null.
  function priceCheckTagEl(row) {
    const tag = row.priceCheckTag;
    return tag ? tagEl(PRICE_CHECK_TAG_CLASS[tag.kind] ?? "", tag.label) : null;
  }

  // The ticket's "Others" row: every other book at the number, best first, yours in its place.
  function renderTicketBooks(row) {
    const check = row.priceCheck;
    const tag = priceCheckTagEl(row);
    view.ticketBooksLabel.hidden = !tag;
    view.ticketBooks.hidden = !tag;
    view.ticketBooks.replaceChildren();
    if (!tag) return;
    const lineEl = (name, price, className) => {
      const div = makeEl("div", className);
      div.append(makeEl("span", null, name), makeEl("span", null, phoneView.fmtAmerican(price)));
      return div;
    };
    const list = makeEl("div", "others-list");
    for (const entry of check.better) list.append(lineEl(entry.bookName, entry.price, "better"));
    list.append(lineEl(`${row.book.name} (you)${check.same ? ` +${check.same} same` : ""}`, row.price, "yours"));
    for (const entry of check.worse.slice(0, OTHERS_WORSE_SHOWN)) list.append(lineEl(entry.bookName, entry.price, "worse"));
    if (check.worse.length > OTHERS_WORSE_SHOWN) list.append(makeEl("div", "worse", `${check.worse.length - OTHERS_WORSE_SHOWN} more, worse`));
    view.ticketBooks.append(tag, list);
  }

  function edgeItem(row, card, now) {
    const li = makeEl("li", `edge-row tier-${row.edgeTier}${row.isBlurred ? " blurred" : ""}`);
    li.dataset.key = row.key;
    if (card) li.dataset.card = card.key;
    li.append(mainEl(row, now), railEl(row, true));
    const related = relatedBlock(row);
    if (related) li.append(related);
    if (card && card.lineCount > 1) {
      li.append(makeEl("div", "card-more", `${card.bookCount} book${card.bookCount === 1 ? "" : "s"} · ${card.lineCount} lines · tap for the others`));
    }
    return li;
  }

  function renderEdges(now) {
    const payload = state.edges.payload;
    if (!payload) {
      view.edgesStatus.textContent = state.edges.error ? "No Edges list yet." : "Loading the Edges list…";
      view.edgesBooks.textContent = "";
      view.edgesList.replaceChildren();
      view.edgesCount.hidden = true;
      return;
    }
    const settings = payload.settings.edges;
    view.edgesStatus.textContent = phoneView.scannerText(payload.scanner, now);
    view.edgesBooks.textContent = [
      phoneView.booksText(payload.books),
      `≥${settings.minEdgePct}%`,
      settings.includeAlts ? "+alts" : null,
      `sorted by ${SORT_LABELS[settings.sortBy] || settings.sortBy}`,
    ].filter(Boolean).join(" · ");
    view.edgesTailFlex.textContent = payload.tailFlex || "";
    view.edgesTailFlex.hidden = !payload.tailFlex;
    view.edgesPriceCheck.textContent = payload.priceCheck || "";
    view.edgesPriceCheck.hidden = !payload.priceCheck;
    const items = payload.items.map((item) => (payload.grouped ? edgeItem(item.best, item, now) : edgeItem(item, null, now)));
    view.edgesList.replaceChildren(...items);
    view.edgesCount.hidden = payload.total === 0;
    view.edgesCount.textContent = String(payload.total);
    view.edgesEmpty.hidden = payload.total > 0 && payload.total <= payload.maxShown;
    view.edgesEmpty.textContent = payload.total === 0
      ? (payload.scanner.phase === "live" ? `No line at or above ${settings.minEdgePct}% edge right now.` : "Waiting for the runner's first snapshot…")
      : `Showing the top ${payload.maxShown} of ${payload.total} ${payload.unit}; raise the minimum edge to see fewer.`;
  }

  // ---- the ticket sheet ---------------------------------------------------------------

  // The row and its card in the current payload, by key; null when it is no longer listed.
  function findTicketRow(ticket) {
    const payload = state.edges.payload;
    if (!payload) return null;
    for (const item of payload.items) {
      const rows = payload.grouped ? [item.best, ...item.others] : [item];
      const row = rows.find((candidate) => candidate.key === ticket.rowKey);
      if (row) return { row, card: payload.grouped ? item : null };
    }
    return null;
  }

  function openTicket(rowKey, cardKey) {
    state.ticket = { rowKey, cardKey, last: null, opened: null };
    view.sheet.hidden = false;
    document.body.classList.add("sheet-open");
    renderTicket(Date.now());
  }

  function closeTicket() {
    state.ticket = null;
    view.sheet.hidden = true;
    document.body.classList.remove("sheet-open");
  }

  function renderTicket(now) {
    if (!state.ticket) return;
    const found = findTicketRow(state.ticket);
    if (found) state.ticket.last = found;
    view.ticketGone.hidden = Boolean(found);
    const shown = found || state.ticket.last;
    if (!shown) {
      closeTicket();
      return;
    }
    const { row, card } = shown;
    if (!state.ticket.opened) state.ticket.opened = { price: row.price, points: row.points };
    const moved = found ? phoneView.lineMovedText(state.ticket.opened, row) : null;
    view.ticketMoved.hidden = !moved;
    view.ticketMoved.textContent = moved || "";
    const ticket = phoneView.ticketView(row, now);
    view.ticketSide.textContent = ticket.sideLabel;
    view.ticketBadges.replaceChildren(...(row.badges || []).map((badge) => tagEl(badge.kind, badge.text, badge.title)));
    view.ticketBetLine.textContent = ticket.betLine;
    view.ticketMatchup.textContent = ticket.matchup;
    view.ticketStart.textContent = ticket.start;
    view.ticketBook.replaceChildren(ticket.book, makeEl("span", "sub", ticket.lineAge));
    view.ticketPrice.textContent = ticket.price;
    view.ticketFair.textContent = ticket.fair;
    view.ticketEdge.replaceChildren(makeEl("span", `edge-pct tier-${ticket.edgeTier}`, ticket.edge));
    if (row.move) {
      view.ticketEdge.append(tagEl(MOVE_TAG_CLASS[row.move.kind] || "", row.move.label), makeEl("span", "sub", row.move.detail));
    }
    renderTicketBooks(row);

    const block = ticket.stakeBlock;
    view.ticketStakeLabel.textContent = block.label;
    view.ticketStake.textContent = block.stake;
    view.ticketStake.classList.toggle("no-edge", block.noEdge);
    view.ticketContracts.hidden = !block.contracts;
    view.ticketContracts.classList.toggle("under", Boolean(block.contracts && block.contracts.under));
    view.ticketContracts.replaceChildren();
    if (block.contracts) {
      view.ticketContracts.append(block.contracts.text);
      if (block.contracts.cost) view.ticketContracts.append(makeEl("span", "cost", ` · ${block.contracts.cost}`));
    }
    view.ticketExposure.hidden = !block.position;
    view.ticketExposure.textContent = block.position;
    view.ticketExposure.classList.toggle("against", block.positionAgainst);
    view.ticketUncapped.hidden = !block.uncapped;
    view.ticketUncapped.textContent = block.uncapped || "";
    view.ticketPayoutRow.hidden = block.toWin == null;
    view.ticketToWin.textContent = block.toWin || "";
    view.ticketPayout.textContent = block.payout || "";

    const related = row.related || [];
    view.ticketRelatedBlock.hidden = related.length === 0;
    view.ticketRelated.replaceChildren(...related.map((line) => relatedLineEl(line, "bet-match")));

    const others = card ? [card.best, ...card.others].filter((other) => other.key !== row.key) : [];
    view.ticketOthersBlock.hidden = others.length === 0;
    view.ticketOthersCount.textContent = others.length ? String(others.length) : "";
    view.ticketOthers.replaceChildren(...others.map((other) => otherLineEl(other, row.sideLabel, card.key, now)));
  }

  // Another line of the same card: rung, book and price, its edge and standalone stake. Tap to open it.
  function otherLineEl(other, sideLabel, cardKey, now) {
    const li = makeEl("li", "other-line");
    li.dataset.key = other.key;
    li.dataset.card = cardKey;
    const main = makeEl("div");
    const rung = other.sideLabel !== sideLabel ? `${other.sideLabel} · ` : "";
    const alt = other.isAlt ? ` · alt of ${phoneView.fmtPoints(other.mainPoints)}` : "";
    main.append(makeEl("span", "ol-price", `${rung}${phoneView.bookPriceText(other)}${alt}`), makeEl("span", "ol-age", phoneView.lineAgeText(other, now)));
    const rail = makeEl("div", "ol-rail");
    rail.append(makeEl("span", "ol-edge", phoneView.fmtEdgePct(other.edgePct)));
    if (other.move) rail.append(tagEl(MOVE_TAG_CLASS[other.move.kind] || "", other.move.label, other.move.detail));
    rail.append(makeEl("span", "ol-stake", other.stake == null ? "—" : phoneView.fmtDollars(other.stake)));
    li.append(main, rail);
    return li;
  }

  // ---- Teasers --------------------------------------------------------------------------

  function summaryCellEl(cell) {
    const div = makeEl("div");
    const value = makeEl("div", "payout-value", cell.value);
    if (cell.small) value.append(" ", makeEl("small", null, cell.small));
    div.append(makeEl("div", "payout-label", cell.label), value);
    return div;
  }

  function teaserTicketEl(ticket) {
    const card = makeEl("li", "ticket-card");
    const head = makeEl("div", "tk-head");
    const rail = makeEl("span", "tk-rail");
    rail.append(makeEl("span", "tk-stake", ticket.stake), makeEl("span", "tk-ev", ticket.ev));
    head.append(makeEl("span", "tk-num", `#${ticket.number}`), makeEl("span", "tk-size", ticket.size), rail);
    const legs = makeEl("ol", "tk-legs");
    for (const leg of ticket.legs) {
      const row = makeEl("li");
      row.append(makeEl("span", "tl-leg", leg.label), makeEl("span", "tl-win", leg.win), makeEl("span", "tl-from", leg.from));
      legs.append(row);
    }
    const actions = makeEl("div", "tk-actions");
    appendPlaceControls(actions, ticket);
    card.append(head, legs, actions);
    return card;
  }

  function placeButton(className, text, action, signature) {
    const button = makeEl("button", className, text);
    button.type = "button";
    Object.assign(button.dataset, { placeAction: action, placeTicket: signature });
    return button;
  }

  // The panel's Place row: Place -> Cancel / Bet $X at BFA -> Placing… ->
  // placed, or why not. An unconfirmed ticket shows no Place: it may be booked.
  function appendPlaceControls(actions, ticket) {
    const placeState = state.teaserPlace.get(ticket.signature) || null;
    const phase = placeState ? placeState.phase : null;
    if (placeState && placeState.message) {
      // A placed ticket whose legs BFA shows differently is red: it is booked, and wrong.
      const tone = phase === "placed" && placeState.legsMatch !== false ? "ok" : phase === "placing" ? "wait" : "bad";
      actions.append(makeEl("span", `tk-place-msg ${tone}`, placeState.message));
    }
    if (phase === "confirm") {
      actions.append(placeButton("btn sm", "Cancel", "cancel", ticket.signature),
        placeButton("btn sm primary", `Bet ${ticket.stake} at BFA`, "confirm", ticket.signature));
    } else if (phase === null || phase === "refused") {
      actions.append(placeButton("btn sm", "Place", "arm", ticket.signature));
    }
  }

  // POST the ticket's request (the runner built it with teaser.placeRequestOf)
  // to the bets service, which bets it at BFA and answers placed / refused /
  // unconfirmed. An HTTP error is a refusal before BFA; no answer at all may
  // still have reached BFA, so it reads unconfirmed.
  async function placeTeaser(signature) {
    const ticket = (state.teasers.payload ? state.teasers.payload.tickets : []).find((shown) => shown.signature === signature);
    const request = ticket ? ticket.place : { error: "the ticket is no longer on the list" };
    if (request.error) {
      state.teaserPlace.set(signature, { phase: "refused", message: `Not placed: ${request.error}` });
      renderAll();
      return;
    }
    state.teaserPlace.set(signature, { phase: "placing", message: "Placing at BFA…" });
    renderAll();
    let next;
    try {
      const reply = await requestJson("POST", "/place_teaser.json", request.body, PLACE_TIMEOUT_MS);
      next = { phase: reply.status, message: reply.message, legsMatch: reply.legsMatch };
    } catch (error) {
      next = error.status
        ? { phase: "refused", message: `Not placed: ${error.message}` }
        : { phase: "unconfirmed", message: `The bets service did not answer (${error.message}). Check BFA's open bets before placing it again.` };
    }
    state.teaserPlace.set(signature, next);
    renderAll();
    if (next.phase === "placed") await Promise.all([pollBets(), pollTeasers(true)]);
  }

  function onTeaserPlaceClick(button, event) {
    const signature = button.dataset.placeTicket;
    const action = button.dataset.placeAction;
    if (action === "arm") {
      const armedAt = Date.now();
      state.teaserPlace.set(signature, { phase: "confirm", message: null, armedAt });
      setTimeout(() => {
        const held = state.teaserPlace.get(signature);
        if (held && held.phase === "confirm" && held.armedAt === armedAt) {
          state.teaserPlace.delete(signature);
          renderAll();
        }
      }, TEASER_CONFIRM_MS);
      renderAll();
      return;
    }
    if (action === "cancel") {
      state.teaserPlace.delete(signature);
      renderAll();
      return;
    }
    const held = state.teaserPlace.get(signature);
    if (action !== "confirm" || !held || held.phase !== "confirm") return;
    if (event.detail > 1 || Date.now() - held.armedAt < TEASER_CONFIRM_MIN_MS) return;
    placeTeaser(signature).catch((error) => console.error("[unabated-ticket] place teaser failed", error));
  }

  function openTeaserEl(ticket) {
    const card = makeEl("li", `ticket-card placed${ticket.inPlay ? "" : " done"}`);
    const head = makeEl("div", "tk-head");
    const rail = makeEl("span", "tk-rail");
    if (ticket.ev) rail.append(makeEl("span", "tk-ev", ticket.ev));
    head.append(makeEl("span", "tk-size", ticket.size), rail);
    const legs = makeEl("ol", "tk-legs");
    for (const leg of ticket.legs) {
      const row = makeEl("li", leg.counted ? "counted" : null);
      row.append(makeEl("span", "tl-leg", leg.label), makeEl("span", "tl-win", leg.win), makeEl("span", "tl-from", leg.note));
      legs.append(row);
    }
    const actions = makeEl("div", "tk-actions");
    if (ticket.outOfMath) actions.append(makeEl("span", null, ticket.outOfMath));
    actions.append(tagEl("held", ticket.placedTag));
    card.append(head, legs, actions);
    return card;
  }

  function teaserLegEl(row, now) {
    const item = makeEl("li", `edge-row leg-row ${row.pool ? "tier-hot" : "tier-thin out"}`);
    const main = makeEl("div");
    const side = makeEl("div", "edge-side", row.side);
    if (row.blocked) side.append(tagEl("dismissed", "can't tease"));
    const meta = makeEl("div", "edge-meta", phoneView.teaserLegMeta(row));
    meta.append(makeEl("span", phoneView.untilLevel(row.eventStartMs, now), phoneView.fmtUntil(row.eventStartMs, now)));
    const book = makeEl("div", "edge-book");
    book.append(makeEl("span", "price", row.book), makeEl("span", "age", ` · ${row.teased} · ${phoneView.fmtLineAge(row.modifiedMs, now)}`));
    main.append(side, meta, book);
    const rail = makeEl("div", "edge-rail");
    rail.append(makeEl("span", `edge-pct tier-${row.pool ? "hot" : "thin"}`, row.win), makeEl("span", "leg-in", row.standing));
    for (const tag of row.tags) rail.append(tagEl(tag.kind, tag.text, tag.title));
    if (row.control) {
      const saving = state.teaserMarksSaving.has(row.marketKey);
      const restore = row.control.action === "restore";
      const button = makeEl("button", restore ? "link-btn" : "chip-btn", saving ? "Saving…" : restore ? "Restore" : "Can't tease");
      button.type = "button";
      button.title = row.control.title;
      button.disabled = saving;
      Object.assign(button.dataset, { teaserMark: restore ? "restore" : "block", marketKey: row.marketKey, eventStartMs: String(row.eventStartMs) });
      rail.append(button);
    }
    item.append(main, rail);
    return item;
  }

  function breakEvenEl(text) {
    const divider = makeEl("li", "be-divider");
    divider.setAttribute("role", "separator");
    divider.append(makeEl("span", null, text));
    return divider;
  }

  function renderTeasers(now) {
    const payload = state.teasers.payload;
    const count = payload ? payload.ticketCount : 0;
    view.teasersCount.hidden = count === 0;
    view.teasersCount.textContent = String(count);
    view.teasersError.hidden = !(payload && payload.error);
    view.teasersError.textContent = payload && payload.error ? payload.error : "";
    const warning = phoneView.teasersWarningText(phoneView.bfaRowOf(state.bets, now));
    view.teasersWarning.hidden = !warning;
    view.teasersWarning.textContent = warning || "";
    if (!payload) {
      view.teasersStatus.textContent = state.teasers.error ? "No Teasers list yet." : "Loading the Teasers list…";
      return;
    }
    view.teasersStatus.textContent = payload.status;
    view.teasersEmpty.hidden = !payload.empty;
    view.teasersEmpty.textContent = payload.empty || "";

    const summary = payload.summary;
    view.teasersSummary.hidden = !summary;
    if (summary) {
      view.teasersSummaryLabel.textContent = summary.label;
      view.teasersSummaryStake.textContent = summary.stake;
      view.teasersSummaryCells.replaceChildren(...summary.cells.map(summaryCellEl));
      view.teasersSummaryNote.textContent = summary.note;
    }

    view.teasersOpen.hidden = !payload.ready;
    const openBlock = teaserView.openBlockView(payload.open, phoneView.bfaRowOf(state.bets, now));
    view.teasersOpenCount.textContent = openBlock.count;
    view.teasersOpenNote.textContent = openBlock.note;
    view.teasersOpenList.replaceChildren(...payload.open.map(openTeaserEl));

    const tickets = payload.tickets;
    view.teasersListLabel.hidden = tickets.length === 0;
    view.teasersListCount.textContent = tickets.length ? String(tickets.length) : "";
    view.teasersList.replaceChildren(...tickets.slice(0, teaserView.TICKETS_SHOWN).map(teaserTicketEl));
    const rest = tickets.slice(teaserView.TICKETS_SHOWN);
    view.teasersMore.hidden = rest.length === 0;
    view.teasersMoreLabel.textContent = payload.moreLabel;
    view.teasersMoreList.replaceChildren(...rest.map(teaserTicketEl));
    // A ticket gone from the list (placed and now open, or rebuilt away) takes its button state with it.
    const shownSignatures = new Set(tickets.map((ticket) => ticket.signature));
    for (const [signature, placeState] of state.teaserPlace) {
      if (!shownSignatures.has(signature) && placeState.phase !== "placing") state.teaserPlace.delete(signature);
    }

    const legs = payload.legs;
    view.teasersLegsLabel.hidden = legs.items.length + legs.folded.length === 0;
    view.teasersLegsCount.textContent = legs.gameCount ? String(legs.gameCount) : "";
    view.teasersLegs.replaceChildren(...legs.items.map((item) => (item.divider ? breakEvenEl(payload.breakEven) : teaserLegEl(item.row, now))));
    view.teasersLegsMore.hidden = legs.folded.length === 0;
    view.teasersLegsMoreLabel.textContent = legs.foldLabel;
    view.teasersLegsMoreList.replaceChildren(...legs.folded.map((row) => teaserLegEl(row, now)));
  }

  // Can't tease (both sides of the game's spread or total, until it starts) or Restore.
  async function setTeaserMark(button) {
    const { marketKey } = button.dataset;
    const blocked = button.dataset.teaserMark === "block";
    state.teaserMarksSaving.add(marketKey);
    renderAll();
    try {
      await requestJson("POST", "/teaser_blocks.json", { marketKey, eventStartMs: Number(button.dataset.eventStartMs), blocked });
      await pollTeasers(true);
    } catch (error) {
      if (state.teasers.payload) state.teasers.payload = { ...state.teasers.payload, error: `${blocked ? "Can't tease" : "Restore"} failed: ${error.message}` };
    } finally {
      state.teaserMarksSaving.delete(marketKey);
      renderAll();
    }
  }

  // ---- Bets -----------------------------------------------------------------------------

  function actionButton(className, text, data, title) {
    const button = makeEl("button", className, text);
    button.type = "button";
    Object.assign(button.dataset, data);
    if (title) button.title = title;
    return button;
  }

  // One bet. options {unmatched, reason, quiet, flagged, attachable, dismissed}:
  // a flagged row offers Dismiss, an attachable one Attach (its panel opens
  // under it), a dismissed one Restore; a matched bet attached by hand, Undo.
  function betEl(entry, options) {
    const { unmatched = false, reason = null, quiet = false, flagged = false, attachable = false, dismissed = false } = options || {};
    const item = entry.item;
    const saving = state.betsSaving.has(item.id);
    const li = makeEl("li", unmatched ? "unmatched" : "matched");
    const main = makeEl("div");
    const what = makeEl("div", "bet-what", item.what);
    const pinned = item.pinned && !unmatched;
    if (pinned) what.append(" ", tagEl("held", "attached"));
    if (dismissed && reason) what.append(" ", tagEl("dismissed", "dismissed"));
    const meta = makeEl("div", "bet-meta", item.meta);
    if (pinned) meta.append(" · ", actionButton("link-btn", saving ? "Saving…" : "Undo", { betAction: "undo", betId: item.id }, "Remove this attach and the names it taught."));
    if (dismissed && reason) meta.append(" · ", actionButton("link-btn", saving ? "Saving…" : "Restore", { betAction: "restore", betId: item.id }, "Flag this bet again."));
    main.append(what, meta);
    const rail = makeEl("div");
    rail.append(makeEl("span", "bet-stake", item.stake));
    if (item.when) rail.append(makeEl("small", "bet-when", item.when));
    li.append(main, rail);
    if (!reason) return li;
    const actions = makeEl("div", "bet-actions");
    actions.append(makeEl("span", `bet-reason${quiet ? " quiet" : ""}`, reason));
    if (flagged) {
      actions.append(actionButton("chip-btn", saving ? "Saving…" : "Dismiss", { betAction: "dismiss", betId: item.id },
        "Stop flagging this bet. It moves to Not on the board and still does not size your next bet."));
    }
    const attachOpen = state.attach && state.attach.betId === item.id;
    if (attachable) actions.append(actionButton(`chip-btn${attachOpen ? "" : " primary"}`, attachOpen ? "Cancel" : "Attach", { betAction: attachOpen ? "attach-cancel" : "attach", betId: item.id }));
    li.append(actions);
    if (attachable && attachOpen) li.append(attachPanelEl(item.id));
    return li;
  }

  // ---- Attach (attach.js, run by the runner on its board) ----

  function attachPanelEl(betId) {
    const panel = makeEl("div", "attach-panel");
    panel.dataset.attachFor = betId;
    if (state.attach.step === "confirm") attachConfirmStep(panel);
    else attachPickStep(panel);
    if (state.attach.error) panel.append(makeEl("div", "attach-error", state.attach.error));
    return panel;
  }

  function attachPickStep(panel) {
    const picks = state.attach.picks;
    const head = makeEl("div", "attach-head", "Which game?");
    head.append(makeEl("span", "scope", picks ? picks.scope : "Loading the board…"));
    const search = makeEl("input", "attach-search");
    Object.assign(search, { type: "search", placeholder: "Search a team", value: state.attach.query, autocomplete: "off" });
    search.dataset.attachSearch = "1";
    const list = makeEl("ol", "candidates");
    for (const event of picks ? picks.events : []) {
      const item = makeEl("li", event.why ? "candidate best" : "candidate");
      item.dataset.attachPick = String(event.eventId);
      item.append(makeEl("div", "c-game", event.label), makeEl("div", "c-meta", event.meta));
      if (event.why) item.append(tagEl("held", event.why));
      list.append(item);
    }
    panel.append(head, search, list);
    if (picks && !picks.events.length) panel.append(makeEl("div", "muted", "No game on the board fits. Search a team by name."));
    if (picks && picks.more) panel.append(makeEl("div", "muted", `${picks.more} more: type to narrow the list.`));
    panel.append(makeEl("div", "muted", "Not listed? The game may not be on the board yet. The bet stays flagged until you attach it."));
  }

  function attachConfirmStep(panel) {
    const plan = state.attach.plan;
    const gameHead = makeEl("div", "attach-head", "Game");
    gameHead.append(actionButton("link-btn", "Change", { betAction: "attach-change", betId: state.attach.betId }));
    panel.append(gameHead);
    if (!plan) {
      panel.append(makeEl("div", "muted", "Loading…"));
      return;
    }
    const chosen = makeEl("div", "c-game", plan.label);
    chosen.append(makeEl("span", "muted", ` · ${plan.meta}`));
    panel.append(chosen);
    if (plan.names.length) {
      const namesHead = makeEl("div", "attach-head", `${plan.venueLabel} calls them`);
      namesHead.append(actionButton("link-btn", "Swap", { betAction: "attach-swap", betId: state.attach.betId }));
      panel.append(namesHead);
      for (const name of plan.names) {
        const row = makeEl("div", "map-row");
        const target = makeEl("span", null, `${name.unabatedTeamName} `);
        target.append(makeEl("span", "muted", name.was ? `${name.eventSide}, was ${name.was}` : name.eventSide));
        row.append(makeEl("span", null, name.venueTeamName), makeEl("span", "muted", "→"), target, tagEl(name.status === "known" ? "" : name.status, name.status));
        panel.append(row);
      }
    } else {
      panel.append(makeEl("div", "muted", "The bet names no team, so nothing is learned. The attach pins the game."));
    }
    const yourBet = makeEl("div");
    yourBet.append(makeEl("span", "muted", "Your bet "), plan.betOnGame);
    const submit = actionButton("btn primary", state.attach.busy ? "Saving…" : plan.learnCount ? "Attach and learn" : "Attach", { betAction: "attach-submit", betId: state.attach.betId });
    submit.disabled = state.attach.busy;
    const learned = plan.learnedNames.map((name) => `"${name}"`).join(" and ");
    panel.append(yourBet, submit, makeEl("div", "muted", plan.learnCount
      ? `Next time ${plan.venueLabel} writes ${learned}, it matches on its own.`
      : "Nothing to learn: every name already matches. This pins the bet to this game."));
  }

  let attachSearchTimer = null;

  async function loadAttachPicks() {
    const { betId, query } = state.attach;
    try {
      const picks = await requestJson("GET", `/attach.json?betId=${encodeURIComponent(betId)}&query=${encodeURIComponent(query)}`);
      if (state.attach && state.attach.betId === betId && state.attach.query === query) state.attach = { ...state.attach, picks, error: null };
    } catch (error) {
      if (state.attach && state.attach.betId === betId) state.attach = { ...state.attach, error: `Could not load the games: ${error.message}` };
    }
    renderBetsKeepingSearch();
  }

  async function loadAttachPlan() {
    const { betId, eventId, swapped } = state.attach;
    try {
      const plan = await requestJson("GET", `/attach.json?betId=${encodeURIComponent(betId)}&eventId=${encodeURIComponent(eventId)}&swapped=${swapped ? 1 : 0}`);
      if (state.attach && state.attach.eventId === eventId && state.attach.swapped === swapped) state.attach = { ...state.attach, plan, error: null };
    } catch (error) {
      if (state.attach && state.attach.betId === betId) state.attach = { ...state.attach, error: `Could not load that game: ${error.message}` };
    }
    renderAll();
  }

  async function submitAttach() {
    const plan = state.attach.plan;
    if (!plan || state.attach.busy) return;
    state.attach = { ...state.attach, busy: true, error: null };
    renderAll();
    try {
      await requestJson("POST", "/pins.json", plan.pinRequest);
      state.attach = null;
      await refreshAfterWrite();
    } catch (error) {
      if (state.attach) state.attach = { ...state.attach, busy: false, error: `Attach failed: ${error.message}` };
      renderAll();
    }
  }

  // The Bets view is redrawn on every poll; keep Attach's search box focused where the caret was.
  function renderBetsKeepingSearch() {
    const active = document.activeElement;
    const caret = active && active.dataset && active.dataset.attachSearch ? active.selectionStart : null;
    renderAll();
    if (caret == null) return;
    const search = view.viewBets.querySelector("[data-attach-search]");
    if (!search) return;
    search.focus();
    search.setSelectionRange(caret, caret);
  }

  // Dismiss / Restore (bets.duckdb::bet_dismissals) and Undo (DELETE /pins.json).
  async function onBetAction(button) {
    const { betAction, betId } = button.dataset;
    if (betAction === "attach") {
      state.attach = { betId, step: "pick", query: "", picks: null, eventId: null, swapped: false, plan: null, busy: false, error: null };
      renderAll();
      loadAttachPicks();
      return;
    }
    if (betAction === "attach-cancel") {
      state.attach = null;
      renderAll();
      return;
    }
    if (betAction === "attach-change") {
      state.attach = { ...state.attach, step: "pick", eventId: null, plan: null, error: null };
      renderAll();
      return;
    }
    if (betAction === "attach-swap") {
      state.attach = { ...state.attach, swapped: !state.attach.swapped, plan: null };
      renderAll();
      loadAttachPlan();
      return;
    }
    if (betAction === "attach-submit") {
      submitAttach();
      return;
    }
    await saveBetMark(betId, betAction);
  }

  async function saveBetMark(betId, betAction) {
    state.betsSaving.add(betId);
    state.betsError = null;
    renderAll();
    try {
      if (betAction === "undo") {
        await requestJson("DELETE", `/pins.json?betId=${encodeURIComponent(betId)}`);
      } else {
        const reply = await requestJson("POST", "/dismissals.json", { betIds: [betId], dismissed: betAction === "dismiss" });
        state.dismissedIds = reply.dismissals.map((row) => row.betId);
      }
    } catch (error) {
      state.betsError = `${{ undo: "Undo", dismiss: "Dismiss", restore: "Restore" }[betAction]} failed: ${error.message}`;
    } finally {
      state.betsSaving.delete(betId);
    }
    await refreshAfterWrite();
  }

  // ---- learned team names ----

  function crosswalkArmed() {
    return Date.now() < state.crosswalkArmedUntil;
  }

  function renderCrosswalk() {
    const rows = betsView.crosswalkRows(state.held.crosswalk);
    view.betsCrosswalkCount.textContent = rows.length ? String(rows.length) : "";
    view.betsCrosswalkClear.disabled = rows.length === 0;
    view.betsCrosswalkClear.classList.toggle("armed", crosswalkArmed());
    view.betsCrosswalkClear.textContent = crosswalkArmed() ? `Clear ${rows.length} row${rows.length === 1 ? "" : "s"}?` : "Clear";
    view.betsCrosswalk.replaceChildren(...rows.map((row) => {
      const li = makeEl("li", "matched");
      li.title = row.title;
      const main = makeEl("div");
      main.append(makeEl("div", "bet-what", row.what), makeEl("div", "bet-meta", row.meta));
      li.append(main);
      return li;
    }));
    view.betsCrosswalkEmpty.hidden = rows.length > 0;
    view.betsCrosswalkEmpty.textContent = "Nothing learned yet. Rows appear when an open bet joins the board by its Kalshi event or Novig outcome id, or when you attach one.";
  }

  async function clearCrosswalk() {
    if (!crosswalkArmed()) {
      state.crosswalkArmedUntil = Date.now() + CROSSWALK_CLEAR_ARM_MS;
      setTimeout(renderAll, CROSSWALK_CLEAR_ARM_MS);
      renderAll();
      return;
    }
    state.crosswalkArmedUntil = 0;
    view.betsCrosswalkClear.disabled = true;
    try {
      await requestJson("DELETE", "/crosswalk.json");
      view.betsCrosswalkError.textContent = "";
    } catch (error) {
      view.betsCrosswalkError.textContent = `Could not clear: ${error.message}`;
    }
    await refreshAfterWrite();
  }

  function venueEl(row) {
    const div = makeEl("div", `venue fresh-${row.level}`);
    const bets = row.count == null ? null : `${row.count} bets`;
    const trouble = row.error || row.note;
    const note = trouble ? [bets, trouble].filter(Boolean).join(" · ") : bets || "no bets";
    const name = makeEl("span");
    name.append(makeEl("span", "vname", betsLib.venueLabel(row.venue)), makeEl("span", "vnote", note));
    div.append(makeEl("span", "vdot"), name, makeEl("span", "vage", row.configured ? row.ageText : "—"));
    return div;
  }

  function renderBets(now) {
    const edgesPayload = state.edges.payload;
    const service = edgesPayload ? edgesPayload.betsService : null;
    const runnerUnmatched = service && Array.isArray(service.unmatched) ? service.unmatched : null;
    const unmatched = phoneView.applyDismissals(runnerUnmatched, state.dismissedIds);
    const tab = phoneView.betsTabView(state.held.records, unmatched, service ? service.boardLineCount : 0, state.bets.payload, now);
    // An Attach whose bet matched, settled or is no longer attachable closes (never mid-save).
    const attachableIds = new Set([...tab.needsGame, ...tab.offBoard].filter((entry) => entry.attachable).map((entry) => entry.betId));
    if (state.attach && !state.attach.busy && !attachableIds.has(state.attach.betId)) state.attach = null;
    view.betsRisk.textContent = tab.atRisk;
    view.betsRiskCaption.textContent = state.betsError ? `${tab.caption} · ${state.betsError}` : tab.caption;
    view.betsNeedsBlock.hidden = tab.needsGame.length === 0;
    view.betsNeedsCount.textContent = String(tab.needsGame.length);
    view.betsNeeds.replaceChildren(...tab.needsGame.map((entry) => betEl(entry, { unmatched: true, reason: entry.reason, flagged: true, attachable: entry.attachable })));
    view.betsFixBlock.hidden = tab.needsFix.length === 0;
    view.betsFixCount.textContent = String(tab.needsFix.length);
    view.betsFix.replaceChildren(...tab.needsFix.map((entry) => betEl(entry, { unmatched: true, reason: entry.reason, flagged: true })));
    view.betsSources.replaceChildren(...tab.venueRows.map(venueEl));
    view.betsOpenCount.textContent = tab.open.length ? String(tab.open.length) : "";
    view.betsOpen.replaceChildren(...tab.open.map((entry) => betEl(entry, { unmatched: entry.unmatched })));
    view.betsOpenEmpty.hidden = tab.open.length > 0;
    view.betsOpenEmpty.textContent = state.bets.okAt != null ? "No open bets." : "No bets loaded yet.";
    view.betsMatchNote.hidden = !tab.matchNote;
    view.betsMatchNote.textContent = tab.matchNote || "";
    view.betsOffboard.hidden = tab.offBoard.length === 0;
    view.betsOffboardCount.textContent = String(tab.offBoard.length);
    view.betsOffboardList.replaceChildren(...tab.offBoard.map((entry) => betEl(entry, {
      unmatched: true, reason: entry.reason, quiet: true, attachable: entry.attachable, dismissed: entry.dismissed,
    })));
    if (state.attach && tab.offBoard.some((entry) => entry.betId === state.attach.betId)) view.betsOffboard.open = true;
    renderCrosswalk();

    const flagged = tab.needsGame.length + tab.needsFix.length;
    view.betsCount.hidden = tab.open.length === 0;
    view.betsCount.textContent = String(tab.open.length);
    view.betsAlert.hidden = flagged === 0;
    view.betsAlert.textContent = String(flagged);
    view.betsTab.classList.toggle("flagged", flagged > 0);
  }

  // ---- Settings ---------------------------------------------------------------------------

  function setSettingsStatus(text, level) {
    view.settingsStatus.textContent = text;
    view.settingsStatus.className = `settings-status${level ? ` ${level}` : ""}`;
  }

  // The hint under a field: whether it is your value or the default, and Reset for yours.
  function hintEl(field, heldValue) {
    const hint = makeEl("div", "field-hint");
    if (heldValue == null) {
      hint.append(phoneView.defaultText(field));
      return hint;
    }
    hint.append(makeEl("span", "held", "your value"), ` · ${phoneView.defaultText(field)}`);
    const reset = makeEl("button", "reset-btn", "Reset to default");
    reset.type = "button";
    reset.dataset.reset = field;
    hint.append(reset);
    return hint;
  }

  function numberField(field, label, value, heldValue, step) {
    const row = makeEl("div", "field");
    const id = `setting-${field}`;
    const labelEl = makeEl("label", "field-label", label);
    labelEl.htmlFor = id;
    const wrap = makeEl("div", "field-input");
    const input = makeEl("input");
    Object.assign(input, { id, type: "number", inputMode: "decimal", step: String(step), min: "0", value: String(value) });
    input.dataset.field = field;
    wrap.append(input);
    row.append(labelEl, wrap, hintEl(field, heldValue));
    return row;
  }

  function toggleField(field, label, value, heldValue) {
    const row = makeEl("div", "field");
    const id = `setting-${field}`;
    const labelEl = makeEl("label", "field-label", label);
    labelEl.htmlFor = id;
    const wrap = makeEl("div", "field-input toggle");
    const input = makeEl("input");
    Object.assign(input, { id, type: "checkbox", checked: Boolean(value) });
    input.dataset.field = field;
    wrap.append(input);
    row.append(labelEl, wrap, hintEl(field, heldValue));
    return row;
  }

  function selectField(field, label, value, heldValue, options) {
    const row = makeEl("div", "field");
    const id = `setting-${field}`;
    const labelEl = makeEl("label", "field-label", label);
    labelEl.htmlFor = id;
    const wrap = makeEl("div", "field-input");
    const select = makeEl("select");
    select.id = id;
    select.dataset.field = field;
    for (const [optionValue, optionLabel] of Object.entries(options)) {
      const option = makeEl("option", null, optionLabel);
      option.value = optionValue;
      option.selected = optionValue === value;
      select.append(option);
    }
    wrap.append(select);
    row.append(labelEl, wrap, hintEl(field, heldValue));
    return row;
  }

  // A row of checkboxes for a list field; `choices` is [[value, label]].
  function checksField(field, label, values, heldValue, choices, extraClass) {
    const row = makeEl("div", "field");
    row.append(makeEl("div", "field-label", label));
    const checks = makeEl("div", `checks${extraClass ? ` ${extraClass}` : ""}`);
    checks.dataset.group = field;
    for (const [choiceValue, choiceLabel] of choices) {
      const labelEl = makeEl("label", "check");
      const input = makeEl("input");
      Object.assign(input, { type: "checkbox", checked: values.includes(choiceValue), value: String(choiceValue) });
      labelEl.append(input, choiceLabel);
      checks.append(labelEl);
    }
    row.append(makeEl("div"), checks, hintEl(field, heldValue));
    return row;
  }

  function section(title) {
    return makeEl("h2", "section-label", title);
  }

  // Every live book the runner sees, plus any stored id it does not (so a tick is never silently dropped).
  function bookChoices(selectedIds) {
    const payload = state.edges.payload;
    const live = payload && payload.books && Array.isArray(payload.books.live) ? payload.books.live : [];
    const choices = live.map((book) => [book.id, book.name]);
    for (const id of selectedIds || []) if (!choices.some(([choiceId]) => choiceId === id)) choices.push([id, `book ${id}`]);
    return choices;
  }

  // The form from the stored settings; `formValues` (what is typed so far)
  // overrides them when the form is redrawn mid-edit (a book mode change).
  function renderSettings(formValues) {
    const held = state.settingsHeld ? state.settingsHeld.settings : null;
    const stored = (field) => (held ? held[field] : null);
    const current = { ...phoneView.effectiveSettings(held), ...(formValues || {}) };
    const periods = Array.from(new Set([...OFFERED_PERIODS, ...current.periods]));
    const form = view.settingsForm;
    form.replaceChildren(
      section("Sizing"),
      numberField("bankroll", "Bankroll $", current.bankroll, stored("bankroll"), 1),
      numberField("multiplier", "Kelly multiplier", current.multiplier, stored("multiplier"), 0.01),
      section("Markets"),
      checksField("leagues", "Sports", phoneView.sportsOfLeagues(current.leagues), stored("leagues"), Object.entries(feed.SPORTS)),
      checksField("periods", "Periods", current.periods, stored("periods"), periods.map((id) => [id, feed.PERIODS[id]])),
      checksField("betTypes", "Bets", current.betTypes, stored("betTypes"), Object.entries(feed.BET_TYPES).map(([id, name]) => [Number(id), name])),
      section("Books"),
      selectField("bookMode", "Books", current.bookMode, stored("bookMode"), BOOK_MODE_LABELS),
    );
    if (current.bookMode === "custom") {
      const books = checksField("bookIds", "Your books", current.bookIds || [], null, bookChoices(current.bookIds), "books-checks");
      books.querySelector(".field-hint").textContent = bookChoices([]).length ? "" : "The runner has no live book list yet.";
      form.append(books);
    }
    form.append(
      section("What is listed"),
      numberField("minEdgePct", "Min edge %", current.minEdgePct, stored("minEdgePct"), 0.1),
      numberField("minStake", "Min suggested bet $", current.minStake, stored("minStake"), 1),
      numberField("maxLineAgeHours", "Max line age h", current.maxLineAgeHours, stored("maxLineAgeHours"), 1),
      numberField("minLiquidityToWin", "Min liq to win $", current.minLiquidityToWin, stored("minLiquidityToWin"), 10),
      toggleField("includeAlts", "Include alt lines", current.includeAlts, stored("includeAlts")),
      selectField("sortBy", "Sort", current.sortBy, stored("sortBy"), SORT_LABELS),
      toggleField("groupByMarket", "Group by market", current.groupByMarket, stored("groupByMarket")),
    );
    if (state.settingsHeld && state.settingsHeld.updatedAt && !state.settingsBusy && !state.settingsDirty) {
      setSettingsStatus(`Saved ${new Date(state.settingsHeld.updatedAt).toLocaleString([], { month: "short", day: "numeric", hour: "numeric", minute: "2-digit" })}`, null);
    }
  }

  const NUMBER_FIELD_LABELS = {
    bankroll: "Bankroll", multiplier: "Kelly multiplier", minEdgePct: "Min edge", minStake: "Min suggested bet",
    maxLineAgeHours: "Max line age", minLiquidityToWin: "Min liq to win",
  };

  // The form as {field: value}: an empty number box is "not given" (left
  // out); a box that is not a number is an error, never sent (JSON would turn
  // NaN into null, a reset).
  function readSettingsForm() {
    const form = view.settingsForm;
    const values = {};
    for (const [field, label] of Object.entries(NUMBER_FIELD_LABELS)) {
      const raw = form.querySelector(`[data-field="${field}"]`).value.trim();
      if (raw === "") continue;
      const number = Number(raw);
      if (!Number.isFinite(number)) return { error: `${label} must be a number, got "${raw}"` };
      values[field] = number;
    }
    const checked = (group) => Array.from(form.querySelectorAll(`[data-group="${group}"] input:checked`)).map((input) => input.value);
    const held = state.settingsHeld ? state.settingsHeld.settings : null;
    values.leagues = phoneView.leaguesForSports(checked("leagues"), phoneView.effectiveSettings(held).leagues);
    values.periods = checked("periods").map(Number);
    values.betTypes = checked("betTypes").map(Number);
    values.bookMode = form.querySelector('[data-field="bookMode"]').value;
    const bookBox = form.querySelector('[data-group="bookIds"]');
    if (bookBox) values.bookIds = checked("bookIds").map(Number);
    values.includeAlts = form.querySelector('[data-field="includeAlts"]').checked;
    values.groupByMarket = form.querySelector('[data-field="groupByMarket"]').checked;
    values.sortBy = form.querySelector('[data-field="sortBy"]').value;
    return { values };
  }

  async function putSettings(update, doneText) {
    state.settingsBusy = true;
    view.settingsSave.disabled = true;
    setSettingsStatus("Saving…", null);
    try {
      const reply = await requestJson("PUT", "/settings.json", { settings: update });
      state.settingsHeld = { settings: reply.settings, updatedAt: reply.updatedAt };
      state.settingsDirty = false;
      state.settingsBusy = false;
      renderSettings();
      setSettingsStatus(`${doneText} The runner picks it up within 10 s.`, "ok");
      pollEdges();
      setTimeout(pollEdges, RUNNER_SETTINGS_LAG_MS);
    } catch (error) {
      // A 400 names the field and value the service refused.
      setSettingsStatus(`Not saved: ${error.message}`, "bad");
    } finally {
      state.settingsBusy = false;
      view.settingsSave.disabled = false;
    }
  }

  function onSave() {
    const parsed = readSettingsForm();
    if (parsed.error) {
      setSettingsStatus(parsed.error, "bad");
      return;
    }
    const held = state.settingsHeld ? state.settingsHeld.settings : null;
    const update = phoneView.settingsUpdate(held, parsed.values);
    if (Object.keys(update).length === 0) {
      setSettingsStatus("Nothing changed.", null);
      state.settingsDirty = false;
      return;
    }
    putSettings(update, "Saved.");
  }

  // ---- views and events ----------------------------------------------------------------------

  function renderAll() {
    const now = Date.now();
    renderTop(now);
    renderEdges(now);
    renderTeasers(now);
    renderBets(now);
    renderTicket(now);
  }

  function showView(name) {
    state.activeView = name;
    for (const button of document.querySelectorAll(".tabbar button")) {
      const active = button.dataset.view === name;
      button.classList.toggle("active", active);
      button.setAttribute("aria-selected", String(active));
    }
    for (const viewName of ["edges", "teasers", "bets", "settings"]) el(`view-${viewName}`).hidden = viewName !== name;
    if (name === "settings" && !state.settingsDirty) loadSettings();
    window.scrollTo(0, 0);
  }

  document.querySelector(".tabbar").addEventListener("click", (event) => {
    const button = event.target.closest("button[data-view]");
    if (button) showView(button.dataset.view);
  });

  view.edgesList.addEventListener("click", (event) => {
    const row = event.target.closest(".edge-row");
    if (row) openTicket(row.dataset.key, row.dataset.card || null);
  });

  view.ticketOthers.addEventListener("click", (event) => {
    const line = event.target.closest(".other-line");
    if (!line) return;
    openTicket(line.dataset.key, line.dataset.card || null);
    view.sheet.querySelector(".sheet-panel").scrollTop = 0;
  });

  view.viewTeasers.addEventListener("click", (event) => {
    const placeButtonHit = event.target.closest("button[data-place-action]");
    if (placeButtonHit) {
      onTeaserPlaceClick(placeButtonHit, event);
      return;
    }
    const button = event.target.closest("button[data-teaser-mark]");
    if (button && !button.disabled) setTeaserMark(button);
  });

  view.viewBets.addEventListener("click", (event) => {
    const button = event.target.closest("button[data-bet-action]");
    if (button && !button.disabled) {
      onBetAction(button);
      return;
    }
    const pick = event.target.closest("[data-attach-pick]");
    if (pick && state.attach) {
      state.attach = { ...state.attach, step: "confirm", eventId: pick.dataset.attachPick, swapped: false, plan: null, error: null };
      renderAll();
      loadAttachPlan();
    }
  });
  view.viewBets.addEventListener("input", (event) => {
    if (!event.target.dataset.attachSearch || !state.attach) return;
    state.attach = { ...state.attach, query: event.target.value };
    clearTimeout(attachSearchTimer);
    attachSearchTimer = setTimeout(loadAttachPicks, ATTACH_SEARCH_DELAY_MS);
  });
  view.betsCrosswalkClear.addEventListener("click", clearCrosswalk);

  view.sheetClose.addEventListener("click", closeTicket);
  view.sheetBackdrop.addEventListener("click", closeTicket);
  document.addEventListener("keydown", (event) => {
    if (event.key === "Escape" && state.ticket) closeTicket();
  });

  view.settingsForm.addEventListener("input", () => {
    state.settingsDirty = true;
    setSettingsStatus("Unsaved changes.", null);
  });
  // The book list appears and goes with the mode, so the mode re-renders the form around what is typed.
  view.settingsForm.addEventListener("change", (event) => {
    if (event.target.dataset.field !== "bookMode") return;
    const parsed = readSettingsForm();
    if (parsed.error) return;
    const values = { ...parsed.values };
    // Picking books starts from the books the list uses now.
    if (values.bookMode === "custom" && !Array.isArray(values.bookIds)) values.bookIds = (state.edges.payload && state.edges.payload.books.ids) || [];
    renderSettings(values);
    state.settingsDirty = true;
    setSettingsStatus("Unsaved changes.", null);
  });
  view.settingsForm.addEventListener("click", (event) => {
    const reset = event.target.closest("[data-reset]");
    if (reset) putSettings(phoneView.resetUpdate(reset.dataset.reset), "Reset to the default.");
  });
  view.settingsSave.addEventListener("click", onSave);

  document.addEventListener("visibilitychange", onVisibility);
  // Relative times ("4s ago", "in 2h 10m") keep moving between polls.
  setInterval(() => {
    if (document.visibilityState === "visible") renderTop(Date.now());
  }, 5000);

  loadSettings();
  onVisibility();
})();
