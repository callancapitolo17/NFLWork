// Unabated Ticket phone page — fetches ./edges.json from the runner that
// served this page, every 30 s and on returning to the tab, and draws what
// phone_view.js returns. Read-only: no settings, no bet buttons, no writes.
// Every string from the feed goes in through textContent, never parsed as markup.

(function () {
  "use strict";

  const REFRESH_MS = 30 * 1000;
  const view = window.UnabatedPhoneView;
  let lastPayload = null;
  let lastError = null;

  function el(tag, className, text) {
    const node = document.createElement(tag);
    if (className) node.className = className;
    if (text != null) node.textContent = text;
    return node;
  }

  function chipNode(chip) {
    return el("span", `chip tone-${chip.tone}`, chip.text);
  }

  function cardNode(card) {
    const article = el("article", "card");

    const top = el("div", "card-top");
    top.append(el("span", null, card.header), el("span", null, card.when));
    article.append(top);

    const matchup = el("div", "matchup");
    matchup.append(el("span", null, card.matchup));
    for (const badge of card.badges) matchup.append(el("span", `badge tone-${badge.tone}`, badge.text));
    article.append(matchup);

    const main = el("div", "card-main");
    const left = el("div", "side-block");
    left.append(el("div", "side", card.side), el("div", "price mono", card.priceLine));
    if (card.liquidity) left.append(el("div", "muted small", card.liquidity));
    const right = el("div", "edge-block");
    right.append(el("div", `edge mono tier-${card.edgeTier}`, card.edge), el("div", "muted tiny", "edge"));
    if (card.moveLabel) right.append(el("div", "move tiny", card.moveLabel));
    main.append(left, right);
    article.append(main);

    if (card.related.length) {
      const related = el("div", "related");
      for (const text of card.related) related.append(el("div", null, text));
      article.append(related);
    }

    const rail = el("div", "rail");
    const stake = el("div", "stake-block");
    stake.append(el("span", `stake mono${card.atSize ? " at-size" : ""}`, card.stake));
    if (card.stakeNote) stake.append(el("span", "muted small", card.stakeNote));
    rail.append(stake, el("span", "muted small", card.othersText));
    article.append(rail);
    return article;
  }

  function render() {
    const now = Date.now();
    const updated = document.getElementById("updated");
    const chips = document.getElementById("chips");
    const tailFlex = document.getElementById("tail-flex");
    const summary = document.getElementById("summary");
    const list = document.getElementById("cards");
    const errorBox = document.getElementById("error");

    errorBox.hidden = !lastError;
    errorBox.textContent = lastError ? `Can't reach the runner: ${lastError}. Showing the last list it sent.` : "";
    if (!lastPayload) {
      updated.textContent = lastError ? "no data" : "loading…";
      return;
    }
    const page = view.pageView(lastPayload, now);
    updated.textContent = page.updated;
    updated.classList.toggle("tone-bad", page.stale);
    chips.replaceChildren(...page.chips.map(chipNode));
    tailFlex.textContent = page.tailFlex;
    summary.textContent = page.summary;
    if (page.empty) {
      list.replaceChildren(el("p", "muted empty", page.empty));
      return;
    }
    list.replaceChildren(...page.cards.map(cardNode));
  }

  async function refresh() {
    try {
      const response = await fetch("edges.json", { cache: "no-store" });
      if (!response.ok) throw new Error(`HTTP ${response.status}`);
      lastPayload = await response.json();
      lastError = null;
    } catch (error) {
      lastError = error.message;
    }
    render();
  }

  refresh();
  setInterval(refresh, REFRESH_MS);
  // The "updated 12s ago" line ages between fetches.
  setInterval(render, 5 * 1000);
  document.addEventListener("visibilitychange", () => {
    if (document.visibilityState === "visible") refresh();
  });
})();
