// Unabated Ticket — "take me to this line" plumbing shared by the side panel
// (row click) and the service worker (notification click).
//
// Writes chrome.storage.local {locate} — the request content.js forwards to
// page.js, which scrolls the grid to the row and flashes the cell — then
// focuses the Unabated tab for that league, navigating or opening one when
// none shows it. The click on the price itself stays yours, so the book's
// deeplink is a real user gesture and never popup-blocked.
//
// Loaded as a plain <script> in panel.html and via importScripts in
// background.js (globalThis.UnabatedLocate).

(function (root) {
  "use strict";

  const ODDS_URL = (league) => `https://tools.unabated.com/${league}/odds`;
  const UNABATED_TAB_PATTERN = "https://tools.unabated.com/*";

  function tabShowsLeague(tab, league) {
    try {
      return new URL(tab.url || "").pathname.startsWith(`/${league}/`);
    } catch (_error) {
      return false;
    }
  }

  // Store the request first so a tab that has to (re)load finds it on start.
  async function locateLine(request) {
    const locate = { ...request, at: Date.now() };
    await chrome.storage.local.set({ locate, locateResult: null });
    const tabs = await chrome.tabs.query({ url: UNABATED_TAB_PATTERN });
    const showing = tabs.find((tab) => tabShowsLeague(tab, request.league));
    const tab = showing || tabs[0];
    if (!tab) {
      const created = await chrome.tabs.create({ url: ODDS_URL(request.league), active: true });
      return { tabId: created.id, navigated: true, opened: true };
    }
    const update = { active: true };
    if (!showing) update.url = ODDS_URL(request.league);
    await chrome.tabs.update(tab.id, update);
    if (tab.windowId != null) await chrome.windows.update(tab.windowId, { focused: true });
    return { tabId: tab.id, navigated: !showing, opened: false };
  }

  const api = { locateLine, ODDS_URL, UNABATED_TAB_PATTERN };
  if (typeof module !== "undefined" && module.exports) {
    module.exports = api;
  } else {
    root.UnabatedLocate = api;
  }
})(typeof globalThis !== "undefined" ? globalThis : this);
