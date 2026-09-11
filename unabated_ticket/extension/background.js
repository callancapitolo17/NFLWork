// Unabated Ticket — service worker.
//
// Deliberately minimal and OFF the hot path: content.js writes the ticket to
// chrome.storage.local directly and the panel reads it, so a sleeping or
// crashed worker cannot stop a click from reaching the panel. This worker
// sets the panel-open-on-click behavior and handles clicks on edge
// notifications (the panel creates them; a click focuses the Unabated tab
// and points page.js at the row through locate.js, same path as a row click).
//
// Side effects: writes chrome.storage.local {locate, locateResult} via
// locate.js; clears the clicked notification; on install/reload re-injects
// page.js + content.js into Unabated tabs that are already open.

importScripts("locate.js");

// Chrome injects content scripts only into pages loaded after the extension
// (re)loads; a tab that was already open keeps a dead copy whose storage
// writes fail silently, so clicks stop reaching the panel until the tab is
// reloaded. Re-inject instead; the new page.js tells the old one to retire.
const UNABATED_TAB_PATTERN = "https://tools.unabated.com/*";

async function reinjectIntoOpenTabs() {
  const tabs = await chrome.tabs.query({ url: UNABATED_TAB_PATTERN });
  for (const tab of tabs) {
    try {
      await chrome.scripting.executeScript({ target: { tabId: tab.id }, files: ["page.js"], world: "MAIN" });
      await chrome.scripting.executeScript({ target: { tabId: tab.id }, files: ["content.js"], world: "ISOLATED" });
      console.info("[unabated-ticket] re-injected into", tab.url);
    } catch (error) {
      console.error("[unabated-ticket] re-inject failed for", tab.url, error);
    }
  }
}

chrome.runtime.onInstalled.addListener(() => {
  reinjectIntoOpenTabs().catch((error) => console.error("[unabated-ticket] re-inject failed", error));
});

// Clicking the toolbar icon opens the side panel.
chrome.sidePanel
  .setPanelBehavior({ openPanelOnActionClick: true })
  .catch((error) => console.error("[unabated-ticket] setPanelBehavior failed", error));


chrome.notifications.onClicked.addListener(async (notificationId) => {
  chrome.notifications.clear(notificationId);
  try {
    const stored = await chrome.storage.local.get("alertTargets");
    const target = stored.alertTargets && stored.alertTargets[notificationId];
    if (!target) {
      console.info("[unabated-ticket] notification has no stored target (panel closed since?)", notificationId);
      return;
    }
    await globalThis.UnabatedLocate.locateLine(target);
  } catch (error) {
    console.error("[unabated-ticket] notification click failed", error);
  }
});
