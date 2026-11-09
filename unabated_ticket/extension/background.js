// Unabated Ticket — service worker.
//
// Deliberately minimal and OFF the hot path: content.js writes the ticket to
// chrome.storage.local directly and the panel reads it, so a sleeping or
// crashed worker cannot stop a click from reaching the panel. This worker
// sets the panel-open-on-click behavior and handles clicks on edge
// notifications (the panel creates them; a click focuses the Unabated tab
// and points page.js at the row through locate.js, same path as a row click).
//
// It also reads Cal's Underdog slips every 2 minutes (2026-11-09, below):
// in his own app.underdogsports.com tab, on the app's own login, POSTed to the
// bets service — here rather than in the panel so it runs with the panel
// closed (Cal uses the phone page; Chrome on the Mac just has to be open).
//
// Side effects: writes chrome.storage.local {locate, locateResult} via
// locate.js; clears the clicked notification; on install/reload re-injects
// page.js + content.js into Unabated tabs that are already open; POSTs
// /underdog.json to the bets service and may reload an idle Underdog tab.

importScripts("locate.js", "underdog.js");

// Chrome injects content scripts only into pages loaded after the extension
// (re)loads; a tab that was already open keeps a dead copy whose storage
// writes fail silently, so clicks stop reaching the panel until the tab is
// reloaded. Re-inject instead. Both re-injected copies replace their
// predecessor: the new page.js posts a `takeover` the old one retires on
// (they are separate MAIN-world IIFEs), and the new content.js calls the
// retire hook the old one published (they share the isolated world's
// `window`). A guard that let the NEW copy return instead would leave the
// dead one holding the page's messages — the bug this pair of handoffs
// closes; page.js alone was not enough, since content.js is the half that
// writes storage.
const UNABATED_TAB_PATTERN = "https://tools.unabated.com/*";

async function reinjectIntoOpenTabs() {
  for (const tab of await chrome.tabs.query({ url: UNABATED_TAB_PATTERN })) {
    try {
      await chrome.scripting.executeScript({ target: { tabId: tab.id }, files: ["page.js"], world: "MAIN" });
      await chrome.scripting.executeScript({ target: { tabId: tab.id }, files: ["content.js"], world: "ISOLATED" });
      console.info("[unabated-ticket] re-injected into", tab.url);
    } catch (error) {
      console.error("[unabated-ticket] re-inject failed for", tab.url, error);
    }
  }
  // The Underdog reader too, so a tab open before a reload keeps reading
  // without being reloaded (underdogpage.js leaves an earlier copy in place).
  for (const tab of await chrome.tabs.query({ url: underdog.APP_TAB_PATTERN })) {
    if (tab.discarded) continue;
    try {
      await chrome.scripting.executeScript({ target: { tabId: tab.id }, files: ["underdog.js", "underdogpage.js"], world: "MAIN" });
    } catch (error) {
      console.error("[unabated-ticket] underdog re-inject failed for", tab.url, error);
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


// ---- Underdog (read in Cal's tab, stored by the service) ---------------------
//
// On an alarm every underdog.POLL_MINUTES: run the read in the first
// app.underdogsports.com tab (underdogpage.js, MAIN world) and POST what it
// returns to the bets service. The read reuses the headers the app last sent;
// when they have expired (the app refreshes only when it makes a request, so
// an idle tab's 10-minute token lapses) the tab is reloaded so the app logs
// itself back in — at most every RELOAD_MIN_MS and never while Cal is looking
// at it — then read again. A failed read is POSTed as an error so the Bets
// tab's Underdog row turns red with the fix; nothing half-read is pushed.
const underdog = globalThis.UnabatedUnderdog;
const UNDERDOG_ALARM = "underdog-read";
const UNDERDOG_RETRY_EVERY_MS = 3 * 1000;
let underdogBusy = false;

function ensureUnderdogAlarm() {
  chrome.alarms.get(UNDERDOG_ALARM).then((alarm) => {
    if (!alarm) chrome.alarms.create(UNDERDOG_ALARM, { periodInMinutes: underdog.POLL_MINUTES, delayInMinutes: 0.1 });
  }).catch((error) => console.error("[unabated-ticket] underdog: alarm setup failed", error));
}
ensureUnderdogAlarm();
chrome.runtime.onStartup.addListener(ensureUnderdogAlarm);

chrome.alarms.onAlarm.addListener((alarm) => {
  if (alarm.name !== UNDERDOG_ALARM) return;
  pollUnderdog().catch((error) => console.error("[unabated-ticket] underdog poll failed", error));
});

async function postUnderdog(body) {
  const stored = await chrome.storage.local.get("betsSettings");
  const response = await fetch(`${underdog.serviceUrlOf(stored.betsSettings)}${underdog.SERVICE_PATH}`, {
    method: "POST", cache: "no-store", headers: { "Content-Type": "application/json" }, body: JSON.stringify(body),
    signal: AbortSignal.timeout(underdog.FETCH_TIMEOUT_MS),
  });
  const reply = await response.json().catch(() => null);
  if (!response.ok) throw new Error(`bets service answered HTTP ${response.status}${reply && reply.error ? `: ${reply.error}` : ""}`);
  return reply;
}

// What the tab's read returned; {missing: true} when underdogpage.js is not in
// the tab, {unreachable: true} when the tab cannot run a script (discarded by
// Memory Saver, or mid-navigation).
async function readInTab(tab) {
  if (tab.discarded) return { unreachable: true, error: "the Underdog tab was put to sleep by Chrome" };
  try {
    const [injection] = await chrome.scripting.executeScript({
      target: { tabId: tab.id }, world: "MAIN",
      func: () => (typeof globalThis.__unabatedUnderdogRead === "function" ? globalThis.__unabatedUnderdogRead() : { missing: true }),
    });
    return (injection && injection.result) || { error: "the Underdog tab returned nothing" };
  } catch (error) {
    return { unreachable: true, error: `the Underdog tab could not be read: ${error.message}` };
  }
}

async function reloadAndRead(tab) {
  await chrome.storage.session.set({ underdogLastReloadAt: Date.now() });
  await chrome.tabs.reload(tab.id);
  const deadline = Date.now() + underdog.RELOAD_WAIT_MS;
  let result = {};
  while (Date.now() < deadline) {
    await new Promise((resolve) => setTimeout(resolve, UNDERDOG_RETRY_EVERY_MS));
    result = await readInTab(tab);
    if (underdog.stepAfterRead(result) !== "reload") return result;
  }
  return { error: underdog.afterReloadError(result) };
}

async function pollUnderdog() {
  if (underdogBusy) return;
  underdogBusy = true;
  try {
    const tabs = await chrome.tabs.query({ url: underdog.APP_TAB_PATTERN });
    if (!tabs.length) {
      await postUnderdog(underdog.errorBody(underdog.NO_TAB));
      return;
    }
    // The first tab that reads wins; else the first tab is the one to revive.
    let first = null;
    for (const tab of tabs) {
      const result = await readInTab(tab);
      if (result.push) {
        await postUnderdog(result.push);
        return;
      }
      if (!first) first = { tab, result };
    }
    let { result } = first;
    if (underdog.stepAfterRead(result) === "reload") {
      const appWindow = await chrome.windows.get(first.tab.windowId).catch(() => null);
      const { underdogLastReloadAt } = await chrome.storage.session.get("underdogLastReloadAt");
      if (underdog.mayReload(first.tab, Boolean(appWindow && appWindow.focused), underdogLastReloadAt || 0, Date.now())) {
        result = await reloadAndRead(first.tab);
      } else if (result.missing) {
        result = { error: "the Underdog tab was open before the extension loaded — reload it" };
      }
    }
    await postUnderdog(result.push || underdog.errorBody(result.error || underdog.STALE));
  } catch (error) {
    console.warn("[unabated-ticket] underdog:", error.message);
    await postUnderdog(underdog.errorBody(error.message)).catch(() => {});
  } finally {
    underdogBusy = false;
  }
}
