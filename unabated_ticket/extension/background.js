// Unabated Ticket — service worker.
//
// Owns every write to chrome.storage.session (the one live ticket) and opens
// the side panel. The panel only reads storage and settings.
//
// Session keys: `ticket` (see the plan's ticket contract), `error`
// ({message, at}) when the last click could not be read, `watchStatus`
// ({capturedAt, seenAt, error}) heartbeat from the line watcher.

const PANEL_ON_ACTION_CLICK = { openPanelOnActionClick: true };

// Runs on every service-worker start, which covers install and update.
chrome.sidePanel.setPanelBehavior(PANEL_ON_ACTION_CLICK).catch(() => {});

async function openPanelForSender(sender) {
  if (!sender.tab || sender.tab.windowId == null) return;
  try {
    // windowId (not tabId) so the panel stays open when the book tab takes focus.
    await chrome.sidePanel.open({ windowId: sender.tab.windowId });
  } catch (_error) {
    // Chrome refuses when no user gesture is attached; the toolbar button still opens it.
  }
}

function linesDiffer(ticket, line) {
  return ticket.price !== line.price || ticket.points !== line.points;
}

function sameCurrent(a, b) {
  if (a == null || b == null) return a == null && b == null;
  return a.price === b.price && a.points === b.points && a.fair === b.fair && a.offBoard === b.offBoard;
}

async function handleTicket(ticket, sender) {
  await chrome.storage.session.set({ ticket, error: null, watchStatus: null });
  await openPanelForSender(sender);
}

async function handleError(payload, sender) {
  await chrome.storage.session.set({ ticket: null, error: payload, watchStatus: null });
  await openPanelForSender(sender);
}

async function handleReady(payload) {
  await chrome.storage.session.set({ pageReady: { url: payload.url, at: payload.at } });
}

async function handleWatch(payload) {
  const { ticket } = await chrome.storage.session.get("ticket");
  if (!ticket || ticket.capturedAt !== payload.capturedAt) return; // stale watcher
  if (payload.error) {
    await chrome.storage.session.set({ watchStatus: { capturedAt: ticket.capturedAt, seenAt: Date.now(), error: payload.error } });
    return;
  }
  const line = payload.line;
  const current = linesDiffer(ticket, line) || line.offBoard ? line : null;
  const updates = { watchStatus: { capturedAt: ticket.capturedAt, seenAt: line.seenAt, error: null } };
  // Only rewrite the ticket when `current` actually changes, so the panel does not re-render every 5s.
  if (!sameCurrent(ticket.current, current)) updates.ticket = { ...ticket, current };
  await chrome.storage.session.set(updates);
}

chrome.runtime.onMessage.addListener((message, sender, sendResponse) => {
  if (!message || typeof message.type !== "string") return false;
  console.info("[unabated-ticket] message", message.type, "from", sender.tab ? sender.tab.url : "no tab");
  const handlers = { ticket: handleTicket, error: handleError, watch: handleWatch, ready: handleReady };
  const handler = handlers[message.type];
  if (!handler) return false;
  handler(message.payload, sender)
    .catch((error) => console.error("[unabated-ticket]", message.type, error))
    .finally(() => sendResponse({ ok: true }));
  return true;
});
