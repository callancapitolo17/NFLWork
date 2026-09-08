// Unabated Ticket — service worker.
//
// Deliberately minimal and OFF the hot path: content.js writes the ticket to
// chrome.storage.session directly and the panel reads it, so a sleeping or
// crashed worker cannot stop a click from reaching the panel. This worker only
// runs two one-time setup calls on each start.

// Clicking the toolbar icon opens the side panel.
chrome.sidePanel
  .setPanelBehavior({ openPanelOnActionClick: true })
  .catch((error) => console.error("[unabated-ticket] setPanelBehavior failed", error));

// content.js is an untrusted context; session storage is TRUSTED_CONTEXTS-only
// by default, so open it to content scripts or their writes silently no-op.
chrome.storage.session
  .setAccessLevel({ accessLevel: "TRUSTED_AND_UNTRUSTED_CONTEXTS" })
  .catch((error) => console.error("[unabated-ticket] setAccessLevel failed", error));
