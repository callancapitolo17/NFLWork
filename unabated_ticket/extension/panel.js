// Unabated Ticket — side panel.
//
// Reads: chrome.storage.session {ticket, error, watchStatus} (written by
// background.js) and chrome.storage.local {bankroll, multiplier}.
// Writes: chrome.storage.local settings only. Re-renders on storage.onChanged.

(function () {
  "use strict";

  const DEFAULT_SETTINGS = { bankroll: 30000, multiplier: 0.25 };
  // Watcher heartbeats every 5s; past this with no heartbeat, the Unabated tab is gone.
  const WATCH_STALE_MS = 15000;
  const kelly = globalThis.UnabatedKelly;

  const el = (id) => document.getElementById(id);
  const view = {
    ticket: el("ticket"), error: el("error"), empty: el("empty"),
    warning: el("warning"), sideLabel: el("side-label"), betLine: el("bet-line"),
    eventLine: el("event-line"), startLine: el("start-line"),
    book: el("book"), price: el("price"), fair: el("fair"), edge: el("edge"),
    stake: el("stake"), fullKelly: el("full-kelly"),
    copy: el("copy"), copyStatus: el("copy-status"), errorDetail: el("error-detail"),
    bankroll: el("bankroll"), multiplier: el("multiplier"), settingsError: el("settings-error"),
  };

  let state = { ticket: null, error: null, watchStatus: null, settings: { ...DEFAULT_SETTINGS } };
  let lastCopyText = "";

  // ---- formatting ----------------------------------------------------------

  function fmtAmerican(price) {
    return price > 0 ? `+${price}` : `${price}`;
  }

  function fmtDollars(value) {
    return value.toLocaleString("en-US", { style: "currency", currency: "USD", minimumFractionDigits: 2, maximumFractionDigits: 2 });
  }

  function fmtPct(fraction) {
    const pct = fraction * 100;
    return `${pct > 0 ? "+" : ""}${pct.toFixed(2)}%`;
  }

  function fmtStart(eventStart) {
    if (!eventStart) return "";
    // Unabated's eventStart is naive UTC ("2026-09-12T23:30:00").
    const date = new Date(eventStart.endsWith("Z") ? eventStart : `${eventStart}Z`);
    if (Number.isNaN(date.getTime())) return eventStart;
    return date.toLocaleString([], { weekday: "short", month: "short", day: "numeric", hour: "numeric", minute: "2-digit" });
  }

  // Edge from our fair/price math; if Unabated's own edge % (display only) disagrees, say so.
  const EDGE_MISMATCH_PCT = 0.05;
  function fmtEdge(edgeFraction, ticket) {
    const ours = fmtPct(edgeFraction);
    if (ticket.current || ticket.edgePct == null) return ours;
    const unabated = ticket.edgePct;
    if (Math.abs(unabated - edgeFraction * 100) <= EDGE_MISMATCH_PCT) return ours;
    return `${ours} (Unabated shows ${unabated > 0 ? "+" : ""}${unabated.toFixed(2)}%)`;
  }

  function fmtPoints(points) {
    return points == null ? "" : `${points > 0 ? "+" : ""}${points}`;
  }

  // ---- pricing -------------------------------------------------------------

  // The line the stake is computed from: the current line if it moved, else the captured one.
  function pricedLine(ticket) {
    const current = ticket.current;
    if (!current) return { price: ticket.price, fair: ticket.fair, points: ticket.points, moved: false };
    return { price: current.price, fair: current.fair, points: current.points, moved: true };
  }

  function computeStake(ticket, settings) {
    const line = pricedLine(ticket);
    if (line.fair == null) return { line, result: null, reason: "no Unabated fair at the new line" };
    try {
      const result = kelly.kellyStake({ bookPrice: line.price, fairPrice: line.fair, bankroll: settings.bankroll, multiplier: settings.multiplier });
      return { line, result, reason: null };
    } catch (error) {
      return { line, result: null, reason: error.message };
    }
  }

  // ---- rendering -----------------------------------------------------------

  function show(which) {
    view.ticket.hidden = which !== "ticket";
    view.error.hidden = which !== "error";
    view.empty.hidden = which !== "empty";
  }

  function watcherIsLive(ticket, watchStatus) {
    if (!watchStatus || watchStatus.capturedAt !== ticket.capturedAt) {
      // No heartbeat yet: live only during the first interval after capture.
      return Date.now() - ticket.capturedAt < WATCH_STALE_MS;
    }
    return !watchStatus.error && Date.now() - watchStatus.seenAt < WATCH_STALE_MS;
  }

  function renderWarning(ticket, watchStatus) {
    const messages = [];
    let bad = false;
    if (ticket.current && ticket.current.offBoard) {
      messages.push("Off the board at this book.");
      bad = true;
    } else if (ticket.current) {
      const pts = ticket.current.points != null ? ` at ${fmtPoints(ticket.current.points)}` : "";
      messages.push(`Line moved: now ${fmtAmerican(ticket.current.price)}${pts} (captured ${fmtAmerican(ticket.price)}${ticket.points != null ? ` at ${fmtPoints(ticket.points)}` : ""}). Stake re-sized.`);
    }
    if (!watcherIsLive(ticket, watchStatus)) {
      const why = watchStatus && watchStatus.error ? `: ${watchStatus.error}` : "";
      messages.push(`Not watching the line${why}. Showing the captured price.`);
    }
    view.warning.hidden = messages.length === 0;
    view.warning.textContent = messages.join(" ");
    view.warning.classList.toggle("bad", bad);
  }

  function renderTicket() {
    const { ticket, settings, watchStatus } = state;
    renderWarning(ticket, watchStatus);

    view.sideLabel.textContent = ticket.sideLabel;
    const rotation = ticket.rotation != null ? ` · rot ${ticket.rotation}` : "";
    view.betLine.textContent = `${ticket.betType}${ticket.points != null ? ` ${fmtPoints(ticket.points)}` : ""} · ${(ticket.league || "").toUpperCase()}${rotation}`;
    view.eventLine.textContent = ticket.eventName || "";
    view.startLine.textContent = fmtStart(ticket.eventStart);

    const { line, result, reason } = computeStake(ticket, settings);
    view.book.textContent = ticket.book.name;
    view.price.textContent = fmtAmerican(line.price);
    view.fair.textContent = line.fair == null ? "unknown" : fmtAmerican(line.fair);
    view.edge.textContent = result ? fmtEdge(result.edge, ticket) : "—";

    view.stake.classList.remove("no-edge");
    if (!result) {
      view.stake.textContent = "—";
      view.stake.classList.add("no-edge");
      view.fullKelly.textContent = `Cannot size: ${reason}`;
    } else if (result.stake <= 0) {
      view.stake.textContent = fmtDollars(0);
      view.stake.classList.add("no-edge");
      view.fullKelly.textContent = "No edge at this price.";
    } else {
      view.stake.textContent = fmtDollars(result.stake);
      view.fullKelly.textContent = `Full Kelly ${fmtDollars(result.fullKellyStake)} (${fmtPct(result.fullKellyFraction)} of bankroll) × ${settings.multiplier}`;
    }

    const stakeText = result ? result.stake.toFixed(2) : "n/a";
    lastCopyText = `${ticket.sideLabel} ${fmtAmerican(line.price)} @ ${ticket.book.name} | fair ${line.fair == null ? "?" : fmtAmerican(line.fair)} | edge ${result ? fmtPct(result.edge) : "?"} | stake $${stakeText} | ${ticket.eventName || ""}`;
    view.copyStatus.textContent = "";
    show("ticket");
  }

  function render() {
    if (state.error) {
      view.errorDetail.textContent = `${state.error.message} (${new Date(state.error.at).toLocaleTimeString()})`;
      show("error");
      return;
    }
    if (!state.ticket) {
      show("empty");
      return;
    }
    renderTicket();
  }

  // ---- settings ------------------------------------------------------------

  function readSettingInputs() {
    const bankroll = Number(view.bankroll.value);
    const multiplier = Number(view.multiplier.value);
    if (!Number.isFinite(bankroll) || bankroll <= 0) return { error: "Bankroll must be a positive number." };
    if (!Number.isFinite(multiplier) || multiplier <= 0 || multiplier > 1) return { error: "Multiplier must be between 0 and 1." };
    return { settings: { bankroll, multiplier } };
  }

  function onSettingsInput() {
    const parsed = readSettingInputs();
    view.settingsError.textContent = parsed.error || "";
    if (parsed.error) return;
    state.settings = parsed.settings;
    chrome.storage.local.set(parsed.settings);
    render();
  }

  function fillSettingInputs() {
    view.bankroll.value = state.settings.bankroll;
    view.multiplier.value = state.settings.multiplier;
  }

  // ---- wiring --------------------------------------------------------------

  async function load() {
    const local = await chrome.storage.local.get(DEFAULT_SETTINGS);
    state.settings = { bankroll: Number(local.bankroll) || DEFAULT_SETTINGS.bankroll, multiplier: Number(local.multiplier) || DEFAULT_SETTINGS.multiplier };
    fillSettingInputs();
    const session = await chrome.storage.session.get(["ticket", "error", "watchStatus"]);
    state.ticket = session.ticket || null;
    state.error = session.error || null;
    state.watchStatus = session.watchStatus || null;
    render();
  }

  chrome.storage.onChanged.addListener((changes, area) => {
    if (area === "session") {
      if ("ticket" in changes) state.ticket = changes.ticket.newValue || null;
      if ("error" in changes) state.error = changes.error.newValue || null;
      if ("watchStatus" in changes) state.watchStatus = changes.watchStatus.newValue || null;
      render();
    }
  });

  view.bankroll.addEventListener("input", onSettingsInput);
  view.multiplier.addEventListener("input", onSettingsInput);

  view.copy.addEventListener("click", async () => {
    try {
      await navigator.clipboard.writeText(lastCopyText);
      view.copyStatus.textContent = "Copied";
    } catch (error) {
      view.copyStatus.textContent = `Copy failed: ${error.message}`;
    }
  });

  // Re-evaluate the "not watching" state even when no storage event arrives.
  setInterval(() => { if (state.ticket && !state.error) renderTicket(); }, 5000);

  load().catch((error) => {
    view.errorDetail.textContent = error.message;
    show("error");
  });
})();
