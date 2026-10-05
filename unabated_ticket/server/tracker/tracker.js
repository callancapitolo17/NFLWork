// Bet Tracker page: fetches the bets service's /bets.json (all history),
// builds tickets with trackerstats.js and renders the Overview and Analysis
// views. Read-only: no writes to the service; the only state it keeps is the
// viewer's own display choices (range, $ / units, unit size, view) in
// localStorage, which may be unavailable and is never required.
//
// The service's Content-Security-Policy blocks inline style attributes, so
// every per-element value is set through element.style (the CSSOM), never
// through markup.

(function () {
  "use strict";

  const stats = globalThis.UnabatedTrackerStats;
  const BETS_URL = "/bets.json?days=3650";
  const POLL_MS = 60 * 1000;
  const PREFS_KEY = "betTracker.prefs";
  const RANGES = ["7D", "30D", "90D", "YTD", "All"];
  const RANGE_DAYS = { "7D": 7, "30D": 30, "90D": 90 };
  const KINDS = ["All", "Straight", "Parlay", "Teaser", "Kalshi combo"];
  const DAILY_ROWS = 14;
  const LOG_PAGE = 50;
  const CALENDAR_WEEKS = 6;
  const CI_SCALE = 0.3;
  const CHART_W = 1000;
  const CHART_H = 240;
  const BAR_HALF_PX = 34;
  const CAL_FULL_COLOR_PNL_UNITS = 15;
  const SVG_NS = "http://www.w3.org/2000/svg";
  const COLORS = { pos: "#3dd68c", neg: "#ff6b6b", muted: "#8b96a5", exp: "#6ea8fe", text: "#e7ecf2", warn: "#f5b74f", dim: "#4e5866" };

  const state = Object.assign({
    view: "overview", range: "30D", units: false, unitSize: 100,
    groupBy: "venue", kind: "All", offVenues: [], offLeagues: [], query: "", logLimit: LOG_PAGE,
  }, loadPrefs(), { view: location.hash === "#analysis" ? "analysis" : "overview" });
  let payload = null;
  let tickets = [];

  // ---- prefs ----------------------------------------------------------------

  function loadPrefs() {
    try {
      const saved = JSON.parse(localStorage.getItem(PREFS_KEY) || "{}");
      const prefs = {};
      if (RANGES.includes(saved.range)) prefs.range = saved.range;
      if (typeof saved.units === "boolean") prefs.units = saved.units;
      if (Number.isFinite(saved.unitSize) && saved.unitSize > 0) prefs.unitSize = saved.unitSize;
      return prefs;
    } catch (_error) {
      return {};
    }
  }

  function savePrefs() {
    try {
      localStorage.setItem(PREFS_KEY, JSON.stringify({ range: state.range, units: state.units, unitSize: state.unitSize }));
    } catch (_error) {
      // Private window or blocked storage: the choice lasts until reload.
    }
  }

  // ---- formatting -----------------------------------------------------------

  function sign(value, signed) {
    return value < 0 ? "−" : signed && value > 0 ? "+" : "";
  }

  /** A value that displays as zero carries no sign ("$0", never "−$0"). */
  function money(value, signed) {
    if (state.units) {
      const units = (Math.abs(value) / state.unitSize).toFixed(2);
      return (Number(units) ? sign(value, signed) : "") + units + "u";
    }
    const dollars = Math.round(Math.abs(value));
    return (dollars ? sign(value, signed) : "") + "$" + dollars.toLocaleString("en-US");
  }

  function moneyShort(value) {
    const abs = Math.abs(value);
    if (state.units) return sign(value, true) + (abs / state.unitSize).toFixed(1) + "u";
    if (Math.round(abs) === 0) return "0";
    return sign(value, true) + (abs >= 1000 ? (abs / 1000).toFixed(1) + "k" : String(Math.round(abs)));
  }

  function pct(value, signed) {
    if (value === null || !Number.isFinite(value)) return "—";
    return sign(value, signed) + Math.abs(value * 100).toFixed(1) + "%";
  }

  function american(price) {
    if (!Number.isFinite(price)) return "—";
    return price > 0 ? "+" + price : "−" + Math.abs(price);
  }

  function toneClass(value) {
    return value > 0 ? "pos" : value < 0 ? "neg" : "muted";
  }

  function toneColor(value) {
    return value > 0 ? COLORS.pos : value < 0 ? COLORS.neg : COLORS.muted;
  }

  function dayLabel(dayKey, withWeekday) {
    const options = { month: "short", day: "numeric", timeZone: "UTC" };
    if (withWeekday) options.weekday = "short";
    return new Date(stats.dayKeyToUtc(dayKey)).toLocaleDateString("en-US", options);
  }

  function startLabel(iso) {
    if (!iso) return "—";
    const ms = Date.parse(iso);
    if (!Number.isFinite(ms)) return "—";
    const today = stats.pacificDay(Date.now());
    const day = stats.pacificDay(ms);
    const time = new Date(ms).toLocaleTimeString("en-US", { hour: "numeric", minute: "2-digit", timeZone: stats.PACIFIC_TZ });
    if (day === today) return "Today " + time;
    return new Date(ms).toLocaleDateString("en-US", { weekday: "short", month: "short", day: "numeric", timeZone: stats.PACIFIC_TZ }) + " " + time;
  }

  function ago(iso) {
    const ms = Date.parse(iso);
    if (!Number.isFinite(ms)) return "never";
    const minutes = Math.max(0, Math.round((Date.now() - ms) / 60000));
    if (minutes < 1) return "just now";
    if (minutes < 60) return minutes + " min ago";
    return Math.round(minutes / 60) + " h ago";
  }

  // ---- DOM helpers ----------------------------------------------------------

  function el(tag, options, children) {
    const node = document.createElement(tag);
    const opts = options || {};
    if (opts.className) node.className = opts.className;
    if (opts.text !== undefined) node.textContent = String(opts.text);
    if (opts.title) node.title = opts.title;
    if (opts.attrs) for (const [name, value] of Object.entries(opts.attrs)) node.setAttribute(name, value);
    if (opts.style) Object.assign(node.style, opts.style);
    if (opts.onClick) node.addEventListener("click", opts.onClick);
    for (const child of children || []) if (child) node.appendChild(child);
    return node;
  }

  function svgEl(tag, attrs) {
    const node = document.createElementNS(SVG_NS, tag);
    for (const [name, value] of Object.entries(attrs || {})) node.setAttribute(name, value);
    return node;
  }

  function fill(id, ...nodes) {
    const target = document.getElementById(id);
    target.replaceChildren(...nodes.filter(Boolean));
    return target;
  }

  function setText(id, text) {
    document.getElementById(id).textContent = text;
  }

  function emptyNote(text) {
    return el("div", { className: "empty", text });
  }

  function kpiTiles(id, tiles) {
    fill(id, ...tiles.map((tile) => el("div", { className: "panel kpi" }, [
      el("span", { className: "lbl", text: tile.label }),
      el("span", { className: "value", text: tile.value, style: { color: tile.color || COLORS.text } }),
      el("span", { className: "sub", text: tile.sub }),
    ])));
  }

  /** A table; each column {label, right?, cell(row) -> string | Node, className?(row)}. */
  function table(columns, rows) {
    const head = el("tr", null, columns.map((col) => el("th", { className: col.right ? "r" : "", text: col.label })));
    const body = rows.map((row) => el("tr", null, columns.map((col) => {
      const value = col.cell(row);
      const classes = [col.right ? "r" : "", col.num ? "num" : "", col.className ? col.className(row) : ""].filter(Boolean).join(" ");
      const td = el("td", { className: classes });
      if (value instanceof Node) td.appendChild(value); else td.textContent = value;
      return td;
    })));
    return el("table", { className: "tbl" }, [el("thead", null, [head]), el("tbody", null, body)]);
  }

  function segButtons(id, options, isActive, onPick) {
    fill(id, ...options.map((option) => el("button", {
      text: option.label || option, attrs: { type: "button", "aria-pressed": String(isActive(option)) },
      onClick: () => onPick(option),
    })));
  }

  // ---- range ----------------------------------------------------------------

  function rangeDays() {
    const today = stats.pacificDay(Date.now());
    if (RANGE_DAYS[state.range]) return { first: stats.addDays(today, 1 - RANGE_DAYS[state.range]), last: today };
    if (state.range === "YTD") return { first: today.slice(0, 4) + "-01-01", last: today };
    return { first: stats.firstSettledDay(tickets) || today, last: today };
  }

  function rangeCaption(first, last) {
    return dayLabel(first) + " to " + dayLabel(last) + " · Pacific time";
  }

  function noResultCount(first, last) {
    return tickets.filter((t) => t.status !== "open" && t.pnl === null && t.closedMs !== null
      && stats.pacificDay(t.closedMs) >= first && stats.pacificDay(t.closedMs) <= last).length;
  }

  // ---- overview -------------------------------------------------------------

  function renderOverview() {
    const { first, last } = rangeDays();
    const settled = stats.inDayRange(tickets, first, last);
    const total = stats.summarize(settled);
    const open = tickets.filter((t) => t.status === "open");
    const openStake = open.reduce((sum, t) => sum + t.stake, 0);
    const openWithFair = open.filter((t) => t.expected !== null);
    const openEv = openWithFair.reduce((sum, t) => sum + t.expected, 0);
    const missing = noResultCount(first, last);
    setText("ov-caption", rangeCaption(first, last) + (state.units ? " · 1u = $" + state.unitSize : ""));

    const luck = total.withFair ? total.fairPnl - total.expected : null;
    kpiTiles("ov-kpis", [
      { label: "Net P&L", value: money(total.pnl, true), color: toneColor(total.pnl),
        sub: total.bets + " settled bets" + (missing ? " · " + missing + " without a result" : "") },
      { label: "ROI", value: pct(total.roi, true), color: toneColor(total.pnl), sub: "on " + money(total.handle) + " handle" },
      { label: "Expected P&L", value: total.withFair ? money(total.expected, true) : "—", color: COLORS.exp,
        sub: total.withFair + " of " + total.bets + " bets have a saved fair" },
      { label: "Actual vs expected", value: luck === null ? "—" : money(luck, true), color: luck === null ? COLORS.muted : toneColor(luck),
        sub: total.z === null ? "needs bets with a saved fair" : "z = " + total.z.toFixed(2) + (Math.abs(total.z) < 1.96 ? ", within noise" : ", outside the 95% band") },
      { label: "Record", value: total.wins + "-" + total.losses + "-" + total.pushes, sub: "win rate " + pct(total.winRate) },
      { label: "Open risk", value: money(openStake), sub: open.length + " open bets" + (openWithFair.length ? " · EV " + money(openEv, true) : "") },
    ]);

    const series = stats.dailySeries(tickets, first, last);
    renderChart(series);
    renderCalendar(last);
    renderDaily(series);
    renderVenues(settled);
    renderOpen(open, openStake);
  }

  function renderChart(series) {
    const withBets = series.filter((d) => d.bets > 0);
    if (!withBets.length) {
      fill("ov-stats");
      fill("ov-chart", emptyNote("No settled bets in this range."));
      return;
    }
    const best = withBets.reduce((a, b) => (b.pnl > a.pnl ? b : a));
    const worst = withBets.reduce((a, b) => (b.pnl < a.pnl ? b : a));
    const stat = (label, text, color) => el("div", null, [el("span", { className: "lbl", text: label }), el("span", { className: "num", text, style: { color } })]);
    fill("ov-stats",
      stat("Best day", money(best.pnl, true) + " " + dayLabel(best.day), COLORS.pos),
      stat("Worst day", money(worst.pnl, true) + " " + dayLabel(worst.day), COLORS.neg),
      stat("Up days", withBets.filter((d) => d.pnl > 0).length + " / " + withBets.length, COLORS.text));

    let actual = 0; let expected = 0;
    const points = [{ actual: 0, expected: 0 }].concat(series.map((d) => {
      actual += d.pnl; expected += d.expected;
      return { actual, expected };
    }));
    const values = points.flatMap((p) => [p.actual, p.expected]);
    const top = Math.max(0, ...values) * 1.08 || 1;
    const bottom = Math.min(0, ...values) * 1.08;
    const span = top - bottom || 1;
    const x = (i) => ((i / Math.max(1, points.length - 1)) * CHART_W).toFixed(1);
    const y = (v) => (CHART_H - ((v - bottom) / span) * CHART_H).toFixed(1);
    const path = (key) => points.map((p, i) => (i ? "L" : "M") + x(i) + " " + y(p[key])).join(" ");
    const line = { "vector-effect": "non-scaling-stroke", fill: "none" };
    const svg = svgEl("svg", { viewBox: "0 0 " + CHART_W + " " + CHART_H, preserveAspectRatio: "none", height: String(CHART_H), role: "img", "aria-label": "Cumulative actual and expected profit" });
    svg.append(
      svgEl("path", Object.assign({ d: "M0 0 H1000 M0 80 H1000 M0 160 H1000 M0 240 H1000", stroke: "#1c232d", "stroke-width": "1" }, line)),
      svgEl("path", Object.assign({ d: "M0 " + y(0) + " H1000", stroke: "#3a4452", "stroke-width": "1" }, line)),
      svgEl("path", { d: path("actual") + " L1000 " + y(0) + " L0 " + y(0) + " Z", fill: "rgba(61,214,140,0.09)" }),
      svgEl("path", Object.assign({ d: path("expected"), stroke: COLORS.exp, "stroke-width": "2", "stroke-dasharray": "6 5" }, line)),
      svgEl("path", Object.assign({ d: path("actual"), stroke: COLORS.pos, "stroke-width": "2.25", "stroke-linejoin": "round" }, line)),
    );
    const yAxis = el("div", { className: "yaxis" }, [0, 1, 2, 3].map((i) => el("span", { text: money(top - (span * i) / 3, true) })));
    yAxis.style.height = CHART_H + "px";

    const maxAbs = Math.max(...series.map((d) => Math.abs(d.pnl))) || 1;
    const bars = el("div", { className: "bars" }, series.map((d) => {
      const height = (value) => (value ? Math.max(2, Math.round((Math.abs(value) / maxAbs) * BAR_HALF_PX)) : 0) + "px";
      return el("div", { className: "bar", title: dayLabel(d.day, true) + ": " + money(d.pnl, true) + " on " + d.bets + " bets" }, [
        el("div", { className: "up" }, [el("div", { style: { height: d.pnl > 0 ? height(d.pnl) : "0px" } })]),
        el("div", { className: "down" }, [el("div", { style: { height: d.pnl < 0 ? height(d.pnl) : "0px" } })]),
      ]);
    }));
    const mid = series[Math.floor((series.length - 1) / 2)];
    const xAxis = el("div", { className: "xaxis" }, [series[0], mid, series[series.length - 1]].map((d) => el("span", { text: dayLabel(d.day) })));
    fill("ov-chart", el("div", { className: "plot" }, [yAxis, svg]), bars, xAxis);
  }

  function renderCalendar(today) {
    const todayIndex = stats.WEEKDAYS.indexOf(stats.weekdayOf(today));
    const start = stats.addDays(today, -todayIndex - (CALENDAR_WEEKS - 1) * 7);
    const end = stats.addDays(start, CALENDAR_WEEKS * 7 - 1);
    const byDay = new Map(stats.dailySeries(tickets, start, today).map((d) => [d.day, d]));
    setText("cal-caption", dayLabel(start) + " to " + dayLabel(end));
    const fullColorDollars = CAL_FULL_COLOR_PNL_UNITS * state.unitSize;
    const cells = stats.WEEKDAYS.map((name) => el("span", { className: "lbl wd", text: name }));
    for (let day = start; day <= end; day = stats.addDays(day, 1)) {
      const totals = byDay.get(day);
      const cell = el("div", { className: "cell" + (day === today ? " today" : "") });
      const dayNumber = el("span", { className: "d", text: Number(day.slice(8)) });
      const value = el("span", { className: "v" });
      if (totals && totals.bets) {
        const strength = Math.min(1, Math.abs(totals.pnl) / fullColorDollars);
        const tint = totals.pnl >= 0 ? "61,214,140" : "255,107,107";
        cell.style.background = "rgba(" + tint + "," + (0.1 + strength * 0.55).toFixed(2) + ")";
        value.textContent = moneyShort(totals.pnl);
        value.style.color = strength > 0.55 ? "#0b0e13" : COLORS.text;
        dayNumber.style.color = strength > 0.55 ? "rgba(11,14,19,0.7)" : COLORS.muted;
        cell.title = dayLabel(day, true) + ": " + money(totals.pnl, true) + " on " + totals.bets + " bets";
      } else {
        cell.title = dayLabel(day, true) + (day === today ? ": nothing settled yet" : day > today ? "" : ": no bets settled");
        if (day === today) value.textContent = "Today";
      }
      cell.append(dayNumber, value);
      cells.push(cell);
    }
    fill("ov-calendar", ...cells);
  }

  function renderDaily(series) {
    const rows = series.filter((d) => d.bets > 0).slice(-DAILY_ROWS).reverse();
    if (!rows.length) { fill("ov-daily", emptyNote("No settled bets in this range.")); return; }
    fill("ov-daily", table([
      { label: "Day", cell: (d) => dayLabel(d.day, true) },
      { label: "Bets", right: true, num: true, cell: (d) => String(d.bets) },
      { label: "W-L-P", right: true, num: true, className: () => "muted", cell: (d) => d.wins + "-" + d.losses + "-" + d.pushes },
      { label: "Handle", right: true, num: true, cell: (d) => money(d.handle) },
      { label: "P&L", right: true, num: true, className: (d) => toneClass(d.pnl), cell: (d) => money(d.pnl, true) },
      { label: "ROI", right: true, num: true, className: (d) => toneClass(d.pnl), cell: (d) => pct(d.roi, true) },
      { label: "Expected", right: true, num: true, className: () => "exp", cell: (d) => (d.withFair ? money(d.expected, true) : "—") },
    ], rows));
  }

  function renderVenues(settled) {
    const rows = stats.groupBy(settled, "venue");
    if (!rows.length) { fill("ov-venues", emptyNote("No settled bets in this range.")); return; }
    const maxAbs = Math.max(...rows.map((r) => Math.abs(r.pnl))) || 1;
    fill("ov-venues", el("div", { className: "venues" }, rows.map((row) => {
      const width = (Math.abs(row.pnl) / maxAbs) * 50;
      const bar = el("div", { className: "fill", style: { width: width.toFixed(1) + "%", left: (row.pnl >= 0 ? 50 : 50 - width).toFixed(1) + "%", background: toneColor(row.pnl) } });
      return el("div", { className: "stack tight" }, [
        el("div", { className: "venue-top" }, [el("span", { text: row.label }), el("span", { className: "num " + toneClass(row.pnl), text: money(row.pnl, true) })]),
        el("div", { className: "venue-bottom" }, [
          el("div", { className: "track" }, [el("div", { className: "mid" }), bar]),
          el("span", { className: "num", text: row.bets + " bets · " + money(row.handle) + " · " + pct(row.roi, true) }),
        ]),
      ]);
    })));
  }

  function renderOpen(open, openStake) {
    const toWin = open.reduce((sum, t) => sum + (Number.isFinite(t.toWin) ? t.toWin : 0), 0);
    setText("open-caption", open.length ? "Risking " + money(openStake) + " to win " + money(toWin) + " · fair is Unabated's at fill" : "");
    if (!open.length) { fill("ov-open", emptyNote("No open bets.")); return; }
    const sorted = open.slice().sort((a, b) => (Date.parse(a.eventStart) || Infinity) - (Date.parse(b.eventStart) || Infinity));
    fill("ov-open", table([
      { label: "Starts", className: () => "muted", cell: (t) => startLabel(t.eventStart) },
      { label: "Venue", cell: (t) => t.venue },
      { label: "League", cell: (t) => el("span", { className: "tag", text: t.league }) },
      { label: "Event", cell: (t) => t.event || t.kind },
      { label: "Bet", className: () => "wrap", cell: (t) => t.selection },
      { label: "Price", right: true, num: true, cell: (t) => american(t.displayPrice) },
      { label: "Fair", right: true, num: true, className: () => "exp", cell: (t) => american(t.fairAmerican) },
      { label: "Edge", right: true, num: true, className: (t) => (t.edge === null ? "muted" : toneClass(t.edge)), cell: (t) => pct(t.edge, true) },
      { label: "Stake", right: true, num: true, cell: (t) => money(t.stake) },
      { label: "To win", right: true, num: true, cell: (t) => (Number.isFinite(t.toWin) ? money(t.toWin) : "—") },
    ], sorted));
  }

  // ---- analysis -------------------------------------------------------------

  function filteredSettled() {
    const { first, last } = rangeDays();
    return stats.inDayRange(tickets, first, last).filter((t) => !state.offVenues.includes(t.venue)
      && !state.offLeagues.includes(t.league) && (state.kind === "All" || t.kind === state.kind));
  }

  function valuesByCount(key) {
    const counts = new Map();
    for (const ticket of tickets) counts.set(ticket[key], (counts.get(ticket[key]) || 0) + 1);
    return [...counts].sort((a, b) => b[1] - a[1]).map(([value]) => value);
  }

  function toggleIn(listKey, value) {
    const list = state[listKey];
    state[listKey] = list.includes(value) ? list.filter((v) => v !== value) : list.concat(value);
    state.logLimit = LOG_PAGE;
    render();
  }

  function chips(id, values, offKey) {
    fill(id, ...values.map((value) => el("button", {
      className: "chip", text: value, attrs: { type: "button", "aria-pressed": String(!state[offKey].includes(value)) },
      onClick: () => toggleIn(offKey, value),
    })));
  }

  function renderAnalysis() {
    const { first, last } = rangeDays();
    chips("f-venues", valuesByCount("venue"), "offVenues");
    chips("f-leagues", valuesByCount("league"), "offLeagues");
    segButtons("f-kinds", KINDS, (k) => k === state.kind, (k) => { state.kind = k; state.logLimit = LOG_PAGE; render(); });
    segButtons("group-tabs", stats.GROUPS, (g) => g.key === state.groupBy, (g) => { state.groupBy = g.key; render(); });

    const settled = filteredSettled();
    const total = stats.summarize(settled);
    setText("an-caption", rangeCaption(first, last) + " · " + total.bets + " settled bets");
    kpiTiles("an-kpis", [
      { label: "Bets", value: String(total.bets), sub: total.wins + "-" + total.losses + "-" + total.pushes + " W-L-P" },
      { label: "Handle", value: money(total.handle), sub: "avg stake " + money(total.bets ? total.handle / total.bets : 0) },
      { label: "P&L", value: money(total.pnl, true), color: toneColor(total.pnl), sub: total.withFair ? "expected " + money(total.expected, true) + " on " + total.withFair + " with a fair" : "no saved fairs" },
      { label: "ROI", value: pct(total.roi, true), color: toneColor(total.roi), sub: "95%: " + pct(total.roi - total.ciHalf, true) + " to " + pct(total.roi + total.ciHalf, true) },
      { label: "Avg edge at fill", value: pct(total.expRoi, true), color: COLORS.exp, sub: "stake-weighted, Unabated fair" },
      { label: "Luck (z)", value: total.z === null ? "—" : total.z.toFixed(2), color: total.z !== null && Math.abs(total.z) >= 1.96 ? COLORS.warn : COLORS.text,
        sub: total.z === null ? "needs bets with a saved fair" : Math.abs(total.z) < 1.96 ? "within normal variance" : "outside the 95% band" },
    ]);
    renderGroups(settled);
    renderCalibration(settled);
    renderLog(settled);
  }

  function ciBar(row) {
    const scale = (roi) => Math.max(0, Math.min(100, 50 + (roi / CI_SCALE) * 50));
    const low = scale(row.roi - row.ciHalf); const high = scale(row.roi + row.ciHalf);
    const proven = row.roi - row.ciHalf > 0 || row.roi + row.ciHalf < 0;
    const color = proven ? toneColor(row.roi) : COLORS.dim;
    return el("div", { className: "ci", title: "95% interval " + pct(row.roi - row.ciHalf, true) + " to " + pct(row.roi + row.ciHalf, true) }, [
      el("div", { className: "base" }), el("div", { className: "zero" }),
      el("div", { className: "range", style: { left: low + "%", width: (high - low) + "%", background: color } }),
      el("div", { className: "cap", style: { left: low + "%", background: color } }),
      el("div", { className: "cap", style: { left: "calc(" + high + "% - 2px)", background: color } }),
      el("div", { className: "pt", style: { left: scale(row.roi) + "%", background: toneColor(row.roi) } }),
    ]);
  }

  function renderGroups(settled) {
    const rows = stats.groupBy(settled, state.groupBy);
    if (!rows.length) { fill("an-groups", emptyNote("No settled bets match these filters.")); return; }
    const label = stats.GROUPS.find((g) => g.key === state.groupBy).label;
    fill("an-groups", table([
      { label, cell: (r) => r.label },
      { label: "Bets", right: true, num: true, cell: (r) => String(r.bets) },
      { label: "W-L-P", right: true, num: true, className: () => "muted", cell: (r) => r.wins + "-" + r.losses + "-" + r.pushes },
      { label: "Handle", right: true, num: true, cell: (r) => money(r.handle) },
      { label: "P&L", right: true, num: true, className: (r) => toneClass(r.pnl), cell: (r) => money(r.pnl, true) },
      { label: "ROI", right: true, num: true, className: (r) => toneClass(r.pnl), cell: (r) => pct(r.roi, true) },
      { label: "ROI, 95% interval", cell: ciBar },
      { label: "Expected ROI", right: true, num: true, className: () => "exp", cell: (r) => pct(r.expRoi, true) },
      { label: "vs expected", right: true, num: true, className: (r) => (r.withFair ? toneClass(r.fairPnl - r.expected) : "muted"), cell: (r) => (r.withFair ? money(r.fairPnl - r.expected, true) : "—") },
      { label: "z", right: true, num: true, className: (r) => (r.z !== null && Math.abs(r.z) >= 1.96 ? "warn" : "muted"), cell: (r) => (r.z === null ? "—" : r.z.toFixed(2)) },
    ], rows));
  }

  function renderCalibration(settled) {
    const bins = stats.calibration(settled);
    if (!bins.length) {
      fill("an-calibration", emptyNote("No settled straight bets with a saved fair yet."));
      fill("an-cal-table");
      return;
    }
    // One square scale for both axes, wide enough for every bin, its win rate and its interval.
    const extents = bins.flatMap((b) => [b.low, b.high, b.actual - b.ciHalf, b.actual + b.ciHalf]);
    const low = Math.max(0, Math.floor(Math.min(...extents) * 10) / 10);
    const high = Math.min(1, Math.ceil(Math.max(...extents) * 10) / 10);
    const span = high - low || 1;
    const x = (p) => (((p - low) / span) * 400).toFixed(1);
    const y = (p) => (300 - ((Math.max(low, Math.min(high, p)) - low) / span) * 300).toFixed(1);
    const line = { "vector-effect": "non-scaling-stroke", fill: "none" };
    const svg = svgEl("svg", { viewBox: "0 0 400 300", role: "img", "aria-label": "Fair win probability against actual win rate" });
    svg.append(
      svgEl("path", Object.assign({ d: "M0 0 H400 M0 150 H400 M0 300 H400 M0 0 V300 M200 0 V300 M400 0 V300", stroke: "#1c232d", "stroke-width": "1" }, line)),
      svgEl("path", Object.assign({ d: "M0 300 L400 0", stroke: "#3a4452", "stroke-width": "1.5", "stroke-dasharray": "5 5" }, line)),
      svgEl("path", Object.assign({ d: bins.map((b) => "M" + x(b.expected) + " " + y(b.actual + b.ciHalf) + " V" + y(b.actual - b.ciHalf)).join(" "), stroke: COLORS.exp, "stroke-width": "1.5" }, line)),
    );
    for (const bin of bins) {
      const dot = svgEl("circle", { cx: x(bin.expected), cy: y(bin.actual), r: "6", fill: COLORS.exp, stroke: "#12161d", "stroke-width": "2", "vector-effect": "non-scaling-stroke" });
      const tip = svgEl("title");
      tip.textContent = bin.bets + " bets: expected " + pct(bin.expected) + ", won " + pct(bin.actual);
      dot.appendChild(tip);
      svg.appendChild(dot);
    }
    const axisLabels = [high, (high + low) / 2, low].map((p) => el("span", { text: pct(p) }));
    // Square-ish plot kept at the viewBox's aspect so the dots stay round; the
    // y labels stretch to whatever height that gives.
    const yAxis = el("div", { className: "yaxis" }, axisLabels);
    const xAxis = el("div", { className: "xaxis" }, [low, (high + low) / 2, high].map((p) => el("span", { text: pct(p) })));
    fill("an-calibration", el("div", { className: "cal-plot" }, [el("div", { className: "plot" }, [yAxis, svg]), xAxis]),
      el("div", { className: "muted", text: "Fair win probability at fill (x) against actual win rate (y); bars are 95% intervals." }));
    fill("an-cal-table", table([
      { label: "Fair prob.", cell: (b) => Math.round(b.low * 100) + " to " + Math.round(b.high * 100) + "%" },
      { label: "Bets", right: true, num: true, cell: (b) => String(b.bets) },
      { label: "Expected win", right: true, num: true, className: () => "exp", cell: (b) => pct(b.expected) },
      { label: "Actual win", right: true, num: true, cell: (b) => pct(b.actual) },
      { label: "Diff", right: true, num: true, className: (b) => (Math.abs(b.actual - b.expected) > b.ciHalf ? "warn" : "muted"), cell: (b) => pct(b.actual - b.expected, true) },
    ], bins));
  }

  function renderLog(settled) {
    const query = state.query.trim().toLowerCase();
    const matching = query
      ? settled.filter((t) => [t.event, t.selection, t.venue, t.league, t.kind].join(" ").toLowerCase().includes(query))
      : settled;
    const shown = matching.slice(0, state.logLimit);
    setText("log-caption", "Showing " + shown.length + " of " + matching.length + " settled bets, newest first");
    document.getElementById("log-more").hidden = matching.length <= shown.length;
    if (!shown.length) { fill("an-log", emptyNote(query ? "No bets match that search." : "No settled bets match these filters.")); return; }
    fill("an-log", table([
      { label: "Settled", className: () => "muted", cell: (t) => dayLabel(t.settledDay, true) },
      { label: "Venue", cell: (t) => t.venue },
      { label: "League", cell: (t) => el("span", { className: "tag", text: t.league }) },
      { label: "Event", cell: (t) => t.event || t.kind },
      { label: "Bet", className: () => "wrap", cell: (t) => t.selection },
      { label: "Price", right: true, num: true, cell: (t) => american(t.displayPrice) },
      { label: "Fair", right: true, num: true, className: () => "exp", cell: (t) => american(t.fairAmerican) },
      { label: "Edge", right: true, num: true, className: (t) => (t.edge === null ? "muted" : toneClass(t.edge)), cell: (t) => pct(t.edge, true) },
      { label: "Stake", right: true, num: true, cell: (t) => money(t.stake) },
      { label: "Result", cell: (t) => el("span", { className: "result " + t.status, text: t.status[0].toUpperCase() + t.status.slice(1) }) },
      { label: "P&L", right: true, num: true, className: (t) => toneClass(t.pnl), cell: (t) => money(t.pnl, true) },
    ], shown));
  }

  // ---- shell ----------------------------------------------------------------

  function renderHeader() {
    for (const button of document.querySelectorAll(".nav button")) {
      if (button.dataset.view === state.view) button.setAttribute("aria-current", "page");
      else button.removeAttribute("aria-current");
    }
    document.getElementById("view-overview").hidden = state.view !== "overview";
    document.getElementById("view-analysis").hidden = state.view !== "analysis";
    segButtons("ranges", RANGES, (r) => r === state.range, (r) => { state.range = r; state.logLimit = LOG_PAGE; savePrefs(); render(); });
    document.getElementById("show-dollars").setAttribute("aria-pressed", String(!state.units));
    document.getElementById("show-units").setAttribute("aria-pressed", String(state.units));
    document.getElementById("unit-size-label").hidden = !state.units;
    const unitInput = document.getElementById("unit-size");
    if (document.activeElement !== unitInput) unitInput.value = String(state.unitSize);
  }

  function renderSync() {
    const dot = document.getElementById("sync-dot");
    const sources = Object.entries((payload && payload.sources) || {});
    if (!sources.length) { setText("sync-text", payload ? "No venues reporting" : "Loading…"); return; }
    const failing = sources.filter(([, s]) => !s.ok);
    dot.classList.toggle("bad", failing.length > 0);
    const okCount = sources.length - failing.length;
    setText("sync-text", okCount + " of " + sources.length + " venues ok · " + ago(payload.generatedAt));
    document.getElementById("sync").title = failing.length
      ? failing.map(([name, s]) => stats.venueName(name) + ": " + (s.error || "no completed poll yet")).join("\n")
      : "Every venue's last poll succeeded";
  }

  function render() {
    renderHeader();
    renderSync();
    if (!payload) return;
    if (state.view === "overview") renderOverview(); else renderAnalysis();
  }

  function showBanner(text) {
    const banner = document.getElementById("banner");
    banner.hidden = !text;
    banner.textContent = text || "";
  }

  async function refresh() {
    try {
      const response = await fetch(BETS_URL, { cache: "no-store" });
      if (!response.ok) throw new Error("the bets service answered HTTP " + response.status);
      payload = await response.json();
      tickets = stats.buildTickets(payload.bets, payload.fillFairs);
      showBanner(null);
    } catch (error) {
      showBanner("Could not load bets: " + error.message + (payload ? ". Showing the last load." : "."));
    }
    render();
  }

  function wire() {
    for (const button of document.querySelectorAll(".nav button")) {
      button.addEventListener("click", () => {
        state.view = button.dataset.view;
        history.replaceState(null, "", state.view === "analysis" ? "#analysis" : "#overview");
        render();
      });
    }
    document.getElementById("show-dollars").addEventListener("click", () => { state.units = false; savePrefs(); render(); });
    document.getElementById("show-units").addEventListener("click", () => { state.units = true; savePrefs(); render(); });
    document.getElementById("unit-size").addEventListener("change", (event) => {
      const size = Number(event.target.value);
      if (Number.isFinite(size) && size > 0) { state.unitSize = size; savePrefs(); }
      render();
    });
    document.getElementById("log-search").addEventListener("input", (event) => {
      state.query = event.target.value; state.logLimit = LOG_PAGE; render();
    });
    document.getElementById("log-more").addEventListener("click", () => { state.logLimit += LOG_PAGE; render(); });
    window.addEventListener("hashchange", () => { state.view = location.hash === "#analysis" ? "analysis" : "overview"; render(); });
    document.addEventListener("visibilitychange", () => { if (!document.hidden) refresh(); });
    setInterval(() => { if (!document.hidden) refresh(); }, POLL_MS);
  }

  wire();
  render();
  refresh();
})();
