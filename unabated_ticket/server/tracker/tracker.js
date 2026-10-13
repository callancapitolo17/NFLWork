// Bet Tracker page: fetches the bets service's /bets.json (all history) and
// /scenarios.json (the server runner's live game cards, passed through by the
// bets service), builds tickets with trackerstats.js and renders the
// Overview, Live, Open, Analysis and Bets views. Its one write is the Bets view's Remove / Restore: POST
// /exclusions.json, which marks a BFA or Wagerzon bet as not Cal's
// (bets.duckdb::bet_exclusions; the bet itself stays) and takes it out of
// every number here. Otherwise it keeps only the viewer's own display choices
// (range, $ / units, unit size) in localStorage, which may be unavailable and
// is never required.
//
// The service's Content-Security-Policy blocks inline style attributes, so
// every per-element value is set through element.style (the CSSOM), never
// through markup.

(function () {
  "use strict";

  const stats = globalThis.UnabatedTrackerStats;
  const BETS_URL = "/bets.json?days=3650";
  const EXCLUSIONS_URL = "/exclusions.json";
  const SCENARIOS_URL = "/scenarios.json";
  const VIEWS = ["overview", "live", "open", "analysis", "bets"];
  // Cal's ask (2026-10-05): only these books carry bets that are not his (service.EXCLUDABLE_VENUES).
  const REMOVABLE_VENUES = ["BFA", "Wagerzon"];
  const BET_FILTERS = ["All", "Counted", "Removed"];
  const POLL_MS = 60 * 1000;
  const PREFS_KEY = "betTracker.prefs";
  const RANGES = ["Today", "Yesterday", "7D", "30D", "90D", "YTD", "All", "Custom"];
  // The Analysis filter bar, one chip row per entry (Cal, 2026-10-12): every
  // value a ticket carries is a chip, all on by default; state.off[key] holds
  // the ones switched off, so a value new to the data shows up on. Keys match
  // trackerstats.GROUPS so a clicked breakdown row or pivot cell becomes a
  // filter (Cal, 2026-10-13: drill down like a BI tool); `fromGroup` maps a
  // group label to the filter's value where the two differ. Rows without
  // `primary` stay folded under More filters until one is in use.
  // A ticket's value for a group key, labelled exactly as trackerstats.groupBy
  // labels it (missing -> "Unknown"), so a clicked "Unknown" row filters to it.
  function groupValue(key) {
    return (ticket) => ticket[key] || "Unknown";
  }
  const FILTERS = [
    { key: "venue", label: "Venue", primary: true, value: groupValue("venue") },
    { key: "league", label: "League", primary: true, value: groupValue("league") },
    { key: "kind", label: "Type", primary: true, value: groupValue("kind") },
    { key: "market", label: "Market", primary: true, value: groupValue("market") },
    { key: "period", label: "Period", primary: true, value: (t) => periodFilterValue(t.period), fromGroup: periodFilterValue },
    { key: "side", label: "Side", primary: true, value: groupValue("side") },
    { key: "closeResult", label: "Close", primary: true, value: groupValue("closeResult") },
    { key: "favDog", label: "Fav / Dog", value: groupValue("favDog") },
    { key: "oddsBucket", label: "Odds", value: groupValue("oddsBucket") },
    { key: "edgeBucket", label: "Edge", value: groupValue("edgeBucket") },
    { key: "timing", label: "Timing", value: groupValue("timing") },
    { key: "weekday", label: "Weekday", value: groupValue("weekday") },
    { key: "stakeBucket", label: "Stake", value: groupValue("stakeBucket") },
    { key: "number", label: "Number", value: groupValue("number") },
  ];
  const PIVOT_METRICS = [
    { key: "roi", label: "ROI" }, { key: "pnl", label: "P&L" }, { key: "clv", label: "CLV" }, { key: "bets", label: "Bets" },
  ];
  // ROI at which a pivot cell is fully shaded; smaller ROIs fade toward clear.
  const PIVOT_FULL_SHADE_ROI = 0.3;
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
    view: "overview", range: "30D", customFirst: null, customLast: null, units: false, unitSize: 100,
    groupBy: "venue", off: Object.fromEntries(FILTERS.map((f) => [f.key, []])), query: "", logLimit: LOG_PAGE,
    // The drill path: filter keys in the order they came into use, so the
    // breadcrumb reads like the clicks that built it. oddsMin/oddsMax are the
    // typed American odds range (null = open end).
    drillOrder: [], oddsMin: null, oddsMax: null, moreFilters: false,
    pivotRows: "side", pivotCols: "timing", pivotMetric: "roi",
    betsQuery: "", betsVenue: "All", betsFilter: "All", betsLimit: LOG_PAGE, saving: false,
    // Per table: {col: column index, dir: "asc" | "desc"}; absent = the table's own order.
    sorts: {},
    settledLimit: LOG_PAGE,
    // A calendar-day click (Cal, 2026-10-07): { day, anchor, range, customFirst,
    // customLast } so a second click restores the range it replaced and the
    // calendar keeps its weeks instead of jumping to end on the picked day.
    calendarPick: null,
    // Cal's decision 2026-10-05: the Kalshi bots' combo fills count in every
    // total by default; the header toggle hides them.
    includeBotCombos: true,
  }, loadPrefs(), { view: viewOfHash() });
  const BOT_COMBO_KIND = "Kalshi combo";
  let payload = null;
  // The runner's live game cards ({games, coveredBetIds, ...}) and why the
  // last read failed; the Live tab falls back to a plain list without them.
  let scenarios = null;
  let scenariosError = null;
  let allTickets = [];
  let tickets = [];

  function viewOfHash() {
    const view = location.hash.slice(1);
    return VIEWS.includes(view) ? view : "overview";
  }

  // ---- prefs ----------------------------------------------------------------

  function loadPrefs() {
    try {
      const saved = JSON.parse(localStorage.getItem(PREFS_KEY) || "{}");
      const prefs = {};
      if (RANGES.includes(saved.range)) prefs.range = saved.range;
      if (stats.isDayKey(saved.customFirst)) prefs.customFirst = saved.customFirst;
      if (stats.isDayKey(saved.customLast)) prefs.customLast = saved.customLast;
      if (typeof saved.units === "boolean") prefs.units = saved.units;
      if (typeof saved.includeBotCombos === "boolean") prefs.includeBotCombos = saved.includeBotCombos;
      if (Number.isFinite(saved.unitSize) && saved.unitSize > 0) prefs.unitSize = saved.unitSize;
      const pick = saved.calendarPick;
      if (pick && stats.isDayKey(pick.day) && stats.isDayKey(pick.anchor) && RANGES.includes(pick.range)) prefs.calendarPick = pick;
      return prefs;
    } catch (_error) {
      return {};
    }
  }

  function savePrefs() {
    try {
      localStorage.setItem(PREFS_KEY, JSON.stringify({ range: state.range, customFirst: state.customFirst, customLast: state.customLast, units: state.units, unitSize: state.unitSize, includeBotCombos: state.includeBotCombos, calendarPick: state.calendarPick }));
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

  /** A calendar cell's P&L: thousands as 3.21K (Cal, 2026-10-11), under $1K in whole dollars. */
  function moneyShort(value) {
    const abs = Math.abs(value);
    if (state.units) return sign(value, true) + (abs / state.unitSize).toFixed(1) + "u";
    if (Math.round(abs) === 0) return "0";
    return sign(value, true) + (abs >= 999.5 ? (abs / 1000).toFixed(2) + "K" : String(Math.round(abs)));
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

  /** Tiles {label, value, sub, color?, onClick?}; a tile with onClick is a button-like link. */
  function kpiTiles(id, tiles) {
    fill(id, ...tiles.map((tile) => {
      const link = typeof tile.onClick === "function";
      const node = el("div", {
        className: "panel kpi" + (link ? " link" : ""),
        attrs: link ? { role: "link", tabindex: "0" } : null,
        onClick: link ? tile.onClick : null,
      }, [
        el("span", { className: "lbl", text: tile.label }),
        el("span", { className: "value", text: tile.value, style: { color: tile.color || COLORS.text } }),
        el("span", { className: "sub", text: tile.sub }),
      ]);
      if (link) node.addEventListener("keydown", (event) => { if (event.key === "Enter") tile.onClick(); });
      return node;
    }));
  }

  /**
   * A table; each column {label, right?, num?, cell(row) -> string | Node, className?(row),
   * sort?(row) -> number | string | null, highFirst? (a left-aligned column whose first click sorts high to low)}.
   * opts.sortKey names the table in state.sorts: columns with a `sort` get a
   * clickable header (click sorts, click again reverses) and the rows are
   * ordered by it. opts.rowClass?(row) -> class. opts.onRowClick?(row) makes rows clickable.
   */
  function table(columns, rows, opts) {
    const { sortKey, rowClass, onRowClick } = opts || {};
    const current = sortKey ? state.sorts[sortKey] : null;
    const head = el("tr", null, columns.map((col, index) => headerCell(col, index, sortKey, current)));
    const body = sortedRows(sortKey, columns, rows).map((row) => el("tr", {
      className: rowClass ? rowClass(row) : "", title: onRowClick ? "Click to drill in" : "", onClick: onRowClick ? () => onRowClick(row) : null,
    }, columns.map((col) => {
      const value = col.cell(row);
      const classes = [col.right ? "r" : "", col.num ? "num" : "", col.className ? col.className(row) : ""].filter(Boolean).join(" ");
      const td = el("td", { className: classes });
      if (value instanceof Node) td.appendChild(value); else td.textContent = value;
      return td;
    })));
    return el("table", { className: "tbl" }, [el("thead", null, [head]), el("tbody", null, body)]);
  }

  function headerCell(col, index, sortKey, current) {
    const className = col.right ? "r" : "";
    if (!sortKey || !col.sort) return el("th", { className, text: col.label });
    const active = current && current.col === index;
    const arrow = active ? (current.dir === "asc" ? " ▲" : " ▼") : "";
    const button = el("button", {
      className: "sort" + (active ? " on" : ""), text: col.label + arrow, title: "Sort by " + col.label.toLowerCase(),
      attrs: { type: "button" }, onClick: () => toggleSort(sortKey, index, col),
    });
    const attrs = { "aria-sort": active ? (current.dir === "asc" ? "ascending" : "descending") : "none" };
    return el("th", { className, attrs }, [button]);
  }

  function sortsHighFirst(col) {
    return Boolean(col.right || col.highFirst);
  }

  // A number column sorts high to low first, a text column A to Z; the next click reverses.
  function toggleSort(sortKey, index, col) {
    const current = state.sorts[sortKey];
    const dir = current && current.col === index ? (current.dir === "asc" ? "desc" : "asc") : (sortsHighFirst(col) ? "desc" : "asc");
    state.sorts = Object.assign({}, state.sorts, { [sortKey]: { col: index, dir } });
    render();
  }

  /** `rows` in the table's chosen order, or unchanged when it has none. Paginated tables call this before slicing. */
  function sortedRows(sortKey, columns, rows) {
    const current = sortKey ? state.sorts[sortKey] : null;
    const col = current ? columns[current.col] : null;
    if (!col || !col.sort) return rows;
    return stats.sortRows(rows, col.sort, current.dir);
  }

  /** "sorted by P&L, high to low", or the table's own order when unsorted. */
  function orderCaption(sortKey, columns, defaultOrder) {
    const current = state.sorts[sortKey];
    const col = current ? columns[current.col] : null;
    if (!col || !col.sort) return defaultOrder;
    const words = sortsHighFirst(col) ? { asc: "low to high", desc: "high to low" } : { asc: "A to Z", desc: "Z to A" };
    return "sorted by " + col.label + ", " + words[current.dir];
  }

  function msOf(iso) {
    const ms = iso ? Date.parse(iso) : NaN;
    return Number.isFinite(ms) ? ms : null;
  }

  function segButtons(id, options, isActive, onPick) {
    fill(id, ...options.map((option) => el("button", {
      text: option.label || option, attrs: { type: "button", "aria-pressed": String(isActive(option)) },
      onClick: () => onPick(option),
    })));
  }

  // ---- range ----------------------------------------------------------------

  function rangeDays() {
    return stats.rangeBounds(state.range, stats.pacificDay(Date.now()),
      { first: state.customFirst, last: state.customLast }, stats.firstSettledDay(tickets));
  }

  function rangeCaption(first, last) {
    return (first === last ? dayLabel(first) : dayLabel(first) + " to " + dayLabel(last)) + " · Pacific time";
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
    const open = openTickets();
    const openStake = open.reduce((sum, t) => sum + t.stake, 0);
    const openWithFair = open.filter((t) => t.expected !== null);
    const openEv = openWithFair.reduce((sum, t) => sum + t.expected, 0);
    const missing = noResultCount(first, last);
    setText("ov-caption", rangeCaption(first, last) + (state.units ? " · 1u = $" + state.unitSize : ""));

    const luck = total.withFair ? total.fairPnl - total.expected : null;
    kpiTiles("ov-kpis", [
      { label: "Net P&L", value: money(total.pnl, true), color: toneColor(total.pnl),
        sub: total.bets + " settled · " + total.wins + "-" + total.losses + "-" + total.pushes
          + (total.winRate === null ? "" : " (" + pct(total.winRate) + ")") + (missing ? " · " + missing + " without a result" : "") },
      { label: "ROI", value: pct(total.roi, true), color: toneColor(total.pnl), sub: "on " + money(total.handle) + " handle" },
      { label: "Expected P&L", value: total.withFair ? money(total.expected, true) : "—", color: COLORS.exp,
        sub: total.withFair + " of " + total.bets + " bets have a fair (fill, else close)" },
      { label: "Actual vs expected", value: luck === null ? "—" : money(luck, true), color: luck === null ? COLORS.muted : toneColor(luck),
        sub: total.z === null ? "needs bets with a fair" : "z = " + total.z.toFixed(2) + (Math.abs(total.z) < 1.96 ? ", within noise" : ", outside the 95% band") },
      clvTile(total),
      { label: "Open risk", value: money(openStake), sub: open.length + " open bets" + (openWithFair.length ? " · EV " + money(openEv, true) : "") + " · view →",
        onClick: () => showView("open") },
    ]);

    const series = stats.dailySeries(tickets, first, last);
    renderChart(series, "ov-chart", "ov-stats", "No settled bets in this range.");
    renderCalendar(pickedCalendarDay() ? state.calendarPick.anchor : last);
    renderDaily(series);
    renderVenues(settled);
    renderSettledBets(settled, total, missing);
  }

  /** CLV: expected return at Unabated's closing fair, stake-weighted, and how often the close was beaten. */
  function clvTile(total) {
    return { label: "CLV", value: pct(total.clvRoi, true), color: total.clvRoi === null ? COLORS.muted : toneColor(total.clvRoi),
      sub: total.withClose ? "beat the close " + pct(total.beatRate) + " · " + total.withClose + " of " + total.bets + " bets" : "no closing fairs yet" };
  }

  function clvCell(ticket) {
    return pct(ticket.clv, true);
  }

  function clvClass(ticket) {
    return ticket.clv === null ? "muted" : toneClass(ticket.clv);
  }

  /** Cumulative actual vs expected P&L line with daily bars; drawn into chartId, best/worst/up-days into statsId. */
  function renderChart(series, chartId, statsId, emptyText) {
    const withBets = series.filter((d) => d.bets > 0);
    if (!withBets.length) {
      fill(statsId);
      fill(chartId, emptyNote(emptyText));
      return;
    }
    const best = withBets.reduce((a, b) => (b.pnl > a.pnl ? b : a));
    const worst = withBets.reduce((a, b) => (b.pnl < a.pnl ? b : a));
    const stat = (label, text, color) => el("div", null, [el("span", { className: "lbl", text: label }), el("span", { className: "num", text, style: { color } })]);
    fill(statsId,
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
    const axisDays = [...new Set([series[0], mid, series[series.length - 1]].map((d) => d.day))];
    const xAxis = el("div", { className: "xaxis" }, axisDays.map((day) => el("span", { text: dayLabel(day) })));
    fill(chartId, el("div", { className: "plot" }, [yAxis, svg]), bars, xAxis);
  }

  function renderCalendar(today) {
    const todayIndex = stats.WEEKDAYS.indexOf(stats.weekdayOf(today));
    const start = stats.addDays(today, -todayIndex - (CALENDAR_WEEKS - 1) * 7);
    const end = stats.addDays(start, CALENDAR_WEEKS * 7 - 1);
    const byDay = new Map(stats.dailySeries(tickets, start, today).map((d) => [d.day, d]));
    setText("cal-caption", dayLabel(start) + " to " + dayLabel(end));
    const fullColorDollars = CAL_FULL_COLOR_PNL_UNITS * state.unitSize;
    const cells = stats.WEEKDAYS.map((name) => el("span", { className: "lbl wd", text: name }));
    // Custom ranges clamp to [first settled day, today], so a click outside it
    // would quietly filter to a different day; those days aren't clickable.
    const realToday = stats.pacificDay(Date.now());
    const firstPickable = stats.firstSettledDay(tickets) || realToday;
    for (let day = start; day <= end; day = stats.addDays(day, 1)) {
      const totals = byDay.get(day);
      const isPicked = day === pickedCalendarDay();
      const isPickable = day >= firstPickable && day <= realToday;
      const cell = el("div", {
        className: "cell" + (day === today ? " today" : "") + (isPicked ? " picked" : "") + (isPickable ? " pickable" : ""),
        onClick: isPickable ? () => toggleCalendarDay(day, today) : null,
      });
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
      if (isPickable) cell.title += isPicked ? " · click to clear" : " · click to filter to this day";
      cell.append(dayNumber, value);
      cells.push(cell);
    }
    fill("ov-calendar", ...cells);
  }

  /** The day a calendar click is filtering to, or null once the range has moved off it. */
  function pickedCalendarDay() {
    const pick = state.calendarPick;
    if (!pick || state.range !== "Custom") return null;
    const { first, last } = rangeDays();
    return first === pick.day && last === pick.day ? pick.day : null;
  }

  /** Click a day: filter every Overview panel to it. Click it again: back to the range it replaced. */
  function toggleCalendarDay(day, anchor) {
    const pick = state.calendarPick;
    if (day === pickedCalendarDay()) {
      state.range = pick.range; state.customFirst = pick.customFirst; state.customLast = pick.customLast;
      state.calendarPick = null;
    } else {
      const previous = pickedCalendarDay() ? pick
        : { range: state.range, customFirst: state.customFirst, customLast: state.customLast };
      state.calendarPick = { day, anchor, range: previous.range, customFirst: previous.customFirst, customLast: previous.customLast };
      state.range = "Custom"; state.customFirst = day; state.customLast = day;
    }
    state.logLimit = LOG_PAGE; state.settledLimit = LOG_PAGE; savePrefs(); render();
  }

  function renderDaily(series) {
    const rows = series.filter((d) => d.bets > 0).slice(-DAILY_ROWS).reverse();
    if (!rows.length) { fill("ov-daily", emptyNote("No settled bets in this range.")); return; }
    fill("ov-daily", table([
      { label: "Day", cell: (d) => dayLabel(d.day, true), sort: (d) => d.day },
      { label: "Bets", right: true, num: true, cell: (d) => String(d.bets), sort: (d) => d.bets },
      { label: "W-L-P", right: true, num: true, className: () => "muted", cell: (d) => d.wins + "-" + d.losses + "-" + d.pushes, sort: (d) => d.wins },
      { label: "Handle", right: true, num: true, cell: (d) => money(d.handle), sort: (d) => d.handle },
      { label: "P&L", right: true, num: true, className: (d) => toneClass(d.pnl), cell: (d) => money(d.pnl, true), sort: (d) => d.pnl },
      { label: "ROI", right: true, num: true, className: (d) => toneClass(d.pnl), cell: (d) => pct(d.roi, true), sort: (d) => d.roi },
      { label: "Expected", right: true, num: true, className: () => "exp", cell: (d) => (d.withFair ? money(d.expected, true) : "—"), sort: (d) => (d.withFair ? d.expected : null) },
    ], rows, { sortKey: "daily" }));
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

  function timeLabel(ms) {
    return new Date(ms).toLocaleTimeString("en-US", { hour: "numeric", minute: "2-digit", timeZone: stats.PACIFIC_TZ });
  }

  /**
   * Overview's bottom panel: every bet settled in the header range (Cal,
   * 2026-10-07), newest first, paged like the bet log.
   */
  function renderSettledBets(settled, total, missing) {
    fillSettledSummary("settled-sum", settled, total, missing);
    const columns = settledColumns((t) => dayLabel(t.settledDay, true) + " " + timeLabel(t.closedMs));
    const newestFirst = settled.slice().sort((a, b) => b.closedMs - a.closedMs);
    const shown = sortedRows("settled", columns, newestFirst).slice(0, state.settledLimit);
    setText("settled-caption", settled.length
      ? "Showing " + shown.length + " of " + settled.length + ", " + orderCaption("settled", columns, "newest first")
      : "");
    document.getElementById("settled-more").hidden = settled.length <= shown.length;
    if (!shown.length) { fill("ov-settled", emptyNote("No settled bets in this range.")); return; }
    fill("ov-settled", table(columns, shown, { sortKey: "settled" }));
  }

  /** The Open page's Settled today panel: bets settled on the current Pacific day, whatever the header range. */
  function renderSettledToday() {
    const today = stats.pacificDay(Date.now());
    const settled = stats.inDayRange(tickets, today, today);
    fillSettledSummary("today-sum", settled, stats.summarize(settled), noResultCount(today, today));
    setText("today-caption", dayLabel(today, true) + " · Pacific settle day");
    if (!settled.length) { fill("op-today", emptyNote("Nothing has settled yet today.")); return; }
    const newestFirst = settled.slice().sort((a, b) => b.closedMs - a.closedMs);
    fill("op-today", table(settledColumns((t) => timeLabel(t.closedMs)), newestFirst, { sortKey: "today" }));
  }

  /** P&L, record, expected and handle above a settled-bets table. */
  function fillSettledSummary(id, settled, total, missing) {
    const part = (label, value, className) => el("span", null, [document.createTextNode(label + " "), el("b", { className: className || "", text: value })]);
    fill(id, ...(settled.length ? [
      part("P&L", money(total.pnl, true), toneClass(total.pnl)),
      part("Record", total.wins + "-" + total.losses + "-" + total.pushes),
      total.withFair ? part("Expected", money(total.expected, true), "exp") : null,
      part("Handle", money(total.handle)),
      missing ? part("Without a result", String(missing)) : null,
    ] : []));
  }

  /** Columns of a settled-bets table; settledCell formats the settle time. */
  function settledColumns(settledCell) {
    return [
      { label: "Settled", className: () => "muted", cell: settledCell, sort: (t) => t.closedMs },
      { label: "Venue", cell: (t) => t.venue, sort: (t) => t.venue },
      { label: "League", cell: (t) => el("span", { className: "tag", text: t.league }), sort: (t) => t.league },
      { label: "Event", cell: (t) => t.event || t.kind, sort: (t) => t.event || t.kind },
      { label: "Bet", className: () => "wrap", cell: (t) => t.selection, sort: (t) => t.selection },
      { label: "Price", right: true, num: true, cell: (t) => american(t.displayPrice), sort: (t) => t.displayPrice },
      { label: "Fair", right: true, num: true, className: () => "exp", cell: (t) => american(t.fairAmerican), sort: (t) => t.fairAmerican },
      { label: "Close", right: true, num: true, className: () => "exp", cell: (t) => american(t.closeAmerican), sort: (t) => t.closeAmerican },
      { label: "CLV", right: true, num: true, className: clvClass, cell: clvCell, sort: (t) => t.clv },
      { label: "Stake", right: true, num: true, cell: (t) => money(t.stake), sort: (t) => t.stake },
      { label: "Result", cell: resultTag, sort: (t) => t.status },
      { label: "P&L", right: true, num: true, className: (t) => toneClass(t.pnl), cell: (t) => money(t.pnl, true), sort: (t) => t.pnl },
    ];
  }

  // ---- open -----------------------------------------------------------------

  function openTickets() {
    return tickets.filter((t) => t.status === "open");
  }

  function renderOpenView() {
    const open = openTickets();
    const openStake = open.reduce((sum, t) => sum + t.stake, 0);
    const toWin = open.reduce((sum, t) => sum + (Number.isFinite(t.toWin) ? t.toWin : 0), 0);
    const withFair = open.filter((t) => t.expected !== null);
    const openEv = withFair.reduce((sum, t) => sum + t.expected, 0);
    const { upcoming, noStart } = liveSplit();
    const live = liveTickets();
    const liveStake = live.reduce((sum, t) => sum + t.stake, 0);
    kpiTiles("op-kpis", [
      { label: "Open bets", value: String(open.length),
        sub: live.length + " live · " + upcoming.length + " upcoming" + (noStart.length ? " · " + noStart.length + " no start time" : "") },
      { label: "At risk", value: money(openStake), sub: "to win " + money(toWin) },
      { label: "Open EV", value: withFair.length ? money(openEv, true) : "—", color: withFair.length ? COLORS.exp : COLORS.muted,
        sub: withFair.length + " of " + open.length + " have a saved fair" },
      { label: "Live now", value: money(liveStake), sub: live.length + (live.length === 1 ? " bet" : " bets") + " in progress · view →",
        onClick: () => showView("live") },
    ]);
    renderOpen(upcoming.concat(noStart), { upcoming, noStart });
    renderSettledToday();
  }

  function renderOpen(notLive, split) {
    const stake = notLive.reduce((sum, t) => sum + t.stake, 0);
    const toWin = notLive.reduce((sum, t) => sum + (Number.isFinite(t.toWin) ? t.toWin : 0), 0);
    setText("open-caption", notLive.length ? "Risking " + money(stake) + " to win " + money(toWin) + " · live bets are on the Live tab" : "");
    if (!notLive.length) { fill("op-groups", emptyNote(liveTickets().length ? "Everything open is live. See the Live tab." : "No open bets.")); return; }
    const { upcoming, noStart } = split;
    fill("op-groups",
      openGroup("Upcoming", upcoming, "Starts", "Nothing else open."),
      noStart.length ? openGroup("No start time", noStart, "Starts", "",
        "The venue sends no game time (BetOnline; Kalshi NFL and CFB), or it is a future or combo") : null);
  }

  /** One titled block of the Open bets panel: a count and stake line, then its table. */
  function openGroup(title, group, startLabelText, emptyText, why, sortKey) {
    const stake = group.reduce((sum, t) => sum + t.stake, 0);
    const head = el("div", { className: "open-group-head" }, [
      el("span", { className: "open-group-title", text: title }),
      el("span", { className: "note", text: group.length ? group.length + " · " + money(stake) + " at risk" : "" }),
      why ? el("span", { className: "note muted", text: why }) : null,
    ]);
    if (!group.length) return el("div", { className: "open-group" }, [head, emptyNote(emptyText)]);
    return el("div", { className: "open-group" }, [head, el("div", { className: "scroll" }, [table([
      { label: startLabelText, className: () => "muted", cell: (t) => startLabel(t.eventStart), sort: (t) => msOf(t.eventStart) },
      { label: "Venue", cell: (t) => t.venue, sort: (t) => t.venue },
      { label: "League", cell: (t) => el("span", { className: "tag", text: t.league }), sort: (t) => t.league },
      { label: "Event", cell: (t) => t.event || t.kind, sort: (t) => t.event || t.kind },
      { label: "Bet", className: () => "wrap", cell: (t) => t.selection, sort: (t) => t.selection },
      { label: "Price", right: true, num: true, cell: (t) => american(t.displayPrice), sort: (t) => t.displayPrice },
      { label: "Fair", right: true, num: true, className: () => "exp", cell: (t) => american(t.fairAmerican), sort: (t) => t.fairAmerican },
      { label: "Edge", right: true, num: true, className: (t) => (t.edge === null ? "muted" : toneClass(t.edge)), cell: (t) => pct(t.edge, true), sort: (t) => t.edge },
      { label: "CLV", right: true, num: true, className: clvClass, cell: clvCell, sort: (t) => t.clv },
      { label: "Stake", right: true, num: true, cell: (t) => money(t.stake), sort: (t) => t.stake },
      { label: "To win", right: true, num: true, cell: (t) => (Number.isFinite(t.toWin) ? money(t.toWin) : "—"), sort: (t) => (Number.isFinite(t.toWin) ? t.toWin : null) },
    ], group, { sortKey: sortKey || "open" })])]);
  }

  // ---- live -----------------------------------------------------------------

  /** Open tickets split for the Live and Open tabs off the runner's covered bet ids. */
  function liveSplit() {
    return stats.splitOpenForLive(openTickets(), scenarios ? scenarios.coveredBetIds : [], Date.now());
  }

  function liveTickets() {
    const split = liveSplit();
    return split.onCard.concat(split.tickets, split.offCard);
  }

  function renderLiveView() {
    const split = liveSplit();
    const games = scenarios ? scenarios.games : [];
    const live = split.onCard.concat(split.tickets, split.offCard);
    const stake = live.reduce((sum, t) => sum + t.stake, 0);
    const pricedGroups = games.flatMap((game) => game.groups).filter((group) => group.ev !== null);
    const ev = pricedGroups.reduce((sum, group) => sum + group.ev, 0);
    const best = games.flatMap((game) => game.groups).reduce((sum, group) => sum + group.best, 0);
    const worst = worstCase(games.flatMap((game) => game.groups));
    setText("live-caption", scenarios ? "Updated " + ago(scenarios.generatedAt) : "");
    kpiTiles("lv-kpis", [
      { label: "Live bets", value: String(live.length), sub: games.length + (games.length === 1 ? " game" : " games") + " in progress · " + money(stake) + " at risk" },
      { label: "Best case", value: games.length ? money(best, true) : "—", color: games.length ? COLORS.pos : COLORS.muted, sub: "every game lands your way (straight bets)" },
      { label: "Worst case", value: games.length ? money(worst, true) : "—", color: games.length ? COLORS.neg : COLORS.muted, sub: "nothing lands, every live teaser dies" },
      { label: "EV at kickoff", value: pricedGroups.length ? money(ev, true) : "—", color: pricedGroups.length ? COLORS.exp : COLORS.muted,
        sub: "Unabated's fair before the start (straight bets)" },
    ]);
    const problem = scenariosError
      ? el("div", { className: "banner", text: "Game cards are unavailable: " + scenariosError + ". Live bets are listed below without them." })
      : null;
    fill("lv-cards", problem, ...(games.length ? games.map(gameCard) : [problem ? null : emptyNote(live.length ? "No live bet could be placed on a game card." : "No games in progress.")]));
    const ticketRows = split.tickets;
    document.getElementById("lv-tickets-panel").hidden = !ticketRows.length;
    fill("lv-tickets", ticketRows.length ? openGroup("", ticketRows, "First leg", "", null, "liveTickets") : null);
    document.getElementById("lv-off-panel").hidden = !split.offCard.length;
    fill("lv-off", split.offCard.length ? openGroup("", split.offCard, "Started", "", null, "liveOff") : null);
  }

  /** Each group's worst result added up, a ticket the worst results of several groups kill counted once. */
  function worstCase(groups) {
    let total = 0;
    const killsByTicket = new Map();
    for (const group of groups) {
      total += group.worst;
      const worstBand = group.bands.find((band) => band.pnl === group.worst);
      const killed = new Map(worstBand.legs.filter((leg) => leg.result === "lost" && Number.isFinite(leg.ticketStake)).map((leg) => [leg.parlayId, leg.ticketStake]));
      for (const [ticket, stake] of killed) killsByTicket.set(ticket, { stake, count: (killsByTicket.get(ticket)?.count || 0) + 1 });
    }
    for (const { stake, count } of killsByTicket.values()) total += stake * (count - 1);
    return total;
  }

  /** Legs with the same label, kind and result, one entry each with their distinct tickets' count, stake and to-win. */
  function legsByTicket(legs, keyOf) {
    const byKey = new Map();
    for (const leg of legs) {
      const key = keyOf(leg);
      if (!byKey.has(key)) byKey.set(key, { leg, tickets: new Map() });
      byKey.get(key).tickets.set(leg.parlayId, leg);
    }
    return Array.from(byKey.values()).map(({ leg, tickets }) => {
      const all = Array.from(tickets.values());
      const total = (field) => all.reduce((sum, ticket) => sum + (Number.isFinite(ticket[field]) ? ticket[field] : 0), 0);
      return { leg, count: all.length, stake: total("ticketStake"), toWin: total("ticketToWin") };
    });
  }

  function ticketWord(kind, count) {
    return count + " " + kind + (count === 1 ? "" : "s");
  }

  /** One game in progress: its bets, then one outcome ladder per market and period. */
  function gameCard(game) {
    const placedBets = game.bets.map((bet) => el("span", { className: "bet-chip" + (bet.reason ? " off" : ""), title: bet.reason ? "Not on the ladder: " + bet.reason : "" }, [
      el("span", { text: bet.label }),
      el("span", { className: "muted", text: " · " + bet.venue + (Number.isFinite(bet.stake) ? " · " + money(bet.stake) : "") + (bet.reason ? " · " + bet.reason : "") }),
    ]));
    const legEntries = legsByTicket(game.legs, (leg) => [leg.label, leg.kind, leg.venue, leg.reason].join("|"));
    const legChips = legEntries.map(({ leg, count, stake }) => el("span", { className: "bet-chip leg", title: leg.reason ? "Not on the ladder: " + leg.reason : "" }, [
      el("span", { text: leg.label }),
      el("span", { className: "muted", text: " · in " + ticketWord(leg.kind, count) + " · " + leg.venue + (stake ? " · " + money(stake) : "") + (leg.reason ? " · " + leg.reason : "") }),
    ]));
    const best = game.groups.reduce((sum, group) => sum + group.best, 0);
    const worst = worstCase(game.groups);
    const priced = game.groups.filter((group) => group.ev !== null);
    const ev = priced.reduce((sum, group) => sum + group.ev, 0);
    const summary = el("div", { className: "note" }, [
      el("span", { text: "Range " }), el("span", { className: "num " + toneClass(worst), text: money(worst, true) }),
      el("span", { text: " to " }), el("span", { className: "num " + toneClass(best), text: money(best, true) }),
      priced.length === game.groups.length && priced.length ? el("span", { text: " · EV at kickoff (straight bets) " }) : null,
      priced.length === game.groups.length && priced.length ? el("span", { className: "num exp", text: money(ev, true) }) : null,
      game.oddsAt === null ? el("span", { text: " · no kickoff odds (the runner started after this game did)" }) : null,
    ]);
    return el("section", { className: "panel game-card", attrs: { "aria-label": game.awayTeam + " at " + game.homeTeam } }, [
      el("div", { className: "game-head" }, [
        el("div", null, [
          el("span", { className: "game-teams", text: (game.awayTeam || "Away") + " @ " + (game.homeTeam || "Home") }),
          el("span", { className: "tag", text: String(game.league).toUpperCase() }),
        ]),
        el("span", { className: "note", text: "Started " + startLabel(new Date(game.startMs).toISOString()) }),
      ]),
      el("div", { className: "game-body" }, [
        el("div", { className: "bet-chips" }, placedBets.concat(legChips)),
        el("div", { className: "ladders" }, game.groups.map(outcomeLadder)),
        game.groups.length ? summary : null,
      ]),
    ]);
  }

  /** One market's results: what has to happen, its chance at kickoff, and the P&L. */
  function outcomeLadder(group) {
    const rows = group.bands.map((band) => {
      const bar = el("span", { className: "prob" }, [el("i", { style: { width: band.prob === null ? "0%" : (band.prob * 100).toFixed(1) + "%" } })]);
      const legTags = legsByTicket(band.legs, (leg) => [leg.label, leg.kind, leg.result].join("|")).map(({ leg, count, stake, toWin }) => {
        const tickets = ticketWord(leg.kind, count);
        if (leg.result === "won") return el("span", { className: "leg-tag alive", text: leg.label + " holds: " + tickets + " alive, " + money(toWin) + " to win" });
        if (leg.result === "push") return el("span", { className: "leg-tag push", text: leg.label + " pushes: " + tickets });
        return el("span", { className: "leg-tag", text: leg.label + " loses: " + tickets + ", " + money(-stake, true) });
      });
      return el("div", { className: "band" }, [
        el("span", { className: "band-what" }, [el("span", { text: band.label })].concat(legTags)),
        bar,
        el("span", { className: "num muted band-pct", text: band.prob === null ? "—" : Math.round(band.prob * 100) + "%" }),
        el("span", { className: "num band-pnl " + toneClass(band.pnl), text: money(band.pnl, true) }),
      ]);
    });
    return el("div", { className: "ladder" }, [
      el("div", { className: "ladder-head" }, [el("span", { className: "lbl", text: group.title }), el("span", { className: "lbl", text: "Chance · P&L" })]),
    ].concat(rows));
  }

  // ---- analysis -------------------------------------------------------------

  function periodFilterValue(period) {
    return period === "FG" ? "Full game" : "Partial";
  }

  function filterSpec(key) {
    return FILTERS.find((f) => f.key === key);
  }

  function inOddsRange(ticket) {
    if (state.oddsMin === null && state.oddsMax === null) return true;
    const price = ticket.displayPrice;
    if (!Number.isFinite(price)) return false;
    return (state.oddsMin === null || price >= state.oddsMin) && (state.oddsMax === null || price <= state.oddsMax);
  }

  function passesFilters(ticket) {
    return inOddsRange(ticket) && FILTERS.every((f) => !state.off[f.key].includes(f.value(ticket)));
  }

  function filteredSettled() {
    const { first, last } = rangeDays();
    return stats.inDayRange(tickets, first, last).filter(passesFilters);
  }

  /** Every value of a filter across all tickets, most bets first. */
  function allFilterValues(filter) {
    return valuesByCount(tickets, filter.value).map(([value]) => value);
  }

  /** Filter keys switched on for a value still on the page, in drill order. */
  function activeFilterKeys() {
    const active = FILTERS.filter((f) => tickets.some((t) => state.off[f.key].includes(f.value(t)))).map((f) => f.key);
    if (state.oddsMin !== null || state.oddsMax !== null) active.push("odds");
    const rank = (key) => (state.drillOrder.includes(key) ? state.drillOrder.indexOf(key) : Infinity);
    return active.sort((a, b) => rank(a) - rank(b));
  }

  function noteDrillKey(key, inUse) {
    state.drillOrder = state.drillOrder.filter((k) => k !== key);
    if (inUse) state.drillOrder.push(key);
  }

  function resetPaging() {
    state.logLimit = LOG_PAGE; state.settledLimit = LOG_PAGE;
  }

  function setOff(key, values) {
    state.off[key] = values;
    noteDrillKey(key, values.length > 0);
    resetPaging();
    render();
  }

  function clearFilter(key) {
    if (key === "odds") { state.oddsMin = null; state.oddsMax = null; } else state.off[key] = [];
    noteDrillKey(key, false);
  }

  /** Keep only these values: [[groupKey, groupLabel], ...] from a clicked row or cell. */
  function drillTo(pairs) {
    for (const [key, label] of pairs) {
      const filter = filterSpec(key);
      if (!filter) continue;
      const value = filter.fromGroup ? filter.fromGroup(label) : label;
      state.off[key] = allFilterValues(filter).filter((v) => v !== value);
      noteDrillKey(key, true);
    }
    resetPaging();
    render();
  }

  /** Back up the path: keep the crumbs up to and including `key`, drop the rest. */
  function backTo(key) {
    const keys = activeFilterKeys();
    for (const later of keys.slice(keys.indexOf(key) + 1)) clearFilter(later);
    resetPaging();
    render();
  }

  function resetFilters() {
    for (const key of activeFilterKeys()) clearFilter(key);
    state.drillOrder = [];
    resetPaging();
    render();
  }

  /** "NFL", "Over, Under", "4 of 7"; the odds crumb shows its range. */
  function crumbText(key) {
    if (key === "odds") {
      return "Odds: " + (state.oddsMin === null ? "any" : american(state.oddsMin)) + " to " + (state.oddsMax === null ? "any" : american(state.oddsMax));
    }
    const filter = filterSpec(key);
    const values = allFilterValues(filter);
    const on = values.filter((v) => !state.off[key].includes(v));
    const shown = on.length === 0 ? "none" : on.length <= 2 ? on.join(", ") : on.length + " of " + values.length;
    return filter.label + ": " + shown;
  }

  function renderCrumbs() {
    const keys = activeFilterKeys();
    const nodes = [el("button", { className: "crumb-all", text: "All bets", attrs: { type: "button" }, onClick: resetFilters })];
    for (const key of keys) {
      nodes.push(el("span", { className: "crumb-sep", text: "›" }));
      nodes.push(el("span", { className: "crumb" }, [
        el("button", { className: "crumb-go", text: crumbText(key), title: "Back up to here", attrs: { type: "button" }, onClick: () => backTo(key) }),
        el("button", { className: "crumb-x", text: "×", title: "Drop this filter", attrs: { type: "button", "aria-label": "Drop " + crumbText(key) },
          onClick: () => { clearFilter(key); resetPaging(); render(); } }),
      ]));
    }
    if (!keys.length) nodes.push(el("span", { className: "note", text: "Click any breakdown row or pivot cell to drill in." }));
    fill("f-crumbs", ...nodes);
  }

  /** [value, count] pairs over these tickets, most bets first. */
  function valuesByCount(list, valueOf) {
    const counts = new Map();
    for (const ticket of list) counts.set(valueOf(ticket), (counts.get(valueOf(ticket)) || 0) + 1);
    return [...counts].sort((a, b) => b[1] - a[1]);
  }

  /** A filter's values in its group's order where it has one (odds, timing...), else most bets first. */
  function orderedFilterValues(filter) {
    const values = allFilterValues(filter);
    const group = stats.GROUPS.find((g) => g.key === filter.key);
    if (group && group.order) return values.sort((a, b) => rankIn(group.order, a) - rankIn(group.order, b));
    if (group && group.numeric) return values.sort((a, b) => numberOrInfinity(a) - numberOrInfinity(b));
    return values;
  }

  function rankIn(order, value) {
    return order.includes(value) ? order.indexOf(value) : order.length;
  }

  function numberOrInfinity(label) {
    const value = Number(String(label).replace("−", "-").replace("+", ""));
    return Number.isNaN(value) ? Infinity : value;
  }

  /** One filter row: label, All / None, then a chip per value with its count in the range. */
  function filterRow(filter, inRange) {
    const allValues = orderedFilterValues(filter);
    const inRangeCounts = new Map(valuesByCount(inRange, filter.value));
    const off = state.off[filter.key];
    const quick = el("span", { className: "quick" }, [
      el("button", { text: "All", attrs: { type: "button", "aria-pressed": String(off.length === 0) }, onClick: () => setOff(filter.key, []) }),
      el("button", { text: "None", attrs: { type: "button", "aria-pressed": String(allValues.every((v) => off.includes(v))) }, onClick: () => setOff(filter.key, allValues.slice()) }),
    ]);
    const chipNodes = allValues.map((value) => el("button", {
      className: "fchip", attrs: { type: "button", "aria-pressed": String(!off.includes(value)) },
      onClick: () => setOff(filter.key, off.includes(value) ? off.filter((v) => v !== value) : off.concat(value)),
    }, [document.createTextNode(value), el("span", { className: "n", text: String(inRangeCounts.get(value) || 0) })]));
    return el("div", { className: "filter-row" }, [el("span", { className: "lbl", text: filter.label }), quick, ...chipNodes]);
  }

  /** Parse a typed American price ("-120", "+120", "120"); null when blank, NaN when not a price. */
  function parseAmerican(text) {
    const trimmed = String(text).trim().replace("−", "-");
    if (!trimmed) return null;
    const value = Number(trimmed);
    return Number.isFinite(value) && Math.abs(value) >= 100 ? value : NaN;
  }

  function oddsInput(id, label, current, onSet) {
    const input = el("input", { className: "odds-in num", attrs: { id, type: "text", inputmode: "numeric", placeholder: "any", "aria-label": label } });
    input.value = current === null ? "" : american(current);
    input.addEventListener("change", () => {
      const value = parseAmerican(input.value);
      if (Number.isNaN(value)) { input.setAttribute("aria-invalid", "true"); input.title = "Type an American price like -120 or +120"; return; }
      onSet(value);
      noteDrillKey("odds", state.oddsMin !== null || state.oddsMax !== null);
      resetPaging();
      render();
    });
    return input;
  }

  function oddsRangeRow() {
    return el("div", { className: "filter-row" }, [
      el("span", { className: "lbl", text: "Price" }),
      el("label", { className: "odds-range" }, [document.createTextNode("from"), oddsInput("odds-min", "Lowest price", state.oddsMin, (v) => { state.oddsMin = v; })]),
      el("label", { className: "odds-range" }, [document.createTextNode("to"), oddsInput("odds-max", "Highest price", state.oddsMax, (v) => { state.oddsMax = v; })]),
      el("span", { className: "note", text: "American odds, e.g. −120 to +120" }),
    ]);
  }

  function renderFilters(first, last) {
    const inRange = stats.inDayRange(tickets, first, last);
    const active = activeFilterKeys();
    const shown = FILTERS.filter((f) => f.primary || state.moreFilters || active.includes(f.key));
    const hiddenCount = FILTERS.length - shown.length;
    const toggle = el("button", {
      className: "reset", attrs: { type: "button" },
      text: state.moreFilters ? "Fewer filters" : "More filters (" + hiddenCount + ")",
      onClick: () => { state.moreFilters = !state.moreFilters; render(); },
    });
    fill("f-rows", ...shown.map((filter) => filterRow(filter, inRange)), oddsRangeRow(), (state.moreFilters || hiddenCount) ? toggle : null);
    document.getElementById("f-reset").hidden = !active.length;
    renderCrumbs();
  }

  function renderAnalysis() {
    const { first, last } = rangeDays();
    renderFilters(first, last);
    segButtons("group-tabs", stats.GROUPS, (g) => g.key === state.groupBy, (g) => { state.groupBy = g.key; render(); });

    const settled = filteredSettled();
    const total = stats.summarize(settled);
    setText("an-caption", rangeCaption(first, last) + " · " + total.bets + " settled bets");
    kpiTiles("an-kpis", [
      { label: "Bets", value: String(total.bets), sub: total.wins + "-" + total.losses + "-" + total.pushes + " W-L-P · " + money(total.handle) + " handle" },
      { label: "P&L", value: money(total.pnl, true), color: toneColor(total.pnl), sub: total.withFair ? "expected " + money(total.expected, true) + " on " + total.withFair + " with a fair" : "no fairs" },
      { label: "ROI", value: pct(total.roi, true), color: toneColor(total.roi), sub: "95%: " + pct(total.roi - total.ciHalf, true) + " to " + pct(total.roi + total.ciHalf, true) },
      { label: "Avg edge", value: pct(total.expRoi, true), color: COLORS.exp, sub: "stake-weighted, fill fair else close" },
      clvTile(total),
      { label: "Luck (z)", value: total.z === null ? "—" : total.z.toFixed(2), color: total.z !== null && Math.abs(total.z) >= 1.96 ? COLORS.warn : COLORS.text,
        sub: total.z === null ? "needs bets with a fair" : Math.abs(total.z) < 1.96 ? "within normal variance" : "outside the 95% band" },
    ]);
    renderChart(stats.dailySeries(settled, first, last), "an-chart", "an-stats", "No settled bets match these filters.");
    renderGroups(settled);
    renderPivot(settled);
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
    const group = stats.GROUPS.find((g) => g.key === state.groupBy);
    const rows = stats.groupBy(settled, state.groupBy);
    if (!rows.length) { fill("an-groups", emptyNote("No settled bets match these filters.")); return; }
    // Bucketed groups (odds, edge, timing, weekday, stake) sort in their own order, numbers by value, the rest A to Z.
    const groupSortKey = group.order ? (r) => group.order.indexOf(r.label) : group.numeric ? (r) => numberOrInfinity(r.label) : (r) => r.label;
    fill("an-groups", groupTable(group.label, groupSortKey, rows, (row) => drillTo([[group.key, row.label]])));
  }

  function pivotValue(cell, metric) {
    if (metric === "pnl") return money(cell.pnl, true);
    if (metric === "clv") return pct(cell.clvRoi, true);
    if (metric === "bets") return String(cell.bets);
    return pct(cell.roi, true);
  }

  /** Green for a winning cell, red for a losing one, stronger with ROI; CLV colours by CLV. */
  function pivotShade(cell, metric) {
    const value = metric === "clv" ? cell.clvRoi : cell.roi;
    if (value === null || metric === "bets") return null;
    const strength = Math.min(1, Math.abs(value) / PIVOT_FULL_SHADE_ROI);
    const alpha = (0.08 + 0.32 * strength).toFixed(2);
    return value >= 0 ? "rgba(61, 214, 140, " + alpha + ")" : "rgba(255, 107, 107, " + alpha + ")";
  }

  function pivotSelect(id, current, onPick) {
    const select = el("select", { attrs: { id, "aria-label": id === "pivot-rows" ? "Pivot rows" : "Pivot columns" } },
      stats.GROUPS.map((g) => el("option", { text: g.label, attrs: { value: g.key } })));
    select.value = current;
    select.addEventListener("change", () => onPick(select.value));
    return select;
  }

  function renderPivot(settled) {
    fill("pivot-pick",
      el("label", { className: "pivot-axis" }, [document.createTextNode("Rows"), pivotSelect("pivot-rows", state.pivotRows, (key) => { state.pivotRows = key; render(); })]),
      el("label", { className: "pivot-axis" }, [document.createTextNode("Columns"), pivotSelect("pivot-cols", state.pivotCols, (key) => { state.pivotCols = key; render(); })]));
    segButtons("pivot-metric", PIVOT_METRICS, (m) => m.key === state.pivotMetric, (m) => { state.pivotMetric = m.key; render(); });
    if (state.pivotRows === state.pivotCols) { fill("an-pivot", emptyNote("Pick two different dimensions.")); return; }
    const grid = stats.pivot(settled, state.pivotRows, state.pivotCols);
    if (!grid.rows.length) { fill("an-pivot", emptyNote("No settled bets match these filters.")); return; }
    const metric = state.pivotMetric;
    const cellNode = (cell, pairs) => {
      if (!cell) return el("td", { className: "pv-empty", text: "·" });
      const shade = pivotShade(cell, metric);
      const button = el("button", {
        className: "pv-cell num", attrs: { type: "button" },
        title: cell.bets + " bets · " + money(cell.pnl, true) + " · ROI " + pct(cell.roi, true) + " · CLV " + pct(cell.clvRoi, true) + " · click to drill in",
        style: shade ? { background: shade } : null, onClick: () => drillTo(pairs),
      }, [el("span", { text: pivotValue(cell, metric) }), el("small", { text: String(cell.bets) })]);
      return el("td", null, [button]);
    };
    const head = el("tr", null, [el("th", { text: "" }), ...grid.cols.map((col) => el("th", { className: "r" }, [
      el("button", { className: "pv-head", text: col, attrs: { type: "button" }, title: "Drill into " + col, onClick: () => drillTo([[state.pivotCols, col]]) }),
    ])), el("th", { className: "r", text: "Total" })]);
    const body = grid.rows.map((row) => el("tr", null, [
      el("th", { className: "pv-row" }, [el("button", { className: "pv-head", text: row.label, attrs: { type: "button" }, title: "Drill into " + row.label, onClick: () => drillTo([[state.pivotRows, row.label]]) })]),
      ...row.cells.map((cell, index) => cellNode(cell, [[state.pivotRows, row.label], [state.pivotCols, grid.cols[index]]])),
      cellNode(row.total, [[state.pivotRows, row.label]]),
    ]));
    fill("an-pivot", el("table", { className: "tbl pivot" }, [el("thead", null, [head]), el("tbody", null, body)]));
  }

  /** One breakdown table: a row per group with the Analysis metrics. */
  function groupTable(firstLabel, firstSort, rows, onRowClick) {
    return table([
      { label: firstLabel, cell: (r) => r.label, sort: firstSort },
      { label: "Bets", right: true, num: true, cell: (r) => String(r.bets), sort: (r) => r.bets },
      { label: "W-L-P", right: true, num: true, className: () => "muted", cell: (r) => r.wins + "-" + r.losses + "-" + r.pushes, sort: (r) => r.wins },
      { label: "Handle", right: true, num: true, cell: (r) => money(r.handle), sort: (r) => r.handle },
      { label: "P&L", right: true, num: true, className: (r) => toneClass(r.pnl), cell: (r) => money(r.pnl, true), sort: (r) => r.pnl },
      { label: "ROI", right: true, num: true, className: (r) => toneClass(r.pnl), cell: (r) => pct(r.roi, true), sort: (r) => r.roi },
      // The interval sorts by its low end: high to low puts the most surely winning groups first.
      { label: "ROI, 95% interval", highFirst: true, cell: ciBar, sort: (r) => r.roi - r.ciHalf },
      { label: "Expected ROI", right: true, num: true, className: () => "exp", cell: (r) => pct(r.expRoi, true), sort: (r) => r.expRoi },
      { label: "vs expected", right: true, num: true, className: (r) => (r.withFair ? toneClass(r.fairPnl - r.expected) : "muted"), cell: (r) => (r.withFair ? money(r.fairPnl - r.expected, true) : "—"), sort: (r) => (r.withFair ? r.fairPnl - r.expected : null) },
      { label: "z", right: true, num: true, className: (r) => (r.z !== null && Math.abs(r.z) >= 1.96 ? "warn" : "muted"), cell: (r) => (r.z === null ? "—" : r.z.toFixed(2)), sort: (r) => r.z },
      { label: "CLV", right: true, num: true, className: (r) => (r.clvRoi === null ? "muted" : toneClass(r.clvRoi)), cell: (r) => pct(r.clvRoi, true), sort: (r) => r.clvRoi },
      { label: "Beat close", right: true, num: true, className: () => "muted", cell: (r) => (r.withClose ? pct(r.beatRate) + " of " + r.withClose : "—"), sort: (r) => (r.withClose ? r.beatRate : null) },
    ], rows, { sortKey: "groups", onRowClick, rowClass: () => "drill" });
  }

  function renderCalibration(settled) {
    const bins = stats.calibration(settled);
    if (!bins.length) {
      fill("an-calibration", emptyNote("No settled straight bets with a fair yet."));
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
      el("div", { className: "muted", text: "Fair win probability, at fill or else at close (x), against actual win rate (y); bars are 95% intervals." }));
    fill("an-cal-table", table([
      { label: "Fair prob.", cell: (b) => Math.round(b.low * 100) + " to " + Math.round(b.high * 100) + "%", sort: (b) => b.low },
      { label: "Bets", right: true, num: true, cell: (b) => String(b.bets), sort: (b) => b.bets },
      { label: "Expected win", right: true, num: true, className: () => "exp", cell: (b) => pct(b.expected), sort: (b) => b.expected },
      { label: "Actual win", right: true, num: true, cell: (b) => pct(b.actual), sort: (b) => b.actual },
      { label: "Diff", right: true, num: true, className: (b) => (Math.abs(b.actual - b.expected) > b.ciHalf ? "warn" : "muted"), cell: (b) => pct(b.actual - b.expected, true), sort: (b) => b.actual - b.expected },
    ], bins, { sortKey: "calibration" }));
  }

  function renderLog(settled) {
    const query = state.query.trim().toLowerCase();
    const matching = query
      ? settled.filter((t) => [t.event, t.selection, t.venue, t.league, t.kind].join(" ").toLowerCase().includes(query))
      : settled;
    const columns = [
      { label: "Settled", className: () => "muted", cell: (t) => dayLabel(t.settledDay, true), sort: (t) => t.closedMs },
      { label: "Venue", cell: (t) => t.venue, sort: (t) => t.venue },
      { label: "League", cell: (t) => el("span", { className: "tag", text: t.league }), sort: (t) => t.league },
      { label: "Event", cell: (t) => t.event || t.kind, sort: (t) => t.event || t.kind },
      { label: "Bet", className: () => "wrap", cell: (t) => t.selection, sort: (t) => t.selection },
      { label: "Price", right: true, num: true, cell: (t) => american(t.displayPrice), sort: (t) => t.displayPrice },
      { label: "Fair", right: true, num: true, className: () => "exp", cell: (t) => american(t.fairAmerican), sort: (t) => t.fairAmerican },
      { label: "Edge", right: true, num: true, className: (t) => (t.edge === null ? "muted" : toneClass(t.edge)), cell: (t) => pct(t.edge, true), sort: (t) => t.edge },
      { label: "Stake", right: true, num: true, cell: (t) => money(t.stake), sort: (t) => t.stake },
      { label: "Result", cell: (t) => el("span", { className: "result " + t.status, text: t.status[0].toUpperCase() + t.status.slice(1) }), sort: (t) => t.status },
      { label: "P&L", right: true, num: true, className: (t) => toneClass(t.pnl), cell: (t) => money(t.pnl, true), sort: (t) => t.pnl },
    ];
    const shown = sortedRows("log", columns, matching).slice(0, state.logLimit);
    setText("log-caption", "Showing " + shown.length + " of " + matching.length + " settled bets, " + orderCaption("log", columns, "newest first"));
    document.getElementById("log-more").hidden = matching.length <= shown.length;
    if (!shown.length) { fill("an-log", emptyNote(query ? "No bets match that search." : "No settled bets match these filters.")); return; }
    fill("an-log", table(columns, shown, { sortKey: "log" }));
  }

  // ---- bets -----------------------------------------------------------------

  function resultTag(ticket) {
    return el("span", { className: "result " + ticket.status, text: ticket.status[0].toUpperCase() + ticket.status.slice(1) });
  }

  function betsShown() {
    const query = state.betsQuery.trim().toLowerCase();
    return allTickets.filter((t) => REMOVABLE_VENUES.includes(t.venue)
      && (state.betsVenue === "All" || t.venue === state.betsVenue)
      && (state.betsFilter === "All" || (state.betsFilter === "Removed") === t.excluded)
      && (!query || [t.event, t.selection, t.league, t.kind, t.id].join(" ").toLowerCase().includes(query)));
  }

  // POST the ticket's record ids (every leg of a parlay), then rebuild from
  // the service's answer so the page shows what is stored, not what was asked.
  async function setExcluded(ticket, excluded) {
    if (state.saving) return;
    state.saving = true;
    render();
    try {
      const response = await fetch(EXCLUSIONS_URL, {
        method: "POST", headers: { "Content-Type": "application/json" },
        body: JSON.stringify({ betIds: ticket.betIds, excluded }),
      });
      const body = await response.json().catch(() => ({}));
      if (!response.ok) throw new Error(body.error || "the bets service answered HTTP " + response.status);
      payload = Object.assign({}, payload, { exclusions: body.exclusions });
      rebuildTickets();
      showBanner(null);
    } catch (error) {
      showBanner("Could not " + (excluded ? "remove" : "restore") + " the bet: " + error.message);
    } finally {
      state.saving = false;
      render();
    }
  }

  function renderBets() {
    segButtons("b-venues", ["All"].concat(REMOVABLE_VENUES), (v) => v === state.betsVenue, (v) => { state.betsVenue = v; state.betsLimit = LOG_PAGE; render(); });
    segButtons("b-filter", BET_FILTERS, (f) => f === state.betsFilter, (f) => { state.betsFilter = f; state.betsLimit = LOG_PAGE; render(); });
    const all = allTickets.filter((t) => REMOVABLE_VENUES.includes(t.venue));
    const removed = all.filter((t) => t.excluded);
    setText("b-caption", all.length + " BFA and Wagerzon bets · " + removed.length + " removed"
      + (removed.length ? " (" + money(removed.reduce((sum, t) => sum + (t.pnl || 0), 0), true) + " P&L left out)" : ""));
    const matching = betsShown();
    // The button and the bet lead the row so a phone shows both without scrolling the table sideways.
    const columns = [
      { label: "", cell: (t) => el("button", {
        className: "chip" + (t.excluded ? "" : " danger"), text: t.excluded ? "Restore" : "Remove",
        title: t.excluded ? "Count this bet again" : "Not my bet: leave it out of every number",
        attrs: Object.assign({ type: "button" }, state.saving ? { disabled: "" } : {}),
        onClick: () => setExcluded(t, !t.excluded),
      }) },
      { label: "Bet", className: () => "wrap", cell: (t) => t.selection, sort: (t) => t.selection },
      { label: "Event", cell: (t) => t.event || t.kind, sort: (t) => t.event || t.kind },
      { label: "Price", right: true, num: true, cell: (t) => american(t.displayPrice), sort: (t) => t.displayPrice },
      { label: "Stake", right: true, num: true, cell: (t) => money(t.stake), sort: (t) => t.stake },
      { label: "Result", cell: resultTag, sort: (t) => t.status },
      { label: "P&L", right: true, num: true, className: (t) => toneClass(t.pnl), cell: (t) => (t.pnl === null ? "—" : money(t.pnl, true)), sort: (t) => t.pnl },
      { label: "Placed", className: () => "muted", cell: (t) => (t.placedMs ? dayLabel(stats.pacificDay(t.placedMs), true) : "—"), sort: (t) => t.placedMs },
      { label: "Venue", cell: (t) => t.venue, sort: (t) => t.venue },
      { label: "League", cell: (t) => el("span", { className: "tag", text: t.league }), sort: (t) => t.league },
    ];
    const shown = sortedRows("bets", columns, matching).slice(0, state.betsLimit);
    setText("bets-caption", "Showing " + shown.length + " of " + matching.length + ", " + orderCaption("bets", columns, "newest first"));
    document.getElementById("bets-more").hidden = matching.length <= shown.length;
    if (!shown.length) { fill("b-list", emptyNote(state.betsQuery ? "No bets match that search." : "No bets here.")); return; }
    fill("b-list", table(columns, shown, { sortKey: "bets", rowClass: (t) => (t.excluded ? "removed" : "") }));
  }

  // ---- shell ----------------------------------------------------------------

  function showView(view) {
    state.view = view;
    history.replaceState(null, "", "#" + view);
    window.scrollTo(0, 0);
    render();
  }

  function renderHeader() {
    for (const button of document.querySelectorAll(".nav button")) {
      if (button.dataset.view === state.view) button.setAttribute("aria-current", "page");
      else button.removeAttribute("aria-current");
    }
    for (const view of VIEWS) document.getElementById("view-" + view).hidden = state.view !== view;
    const liveCount = payload ? liveTickets().length : 0;
    const liveBadge = document.getElementById("live-count");
    liveBadge.hidden = !liveCount;
    liveBadge.textContent = String(liveCount);
    const openCount = payload ? openTickets().length - liveCount : 0;
    const badge = document.getElementById("open-count");
    badge.hidden = !openCount;
    badge.textContent = String(openCount);
    segButtons("ranges", RANGES, (r) => r === state.range, pickRange);
    renderCustomRange();
    document.getElementById("show-dollars").setAttribute("aria-pressed", String(!state.units));
    document.getElementById("show-units").setAttribute("aria-pressed", String(state.units));
    document.getElementById("unit-size-label").hidden = !state.units;
    document.getElementById("bot-combos").setAttribute("aria-pressed", String(state.includeBotCombos));
    const unitInput = document.getElementById("unit-size");
    if (document.activeElement !== unitInput) unitInput.value = String(state.unitSize);
  }

  /** Picking Custom starts from the range on screen, so the dates are never blank. */
  function pickRange(range) {
    if (range === "Custom" && state.range !== "Custom") {
      const { first, last } = rangeDays();
      state.customFirst = first; state.customLast = last;
    }
    state.calendarPick = null;
    state.range = range; state.logLimit = LOG_PAGE; state.settledLimit = LOG_PAGE; savePrefs(); render();
  }

  function renderCustomRange() {
    document.getElementById("custom-range").hidden = state.range !== "Custom";
    if (state.range !== "Custom") return;
    const { first, last } = rangeDays();
    const today = stats.pacificDay(Date.now());
    for (const [id, value] of [["custom-first", first], ["custom-last", last]]) {
      const input = document.getElementById(id);
      input.max = today;
      if (document.activeElement !== input) input.value = value;
    }
  }

  function renderSync() {
    const dot = document.getElementById("sync-dot");
    const sources = Object.entries((payload && payload.sources) || {});
    if (!sources.length) { setText("sync-text", payload ? "No venues reporting" : "Loading…"); return; }
    const failing = sources.filter(([, s]) => !s.ok);
    dot.classList.toggle("bad", failing.length > 0);
    const okCount = sources.length - failing.length;
    const failingNames = failing.map(([name]) => stats.venueName(name)).join(", ");
    setText("sync-text", okCount + " of " + sources.length + " venues ok · "
      + (failing.length ? failingNames + " down · " : "") + ago(payload.generatedAt));
    document.getElementById("sync").title = failing.length
      ? failing.map(([name, s]) => stats.venueName(name) + ": " + (s.error || "no completed poll yet")).join("\n")
      : "Every venue's last poll succeeded";
  }

  // Bets no rule could give a market: hidden at zero, else the count and venues.
  function renderUntyped() {
    const element = document.getElementById("untyped");
    const summary = stats.untypedSummary(allTickets);
    element.hidden = summary.count === 0;
    if (!summary.count) return;
    setText("untyped-text", summary.count + " untyped · " + summary.venues.map((v) => v.venue).join(", "));
    element.title = "Bets with no market type (counted as Other); the bet type rules in trackerstats.js need a case for them:\n"
      + summary.examples.join("\n") + (summary.count > summary.examples.length ? "\n…" : "");
  }

  function render() {
    renderHeader();
    renderSync();
    renderUntyped();
    if (!payload) return;
    if (state.view === "overview") renderOverview();
    else if (state.view === "live") renderLiveView();
    else if (state.view === "open") renderOpenView();
    else if (state.view === "analysis") renderAnalysis();
    else renderBets();
  }

  function showBanner(text) {
    const banner = document.getElementById("banner");
    banner.hidden = !text;
    banner.textContent = text || "";
  }

  // Every number leaves out removed bets; the bot-combo toggle may leave out those too.
  function applyBotComboFilter() {
    tickets = allTickets.filter((t) => !t.excluded && (state.includeBotCombos || t.kind !== BOT_COMBO_KIND));
  }

  function rebuildTickets() {
    allTickets = stats.buildTickets(payload.bets, payload.fillFairs, payload.closingFairs, payload.exclusions);
    applyBotComboFilter();
  }

  // The live game cards; a failure keeps the last cards read and says why.
  async function refreshScenarios() {
    try {
      const response = await fetch(SCENARIOS_URL, { cache: "no-store" });
      const body = await response.json().catch(() => null);
      if (!response.ok) throw new Error((body && body.error) || "the bets service answered HTTP " + response.status);
      if (!body || !Array.isArray(body.games) || !Array.isArray(body.coveredBetIds)) throw new Error("scenarios.json has no games list");
      scenarios = body;
      scenariosError = null;
    } catch (error) {
      scenariosError = error.message;
    }
  }

  async function refresh() {
    const scenariosRead = refreshScenarios();
    try {
      const response = await fetch(BETS_URL, { cache: "no-store" });
      if (!response.ok) throw new Error("the bets service answered HTTP " + response.status);
      payload = await response.json();
      rebuildTickets();
      showBanner(null);
    } catch (error) {
      showBanner("Could not load bets: " + error.message + (payload ? ". Showing the last load." : "."));
    }
    await scenariosRead;
    render();
  }

  function wire() {
    for (const button of document.querySelectorAll(".nav button")) {
      button.addEventListener("click", () => showView(button.dataset.view));
    }
    document.getElementById("settled-more").addEventListener("click", () => { state.settledLimit += LOG_PAGE; render(); });
    document.getElementById("f-reset").addEventListener("click", resetFilters);
    document.getElementById("show-dollars").addEventListener("click", () => { state.units = false; savePrefs(); render(); });
    document.getElementById("show-units").addEventListener("click", () => { state.units = true; savePrefs(); render(); });
    document.getElementById("bot-combos").addEventListener("click", () => {
      state.includeBotCombos = !state.includeBotCombos;
      state.logLimit = LOG_PAGE; state.settledLimit = LOG_PAGE;
      savePrefs(); applyBotComboFilter(); render();
    });
    document.getElementById("unit-size").addEventListener("change", (event) => {
      const size = Number(event.target.value);
      if (Number.isFinite(size) && size > 0) { state.unitSize = size; savePrefs(); }
      render();
    });
    for (const [id, key] of [["custom-first", "customFirst"], ["custom-last", "customLast"]]) {
      document.getElementById(id).addEventListener("change", (event) => {
        if (!stats.isDayKey(event.target.value)) return;
        state.calendarPick = null;
        state[key] = event.target.value; state.logLimit = LOG_PAGE; state.settledLimit = LOG_PAGE; savePrefs(); render();
      });
    }
    document.getElementById("log-search").addEventListener("input", (event) => {
      state.query = event.target.value; state.logLimit = LOG_PAGE; state.settledLimit = LOG_PAGE; render();
    });
    document.getElementById("log-more").addEventListener("click", () => { state.logLimit += LOG_PAGE; render(); });
    document.getElementById("bets-search").addEventListener("input", (event) => {
      state.betsQuery = event.target.value; state.betsLimit = LOG_PAGE; render();
    });
    document.getElementById("bets-more").addEventListener("click", () => { state.betsLimit += LOG_PAGE; render(); });
    window.addEventListener("hashchange", () => { state.view = viewOfHash(); render(); });
    document.addEventListener("visibilitychange", () => { if (!document.hidden) refresh(); });
    setInterval(() => { if (!document.hidden) refresh(); }, POLL_MS);
  }

  wire();
  render();
  refresh();
})();
