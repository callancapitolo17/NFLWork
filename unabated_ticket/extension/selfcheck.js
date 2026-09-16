// Unabated Ticket — bundle self-check (#123).
//
// page.js reads Unabated's React fiber props and AG Grid rows. When Unabated
// ships a new bundle those shapes move and the result is silent: an empty
// edges list, a ticket that captures nothing, "Not watching the line" — none
// of it distinguishable from a quiet board. This module turns one raw sample
// of the page (gathered by page.js, which is the only side that may touch the
// DOM) into a per-probe pass/fail verdict the panel can name.
//
// Pure functions: no DOM, no fetch, no chrome.* — loaded as a plain <script>
// in panel.html and as a MAIN-world content script beside page.js (both
// expose globalThis.UnabatedSelfCheck), and via require() in
// tests/page_rows.test.js.
//
// The sample page.js hands `run()` (raw reads, no decisions):
//   {
//     at, url, build,                 // when, where, which page.js build
//     path, sinceLoadMs,              // location.pathname, age of this page.js copy
//     shells,                         // odds-cell elements rendered (int)
//     rowElements, gridRoots,         // .ag-row and .ag-root nodes (int)
//     rows,                           // grid rows carrying `sides` (int)
//     lines: [line, ...],             // book entries harvested off those rows
//     cell: null | {                  // one rendered odds cell, as capture reads it
//       hasFiber, hasLineProps, hasGridProps,
//       marketLine, sideIndex, rowData, bookId,
//     },
//     altsExpandable,                 // bool|null: setExpanded on a top-level row node
//     bookSelection: { bookIds, error, entryCount },  // userSettings.gameOdds
//     locate: null | { rowElement, bookCell },  // the row/cell selectors locate uses
//   }
//
// Shape, not value: a probe passes when the FIELD is still there on some line,
// even if today's value is null. A book that stops publishing a fair is not a
// bundle change, and treating it as one would cry wolf on every quiet board.
// Value counts go in each probe's detail line instead.
//
// Verdicts are "pass" / "fail" / "unknown" (nothing on the page to judge yet).
// `status` is "waiting" when no odds cell is rendered at all (a loading tab, or
// a tools.unabated.com page that is not the odds screen — never an alarm),
// "changed" when any probe failed, else "ok".

(function (root) {
  "use strict";

  const PASS = "pass";
  const FAIL = "fail";
  const UNKNOWN = "unknown";

  function isFiniteNumber(value) {
    return typeof value === "number" && Number.isFinite(value);
  }

  function linesOf(sample) {
    return Array.isArray(sample.lines) ? sample.lines.filter((line) => line && typeof line === "object") : [];
  }

  function countLines(sample, has) {
    const lines = linesOf(sample);
    return { total: lines.length, hits: lines.filter(has).length };
  }

  // A line-shape probe: `carries` is the field being read, `valued` the
  // values that are actually usable. Unknown (not failed) with no lines to
  // look at — an empty board proves nothing about the bundle.
  function lineProbe(carries, valued) {
    return (sample) => {
      const field = countLines(sample, carries);
      if (field.total === 0) return { state: UNKNOWN, detail: "no book lines on the board to read" };
      const usable = countLines(sample, (line) => carries(line) && valued(line));
      return {
        state: field.hits > 0 ? PASS : FAIL,
        detail: `${field.hits} of ${field.total} sampled lines carry the field (${usable.hits} with a usable value)`,
      };
    };
  }

  function cellOf(sample) {
    return sample.cell && typeof sample.cell === "object" ? sample.cell : null;
  }

  // "/nfl/odds", "/cfb/odds/..." — the screen this extension reads. A tab
  // sitting on one with nothing rendered long after load is a break, not a load.
  function isOddsScreen(path) {
    return typeof path === "string" && /\/odds(\/|$)/.test(path);
  }

  const GRID_GRACE_MS = 60000;

  // scope: which panel views a failure makes wrong — "ticket" (a capture),
  // "edges" (the Edges list and its book filter) or "both".
  const PROBES = [
    {
      id: "cells",
      label: "odds cells on the page",
      scope: "both",
      check: (sample) => {
        if (sample.shells > 0) return { state: PASS, detail: `${sample.shells} odds cells, ${sample.gridRoots} grid roots rendered` };
        // Rows rendered with no odds cell in any of them: a renamed cell class,
        // which is exactly the silent break this file exists to catch.
        if (sample.rowElements > 0) return { state: FAIL, detail: `${sample.rowElements} grid rows rendered but no odds cells matched` };
        // A grid with no rows is an empty board or the user's own filters.
        if (sample.gridRoots > 0) return { state: UNKNOWN, detail: "the odds grid is up with no rows" };
        if (isOddsScreen(sample.path) && (sample.sinceLoadMs ?? 0) >= GRID_GRACE_MS) {
          return { state: FAIL, detail: `nothing rendered on an odds screen ${Math.round((sample.sinceLoadMs ?? 0) / 1000)}s after load` };
        }
        return { state: UNKNOWN, detail: "no odds grid on this page yet" };
      },
    },
    {
      id: "fiber",
      label: "the cell's React props",
      scope: "both",
      check: (sample) => {
        const cell = cellOf(sample);
        if (!cell) return { state: FAIL, detail: "no odds cell exposed React props" };
        const missing = [];
        if (!cell.hasFiber) missing.push("no __reactFiber$ on the cell");
        if (!cell.hasLineProps) missing.push("no marketLine + sideIndex props above it");
        if (!cell.marketLine) missing.push("marketLine is empty");
        if (!Number.isInteger(cell.sideIndex)) missing.push(`sideIndex is ${JSON.stringify(cell.sideIndex)}`);
        return missing.length
          ? { state: FAIL, detail: missing.join("; ") }
          : { state: PASS, detail: `marketLine + sideIndex ${cell.sideIndex} read off a cell` };
      },
    },
    {
      id: "grid",
      label: "the grid row behind a cell",
      scope: "both",
      check: (sample) => {
        const cell = cellOf(sample);
        const rowData = cell && cell.rowData;
        if (!cell || !cell.hasGridProps) return { state: FAIL, detail: "no AG Grid cell props (api + node.data) above the cell" };
        if (!rowData || !rowData.sides || typeof rowData.sides !== "object") {
          return { state: FAIL, detail: `the row carries no sides object (keys: ${Object.keys(rowData || {}).slice(0, 10).join(",") || "none"})` };
        }
        const sideKeys = Object.keys(rowData.sides);
        const keyed = sideKeys.filter((key) => /^si\d+:/.test(key)).length;
        if (!keyed) return { state: FAIL, detail: `sides keys are not "si<n>:..." (${sideKeys.slice(0, 6).join(",") || "none"})` };
        return { state: PASS, detail: `${sample.rows} rows with sides; this one keyed ${sideKeys.slice(0, 4).join(",")}` };
      },
    },
    {
      id: "edge",
      label: "Unabated's edge on a line",
      scope: "ticket",
      check: lineProbe(
        (line) => "edge" in line || "ge" in line,
        (line) => isFiniteNumber(line.edge && line.edge.edge) || isFiniteNumber(line.ge),
      ),
    },
    {
      id: "fair",
      label: "Unabated's fair price (bacr)",
      scope: "ticket",
      check: lineProbe((line) => "bacr" in line, (line) => isFiniteNumber(line.bacr)),
    },
    {
      id: "sourcePrice",
      label: "the exchange source price",
      scope: "ticket",
      check: lineProbe(
        (line) => "sourcePrice" in line && "sourceFormat" in line,
        (line) => (line.sourceFormat === 2 || line.sourceFormat === 4) && isFiniteNumber(line.sourcePrice) && line.sourcePrice > 0,
      ),
    },
    {
      id: "books",
      label: "your Unabated book selection",
      scope: "edges",
      check: (sample) => {
        const selection = sample.bookSelection || {};
        const ids = Array.isArray(selection.bookIds) ? selection.bookIds : null;
        if (ids && ids.length) return { state: PASS, detail: `${ids.length} books enabled under userSettings.gameOdds` };
        // The shape is what is probed: entries found and none of them enabled
        // is a user who has unticked every book, not a moved field.
        if (selection.entryCount > 0) {
          return { state: PASS, detail: `${selection.entryCount} book entries, none of them enabled (your selection, not a change)` };
        }
        return { state: FAIL, detail: selection.error || "userSettings.gameOdds yielded no book entries" };
      },
    },
    {
      id: "alts",
      label: "the Alts expander",
      scope: "both",
      check: (sample) => {
        if (sample.altsExpandable == null) return { state: UNKNOWN, detail: "no top-level grid row to try" };
        return sample.altsExpandable
          ? { state: PASS, detail: "setExpanded is on the row node" }
          : { state: FAIL, detail: "no setExpanded on a top-level row node; alt lines cannot be opened or located" };
      },
    },
    {
      id: "cellLookup",
      label: "finding a cell by number, side, book and game",
      scope: "both",
      check: (sample) => {
        const cell = cellOf(sample);
        if (!cell || !cell.marketLine || !cell.rowData) return { state: UNKNOWN, detail: "no cell to match against" };
        const row = cell.rowData;
        const missing = [];
        // A moneyline has no number; every other market must carry one.
        if (row.betTypeId !== 1 && !isFiniteNumber(cell.marketLine.points)) missing.push("the cell carries no points");
        if (!Number.isInteger(cell.sideIndex)) missing.push("no side index on the cell");
        if (cell.bookId == null) missing.push("the book could not be identified");
        if (row.eventId == null) missing.push("no eventId on the row");
        if (row.betTypeId == null) missing.push("no betTypeId on the row");
        const locate = sample.locate || {};
        if (locate.rowElement === false) missing.push('no .ag-row[row-id] for the row');
        if (locate.rowElement === true && locate.bookCell === false) missing.push("the row has no cell in the book's column");
        return missing.length
          ? { state: FAIL, detail: missing.join("; ") }
          : { state: PASS, detail: `book ${cell.bookId}, side ${cell.sideIndex}, event ${row.eventId} bt${row.betTypeId}` };
      },
    },
  ];

  // Every other probe reads through a cell or the grid behind one, so when
  // there are no cells they can say nothing — unknown, not failed. Whether
  // THAT is an alarm is the cells probe's own call: unknown while a tab is
  // still loading (status "waiting", silent), failed when the grid is up and
  // the cells are not (status "changed", loud).
  const GATE_PROBE_ID = "cells";

  function run(sample) {
    const input = sample && typeof sample === "object" ? sample : {};
    const gate = safeCheck(PROBES.find((probe) => probe.id === GATE_PROBE_ID), input);
    const probes = PROBES.map((probe) => {
      const verdict = probe.id === GATE_PROBE_ID
        ? gate
        : gate.state !== PASS
          ? { state: UNKNOWN, detail: "waiting for the odds grid" }
          : safeCheck(probe, input);
      return { id: probe.id, label: probe.label, scope: probe.scope, state: verdict.state, detail: verdict.detail };
    });
    const failed = probes.filter((probe) => probe.state === FAIL);
    return {
      at: input.at ?? null,
      url: input.url ?? null,
      build: input.build ?? null,
      status: failed.length ? "changed" : gate.state === PASS ? "ok" : "waiting",
      probes,
      failed: failed.map((probe) => probe.id),
    };
  }

  // A probe that throws IS a change (the shape it reached into moved), but it
  // must not take the rest of the self-check down with it.
  function safeCheck(probe, sample) {
    try {
      return probe.check(sample);
    } catch (error) {
      return { state: FAIL, detail: `probe threw: ${error && error.message}` };
    }
  }

  // page.js could not even gather the sample: report it as a change rather
  // than staying silent, which is the failure mode this whole file exists for.
  function errorReport({ at, url, build, message }) {
    return {
      at: at ?? null,
      url: url ?? null,
      build: build ?? null,
      status: "changed",
      probes: [{ id: "sample", label: "reading the page", scope: "both", state: FAIL, detail: String(message) }],
      failed: ["sample"],
    };
  }

  function failedProbes(report) {
    if (!report || report.status !== "changed" || !Array.isArray(report.probes)) return [];
    return report.probes.filter((probe) => probe.state === FAIL);
  }

  // The panel's one loud line, or null when there is nothing to shout about.
  function summaryOf(report) {
    const failed = failedProbes(report);
    if (!failed.length) return null;
    return `Unabated changed: ${failed.map((probe) => probe.label).join(", ")}`;
  }

  // One line per failed probe, for the detail under the banner.
  function detailOf(report) {
    return failedProbes(report).map((probe) => `${probe.label} — ${probe.detail}`).join("\n");
  }

  // What a view should say instead of rendering empty, or null when this
  // view's own reads are unaffected by what failed.
  const VIEW_NOTES = {
    ticket: "Unabated's page changed, so clicking a price may capture nothing",
    edges: "Unabated's page changed, so the book selection read from the tab may be wrong",
  };

  function noteFor(report, where) {
    const relevant = failedProbes(report).filter((probe) => probe.scope === where || probe.scope === "both");
    if (!relevant.length || !VIEW_NOTES[where]) return null;
    return `${VIEW_NOTES[where]}: ${relevant.map((probe) => probe.label).join(", ")}.`;
  }

  const api = { PROBES, PASS, FAIL, UNKNOWN, run, errorReport, summaryOf, detailOf, noteFor };

  if (typeof module !== "undefined" && module.exports) {
    module.exports = api;
  } else {
    root.UnabatedSelfCheck = api;
  }
})(typeof globalThis !== "undefined" ? globalThis : this);
