# MLB Dashboard — Header Pill Row Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Replace the standalone Wagerzon multi-account header bar with a click-to-select pill row living inside the existing dashboard header. Fix the "stale 0s ago" rendering bug as part of the pill-renderer rewrite.

**Architecture:** Layout-only change. Single file modified (`mlb_dashboard.R`) plus a one-line README touch-up. The standalone `<div id="wz-account-bar">` above `.container` is removed; pill row + ↻ button are placed inside `.header` as a second row. The `<select>` dropdown is replaced by clickable pills using new CSS classes (`.wz-pill`, `.wz-pill.selected`, etc.). No backend, schema, or API contract changes.

**Tech Stack:** R Shiny (`mlb_dashboard.R` generates `report.html`), vanilla JavaScript (embedded `tags$script`), Flask server (`mlb_dashboard_server.py`, untouched), DuckDB (`mlb_dashboard.duckdb`, untouched).

**Worktree:** `.worktrees/mlb-dash-header-merge` on branch `feature/mlb-dashboard-header-merge`. Spec at `docs/superpowers/specs/2026-05-02-mlb-dashboard-header-merge-design.md` (committed `1b55f96`).

---

## Pre-flight

These steps confirm you're in the right place. Run them once at the start.

- [ ] **Step 0.1: Confirm worktree + branch**

```bash
cd /Users/callancapitolo/NFLWork/.worktrees/mlb-dash-header-merge
pwd
git branch --show-current
git log --oneline -1
```

Expected: working dir is the worktree path, branch is `feature/mlb-dashboard-header-merge`, latest commit is `1b55f96 docs(mlb-dashboard): design — merge WZ account bar into header pill row`.

- [ ] **Step 0.2: Confirm dashboard is running and reachable**

```bash
curl -sI http://localhost:8083/ | head -1
curl -s http://localhost:8083/api/wagerzon/balances | head -c 200
```

Expected: `HTTP/1.1 200 OK` and a JSON `{"balances":[...]}` payload (the current `wz_error` state is fine — UI will degrade gracefully).

If the dashboard isn't running, ask the user to start it (`bash "Answer Keys/MLB Dashboard/run.sh"` from the **main** repo, not the worktree — the running server is shared).

---

## Task 1: Add new CSS classes

**Files:**
- Modify: `Answer Keys/MLB Dashboard/mlb_dashboard.R` (the embedded `<style>` block, current line ~1759 region — search for the existing `.pill {` rule that styles the books strip; the new classes go in their own block above the existing parlay-tab section so they're easy to find).

These classes are pure additions. The page still renders identically until Task 2 starts using them.

- [ ] **Step 1.1: Locate the existing `<style>` block**

```bash
grep -n "Parlay tab — books strip" "Answer Keys/MLB Dashboard/mlb_dashboard.R"
```

Expected: one line, currently `1759:        /* Parlay tab — books strip (M / DK / FD / PX / NV / Cons pill row) */`.

- [ ] **Step 1.2: Insert the new CSS block immediately *before* that comment**

Use the Edit tool. `old_string` is the exact existing comment line plus enough surrounding context to be unique. `new_string` is the new CSS block followed by the original comment line.

```r
        /* === Wagerzon multi-account header pill row (Phase 7 layout merge) ===
           Lives inside .header as a second row. Replaces the old
           full-bleed wz-account-bar that sat above .container. */
        .header-row-top {
          display: flex;
          justify-content: space-between;
          align-items: center;
          padding-bottom: 10px;
        }
        .header-row-accounts {
          display: flex;
          align-items: center;
          gap: 8px;
          padding: 10px 0 4px 0;
        }
        .header-label {
          font-size: 11px;
          color: #8b949e;
          text-transform: uppercase;
          letter-spacing: 0.5px;
          margin-right: 2px;
        }
        .wz-pill {
          padding: 4px 10px;
          border-radius: 14px;
          background: #21262d;
          border: 1px solid #30363d;
          color: #c9d1d9;
          font-size: 13px;
          cursor: pointer;
          user-select: none;
          transition: border-color 0.12s, background 0.12s;
        }
        .wz-pill:hover { border-color: #58a6ff; }
        .wz-pill.selected {
          background: #1f6feb;
          border-color: #1f6feb;
          color: #ffffff;
          font-weight: 600;
        }
        .wz-pill.stale {
          background: #3a1d1d;
          border-color: #4a2a2a;
          color: #ffa198;
        }
        .wz-pill.empty {
          background: transparent;
          border-style: dashed;
          color: #6e7681;
          cursor: default;
        }
        .wz-pill.empty:hover { border-color: #30363d; }
        .wz-icon-btn {
          background: transparent;
          border: 1px solid #30363d;
          color: #8b949e;
          width: 28px;
          height: 28px;
          border-radius: 6px;
          cursor: pointer;
          display: inline-flex;
          align-items: center;
          justify-content: center;
          font-size: 14px;
          padding: 0;
        }
        .wz-icon-btn:hover { color: #c9d1d9; border-color: #58a6ff; }

        /* Parlay tab — books strip (M / DK / FD / PX / NV / Cons pill row) */
```

The trailing `/* Parlay tab — books strip ... */` line is the existing comment we matched on; preserving it is what keeps the Edit tool's `old_string` matched.

- [ ] **Step 1.3: Regenerate the dashboard HTML**

```bash
cd /Users/callancapitolo/NFLWork
Rscript "Answer Keys/MLB Dashboard/mlb_dashboard.R"
```

Expected: completes without R errors. (`Warning` lines about deprecated reactable args are pre-existing and fine.)

- [ ] **Step 1.4: Confirm new classes are present in the rendered HTML**

```bash
grep -c "\.wz-pill" "Answer Keys/MLB Dashboard/report.html"
grep -c "\.header-row-accounts" "Answer Keys/MLB Dashboard/report.html"
```

Expected: both > 0. The classes are now in the served CSS even though no DOM uses them yet.

- [ ] **Step 1.5: Sanity-check the live page**

Refresh http://localhost:8083 in the browser. The page should look **identical** to before — the new CSS is dormant until DOM/JS use it.

If anything visually changed: revert the Edit and re-apply, checking you didn't accidentally consume an unrelated CSS rule.

**Do NOT commit yet.** Task 1 leaves the page in a half-state (CSS without DOM); commits happen at the end after the full feature works.

---

## Task 2: Restructure the header DOM

**Files:**
- Modify: `Answer Keys/MLB Dashboard/mlb_dashboard.R` lines ~1992–2025 (the `tags$body(...)` opener through the existing `.header` block).

Two edits:

1. **Remove** the entire `tags$div(id = "wz-account-bar", ...)` block above `.container`.
2. **Restructure** the `.header` div inside `.container` from a single flex row into two rows: `.header-row-top` (existing title + Refresh) and `.header-row-accounts` (new pill row).

- [ ] **Step 2.1: Remove the standalone account bar**

Use the Edit tool. Match the entire `wz-account-bar` block plus the `tags$div(class = "container",` opener that follows it, so the replacement keeps the container intact.

`old_string` is the block from `      # Wagerzon multi-account header bar (Phase 6).` through and including the `tags$select(id = "wz-account-select", style = "padding:4px 8px; font-size:14px;")` and the closing `)` and `,` of that outer div, ending just before `      tags$div(class = "container",`.

`new_string` is empty — drop the entire block. The `tags$div(class = "container", ...)` line that followed becomes the first child of `tags$body(`.

After the edit, the structure should be:

```r
    tags$body(
      tags$div(class = "container",
        # Header
        tags$div(class = "header",
          ...
        ),
        ...
      )
    )
```

- [ ] **Step 2.2: Restructure `.header` into two rows**

Locate the existing `.header` block:

```r
        # Header
        tags$div(class = "header",
          tags$div(
            tags$h1("MLB Answer Key Dashboard"),
            tags$div(class = "subtitle", paste("Updated", timestamp))
          ),
          tags$button(class = "refresh-btn", onclick = "refreshData()", "Refresh")
        ),
```

Replace it with:

```r
        # Header — two rows.
        # Row 1: title + subtitle (left) | "Refresh" data button (right).
        # Row 2: "Placing on" caption + clickable Wagerzon account pills +
        #        balance-refresh icon. Replaces the old full-bleed
        #        wz-account-bar that lived above .container.
        tags$div(class = "header",
          tags$div(class = "header-row-top",
            tags$div(
              tags$h1("MLB Answer Key Dashboard"),
              tags$div(class = "subtitle", paste("Updated", timestamp))
            ),
            tags$button(class = "refresh-btn", onclick = "refreshData()", "Refresh")
          ),
          tags$div(class = "header-row-accounts", id = "wz-account-row",
            tags$span(class = "header-label", "Placing on"),
            tags$div(id = "wz-account-pills", style = "display:flex; gap:6px; align-items:center;"),
            tags$button(
              id = "wz-refresh-btn", type = "button", class = "wz-icon-btn",
              title = "Refresh balances",
              HTML("&#x21bb;")
            )
          )
        ),
```

Notes:
- `#wz-account-pills` keeps its old ID so the existing JS controller's `BAR_ID = 'wz-account-pills'` keeps matching.
- `#wz-refresh-btn` keeps its old ID.
- The `<select id="wz-account-select">` is gone — JS for it is removed in Task 3.
- The pill container's inline `style` (`display:flex; gap:6px; align-items:center;`) keeps pills laid out correctly even if a future class change forgets it.

- [ ] **Step 2.3: Regenerate and inspect HTML**

```bash
Rscript "Answer Keys/MLB Dashboard/mlb_dashboard.R"
grep -c "wz-account-bar" "Answer Keys/MLB Dashboard/report.html"
grep -c "wz-account-row" "Answer Keys/MLB Dashboard/report.html"
grep -c "wz-account-select" "Answer Keys/MLB Dashboard/report.html"
```

Expected:
- `wz-account-bar` → `0` (removed)
- `wz-account-row` → `1` (new)
- `wz-account-select` → `0` (removed)

- [ ] **Step 2.4: Reload http://localhost:8083 — expect a transient broken state**

The page renders with an empty pill row (no pills, no caption interaction). The browser console will show errors from the existing JS controller calling `document.getElementById('wz-account-select')` — that's expected; Task 3 fixes it.

Do **not** stop here. Move directly to Task 3.

---

## Task 3: Rewrite the JS controller

**Files:**
- Modify: `Answer Keys/MLB Dashboard/mlb_dashboard.R` lines ~4108–4297 (the `tags$script(HTML(r"(...)"))` block tagged `Wagerzon multi-account header bar — JS controller (Phase 6).`).

Five behavioral changes inside the IIFE:

1. Remove all `<select>`-related code (`SELECT_ID` constant, `renderSelect()` function, `sel.addEventListener('change', ...)`).
2. Rewrite `renderPills()` to emit class-based pills with click handlers. Drop the old inline `style.cssText` blob.
3. After balances load, if `WZ_SELECTED_ACCOUNT` is `null`, default to the first label and persist via `POST /api/wagerzon/last-used`. (The old `<select>` got this for free; the new div-based renderer must do it explicitly.)
4. Fix the "stale 0s ago" bug: the stale suffix only renders when `stale_seconds >= 60`.
5. Render an empty-state pill (`.wz-pill.empty` with text `No Wagerzon accounts configured`) when `orderedLabels.length === 0`.

- [ ] **Step 3.1: Locate the JS controller block**

```bash
grep -n "Wagerzon multi-account header bar — JS controller" "Answer Keys/MLB Dashboard/mlb_dashboard.R"
```

Expected: one line in the comment header, around line 4109.

- [ ] **Step 3.2: Replace the entire IIFE with the rewritten version**

Use the Edit tool with a large `old_string` covering everything from `tags$script(HTML(r"(` through the closing `)"))` (the line right before `    )` and `  )` that close out the page).

Replace with:

```r
      tags$script(HTML(r"(
(function() {
  // Wagerzon multi-account header pill row controller (Phase 7).
  // Owns:
  //   * window.WZ_SELECTED_ACCOUNT (string, source-of-truth selection)
  //   * window.WZ_BALANCES         (label -> snapshot from /api/wagerzon/balances)
  //   * window.wzApplyBalanceAfter (called by placeParlay() after a successful
  //                                 placement to refresh the pill without a full GET)
  //   * window._wzRecomputeWarnings (sweeps every Place button with data-risk and
  //                                  populates the sibling .wz-insufficient-warning)
  var PILLS_ID  = 'wz-account-pills';
  var REFRESH_BTN_ID = 'wz-refresh-btn';

  window.WZ_SELECTED_ACCOUNT = null;
  window.WZ_BALANCES = {};   // label -> snapshot

  function fmtMoney(n) {
    if (n === null || n === undefined) return '—';  // em dash
    return '$' + Number(n).toLocaleString('en-US', {minimumFractionDigits:2, maximumFractionDigits:2});
  }

  function persistSelection(label) {
    return fetch('/api/wagerzon/last-used', {
      method: 'POST',
      headers: {'Content-Type': 'application/json'},
      body: JSON.stringify({label: label})
    });
  }

  function pillTextFor(label, snap) {
    // Healthy: "Label · $1,234.56"
    // Errored, fresh:  "Label · — ⚠"
    // Errored, stale:  "Label · — ⚠ (stale Nm ago)"
    var text = label + ' · ' + fmtMoney(snap.available);
    if (snap.error) {
      text += ' ⚠';
      // Only label as "stale" when it has been at least a minute since the
      // last successful fetch. Avoids the contradictory "stale 0s ago" text
      // we used to render on every fresh-but-errored fetch.
      if (snap.stale_seconds >= 60) {
        var sub = Math.floor(snap.stale_seconds / 60) + 'm';
        text += ' (stale ' + sub + ' ago)';
      }
    }
    return text;
  }

  function isStale(snap) {
    return snap.error && snap.stale_seconds >= 60;
  }

  function renderEmpty(bar) {
    var pill = document.createElement('span');
    pill.className = 'wz-pill empty';
    pill.textContent = 'No Wagerzon accounts configured';
    bar.appendChild(pill);
  }

  function renderPills(orderedLabels) {
    var bar = document.getElementById(PILLS_ID);
    if (!bar) return;
    bar.innerHTML = '';

    if (orderedLabels.length === 0) {
      renderEmpty(bar);
      return;
    }

    orderedLabels.forEach(function(label) {
      var snap = WZ_BALANCES[label] || {available: null, error: 'no_data', stale_seconds: 0};
      var pill = document.createElement('span');
      var classes = ['wz-pill'];
      if (label === window.WZ_SELECTED_ACCOUNT) classes.push('selected');
      if (isStale(snap)) classes.push('stale');
      pill.className = classes.join(' ');
      pill.dataset.label = label;
      pill.textContent = pillTextFor(label, snap);

      pill.addEventListener('click', function() {
        if (label === window.WZ_SELECTED_ACCOUNT) return;
        window.WZ_SELECTED_ACCOUNT = label;
        persistSelection(label).catch(function() {
          // Network blip on persist shouldn't block UI; selection still
          // applies for this session and will retry on next click.
        });
        renderPills(orderedLabels);
        recomputeAllInsufficiencyWarnings();
      });

      bar.appendChild(pill);
    });
  }

  function refreshBalances() {
    return fetch('/api/wagerzon/balances')
      .then(function(r) { return r.json(); })
      .then(function(payload) {
        var orderedLabels = [];
        WZ_BALANCES = {};
        (payload.balances || []).forEach(function(snap) {
          WZ_BALANCES[snap.label] = snap;
          orderedLabels.push(snap.label);
        });

        // Default-to-first behaviour. The old <select> got this for free
        // via browser default-first-option; the div-based pill row needs
        // an explicit default + persist so the next page load is stable.
        if (!window.WZ_SELECTED_ACCOUNT && orderedLabels.length > 0) {
          window.WZ_SELECTED_ACCOUNT = orderedLabels[0];
          persistSelection(orderedLabels[0]).catch(function() {});
        }

        renderPills(orderedLabels);
        recomputeAllInsufficiencyWarnings();
      });
  }

  function loadLastUsed() {
    return fetch('/api/wagerzon/last-used')
      .then(function(r) { return r.json(); })
      .then(function(payload) {
        window.WZ_SELECTED_ACCOUNT = payload.label || null;
      });
  }

  function recomputeAllInsufficiencyWarnings() {
    if (typeof window._wzRecomputeWarnings === 'function') {
      window._wzRecomputeWarnings();
    }
  }

  // Public hook used by placeParlay() after a successful placement
  // to update the pill without a full re-fetch.
  window.wzApplyBalanceAfter = function(label, snap) {
    if (!snap) return;
    WZ_BALANCES[label] = snap;
    // Caller doesn't have orderedLabels — derive from the current DOM
    // so render order is preserved.
    var bar = document.getElementById(PILLS_ID);
    var labels = [];
    if (bar) {
      bar.querySelectorAll('.wz-pill[data-label]').forEach(function(el) {
        labels.push(el.dataset.label);
      });
    }
    if (labels.length === 0) labels = Object.keys(WZ_BALANCES);
    renderPills(labels);
    recomputeAllInsufficiencyWarnings();
  };

  // Per-parlay insufficient-balance warning. Reads each Place button's
  // data-risk against the currently-selected account's available balance
  // and writes a short message into the sibling .wz-insufficient-warning
  // span. Suppresses output when the account has no successful fetch yet
  // (available === null) to avoid showing misleading text.
  window._wzRecomputeWarnings = function() {
    var label = window.WZ_SELECTED_ACCOUNT;
    var snap  = label ? WZ_BALANCES[label] : null;
    var available = snap ? snap.available : null;

    document.querySelectorAll('button[data-hash][data-risk]').forEach(function(btn) {
      var risk = parseFloat(btn.dataset.risk);
      var warn = btn.parentElement.querySelector('.wz-insufficient-warning');
      if (!warn) return;

      if (available === null || isNaN(risk)) {
        warn.textContent = '';
        return;
      }
      if (risk > available) {
        warn.textContent = '⚠ insufficient on ' + label +
          ' (' + fmtMoney(available) + ' < ' + fmtMoney(risk) + ' risk)';
      } else {
        warn.textContent = '';
      }
    });
  };

  document.addEventListener('DOMContentLoaded', function() {
    var refreshBtn = document.getElementById(REFRESH_BTN_ID);
    if (refreshBtn) {
      refreshBtn.addEventListener('click', refreshBalances);
    }

    loadLastUsed().catch(function() {
      // Network blip on the last-used GET shouldn't block the rest of the
      // initial render. Pills will still populate; default-to-first kicks
      // in inside refreshBalances if WZ_SELECTED_ACCOUNT is still null.
    }).then(refreshBalances);

    // Watch for parlay-table re-renders (reactable pagination, hot-swap
    // after combined placement, etc.) so the warning is recomputed for
    // the new rows. Same pattern as the existing same-game observer on
    // bets-table-container.
    var parlayContainer = document.getElementById('parlays-table-container');
    if (parlayContainer && typeof MutationObserver === 'function') {
      var _wzWarnDebounce = null;
      var obs = new MutationObserver(function() {
        if (_wzWarnDebounce) clearTimeout(_wzWarnDebounce);
        _wzWarnDebounce = setTimeout(function() {
          _wzWarnDebounce = null;
          recomputeAllInsufficiencyWarnings();
        }, 100);
      });
      obs.observe(parlayContainer, {childList: true, subtree: true});
    }
  });
})();
)"))
```

Key changes vs. the previous controller:
- `SELECT_ID` constant removed.
- `renderSelect()` removed entirely.
- `renderPills()` no longer takes a separate `orderedLabels` parameter from `refreshBalances` only; the new signature `renderPills(orderedLabels)` is called from both `refreshBalances` and `wzApplyBalanceAfter`. Click handlers are added inline.
- `pillTextFor()` and `isStale()` extracted as named helpers — the stale-suffix gate `stale_seconds >= 60` is now in one place.
- `wzApplyBalanceAfter` derives labels from current DOM (so render order is stable across in-place updates).
- `loadLastUsed()` failure no longer leaves the dashboard with no selection — `refreshBalances` defaults to the first label and persists.

- [ ] **Step 3.3: Regenerate and inspect HTML**

```bash
Rscript "Answer Keys/MLB Dashboard/mlb_dashboard.R"
grep -c "renderSelect" "Answer Keys/MLB Dashboard/report.html"
grep -c "wz-account-select" "Answer Keys/MLB Dashboard/report.html"
grep -c "stale 0s ago" "Answer Keys/MLB Dashboard/report.html"
grep -c "stale_seconds >= 60" "Answer Keys/MLB Dashboard/report.html"
```

Expected:
- `renderSelect` → `0`
- `wz-account-select` → `0`
- `stale 0s ago` → `0` (the literal old text is gone — we now compose the suffix conditionally)
- `stale_seconds >= 60` → `2` (one in `pillTextFor`, one in `isStale`)

---

## Task 4: Manual visual verification

This is the gate before commit. The R/CSS/JS edits are in place; now we drive the page in a browser and confirm each manual test case from the spec.

- [ ] **Step 4.1: Hard-refresh the dashboard**

In the browser at http://localhost:8083, press `Cmd+Shift+R` (Mac) to force a fresh fetch of CSS/JS. The page should load without console errors.

- [ ] **Step 4.2: Verify alignment**

The "Placing on" caption + pills should sit directly under "Updated …", aligned with the same left edge as `MLB Answer Key Dashboard`. The green "Refresh" button is on the right of the title row. The tab bar (`Bets / Parlays / Trifectas`) sits below the pill row. No full-bleed bar above the title.

- [ ] **Step 4.3: Verify error rendering (current live state)**

Each pill should read `Label · — ⚠` with no "(stale 0s ago)" suffix. If you wait long enough that `stale_seconds >= 60`, the suffix appears as `(stale 1m ago)` and the pill goes red-tinted (`.wz-pill.stale`). The selected pill is filled blue.

- [ ] **Step 4.4: Verify click-to-select**

Click an inactive pill. It should:
- become filled blue
- prior selected pill returns to the inactive (dark grey) state
- a network request fires to `POST /api/wagerzon/last-used` (visible in DevTools Network tab) with the new label

Refresh the page. The clicked pill is still filled blue (persistence works).

- [ ] **Step 4.5: Verify the ↻ refresh button**

Click the small ↻ icon next to the rightmost pill. A `GET /api/wagerzon/balances` fires. Pills re-render with the same content (since the underlying `wz_error` state is unchanged).

- [ ] **Step 4.6: Verify the green "Refresh" button still works**

Click the dashboard-wide green "Refresh" button on the right of the title row. The dashboard data refreshes (pipeline run kicks off — usual `refreshData()` behavior).

- [ ] **Step 4.7: Verify the placeParlay flow**

On the Parlays tab, find any row with a Place button. The button's `data-risk` value vs. the active account's available balance should still produce or clear the `.wz-insufficient-warning` text. (Underlying balances are `null` due to `wz_error`, so the helper short-circuits — that's expected behaviour, not a regression.) If a parlay-table hot-swap happens after a combined placement, warnings still recompute on the new rows.

- [ ] **Step 4.8: Console check**

DevTools console should be clean — no `Cannot read properties of null` errors from the removed `<select>` references.

If any of 4.1–4.8 fail, fix in the R file, re-run `Rscript`, hard-refresh, retest. Do not commit until every case passes.

---

## Task 5: Update README

**Files:**
- Modify: `Answer Keys/MLB Dashboard/README.md`

The README currently has no section describing the multi-account header bar. Add a short new section between **Auto-Placement of Correlated Parlays** and **Troubleshooting** describing the in-header pill row. CLAUDE.md requires docs to ship in the same commit as the feature.

- [ ] **Step 5.1: Insert the new section**

Use the Edit tool. `old_string` matches the line right before the Troubleshooting heading; `new_string` adds the new section followed by the original anchor line.

`old_string`:

```
## Troubleshooting
```

`new_string`:

```
## Wagerzon account selector

The dashboard header includes a row of pills — one per configured Wagerzon
account (discovered by `wagerzon_odds/wagerzon_accounts.py`). Each pill shows
the account label and current available balance. Click a pill to switch the
active placement account; the selection is persisted to
`dashboard_settings.wagerzon_last_used` via `POST /api/wagerzon/last-used`
and used by `POST /api/place-parlay`.

- The selected pill is filled blue. Click an inactive pill (dark grey) to
  switch.
- A pill rendered as `Label · — ⚠` indicates the latest balance fetch
  failed but is still considered fresh (under one minute old). After one
  minute, it gains a `(stale Nm ago)` suffix and turns red-tinted.
- The small `↻` icon next to the pills refetches balances on demand. The
  green **Refresh** button on the right of the title row re-runs the
  dashboard pipeline.
- With zero accounts configured, the row renders a single dashed
  "No Wagerzon accounts configured" pill. Placement is disabled.

## Troubleshooting
```

- [ ] **Step 5.2: Confirm the edit**

```bash
grep -n "Wagerzon account selector" "Answer Keys/MLB Dashboard/README.md"
grep -n "## Troubleshooting" "Answer Keys/MLB Dashboard/README.md"
```

Expected: the new section header appears immediately before the Troubleshooting header.

---

## Task 6: Commit

**Files:**
- All changes from Tasks 1–3 (`Answer Keys/MLB Dashboard/mlb_dashboard.R`)
- README update from Task 5 (`Answer Keys/MLB Dashboard/README.md`)

One commit per CLAUDE.md (docs ship with the feature).

- [ ] **Step 6.1: Confirm the diff is what you expect**

```bash
git -C /Users/callancapitolo/NFLWork/.worktrees/mlb-dash-header-merge status --short
git -C /Users/callancapitolo/NFLWork/.worktrees/mlb-dash-header-merge diff --stat
```

Expected: only `Answer Keys/MLB Dashboard/mlb_dashboard.R` and `Answer Keys/MLB Dashboard/README.md` modified. No stray files.

- [ ] **Step 6.2: Stage and commit**

```bash
cd /Users/callancapitolo/NFLWork/.worktrees/mlb-dash-header-merge
git add "Answer Keys/MLB Dashboard/mlb_dashboard.R" "Answer Keys/MLB Dashboard/README.md"
git commit -m "$(cat <<'EOF'
feat(mlb-dashboard): merge WZ account bar into header pill row

Replace the standalone Wagerzon multi-account header bar (full-bleed,
sat above .container) with a click-to-select pill row inside the
existing dashboard header. The dropdown selector is removed; pills
double as the selector and use real CSS classes (.wz-pill,
.wz-pill.selected, .wz-pill.stale, .wz-pill.empty) instead of
inline-styled spans. The "(stale 0s ago)" rendering bug is fixed —
the stale suffix only renders when stale_seconds >= 60.

Layout-only change. No backend, schema, or API contract changes.
README updated to reference the new in-header pill row.

See docs/superpowers/specs/2026-05-02-mlb-dashboard-header-merge-design.md.

Co-Authored-By: Claude Opus 4.7 (1M context) <noreply@anthropic.com>
EOF
)"
```

- [ ] **Step 6.3: Verify commit landed**

```bash
git -C /Users/callancapitolo/NFLWork/.worktrees/mlb-dash-header-merge log --oneline -3
```

Expected: top commit is the new `feat(mlb-dashboard): merge WZ account bar into header pill row`. Below it: `1b55f96 docs(mlb-dashboard): design — merge WZ account bar into header pill row`.

---

## Task 7: Pre-merge review (per CLAUDE.md)

Before asking for merge approval, run the pre-merge checklist from `CLAUDE.md`. This is review, not implementation — read the diff against the checklist and report findings.

- [ ] **Step 7.1: Generate the full diff**

```bash
git -C /Users/callancapitolo/NFLWork/.worktrees/mlb-dash-header-merge diff main..HEAD -- "Answer Keys/MLB Dashboard/" > /tmp/header-merge.diff
wc -l /tmp/header-merge.diff
```

- [ ] **Step 7.2: Walk the checklist**

For each item, write either ✅ (no issue) or ⚠️ (issue + line numbers in the diff). Items:

- **Data integrity**: This change does not write to any DB. No risk of duplicate writes / dedup issues.
- **Resource safety**: No new DB connections. Existing `on.exit(dbDisconnect(...))` calls in `mlb_dashboard.R` are untouched.
- **Edge cases**: 0 accounts, 1 account, all errored, fresh-but-errored, and stale-errored states are all enumerated in the spec and exercised in Task 4 manual tests.
- **Dead code**: `renderSelect()`, `SELECT_ID`, the `<select>` change-listener — confirm all three are gone (`grep -n` for each in the diff).
- **Log/disk hygiene**: Layout change only; no new logs or disk writes.
- **Security**: No secrets or API keys. `placeParlay` still uses the same `account` field; no exposure changes.

- [ ] **Step 7.3: Document findings**

Post the checklist results inline (or as a comment on the eventual PR/merge approval). If any ⚠️ items: fix and re-run Task 6 + 7. If all ✅: proceed to Task 8.

---

## Task 8: Merge (only with explicit user approval)

Per CLAUDE.md: "Never merge to `main` or push to remote without explicit user approval, even if tests pass." Ask the user "Pre-merge review clean. Merge to main now?" and wait for an explicit "yes" before running these steps.

- [ ] **Step 8.1: Switch back to main and merge**

```bash
cd /Users/callancapitolo/NFLWork
git checkout main
git pull --ff-only origin main 2>/dev/null || true   # optional, only if remote exists
git merge --no-ff feature/mlb-dashboard-header-merge
```

The `--no-ff` matches the project's recent merge style (e.g. `b31f7fb`, `864542b`).

- [ ] **Step 8.2: Regenerate dashboard from main and re-verify in browser**

```bash
Rscript "Answer Keys/MLB Dashboard/mlb_dashboard.R"
```

Hard-refresh http://localhost:8083 — confirm the merged page renders identically to what you tested in the worktree.

- [ ] **Step 8.3: Clean up worktree + branch**

```bash
git worktree remove .worktrees/mlb-dash-header-merge
git branch -d feature/mlb-dashboard-header-merge
git worktree list
```

Expected: the worktree is gone; only `main` (and any other active worktrees) remain.

---

## Self-review notes

Before handing this plan to an executor:

**Spec coverage** — every section of the spec has a task:
- Problem statement → addressed by the solution as a whole
- DOM layout → Task 2
- CSS classes → Task 1
- Pill behaviour table → Task 3 (renderer + click handler)
- "stale 0s ago" fix → Task 3 (`pillTextFor` + `isStale`)
- Default-to-first-account → Task 3 (`refreshBalances`)
- Manual test plan (9 cases) → Task 4 (8 sub-steps; cases collapsed where overlap exists)
- Files Affected (`mlb_dashboard.R`, `README.md`) → Tasks 1–3 + Task 5
- Documentation discipline → Task 5 (same commit as feature in Task 6)
- Worktree lifecycle → Task 8

**Placeholder scan** — no `TBD`, `TODO`, "implement later", or vague "handle edge cases" instructions. Every code step has the actual code. Every command has expected output.

**Type/identifier consistency** — `PILLS_ID = 'wz-account-pills'`, `REFRESH_BTN_ID = 'wz-refresh-btn'`, classes `.wz-pill`, `.wz-pill.selected`, `.wz-pill.stale`, `.wz-pill.empty`, `.wz-icon-btn`, `.header-row-top`, `.header-row-accounts`, `.header-label` — all defined in Task 1, used consistently in Tasks 2–3.
