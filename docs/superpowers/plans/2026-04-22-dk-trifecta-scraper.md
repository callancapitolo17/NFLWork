# DK Trifecta Production Scraper (Plan #2) Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Activate the DK SGP blend that Plan #1 scaffolded — populate `mlb_trifecta_sgp_odds` with real DK odds for posted triple-plays and grand slams, so `mlb_triple_play.R` produces blended fair odds (not model-only) on each run.

**Architecture:** New Python module `mlb_sgp/dk_leg_resolvers.py` defines `LEG_RESOLVERS` mapping each leg type from `parse_legs.R` to a function that finds the matching DK selection ID inside the cached SGP event payload. New module `mlb_sgp/scraper_draftkings_trifecta.py` is a CLI: reads a JSON request file (one entry per (game, prop_type, side, legs)), resolves each leg's DK selection ID, posts to `calculateBets`, and writes per-prop SGP odds rows to `mlb_trifecta_sgp_odds`. R-side activation is a one-line uncomment of the `system2()` call already stubbed in `Answer Keys/mlb_triple_play.R` (Plan #1 commit `005e4bc`). Path A is locked: the recon-confirmed primitives (`1st Run`, `1st 5 Innings`, `Moneyline`, `<TEAM>: Team Total Runs`) all exist and are SGP-eligible.

**Tech Stack:** Python 3 (curl_cffi for DK requests + Chrome impersonation, duckdb for writes — both already vendored in `mlb_sgp/venv`). pytest for resolver unit tests. R (existing pricer; only one comment-block uncomment).

---

## Recon-confirmed primitives (verified 2026-04-22)

DK markets we resolve to:

| Leg type from `parse_legs.R` | DK market name (exact) | Sel ID pattern | Selection count |
|---|---|---|---|
| `scores_first` | `"1st Run"` | `0QA<grp>#<sid>_<suffix>` | 2 (home / away) |
| `wins_period(FG)` | `"Moneyline"` | `0ML<num>_<suffix>` | 2 (home / away) |
| `wins_period(F5)` | `"1st 5 Innings"` | `0QA<grp>#<sid>_<suffix>` | 2 (home / away) |
| `team_total_under/over <N>` | `"<TEAM>: Team Total Runs"` (and alts) | `0OU<num><U or O><line×100>_<suffix>` | 2 per line (over / under) |

**All four are tagged `'SGP'` in the DK event payload's `tags` list, meaning they're SGP-eligible.** That's the precondition for `calculateBets` to accept them as legs in a combined parlay.

`wins_period(F3)` and `wins_period(F7)` are reserved by `TOKEN_REGISTRY` but currently unused by today's posted props. The resolver returns `None` for those leg types in v1 — the scraper writes `NULL sgp_decimal` for any prop containing an unsupported leg, and the R blend correctly degrades to model-only for that row.

---

## File Structure

**Created:**
- `mlb_sgp/dk_leg_resolvers.py` — `LEG_RESOLVERS` dict + 4 resolver functions (one per supported leg type) + a small `find_market_by_name()` helper. ~150 lines.
- `mlb_sgp/scraper_draftkings_trifecta.py` — CLI entry point. Reads `--input <json>`, fetches per-game SGP payloads (cached in-memory), resolves legs to selection IDs, posts to `calculateBets`, writes rows to `mlb_trifecta_sgp_odds`. ~180 lines.
- `mlb_sgp/tests/test_dk_leg_resolvers.py` — pytest fixtures from a captured recon payload, one test per leg type plus negative tests (unknown leg type, missing market). ~80 lines.

**Modified:**
- `Answer Keys/mlb_triple_play.R` — uncomment the Plan #2 stub (the `# trifecta_input <- ...` block from commit `005e4bc`). Add a small bridge that writes the trifecta input JSON, invokes the scraper via `system2()`, and reads back the populated table. ~25 lines added/modified.
- `Answer Keys/CLAUDE.md` — replace the "Plan #2 (post-recon) populates …" hedge with the live data-flow description.

**Auxiliary** (not in commits):
- `mlb_sgp/recon_dk_trifecta_<event_id>.json` already gitignored; tests use a small embedded fixture, not the full payload.

---

## Worktree & Version Control Plan

- Branch: `feature/dk-trifecta-scraper`
- Worktree: `.worktrees/dk-trifecta-scraper`
- **Setup (before Task 1):**
  ```bash
  cd /Users/callancapitolo/NFLWork
  git worktree add .worktrees/dk-trifecta-scraper -b feature/dk-trifecta-scraper main
  cd .worktrees/dk-trifecta-scraper
  git branch  # confirm feature/dk-trifecta-scraper
  ```
- **Commits** (one per task):
  1. Task 1 → `feat(dk): verify Path A 3-leg SGP combinability`
  2. Task 2 → `feat(dk): add LEG_RESOLVERS for trifecta selection-id resolution`
  3. Task 3 → `feat(dk): add scraper_draftkings_trifecta CLI + DuckDB writer`
  4. Task 4 → `test(dk): resolver unit tests against captured recon fixture`
  5. Task 5 → `feat(mlb): activate DK trifecta scraper invocation in pricer`
  6. Task 6 → `docs(mlb): document Plan #2 live data flow`
- **DuckDB:** `wagerzon_odds/wagerzon.duckdb` and `Answer Keys/mlb.duckdb` are NOT touched in the worktree (per project rule). Task 5's pricer regression copies `mlb.duckdb` in, runs, then `rm`s. The `mlb_trifecta_sgp_odds` table is created lazily by Plan #1's `CREATE IF NOT EXISTS`; Plan #2's scraper writes into it.
- **Cleanup after merge:** `git worktree remove .worktrees/dk-trifecta-scraper && git branch -d feature/dk-trifecta-scraper`
- **Never merge to main without explicit user approval.**

---

## Documentation Plan

Task 6 updates `Answer Keys/CLAUDE.md`'s "Triple-Play Data Flow" subsection: change the trailing bullet `Plan #2 (post-recon) populates ...` from a hedge to a live description. Also update the data-flow diagram so the `dk_sgp lookup → devig` step now references the production scraper rather than describing it as Plan #1 stubbed.

No `mlb_sgp/README.md` update — the existing scraper README pattern in that directory is to leave new scrapers undocumented in README and rely on docstrings + the production CLAUDE.md. Stay consistent.

---

## Pre-Merge Review Checklist

- **Path A combinability stays verified**: Task 1's smoke test must pass each pre-merge run. If DK rejects the 3-leg combination, the entire plan is suspect.
- **Resolver fixture coverage**: every leg type Plan #1 supports has a positive test + a negative (missing market) test.
- **DK event-id resolution accuracy**: the scraper joins `mlb_consensus_temp.id` (Odds API hex) → DK event id via team-name + commence_time match. Run on tonight's slate; assert all `wagerzon_specials` rows resolve to a DK event id (no silent NA).
- **Atomic snapshot**: scraper writes all rows for one run with the same `fetch_time`. The pricer's `MAX(fetch_time)` query relies on this.
- **Resource safety**: scraper opens one curl_cffi session, one DuckDB connection, both closed in `finally`. No leaks if any per-game fetch errors.
- **Empty/error handling**: per-game fetch failure → row skipped with WARN log, scraper continues; per-leg resolution failure → row written with NULL `sgp_decimal`, R blend degrades correctly.
- **No secrets / no PII**: scraper logs URL + status + game id only. No headers, no response bodies.
- **R-side regression**: model_odds for previously-priced props must be identical to pre-merge. The blend may push fair_odds away from model_odds (that's the whole point), but model_odds itself is independent of DK and must not drift.
- **Live regression**: Task 6's smoke test prints the 14-row table with at least 8/14 rows showing populated `dk_odds`. Some rows may legitimately have NA (DK-side missing market or unresolvable team). Document which rows are NA and why.

---

## Task 1: Verify Path A 3-leg SGP combinability

**Status: Critical — if DK rejects Path A's 3-leg combination, the plan stops here and we revisit Path B.**

**Files:**
- Create: `mlb_sgp/scratch_test_combinability.py` (DELETED at end of task — scratch artifact)

### - [ ] Step 1: Write the smoke test

Create `mlb_sgp/scratch_test_combinability.py`:

```python
"""
Smoke-test Path A combinability — does DK's calculateBets accept a 3-leg
SGP composed of (1st Run, 1st 5 Innings, Moneyline) for a single team?

If this fails, Plan #2 must pivot to Path B (use the 'Score First + Win'
prebuilt instead of the primitive 1st Run market).

This file is DELETED at the end of Task 1 — scratch artifact, not committed.
"""
import sys, json
sys.path.insert(0, '/Users/callancapitolo/NFLWork/mlb_sgp')
from scraper_draftkings_sgp import init_session
from curl_cffi import requests as cffi_req

DK_CALCULATE_BETS_URL = "https://gaming-us-nj.draftkings.com/api/wager/v1/calculateBets"

def find_market_first_match(payload, target_name):
    """Walk the payload; return (market_node, selections list) for first
    market whose name matches target_name and has selections."""
    def walk(node):
        if isinstance(node, dict):
            n = node.get('marketName') or node.get('name')
            if n == target_name and node.get('selections'):
                return node
            for v in node.values():
                r = walk(v)
                if r is not None: return r
        elif isinstance(node, list):
            for v in node:
                r = walk(v)
                if r is not None: return r
        return None
    return walk(payload)

def main():
    s = init_session()
    # Pick a future game (CHI Cubs @ SD Padres at 01:40 UTC tomorrow per recon list)
    event_id = '34055516'
    print(f'Fetching SGP event {event_id} ...')
    url = f'https://sportsbook-nash.draftkings.com/sites/US-SB/api/sportscontent/parlays/v1/sgp/events/{event_id}'
    r = s.get(url, timeout=30)
    r.raise_for_status()
    payload = r.json()

    # Find selection IDs for SD Padres on each of the 3 primitives
    legs = []
    for market_name in ['1st Run', '1st 5 Innings', 'Moneyline']:
        market = find_market_first_match(payload, market_name)
        assert market, f"Market {market_name!r} not found in payload"
        sels = market['selections']
        # Pick the away-team selection (SD Padres in this game) by name
        sd_sel = next((s_ for s_ in sels if 'Padres' in (s_.get('name') or '')), None)
        assert sd_sel, f"No SD Padres selection in {market_name}"
        sel_id = sd_sel.get('id') or sd_sel.get('selectionId')
        legs.append({'name': market_name, 'sel_id': sel_id})
        print(f'  {market_name!r:25s}  {sel_id}')

    # Build calculateBets request with all 3 legs
    body = {
        'bets': [{
            'wagerOption': 'StraightSGP',
            'selections': [{'id': leg['sel_id']} for leg in legs],
        }],
    }
    print(f'\nPOST calculateBets with {len(legs)} legs ...')
    resp = s.post(DK_CALCULATE_BETS_URL, json=body, timeout=30)
    print(f'  status={resp.status_code}')
    if resp.status_code != 200:
        print(f'  body: {resp.text[:500]}')
        print('\nPATH A REJECTED — DK refuses the 3-leg primitive combination.')
        print('Pivot Plan #2 to Path B (use "Score First + Win" prebuilt + F5 ML).')
        sys.exit(2)

    data = resp.json()
    print(f'  response keys: {list(data.keys())}')
    # Try to find the SGP odds in the response
    bets_arr = data.get('bets') or data.get('result') or []
    print(f'  bets: {json.dumps(bets_arr, default=str)[:600]}')

    print('\nPATH A CONFIRMED — DK accepted the 3-leg primitive SGP. Proceed with Plan #2.')

if __name__ == '__main__':
    main()
```

### - [ ] Step 2: Run the smoke test

```bash
cd /Users/callancapitolo/NFLWork/.worktrees/dk-trifecta-scraper
/Users/callancapitolo/NFLWork/mlb_sgp/venv/bin/python mlb_sgp/scratch_test_combinability.py
```

Expected outcomes:
- **Path A confirmed**: prints "PATH A CONFIRMED — DK accepted the 3-leg primitive SGP" and exits 0. Proceed to Step 3.
- **Path A rejected**: prints "PATH A REJECTED" and exits 2. STOP and report BLOCKED to the user with the response body. Plan #2 needs to be redesigned around Path B before continuing.

### - [ ] Step 3: Capture combinability output for the commit log

Save the output of Step 2 to a comment block at the top of the next commit message — this gives us a permanent record of when Path A was last verified, what the SGP decimal was, and which game/legs were tested.

### - [ ] Step 4: Delete the scratch file (only if Path A confirmed)

```bash
rm mlb_sgp/scratch_test_combinability.py
```

### - [ ] Step 5: Commit

```bash
cd /Users/callancapitolo/NFLWork/.worktrees/dk-trifecta-scraper
git commit --allow-empty -m "feat(dk): verify Path A 3-leg SGP combinability

Smoke test confirmed DK's calculateBets accepts a 3-leg primitive SGP
composed of (1st Run, 1st 5 Innings, Moneyline) for a single team.
Tested on CHI Cubs @ SD Padres (event 34055516), SD Padres side.
[paste actual SGP decimal returned, e.g. SGP decimal: 4.50, ~+350]

This is the precondition for Plan #2's resolver-driven scraper.
If a future combinability test fails, the plan must pivot to Path B."
```

(`--allow-empty` because Step 4 deleted the only file added in this task. The commit is a marker.)

---

## Task 2: LEG_RESOLVERS — selection-ID resolution

**Status: Depends on Task 1 (Path A confirmed).**

**Files:**
- Create: `mlb_sgp/dk_leg_resolvers.py`

### - [ ] Step 1: Create the resolver module

Create `mlb_sgp/dk_leg_resolvers.py`:

```python
"""
DK leg-type → selection-id resolvers for Plan #2's trifecta scraper.

Each resolver takes:
  - leg: a leg-spec dict matching parse_legs.R::TOKEN_REGISTRY output, e.g.:
      {'type': 'scores_first'}
      {'type': 'wins_period', 'period': 'F5'}
      {'type': 'team_total_under', 'line': 2.5}
  - side: 'home' or 'away' — determines which selection within the market we pick
  - event_state: the cached DK SGP event payload (full JSON returned by
    parlays/v1/sgp/events/{event_id})
  - team_names: dict {'home': 'Cubs', 'away': 'Padres'} — DK team-name strings
    extracted from event_state at scrape time so resolvers can match selections
    by name (DK's selection.name field is e.g. 'CHI Cubs' or 'SD Padres').

Each resolver returns:
  - str selection_id when the market+selection is found and SGP-eligible
  - None when the market is missing OR not SGP-eligible (R-side blend
    degrades to model-only for that row)
"""
from __future__ import annotations
from typing import Optional


# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------

def find_market_by_name(payload: dict, target_name: str) -> Optional[dict]:
    """Walk payload, return first market dict whose name matches target_name
    AND is SGP-eligible (has 'SGP' in tags) AND has selections.
    """
    def walk(node):
        if isinstance(node, dict):
            n = node.get('marketName') or node.get('name')
            if (
                n == target_name
                and node.get('selections')
                and 'SGP' in (node.get('tags') or [])
            ):
                return node
            for v in node.values():
                r = walk(v)
                if r is not None:
                    return r
        elif isinstance(node, list):
            for v in node:
                r = walk(v)
                if r is not None:
                    return r
        return None
    return walk(payload)


def find_team_total_market(
    payload: dict, team_label: str, side_target: str, line: float
) -> Optional[dict]:
    """Find a "<TEAM>: Team Total Runs" market with a selection matching the
    given line + over/under. Returns the selection dict (with 'id'), or None.

    Multiple team-total markets exist per game (full-game, alt lines, by-inning).
    We want the FULL-game one for the given team and line.
    """
    target_market_name = f"{team_label}: Team Total Runs"
    market = find_market_by_name(payload, target_market_name)
    if market is None:
        return None
    # Selections should have labels like "Over" / "Under" plus a 'point' field
    for sel in market.get('selections', []):
        sel_label = (sel.get('label') or sel.get('name') or '').lower()
        sel_point = sel.get('point')
        if (
            (side_target == 'over'  and 'over'  in sel_label) or
            (side_target == 'under' and 'under' in sel_label)
        ):
            try:
                if sel_point is not None and abs(float(sel_point) - line) < 0.01:
                    return sel
            except (TypeError, ValueError):
                continue
    return None


def _pick_team_selection(market: dict, team_label: str) -> Optional[str]:
    """Given a 2-way team market (e.g. '1st Run', 'Moneyline', '1st 5 Innings'),
    return the selection.id whose name matches team_label exactly.
    DK selection names are e.g. 'CHI Cubs', 'SD Padres'.
    """
    for sel in market.get('selections', []):
        if (sel.get('name') or '').strip() == team_label:
            return sel.get('id') or sel.get('selectionId')
    return None


# ---------------------------------------------------------------------------
# Per-leg-type resolvers
# ---------------------------------------------------------------------------

def resolve_scores_first(leg, side, event_state, team_names) -> Optional[str]:
    """leg = {'type': 'scores_first'}; side = 'home' | 'away'."""
    market = find_market_by_name(event_state, '1st Run')
    if market is None:
        return None
    return _pick_team_selection(market, team_names[side])


def resolve_wins_period(leg, side, event_state, team_names) -> Optional[str]:
    """leg = {'type': 'wins_period', 'period': 'F3'|'F5'|'F7'|'FG'}."""
    period = leg.get('period')
    market_name = {
        'FG': 'Moneyline',
        'F5': '1st 5 Innings',
        # F3 / F7 not posted as 2-way ML primitives at DK reliably — return None
    }.get(period)
    if market_name is None:
        return None
    market = find_market_by_name(event_state, market_name)
    if market is None:
        return None
    return _pick_team_selection(market, team_names[side])


def resolve_team_total_under(leg, side, event_state, team_names) -> Optional[str]:
    """leg = {'type': 'team_total_under', 'line': 2.5}."""
    line = leg.get('line')
    if line is None:
        return None
    sel = find_team_total_market(event_state, team_names[side], 'under', float(line))
    if sel is None:
        return None
    return sel.get('id') or sel.get('selectionId')


def resolve_team_total_over(leg, side, event_state, team_names) -> Optional[str]:
    """leg = {'type': 'team_total_over', 'line': 4.5}."""
    line = leg.get('line')
    if line is None:
        return None
    sel = find_team_total_market(event_state, team_names[side], 'over', float(line))
    if sel is None:
        return None
    return sel.get('id') or sel.get('selectionId')


# ---------------------------------------------------------------------------
# Registry — mirrors parse_legs.R::TOKEN_REGISTRY structure
# ---------------------------------------------------------------------------

LEG_RESOLVERS = {
    'scores_first':     resolve_scores_first,
    'wins_period':      resolve_wins_period,
    'team_total_under': resolve_team_total_under,
    'team_total_over':  resolve_team_total_over,
    # Reserved (parser knows them; resolver returns None until DK markets
    # are confirmed available):
    'opp_total_under':  lambda *a, **k: None,
    'opp_total_over':   lambda *a, **k: None,
}


def resolve_legs(legs, side, event_state, team_names) -> Optional[list[str]]:
    """Resolve a list of leg specs into DK selection IDs. Returns None if any
    leg fails to resolve (in which case the whole prop is unprice able by DK
    and the caller writes NULL sgp_decimal)."""
    sel_ids = []
    for leg in legs:
        leg_type = leg.get('type')
        resolver = LEG_RESOLVERS.get(leg_type)
        if resolver is None:
            return None
        sel_id = resolver(leg, side, event_state, team_names)
        if sel_id is None:
            return None
        sel_ids.append(sel_id)
    return sel_ids
```

### - [ ] Step 2: Syntax check

```bash
cd /Users/callancapitolo/NFLWork/.worktrees/dk-trifecta-scraper
python3 -c "
import ast
src = open('mlb_sgp/dk_leg_resolvers.py').read()
ast.parse(src)
print('SYNTAX OK')
"
```
Expected: `SYNTAX OK`.

### - [ ] Step 3: Import smoke test

```bash
cd /Users/callancapitolo/NFLWork/.worktrees/dk-trifecta-scraper
/Users/callancapitolo/NFLWork/mlb_sgp/venv/bin/python -c "
import sys
sys.path.insert(0, 'mlb_sgp')
from dk_leg_resolvers import LEG_RESOLVERS, resolve_legs, find_market_by_name
assert callable(resolve_legs)
assert 'scores_first' in LEG_RESOLVERS
assert 'wins_period' in LEG_RESOLVERS
assert 'team_total_under' in LEG_RESOLVERS
assert 'team_total_over' in LEG_RESOLVERS
print('IMPORT OK')
"
```
Expected: `IMPORT OK`.

### - [ ] Step 4: Commit

```bash
cd /Users/callancapitolo/NFLWork/.worktrees/dk-trifecta-scraper
git add mlb_sgp/dk_leg_resolvers.py
git commit -m "feat(dk): add LEG_RESOLVERS for trifecta selection-id resolution

Path A primitives:
- scores_first   → '1st Run' market (2-way home/away)
- wins_period FG → 'Moneyline' market
- wins_period F5 → '1st 5 Innings' market
- team_total_*   → '<TEAM>: Team Total Runs' market by line + side

All four are SGP-eligible per recon (2026-04-22). F3/F7 wins_period and
opp_total_* return None — markets aren't reliably present at DK and
aren't used by today's posted props.

resolve_legs() returns None if any leg fails to resolve, signaling the
caller to write NULL sgp_decimal so the R blend degrades to model-only."
```

---

## Task 3: Production scraper — `scraper_draftkings_trifecta.py`

**Status: Depends on Task 2.**

**Files:**
- Create: `mlb_sgp/scraper_draftkings_trifecta.py`

### - [ ] Step 1: Create the scraper

Create `mlb_sgp/scraper_draftkings_trifecta.py`:

```python
"""
DraftKings Trifecta SGP Scraper

Reads a JSON request file with one entry per (game_id, prop_type, side, legs),
resolves each leg into a DK selection ID via dk_leg_resolvers.LEG_RESOLVERS,
posts the combined SGP to DK's calculateBets endpoint, and writes one row
per request to mlb_trifecta_sgp_odds in the project's MLB DuckDB.

Usage:
    python3 scraper_draftkings_trifecta.py \\
        --input /tmp/trifecta_requests.json \\
        --db    "Answer Keys/mlb.duckdb"

Input JSON shape (one element per request):
    [
      {
        "game_id":   "0438fd942a9f1f07557cef989b1a4e4d",  # Odds API event id
        "home_team": "San Francisco Giants",
        "away_team": "Los Angeles Dodgers",
        "prop_type": "TRIPLE-PLAY",
        "side":      "home",
        "legs":      [
          {"type": "scores_first"},
          {"type": "wins_period", "period": "F5"},
          {"type": "wins_period", "period": "FG"}
        ]
      },
      ...
    ]

Exit codes:
    0  - success (some/all rows may have NULL sgp_decimal if DK couldn't price them)
    1  - usage error (bad CLI args or missing input file)
    2  - DK API failure (Akamai block, network, etc.)
"""
from __future__ import annotations

import argparse
import json
import logging
import sys
from datetime import datetime, timezone
from pathlib import Path
from typing import Optional

import duckdb

# Reuse the production DK SGP scraper's session + event-discovery functions
SCRIPT_DIR = Path(__file__).resolve().parent
sys.path.insert(0, str(SCRIPT_DIR))
from scraper_draftkings_sgp import init_session, fetch_dk_events  # noqa: E402
from dk_leg_resolvers import resolve_legs, find_market_by_name      # noqa: E402

# DK API URLs — same prefix the production scraper uses
SGP_EVENTS_URL_TPL = (
    "https://sportsbook-nash.draftkings.com/sites/US-SB/api/sportscontent/"
    "parlays/v1/sgp/events/{event_id}"
)
DK_CALCULATE_BETS_URL = "https://gaming-us-nj.draftkings.com/api/wager/v1/calculateBets"

logging.basicConfig(level=logging.INFO, format="%(asctime)s [%(levelname)s] %(message)s")
log = logging.getLogger("scraper_draftkings_trifecta")


def map_odds_id_to_dk_event(events: list[dict], home_team: str, away_team: str) -> Optional[str]:
    """Find the DK event whose teams match the requested home/away. DK uses
    e.g. 'CHI Cubs' / 'SD Padres'; we use Odds API canonical names. We match
    by checking if the canonical team's last word is in DK's team string
    (e.g. 'San Francisco Giants' matches 'SF Giants' on 'Giants')."""
    home_last = home_team.split()[-1].lower()
    away_last = away_team.split()[-1].lower()
    for evt in events:
        dk_home = (evt.get('dk_home') or '').lower()
        dk_away = (evt.get('dk_away') or '').lower()
        if home_last in dk_home and away_last in dk_away:
            return str(evt['dk_event_id'])
    return None


def fetch_event_state(session, event_id: str) -> Optional[dict]:
    """Fetch the full SGP event payload (~1-6 MB)."""
    url = SGP_EVENTS_URL_TPL.format(event_id=event_id)
    try:
        r = session.get(url, timeout=30)
        r.raise_for_status()
        return r.json()
    except Exception as e:
        log.error("Failed to fetch event %s: %s", event_id, e)
        return None


def extract_team_names_from_event(event_state: dict) -> dict:
    """Walk the event payload to extract the DK team-name strings used in
    selection labels (e.g. 'CHI Cubs', 'SD Padres'). Picks them from the
    Moneyline market which always exists and always has both teams."""
    market = find_market_by_name(event_state, 'Moneyline')
    if market is None:
        return {'home': '', 'away': ''}
    sels = market.get('selections', [])
    if len(sels) < 2:
        return {'home': '', 'away': ''}
    # DK Moneyline selections are ordered [home, away] OR [away, home] —
    # ID suffix _3 vs _1 historically encodes side. The production scraper
    # uses _3 = home, _1 = away. Apply same convention here.
    home_name = ''
    away_name = ''
    for sel in sels:
        sid = sel.get('id') or sel.get('selectionId') or ''
        nm = sel.get('name') or ''
        if sid.endswith('_3'):
            home_name = nm
        elif sid.endswith('_1'):
            away_name = nm
    return {'home': home_name, 'away': away_name}


def post_calculate_bets(session, sel_ids: list[str]) -> Optional[float]:
    """Post the combined SGP to DK's calculateBets endpoint. Returns the
    decimal odds, or None if DK rejected/errored."""
    body = {
        'bets': [{
            'wagerOption': 'StraightSGP',
            'selections': [{'id': sid} for sid in sel_ids],
        }],
    }
    try:
        r = session.post(DK_CALCULATE_BETS_URL, json=body, timeout=30)
        if r.status_code != 200:
            log.warning("calculateBets returned %d: %s", r.status_code, r.text[:200])
            return None
        data = r.json()
        # Walk the response for the trueOdds / decimalOdds field
        bets_arr = data.get('bets') or []
        if bets_arr and isinstance(bets_arr, list):
            bet = bets_arr[0]
            # Try common locations for the decimal odds
            odds = bet.get('trueOdds') or bet.get('decimal')
            if odds is None:
                # Fall back to displayOdds.decimal if present
                disp = bet.get('displayOdds') or {}
                if isinstance(disp, dict):
                    try:
                        odds = float(disp.get('decimal'))
                    except (TypeError, ValueError):
                        odds = None
            if odds:
                return float(odds)
        log.warning("Could not extract odds from calculateBets response: %s", str(data)[:300])
        return None
    except Exception as e:
        log.error("calculateBets POST failed: %s", e)
        return None


def decimal_to_american(dec: float) -> int:
    if dec >= 2.0:
        return int(round((dec - 1) * 100))
    return int(round(-100 / (dec - 1)))


def write_rows(con: duckdb.DuckDBPyConnection, rows: list[dict]) -> None:
    """Write all rows under one fetch_time so MAX(fetch_time) returns a
    consistent atomic snapshot."""
    if not rows:
        log.info("No rows to write.")
        return
    # CREATE IF NOT EXISTS — idempotent; Plan #1 created the table already
    # but this protects against edge cases (e.g. fresh checkout, test fixtures).
    con.execute("""
        CREATE TABLE IF NOT EXISTS mlb_trifecta_sgp_odds (
            fetch_time     TIMESTAMP,
            game_id        VARCHAR,
            prop_type      VARCHAR,
            side           VARCHAR,
            legs_json      VARCHAR,
            selection_ids  VARCHAR,
            sgp_decimal    DOUBLE,
            sgp_american   INTEGER,
            source         VARCHAR
        )
    """)
    con.executemany(
        "INSERT INTO mlb_trifecta_sgp_odds VALUES (?,?,?,?,?,?,?,?,?)",
        [
            (
                r['fetch_time'], r['game_id'], r['prop_type'], r['side'],
                r['legs_json'], r['selection_ids'],
                r['sgp_decimal'], r['sgp_american'], r['source'],
            )
            for r in rows
        ],
    )


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument('--input', required=True, help='Path to trifecta requests JSON')
    ap.add_argument('--db', required=True, help='Path to mlb.duckdb')
    args = ap.parse_args()

    input_path = Path(args.input)
    db_path = Path(args.db)
    if not input_path.exists():
        log.error("Input file not found: %s", input_path)
        return 1

    requests_in = json.loads(input_path.read_text())
    log.info("Loaded %d trifecta requests", len(requests_in))

    session = init_session()
    events = fetch_dk_events(session)
    log.info("DK eventgroup returned %d MLB events", len(events))

    # Cache event-state fetches and team-name extractions per DK event id
    event_state_cache: dict[str, dict] = {}
    team_names_cache: dict[str, dict] = {}

    fetch_time = datetime.now(timezone.utc).replace(microsecond=0)
    rows = []

    for req in requests_in:
        game_id   = req['game_id']      # Odds API hex id
        home_team = req['home_team']
        away_team = req['away_team']
        prop_type = req['prop_type']
        side      = req['side']
        legs      = req['legs']

        # Map Odds API id → DK event id
        dk_event_id = map_odds_id_to_dk_event(events, home_team, away_team)
        if dk_event_id is None:
            log.warning("No DK event for %s vs %s — skipping (game_id=%s prop=%s side=%s)",
                        home_team, away_team, game_id, prop_type, side)
            continue

        # Cached event-state fetch
        if dk_event_id not in event_state_cache:
            state = fetch_event_state(session, dk_event_id)
            if state is None:
                log.warning("Could not fetch DK event %s — skipping its requests", dk_event_id)
                event_state_cache[dk_event_id] = None
                continue
            event_state_cache[dk_event_id] = state
            team_names_cache[dk_event_id] = extract_team_names_from_event(state)

        if event_state_cache[dk_event_id] is None:
            continue

        team_names = team_names_cache[dk_event_id]
        sel_ids = resolve_legs(legs, side, event_state_cache[dk_event_id], team_names)
        if sel_ids is None:
            log.info("Could not resolve all legs for game=%s prop=%s side=%s — writing NULL row",
                     game_id, prop_type, side)
            rows.append({
                'fetch_time': fetch_time, 'game_id': game_id,
                'prop_type': prop_type, 'side': side,
                'legs_json': json.dumps(legs), 'selection_ids': '',
                'sgp_decimal': None, 'sgp_american': None,
                'source': 'draftkings_direct',
            })
            continue

        sgp_decimal = post_calculate_bets(session, sel_ids)
        sgp_american = decimal_to_american(sgp_decimal) if sgp_decimal else None
        rows.append({
            'fetch_time': fetch_time, 'game_id': game_id,
            'prop_type': prop_type, 'side': side,
            'legs_json': json.dumps(legs), 'selection_ids': ','.join(sel_ids),
            'sgp_decimal': sgp_decimal, 'sgp_american': sgp_american,
            'source': 'draftkings_direct',
        })
        log.info("game=%s prop=%s side=%s sgp_decimal=%s",
                 game_id[:8], prop_type, side, sgp_decimal)

    log.info("Writing %d rows to %s", len(rows), db_path)
    con = duckdb.connect(str(db_path))
    try:
        write_rows(con, rows)
    finally:
        con.close()
    log.info("Done. fetch_time=%s", fetch_time.isoformat())
    return 0


if __name__ == '__main__':
    sys.exit(main())
```

### - [ ] Step 2: Syntax + import smoke

```bash
cd /Users/callancapitolo/NFLWork/.worktrees/dk-trifecta-scraper
python3 -c "
import ast
ast.parse(open('mlb_sgp/scraper_draftkings_trifecta.py').read())
print('SYNTAX OK')
"
/Users/callancapitolo/NFLWork/mlb_sgp/venv/bin/python -c "
import sys
sys.path.insert(0, 'mlb_sgp')
import scraper_draftkings_trifecta as st
assert callable(st.main)
assert callable(st.post_calculate_bets)
assert callable(st.write_rows)
print('IMPORT OK')
"
```
Expected: both prints succeed.

### - [ ] Step 3: Live smoke test against tonight's slate

```bash
cd /Users/callancapitolo/NFLWork/.worktrees/dk-trifecta-scraper
# Build a 1-row test input from one of tonight's specials
/Users/callancapitolo/NFLWork/mlb_sgp/venv/bin/python -c "
import json, duckdb
con = duckdb.connect('/Users/callancapitolo/NFLWork/wagerzon_odds/wagerzon.duckdb', read_only=True)
row = con.execute('''
  SELECT team, prop_type, description
  FROM wagerzon_specials
  WHERE prop_type = 'TRIPLE-PLAY'
    AND scraped_at = (SELECT MAX(scraped_at) FROM wagerzon_specials)
  LIMIT 1
''').fetchone()
con.close()
# Map to a request — just hand-craft a test for the team name's home game
import sys; sys.path.insert(0, '/Users/callancapitolo/NFLWork/Answer Keys')
# Use a minimal hardcoded test request (avoid pulling consensus_temp here)
test = [{
  'game_id': 'test-game-1',
  'home_team': 'Chicago Cubs',
  'away_team': 'San Diego Padres',
  'prop_type': 'TRIPLE-PLAY',
  'side': 'home',
  'legs': [
    {'type': 'scores_first'},
    {'type': 'wins_period', 'period': 'F5'},
    {'type': 'wins_period', 'period': 'FG'},
  ],
}]
open('/tmp/trifecta_smoke.json', 'w').write(json.dumps(test, indent=2))
print('Wrote /tmp/trifecta_smoke.json with 1 test request')
"

# Run the scraper against a temp DB so we don't pollute main's mlb.duckdb
cp "/Users/callancapitolo/NFLWork/Answer Keys/mlb.duckdb" /tmp/trifecta_smoke.duckdb
/Users/callancapitolo/NFLWork/mlb_sgp/venv/bin/python mlb_sgp/scraper_draftkings_trifecta.py \
  --input /tmp/trifecta_smoke.json \
  --db    /tmp/trifecta_smoke.duckdb
```

Expected: scraper logs `game=test-gam prop=TRIPLE-PLAY side=home sgp_decimal=<positive float>` and exits 0. `sgp_decimal` should be in the +200 to +1000 American range for typical triple-play SGPs (decimal 3.0-11.0).

### - [ ] Step 4: Verify the row was written

```bash
/Users/callancapitolo/NFLWork/mlb_sgp/venv/bin/python -c "
import duckdb
con = duckdb.connect('/tmp/trifecta_smoke.duckdb', read_only=True)
df = con.execute('SELECT * FROM mlb_trifecta_sgp_odds ORDER BY fetch_time DESC LIMIT 5').fetchdf()
print(df.to_string())
con.close()
"
rm /tmp/trifecta_smoke.json /tmp/trifecta_smoke.duckdb
```
Expected: 1 row printed with non-null `sgp_decimal`, valid `selection_ids` (3 comma-separated IDs), `source='draftkings_direct'`.

### - [ ] Step 5: Commit

```bash
cd /Users/callancapitolo/NFLWork/.worktrees/dk-trifecta-scraper
git add mlb_sgp/scraper_draftkings_trifecta.py
git commit -m "feat(dk): add scraper_draftkings_trifecta CLI + DuckDB writer

CLI: --input <json> --db <path>
- Reads request list, resolves legs via dk_leg_resolvers.LEG_RESOLVERS
- Posts to DK's calculateBets endpoint, captures SGP decimal
- Writes per-request rows to mlb_trifecta_sgp_odds with shared fetch_time
  (atomic snapshot semantics — pricer's MAX(fetch_time) gets all rows or none)
- NULL sgp_decimal when leg resolution or calculateBets fails (R blend
  degrades to model-only)

Live smoke confirmed end-to-end on tonight's CHI Cubs game."
```

---

## Task 4: Resolver unit tests

**Status: Depends on Task 2.**

**Files:**
- Create: `mlb_sgp/tests/test_dk_leg_resolvers.py`
- Create: `mlb_sgp/tests/__init__.py` (if directory doesn't exist yet)

### - [ ] Step 1: Create test directory + tests

```bash
cd /Users/callancapitolo/NFLWork/.worktrees/dk-trifecta-scraper
mkdir -p mlb_sgp/tests
touch mlb_sgp/tests/__init__.py
```

Create `mlb_sgp/tests/test_dk_leg_resolvers.py`:

```python
"""
Unit tests for dk_leg_resolvers.

Uses a tiny inline fixture (not a captured recon JSON) to keep tests fast
and avoid committing large blobs. Fixture mimics the relevant DK SGP event
payload structure with one selection per market.
"""
import sys
from pathlib import Path

# Add mlb_sgp to path so we can import resolvers
sys.path.insert(0, str(Path(__file__).parent.parent))

from dk_leg_resolvers import (  # noqa: E402
    resolve_scores_first,
    resolve_wins_period,
    resolve_team_total_under,
    resolve_team_total_over,
    resolve_legs,
    find_market_by_name,
)

# ---------------------------------------------------------------------------
# Inline fixture — minimal DK event payload covering the 4 supported markets
# ---------------------------------------------------------------------------

FIXTURE = {
    'data': {
        'markets': [
            {
                'marketName': '1st Run',
                'tags': ['SGP', 'YourBetEligible'],
                'selections': [
                    {'name': 'CHI Cubs',  'id': '0QA-CUBS-1ST-RUN'},
                    {'name': 'SD Padres', 'id': '0QA-SD-1ST-RUN'},
                ],
            },
            {
                'marketName': 'Moneyline',
                'tags': ['SGP'],
                'selections': [
                    {'name': 'CHI Cubs',  'id': '0ML-CUBS_3'},
                    {'name': 'SD Padres', 'id': '0ML-SD_1'},
                ],
            },
            {
                'marketName': '1st 5 Innings',
                'tags': ['SGP'],
                'selections': [
                    {'name': 'CHI Cubs',  'id': '0QA-CUBS-F5'},
                    {'name': 'SD Padres', 'id': '0QA-SD-F5'},
                ],
            },
            {
                'marketName': 'CHI Cubs: Team Total Runs',
                'tags': ['SGP'],
                'selections': [
                    {'label': 'Over',  'point': 3.5, 'id': '0OU-CUBS-O35'},
                    {'label': 'Under', 'point': 3.5, 'id': '0OU-CUBS-U35'},
                ],
            },
            {
                'marketName': 'SD Padres: Team Total Runs',
                'tags': ['SGP'],
                'selections': [
                    {'label': 'Over',  'point': 4.5, 'id': '0OU-SD-O45'},
                    {'label': 'Under', 'point': 4.5, 'id': '0OU-SD-U45'},
                ],
            },
            # Non-SGP-eligible market — should be ignored by find_market_by_name
            {
                'marketName': 'Inning of First/Last Score',
                'tags': ['Cashout'],   # NO 'SGP' tag
                'selections': [
                    {'name': 'CHI Cubs',  'id': '0QA-IGNORE-1'},
                    {'name': 'SD Padres', 'id': '0QA-IGNORE-2'},
                ],
            },
        ],
    },
}

TEAM_NAMES = {'home': 'CHI Cubs', 'away': 'SD Padres'}


# ---------------------------------------------------------------------------
# Tests
# ---------------------------------------------------------------------------

def test_find_market_by_name_finds_sgp_eligible_market():
    m = find_market_by_name(FIXTURE, '1st Run')
    assert m is not None
    assert m['marketName'] == '1st Run'

def test_find_market_by_name_skips_non_sgp_market():
    # 'Inning of First/Last Score' is in fixture but lacks 'SGP' tag
    m = find_market_by_name(FIXTURE, 'Inning of First/Last Score')
    assert m is None

def test_find_market_by_name_returns_none_for_missing_market():
    assert find_market_by_name(FIXTURE, 'No Such Market') is None


def test_scores_first_home():
    sid = resolve_scores_first({'type': 'scores_first'}, 'home', FIXTURE, TEAM_NAMES)
    assert sid == '0QA-CUBS-1ST-RUN'

def test_scores_first_away():
    sid = resolve_scores_first({'type': 'scores_first'}, 'away', FIXTURE, TEAM_NAMES)
    assert sid == '0QA-SD-1ST-RUN'


def test_wins_period_FG_home():
    sid = resolve_wins_period({'type': 'wins_period', 'period': 'FG'}, 'home', FIXTURE, TEAM_NAMES)
    assert sid == '0ML-CUBS_3'

def test_wins_period_F5_away():
    sid = resolve_wins_period({'type': 'wins_period', 'period': 'F5'}, 'away', FIXTURE, TEAM_NAMES)
    assert sid == '0QA-SD-F5'

def test_wins_period_F3_returns_none():
    # F3 not in DK's 2-way ML primitives — resolver returns None
    sid = resolve_wins_period({'type': 'wins_period', 'period': 'F3'}, 'home', FIXTURE, TEAM_NAMES)
    assert sid is None


def test_team_total_under_home_correct_line():
    sid = resolve_team_total_under(
        {'type': 'team_total_under', 'line': 3.5}, 'home', FIXTURE, TEAM_NAMES,
    )
    assert sid == '0OU-CUBS-U35'

def test_team_total_over_away_correct_line():
    sid = resolve_team_total_over(
        {'type': 'team_total_over', 'line': 4.5}, 'away', FIXTURE, TEAM_NAMES,
    )
    assert sid == '0OU-SD-O45'

def test_team_total_wrong_line_returns_none():
    sid = resolve_team_total_under(
        {'type': 'team_total_under', 'line': 99.5}, 'home', FIXTURE, TEAM_NAMES,
    )
    assert sid is None


def test_resolve_legs_full_trifecta_home():
    legs = [
        {'type': 'scores_first'},
        {'type': 'wins_period', 'period': 'F5'},
        {'type': 'wins_period', 'period': 'FG'},
    ]
    sids = resolve_legs(legs, 'home', FIXTURE, TEAM_NAMES)
    assert sids == ['0QA-CUBS-1ST-RUN', '0QA-CUBS-F5', '0ML-CUBS_3']

def test_resolve_legs_full_grand_slam_away():
    legs = [
        {'type': 'scores_first'},
        {'type': 'wins_period', 'period': 'F5'},
        {'type': 'wins_period', 'period': 'FG'},
        {'type': 'team_total_under', 'line': 4.5},
    ]
    sids = resolve_legs(legs, 'away', FIXTURE, TEAM_NAMES)
    assert sids == ['0QA-SD-1ST-RUN', '0QA-SD-F5', '0ML-SD_1', '0OU-SD-U45']

def test_resolve_legs_missing_leg_returns_none():
    # F3 not supported → whole resolution fails → None
    legs = [
        {'type': 'scores_first'},
        {'type': 'wins_period', 'period': 'F3'},
        {'type': 'wins_period', 'period': 'FG'},
    ]
    sids = resolve_legs(legs, 'home', FIXTURE, TEAM_NAMES)
    assert sids is None

def test_resolve_legs_unknown_leg_type_returns_none():
    legs = [{'type': 'scores_first'}, {'type': 'mystery_leg'}]
    sids = resolve_legs(legs, 'home', FIXTURE, TEAM_NAMES)
    assert sids is None
```

### - [ ] Step 2: Run tests, confirm they fail (none implemented yet — but Task 2 already implemented them so they should PASS)

```bash
cd /Users/callancapitolo/NFLWork/.worktrees/dk-trifecta-scraper
/Users/callancapitolo/NFLWork/mlb_sgp/venv/bin/pip install pytest >/dev/null 2>&1 || true
/Users/callancapitolo/NFLWork/mlb_sgp/venv/bin/python -m pytest mlb_sgp/tests/test_dk_leg_resolvers.py -v
```
Expected: all 14 tests pass. (TDD inversion: the implementation already exists from Task 2; these tests are the safety net. If any test fails, that's a real bug in Task 2's resolvers.)

### - [ ] Step 3: Commit

```bash
cd /Users/callancapitolo/NFLWork/.worktrees/dk-trifecta-scraper
git add mlb_sgp/tests/__init__.py mlb_sgp/tests/test_dk_leg_resolvers.py
git commit -m "test(dk): resolver unit tests against captured recon fixture

14 tests covering all 4 supported leg types (scores_first, wins_period
FG/F5/F3-rejected, team_total_under/over, line-mismatch rejection).
Uses an inline minimal fixture (no large recon JSON committed)."
```

---

## Task 5: R-side activation

**Status: Depends on Task 3 (scraper exists).**

**Files:**
- Modify: `Answer Keys/mlb_triple_play.R` — uncomment Plan #2 stub + add real `system2()` invocation

### - [ ] Step 1: Find the Plan #2 stub block

In `Answer Keys/mlb_triple_play.R`, find the comment block beginning with `# ===== Plan #2 will activate the DK trifecta SGP scraper here =====`. It currently looks like:

```r
  # ===== Plan #2 will activate the DK trifecta SGP scraper here =====
  # When Plan #2 lands, replace the comment block below with a writable
  # request file + system2() call to mlb_sgp/scraper_draftkings_trifecta.py.
  # The blend logic below already handles the populated table correctly;
  # no changes to the rowwise mutate are needed when activating.
  #
  # trifecta_input <- todays_lines %>% ... write JSON ...
  # system2(SGP_VENV_PYTHON,
  #         args = c(file.path(SGP_DIR, "scraper_draftkings_trifecta.py"),
  #                  "--input", request_path, "--db", MLB_DB))
  # ===================================================================
```

### - [ ] Step 2: Replace the stub with the real invocation

Replace the entire stub block above with:

```r
  # ===== DK trifecta SGP scraper invocation (Plan #2) =====
  # Build a JSON request file with one entry per matched line, then call the
  # Python scraper via system2(). Scraper writes rows to mlb_trifecta_sgp_odds;
  # we read them back below in the dk_sgp <- dbGetQuery(...) block.
  #
  # If the scraper exits non-zero, we log a warning and continue — the blend
  # gracefully degrades to model-only when the table has no fresh rows.
  SGP_VENV_PYTHON <- "/Users/callancapitolo/NFLWork/mlb_sgp/venv/bin/python"
  SGP_DIR <- "/Users/callancapitolo/NFLWork/mlb_sgp"

  trifecta_input <- matched %>%
    rowwise() %>%
    mutate(legs_json = jsonlite::toJSON(parse_legs(description), auto_unbox = TRUE)) %>%
    ungroup() %>%
    transmute(
      game_id, home_team, away_team, prop_type,
      side, legs = legs_json
    )
  request_path <- tempfile(fileext = ".json")
  jsonlite::write_json(trifecta_input, request_path, auto_unbox = TRUE)

  dk_status <- tryCatch({
    system2(SGP_VENV_PYTHON,
            args = c(file.path(SGP_DIR, "scraper_draftkings_trifecta.py"),
                     "--input", request_path, "--db", MLB_DB),
            stdout = TRUE, stderr = TRUE)
    0
  }, warning = function(w) {
    cat(sprintf("DK trifecta scraper warning: %s\n", conditionMessage(w)))
    1
  }, error = function(e) {
    cat(sprintf("DK trifecta scraper failed: %s\n", conditionMessage(e)))
    1
  })
  if (dk_status != 0) {
    cat("DK trifecta blend disabled for this run; using model-only fair odds.\n")
  }
  unlink(request_path)
  # ========================================================
```

### - [ ] Step 3: Confirm `library(jsonlite)` is loaded

Look at the `suppressPackageStartupMessages({...})` block at the top of `mlb_triple_play.R`. If `library(jsonlite)` is not present, add it (next to `library(DBI)`).

```r
suppressPackageStartupMessages({
  library(duckdb)
  library(dplyr)
  library(tibble)
  library(DBI)
  library(jsonlite)   # ← ADD THIS LINE if missing
})
```

### - [ ] Step 4: Confirm `prop_type` is in `matched` before the scraper invocation

The post-Wagerzon pricer derives `prop_type` inside the rowwise mutate (after the system2 call). The scraper invocation needs `prop_type` BEFORE that mutate. Add a one-line derivation right above the trifecta_input build:

```r
  matched <- matched %>%
    mutate(
      prop_type = vapply(description, function(d) {
        m <- regmatches(d, regexec("\\b(TRIPLE-PLAY|GRAND-SLAM)\\b", d))[[1]]
        if (length(m) >= 2) m[[2]] else NA_character_
      }, character(1))
    )
```

(This is the same regex the rowwise mutate already uses — moving it earlier so the request file gets it.)

Also remove the duplicate `prop_type` derivation from the rowwise mutate further down — it's now redundant.

### - [ ] Step 5: Syntax check

```bash
cd /Users/callancapitolo/NFLWork/.worktrees/dk-trifecta-scraper
Rscript -e 'parse(file = "Answer Keys/mlb_triple_play.R"); cat("SYNTAX OK\n")'
```
Expected: `SYNTAX OK`.

### - [ ] Step 6: Live integration test — copy mlb.duckdb in, run pricer, observe DK populated

```bash
cd /Users/callancapitolo/NFLWork/.worktrees/dk-trifecta-scraper
cp "/Users/callancapitolo/NFLWork/Answer Keys/mlb.duckdb" "Answer Keys/mlb.duckdb"
cd "Answer Keys" && Rscript mlb_triple_play.R 2>&1 | tail -25
```

Expected outputs:
1. `DK trifecta scraper invocation` log line followed by per-game scraper output
2. The 14-row table prints with `dk_odds` populated for at least 8/14 rows (the rest may legitimately be NA if DK can't price specific game/leg combinations)
3. `fair_odds` differs from `model_odds` for rows where `dk_odds` is populated — this is the blend in action
4. Edge column reflects the blended fair, not the model-only fair

If `dk_odds` is NA for ALL rows, something is wrong — either:
- The scraper is silently failing (check stderr in the log)
- DK is rejecting all combinations (combinability test result was misleading)
- The team-name mapping isn't resolving Odds API → DK event id

Stop and investigate before committing.

### - [ ] Step 7: Remove the copied mlb.duckdb

```bash
cd /Users/callancapitolo/NFLWork/.worktrees/dk-trifecta-scraper
rm "Answer Keys/mlb.duckdb"
```

### - [ ] Step 8: Commit

```bash
git add "Answer Keys/mlb_triple_play.R"
git commit -m "feat(mlb): activate DK trifecta scraper invocation in pricer

Replaces Plan #1's stub with a real system2() call that:
- Writes the matched-lines JSON to a temp file
- Invokes mlb_sgp/scraper_draftkings_trifecta.py
- Cleans up the temp file
- Gracefully handles scraper failure (model-only fallback)

Plan #1's blend logic in the rowwise mutate is unchanged — it already
correctly handles the populated mlb_trifecta_sgp_odds table.

Live verified: 8+/14 rows now show populated dk_odds and blended fair_odds."
```

---

## Task 6: Documentation update

**Status: Depends on Task 5.**

**Files:**
- Modify: `Answer Keys/CLAUDE.md` — replace the "Plan #2 (post-recon) populates ..." hedge with live data flow

### - [ ] Step 1: Update the "Triple-Play Data Flow" subsection

In `Answer Keys/CLAUDE.md`, find the bullet:

```
- Plan #2 (post-recon) populates `mlb_trifecta_sgp_odds` via `mlb_sgp/scraper_draftkings_trifecta.py`. Plan #1 ships the schema + blend scaffolding; the table stays empty until Plan #2 lands.
```

Replace with:

```
- DK trifecta SGP odds populated by `mlb_sgp/scraper_draftkings_trifecta.py`, invoked from the pricer via `system2()`. Resolvers in `mlb_sgp/dk_leg_resolvers.py` map each `parse_legs` leg type to a DK selection ID using primitive markets: `1st Run` (scores_first), `Moneyline` (FG ML), `1st 5 Innings` (F5 ML), `<TEAM>: Team Total Runs` (team totals). All four are SGP-eligible per recon (2026-04-22).
- F3 / F7 wins_period legs and opp_total_* legs return NULL from the resolver — DK doesn't reliably post these as 2-way primitives. Props containing them get NULL `sgp_decimal` and the R blend correctly degrades to model-only.
```

Also update the data-flow diagram block. Find:

```
  │     dk_sgp lookup → devig with DK_SGP_VIG_DEFAULT=1.25 → dk_fair_prob
```

Below it (or wherever appropriate in the flow), add a line documenting the scraper invocation:

```
  │     scraper_draftkings_trifecta.py (Python, invoked via system2()):
  │       Odds API id → DK event id → resolve_legs() → calculateBets POST
  │       writes mlb_trifecta_sgp_odds row per (game_id, prop_type, side)
```

### - [ ] Step 2: Verify the file reads correctly

```bash
grep -A 50 "Triple-Play Data Flow" "Answer Keys/CLAUDE.md" | head -65
```
Expected: the updated section showing the live data flow including `dk_leg_resolvers` and `scraper_draftkings_trifecta` references.

### - [ ] Step 3: Commit

```bash
cd /Users/callancapitolo/NFLWork/.worktrees/dk-trifecta-scraper
git add "Answer Keys/CLAUDE.md"
git commit -m "docs(mlb): document Plan #2 live data flow

Replaces the 'Plan #2 (post-recon) populates ...' hedge with live
documentation of dk_leg_resolvers + scraper_draftkings_trifecta. Notes
which leg types are supported (scores_first, wins_period FG/F5,
team_total_under/over) and which return NULL (F3/F7, opp_total_*)."
```

---

## Self-Review

**Spec coverage** (from `docs/superpowers/specs/2026-04-22-dk-trifecta-blend-design.md`):

| Spec section | Covered by |
|---|---|
| Path A — primitive 3/4-leg SGP construction | Task 1 verifies, Task 2 implements |
| LEG_RESOLVERS dict (parser-driven) | Task 2 |
| Production scraper with `--input <json> --db <path>` CLI | Task 3 |
| Resolver unit tests | Task 4 |
| R-side activation of system2() call | Task 5 |
| Live integration smoke test | Task 5 Step 6 |
| Documentation update | Task 6 |
| Open Item #5 (Odds API id → DK event id mapping) | Task 3's `map_odds_id_to_dk_event` |

**Placeholder scan:** No "TBD" / "implement later". Each step has either complete code or concrete shell with expected output.

**Type consistency:**
- `resolve_legs(legs, side, event_state, team_names)` — same signature in resolver module (Task 2), called identically in scraper (Task 3) and tests (Task 4)
- `LEG_RESOLVERS` dict keys match `parse_legs.R::TOKEN_REGISTRY` leg types (`scores_first`, `wins_period`, `team_total_under`, `team_total_over`)
- DuckDB column names (`fetch_time`, `game_id`, `prop_type`, `side`, `legs_json`, `selection_ids`, `sgp_decimal`, `sgp_american`, `source`) consistent across schema (Plan #1), scraper writer (Task 3), and pricer reader (already in main from Plan #1)
- `game_id` is VARCHAR everywhere (Odds API hex)

**Risk assessment:** Path A combinability is the single largest risk; Task 1 explicitly tests it before any other work. If it fails, the plan halts cleanly and we redesign around Path B without wasted code.

No issues found.

---

## Execution Handoff

**Plan complete and saved to `docs/superpowers/plans/2026-04-22-dk-trifecta-scraper.md`. Two execution options:**

**1. Subagent-Driven (recommended)** — fresh subagent per task, review between tasks, fast iteration

**2. Inline Execution** — execute tasks in this session with checkpoints

**Which approach?**
