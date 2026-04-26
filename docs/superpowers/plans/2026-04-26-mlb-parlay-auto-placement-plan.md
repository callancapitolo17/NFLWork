# MLB Parlay Auto-Placement Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Add one-click MLB correlated-parlay placement to the MLB answer key dashboard via Wagerzon's REST API, with zero-tolerance price-drift abort, in-process placement (no browser), and a dry-run validation toggle.

**Architecture:** New Python module `wagerzon_odds/parlay_placer.py` reuses the existing `parlay_pricer.py` session/login + leg encoding to call `ConfirmWagerHelper` (price preflight) then `PostWagerMultipleHelper` (placement). A new `/api/place-parlay` endpoint on `mlb_dashboard_server.py` wires the dashboard's Place button to the placer. State lives in `placed_parlays` and `placement_orphans` tables in `mlb_dashboard.duckdb`.

**Tech Stack:** Python 3 (`requests`), R/Shiny (existing dashboard), DuckDB.

**Spec reference:** `docs/superpowers/specs/2026-04-26-mlb-parlay-auto-placement-design.md`

---

## Pre-flight checklist (read once, then proceed)

1. **Branch hygiene** — All code work happens on `feature/mlb-parlay-auto-place` in a worktree. Per repo CLAUDE.md, never write code on `main`.
2. **Never symlink DuckDB** — `cp` databases into the worktree if needed; better yet, test from `main` after merge.
3. **Recon JSON** — Already scrubbed (password replaced with `REDACTED`, gitignored at `.gitignore:38`). Don't re-record placement traffic without scrubbing.
4. **Real money risk** — Tasks 13 and 14 place real bets. Use the $15 minimum. Verify state against Wagerzon's bet history before declaring success.
5. **Pre-merge approval gate** — Per repo CLAUDE.md, never merge to `main` without explicit user approval. Task 17 is the gate.

---

## Task 1: Set up worktree, branch, and environment

**Files:** none modified — repo-level setup only.

- [ ] **Step 1: Create worktree and feature branch**

  Run from main repo (`/Users/callancapitolo/NFLWork`):
  ```bash
  git worktree add -b feature/mlb-parlay-auto-place ../NFLWork-mlb-parlay-place main
  cd ../NFLWork-mlb-parlay-place
  git branch --show-current
  ```
  Expected: prints `feature/mlb-parlay-auto-place`.

- [ ] **Step 2: Verify wagerzon_odds venv works inside the worktree**

  ```bash
  source wagerzon_odds/venv/bin/activate
  python3 -c "import requests, playwright; print('OK')"
  deactivate
  ```
  Expected: `OK`. If venv is missing in the worktree (it's a directory, not a symlink, so it should be present), recreate per `wagerzon_odds/README.md`.

- [ ] **Step 3: Confirm `.env` is accessible (read-only)**

  ```bash
  test -f bet_logger/.env && echo "env present" || echo "MISSING"
  ```
  Expected: `env present`. The `.env` is shared across scrapers, placer, etc. It must contain `WAGERZON_USERNAME` and `WAGERZON_PASSWORD`.

- [ ] **Step 4: Make a probe commit on the new branch**

  Don't create a placeholder file. Just verify the branch is live:
  ```bash
  git log --oneline -3
  ```
  Expected: most recent commit is the spec on main, history shared.

---

## Task 2: Resolve the `idgm` open question

The placer needs Wagerzon's internal game ID (`idgm`) to build `sel`/`detailData` payloads. Spec §10 marks this as the first thing to resolve.

**Files:**
- Read: `wagerzon_odds/parlay_pricer.py` (no changes — discovery only)
- Read: `Answer Keys/mlb_correlated_parlay.R` (no changes — discovery only)
- Possible modify: `Answer Keys/mlb_correlated_parlay.R` (only if `idgm` needs to be added to `mlb_parlay_opportunities`)

- [ ] **Step 1: Read the pricer to find where `idgm` comes from**

  Open `wagerzon_odds/parlay_pricer.py`. Search for: `idgm`, `IdGame`, `NewScheduleHelper`. Identify the function/variable that produces `idgm` for each parlay leg before calling `ConfirmWagerHelper`.

- [ ] **Step 2: Trace how the R parlay generator uses (or doesn't use) `idgm`**

  Open `Answer Keys/mlb_correlated_parlay.R`. Search for: `idgm`, `IdGame`. Determine whether `idgm` is:
  - (a) Already a column in `mlb_parlay_opportunities` (read it and you're done)
  - (b) Computed at pricing time but not persisted (need to persist it)
  - (c) Looked up via a separate join (placer must do the same join)

- [ ] **Step 3: Verify by querying the DB**

  ```bash
  python3 -c "
  import duckdb
  con = duckdb.connect('Answer Keys/mlb.duckdb', read_only=True)
  cols = [r[1] for r in con.execute('PRAGMA table_info(mlb_parlay_opportunities)').fetchall()]
  print('idgm-bearing cols:', [c for c in cols if 'idgm' in c.lower() or 'gameid' in c.lower().replace('_','')])
  print('all cols:', cols)
  con.close()
  "
  ```
  Expected: prints column list. If a column with `idgm`/`IdGame`/`game_id_wz` exists → case (a). Otherwise → case (b) or (c).

- [ ] **Step 4: Document the finding in this plan**

  Edit this file: under this task, add a "Finding" subsection with one paragraph describing where `idgm` comes from and which option (a/b/c) we landed on. This becomes the source of truth for Task 4 onward.

  ```bash
  $EDITOR docs/superpowers/plans/2026-04-26-mlb-parlay-auto-placement-plan.md
  # add finding paragraph
  ```

- [ ] **Step 5: If case (b) — persist `idgm` to `mlb_parlay_opportunities`**

  Modify `Answer Keys/mlb_correlated_parlay.R` to write `idgm_home` and `idgm_away` (or just `idgm` if there's one per parlay row) to the table. Run the pricer once to refresh the table. Spec §2 said "no schema changes" for this PR — this exception is justified because the schema is missing data the placer cannot synthesize.

  If case (a) or (c), skip this step.

- [ ] **Step 6: Commit the finding (and any pricer change)**

  ```bash
  git add docs/superpowers/plans/2026-04-26-mlb-parlay-auto-placement-plan.md
  # if step 5 ran:
  git add "Answer Keys/mlb_correlated_parlay.R"
  git commit -m "docs(plans): document idgm sourcing for parlay placer"
  ```

## Finding

**Case (a) — `idgm` is already a column in `mlb_parlay_opportunities`.** The Wagerzon scraper populates `mlb_odds.idgm` (sourced from the `IdGame` field in Wagerzon's schedule/odds response), and `parlay_pricer.py::_price_mlb_parlays_inner` reads it via `SELECT DISTINCT home_team, away_team, idgm, ...` from `mlb_odds`. The R generator (`mlb_correlated_parlay.R`) then carries `idgm` forward from the `wz` (Wagerzon odds) data frame into the per-combo result tibble at line 744: `idgm = if ("idgm" %in% names(row)) as.integer(row$idgm) else NA_integer_`. Since all four combos for a game share the same game, a single `idgm` integer per parlay row is sufficient — no per-leg fan-out is needed. The column is confirmed populated with real integer values in the live DB (e.g., `5628895`, `5631648`). **No R or Python changes are required.** The placer (Task 9 `_build_spec_from_row`) can read `idgm` directly from the DB row as `int(row["idgm"])` — the `TODO_T2` placeholder in the plan's Task 9 pseudocode is already correct as written.

---

## Task 3: Database schema migration

**Files:**
- Create: `wagerzon_odds/migrate_placed_parlays.py` (one-shot, will be run once and committed)

- [ ] **Step 1: Write the migration script**

  ```python
  # wagerzon_odds/migrate_placed_parlays.py
  """One-shot migration: ensures placed_parlays has all columns required by
  parlay_placer + creates placement_orphans. Idempotent (safe to re-run)."""
  import os
  import duckdb

  DB_PATH = os.path.join(
      os.path.dirname(__file__), "..", "Answer Keys", "MLB Dashboard",
      "mlb_dashboard.duckdb"
  )

  REQUIRED_COLS = {
      "parlay_hash":      "VARCHAR PRIMARY KEY",
      "status":           "VARCHAR NOT NULL",
      "combo":            "VARCHAR",
      "game_id":          "VARCHAR",
      "game_time":        "TIMESTAMP",
      "recommended_size": "DOUBLE",
      "expected_odds":    "INTEGER",
      "expected_win":     "DOUBLE",
      "actual_size":      "DOUBLE",
      "actual_win":       "DOUBLE",
      "ticket_number":    "VARCHAR",
      "idwt":             "BIGINT",
      "legs_json":        "VARCHAR",
      "error_msg":        "VARCHAR",
      "error_msg_key":    "VARCHAR",
      "placed_at":        "TIMESTAMP DEFAULT CURRENT_TIMESTAMP",
      "updated_at":       "TIMESTAMP DEFAULT CURRENT_TIMESTAMP",
  }

  def main():
      con = duckdb.connect(DB_PATH)
      try:
          # Create placed_parlays if missing
          con.execute(f"""
              CREATE TABLE IF NOT EXISTS placed_parlays (
                  parlay_hash VARCHAR PRIMARY KEY,
                  status      VARCHAR NOT NULL,
                  placed_at   TIMESTAMP DEFAULT CURRENT_TIMESTAMP
              )
          """)

          # Add any missing columns
          existing = {r[1] for r in con.execute(
              "PRAGMA table_info(placed_parlays)"
          ).fetchall()}
          for col, ddl in REQUIRED_COLS.items():
              if col in existing:
                  continue
              # PK / NOT NULL constraints can't be added via ALTER. If they're
              # missing on an existing table, that's fine (the create-table
              # above guarantees PK; NOT NULL is enforced application-side).
              ddl_alter = ddl.replace("PRIMARY KEY", "").replace("NOT NULL", "").strip()
              con.execute(f"ALTER TABLE placed_parlays ADD COLUMN IF NOT EXISTS {col} {ddl_alter}")
              print(f"  added column placed_parlays.{col}")

          # Create placement_orphans
          con.execute("""
              CREATE TABLE IF NOT EXISTS placement_orphans (
                  idwt          BIGINT PRIMARY KEY,
                  ticket_number VARCHAR,
                  parlay_hash   VARCHAR,
                  raw_response  VARCHAR,
                  error         VARCHAR,
                  created_at    TIMESTAMP DEFAULT CURRENT_TIMESTAMP
              )
          """)
          print("placement_orphans ready")
      finally:
          con.close()

  if __name__ == "__main__":
      main()
  ```

- [ ] **Step 2: Verify the script is idempotent**

  Run twice in a row:
  ```bash
  source wagerzon_odds/venv/bin/activate
  python3 wagerzon_odds/migrate_placed_parlays.py
  python3 wagerzon_odds/migrate_placed_parlays.py
  deactivate
  ```
  Expected: first run prints any added columns; second run is silent (no errors, no duplicate column adds).

- [ ] **Step 3: Inspect resulting schema**

  ```bash
  python3 -c "
  import duckdb
  con = duckdb.connect('Answer Keys/MLB Dashboard/mlb_dashboard.duckdb', read_only=True)
  for t in ('placed_parlays', 'placement_orphans'):
      print(f'\\n=== {t} ===')
      for r in con.execute(f'PRAGMA table_info({t})').fetchall():
          print(' ', r[1], r[2])
  con.close()
  "
  ```
  Expected: all 17 columns of `placed_parlays`, all 6 columns of `placement_orphans`.

- [ ] **Step 4: Commit**

  ```bash
  git add wagerzon_odds/migrate_placed_parlays.py
  git commit -m "feat: add placed_parlays + placement_orphans schema migration"
  ```

---

## Task 4: Placer module — dataclasses + leg encoding helpers

**Files:**
- Create: `wagerzon_odds/parlay_placer.py`
- Create: `wagerzon_odds/test_parlay_placer.py`

- [ ] **Step 1: Write failing test for dataclasses + encoding**

  ```python
  # wagerzon_odds/test_parlay_placer.py
  """Unit tests for parlay_placer. All Wagerzon HTTP calls are mocked."""
  import pytest
  from parlay_placer import Leg, ParlaySpec, encode_sel, encode_detail_data


  def test_leg_dataclass_basic():
      leg = Leg(idgm=5632938, play=0, points=-1.5, odds=117, pitcher=0)
      assert leg.idgm == 5632938
      assert leg.play == 0


  def test_parlay_spec_dataclass_basic():
      spec = ParlaySpec(
          parlay_hash="abc123",
          legs=[
              Leg(idgm=5632938, play=1, points=-1.5, odds=117),
              Leg(idgm=5632938, play=2, points=-7.5, odds=105),
          ],
          amount=15.0,
          expected_win=12.50,
          expected_risk=15.0,
      )
      assert len(spec.legs) == 2


  def test_encode_sel_two_legs():
      legs = [
          Leg(idgm=5632938, play=1, points=-1.5, odds=117),
          Leg(idgm=5632938, play=2, points=-7.5, odds=105),
      ]
      assert encode_sel(legs) == "1_5632938_-1.5_117,2_5632938_-7.5_105"


  def test_encode_sel_negative_odds():
      legs = [Leg(idgm=5632938, play=5, points=0, odds=-140)]
      assert encode_sel(legs) == "5_5632938_0_-140"


  def test_encode_detail_data_shape():
      legs = [Leg(idgm=5632938, play=1, points=-1.5, odds=117)]
      out = encode_detail_data(legs, amount=15.0)
      # detail_data is a list of dicts ready for json.dumps
      assert out[0]["IdGame"] == 5632938
      assert out[0]["Play"] == 1
      assert out[0]["Amount"] == "15"
      assert out[0]["Points"]["selected"] is True
  ```

- [ ] **Step 2: Run test to verify failure**

  ```bash
  cd wagerzon_odds && source venv/bin/activate
  pytest test_parlay_placer.py -v
  ```
  Expected: `ModuleNotFoundError: No module named 'parlay_placer'` (or `ImportError`).

- [ ] **Step 3: Create the parlay_placer module skeleton**

  ```python
  # wagerzon_odds/parlay_placer.py
  """Wagerzon parlay placer — pure REST, no browser.

  See docs/superpowers/specs/2026-04-26-mlb-parlay-auto-placement-design.md
  for the full design. This module provides:
      - ParlaySpec / Leg / PlacementResult dataclasses
      - encode_sel / encode_detail_data leg-encoding helpers
      - place_parlays(specs, dry_run=False) — top-level entry point
  """
  from __future__ import annotations
  from dataclasses import dataclass, field
  from typing import Optional


  @dataclass
  class Leg:
      """One leg of a parlay.

      play codes: 0=away spread, 1=home spread, 2=over, 3=under,
                  4=away ML, 5=home ML
      points: signed line value (negative for fav spread / over)
      odds: integer American odds (e.g. 117 = +117, -130 = -130)
      pitcher: 0 = action / 3 = listed (mirrors what ConfirmWagerHelper expects)
      """
      idgm: int
      play: int
      points: float
      odds: int
      pitcher: int = 0


  @dataclass
  class ParlaySpec:
      parlay_hash: str
      legs: list[Leg]
      amount: float
      expected_win: float
      expected_risk: float


  @dataclass
  class PlacementResult:
      parlay_hash: str
      status: str
      ticket_number: Optional[str] = None
      idwt: Optional[int] = None
      actual_win: Optional[float] = None
      actual_risk: Optional[float] = None
      error_msg: str = ""
      error_msg_key: str = ""
      raw_response: Optional[str] = None


  def encode_sel(legs: list[Leg]) -> str:
      """Encode legs as Wagerzon's `sel` parameter: play_idgm_points_odds, ..."""
      parts = []
      for leg in legs:
          # points printed without trailing .0 if integer
          pts = int(leg.points) if leg.points == int(leg.points) else leg.points
          parts.append(f"{leg.play}_{leg.idgm}_{pts}_{leg.odds}")
      return ",".join(parts)


  def encode_detail_data(legs: list[Leg], amount: float) -> list[dict]:
      """Build the `detailData` JSON array for ConfirmWagerHelper / PostWager."""
      return [
          {
              "Amount": str(int(amount)) if amount == int(amount) else str(amount),
              "RiskWin": 0,
              "TeaserPointsPurchased": 0,
              "IdGame": leg.idgm,
              "Play": leg.play,
              "Pitcher": leg.pitcher,
              "Points": {
                  "BuyPoints": 0,
                  "BuyPointsDesc": "",
                  "LineDesc": "",
                  "selected": True,
              },
          }
          for leg in legs
      ]
  ```

- [ ] **Step 4: Run tests, verify pass**

  ```bash
  pytest test_parlay_placer.py -v
  ```
  Expected: all 5 tests pass.

- [ ] **Step 5: Commit**

  ```bash
  git add wagerzon_odds/parlay_placer.py wagerzon_odds/test_parlay_placer.py
  git commit -m "feat(parlay_placer): dataclasses + sel/detailData encoding"
  ```

---

## Task 5: Placer module — login and session management

**Files:**
- Modify: `wagerzon_odds/parlay_placer.py` (add `_get_session`)
- Modify: `wagerzon_odds/test_parlay_placer.py`

- [ ] **Step 1: Read the existing login pattern**

  Open `wagerzon_odds/parlay_pricer.py`. Locate the function that establishes a logged-in `requests.Session` (likely `_api_login` or similar — the same pattern is also in `bet_placer/navigator_wagerzon.py:646-669`). Note: it parses `__VIEWSTATE` etc. from an HTML response, then form-POSTs `Account` + `Password`.

- [ ] **Step 2: Write a failing test for `_get_session`**

  ```python
  # additions to test_parlay_placer.py
  from unittest.mock import patch, MagicMock
  import parlay_placer


  def test_get_session_caches_within_call(monkeypatch):
      """Two calls to _get_session in the same place_parlays() invocation
      should reuse one session, not log in twice."""
      monkeypatch.setenv("WAGERZON_USERNAME", "user")
      monkeypatch.setenv("WAGERZON_PASSWORD", "pass")

      with patch("parlay_placer.requests.Session") as MockSession:
          mock_session = MagicMock()
          # Simulate already-authenticated GET (URL contains NewSchedule)
          response = MagicMock(url="https://backend.wagerzon.com/wager/NewSchedule.aspx")
          mock_session.get.return_value = response
          MockSession.return_value = mock_session

          s1 = parlay_placer._get_session()
          s2 = parlay_placer._get_session()
          assert s1 is s2
          # Only one Session() instantiation
          assert MockSession.call_count == 1
  ```

- [ ] **Step 3: Run, verify fail**

  ```bash
  pytest test_parlay_placer.py::test_get_session_caches_within_call -v
  ```
  Expected: FAIL — `_get_session` not defined.

- [ ] **Step 4: Implement `_get_session` in parlay_placer.py**

  Add to `parlay_placer.py`:
  ```python
  import os
  import re
  import requests

  WAGERZON_BASE_URL = "https://backend.wagerzon.com"

  _CACHED_SESSION: Optional[requests.Session] = None


  def _get_session() -> requests.Session:
      """Return a logged-in Wagerzon session, reusing a cached one if present.

      Mirrors the login pattern in wagerzon_odds/parlay_pricer.py: GET base URL,
      if not already authenticated, parse __VIEWSTATE etc. and form-POST
      Account/Password.

      Resets via _clear_session() (called from auth_error retry path).
      """
      global _CACHED_SESSION
      if _CACHED_SESSION is not None:
          return _CACHED_SESSION

      session = requests.Session()
      resp = session.get(WAGERZON_BASE_URL, timeout=15)
      resp.raise_for_status()
      if "NewSchedule" in resp.url or "Welcome" in resp.url:
          _CACHED_SESSION = session
          return session

      html = resp.text
      fields = {}
      for name in ("__VIEWSTATE", "__VIEWSTATEGENERATOR", "__EVENTVALIDATION",
                   "__EVENTTARGET", "__EVENTARGUMENT"):
          m = re.search(rf'(?:name|id)="{name}"[^>]*value="([^"]*)"', html)
          if m:
              fields[name] = m.group(1)
      fields["Account"] = os.environ["WAGERZON_USERNAME"]
      fields["Password"] = os.environ["WAGERZON_PASSWORD"]
      fields["BtnSubmit"] = ""

      resp = session.post(WAGERZON_BASE_URL, data=fields, timeout=15)
      resp.raise_for_status()
      _CACHED_SESSION = session
      return session


  def _clear_session() -> None:
      """Force the next _get_session() call to re-login."""
      global _CACHED_SESSION
      _CACHED_SESSION = None
  ```

- [ ] **Step 5: Run, verify pass**

  ```bash
  pytest test_parlay_placer.py -v
  ```
  Expected: all tests pass.

- [ ] **Step 6: Add a fixture-reset helper to the test file**

  At top of `test_parlay_placer.py`, add:
  ```python
  @pytest.fixture(autouse=True)
  def reset_session():
      """Clear cached session before/after each test to avoid leakage."""
      import parlay_placer
      parlay_placer._clear_session()
      yield
      parlay_placer._clear_session()
  ```

- [ ] **Step 7: Re-run tests to confirm fixture works**

  ```bash
  pytest test_parlay_placer.py -v
  ```
  Expected: all pass.

- [ ] **Step 8: Commit**

  ```bash
  git add wagerzon_odds/parlay_placer.py wagerzon_odds/test_parlay_placer.py
  git commit -m "feat(parlay_placer): session login + caching with reset"
  ```

---

## Task 6: Placer module — ConfirmWagerHelper preflight + drift check

**Files:**
- Modify: `wagerzon_odds/parlay_placer.py`
- Modify: `wagerzon_odds/test_parlay_placer.py`

- [ ] **Step 1: Write failing tests for preflight**

  Add to `test_parlay_placer.py`:
  ```python
  import json


  CONFIRM_OK_RESPONSE = {
      "result": {
          "details": [{"Risk": 15.0, "Win": 30.0, "WagerType": 1,
                       "WagerTypeDesc": "PARLAY (2 TEAMS)"}],
          "Confirm": True,
      }
  }


  def _make_session_with_post(json_response):
      """Build a mocked requests.Session whose .post returns the given JSON."""
      sess = MagicMock()
      r = MagicMock()
      r.json.return_value = json_response
      r.headers = {"content-type": "application/json"}
      r.status_code = 200
      r.text = json.dumps(json_response)
      sess.post.return_value = r
      return sess


  def test_confirm_preflight_returns_win_risk(monkeypatch):
      monkeypatch.setattr(parlay_placer, "_get_session",
                          lambda: _make_session_with_post(CONFIRM_OK_RESPONSE))
      legs = [Leg(idgm=5632938, play=1, points=-1.5, odds=117)]
      win, risk = parlay_placer._confirm_preflight(legs, amount=15.0)
      assert win == 30.0
      assert risk == 15.0


  def test_drift_check_within_penny_passes():
      assert parlay_placer._drift_ok(expected=30.00, actual=30.005) is True
      assert parlay_placer._drift_ok(expected=30.00, actual=29.995) is True


  def test_drift_check_beyond_penny_fails():
      assert parlay_placer._drift_ok(expected=30.00, actual=29.50) is False
      assert parlay_placer._drift_ok(expected=30.00, actual=28.50) is False
  ```

- [ ] **Step 2: Run, verify fail**

  ```bash
  pytest test_parlay_placer.py -v -k "preflight or drift"
  ```
  Expected: FAIL — `_confirm_preflight` and `_drift_ok` not defined.

- [ ] **Step 3: Implement preflight + drift check**

  Add to `parlay_placer.py`:
  ```python
  import json as _json

  CONFIRM_URL = f"{WAGERZON_BASE_URL}/wager/ConfirmWagerHelper.aspx"
  DRIFT_TOLERANCE_USD = 0.01


  def _drift_ok(expected: float, actual: float) -> bool:
      """True if returned Win matches expected within $0.01."""
      return abs(actual - expected) <= DRIFT_TOLERANCE_USD


  def _confirm_preflight(legs: list[Leg], amount: float) -> tuple[float, float]:
      """Call ConfirmWagerHelper, return (Win, Risk).

      Raises ValueError on malformed response. Auth/HTML responses raise
      AuthExpired (handled by caller).
      """
      session = _get_session()
      detail_data = encode_detail_data(legs, amount)
      params = {
          "IDWT": "0",
          "WT": "1",
          "amountType": "0",
          "open": "0",
          "sameAmount": "false",
          "sameAmountNumber": str(int(amount)),
          "useFreePlayAmount": "false",
          "sel": encode_sel(legs),
          "detailData": _json.dumps(detail_data),
      }
      resp = session.post(CONFIRM_URL, data=params, timeout=15,
                          headers={"Accept": "application/json"})
      _raise_if_html(resp)
      data = resp.json()
      details = data["result"]["details"][0]
      return float(details["Win"]), float(details["Risk"])


  class AuthExpired(Exception):
      """Wagerzon returned an HTML login page instead of JSON."""


  def _raise_if_html(resp) -> None:
      ct = resp.headers.get("content-type", "")
      if "json" not in ct:
          raise AuthExpired(f"Non-JSON response: content-type={ct!r}")
  ```

- [ ] **Step 4: Run, verify pass**

  ```bash
  pytest test_parlay_placer.py -v
  ```
  Expected: all tests pass.

- [ ] **Step 5: Commit**

  ```bash
  git add wagerzon_odds/parlay_placer.py wagerzon_odds/test_parlay_placer.py
  git commit -m "feat(parlay_placer): ConfirmWagerHelper preflight + drift check"
  ```

---

## Task 7: Placer module — PostWagerMultipleHelper batch placement

**Files:**
- Modify: `wagerzon_odds/parlay_placer.py`
- Modify: `wagerzon_odds/test_parlay_placer.py`

- [ ] **Step 1: Write failing tests for `_post_wagers`**

  Add to `test_parlay_placer.py`:
  ```python
  POST_OK_RESPONSE = {
      "result": [{
          "WagerPostResult": {
              "details": [],
              "IDWT": 438173,
              "TicketNumber": "212147131",
              "Risk": 15.0,
              "Win": 30.0,
              "Confirm": True,
              "ErrorMsg": "",
              "ErrorMsgKey": "",
              "ErrorCode": {},
          }
      }]
  }

  POST_REJECTED_RESPONSE = {
      "result": [{
          "WagerPostResult": {
              "details": [],
              "Confirm": False,
              "ErrorMsg": "Insufficient funds",
              "ErrorMsgKey": "insufficient_funds",
              "ErrorCode": {"code": 42},
          }
      }]
  }


  def test_post_wagers_success(monkeypatch):
      monkeypatch.setattr(parlay_placer, "_get_session",
                          lambda: _make_session_with_post(POST_OK_RESPONSE))
      monkeypatch.setenv("WAGERZON_PASSWORD", "secret")
      specs = [ParlaySpec(
          parlay_hash="h1",
          legs=[Leg(idgm=5632938, play=1, points=-1.5, odds=117)],
          amount=15.0, expected_win=30.0, expected_risk=15.0,
      )]
      results = parlay_placer._post_wagers(specs)
      assert len(results) == 1
      r = results[0]
      assert r.status == "placed"
      assert r.ticket_number == "212147131"
      assert r.idwt == 438173
      assert r.actual_win == 30.0


  def test_post_wagers_rejected(monkeypatch):
      monkeypatch.setattr(parlay_placer, "_get_session",
                          lambda: _make_session_with_post(POST_REJECTED_RESPONSE))
      monkeypatch.setenv("WAGERZON_PASSWORD", "secret")
      specs = [ParlaySpec(
          parlay_hash="h1",
          legs=[Leg(idgm=5632938, play=1, points=-1.5, odds=117)],
          amount=15.0, expected_win=30.0, expected_risk=15.0,
      )]
      results = parlay_placer._post_wagers(specs)
      assert results[0].status == "rejected"
      assert results[0].error_msg_key == "insufficient_funds"
  ```

- [ ] **Step 2: Run, verify fail**

  ```bash
  pytest test_parlay_placer.py -v -k "post_wagers"
  ```
  Expected: FAIL — `_post_wagers` not defined.

- [ ] **Step 3: Implement `_post_wagers`**

  Add to `parlay_placer.py`:
  ```python
  POST_URL = f"{WAGERZON_BASE_URL}/wager/PostWagerMultipleHelper.aspx"

  ERROR_KEY_MAP = {
      "insufficient_funds": "insufficient balance",
      "bet_too_large":      "exceeds limit",
      "line_unavailable":   "line pulled",
  }


  def _build_post_request(spec: ParlaySpec, password: str) -> dict:
      detail_data = encode_detail_data(spec.legs, spec.amount)
      return {
          "WT": "1",
          "open": 0,
          "IDWT": "0",
          "sel": encode_sel(spec.legs),
          "sameAmount": False,
          "amountType": "0",
          "detailData": _json.dumps(detail_data),
          "confirmPassword": password,
          "sameAmountNumber": str(int(spec.amount)) if spec.amount == int(spec.amount) else str(spec.amount),
          "useFreePlayAmount": False,
          "roundRobinCombinations": "",
      }


  def _post_wagers(specs: list[ParlaySpec]) -> list[PlacementResult]:
      """POST one bulk request to PostWagerMultipleHelper. Returns one
      PlacementResult per input spec, in order.

      Raises AuthExpired on HTML response.
      """
      session = _get_session()
      password = os.environ["WAGERZON_PASSWORD"]
      payload = [_build_post_request(s, password) for s in specs]
      body = {"postWagerRequests": _json.dumps(payload)}

      resp = session.post(POST_URL, data=body, timeout=30,
                          headers={"Accept": "application/json"})
      _raise_if_html(resp)
      data = resp.json()
      raw = resp.text

      results = []
      for spec, item in zip(specs, data["result"]):
          wpr = item["WagerPostResult"]
          if wpr.get("Confirm") and not wpr.get("ErrorMsgKey"):
              results.append(PlacementResult(
                  parlay_hash=spec.parlay_hash,
                  status="placed",
                  ticket_number=str(wpr.get("TicketNumber")) if wpr.get("TicketNumber") else None,
                  idwt=int(wpr.get("IDWT")) if wpr.get("IDWT") else None,
                  actual_win=float(wpr.get("Win", 0)),
                  actual_risk=float(wpr.get("Risk", 0)),
                  raw_response=raw,
              ))
          else:
              key = wpr.get("ErrorMsgKey", "") or ""
              friendly = ERROR_KEY_MAP.get(key, key) or "unknown error"
              results.append(PlacementResult(
                  parlay_hash=spec.parlay_hash,
                  status="rejected",
                  error_msg=f"rejected: {friendly}" if key else (wpr.get("ErrorMsg") or "rejected"),
                  error_msg_key=key,
                  raw_response=raw,
              ))
      return results
  ```

- [ ] **Step 4: Run, verify pass**

  ```bash
  pytest test_parlay_placer.py -v
  ```
  Expected: all tests pass.

- [ ] **Step 5: Commit**

  ```bash
  git add wagerzon_odds/parlay_placer.py wagerzon_odds/test_parlay_placer.py
  git commit -m "feat(parlay_placer): PostWagerMultipleHelper batch placement"
  ```

---

## Task 8: Placer module — top-level `place_parlays`, retry, dry-run

**Files:**
- Modify: `wagerzon_odds/parlay_placer.py`
- Modify: `wagerzon_odds/test_parlay_placer.py`

- [ ] **Step 1: Write failing integration tests**

  Add to `test_parlay_placer.py`:
  ```python
  def _spec(hash="h", win=30.0, risk=15.0):
      return ParlaySpec(
          parlay_hash=hash,
          legs=[Leg(idgm=5632938, play=1, points=-1.5, odds=117)],
          amount=risk, expected_win=win, expected_risk=risk,
      )


  def _session_for_calls(*responses_per_call):
      """Return a session whose .post returns these responses in order."""
      sess = MagicMock()
      mocks = []
      for r in responses_per_call:
          mr = MagicMock()
          mr.json.return_value = r
          mr.headers = {"content-type": "application/json"}
          mr.text = json.dumps(r)
          mocks.append(mr)
      sess.post.side_effect = mocks
      return sess


  def test_place_parlays_happy_path(monkeypatch):
      monkeypatch.setenv("WAGERZON_PASSWORD", "secret")
      sess = _session_for_calls(CONFIRM_OK_RESPONSE, POST_OK_RESPONSE)
      monkeypatch.setattr(parlay_placer, "_get_session", lambda: sess)
      results = parlay_placer.place_parlays([_spec()])
      assert results[0].status == "placed"


  def test_place_parlays_price_drift_aborts(monkeypatch):
      monkeypatch.setenv("WAGERZON_PASSWORD", "secret")
      drifted_response = {
          "result": {"details": [{"Win": 28.50, "Risk": 15.0}], "Confirm": True}
      }
      sess = _session_for_calls(drifted_response)  # only the preflight call expected
      monkeypatch.setattr(parlay_placer, "_get_session", lambda: sess)
      results = parlay_placer.place_parlays([_spec(win=30.0)])
      assert results[0].status == "price_moved"
      assert "28.5" in results[0].error_msg
      # PostWager should NOT have been called
      assert sess.post.call_count == 1


  def test_place_parlays_dry_run_skips_post(monkeypatch):
      monkeypatch.setenv("WAGERZON_PASSWORD", "secret")
      sess = _session_for_calls(CONFIRM_OK_RESPONSE)  # only preflight
      monkeypatch.setattr(parlay_placer, "_get_session", lambda: sess)
      results = parlay_placer.place_parlays([_spec()], dry_run=True)
      assert results[0].status == "would_place"
      assert results[0].actual_win == 30.0
      assert sess.post.call_count == 1


  def test_place_parlays_auth_retry_succeeds(monkeypatch):
      """First preflight returns HTML (session expired) → re-login → retry → OK."""
      monkeypatch.setenv("WAGERZON_USERNAME", "u")
      monkeypatch.setenv("WAGERZON_PASSWORD", "p")
      html_resp = MagicMock()
      html_resp.headers = {"content-type": "text/html"}
      html_resp.text = "<html>login</html>"
      ok_resp = MagicMock()
      ok_resp.headers = {"content-type": "application/json"}
      ok_resp.json.return_value = CONFIRM_OK_RESPONSE
      ok_resp.text = json.dumps(CONFIRM_OK_RESPONSE)
      post_ok = MagicMock()
      post_ok.headers = {"content-type": "application/json"}
      post_ok.json.return_value = POST_OK_RESPONSE
      post_ok.text = json.dumps(POST_OK_RESPONSE)

      session1 = MagicMock(); session1.post.return_value = html_resp
      session2 = MagicMock(); session2.post.side_effect = [ok_resp, post_ok]
      session_iter = iter([session1, session2])
      monkeypatch.setattr(parlay_placer, "_get_session",
                          lambda: next(session_iter))
      monkeypatch.setattr(parlay_placer, "_clear_session", lambda: None)
      results = parlay_placer.place_parlays([_spec()])
      assert results[0].status == "placed"


  def test_place_parlays_network_error_no_retry(monkeypatch):
      sess = MagicMock()
      sess.post.side_effect = parlay_placer.requests.exceptions.Timeout()
      monkeypatch.setattr(parlay_placer, "_get_session", lambda: sess)
      monkeypatch.setenv("WAGERZON_PASSWORD", "p")
      results = parlay_placer.place_parlays([_spec()])
      assert results[0].status == "network_error"
      assert sess.post.call_count == 1  # no retry
  ```

- [ ] **Step 2: Run, verify fail**

  ```bash
  pytest test_parlay_placer.py -v -k "place_parlays"
  ```
  Expected: FAIL — `place_parlays` not defined.

- [ ] **Step 3: Implement `place_parlays` with retry, drift, dry-run, network handling**

  Add to `parlay_placer.py`:
  ```python
  def _drift_error_msg(expected: float, actual: float) -> str:
      return f"Expected ${expected:.2f}, Wagerzon offered ${actual:.2f}"


  def _preflight_with_retry(spec: ParlaySpec) -> tuple[float, float]:
      """ConfirmWagerHelper with one re-login retry on AuthExpired."""
      try:
          return _confirm_preflight(spec.legs, spec.amount)
      except AuthExpired:
          _clear_session()
          return _confirm_preflight(spec.legs, spec.amount)


  def _post_with_retry(specs: list[ParlaySpec]) -> list[PlacementResult]:
      """PostWagerMultipleHelper with one re-login + re-preflight retry."""
      try:
          return _post_wagers(specs)
      except AuthExpired:
          _clear_session()
          # Defensive: re-run preflight on each spec to reconfirm price
          for s in specs:
              win, _ = _confirm_preflight(s.legs, s.amount)
              if not _drift_ok(s.expected_win, win):
                  # Drift appeared during the retry window
                  return [
                      PlacementResult(
                          parlay_hash=s.parlay_hash, status="price_moved",
                          error_msg=_drift_error_msg(s.expected_win, win),
                      )
                      for s in specs
                  ]
          return _post_wagers(specs)


  def place_parlays(specs: list[ParlaySpec], dry_run: bool = False) -> list[PlacementResult]:
      """Place a batch of parlays at Wagerzon.

      For each spec:
        1) ConfirmWagerHelper preflight (with one auth-retry)
        2) Drift check vs expected_win — abort that spec on drift > $0.01
      Then for surviving specs:
        3) If dry_run: return would_place results
        4) Else: PostWagerMultipleHelper (with one auth-retry that re-runs
           preflight defensively before placing).
      """
      results_by_hash: dict[str, PlacementResult] = {}
      survivors: list[ParlaySpec] = []

      # Preflight + drift
      for spec in specs:
          try:
              win, risk = _preflight_with_retry(spec)
          except AuthExpired:
              results_by_hash[spec.parlay_hash] = PlacementResult(
                  parlay_hash=spec.parlay_hash, status="auth_error",
                  error_msg="auth_error: session expired",
              )
              continue
          except requests.exceptions.RequestException as e:
              results_by_hash[spec.parlay_hash] = PlacementResult(
                  parlay_hash=spec.parlay_hash, status="network_error",
                  error_msg=f"network_error: {type(e).__name__}",
              )
              continue
          if not _drift_ok(spec.expected_win, win):
              results_by_hash[spec.parlay_hash] = PlacementResult(
                  parlay_hash=spec.parlay_hash, status="price_moved",
                  error_msg=_drift_error_msg(spec.expected_win, win),
                  actual_win=win, actual_risk=risk,
              )
              continue
          if dry_run:
              results_by_hash[spec.parlay_hash] = PlacementResult(
                  parlay_hash=spec.parlay_hash, status="would_place",
                  actual_win=win, actual_risk=risk,
              )
              continue
          survivors.append(spec)

      # Live placement
      if survivors:
          try:
              live_results = _post_with_retry(survivors)
          except AuthExpired:
              live_results = [
                  PlacementResult(parlay_hash=s.parlay_hash, status="auth_error",
                                  error_msg="auth_error: re-login failed")
                  for s in survivors
              ]
          except requests.exceptions.RequestException as e:
              live_results = [
                  PlacementResult(parlay_hash=s.parlay_hash, status="network_error",
                                  error_msg=f"network_error: {type(e).__name__}")
                  for s in survivors
              ]
          for r in live_results:
              results_by_hash[r.parlay_hash] = r

      # Preserve input order
      return [results_by_hash[s.parlay_hash] for s in specs]
  ```

- [ ] **Step 4: Run all tests**

  ```bash
  pytest test_parlay_placer.py -v
  ```
  Expected: all tests pass.

- [ ] **Step 5: Commit**

  ```bash
  git add wagerzon_odds/parlay_placer.py wagerzon_odds/test_parlay_placer.py
  git commit -m "feat(parlay_placer): top-level place_parlays with retry + dry-run"
  ```

---

## Task 9: Dashboard server endpoint — `/api/place-parlay` (live mode)

**Files:**
- Modify: `Answer Keys/MLB Dashboard/mlb_dashboard_server.py`
- Read for reference: `Answer Keys/CBB Dashboard/cbb_dashboard_server.py:417-510` (the `/api/place-bet` pattern)

- [ ] **Step 1: Sketch the route handler structure**

  Open `Answer Keys/MLB Dashboard/mlb_dashboard_server.py`. Find where existing routes are registered (Flask `@app.route(...)` decorators). Add a new route `POST /api/place-parlay` accepting JSON `{parlay_hash: str, dry_run: bool = false}`.

- [ ] **Step 2: Implement the live-mode handler**

  Add to `mlb_dashboard_server.py`:
  ```python
  import sys
  import json as _json
  from pathlib import Path

  # Make wagerzon_odds importable
  REPO_ROOT = Path(__file__).resolve().parents[2]
  sys.path.insert(0, str(REPO_ROOT / "wagerzon_odds"))
  import parlay_placer  # noqa: E402

  MLB_DB = REPO_ROOT / "Answer Keys" / "mlb.duckdb"
  DASHBOARD_DB = REPO_ROOT / "Answer Keys" / "MLB Dashboard" / "mlb_dashboard.duckdb"


  def _load_parlay_row(parlay_hash: str) -> dict | None:
      """Load one parlay opportunity row from mlb.duckdb (read-only)."""
      con = duckdb.connect(str(MLB_DB), read_only=True)
      try:
          row = con.execute(
              "SELECT * FROM mlb_parlay_opportunities WHERE parlay_hash = ?",
              [parlay_hash],
          ).fetchone()
          if not row:
              return None
          cols = [d[0] for d in con.description]
          return dict(zip(cols, row))
      finally:
          con.close()


  def _build_spec_from_row(row: dict) -> parlay_placer.ParlaySpec:
      """Translate a mlb_parlay_opportunities row into a ParlaySpec.

      NOTE: depends on Task 2's idgm-resolution finding. Adjust the column
      names below based on how idgm is stored.
      """
      # TODO_T2: replace these with the actual columns identified in Task 2
      legs = [
          parlay_placer.Leg(
              idgm=int(row["idgm"]),
              play=int(row["spread_play"]),  # 0=away, 1=home
              points=float(row["spread_line"]),
              odds=int(row["spread_price"]),
          ),
          parlay_placer.Leg(
              idgm=int(row["idgm"]),
              play=2 if row["total_side"] == "Over" else 3,
              points=-float(row["total_line"]) if row["total_side"] == "Over" else float(row["total_line"]),
              odds=int(row["total_price"]),
          ),
      ]
      amount = float(row["kelly_bet"])
      # Compute expected_win from wz_odds (American → decimal × amount)
      wz = int(row["wz_odds"])
      decimal = (wz / 100 + 1) if wz > 0 else (100 / -wz + 1)
      expected_win = round(amount * (decimal - 1), 2)
      return parlay_placer.ParlaySpec(
          parlay_hash=row["parlay_hash"],
          legs=legs,
          amount=amount,
          expected_win=expected_win,
          expected_risk=amount,
      )


  def _upsert_placed_parlay(result: parlay_placer.PlacementResult, row: dict) -> None:
      """Insert or update placed_parlays. Idempotent on parlay_hash."""
      con = duckdb.connect(str(DASHBOARD_DB))
      try:
          con.execute("""
              INSERT INTO placed_parlays (
                  parlay_hash, status, combo, game_id, game_time,
                  recommended_size, expected_odds, expected_win,
                  actual_size, actual_win, ticket_number, idwt,
                  legs_json, error_msg, error_msg_key, updated_at
              ) VALUES (?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, CURRENT_TIMESTAMP)
              ON CONFLICT (parlay_hash) DO UPDATE SET
                  status        = EXCLUDED.status,
                  actual_size   = EXCLUDED.actual_size,
                  actual_win    = EXCLUDED.actual_win,
                  ticket_number = EXCLUDED.ticket_number,
                  idwt          = EXCLUDED.idwt,
                  error_msg     = EXCLUDED.error_msg,
                  error_msg_key = EXCLUDED.error_msg_key,
                  updated_at    = CURRENT_TIMESTAMP
          """, [
              result.parlay_hash, result.status, row.get("combo"),
              row.get("game_id"), row.get("game_time"),
              float(row["kelly_bet"]), int(row["wz_odds"]), float(row.get("expected_win") or 0),
              result.actual_risk, result.actual_win, result.ticket_number, result.idwt,
              _json.dumps([{"play": l.play, "idgm": l.idgm, "points": l.points,
                            "odds": l.odds} for l in _build_spec_from_row(row).legs]),
              result.error_msg, result.error_msg_key,
          ])
      finally:
          con.close()


  @app.route("/api/place-parlay", methods=["POST"])
  def api_place_parlay():
      payload = request.get_json(force=True) or {}
      parlay_hash = payload.get("parlay_hash")
      dry_run = bool(payload.get("dry_run", False))
      if not parlay_hash:
          return jsonify({"status": "error", "error_msg": "missing parlay_hash"}), 400

      # Idempotency check (live only)
      if not dry_run:
          con = duckdb.connect(str(DASHBOARD_DB), read_only=True)
          try:
              existing = con.execute(
                  "SELECT status, ticket_number, error_msg FROM placed_parlays "
                  "WHERE parlay_hash = ? AND status IN ('placing', 'placed')",
                  [parlay_hash],
              ).fetchone()
          finally:
              con.close()
          if existing:
              return jsonify({"status": existing[0], "ticket_number": existing[1],
                              "error_msg": existing[2] or ""})

      row = _load_parlay_row(parlay_hash)
      if not row:
          return jsonify({"status": "error", "error_msg": "parlay not found"}), 404

      # Mark placing (live mode only)
      if not dry_run:
          con = duckdb.connect(str(DASHBOARD_DB))
          try:
              con.execute("""
                  INSERT INTO placed_parlays (parlay_hash, status, recommended_size,
                                              expected_odds, expected_win, combo,
                                              game_id, game_time, updated_at)
                  VALUES (?, 'placing', ?, ?, ?, ?, ?, ?, CURRENT_TIMESTAMP)
                  ON CONFLICT (parlay_hash) DO NOTHING
              """, [parlay_hash, float(row["kelly_bet"]), int(row["wz_odds"]),
                    None, row.get("combo"), row.get("game_id"), row.get("game_time")])
          finally:
              con.close()

      spec = _build_spec_from_row(row)
      results = parlay_placer.place_parlays([spec], dry_run=dry_run)
      result = results[0]

      if dry_run:
          # Skip DB write entirely in dry-run mode
          return jsonify({"status": result.status,
                          "ticket_number": None,
                          "error_msg": result.error_msg,
                          "actual_win": result.actual_win})

      try:
          _upsert_placed_parlay(result, row)
      except Exception as e:
          # Local DB write failed AFTER Wagerzon placement
          if result.status == "placed":
              _record_orphan(result, parlay_hash, e)
              return jsonify({"status": "orphaned",
                              "ticket_number": result.ticket_number,
                              "error_msg": f"orphan: {e}"})
          raise

      return jsonify({"status": result.status,
                      "ticket_number": result.ticket_number,
                      "error_msg": result.error_msg})
  ```

- [ ] **Step 3: Replace the `TODO_T2` placeholder with real column names**

  Based on the finding from Task 2, edit `_build_spec_from_row` to use the actual `idgm` source. If Task 2 step 5 ran (idgm added to the table), use `row["idgm"]` directly. If Task 2 found case (c) (separate join), implement the join here.

- [ ] **Step 4: Smoke-test the import path**

  ```bash
  cd "Answer Keys/MLB Dashboard"
  python3 -c "
  import sys
  sys.path.insert(0, '.')
  # Just verify the file imports cleanly (don't start the server)
  import importlib.util
  spec = importlib.util.spec_from_file_location('m', 'mlb_dashboard_server.py')
  mod = importlib.util.module_from_spec(spec)
  spec.loader.exec_module(mod)
  print('imports OK')
  "
  ```
  Expected: `imports OK`.

- [ ] **Step 5: Commit**

  ```bash
  git add "Answer Keys/MLB Dashboard/mlb_dashboard_server.py"
  git commit -m "feat: /api/place-parlay endpoint (live + dry-run modes)"
  ```

---

## Task 10: Dashboard server — orphan handling

**Files:**
- Modify: `Answer Keys/MLB Dashboard/mlb_dashboard_server.py`

- [ ] **Step 1: Implement `_record_orphan`**

  Add to `mlb_dashboard_server.py`:
  ```python
  def _record_orphan(result: parlay_placer.PlacementResult,
                     parlay_hash: str, db_error: Exception) -> None:
      """Wagerzon confirmed placement but local placed_parlays write failed.
      Insert a forensics row into placement_orphans + log loudly."""
      try:
          con = duckdb.connect(str(DASHBOARD_DB))
          try:
              con.execute("""
                  INSERT INTO placement_orphans
                      (idwt, ticket_number, parlay_hash, raw_response, error)
                  VALUES (?, ?, ?, ?, ?)
                  ON CONFLICT (idwt) DO NOTHING
              """, [result.idwt, result.ticket_number, parlay_hash,
                    result.raw_response, str(db_error)])
          finally:
              con.close()
      except Exception as inner:
          # Last resort: shout to stderr
          print(
              f"!! ORPHAN UNRECORDED: idwt={result.idwt} "
              f"ticket={result.ticket_number} parlay_hash={parlay_hash} "
              f"raw_response={result.raw_response!r} "
              f"original_error={db_error!r} orphan_write_error={inner!r}",
              file=sys.stderr, flush=True,
          )
          raise
      print(
          f"!! ORPHAN RECORDED: idwt={result.idwt} ticket={result.ticket_number} "
          f"parlay_hash={parlay_hash} reason={db_error!r}",
          file=sys.stderr, flush=True,
      )
  ```

- [ ] **Step 2: Re-import smoke test**

  ```bash
  cd "Answer Keys/MLB Dashboard"
  python3 -c "
  import importlib.util, sys
  sys.path.insert(0, '.')
  spec = importlib.util.spec_from_file_location('m', 'mlb_dashboard_server.py')
  mod = importlib.util.module_from_spec(spec)
  spec.loader.exec_module(mod)
  print('imports OK')
  "
  ```
  Expected: `imports OK`.

- [ ] **Step 3: Commit**

  ```bash
  git add "Answer Keys/MLB Dashboard/mlb_dashboard_server.py"
  git commit -m "feat: orphan recording on local DB write failure"
  ```

---

## Task 11: Dashboard UI — Place column + status rendering

**Files:**
- Modify: `Answer Keys/MLB Dashboard/mlb_dashboard.R` (specifically `create_parlays_table()` around lines 246-330 and the loader at 3105-3142 per the explorer)

- [ ] **Step 1: Locate the parlays table renderer and the data loader**

  Open `Answer Keys/MLB Dashboard/mlb_dashboard.R`. Find:
  - The function `create_parlays_table()` (around line 246) — the renderer.
  - The block that loads `mlb_parlay_opportunities` (around line 3105) — augment to also LEFT JOIN `placed_parlays`.

- [ ] **Step 2: Augment the loader to join placement state**

  Change the parlay-load query so each row also includes `placement_status`, `ticket_number`, and `placement_error` from `placed_parlays`:
  ```r
  parlay_query <- "
    ATTACH 'Answer Keys/MLB Dashboard/mlb_dashboard.duckdb' AS dash (READ_ONLY);
    SELECT op.*, pp.status AS placement_status,
           pp.ticket_number AS ticket_number,
           pp.error_msg AS placement_error
    FROM mlb_parlay_opportunities op
    LEFT JOIN dash.placed_parlays pp ON pp.parlay_hash = op.parlay_hash
  "
  parlay_opps <- dbGetQuery(con, parlay_query)
  ```
  (Adjust to match the existing connection pattern; use `tryCatch` + `on.exit(dbDisconnect)` per repo CLAUDE.md.)

- [ ] **Step 3: Add a Place column to `create_parlays_table()`**

  Inside the function that builds the table (around line 246), add a new column. Pseudocode for what should render per row:
  - If `placement_status == "placed"` → `placed · #<ticket_number>` as muted text.
  - Else if `placement_status %in% c("price_moved","rejected","auth_error","network_error","orphaned")` → red pill showing the `placement_error`.
  - Else (no row OR `would_place`) → an actual button labeled "Place" (or "Dry run" when toggle is on — Task 12 wires this) with `onclick` handler `placeParlay('<parlay_hash>')`.

  Example R / DT cell formatter:
  ```r
  render_place_cell <- function(status, ticket, err, parlay_hash, dry_run) {
    if (!is.na(status) && status == "placed") {
      return(paste0('<span class="placed">placed · #', ticket, '</span>'))
    }
    if (!is.na(status) && status %in% c("price_moved","rejected","auth_error",
                                         "network_error","orphaned")) {
      return(paste0('<span class="pill error">', htmltools::htmlEscape(err), '</span>'))
    }
    label <- if (isTRUE(dry_run)) "Dry run" else "Place"
    paste0('<button class="place-btn" onclick="placeParlay(\'',
           parlay_hash, '\', ', tolower(as.character(isTRUE(dry_run))),
           ')">', label, '</button>')
  }
  ```

- [ ] **Step 4: Add the JS handler for the Place button**

  In the dashboard's HTML head (or a `tags$script(...)` block), add:
  ```html
  <script>
  function placeParlay(parlayHash, dryRun) {
    fetch('/api/place-parlay', {
      method: 'POST',
      headers: {'Content-Type': 'application/json'},
      body: JSON.stringify({parlay_hash: parlayHash, dry_run: !!dryRun})
    })
    .then(r => r.json())
    .then(data => {
      // Trigger a Shiny refresh of the parlay table
      Shiny.setInputValue('parlay_placement_done',
                          {hash: parlayHash, status: data.status, ts: Date.now()});
      console.log('place result:', data);
    })
    .catch(err => console.error('place failed:', err));
  }
  </script>
  ```
  And on the server side, observe `input$parlay_placement_done` to re-run the parlay loader.

- [ ] **Step 5: Smoke-test the dashboard**

  Start the dashboard:
  ```bash
  Rscript "Answer Keys/MLB Dashboard/mlb_dashboard.R"
  ```
  Open in a browser. Confirm:
  - The parlays table loads
  - The Place column renders with buttons for un-placed parlays
  - Hovering the buttons shows clickable cursor (no JS errors in browser console)
  Don't actually click yet — that's Task 13/14.

- [ ] **Step 6: Commit**

  ```bash
  git add "Answer Keys/MLB Dashboard/mlb_dashboard.R"
  git commit -m "feat(dashboard): add Place column + status rendering to parlays table"
  ```

---

## Task 12: Dashboard UI — Dry run toggle

**Files:**
- Modify: `Answer Keys/MLB Dashboard/mlb_dashboard.R`

- [ ] **Step 1: Add the toggle near the parlays table**

  Add a `checkboxInput("parlay_dry_run", "Dry run (no real bet)", value = TRUE)` somewhere visible above the parlays table (e.g. in the same row as the edge filter). Default to `TRUE` for safety.

- [ ] **Step 2: Pass `dry_run` through to the renderer**

  In the renderer call (Task 11 step 3), use `input$parlay_dry_run` as the `dry_run` arg of `render_place_cell`. The button label and JS call automatically reflect the current toggle state.

- [ ] **Step 3: Smoke-test toggle**

  Restart the dashboard. Toggle the checkbox. Confirm:
  - Buttons relabel "Place" ↔ "Dry run"
  - The button's `onclick` arg flips between `true`/`false`
  - No JS console errors

- [ ] **Step 4: Commit**

  ```bash
  git add "Answer Keys/MLB Dashboard/mlb_dashboard.R"
  git commit -m "feat(dashboard): dry-run toggle for parlay placement"
  ```

---

## Task 13: Validation — Stage 2 dry-run via dashboard

This is a **manual smoke test**, not a code commit. Document observations in the worktree.

- [ ] **Step 1: Start the dashboard**

  ```bash
  Rscript "Answer Keys/MLB Dashboard/mlb_dashboard.R"
  ```

- [ ] **Step 2: Toggle "Dry run" ON, click Place on a real +EV parlay**

  Pick a parlay with positive `edge_pct`. Click. Watch:
  - Browser console shows the `/api/place-parlay` round-trip
  - Response includes `status: "would_place"` and an `actual_win` value
  - `placed_parlays` table is **NOT** updated (verify via DuckDB)

- [ ] **Step 3: Verify drift detection on a moving line**

  Find a parlay whose price has shifted since the last pricer run (or wait for a price change). Click Place (still in dry-run). Verify:
  - Response status is `would_place` BUT actual_win differs from displayed wz_odds — OR —
  - Response status is `price_moved` with the expected/actual error message

- [ ] **Step 4: Verify Wagerzon shows NO bets in history**

  Log into Wagerzon manually. Confirm no new bets were placed during this stage.

- [ ] **Step 5: Document any anomalies**

  Append observations to a checklist file or commit message. Examples:
  - Did the response shape match what tests expected?
  - Were there F5 parlays that broke `_build_spec_from_row`? (Spec §10 risk.)
  - Were there error paths not covered by tests?

  If anomalies found: fix them in a new task before proceeding.

- [ ] **Step 6: No commit unless code changes**

  This task may produce no commit if all paths work. If tweaks were needed, commit them with `fix:` prefix.

---

## Task 14: Validation — Stage 3 first $15 real bet via dashboard

**Real money. $15 minimum. Do not skip.**

- [ ] **Step 1: Verify Wagerzon balance has at least $15**

  Log into Wagerzon manually. Confirm balance.

- [ ] **Step 2: Pick a parlay** with `kelly_bet >= 15` and click Place with **Dry run OFF**

  If no parlay has a Kelly recommendation ≥ $15, you may need to (a) wait for a bigger one, or (b) for this validation only, temporarily edit the dashboard server to override `amount = 15.0` regardless of `kelly_bet`. Revert that temporary edit before merging.

- [ ] **Step 3: Watch the response**

  - Browser console shows `status: "placed"`, a real `ticket_number`, `actual_win` matches expected
  - Dashboard refreshes the row to show `placed · #<ticket>`

- [ ] **Step 4: Verify in Wagerzon**

  Log into Wagerzon. Confirm the bet appears in their bet history with the correct legs and amount.

- [ ] **Step 5: Verify `placed_parlays` row**

  ```bash
  python3 -c "
  import duckdb
  con = duckdb.connect('Answer Keys/MLB Dashboard/mlb_dashboard.duckdb', read_only=True)
  print(con.execute('SELECT parlay_hash, status, ticket_number, idwt, actual_size, actual_win, placed_at FROM placed_parlays ORDER BY placed_at DESC LIMIT 1').fetchall())
  con.close()
  "
  ```
  Expected: row with `status='placed'`, ticket number, IDWT, correct amounts.

- [ ] **Step 6: Wait for `bet_logger/scraper_wagerzon.py` to run and verify Sheets update**

  Run it manually if you can't wait for cron:
  ```bash
  source bet_logger/venv/bin/activate
  python3 bet_logger/scraper_wagerzon.py
  deactivate
  ```
  Open the Google Sheet. Confirm the new bet appears in the MLB Summary tab (or wherever Wagerzon parlays land).

- [ ] **Step 7: Document validation**

  Append a "Stage 3 validation" line to this plan with the parlay_hash, ticket_number, and timestamp. Commit:
  ```bash
  git add docs/superpowers/plans/2026-04-26-mlb-parlay-auto-placement-plan.md
  git commit -m "docs(plans): record stage 3 first-bet validation"
  ```

---

## Task 15: Documentation

**Files:**
- Modify: `wagerzon_odds/README.md`
- Modify or create: `Answer Keys/MLB Dashboard/README.md`

- [ ] **Step 1: Update `wagerzon_odds/README.md`**

  Add sections documenting:
  - `parlay_placer.py` — what it does, public API (`place_parlays(specs, dry_run)`), how to run tests (`pytest test_parlay_placer.py -v`)
  - `recon_place_parlay.py` — when to use it (capture new placement endpoints), how to run, security note (scrub passwords)
  - `migrate_placed_parlays.py` — what it does, when to run

- [ ] **Step 2: Update or create `Answer Keys/MLB Dashboard/README.md`**

  Add section: "Auto-placement of correlated parlays" with:
  - How to use the Place button + Dry run toggle
  - Status meanings (placed / price_moved / rejected / auth_error / network_error / orphaned)
  - What to do on failure (especially the network_error ambiguity warning from spec §6.4)
  - Where state lives (`placed_parlays`, `placement_orphans` in `mlb_dashboard.duckdb`)

- [ ] **Step 3: Commit**

  ```bash
  git add wagerzon_odds/README.md "Answer Keys/MLB Dashboard/README.md"
  git commit -m "docs: parlay_placer + dashboard auto-placement usage"
  ```

---

## Task 16: Pre-merge review (REQUIRED per repo CLAUDE.md)

**Files:** none modified — review only.

- [ ] **Step 1: Generate the full diff**

  ```bash
  git diff main..HEAD --stat
  git diff main..HEAD > /tmp/mlb-parlay-auto-place.diff
  wc -l /tmp/mlb-parlay-auto-place.diff
  ```

- [ ] **Step 2: Walk the review checklist from repo CLAUDE.md**

  For each item, write a one-line note (PASS / fix needed):
  - **Data integrity**: Does `placed_parlays` get duplicate writes? Does idempotency actually fire? Are in-progress (`status='placing'`) records filtered from "needs placement" queries?
  - **Resource safety**: Every `duckdb.connect(...)` uses try/finally close? No file lock leaks?
  - **Edge cases**: What happens if `mlb_parlay_opportunities` is empty? If a parlay has no `idgm`? If two `/api/place-parlay` requests arrive concurrently?
  - **Dead code**: Any unused imports / commented-out blocks left in?
  - **Log/disk hygiene**: Are placement responses written to logs? (They shouldn't be — they contain the password.)
  - **Security**: `WAGERZON_PASSWORD` never logged? Recon JSON not committed?

- [ ] **Step 3: Document findings**

  Append the review notes to a new section "Pre-merge review" at the bottom of this plan file. Commit:
  ```bash
  git add docs/superpowers/plans/2026-04-26-mlb-parlay-auto-placement-plan.md
  git commit -m "docs(plans): pre-merge review notes"
  ```

- [ ] **Step 4: Address every fix-needed item**

  Each fix is its own commit with `fix:` prefix. Re-run the review until all items PASS.

---

## Task 17: Merge approval and cleanup

**Files:** none modified.

- [ ] **Step 1: Stop and ask the user for explicit merge approval**

  Per repo CLAUDE.md: "Approval required: Never merge to `main` or push to remote without explicit user approval."

  Show:
  - The full diff stat
  - A summary of stages 2 and 3 validation results (parlay_hash, ticket_number)
  - Any open concerns from the pre-merge review

  Wait for user "yes."

- [ ] **Step 2: Merge feature branch into main**

  From the worktree:
  ```bash
  git checkout main
  git merge --no-ff feature/mlb-parlay-auto-place -m "feat: MLB parlay auto-placement (Wagerzon, pure-API)"
  git log --oneline -5
  ```

- [ ] **Step 3: Clean up worktree and feature branch**

  ```bash
  cd /Users/callancapitolo/NFLWork
  git worktree remove ../NFLWork-mlb-parlay-place
  git branch -d feature/mlb-parlay-auto-place
  git worktree list
  ```
  Expected: only the main worktree remains.

- [ ] **Step 4: Do not push to remote unless user asks**

  Per repo CLAUDE.md: "DO NOT push to the remote repository unless the user explicitly asks you to do so."

---

## Self-review (post-write)

Spec coverage check (each spec section → which task implements it):

- §1 Mission & Problem — context only, no task needed
- §2 Scope (in scope) — Task 4-12 collectively cover all in-scope items
- §2 Scope (out of scope) — explicitly NOT in any task ✓
- §3 Architecture — Task 4-10 implement
- §4.1 parlay_placer.py — Tasks 4, 5, 6, 7, 8 ✓
- §4.2 mlb_dashboard_server.py — Tasks 9, 10 ✓
- §4.3 mlb_dashboard.R — Tasks 11, 12 ✓
- §5 Data model — Task 3 (migration) ✓
- §6 Failure handling — covered across Tasks 6, 7, 8, 10 ✓
- §7.1 Pre-flight (recon scrub, F5 verification) — recon already scrubbed; F5 verification is folded into Task 13 step 5 ✓
- §7.2 Test stages — Stage 1 (Tasks 4-8 unit tests), Stage 2 (Task 13), Stage 3 (Task 14) ✓
- §7.3 Rollout sequence — Tasks 13-17 in order ✓
- §8 Version control & docs — Task 1 (worktree/branch), Task 15 (docs) ✓
- §9 Follow-ups — explicitly out of scope, no tasks (correct) ✓
- §10 Open questions — `idgm` resolution is Task 2 ✓; placed_parlays schema verification is Task 3 ✓; ErrorMsgKey discovery is implicit in §6.2 raw passthrough; polling cadence revisited only if observed in Task 13/14

Placeholder scan: one intentional `TODO_T2` marker in Task 9 step 2 — flagged because Task 9 step 3 explicitly replaces it with the Task 2 finding. Documented in-line. Not a real placeholder.

Type consistency:
- `Leg.idgm: int`, `play: int`, `points: float`, `odds: int`, `pitcher: int` — used consistently across Tasks 4, 6, 7, 9.
- `ParlaySpec.amount: float`, `expected_win: float`, `expected_risk: float` — consistent.
- `PlacementResult.status: str` values are: `placed | price_moved | rejected | auth_error | network_error | orphaned | would_place` — matches spec §5.1 state machine plus the dry-run `would_place`. ✓
- Function names consistent: `_get_session`, `_clear_session`, `_confirm_preflight`, `_post_wagers`, `_preflight_with_retry`, `_post_with_retry`, `place_parlays`. No drift between tasks.

No gaps found.
