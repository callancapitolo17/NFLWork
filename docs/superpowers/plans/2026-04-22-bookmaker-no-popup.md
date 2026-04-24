# Bookmaker Scraper — Stop Recurring Chrome Popup Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Stop `bookmaker_odds/scraper.py` from launching the interactive `recon_bookmaker.py` Chrome browser (a) when the session is actually fine and Bookmaker simply has no games posted, and (b) when the scraper is running inside a piped subprocess that can't handle interactive input (the MLB/CBB answer key pipeline).

**Architecture:** Three surgical changes inside `bookmaker_odds/scraper.py`:
1. After `session.get(SITE_URL)` returns 200 and `login()` returns 200, treat the session as healthy — empty schedule means no games posted, not a broken session. Save empty and return instead of escalating to recon.
2. Gate the interactive `refresh_cookies()` call behind `sys.stdin.isatty()`. When run from `run.py` (async subprocess with piped stdio), the scraper now logs a clear message and exits 0 with an empty table instead of hanging forever on a blocking `input()` inside Playwright.
3. Rate-limit recon launches with a `.last_recon_attempt` sentinel file (1h minimum gap). Prevents rapid repeated browser popups if the user re-runs the pipeline after a real CF 403.

No changes to `recon_bookmaker.py`, `run.py`, or any other scraper. All behavior except the popup itself is preserved: cookies, login, schedule fetch, DuckDB write.

**Tech Stack:** Python 3.11, curl_cffi, DuckDB, existing `bookmaker_odds/venv`.

---

## File Structure

**Modified:**
- `bookmaker_odds/scraper.py` — 3 small helpers + rewire the recovery branch of `scrape_bookmaker()`. Single file, single responsibility preserved.
- `bookmaker_odds/README.md` — one-paragraph note on the new "fail quietly when non-interactive" behavior.

**Created:**
- `bookmaker_odds/tests/__init__.py` — empty, marks test dir.
- `bookmaker_odds/tests/test_scraper_helpers.py` — unit tests for the three new helpers. These are pure-Python, no network/mocking needed.

**Not touched:**
- `recon_bookmaker.py` — remains interactive. That's fine; it's only invoked now when the user explicitly sits at a terminal.
- `run.py` — doesn't need to know anything. The scraper now cooperates by exiting 0 with empty data in bad conditions.

Cookie file (`.bookmaker_cookies.json`) and sentinel (`.last_recon_attempt`) stay inside `bookmaker_odds/`, both already gitignored by the broad `*.json` and hidden-file rules in that dir's conventions (verified: `recon_bookmaker_cookies.json` is not tracked).

## Version Control

- **Branch:** `fix/bookmaker-no-popup`
- **Worktree path:** `~/NFLWork-worktrees/bookmaker-no-popup`
- **Commit shape:** Four commits: (1) tests + helpers, (2) wire helpers into scrape_bookmaker, (3) README, (4) any fixup from manual validation.
- **Merge:** Fast-forward or merge commit to `main` after manual validation + user approval. Delete branch and worktree immediately after merge.

## Documentation

- `bookmaker_odds/README.md` — add a short "Non-interactive behavior" subsection under Auth explaining: scraper never pops a browser when stdin isn't a TTY, rate-limits recon to once per hour, and no-games-posted is not treated as an auth failure.
- `CLAUDE.md` (root) — no update needed; this is a scraper-local behavior change with no cross-cutting implications.
- `Answer Keys/CLAUDE.md` — no update needed; the R side and `run.py` contract is unchanged.

---

### Task 1: Unit tests for the three new helpers

**Files:**
- Create: `bookmaker_odds/tests/__init__.py`
- Create: `bookmaker_odds/tests/test_scraper_helpers.py`
- Test: same file (it *is* the test)

- [ ] **Step 1: Create the empty package marker**

```bash
touch /Users/callancapitolo/NFLWork/bookmaker_odds/tests/__init__.py
```

- [ ] **Step 2: Write the failing tests**

Create `bookmaker_odds/tests/test_scraper_helpers.py`:

```python
"""Unit tests for non-network helpers in scraper.py.

These tests are pure-Python — no curl_cffi, no subprocess, no DuckDB. They
pin down the three behaviors the 'stop the popup' fix relies on:

  1. _session_looks_healthy()        — distinguishes "no games posted" from
                                       "session is broken"
  2. _can_launch_interactive_recon() — refuses to launch when stdin isn't a
                                       TTY, so piped subprocesses can't hang
  3. _recon_rate_limited()           — skips recon if another attempt
                                       happened < 1h ago, so a misbehaving
                                       caller can't storm the user with
                                       Chrome popups
"""
import io
import os
import sys
import time
from pathlib import Path

import pytest

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))
import scraper  # noqa: E402


# ---- _session_looks_healthy ------------------------------------------------

def test_session_healthy_when_not_blocked_and_login_ok():
    assert scraper._session_looks_healthy(blocked=False, login_ok=True) is True


def test_session_unhealthy_when_blocked():
    assert scraper._session_looks_healthy(blocked=True, login_ok=True) is False


def test_session_unhealthy_when_login_failed():
    assert scraper._session_looks_healthy(blocked=False, login_ok=False) is False


# ---- _can_launch_interactive_recon ----------------------------------------

def test_recon_blocked_when_stdin_not_tty(monkeypatch):
    fake_stdin = io.StringIO("")  # StringIO is not a tty
    monkeypatch.setattr(sys, "stdin", fake_stdin)
    assert scraper._can_launch_interactive_recon() is False


def test_recon_allowed_when_stdin_is_tty(monkeypatch):
    class FakeTTY:
        def isatty(self):
            return True
    monkeypatch.setattr(sys, "stdin", FakeTTY())
    assert scraper._can_launch_interactive_recon() is True


# ---- _recon_rate_limited ---------------------------------------------------

def test_rate_limit_false_when_sentinel_missing(tmp_path):
    sentinel = tmp_path / "never_existed"
    assert scraper._recon_rate_limited(sentinel, min_gap_sec=3600) is False


def test_rate_limit_true_when_sentinel_fresh(tmp_path):
    sentinel = tmp_path / "fresh"
    sentinel.touch()
    # Just now — less than an hour ago
    assert scraper._recon_rate_limited(sentinel, min_gap_sec=3600) is True


def test_rate_limit_false_when_sentinel_old(tmp_path):
    sentinel = tmp_path / "old"
    sentinel.touch()
    # Rewind mtime by 2 hours
    old = time.time() - 7200
    os.utime(sentinel, (old, old))
    assert scraper._recon_rate_limited(sentinel, min_gap_sec=3600) is False
```

- [ ] **Step 3: Run tests to verify they fail (helpers don't exist yet)**

```bash
cd /Users/callancapitolo/NFLWork/bookmaker_odds
./venv/bin/python -m pytest tests/test_scraper_helpers.py -v
```

Expected: all tests fail with `AttributeError: module 'scraper' has no attribute '_session_looks_healthy'` (and similar for the other two).

- [ ] **Step 4: Implement the three helpers in scraper.py**

Add these three functions near the other module-level helpers in `bookmaker_odds/scraper.py`, directly above the `def refresh_cookies():` definition (currently line 172):

```python
def _session_looks_healthy(*, blocked: bool, login_ok: bool) -> bool:
    """Did we pass Cloudflare AND authenticate successfully?

    If both are true, an empty schedule means no games posted right now —
    not a broken session. Caller should save empty and exit cleanly rather
    than escalating to the interactive recon browser.
    """
    return (not blocked) and login_ok


def _can_launch_interactive_recon() -> bool:
    """Only allow the Playwright browser popup when a human is at the keyboard.

    recon_bookmaker.py has three blocking input() calls. When the scraper is
    running as a piped subprocess of run.py, stdin is not a TTY and those
    prompts would hang the whole MLB/CBB pipeline indefinitely.
    """
    try:
        return sys.stdin.isatty()
    except (AttributeError, ValueError):
        return False


def _recon_rate_limited(sentinel: Path, *, min_gap_sec: int = 3600) -> bool:
    """Return True if recon was attempted less than `min_gap_sec` ago.

    Prevents a misbehaving caller (e.g. a user hammering ./run.sh after a
    real CF 403) from spawning back-to-back Chrome windows.
    """
    if not sentinel.exists():
        return False
    age = time.time() - sentinel.stat().st_mtime
    return age < min_gap_sec
```

At the top of the file, add `import time` to the existing import block if it isn't already there (check — `datetime` is imported but `time` may not be). Leave `Path` alone; it's already imported.

- [ ] **Step 5: Run tests to verify they pass**

```bash
cd /Users/callancapitolo/NFLWork/bookmaker_odds
./venv/bin/python -m pytest tests/test_scraper_helpers.py -v
```

Expected: 7 passed.

- [ ] **Step 6: Commit**

```bash
git add bookmaker_odds/tests/__init__.py \
        bookmaker_odds/tests/test_scraper_helpers.py \
        bookmaker_odds/scraper.py
git commit -m "$(cat <<'EOF'
test(bookmaker): helpers to suppress Chrome popup storms

Three pure helpers in scraper.py, unit-tested without network:
- _session_looks_healthy: empty schedule after CF-OK + login-OK is just "no
  games", not a broken session.
- _can_launch_interactive_recon: stdin-isatty gate so piped subprocesses
  can't hang on recon_bookmaker.py's blocking input() calls.
- _recon_rate_limited: 1h gap between recon attempts to prevent repeated
  browser popups.

Wiring into scrape_bookmaker() comes in the next commit.

Co-Authored-By: Claude Opus 4.7 (1M context) <noreply@anthropic.com>
EOF
)"
```

---

### Task 2: Wire the helpers into scrape_bookmaker

**Files:**
- Modify: `bookmaker_odds/scraper.py:412-469` (the `scrape_bookmaker` recovery branch)

- [ ] **Step 1: Replace the recovery branch**

In `bookmaker_odds/scraper.py`, locate the current recovery block at lines 439-469 (starts with `# Test with first league` and ends with `test_data = fetch_schedule(session, first_league["id"])`/`return []`). Replace it with this block.

Find exactly this existing block:

```python
    # Test with first league
    first_league = config["leagues"][0]
    test_data = None if blocked else fetch_schedule(session, first_league["id"])

    # If blocked or no data, try refreshing cookies via recon
    if not _has_games(test_data):
        if blocked:
            print("Cloudflare blocked request — cookies expired.")
        else:
            # Not blocked but no games — try login first
            if BOOKMAKER_USERNAME and BOOKMAKER_PASSWORD:
                print("No games returned, logging in...")
                if login(session, BOOKMAKER_USERNAME, BOOKMAKER_PASSWORD):
                    test_data = fetch_schedule(session, first_league["id"])

        # Still no games? Refresh cookies via recon
        if not _has_games(test_data):
            refresh_cookies()
            session = _create_session()
            resp = session.get(SITE_URL, timeout=15)
            if resp.status_code == 403:
                print("ERROR: Still blocked after recon. Check Cloudflare manually. Clearing stale data.")
                save_to_database(sport, [])
                return []
            if BOOKMAKER_USERNAME and BOOKMAKER_PASSWORD:
                login(session, BOOKMAKER_USERNAME, BOOKMAKER_PASSWORD)
            test_data = fetch_schedule(session, first_league["id"])
            if not _has_games(test_data):
                print("ERROR: No games returned after recon + login. Clearing stale data.")
                save_to_database(sport, [])
                return []
```

Replace it with:

```python
    # Test with first league
    first_league = config["leagues"][0]
    test_data = None if blocked else fetch_schedule(session, first_league["id"])

    login_ok = True  # Default: assume prior-session cookies are still valid

    # If we didn't find games on the first probe, try login before anything drastic
    if not _has_games(test_data):
        login_ok = False
        if not blocked and BOOKMAKER_USERNAME and BOOKMAKER_PASSWORD:
            print("No games returned, logging in...")
            login_ok = login(session, BOOKMAKER_USERNAME, BOOKMAKER_PASSWORD)
            if login_ok:
                test_data = fetch_schedule(session, first_league["id"])

    # If still no games, decide WHY before launching a browser popup.
    if not _has_games(test_data):
        # Case A: session is healthy, BM just has nothing posted right now.
        # Happens at odd hours / off days. No recon needed. Save empty so
        # any previous scrape's stale data gets cleared.
        if _session_looks_healthy(blocked=blocked, login_ok=login_ok):
            print("Session healthy but no games posted — saving empty.")
            save_to_database(sport, [])
            return []

        # Case B: session is broken. Only allow the interactive recon popup
        # when a human is at the keyboard AND we haven't attempted recon in
        # the last hour. Otherwise log a clear message and exit cleanly —
        # the pipeline cannot hang on a blocking input().
        sentinel = Path(__file__).parent / ".last_recon_attempt"
        if not _can_launch_interactive_recon():
            print("Cookies appear stale but stdin is not a TTY.")
            print("Skipping interactive recon. Run manually:")
            print(f"  cd {Path(__file__).parent} && ./venv/bin/python recon_bookmaker.py")
            save_to_database(sport, [])
            return []
        if _recon_rate_limited(sentinel):
            print("Recon was attempted recently (< 1h ago). Skipping to avoid popup loop.")
            save_to_database(sport, [])
            return []

        # Green light: touch the sentinel, run recon, retry once.
        sentinel.touch()
        refresh_cookies()
        session = _create_session()
        resp = session.get(SITE_URL, timeout=15)
        if resp.status_code == 403:
            print("ERROR: Still blocked after recon. Check Cloudflare manually. Clearing stale data.")
            save_to_database(sport, [])
            return []
        if BOOKMAKER_USERNAME and BOOKMAKER_PASSWORD:
            login(session, BOOKMAKER_USERNAME, BOOKMAKER_PASSWORD)
        test_data = fetch_schedule(session, first_league["id"])
        if not _has_games(test_data):
            print("ERROR: No games returned after recon + login. Clearing stale data.")
            save_to_database(sport, [])
            return []
```

- [ ] **Step 2: Smoke-test the standalone scraper from a terminal (Case A — no games posted)**

Make sure the fresh cookies you promoted earlier are still in place:

```bash
cd /Users/callancapitolo/NFLWork/bookmaker_odds
./venv/bin/python scraper.py mlb
```

Expected: either "Session healthy but no games posted — saving empty." (if BM is currently showing 0 games) OR a normal scrape with some records. Either way, **no Chrome window opens**. Exit code 0.

Verify exit code:

```bash
echo $?
```

Expected: `0`.

- [ ] **Step 3: Smoke-test the piped-subprocess path (Case "pretend stdin is not a TTY")**

This forces the non-TTY branch even if stdin is a real terminal, so we exercise the piped-subprocess behavior without having to run the full MLB dashboard:

```bash
cd /Users/callancapitolo/NFLWork/bookmaker_odds
# Simulate stale cookies by moving the good ones out of the way temporarily
mv .bookmaker_cookies.json .bookmaker_cookies.json.bak
# Pipe stdin from /dev/null so sys.stdin.isatty() returns False
./venv/bin/python scraper.py mlb < /dev/null
echo "EXIT=$?"
# Restore real cookies immediately
mv .bookmaker_cookies.json.bak .bookmaker_cookies.json
```

Expected output includes:
```
Cookies appear stale but stdin is not a TTY.
Skipping interactive recon. Run manually:
  cd .../bookmaker_odds && ./venv/bin/python recon_bookmaker.py
```
and `EXIT=0`. **No Chrome window opens.** No hang.

- [ ] **Step 4: Verify unit tests still pass**

```bash
cd /Users/callancapitolo/NFLWork/bookmaker_odds
./venv/bin/python -m pytest tests/test_scraper_helpers.py -v
```

Expected: 7 passed.

- [ ] **Step 5: Commit**

```bash
git add bookmaker_odds/scraper.py
git commit -m "$(cat <<'EOF'
fix(bookmaker): never pop Chrome when session is healthy or stdin is piped

Three behavior changes in scrape_bookmaker():
- Empty schedule after CF-OK + login-OK is now "no games posted"; save empty
  and return. Previously this was treated as a broken session and launched
  recon_bookmaker.py's interactive Playwright browser.
- Interactive recon is gated behind sys.stdin.isatty(). When run.py spawns
  the scraper as a piped subprocess, the scraper now logs a clear message
  and exits 0 instead of hanging on a blocking input() inside the recon
  script.
- Recon is rate-limited to once per hour via .last_recon_attempt sentinel
  so rapid re-runs can't open back-to-back Chrome windows.

Fixes the MLB Dashboard pipeline hang where bookmaker's recon subprocess
blocked on an input() whose prompt never reached the terminal.

Co-Authored-By: Claude Opus 4.7 (1M context) <noreply@anthropic.com>
EOF
)"
```

---

### Task 3: README update

**Files:**
- Modify: `bookmaker_odds/README.md`

- [ ] **Step 1: Replace the Auth section**

Find this block at the end of the README:

```markdown
## Auth

Requires in `.env`:
- `BOOKMAKER_USERNAME`
- `BOOKMAKER_PASSWORD`

Cookies cached in `.bookmaker_cookies.json`. If blocked/expired, runs `recon_bookmaker.py` to refresh.
```

Replace with:

```markdown
## Auth

Requires in `.env`:
- `BOOKMAKER_USERNAME`
- `BOOKMAKER_PASSWORD`

Cookies cached in `.bookmaker_cookies.json`. If Cloudflare blocks the request AND stdin is a TTY AND the last recon attempt was > 1h ago, the scraper launches `recon_bookmaker.py` in a real Chrome window to refresh `cf_clearance`.

Otherwise (piped subprocess from `run.py`, recent recon attempt, or healthy session with no games posted) the scraper clears stale data, logs the reason, and exits 0 without opening a browser. To refresh cookies manually when the pipeline has skipped it:

```bash
cd bookmaker_odds
./venv/bin/python recon_bookmaker.py
```
```

- [ ] **Step 2: Commit**

```bash
git add bookmaker_odds/README.md
git commit -m "$(cat <<'EOF'
docs(bookmaker): document non-interactive fail-quiet behavior

Co-Authored-By: Claude Opus 4.7 (1M context) <noreply@anthropic.com>
EOF
)"
```

---

### Task 4: Full-pipeline validation + merge

**Files:**
- None modified. Validates the feature end-to-end and cleans up.

- [ ] **Step 1: Run the MLB Dashboard pipeline end-to-end**

```bash
cd "/Users/callancapitolo/NFLWork/Answer Keys/MLB Dashboard"
./run.sh
```

Expected: `[scraper] ✓ bookmaker complete (Ns)` with no hang, no Chrome popup. If BM has no games posted for MLB right now, the scraper prints "Session healthy but no games posted — saving empty." and proceeds. The dashboard finishes building.

- [ ] **Step 2: Pre-merge review**

Per `CLAUDE.md` pre-merge checklist, from the worktree:

```bash
git diff main..HEAD
```

Scan for:
- Data integrity: `save_to_database(sport, [])` is still called in every bail path (✓ verify by reading the new block).
- Resource safety: no new DB connections added; existing `conn.close()` in `save_to_database` still fires.
- Edge cases: three branches covered (healthy-empty, non-TTY, rate-limited). Confirm each has an explicit test or smoke-test above.
- Dead code: no unused imports or flags introduced.
- Log hygiene: no secrets printed. `BOOKMAKER_USERNAME/PASSWORD` are never logged.
- Documentation: README updated in the same branch.

Document findings inline before merging. Ask the user for explicit approval to merge.

- [ ] **Step 3: Merge to main (ONLY after user approval)**

```bash
cd /Users/callancapitolo/NFLWork
git checkout main
git merge --no-ff fix/bookmaker-no-popup
```

- [ ] **Step 4: Clean up the worktree and branch**

```bash
git worktree remove ~/NFLWork-worktrees/bookmaker-no-popup
git branch -d fix/bookmaker-no-popup
git worktree list
```

Expected: only the main working tree is listed.

- [ ] **Step 5: Final commit hygiene check**

```bash
git status
git log --oneline -5
```

Expected: clean working tree, last 3 commits are the test/helper, fix, and docs commits from this plan.
