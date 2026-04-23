"""Unit tests: _run_scrape_once must run in-process in a daemon thread,
not via a subprocess. Under the in-process writer consolidation, any
subprocess scrape would re-introduce cross-process DuckDB lock races."""
import importlib
import subprocess
import sys
import threading
import time
from pathlib import Path

# kalshi_draft/app.py uses sibling imports (`import db`), so its directory
# must be on sys.path before import. It is not a package (no __init__.py).
_REPO_ROOT = Path(__file__).resolve().parents[3]
_KALSHI_DIR = _REPO_ROOT / "kalshi_draft"


def _import_app(monkeypatch):
    """Import kalshi_draft/app.py with the correct sys.path so sibling
    imports resolve. Uses the same pattern as test_app_helpers.py."""
    monkeypatch.syspath_prepend(str(_KALSHI_DIR))
    # Drop any cached module so patches take full effect on each test.
    sys.modules.pop("app", None)
    return importlib.import_module("app")


def test_run_scrape_once_does_not_spawn_subprocess(monkeypatch):
    """If _run_scrape_once ever falls back to subprocess, the patched
    Popen/run will raise and this test will fail — acting as a tripwire."""
    app_mod = _import_app(monkeypatch)

    def _boom(*a, **kw):
        raise AssertionError("subprocess usage is forbidden under in-process consolidation")
    monkeypatch.setattr(subprocess, "Popen", _boom)
    monkeypatch.setattr(subprocess, "run", _boom)

    # Also stub run_scrape so it doesn't actually hit the network.
    import nfl_draft.run as nfl_run
    monkeypatch.setattr(nfl_run, "run_scrape", lambda book: None)

    app_mod._run_scrape_once()
    # Give the daemon thread a moment to actually call run_scrape.
    time.sleep(0.2)


def test_run_scrape_once_invokes_run_scrape_in_thread(monkeypatch):
    """The in-process refactor must delegate to nfl_draft.run.run_scrape('all')."""
    app_mod = _import_app(monkeypatch)

    called = {"count": 0, "arg": None, "thread_name": None}
    done = threading.Event()

    def _fake_run_scrape(book):
        called["count"] += 1
        called["arg"] = book
        called["thread_name"] = threading.current_thread().name
        done.set()

    import nfl_draft.run as nfl_run
    monkeypatch.setattr(nfl_run, "run_scrape", _fake_run_scrape)

    app_mod._run_scrape_once()

    assert done.wait(timeout=2.0), "scrape thread did not invoke run_scrape within 2s"
    assert called["count"] == 1
    assert called["arg"] == "all"
    # Must run in a NON-main thread so it doesn't block the caller.
    assert called["thread_name"] != "MainThread"


def test_run_scrape_once_survives_run_scrape_exception(monkeypatch):
    """If run_scrape raises, _run_scrape_once must not propagate the
    exception into the caller's thread (the _periodic_scrape_loop)."""
    app_mod = _import_app(monkeypatch)

    def _boom(book):
        raise RuntimeError("simulated scrape failure")

    import nfl_draft.run as nfl_run
    monkeypatch.setattr(nfl_run, "run_scrape", _boom)

    # Must not raise:
    app_mod._run_scrape_once()
    time.sleep(0.2)
