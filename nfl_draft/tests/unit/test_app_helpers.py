"""Dashboard startup/background-scrape helpers."""
import subprocess
import sys
from pathlib import Path


def test_run_scrape_once_uses_detached_popen(monkeypatch):
    """_run_scrape_once must spawn a detached subprocess that doesn't block."""
    calls = []

    class FakePopen:
        def __init__(self, *args, **kwargs):
            calls.append((args, kwargs))

    monkeypatch.setattr(subprocess, "Popen", FakePopen)

    # kalshi_draft/app.py uses sibling imports (`import db`), so its directory
    # must be on sys.path before import. It's not a real package.
    repo_root = Path(__file__).resolve().parents[3]
    kalshi_dir = repo_root / "kalshi_draft"
    monkeypatch.syspath_prepend(str(kalshi_dir))

    # Drop any cached app module so the patched Popen takes effect on reimport.
    sys.modules.pop("app", None)
    import importlib
    app_module = importlib.import_module("app")

    app_module._run_scrape_once()

    assert len(calls) == 1
    args, kwargs = calls[0]
    # Must use start_new_session to detach
    assert kwargs.get("start_new_session") is True
    # Must invoke run.py with --mode scrape --book all
    cmd = args[0]
    assert "--mode" in cmd and "scrape" in cmd and "--book" in cmd and "all" in cmd
