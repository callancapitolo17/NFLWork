"""Runner-log retention on the legacy subprocess path (issue #38, item 4).

``run_scrapers`` opened each per-scraper log with mode "w", so every cycle
erased the previous cycle's output. A book that died at 14:02 left no trace
by 14:03 — part of why four weeks of rot went unnoticed.

This path is LEGACY (dashboard-only; the bots price in-process). The fix is
deliberately minimal: append instead of truncate, with a size cap so an
unattended dashboard cannot fill the disk. No rotation, no rollover files.
"""
from __future__ import annotations

import sys
from pathlib import Path

from kalshi_common import sgp_runner


def _fake_scraper(tmp_path: Path, name: str, text: str) -> None:
    (tmp_path / name).write_text(
        f"import sys\nsys.stdout.write({text!r})\n")


def test_runner_log_appends_across_cycles(tmp_path):
    _fake_scraper(tmp_path, "s.py", "cycle-output\n")
    for _ in range(2):
        rcs = sgp_runner.run_scrapers(
            scraper_dir=str(tmp_path), scraper_names=["s.py"],
            venv_python=sys.executable, timeout_sec=30)
        assert rcs == {"s.py": 0}

    body = (tmp_path / "s.py.runner.log").read_text()
    assert body.count("cycle-output") == 2, (
        "previous cycle's output was erased — the exact blindness #38 fixes")


def test_runner_log_is_truncated_past_the_size_cap(tmp_path, monkeypatch):
    monkeypatch.setattr(sgp_runner, "RUNNER_LOG_MAX_BYTES", 256)
    log_path = tmp_path / "s.py.runner.log"
    log_path.write_text("x" * 1024)
    _fake_scraper(tmp_path, "s.py", "fresh\n")

    sgp_runner.run_scrapers(scraper_dir=str(tmp_path), scraper_names=["s.py"],
                            venv_python=sys.executable, timeout_sec=30)

    body = log_path.read_text()
    assert "x" * 100 not in body, "oversized log was not reset"
    assert "fresh" in body


def test_runner_log_under_the_cap_is_kept(tmp_path, monkeypatch):
    monkeypatch.setattr(sgp_runner, "RUNNER_LOG_MAX_BYTES", 10_000)
    log_path = tmp_path / "s.py.runner.log"
    log_path.write_text("older-cycle\n")
    _fake_scraper(tmp_path, "s.py", "fresh\n")

    sgp_runner.run_scrapers(scraper_dir=str(tmp_path), scraper_names=["s.py"],
                            venv_python=sys.executable, timeout_sec=30)

    body = log_path.read_text()
    assert "older-cycle" in body and "fresh" in body
