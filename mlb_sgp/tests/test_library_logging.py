"""Library code logs, it does not print (issue #36).

Both bots run these scrapers IN-PROCESS via ``kalshi_common.sgp_service``.
A ``print()`` inside library code bypasses the bot's log handlers and levels
entirely, so sanity drops and event-match misses never reach ``bot.log``.

Scope note: this scan covers the modules the BOTS import. The ``recon_*.py``
diagnostics, ``quick_recon.py``, ``probe_concurrency.py`` and the Pikkit
modules are standalone command-line tools that no bot imports; they keep
their prints deliberately.
"""
from __future__ import annotations

import ast
import logging
from pathlib import Path

import pytest

MLB_SGP_DIR = Path(__file__).resolve().parents[1]

# Every module reachable from an in-process bot run.
BOT_PATH_MODULES = (
    # orchestrators
    "draftkings.py", "fanduel.py", "prophetx.py", "novig.py",
    "betmgm.py", "caesars.py",
    # clients
    "dk_client.py", "fd_client.py", "prophetx_client.py",
    "novig_client.py", "betmgm_client.py", "caesars_client.py",
    # support
    "caesars_waf.py", "db.py", "_shared.py", "dk_leg_resolvers.py",
    "integer_line_derivation.py",
    # CLI shims — imported as libraries by the orchestrators / SGPService
    "scraper_draftkings_sgp.py", "scraper_fanduel_sgp.py",
    "scraper_prophetx_sgp.py", "scraper_novig_sgp.py",
    "scraper_betmgm_sgp.py", "scraper_caesars_sgp.py",
    # singles scrapers (dashboard pill data; same in-process concerns)
    "scraper_draftkings_singles.py", "scraper_fanduel_singles.py",
)


def _main_guard_line_ranges(tree: ast.Module) -> list[tuple[int, int]]:
    """Line spans of every ``if __name__ == "__main__":`` block."""
    spans = []
    for node in ast.walk(tree):
        if not isinstance(node, ast.If):
            continue
        test = node.test
        if (isinstance(test, ast.Compare)
                and isinstance(test.left, ast.Name)
                and test.left.id == "__name__"):
            spans.append((node.lineno, max(
                getattr(n, "end_lineno", node.lineno) or node.lineno
                for n in ast.walk(node))))
    return spans


def _prints_outside_main(path: Path) -> list[int]:
    tree = ast.parse(path.read_text())
    spans = _main_guard_line_ranges(tree)
    out = []
    for node in ast.walk(tree):
        if (isinstance(node, ast.Call) and isinstance(node.func, ast.Name)
                and node.func.id == "print"):
            if not any(lo <= node.lineno <= hi for lo, hi in spans):
                out.append(node.lineno)
    return sorted(out)


@pytest.mark.parametrize("module_name", BOT_PATH_MODULES)
def test_no_print_outside_main_guard(module_name):
    path = MLB_SGP_DIR / module_name
    assert path.exists(), f"{module_name} moved — update BOT_PATH_MODULES"
    offenders = _prints_outside_main(path)
    assert offenders == [], (
        f"{module_name} print()s at lines {offenders} bypass the bots' log "
        f"handlers — use module logger instead (issue #36)")


@pytest.mark.parametrize("module_name", BOT_PATH_MODULES)
def test_module_defines_a_logger(module_name):
    """Every bot-path module owns a ``logging.getLogger(__name__)``, so the
    host process can filter per module."""
    source = (MLB_SGP_DIR / module_name).read_text()
    assert "logging.getLogger(__name__)" in source, (
        f"{module_name} has no module logger (issue #36)")


def test_library_modules_attach_no_handlers():
    """Library code must NOT configure logging — the host process decides
    where logs go (same principle as R's ``message()`` vs ``cat()``).
    ``basicConfig`` is allowed only inside a ``__main__`` guard."""
    for module_name in BOT_PATH_MODULES:
        path = MLB_SGP_DIR / module_name
        tree = ast.parse(path.read_text())
        spans = _main_guard_line_ranges(tree)
        for node in ast.walk(tree):
            if not isinstance(node, ast.Call):
                continue
            func = node.func
            name = getattr(func, "attr", None)
            if name in ("basicConfig", "addHandler"):
                assert any(lo <= node.lineno <= hi for lo, hi in spans), (
                    f"{module_name}:{node.lineno} configures logging at "
                    f"import time — only a __main__ block may do that")


def test_sanity_drop_reaches_a_standard_logging_handler(caplog, monkeypatch):
    """#36 acceptance criterion: a sanity-drop event emitted during an
    in-process run is captured by a standard logging handler (it used to go
    to stdout, where a bot's log handlers never saw it)."""
    from mlb_sgp.tests.sweep_drives import drive_caesars

    with caplog.at_level(logging.DEBUG, logger="mlb_sgp.caesars"):
        # naive = 1.9 * 1.9 = 3.61; 99.0 blows past SANITY_MULT_RATIO * naive.
        rows = drive_caesars(monkeypatch, healthy=True, combo_decimal=99.0,
                             leg_decimal=1.9)

    assert rows == []
    drops = [r for r in caplog.records if "SANITY-DROP" in r.getMessage()]
    assert drops, "sanity drop never reached a logging handler"
    assert all(r.name == "mlb_sgp.caesars" for r in drops)


def test_in_process_run_writes_nothing_to_stdout(capsys, monkeypatch):
    """The bots run these scrapers in-process; library stdout would bypass
    their handlers entirely."""
    from mlb_sgp.tests.sweep_drives import drive_caesars

    drive_caesars(monkeypatch, healthy=True)
    assert capsys.readouterr().out == ""
