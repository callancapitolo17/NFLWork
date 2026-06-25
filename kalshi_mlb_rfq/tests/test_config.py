"""Tests for kalshi_mlb_rfq/config.py — worktree path + new env knobs."""
from pathlib import Path

import kalshi_mlb_rfq.config as config_mod


def test_project_root_strips_worktree_suffix():
    """When config.py lives in a .worktrees/ subdir, PROJECT_ROOT should
    resolve to the main repo, not to the worktree itself."""
    pkg_dir = Path(config_mod.PKG_DIR)
    if ".worktrees" in str(pkg_dir):
        assert ".worktrees" not in str(config_mod.PROJECT_ROOT), \
            f"PROJECT_ROOT={config_mod.PROJECT_ROOT} still contains worktree path"
        assert str(config_mod.PROJECT_ROOT).endswith("NFLWork"), \
            f"PROJECT_ROOT should end with NFLWork, got {config_mod.PROJECT_ROOT}"


def test_new_env_knobs_have_defaults():
    assert config_mod.SGP_REFRESH_SEC == 60
    assert config_mod.SGP_SCRAPER_TIMEOUT_SEC == 90
    assert config_mod.MIN_BOOK_COUNT_FOR_BLEND == 2
    assert isinstance(config_mod.BOT_MARKET_DB, Path)
    assert str(config_mod.BOT_MARKET_DB).endswith("kalshi_mlb_rfq_market.duckdb")
    assert isinstance(config_mod.MLB_SGP_DIR, Path)


def test_mlb_sgp_dir_resolves_to_main_repo():
    assert config_mod.MLB_SGP_DIR.name == "mlb_sgp"
    assert ".worktrees" not in str(config_mod.MLB_SGP_DIR)


def test_logging_knobs_have_defaults():
    assert config_mod.LOG_MAX_BYTES == 50 * 1024 * 1024
    assert config_mod.LOG_BACKUP_COUNT == 5
    assert config_mod.LOG_LEVEL == "INFO"


def test_research_knobs_have_defaults():
    assert str(config_mod.RESEARCH_DB_PATH).endswith(
        "kalshi_mlb_rfq_research.duckdb")
    assert config_mod.RESEARCH_CANDIDATE_SAMPLING == 1.0
    assert config_mod.RESEARCH_BUFFER_MAX == 50000


def test_use_model_defaults_off(monkeypatch):
    import importlib
    # Unset any ambient USE_MODEL env var so we test the true default (False).
    monkeypatch.delenv("USE_MODEL", raising=False)
    import kalshi_mlb_rfq.config as c
    importlib.reload(c)
    assert c.USE_MODEL is False


def test_use_model_parses_truthy(monkeypatch):
    import importlib
    monkeypatch.setenv("USE_MODEL", "true")
    import kalshi_mlb_rfq.config as c
    importlib.reload(c)
    assert c.USE_MODEL is True
    monkeypatch.setenv("USE_MODEL", "0")
    importlib.reload(c)
    assert c.USE_MODEL is False
