import pytest
from nfl_draft.lib import db as db_module
from nfl_draft.lib import seed
from nfl_draft.lib.normalize import resolve_player, resolve_team


@pytest.fixture
def seeded(monkeypatch, tmp_path):
    monkeypatch.setattr(db_module, "DB_PATH", tmp_path / "test.duckdb")
    db_module.init_schema()
    seed.run()


def test_resolve_player_known_canonical(seeded):
    assert resolve_player("Cam Ward") == "Cam Ward"


def test_resolve_player_alias(seeded):
    assert resolve_player("Cameron Ward") == "Cam Ward"


def test_resolve_player_case_insensitive(seeded):
    assert resolve_player("cAm WaRd") == "Cam Ward"


def test_resolve_player_whitespace_tolerant(seeded):
    assert resolve_player("  Cam Ward  ") == "Cam Ward"


def test_resolve_player_unknown_returns_none(seeded):
    assert resolve_player("Joe Nobody") is None


def test_resolve_team_canonical(seeded):
    assert resolve_team("CHI") == "CHI"


def test_resolve_team_alias(seeded):
    assert resolve_team("Chicago Bears") == "CHI"


def test_resolve_team_unknown_returns_none(seeded):
    assert resolve_team("Detroit Pistons") is None
