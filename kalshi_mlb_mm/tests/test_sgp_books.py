"""SGP_BOOKS: the one list that decides which books the maker's SGP service
prices (on-demand flights + structure warming). ProphetX is out by default
(2026-09-15): its warming requests were 100% 403s and funded the WAF block."""
import importlib
import inspect


def test_default_excludes_prophetx_and_keeps_the_other_five():
    import kalshi_mlb_mm.config as cfg
    assert "prophetx" not in cfg.SGP_BOOKS
    assert set(cfg.SGP_BOOKS) == {"draftkings", "fanduel", "novig", "betmgm", "caesars"}


def test_env_override_parses_and_can_re_add_prophetx(monkeypatch):
    import kalshi_mlb_mm.config as cfg
    monkeypatch.setenv("SGP_BOOKS", " fanduel, prophetx ,novig ")
    reloaded = importlib.reload(cfg)
    try:
        assert reloaded.SGP_BOOKS == ("fanduel", "prophetx", "novig")
    finally:
        monkeypatch.delenv("SGP_BOOKS")
        importlib.reload(cfg)


def test_main_passes_sgp_books_to_the_service():
    from kalshi_mlb_mm import main
    assert "books=config.SGP_BOOKS" in inspect.getsource(main)
