from coverage_audit.registry import get_books, BOOK_REGISTRY

def test_registry_has_nine_books():
    names = {b.name for b in get_books()}
    assert names == {
        "wagerzon", "wagerzon_specials", "hoop88", "bfa", "bookmaker",
        "bet105", "kalshi", "draftkings_singles", "fanduel_singles",
    }

def test_only_five_books_expected_on_screen():
    on_screen = {b.screen_name for b in get_books() if b.expected_on_screen}
    assert on_screen == {"wagerzon", "bookmaker", "bet105", "draftkings", "fanduel"}

def test_db_paths_are_absolute_and_named_correctly():
    by = {b.name: b for b in get_books()}
    assert by["fanduel_singles"].db_path.name == "fd.duckdb"
    assert by["draftkings_singles"].db_path.name == "dk.duckdb"
    assert by["wagerzon"].db_path.is_absolute()
    # wagerzon_specials reads a different table in the same DB
    assert by["wagerzon_specials"].db_path.name == "wagerzon.duckdb"
    assert by["wagerzon_specials"].table == "wagerzon_specials"
