from dataclasses import dataclass
from pathlib import Path

# repo root = parent of the coverage_audit/ package dir
# (coverage_audit/registry.py -> coverage_audit/ -> repo root)
REPO_ROOT = Path(__file__).resolve().parents[1]

@dataclass(frozen=True)
class Book:
    name: str
    db_path: Path
    table: str
    screen_name: str | None      # bookmaker label in mlb_bets_book_prices, or None
    expected_on_screen: bool
    recon_hint: str | None       # how the agent refreshes auth for this book

def _db(rel: str) -> Path:
    return (REPO_ROOT / rel).resolve()

BOOK_REGISTRY: list[Book] = [
    Book("wagerzon", _db("wagerzon_odds/wagerzon.duckdb"), "mlb_odds", "wagerzon", True, "wagerzon_odds/recon_wagerzon.py"),
    Book("wagerzon_specials", _db("wagerzon_odds/wagerzon.duckdb"), "wagerzon_specials", None, False, "wagerzon_odds/recon_wagerzon.py"),
    Book("hoop88", _db("hoop88_odds/hoop88.duckdb"), "mlb_odds", None, False, "hoop88_odds/recon_hoop88.py"),
    Book("bfa", _db("bfa_odds/bfa.duckdb"), "mlb_odds", None, False, "bfa_odds (single API call)"),
    Book("bookmaker", _db("bookmaker_odds/bookmaker.duckdb"), "mlb_odds", "bookmaker", True, "bookmaker_odds/recon_bookmaker.py"),
    Book("bet105", _db("bet105_odds/bet105.duckdb"), "mlb_odds", "bet105", True, "bet105_odds/recon_bet105.py"),
    Book("kalshi", _db("kalshi_odds/kalshi.duckdb"), "mlb_odds", None, False, None),
    Book("draftkings_singles", _db("dk_odds/dk.duckdb"), "mlb_odds", "draftkings", True, None),
    Book("fanduel_singles", _db("fd_odds/fd.duckdb"), "mlb_odds", "fanduel", True, None),
]

def get_books() -> list[Book]:
    return list(BOOK_REGISTRY)
