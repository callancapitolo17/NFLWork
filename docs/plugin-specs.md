# Plugin Specifications

Standalone tools that should be built as MCP servers or utility scripts to eliminate recurring manual work.

---

## 1. Team Name Normalizer

**Problem**: Team name mismatches across sportsbooks are the #1 source of bugs (99 debugging sessions). Each book formats names differently and there's no single source of truth beyond The Odds API.

**Solution**: A shared Python module + DuckDB lookup table that all scrapers import.

**Design**:
```
team_normalizer/
├── normalizer.py        # Core logic
├── mappings.duckdb      # Canonical mappings table
├── learn.py             # Auto-learn new mappings from Odds API
└── test_normalizer.py
```

**Core table** (`mappings.duckdb`):
```sql
CREATE TABLE team_mappings (
    raw_name TEXT,          -- "SEA SEAHAWKS", "Seattle", "SEATTLE SEAHAWKS"
    book TEXT,              -- "wagerzon", "hoop88", "bfa", etc.
    sport TEXT,             -- "nfl", "cbb", "nba"
    canonical_name TEXT,    -- "Seattle Seahawks" (from Odds API)
    PRIMARY KEY (raw_name, book, sport)
);
```

**API**:
```python
from team_normalizer import normalize
name = normalize("SEA SEAHAWKS", book="wagerzon", sport="nfl")
# Returns: "Seattle Seahawks"

# Auto-learn: when a new name appears, fuzzy-match against canonical list
# Log unmatched names for manual review
```

**Priority**: HIGH — addresses the single most common debugging issue.

---

## 2. Sportsbook Auth Manager

**Problem**: Auth tokens expire constantly (54 debugging sessions). Pipeline fails silently when tokens are stale.

**Solution**: A pre-pipeline check that validates tokens and rotates them when possible.

**Design**:
```python
# auth_manager.py
def check_all_tokens() -> dict:
    """Returns {book: {status, token_age, expires_at}}"""

def refresh_token(book: str) -> bool:
    """Attempt to refresh. Returns True if successful."""

def is_healthy(book: str) -> bool:
    """Quick connectivity check."""
```

**Integration**: Called at the start of `run.py` before launching scrapers. If a token is expired and can be auto-refreshed, do it. If not, skip that scraper and log a warning (don't fail the whole pipeline).

**Priority**: HIGH — eliminates the #2 failure mode.

---

## 3. Odds Dedup Validator

**Problem**: Duplicate odds entries in DuckDB and duplicate bets in Google Sheets (35 debugging sessions).

**Solution**: A validation script that runs post-scraper to catch duplicates before they propagate.

**Design**:
```python
# dedup_validator.py
def validate_odds_db(db_path: str, sport: str) -> Report:
    """Check for duplicate (game, market, book, period) entries."""

def validate_bet_log(sheet_id: str) -> Report:
    """Check Google Sheet for duplicate bet entries."""

def auto_fix(db_path: str, strategy: str = "keep_latest") -> int:
    """Remove duplicates, keeping latest entry. Returns count removed."""
```

**Priority**: MEDIUM — annoying but not blocking.

---

## 4. Dashboard Health Monitor

**Problem**: After code changes, user manually checks dashboards by pasting screenshots (56 times).

**Solution**: A script that screenshots all running dashboards and reports status.

**Design**:
```python
# dashboard_monitor.py
DASHBOARDS = {
    "nfl": "http://localhost:8081",
    "cbb": "http://localhost:8082",
    "coaching": "http://localhost:8080",
    "draft": "http://localhost:8083",
}

def health_check() -> dict:
    """Ping all dashboards, return status."""

def screenshot_all(output_dir: str) -> list:
    """Screenshot each dashboard, return paths."""

def compare_to_baseline(current: str, baseline: str) -> float:
    """Visual diff score (0=identical, 1=completely different)."""
```

**Priority**: LOW — nice to have but not critical.
