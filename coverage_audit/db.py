from pathlib import Path
import duckdb
from coverage_audit.registry import REPO_ROOT

COVERAGE_DB_PATH = (REPO_ROOT / "coverage_audit" / "coverage.duckdb").resolve()

def connect_readonly(path: Path):
    """Open a book DB read-only. Returns None on any failure (missing, locked, corrupt)."""
    if not Path(path).exists():
        return None
    try:
        return duckdb.connect(str(path), read_only=True)
    except Exception:
        return None

def connect_coverage():
    """Open (creating if needed) the coverage output DB read-write."""
    # Read the module attribute at call time (not a local import) so tests can
    # monkeypatch coverage_audit.db.COVERAGE_DB_PATH to redirect to a temp DB.
    COVERAGE_DB_PATH.parent.mkdir(parents=True, exist_ok=True)
    return duckdb.connect(str(COVERAGE_DB_PATH))
