"""Store under the app's threads — one DuckDB file, many callers at once."""
import threading
from datetime import datetime, timezone

from nfl_specials.store import Store


def test_threads_reading_and_writing_at_once_all_succeed(tmp_path):
    """2026-11-13: a pricing lane's append and the page's read attached the
    file at the same moment and DuckDB refused one ("Unique file handle
    conflict"), which failed that refresh's whole fecta board."""
    store = Store(tmp_path / "state.duckdb")
    errors = []

    def busy() -> None:
        for _ in range(10):
            try:
                store.save_settings(5000.0, 0.25, False)
                store.placements_since(datetime(2026, 1, 1, tzinfo=timezone.utc))
            except Exception as exc:  # the regression: any DuckDB error here
                errors.append(repr(exc))

    threads = [threading.Thread(target=busy) for _ in range(4)]
    for thread in threads:
        thread.start()
    for thread in threads:
        thread.join()
    assert errors == []
    assert store.settings()["bankroll"] == 5000.0
