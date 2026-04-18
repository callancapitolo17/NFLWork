"""Stress test: writer + reader threads, asserts no lock errors."""
import threading
import time

from nfl_draft.lib import db as db_module


def test_writer_and_reader_coexist(monkeypatch, tmp_path):
    monkeypatch.setattr(db_module, "DB_PATH", tmp_path / "stress.duckdb")
    db_module.init_schema()
    errors = []
    stop = threading.Event()

    def writer():
        from datetime import datetime

        from nfl_draft.lib.db import write_connection

        i = 0
        while not stop.is_set():
            try:
                with write_connection() as con:
                    con.execute(
                        "INSERT INTO draft_odds VALUES (?, 'kalshi', 100, 0.5, 0.5, ?)",
                        [f"m{i}", datetime.now()],
                    )
                i += 1
            except Exception as e:
                errors.append(("writer", str(e)))

    def reader():
        from nfl_draft.lib.db import read_connection

        while not stop.is_set():
            try:
                with read_connection() as con:
                    con.execute("SELECT COUNT(*) FROM draft_odds").fetchone()
            except Exception as e:
                errors.append(("reader", str(e)))

    # Seed at least one row first so the reader has something to query
    from datetime import datetime

    from nfl_draft.lib.db import write_connection

    with write_connection() as con:
        con.execute(
            "INSERT INTO draft_odds VALUES ('seed', 'kalshi', 100, 0.5, 0.5, ?)",
            [datetime.now()],
        )

    w = threading.Thread(target=writer)
    r = threading.Thread(target=reader)
    w.start()
    r.start()
    time.sleep(5)  # 5-second stress
    stop.set()
    w.join()
    r.join()
    assert not errors, f"Lock errors observed: {errors[:5]}"
