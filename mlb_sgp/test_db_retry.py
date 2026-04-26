"""Unit tests for db._connect_with_retry."""
import unittest
from unittest.mock import patch, MagicMock

import duckdb

import db


class TestConnectWithRetry(unittest.TestCase):
    def test_succeeds_on_first_try_when_no_lock(self):
        fake_con = MagicMock()
        with patch.object(db.duckdb, "connect", return_value=fake_con) as m:
            con = db._connect_with_retry("/tmp/test_sgp.duckdb")
        self.assertIs(con, fake_con)
        self.assertEqual(m.call_count, 1)

    def test_retries_then_succeeds_on_lock_error(self):
        fake_con = MagicMock()
        attempts = {"n": 0}

        def fake_connect(*args, **kwargs):
            attempts["n"] += 1
            if attempts["n"] < 3:
                raise duckdb.IOException("Could not set lock on file")
            return fake_con

        with patch.object(db.duckdb, "connect", side_effect=fake_connect):
            con = db._connect_with_retry(
                "/tmp/test_sgp.duckdb", base_delay=0.001
            )
        self.assertIs(con, fake_con)
        self.assertEqual(attempts["n"], 3)

    def test_raises_on_non_lock_error(self):
        with patch.object(db.duckdb, "connect",
                          side_effect=ValueError("bad path")):
            with self.assertRaises(ValueError):
                db._connect_with_retry("/bad/path", base_delay=0.001)

    def test_gives_up_after_max_attempts(self):
        with patch.object(db.duckdb, "connect",
                          side_effect=duckdb.IOException("lock")):
            with self.assertRaises(duckdb.IOException):
                db._connect_with_retry(
                    "/tmp/test_sgp.duckdb",
                    max_attempts=3,
                    base_delay=0.001,
                )


if __name__ == "__main__":
    unittest.main()
