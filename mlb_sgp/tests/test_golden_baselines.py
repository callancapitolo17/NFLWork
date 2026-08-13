"""Golden baselines for the six SGP orchestrators (issue #42).

Replays each book's frozen payloads through its real ``price_sgps`` and pins
the resulting rows. Fully offline: no live call, no DuckDB, no shared state.

These replace the two live-API baselines in ``test_sgp_regression.py``, which
re-ran the real scraper against the book, wrote into the shared production
DuckDB, and compared against a CSV keyed by ``game_id`` — so they went stale
within hours of capture and could never pass on a normal day.

Recapture with ``python3 mlb_sgp/tests/capture_goldens.py`` (see that file
for when that is and is not appropriate).
"""
from __future__ import annotations

import sys
from datetime import datetime
from pathlib import Path

import pytest

MLB_SGP_DIR = Path(__file__).resolve().parents[1]
if str(MLB_SGP_DIR) not in sys.path:
    sys.path.insert(0, str(MLB_SGP_DIR))

from mlb_sgp._shared import TargetLine                 # noqa: E402
from mlb_sgp.tests import golden_replay as gr          # noqa: E402
from mlb_sgp.verify_books import ALL_BOOKS             # noqa: E402


def _targets(blob: dict) -> list[TargetLine]:
    return [TargetLine(game_id=t["game_id"], home_team=t["home_team"],
                       away_team=t["away_team"],
                       commence_time=datetime.fromisoformat(
                           t["commence_time"]),
                       period=t["period"], spread=t["spread"],
                       total=t["total"])
            for t in blob["targets"]]


@pytest.mark.parametrize("book", ALL_BOOKS)
def test_golden_rows_reproduce_offline(book, monkeypatch):
    """The orchestrator, replayed on frozen payloads, yields the pinned rows.

    A book with no fixture is skipped rather than failed — fixtures are
    captured from a live slate and a book may legitimately not have been
    capturable when they were taken. ``test_every_book_has_a_golden`` is what
    notices a missing one.
    """
    blob = gr.load_golden(book)
    if blob is None:
        pytest.skip(f"no golden fixture captured for {book}")

    client, fetchers = gr.replay(book, monkeypatch, blob["payloads"])
    rows = gr.price_with(book, _targets(blob), client, fetchers)
    actual = gr.sort_rows(gr.rows_to_json(rows))

    assert actual == blob["rows"], (
        f"{book}: replayed output drifted from the golden. "
        f"{len(actual)} rows now vs {len(blob['rows'])} pinned.")


@pytest.mark.parametrize("book", ALL_BOOKS)
def test_golden_replay_touches_no_network(book, monkeypatch):
    """Replay must be hermetic — a golden that silently reaches the wire is
    not a golden. Every transport the six clients use is poisoned first."""
    blob = gr.load_golden(book)
    if blob is None:
        pytest.skip(f"no golden fixture captured for {book}")

    def _boom(*_a, **_k):
        raise AssertionError(f"{book}: golden replay attempted a live call")

    import requests
    monkeypatch.setattr(requests.Session, "request", _boom, raising=False)
    monkeypatch.setattr(requests, "get", _boom, raising=False)
    monkeypatch.setattr(requests, "post", _boom, raising=False)
    try:
        from curl_cffi import requests as cffi_requests
        monkeypatch.setattr(cffi_requests.Session, "request", _boom,
                            raising=False)
    except ImportError:                       # pragma: no cover
        pass

    client, fetchers = gr.replay(book, monkeypatch, blob["payloads"])
    rows = gr.price_with(book, _targets(blob), client, fetchers)
    assert rows, f"{book}: replay produced no rows"


def test_every_book_has_a_golden():
    """The gate (#42) requires a golden per ACTIVE book.

    Listed explicitly rather than derived, so deleting a fixture is a test
    failure instead of a silently smaller matrix.
    """
    missing = [b for b in ALL_BOOKS if gr.load_golden(b) is None]
    assert not missing, (
        f"no golden fixture for: {', '.join(missing)}. Capture one with "
        f"mlb_sgp/tests/capture_goldens.py --books {','.join(missing)} "
        f"from a slate where those books are healthy.")


@pytest.mark.parametrize("book", ALL_BOOKS)
def test_golden_rows_are_self_consistent(book):
    """Pinned rows must look like real prices, not placeholders.

    Cheap arithmetic guards so a fixture captured from a half-broken book
    (all-equal decimals, non-positive prices) cannot become the baseline.
    """
    blob = gr.load_golden(book)
    if blob is None:
        pytest.skip(f"no golden fixture captured for {book}")
    rows = blob["rows"]
    assert rows, f"{book}: golden pins zero rows"
    for row in rows:
        assert row["sgp_decimal"] > 1.0, (
            f"{book}: decimal {row['sgp_decimal']} is not a payout")
        assert row["bookmaker"], f"{book}: row has no bookmaker"

    # Issue #64's failure mode, generalised: rows carrying different line
    # labels must not all share one price.
    by_line = {}
    for row in rows:
        by_line.setdefault((row["spread_line"], row["total_line"]),
                           set()).add(row["sgp_decimal"])
    if len(by_line) > 1:
        vectors = {frozenset(v) for v in by_line.values()}
        assert len(vectors) > 1, (
            f"{book}: every line carries an identical price vector — the "
            f"orchestrator is pricing one leg set under many labels")
