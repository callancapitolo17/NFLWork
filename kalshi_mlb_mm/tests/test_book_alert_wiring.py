"""Both bots are actually wired to the #37 alerter — and to the RIGHT floor.

The alerter itself is covered in kalshi_common/tests/test_book_health_alerts.py.
What can silently rot instead is the wiring: a main_loop refactor that drops
the per-tick dispatch, or a Rule B floor that stops matching the gate it is
supposed to be warning about. Neither shows up as a test failure anywhere
else, and both fail exactly the way this whole ticket exists to prevent —
silently.
"""
from __future__ import annotations

import re
from pathlib import Path

from kalshi_common.book_health import BookHealthAlerter, build_alerter
from kalshi_mlb_mm import config as mm_config
from kalshi_mlb_rfq import config as rfq_config

REPO = Path(__file__).resolve().parents[2]


def _source(rel: str) -> str:
    return (REPO / rel).read_text()


# ------------------------------------------------------------------ #
# Config defaults                                                     #
# ------------------------------------------------------------------ #

def test_maker_book_alert_knobs_have_defaults():
    assert mm_config.BOOK_ALERT_ENABLED is True
    assert mm_config.BOOK_ALERT_STREAK == 3
    assert mm_config.BOOK_ALERT_PATHS == ("on_demand",)


def test_paths_knob_parses_a_single_path():
    """#57 sets BOOK_ALERT_PATHS=on_demand. If that parsed to a string, the
    `path not in self.paths` check would match SUBSTRINGS and quietly keep
    alerting on the sweep."""
    raw = "on_demand"
    parsed = tuple(p.strip() for p in raw.split(",") if p.strip())
    assert parsed == ("on_demand",)
    a = BookHealthAlerter(label="MM", paths=parsed, notify_async=False)
    a.observe(book="dk", path="sweep", outcome="transport_error")
    assert a.snapshot()["streaks"] == {}


# ------------------------------------------------------------------ #
# build_alerter                                                       #
# ------------------------------------------------------------------ #

def test_build_alerter_disabled_returns_none():
    """Disabled must mean NO OBJECT, so SGPService takes its inert path
    rather than carrying a live-looking alerter that never fires."""
    assert build_alerter(label="MM", enabled=False, streak_threshold=3,
                         min_healthy_books=2, paths=("sweep",)) is None


def test_build_alerter_carries_config_through():
    a = build_alerter(label="RFQ", enabled=True, streak_threshold=5,
                      min_healthy_books=3, paths=("on_demand",))
    assert (a.label, a.streak_threshold, a.min_healthy_books, a.paths) == (
        "RFQ", 5, 3, ("on_demand",))


# ------------------------------------------------------------------ #
# Tick-loop wiring                                                    #
# ------------------------------------------------------------------ #

def test_both_tick_loops_dispatch_alerts_every_tick():
    """An alerter nobody dispatches is a silent failure detector — the exact
    shape of bug #37 exists to remove."""
    for rel in ("kalshi_mlb_mm/main.py", "kalshi_mlb_rfq/main.py"):
        src = _source(rel)
        assert "sgp_service.check_book_health()" in src, rel
        # Must sit next to the health flush, i.e. inside the tick loop and
        # not in some startup-only branch.
        assert re.search(r"sgp_service\.flush_health\(\)[^\n]*\n\s*"
                         r"sgp_service\.check_book_health\(\)", src), rel


def test_each_bot_alerts_at_its_own_consensus_floor():
    """Rule B must fire at the count the bot's REAL gate uses. A separate
    knob could drift; these two assertions are what keep it honest."""
    assert "min_healthy_books=config.MIN_AGREEING_BOOKS" in _source(
        "kalshi_mlb_mm/main.py")
    assert "min_healthy_books=config.MIN_BOOK_COUNT_FOR_BLEND" in _source(
        "kalshi_mlb_rfq/main.py")
    # And the two gates are the values we think they are.
    assert mm_config.MIN_AGREEING_BOOKS == 2
    assert rfq_config.MIN_BOOK_COUNT_FOR_BLEND == 2
