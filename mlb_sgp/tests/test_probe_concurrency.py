"""Unit tests for the pure planning/summary helpers in probe_concurrency."""
from datetime import datetime, timezone

from mlb_sgp._shared import TargetLine
from mlb_sgp.probe_concurrency import (
    BOOK_PLANS,
    pick_probe_targets,
    summarize_level,
)

CT = datetime(2026, 6, 8, 23, 0, tzinfo=timezone.utc)


def _t(gid, spread=-1.5, total=8.5):
    return TargetLine(game_id=gid, home_team="H", away_team="A",
                      commence_time=CT, period="FG", spread=spread, total=total)


def test_pick_probe_targets_one_per_game_capped():
    targets = [_t("g1"), _t("g1", -2.5), _t("g2"), _t("g3"), _t("g3", -3.5)]
    picked = pick_probe_targets(targets, n_games=2)
    assert len(picked) == 2
    assert [p.game_id for p in picked] == ["g1", "g2"]


def test_pick_probe_targets_prefers_main_spread():
    targets = [_t("g1", -4.5), _t("g1", -1.5), _t("g1", -2.5)]
    picked = pick_probe_targets(targets, n_games=1)
    assert picked[0].spread == -1.5


def test_summarize_level_computes_rates():
    s = summarize_level(level=8, n_targets=6, n_rows=20, elapsed_sec=4.0)
    assert s["level"] == 8
    assert s["expected_rows"] == 24
    assert s["rows_per_sec"] == 5.0
    assert abs(s["miss_pct"] - (4 / 24 * 100)) < 1e-9


def test_book_plans_respect_budgets():
    # call budget = sum over levels of (targets_per_level * 4 combos)
    # + health check (4 calls). Spec: DK <= 150, PX <= 40.
    for book, plan in BOOK_PLANS.items():
        calls = sum(plan["targets_per_level"] * 4 for _ in plan["levels"]) + 4
        assert calls <= plan["budget"], f"{book} plan exceeds budget"
