"""Tests for _recon_util helpers that must be worktree-safe."""

from nfl_draft.scrapers._recon_util import _main_repo_root


def test_main_repo_root_resolves_to_main_working_tree():
    """_main_repo_root must resolve to the main repo even from a worktree.

    Heuristic: bet_logger directory should exist under the returned root
    (it lives in the main repo, not a worktree checkout).
    """
    root = _main_repo_root()
    # In our repo, bet_logger/ always exists at the main repo root
    assert (root / "bet_logger").is_dir(), f"bet_logger not found under {root}"
