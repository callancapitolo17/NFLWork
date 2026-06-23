"""Per-bot adapter config for the Kalshi MLB monitor.

The two bots store the same *concepts* under different table/column names. Rather
than write every query twice, we describe each bot once here and let
``queries.py`` build SQL from these maps. This is the only place that knows about
schema differences between the maker and the taker.

Paths point at the *canonical live* bot DBs under ``~/NFLWork`` (override with
``KALSHI_MLB_MONITOR_ROOT``). The monitor always watches the running bots, so it
reads the real DBs even when this code is run from a git worktree.
"""

from __future__ import annotations

import os

# Repo root that holds the two bot directories. Defaults to the user's checkout
# so the dashboard reads live data regardless of where this file is executed
# from (e.g. a worktree). Override for tests / alternate layouts.
ROOT = os.environ.get(
    "KALSHI_MLB_MONITOR_ROOT",
    os.path.expanduser("~/NFLWork"),
)


def _p(*parts: str) -> str:
    return os.path.join(ROOT, *parts)


# Each bot adapter is a plain dict. ``kind`` drives the handful of places where
# the maker and taker genuinely diverge (funnel stages, fill columns); the rest
# of the fields are looked up generically.
BOTS = {
    "maker": {
        "kind": "maker",
        "label": "Maker",
        "subtitle": "kalshi_mlb_mm · provides liquidity by quoting others' RFQs",
        "proc_match": "kalshi_mlb_mm.main",  # pgrep pattern to detect if live
        "state_db": _p("kalshi_mlb_mm", "kalshi_mlb_mm.duckdb"),
        "research_db": _p("kalshi_mlb_mm", "kalshi_mlb_mm_research.duckdb"),
        # RFQ ledger: one row per RFQ the maker *observed*.
        "rfq_table": "seen_rfqs",
        "rfq_ts": "first_seen_at",
        # Decision ledger: one row per quote-eval decision.
        "decision_table": "quote_decisions",
        "decision_col": "decision",
        "reason_col": "reason",
        "decision_ts": "observed_at",
        # Fills column map (logical -> real).
        "fills_table": "fills",
        "fill_side": "side_held",
        "fill_price": "price",
        "fill_fee": "fee",
        "fill_ts": "filled_at",
        "fill_fair_quote": "blended_fair_at_quote",
        "fill_fair_confirm": "fair_at_confirm",
        "fill_pnl": "realized_pnl",
        # Open working orders the maker has resting.
        "open_table": "live_quotes",
        "open_ts": "submitted_at",
        # Positions.
        "positions_table": "positions",
        "has_position_legs": False,
    },
    "taker": {
        "kind": "taker",
        "label": "Taker",
        "subtitle": "kalshi_mlb_rfq · sends RFQs to take liquidity",
        "proc_match": "kalshi_mlb_rfq.main",
        "state_db": _p("kalshi_mlb_rfq", "kalshi_mlb_rfq.duckdb"),
        "research_db": _p("kalshi_mlb_rfq", "kalshi_mlb_rfq_research.duckdb"),
        # RFQ ledger: one row per RFQ the taker *sent*.
        "rfq_table": "live_rfqs",
        "rfq_ts": "submitted_at",
        # Decision ledger: one row per quote received & evaluated.
        "decision_table": "quote_log",
        "decision_col": "decision",
        "reason_col": "reason_detail",
        "decision_ts": "observed_at",
        # Fills column map.
        "fills_table": "fills",
        "fill_side": "side",
        "fill_price": "price_dollars",
        "fill_fee": "fee_dollars",
        "fill_ts": "filled_at",
        "fill_fair_quote": "blended_fair_at_fill",
        "fill_fair_confirm": None,  # taker has no confirm-time re-price
        "fill_pnl": None,           # taker tracks no per-fill realized pnl column
        # Open working orders: RFQs still resting on the book.
        "open_table": "live_rfqs",
        "open_ts": "submitted_at",
        # Positions.
        "positions_table": "positions",
        "has_position_legs": True,
    },
}

# Plain-language gloss for the decision/reason codes, so the dashboard can show a
# human legend next to the raw strings. Unknown codes simply fall through to "".
REASON_GLOSS = {
    # ---- shared / maker ----
    "out_of_scope": "Not a spread×total 2-leg combo — outside what the bot quotes",
    "no_fair": "No blended fair (too few fresh books, or fair out of bounds)",
    "size_gate": "Worst-case exposure exceeds the per-fill contract cap",
    "size_gate_dollars": "Worst-case exposure exceeds the per-fill dollar cap",
    "size_unknown": "RFQ had no size and no target cost — can't suggest a quote",
    "no_game": "Couldn't resolve a game from the market ticker",
    "no_mapping": "A leg failed to type-check against the MLB event mapping",
    "daily_cap": "Today's total exposure already at the daily cap",
    "per_game_cap": "Today's exposure on this game already at the per-game cap",
    "per_combo_cap": "Concentration on this combo already at the per-combo cap",
    "in_cooldown": "Combo is in post-fill cooldown",
    "tipoff": "Game starts inside the tipoff blackout window",
    "unpriceable": "Quote pricing failed",
    "skipped_creator_halt": "Per-creator fill halt active",
    "circuit_breaker": "Book moved past the circuit-breaker threshold — quotes pulled",
    "quoted": "Quote submitted",
    "skipped": "Did not quote this RFQ",
    # ---- taker decisions ----
    "accepted": "Filled — quote accepted",
    "declined_ev": "Quote's post-fee EV below the minimum threshold",
    "declined_stale_predictions": "Model samples too old to trust",
    "halted_low_fill_ratio": "Fill ratio dropped — trading halted to investigate",
    "failed_quote_walked": "Tried to accept but the quote walked (race/expired)",
    "declined_sanity": "Quote failed a sanity check",
    "declined_tipoff": "Inside the tipoff blackout window",
    "declined_inverse_lock": "Opposite side already held — avoid a hedge trap",
    "declined_kelly_zero": "Kelly sizing returned 0 contracts",
    "declined_side_mismatch": "Chosen side didn't match the RFQ's intended side",
    "declined_per_game_cap": "Per-game exposure cap hit",
    "declined_daily_cap": "Daily exposure cap hit",
    "declined_dry_run": "Would have accepted, but running in dry-run",
    "declined_positions_unhealthy": "Positions API unhealthy — trading paused",
    "declined_killswitch": "Manual kill-switch file present",
    "declined_cooldown": "Combo+side in post-fill cooldown",
    "declined_math_invariant": "Both sides looked +EV (impossible) — guard tripped",
    # ---- maker void/last-look ----
    "voided_last_look": "Accepted then voided — fair drifted past tolerance",
    "voided_no_legs": "Accepted then voided — leg data missing on confirm",
    "voided_no_fresh_books": "Accepted then voided — no fresh books to re-price",
    "voided_blend_failed": "Accepted then voided — blend recompute failed",
    "halted_high_void_rate": "Void rate too high — quoting halted",
}


def gloss(code: str) -> str:
    """Human-readable description for a decision/reason code (or '' if unknown)."""
    if code is None:
        return ""
    return REASON_GLOSS.get(code, "")
