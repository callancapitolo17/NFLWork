"""Fair-value assembly for the maker — reuses kalshi_common math against the
maker's own caches. Mirrors the taker's _fresh_blended_fair, minus RFQ creation."""
import statistics

import pandas as pd

from kalshi_common import fair_value
from kalshi_common.leg_types import _leg_dict_to_typed


def blended_fair(legs: list[dict], game_id: str, samples: pd.DataFrame | None,
                 book_fairs: dict[str, float]) -> tuple[float | None, float | None, float | None]:
    """Returns (model_fair, book_fair_median, blended) — each None if unavailable.

    book_fair_median is the median of the per-book devigged fairs (for logging).
    blended is None when fewer than 2 sources (model + books) are present —
    fair_value.blend enforces that 2-source gate.
    """
    if samples is None or samples.empty:
        return None, None, None
    typed = [_leg_dict_to_typed(l, game_id) for l in legs]
    if any(t is None for t in typed):
        return None, None, None
    model = fair_value.model_fair(samples, typed)
    blended = fair_value.blend(model, book_fairs)
    book_med = statistics.median(list(book_fairs.values())) if book_fairs else None
    return model, book_med, blended
