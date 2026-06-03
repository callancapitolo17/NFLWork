"""Fair-value assembly for the maker — pure-book consensus blend.

v1 hardening: the model (kalshi_common.fair_value.model_fair against
mlb_game_samples) was removed from the blend because it was being medianed-out
in the median(model, books) combination (so it barely moved the answer in
typical configurations), carried documented bias on certain combo families
(`mlb_parlay_edge_overestimation`), and created a soft dependency on the R
answer-key pipeline writing fresh samples. The blend is now simply the median
of the input book devigged fairs (which the caller has already filtered through
the book-consensus-band gate in `main._consensus_filter`).
"""
import statistics


def blended_fair(legs: list[dict], game_id: str,
                 book_fairs: dict[str, float]) -> tuple[float | None, float | None]:
    """Returns (book_fair_median, blended) — each None if no books supplied.

    Both values are the median of `book_fairs.values()`. The caller is expected
    to have already applied any consensus gating (outlier rejection,
    minimum-book-count) upstream; this function takes the surviving books at
    face value.
    """
    if not book_fairs:
        return None, None
    med = statistics.median(list(book_fairs.values()))
    return med, med
