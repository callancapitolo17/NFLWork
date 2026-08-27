"""Cached single-leg book fairs — the maker's leg surface (epic #94, #96).

Two decoupled loops. The INGEST loop (this package) refreshes single-leg book
prices on our schedule, so its work is O(games) and burst-immune. The QUOTE
loop (issue #98) answers a CROSS-GAME RFQ by reading the in-memory store here
and multiplying per-game fairs — zero network I/O. Same-game combos are NOT in
scope and keep the live on-demand path unchanged.

Public surface:
    SurfaceKey / SurfaceRow  — one leg's devigged fair at one book
    LegSurface               — the thread-safe in-memory store #98 reads
    SurfaceIngest            — the ingest loop (runner.py)

Side effects: live HTTP to the books, and writes to
``kalshi_mlb_mm/kalshi_mlb_mm_surface.duckdb`` (its own file and write lock —
never the market DB the pricing path reads).
"""
from kalshi_mlb_mm.leg_surface.store import (LegSurface, SurfaceKey,
                                             SurfaceRow, surface_key)

__all__ = ["LegSurface", "SurfaceKey", "SurfaceRow", "surface_key"]
