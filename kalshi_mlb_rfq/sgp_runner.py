"""Backward-compat shim. Logic lives in kalshi_common.sgp_runner."""
from kalshi_common.sgp_runner import *  # noqa: F401,F403
from kalshi_common.sgp_runner import (
    should_scrape, latest_sgp_fetch_time, enumerate_kalshi_targets,
    write_target_lines, run_scrapers, read_priced_rows, sgp_cycle,
    SCRAPER_NAMES,
)
