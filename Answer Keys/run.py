#!/usr/bin/env python3
"""
Sports betting orchestrator
Runs scrapers AND R sample generation in parallel, then combines
"""
import asyncio
import json
import os
import subprocess
import sys
from pathlib import Path

# Force unbuffered output so prints appear in correct order with subprocess output
sys.stdout.reconfigure(line_buffering=True)

# Scraper configurations - add new books here
# Paths relative to project root (NFLWork/)
SCRAPER_CONFIGS = {
    "wagerzon": {
        "script": "../wagerzon_odds/scraper_v2.py",
        "sports": ["nfl", "nba", "cbb"]
    },
    "hoop88": {
        "script": "../hoop88_odds/scraper.py",
        "sports": ["nfl", "ncaaf", "cbb"]
    },
    "bfa": {
        "script": "../bfa_odds/scraper.py",
        "sports": ["nfl", "nba", "cbb"]
    },
}

# Sport key mapping for Odds API
SPORT_KEYS = {
    "nfl": "americanfootball_nfl",
    "ncaaf": "americanfootball_ncaaf",
    "cbb": "basketball_ncaab",
    "nba": "basketball_nba",
}

# R script configurations - add new sports here
# Paths relative to Answer Keys/
R_SCRIPTS = {
    "nfl": {
        "prepare": "NFL Answer Key/NFLPrepare.R",
        "combine": "NFL Answer Key/NFLCombine.R",
    },
    "cbb": {
        "prepare": "CBB Answer Key/CBBPrepare.R",
        "combine": "CBB Answer Key/CBBCombine.R",
    },
    # "nba": {
    #     "prepare": "NBA Answer Key/NBAPrepare.R",
    #     "combine": "NBA Answer Key/NBACombine.R",
    # },
}


async def run_scraper(name: str, script: str, sport: str) -> dict:
    """Run a single scraper asynchronously."""
    print(f"  [scraper] Starting {name}...", flush=True)

    script_path = Path(__file__).parent / script
    script_dir = script_path.parent

    # Use venv if it exists in the scraper directory
    venv_python = script_dir / "venv" / "bin" / "python"
    if venv_python.exists():
        python_exe = str(venv_python)
    else:
        python_exe = sys.executable

    proc = await asyncio.create_subprocess_exec(
        python_exe, str(script_path), sport,
        stdout=asyncio.subprocess.PIPE,
        stderr=asyncio.subprocess.PIPE,
        cwd=str(script_dir)
    )
    stdout, stderr = await proc.communicate()

    success = proc.returncode == 0
    status = "✓" if success else "✗"
    print(f"  [scraper] {status} {name} complete", flush=True)

    if not success:
        print(f"  [scraper] Error: {stderr.decode()[-2000:]}", flush=True)

    return {"name": name, "success": success, "type": "scraper"}


async def run_r_prepare(sport: str) -> dict:
    """Run R Part 1: load data, generate samples, save to DuckDB."""
    print(f"  [R] Starting sample generation...", flush=True)

    script_path = Path(__file__).parent / R_SCRIPTS[sport]["prepare"]

    proc = await asyncio.create_subprocess_exec(
        "Rscript", str(script_path),
        stdout=asyncio.subprocess.PIPE,
        stderr=asyncio.subprocess.PIPE,
        cwd=str(script_path.parent)
    )
    stdout, stderr = await proc.communicate()

    success = proc.returncode == 0
    status = "✓" if success else "✗"
    print(f"  [R] {status} Sample generation complete", flush=True)

    if not success:
        print(f"  [R] Error: {stderr.decode()[:500]}", flush=True)

    return {"name": "r_prepare", "success": success, "type": "r"}


async def run_parallel_phase(sport: str) -> list:
    """Run scrapers AND R preparation in parallel."""
    tasks = []

    # Add all applicable scrapers
    for name, config in SCRAPER_CONFIGS.items():
        if sport in config["sports"]:
            tasks.append(run_scraper(name, config["script"], sport))

    # Add R sample generation
    if sport in R_SCRIPTS:
        tasks.append(run_r_prepare(sport))
    else:
        print(f"Warning: No R scripts configured for {sport}", flush=True)

    if not tasks:
        print(f"No tasks configured for {sport}", flush=True)
        return []

    print(f"Running {len(tasks)} tasks in parallel...", flush=True)
    results = await asyncio.gather(*tasks, return_exceptions=True)
    return results


def run_r_combine(sport: str) -> int:
    """Run R Part 2: load samples + scraped odds, find edge."""
    script_path = Path(__file__).parent / R_SCRIPTS[sport]["combine"]
    print(f"\n[R] Combining odds and finding edge...", flush=True)

    result = subprocess.run(
        ["Rscript", str(script_path)],
        cwd=str(script_path.parent)
    )
    return result.returncode


def fetch_canonical_games(sport: str):
    """Fetch today's Odds API games and save as JSON for scrapers to use."""
    try:
        import requests
    except ImportError:
        print("Warning: 'requests' not installed, skipping canonical games fetch", flush=True)
        return

    api_key = os.environ.get("ODDS_API_KEY")
    sport_key = SPORT_KEYS.get(sport)
    if not api_key or not sport_key:
        return

    try:
        resp = requests.get(
            f"https://api.the-odds-api.com/v4/sports/{sport_key}/odds",
            params={"apiKey": api_key, "regions": "us", "markets": "spreads",
                    "oddsFormat": "american", "dateFormat": "iso"},
            timeout=30,
        )
        if resp.status_code != 200:
            print(f"Warning: Could not fetch canonical games (HTTP {resp.status_code})", flush=True)
            return

        games = [{"home_team": g["home_team"], "away_team": g["away_team"]}
                 for g in resp.json()]

        out_path = Path(__file__).parent / f".canonical_games_{sport}.json"
        with open(out_path, "w") as f:
            json.dump(games, f)
        print(f"Fetched {len(games)} canonical games for {sport}", flush=True)
    except Exception as e:
        print(f"Warning: Could not fetch canonical games: {e}", flush=True)


def cleanup_canonical_games(sport: str):
    """Remove temporary canonical game files."""
    path = Path(__file__).parent / f".canonical_games_{sport}.json"
    if path.exists():
        path.unlink()


def main():
    sport = sys.argv[1] if len(sys.argv) > 1 else "nfl"
    print(f"=== {sport.upper()} Answer Key ===\n", flush=True)

    if sport not in R_SCRIPTS:
        print(f"Error: Sport '{sport}' not configured. Available: {list(R_SCRIPTS.keys())}", flush=True)
        return 1

    # Phase 0: Fetch canonical game list for team name resolution
    fetch_canonical_games(sport)

    # Phase 1: Run scrapers + R samples in parallel
    results = asyncio.run(run_parallel_phase(sport))

    # Check for failures
    failures = [r for r in results if isinstance(r, dict) and not r.get("success")]
    if failures:
        print(f"\nWarning: {len(failures)} tasks failed:", flush=True)
        for f in failures:
            print(f"  - {f.get('name', 'unknown')}", flush=True)

    # Phase 2: Combine and find edge
    rc = run_r_combine(sport)

    # Cleanup temp files
    cleanup_canonical_games(sport)

    return rc


if __name__ == "__main__":
    sys.exit(main())
