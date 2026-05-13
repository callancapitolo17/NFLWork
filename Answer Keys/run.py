#!/usr/bin/env python3
"""
Sports betting orchestrator
Runs scrapers AND R answer key in parallel, with sentinel signaling
"""
import asyncio
import json
import os
import sys
import time
from pathlib import Path

# Force unbuffered output so prints appear in correct order with subprocess output
sys.stdout.reconfigure(line_buffering=True)

# Scraper configurations - add new books here
# Paths relative to project root (NFLWork/)
SCRAPER_CONFIGS = {
    "wagerzon": {
        "script": "../wagerzon_odds/scraper_v2.py",
        "sports": ["nfl", "nba", "cbb", "college_baseball", "mlb"]
    },
    "wagerzon_specials": {
        "script": "../wagerzon_odds/scraper_specials.py",
        "sports": ["mlb"]
    },
    "hoop88": {
        "script": "../hoop88_odds/scraper.py",
        "sports": ["nfl", "ncaaf", "cbb", "college_baseball", "mlb"]
    },
    "bfa": {
        "script": "../bfa_odds/scraper.py",
        "sports": ["nfl", "nba", "cbb", "mlb"]
    },
    "bookmaker": {
        "script": "../bookmaker_odds/scraper.py",
        "sports": ["cbb", "mlb"]
    },
    "bet105": {
        "script": "../bet105_odds/scraper.py",
        "sports": ["cbb", "mlb"]
    },
    "kalshi": {
        "script": "../kalshi_odds/scraper.py",
        "sports": ["cbb", "mlb"]
    },
    "draftkings_singles": {
        "script": "../mlb_sgp/scraper_draftkings_singles.py",
        "sports": ["mlb"]
    },
    "fanduel_singles": {
        "script": "../mlb_sgp/scraper_fanduel_singles.py",
        "sports": ["mlb"]
    },
}

# Sport key mapping for Odds API
SPORT_KEYS = {
    "nfl": "americanfootball_nfl",
    "ncaaf": "americanfootball_ncaaf",
    "cbb": "basketball_ncaab",
    "nba": "basketball_nba",
    "college_baseball": "baseball_ncaa",
    "mlb": "baseball_mlb",
}

# R script configurations - add new sports here
# Paths relative to Answer Keys/
R_SCRIPTS = {
    "nfl": {
        "prepare": "NFL Answer Key/NFLPrepare.R",
        "combine": "NFL Answer Key/NFLCombine.R",
    },
    "cbb": {
        "script": "CBB Answer Key/CBB.R",
    },
    "college_baseball": {
        "script": "College Baseball Answer Key/CollegeBaseball.R",
    },
    "mlb": {
        "script": "MLB Answer Key/MLB.R",
    },
    # "nba": {
    #     "script": "NBA Answer Key/NBA.R",
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

    t0 = time.time()
    proc = await asyncio.create_subprocess_exec(
        python_exe, str(script_path), sport,
        stdout=asyncio.subprocess.PIPE,
        stderr=asyncio.subprocess.PIPE,
        cwd=str(script_dir)
    )
    stdout, stderr = await proc.communicate()
    elapsed = time.time() - t0

    success = proc.returncode == 0
    status = "✓" if success else "✗"
    print(f"  [scraper] {status} {name} complete ({elapsed:.1f}s)", flush=True)

    if not success:
        print(f"  [scraper] Error: {stderr.decode()[-2000:]}", flush=True)

    return {"name": name, "success": success, "type": "scraper", "elapsed": elapsed}


async def run_r(sport: str, script_path: Path) -> dict:
    """Run R answer key script asynchronously."""
    print(f"  [R] Starting answer key...", flush=True)

    t0 = time.time()
    proc = await asyncio.create_subprocess_exec(
        "Rscript", str(script_path),
        stdout=asyncio.subprocess.PIPE,
        stderr=asyncio.subprocess.PIPE,
        cwd=str(script_path.parent)
    )
    stdout, stderr = await proc.communicate()
    elapsed = time.time() - t0

    success = proc.returncode == 0
    status = "✓" if success else "✗"
    print(f"  [R] {status} Answer key complete ({elapsed:.1f}s)", flush=True)

    # Print R output (contains edge summary)
    r_output = stdout.decode()
    if r_output.strip():
        print(r_output, flush=True)

    if not success:
        print(f"  [R] Error: {stderr.decode()[:500]}", flush=True)

    return {"name": "r_answer_key", "success": success, "type": "r", "elapsed": elapsed}


async def run_r_two_phase(sport: str) -> dict:
    """Run R prepare + combine as two separate processes (legacy sports)."""
    config = R_SCRIPTS[sport]
    prepare_path = Path(__file__).parent / config["prepare"]
    combine_path = Path(__file__).parent / config["combine"]

    # Phase 1: prepare
    print(f"  [R] Starting sample generation...", flush=True)
    t0 = time.time()
    proc = await asyncio.create_subprocess_exec(
        "Rscript", str(prepare_path),
        stdout=asyncio.subprocess.PIPE,
        stderr=asyncio.subprocess.PIPE,
        cwd=str(prepare_path.parent)
    )
    stdout, stderr = await proc.communicate()
    elapsed = time.time() - t0

    success = proc.returncode == 0
    status = "✓" if success else "✗"
    print(f"  [R] {status} Sample generation complete ({elapsed:.1f}s)", flush=True)

    if not success:
        print(f"  [R] Error: {stderr.decode()[:500]}", flush=True)

    return {"name": "r_prepare", "success": success, "type": "r", "elapsed": elapsed}


async def run_pipeline(sport: str) -> tuple[list, int]:
    """Run scrapers + R in parallel with sentinel signaling."""
    config = R_SCRIPTS[sport]
    is_merged = "script" in config  # merged single-script sports

    # Sharp scrapers (bookmaker, bet105) must finish BEFORE R starts Phase 2
    # so their data is in DuckDB when CBB.R reads it for consensus.
    SHARP_SCRAPERS = {"bookmaker", "bet105"}

    sharp_tasks = []
    other_tasks = []
    for name, scraper_config in SCRAPER_CONFIGS.items():
        if sport in scraper_config["sports"]:
            task = asyncio.create_task(
                run_scraper(name, scraper_config["script"], sport)
            )
            if name in SHARP_SCRAPERS:
                sharp_tasks.append(task)
            else:
                other_tasks.append(task)

    # Run sharp scrapers first (in parallel with each other)
    if sharp_tasks:
        print(f"Running {len(sharp_tasks)} sharp scrapers first...", flush=True)
        sharp_results = await asyncio.gather(*sharp_tasks, return_exceptions=True)
    else:
        sharp_results = []

    # Now start R (which reads sharp scraper DBs in Phase 2) + remaining scrapers
    if is_merged:
        script_path = Path(__file__).parent / config["script"]
        r_task = asyncio.create_task(run_r(sport, script_path))
    else:
        r_task = asyncio.create_task(run_r_two_phase(sport))

    n_tasks = len(sharp_tasks) + len(other_tasks) + 1
    print(f"Running {n_tasks} tasks total ({len(other_tasks)} scrapers + R in parallel)...", flush=True)

    # Wait for remaining scrapers
    other_results = await asyncio.gather(*other_tasks, return_exceptions=True)
    scraper_results = list(sharp_results) + list(other_results)

    if is_merged:
        # Write sentinel file so merged R script knows scraper data is ready
        sentinel = Path(__file__).parent / f".scrapers_done_{sport}"
        sentinel.touch()

    # Wait for R to finish (it detects sentinel and does combine phase)
    r_result = await r_task

    # Cleanup sentinel
    if is_merged:
        sentinel.unlink(missing_ok=True)

    results = list(scraper_results) + [r_result]

    # For legacy two-phase sports, run combine after everything
    rc = 0
    if not is_merged and r_result.get("success"):
        combine_path = Path(__file__).parent / config["combine"]
        print(f"\n[R] Combining odds and finding edge...", flush=True)
        t0 = time.time()
        import subprocess
        result = subprocess.run(
            ["Rscript", str(combine_path)],
            cwd=str(combine_path.parent)
        )
        combine_elapsed = time.time() - t0
        rc = result.returncode
        results.append({
            "name": "r_combine", "success": rc == 0,
            "type": "r", "elapsed": combine_elapsed
        })

    # Check R exit code for merged scripts
    if is_merged and not r_result.get("success"):
        rc = 1

    return results, rc


def fetch_canonical_games(sport: str):
    """Fetch today's Odds API games and save as JSON for scrapers to use."""
    try:
        import requests
    except ImportError:
        print("Warning: 'requests' not installed, skipping canonical games fetch", flush=True)
        return

    api_key = os.environ.get("ODDS_API_KEY")
    # Fallback: read from ~/.Renviron (where R stores it)
    if not api_key:
        renviron = Path.home() / ".Renviron"
        if renviron.exists():
            for line in renviron.read_text().splitlines():
                if line.strip().startswith("ODDS_API_KEY"):
                    api_key = line.split("=", 1)[1].strip()
                    break
    sport_key = SPORT_KEYS.get(sport)
    if not api_key or not sport_key:
        return

    out_path = Path(__file__).parent / f".canonical_games_{sport}.json"

    # Reuse cached canonical games if <2 hours old (saves Odds API credits)
    if out_path.exists():
        age_sec = time.time() - out_path.stat().st_mtime
        if age_sec < 7200:
            with open(out_path) as f:
                games = json.load(f)
            print(f"Using cached canonical games for {sport} ({len(games)} games, {age_sec/60:.0f}m old)", flush=True)
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

        with open(out_path, "w") as f:
            json.dump(games, f)
        print(f"Fetched {len(games)} canonical games for {sport}", flush=True)
    except Exception as e:
        print(f"Warning: Could not fetch canonical games: {e}", flush=True)


def cleanup_canonical_games(sport: str):
    """Keep canonical games cache for reuse (saves Odds API credits)."""
    pass


def main():
    sport = sys.argv[1] if len(sys.argv) > 1 else "nfl"
    print(f"=== {sport.upper()} Answer Key ===\n", flush=True)

    if sport not in R_SCRIPTS:
        print(f"Error: Sport '{sport}' not configured. Available: {list(R_SCRIPTS.keys())}", flush=True)
        return 1

    timings = {}
    t_total = time.time()

    # Phase 0: Fetch canonical game list for team name resolution
    t0 = time.time()
    fetch_canonical_games(sport)
    timings["phase0"] = time.time() - t0

    # Phase 1: Run scrapers + R in parallel (with sentinel signaling for merged scripts)
    t0 = time.time()
    results, rc = asyncio.run(run_pipeline(sport))
    timings["phase1"] = time.time() - t0

    # Extract per-task timings
    for r in results:
        if isinstance(r, dict) and "elapsed" in r:
            timings[r["name"]] = r["elapsed"]

    # Check for failures
    failures = [r for r in results if isinstance(r, dict) and not r.get("success")]
    if failures:
        print(f"\nWarning: {len(failures)} tasks failed:", flush=True)
        for f in failures:
            print(f"  - {f.get('name', 'unknown')}", flush=True)

    # Cleanup temp files
    cleanup_canonical_games(sport)

    timings["total"] = time.time() - t_total

    # Print timing summary
    print(f"\n=== PIPELINE TIMING ===", flush=True)
    print(f"  Phase 0 (canonical games): {timings.get('phase0', 0):5.1f}s", flush=True)
    print(f"  Phase 1 (parallel):        {timings.get('phase1', 0):5.1f}s", flush=True)
    for name in ["bookmaker", "bet105", "wagerzon", "hoop88", "bfa", "kalshi", "r_answer_key", "r_prepare", "r_combine"]:
        if name in timings:
            print(f"    - {name:20s}  {timings[name]:5.1f}s", flush=True)
    print(f"  Total:                     {timings.get('total', 0):5.1f}s", flush=True)

    return rc


if __name__ == "__main__":
    sys.exit(main())
