#!/usr/bin/env python3
"""
Acquire CBB game scores using CBBpy - Simple sequential version.
Run overnight: nohup python -u acquire_cbb_pbp.py --all > acquire.log 2>&1 &
"""

import cbbpy.mens_scraper as scraper
import pandas as pd
import duckdb
from datetime import datetime, timedelta
import time
import sys
import argparse

DB_PATH = "../cbb.duckdb"
TABLE_NAME = "cbb_pbp_v2"


def extract_game_scores(game_id, game_date):
    """Extract scores from PBP data"""
    try:
        pbp = scraper.get_game_pbp(game_id)
        if pbp is None or len(pbp) == 0:
            return None

        home = pbp['home_team'].iloc[0]
        away = pbp['away_team'].iloc[0]

        h1 = pbp[pbp['half'] == 1]
        h2 = pbp[pbp['half'] == 2]

        if len(h1) == 0 or len(h2) == 0:
            return None

        h1_home = int(h1['home_score'].max())
        h1_away = int(h1['away_score'].max())
        h2_end_home = int(h2['home_score'].max())
        h2_end_away = int(h2['away_score'].max())

        final_home = int(pbp['home_score'].max())
        final_away = int(pbp['away_score'].max())

        ot = pbp[pbp['half'] > 2]
        ot_home = final_home - h2_end_home if len(ot) > 0 else 0
        ot_away = final_away - h2_end_away if len(ot) > 0 else 0

        return {
            'game_id': game_id,
            'game_date': game_date,
            'home_team': home,
            'away_team': away,
            'home_h1_score': h1_home,
            'away_h1_score': h1_away,
            'game_home_margin_h1': h1_home - h1_away,
            'game_total_h1': h1_home + h1_away,
            'home_h2_score': h2_end_home - h1_home,
            'away_h2_score': h2_end_away - h1_away,
            'game_home_margin_h2': (h2_end_home - h1_home) - (h2_end_away - h1_away),
            'game_total_h2': (h2_end_home - h1_home) + (h2_end_away - h1_away),
            'home_ot_score': ot_home,
            'away_ot_score': ot_away,
            'game_home_margin_ot': ot_home - ot_away,
            'game_total_ot': ot_home + ot_away,
            'home_final_score': final_home,
            'away_final_score': final_away,
            'game_home_margin_fg': final_home - final_away,
            'game_total_fg': final_home + final_away,
            'home_winner': 1 if final_home > final_away else 0,
            'went_to_ot': 1 if len(ot) > 0 else 0
        }
    except:
        return None


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--season', type=int)
    parser.add_argument('--all', action='store_true')
    parser.add_argument('--db', type=str, help='Custom database path')
    args = parser.parse_args()

    if not args.season and not args.all:
        parser.print_help()
        sys.exit(1)

    seasons = [args.season] if args.season else [2021, 2022, 2023, 2024, 2025, 2026]
    db_path = args.db if args.db else DB_PATH

    print(f"CBB PBP Acquisition - Season {args.season if args.season else 'ALL'}")
    print("=" * 50)

    con = duckdb.connect(db_path)
    con.execute(f"""
        CREATE TABLE IF NOT EXISTS {TABLE_NAME} (
            game_id VARCHAR PRIMARY KEY,
            game_date VARCHAR,
            home_team VARCHAR,
            away_team VARCHAR,
            home_h1_score INTEGER,
            away_h1_score INTEGER,
            game_home_margin_h1 INTEGER,
            game_total_h1 INTEGER,
            home_h2_score INTEGER,
            away_h2_score INTEGER,
            game_home_margin_h2 INTEGER,
            game_total_h2 INTEGER,
            home_ot_score INTEGER,
            away_ot_score INTEGER,
            game_home_margin_ot INTEGER,
            game_total_ot INTEGER,
            home_final_score INTEGER,
            away_final_score INTEGER,
            game_home_margin_fg INTEGER,
            game_total_fg INTEGER,
            home_winner INTEGER,
            went_to_ot INTEGER
        )
    """)

    try:
        existing = set(con.execute(f"SELECT game_id FROM {TABLE_NAME}").fetchdf()['game_id'].tolist())
    except:
        existing = set()

    print(f"Existing: {len(existing)} games")

    for season in seasons:
        start_date = datetime(season - 1, 11, 1)
        end_date = datetime(season, 4, 15)
        current = start_date

        season_games = []
        season_count = 0

        print(f"\n{'='*50}")
        print(f"Season {season-1}-{str(season)[2:]}")
        print(f"{'='*50}", flush=True)

        while current <= end_date:
            date_str = current.strftime("%Y-%m-%d")

            try:
                game_ids = scraper.get_game_ids(date_str)
                if game_ids:
                    new_ids = [g for g in game_ids if g not in existing]

                    for i, gid in enumerate(new_ids):
                        result = extract_game_scores(gid, date_str)
                        if result:
                            season_games.append(result)
                            existing.add(gid)
                            season_count += 1
                            print(f"    {gid}: OK", flush=True)
                        else:
                            print(f"    {gid}: skip", flush=True)
                        time.sleep(0.05)

                    # Save every 500 games
                    if len(season_games) >= 500:
                        df = pd.DataFrame(season_games)
                        con.execute(f"INSERT INTO {TABLE_NAME} SELECT * FROM df")
                        print(f"  Checkpoint: {len(season_games)} games saved (total: {season_count})", flush=True)
                        season_games = []

                    # Progress every 50 games
                    if season_count > 0 and season_count % 50 == 0:
                        print(f"  {date_str}: {season_count} games", flush=True)

            except Exception as e:
                print(f"  {date_str}: Error - {str(e)[:40]}", flush=True)

            current += timedelta(days=1)

        # Save remaining
        if season_games:
            df = pd.DataFrame(season_games)
            con.execute(f"INSERT INTO {TABLE_NAME} SELECT * FROM df")

        print(f"Season complete: {season_count} games")

    total = con.execute(f"SELECT COUNT(*) FROM {TABLE_NAME}").fetchone()[0]
    print(f"\n{'='*50}")
    print(f"Done! Total: {total} games")
    con.close()


if __name__ == "__main__":
    main()
