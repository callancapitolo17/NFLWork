"""One-shot seeder: dict literals in config/*.py -> DuckDB lookup tables.

Idempotent: TRUNCATE-and-INSERT inside a single transaction, so re-runs
produce identical state. Called automatically at the start of every
run.py invocation.
"""

from nfl_draft.lib.db import write_connection
from nfl_draft.config.players import PLAYERS
from nfl_draft.config.teams import TEAMS
from nfl_draft.config.markets import MARKET_MAP


def run() -> None:
    """Truncate + repopulate all 5 lookup tables in one transaction."""
    with write_connection() as con:
        con.execute("BEGIN TRANSACTION")
        try:
            # players
            con.execute("DELETE FROM players")
            con.execute("DELETE FROM player_aliases")
            for canonical, info in PLAYERS.items():
                con.execute(
                    "INSERT INTO players VALUES (?, ?, ?)",
                    [canonical, info["position"], info.get("college")],
                )
                for alias in info.get("aliases", []):
                    con.execute(
                        "INSERT INTO player_aliases VALUES (?, ?)",
                        [alias.lower().strip(), canonical],
                    )
            # teams
            con.execute("DELETE FROM teams")
            con.execute("DELETE FROM team_aliases")
            for abbr, info in TEAMS.items():
                con.execute(
                    "INSERT INTO teams VALUES (?, ?)",
                    [abbr, info["full_name"]],
                )
                for alias in info.get("aliases", []):
                    con.execute(
                        "INSERT INTO team_aliases VALUES (?, ?)",
                        [alias.lower().strip(), abbr],
                    )
            # market_map
            con.execute("DELETE FROM market_map")
            for book, label, subject, market_id in MARKET_MAP:
                con.execute(
                    "INSERT INTO market_map VALUES (?, ?, ?, ?)",
                    [book, label, subject, market_id],
                )
            con.execute("COMMIT")
        except Exception:
            con.execute("ROLLBACK")
            raise


if __name__ == "__main__":
    run()
    print("Seed complete.")
