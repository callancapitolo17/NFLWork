#!/usr/bin/env python3
"""
NFL +EV Betting Dashboard Server

Flask server that:
- Serves the dashboard HTML
- Provides API for marking bets as placed/removed
- Handles refresh to regenerate predictions

Run with: python nfl_dashboard_server.py
Then open: http://localhost:8081
"""

import subprocess
import json
import hashlib
import mimetypes
from datetime import datetime
from pathlib import Path

import duckdb
from flask import Flask, Response, jsonify, request

BASE_DIR = Path(__file__).parent.resolve()
DB_PATH = BASE_DIR / "nfl_dashboard.duckdb"
app = Flask(__name__)


# =============================================================================
# DATABASE INITIALIZATION
# =============================================================================

def init_db():
    """Create tables if they don't exist."""
    con = duckdb.connect(str(DB_PATH))

    con.execute("""
        CREATE TABLE IF NOT EXISTS placed_bets (
            bet_hash TEXT PRIMARY KEY,
            game_id TEXT NOT NULL,
            home_team TEXT NOT NULL,
            away_team TEXT NOT NULL,
            game_time TIMESTAMP NOT NULL,
            market TEXT NOT NULL,
            bet_on TEXT NOT NULL,
            line REAL,
            model_prob REAL NOT NULL,
            model_ev REAL NOT NULL,
            recommended_size REAL NOT NULL,
            actual_size REAL,
            odds INTEGER NOT NULL,
            bookmaker TEXT NOT NULL,
            placed_at TIMESTAMP DEFAULT CURRENT_TIMESTAMP,
            status TEXT DEFAULT 'pending'
        )
    """)

    con.close()
    print(f"Database initialized at: {DB_PATH}")


# =============================================================================
# ROUTES
# =============================================================================

@app.route("/")
def index():
    """Serve the dashboard HTML."""
    html_path = BASE_DIR / "report.html"
    if not html_path.exists():
        return Response(
            "<h1>Dashboard not generated yet</h1>"
            "<p>Click the refresh button or run: <code>Rscript nfl_dashboard.R</code></p>",
            mimetype="text/html"
        )
    with open(html_path, "r", encoding="utf-8") as f:
        content = f.read()
    return Response(content, mimetype="text/html")


@app.route("/lib/<path:filename>")
def serve_lib(filename):
    """Serve static library files."""
    file_path = BASE_DIR / "lib" / filename
    if not file_path.exists():
        return "File not found", 404

    mimetype, _ = mimetypes.guess_type(str(file_path))
    with open(file_path, "rb") as f:
        content = f.read()
    return Response(content, mimetype=mimetype or "application/octet-stream")


@app.route("/api/place-bet", methods=["POST"])
def place_bet():
    """Mark a bet as placed in the database."""
    data = request.json

    required = ["bet_hash", "game_id", "home_team", "away_team", "game_time",
                "market", "bet_on", "model_prob", "model_ev", "recommended_size",
                "odds", "bookmaker"]

    missing = [k for k in required if k not in data]
    if missing:
        return jsonify({"success": False, "error": f"Missing fields: {missing}"}), 400

    try:
        con = duckdb.connect(str(DB_PATH))

        # Check if already placed
        existing = con.execute(
            "SELECT bet_hash FROM placed_bets WHERE bet_hash = ?",
            [data["bet_hash"]]
        ).fetchone()

        if existing:
            con.close()
            return jsonify({"success": False, "error": "Bet already placed"}), 409

        # Insert new placed bet
        con.execute("""
            INSERT INTO placed_bets (
                bet_hash, game_id, home_team, away_team, game_time,
                market, bet_on, line, model_prob, model_ev, recommended_size,
                actual_size, odds, bookmaker, placed_at, status
            ) VALUES (?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, 'pending')
        """, [
            data["bet_hash"],
            data["game_id"],
            data["home_team"],
            data["away_team"],
            data["game_time"],
            data["market"],
            data["bet_on"],
            data.get("line"),
            data["model_prob"],
            data["model_ev"],
            data["recommended_size"],
            data.get("actual_size", data["recommended_size"]),
            data["odds"],
            data["bookmaker"],
            datetime.now().isoformat()
        ])

        con.close()
        return jsonify({"success": True, "message": "Bet placed successfully"})

    except Exception as e:
        return jsonify({"success": False, "error": str(e)}), 500


@app.route("/api/remove-bet", methods=["POST"])
def remove_bet():
    """Remove a bet from placed status."""
    data = request.json
    bet_hash = data.get("bet_hash")

    if not bet_hash:
        return jsonify({"success": False, "error": "bet_hash required"}), 400

    try:
        con = duckdb.connect(str(DB_PATH))
        result = con.execute(
            "DELETE FROM placed_bets WHERE bet_hash = ? RETURNING bet_hash",
            [bet_hash]
        ).fetchone()
        con.close()

        if result:
            return jsonify({"success": True, "message": "Bet removed"})
        else:
            return jsonify({"success": False, "error": "Bet not found"}), 404

    except Exception as e:
        return jsonify({"success": False, "error": str(e)}), 500


@app.route("/api/placed-bets", methods=["GET"])
def get_placed_bets():
    """Get all placed bets."""
    try:
        con = duckdb.connect(str(DB_PATH))
        bets = con.execute("""
            SELECT * FROM placed_bets
            WHERE status = 'pending'
            ORDER BY placed_at DESC
        """).fetchdf()
        con.close()

        return jsonify(bets.to_dict(orient="records"))

    except Exception as e:
        return jsonify({"success": False, "error": str(e)}), 500


@app.route("/api/exposure", methods=["GET"])
def get_exposure():
    """Get current exposure summary by game."""
    try:
        con = duckdb.connect(str(DB_PATH))

        exposure = con.execute("""
            SELECT
                game_id,
                home_team || ' vs ' || away_team as game,
                COUNT(*) as num_bets,
                SUM(COALESCE(actual_size, recommended_size)) as total_exposure,
                STRING_AGG(DISTINCT market, ', ') as markets
            FROM placed_bets
            WHERE status = 'pending'
            GROUP BY game_id, home_team, away_team
            ORDER BY total_exposure DESC
        """).fetchdf()

        con.close()
        return jsonify(exposure.to_dict(orient="records"))

    except Exception as e:
        return jsonify({"success": False, "error": str(e)}), 500


# =============================================================================
# PARLAY & PROPS API
# =============================================================================

@app.route("/api/parlay", methods=["POST"])
def calculate_parlay():
    """Calculate fair odds for a parlay."""
    data = request.json
    legs = data.get("legs", [])

    if not legs or len(legs) < 2:
        return jsonify({"success": False, "error": "At least 2 legs required"}), 400

    try:
        answer_keys_dir = BASE_DIR.parent

        # Build command: Rscript parlay.R "1H home -3" "1H under 22.5" ...
        cmd = ["Rscript", str(answer_keys_dir / "parlay.R")] + legs

        result = subprocess.run(
            cmd,
            capture_output=True,
            text=True,
            cwd=str(answer_keys_dir),
            timeout=60
        )

        if result.returncode != 0:
            return jsonify({
                "success": False,
                "error": f"Parlay calculation failed: {result.stderr[:500]}"
            }), 500

        # Parse the output to extract results
        output = result.stdout
        return jsonify({
            "success": True,
            "output": output,
            "legs": legs
        })

    except subprocess.TimeoutExpired:
        return jsonify({"success": False, "error": "Calculation timed out"}), 500
    except Exception as e:
        return jsonify({"success": False, "error": str(e)}), 500


@app.route("/api/props", methods=["POST"])
def calculate_props():
    """Calculate fair odds for a prop bet."""
    data = request.json

    column = data.get("column")
    condition = data.get("condition")
    period = data.get("period")  # Optional
    team = data.get("team")  # Optional

    if not column or not condition:
        return jsonify({"success": False, "error": "column and condition required"}), 400

    try:
        answer_keys_dir = BASE_DIR.parent
        pbp_db = answer_keys_dir / "pbp.duckdb"

        # Build R code to run props calculation
        period_arg = f'"{period}"' if period else "NULL"
        team_arg = f'"{team}"' if team else "NULL"

        r_code = f'''
        suppressPackageStartupMessages({{
            setwd("{answer_keys_dir}")
            source("props.R")
        }})

        result <- prop_odds("{column}", "{condition}", period = {period_arg}, team = {team_arg})

        # Build histogram data (binned counts)
        hist_data <- tryCatch({{
            vals <- result$raw_values
            if (length(unique(vals)) <= 10) {{
                # Discrete data - count each value
                tbl <- table(vals)
                list(bins = as.numeric(names(tbl)), counts = as.numeric(tbl), type = "discrete")
            }} else {{
                # Continuous data - use histogram bins
                h <- hist(vals, breaks = "Sturges", plot = FALSE)
                list(bins = h$mids, counts = h$counts, type = "continuous")
            }}
        }}, error = function(e) NULL)

        # Output as JSON
        cat("---JSON_START---\\n")
        cat(jsonlite::toJSON(list(
            column = result$column,
            condition = result$condition,
            period = result$period,
            team = result$team,
            prob = result$prob,
            fair_american_odds = result$fair_american_odds,
            fair_decimal_odds = result$fair_decimal_odds,
            n_games = result$n_games,
            n_hits = result$n_hits,
            stat_min = as.numeric(result$stat_distribution["Min."]),
            stat_q1 = as.numeric(result$stat_distribution["1st Qu."]),
            stat_median = as.numeric(result$stat_distribution["Median"]),
            stat_mean = as.numeric(result$stat_distribution["Mean"]),
            stat_q3 = as.numeric(result$stat_distribution["3rd Qu."]),
            stat_max = as.numeric(result$stat_distribution["Max."]),
            histogram = hist_data
        ), auto_unbox = TRUE))
        cat("\\n---JSON_END---\\n")
        '''

        result = subprocess.run(
            ["Rscript", "-e", r_code],
            capture_output=True,
            text=True,
            cwd=str(answer_keys_dir),
            timeout=60
        )

        # Check if we got valid output (even if there were warnings)
        output = result.stdout
        if "---JSON_START---" in output and "---JSON_END---" in output:
            json_str = output.split("---JSON_START---")[1].split("---JSON_END---")[0].strip()
            props_result = json.loads(json_str)
            return jsonify({"success": True, "result": props_result})
        else:
            # No JSON output - actual error
            error_msg = result.stderr if result.stderr else "Unknown error - no output"
            # Clean up common R noise
            error_msg = error_msg.replace("Warning message:", "").strip()
            return jsonify({
                "success": False,
                "error": f"Props calculation failed: {error_msg[:300]}"
            }), 500

    except subprocess.TimeoutExpired:
        return jsonify({"success": False, "error": "Calculation timed out"}), 500
    except Exception as e:
        return jsonify({"success": False, "error": str(e)}), 500


@app.route("/api/props-custom", methods=["POST"])
def calculate_props_custom():
    """Calculate fair odds using a custom SQL WHERE clause."""
    data = request.json
    query = data.get("query", "").strip()

    if not query:
        return jsonify({"success": False, "error": "Query is required"}), 400

    # Basic SQL injection prevention - only allow SELECT-like conditions
    forbidden = ["drop", "delete", "insert", "update", "alter", "create", "truncate", ";", "--"]
    query_lower = query.lower()
    for word in forbidden:
        if word in query_lower:
            return jsonify({"success": False, "error": f"Forbidden keyword: {word}"}), 400

    try:
        answer_keys_dir = BASE_DIR.parent
        pbp_db = answer_keys_dir / "pbp.duckdb"

        con = duckdb.connect(str(pbp_db), read_only=True)

        # Get sample game IDs from the samples table
        samples = con.execute("""
            SELECT DISTINCT b.game_id
            FROM nfl_samples_temp s
            JOIN (
                SELECT game_id, game_date, home_team, away_team FROM nfl_betting_pbp
                UNION ALL
                SELECT game_id, game_date, home_team, away_team FROM nfl_pre_20_betting_history
            ) b ON s.game_date = b.game_date
               AND s.home_team = b.home_team
               AND s.away_team = b.away_team
        """).fetchall()

        if not samples:
            con.close()
            return jsonify({"success": False, "error": "No sample games found"}), 500

        game_ids = [s[0] for s in samples]
        game_ids_str = "','".join(game_ids)

        # Count games where condition is met
        count_query = f"""
            SELECT COUNT(DISTINCT game_id) as hits
            FROM nfl_pbp
            WHERE game_id IN ('{game_ids_str}')
            AND ({query})
        """

        hits = con.execute(count_query).fetchone()[0]
        total = len(game_ids)
        con.close()

        prob = hits / total if total > 0 else 0

        # Convert to odds
        if prob > 0 and prob < 1:
            fair_decimal = 1 / prob
            if prob >= 0.5:
                fair_american = -round(prob / (1 - prob) * 100)
            else:
                fair_american = round((1 - prob) / prob * 100)
        else:
            fair_decimal = float('inf') if prob == 0 else 1
            fair_american = None

        return jsonify({
            "success": True,
            "result": {
                "query": query,
                "prob": prob,
                "fair_american_odds": fair_american,
                "fair_decimal_odds": round(fair_decimal, 2) if fair_decimal != float('inf') else None,
                "n_games": total,
                "n_hits": hits
            }
        })

    except Exception as e:
        return jsonify({"success": False, "error": str(e)}), 500


@app.route("/api/columns", methods=["GET"])
def get_columns():
    """Get available columns for props calculation."""
    try:
        answer_keys_dir = BASE_DIR.parent
        pbp_db = answer_keys_dir / "pbp.duckdb"

        con = duckdb.connect(str(pbp_db), read_only=True)
        columns = con.execute("""
            SELECT column_name
            FROM information_schema.columns
            WHERE table_name = 'nfl_pbp'
            ORDER BY column_name
        """).fetchall()
        con.close()

        column_list = [c[0] for c in columns]

        # Highlight common prop columns
        common = [
            "touchdown", "pass_touchdown", "rush_touchdown", "return_touchdown",
            "field_goal_attempt", "field_goal_result", "safety",
            "extra_point_attempt", "two_point_attempt", "interception",
            "fumble", "sack", "penalty"
        ]

        return jsonify({
            "success": True,
            "columns": column_list,
            "common": [c for c in common if c in column_list]
        })

    except Exception as e:
        return jsonify({"success": False, "error": str(e)}), 500


@app.route("/refresh", methods=["POST"])
def refresh():
    """Run full pipeline (scrapers + predictions) and regenerate dashboard."""
    try:
        nfl_work_dir = BASE_DIR.parent.parent  # NFLWork directory
        answer_keys_dir = BASE_DIR.parent  # Answer Keys directory

        # Step 1: Run the full pipeline (scrapers + R prepare + R combine)
        print("Running full pipeline (scrapers + predictions)...")
        result = subprocess.run(
            ["python3", str(answer_keys_dir / "run.py"), "nfl"],
            capture_output=True,
            text=True,
            cwd=nfl_work_dir
        )

        if result.returncode != 0:
            print(f"Pipeline stderr: {result.stderr}")
            return jsonify({
                "success": False,
                "error": f"Pipeline failed: {result.stderr[:500]}"
            }), 500

        # Step 2: Generate dashboard HTML from the saved data
        print("Generating dashboard HTML...")
        result = subprocess.run(
            ["Rscript", str(BASE_DIR / "nfl_dashboard.R")],
            capture_output=True,
            text=True,
            cwd=nfl_work_dir
        )

        if result.returncode != 0:
            print(f"Dashboard stderr: {result.stderr}")
            return jsonify({
                "success": False,
                "error": f"Dashboard generation failed: {result.stderr[:500]}"
            }), 500

        print("Refresh complete!")
        return jsonify({"success": True, "message": "Dashboard refreshed successfully"})

    except Exception as e:
        return jsonify({"success": False, "error": str(e)}), 500


# =============================================================================
# MAIN
# =============================================================================

if __name__ == "__main__":
    init_db()
    print("\n" + "=" * 50)
    print("NFL +EV Betting Dashboard")
    print("=" * 50)
    print(f"\nOpen in browser: http://127.0.0.1:8081")
    print("Press Ctrl+C to stop\n")
    app.run(host="127.0.0.1", port=8081, debug=False)
