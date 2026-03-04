#!/usr/bin/env python3
"""
CBB +EV Betting Dashboard Server

Flask server that:
- Serves the dashboard HTML
- Provides API for marking bets as placed/removed
- Handles refresh to regenerate predictions

Run with: python cbb_dashboard_server.py
Then open: http://localhost:8082
"""

import subprocess
import sys
import json
import hashlib
import mimetypes
from datetime import datetime
from pathlib import Path

import duckdb
from flask import Flask, Response, jsonify, request

BASE_DIR = Path(__file__).parent.resolve()
DB_PATH = BASE_DIR / "cbb_dashboard.duckdb"
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
            game_time TIMESTAMP,
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

    con.execute("""
        CREATE TABLE IF NOT EXISTS book_settings (
            bookmaker_key TEXT PRIMARY KEY,
            enabled BOOLEAN DEFAULT FALSE,
            discovered_at TIMESTAMP DEFAULT CURRENT_TIMESTAMP
        )
    """)

    con.execute("""
        CREATE TABLE IF NOT EXISTS sizing_settings (
            param TEXT PRIMARY KEY,
            value REAL NOT NULL
        )
    """)
    # Seed defaults if empty
    con.execute("""
        INSERT INTO sizing_settings (param, value) VALUES ('bankroll', 100)
        ON CONFLICT (param) DO NOTHING
    """)
    con.execute("""
        INSERT INTO sizing_settings (param, value) VALUES ('kelly_mult', 0.25)
        ON CONFLICT (param) DO NOTHING
    """)

    con.execute("""
        CREATE TABLE IF NOT EXISTS filter_settings (
            filter_type TEXT PRIMARY KEY,
            selected_values TEXT NOT NULL
        )
    """)
    # Default: Status shows only "Not Placed" (hide already-placed bets)
    con.execute("""
        INSERT INTO filter_settings (filter_type, selected_values)
        VALUES ('status', '["Not Placed"]')
        ON CONFLICT (filter_type) DO NOTHING
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
            "<p>Click the refresh button or run: <code>Rscript cbb_dashboard.R</code></p>",
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
            None if data["game_time"] in ("NA", "", None) else data["game_time"],
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


@app.route("/api/update-bet", methods=["POST"])
def update_bet():
    """Update a placed bet's actual size."""
    data = request.json
    bet_hash = data.get("bet_hash")
    actual_size = data.get("actual_size")

    if not bet_hash:
        return jsonify({"success": False, "error": "bet_hash required"}), 400
    if actual_size is None:
        return jsonify({"success": False, "error": "actual_size required"}), 400

    try:
        con = duckdb.connect(str(DB_PATH))
        result = con.execute(
            "UPDATE placed_bets SET actual_size = ? WHERE bet_hash = ? RETURNING bet_hash",
            [float(actual_size), bet_hash]
        ).fetchone()
        con.close()

        if result:
            return jsonify({"success": True, "message": "Bet updated"})
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


@app.route("/api/book-settings", methods=["GET"])
def get_book_settings():
    """Get all bookmaker enabled/disabled settings."""
    try:
        con = duckdb.connect(str(DB_PATH))
        rows = con.execute("SELECT bookmaker_key, enabled FROM book_settings").fetchall()
        con.close()
        return jsonify({row[0]: row[1] for row in rows})
    except Exception as e:
        return jsonify({}), 200


@app.route("/api/book-settings", methods=["POST"])
def update_book_setting():
    """Update a single bookmaker's enabled/disabled state."""
    data = request.json
    book = data.get("book")
    enabled = data.get("enabled", False)

    if not book:
        return jsonify({"success": False, "error": "book required"}), 400

    try:
        con = duckdb.connect(str(DB_PATH))
        con.execute("""
            INSERT INTO book_settings (bookmaker_key, enabled)
            VALUES (?, ?)
            ON CONFLICT (bookmaker_key) DO UPDATE SET enabled = ?
        """, [book, enabled, enabled])
        con.close()
        return jsonify({"success": True})
    except Exception as e:
        return jsonify({"success": False, "error": str(e)}), 500


@app.route("/api/book-settings/bulk", methods=["POST"])
def bulk_add_books():
    """Register newly discovered books as disabled by default."""
    data = request.json
    books = data.get("books", [])

    try:
        con = duckdb.connect(str(DB_PATH))
        for book in books:
            con.execute("""
                INSERT INTO book_settings (bookmaker_key, enabled)
                VALUES (?, FALSE)
                ON CONFLICT (bookmaker_key) DO NOTHING
            """, [book])
        con.close()
        return jsonify({"success": True, "added": len(books)})
    except Exception as e:
        return jsonify({"success": False, "error": str(e)}), 500


@app.route("/api/sizing-settings", methods=["GET"])
def get_sizing_settings():
    """Get bankroll and kelly multiplier settings."""
    try:
        con = duckdb.connect(str(DB_PATH))
        rows = con.execute("SELECT param, value FROM sizing_settings").fetchall()
        con.close()
        return jsonify({row[0]: row[1] for row in rows})
    except Exception as e:
        return jsonify({"bankroll": 100, "kelly_mult": 0.25}), 200


@app.route("/api/sizing-settings", methods=["POST"])
def update_sizing_settings():
    """Update bankroll and/or kelly multiplier."""
    data = request.json
    try:
        con = duckdb.connect(str(DB_PATH))
        for param in ("bankroll", "kelly_mult"):
            if param in data:
                con.execute("""
                    INSERT INTO sizing_settings (param, value) VALUES (?, ?)
                    ON CONFLICT (param) DO UPDATE SET value = ?
                """, [param, float(data[param]), float(data[param])])
        con.close()
        return jsonify({"success": True})
    except Exception as e:
        return jsonify({"success": False, "error": str(e)}), 500


@app.route("/api/filter-settings", methods=["GET"])
def get_filter_settings():
    """Get persisted filter states (market, correlation, status)."""
    try:
        con = duckdb.connect(str(DB_PATH))
        rows = con.execute("SELECT filter_type, selected_values FROM filter_settings").fetchall()
        con.close()
        return jsonify({row[0]: json.loads(row[1]) for row in rows})
    except Exception as e:
        return jsonify({}), 200


@app.route("/api/filter-settings", methods=["POST"])
def update_filter_settings():
    """Save a filter type's selected values."""
    data = request.json
    filter_type = data.get("filter_type")
    selected_values = data.get("selected_values")

    if not filter_type or filter_type not in ("market", "correlation", "status"):
        return jsonify({"success": False, "error": "Invalid filter_type"}), 400

    try:
        con = duckdb.connect(str(DB_PATH))
        con.execute("""
            INSERT INTO filter_settings (filter_type, selected_values)
            VALUES (?, ?)
            ON CONFLICT (filter_type) DO UPDATE SET selected_values = ?
        """, [filter_type, json.dumps(selected_values), json.dumps(selected_values)])
        con.close()
        return jsonify({"success": True})
    except Exception as e:
        return jsonify({"success": False, "error": str(e)}), 500


def run_pipeline():
    """Run full pipeline (scrapers + predictions) and regenerate dashboard."""
    nfl_work_dir = BASE_DIR.parent.parent  # NFLWork directory
    answer_keys_dir = BASE_DIR.parent  # Answer Keys directory

    # Step 1: Run the full pipeline (scrapers + R prepare + R combine)
    print("Running full CBB pipeline (scrapers + predictions)...")
    result = subprocess.run(
        [sys.executable, str(answer_keys_dir / "run.py"), "cbb"],
        capture_output=True,
        text=True,
        cwd=nfl_work_dir
    )

    if result.returncode != 0:
        print(f"Pipeline stderr: {result.stderr}")
        return False, f"Pipeline failed: {result.stderr[:500]}"

    # Step 2: Generate dashboard HTML from the saved data
    print("Generating CBB dashboard HTML...")
    result = subprocess.run(
        ["Rscript", str(BASE_DIR / "cbb_dashboard.R")],
        capture_output=True,
        text=True,
        cwd=nfl_work_dir
    )

    if result.returncode != 0:
        print(f"Dashboard stderr: {result.stderr}")
        return False, f"Dashboard generation failed: {result.stderr[:500]}"

    print("CBB refresh complete!")
    return True, "Dashboard refreshed successfully"


@app.route("/refresh", methods=["POST"])
def refresh():
    """Run full pipeline and regenerate dashboard."""
    try:
        success, message = run_pipeline()
        if success:
            return jsonify({"success": True, "message": message})
        else:
            return jsonify({"success": False, "error": message}), 500
    except Exception as e:
        return jsonify({"success": False, "error": str(e)}), 500


# =============================================================================
# MAIN
# =============================================================================

if __name__ == "__main__":
    init_db()

    print("\n" + "=" * 50)
    print("CBB +EV Betting Dashboard")
    print("=" * 50)
    print(f"\nOpen in browser: http://localhost:8082")
    print("Press Ctrl+C to stop\n")
    app.run(host="0.0.0.0", port=8082, debug=False)
