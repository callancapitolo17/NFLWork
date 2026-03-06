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
import logging
import mimetypes
import threading
from datetime import datetime, timedelta, timezone
from pathlib import Path

import duckdb
from flask import Flask, Response, jsonify, request

BASE_DIR = Path(__file__).parent.resolve()
PROJECT_ROOT = BASE_DIR.parent.parent  # NFLWork/
DB_PATH = BASE_DIR / "cbb_dashboard.duckdb"
app = Flask(__name__)
log = logging.getLogger("clv")

# Scraper configs for closing-odds capture (paths relative to PROJECT_ROOT)
OFFSHORE_SCRAPERS = {
    "wagerzon": {
        "script": "wagerzon_odds/scraper_v2.py",
        "db": "wagerzon_odds/wagerzon.duckdb",
        "table": "cbb_odds",
    },
    "hoop88": {
        "script": "hoop88_odds/scraper.py",
        "db": "hoop88_odds/hoop88.duckdb",
        "table": "cbb_odds",
    },
    "bfa": {
        "script": "bfa_odds/scraper.py",
        "db": "bfa_odds/bfa.duckdb",
        "table": "cbb_odds",
    },
    "bookmaker": {
        "script": "bookmaker_odds/scraper.py",
        "db": "bookmaker_odds/bookmaker.duckdb",
        "table": "cbb_odds",
    },
}

# Active capture timers: game_id -> threading.Timer
scheduled_captures = {}


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

    # CLV tracking tables
    con.execute("""
        CREATE TABLE IF NOT EXISTS closing_snapshots (
            snapshot_time TIMESTAMP NOT NULL,
            bookmaker TEXT NOT NULL,
            game_id TEXT,
            home_team TEXT NOT NULL,
            away_team TEXT NOT NULL,
            market TEXT NOT NULL,
            bet_on TEXT NOT NULL,
            line REAL,
            odds INTEGER NOT NULL,
            counter_odds INTEGER
        )
    """)

    con.execute("""
        CREATE TABLE IF NOT EXISTS bet_clv (
            bet_hash TEXT PRIMARY KEY,
            game_id TEXT,
            game_time TIMESTAMP,
            bookmaker TEXT,
            market TEXT,
            bet_on TEXT,
            placement_line REAL,
            placement_odds INTEGER,
            placement_novig_prob REAL,
            market_closing_line REAL,
            market_closing_odds INTEGER,
            market_closing_counter_odds INTEGER,
            market_closing_novig_prob REAL,
            market_clv REAL,
            book_closing_line REAL,
            book_closing_odds INTEGER,
            book_closing_counter_odds INTEGER,
            book_closing_novig_prob REAL,
            book_clv REAL,
            line_moved BOOLEAN,
            clv_method TEXT,
            sigma_used REAL,
            computed_at TIMESTAMP DEFAULT CURRENT_TIMESTAMP
        )
    """)

    con.close()
    print(f"Database initialized at: {DB_PATH}")


# =============================================================================
# CLV CLOSING-ODDS CAPTURE
# =============================================================================

def run_scraper(name: str, config: dict):
    """Run a single offshore scraper and return success status."""
    script = PROJECT_ROOT / config["script"]
    venv_python = script.parent / "venv" / "bin" / "python"
    python_exe = str(venv_python) if venv_python.exists() else sys.executable

    result = subprocess.run(
        [python_exe, str(script), "cbb"],
        capture_output=True, text=True,
        cwd=str(script.parent),
        timeout=120,
    )
    if result.returncode != 0:
        log.warning("Scraper %s failed: %s", name, result.stderr[-500:])
    return result.returncode == 0


def snapshot_book_odds(bookmaker: str, config: dict, game_teams: list[tuple[str, str]]):
    """Read current odds from a scraper's DuckDB and save to closing_snapshots."""
    db_path = PROJECT_ROOT / config["db"]
    if not db_path.exists():
        log.warning("DB not found for %s: %s", bookmaker, db_path)
        return 0

    now = datetime.now()
    rows = []

    try:
        src = duckdb.connect(str(db_path), read_only=True)
        odds = src.execute(f"SELECT * FROM {config['table']}").fetchall()
        cols = [d[0] for d in src.description]
        src.close()
    except Exception as e:
        log.warning("Failed to read %s DB: %s", bookmaker, e)
        return 0

    for row in odds:
        r = dict(zip(cols, row))
        home = r.get("home_team", "")
        away = r.get("away_team", "")
        market = r.get("market", "")
        game_id = r.get("game_id", "")

        # Spread records
        if r.get("away_spread") is not None and r.get("away_spread_price") is not None:
            rows.append((now, bookmaker, game_id, home, away, market,
                         "away", r["away_spread"], r["away_spread_price"], r.get("home_spread_price")))
            rows.append((now, bookmaker, game_id, home, away, market,
                         "home", r["home_spread"], r["home_spread_price"], r.get("away_spread_price")))

        # Total records
        if r.get("total") is not None and r.get("over_price") is not None:
            rows.append((now, bookmaker, game_id, home, away, market,
                         "over", r["total"], r["over_price"], r.get("under_price")))
            rows.append((now, bookmaker, game_id, home, away, market,
                         "under", r["total"], r["under_price"], r.get("over_price")))

    if not rows:
        return 0

    try:
        dst = duckdb.connect(str(DB_PATH))
        dst.executemany("""
            INSERT INTO closing_snapshots
                (snapshot_time, bookmaker, game_id, home_team, away_team,
                 market, bet_on, line, odds, counter_odds)
            VALUES (?, ?, ?, ?, ?, ?, ?, ?, ?, ?)
        """, rows)
        dst.close()
    except Exception as e:
        log.warning("Failed to write closing_snapshots for %s: %s", bookmaker, e)
        return 0

    return len(rows)


def run_closing_capture(game_id: str, bookmakers: list[str]):
    """Run offshore scrapers and snapshot odds for a specific game."""
    log.info("Capturing closing odds for game %s (books: %s)", game_id, bookmakers)

    for book in bookmakers:
        config = OFFSHORE_SCRAPERS.get(book)
        if not config:
            continue
        run_scraper(book, config)

        # Get teams for this game from placed_bets
        try:
            con = duckdb.connect(str(DB_PATH), read_only=True)
            teams = con.execute(
                "SELECT DISTINCT home_team, away_team FROM placed_bets WHERE game_id = ?",
                [game_id]
            ).fetchall()
            con.close()
        except Exception:
            teams = []

        n = snapshot_book_odds(book, config, teams)
        log.info("Snapshotted %d odds rows for %s", n, book)

    # Cleanup timer reference
    scheduled_captures.pop(game_id, None)


def schedule_capture(game_id: str, game_time_str: str, bookmaker: str):
    """Schedule a closing-odds capture for 15 min before game time."""
    if game_id in scheduled_captures:
        # Already scheduled — just note the additional bookmaker
        return

    try:
        # game_time from dashboard is UTC (formatted with trailing Z)
        raw = str(game_time_str).replace("Z", "+00:00")
        game_time = datetime.fromisoformat(raw)
        if game_time.tzinfo is None:
            game_time = game_time.replace(tzinfo=timezone.utc)
    except (ValueError, TypeError):
        log.warning("Cannot parse game_time for CLV capture: %s", game_time_str)
        return

    capture_at = game_time - timedelta(minutes=15)
    delay = (capture_at - datetime.now(timezone.utc)).total_seconds()

    if delay <= 0:
        log.info("Game %s already within 15 min of tipoff, capturing now", game_id)
        threading.Thread(
            target=run_closing_capture, args=(game_id, [bookmaker]), daemon=True
        ).start()
        return

    # Gather all bookmakers with bets on this game
    try:
        con = duckdb.connect(str(DB_PATH), read_only=True)
        books = [r[0] for r in con.execute(
            "SELECT DISTINCT bookmaker FROM placed_bets WHERE game_id = ? AND status = 'pending'",
            [game_id]
        ).fetchall()]
        con.close()
    except Exception:
        books = [bookmaker]

    timer = threading.Timer(delay, run_closing_capture, args=[game_id, books])
    timer.daemon = True
    timer.start()
    scheduled_captures[game_id] = timer

    log.info("Scheduled CLV capture for game %s at %s (in %.0f min)",
             game_id, capture_at.strftime("%H:%M"), delay / 60)


def schedule_pending_captures():
    """On startup, schedule captures for any upcoming games with placed bets."""
    try:
        con = duckdb.connect(str(DB_PATH), read_only=True)
        upcoming = con.execute("""
            SELECT game_id, game_time, STRING_AGG(DISTINCT bookmaker, ',') as books
            FROM placed_bets
            WHERE status = 'pending'
              AND game_time > CURRENT_TIMESTAMP
            GROUP BY game_id, game_time
        """).fetchall()
        con.close()
    except Exception as e:
        log.warning("Failed to load pending captures on startup: %s", e)
        return

    for game_id, game_time, books_csv in upcoming:
        books = books_csv.split(",")
        # game_time from DB is naive UTC; compare against UTC now
        game_time_utc = game_time.replace(tzinfo=timezone.utc) if game_time.tzinfo is None else game_time
        capture_at = game_time_utc - timedelta(minutes=15)
        delay = (capture_at - datetime.now(timezone.utc)).total_seconds()

        if delay <= 0:
            continue  # Game already started or too close

        timer = threading.Timer(delay, run_closing_capture, args=[game_id, books])
        timer.daemon = True
        timer.start()
        scheduled_captures[game_id] = timer
        log.info("Startup: scheduled CLV capture for game %s at %s",
                 game_id, capture_at.strftime("%H:%M"))

    if upcoming:
        log.info("Scheduled %d CLV captures on startup", len(scheduled_captures))


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

        # Schedule CLV closing-odds capture for this game (best-effort, don't fail the bet)
        try:
            if data.get("game_time") not in ("NA", "", None):
                schedule_capture(data["game_id"], data["game_time"], data["bookmaker"])
        except Exception as e:
            log.warning("Failed to schedule CLV capture: %s", e)

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

    if not filter_type or filter_type not in ("market", "correlation", "size", "status"):
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


# =============================================================================
# CLV API ENDPOINTS
# =============================================================================

@app.route("/api/clv-summary", methods=["GET"])
def clv_summary():
    """Aggregate CLV stats: overall, by market type, and by book."""
    try:
        con = duckdb.connect(str(DB_PATH), read_only=True)

        overall = con.execute("""
            SELECT
                COUNT(*) as total_bets,
                COUNT(market_clv) as market_clv_count,
                COUNT(book_clv) as book_clv_count,
                AVG(market_clv) as avg_market_clv,
                AVG(book_clv) as avg_book_clv,
                COUNT(CASE WHEN market_clv > 0 THEN 1 END) as market_clv_positive,
                COUNT(CASE WHEN book_clv > 0 THEN 1 END) as book_clv_positive
            FROM bet_clv
        """).fetchdf().to_dict(orient="records")[0]

        by_market = con.execute("""
            SELECT
                market,
                COUNT(*) as n,
                AVG(market_clv) as avg_market_clv,
                AVG(book_clv) as avg_book_clv,
                COUNT(CASE WHEN market_clv > 0 THEN 1 END) as positive_pct
            FROM bet_clv
            WHERE market_clv IS NOT NULL
            GROUP BY market
            ORDER BY n DESC
        """).fetchdf().to_dict(orient="records")

        by_book = con.execute("""
            SELECT
                bookmaker,
                COUNT(*) as n,
                AVG(market_clv) as avg_market_clv,
                AVG(book_clv) as avg_book_clv
            FROM bet_clv
            WHERE market_clv IS NOT NULL
            GROUP BY bookmaker
            ORDER BY n DESC
        """).fetchdf().to_dict(orient="records")

        con.close()
        return jsonify({
            "overall": overall,
            "by_market": by_market,
            "by_book": by_book,
        })
    except Exception as e:
        return jsonify({"success": False, "error": str(e)}), 500


@app.route("/api/clv-details", methods=["GET"])
def clv_details():
    """Per-bet CLV data for inspection."""
    try:
        con = duckdb.connect(str(DB_PATH), read_only=True)
        rows = con.execute("""
            SELECT c.*, p.home_team, p.away_team, p.model_ev
            FROM bet_clv c
            JOIN placed_bets p ON c.bet_hash = p.bet_hash
            ORDER BY c.game_time DESC
        """).fetchdf()
        con.close()

        return jsonify(rows.to_dict(orient="records"))
    except Exception as e:
        return jsonify({"success": False, "error": str(e)}), 500


@app.route("/api/clv-compute", methods=["POST"])
def trigger_clv_compute():
    """Trigger CLV computation for completed games."""
    try:
        script = BASE_DIR.parent / "CBB Answer Key" / "clv_compute.py"
        result = subprocess.run(
            [sys.executable, str(script)],
            capture_output=True, text=True,
            cwd=str(script.parent),
            timeout=300,
        )
        if result.returncode == 0:
            return jsonify({"success": True, "output": result.stdout[-1000:]})
        else:
            return jsonify({"success": False, "error": result.stderr[-500:]}), 500
    except Exception as e:
        return jsonify({"success": False, "error": str(e)}), 500


@app.route("/api/scheduled-captures", methods=["GET"])
def get_scheduled_captures():
    """Show currently scheduled CLV capture timers."""
    captures = []
    for game_id, timer in scheduled_captures.items():
        captures.append({
            "game_id": game_id,
            "active": timer.is_alive() if timer else False,
        })
    return jsonify(captures)



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
    logging.basicConfig(
        level=logging.INFO,
        format="%(asctime)s [%(name)s] %(message)s",
        datefmt="%H:%M:%S",
    )

    init_db()
    schedule_pending_captures()

    print("\n" + "=" * 50)
    print("CBB +EV Betting Dashboard")
    print("=" * 50)
    print(f"\nOpen in browser: http://localhost:8082")
    print("Press Ctrl+C to stop\n")
    app.run(host="0.0.0.0", port=8082, debug=False)
