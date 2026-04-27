#!/usr/bin/env python3
"""
MLB +EV Betting Dashboard Server

Flask server that:
- Serves the dashboard HTML
- Provides API for marking bets as placed/removed
- Handles refresh to regenerate predictions

Run with: python mlb_dashboard_server.py
Then open: http://localhost:8083
"""

import subprocess
import sys
import json
import json as _json
import logging
import mimetypes
import threading
import time
from datetime import datetime, timedelta, timezone
from pathlib import Path

import duckdb
from flask import Flask, Response, jsonify, request

BASE_DIR = Path(__file__).parent.resolve()
PROJECT_ROOT = BASE_DIR.parent.parent  # NFLWork/
DB_PATH = BASE_DIR / "mlb_dashboard.duckdb"
MLB_DB_PATH = BASE_DIR.parent / "mlb.duckdb"

# Combined parlay pricing imports — pull from wagerzon_odds/ + local helper
sys.path.insert(0, str(PROJECT_ROOT / "wagerzon_odds"))
from combined_parlay import joint_pricing
from parlay_pricer import (
    get_combined_parlay_price as wz_get_combined_parlay_price,
    get_wz_session as _get_wz_session,
)

# Combined parlay pricing cache: keyed on sorted (hash_a, hash_b) tuple.
# TTL prevents stale prices when WZ moves their lines. 60s is a safe default.
_COMBO_PRICE_CACHE: dict[tuple[str, str], tuple[float, dict]] = {}
_COMBO_PRICE_TTL_SECONDS = 60
app = Flask(__name__)

# Make wagerzon_odds importable for parlay placement
REPO_ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(REPO_ROOT / "wagerzon_odds"))
import parlay_placer  # noqa: E402

MLB_DB = REPO_ROOT / "Answer Keys" / "mlb.duckdb"
DASHBOARD_DB = DB_PATH
log = logging.getLogger("clv")

# Books with working navigators for auto-queue
SUPPORTED_AUTO_BOOKS = ("wagerzon", "hoop88", "bfa", "betonlineag")

# DuckDB with pipeline bets — always in main repo (not worktree)
_REPO_ROOT = PROJECT_ROOT
if ".claude/worktrees" in str(_REPO_ROOT):
    _REPO_ROOT = Path(str(_REPO_ROOT).split(".claude/worktrees")[0])
# Scraper configs for closing-odds capture (paths relative to PROJECT_ROOT)
OFFSHORE_SCRAPERS = {
    "wagerzon": {
        "script": "wagerzon_odds/scraper_v2.py",
        "db": "wagerzon_odds/wagerzon.duckdb",
        "table": "mlb_odds",
    },
    "hoop88": {
        "script": "hoop88_odds/scraper.py",
        "db": "hoop88_odds/hoop88.duckdb",
        "table": "mlb_odds",
    },
    "bfa": {
        "script": "bfa_odds/scraper.py",
        "db": "bfa_odds/bfa.duckdb",
        "table": "mlb_odds",
    },
    "bookmaker": {
        "script": "bookmaker_odds/scraper.py",
        "db": "bookmaker_odds/bookmaker.duckdb",
        "table": "mlb_odds",
    },
}

# Active capture timers: game_id -> threading.Timer
scheduled_captures = {}

# Mutex to prevent concurrent pipeline/refresh runs
_refresh_lock = threading.Lock()


# =============================================================================
# DATABASE INITIALIZATION
# =============================================================================

def init_db():
    """Create tables if they don't exist."""
    con = duckdb.connect(str(DB_PATH))
    try:
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
            INSERT INTO sizing_settings (param, value) VALUES ('parlay_bankroll', 100)
            ON CONFLICT (param) DO NOTHING
        """)
        con.execute("""
            INSERT INTO sizing_settings (param, value) VALUES ('parlay_kelly_mult', 0.25)
            ON CONFLICT (param) DO NOTHING
        """)
        con.execute("""
            INSERT INTO sizing_settings (param, value) VALUES ('parlay_min_edge', 0)
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

        # Parlay bet tracking
        con.execute("""
            CREATE TABLE IF NOT EXISTS placed_parlays (
                parlay_hash TEXT PRIMARY KEY,
                game_id TEXT,
                home_team TEXT,
                away_team TEXT,
                game_time TIMESTAMP,
                combo TEXT,
                spread_line REAL,
                total_line REAL,
                fair_odds INTEGER,
                wz_odds INTEGER,
                edge_pct REAL,
                kelly_bet REAL,
                actual_size REAL,
                placed_at TIMESTAMP DEFAULT CURRENT_TIMESTAMP,
                status TEXT DEFAULT 'pending',
                -- columns added for auto-placement (Task 9)
                recommended_size REAL,
                expected_odds INTEGER,
                expected_win REAL,
                actual_win REAL,
                ticket_number TEXT,
                idwt INTEGER,
                legs_json TEXT,
                error_msg TEXT,
                error_msg_key TEXT,
                updated_at TIMESTAMP
            )
        """)
        # Migration: add auto-placement columns to existing placed_parlays tables
        # DuckDB ALTER TABLE ADD COLUMN is idempotent-ish but raises if col exists.
        # We check information_schema first to make it truly safe.
        existing_cols = {
            r[0] for r in con.execute(
                "SELECT column_name FROM information_schema.columns "
                "WHERE table_name = 'placed_parlays'"
            ).fetchall()
        }
        migration_cols = {
            "recommended_size": "REAL",
            "expected_odds": "INTEGER",
            "expected_win": "REAL",
            "actual_win": "REAL",
            "ticket_number": "TEXT",
            "idwt": "INTEGER",
            "legs_json": "TEXT",
            "error_msg": "TEXT",
            "error_msg_key": "TEXT",
            "updated_at": "TIMESTAMP",
        }
        for col, dtype in migration_cols.items():
            if col not in existing_cols:
                con.execute(f"ALTER TABLE placed_parlays ADD COLUMN {col} {dtype}")

        # Orphan forensics: rare case where Wagerzon confirms placement but
        # the local placed_parlays write fails. Created here AND by
        # wagerzon_odds/migrate_placed_parlays.py so a fresh server start
        # without running the migration script still has the safety net.
        con.execute("""
            CREATE TABLE IF NOT EXISTS placement_orphans (
                idwt          BIGINT PRIMARY KEY,
                ticket_number TEXT,
                parlay_hash   TEXT,
                raw_response  TEXT,
                error         TEXT,
                created_at    TIMESTAMP DEFAULT CURRENT_TIMESTAMP
            )
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
    finally:
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
        [python_exe, str(script), "mlb"],
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
        try:
            odds = src.execute(f"SELECT * FROM {config['table']}").fetchall()
            cols = [d[0] for d in src.description]
        finally:
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
        try:
            dst.executemany("""
                INSERT INTO closing_snapshots
                    (snapshot_time, bookmaker, game_id, home_team, away_team,
                     market, bet_on, line, odds, counter_odds)
                VALUES (?, ?, ?, ?, ?, ?, ?, ?, ?, ?)
            """, rows)
        finally:
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
            try:
                teams = con.execute(
                    "SELECT DISTINCT home_team, away_team FROM placed_bets WHERE game_id = ?",
                    [game_id]
                ).fetchall()
            finally:
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
        log.info("Game %s already within 15 min of first pitch, capturing now", game_id)
        threading.Thread(
            target=run_closing_capture, args=(game_id, [bookmaker]), daemon=True
        ).start()
        return

    # Gather all bookmakers with bets on this game
    try:
        con = duckdb.connect(str(DB_PATH), read_only=True)
        try:
            books = [r[0] for r in con.execute(
                "SELECT DISTINCT bookmaker FROM placed_bets WHERE game_id = ? AND status = 'pending'",
                [game_id]
            ).fetchall()]
        finally:
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
        try:
            upcoming = con.execute("""
                SELECT game_id, game_time, STRING_AGG(DISTINCT bookmaker, ',') as books
                FROM placed_bets
                WHERE status = 'pending'
                  AND game_time > CURRENT_TIMESTAMP
                GROUP BY game_id, game_time
            """).fetchall()
        finally:
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
            "<p>Click the refresh button or run: <code>Rscript mlb_dashboard.R</code></p>",
            mimetype="text/html"
        )
    # R may write report.html in the platform default encoding (Latin-1 on
    # macOS) when the source contains characters outside ASCII. Read with
    # errors="replace" so a stray byte never 500s the dashboard.
    with open(html_path, "r", encoding="utf-8", errors="replace") as f:
        content = f.read()
    return Response(content, mimetype="text/html; charset=utf-8")


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
        try:
            # Check if already placed — allow overwrite if status is not 'pending' (i.e. queued/nav_error/ready_to_confirm)
            existing = con.execute(
                "SELECT bet_hash, status FROM placed_bets WHERE bet_hash = ?",
                [data["bet_hash"]]
            ).fetchone()

            if existing and existing[1] == "pending":
                return jsonify({"success": False, "error": "Bet already placed"}), 409

            if existing:
                # Overwrite non-confirmed bet (queued, nav_error, ready_to_confirm)
                con.execute("""
                    UPDATE placed_bets SET
                        actual_size = ?, recommended_size = ?, odds = ?,
                        placed_at = ?, status = 'pending'
                    WHERE bet_hash = ?
                """, [
                    float(data.get("actual_size", data["recommended_size"])),
                    float(data["recommended_size"]),
                    int(data["odds"]),
                    datetime.now().isoformat(),
                    data["bet_hash"],
                ])
            else:
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
        finally:
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
        try:
            result = con.execute(
                "DELETE FROM placed_bets WHERE bet_hash = ? RETURNING bet_hash",
                [bet_hash]
            ).fetchone()
        finally:
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
        try:
            result = con.execute(
                "UPDATE placed_bets SET actual_size = ? WHERE bet_hash = ? RETURNING bet_hash",
                [float(actual_size), bet_hash]
            ).fetchone()
        finally:
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
        try:
            bets = con.execute("""
                SELECT * FROM placed_bets
                WHERE status = 'pending'
                ORDER BY placed_at DESC
            """).fetchdf()
        finally:
            con.close()

        return jsonify(bets.to_dict(orient="records"))

    except Exception as e:
        return jsonify({"success": False, "error": str(e)}), 500


# =============================================================================
# PARLAY BET PLACEMENT
# =============================================================================


def _spread_play_from_combo(combo: str) -> int:
    """1 = home spread, 0 = away spread. Derived from the combo description.

    Combo strings start with the side bet on, e.g. 'Home -1.5 + Over 9.5',
    'Away +1.5 + Under 7.5', or 'F5 Home Spread + Over'. We can't infer the
    side from the line sign because home can be the underdog.
    """
    # Match 'Home' or 'Away' as a standalone word at the start or after 'F5 '
    if "Home Spread" in combo or combo.startswith("Home"):
        return 1
    if "Away Spread" in combo or combo.startswith("Away"):
        return 0
    raise ValueError(f"Cannot determine spread side from combo: {combo!r}")


def _total_play_from_combo(combo: str) -> int:
    """2 = over, 3 = under. Looks for the literal '+ Over' / '+ Under' segment."""
    if "+ Over" in combo:
        return 2
    if "+ Under" in combo:
        return 3
    raise ValueError(f"Cannot determine total side from combo: {combo!r}")


def _total_points_for_play(total_line: float, play: int) -> str:
    """Wagerzon expects negative points for over (play=2), positive for under (play=3)."""
    return str(-abs(total_line)) if play == 2 else str(total_line)


@app.route("/api/price-combined-parlay", methods=["POST"])
def price_combined_parlay():
    """Compute joint pricing + recommended Kelly stake for two parlay rows.

    Body: {"parlay_hash_a": "...", "parlay_hash_b": "..."}
    Returns: joint_fair_dec, joint_fair_prob, joint_edge, wz_dec, wz_win,
             kelly_stake, and the per-leg combo descriptions.
    """
    data = request.get_json(silent=True) or {}
    hash_a = data.get("parlay_hash_a")
    hash_b = data.get("parlay_hash_b")
    if not hash_a or not hash_b:
        return jsonify({"success": False, "error": "Missing parlay_hash_a or parlay_hash_b"}), 400
    if hash_a == hash_b:
        return jsonify({"success": False, "error": "Cannot combine a row with itself"}), 400

    # Pull the two source rows from mlb.duckdb
    mcon = duckdb.connect(str(MLB_DB_PATH))
    try:
        rows = mcon.execute("""
            SELECT parlay_hash, fair_dec, wz_dec, idgm,
                   spread_line, total_line, spread_price, total_price, combo, game_id
            FROM mlb_parlay_opportunities
            WHERE parlay_hash IN (?, ?)
        """, [hash_a, hash_b]).fetchall()
    finally:
        mcon.close()

    if len(rows) != 2:
        return jsonify({"success": False, "error": "One or both parlays not found"}), 404

    # Reject same-game combos
    if rows[0][9] == rows[1][9]:
        return jsonify({"success": False, "error": "Same-game combos not supported"}), 400

    by_hash = {r[0]: r for r in rows}
    a = by_hash[hash_a]; b = by_hash[hash_b]

    # Cache check (keyed on sorted hashes so order doesn't matter)
    cache_key = tuple(sorted([hash_a, hash_b]))
    now = time.time()
    cached = _COMBO_PRICE_CACHE.get(cache_key)
    if cached and (now - cached[0]) < _COMBO_PRICE_TTL_SECONDS:
        return jsonify(cached[1])

    # Pull parlay sizing settings from dashboard DB
    dcon = duckdb.connect(str(DB_PATH))
    try:
        bankroll_row = dcon.execute(
            "SELECT value FROM sizing_settings WHERE param = 'parlay_bankroll'"
        ).fetchone()
        kmult_row = dcon.execute(
            "SELECT value FROM sizing_settings WHERE param = 'parlay_kelly_mult'"
        ).fetchone()
    finally:
        dcon.close()
    if not bankroll_row or not kmult_row:
        return jsonify({"success": False, "error": "Parlay sizing settings missing"}), 500
    bankroll = bankroll_row[0]
    kmult = kmult_row[0]

    # Build legs for WZ pricing
    a_total_play = _total_play_from_combo(a[8])
    b_total_play = _total_play_from_combo(b[8])
    legs = [
        {"idgm": a[3], "play": _spread_play_from_combo(a[8]), "points": str(a[4]),                                 "odds": int(a[6])},
        {"idgm": a[3], "play": a_total_play,                  "points": _total_points_for_play(a[5], a_total_play), "odds": int(a[7])},
        {"idgm": b[3], "play": _spread_play_from_combo(b[8]), "points": str(b[4]),                                 "odds": int(b[6])},
        {"idgm": b[3], "play": b_total_play,                  "points": _total_points_for_play(b[5], b_total_play), "odds": int(b[7])},
    ]

    # Live WZ price (10000 first, fallback to 100)
    try:
        wz_session = _get_wz_session()
    except RuntimeError as e:
        return jsonify({"success": False, "error": f"Wagerzon auth failed: {e}"}), 502
    wz = wz_get_combined_parlay_price(wz_session, legs, amount=10000)
    if wz is None:
        wz = wz_get_combined_parlay_price(wz_session, legs, amount=100)
    if wz is None:
        return jsonify({"success": False, "error": "WZ pricing unavailable"}), 502

    pricing = joint_pricing(
        fair_dec_a=float(a[1]), fair_dec_b=float(b[1]), wz_dec=float(wz["decimal"]),
        bankroll=bankroll, kelly_mult=kmult,
    )

    response = {
        "success": True,
        "joint_fair_dec": pricing["joint_fair_dec"],
        "joint_fair_prob": pricing["joint_fair_prob"],
        "joint_edge": pricing["joint_edge"],
        "kelly_stake": pricing["kelly_stake"],
        "wz_dec": wz["decimal"],
        "wz_win": wz["win"],
        "amount": wz["amount"],
        "leg_a_combo": a[8],
        "leg_b_combo": b[8],
    }
    _COMBO_PRICE_CACHE[cache_key] = (now, response)
    return jsonify(response)


# NOTE: main previously defined a manual /api/place-parlay endpoint here that
# took game_id/home_team/away_team/combo/wz_odds/kelly_bet/etc. and just
# recorded a "pending" row in placed_parlays. That route is now superseded
# by api_place_parlay() below, which actually places the bet at Wagerzon
# via the REST API and records the resulting ticket_number. The manual
# route is intentionally NOT carried forward to avoid a Flask duplicate-
# endpoint error and to keep the placement contract single.


def _load_parlay_row(parlay_hash: str) -> dict | None:
    """Load one parlay opportunity row from mlb.duckdb (read-only)."""
    con = duckdb.connect(str(MLB_DB), read_only=True)
    try:
        row = con.execute(
            "SELECT * FROM mlb_parlay_opportunities WHERE parlay_hash = ?",
            [parlay_hash],
        ).fetchone()
        if not row:
            return None
        cols = [d[0] for d in con.description]
        return dict(zip(cols, row))
    finally:
        con.close()


def _build_spec_from_row(row: dict) -> parlay_placer.ParlaySpec:
    """Translate a mlb_parlay_opportunities row into a ParlaySpec.

    The combo column encodes both leg directions, e.g.:
      'Away Spread + Over'  → spread play=0 (away), total play=2 (over)
      'Home Spread + Over'  → spread play=1 (home), total play=2 (over)
      'Home Spread + Under' → spread play=1 (home), total play=3 (under)
      'Away Spread + Under' → spread play=0 (away), total play=3 (under)

    spread_line is already signed correctly in the DB for the picked side
    (favorite negative, dog positive) and is passed through as-is.

    total_line is stored positive in the DB, but Wagerzon's sel encoding
    requires NEGATIVE points for Over and POSITIVE for Under. This matches
    the existing convention in wagerzon_odds/parlay_pricer.py:331
    ("WZ expects -total for over, +total for under (recon-confirmed)").
    """
    idgm = int(row["idgm"])
    combo = row.get("combo", "")

    spread_play = 0 if "Away Spread" in combo else 1
    is_over = "Over" in combo
    total_play = 2 if is_over else 3
    total_line = float(row["total_line"])
    total_points = -abs(total_line) if is_over else total_line

    legs = [
        parlay_placer.Leg(
            idgm=idgm,
            play=spread_play,
            points=float(row["spread_line"]),
            odds=int(row["spread_price"]),
        ),
        parlay_placer.Leg(
            idgm=idgm,
            play=total_play,
            points=total_points,
            odds=int(row["total_price"]),
        ),
    ]
    amount, expected_win = _resolve_amount_and_win(row)
    return parlay_placer.ParlaySpec(
        parlay_hash=row["parlay_hash"],
        legs=legs,
        amount=amount,
        expected_win=expected_win,
        expected_risk=amount,
    )


def _resolve_amount_and_win(row: dict) -> tuple[float, float]:
    """Pick the placement stake AND the expected win in one place.

    The pricer (wagerzon_odds/parlay_pricer.py --exact-payouts) calls
    Wagerzon's ConfirmWagerHelper at an integer stake near Kelly and
    stores both (exact_wager, exact_to_win). The dashboard's SIZE / TO
    WIN columns show those numbers — so the user's mental "this is what
    I'm placing" matches the pricer's pre-confirmed pair.

    Source of truth: exact_wager + exact_to_win. They came from the same
    ConfirmWagerHelper call so they match each other to the cent and will
    match the placer's preflight again unless the line genuinely moves.

    Fallback (rare — the pricer hasn't run --exact-payouts yet for this
    row): use kelly_bet rounded to integer, expected_win from wz_odds
    math. The drift check still catches real moves; the only cost is
    a few cents of slack from Wagerzon's per-stake rounding.
    """
    exact_wager = row.get("exact_wager")
    exact_to_win = row.get("exact_to_win")

    if exact_wager is not None and exact_to_win is not None:
        return float(exact_wager), float(exact_to_win)

    # Pricer's exact-payouts step hasn't populated this row yet.
    amount = float(int(round(float(row["kelly_bet"]))))
    wz = int(row["wz_odds"])
    decimal = (wz / 100 + 1) if wz > 0 else (100 / -wz + 1)
    return amount, round(amount * (decimal - 1), 2)


def _upsert_placed_parlay(result: parlay_placer.PlacementResult, row: dict) -> None:
    """Insert or update placed_parlays in mlb_dashboard.duckdb. Idempotent on parlay_hash."""
    spec = _build_spec_from_row(row)
    legs_json = _json.dumps([
        {"play": l.play, "idgm": l.idgm, "points": l.points, "odds": l.odds}
        for l in spec.legs
    ])
    con = duckdb.connect(str(DASHBOARD_DB))
    try:
        con.execute("""
            INSERT INTO placed_parlays (
                parlay_hash, status, combo, game_id, game_time,
                recommended_size, expected_odds, expected_win,
                actual_size, actual_win, ticket_number, idwt,
                legs_json, error_msg, error_msg_key, updated_at
            ) VALUES (?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, CURRENT_TIMESTAMP)
            ON CONFLICT (parlay_hash) DO UPDATE SET
                status        = EXCLUDED.status,
                actual_size   = EXCLUDED.actual_size,
                actual_win    = EXCLUDED.actual_win,
                ticket_number = EXCLUDED.ticket_number,
                idwt          = EXCLUDED.idwt,
                error_msg     = EXCLUDED.error_msg,
                error_msg_key = EXCLUDED.error_msg_key,
                updated_at    = CURRENT_TIMESTAMP
        """, [
            result.parlay_hash,
            result.status,
            row.get("combo"),
            row.get("game_id"),
            row.get("game_time"),
            float(row["kelly_bet"]),
            int(row["wz_odds"]),
            float(row.get("exact_to_win") or 0) or round(float(row["kelly_bet"]) * (
                (int(row["wz_odds"]) / 100) if int(row["wz_odds"]) > 0
                else (100 / -int(row["wz_odds"]))
            ), 2),
            result.actual_risk,
            result.actual_win,
            result.ticket_number,
            result.idwt,
            legs_json,
            result.error_msg,
            result.error_msg_key,
        ])
    finally:
        con.close()


def _record_orphan(result: parlay_placer.PlacementResult,
                   parlay_hash: str, db_error: Exception) -> None:
    """Wagerzon confirmed placement but local placed_parlays write failed.
    Insert a forensics row into placement_orphans + log loudly."""
    try:
        con = duckdb.connect(str(DASHBOARD_DB))
        try:
            con.execute("""
                INSERT INTO placement_orphans
                    (idwt, ticket_number, parlay_hash, raw_response, error)
                VALUES (?, ?, ?, ?, ?)
                ON CONFLICT (idwt) DO NOTHING
            """, [result.idwt, result.ticket_number, parlay_hash,
                  result.raw_response, str(db_error)])
        finally:
            con.close()
    except Exception as inner:
        # Last resort: shout to stderr
        print(
            f"!! ORPHAN UNRECORDED: idwt={result.idwt} "
            f"ticket={result.ticket_number} parlay_hash={parlay_hash} "
            f"raw_response={result.raw_response!r} "
            f"original_error={db_error!r} orphan_write_error={inner!r}",
            file=sys.stderr, flush=True,
        )
        raise
    print(
        f"!! ORPHAN RECORDED: idwt={result.idwt} ticket={result.ticket_number} "
        f"parlay_hash={parlay_hash} reason={db_error!r}",
        file=sys.stderr, flush=True,
    )


@app.route("/api/place-parlay", methods=["POST"])
def api_place_parlay():
    """Auto-place a parlay at Wagerzon via the REST API (live or dry-run mode).

    Accepts JSON: {parlay_hash: str, dry_run: bool = false}
    Returns JSON: {status, ticket_number, error_msg, actual_win}
    """
    payload = request.get_json(force=True) or {}
    parlay_hash = payload.get("parlay_hash")
    dry_run = bool(payload.get("dry_run", False))

    if not parlay_hash:
        return jsonify({"status": "error", "error_msg": "missing parlay_hash"}), 400

    # Idempotency check: skip if already in a terminal placed/placing state (live only)
    if not dry_run:
        con = duckdb.connect(str(DASHBOARD_DB), read_only=True)
        try:
            existing = con.execute(
                "SELECT status, ticket_number, error_msg FROM placed_parlays "
                "WHERE parlay_hash = ? AND status IN ('placing', 'placed')",
                [parlay_hash],
            ).fetchone()
        finally:
            con.close()
        if existing:
            return jsonify({
                "status": existing[0],
                "ticket_number": existing[1],
                "error_msg": existing[2] or "",
            })

    row = _load_parlay_row(parlay_hash)
    if not row:
        return jsonify({"status": "error", "error_msg": "parlay not found"}), 404

    # Mark as 'placing' before network call so a crash leaves a breadcrumb
    if not dry_run:
        con = duckdb.connect(str(DASHBOARD_DB))
        try:
            con.execute("""
                INSERT INTO placed_parlays (
                    parlay_hash, status, recommended_size, expected_odds,
                    combo, game_id, game_time, updated_at
                ) VALUES (?, 'placing', ?, ?, ?, ?, ?, CURRENT_TIMESTAMP)
                ON CONFLICT (parlay_hash) DO NOTHING
            """, [
                parlay_hash,
                float(row["kelly_bet"]),
                int(row["wz_odds"]),
                row.get("combo"),
                row.get("game_id"),
                row.get("game_time"),
            ])
        finally:
            con.close()

    spec = _build_spec_from_row(row)
    results = parlay_placer.place_parlays([spec], dry_run=dry_run)
    result = results[0]

    # Dry-run: skip DB write entirely, return what would have happened
    if dry_run:
        return jsonify({
            "status": result.status,
            "short_label": _short_label_for(result.status, result.error_msg_key, result.error_msg),
            "ticket_number": None,
            "error_msg": result.error_msg,
            "actual_win": result.actual_win,
        })

    # Persist result; handle orphan case (placed at Wagerzon but local write failed)
    try:
        _upsert_placed_parlay(result, row)
    except Exception as e:
        if result.status == "placed":
            _record_orphan(result, parlay_hash, e)
            return jsonify({
                "status": "orphaned",
                "short_label": _short_label_for("orphaned", "", f"orphan: {e}"),
                "ticket_number": result.ticket_number,
                "error_msg": f"orphan: {e}",
            })
        raise

    return jsonify({
        "status": result.status,
        "short_label": _short_label_for(result.status, result.error_msg_key, result.error_msg),
        "ticket_number": result.ticket_number,
        "error_msg": result.error_msg,
    })


# Compact labels for the dashboard's status pill. Full text remains in
# error_msg (toast + hover tooltip). Keep all under ~15 chars.
_SHORT_LABEL_BY_STATUS = {
    "placed":         "placed",
    "would_place":    "would place",
    "price_moved":    "price moved",
    "auth_error":     "auth fail",
    "network_error":  "network err",
    "orphaned":       "orphaned",
}
_SHORT_LABEL_BY_REJECT_KEY = {
    "insufficient_funds": "insufficient $",
    "bet_too_large":      "exceeds limit",
    "line_unavailable":   "line pulled",
}


def _short_label_for(status: str, error_msg_key: str, error_msg: str) -> str:
    if status == "rejected":
        return _SHORT_LABEL_BY_REJECT_KEY.get(error_msg_key or "", "rejected")
    return _SHORT_LABEL_BY_STATUS.get(status, status)


@app.route("/api/place-combined-parlay", methods=["POST"])
def place_combined_parlay():
    """Record a combined (cross-game) parlay placement.

    Inserts a single row into placed_parlays with is_combo=TRUE and
    combo_leg_ids = JSON list of the two source parlay_hashes.
    """
    data = request.get_json(silent=True) or {}
    required = ["combo_hash", "parlay_hash_a", "parlay_hash_b",
                "wz_odds", "kelly_bet", "actual_size", "combo_label"]
    missing = [k for k in required if k not in data]
    if missing:
        return jsonify({"success": False, "error": f"Missing fields: {missing}"}), 400

    leg_ids_json = json.dumps([data["parlay_hash_a"], data["parlay_hash_b"]])

    try:
        con = duckdb.connect(str(DB_PATH))
        try:
            existing = con.execute(
                "SELECT parlay_hash, status FROM placed_parlays WHERE parlay_hash = ?",
                [data["combo_hash"]]
            ).fetchone()
            if existing and existing[1] == "pending":
                return jsonify({"success": False, "error": "Combo already placed"}), 409

            con.execute("""
                INSERT INTO placed_parlays (
                    parlay_hash, game_id, home_team, away_team, combo,
                    wz_odds, kelly_bet, actual_size, placed_at, status,
                    is_combo, combo_leg_ids
                ) VALUES (?, ?, ?, ?, ?, ?, ?, ?, ?, 'pending', TRUE, ?)
            """, [
                data["combo_hash"],
                "COMBO",  # synthetic game_id for combo rows
                "(combined)", "(combined)",
                data["combo_label"],
                int(data["wz_odds"]),
                float(data["kelly_bet"]),
                float(data["actual_size"]),
                datetime.now().isoformat(),
                leg_ids_json,
            ])
        finally:
            con.close()
        return jsonify({"success": True, "message": "Combined parlay placed"})
    except Exception as e:
        return jsonify({"success": False, "error": str(e)}), 500


@app.route("/api/remove-parlay", methods=["POST"])
def remove_parlay():
    """Remove a parlay from placed status."""
    data = request.json
    parlay_hash = data.get("parlay_hash")

    if not parlay_hash:
        return jsonify({"success": False, "error": "parlay_hash required"}), 400

    try:
        con = duckdb.connect(str(DB_PATH))
        try:
            result = con.execute(
                "DELETE FROM placed_parlays WHERE parlay_hash = ? RETURNING parlay_hash",
                [parlay_hash]
            ).fetchone()
        finally:
            con.close()

        if result:
            return jsonify({"success": True, "message": "Parlay removed"})
        else:
            return jsonify({"success": False, "error": "Parlay not found"}), 404

    except Exception as e:
        return jsonify({"success": False, "error": str(e)}), 500


@app.route("/api/exposure", methods=["GET"])
def get_exposure():
    """Get current exposure summary by game."""
    try:
        con = duckdb.connect(str(DB_PATH))
        try:
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
        finally:
            con.close()
        return jsonify(exposure.to_dict(orient="records"))

    except Exception as e:
        return jsonify({"success": False, "error": str(e)}), 500


@app.route("/api/book-settings", methods=["GET"])
def get_book_settings():
    """Get all bookmaker enabled/disabled settings."""
    try:
        con = duckdb.connect(str(DB_PATH))
        try:
            rows = con.execute("SELECT bookmaker_key, enabled FROM book_settings").fetchall()
        finally:
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
        try:
            con.execute("""
                INSERT INTO book_settings (bookmaker_key, enabled)
                VALUES (?, ?)
                ON CONFLICT (bookmaker_key) DO UPDATE SET enabled = ?
            """, [book, enabled, enabled])
        finally:
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
        try:
            for book in books:
                con.execute("""
                    INSERT INTO book_settings (bookmaker_key, enabled)
                    VALUES (?, FALSE)
                    ON CONFLICT (bookmaker_key) DO NOTHING
                """, [book])
        finally:
            con.close()
        return jsonify({"success": True, "added": len(books)})
    except Exception as e:
        return jsonify({"success": False, "error": str(e)}), 500


@app.route("/api/sizing-settings", methods=["GET"])
def get_sizing_settings():
    """Get bankroll and kelly multiplier settings."""
    try:
        con = duckdb.connect(str(DB_PATH))
        try:
            rows = con.execute("SELECT param, value FROM sizing_settings").fetchall()
        finally:
            con.close()
        return jsonify({row[0]: row[1] for row in rows})
    except Exception as e:
        return jsonify({"bankroll": 100, "kelly_mult": 0.25}), 200


@app.route("/api/sizing-settings", methods=["POST"])
def update_sizing_settings():
    """Update bankroll and/or kelly multiplier settings."""
    data = request.json
    allowed = ("bankroll", "kelly_mult", "parlay_bankroll", "parlay_kelly_mult", "parlay_min_edge")
    try:
        con = duckdb.connect(str(DB_PATH))
        try:
            # Support both {param: "name", value: X} and {bankroll: X, kelly_mult: Y} formats
            if "param" in data and "value" in data and data["param"] in allowed:
                con.execute("""
                    INSERT INTO sizing_settings (param, value) VALUES (?, ?)
                    ON CONFLICT (param) DO UPDATE SET value = ?
                """, [data["param"], float(data["value"]), float(data["value"])])
            else:
                for param in allowed:
                    if param in data:
                        con.execute("""
                            INSERT INTO sizing_settings (param, value) VALUES (?, ?)
                            ON CONFLICT (param) DO UPDATE SET value = ?
                        """, [param, float(data[param]), float(data[param])])
        finally:
            con.close()
        return jsonify({"success": True})
    except Exception as e:
        return jsonify({"success": False, "error": str(e)}), 500


@app.route("/api/filter-settings", methods=["GET"])
def get_filter_settings():
    """Get persisted filter states (market, correlation, status)."""
    try:
        con = duckdb.connect(str(DB_PATH))
        try:
            rows = con.execute("SELECT filter_type, selected_values FROM filter_settings").fetchall()
        finally:
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

    if not filter_type or filter_type not in ("market", "correlation", "size", "status", "parlay_status", "parlay_size"):
        return jsonify({"success": False, "error": "Invalid filter_type"}), 400

    try:
        con = duckdb.connect(str(DB_PATH))
        try:
            con.execute("""
                INSERT INTO filter_settings (filter_type, selected_values)
                VALUES (?, ?)
                ON CONFLICT (filter_type) DO UPDATE SET selected_values = ?
            """, [filter_type, json.dumps(selected_values), json.dumps(selected_values)])
        finally:
            con.close()
        return jsonify({"success": True})
    except Exception as e:
        return jsonify({"success": False, "error": str(e)}), 500


@app.route("/api/auto-place", methods=["POST"])
def auto_place():
    """Launch Playwright to navigate to book and pre-fill bet.

    Spawns placer.py as a non-blocking subprocess with the bet data as JSON arg.
    The subprocess opens a visible browser, logs in, navigates to the game,
    and pre-fills the bet slip. User confirms manually on the book's site.
    """
    data = request.json
    bookmaker = data.get("bookmaker")

    if bookmaker not in SUPPORTED_AUTO_BOOKS:
        return jsonify({
            "success": False,
            "error": f"Auto-place not supported for '{bookmaker}'. Supported: {SUPPORTED_AUTO_BOOKS}"
        }), 400

    try:
        local_root = PROJECT_ROOT
        repo_root = _REPO_ROOT

        placer_script = str(local_root / "bet_placer" / "placer.py")
        bet_json = json.dumps(data)

        # Use wagerzon_odds venv (has playwright + duckdb + dotenv)
        placer_python = str(repo_root / "wagerzon_odds" / "venv" / "bin" / "python3")
        if not Path(placer_python).exists():
            placer_python = sys.executable  # fallback

        log_path = str(local_root / "bet_placer" / "placer.log")
        log_file = open(log_path, "a")
        subprocess.Popen(
            [placer_python, placer_script, bet_json],
            cwd=str(local_root / "bet_placer"),
            stdout=log_file,
            stderr=log_file,
            start_new_session=True,
        )
        log_file.close()

        return jsonify({"success": True, "message": f"Browser launching for {bookmaker}..."})

    except Exception as e:
        return jsonify({"success": False, "error": str(e)}), 500


@app.route("/api/nav-status/<bet_hash>", methods=["GET"])
def nav_status(bet_hash):
    """Get current navigator status for a bet (used by dashboard polling)."""
    try:
        con = duckdb.connect(str(DB_PATH), read_only=True)
        try:
            row = con.execute(
                "SELECT status FROM placed_bets WHERE bet_hash = ?", [bet_hash]
            ).fetchone()
        finally:
            con.close()
        if row:
            return jsonify({"status": row[0]})
        return jsonify({"status": "unknown"}), 404
    except Exception as e:
        return jsonify({"status": "error", "error": str(e)}), 500


# =============================================================================
# CLV API ENDPOINTS
# =============================================================================

@app.route("/api/clv-summary", methods=["GET"])
def clv_summary():
    """Aggregate CLV stats: overall, by market type, and by book."""
    try:
        con = duckdb.connect(str(DB_PATH), read_only=True)
        try:
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
        finally:
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
        try:
            rows = con.execute("""
                SELECT c.*, p.home_team, p.away_team, p.model_ev
                FROM bet_clv c
                JOIN placed_bets p ON c.bet_hash = p.bet_hash
                ORDER BY c.game_time DESC
            """).fetchdf()
        finally:
            con.close()

        return jsonify(rows.to_dict(orient="records"))
    except Exception as e:
        return jsonify({"success": False, "error": str(e)}), 500


@app.route("/api/clv-compute", methods=["POST"])
def trigger_clv_compute():
    """Trigger CLV computation for completed games."""
    try:
        script = BASE_DIR.parent / "MLB Answer Key" / "clv_compute.py"
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
    if not _refresh_lock.acquire(blocking=False):
        return False, "Refresh already in progress — please wait"

    try:
        nfl_work_dir = BASE_DIR.parent.parent  # NFLWork directory
        answer_keys_dir = BASE_DIR.parent  # Answer Keys directory

        # Step 1: Run the full pipeline (scrapers + R prepare + R combine)
        print("Running full MLB pipeline (scrapers + predictions)...")
        result = subprocess.run(
            [sys.executable, str(answer_keys_dir / "run.py"), "mlb"],
            capture_output=True,
            text=True,
            cwd=nfl_work_dir
        )

        if result.returncode != 0:
            # run.py prints errors to stdout, not stderr — check both
            error_output = (result.stderr or result.stdout or "Unknown error")[-500:]
            print(f"Pipeline failed: {error_output}")
            return False, f"Pipeline failed: {error_output}"

        # Step 1.5: Price parlays via Wagerzon ConfirmWagerHelper (non-fatal)
        # Must run AFTER scrapers (reads mlb_odds) and BEFORE correlated parlay
        # finder (reads mlb_parlay_prices). Without this, stale exact prices from
        # a previous manual run get used, causing WZ odds / payout mismatches.
        print("Pricing parlays on Wagerzon...")
        pricing_result = subprocess.run(
            [sys.executable, str(nfl_work_dir / "wagerzon_odds" / "parlay_pricer.py"), "mlb"],
            capture_output=True,
            text=True,
            cwd=str(nfl_work_dir)
        )
        if pricing_result.returncode != 0:
            print(f"Parlay pricing warning: {(pricing_result.stderr or pricing_result.stdout or '')[-300:]}")

        # Step 2: Find correlated parlay opportunities (non-fatal)
        print("Finding parlay opportunities...")
        parlay_result = subprocess.run(
            ["Rscript", str(answer_keys_dir / "mlb_correlated_parlay.R")],
            capture_output=True,
            text=True,
            cwd=str(nfl_work_dir)
        )
        if parlay_result.returncode != 0:
            print(f"Parlay finder warning: {(parlay_result.stderr or '')[-300:]}")

        # Step 2.5: Empirical nudge + exact payout per sized parlay (non-fatal).
        # Queries ConfirmWagerHelper at stake ± NUDGE_RANGE around each Kelly-ideal
        # wager and stores the best (exact_wager, exact_to_win) on the opportunities
        # row. Dashboard renders those directly so "To Win" matches the WZ slip.
        print("Computing exact payouts at Kelly stakes...")
        exact_result = subprocess.run(
            [sys.executable, str(nfl_work_dir / "wagerzon_odds" / "parlay_pricer.py"),
             "mlb", "--exact-payouts"],
            capture_output=True,
            text=True,
            cwd=str(nfl_work_dir)
        )
        if exact_result.returncode != 0:
            print(f"Exact payout warning: "
                  f"{(exact_result.stderr or exact_result.stdout or '')[-300:]}")

        # Step 3: Generate dashboard HTML from the saved data
        print("Generating MLB dashboard HTML...")
        result = subprocess.run(
            ["Rscript", str(BASE_DIR / "mlb_dashboard.R")],
            capture_output=True,
            text=True,
            cwd=nfl_work_dir
        )

        if result.returncode != 0:
            error_output = (result.stderr or result.stdout or "Unknown error")[-500:]
            print(f"Dashboard failed: {error_output}")
            return False, f"Dashboard generation failed: {error_output}"

        print("MLB refresh complete!")
        return True, "Dashboard refreshed successfully"
    finally:
        _refresh_lock.release()


@app.route("/refresh", methods=["POST"])
def refresh():
    """Run full pipeline and regenerate dashboard."""
    try:
        success, message = run_pipeline()
        if not success:
            return jsonify({"success": False, "error": message}), 500

        return jsonify({
            "success": True,
            "message": message,
        })
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
    print("MLB +EV Betting Dashboard")
    print("=" * 50)
    print(f"\nOpen in browser: http://localhost:8083")
    print("Press Ctrl+C to stop\n")
    app.run(host="0.0.0.0", port=8083, debug=False)
