#!/bin/bash
# NFL +EV Dashboard Launcher
# Starts the Flask server and opens the dashboard in browser

cd "$(dirname "$0")"
SCRIPT_DIR=$(pwd)

echo "=================================================="
echo "NFL +EV Dashboard"
echo "=================================================="
echo ""

# Activate virtual environment
if [ -d "$SCRIPT_DIR/venv" ]; then
    source "$SCRIPT_DIR/venv/bin/activate"
else
    echo "Creating virtual environment..."
    python3 -m venv "$SCRIPT_DIR/venv"
    source "$SCRIPT_DIR/venv/bin/activate"
    pip install flask duckdb
fi

cd ../..  # Go to NFLWork directory

# Generate initial dashboard if report.html doesn't exist
if [ ! -f "$SCRIPT_DIR/report.html" ]; then
    echo "Generating initial dashboard..."
    Rscript "$SCRIPT_DIR/nfl_dashboard.R"
    if [ $? -ne 0 ]; then
        echo "Warning: Failed to generate initial dashboard."
        echo "The server will start anyway - use the Refresh button to generate."
    fi
fi

echo ""
echo "Starting server..."
echo "Dashboard will be available at: http://127.0.0.1:8081"
echo ""
echo "Press Ctrl+C to stop the server"
echo ""

# Open browser after a short delay (in background)
(sleep 2 && open "http://127.0.0.1:8081") &

# Start the server
python3 "$SCRIPT_DIR/nfl_dashboard_server.py"
