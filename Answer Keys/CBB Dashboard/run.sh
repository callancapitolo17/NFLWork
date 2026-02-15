#!/bin/bash
# CBB +EV Dashboard Launcher
# Starts the Flask server and opens the dashboard in browser

cd "$(dirname "$0")"
SCRIPT_DIR=$(pwd)

echo "=================================================="
echo "CBB +EV Dashboard"
echo "=================================================="
echo ""

# Activate virtual environment (shared with NFL Dashboard if it exists)
NFL_VENV="$SCRIPT_DIR/../NFL Dashboard/venv"
if [ -d "$NFL_VENV" ]; then
    echo "Using shared venv from NFL Dashboard..."
    source "$NFL_VENV/bin/activate"
elif [ -d "$SCRIPT_DIR/venv" ]; then
    source "$SCRIPT_DIR/venv/bin/activate"
else
    echo "Creating virtual environment..."
    python3 -m venv "$SCRIPT_DIR/venv"
    source "$SCRIPT_DIR/venv/bin/activate"
    pip install flask duckdb pandas numpy
fi

cd ../..  # Go to NFLWork directory

# Kill any existing server on port 8082
lsof -ti:8082 | xargs kill 2>/dev/null
sleep 1

# Always run full pipeline on launch
echo "Running full CBB pipeline (scrapers + predictions)..."
python3 "$SCRIPT_DIR/../run.py" cbb
if [ $? -ne 0 ]; then
    echo "Warning: Pipeline had errors, but continuing..."
fi

echo "Generating dashboard HTML..."
Rscript "$SCRIPT_DIR/cbb_dashboard.R"
if [ $? -ne 0 ]; then
    echo "Warning: Failed to generate dashboard."
    echo "The server will start anyway - use the Refresh button to try again."
fi

echo ""
echo "Starting server..."
echo "Dashboard will be available at: http://127.0.0.1:8082"
echo ""
echo "Press Ctrl+C to stop the server"
echo ""

# Open browser after a short delay (in background)
(sleep 2 && open "http://127.0.0.1:8082") &

# Start the server
python3 "$SCRIPT_DIR/cbb_dashboard_server.py"
