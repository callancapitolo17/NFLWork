#!/bin/bash
# NFL Draft Dashboard - Entry point
# Usage: ./run.sh [--fetch-only]

set -e
cd "$(dirname "$0")"

VENV="./venv/bin/python"

# Kill existing server on port 8083
lsof -ti:8083 | xargs kill -9 2>/dev/null || true

echo "=== NFL Draft Dashboard ==="

# Step 1: Fetch fresh data
echo "[1/3] Fetching draft odds from Kalshi..."
$VENV fetcher.py

# Step 2: Run edge detection
echo "[2/3] Running edge detection..."
$VENV edge_detector.py

# Step 3: Scrape consensus
echo "[3/3] Scraping mock draft consensus..."
$VENV consensus.py

if [ "$1" = "--fetch-only" ]; then
    echo "Done (fetch-only mode)."
    exit 0
fi

# Step 4: Start Dash server
echo ""
echo "Starting dashboard on http://127.0.0.1:8083"
$VENV app.py &
SERVER_PID=$!

# Wait for server to start
sleep 2

# Open browser
open "http://127.0.0.1:8083" 2>/dev/null || true

echo "Dashboard running (PID: $SERVER_PID). Press Ctrl+C to stop."
wait $SERVER_PID
