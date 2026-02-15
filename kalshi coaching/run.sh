#!/bin/bash
# Kalshi NFL Coaching Odds - Fetch and Display
#
# Usage:
#   ./run.sh          - Fetch data, generate report, start server
#   ./run.sh --static - Just generate the HTML (no server)

set -e  # Exit on error

SCRIPT_DIR="$(cd "$(dirname "$0")" && pwd)"
PROJECT_DIR="$(dirname "$SCRIPT_DIR")"

cd "$PROJECT_DIR"

echo "=================================="
echo "Kalshi NFL Coaching Odds"
echo "=================================="
echo ""

# Step 1: Fetch latest data
echo "Step 1: Fetching latest odds from Kalshi..."
python3 "$SCRIPT_DIR/kalshi_coaching.py"
echo ""

# Step 2: Generate visualizations
echo "Step 2: Generating heatmap and table..."
Rscript "$SCRIPT_DIR/kalshi_coaching_display.R"
echo ""

# Check for --static flag
if [ "$1" = "--static" ]; then
    echo "Step 3: Opening static dashboard..."
    open "$SCRIPT_DIR/report.html"
    echo ""
    echo "Done! (static mode - refresh button won't work)"
else
    # Step 3: Start server
    echo "Step 3: Starting dashboard server..."
    echo ""

    # Install flask if needed
    pip3 install flask -q 2>/dev/null || true

    # Open browser after short delay
    (sleep 1 && open "http://127.0.0.1:8080") &

    # Start server
    python3 "$SCRIPT_DIR/server.py"
fi
