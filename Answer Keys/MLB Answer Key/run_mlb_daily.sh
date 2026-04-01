#!/bin/bash
# MLB Daily Data Acquisition
# Runs Acquire New MLB Data.R to fetch closing odds and PBP data.
# Schedule via launchd during MLB season (April-October).

set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "$0")" && pwd)"
LOG_DIR="$HOME/Library/Logs/mlb-daily-acquire"
mkdir -p "$LOG_DIR"

LOGFILE="$LOG_DIR/run_$(date +%Y-%m-%d_%H-%M-%S).log"

# Rotate logs older than 30 days
find "$LOG_DIR" -name "run_*.log" -mtime +30 -delete 2>/dev/null || true

echo "=== MLB Daily Acquire $(date) ===" | tee "$LOGFILE"

cd "$SCRIPT_DIR"

# Step 1: Acquire closing odds from Odds API
echo "Step 1: Acquiring closing odds..." | tee -a "$LOGFILE"
Rscript "Acquire New MLB Data.R" --daily >> "$LOGFILE" 2>&1

# Step 2: Acquire PBP data from baseballr
echo "Step 2: Acquiring PBP data..." | tee -a "$LOGFILE"
Rscript "Acquire New MLB Data.R" --daily-pbp >> "$LOGFILE" 2>&1

echo "=== Done $(date) ===" | tee -a "$LOGFILE"

# Desktop notification
echo "MLB daily acquire completed at $(date)" > "$HOME/Desktop/MLB Daily Acquire - last run.txt"
