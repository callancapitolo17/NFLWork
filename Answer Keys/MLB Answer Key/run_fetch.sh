#!/bin/zsh
# MLB Daily Data Acquisition — canonical entry point for launchd
# Runs three steps: (1) closing odds, (2) play-by-play, (3) consensus rebuild.
# Idempotent at the data layer — dedup prevents duplicate writes if re-run.
#
# Invoked by launchd via ~/Library/LaunchAgents/com.callan.fetch-mlb-odds.plist
# Manual invocation: bash run_fetch.sh

export PATH="/opt/homebrew/bin:/usr/local/bin:/usr/bin:/bin"
export TZ="America/Los_Angeles"
export R_ENVIRON_USER="/Users/callancapitolo/.Renviron"

SCRIPT_DIR="$(cd "$(dirname "$0")" && pwd)"
LOGDIR="/Users/callancapitolo/Library/Logs/fetch-mlb-odds"
mkdir -p "$LOGDIR"
STAMP=$(date +"%Y-%m-%d_%H-%M-%S")
LOGFILE="$LOGDIR/run_$STAMP.log"

# Rotate logs older than 30 days
find "$LOGDIR" -name "run_*.log" -mtime +30 -delete 2>/dev/null || true
find "$LOGDIR" -name "run_*.out" -mtime +30 -delete 2>/dev/null || true
find "$LOGDIR" -name "run_*.err" -mtime +30 -delete 2>/dev/null || true

cd "$SCRIPT_DIR"

echo "=== MLB Daily Acquire $(date) ===" | tee "$LOGFILE"

# Step 1: Acquire closing odds from Odds API
echo "Step 1: Acquiring closing odds..." | tee -a "$LOGFILE"
Rscript "Acquire New MLB Data.R" --daily >> "$LOGFILE" 2>&1
STEP1=$?

# Step 2: Acquire PBP data from baseballr
echo "Step 2: Acquiring PBP data..." | tee -a "$LOGFILE"
Rscript "Acquire New MLB Data.R" --daily-pbp >> "$LOGFILE" 2>&1
STEP2=$?

# Step 3: Rebuild mlb_betting_pbp (joins odds + PBP into the table MLB.R consumes)
echo "Step 3: Rebuilding betting PBP consensus..." | tee -a "$LOGFILE"
Rscript "Consensus Betting History.R" >> "$LOGFILE" 2>&1
STEP3=$?

echo "=== Done $(date) ===" | tee -a "$LOGFILE"

# Summary + macOS notification (only on failure or new data)
if [ $STEP1 -eq 0 ] && [ $STEP2 -eq 0 ] && [ $STEP3 -eq 0 ]; then
  # Check if new odds or games were actually added (R scripts print these on success)
  NEW_ODDS=$(grep -o "Inserted [0-9]* new odds rows" "$LOGFILE" | head -1)
  NEW_PBP=$(grep -o "Inserted [0-9]* PBP rows for [0-9]* games" "$LOGFILE" | head -1)

  if [ -n "$NEW_ODDS" ] || [ -n "$NEW_PBP" ]; then
    DETAILS=""
    [ -n "$NEW_ODDS" ] && DETAILS="$NEW_ODDS"
    [ -n "$NEW_PBP" ] && DETAILS="${DETAILS:+$DETAILS, }$NEW_PBP"
    echo "$(date) - Success with new data: $DETAILS" >> "$LOGFILE"
    osascript -e "display notification \"$DETAILS\" with title \"MLB Daily Acquire\""
  else
    echo "$(date) - All steps passed, no new data" >> "$LOGFILE"
  fi
  exit 0
else
  echo "$(date) - FAILED: odds=$STEP1 pbp=$STEP2 consensus=$STEP3" >> "$LOGFILE"
  osascript -e "display notification \"FAILED — odds=$STEP1 pbp=$STEP2 consensus=$STEP3\" with title \"MLB Daily Acquire\""
  exit 1
fi
