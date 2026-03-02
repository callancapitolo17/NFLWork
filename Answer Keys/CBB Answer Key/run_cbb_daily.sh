#!/bin/zsh
# run_cbb_daily.sh - Daily CBB historical data acquisition
# Called by launchd at 5 AM daily. Each step dynamically fills
# from the last acquired date to yesterday (no hardcoded lookback).

export PATH="/opt/homebrew/bin:/usr/local/bin:/usr/bin:/bin"
export TZ="America/Los_Angeles"
export R_ENVIRON_USER="/Users/callancapitolo/.Renviron"

SCRIPT_DIR="$(cd "$(dirname "$0")" && pwd)"
ANSWER_KEYS_DIR="$(dirname "$SCRIPT_DIR")"
LOGDIR="/Users/callancapitolo/Library/Logs/cbb-daily-acquire"
mkdir -p "$LOGDIR"

STAMP=$(date +"%Y-%m-%d_%H-%M-%S")
LOGFILE="$LOGDIR/run_${STAMP}.log"

echo "$(date) - Starting CBB daily acquisition" >> "$LOGFILE"

# Step 1: Acquire closing odds (dynamic gap fill)
echo "$(date) - [Step 1/3] Acquiring closing odds..." >> "$LOGFILE"
cd "$ANSWER_KEYS_DIR"
Rscript "$SCRIPT_DIR/Acquire CBB Data.R" --daily >> "$LOGFILE" 2>&1
STEP1=$?
echo "$(date) - [Step 1/3] Exit code: $STEP1" >> "$LOGFILE"

# Step 2: Acquire PBP scores (dynamic gap fill)
echo "$(date) - [Step 2/3] Acquiring PBP scores..." >> "$LOGFILE"
cd "$SCRIPT_DIR"
python3 "$SCRIPT_DIR/acquire_cbb_pbp.py" --daily >> "$LOGFILE" 2>&1
STEP2=$?
echo "$(date) - [Step 2/3] Exit code: $STEP2" >> "$LOGFILE"

# Step 3: Rebuild betting_pbp join table
echo "$(date) - [Step 3/3] Building betting_pbp table..." >> "$LOGFILE"
cd "$ANSWER_KEYS_DIR"
Rscript "$SCRIPT_DIR/build_betting_pbp.R" >> "$LOGFILE" 2>&1
STEP3=$?
echo "$(date) - [Step 3/3] Exit code: $STEP3" >> "$LOGFILE"

# Summary
if [ $STEP1 -eq 0 ] && [ $STEP2 -eq 0 ] && [ $STEP3 -eq 0 ]; then
  echo "$(date) - All steps completed successfully" >> "$LOGFILE"
  exit 0
else
  echo "$(date) - FAILED: odds=$STEP1 pbp=$STEP2 build=$STEP3" >> "$LOGFILE"
  exit 1
fi
