#!/bin/zsh
# run_college_baseball.sh - College baseball correlated parlay edge finder
# Runs scrapers (Wagerzon + Hoop88) + R answer key via run.py orchestrator.
# Called by launchd at scheduled intervals during baseball season (Feb-Jun).

export PATH="/opt/homebrew/bin:/usr/local/bin:/usr/bin:/bin"
export TZ="America/Los_Angeles"
export R_ENVIRON_USER="/Users/callancapitolo/.Renviron"

SCRIPT_DIR="$(cd "$(dirname "$0")" && pwd)"
ANSWER_KEYS_DIR="$(dirname "$SCRIPT_DIR")"
LOGDIR="/Users/callancapitolo/Library/Logs/college-baseball"
mkdir -p "$LOGDIR"

STAMP=$(date +"%Y-%m-%d_%H-%M-%S")
LOGFILE="$LOGDIR/run_${STAMP}.log"

echo "$(date) - Starting college baseball parlay edge finder" >> "$LOGFILE"

# Run pipeline via orchestrator (scrapers + R in parallel)
cd "$ANSWER_KEYS_DIR"
python3 run.py college_baseball >> "$LOGFILE" 2>&1
RC=$?

echo "$(date) - Pipeline exit code: $RC" >> "$LOGFILE"

# Clean up logs older than 30 days
find "$LOGDIR" -name "run_*.log" -mtime +30 -delete 2>/dev/null

# macOS notification (click opens log file)
if [ $RC -eq 0 ]; then
  EDGES=$(grep "+EV parlays" "$LOGFILE" | tail -1)
  terminal-notifier \
    -title "College Baseball AK" \
    -message "${EDGES:-Complete — no edges}" \
    -sound Glass \
    -execute "open -a Console '$LOGFILE'"
else
  terminal-notifier \
    -title "College Baseball AK" \
    -message "Pipeline failed — click to view log" \
    -sound Basso \
    -execute "open -a Console '$LOGFILE'"
fi

exit $RC
