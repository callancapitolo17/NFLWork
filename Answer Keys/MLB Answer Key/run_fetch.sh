#!/bin/zsh
export PATH="/opt/homebrew/bin:/usr/local/bin:/usr/bin:/bin"
export TZ="America/Los_Angeles"
export R_ENVIRON_USER="/Users/callancapitolo/.Renviron"

LOGDIR="/Users/callancapitolo/Library/Logs/fetch-mlb-odds"
mkdir -p "$LOGDIR"
STAMP=$(date +"%Y-%m-%d_%H-%M-%S")

# --- Friday-only + once-per-day gate (safe to keep forever) ---
MARKER="$LOGDIR/last_run_date.txt"
TODAY=$(date +%F)
WEEKDAY=$(date +%u)  # 5 = Friday
if [ "${FORCE_RUN:-0}" != "1" ]; then
  if [ "$WEEKDAY" != "5" ]; then
    exit 0
  fi
  if [ -f "$MARKER" ] && [ "$(cat "$MARKER" 2>/dev/null)" = "$TODAY" ]; then
    exit 0
  fi
fi
# ---------------------------------------------------------------

/usr/bin/env Rscript "/Users/callancapitolo/NFLWork/Answer Keys/MLB Answer Key/Acquire New MLB Data.R" \
  >> "$LOGDIR/run_$STAMP.out" 2>> "$LOGDIR/run_$STAMP.err"

STATUS=$?
echo "$TODAY" > "$MARKER"

DESKTOP_FILE="/Users/callancapitolo/Desktop/Fetch MLB Odds - last run.txt"
if [ $STATUS -eq 0 ]; then
  echo "$(date) — finished successfully" > "$DESKTOP_FILE"
  open "$DESKTOP_FILE"
else
  echo "$(date) — FAILED (see logs at $LOGDIR)" > "$DESKTOP_FILE"
  open "$DESKTOP_FILE"
fi

exit $STATUS
