#!/bin/bash

# Run All Bet Scrapers
# Runs all betting platform scrapers in sequence.
# Continues to the next scraper even if one fails.

cd "/Users/callancapitolo/NFLWork/bet_logger"

FAILED=0

echo "========================================"
echo "MULTI-PLATFORM BET SCRAPER"
echo "Started: $(date '+%Y-%m-%d %H:%M:%S')"
echo "========================================"
echo ""

# Run Wagerzon scraper
echo "[$(date '+%H:%M:%S')] Running Wagerzon scraper..."
echo "----------------------------------------"
if ./venv/bin/python3 scraper_wagerzon.py; then
    echo "[$(date '+%H:%M:%S')] Wagerzon: done"
else
    echo "[$(date '+%H:%M:%S')] Wagerzon: FAILED (exit $?)"
    FAILED=$((FAILED + 1))
fi
echo ""

# Run Hoop88 scraper
echo "[$(date '+%H:%M:%S')] Running Hoop88 scraper..."
echo "----------------------------------------"
if ./venv/bin/python3 scraper_hoop88.py; then
    echo "[$(date '+%H:%M:%S')] Hoop88: done"
else
    echo "[$(date '+%H:%M:%S')] Hoop88: FAILED (exit $?)"
    FAILED=$((FAILED + 1))
fi
echo ""

# Run BFA Gaming scraper (recon first to get fresh auth token)
echo "[$(date '+%H:%M:%S')] Running BFA Gaming recon + scraper..."
echo "----------------------------------------"
./venv/bin/python3 recon_bfa.py
if ./venv/bin/python3 scraper_bfa.py; then
    echo "[$(date '+%H:%M:%S')] BFA: done"
else
    echo "[$(date '+%H:%M:%S')] BFA: FAILED (exit $?)"
    FAILED=$((FAILED + 1))
fi
echo ""

# Run BetOnline scraper
echo "[$(date '+%H:%M:%S')] Running BetOnline scraper..."
echo "----------------------------------------"
if ./venv/bin/python3 scraper_betonline.py --since-last; then
    echo "[$(date '+%H:%M:%S')] BetOnline: done"
else
    echo "[$(date '+%H:%M:%S')] BetOnline: FAILED (exit $?)"
    FAILED=$((FAILED + 1))
fi
echo ""

echo "========================================"
if [ $FAILED -eq 0 ]; then
    MSG="All 4 scrapers completed successfully."
    echo "$MSG"
    osascript -e "display notification \"$MSG\" with title \"Bet Logger\" sound name \"Glass\""
else
    MSG="Completed with $FAILED scraper(s) failed. Check logs."
    echo "$MSG"
    osascript -e "display notification \"$MSG\" with title \"Bet Logger\" sound name \"Basso\""
fi
echo "Finished: $(date '+%Y-%m-%d %H:%M:%S')"
echo "========================================"

exit $FAILED
