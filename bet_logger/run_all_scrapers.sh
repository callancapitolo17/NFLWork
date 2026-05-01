#!/bin/bash

# Run All Bet Scrapers
# Runs all betting platform scrapers in sequence.
# Continues to the next scraper even if one fails.

cd "/Users/callancapitolo/NFLWork/bet_logger"

FAILED=0
FAILED_NAMES=""

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
    FAILED_NAMES="${FAILED_NAMES}Wagerzon, "
fi
echo ""

# Run Wagerzon scraper — WagerzonJ account
echo "[$(date '+%H:%M:%S')] Running Wagerzon scraper (WagerzonJ)..."
echo "----------------------------------------"
if ./venv/bin/python3 scraper_wagerzon.py --account j; then
    echo "[$(date '+%H:%M:%S')] WagerzonJ: done"
else
    echo "[$(date '+%H:%M:%S')] WagerzonJ: FAILED (exit $?)"
    FAILED=$((FAILED + 1))
    FAILED_NAMES="${FAILED_NAMES}WagerzonJ, "
fi
echo ""

# Run Wagerzon scraper — WagerzonC account
echo "[$(date '+%H:%M:%S')] Running Wagerzon scraper (WagerzonC)..."
echo "----------------------------------------"
if ./venv/bin/python3 scraper_wagerzon.py --account c; then
    echo "[$(date '+%H:%M:%S')] WagerzonC: done"
else
    echo "[$(date '+%H:%M:%S')] WagerzonC: FAILED (exit $?)"
    FAILED=$((FAILED + 1))
    FAILED_NAMES="${FAILED_NAMES}WagerzonC, "
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
    FAILED_NAMES="${FAILED_NAMES}Hoop88, "
fi
echo ""

# Run BFA Gaming scraper — primary account (recon first to get fresh auth token)
echo "[$(date '+%H:%M:%S')] Running BFA Gaming recon + scraper (primary)..."
echo "----------------------------------------"
./venv/bin/python3 recon_bfa.py
if ./venv/bin/python3 scraper_bfa.py; then
    echo "[$(date '+%H:%M:%S')] BFA primary: done"
else
    echo "[$(date '+%H:%M:%S')] BFA primary: FAILED (exit $?)"
    FAILED=$((FAILED + 1))
    FAILED_NAMES="${FAILED_NAMES}BFA, "
fi
echo ""

# Run BFA Gaming scraper — BFAJ account
echo "[$(date '+%H:%M:%S')] Running BFA Gaming recon + scraper (BFAJ)..."
echo "----------------------------------------"
./venv/bin/python3 recon_bfa.py --account j
if ./venv/bin/python3 scraper_bfa.py --account j; then
    echo "[$(date '+%H:%M:%S')] BFAJ: done"
else
    echo "[$(date '+%H:%M:%S')] BFAJ: FAILED (exit $?)"
    FAILED=$((FAILED + 1))
    FAILED_NAMES="${FAILED_NAMES}BFAJ, "
fi
echo ""

# Run BetOnline scraper — if token refresh fails, auto-recon and retry
echo "[$(date '+%H:%M:%S')] Running BetOnline scraper..."
echo "----------------------------------------"
if ./venv/bin/python3 scraper_betonline.py --since-last; then
    echo "[$(date '+%H:%M:%S')] BetOnline: done"
else
    echo "[$(date '+%H:%M:%S')] BetOnline: token may be expired, running auto-recon..."
    if ./venv/bin/python3 recon_betonline.py; then
        echo "[$(date '+%H:%M:%S')] Recon succeeded, retrying scraper..."
        if ./venv/bin/python3 scraper_betonline.py --since-last; then
            echo "[$(date '+%H:%M:%S')] BetOnline: done (after recon)"
        else
            echo "[$(date '+%H:%M:%S')] BetOnline: FAILED after recon (exit $?)"
            FAILED=$((FAILED + 1))
            FAILED_NAMES="${FAILED_NAMES}BetOnline, "
        fi
    else
        echo "[$(date '+%H:%M:%S')] BetOnline: FAILED — recon also failed (exit $?)"
        FAILED=$((FAILED + 1))
        FAILED_NAMES="${FAILED_NAMES}BetOnline, "
    fi
fi
echo ""

echo "========================================"
if [ $FAILED -eq 0 ]; then
    MSG="All 7 scrapers completed successfully."
    echo "$MSG"
    /usr/local/bin/terminal-notifier -title "Bet Logger ✓" -message "$MSG" -sound Glass -group betlogger
else
    # Trim trailing ", " from the list of failed scraper names
    FAILED_NAMES="${FAILED_NAMES%, }"
    MSG="Failed: ${FAILED_NAMES}"
    echo "$MSG"
    /usr/local/bin/terminal-notifier -title "Bet Logger — ${FAILED} failed" -message "$MSG" -sound Basso -group betlogger
fi
echo "Finished: $(date '+%Y-%m-%d %H:%M:%S')"
echo "========================================"

exit $FAILED
