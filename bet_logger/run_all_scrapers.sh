#!/bin/bash

# Run All Bet Scrapers
# This script runs all betting platform scrapers in sequence

cd "/Users/callancapitolo/NFLWork/bet_logger"

echo "========================================"
echo "MULTI-PLATFORM BET SCRAPER"
echo "========================================"
echo ""

# Run Wagerzon scraper
echo "🎯 Running Wagerzon scraper..."
echo "----------------------------------------"
./venv/bin/python3 scraper.py
echo ""

# Run Hoop88 scraper
echo "🎯 Running Hoop88 scraper..."
echo "----------------------------------------"
./venv/bin/python3 scraper_hoop88.py
echo ""

echo "========================================"
echo "✅ All scrapers completed!"
echo "========================================"
