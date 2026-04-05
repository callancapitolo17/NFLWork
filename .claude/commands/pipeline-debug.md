# Pipeline Debug

Run the pipeline and diagnose any issues. Sport: $ARGUMENTS (default: cbb)

## Steps

1. **Pre-flight checks**
   - Verify `.env` files exist for all scraper directories
   - Check that DuckDB files are not locked (no stale WAL files)
   - Confirm we're on the correct branch (`git branch`)

2. **Run the pipeline**
   ```bash
   cd "/Users/callancapitolo/NFLWork/Answer Keys" && python run.py --sport <sport>
   ```

3. **Analyze output**
   - Check exit code (0 = success)
   - Review timing breakdown for slow scrapers
   - Look for "Warning: Unknown team name" messages (team mapping gaps)
   - Check for empty results from any scraper
   - Verify R script completed without errors

4. **If scraper failed**: Check auth tokens
   - Run the relevant recon script to test connectivity
   - Check if tokens in `.env` are expired
   - For Playwright scrapers: check if site structure changed

5. **If R script failed**: Check data
   - Look at DuckDB tables for the sport — are there fresh odds?
   - Check for team name mismatches between Odds API and scrapers
   - Look for NA/NULL values in critical columns

6. **Verify dashboard**
   - Check if `report.html` was regenerated (modification time)
   - Start/restart the Flask server if needed
   - Open dashboard and verify data looks correct

7. **Report**
   - Summarize: which scrapers succeeded/failed, how many games/markets found
   - Flag any team name mapping gaps
   - Flag any data quality issues (stale odds, missing markets)
