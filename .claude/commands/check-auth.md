# Check Sportsbook Auth Tokens

Test connectivity and auth for all configured sportsbooks (or specific book: $ARGUMENTS).

## Steps

1. **Identify books to check**
   - If argument provided, check only that book
   - Otherwise check all: wagerzon, hoop88, bfa, bookmaker, bet105

2. **For each book**, run its recon script:
   - `wagerzon_odds/recon_wagerzon.py`
   - `hoop88_odds/recon_hoop88.py`
   - `bookmaker_odds/recon_bookmaker.py`
   - `bet105_odds/recon_bet105.py`
   - `bfa_odds/` — test with a single API call

3. **Report status**
   ```
   Book        | Status | Token Age | Notes
   ------------|--------|-----------|------
   wagerzon    | OK     | 2h        |
   hoop88      | EXPIRED| 26h       | Needs manual login
   bfa         | OK     | 1h        |
   bookmaker   | OK     | 4h        |
   bet105      | OK     | N/A       | API key (no expiry)
   ```

4. **For expired tokens**: Provide instructions on how to refresh
   - Playwright books: May need manual browser login to refresh session
   - REST API books: Run recon script to get new token
   - API key books: Check if key is still valid

## Notes
- Auth token expiration is the #2 source of pipeline failures
- BFA tokens expire frequently — always check BFA first
- Bookmaker uses curl_cffi with specific TLS fingerprinting — if recon fails, the site may have updated
