# Wagerzon Bet Logger

Automatically scrapes bet history from Wagerzon and uploads to Google Sheets.

## Features

- Auto-login to Wagerzon
- Scrapes bet history (defaults to "Last Week")
- Automatic duplicate detection (date + description + bet amount)
- Dynamic row detection (appends to next empty row)
- Calculates American odds and decimal odds

## Setup Instructions

### 1. Install Python Dependencies

```bash
cd "/Users/callancapitolo/Sports Stuff/bet_logger"
python3 -m venv venv
source venv/bin/activate
pip install -r requirements.txt
playwright install chromium
```

### 2. Set Up Google Cloud Credentials

1. Go to [Google Cloud Console](https://console.cloud.google.com/)
2. Create a new project (or select existing)
3. Enable the **Google Sheets API**:
   - Go to "APIs & Services" > "Library"
   - Search for "Google Sheets API"
   - Click "Enable"
4. Create a Service Account:
   - Go to "APIs & Services" > "Credentials"
   - Click "Create Credentials" > "Service Account"
   - Name it something like "bet-logger"
   - Click "Done"
5. Create a key for the service account:
   - Click on the service account you just created
   - Go to "Keys" tab
   - Click "Add Key" > "Create new key" > "JSON"
   - Save the downloaded file as `credentials.json` in this folder

### 3. Share Your Google Sheet

1. Open your betting spreadsheet
2. Click "Share" button
3. Add the service account email (looks like `bet-logger@your-project.iam.gserviceaccount.com`)
4. Give it "Editor" access

### 4. Configure Environment Variables

```bash
cp .env.example .env
```

Edit `.env` with your actual values:
- `WAGERZON_USERNAME`: Your Wagerzon username
- `WAGERZON_PASSWORD`: Your Wagerzon password
- `WAGERZON_URL`: The Wagerzon website URL

## Usage

### Run Full Scraper

```bash
cd "/Users/callancapitolo/Sports Stuff/bet_logger" && source venv/bin/activate && python scraper.py
```

This will:
1. Open a browser window
2. Navigate to Wagerzon and log in
3. Go to bet history
4. Scrape all bets from "Last Week"
5. Skip any duplicates
6. Upload new bets to your Google Sheet

### Test the Parser (without logging in)

```bash
cd "/Users/callancapitolo/Sports Stuff/bet_logger" && source venv/bin/activate && python scraper.py --test
```

### Test Google Sheets Connection

```bash
cd "/Users/callancapitolo/Sports Stuff/bet_logger" && source venv/bin/activate && python sheets.py --test
```

## Output Columns

The scraper fills these columns in your sheet:
- A: Date (e.g., "1/6/26")
- B: Platform ("Wagerzon")
- C: Sports (NFL, NHL, etc.)
- D: Bet Description
- E: Bet Type (Parlay, Straight, Prop, etc.)
- F: Line
- G: Odds (American format, e.g., "+450", "-110")
- H: Bet Amount
- I: Decimal Odds
- J: Result (win/loss)

## Troubleshooting

### Auto-login doesn't work
The script will wait 60 seconds for you to log in manually if auto-login fails.

### Permission denied on Google Sheets
Make sure you shared the spreadsheet with your service account email address.

### Duplicates still appearing
Duplicate detection matches on date + description + bet amount. If any of these differ slightly, it won't be detected as a duplicate.
