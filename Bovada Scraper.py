from selenium import webdriver
from selenium.webdriver.common.by import By
from selenium.webdriver.chrome.options import Options
import pandas as pd
import time

# Setup Chrome
options = Options()
options.add_argument("--start-maximized")
driver = webdriver.Chrome(options=options)

# Step 1: Open Bovada and login manually
driver.get("https://www.bovada.lv/")
input("Log in to your account manually, then press Enter here to continue...")

# Step 2: Navigate to bet history
driver.get("https://www.bovada.lv/account/bets")

time.sleep(5)  # Adjust this based on your internet speed

# Step 3: Extract data
bets = driver.find_elements(By.CLASS_NAME, "bet-history-item")
data = []

for bet in bets:
  try:
  date = bet.find_element(By.CLASS_NAME, "bet-history-date").text
details = bet.find_element(By.CLASS_NAME, "bet-history-details").text
outcome = bet.find_element(By.CLASS_NAME, "bet-history-status").text
payout = bet.find_element(By.CLASS_NAME, "bet-history-payout").text

data.append({
  "date": date,
  "details": details,
  "outcome": outcome,
  "payout": payout
})
except Exception as e:
  print("Error processing a bet:", e)

# Step 4: Save to CSV
df = pd.DataFrame(data)
df.to_csv("bovada_bet_history.csv", index=False)

print("Saved bet history to bovada_bet_history.csv")
driver.quit()
