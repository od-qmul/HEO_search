import requests
import datetime
import csv
import os

# --- CONFIG ---
OWNER = "od-qmul"         # e.g. "torvalds"
REPO = "HEO_search"              # e.g. "linux"
OUTPUT_FILE = "daily_downloads.csv"

# --- AUTH (use GitHub token if available) ---
headers = {}
token = os.getenv("GITHUB_TOKEN")
if token:
    headers["Authorization"] = f"token {token}"

# --- SCRIPT ---
url = f"https://api.github.com/repos/{OWNER}/{REPO}/releases"
response = requests.get(url, headers=headers)

if response.status_code != 200:
    raise Exception(f"GitHub API error: {response.status_code} {response.text}")

releases = response.json()
today = datetime.date.today().isoformat()

# Load previous totals (to compute daily increments)
previous_totals = {}
if os.path.exists(OUTPUT_FILE):
    with open(OUTPUT_FILE, "r") as f:
        reader = csv.DictReader(f)
        for row in reader:
            key = (row["release_tag"], row["asset_name"])
            previous_totals[key] = int(row["cumulative_downloads"])

# Prepare new rows for today
rows = []
for release in releases:
    tag = release["tag_name"]
    for asset in release.get("assets", []):
        asset_name = asset["name"]
        total = asset["download_count"]

        key = (tag, asset_name)
        prev_total = previous_totals.get(key, 0)
        daily = total - prev_total

        rows.append({
            "date": today,
            "release_tag": tag,
            "asset_name": asset_name,
            "daily_downloads": daily,
            "cumulative_downloads": total
        })

# Write/append to CSV
write_header = not os.path.exists(OUTPUT_FILE)
with open(OUTPUT_FILE, "a", newline="") as f:
    writer = csv.DictWriter(f, fieldnames=["date", "release_tag", "asset_name", "daily_downloads", "cumulative_downloads"])
    if write_header:
        writer.writeheader()
    writer.writerows(rows)

print(f"Saved daily downloads for {today}")
