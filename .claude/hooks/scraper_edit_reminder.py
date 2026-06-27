#!/usr/bin/env python3
"""PostToolUse hook — scraper edit reminder.

When Claude edits a scraper file, inject a reminder to run the timezone
parity regression gate and keep game_start_time as TIMESTAMPTZ UTC
(NFLWork CLAUDE.md housekeeping rule #6).

Reads the PostToolUse event JSON from stdin. ALWAYS exits 0 so it can never
block an edit — worst case it is a silent no-op.
"""
import json
import sys


def main() -> None:
    try:
        data = json.load(sys.stdin)
    except Exception:
        return  # malformed input -> no-op

    file_path = ""
    try:
        file_path = (data.get("tool_input") or {}).get("file_path", "") or ""
    except Exception:
        return

    low = file_path.lower()
    is_scraper = (
        "scraper" in low
        or (low.endswith(".py") and "_odds/" in low)
        or (low.endswith(".py") and "mlb_sgp/" in low)
    )
    if not is_scraper:
        return

    msg = (
        f"⚠️ Scraper file edited: {file_path}\n"
        "Before committing, run the timezone regression gate:\n"
        "    python3 tests/timezone_parity_test.py\n"
        "and confirm game_start_time is written as TIMESTAMPTZ in UTC "
        "(CLAUDE.md housekeeping #6)."
    )
    print(json.dumps({
        "hookSpecificOutput": {
            "hookEventName": "PostToolUse",
            "additionalContext": msg,
        }
    }))


if __name__ == "__main__":
    try:
        main()
    except Exception:
        pass  # never fail the tool call
    sys.exit(0)
