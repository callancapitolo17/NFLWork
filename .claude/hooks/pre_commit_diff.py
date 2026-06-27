#!/usr/bin/env python3
"""PreToolUse hook — surface the staged diff before a git commit.

NFLWork CLAUDE.md mandates "always git diff --stat before commit to confirm
all intended files are included." This hook injects the staged stat as context
so the commit is never blind.

Reads the PreToolUse event JSON from stdin. ALWAYS exits 0 (default-allow) so it
can never block a commit — worst case it is a silent no-op.
"""
import json
import subprocess
import sys


def main() -> None:
    try:
        data = json.load(sys.stdin)
    except Exception:
        return

    command = ""
    try:
        command = (data.get("tool_input") or {}).get("command", "") or ""
    except Exception:
        return

    if "git commit" not in command:
        return

    try:
        stat = subprocess.run(
            ["git", "diff", "--cached", "--stat"],
            capture_output=True, text=True, timeout=10,
        ).stdout.strip()
    except Exception:
        stat = ""

    if not stat:
        stat = "(nothing staged — did you forget to `git add`?)"

    msg = "Staged changes about to be committed:\n" + stat
    print(json.dumps({
        "hookSpecificOutput": {
            "hookEventName": "PreToolUse",
            "additionalContext": msg,
        }
    }))


if __name__ == "__main__":
    try:
        main()
    except Exception:
        pass
    sys.exit(0)
