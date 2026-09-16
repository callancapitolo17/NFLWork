#!/usr/bin/env bash
# The one pre-merge check for unabated_ticket/: ESLint, the node test suite,
# and the bets service's pytest suite. Exits non-zero if any of them fails,
# and still runs all three so one run shows every failure.
#
#   ./unabated_ticket/check.sh          # from anywhere: repo root, a worktree, this dir
#
# Python resolves the way bets_service/run.sh does: the kalshi_draft venv in
# the MAIN checkout (it has duckdb + cryptography; a Claude worktree has no
# venv of its own), else python3. `npm install` once per checkout puts eslint
# in node_modules/ (gitignored).
set -uo pipefail
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
REPO_ROOT="$(cd "$SCRIPT_DIR/.." && pwd)"
MAIN_ROOT="${REPO_ROOT%%/.claude/worktrees/*}"  # Claude desktop worktrees
MAIN_ROOT="${MAIN_ROOT%%/.worktrees/*}"
PYTHON="$MAIN_ROOT/kalshi_draft/venv/bin/python3"
if [ ! -x "$PYTHON" ]; then PYTHON=python3; fi

if [ ! -x "$SCRIPT_DIR/node_modules/.bin/eslint" ]; then
  echo "check.sh: eslint is not installed — run: (cd \"$SCRIPT_DIR\" && npm install)" >&2
  exit 2
fi

failed=0
run_step() {
  local name="$1"; shift
  echo "==> $name"
  if "$@"; then
    echo "==> $name: ok"
  else
    echo "==> $name: FAILED" >&2
    failed=1
  fi
  echo
}

cd "$SCRIPT_DIR"
run_step "eslint"      npm run --silent lint
run_step "node tests"  npm test --silent
cd "$REPO_ROOT"
run_step "pytest"      "$PYTHON" -m pytest -q unabated_ticket/bets_service/tests

if [ "$failed" -ne 0 ]; then
  echo "check.sh: FAILED" >&2
  exit 1
fi
echo "check.sh: all checks passed"
