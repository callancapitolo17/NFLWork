#!/usr/bin/env bash
# Update and (re)start Unabated Ticket on the server, then check it answers.
#
#   ~/NFLWork/unabated_ticket/deploy/deploy.sh             # pull main, rebuild, restart, check
#   ~/NFLWork/unabated_ticket/deploy/deploy.sh --no-pull   # deploy the checkout as it is
#
# Inputs:  unabated_ticket/deploy/.env (from .env.example: BETS_EXTRA_ALLOWED_HOSTS,
#          UNABATED_DATA_DIR, KALSHI_PEM_PATH, DEPLOY_UID/GID) — paths and names only.
# Effects: `git pull --ff-only origin main` in the repo (refuses any other branch
#          unless --no-pull); creates UNABATED_DATA_DIR/logs (mode 700); rebuilds
#          the bets image and recreates both containers (compose.yaml), so the
#          bind-mounted code is re-read. bets.duckdb and the Novig token in the data
#          dir are kept. Prints PASS/FAIL per check and exits 1 on any FAIL.
# Never prints a credential: the checks print status codes, not bodies.
set -euo pipefail

BETS_URL="http://127.0.0.1:8094"
LOOPBACK_HOST="127.0.0.1:8094"
FOREIGN_HOST="evil.example"
WAIT_SEC=90
POLL_SEC=3

fail() {
  echo "deploy: FAIL — $*" >&2
  exit 1
}

# One KEY=value from deploy/.env, quotes stripped; empty when absent.
env_value() {
  local env_file="$1" key="$2"
  sed -n "s/^${key}=//p" "$env_file" | tail -n 1 | sed -e 's/^["'\'']//' -e 's/["'\'']$//'
}

# The status code of GET $BETS_URL$path sent with this Host header ("000" when
# nothing answers). --noproxy: the request must reach the VM's own loopback.
status_of() {
  local path="$1" host="$2"
  curl --noproxy '*' -sS -o /dev/null -w '%{http_code}' --max-time 5 \
    -H "Host: ${host}" "${BETS_URL}${path}" 2>/dev/null || true
}

# Polls until GET path with this Host returns the expected status or WAIT_SEC
# passes (the containers need a few seconds to bind). Prints PASS/FAIL.
check() {
  local label="$1" path="$2" host="$3" expected="$4"
  local waited=0 status
  while :; do
    status="$(status_of "$path" "$host")"
    if [ "$status" = "$expected" ]; then
      echo "PASS  ${label}: GET ${path} Host ${host} -> ${status}"
      return 0
    fi
    if [ "$waited" -ge "$WAIT_SEC" ]; then
      echo "FAIL  ${label}: GET ${path} Host ${host} -> ${status}, expected ${expected}" >&2
      return 1
    fi
    sleep "$POLL_SEC"
    waited=$((waited + POLL_SEC))
  done
}

main() {
  local pull=1
  case "${1:-}" in
    "") ;;
    --no-pull) pull=0 ;;
    *) fail "unknown argument ${1}; expected nothing or --no-pull" ;;
  esac

  local deploy_dir repo_root env_file
  deploy_dir="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
  repo_root="$(cd "$deploy_dir/../.." && pwd)"
  env_file="$deploy_dir/.env"

  command -v curl >/dev/null || fail "curl is not installed (sudo apt install -y curl)"
  docker info >/dev/null 2>&1 \
    || fail "cannot reach Docker as $(id -un); run: sudo usermod -aG docker $(id -un), then log out and back in"
  [ -f "$env_file" ] || fail "missing ${env_file}; copy deploy/.env.example to it and fill it in"

  if [ "$pull" -eq 1 ]; then
    local branch
    branch="$(git -C "$repo_root" rev-parse --abbrev-ref HEAD)"
    [ "$branch" = "main" ] || fail "the checkout is on ${branch}, not main; check out main or pass --no-pull"
    echo "==> git pull --ff-only origin main"
    git -C "$repo_root" pull --ff-only origin main
  fi
  echo "==> deploying $(git -C "$repo_root" log -1 --format='%h %s')"

  local data_dir pem_path tailnet_hosts tailnet_host
  data_dir="$(env_value "$env_file" UNABATED_DATA_DIR)"
  pem_path="$(env_value "$env_file" KALSHI_PEM_PATH)"
  tailnet_hosts="$(env_value "$env_file" BETS_EXTRA_ALLOWED_HOSTS)"
  [ -n "$data_dir" ] || fail "UNABATED_DATA_DIR is empty in ${env_file}"
  [ -n "$tailnet_hosts" ] || fail "BETS_EXTRA_ALLOWED_HOSTS is empty in ${env_file}"
  [ -f "$pem_path" ] || fail "KALSHI_PEM_PATH (${pem_path}) is not a file"
  tailnet_host="${tailnet_hosts%%,*}"
  tailnet_host="${tailnet_host// /}"

  # Docker would create a missing bind source as root; the service runs as DEPLOY_UID.
  mkdir -p "$data_dir/logs"
  chmod 700 "$data_dir"
  if [ -f "$repo_root/bet_logger/.env" ]; then
    echo "WARN  ${repo_root}/bet_logger/.env exists: the bets service reads it, so BFA/Wagerzon" \
         "would log in from this VM. Remove it unless you meant to move those books here." >&2
  fi

  local compose=(docker compose --project-directory "$deploy_dir" -f "$deploy_dir/compose.yaml")
  echo "==> docker compose build"
  "${compose[@]}" build --pull
  echo "==> docker compose up"
  "${compose[@]}" up -d --force-recreate --remove-orphans

  echo "==> health checks (up to ${WAIT_SEC}s each)"
  local failed=0
  check "bets service"             /health     "$LOOPBACK_HOST" 200 || failed=1
  check "tailnet name allowed"     /health     "$tailnet_host"  200 || failed=1
  check "phone page"               /           "$tailnet_host"  200 || failed=1
  check "edges via the runner"     /edges.json "$tailnet_host"  200 || failed=1
  check "foreign Host refused"     /health     "$FOREIGN_HOST"  403 || failed=1

  if [ "$failed" -ne 0 ]; then
    echo "deploy: FAIL — see: docker compose --project-directory ${deploy_dir} logs --tail 100" >&2
    exit 1
  fi
  echo "deploy: PASS — open https://${tailnet_host}/ on a device in the tailnet"
}

# Everything runs inside main(), so bash has parsed the whole file before the
# git pull can rewrite it.
main "$@"
