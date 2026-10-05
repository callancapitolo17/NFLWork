# Unabated Ticket on the Oracle VM (phone page plan step 4)

Two containers on the VM `mlb-stack` (Ubuntu 22.04, ARM64), both bound to
**127.0.0.1 only**: the bets service (:8094 — phone page, `/bets.json`,
`/settings.json`, `/edges.json` proxy) and the Edges runner (:8095). The
phone reaches :8094 through `tailscale serve` on the host, over the tailnet.
Nothing listens on a public interface; Oracle's security list keeps only SSH.

| File | What it is |
|---|---|
| `compose.yaml` | both services, `network_mode: host`, `restart: unless-stopped`, json-file logs rotated at 10 MB × 3 |
| `bets_service.Dockerfile` | `python:3.12-slim` + `bets_service/requirements.txt`, wheels only (all have aarch64 wheels); no code, no secrets |
| `.env.example` | copy to `deploy/.env` (gitignored): tailnet name, data dir, Kalshi `.pem` path, uid/gid |
| `deploy.sh` | pull `main`, rebuild, recreate both containers, health-check; prints PASS/FAIL |

Where things live on the VM:

| Host path | In the container | Mode |
|---|---|---|
| `~/NFLWork` (the repo) | `/app` | read-only (code; an update is pull + restart) |
| `~/NFLWork/unabated_ticket/bets_service/.env` | inside `/app` | read-only: `KALSHI_API_KEY_ID`, `POLYMARKET_US_KEY_ID`, `POLYMARKET_US_SECRET_KEY`, and the book logins `BFA_USERNAME`, `BFA_PASSWORD`, `WAGERZONC_USERNAME`, `WAGERZONC_PASSWORD` (step 10) |
| `KALSHI_PEM_PATH` (e.g. `~/unabated-secrets/kalshi.pem`) | `/run/secrets/kalshi.pem` | read-only |
| `UNABATED_DATA_DIR` (e.g. `~/unabated-data`) | `/data` | read-write: `bets.duckdb`, `novig_token.json` (rotates; rewritten in place), `recon_betonline_cookies.json` (BetOnline's refresh token; rotates, step 10), `logs/bets_service.log` |

`bet_logger/.env` is never put on the VM: the BFA and Wagerzon logins go in
`bets_service/.env` instead, and BetOnline's token file in the data dir
(step 10). A book with neither reports "not configured".

## One-time setup

Run as `ubuntu` on the VM (`ssh -i ~/.ssh/oracle_mlb.key ubuntu@<public ip>`).

1. **Old MLB containers stay stopped.** `docker ps -a`; for each old MLB
   container, `docker update --restart=no <name>` so a reboot cannot revive
   it. Do not start them.
2. **Docker without sudo.** `docker info` must work as `ubuntu`; if not:
   `sudo usermod -aG docker ubuntu`, log out and back in. `sudo systemctl
   enable docker` so the containers come back after a reboot.
3. **Repo.** If `~/NFLWork` is not there yet, clone it with the read-only
   deploy key (main README § Running on a server, step 0 item 2). Then let
   `git pull` use that key: `git -C ~/NFLWork config core.sshCommand "ssh -i ~/.ssh/nflwork"`.
4. **Credentials out of `~/betscheck`.** Step 0's copy holds the keys; move
   them, then delete it:
   ```bash
   ls -la ~/betscheck                       # see what is there
   mkdir -p ~/unabated-secrets ~/unabated-data && chmod 700 ~/unabated-secrets ~/unabated-data
   mv ~/betscheck/<the kalshi .pem> ~/unabated-secrets/kalshi.pem && chmod 600 ~/unabated-secrets/kalshi.pem
   # KALSHI_API_KEY_ID and the two POLYMARKET_US_* lines go here (an editor, not echo — keeps them out of shell history):
   nano ~/NFLWork/unabated_ticket/bets_service/.env && chmod 600 ~/NFLWork/unabated_ticket/bets_service/.env
   # Only if step 0 made a Novig token ON THE VM: move it, never copy (Auth0 revokes a reused refresh token).
   mv ~/betscheck/<path>/novig_token.json ~/unabated-data/novig_token.json
   rm -rf ~/betscheck                       # after deploy.sh says PASS
   ```
   A `KALSHI_PRIVATE_KEY_PATH` line in that `.env` is harmless: compose
   overrides it with `/run/secrets/kalshi.pem`.
5. **Tailscale** (host, not a container):
   ```bash
   curl -fsSL https://tailscale.com/install.sh | sh
   sudo tailscale up                        # prints a login URL: open it, log in with your account
   tailscale status --json | python3 -c 'import json,sys; print(json.load(sys.stdin)["Self"]["DNSName"].rstrip("."))'
   ```
   The last line prints the VM's MagicDNS name, e.g. `mlb-stack.tail1234.ts.net`.
   In the Tailscale admin console (DNS page) MagicDNS and **HTTPS
   Certificates** must be on. Then:
   ```bash
   sudo tailscale serve --bg 8094           # https://<name>/ -> http://127.0.0.1:8094, survives reboots
   tailscale serve status                   # must show https://<name> proxy http://127.0.0.1:8094
   ```
   If `serve` prints a link instead (HTTPS not enabled for the tailnet), open
   it, enable, and rerun. The first HTTPS request fetches the certificate
   (a few seconds). The certificate puts the name in public CT logs; the
   name is still reachable only from inside the tailnet.
6. **deploy/.env.** `cp ~/NFLWork/unabated_ticket/deploy/.env.example ~/NFLWork/unabated_ticket/deploy/.env`,
   set `BETS_EXTRA_ALLOWED_HOSTS` to the name from step 5 (lowercase, no
   `https://`, no port), the two paths, and `DEPLOY_UID`/`DEPLOY_GID` to
   `id -u` / `id -g`. The bets service refuses to start on a malformed name.
   Why it is needed: `tailscale serve` forwards the browser's `Host` header
   (the tailnet name), and the #125 DNS-rebinding guard 403s any name it does
   not list.
7. **First deploy:** `~/NFLWork/unabated_ticket/deploy/deploy.sh`.
8. **Novig** (if no token was moved in step 4): its own login on the VM,
   with the service stopped so no second process rotates the token:
   ```bash
   cd ~/NFLWork/unabated_ticket/deploy
   docker compose stop bets
   docker compose run --rm bets python -m unabated_ticket.bets_service.sources.novig_auth connect
   docker compose start bets
   ```
9. **Phone:** install the Tailscale app, log in with the same account, open
   `https://<name>/` (Edges) or `https://<name>/tracker` (Bet Tracker).
10. **Book logins: BFA, Wagerzon, BetOnline** (2026-10-05, Cal's call: these
   books don't mind a data-center login). With all three on the VM the
   tracker covers every polled book with the Mac off. **Only when the VM
   becomes the primary host:** while the Mac is (the current setup, main
   README § Bet Tracker), skip this step, because moving BetOnline's token
   takes it away from the Mac.
   - **BFA and Wagerzon** are plain password logins, so the Mac can keep its
     own. Add four lines to `~/NFLWork/unabated_ticket/bets_service/.env`
     with an editor: `BFA_USERNAME`, `BFA_PASSWORD`, `WAGERZONC_USERNAME`,
     `WAGERZONC_PASSWORD` (the values from the Mac's `bet_logger/.env`; the
     Wagerzon C account, or `WAGERZON_*` if that is where its login sits).
   - **BetOnline: move the token, never copy it.** Every refresh rotates its
     Keycloak token and a reused one kills the chain, so only one machine
     can hold it. On the Mac first stop everything that refreshes it: the
     bets service (`launchctl bootout gui/$(id -u)/com.nflwork.bets-service`)
     and the weekly sheet scraper (`launchctl unload ~/Library/LaunchAgents/com.callancapitolo.betlogger.plist`;
     load it again afterwards if you still want the sheet). Then:
     ```bash
     scp ~/NFLWork/bet_logger/recon_betonline_cookies.json ubuntu@<vm>:~/unabated-data/
     mv ~/NFLWork/bet_logger/recon_betonline_cookies.json ~/NFLWork/bet_logger/recon_betonline_cookies.json.moved-to-vm
     ```
     On the VM `chmod 600 ~/unabated-data/recon_betonline_cookies.json`. The
     Mac's bets service then reports BetOnline "not configured" and can be
     started again; the weekly sheet scraper's BetOnline step fails until the
     file comes back.
   - **Check, then restart:** the login check below. `ok` for all three means
     done. If BetOnline fails on Cloudflare (its cookies were issued to the
     Mac's IP), move the file back to the Mac the same way.

## Day to day

```bash
cd ~/NFLWork/unabated_ticket/deploy
./deploy.sh                          # update: pull main, rebuild, restart, PASS/FAIL checks
docker compose ps                    # state
docker compose logs -f --tail 100 bets      # or: runner
tail -f ~/unabated-data/logs/bets_service.log
docker compose restart bets          # e.g. after editing bets_service/.env
docker compose down                  # stop both (data kept); `up -d` or deploy.sh starts again
sudo tailscale serve reset           # take the page off the tailnet
```

Login check on the VM (step 10 too), as step 0 but in the container — stop the service
first, because a Novig check rotates the same refresh token the running
service holds:

```bash
docker compose stop bets && docker compose run --rm bets python -m unabated_ticket.bets_service.check_sources; docker compose start bets
```

## Nothing public listens

```bash
sudo ss -tlnp
```

Expect `127.0.0.1:8094` (python), `127.0.0.1:8095` (node), `sshd` on `:22`,
and local system ones (`127.0.0.53:53`). `tailscaled` may hold the tailnet
address (`100.x.y.z`) only. Anything else on `0.0.0.0`, `[::]` or the
public IP is wrong — `docker compose down` and look.

## What `deploy.sh` checks

Against `http://127.0.0.1:8094`, polling up to 90 s each: `/health` with the
loopback `Host` → 200; `/health`, `/` and `/edges.json` with the tailnet
name as `Host` (what `tailscale serve` forwards) → 200 (`/edges.json`
proves the runner is up); `/health` with `Host: evil.example` → 403 (the
guard is still on). It prints status codes, never bodies. On FAIL read
`docker compose logs --tail 100`. Run it with `--no-pull` to deploy the
checkout as it is (e.g. trying a branch).
