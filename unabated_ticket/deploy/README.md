# Unabated Ticket on the Oracle VM (phone page plan step 4)

Three containers on the VM `mlb-stack` (Ubuntu 22.04, ARM64), all bound to
**127.0.0.1 only**: the bets service (:8094, phone page, `/bets.json`,
`/settings.json`, `/edges.json` proxy), the Edges runner (:8095) and the NFL
fecta pricer `fectas` (:8096, `nfl_specials/`). The phone reaches :8094 and
:8096 through `tailscale serve` on the host, over the tailnet.
Nothing listens on a public interface; Oracle's security list keeps only SSH.

| File | What it is |
|---|---|
| `compose.yaml` | all three services, `network_mode: host`, `restart: unless-stopped`, json-file logs rotated at 10 MB × 3 |
| `bets_service.Dockerfile` | `python:3.12-slim` + `bets_service/requirements.txt`, wheels only (all have aarch64 wheels); no code, no secrets |
| `nfl_specials.Dockerfile` | the same for the fecta pricer: `python:3.12-slim` + `nfl_specials/requirements.txt` |
| `.env.example` | copy to `deploy/.env` (gitignored): tailnet name, data dir, Kalshi `.pem` path, uid/gid, the fecta pricer's DraftKings relay URL |
| `deploy.sh` | pull `main`, rebuild, recreate all three containers, health-check; prints PASS/FAIL |

Where things live on the VM:

| Host path | In the container | Mode |
|---|---|---|
| `~/NFLWork` (the repo) | `/app` | read-only (code; an update is pull + restart) |
| `~/NFLWork/unabated_ticket/bets_service/.env` | inside `/app` | read-only: `KALSHI_API_KEY_ID`, `POLYMARKET_US_KEY_ID`, `POLYMARKET_US_SECRET_KEY`, and the book logins `BFA_USERNAME`, `BFA_PASSWORD`, `WAGERZONC_USERNAME`, `WAGERZONC_PASSWORD` (step 10) |
| `KALSHI_PEM_PATH` (e.g. `~/unabated-secrets/kalshi.pem`) | `/run/secrets/kalshi.pem` | read-only |
| `UNABATED_DATA_DIR` (e.g. `~/unabated-data`) | `/data` | read-write: `bets.duckdb`, `novig_token.json` (rotates; rewritten in place), `recon_betonline_cookies.json` (BetOnline's refresh token; rotates, step 10), `logs/bets_service.log`, and the fecta pricer's `nfl_specials.duckdb` and `logs/.placement_debug_*.json` (`WAGERZON_DEBUG_DIR`: Wagerzon's answer to each Place that did not place) |

`bet_logger/.env` is never put on the VM: the BFA and Wagerzon logins go in
`bets_service/.env` instead, and BetOnline's token file in the data dir
(step 10). A book with neither reports "not configured". The fecta pricer
reads only the `WAGERZON*_USERNAME` / `_PASSWORD` lines of that same file.

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
   `id -u` / `id -g`, and `NFL_SPECIALS_DK_PROXY` (step 12; deploy.sh refuses
   to run without it). The bets service refuses to start on a malformed name.
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
     started again. The weekly sheet scraper (`bet_logger/scraper_betonline.py`)
     reads BetOnline from this service's `/bets.json` over Tailscale, so it
     keeps working; never run its old login on the Mac again.
   - **Check, then restart:** the login check below. `ok` for all three means
     done. If BetOnline fails on Cloudflare (its cookies were issued to the
     Mac's IP), move the file back to the Mac the same way.

11. **Automated taker: its own subaccount and key** (2026-10-16; the bot
    places real Kalshi orders once this is in). Public SSH stays off: SSH over
    Tailscale only (the key can trade). The setup must sign with a Kalshi key
    that has **write** scope (an unrestricted key: it makes a subaccount,
    moves money and mints a key). The service's own `KALSHI_API_KEY_ID` is
    read-only, and a read-only key stops at once with "the signing key ... is
    read-only", nothing changed at Kalshi. So override both settings for this
    one run with `-e`, pointing at a **copy** of a write key's file in
    `~/unabated-data/` (the setup overwrites `kalshi_taker.pem`, so never sign
    from that path itself). On the VM, from `~/NFLWork/unabated_ticket/deploy`:
    ```bash
    cp -p <write key file> ~/unabated-data/kalshi_setup_signing.pem && chmod 600 ~/unabated-data/kalshi_setup_signing.pem
    docker compose exec -T -e KALSHI_API_KEY_ID=<write key id> \
      -e KALSHI_PRIVATE_KEY_PATH=/data/kalshi_setup_signing.pem \
      bets python -m unabated_ticket.bets_service.taker_subaccount_setup
    rm ~/unabated-data/kalshi_setup_signing.pem   # afterwards, Done or STOPPED
    ```
    Use a main-account key whose scopes, in Kalshi's `GET /api_keys` listing,
    include `read` and `write` (a key made at kalshi.com listed as both on
    2026-10-16); check before using the Place key of step 12. (2026-10-16 the write key
    was the earlier bot key, then at `~/unabated-data/kalshi_taker.pem`; `-e`
    on `exec` wins over `.env`, which stays untouched.) It puts the stop file (`~/unabated-data/taker.stop`) in
    place first, so nothing trades until you `rm` it; that holds on a STOPPED
    run too, and on a VM where the bot already trades it pauses the bot. Then
    it makes subaccount 1 (asking for Kalshi's Advanced API tier if needed),
    moves $500 into it once and checks it landed (waiting up to 30 s for
    Kalshi's balance to show it), mints a key locked to it
    with no transfer permission, writes it to `/data/kalshi_taker.pem` (mode
    0600, replacing the old bot key file; never printed) and proves the lock:
    Kalshi's key listing must show it on subaccount 1 with only `read` and
    `write::trade`, and it must read that balance and be refused the primary's.
    A key that fails is deleted from disk and the output names it for
    revoking. Re-running it is safe (it reuses the subaccount, transfers
    nothing while the subaccount holds money, and keeps a key that still
    passes); the one case it stops on is a subaccount back at $0 after the
    $500 went in, which is a top-up you make yourself as a transfer. Then set
    the two lines it prints in `~/NFLWork/unabated_ticket/bets_service/.env`
    (`KALSHI_TAKER_API_KEY_ID=<new id>`, `KALSHI_TAKER_SUBACCOUNT=1`) and:
    ```bash
    docker compose restart bets   # the log says "taker: ARMED on its own key, locked to Kalshi subaccount 1"
    ```
    Revoke any earlier bot key at kalshi.com > API keys. Stop the bot:
    `touch ~/unabated-data/taker.stop` (`rm` to resume). For good: revoke the
    key. State: `curl -s 127.0.0.1:8095/taker.json`.
12. **The Ticket's Place key** (2026-10-16, Cal: its own key, "Unabated Place").
    At kalshi.com > Account & security > API keys, make an RSA key on the MAIN
    account (no subaccount) with read + trade, and save its .pem. From the Mac:
    `scp -i ~/.ssh/oracle_mlb.key ~/Downloads/<file>.pem
    ubuntu@100.66.14.98:unabated-data/kalshi_place.pem`, then on the VM `chmod
    600 ~/unabated-data/kalshi_place.pem`, add `KALSHI_PLACE_API_KEY_ID=<its id>`
    to `~/NFLWork/unabated_ticket/bets_service/.env` and `docker compose restart
    bets` (this restarts the automated taker too; it re-arms by itself). The log
    says "kalshi place: on, signing with Place key". Without it Place answers
    503 "no Kalshi trade key for Place": the read key is read-only.
13. **Automated taker on Novig: its own subaccount** (2026-10-16, Cal: $500 of
    its own; the bot places real Novig orders once this is in). Novig's
    management key (`~/.novig/novig_mgmt.pem` + `novig_mgmt.key_id`) stays on
    the Mac; the bot's trading key is made on the VM and its private half
    never moves. One line on the Mac, from `~/NFLWork`:
    ```bash
    ssh -i ~/.ssh/oracle_mlb.key ubuntu@100.66.14.98 \
      'cd ~/NFLWork/unabated_ticket/deploy && docker compose exec -T bets python -m unabated_ticket.bets_service.novig_subaccount_setup keygen --name bot' \
      | kalshi_draft/venv/bin/python3 -m unabated_ticket.bets_service.novig_subaccount_setup open --name bot --fund 500
    ```
    The VM half, making a NEW key, puts the Novig stop file
    (`~/unabated-data/novig_taker.stop`) in place first (a kept key, a
    top-up, leaves a trading bot running), writes the key to `/data/novig_bot.pem` (mode 0600; an
    existing usable key is kept) and prints only its public half. The Mac half
    opens the subaccount "unabated-bot" with it (or finds it), checks Novig's
    record of the key (scope trading, that subaccount, the VM key's
    fingerprint), moves $500 in from the app's cash and waits for Novig to
    apply it, then prints `NOVIG_TAKER_KEY_ID=<id>`. Re-running is safe: the
    transfer id is fixed per run (`--run-id`, "initial" by default), so a
    retry never moves money twice; a top-up is the same line with
    `--run-id topup-1` and the amount. Set that line in
    `~/NFLWork/unabated_ticket/bets_service/.env`, `docker compose up -d`
    (the log must say "novig taker: ARMED on its own trading key"; else it
    names why the desk is off), then `rm ~/unabated-data/novig_taker.stop` to
    start trading. `curl -s 127.0.0.1:8095/taker.json`: `pass.novigSkipped`
    says why a pass decided nothing on Novig, `novigDesk` what the desk last
    answered (`armed`, `stopped: ...`, `off: ...`, `geo_blocked: ...`; `not
    asked yet` until a Novig bet qualifies). Stop it: `touch` that
    file (or `taker.stop`, which stops both bots). A 451 means Novig wants the
    phone app opened (location check every 3 days): `novigDesk` reads
    `geo_blocked` and nothing is placed until then. The Place button's
    subaccount is the same two halves with `--name place` and its own amount
    (its line is `NOVIG_PLACE_KEY_ID`; top it up by re-running `open --name
    place --fund <dollars>` with a new `--run-id`). After `docker compose up -d` the log must say "novig place: on,
    signing with Place key"; it has no stop file, every order is a click.
14. **The Ticket's Polymarket US Place key** (2026-10-21, Cal: option A, its
    own key on the MAIN account). At polymarket.us/developer make a new API
    key and copy its Secret Key (shown once) into a file on the Mac, e.g.
    `~/Downloads/polymarket_place.key` (the secret alone, one line). Don't
    upload it anywhere. From the Mac:
    `scp -i ~/.ssh/oracle_mlb.key ~/Downloads/polymarket_place.key
    ubuntu@100.66.14.98:unabated-data/polymarket_place.key`, then on the VM
    `chmod 600 ~/unabated-data/polymarket_place.key`, add
    `POLYMARKET_PLACE_KEY_ID=<its key id>` to
    `~/NFLWork/unabated_ticket/bets_service/.env`, `docker compose up -d`, and
    delete the Mac copy. The log must say "polymarket place: on, signing with
    Place key"; else it names why (no key, an unreadable file, Polymarket US
    refusing the balance read). It has no stop file: every order is a click.

12. **NFL fecta pricer** (`fectas`; its own README, Running on the VM):
    - On the Mac, install the DraftKings relay agent (`nfl_specials/README.md`,
      Running on the VM) and read its address: `tailscale ip -4`.
    - In `deploy/.env` set `NFL_SPECIALS_DK_PROXY=http://<that address>:8097`.
      Every deploy needs it (compose refuses to start without it).
    - Publish the page on the tailnet at port 8443:
      ```bash
      sudo tailscale serve --bg --https=8443 http://127.0.0.1:8096
      tailscale serve status        # must also show https://<name>:8443 proxy http://127.0.0.1:8096
      ```
    - Once, before the first deploy and with the Mac's page stopped, move its
      history over, then retire the Mac's copy as step 10 does:
      `scp ~/NFLWork/nfl_specials/nfl_specials.duckdb ubuntu@<vm>:~/unabated-data/`, then
      `mv ~/NFLWork/nfl_specials/nfl_specials.duckdb ~/NFLWork/nfl_specials/nfl_specials.duckdb.moved-to-vm`.
    - `./deploy.sh`, then open `https://<name>:8443/` on the phone.

## Day to day

```bash
cd ~/NFLWork/unabated_ticket/deploy
./deploy.sh                          # update: pull main, rebuild, restart, PASS/FAIL checks
docker compose ps                    # state
docker compose logs -f --tail 100 bets      # or: runner, fectas
tail -f ~/unabated-data/logs/bets_service.log
docker compose restart bets fectas   # e.g. after editing bets_service/.env
docker compose down                  # stop all three (data kept); `up -d` or deploy.sh starts again
sudo tailscale serve reset           # take both pages off the tailnet
sudo tailscale serve --https=8443 off   # take only the fecta pricer off
```

What is freezing the runner, if anything: a 3-minute CPU profile of the
live runner (no restart; it opens Node's inspector on the VM's loopback,
samples, then closes it), printing each busy stretch over 300 ms with its
call path. Piped in, nothing written on the VM:

```bash
docker compose exec -T runner node - < ~/NFLWork/unabated_ticket/server/profile_stalls.js
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

Expect `127.0.0.1:8094` (python), `127.0.0.1:8095` (node), `127.0.0.1:8096` (python), `sshd` on `:22`,
and local system ones (`127.0.0.53:53`). `tailscaled` may hold the tailnet
address (`100.x.y.z`) only. Anything else on `0.0.0.0`, `[::]` or the
public IP is wrong — `docker compose down` and look.

## What `deploy.sh` checks

Against `http://127.0.0.1:8094`, polling up to 90 s each: `/health` with the
loopback `Host` → 200; `/health`, `/` and `/edges.json` with the tailnet
name as `Host` (what `tailscale serve` forwards) → 200 (`/edges.json`
proves the runner is up); `/health` with `Host: evil.example` → 403 (the
guard is still on). Then against `http://127.0.0.1:8096`: `/` with the
loopback `Host` and with `<tailnet name>:8443` → 200, `/api/board` (reads the
pricer's database) → 200, `/` with `evil.example` → 403. It prints status codes, never bodies. On FAIL read
`docker compose logs --tail 100`. Run it with `--no-pull` to deploy the
checkout as it is (e.g. trying a branch).
