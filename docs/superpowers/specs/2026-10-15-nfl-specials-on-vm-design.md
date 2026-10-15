# NFL fecta pricer on the Oracle VM

**Problem.** The pricer (`nfl_specials`, :8096) only runs on the laptop.

**Approach.**
1. Run it as a third container in the VM's existing stack (`unabated_ticket/deploy/compose.yaml`): repo mounted read-only, `nfl_specials.duckdb` in the data dir (the Mac's copy moved over once), a Wagerzon-only secrets file. Phone URL `https://mlb-stack.tail87ab2d.ts.net:8443` via `tailscale serve`; `deploy.sh` checks it answers.
2. Its Host check accepts the VM's tailnet name (the same `deploy/.env` value the bets service uses).
3. DraftKings 403s the VM (Akamai blocks Oracle's IPs, rechecked 2026-10-15), so DK requests go through a small proxy on the Mac over Tailscale that only reaches DraftKings. Mac off: superfectas show unpriced, trifectas still price.

**Out of scope.** Superfectas while the Mac is off; MLB tools.

**Done when.** From the phone over Tailscale the board loads, FanDuel, BetMGM and DraftKings (via the Mac) all price, and a Wagerzon preview (no bet) succeeds from the VM.
