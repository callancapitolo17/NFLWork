#!/usr/bin/env python3
"""
BetMGM (Entain "CDS" platform) MLB SGP Recon Tool
==================================================

Goal of this recon
------------------
Determine whether we can programmatically, against BetMGM's US sportsbook:
  (a) list today's MLB events,
  (b) discover selection/option IDs for a run-line (spread) leg and a
      total (over/under) leg,
  (c) POST a priced SGP slip (spread + total combo) and get correlated odds.

What was confirmed during recon (2026-06-10, see module docstring + REPORT)
---------------------------------------------------------------------------
BetMGM US runs the Entain "CDS" read API. The host that actually serves
fixtures/markets/options/prices for US (NJ) is:

    https://sports.<state>.betmgm.com/cds-api/...        (e.g. state=nj)

NOT the EU-style host `cds-api.<state>.betmgm.com` (that does not resolve
for US). Confirmed live:

    GET https://sports.nj.betmgm.com/cds-api/bettingoffer/fixtures?x-bwin-accessid=&...
    -> HTTP 400 {"message":"Access id missing", ...}

So the endpoint exists and is reachable WITHOUT any cookie/login. It is
gated only by a query param `x-bwin-accessid` (a base64 client/region
identifier — NOT a user-auth token). There is ALSO a documented read API
mirror with swagger:

    https://sportsapi.<state>.betmgm.com/restapi/swagger.html
    base path: /offer/api/{sportId}/{country}/fixtures  (GET only)

The HARD parts:
  1. Obtaining a valid `x-bwin-accessid`. The SPA mints it at runtime
     AFTER it resolves the user's state from IP geolocation. On a network
     whose IP does NOT geolocate to a BetMGM-legal sports state, the SPA
     redirects sports.<state>.betmgm.com -> www.betmgm.com and never
     bootstraps, so no accessid is ever generated. (This is exactly what
     happened on the recon machine — every load parked on www.betmgm.com.)
     => You need to run the Playwright harvest step from an IP that
        geolocates to a legal state (NJ/PA/MI/etc.) — typically a VPN or
        a box physically in-state.
  2. The priced-slip / bet-builder (SGP) endpoint. BetMGM's correlated
     SGP pricing comes from the Angstrom engine and is a POST that the
     public CDS read swagger does NOT document. This script captures it
     live (Phase B) so we can see the exact path + body shape.

How to use
----------
This script has THREE modes:

  python recon_betmgm_sgp.py probe
      Pure-HTTP probe (no browser). Confirms the CDS host is reachable
      and reports the "Access id missing" gate. Runs anywhere. Fast.

  python recon_betmgm_sgp.py harvest [--state nj] [--headful]
      Launches Playwright, loads the live sportsbook, and intercepts the
      `x-bwin-accessid` the SPA uses + every cds-api / bettingoffer /
      betslip / price call (this is where the SGP POST endpoint reveals
      itself). MUST be run from a legal-state IP or it will just park on
      www.betmgm.com (the script detects and reports that).

  python recon_betmgm_sgp.py pull --accessid <AID> [--state nj]
      Given an accessid (from the harvest step), exercises the read API:
      MLB fixtures -> one fixture's spread + total markets -> option IDs,
      then ATTEMPTS the priced-slip combo endpoints with those options.

  python recon_betmgm_sgp.py mgm [--state pa]
      *** THE WORKING END-TO-END PATH (no browser, runs from California). ***
      1) Harvests a brand-correct accessid headlessly from the clientconfig
         endpoint (the `x-bwin-sports-api: prod` header is what unlocks the
         sports-app config block carrying `msApp.publicAccessId`).
      2) Pulls real MLB fixtures + spread/total option IDs.
      3) POSTs the bet-builder combo to /cds-api/bettingoffer/picks and reads
         the correlated Angstrom price out of `betBuilderPricingGroups`.

SOLVED (2026-06-15) — feasibility GREEN from a California IP
-----------------------------------------------------------
The prior recon failed because it tried to mint the accessid via the SPA,
which IP-geo-redirects sports.<state>.betmgm.com -> www.betmgm.com and never
bootstraps. The accessid does NOT actually require an in-state browser.

KEY UNLOCK: the clientconfig endpoint serves the full sports-app config
(with the accessid) to ANY IP, but ONLY when the request carries the header
`x-bwin-sports-api: prod`. Without it you get the 11 KB host-app shell
(no accessid); with it you get the ~97 KB sports config containing
`msApp.publicAccessId` (and `msConnection.publicAccessId`).

  GET https://www.<state>.betmgm.com/en/api/clientconfig
        ?browserUrl=<urlenc>&x-from-product=host-app
      headers: x-bwin-sports-api: prod, x-bwin-browser-url=<urlenc>
  -> JSON; grep "publicAccessId" -> base64(UUID), e.g. PA:
     YWIzOGYzMjgtNzU3OS00NjU1LTk1MjUtZjQ4Y2UxODQyOTY0

The accessid is brand-partitioned: a bwin.com accessid is rejected by the
betmgm host with {"message":"Access id not allowed for application"}.
It is STATIC/long-lived (identical across re-harvests) and per-state. NJ
omits the pubid in config; PA/MI/CO/VA/TN all return clean pubids.

SGP PRICE ENDPOINT (Angstrom bet-builder), confirmed live from CA:
  POST https://www.<state>.betmgm.com/cds-api/bettingoffer/picks
        ?x-bwin-accessid=<AID>&lang=en&country=US&userCountry=US
  body: {"tv1Picks":[{"fixtureId":"<fxId>","gameId":<marketId>,
          "resultId":<optionId>,"useLiveFallback":false,
          "pickGroupId":"<shared-uuid>"}, ...],"tv2Picks":[]}
  (legs sharing one pickGroupId = one same-game combo)
  -> resp.betBuilderPricingGroups[<gid>].odds.{odds,americanOdds}
     is the TRUE correlated price (providerId "Angstrom").

Reference: pretrehr/Sports-betting bookmakers/bwin.py (EU sibling, same
CDS shape — fixtures endpoint + accessid param). Adapted for US host.
"""

import argparse
import json
import re
import sys

from curl_cffi import requests

# ---------------------------------------------------------------------------
# Config
# ---------------------------------------------------------------------------
# Entain CDS sport id for baseball is 23; US country code is the CDS "US".
MLB_SPORT_ID = "23"
COUNTRY = "US"

# Keywords that flag a request as SGP / bet-builder / pricing related.
SGP_KEYWORDS = [
    "betslip", "bet-slip", "bettingslip", "samegame", "same-game", "sgp",
    "bet-builder", "betbuilder", "combo", "combination", "price", "pricing",
    "quickbet", "selections", "parlay", "accumulator", "angstrom", "calculate",
]


def cds_base(state: str) -> str:
    return f"https://sports.{state}.betmgm.com/cds-api"


def fixtures_url(state: str, accessid: str, sport_id: str = MLB_SPORT_ID) -> str:
    """The CDS fixtures read endpoint (mirrors bwin.py, US host)."""
    return (
        f"{cds_base(state)}/bettingoffer/fixtures"
        f"?x-bwin-accessid={accessid}"
        f"&lang=en&country={COUNTRY}&userCountry={COUNTRY}"
        f"&fixtureTypes=Standard&state=Latest"
        f"&offerMapping=Filtered&offerCategories=Gridable"
        f"&fixtureCategories=Gridable"
        f"&sportIds={sport_id}"
        f"&skip=0&take=50&sortBy=Tags"
    )


def _session() -> requests.Session:
    s = requests.Session(impersonate="chrome")
    s.headers.update({
        "Accept": "application/json, text/plain, */*",
        "Accept-Language": "en-US,en;q=0.9",
    })
    return s


# ---------------------------------------------------------------------------
# Mode 1: pure-HTTP probe (runs anywhere; no accessid needed)
# ---------------------------------------------------------------------------
def mode_probe(state: str):
    print(f"=== PROBE: CDS reachability for state={state} ===\n")
    s = _session()
    s.headers["Referer"] = f"https://sports.{state}.betmgm.com/"

    # 1) fixtures with EMPTY accessid -> should return the structured gate error
    url = fixtures_url(state, accessid="")
    print(f"[1] GET fixtures (no accessid)\n    {url[:110]}...")
    try:
        r = s.get(url, timeout=20)
        print(f"    -> HTTP {r.status_code}  ({len(r.text)} bytes)")
        print(f"    -> body: {r.text[:300]}")
        if r.status_code == 400 and "ccess id" in r.text.lower():
            print("    => CONFIRMED: endpoint live, gated only by x-bwin-accessid\n")
        elif r.status_code == 200:
            print("    => fixtures returned WITHOUT accessid (unexpected — inspect body)\n")
    except Exception as e:
        print(f"    -> ERROR: {e!r}\n")

    # 2) swagger read-API mirror reachability
    sw = f"https://sportsapi.{state}.betmgm.com/restapi/swagger.html"
    print(f"[2] GET swagger mirror\n    {sw}")
    try:
        r = s.get(sw, timeout=20)
        print(f"    -> HTTP {r.status_code} ({len(r.text)} bytes) "
              f"{'(docs reachable)' if r.status_code == 200 else ''}\n")
    except Exception as e:
        print(f"    -> ERROR: {e!r}\n")

    print("NEXT: run `harvest` from a legal-state IP to grab a live accessid.")


# ---------------------------------------------------------------------------
# Mode 2: Playwright harvest of accessid + SGP/price endpoints
# ---------------------------------------------------------------------------
def mode_harvest(state: str, headful: bool):
    try:
        from playwright.sync_api import sync_playwright
    except Exception as e:
        print(f"playwright not importable: {e!r}")
        sys.exit(1)

    accessids = set()
    cds_calls = []           # (method, url, post_data) for cds-api/bettingoffer
    flagged = []             # SGP/price-keyword requests

    def on_request(req):
        u = req.url
        # capture accessid wherever it appears
        if "x-bwin-accessid=" in u:
            m = re.search(r"x-bwin-accessid=([^&]+)", u)
            if m and m.group(1):
                accessids.add(m.group(1))
        h = req.headers
        if "x-bwin-accessid" in h and h["x-bwin-accessid"]:
            accessids.add(h["x-bwin-accessid"])

        low = u.lower()
        if "cds-api" in low or "bettingoffer" in low:
            cds_calls.append((req.method, u[:240], (req.post_data or "")[:300]))
        if any(k in low for k in SGP_KEYWORDS):
            flagged.append((req.method, u[:240], (req.post_data or "")[:400]))

    print(f"=== HARVEST: launching browser for state={state} (headful={headful}) ===\n")
    with sync_playwright() as p:
        b = p.chromium.launch(
            headless=not headful,
            args=["--disable-blink-features=AutomationControlled"],
        )
        # NJ centroid geolocation; pair with an in-state IP/VPN for it to matter.
        ctx = b.new_context(
            locale="en-US",
            geolocation={"latitude": 40.0583, "longitude": -74.4057},
            permissions=["geolocation"],
        )
        pg = ctx.new_page()
        pg.on("request", on_request)

        targets = [
            f"https://sports.{state}.betmgm.com/en/sports/baseball-23",
            f"https://sports.{state}.betmgm.com/en/sports",
        ]
        for url in targets:
            try:
                pg.goto(url, wait_until="domcontentloaded", timeout=60000)
            except Exception as e:
                print(f"  nav {url[-25:]} slow/err: {e!r}")
            pg.wait_for_timeout(12000)
            final = pg.frames[0].url
            print(f"  loaded {url[-30:]}  ->  landed on {final[:60]}")
            if "www.betmgm.com" in final and accessids == set():
                print("  !! Redirected to www.betmgm.com (geo router). This IP is NOT")
                print("     geolocating to a legal sports state — SPA won't mint an")
                print("     accessid here. Re-run from an in-state IP / VPN.")

        print(f"\n  accessids harvested: {list(accessids)[:3] or '(none)'}")
        print(f"  cds-api/bettingoffer calls seen: {len(cds_calls)}")
        for m, u, body in cds_calls[:15]:
            print(f"    {m} {u}")
            if body:
                print(f"        POST: {body}")
        print(f"\n  SGP/price-flagged calls: {len(flagged)}")
        for m, u, body in flagged[:25]:
            print(f"    *** {m} {u}")
            if body:
                print(f"        POST: {body}")

        out = {
            "accessids": list(accessids),
            "cds_calls": cds_calls,
            "flagged": flagged,
        }
        with open("recon_betmgm_sgp.json", "w") as f:
            json.dump(out, f, indent=2, default=str)
        print("\n  saved -> recon_betmgm_sgp.json")
        if accessids:
            aid = list(accessids)[0]
            print("\n  NEXT: python recon_betmgm_sgp.py pull "
                  f"--accessid '{aid}' --state {state}")
        b.close()


# ---------------------------------------------------------------------------
# Mode 3: exercise the read API + attempt the priced-slip combo endpoint
# ---------------------------------------------------------------------------
def _walk_options(fixture):
    """Yield (market_name, option_name, option_id, price_decimal) tuples."""
    for g in fixture.get("games", []):
        gname = g.get("name", {}).get("value", "") if isinstance(g.get("name"), dict) else g.get("name", "")
        for r in g.get("results", []):
            oid = r.get("id")
            oname = r.get("name", {}).get("value", "") if isinstance(r.get("name"), dict) else r.get("name", "")
            price = None
            odds = r.get("odds")
            if isinstance(odds, (int, float)):
                price = odds
            yield gname, oname, oid, price


def mode_pull(state: str, accessid: str):
    s = _session()
    s.headers["Referer"] = f"https://sports.{state}.betmgm.com/"

    print(f"=== PULL: read API + priced-slip attempt (state={state}) ===\n")
    url = fixtures_url(state, accessid)
    print(f"[a] GET MLB fixtures\n    {url[:110]}...")
    r = s.get(url, timeout=30)
    print(f"    -> HTTP {r.status_code} ({len(r.text)} bytes)")
    if r.status_code != 200:
        print(f"    -> body: {r.text[:300]}")
        print("    accessid likely invalid/expired — re-harvest. Aborting.")
        return
    try:
        data = r.json()
    except Exception:
        print("    -> non-JSON; aborting.")
        return

    fixtures = data.get("fixtures") or data.get("items") or []
    print(f"    -> {len(fixtures)} fixtures parsed\n")
    if not fixtures:
        print("    No fixtures (off-day or wrong sportId). Aborting.")
        return

    # Pick first fixture, find a spread (run line) option + a total option
    fx = fixtures[0]
    parts = fx.get("participants", [])
    name = " vs ".join(
        (pp.get("name", {}).get("value", "") if isinstance(pp.get("name"), dict) else str(pp.get("name", "")))
        for pp in parts
    ) or fx.get("name", {}).get("value", "?")
    fx_id = fx.get("id")
    print(f"[b] Fixture: {name}  (id={fx_id})")

    spread_opt = total_opt = None
    for gname, oname, oid, price in _walk_options(fx):
        gl = gname.lower()
        if spread_opt is None and ("run line" in gl or "handicap" in gl or "spread" in gl):
            spread_opt = (gname, oname, oid, price)
        if total_opt is None and ("total" in gl or "over/under" in gl or "over under" in gl):
            total_opt = (gname, oname, oid, price)
    print(f"    spread leg: {spread_opt}")
    print(f"    total  leg: {total_opt}\n")

    if not (spread_opt and total_opt):
        print("    Could not locate both a spread and a total option in fixture[0].")
        print("    (Inspect recon_betmgm_sgp.json fixtures dump to map market names.)")
        with open("recon_betmgm_sgp.json", "w") as f:
            json.dump(data, f, indent=2, default=str)
        return

    option_ids = [spread_opt[2], total_opt[2]]
    print(f"[c] Attempting priced-slip combo with optionIds={option_ids}\n")

    # Candidate bet-builder / priced-slip POST endpoints (CDS family + Angstrom).
    # The exact one is confirmed by the `harvest` capture; we try the known
    # shapes here so a valid accessid + in-state run nails it down.
    candidates = [
        (f"{cds_base(state)}/betslip/v1/quickbets",
         {"accessId": accessid, "selectionIds": option_ids, "wantsBetBuilder": True}),
        (f"{cds_base(state)}/bettingoffer/betbuilder/price",
         {"optionIds": option_ids}),
        (f"{cds_base(state)}/betslip/api/betbuilder",
         {"selections": [{"optionId": o} for o in option_ids]}),
        (f"{cds_base(state)}/sgp/v1/price",
         {"fixtureId": fx_id, "optionIds": option_ids}),
    ]
    for u, body in candidates:
        try:
            rr = s.post(u, json=body, timeout=20)
            tag = ""
            if rr.status_code == 200:
                tag = "  <<< SUCCESS — inspect body"
            print(f"    POST {u[len('https://'):]}\n      -> HTTP {rr.status_code} "
                  f"({len(rr.text)} bytes){tag}")
            if rr.text:
                print(f"      body: {rr.text[:240]}")
        except Exception as e:
            print(f"    POST {u[len('https://'):]} -> ERR {e!r}")

    print("\n  If all candidates 404/400: open recon_betmgm_sgp.json from the")
    print("  `harvest` run — the real SGP POST path is captured there under")
    print("  'flagged' (look for betslip / price / samegame).")


# ---------------------------------------------------------------------------
# Mode 4: THE WORKING PATH — headless accessid harvest + fixtures + SGP price
# ---------------------------------------------------------------------------
def harvest_accessid(state: str) -> str:
    """Pull the brand-correct accessid from clientconfig — no browser.

    The `x-bwin-sports-api: prod` header is the unlock: without it the
    endpoint returns the host-app shell (no accessid) even from CA.
    """
    host = f"www.{state}.betmgm.com"
    burl = f"http%3A%2F%2F{host}%2Fen%2Fsports%2Fbaseball-23"
    url = (f"https://{host}/en/api/clientconfig"
           f"?browserUrl={burl}&x-from-product=host-app")
    s = _session()
    r = s.get(url, headers={
        "Referer": f"https://{host}/en/sports/baseball-23",
        "x-bwin-browser-url": burl,
        "x-from-product": "host-app",
        "x-bwin-sports-api": "prod",
    }, timeout=25)
    ids = re.findall(r'"publicAccessId":"([^"]+)"', r.text)
    return ids[0] if ids else ""


def mode_mgm(state: str):
    import uuid
    print(f"=== MGM end-to-end (state={state}, from California) ===\n")

    aid = harvest_accessid(state)
    print(f"[1] accessid (headless harvest): {aid or '(EMPTY — try pa/mi/co/va/tn)'}")
    if not aid:
        print("    NJ omits pubid in config; use a state that returns one.")
        return

    s = _session()
    s.headers["Referer"] = f"https://www.{state}.betmgm.com/"
    base = f"https://www.{state}.betmgm.com/cds-api"
    fx_url = (f"{base}/bettingoffer/fixtures?x-bwin-accessid={aid}"
              f"&lang=en&country=US&userCountry=US&fixtureTypes=Standard&state=Latest"
              f"&offerMapping=Filtered&offerCategories=Gridable"
              f"&fixtureCategories=Gridable&sportIds={MLB_SPORT_ID}"
              f"&skip=0&take=10&sortBy=Tags")
    r = s.get(fx_url, timeout=30)
    print(f"[2] fixtures -> HTTP {r.status_code} ({len(r.text)} bytes)")
    if r.status_code != 200:
        print(f"    body: {r.text[:160]}"); return
    fixtures = r.json().get("fixtures", [])
    if not fixtures:
        print("    no MLB fixtures (off-day?)"); return

    # find a fixture exposing both a run-line spread and a totals market
    def markets(fx):
        out = {}
        for om in fx.get("optionMarkets", []):
            nm = om.get("name", {}).get("value", "")
            if nm == "Run Line Spread" and "spread" not in out:
                out["spread"] = om
            if nm == "Totals" and "total" not in out:
                out["total"] = om
        return out

    fx = sp = tot = None
    for cand in fixtures:
        m = markets(cand)
        if "spread" in m and "total" in m:
            fx, sp, tot = cand, m["spread"], m["total"]; break
    if not fx:
        print("    no fixture with both spread+total markets"); return

    name = fx.get("name", {}).get("value", "?")
    fx_id = str(fx.get("id"))
    sp_opt = sp["options"][0]   # first spread side
    tot_opt = tot["options"][0]  # Over
    print(f"[3] fixture: {name} (id={fx_id})")
    print(f"    spread leg: market={sp['id']} option={sp_opt['id']} "
          f"'{sp_opt['name']['value']}' @ {sp_opt['price']['odds']}")
    print(f"    total  leg: market={tot['id']} option={tot_opt['id']} "
          f"'{tot_opt['name']['value']}' @ {tot_opt['price']['odds']}\n")

    pg = str(uuid.uuid4())
    body = {"tv1Picks": [
        {"fixtureId": fx_id, "gameId": sp["id"], "resultId": sp_opt["id"],
         "useLiveFallback": False, "pickGroupId": pg},
        {"fixtureId": fx_id, "gameId": tot["id"], "resultId": tot_opt["id"],
         "useLiveFallback": False, "pickGroupId": pg},
    ], "tv2Picks": []}
    pick_url = (f"{base}/bettingoffer/picks?x-bwin-accessid={aid}"
                f"&lang=en&country=US&userCountry=US")
    rr = s.post(pick_url, json=body, timeout=30)
    print(f"[4] POST bettingoffer/picks -> HTTP {rr.status_code}")
    groups = rr.json().get("betBuilderPricingGroups", {}) if rr.status_code == 200 else {}
    if not groups:
        print(f"    no pricing group returned. body: {rr.text[:200]}"); return
    grp = next(iter(groups.values()))
    combo = grp["odds"]["odds"]
    naive = round(sp_opt["price"]["odds"] * tot_opt["price"]["odds"], 3)
    print(f"    CORRELATED SGP price: {combo}  (americanOdds {grp['odds']['americanOdds']}, "
          f"provider {grp.get('providerId')})")
    print(f"    naive product of legs: {naive}  -> "
          f"correlation diff {(combo/naive - 1) * 100:+.1f}%")
    print("\n  => GREEN: fixtures + correlated SGP price pulled from California.")


# ---------------------------------------------------------------------------
def main():
    ap = argparse.ArgumentParser(description="BetMGM MLB SGP recon")
    ap.add_argument("mode", choices=["probe", "harvest", "pull", "mgm"])
    ap.add_argument("--state", default="nj", help="betmgm state subdomain (nj, pa, mi, ...)")
    ap.add_argument("--accessid", default="", help="x-bwin-accessid for pull mode")
    ap.add_argument("--headful", action="store_true", help="show the browser (harvest mode)")
    args = ap.parse_args()

    if args.mode == "probe":
        mode_probe(args.state)
    elif args.mode == "harvest":
        mode_harvest(args.state, args.headful)
    elif args.mode == "pull":
        if not args.accessid:
            print("pull mode requires --accessid (get one from `harvest`).")
            sys.exit(1)
        mode_pull(args.state, args.accessid)
    elif args.mode == "mgm":
        mode_mgm("pa" if args.state == "nj" else args.state)


if __name__ == "__main__":
    main()
