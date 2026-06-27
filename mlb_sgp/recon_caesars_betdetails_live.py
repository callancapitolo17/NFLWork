#!/usr/bin/env python3
"""
LIVE TEST of the Caesars BYO/SGP pricer: POST /v2/bets/details.

From bundle static analysis the combined SGP price is an RTK-Query mutation
`getBetDetails` -> service `getBetDetails`:
    method POST  url /v2/bets/details
    headers {PlayerId, SessionId}   (logged-out => empty)
    body {legs:[...], combinationSelections:[], channel, channelDetail:"cordova-desktop",
          includeMultiLine, walletType?}
    skipContentType:true
The base host is api.americanwagering.com brand path (same as the event feed).
This is NOT the be-push socket.io (that socket is bet-referral '{acct}:RefBet',
login-only). So arbitrary-combo pricing should be a pure REST call.

This script:
  1. warm browser -> mint WAF token + cookies + device id; capture the REAL
     /v2/bets/details (or /bets/details) request the SPA fires (url+method+
     headers+postData) if any betslip activity occurs.
  2. harvest a live MLB event -> run-line selection + total selection (ids/lines/d).
  3. POST /v2/bets/details over curl_cffi (token + app headers + real UA),
     logged-out, building the leg payload per the bundle shape. Try the BYO
     `type:"BYO"` single combinationSelection AND the plain 2-leg form.
  4. read the combined decimal from the response; compare to naive product.

Output -> recon_caesars_betdetails_live.json (+ _czr_bd_*.json bodies)
"""
import json, os, re

STATE = os.environ.get("CAESARS_STATE", "nj")
PROFILE = os.environ.get("CAESARS_PROFILE", "/tmp/czr_prof_bd")
API = "https://api.americanwagering.com"
BASE = f"{API}/regions/us/locations/{STATE}/brands/czr"
SB = f"{BASE}/sb"
OUT = os.path.join(os.path.dirname(__file__), "recon_caesars_betdetails_live.json")
REAL_UA = ("Mozilla/5.0 (Macintosh; Intel Mac OS X 10_15_7) AppleWebKit/537.36 "
           "(KHTML, like Gecko) Chrome/149.0.0.0 Safari/537.36")
out = {"state": STATE}


def warm_and_sniff():
    """Mint token+cookies+device; sniff any real betslip /bets/details call."""
    from playwright.sync_api import sync_playwright
    sniff = {"calls": []}
    dev = {"v": None}
    with sync_playwright() as pw:
        ctx = pw.chromium.launch_persistent_context(
            user_data_dir=PROFILE, channel="chrome", headless=True,
            viewport={"width": 1400, "height": 900}, user_agent=REAL_UA,
            args=["--disable-blink-features=AutomationControlled", f"--user-agent={REAL_UA}"])
        ctx.add_init_script("Object.defineProperty(navigator,'webdriver',{get:()=>undefined});")
        page = ctx.pages[0] if ctx.pages else ctx.new_page()

        def on_req(req):
            u = req.url
            if "americanwagering.com" in u:
                d = req.headers.get("x-unique-device-id")
                if d:
                    dev["v"] = d
                if "bets/details" in u or "/bets/transform" in u:
                    pd = None
                    try:
                        pd = req.post_data
                    except Exception:
                        pass
                    sniff["calls"].append({
                        "url": u, "method": req.method,
                        "headers": {k: v for k, v in req.headers.items()
                                    if k.lower() in ("playerid", "sessionid", "x-platform",
                                                     "x-app-version", "content-type",
                                                     "x-aws-waf-token", "x-unique-device-id")},
                        "post_data": pd[:1500] if pd else None})
        page.on("request", on_req)
        page.goto(f"https://sportsbook.caesars.com/us/{STATE}/bet/baseball/competition/mlb-04f8d",
                  wait_until="commit", timeout=90000)
        tok = ""
        for _ in range(16):
            page.wait_for_timeout(2500)
            c = next((c for c in ctx.cookies() if c["name"] == "aws-waf-token"), None)
            if c:
                tok = c["value"]
        cookies = {c["name"]: c["value"] for c in ctx.cookies()
                   if "americanwagering" in c["domain"] or c["name"] == "aws-waf-token"}
        ctx.close()
    return tok, cookies, dev["v"], sniff["calls"]


def harvest_legs(get):
    """Return seed dict with run-line + total selection ids/lines/decimals."""
    st, txt = get(f"{SB}/v4/sports/baseball/tabs/SCHEDULE%7CGames%20%E2%9A%BE")
    if st != 200 or not txt.startswith("{"):
        return None, st
    tj = json.loads(txt)
    for comp in tj.get("competitions", []):
        cid = comp.get("id")
        for ev in comp.get("events", []):
            eid = ev.get("id")
            rl = tot = None
            for kmg in ev.get("keyMarketGroups", []):
                for m in kmg.get("markets", []):
                    nm = (m.get("name") or "").strip("|").lower()
                    sels = m.get("selections") or []
                    if not sels:
                        continue
                    if ("run line" in nm or "handicap" in nm) and not rl:
                        rl = (m, sels[0])
                    if "total" in nm and "team" not in nm and not tot:
                        tot = (m, sels[0])
            if eid and rl and tot:
                md = ev.get("metadata", {})
                def leg(mm, ss, signflip=False):
                    line = mm.get("line")
                    return {"selectionId": ss["id"], "marketId": mm["id"], "eventId": eid,
                            "competitionId": cid, "line": line,
                            "name": ss.get("name"), "selectionType": ss.get("type") or "",
                            "price": ss.get("price") or {},
                            "d": (ss.get("price") or {}).get("d")}
                return {"competitionId": cid, "eventId": eid,
                        "game": f'{md.get("awayTeamName")} @ {md.get("homeTeamName")}',
                        "rl": leg(rl[0], rl[1]), "tot": leg(tot[0], tot[1])}, st
    return None, st


def build_legs(seed):
    """Build the legs[] array in the bundle's shape (priceType fp)."""
    legs = []
    for k in ("rl", "tot"):
        s = seed[k]
        legs.append({
            "selectionId": s["selectionId"], "eventId": s["eventId"],
            "marketId": s["marketId"], "competitionId": s["competitionId"],
            "priceType": "fp", "related": True, "eachWay": False,
            "stakePerLine": 0, "line": s["line"],
            "selectionType": s["selectionType"], "name": s["name"],
            "price": s["price"], "origPrice": s["price"],
        })
    return legs


def main():
    from curl_cffi import requests as cr
    tok, cookies, dev, sniffed = warm_and_sniff()
    out["token_len"] = len(tok or "")
    out["device_id"] = dev
    out["sniffed_betdetails_calls"] = sniffed
    print(f"token_len={len(tok or '')} dev={dev} sniffed={len(sniffed)}")
    for c in sniffed:
        print("  SNIFF", c["method"], c["url"][:80], "PlayerId=", c["headers"].get("playerid"))

    if not tok:
        json.dump(out, open(OUT, "w"), indent=2, default=str); print("no token"); return

    hdrs = {"User-Agent": REAL_UA, "Accept": "application/json, text/plain, */*",
            "Content-Type": "application/json",
            "Origin": "https://sportsbook.caesars.com",
            "Referer": "https://sportsbook.caesars.com/",
            "X-Aws-Waf-Token": tok, "x-platform": "cordova-desktop",
            "x-app-version": "7.49.0"}
    if dev:
        hdrs["x-unique-device-id"] = dev

    def get(url):
        try:
            r = cr.get(url, headers={**hdrs, "Content-Type": "application/json"},
                       cookies=cookies, impersonate="chrome", timeout=25)
            return r.status_code, r.text
        except Exception as e:
            return f"ERR:{e}", ""

    def post(url, body):
        try:
            r = cr.post(url, headers=hdrs, cookies=cookies, json=body,
                        impersonate="chrome", timeout=30)
            return r.status_code, r.text
        except Exception as e:
            return f"ERR:{e}", ""

    seed, tabst = harvest_legs(get)
    out["tabs_status"] = tabst
    out["seed"] = seed
    if not seed:
        json.dump(out, open(OUT, "w"), indent=2, default=str); print("no seed"); return
    print("GAME:", seed["game"])
    print("  RL :", seed["rl"]["name"], seed["rl"]["line"], "@", seed["rl"]["d"])
    print("  TOT:", seed["tot"]["name"], seed["tot"]["line"], "@", seed["tot"]["d"])
    naive = None
    if seed["rl"]["d"] and seed["tot"]["d"]:
        naive = round(seed["rl"]["d"] * seed["tot"]["d"], 4)
        out["naive_product"] = naive
        print("  naive RLxTOT =", naive)

    legs = build_legs(seed)
    # combinationSelections: a single BYO across the 2 legs
    byo_combo = {"type": "BYO", "combinationType": "BYO",
                 "legs": [{"selectionId": l["selectionId"], "eventId": l["eventId"],
                           "marketId": l["marketId"]} for l in legs],
                 "stakePerLine": 0}
    common = {"channel": "desktop", "channelDetail": "cordova-desktop",
              "includeMultiLine": False}

    bodies = {
        # A. v2 with BYO combinationSelection
        "A_v2_byo": {**common, "legs": legs, "combinationSelections": [byo_combo]},
        # B. v2 plain legs, empty combinationSelections (server may infer SGP)
        "B_v2_plain": {**common, "legs": legs, "combinationSelections": []},
        # C. v1 bets/details with BYO
        "C_v1_byo": {**common, "legs": legs, "combinationSelections": [byo_combo]},
        # D. minimal selectionIds-only (some BE accept thin form)
        "D_v2_thin": {**common,
                      "legs": [{"selectionId": l["selectionId"], "eventId": l["eventId"],
                                "marketId": l["marketId"], "priceType": "fp",
                                "line": l["line"]} for l in legs],
                      "combinationSelections": [byo_combo]},
    }
    urls = {"A_v2_byo": f"{BASE}/{ 'sb/v2/bets/details' }",
            "B_v2_plain": f"{BASE}/sb/v2/bets/details",
            "C_v1_byo": f"{BASE}/sb/bets/details",
            "D_v2_thin": f"{BASE}/sb/v2/bets/details"}
    # also try without the /sb segment (path map uses brand base directly)
    alt_urls = {k: v.replace("/sb/", "/") for k, v in urls.items()}

    out["attempts"] = []
    for name, body in bodies.items():
        for label, u in ((name, urls[name]), (name + "_noSB", alt_urls[name])):
            st, txt = post(u, body)
            rec = {"name": label, "url": u, "status": st, "len": len(txt),
                   "sample": txt[:500]}
            if isinstance(st, int) and st == 200 and txt.strip().startswith(("{", "[")):
                open(os.path.join(os.path.dirname(__file__), f"_czr_bd_{label}.json"), "w").write(txt)
                rec["saved"] = f"_czr_bd_{label}.json"
                rec["combined_d"] = find_combined_decimal(json.loads(txt))
            out["attempts"].append(rec)
            cd = rec.get("combined_d")
            print(f"[{label}] {st} {len(txt)}B combined_d={cd}")

    json.dump(out, open(OUT, "w"), indent=2, default=str)
    print(f"\nsaved -> {OUT}")


def find_combined_decimal(j):
    """Walk the betDetails response for the combined parlay/BYO decimal price."""
    found = []

    def walk(o, path=""):
        if isinstance(o, dict):
            # price.d or price.decimal at a parlay/byo node
            p = o.get("price")
            d = None
            if isinstance(p, dict):
                d = p.get("d") if p.get("d") is not None else p.get("decimal")
            t = o.get("type") or o.get("combinationType")
            legs = o.get("legs") or o.get("betLegs")
            if d and isinstance(legs, list) and len(legs) >= 2:
                found.append({"path": path, "type": t, "d": d, "n_legs": len(legs)})
            for k, v in o.items():
                walk(v, f"{path}.{k}")
        elif isinstance(o, list):
            for i, v in enumerate(o):
                walk(v, f"{path}[{i}]")
    walk(j)
    return found


if __name__ == "__main__":
    main()
