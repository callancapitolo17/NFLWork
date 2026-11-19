"""One-off READ-ONLY recon of FanDuel NFL markets for the SGP compare tool.

Run from the repo root on a machine that can reach FanDuel:
    python sgp_compare/recon/fanduel/recon_fd_nfl.py

Calls: event scan (NFL list), event-page for the soonest game (default tab, then
every tab the layout lists), and a handful of implyBets PRICE requests
(implyBets only quotes; it cannot place a bet). Writes gzipped raw JSON and a
summary.txt into sgp_compare/recon/fanduel/out/. No DB writes, no login.
"""
from __future__ import annotations

import gzip
import json
import sys
import time
from collections import defaultdict
from pathlib import Path

REPO = Path(__file__).resolve().parents[3]
sys.path.insert(0, str(REPO / "mlb_sgp"))
sys.path.insert(0, str(REPO))
from curl_cffi import requests as cffi_requests  # noqa: E402
from scraper_fanduel_sgp import FD_AK, FD_EVENT_PAGE_URL, FD_HEADERS, FD_IMPLY_BETS_URL, FD_SCAN_URL  # noqa: E402

OUT = Path(__file__).resolve().parent / "out"
OUT.mkdir(exist_ok=True)
NFL_COMPETITION_ID = 12282733
S = cffi_requests.Session(impersonate="chrome")
summary: list[str] = []


def note(line: str) -> None:
    print(line)
    summary.append(line)


def save(name: str, data) -> None:
    with gzip.open(OUT / f"{name}.json.gz", "wt") as f:
        json.dump(data, f)


def walk(node, keep):
    found, stack = [], [node]
    while stack:
        item = stack.pop()
        if isinstance(item, dict):
            if keep(item):
                found.append(item)
            stack.extend(item.values())
        elif isinstance(item, list):
            stack.extend(item)
    return found


def scan():
    body = {"filter": {"competitionIds": [NFL_COMPETITION_ID],
                       "contentGroup": {"language": "en", "regionCode": "NAMERICA"},
                       "marketLevels": ["AVB_EVENT"], "maxResults": 100,
                       "productTypes": ["SPORTSBOOK"], "selectBy": "FIRST_TO_START"},
            "facets": [{"type": "EVENT"}], "currencyCode": "USD"}
    r = S.post(FD_SCAN_URL, headers={**FD_HEADERS, "Content-Type": "application/json"}, data=json.dumps(body), timeout=30)
    note(f"scan HTTP {r.status_code}")
    data = r.json()
    save("scan", data)
    return data


def event_page(event_id: str, tab: str | None):
    url = f"{FD_EVENT_PAGE_URL}?_ak={FD_AK}&eventId={event_id}&useCombinedTouchdownsVirtualMarket=true&useQuickBets=true"
    if tab:
        url += f"&tab={tab}"
    r = S.get(url, headers=FD_HEADERS, timeout=30)
    note(f"event-page tab={tab!r} HTTP {r.status_code} bytes={len(r.content)}")
    return r.json() if r.status_code == 200 else None


def imply(refs, label):
    body = {"betLegs": [{"legType": "SIMPLE_SELECTION",
                         "betRunners": [{"runner": {"marketId": m, "selectionId": s}}]} for m, s in refs]}
    r = S.post(FD_IMPLY_BETS_URL, headers={**FD_HEADERS, "Content-Type": "application/json"}, data=json.dumps(body), timeout=20)
    data = r.json() if r.headers.get("content-type", "").startswith("application/json") else {"text": r.text[:2000]}
    save(f"imply_{label}", {"request": body, "status": r.status_code, "response": data})
    sgm = [c for c in data.get("betCombinations", []) if c.get("isSGM")] if isinstance(data, dict) else []
    note(f"implyBets {label}: HTTP {r.status_code} isSGM={bool(sgm)} keys={list(data)[:6] if isinstance(data, dict) else '?'}")
    if sgm:
        note(f"   SGM combo: {json.dumps(sgm[0])[:600]}")
    else:
        note(f"   body: {json.dumps(data)[:600]}")
    time.sleep(1)


def main():
    data = scan()
    events = walk(data, lambda d: "eventId" in d and "openDate" in d and "name" in d)
    events = sorted({e["eventId"]: e for e in events if " @ " in e["name"]}.values(), key=lambda e: e["openDate"])
    for e in events[:20]:
        note(f"  event {e['eventId']} {e['name']} {e['openDate']}")
    if not events:
        note("no NFL events found")
        return
    ev = events[0]
    eid = str(ev["eventId"])
    note(f"chosen event {eid} {ev['name']}")

    default = event_page(eid, None)
    save("event_default", default)
    tabs = ((default or {}).get("layout") or {}).get("tabs") or {}
    tab_keys = []
    for tid, t in tabs.items():
        note(f"  tab id={tid} title={t.get('title')!r} keys={list(t)[:8]}")
        tab_keys.append(t)
    # FD's tab query param is the title slug ('same-game-parlay-' style). Try slugs.
    slugs = []
    for t in tab_keys:
        title = (t.get("title") or "").strip()
        slug = "".join(ch.lower() if ch.isalnum() else "-" for ch in title)
        slugs.append(slug)
        if not slug.endswith("-"):
            slugs.append(slug + "-")
    slugs = list(dict.fromkeys(["same-game-parlay-", *slugs]))

    markets_all: dict[str, dict] = {}
    tab_of_market: dict[str, set] = defaultdict(set)
    for slug in slugs:
        page = event_page(eid, slug)
        time.sleep(0.7)
        if not page:
            continue
        ms = walk(page, lambda d: "marketId" in d and "runners" in d and "marketName" in d)
        if not ms:
            continue
        save(f"tab_{slug.strip('-') or 'x'}", page)
        note(f"  tab {slug}: {len(ms)} markets, {sum(1 for m in ms if m.get('sgmMarket'))} sgm")
        for m in ms:
            markets_all.setdefault(m["marketId"], m)
            tab_of_market[m["marketId"]].add(slug)

    note("\n=== MARKETS (sgm, status, type, #runners, tabs, sample runners) ===")
    by_type = defaultdict(list)
    for m in markets_all.values():
        by_type[m.get("marketType")].append(m)
    for mtype, ms in sorted(by_type.items(), key=lambda kv: str(kv[0])):
        m = ms[0]
        runners = m.get("runners", [])
        sample = [(r.get("runnerName"), r.get("handicap"), r.get("runnerStatus"),
                   ((r.get("winRunnerOdds") or {}).get("americanDisplayOdds") or {}).get("americanOdds"))
                  for r in runners[:4]]
        note(f"{mtype} x{len(ms)} | e.g. {m.get('marketName')!r} sgm={m.get('sgmMarket')} "
             f"status={m.get('marketStatus')} runners={len(runners)} tabs={sorted(tab_of_market[m['marketId']])} {sample}")
    note("\nmarket keys example: " + json.dumps({k: v for k, v in next(iter(markets_all.values())).items() if k != 'runners'})[:1500])
    note("runner keys example: " + json.dumps(next(iter(markets_all.values()))["runners"][0])[:1500])

    # --- price requests (quotes only) ---
    def first(pred):
        return next((m for m in markets_all.values() if m.get("sgmMarket") and m.get("marketStatus") == "OPEN" and pred(m)), None)

    def ref(m, i=0):
        rs = [r for r in m["runners"] if r.get("runnerStatus") == "ACTIVE"]
        return (m["marketId"], rs[i]["selectionId"]) if len(rs) > i else None

    ml = first(lambda m: m.get("marketName") == "Moneyline")
    tot = first(lambda m: m.get("marketName") == "Total Match Points" or m.get("marketType") == "TOTAL_POINTS_(OVER/UNDER)")
    spr = first(lambda m: m.get("marketName") == "Spread")
    prop_ou = first(lambda m: "Receiving Yds" in m.get("marketName", "") or "Receiving Yards" in m.get("marketName", ""))
    ladder = first(lambda m: "Alt Receiving" in m.get("marketName", "") or "Receiving Yds" in m.get("marketName", "") and len(m["runners"]) > 2)
    td = first(lambda m: "Touchdown" in m.get("marketName", "") and len(m["runners"]) > 2)
    note(f"\npicked: ml={ml and ml['marketName']} tot={tot and tot['marketName']} spr={spr and spr['marketName']} "
         f"prop={prop_ou and prop_ou['marketName']} ladder={ladder and ladder['marketName']} td={td and td['marketName']}")
    if ml and tot:
        imply([ref(ml, 0), ref(tot, 0)], "ml_total")
        imply([ref(ml, 0), ref(tot, 1)], "ml_under")
    if ml and spr:
        imply([ref(ml, 0), ref(spr, 0)], "ml_spread_same_team")
    if ml and tot and prop_ou:
        imply([ref(ml, 0), ref(tot, 0), ref(prop_ou, 0)], "ml_total_prop")
    if ml and td:
        imply([ref(ml, 0), ref(td, 0)], "ml_td")
    if ladder and ml:
        imply([ref(ml, 0), ref(ladder, 1)], "ml_ladder")
    if ml:
        imply([ref(ml, 0), ref(ml, 1)], "ml_both_sides_expect_decline")
    (OUT / "summary.txt").write_text("\n".join(summary))


if __name__ == "__main__":
    main()
