#!/usr/bin/env python3
"""
Kalshi CBB 1H Market Maker — Performance Analysis

Pulls data from the Kalshi API (no DuckDB dependency) and prints a
formatted report covering:
  1. Settled P&L
  2. Closing Line Value (CLV)
  3. Maker vs Taker edge
  4. Fill pattern analysis
  5. Position concentration & risk

Usage:
    python3 analyze_performance.py
"""

import sys
import re
from pathlib import Path
from collections import defaultdict
from datetime import datetime

# Auth lives in kalshi_draft/ — resolve to the actual repo root
# (worktrees may not contain all directories)
SCRIPT_DIR = Path(__file__).parent
PROJECT_ROOT = SCRIPT_DIR.parent
_kalshi_draft = PROJECT_ROOT / "kalshi_draft"
if not _kalshi_draft.exists():
    import subprocess
    _main_root = Path(subprocess.check_output(
        ["git", "-C", str(SCRIPT_DIR), "rev-parse", "--path-format=absolute",
         "--git-common-dir"], text=True).strip()).parent
    _kalshi_draft = _main_root / "kalshi_draft"
sys.path.insert(0, str(_kalshi_draft))
from auth import load_credentials, authenticated_request, public_request

# ── credentials ──────────────────────────────────────────────────────
MM_DIR = SCRIPT_DIR
creds = load_credentials(env_dir=MM_DIR)
if not creds["api_key_id"]:
    creds = load_credentials(env_dir=_kalshi_draft.parent / "kalshi_mm")
API_KEY = creds["api_key_id"]
PK_PATH = creds["private_key_path"]
if not API_KEY or not PK_PATH:
    print("ERROR: Kalshi credentials not found.")
    print("  Checked .env in:", MM_DIR)
    print("  Set KALSHI_API_KEY_ID and KALSHI_PRIVATE_KEY_PATH in kalshi_mm/.env")
    sys.exit(1)

CBB_PREFIX = "KXNCAAMB"


# ── helpers ──────────────────────────────────────────────────────────
def api(path):
    return authenticated_request("GET", path, API_KEY, PK_PATH)


def paginate_auth(base_path):
    """Paginate an authenticated Kalshi list endpoint."""
    cursor = None
    sep = "&" if "?" in base_path else "?"
    while True:
        url = f"{base_path}{sep}cursor={cursor}" if cursor else base_path
        data = api(url)
        if not data:
            break
        items = None
        for key in ("fills", "market_positions", "event_positions", "markets",
                     "orders", "settlements"):
            if key in data:
                items = data[key]
                break
        if items is None:
            for v in data.values():
                if isinstance(v, list):
                    items = v
                    break
        if not items:
            break
        yield from items
        cursor = data.get("cursor")
        if not cursor:
            break


def parse_ticker(ticker):
    """Extract market_type, game_slug, and strike from a ticker."""
    parts = ticker.split("-", 2)
    prefix = parts[0]
    game = parts[1] if len(parts) > 1 else ""
    strike = parts[2] if len(parts) > 2 else ""
    if "SPREAD" in prefix:
        mtype = "spread"
    elif "TOTAL" in prefix:
        mtype = "total"
    elif "WINNER" in prefix:
        mtype = "moneyline"
    else:
        mtype = "other"
    return mtype, game, strike


def parse_game_teams(game_slug):
    """Extract a readable game label from the slug."""
    m = re.match(r"\d+[A-Z]{3}\d+(.+)", game_slug)
    return m.group(1) if m else game_slug


def dollars(val):
    if isinstance(val, str):
        return float(val)
    return float(val) if val else 0.0


def fmt_dollars(val):
    return f"${val:,.2f}" if val >= 0 else f"-${abs(val):,.2f}"


def pct(num, denom):
    return num / denom * 100 if denom else 0.0


def contracts(fill):
    return int(float(fill.get("count_fp", 1)))


def yes_price_cents(fill):
    return dollars(fill.get("yes_price_dollars", 0)) * 100


def fill_cost_cents(fill):
    """Actual cost paid for the fill (yes_price for YES buys, no_price for NO buys)."""
    if fill["side"] == "yes":
        return dollars(fill.get("yes_price_dollars", 0)) * 100
    else:
        return dollars(fill.get("no_price_dollars", 0)) * 100


def fill_pnl_cents(fill, settle_result):
    """P&L per contract in cents given settlement result ('yes' or 'no').
    All fills are action='buy'. Side determines what was bought.
    """
    settle_val = 100 if settle_result == "yes" else 0
    if fill["side"] == "yes":
        # Paid yes_price, receive 100 if YES, 0 if NO
        return settle_val - dollars(fill["yes_price_dollars"]) * 100
    else:
        # Paid no_price, receive 100 if NO, 0 if YES
        no_price = dollars(fill["no_price_dollars"]) * 100
        return (100 - settle_val) - no_price


# ── data pull ────────────────────────────────────────────────────────
def pull_data():
    print("Pulling data from Kalshi API...")

    balance = api("/portfolio/balance") or {}
    bal_cash = balance.get("balance", 0) / 100
    bal_portfolio = balance.get("portfolio_value", 0) / 100

    # All fills
    all_fills = list(paginate_auth("/portfolio/fills?limit=200"))
    cbb_fills = [f for f in all_fills if CBB_PREFIX in f.get("ticker", "")]

    # All positions (single pull)
    all_market_pos = []
    all_event_pos = []
    data = api("/portfolio/positions?limit=200")
    if data:
        all_market_pos = data.get("market_positions", [])
        all_event_pos = data.get("event_positions", [])
        cursor = data.get("cursor")
        while cursor:
            data2 = api(f"/portfolio/positions?limit=200&cursor={cursor}")
            if not data2:
                break
            all_market_pos.extend(data2.get("market_positions", []))
            all_event_pos.extend(data2.get("event_positions", []))
            cursor = data2.get("cursor")

    all_market_pos = [p for p in all_market_pos if CBB_PREFIX in p.get("ticker", "")]
    all_event_pos = [p for p in all_event_pos if CBB_PREFIX in p.get("event_ticker", "")]

    # Check settlement status for each unique market ticker
    unique_tickers = {f["ticker"] for f in cbb_fills}
    unique_tickers.update(p["ticker"] for p in all_market_pos)
    print(f"  {len(cbb_fills)} CBB fills, {len(all_market_pos)} market positions, "
          f"{len(all_event_pos)} events, {len(unique_tickers)} unique tickers")

    print("  Checking market settlement status...")
    market_results = {}  # ticker -> "yes" | "no" | None
    checked = 0
    for ticker in unique_tickers:
        mdata = public_request(f"/markets/{ticker}")
        if mdata and mdata.get("market"):
            mk = mdata["market"]
            result = mk.get("result") or None
            market_results[ticker] = result
        checked += 1
        if checked % 50 == 0:
            print(f"    ...checked {checked}/{len(unique_tickers)} markets")

    settled_count = sum(1 for r in market_results.values() if r is not None)
    print(f"  Settlement: {settled_count} settled, "
          f"{len(market_results) - settled_count} still open")

    return {
        "balance_cash": bal_cash,
        "balance_portfolio": bal_portfolio,
        "cbb_fills": cbb_fills,
        "all_market_pos": all_market_pos,
        "all_event_pos": all_event_pos,
        "market_results": market_results,
    }


# ── 1. Settled P&L ──────────────────────────────────────────────────
def report_settled_pnl(data):
    print("\n" + "=" * 60)
    print("  1. P&L REPORT")
    print("=" * 60)

    results = data["market_results"]
    fills = data["cbb_fills"]

    settled_tickers = {t for t, r in results.items() if r is not None}
    unsettled_fills = [f for f in fills if f["ticker"] not in settled_tickers]

    if not settled_tickers:
        print("\n  No markets have settled yet.")
        # Show open position summary
        _show_open_summary(data, fills)
        return

    # Compute P&L from fills on settled markets
    by_type = defaultdict(lambda: {"wins": 0, "losses": 0, "profit": 0,
                                    "cost": 0, "fees": 0, "contracts": 0})
    by_event = defaultdict(lambda: {"profit": 0, "cost": 0, "fees": 0})

    for f in fills:
        ticker = f["ticker"]
        if ticker not in settled_tickers:
            continue
        result = results[ticker]
        cnt = contracts(f)
        fee = dollars(f.get("fee_cost", 0))
        mtype, game, _ = parse_ticker(ticker)

        pnl_per = fill_pnl_cents(f, result)
        cost_per = fill_cost_cents(f)
        profit_dollars = pnl_per * cnt / 100
        cost_dollars = cost_per * cnt / 100

        d = by_type[mtype]
        d["profit"] += profit_dollars
        d["fees"] += fee
        d["cost"] += cost_dollars
        d["contracts"] += cnt
        if pnl_per > 0:
            d["wins"] += 1
        elif pnl_per < 0:
            d["losses"] += 1

        teams = parse_game_teams(game)
        ek = f"{teams} ({mtype})"
        by_event[ek]["profit"] += profit_dollars
        by_event[ek]["cost"] += cost_dollars
        by_event[ek]["fees"] += fee

    # Count unique settled events (games)
    settled_events = set()
    for f in fills:
        if f["ticker"] in settled_tickers:
            _, game, _ = parse_ticker(f["ticker"])
            settled_events.add(game)

    print(f"\n  SETTLED ({len(settled_tickers)} markets across {len(settled_events)} games)")
    print(f"  {'Type':<12} {'Fills':>6} {'Win':>5} {'Loss':>5} {'Cost':>10} "
          f"{'P&L':>10} {'Fees':>8} {'Net':>10} {'ROI':>8}")
    print("  " + "-" * 78)
    total = {"profit": 0, "fees": 0, "cost": 0, "fills": 0}
    for mtype in ("spread", "total", "moneyline", "other"):
        d = by_type.get(mtype)
        if not d:
            continue
        n = d["wins"] + d["losses"]
        net = d["profit"] - d["fees"]
        roi = pct(net, d["cost"]) if d["cost"] else 0
        print(f"  {mtype:<12} {n:>6} {d['wins']:>5} {d['losses']:>5} "
              f"{fmt_dollars(d['cost']):>10} {fmt_dollars(d['profit']):>10} "
              f"{fmt_dollars(d['fees']):>8} {fmt_dollars(net):>10} {roi:>7.1f}%")
        total["profit"] += d["profit"]
        total["fees"] += d["fees"]
        total["cost"] += d["cost"]
        total["fills"] += n
    net_total = total["profit"] - total["fees"]
    roi_total = pct(net_total, total["cost"]) if total["cost"] else 0
    print("  " + "-" * 78)
    print(f"  {'TOTAL':<12} {total['fills']:>6} {'':>5} {'':>5} "
          f"{fmt_dollars(total['cost']):>10} {fmt_dollars(total['profit']):>10} "
          f"{fmt_dollars(total['fees']):>8} {fmt_dollars(net_total):>10} {roi_total:>7.1f}%")

    # Best/worst events
    sorted_events = sorted(by_event.items(), key=lambda x: x[1]["profit"] - x[1]["fees"])
    if len(sorted_events) >= 2:
        best = sorted_events[-1]
        worst = sorted_events[0]
        print(f"\n  Best:  {best[0]} → net {fmt_dollars(best[1]['profit'] - best[1]['fees'])}")
        print(f"  Worst: {worst[0]} → net {fmt_dollars(worst[1]['profit'] - worst[1]['fees'])}")

    # Open positions summary
    if unsettled_fills:
        open_cost = sum(fill_cost_cents(f) * contracts(f) / 100 for f in unsettled_fills)
        print(f"\n  OPEN POSITIONS (unsettled)")
        print(f"    {len(unsettled_fills)} fills, {fmt_dollars(open_cost)} cost basis")


def _show_open_summary(data, fills):
    """Fallback when nothing has settled."""
    event_pos = data["all_event_pos"]
    total_exposure = sum(dollars(e.get("event_exposure_dollars", 0)) for e in event_pos)
    print(f"  Open exposure: {fmt_dollars(total_exposure)} across {len(event_pos)} events")


# ── 2. CLV Analysis ─────────────────────────────────────────────────
def report_clv(data):
    print("\n" + "=" * 60)
    print("  2. SETTLEMENT P&L BY MAKER / TAKER")
    print("=" * 60)

    fills = data["cbb_fills"]
    results = data["market_results"]
    settled_tickers = {t for t, r in results.items() if r is not None}

    if not settled_tickers:
        print("\n  No markets settled yet — true CLV requires settlement outcomes.")
        print("  Once games finish, this section will show:")
        print("    - Win rate by maker/taker and market type")
        print("    - Average profit per contract")
        print("    - CLV distribution")
        print()

        # Show what we CAN report: entry price distribution
        maker_prices = [yes_price_cents(f) for f in fills if not f.get("is_taker", False)]
        taker_prices = [yes_price_cents(f) for f in fills if f.get("is_taker", False)]
        if maker_prices:
            avg_m = sum(maker_prices) / len(maker_prices)
            print(f"  Entry price snapshot (pre-settlement):")
            print(f"    Maker avg YES price: {avg_m:.1f}c ({len(maker_prices)} fills)")
        if taker_prices:
            avg_t = sum(taker_prices) / len(taker_prices)
            print(f"    Taker avg YES price: {avg_t:.1f}c ({len(taker_prices)} fills)")
        return

    # Compute CLV from settled markets
    clv_data = defaultdict(list)

    for f in fills:
        ticker = f["ticker"]
        if ticker not in settled_tickers:
            continue
        result = results[ticker]
        is_taker = f.get("is_taker", False)
        cnt = contracts(f)
        mtype, _, _ = parse_ticker(ticker)

        pnl_per = fill_pnl_cents(f, result)
        cost_per = fill_cost_cents(f)
        fee_dollars = dollars(f.get("fee_cost", 0))

        entry = {"profit_cents": pnl_per, "count": cnt,
                 "cost_cents": cost_per, "fee_dollars": fee_dollars,
                 "side": f["side"]}
        clv_data["maker" if not is_taker else "taker"].append(entry)
        clv_data[mtype].append(entry)
        clv_data["all"].append(entry)

    def clv_summary(label, entries):
        if not entries:
            return
        wins = sum(1 for e in entries if e["profit_cents"] > 0)
        total = len(entries)
        total_ct = sum(e["count"] for e in entries)
        weighted_profit = sum(e["profit_cents"] * e["count"] for e in entries)
        avg_profit = weighted_profit / max(total_ct, 1)
        total_profit = weighted_profit / 100
        total_fees = sum(e["fee_dollars"] for e in entries)
        net = total_profit - total_fees
        print(f"  {label:<12} {total:>5} fills  {total_ct:>6} cts  "
              f"Win: {pct(wins, total):>5.1f}%  "
              f"Avg: {avg_profit:>+6.1f}c/ct  "
              f"P&L: {fmt_dollars(total_profit):>10}  "
              f"Fees: {fmt_dollars(total_fees):>8}  "
              f"Net: {fmt_dollars(net):>10}")

    if not clv_data["all"]:
        print("\n  No fills matched to settled markets.\n")
        return

    print()
    clv_summary("ALL", clv_data["all"])
    print("  " + "-" * 80)
    clv_summary("Maker", clv_data["maker"])
    clv_summary("Taker", clv_data["taker"])
    print("  " + "-" * 80)
    clv_summary("Spread", clv_data["spread"])
    clv_summary("Total", clv_data["total"])
    clv_summary("Moneyline", clv_data["moneyline"])


# ── 3. Maker vs Taker ───────────────────────────────────────────────
def report_maker_taker(data):
    print("\n" + "=" * 60)
    print("  3. MAKER vs TAKER BREAKDOWN")
    print("=" * 60)

    fills = data["cbb_fills"]
    maker_fills = [f for f in fills if not f.get("is_taker", False)]
    taker_fills = [f for f in fills if f.get("is_taker", False)]

    def fill_stats(label, group):
        if not group:
            print(f"\n  {label}: No fills")
            return
        total_ct = sum(contracts(f) for f in group)
        total_cost = sum(fill_cost_cents(f) * contracts(f) / 100 for f in group)
        total_fees = sum(dollars(f.get("fee_cost", 0)) for f in group)
        avg_cost = total_cost / total_ct * 100 if total_ct else 0  # cents
        yes_count = sum(1 for f in group if f["side"] == "yes")
        no_count = len(group) - yes_count

        by_type = defaultdict(int)
        for f in group:
            mtype, _, _ = parse_ticker(f["ticker"])
            by_type[mtype] += 1

        print(f"\n  {label}:")
        print(f"    Fills: {len(group):,}  |  Contracts: {total_ct:,}  |  "
              f"Cost: {fmt_dollars(total_cost)}  |  Fees: {fmt_dollars(total_fees)}")
        print(f"    Avg cost/contract: {avg_cost:.1f}c  |  YES/NO split: {yes_count}/{no_count}")
        print(f"    By type: " + ", ".join(f"{k}: {v}" for k, v in sorted(by_type.items())))
        print(f"    Avg contracts/fill: {total_ct/len(group):.1f}")

    fill_stats("MAKER", maker_fills)
    fill_stats("TAKER", taker_fills)

    if maker_fills and taker_fills:
        maker_cost = sum(fill_cost_cents(f) * contracts(f) / 100 for f in maker_fills)
        taker_cost = sum(fill_cost_cents(f) * contracts(f) / 100 for f in taker_fills)
        print(f"\n  Maker/Taker cost ratio: {pct(maker_cost, maker_cost + taker_cost):.0f}% / "
              f"{pct(taker_cost, maker_cost + taker_cost):.0f}%")


# ── 4. Fill Patterns ────────────────────────────────────────────────
def report_fill_patterns(data):
    print("\n" + "=" * 60)
    print("  4. FILL PATTERN ANALYSIS")
    print("=" * 60)

    fills = data["cbb_fills"]

    maker_by_order = defaultdict(list)
    taker_by_order = defaultdict(list)
    for f in fills:
        if f.get("is_taker", False):
            taker_by_order[f["order_id"]].append(f)
        else:
            maker_by_order[f["order_id"]].append(f)

    def order_stats(label, orders_dict):
        if not orders_dict:
            print(f"\n  {label}: No orders")
            return

        fills_per_order = [len(v) for v in orders_dict.values()]
        contracts_per_order = [sum(contracts(f) for f in v) for v in orders_dict.values()]

        time_spreads = []
        for fills_list in orders_dict.values():
            if len(fills_list) < 2:
                continue
            times = sorted(f["created_time"] for f in fills_list)
            t0 = datetime.fromisoformat(times[0].replace("Z", "+00:00"))
            t1 = datetime.fromisoformat(times[-1].replace("Z", "+00:00"))
            time_spreads.append((t1 - t0).total_seconds())

        avg_fills = sum(fills_per_order) / len(fills_per_order)
        avg_ct = sum(contracts_per_order) / len(contracts_per_order)
        sorted_ct = sorted(contracts_per_order)
        median_ct = sorted_ct[len(sorted_ct) // 2]

        single_fill = sum(1 for x in fills_per_order if x == 1)
        multi_fill = len(fills_per_order) - single_fill

        print(f"\n  {label}:")
        print(f"    Unique orders: {len(orders_dict):,}")
        print(f"    Avg fills/order: {avg_fills:.1f}  |  "
              f"Single-fill: {single_fill} ({pct(single_fill, len(orders_dict)):.0f}%)  |  "
              f"Multi-fill: {multi_fill} ({pct(multi_fill, len(orders_dict)):.0f}%)")
        print(f"    Contracts/order — avg: {avg_ct:.1f}  |  median: {median_ct}")

        if time_spreads:
            avg_spread = sum(time_spreads) / len(time_spreads)
            max_spread = max(time_spreads)
            min_spread = min(time_spreads)
            print(f"    Multi-fill time spread: avg {avg_spread:.0f}s  |  "
                  f"min {min_spread:.0f}s  |  max {max_spread:.0f}s")

        size_buckets = {"1": 0, "2-5": 0, "6-10": 0, "11-25": 0, "26-50": 0, "50+": 0}
        for c in contracts_per_order:
            if c == 1:
                size_buckets["1"] += 1
            elif c <= 5:
                size_buckets["2-5"] += 1
            elif c <= 10:
                size_buckets["6-10"] += 1
            elif c <= 25:
                size_buckets["11-25"] += 1
            elif c <= 50:
                size_buckets["26-50"] += 1
            else:
                size_buckets["50+"] += 1
        print(f"    Size distribution: " +
              "  ".join(f"{k}: {v}" for k, v in size_buckets.items() if v > 0))

    order_stats("MAKER orders", maker_by_order)
    order_stats("TAKER orders", taker_by_order)

    # Red flags: large single-fill maker orders (above 75th percentile)
    print("\n  Adverse selection check:")
    maker_order_sizes = [sum(contracts(f) for f in v) for v in maker_by_order.values()]
    if maker_order_sizes:
        sorted_sizes = sorted(maker_order_sizes)
        p75 = sorted_sizes[int(len(sorted_sizes) * 0.75)]
        large_instant = 0
        large_instant_tickers = defaultdict(int)
        for oid, fills_list in maker_by_order.items():
            ct = sum(contracts(f) for f in fills_list)
            if ct > p75 and len(fills_list) == 1:
                large_instant += 1
                large_instant_tickers[fills_list[0]["ticker"]] += 1

        if large_instant:
            print(f"    {large_instant} maker orders above p75 ({p75} cts) filled in single trade")
            # Show top affected tickers
            top_tickers = sorted(large_instant_tickers.items(), key=lambda x: -x[1])[:5]
            for ticker, cnt in top_tickers:
                mtype, game, strike = parse_ticker(ticker)
                teams = parse_game_teams(game)
                print(f"      {teams} {mtype} {strike}: {cnt} times")
        else:
            print(f"    No large single-fill maker orders above p75 ({p75} cts)")


# ── 5. Position Concentration ───────────────────────────────────────
def report_concentration(data):
    print("\n" + "=" * 60)
    print("  5. POSITION CONCENTRATION & RISK")
    print("=" * 60)

    event_pos = data["all_event_pos"]
    market_pos = data["all_market_pos"]

    if not event_pos:
        print("\n  No event positions found.\n")
        return

    events = []
    for e in event_pos:
        cost = dollars(e.get("total_cost_dollars", 0))
        exposure = dollars(e.get("event_exposure_dollars", 0))
        rpnl = dollars(e.get("realized_pnl_dollars", 0))
        fees = dollars(e.get("fees_paid_dollars", 0))
        mtype, game, _ = parse_ticker(e["event_ticker"])
        teams = parse_game_teams(game)
        events.append({
            "event": e["event_ticker"], "teams": teams, "mtype": mtype,
            "cost": cost, "exposure": exposure, "rpnl": rpnl, "fees": fees,
        })

    events.sort(key=lambda x: x["exposure"], reverse=True)

    print(f"\n  Top 10 positions by exposure:")
    print(f"  {'Game':<25} {'Type':<10} {'Exposure':>10} {'Cost':>10} {'P&L':>10}")
    print("  " + "-" * 67)
    for e in events[:10]:
        print(f"  {e['teams']:<25} {e['mtype']:<10} {fmt_dollars(e['exposure']):>10} "
              f"{fmt_dollars(e['cost']):>10} {fmt_dollars(e['rpnl']):>10}")

    # By market type
    print(f"\n  Allocation by market type:")
    by_type = defaultdict(lambda: {"cost": 0, "exposure": 0, "count": 0})
    total_cost = total_exposure = 0
    for e in events:
        by_type[e["mtype"]]["cost"] += e["cost"]
        by_type[e["mtype"]]["exposure"] += e["exposure"]
        by_type[e["mtype"]]["count"] += 1
        total_cost += e["cost"]
        total_exposure += e["exposure"]

    print(f"  {'Type':<12} {'Events':>6} {'Exposure':>10} {'% Exp':>7} {'Cost':>10} {'% Cost':>7}")
    print("  " + "-" * 55)
    for mtype in ("spread", "total", "moneyline", "other"):
        d = by_type.get(mtype)
        if not d:
            continue
        print(f"  {mtype:<12} {d['count']:>6} {fmt_dollars(d['exposure']):>10} "
              f"{pct(d['exposure'], total_exposure):>6.1f}% "
              f"{fmt_dollars(d['cost']):>10} {pct(d['cost'], total_cost):>6.1f}%")
    print("  " + "-" * 55)
    print(f"  {'TOTAL':<12} {len(events):>6} {fmt_dollars(total_exposure):>10} "
          f"{'100.0%':>7} {fmt_dollars(total_cost):>10} {'100.0%':>7}")

    # Concentration
    exposures = sorted([e["exposure"] for e in events if e["exposure"] > 0], reverse=True)
    if exposures:
        print(f"\n  Concentration:")
        print(f"    Largest single event: {fmt_dollars(exposures[0])} "
              f"({pct(exposures[0], total_exposure):.1f}% of total)")
        if len(exposures) >= 3:
            top3 = sum(exposures[:3])
            print(f"    Top 3 events: {fmt_dollars(top3)} "
                  f"({pct(top3, total_exposure):.1f}% of total)")
        print(f"    Events with exposure: {len(exposures)} of {len(events)}")

    # Directional correlation
    print(f"\n  Directional correlation (spread + ML on same game):")
    games = defaultdict(lambda: {"spread": 0, "total": 0, "moneyline": 0})
    for p in market_pos:
        pos = float(p.get("position_fp", 0))
        if pos == 0:
            continue
        mtype, game, _ = parse_ticker(p["ticker"])
        games[game][mtype] += pos

    correlated = []
    for game, positions in games.items():
        sp, ml = positions["spread"], positions["moneyline"]
        if sp != 0 and ml != 0:
            same = (sp > 0) == (ml > 0)
            teams = parse_game_teams(game)
            correlated.append((teams, sp, ml, same))

    if correlated:
        for teams, sp, ml, same in correlated:
            flag = "SAME dir" if same else "opposite"
            print(f"    {teams}: spread={sp:+.0f} ML={ml:+.0f} ({flag})")
    else:
        print(f"    No games with both spread + ML positions")

    # Worst case
    if events:
        worst = max(events, key=lambda x: x["exposure"])
        print(f"\n  Worst-case single-event loss: {fmt_dollars(worst['exposure'])} "
              f"({worst['teams']} {worst['mtype']})")


# ── Header ───────────────────────────────────────────────────────────
def print_header(data):
    print("\n" + "=" * 60)
    print("  KALSHI CBB 1H MARKET MAKER — PERFORMANCE REPORT")
    print(f"  Generated: {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}")
    print("=" * 60)

    cash = data["balance_cash"]
    portfolio = data["balance_portfolio"]
    print(f"\n  Account: {fmt_dollars(cash)} cash + {fmt_dollars(portfolio)} portfolio "
          f"= {fmt_dollars(cash + portfolio)} total")

    n_fills = len(data["cbb_fills"])
    n_events = len(data["all_event_pos"])
    n_open = sum(1 for p in data["all_market_pos"] if float(p.get("position_fp", 0)) != 0)
    total_cost = sum(fill_cost_cents(f) * contracts(f) / 100 for f in data["cbb_fills"])
    total_fees = sum(dollars(f.get("fee_cost", 0)) for f in data["cbb_fills"])
    total_contracts = sum(contracts(f) for f in data["cbb_fills"])

    print(f"  Fills: {n_fills:,}  |  Contracts: {total_contracts:,}  |  "
          f"Open events: {n_events}  |  Open positions: {n_open}")
    print(f"  Total cost: {fmt_dollars(total_cost)}  |  Fees: {fmt_dollars(total_fees)}")


# ── main ─────────────────────────────────────────────────────────────
def main():
    data = pull_data()
    print_header(data)
    report_settled_pnl(data)
    report_clv(data)
    report_maker_taker(data)
    report_fill_patterns(data)
    report_concentration(data)
    print("\n" + "=" * 60)
    print("  END OF REPORT")
    print("=" * 60 + "\n")


if __name__ == "__main__":
    main()
