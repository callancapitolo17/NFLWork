#!/usr/bin/env python3
"""
Kalshi CBB 1H Market Maker — Performance Dashboard

Dash web app visualizing the same 5 analyses as analyze_performance.py.
Accessible from phone on same Wi-Fi at http://<mac-ip>:8084

Usage:
    python3 dashboard.py
"""

import dash
from dash import dcc, html, dash_table, callback, Input, Output, ctx
# Format objects don't serialize in Dash 4.0.0 — use raw dicts instead
FMT_2F = {"specifier": ".2f"}  # 2 decimal fixed
FMT_1F = {"specifier": ".1f"}  # 1 decimal fixed
import plotly.graph_objects as go
from collections import defaultdict
from datetime import datetime
import traceback

# Reuse all helpers from analyze_performance
from analyze_performance import (
    pull_data, parse_ticker, parse_game_teams, dollars, contracts,
    fill_cost_cents, fill_pnl_cents, pct, fmt_dollars,
)

# ── Theme (GitHub dark, matching kalshi_draft) ───────────────────────
COLORS = {
    "bg": "#0d1117",
    "card": "#161b22",
    "card_border": "#30363d",
    "text": "#c9d1d9",
    "text_muted": "#8b949e",
    "accent": "#00cec9",
    "accent2": "#0984e3",
    "green": "#3fb950",
    "red": "#f85149",
    "yellow": "#f0ad4e",
    "header_bg": "linear-gradient(135deg, #00b894 0%, #00cec9 50%, #0984e3 100%)",
}

CARD_STYLE = {
    "backgroundColor": COLORS["card"],
    "border": f"1px solid {COLORS['card_border']}",
    "borderRadius": "12px",
    "padding": "20px",
    "marginBottom": "16px",
}

TABLE_HEADER = {
    "backgroundColor": "#0d1b2a",
    "color": COLORS["accent"],
    "fontWeight": "700",
    "fontSize": "0.85em",
    "textTransform": "uppercase",
    "border": "none",
    "borderBottom": f"2px solid {COLORS['accent']}",
}

TABLE_DATA = {
    "backgroundColor": COLORS["card"],
    "color": COLORS["text"],
    "border": "none",
    "borderBottom": f"1px solid {COLORS['card_border']}",
}

TABLE_COND = [
    {"if": {"row_index": "odd"}, "backgroundColor": "#1c2333"},
]

GRAPH_LAYOUT = dict(
    template="plotly_dark",
    paper_bgcolor=COLORS["card"],
    plot_bgcolor=COLORS["card"],
    margin=dict(l=40, r=20, t=40, b=40),
    font=dict(color=COLORS["text"]),
)

# Pie chart color palette (enough for any market type count)
PIE_COLORS = [COLORS["accent"], COLORS["accent2"], COLORS["yellow"],
              COLORS["text_muted"], COLORS["green"]]

# CSS for dark-themed dropdown
DROPDOWN_CSS = """
#event-dropdown [class*="control"] {
    background-color: #161b22 !important;
    border-color: #30363d !important;
}
#event-dropdown [class*="menu"] {
    background-color: #161b22 !important;
    border-color: #30363d !important;
}
#event-dropdown [class*="option"] {
    background-color: #161b22 !important;
    color: #c9d1d9 !important;
}
#event-dropdown [class*="option"]:hover,
#event-dropdown [class*="option--is-focused"] {
    background-color: #30363d !important;
}
#event-dropdown [class*="singleValue"],
#event-dropdown [class*="input"] input {
    color: #c9d1d9 !important;
}
#event-dropdown [class*="placeholder"] {
    color: #8b949e !important;
}
#event-dropdown [class*="indicatorSeparator"] {
    background-color: #30363d !important;
}
"""


# ── Data store (refreshed on demand) ────────────────────────────────
DATA = {"loaded": False, "error": None}


def load_data():
    try:
        DATA.update(pull_data())
        DATA["loaded"] = True
        DATA["ts"] = datetime.now().strftime("%H:%M:%S")
        DATA["error"] = None
        _compute_all(DATA)
    except Exception as e:
        DATA["error"] = f"{type(e).__name__}: {e}"
        DATA["ts"] = datetime.now().strftime("%H:%M:%S")
        print(f"ERROR loading data: {DATA['error']}")
        traceback.print_exc()


def _compute_all(data):
    """Pre-compute aggregates for all tabs."""
    results = data.get("market_results", {})
    fills = data.get("cbb_fills", [])
    settled_tickers = {t for t, r in results.items() if r is not None}

    # --- Settled P&L by type ---
    by_type = defaultdict(lambda: {"wins": 0, "losses": 0, "pnl": 0.0,
                                    "cost": 0.0, "fees": 0.0, "fills": 0})
    by_event = defaultdict(lambda: {"pnl": 0.0, "cost": 0.0, "fees": 0.0, "fills": 0})

    for f in fills:
        t = f["ticker"]
        if t not in settled_tickers:
            continue
        result = results[t]
        cnt = contracts(f)
        fee = dollars(f.get("fee_cost", 0))
        mtype, game, _ = parse_ticker(t)
        pnl_per = fill_pnl_cents(f, result)
        cost_per = fill_cost_cents(f)

        d = by_type[mtype]
        d["pnl"] += pnl_per * cnt / 100
        d["cost"] += cost_per * cnt / 100
        d["fees"] += fee
        d["fills"] += 1
        if pnl_per > 0:
            d["wins"] += 1
        elif pnl_per < 0:
            d["losses"] += 1

        teams = parse_game_teams(game)
        ek = f"{teams} ({mtype})"
        by_event[ek]["pnl"] += pnl_per * cnt / 100
        by_event[ek]["cost"] += cost_per * cnt / 100
        by_event[ek]["fees"] += fee
        by_event[ek]["fills"] += 1

    data["pnl_by_type"] = [
        {"type": mt, **by_type[mt],
         "net": by_type[mt]["pnl"] - by_type[mt]["fees"],
         "roi": pct(by_type[mt]["pnl"] - by_type[mt]["fees"], by_type[mt]["cost"])
                if by_type[mt]["cost"] else 0}
        for mt in ("spread", "total", "moneyline", "other") if mt in by_type
    ]

    data["pnl_by_event"] = sorted([
        {"event": k, "pnl": v["pnl"], "cost": v["cost"], "fees": v["fees"],
         "net": v["pnl"] - v["fees"], "fills": v["fills"],
         "roi": pct(v["pnl"] - v["fees"], v["cost"]) if v["cost"] else 0}
        for k, v in by_event.items()
    ], key=lambda x: x["net"], reverse=True)

    # --- Maker vs Taker ---
    for role in ("maker", "taker"):
        group = [f for f in fills
                 if (not f.get("is_taker", False)) == (role == "maker")]
        settled = [f for f in group if f["ticker"] in settled_tickers]
        total_ct = sum(contracts(f) for f in group)
        total_cost = sum(fill_cost_cents(f) * contracts(f) / 100 for f in group)
        total_fees = sum(dollars(f.get("fee_cost", 0)) for f in group)
        settled_pnl = sum(fill_pnl_cents(f, results[f["ticker"]]) * contracts(f) / 100
                          for f in settled)
        settled_cost = sum(fill_cost_cents(f) * contracts(f) / 100 for f in settled)
        settled_fees = sum(dollars(f.get("fee_cost", 0)) for f in settled)
        settled_net = settled_pnl - settled_fees

        data[f"{role}_stats"] = {
            "fills": len(group), "contracts": total_ct,
            "cost": total_cost, "fees": total_fees,
            "settled_pnl": settled_pnl, "settled_cost": settled_cost,
            "settled_net": settled_net,
            "roi": pct(settled_net, settled_cost) if settled_cost else 0,
            "win_pct": pct(
                sum(1 for f in settled if fill_pnl_cents(f, results[f["ticker"]]) > 0),
                len(settled)) if settled else 0,
        }

    # Maker/taker by market type
    mt_rows = []
    for role in ("maker", "taker"):
        group = [f for f in fills
                 if (not f.get("is_taker", False)) == (role == "maker")
                 and f["ticker"] in settled_tickers]
        by_mt = defaultdict(lambda: {"pnl": 0, "cost": 0, "fees": 0, "fills": 0})
        for f in group:
            mtype, _, _ = parse_ticker(f["ticker"])
            cnt = contracts(f)
            by_mt[mtype]["pnl"] += fill_pnl_cents(f, results[f["ticker"]]) * cnt / 100
            by_mt[mtype]["cost"] += fill_cost_cents(f) * cnt / 100
            by_mt[mtype]["fees"] += dollars(f.get("fee_cost", 0))
            by_mt[mtype]["fills"] += 1
        for mt in ("spread", "total", "moneyline"):
            if mt in by_mt:
                d = by_mt[mt]
                net = d["pnl"] - d["fees"]
                mt_rows.append({"role": role.title(), "type": mt, "fills": d["fills"],
                                "cost": d["cost"], "net": net,
                                "roi": pct(net, d["cost"]) if d["cost"] else 0})
    data["maker_taker_by_type"] = mt_rows

    # --- Fill Patterns ---
    for role in ("maker", "taker"):
        group = [f for f in fills
                 if (not f.get("is_taker", False)) == (role == "maker")]
        by_order = defaultdict(list)
        for f in group:
            by_order[f["order_id"]].append(f)

        order_sizes = [sum(contracts(f) for f in v) for v in by_order.values()]
        single = sum(1 for s in [len(v) for v in by_order.values()] if s == 1)

        data[f"{role}_fill_stats"] = {
            "orders": len(by_order),
            "avg_size": sum(order_sizes) / len(order_sizes) if order_sizes else 0,
            "median_size": sorted(order_sizes)[len(order_sizes) // 2] if order_sizes else 0,
            "single_pct": pct(single, len(by_order)) if by_order else 0,
            "order_sizes": order_sizes,
        }

    # Adverse selection
    maker_orders = defaultdict(list)
    for f in fills:
        if not f.get("is_taker", False):
            maker_orders[f["order_id"]].append(f)
    order_sizes = [sum(contracts(f) for f in v) for v in maker_orders.values()]
    if order_sizes:
        p75 = sorted(order_sizes)[int(len(order_sizes) * 0.75)]
        adverse = []
        for oid, fl in maker_orders.items():
            ct = sum(contracts(f) for f in fl)
            if ct > p75 and len(fl) == 1:
                mtype, game, strike = parse_ticker(fl[0]["ticker"])
                is_settled = fl[0]["ticker"] in settled_tickers
                pnl = (fill_pnl_cents(fl[0], results.get(fl[0]["ticker"], ""))
                       * ct / 100 if is_settled else None)
                adverse.append({"game": parse_game_teams(game), "type": mtype,
                                "strike": strike, "contracts": ct,
                                "pnl": pnl,
                                "status": "Settled" if is_settled else "Open"})
        data["adverse_orders"] = sorted(adverse, key=lambda x: x["contracts"],
                                        reverse=True)
        data["adverse_p75"] = p75
    else:
        data["adverse_orders"] = []
        data["adverse_p75"] = 0

    # --- Fill Rate (maker orders) ---
    orders = data.get("cbb_orders", [])
    order_lookup = {o["order_id"]: o for o in orders}
    maker_order_ids = set(maker_orders.keys()) if order_sizes else set()

    fill_rate_rows = []
    for oid in maker_order_ids:
        order = order_lookup.get(oid)
        if not order:
            continue
        initial = int(float(order.get("initial_count_fp", 0)))
        filled = int(float(order.get("fill_count_fp", 0)))
        if initial <= 0 and filled <= 0:
            continue
        # initial_count_fp reflects last-amended size, not original.
        # Use max(initial, filled) as the effective order size so
        # amended-down orders cap at 100% instead of >100%.
        effective_size = max(initial, filled)
        fill_pct_val = min(filled / effective_size * 100, 100.0)
        mtype, game, _ = parse_ticker(order.get("ticker", ""))
        fill_rate_rows.append({
            "fill_pct": fill_pct_val,
            "original": effective_size,
            "filled": filled,
            "type": mtype,
        })

    # Per-fill bite size: each individual fill as % of resting order
    bite_pcts = []
    bite_by_type = defaultdict(list)
    for f in fills:
        if f.get("is_taker", False):
            continue
        order = order_lookup.get(f["order_id"])
        if not order:
            continue
        initial = int(float(order.get("initial_count_fp", 0)))
        total_filled = int(float(order.get("fill_count_fp", 0)))
        effective_size = max(initial, total_filled)
        if effective_size <= 0:
            continue
        bite = contracts(f) / effective_size * 100
        bite_pcts.append(min(bite, 100.0))
        mtype, _, _ = parse_ticker(f["ticker"])
        bite_by_type[mtype].append(min(bite, 100.0))

    if bite_pcts:
        sorted_bites = sorted(bite_pcts)
        data["bite_stats"] = {
            "avg": sum(bite_pcts) / len(bite_pcts),
            "median": sorted_bites[len(sorted_bites) // 2],
            "count": len(bite_pcts),
        }
        data["bite_by_type"] = [
            {"type": mt, "avg": sum(v) / len(v), "count": len(v)}
            for mt, v in sorted(bite_by_type.items())
        ]
    else:
        data["bite_stats"] = {"avg": 0, "median": 0, "count": 0}
        data["bite_by_type"] = []
    data["bite_pcts"] = bite_pcts

    data["fill_rate_rows"] = fill_rate_rows
    if fill_rate_rows:
        pcts = [r["fill_pct"] for r in fill_rate_rows]
        sorted_pcts = sorted(pcts)
        data["fill_rate_stats"] = {
            "avg": sum(pcts) / len(pcts),
            "median": sorted_pcts[len(sorted_pcts) // 2],
            "fully_filled": pct(sum(1 for p in pcts if p >= 100), len(pcts)),
            "partial": pct(sum(1 for p in pcts if 0 < p < 100), len(pcts)),
            "unfilled": pct(sum(1 for p in pcts if p == 0), len(pcts)),
            "count": len(pcts),
        }
        # By market type
        by_mt = defaultdict(list)
        for r in fill_rate_rows:
            by_mt[r["type"]].append(r["fill_pct"])
        data["fill_rate_by_type"] = [
            {"type": mt, "avg": sum(v) / len(v), "count": len(v),
             "fully_filled": pct(sum(1 for p in v if p >= 100), len(v))}
            for mt, v in sorted(by_mt.items())
        ]
    else:
        data["fill_rate_stats"] = {
            "avg": 0, "median": 0, "fully_filled": 0,
            "partial": 0, "unfilled": 0, "count": 0,
        }
        data["fill_rate_by_type"] = []

    # --- Fill Timing (time before market close) ---
    close_times = data.get("market_close_times", {})
    timing_buckets = defaultdict(lambda: {"fills": 0, "contracts": 0,
                                           "maker_fills": 0, "pnl": 0.0,
                                           "settled_fills": 0})
    BUCKET_EDGES = [0.25, 0.5, 1, 2, 4, 8, 24, float("inf")]
    BUCKET_LABELS = ["<15m", "15-30m", "30m-1h", "1-2h",
                     "2-4h", "4-8h", "8-24h", ">24h"]
    timing_hours = []

    for f in fills:
        t = f["ticker"]
        ct_str = close_times.get(t)
        if not ct_str:
            continue
        try:
            fill_dt = datetime.fromisoformat(
                f["created_time"].replace("Z", "+00:00"))
            close_dt = datetime.fromisoformat(ct_str.replace("Z", "+00:00"))
            hours_before = (close_dt - fill_dt).total_seconds() / 3600
        except (ValueError, TypeError):
            continue
        if hours_before < 0:
            continue
        cnt = contracts(f)
        is_maker = not f.get("is_taker", False)
        is_settled = t in settled_tickers
        fill_pnl = (fill_pnl_cents(f, results.get(t, "")) * cnt / 100
                     if is_settled else None)

        timing_hours.append(hours_before)

        # Find bucket
        for i, edge in enumerate(BUCKET_EDGES):
            if hours_before < edge:
                b = timing_buckets[BUCKET_LABELS[i]]
                b["fills"] += 1
                b["contracts"] += cnt
                if is_maker:
                    b["maker_fills"] += 1
                if fill_pnl is not None:
                    b["pnl"] += fill_pnl
                    b["settled_fills"] += 1
                break

    data["timing_buckets"] = {k: timing_buckets[k] for k in BUCKET_LABELS}
    data["timing_bucket_labels"] = BUCKET_LABELS
    data["timing_hours"] = timing_hours

    # --- Concentration ---
    event_pos = data.get("all_event_pos", [])
    conc_rows = []
    for e in event_pos:
        mtype, game, _ = parse_ticker(e["event_ticker"])
        conc_rows.append({
            "game": parse_game_teams(game), "type": mtype,
            "exposure": dollars(e.get("event_exposure_dollars", 0)),
            "cost": dollars(e.get("total_cost_dollars", 0)),
        })
    data["concentration"] = sorted(conc_rows, key=lambda x: x["exposure"],
                                   reverse=True)

    # Directional correlation
    market_pos = data.get("all_market_pos", [])
    games = defaultdict(lambda: {"spread": 0, "total": 0, "moneyline": 0})
    for p in market_pos:
        pos = float(p.get("position_fp", 0))
        if pos == 0:
            continue
        mtype, game, _ = parse_ticker(p["ticker"])
        games[game][mtype] += pos
    corr_rows = []
    for game, pos in games.items():
        if pos["spread"] != 0 and pos["moneyline"] != 0:
            same = (pos["spread"] > 0) == (pos["moneyline"] > 0)
            corr_rows.append({
                "game": parse_game_teams(game),
                "spread": pos["spread"], "ml": pos["moneyline"],
                "direction": "SAME" if same else "Opposite",
            })
    data["correlation"] = corr_rows

    # --- Event detail (all fills grouped by event) ---
    events_map = defaultdict(list)
    for f in fills:
        mtype, game, _ = parse_ticker(f["ticker"])
        events_map[f"{parse_game_teams(game)} ({mtype})"].append(f)
    data["events_map"] = dict(events_map)


# ── UI helpers ───────────────────────────────────────────────────────
def stat_card(label, value, color=None):
    return html.Div(style={**CARD_STYLE, "flex": "1", "minWidth": "140px",
                           "textAlign": "center"}, children=[
        html.Div(label, style={"color": COLORS["text_muted"], "fontSize": "0.75em",
                                "textTransform": "uppercase", "letterSpacing": "1px"}),
        html.Div(value, style={"color": color or COLORS["accent"],
                                "fontSize": "1.5em", "fontWeight": "700",
                                "marginTop": "4px"}),
    ])


def make_table(data_records, columns, id_prefix, conditional=None, page_size=20):
    cond = TABLE_COND + (conditional or [])
    return dash_table.DataTable(
        data=data_records,
        columns=columns,
        sort_action="native",
        sort_mode="multi",
        page_size=page_size,
        style_header=TABLE_HEADER,
        style_data=TABLE_DATA,
        style_data_conditional=cond,
        style_table={"overflowX": "auto"},
        style_cell={"textAlign": "right", "padding": "8px 12px",
                     "fontSize": "0.9em"},
        style_cell_conditional=[
            {"if": {"column_id": c}, "textAlign": "left"}
            for c in ("type", "event", "game", "role", "strike",
                       "direction", "status", "bucket")
        ],
        id=f"table-{id_prefix}",
    )


def pnl_color(val):
    """Return green/red/neutral color for a P&L value."""
    if val > 0:
        return COLORS["green"]
    elif val < 0:
        return COLORS["red"]
    return COLORS["text_muted"]


def _truncate(text, max_len=18):
    """Truncate long chart labels."""
    return text[:max_len] + "\u2026" if len(text) > max_len else text


PNL_COND = [
    {"if": {"filter_query": "{net} > 0", "column_id": "net"},
     "color": COLORS["green"], "fontWeight": "bold"},
    {"if": {"filter_query": "{net} < 0", "column_id": "net"},
     "color": COLORS["red"], "fontWeight": "bold"},
    {"if": {"filter_query": "{roi} > 0", "column_id": "roi"},
     "color": COLORS["green"]},
    {"if": {"filter_query": "{roi} < 0", "column_id": "roi"},
     "color": COLORS["red"]},
    {"if": {"filter_query": "{pnl} > 0", "column_id": "pnl"},
     "color": COLORS["green"]},
    {"if": {"filter_query": "{pnl} < 0", "column_id": "pnl"},
     "color": COLORS["red"]},
]


# ── Tab renderers ────────────────────────────────────────────────────
def render_pnl():
    d = DATA
    results = d.get("market_results", {})
    settled_n = sum(1 for r in results.values() if r is not None)
    total_pnl = sum(r["net"] for r in d.get("pnl_by_type", []))
    total_cost = sum(r["cost"] for r in d.get("pnl_by_type", []))
    total_roi = pct(total_pnl, total_cost) if total_cost else 0
    open_cost = sum(r["cost"] for r in d.get("concentration", []))

    cards = html.Div(style={"display": "flex", "flexWrap": "wrap", "gap": "12px",
                             "marginBottom": "16px"}, children=[
        stat_card("Account", fmt_dollars(
            d.get("balance_cash", 0) + d.get("balance_portfolio", 0))),
        stat_card("Settled P&L", fmt_dollars(total_pnl), pnl_color(total_pnl)),
        stat_card("ROI", f"{total_roi:+.1f}%", pnl_color(total_roi)),
        stat_card("Settled Mkts", str(settled_n)),
        stat_card("Open Exposure", fmt_dollars(open_cost)),
    ])

    type_table = make_table(
        [{**r, "net": round(r["net"], 2), "roi": round(r["roi"], 1),
          "pnl": round(r["pnl"], 2), "cost": round(r["cost"], 2),
          "fees": round(r["fees"], 2)}
         for r in d.get("pnl_by_type", [])],
        [{"name": "Type", "id": "type"},
         {"name": "Fills", "id": "fills", "type": "numeric"},
         {"name": "Win", "id": "wins", "type": "numeric"},
         {"name": "Loss", "id": "losses", "type": "numeric"},
         {"name": "Cost", "id": "cost", "type": "numeric",
          "format": FMT_2F},
         {"name": "P&L", "id": "pnl", "type": "numeric",
          "format": FMT_2F},
         {"name": "Fees", "id": "fees", "type": "numeric",
          "format": FMT_2F},
         {"name": "Net", "id": "net", "type": "numeric",
          "format": FMT_2F},
         {"name": "ROI %", "id": "roi", "type": "numeric",
          "format": FMT_1F}],
        "pnl-type", PNL_COND,
    )

    # Bar chart by event — truncate labels, add hover for full name
    events = d.get("pnl_by_event", [])[:20]
    bar_fig = go.Figure()
    bar_fig.add_bar(
        x=[_truncate(e["event"]) for e in events],
        y=[e["net"] for e in events],
        marker_color=[COLORS["green"] if e["net"] >= 0 else COLORS["red"]
                      for e in events],
        text=[f"${e['net']:+,.0f}" for e in events],
        textposition="outside",
        hovertext=[f"{e['event']}<br>Net: ${e['net']:+,.2f}"
                   f"<br>ROI: {e['roi']:+.1f}%" for e in events],
        hoverinfo="text",
    )
    bar_fig.update_layout(**GRAPH_LAYOUT, height=420, xaxis_tickangle=-45,
                          yaxis_title="Net P&L ($)",
                          title=dict(text="P&L by Game Event (top 20)", x=0.5))

    event_table = make_table(
        [{**r, "net": round(r["net"], 2), "roi": round(r["roi"], 1),
          "pnl": round(r["pnl"], 2), "cost": round(r["cost"], 2)}
         for r in d.get("pnl_by_event", [])],
        [{"name": "Event", "id": "event"},
         {"name": "Fills", "id": "fills", "type": "numeric"},
         {"name": "Cost", "id": "cost", "type": "numeric",
          "format": FMT_2F},
         {"name": "Net", "id": "net", "type": "numeric",
          "format": FMT_2F},
         {"name": "ROI %", "id": "roi", "type": "numeric",
          "format": FMT_1F}],
        "pnl-event", PNL_COND,
    )

    return html.Div([
        cards,
        html.Div(style=CARD_STYLE, children=[
            html.H3("By Market Type",
                     style={"color": COLORS["text"], "marginTop": 0}),
            type_table]),
        html.Div(style=CARD_STYLE, children=[dcc.Graph(figure=bar_fig)]),
        html.Div(style=CARD_STYLE, children=[
            html.H3("By Game Event",
                     style={"color": COLORS["text"], "marginTop": 0}),
            event_table]),
    ])


def render_maker_taker():
    d = DATA
    ms = d.get("maker_stats", {})
    ts = d.get("taker_stats", {})

    cards = html.Div(style={"display": "flex", "flexWrap": "wrap", "gap": "12px",
                             "marginBottom": "16px"}, children=[
        stat_card("Maker Fills", f"{ms.get('fills', 0):,}"),
        stat_card("Maker ROI", f"{ms.get('roi', 0):+.1f}%",
                  pnl_color(ms.get("roi", 0))),
        stat_card("Maker Win%", f"{ms.get('win_pct', 0):.1f}%"),
        stat_card("Taker Fills", f"{ts.get('fills', 0):,}"),
        stat_card("Taker ROI", f"{ts.get('roi', 0):+.1f}%",
                  pnl_color(ts.get("roi", 0))),
        stat_card("Taker Win%", f"{ts.get('win_pct', 0):.1f}%"),
    ])

    # ROI bar chart by role + type
    mt_data = d.get("maker_taker_by_type", [])
    bar_fig = go.Figure()
    for role in ("Maker", "Taker"):
        rows = [r for r in mt_data if r["role"] == role]
        bar_fig.add_bar(
            name=role,
            x=[r["type"] for r in rows],
            y=[r["roi"] for r in rows],
            text=[f"{r['roi']:+.1f}%" for r in rows],
            textposition="outside",
            marker_color=(COLORS["accent"] if role == "Maker"
                          else COLORS["accent2"]),
        )
    bar_fig.update_layout(**GRAPH_LAYOUT, barmode="group", height=350,
                          yaxis_title="ROI %",
                          title=dict(text="ROI by Role & Market Type", x=0.5))

    mt_table = make_table(
        [{**r, "net": round(r["net"], 2), "roi": round(r["roi"], 1),
          "cost": round(r["cost"], 2)}
         for r in mt_data],
        [{"name": "Role", "id": "role"},
         {"name": "Type", "id": "type"},
         {"name": "Fills", "id": "fills", "type": "numeric"},
         {"name": "Cost", "id": "cost", "type": "numeric",
          "format": FMT_2F},
         {"name": "Net", "id": "net", "type": "numeric",
          "format": FMT_2F},
         {"name": "ROI %", "id": "roi", "type": "numeric",
          "format": FMT_1F}],
        "mt-detail", PNL_COND,
    )

    # Pie chart
    pie_fig = go.Figure(data=[go.Pie(
        labels=["Maker", "Taker"],
        values=[ms.get("cost", 0), ts.get("cost", 0)],
        marker=dict(colors=[COLORS["accent"], COLORS["accent2"]]),
        textinfo="label+percent", hole=0.4,
    )])
    pie_fig.update_layout(**GRAPH_LAYOUT, height=300,
                          title=dict(text="Cost Allocation", x=0.5))

    return html.Div([
        cards,
        html.Div(style={"display": "flex", "flexWrap": "wrap", "gap": "16px"},
                 children=[
            html.Div(style={**CARD_STYLE, "flex": "1", "minWidth": "300px"},
                     children=[dcc.Graph(figure=bar_fig)]),
            html.Div(style={**CARD_STYLE, "flex": "0 0 300px"},
                     children=[dcc.Graph(figure=pie_fig)]),
        ]),
        html.Div(style=CARD_STYLE, children=[
            html.H3("Detail",
                     style={"color": COLORS["text"], "marginTop": 0}),
            mt_table]),
    ])


def render_fill_patterns():
    d = DATA
    ms = d.get("maker_fill_stats", {})
    ts = d.get("taker_fill_stats", {})
    fr = d.get("fill_rate_stats", {})

    # ── Section 1: Order size stats ──
    size_cards = html.Div(style={"display": "flex", "flexWrap": "wrap",
                                  "gap": "12px", "marginBottom": "16px"},
                          children=[
        stat_card("Maker Orders", f"{ms.get('orders', 0):,}"),
        stat_card("Maker Avg Size", f"{ms.get('avg_size', 0):.0f} cts"),
        stat_card("Maker Single-Fill", f"{ms.get('single_pct', 0):.0f}%"),
        stat_card("Taker Orders", f"{ts.get('orders', 0):,}"),
        stat_card("Taker Avg Size", f"{ts.get('avg_size', 0):.0f} cts"),
        stat_card("Adverse Flags", str(len(d.get("adverse_orders", []))),
                  COLORS["yellow"] if d.get("adverse_orders")
                  else COLORS["green"]),
    ])

    # Order size histogram
    hist_fig = go.Figure()
    if ms.get("order_sizes"):
        hist_fig.add_histogram(x=ms["order_sizes"], name="Maker",
                               marker_color=COLORS["accent"], opacity=0.7,
                               xbins=dict(size=5))
    if ts.get("order_sizes"):
        hist_fig.add_histogram(x=ts["order_sizes"], name="Taker",
                               marker_color=COLORS["accent2"], opacity=0.7,
                               xbins=dict(size=5))
    hist_fig.update_layout(**GRAPH_LAYOUT, barmode="overlay", height=350,
                           xaxis_title="Contracts per Order",
                           yaxis_title="Count",
                           title=dict(text="Order Size Distribution", x=0.5))

    # Adverse selection table
    adverse = d.get("adverse_orders", [])[:20]
    adverse_table = make_table(
        [{"game": r["game"], "type": r["type"], "strike": r["strike"],
          "contracts": r["contracts"], "status": r["status"],
          "pnl": round(r["pnl"], 2) if r["pnl"] is not None else None}
         for r in adverse],
        [{"name": "Game", "id": "game"},
         {"name": "Type", "id": "type"},
         {"name": "Strike", "id": "strike"},
         {"name": "Contracts", "id": "contracts", "type": "numeric"},
         {"name": "Status", "id": "status"},
         {"name": "P&L", "id": "pnl", "type": "numeric",
          "format": FMT_2F}],
        "adverse",
        [{"if": {"filter_query": "{pnl} < 0", "column_id": "pnl"},
          "color": COLORS["red"]},
         {"if": {"filter_query": "{pnl} > 0", "column_id": "pnl"},
          "color": COLORS["green"]},
         {"if": {"filter_query": '{status} = "Open"', "column_id": "status"},
          "color": COLORS["yellow"]}],
    )

    # ── Section 2: Maker Fill Rate ──
    bs = d.get("bite_stats", {})
    fill_rate_cards = html.Div(
        style={"display": "flex", "flexWrap": "wrap", "gap": "12px",
               "marginBottom": "16px"},
        children=[
            stat_card("Avg Bite Size", f"{bs.get('avg', 0):.1f}%"),
            stat_card("Median Bite Size", f"{bs.get('median', 0):.1f}%"),
            stat_card("Order Fill Rate", f"{fr.get('avg', 0):.1f}%"),
            stat_card("Fully Filled", f"{fr.get('fully_filled', 0):.1f}%"),
            stat_card("Partial Fill", f"{fr.get('partial', 0):.1f}%",
                      COLORS["yellow"] if fr.get("partial", 0) > 30
                      else COLORS["green"]),
        ])

    # Bite size histogram (per-fill % of order)
    bite_fig = go.Figure()
    bite_pcts = d.get("bite_pcts", [])
    if bite_pcts:
        bite_fig.add_histogram(
            x=bite_pcts, marker_color=COLORS["accent"],
            xbins=dict(start=0, end=105, size=5),
        )
    bite_fig.update_layout(
        **GRAPH_LAYOUT, height=300,
        xaxis_title="Fill Size as % of Resting Order",
        yaxis_title="Fills",
        title=dict(text="Maker Bite Size Distribution", x=0.5))

    # Order fill rate histogram
    fill_rate_fig = go.Figure()
    fill_rates = [r["fill_pct"] for r in d.get("fill_rate_rows", [])]
    if fill_rates:
        fill_rate_fig.add_histogram(
            x=fill_rates, marker_color=COLORS["accent2"],
            xbins=dict(start=0, end=105, size=10),
        )
    fill_rate_fig.update_layout(
        **GRAPH_LAYOUT, height=300,
        xaxis_title="Total Filled % of Order Size",
        yaxis_title="Orders",
        title=dict(text="Maker Order Fill Rate Distribution", x=0.5))

    # Fill rate by market type table (includes bite size)
    fr_by_type = d.get("fill_rate_by_type", [])
    bite_by_type = {r["type"]: r for r in d.get("bite_by_type", [])}
    fr_type_table = make_table(
        [{"type": r["type"], "orders": r["count"],
          "avg_fill": round(r["avg"], 1),
          "fully_filled": round(r["fully_filled"], 1),
          "avg_bite": round(bite_by_type.get(r["type"], {}).get("avg", 0), 1),
          "bite_fills": bite_by_type.get(r["type"], {}).get("count", 0)}
         for r in fr_by_type],
        [{"name": "Type", "id": "type"},
         {"name": "Orders", "id": "orders", "type": "numeric"},
         {"name": "Avg Fill %", "id": "avg_fill", "type": "numeric",
          "format": FMT_1F},
         {"name": "Fully Filled %", "id": "fully_filled", "type": "numeric",
          "format": FMT_1F},
         {"name": "Avg Bite %", "id": "avg_bite", "type": "numeric",
          "format": FMT_1F},
         {"name": "Fills", "id": "bite_fills", "type": "numeric"}],
        "fr-type",
    )

    # ── Section 3: Fill Timing ──
    buckets = d.get("timing_buckets", {})
    labels = d.get("timing_bucket_labels", [])
    timing_hours = d.get("timing_hours", [])

    total_timing_fills = sum(b["fills"] for b in buckets.values())
    last_hour = sum(buckets.get(l, {}).get("fills", 0)
                    for l in ["<15m", "15-30m", "30m-1h"])
    avg_hours = (sum(timing_hours) / len(timing_hours)) if timing_hours else 0

    timing_cards = html.Div(
        style={"display": "flex", "flexWrap": "wrap", "gap": "12px",
               "marginBottom": "16px"},
        children=[
            stat_card("Avg Time Before Close",
                      f"{avg_hours:.1f}h" if avg_hours < 24
                      else f"{avg_hours:.0f}h"),
            stat_card("Fills < 1h Before Close",
                      f"{pct(last_hour, total_timing_fills):.1f}%"
                      if total_timing_fills else "N/A",
                      COLORS["yellow"]
                      if total_timing_fills and
                      pct(last_hour, total_timing_fills) > 30
                      else COLORS["green"]),
            stat_card("Total Timed Fills", f"{total_timing_fills:,}"),
        ])

    # Timing bar chart (fills by bucket)
    bucket_fills = [buckets.get(l, {}).get("fills", 0) for l in labels]
    bucket_maker = [buckets.get(l, {}).get("maker_fills", 0) for l in labels]
    bucket_taker = [f - m for f, m in zip(bucket_fills, bucket_maker)]

    timing_fig = go.Figure()
    timing_fig.add_bar(name="Maker", x=labels, y=bucket_maker,
                       marker_color=COLORS["accent"])
    timing_fig.add_bar(name="Taker", x=labels, y=bucket_taker,
                       marker_color=COLORS["accent2"])
    timing_fig.update_layout(**GRAPH_LAYOUT, barmode="stack", height=350,
                             xaxis_title="Time Before Market Close",
                             yaxis_title="Fills",
                             title=dict(text="Fill Timing Distribution",
                                        x=0.5))

    # P&L by timing bucket
    timing_pnl_fig = go.Figure()
    bucket_pnl = [buckets.get(l, {}).get("pnl", 0) for l in labels]
    timing_pnl_fig.add_bar(
        x=labels, y=bucket_pnl,
        marker_color=[COLORS["green"] if p >= 0 else COLORS["red"]
                      for p in bucket_pnl],
        text=[f"${p:+,.0f}" for p in bucket_pnl],
        textposition="outside",
    )
    timing_pnl_fig.update_layout(
        **GRAPH_LAYOUT, height=350,
        xaxis_title="Time Before Market Close",
        yaxis_title="Settled P&L ($)",
        title=dict(text="P&L by Fill Timing (settled only)", x=0.5))

    # Timing detail table — add sort key for chronological order
    timing_table_data = []
    for i, l in enumerate(labels):
        b = buckets.get(l, {})
        if b.get("fills", 0) == 0:
            continue
        avg_pnl = (b["pnl"] / b["settled_fills"]
                   if b["settled_fills"] else None)
        timing_table_data.append({
            "sort_key": i,
            "bucket": l,
            "fills": b["fills"],
            "contracts": b["contracts"],
            "maker_pct": round(pct(b["maker_fills"], b["fills"]), 1),
            "pnl": round(b["pnl"], 2) if b["settled_fills"] else None,
            "avg_pnl": round(avg_pnl, 2) if avg_pnl is not None else None,
        })

    timing_table = make_table(
        timing_table_data,
        [{"name": "#", "id": "sort_key", "type": "numeric"},
         {"name": "Window", "id": "bucket"},
         {"name": "Fills", "id": "fills", "type": "numeric"},
         {"name": "Contracts", "id": "contracts", "type": "numeric"},
         {"name": "Maker %", "id": "maker_pct", "type": "numeric",
          "format": FMT_1F},
         {"name": "Net P&L", "id": "pnl", "type": "numeric",
          "format": FMT_2F},
         {"name": "Avg P&L/Fill", "id": "avg_pnl", "type": "numeric",
          "format": FMT_2F}],
        "timing",
        [{"if": {"filter_query": "{pnl} < 0", "column_id": "pnl"},
          "color": COLORS["red"]},
         {"if": {"filter_query": "{pnl} > 0", "column_id": "pnl"},
          "color": COLORS["green"]},
         {"if": {"filter_query": "{avg_pnl} < 0", "column_id": "avg_pnl"},
          "color": COLORS["red"]},
         {"if": {"filter_query": "{avg_pnl} > 0", "column_id": "avg_pnl"},
          "color": COLORS["green"]}],
    )

    # ── Assemble ──
    return html.Div([
        # Order sizes
        size_cards,
        html.Div(style=CARD_STYLE, children=[dcc.Graph(figure=hist_fig)]),
        html.Div(style=CARD_STYLE, children=[
            html.H3(f"Adverse Selection \u2014 Maker orders "
                     f">{d.get('adverse_p75', 0)} cts, single fill",
                     style={"color": COLORS["text"], "marginTop": 0,
                            "fontSize": "1em"}),
            adverse_table]),

        # Fill rate
        html.H2("Maker Fill Rate",
                 style={"color": COLORS["accent"], "marginTop": "28px",
                        "marginBottom": "12px", "fontSize": "1.2em",
                        "borderBottom": f"1px solid {COLORS['card_border']}",
                        "paddingBottom": "8px"}),
        fill_rate_cards,
        html.Div(style={"display": "flex", "flexWrap": "wrap",
                         "gap": "16px"}, children=[
            html.Div(style={**CARD_STYLE, "flex": "1", "minWidth": "300px"},
                     children=[dcc.Graph(figure=bite_fig)]),
            html.Div(style={**CARD_STYLE, "flex": "1", "minWidth": "300px"},
                     children=[dcc.Graph(figure=fill_rate_fig)]),
        ]),
        html.Div(style=CARD_STYLE, children=[
            html.H3("By Market Type",
                     style={"color": COLORS["text"], "marginTop": 0}),
            fr_type_table]),

        # Timing
        html.H2("Fill Timing",
                 style={"color": COLORS["accent"], "marginTop": "28px",
                        "marginBottom": "12px", "fontSize": "1.2em",
                        "borderBottom": f"1px solid {COLORS['card_border']}",
                        "paddingBottom": "8px"}),
        timing_cards,
        html.Div(style={"display": "flex", "flexWrap": "wrap",
                         "gap": "16px"}, children=[
            html.Div(style={**CARD_STYLE, "flex": "1", "minWidth": "300px"},
                     children=[dcc.Graph(figure=timing_fig)]),
            html.Div(style={**CARD_STYLE, "flex": "1", "minWidth": "300px"},
                     children=[dcc.Graph(figure=timing_pnl_fig)]),
        ]),
        html.Div(style=CARD_STYLE, children=[
            html.H3("Timing Detail",
                     style={"color": COLORS["text"], "marginTop": 0}),
            timing_table]),
    ])


def render_concentration():
    d = DATA
    conc = d.get("concentration", [])
    total_exp = sum(r["exposure"] for r in conc)
    top1 = conc[0]["exposure"] if conc else 0
    top3 = (sum(r["exposure"] for r in conc[:3])
            if len(conc) >= 3 else total_exp)

    cards = html.Div(style={"display": "flex", "flexWrap": "wrap", "gap": "12px",
                             "marginBottom": "16px"}, children=[
        stat_card("Total Exposure", fmt_dollars(total_exp)),
        stat_card("Largest Event", fmt_dollars(top1)),
        stat_card("Top 1 %", f"{pct(top1, total_exp):.1f}%",
                  COLORS["yellow"] if total_exp and pct(top1, total_exp) > 20
                  else COLORS["green"]),
        stat_card("Top 3 %", f"{pct(top3, total_exp):.1f}%",
                  COLORS["yellow"] if total_exp and pct(top3, total_exp) > 50
                  else COLORS["green"]),
        stat_card("Events", str(sum(1 for r in conc if r["exposure"] > 0))),
    ])

    # Pie by market type — guard against empty
    by_type = defaultdict(float)
    for r in conc:
        by_type[r["type"]] += r["exposure"]

    if by_type and any(v > 0 for v in by_type.values()):
        pie_fig = go.Figure(data=[go.Pie(
            labels=list(by_type.keys()),
            values=list(by_type.values()),
            marker=dict(colors=PIE_COLORS[:len(by_type)]),
            textinfo="label+percent", hole=0.4,
        )])
    else:
        pie_fig = go.Figure()
        pie_fig.add_annotation(text="No positions", showarrow=False,
                               font=dict(size=16, color=COLORS["text_muted"]))
    pie_fig.update_layout(**GRAPH_LAYOUT, height=300,
                          title=dict(text="Exposure by Market Type", x=0.5))

    conc_table = make_table(
        [{"game": r["game"], "type": r["type"],
          "exposure": round(r["exposure"], 2), "cost": round(r["cost"], 2),
          "pct": round(pct(r["exposure"], total_exp), 1) if total_exp else 0}
         for r in conc if r["exposure"] > 0],
        [{"name": "Game", "id": "game"},
         {"name": "Type", "id": "type"},
         {"name": "Exposure", "id": "exposure", "type": "numeric",
          "format": FMT_2F},
         {"name": "Cost", "id": "cost", "type": "numeric",
          "format": FMT_2F},
         {"name": "% of Total", "id": "pct", "type": "numeric",
          "format": FMT_1F}],
        "conc",
    )

    # Correlation table
    corr = d.get("correlation", [])
    corr_table = make_table(
        [{"game": r["game"], "spread": round(r["spread"]),
          "ml": round(r["ml"]), "direction": r["direction"]}
         for r in corr],
        [{"name": "Game", "id": "game"},
         {"name": "Spread Pos", "id": "spread", "type": "numeric"},
         {"name": "ML Pos", "id": "ml", "type": "numeric"},
         {"name": "Direction", "id": "direction"}],
        "corr",
        [{"if": {"filter_query": '{direction} = "SAME"',
                 "column_id": "direction"},
          "color": COLORS["yellow"], "fontWeight": "bold"}],
    )

    return html.Div([
        cards,
        html.Div(style={"display": "flex", "flexWrap": "wrap", "gap": "16px"},
                 children=[
            html.Div(style={**CARD_STYLE, "flex": "1", "minWidth": "300px"},
                     children=[
                html.H3("Positions",
                         style={"color": COLORS["text"], "marginTop": 0}),
                conc_table]),
            html.Div(style={**CARD_STYLE, "flex": "0 0 320px"}, children=[
                dcc.Graph(figure=pie_fig)]),
        ]),
        html.Div(style=CARD_STYLE, children=[
            html.H3("Directional Correlation (Spread + ML)",
                     style={"color": COLORS["text"], "marginTop": 0}),
            html.P("Open positions only \u2014 settled games are zeroed out "
                   "by Kalshi.",
                   style={"color": COLORS["text_muted"], "fontSize": "0.8em",
                           "marginTop": "-8px", "marginBottom": "12px"}),
            corr_table]),
    ])


def render_event_detail():
    d = DATA
    events_map = d.get("events_map", {})
    options = sorted(events_map.keys())

    return html.Div([
        html.Div(style=CARD_STYLE, children=[
            html.Label("Select Event:",
                        style={"color": COLORS["text_muted"],
                               "marginRight": "12px"}),
            dcc.Dropdown(
                id="event-dropdown",
                options=[{"label": e, "value": e} for e in options],
                value=options[0] if options else None,
                style={"minWidth": "300px"},
            ),
        ]),
        html.Div(id="event-detail-content"),
    ])


# ── App setup ────────────────────────────────────────────────────────
app = dash.Dash(__name__, title="MM Performance",
                suppress_callback_exceptions=True)

app.index_string = (
    '<!DOCTYPE html><html><head>'
    '{%metas%}<title>{%title%}</title>{%favicon%}{%css%}'
    '<style>' + DROPDOWN_CSS + '</style>'
    '</head><body>{%app_entry%}<footer>'
    '{%config%}{%scripts%}{%renderer%}'
    '</footer></body></html>'
)

TAB_STYLE = {"color": COLORS["text_muted"],
             "backgroundColor": COLORS["card"],
             "border": "none", "padding": "12px 20px"}
TAB_SELECTED = {"color": COLORS["accent"],
                "backgroundColor": COLORS["bg"],
                "borderTop": f"2px solid {COLORS['accent']}",
                "padding": "12px 20px"}

app.layout = html.Div(
    style={"backgroundColor": COLORS["bg"], "minHeight": "100vh",
           "padding": "16px",
           "fontFamily": ("-apple-system, BlinkMacSystemFont, "
                          "'Segoe UI', sans-serif"),
           "color": COLORS["text"]},
    children=[
        # Header
        html.Div(style={"background": COLORS["header_bg"],
                         "padding": "20px 28px",
                         "borderRadius": "12px", "marginBottom": "16px",
                         "display": "flex", "justifyContent": "space-between",
                         "alignItems": "center", "flexWrap": "wrap"},
                 children=[
            html.Div([
                html.H1("CBB 1H Market Maker",
                         style={"margin": 0, "color": "white",
                                "fontSize": "1.6em"}),
                html.P(id="last-updated",
                       style={"margin": "4px 0 0 0",
                              "color": "rgba(255,255,255,0.85)"}),
            ]),
            html.Button("Refresh", id="refresh-btn", n_clicks=0, style={
                "background": "rgba(255,255,255,0.2)",
                "border": "1px solid rgba(255,255,255,0.3)",
                "color": "white", "padding": "10px 20px",
                "borderRadius": "8px", "fontSize": "0.95em",
                "fontWeight": "600", "cursor": "pointer"}),
        ]),
        # Tabs
        dcc.Tabs(id="main-tabs", value="pnl", colors={
            "border": COLORS["card_border"],
            "primary": COLORS["accent"],
            "background": COLORS["card"]}, children=[
            dcc.Tab(label="P&L", value="pnl",
                    style=TAB_STYLE, selected_style=TAB_SELECTED),
            dcc.Tab(label="Maker / Taker", value="mt",
                    style=TAB_STYLE, selected_style=TAB_SELECTED),
            dcc.Tab(label="Fill Patterns", value="fills",
                    style=TAB_STYLE, selected_style=TAB_SELECTED),
            dcc.Tab(label="Concentration", value="conc",
                    style=TAB_STYLE, selected_style=TAB_SELECTED),
            dcc.Tab(label="Event Detail", value="detail",
                    style=TAB_STYLE, selected_style=TAB_SELECTED),
        ]),
        # Loading spinner wraps tab content (visible during load + refresh)
        dcc.Loading(id="loading-tabs", type="circle",
                    color=COLORS["accent"],
                    children=[html.Div(id="tab-content")]),
    ],
)


# Single callback handles both tab switches AND refresh clicks
@callback(Output("tab-content", "children"),
          Output("last-updated", "children"),
          Input("main-tabs", "value"),
          Input("refresh-btn", "n_clicks"))
def render_tab(tab, n_clicks):
    triggered = ctx.triggered_id

    if triggered == "refresh-btn":
        load_data()
    elif not DATA.get("loaded"):
        load_data()

    # Handle total failure (no data at all)
    if DATA.get("error") and not DATA.get("loaded"):
        return (html.Div(style=CARD_STYLE, children=[
            html.H3("Error Loading Data",
                     style={"color": COLORS["red"]}),
            html.P(DATA["error"], style={"color": COLORS["text"]}),
            html.P("Check API credentials and try Refresh.",
                   style={"color": COLORS["text_muted"]}),
        ]), f"Error: {DATA.get('ts', '?')}")

    ts_text = f"Updated: {DATA.get('ts', '?')}"
    if DATA.get("error"):
        ts_text += f"  (refresh failed: {DATA['error']})"

    renderers = {
        "pnl": render_pnl,
        "mt": render_maker_taker,
        "fills": render_fill_patterns,
        "conc": render_concentration,
        "detail": render_event_detail,
    }
    content = renderers.get(tab, lambda: html.Div("Unknown tab"))()
    return content, ts_text


# Event detail — no prevent_initial_call so first event renders immediately
@callback(Output("event-detail-content", "children"),
          Input("event-dropdown", "value"))
def render_event(event_name):
    if not event_name:
        return html.Div()
    fills = DATA.get("events_map", {}).get(event_name, [])
    results = DATA.get("market_results", {})
    rows = []
    for f in sorted(fills, key=lambda x: x["created_time"]):
        t = f["ticker"]
        _, _, strike = parse_ticker(t)
        cnt = contracts(f)
        cost = fill_cost_cents(f) * cnt / 100
        settled = t in results and results[t] is not None
        pnl = (fill_pnl_cents(f, results.get(t, "")) * cnt / 100
               if settled else None)
        rows.append({
            "time": f["created_time"][:19].replace("T", " "),
            "strike": strike, "side": f["side"], "contracts": cnt,
            "cost": round(cost, 2),
            "taker": "Yes" if f.get("is_taker") else "",
            "result": results.get(t, "") or "open",
            "pnl": round(pnl, 2) if pnl is not None else None,
        })
    total_cost = sum(r["cost"] for r in rows)
    total_pnl = sum(r["pnl"] for r in rows
                    if isinstance(r["pnl"], (int, float)))
    has_settled = any(isinstance(r["pnl"], (int, float)) for r in rows)

    pnl_display = fmt_dollars(total_pnl) if has_settled else "open"
    pnl_card_color = pnl_color(total_pnl) if has_settled else COLORS["text_muted"]

    return html.Div(style=CARD_STYLE, children=[
        html.Div(style={"display": "flex", "gap": "12px",
                         "marginBottom": "12px"}, children=[
            stat_card("Fills", str(len(rows))),
            stat_card("Cost", fmt_dollars(total_cost)),
            stat_card("P&L", pnl_display, pnl_card_color),
        ]),
        make_table(rows,
            [{"name": "Time", "id": "time"},
             {"name": "Strike", "id": "strike"},
             {"name": "Side", "id": "side"},
             {"name": "Cts", "id": "contracts", "type": "numeric"},
             {"name": "Cost", "id": "cost", "type": "numeric",
              "format": FMT_2F},
             {"name": "Taker", "id": "taker"},
             {"name": "Result", "id": "result"},
             {"name": "P&L", "id": "pnl", "type": "numeric",
              "format": FMT_2F}],
            "event-fills",
            [{"if": {"filter_query": "{pnl} < 0", "column_id": "pnl"},
              "color": COLORS["red"]},
             {"if": {"filter_query": "{pnl} > 0", "column_id": "pnl"},
              "color": COLORS["green"]},
             {"if": {"filter_query": '{result} = "yes"',
                     "column_id": "result"},
              "color": COLORS["green"]},
             {"if": {"filter_query": '{result} = "no"',
                     "column_id": "result"},
              "color": COLORS["red"]}],
            page_size=50),
    ])


if __name__ == "__main__":
    print("Loading data from Kalshi API...")
    load_data()
    print(f"Loaded {len(DATA.get('cbb_fills', []))} fills. "
          f"Starting dashboard...")
    app.run(debug=False, host="0.0.0.0", port=8084, threaded=True)
