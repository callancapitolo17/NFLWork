#!/usr/bin/env python3
"""
Kalshi CBB 1H Market Maker — Performance Dashboard

Dash web app visualizing the same 5 analyses as analyze_performance.py.
Accessible from phone on same Wi-Fi at http://<mac-ip>:8084

Usage:
    python3 dashboard.py
"""

import dash
from dash import dcc, html, dash_table, callback, Input, Output
from dash.dash_table.Format import Format, Scheme
import plotly.graph_objects as go
from collections import defaultdict
from datetime import datetime

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


# ── Data store (refreshed on demand) ────────────────────────────────
DATA = {"loaded": False}


def load_data():
    DATA.update(pull_data())
    DATA["loaded"] = True
    DATA["ts"] = datetime.now().strftime("%H:%M:%S")
    # Pre-compute aggregates
    _compute_all(DATA)


def _compute_all(data):
    """Pre-compute DataFrames for all tabs."""
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
                settled = fl[0]["ticker"] in settled_tickers
                pnl = fill_pnl_cents(fl[0], results.get(fl[0]["ticker"], "")) * ct / 100 if settled else None
                adverse.append({"game": parse_game_teams(game), "type": mtype,
                                "strike": strike, "contracts": ct,
                                "pnl": pnl, "settled": settled})
        data["adverse_orders"] = sorted(adverse, key=lambda x: x["contracts"], reverse=True)
        data["adverse_p75"] = p75
    else:
        data["adverse_orders"] = []
        data["adverse_p75"] = 0

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
    data["concentration"] = sorted(conc_rows, key=lambda x: x["exposure"], reverse=True)

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
            for c in ("type", "event", "game", "role", "strike", "direction")
        ],
        id=f"table-{id_prefix}",
    )


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
        stat_card("Account", fmt_dollars(d.get("balance_cash", 0) + d.get("balance_portfolio", 0))),
        stat_card("Settled P&L", fmt_dollars(total_pnl),
                  COLORS["green"] if total_pnl > 0 else COLORS["red"]),
        stat_card("ROI", f"{total_roi:+.1f}%",
                  COLORS["green"] if total_roi > 0 else COLORS["red"]),
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
          "format": Format(precision=2, scheme=Scheme.fixed)},
         {"name": "P&L", "id": "pnl", "type": "numeric",
          "format": Format(precision=2, scheme=Scheme.fixed)},
         {"name": "Fees", "id": "fees", "type": "numeric",
          "format": Format(precision=2, scheme=Scheme.fixed)},
         {"name": "Net", "id": "net", "type": "numeric",
          "format": Format(precision=2, scheme=Scheme.fixed)},
         {"name": "ROI %", "id": "roi", "type": "numeric",
          "format": Format(precision=1, scheme=Scheme.fixed)}],
        "pnl-type", PNL_COND,
    )

    # Bar chart by event
    events = d.get("pnl_by_event", [])[:30]
    bar_fig = go.Figure()
    bar_fig.add_bar(
        x=[e["event"] for e in events],
        y=[e["net"] for e in events],
        marker_color=[COLORS["green"] if e["net"] >= 0 else COLORS["red"] for e in events],
        text=[f"${e['net']:+,.0f}" for e in events],
        textposition="outside",
    )
    bar_fig.update_layout(**GRAPH_LAYOUT, height=400, xaxis_tickangle=-45,
                          yaxis_title="Net P&L ($)",
                          title=dict(text="P&L by Game Event", x=0.5))

    event_table = make_table(
        [{**r, "net": round(r["net"], 2), "roi": round(r["roi"], 1),
          "pnl": round(r["pnl"], 2), "cost": round(r["cost"], 2)}
         for r in d.get("pnl_by_event", [])],
        [{"name": "Event", "id": "event"},
         {"name": "Fills", "id": "fills", "type": "numeric"},
         {"name": "Cost", "id": "cost", "type": "numeric",
          "format": Format(precision=2, scheme=Scheme.fixed)},
         {"name": "Net", "id": "net", "type": "numeric",
          "format": Format(precision=2, scheme=Scheme.fixed)},
         {"name": "ROI %", "id": "roi", "type": "numeric",
          "format": Format(precision=1, scheme=Scheme.fixed)}],
        "pnl-event", PNL_COND,
    )

    return html.Div([
        cards,
        html.Div(style=CARD_STYLE, children=[
            html.H3("By Market Type", style={"color": COLORS["text"], "marginTop": 0}),
            type_table]),
        html.Div(style=CARD_STYLE, children=[dcc.Graph(figure=bar_fig)]),
        html.Div(style=CARD_STYLE, children=[
            html.H3("By Game Event", style={"color": COLORS["text"], "marginTop": 0}),
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
                  COLORS["green"] if ms.get("roi", 0) > 0 else COLORS["red"]),
        stat_card("Maker Win%", f"{ms.get('win_pct', 0):.1f}%"),
        stat_card("Taker Fills", f"{ts.get('fills', 0):,}"),
        stat_card("Taker ROI", f"{ts.get('roi', 0):+.1f}%",
                  COLORS["green"] if ts.get("roi", 0) > 0 else COLORS["red"]),
        stat_card("Taker Win%", f"{ts.get('win_pct', 0):.1f}%"),
    ])

    # ROI bar chart by role + type
    mt_data = d.get("maker_taker_by_type", [])
    bar_fig = go.Figure()
    for role in ("Maker", "Taker"):
        rows = [r for r in mt_data if r["role"] == role]
        bar_fig.add_bar(
            name=role, x=[r["type"] for r in rows], y=[r["roi"] for r in rows],
            text=[f"{r['roi']:+.1f}%" for r in rows], textposition="outside",
            marker_color=COLORS["accent"] if role == "Maker" else COLORS["accent2"],
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
          "format": Format(precision=2, scheme=Scheme.fixed)},
         {"name": "Net", "id": "net", "type": "numeric",
          "format": Format(precision=2, scheme=Scheme.fixed)},
         {"name": "ROI %", "id": "roi", "type": "numeric",
          "format": Format(precision=1, scheme=Scheme.fixed)}],
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
        html.Div(style={"display": "flex", "flexWrap": "wrap", "gap": "16px"}, children=[
            html.Div(style={**CARD_STYLE, "flex": "1", "minWidth": "300px"},
                     children=[dcc.Graph(figure=bar_fig)]),
            html.Div(style={**CARD_STYLE, "flex": "0 0 300px"},
                     children=[dcc.Graph(figure=pie_fig)]),
        ]),
        html.Div(style=CARD_STYLE, children=[
            html.H3("Detail", style={"color": COLORS["text"], "marginTop": 0}),
            mt_table]),
    ])


def render_fill_patterns():
    d = DATA
    ms = d.get("maker_fill_stats", {})
    ts = d.get("taker_fill_stats", {})

    cards = html.Div(style={"display": "flex", "flexWrap": "wrap", "gap": "12px",
                             "marginBottom": "16px"}, children=[
        stat_card("Maker Orders", f"{ms.get('orders', 0):,}"),
        stat_card("Maker Avg Size", f"{ms.get('avg_size', 0):.0f} cts"),
        stat_card("Maker Single-Fill", f"{ms.get('single_pct', 0):.0f}%"),
        stat_card("Taker Orders", f"{ts.get('orders', 0):,}"),
        stat_card("Taker Avg Size", f"{ts.get('avg_size', 0):.0f} cts"),
        stat_card("Adverse Flags", str(len(d.get("adverse_orders", []))),
                  COLORS["yellow"] if d.get("adverse_orders") else COLORS["green"]),
    ])

    # Histogram
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
        [{**r, "pnl": round(r["pnl"], 2) if r["pnl"] is not None else "open"}
         for r in adverse],
        [{"name": "Game", "id": "game"},
         {"name": "Type", "id": "type"},
         {"name": "Strike", "id": "strike"},
         {"name": "Contracts", "id": "contracts", "type": "numeric"},
         {"name": "P&L", "id": "pnl"}],
        "adverse",
        [{"if": {"filter_query": '{pnl} contains "-"', "column_id": "pnl"},
          "color": COLORS["red"]},
         {"if": {"filter_query": '{pnl} > 0', "column_id": "pnl"},
          "color": COLORS["green"]}],
    )

    return html.Div([
        cards,
        html.Div(style=CARD_STYLE, children=[dcc.Graph(figure=hist_fig)]),
        html.Div(style=CARD_STYLE, children=[
            html.H3(f"Adverse Selection — Maker orders >{d.get('adverse_p75', 0)} cts, single fill",
                     style={"color": COLORS["text"], "marginTop": 0, "fontSize": "1em"}),
            adverse_table]),
    ])


def render_concentration():
    d = DATA
    conc = d.get("concentration", [])
    total_exp = sum(r["exposure"] for r in conc)
    top1 = conc[0]["exposure"] if conc else 0
    top3 = sum(r["exposure"] for r in conc[:3]) if len(conc) >= 3 else total_exp

    cards = html.Div(style={"display": "flex", "flexWrap": "wrap", "gap": "12px",
                             "marginBottom": "16px"}, children=[
        stat_card("Total Exposure", fmt_dollars(total_exp)),
        stat_card("Largest Event", fmt_dollars(top1)),
        stat_card("Top 1 %", f"{pct(top1, total_exp):.1f}%",
                  COLORS["yellow"] if pct(top1, total_exp) > 20 else COLORS["green"]),
        stat_card("Top 3 %", f"{pct(top3, total_exp):.1f}%",
                  COLORS["yellow"] if pct(top3, total_exp) > 50 else COLORS["green"]),
        stat_card("Events", str(sum(1 for r in conc if r["exposure"] > 0))),
    ])

    # Pie by market type
    by_type = defaultdict(float)
    for r in conc:
        by_type[r["type"]] += r["exposure"]
    pie_fig = go.Figure(data=[go.Pie(
        labels=list(by_type.keys()), values=list(by_type.values()),
        marker=dict(colors=[COLORS["accent"], COLORS["accent2"],
                            COLORS["yellow"], COLORS["text_muted"]]),
        textinfo="label+percent", hole=0.4,
    )])
    pie_fig.update_layout(**GRAPH_LAYOUT, height=300,
                          title=dict(text="Exposure by Market Type", x=0.5))

    conc_table = make_table(
        [{"game": r["game"], "type": r["type"],
          "exposure": round(r["exposure"], 2), "cost": round(r["cost"], 2),
          "pct": round(pct(r["exposure"], total_exp), 1)}
         for r in conc if r["exposure"] > 0],
        [{"name": "Game", "id": "game"},
         {"name": "Type", "id": "type"},
         {"name": "Exposure", "id": "exposure", "type": "numeric",
          "format": Format(precision=2, scheme=Scheme.fixed)},
         {"name": "Cost", "id": "cost", "type": "numeric",
          "format": Format(precision=2, scheme=Scheme.fixed)},
         {"name": "% of Total", "id": "pct", "type": "numeric",
          "format": Format(precision=1, scheme=Scheme.fixed)}],
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
        [{"if": {"filter_query": '{direction} = "SAME"', "column_id": "direction"},
          "color": COLORS["yellow"], "fontWeight": "bold"}],
    )

    return html.Div([
        cards,
        html.Div(style={"display": "flex", "flexWrap": "wrap", "gap": "16px"}, children=[
            html.Div(style={**CARD_STYLE, "flex": "1", "minWidth": "300px"},
                     children=[
                html.H3("Positions", style={"color": COLORS["text"], "marginTop": 0}),
                conc_table]),
            html.Div(style={**CARD_STYLE, "flex": "0 0 320px"}, children=[
                dcc.Graph(figure=pie_fig)]),
        ]),
        html.Div(style=CARD_STYLE, children=[
            html.H3("Directional Correlation (Spread + ML)", style={"color": COLORS["text"], "marginTop": 0}),
            corr_table]),
    ])


def render_event_detail():
    d = DATA
    events_map = d.get("events_map", {})
    results = d.get("market_results", {})
    options = sorted(events_map.keys())

    return html.Div([
        html.Div(style=CARD_STYLE, children=[
            html.Label("Select Event:", style={"color": COLORS["text_muted"],
                                                "marginRight": "12px"}),
            dcc.Dropdown(
                id="event-dropdown",
                options=[{"label": e, "value": e} for e in options],
                value=options[0] if options else None,
                style={"backgroundColor": COLORS["card"], "color": "#000",
                       "minWidth": "300px"},
            ),
        ]),
        html.Div(id="event-detail-content"),
    ])


# ── App setup ────────────────────────────────────────────────────────
app = dash.Dash(__name__, title="MM Performance", suppress_callback_exceptions=True)

TAB_STYLE = {"color": COLORS["text_muted"], "backgroundColor": COLORS["card"],
             "border": "none", "padding": "12px 20px"}
TAB_SELECTED = {"color": COLORS["accent"], "backgroundColor": COLORS["bg"],
                "borderTop": f"2px solid {COLORS['accent']}", "padding": "12px 20px"}

app.layout = html.Div(
    style={"backgroundColor": COLORS["bg"], "minHeight": "100vh", "padding": "16px",
           "fontFamily": "-apple-system, BlinkMacSystemFont, 'Segoe UI', sans-serif",
           "color": COLORS["text"]},
    children=[
        # Header
        html.Div(style={"background": COLORS["header_bg"], "padding": "20px 28px",
                         "borderRadius": "12px", "marginBottom": "16px",
                         "display": "flex", "justifyContent": "space-between",
                         "alignItems": "center", "flexWrap": "wrap"}, children=[
            html.Div([
                html.H1("CBB 1H Market Maker", style={"margin": 0, "color": "white",
                                                        "fontSize": "1.6em"}),
                html.P(id="last-updated", style={"margin": "4px 0 0 0",
                                                   "color": "rgba(255,255,255,0.85)"}),
            ]),
            html.Button("Refresh", id="refresh-btn", n_clicks=0, style={
                "background": "rgba(255,255,255,0.2)", "border": "1px solid rgba(255,255,255,0.3)",
                "color": "white", "padding": "10px 20px", "borderRadius": "8px",
                "fontSize": "0.95em", "fontWeight": "600", "cursor": "pointer"}),
        ]),
        # Tabs
        dcc.Tabs(id="main-tabs", value="pnl", colors={
            "border": COLORS["card_border"], "primary": COLORS["accent"],
            "background": COLORS["card"]}, children=[
            dcc.Tab(label="P&L", value="pnl", style=TAB_STYLE, selected_style=TAB_SELECTED),
            dcc.Tab(label="Maker / Taker", value="mt", style=TAB_STYLE, selected_style=TAB_SELECTED),
            dcc.Tab(label="Fill Patterns", value="fills", style=TAB_STYLE, selected_style=TAB_SELECTED),
            dcc.Tab(label="Concentration", value="conc", style=TAB_STYLE, selected_style=TAB_SELECTED),
            dcc.Tab(label="Event Detail", value="detail", style=TAB_STYLE, selected_style=TAB_SELECTED),
        ]),
        html.Div(id="tab-content"),
        dcc.Loading(id="loading", type="circle", children=[
            html.Div(id="refresh-output", style={"display": "none"})]),
    ],
)


@callback(Output("tab-content", "children"), Input("main-tabs", "value"))
def render_tab(tab):
    if not DATA.get("loaded"):
        load_data()
    if tab == "pnl":
        return render_pnl()
    elif tab == "mt":
        return render_maker_taker()
    elif tab == "fills":
        return render_fill_patterns()
    elif tab == "conc":
        return render_concentration()
    elif tab == "detail":
        return render_event_detail()
    return html.Div("Unknown tab")


@callback(Output("refresh-output", "children"), Output("last-updated", "children"),
          Input("refresh-btn", "n_clicks"), prevent_initial_call=True)
def refresh(n):
    load_data()
    return "", f"Updated: {DATA.get('ts', '?')}"


@callback(Output("event-detail-content", "children"),
          Input("event-dropdown", "value"), prevent_initial_call=True)
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
        pnl = fill_pnl_cents(f, results.get(t, "")) * cnt / 100 if settled else None
        rows.append({
            "time": f["created_time"][:19].replace("T", " "),
            "strike": strike, "side": f["side"], "contracts": cnt,
            "cost": round(cost, 2),
            "taker": "Yes" if f.get("is_taker") else "",
            "result": results.get(t, "") or "open",
            "pnl": round(pnl, 2) if pnl is not None else "open",
        })
    total_cost = sum(r["cost"] for r in rows)
    total_pnl = sum(r["pnl"] for r in rows if isinstance(r["pnl"], (int, float)))
    return html.Div(style=CARD_STYLE, children=[
        html.Div(style={"display": "flex", "gap": "12px", "marginBottom": "12px"}, children=[
            stat_card("Fills", str(len(rows))),
            stat_card("Cost", fmt_dollars(total_cost)),
            stat_card("P&L", fmt_dollars(total_pnl) if total_pnl else "open",
                      COLORS["green"] if total_pnl and total_pnl > 0 else COLORS["red"]),
        ]),
        make_table(rows,
            [{"name": "Time", "id": "time"},
             {"name": "Strike", "id": "strike"},
             {"name": "Side", "id": "side"},
             {"name": "Cts", "id": "contracts", "type": "numeric"},
             {"name": "Cost", "id": "cost", "type": "numeric",
              "format": Format(precision=2, scheme=Scheme.fixed)},
             {"name": "Taker", "id": "taker"},
             {"name": "Result", "id": "result"},
             {"name": "P&L", "id": "pnl"}],
            "event-fills",
            [{"if": {"filter_query": '{pnl} contains "-"', "column_id": "pnl"},
              "color": COLORS["red"]},
             {"if": {"filter_query": '{result} = "yes"', "column_id": "result"},
              "color": COLORS["green"]},
             {"if": {"filter_query": '{result} = "no"', "column_id": "result"},
              "color": COLORS["red"]}],
            page_size=50),
    ])


if __name__ == "__main__":
    print("Loading data from Kalshi API...")
    load_data()
    print(f"Loaded {len(DATA.get('cbb_fills', []))} fills. Starting dashboard...")
    app.run(debug=False, host="0.0.0.0", port=8084)
