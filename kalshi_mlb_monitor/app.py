"""Kalshi MLB Bots Monitor — read-only Dash dashboard over both bots.

Run: ``python -m kalshi_mlb_monitor.app`` (or ``kalshi_mlb_monitor/run.sh``).
Opens at http://127.0.0.1:8092 (override KALSHI_MLB_MONITOR_PORT).

Read-only: never writes to any bot DB, imports no bot code. Each callback opens a
short-lived read-only DuckDB handle via ``queries.py`` and degrades gracefully if
the live bot briefly holds the write lock.
"""

from __future__ import annotations

import os
from datetime import datetime

import pandas as pd
import plotly.graph_objects as go
from dash import (Dash, Input, Output, State, callback, ctx, dash_table, dcc, html,
                  no_update)
from dash.exceptions import PreventUpdate

from . import queries as Q
from .bots import BOTS, gloss

# --------------------------------------------------------------------------- #
# Theme (mirrors kalshi_draft GitHub-dark palette)
# --------------------------------------------------------------------------- #
C = {
    "bg": "#0d1117", "card": "#161b22", "border": "#30363d",
    "text": "#c9d1d9", "muted": "#8b949e", "accent": "#00cec9",
    "accent2": "#0984e3", "green": "#3fb950", "red": "#f85149",
    "amber": "#d29922", "header": "linear-gradient(135deg,#00b894 0%,#00cec9 50%,#0984e3 100%)",
}
CARD = {"backgroundColor": C["card"], "border": f"1px solid {C['border']}",
        "borderRadius": "12px", "padding": "18px", "marginBottom": "14px"}
TAB = {"color": C["muted"], "backgroundColor": C["card"], "border": "none", "padding": "10px 18px"}
TAB_SEL = {"color": C["accent"], "backgroundColor": C["bg"],
           "borderTop": f"2px solid {C['accent']}", "padding": "10px 18px"}

DTABLE_KW = dict(
    style_header={"backgroundColor": "#0d1b2a", "color": C["accent"], "fontWeight": "700",
                  "fontSize": "0.8em", "textTransform": "uppercase", "border": "none",
                  "borderBottom": f"2px solid {C['accent']}"},
    style_cell={"backgroundColor": C["card"], "color": C["text"], "border": "none",
                "borderBottom": f"1px solid {C['border']}", "fontSize": "0.85em",
                "padding": "6px 10px", "textAlign": "left",
                "fontFamily": "ui-monospace,SFMono-Regular,Menlo,monospace",
                "maxWidth": "320px", "overflow": "hidden", "textOverflow": "ellipsis"},
    style_data_conditional=[{"if": {"row_index": "odd"}, "backgroundColor": "#1c2333"}],
    page_size=25, sort_action="native",
)


# --------------------------------------------------------------------------- #
# Small UI helpers
# --------------------------------------------------------------------------- #
def _fmt_int(v):
    return "—" if v is None or (isinstance(v, float) and pd.isna(v)) else f"{int(v):,}"


def _fmt_money(v):
    return "—" if v is None or (isinstance(v, float) and pd.isna(v)) else f"${float(v):,.2f}"


def _fmt_pct(v):
    return "—" if v is None or (isinstance(v, float) and pd.isna(v)) else f"{float(v)*100:.2f}%"


def _age_str(ts):
    if ts is None:
        return "never"
    secs = (datetime.now() - ts).total_seconds()
    if secs < 0:
        secs = 0
    if secs < 90:
        return f"{int(secs)}s ago"
    if secs < 5400:
        return f"{int(secs/60)}m ago"
    if secs < 172800:
        return f"{int(secs/3600)}h ago"
    return f"{int(secs/86400)}d ago"


def kpi_card(label, value, sub=None, color=None):
    return html.Div(style={**CARD, "flex": "1", "minWidth": "150px", "marginBottom": "0"}, children=[
        html.Div(label, style={"color": C["muted"], "fontSize": "0.72em",
                                "textTransform": "uppercase", "letterSpacing": "0.04em"}),
        html.Div(value, style={"color": color or C["accent"], "fontSize": "1.7em", "fontWeight": "700"}),
        html.Div(sub or "", style={"color": C["muted"], "fontSize": "0.72em"}),
    ])


def section(title, *children):
    return html.Div(style=CARD, children=[
        html.Div(title, style={"color": C["text"], "fontSize": "1.05em", "fontWeight": "700",
                               "marginBottom": "10px"}),
        *children,
    ])


def empty_state(msg):
    return html.Div(msg, style={"color": C["muted"], "padding": "26px", "textAlign": "center",
                                "fontStyle": "italic"})


def no_activity_note(bot_key, window):
    """Callout when the chosen window has no activity but the bot has history —
    so a dormant bot (e.g. the taker, down for weeks) reads as 'dormant', not
    'broken'. Returns None if 'all' is selected (nothing more to suggest)."""
    if window == "all":
        return None
    last = Q.status(bot_key)["last_activity"]
    return html.Div(style={**CARD, "borderColor": C["amber"]}, children=[
        html.Span("No activity in the last " + window + ". ",
                  style={"color": C["amber"], "fontWeight": "700"}),
        html.Span((f"This bot last traded {_age_str(last)}. " if last else ""),
                  style={"color": C["text"]}),
        html.Span("Switch the Window to ‘all’ to see its full history. "
                  "(Open positions and exposure below are all-time, not windowed.)",
                  style={"color": C["muted"]}),
    ])


def locked_state():
    return html.Div("⏳ data temporarily locked by the bot — showing nothing this tick",
                    style={"color": C["amber"], "padding": "20px", "textAlign": "center"})


def df_table(df, **kw):
    opts = {**DTABLE_KW, **kw}  # caller kwargs override defaults (no dup keys)
    return dash_table.DataTable(
        data=df.to_dict("records"),
        columns=[{"name": c, "id": c} for c in df.columns],
        **opts)


def _dark(fig, height=340):
    fig.update_layout(template="plotly_dark", paper_bgcolor=C["card"], plot_bgcolor=C["card"],
                      font_color=C["text"], height=height, margin=dict(l=10, r=10, t=30, b=10),
                      legend=dict(font=dict(size=10)))
    return fig


# --------------------------------------------------------------------------- #
# App
# --------------------------------------------------------------------------- #
app = Dash(__name__, title="Kalshi MLB Bots Monitor", update_title=None)
server = app.server

app.layout = html.Div(style={"backgroundColor": C["bg"], "minHeight": "100vh",
                             "fontFamily": "system-ui,-apple-system,sans-serif", "color": C["text"]}, children=[
    dcc.Store(id="bot", data="maker"),
    dcc.Store(id="window", data="24h"),
    dcc.Store(id="render-sig", data=None),  # last-rendered data fingerprint (poll guard)
    dcc.Interval(id="tick", interval=30_000, n_intervals=0),

    # Header
    html.Div(style={"background": C["header"], "padding": "16px 24px"}, children=[
        html.Div(style={"display": "flex", "alignItems": "center", "justifyContent": "space-between",
                        "flexWrap": "wrap", "gap": "12px"}, children=[
            html.Div([
                html.Span("⚾ Kalshi MLB Bots Monitor", style={"fontSize": "1.4em", "fontWeight": "800",
                                                              "color": "#fff"}),
                html.Span(" · read-only", style={"color": "rgba(255,255,255,.8)", "fontSize": "0.85em"}),
            ]),
            html.Div(style={"display": "flex", "gap": "8px", "alignItems": "center"}, children=[
                dcc.RadioItems(id="bot-select",
                               options=[{"label": f"  {b['label']}", "value": k} for k, b in BOTS.items()],
                               value="maker", inline=True,
                               style={"color": "#fff", "fontWeight": "700"},
                               inputStyle={"marginLeft": "12px", "marginRight": "4px"}),
            ]),
        ]),
    ]),

    # Status strip
    html.Div(id="status-strip", style={"padding": "10px 24px", "backgroundColor": C["card"],
                                        "borderBottom": f"1px solid {C['border']}", "display": "flex",
                                        "gap": "22px", "flexWrap": "wrap", "alignItems": "center",
                                        "fontSize": "0.85em"}),

    # Controls
    html.Div(style={"padding": "12px 24px", "display": "flex", "gap": "18px", "alignItems": "center",
                    "flexWrap": "wrap"}, children=[
        html.Span("Window:", style={"color": C["muted"]}),
        dcc.RadioItems(id="window-select",
                       options=[{"label": w, "value": w} for w in ["1h", "24h", "7d", "all"]],
                       value="24h", inline=True, style={"color": C["text"]},
                       inputStyle={"marginLeft": "10px", "marginRight": "4px"}),
        html.Span("Refresh:", style={"color": C["muted"], "marginLeft": "12px"}),
        dcc.RadioItems(id="refresh-select",
                       options=[{"label": "10s", "value": 10}, {"label": "30s", "value": 30},
                                {"label": "off", "value": 0}],
                       value=30, inline=True, style={"color": C["text"]},
                       inputStyle={"marginLeft": "10px", "marginRight": "4px"}),
        html.Button("↻ Refresh now", id="refresh-now", n_clicks=0,
                    style={"backgroundColor": C["card"], "color": C["accent"],
                           "border": f"1px solid {C['accent']}", "borderRadius": "8px",
                           "padding": "5px 12px", "cursor": "pointer"}),
        html.Span(id="updated-at", style={"color": C["muted"], "fontSize": "0.78em",
                                          "marginLeft": "auto"}),
    ]),

    # Tabs
    dcc.Tabs(id="tabs", value="overview", children=[
        dcc.Tab(label="Overview", value="overview", style=TAB, selected_style=TAB_SEL),
        dcc.Tab(label="Why Not Filled ★", value="why", style=TAB, selected_style=TAB_SEL),
        dcc.Tab(label="Fills & P&L", value="fills", style=TAB, selected_style=TAB_SEL),
        dcc.Tab(label="Positions & Exposure", value="positions", style=TAB, selected_style=TAB_SEL),
        dcc.Tab(label="Adverse Selection", value="adverse", style=TAB, selected_style=TAB_SEL),
    ]),
    # No dcc.Loading here on purpose: it flashed a spinner on every 30s refresh.
    # The poll guard below skips rebuilds when nothing changed, so updates are quiet.
    html.Div(id="tab-content", style={"padding": "16px 24px"}),
])


# --------------------------------------------------------------------------- #
# Control wiring
# --------------------------------------------------------------------------- #
@callback(Output("bot", "data"), Input("bot-select", "value"))
def _set_bot(v):
    return v


@callback(Output("window", "data"), Input("window-select", "value"))
def _set_window(v):
    return v


@callback(Output("tick", "interval"), Input("refresh-select", "value"))
def _set_interval(sec):
    return (sec or 0) * 1000 if sec else 24 * 3600 * 1000  # 0 => effectively off


# --------------------------------------------------------------------------- #
# Status strip
# --------------------------------------------------------------------------- #
@callback(Output("status-strip", "children"),
          Input("bot", "data"), Input("tick", "n_intervals"), Input("refresh-now", "n_clicks"))
def _status(bot_key, _n, _c):
    s = Q.status(bot_key)
    items = []
    # live dot
    if s["running"]:
        dot, txt, col = "●", "RUNNING", C["green"]
    else:
        dot, txt, col = "●", "NOT RUNNING", C["red"]
    items.append(html.Span([html.Span(dot + " ", style={"color": col}),
                            html.Span(txt, style={"color": col, "fontWeight": "700"})]))

    last = s["last_activity"]
    age = _age_str(last)
    if last is None:
        fcol = C["muted"]
    else:
        secs = (datetime.now() - last).total_seconds()
        fcol = C["green"] if secs < 120 else (C["amber"] if secs < 1800 else C["red"])
    items.append(html.Span([html.Span("last activity: ", style={"color": C["muted"]}),
                            html.Span(age, style={"color": fcol, "fontWeight": "700"})]))

    # "mode:" prefix so this never reads as a contradiction next to NOT RUNNING —
    # it describes dry-run vs live trading, not whether the process is up.
    if s["dry_run"] is True:
        items.append(html.Span("mode: DRY-RUN", style={"color": C["amber"], "fontWeight": "700"}))
    elif s["dry_run"] is False:
        items.append(html.Span("mode: LIVE", style={"color": C["green"], "fontWeight": "700"}))

    if s["session_start"]:
        items.append(html.Span([html.Span("session since: ", style={"color": C["muted"]}),
                                html.Span(s["session_start"].strftime("%b %d %H:%M"),
                                          style={"color": C["text"]})]))

    # Stale banner: not running AND last activity old
    stale = (not s["running"]) and last is not None and (datetime.now() - last).total_seconds() > 3600
    if stale:
        items.append(html.Span(f"⚠ STALE — bot down, last write {age}",
                               style={"color": C["red"], "fontWeight": "700",
                                      "backgroundColor": "rgba(248,81,73,.12)",
                                      "padding": "2px 10px", "borderRadius": "6px"}))
    return items


# --------------------------------------------------------------------------- #
# Tab renderer
# --------------------------------------------------------------------------- #
@callback(Output("tab-content", "children"), Output("render-sig", "data"),
          Output("updated-at", "children"),
          Input("tabs", "value"), Input("bot", "data"), Input("window", "data"),
          Input("tick", "n_intervals"), Input("refresh-now", "n_clicks"),
          State("render-sig", "data"))
def _render(tab, bot_key, window, _n, _c, prev_sig):
    # If the bot DB can't be read right now (genuinely stuck lock after retries),
    # don't render dashes — that would make a live bot look dormant. Keep the last
    # good screen if we have one; otherwise show a transient "reading…" placeholder.
    sig_data = Q.signature(bot_key)
    if sig_data is None:
        if prev_sig is not None:
            raise PreventUpdate
        return (html.Div("⏳ reading… (bot is writing; will refresh momentarily)",
                         style={"color": C["amber"], "padding": "40px", "textAlign": "center"}),
                no_update, "reading…")

    # Poll guard: if this fire came from the auto-refresh interval and the bot's
    # data fingerprint is unchanged since the last render, don't rebuild — that
    # rebuild is exactly what made the page look like it was constantly buffering.
    # User actions (switching tab/bot/window, clicking Refresh now) always render.
    sig = f"{tab}|{bot_key}|{window}|{sig_data}"
    if ctx.triggered_id == "tick" and prev_sig is not None and sig == prev_sig:
        raise PreventUpdate

    stamp = "updated " + datetime.now().strftime("%H:%M:%S")
    cutoff = Q.window_cutoff(window)
    builders = {
        "overview": lambda: _tab_overview(bot_key, window, cutoff),
        "why": lambda: _tab_why(bot_key, window, cutoff),
        "fills": lambda: _tab_fills(bot_key),
        "positions": lambda: _tab_positions(bot_key),
        "adverse": lambda: _tab_adverse(bot_key, cutoff),
    }
    build = builders.get(tab)
    if build is None:
        raise PreventUpdate
    return build(), sig, stamp


# ---- Tab 1: Overview ---- #
def _tab_overview(bot_key, window, cutoff):
    bot = BOTS[bot_key]
    k = Q.kpis(bot_key, cutoff)
    if bot["kind"] == "maker":
        cards = [
            kpi_card("RFQs seen", _fmt_int(k["rfqs"]), f"in {window}"),
            kpi_card("In scope", _fmt_int(k["in_scope"]),
                     _fmt_pct((k["in_scope"] / k["rfqs"]) if k["rfqs"] else None) + " of seen"),
            kpi_card("Quoted", _fmt_int(k["quotes"])),
            kpi_card("Fills", _fmt_int(k["fills"]),
                     "fill rate " + _fmt_pct(k["fill_rate"]),
                     color=C["green"] if (k["fills"] or 0) > 0 else C["red"]),
            kpi_card("Open positions", _fmt_int(k["open_positions"]),
                     _fmt_money(k["exposure"]) + " exposure"),
        ]
    else:
        cards = [
            kpi_card("RFQs sent", _fmt_int(k["rfqs"]), f"in {window}"),
            kpi_card("Quotes evaluated", _fmt_int(k["in_scope"])),
            kpi_card("Fills", _fmt_int(k["fills"]),
                     "fill rate " + _fmt_pct(k["fill_rate"]),
                     color=C["green"] if (k["fills"] or 0) > 0 else C["muted"]),
            kpi_card("Staked", _fmt_money(k["staked"])),
            kpi_card("Open positions", _fmt_int(k["open_positions"]),
                     _fmt_money(k["exposure"]) + " exposure"),
        ]
    kpi_row = html.Div(cards, style={"display": "flex", "gap": "12px", "flexWrap": "wrap",
                                     "marginBottom": "14px"})

    # Dormant-in-window guard: if nothing happened in the window but the bot has
    # history, lead with a note so empty KPIs don't read as a broken dashboard.
    dormant = no_activity_note(bot_key, window) if not (k["rfqs"] or 0) else None

    # Funnel
    stages = Q.funnel(bot_key, cutoff)
    labels = [s[0] for s in stages]
    vals = [int(s[1]) if s[1] is not None else 0 for s in stages]
    fig = go.Figure(go.Funnel(y=labels, x=vals, textinfo="value+percent initial",
                              marker={"color": [C["accent"], C["accent2"], "#6f42c1", C["green"]][:len(vals)]}))
    funnel_card = section("Conversion funnel — where RFQs drop off", dcc.Graph(figure=_dark(fig, 320)))

    note = ""
    if bot["kind"] == "maker" and (k["fills"] or 0) == 0:
        note = section("Reading this",
                       html.Div([
                           "The maker has seen ", html.B(_fmt_int(k["rfqs"])),
                           " RFQs but only ", html.B(_fmt_int(k["in_scope"])),
                           " were in scope (spread×total combos) and ", html.B("0"),
                           " filled. The bottleneck is at the very first stage — almost all "
                           "incoming RFQ flow isn't the kind of market this bot quotes. "
                           "See the ", html.B("Why Not Filled"), " tab for the full breakdown."],
                           style={"color": C["text"], "lineHeight": "1.6"}))
    children = [c for c in (dormant, kpi_row, funnel_card, note) if c]
    return html.Div(children)


# ---- Tab 2: Why Not Filled ---- #
def _tab_why(bot_key, window, cutoff):
    bd = Q.decision_breakdown(bot_key, cutoff)
    if bd is Q.LOCKED:
        return locked_state()
    if not len(bd):
        note = no_activity_note(bot_key, window)
        return note or empty_state("No decisions recorded.")

    total = bd["n"].sum()
    bd = bd.copy()
    bd["pct"] = bd["n"] / total
    bd["what it means"] = bd["category"].map(gloss)

    bar = go.Figure(go.Bar(x=bd["n"], y=bd["category"], orientation="h",
                           marker_color=C["accent"], text=bd["n"], textposition="auto"))
    bar.update_layout(yaxis={"autorange": "reversed"}, xaxis_title="count")
    bar_card = section(f"Why not filled — decision breakdown ({window}, {int(total):,} decisions)",
                       dcc.Graph(figure=_dark(bar, max(260, 40 * len(bd)))))

    tbl = bd.assign(count=bd["n"].map(lambda v: f"{int(v):,}"),
                    share=bd["pct"].map(lambda v: f"{v*100:.1f}%"))[
        ["category", "count", "share", "what it means"]]
    tbl_card = section("Reason legend", df_table(tbl, page_size=20))

    # time-series stacked
    ts = Q.decision_timeseries(bot_key, cutoff, window)
    if ts is Q.LOCKED or not len(ts):
        ts_card = section("Decisions over time", empty_state("No data in window."))
    else:
        top = list(bd["category"].head(6))
        fig = go.Figure()
        for cat in top:
            sub = ts[ts["category"] == cat]
            fig.add_trace(go.Scatter(x=sub["bucket"], y=sub["n"], name=cat,
                                     stackgroup="one", mode="lines"))
        ts_card = section("Decisions over time (top reasons)", dcc.Graph(figure=_dark(fig, 320)))

    return html.Div([bar_card, ts_card, tbl_card])


# ---- Tab 3: Fills & P&L ---- #
def _tab_fills(bot_key):
    bot = BOTS[bot_key]
    rf = Q.recent_fills(bot_key, 300)
    if rf is Q.LOCKED:
        return locked_state()
    if not len(rf):
        return section("Fills", empty_state(
            "0 fills recorded. " + ("This maker has never been filled — see the "
                                    "Why Not Filled tab for why." if bot["kind"] == "maker" else "")))

    cum = Q.fills_cumulative(bot_key)
    charts = []
    if cum is not Q.LOCKED and len(cum):
        f = go.Figure()
        f.add_trace(go.Scatter(x=cum["filled_at"], y=cum["cum_fills"], name="cumulative fills",
                               mode="lines", line=dict(color=C["accent"])))
        f.add_trace(go.Scatter(x=cum["filled_at"], y=cum["cum_stake"], name="cumulative stake $",
                               mode="lines", line=dict(color=C["green"]), yaxis="y2"))
        f.update_layout(yaxis2=dict(overlaying="y", side="right", title="stake $"),
                        yaxis_title="fills")
        charts.append(section("Cumulative fills & stake", dcc.Graph(figure=_dark(f, 320))))

    disp = rf.copy()
    if "filled_at" in disp:
        disp["filled_at"] = pd.to_datetime(disp["filled_at"]).dt.strftime("%m-%d %H:%M:%S")
    for col in ("price", "fee", "fair_at_quote", "realized_pnl"):
        if col in disp:
            disp[col] = disp[col].map(lambda v: "" if pd.isna(v) else f"{float(v):.3f}")
    tbl = section(f"Recent fills ({len(rf)})", df_table(disp))
    return html.Div(charts + [tbl])


# ---- Tab 4: Positions & Exposure ---- #
def _tab_positions(bot_key):
    bot = BOTS[bot_key]
    pos = Q.positions(bot_key)
    summ = Q.open_orders_summary(bot_key)
    oo = Q.open_orders(bot_key)
    blocks = []

    # Open working orders summary + orphan flag
    if summ is not Q.LOCKED and len(summ):
        chips = []
        for _, r in summ.iterrows():
            col = C["amber"] if r["status"] == "open" and not Q.is_running(bot_key) else C["text"]
            chips.append(html.Span(f"{r['status']}: {int(r['n'])}",
                                   style={"color": col, "marginRight": "16px", "fontWeight": "600"}))
        orphan = (not Q.is_running(bot_key)) and (summ["status"] == "open").any()
        if orphan:
            n_open = int(summ.loc[summ["status"] == "open", "n"].iloc[0])
            chips.append(html.Span(f"⚠ {n_open} orphaned 'open' (bot not running)",
                                   style={"color": C["red"], "fontWeight": "700"}))
        blocks.append(section("Working orders", html.Div(chips)))

    if pos is Q.LOCKED:
        blocks.append(locked_state())
    elif not len(pos):
        blocks.append(section("Open positions", empty_state("No open positions.")))
    else:
        total_exp = pos["exposure"].sum()
        disp = pos.copy()
        if "updated_at" in disp:
            disp["updated_at"] = pd.to_datetime(disp["updated_at"]).dt.strftime("%m-%d %H:%M")
        for col in ("weighted_price", "exposure", "net_contracts"):
            if col in disp:
                disp[col] = disp[col].map(lambda v: "" if pd.isna(v) else f"{float(v):.3f}")
        if "legs_json" in disp:
            disp["legs_json"] = disp["legs_json"].astype(str).str.slice(0, 60)
        blocks.append(section(f"Open positions ({len(pos)}) — {_fmt_money(total_exp)} total exposure",
                              df_table(disp)))

    # detail of open orders
    if oo is not Q.LOCKED and len(oo):
        disp = oo.copy()
        for c in disp.columns:
            if "_at" in c:
                disp[c] = pd.to_datetime(disp[c], errors="coerce").dt.strftime("%m-%d %H:%M")
        blocks.append(section("Working-order detail", df_table(disp)))
    return html.Div(blocks)


# ---- Tab 5: Adverse Selection ---- #
def _tab_adverse(bot_key, cutoff):
    bot = BOTS[bot_key]
    if bot["kind"] == "maker":
        adv = Q.adverse(bot_key)
        if adv is Q.LOCKED:
            return locked_state()
        if not len(adv):
            return section("Adverse selection (measurement scaffold)", empty_state(
                "No fills yet, so there's nothing to measure. Once the maker fills, this shows "
                "quoted margin vs realized fair drift (fair_at_confirm − fair_at_quote) per fill — "
                "the core v1 question: does the 5% margin survive adverse selection?"))
        adv = adv.copy()
        f = go.Figure(go.Scatter(x=adv["filled_at"], y=adv["fair_drift"], mode="markers",
                                 marker=dict(color=C["accent"])))
        f.update_layout(yaxis_title="fair drift (confirm − quote)")
        return section("Fair drift per fill (adverse-selection signal)", dcc.Graph(figure=_dark(f, 340)))

    # taker
    qq = Q.taker_quote_quality(bot_key, cutoff)
    if qq is Q.LOCKED:
        return locked_state()
    if not len(qq):
        return section("Fill quality", empty_state("No accept/walk/halt events in window."))
    f = go.Figure()
    for dec, col in [("accepted", C["green"]), ("failed_quote_walked", C["red"]),
                     ("halted_low_fill_ratio", C["amber"])]:
        sub = qq[qq["decision"] == dec]
        if len(sub):
            f.add_trace(go.Bar(x=sub["bucket"], y=sub["n"], name=dec, marker_color=col))
    f.update_layout(barmode="stack", yaxis_title="count/hour")
    note = html.Div("Walks = quotes that walked when we tried to accept (we lost the race). "
                    "Halts = trading paused because the fill ratio dropped.",
                    style={"color": C["muted"], "fontSize": "0.85em", "marginTop": "8px"})
    return section("Accept vs walk vs halt over time", dcc.Graph(figure=_dark(f, 340)), note)


def main():
    port = int(os.environ.get("KALSHI_MLB_MONITOR_PORT", "8092"))
    app.run(host="127.0.0.1", port=port, debug=False)


if __name__ == "__main__":
    main()
