"""
NFL Draft Prediction Market Dashboard
Dash application with tabs for market overview, price history,
edge detection, consensus comparison, and portfolio.
"""

import subprocess
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))
from nfl_draft.lib import db as nfl_db
from nfl_draft.lib import queries as nfl_queries
from nfl_draft.lib.queries import QueryLocked

import dash
from dash import dcc, html, dash_table, callback, Input, Output, State
from dash.dash_table.Format import Format
from dash.exceptions import PreventUpdate
import plotly.express as px
import plotly.graph_objects as go
import pandas as pd
import uuid
from datetime import datetime

import db

# ---------------------------------------------------------------------------
# Series display names
# ---------------------------------------------------------------------------
SERIES_DISPLAY = {
    "KXNFLDRAFT1": "#1 Overall Pick",
    "KXNFLDRAFT1ST": "Team with #1 Pick",
    "KXNFLDRAFTPICK": "Draft Pick Position",
    "KXNFLDRAFTTOP": "Top N Draft Range",
    "KXNFLFIRSTPICK": "First Pick",
}


def display_series(ticker):
    """Convert series ticker to human-readable name."""
    return SERIES_DISPLAY.get(ticker, ticker)


# ---------------------------------------------------------------------------
# App init
# ---------------------------------------------------------------------------
app = dash.Dash(
    __name__,
    title="NFL Draft Dashboard",
    suppress_callback_exceptions=True,
)

# ---------------------------------------------------------------------------
# Dark theme colors
# ---------------------------------------------------------------------------
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
    "header_bg": "linear-gradient(135deg, #00b894 0%, #00cec9 50%, #0984e3 100%)",
}

CARD_STYLE = {
    "backgroundColor": COLORS["card"],
    "border": f"1px solid {COLORS['card_border']}",
    "borderRadius": "12px",
    "padding": "20px",
    "marginBottom": "16px",
}

TABLE_STYLE_HEADER = {
    "backgroundColor": "#0d1b2a",
    "color": COLORS["accent"],
    "fontWeight": "700",
    "fontSize": "0.85em",
    "textTransform": "uppercase",
    "border": "none",
    "borderBottom": f"2px solid {COLORS['accent']}",
}

TABLE_STYLE_DATA = {
    "backgroundColor": COLORS["card"],
    "color": COLORS["text"],
    "border": "none",
    "borderBottom": f"1px solid {COLORS['card_border']}",
}

TABLE_STYLE_DATA_CONDITIONAL = [
    {"if": {"row_index": "odd"}, "backgroundColor": "#1c2333"},
]


# ---------------------------------------------------------------------------
# Layout
# ---------------------------------------------------------------------------
def make_header():
    snapshot = db.get_snapshot_count()
    last_fetch = snapshot[2] if snapshot and snapshot[2] else "Never"
    return html.Div(
        style={
            "background": COLORS["header_bg"],
            "padding": "24px 32px",
            "borderRadius": "12px",
            "marginBottom": "20px",
            "display": "flex",
            "justifyContent": "space-between",
            "alignItems": "center",
        },
        children=[
            html.Div([
                html.H1("NFL Draft Dashboard", style={
                    "margin": "0", "color": "white", "fontSize": "1.8em",
                }),
                html.P(f"Updated: {last_fetch}", id="last-updated", style={
                    "margin": "4px 0 0 0", "color": "rgba(255,255,255,0.85)",
                }),
            ]),
            html.Button("Refresh Data", id="refresh-btn", n_clicks=0, style={
                "background": "rgba(255,255,255,0.2)",
                "border": "1px solid rgba(255,255,255,0.3)",
                "color": "white",
                "padding": "10px 20px",
                "borderRadius": "8px",
                "fontSize": "0.95em",
                "fontWeight": "600",
                "cursor": "pointer",
            }),
        ],
    )


def make_stats_cards(df):
    """Summary stat cards from latest odds."""
    if df is None or df.empty:
        return html.Div("No data yet. Click Refresh.", style={"color": COLORS["text_muted"]})

    n_series = df["series_ticker"].nunique()
    n_markets = len(df)
    total_vol = df["volume"].sum()
    total_liq = df["liquidity"].sum()

    cards_data = [
        ("Series", str(n_series)),
        ("Markets", f"{n_markets:,}"),
        ("Total Volume", f"{total_vol:,.0f}"),
        ("Total Liquidity", f"${total_liq:,.0f}"),
    ]

    return html.Div(
        style={"display": "flex", "gap": "16px", "marginBottom": "16px"},
        children=[
            html.Div(style={
                **CARD_STYLE,
                "flex": "1",
                "textAlign": "center",
            }, children=[
                html.Div(label, style={"color": COLORS["text_muted"], "fontSize": "0.8em", "textTransform": "uppercase"}),
                html.Div(value, style={"color": COLORS["accent"], "fontSize": "1.6em", "fontWeight": "700"}),
            ])
            for label, value in cards_data
        ],
    )


# Common styling helpers for tabs to stay DRY
_TAB_STYLE = {"color": COLORS["text_muted"], "backgroundColor": COLORS["card"], "border": "none", "padding": "12px 20px"}
_TAB_SELECTED = {"color": COLORS["accent"], "backgroundColor": COLORS["bg"], "borderTop": f"2px solid {COLORS['accent']}", "padding": "12px 20px"}


def _make_header_with_mode_toggle():
    """Header with mode toggle (pre-draft = 60s poll, draft-day = 15s poll)."""
    snapshot = db.get_snapshot_count()
    last_fetch = snapshot[2] if snapshot and snapshot[2] else "Never"
    return html.Div(
        style={
            "background": COLORS["header_bg"],
            "padding": "24px 32px",
            "borderRadius": "12px",
            "marginBottom": "20px",
            "display": "flex",
            "justifyContent": "space-between",
            "alignItems": "center",
        },
        children=[
            html.Div([
                html.H1("NFL Draft Dashboard", style={
                    "margin": "0", "color": "white", "fontSize": "1.8em",
                }),
                html.P(f"Updated: {last_fetch}", id="last-updated", style={
                    "margin": "4px 0 0 0", "color": "rgba(255,255,255,0.85)",
                }),
            ]),
            html.Div([
                dcc.RadioItems(
                    id="nfl_draft__mode_toggle_radio",
                    options=[
                        {"label": " Pre-draft (60s) ", "value": "pre"},
                        {"label": " Draft-day (15s) ", "value": "draft"},
                    ],
                    value="pre",
                    inline=True,
                    style={"color": "white", "marginRight": "16px"},
                ),
                html.Button("Refresh Data", id="refresh-btn", n_clicks=0, style={
                    "background": "rgba(255,255,255,0.2)",
                    "border": "1px solid rgba(255,255,255,0.3)",
                    "color": "white",
                    "padding": "10px 20px",
                    "borderRadius": "8px",
                    "fontSize": "0.95em",
                    "fontWeight": "600",
                    "cursor": "pointer",
                }),
            ], style={"display": "flex", "alignItems": "center"}),
        ],
    )


def _portal_tabs():
    """Inner tabs under the 'Portal' section — 4 new cross-venue tabs."""
    return dcc.Tabs(
        id="portal-tabs",
        value="crossbook",
        colors={
            "border": COLORS["card_border"],
            "primary": COLORS["accent"],
            "background": COLORS["card"],
        },
        style={"marginBottom": "16px"},
        children=[
            dcc.Tab(label="Cross-Book Grid", value="crossbook", style=_TAB_STYLE, selected_style=_TAB_SELECTED),
            dcc.Tab(label="+EV Candidates", value="ev", style=_TAB_STYLE, selected_style=_TAB_SELECTED),
            dcc.Tab(label="Trade Tape", value="tape", style=_TAB_STYLE, selected_style=_TAB_SELECTED),
            dcc.Tab(label="Bet Log", value="betlog", style=_TAB_STYLE, selected_style=_TAB_SELECTED),
        ],
    )


def _legacy_tabs():
    """Inner tabs under the 'Kalshi-only legacy' section — the original 5 tabs."""
    return dcc.Tabs(
        id="main-tabs",
        value="overview",
        colors={
            "border": COLORS["card_border"],
            "primary": COLORS["accent"],
            "background": COLORS["card"],
        },
        style={"marginBottom": "16px"},
        children=[
            dcc.Tab(label="Market Overview", value="overview", style=_TAB_STYLE, selected_style=_TAB_SELECTED),
            dcc.Tab(label="Price History", value="history", style=_TAB_STYLE, selected_style=_TAB_SELECTED),
            dcc.Tab(label="Edge Detection", value="edges", style=_TAB_STYLE, selected_style=_TAB_SELECTED),
            dcc.Tab(label="Consensus", value="consensus", style=_TAB_STYLE, selected_style=_TAB_SELECTED),
            dcc.Tab(label="Portfolio", value="portfolio", style=_TAB_STYLE, selected_style=_TAB_SELECTED),
        ],
    )


app.layout = html.Div(
    style={
        "backgroundColor": COLORS["bg"],
        "minHeight": "100vh",
        "padding": "20px",
        "fontFamily": "-apple-system, BlinkMacSystemFont, 'Segoe UI', Roboto, sans-serif",
        "color": COLORS["text"],
    },
    children=[
        # dcc.Store registry — prefixed with 'nfl_draft.' per spec
        dcc.Store(id="nfl_draft__mode_toggle", storage_type="local", data="pre"),
        dcc.Store(id="nfl_draft__subnav", storage_type="local"),
        dcc.Store(id="nfl_draft__last_fetched_odds", storage_type="local"),
        dcc.Store(id="nfl_draft__last_fetched_trades", storage_type="local"),
        dcc.Store(id="nfl_draft__bet_log_prefill", storage_type="session"),
        # Auto-refresh interval. 60s default (pre-draft), switched to 15s by mode toggle.
        dcc.Interval(id="nfl_draft__interval", interval=60_000, n_intervals=0),

        _make_header_with_mode_toggle(),

        # Outer Portal / Legacy section selector
        dcc.Tabs(
            id="section-tabs",
            value="portal",
            colors={
                "border": COLORS["card_border"],
                "primary": COLORS["accent2"],
                "background": COLORS["card"],
            },
            style={"marginBottom": "12px"},
            children=[
                dcc.Tab(label="Portal (Cross-Venue)", value="portal",
                        style=_TAB_STYLE, selected_style=_TAB_SELECTED),
                dcc.Tab(label="Kalshi-only (Legacy)", value="legacy",
                        style=_TAB_STYLE, selected_style=_TAB_SELECTED),
            ],
        ),
        html.Div(id="section-content"),
        html.Div(id="tab-content"),
        # Hidden div for refresh callback
        html.Div(id="refresh-output", style={"display": "none"}),
        # Toast notification area
        html.Div(id="toast-area", style={
            "position": "fixed", "bottom": "30px", "right": "30px", "zIndex": "2000",
        }),
    ],
)


# ---------------------------------------------------------------------------
# Tab content renderers
# ---------------------------------------------------------------------------

def render_overview():
    df = db.get_latest_odds()
    if df is None or df.empty:
        return html.Div("No data. Click Refresh to fetch.", style={"color": COLORS["text_muted"], "padding": "40px"})

    stats = make_stats_cards(df)

    # --- Heatmap: Players × Draft Position (1st through 5th) ---
    # Combine KXNFLDRAFT1 (pick 1) with KXNFLDRAFTPICK (picks 2-5)
    pick_map = {
        "KXNFLDRAFT1-26": "1st",
        "KXNFLDRAFTPICK-26-2": "2nd",
        "KXNFLDRAFTPICK-26-3": "3rd",
        "KXNFLDRAFTPICK-26-4": "4th",
        "KXNFLDRAFTPICK-26-5": "5th",
    }
    pick_events = set(pick_map.keys())
    heatmap_df = df[df["event_ticker"].isin(pick_events)].copy()
    heatmap_df["pick"] = heatmap_df["event_ticker"].map(pick_map)

    if not heatmap_df.empty:
        # Pivot: players as rows, pick positions as columns
        # Build separate pivots for bid, ask, and midpoint (for color)
        pivot_bid = heatmap_df.pivot_table(
            index="candidate", columns="pick",
            values="yes_bid", aggfunc="first"
        ).fillna(0)
        pivot_ask = heatmap_df.pivot_table(
            index="candidate", columns="pick",
            values="yes_ask", aggfunc="first"
        ).fillna(0)

        # Reorder columns
        col_order = [c for c in ["1st", "2nd", "3rd", "4th", "5th"] if c in pivot_bid.columns]
        pivot_bid = pivot_bid[col_order]
        pivot_ask = pivot_ask.reindex(columns=col_order, fill_value=0)

        # Midpoint for color scale
        pivot_mid = (pivot_bid + pivot_ask) / 2

        # Sort players by highest midpoint in any pick slot
        pivot_mid["_max"] = pivot_mid.max(axis=1)
        sort_order = pivot_mid.sort_values("_max", ascending=True).drop("_max", axis=1).index
        pivot_bid = pivot_bid.loc[sort_order]
        pivot_ask = pivot_ask.loc[sort_order]
        pivot_mid = pivot_mid.drop("_max", axis=1).loc[sort_order]

        # Only show players with >= 3% midpoint somewhere
        mask = pivot_mid.max(axis=1) >= 3
        pivot_bid = pivot_bid[mask]
        pivot_ask = pivot_ask[mask]
        pivot_mid = pivot_mid[mask]

        if not pivot_mid.empty:
            # Build bid/ask text annotations
            text_labels = []
            for i in range(len(pivot_bid)):
                row_labels = []
                for j in range(len(col_order)):
                    b = int(pivot_bid.iloc[i, j])
                    a = int(pivot_ask.iloc[i, j])
                    if b == 0 and a == 0:
                        row_labels.append("")
                    else:
                        row_labels.append(f"{b}/{a}")
                text_labels.append(row_labels)

            heatmap_fig = go.Figure(data=go.Heatmap(
                z=pivot_mid.values,
                x=pivot_mid.columns.tolist(),
                y=pivot_mid.index.tolist(),
                colorscale=[
                    [0, "#0d1b2a"],
                    [0.01, "#1a472a"],
                    [0.15, "#2ecc71"],
                    [0.5, "#f1c40f"],
                    [1.0, "#e74c3c"],
                ],
                zmin=0,
                zmax=100,
                text=text_labels,
                texttemplate="%{text}",
                textfont={"size": 16, "color": "white"},
                hovertemplate="<b>%{y}</b><br>Pick: %{x}<br>Bid/Ask: %{text}<br>Midpoint: %{z:.0f}%<extra></extra>",
                colorbar=dict(title="Mid %", ticksuffix="%"),
            ))
            heatmap_fig.update_layout(
                template="plotly_dark",
                paper_bgcolor=COLORS["card"],
                plot_bgcolor=COLORS["card"],
                height=max(500, len(pivot_mid) * 40 + 120),
                margin=dict(l=220, r=40, t=40, b=60),
                xaxis=dict(side="top", title="Draft Position"),
                yaxis=dict(title=None),
            )
            heatmap = dcc.Graph(figure=heatmap_fig)
        else:
            heatmap = html.Div("No candidates above threshold", style={"color": COLORS["text_muted"]})
    else:
        heatmap = html.Div("No draft position data for heatmap", style={"color": COLORS["text_muted"]})

    # --- Odds table ---
    table_df = df[["series_ticker", "candidate", "last_price", "yes_bid", "yes_ask",
                    "volume", "liquidity", "open_interest"]].copy()
    table_df["series"] = table_df["series_ticker"].map(display_series)
    table_df["bid_ask"] = table_df["yes_bid"].astype(str) + "/" + table_df["yes_ask"].astype(str)
    table_df = table_df.sort_values(["series_ticker", "last_price"], ascending=[True, False])

    odds_table = dash_table.DataTable(
        data=table_df.to_dict("records"),
        columns=[
            {"name": "Series", "id": "series"},
            {"name": "Candidate", "id": "candidate"},
            {"name": "Prob %", "id": "last_price", "type": "numeric"},
            {"name": "Bid/Ask", "id": "bid_ask"},
            {"name": "Volume", "id": "volume", "type": "numeric", "format": Format().group(True)},
            {"name": "Liquidity", "id": "liquidity", "type": "numeric", "format": dash_table.FormatTemplate.money(0)},
            {"name": "Open Int", "id": "open_interest", "type": "numeric", "format": Format().group(True)},
        ],
        filter_action="native",
        sort_action="native",
        sort_mode="multi",
        page_size=25,
        style_header=TABLE_STYLE_HEADER,
        style_data=TABLE_STYLE_DATA,
        style_data_conditional=TABLE_STYLE_DATA_CONDITIONAL,
        style_filter={"backgroundColor": "#0d1b2a", "color": COLORS["text"]},
        style_table={"overflowX": "auto"},
    )

    return html.Div([
        stats,
        html.Div(style=CARD_STYLE, children=[
            html.H3("Draft Position Probability", style={"color": COLORS["accent"], "marginTop": 0}),
            html.P("Probability of each player being selected at each draft position (1st-5th). Hover for details.", style={"color": COLORS["text_muted"], "fontSize": "0.85em"}),
            heatmap,
        ]),
        html.Div(style=CARD_STYLE, children=[
            html.H3("All Markets", style={"color": COLORS["accent"], "marginTop": 0}),
            odds_table,
        ]),
    ])


def render_history():
    df = db.get_latest_odds()
    if df is None or df.empty:
        return html.Div("No data yet.", style={"color": COLORS["text_muted"], "padding": "40px"})

    # Build options from all candidates with > 0% probability
    candidates = df[df["last_price"] > 0].sort_values("last_price", ascending=False)
    options = [
        {"label": f"{row['candidate']} ({display_series(row['series_ticker'])}: {row['last_price']}%)", "value": row["ticker"]}
        for _, row in candidates.iterrows()
    ]

    return html.Div([
        html.Div(style=CARD_STYLE, children=[
            html.H3("Price History", style={"color": COLORS["accent"], "marginTop": 0}),
            html.P("Select markets to compare price movement over time.", style={"color": COLORS["text_muted"], "fontSize": "0.85em"}),
            dcc.Dropdown(
                id="history-selector",
                options=options,
                multi=True,
                placeholder="Select players/markets...",
                value=[options[0]["value"]] if options else [],
                style={"backgroundColor": "#0d1b2a", "color": COLORS["text"], "marginBottom": "16px"},
            ),
            dcc.Graph(id="history-chart"),
        ]),
    ])


def render_edges():
    edges_df = db.get_latest_edges()

    if edges_df is None or edges_df.empty:
        return html.Div(style=CARD_STYLE, children=[
            html.H3("Edge Detection", style={"color": COLORS["accent"], "marginTop": 0}),
            html.P("No edges detected yet. Click Refresh to run edge detection.", style={"color": COLORS["text_muted"]}),
        ])

    # Color code by confidence
    style_conditional = [
        {"if": {"filter_query": '{confidence} = "high"', "column_id": "confidence"},
         "color": COLORS["green"], "fontWeight": "bold"},
        {"if": {"filter_query": '{confidence} = "medium"', "column_id": "confidence"},
         "color": "#f0ad4e"},
        {"if": {"filter_query": '{confidence} = "low"', "column_id": "confidence"},
         "color": COLORS["text_muted"]},
    ] + TABLE_STYLE_DATA_CONDITIONAL

    edges_table = dash_table.DataTable(
        data=edges_df.to_dict("records"),
        columns=[
            {"name": "Type", "id": "edge_type"},
            {"name": "Description", "id": "description"},
            {"name": "Market A", "id": "market_a"},
            {"name": "Market B", "id": "market_b"},
            {"name": "Price A", "id": "price_a", "type": "numeric"},
            {"name": "Price B", "id": "price_b", "type": "numeric"},
            {"name": "Edge", "id": "implied_edge", "type": "numeric"},
            {"name": "Confidence", "id": "confidence"},
        ],
        sort_action="native",
        sort_mode="multi",
        style_header=TABLE_STYLE_HEADER,
        style_data=TABLE_STYLE_DATA,
        style_data_conditional=style_conditional,
        style_table={"overflowX": "auto"},
    )

    # Sum-to-1 chart
    latest = db.get_latest_odds()
    sum_chart = None
    if latest is not None and not latest.empty:
        sums = latest.groupby("series_ticker")["last_price"].sum().reset_index()
        sums.columns = ["series_ticker", "total_implied"]
        sums["series"] = sums["series_ticker"].map(display_series)
        sums["overround"] = sums["total_implied"] - 100

        fig = go.Figure()
        fig.add_bar(
            x=sums["series"], y=sums["overround"],
            marker_color=[COLORS["green"] if v < 0 else COLORS["red"] for v in sums["overround"]],
            text=[f"{v:+.1f}%" for v in sums["overround"]],
            textposition="outside",
        )
        fig.add_hline(y=0, line_dash="dash", line_color=COLORS["text_muted"])
        fig.update_layout(
            template="plotly_dark",
            paper_bgcolor=COLORS["card"],
            plot_bgcolor=COLORS["card"],
            title="Sum-to-1 Overround by Series",
            yaxis_title="Overround %",
            height=350,
            margin=dict(t=60, b=40),
        )
        sum_chart = dcc.Graph(figure=fig)

    return html.Div([
        html.Div(style=CARD_STYLE, children=[
            html.H3("Detected Edges", style={"color": COLORS["accent"], "marginTop": 0}),
            html.P(f"{len(edges_df)} edges found. Sorted by edge magnitude.", style={"color": COLORS["text_muted"], "fontSize": "0.85em"}),
            edges_table,
        ]),
        html.Div(style=CARD_STYLE, children=[
            html.H3("Series Overround", style={"color": COLORS["accent"], "marginTop": 0}),
            html.P("Negative = potential arbitrage. Positive = normal vig.", style={"color": COLORS["text_muted"], "fontSize": "0.85em"}),
            sum_chart or html.Div("No data"),
        ]),
    ])


def render_consensus():
    consensus_df = db.get_latest_consensus()
    latest = db.get_latest_odds()

    if consensus_df is None or consensus_df.empty:
        return html.Div(style=CARD_STYLE, children=[
            html.H3("Consensus Comparison", style={"color": COLORS["accent"], "marginTop": 0}),
            html.P("No consensus data yet. Click Refresh to scrape mock draft data.", style={"color": COLORS["text_muted"]}),
        ])

    # Build comparison: consensus rank vs Kalshi #1 pick probability
    draft1 = latest[latest["series_ticker"] == "KXNFLDRAFT1"] if latest is not None else pd.DataFrame()

    if not draft1.empty:
        # Fuzzy join on player name
        from consensus import fuzzy_match_name
        consensus_df["kalshi_prob"] = consensus_df["player_name"].apply(
            lambda name: _find_kalshi_prob(name, draft1)
        )
        comparison = consensus_df[consensus_df["kalshi_prob"].notna()].copy()

        if not comparison.empty:
            # Scatter: rank vs probability
            fig = px.scatter(
                comparison, x="rank", y="kalshi_prob",
                text="player_name", size_max=12,
                labels={"rank": "Consensus Board Rank", "kalshi_prob": "Kalshi #1 Pick Prob (%)"},
            )
            fig.update_traces(textposition="top center", marker=dict(size=10, color=COLORS["accent"]))
            fig.update_layout(
                template="plotly_dark",
                paper_bgcolor=COLORS["card"],
                plot_bgcolor=COLORS["card"],
                height=450,
                margin=dict(t=40, b=40),
            )
            scatter = dcc.Graph(figure=fig)
        else:
            scatter = html.Div("No matched players between consensus and Kalshi.", style={"color": COLORS["text_muted"]})
    else:
        scatter = html.Div("No #1 pick series data to compare.", style={"color": COLORS["text_muted"]})
        comparison = pd.DataFrame()

    # Consensus board table
    board_table = dash_table.DataTable(
        data=consensus_df.to_dict("records"),
        columns=[
            {"name": "Rank", "id": "rank", "type": "numeric"},
            {"name": "Player", "id": "player_name"},
            {"name": "Pos", "id": "position"},
            {"name": "School", "id": "school"},
            {"name": "Kalshi #1 %", "id": "kalshi_prob", "type": "numeric"},
        ],
        sort_action="native",
        page_size=30,
        style_header=TABLE_STYLE_HEADER,
        style_data=TABLE_STYLE_DATA,
        style_data_conditional=TABLE_STYLE_DATA_CONDITIONAL,
        style_table={"overflowX": "auto"},
    )

    return html.Div([
        html.Div(style=CARD_STYLE, children=[
            html.H3("Consensus vs Market", style={"color": COLORS["accent"], "marginTop": 0}),
            html.P("Consensus big board rank vs Kalshi implied probability.", style={"color": COLORS["text_muted"], "fontSize": "0.85em"}),
            scatter,
        ]),
        html.Div(style=CARD_STYLE, children=[
            html.H3("Consensus Big Board", style={"color": COLORS["accent"], "marginTop": 0}),
            board_table,
        ]),
    ])


def _find_kalshi_prob(name, draft1_df):
    """Find Kalshi probability for a consensus player name."""
    try:
        from consensus import fuzzy_match_name
        for _, row in draft1_df.iterrows():
            if fuzzy_match_name(name, row["candidate"]):
                return row["last_price"]
    except Exception:
        pass
    return None


def render_portfolio():
    positions, orders = db.get_portfolio()
    changes = db.get_position_changes()

    sections = []

    # Positions
    if positions is not None and not positions.empty:
        pos_df = positions.copy()
        pos_df["side"] = pos_df["position"].apply(lambda x: "YES" if x > 0 else "NO")
        pos_df["contracts"] = pos_df["position"].abs()

        sections.append(html.Div(style=CARD_STYLE, children=[
            html.H3(f"Open Positions ({len(pos_df)})", style={"color": COLORS["accent"], "marginTop": 0}),
            dash_table.DataTable(
                data=pos_df.to_dict("records"),
                columns=[
                    {"name": "Market", "id": "market_name"},
                    {"name": "Contracts", "id": "contracts", "type": "numeric"},
                    {"name": "Side", "id": "side"},
                    {"name": "Exposure", "id": "market_exposure", "type": "numeric", "format": dash_table.FormatTemplate.money(2)},
                    {"name": "P&L", "id": "realized_pnl", "type": "numeric", "format": dash_table.FormatTemplate.money(2)},
                ],
                sort_action="native",
                style_header=TABLE_STYLE_HEADER,
                style_data=TABLE_STYLE_DATA,
                style_data_conditional=TABLE_STYLE_DATA_CONDITIONAL + [
                    {"if": {"filter_query": "{realized_pnl} > 0", "column_id": "realized_pnl"},
                     "color": COLORS["green"], "fontWeight": "bold"},
                    {"if": {"filter_query": "{realized_pnl} < 0", "column_id": "realized_pnl"},
                     "color": COLORS["red"], "fontWeight": "bold"},
                ],
            ),
        ]))
    else:
        sections.append(html.Div(style=CARD_STYLE, children=[
            html.H3("Open Positions", style={"color": COLORS["accent"], "marginTop": 0}),
            html.P("No open positions.", style={"color": COLORS["text_muted"]}),
        ]))

    # Orders
    if orders is not None and not orders.empty:
        sections.append(html.Div(style=CARD_STYLE, children=[
            html.H3(f"Resting Orders ({len(orders)})", style={"color": COLORS["accent"], "marginTop": 0}),
            dash_table.DataTable(
                data=orders.to_dict("records"),
                columns=[
                    {"name": "Market", "id": "market_name"},
                    {"name": "Side", "id": "side"},
                    {"name": "Price", "id": "yes_price", "type": "numeric"},
                    {"name": "Contracts", "id": "remaining_count", "type": "numeric"},
                ],
                sort_action="native",
                style_header=TABLE_STYLE_HEADER,
                style_data=TABLE_STYLE_DATA,
                style_data_conditional=TABLE_STYLE_DATA_CONDITIONAL,
            ),
        ]))

    # Changes
    if changes is not None and not changes.empty:
        sections.append(html.Div(style=CARD_STYLE, children=[
            html.H3("Changes Since Last Run", style={"color": COLORS["accent"], "marginTop": 0}),
            dash_table.DataTable(
                data=changes.to_dict("records"),
                columns=[
                    {"name": "Market", "id": "market_name"},
                    {"name": "Type", "id": "change_type"},
                    {"name": "Pos Change", "id": "pos_change", "type": "numeric"},
                    {"name": "Current", "id": "current_pos", "type": "numeric"},
                ],
                sort_action="native",
                style_header=TABLE_STYLE_HEADER,
                style_data=TABLE_STYLE_DATA,
                style_data_conditional=TABLE_STYLE_DATA_CONDITIONAL,
            ),
        ]))

    if not sections:
        return html.Div(style=CARD_STYLE, children=[
            html.H3("Portfolio", style={"color": COLORS["accent"], "marginTop": 0}),
            html.P("No portfolio data. Set up .env with Kalshi credentials.", style={"color": COLORS["text_muted"]}),
        ])

    return html.Div(sections)


# ---------------------------------------------------------------------------
# Portal tab renderers (Tasks 22-23)
# ---------------------------------------------------------------------------

# All 6 venues — order of columns in the cross-book grid
VENUES = ["kalshi", "draftkings", "fanduel", "bookmaker", "wagerzon", "hoop88"]


def render_crossbook_grid():
    """Markets on rows, venues on columns. Threshold slider controls outlier flags.

    Cheap-poll guard: the callback that updates this tab only re-renders when
    MAX(fetched_at) in draft_odds changes OR the user moved the threshold slider.
    """
    return html.Div([
        html.Div(style=CARD_STYLE, children=[
            html.H3("Cross-Book Grid", style={"color": COLORS["accent"], "marginTop": 0}),
            html.P(
                "Devigged probability per venue. Flagged cells (⚑) differ from the "
                "cross-venue median by at least the threshold.",
                style={"color": COLORS["text_muted"], "fontSize": "0.85em"},
            ),
            html.Div([
                html.Label("Outlier threshold (percentage points):",
                           style={"color": COLORS["text_muted"], "marginRight": "12px"}),
                dcc.Slider(
                    id="crossbook-threshold",
                    min=0, max=25, step=1, value=10,
                    marks={i: str(i) for i in range(0, 26, 5)},
                ),
            ], style={"marginBottom": "16px"}),
            html.Div(id="crossbook-table-wrap"),
        ]),
    ])


def render_ev_candidates():
    """Flat list of flagged (market, venue) outliers, sorted by |delta| desc."""
    return html.Div([
        html.Div(style=CARD_STYLE, children=[
            html.H3("+EV Candidates", style={"color": COLORS["accent"], "marginTop": 0}),
            html.P(
                "One row per outlier. Negative delta = venue is lower than consensus "
                "(bet YES); positive delta = venue is higher (bet NO).",
                style={"color": COLORS["text_muted"], "fontSize": "0.85em"},
            ),
            html.Div([
                html.Label("Threshold (pp):",
                           style={"color": COLORS["text_muted"], "marginRight": "12px"}),
                dcc.Slider(
                    id="ev-threshold",
                    min=0, max=25, step=1, value=10,
                    marks={i: str(i) for i in range(0, 26, 5)},
                ),
            ], style={"marginBottom": "16px"}),
            html.Div(id="ev-table-wrap"),
        ]),
    ])


def render_trade_tape():
    """Recent Kalshi trades with large-fill highlighting."""
    return html.Div([
        html.Div(style=CARD_STYLE, children=[
            html.H3("Trade Tape", style={"color": COLORS["accent"], "marginTop": 0}),
            html.P(
                "Most-recent 200 Kalshi trades. Large fills highlighted.",
                style={"color": COLORS["text_muted"], "fontSize": "0.85em"},
            ),
            html.Div([
                html.Label("Large fill threshold (USD):",
                           style={"color": COLORS["text_muted"], "marginRight": "12px"}),
                dcc.Slider(
                    id="tape-threshold",
                    min=0, max=5000, step=100, value=500,
                    marks={i: f"${i}" for i in range(0, 5001, 1000)},
                ),
            ], style={"marginBottom": "16px"}),
            dcc.Input(
                id="tape-ticker-filter",
                placeholder="Filter by ticker substring (optional)...",
                style={
                    "width": "100%", "padding": "8px", "marginBottom": "12px",
                    "backgroundColor": "#0d1b2a", "color": COLORS["text"],
                    "border": f"1px solid {COLORS['card_border']}", "borderRadius": "6px",
                },
            ),
            html.Div(id="tape-table-wrap"),
        ]),
    ])


def render_bet_log():
    """Bet entry form + table of past bets.

    Market dropdown is searchable so the user can type any market_id and
    narrow it down. The form reads nfl_draft.bet_log_prefill store on render
    to pre-fill fields when the user clicked 'Log this bet' on +EV.
    """
    market_ids = nfl_queries.all_market_ids()
    if isinstance(market_ids, QueryLocked):
        # Lock contention during a tab-switch render: show an empty dropdown;
        # the interval tick on the new tab will repopulate on next render.
        market_ids = []
    return html.Div([
        html.Div(style=CARD_STYLE, children=[
            html.H3("Log a Bet", style={"color": COLORS["accent"], "marginTop": 0}),
            html.Div([
                html.Label("Market", style={"color": COLORS["text_muted"]}),
                dcc.Dropdown(
                    id="betlog-market",
                    searchable=True,
                    options=[{"label": m, "value": m} for m in market_ids],
                    placeholder="Start typing a market_id...",
                    style={"backgroundColor": "#0d1b2a", "color": COLORS["text"],
                           "marginBottom": "10px"},
                ),
                html.Label("Book", style={"color": COLORS["text_muted"]}),
                dcc.Dropdown(
                    id="betlog-book",
                    options=[{"label": v, "value": v} for v in VENUES],
                    style={"backgroundColor": "#0d1b2a", "color": COLORS["text"],
                           "marginBottom": "10px"},
                ),
                html.Label("Side (yes/no)", style={"color": COLORS["text_muted"]}),
                dcc.Dropdown(
                    id="betlog-side",
                    options=[{"label": "YES", "value": "yes"}, {"label": "NO", "value": "no"}],
                    value="yes",
                    style={"backgroundColor": "#0d1b2a", "color": COLORS["text"],
                           "marginBottom": "10px"},
                ),
                html.Label("American Odds", style={"color": COLORS["text_muted"]}),
                dcc.Input(
                    id="betlog-odds",
                    type="number",
                    placeholder="e.g. -110 or +150",
                    style={"width": "100%", "padding": "8px", "marginBottom": "10px",
                           "backgroundColor": "#0d1b2a", "color": COLORS["text"],
                           "border": f"1px solid {COLORS['card_border']}", "borderRadius": "6px"},
                ),
                html.Label("Stake (USD)", style={"color": COLORS["text_muted"]}),
                dcc.Input(
                    id="betlog-stake",
                    type="number",
                    placeholder="50",
                    style={"width": "100%", "padding": "8px", "marginBottom": "10px",
                           "backgroundColor": "#0d1b2a", "color": COLORS["text"],
                           "border": f"1px solid {COLORS['card_border']}", "borderRadius": "6px"},
                ),
                html.Label("Note", style={"color": COLORS["text_muted"]}),
                dcc.Textarea(
                    id="betlog-note",
                    placeholder="Rationale, stale line on book X, etc.",
                    style={"width": "100%", "padding": "8px", "marginBottom": "10px",
                           "backgroundColor": "#0d1b2a", "color": COLORS["text"],
                           "border": f"1px solid {COLORS['card_border']}", "borderRadius": "6px",
                           "minHeight": "60px"},
                ),
                html.Button("Log Bet", id="betlog-submit", n_clicks=0, style={
                    "background": COLORS["accent"], "color": "white",
                    "padding": "10px 20px", "borderRadius": "8px", "border": "none",
                    "fontWeight": "600", "cursor": "pointer",
                }),
                html.Div(id="betlog-submit-msg", style={"marginTop": "10px", "color": COLORS["green"]}),
            ]),
        ]),
        html.Div(style=CARD_STYLE, children=[
            html.H3("Past Bets", style={"color": COLORS["accent"], "marginTop": 0}),
            html.Div(id="betlog-table-wrap"),
        ]),
    ])


# ---------------------------------------------------------------------------
# Callbacks
# ---------------------------------------------------------------------------


@callback(Output("section-content", "children"), Input("section-tabs", "value"))
def render_section(section):
    """Outer section picker: which group of tabs to show."""
    if section == "portal":
        return _portal_tabs()
    elif section == "legacy":
        return _legacy_tabs()
    return html.Div("Unknown section")


@callback(Output("tab-content", "children"),
          Input("portal-tabs", "value"),
          prevent_initial_call=False)
def render_portal_tab(portal_tab):
    """Render the active Portal-section inner tab. Fires only when portal-tabs
    mounts (after section-content is populated with _portal_tabs())."""
    if portal_tab == "crossbook":
        return render_crossbook_grid()
    elif portal_tab == "ev":
        return render_ev_candidates()
    elif portal_tab == "tape":
        return render_trade_tape()
    elif portal_tab == "betlog":
        return render_bet_log()
    return html.Div("Unknown Portal tab")


@callback(Output("tab-content", "children", allow_duplicate=True),
          Input("main-tabs", "value"),
          prevent_initial_call=True)
def render_legacy_tab(legacy_tab):
    """Render the active Kalshi-legacy inner tab. Fires only when main-tabs
    mounts (after section-content is populated with _legacy_tabs())."""
    if legacy_tab == "overview":
        return render_overview()
    elif legacy_tab == "history":
        return render_history()
    elif legacy_tab == "edges":
        return render_edges()
    elif legacy_tab == "consensus":
        return render_consensus()
    elif legacy_tab == "portfolio":
        return render_portfolio()
    return html.Div("Unknown tab")


# --- Mode toggle: sets interval ms and persists to Store ---
@callback(
    Output("nfl_draft__interval", "interval"),
    Output("nfl_draft__mode_toggle", "data"),
    Input("nfl_draft__mode_toggle_radio", "value"),
)
def _apply_mode(mode):
    # Pre-draft (60s) vs Draft-day (15s) — just controls auto-poll cadence
    return (15_000 if mode == "draft" else 60_000), mode


# --- Cross-Book Grid callback with cheap-poll guard ---
@callback(
    Output("crossbook-table-wrap", "children"),
    Output("nfl_draft__last_fetched_odds", "data"),
    Input("crossbook-threshold", "value"),
    Input("nfl_draft__interval", "n_intervals"),
    State("nfl_draft__last_fetched_odds", "data"),
    prevent_initial_call=True,
)
def _update_crossbook(threshold_pp, _n_intervals, last_seen):
    """Cheap-poll guard: on interval tick, skip render if MAX(fetched_at)
    hasn't changed. User-triggered threshold changes always re-render.

    If DuckDB is locked by the cron writer, skip this render entirely —
    the next interval tick (seconds away) will retry cleanly.
    """
    ctx = dash.callback_context
    triggered_by_interval = (
        ctx.triggered and ctx.triggered[0]["prop_id"].startswith("nfl_draft__interval")
    )
    latest = nfl_queries.latest_max_fetched_at("draft_odds")
    if isinstance(latest, QueryLocked):
        raise PreventUpdate
    latest_iso = latest.isoformat() if latest else None
    if triggered_by_interval and latest_iso == last_seen:
        # Nothing new on the wire since last poll — don't waste cycles
        raise PreventUpdate

    grid = nfl_queries.cross_book_grid(threshold_pp=threshold_pp or 0)
    if isinstance(grid, QueryLocked):
        raise PreventUpdate
    rows = []
    for m in grid:
        row = {"market_id": m["market_id"]}
        for venue in VENUES:
            prob = m["books"].get(venue)
            flagged = m["flags"].get(venue, False)
            if prob is None:
                row[venue] = ""
            else:
                row[venue] = f"{prob*100:.1f}%" + (" \u2691" if flagged else "")
        row["median"] = f"{m['median']*100:.1f}%" if m["median"] is not None else ""
        row["outliers"] = m["outlier_count"]
        rows.append(row)

    table = dash_table.DataTable(
        data=rows,
        columns=[{"name": "Market", "id": "market_id"}]
        + [{"name": v.capitalize(), "id": v} for v in VENUES]
        + [{"name": "Median", "id": "median"}, {"name": "Outliers", "id": "outliers", "type": "numeric"}],
        filter_action="native",
        sort_action="native",
        page_size=50,
        style_header=TABLE_STYLE_HEADER,
        style_data=TABLE_STYLE_DATA,
        style_data_conditional=TABLE_STYLE_DATA_CONDITIONAL + [
            # Highlight any cell containing the flag glyph
            {"if": {"filter_query": f"{{{v}}} contains '\u2691'", "column_id": v},
             "backgroundColor": "#4d2d15", "color": COLORS["red"], "fontWeight": "bold"}
            for v in VENUES
        ],
        style_table={"overflowX": "auto"},
        style_filter={"backgroundColor": "#0d1b2a", "color": COLORS["text"]},
    )
    return table, latest_iso


# --- +EV Candidates callback ---
@callback(
    Output("ev-table-wrap", "children"),
    Input("ev-threshold", "value"),
    Input("nfl_draft__interval", "n_intervals"),
)
def _update_ev(threshold_pp, _n_intervals):
    rows_raw = nfl_queries.ev_candidates(threshold_pp=threshold_pp or 0)
    if isinstance(rows_raw, QueryLocked):
        raise PreventUpdate
    rows = []
    for r in rows_raw:
        rows.append({
            "market_id": r["market_id"],
            "book": r["book"],
            "book_prob": f"{r['book_prob']*100:.1f}%",
            "median": f"{r['median']*100:.1f}%",
            "delta": f"{r['delta']*100:+.1f}pp",
            "direction": "Bet YES (book low)" if r["delta"] < 0 else "Bet NO (book high)",
            # Hidden raw fields used by the log-this-bet action
            "_delta_raw": r["delta"],
        })

    table = dash_table.DataTable(
        id="ev-table",
        data=rows,
        columns=[
            {"name": "Market", "id": "market_id"},
            {"name": "Book", "id": "book"},
            {"name": "Book Prob", "id": "book_prob"},
            {"name": "Median", "id": "median"},
            {"name": "Delta", "id": "delta"},
            {"name": "Direction", "id": "direction"},
        ],
        row_selectable="single",
        selected_rows=[],
        filter_action="native",
        sort_action="native",
        page_size=50,
        style_header=TABLE_STYLE_HEADER,
        style_data=TABLE_STYLE_DATA,
        style_data_conditional=TABLE_STYLE_DATA_CONDITIONAL,
        style_table={"overflowX": "auto"},
        style_filter={"backgroundColor": "#0d1b2a", "color": COLORS["text"]},
    )
    hint = html.P(
        "Select a row to pre-fill the Bet Log with this (market, book). Switch to the Bet Log tab to confirm.",
        style={"color": COLORS["text_muted"], "fontSize": "0.85em", "marginTop": "8px"},
    )
    return html.Div([table, hint])


# --- "Log this bet" hand-off: EV-candidate row selection -> prefill store ---
@callback(
    Output("nfl_draft__bet_log_prefill", "data"),
    Input("ev-table", "selected_rows"),
    State("ev-table", "data"),
    prevent_initial_call=True,
)
def _ev_to_prefill(selected_rows, data):
    if not selected_rows or not data:
        raise PreventUpdate
    row = data[selected_rows[0]]
    return {"market_id": row["market_id"], "book": row["book"]}


# --- Trade Tape callback ---
@callback(
    Output("tape-table-wrap", "children"),
    Output("nfl_draft__last_fetched_trades", "data"),
    Input("tape-threshold", "value"),
    Input("tape-ticker-filter", "value"),
    Input("nfl_draft__interval", "n_intervals"),
    State("nfl_draft__last_fetched_trades", "data"),
    prevent_initial_call=True,
)
def _update_tape(threshold_usd, ticker_filter, _n_intervals, last_seen):
    ctx = dash.callback_context
    triggered_by_interval = (
        ctx.triggered and ctx.triggered[0]["prop_id"].startswith("nfl_draft__interval")
    )
    latest = nfl_queries.latest_max_fetched_at("kalshi_trades")
    if isinstance(latest, QueryLocked):
        raise PreventUpdate
    latest_iso = latest.isoformat() if latest else None
    if triggered_by_interval and latest_iso == last_seen:
        raise PreventUpdate

    trades = nfl_queries.trade_tape(limit=200, large_threshold_usd=threshold_usd or 0)
    if isinstance(trades, QueryLocked):
        raise PreventUpdate
    if ticker_filter:
        needle = ticker_filter.strip().upper()
        trades = [t for t in trades if needle in (t.get("ticker") or "").upper()]

    # Convert to display-friendly rows
    for t in trades:
        t["notional_usd"] = round(t["notional_usd"], 2) if t["notional_usd"] else 0
        t["traded_at"] = str(t["traded_at"])

    table = dash_table.DataTable(
        data=trades,
        columns=[
            {"name": "Traded At", "id": "traded_at"},
            {"name": "Ticker", "id": "ticker"},
            {"name": "Side", "id": "side"},
            {"name": "Price (¢)", "id": "price_cents", "type": "numeric"},
            {"name": "Count", "id": "count", "type": "numeric"},
            {"name": "Notional", "id": "notional_usd", "type": "numeric",
             "format": dash_table.FormatTemplate.money(2)},
        ],
        filter_action="native",
        sort_action="native",
        page_size=50,
        style_header=TABLE_STYLE_HEADER,
        style_data=TABLE_STYLE_DATA,
        style_data_conditional=TABLE_STYLE_DATA_CONDITIONAL + [
            {"if": {"filter_query": "{is_large} = True"},
             "backgroundColor": "#1a3b2a", "fontWeight": "600"},
        ],
        style_table={"overflowX": "auto"},
        style_filter={"backgroundColor": "#0d1b2a", "color": COLORS["text"]},
    )
    return table, latest_iso


# --- Bet Log: apply prefill and render past bets ---
@callback(
    Output("betlog-market", "value"),
    Output("betlog-book", "value"),
    Input("nfl_draft__bet_log_prefill", "data"),
)
def _apply_prefill(prefill):
    if not prefill:
        raise PreventUpdate
    return prefill.get("market_id"), prefill.get("book")


@callback(
    Output("betlog-submit-msg", "children"),
    Output("betlog-table-wrap", "children"),
    Input("betlog-submit", "n_clicks"),
    Input("nfl_draft__interval", "n_intervals"),
    State("betlog-market", "value"),
    State("betlog-book", "value"),
    State("betlog-side", "value"),
    State("betlog-odds", "value"),
    State("betlog-stake", "value"),
    State("betlog-note", "value"),
)
def _log_bet_and_render(n_clicks, _n_intervals, market_id, book, side, odds, stake, note):
    """Insert a new row on submit; always re-render the history table."""
    ctx = dash.callback_context
    msg = ""
    if ctx.triggered and ctx.triggered[0]["prop_id"].startswith("betlog-submit"):
        if market_id and book and odds is not None and stake is not None:
            try:
                with nfl_db.write_connection() as con:
                    con.execute(
                        "INSERT INTO draft_bets VALUES (?, ?, ?, ?, ?, ?, ?, ?)",
                        [str(uuid.uuid4()), market_id, book, side or "yes",
                         int(odds), float(stake), datetime.now(), note or ""],
                    )
                msg = f"Logged bet on {market_id} at {book} ({odds})."
            except Exception as e:  # noqa: BLE001 — surface any DB error to the UI
                msg = f"Error: {e}"
        else:
            msg = "Missing required fields (market, book, odds, stake)."

    rows = nfl_queries.bet_log_rows()
    if isinstance(rows, QueryLocked):
        # Cron writer holds the lock — let the last rendered table stand; the
        # next interval tick (seconds away) will refresh cleanly.
        raise PreventUpdate
    for r in rows:
        r["taken_at"] = str(r["taken_at"])
    table = dash_table.DataTable(
        data=rows,
        columns=[
            {"name": "Taken At", "id": "taken_at"},
            {"name": "Market", "id": "market_id"},
            {"name": "Book", "id": "book"},
            {"name": "Side", "id": "side"},
            {"name": "Odds", "id": "american_odds", "type": "numeric"},
            {"name": "Stake", "id": "stake_usd", "type": "numeric",
             "format": dash_table.FormatTemplate.money(2)},
            {"name": "Note", "id": "note"},
        ],
        sort_action="native",
        page_size=25,
        style_header=TABLE_STYLE_HEADER,
        style_data=TABLE_STYLE_DATA,
        style_data_conditional=TABLE_STYLE_DATA_CONDITIONAL,
        style_table={"overflowX": "auto"},
    )
    return msg, table


@callback(
    Output("history-chart", "figure"),
    Input("history-selector", "value"),
    prevent_initial_call=True,
)
def update_history_chart(selected_tickers):
    if not selected_tickers:
        fig = go.Figure()
        fig.update_layout(
            template="plotly_dark",
            paper_bgcolor=COLORS["card"],
            plot_bgcolor=COLORS["card"],
            annotations=[dict(text="Select markets above", showarrow=False, font=dict(size=16, color=COLORS["text_muted"]))],
        )
        return fig

    history = db.get_price_history(tickers=selected_tickers, days=90)
    if history is None or history.empty:
        fig = go.Figure()
        fig.update_layout(
            template="plotly_dark",
            paper_bgcolor=COLORS["card"],
            plot_bgcolor=COLORS["card"],
            annotations=[dict(text="Not enough history yet. Run fetcher multiple times.", showarrow=False, font=dict(size=14, color=COLORS["text_muted"]))],
        )
        return fig

    fig = go.Figure()
    for ticker in selected_tickers:
        t_data = history[history["ticker"] == ticker]
        if t_data.empty:
            continue
        label = t_data.iloc[0]["candidate"]
        fig.add_scatter(
            x=t_data["fetch_time"], y=t_data["last_price"],
            name=label, mode="lines+markers",
            hovertemplate=f"<b>{label}</b><br>%{{x}}<br>%{{y}}%<extra></extra>",
        )

    fig.update_layout(
        template="plotly_dark",
        paper_bgcolor=COLORS["card"],
        plot_bgcolor=COLORS["card"],
        yaxis_title="Probability (%)",
        xaxis_title="Time",
        height=500,
        margin=dict(t=20, b=40),
        legend=dict(orientation="h", yanchor="bottom", y=1.02, xanchor="right", x=1),
    )
    return fig


@callback(
    Output("refresh-output", "children"),
    Output("toast-area", "children"),
    Input("refresh-btn", "n_clicks"),
    prevent_initial_call=True,
)
def refresh_data(n_clicks):
    """Run the full fetch + edge detection + consensus pipeline."""
    try:
        venv_python = str(Path(__file__).parent / "venv" / "bin" / "python")
        script_dir = str(Path(__file__).parent)

        # Run fetcher
        subprocess.run(
            [venv_python, "fetcher.py"],
            cwd=script_dir, capture_output=True, text=True, timeout=120,
        )

        # Run edge detector
        subprocess.run(
            [venv_python, "edge_detector.py"],
            cwd=script_dir, capture_output=True, text=True, timeout=60,
        )

        # Run consensus scraper
        subprocess.run(
            [venv_python, "consensus.py"],
            cwd=script_dir, capture_output=True, text=True, timeout=60,
        )

        toast = html.Div(
            "Data refreshed successfully!",
            style={
                "background": "linear-gradient(135deg, #00b894, #00cec9)",
                "color": "white", "padding": "16px 24px", "borderRadius": "10px",
                "fontWeight": "500", "boxShadow": "0 8px 30px rgba(0,0,0,0.3)",
            },
        )
        return "", toast

    except Exception as e:
        toast = html.Div(
            f"Refresh failed: {str(e)}",
            style={
                "background": "linear-gradient(135deg, #e74c3c, #c0392b)",
                "color": "white", "padding": "16px 24px", "borderRadius": "10px",
                "fontWeight": "500",
            },
        )
        return "", toast


# ---------------------------------------------------------------------------
# Run
# ---------------------------------------------------------------------------
if __name__ == "__main__":
    import os
    nfl_db.init_schema()
    port = int(os.environ.get("NFL_DRAFT_DASHBOARD_PORT", "8090"))
    app.run(debug=False, host="127.0.0.1", port=port)
