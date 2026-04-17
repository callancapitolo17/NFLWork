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

import dash
from dash import dcc, html, dash_table, callback, Input, Output, State
from dash.dash_table.Format import Format
import plotly.express as px
import plotly.graph_objects as go
import pandas as pd

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


app.layout = html.Div(
    style={
        "backgroundColor": COLORS["bg"],
        "minHeight": "100vh",
        "padding": "20px",
        "fontFamily": "-apple-system, BlinkMacSystemFont, 'Segoe UI', Roboto, sans-serif",
        "color": COLORS["text"],
    },
    children=[
        make_header(),
        dcc.Tabs(
            id="main-tabs",
            value="overview",
            colors={
                "border": COLORS["card_border"],
                "primary": COLORS["accent"],
                "background": COLORS["card"],
            },
            style={"marginBottom": "16px"},
            children=[
                dcc.Tab(label="Market Overview", value="overview",
                        style={"color": COLORS["text_muted"], "backgroundColor": COLORS["card"], "border": "none", "padding": "12px 20px"},
                        selected_style={"color": COLORS["accent"], "backgroundColor": COLORS["bg"], "borderTop": f"2px solid {COLORS['accent']}", "padding": "12px 20px"}),
                dcc.Tab(label="Price History", value="history",
                        style={"color": COLORS["text_muted"], "backgroundColor": COLORS["card"], "border": "none", "padding": "12px 20px"},
                        selected_style={"color": COLORS["accent"], "backgroundColor": COLORS["bg"], "borderTop": f"2px solid {COLORS['accent']}", "padding": "12px 20px"}),
                dcc.Tab(label="Edge Detection", value="edges",
                        style={"color": COLORS["text_muted"], "backgroundColor": COLORS["card"], "border": "none", "padding": "12px 20px"},
                        selected_style={"color": COLORS["accent"], "backgroundColor": COLORS["bg"], "borderTop": f"2px solid {COLORS['accent']}", "padding": "12px 20px"}),
                dcc.Tab(label="Consensus", value="consensus",
                        style={"color": COLORS["text_muted"], "backgroundColor": COLORS["card"], "border": "none", "padding": "12px 20px"},
                        selected_style={"color": COLORS["accent"], "backgroundColor": COLORS["bg"], "borderTop": f"2px solid {COLORS['accent']}", "padding": "12px 20px"}),
                dcc.Tab(label="Portfolio", value="portfolio",
                        style={"color": COLORS["text_muted"], "backgroundColor": COLORS["card"], "border": "none", "padding": "12px 20px"},
                        selected_style={"color": COLORS["accent"], "backgroundColor": COLORS["bg"], "borderTop": f"2px solid {COLORS['accent']}", "padding": "12px 20px"}),
            ],
        ),
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
                height=max(500, len(pivot) * 40 + 120),
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
# Callbacks
# ---------------------------------------------------------------------------

@callback(Output("tab-content", "children"), Input("main-tabs", "value"))
def render_tab(tab):
    if tab == "overview":
        return render_overview()
    elif tab == "history":
        return render_history()
    elif tab == "edges":
        return render_edges()
    elif tab == "consensus":
        return render_consensus()
    elif tab == "portfolio":
        return render_portfolio()
    return html.Div("Unknown tab")


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
    nfl_db.init_schema()
    app.run(debug=False, host="127.0.0.1", port=8083)
