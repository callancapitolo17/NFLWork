"""
Cross-Market Edge Detection for NFL Draft Prediction Markets.
Identifies mispricing, arbitrage, and inconsistencies across series.
"""

from datetime import datetime, timezone
import db


def check_sum_violations(latest_odds):
    """Check if YES prices within an event sum to ~100%.

    Groups by event_ticker (not series_ticker) since series like KXNFLDRAFTTOP
    contain multiple events (top 5, top 10, etc.) that each independently sum to 100.

    Overround > 0 means normal vig (books profit).
    Overround < 0 means potential arbitrage (free money).
    """
    edges = []

    # Only check series where outcomes are mutually exclusive
    # KXNFLDRAFT1 = who goes #1 (only one player), KXNFLDRAFT1ST = which team picks #1
    # KXNFLDRAFTTOP/KXNFLDRAFTPICK = NOT mutually exclusive (multiple players in top 10)
    mutually_exclusive_series = {"KXNFLDRAFT1", "KXNFLDRAFT1ST"}

    for event in latest_odds["event_ticker"].unique():
        if not event:
            continue

        event_data = latest_odds[latest_odds["event_ticker"] == event]
        series = event_data["series_ticker"].iloc[0]

        if series not in mutually_exclusive_series:
            continue

        # Use last_price for implied probability
        total = int(event_data["last_price"].sum())
        overround = total - 100

        # Use midpoint for sharper estimate
        event_data = event_data.copy()
        mid_mask = (event_data["yes_bid"] > 0) & (event_data["yes_ask"] > 0)
        if mid_mask.any():
            event_data.loc[mid_mask, "midpoint"] = (
                event_data.loc[mid_mask, "yes_bid"] + event_data.loc[mid_mask, "yes_ask"]
            ) / 2

        # Only flag meaningful deviations (> 5%) for events with multiple markets
        if len(event_data) < 2:
            continue

        if abs(overround) > 5:
            if overround < -5:
                desc = f"{event} ({series}): sum={total}%, underround by {abs(overround):.1f}% (potential arb)"
                confidence = "high" if overround < -10 else "medium"
            else:
                desc = f"{event} ({series}): sum={total}%, overround by {overround:.1f}%"
                confidence = "low"

            edges.append({
                "edge_type": "sum_violation",
                "description": desc,
                "market_a": event,
                "market_b": f"sum={total}%",
                "price_a": float(total),
                "price_b": 100.0,
                "implied_edge": round(float(overround), 2),
                "confidence": confidence,
            })

    return edges


def check_cross_market_consistency(latest_odds):
    """Check consistency between related draft series.

    Example: If player X is 95% to go #1 (KXNFLDRAFT1) and team Y is 93%
    to pick #1 (KXNFLDRAFT1ST), then player X to team Y should be ~88%.
    """
    edges = []

    draft1 = latest_odds[latest_odds["series_ticker"] == "KXNFLDRAFT1"]
    draft1st = latest_odds[latest_odds["series_ticker"] == "KXNFLDRAFT1ST"]
    drafttop = latest_odds[latest_odds["series_ticker"] == "KXNFLDRAFTTOP"]
    draftpick = latest_odds[latest_odds["series_ticker"] == "KXNFLDRAFTPICK"]

    # Check 1: Top player in DRAFT1 vs top team in DRAFT1ST
    if not draft1.empty and not draft1st.empty:
        top_player = draft1.loc[draft1["last_price"].idxmax()]
        top_team = draft1st.loc[draft1st["last_price"].idxmax()]

        player_prob = top_player["last_price"]
        team_prob = top_team["last_price"]

        # If both are very high, they should be consistent
        # P(player goes #1) * P(team picks #1) <= min(player, team) if independent
        # But they're NOT independent - if the team has the pick, they pick the player
        if player_prob > 50 and team_prob > 50:
            implied_joint = min(player_prob, team_prob)
            gap = abs(player_prob - team_prob)

            if gap > 5:
                edges.append({
                    "edge_type": "cross_market",
                    "description": f"#1 pick: {top_player['candidate']} ({player_prob}%) vs {top_team['candidate']} ({team_prob}%) - gap {gap}%",
                    "market_a": top_player["candidate"],
                    "market_b": top_team["candidate"],
                    "price_a": float(player_prob),
                    "price_b": float(team_prob),
                    "implied_edge": round(gap, 2),
                    "confidence": "medium" if gap > 10 else "low",
                })

    # Check 2: DRAFTTOP consistency with DRAFT1
    # If a player is X% for #1, they should be >= X% for "top N" markets
    # (being #1 implies being in top 5, top 10, etc.)
    if not draft1.empty and not drafttop.empty:
        for _, player in draft1[draft1["last_price"] >= 5].iterrows():
            name = player["candidate"]
            prob_1 = player["last_price"]

            top_matches = drafttop[
                drafttop["candidate"].str.lower() == name.lower()
            ]

            for _, top_row in top_matches.iterrows():
                prob_top = top_row["last_price"]
                title = top_row["market_title"].lower()

                # Only flag if top-N probability is LESS than #1 probability
                # (should be impossible: #1 pick subset of top-N)
                if "top" in title and prob_1 > prob_top + 3:
                    edges.append({
                        "edge_type": "cross_market",
                        "description": f"{name}: {prob_1}% for #1 but only {prob_top}% for {top_row['market_title']}",
                        "market_a": player["ticker"],
                        "market_b": top_row["ticker"],
                        "price_a": float(prob_1),
                        "price_b": float(prob_top),
                        "implied_edge": round(float(prob_1 - prob_top), 2),
                        "confidence": "high",
                    })

    return edges


def check_spread_opportunities(latest_odds):
    """Find markets with wide bid/ask spreads or stale last_price."""
    edges = []

    for _, row in latest_odds.iterrows():
        bid = row["yes_bid"]
        ask = row["yes_ask"]
        last = row["last_price"]
        spread = ask - bid

        if spread <= 0 or bid == 0:
            continue

        midpoint = (bid + ask) / 2.0

        # Wide spread with stale last price
        if spread >= 5 and abs(last - midpoint) >= 3:
            direction = "above" if last > midpoint else "below"
            edges.append({
                "edge_type": "spread_opportunity",
                "description": f"{row['candidate']}: spread={spread}c, last={last}c {direction} mid={midpoint:.0f}c",
                "market_a": row["ticker"],
                "market_b": f"spread={spread}c",
                "price_a": float(last),
                "price_b": float(midpoint),
                "implied_edge": round(abs(last - midpoint), 2),
                "confidence": "low",
            })

    return edges


def detect_all_edges():
    """Run all edge detection checks and save results."""
    latest = db.get_latest_odds()
    if latest is None or latest.empty:
        print("No odds data for edge detection")
        return []

    print("Running edge detection...")
    all_edges = []
    all_edges.extend(check_sum_violations(latest))
    all_edges.extend(check_cross_market_consistency(latest))
    all_edges.extend(check_spread_opportunities(latest))

    # Save to DB
    fetch_time = datetime.now(timezone.utc).strftime("%Y-%m-%d %H:%M:%S")
    con = db.get_connection()

    # Clear old edges for this run
    con.execute("DELETE FROM detected_edges WHERE fetch_time = ?", [fetch_time])

    if all_edges:
        con.executemany("""
            INSERT INTO detected_edges VALUES (?, ?, ?, ?, ?, ?, ?, ?, ?)
        """, [
            (
                fetch_time,
                str(e["edge_type"]), str(e["description"]),
                str(e["market_a"]), str(e["market_b"]),
                float(e["price_a"]), float(e["price_b"]),
                float(e["implied_edge"]), str(e["confidence"])
            )
            for e in all_edges
        ])

    con.close()

    print(f"  Found {len(all_edges)} edges:")
    for e in sorted(all_edges, key=lambda x: abs(x["implied_edge"]), reverse=True)[:5]:
        print(f"    [{e['confidence']}] {e['description']}")

    return all_edges


if __name__ == "__main__":
    db.init_schema()
    detect_all_edges()
