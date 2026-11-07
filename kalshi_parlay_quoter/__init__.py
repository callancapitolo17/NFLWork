"""Kalshi parlay RFQ quoter (epic #281). Phase 1 (#283) is dry-run only:
it watches the RFQ feed, prices in-scope parlays through the runner's
/rfq_price.json and logs a decision per RFQ. It never sends a quote."""
