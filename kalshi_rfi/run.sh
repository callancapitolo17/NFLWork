#!/bin/bash
# Launch the one-sided Kalshi YRFI maker from the repo root.
# Shadow (default): ./kalshi_rfi/run.sh
# Live:             RFI_MODE=live RFI_LIVE_ACK=1 ./kalshi_rfi/run.sh
cd "$(dirname "$0")/.." || exit 1
exec python3 -m kalshi_rfi.main "$@"
