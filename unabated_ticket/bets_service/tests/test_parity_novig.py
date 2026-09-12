"""normalize_novig and extension/novig_bets.js normalizeNovig must agree byte
for byte on the shared fixture (test_parity.py's rule for the Kalshi pair)."""
import json
import shutil
import subprocess
from pathlib import Path

import pytest

from unabated_ticket.bets_service.sources.novig import normalize_novig
from unabated_ticket.bets_service.tests.conftest import NOVIG_FIXTURE_PATH, NOVIG_READ_AT
from unabated_ticket.bets_service.tests.test_parity import sort_keys

EXTENSION_DIR = Path(__file__).parents[2] / "extension"

NODE_SCRIPT = """
const fs = require("fs");
const novig = require(process.argv[1] + "/novig_bets.js");
const fx = JSON.parse(fs.readFileSync(process.argv[2], "utf8"));
const orders = fx.responses[0].data.ActivePortfolioOrders_Query.concat(fx.responses[1].data.SettledPortfolioOrders_Query);
const parlays = fx.responses[2].data.parlay.concat(fx.responses[3].data.parlay);
const sortKeys = (o) => Array.isArray(o) ? o.map(sortKeys) : (o && typeof o === "object") ? Object.fromEntries(Object.keys(o).sort().map((k) => [k, sortKeys(o[k])])) : o;
console.log(JSON.stringify(sortKeys(novig.normalizeNovig({ orders, parlays, readAt: process.argv[3] }))));
"""


def test_python_records_are_byte_equivalent_to_the_node_normaliser(novig_rows):
    node = shutil.which("node")
    if node is None:
        pytest.skip("node not on PATH — parity with extension/novig_bets.js not checked")
    orders, parlays = novig_rows
    python_records = normalize_novig(orders, parlays, NOVIG_READ_AT)
    assert len(python_records) == 19
    python_json = json.dumps(sort_keys(python_records), separators=(",", ":"), ensure_ascii=False)
    result = subprocess.run([node, "-e", NODE_SCRIPT, str(EXTENSION_DIR), str(NOVIG_FIXTURE_PATH), NOVIG_READ_AT],
                            capture_output=True, text=True, check=True, cwd=EXTENSION_DIR)
    assert python_json == result.stdout.strip()
