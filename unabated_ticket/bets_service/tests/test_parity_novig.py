"""normalize_novig and extension/novig_bets.js normalizeNovig must agree byte
for byte on the shared fixture (test_parity.py's rule for the Kalshi pair)."""
import copy
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
// Extra rows on stdin (JSON {orders, parlays}) are normalised after the fixture's.
const extra = JSON.parse(fs.readFileSync(0, "utf8") || "{}");
orders.push(...(extra.orders || []));
parlays.push(...(extra.parlays || []));
const sortKeys = (o) => Array.isArray(o) ? o.map(sortKeys) : (o && typeof o === "object") ? Object.fromEntries(Object.keys(o).sort().map((k) => [k, sortKeys(o[k])])) : o;
console.log(JSON.stringify(sortKeys(novig.normalizeNovig({ orders, parlays, readAt: process.argv[3] }))));
"""


def _node_records(orders_extra: list[dict], parlays_extra: list[dict]) -> str:
    node = shutil.which("node")
    if node is None:
        pytest.skip("node not on PATH — parity with extension/novig_bets.js not checked")
    result = subprocess.run([node, "-e", NODE_SCRIPT, str(EXTENSION_DIR), str(NOVIG_FIXTURE_PATH), NOVIG_READ_AT],
                            input=json.dumps({"orders": orders_extra, "parlays": parlays_extra}),
                            capture_output=True, text=True, check=True, cwd=EXTENSION_DIR)
    return result.stdout.strip()


def _python_records(orders: list[dict], parlays: list[dict]) -> tuple[list[dict], str]:
    python_records = normalize_novig(orders, parlays, NOVIG_READ_AT)
    return python_records, json.dumps(sort_keys(python_records), separators=(",", ":"), ensure_ascii=False)


def test_python_records_are_byte_equivalent_to_the_node_normaliser(novig_rows):
    orders, parlays = novig_rows
    python_records, python_json = _python_records(orders, parlays)
    assert len(python_records) == 19
    assert python_json == _node_records([], [])


def test_degraded_venue_ids_stay_byte_equivalent(novig_rows):
    """Ids the venue sends malformed or not at all (#118) must null out the same way in both."""
    orders, parlays = novig_rows
    odd = copy.deepcopy(next(o for o in orders if o["id"] == "o-cfb"))
    odd.update({"id": "o-odd-ids"})
    odd["market"]["id"] = 42
    odd["outcome"]["id"] = ""
    del odd["market"]["event"]["id"]
    odd["market"]["event"]["game"]["id"] = None
    odd["market"]["event"]["game"]["awayTeam"] = {"name": "Chattanooga", "symbol": ""}
    leg_parlay = copy.deepcopy(parlays[0])
    leg_parlay["id"] = "pl-odd-ids"
    leg_parlay["legs"][0]["outcome"]["market"]["event"] = None
    extra_orders, extra_parlays = [odd, {"id": "o-bare"}], [leg_parlay]
    python_records, python_json = _python_records(orders + extra_orders, parlays + extra_parlays)
    assert len(python_records) == 19 + 2 + 2
    assert python_json == _node_records(extra_orders, extra_parlays)
