"""The Python normaliser and extension/bets.js normalizeKalshi must agree byte
for byte on the shared fixture: the node tests pin the JS semantics, this test
pins the port to them. Runs the JS through `node -e`; skipped (loudly) when
node is not on PATH."""
import json
import shutil
import subprocess
from pathlib import Path

import pytest

from unabated_ticket.bets_service.sources.kalshi import normalize_kalshi
from unabated_ticket.bets_service.tests.conftest import FETCHED_AT, FIXTURE_PATH

EXTENSION_DIR = Path(__file__).parents[2] / "extension"

# Prints two JSON documents on one line each: the JS normaliser's records with
# awayKey/homeKey nulled (what the service emits), and resolveTeamKeys() over
# the Python records (what the panel holds after load) next to the JS records.
NODE_SCRIPT = """
const fs = require("fs");
const bets = require(process.argv[1] + "/bets.js");
const fx = JSON.parse(fs.readFileSync(process.argv[2], "utf8"));
const py = JSON.parse(fs.readFileSync(process.argv[3], "utf8"));
const js = bets.normalizeKalshi({ fills: fx.fills, positions: fx.positions, markets: fx.markets, events: fx.events, fetchedAt: process.argv[4] });
const sortKeys = (o) => Array.isArray(o) ? o.map(sortKeys) : (o && typeof o === "object") ? Object.fromEntries(Object.keys(o).sort().map((k) => [k, sortKeys(o[k])])) : o;
const serviceShape = js.map((r) => Object.assign({}, r, { awayKey: null, homeKey: null }));
console.log(JSON.stringify(sortKeys(serviceShape)));
console.log(JSON.stringify(sortKeys(bets.resolveTeamKeys(py))));
console.log(JSON.stringify(sortKeys(js)));
"""


def sort_keys(value):
    if isinstance(value, list):
        return [sort_keys(item) for item in value]
    if isinstance(value, dict):
        return {key: sort_keys(value[key]) for key in sorted(value)}
    return value


def test_python_records_are_byte_equivalent_to_the_node_normaliser(tmp_path):
    node = shutil.which("node")
    if node is None:
        pytest.skip("node not on PATH — parity with extension/bets.js not checked")
    fixture = json.loads(FIXTURE_PATH.read_text())
    python_records = normalize_kalshi(fixture["fills"], fixture["positions"], fixture["markets"],
                                      fixture["events"], FETCHED_AT)
    assert len(python_records) == 12
    python_json = json.dumps(sort_keys(python_records), separators=(",", ":"))
    python_path = tmp_path / "python_records.json"
    python_path.write_text(python_json)
    result = subprocess.run(
        [node, "-e", NODE_SCRIPT, str(EXTENSION_DIR), str(FIXTURE_PATH), str(python_path), FETCHED_AT],
        capture_output=True, text=True, check=True, cwd=EXTENSION_DIR)
    js_service_shape, resolved_python, js_full = result.stdout.strip().split("\n")
    assert python_json == js_service_shape  # byte-equivalent as the service emits them
    assert resolved_python == js_full  # and identical to the JS once the panel resolves keys
