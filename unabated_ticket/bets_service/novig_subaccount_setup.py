"""Setup of a Novig subaccount and its trading key: the automated taker's
("bot", Cal 2026-10-16: $500 of its own) or the Ticket's Novig Place button's
("place", funded when Cal wants). Plan: automated_taker_plan.md § Novig.

A subaccount is a separate Novig wallet; its one trading key places and
cancels there and can never move money. Only Cal's management key opens,
funds and lists subaccounts, and it never leaves his Mac. So the setup runs
in two halves, and the trading key's private half never moves:

  1. on the VM, inside the bets container (`keygen`): makes the trading
     keypair, writes the private half to the name's key file at mode 0600
     (an existing usable key is kept), and prints ONLY its public half on
     stdout. For a NEW "bot" key it first puts the Novig stop file
     (NOVIG_TAKER_STOP_PATH) in place, so the restart that arms the desk
     cannot trade until Cal removes it (a kept key, a top-up, leaves it be);
  2. on the Mac (`open`), signed by the management key in ~/.novig: reads that
     public half on stdin, opens the subaccount labelled "unabated-<name>"
     with it (or finds the one already open), checks Novig's record of the
     trading key (scope "trading", that label, the same fingerprint), then
     moves `--fund` dollars from the app's cash wallet into it.

One Mac line runs both (deploy/README.md step 13):

    ssh ... 'docker compose exec -T bets python -m unabated_ticket.bets_service.novig_subaccount_setup keygen --name bot' \
      | python -m unabated_ticket.bets_service.novig_subaccount_setup open --name bot --fund 500

Funding is safe to re-run: the transfer's clientTransferId is
"unabated-<name>-<run id>" (--run-id, "initial" by default), and Novig answers
a repeat id with the first request's reply and moves nothing. A top-up is the
same line with a new --run-id. The transfer is polled to Applied or Rejected.

Side effects: keygen writes the key file and (bot) the stop file. open makes,
on Novig, at most one subaccount with its trading key and one transfer per run
id. Neither reads nor writes any DuckDB, and no private key is ever printed.
"""
import argparse
import hashlib
import os
import sys
import time
from collections.abc import Callable
from pathlib import Path
from typing import TextIO

from unabated_ticket.bets_service import config
from unabated_ticket.bets_service.novig_orders import KEYS_PATH, NovigError, NovigOrderClient, public_key_pem_of
from unabated_ticket.bets_service.sources.novig_subaccount import NAMES

LABEL_PREFIX = "unabated-"
SUBACCOUNTS_PATH = "/v3/account/subaccounts"
TRADING_SCOPE = "trading"
DEFAULT_RUN_ID = "initial"
MAX_TRANSFER_ID_CHARS = 64
TRANSFER_APPLIED = "Applied"
TRANSFER_REJECTED = "Rejected"
TRANSFER_POLL_ATTEMPTS = 10
TRANSFER_POLL_SEC = 2.0
DEFAULT_MANAGEMENT_DIR = Path.home() / ".novig"


class SetupError(RuntimeError):
    """A step Novig refused or the state rules out; the message says which."""


def label_of(name: str) -> str:
    return f"{LABEL_PREFIX}{name}"


def key_path_of(name: str) -> Path:
    """Where the name's trading key lives on the VM (compose sets both paths under /data)."""
    raw = config.NOVIG_TAKER_PRIVATE_KEY_PATH if name == "bot" else config.NOVIG_PLACE_PRIVATE_KEY_PATH
    if not raw:
        variable = "NOVIG_TAKER_PRIVATE_KEY_PATH" if name == "bot" else "NOVIG_PLACE_PRIVATE_KEY_PATH"
        raise SetupError(f"{variable} (where the {name} trading key goes) is not set")
    return Path(raw)


def fingerprint_of(public_key_pem: str) -> str:
    """Novig's key fingerprint: "sha256:" + hex SHA-256 of the SPKI DER."""
    from cryptography.hazmat.primitives import serialization
    public_key = serialization.load_pem_public_key(public_key_pem.encode())
    der = public_key.public_bytes(serialization.Encoding.DER, serialization.PublicFormat.SubjectPublicKeyInfo)
    return f"sha256:{hashlib.sha256(der).hexdigest()}"


def algorithm_of(public_key_pem: str) -> str:
    from cryptography.hazmat.primitives import serialization
    from cryptography.hazmat.primitives.asymmetric import ec, ed25519
    public_key = serialization.load_pem_public_key(public_key_pem.encode())
    if isinstance(public_key, ed25519.Ed25519PublicKey):
        return "Ed25519"
    if isinstance(public_key, ec.EllipticCurvePublicKey) and public_key.curve.name == "secp256r1":
        return "P-256"
    raise SetupError(f"expected an Ed25519 or P-256 public key, got {type(public_key).__name__}")


# ---- half 1: keygen, on the VM ---------------------------------------------------------------

def write_new_private_key(path: Path) -> None:
    """A fresh Ed25519 key at `path`, mode 0600, swapped in whole."""
    from cryptography.hazmat.primitives import serialization
    from cryptography.hazmat.primitives.asymmetric import ed25519
    pem = ed25519.Ed25519PrivateKey.generate().private_bytes(
        serialization.Encoding.PEM, serialization.PrivateFormat.PKCS8, serialization.NoEncryption())
    staging = path.with_name(path.name + ".new")
    try:
        descriptor = os.open(staging, os.O_WRONLY | os.O_CREAT | os.O_TRUNC, 0o600)
        with os.fdopen(descriptor, "wb") as key_file:
            key_file.write(pem)
            key_file.flush()
            os.fsync(key_file.fileno())
        os.replace(staging, path)
    finally:
        staging.unlink(missing_ok=True)


def ensure_novig_stopped(stop_path: Path, err: TextIO) -> None:
    if stop_path.exists():
        err.write(f"Stop file {stop_path} already in place.\n")
        return
    stop_path.parent.mkdir(parents=True, exist_ok=True)
    stop_path.touch()
    err.write(f"Put the Novig stop file in place: {stop_path} (the Novig desk cannot trade until you rm it).\n")


def run_keygen(name: str, key_path: Path, stop_path: Path, out: TextIO, err: TextIO) -> int:
    """Half 1. Prints the public key on `out` and nothing else there."""
    try:
        if key_path.exists():
            try:
                public_pem = public_key_pem_of(key_path.read_bytes())
            except (OSError, ValueError) as error:
                raise SetupError(f"{key_path} exists but is not a usable private key ({type(error).__name__}); "
                                 f"move it aside and run this again") from error
            err.write(f"Kept the existing {name} trading key {key_path}.\n")
        else:
            # A new bot key is a first setup: stop the desk before the restart that arms it. A kept
            # key is a re-run or a top-up, which must not stop a bot already trading.
            if name == "bot":
                ensure_novig_stopped(stop_path, err)
            key_path.parent.mkdir(parents=True, exist_ok=True)
            write_new_private_key(key_path)
            public_pem = public_key_pem_of(key_path.read_bytes())
            err.write(f"Made a new {name} trading key at {key_path} (mode 0600).\n")
    except (SetupError, OSError) as error:
        err.write(f"STOPPED: {error}\n")
        return 1
    out.write(public_pem)
    return 0


# ---- half 2: open and fund, on the Mac --------------------------------------------------------

def _signed(client: NovigOrderClient, method: str, path: str, body: dict | None = None) -> tuple[int, object]:
    try:
        return client.request(method, path, body=body)
    except OSError as error:
        raise SetupError(f"{method} {path}: no answer from Novig ({error})") from error


def find_subaccount(client: NovigOrderClient, label: str) -> dict | None:
    """The live subaccount row with this label ({keyId, label, balance}), or None."""
    status, reply = _signed(client, "GET", SUBACCOUNTS_PATH)
    rows = reply.get("items") if isinstance(reply, dict) else reply
    if status != 200 or not isinstance(rows, list):
        raise SetupError(f"GET {SUBACCOUNTS_PATH}: expected 200 with a list, got HTTP {status} {str(reply)[:300]}")
    matches = [row for row in rows if isinstance(row, dict) and row.get("label") == label]
    if len(matches) > 1:
        raise SetupError(f"{len(matches)} subaccounts are labelled {label}; expected at most one")
    return matches[0] if matches else None


def open_subaccount(client: NovigOrderClient, label: str, public_key_pem: str, err: TextIO) -> str:
    """The trading key id of the labelled subaccount, opened with this public key when absent."""
    existing = find_subaccount(client, label)
    if existing:
        err.write(f"Subaccount {label} is already open (trading key {existing.get('keyId')}).\n")
        return existing["keyId"]
    body = {"label": label, "publicKey": public_key_pem, "algorithm": algorithm_of(public_key_pem)}
    status, reply = _signed(client, "POST", SUBACCOUNTS_PATH, body)
    key_id = reply.get("keyId") if isinstance(reply, dict) else None
    if status != 201 or not isinstance(key_id, str) or not key_id:
        raise SetupError(f"POST {SUBACCOUNTS_PATH} was refused: HTTP {status} {str(reply)[:300]}")
    err.write(f"Opened subaccount {label}; its trading key is {key_id}.\n")
    return key_id


def trading_key_problem(client: NovigOrderClient, key_id: str, label: str, public_key_pem: str) -> str | None:
    """Why Novig's record of the trading key is not the VM's key for this
    subaccount, or None: scope trading, restricted to the label, same fingerprint."""
    status, record = _signed(client, "GET", f"{KEYS_PATH}/{key_id}")
    if status != 200 or not isinstance(record, dict):
        return f"GET {KEYS_PATH}/{key_id}: HTTP {status} {str(record)[:300]}"
    if record.get("scope") != TRADING_SCOPE:
        return f"key {key_id} has scope {record.get('scope')!r}, not {TRADING_SCOPE!r}"
    if record.get("restriction") != label:
        return f"key {key_id} is restricted to {record.get('restriction')!r}, not {label!r}"
    if record.get("fingerprint") != fingerprint_of(public_key_pem):
        return (f"key {key_id} is not the VM's key (fingerprint {record.get('fingerprint')}, the VM's is "
                f"{fingerprint_of(public_key_pem)}): the subaccount was opened with another key")
    return None


def transfer_id_of(name: str, run_id: str) -> str:
    transfer_id = f"{LABEL_PREFIX}{name}-{run_id}"
    if len(transfer_id) > MAX_TRANSFER_ID_CHARS:
        raise SetupError(f"transfer id {transfer_id!r} is over {MAX_TRANSFER_ID_CHARS} characters; shorten --run-id")
    return transfer_id


def fund(client: NovigOrderClient, key_id: str, amount: float, transfer_id: str, err: TextIO,
         sleep: Callable[[float], None]) -> str:
    """Moves `amount` from the cash wallet into the subaccount under
    `transfer_id` and waits for Novig to apply or reject it; returns the status."""
    body = {"direction": "fund", "amount": f"{amount:.5f}", "clientTransferId": transfer_id}
    status, reply = _signed(client, "POST", f"{SUBACCOUNTS_PATH}/{key_id}/transfer", body)
    transfer = reply if isinstance(reply, dict) else {}
    if status not in (200, 202) or not transfer.get("transferId"):
        raise SetupError(f"the transfer was refused: HTTP {status} {str(reply)[:300]}")
    for attempt in range(TRANSFER_POLL_ATTEMPTS):
        if transfer.get("status") in (TRANSFER_APPLIED, TRANSFER_REJECTED):
            break
        if attempt:
            sleep(TRANSFER_POLL_SEC)
        status, reply = _signed(client, "GET", f"{SUBACCOUNTS_PATH}/{key_id}/transfers/{transfer['transferId']}")
        if status == 200 and isinstance(reply, dict):
            transfer = reply
    final = transfer.get("status")
    if final == TRANSFER_REJECTED:
        raise SetupError(f"Novig rejected the transfer {transfer_id}: the cash wallet held {transfer.get('actualBalance')}")
    if final != TRANSFER_APPLIED:
        raise SetupError(f"the transfer {transfer_id} is still {final!r}; do not re-send it under a new id "
                         f"(Novig: contact support if it stays Requested)")
    err.write(f"Moved ${float(transfer['amount']):,.2f} into the subaccount (transfer {transfer_id}, Applied).\n")
    return final


def subaccount_balance(client: NovigOrderClient, key_id: str) -> float:
    status, reply = _signed(client, "GET", f"{SUBACCOUNTS_PATH}/{key_id}/balance")
    balance = reply.get("balance") if isinstance(reply, dict) else None
    if status != 200 or not isinstance(balance, str):
        raise SetupError(f"the subaccount's balance could not be read: HTTP {status} {str(reply)[:300]}")
    return float(balance)


def run_open(client: NovigOrderClient, name: str, public_key_pem: str, amount: float | None, run_id: str,
             out: TextIO, sleep: Callable[[float], None] = time.sleep) -> int:
    """Half 2. Returns the process exit code."""
    label = label_of(name)
    try:
        if "PRIVATE KEY" in public_key_pem:
            raise SetupError("stdin holds a PRIVATE key; pipe in only what keygen prints")
        algorithm_of(public_key_pem)
        if amount is not None and not amount > 0:
            raise SetupError(f"--fund must be a positive number of dollars, got {amount}")
        key_id = open_subaccount(client, label, public_key_pem, out)
        problem = trading_key_problem(client, key_id, label, public_key_pem)
        if problem:
            raise SetupError(problem)
        if amount is not None:
            fund(client, key_id, amount, transfer_id_of(name, run_id), out, sleep)
        balance = subaccount_balance(client, key_id)
    except (SetupError, NovigError, ValueError) as error:
        out.write(f"STOPPED: {error}\n")
        return 1
    variable = "NOVIG_TAKER_KEY_ID" if name == "bot" else "NOVIG_PLACE_KEY_ID"
    out.write(f"\nDone. Subaccount {label} holds ${balance:,.2f}; its trading key {key_id} trades only there "
              f"and cannot move money.\nSet this in unabated_ticket/bets_service/.env on the VM, then restart the bets "
              f"service:\n  {variable}={key_id}\n")
    if name == "bot":
        out.write("If this made a new bot key, the Novig stop file stays until you rm it on the VM.\n")
    return 0


def _management_client(directory: Path) -> NovigOrderClient:
    key_id_path = directory / "novig_mgmt.key_id"
    pem_path = directory / "novig_mgmt.pem"
    try:
        key_id = key_id_path.read_text().strip()
        return NovigOrderClient(key_id, pem_path, config.NOVIG_API_URL)
    except (OSError, ValueError, NovigError) as error:
        raise SetupError(f"the management key in {directory} could not be read ({type(error).__name__})") from error


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description="Set up a Novig subaccount and its trading key.")
    commands = parser.add_subparsers(dest="command", required=True)
    keygen = commands.add_parser("keygen", help="on the VM: make the trading key, print its public half")
    keygen.add_argument("--name", choices=NAMES, required=True)
    opener = commands.add_parser("open", help="on the Mac: open and fund the subaccount (public key on stdin)")
    opener.add_argument("--name", choices=NAMES, required=True)
    opener.add_argument("--fund", type=float, default=None, help="dollars to move in from the cash wallet")
    opener.add_argument("--run-id", default=DEFAULT_RUN_ID, help="names this funding; a new one for each top-up")
    opener.add_argument("--management-dir", type=Path, default=DEFAULT_MANAGEMENT_DIR)
    args = parser.parse_args(argv)
    if args.command == "keygen":
        try:
            key_path = key_path_of(args.name)
        except SetupError as error:
            print(f"STOPPED: {error}", file=sys.stderr)
            return 1
        return run_keygen(args.name, key_path, config.NOVIG_TAKER_STOP_PATH, sys.stdout, sys.stderr)
    try:
        client = _management_client(args.management_dir)
    except SetupError as error:
        print(f"STOPPED: {error}", file=sys.stderr)
        return 1
    return run_open(client, args.name, sys.stdin.read(), args.fund, args.run_id, sys.stdout)


if __name__ == "__main__":
    sys.exit(main())
