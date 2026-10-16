"""One-time setup of the automated taker's own Kalshi subaccount (Cal, 2026-10-16).

The bot trades only from a Kalshi subaccount holding its money, on a key
locked to that subaccount: the money in it is all the bot can lose, and the
key cannot move money or touch Cal's primary account. Kalshi's website and app
cannot make subaccounts yet; this does, through the API.

Run on the VM inside the bets container (it reads the same config as the
service; /data is where the bot's key file lives):

    docker compose exec bets python -m unabated_ticket.bets_service.taker_subaccount_setup

It signs with Cal's MAIN key (KALSHI_API_KEY_ID / KALSHI_PRIVATE_KEY_PATH),
which must be allowed to make subaccounts, move money and make keys. Steps,
each safe to re-run:
  0. the stop file (TAKER_STOP_PATH) is put in place if absent, so the restart
     that arms the bot cannot trade until Cal removes it; the bot key file's
     directory is checked writable before anything is asked of Kalshi;
  1. the subaccount: --subaccount N, else the account's one numbered
     subaccount, else a new one (POST /portfolio/subaccounts). Kalshi makes
     subaccounts on its Advanced API tier only; on a refusal the tier is asked
     for once (POST /account/api_usage_level/upgrade, granted when one of the
     account's last 100 orders came by API) and the create retried;
  2. the money: BOT_BANKROLL ($500) moved from the primary account, only while
     the subaccount holds $0, and checked to have landed. The transfer id is
     fixed per subaccount, so the setup funds a subaccount once at most; a
     later top-up is Cal's own transfer;
  3. the bot's key: POST /api_keys/generate locked to the subaccount, scopes
     read + write::trade (no transfer). Its private key is written to
     KALSHI_TAKER_PRIVATE_KEY_PATH at mode 0600, replacing the file there,
     and is never printed. Skipped when KALSHI_TAKER_SUBACCOUNT already names
     this subaccount and the configured key passes step 4;
  4. two checks of the lock: Kalshi's own listing of the key (GET /api_keys on
     the main key: `subaccount` = N, scopes exactly read + write::trade, no
     FCM binding), then the check the desk runs before it arms
     (taker_orders.key_lock_problem) on the new key itself. A key that fails
     either is deleted from disk and named for revoking.
It prints the key id and the .env lines to set; it edits no .env.

Side effects: on Kalshi, a subaccount, one transfer of money and one API key;
on disk, the bot's key file and the stop file. Reads and writes no DuckDB.
"""
import argparse
import os
import sys
import time
import uuid
from collections.abc import Callable
from datetime import datetime, timezone
from pathlib import Path
from typing import TextIO

from unabated_ticket.bets_service import config
from unabated_ticket.bets_service.sources.kalshi import PRIMARY_SUBACCOUNT
from unabated_ticket.bets_service.taker_orders import BOT_BANKROLL, KalshiOrderClient, key_lock_problem

SUBACCOUNTS_PATH = "/portfolio/subaccounts"
SUBACCOUNT_BALANCES_PATH = "/portfolio/subaccounts/balances"
TRANSFER_PATH = "/portfolio/subaccounts/transfer"
UPGRADE_PATH = "/account/api_usage_level/upgrade"
GENERATE_KEY_PATH = "/api_keys/generate"
API_KEYS_PATH = "/api_keys"
# Read and trade in its own subaccount; never write::transfer.
BOT_KEY_SCOPES = ("read", "write::trade")
# The setup's transfer id is uuid5 of this and the subaccount: one funding per subaccount, ever.
FUNDING_ID_NAMESPACE = uuid.UUID("6f1d2b7e-3c4a-4e59-9d0b-6a2f5c8e1b37")
# Kalshi reports dollars to the cent; a landed transfer may show a hair under.
FUNDING_TOLERANCE_DOLLARS = 0.01
# A key Kalshi just made may take a moment to be accepted.
NEW_KEY_CHECK_ATTEMPTS = 3
NEW_KEY_CHECK_RETRY_SEC = 3.0


class SetupError(RuntimeError):
    """A step Kalshi refused or the account's state rules out; the message says which."""


def funding_transfer_id(subaccount: int) -> str:
    return str(uuid.uuid5(FUNDING_ID_NAMESPACE, f"unabated-taker-funding-subaccount-{subaccount}"))


def _signed(client: KalshiOrderClient, method: str, path: str, body: dict | None = None) -> tuple[int, object]:
    try:
        return client.request(method, path, body)
    except OSError as error:
        raise SetupError(f"{method} {path}: no answer from Kalshi ({error})") from error


# ---- step 0: nothing can trade, and the key has somewhere to go ---------------------

def ensure_stopped(stop_path: Path, out: TextIO) -> None:
    """The stop file in place: the restart that arms the bot sends no order."""
    if stop_path.exists():
        out.write(f"Stop file {stop_path} is in place; the bot sends no order until you rm it.\n")
        return
    try:
        stop_path.touch()
    except OSError as error:
        raise SetupError(
            f"could not create the stop file {stop_path} ({type(error).__name__}): nothing done"
        ) from error
    out.write(f"Put the stop file {stop_path} in place; the bot sends no order until you rm it.\n")


def ensure_key_dir_writable(key_path: Path) -> None:
    """A key Kalshi mints can be saved only once: find out first that it can be."""
    probe = key_path.with_name(key_path.name + ".probe")
    try:
        descriptor = os.open(probe, os.O_WRONLY | os.O_CREAT | os.O_EXCL, 0o600)
        os.close(descriptor)
        probe.unlink()
    except OSError as error:
        raise SetupError(f"cannot write the bot key file next to {key_path} ({type(error).__name__}): "
                         f"nothing asked of Kalshi") from error


# ---- step 1: the subaccount --------------------------------------------------------

def subaccount_balances(client: KalshiOrderClient) -> dict[int, float]:
    """Dollars per subaccount number (summed over exchange shards), the primary included."""
    status, reply = _signed(client, "GET", SUBACCOUNT_BALANCES_PATH)
    rows = reply.get("subaccount_balances") if isinstance(reply, dict) else None
    if status != 200 or not isinstance(rows, list):
        raise SetupError(f"GET {SUBACCOUNT_BALANCES_PATH}: expected 200 with subaccount_balances, "
                         f"got HTTP {status} {str(reply)[:300]}")
    balances: dict[int, float] = {}
    for row in rows:
        number = row.get("subaccount_number") if isinstance(row, dict) else None
        if isinstance(number, bool) or not isinstance(number, int):
            raise SetupError(f"unreadable subaccount balance row {row!r}")
        try:
            dollars = float(row.get("balance"))
        except (TypeError, ValueError) as error:
            raise SetupError(f"unreadable balance on subaccount {number}: {row.get('balance')!r}") from error
        balances[number] = balances.get(number, 0.0) + dollars
    return balances


def create_subaccount(client: KalshiOrderClient, out: TextIO) -> int:
    """A new numbered subaccount; asks for the Advanced API tier once if Kalshi refuses."""
    status, reply = _signed(client, "POST", SUBACCOUNTS_PATH, {})
    if status != 201:
        out.write(f"Kalshi refused a new subaccount (HTTP {status}); asking for the Advanced API tier it needs.\n")
        upgrade_status, upgrade_reply = _signed(client, "POST", UPGRADE_PATH, {})
        if upgrade_status != 201:
            raise SetupError(f"the Advanced API tier was refused: HTTP {upgrade_status} {str(upgrade_reply)[:300]} "
                             f"(the create answered HTTP {status} {str(reply)[:300]})")
        status, reply = _signed(client, "POST", SUBACCOUNTS_PATH, {})
    number = reply.get("subaccount_number") if isinstance(reply, dict) else None
    if status != 201 or isinstance(number, bool) or not isinstance(number, int) or number < 1:
        raise SetupError(f"POST {SUBACCOUNTS_PATH}: expected 201 with a subaccount_number, "
                         f"got HTTP {status} {str(reply)[:300]}")
    out.write(f"Made subaccount {number}.\n")
    return number


def choose_subaccount(client: KalshiOrderClient, wanted: int | None, out: TextIO) -> int:
    """`wanted` when Kalshi lists it, else the one numbered subaccount, else a new one."""
    numbered = sorted(number for number in subaccount_balances(client) if number != PRIMARY_SUBACCOUNT)
    if wanted is not None:
        if wanted not in numbered:
            raise SetupError(f"Kalshi lists no subaccount {wanted} on this account (it lists {numbered or 'none'})")
        return wanted
    if len(numbered) > 1:
        raise SetupError(f"the account already has subaccounts {numbered}: name the bot's with --subaccount N")
    if numbered:
        out.write(f"Using the account's subaccount {numbered[0]}.\n")
        return numbered[0]
    return create_subaccount(client, out)


# ---- step 2: the money -------------------------------------------------------------

def fund_once(client: KalshiOrderClient, subaccount: int, out: TextIO) -> float:
    """Moves BOT_BANKROLL from the primary account while the subaccount holds
    $0, and checks it landed; returns what the subaccount holds after."""
    held = subaccount_balances(client).get(subaccount, 0.0)
    if held > 0:
        out.write(f"Subaccount {subaccount} already holds ${held:,.2f}; no transfer.\n")
        return held
    body = {"client_transfer_id": funding_transfer_id(subaccount), "from_subaccount": PRIMARY_SUBACCOUNT,
            "to_subaccount": subaccount, "amount_cents": int(round(BOT_BANKROLL * 100))}
    status, reply = _signed(client, "POST", TRANSFER_PATH, body)
    if status != 200:
        raise SetupError(f"the ${BOT_BANKROLL:,.0f} transfer to subaccount {subaccount} was refused: HTTP {status} "
                         f"{str(reply)[:300]} (the setup funds a subaccount once; a top-up is your own transfer)")
    held = subaccount_balances(client).get(subaccount, 0.0)
    if held < BOT_BANKROLL - FUNDING_TOLERANCE_DOLLARS:
        raise SetupError(f"Kalshi answered 200 to the ${BOT_BANKROLL:,.0f} transfer but subaccount {subaccount} "
                         f"holds ${held:,.2f}: check your Kalshi balances before running this again")
    out.write(f"Moved ${BOT_BANKROLL:,.2f} from your primary account; subaccount {subaccount} holds ${held:,.2f}.\n")
    return held


# ---- step 3: the key ---------------------------------------------------------------

def write_private_key(path: Path, pem: str) -> None:
    """The key file at mode 0600, swapped in whole (never a half-written key)."""
    from cryptography.hazmat.primitives import serialization
    serialization.load_pem_private_key(pem.encode(), password=None)  # raises on anything but a usable key
    staging = path.with_name(path.name + ".new")
    try:
        descriptor = os.open(staging, os.O_WRONLY | os.O_CREAT | os.O_TRUNC, 0o600)
        with os.fdopen(descriptor, "w") as key_file:
            key_file.write(pem)
            key_file.flush()
            os.fsync(key_file.fileno())
        os.chmod(staging, 0o600)
        os.replace(staging, path)
    finally:
        # A failed write leaves no stray copy of the key behind.
        staging.unlink(missing_ok=True)


def mint_locked_key(client: KalshiOrderClient, subaccount: int, key_path: Path, now: datetime) -> str:
    """A new key locked to the subaccount, its private key written to
    `key_path`; returns its id. The reply holds the private key, so no error
    here ever echoes a 201 reply."""
    body = {"name": f"unabated-taker-subaccount-{subaccount}-{now:%Y%m%d}", "key_type": "rsa",
            "scopes": list(BOT_KEY_SCOPES), "subaccount": subaccount}
    status, reply = _signed(client, "POST", GENERATE_KEY_PATH, body)
    if status != 201:
        raise SetupError(f"POST {GENERATE_KEY_PATH} was refused: HTTP {status} {str(reply)[:300]}")
    key_id = reply.get("api_key_id") if isinstance(reply, dict) else None
    pem = reply.get("private_key") if isinstance(reply, dict) else None
    if not isinstance(key_id, str) or not key_id or not isinstance(pem, str):
        fields = sorted(reply) if isinstance(reply, dict) else type(reply).__name__
        raise SetupError(f"POST {GENERATE_KEY_PATH}: expected api_key_id and private_key, got fields {fields}")
    try:
        write_private_key(key_path, pem)
    except (OSError, ValueError) as error:
        raise SetupError(f"Kalshi made key {key_id} but its private key could not be saved to {key_path} "
                         f"({type(error).__name__}); revoke {key_id} at kalshi.com and run this again") from error
    return key_id


# ---- step 4: the lock, two ways ----------------------------------------------------

def listed_key(client: KalshiOrderClient, key_id: str) -> dict | None:
    """Kalshi's own record of the key (GET /api_keys on the main key), or None
    when the listing does not carry it."""
    status, reply = _signed(client, "GET", API_KEYS_PATH)
    keys = reply.get("api_keys") if isinstance(reply, dict) else None
    if status != 200 or not isinstance(keys, list):
        raise SetupError(f"GET {API_KEYS_PATH}: expected 200 with api_keys, got HTTP {status} {str(reply)[:300]}")
    for key in keys:
        if isinstance(key, dict) and key.get("api_key_id") == key_id:
            return key
    return None


def key_restriction_problem(record: dict | None, key_id: str, subaccount: int) -> str | None:
    """Why Kalshi's record of the key is not the lock the bot needs, or None:
    restricted to `subaccount`, scopes exactly BOT_KEY_SCOPES, no FCM binding."""
    if record is None:
        return f"Kalshi's key listing does not carry {key_id}"
    if record.get("subaccount") != subaccount:
        return f"Kalshi lists {key_id} with subaccount {record.get('subaccount')!r}, not {subaccount}"
    scopes = record.get("scopes")
    if not isinstance(scopes, list) or set(scopes) != set(BOT_KEY_SCOPES) or len(scopes) != len(BOT_KEY_SCOPES):
        return f"Kalshi lists {key_id} with scopes {scopes!r}, not {list(BOT_KEY_SCOPES)}"
    if record.get("fcm_subtrader_id"):
        return f"Kalshi lists {key_id} bound to an FCM subtrader"
    return None


def confirmed_lock(make_bot_client: Callable[[str, Path], KalshiOrderClient], key_id: str, key_path: Path,
                   subaccount: int, attempts: int, sleep: Callable[[float], None]) -> str | None:
    """None once the key passes the desk's lock check, else the last reason."""
    try:
        client = make_bot_client(key_id, key_path)
    except (OSError, ValueError) as error:
        return f"the key file {key_path} could not be read ({type(error).__name__})"
    problem = "not checked"
    for attempt in range(attempts):
        if attempt:
            sleep(NEW_KEY_CHECK_RETRY_SEC)
        problem, _balance = key_lock_problem(client, subaccount, sleep=sleep)
        if problem is None:
            return None
    return problem


def _configured_key_stands(main_client: KalshiOrderClient, make_bot_client: Callable[[str, Path], KalshiOrderClient],
                           key_id: str | None, key_path: Path, configured_subaccount: int | None, subaccount: int,
                           sleep: Callable[[float], None]) -> bool:
    """The key the service is configured with is already the locked key for this subaccount."""
    if configured_subaccount != subaccount or not key_id or not key_path.exists():
        return False
    if key_restriction_problem(listed_key(main_client, key_id), key_id, subaccount) is not None:
        return False
    return confirmed_lock(make_bot_client, key_id, key_path, subaccount, 1, sleep) is None


def run_setup(main_client: KalshiOrderClient, make_bot_client: Callable[[str, Path], KalshiOrderClient],
              key_path: Path, stop_path: Path, wanted_subaccount: int | None, configured_key_id: str | None,
              configured_subaccount: int | None, out: TextIO, now: datetime,
              sleep: Callable[[float], None] = time.sleep) -> int:
    """The steps of the module docstring, in order. Returns the process exit code."""
    try:
        ensure_stopped(stop_path, out)
        ensure_key_dir_writable(key_path)
        subaccount = choose_subaccount(main_client, wanted_subaccount, out)
        held = fund_once(main_client, subaccount, out)
        if _configured_key_stands(main_client, make_bot_client, configured_key_id, key_path, configured_subaccount,
                                  subaccount, sleep):
            out.write(f"The configured bot key {configured_key_id} is already locked to subaccount {subaccount}; "
                      f"no new key.\n")
            return 0
        key_id = mint_locked_key(main_client, subaccount, key_path, now)
        out.write(f"Made the bot's key {key_id}, locked to subaccount {subaccount}; saved to {key_path} (mode 0600).\n")
        problem = (key_restriction_problem(listed_key(main_client, key_id), key_id, subaccount)
                   or confirmed_lock(make_bot_client, key_id, key_path, subaccount, NEW_KEY_CHECK_ATTEMPTS, sleep))
        if problem is not None:
            # Never leave a key that may reach the primary account where the bot's key goes.
            key_path.unlink(missing_ok=True)
            raise SetupError(f"the new key {key_id} did not pass the lock check ({problem}); its file is deleted. "
                             f"Revoke {key_id} at kalshi.com")
    except SetupError as error:
        out.write(f"STOPPED: {error}\n")
        return 1
    out.write(
        f"\nDone. Subaccount {subaccount} holds ${held:,.2f}; key {key_id} trades only there and cannot move money.\n"
        f"Set these in unabated_ticket/bets_service/.env on the host, then restart the bets service:\n"
        f"  KALSHI_TAKER_API_KEY_ID={key_id}\n"
        f"  KALSHI_TAKER_SUBACCOUNT={subaccount}\n"
        f"Then revoke the bot's previous key at kalshi.com. The stop file {stop_path} stays until you rm it.\n")
    return 0


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description="Set up the automated taker's own Kalshi subaccount and locked key.")
    parser.add_argument("--subaccount", type=int, choices=range(1, 64), metavar="N",
                        help="use this existing subaccount (1-63) instead of finding or making one")
    args = parser.parse_args(argv)
    if not config.KALSHI_API_KEY_ID or not config.KALSHI_PRIVATE_KEY_PATH:
        print("STOPPED: KALSHI_API_KEY_ID / KALSHI_PRIVATE_KEY_PATH (your main key) are not set", file=sys.stderr)
        return 1
    if not config.KALSHI_TAKER_PRIVATE_KEY_PATH:
        print("STOPPED: KALSHI_TAKER_PRIVATE_KEY_PATH (where the bot's key goes) is not set", file=sys.stderr)
        return 1
    try:
        main_client = KalshiOrderClient(config.KALSHI_API_KEY_ID, config.KALSHI_PRIVATE_KEY_PATH,
                                        config.KALSHI_BASE_URL)
    except (OSError, ValueError) as error:
        print(f"STOPPED: your main key {config.KALSHI_PRIVATE_KEY_PATH} could not be read ({type(error).__name__})",
              file=sys.stderr)
        return 1

    def make_bot_client(key_id: str, key_path: Path) -> KalshiOrderClient:
        return KalshiOrderClient(key_id, str(key_path), config.KALSHI_BASE_URL)

    return run_setup(main_client, make_bot_client, Path(config.KALSHI_TAKER_PRIVATE_KEY_PATH),
                     config.TAKER_STOP_PATH, args.subaccount, config.KALSHI_TAKER_API_KEY_ID,
                     config.KALSHI_TAKER_SUBACCOUNT, sys.stdout, datetime.now(timezone.utc))


if __name__ == "__main__":
    sys.exit(main())
