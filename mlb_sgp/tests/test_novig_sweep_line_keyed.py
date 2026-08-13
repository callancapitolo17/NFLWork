"""Issue #64: the Novig SWEEP must price each target at ITS OWN lines.

The defect these tests pin: ``novig._price_sgps`` used to collapse every
TargetLine of a game into one ``fg_spread_line`` / ``fg_total_line``,
resolve the leg UUIDs ONCE per game, then price all of that game's
targets with those same legs while stamping each row with its own
``spread_line`` / ``total_line`` — N rows carrying the SAME price under
DIFFERENT labels, which poisons the maker's cross-book median.

These drive the REAL ``price_sgps`` orchestration with the network faked
at the two legacy seams it imports (``fetch_event_legs`` for the raw
market tree, ``submit_parlay`` for pricing), so the control flow under
test is production's.

The fake pricer returns a decimal that is a deterministic function of the
outcome-id pair it was handed. That is the whole trick: if two targets are
priced with the same legs, they get the same decimal, and the "different
totals must give different prices" assertions fail — which is exactly what
happened before the fix.
"""
from __future__ import annotations

import inspect
import sys
import threading
from datetime import datetime, timezone
from pathlib import Path
from types import SimpleNamespace
from unittest.mock import MagicMock

MLB_SGP_DIR = Path(__file__).resolve().parents[1]
if str(MLB_SGP_DIR) not in sys.path:
    sys.path.insert(0, str(MLB_SGP_DIR))

from mlb_sgp import novig  # noqa: E402
from mlb_sgp._shared import TargetLine  # noqa: E402

CT = datetime(2026, 8, 3, 23, 5, tzinfo=timezone.utc)
GAME_ID = "g1"
HOME_SYM, AWAY_SYM = "NYY", "BOS"
SPREAD = -1.5


# --------------------------------------------------------------------- #
# Fake Novig market tree — ids encode the line so an assertion can read  #
# straight off a submitted parlay which line was actually priced.        #
# --------------------------------------------------------------------- #

def _spread_market(strike: float, mtype: str = "SPREAD") -> dict:
    """One SPREAD market. ``strike`` is HOME-perspective, matching Novig."""
    return {
        "type": mtype, "strike": strike, "is_consensus": True,
        "outcomes": [
            {"id": f"home_spread@{strike}", "available": 0.5,
             "competitor": {"symbol": HOME_SYM}},
            {"id": f"away_spread@{strike}", "available": 0.5,
             "competitor": {"symbol": AWAY_SYM}},
        ],
    }


def _total_market(strike: float, mtype: str = "TOTAL") -> dict:
    return {
        "type": mtype, "strike": strike, "is_consensus": True,
        "outcomes": [
            {"id": f"over@{strike}", "available": 0.5,
             "description": f"Over {strike}"},
            {"id": f"under@{strike}", "available": 0.5,
             "description": f"Under {strike}"},
        ],
    }


def _money_market(mtype: str = "MONEY") -> dict:
    return {
        "type": mtype, "strike": None, "is_consensus": True,
        "outcomes": [
            {"id": "home_ml", "available": 0.55,
             "competitor": {"symbol": HOME_SYM}},
            {"id": "away_ml", "available": 0.45,
             "competitor": {"symbol": AWAY_SYM}},
        ],
    }


def _markets(totals=(7.5, 8.5, 9.5), spreads=(SPREAD,), money=True) -> list:
    out = [_spread_market(s) for s in spreads]
    out += [_total_market(t) for t in totals]
    if money:
        out.append(_money_market())
    return out


def _line_of(outcome_id: str) -> float | None:
    """Recover the line an id was minted at ('over@8.5' -> 8.5)."""
    if "@" not in outcome_id:
        return None
    return float(outcome_id.split("@", 1)[1])


class _Recorder:
    """Stands in for ``scraper_novig_sgp.submit_parlay``: remembers every
    outcome-id pair it was asked to price, and answers each DISTINCT pair
    with its own decimal.

    One decimal per leg set is the whole mechanism — a row's price is a
    fingerprint of the legs that produced it, so ``ids_for`` can read back
    which lines were really priced, and two targets sharing one leg set
    necessarily share a price (which is the defect).

    Decimals land in [2.0, 2.01): both legs carry available=0.5, so the
    naive-multiply ceiling is 4.0 * SANITY_MULT_RATIO = 6.0 and nothing is
    dropped by the sanity filter for the wrong reason.
    """

    def __init__(self, price_fn=None):
        self.calls: list[tuple[str, ...]] = []
        self._decimal_by_legs: dict[frozenset, float] = {}
        self._legs_by_decimal: dict[float, tuple[str, ...]] = {}
        # price_sgps prices combos on a thread pool.
        self._lock = threading.Lock()
        self._price_fn = price_fn

    def __call__(self, session, outcome_ids, verbose=False, **kwargs):
        ids = tuple(outcome_ids)
        key = frozenset(ids)
        with self._lock:
            self.calls.append(ids)
            dec = self._decimal_by_legs.get(key)
            if dec is None:
                dec = (self._price_fn(ids) if self._price_fn
                       else 2.0 + 0.0001 * (len(self._decimal_by_legs) + 1))
                self._decimal_by_legs[key] = dec
                self._legs_by_decimal.setdefault(round(dec, 4), ids)
        return {"decimal": dec, "american": 150}, False

    def ids_for(self, decimal: float) -> tuple[str, ...]:
        """The leg pair that produced this row's decimal."""
        ids = self._legs_by_decimal.get(round(decimal, 4))
        if ids is None:
            raise AssertionError(f"no recorded price call produced {decimal}")
        return ids


# A coherent 4-cell joint distribution at each adjacent alt, quoted at an
# 18% hold (Novig's observed hold). Over mass shifts down as the total rises,
# which is what gives the integer-line derivation its push mass to solve for.
_ALT_JOINT = {
    7.5: {"home_over": 0.25, "home_under": 0.20,
          "away_over": 0.28, "away_under": 0.27},
    8.5: {"home_over": 0.21, "home_under": 0.24,
          "away_over": 0.24, "away_under": 0.31},
}
_ALT_HOLD = 1.18


def _coherent_alt_price(outcome_ids) -> float:
    """Price a (spread, total) pair off ``_ALT_JOINT`` instead of the default
    fingerprint decimals — used where real derivation math must succeed."""
    total_leg = next(i for i in outcome_ids
                     if i.startswith("over@") or i.startswith("under@"))
    spread_leg = next(i for i in outcome_ids if "spread@" in i)
    cell = f"{_side_of(spread_leg)}_{'over' if 'over@' in total_leg else 'under'}"
    return 1.0 / (_ALT_JOINT[_line_of(total_leg)][cell] * _ALT_HOLD)


def _side_of(outcome_id: str) -> str:
    return "home" if outcome_id.startswith("home") else "away"


def _assert_legs_match_label(row, ids):
    """A row's priced legs must carry its own labelled lines and sides.

    Novig models a spread as ONE market with a single HOME-perspective
    strike and two outcomes, so both sides' ids encode the home-perspective
    line — the side itself is carried by the outcome, which is what the
    ``home``/``away`` prefix stands in for here.
    """
    total_leg = next(i for i in ids
                     if i.startswith("over@") or i.startswith("under@"))
    assert _line_of(total_leg) == row.total_line, (
        f"row labelled total={row.total_line} was priced on {total_leg}")
    assert ("over@" in total_leg) == ("+ Over" in row.combo), (
        f"row {row.combo!r} was priced on {total_leg}")

    if row.spread_line is None:            # ML x total row
        ml_leg = next(i for i in ids if i.endswith("_ml"))
        assert _side_of(ml_leg) == ("home" if "Home" in row.combo else "away")
        return
    spread_leg = next(i for i in ids if "spread@" in i)
    assert _line_of(spread_leg) == row.spread_line, (
        f"row labelled spread={row.spread_line} was priced on {spread_leg}")
    assert _side_of(spread_leg) == ("home" if "Home" in row.combo else "away")


def _target(total: float, spread: float = SPREAD,
            period: str = "FG") -> TargetLine:
    return TargetLine(game_id=GAME_ID, home_team="New York Yankees",
                      away_team="Boston Red Sox", commence_time=CT,
                      period=period, spread=spread, total=total)


def _run(monkeypatch, targets, markets, periods=("FG",), price_fn=None):
    """Drive the real price_sgps with the two legacy seams faked.

    Returns (rows, recorder). ``match_events`` is faked because it only
    does team-name/commence-time resolution, which is not what is under
    test; everything downstream of it is production code.
    """
    import scraper_novig_sgp as legacy

    matched = {
        "game_id": GAME_ID, "nv_event_id": "e1",
        "nv_home": "New York Yankees", "nv_away": "Boston Red Sox",
        "nv_home_sym": HOME_SYM, "nv_away_sym": AWAY_SYM,
        # Whatever the collapsed grouping produced upstream; the ML phase
        # reads these, the spread/total phase must NOT.
        "fg_spread_line": targets[-1].spread, "fg_total_line": targets[-1].total,
        "f5_spread_line": None, "f5_total_line": None,
    }
    monkeypatch.setattr(legacy, "match_events", lambda evs, lines: [matched])

    # fetch_event_legs' FIRST element is the single-line legs dict the sweep
    # used to price from. Hand back the legs for the LAST target's line —
    # the collapse the defect was built on — so a fixed sweep must be
    # ignoring it. Its SECOND element (the raw market tree) is the input
    # the fix reads.
    last = targets[-1]
    collapsed = {
        "fg": {"home_spread": {"id": f"home_spread@{last.spread}",
                               "available": 0.5},
               "away_spread": {"id": f"away_spread@{-last.spread}",
                               "available": 0.5},
               "over": {"id": f"over@{last.total}", "available": 0.5},
               "under": {"id": f"under@{last.total}", "available": 0.5},
               "home_ml": {"id": "home_ml", "available": 0.55},
               "away_ml": {"id": "away_ml", "available": 0.45}},
        "f5": {"home_spread": None, "away_spread": None, "over": None,
               "under": None, "home_ml": None, "away_ml": None},
    }
    monkeypatch.setattr(legacy, "fetch_event_legs",
                        lambda session, game, verbose=False, **kw:
                        (collapsed, markets))

    recorder = _Recorder(price_fn=price_fn)
    monkeypatch.setattr(legacy, "submit_parlay", recorder)

    client = MagicMock()
    client.list_events.return_value = [
        SimpleNamespace(event_id="e1", home_team="New York Yankees",
                        away_team="Boston Red Sox", home_sym=HOME_SYM,
                        away_sym=AWAY_SYM, start_time=CT)]
    rows = novig.price_sgps(targets, periods=periods, client=client)
    return rows, recorder


def _spread_total_rows(rows):
    """Rows from the spread x total grid (drops the ML x total phase)."""
    return [r for r in rows if r.spread_line is not None]


# --------------------------------------------------------------------- #
# The defect                                                            #
# --------------------------------------------------------------------- #

def test_three_totals_request_three_distinct_leg_sets(monkeypatch):
    """3 targets differing only in total must resolve 3 different leg sets.

    Pre-fix this produced 4 distinct leg sets (one per combo flavor,
    shared across every total) instead of 12.
    """
    totals = [7.5, 8.5, 9.5]
    rows, rec = _run(monkeypatch, [_target(t) for t in totals], _markets())

    spread_total_calls = [ids for ids in rec.calls
                          if any(i.startswith("home_spread@")
                                 or i.startswith("away_spread@") for i in ids)]
    assert len({frozenset(ids) for ids in spread_total_calls}) == 12, (
        "each (total x combo) must be its own leg set; got "
        f"{len({frozenset(i) for i in spread_total_calls})}")

    over_lines = {_line_of(i) for ids in spread_total_calls for i in ids
                  if i.startswith("over@") or i.startswith("under@")}
    assert over_lines == set(totals)


def test_rows_at_different_totals_carry_different_prices(monkeypatch):
    """The acceptance criterion: never the same price under different labels."""
    totals = [7.5, 8.5, 9.5]
    rows, _ = _run(monkeypatch, [_target(t) for t in totals], _markets())

    for combo in ("Home Spread + Over", "Home Spread + Under",
                  "Away Spread + Over", "Away Spread + Under"):
        decs = {r.total_line: r.sgp_decimal
                for r in _spread_total_rows(rows) if r.combo == combo}
        assert set(decs) == set(totals), f"{combo}: missing a total"
        assert len(set(decs.values())) == 3, (
            f"{combo}: same price under different total labels -> {decs}")


def test_every_row_was_priced_with_the_legs_its_label_claims(monkeypatch):
    """Strongest form: label and priced legs must agree, row by row."""
    targets = [_target(t) for t in (7.5, 8.5, 9.5)]
    rows, rec = _run(monkeypatch, targets, _markets())

    assert _spread_total_rows(rows)
    for r in _spread_total_rows(rows):
        _assert_legs_match_label(r, rec.ids_for(r.sgp_decimal))


def test_alt_spreads_resolve_per_target_too(monkeypatch):
    """The same collapse applied to the spread axis — pin both axes."""
    spreads = [-1.5, -2.5]
    targets = [_target(8.5, spread=s) for s in spreads]
    rows, rec = _run(monkeypatch, targets,
                     _markets(totals=(8.5,), spreads=tuple(spreads)))

    assert {r.spread_line for r in _spread_total_rows(rows)} == set(spreads)
    for r in _spread_total_rows(rows):
        _assert_legs_match_label(r, rec.ids_for(r.sgp_decimal))


def test_f5_targets_resolve_against_f5_markets(monkeypatch):
    """F5 targets must read the *_1H market tree at their own lines."""
    targets = [_target(4.5, spread=-0.5, period="F5"),
               _target(5.5, spread=-0.5, period="F5")]
    markets = ([_spread_market(-0.5, "SPREAD_1H")]
               + [_total_market(t, "TOTAL_1H") for t in (4.5, 5.5)]
               + [_money_market("MONEY_1H")]
               # FG markets at unrelated lines must never be picked up.
               + _markets(totals=(9.5,), spreads=(-3.5,)))
    rows, rec = _run(monkeypatch, targets, markets, periods=("F5",))

    f5_rows = _spread_total_rows(rows)
    assert {r.total_line for r in f5_rows} == {4.5, 5.5}
    assert all(r.combo.startswith("F5 ") for r in f5_rows)
    for r in f5_rows:
        ids = rec.ids_for(r.sgp_decimal)
        total_leg = next(i for i in ids
                         if i.startswith("over@") or i.startswith("under@"))
        assert _line_of(total_leg) == r.total_line


# --------------------------------------------------------------------- #
# Guards: the fix must not invent rows, drop the fallback, or add calls  #
# --------------------------------------------------------------------- #

def test_unoffered_total_produces_no_row(monkeypatch):
    """A half-point total the book doesn't quote yields nothing — never an
    interpolated row invented to fill the gap."""
    targets = [_target(8.5), _target(12.5)]
    rows, _ = _run(monkeypatch, targets, _markets(totals=(7.5, 8.5, 9.5)))

    assert {r.total_line for r in _spread_total_rows(rows)} == {8.5}
    assert not [r for r in rows if r.source == novig.SOURCE_LABEL_FALLBACK]


def test_integer_total_still_routes_to_the_interpolated_fallback(monkeypatch):
    """try_integer_fallback_nv was already per-target correct — keep it.

    Target total 8.0 with 7.5 / 8.5 both offered must still derive an
    integer-line price and be labelled novig_interpolated.

    The derivation applies real bounds checks (per-alt vig in [1.05, 1.30],
    push mass in [0.03, 0.18]), so this test prices off a coherent joint
    distribution rather than the default fingerprint decimals.
    """
    targets = [_target(8.0)]
    rows, _ = _run(monkeypatch, targets, _markets(totals=(7.5, 8.5)),
                   price_fn=_coherent_alt_price)

    grid = _spread_total_rows(rows)
    assert grid, "integer fallback produced no rows"
    assert {r.source for r in grid} == {novig.SOURCE_LABEL_FALLBACK}
    assert {r.total_line for r in grid} == {8.0}


def test_price_call_count_does_not_exceed_four_per_target(monkeypatch):
    """Acceptance criterion: no increase in Novig's wire-call count.

    Four combos per resolvable target plus one four-call ML x total phase
    per game is the pre-fix budget; the fix must stay inside it.
    """
    totals = [7.5, 8.5, 9.5]
    targets = [_target(t) for t in totals]
    _rows, rec = _run(monkeypatch, targets, _markets(totals=tuple(totals)))

    assert len(rec.calls) <= 4 * len(targets) + 4


def test_single_target_per_game_matches_the_line_it_asked_for(monkeypatch):
    """The dashboard shape (one target per game) must be unchanged."""
    rows, rec = _run(monkeypatch, [_target(8.5)], _markets())

    grid = _spread_total_rows(rows)
    assert len(grid) == 4
    assert {r.total_line for r in grid} == {8.5}
    assert {r.spread_line for r in grid} == {SPREAD}
    for ids in rec.calls:
        for i in ids:
            line = _line_of(i)
            if i.startswith("over@") or i.startswith("under@"):
                assert line == 8.5


def test_moneyline_rows_are_labelled_with_the_total_they_priced(monkeypatch):
    """ML x total is priced once per game; its label must match its legs."""
    targets = [_target(7.5), _target(8.5)]
    rows, rec = _run(monkeypatch, targets, _markets())

    ml_rows = [r for r in rows if r.spread_line is None]
    assert ml_rows, "ML x total phase produced no rows"
    for r in ml_rows:
        _assert_legs_match_label(r, rec.ids_for(r.sgp_decimal))


# --------------------------------------------------------------------- #
# Signature contract                                                     #
# --------------------------------------------------------------------- #

def test_build_line_structure_period_is_keyword_only():
    """A new param must be keyword-only so no positional caller shifts."""
    sig = inspect.signature(novig.build_line_structure)
    sig.bind([], HOME_SYM, AWAY_SYM)                    # on-demand call site
    sig.bind([], HOME_SYM, AWAY_SYM, period="f5")
    try:
        sig.bind([], HOME_SYM, AWAY_SYM, "f5")
    except TypeError:
        return
    raise AssertionError("period must be keyword-only")


def test_build_line_structure_reads_the_f5_market_types():
    """F5 bucket comes from *_1H markets, and never from the FG ones."""
    markets = [_spread_market(-0.5, "SPREAD_1H"), _total_market(4.5, "TOTAL_1H"),
               _money_market("MONEY_1H"),
               _spread_market(-1.5), _total_market(8.5), _money_market()]

    f5 = novig.build_line_structure(markets, HOME_SYM, AWAY_SYM, period="f5")
    assert set(f5["over"]) == {4.5}
    assert set(f5["home_spread"]) == {-0.5}

    fg = novig.build_line_structure(markets, HOME_SYM, AWAY_SYM)
    assert set(fg["over"]) == {8.5}
    assert set(fg["home_spread"]) == {-1.5}
