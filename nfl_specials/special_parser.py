"""Parse a Wagerzon NFL trifecta / superfecta description into legs.

    parse_fecta("SEAHAWKS TRIFECTA (1Q -½, 1H -4½ & GM -7½)")
      -> Fecta(team='SEA', prop_type='TRIFECTA',
               legs=(Leg('spread','Q1',-0.5), Leg('spread','H1',-4.5), Leg('spread','GM',-7.5)))

Leg grammar inside the parentheses (items split on ',' and '&'):
    SCR 1ST         team makes the game's first score
    1Q / 1H / GM    team WINS that period outright. Quarter and half wins are
                    3-way: a tied quarter or half LOSES (Wagerzon grading).
    1Q -½ / GM +3½  team covers that spread in that period

Anything else returns None with a reason — a guessed leg would print a fake
edge, a skipped special costs nothing.
"""
from __future__ import annotations

import re
from dataclasses import dataclass

from nfl_specials.teams import resolve_wz_team

FECTA_RE = re.compile(r"^(\S+) (TRIFECTA|SUPERFECTA) \((.+)\)$")
PERIOD_TOKENS = {"1Q": "Q1", "1H": "H1", "GM": "GM"}
LEG_RE = re.compile(r"(1Q|1H|GM)(?: ([+-]\d+(?:\.\d+)?))?")


@dataclass(frozen=True)
class Leg:
    kind: str                 # 'scores_first' | 'win' | 'spread'
    period: str = "GM"        # 'Q1' | 'H1' | 'GM'
    line: float | None = None  # team-perspective spread, spread legs only

    def describe(self, team: str) -> str:
        if self.kind == "scores_first":
            return f"{team} scores first"
        if self.kind == "win":
            return f"{team} wins {self.period}"
        return f"{team} {self.line:+g} {self.period}"


@dataclass(frozen=True)
class Fecta:
    team: str                 # nflverse abbreviation, e.g. 'SEA'
    prop_type: str            # 'TRIFECTA' | 'SUPERFECTA'
    legs: tuple[Leg, ...]


@dataclass(frozen=True)
class ParseFailure:
    reason: str


def normalize(text: str) -> str:
    text = text.upper().replace("½", ".5")
    text = re.sub(r"(?<![\d.])\.5", "0.5", text)   # "-.5" -> "-0.5"
    return re.sub(r"\s+", " ", text).strip()


def parse_fecta(description: str) -> Fecta | ParseFailure:
    text = normalize(description)
    match = FECTA_RE.match(text)
    if not match:
        return ParseFailure("not a trifecta/superfecta")
    team = resolve_wz_team(match.group(1))
    if team is None:
        return ParseFailure(f"unknown team '{match.group(1)}'")

    legs = []
    for item in re.split(r"\s*(?:,|&)\s*", match.group(3)):
        if item == "SCR 1ST":
            legs.append(Leg(kind="scores_first"))
            continue
        leg_match = LEG_RE.fullmatch(item)
        if not leg_match:
            return ParseFailure(f"unrecognized leg '{item}'")
        period = PERIOD_TOKENS[leg_match.group(1)]
        if leg_match.group(2) is None:
            legs.append(Leg(kind="win", period=period))
        else:
            legs.append(Leg(kind="spread", period=period, line=float(leg_match.group(2))))

    expected_legs = 3 if match.group(2) == "TRIFECTA" else 4
    if len(legs) != expected_legs:
        return ParseFailure(f"{match.group(2)} with {len(legs)} legs, expected {expected_legs}")
    return Fecta(team=team, prop_type=match.group(2), legs=tuple(legs))
