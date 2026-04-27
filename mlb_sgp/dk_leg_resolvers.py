"""
DK leg-type → selection-id resolvers for Plan #2's trifecta scraper.

Each resolver takes:
  - leg: a leg-spec dict matching parse_legs.R::TOKEN_REGISTRY output, e.g.:
      {'type': 'scores_first'}
      {'type': 'wins_period', 'period': 'F5'}
      {'type': 'team_total_under', 'line': 2.5}
  - side: 'home' or 'away' — determines which selection within the market we pick
  - event_state: the cached DK SGP event payload (full JSON returned by
    parlays/v1/sgp/events/{event_id})
  - team_names: dict {'home': 'Cubs', 'away': 'Padres'} — DK team-name strings
    extracted from event_state at scrape time so resolvers can match selections
    by name (DK's selection.name field is e.g. 'CHI Cubs' or 'SD Padres').

Each resolver returns:
  - str selection_id when the market+selection is found and SGP-eligible
  - None when the market is missing OR not SGP-eligible (R-side blend
    degrades to model-only for that row)
"""
from __future__ import annotations
from typing import Optional


# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------

def find_market_by_name(payload: dict, target_name: str) -> Optional[dict]:
    """Walk payload, return first market dict whose name matches target_name
    AND is SGP-eligible (has 'SGP' in tags) AND has selections.
    """
    def walk(node):
        if isinstance(node, dict):
            n = node.get('marketName') or node.get('name')
            if (
                n == target_name
                and node.get('selections')
                and 'SGP' in (node.get('tags') or [])
            ):
                return node
            for v in node.values():
                r = walk(v)
                if r is not None:
                    return r
        elif isinstance(node, list):
            for v in node:
                r = walk(v)
                if r is not None:
                    return r
        return None
    return walk(payload)


def find_team_total_market(
    payload: dict, team_label: str, side_target: str, line: float
) -> Optional[dict]:
    """Find a "<TEAM>: Team Total Runs" market with a selection matching the
    given line + over/under. Returns the selection dict (with 'id'), or None.

    Multiple team-total markets exist per game (full-game, alt lines, by-inning).
    We want the FULL-game one for the given team and line.
    """
    target_market_name = f"{team_label}: Team Total Runs"
    market = find_market_by_name(payload, target_market_name)
    if market is None:
        return None
    # Selections should have labels like "Over" / "Under" plus a 'point' field
    for sel in market.get('selections', []):
        sel_label = (sel.get('label') or sel.get('name') or '').lower()
        sel_point = sel.get('point')
        if (
            (side_target == 'over'  and 'over'  in sel_label) or
            (side_target == 'under' and 'under' in sel_label)
        ):
            try:
                if sel_point is not None and abs(float(sel_point) - line) < 0.01:
                    return sel
            except (TypeError, ValueError):
                continue
    return None


def _pick_team_selection(market: dict, team_label: str) -> Optional[str]:
    """Given a 2-way team market (e.g. '1st Run', 'Moneyline', '1st 5 Innings'),
    return the selection.id whose name matches team_label exactly.
    DK selection names are e.g. 'CHI Cubs', 'SD Padres'.
    """
    for sel in market.get('selections', []):
        if (sel.get('name') or '').strip() == team_label:
            return sel.get('id') or sel.get('selectionId')
    return None


# ---------------------------------------------------------------------------
# Per-leg-type resolvers
# ---------------------------------------------------------------------------

def resolve_scores_first(leg, side, event_state, team_names) -> Optional[str]:
    """leg = {'type': 'scores_first'}; side = 'home' | 'away'."""
    market = find_market_by_name(event_state, '1st Run')
    if market is None:
        return None
    return _pick_team_selection(market, team_names[side])


def resolve_wins_period(leg, side, event_state, team_names) -> Optional[str]:
    """leg = {'type': 'wins_period', 'period': 'F3'|'F5'|'F7'|'FG'}."""
    period = leg.get('period')
    market_name = {
        'FG': 'Moneyline',
        'F5': '1st 5 Innings',
        # F3 / F7 not posted as 2-way ML primitives at DK reliably — return None
    }.get(period)
    if market_name is None:
        return None
    market = find_market_by_name(event_state, market_name)
    if market is None:
        return None
    return _pick_team_selection(market, team_names[side])


def resolve_team_total_under(leg, side, event_state, team_names) -> Optional[str]:
    """leg = {'type': 'team_total_under', 'line': 2.5}."""
    line = leg.get('line')
    if line is None:
        return None
    sel = find_team_total_market(event_state, team_names[side], 'under', float(line))
    if sel is None:
        return None
    return sel.get('id') or sel.get('selectionId')


def resolve_team_total_over(leg, side, event_state, team_names) -> Optional[str]:
    """leg = {'type': 'team_total_over', 'line': 4.5}."""
    line = leg.get('line')
    if line is None:
        return None
    sel = find_team_total_market(event_state, team_names[side], 'over', float(line))
    if sel is None:
        return None
    return sel.get('id') or sel.get('selectionId')


# ---------------------------------------------------------------------------
# Registry — mirrors parse_legs.R::TOKEN_REGISTRY structure
# ---------------------------------------------------------------------------

LEG_RESOLVERS = {
    'scores_first':     resolve_scores_first,
    'wins_period':      resolve_wins_period,
    'team_total_under': resolve_team_total_under,
    'team_total_over':  resolve_team_total_over,
    # Reserved (parser knows them; resolver returns None until DK markets
    # are confirmed available):
    'opp_total_under':  lambda *a, **k: None,
    'opp_total_over':   lambda *a, **k: None,
}


def resolve_legs(legs, side, event_state, team_names) -> Optional[list[str]]:
    """Resolve a list of leg specs into DK selection IDs. Returns None if any
    leg fails to resolve (in which case the whole prop is unprice able by DK
    and the caller writes NULL sgp_decimal)."""
    sel_ids = []
    for leg in legs:
        leg_type = leg.get('type')
        resolver = LEG_RESOLVERS.get(leg_type)
        if resolver is None:
            return None
        sel_id = resolver(leg, side, event_state, team_names)
        if sel_id is None:
            return None
        sel_ids.append(sel_id)
    return sel_ids
