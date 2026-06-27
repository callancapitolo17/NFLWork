from unabated_edge.sports.soccer import Soccer
ADAPTERS = [Soccer()]
def league_prefixes(): return {a.league_prefix for a in ADAPTERS}
def by_league(prefix): return next((a for a in ADAPTERS if a.league_prefix == prefix), None)
