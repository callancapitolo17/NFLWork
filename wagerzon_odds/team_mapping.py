"""
Team name mappings for normalizing Wagerzon team names to standard names
"""

# Wagerzon uses abbreviated city + team name (e.g., "SEA SEAHAWKS")
# Map to full names matching The Odds API format

NFL_TEAMS = {
    # AFC East
    "BUF BILLS": "Buffalo Bills",
    "BUFFALO BILLS": "Buffalo Bills",
    "MIA DOLPHINS": "Miami Dolphins",
    "MIAMI DOLPHINS": "Miami Dolphins",
    "NE PATRIOTS": "New England Patriots",
    "NEW ENGLAND PATRIOTS": "New England Patriots",
    "NY JETS": "New York Jets",
    "NEW YORK JETS": "New York Jets",
    "NYJ JETS": "New York Jets",

    # AFC North
    "BAL RAVENS": "Baltimore Ravens",
    "BALTIMORE RAVENS": "Baltimore Ravens",
    "CIN BENGALS": "Cincinnati Bengals",
    "CINCINNATI BENGALS": "Cincinnati Bengals",
    "CLE BROWNS": "Cleveland Browns",
    "CLEVELAND BROWNS": "Cleveland Browns",
    "PIT STEELERS": "Pittsburgh Steelers",
    "PITTSBURGH STEELERS": "Pittsburgh Steelers",

    # AFC South
    "HOU TEXANS": "Houston Texans",
    "HOUSTON TEXANS": "Houston Texans",
    "IND COLTS": "Indianapolis Colts",
    "INDIANAPOLIS COLTS": "Indianapolis Colts",
    "JAC JAGUARS": "Jacksonville Jaguars",
    "JAX JAGUARS": "Jacksonville Jaguars",
    "JACKSONVILLE JAGUARS": "Jacksonville Jaguars",
    "TEN TITANS": "Tennessee Titans",
    "TENNESSEE TITANS": "Tennessee Titans",

    # AFC West
    "DEN BRONCOS": "Denver Broncos",
    "DENVER BRONCOS": "Denver Broncos",
    "KC CHIEFS": "Kansas City Chiefs",
    "KANSAS CITY CHIEFS": "Kansas City Chiefs",
    "LV RAIDERS": "Las Vegas Raiders",
    "LAS VEGAS RAIDERS": "Las Vegas Raiders",
    "LAC CHARGERS": "Los Angeles Chargers",
    "LA CHARGERS": "Los Angeles Chargers",
    "LOS ANGELES CHARGERS": "Los Angeles Chargers",

    # NFC East
    "DAL COWBOYS": "Dallas Cowboys",
    "DALLAS COWBOYS": "Dallas Cowboys",
    "NYG GIANTS": "New York Giants",
    "NY GIANTS": "New York Giants",
    "NEW YORK GIANTS": "New York Giants",
    "PHI EAGLES": "Philadelphia Eagles",
    "PHILADELPHIA EAGLES": "Philadelphia Eagles",
    "WAS COMMANDERS": "Washington Commanders",
    "WASHINGTON COMMANDERS": "Washington Commanders",
    "WSH COMMANDERS": "Washington Commanders",

    # NFC North
    "CHI BEARS": "Chicago Bears",
    "CHICAGO BEARS": "Chicago Bears",
    "DET LIONS": "Detroit Lions",
    "DETROIT LIONS": "Detroit Lions",
    "GB PACKERS": "Green Bay Packers",
    "GREEN BAY PACKERS": "Green Bay Packers",
    "MIN VIKINGS": "Minnesota Vikings",
    "MINNESOTA VIKINGS": "Minnesota Vikings",

    # NFC South
    "ATL FALCONS": "Atlanta Falcons",
    "ATLANTA FALCONS": "Atlanta Falcons",
    "CAR PANTHERS": "Carolina Panthers",
    "CAROLINA PANTHERS": "Carolina Panthers",
    "NO SAINTS": "New Orleans Saints",
    "NEW ORLEANS SAINTS": "New Orleans Saints",
    "TB BUCCANEERS": "Tampa Bay Buccaneers",
    "TAMPA BAY BUCCANEERS": "Tampa Bay Buccaneers",

    # NFC West
    "ARI CARDINALS": "Arizona Cardinals",
    "ARIZONA CARDINALS": "Arizona Cardinals",
    "LA RAMS": "Los Angeles Rams",
    "LAR RAMS": "Los Angeles Rams",
    "LOS ANGELES RAMS": "Los Angeles Rams",
    "SF 49ERS": "San Francisco 49ers",
    "SAN FRANCISCO 49ERS": "San Francisco 49ers",
    "SEA SEAHAWKS": "Seattle Seahawks",
    "SEATTLE SEAHAWKS": "Seattle Seahawks",
}

# NBA teams (placeholder for future)
NBA_TEAMS = {
    # Add as needed
}

# College Basketball - handled differently due to variety
CBB_TEAMS = {
    # Would need a more comprehensive mapping or fuzzy matching
}


def normalize_team_name(team: str, sport: str = "nfl") -> str:
    """
    Normalize team name from Wagerzon format to standard format.

    Args:
        team: Raw team name from Wagerzon (e.g., "SEA SEAHAWKS")
        sport: Sport key (nfl, nba, cbb)

    Returns:
        Normalized team name (e.g., "Seattle Seahawks")
    """
    if not team:
        return team

    # Clean up the team name
    team_upper = team.upper().strip()

    # Get the appropriate mapping
    if sport == "nfl":
        mapping = NFL_TEAMS
    elif sport == "nba":
        mapping = NBA_TEAMS
    else:
        mapping = {}

    # Look up the team
    if team_upper in mapping:
        return mapping[team_upper]

    # If not found, return original (might need to be added to mapping)
    print(f"Warning: Unknown team name '{team}' for sport '{sport}'")
    return team
