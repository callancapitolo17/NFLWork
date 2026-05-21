"""
Wagerzon Scraper Configuration

Each sport declares its wanted markets by WZ's human-readable `Description`
label (the same string visible on the betslip), NOT by a hardcoded `lg=` ID.
At scrape time, `scraper_v2.resolve_leagues()` matches each market against
WZ's live `ActiveLeaguesHelper.aspx` catalog and resolves the current
`IdLeague`. This prevents the silent-mislabel failure mode (BKM's league
503 was labeled "1st 3 Innings" but is actually "2ND HALVES" — same trap
existed here with 25 MLB IDs, most undocumented).

To add a new market: hit /wager/ActiveLeaguesHelper.aspx?WT=0, find the
league's `Description` + `IndexName`, and append an entry to the sport's
`markets` list.

`period` is informational — the parser still detects period per response
league from substring matching the `Description`, so the value here is
mostly documentation. `kind` is also informational today; reserved for
future use when the parser branches more explicitly on market kind.
"""

# Sport configurations
SPORTS = {
    "mlb": {
        "sport_key": "baseball_mlb",
        "table_name": "mlb_odds",
        "wz_index_name": "MLB",
        "markets": [
            # --- Core lines (parser fully supports these) ---
            {"description": "MLB - GAME LINES",                  "period": "fg", "kind": "lines"},
            {"description": "MLB - 1ST 5 FULL INNINGS",          "period": "F5", "kind": "lines"},
            {"description": "MLB - 3 INNINGS LINE",              "period": "F3", "kind": "lines"},
            {"description": "MLB - 7 INNINGS LINE",              "period": "F7", "kind": "lines"},
            {"description": "MLB - ALTERNATE RUNLINES & TOTALS", "period": "fg", "kind": "alts"},

            # --- 3-way moneylines (routed via parse_3way_line) ---
            {"description": "MLB - 1ST INNING WINNER (3 WAY)",            "period": "f1", "kind": "3way"},
            {"description": "MLB - 1ST 5 INN WINNER (3-WAY)",             "period": "F5", "kind": "3way"},
            {"description": "MLB - TEAM WITH THE HIGHEST SCORING INNING", "period": "fg", "kind": "3way"},

            # --- Team totals ---
            {"description": "MLB - TEAM TOTALS",                 "period": "fg",    "kind": "team_total"},
            {"description": "MLB - 1ST HALF TEAM TOTALS",        "period": "Half1", "kind": "team_total"},  # NEW 2026-05-20

            # --- Inning / score-first / odd-even moneylines ---
            {"description": "MLB - TEAM TO SCORE 1ST",                  "period": "fg", "kind": "ml_only"},
            {"description": "MLB - TEAM TO SCORE FIRST WINS THE GAME",  "period": "fg", "kind": "ml_only"},
            {"description": "MLB - SCORE 1ST INNING",                   "period": "f1", "kind": "ml_only"},
            {"description": "MLB - TOTAL RUNS ODD/EVEN",                "period": "fg", "kind": "ml_only"},

            # --- Game-level totals + result props ---
            {"description": "MLB - DOUBLE RESULT",   "period": "fg", "kind": "double_result"},
            {"description": "MLB - WINNING MARGIN",  "period": "fg", "kind": "winning_margin"},
            {"description": "MLB - HITS +RUNS +ERRORS", "period": "fg", "kind": "hre"},
            {"description": "MLB - TOTAL HITS",      "period": "fg", "kind": "total_only"},
            {"description": "MLB - TOTAL BASES",     "period": "fg", "kind": "total_only"},

            # --- Pitcher props (parser handles via idgmtyp=30) ---
            {"description": "MLB - PITCHER STRIKEOUTS",     "period": "fg", "kind": "pitcher_prop"},
            {"description": "MLB - PITCHER HITS ALLOWED",   "period": "fg", "kind": "pitcher_prop"},
            {"description": "MLB - PITCHER WALKS ALLOWED",  "period": "fg", "kind": "pitcher_prop"},
            {"description": "MLB - PITCHER TOTAL OUTS",     "period": "fg", "kind": "pitcher_prop"},

            # --- Player props + exact-result markets ---
            # NOTE (2026-05-20): the parser does NOT currently handle the
            # player-prop response shape (idgmtyp=30 at parent level, with
            # vtm/htm as labels/player names and GameLines keyed by
            # `bml/odds/oddsh/tmname/tmnum`). These markets resolve and get
            # fetched but emit NULL-only rows — no downstream consumer
            # reads them. Listed here for parity with the pre-refactor
            # config; building real player-prop parsing is a separate task.
            # MLB - PLAYER TO HIT 1ST HOME RUN (lg=4038) intentionally
            # omitted — same parser limitation, no point adding more
            # NULL-row noise. Add once real player-prop parsing exists.
            {"description": "MLB - HOME RUN MARKET",                     "period": "fg", "kind": "player_prop_unparsed"},
            {"description": "MLB - PLAYER TO HIT A HOME RUN",            "period": "fg", "kind": "player_prop_unparsed"},
            {"description": "MLB - 1ST PLATE APPEARANCE - EXACT RESULT", "period": "fg", "kind": "player_prop_unparsed"},
            {"description": "MLB - 1ST INNING EXACT HITS",               "period": "fg", "kind": "player_prop_unparsed"},
        ],
    },
    "nba": {
        "sport_key": "basketball_nba",
        "table_name": "nba_odds",
        "wz_index_name": "NBA",
        "markets": [
            {"description": "NBA - GAME LINES",     "period": "fg",    "kind": "lines"},
            {"description": "NBA - 1ST HALF LINES", "period": "Half1", "kind": "lines"},
            {"description": "NBA - TEAM TOTALS",    "period": "fg",    "kind": "team_total"},
            {"description": "NBA - 1H TEAM TOTALS", "period": "Half1", "kind": "team_total"},
            {"description": "NBA - ADJUSTED LINES", "period": "fg",    "kind": "alts"},
        ],
    },
    "nfl": {
        "sport_key": "americanfootball_nfl",
        "table_name": "nfl_odds",
        "wz_index_name": "NFL",
        # UNVERIFIED (off-season at refactor time, 2026-05-20). The previous
        # config hardcoded lg IDs (4029, 430, 432, 433, 2239, 2140, 98,
        # 2767, 124) without comments. None of those except 124 (Super Bowl
        # coin toss) appear in the current catalog because NFL is between
        # seasons. First in-season scrape will WARN-and-skip any patterns
        # below that don't match — the diagnostic prints the actual
        # Descriptions WZ posts. Update this list to match at that point.
        "markets": [
            {"description": "NFL - GAME LINES",       "period": "fg",    "kind": "lines"},
            {"description": "NFL - 1ST HALF LINES",   "period": "Half1", "kind": "lines"},
            {"description": "NFL - 2ND HALF LINES",   "period": "Half2", "kind": "lines"},
            {"description": "NFL - QUARTER LINES",    "period": "Q",     "kind": "lines"},
            {"description": "NFL - TEAM TOTALS",      "period": "fg",    "kind": "team_total"},
            {"description": "NFL - ALTERNATE LINES",  "period": "fg",    "kind": "alts"},
        ],
    },
    "cbb": {
        "sport_key": "basketball_ncaab",
        "table_name": "cbb_odds",
        "wz_index_name": "CBB",
        # UNVERIFIED (off-season at refactor time, 2026-05-20). Previous
        # config used lg=43,403,45 — none currently in the catalog. First
        # in-season scrape will WARN-and-skip; diagnostic names what
        # IS available. Race-to-X props remain hardcoded below because
        # they use a different URL parameter shape.
        "markets": [
            {"description": "CBB - GAME LINES",     "period": "fg",    "kind": "lines"},
            {"description": "CBB - 1ST HALF LINES", "period": "Half1", "kind": "lines"},
        ],
        # Race-to-X props are fetched as standalone leagues via separate
        # URL params, not as derivatives of the main games. Keep the
        # hardcoded `lg=` IDs here because they're stable, niche markets
        # that the parser handles via lg-specific logic at
        # scraper_v2.py:591+.
        "prop_params": {
            "race_to_10": "lg=1852",
            "race_to_20": "lg=1089",
            "race_to_40": "lg=4139",
        },
    },
    "college_baseball": {
        "sport_key": "baseball_ncaa",
        "table_name": "college_baseball_odds",
        # NOTE: WZ uses IndexName "NCAA BASEBALL" for college baseball
        # (the only sport in this codebase where IndexName differs from
        # the wz-config sport key).
        "wz_index_name": "NCAA BASEBALL",
        "markets": [
            {"description": "COLLEGE BASEBALL",     "period": "fg", "kind": "lines"},
            # NCAA GAME + 1ST 5 INNINGS markets were in previous config
            # (lg=1554, 4321) but aren't currently in the catalog.
            # Late-season slate is thin; these may reappear next season.
            {"description": "NCAA GAME",            "period": "fg", "kind": "lines"},
            {"description": "1ST 5 INNINGS",        "period": "F5", "kind": "lines"},
        ],
    },
}

# Wagerzon backend URLs
WAGERZON_BASE_URL = "https://backend.wagerzon.com"
WAGERZON_HELPER_URL = f"{WAGERZON_BASE_URL}/wager/NewScheduleHelper.aspx"
WAGERZON_SCHEDULE_URL = f"{WAGERZON_BASE_URL}/wager/NewSchedule.aspx?WT=0&"


def get_sport_url(sport: str) -> str:
    """Deprecated — league IDs are now resolved at runtime.

    Pre-refactor this built a static URL from `url_params`. Now that
    league IDs come from `ActiveLeaguesHelper` at scrape time, this
    helper has nothing useful to return. Use
    `scraper_v2.fetch_odds_json()` instead.

    Raises to surface any external caller that still depends on it.
    """
    if sport not in SPORTS:
        raise ValueError(f"Unknown sport: {sport}. Available: {list(SPORTS.keys())}")
    raise RuntimeError(
        "get_sport_url() is deprecated — league IDs are resolved at "
        "runtime via ActiveLeaguesHelper. Use "
        "scraper_v2.fetch_odds_json() instead."
    )


def get_sport_config(sport: str) -> dict:
    """Get configuration for a sport."""
    if sport not in SPORTS:
        raise ValueError(f"Unknown sport: {sport}. Available: {list(SPORTS.keys())}")
    return SPORTS[sport]
