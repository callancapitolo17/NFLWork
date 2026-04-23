"""End-to-end test: unmapped player + unmapped market land in quarantine, NOT draft_odds."""
import pytest
from datetime import datetime
from nfl_draft.lib import db as db_module
from nfl_draft.lib import seed
from nfl_draft.lib.quarantine import write_or_quarantine
from nfl_draft.scrapers._base import OddsRow


@pytest.fixture
def seeded(monkeypatch, tmp_path):
    monkeypatch.setattr(db_module, "DB_PATH", tmp_path / "test.duckdb")
    db_module.init_schema()
    seed.run()


def test_unmapped_player_lands_in_quarantine(seeded):
    rows = [OddsRow(
        book="draftkings", book_label="First QB", book_subject="Made Up Player",
        american_odds=200, fetched_at=datetime.now(),
    )]
    write_or_quarantine(rows)
    with db_module.read_connection() as con:
        odds_count = con.execute("SELECT COUNT(*) FROM draft_odds").fetchone()[0]
        quarantine_count = con.execute("SELECT COUNT(*) FROM draft_odds_unmapped").fetchone()[0]
    assert odds_count == 0
    assert quarantine_count == 1


def test_mapped_row_lands_in_draft_odds(seeded):
    from nfl_draft.lib.db import write_connection
    with write_connection() as con:
        # Insert a market row first to satisfy any FK
        con.execute("INSERT INTO draft_markets (market_id, market_type, subject_player, position) VALUES ('first_qb_cam-ward', 'first_at_position', 'Cam Ward', 'QB')")
        con.execute("INSERT INTO market_map VALUES ('draftkings', 'First QB', 'Cam Ward', 'first_qb_cam-ward')")
    rows = [OddsRow(
        book="draftkings", book_label="First QB", book_subject="Cam Ward",
        american_odds=200, fetched_at=datetime.now(),
    )]
    write_or_quarantine(rows)
    with db_module.read_connection() as con:
        odds_count = con.execute("SELECT COUNT(*) FROM draft_odds").fetchone()[0]
        quarantine_count = con.execute("SELECT COUNT(*) FROM draft_odds_unmapped").fetchone()[0]
    assert odds_count == 1
    assert quarantine_count == 0


def test_write_or_quarantine_perf_1000_rows(seeded):
    """1000 rows should complete in < 2s (pre-fetch lookup avoids per-row read)."""
    import time
    rows = [
        OddsRow(
            book="draftkings",
            book_label=f"label_{i}",
            book_subject=f"player_{i}",
            american_odds=100,
            fetched_at=datetime.now(),
        )
        for i in range(1000)
    ]
    t0 = time.time()
    write_or_quarantine(rows)
    elapsed = time.time() - t0
    assert elapsed < 2.0, f"write_or_quarantine took {elapsed:.2f}s for 1000 rows (expected < 2s)"


def _seed_market_map_for_first_wr(con):
    """Helper: seed draft_markets + market_map for a 3-candidate first_wr outright.

    Uses DK as the only book so we can observe devig behavior in isolation.
    """
    con.execute("INSERT INTO draft_markets (market_id, market_type, subject_player, position) "
                "VALUES ('first_wr_hunter', 'first_at_position', 'Hunter', 'WR')")
    con.execute("INSERT INTO draft_markets (market_id, market_type, subject_player, position) "
                "VALUES ('first_wr_tyson', 'first_at_position', 'Tyson', 'WR')")
    con.execute("INSERT INTO draft_markets (market_id, market_type, subject_player, position) "
                "VALUES ('first_wr_tate', 'first_at_position', 'Tate', 'WR')")
    con.execute("INSERT INTO market_map VALUES ('draftkings', '1st WR', 'Hunter', 'first_wr_hunter')")
    con.execute("INSERT INTO market_map VALUES ('draftkings', '1st WR', 'Tyson',  'first_wr_tyson')")
    con.execute("INSERT INTO market_map VALUES ('draftkings', '1st WR', 'Tate',   'first_wr_tate')")


def test_write_or_quarantine_devigs_outright_group(seeded):
    """3-candidate first_wr outright at DK gets proportional-devigged."""
    from nfl_draft.lib.db import write_connection, read_connection
    with write_connection() as con:
        _seed_market_map_for_first_wr(con)
    now = datetime.now()
    rows = [
        OddsRow(book="draftkings", book_label="1st WR", book_subject="Hunter",
                american_odds=+110, fetched_at=now, market_group="first_at_position"),
        OddsRow(book="draftkings", book_label="1st WR", book_subject="Tyson",
                american_odds=+200, fetched_at=now, market_group="first_at_position"),
        OddsRow(book="draftkings", book_label="1st WR", book_subject="Tate",
                american_odds=+400, fetched_at=now, market_group="first_at_position"),
    ]
    write_or_quarantine(rows)
    with read_connection() as con:
        result = con.execute(
            "SELECT market_id, implied_prob, devig_prob "
            "FROM draft_odds ORDER BY market_id"
        ).fetchall()
    # Sum of raw implieds: 100/210 + 100/300 + 100/500 = 0.4762 + 0.3333 + 0.2000 = 1.0095
    # Devigged: each raw / 1.0095
    assert len(result) == 3
    total_devig = sum(r[2] for r in result)
    assert abs(total_devig - 1.0) < 1e-6, f"devig should sum to 1.0, got {total_devig}"
    # Each implied_prob should NOT equal devig_prob (would indicate no devig ran)
    for market_id, implied, devig in result:
        assert abs(implied - devig) > 1e-4, f"{market_id}: implied={implied} equals devig={devig}"


def test_write_or_quarantine_single_row_bucket_uses_implied(seeded):
    """A prop with no sibling (group=None) gets devig_prob = implied_prob."""
    from nfl_draft.lib.db import write_connection, read_connection
    with write_connection() as con:
        con.execute("INSERT INTO draft_markets (market_id, market_type) "
                    "VALUES ('prop_over_9p5_total_qbs_drafted_round_1', 'prop')")
        con.execute("INSERT INTO market_map VALUES "
                    "('draftkings', 'Total QBs Round 1', 'Over 9.5', 'prop_over_9p5_total_qbs_drafted_round_1')")
    rows = [OddsRow(
        book="draftkings", book_label="Total QBs Round 1", book_subject="Over 9.5",
        american_odds=-120, fetched_at=datetime.now(),
        market_group="prop_first_round_total_ou",
    )]
    write_or_quarantine(rows)
    with read_connection() as con:
        implied, devig = con.execute(
            "SELECT implied_prob, devig_prob FROM draft_odds"
        ).fetchone()
    # -120 -> implied = 120/220 = 0.5454...
    assert abs(implied - 120/220) < 1e-6
    assert abs(devig - implied) < 1e-9, "1-row group should pass implied through as devig"


def test_write_or_quarantine_respects_kalshi_pre_set_devig(seeded):
    """Kalshi rows arrive with devig_prob already set to mid; must not be overwritten."""
    from nfl_draft.lib.db import write_connection, read_connection
    with write_connection() as con:
        con.execute("INSERT INTO draft_markets (market_id, market_type, subject_player, position) "
                    "VALUES ('first_wr_hunter', 'first_at_position', 'Hunter', 'WR')")
        con.execute("INSERT INTO market_map VALUES "
                    "('kalshi', 'KXNFLDRAFT-FIRST-WR', 'Hunter', 'first_wr_hunter')")
    rows = [OddsRow(
        book="kalshi", book_label="KXNFLDRAFT-FIRST-WR", book_subject="Hunter",
        american_odds=+138, fetched_at=datetime.now(),
        market_group="",  # Kalshi doesn't set market_group
        implied_prob=0.42, devig_prob=0.40,  # take vs mid
    )]
    write_or_quarantine(rows)
    with read_connection() as con:
        implied, devig = con.execute(
            "SELECT implied_prob, devig_prob FROM draft_odds"
        ).fetchone()
    assert abs(implied - 0.42) < 1e-9
    assert abs(devig - 0.40) < 1e-9, "Kalshi pre-set devig_prob must be preserved"


def _seed_market_map_for_top_10(con, players):
    """Seed draft_markets + DK market_map for a list of top_10 candidates."""
    for p in players:
        slug_ = p.lower().replace(" ", "-")
        con.execute(
            "INSERT INTO draft_markets (market_id, market_type, subject_player) "
            "VALUES (?, 'top_n_range', ?)",
            [f"top_10_{slug_}", p],
        )
        con.execute(
            "INSERT INTO market_map VALUES ('draftkings', 'Top 10 Pick', ?, ?)",
            [p, f"top_10_{slug_}"],
        )


def test_write_or_quarantine_pool_devigs_top_n_range(seeded):
    """top_10_range rows at DK get pool-normalized to sum to 10, not 1.

    15 candidates at -150 (implied 0.60 each) -> sum(implieds)=9.0, which
    clears the 0.9*10=9.0 coverage floor, so devig_pool(odds, 10) fires
    and each fair scales by 10/9.
    """
    from nfl_draft.lib.db import write_connection, read_connection
    players = [f"Player {i:02d}" for i in range(15)]
    with write_connection() as con:
        _seed_market_map_for_top_10(con, players)
    now = datetime.now()
    rows = [
        OddsRow(book="draftkings", book_label="Top 10 Pick", book_subject=p,
                american_odds=-150, fetched_at=now, market_group="top_10_range")
        for p in players
    ]
    write_or_quarantine(rows)
    with read_connection() as con:
        result = con.execute(
            "SELECT implied_prob, devig_prob FROM draft_odds"
        ).fetchall()
    assert len(result) == 15
    total_devig = sum(r[1] for r in result)
    assert abs(total_devig - 10.0) < 1e-6, f"pool devig should sum to 10.0, got {total_devig}"
    # Each fair = raw_implied * 10/9 = 0.60 * 10/9 ≈ 0.6667
    for implied, devig in result:
        assert abs(implied - 0.60) < 1e-6
        assert abs(devig - 0.60 * 10 / 9) < 1e-6


def test_write_or_quarantine_pool_devig_sparse_coverage_falls_back(seeded):
    """Sparse top_10 board (coverage < 0.9*N) must NOT pool-normalize.

    3 candidates at -150 sum to 1.8 in implied; pool-normalizing to 10
    would inflate each fair to ~3.33 (nonsense). Guardrail must pass raw
    implieds through instead.
    """
    from nfl_draft.lib.db import write_connection, read_connection
    players = ["Player A", "Player B", "Player C"]
    with write_connection() as con:
        _seed_market_map_for_top_10(con, players)
    now = datetime.now()
    rows = [
        OddsRow(book="draftkings", book_label="Top 10 Pick", book_subject=p,
                american_odds=-150, fetched_at=now, market_group="top_10_range")
        for p in players
    ]
    write_or_quarantine(rows)
    with read_connection() as con:
        result = con.execute(
            "SELECT implied_prob, devig_prob FROM draft_odds"
        ).fetchall()
    assert len(result) == 3
    # Raw implied is 0.60 for -150; guardrail kicks in (1.8 < 9.0), so devig
    # = implied. No inflation.
    for implied, devig in result:
        assert abs(implied - 0.60) < 1e-6
        assert abs(devig - implied) < 1e-9, \
            f"sparse top_N bucket should pass implied through, got devig={devig}"


def test_write_or_quarantine_groups_only_within_book(seeded):
    """DK's first_wr devig must not be polluted by another book's first_wr rows."""
    from nfl_draft.lib.db import write_connection, read_connection
    with write_connection() as con:
        _seed_market_map_for_first_wr(con)
        # Add bookmaker as second book on the same outright
        con.execute("INSERT INTO market_map VALUES "
                    "('bookmaker', '1st WR', 'Hunter', 'first_wr_hunter')")
        con.execute("INSERT INTO market_map VALUES "
                    "('bookmaker', '1st WR', 'Tyson',  'first_wr_tyson')")
    now = datetime.now()
    rows = [
        # DK: 2-candidate subset
        OddsRow(book="draftkings", book_label="1st WR", book_subject="Hunter",
                american_odds=+100, fetched_at=now, market_group="first_at_position"),
        OddsRow(book="draftkings", book_label="1st WR", book_subject="Tyson",
                american_odds=+100, fetched_at=now, market_group="first_at_position"),
        # Bookmaker: 2-candidate subset with different numbers
        OddsRow(book="bookmaker", book_label="1st WR", book_subject="Hunter",
                american_odds=-150, fetched_at=now, market_group="first_at_position"),
        OddsRow(book="bookmaker", book_label="1st WR", book_subject="Tyson",
                american_odds=+200, fetched_at=now, market_group="first_at_position"),
    ]
    write_or_quarantine(rows)
    with read_connection() as con:
        dk_sum = con.execute(
            "SELECT SUM(devig_prob) FROM draft_odds WHERE book='draftkings'"
        ).fetchone()[0]
        bm_sum = con.execute(
            "SELECT SUM(devig_prob) FROM draft_odds WHERE book='bookmaker'"
        ).fetchone()[0]
    assert abs(dk_sum - 1.0) < 1e-6, f"DK devig should sum to 1.0 within book, got {dk_sum}"
    assert abs(bm_sum - 1.0) < 1e-6, f"Bookmaker devig should sum to 1.0 within book, got {bm_sum}"
