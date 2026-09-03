"""Unit tests for KXMLBRFI market parsing (no network)."""
import pytest

from kalshi_rfi.discovery import (drop_doubleheaders, parse_market,
                                  parse_suffix_start_utc)


def _market(ticker, yes_bid=44, yes_ask=46, status="active"):
    return {"ticker": ticker, "yes_bid": yes_bid, "yes_ask": yes_ask,
            "status": status}


class TestSuffixTime:
    def test_et_to_utc(self):
        # 26AUG262105 = 2026-08-26 21:05 ET (EDT, UTC-4) -> 2026-08-27 01:05 UTC
        dt = parse_suffix_start_utc("26AUG262105MINATH")
        assert (dt.year, dt.month, dt.day, dt.hour, dt.minute) == (2026, 8, 27, 1, 5)
        assert dt.tzinfo is None    # naive UTC, GameRef convention

    def test_winter_offset(self):
        # January is EST (UTC-5)
        dt = parse_suffix_start_utc("26JAN151905MINATH")
        assert (dt.day, dt.hour) == (16, 0)

    def test_garbage_returns_none(self):
        assert parse_suffix_start_utc("short") is None
        assert parse_suffix_start_utc("26XXX262105MINATH") is None


class TestParseMarket:
    def test_happy_path(self):
        g = parse_market(_market("KXMLBRFI-26AUG262105MINATH"))
        assert g.home_team == "Athletics"
        assert g.away_team == "Minnesota Twins"
        assert g.suffix == "26AUG262105MINATH"
        assert (g.yes_bid_cents, g.yes_ask_cents) == (44, 46)

    def test_two_letter_home_code(self):
        g = parse_market(_market("KXMLBRFI-26AUG261610PITSD"))
        assert g.home_team == "San Diego Padres"
        assert g.away_team == "Pittsburgh Pirates"

    def test_other_series_rejected(self):
        assert parse_market(_market("KXMLBTOTAL-26AUG262105MINATH-8")) is None

    def test_unknown_team_code_rejected(self):
        assert parse_market(_market("KXMLBRFI-26AUG262105XXXYYY")) is None

    def test_exchange_index_is_read_from_the_payload(self):
        # Sharding (2026-08-24): baseball is 3 today, NFL is 0 — cancels
        # must target the market's own shard, so read it, never assume it.
        m = _market("KXMLBRFI-26AUG262105MINATH")
        assert parse_market(m).exchange_index is None
        assert parse_market({**m, "exchange_index": 3}).exchange_index == 3
        assert parse_market({**m, "exchange_index": "x"}).exchange_index is None


class TestDoubleheaders:
    def test_same_day_pair_dropped(self):
        g1 = parse_market(_market("KXMLBRFI-26AUG261305BALSTL"))
        g2 = parse_market(_market("KXMLBRFI-26AUG261945BALSTL"))
        other = parse_market(_market("KXMLBRFI-26AUG262105MINATH"))
        kept = drop_doubleheaders([g1, g2, other])
        assert kept == [other]

    def test_same_pair_across_days_kept(self):
        g1 = parse_market(_market("KXMLBRFI-26AUG261905HOUNYY"))
        g2 = parse_market(_market("KXMLBRFI-26AUG271905HOUNYY"))
        assert drop_doubleheaders([g1, g2]) == [g1, g2]


if __name__ == "__main__":
    raise SystemExit(pytest.main([__file__, "-v"]))
