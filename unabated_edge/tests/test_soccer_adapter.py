import datetime
from unabated_edge import feed
from unabated_edge.sports.soccer import Soccer


def _state_with_bt3(eid=9, ms=7, over_price=115, under_price=-140, line=2.5):
    """FeedState with bt3 total over/under lines for one event."""
    return feed.parse_snapshot({
        "marketSources": [{"id": ms, "name": "S"}],
        "teams": {"1": {"name": "Colombia"}, "2": {"name": "Ghana"}},
        "gameOddsEvents": {
            "lg21:pt1:pregame": [{
                "eventId": eid,
                "eventStart": "2026-06-22T17:00:00+00:00",
                "eventTeams": {"1": {"id": 1}, "0": {"id": 2}},
                "gameOddsMarketSourcesLines": {
                    f"si0:ms{ms}:an0": {"bt3": {"price": over_price, "points": line}},  # Over
                    f"si1:ms{ms}:an0": {"bt3": {"price": under_price, "points": line}},  # Under
                }
            }]
        }
    }, {"lg21"})


def test_canon_alias():
    assert Soccer().canon_team("Korea Republic") == "South Korea"


def test_canon_alias_expanded_collapses_both_spellings():
    s = Soccer()
    assert s.canon_team("Côte d'Ivoire") == "Ivory Coast"
    assert s.canon_team("Ivory Coast") == "Ivory Coast"      # plain form falls through to itself
    assert s.canon_team("Czechia") == s.canon_team("Czech Republic")
    assert s.canon_team("Türkiye") == "Turkey"


def test_anchor_totals_returns_ladder():
    s = Soccer()
    st = _state_with_bt3(eid=9, ms=7, over_price=115, under_price=-140, line=2.5)
    ladder = s._anchor_totals(st, 9)
    assert set(ladder) == {2.5}
    assert 0 < ladder[2.5] < 1


def test_anchor_totals_devig_reflects_prices():
    """Devig must produce a sane prob that reflects the input prices, not a tautology.
    Under is favored (-140) vs Over (+115), so devigged P(over) < 0.5 and below the
    raw vig-inflated over-implied (vig was removed)."""
    s = Soccer()
    st = _state_with_bt3(eid=9, ms=7, over_price=115, under_price=-140, line=2.5)
    p_over = s._anchor_totals(st, 9)[2.5]
    raw_over_implied = 100 / 215            # american_to_prob(+115) ≈ 0.465 (still has vig)
    assert p_over < 0.5                     # -140 under is favored
    assert p_over < raw_over_implied        # devig shrank the inflated implied
    assert 0.40 < p_over < 0.47             # sane devigged range


def test_anchor_totals_ladder_from_alternate_lines():
    """Main line 2.75 (quarter-goal, no Kalshi rung) + alternateLines at 2.5/3.5 on
    BOTH sides -> ladder includes the half-goal rungs; price_event matches them."""
    s = Soccer()
    st = feed.parse_snapshot({
        "marketSources": [{"id": 7, "name": "S"}],
        "teams": {"1": {"name": "England"}, "2": {"name": "Norway"}},
        "gameOddsEvents": {"lg21:pt1:pregame": [{
            "eventId": 9, "eventStart": "2026-07-11T21:00:00+00:00",
            "eventTeams": {"1": {"id": 1}, "0": {"id": 2}},
            "gameOddsMarketSourcesLines": {
                "si0:ms7:an0": {"bt3": {"price": -107, "points": 2.75, "alternateLines": [
                    {"points": 2.5, "price": -140}, {"points": 3.5, "price": 254}]}},
                "si1:ms7:an0": {"bt3": {"price": -105, "points": 2.75, "alternateLines": [
                    {"points": 2.5, "price": 120}, {"points": 3.5, "price": -347}]}},
            }}]}}, {"lg21"})
    ladder = s._anchor_totals(st, 9)
    assert set(ladder) == {2.5, 2.75, 3.5}
    assert ladder[2.5] > ladder[2.75] > ladder[3.5]     # P(over) monotone falling in line
    kev = {"title": "England vs Norway: Regulation Time Total Goals", "markets": [
        {"ticker": "T-O25", "strike_type": "greater", "floor_strike": 2.5},
        {"ticker": "T-O35", "strike_type": "greater", "floor_strike": 3.5},
        {"ticker": "T-O45", "strike_type": "greater", "floor_strike": 4.5},   # no anchor line -> skipped
    ]}
    cands = s.price_event(st, st.events[9], kev)
    assert len(cands) == 4                              # 2 rungs x (yes+no); 4.5 skipped
    assert {c.market_ticker for c in cands} == {"T-O25", "T-O35"}


def test_side_prices_main_wins_over_alt_collision():
    """An alternateLines entry at the SAME points as the main quote must not
    overwrite the main (raw book) price."""
    s = Soccer()
    prices = s._side_prices({"price": -140, "points": 2.5,
                             "alternateLines": [{"points": 2.5, "price": -999},
                                                {"points": 3.5, "price": 254}]})
    assert prices[2.5] == (-140, False)          # main survived, marked non-alt
    assert prices[3.5] == (254, True)            # genuine alt marked alt


def test_candidates_carry_provenance_meta():
    """Every candidate carries book/alt/overround so the research firehose can
    audit whether alt rungs behave like real quotes or derived numbers."""
    s = Soccer()
    st = _state_with_bt3(eid=9, ms=7, over_price=115, under_price=-140, line=2.5)
    kev = {"title": "Colombia vs Ghana: Regulation Time Total Goals",
           "markets": [{"ticker": "T-O25", "strike_type": "greater", "floor_strike": 2.5}]}
    cands = s.price_event(st, st.events[9], kev)
    for c in cands:
        assert c.meta["book"] == 7
        assert c.meta["alt"] is False            # built from the main quote
        assert 1.0 < c.meta["overround"] < 1.2   # raw implied sum carries the vig


def test_anchor_totals_same_book_only():
    """Book 7 incomplete (over side only) -> fall through to book 6's complete
    ladder; never mix book 7's over with book 6's under."""
    s = Soccer()
    st = feed.parse_snapshot({
        "marketSources": [{"id": 7, "name": "S7"}, {"id": 6, "name": "S6"}],
        "teams": {"1": {"name": "A"}, "2": {"name": "B"}},
        "gameOddsEvents": {"lg21:pt1:pregame": [{
            "eventId": 9, "eventStart": "2026-07-11T21:00:00+00:00",
            "eventTeams": {"1": {"id": 1}, "0": {"id": 2}},
            "gameOddsMarketSourcesLines": {
                "si0:ms7:an0": {"bt3": {"price": 115, "points": 2.5}},   # book 7: over only
                "si0:ms6:an0": {"bt3": {"price": 100, "points": 2.5}},
                "si1:ms6:an0": {"bt3": {"price": -120, "points": 2.5}},  # book 6: complete
            }}]}}, {"lg21"})
    ladder = s._anchor_totals(st, 9)
    assert set(ladder) == {2.5}
    # devig of book 6's 100/-120, NOT book 7's 115 mixed with book 6's -120:
    from unabated_edge import pricing
    expected, _ = pricing.devig([pricing.american_to_prob(100), pricing.american_to_prob(-120)])
    assert abs(ladder[2.5] - expected) < 1e-9


def test_candidate_side_label_binding():
    """Pin the money-critical binding: yes→over→P(over), no→under→P(under).
    A yes↔under swap (or label/prob mix-up) must fail this."""
    s = Soccer()
    st = _state_with_bt3(eid=9, ms=7, over_price=115, under_price=-140, line=2.5)  # under favored
    p_over = s._anchor_totals(st, 9)[2.5]
    kev = {
        "title": "Colombia vs Ghana: Regulation Time Total Goals",
        "markets": [{"ticker": "T-O25", "strike_type": "greater", "floor_strike": 2.5,
                     "yes_sub_title": "Reg Time: Over 2.5 goals scored"}],
    }
    cands = s.price_event(st, st.events[9], kev)
    yes = next(c for c in cands if c.side == "yes")
    no = next(c for c in cands if c.side == "no")
    assert yes.label == "over_2.5" and abs(yes.fair_prob - p_over) < 1e-9
    assert no.label == "under_2.5" and abs(no.fair_prob - (1 - p_over)) < 1e-5  # code rounds to 6dp
    assert yes.fair_prob < no.fair_prob     # over (underdog) < under (favored)


def test_price_event_returns_yes_and_no_candidates():
    s = Soccer()
    st = _state_with_bt3(eid=9, ms=7, over_price=115, under_price=-140, line=2.5)
    ev_meta = st.events[9]
    kalshi_event = {
        "title": "Colombia vs Ghana: Regulation Time Total Goals",
        "markets": [
            {"ticker": "T-O25", "strike_type": "greater", "floor_strike": 2.5,
             "yes_sub_title": "Reg Time: Over 2.5 goals scored"},
        ]
    }
    candidates = s.price_event(st, ev_meta, kalshi_event)
    assert len(candidates) == 2
    assert all(c.market_ticker == "T-O25" for c in candidates)
    sides = {c.side for c in candidates}
    assert sides == {"yes", "no"}
    labels = {c.label for c in candidates}
    assert "over_2.5" in labels
    assert "under_2.5" in labels
    # fair probs sum to 1
    total_prob = sum(c.fair_prob for c in candidates)
    assert abs(total_prob - 1.0) < 1e-5


def test_event_teams_parsed_from_title():
    s = Soccer()
    result = s.event_teams({"title": "South Africa vs Canada: Regulation Time Total Goals"})
    assert result == frozenset({s.canon_team("South Africa"), s.canon_team("Canada")})


def test_event_teams_with_alias():
    """Aliases are applied during event_teams parsing."""
    s = Soccer()
    result = s.event_teams({"title": "Korea Republic vs Mexico: Regulation Time Total Goals"})
    assert "South Korea" in result
    assert "Mexico" in result


def test_fail_closed_no_matching_floor_strike():
    """When no Kalshi market matches the anchor line, price_event returns []."""
    s = Soccer()
    st = _state_with_bt3(eid=9, ms=7, over_price=115, under_price=-140, line=2.5)
    ev_meta = st.events[9]
    # Market at 1.5, but anchor quoted 2.5 → no match
    kalshi_event = {
        "title": "Colombia vs Ghana: Regulation Time Total Goals",
        "markets": [
            {"ticker": "T-O15", "strike_type": "greater", "floor_strike": 1.5,
             "yes_sub_title": "Reg Time: Over 1.5 goals scored"},
        ]
    }
    candidates = s.price_event(st, ev_meta, kalshi_event)
    assert candidates == []


def test_fail_closed_anchor_missing_under():
    """When the Under side is missing, _anchor_totals returns {} → price_event []."""
    s = Soccer()
    st = feed.parse_snapshot({
        "marketSources": [{"id": 7, "name": "S"}],
        "teams": {"1": {"name": "Colombia"}, "2": {"name": "Ghana"}},
        "gameOddsEvents": {
            "lg21:pt1:pregame": [{
                "eventId": 9,
                "eventStart": "2026-06-22T17:00:00+00:00",
                "eventTeams": {"1": {"id": 1}, "0": {"id": 2}},
                "gameOddsMarketSourcesLines": {
                    "si0:ms7:an0": {"bt3": {"price": 115, "points": 2.5}},
                    # si1 (Under) deliberately absent
                }
            }]
        }
    }, {"lg21"})
    assert s._anchor_totals(st, 9) == {}
    ev_meta = st.events[9]
    kalshi_event = {
        "title": "Colombia vs Ghana: Regulation Time Total Goals",
        "markets": [{"ticker": "T-O25", "strike_type": "greater", "floor_strike": 2.5,
                     "yes_sub_title": "Reg Time: Over 2.5 goals scored"}]
    }
    candidates = s.price_event(st, ev_meta, kalshi_event)
    assert candidates == []


def test_fair_ladder_exposes_anchor_ladder():
    from unabated_edge.tests.test_runner_tick import _state, _NOW  # reuse fixture helpers
    st = _state()
    em = next(iter(st.events.values()))
    ladder = Soccer().fair_ladder(st, em)
    assert 2.5 in ladder and 0 < ladder[2.5]["p_over"] < 1


def test_fair_ladder_none_when_no_anchor():
    from unabated_edge.tests.test_runner_tick import _state
    st = _state()
    st.lines.clear()
    em = next(iter(st.events.values()))
    assert Soccer().fair_ladder(st, em) is None


def test_ladder_rung_carries_older_modified_on():
    """Each rung's modified_on is the OLDER of the two injected side modifiedOn
    strings (the staleness gate must trust the stalest half of the rung)."""
    s = Soccer()
    st = feed.parse_snapshot({
        "marketSources": [{"id": 7, "name": "S"}],
        "teams": {"1": {"name": "Colombia"}, "2": {"name": "Ghana"}},
        "gameOddsEvents": {
            "lg21:pt1:pregame": [{
                "eventId": 9,
                "eventStart": "2026-06-22T17:00:00+00:00",
                "eventTeams": {"1": {"id": 1}, "0": {"id": 2}},
                "gameOddsMarketSourcesLines": {
                    "si0:ms7:an0": {"bt3": {"price": 115, "points": 2.5,
                                            "modifiedOn": "2026-07-14T10:00:00"}},
                    "si1:ms7:an0": {"bt3": {"price": -140, "points": 2.5,
                                            "modifiedOn": "2026-07-14T09:00:00"}},
                }
            }]
        }
    }, {"lg21"})
    ladder = s._anchor_ladder(st, 9)
    assert ladder[2.5]["modified_on"] == "2026-07-14T09:00:00"


def test_older_mo_helper():
    from unabated_edge.sports.soccer import _older_mo
    assert _older_mo("2026-07-14T10:00:00", "2026-07-14T09:00:00") == "2026-07-14T09:00:00"
    assert _older_mo(None, "x") == "x"
    assert _older_mo(None, None) is None


# ---------- overround feed-integrity gate (issue #73) ----------

def test_crossed_main_rung_rejected():
    """MAIN rung whose raw sides sum < 1 (one side stale mid-refresh) must be
    refused before devig — marked rejected, never priced. Fails on pre-#73 main:
    devig used to clip-and-solve the crossed pair into a plausible p_over."""
    s = Soccer()
    st = _state_with_bt3(eid=9, ms=7, over_price=150, under_price=150, line=2.5)  # 0.4+0.4=0.8
    ladder = s._anchor_ladder(st, 9)
    assert ladder == {} or ladder.get(2.5, {}).get("reject") == "anchor_crossed"
    assert s._anchor_totals(st, 9) == {}
    kev = {"title": "Colombia vs Ghana: Regulation Time Total Goals",
           "markets": [{"ticker": "T-O25", "strike_type": "greater", "floor_strike": 2.5}]}
    assert s.price_event(st, st.events[9], kev) == []


def test_blown_vig_main_rung_rejected():
    """MAIN rung with implausible vig (sum 1.33) is refused with the distinct
    band reason, not silently devigged."""
    s = Soccer()
    st = _state_with_bt3(eid=9, ms=7, over_price=-200, under_price=-200, line=2.5)  # 0.667*2
    assert s._anchor_totals(st, 9) == {}
    kev = {"title": "Colombia vs Ghana: Regulation Time Total Goals",
           "markets": [{"ticker": "T-O25", "strike_type": "greater", "floor_strike": 2.5}]}
    assert s.price_event(st, st.events[9], kev) == []


def _poisoned_transient_state():
    """The ticket's one-sided-refresh shape in ONE ladder: main 2.5 healthy
    (sum 1.048), alt 1.5 crossed (sum ~0.82), alt 3.5 blown (sum ~1.35)."""
    return feed.parse_snapshot({
        "marketSources": [{"id": 7, "name": "S"}],
        "teams": {"1": {"name": "Colombia"}, "2": {"name": "Ghana"}},
        "gameOddsEvents": {"lg21:pt1:pregame": [{
            "eventId": 9, "eventStart": "2026-06-22T17:00:00+00:00",
            "eventTeams": {"1": {"id": 1}, "0": {"id": 2}},
            "gameOddsMarketSourcesLines": {
                "si0:ms7:an0": {"bt3": {"price": 115, "points": 2.5, "alternateLines": [
                    {"points": 1.5, "price": 144}, {"points": 3.5, "price": -208}]}},
                "si1:ms7:an0": {"bt3": {"price": -140, "points": 2.5, "alternateLines": [
                    {"points": 1.5, "price": 144}, {"points": 3.5, "price": -208}]}},
            }}]}}, {"lg21"})


def test_transient_one_sided_refresh_only_healthy_rung_prices():
    """Poisoned rungs are marked with distinct reasons + full provenance; the
    healthy rung in the SAME ladder still prices (the gate is per-rung, not
    per-ladder)."""
    s = Soccer()
    st = _poisoned_transient_state()
    ladder = s._anchor_ladder(st, 9)
    assert set(ladder) == {1.5, 2.5, 3.5}
    assert ladder[1.5]["reject"] == "anchor_crossed" and ladder[1.5]["p_over"] is None
    assert ladder[3.5]["reject"] == "anchor_overround" and ladder[3.5]["p_over"] is None
    assert ladder[2.5].get("reject") is None and 0 < ladder[2.5]["p_over"] < 1
    # rejected rungs keep their research provenance
    assert ladder[1.5]["book"] == 7 and ladder[1.5]["alt"] is True
    assert 0.81 < ladder[1.5]["overround"] < 0.83
    kev = {"title": "Colombia vs Ghana: Regulation Time Total Goals", "markets": [
        {"ticker": "T-O15", "strike_type": "greater", "floor_strike": 1.5},
        {"ticker": "T-O25", "strike_type": "greater", "floor_strike": 2.5},
        {"ticker": "T-O35", "strike_type": "greater", "floor_strike": 3.5}]}
    cands = s.price_event(st, st.events[9], kev)
    assert {c.market_ticker for c in cands} == {"T-O25"}      # only the healthy rung


def test_rejected_rung_logs_warning(caplog):
    """Rejects must be LOUD: a WARNING per rejected rung with reason + provenance."""
    import logging
    with caplog.at_level(logging.WARNING, logger="unabated_edge"):
        Soccer()._anchor_ladder(_poisoned_transient_state(), 9)
    msgs = [r.getMessage() for r in caplog.records]
    assert any("anchor_crossed" in m and "line=1.5" in m for m in msgs)
    assert any("anchor_overround" in m and "line=3.5" in m for m in msgs)


def test_all_rejected_book_falls_through_to_next_anchor():
    """A book whose every rung fails the gate is not 'found' — fall through to
    the next anchor source's healthy ladder (same-book-only preserved)."""
    s = Soccer()
    st = feed.parse_snapshot({
        "marketSources": [{"id": 7, "name": "S7"}, {"id": 6, "name": "S6"}],
        "teams": {"1": {"name": "A"}, "2": {"name": "B"}},
        "gameOddsEvents": {"lg21:pt1:pregame": [{
            "eventId": 9, "eventStart": "2026-07-11T21:00:00+00:00",
            "eventTeams": {"1": {"id": 1}, "0": {"id": 2}},
            "gameOddsMarketSourcesLines": {
                "si0:ms7:an0": {"bt3": {"price": 150, "points": 2.5}},   # book 7: crossed
                "si1:ms7:an0": {"bt3": {"price": 150, "points": 2.5}},
                "si0:ms6:an0": {"bt3": {"price": 100, "points": 2.5}},   # book 6: healthy
                "si1:ms6:an0": {"bt3": {"price": -120, "points": 2.5}},
            }}]}}, {"lg21"})
    ladder = s._anchor_ladder(st, 9)
    assert ladder[2.5]["book"] == 6 and ladder[2.5].get("reject") is None


def test_fail_closed_mismatched_points():
    """Over and Under at different `points` share no common line → empty ladder."""
    s = Soccer()
    st = feed.parse_snapshot({
        "marketSources": [{"id": 7, "name": "S"}],
        "teams": {"1": {"name": "Colombia"}, "2": {"name": "Ghana"}},
        "gameOddsEvents": {
            "lg21:pt1:pregame": [{
                "eventId": 9,
                "eventStart": "2026-06-22T17:00:00+00:00",
                "eventTeams": {"1": {"id": 1}, "0": {"id": 2}},
                "gameOddsMarketSourcesLines": {
                    "si0:ms7:an0": {"bt3": {"price": 115, "points": 2.5}},
                    "si1:ms7:an0": {"bt3": {"price": -140, "points": 3.5}},  # different line!
                }
            }]
        }
    }, {"lg21"})
    assert s._anchor_totals(st, 9) == {}
