from unabated_edge import mapping
from unabated_edge.sports.soccer import Soccer


def test_validate_rejects_incomplete():
    a = Soccer()
    p = mapping.Pairing(None, {"home": "H", "draw": "D"})       # missing away
    assert mapping.validate(a, p, {"home": .4, "draw": .3, "away": .3}) is False


def test_validate_accepts_complete():
    a = Soccer()
    p = mapping.Pairing(None, {"home": "H", "draw": "D", "away": "A"})
    assert mapping.validate(a, p, {"home": .4, "draw": .3, "away": .3}) is True
