from unabated_edge import ev


def test_positive_and_negative():
    assert ev.edge_for_yes(0.55, 0.45)[0] > 0
    assert ev.edge_for_yes(0.40, 0.50)[0] < 0
