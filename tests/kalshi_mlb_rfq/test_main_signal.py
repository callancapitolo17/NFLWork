import threading

from kalshi_mlb_rfq import main


def test_signal_handler_sets_running_false():
    main._running.set()
    main._signal_handler(15, None)  # SIGTERM
    assert not main._running.is_set()
    main._running.set()  # restore for other tests


def test_accept_lock_serializes():
    """Two threads attempting accept simultaneously should serialize."""
    counter = {"value": 0}
    snapshots = []

    def attempt():
        with main.ACCEPT_LOCK:
            v = counter["value"]
            snapshots.append(v)
            counter["value"] = v + 1

    threads = [threading.Thread(target=attempt) for _ in range(20)]
    for t in threads: t.start()
    for t in threads: t.join()
    assert sorted(snapshots) == list(range(20))
