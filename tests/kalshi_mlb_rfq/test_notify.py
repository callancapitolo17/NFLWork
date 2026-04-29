import json
from unittest.mock import patch, MagicMock

from kalshi_mlb_rfq import notify


def test_fill_appends_log_line(tmp_path, monkeypatch):
    log = tmp_path / "bot.log"
    monkeypatch.setattr(notify, "LOG_PATH", log)
    notify.fill(rfq_id="r1", combo_market_ticker="MT", contracts=10,
                price=0.30, ev_pct=0.06)
    assert "[FILL]" in log.read_text()
    assert "r1" in log.read_text()


def test_halt_appends_log_line(tmp_path, monkeypatch):
    log = tmp_path / "bot.log"
    monkeypatch.setattr(notify, "LOG_PATH", log)
    notify.halt(reason="kill_switch", detail="user request")
    assert "[HALT]" in log.read_text()
    assert "kill_switch" in log.read_text()


def test_webhook_post_when_url_set(tmp_path, monkeypatch):
    monkeypatch.setattr(notify, "LOG_PATH", tmp_path / "bot.log")
    monkeypatch.setattr(notify, "NOTIFY_WEBHOOK_URL", "https://hook.example.com/x")
    fake_resp = MagicMock()
    fake_resp.__enter__.return_value = fake_resp
    fake_resp.__exit__.return_value = False
    fake_resp.status = 200
    with patch("kalshi_mlb_rfq.notify.urllib.request.urlopen",
               return_value=fake_resp) as m:
        notify.fill(rfq_id="r1", combo_market_ticker="MT", contracts=1,
                    price=0.30, ev_pct=0.06)
    assert m.called
    sent_data = m.call_args[0][0].data.decode()
    assert json.loads(sent_data)["event"] == "fill"
