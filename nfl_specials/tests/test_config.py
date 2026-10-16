import pytest

from nfl_specials import config


def test_extra_allowed_hosts_are_trimmed_and_deduplicated():
    raw = " mlb-stack.tail87ab2d.ts.net, ,mlb-stack.tail87ab2d.ts.net,other.tail87ab2d.ts.net "
    assert config.parse_extra_allowed_hosts(raw) == ("mlb-stack.tail87ab2d.ts.net", "other.tail87ab2d.ts.net")


def test_no_extra_allowed_hosts_is_empty():
    assert config.parse_extra_allowed_hosts(None) == ()
    assert config.parse_extra_allowed_hosts("") == ()


def test_dk_proxy_unset_means_direct():
    assert config.parse_dk_proxy(None) is None
    assert config.parse_dk_proxy("  ") is None


def test_dk_proxy_is_an_http_host_and_port():
    assert config.parse_dk_proxy(" http://100.120.136.59:8097 ") == "http://100.120.136.59:8097"


@pytest.mark.parametrize("bad", ["100.120.136.59:8097", "https://100.120.136.59:8097", "http://100.120.136.59",
                                 "http://100.120.136.59:8097/"])
def test_a_misread_dk_proxy_refuses_to_start_and_names_the_value(bad):
    with pytest.raises(ValueError) as raised:
        config.parse_dk_proxy(bad)
    assert str(raised.value) == f"NFL_SPECIALS_DK_PROXY: expected http://<host>:<port>, got {bad!r}"
