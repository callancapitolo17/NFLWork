"""NovigAuth (refresh, rotation, connect parsing) and NovigSource (cursor
pagination of the Portfolio feed, failures) against fakes — no network."""
import base64
import json
from urllib.parse import parse_qs, urlparse

import pytest

from unabated_ticket.bets_service.sources import novig, novig_auth
from unabated_ticket.bets_service.sources.novig import NovigSource, portfolio_url, source_if_connected
from unabated_ticket.bets_service.sources.novig_auth import NovigAuth, NovigAuthError, code_from_redirect, jwt_subject


def fake_jwt(sub: str) -> str:
    payload = base64.urlsafe_b64encode(json.dumps({"sub": sub}).encode()).rstrip(b"=").decode()
    return f"header.{payload}.sig"


class FakeTokenEndpoint:
    def __init__(self, rotate: bool = True, fail: bool = False):
        self.calls: list[dict] = []
        self.rotate = rotate
        self.fail = fail

    def __call__(self, body: dict) -> dict:
        self.calls.append(body)
        if self.fail:
            raise NovigAuthError("Auth0 refresh_token failed (403): invalid_grant")
        payload = {"access_token": fake_jwt("auth0|user1"), "expires_in": 1800}
        if self.rotate:
            payload["refresh_token"] = f"rt-{len(self.calls)}"
        return payload


@pytest.fixture
def token_file(tmp_path):
    path = tmp_path / "novig_token.json"
    novig_auth.save_token_file(path, "rt-0", "auth0|user1")
    return path


def test_jwt_subject_reads_sub_without_verifying():
    assert jwt_subject(fake_jwt("auth0|abc")) == "auth0|abc"
    with pytest.raises(NovigAuthError):
        jwt_subject("not-a-jwt")


def test_auth_refreshes_once_then_caches_until_the_margin(token_file):
    endpoint = FakeTokenEndpoint(rotate=False)
    now = [1_000.0]
    auth = NovigAuth(token_file, request=endpoint, clock=lambda: now[0])
    first = auth.token()
    assert first.auth_id == "auth0|user1"
    assert auth.token().token == first.token and len(endpoint.calls) == 1
    now[0] += 1800 - 60  # inside REFRESH_MARGIN_SEC of expiry
    auth.token()
    assert len(endpoint.calls) == 2
    assert endpoint.calls[0]["grant_type"] == "refresh_token" and endpoint.calls[0]["refresh_token"] == "rt-0"


def test_a_rotated_refresh_token_is_persisted_at_once(token_file):
    auth = NovigAuth(token_file, request=FakeTokenEndpoint(rotate=True), clock=lambda: 0.0)
    auth.token()
    assert json.loads(token_file.read_text())["refresh_token"] == "rt-1"
    assert oct(token_file.stat().st_mode & 0o777) == "0o600"


def test_missing_or_broken_token_file(tmp_path):
    with pytest.raises(NovigAuthError, match="connect"):
        NovigAuth(tmp_path / "none.json", request=FakeTokenEndpoint()).token()
    broken = tmp_path / "broken.json"
    broken.write_text("{}")
    with pytest.raises(NovigAuthError, match="no refresh_token"):
        NovigAuth(broken, request=FakeTokenEndpoint()).token()


def test_code_from_redirect_checks_state_and_errors():
    assert code_from_redirect("https://app.novig.us/?code=abc&state=s1", "s1") == "abc"
    with pytest.raises(NovigAuthError, match="state mismatch"):
        code_from_redirect("https://app.novig.us/?code=abc&state=old", "s1")
    with pytest.raises(NovigAuthError, match="no \\?code="):
        code_from_redirect("https://app.novig.us/?state=s1", "s1")
    with pytest.raises(NovigAuthError, match="login refused"):
        code_from_redirect("https://app.novig.us/?error=access_denied&error_description=nope&state=s1", "s1")


class FakePortfolio:
    """Serves the fixture's two lists in cursor pages, records every GET."""

    def __init__(self, active, settled, page=3, items_key="items"):
        self.lists = {"active": active, "settled": settled}
        self.page = page
        self.items_key = items_key
        self.calls: list[tuple[str, str]] = []

    def __call__(self, url, bearer):
        self.calls.append((url, bearer))
        parsed = urlparse(url)
        tab = parsed.path.rsplit("/", 1)[1]
        query = parse_qs(parsed.query)
        offset = int(query.get("cursor", ["0"])[0])
        rows = self.lists[tab]
        page = rows[offset:offset + self.page]
        next_offset = offset + self.page
        return {self.items_key: page, "nextCursor": str(next_offset) if next_offset < len(rows) else None}


def make_source(token_file, http_get):
    auth = NovigAuth(token_file, request=FakeTokenEndpoint(rotate=False), clock=lambda: 0.0)
    return NovigSource(auth=auth, http_get=http_get, poll_sec=1, base_url="https://api.example/nbx/v1/")


def test_fetch_pages_both_lists_to_the_end_and_returns_normalised_records(token_file, novig_fixture):
    portfolio = FakePortfolio(novig_fixture["active"], novig_fixture["settled"])
    records = make_source(token_file, portfolio).fetch()
    assert len(records) == 31
    assert {r["venue"] for r in records} == {"novig"}
    assert any(r["id"] == "novig:01a0cabf-8ef6-7cc2-8bde-8451ccd8ca33" and r["side"] == "home" for r in records)
    urls = [url for url, _bearer in portfolio.calls]
    assert urls[0] == "https://api.example/nbx/v1/portfolio/active?currency=CASH&sort=active_recency&limit=50"
    assert urls[1] == "https://api.example/nbx/v1/portfolio/active?currency=CASH&sort=active_recency&limit=50&cursor=3"
    assert [u.rsplit("/", 1)[1].split("?")[0] for u in urls] == ["active", "active", "settled", "settled", "settled", "settled"]  # 4 and 11 cards, pages of 3
    assert all(bearer == fake_jwt("auth0|user1") for _url, bearer in portfolio.calls)


def test_portfolio_url():
    assert portfolio_url("https://api.novig.us/nbx/v1", "settled", "settled_recency", 50, None) == \
        "https://api.novig.us/nbx/v1/portfolio/settled?currency=CASH&sort=settled_recency&limit=50"
    assert portfolio_url("https://api.novig.us/nbx/v1", "active", "active_recency", 10, "abc").endswith("&limit=10&cursor=abc")


def test_fetch_fails_loudly_on_a_bad_page_an_unended_list_auth_failure_or_http_error(token_file, novig_fixture, monkeypatch):
    with pytest.raises(RuntimeError, match="expected an `items` list"):
        make_source(token_file, FakePortfolio(novig_fixture["active"], novig_fixture["settled"], items_key="cards")).fetch()

    def never_ends(url, bearer):
        return {"items": [], "nextCursor": "again"}
    monkeypatch.setattr(novig, "MAX_PAGES", 3)
    with pytest.raises(RuntimeError, match="did not end within 3 pages"):
        make_source(token_file, never_ends).fetch()

    def broken(url, bearer):
        raise RuntimeError("novig portfolio HTTP 503: upstream")
    with pytest.raises(RuntimeError, match="HTTP 503"):
        make_source(token_file, broken).fetch()
    dead_auth = NovigAuth(token_file, request=FakeTokenEndpoint(fail=True), clock=lambda: 0.0)
    with pytest.raises(RuntimeError, match="novig auth"):
        NovigSource(auth=dead_auth, http_get=FakePortfolio([], []), poll_sec=1).fetch()


def test_source_if_connected_is_none_without_a_token_file(tmp_path):
    assert source_if_connected(tmp_path / "none.json") is None
