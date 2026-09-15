"""NovigAuth (refresh, rotation, connect parsing) and NovigSource (trader
resolution, pagination, failures) against fakes — no network."""
import base64
import json

import pytest

from unabated_ticket.bets_service.sources import novig, novig_auth
from unabated_ticket.bets_service.sources.novig import NovigSource, source_if_connected
from unabated_ticket.bets_service.sources.novig_auth import NovigAuth, NovigAuthError, code_from_redirect, jwt_subject
from unabated_ticket.bets_service.tests.conftest import NOVIG_READ_AT


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


class FakeGraphql:
    """Serves the fixture rows in pages, records every call."""

    def __init__(self, orders, parlays, users=None, page_limit=None):
        self.orders, self.parlays = orders, parlays
        self.users = users if users is not None else [{"id": "u1", "trader_id": "t1"}]
        self.calls: list[tuple[str, dict, str]] = []
        self.page_limit = page_limit

    def __call__(self, query, variables, bearer):
        self.calls.append((query, variables, bearer))
        if "BetsService_User" in query:
            return {"user": self.users}
        rows = self.orders if "BetsService_Orders" in query else self.parlays
        limit = self.page_limit or variables["limit"]
        key = "order" if "BetsService_Orders" in query else "parlay"
        return {key: rows[variables["offset"]:variables["offset"] + limit]}


def make_source(token_file, graphql):
    auth = NovigAuth(token_file, request=FakeTokenEndpoint(rotate=False), clock=lambda: 0.0)
    return NovigSource(auth=auth, graphql=graphql, poll_sec=1, retention_days=30)


def test_fetch_resolves_the_trader_once_and_returns_normalised_records(token_file, novig_rows):
    orders, parlays = novig_rows
    graphql = FakeGraphql(orders, parlays)
    source = make_source(token_file, graphql)
    records = source.fetch()
    assert len(records) == 19
    assert {r["venue"] for r in records} == {"novig"}
    assert any(r["id"] == "novig:o-ml-car" and r["side"] == "home" for r in records)
    user_calls = [c for c in graphql.calls if "BetsService_User" in c[0]]
    assert len(user_calls) == 1 and user_calls[0][1] == {"auth_id": "auth0|user1"}
    assert all(c[2] == fake_jwt("auth0|user1") for c in graphql.calls)
    source.fetch()
    assert len([c for c in graphql.calls if "BetsService_User" in c[0]]) == 1  # cached trader id


def test_fetch_paginates_to_a_short_page(token_file, novig_rows):
    orders, parlays = novig_rows
    graphql = FakeGraphql(orders * 30, parlays)  # 450 orders: pages of 100 + a short one
    records = make_source(token_file, graphql).fetch()
    order_calls = [c for c in graphql.calls if "BetsService_Orders" in c[0]]
    assert [c[1]["offset"] for c in order_calls] == [0, 100, 200, 300, 400]
    assert len(records) == 19  # duplicates collapse on native id


def test_fetch_fails_loudly_on_no_user_auth_failure_or_graphql_error(token_file, novig_rows, tmp_path):
    orders, parlays = novig_rows
    with pytest.raises(RuntimeError, match="exactly one Novig user"):
        make_source(token_file, FakeGraphql(orders, parlays, users=[])).fetch()

    def broken(query, variables, bearer):
        raise RuntimeError("novig graphql error: field 'trader' not found")
    with pytest.raises(RuntimeError, match="novig graphql error"):
        make_source(token_file, broken).fetch()
    dead_auth = NovigAuth(token_file, request=FakeTokenEndpoint(fail=True), clock=lambda: 0.0)
    with pytest.raises(RuntimeError, match="novig auth"):
        NovigSource(auth=dead_auth, graphql=FakeGraphql(orders, parlays), poll_sec=1, retention_days=30).fetch()


def test_source_if_connected_is_none_without_a_token_file(tmp_path):
    assert source_if_connected(tmp_path / "none.json") is None
