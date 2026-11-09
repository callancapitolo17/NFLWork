// Run: node --test unabated_ticket/tests
// underdog.js: the Underdog read run in Cal's own app.underdogsports.com tab on
// the headers the app just sent, and the push body it hands the bets service
// (2026-11-09). The slip grammar is pinned by the Python side
// (bets_service/tests/test_underdog.py); these pin the headers, the trim, the
// settled walk and the reload rule.
const test = require("node:test");
const assert = require("node:assert/strict");
const underdog = require("../extension/underdog.js");
const betsView = require("../extension/betsview.js");

const NOW = Date.parse("2026-11-09T19:30:00Z");
const DAY = 24 * 60 * 60 * 1000;

function jwtWithExp(expSeconds) {
  const part = (object) => Buffer.from(JSON.stringify(object)).toString("base64url");
  return `${part({ alg: "RS256" })}.${part({ exp: expSeconds, sub: "auth0|cal" })}.sig`;
}

function slip(id, product, createdAtMs, optionId = "o-1", matchId = 1) {
  return { id, product, created_at: new Date(createdAtMs).toISOString(),
    selection_groups: [{ match_id: matchId, match_type: "Game", selections: [{ option_id: optionId }] }] };
}

function answer(slips, next = null) {
  return {
    data: {
      entry_slips: slips,
      over_under_options: [{ id: "o-1", over_under_line_id: "l-1" }, { id: "o-pick", over_under_line_id: "l-pick" }],
      over_under_lines: [{ id: "l-1", over_under_id: "ou-1" }, { id: "l-pick", over_under_id: "ou-pick" }],
      over_unders: [{ id: "ou-1", title: "LV @ SF Moneyline" }, { id: "ou-pick", title: "Player points" }],
      games: [{ id: 1, away_team_id: "t-lv", home_team_id: "t-sf" }, { id: 2, away_team_id: "t-a", home_team_id: "t-b" }],
      teams: [{ id: "t-lv", name: "Las Vegas Raiders" }, { id: "t-sf", name: "San Francisco 49ers" }, { id: "t-a" }, { id: "t-b" }],
      players: [{ id: "p-1" }],
      appearances: [{ id: "a-1" }],
    },
    meta: { next, page: 1 },
  };
}

test("readHeadersOf keeps the read headers off the app's request and drops location; none without a token", () => {
  const sent = { Accept: "application/json", "Client-Type": "web", "Client-Version": "20261106182101",
    "Client-Device-Id": "dev-1", "Client-Request-Id": "req-1", "User-Latitude": "37.7", "User-Longitude": "-122.4",
    "User-Location-Token": "loc", authorization: "eyJ.x.y" };
  assert.deepEqual(underdog.readHeadersOf(sent), { Accept: "application/json", "Client-Type": "web",
    "Client-Version": "20261106182101", "Client-Device-Id": "dev-1", Authorization: "eyJ.x.y" });
  assert.equal(underdog.readHeadersOf({ Accept: "application/json" }), null);
  assert.equal(underdog.readHeadersOf(null), null);
});

test("tokenExpiryMs reads exp, with or without Bearer; isFresh wants a minute left", () => {
  const decode = (b64) => Buffer.from(b64, "base64").toString();
  const exp = NOW / 1000 + 600;
  assert.equal(underdog.tokenExpiryMs(jwtWithExp(exp), decode), exp * 1000);
  assert.equal(underdog.tokenExpiryMs(`Bearer ${jwtWithExp(exp)}`, decode), exp * 1000);
  assert.equal(underdog.tokenExpiryMs("not-a-jwt", decode), null);
  const headers = { Authorization: "x" };
  assert.equal(underdog.isFresh({ headers, expiresAtMs: NOW + 61 * 1000 }, NOW), true);
  assert.equal(underdog.isFresh({ headers, expiresAtMs: NOW + 59 * 1000 }, NOW), false);
  assert.equal(underdog.isFresh({ headers, expiresAtMs: null }, NOW), false);
  assert.equal(underdog.isFresh(null, NOW), false);
});

test("trimPage keeps prediction slips and only the lookups they name", () => {
  const page = answer([slip("p", "prediction_markets", NOW), slip("f", "fantasy", NOW, "o-pick", 2)], 2);
  const trimmed = underdog.trimPage(page);
  assert.deepEqual(trimmed.data.entry_slips.map((s) => s.id), ["p"]);
  assert.deepEqual(trimmed.data.over_under_options.map((o) => o.id), ["o-1"]);
  assert.deepEqual(trimmed.data.over_under_lines.map((l) => l.id), ["l-1"]);
  assert.deepEqual(trimmed.data.over_unders.map((o) => o.id), ["ou-1"]);
  assert.deepEqual(trimmed.data.games.map((g) => g.id), [1]);
  assert.deepEqual(trimmed.data.teams.map((t) => t.name), ["Las Vegas Raiders", "San Francisco 49ers"]);
  assert.equal(trimmed.data.players, undefined);
  assert.deepEqual(trimmed.meta, { next: 2, page: 1 });
});

function fakeFetch(answers) {
  const calls = [];
  const fetchImpl = async (url, init) => {
    calls.push({ url, init });
    const next = answers.shift();
    if (!next) throw new Error(`unexpected request ${url}`);
    return { status: next.status || 200, json: async () => next.body };
  };
  return { fetchImpl, calls };
}

test("readAccount reads the open slips then settled pages until one reaches past 31 days", async () => {
  const { fetchImpl, calls } = fakeFetch([
    { body: answer([slip("open", "prediction_markets", NOW - DAY)]) },
    { body: answer([slip("s1", "prediction_markets", NOW - 2 * DAY)], 2) },
    { body: answer([slip("s2", "fantasy", NOW - 40 * DAY)], 3) },
  ]);
  const headers = { Authorization: "eyJ", "Client-Type": "web" };
  const { push } = await underdog.readAccount({ fetchImpl, headers, nowMs: NOW });
  assert.deepEqual(calls.map((call) => call.url), [underdog.ACTIVE_URL, `${underdog.SETTLED_URL}?page=1`, `${underdog.SETTLED_URL}?page=2`]);
  assert.equal(calls[0].init.headers, headers);
  assert.equal(calls[0].init.credentials, "omit");
  assert.equal(push.fetchedAt, new Date(NOW).toISOString());
  assert.equal(push.complete, false);
  assert.deepEqual(push.settled.map((page) => page.data.entry_slips.map((s) => s.id)), [["s1"], []]);
  assert.deepEqual(push.active.data.entry_slips.map((s) => s.id), ["open"]);
});

test("readAccount: the list ending is complete; the page cap stops a list that never ends", async () => {
  let { fetchImpl } = fakeFetch([{ body: answer([]) }, { body: answer([slip("s1", "fantasy", NOW)]) }]);
  assert.equal((await underdog.readAccount({ fetchImpl, headers: {}, nowMs: NOW })).push.complete, true);
  const endless = [{ body: answer([]) }];
  for (let n = 0; n < underdog.MAX_SETTLED_PAGES; n += 1) endless.push({ body: answer([slip(`s${n}`, "fantasy", NOW)], n + 2) });
  ({ fetchImpl } = fakeFetch(endless));
  const { push } = await underdog.readAccount({ fetchImpl, headers: {}, nowMs: NOW });
  assert.equal(push.settled.length, underdog.MAX_SETTLED_PAGES);
  assert.equal(push.complete, false);
});

test("readAccount: a refused token is stale, an old app version names the reload, a bad answer is an error", async () => {
  let { fetchImpl } = fakeFetch([{ status: 401, body: {} }]);
  assert.deepEqual(await underdog.readAccount({ fetchImpl, headers: {}, nowMs: NOW }), { error: underdog.STALE, stale: true });
  ({ fetchImpl } = fakeFetch([{ body: answer([]) }, { status: 426, body: {} }]));
  const old = await underdog.readAccount({ fetchImpl, headers: {}, nowMs: NOW });
  assert.match(old.error, /settled page 1 wants a newer app version \(HTTP 426\)/);
  assert.equal(old.stale, false);
  ({ fetchImpl } = fakeFetch([{ body: { data: {} } }]));
  assert.equal((await underdog.readAccount({ fetchImpl, headers: {}, nowMs: NOW })).error,
    "Underdog open slips answered without data.entry_slips");
});

test("mayReload: never the tab Cal is looking at, at most every five minutes", () => {
  const idle = { active: false };
  assert.equal(underdog.mayReload(idle, true, 0, NOW), true);
  assert.equal(underdog.mayReload({ active: true }, true, 0, NOW), false);
  assert.equal(underdog.mayReload({ active: true }, false, 0, NOW), true);
  assert.equal(underdog.mayReload(idle, true, NOW - underdog.RELOAD_MIN_MS + 1, NOW), false);
  assert.equal(underdog.mayReload(idle, true, NOW - underdog.RELOAD_MIN_MS, NOW), true);
});

test("serviceUrlOf follows the panel's own settings rule", () => {
  for (const stored of [undefined, null, {}, { serviceUrl: "ftp://x" }, { serviceUrl: "https://mlb-stack.tail87ab2d.ts.net/" },
    { serviceUrl: "http://127.0.0.1:8094" }]) {
    assert.equal(underdog.serviceUrlOf(stored), betsView.sanitizeBetsSettings(stored).serviceUrl);
  }
});
