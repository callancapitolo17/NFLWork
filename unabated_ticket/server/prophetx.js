// Unabated Ticket server — the Ticket's Place button at ProphetX (Cal
// 2026-10-31, option A): which board line Place buys and the ProphetX strike
// ids Unabated's feed carries for its game. The runner answers GET
// /prophetx_quote.json with these; the bets service (bets_service/
// prophetx_place.py) reads ProphetX on its key pair, finds the strike and the
// live price, and sends the order. ProphetX's Trading API needs that key for
// every read, so nothing here reads ProphetX.
//
// Unabated's ProphetX alt rungs carry ProphetX's strike id as `sourceData`
// (32 hex characters, README "Venue ids on the rungs"); main lines carry none.
// Any one of them names the ProphetX game (GET /v4/mm/get_strikes), so the
// quote sends every one the game has, the line's own first.
//
// Pure: no network, no clock. Reused by runner.js only.

// Unabated's book id for Prophet Exchange (fillfair.js BOOK_ID_OF_VENUE).
const PROPHETX_BOOK_ID = 66;
const STRIKE_ID_RE = /^[0-9a-f]{32}$/i;
// The bets service reads at most this many strikes in one get_strikes call.
const MAX_STRIKE_IDS = 20;

function isProphetxStrikeId(sourceData) {
  return typeof sourceData === "string" && STRIKE_ID_RE.test(sourceData);
}

// {lineStrikeId, strikeIds}: the line's own strike id (null on a main line)
// and the game's strike ids, the line's own first, at most MAX_STRIKE_IDS.
//   line        the board line Place buys
//   gameLines   every board line of the game
function strikeIdsOf(line, gameLines) {
  const lineStrikeId = isProphetxStrikeId(line.sourceData) ? line.sourceData : null;
  const ids = new Set(lineStrikeId ? [lineStrikeId] : []);
  for (const other of gameLines) {
    if (ids.size >= MAX_STRIKE_IDS) break;
    if (other.bookId === PROPHETX_BOOK_ID && isProphetxStrikeId(other.sourceData)) ids.add(other.sourceData);
  }
  return { lineStrikeId, strikeIds: [...ids] };
}

module.exports = { PROPHETX_BOOK_ID, MAX_STRIKE_IDS, isProphetxStrikeId, strikeIdsOf };
