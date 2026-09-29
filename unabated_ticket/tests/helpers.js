// Shared test helpers (not a test file: package.json runs tests/*.test.js).

// Simulate a book moving one line between snapshots: overwrite `fields` on
// the held line in place, the way a fresh snapshot would. Throws when the key
// is not in the state, so a typo cannot pass as "nothing moved".
function moveLine(state, key, fields) {
  const held = state.lines[key];
  if (!held) throw new Error(`moveLine: expected line ${key} in the state, found none`);
  state.lines[key] = { ...held, ...fields };
  return state.lines[key];
}

module.exports = { moveLine };
