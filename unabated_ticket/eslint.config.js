// ESLint flat config for the Unabated Ticket extension and its node tests.
//
// Errors only for real defects — undefined variables, unused variables,
// unreachable code and the rest of eslint's "recommended" set. No style
// rules: the extension is plain JS loaded unpacked, and formatting churn
// would bury the one diff line that matters in a pre-merge review.
//
// The pure modules (feed.js, bets.js, teams.js, kelly.js, edgemove.js,
// betsview.js, scanner.js, locate.js) are dual-loaded: as a plain <script>
// in panel.html, where they publish onto globalThis.UnabatedX, and via
// require() in tests/. That is why extension files get the CommonJS globals
// (module, require) next to the browser + chrome.* ones — one config for
// the pattern, no per-file disables.
//
// Run: npm run lint   (or ./check.sh for lint + node tests + pytest)
const js = require("@eslint/js");
const globals = require("globals");

module.exports = [
  { ignores: ["node_modules/**"] },
  js.configs.recommended,
  {
    rules: {
      "no-unused-vars": ["error", {
        // `catch (_error)` is the codebase's own way of saying "ignored on
        // purpose"; a rest-destructure (`const { x, ...rest } = o`) names x
        // only to leave it out of rest.
        caughtErrorsIgnorePattern: "^_",
        ignoreRestSiblings: true,
      }],
    },
  },
  {
    // Extension scripts: side panel, service worker, content scripts (MAIN
    // and ISOLATED world). All classic scripts, never ES modules.
    files: ["extension/**/*.js"],
    languageOptions: {
      ecmaVersion: "latest",
      sourceType: "script",
      globals: {
        ...globals.browser,
        ...globals.serviceworker,
        ...globals.webextensions,
        ...globals.commonjs,
      },
    },
  },
  {
    // node --test suites and this config file.
    files: ["tests/**/*.js", "eslint.config.js"],
    languageOptions: {
      ecmaVersion: "latest",
      sourceType: "commonjs",
      globals: { ...globals.node },
    },
  },
];
