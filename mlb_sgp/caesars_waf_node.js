// Browser-free AWS WAF token mint for Caesars (api.americanwagering.com).
//
// Runs AWS WAF's real challenge.js under node with minimal DOM shims and calls
// window.AwsWafIntegration.forceRefreshToken() -> getToken(). The SDK performs
// the NetworkBandwidth challenge itself (GET /inputs, multipart mp_verify POSTs)
// using node's native fetch / crypto.subtle / Blob / FormData / performance.
// The only browser APIs the SDK needs that node lacks are shimmed below
// (document scripts-scan, navigator/screen fingerprint, FileReader for the
// bandwidth blob). No Chromium.
//
// Usage:  node recon_caesars_waf_node.js [issuerHost] [issuerKey] [UA]
// Emits a single JSON line on stdout: {"token": "...", "ok": true}
// Diagnostics go to stderr (suppressed unless CZR_WAF_DEBUG=1).
'use strict';

const ISSUER = process.argv[2] || '4ad3fec456d9.edge.sdk.awswaf.com';
const KEY = process.argv[3] || '4ad3fec456d9';
const UA = process.argv[4] ||
  'Mozilla/5.0 (Macintosh; Intel Mac OS X 10_15_7) AppleWebKit/537.36 ' +
  '(KHTML, like Gecko) Chrome/149.0.0.0 Safari/537.36';
const ORIGIN = 'https://sportsbook.caesars.com';
const DEBUG = process.env.CZR_WAF_DEBUG === '1';
const dlog = (...a) => { if (DEBUG) process.stderr.write(a.join(' ') + '\n'); };

// ---------- browser-global shims ----------
const store = {};
const localStorage = {
  getItem: (k) => (k in store ? store[k] : null),
  setItem: (k, v) => { store[k] = String(v); },
  removeItem: (k) => { delete store[k]; },
};
function fakeWebGL() {
  const exts = ['ANGLE_instanced_arrays', 'EXT_blend_minmax', 'OES_texture_float',
    'WEBGL_debug_renderer_info', 'WEBGL_lose_context'];
  return {
    getParameter(p) {
      if (p === 37445) return 'Google Inc. (Apple)';
      if (p === 37446) return 'ANGLE (Apple, ANGLE Metal Renderer: Apple M1, Unspecified Version)';
      return 'WebGL';
    },
    getExtension(n) { return n === 'WEBGL_debug_renderer_info' ? { UNMASKED_VENDOR_WEBGL: 37445, UNMASKED_RENDERER_WEBGL: 37446 } : {}; },
    getSupportedExtensions() { return exts; },
    getContextAttributes() { return {}; },
    getShaderPrecisionFormat() { return { precision: 23, rangeMin: 127, rangeMax: 127 }; },
  };
}
function fake2d() {
  return {
    fillRect() {}, fillText() {}, getImageData() { return { data: new Uint8ClampedArray(4) }; },
    measureText() { return { width: 100 }; }, beginPath() {}, arc() {}, fill() {},
    save() {}, restore() {}, translate() {}, rotate() {}, bezierCurveTo() {}, rect() {},
    closePath() {}, stroke() {}, font: '', fillStyle: '', strokeStyle: '',
  };
}
function makeEl() {
  return {
    style: {}, setAttribute() {}, getAttribute() { return null; }, appendChild() {},
    remove() {}, addEventListener() {}, removeEventListener() {},
    getContext(t) { return (t && t.indexOf('webgl') >= 0) ? fakeWebGL() : fake2d(); },
    toDataURL() { return 'data:image/png;base64,iVBORw0KGgoAAAANSUhEUg=='; },
    width: 300, height: 150, classList: { add() {}, remove() {}, contains() { return false; } },
  };
}
const documentShim = {
  cookie: '', referrer: '', title: 'Caesars Sportsbook', visibilityState: 'visible',
  hidden: false, readyState: 'complete', documentElement: makeEl(), body: makeEl(),
  scripts: [{ innerHTML: '', src: '', textContent: '' }],
  forms: [], images: [], links: [], all: { length: 0 },
  createElement: () => makeEl(), getElementsByTagName: () => [], getElementById: () => null,
  querySelector: () => null, querySelectorAll: () => [],
  addEventListener: () => {}, removeEventListener: () => {},
  createEvent: () => ({ initEvent() {} }),
};
const _plugins = [
  { name: 'PDF Viewer', description: '', filename: 'internal-pdf-viewer' },
  { name: 'Chrome PDF Viewer', description: '', filename: 'internal-pdf-viewer' },
  { name: 'Chromium PDF Viewer', description: '', filename: 'internal-pdf-viewer' },
  { name: 'Microsoft Edge PDF Viewer', description: '', filename: 'internal-pdf-viewer' },
  { name: 'WebKit built-in PDF', description: '', filename: 'internal-pdf-viewer' },
];
_plugins.item = (i) => _plugins[i];
_plugins.namedItem = (n) => _plugins.find((p) => p.name === n) || null;
const _mimes = [{ type: 'application/pdf', suffixes: 'pdf' }, { type: 'text/pdf', suffixes: 'pdf' }];
_mimes.item = (i) => _mimes[i];
const navigatorShim = {
  userAgent: UA, appVersion: UA.replace('Mozilla/', ''), appName: 'Netscape',
  appCodeName: 'Mozilla', product: 'Gecko', productSub: '20030107', platform: 'MacIntel',
  language: 'en-US', languages: ['en-US', 'en'], hardwareConcurrency: 8, deviceMemory: 8,
  maxTouchPoints: 0, webdriver: false, vendor: 'Google Inc.', vendorSub: '',
  plugins: _plugins, mimeTypes: _mimes, doNotTrack: null, cookieEnabled: true, onLine: true,
  connection: { effectiveType: '4g', downlink: 10, rtt: 50, saveData: false },
  permissions: { query: async () => ({ state: 'granted' }) },
  sendBeacon: () => true, getBattery: async () => ({ level: 1, charging: true }), credentials: {},
  userAgentData: { brands: [{ brand: 'Chromium', version: '149' }], mobile: false, platform: 'macOS' },
};
const screenShim = { width: 1920, height: 1080, availWidth: 1920, availHeight: 1032, colorDepth: 24, pixelDepth: 24 };
const locationShim = {
  href: ORIGIN + '/us/nj/bet/baseball', origin: ORIGIN, protocol: 'https:',
  host: 'sportsbook.caesars.com', hostname: 'sportsbook.caesars.com',
  pathname: '/us/nj/bet/baseball', search: '', hash: '',
  reload: () => {}, assign: () => {}, replace: () => {},
};
const win = globalThis;
// node v24 exposes a read-only `navigator` global getter — override via
// defineProperty (plain assignment throws "has only a getter").
function setGlobal(name, value) {
  try { win[name] = value; } catch (e) {
    Object.defineProperty(win, name, { value, writable: true, configurable: true });
  }
}
win.window = win; win.self = win; win.top = win;
setGlobal('document', documentShim); setGlobal('navigator', navigatorShim);
setGlobal('screen', screenShim);
win.location = locationShim; win.localStorage = localStorage; win.sessionStorage = localStorage;
win.origin = ORIGIN; win.gokuProps = { key: '', iv: '', context: '' };
// pre-stub so challenge.js auto-init (reads window.AwsWafIntegration.checkForceRefresh)
// doesn't crash before the real object is assigned at end of the IIFE.
win.AwsWafIntegration = {
  checkForceRefresh: () => Promise.resolve(false), getToken: () => Promise.resolve(''),
  forceRefreshToken: () => Promise.resolve(), hasToken: () => false, saveReferrer: () => {},
};
win.addEventListener = () => {}; win.removeEventListener = () => {};
win.requestAnimationFrame = (cb) => setTimeout(() => cb(Date.now()), 16);
win.matchMedia = () => ({ matches: false, addListener() {}, removeListener() {} });
win.Worker = class { postMessage() {} terminate() {} addEventListener() {} };
win.WebGLRenderingContext = function () {}; win.WebGL2RenderingContext = function () {};
win.HTMLCanvasElement = function () {}; win.HTMLElement = function () {};
win.Image = class { constructor() { this.addEventListener = () => {}; } };
win.OffscreenCanvas = class { getContext(t) { return (t && t.indexOf('webgl') >= 0) ? fakeWebGL() : fake2d(); } };
if (!performance.getEntriesByType) performance.getEntriesByType = () => [];
if (!performance.getEntries) performance.getEntries = () => [];
win.performance = performance;
win.PerformanceObserver = class { observe() {} disconnect() {} };
win.devicePixelRatio = 2;
// FileReader: the NetworkBandwidth challenge reads its measurement Blob via FileReader.
win.FileReader = class {
  constructor() { this.onload = null; this.onloadend = null; this.onerror = null; this.result = null; this.readyState = 0; }
  _emit() {
    this.readyState = 2; const ev = { target: this };
    if (this.onload) this.onload(ev);
    if (this.onloadend) this.onloadend(ev);
    if (this._l) (this._l.load || []).forEach((f) => f(ev));
    if (this._l) (this._l.loadend || []).forEach((f) => f(ev));
  }
  addEventListener(t, f) { (this._l = this._l || {})[t] = (this._l[t] || []).concat(f); }
  removeEventListener() {}
  async readAsArrayBuffer(b) { this.result = await b.arrayBuffer(); this._emit(); }
  async readAsText(b) { this.result = await b.text(); this._emit(); }
  async readAsDataURL(b) { const buf = Buffer.from(await b.arrayBuffer()); this.result = 'data:' + (b.type || '') + ';base64,' + buf.toString('base64'); this._emit(); }
  async readAsBinaryString(b) { const buf = Buffer.from(await b.arrayBuffer()); this.result = buf.toString('binary'); this._emit(); }
  abort() {}
};
win.File = class extends Blob { constructor(parts, name, opts) { super(parts, opts); this.name = name; } };
win.URL = URL;
if (!URL.createObjectURL) URL.createObjectURL = () => 'blob:czr/' + Math.random().toString(36).slice(2);
if (!URL.revokeObjectURL) URL.revokeObjectURL = () => {};
win.XMLHttpRequest = undefined; // force SDK to use fetch

// inject Origin/Referer/UA on the SDK's outbound WAF requests
const _origFetch = globalThis.fetch;
globalThis.fetch = (url, opts = {}) => {
  opts.headers = Object.assign({ Origin: ORIGIN, Referer: ORIGIN + '/', 'User-Agent': UA }, opts.headers || {});
  if (DEBUG) {
    const u = typeof url === 'string' ? url : (url && url.url) || String(url);
    dlog('[fetch]', opts.method || 'GET', u);
  }
  return _origFetch(url, opts);
};

function fail(msg) { process.stdout.write(JSON.stringify({ ok: false, error: msg }) + '\n'); process.exit(1); }
process.on('uncaughtException', (e) => fail('uncaught:' + (e && e.message || e)));
process.on('unhandledRejection', (e) => fail('unhandled:' + (e && e.message || e)));

(async () => {
  const r = await _origFetch(`https://${ISSUER}/${KEY}/challenge.js`, {
    headers: { Origin: ORIGIN, Referer: ORIGIN + '/', 'User-Agent': UA },
  });
  if (!r.ok) return fail('challenge.js ' + r.status);
  const js = await r.text();
  // eslint-disable-next-line no-eval
  (0, eval)(js); // runs IIFE -> assigns the real window.AwsWafIntegration
  const integ = win.AwsWafIntegration;
  if (!integ || typeof integ.forceRefreshToken !== 'function') return fail('no integration');
  await integ.forceRefreshToken();
  const token = await integ.getToken();
  if (!token) return fail('empty token');
  process.stdout.write(JSON.stringify({ ok: true, token }) + '\n');
  process.exit(0);
})();
