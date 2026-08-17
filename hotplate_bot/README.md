# hotplate_bot

Races the add-to-cart click at the end of a Hotplate waiting-room countdown.

Unrelated to the betting stack in the rest of this repo — it lives here only because
this is the repo the session had access to. It reads and writes nothing outside
`hotplate_bot/`, and touches no DuckDB.

## What it actually does, and what it can't

Hotplate drops open with a ~5 minute waiting-room countdown; at zero, the add-to-cart
buttons enable for everyone at once. **Adding to cart is what reserves the item**, so
that single click is the entire contest. Checkout afterwards runs against a cart-hold
timer that is generous by comparison.

So this tool optimises exactly one thing: the latency between "button becomes
clickable" and "button is clicked."

| | Latency to click |
|---|---|
| Human (reaction + mouse travel) | ~250–400ms |
| This bot (measured, `test_sniper.py`) | **1–2ms** |

That gap is the whole edge. It's achieved by running the detect-and-click loop
*inside the page* rather than from Python — a Python-side polling loop pays a
CDP round-trip plus interpreter jitter on every check, which costs more than the
entire margin being competed for.

**It cannot help you if:**
- The drop sells out to people ahead of you in a server-side queue. Speed only wins
  the client-side race; it can't reorder a queue.
- Hotplate randomises queue position after the countdown (some drops do this — if
  yours does, this tool is pointless and you should stop here).
- You're not signed in with a saved card before the drop. That's most of the
  realised losses, and no amount of speed fixes it.

## Setup

```bash
cd hotplate_bot
pip install -r requirements.txt
playwright install chromium

cp config.example.json config.json   # then edit it
```

Edit `config.json`:
- `event_url` — the drop link you were sent.
- `target_items` — **priority order**. The first item whose button becomes clickable
  is the one clicked. Partial matches, case-insensitive.

## Usage

```bash
# 1. One time: sign in by hand, session is saved to storage_state.json.
#    Add your card and address in this window too.
python hotplate_snipe.py login

# 2. Before the drop: confirm the bot can see your items and their buttons.
python hotplate_snipe.py recon

# 3. At the drop, several minutes early: arm and leave it alone.
python hotplate_snipe.py arm
```

### Read the recon output

`recon` is the step people skip and then lose the drop. It prints every button on
the page and whether your target items were found. During a waiting room you want
to see your items listed and their buttons `DISABLED` — that's proof the matcher is
locked onto the right elements and is simply waiting. If it prints
`Items found: NONE`, your `target_items` don't match the page text and `arm` will
sit there and time out.

If no button matches, widen `add_to_cart_pattern` using the button texts recon dumps.

### Console recon (no setup at all)

The fastest check, and the one to use during a live waiting room. Open the drop page,
DevTools → Console, paste all of `console_recon.js`, Enter:

```
=== HOTPLATE RECON ===
┌─────────┬────────────────────────┬───────────────┬───────────────┐
│ (index) │        itemName        │  buttonLabel  │ clickableNow  │
├─────────┼────────────────────────┼───────────────┼───────────────┤
│    0    │ 'Pistachio Croissant'  │ 'Add to Cart' │     false     │
│    1    │    'Cinnamon Roll'     │ 'Add to Cart' │     false     │
└─────────┴────────────────────────┴───────────────┴───────────────┘
All buttons currently locked — expected during a waiting room. Matcher is working.
Report copied to clipboard.
```

It reads the DOM and clicks nothing. Use the exact `itemName` strings it prints as your
`target_items` — those are the strings the sniper matches against, so copying them
removes the guesswork entirely.

It also pairs each button back to its item name, which is the check that matters: the
sniper finds an item by name then walks *up* to the nearest ancestor containing a
button. A card whose name and button don't share an ancestor within 8 levels is
invisible to it, and `(no name found)` in this table is where you'd catch that.

### Offline recon

The probe only needs the DOM, so it can run against a saved copy of the page with no
network at all:

```bash
python hotplate_snipe.py recon --html-file ~/Downloads/drop.html
```

Save the page from the browser (`Cmd+S`, or DevTools → Elements → right-click `<html>`
→ Copy → Copy outerHTML, into a file). Useful when the machine holding the page can't
reach the machine checking the config — including handing the HTML to someone else to
sanity-check your `target_items` against.

### What `arm` does

Loads the page, warns you if your clock is skewed from Hotplate's server, then waits.
On open it clicks add-to-cart for the highest-priority available item, screenshots,
plays a sound, and **hands the browser to you** to finish quantity, pickup slot, and
payment. It does not submit payment — an unattended mis-click on a real card is not
worth the seconds it would save.

Leave the process running after the click; closing it closes the browser.

## Tests

```bash
python test_sniper.py
```

Runs the sniper against a mock drop page — buttons disabled through a waiting room,
then enabled. Covers the ordinary enable-in-place case, a full framework re-render at
open (which detaches cached elements), falling through to a second choice when the
top item is sold out, not clicking while still disabled, and item targeting. Asserts
click latency measured *inside the page*, so a regression back to Python-side polling
fails the suite.

## Deliberately not implemented

Multi-account entry stacking, proxy rotation, CAPTCHA solving, browser-fingerprint
spoofing. This drives one ordinary logged-in browser on your own machine, faster than
your hand can — it doesn't pretend to be several people or hide that it's automated.

## Worth knowing before you run it

Automating orders almost certainly violates Hotplate's terms of service, and a ban
costs you every future drop from this baker. Against a single order's upside that's a
bad trade on a repeated game — which is a reason to keep this to one account and not
hammer the site, not a reason it won't work.

`storage_state.json` holds live session cookies. It's gitignored and written mode
600; treat it like a password.
