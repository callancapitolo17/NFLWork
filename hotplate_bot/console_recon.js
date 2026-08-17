/*
 * Paste-into-DevTools recon. No Python, no saving files, no setup.
 *
 * Open the live drop page (during the waiting room is the useful moment), open
 * DevTools -> Console, paste this whole file, press Enter. It prints every
 * item->button pair the sniper would consider and copies a compact report to
 * your clipboard.
 *
 * Read-only: it inspects the DOM and clicks nothing.
 *
 * The item->button pairing is what matters. hotplate_snipe.py finds an item by
 * its name text, then walks UP to the nearest ancestor that also contains a
 * button — so a card whose name and button do not share an ancestor within 8
 * levels is invisible to the sniper, and this report is where you would see that.
 */
(() => {
  const BUTTON_PATTERN = /add to cart|add|order now|reserve|sold out/i;
  const MAX_CARD_DEPTH = 8;

  const readText = (element) => ((element.innerText || element.textContent || "").trim());

  const isClickable = (element) => {
    if (!element || !element.isConnected || element.disabled) return false;
    if (element.getAttribute("aria-disabled") === "true") return false;
    const rect = element.getBoundingClientRect();
    if (rect.width === 0 || rect.height === 0) return false;
    const style = window.getComputedStyle(element);
    return style.visibility !== "hidden" && style.display !== "none" && style.pointerEvents !== "none";
  };

  // Mirror of the sniper's walk, in reverse: from a button, find the item name
  // that shares its card. If this comes back "(no name found)", the sniper will
  // not be able to target that button by item name either.
  const itemNameForButton = (button) => {
    let current = button.parentElement;
    for (let depth = 0; depth < MAX_CARD_DEPTH && current; depth += 1) {
      const heading = current.querySelector("h1,h2,h3,h4,h5,h6");
      if (heading && readText(heading)) return readText(heading);
      current = current.parentElement;
    }
    // No heading: fall back to the card's own text minus the button label.
    let container = button.parentElement;
    for (let depth = 0; depth < MAX_CARD_DEPTH && container; depth += 1) {
      const text = readText(container).replace(readText(button), "").trim();
      if (text.length > 2) return text.split("\n")[0].slice(0, 60);
      container = container.parentElement;
    }
    return "(no name found)";
  };

  const buttons = [...document.querySelectorAll("button, [role=button], a")]
    .filter((element) => BUTTON_PATTERN.test(readText(element)));

  const report = {
    url: location.href,
    pageTitle: document.title,
    capturedAt: new Date().toISOString(),
    items: buttons.map((button) => ({
      itemName: itemNameForButton(button),
      buttonLabel: readText(button).slice(0, 40),
      clickableNow: isClickable(button),
    })),
    allButtonTexts: [...new Set(
      [...document.querySelectorAll("button, [role=button]")].map(readText).filter(Boolean)
    )].slice(0, 40),
  };

  console.log("%c=== HOTPLATE RECON ===", "font-weight:bold;font-size:14px");
  console.table(report.items);
  console.log("All button texts:", report.allButtonTexts);

  if (report.items.length === 0) {
    console.warn("No add-to-cart-ish buttons found. Copy `allButtonTexts` above — the pattern needs widening.");
  } else if (report.items.every((entry) => !entry.clickableNow)) {
    console.log("%cAll buttons currently locked — expected during a waiting room. Matcher is working.", "color:green");
  } else {
    console.log("%cSome buttons are already clickable — the drop is OPEN right now.", "color:orange;font-weight:bold");
  }

  const asText = JSON.stringify(report, null, 2);
  if (navigator.clipboard) {
    navigator.clipboard.writeText(asText)
      .then(() => console.log("%cReport copied to clipboard — paste it to Claude.", "color:green;font-weight:bold"))
      .catch(() => console.log("Clipboard blocked; copy the object below by hand.\n" + asText));
  } else {
    console.log(asText);
  }
  return report;
})();
