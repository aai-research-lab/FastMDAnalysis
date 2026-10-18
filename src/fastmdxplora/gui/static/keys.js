/* The GUI's keys, in one place, and the sheet that lists every one of them
 * (?, anywhere but in a field).
 *
 * Keys that act on the whole GUI are here: the study's pages, the sidebar,
 * the side panel, a file by name. The palette (palette.js, ⌘K), the
 * Agent (agent-beside.js, ⌘J) and the Viewer's own (molecule-viewer.js,
 * viewer-rail.js) keep their handlers and are listed here, so the sheet and
 * the palette say each key once (FastMDXKeys.of(name) for a command's key).
 *
 * ⌘ on a Mac, Ctrl elsewhere. A browser keeps ⌘1 to ⌘9 (Ctrl elsewhere) for
 * its own tabs, so the pages are also on ⌥1 to ⌥5 (Alt elsewhere), which a
 * browser leaves to the page; an installed window has both.
 */
(function () {
  "use strict";

  var MAC = /Mac|iPhone|iPad/.test(navigator.platform || navigator.userAgent || "");
  var PAGES = [["overview", "Overview"], ["viewer", "Viewer"], ["analysis", "Analysis"],
               ["report", "Report"], ["files", "Files"]];

  function dashboard() { return window.FastMDXDashboard; }
  function go(page) { return function () { if (dashboard()) dashboard().navigate(page); }; }
  function clickId(id) { var b = document.getElementById(id); if (b && !b.hidden) b.click(); }

  function toggleSidebar() {
    clickId(document.body.classList.contains("sidebar-collapsed") ? "sidebar-expand" : "sidebar-collapse");
  }
  function togglePanel() {
    if (!window.FastMDXFrame) return;
    window.FastMDXFrame.setCollapsed(!document.body.classList.contains("panel-collapsed"));
  }
  function openSheet() { if (window.FastMDXDialog) window.FastMDXDialog.open("viewer-help"); }

  /* Each key: what it is (`mod` is ⌘ or Ctrl), what it does, and the
   * command's name in the palette where it has one. `here` is false for a
   * key another script answers. */
  var KEYS = [
    { keys: ["mod", "K"], does: "Go to anything: a page, a study, a file, an analysis, a setting",
      name: "palette", here: false },
    { keys: ["mod", "P"], does: "Go to a file of the open study", name: "file",
      run: function () { if (window.FastMDXPalette) window.FastMDXPalette.open("files"); } },
    { keys: ["mod", "J"], does: "The Agent beside the page", name: "agent-beside", here: false },
  ].concat(PAGES.map(function (page, i) {
    return { keys: ["mod", String(i + 1)], also: ["alt", String(i + 1)], does: "The study's " + page[1],
             name: "page:" + page[0], run: go(page[0]) };
  })).concat([
    { keys: ["mod", "B"], does: "Fold or show the sidebar", name: "sidebar", run: toggleSidebar },
    { keys: ["mod", "alt", "B"], does: "Fold or show the side panel", name: "panel", run: togglePanel },
    { keys: ["?"], does: "Every key (this sheet)", name: "keys", run: openSheet, bare: true },
  ]);

  var WORD = { mod: MAC ? "⌘" : "Ctrl", alt: MAC ? "⌥" : "Alt", shift: MAC ? "⇧" : "Shift" };

  function said(keys) { return keys.map(function (k) { return WORD[k] || k; }); }

  /* In a browser's tab its own keys come first (⌘1 is its first tab); an
   * installed window leaves them all to the page. */
  var WINDOW = window.matchMedia && (window.matchMedia("(display-mode: standalone)").matches
    || window.matchMedia("(display-mode: window-controls-overlay)").matches);

  /* The key of a command, as words for <kbd>s ([] where it has none): the
   * one that works where the page is. */
  function of(name) {
    var key = KEYS.find(function (k) { return k.name === name; });
    if (!key) return [];
    return said(key.also && !WINDOW ? key.also : key.keys);
  }

  function matches(keys, event) {
    var mod = keys.indexOf("mod") >= 0, alt = keys.indexOf("alt") >= 0, shift = keys.indexOf("shift") >= 0;
    var main = keys[keys.length - 1];
    if (mod !== (MAC ? event.metaKey : event.ctrlKey)) return false;
    if (!MAC && event.metaKey) return false;
    if (MAC && event.ctrlKey) return false;
    if (alt !== event.altKey || shift !== event.shiftKey) return false;
    if (/^\d$/.test(main)) return event.code === "Digit" + main;
    if (/^[A-Z]$/.test(main)) return event.code === "Key" + main;
    return false;
  }

  function typing(target) {
    return target && (target.isContentEditable || /^(INPUT|TEXTAREA|SELECT)$/.test(target.tagName));
  }

  document.addEventListener("keydown", function (event) {
    if (event.isComposing || event.defaultPrevented) return;
    // A dialog over the page has the keyboard.
    if (document.querySelector(".agent-dialog:not([hidden])")) return;
    for (var i = 0; i < KEYS.length; i += 1) {
      var key = KEYS[i];
      if (key.here === false || !key.run) continue;
      if (key.bare) {
        // The Viewer answers its own ? (viewer-rail.js): the same sheet.
        if (event.key !== "?" || event.ctrlKey || event.metaKey || event.altKey) continue;
        if (typing(event.target) || document.documentElement.dataset.page === "viewer") continue;
      } else if (!matches(key.keys, event) && !(key.also && matches(key.also, event))) {
        continue;
      }
      if (!key.bare && event.target && event.target.isContentEditable) continue;
      event.preventDefault();
      if (event.repeat) return;
      key.run();
      return;
    }
  });

  /* The sheet: every key above, then the Viewer's. */
  function fillSheet() {
    var host = document.getElementById("every-key");
    if (!host) return;
    host.innerHTML = "";
    KEYS.forEach(function (key) {
      var row = document.createElement("div");
      var dt = document.createElement("dt");
      said(key.keys).forEach(function (word) {
        var kbd = document.createElement("kbd");
        kbd.textContent = word;
        dt.appendChild(kbd);
        dt.appendChild(document.createTextNode(" "));
      });
      var dd = document.createElement("dd");
      dd.textContent = key.does + (key.also ? " (also " + said(key.also).join(" ") + ")" : "");
      row.appendChild(dt);
      row.appendChild(dd);
      host.appendChild(row);
    });
  }

  if (document.readyState === "loading") document.addEventListener("DOMContentLoaded", fillSheet);
  else fillSheet();

  window.FastMDXKeys = { of: of, list: function () { return KEYS.slice(); }, mac: MAC };
}());
