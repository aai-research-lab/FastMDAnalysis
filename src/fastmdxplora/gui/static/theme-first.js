/* The scheme before the first paint, as frame.js chooses it: Light or
 * Dark where one was chosen (Graphite and Ink, chosen before, are Dark;
 * Paper is Light), otherwise the computer's. Loaded in the page's head, so
 * a light computer's page does not open dark and turn light.
 *
 * The accent too: the person's own (`fmx.accent`), else the one the
 * service hosting the GUI chose (`data-accent-default` on the root
 * element, `fastmdx gui --accent`), else Cyan, so the page never shows
 * Cyan first and turns. */
(function () {
  var kept = null, former = { graphite: "dark", ink: "dark", paper: "light" };
  try { kept = localStorage.getItem("fmx.theme"); } catch (e) { /* private mode */ }
  kept = former[kept] || kept;
  if (kept !== "light" && kept !== "dark") {
    try { kept = window.matchMedia("(prefers-color-scheme: light)").matches ? "light" : "dark"; }
    catch (e) { kept = "dark"; }
  }
  var root = document.documentElement;
  root.dataset.theme = kept;
  var accents = ["cyan", "blue", "violet", "magenta", "graphite"];
  var accent = null;
  try { accent = localStorage.getItem("fmx.accent"); } catch (e) { /* private mode */ }
  if (accents.indexOf(accent) < 0) accent = root.getAttribute("data-accent-default");
  if (accents.indexOf(accent) >= 0 && accent !== "cyan") root.dataset.accent = accent;
  else root.removeAttribute("data-accent");
})();

/* The Display size, before the first paint (IMDE 2.9; the ideas session's
 * large-screen entry, Auto by default, user 10-21): everything at once,
 * text included, as CSS `zoom` on the page. Auto picks from the screen the
 * browser reports, in CSS pixels, so a system's own scaling is already in
 * it and never doubled: the GUI was drawn for about 1920 by 1080, and Auto
 * is the largest of 100%, 110%, 125% and 150% that the screen holds as many
 * times over (100% to about 1920 by 1080, 125% at 2560 by 1440 or an
 * ultrawide 3440 by 1440, 150% at 3840 by 2160). It is chosen as the page
 * opens and again only when the window is moved to another screen, never
 * while it is resized. A size chosen in Settings (`fmx.displaySize`, 90 to
 * 150) is used as it is.
 *
 * A page zoomed this way reports what is on the screen in zoomed pixels
 * (an element's box, where the pointer is, the window's width) while the
 * page lays itself out and is positioned in its own: every popover, menu
 * and tooltip placed from a box would land a quarter off at 125%, and the
 * Viewer's canvas would draw at the page's pixels, not the screen's. While
 * the page is zoomed, what it is told of the screen (a box, the pointer,
 * the window's width and height however read) is said in its own pixels,
 * as at 100%, and the screen's pixels per page pixel include the zoom; at
 * 100% nothing is changed. Where the browser does not zoom so (no
 * `currentCSSZoom`), the page stays at 100% and Settings says to use the
 * browser's own zoom. */
(function () {
  var root = document.documentElement;
  var STEPS = [1, 1.1, 1.25, 1.5];
  var CHOICES = ["auto", "90", "100", "110", "125", "150"];
  var works = typeof Element !== "undefined" && "currentCSSZoom" in Element.prototype;
  var zoom = 1;
  var shimmed = false;
  var screenWas = "";

  function chosen() {
    var saved = null;
    try { saved = localStorage.getItem("fmx.displaySize"); } catch (e) { /* private mode */ }
    return CHOICES.indexOf(saved) >= 0 ? saved : "auto";
  }

  function auto() {
    var w = (window.screen && screen.width) || 1920;
    var h = (window.screen && screen.height) || 1080;
    var room = Math.min(w / 1920, h / 1080);
    var best = 1;
    STEPS.forEach(function (s) { if (s <= room + 0.001) best = s; });
    return best;
  }

  function wanted() {
    var c = chosen();
    return c === "auto" ? auto() : Number(c) / 100;
  }

  function apply(z) {
    zoom = works ? z : 1;
    if (zoom === 1) {
      root.style.removeProperty("zoom");
      root.removeAttribute("data-display-size");
    } else {
      root.style.zoom = String(zoom);
      root.setAttribute("data-display-size", String(Math.round(zoom * 100)));
      shim();
    }
    root.style.setProperty("--fmx-zoom", String(zoom));
    try { window.dispatchEvent(new CustomEvent("fmx:display-size", { detail: { zoom: zoom } })); }
    catch (e) { /* an old browser: the size still applies */ }
  }

  /* What the page is told of the screen, in its own pixels while zoomed. */
  function shim() {
    if (shimmed) return;
    shimmed = true;
    function scaled(r) {
      return zoom === 1 ? r : new DOMRect(r.x / zoom, r.y / zoom, r.width / zoom, r.height / zoom);
    }
    function boxes(proto) {
      var rect = proto.getBoundingClientRect, rects = proto.getClientRects;
      if (rect) proto.getBoundingClientRect = function () { return scaled(rect.call(this)); };
      if (rects) proto.getClientRects = function () {
        var all = rects.call(this);
        if (zoom === 1) return all;
        return Array.prototype.map.call(all, scaled);
      };
    }
    boxes(Element.prototype);
    if (window.Range) boxes(Range.prototype);
    function divided(proto, names) {
      names.forEach(function (name) {
        var d = Object.getOwnPropertyDescriptor(proto, name);
        if (!d || !d.get || !d.configurable) return;
        Object.defineProperty(proto, name, { configurable: true, enumerable: d.enumerable,
          get: function () { return d.get.call(this) / zoom; } });
      });
    }
    divided(MouseEvent.prototype, ["clientX", "clientY", "x", "y", "offsetX", "offsetY"]);
    if (window.Touch) divided(Touch.prototype, ["clientX", "clientY"]);
    ["pageX", "pageY"].forEach(function (name) {
      var d = Object.getOwnPropertyDescriptor(MouseEvent.prototype, name);
      var client = name === "pageX" ? "clientX" : "clientY";
      var scroll = name === "pageX" ? "scrollX" : "scrollY";
      if (!d || !d.configurable) return;
      Object.defineProperty(MouseEvent.prototype, name, { configurable: true, enumerable: d.enumerable,
        get: function () { return zoom === 1 ? d.get.call(this) : this[client] + window[scroll]; } });
    });
    ["innerWidth", "innerHeight"].forEach(function (name) {
      var d = Object.getOwnPropertyDescriptor(window, name);
      if (!d || !d.configurable || !d.get) return;
      Object.defineProperty(window, name, { configurable: true, enumerable: d.enumerable,
        get: function () { return d.get.call(window) / zoom; } });
    });
    // The window's size read from the root element (as tooltips.js reads
    // it): the root's own box is the page's, and only the root's is the
    // window's.
    ["clientWidth", "clientHeight", "scrollWidth", "scrollHeight"].forEach(function (name) {
      var d = Object.getOwnPropertyDescriptor(Element.prototype, name);
      if (!d || !d.configurable || !d.get) return;
      Object.defineProperty(Element.prototype, name, { configurable: true, enumerable: d.enumerable,
        get: function () {
          var value = d.get.call(this);
          return this === document.documentElement ? value / zoom : value;
        } });
    });
    var ratio = Object.getOwnPropertyDescriptor(window, "devicePixelRatio");
    if (ratio && ratio.configurable && ratio.get) {
      Object.defineProperty(window, "devicePixelRatio", { configurable: true, enumerable: ratio.enumerable,
        get: function () { return ratio.get.call(window) * zoom; } });
    }
    ["elementFromPoint", "elementsFromPoint", "caretPositionFromPoint", "caretRangeFromPoint"].forEach(function (name) {
      var at = Document.prototype[name];
      if (typeof at !== "function") return;
      Document.prototype[name] = function (x, y) { return at.call(this, x * zoom, y * zoom); };
    });
  }

  /* Rules sized by the window (vh, vw) sized by the page's own window:
   * zoomed, 100vh is the window's height in screen pixels and would run
   * a quarter past its foot at 125%. Copied, after the page's own rules,
   * divided by the zoom, for while it is zoomed. A property sized so
   * anywhere is copied from every rule that sets it, in the page's order,
   * so a later rule's own value still wins over an earlier one sized by
   * the window: copied alone, `.app-shell`'s first `min-height: 100vh` came
   * after its later `min-height: 0`, and the shell ran under the status
   * bar. */
  var VIEWPORT = /(-?\d*\.?\d+)(dvh|svh|lvh|vh|dvw|svw|lvw|vw|vmin|vmax)\b/g;
  function sizedByTheWindow() {
    if (document.getElementById("fmx-display-size-rules")) return;
    function scope(selector) {
      return selector.split(",").map(function (one) {
        var s = one.trim();
        if (/^:root\b/.test(s)) return s.replace(/^:root/, ":root[data-display-size]");
        if (/^html\b/.test(s)) return s.replace(/^html/, "html[data-display-size]");
        return ":where(html[data-display-size]) " + s;
      }).join(", ");
    }
    function each(rules, one) {
      Array.prototype.forEach.call(rules, function (rule) {
        if (rule.style && rule.selectorText) {
          for (var i = 0; i < rule.style.length; i += 1) one(rule, rule.style[i]);
        } else if (rule.cssRules) each(rule.cssRules, one);
      });
    }
    var sheets = [];
    Array.prototype.forEach.call(document.styleSheets, function (sheet) {
      try { sheets.push(sheet.cssRules); } catch (e) { /* another origin's sheet: none here */ }
    });
    var sized = {};
    sheets.forEach(function (rules) {
      each(rules, function (rule, name) {
        VIEWPORT.lastIndex = 0;
        if (VIEWPORT.test(rule.style.getPropertyValue(name))) sized[name] = true;
      });
    });
    function walk(rules, into) {
      Array.prototype.forEach.call(rules, function (rule) {
        if (rule.cssRules && rule.media) {
          var inner = [];
          walk(rule.cssRules, inner);
          if (inner.length) into.push("@media " + rule.media.mediaText + " {\n" + inner.join("\n") + "\n}");
          return;
        }
        if (rule.cssRules && rule.conditionText !== undefined && rule.constructor.name === "CSSSupportsRule") {
          var held = [];
          walk(rule.cssRules, held);
          if (held.length) into.push("@supports " + rule.conditionText + " {\n" + held.join("\n") + "\n}");
          return;
        }
        if (!rule.style || !rule.selectorText) return;
        var props = [];
        for (var i = 0; i < rule.style.length; i += 1) {
          var name = rule.style[i];
          var value = rule.style.getPropertyValue(name);
          if (!sized[name] || !value) continue;
          var priority = rule.style.getPropertyPriority(name);
          props.push(name + ": " + value.replace(VIEWPORT, "calc($1$2 / var(--fmx-zoom, 1))") +
                     (priority ? " !" + priority : "") + ";");
        }
        if (props.length) into.push(scope(rule.selectorText) + " { " + props.join(" ") + " }");
      });
    }
    var out = [];
    sheets.forEach(function (rules) { walk(rules, out); });
    var style = document.createElement("style");
    style.id = "fmx-display-size-rules";
    style.textContent = out.join("\n");
    document.head.appendChild(style);
  }

  apply(wanted());
  try { screenWas = screen.width + "x" + screen.height; } catch (e) { /* none */ }
  if (document.readyState === "loading") document.addEventListener("DOMContentLoaded", sizedByTheWindow);
  else sizedByTheWindow();
  // Moved to another screen: Auto chooses again. Resized on one: nothing.
  window.addEventListener("resize", function () {
    var now = "";
    try { now = screen.width + "x" + screen.height; } catch (e) { return; }
    if (now === screenWas) return;
    screenWas = now;
    if (chosen() === "auto") apply(auto());
  });

  window.FastMDXDisplay = {
    works: works,
    choices: CHOICES,
    chosen: chosen,
    auto: auto,
    zoom: function () { return zoom; },
    /* "auto", or "90" to "150": saved, and applied at once. */
    choose: function (value) {
      if (CHOICES.indexOf(value) < 0) return;
      try { localStorage.setItem("fmx.displaySize", value); } catch (e) { /* this page only */ }
      apply(wanted());
    },
  };
})();
