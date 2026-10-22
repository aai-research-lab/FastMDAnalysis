/* New study and the workspace's defaults file (IMDE 2.9). A file that is
 * refused (a setting of one study in it, a value the schema refuses) fills
 * nothing, and the form offers FastMDXplora's own values meanwhile. Only
 * Settings' This folder said so: New study showed FastMDXplora's values
 * with nothing said, so a lab's file refused for one setting dropped all
 * its values silently, and the study was refused only once it started.
 * The refusal `/api/schema` gives (`defaults_refused`) is said above the
 * form each time New study is opened, with the way to This folder. */
(function () {
  "use strict";

  var asking = null;

  function el(id) { return document.getElementById(id); }

  function make(tag, className, text) {
    var node = document.createElement(tag);
    if (className) node.className = className;
    if (text !== undefined) node.textContent = text;
    return node;
  }

  /* Words with `names` in them, the names set as code. */
  function words(into, text) {
    String(text).split("`").forEach(function (part, i) {
      if (part) into.appendChild(i % 2 ? make("code", "mono", part) : document.createTextNode(part));
    });
  }

  function say(note, refused) {
    note.textContent = "";
    if (!refused) { note.hidden = true; return; }
    note.appendChild(make("b", "", "Your defaults are not used. "));
    words(note, refused);
    note.appendChild(document.createTextNode(
      " This form offers FastMDXplora's own values until the file is put right. "));
    var link = make("a", "", "This folder in Settings");
    link.href = "#settings";
    link.addEventListener("click", function (e) {
      if (!window.FastMDXSettings) return;
      e.preventDefault();
      window.FastMDXSettings.scope("folder");
    });
    note.appendChild(link);
    note.hidden = false;
  }

  function read() {
    var note = el("run-defaults-refused");
    if (!note || asking) return;
    asking = fetch("/api/schema").then(function (response) { return response.json(); })
      .then(function (schema) { say(note, schema && schema.defaults_refused); })
      .catch(function () { /* the form says it could not read the settings */ })
      .then(function () { asking = null; });
  }

  function start() {
    if (window.FastMDXDashboard) {
      window.FastMDXDashboard.on("navigate", function (detail) {
        if (detail && detail.page === "run") read();
      });
    }
    if (document.documentElement.dataset.page === "run") read();
  }

  if (document.readyState === "loading") document.addEventListener("DOMContentLoaded", start);
  else start();
}());
