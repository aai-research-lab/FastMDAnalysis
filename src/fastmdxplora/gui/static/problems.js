/* Recorded problems in the study on screen. This panel only reads; the
 * existing Overview keeps the confirmed actions for running a remedy. */
(function () {
  "use strict";

  var active = null;
  var study = "";
  var generation = 0;
  var pending = null;
  var lastRead = 0;
  var rendered = "";
  var RECHECK_MS = 15000;
  // Read for the status bar's count while the tab is not shown, less often.
  var COUNTED_MS = 60000;

  function byId(id) { return document.getElementById(id); }
  function element(tag, cls, text) {
    var node = document.createElement(tag);
    if (cls) node.className = cls;
    if (text != null) node.textContent = text;
    return node;
  }
  function visible() { return byId("problems-pane") && !byId("problems-pane").hidden; }
  /* The count, for the status bar (statusbar.js); null where not known. */
  function counted(count) {
    window.dispatchEvent(new CustomEvent("dashboard:problems-counted", { detail: { count: count } }));
  }
  /* Where the tab stands, in its summary, and below it the page's one
   * block (page-states.js): reading, nothing open, or not read with Try
   * again. */
  function message(words, kind) {
    counted(null);
    byId("problems-summary").textContent = words;
    var body = byId("problems-body");
    body.replaceChildren();
    if (kind && window.FastMDXPageState) {
      body.appendChild(window.FastMDXPageState.block(kind === "error" ? {
        kind: "error", title: "The problems could not be read", detail: words,
        action: { label: "Try again", run: function () { load(true); } },
      } : kind === "loading" ? { kind: "loading", title: words } : { kind: "empty", title: words }));
    }
    rendered = "";
  }
  function clear() {
    generation += 1;
    if (pending) pending.abort();
    pending = null;
    byId("problems-refresh").disabled = false;
    study = "";
    lastRead = 0;
    message(active ? "Reading the study's records…" : "Open a study to read its problems.",
            active ? "loading" : "empty");
  }
  function render(data) {
    counted(data.available === false ? null : (data.items || []).length);
    var key = JSON.stringify(data);
    if (key === rendered) return;
    rendered = key;
    var body = byId("problems-body");
    body.replaceChildren();
    var items = data.items || [];
    byId("problems-summary").textContent = items.length
      ? items.length + (items.length === 1 ? " recorded problem" : " recorded problems")
      : (data.available ? "No problems listed in the available records." : "Records unavailable.");
    items.forEach(function (item) {
      var card = element("article", "problem-item");
      card.dataset.problemId = item.id;
      card.appendChild(element("h3", "", item.title));
      card.appendChild(element("p", "", item.reason));
      if (item.code) card.appendChild(element("code", "problem-code", item.code));
      card.appendChild(element("h4", "", "What would fix it"));
      card.appendChild(element("p", "", item.fix || "Review the recorded reason before choosing a change."));
      if (item.file) {
        var file = element("div", "problem-file");
        file.appendChild(element("code", "", item.file));
        var open = window.FastMDXIcons.button("eye", "Open " + item.file + " in Files");
        open.type = "button";
        open.addEventListener("click", function () {
          window.FastMDXFrame.previewPath(item.file);
        });
        file.appendChild(open);
        card.appendChild(file);
      }
      body.appendChild(card);
    });
    (data.notes || []).forEach(function (note) {
      body.appendChild(element("p", "problems-note", note));
    });
  }
  function load(force) {
    if (!active || !study || pending) return;
    if (!force && Date.now() - lastRead < (visible() ? RECHECK_MS : COUNTED_MS)) return;
    var at = generation;
    var wanted = study;
    var root = active;
    var controller = new AbortController();
    pending = controller;
    lastRead = Date.now();
    byId("problems-refresh").disabled = true;
    fetch("/api/problems?study=" + encodeURIComponent(wanted), {
      cache: "no-store", signal: controller.signal
    }).then(function (response) {
      if (!response.ok) throw new Error("Could not read the study's problems.");
      return response.json();
    }).then(function (data) {
      if (at !== generation || root !== active) return;
      if (data.ok && data.active_run === null) {
        message("Open a study to read its problems.", "empty");
        return;
      }
      if (data.study !== wanted || data.active_run !== root) {
        message("The open study changed. Waiting for its records…", "loading");
        return;
      }
      if (!data.ok) throw new Error(data.reason || "Could not read the study's problems.");
      render(data);
    }).catch(function (error) {
      if (at !== generation || error.name === "AbortError") return;
      message("Could not read problems.", "error");
    }).finally(function () {
      if (pending !== controller) return;
      pending = null;
      byId("problems-refresh").disabled = false;
    });
  }
  function boot() {
    if (!byId("problems-pane")) return;
    var panel = byId("side-panel");
    var home = panel.parentNode;
    var next = panel.nextSibling;
    var dialog = byId("problems-mobile");
    var phone = byId("problems-phone");
    phone.addEventListener("click", function () {
      dialog.appendChild(panel);
      dialog.showModal();
      window.FastMDXFrame.showTab("problems");
      byId("problems-tab").focus();
      load(true);
    });
    dialog.addEventListener("close", function () {
      home.insertBefore(panel, next);
      phone.focus();
    });
    window.matchMedia("(max-width: 800px)").addEventListener("change", function (event) {
      if (!event.matches && dialog.open) dialog.close();
    });
    byId("problems-tab").addEventListener("click", function () { load(true); });
    byId("problems-refresh").addEventListener("click", function () { load(true); });
    window.FastMDXDashboard.on("app-state", function (state) {
      var next = state.active_run || null;
      if (next === active) return;
      active = next;
      clear();
    });
    window.FastMDXDashboard.on("status-updated", function (update) {
      if (!active) return;
      if (study !== update.study) {
        clear();
        study = update.study || "";
      }
      load(false);
    });
  }
  document.addEventListener("DOMContentLoaded", boot);
}());
