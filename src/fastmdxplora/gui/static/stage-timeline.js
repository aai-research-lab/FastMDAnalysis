/* Retained stage structures have their own frozen canvas. The study's live
 * Viewer and its scientific results keep their own atom identities. */
(function () {
  "use strict";
  var api, host, list, note, dialog, canvas;
  var study = "", run = "", runName = "", generation = 0, view = 0, lastRead = 0;
  var reading = null, loading = null, engine = null, selected = null, signature = "";
  var disposed = new WeakSet();
  var ids = ["input", "prepared", "solvated", "minimized", "equilibrated", "production"];
  var words = {retained: "Structure retained", completed: "Completed", current: "In progress",
    failed: "Failed", error: "Failed", stopped: "Stopped", interrupted: "Interrupted",
    skipped: "Skipped", unknown: "Not recorded"};
  function el(id) { return document.getElementById(id); }
  function drop(item) { if (item && !disposed.has(item)) { disposed.add(item); item.dispose(); } }
  function closeView() {
    var returnStage = selected?.stage.id;
    var focus = document.activeElement;
    if (dialog.hidden && returnStage && (focus === document.body || dialog.contains(focus))) {
      list.querySelector('[data-stage="' + returnStage + '"]')?.focus({preventScroll: true});
    }
    view += 1;
    if (loading) loading.abort();
    loading = null; selected = null;
    drop(engine); engine = null;
    canvas.replaceChildren(); canvas.removeAttribute("data-ready");
  }
  function clear() {
    generation += 1;
    if (reading) reading.abort();
    reading = null; lastRead = 0; signature = ""; run = ""; runName = "";
    window.FastMDXDialog.close("stage-structure-dialog"); closeView();
    list.replaceChildren(); note.textContent = "Reading retained stage records...";
    host.hidden = !study;
  }
  function visible() {
    var page = api.state.activePage;
    return !!study && page !== "run" && page !== "studies";
  }
  function render(data) {
    if (data.study !== study) return;
    host.hidden = !visible();
    var sig = JSON.stringify(data);
    if (sig === signature) return;
    signature = sig; run = data.run || ""; runName = data.run_name || "";
    var focused = document.activeElement?.dataset.stage;
    var scroll = list.scrollLeft;
    list.replaceChildren();
    (data.stages || []).forEach(function (stage) {
      if (ids.indexOf(stage.id) < 0) return;
      var item = document.createElement("li"), button = document.createElement("button");
      button.type = "button"; button.dataset.stage = stage.id;
      button.dataset.state = stage.state;
      var label = document.createElement("span"); label.textContent = stage.label;
      var state = document.createElement("small"); state.textContent = words[stage.state] || words.unknown;
      button.append(label, state);
      button.title = stage.available ? "View retained " + stage.label.toLowerCase() + " coordinates"
        : stage.reason;
      if (!stage.available) {
        button.classList.add("stage-unavailable");
        button.setAttribute("aria-description", stage.reason || "No retained coordinates");
      }
      button.onclick = function () { open(stage, button); };
      item.appendChild(button); list.appendChild(item);
    });
    list.scrollLeft = scroll;
    if (focused) list.querySelector('[data-stage="' + focused + '"]')?.focus({preventScroll: true});
    note.textContent = data.reason || (runName ? "For " + runName + ". " : "")
      + "Completion and retained coordinates are separate records.";
    if (!(data.stages || []).length && data.available) host.hidden = true;
  }
  function load(force) {
    if (!study || reading || (!force && Date.now() - lastRead < 15000)) return;
    var at = generation, key = study, controller = reading = new AbortController();
    lastRead = Date.now();
    fetch("/api/stage-timeline?study=" + encodeURIComponent(key), {
      cache: "no-store", signal: controller.signal
    }).then(function (response) {
      if (!response.ok) throw new Error("The timeline could not be read. Refresh to try again.");
      return response.json();
    }).then(function (data) { if (at === generation && key === study) render(data); })
      .catch(function (error) {
        if (at !== generation || error.name === "AbortError") return;
        note.textContent = error.message; host.hidden = !visible();
      }).finally(function () { if (reading === controller) reading = null; });
  }
  function open(stage, opener) {
    closeView();
    selected = {stage: stage, study: study, run: run};
    el("stage-structure-title").textContent = stage.label + ": retained structure";
    el("stage-structure-source").textContent = stage.available
      ? (runName ? runName + "/" : "") + stage.source
        + (stage.topology ? " + " + (runName ? runName + "/" : "") + stage.topology : "")
      : "Coordinates unavailable";
    el("stage-structure-solvent").checked = false;
    el("stage-structure-solvent").disabled = !stage.available;
    window.FastMDXDialog.open("stage-structure-dialog", opener);
    if (stage.available) show();
    else { el("stage-structure-note").textContent = stage.reason; canvas.replaceChildren(); canvas.hidden = true; }
  }
  async function show() {
    var choice = selected;
    if (!choice) return;
    var ticket = ++view, at = generation;
    if (loading) loading.abort();
    var controller = loading = new AbortController(), created = null;
    drop(engine); engine = null; canvas.replaceChildren(); canvas.removeAttribute("data-ready");
    canvas.hidden = false;
    el("stage-structure-note").textContent = "Reading the saved coordinates...";
    var stale = function () { return ticket !== view || at !== generation || dialog.hidden || choice.study !== study; };
    try {
      var query = new URLSearchParams({study: choice.study, run: choice.run, stage: choice.stage.id,
        version: choice.stage.version, solvent: el("stage-structure-solvent").checked ? "1" : "0"});
      var response = await fetch("/structure/stage.pdb?" + query, {cache: "no-store", signal: controller.signal});
      if (!response.ok) {
        var refused = await response.json();
        throw new Error(refused.reason || "The saved coordinates could not be read.");
      }
      var pdb = await response.text();
      if (stale()) return;
      var light = getComputedStyle(document.documentElement).getPropertyValue("--background-primary").trim();
      var colour = /^#[0-9a-f]{6}$/i.test(light) ? parseInt(light.slice(1), 16) : 0x1b1f23;
      created = await window.FastMDXViewerEngine.create(canvas, {background: colour, axes: true});
      if (stale()) { drop(created); return; }
      engine = created;
      created.representation = "sticks";
      created.scene = {allAtoms: true, show: {hydrogens: true}};
      var result = await created.loadStructure({text: pdb, label: choice.stage.label});
      if (stale()) { drop(created); return; }
      if (!result.atoms) throw new Error("No atoms could be rendered from these coordinates.");
      canvas.dataset.ready = "true";
      el("stage-structure-note").textContent = result.atoms + " atoms. Frozen saved coordinates; live updates stay in the study view.";
    } catch (error) {
      if (!stale() && error.name !== "AbortError") {
        el("stage-structure-note").textContent = error.message || "The saved structure could not be displayed.";
        canvas.hidden = true;
      }
      if (created) drop(created);
    } finally { if (loading === controller) loading = null; }
  }
  function boot() {
    api = window.FastMDXDashboard; host = el("stage-timeline"); list = el("stage-timeline-list");
    note = el("stage-timeline-note"); dialog = el("stage-structure-dialog"); canvas = el("stage-structure-canvas");
    if (!api || !host || !window.FastMDXDialog) return;
    el("stage-timeline-refresh").onclick = function () { load(true); };
    el("stage-structure-solvent").onchange = show;
    new MutationObserver(function () { if (dialog.hidden) closeView(); }).observe(dialog, {attributes: true, attributeFilter: ["hidden"]});
    api.on("run-changed", function () { clear(); });
    api.on("study-addressed", function (data) {
      if (study !== data.key) { study = data.key || ""; clear(); }
      host.hidden = !visible(); load(false);
    });
    api.on("navigate", function () { host.hidden = !visible() || !list.children.length; });
    window.FastMDXStageTimeline = {refresh: function () { load(true); }, get engine() { return engine; }};
  }
  document.addEventListener("DOMContentLoaded", boot);
})();
