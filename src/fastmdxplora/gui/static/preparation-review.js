/* Preparation observations remain separate from requested choices. Every
 * atom highlight is mapped against the actual retained structure by identity. */
(function () {
  "use strict";
  let api, dialog, opener, canvas, study = "", generation = 0, view = 0;
  let payload = null, change = null, engine = null, reading = null, loading = null, accepting = false;
  let outcome = null;
  const disposed = new WeakSet();
  function dispose(item) { if (item && !disposed.has(item)) { disposed.add(item); item.dispose(); } }
  const el = (id) => document.getElementById(id);
  const text = (id, value) => { el(id).textContent = value || ""; };
  function node(tag, value) { const item = document.createElement(tag); item.textContent = value || ""; return item; }
  function drop() {
    view += 1;
    if (loading) loading.abort();
    loading = null;
    dispose(engine);
    engine = null;
    canvas.replaceChildren(); canvas.hidden = true; canvas.removeAttribute("data-ready");
  }
  function clear() {
    generation += 1;
    if (reading) reading.abort();
    reading = null; payload = null; change = null; outcome = null; opener.hidden = true;
    window.FastMDXDialog.close("preparation-review-dialog"); drop();
  }
  function render() {
    if (!payload) return;
    opener.hidden = !payload.available && !payload.gate?.requested;
    const gate = payload.gate || {};
    if (outcome && outcome.signature !== signature()) outcome = null;
    text("preparation-review-notice", payload.notice || payload.reason);
    text("preparation-review-status", outcome?.message || gate.reason || "Inspect recorded observations. Accepting setup is available only for a local study launched for review.");
    text("preparation-review-plan", gate.simulation ? JSON.stringify(gate.simulation, null, 2) : "No pending local simulation plan is available.");
    const accept = el("preparation-accept");
    accept.hidden = !payload.control || !gate.requested;
    accept.disabled = accepting || !!outcome?.uncertain || !gate.can_accept;
    const events = el("preparation-review-events"), choices = el("preparation-review-choices");
    const focused = document.activeElement?.dataset.preparationEvent;
    events.replaceChildren(); choices.replaceChildren();
    (payload.warnings || []).forEach((warning) => choices.append(node("p", warning)));
    (payload.events || []).forEach((event) => {
      if (event.operation === "recorded_setup_choices") {
        const detail = node("details"), summary = node("summary", "Requested choices and recorded decisions");
        detail.append(summary, node("p", event.reason), node("pre", JSON.stringify(event.details, null, 2)));
        choices.append(detail);
      } else {
        const item = node("li"), button = node("button", event.operation.replace(/_/g, " "));
        button.type = "button"; button.className = "ghost-btn"; button.dataset.preparationEvent = event.id;
        button.onclick = () => select(event);
        item.append(button); events.append(item);
      }
    });
    if (focused) events.querySelector('[data-preparation-event="' + focused + '"]')?.focus({preventScroll: true});
  }
  function signature() {
    const gate = payload?.gate || {};
    return JSON.stringify([study, payload?.version, gate.revision, gate.status, gate.reason, gate.can_accept]);
  }
  async function load(manual = false) {
    if (manual) { outcome = null; render(); }
    if (!study || reading) return;
    const at = generation, key = study, controller = reading = new AbortController();
    try {
      const response = await fetch("/api/preparation-review?study=" + encodeURIComponent(key), {cache: "no-store", signal: controller.signal});
      if (!response.ok) throw new Error("The preparation record could not be read.");
      const answer = await response.json();
      if (at !== generation || answer.study !== study) return;
      if (payload && payload.version !== answer.version) {
        drop(); change = null;
        el("preparation-change-atoms").replaceChildren();
        text("preparation-change-note", "The preparation record changed. Choose an operation from the refreshed record.");
      }
      payload = answer; render();
    } catch (error) {
      if (at === generation && error.name !== "AbortError") text("preparation-review-status", error.message);
    } finally { if (reading === controller) reading = null; }
  }
  async function select(event) {
    drop(); change = null;
    const ticket = view, at = generation, key = study, version = payload?.version;
    const controller = loading = new AbortController();
    text("preparation-change-title", event.operation.replace(/_/g, " "));
    text("preparation-change-details", [event.reason, JSON.stringify(event.details || {}, null, 2)].filter(Boolean).join("\n"));
    text("preparation-change-note", "Reading retained observations...");
    el("preparation-change-atoms").replaceChildren();
    el("preparation-before").disabled = el("preparation-after").disabled = true;
    try {
      const query = new URLSearchParams({study: key, version, event: event.id});
      const response = await fetch("/api/preparation-change?" + query, {cache: "no-store", signal: controller.signal});
      if (!response.ok) throw new Error("The preparation observation could not be read.");
      const answer = await response.json();
      if (at !== generation || ticket !== view || dialog.hidden || answer.study !== study) return;
      change = answer;
      el("preparation-before").disabled = !answer.before;
      el("preparation-after").disabled = !answer.after;
      if (answer.after || answer.before) await show(answer.after ? "after" : "before");
      else text("preparation-change-note", answer.reason || answer.after_reason || answer.before_reason || "No supported retained coordinates for this operation.");
      if (at !== generation || dialog.hidden || change !== answer || answer.study !== study) return;
      if (!answer.compared) text("preparation-change-note", (el("preparation-change-note").textContent + " "
        + (answer.before_reason || answer.after_reason || "Atom comparison is unavailable.")).trim());
      if (answer.compared) {
        const host = el("preparation-change-atoms");
        [["added", "after", "Added"], ["removed", "before", "Removed"]].forEach(([kind, side, label]) => {
          (answer[kind] || []).slice(0, 200).forEach((identity) => {
            const button = node("button", label + " " + identity.join(" / "));
            button.className = "ghost-btn"; button.type = "button";
            button.onclick = () => show(side, identity); host.append(button);
          });
          if ((answer[kind] || []).length > 200) host.append(node("p", label + ": first 200 atom controls shown; all changes remain in the audit."));
        });
        if (!(answer.added?.length || answer.removed?.length)) host.append(node("p", "No atom identities were added or removed between these snapshots."));
      }
    } catch (error) {
      if (at === generation && ticket === view && error.name !== "AbortError") text("preparation-change-note", error.message);
    } finally { if (loading === controller) loading = null; }
  }
  async function show(side, identity) {
    const selected = change, source = selected?.[side];
    if (!source) return;
    drop(); const ticket = view, at = generation;
    const stale = () => ticket !== view || at !== generation || dialog.hidden || selected !== change;
    let made = null;
    try {
      canvas.hidden = false;
      const background = getComputedStyle(document.documentElement).getPropertyValue("--background-primary").trim();
      const colour = /^#[0-9a-f]{6}$/i.test(background) ? parseInt(background.slice(1), 16) : 0x1b1f23;
      made = await window.FastMDXViewerEngine.create(canvas, {axes: true, background: colour});
      if (stale()) { dispose(made); return; }
      engine = made; made.representation = "sticks";
      made.scene = {allAtoms: true, show: {hydrogens: true}};
      await made.loadStructure({text: source.pdb, label: "Preparation " + side});
      if (stale()) { dispose(made); return; }
      canvas.dataset.ready = "true";
      text("preparation-change-note", (side === "before" ? "Before" : "After") + ": " + source.source + ". Frozen retained coordinates.");
      if (identity) {
        const matching = made.find({}).filter((atom) => JSON.stringify([String(atom.chain).trim(), String(atom.resi),
          String(atom.icode || "").trim(), atom.resn, atom.atom]) === JSON.stringify(identity));
        const evidence = source.atoms.filter((atom) => JSON.stringify(atom.identity) === JSON.stringify(identity));
        if (matching.length !== 1 || evidence.length !== 1 || ["x", "y", "z"].some((axis, i) =>
            Math.abs(matching[0][axis] - evidence[0].xyz[i]) > 0.005)) {
          text("preparation-change-note", "This atom cannot be mapped uniquely to the rendered snapshot. No highlight was guessed.");
          return;
        }
        await made.showSelected([matching[0].index]);
        if (!stale()) made.focus([matching[0].index]);
      }
    } catch (error) {
      if (!stale()) text("preparation-change-note", error.message || "The retained structure could not be rendered.");
      if (engine === made) { engine = null; dispose(made); }
    }
  }
  async function accept() {
    if (accepting || outcome?.uncertain || !payload?.gate?.can_accept) return;
    outcome = null; accepting = true; render();
    const key = study, revision = payload.gate.revision;
    try {
      const response = await fetch("/api/accept-preparation", {method: "POST", headers: {"Content-Type": "application/json"},
        body: JSON.stringify({study: key, revision})});
      const answer = await response.json();
      if (key !== study) return;
      if (answer.ok) { payload.gate.can_accept = false; payload.gate.status = "accepted"; }
      outcome = {signature: signature(), message: answer.ok ? "Setup accepted. Simulation is starting on this prepared study." : answer.error || "Simulation could not start."};
    } catch (error) {
      if (key === study) outcome = {signature: signature(), uncertain: true,
        message: "The server did not answer. Refresh the study before retrying; simulation may have started."};
    } finally { accepting = false; render(); }
  }
  document.addEventListener("DOMContentLoaded", function () {
    api = window.FastMDXDashboard; dialog = el("preparation-review-dialog");
    opener = el("preparation-review-open"); canvas = el("preparation-change-canvas");
    if (!api || !dialog || !window.FastMDXDialog) return;
    opener.onclick = () => { window.FastMDXDialog.open("preparation-review-dialog", opener); load(); };
    el("preparation-review-refresh").onclick = () => load(true);
    el("preparation-before").onclick = () => show("before");
    el("preparation-after").onclick = () => show("after");
    el("preparation-accept").onclick = accept;
    new MutationObserver(() => { if (dialog.hidden) drop(); }).observe(dialog, {attributes: true, attributeFilter: ["hidden"]});
    api.on("run-changed", clear);
    api.on("study-addressed", (data) => { if (study !== data.key) { study = data.key || ""; clear(); } load(); });
    api.on("app-state", () => { if (!dialog.hidden) load(); });
    window.FastMDXPreparationReview = {refresh: load, get engine() { return engine; },
      beforeStudySwitch: () => accepting ? "Wait for setup acceptance to finish before switching studies." : ""};
  });
})();
