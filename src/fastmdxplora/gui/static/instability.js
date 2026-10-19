/* Frozen force observations belong to one study, run and snapshot. */
(function () {
  "use strict";
  let api, dialog, canvas, study = "", serial = 0, request = null, engine = null, evidence = null, proposal = null;
  const disposed = new WeakSet();
  const el = id => document.getElementById(id);
  const text = (id, value) => { el(id).textContent = value; };
  function drop(item) { if (item && !disposed.has(item)) { disposed.add(item); item.dispose(); } }
  function clear() {
    serial += 1;
    if (request) request.abort();
    request = null; drop(engine); engine = null; evidence = null; proposal = null;
    canvas.replaceChildren(); canvas.removeAttribute("data-ready"); canvas.hidden = true;
    el("instability-ranking").replaceChildren(); el("instability-proposal").hidden = true;
    el("instability-propose").disabled = true; el("instability-download").disabled = true;
    text("instability-atom-note", ""); text("instability-source", "");
  }
  async function load(run) {
    clear(); const ticket = serial, key = study;
    const controller = request = new AbortController();
    const stale = () => ticket !== serial || key !== study || dialog.hidden;
    let made = null;
    text("instability-status", "Reading retained force evidence...");
    try {
      const query = new URLSearchParams({study: key}); if (run) query.set("run", run);
      const response = await fetch("/api/instability?" + query, {cache: "no-store", signal: controller.signal});
      const data = await response.json();
      if (stale()) return;
      if (!response.ok || data.study !== key) throw new Error(data.reason || "The open study changed.");
      const select = el("instability-run"); select.replaceChildren();
      (data.runs || []).forEach(item => { const option = document.createElement("option"); option.value = item.key; option.textContent = item.label; select.append(option); });
      select.value = data.run || ""; select.disabled = select.options.length < 2;
      if (!data.available) { text("instability-status", data.reason || "No retained finite force evidence is available."); return; }
      evidence = data;
      text("instability-status", data.notice + " " + data.failure_step_note);
      const cadence = data.stage === "minimization" ? "Observed at the start of minimization; intermediate iterations were not sampled. "
        : "Observations every " + data.interval_steps + " steps. ";
      text("instability-source", data.run_name + ": " + data.stage + ", step " + data.step + ", " + data.time_ps + " ps. " + cadence
        + "Failure reported at step " + data.failure_reported_step + " (" + data.reported_step_gap + " steps after this observation on that counter). "
        + data.force_scope + "\n" + data.source + "\nSHA-256: " + data.snapshot_sha256);
      el("instability-propose").disabled = !data.can_propose;
      canvas.hidden = false;
      const background = getComputedStyle(document.documentElement).getPropertyValue("--background-primary").trim();
      made = await window.FastMDXViewerEngine.create(canvas, {axes: true,
        background: /^#[0-9a-f]{6}$/i.test(background) ? parseInt(background.slice(1), 16) : 0x1b1f23});
      if (stale()) { drop(made); return; }
      engine = made; made.representation = "sticks"; made.scene = {allAtoms: true, show: {hydrogens: true}};
      const result = await made.loadStructure({text: data.pdb, label: "Last recorded finite observation"});
      if (stale()) { drop(made); return; }
      if (result.atoms !== data.atom_count) throw new Error("The rendered atom count does not match the full-system snapshot. No highlight was guessed.");
      canvas.dataset.ready = "true";
      const atoms = made.find({});
      let pickNumber = 0, picks = Promise.resolve();
      data.ranking.forEach(row => {
        const li = document.createElement("li"), button = document.createElement("button");
        button.type = "button"; button.className = "ghost-btn";
        button.textContent = (row.identity[0] || "(no chain)") + ":" + row.identity[1] + row.identity[2] + " " + row.identity[3] + " " + row.identity[4]
          + " — " + row.force_kj_mol_nm.toPrecision(5) + " kJ mol⁻¹ nm⁻¹";
        button.onclick = () => {
          const picked = ++pickNumber;
          picks = picks.catch(() => {}).then(async () => {
          if (stale() || engine !== made || picked !== pickNumber) return;
          const matching = atoms.filter(atom => JSON.stringify([String(atom.chain).trim(), String(atom.resi), String(atom.icode || "").trim(), atom.resn, atom.atom]) === JSON.stringify(row.identity));
          if (matching.length !== 1 || ["x", "y", "z"].some((axis, i) => Math.abs(matching[0][axis] - row.position_nm[i] * 10) > .005)) {
            text("instability-atom-note", "This atom cannot be mapped uniquely to the frozen coordinates. No highlight was guessed."); return;
          }
          await made.showSelected([matching[0].index], {hydrogens: true, emphasize: true});
          if (stale() || picked !== pickNumber) return;
          made.focus([matching[0].index]);
          text("instability-atom-note", button.textContent + ". Highlighted in this frozen full-system snapshot.");
          }).catch(error => { if (!stale() && picked === pickNumber) text("instability-atom-note", error.message || "The atom could not be highlighted."); });
        };
        li.append(button); el("instability-ranking").append(li);
      });
    } catch (error) {
      if (!stale() && error.name !== "AbortError") {
        text("instability-status", error.message || "The saved evidence could not be rendered."); canvas.hidden = true;
      }
      if (made) drop(made);
    } finally { if (request === controller) request = null; }
  }
  async function propose(download) {
    const data = evidence;
    if (!data || !data.can_propose || request) return;
    const ticket = serial, old = proposal, controller = request = new AbortController();
    el("instability-propose").disabled = true; el("instability-download").disabled = true;
    const stale = () => ticket !== serial || data !== evidence || dialog.hidden;
    try {
      const query = new URLSearchParams({study, run: data.run, version: data.version});
      const response = await fetch("/api/instability-proposal?" + query, {cache: "no-store", signal: controller.signal});
      const answer = await response.json(); if (stale()) return;
      if (!response.ok || !answer.ok || answer.study !== study || answer.run !== data.run) throw new Error(answer.reason || "A proposal could not be checked.");
      proposal = answer; el("instability-proposal").hidden = false;
      text("instability-proposal-note", answer.reason);
      text("instability-diff", answer.changes.map(change => change.setting + ": " + change.before + " → " + change.after).join("\n"));
      text("instability-yaml", answer.yaml);
      if (download && old && old.yaml === answer.yaml && old.config_sha256 === answer.config_sha256) {
        const url = URL.createObjectURL(new Blob([answer.yaml], {type: "text/yaml"}));
        const link = document.createElement("a"); link.href = url; link.download = "fastmdxplora_instability_proposal.yml"; link.click();
        setTimeout(() => URL.revokeObjectURL(url), 1000);
        text("instability-proposal-note", "Reviewed Config downloaded. The failed run is preserved. No run was started.");
      } else if (download) text("instability-proposal-note", "The source Config changed. Review the updated proposal before downloading it.");
      el("instability-download").disabled = false;
    } catch (error) {
      if (!stale() && error.name !== "AbortError") { text("instability-proposal-note", error.message); el("instability-proposal").hidden = false; }
    }
    finally { if (!stale()) el("instability-propose").disabled = false; if (request === controller) request = null; }
  }
  function open(opener) {
    // A phone's Problems pane is a native dialog. Close it before opening
    // the shared workspace dialog and return focus to its visible opener.
    const phone = el("problems-mobile");
    if (phone?.open) { phone.close(); opener = el("problems-phone"); }
    window.FastMDXDialog.open("instability-dialog", opener); load();
  }
  function boot() {
    api = window.FastMDXDashboard; dialog = el("instability-dialog"); canvas = el("instability-canvas");
    if (!api || !dialog) return;
    el("instability-open").onclick = () => open(el("instability-open"));
    el("instability-run").onchange = () => load(el("instability-run").value);
    el("instability-refresh").onclick = () => load(el("instability-run").value);
    el("instability-propose").onclick = () => propose(false);
    el("instability-download").onclick = () => propose(true);
    new MutationObserver(() => { if (dialog.hidden) clear(); }).observe(dialog, {attributes: true, attributeFilter: ["hidden"]});
    api.on("study-addressed", data => { if (study !== (data.key || "")) { study = data.key || ""; clear(); window.FastMDXDialog.close("instability-dialog"); } });
    api.on("run-changed", () => { clear(); window.FastMDXDialog.close("instability-dialog"); });
    window.FastMDXInstability = {open, get engine() { return engine; }};
  }
  document.addEventListener("DOMContentLoaded", boot);
})();
