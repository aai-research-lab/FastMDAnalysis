/* FastMDXplora dashboard-first simulation builder.
 *
 * This module only collects and validates options.  The Python backend
 * launches the canonical CLI in a separate process, so the scientific
 * workflow remains identical to a normal `fastmdx explore` command.
 */
(function () {
  "use strict";

  const byId = (id) => document.getElementById(id);
  const $$ = (selector, root) => Array.from((root || document).querySelectorAll(selector));

  const builder = {
    defaults: null,
    appState: {},
    initialized: false,
    lastValidation: null,
  };

  document.addEventListener("DOMContentLoaded", init);

  async function init() {
    const form = byId("simulation-builder-form");
    if (!form) return;
    wireForm();
    try {
      builder.defaults = await requestJSON("/api/launcher/defaults");
      populateDefaults(builder.defaults);
      builder.initialized = true;
      updateSummary();
      await refreshAppState();
    } catch (error) {
      setMessage(`Could not load simulation defaults: ${error.message}`, "error");
    }

    window.FastMDXDashboard?.on("app-state", (payload) => applyAppState(payload || {}));
    window.setInterval(refreshAppState, 2000);
  }

  function wireForm() {
    const pairs = [
      ["builder-ph-range", "builder-ph"],
      ["builder-ion-range", "builder-ion"],
      ["builder-padding-range", "builder-padding"],
      ["builder-nvt-range", "builder-nvt"],
      ["builder-npt-range", "builder-npt"],
      ["builder-production-range", "builder-production"],
      ["builder-temperature-range", "builder-temperature"],
      ["builder-timestep-range", "builder-timestep"],
      ["builder-friction-range", "builder-friction"],
    ];
    pairs.forEach(([rangeId, numberId]) => syncPair(rangeId, numberId));

    byId("simulation-builder-form")?.addEventListener("input", () => {
      clearFieldErrors();
      updateSummary();
    });
    byId("simulation-builder-form")?.addEventListener("change", updateSummary);
    byId("simulation-builder-form")?.addEventListener("submit", async (event) => {
      event.preventDefault();
      await launchSimulation();
    });
    byId("builder-validate")?.addEventListener("click", validateForm);
    byId("builder-reset")?.addEventListener("click", () => {
      if (builder.defaults) populateDefaults(builder.defaults);
      clearFieldErrors();
      setMessage("Defaults restored.", "ok");
    });
    byId("builder-stop")?.addEventListener("click", stopWorkflow);
  }

  function syncPair(rangeId, numberId) {
    const range = byId(rangeId);
    const number = byId(numberId);
    if (!range || !number) return;
    range.addEventListener("input", () => {
      number.value = range.value;
      updateSummary();
    });
    number.addEventListener("input", () => {
      const numeric = Number(number.value);
      const min = Number(range.min);
      const max = Number(range.max);
      if (Number.isFinite(numeric)) {
        range.value = String(Math.max(min, Math.min(max, numeric)));
      }
      updateSummary();
    });
  }

  function populateDefaults(payload) {
    const setup = payload.setup || {};
    const sim = payload.simulation || {};
    const workflow = payload.workflow || {};
    const choices = payload.choices || {};

    byId("builder-system").value = payload.system || "1L2Y";
    byId("builder-run-name").value = payload.run_name || "";
    fillSelect("builder-forcefield", choices.forcefields || [], setup.forcefield);
    fillSelect("builder-integrator", choices.integrators || [], sim.integrator);
    fillSelect("builder-platform", choices.platforms || [], sim.platform);
    fillSelect("builder-precision", choices.precisions || [], sim.precision);

    setPair("builder-ph-range", "builder-ph", setup.ph);
    setPair("builder-ion-range", "builder-ion", setup.ion_concentration_M);
    setPair("builder-padding-range", "builder-padding", setup.solvent_padding_nm);
    setPair("builder-nvt-range", "builder-nvt", sim.nvt_steps);
    setPair("builder-npt-range", "builder-npt", sim.npt_steps);
    setPair("builder-production-range", "builder-production", sim.production_steps);
    setPair("builder-temperature-range", "builder-temperature", sim.temperature_K);
    setPair("builder-timestep-range", "builder-timestep", sim.timestep_fs);
    setPair("builder-friction-range", "builder-friction", sim.friction_per_ps);

    byId("builder-water-model").value = setup.water_model || "auto";
    byId("builder-keep-heterogens").checked = !!setup.keep_heterogens;
    byId("builder-keep-water").checked = !!setup.keep_water;
    byId("builder-minimize").checked = sim.minimize !== false;
    byId("builder-trajectory-interval").value = sim.trajectory_interval_steps;
    byId("builder-telemetry-interval").value = sim.telemetry_interval;
    byId("builder-checkpoint-interval").value = sim.checkpoint_interval_steps;

    byId("builder-run-analysis").checked = workflow.run_analysis !== false;
    byId("builder-run-report").checked = workflow.run_report !== false;
    byId("builder-report-document").checked = workflow.report_document !== false;
    byId("builder-report-slides").checked = workflow.report_slides !== false;
    byId("builder-report-bundle").checked = workflow.report_bundle !== false;
    renderAnalysisChoices(choices.analyses || [], workflow.analyses || []);
    updateSummary();
  }

  function fillSelect(id, values, selected) {
    const element = byId(id);
    if (!element) return;
    element.innerHTML = "";
    values.forEach((value) => {
      const option = document.createElement("option");
      option.value = value;
      option.textContent = humanize(value);
      option.selected = value === selected;
      element.appendChild(option);
    });
  }

  function renderAnalysisChoices(analyses, selected) {
    const root = byId("builder-analysis-choices");
    if (!root) return;
    const chosen = new Set(selected);
    root.innerHTML = analyses.map((name) => `
      <label class="analysis-choice">
        <input type="checkbox" value="${escapeAttr(name)}" ${chosen.has(name) ? "checked" : ""}>
        <span>${escapeHTML(humanize(name))}</span>
      </label>`).join("");
  }

  function setPair(rangeId, numberId, value) {
    const range = byId(rangeId);
    const number = byId(numberId);
    if (range) {
      const min = Number(range.min);
      const max = Number(range.max);
      range.value = String(Math.max(min, Math.min(max, Number(value))));
    }
    if (number) number.value = String(value);
  }

  function collectPayload() {
    return {
      system: byId("builder-system")?.value.trim() || "",
      run_name: byId("builder-run-name")?.value.trim() || "",
      setup: {
        forcefield: byId("builder-forcefield")?.value || "charmm36",
        water_model: byId("builder-water-model")?.value.trim() || "auto",
        ph: numberValue("builder-ph", 7),
        ion_concentration_M: numberValue("builder-ion", 0.15),
        solvent_padding_nm: numberValue("builder-padding", 1),
        keep_heterogens: !!byId("builder-keep-heterogens")?.checked,
        keep_water: !!byId("builder-keep-water")?.checked,
      },
      simulation: {
        minimize: !!byId("builder-minimize")?.checked,
        nvt_steps: integerValue("builder-nvt", 250000),
        npt_steps: integerValue("builder-npt", 500000),
        production_steps: integerValue("builder-production", 1000000),
        timestep_fs: numberValue("builder-timestep", 2),
        temperature_K: numberValue("builder-temperature", 300),
        friction_per_ps: numberValue("builder-friction", 1),
        integrator: byId("builder-integrator")?.value || "langevin_middle",
        platform: byId("builder-platform")?.value || "auto",
        precision: byId("builder-precision")?.value || "mixed",
        trajectory_interval_steps: integerValue("builder-trajectory-interval", 1000),
        telemetry_interval: integerValue("builder-telemetry-interval", 1000),
        checkpoint_interval_steps: integerValue("builder-checkpoint-interval", 10000),
      },
      workflow: {
        run_analysis: !!byId("builder-run-analysis")?.checked,
        run_report: !!byId("builder-run-report")?.checked,
        report_document: !!byId("builder-report-document")?.checked,
        report_slides: !!byId("builder-report-slides")?.checked,
        report_bundle: !!byId("builder-report-bundle")?.checked,
        analyses: $$("#builder-analysis-choices input:checked").map((input) => input.value),
      },
    };
  }

  function updateSummary() {
    const payload = collectPayload();
    const sim = payload.simulation;
    const timestep = Number(sim.timestep_fs) || 0;
    const duration = (steps) => Number(steps || 0) * timestep / 1_000_000;
    const nvt = duration(sim.nvt_steps);
    const npt = duration(sim.npt_steps);
    const production = duration(sim.production_steps);
    const total = nvt + npt + production;
    const interval = Math.max(1, Number(sim.trajectory_interval_steps) || 1);
    const frames = Math.ceil(Math.max(0, sim.production_steps) / interval);

    setText("builder-nvt-duration", `${formatNumber(nvt)} ns`);
    setText("builder-npt-duration", `${formatNumber(npt)} ns`);
    setText("builder-production-duration", `${formatNumber(production)} ns`);
    setText("builder-summary-system", payload.system || "—");
    setText("builder-summary-output", payload.run_name || "auto-generated");
    setText("builder-summary-nvt", `${formatInteger(sim.nvt_steps)} steps · ${formatNumber(nvt)} ns`);
    setText("builder-summary-npt", `${formatInteger(sim.npt_steps)} steps · ${formatNumber(npt)} ns`);
    setText("builder-summary-production", `${formatInteger(sim.production_steps)} steps · ${formatNumber(production)} ns`);
    setText("builder-summary-total", `${formatNumber(total)} ns`);
    setText("builder-summary-frames", formatInteger(frames));
    setText("builder-summary-platform", `${sim.platform} · ${sim.precision}`);
  }

  async function validateForm() {
    clearFieldErrors();
    setMessage("Validating configuration…", "working");
    try {
      const result = await postJSON("/api/launcher/validate", collectPayload());
      builder.lastValidation = result;
      applyValidation(result);
      if (result.valid) setMessage("Configuration is valid and ready to launch.", "ok");
      return result;
    } catch (error) {
      const payload = error.payload || {};
      applyValidation(payload);
      setMessage(payload.error || "Please correct the highlighted fields.", "error");
      return payload;
    }
  }

  function applyValidation(result) {
    clearFieldErrors();
    Object.entries(result.errors || {}).forEach(([key, message]) => {
      const target = document.querySelector(`[data-error-for="${cssEscape(key)}"]`);
      if (target) target.textContent = message;
    });
    const warnings = result.warnings || [];
    const warningRoot = byId("builder-warning-list");
    if (warningRoot) {
      warningRoot.hidden = warnings.length === 0;
      warningRoot.innerHTML = warnings.map((warning) => `<div>${escapeHTML(warning)}</div>`).join("");
    }
    if (result.output) setText("builder-summary-output", result.output);
    if (Array.isArray(result.command)) {
      setText("builder-command-preview", formatCommand(result.command));
    }
  }

  async function launchSimulation() {
    if (builder.appState.process_running) {
      setMessage("A workflow is already running.", "error");
      return;
    }
    clearFieldErrors();
    setBusy(true);
    setMessage("Launching FastMDXplora…", "working");
    try {
      const result = await postJSON("/api/launcher/launch", collectPayload());
      applyValidation(result);
      if (!result.launched) throw makeRequestError("Launch was not accepted.", result);
      setMessage(`Workflow launched in ${result.output}`, "ok");
      applyAppState(result.state || {});
      window.FastMDXDashboard?.navigate("overview");
    } catch (error) {
      const payload = error.payload || {};
      applyValidation(payload);
      setMessage(payload.error || error.message || "Could not launch the workflow.", "error");
    } finally {
      setBusy(false);
    }
  }

  async function stopWorkflow() {
    if (!window.confirm("Request termination of the running FastMDXplora workflow?")) return;
    try {
      const result = await postJSON("/api/launcher/stop", {});
      applyAppState(result.state || {});
      setMessage(result.detail || "Termination requested.", result.stopped ? "warning" : "error");
    } catch (error) {
      setMessage(error.message || "Could not stop the workflow.", "error");
    }
  }

  async function refreshAppState() {
    try {
      applyAppState(await requestJSON("/api/app-state"));
    } catch (_) {
      // The main dashboard controller owns global connection reporting.
    }
  }

  function applyAppState(payload) {
    builder.appState = payload || {};
    const running = !!payload.process_running;
    const status = payload.status || (payload.active_run ? "run" : "idle");
    const stateRoot = byId("builder-run-state");
    if (stateRoot) stateRoot.dataset.state = status;
    setText("builder-state-text", stateLabel(payload));
    const dot = byId("builder-state-dot");
    if (dot) {
      dot.className = `status-dot ${running ? "status-dot-live" : status === "failed" ? "status-dot-error" : status === "completed" ? "status-dot-completed" : "status-dot-waiting"}`;
    }
    byId("builder-stop")?.toggleAttribute("hidden", !running);
    byId("builder-launch")?.toggleAttribute("disabled", running);
    if (payload.active_run) setText("builder-summary-output", payload.active_run);
  }

  function stateLabel(payload) {
    if (payload.process_running) return "Workflow running";
    if (payload.status === "completed") return "Workflow completed";
    if (payload.status === "failed") return `Workflow failed${payload.returncode != null ? ` (code ${payload.returncode})` : ""}`;
    if (payload.active_run) return "Run selected";
    return "No active workflow";
  }

  function setBusy(busy) {
    byId("builder-launch")?.toggleAttribute("disabled", busy || builder.appState.process_running);
    byId("builder-validate")?.toggleAttribute("disabled", busy);
  }

  function clearFieldErrors() {
    $$(".field-error").forEach((element) => { element.textContent = ""; });
  }

  function setMessage(message, kind) {
    const element = byId("builder-message");
    if (!element) return;
    element.textContent = message || "";
    element.dataset.kind = kind || "";
  }

  async function requestJSON(url) {
    const response = await fetch(url, {cache: "no-store", headers: {"X-FastMDX": "launcher"}});
    const payload = await response.json().catch(() => ({}));
    if (!response.ok) throw makeRequestError(payload.error || `HTTP ${response.status}`, payload);
    return payload;
  }

  async function postJSON(url, data) {
    const response = await fetch(url, {
      method: "POST",
      cache: "no-store",
      headers: {"Content-Type": "application/json", "X-FastMDX": "launcher"},
      body: JSON.stringify(data),
    });
    const payload = await response.json().catch(() => ({}));
    if (!response.ok) throw makeRequestError(payload.error || `HTTP ${response.status}`, payload);
    return payload;
  }

  function makeRequestError(message, payload) {
    const error = new Error(message);
    error.payload = payload;
    return error;
  }

  function formatCommand(command) {
    return command.map((part) => {
      const value = String(part);
      return /[\s"']/u.test(value) ? `"${value.replace(/"/g, '\\"')}"` : value;
    }).join(" ");
  }

  function numberValue(id, fallback) {
    const value = Number(byId(id)?.value);
    return Number.isFinite(value) ? value : fallback;
  }

  function integerValue(id, fallback) {
    const value = parseInt(byId(id)?.value, 10);
    return Number.isFinite(value) ? value : fallback;
  }

  function formatNumber(value) {
    const number = Number(value);
    if (!Number.isFinite(number)) return "—";
    if (number === 0) return "0";
    if (Math.abs(number) < 0.001) return number.toExponential(3);
    return number.toLocaleString(undefined, {maximumFractionDigits: 6});
  }

  function formatInteger(value) {
    const number = Number(value);
    return Number.isFinite(number) ? Math.round(number).toLocaleString() : "—";
  }

  function humanize(value) {
    return String(value || "").replace(/[_-]+/g, " ").replace(/\b\w/g, (letter) => letter.toUpperCase());
  }

  function setText(id, value) {
    const element = byId(id);
    if (element) element.textContent = value == null ? "" : String(value);
  }

  function escapeHTML(value) {
    return String(value == null ? "" : value).replace(/[&<>"']/g, (character) => ({
      "&": "&amp;", "<": "&lt;", ">": "&gt;", '"': "&quot;", "'": "&#39;",
    })[character]);
  }

  function escapeAttr(value) {
    return escapeHTML(value);
  }

  function cssEscape(value) {
    return String(value).replace(/["\\]/g, "\\$&");
  }

  window.FastMDXBuilder = {
    collectPayload,
    validate: validateForm,
    launch: launchSimulation,
    refresh: refreshAppState,
  };
}());
