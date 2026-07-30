/* FastMDXplora Live Dashboard — application controller.
 *
 * Framework-free and fully offline.  This module owns navigation, API
 * polling, status/results/file rendering, and the shared event bus used by
 * charts.js and molecule-viewer.js.  Every API renderer is isolated so one
 * malformed payload can never prevent the other dashboard sections from
 * updating.
 */

(function () {
  "use strict";

  const $ = (selector, root) => (root || document).querySelector(selector);
  const $$ = (selector, root) => Array.from((root || document).querySelectorAll(selector));
  const byId = (id) => document.getElementById(id);

  const state = {
    outputDir: "",
    paused: false,
    pollIntervalMs: 3000,
    lastUpdateMs: 0,
    refreshTimer: null,
    pages: ["builder", "overview", "live", "viewer", "analysis", "files", "settings"],
    activePage: "overview",
    appState: {},
    initialRouteResolved: false,
    ligandResname: null,
    bindingPocketCutoff: 5,
    playbackAvailable: false,
    playbackFrames: 0,
    playbackTotalFrames: 0,
    playbackFrameTimes: [],
    playbackSignature: null,
    metrics: [],
    status: {},
    health: {},
    results: {},
    structureInfo: null,
    setupManifest: {},
    simManifest: {},
    runId: "",
    runTitle: "FastMDXplora Live",
    apiErrors: {},
  };

  let chartHistorySamples = 600;

  document.addEventListener("DOMContentLoaded", boot);

  function boot() {
    const configuredRefresh = Number(document.body?.dataset.refreshSeconds || 3);
    state.pollIntervalMs = clampInt(configuredRefresh * 1000, 1000, 60000, 3000);
    wireNavigation();
    wireTopBar();
    wireViewerToggle();
    wireInfoTabs();
    wireTrajectoryControls();
    wireSettings();
    navigate(location.hash.replace(/^#/, "") || "overview", {updateHash: false});
    startLoadingChecklist();
    schedulePoll(0);
  }

  /* ------------------------------------------------------------------ */
  /* Navigation                                                          */
  /* ------------------------------------------------------------------ */
  function wireNavigation() {
    $$('[data-view-link]').forEach((element) => {
      element.addEventListener("click", (event) => {
        event.preventDefault();
        navigate(element.getAttribute("data-view-link"));
      });
    });
    window.addEventListener("hashchange", () => {
      const page = location.hash.replace(/^#/, "");
      if (state.pages.includes(page)) navigate(page, {updateHash: false});
    });
  }

  function navigate(page, options) {
    const opts = options || {};
    if (!state.pages.includes(page)) page = "overview";
    state.activePage = page;
    $$('.page').forEach((element) => {
      const hidden = element.getAttribute("data-page") !== page;
      element.hidden = hidden;
    });
    $$('[data-view-link]').forEach((element) => {
      element.classList.toggle(
        "active",
        element.getAttribute("data-view-link") === page
      );
    });
    document.documentElement.setAttribute("data-page", page);
    if (opts.updateHash !== false && location.hash !== `#${page}`) {
      history.replaceState(null, "", `#${page}`);
    }
    requestAnimationFrame(() => {
      if (page === "live") emit("live-page-opened", {});
      if (page === "viewer") emit("viewer-page-opened", {});
      if (page === "analysis" || page === "files") emit("results-page-opened", {});
    });
  }

  /* ------------------------------------------------------------------ */
  /* Loading screen                                                      */
  /* ------------------------------------------------------------------ */
  function startLoadingChecklist() {
    const checks = [
      ["telemetry", "/api/status"],
      ["structure", "/api/structure-info"],
      ["metrics", "/api/metrics"],
    ].map(([name, url]) => {
      setLoadingStep(name, "active");
      return fetchJSON(url)
        .then(() => setLoadingStep(name, "done"))
        .catch(() => setLoadingStep(name, "error"));
    });
    Promise.allSettled(checks).finally(() => {
      window.setTimeout(() => {
        document.body.classList.remove("state-loading");
        document.body.classList.add("state-ready");
        const loading = byId("loading-screen");
        if (loading) loading.setAttribute("aria-hidden", "true");
      }, 300);
    });
  }

  function setLoadingStep(name, status) {
    const element = $(`.loading-step[data-step="${name}"]`);
    if (element) element.setAttribute("data-state", status);
  }

  /* ------------------------------------------------------------------ */
  /* Controls                                                            */
  /* ------------------------------------------------------------------ */
  function wireTopBar() {
    byId("pause-toggle")?.addEventListener("click", () => {
      state.paused = !state.paused;
      byId("pause-toggle")?.setAttribute("aria-pressed", String(state.paused));
      setText("pause-label", state.paused ? "Resume Updates" : "Pause Updates");
      showToast(
        state.paused
          ? "Browser updates paused. The OpenMM simulation is still running."
          : "Browser updates resumed."
      );
      if (!state.paused) schedulePoll(0);
    });
    byId("refresh-now")?.addEventListener("click", () => schedulePoll(0));
    byId("open-output")?.addEventListener("click", async () => {
      try {
        const payload = await fetchJSON("/api/open-output");
        if (payload.opened) {
          showToast(`Opened output folder: ${payload.path}`);
        } else {
          await copyText(payload.path || state.outputDir || "");
          showToast("Could not open the folder automatically; its path was copied.", "warning");
        }
      } catch (error) {
        if (state.outputDir) await copyText(state.outputDir);
        showToast("Could not open the output folder; its path was copied.", "warning");
      }
    });
  }

  function wireViewerToggle() {
    byId("mini-preview-frame")?.addEventListener("click", () => navigate("viewer"));
  }

  function wireInfoTabs() {
    $$('.info-tab').forEach((tab) => {
      tab.addEventListener("click", () => {
        const name = tab.getAttribute("data-tab");
        $$('.info-tab').forEach((item) => item.classList.toggle("active", item === tab));
        $$('.info-pane').forEach((pane) => {
          const active = pane.getAttribute("data-tab") === name;
          pane.hidden = !active;
          pane.classList.toggle("active", active);
        });
      });
    });
  }

  function wireTrajectoryControls() {
    byId("traj-slider")?.addEventListener("input", (event) => {
      emit("trajectory-seek", {frame: parseInt(event.target.value, 10) || 0});
    });
    $$('[data-traj]').forEach((button) => {
      button.addEventListener("click", () => {
        emit("trajectory-action", {action: button.getAttribute("data-traj")});
      });
    });
  }

  function wireSettings() {
    const ids = [
      "setting-protein-rep", "setting-ligand-rep", "setting-background",
      "setting-show-water", "setting-show-ions", "setting-spin", "setting-fog",
      "setting-preserve-camera", "setting-refresh-seconds", "setting-chart-history",
      "setting-compact", "setting-reduced-motion", "setting-advanced-metrics",
      "setting-time-format", "setting-run-name", "setting-ligand-resname",
      "setting-pocket-cutoff", "setting-scinote",
    ];
    ids.forEach((id) => byId(id)?.addEventListener("change", onSettingsChanged));
    byId("pocket-cutoff")?.addEventListener("change", (event) => {
      const value = clampFloat(parseFloat(event.target.value), 3, 15, 5);
      state.bindingPocketCutoff = value;
      if (byId("setting-pocket-cutoff")) byId("setting-pocket-cutoff").value = String(value);
      onSettingsChanged();
    });
  }

  function onSettingsChanged() {
    state.pollIntervalMs = clampInt(
      parseFloat(byId("setting-refresh-seconds")?.value) * 1000,
      1000,
      60000,
      3000
    );
    state.bindingPocketCutoff = clampFloat(
      parseFloat(byId("setting-pocket-cutoff")?.value),
      3,
      15,
      5
    );
    const ligand = (byId("setting-ligand-resname")?.value || "").trim().toUpperCase();
    state.ligandResname = ligand || null;
    chartHistorySamples = clampInt(
      parseInt(byId("setting-chart-history")?.value, 10),
      60,
      5000,
      600
    );
    const customRunName = (byId("setting-run-name")?.value || "").trim();
    if (customRunName) state.runTitle = customRunName;

    document.body.classList.toggle("compact-mode", !!byId("setting-compact")?.checked);
    document.body.classList.toggle("reduced-motion", !!byId("setting-reduced-motion")?.checked);
    document.body.classList.toggle(
      "advanced-metrics",
      !!byId("setting-advanced-metrics")?.checked
    );

    renderTopBar(state.status, state.health);
    emit("settings-updated", {
      ligand: state.ligandResname,
      pocketCutoff: state.bindingPocketCutoff,
      chartHistory: chartHistorySamples,
      proteinRepresentation: byId("setting-protein-rep")?.value || "cartoon",
      ligandRepresentation: byId("setting-ligand-rep")?.value || "sticks",
      background: byId("setting-background")?.value || "matte-black",
      showWater: !!byId("setting-show-water")?.checked,
      showIons: !!byId("setting-show-ions")?.checked,
      spin: !!byId("setting-spin")?.checked,
      fog: !!byId("setting-fog")?.checked,
      preserveCamera: byId("setting-preserve-camera")?.checked !== false,
    });
    schedulePoll(0);
  }

  /* ------------------------------------------------------------------ */
  /* Polling                                                             */
  /* ------------------------------------------------------------------ */
  async function fetchJSON(url) {
    const response = await fetch(url, {
      cache: "no-store",
      headers: {"X-FastMDX": "live"},
    });
    if (!response.ok) throw new Error(`HTTP ${response.status} on ${url}`);
    return response.json();
  }

  function schedulePoll(delayMs) {
    if (state.refreshTimer) window.clearTimeout(state.refreshTimer);
    state.refreshTimer = window.setTimeout(async () => {
      if (!state.paused) {
        try {
          await pollNow();
        } catch (error) {
          console.warn("dashboard poll error", error);
        }
      }
      schedulePoll(state.pollIntervalMs);
    }, Math.max(0, delayMs));
  }

  async function pollNow() {
    const names = ["app", "status", "metrics", "events", "results", "structure", "playback"];
    const urls = [
      "/api/app-state", "/api/status", "/api/metrics", "/api/events", "/api/results",
      "/api/structure-info", "/api/playback-info",
    ];
    const settled = await Promise.allSettled(urls.map(fetchJSON));
    let successes = 0;

    settled.forEach((result, index) => {
      const name = names[index];
      if (result.status === "rejected") {
        state.apiErrors[name] = String(result.reason || "request failed");
        console.warn(`dashboard ${name} request failed`, result.reason);
        return;
      }
      delete state.apiErrors[name];
      successes += 1;
      const payload = result.value;
      if (name === "app") safeApply(name, () => applyAppState(payload));
      if (name === "status") safeApply(name, () => applyStatus(payload));
      if (name === "metrics") safeApply(name, () => applyMetrics(payload.metrics || []));
      if (name === "events") safeApply(name, () => renderEvents(payload.events || []));
      if (name === "results") safeApply(name, () => applyResults(payload));
      if (name === "structure") safeApply(name, () => applyStructure(payload));
      if (name === "playback") safeApply(name, () => applyPlayback(payload));
    });

    if (successes) {
      state.lastUpdateMs = Date.now();
      updateRefreshedAt();
    }
    updateConnectionState();
  }

  function applyAppState(payload) {
    state.appState = payload || {};
    const activeRun = payload?.active_run || "";
    if (activeRun) state.outputDir = activeRun;

    if (!state.initialRouteResolved) {
      state.initialRouteResolved = true;
      const requested = location.hash.replace(/^#/, "");
      if (!requested && !activeRun) navigate("builder");
    }

    if (!activeRun) {
      setText("topbar-run-id", "workspace");
      setText("topbar-run-title", "No active run");
      setText("topbar-stage", "configure a simulation");
      setText("sidebar-run-name", "No active run");
      setText("sidebar-platform", "—");
    }
    emit("app-state", payload || {});
  }

  function safeApply(name, callback) {
    try {
      callback();
    } catch (error) {
      state.apiErrors[`render-${name}`] = String(error);
      console.error(`dashboard ${name} renderer failed`, error);
      showToast(`${humanise(name)} data could not be rendered. See the browser console.`, "warning");
    }
  }

  function updateConnectionState() {
    const ageSeconds = state.lastUpdateMs ? (Date.now() - state.lastUpdateMs) / 1000 : Infinity;
    const requestFailed = Object.keys(state.apiErrors).length > 0;
    if (ageSeconds > 30 || requestFailed) {
      byId("topbar-status-dot")?.classList.add("status-dot-stale");
    }
  }

  function updateRefreshedAt() {
    const date = new Date(state.lastUpdateMs);
    const use24 = byId("setting-time-format")?.value !== "12h";
    const text = use24
      ? `${date.toLocaleDateString()} ${date.toLocaleTimeString([], {hour12: false})}`
      : date.toLocaleTimeString();
    setText("refreshed-at", text);
  }

  /* ------------------------------------------------------------------ */
  /* Status                                                              */
  /* ------------------------------------------------------------------ */
  function applyStatus(payload) {
    const status = payload.status || {};
    const health = payload.health || {};
    state.status = status;
    state.health = health;
    renderTopBar(status, health);
    renderHero(status);
    renderHealth(health);
    renderStageTimeline(status);
    renderLiveProgress(status);
    emit("status-updated", {status, health});
  }

  function renderTopBar(status, health) {
    if (!state.appState?.active_run) {
      setClassName("topbar-status-dot", "status-dot status-dot-waiting");
      setClassName("sidebar-status-dot", "status-dot status-dot-waiting");
      setText("topbar-status-text", "ready");
      setText("topbar-stage", "configure a simulation");
      setText("topbar-step", "—");
      setText("topbar-total", "—");
      setText("topbar-platform", "—");
      setText("topbar-temperature", "—");
      setText("sidebar-connection-state", "ready");
      setText("sidebar-platform", "—");
      setText("sidebar-run-name", "No active run");
      setText("topbar-run-id", "workspace");
      setText("topbar-run-title", "FastMDXplora Dashboard");
      return;
    }
    const statusName = String(health.state || status.status || "waiting").toLowerCase();
    const dotClass = stateDotClass(statusName);
    setClassName("topbar-status-dot", `status-dot ${dotClass}`);
    setClassName("sidebar-status-dot", `status-dot ${dotClass}`);
    setText("topbar-status-text", statusName);
    setText("topbar-stage", status.stage || "not available");
    setText("topbar-step", valueOrDash(status.current_step));
    setText("topbar-total", valueOrDash(status.total_planned_steps));
    setText("topbar-platform", status.platform || state.simManifest.platform || "—");

    const latest = state.metrics[state.metrics.length - 1] || {};
    const temperature = firstPresent(
      latest.temperature,
      status.target_temperature_K,
      state.simManifest.temperature_K
    );
    setText("topbar-temperature", temperature != null ? `${formatNumber(temperature, 1)} K` : "—");
    setText("sidebar-connection-state", statusName);
    setText("sidebar-platform", status.platform || state.simManifest.platform || "not available");
    setText("sidebar-run-name", state.runTitle);

    state.runId = status.system_id || state.results?.system?.system || state.runId;
    setText("topbar-run-id", state.runId || "system");
    setText("topbar-run-title", state.runTitle);
  }

  function stateDotClass(value) {
    if (["ok", "live", "completed"].includes(value)) return "status-dot-live";
    if (["failed", "error", "critical"].includes(value)) return "status-dot-error";
    return "status-dot-stale";
  }

  function renderHero(status) {
    const card = byId("hero-card");
    if (card) card.setAttribute("data-state", String(status.status || "running").toLowerCase());
    setText("hero-status-text", humanise(status.status || status.stage || "waiting"));
    setText("hero-stage", status.stage || "—");
    const pct = progressPercent(status);
    setWidth("hero-progress-fill", pct);
    setText("hero-progress-pct", pct != null ? pct.toFixed(1) : "—");
    setText(
      "hero-sim-time",
      status.simulation_time_completed_ns != null
        ? formatNumber(status.simulation_time_completed_ns, 3)
        : "—"
    );
    setText(
      "hero-elapsed",
      status.elapsed_wall_time_s != null ? fmtDuration(status.elapsed_wall_time_s) : "—"
    );
    setText("hero-eta", computeETA(status));
    setText(
      "hero-step",
      status.current_step != null
        ? `${status.current_step} / ${status.total_planned_steps ?? "—"}`
        : "—"
    );
  }

  function renderHealth(health) {
    const stateName = String(health.state || "unknown").toLowerCase();
    byId("hero-health")?.setAttribute("data-state", stateName);
    setText("health-headline", health.message || humanise(stateName));
    setText("health-explanation", health.explanation || "");
    setText("health-pill", stateName);
    byId("health-pill")?.setAttribute("data-state", stateName);

    const list = byId("health-list");
    if (list) {
      const items = Array.isArray(health.items) ? health.items.slice(0, 4) : [];
      list.innerHTML = items.map((item) => `
        <li data-state="${escapeAttr(item.severity || "ok")}">
          <strong>${escapeHTML(item.title || item.severity || "info")}</strong>
          <span class="muted small"> — ${escapeHTML(item.detail || "")}</span>
        </li>
      `).join("");
    }
    setText("live-health-pill", stateName);
    byId("live-health-pill")?.setAttribute("data-state", stateName);
    setText("live-health-message", health.message || "not available");
    setText("live-health-explanation", health.explanation || "");
  }

  function renderStageTimeline(status) {
    const order = ["setup", "minimization", "nvt", "npt", "production", "analysis", "report"];
    const phaseMap = {};
    (state.results.phases || []).forEach((phase) => {
      phaseMap[String(phase.name || "").toLowerCase()] = String(phase.status || "").toLowerCase();
    });
    const liveStates = status?.stage_states && typeof status.stage_states === "object"
      ? status.stage_states : {};
    const current = normaliseStage(status.stage);
    const currentIndex = order.indexOf(current);
    const simulationDone = isPhaseDone(phaseMap.simulation);

    $$('.stage-step').forEach((element) => {
      const stage = element.getAttribute("data-stage");
      let stageState = phaseVisualState(liveStates[stage]);

      // Older/completed runs may not have the new live stage map, so retain
      // manifest-based fallback without overwriting a real current/failed state.
      if (stageState === "waiting") {
        if (stage === "setup") stageState = phaseVisualState(phaseMap.setup);
        if (["minimization", "nvt", "npt", "production"].includes(stage) && simulationDone) {
          stageState = "completed";
        }
        if (stage === "analysis") stageState = phaseVisualState(phaseMap.analysis);
        if (stage === "report") stageState = phaseVisualState(phaseMap.report);
      }

      if (currentIndex >= 0 && stageState === "waiting") {
        const index = order.indexOf(stage);
        if (index < currentIndex) stageState = "completed";
        if (index === currentIndex) stageState = "current";
      }
      if (stage === current && phaseVisualState(liveStates[stage]) === "current") {
        stageState = "current";
      }
      element.setAttribute("data-state", stageState);
    });
  }

  function renderLiveProgress(status) {
    const pct = progressPercent(status);
    setWidth("live-progress-fill", pct);
    setText("live-progress-pct", pct != null ? pct.toFixed(1) : "—");
    setText(
      "live-sim-time",
      status.simulation_time_completed_ns != null
        ? formatNumber(status.simulation_time_completed_ns, 3)
        : "—"
    );
    setText("live-stage-cell", status.stage || "not available");
    setText("live-step-cell", present(status.current_step));
    setText("live-total-cell", present(status.total_planned_steps));
    setText(
      "live-frames-cell",
      status.current_frame_count != null
        ? `${status.current_frame_count}${status.planned_frame_count != null ? ` / ${status.planned_frame_count}` : ""}`
        : "not available"
    );
    setText(
      "live-simtime-cell",
      status.simulation_time_completed_ns != null
        ? `${formatNumber(status.simulation_time_completed_ns, 6)} ns`
        : "not available"
    );
    setText(
      "live-elapsed-cell",
      status.elapsed_wall_time_s != null ? fmtDuration(status.elapsed_wall_time_s) : "not available"
    );
    setText("live-eta-cell", computeETA(status));
    setText("live-checkpoint-cell", status.current_checkpoint_path || "not available");
    setText("live-lastupdate-cell", formatTimestamp(status.last_update_timestamp));
    setText("live-card-step", valueOrDash(status.current_step));
  }

  /* ------------------------------------------------------------------ */
  /* Metrics and events                                                  */
  /* ------------------------------------------------------------------ */
  function applyMetrics(metrics) {
    state.metrics = (Array.isArray(metrics) ? metrics : [])
      .map(normaliseMetricRow)
      .slice(-chartHistorySamples);
    renderMetricCards();
    if (window.FastMDXCharts) window.FastMDXCharts.update(state.metrics);
    emit("metrics-updated", {metrics: state.metrics});
  }

  function normaliseMetricRow(row) {
    const output = Object.assign({}, row || {});
    const aliases = {
      potentialEnergy: "potential_energy",
      kineticEnergy: "kinetic_energy",
      totalEnergy: "total_energy",
      simulationSpeed: "speed",
      frame: "current_frame_count",
      frames: "current_frame_count",
    };
    Object.keys(aliases).forEach((key) => {
      if (output[aliases[key]] == null && output[key] != null) output[aliases[key]] = output[key];
    });
    return output;
  }

  function renderMetricCards() {
    const latest = state.metrics[state.metrics.length - 1] || {};
    const units = {
      potential_energy: "kJ/mol",
      temperature: "K",
      density: "g/mL",
      speed: "ns/day",
      frames: "frames",
      pressure: "bar",
    };
    $$('.metric-card').forEach((card) => {
      const key = card.getAttribute("data-metric");
      const derived = deriveMetric(key, latest);
      const value = card.querySelector("[data-value]");
      const unit = card.querySelector(".metric-card-unit");
      if (value) value.textContent = derived != null ? derived : "—";
      if (unit && units[key]) unit.textContent = units[key];
      const previous = card.getAttribute("data-last-value");
      if (derived != null && previous && previous !== String(derived)) {
        card.setAttribute("data-pulse", "");
        window.setTimeout(() => card.removeAttribute("data-pulse"), 900);
      }
      card.setAttribute("data-last-value", derived != null ? String(derived) : "");
    });
  }

  function deriveMetric(key, latest) {
    if (key === "frames") {
      const value = firstPresent(
        latest.current_frame_count,
        state.status.current_frame_count,
        state.simManifest.n_production_frames
      );
      return value != null ? String(value) : null;
    }
    const value = latest[key];
    const number = Number(value);
    if (!Number.isFinite(number)) return null;
    if (key === "speed") return number.toFixed(2);
    if (key === "density") return number.toFixed(4);
    return number.toFixed(2);
  }

  function renderEvents(events) {
    const list = byId("events-list");
    if (!list) return;
    const items = (Array.isArray(events) ? events : []).slice(-50).reverse();
    list.innerHTML = items.length ? items.map((event) => `
      <li data-level="${escapeAttr(event.level || "info")}">
        <span class="level-badge">${escapeHTML(event.level || "info")}</span>
        <span class="event-message">${escapeHTML(event.message || "")}</span>
        <span class="event-time">${escapeHTML(formatEventTime(event.timestamp))}</span>
      </li>
    `).join("") : `
      <li data-level="info">
        <span class="level-badge">info</span>
        <span class="event-message">No telemetry events were recorded.</span>
        <span class="event-time"></span>
      </li>`;
  }

  /* ------------------------------------------------------------------ */
  /* Results / analyses / files                                          */
  /* ------------------------------------------------------------------ */
  function applyResults(payload) {
    state.results = payload || {};
    state.outputDir = payload.output_dir || state.outputDir;
    state.setupManifest = payload.setup || {};
    state.simManifest = payload.simulation || {};
    state.runTitle = (byId("setting-run-name")?.value || "").trim()
      || payload.run_title
      || state.runTitle;
    state.runId = payload.system?.system || state.runId;

    renderTopBar(state.status, state.health);
    renderStageTimeline(state.status);
    renderAnalysis(payload);
    renderFiles(payload);
    renderRunSummary(payload);
    renderSimulationTab(state.structureInfo || {});
    emit("results-updated", payload);
  }

  function renderAnalysis(payload) {
    const grid = byId("analysis-grid");
    const empty = byId("analysis-empty");
    if (!grid || !empty) return;

    const analyses = Array.isArray(payload.analyses) ? payload.analyses : [];
    const plots = Array.isArray(payload.plots) ? payload.plots : [];
    const cards = [];
    const usedPlotPaths = new Set();

    analyses.forEach((analysis) => {
      const plot = analysis.plot || null;
      if (plot?.path) usedPlotPaths.add(plot.path);
      cards.push(analysisCardHtml(analysis, plot));
    });
    plots.forEach((plot) => {
      if (usedPlotPaths.has(plot.path)) return;
      cards.push(analysisCardHtml({
        name: plot.title,
        title: plot.title,
        status: "complete",
        message: plot.category || "",
        artifacts: [plot],
      }, plot));
    });

    grid.innerHTML = cards.join("");
    const svgBundle = byId("download-all-svg");
    const svgCount = Number(payload.svg_figure_count || 0);
    if (svgBundle) {
      svgBundle.hidden = svgCount < 1;
      svgBundle.href = payload.svg_bundle_href || "/analysis-figures-svg.zip";
      svgBundle.textContent = svgCount === 1
        ? "Download SVG figure"
        : `Download all ${svgCount} SVG figures`;
    }
    const count = analyses.length || plots.length;
    setText(
      "analysis-meta",
      count ? `${count} analysis output${count === 1 ? "" : "s"}` : "no analyses yet"
    );
    if (cards.length) empty.setAttribute("hidden", "");
    else empty.removeAttribute("hidden");
  }

  function analysisCardHtml(analysis, plot) {
    const status = String(analysis.status || "unknown").toLowerCase();
    const title = analysis.title || humanise(analysis.name || "Analysis");
    const artifacts = Array.isArray(analysis.artifacts) ? analysis.artifacts : [];
    const imageArtifact = artifacts.find((item) => /\.(png|jpe?g|gif)$/i.test(item?.path || item?.name || ""));
    const svgArtifact = artifacts.find((item) => /\.svg$/i.test(item?.path || item?.name || ""));
    const displayPlot = plot?.href ? plot : imageArtifact || svgArtifact || null;
    const vectorPlot = plot?.svg_href
      ? {href: plot.svg_href, download_href: plot.svg_download_href}
      : svgArtifact;
    const dataArtifact = artifacts.find((item) => !/\.(png|jpe?g|gif|svg)$/i.test(item?.path || item?.name || ""));
    const primary = dataArtifact || artifacts[0] || displayPlot;
    const image = displayPlot?.href
      ? `<div class="ac-frame"><img src="${escapeAttr(displayPlot.href)}" alt="${escapeAttr(title)}" loading="lazy"></div>`
      : `<div class="ac-frame ac-no-plot"><div class="muted">No saved figure file was found for this analysis. Re-run the analysis with this revised version to generate PNG and SVG figures.</div></div>`;
    const links = [];
    if (displayPlot?.href) {
      links.push(`<a class="file-action" href="${escapeAttr(displayPlot.href)}" target="_blank" rel="noopener">Open figure</a>`);
      if (!/\.svg(?:$|\?)/i.test(displayPlot.href)) {
        links.push(`<a class="file-action" href="${escapeAttr(displayPlot.download_href || `${displayPlot.href}${displayPlot.href.includes("?") ? "&" : "?"}download=1`)}" download>Download PNG</a>`);
      }
    }
    if (vectorPlot?.href) {
      links.push(`<a class="file-action" href="${escapeAttr(vectorPlot.download_href || vectorPlot.href)}" download>Download SVG</a>`);
    }
    if (primary?.href && primary?.href !== displayPlot?.href && primary?.href !== vectorPlot?.href) {
      links.push(`<a class="file-action" href="${escapeAttr(primary.href)}" target="_blank" rel="noopener">Open data</a>`);
    }
    return `
      <article class="analysis-card" data-state="${escapeAttr(status)}">
        <div class="ac-header">
          <div class="ac-title">${escapeHTML(title)}</div>
          <div class="ac-status">${escapeHTML(status)}</div>
        </div>
        ${image}
        <div class="ac-body">${escapeHTML(analysis.message || plot?.category || "")}</div>
        <div class="ac-footer">${links.join("") || '<span class="muted small">Artifacts will appear when available.</span>'}</div>
      </article>`;
  }

  function renderFiles(payload) {
    const artifacts = Array.isArray(payload.artifacts) ? payload.artifacts : [];
    const groups = {
      "reports-files": Array.isArray(payload.reports) ? payload.reports : artifacts.filter((item) => item.path.startsWith("report/")),
      "simulation-files": artifacts.filter((item) => item.path.startsWith("simulation/")),
      "analysis-files": artifacts.filter((item) => item.path.startsWith("analysis/")),
    };
    Object.entries(groups).forEach(([id, files]) => {
      const root = byId(id);
      if (!root) return;
      root.innerHTML = files.length
        ? files.map(fileRowHtml).join("")
        : '<div class="muted small">No files are available in this section yet.</div>';
    });
    wireCopyActions();
  }

  function fileRowHtml(file) {
    const title = file.name || file.path || "Artifact";
    const size = file.size != null ? humanSize(parseInt(file.size, 10)) : "—";
    const mtime = file.mtime != null
      ? new Date(parseFloat(file.mtime) * 1000).toLocaleString()
      : "—";
    const href = file.href || `/artifacts/${encodeURI(file.path || "")}`;
    const downloadHref = file.download_href || `${href}${href.includes("?") ? "&" : "?"}download=1`;
    return `
      <div class="file-row" data-path="${escapeAttr(file.absolute_path || file.path || "")}">
        <div class="file-title" title="${escapeAttr(file.path || "")}">${escapeHTML(title)}</div>
        <div class="file-meta">
          <span>${escapeHTML(size)}</span>
          <span class="muted">${escapeHTML(mtime)}</span>
          <div class="file-actions">
            <a class="file-action" href="${escapeAttr(href)}" target="_blank" rel="noopener">Open</a>
            <a class="file-action" href="${escapeAttr(downloadHref)}" download>Download</a>
            <button class="file-action" type="button" data-copy-path>Copy path</button>
          </div>
        </div>
      </div>`;
  }

  function wireCopyActions() {
    $$('[data-copy-path]').forEach((button) => {
      button.addEventListener("click", async () => {
        const row = button.closest(".file-row");
        const path = row?.getAttribute("data-path") || "";
        const copied = await copyText(path);
        showToast(copied ? "File path copied." : "Could not copy the file path.", copied ? "ok" : "warning");
      });
    });
  }

  function renderRunSummary(payload) {
    const card = byId("summary-card");
    if (!card) return;
    const analyses = Array.isArray(payload.analyses) ? payload.analyses : [];
    const plots = Array.isArray(payload.plots) ? payload.plots : [];
    const reports = Array.isArray(payload.reports) ? payload.reports : [];
    const hasResults = payload.has_report || payload.has_analysis || reports.length || analyses.length;
    card.hidden = !hasResults;
    if (!hasResults) return;
    const stats = [
      {label: "Analyses", value: String(analyses.length)},
      {label: "Figures", value: String(plots.length)},
      {label: "Report formats", value: String(reports.length)},
      {label: "Result bundle", value: reports.some((item) => item.path?.endsWith(".zip")) ? "ready" : "—"},
    ];
    const grid = byId("summary-grid");
    if (grid) grid.innerHTML = stats.map((item) => `
      <div class="summary-stat">
        <span class="summary-stat-label">${escapeHTML(item.label)}</span>
        <span class="summary-stat-value">${escapeHTML(item.value)}</span>
      </div>`).join("");
  }

  /* ------------------------------------------------------------------ */
  /* Structure and playback                                              */
  /* ------------------------------------------------------------------ */
  function applyStructure(info) {
    state.structureInfo = info || {};
    renderStructureTab(state.structureInfo);
    if (state.structureInfo.ligand_resnames?.length && !state.ligandResname) {
      state.ligandResname = state.structureInfo.ligand_resnames[0];
    }
    renderLigandTab(state.structureInfo);
    renderSimulationTab(state.structureInfo);
    emit("structure-updated", state.structureInfo);
  }

  function renderStructureTab(info) {
    const body = byId("structure-tab-tbody");
    if (!body) return;
    if (!info?.valid) {
      const reason = info?.reason ? ` (${humanise(info.reason)})` : "";
      body.innerHTML = `<tr><td colspan="2" class="muted">Structure metadata is not available${escapeHTML(reason)}.</td></tr>`;
      return;
    }
    body.innerHTML = `
      <tr><th>Protein chains</th><td>${escapeHTML(info.n_chains)}</td></tr>
      <tr><th>Protein residues</th><td>${escapeHTML(info.protein_residues)}</td></tr>
      <tr><th>Protein atoms</th><td>${escapeHTML(info.protein_atoms)}</td></tr>
      <tr><th>Ligands</th><td>${escapeHTML((info.ligand_resnames || []).join(", ") || "none")}</td></tr>
      <tr><th>Water</th><td>${escapeHTML(info.water_residues)}</td></tr>
      <tr><th>Ions</th><td>${escapeHTML(info.ions)}</td></tr>`;
  }

  function renderSimulationTab() {
    const body = byId("simulation-tab-tbody");
    if (!body) return;
    const setup = state.setupManifest || {};
    const simulation = state.simManifest || {};
    const status = state.status || {};
    body.innerHTML = `
      <tr><th>Force field</th><td>${escapeHTML(setup.force_field || status.force_field || "—")}</td></tr>
      <tr><th>Water model</th><td>${escapeHTML(setup.water_model || status.water_model || "—")}</td></tr>
      <tr><th>pH</th><td>${escapeHTML(setup.ph ?? "—")}</td></tr>
      <tr><th>Ion concentration</th><td>${escapeHTML(setup.ion_concentration_M != null ? `${setup.ion_concentration_M} M` : "—")}</td></tr>
      <tr><th>Temperature</th><td>${escapeHTML(firstPresent(status.target_temperature_K, simulation.temperature_K) != null ? `${firstPresent(status.target_temperature_K, simulation.temperature_K)} K` : "—")}</td></tr>
      <tr><th>Timestep</th><td>${escapeHTML(firstPresent(status.timestep_fs, simulation.timestep_fs) != null ? `${firstPresent(status.timestep_fs, simulation.timestep_fs)} fs` : "—")}</td></tr>
      <tr><th>Precision</th><td>${escapeHTML(simulation.precision || status.precision || "—")}</td></tr>
      <tr><th>Platform</th><td>${escapeHTML(status.platform || simulation.platform || "—")}</td></tr>`;
  }

  function renderLigandTab(info) {
    const tools = byId("ligand-tools");
    const meta = byId("ligand-meta");
    const body = byId("ligand-tab-tbody");
    if (!tools || !meta || !body) return;
    const instances = Array.isArray(info?.ligand_instances) ? info.ligand_instances : [];
    const instance = instances.find((item) => item.resname === state.ligandResname) || instances[0];
    if (!instance) {
      tools.hidden = true;
      meta.textContent = "not available";
      body.innerHTML = '<tr><td colspan="2" class="muted">No ligand was detected.</td></tr>';
      return;
    }
    tools.hidden = false;
    const residueName = instance.resname || state.ligandResname;
    const atoms = info?.atoms_by_resname?.[residueName] ?? "—";
    state.ligandResname = residueName;
    meta.textContent = `${residueName} · chain ${instance.chain || "—"} · resi ${instance.resi || "—"} · ${atoms} atoms`;
    body.innerHTML = `
      <tr><th>Ligand</th><td>${escapeHTML(residueName)}</td></tr>
      <tr><th>Chain</th><td>${escapeHTML(instance.chain || "—")}</td></tr>
      <tr><th>Residue ID</th><td>${escapeHTML(instance.resi || "—")}</td></tr>
      <tr><th>Atom count</th><td>${escapeHTML(atoms)}</td></tr>
      <tr><th>Nearby residues</th><td>Use “Show pocket residues”</td></tr>
      <tr><th>Pocket distance</th><td>${escapeHTML(state.bindingPocketCutoff)} Å cutoff</td></tr>
      <tr><th>H-bonds</th><td>Requires analysis output</td></tr>
      <tr><th>Hydrophobic contacts</th><td>Requires analysis output</td></tr>
      <tr><th>Salt bridges</th><td>Requires analysis output</td></tr>`;
  }

  function applyPlayback(payload) {
    const previousSignature = state.playbackSignature;
    const signature = payload?.source_signature || payload?.compiled_at || null;
    const previousFrame = parseInt(byId("traj-slider")?.value || "0", 10) || 0;

    state.playbackAvailable = !!payload?.playback_available;
    state.playbackFrames = Number(payload?.n_frames_browser || 0);
    state.playbackTotalFrames = Number(payload?.n_frames_total || 0);
    state.playbackFrameTimes = Array.isArray(payload?.frame_times_ns) ? payload.frame_times_ns : [];
    state.playbackSignature = signature;

    const slider = byId("traj-slider");
    const row = byId("trajectory-row");
    if (!slider || !row) return;
    if (!state.playbackAvailable) {
      row.hidden = true;
      emit("playback-ready", payload || {});
      return;
    }

    const maxFrame = Math.max(0, state.playbackFrames - 1);
    const frame = Math.min(previousFrame, maxFrame);
    slider.min = "0";
    slider.max = String(maxFrame);
    // Do not reset the user's scrubber on every three-second API poll.
    slider.value = String(frame);
    setText("traj-current", String(frame));
    setText("traj-total", String(state.playbackFrames));
    setText(
      "traj-simtime",
      state.playbackFrameTimes[frame] != null
        ? formatNumber(state.playbackFrameTimes[frame], 3) : "—"
    );
    row.hidden = false;
    emit("playback-ready", Object.assign({}, payload, {
      source_changed: previousSignature !== null && previousSignature !== signature,
    }));
  }

  /* ------------------------------------------------------------------ */
  /* Helpers                                                             */
  /* ------------------------------------------------------------------ */
  function progressPercent(status) {
    if (Number.isFinite(Number(status.progress_percent))) {
      return clampFloat(Number(status.progress_percent), 0, 100, null);
    }
    if (status.current_step != null && status.total_planned_steps != null && Number(status.total_planned_steps) > 0) {
      return clampFloat(Number(status.current_step) / Number(status.total_planned_steps) * 100, 0, 100, null);
    }
    return null;
  }

  function fmtDuration(seconds) {
    const total = Math.max(0, Math.floor(Number(seconds) || 0));
    const hours = Math.floor(total / 3600);
    const minutes = Math.floor((total % 3600) / 60);
    const secs = total % 60;
    return hours ? `${hours}h ${minutes}m` : `${minutes}m ${secs}s`;
  }

  function computeETA(status) {
    const step = Number(status.current_step);
    const total = Number(status.total_planned_steps);
    const elapsed = Number(status.elapsed_wall_time_s);
    if (!(step > 0) || !(total > 0) || !(elapsed >= 0) || step >= total) return step >= total ? "0m 0s" : "—";
    return fmtDuration(elapsed * (total / step - 1));
  }

  function normaliseStage(value) {
    const stage = String(value || "").toLowerCase();
    if (stage.includes("minim")) return "minimization";
    if (stage.includes("nvt")) return "nvt";
    if (stage.includes("npt")) return "npt";
    if (stage.includes("production")) return "production";
    if (stage.includes("analysis")) return "analysis";
    if (stage.includes("report")) return "report";
    if (stage.includes("setup") || stage.includes("loading")) return "setup";
    return stage;
  }

  function phaseVisualState(value) {
    const status = String(value || "").toLowerCase();
    if (["ok", "complete", "completed", "success", "succeeded"].includes(status)) return "completed";
    if (["error", "failed"].includes(status)) return "failed";
    if (["skipped", "not run"].includes(status)) return "skipped";
    if (["running", "active", "current"].includes(status)) return "current";
    return "waiting";
  }

  function isPhaseDone(value) {
    return ["ok", "complete", "completed", "success"].includes(String(value || "").toLowerCase());
  }

  function formatTimestamp(value) {
    if (!value) return "not available";
    const date = new Date(value);
    return Number.isNaN(date.getTime()) ? String(value) : date.toLocaleString();
  }

  function formatEventTime(value) {
    if (!value) return "";
    const date = new Date(value);
    return Number.isNaN(date.getTime())
      ? String(value)
      : date.toLocaleTimeString([], {hour12: false});
  }

  function formatNumber(value, digits) {
    const number = Number(value);
    return Number.isFinite(number) ? number.toFixed(digits) : String(value ?? "—");
  }

  function firstPresent(...values) {
    return values.find((value) => value !== null && value !== undefined && value !== "");
  }

  function present(value) {
    return value !== null && value !== undefined && value !== "" ? String(value) : "not available";
  }

  function valueOrDash(value) {
    return value !== null && value !== undefined && value !== "" ? String(value) : "—";
  }

  function setText(id, value) {
    const element = byId(id);
    if (element) element.textContent = value == null ? "" : String(value);
  }

  function setClassName(id, value) {
    const element = byId(id);
    if (element) element.className = value;
  }

  function setWidth(id, value) {
    const element = byId(id);
    if (element) element.style.width = value == null ? "0%" : `${Math.max(0, Math.min(100, value))}%`;
  }

  function humanSize(bytes) {
    if (!Number.isFinite(bytes)) return "—";
    if (bytes < 1024) return `${bytes} B`;
    if (bytes < 1024 ** 2) return `${(bytes / 1024).toFixed(1)} KB`;
    if (bytes < 1024 ** 3) return `${(bytes / 1024 ** 2).toFixed(2)} MB`;
    return `${(bytes / 1024 ** 3).toFixed(2)} GB`;
  }

  function humanise(value) {
    return String(value || "")
      .replace(/[_-]+/g, " ")
      .replace(/\b\w/g, (character) => character.toUpperCase());
  }

  async function copyText(text) {
    if (!text) return false;
    try {
      if (navigator.clipboard?.writeText) {
        await navigator.clipboard.writeText(text);
        return true;
      }
      const textarea = document.createElement("textarea");
      textarea.value = text;
      textarea.style.position = "fixed";
      textarea.style.opacity = "0";
      document.body.appendChild(textarea);
      textarea.select();
      const ok = document.execCommand("copy");
      textarea.remove();
      return ok;
    } catch (error) {
      return false;
    }
  }

  function showToast(message, kind) {
    let toast = byId("dashboard-toast");
    if (!toast) {
      toast = document.createElement("div");
      toast.id = "dashboard-toast";
      toast.className = "dashboard-toast";
      toast.setAttribute("role", "status");
      document.body.appendChild(toast);
    }
    toast.textContent = message;
    toast.setAttribute("data-kind", kind || "ok");
    toast.classList.add("show");
    window.clearTimeout(showToast.timer);
    showToast.timer = window.setTimeout(() => toast.classList.remove("show"), 3500);
  }

  function escapeHTML(value) {
    return String(value == null ? "" : value).replace(/[&<>"']/g, (character) => ({
      "&": "&amp;", "<": "&lt;", ">": "&gt;", '"': "&quot;", "'": "&#39;",
    })[character]);
  }

  function escapeAttr(value) {
    return escapeHTML(value);
  }

  function clampInt(value, low, high, fallback) {
    const number = parseInt(value, 10);
    if (!Number.isFinite(number)) return fallback;
    return Math.max(low, Math.min(high, number));
  }

  function clampFloat(value, low, high, fallback) {
    const number = Number(value);
    if (!Number.isFinite(number)) return fallback;
    return Math.max(low, Math.min(high, number));
  }

  function emit(name, detail) {
    window.dispatchEvent(new CustomEvent(`dashboard:${name}`, {detail}));
  }

  window.FastMDXDashboard = {
    get state() { return JSON.parse(JSON.stringify(state)); },
    navigate,
    applyStatus,
    applyAppState,
    applyMetrics,
    applyStructure,
    applyPlayback,
    applyResults,
    on(eventName, handler) {
      window.addEventListener(`dashboard:${eventName}`, (event) => handler(event.detail));
    },
  };
}());
