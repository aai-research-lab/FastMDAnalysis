/* Building a run.
 *
 * There was a page for starting from a protein and a page for starting from a
 * trajectory, and they were the same thing: both write a config and run it.
 * They differed only in which phases the config named -- and having been built
 * separately, one of them offered eleven of the eighty-three settings that
 * exist while the other offered all of them.
 *
 * So there is one page, and it asks three questions in the order somebody
 * actually answers them: what have you got, what should happen to it, and is
 * there anything you want to change. The first answer sets the second, the
 * second decides which settings are worth showing, and nothing is drawn from a
 * list kept here -- it all comes from the schema, which is read from the
 * software itself.
 */

(function () {
  "use strict";

  const PHASES = [
    {
      name: "setup",
      label: "Setup",
      blurb: "Protonate, solvate, add ions, and assign a force field.",
    },
    {
      name: "simulation",
      label: "Simulate",
      blurb: "Minimise, equilibrate, and run production dynamics.",
    },
    {
      name: "analysis",
      label: "Analyze",
      blurb: "Measure the trajectory: structure, flexibility, interactions.",
    },
    {
      name: "report",
      label: "Report",
      blurb: "Collect the figures and numbers into a document.",
    },
  ];

  /* What each starting point implies. Somebody with a protein wants the whole
   * thing; somebody with a trajectory has already done the expensive part. */
  const STARTING_POINTS = {
    structure: {
      label: "A structure",
      detail:
        "A four-character PDB identifier, or a path to a PDB or CIF file. " +
        "Everything is built from it.",
      phases: ["setup", "simulation", "analysis", "report"],
      // Sequence-to-structure is not implemented -- the setup phase refuses
      // it -- so it is not offered here.
      offers: ["setup", "simulation", "analysis", "report"],
    },
    trajectory: {
      label: "A trajectory",
      detail:
        "A simulation you already have, from here or anywhere else. It is " +
        "measured as it stands.",
      phases: ["analysis", "report"],
      // Setup and simulation have nothing to act on: there is no supported
      // way to continue a run from a trajectory, and re-preparing the
      // structure would not connect to the frames already recorded.
      offers: ["analysis", "report"],
    },
    config: {
      label: "A config I already have",
      detail:
        "One written earlier, edited by hand, or brought back from a " +
        "cluster. It is checked before anything runs.",
      phases: [],
      offers: [],
    },
  };

  const state = {
    schema: null,
    start: null,
    phases: new Set(),
    values: {},            // phase -> { field: value }
    analyses: new Set(),
    analysisOptions: {},   // analysis -> { option: value }
    open: new Set(),
    found: null,
    configVerdict: null,
    loadedFrom: null,
  };

  const el = (id) => document.getElementById(id);
  const text = (node, value) => { if (node) node.textContent = value; };

  // ------------------------------------------------------------------ schema

  async function loadSchema() {
    if (state.schema) return state.schema;
    const response = await fetch("/api/schema");
    state.schema = await response.json();
    return state.schema;
  }

  function fieldsFor(phase) {
    const block = state.schema.phases[phase];
    return block ? block.fields : [];
  }

  /* Settings a run supplies for itself. Asking somebody to type the trajectory
   * path in the settings when they chose it two questions ago is the kind of
   * thing that makes a form feel like paperwork. */
  const SUPPLIED = new Set(["trajectory", "topology"]);

  function settingsFor(phase) {
    return fieldsFor(phase).filter((field) => !SUPPLIED.has(field.name));
  }

  // ------------------------------------------------------- what have you got

  function renderStart() {
    const select = el("run-start");
    if (!select) return;
    if (!select.options.length) {
      const blank = document.createElement("option");
      blank.value = "";
      blank.textContent = "Choose…";
      select.appendChild(blank);
      Object.keys(STARTING_POINTS).forEach((key) => {
        const option = document.createElement("option");
        option.value = key;
        option.textContent = STARTING_POINTS[key].label;
        select.appendChild(option);
      });
      select.addEventListener("change", () => {
        state.start = select.value || null;
        const choice = STARTING_POINTS[state.start] || { phases: [] };
        state.phases = new Set(choice.phases);
        renderAll();
      });
    }
    select.value = state.start || "";

    // The description belongs with the control, not on a card of its own:
    // it is what the option means, and it changes as the option does.
    text(
      el("run-start-detail"),
      state.start ? STARTING_POINTS[state.start].detail : ""
    );

    const started = Boolean(state.start);
    // Once a config has been opened into the form, the form is what is being
    // edited -- the starting point is whatever the config described.
    const fromConfig = state.start === "config";
    el("run-input-card").hidden = !started;
    el("run-phases-card").hidden = !started || fromConfig;
    el("run-settings-card").hidden = !started || fromConfig;
    el("run-actions-card").hidden = !started;

    el("run-structure-field").hidden = state.start !== "structure";
    el("run-trajectory-field").hidden = state.start !== "trajectory";
    el("run-config-field").hidden = !fromConfig;
    el("run-output-field").hidden = fromConfig;
  }

  // ------------------------------------------------ what should happen to it

  function renderPhases() {
    const host = el("run-phases");
    if (!host) return;
    host.innerHTML = "";
    PHASES.forEach((phase) => {
      const row = document.createElement("label");
      row.className = "run-phase";
      row.dataset.chosen = String(state.phases.has(phase.name));

      const offered = STARTING_POINTS[state.start]
        ? STARTING_POINTS[state.start].offers
        : [];
      const available = offered.includes(phase.name);
      row.dataset.available = String(available);

      const box = document.createElement("input");
      box.type = "checkbox";
      box.checked = available && state.phases.has(phase.name);
      box.disabled = !available;
      box.addEventListener("change", () => {
        if (box.checked) state.phases.add(phase.name);
        else state.phases.delete(phase.name);
        renderAll();
      });

      const naming = document.createElement("span");
      naming.className = "run-phase-naming";
      const label = document.createElement("span");
      label.className = "run-phase-label";
      label.textContent = phase.label;
      const blurb = document.createElement("span");
      blurb.className = "run-phase-blurb";
      blurb.textContent = available
        ? phase.blurb
        : "Nothing to act on: a trajectory is already the result of this.";
      naming.appendChild(label);
      naming.appendChild(blurb);

      row.appendChild(box);
      row.appendChild(naming);
      host.appendChild(row);
    });
  }

  // ------------------------------------------- anything you want to change

  function renderSettings() {
    const host = el("run-settings");
    if (!host) return;
    host.innerHTML = "";

    const chosen = PHASES.filter((phase) => state.phases.has(phase.name));
    if (!chosen.length) {
      host.innerHTML =
        '<div class="empty-detail muted">Choose at least one thing to do.</div>';
      return;
    }

    chosen.forEach((phase) => {
      host.appendChild(settingsSection(phase));
      if (phase.name === "analysis") host.appendChild(analysisSection());
    });
  }

  function settingsSection(phase) {
    const section = document.createElement("div");
    section.className = "run-section";

    const settings = settingsFor(phase.name);
    const changed = Object.keys(state.values[phase.name] || {}).length;

    const head = document.createElement("button");
    head.type = "button";
    head.className = "run-section-head";
    head.setAttribute("aria-expanded", String(state.open.has(phase.name)));
    head.innerHTML =
      `<span class="run-section-name">${phase.label}</span>` +
      `<span class="run-section-count">${
        changed ? `${changed} changed` : `${settings.length} settings`
      }</span>`;
    head.addEventListener("click", () => {
      if (state.open.has(phase.name)) state.open.delete(phase.name);
      else state.open.add(phase.name);
      renderSettings();
    });
    section.appendChild(head);

    if (state.open.has(phase.name)) {
      const grid = document.createElement("div");
      grid.className = "run-section-body";
      settings.forEach((field) => grid.appendChild(control(phase.name, field)));
      section.appendChild(grid);
    }
    return section;
  }

  function control(phase, field) {
    const wrap = document.createElement("label");
    wrap.className = "builder-field";

    const label = document.createElement("span");
    label.className = "builder-label";
    label.textContent = field.name.replace(/_/g, " ");
    wrap.appendChild(label);

    let input;
    if (field.choices && field.choices.length) {
      input = document.createElement("select");
      field.choices.forEach((choice) => {
        const option = document.createElement("option");
        option.value = choice;
        option.textContent = choice;
        if (choice === field.default) option.selected = true;
        input.appendChild(option);
      });
    } else if (field.control === "checkbox") {
      input = document.createElement("input");
      input.type = "checkbox";
      input.checked = Boolean(field.default);
    } else {
      input = document.createElement("input");
      input.type = field.control === "number" ? "number" : "text";
      if (field.control === "number") input.step = "any";
      // Shown rather than filled in, so the config records what was decided.
      input.placeholder =
        field.default === null || field.default === undefined
          ? (field.example === null ? "" : String(field.example ?? ""))
          : String(field.default);
    }
    const current = (state.values[phase] || {})[field.name];
    if (current !== undefined) {
      if (input.type === "checkbox") input.checked = Boolean(current);
      else input.value = current;
    }

    input.addEventListener("change", () => {
      const value = input.type === "checkbox" ? input.checked : input.value;
      state.values[phase] = state.values[phase] || {};
      if (value === "" || value === null) delete state.values[phase][field.name];
      else state.values[phase][field.name] = value;
      renderSettings();
      updateSummary();
    });
    wrap.appendChild(input);

    if (field.help) {
      const help = document.createElement("span");
      help.className = "builder-card-note";
      help.textContent = field.help;
      wrap.appendChild(help);
    }
    return wrap;
  }

  // --------------------------------------------------- which measurements

  function analysisSection() {
    const section = document.createElement("div");
    section.className = "run-section run-section-open";

    const head = document.createElement("div");
    head.className = "run-section-head run-section-head-static";
    head.innerHTML =
      '<span class="run-section-name">Measurements</span>' +
      `<span class="run-section-count">${
        state.analyses.size ? `${state.analyses.size} chosen` : "all of them"
      }</span>`;
    section.appendChild(head);

    const note = document.createElement("p");
    note.className = "builder-card-note";
    note.textContent =
      "Choose nothing and every analysis that applies is run.";
    section.appendChild(note);

    const options = state.schema.analysis_options;
    const body = document.createElement("div");
    body.className = "run-analyses";
    if (!options || !options.available) {
      body.innerHTML = `<div class="empty-detail muted">${
        (options && options.reason) || "No analyses available."
      }</div>`;
      section.appendChild(body);
      return section;
    }

    Object.keys(options.analyses).sort().forEach((name) => {
      body.appendChild(analysisRow(name, options));
    });
    section.appendChild(body);
    return section;
  }

  function analysisRow(name, options) {
    const explanation = (options.explanations || {})[name] || {};
    const settings = (options.analyses[name] || []).filter((o) => !o.shared);

    const row = document.createElement("div");
    row.className = "analyse-row";
    row.dataset.chosen = String(state.analyses.has(name));

    const head = document.createElement("div");
    head.className = "analyse-row-head";

    const label = document.createElement("label");
    label.className = "analyse-choice";
    const box = document.createElement("input");
    box.type = "checkbox";
    box.checked = state.analyses.has(name);
    box.addEventListener("change", () => {
      if (box.checked) state.analyses.add(name);
      else state.analyses.delete(name);
      renderSettings();
      updateSummary();
    });
    const naming = document.createElement("span");
    naming.className = "analyse-naming";
    const title = document.createElement("span");
    title.className = "analyse-name";
    title.textContent = name;
    naming.appendChild(title);
    if (explanation.title) {
      const what = document.createElement("span");
      what.className = "analyse-what";
      what.textContent = explanation.title;
      naming.appendChild(what);
    }
    label.appendChild(box);
    label.appendChild(naming);

    const toggle = document.createElement("button");
    toggle.type = "button";
    toggle.className = "analyse-toggle";
    toggle.textContent = settings.length
      ? `${settings.length} setting${settings.length === 1 ? "" : "s"}`
      : "no settings";
    toggle.disabled = !settings.length;
    const key = `analysis:${name}`;
    toggle.setAttribute("aria-expanded", String(state.open.has(key)));
    toggle.addEventListener("click", () => {
      if (state.open.has(key)) state.open.delete(key);
      else state.open.add(key);
      renderSettings();
    });

    head.appendChild(label);
    head.appendChild(toggle);
    row.appendChild(head);

    if (explanation.detail) {
      const explains = document.createElement("p");
      explains.className = "analyse-explains";
      explains.textContent = explanation.detail;
      row.appendChild(explains);
    }

    if (state.open.has(key)) {
      const body = document.createElement("div");
      body.className = "analyse-settings";
      settings.forEach((option) => body.appendChild(analysisControl(name, option)));
      row.appendChild(body);
    }
    return row;
  }

  function analysisControl(analysis, option) {
    const wrap = document.createElement("label");
    wrap.className = "builder-field";
    const label = document.createElement("span");
    label.className = "builder-label";
    label.textContent = option.name.replace(/_/g, " ");
    wrap.appendChild(label);

    let input;
    if (option.choices && option.choices.length) {
      input = document.createElement("select");
      option.choices.forEach((choice) => {
        const item = document.createElement("option");
        item.value = choice;
        item.textContent = choice;
        if (choice === option.default) item.selected = true;
        input.appendChild(item);
      });
    } else if (typeof option.default === "boolean") {
      input = document.createElement("input");
      input.type = "checkbox";
      input.checked = option.default;
    } else {
      input = document.createElement("input");
      input.type = typeof option.default === "number" ? "number" : "text";
      if (input.type === "number") input.step = "any";
      input.placeholder =
        option.default === null || option.default === undefined
          ? ""
          : String(option.default);
    }
    const current = (state.analysisOptions[analysis] || {})[option.name];
    if (current !== undefined) {
      if (input.type === "checkbox") input.checked = Boolean(current);
      else input.value = current;
    }
    input.addEventListener("change", () => {
      const value = input.type === "checkbox" ? input.checked : input.value;
      state.analysisOptions[analysis] = state.analysisOptions[analysis] || {};
      if (value === "" || value === null) {
        delete state.analysisOptions[analysis][option.name];
      } else {
        state.analysisOptions[analysis][option.name] = value;
      }
      updateSummary();
    });
    wrap.appendChild(input);

    if (option.help) {
      const help = document.createElement("span");
      help.className = "builder-card-note";
      help.textContent = option.help;
      wrap.appendChild(help);
    }
    return wrap;
  }

  // ------------------------------------------------------------- the config

  function currentState() {
    const config = {
      output: el("run-output").value.trim() || "fastmdxplora_output",
      include: PHASES.filter((p) => state.phases.has(p.name)).map((p) => p.name),
    };

    if (state.start === "structure") {
      config.system = el("run-system").value.trim();
    } else {
      // The structure a trajectory refers to is the system, and the
      // trajectory itself is where the analysis looks.
      config.system = el("run-topology").value.trim();
      config.analysis = {
        trajectory: el("run-trajectory").value.trim(),
        topology: el("run-topology").value.trim(),
      };
    }

    PHASES.forEach((phase) => {
      if (!state.phases.has(phase.name)) return;
      const chosen = state.values[phase.name];
      if (chosen && Object.keys(chosen).length) {
        config[phase.name] = Object.assign(config[phase.name] || {}, chosen);
      }
    });

    if (state.phases.has("analysis")) {
      const analysis = config.analysis || {};
      if (state.analyses.size) {
        analysis.include = Array.from(state.analyses).join(", ");
      }
      const nested = {};
      Object.keys(state.analysisOptions).forEach((name) => {
        const values = state.analysisOptions[name];
        const wanted = !state.analyses.size || state.analyses.has(name);
        if (wanted && values && Object.keys(values).length) nested[name] = values;
      });
      if (Object.keys(nested).length) analysis.options = nested;
      config.analysis = Object.assign(analysis, config.analysis || {});
    }

    const everything = el("run-full-config");
    config.full = Boolean(everything && everything.checked);
    return config;
  }

  function ready() {
    if (!state.start) return false;
    if (state.start === "config") {
      // Nothing is offered until the file has been read and found sound.
      return Boolean(state.configVerdict && state.configVerdict.ok);
    }
    if (!state.phases.size) return false;
    if (state.start === "structure") return Boolean(el("run-system").value.trim());
    return Boolean(el("run-trajectory").value.trim() && el("run-topology").value.trim());
  }

  function updateSummary() {
    if (state.start === "config") {
      const verdict = state.configVerdict;
      text(
        el("run-summary"),
        verdict && verdict.ok ? verdict.phases.join(" → ") : "No config checked yet"
      );
    } else {
      const doing = PHASES.filter((p) => state.phases.has(p.name))
        .map((p) => p.label);
      text(
        el("run-summary"),
        doing.length ? doing.join(" → ") : "Nothing chosen yet"
      );
    }
    const can = ready();
    ["run-start-button", "run-download"].forEach((id) => {
      const button = el(id);
      if (button) button.disabled = !can;
    });
  }

  async function fetchConfig() {
    const response = await fetch("/api/config", {
      method: "POST",
      headers: { "Content-Type": "application/json" },
      body: JSON.stringify(currentState()),
    });
    return response.json();
  }

  async function showConfig() {
    const box = el("run-config-preview");
    const button = el("run-preview");
    if (!box) return;
    if (!box.hidden) {
      box.hidden = true;
      text(button, "Show the config");
      button.setAttribute("aria-expanded", "false");
      return;
    }
    const built = await fetchConfig();
    box.textContent = built.ok ? built.yaml : built.error;
    box.hidden = false;
    text(button, "Hide the config");
    button.setAttribute("aria-expanded", "true");
  }

  async function download() {
    const built = await fetchConfig();
    if (!built.ok) {
      text(el("run-note"), built.error);
      return;
    }
    text(el("run-note"), `${built.settings_changed} setting(s) written.`);
    const blob = new Blob([built.yaml], { type: "text/yaml" });
    const link = document.createElement("a");
    link.href = URL.createObjectURL(blob);
    link.download = "fastmdxplora.yml";
    link.click();
    URL.revokeObjectURL(link.href);
  }

  /* A config that has been elsewhere can be run as it stands, or opened and
   * changed. Opening it never writes to it: the file may be committed beside
   * a paper or be the record of a run that already happened, so anything
   * altered is saved as a new file and the original is left alone. */
  async function loadConfigIntoForm() {
    const path = el("run-config-path").value.trim();
    if (!path) return;
    const response = await fetch("/api/load-config", {
      method: "POST",
      headers: { "Content-Type": "application/json" },
      body: JSON.stringify({ path }),
    });
    const loaded = await response.json();
    if (!loaded.ok) {
      text(el("run-config-verdict"), loaded.error);
      return;
    }

    const from = loaded.state;
    state.start = from.start;
    state.phases = new Set(from.include);
    state.values = from.phases || {};
    state.analyses = new Set(from.analyses || []);
    state.analysisOptions = from.analysis_options || {};
    state.loadedFrom = path;

    renderAll();
    // The fields the form owns rather than the settings tables.
    if (el("run-system")) el("run-system").value = from.system || "";
    if (el("run-trajectory")) el("run-trajectory").value = from.trajectory || "";
    if (el("run-topology")) el("run-topology").value = from.topology || "";
    if (el("run-output")) el("run-output").value = from.output || "";
    updateSummary();
    text(
      el("run-note"),
      `Opened ${path}. Changes here are saved as a new file; that one is ` +
      "left as it is."
    );
  }

  async function runAsItStands() {
    const path = el("run-config-path").value.trim();
    if (!path) return;
    const button = el("run-as-is");
    button.disabled = true;
    text(el("run-note"), "Starting\u2026");
    const response = await fetch("/api/run-config", {
      method: "POST",
      headers: { "Content-Type": "application/json" },
      body: JSON.stringify({ path, output: el("run-output").value.trim() }),
    });
    const started = await response.json();
    if (!started.ok) {
      text(el("run-note"), started.error || "Could not start.");
      button.disabled = false;
      return;
    }
    text(el("run-note"), `Running ${started.config_path} as it stands.`);
    if (window.FastMDXDashboard && window.FastMDXDashboard.navigate) {
      window.FastMDXDashboard.navigate("overview");
    }
  }

  async function checkConfigFile() {
    const path = el("run-config-path").value.trim();
    const note = el("run-config-verdict");
    if (!path) {
      text(note, "");
      return null;
    }
    const response = await fetch("/api/check-config", {
      method: "POST",
      headers: { "Content-Type": "application/json" },
      body: JSON.stringify({ path }),
    });
    const verdict = await response.json();
    text(
      note,
      verdict.ok
        ? `Runs. ${verdict.phases.join(" → ")}, ` +
          `${verdict.systems} system(s), ${verdict.settings_named} setting(s) named.`
        : verdict.error
    );
    note.dataset.ok = String(Boolean(verdict.ok));
    state.configVerdict = verdict;
    // Opening and running are only offered once the file is known to be sound.
    ["run-open-config", "run-as-is"].forEach((id) => {
      const button = el(id);
      if (button) button.disabled = !verdict.ok;
    });
    updateSummary();
    return verdict;
  }

  /* Stopping a run. This lived only on the page being retired, so deleting
   * that page without bringing it across would have left the browser able to
   * start something it could not stop. */
  async function stopRunning() {
    const button = el("run-stop");
    button.disabled = true;
    text(el("run-note"), "Stopping\u2026");
    try {
      const response = await fetch("/api/explore/stop", {
        method: "POST",
        headers: { "Content-Type": "application/json" },
        body: JSON.stringify({}),
      });
      const stopped = await response.json();
      text(el("run-note"), stopped.ok ? "Stopped." : (stopped.error || "Could not stop."));
    } catch (error) {
      text(el("run-note"), "Could not reach the server.");
    }
    button.disabled = false;
  }

  /* The button appears only while something is running, and the dashboard is
   * what knows. */
  function watchForARun() {
    if (!window.FastMDXDashboard || !window.FastMDXDashboard.on) return;
    // The handler is given the detail itself, not the event carrying it.
    window.FastMDXDashboard.on("app-state", (detail) => {
      const running = Boolean(detail && detail.active_run);
      const stop = el("run-stop");
      if (stop) stop.hidden = !running;
    });
  }

  async function start() {
    const button = el("run-start-button");
    button.disabled = true;
    text(el("run-note"), "Starting\u2026");
    let started;
    try {
      const response = await fetch("/api/run", {
        method: "POST",
        headers: { "Content-Type": "application/json" },
        body: JSON.stringify(currentState()),
      });
      started = await response.json();
    } catch (error) {
      text(el("run-note"), "Could not reach the server.");
      button.disabled = false;
      return;
    }
    if (!started.ok) {
      text(el("run-note"), started.error || "Could not start.");
      button.disabled = false;
      return;
    }
    text(el("run-note"), `Running. Config written to ${started.config_path}`);
    if (window.FastMDXDashboard && window.FastMDXDashboard.navigate) {
      window.FastMDXDashboard.navigate("overview");
    }
  }

  function resetEverything() {
    state.phases = state.start
      ? new Set(STARTING_POINTS[state.start].phases)
      : new Set();
    state.values = {};
    state.analyses.clear();
    state.analysisOptions = {};
    state.open.clear();
    const output = el("run-output");
    if (output) output.value = "";
    const everything = el("run-full-config");
    if (everything) everything.checked = false;
    const box = el("run-config-preview");
    if (box) box.hidden = true;
    text(el("run-preview"), "Show the config");
    text(el("run-note"), "Everything is back to its default.");
    renderAll();
  }

  // ------------------------------------------------------------------- wire

  function renderAll() {
    renderStart();
    renderPhases();
    renderSettings();
    updateSummary();
  }


  /* Where the results will actually be written, for whatever is in the box.
   * Saying "a name lands in /somewhere" left the reader to do the joining;
   * this says the path. */
  function describeOutput() {
    const box = el("run-output");
    const where = (state.schema && state.schema.workspace) || "";
    if (!box) return;
    const typed = box.value.trim();
    const name = typed || "fastmdxplora_output";
    const absolute = name.startsWith("/") || name.startsWith("~");
    const full = absolute ? name : (where ? `${where}/${name}` : name);
    text(el("run-output-note"), `Results will be saved in ${full}`);
  }

  /* A trajectory usually sits beside the structure it moves. Choosing one and
   * then being asked to find the other in the same folder is a step the page
   * can take itself, and the button that used to do it said "Find beside it",
   * which explained nothing. */
  async function findStructureBeside() {
    const trajectory = el("run-trajectory").value.trim();
    const topology = el("run-topology");
    if (!trajectory || !topology || topology.value.trim()) return;

    const folder = trajectory.replace(/[/\\][^/\\]*$/, "");
    if (!folder) return;
    try {
      const response = await fetch(
        "/api/inspect-directory?path=" + encodeURIComponent(folder)
      );
      const found = await response.json();
      if (found.ok && found.suggestion && found.suggestion.topology) {
        topology.value = found.suggestion.topology;
        updateSummary();
      }
    } catch (error) {
      // Not finding one is not a failure: the field is still there to type in.
    }
  }

  function attach() {
    if (!el("run-start")) return;

    ["run-system", "run-trajectory", "run-topology", "run-output"].forEach((id) => {
      const input = el(id);
      if (input) input.addEventListener("input", updateSummary);
    });
    const output = el("run-output");
    if (output) output.addEventListener("input", describeOutput);
    const trajectory = el("run-trajectory");
    if (trajectory) trajectory.addEventListener("change", findStructureBeside);

    const configPath = el("run-config-path");
    if (configPath) configPath.addEventListener("change", checkConfigFile);
    const checkButton = el("run-check-config");
    if (checkButton) checkButton.addEventListener("click", checkConfigFile);
    const openButton = el("run-open-config");
    if (openButton) openButton.addEventListener("click", loadConfigIntoForm);
    const asIsButton = el("run-as-is");
    if (asIsButton) asIsButton.addEventListener("click", runAsItStands);

    const preview = el("run-preview");
    if (preview) preview.addEventListener("click", showConfig);
    const downloadButton = el("run-download");
    if (downloadButton) downloadButton.addEventListener("click", download);
    const startButton = el("run-start-button");
    if (startButton) startButton.addEventListener("click", start);
    const reset = el("run-reset");
    if (reset) reset.addEventListener("click", resetEverything);
    const stop = el("run-stop");
    if (stop) stop.addEventListener("click", stopRunning);
    watchForARun();

    loadSchema().then(() => {
      describeOutput();
      renderAll();
    }).catch(() => {
      text(el("run-summary"), "Could not read the settings.");
    });
  }

  if (document.readyState === "loading") {
    document.addEventListener("DOMContentLoaded", attach);
  } else {
    attach();
  }

  window.FastMDXRun = { state, currentState, PHASES, STARTING_POINTS };
})();
