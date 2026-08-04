/* Analysing data you already have.
 *
 * The other way in. "New Exploration" is for someone starting from a protein;
 * this is for someone who already has a trajectory -- from GROMACS, from
 * AMBER, from a run last week -- and wants it analysed without reproducing the
 * simulation first. That is most people.
 *
 * Nothing here holds a list of what can be set. The analyses and their
 * settings come from /api/schema, which reads them from the analyses
 * themselves, so an analysis gaining an option puts a control on this page and
 * there is no second list to keep in step. The form that used to be written by
 * hand offered eleven of eighty-three settings, which is what happens to a
 * list somebody has to remember to update.
 */

(function () {
  "use strict";

  const state = {
    schema: null,
    found: null,
    chosen: new Set(),
    options: {},          // analysis name -> { option name: value }
    expanded: new Set(),
  };

  const el = (id) => document.getElementById(id);

  function text(node, value) {
    if (node) node.textContent = value;
  }

  // ---------------------------------------------------------------- helpers

  /* Nothing is offered before it is known to exist, so the page is asked to
   * draw only after the schema has arrived. */
  async function loadSchema() {
    if (state.schema) return state.schema;
    const response = await fetch("/api/schema");
    state.schema = await response.json();
    return state.schema;
  }

  function analysisNames() {
    const options = state.schema && state.schema.analysis_options;
    if (!options || !options.available) return [];
    return Object.keys(options.analyses).sort();
  }

  function optionsFor(name) {
    const all = state.schema.analysis_options.analyses[name] || [];
    // Settings every analysis shares -- figure size, titles -- are real but
    // they are not what distinguishes one analysis from another, so they do
    // not compete for attention with the ones that do.
    return all.filter((option) => !option.shared);
  }

  // ------------------------------------------------------------- the folder

  async function inspect() {
    const path = el("analyse-folder").value.trim();
    const summary = el("analyse-found");
    if (!path) {
      summary.hidden = true;
      return;
    }
    text(el("analyse-folder-note"), "Looking…");
    let found;
    try {
      const response = await fetch(
        "/api/inspect-directory?path=" + encodeURIComponent(path)
      );
      found = await response.json();
    } catch (error) {
      text(el("analyse-folder-note"), "Could not read that folder.");
      return;
    }

    state.found = found;
    if (!found.ok) {
      text(el("analyse-folder-note"), found.error);
      summary.hidden = true;
      return;
    }

    if (!found.can_analyse) {
      text(
        el("analyse-folder-note"),
        found.can_set_up
          ? "A structure, but no trajectory. Start from New Exploration to simulate it."
          : "Nothing to analyse here."
      );
      summary.hidden = true;
      return;
    }

    text(
      el("analyse-folder-note"),
      found.is_previous_run
        ? "A previous run of this software."
        : `${found.trajectories.length} trajectory file(s) found.`
    );

    fillFileChoice("analyse-trajectory", found.trajectories,
                   found.suggestion && found.suggestion.trajectory,
                   found.directory);
    fillFileChoice("analyse-topology", found.topologies,
                   found.suggestion && found.suggestion.topology,
                   found.directory);
    text(
      el("analyse-trajectory-note"),
      found.trajectories.length > 1
        ? `${found.trajectories.length} found. The largest is selected; change it if the production run is not the longest.`
        : "The trajectory to measure."
    );
    summary.hidden = false;
    render();
  }

  function fillFileChoice(id, entries, preferred, base) {
    const select = el(id);
    if (!select) return;
    select.innerHTML = "";
    entries.forEach((entry) => {
      const option = document.createElement("option");
      option.value = entry.path;
      // Named by where it sits, not just what it is called: a run with
      // replicates has a production.dcd in each of them, and a list of
      // identical names is no choice at all. The size is how somebody tells
      // the production run from the equilibration without opening either.
      let label = entry.path;
      if (base && label.startsWith(base)) {
        label = label.slice(base.length).replace(/^[\/\\]/, "");
      }
      option.textContent = `${label} — ${entry.size_mb} MB`;
      if (entry.path === preferred) option.selected = true;
      select.appendChild(option);
    });
  }

  // ------------------------------------------------------------ the choices

  function render() {
    const host = el("analyse-choices");
    if (!host) return;
    host.innerHTML = "";

    const names = analysisNames();
    if (!names.length) {
      const reason =
        (state.schema.analysis_options &&
          state.schema.analysis_options.reason) ||
        "No analyses are available.";
      host.innerHTML = `<div class="empty-detail muted">${reason}</div>`;
      return;
    }

    names.forEach((name) => host.appendChild(analysisRow(name)));
    updateSummary();
  }

  function analysisRow(name) {
    const row = document.createElement("div");
    row.className = "analyse-row";

    const head = document.createElement("div");
    head.className = "analyse-row-head";

    const label = document.createElement("label");
    label.className = "analyse-choice";
    const box = document.createElement("input");
    box.type = "checkbox";
    box.checked = state.chosen.has(name);
    box.addEventListener("change", () => {
      if (box.checked) state.chosen.add(name);
      else state.chosen.delete(name);
      row.dataset.chosen = box.checked ? "true" : "false";
      updateSummary();
    });
    const explanation =
      (state.schema.analysis_options.explanations || {})[name] || {};

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

    const settings = optionsFor(name);
    const toggle = document.createElement("button");
    toggle.type = "button";
    toggle.className = "analyse-toggle";
    toggle.textContent = settings.length
      ? `${settings.length} setting${settings.length === 1 ? "" : "s"}`
      : "no settings";
    toggle.disabled = !settings.length;
    toggle.addEventListener("click", () => {
      if (state.expanded.has(name)) state.expanded.delete(name);
      else state.expanded.add(name);
      body.hidden = !state.expanded.has(name);
      toggle.setAttribute("aria-expanded", String(!body.hidden));
    });
    toggle.setAttribute("aria-expanded", String(state.expanded.has(name)));

    head.appendChild(label);
    head.appendChild(toggle);
    row.appendChild(head);

    /* What the measurement means and how to read it -- already written at the
     * top of every analysis module, and until now only visible to somebody
     * reading the source. It belongs next to the decision it informs. */
    if (explanation.detail) {
      const note = document.createElement("p");
      note.className = "analyse-explains";
      note.textContent = explanation.detail;
      row.appendChild(note);
    }

    const body = document.createElement("div");
    body.className = "analyse-settings";
    body.hidden = !state.expanded.has(name);
    settings.forEach((option) => body.appendChild(optionField(name, option)));
    row.appendChild(body);

    row.dataset.chosen = box.checked ? "true" : "false";
    return row;
  }

  function optionField(analysis, option) {
    const field = document.createElement("label");
    field.className = "builder-field";

    const label = document.createElement("span");
    label.className = "builder-label";
    label.textContent = option.name;
    field.appendChild(label);

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
    } else if (option.control === "checkbox") {
      input = document.createElement("input");
      input.type = "checkbox";
      input.checked = Boolean(option.default);
    } else {
      input = document.createElement("input");
      input.type = option.control === "number" ? "number" : "text";
      input.step = "any";
      // The default is shown rather than filled in, so what is written here
      // is what somebody decided and the file records only that.
      input.placeholder =
        option.default === null || option.default === undefined
          ? ""
          : String(option.default);
    }

    input.addEventListener("change", () => {
      const value = input.type === "checkbox" ? input.checked : input.value;
      state.options[analysis] = state.options[analysis] || {};
      if (value === "" || value === null) delete state.options[analysis][option.name];
      else state.options[analysis][option.name] = value;
      updateSummary();
    });
    field.appendChild(input);

    if (option.help) {
      const help = document.createElement("span");
      help.className = "builder-card-note";
      help.textContent = option.help;
      field.appendChild(help);
    }
    return field;
  }



  // -------------------------------------------------------------- the reset

  /* Putting everything back is a real need, not a nicety: somebody exploring
   * what the settings do will change several, lose track of which, and have
   * no way to tell a considered choice from a leftover. The config records
   * only what differs from the default, so a forgotten change is invisible in
   * the file and shows up in the results. */
  function resetEverything() {
    state.chosen.clear();
    state.options = {};
    state.expanded.clear();

    const scope = el("analyse-scope");
    if (scope) {
      const fromSchema = (state.schema.phases.analysis.fields || [])
        .find((field) => field.name === "scope");
      scope.value = (fromSchema && fromSchema.default) || "solute";
    }
    const output = el("analyse-output");
    if (output) output.value = "";
    const everything = el("analyse-full-config");
    if (everything) everything.checked = false;

    const box = el("analyse-config-preview");
    if (box) box.hidden = true;
    const previewButton = el("analyse-preview");
    if (previewButton) {
      previewButton.textContent = "Show the config";
      previewButton.setAttribute("aria-expanded", "false");
    }
    text(el("analyse-config-note"), "Everything is back to its default.");

    render();
  }

  // ------------------------------------------------------------- the config

  function currentState() {
    const analysis = {
      trajectory: el("analyse-trajectory").value,
      topology: el("analyse-topology").value,
    };
    if (state.chosen.size) analysis.include = Array.from(state.chosen).join(", ");
    const scope = el("analyse-scope");
    if (scope && scope.value) analysis.scope = scope.value;

    const nested = {};
    Object.keys(state.options).forEach((name) => {
      if (state.chosen.has(name) && Object.keys(state.options[name]).length) {
        nested[name] = state.options[name];
      }
    });
    if (Object.keys(nested).length) analysis.options = nested;

    return {
      output: el("analyse-output").value.trim() || "analysis_output",
      include: ["analysis"],
      systems: [{ system: analysis.topology || "", id: "analysed" }],
      analysis,
    };
  }

  function updateSummary() {
    const chosen = state.chosen.size;
    text(
      el("analyse-summary"),
      chosen
        ? `${chosen} analysis${chosen === 1 ? "" : "es"} selected`
        : "Nothing selected yet"
    );
    const ready = chosen > 0 && state.found && state.found.can_analyse;
    ["analyse-run", "analyse-download"].forEach((id) => {
      const button = el(id);
      if (button) button.disabled = !ready;
    });
  }

  async function fetchConfig() {
    const everything = el("analyse-full-config");
    const body = currentState();
    // Two files answer two questions: what did I decide, and what did the run
    // actually use. Both are read back the same way.
    body.full = Boolean(everything && everything.checked);
    const response = await fetch("/api/config", {
      method: "POST",
      headers: { "Content-Type": "application/json" },
      body: JSON.stringify(body),
    });
    return response.json();
  }

  async function download() {
    const built = await fetchConfig();
    const note = el("analyse-config-note");
    if (!built.ok) {
      // The settings would be refused, and saying so here is the whole point
      // of checking before the file leaves: the form that produced it is
      // still on the screen.
      text(note, built.error);
      return;
    }
    text(note, `${built.settings_changed} setting(s) written.`);
    const blob = new Blob([built.yaml], { type: "text/yaml" });
    const link = document.createElement("a");
    link.href = URL.createObjectURL(blob);
    link.download = "analysis.yml";
    link.click();
    URL.revokeObjectURL(link.href);
  }

  async function run() {
    const note = el("analyse-config-note");
    const button = el("analyse-run");
    button.disabled = true;
    text(note, "Starting\u2026");
    let started;
    try {
      const response = await fetch("/api/run", {
        method: "POST",
        headers: { "Content-Type": "application/json" },
        body: JSON.stringify(currentState()),
      });
      started = await response.json();
    } catch (error) {
      text(note, "Could not reach the server.");
      button.disabled = false;
      return;
    }
    if (!started.ok) {
      text(note, started.error || "Could not start.");
      button.disabled = false;
      return;
    }
    // The config it ran from is written beside the results, so the run can be
    // repeated later without anybody remembering what was clicked.
    text(note, `Running. Config written to ${started.config_path}`);
    // Straight to the overview, where the run reports itself. The dashboard
    // module owns navigation, so it is asked rather than the page reloaded.
    if (window.FastMDXDashboard && window.FastMDXDashboard.navigate) {
      window.FastMDXDashboard.navigate("overview");
    } else {
      location.hash = "#overview";
      location.reload();
    }
  }

  async function preview() {
    const box = el("analyse-config-preview");
    const button = el("analyse-preview");
    if (!box) return;
    // Toggling, because the same button that opened it is where somebody
    // will look to close it again.
    if (!box.hidden) {
      box.hidden = true;
      if (button) {
        button.textContent = "Show the config";
        button.setAttribute("aria-expanded", "false");
      }
      return;
    }
    const built = await fetchConfig();
    box.textContent = built.ok ? built.yaml : built.error;
    box.hidden = false;
    if (button) {
      button.textContent = "Hide the config";
      button.setAttribute("aria-expanded", "true");
    }
  }

  // ------------------------------------------------------------------- wire

  function attach() {
    const folder = el("analyse-folder");
    if (!folder) return;   // this page is not on the document

    folder.addEventListener("change", inspect);
    const look = el("analyse-look");
    if (look) look.addEventListener("click", inspect);

    const everything = el("analyse-full-config");
    if (everything) {
      everything.addEventListener("change", () => {
        // If the config is on screen, it should be the one just asked for.
        const box = el("analyse-config-preview");
        if (box && !box.hidden) {
          box.hidden = true;
          preview();
        }
      });
    }


    const downloadButton = el("analyse-download");
    if (downloadButton) downloadButton.addEventListener("click", download);

    const previewButton = el("analyse-preview");
    if (previewButton) previewButton.addEventListener("click", preview);

    const runButton = el("analyse-run");
    if (runButton) runButton.addEventListener("click", run);

    const resetButton = el("analyse-reset");
    if (resetButton) resetButton.addEventListener("click", resetEverything);

    loadSchema().then((schema) => {
      // Say where a plain name lands, rather than leaving somebody to guess
      // which directory "analysis_output" is relative to.
      const where = schema && schema.workspace;
      if (where) {
        text(
          el("analyse-output-note"),
          `A name lands in ${where}. A full path is used as given.`
        );
      }
      render();
    }).catch(() => {
      text(el("analyse-summary"), "Could not read the settings.");
    });
  }

  if (document.readyState === "loading") {
    document.addEventListener("DOMContentLoaded", attach);
  } else {
    attach();
  }

  // Exposed so the page's tests can drive it without a browser.
  window.FastMDXAnalyse = { state, currentState, optionsFor };
})();
