/* A YAML draft never becomes the runnable Builder Config without explicit
 * validation and an independently checked conversion through its real mapper. */
(function () {
  "use strict";
  let api, area, dialog, status, opener, key = "", seen = null, drafted = false;
  let serial = 0, checked = null, busy = false, controller = null, diagnostics = [];
  const el = id => document.getElementById(id);
  const say = words => { status.textContent = words; };
  const textNode = (tag, words) => { const n = document.createElement(tag); n.textContent = words; return n; };
  function storageKey() {
    const workspace = window.FastMDXRun?.state?.schema?.workspace || api.state.app?.workspace || "";
    return "fmx.yamlDraft.v1:" + workspace + ":" + (api.state.studyShown || "new");
  }
  function persist() {
    if (!drafted || !key) return true;
    try {
      const current = localStorage.getItem(key);
      if (current !== seen) { say("Another window changed this YAML draft. Keep this window open or download your checked text before switching studies."); return false; }
      const value = JSON.stringify({text: area.value});
      localStorage.setItem(key, value); seen = value;
      return true;
    } catch (error) { say("This browser could not save the YAML draft. Keep this window open before switching studies."); return false; }
  }
  function invalidate() {
    serial += 1; checked = null; busy = false;
    controller?.abort(); controller = null;
    diagnostics = []; el("config-yaml-diagnostics").replaceChildren(); buttons(); gutter();
  }
  function hydrate() {
    const next = storageKey();
    if (next === key) return true;
    if (!persist()) return false;
    invalidate(); key = next; drafted = false; area.value = ""; seen = null;
    try {
      seen = localStorage.getItem(key);
      const saved = JSON.parse(seen || "null");
      if (saved && typeof saved.text === "string") { area.value = saved.text; drafted = true; }
    } catch (error) { say("The saved YAML draft could not be read. It has not been replaced."); }
    gutter(); return true;
  }
  function buttons() {
    const same = checked && checked.text === area.value;
    el("config-yaml-validate").disabled = busy;
    el("config-yaml-apply").disabled = busy || !same || !checked.answer.apply_ok;
    el("config-yaml-download").disabled = busy || !same;
  }
  function gutter() {
    if (!area) return;
    const top = Math.max(0, Math.floor(area.scrollTop / 22));
    const count = area.value.split("\n").length;
    const end = Math.min(count, top + Math.ceil(area.clientHeight / 22) + 2);
    const marked = new Set(diagnostics.flatMap(d => (d.locations || []).map(p => p.line)));
    const lines = el("config-yaml-lines"); lines.replaceChildren();
    for (let i = top; i < end; i += 1) {
      const n = textNode("span", String(i + 1));
      if (marked.has(i + 1)) n.className = "yaml-line-refused";
      lines.append(n);
    }
    lines.style.transform = "translateY(" + (12 - area.scrollTop % 22) + "px)";
  }
  function jump(place) {
    const rows = area.value.split("\n");
    const line = Math.min(rows.length, Math.max(1, place.line));
    const start = rows.slice(0, line - 1).reduce((n, row) => n + row.length + 1, 0);
    area.focus(); area.setSelectionRange(start, start + rows[line - 1].length);
    area.scrollTop = Math.max(0, (line - 3) * 22); gutter();
  }
  function showDiagnostics(items) {
    diagnostics = items || [];
    const host = el("config-yaml-diagnostics"); host.replaceChildren();
    diagnostics.forEach(item => {
      const row = textNode("li", "");
      (item.locations || []).forEach(place => {
        const button = textNode("button", "Line " + place.line);
        button.type = "button"; button.className = "ghost-btn";
        button.onclick = () => jump(place); row.append(button);
      });
      row.append(textNode("span", (item.code ? item.code + ": " : "") + item.message));
      if (!(item.locations || []).length) row.append(textNode("span", " (Document-level; no unique source line.)"));
      host.append(row);
    }); gutter();
  }
  async function fromBuilder() {
    invalidate(); const at = serial, where = key;
    try {
      const built = await window.FastMDXRun.fetchConfig(null, false);
      if (at !== serial || where !== key || dialog.hidden) return;
      area.value = built.ok ? built.yaml : ""; drafted = true;
      if (persist()) say(built.ok ? "Builder Config copied into a separate YAML draft. Validate after editing." : "The Builder is incomplete. Write a complete Config here, then validate it.");
      gutter();
    } catch (error) { if (at === serial) say("The Builder Config could not be read. Your saved draft is unchanged."); }
  }
  async function open() {
    if (!window.FastMDXRun?.ready) return;
    if (!hydrate()) return;
    window.FastMDXDialog.open("config-yaml-dialog", opener);
    el("config-yaml-replace-ask").hidden = true;
    if (!drafted) await fromBuilder();
    else say("Your separate YAML draft is restored. Validate this text before applying it.");
    gutter(); area.focus();
  }
  async function validate() {
    if (busy) return;
    invalidate(); busy = true; buttons();
    const at = serial, where = key, text = area.value;
    const request = controller = new AbortController();
    const stale = () => at !== serial || where !== key || area.value !== text || dialog.hidden;
    say("Checking YAML with the Config validator...");
    try {
      const ask = async body => {
        const response = await fetch("/api/config-yaml", {method: "POST", headers: {"Content-Type": "application/json"}, body: JSON.stringify(body), signal: request.signal});
        if (!response.ok) throw new Error("YAML validation is unavailable in this view.");
        return response.json();
      };
      let answer = await ask({text});
      if (stale()) return;
      if (answer.ok) {
        const run = window.FastMDXRun;
        const form = run.bodyOf(run.snapshotOfLoaded(answer.state));
        answer = await ask({text, form});
        if (stale()) return;
      }
      showDiagnostics(answer.diagnostics);
      if (answer.ok) {
        checked = {text, answer};
        say(answer.apply_ok ? "Valid YAML. Apply to Builder when ready; no run will start." : answer.mapping_reason || "Valid YAML can be downloaded; the Builder conversion could not be verified.");
      } else say(answer.error || "The Config validator refused this YAML.");
    } catch (error) { if (!stale() && error.name !== "AbortError") say(error.message); }
    finally { if (at === serial) { busy = false; controller = null; buttons(); } }
  }
  function attach() {
    api = window.FastMDXDashboard; area = el("config-yaml-text"); dialog = el("config-yaml-dialog");
    status = el("config-yaml-status"); opener = el("config-yaml-open");
    if (!api || !area || !window.FastMDXDialog) return;
    opener.onclick = open;
    area.addEventListener("input", () => { drafted = true; invalidate(); if (persist()) say("YAML draft saved. Validate the changed text before applying."); });
    area.addEventListener("scroll", gutter);
    new ResizeObserver(gutter).observe(area);
    el("config-yaml-validate").onclick = validate;
    el("config-yaml-apply").onclick = () => {
      if (!checked || !checked.answer.apply_ok || checked.text !== area.value || busy) return;
      window.FastMDXRun.applyLoadedState(checked.answer.state, {note: "Validated YAML applied to the Builder. No run was started."});
      window.FastMDXDialog.close("config-yaml-dialog");
    };
    el("config-yaml-download").onclick = () => {
      if (!checked || checked.text !== area.value || busy) return;
      const url = URL.createObjectURL(new Blob([checked.text], {type: "text/yaml"}));
      const a = document.createElement("a"); a.href = url; a.download = "fastmdxplora_edited.yml"; a.click();
      setTimeout(() => URL.revokeObjectURL(url), 1000);
    };
    el("config-yaml-replace").onclick = () => { el("config-yaml-replace-ask").hidden = false; el("config-yaml-replace-yes").focus(); };
    el("config-yaml-replace-no").onclick = () => { el("config-yaml-replace-ask").hidden = true; el("config-yaml-replace").focus(); };
    el("config-yaml-replace-yes").onclick = () => { el("config-yaml-replace-ask").hidden = true; fromBuilder(); };
    new MutationObserver(() => { if (dialog.hidden) { persist(); invalidate(); } }).observe(dialog, {attributes: true, attributeFilter: ["hidden"]});
    api.on("study-addressed", () => { if (key && storageKey() !== key) { window.FastMDXDialog.close("config-yaml-dialog"); hydrate(); } });
    window.FastMDXYamlEditor = {beforeStudySwitch: () => persist() ? "" : "Save or keep the YAML draft before switching studies."};
  }
  document.addEventListener("DOMContentLoaded", attach);
})();
