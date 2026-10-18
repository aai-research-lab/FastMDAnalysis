/* Recorded problems in the study on screen. This panel only reads; the
 * existing Overview keeps the confirmed actions for running a remedy.
 *
 * Read when the study opens, when the tab is shown, with Refresh, and again
 * only when the run writes something new: its live record's update time or
 * its events log changed, as the page learns them from its status poll or
 * event stream (user, 10-18). No timer reads it in the background; the
 * status bar's count is this read's. Refresh reads it at once. */
(function () {
  "use strict";

  var active = null;
  var study = "";
  var generation = 0;
  var pending = null;
  var rendered = "";
  // What the run had written when last read: its live record's update
  // time, and its events log's newest line. null until first seen.
  var liveSeen = null;
  var eventsSeen = null;

  function byId(id) { return document.getElementById(id); }
  function element(tag, cls, text) {
    var node = document.createElement(tag);
    if (cls) node.className = cls;
    if (text != null) node.textContent = text;
    return node;
  }
  function visible() { return byId("problems-pane") && !byId("problems-pane").hidden; }
  /* The count, for the status bar (statusbar.js); null where not known. */
  function counted(count, partial) {
    window.dispatchEvent(new CustomEvent("dashboard:problems-counted",
      { detail: { count: count, partial: Boolean(partial) } }));
  }
  function message(words) {
    counted(null);
    byId("problems-summary").textContent = words;
    byId("problems-body").replaceChildren();
    rendered = "";
  }
  function clear() {
    generation += 1;
    if (pending) pending.abort();
    pending = null;
    again = false;
    byId("problems-refresh").disabled = false;
    study = "";
    liveSeen = null;
    eventsSeen = null;
    message(active ? "Reading the study's records..." : "Open a study to read its problems.");
  }
  function render(data) {
    var key = JSON.stringify(data);
    if (key === rendered) return;
    rendered = key;
    var body = byId("problems-body");
    body.replaceChildren();
    var items = data.items || [];
    var read = data.read || [];
    var missing = (data.unavailable || []).length > 0;
    /* For the status bar: none is said only of every record read, as the
     * headline says it; where a record could not be read, those found. */
    var partial = !data.available || missing;
    counted(partial && !items.length ? null : items.length, partial);
    var found = items.length + (items.length === 1 ? " recorded problem" : " recorded problems");
    /* "No problems" is said only of every record read. Where one could
     * not be, the headline says so, and the notes name it. */
    byId("problems-summary").textContent = !data.available
      ? (items.length ? found + "; other records are unavailable." : "Records unavailable.")
      : items.length
        ? found + (missing ? "; some records are unavailable." : "")
        : missing
          ? "Some records are unavailable. The " + read.length +
            (read.length === 1 ? " record read lists" : " records read list") + " no problems."
          : "No problems listed in the records read.";
    items.forEach(function (item) {
      var card = element("article", "problem-item");
      card.dataset.problemId = item.id;
      card.appendChild(element("h3", "", item.title));
      card.appendChild(element("p", "", item.reason));
      if (item.code) card.appendChild(element("code", "problem-code", item.code));
      card.appendChild(element("h4", "", "What would fix it"));
      card.appendChild(element("p", "", item.fix || "No remedy is recorded."));
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
    if (read.length) {
      body.appendChild(element("p", "problems-note problems-read", "Records read: " + read.join(", ") + "."));
    }
  }
  // Asked while a read was under way: read once more after it, since the
  // read under way may have begun before what is new was written.
  var again = false;
  function load() {
    if (!active || !study) return;
    if (pending) {
      again = true;
      return;
    }
    again = false;
    var at = generation;
    var wanted = study;
    var root = active;
    var controller = new AbortController();
    pending = controller;
    byId("problems-refresh").disabled = true;
    fetch("/api/problems?study=" + encodeURIComponent(wanted), {
      cache: "no-store", signal: controller.signal
    }).then(function (response) {
      if (!response.ok) throw new Error("Could not read the study's problems.");
      return response.json();
    }).then(function (data) {
      if (at !== generation || root !== active) return;
      /* No study open has no key. A folder open that is not yet a study
       * (one a run has just begun in) has its key, and no run named. */
      if (data.ok && data.active_run === null && !data.study) {
        message("Open a study to read its problems.");
        return;
      }
      if (data.study !== wanted || (data.active_run !== null && data.active_run !== root)) {
        message("The open study changed. Waiting for its records...");
        return;
      }
      if (!data.ok) throw new Error(data.reason || "Could not read the study's problems.");
      if (data.state === "gone" || data.state === "unreadable" || data.state === "empty") {
        /* Refresh cannot bring a folder back: said once, without it. A
         * folder with no records yet is said as that, not as missing ones. */
        message((data.notes || [])[0] || "The study's records cannot be read.");
        return;
      }
      render(data);
    }).catch(function (error) {
      if (at !== generation || error.name === "AbortError") return;
      message("Could not read problems. Use Refresh to try again.");
    }).finally(function () {
      if (pending !== controller) return;
      pending = null;
      byId("problems-refresh").disabled = false;
      if (again && at === generation) load();
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
      load();
    });
    dialog.addEventListener("close", function () {
      home.insertBefore(panel, next);
      phone.focus();
    });
    window.matchMedia("(max-width: 800px)").addEventListener("change", function (event) {
      if (!event.matches && dialog.open) dialog.close();
    });
    byId("problems-refresh").addEventListener("click", load);
    /* Shown, however it was (its tab clicked or chosen by the keys, the
     * status bar's count, a phone's button): read again. */
    var pane = byId("problems-pane");
    var wasShown = visible();
    new MutationObserver(function () {
      var shown = visible();
      if (shown && !wasShown) load();
      wasShown = shown;
    }).observe(pane, { attributes: true, attributeFilter: ["hidden"] });
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
      var status = (update && update.status) || {};
      var live = String(status.last_update_timestamp || "") + "|" + String(status.status || "");
      if (live === liveSeen) return;
      liveSeen = live;
      load();
    });
    window.FastMDXDashboard.on("events-read", function (detail) {
      if (!active || !study) return;
      var events = (detail && detail.events) || [];
      var last = events.length ? events[events.length - 1] : {};
      var key = events.length + "|" + String(last.timestamp || "") + "|" + String(last.message || "");
      // The first seen of a study is read with it; only a change after reads again.
      var changed = eventsSeen !== null && key !== eventsSeen;
      eventsSeen = key;
      if (changed) load();
    });
  }
  document.addEventListener("DOMContentLoaded", boot);
}());
