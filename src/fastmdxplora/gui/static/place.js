/* Where you are, at the head of every page: `villin › seed 1 › Analysis`.
 *
 * The study, then the run of a study of several (a switcher: its other runs
 * and the study as a whole, each with its state), then the page, whose
 * heading is the last part. A page of no study (All studies, New study, the
 * Agent) is its name alone. The sidebar's list of runs and its "Run of"
 * card were the same thing in two places, the run names cut to
 * "random_see…" there; they are here, once.
 *
 * Read from the app's state (exploration.py: `runs` at a study of several,
 * `run_of` in one of its runs, each run named by `exploration.run_label`).
 */
(function () {
  "use strict";

  var STUDY_PAGES = ["overview", "viewer", "analysis", "report", "files"];
  var DOT = { running: "status-dot-live", completed: "status-dot-completed",
              failed: "status-dot-error", stopped: "status-dot-waiting",
              interrupted: "status-dot-waiting", waiting: "status-dot-waiting" };
  var app = null;
  var shownKey = "";

  function el(tag, cls, text) {
    var node = document.createElement(tag);
    if (cls) node.className = cls;
    if (text !== undefined && text !== null) node.textContent = text;
    return node;
  }

  function page() { return document.documentElement.dataset.page || "overview"; }

  /* The open study's name as the sidebar's card says it. */
  function studyName() {
    var title = document.getElementById("topbar-run-title");
    var said = title ? title.textContent.trim() : "";
    return said && said !== "No active study" ? said : "";
  }

  /* What the place says now: the study, and the run where there are runs. */
  function place() {
    if (!app || !app.active_run || STUDY_PAGES.indexOf(page()) < 0) return null;
    var of = app.run_of;
    if (of) {
      return { study: of.study, studyPath: of.path, run: of.label || of.run_id,
               runs: of.runs || [], current: of.run_id, studyTitle: of.path };
    }
    var name = studyName();
    if (!name) return null;
    if (Array.isArray(app.runs)) {
      return { study: name, studyPath: app.active_run, run: "All runs", runs: app.runs,
               current: "", studyTitle: app.active_run };
    }
    return { study: name, studyPath: app.active_run, run: "", runs: null, current: "",
             studyTitle: app.active_run };
  }

  function stateSaid(run) {
    if (run.state === "running" && typeof run.fraction === "number") {
      return (100 * run.fraction).toFixed(0) + "%";
    }
    if (run.state === "waiting" && run.so_far) {
      return "next round" + (typeof run.fraction === "number"
        ? ", " + (100 * run.fraction).toFixed(0) + "%" : "");
    }
    return run.state || "";
  }

  /* Open a study's folder on the page shown, or the one given (a step of
   * its own). */
  function open(folder, onPage) {
    return fetch("/api/explore/switch", {
      method: "POST", headers: { "content-type": "application/json" },
      body: JSON.stringify({ folder: folder })
    }).then(function (r) { return r.json(); }).then(function (d) {
      if (d && d.ok && window.FastMDXDashboard) {
        window.FastMDXDashboard.studyOpened(d, onPage || page());
      }
      else if (window.FastMDXDashboard) {
        window.FastMDXDashboard.toast((d && d.error) || "That study could not be opened.", "error");
      }
    });
  }

  /* -------------------------------------------------------------------- */
  /* The run switcher's menu                                               */
  /* -------------------------------------------------------------------- */

  var menu = null;
  var menuFor = null;

  function closeMenu(focusBack) {
    if (!menu) return;
    menu.remove();
    menu = null;
    if (menuFor) {
      menuFor.setAttribute("aria-expanded", "false");
      if (focusBack) menuFor.focus();
    }
    menuFor = null;
  }

  function item(label, said, state, chosen, onChoose, title) {
    var row = el("button", "place-menu-item");
    row.type = "button";
    row.setAttribute("role", "menuitemradio");
    row.setAttribute("aria-checked", chosen ? "true" : "false");
    if (state) {
      var dot = el("span", "status-dot " + (DOT[state] || "status-dot-stale"));
      dot.setAttribute("aria-hidden", "true");
      row.appendChild(dot);
    } else {
      row.appendChild(el("span", "place-menu-blank"));
    }
    row.appendChild(el("span", "place-menu-label", label));
    row.appendChild(el("span", "place-menu-said mono", said));
    row.title = title || label;
    row.addEventListener("click", function () {
      closeMenu(false);
      if (!chosen) onChoose();
    });
    return row;
  }

  function openMenu(button, where) {
    closeMenu(false);
    menuFor = button;
    button.setAttribute("aria-expanded", "true");
    menu = el("div", "place-menu");
    menu.id = "place-menu";
    menu.setAttribute("role", "menu");
    menu.setAttribute("aria-label", "Runs of " + where.study);
    var runs = where.runs || [];
    var done = runs.filter(function (r) { return r.state === "completed"; }).length;
    menu.appendChild(item("All runs", done + " of " + runs.length + " completed", "",
      !where.current, function () { open(where.studyPath); }, "The study as a whole: its runs compared"));
    runs.forEach(function (run) {
      menu.appendChild(item(run.label || run.run_id, stateSaid(run), run.state,
        run.run_id === where.current, function () { open(run.path); },
        (run.label || run.run_id) + ", " + stateSaid(run) + "\n" + run.run_id));
    });
    document.body.appendChild(menu);
    var box = button.getBoundingClientRect();
    menu.style.left = Math.max(8, Math.min(box.left, innerWidth - menu.offsetWidth - 8)) + "px";
    menu.style.top = (box.bottom + 4) + "px";
    var items = Array.prototype.slice.call(menu.querySelectorAll(".place-menu-item"));
    var at = Math.max(0, items.findIndex(function (i) { return i.getAttribute("aria-checked") === "true"; }));
    items[at].focus();
    menu.addEventListener("keydown", function (e) {
      var now = items.indexOf(document.activeElement);
      if (e.key === "ArrowDown") { e.preventDefault(); items[(now + 1) % items.length].focus(); }
      else if (e.key === "ArrowUp") { e.preventDefault(); items[(now - 1 + items.length) % items.length].focus(); }
      else if (e.key === "Home") { e.preventDefault(); items[0].focus(); }
      else if (e.key === "End") { e.preventDefault(); items[items.length - 1].focus(); }
      else if (e.key === "Escape") { e.preventDefault(); closeMenu(true); }
      else if (e.key === "Tab") { closeMenu(false); }
    });
  }

  document.addEventListener("mousedown", function (e) {
    if (menu && !menu.contains(e.target) && e.target !== menuFor
        && !(menuFor && menuFor.contains(e.target))) closeMenu(false);
  });
  window.addEventListener("resize", function () { closeMenu(false); });

  /* -------------------------------------------------------------------- */
  /* The place, in each page's head                                        */
  /* -------------------------------------------------------------------- */

  function crumb(text, cls, title) {
    var b = el("button", "place-crumb " + cls, text);
    b.type = "button";
    if (title) b.title = title;
    return b;
  }

  function sep() {
    var s = el("span", "place-sep", "›");
    s.setAttribute("aria-hidden", "true");
    return s;
  }

  function fill(nav, where) {
    nav.innerHTML = "";
    var study = crumb(where.study, "place-study", "The study's Overview\n" + where.studyTitle);
    study.addEventListener("click", function () {
      if (where.current) open(where.studyPath, "overview");
      else if (window.FastMDXDashboard) window.FastMDXDashboard.navigate("overview");
    });
    nav.appendChild(study);
    nav.appendChild(sep());
    if (where.runs) {
      var run = crumb(where.run, "place-run", "Another run of this study");
      run.setAttribute("aria-haspopup", "menu");
      run.setAttribute("aria-expanded", "false");
      if (window.FastMDXIcons) {
        run.insertAdjacentHTML("beforeend", window.FastMDXIcons.svg("chevron", "place-chevron"));
      }
      run.addEventListener("click", function () {
        if (menu && menuFor === run) closeMenu(true);
        else openMenu(run, where);
      });
      run.addEventListener("keydown", function (e) {
        if (e.key === "ArrowDown") { e.preventDefault(); openMenu(run, where); }
      });
      nav.appendChild(run);
      nav.appendChild(sep());
    }
  }

  function render() {
    var where = place();
    var key = JSON.stringify([page(), where && [where.study, where.run, where.current,
      (where.runs || []).map(function (r) { return [r.run_id, r.label, r.state, r.fraction]; })]]);
    if (key === shownKey) return;
    shownKey = key;
    Array.prototype.forEach.call(document.querySelectorAll(".page[data-page] .page-header"), function (head) {
      var lead = head.firstElementChild;
      if (!lead) return;
      var nav = lead.querySelector(":scope > .place");
      var name = head.closest(".page").getAttribute("data-page");
      if (!where || STUDY_PAGES.indexOf(name) < 0) {
        if (nav) nav.hidden = true;
        return;
      }
      if (!nav) {
        nav = el("nav", "place");
        nav.setAttribute("aria-label", "Where you are");
        lead.insertBefore(nav, lead.firstChild);
      }
      nav.hidden = false;
      // The menu belongs to the place it was opened from; a new place
      // closes it, as the list it showed has changed.
      if (menu && nav.contains(menuFor)) closeMenu(false);
      fill(nav, where);
    });
  }

  function start() {
    if (!window.FastMDXDashboard) return;
    window.FastMDXDashboard.on("app-state", function (s) { app = s; render(); });
    window.FastMDXDashboard.on("status-updated", render);
    window.FastMDXDashboard.on("navigate", render);
  }

  if (document.readyState === "loading") document.addEventListener("DOMContentLoaded", start);
  else start();

  window.FastMDXPlace = { render: render, menuOpen: function () { return !!menu; } };
}());
