/* Notices saved: every notice the page shows, every run that ended and
 * every refusal, in a list behind the bell in the status bar (IMDE 2.8).
 *
 * The toast is a notice's first sight and goes in seconds; what it said
 * was then lost, and a run that ended while the person was elsewhere said
 * so only in the tab's title. Each is saved now on this computer
 * (localStorage `fmx.notices`, the last 200, newest first): the toast as it
 * is shown (read from `#dashboard-toast`, so every caller's notice is
 * saved without each saving it), a run's end as `run-notice.js` sees it
 * (`fmx:run-ended`), and a refusal as the Log reads it from the study's
 * events (`fmx:events`), only those written while the page is open. The
 * bell says how many are new; its list has All, Runs and Refusals and
 * warnings; a notice keeps the study it was about and, for a run's end,
 * opens that study's Overview.
 */
(function () {
  "use strict";

  var SAVED = "fmx.notices";
  var SEEN = "fmx.noticesSeen";
  var MOST = 200;
  var notices = [];
  var seenAt = 0;
  var filter = "all";
  var cleared = null;

  function el(id) { return document.getElementById(id); }
  function make(tag, cls, text) {
    var node = document.createElement(tag);
    if (cls) node.className = cls;
    if (text != null) node.textContent = text;
    return node;
  }

  function read() {
    try {
      var kept = JSON.parse(localStorage.getItem(SAVED) || "[]");
      notices = Array.isArray(kept) ? kept.filter(function (n) { return n && typeof n.said === "string"; }) : [];
      seenAt = Number(localStorage.getItem(SEEN)) || 0;
    } catch (e) { notices = []; }
  }

  function write() {
    try {
      localStorage.setItem(SAVED, JSON.stringify(notices.slice(0, MOST)));
      localStorage.setItem(SEEN, String(seenAt));
    } catch (e) { /* a browser that saves nothing: this page's list still holds them */ }
  }

  function study() {
    var board = window.FastMDXDashboard;
    var folder = board && board.state ? String(board.state.outputDir || "") : "";
    return { folder: folder, name: folder.replace(/[\\/]+$/, "").split(/[\\/]/).pop() || "" };
  }

  /* A notice saved: the same words of the same kind within two seconds are
   * one (a toast shown again as a page is drawn again). */
  function add(kind, said, about, act) {
    said = String(said || "").replace(/\s+/g, " ").trim();
    if (!said) return;
    var now = Date.now();
    var last = notices[0];
    if (last && last.said === said && last.kind === kind && now - last.at < 2000) return;
    var where = about || study();
    notices.unshift({ at: now, kind: kind, said: said, folder: where.folder || "", name: where.name || "",
                      act: act || "" });
    if (notices.length > MOST) notices.length = MOST;
    write();
    drawBell();
    if (!panel().hidden) drawList();
  }

  /* ---- What is saved ---- */
  // Every toast, as it is shown.
  function watchToast() {
    var watched = null;
    var lastWords = null;
    function readToast() {
      var toast = el("dashboard-toast");
      var said = toast && toast.querySelector(".dashboard-toast-said");
      if (!said || !toast.classList.contains("show")) return;
      // Its words are set anew each time it is shown; anything else that
      // changes in it (its Undo going) is not a notice.
      if (said.firstChild === lastWords) return;
      lastWords = said.firstChild;
      var kind = toast.getAttribute("data-kind") || "ok";
      add(kind === "warning" || kind === "error" ? "warning" : "ok", said.textContent);
    }
    var reading = new MutationObserver(readToast);
    function attach() {
      var toast = el("dashboard-toast");
      if (!toast || toast === watched) return;
      watched = toast;
      reading.observe(toast, { childList: true, subtree: true, characterData: true,
                               attributes: true, attributeFilter: ["class"] });
      // Made for its first notice, already shown.
      readToast();
    }
    new MutationObserver(attach).observe(document.body, { childList: true });
    attach();
  }

  // A run that ended, as run-notice.js sees it.
  window.addEventListener("fmx:run-ended", function (event) {
    var d = event.detail || {};
    var ended = { failed: "failed", stopped: "stopped", interrupted: "was interrupted" }[d.run] || "completed";
    var said = (d.name || "The run") + " " + ended + (d.stage ? ": " + d.stage : ".");
    add(d.run === "failed" || d.run === "interrupted" ? "refused" : "ended", said,
        { folder: d.study || study().folder, name: d.name || study().name }, "overview");
  });

  // A refusal the Log reads, written while the page is open: the events a
  // study had when it was first read are its past, the Log's to show.
  var since = {};
  window.addEventListener("fmx:events", function (event) {
    var d = event.detail || {};
    var where = study();
    var events = Array.isArray(d.events) ? d.events : [];
    var times = events.map(function (e) { return Date.parse(e.timestamp || ""); }).filter(isFinite);
    var newest = times.length ? Math.max.apply(null, times) : 0;
    if (!(where.folder in since)) { since[where.folder] = newest; return; }
    var from = since[where.folder];
    events.forEach(function (e) {
      var at = Date.parse(e.timestamp || "");
      if (!isFinite(at) || at <= from || typeof d.refused !== "function" || !d.refused(e)) return;
      add("refused", e.message, where);
    });
    since[where.folder] = Math.max(from, newest);
  });

  /* ---- The bell and its list ---- */
  function unread() {
    return notices.filter(function (n) { return n.at > seenAt; }).length;
  }

  function drawBell() {
    var bell = el("sb-notices");
    if (!bell) return;
    var n = unread();
    var count = bell.querySelector(".sb-count");
    count.textContent = n > 99 ? "99+" : String(n);
    count.hidden = !n;
    bell.setAttribute("aria-label", "Notices" + (n ? ", " + n + " new" : ""));
    bell.title = "Every notice, refusal and run that ended" + (n ? ": " + n + " new" : "");
  }

  function panel() {
    var box = el("notices");
    if (box) return box;
    box = make("div", "notices");
    box.id = "notices";
    box.hidden = true;
    box.setAttribute("role", "dialog");
    box.setAttribute("aria-label", "Notices");
    box.addEventListener("keydown", function (event) {
      if (event.key === "Escape") { event.preventDefault(); close(true); }
    });
    document.body.appendChild(box);
    return box;
  }

  var KINDS = { all: "All", runs: "Runs", refused: "Refusals and warnings" };
  function shownNow() {
    return notices.filter(function (n) {
      return filter === "all" || (filter === "runs" ? n.act === "overview"
        : n.kind === "refused" || n.kind === "warning");
    });
  }

  function dayOf(at) {
    var d = new Date(at);
    var today = new Date();
    var yesterday = new Date(today.getTime() - 86400000);
    if (d.toDateString() === today.toDateString()) return "Today";
    if (d.toDateString() === yesterday.toDateString()) return "Yesterday";
    return d.toLocaleDateString(undefined, { weekday: "short", day: "numeric", month: "short" });
  }

  function timeOf(at) {
    var d = new Date(at);
    return String(d.getHours()).padStart(2, "0") + ":" + String(d.getMinutes()).padStart(2, "0");
  }

  function drawList() {
    var box = panel();
    box.replaceChildren();
    var head = make("div", "notices-head");
    head.appendChild(make("b", "", "Notices"));
    head.appendChild(make("span", "notices-n", notices.length + " saved"));
    box.appendChild(head);
    var chips = make("div", "notices-chips");
    chips.setAttribute("role", "group");
    chips.setAttribute("aria-label", "Show");
    Object.keys(KINDS).forEach(function (k) {
      var chip = make("button", "notices-chip", KINDS[k]);
      chip.type = "button";
      chip.dataset.show = k;
      chip.setAttribute("aria-pressed", String(filter === k));
      chip.addEventListener("click", function () { filter = k; drawList(); refocus("[data-show='" + k + "']"); });
      chips.appendChild(chip);
    });
    box.appendChild(chips);
    var list = make("div", "notices-list");
    var shown = shownNow();
    var day = "";
    shown.forEach(function (n) {
      var on = dayOf(n.at);
      if (on !== day) { day = on; list.appendChild(make("div", "notices-day", on)); }
      var row = make("div", "notice" + (n.at > seenAt ? " is-new" : ""));
      row.dataset.kind = n.kind;
      row.appendChild(make("span", "notice-dot"));
      row.appendChild(make("span", "notice-said", n.said));
      row.appendChild(make("span", "notice-when", timeOf(n.at)));
      var about = make("span", "notice-about");
      if (n.name) about.appendChild(make("span", "", n.name));
      if (n.act === "overview" && n.folder) {
        var go = make("button", "notice-act", "Open its Overview");
        go.type = "button";
        go.addEventListener("click", function () { openOverview(n.folder, go); });
        about.appendChild(go);
      }
      if (about.childNodes.length) row.appendChild(about);
      list.appendChild(row);
    });
    if (!shown.length) {
      list.appendChild(make("p", "notices-none", notices.length ? "None of these yet." :
        "Nothing yet. What the page tells you, each run that ends and each refusal are saved here."));
    }
    box.appendChild(list);
    var foot = make("div", "notices-foot");
    if (cleared) {
      foot.appendChild(make("span", "", "Cleared " + cleared.length + "."));
      var undo = make("button", "notices-undo", "Undo");
      undo.type = "button";
      undo.addEventListener("click", function () {
        notices = cleared.concat(notices).slice(0, MOST);
        cleared = null;
        write();
        drawBell();
        drawList();
        refocus(".notices-clear");
      });
      foot.appendChild(undo);
    } else {
      foot.appendChild(make("span", "", "Saved on this computer, the last " + MOST));
      var clear = make("button", "notices-clear", "Clear the list");
      clear.type = "button";
      clear.disabled = !notices.length;
      clear.addEventListener("click", function () {
        cleared = notices;
        notices = [];
        write();
        drawBell();
        drawList();
        refocus(".notices-undo");
      });
      foot.appendChild(clear);
    }
    box.appendChild(foot);
  }

  function refocus(selector) {
    var node = panel().querySelector(selector);
    if (node) node.focus();
  }

  function openOverview(folder, button) {
    var board = window.FastMDXDashboard;
    if (!board) return;
    close(false);
    if (study().folder === folder) { board.navigate("overview"); return; }
    button.disabled = true;
    fetch("/api/explore/switch", {
      method: "POST", headers: { "content-type": "application/json" },
      body: JSON.stringify({ folder: folder }),
    }).then(function (r) { return r.json(); }).then(function (d) {
      button.disabled = false;
      if (d && d.ok) board.studyOpened(d, "overview");
      else if (board.toast) board.toast((d && (d.error || d.reason)) || "The study could not open.", "warning");
    }).catch(function () { button.disabled = false; });
  }

  function open() {
    var box = panel();
    cleared = null;
    drawList();
    box.hidden = false;
    var bell = el("sb-notices");
    if (bell) bell.setAttribute("aria-expanded", "true");
    // Read once seen: the bell's count goes, the new ones stay marked in
    // the list while it is open.
    seenAt = notices.length ? notices[0].at : Date.now();
    write();
    drawBell();
    var first = box.querySelector(".notices-chip[aria-pressed='true']");
    if (first) first.focus();
  }

  function close(back) {
    var box = el("notices");
    if (!box || box.hidden) return;
    box.hidden = true;
    var bell = el("sb-notices");
    if (bell) {
      bell.setAttribute("aria-expanded", "false");
      if (back) bell.focus();
    }
  }

  function attach() {
    read();
    var bell = el("sb-notices");
    if (bell) {
      bell.hidden = false;
      bell.addEventListener("click", function () {
        if (panel().hidden) open(); else close(false);
      });
    }
    document.addEventListener("click", function (event) {
      var box = el("notices");
      // A click in the list that drew it again took its target out of the
      // page: still in the list.
      if (!event.target.isConnected) return;
      if (box && !box.hidden && !box.contains(event.target) && !event.target.closest("#sb-notices")) close(false);
    });
    // Another window of the GUI saved one.
    window.addEventListener("storage", function (event) {
      if (event.key !== SAVED && event.key !== SEEN) return;
      read();
      drawBell();
      if (!panel().hidden) drawList();
    });
    watchToast();
    drawBell();
  }

  if (document.readyState === "loading") document.addEventListener("DOMContentLoaded", attach);
  else attach();

  window.FastMDXNotices = { add: add, open: open, close: close, get all() { return notices.slice(); } };
}());
