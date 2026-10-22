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

  function newestFirst(a, b) { return b.at - a.at; }

  /* The list as saved on this computer, newest first, or null where the
   * browser saves nothing. Only a run's end keeps its study's key (to open
   * it); nothing keeps a folder. */
  function stored() {
    try {
      var kept = JSON.parse(localStorage.getItem(SAVED) || "[]");
      var list = Array.isArray(kept) ? kept.filter(function (n) {
        return n && typeof n.said === "string" && isFinite(n.at);
      }) : [];
      list.forEach(function (n) {
        delete n.folder;
        if (n.act !== "overview") n.key = "";
      });
      return { list: list.sort(newestFirst), seen: Number(localStorage.getItem(SEEN)) || 0 };
    } catch (e) { return null; }
  }

  function read() {
    var kept = stored();
    if (!kept) return;
    notices = kept.list;
    seenAt = kept.seen;
  }

  function write() {
    try {
      localStorage.setItem(SAVED, JSON.stringify(notices.slice(0, MOST)));
      localStorage.setItem(SEEN, String(seenAt));
    } catch (e) { /* a browser that saves nothing: this page's list still holds them */ }
  }

  function board() { return window.FastMDXDashboard; }

  /* The open study: the server's key for it and its folder's name, as the
   * sidebar names it. */
  function study() {
    var state = board() && board().state ? board().state : {};
    var named = el("study-folder");
    var folder = String(state.outputDir || "");
    return { key: String(state.studyShown || ""),
             name: (named && named.dataset.name) || folder.replace(/[\\/]+$/, "").split(/[\\/]/).pop() || "" };
  }

  /* The study a toast is about: the open one on its pages; on a page of no
   * study (All studies, Compare, New study, the Agent) none, as a toast
   * there is about another study or none. */
  var STUDY_PAGES = ["overview", "viewer", "analysis", "report", "files"];
  function studyOfThePage() {
    var page = document.documentElement.dataset.page || "overview";
    return STUDY_PAGES.indexOf(page) >= 0 ? study() : { key: "", name: "" };
  }

  /* A notice saved, unless it is saved already: one with the same identity
   * (a run's end, a refusal of the Log, as another window of the GUI saved
   * it), or the same words of the same kind about the same study within two
   * seconds (a toast shown again as a page is drawn again). The list is
   * read from this computer first, where another window may have saved it.
   * Saved is what was said, the study's name and, only for what opens it
   * (a run's end), the study's key: nothing of the study the notice did not
   * show. */
  function add(kind, said, about, act, id) {
    said = String(said || "").replace(/\s+/g, " ").trim();
    if (!said) return;
    read();
    var now = Date.now();
    var where = about || studyOfThePage();
    var name = where.name || "";
    var already = notices.some(function (n) {
      if (id) return n.id === id && Math.abs(now - n.at) < 300000;
      return n.kind === kind && n.said === said && (n.name || "") === name && Math.abs(now - n.at) < 2000;
    });
    if (already) return;
    var notice = { at: now, kind: kind, said: said, name: name, act: act || "",
                   key: act === "overview" ? where.key || "" : "" };
    if (id) notice.id = id;
    notices.unshift(notice);
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

  // A run that ended, as run-notice.js sees it: how it ended, as the
  // browser's notice says it; the study by the server's key for it, named
  // by its folder.
  window.addEventListener("fmx:run-ended", function (event) {
    var d = event.detail || {};
    var ended = { failed: "failed", stopped: "stopped", interrupted: "was interrupted" }[d.run] || "completed";
    var said = (d.name || "The run") + " " + ended + (d.stage ? ": " + d.stage : ".");
    var here = study();
    var about = { key: d.study || "", name: d.folder || (d.study && d.study === here.key ? here.name : "") };
    add(d.run === "failed" || d.run === "interrupted" ? "refused" : "ended", said, about, "overview",
        "run:" + (d.study || "") + ":" + (d.started || "") + ":" + ended);
  });

  /* A refusal as a notice says it: its first line, at most 200
   * characters; what follows (a traceback, the paths in it) stays in the
   * Log. */
  function firstLine(text) {
    var line = String(text || "").split(/\r?\n/)[0].trim();
    return line.length > 200 ? line.slice(0, 199) + "\u2026" : line;
  }

  // A refusal the Log reads, written while the page is open: the events a
  // study had when it was first read are its past, the Log's to show.
  var since = {};
  window.addEventListener("fmx:events", function (event) {
    var d = event.detail || {};
    var where = study();
    var events = Array.isArray(d.events) ? d.events : [];
    var times = events.map(function (e) { return Date.parse(e.timestamp || ""); }).filter(isFinite);
    var newest = times.length ? Math.max.apply(null, times) : 0;
    if (!(where.key in since)) { since[where.key] = newest; return; }
    var from = since[where.key];
    events.forEach(function (e) {
      var at = Date.parse(e.timestamp || "");
      if (!isFinite(at) || at <= from || typeof d.refused !== "function" || !d.refused(e)) return;
      var line = firstLine(e.message);
      add("refused", line, where, "", "event:" + where.key + ":" + e.timestamp + ":" + line);
    });
    since[where.key] = Math.max(from, newest);
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
    // Keys pressed in the list stay in it: Escape closes the list alone,
    // and a page's own keys (All studies' c, j, k, t, x) do nothing under it.
    box.addEventListener("keydown", function (event) {
      if (event.ctrlKey || event.metaKey || event.altKey) return;
      event.stopPropagation();
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
      // Another study opens by a POST, refused to a page bound beyond
      // this computer: not offered there.
      if (n.act === "overview" && n.key && (n.key === study().key || !beyondLoopback())) {
        var go = make("button", "notice-act", "Open its Overview");
        go.type = "button";
        go.addEventListener("click", function () { openOverview(n.key); });
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
        // Newest first, with any saved since the list was cleared.
        notices = cleared.concat(notices).sort(newestFirst).slice(0, MOST);
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

  function beyondLoopback() { return document.body.dataset.beyondLoopback === "true"; }

  /* Another study opens by its key as the strip's tabs open one
   * (dashboard.js): asked first whether the switch would lose work here,
   * one switch at a time; why it did not open is said. */
  function openOverview(key) {
    var dashboard = board();
    if (!dashboard) return;
    close(false);
    if (study().key === key) { dashboard.navigate("overview"); return; }
    if (!dashboard.openStudyTab) return;
    dashboard.openStudyTab(key, "overview").then(function (answer) {
      if (answer && !answer.ok && dashboard.toast) {
        dashboard.toast(answer.error || "The study could not be opened.", "warning");
      }
    });
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
    seenAt = notices.reduce(function (most, n) { return Math.max(most, n.at); }, 0) || Date.now();
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
