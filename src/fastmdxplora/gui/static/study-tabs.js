/* Open study documents in this browser window. Only keys and last pages
 * are kept here; records still come from the one active study on the
 * server. Closing a document never changes files or stops its run. */
(function () {
  "use strict";
  var api, host, list, note, storageKey = "", tabs = [], closed = [];
  var active = "", pending = false, hydrated = false;
  var lastRender = "", explicitKey = "";
  var focusOnOpen = "";
  var MAX_TABS = 100;
  function validKey(key) { return typeof key === "string" && /^[0-9a-f]{16}$/.test(key); }
  function safeRoute(route) {
    var match = /^#?([a-z-]+)(?:\/([a-zA-Z0-9_-]{1,120}))?$/.exec(String(route || ""));
    if (!match || api.state.pages.indexOf(match[1]) < 0) return "overview";
    return match[1] + (match[1] === "analysis" && match[2] ? "/" + match[2] : "");
  }
  function pageNow() {
    var state = api.state;
    return safeRoute(state.activePage + (state.activePlace ? "/" + state.activePlace : ""));
  }
  function studyPage() { return ["studies", "run", "compare"].indexOf(api.state.activePage) < 0; }
  function save() {
    if (!storageKey) return;
    try { sessionStorage.setItem(storageKey, JSON.stringify({tabs: tabs, closed: closed, focus: focusOnOpen})); }
    catch (error) { /* Browsing still works without storage. */ }
  }
  function hydrate(app) {
    var workspace = app && (app.exploration_root || app.workspace);
    if (!workspace) return;
    var key = "fmx.studyTabs:" + workspace;
    if (key === storageKey) return;
    storageKey = key;
    tabs = []; closed = []; focusOnOpen = ""; hydrated = true;
    try {
      var stored = JSON.parse(sessionStorage.getItem(key) || "{}");
      if (Array.isArray(stored.tabs)) stored.tabs.slice(0, MAX_TABS).forEach(function (item) {
        if (!item || !validKey(item.key) || tabs.some(function (t) { return t.key === item.key; })) return;
        tabs.push({key: item.key, title: typeof item.title === "string" ? item.title.slice(0, 200) : "Study",
          page: safeRoute(item.page)});
      });
      if (Array.isArray(stored.closed)) closed = stored.closed.filter(validKey).slice(-MAX_TABS);
      if (validKey(stored.focus)) focusOnOpen = stored.focus;
    } catch (error) { /* Invalid stored state is ignored. */ }
    render();
  }
  function remember(key, title, explicit) {
    if (!hydrated || pending || !validKey(key)) return;
    active = key;
    if (explicit || key === explicitKey) closed = closed.filter(function (k) { return k !== key; });
    explicitKey = "";
    if (closed.indexOf(key) >= 0) { render(); return; }
    var tab = tabs.find(function (t) { return t.key === key; });
    if (!tab) {
      tab = {key: key, title: "Study", page: "overview"};
      tabs.push(tab);
      if (tabs.length > MAX_TABS) tabs.shift();
    }
    if (typeof title === "string" && title.trim()) tab.title = title.slice(0, 200);
    if (studyPage()) tab.page = pageNow();
    save(); render();
    if (focusOnOpen === key) { focusOnOpen = ""; save(); focusTab(key); }
  }
  function rememberShown(explicit) {
    var state = api.state;
    if (state.studyIsOne === false) return;
    remember(state.studyShown, document.getElementById("topbar-run-title")?.textContent, explicit);
  }
  function hrefFor(tab) {
    var url = new URL(location.href);
    url.searchParams.set("study", tab.key);
    url.hash = safeRoute(tab.page);
    return url.pathname + url.search + url.hash;
  }
  function message(words) { note.textContent = words || ""; note.hidden = !words; }
  function forget(key) {
    tabs = tabs.filter(function (t) { return t.key !== key; });
    closed = closed.filter(function (k) { return k !== key; });
    closed.push(key);
    closed = closed.slice(-MAX_TABS);
    save();
  }
  function focusTab(key) {
    var links = Array.from(list.querySelectorAll(".study-tab-link"));
    var link = links.find(function (a) { return a.dataset.studyKey === key; }) || links[0];
    if (link) { link.focus(); link.scrollIntoView({block: "nearest", inline: "nearest"}); }
    else document.getElementById("study-tabs-open").focus();
  }
  async function activate(tab, closing) {
    if (pending) return;
    // Navigation announces its place on the next frame. Flush the current
    // place before pending suppresses that announcement.
    rememberShown(false);
    if (tab.key === active && !closing) { api.navigate(safeRoute(tab.page)); return; }
    pending = true; message("Opening " + tab.title + "..."); render();
    try {
      var answer = await api.openStudyTab(tab.key, safeRoute(tab.page), function () {
        if (closing) forget(closing);
        closed = closed.filter(function (key) { return key !== tab.key; });
        focusOnOpen = tab.key;
        save();
      });
      if (!answer || !answer.ok) {
        pending = false;
        message((answer && answer.error) || "Could not open this study. The current study stays open.");
        render(); focusTab(closing || tab.key);
      }
    } catch (error) {
      pending = false;
      message("Could not open this study. The current study stays open.");
      render(); focusTab(closing || tab.key);
    }
  }
  function close(tab) {
    if (pending) return;
    var index = tabs.indexOf(tab);
    var next = tabs[index + 1] || tabs[index - 1];
    if (tab.key !== active) {
      forget(tab.key); render(); focusTab(next && next.key); return;
    }
    if (next) { activate(next, tab.key); return; }
    // No study switch is needed to show All studies. The active study's
    // records remain open on the server, so a running simulation goes on.
    forget(tab.key); render(); api.navigate("studies"); focusTab("");
  }
  function render() {
    if (!host || !hydrated) return;
    var signature = JSON.stringify([tabs, active, pending]);
    if (signature === lastRender) return;
    lastRender = signature;
    var scroll = list.scrollLeft;
    var focused = document.activeElement;
    var focusKey = focused && focused.dataset.studyKey;
    var focusClose = focused && focused.classList.contains("study-tab-close");
    list.replaceChildren();
    tabs.forEach(function (tab) {
      var item = document.createElement("li");
      item.className = "study-tab";
      if (tab.key === active) item.classList.add("is-current");
      var link = document.createElement("a");
      link.className = "study-tab-link";
      link.dataset.studyKey = tab.key;
      link.href = hrefFor(tab);
      link.textContent = tab.title;
      link.title = tab.title;
      if (tab.key === active) link.setAttribute("aria-current", "page");
      if (pending) link.setAttribute("aria-disabled", "true");
      link.onclick = function (event) {
        if (pending) { event.preventDefault(); return; }
        if (event.metaKey || event.ctrlKey || event.shiftKey || event.altKey || event.button !== 0) return;
        event.preventDefault(); activate(tab);
      };
      var remove = window.FastMDXIcons.button("close", "Close " + tab.title, "study-tab-close");
      remove.dataset.studyKey = tab.key;
      remove.disabled = pending;
      remove.onclick = function () { close(tab); };
      item.append(link, remove); list.appendChild(item);
    });
    host.hidden = false;
    host.setAttribute("aria-busy", String(pending));
    list.scrollLeft = scroll;
    if (focusKey) {
      var selector = focusClose ? ".study-tab-close" : ".study-tab-link";
      var restored = Array.from(list.querySelectorAll(selector)).find(function (el) { return el.dataset.studyKey === focusKey; });
      if (restored) restored.focus({preventScroll: true});
    }
  }
  function boot() {
    api = window.FastMDXDashboard;
    host = document.getElementById("study-tabs");
    list = document.getElementById("study-tabs-list");
    note = document.getElementById("study-tabs-note");
    if (!api || !host || !list) return;
    if (performance.getEntriesByType?.("navigation")?.[0]?.type !== "reload") {
      explicitKey = new URL(location.href).searchParams.get("study") || "";
    }
    document.getElementById("study-tabs-open").onclick = function () { if (!pending) api.navigate("studies"); };
    list.onkeydown = function (event) {
      var key = event.target.dataset.studyKey;
      var tab = tabs.find(function (t) { return t.key === key; });
      if (!tab || pending) return;
      if (event.key === "Delete") { event.preventDefault(); close(tab); return; }
      var index = tabs.indexOf(tab), at = null;
      if (event.key === "ArrowRight") at = (index + 1) % tabs.length;
      if (event.key === "ArrowLeft") at = (index - 1 + tabs.length) % tabs.length;
      if (event.key === "Home") at = 0;
      if (event.key === "End") at = tabs.length - 1;
      if (at !== null) { event.preventDefault(); focusTab(tabs[at].key); }
    };
    api.on("app-state", function (app) {
      hydrate(app);
      if (!app.active_run && !pending) { active = ""; render(); }
    });
    api.on("study-addressed", function (study) {
      if (!study.key && !pending) { active = ""; render(); }
      else remember(study.key, study.title, false);
    });
    api.on("results-updated", function () { rememberShown(false); });
    api.on("study-opening", function (study) {
      if (validKey(study.key)) closed = closed.filter(function (key) { return key !== study.key; });
      save();
    });
    api.on("navigate", function () {
      if (studyPage()) rememberShown(true);
    });
    window.addEventListener("pagehide", save);
    window.addEventListener("pageshow", function (event) {
      if (event.persisted) {
        pending = false; message(""); storageKey = ""; lastRender = "";
        explicitKey = new URL(location.href).searchParams.get("study") || "";
        hydrate(api.state.appState);
      }
    });
    hydrate(api.state.appState);
  }
  document.addEventListener("DOMContentLoaded", boot);
})();
