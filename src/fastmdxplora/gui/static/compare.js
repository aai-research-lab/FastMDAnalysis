/* Compare: two studies, or two runs, side by side (compare_view.py).
 *
 * One view for what was in three places: the settings in which the two
 * differ, the means each recorded with the difference marked only where it
 * is resolved, each series both recorded over time overlaid, and the Viewer
 * playing the second beside the first. Opened from All studies (two
 * chosen), from a study of several runs (two of its runs chosen on its
 * Overview) and from the palette (Compare with...). Its address names the
 * two, so Back, Forward and a reload keep them: #compare/<first>|<second>.
 */
(function () {
  "use strict";

  /* The settings shown before "Show the other N", as All studies did. */
  var SHOWN_FIRST = 20;

  var shown = null;

  function byId(id) { return document.getElementById(id); }

  function make(tag, cls, text) {
    var node = document.createElement(tag);
    if (cls) node.className = cls;
    if (text !== undefined && text !== null) node.textContent = text;
    return node;
  }

  /* A mean with its error, to the error's second significant figure. */
  function withError(mean, error, unit) {
    if (mean === null || mean === undefined) return "not recorded";
    var said;
    if (error !== null && error !== undefined && error > 0) {
      var places = Math.max(0, Math.min(6, 1 - Math.floor(Math.log10(error))));
      said = mean.toFixed(places) + " \u00b1 " + error.toFixed(places);
    } else {
      said = Number(mean.toPrecision(4)).toString();
    }
    return said + (unit ? " " + unit : "");
  }

  function meanSaid(record) {
    if (!record) return "not recorded";
    if (record.withheld) return withError(record.mean, null, record.unit) + ", not determined";
    return withError(record.mean, record.error, record.unit);
  }

  function settingSaid(value) {
    return value == null ? "not set" : JSON.stringify(value);
  }

  /* The two the address names, or none. */
  function pairOf(hash) {
    var place = String(hash || "").replace(/^#/, "").split(/\/(.*)/)[1] || "";
    var parts = place.split("|");
    if (parts.length !== 2 || !parts[0] || !parts[1]) return null;
    try {
      return [decodeURIComponent(parts[0]), decodeURIComponent(parts[1])];
    } catch (error) {
      return null;
    }
  }

  function placeOf(first, second) {
    return encodeURIComponent(first) + "|" + encodeURIComponent(second);
  }

  /* Open the page on two folders (a step of its own). */
  function open(first, second) {
    if (!first || !second || !window.FastMDXDashboard) return;
    window.FastMDXDashboard.navigate("compare", { place: placeOf(first, second) });
    // What opened it is on a page now hidden: the keyboard and a screen
    // reader go on from the page's heading, as Back takes them there.
    var heading = document.querySelector('.page[data-page="compare"] .page-title');
    if (heading) {
      if (!heading.hasAttribute("tabindex")) heading.setAttribute("tabindex", "-1");
      heading.focus({ preventScroll: true });
    }
  }

  /* Nothing to compare, comparing, or not compared: one block, as every
   * page says it (page-states.js). */
  function say(state) {
    var empty = byId("compare-empty");
    if (!empty) return;
    if (!state) { empty.hidden = true; empty.replaceChildren(); return; }
    if (window.FastMDXPageState) window.FastMDXPageState.show(empty, state);
    else { empty.hidden = false; empty.textContent = state.title + (state.detail ? " " + state.detail : ""); }
  }

  function clear() {
    ["compare-heads", "compare-settings", "compare-means", "compare-series"].forEach(function (id) {
      var node = byId(id);
      if (node) { node.replaceChildren(); node.hidden = true; }
    });
    byId("compare-swap").hidden = true;
    byId("compare-beside").hidden = true;
    var page = document.querySelector('.page[data-page="compare"]');
    if (page) page.removeAttribute("data-compared");
  }

  var asking = "";

  function load() {
    var pair = pairOf(location.hash);
    // The page opened again on the two it shows, or is asking for: nothing
    // to read again.
    if (pair && (asking === location.hash || (shown && shown.hash === location.hash))) {
      return Promise.resolve(shown && shown.data);
    }
    clear();
    shown = null;
    if (!pair) {
      byId("compare-said").textContent = "Two studies, or two runs, side by side.";
      say({ kind: "empty", title: "Nothing to compare yet",
            detail: "Choose two studies in All studies, two runs on a study's Overview, or Compare " +
                    "with\u2026 in the palette (\u2318K)." });
      return Promise.resolve(null);
    }
    say({ kind: "loading", title: "Comparing\u2026" });
    var asked = location.hash;
    asking = asked;
    return fetch("/api/compared?" + new URLSearchParams({ a: pair[0], b: pair[1] }),
                 { cache: "no-store" })
      .then(function (r) { return r.json(); })
      .catch(function () { return { ok: false, reason: "The server did not answer." }; })
      .then(function (data) {
        if (asking === asked) asking = "";
        if (location.hash !== asked) return null;
        draw(pair, data);
        return data;
      });
  }

  function draw(pair, data) {
    var page = document.querySelector('.page[data-page="compare"]');
    if (!data || !data.ok) {
      say({ kind: "error", title: "They could not be compared",
            detail: (data && (data.reason || data.error)) || "The server gave no reason.",
            action: { label: "Try again", run: function () { shown = null; asking = ""; load(); } } });
      if (page) page.setAttribute("data-compared", "no");
      return;
    }
    shown = { pair: pair, data: data, hash: location.hash };
    say(null);
    byId("compare-said").textContent = data.first.said + " and " + data.second.said;
    byId("compare-swap").hidden = false;
    byId("compare-beside").hidden = false;
    heads(data);
    settings(data);
    means(data);
    series(data);
    if (page) page.setAttribute("data-compared", "yes");
  }

  /* Each of the two: its name, its system and state, and the way to open
   * it. */
  function heads(data) {
    var host = byId("compare-heads");
    host.hidden = false;
    [data.first, data.second].forEach(function (side, i) {
      var card = make("div", "card compare-head");
      card.dataset.side = i ? "second" : "first";
      var swatch = make("span", "compare-swatch");
      swatch.style.background = colourOf(i);
      var title = make("div", "compare-head-title");
      title.appendChild(swatch);
      title.appendChild(make("strong", "compare-head-name", side.said));
      card.appendChild(title);
      var facts = [side.system, side.kind, side.state].filter(Boolean).join(" \u00b7 ");
      if (facts) card.appendChild(make("div", "muted small", facts));
      if (side.run_of) {
        card.appendChild(make("div", "muted small", "A run of " + side.run_of.name + ": " +
                                                    side.run_of.label));
      }
      var go = make("button", "chip-btn compare-open", "Open");
      go.type = "button";
      go.title = "Open " + side.said + " on its Overview";
      go.addEventListener("click", function () { switchTo(side.path, "overview"); });
      card.appendChild(go);
      host.appendChild(card);
    });
  }

  function settings(data) {
    var host = byId("compare-settings");
    host.hidden = false;
    host.appendChild(make("h2", "card-title", "Settings that differ"));
    var found = data.settings || [];
    if (!found.length) {
      host.appendChild(make("p", "muted small", "None: they asked for the same study."));
    } else {
      var table = make("table", "kv-table compare-settings-table");
      var top = make("tr");
      ["Setting", data.first.said, data.second.said].forEach(function (label) {
        top.appendChild(make("th", "", label));
      });
      table.appendChild(top);
      found.forEach(function (d, i) {
        var row = make("tr");
        if (i >= SHOWN_FIRST) row.hidden = true;
        row.appendChild(make("td", "mono", d.setting));
        row.appendChild(make("td", "mono", settingSaid(d.first)));
        row.appendChild(make("td", "mono", settingSaid(d.second)));
        table.appendChild(row);
      });
      host.appendChild(table);
      if (found.length > SHOWN_FIRST) {
        var more = make("button", "file-action", "Show the other " + (found.length - SHOWN_FIRST));
        more.type = "button";
        more.addEventListener("click", function () {
          table.querySelectorAll("tr[hidden]").forEach(function (row) { row.hidden = false; });
          more.remove();
        });
        host.appendChild(more);
      }
    }
    var unrecorded = data.unrecorded || [];
    if (unrecorded.length) {
      host.appendChild(make("p", "muted small",
        "Not compared: " + unrecorded.join(", ") + ". The run works " +
        (unrecorded.length === 1 ? "it" : "these") + " out from the other settings, " +
        "and one of the two never recorded " + (unrecorded.length === 1 ? "it" : "them") +
        " (it stopped or failed before them, or analysed a trajectory and ran no simulation)."));
    }
  }

  function means(data) {
    var host = byId("compare-means");
    host.hidden = false;
    host.appendChild(make("h2", "card-title", "What they recorded"));
    var found = data.measures || [];
    if (!found.length) {
      host.appendChild(make("p", "muted small", "Neither recorded a mean."));
      return;
    }
    var table = make("table", "kv-table compare-means-table");
    var top = make("tr");
    ["Quantity", data.first.said, data.second.said, "Difference"].forEach(function (label) {
      top.appendChild(make("th", "", label));
    });
    table.appendChild(top);
    found.forEach(function (m) {
      var row = make("tr");
      row.dataset.analysis = m.analysis;
      row.appendChild(make("td", "", m.label || m.analysis));
      row.appendChild(make("td", "mono", meanSaid(m.first)));
      row.appendChild(make("td", "mono", meanSaid(m.second)));
      var cell = make("td", "mono");
      if (m.versus) {
        cell.textContent = withError(m.versus.difference, m.versus.error, m.unit) +
          (m.versus.resolved ? ", resolved" : ", within the errors");
        if (m.versus.resolved) row.classList.add("is-resolved");
      } else {
        cell.textContent = "not compared";
      }
      row.appendChild(cell);
      table.appendChild(row);
    });
    host.appendChild(table);
    host.appendChild(make("p", "muted small",
      "The difference is the second's mean less the first's. It is resolved past " +
      data.resolved_at + " times the two's combined standard error. A mean an analysis " +
      "did not stand behind is not compared."));
  }

  function colourOf(i) {
    var runs = window.FastMDXRunsCompared;
    return runs && runs.colourOf ? runs.colourOf(i) : (i ? "#E69F00" : "#0072B2");
  }

  /* Each series both recorded, the two overlaid on one clock. */
  function series(data) {
    var host = byId("compare-series");
    host.hidden = false;
    host.appendChild(make("h2", "card-title", "Over time"));
    var found = data.series || [];
    if (!found.length) {
      host.appendChild(make("p", "muted small", "They recorded no series over time in common."));
      return;
    }
    var runs = window.FastMDXRunsCompared;
    found.forEach(function (s) {
      var block = make("div", "compare-series-block");
      block.dataset.analysis = s.analysis;
      var chart = make("div", "compare-series-chart");
      var legend = make("div", "runs-overlay-legend compare-series-legend");
      block.appendChild(chart);
      block.appendChild(legend);
      host.appendChild(block);
      if (s.first.x_label !== s.second.x_label) {
        // One on production's clock, the other on frame numbers: drawn on
        // one axis, they would say the frames line up in time.
        chart.appendChild(make("p", "muted small", (s.label || s.analysis) + ": not overlaid. " +
          "One is recorded against " + s.first.x_label.toLowerCase() + ", the other against " +
          s.second.x_label.toLowerCase() + "."));
        return;
      }
      if (!runs || !runs.drawOverlay) return;
      var drawn = [data.first, data.second].map(function (side, i) {
        var line = i ? s.second : s.first;
        return { r: { run: { run_id: i ? "second" : "first", label: side.said }, index: i },
                 series: { x: line.x, y: line.y, x_label: line.x_label } };
      });
      runs.drawOverlay(chart, legend, drawn, { label: s.label, unit: s.unit }, 0);
      chart.setAttribute("data-drawn", "2");
    });
  }

  /* Open a folder on a page: the server's study switched, the page loaded
   * again on it. */
  function switchTo(folder, page) {
    return fetch("/api/explore/switch", {
      method: "POST", headers: { "content-type": "application/json" },
      body: JSON.stringify({ folder: folder })
    }).then(function (r) { return r.json(); }).then(function (d) {
      if (d && d.ok && window.FastMDXDashboard) window.FastMDXDashboard.studyOpened(d, page);
      else if (window.FastMDXDashboard) {
        window.FastMDXDashboard.toast((d && d.error) || "That study could not be opened.", "error");
      }
      return d;
    });
  }

  /* The first in the Viewer, the second played beside it (viewer-beside.js).
   * Opening the first loads the page again; what to show beside it is kept
   * for the Viewer to take up then. */
  function playBeside() {
    if (!shown) return;
    var first = shown.pair[0], second = shown.pair[1];
    var beside = window.FastMDXBeside;
    var here = "";
    try {
      here = window.FastMDXDashboard.state.appState.active_run || "";
    } catch (error) {
      here = "";
    }
    if (beside && beside.next) beside.next(second, shown.data.second.said);
    if (here === first) window.FastMDXDashboard.navigate("viewer");
    else switchTo(first, "viewer");
  }

  function swap() {
    if (shown) open(shown.pair[1], shown.pair[0]);
  }

  function wire() {
    if (!byId("compare-heads")) return;
    byId("compare-swap").addEventListener("click", swap);
    byId("compare-beside").addEventListener("click", playBeside);
    window.addEventListener("dashboard:navigate", function (event) {
      if (event.detail && event.detail.page === "compare") load();
    });
  }

  if (document.readyState === "loading") document.addEventListener("DOMContentLoaded", wire);
  else wire();

  window.FastMDXCompare = { open: open, load: load, pairOf: pairOf };
}());
