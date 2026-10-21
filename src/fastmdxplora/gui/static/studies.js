/* Every study in the workspace, as cards, and two of them compared
 * (workspace.py).
 *
 * Each card says the study's structure, what kind of study it is, where it
 * stands, when it began, the means it recorded with their errors and a
 * figure it plotted. Search narrows the cards as typed. Two chosen open
 * on the Compare page (compare.js). Open makes a study the one on
 * screen. Everything is built as elements: a study's name and its
 * records are data. */
(function () {
  "use strict";

  var studies = [];
  /* The page's own line icons (icons.js, from sidebar_icons.py). */
  var ICONS = window.FastMDXIcons;
  var chosen = [];
  var lookedIn = "";
  var loading = null;
  // Tags: the one the cards are narrowed to, the card being tagged, and
  // the tags the studies here carry (offered as one is typed).
  var tagged = null;
  var editing = null;
  var tagsUsed = [];
  var TAG_LENGTH = 40;
  var NOTE_LENGTH = 200;

  function el(id) { return document.getElementById(id); }

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
      said = mean.toFixed(places) + " ± " + error.toFixed(places);
    } else {
      // Four figures kept as figures: 1.000 nm, not 1 nm, which said less
      // than the record holds.
      said = Math.abs(mean) >= 1e4 ? Math.round(mean).toLocaleString("en-US") : mean.toPrecision(4);
      if (/e/.test(said)) said = Number(said).toString();
    }
    return said + (unit ? " " + unit : "");
  }

  /* A time within the day said as how long ago; earlier, its date. */
  function timeSaid(when) {
    var date = new Date(when);
    if (isNaN(date.getTime())) return "";
    var minutes = Math.round((Date.now() - date.getTime()) / 60000);
    if (minutes < 1) return "Just now";
    if (minutes < 60) return minutes + " min ago";
    if (minutes < 24 * 60) return Math.round(minutes / 60) + " h ago";
    return dateSaid(when);
  }

  /* Tags and a note of the person's own, changed on the card or under
   * the table's row. */
  /* A page bound beyond loopback (the page's bind flag): the server
   * refuses every POST there but the Config's, so Tag, Archive, Open and
   * Files, which POST, and Compare, which is read on loopback only, are
   * not offered. */
  function beyondLoopback() {
    return document.body.dataset.beyondLoopback === "true";
  }

  function tagButton(study) {
    if (beyondLoopback()) return null;
    var count = (study.tags || []).length;
    var some = count || study.note;
    // Tags once it has some, marked so on the icon; the count to a reader.
    var tag = ICONS.button("tag", some ? "Tags" + (count ? " (" + count + ")" : "") : "Tag",
                           "study-tag-edit" + (some ? " has-tags" : ""));
    tag.title = (some ? "Tags" : "Tag") + ": words and a note of your own, kept in the study's folder";
    tag.addEventListener("click", function () {
      editing = study.path;
      draw();
      var input = document.querySelector(".study-tag-editor .study-tag-input");
      if (input) input.focus();
    });
    return tag;
  }

  /* A study's system as it is shown: the ID the person gave it in this
   * browser (preferences.js, by its folder), else its own. */
  function shownId(study) {
    try {
      var kept = JSON.parse(localStorage.getItem("fmx.study:" + study.path) || "null");
      if (kept && /^[A-Z0-9]{1,4}$/.test(kept.name || "")) return kept.name;
    } catch (e) { /* not kept */ }
    return study.system || "";
  }

  function dateSaid(when) {
    var date = new Date(when);
    return isNaN(date.getTime()) ? "" : date.toLocaleDateString(undefined,
      { year: "numeric", month: "short", day: "numeric" });
  }

  function searchText(study) {
    return [study.name, shownId(study), study.system, study.structure || "", study.kind,
      study.state, study.forcefield || "",
      (study.tags || []).join(" "), study.note || ""].join(" ").toLowerCase();
  }

  function hasTag(study, tag) {
    var wanted = tag.toLowerCase();
    return (study.tags || []).some(function (t) { return t.toLowerCase() === wanted; });
  }

  /* The tags a card carries, each narrowing the cards to it, and its note. */
  function tagsShown(study) {
    var holder = make("div", "study-tags");
    (study.tags || []).forEach(function (tag) {
      var chip = make("button", "study-tag", tag);
      chip.type = "button";
      chip.title = "Show only the studies tagged " + tag;
      chip.addEventListener("click", function () { narrowTo(tag); });
      holder.appendChild(chip);
    });
    return holder;
  }

  function narrowTo(tag) {
    tagged = tag;
    draw();
  }

  /* The card's tags and note, changed on it and kept in the study's folder. */
  function tagEditor(study) {
    var form = make("form", "study-tag-editor");
    form.setAttribute("aria-label", "Tags of " + study.name);
    var tags = (study.tags || []).slice();
    var list = make("div", "study-tags");
    function drawTags() {
      list.replaceChildren.apply(list, tags.map(function (tag, i) {
        var chip = make("span", "study-tag is-editing", tag);
        var drop = make("button", "study-tag-drop", "\u00d7");
        drop.type = "button";
        drop.setAttribute("aria-label", "Remove the tag " + tag);
        drop.addEventListener("click", function () {
          tags.splice(i, 1);
          drawTags();
        });
        chip.appendChild(drop);
        return chip;
      }));
    }
    drawTags();
    var adding = make("input", "study-tag-input");
    adding.type = "text";
    adding.maxLength = TAG_LENGTH;
    adding.setAttribute("list", "studies-tags-used");
    adding.placeholder = "Add a tag";
    adding.setAttribute("aria-label", "A tag to add");
    function add() {
      var tag = adding.value.replace(/\s+/g, " ").trim();
      if (tag && !tags.some(function (t) { return t.toLowerCase() === tag.toLowerCase(); })) {
        tags.push(tag);
        drawTags();
      }
      adding.value = "";
    }
    adding.addEventListener("keydown", function (event) {
      if (event.key === "Enter") {
        event.preventDefault();
        add();
      }
    });
    var addButton = make("button", "file-action", "Add");
    addButton.type = "button";
    addButton.addEventListener("click", add);
    var note = make("input", "study-note-input");
    note.type = "text";
    note.maxLength = NOTE_LENGTH;
    note.value = study.note || "";
    note.placeholder = "A note: one line";
    note.setAttribute("aria-label", "A note on the study");
    var said = make("p", "muted small study-tag-said");
    said.setAttribute("role", "status");
    var save = make("button", "file-action", "Save");
    save.type = "submit";
    var cancel = make("button", "file-action", "Cancel");
    cancel.type = "button";
    cancel.addEventListener("click", function () {
      editing = null;
      draw();
    });
    form.addEventListener("submit", function (event) {
      event.preventDefault();
      if (adding.value.trim()) add();
      save.disabled = true;
      var before = { tags: (study.tags || []).slice(), note: study.note || "" };
      // The change, not the whole list: what was put on and taken off
      // here, applied by the server to the tags as kept, so a tag an AI
      // app or another window kept since this list was read stays.
      var change = tagsChanged(before.tags, tags);
      before.change = { add: change.remove, remove: change.add };
      fetch("/api/study-tags", {
        method: "POST",
        headers: { "content-type": "application/json" },
        body: JSON.stringify({ path: study.path, add: change.add, remove: change.remove,
                               note: note.value })
      }).then(function (r) { return r.json(); }).then(function (d) {
        save.disabled = false;
        if (!d || !d.ok) {
          said.textContent = (d && (d.reason || d.error)) || "The tags could not be kept.";
          return;
        }
        study.tags = d.tags;
        study.note = d.note;
        edited += 1;
        editing = null;
        tagsUsed = countTags();
        offerTags();
        draw();
        tagsOfFocused(study);
        offerUndoOfTags(study, before);
      }).catch(function () {
        save.disabled = false;
        said.textContent = "The server did not answer.";
      });
    });
    var row = make("div", "study-tag-row");
    row.append(adding, addButton);
    var buttons = make("div", "study-tag-row");
    buttons.append(save, cancel);
    form.append(list, row, note, buttons, said);
    return form;
  }

  /* Tags and a note kept, said in the page's notice with Undo, which keeps
   * them as they were. */
  // The keyboard on a card's Tags, drawn anew (where the person was, and
  // where the notice's Undo brings it back to).
  function tagsOfFocused(study) {
    Array.prototype.forEach.call(document.querySelectorAll(".study-card"), function (card) {
      if (card.dataset.path === study.path) {
        var edit = card.querySelector(".study-tag-edit");
        if (edit) edit.focus({ preventScroll: true });
      }
    });
  }

  /* What was put on (``add``) and taken off (``remove``) between two lists
   * of tags, a tag whatever its case. */
  function tagsChanged(was, now) {
    var key = function (t) { return String(t).toLowerCase(); };
    var had = was.map(key), has = now.map(key);
    return { add: now.filter(function (t) { return had.indexOf(key(t)) < 0; }),
             remove: was.filter(function (t) { return has.indexOf(key(t)) < 0; }) };
  }

  function offerUndoOfTags(study, before) {
    var board = window.FastMDXDashboard;
    if (!board || !board.toast) return;
    board.toast("Kept the tags and note of " + (study.name || "the study") + ".", "ok", { label: "Undo", run: function () {
      // Settled once its answer is said: a notice waiting on the Undo comes after it.
      return fetch("/api/study-tags", {
        method: "POST",
        headers: { "content-type": "application/json" },
        body: JSON.stringify({ path: study.path, add: before.change.add,
                               remove: before.change.remove, note: before.note })
      }).then(function (r) { return r.json(); }).then(function (d) {
        if (!d || !d.ok) {
          board.toast((d && (d.reason || d.error)) || "It could not be undone.", "warning");
          return;
        }
        // The card as listed now: the list is built again on each visit,
        // and the old one's tags edited next would be saved back.
        studies.filter(function (one) { return one.path === study.path; })
          .concat([study]).forEach(function (one) {
            one.tags = d.tags;
            one.note = d.note;
          });
        edited += 1;
        tagsUsed = countTags();
        offerTags();
        draw();
        // The card drawn anew took the keyboard with it: back on its Tags.
        if (document.activeElement === document.body) tagsOfFocused(study);
        board.toast("The tags and note are back as they were.", "ok");
      }).catch(function () { board.toast("The server did not answer.", "warning"); });
    } });
  }

  function countTags() {
    var counted = {};
    studies.forEach(function (study) {
      (study.tags || []).forEach(function (tag) {
        var key = tag.toLowerCase();
        if (!counted[key]) counted[key] = { tag: tag, studies: 0 };
        counted[key].studies += 1;
      });
    });
    return Object.keys(counted).map(function (key) { return counted[key]; })
      .sort(function (a, b) { return b.studies - a.studies || a.tag.localeCompare(b.tag); });
  }

  /* The tags used here, offered as a tag is typed. */
  function offerTags() {
    var listed = el("studies-tags-used");
    if (!listed) return;
    listed.replaceChildren.apply(listed, tagsUsed.map(function (used) {
      var option = make("option");
      option.value = used.tag;
      return option;
    }));
  }

  /* The page's scheme, for a picture drawn in it. */
  function scheme() {
    return document.documentElement.dataset.theme === "light" ? "light" : "dark";
  }

  var SVG = "http://www.w3.org/2000/svg";
  function svg(tag, attrs, parent) {
    var node = document.createElementNS(SVG, tag);
    Object.keys(attrs).forEach(function (key) { node.setAttribute(key, attrs[key]); });
    if (parent) parent.appendChild(node);
    return node;
  }

  /* ---- The molecule ------------------------------------------------- */
  /* The study's backbone (workspace.py `molecule`, backbone_picture.py
   * `backbone_points`) drawn at the size asked: each piece's points in
   * the plane of its two widest axes, coloured from the N terminus to the
   * C terminus in viridis (its lightest fifth left out on paper, its
   * darkest on the dark scheme), the nearer part wider and stronger. */
  var VIRIDIS = [[0, [0x44, 0x01, 0x54]], [0.2, [0x41, 0x44, 0x87]], [0.4, [0x2a, 0x78, 0x8e]],
                 [0.6, [0x22, 0xa8, 0x84]], [0.8, [0x7a, 0xd1, 0x51]], [1, [0xfd, 0xe7, 0x25]]];

  function viridis(t) {
    var range = scheme() === "light" ? [0, 0.8] : [0.2, 1];
    t = range[0] + Math.min(1, Math.max(0, t)) * (range[1] - range[0]);
    for (var i = 1; i < VIRIDIS.length; i += 1) {
      if (t <= VIRIDIS[i][0]) {
        var a = VIRIDIS[i - 1], b = VIRIDIS[i], f = (t - a[0]) / (b[0] - a[0]);
        return "rgb(" + a[1].map(function (lo, k) { return Math.round(lo + f * (b[1][k] - lo)); }).join(",") + ")";
      }
    }
    return "rgb(253,231,37)";
  }

  var FROM_SAID = { "latest frame": "the latest frame its run wrote",
                    structure: "the structure it was set up from",
                    start: "the structure it is to start from" };

  function moleculeSaid(m) {
    return "Its backbone (" + (m.residues || 0) + " residues, " + (FROM_SAID[m.from] || "its structure") +
      "), coloured from the N terminus (purple) to the C terminus (" +
      (scheme() === "light" ? "green" : "yellow") + "), the nearer part stronger";
  }

  function moleculeSvg(m, width, height, pad) {
    var plot = svg("svg", { viewBox: "0 0 " + width + " " + height, width: width, height: height,
                            class: "study-molecule", role: "img", "aria-label": moleculeSaid(m) }, null);
    var aspect = m.aspect > 0 ? m.aspect : 1;
    var scale = Math.min((width - 2 * pad) / aspect, height - 2 * pad);
    var ox = (width - scale * aspect) / 2, oy = (height - scale) / 2;
    var thick = Math.max(0.9, Math.min(3.2, Math.min(width, height) / 13 *
      Math.sqrt(40 / Math.max(40, m.residues || 40))));
    var segments = [];
    (m.pieces || []).forEach(function (piece) {
      for (var k = 0; k + 1 < piece.length; k += 1) {
        segments.push([(piece[k][2] + piece[k + 1][2]) / 2, piece[k], piece[k + 1]]);
      }
    });
    segments.sort(function (a, b) { return a[0] - b[0]; });
    segments.forEach(function (s) {
      svg("line", { x1: (ox + s[1][0] * scale).toFixed(1), y1: (oy + s[1][1] * scale).toFixed(1),
                    x2: (ox + s[2][0] * scale).toFixed(1), y2: (oy + s[2][1] * scale).toFixed(1),
                    stroke: viridis((s[1][3] + s[2][3]) / 2), "stroke-opacity": (0.35 + 0.65 * s[0]).toFixed(2),
                    "stroke-width": (thick * (0.7 + 0.6 * s[0])).toFixed(2), "stroke-linecap": "round" }, plot);
    });
    return plot;
  }

  /* The molecule's frame: the picture, else the study's ID or name, with
   * why: nothing to draw, or its structure not read (linked out of the
   * study, or outside it), which is not a study with none. */
  function moleculeFrame(study, cls, width, height, pad) {
    var frame = make("div", cls);
    if (study.molecule && (study.molecule.pieces || []).length) {
      frame.appendChild(moleculeSvg(study.molecule, width, height, pad));
      frame.title = moleculeSaid(study.molecule);
      frame.dataset.picture = "molecule";
    } else {
      frame.appendChild(make("span", "study-thumb-none", shownId(study) || study.name));
      frame.dataset.picture = "none";
      frame.title = unread(study, "molecule")
        ? "Not drawn: its structure was not read. " + (study.not_read || "")
        : study.molecule_outside
          ? "Not drawn: the structure it is to start from, " + study.molecule_outside +
            ", is outside its folder, and only the study's own files are read"
          : "Nothing to draw yet: no structure of it is here";
    }
    return frame;
  }

  /* Whether a record not read leaves this part of the study unknown
   * (workspace.py `card_of`: state, length, means, molecule). */
  function unread(study, part) {
    return (study.unread || []).indexOf(part) >= 0;
  }

  /* The records of a study not read, named, the why on hover. */
  function notReadShown(study, cls) {
    var said = make("span", cls, "Not read: " + (study.unavailable || []).join(", "));
    said.title = study.not_read || "";
    if (!(study.unavailable || []).length) said.textContent = study.not_read || "Not read";
    return said;
  }

  /* Drawn again in the new scheme's colours. */
  document.addEventListener("fmx:theme", function () { if (studies.length) draw(); });

  /* ---- What a row and a card say ------------------------------------ */
  /* A length of simulation in the unit it reads in: under a nanosecond in
   * picoseconds, as every plot of the page reads it. */
  function lengthSaid(ns) {
    if (ns == null || !isFinite(ns)) return "";
    if (ns === 0) return "0 ns";
    if (Math.abs(ns) < 1) return Number((ns * 1000).toPrecision(3)) + " ps";
    return Number(ns.toPrecision(4)) + " ns";
  }

  function durationSaid(seconds) {
    if (seconds == null || !isFinite(seconds)) return "";
    if (seconds < 90) return Math.max(1, Math.round(seconds)) + " s";
    if (seconds < 5400) return Math.round(seconds / 60) + " min";
    if (seconds < 129600) return (seconds / 3600).toFixed(1) + " h";
    return (seconds / 86400).toFixed(1) + " days";
  }

  /* A mean as a row says it: with its error where determined; a study of
   * replicas, the mean of their means and their spread. */
  function meanShort(m) {
    if (!m) return "";
    if (m.mean == null) return "not recorded";
    if (m.withheld && m.spread != null && m.spread > 0) {
      // To the spread's second figure, as a mean is to its error's: no
      // more figures than the replicas agree to.
      var places = Math.max(0, Math.min(6, 1 - Math.floor(Math.log10(m.spread))));
      return m.mean.toFixed(places) + (m.unit ? " " + m.unit : "") + ", spread " + m.spread.toFixed(places);
    }
    return withError(m.mean, m.withheld ? null : m.error, m.unit);
  }

  function meanItem(m) {
    var item = make("span", "studies-mean-item");
    item.dataset.determined = String(!m.withheld && m.error != null);
    var dot = make("span", "studies-mean-dot");
    dot.setAttribute("aria-hidden", "true");
    var said = meanShort(m);
    var value = make("span", "studies-mean-v mono");
    // The value and its unit kept together; a replicas' spread may go on
    // the next line where the column is narrow.
    var at = said.indexOf(", spread ");
    if (at < 0) {
      value.appendChild(make("span", "studies-mean-n", said));
    } else {
      value.append(make("span", "studies-mean-n", said.slice(0, at + 1)), " ",
                   make("span", "studies-mean-n", said.slice(at + 2)));
    }
    item.append(dot, make("span", "studies-mean-k", m.label || m.analysis), value);
    item.title = m.withheld ? "Not determined: " + firstSentence(m.withheld) + " (M for why)"
      : m.error != null ? "Determined, with its standard error" : "Recorded without an error";
    return item;
  }

  /* Every mean a study recorded, each with why it is not determined where
   * it is not: opened from +N, or with M on a row. */
  var popover = null;
  function closeMeans() {
    if (popover) { popover.remove(); popover = null; }
  }

  /* Why a mean is not yet determined: its first sentence, the rest
   * opened on request, so ten reasons do not fill the screen. */
  function sentenceEnd(text) { return text.search(/[.)]\s+(?=[A-Z])/); }
  function firstSentence(text) {
    var cut = sentenceEnd(text);
    return cut < 0 ? text : text.slice(0, cut + 1);
  }
  function whySaid(text) {
    var dd = make("dd", "studies-means-pop-why");
    var cut = sentenceEnd(text);
    if (cut < 0) { dd.textContent = text; return dd; }
    var details = make("details");
    details.appendChild(make("summary", "", text.slice(0, cut + 1)));
    details.appendChild(make("p", "studies-means-pop-rest", text.slice(cut + 1).trim()));
    dd.appendChild(details);
    return dd;
  }

  function openMeans(study, anchor) {
    var wasFor = popover && popover.dataset.path === study.path;
    closeMeans();
    if (wasFor) return;
    var box = make("div", "studies-means-pop");
    box.dataset.path = study.path;
    box.setAttribute("role", "dialog");
    box.setAttribute("aria-label", "Every mean " + study.name + " recorded");
    box.appendChild(make("div", "studies-means-pop-head", study.name + ": every mean recorded"));
    var list = make("dl", "studies-means-pop-list");
    (study.means || []).forEach(function (m) {
      var k = make("dt");
      k.appendChild(meanItem(m));
      list.appendChild(k);
      if (m.withheld) list.appendChild(whySaid(m.withheld));
    });
    box.appendChild(list);
    box.appendChild(make("div", "muted small", "A filled dot is a mean determined with its error; a ring, one not yet."));
    document.body.appendChild(box);
    var r = anchor.getBoundingClientRect();
    var w = box.offsetWidth, h = box.offsetHeight;
    var left = Math.max(8, Math.min(r.left, window.innerWidth - w - 8));
    var top = r.bottom + 6 + h > window.innerHeight ? Math.max(8, r.top - h - 6) : r.bottom + 6;
    box.style.left = left + "px";
    box.style.top = top + "px";
    popover = box;
    popoverAt = { anchor: anchor, top: r.top, left: r.left };
  }

  /* Closed once the page moves under it: not by a scroll that put its
   * button in view just before it opened (Tab, or a click that scrolled
   * first), which closed it as it opened. */
  var popoverAt = null;
  function movedUnder() {
    if (!popover || !popoverAt) return;
    var r = popoverAt.anchor.getBoundingClientRect();
    if (Math.abs(r.top - popoverAt.top) > 1 || Math.abs(r.left - popoverAt.left) > 1) closeMeans();
  }

  document.addEventListener("click", function (event) {
    if (popover && !popover.contains(event.target) && !event.target.closest(".studies-means-more")) closeMeans();
  });
  document.addEventListener("scroll", movedUnder, true);

  function meansShown(study, most) {
    var host = make("div", "studies-means");
    var means = study.means || [];
    means.slice(0, most).forEach(function (m) { host.appendChild(meanItem(m)); });
    // The rest, or with none more why those shown are not determined: a
    // button, so the why is reached from the keyboard, not on hover only.
    var why = means.slice(0, most).some(function (m) { return !!m.withheld; });
    if (means.length > most || why) {
      var more = make("button", "studies-means-more",
        means.length > most ? "+" + (means.length - most) : "why");
      more.type = "button";
      more.title = (means.length > most ? "Every mean " + study.name + " recorded"
        : "Why it is not determined") + " (M)";
      more.setAttribute("aria-label", means.length > most
        ? "Every mean recorded, " + means.length + " in all" : "Why it is not determined");
      more.addEventListener("click", function () { openMeans(study, more); });
      host.appendChild(more);
    }
    // A mean's record not read: said, never taken for one not recorded.
    if (unread(study, "means")) host.appendChild(notReadMark(study));
    return host;
  }

  function notReadMark(study) {
    var mark = make("span", "studies-means-unread muted small", "not read");
    mark.title = study.not_read || "A record of its means was not read";
    return mark;
  }

  /* Where it stands: its light and word, and for a run going how far. */
  function stateShown(study) {
    var host = make("div", "studies-state-said");
    var line = make("span", "studies-state-line");
    var light = make("span", "recent-light");
    light.dataset.state = String(study.state || "").toLowerCase();
    light.setAttribute("aria-hidden", "true");
    var word = study.state || "";
    if (study.state === "running" && study.progress != null) word += ", " + Math.round(study.progress * 100) + "%";
    line.append(light, document.createTextNode(" " + word));
    host.appendChild(line);
    if (study.state === "running" && study.progress != null) {
      var bar = make("span", "studies-progress");
      var fill = make("i");
      fill.style.width = Math.round(study.progress * 100) + "%";
      bar.appendChild(fill);
      bar.setAttribute("role", "progressbar");
      bar.setAttribute("aria-valuenow", String(Math.round(study.progress * 100)));
      bar.setAttribute("aria-valuemin", "0");
      bar.setAttribute("aria-valuemax", "100");
      host.appendChild(bar);
      var said = [];
      if (study.planned_ns) said.push(lengthSaid(study.simulated_ns || 0) + " of " + lengthSaid(study.planned_ns));
      if (study.seconds_left != null) said.push(durationSaid(study.seconds_left) + " left");
      if (said.length) host.appendChild(make("span", "studies-state-more muted small", said.join(" · ")));
    }
    return host;
  }

  /* How long it simulated, and for runs side by side, each one's. */
  function lengthShown(study) {
    var host = make("div", "studies-length");
    if (study.simulated_ns == null && unread(study, "length")) {
      host.textContent = "not read";
      host.title = study.not_read || "A record of how long it ran was not read";
      return host;
    }
    if (study.simulated_ns == null) {
      host.textContent = "none";
      host.title = study.state === "analysed" ? "A trajectory analysed: its length is its frames'"
        : "Nothing simulated yet";
      return host;
    }
    host.appendChild(make("span", "mono", lengthSaid(study.simulated_ns)));
    if (study.runs > 1 && study.planned_ns) {
      host.appendChild(make("span", "studies-length-each muted small",
        study.runs + " × " + lengthSaid(study.planned_ns / study.runs)));
    }
    return host;
  }

  /* Its folder, from where the studies were looked for: what tells two of
   * one name apart. */
  function folderSaid(study) {
    var path = String(study.path || "");
    if (lookedIn && path.indexOf(lookedIn) === 0) path = path.slice(lookedIn.length).replace(/^[\\/]/, "");
    return path || study.path;
  }

  function clashes(study) {
    return studies.filter(function (s) { return s.name === study.name; }).length > 1;
  }

  function openButton(study, page) {
    if (beyondLoopback()) return null;
    var button = ICONS.button(page === "files" ? "files" : "open",
      page === "files" ? "Its files" : "Open", page === "files" ? "study-files" : "study-open");
    button.addEventListener("click", function () { openStudy(study.path, button, page); });
    return button;
  }

  function card(study) {
    var item = make("article", "study-card");
    item.setAttribute("role", "listitem");
    item.dataset.path = study.path;
    item.dataset.state = study.state;
    item.appendChild(moleculeFrame(study, "study-thumb", 260, 130, 12));
    var body = make("div", "study-body");
    var head = make("div", "study-head");
    head.appendChild(make("span", "study-name", study.name));
    head.appendChild(make("span", "study-state", study.state === "running" && study.progress != null
      ? "running, " + Math.round(study.progress * 100) + "%" : study.state));
    body.appendChild(head);
    var what = [shownId(study), study.kind, lengthSaid(study.simulated_ns), timeSaid(study.last_active || study.when)];
    body.appendChild(make("div", "study-what", what.filter(Boolean).join(" · ")));
    if (clashes(study)) body.appendChild(make("div", "study-folder-said muted small mono", folderSaid(study)));
    if (study.not_read) body.appendChild(notReadShown(study, "study-not-read small"));
    if (study.state === "running" && study.progress != null) {
      var bar = make("span", "studies-progress");
      var fill = make("i");
      fill.style.width = Math.round(study.progress * 100) + "%";
      bar.appendChild(fill);
      body.appendChild(bar);
    }
    if ((study.means || []).length) body.appendChild(meansShown(study, 2));
    else if (unread(study, "means")) body.appendChild(notReadMark(study));
    if (study.free_energy) {
      body.appendChild(make("div", "study-free-energy",
        study.free_energy.refused ? "Free energy refused: " + study.free_energy.refused
          : "Free energy recombined."));
    }
    if (editing === study.path) {
      body.appendChild(tagEditor(study));
    } else {
      if ((study.tags || []).length) body.appendChild(tagsShown(study));
      if (study.note) body.appendChild(make("div", "study-note", study.note));
    }
    var actions = make("div", "study-actions");
    var pick = make("label", "study-pick");
    var box = make("input");
    box.type = "checkbox";
    box.checked = chosen.indexOf(study.path) >= 0;
    box.dataset.path = study.path;
    box.addEventListener("click", function (event) { picked(study.path, box.checked, event); });
    pick.append(box, document.createTextNode(" Choose"));
    actions.append.apply(actions, [openButton(study), tagButton(study), openButton(study, "files"), pick]
      .filter(Boolean));
    body.appendChild(actions);
    item.appendChild(body);
    return item;
  }

  /* ---- As a list ---------------------------------------------------- */
  /* The studies as rows first (user, 10-18: the table first), each with
   * its molecule, where it stands, how long it ran, every mean and when
   * it was last active; cards as a gallery where chosen. Kept per person
   * in this browser, as the views saved, the grouping, the density and the
   * columns shown. */
  var VIEW_KEPT = "fmx.studiesView";
  var LAYOUT_KEPT = "fmx.studiesLayout";
  var VIEWS_KEPT = "fmx.studiesViews";
  var view = "table";
  var sortBy = { key: "active", down: true };
  var meanSort = "";
  var STATE_ORDER = { running: 0, failed: 1, stopped: 2, interrupted: 3, incomplete: 4,
                      completed: 5, analysed: 6, "not started": 7, "not known": 8 };
  var layout = { group: "study", density: "comfortable", hidden: [] };
  var filters = { states: [], systems: [] };
  var savedViews = [];
  var activeView = "all";
  var expanded = {};
  var runsOf = {};
  var runsNotRead = {};
  var COLUMNS = [["molecule", "Molecule"], ["system", "System"], ["length", "Simulated"],
                 ["means", "Means"], ["active", "Active"]];
  var BUILT_IN = [
    { id: "all", name: "All", test: function () { return true; } },
    { id: "running", name: "Running", test: function (s) { return s.state === "running"; } },
    { id: "sampling", name: "Needs more sampling",
      test: function (s) { return (s.means || []).some(function (m) { return !!m.withheld; }); } },
    { id: "not-started", name: "Not started", test: function (s) { return s.state === "not started"; } },
    // Studies the person archived, left out of every other view; shown
    // only while there are some.
    { id: "archived", name: "Archived", test: function (s) { return !!s.archived; } },
  ];

  function readKept(key, fallback) {
    try {
      var kept = JSON.parse(localStorage.getItem(key) || "null");
      return kept == null ? fallback : kept;
    } catch (e) { return fallback; }
  }

  function keep(key, value) {
    try { localStorage.setItem(key, JSON.stringify(value)); } catch (e) { /* not kept */ }
  }

  function keptView() {
    // A choice kept as words before this page kept JSON.
    try {
      var raw = localStorage.getItem(VIEW_KEPT);
      if (raw === "cards" || raw === "\"cards\"") return "cards";
    } catch (e) { /* not kept */ }
    return "table";
  }

  function showAs(chosenView, remember) {
    view = chosenView === "cards" ? "cards" : "table";
    Array.prototype.forEach.call(document.querySelectorAll("[data-studies-view]"), function (b) {
      var on = b.getAttribute("data-studies-view") === view;
      b.classList.toggle("active", on);
      b.setAttribute("aria-pressed", String(on));
    });
    if (remember) {
      try { localStorage.setItem(VIEW_KEPT, view); } catch (e) { /* not kept */ }
    }
    draw();
  }

  function shownColumn(key) { return layout.hidden.indexOf(key) < 0; }

  function sortValue(study, key) {
    if (key === "name") return String(study.name || "").toLowerCase();
    if (key === "system") return String(shownId(study) || "").toLowerCase();
    if (key === "state") return STATE_ORDER[String(study.state || "").toLowerCase()];
    if (key === "length") return study.simulated_ns == null ? null : study.simulated_ns;
    if (key.indexOf("mean:") === 0) {
      var m = (study.means || []).find(function (x) { return meanName(x) === key.slice(5); });
      return m && typeof m.mean === "number" && isFinite(m.mean) ? m.mean : null;
    }
    if (key === "when" || key === "active") {
      var t = new Date(key === "when" ? study.when : study.last_active || study.when).getTime();
      return isNaN(t) ? null : t;
    }
    return null;
  }

  /* The observables the shown studies have a mean of, the one most of
   * them have first: what the Means column is sorted by. */
  function meanName(m) { return m.label || m.analysis; }
  function meanNames(shown) {
    var count = {};
    var order = [];
    shown.forEach(function (study) {
      (study.means || []).forEach(function (m) {
        var name = meanName(m);
        if (!(name in count)) { count[name] = 0; order.push(name); }
        count[name] += 1;
      });
    });
    return order.sort(function (a, b) { return count[b] - count[a]; });
  }

  function sorted(shown) {
    var key = sortBy.key;
    return shown.slice().sort(function (a, b) {
      var x = sortValue(a, key), y = sortValue(b, key);
      // What has no value goes last, either way.
      if (x == null && y == null) return 0;
      if (x == null) return 1;
      if (y == null) return -1;
      var order = x < y ? -1 : x > y ? 1 : 0;
      return sortBy.down ? -order : order;
    });
  }

  function refocus(selector, key, value) {
    var found = Array.prototype.find.call(document.querySelectorAll(selector), function (node) {
      return node.dataset[key] === value;
    });
    if (found) found.focus();
  }

  /* The groups the rows fall in, as Group by asks: none (one group), the
   * system, a tag, the state or the folder they sit in. By study, replicas
   * fold under their study and nothing else is grouped, but the studies
   * where nothing ran yet, gathered at the foot. */
  function grouped(shown) {
    var by = layout.group;
    if (by === "none") return [{ name: "", rows: shown }];
    if (by === "study") {
      var begun = shown.filter(function (s) { return s.state !== "not started"; });
      var waiting = shown.filter(function (s) { return s.state === "not started"; });
      var out = [{ name: "", rows: begun }];
      if (waiting.length) out.push({ name: "Not started", said: "folders a study was begun in and nothing ran", rows: waiting });
      return out;
    }
    var keyed = {};
    var order = [];
    shown.forEach(function (s) {
      var keys = by === "system" ? [shownId(s) || "No system"]
        : by === "state" ? [s.state || "unknown"]
        : by === "folder" ? [folderSaid(s).split(/[\\/]/).slice(0, -1).join("/") || "Here"]
        : (s.tags || []).length ? s.tags : ["No tag"];
      keys.forEach(function (k) {
        if (!keyed[k]) { keyed[k] = []; order.push(k); }
        keyed[k].push(s);
      });
    });
    if (by === "state") order.sort(function (a, b) { return (STATE_ORDER[a] || 9) - (STATE_ORDER[b] || 9); });
    else order.sort(function (a, b) { return a.localeCompare(b); });
    return order.map(function (k) { return { name: k, rows: keyed[k] }; });
  }

  function cell(cls, child) {
    var td = make("td", cls);
    if (child) td.appendChild(child);
    return td;
  }

  /* One study as a row; a run of a study of several, under it, as one. */
  function row(study, child) {
    var tr = make("tr", "studies-row" + (child ? " studies-child" : ""));
    tr.dataset.path = study.path;
    tr.dataset.state = String(study.state || "").toLowerCase();
    var pick = make("td", "studies-pick");
    var box = make("input");
    box.type = "checkbox";
    box.checked = chosen.indexOf(study.path) >= 0;
    box.setAttribute("aria-label", "Compare " + study.name);
    box.dataset.path = study.path;
    box.addEventListener("click", function (event) { picked(study.path, box.checked, event); });
    pick.appendChild(box);
    tr.appendChild(pick);
    if (shownColumn("molecule")) tr.appendChild(cell("studies-mol", moleculeFrame(study, "studies-mol-frame", 64, 34, 3)));
    var name = make("td", "studies-name");
    var line = make("div", "studies-name-line");
    if (!child && study.runs > 1) {
      var fold = make("button", "studies-fold");
      fold.type = "button";
      var open = !!expanded[study.path];
      fold.setAttribute("aria-expanded", String(open));
      fold.setAttribute("aria-label", (open ? "Fold " : "Show ") + "its " + study.runs + " runs");
      fold.title = (open ? "Fold" : "Show") + " its " + study.runs + " runs";
      fold.textContent = open ? "▾" : "▸";
      fold.addEventListener("click", function () { toggleRuns(study); });
      line.appendChild(fold);
    }
    // The name opens the study, as Open does at the row's end.
    var named = make(beyondLoopback() ? "span" : "button", "studies-name-text",
                     child ? runName(study) : study.name);
    if (beyondLoopback()) {
      // Reached by j and k, as a row is, but not opened here.
      named.tabIndex = 0;
      named.title = study.path;
    } else {
      named.type = "button";
      named.title = "Open " + study.path;
      named.addEventListener("click", function () { openStudy(study.path, named); });
    }
    line.appendChild(named);
    (study.tags || []).forEach(function (t) { line.appendChild(make("span", "study-tag studies-row-tag", t)); });
    name.appendChild(line);
    var sub = [];
    if (!child && clashes(study)) sub.push(folderSaid(study));
    if (!child) sub.push([study.kind, study.forcefield].filter(Boolean).join(" · "));
    if (study.note) sub.push(study.note);
    sub = sub.filter(Boolean);
    if (sub.length) name.appendChild(make("span", "studies-kind muted small", sub.join(" · ")));
    if (study.not_read) name.appendChild(notReadShown(study, "studies-not-read small"));
    tr.appendChild(name);
    tr.appendChild(cell("studies-state", stateShown(study)));
    if (shownColumn("system")) {
      var system = make("td", "studies-system mono", shownId(study));
      if (study.structure && study.structure !== system.textContent) system.title = study.structure;
      tr.appendChild(system);
    }
    if (shownColumn("length")) tr.appendChild(cell("studies-len", lengthShown(study)));
    if (shownColumn("means")) {
      tr.appendChild(cell("studies-means-cell", (study.means || []).length ? meansShown(study, 2)
        : unread(study, "means") ? notReadMark(study)
        : make("span", "muted small", study.state === "not started" ? "none yet" : "no mean recorded")));
    }
    if (shownColumn("active")) {
      var active = make("td", "studies-active muted", timeSaid(study.last_active || study.when));
      if (study.last_active) active.title = "Last active " + new Date(study.last_active).toLocaleString() +
        (study.when ? "; begun " + new Date(study.when).toLocaleString() : "");
      tr.appendChild(active);
    }
    var acts = make("div", "studies-actions");
    acts.append.apply(acts, [openButton(study), tagButton(study), openButton(study, "files")]
      .filter(Boolean));
    tr.appendChild(cell("studies-acts", acts));
    return tr;
  }

  function runName(run) {
    var name = String(run.name || "");
    var seed = /random[-_]seed[-_](\d+)/i.exec(name);
    return seed ? "seed " + seed[1] : name;
  }

  /* A study of several runs unfolded: its runs, listed as All studies
   * lists studies (each run is one), the first time from the server. */
  function toggleRuns(study) {
    expanded[study.path] = !expanded[study.path];
    if (expanded[study.path] && !runsOf[study.path]) {
      runsOf[study.path] = "loading";
      var folder = study.path.replace(/[\\/]+$/, "") + "/runs";
      fetch("/api/studies?" + new URLSearchParams({ path: folder }))
        .then(function (r) { return r.json(); })
        .then(function (d) {
          runsOf[study.path] = d && d.ok ? (d.studies || []).sort(function (a, b) {
            return String(a.name).localeCompare(String(b.name), undefined, { numeric: true });
          }) : [];
          // Its runs not read (its `runs` a link out of it): said so.
          runsNotRead[study.path] = (d && d.not_read) || "";
          draw();
        })
        .catch(function () { runsOf[study.path] = []; draw(); });
    }
    draw();
    refocus(".studies-row .studies-name-text", "path", study.path);
  }

  /* Which mean the Means column is sorted by: chosen beside its head,
   * smallest first. */
  function meanSortChooser(names, now) {
    var pick = make("select", "studies-means-sort");
    pick.setAttribute("aria-label", "The mean the Means column is sorted by");
    names.forEach(function (name) {
      var option = make("option", "", name);
      option.value = name;
      option.selected = name === now;
      pick.appendChild(option);
    });
    pick.addEventListener("change", function () {
      meanSort = pick.value;
      sortBy = { key: "mean:" + pick.value, down: false };
      draw();
      var again = document.querySelector("#studies-table .studies-means-sort");
      if (again) again.focus();
    });
    return pick;
  }

  function table(shown) {
    var host = el("studies-table");
    host.classList.toggle("is-compact", layout.density === "compact");
    var heads = [["", null, "pick"]];
    if (shownColumn("molecule")) heads.push(["", null, "mol"]);
    heads.push(["Study", "name"], ["State", "state"]);
    if (shownColumn("system")) heads.push(["System", "system", "system"]);
    if (shownColumn("length")) heads.push(["Simulated", "length", "length"]);
    var names = meanNames(shown);
    if (meanSort && names.indexOf(meanSort) < 0) meanSort = "";
    var meanKey = "mean:" + (meanSort || names[0] || "");
    if (shownColumn("means")) heads.push(["Means", names.length ? meanKey : null, "means"]);
    if (shownColumn("active")) heads.push(["Active", "active", "active"]);
    heads.push(["", null, "acts"]);
    var head = make("thead");
    var tr = make("tr");
    heads.forEach(function (h) {
      var th = make("th", h[2] ? "studies-th-" + h[2] : "");
      th.scope = "col";
      if (h[1]) {
        var button = make("button", "studies-sort", h[0]);
        button.type = "button";
        button.dataset.key = h[1];
        var on = sortBy.key === h[1];
        th.setAttribute("aria-sort", on ? (sortBy.down ? "descending" : "ascending") : "none");
        if (on) button.dataset.dir = sortBy.down ? "down" : "up";
        button.addEventListener("click", function () {
          // A date, and a length, is sorted largest first to begin with.
          var first = h[1] === "when" || h[1] === "active" || h[1] === "length";
          sortBy = { key: h[1], down: sortBy.key === h[1] ? !sortBy.down : first };
          draw();
          // The table is built again: the keyboard stays where it was.
          refocus(".studies-sort", "key", h[1]);
        });
        if (h[2] === "means") {
          button.title = "Sort by the mean of " + meanKey.slice(5);
          th.appendChild(button);
          th.appendChild(meanSortChooser(names, meanKey.slice(5)));
        } else {
          th.appendChild(button);
        }
        if (h[2] === "means") th.title = "A filled dot is a mean determined with its error; a ring, one not yet";
      } else if (h[2] === "means") {
        th.textContent = h[0];
        th.title = "A filled dot is a mean determined with its error; a ring, one not yet";
      } else if (h[0]) {
        th.textContent = h[0];
      } else if (h[2] === "pick") {
        // Every study shown chosen, or let go.
        var all = make("input", "studies-pick-all");
        all.type = "checkbox";
        all.title = "Choose every study shown";
        all.setAttribute("aria-label", "Choose every study shown");
        all.addEventListener("click", function () { chooseEveryShown(all.checked); });
        th.appendChild(all);
      } else {
        th.appendChild(make("span", "sr-only", h[2] === "mol" ? "Molecule" : "Open"));
      }
      tr.appendChild(th);
    });
    head.appendChild(tr);
    var body = make("tbody");
    grouped(sorted(shown)).forEach(function (group) {
      if (group.name) {
        var g = make("tr", "studies-group");
        var gc = make("td");
        gc.colSpan = heads.length;
        gc.appendChild(make("b", "", group.name));
        gc.appendChild(document.createTextNode(" " + group.rows.length + (group.said ? " · " + group.said : "")));
        g.appendChild(gc);
        body.appendChild(g);
      }
      group.rows.forEach(function (study) {
        body.appendChild(row(study));
        if (editing === study.path) body.appendChild(editingRow(study, heads.length));
        if (expanded[study.path]) {
          var runs = runsOf[study.path];
          if (runs === "loading" || (runs && !runs.length)) {
            var wait = make("tr", "studies-child-wait");
            var wc = make("td", "muted small", runs === "loading" ? "Reading its runs…"
              : runsNotRead[study.path] ? "Not read: " + runsNotRead[study.path]
              : "None of its runs is here to read: its runs folder holds no study.");
            wc.colSpan = heads.length;
            wait.appendChild(wc);
            body.appendChild(wait);
          } else {
            (runs || []).forEach(function (run) {
              body.appendChild(row(run, true));
              if (editing === run.path) body.appendChild(editingRow(run, heads.length));
            });
          }
        }
      });
    });
    host.replaceChildren(head, body);
  }

  function editingRow(study, span) {
    var under = make("tr", "studies-editing");
    var c = make("td");
    c.colSpan = span;
    c.appendChild(tagEditor(study));
    under.appendChild(c);
    return under;
  }

  /* ---- Views, filters and how the list is laid out ------------------- */
  function viewNamed(id) {
    return BUILT_IN.find(function (v) { return v.id === id; }) ||
      savedViews.find(function (v) { return "saved:" + v.name === id; }) || null;
  }

  function inView(study, v) {
    if (!v || v.id === "all") return true;
    if (v.test) return v.test(study);
    return true;
  }

  function filtered() {
    var wanted = (el("studies-search").value || "").trim().toLowerCase().split(/\s+/).filter(Boolean);
    var v = viewNamed(activeView);
    return studies.filter(function (study) {
      var text = searchText(study);
      return (activeView === "archived" ? !!study.archived : !study.archived)
        && inView(study, v)
        && (!tagged || hasTag(study, tagged))
        && (!filters.states.length || filters.states.indexOf(study.state) >= 0)
        && (!filters.systems.length || filters.systems.indexOf(shownId(study) || "") >= 0)
        && wanted.every(function (word) { return text.indexOf(word) >= 0; });
    });
  }

  function drawViews() {
    var host = el("studies-views");
    if (!host) return;
    host.replaceChildren();
    var all = BUILT_IN.map(function (v) { return { id: v.id, name: v.name, test: v.test }; })
      .concat(savedViews.map(function (v) { return { id: "saved:" + v.name, name: v.name, saved: v }; }));
    all.forEach(function (v) {
      var b = make("button", "studies-view-tab" + (v.saved ? " is-saved" : ""));
      b.type = "button";
      b.setAttribute("role", "tab");
      b.dataset.view = v.id;
      b.setAttribute("aria-selected", String(activeView === v.id));
      b.appendChild(make("span", "", v.name));
      var kept = studies.filter(function (s) { return v.id === "archived" ? s.archived : !s.archived; });
      var n = v.saved ? kept.filter(function (s) { return matchesSaved(s, v.saved); }).length
        : kept.filter(function (s) { return v.test(s); }).length;
      if (v.id === "archived" && !n && activeView !== "archived") return;
      b.appendChild(make("span", "studies-view-n", String(n)));
      b.addEventListener("click", function () { chooseView(v.id); });
      if (v.saved) {
        b.title = "A view of your own; right-click or Delete to remove it";
        b.addEventListener("contextmenu", function (event) { event.preventDefault(); dropView(v.saved.name); });
        b.addEventListener("keydown", function (event) {
          if (event.key === "Delete" || event.key === "Backspace") { event.preventDefault(); dropView(v.saved.name); }
        });
      }
      host.appendChild(b);
    });
    if (naming) {
      var input = make("input", "studies-view-name");
      input.type = "text";
      input.maxLength = 40;
      input.placeholder = "A name for this view";
      input.setAttribute("aria-label", "A name for this view: its search, filters and grouping");
      input.addEventListener("keydown", function (event) {
        if (event.key === "Enter") { event.preventDefault(); keepView(input.value); }
        else if (event.key === "Escape") { event.preventDefault(); naming = false; drawViews(); }
      });
      input.addEventListener("blur", function () { if (naming) { naming = false; drawViews(); } });
      host.appendChild(input);
      return;
    }
    var save = make("button", "studies-view-save", "+ Save this view");
    save.type = "button";
    save.title = "Keep the search, the filters and the grouping under a name, in this browser";
    save.addEventListener("click", saveView);
    host.appendChild(save);
  }

  function matchesSaved(study, saved) {
    var text = searchText(study);
    var words = String(saved.search || "").toLowerCase().split(/\s+/).filter(Boolean);
    return (!saved.tag || hasTag(study, saved.tag))
      && (!(saved.states || []).length || saved.states.indexOf(study.state) >= 0)
      && (!(saved.systems || []).length || saved.systems.indexOf(shownId(study) || "") >= 0)
      && words.every(function (w) { return text.indexOf(w) >= 0; });
  }

  function chooseView(id) {
    activeView = id;
    var v = viewNamed(id);
    if (v && v.search !== undefined) {
      // A view of one's own puts back what it was saved with.
      el("studies-search").value = v.search || "";
      tagged = v.tag || null;
      filters = { states: (v.states || []).slice(), systems: (v.systems || []).slice() };
      if (v.group) layout.group = v.group;
      keep(LAYOUT_KEPT, layout);
    } else {
      filters = { states: [], systems: [] };
      tagged = null;
      el("studies-search").value = "";
    }
    draw();
    refocus(".studies-view-tab", "view", id);
  }

  /* Saved under a name typed in place of the button: Enter keeps it,
   * Escape lets it go. */
  var naming = false;
  function saveView() {
    naming = true;
    drawViews();
    var input = document.querySelector(".studies-view-name");
    if (input) input.focus();
  }

  function keepView(name) {
    naming = false;
    name = (name || "").replace(/\s+/g, " ").trim().slice(0, 40);
    if (!name || BUILT_IN.some(function (v) { return v.name.toLowerCase() === name.toLowerCase(); })) {
      drawViews();
      return;
    }
    savedViews = savedViews.filter(function (v) { return v.name !== name; });
    savedViews.push({ name: name, search: el("studies-search").value || "", tag: tagged || "",
                      states: filters.states.slice(), systems: filters.systems.slice(), group: layout.group });
    keep(VIEWS_KEPT, savedViews);
    activeView = "saved:" + name;
    draw();
    refocus(".studies-view-tab", "view", activeView);
  }

  function dropView(name) {
    savedViews = savedViews.filter(function (v) { return v.name !== name; });
    keep(VIEWS_KEPT, savedViews);
    if (activeView === "saved:" + name) activeView = "all";
    draw();
  }

  /* A filter as a chip that opens on its choices, each ticked to keep. */
  function chooser(label, key, values) {
    var box = make("details", "studies-chip");
    var chosenNow = filters[key];
    box.classList.toggle("is-on", chosenNow.length > 0);
    var summary = make("summary", "", chosenNow.length ? label + ": " + chosenNow.join(", ") : label);
    box.appendChild(summary);
    var list = make("div", "studies-chip-list");
    values.forEach(function (value) {
      var item = make("label", "studies-chip-item");
      var input = make("input");
      input.type = "checkbox";
      input.checked = chosenNow.indexOf(value) >= 0;
      input.addEventListener("change", function () {
        filters[key] = filters[key].filter(function (v) { return v !== value; });
        if (input.checked) filters[key].push(value);
        draw();
        var again = document.querySelector('.studies-chip[data-key="' + key + '"]');
        if (again) { again.open = true; var first = again.querySelector("input"); if (first) first.focus(); }
      });
      item.append(input, document.createTextNode(" " + (value || "none")));
      list.appendChild(item);
    });
    if (!values.length) list.appendChild(make("span", "muted small", "None here"));
    box.dataset.key = key;
    box.appendChild(list);
    return box;
  }

  function seg(label, options, chosenNow, onChoose) {
    var host = make("div", "seg studies-seg");
    host.setAttribute("role", "group");
    host.setAttribute("aria-label", label);
    options.forEach(function (o) {
      var b = make("button", "seg-btn" + (o[0] === chosenNow ? " active" : ""), o[1]);
      b.type = "button";
      b.setAttribute("aria-pressed", String(o[0] === chosenNow));
      b.addEventListener("click", function () { onChoose(o[0]); });
      host.appendChild(b);
    });
    return host;
  }

  function drawFilters() {
    var host = el("studies-filters");
    if (!host) return;
    host.replaceChildren();
    var states = [], systems = [];
    studies.forEach(function (s) {
      if (s.state && states.indexOf(s.state) < 0) states.push(s.state);
      var id = shownId(s) || "";
      if (systems.indexOf(id) < 0) systems.push(id);
    });
    states.sort(function (a, b) { return (STATE_ORDER[a] || 9) - (STATE_ORDER[b] || 9); });
    systems.sort();
    host.appendChild(chooser("State", "states", states));
    host.appendChild(chooser("System", "systems", systems));
    if (tagsUsed.length) {
      var tags = make("details", "studies-chip" + (tagged ? " is-on" : ""));
      tags.dataset.key = "tags";
      tags.appendChild(make("summary", "", tagged ? "Tag: " + tagged : "Tags"));
      var tl = make("div", "studies-chip-list");
      tagsUsed.forEach(function (used) {
        var b = make("button", "studies-chip-item studies-chip-pick", used.tag + " (" + used.studies + ")");
        b.type = "button";
        b.addEventListener("click", function () { narrowTo(tagged === used.tag ? null : used.tag); });
        tl.appendChild(b);
      });
      tags.appendChild(tl);
      host.appendChild(tags);
    }
    host.appendChild(make("span", "studies-filters-gap"));
    if (view === "table") {
      var g = make("label", "studies-group-by muted small", "Group by ");
      var select = make("select", "studies-group-select");
      [["study", "Study"], ["system", "System"], ["tag", "Tag"], ["state", "State"], ["folder", "Folder"],
       ["none", "Nothing"]].forEach(function (o) {
        var option = make("option", "", o[1]);
        option.value = o[0];
        option.selected = layout.group === o[0];
        select.appendChild(option);
      });
      select.addEventListener("change", function () {
        layout.group = select.value;
        keep(LAYOUT_KEPT, layout);
        draw();
        var again = document.querySelector(".studies-group-select");
        if (again) again.focus();
      });
      g.appendChild(select);
      host.appendChild(g);
      host.appendChild(seg("Density", [["comfortable", "Comfortable"], ["compact", "Compact"]], layout.density,
        function (d) { layout.density = d; keep(LAYOUT_KEPT, layout); draw(); }));
      var cols = make("details", "studies-chip studies-columns");
      cols.dataset.key = "columns";
      cols.appendChild(make("summary", "", "Columns"));
      var cl = make("div", "studies-chip-list");
      COLUMNS.forEach(function (c) {
        var item = make("label", "studies-chip-item");
        var input = make("input");
        input.type = "checkbox";
        input.checked = shownColumn(c[0]);
        input.addEventListener("change", function () {
          layout.hidden = layout.hidden.filter(function (k) { return k !== c[0]; });
          if (!input.checked) layout.hidden.push(c[0]);
          keep(LAYOUT_KEPT, layout);
          draw();
          var again = document.querySelector('.studies-chip[data-key="columns"]');
          if (again) again.open = true;
        });
        item.append(input, document.createTextNode(" " + c[1]));
        cl.appendChild(item);
      });
      cols.appendChild(cl);
      host.appendChild(cols);
    }
  }

  /* ---- Keys on the list --------------------------------------------- */
  /* / to search, J and K from row to row, Enter to open, X to choose for
   * Compare, C to compare the two chosen, M for every mean of the row. */
  function onStudiesPage() {
    var page = document.querySelector('.page[data-page="studies"]');
    return page && !page.hidden;
  }

  function rowsShown() {
    return Array.prototype.slice.call(document.querySelectorAll(
      view === "table" ? "#studies-table .studies-row" : "#studies-grid .study-card"));
  }

  function rowOf(node) {
    return node && node.closest ? node.closest(view === "table" ? ".studies-row" : ".study-card") : null;
  }

  function focusRow(r) {
    if (!r) return;
    var target = r.querySelector(".studies-name-text, .study-open");
    if (target) { target.focus(); r.scrollIntoView({ block: "nearest" }); }
  }

  document.addEventListener("keydown", function (event) {
    if (!onStudiesPage() || event.ctrlKey || event.metaKey || event.altKey || event.isComposing) return;
    var target = event.target;
    var typing = target && (target.isContentEditable || /^(INPUT|TEXTAREA|SELECT)$/.test(target.tagName)) &&
      !(target.type === "checkbox");
    if (event.key === "Escape" && popover) { closeMeans(); return; }
    if (typing) {
      if (event.key === "Escape" && target === el("studies-search") && target.value) {
        target.value = "";
        draw();
      }
      return;
    }
    if (document.querySelector(".agent-dialog:not([hidden])")) return;
    var rows = rowsShown();
    var here = rowOf(document.activeElement);
    var at = here ? rows.indexOf(here) : -1;
    if (event.key === "/") {
      event.preventDefault();
      el("studies-search").focus();
    } else if (event.key === "j" || event.key === "k") {
      if (!rows.length) return;
      event.preventDefault();
      var next = at < 0 ? 0 : Math.max(0, Math.min(rows.length - 1, at + (event.key === "j" ? 1 : -1)));
      focusRow(rows[next]);
    } else if (event.key === "x" && here) {
      var box = here.querySelector('input[type="checkbox"]');
      if (box) {
        event.preventDefault();
        box.checked = !box.checked;
        choose(box.dataset.path, box.checked, box);
        focusRow(here);
      }
    } else if (beyondLoopback() && (event.key === "c" || event.key === "t")) {
      return;
    } else if (event.key === "c" && chosen.length === 2) {
      event.preventDefault();
      compare();
    } else if (event.key === "t" && chosen.length) {
      event.preventDefault();
      openTagger();
    } else if (event.key === "Escape" && chosen.length && !document.querySelector(".studies-tagger")) {
      event.preventDefault();
      clearChosen(true);
    } else if (event.key === "m" && here) {
      var path = here.dataset.path;
      var study = studies.find(function (s) { return s.path === path; }) ||
        Object.keys(runsOf).reduce(function (found, k) {
          return found || (Array.isArray(runsOf[k]) ? runsOf[k].find(function (s) { return s.path === path; }) : null);
        }, null);
      if (study && (study.means || []).length) {
        event.preventDefault();
        openMeans(study, here.querySelector(".studies-means") || here);
      }
    }
  });

  /* Whether the table runs on under its pinned last column, so its edge
   * can say so: columns scrolled beneath Open and Tag were cut with nothing
   * to show there was more. */
  function moreToTheRight(wrap) {
    if (!wrap) return;
    var said = function () {
      wrap.dataset.more = wrap.scrollLeft + wrap.clientWidth < wrap.scrollWidth - 2 ? "right" : "";
    };
    if (!wrap.dataset.watched) {
      wrap.dataset.watched = "1";
      wrap.addEventListener("scroll", said, { passive: true });
      window.addEventListener("resize", said);
    }
    said();
  }

  /* Each study as two or three lines where the column cannot hold the
   * table (studies-table.css `.is-lines`): the side panel open, or the
   * Agent beside, on a laptop, and a phone. The table is measured as it
   * would be laid out, what it holds and the column's width together; it
   * ran past its edge, its means and buttons under a sideways scroll. */
  function fitTheList(wrap) {
    if (!wrap || wrap.hidden) return;
    wrap.classList.remove("is-lines");
    if (wrap.scrollWidth > wrap.clientWidth + 1) wrap.classList.add("is-lines");
    if (!wrap.dataset.fitted && window.ResizeObserver) {
      wrap.dataset.fitted = "1";
      var width = wrap.clientWidth;
      // Measured again as the column widens or narrows, not as rows are
      // unfolded, which changes only its height.
      new ResizeObserver(function () {
        if (wrap.hidden || wrap.clientWidth === width) return;
        width = wrap.clientWidth;
        fitTheList(wrap);
        moreToTheRight(wrap);
      }).observe(wrap);
    }
  }

  function draw() {
    var grid = el("studies-grid");
    var empty = el("studies-empty");
    if (!grid) return;
    closeMeans();
    var shown = filtered();
    var filter = el("studies-tag-filter");
    if (filter) {
      filter.hidden = !tagged;
      el("studies-tag-filter-said").textContent = tagged ? "Tagged " + tagged : "";
    }
    drawViews();
    drawFilters();
    var wrap = el("studies-table-wrap");
    var asTable = view === "table" && !!wrap;
    grid.hidden = asTable;
    if (wrap) wrap.hidden = !asTable || !shown.length;
    if (asTable) {
      grid.replaceChildren();
      table(shown);
      fitTheList(wrap);
      moreToTheRight(wrap);
    } else {
      grid.replaceChildren.apply(grid, shown.map(card));
    }
    var keysSaid = el("studies-keys");
    if (keysSaid) keysSaid.hidden = !shown.length;
    empty.hidden = shown.length > 0;
    empty.textContent = studies.length
      ? "No study matches that."
      : "No studies here. Look in another folder, or start one from the Agent or the Config page.";
    sayChosen();
  }

  /* The studies chosen, and one bar at the foot of the list while any is:
   * Compare (two), Tag, Files (one), Archive and Clear (user, 10-18: the
   * mock's bar "seems better" than Compare and Clear beside the pointer of
   * 10-09; "good ..... built it"). It is never over a row: the list keeps
   * room for it. Shift and a click choose every study between it and the
   * one chosen last; the box in the list's head, every study shown. */
  var lastPicked = null;

  function choose(path, on) {
    chosen = chosen.filter(function (p) { return p !== path; });
    if (on) chosen.push(path);
    lastPicked = path;
    sayChosen();
  }

  /* The boxes in the order the page shows them. */
  function boxesShown() {
    var selector = view === "table" ? "#studies-table tbody input[type=checkbox]"
      : ".study-card .study-pick input";
    return Array.prototype.filter.call(document.querySelectorAll(selector), function (b) {
      return !!b.dataset.path;
    });
  }

  /* A box clicked: with Shift, every study from the one chosen last to
   * this one is chosen (or let go) with it. */
  function picked(path, on, event) {
    var paths = boxesShown().map(function (b) { return b.dataset.path; });
    var from = lastPicked ? paths.indexOf(lastPicked) : -1;
    var to = paths.indexOf(path);
    if (event && event.shiftKey && from >= 0 && to >= 0) {
      paths.slice(Math.min(from, to), Math.max(from, to) + 1).forEach(function (p) {
        chosen = chosen.filter(function (q) { return q !== p; });
        if (on) chosen.push(p);
      });
      chosen = chosen.filter(function (q) { return q !== path; });
      if (on) chosen.push(path);
      lastPicked = path;
      sayChosen();
      return;
    }
    choose(path, on);
  }

  function chooseEveryShown(on) {
    var paths = boxesShown().map(function (b) { return b.dataset.path; });
    chosen = chosen.filter(function (p) { return paths.indexOf(p) < 0; });
    if (on) chosen = chosen.concat(paths);
    lastPicked = null;
    sayChosen();
  }

  function clearChosen(fromKeys) {
    var last = chosen[chosen.length - 1];
    var inBar = el("studies-bar").contains(document.activeElement);
    chosen = [];
    lastPicked = null;
    sayChosen();
    var heard = el("studies-chosen-heard");
    if (heard) heard.textContent = "Cleared.";
    // The keyboard was in the bar, now hidden: back to the study last
    // chosen, else the search.
    if (inBar || fromKeys) {
      var box = boxOf(last) || el("studies-search");
      if (box) box.focus({ preventScroll: !fromKeys });
    }
  }

  function boxOf(path) {
    return boxesShown().find(function (b) { return b.dataset.path === path; }) || null;
  }

  function nameOf(path) {
    var study = studyAt(path);
    return study ? study.name : "";
  }

  /* A study listed, or a run of one unfolded. */
  function studyAt(path) {
    var found = studies.find(function (s) { return s.path === path; });
    if (found) return found;
    var keys = Object.keys(runsOf);
    for (var k = 0; k < keys.length; k += 1) {
      var runs = runsOf[keys[k]];
      if (Array.isArray(runs)) {
        found = runs.find(function (s) { return s.path === path; });
        if (found) return found;
      }
    }
    return null;
  }

  /* The boxes, the head's box, the bar and what a screen reader hears, as
   * the studies chosen now are. */
  function sayChosen() {
    Array.prototype.forEach.call(document.querySelectorAll(
      ".study-card .study-pick input, #studies-table tbody input[type=checkbox]"), function (b) {
      if (b.dataset.path) b.checked = chosen.indexOf(b.dataset.path) >= 0;
      var holder = b.closest(".studies-row, .study-card");
      if (holder) holder.classList.toggle("is-chosen", b.checked);
    });
    var head = document.querySelector("#studies-table thead .studies-pick-all");
    if (head) {
      var paths = boxesShown().map(function (b) { return b.dataset.path; });
      var on = paths.filter(function (p) { return chosen.indexOf(p) >= 0; }).length;
      head.checked = paths.length > 0 && on === paths.length;
      head.indeterminate = on > 0 && on < paths.length;
    }
    var table = el("studies-table");
    if (table) table.classList.toggle("is-choosing", chosen.length > 0);
    drawBar();
    var heard = el("studies-chosen-heard");
    var n = chosen.length;
    var said = !n ? "" : n === 1 ? (nameOf(chosen[0]) || "One") + " chosen."
      : n === 2 ? chosen.map(nameOf).join(" and ") + " chosen. Compare is in the bar under the list."
      : n + " chosen.";
    if (heard && said && heard.textContent !== said) heard.textContent = said;
  }

  var barSaidFor = 0;
  function drawBar() {
    var bar = el("studies-bar");
    if (!bar) return;
    closeTagger();
    var n = chosen.length;
    bar.hidden = n === 0 && el("studies-bar-said").hidden;
    if (!n) {
      if (bar.contains(document.activeElement) && el("studies-bar-said").hidden) el("studies-search").focus();
      return;
    }
    el("studies-bar-said").hidden = true;
    el("studies-bar-in").hidden = false;
    bar.hidden = false;
    clearTimeout(barSaidFor);
    var count = el("studies-chosen-said");
    count.textContent = n + " chosen";
    count.title = chosen.map(nameOf).filter(Boolean).join(", ");
    var compare = el("studies-compare");
    compare.disabled = n !== 2;
    compare.title = n === 2 ? "Open " + chosen.map(nameOf).join(" and ") + " on the Compare page (C)"
      : "Choose two to compare: " + n + " chosen";
    var files = el("studies-files-chosen");
    files.disabled = n !== 1;
    files.title = n === 1 ? "Open its files" : "Choose one study to open its files";
    var archive = el("studies-archive");
    var undoing = activeView === "archived";
    archive.querySelector("span").textContent = undoing ? "Unarchive" : "Archive";
    archive.title = undoing ? "Show them in All studies again"
      : "Leave them out of All studies; nothing is moved or deleted, and Archived shows them";
    // A run of a study unfolded is not archived alone (study_tags.py says
    // the same): the study is.
    var aRun = chosen.some(function (p) { return !studies.some(function (s) { return s.path === p; }); });
    archive.disabled = aRun;
    if (aRun) archive.title = "A run of a study is archived with its study: archive the study.";
    el("studies-tag-chosen").title = "Add or take off a tag on every study chosen (T)";
  }

  /* What was done, said in the bar for a while, with Undo where it can be. */
  function barSays(text, undo) {
    var bar = el("studies-bar");
    var said = el("studies-bar-said");
    closeTagger();
    said.replaceChildren(make("span", "studies-bar-text", text));
    if (undo) {
      var back = make("button", "studies-bar-btn studies-bar-undo");
      back.type = "button";
      back.innerHTML = ICONS.svg("undo");
      back.appendChild(make("span", "", "Undo"));
      back.addEventListener("click", function () { said.hidden = true; undo(); });
      said.appendChild(back);
    }
    var close = make("button", "studies-bar-btn is-icon");
    close.type = "button";
    close.innerHTML = ICONS.svg("close");
    close.title = "Close";
    close.setAttribute("aria-label", "Close");
    close.addEventListener("click", function () { said.hidden = true; drawBar(); });
    said.appendChild(close);
    el("studies-bar-in").hidden = true;
    said.hidden = false;
    bar.hidden = false;
    clearTimeout(barSaidFor);
    // Gone after a while, but not from under the pointer or the keyboard.
    var wait = undo ? 8000 : 4000;
    var later = function () {
      if (bar.matches(":hover") || bar.contains(document.activeElement)) {
        barSaidFor = setTimeout(later, 1500);
        return;
      }
      said.hidden = true;
      drawBar();
    };
    barSaidFor = setTimeout(later, wait);
  }

  /* The two chosen, on the Compare page (compare.js). */
  function compare() {
    if (chosen.length !== 2 || !window.FastMDXCompare) return;
    window.FastMDXCompare.open(chosen[0], chosen[1]);
  }

  /* Tags of the studies chosen, changed one study at a time: what is put
   * on and what is taken off, applied by the server to the tags as kept
   * (POST /api/study-tags), so a tag kept since the list was read (an AI
   * app, another window) stays; its note and whether it is archived kept. */
  function keepTags(list) {
    return Promise.all(list.map(function (one) {
      return fetch("/api/study-tags", {
        method: "POST",
        headers: { "content-type": "application/json" },
        body: JSON.stringify({ path: one.path, add: one.add, remove: one.remove })
      }).then(function (r) { return r.json(); }).then(function (d) {
        if (!d || !d.ok) return { path: one.path, reason: (d && (d.reason || d.error)) || "not kept" };
        var study = studyAt(one.path);
        if (study) { study.tags = d.tags; study.note = d.note; }
        edited += 1;
        return { path: one.path };
      }).catch(function () { return { path: one.path, reason: "The server did not answer." }; });
    }));
  }

  function afterTags(results, done, undo) {
    var refused = results.filter(function (r) { return r.reason; });
    tagsUsed = countTags();
    offerTags();
    draw();
    if (refused.length) {
      barSays((nameOf(refused[0].path) || "A study") + ": " + refused[0].reason +
        (refused.length > 1 ? " (and " + (refused.length - 1) + " more)" : ""));
    } else {
      barSays(done, undo);
    }
  }

  /* Tag for many: each tag with how many of the chosen carry it. Ticked,
   * it goes on all; cleared, it comes off all; left alone, it stays as it
   * is on each, and a study's other tags are not touched. */
  function closeTagger() {
    var open = document.querySelector(".studies-tagger");
    if (open) open.remove();
  }

  function openTagger() {
    if (document.querySelector(".studies-tagger")) { closeTagger(); return; }
    var picked = chosen.map(studyAt).filter(Boolean);
    var n = picked.length;
    if (!n) return;
    var count = {};
    studies.forEach(function (s) { (s.tags || []).forEach(function (t) { if (!(t in count)) count[t] = 0; }); });
    picked.forEach(function (s) { (s.tags || []).forEach(function (t) { count[t] = (count[t] || 0) + 1; }); });
    var want = {};
    var shows = {};
    var box = make("div", "studies-tagger");
    box.setAttribute("role", "dialog");
    box.setAttribute("aria-label", "Tags of the " + n + " chosen");
    box.appendChild(make("div", "studies-tagger-head", "Tags of the " + n + " chosen"));
    var list = make("div", "studies-tagger-list");
    function item(tag) {
      var row = make("label", "studies-tagger-item");
      var tick = make("input");
      tick.type = "checkbox";
      var said = make("span", "studies-tagger-n");
      function show() {
        var k = want[tag] === "on" ? n : want[tag] === "off" ? 0 : count[tag];
        tick.checked = k === n;
        tick.indeterminate = k > 0 && k < n;
        said.textContent = k === 0 ? "none" : k === n ? (n === 1 ? "on it" : "all " + n) : k + " of " + n;
        row.classList.toggle("is-changed", want[tag] !== undefined);
      }
      tick.addEventListener("click", function (event) {
        event.preventDefault();
        // On some: on all, then off all, then as it was. On all: off all,
        // then as it was. On none: on all, then as it was.
        var now = want[tag];
        var next = count[tag] === n ? (now ? undefined : "off")
          : count[tag] === 0 ? (now ? undefined : "on")
          : now === undefined ? "on" : now === "on" ? "off" : undefined;
        if (next) want[tag] = next; else delete want[tag];
        show();
      });
      row.append(tick, make("span", "studies-tagger-tag", tag), said);
      shows[tag] = show;
      show();
      return row;
    }
    Object.keys(count).sort(function (a, b) { return a.localeCompare(b); })
      .forEach(function (t) { list.appendChild(item(t)); });
    box.appendChild(list);
    var adding = make("input", "studies-tagger-new");
    adding.type = "text";
    adding.maxLength = 40;
    adding.placeholder = "A new tag, then Enter";
    adding.setAttribute("aria-label", "A new tag for every study chosen");
    adding.setAttribute("list", "studies-tags-used");
    adding.addEventListener("keydown", function (event) {
      if (event.key === "Enter") {
        event.preventDefault();
        var tag = adding.value.replace(/\s+/g, " ").trim();
        if (!tag) { apply.click(); return; }
        var known = Object.keys(count).find(function (t) { return t.toLowerCase() === tag.toLowerCase(); });
        if (!known) { count[tag] = 0; list.appendChild(item(tag)); known = tag; }
        if (count[known] < n) want[known] = "on";
        shows[known]();
        adding.value = "";
      } else if (event.key === "Escape") {
        event.preventDefault();
        closeTagger();
        el("studies-tag-chosen").focus();
      }
    });
    box.appendChild(adding);
    box.appendChild(make("p", "studies-tagger-hint", "A tag ticked goes on all " + n +
      "; one cleared comes off all " + n + "; one left alone stays as it is on each."));
    var foot = make("div", "studies-tagger-foot");
    var cancel = make("button", "studies-bar-btn", "Cancel");
    cancel.type = "button";
    cancel.addEventListener("click", function () { closeTagger(); el("studies-tag-chosen").focus(); });
    var apply = make("button", "studies-bar-btn is-primary", "Apply");
    apply.type = "button";
    apply.addEventListener("click", function () {
      // Each study's change: the tags ticked it lacks, the tags cleared it
      // has. Undo is the change the other way.
      var after = picked.map(function (s) {
        var has = function (t) {
          return (s.tags || []).some(function (x) { return x.toLowerCase() === t.toLowerCase(); });
        };
        var keys = Object.keys(want);
        return { path: s.path,
                 add: keys.filter(function (t) { return want[t] === "on" && !has(t); }),
                 remove: keys.filter(function (t) { return want[t] === "off" && has(t); }) };
      }).filter(function (one) { return one.add.length || one.remove.length; });
      closeTagger();
      if (!after.length) { barSays("Nothing changed."); return; }
      apply.disabled = true;
      var was = after.map(function (one) { return { path: one.path, add: one.remove, remove: one.add }; });
      keepTags(after).then(function (results) {
        afterTags(results, "Tags changed on " + after.length + (after.length === 1 ? " study." : " studies."),
          function () {
            keepTags(was).then(function (back) { afterTags(back, "The tags are back as they were."); });
          });
      });
    });
    foot.append(cancel, apply);
    box.appendChild(foot);
    box.addEventListener("keydown", function (event) {
      if (event.key === "Escape") { event.preventDefault(); closeTagger(); el("studies-tag-chosen").focus(); }
    });
    el("studies-bar").appendChild(box);
    var at = el("studies-tag-chosen");
    var left = at.getBoundingClientRect().left - el("studies-bar").getBoundingClientRect().left - 40;
    box.style.left = Math.max(0, left) + "px";
    adding.focus();
  }

  /* Archived: left out of All studies, kept in the study's own
   * study_tags.json; nothing moved or deleted. Archived shows them. */
  function setArchived(on) {
    var list = chosen.map(studyAt).filter(Boolean);
    if (!list.length) return;
    Promise.all(list.map(function (one) {
      return fetch("/api/study-tags", {
        method: "POST",
        headers: { "content-type": "application/json" },
        // Archived alone: the tags as kept are not sent, nor touched.
        body: JSON.stringify({ path: one.path, archived: on })
      }).then(function (r) { return r.json(); }).then(function (d) {
        if (!d || !d.ok) return { path: one.path, reason: (d && (d.reason || d.error)) || "not kept" };
        // A server that does not pass "archived" on keeps the tags alone.
        if (!!d.archived !== on) return { path: one.path, reason: "this GUI's server does not keep Archive yet" };
        one.archived = on;
        edited += 1;
        return { path: one.path };
      }).catch(function () { return { path: one.path, reason: "The server did not answer." }; });
    })).then(function (results) {
      var done = results.filter(function (r) { return !r.reason; }).map(function (r) { return r.path; });
      var refused = results.filter(function (r) { return r.reason; });
      if (refused.length && !done.length) {
        // Nothing kept: still chosen, as they were.
        barSays((nameOf(refused[0].path) || "A study") + ": " + refused[0].reason + ".");
        return;
      }
      // Those kept leave the view; any refused stay chosen.
      chosen = refused.map(function (r) { return r.path; });
      lastPicked = null;
      if (activeView === "archived" && !studies.some(function (s) { return s.archived; })) activeView = "all";
      draw();
      if (refused.length) {
        barSays(done.length + " of " + results.length + " " + (on ? "archived" : "back") + "; " +
          (nameOf(refused[0].path) || "a study") + ": " + refused[0].reason + ".");
        return;
      }
      var said = done.length === 1 ? nameOf(done[0]) : done.length + " studies";
      barSays(on ? said + " archived: left out of All studies, nothing moved."
        : said + " back in All studies.", function () {
        chosen = done.slice();
        setArchived(!on);
      });
    });
  }

  function openStudy(path, button, page) {
    // Work a switch would lose holds it: said beside what was pressed.
    // One that goes on is in flight until answered (the strip is busy).
    var dashboard = window.FastMDXDashboard;
    var held = switchBegun();
    if (held) { notOpened(button, held); return; }
    button.disabled = true;
    fetch("/api/explore/switch", {
      method: "POST",
      headers: { "content-type": "application/json" },
      body: JSON.stringify({ folder: path })
    }).then(function (r) { return r.json(); }).then(function (d) {
      button.disabled = false;
      if (d && d.ok) {
        dashboard.studyOpened(d, page);
      } else {
        switchEnded();
        notOpened(button, (d && d.error) || "Could not open it.");
      }
    }).catch(function () {
      button.disabled = false;
      switchEnded();
      notOpened(button, "The server did not answer.");
    });
  }

  function switchBegun() {
    var dashboard = window.FastMDXDashboard;
    return dashboard && dashboard.beginStudySwitch ? dashboard.beginStudySwitch() : "";
  }

  function switchEnded() {
    var dashboard = window.FastMDXDashboard;
    if (dashboard && dashboard.endStudySwitch) dashboard.endStudySwitch();
  }

  /* Why a study did not open, said beside its button, which is an icon (the
   * page's notice reads it out, so this is not a second one). */
  function notOpened(button, why) {
    ICONS.flash(button, false, why);
    var host = button.parentNode;
    var said = host.querySelector(".study-open-said");
    if (!said) {
      said = make("span", "study-open-said muted small");
      host.appendChild(said);
    }
    said.textContent = why;
  }

  /* Tags or an archive kept from this page, counted: a list asked for
   * before one was kept, and answered after, is not drawn over it. */
  var edited = 0;

  function load(where) {
    if (loading) return loading;
    var query = where ? "?" + new URLSearchParams({ path: where }) : "";
    var began = edited;
    loading = fetch("/api/studies" + query)
      .then(function (r) { return r.json(); })
      .then(function (data) {
        loading = null;
        // Read before tags were kept here: read again rather than draw
        // the tags as they were over those just kept.
        if (began !== edited) return load(where);
        if (!data || !data.ok) {
          studies = [];
          el("studies-where").textContent = (data && (data.reason || data.error)) ||
            "The studies could not be listed.";
          draw();
          return;
        }
        lookedIn = data.root;
        studies = data.studies || [];
        // Listed again: a study's runs may have new numbers since.
        runsOf = {};
        if (!where) counted(data);
        tagsUsed = data.tags_used || [];
        offerTags();
        if (tagged && !studies.some(function (s) { return hasTag(s, tagged); })) tagged = null;
        chosen = chosen.filter(function (path) {
          return studies.some(function (s) { return s.path === path; });
        });
        el("studies-where").textContent = studies.length + " stud" +
          (studies.length === 1 ? "y" : "ies") + " in " + data.root +
          (data.more ? ", the newest first; there are more than this shows." : ".");
        draw();
      })
      .catch(function () {
        loading = null;
        el("studies-where").textContent = "The server did not answer.";
      });
    return loading;
  }

  function attach() {
    if (!el("studies-grid")) return;
    el("studies-search").addEventListener("input", draw);
    view = keptView();
    var kept = readKept(LAYOUT_KEPT, {});
    layout = { group: ["study", "system", "tag", "state", "folder", "none"].indexOf(kept.group) >= 0 ? kept.group : "study",
               density: kept.density === "compact" ? "compact" : "comfortable",
               hidden: Array.isArray(kept.hidden) ? kept.hidden.filter(function (k) { return typeof k === "string"; }) : [] };
    savedViews = (readKept(VIEWS_KEPT, []) || []).filter(function (v) { return v && typeof v.name === "string"; });
    Array.prototype.forEach.call(document.querySelectorAll("[data-studies-view]"), function (b) {
      var on = b.getAttribute("data-studies-view") === view;
      b.classList.toggle("active", on);
      b.setAttribute("aria-pressed", String(on));
      b.addEventListener("click", function () { showAs(b.getAttribute("data-studies-view"), true); });
    });
    el("studies-tag-filter-clear").addEventListener("click", function () { narrowTo(null); });
    el("studies-compare").addEventListener("click", compare);
    el("studies-clear").addEventListener("click", function (event) { clearChosen(event.detail === 0); });
    el("studies-tag-chosen").addEventListener("click", openTagger);
    el("studies-files-chosen").addEventListener("click", function () {
      if (chosen.length === 1) openStudy(chosen[0], el("studies-files-chosen"), "files");
    });
    el("studies-archive").addEventListener("click", function () { setArchived(activeView !== "archived"); });
    // A click outside the tagger closes it, as a menu does.
    document.addEventListener("click", function (event) {
      var open = document.querySelector(".studies-tagger");
      if (open && !open.contains(event.target) && !event.target.closest("#studies-tag-chosen")) closeTagger();
    });
    var path = el("studies-path");
    el("studies-look-in").addEventListener("click", function () {
      if (window.FastMDXPicker) window.FastMDXPicker.open({ into: "studies-path", mode: "folder" });
    });
    path.addEventListener("change", function () {
      if (path.value.trim()) load(path.value.trim());
    });
    wireOpenShared();
    keepThePlace();
    var opened = function (detail) {
      var page = detail && detail.page ? detail.page : detail;
      if (page === "studies") {
        var back = load(lookedIn);
        if (back && back.then) back.then(returnToThePlace, returnToThePlace);
        else returnToThePlace();
      }
    };
    if (window.FastMDXDashboard && window.FastMDXDashboard.on) {
      window.FastMDXDashboard.on("navigate", opened);
    }
    if (location.hash === "#studies") load();
  }

  /* Where All studies was scrolled to, and its search, kept when the page
   * is left (for Compare, or a study opened, which loads the page again)
   * and given back when it is shown again: every page opens at its top
   * (dashboard.js), and the two chosen were below the fold on Back. */
  var PLACE_KEPT = "fmx.studiesPlace";
  var keptScroll = 0;
  var returning = null;

  function column() { return document.querySelector(".main"); }

  function keepThePlace() {
    var page = document.querySelector('.page[data-page="studies"]');
    var search = el("studies-search");
    var main = column();
    if (!page || !search || !main) return;
    try {
      var kept = JSON.parse(sessionStorage.getItem(PLACE_KEPT) || "null");
      if (kept && typeof kept.search === "string" && !search.value) search.value = kept.search;
      if (kept && kept.scroll > 0) { keptScroll = kept.scroll; returning = kept.scroll; }
    } catch (error) { /* kept for this load only */ }
    main.addEventListener("scroll", function () {
      if (page.hidden || returning !== null) return;
      keptScroll = main.scrollTop;
    }, { passive: true });
    // Shown again: the page's own reset to the top (a scroll that comes
    // after) is not where it was left.
    new MutationObserver(function () {
      if (page.hidden) {
        try {
          sessionStorage.setItem(PLACE_KEPT, JSON.stringify(
            { search: search.value || "", scroll: keptScroll }));
        } catch (error) { /* kept in this page only */ }
      } else {
        returning = keptScroll;
      }
    }).observe(page, { attributes: true, attributeFilter: ["hidden"] });
    window.addEventListener("pagehide", function () {
      if (page.hidden) return;
      try {
        sessionStorage.setItem(PLACE_KEPT, JSON.stringify(
          { search: search.value || "", scroll: main.scrollTop }));
      } catch (error) { /* nothing to keep it in */ }
    });
  }

  function returnToThePlace() {
    var main = column();
    if (returning === null || !main) return;
    var to = returning;
    requestAnimationFrame(function () {
      main.scrollTop = to;
      keptScroll = main.scrollTop;
      returning = null;
    });
  }

  /* A shared study, from its DOI, its Zenodo address or a zip on this
   * computer: downloaded, checked file by file and opened in a folder of
   * its own in the workspace (POST /api/open-shared). */
  function wireOpenShared() {
    var toggle = el("studies-open-shared");
    var form = el("studies-open-form");
    if (!toggle || !form) return;
    var input = el("studies-open-source");
    var go = el("studies-open-go");
    var said = el("studies-open-said");
    toggle.addEventListener("click", function () {
      form.hidden = !form.hidden;
      toggle.setAttribute("aria-expanded", String(!form.hidden));
      if (!form.hidden) input.focus();
    });
    form.addEventListener("submit", async function (event) {
      event.preventDefault();
      var source = input.value.trim();
      if (!source) { input.focus(); return; }
      var held = switchBegun();
      if (held) { said.textContent = held; return; }
      go.disabled = true;
      input.disabled = true;
      said.textContent = "Fetching it and checking every file\u2026";
      try {
        var response = await fetch("/api/open-shared", {
          method: "POST", headers: { "content-type": "application/json" },
          body: JSON.stringify({ source: source }),
        });
        var answer = await response.json();
        if (answer && answer.ok) {
          said.textContent = "Opened. Every file matched its SHA-256.";
          window.FastMDXDashboard.studyOpened(answer);
          return;
        }
        said.textContent = (answer && answer.error) || "It could not be opened.";
      } catch (error) {
        said.textContent = "The server did not answer.";
      }
      switchEnded();
      go.disabled = false;
      input.disabled = false;
    });
  }

  /* The sidebar's part: how many studies the workspace holds, beside All
   * studies, and the newest of them under Recent, to switch to. */
  var workspace = null;
  var asked = null;

  function counted(data) {
    workspace = data && data.ok ? data : null;
    var count = el("nav-studies-count");
    if (count) {
      count.textContent = workspace
        ? String(workspace.studies.length) + (workspace.more ? "+" : "") : "";
    }
  }

  /* Asked once and kept, unless `fresh`: the sidebar's Recent asks again
   * when another study is opened, and a while after it last asked, so a
   * study's state and a study begun since are said. */
  var askedAt = 0;
  function askTheWorkspace(fresh) {
    if (fresh) asked = null;
    if (asked) return asked;
    askedAt = Date.now();
    asked = fetch("/api/studies")
      .then(function (r) { return r.json(); })
      .then(function (data) { counted(data); return workspace; })
      .catch(function () { asked = null; return null; });
    return asked;
  }

  /* How many of the workspace's studies share each shown ID: one shared
   * is named by its folder too, in Recent and on the strip of open studies
   * alike. */
  function idsShared(data) {
    var times = {};
    (data && data.studies || []).forEach(function (study) {
      var id = shownId(study) || study.name;
      times[id] = (times[id] || 0) + 1;
    });
    return times;
  }

  function sameFolder(a, b) {
    return String(a || "").replace(/\/+$/, "") === String(b || "").replace(/\/+$/, "");
  }

  /* The sidebar's Recent: the workspace's newest studies, the open one
   * marked, each opened with a click. Folded or not as last left. Each
   * name is shortened in its middle to the width there is, as the active
   * study's folder is, and its state is a light: the word is on hover and
   * read out with the name. */
  var SHOWN_RECENT = 6;
  var STATE_WORDS = { running: "Running", completed: "Completed", stopped: "Stopped",
                      failed: "Failed", interrupted: "Interrupted",
                      incomplete: "Incomplete", "not started": "Not started" };

  function stateWord(state) {
    var key = String(state || "").toLowerCase();
    return STATE_WORDS[key] || (key ? key.charAt(0).toUpperCase() + key.slice(1) : "");
  }

  function fitRecentNames() {
    var dashboard = window.FastMDXDashboard;
    if (!dashboard || !dashboard.fitInTheMiddle) return;
    Array.prototype.forEach.call(document.querySelectorAll(".sidebar-recent-name"),
      function (shown) {
        var room = shown.getBoundingClientRect().width;
        if (room > 0) dashboard.fitInTheMiddle(shown, shown.dataset.name || "", room);
      });
  }

  function showRecent(fresh, evenUnderThePointer) {
    var list = el("sidebar-recent-list");
    if (!list) return Promise.resolve();
    return askTheWorkspace(fresh).then(function (data) {
      // Not built again under somebody's keyboard, nor, unless another
      // study was opened, under their pointer.
      if (list.contains(document.activeElement)) return;
      if (!evenUnderThePointer && list.matches(":hover")) return;
      var open = ((el("sidebar-output-folder") || {}).textContent || "").trim();
      var recent = (data && data.studies || []).slice(0, SHOWN_RECENT);
      if (!recent.length) {
        list.replaceChildren(make("div", "sidebar-recent-empty",
          data ? "No studies in this workspace yet." : "The studies could not be listed."));
        return;
      }
      // Several studies of one structure share its ID: each is named by
      // its folder too, or Recent, and ⌘K from it, list "1L2Y" twice.
      var times = idsShared(data);
      list.replaceChildren.apply(list, recent.map(function (study) {
        var here = sameFolder(study.path, open);
        var word = stateWord(study.state);
        var name = shownId(study) || study.name;
        if (times[name] > 1 && study.name && study.name !== name) name += " \u00b7 " + study.name;
        var item = make("button", "sidebar-recent-item");
        item.type = "button";
        item.title = (word ? word + "\n" : "") + study.path;
        item.dataset.path = study.path;
        if (here) item.setAttribute("aria-current", "true");
        var light = make("span", "recent-light");
        light.dataset.state = String(study.state || "").toLowerCase();
        light.setAttribute("aria-hidden", "true");
        var shown = make("span", "sidebar-recent-name", name);
        shown.dataset.name = name;
        item.append(shown, light);
        // Named in full, as the shortened name on screen is not.
        item.setAttribute("aria-label", name + (word ? ", " + word : ""));
        item.addEventListener("click", function () {
          if (!here) openStudy(study.path, item);
          else if (window.FastMDXDashboard) window.FastMDXDashboard.navigate("overview");
        });
        return item;
      }));
      fitRecentNames();
    });
  }

  /* Fitted again as the sidebar's width changes, and once its typeface
   * has come. */
  function fitRecentAsTheSidebarChanges() {
    var list = el("sidebar-recent-list");
    if (!list) return;
    var width = 0;
    if (window.ResizeObserver) {
      new ResizeObserver(function (seen) {
        var now = seen[0] ? seen[0].contentRect.width : 0;
        // Hidden, nothing is fitted; shown again, it is, at whatever width.
        if (!now) { width = 0; return; }
        if (Math.abs(now - width) > 0.5) { width = now; fitRecentNames(); }
      }).observe(list);
    }
    if (document.fonts && document.fonts.ready) document.fonts.ready.then(fitRecentNames);
  }

  function keepRecentFolded() {
    var fold = el("sidebar-recent");
    if (!fold) return;
    try { if (localStorage.getItem("fmx.recentOpen") === "0") fold.open = false; } catch (e) {}
    fold.addEventListener("toggle", function () {
      try { localStorage.setItem("fmx.recentOpen", fold.open ? "1" : "0"); } catch (e) {}
    });
  }

  /* With the sidebar folded to its strip, Recent and Chats are their icons
   * there, and a click shows the list beside the strip, to open one; a
   * click elsewhere, Escape, or one opened puts it away. The fold kept for
   * the sidebar shown is left as it was. */
  function folded() {
    return document.body.classList.contains("sidebar-collapsed")
      && !document.body.classList.contains("sidebar-peek")
      && window.matchMedia("(min-width: 801px)").matches;
  }

  function besideTheStrip(fold, list, item, onShow) {
    if (!fold || !list) return;
    var head = fold.querySelector("summary");
    var wasOpen = fold.open;

    function shut(back) {
      if (!fold.classList.contains("flyout")) return;
      fold.classList.remove("flyout");
      head.setAttribute("aria-expanded", "false");
      fold.open = wasOpen;
      list.style.top = "";
      if (back) head.focus();
    }

    /* Beside its icon, and kept inside the window: Chats sits low in the
     * strip, and its list would run off the foot. */
    function place() {
      var at = head.getBoundingClientRect();
      var room = window.innerHeight - list.offsetHeight - 8;
      list.style.top = Math.round(Math.max(8, Math.min(at.top, room))) + "px";
    }

    head.addEventListener("click", function (event) {
      if (!folded()) return;
      // The head's own buttons (a new chat) are not the fold.
      if (event.target.closest("button")) return;
      event.preventDefault();
      if (fold.classList.contains("flyout")) { shut(false); return; }
      wasOpen = fold.open;
      fold.classList.add("flyout");
      fold.open = true;
      head.setAttribute("aria-expanded", "true");
      place();
      Promise.resolve(onShow()).then(function () {
        if (!fold.classList.contains("flyout")) return;
        place();
        var first = list.querySelector(item);
        if (first && event.detail === 0) first.focus();
      });
    });
    // Its toggle is not the fold kept for the sidebar shown.
    fold.addEventListener("toggle", function (event) {
      if (fold.classList.contains("flyout")) event.stopImmediatePropagation();
    }, true);
    list.addEventListener("click", function (event) {
      if (event.target.closest(item)) shut(false);
    });
    document.addEventListener("click", function (event) {
      if (!fold.contains(event.target)) shut(false);
    });
    document.addEventListener("keydown", function (event) {
      if (event.key === "Escape" && fold.classList.contains("flyout")) shut(true);
    });
    // The sidebar shown again, or folded: the list goes back in place.
    new MutationObserver(function () { if (!folded()) shut(false); })
      .observe(document.body, { attributes: true, attributeFilter: ["class"] });
  }

  function recentBesideTheStrip() {
    besideTheStrip(el("sidebar-recent"), el("sidebar-recent-list"), ".sidebar-recent-item",
      function () { return showRecent(false, true); });
    besideTheStrip(el("sidebar-chats"), el("sidebar-chats-list"), ".sidebar-conv, .sidebar-conv-new",
      function () { return window.FastMDXChats ? window.FastMDXChats.show() : null; });
  }

  /* Beyond loopback /api/compared is refused (the page's bind flag, as the
   * palette leaves out Compare with...): the way to the Compare page is
   * not offered; nor are the bar's Tag, Files and Archive, which POST. */
  function compareOnLoopbackOnly() {
    if (!beyondLoopback()) return;
    var link = el("studies-compare-page");
    if (link) link.remove();
    // Hidden, not taken out: the page's handlers are given them as ever.
    Array.prototype.forEach.call(document.querySelectorAll(
      "#studies-compare, #studies-tag-chosen, #studies-files-chosen, #studies-archive, " +
      "#studies-bar .studies-bar-sep"), function (node) { node.hidden = true; });
  }

  if (document.readyState === "loading") {
    document.addEventListener("DOMContentLoaded", compareOnLoopbackOnly);
    document.addEventListener("DOMContentLoaded", attach);
    document.addEventListener("DOMContentLoaded", recentBesideTheStrip);
    document.addEventListener("DOMContentLoaded", keepRecentFolded);
    document.addEventListener("DOMContentLoaded", fitRecentAsTheSidebarChanges);
  } else {
    compareOnLoopbackOnly();
    attach();
    recentBesideTheStrip();
    keepRecentFolded();
    fitRecentAsTheSidebarChanges();
  }
  setTimeout(showRecent, 1500);
  // Again once the open study is known, or another is opened.
  var recentFor = null;
  var RECENT_AGAIN_MS = 30000;
  window.addEventListener("dashboard:app-state", function () {
    var open = ((el("sidebar-output-folder") || {}).textContent || "").trim();
    var stale = Date.now() - askedAt > RECENT_AGAIN_MS;
    if (open === recentFor && !stale) return;
    var moved = open !== recentFor;
    recentFor = open;
    showRecent(moved || stale, moved);
  });

  // The open study's ID given or cleared: shown so at once.
  document.addEventListener("change", function (event) {
    if (!event.target || event.target.id !== "setting-run-name") return;
    showRecent(false, true);
    if (studies.length) draw();
  });

  window.FastMDXStudies = { load: load, narrowTo: narrowTo,
                            showRecent: showRecent,
                            idsShared: function () { return askTheWorkspace(false).then(idsShared); } };
}());
