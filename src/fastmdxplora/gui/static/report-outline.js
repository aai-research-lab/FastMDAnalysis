/* The report's outline beside it, as an editor's outline beside a file:
 * its sections numbered as the PDF numbers them (2, 2.1), the one being
 * read marked, each a way to it. And each number the report states, on a
 * click, the record it was traced to (report_numbers.py): the file, the
 * place in it and the value as recorded; or that no record holds it, with
 * the nearest of its kind, and the records changed after the report was
 * written.
 */
(function () {
  "use strict";

  function byId(id) { return document.getElementById(id); }

  function make(tag, cls, text) {
    var node = document.createElement(tag);
    if (cls) node.className = cls;
    if (text !== undefined && text !== null) node.textContent = text;
    return node;
  }

  var numbers = [];
  var changed = [];
  var headings = [];
  var card = null;

  /* Each h2 and h3 numbered: 1, 2, 2.1, as the PDF's stylesheet counts
   * them (report/pdf.py). The title is not a section. */
  function number(doc) {
    var major = 0, minor = 0;
    headings = [];
    Array.prototype.forEach.call(doc.querySelectorAll("h2, h3"), function (heading, i) {
      if (heading.tagName === "H2") { major += 1; minor = 0; } else { minor += 1; }
      var said = heading.tagName === "H2" ? String(major) : major + "." + minor;
      heading.dataset.section = said;
      if (!heading.id) heading.id = "report-section-" + said.replace(".", "-");
      headings.push({ node: heading, said: said, level: heading.tagName === "H2" ? 2 : 3,
                      words: heading.textContent.replace(/\s+/g, " ").trim(), at: i });
    });
  }

  /* The section gone to from the outline: kept marked while its heading
   * is in sight, though the report's end cannot scroll it to the top. */
  var chosen = null;

  function go(heading) {
    chosen = heading;
    heading.scrollIntoView({ block: "start" });
    mark(heading);
  }

  function mark(current) {
    var outline = byId("report-outline");
    if (!outline) return;
    Array.prototype.forEach.call(outline.querySelectorAll(".outline-item"), function (item) {
      var here = current && item.dataset.target === current.id;
      if (here) item.setAttribute("aria-current", "location");
      else item.removeAttribute("aria-current");
    });
  }

  /* The section being read: the last heading in the top third of the
   * window or above it. */
  var asked = false;
  function follow() {
    if (asked) return;
    asked = true;
    requestAnimationFrame(function () {
      asked = false;
      var doc = byId("report-document");
      var page = document.querySelector('.page[data-page="report"]');
      if (!doc || doc.hidden || !page || page.hidden || !headings.length) return;
      if (chosen) {
        var top = chosen.getBoundingClientRect().top;
        if (top >= 0 && top < window.innerHeight - 40) { mark(chosen); return; }
        chosen = null;
      }
      var reading = null;
      var line = Math.max(140, window.innerHeight / 3);
      headings.forEach(function (h) {
        if (h.node.getBoundingClientRect().top <= line) reading = h.node;
      });
      mark(reading || headings[0].node);
    });
  }
  document.addEventListener("scroll", follow, true);

  function counts() {
    var traced = numbers.filter(function (n) { return n.traced; }).length;
    return { traced: traced, missing: numbers.length - traced };
  }

  function render(doc, data) {
    var outline = byId("report-outline");
    numbers = (data && data.numbers) || [];
    changed = (data && data.records_changed) || [];
    close();
    if (!outline) return;
    outline.replaceChildren();
    if (!doc || doc.hidden) { outline.hidden = true; return; }
    number(doc);
    outline.hidden = !headings.length;
    if (!headings.length) return;

    var toggle = make("button", "outline-toggle", "Outline");
    toggle.type = "button";
    toggle.setAttribute("aria-expanded", "false");
    toggle.appendChild(make("span", "muted", " \u00b7 " + headings.filter(function (h) {
      return h.level === 2;
    }).length + " sections"));
    toggle.addEventListener("click", function () {
      var open = outline.dataset.open !== "true";
      outline.dataset.open = String(open);
      toggle.setAttribute("aria-expanded", String(open));
    });
    outline.appendChild(toggle);

    var body = make("div", "outline-body");
    body.appendChild(make("div", "outline-title", "Outline"));
    var list = make("ol", "outline-list");
    var under = list;
    headings.forEach(function (h) {
      var li = make("li", "outline-entry outline-level-" + h.level);
      var item = make("button", "outline-item");
      item.type = "button";
      item.dataset.target = h.node.id;
      item.appendChild(make("span", "outline-number", h.said));
      item.appendChild(make("span", "outline-words", h.words));
      item.title = h.said + " " + h.words;
      item.addEventListener("click", function () {
        go(h.node);
        if (outline.dataset.open === "true") toggle.click();
      });
      li.appendChild(item);
      if (h.level === 2) {
        list.appendChild(li);
        under = make("ol", "outline-list outline-sub");
        li.appendChild(under);
      } else {
        under.appendChild(li);
      }
    });
    Array.prototype.forEach.call(list.querySelectorAll(".outline-sub"), function (sub) {
      if (!sub.childElementCount) sub.remove();
    });
    body.appendChild(list);
    body.appendChild(said(doc));
    outline.appendChild(body);
    follow();
  }

  /* What the report's numbers trace to, in a line under the outline. */
  function said(doc) {
    var box = make("div", "outline-numbers");
    box.id = "report-numbers-said";
    if (!numbers.length) {
      box.appendChild(make("p", "muted", "No number in it to trace to a record."));
      return box;
    }
    var n = counts();
    box.appendChild(make("p", "", n.traced + " of " + numbers.length + " numbers traced to a record; " +
      "click one for it."));
    if (n.missing) {
      var show = make("button", "plot-tool outline-mark",
        "Mark the " + n.missing + " not found");
      show.type = "button";
      show.setAttribute("aria-pressed", "false");
      show.title = "Mark each number no record of this study holds";
      show.addEventListener("click", function () {
        var on = doc.dataset.marked !== "true";
        doc.dataset.marked = String(on);
        show.setAttribute("aria-pressed", String(on));
        show.textContent = on ? "Unmark them" : "Mark the " + n.missing + " not found";
        if (on) {
          var first = doc.querySelector('.report-number[data-traced="false"]');
          if (first) first.scrollIntoView({ block: "center" });
        }
      });
      box.appendChild(show);
    }
    if (changed.length) {
      box.appendChild(make("p", "outline-changed",
        changed.length + (changed.length === 1 ? " record" : " records") +
        " changed after the report was written: " + changed.join(", ") +
        ". Write it again to have it say them."));
    }
    return box;
  }

  function close() {
    if (card) { card.remove(); card = null; }
    Array.prototype.forEach.call(document.querySelectorAll('.report-number[aria-expanded="true"]'),
      function (n) { n.setAttribute("aria-expanded", "false"); });
  }

  function plain(value) {
    var x = Number(value);
    if (!isFinite(x)) return String(value);
    return String(+x.toPrecision(10));
  }

  /* The record's place, as a path into the file: findings.mean.mean. */
  function where(record) {
    var line = make("div", "number-where");
    line.appendChild(make("code", "", record.file));
    line.appendChild(document.createTextNode(" "));
    line.appendChild(make("code", "", record.key + " = " + plain(record.value) +
      (record.unit ? " " + record.unit : "")));
    return line;
  }

  function open(span) {
    var found = numbers[Number(span.dataset.number)];
    if (!found) return;
    var wasOpen = span.getAttribute("aria-expanded") === "true";
    close();
    if (wasOpen) return;
    var doc = byId("report-document");
    card = make("div", "report-number-card");
    card.setAttribute("role", "dialog");
    card.setAttribute("aria-label", "Where " + found.said + " comes from");
    var head = make("div", "number-head");
    head.appendChild(make("b", "", found.said));
    var shut = make("button", "line-btn number-close", "\u00d7");
    shut.type = "button";
    shut.setAttribute("aria-label", "Close");
    shut.addEventListener("click", function () { close(); span.focus(); });
    head.appendChild(shut);
    card.appendChild(head);
    var record = found.record || found.nearest;
    if (found.traced) {
      card.dataset.traced = "true";
      card.appendChild(make("p", "", "From the study's record:"));
      card.appendChild(where(record));
      if (record.unit && record.in && record.unit !== record.in) {
        card.appendChild(make("p", "muted", "That is " + plain(record.as_said) + " " + record.in +
          ", as the report gives it."));
      }
      if (found.others) {
        card.appendChild(make("p", "muted", found.others === 1
          ? "One other record holds a value that rounds to it too."
          : found.others + " other records hold a value that rounds to it too."));
      }
    } else {
      card.dataset.traced = "false";
      card.appendChild(make("p", "", "No record of this study holds it."));
      if (record) {
        card.appendChild(make("p", "", "The nearest of its kind:"));
        card.appendChild(where(record));
      }
      // The records of what it is about, changed since the report said it.
      var stale = changed.filter(function (file) {
        return (record && file === record.file) || (found.subjects || []).some(function (subject) {
          return ("/" + file).indexOf("/" + subject + "/") >= 0;
        });
      });
      if (stale.length) {
        card.appendChild(make("p", "number-stale", (stale.length === 1 ? stale[0] + " was" :
          stale.join(", ") + " were") + " written again after the report: Write it again to have " +
          "the report say what the records hold now."));
      }
      card.appendChild(make("p", "muted",
        "A threshold of the method, a citation, or a number the report works out from more than " +
        "one record has none."));
    }
    if (record && window.FastMDXFrame && window.FastMDXFrame.previewPath) {
      var show = make("button", "plot-tool", "Open the record");
      show.type = "button";
      show.addEventListener("click", function () { window.FastMDXFrame.previewPath(record.file); });
      card.appendChild(show);
    }
    doc.appendChild(card);
    var box = doc.getBoundingClientRect(), at = span.getBoundingClientRect();
    var left = Math.max(8, Math.min(at.left - box.left, box.width - card.offsetWidth - 8));
    card.style.left = left + "px";
    card.style.top = (at.bottom - box.top + 6) + "px";
    span.setAttribute("aria-expanded", "true");
  }

  document.addEventListener("click", function (event) {
    var span = event.target.closest && event.target.closest(".report-number");
    if (span && span.closest("#report-document")) { open(span); return; }
    if (card && !card.contains(event.target)) close();
  });
  document.addEventListener("keydown", function (event) {
    var span = event.target.closest && event.target.closest(".report-number");
    if (span && (event.key === "Enter" || event.key === " ")) {
      event.preventDefault();
      open(span);
    } else if (event.key === "Escape" && card) {
      var from = document.querySelector('.report-number[aria-expanded="true"]');
      close();
      if (from) from.focus();
    }
  });

  window.FastMDXReportOutline = { render: render, close: close };
}());
