/* The scheme before the first paint, as frame.js chooses it: Light or
 * Dark where one was chosen (Graphite and Ink, chosen before, are Dark;
 * Paper is Light), otherwise the computer's. Loaded in the page's head, so
 * a light computer's page does not open dark and turn light.
 *
 * The accent too: the person's own (`fmx.accent`), else the one the
 * service hosting the GUI chose (`data-accent-default` on the root
 * element, `fastmdx gui --accent`), else Cyan, so the page never shows
 * Cyan first and turns. */
(function () {
  var kept = null, former = { graphite: "dark", ink: "dark", paper: "light" };
  try { kept = localStorage.getItem("fmx.theme"); } catch (e) { /* private mode */ }
  kept = former[kept] || kept;
  if (kept !== "light" && kept !== "dark") {
    try { kept = window.matchMedia("(prefers-color-scheme: light)").matches ? "light" : "dark"; }
    catch (e) { kept = "dark"; }
  }
  var root = document.documentElement;
  root.dataset.theme = kept;
  var accents = ["cyan", "blue", "violet", "magenta", "graphite"];
  var accent = null;
  try { accent = localStorage.getItem("fmx.accent"); } catch (e) { /* private mode */ }
  if (accents.indexOf(accent) < 0) accent = root.getAttribute("data-accent-default");
  if (accents.indexOf(accent) >= 0 && accent !== "cyan") root.dataset.accent = accent;
  else root.removeAttribute("data-accent");
})();
