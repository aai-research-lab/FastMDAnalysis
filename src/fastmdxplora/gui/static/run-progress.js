/* How far a run has got, and the time it has left, as the Overview and the
 * status bar say them: the share of its steps taken (progress_percent where
 * its record gives one), and the steps left over the steps a second of its
 * newest samples. The run's whole time so far includes building the system
 * and the first, slow samples, so early in a run it said 33 minutes with
 * about six left. These are dashboard.js's own `progressPercent`,
 * `stepsPerSecond` and `computeETA`, kept here for every page that says a
 * run's progress (All studies' rows), so they are worked out one way.
 *
 * For a study of several runs (studies.js): its steps taken over its steps
 * planned, every run's, begun or not; its time left as it runs them, one
 * after another (the sum) or side by side (each worker given the next run
 * as it frees, the last to end), and not known where how is not known. */
(function () {
  "use strict";

  function clampFloat(value, low, high, fallback) {
    var number = Number(value);
    if (!Number.isFinite(number)) return fallback;
    return Math.max(low, Math.min(high, number));
  }

  function progressPercent(status) {
    if (Number.isFinite(Number(status.progress_percent))) {
      return clampFloat(Number(status.progress_percent), 0, 100, null);
    }
    if (status.current_step != null && status.total_planned_steps != null &&
        Number(status.total_planned_steps) > 0) {
      return clampFloat(Number(status.current_step) / Number(status.total_planned_steps) * 100,
                        0, 100, null);
    }
    return null;
  }

  /* Steps a second from the samples the record saw taken: from where the
   * steps began (the last sample at step 0, minimisation's end), and since
   * the last long pause (a study carried on in pieces has the time between
   * them in its record, and its rate jumped). */
  function stepsPerSecond(rows) {
    var samples = (rows || []).map(function (row) {
      return { step: Number(row.step), at: Date.parse(row.timestamp || "") };
    }).filter(function (s) { return s.step >= 0 && Number.isFinite(s.at); });
    var began = samples.map(function (s) { return s.step; }).lastIndexOf(0);
    if (began > 0) samples = samples.slice(began);
    var recent = samples.slice(-12);
    var gaps = recent.slice(1).map(function (s, i) { return s.at - recent[i].at; });
    var usual = gaps.slice().sort(function (a, b) { return a - b; })[Math.floor(gaps.length / 2)] || 0;
    var pause = gaps.reduce(function (at, gap, i) {
      return gap > Math.max(60000, 5 * usual) ? i + 1 : at;
    }, 0);
    recent = recent.slice(pause);
    if (recent.length < 2) return 0;
    var first = recent[0];
    var last = recent[recent.length - 1];
    var seconds = (last.at - first.at) / 1000;
    return seconds > 0 && last.step > first.step ? (last.step - first.step) / seconds : 0;
  }

  /* Seconds left for one run going; null where nothing says. */
  function secondsLeft(status, rows) {
    var step = Number(status.current_step);
    var total = Number(status.total_planned_steps);
    if (!(step > 0) || !(total > 0) || step >= total) return null;
    var rate = stepsPerSecond(rows);
    return rate > 0 ? (total - step) / rate : null;
  }

  /* A study's progress and time left from its runs (workspace.py
   * `_how_far`'s `pace`): {percent, seconds}, each null where not known. */
  function ofPace(pace) {
    var runs = (pace && pace.runs) || [];
    if (pace && pace.mode === "one") {
      var only = runs[0] || {};
      if (only.state !== "going" || !only.status) return { percent: null, seconds: null };
      return { percent: progressPercent(only.status), seconds: secondsLeft(only.status, only.samples) };
    }
    var taken = 0, planned = 0, steps = true;
    runs.forEach(function (run) {
      if (run.steps == null) { steps = false; return; }
      planned += Number(run.steps);
      taken += run.state === "waiting" ? 0 : Math.min(Number(run.taken) || 0, Number(run.steps));
    });
    var percent = steps && planned > 0 ? clampFloat(taken / planned * 100, 0, 100, null) : null;
    // The rate of the runs going, for the runs not begun.
    var rate = 0;
    runs.forEach(function (run) {
      if (!rate && run.state === "going") rate = stepsPerSecond(run.samples);
    });
    var going = [], waiting = [], known = true;
    runs.forEach(function (run) {
      if (run.state === "going") {
        var left = run.status ? secondsLeft(run.status, run.samples) : null;
        if (left == null) known = false; else going.push(left);
      } else if (run.state === "waiting") {
        if (run.steps == null || !(rate > 0)) known = false; else waiting.push(run.steps / rate);
      } else if (run.state === "not known") {
        known = false;
      }
    });
    var seconds = null;
    if (known && pace.mode === "sequential") {
      seconds = going.concat(waiting).reduce(function (a, b) { return a + b; }, 0);
    } else if (known && pace.mode === "parallel" && pace.workers > 0 && going.length <= pace.workers) {
      var free = going.slice();
      while (free.length < pace.workers) free.push(0);
      waiting.forEach(function (one) {
        var soonest = free.indexOf(Math.min.apply(null, free));
        free[soonest] += one;
      });
      seconds = Math.max.apply(null, free);
    }
    return { percent: percent, seconds: seconds };
  }

  /* A time left as the status bar says it. */
  function durationSaid(seconds) {
    var total = Math.max(0, Math.floor(Number(seconds) || 0));
    var hours = Math.floor(total / 3600);
    var minutes = Math.floor((total % 3600) / 60);
    var secs = total % 60;
    return hours ? hours + "h " + minutes + "m" : minutes + "m " + secs + "s";
  }

  window.FastMDXRunProgress = { progressPercent: progressPercent, stepsPerSecond: stepsPerSecond,
                                secondsLeft: secondsLeft, ofPace: ofPace,
                                durationSaid: durationSaid };
}());
