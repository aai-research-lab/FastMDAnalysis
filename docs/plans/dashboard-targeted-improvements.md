# FastMDXplora dashboard workplan

Proposed work • 15 October 2026 • Prince Otegbulu

Branch: princeote/dashboard-targeted-improvements. Baseline: a3b54aac99dc7cd13170360aa982b136c9f50570. Working checkout: C:/Users/User/OneDrive/Documents/GitHub/FastMDXplora.

Goal: dashboard figures built from their saved data should open the exact corresponding molecular-viewer frame when a point represents a frame. Reports and other exports remain PNG/SVG. Execution status and evidence are maintained in the coverage table below.

## 1. Scope and design rules

Use the current Mol* viewer, chart renderers, recorded analysis data, selection state and Agent current_view. Keep research bookmarks out. No editable PowerPoint chart project, new graph library, replacement viewer, new persistence system or unrelated visual redesign.

Preserve scientific calculations, original result files and full-analysis statistics. Display filtering must not silently recompute an analysis. Reuse existing modules; extract a shared helper only when two real callers need it. Keep cache sizes bounded and heavy optional imports lazy.

## 2. What the current code already does

gui/series.py and gui/static/series-chart.js plot time series and RMSF. gui/figure_data.py and gui/static/figure-chart.js cover secondary structure, clustering, projections, landscapes, density, thermodynamics, inertia and ligand contacts/interactions. Several already open frames. Extend missing behavior rather than rebuilding those charts.

Baseline source issue (addressed in Stage A): series-chart.js playedIndex chose the nearest sampled playback frame. molecule-viewer.js handled trajectory-seek using the playback index, ignoring trajectoryFrame. Stage A now carries the source frame into one shared exact-navigation path.

The prior audit reproduced per-atom RMSF being presented as residues. Graph residue focus also uses focusResidue, while Agent view hints read STATE.selection; reproduce the resulting interaction in a browser before changing it.

analysis/plotting.py already saves SVG companions. Reports/slides use static figures. Preserve that contract. report/summary_figure.py composes PNG panels inside its SVG; do not describe that composite as fully vector-editable or broaden this task into an export rewrite.

## 3. Establish the baseline once

Confirm the branch and imported package path, then run the existing relevant chart, frame, selection and Agent tests. Use the repository virtual environment; browser skips do not count as successful interaction checks.

Make one coverage table in this plan for every current chart/figure type: data source, frame/residue identity, existing click behavior, missing behavior, static export and acceptance test. Include the Analysis page and any chart instances on other dashboard pages. Start from SERIES and FIGURES/PLOTS and reconcile them with figures actually produced; avoid introducing a parallel registry just for planning.

Use small deterministic trajectory fixtures with distinguishable coordinates. Cover nonzero analysis start, stride, sampled playback, duplicate residue numbers across chains, atom selections and a second study. Record existing failures separately.

Follow CONTRIBUTING.md's issue-first workflow for nontrivial changes: prepare concrete issue text from each deliverable below. Submit externally only within the user's authorization. Do not reopen the old PR wholesale.

### Figure and interaction coverage

Baseline was inspected on 15 October 2026. The exact-frame regression was reproduced in Chromium/WebGL with a deterministic 11-frame trajectory and a three-frame browser sample: clicking source frame 7 showed source frame 5 (playback index 1), with the selected CA at x=5 Å instead of x=7 Å. Stage A now passes frame 7 → x=7.00 Å, and its focused implementation/verification evidence is recorded in the coverage table.

| Stage and chart instance | Saved identity and current behavior | Static export | Baseline evidence and next action |
|---|---|---|---|
| A/D — Analysis time series (RMSD, radius of gyration, hydrogen bonds, total SASA, Q-value, ligand RMSD, and other recorded per-frame means) | `analysis/<name>/<name>.dat` plus `analysis_manifest.json` records analysis frame, start, stride, and clock. The chart carries source frame and trajectory signature; the viewer reuses a sampled frame or requests one exact frame and maintains a source-ordered logical playback map. | Keep original PNG/SVG and report figures unchanged. | Baseline Chromium/WebGL repro was frame 7 → 5 (CA x=5 Å). Stage A now opens frame 7 at x=7.00 Å with sample `[0, 5, 10]`; frame 0/final, start/stride, stale/missing sources, interrupted loads, rapid selection and study switch are covered. The 10,000-atom measurement remains bounded to a 27,316-byte one-frame XTC for frame 73 (source DCD 12,122,700 bytes), 0.457 s warm; coordinates matched x=73 Å. Earlier A–D browser/playback group: 64 passed; the Windows default-encoding receipt case passed its `-X utf8` rerun. A transient live-overlay assertion was made atomic with the frame state after a browser polling race; the latest exact-frame/navigation group passed 46/46. `ruff --select F,B`, GUI parsing, the concurrent frame-writer race test and built-wheel asset smoke pass. Stages A–D verified. |
| A/C — Viewer transport frame series | `/api/frames-info` supplies sampled playback frames and their recorded times; clicking the transport series follows the selected loaded source frame. | Viewer-only control; not an analysis export. | Exact-source labels and the viewer's source-ordered playback map are retained. Real Chromium verifies the line follows playback and pointer seeking shows the correct frame; keyboard and cropped-series actions are covered in `test_a_series_follows_the_frame_played.py`. C/D verified. |
| F — Playback counter — Molecular Viewer toolbar and overlay | A live snapshot has a `sequence` (playable snapshot identity), `frame_index` (simulation step) and recorded `simulation_time_ns`; the viewer reads only the loaded frame's position. | Viewer-only control; not an analysis export. | Fixed: before, 200 archived snapshots with steps 163000–362000 returned those steps as `frame_indices`, and the toolbar showed a step against total frames. Now `frame_indices` is the source snapshot sequence, `simulation_steps` and recorded time are separate, and every playback control/readout uses the same zero-based loaded index with a user counter from 1. Real Chromium verified `1 / 200` through `200 / 200`, `step 362,000`, `time 0.724 ns`, slider seek, previous/next, first/last, play/pause, speed, loop, engine bound 199, one-frame playback and empty transport. Playback/frame API/exact-navigation suite: 24 passed; isolated playback suite: 4 passed; latest exact-frame/navigation group: 46 passed. |
| G — Agent Settings — provider authentication and model routing | Saved provider, model, and auth method determine the actual request route and subtitle. API-key fields appear only for API-key methods; OpenAI SIWC is offered for eligible open-source clients. Anthropic Claude.ai subscription sign-in is unavailable; Console OAuth is used only through the official `ant` CLI. | Credential material stays in the OS keyring or the provider CLI's app-owned profile, never ordinary settings or frontend assets. | Official docs: [OpenAI SIWC quickstart](https://developers.openai.com/siwc/quickstart), [OpenAI token sharing and sign-in](https://developers.openai.com/siwc/token-sharing-open-source/sign-in), [OpenAI models and inference](https://developers.openai.com/siwc/token-sharing-open-source/models-and-inference), [OpenAI preview limitations](https://developers.openai.com/siwc/token-sharing-open-source/preview-limitations), [OpenAI API authentication](https://developers.openai.com/api/reference/overview), [Anthropic API authentication](https://platform.claude.com/docs/en/manage-claude/authentication), [Anthropic API keys](https://platform.claude.com/docs/en/get-api-key), [Anthropic CLI authentication](https://platform.claude.com/docs/en/cli-sdks-libraries/cli/authentication), [Anthropic CLI scripting](https://platform.claude.com/docs/en/cli-sdks-libraries/cli/scripting), and [Anthropic API credits](https://platform.claude.com/docs/en/about-claude/api-credits-for-subscribers). Auth suite rerun: 17 passed, including the live browser test; auth plus movie-cancellation rerun: 18 passed. One earlier combined run hit a transient Windows `Permission denied` during the model save; it passed alone and on the full combined rerun. No user provider login or consent was performed. |
| B/D — RMSF profile on Analysis | `analysis/rmsf/rmsf.dat`, saved `options.json`, and the topology actually analyzed determine atom/residue mode and identity. Residue and atom points resolve to explicit chain/residue/insertion/atom identity before actionable; ambiguous mappings stay value-only. The graph action updates shared `STATE.selection` and Agent `current_view`. | Keep original PNG/SVG and report figures unchanged. | B verified by 55 focused tests and real Chromium/WebGL atom and residue selection, frame preservation, Agent context and clearing. Same-number residues on different chains stay value-only unless the saved identity includes chain. D ranges apply/reset without changing original point identity or full-analysis statistics; the C/D browser group is 64 passed plus the Windows encoding-only receipt rerun noted above. |
| D — Temporary ranges on time series and profiles | Range state lives in page memory; empty, reversed and out-of-bounds inputs receive visible errors, reset restores full span, and cropped points retain source identity. Run changes reset the limits. | Analysis files, full-analysis statistics and exported PNG/SVG remain unchanged. | Real Chromium coverage in `test_a_series_follows_the_frame_played.py` applies invalid/valid ranges, seeks a cropped point to its original frame, resets and verifies run-change reset. C/D focused group passed all range-related cases. |
| E — Optional caption in viewer movies | A bounded plain-text caption composes in the existing canvas overlay with recorded-time and legend overlays; empty caption preserves the default. | Movie-only option; static figure exports are unchanged. | `test_a_movie_is_made_in_the_viewer.py::test_a_movie_caption_wraps_inside_the_frame_with_time_kept` decoded the real captioned movie and checked overlay bounds at 1280×720, 1920×1080 and 3840×2160. Caption test plus concurrent movie-writer race regression: 2 passed. |
| A/C — Secondary-structure map | `analysis/ss/ss.dat`: one state string per frame plus residue records. Actionable cells route frame jumps through the shared exact-frame path and focus a residue only where both identities are known. | Keep original PNG/SVG. | C browser coverage verifies exact frame and coordinate result, pointer and keyboard paths; unavailable/aggregate identities remain value-only. |
| A/C — Cluster-label timelines (KMeans, hierarchical, DBSCAN) | Saved label rows and analysis clock retain each represented source frame; the chart uses the same exact-frame route as time series. | Keep original PNG/SVG. | C browser tests verify the selected point's displayed molecular coordinates and keyboard equivalent; no nearest sampled substitute is used. |
| A/C — Cluster-population bars | Saved cluster counts may include a `medoid_frame`; only that saved frame is a defensible representative. | Keep original PNG/SVG. | C routes a bar only to its saved medoid and labels the action representative; no inferred medoid is added. Figure coverage tests pass. |
| C — Cluster dendrogram | Linkage branches aggregate multiple frames and carry no unique source frame. | Keep original PNG/SVG. | Verified as value/source inspection only; no representative frame is invented for an aggregate branch. |
| C — RMSD frame matrix | `cluster_rmsd_matrix.npz` may be block-averaged to fit the display; an aggregated cell represents ranges, not a unique pair. | Keep original PNG/SVG. | Verified singleton cells offer explicit row/column frame actions; aggregated cells retain represented ranges and value inspection. |
| A/C — PCA/MDS/t-SNE projection points | Saved projection points align with `data.frames` and the saved analysis clock. | Keep original PNG/SVG. | C routes actionable points through A, with keyboard access and source frame identity retained; sparse and null data retain their fallback. |
| C — Free-energy landscapes | Binned projection density/free energy; a cell has no inherent frame. | Keep original PNG/SVG. | Verified as aggregate/value and source-data behavior; no frame is guessed. |
| C — Backbone-dihedral density (Ramachandran) | Binned phi/psi density; a cell has no inherent frame. | Keep original PNG/SVG. | Verified as aggregate/value and source-data behavior; no frame is guessed. |
| C — Thermodynamic observables | Aggregate bars and statistics, with no single frame or residue identity. | Keep original PNG/SVG. | Verified as value inspection; the chart does not claim a frame or residue. |
| A/C — Moments of inertia | Per-frame principal moments and the analysis clock. | Keep original PNG/SVG. | C routes frame-bearing points through A and supports keyboard access; data and export behavior stay intact. |
| B/C — Protein-ligand contact-frequency bars | `pl_contacts_per_residue.csv` identifies residues only when saved chain/residue identity matches the analyzed topology. | Keep original PNG/SVG. | B/C tests verify mapped residues select correctly; absent or ambiguous identities remain value-only. |
| B/C — Protein-ligand interaction-type bars | `pl_interactions.dat` may include aggregate-by-type values; optional by-residue and contributing-frame records establish identity. | Keep original PNG/SVG. | Verified residue actions use topology-matched saved identities; aggregate types remain value-only absent provenance. |
| C — Report page and other static figure instances | `report-page.js` renders report HTML/images and provenance chips; no dashboard data-chart instance is present there. | Static images and PNG/SVG downloads remain the report contract. | C preserves static report behavior; representative report image links load in the browser and return PNGs. |
| A–D — Empty, null, thinned, large, resized, theme-switched, and study-switched charts | Actionable points retain original source identity; missing values and aggregates do not create guessed actions. | Original exported image/data stay complete and unchanged by display filters. | A–D regression group covers exact frame, playback mapping, keyboard/pointer paths, temporary filters, null/aggregate fallback, and study/run reset. 64 passed; one Windows default-encoding receipt assertion passed on isolated `-X utf8` rerun. |
| A–D — Export regression gate (PNG/SVG, SVG archive, report images, slides, study bundle) | Exports continue from saved static analysis figures; `report/analysis_summary.svg` may contain embedded raster panels. | Keep PNG/SVG and static report/export behavior. | Representative checks: 12 passed, 1 skipped because `pypdf` is unavailable. PNG preview/SVG download, SVG archive, report image references/browser loading, slide-outline generation and project bundle passed. |
| Full-suite gate | Repository-wide tests include expensive setup/science audits and Windows-specific path/symlink assertions. | No export behavior changes. | Full suite was started but stopped after a failure at 13%; the prior verbose array diff stalled pytest. The failing `tests/test_a_preparation_can_be_repeated.py::test_the_velocities_repeat_in_processes_of_their_own` was rerun alone with plain assertions and failed in 112.50 s because process-generated velocities differ. This preparation reproducibility audit is outside dashboard scope and was already recorded as unresolved; no scientific code was changed. Other previously isolated full-suite issues are Windows path expectations/symlink privilege and an absent `pypdf`; see the execution report before delivery. |

### Local issue drafts — tracker unavailable; local implementation authorized

These drafts document the repository's issue-first step. No matching open issue was found in `princeote/FastMDXplora` or `aai-research-lab/FastMDXplora`; GitHub returned HTTP 410 because Issues is disabled in the approved fork. The user authorized local implementation on the existing branch. No issue has been submitted, and the issue tracker remains an external limitation.

1. **A — Open the exact source frame selected by a dashboard plot.** A reproducible browser case has an 11-frame source trajectory and playback sample `[0, 5, 10]`. Selecting source frame 7 currently opens source frame 5 because the chart substitutes the nearest playback index and the viewer ignores `trajectoryFrame`. Carry source identity through the chart, reuse an already-loaded exact frame, otherwise load one frame or a bounded segment through the existing guarded routes, and keep time/labels/playback consistent. Validate coordinates numerically for frame 0, final frame, nonzero start/stride, absent sampled frame, missing trajectory, interrupted request, rapid clicks, and study switch. Preserve topology filtering/alignment and keep transfer bounded.
2. **B — Preserve RMSF atom/residue identity in viewer selection and Agent context.** Read granularity and identity from saved RMSF options and the topology actually analyzed. Map atoms and residues safely, preserving chain and insertion-code distinctions; do not guess when identity is ambiguous or absent. Route graph focus through the existing shared selection mechanism so highlight, controls, `STATE.selection`, and the next Agent `current_view` agree. Clearing or changing study clears stale context; RMSF focus keeps the current frame.
3. **C — Complete accessible click behavior for data-backed dashboard figures.** Route frame-bearing chart points through exact navigation and identified residue actions through the shared selection path. Cover secondary structure, cluster timelines and saved medoids, projections, inertia, RMSD matrix row/column choices, and residue contacts/interactions. Explain aggregates (including block-averaged matrix cells) without inventing one frame. Give each pointer action an equivalent keyboard action, retain value/source behavior for aggregate figures, and preserve static image fallback. Verify rendered viewer state in a real browser across resize, theme, sparse/null data, and study changes.
4. **D — Add temporary From/To/Reset ranges to suitable dashboard series and profiles.** Keep range state in page memory, validate empty/reversed/out-of-bounds values, preserve each retained point's original frame/residue identity, and reset on study change. Keep full-analysis statistics labeled as full-analysis and leave analysis files and static PNG/SVG exports unchanged. Verify pointer and keyboard use, then click a cropped point and confirm its original source identity is selected.
5. **E — Add an optional bounded plain-text caption to viewer movies.** Reuse the existing canvas overlay and encoder, wrap text within supported output sizes, compose with existing time/legend overlays, and preserve default empty-caption behavior. Verify both captioned and uncaptioned output plus cancellation/error restoration. Keep this last and separable from graph navigation.
6. **F — Keep playback transport, slider and all frame counters on one loaded-frame index.** Reproduce a 200-frame trajectory whose archived simulation steps reached 362,000 while the toolbar compared that step to 200 frames. Use a 1-based user counter over the actual playable frames, preserve simulation step and metadata-derived time separately, and reconcile state on trajectory changes, empty data and a single frame. Verify every control and both overlay/toolbar in a live dashboard.
7. **G — Add only officially documented provider authentication to Agent Settings.** For eligible open-source clients, use OpenAI's Sign in with ChatGPT flow and supported public Responses route. For Anthropic, use API keys and the official app-owned `ant` Console OAuth profile only within its documented local scripting scope; keep Claude.ai subscription sign-in unavailable. Keep credentials out of settings/frontend/logs, route requests by the saved provider/method/model, populate models from documented catalogs, and test cancellation, failure, expiry and disconnect. No login or registration is submitted by the app developer during implementation.

## 4. PR A — Exact graph-to-frame navigation

Files: gui/series.py, gui/figure_data.py, gui/static/series-chart.js, gui/static/figure-chart.js, gui/static/molecule-viewer.js; inspect gui/trajectory_frames.py, gui/static/viewer-engine.js and existing server routes before changing them.

Carry the authoritative source trajectory frame with each selectable point. Keep analysis-row index, source frame and browser playback index distinct. Resolve time from recorded timing; when no clock exists, show frame numbers. Never infer a source frame from the downsampled display index.

Route all frame-bearing charts through one existing navigation entry point. If the frame is loaded, select it. Otherwise use the existing frame/trajectory machinery to load that one frame or a bounded segment. Add a narrowly scoped request only if existing routes cannot do this. Do not increase global playback limits or preload the entire trajectory.

Reuse topology filtering, periodic-boundary handling and alignment so an exact frame matches the ordinary viewer. Define behavior for smoothing explicitly: show the actual frame for a direct jump and clearly indicate any changed smoothing mode. Keep frame/time labels, Agent hints and subsequent play/step operations consistent.

Handle rapid clicks and study changes with request identity checks so a late response cannot show the wrong study/frame. Keep loading and unavailable states visible; never silently substitute a nearby frame.

Done when: a point whose source frame is absent from sampled playback opens coordinates matching that exact source frame; start/stride and recorded-time tests pass; rapid clicks settle on the last selection; missing trajectories produce a useful message; memory stays bounded.

## 5. PR B — Correct RMSF and selection context

Files: gui/series.py, gui/static/series-chart.js, gui/static/molecule-viewer.js; reuse gui/selection.py, gui/viewer_selections.py and viewer-engine.js. Change agent/tools.py only if an existing view-hint field cannot express the verified selection.

Read RMSF granularity from saved options: atom RMSF must be labeled and mapped as atoms; residue RMSF must preserve chain and residue identity, including insertion codes where present. Resolve against the topology actually analyzed and the atoms the viewer displays. A missing or ambiguous mapping gets a clear unavailable action, not a guessed match.

Send graph selections through the existing selection update path so molecular highlighting, selection controls and Agent current_view agree. Avoid a second graph-only selection state. A residue profile does not identify a time: focus that residue on the current frame without inventing a frame jump.

Done when: per-atom/per-residue fixtures identify the right atoms; same-number residues on different chains do not collide; cropped/filtered topology and unavailable identities are handled; graph selection appears in the next Agent context; clearing/switching studies clears stale context.

## 6. PR C — Complete meaningful dashboard click coverage

Files: gui/figure_data.py, gui/static/figure-chart.js and existing dashboard card wiring. Reuse the exact-frame navigation from PR A and selection path from PR B. Add gui/static/report-page.js wiring only for dashboard Report-page instances where the same data-backed chart can be reused.

Time series, inertia, per-frame cluster labels and projection points: open the exact recorded source frame. Secondary-structure cells: open their frame and focus their residue when both identities are known. Cluster population bars: open the saved medoid only when present and label it as a representative frame.

RMSD matrix cells represent two frames: offer explicit row/column frame actions. Block-averaged cells represent frame ranges, so do not pretend a cell is one exact pair. Residue contact/interaction bars focus identified residues; offer frame navigation only if saved records identify the contributing frames.

Landscapes, density maps, thermodynamic summaries and dendrogram aggregate branches have no single exact frame. Show their value/meaning and link to existing source data; add a representative-frame action only when provenance supplies a defined representative. Do not manufacture a frame by guessing.

Support pointer selection and an equivalent keyboard action for actionable points. Reuse existing focus/tooltips and accessible labels. Preserve the original PNG/SVG download and fallback when numbers are missing. Single-point charts, null values, large/downsampled results, resizing and theme changes must remain usable.

Done when: every inventory row has either verified frame/selection behavior or an explicit, justified aggregate/static behavior; supported actions work in a real browser; already-working charts retain their behavior.

## 7. PR D — Temporary graph ranges

Files: gui/static/series-chart.js and the existing dashboard styles; extend figure-chart.js only for chart types sharing the same meaningful axis. Add From/To/Reset controls for time/frame series and appropriate profile ranges.

Keep limits in page memory and reset on study changes. Validate reversed/out-of-bounds/empty ranges. Keep source-frame identity attached after cropping so PR A still opens the correct frame. Full-analysis statistics stay labeled as full-analysis statistics; exports remain the complete original PNG/SVG.

Done when: filtering/reset works by pointer and keyboard; filtered points navigate correctly; changing studies clears the range; data files and analysis statistics are unchanged. No bookmark files or persistent range storage.

## 8. PR E — Optional movie caption

Files: gui/static/viewer-movie.js and its existing form. Reuse the current canvas overlay and movie encoder; touch movies.py/movie_maker.py only if required by the existing option flow.

Add an optional plain-text caption with a practical length limit and predictable wrapping inside the output bounds. Compose it with the existing time and legend overlays. Default empty preserves existing movies.

Done when: captioned and uncaptioned movies render correctly at supported sizes; cancellation restores the viewer; no duplicate export backend or encoder dependency is introduced. This is lower priority than all graph navigation work.

## 9. Static exports are a regression gate

Reports, PDF output, slides and other exports continue using PNG/SVG figures. Do not embed dashboard JavaScript or interactive charts in those deliverables. Keep original figure/data/provenance links available in the dashboard.

Verify existing PNG/SVG downloads, SVG archive, report image references, slides and study bundle on a representative study. Only change analysis/plotting.py, report/document.py, report/slides.py, report/bundle.py or download routes if a demonstrated defect blocks this contract.

Do not rename an embedded bitmap as a vector figure or require all artifacts to become native editable files. The user's clarified requirement is dashboard interaction with static report exports.

## 10. Execution, tests and review

Order: baseline → A → B → C → D → E. Keep one cohesive behavior and its tests per commit. A/B may share one small review only if identity handling genuinely cannot be separated; otherwise review independently. Work remains on the requested princeote branch until a specific delivery requires another topic branch.

After each behavior, run its focused Python and real-browser tests once. Exercise actual click → viewer coordinates/frame label → Agent context, not just emitted events or source strings. Reuse fixtures and avoid broad reruns until integration or a new failure warrants them.

Before each PR: run ruff check src tests --select F,B, repository JavaScript parse checks and the required full pytest suite. Run pytest -m network tests/validation if setup, chemistry or scientific guardrails are touched. Visually inspect light/dark dashboard behavior and exported PNG/SVG/report pages; verify a built-package smoke test can load the changed assets.

Prior audit caveats are not new feature regressions: an Agent receipt test needed UTF-8 mode on Windows; the preparation audit on/off byte-equivalence test failed and its cause remained unresolved. Re-establish these on the actual implementation baseline, report them honestly and do not alter science to silence them.

Measure exact-frame navigation with the same small and larger fixtures: requested frame, shown frame, transfer size and elapsed time. Acceptance requires correct identity, bounded loading, responsive controls and no whole-trajectory download; record measured results rather than promise an unbenchmarked speed.

Each review includes the user-visible before/after behavior, a short demonstration, focused/full test results, export checks, known limitations and its revertable commits. Inspect the final diff for dead code, duplicate parsing, unnecessary endpoints/dependencies and unrelated edits.

Completion: all scoped inventory rows and acceptance tests are resolved; exact frame/selection identity survives dashboard navigation and Agent context; static exports work; optional changes remain separate; contribution history is updated only after implementation and acceptance are verified.

## Source anchors

https://github.com/princeote/FastMDXplora/tree/dashboard-targeted-improvements

https://github.com/aai-research-lab/FastMDXplora/blob/a3b54aac99dc7cd13170360aa982b136c9f50570/src/fastmdxplora/gui/static/series-chart.js

https://github.com/aai-research-lab/FastMDXplora/blob/a3b54aac99dc7cd13170360aa982b136c9f50570/src/fastmdxplora/gui/figure_data.py

https://github.com/aai-research-lab/FastMDXplora/blob/a3b54aac99dc7cd13170360aa982b136c9f50570/CONTRIBUTING.md
