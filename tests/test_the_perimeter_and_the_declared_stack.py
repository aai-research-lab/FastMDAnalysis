"""What a non-loopback dashboard will answer, and what the env file installs.

Two findings that are not about a number being wrong. The dashboard's
workflow-control endpoints were gated on a non-loopback bind and the ones
that walk the filesystem and read a named file were not, while `docs/gui.md`
said there was no network exposure at all. And `environment.yml` -- the
documented development setup -- omitted eight packages while the page
describing it said a clone set up that way "has the whole stack".
"""

from __future__ import annotations

from pathlib import Path

import pytest
import yaml

ROOT = Path(__file__).resolve().parents[1]


class TestARemoteDashboardWillNotWalkTheDisk:
    """Asked of a running server rather than read off its source.

    `_is_loopback_host` is monkeypatched to False so the handler is built the
    way a `--dashboard-host 0.0.0.0` bind builds it, while the socket stays
    on localhost and the test stays local.
    """

    @pytest.fixture()
    def remote(self, tmp_path, monkeypatch):
        import json
        import urllib.error
        import urllib.request

        from fastmdxplora.gui import server

        monkeypatch.setattr(server, "_is_loopback_host", lambda host: False)
        session = server.start_dashboard_session(
            output=tmp_path, host="127.0.0.1", port=0)
        try:
            def ask(path, body=None):
                url = f"{session.url.rstrip('/')}{path}"
                data = json.dumps(body).encode() if body is not None else None
                request = urllib.request.Request(
                    url, data=data,
                    headers={"Content-Type": "application/json"},
                    method="POST" if data is not None else "GET")
                try:
                    with urllib.request.urlopen(request, timeout=10) as reply:
                        return reply.status, reply.read().decode()
                except urllib.error.HTTPError as refused:
                    return refused.code, refused.read().decode()
            yield ask
        finally:
            session.server.shutdown()
            session.server.server_close()

    @pytest.mark.parametrize("path", [
        "/api/browse?path=/etc",
        "/api/inspect-directory?path=/etc",
    ])
    def test_it_will_not_list_a_directory(self, remote, path: str) -> None:
        status, body = remote(path)
        assert status == 403, (
            f"{path} answered {status} to a remote caller: {body[:300]}"
        )

    @pytest.mark.parametrize("path", ["/api/load-config", "/api/check-config"])
    def test_it_will_not_read_a_named_file(self, remote, path, tmp_path) -> None:
        """These quote the file back on a parse failure, which over a network
        is a file-content oracle: a planted token and an AWS key were both
        recovered through it."""
        secret = tmp_path / "creds"
        secret.write_text("AKIAIOSFODNN7EXAMPLE\n", encoding="utf-8")

        status, body = remote(path, {"path": str(secret)})

        assert status == 403
        assert "AKIA" not in body

    def test_the_control_endpoints_are_still_refused(self, remote) -> None:
        status, _ = remote("/api/run", {"system": "1UBQ"})
        assert status == 403

    def test_the_run_is_still_watchable(self, remote) -> None:
        """The narrowing has to leave what a colleague watching a job needs,
        or the flag has no use left."""
        status, _ = remote("/api/status")
        assert status == 200

    def test_the_documentation_no_longer_claims_immunity(self) -> None:
        text = (ROOT / "docs" / "gui.md").read_text(encoding="utf-8")
        assert "There is no authentication because there is no network" not in text
        assert "--dashboard-host" in text


class TestTheEnvironmentFileCarriesWhatThePageClaims:
    """`docs/installation.md` says a clone set up from it has the whole
    stack. That was false for four releases, and nothing checked it: the
    dependency test reads the conda recipe, not this file."""

    def _declared(self) -> set[str]:
        loaded = yaml.safe_load((ROOT / "environment.yml").read_text(
            encoding="utf-8"))
        names: set[str] = set()
        for entry in loaded["dependencies"]:
            if isinstance(entry, str):
                names.add(entry.split(">")[0].split("=")[0].split("<")[0].strip())
        return names

    @pytest.mark.parametrize("package", [
        "rdkit",            # ligand perception
        "propka",           # pKa assignment
        "openff-toolkit",   # ligand parameterisation
        "openmm-plumed",    # every enhanced-sampling method
        "scipy",
        "pillow",
        "netcdf4",          # Amber trajectories
        "weasyprint",       # the PDF report
    ])
    def test_the_stack_is_actually_there(self, package: str) -> None:
        assert package in self._declared(), (
            f"environment.yml does not install {package}, so the documented "
            f"development setup lacks it while docs/installation.md says "
            f"otherwise"
        )

    def test_the_conda_only_ones_are_the_reason_this_file_exists(self) -> None:
        """openmm-plumed has no PyPI distribution at all, which is why there
        is deliberately no `plumed` pip extra. A developer set up any other
        way has no umbrella, metadynamics or steered support."""
        assert "openmm-plumed" in self._declared()

    def test_the_floors_match_pyproject(self) -> None:
        """Covered generally in `test_the_declarations_agree`; asserted here
        for the two that have moved, so this file fails on its own terms."""
        text = (ROOT / "environment.yml").read_text(encoding="utf-8")
        assert "numpy>=2.0" in text
        assert "openmm>=8.2" in text


class TestTheCorpusHasSomewhereToRun:
    """AUD31. 105 network-marked tests -- the structure corpus and the
    method-produces-its-result suite -- ran in no environment at all:
    deselected by pytest.ini, absent from CI, unmentioned in CONTRIBUTING."""

    def _workflow(self) -> dict:
        return yaml.safe_load(
            (ROOT / ".github" / "workflows" / "tests.yml").read_text(
                encoding="utf-8"))

    def test_there_is_a_job_for_it(self) -> None:
        assert "corpus" in self._workflow()["jobs"]

    def test_it_runs_the_marked_tests(self) -> None:
        steps = self._workflow()["jobs"]["corpus"]["steps"]
        run = " ".join(str(step.get("run", "")) for step in steps)
        assert "-m network" in run
        assert "tests/validation" in run

    def test_it_installs_the_extras_no_job_ever_has(self) -> None:
        """`[validation]` brings MDAnalysis and ProLIF. Without it
        `cross_tool.py` -- 853 lines of independent comparison against
        pre-registered tolerances -- skips, which is how it has never run."""
        steps = self._workflow()["jobs"]["corpus"]["steps"]
        run = " ".join(str(step.get("run", "")) for step in steps)
        assert "validation" in run and "ligand" in run and "md" in run

    def test_it_is_scheduled_and_can_be_asked_for(self) -> None:
        workflow = self._workflow()
        triggers = workflow[True] if True in workflow else workflow["on"]
        assert "schedule" in triggers
        assert "workflow_dispatch" in triggers

    def test_it_does_not_gate_a_pull_request(self) -> None:
        """A corpus that fetches from RCSB will go red for reasons that are
        not the contributor's, and a required check that does that gets
        ignored or removed."""
        condition = self._workflow()["jobs"]["corpus"]["if"]
        assert "schedule" in condition or "workflow_dispatch" in condition

    def test_contributing_says_how_to_run_it(self) -> None:
        text = (ROOT / "CONTRIBUTING.md").read_text(encoding="utf-8")
        assert "-m network" in text


class TestTheLinterHasAJobAndAStatedScope:
    """AUD35. CONTRIBUTING required it before there was a job for it."""

    def test_there_is_a_lint_step(self) -> None:
        workflow = yaml.safe_load(
            (ROOT / ".github" / "workflows" / "tests.yml").read_text(
                encoding="utf-8"))
        names = [s.get("name") for s in workflow["jobs"]["test"]["steps"]]
        assert "Lint" in names

    def test_the_rules_it_gates_on_are_clean(self) -> None:
        """The gate is F and B -- the rules that find bugs. Holding those
        hostage to a backlog of import order and line length is how a lint
        job ends up permanently disabled.

        Skipped where ruff is absent rather than erroring. The first version
        of this called `subprocess.run(["ruff", ...])` unguarded, and ruff
        was in the `[dev]` extra and not in `[test]` -- so it raised
        FileNotFoundError in every environment installed the documented way,
        including CI's own test job. A test that shells out to a tool has to
        say what it needs. It is in `[test]` now as well, so the skip should
        be rare; it stays because an environment predating that change is
        not a broken one.
        """
        import shutil
        import subprocess

        if shutil.which("ruff") is None:
            pytest.skip("ruff is not installed; `pip install -e '.[test]'`")

        result = subprocess.run(
            ["ruff", "check", "src", "tests", "--select", "F,B"],
            cwd=ROOT, capture_output=True, text=True)
        assert result.returncode == 0, result.stdout[-2000:]

    def test_contributing_says_which_subset_must_pass(self) -> None:
        text = (ROOT / "CONTRIBUTING.md").read_text(encoding="utf-8")
        assert "--select F,B" in text


class TestThreeTestsThatHadNeverRun:
    """Found by the F811 rule the lint gate now enforces: a function defined
    twice in one file, so the earlier one is replaced before pytest collects
    it. Not in the 39 findings -- this is what turning the rule on found."""

    @pytest.mark.parametrize("path", [
        "tests/test_live_status.py",
        "tests/test_concrete_analyses.py",
    ])
    def test_no_test_is_defined_twice(self, path: str) -> None:
        import ast
        from collections import Counter

        tree = ast.parse((ROOT / path).read_text(encoding="utf-8"))
        repeated: list[str] = []
        for node in ast.walk(tree):
            if not isinstance(node, ast.ClassDef):
                continue
            names = Counter(
                child.name for child in node.body
                if isinstance(child, ast.FunctionDef)
                and child.name.startswith("test_"))
            repeated += [f"{node.name}.{n}" for n, c in names.items() if c > 1]

        assert not repeated, (
            f"defined twice, so the first never runs: {repeated}"
        )


class TestAnUnguardedImportIsACoreDependency:
    """AUD18, and the reason it was invisible.

    `analysis/ligand_chemistry.py` imported RDKit at the top of `_perceive`
    with no guard, and RDKit is in `[ligand]` and `[test]` and not in the
    core dependency list. So on a plain `pip install fastmdxplora` the call
    raised ModuleNotFoundError, escaping `resolve_ligand_chemistry` and its
    carefully written refusal -- the one that names every route tried and
    says to supply an SDF. The user got a traceback instead of a next step,
    on the single path where "it refuses rather than guesses" was not true.

    `test_dependencies_declared` was supposed to catch exactly this and could
    not: it sliced pyproject from `dependencies = [` to `[project.urls]`,
    which spans the whole optional-dependencies table, so anything named in
    any extra counted as a core declaration. That slice is a TOML parse now,
    and the distinction it lost is what this class asserts.
    """

    def _core(self) -> set[str]:
        from tests._toml import load_toml

        project = load_toml(ROOT / "pyproject.toml")["project"]
        import re as _re
        return {_re.split(r"[><=!~;\s\[]", entry, maxsplit=1)[0].strip().lower()
                for entry in project["dependencies"]}

    def test_rdkit_is_not_a_core_dependency(self) -> None:
        """Stated so the test below means something. RDKit belongs in the
        ligand extra -- the point is that the code must cope without it."""
        assert "rdkit" not in self._core()

    def test_perception_returns_none_rather_than_raising(self, monkeypatch) -> None:
        """The behaviour that lets the refusal downstream do its job."""
        import builtins

        from fastmdxplora.analysis import ligand_chemistry

        real_import = builtins.__import__

        def without_rdkit(name, *args, **kwargs):
            if name.startswith("rdkit"):
                raise ImportError("No module named 'rdkit'")
            return real_import(name, *args, **kwargs)

        monkeypatch.setattr(builtins, "__import__", without_rdkit)
        assert ligand_chemistry._perceive(None, [], "LIG", None) is None

    def test_the_refusal_says_rdkit_is_missing(self, monkeypatch) -> None:
        """And names the install line, rather than reading as though
        perception had been tried and had not worked."""
        from fastmdxplora.analysis import ligand_chemistry

        monkeypatch.setattr(ligand_chemistry, "_rdkit_is_absent", lambda: True)
        monkeypatch.setattr(ligand_chemistry, "_perceive",
                            lambda *a, **k: None)

        with pytest.raises(ValueError) as refusal:
            ligand_chemistry.resolve_ligand_chemistry(
                traj=None, atom_indices=[], resname="LIG")

        message = str(refusal.value)
        assert "RDKit is not installed" in message
        assert "conda install" in message
        assert "Supply an SDF" in message
