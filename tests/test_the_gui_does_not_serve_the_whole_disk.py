"""The Run page's output folder is also the artifact server's root.

`launch_from_config` accepts an absolute path on purpose -- browsing to a
folder puts one in the box, and slugging it produced a folder called
`Users_someone_work` inside the launch directory, which is worse. But
`_spawn` then makes that directory `active_root`, and `/artifacts/<path>`
serves files from `active_root`. With no emptiness check, an
unauthenticated POST naming any directory made every file in it readable
over HTTP: an SSH key, a credentials file, whatever was there.

`launch` and `launch_existing_config` both required the folder to be empty
already, for the unrelated reason that a run must not clobber a previous
one. This is the third caller and it did not.

The guard inside `_send_artifact` is correct and is not what these tests
exercise -- it confines a path *within* the root. The root was the part
being chosen.
"""

from __future__ import annotations

from pathlib import Path

import pytest

from fastmdxplora.gui.exploration import DashboardRuntime


@pytest.fixture()
def runtime(tmp_path: Path) -> DashboardRuntime:
    workspace = tmp_path / "workspace"
    launches = tmp_path / "launches"
    workspace.mkdir()
    launches.mkdir()
    return DashboardRuntime(workspace_root=workspace, exploration_root=launches)


@pytest.fixture()
def somebodys_secrets(tmp_path: Path) -> Path:
    """A directory of the kind the exploit named: real files, not a run."""
    secrets = tmp_path / "home" / ".ssh"
    secrets.mkdir(parents=True)
    (secrets / "id_rsa").write_text(
        "-----BEGIN OPENSSH PRIVATE KEY-----\nnot-a-real-key\n",
        encoding="utf-8")
    return secrets


class TestAnOutputFolderMustBeEmpty:

    def test_a_directory_with_files_in_it_is_refused(
        self, runtime, somebodys_secrets
    ) -> None:
        result = runtime.launch_from_config(
            {"system": "1L2Y", "output": str(somebodys_secrets)})

        assert result["ok"] is False
        assert "not empty" in result["error"]

    def test_the_refusal_does_not_make_it_the_served_root(
        self, runtime, somebodys_secrets
    ) -> None:
        """The whole point: a refused launch must not leave the artifact
        server pointed at somebody's home directory."""
        runtime.launch_from_config(
            {"system": "1L2Y", "output": str(somebodys_secrets)})

        assert runtime.active_root != somebodys_secrets
        assert somebodys_secrets not in Path(str(runtime.data_root())).parents
        assert Path(str(runtime.data_root())) != somebodys_secrets

    def test_nothing_in_it_is_touched(self, runtime, somebodys_secrets) -> None:
        before = (somebodys_secrets / "id_rsa").read_text()

        runtime.launch_from_config(
            {"system": "1L2Y", "output": str(somebodys_secrets)})

        assert sorted(p.name for p in somebodys_secrets.iterdir()) == ["id_rsa"]
        assert (somebodys_secrets / "id_rsa").read_text() == before

    def test_an_absolute_path_to_an_empty_folder_is_still_allowed(
        self, runtime, tmp_path
    ) -> None:
        """The feature this guard must not break. Browsing to a folder puts
        an absolute path in the box, and that has to keep working."""
        chosen = tmp_path / "somewhere" / "else"
        chosen.mkdir(parents=True)

        result = runtime.launch_from_config(
            {"system": "1L2Y", "output": str(chosen)})

        # It gets past the folder check; whether it can then launch depends
        # on the chemistry backends, which is a different question.
        assert "not empty" not in str(result.get("error") or "")

    def test_a_folder_that_does_not_exist_yet_is_allowed(
        self, runtime, tmp_path
    ) -> None:
        result = runtime.launch_from_config(
            {"system": "1L2Y", "output": str(tmp_path / "brand" / "new")})

        assert "not empty" not in str(result.get("error") or "")


class TestTheSiblingsAlreadyDidThis:
    """Stated as a test so the three callers cannot drift apart again."""

    def test_launch_existing_config_refuses_a_full_folder(
        self, runtime, somebodys_secrets, tmp_path
    ) -> None:
        config = tmp_path / "study.yml"
        config.write_text("systems:\n  - system: 1L2Y\n", encoding="utf-8")

        result = runtime.launch_existing_config(
            str(config), output=str(somebodys_secrets))

        assert result["ok"] is False
        assert "not empty" in result["error"]
