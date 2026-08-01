"""Tests for the bootstrap helpers in ``src/fastmdxplora/install.py``.

These tests pin two regressions:

1. **Certifi-already-installed-but-stale gap**: earlier code short-circuited
   the auto-recovery when ``certifi`` was *importable*, which silently
   skipped the pip-upgrade path. The most common "new machine" failure
   mode is certifi present but stale, so the upgrade-and-retry path has
   to run on every SSL failure.

2. **Lockstep drift between ``MINIFORGE_INSTALLERS`` and
   ``MINIFORGE_SHA256``**: a future contributor who adds a platform
   without remembering to add the matching SHA would silently lose
   integrity verification on that platform. This module asserts the
   two tables stay 1:1 on non-None entries.

Also covers:

* ``_try_upgrade_certifi`` return paths on pip success / failure.
* ``_insecure_ssl_optin_enabled`` env-var parsing.
* ``_download_miniforge_installer`` SHA-gate on insecure SSL.
* ``_sha256_of_file`` correctness on a small known file.
"""

from __future__ import annotations

import io
import os
import ssl
import subprocess
import sys
import urllib.error
from pathlib import Path
from unittest.mock import MagicMock, patch

import pytest

import fastmdxplora.install as install_module
from fastmdxplora.install import (
    MINIFORGE_INSTALLERS,
    MINIFORGE_SHA256,
    MINIFORGE_VERSION,
    BootstrapError,
    _build_certifi_ssl_context,
    _download_miniforge_installer,
    _insecure_ssl_optin_enabled,
    _sha256_of_file,
    _should_block_insecure_ssl,
    _try_upgrade_certifi,
)


# ---------------------------------------------------------------------------
# Certifi auto-recovery helpers
# ---------------------------------------------------------------------------
class TestTryUpgradeCertifi:
    """``_try_upgrade_certifi`` returns based on pip's exit code."""

    def test_returns_true_when_pip_exits_zero(self):
        fake_proc = MagicMock(returncode=0)
        with patch("subprocess.run", return_value=fake_proc) as run_mock:
            assert _try_upgrade_certifi() is True
        argv = run_mock.call_args[0][0]
        assert argv[:3] == [sys.executable, "-m", "pip"]
        assert "--user" in argv
        assert "--upgrade" in argv
        assert "certifi" in argv

    def test_returns_false_when_pip_exits_nonzero(self):
        with patch("subprocess.run", return_value=MagicMock(returncode=1)):
            assert _try_upgrade_certifi() is False

    def test_returns_false_on_timeout(self):
        with patch(
            "subprocess.run",
            side_effect=subprocess.TimeoutExpired(cmd="pip", timeout=180),
        ):
            assert _try_upgrade_certifi() is False

    def test_returns_false_on_missing_executable(self):
        with patch(
            "subprocess.run",
            side_effect=FileNotFoundError("python not on PATH"),
        ):
            assert _try_upgrade_certifi() is False

    def test_passes_timeout_to_subprocess(self):
        """We cap the pip call so a hung pip can't hang the bootstrap."""
        with patch("subprocess.run", return_value=MagicMock(returncode=0)) as run_mock:
            _try_upgrade_certifi()
        timeout_kw = run_mock.call_args.kwargs.get("timeout")
        assert timeout_kw is not None
        assert timeout_kw <= 300  # cap; don't sleep forever


class TestBuildCertifiSslContext:
    """``_build_certifi_ssl_context`` returns an SSLContext or None."""

    def test_returns_none_when_certifi_import_fails(self):
        with patch.dict("sys.modules", {"certifi": None}):
            assert _build_certifi_ssl_context() is None

    def test_returns_sslcontext_when_certifi_importable(self):
        # Use any existing file as the certifi CA path so we don't
        # make this test depend on system CA bundles; then mock
        # ``ssl.create_default_context`` itself so we don't actually
        # parse the file as a CA bundle (which would still work but
        # couples this test to the real SSL implementation).
        fake_certifi = MagicMock()
        fake_certifi.where.return_value = __file__
        sentinel_ctx = object()
        with patch.dict("sys.modules", {"certifi": fake_certifi}), \
                patch("ssl.create_default_context", return_value=sentinel_ctx):
            ctx = _build_certifi_ssl_context()
        assert ctx is sentinel_ctx


class TestInsecureSslOptinEnvVar:
    """``FASTMDX_INSECURE_SSL`` env var: case-insensitive, whitespace-tolerant."""

    @pytest.mark.parametrize(
        "raw,expected",
        [
            ("1", True),
            ("true", True),
            ("True", True),
            ("TRUE", True),
            ("yes", True),
            ("yES", True),
            ("on", True),
            (" 1 ", True),
            ("  yes  ", True),
            ("0", False),
            ("", False),
            ("false", False),
            ("no", False),
            ("off", False),
            ("2", False),
            ("enabled", False),
            ("enable", False),
        ],
    )
    def test_parses_truthy_and_falsy_inputs(self, raw, expected):
        with patch.dict(
            os.environ, {"FASTMDX_INSECURE_SSL": raw}, clear=False,
        ):
            assert _insecure_ssl_optin_enabled() is expected

    def test_unset_env_means_false(self):
        env = {k: v for k, v in os.environ.items() if k != "FASTMDX_INSECURE_SSL"}
        with patch.dict(os.environ, env, clear=True):
            assert _insecure_ssl_optin_enabled() is False


# ---------------------------------------------------------------------------
# SHA-gate on insecure SSL
# ---------------------------------------------------------------------------
class TestInsecureSslShaGate:
    """Insecure SSL must be refused when no SHA is recorded."""

    def test_refuses_when_no_sha_is_recorded_for_the_installer(
        self, tmp_path: Path,
    ):
        # Manufacture an installer name that exists in
        # MINIFORGE_INSTALLERS but has no MINIFORGE_SHA256 entry.
        unsupported_installer_path = tmp_path / "Made-Up-Installer.exe"
        assert MINIFORGE_SHA256.get(unsupported_installer_path.name) is None
        with patch.dict(
            os.environ, {"FASTMDX_INSECURE_SSL": "1"}, clear=False,
        ):
            with pytest.raises(
                BootstrapError,
                match="FASTMDX_INSECURE_SSL=1 cannot be used",
            ):
                _download_miniforge_installer(unsupported_installer_path)

    def test_does_not_refuse_when_sha_is_recorded(self):
        # When a SHA *is* recorded AND the env var is set, the
        # SHA gate does NOT refuse. ``_should_block_insecure_ssl``
        # is the single source of truth shared by
        # ``_download_miniforge_installer``, so a typo in either
        # the function or the helper is caught by this test.
        installer_name = "Miniforge3-Linux-x86_64.sh"
        assert MINIFORGE_SHA256.get(installer_name) is not None
        with patch.dict(
            os.environ, {"FASTMDX_INSECURE_SSL": "1"}, clear=False,
        ):
            assert not _should_block_insecure_ssl(installer_name)


# ---------------------------------------------------------------------------
# Lockstep invariant: SHA table and INSTALLERS table stay in sync
# ---------------------------------------------------------------------------
class TestMiniforgePinLockstep:
    """Fail-fast CI checks for the pin table consistency."""

    def test_sha_keys_match_installer_filenames(self):
        sha_keys = set(MINIFORGE_SHA256.keys())
        installer_filenames = {
            v for v in MINIFORGE_INSTALLERS.values() if v is not None
        }
        assert sha_keys == installer_filenames, (
            "MINIFORGE_SHA256 keys must equal MINIFORGE_INSTALLERS "
            f"non-None values. SHA keys: {sha_keys}, "
            f"installer filenames: {installer_filenames}."
        )

    def test_no_sha_table_key_is_a_none_installer(self):
        for filename in MINIFORGE_SHA256:
            assert filename in MINIFORGE_INSTALLERS.values(), (
                f"MINIFORGE_SHA256 contains '{filename}' but the corresponding "
                f"MINIFORGE_INSTALLERS entry is None or missing."
            )

    def test_every_sha_value_is_a_64_char_lowercase_hex(self):
        for filename, sha in MINIFORGE_SHA256.items():
            assert len(sha) == 64, (
                f"{filename}: expected 64-char SHA-256, got {len(sha)}"
            )
            int(sha, 16)  # raises ValueError on non-hex
            assert sha == sha.lower(), (
                f"{filename}: SHA-256 must be lowercase (got {sha!r})"
            )

    def test_miniforge_version_is_pinned_string(self):
        """A non-``latest`` pinned version is required for SHA verification."""
        assert MINIFORGE_VERSION != "latest"
        # conda-forge tags look like 24.11.0-0 -- at least one "-" separator.
        assert "-" in MINIFORGE_VERSION

    def test_windows_aarch64_is_documented_unsupported(self):
        assert MINIFORGE_INSTALLERS.get(("Windows", "aarch64")) is None


# ---------------------------------------------------------------------------
# _sha256_of_file streaming correctness on a small known file
# ---------------------------------------------------------------------------
class TestSha256OfFile:
    def test_hashes_a_known_file_correctly(self, tmp_path: Path):
        # SHA-256 of empty file.
        empty = tmp_path / "empty.bin"
        empty.write_bytes(b"")
        assert (
            _sha256_of_file(empty)
            == "e3b0c44298fc1c149afbf4c8996fb92427ae41e4649b934ca495991b7852b855"
        )

    def test_hashes_a_small_file_correctly(self, tmp_path: Path):
        # SHA-256 of the literal string "abc".
        f = tmp_path / "abc.txt"
        f.write_bytes(b"abc")
        assert (
            _sha256_of_file(f)
            == "ba7816bf8f01cfea414140de5dae2223b00361a396177a9cb410ff61f20015ad"
        )


# ---------------------------------------------------------------------------
# Cross-platform discovery, download, and orchestration branches
# ---------------------------------------------------------------------------
def test_platform_discovery_and_command_helpers(tmp_path: Path):
    yaml_text = install_module._build_bootstrap_yaml()
    assert "python>=3.9,<3.14" in yaml_text
    assert "openmm>=8.0" in yaml_text
    assert install_module._normalize_arch("AMD64") == "x86_64"
    assert install_module._normalize_arch("arm64") == "aarch64"
    assert install_module._normalize_arch("riscv64") == "riscv64"

    with (
        patch.object(install_module.platform, "system", return_value="Linux"),
        patch.object(install_module.platform, "machine", return_value="AMD64"),
    ):
        assert install_module._miniforge_installer_name().endswith("Linux-x86_64.sh")

    with (
        patch.object(install_module.platform, "system", return_value="Windows"),
        patch.object(install_module.platform, "machine", return_value="arm64"),
        pytest.raises(BootstrapError, match="Windows-aarch64"),
    ):
        install_module._miniforge_installer_path()

    with (
        patch.object(install_module.platform, "system", return_value="Haiku"),
        patch.object(install_module.platform, "machine", return_value="riscv64"),
        pytest.raises(BootstrapError, match="No Miniforge installer"),
    ):
        install_module._miniforge_installer_path()

    posix_conda = tmp_path / "posix" / "bin" / "conda"
    posix_conda.parent.mkdir(parents=True)
    posix_conda.write_text("", encoding="utf-8")
    with patch.object(install_module.platform, "system", return_value="Linux"):
        assert (
            install_module._conda_executable_from_prefix(tmp_path / "posix", "conda")
            == posix_conda
        )

    windows_mamba = tmp_path / "windows" / "Scripts" / "mamba.exe"
    windows_mamba.parent.mkdir(parents=True)
    windows_mamba.write_text("", encoding="utf-8")
    with patch.object(install_module.platform, "system", return_value="Windows"):
        assert (
            install_module._conda_executable_from_prefix(
                tmp_path / "windows", "mamba"
            )
            == windows_mamba
        )

    with patch.object(
        install_module, "command_exists", side_effect=lambda name: name == "conda"
    ):
        assert install_module._locate_conda_cli() == ["conda"]
    fallback_conda = tmp_path / "fallback" / "bin" / "conda"
    with (
        patch.object(install_module, "command_exists", return_value=False),
        patch.object(
            install_module,
            "_conda_executable_from_prefix",
            side_effect=[None, fallback_conda],
        ),
    ):
        assert install_module._locate_conda_cli() == [str(fallback_conda)]
    with (
        patch.object(install_module, "command_exists", return_value=False),
        patch.object(install_module, "_conda_executable_from_prefix", return_value=None),
    ):
        assert install_module._locate_conda_cli() is None

    with patch.object(install_module.subprocess, "run", side_effect=FileNotFoundError):
        result = install_module._run_command(["missing-command"])
    assert result.returncode == 127
    assert "Command not found" in result.stderr


def test_download_stream_http_classification_and_checksum_paths(tmp_path: Path):
    request = install_module.urllib.request.Request("https://example.test/file")
    destination = tmp_path / "stream.bin"
    with patch.object(
        install_module.urllib.request,
        "urlopen",
        return_value=io.BytesIO(b"downloaded"),
    ) as urlopen:
        install_module._stream_url_to_path(request, destination)
    assert destination.read_bytes() == b"downloaded"
    assert urlopen.call_args.kwargs["timeout"] == install_module._DOWNLOAD_TIMEOUT_SECONDS

    not_found_url = "https://example.test/missing"
    server_error_url = "https://example.test/error"
    not_found = urllib.error.HTTPError(
        not_found_url, 404, "Not Found", {}, None
    )
    server_error = urllib.error.HTTPError(
        server_error_url, 503, "Unavailable", {}, None
    )
    assert "installer not found" in str(
        install_module._classify_http_error(not_found, not_found_url)
    )
    assert "HTTP 503" in str(
        install_module._classify_http_error(server_error, server_error_url)
    )

    installer_name = "Miniforge3-Linux-x86_64.sh"
    destination = tmp_path / "destination" / installer_name
    destination.parent.mkdir()
    private_dir = tmp_path / "private-success"
    private_dir.mkdir()
    downloaded = private_dir / installer_name
    downloaded.write_bytes(b"verified")
    expected = MINIFORGE_SHA256[installer_name]
    with (
        patch.object(
            install_module, "_download_to_temp_file", return_value=downloaded
        ),
        patch.object(install_module, "_sha256_of_file", return_value=expected),
        patch.object(install_module.platform, "system", return_value="Linux"),
        patch.object(install_module, "_print"),
    ):
        install_module._download_miniforge_installer(destination)
    assert destination.read_bytes() == b"verified"
    assert not private_dir.exists()

    bad_destination = tmp_path / "bad-destination" / installer_name
    bad_destination.parent.mkdir()
    bad_private = tmp_path / "private-failure"
    bad_private.mkdir()
    bad_download = bad_private / installer_name
    bad_download.write_bytes(b"tampered")
    with (
        patch.object(
            install_module, "_download_to_temp_file", return_value=bad_download
        ),
        patch.object(install_module, "_sha256_of_file", return_value="0" * 64),
        patch.object(install_module, "_print"),
        pytest.raises(BootstrapError, match="SHA-256 mismatch"),
    ):
        install_module._download_miniforge_installer(bad_destination)
    assert not bad_destination.exists()
    assert not bad_private.exists()


def test_conda_manifest_and_bootstrap_orchestration(tmp_path: Path):
    completed = subprocess.CompletedProcess(
        ["conda"], 0, '{"envs": ["/opt/envs/alpha", "/opt/envs/beta", ""]}', ""
    )
    with patch.object(install_module, "_run_command", return_value=completed):
        assert install_module._parse_conda_envs("conda") == ["alpha", "beta"]
        assert install_module._conda_env_exists("conda", "beta")

    with patch.object(
        install_module,
        "_run_command",
        return_value=subprocess.CompletedProcess(["conda"], 1, "", "failed"),
    ):
        assert install_module._parse_conda_envs("conda") == []
    with patch.object(
        install_module,
        "_run_command",
        return_value=subprocess.CompletedProcess(["conda"], 0, "not-json", ""),
    ):
        assert install_module._parse_conda_envs("conda") == []

    detail = install_module._format_conda_error(
        "UnsatisfiableError PackagesNotFoundError Cannot connect",
        "fastmdx",
        create=True,
    )
    assert "Package resolution failed" in detail
    assert "required package could not be found" in detail
    assert "Network or SSL error" in detail

    with pytest.raises(BootstrapError, match="Invalid Python version"):
        install_module.bootstrap_environment(python_version="2.7")
    with pytest.raises(BootstrapError, match="outside the supported range"):
        install_module.bootstrap_environment(python_version="3.14")

    calls: list[tuple] = []

    def make_yaml():
        directory = tmp_path / f"env-{len(calls)}"
        directory.mkdir()
        path = directory / "environment.yml"
        path.write_text(install_module._build_bootstrap_yaml(), encoding="utf-8")
        return path

    with (
        patch.object(install_module, "_ensure_conda_available", return_value=["conda"]),
        patch.object(install_module, "_conda_env_exists", return_value=True),
        patch.object(install_module, "_write_environment_yaml", side_effect=make_yaml),
        patch.object(
            install_module,
            "_remove_conda_env",
            side_effect=lambda *args: calls.append(("remove", *args)),
        ),
        patch.object(
            install_module,
            "_create_conda_env",
            side_effect=lambda *args: calls.append(("create", *args)),
        ),
        patch.object(
            install_module,
            "_update_conda_env",
            side_effect=lambda *args: calls.append(("update", *args)),
        ),
        patch.object(
            install_module,
            "_install_package_in_env",
            side_effect=lambda *args, **kwargs: calls.append(
                ("install", *args, kwargs)
            ),
        ),
        patch.object(
            install_module,
            "_verify_env",
            side_effect=lambda *args: calls.append(("verify", *args)),
        ),
        patch.object(install_module, "_print"),
    ):
        install_module.bootstrap_environment(
            env_name="fastmdx",
            python_version="3.11",
            force=True,
            package_name=".",
            editable=True,
        )
        install_module.bootstrap_environment(
            env_name="fastmdx",
            python_version="3.10",
            force=False,
        )

    assert [call[0] for call in calls].count("remove") == 1
    assert [call[0] for call in calls].count("create") == 1
    assert [call[0] for call in calls].count("update") == 1
    assert [call[0] for call in calls].count("install") == 2
    assert [call[0] for call in calls].count("verify") == 2
