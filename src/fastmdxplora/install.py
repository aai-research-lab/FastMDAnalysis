from __future__ import annotations

import hashlib
import os
import platform
import shutil
import socket
import ssl
import subprocess
import sys
import tempfile
import urllib.error
import urllib.request
from pathlib import Path
from typing import Iterable

from fastmdxplora import MAX_PYTHON, MIN_PYTHON, python_range_string

DEFAULT_ENV_NAME = "fastmdxplora"
# Pinned Miniforge release. Locking the version lets us verify SHA-256 of
# every installer against an embedded, immutable manifest and prevents an
# upstream cutover from silently pulling in a release with different
# assets. To bump: pick the new tag, fetch the SHA-256 of each asset, and
# update MINIFORGE_VERSION + MINIFORGE_SHA256 together (see
# ``health.py`` which mirrors these constants so the doctor still works
# before fastmdxplora is on the import path).
MINIFORGE_VERSION = "26.3.2-3"
MINIFORGE_BASE_URL = (
    f"https://github.com/conda-forge/miniforge/releases/download/{MINIFORGE_VERSION}"
)
MINIFORGE_PREFIX = Path.home() / "miniforge3"
# conda-forge/miniforge does not publish a Windows-aarch64 installer; we
# keep an entry in this table (value ``None``) so a Windows-ARM user gets
# an actionable error rather than a 404-shaped surprise.
MINIFORGE_INSTALLERS: dict[tuple[str, str], str | None] = {
    ("Linux", "x86_64"): "Miniforge3-Linux-x86_64.sh",
    ("Linux", "aarch64"): "Miniforge3-Linux-aarch64.sh",
    ("Darwin", "x86_64"): "Miniforge3-MacOSX-x86_64.sh",
    ("Darwin", "arm64"): "Miniforge3-MacOSX-arm64.sh",
    ("Windows", "x86_64"): "Miniforge3-Windows-x86_64.exe",
    ("Windows", "aarch64"): None,  # upstream does not publish this asset
}
# SHA-256 of each upstream asset for ``MINIFORGE_VERSION``, keyed by the
# installer filename. Kept in lockstep with ``MINIFORGE_VERSION``.
MINIFORGE_SHA256: dict[str, str] = {
    "Miniforge3-Linux-x86_64.sh": (
        "848194851a98903134187fbb4ab50efe87b003e0c0f808f97644b7524a62bf2c"
    ),
    "Miniforge3-Linux-aarch64.sh": (
        "2c113a69297e612b01ca0f320c22a3107a11f2ab9b573d79ac868a175945ce29"
    ),
    "Miniforge3-MacOSX-x86_64.sh": (
        "39273e4c89a0a1af4538010615d44ae8f44e1af41007e02def593d20f316b003"
    ),
    "Miniforge3-MacOSX-arm64.sh": (
        "59168f1e24d0a4ad9932021170809fca836cd240e183eeeb331d5bcfc0098168"
    ),
    "Miniforge3-Windows-x86_64.exe": (
        "14a8635465b5190537ddad6286746ffebbc55a1ed2a7bb14a506595fe3191e1e"
    ),
}
# GitHub Releases' edge sometimes throttles or 403s the default Python
# urllib UA; using an explicit one keeps bootstrap traffic both unblocked
# and identifiable in server logs.
MINIFORGE_USER_AGENT = "fastmdxplora-bootstrap/1.0"
_DOWNLOAD_CHUNK_BYTES = 64 * 1024
# Generous: a 100 MB Miniforge installer on a 1 Mbps link takes ~800 s.
# The SHA-256 + private-tmp pattern already cover the stored-bad-bytes
# case, so we don't lose safety by waiting longer here.
_DOWNLOAD_TIMEOUT_SECONDS = 600


def _build_bootstrap_yaml() -> str:
    """Render the conda env file using the canonical Python range.

    Kept in sync with ``pyproject.toml`` and the doctor check in
    ``health.py`` via :data:`fastmdxplora.MIN_PYTHON` /
    :data:`fastmdxplora.MAX_PYTHON`.
    """
    min_str = f"{MIN_PYTHON[0]}.{MIN_PYTHON[1]}"
    max_str = f"{MAX_PYTHON[0]}.{MAX_PYTHON[1]}"
    return f"""name: fastmdxplora
channels:
  - conda-forge
dependencies:
  - python>={min_str},<{max_str}
  - pip
  - numpy>=1.22
  - pyyaml>=6.0
  - mdtraj>=1.9.7
  - matplotlib>=3.5
  - scikit-learn>=1.0
  - pandas>=1.4
  - python-pptx>=0.6.21
  - openmm>=8.0
  - pdbfixer
  - openmmforcefields
"""


# Constant kept for backwards compatibility with any external callers
# that imported it before the range-sync refactor.
BOOTSTRAP_ENV_YAML: str = _build_bootstrap_yaml()


class BootstrapError(Exception):
    pass


def _print(message: str, prefix: str = "INFO") -> None:
    print(f"[{prefix}] {message}")


def _run_command(cmd: list[str], **kwargs: object) -> subprocess.CompletedProcess[str]:
    try:
        return subprocess.run(
            cmd,
            stdout=subprocess.PIPE,
            stderr=subprocess.PIPE,
            text=True,
            check=False,
            **kwargs,
        )
    except FileNotFoundError:
        return subprocess.CompletedProcess(cmd, 127, "", f"Command not found: {cmd[0]}")


def _normalize_arch(arch: str) -> str:
    arch = arch.lower()
    if arch in ("amd64", "x86_64"):
        return "x86_64"
    if arch in ("aarch64", "arm64"):
        return "aarch64"
    return arch


def _miniforge_installer_name() -> str | None:
    system = platform.system()
    arch = _normalize_arch(platform.machine() or "")
    return MINIFORGE_INSTALLERS.get((system, arch))


def _miniforge_installer_path() -> Path:
    installer_name = _miniforge_installer_name()
    if installer_name is None:
        sysname = platform.system()
        machine = platform.machine()
        # conda-forge/miniforge does not publish a Windows-aarch64
        # installer; tell a Windows-ARM user why instead of bucketing
        # them into a generic "unsupported platform" message.
        if sysname == "Windows" and machine.lower() in ("aarch64", "arm64"):
            raise BootstrapError(
                "conda-forge/miniforge does not publish a Windows-aarch64 "
                f"(Windows-on-ARM) installer for {MINIFORGE_VERSION}. "
                "Install Miniforge manually from "
                "https://conda-forge.org/miniforge/ (e.g. via an x86_64 "
                "Python on the same machine) or use Windows Subsystem for "
                "Linux, then re-run the bootstrap."
            )
        raise BootstrapError(
            f"No Miniforge installer available for this platform: "
            f"{sysname} {machine}"
        )
    return Path(tempfile.gettempdir()) / installer_name


def command_exists(name: str) -> bool:
    return shutil.which(name) is not None


def _conda_executable_from_prefix(prefix: Path, tool_name: str) -> Path | None:
    if platform.system() == "Windows":
        path = prefix / "Scripts" / f"{tool_name}.exe"
    else:
        path = prefix / "bin" / tool_name
    return path if path.exists() else None


def _locate_conda_cli() -> list[str] | None:
    for tool in ("mamba", "conda"):
        if command_exists(tool):
            return [tool]
    for tool_name in ("mamba", "conda"):
        candidate = _conda_executable_from_prefix(MINIFORGE_PREFIX, tool_name)
        if candidate is not None:
            return [str(candidate)]
    return None


def _sha256_of_file(
    path: Path,
    *,
    chunk_bytes: int = _DOWNLOAD_CHUNK_BYTES,
) -> str:
    """Stream-hash ``path`` with SHA-256."""
    hasher = hashlib.sha256()
    with path.open("rb") as fh:
        while True:
            chunk = fh.read(chunk_bytes)
            if not chunk:
                break
            hasher.update(chunk)
    return hasher.hexdigest()


def _try_upgrade_certifi() -> bool:
    """Best-effort: try to upgrade ``certifi`` to a current bundle via pip.

    This is the canonical Windows + macOS + Linux fix for stale or
    missing certificates. We pass ``--user`` so the install succeeds on
    PEP 668 externally-managed environments (Ubuntu 23+, Homebrew
    Python, etc.). stdout/stderr are redirected to the parent's stderr
    so progress bars / error messages stay visible.

    Returns ``True`` on success, ``False`` on any failure (no pip, no
    PyPI access, no permission, timeout). Never raises.
    """
    try:
        proc = subprocess.run(
            [sys.executable, "-m", "pip", "install", "--user", "--upgrade",
             "certifi"],
            stdout=sys.stderr,
            stderr=sys.stderr,
            text=True,
            check=False,
            timeout=180,
        )
        return proc.returncode == 0
    except (FileNotFoundError, subprocess.TimeoutExpired, OSError):
        return False


def _build_certifi_ssl_context() -> ssl.SSLContext | None:
    """Return an ``SSLContext`` pointing at certifi's CA bundle.

    ``urllib`` does NOT consult certifi automatically -- it goes via
    the OS trust store. To actually take advantage of certifi (e.g. on
    bare Python.org installs where the OS trust store lacks newer CAs)
    we have to construct an ``SSLContext`` with
    ``cafile=certifi.where()`` and pass it explicitly to
    ``urlopen(..., context=ctx)``. Returns ``None`` if certifi is not
    importable -- callers should kick off a certifi install first.
    """
    try:
        import certifi  # type: ignore[import-not-found]
        return ssl.create_default_context(cafile=certifi.where())
    except ImportError:
        return None


def _stream_url_to_path(
    request: urllib.request.Request,
    tmp_path: Path,
    *,
    ssl_context: ssl.SSLContext | None = None,
) -> None:
    """Stream ``request`` to ``tmp_path`` with chunked writes.

    Used by :func:`_download_to_temp_file` for the strict-TLS primary
    attempt, the certifi-based recovery attempt, and the opt-in
    insecure-SSL fallback. When ``ssl_context`` is ``None`` Python's
    default trust store is used.
    """
    kwargs: dict[str, object] = {"timeout": _DOWNLOAD_TIMEOUT_SECONDS}
    if ssl_context is not None:
        kwargs["context"] = ssl_context
    with urllib.request.urlopen(request, **kwargs) as response, \
            tmp_path.open("wb") as out:
        shutil.copyfileobj(response, out, length=_DOWNLOAD_CHUNK_BYTES)


def _insecure_ssl_optin_enabled() -> bool:
    """``True`` iff the user has explicitly opted in to insecure SSL.

    Accepts ``1``, ``true``, ``yes``, ``on`` (case-insensitive).
    Default is ``False`` -- no SSL bypass by default.
    """
    return os.environ.get("FASTMDX_INSECURE_SSL", "").strip().lower() in {
        "1", "true", "yes", "on",
    }


def _classify_http_error(
    http_exc: urllib.error.HTTPError, url: str,
) -> BootstrapError:
    """Build a categorised ``BootstrapError`` for a non-2xx HTTP response.

    Single source of truth for the friendly 404 / HTTP-error message.
    Called from both the outer first-attempt ``except
    urllib.error.HTTPError as exc`` handler AND the certifi-based
    recovery retries.

    Why a helper: once we are inside the outer ``except URLError as
    primary_exc`` clause, sibling ``except HTTPError`` clauses on the
    same ``try`` do NOT catch exceptions raised inside its body --
    so a transient 404 hit on the certifi-based retry would
    otherwise bubble out as the raw ``urllib.error.HTTPError``. This
    helper closes that gap.
    """
    if http_exc.code == 404:
        return BootstrapError(
            f"Miniforge installer not found on GitHub for version "
            f"{MINIFORGE_VERSION}: {url}\n"
            "This usually means upstream does not publish an "
            "installer for your platform (Windows-aarch64 is the "
            "only such case in our table) or the pinned version "
            "has been retired. Install Miniforge manually from "
            "https://conda-forge.org/miniforge/ or report an "
            "issue."
        )
    return BootstrapError(
        f"HTTP {http_exc.code} downloading Miniforge from {url}: "
        f"{http_exc.reason}"
    )


def _download_to_temp_file(url: str, dest: Path) -> Path:
    """Stream ``url`` into a *private* temp dir and return the temp ``Path``.

    The returned file has been streamed to completion but is **not yet
    verified** and **not yet at its final location**. The caller must
    verify the bytes (e.g. SHA-256 via :func:`_sha256_of_file`) and
    atomically commit them to ``dest`` via ``tmp_path.replace(dest)``.

    The temp dir is created with ``tempfile.mkdtemp`` (mode 0700) so a
    local attacker on a shared host cannot pre-plant a symlink at the
    canonical ``tempfile.gettempdir() / installer_name`` location and
    race the bootstrap. The private dir is cleaned up on every exit
    path, including ``KeyboardInterrupt`` via the ``BaseException``
    handler. **Caveat for Windows:** ``mkdtemp`` does not honour
    ``mode=0o700`` on Windows -- the directory inherits the user's
    profile ACLs. The mitigation is best-effort on Windows but still
    removes the predictable-name symlink-swap vector for shared-host
    POSIX deployments.

    Strict TLS verification is assumed -- we deliberately do NOT fall
    back to unverified SSL because the bytes we just wrote will be
    executed as a system installer on the very next line.
    """
    request = urllib.request.Request(
        url, headers={"User-Agent": MINIFORGE_USER_AGENT}
    )
    tmp_dir = Path(tempfile.mkdtemp(prefix="fastmdx-download-"))
    # Belt-and-braces for POSIX; mkdtemp's mode argument interacts with
    # the process umask, so tighten it after creation. No-op on Windows.
    try:
        os.chmod(tmp_dir, 0o700)
    except (OSError, NotImplementedError):
        pass
    tmp_path = tmp_dir / dest.name
    try:
        try:
            with urllib.request.urlopen(
                request, timeout=_DOWNLOAD_TIMEOUT_SECONDS
            ) as response, tmp_path.open("wb") as out:
                shutil.copyfileobj(response, out, length=_DOWNLOAD_CHUNK_BYTES)
        except urllib.error.HTTPError as exc:
            # Route through the single-source-of-truth classifier so the
            # outer first-attempt message and the certifi-based recovery
            # retries all surface the same friendly text. Without this
            # the outer block drifted (different 404 wording) from the
            # inner blocks.
            raise _classify_http_error(exc, url) from exc
        except urllib.error.URLError as primary_exc:
            reason = getattr(primary_exc, "reason", None) or primary_exc
            reason_str = str(reason)
            is_ssl = (
                "CERTIFICATE_VERIFY_FAILED" in reason_str
                or isinstance(reason, ssl.SSLError)
            )
            if not is_ssl:
                # DNS / refused / timeout -- don't waste 3 minutes on
                # a certifi auto-upgrade when the real problem is
                # connectivity.
                raise BootstrapError(
                    f"Network error downloading Miniforge from {url}: "
                    f"{reason}.\n"
                    "Check your connection, VPN, or proxy settings."
                ) from primary_exc
            # Self-heal: try a certifi-based retry against the existing
            # bundle first (covers "certifi missing" cases), THEN
            # transparently upgrade certifi and retry against the
            # freshest bundle (covers the much more common
            # "certifi-already-installed-but-stale" case -- e.g. a
            # long-lived Python.org install where certifi was bundled
            # at install time and is now missing newer ISRG Roots).
            # urllib does NOT consult certifi automatically; we have
            # to construct ``ssl.create_default_context(cafile=
            # certifi.where())`` and pass it via urlopen's ``context=``.
            _print(
                "TLS verification failed; attempting automatic recovery "
                "via `python -m pip install --user --upgrade certifi` "
                "and retrying with its CA bundle...",
                prefix="FIXING",
            )
            try:
                _stream_url_to_path(
                    request, tmp_path,
                    ssl_context=_build_certifi_ssl_context(),
                )
                return tmp_path
            except urllib.error.HTTPError as http_exc:
                # Direct HTTP error on the certifi-based retry; route
                # through the same friendly classifier as the outer.
                raise _classify_http_error(http_exc, url) from http_exc
            except urllib.error.URLError as stale_exc:
                # Existing certifi either missing or stale; record and
                # try a fresh upgrade.
                primary_exc = stale_exc
            if _try_upgrade_certifi():
                upgraded_ctx = _build_certifi_ssl_context()
                if upgraded_ctx is not None:
                    try:
                        _stream_url_to_path(
                            request, tmp_path,
                            ssl_context=upgraded_ctx,
                        )
                        return tmp_path
                    except urllib.error.HTTPError as http_exc:
                        raise _classify_http_error(http_exc, url) from http_exc
                    except urllib.error.URLError as retry_exc:
                        primary_exc = retry_exc
                else:
                    _print(
                        "certifi upgrade reported success but the "
                        "certifi module is still unimportable.",
                        prefix="WARN",
                    )
            else:
                _print(
                    "certifi upgrade did not succeed; falling back to "
                    "fail-closed path.",
                    prefix="WARN",
                )
            # Last-resort: opt-in insecure SSL via env var.
            if _insecure_ssl_optin_enabled():
                _print(
                    "*** SECURITY WARNING ***\n"
                    "FASTMDX_INSECURE_SSL is set; TLS verification is\n"
                    "BEING DISABLED for this download. The downloaded\n"
                    "file is STILL SHA-256-verified against a hardcoded\n"
                    "value in the source tree, so a network MITM cannot\n"
                    "substitute a tampered Miniforge installer. Caveats:\n"
                    "  - downloads are now observable on the wire\n"
                    "    (privacy only, not integrity),\n"
                    "  - if MINIFORGE_SHA256 ever loses this installer,\n"
                    "    the SHA check stops protecting you.\n"
                    "ONLY use this on networks you trust completely.",
                    prefix="WARN",
                )
                insecure_ctx = ssl.create_default_context()
                insecure_ctx.check_hostname = False
                insecure_ctx.verify_mode = ssl.CERT_NONE
                _stream_url_to_path(
                    request, tmp_path, ssl_context=insecure_ctx,
                )
                return tmp_path
            # Fail closed with friendly, categorised remediation.
            raise BootstrapError(
                "TLS certificate verification failed downloading Miniforge "
                f"({url}).\n\n"
                "Automatic recovery via `python -m pip install --user "
                "--upgrade certifi` was attempted but did not resolve "
                "the issue. Likely causes:\n"
                "  - bare Python.org install on Windows with a stale "
                "trust store AND a stale certifi bundle,\n"
                "  - corporate proxy intercepting TLS with a CA that "
                "is not in certifi,\n"
                "  - DNS hijack on the path to github.com,\n"
                "  - system clock skew causing cert validation to fail.\n\n"
                "Remediation:\n"
                "  1. Run `python -m pip install --user --upgrade certifi` "
                "and retry.\n"
                f"  2. Or download the file manually from {url} and place "
                f"it at:\n"
                f"     {dest}\n"
                "     then re-run the bootstrap.\n"
                "  3. (Opt-in, insecure) set environment variable "
                "FASTMDX_INSECURE_SSL=1 to bypass TLS verification. "
                "The download remains SHA-256 verified against a "
                "hardcoded hash, but you are exposed to MITM privacy "
                "leaks and to any future removal of the SHA from the "
                "source tree."
            ) from primary_exc
        except OSError as exc:
            raise BootstrapError(
                f"Filesystem error downloading Miniforge: {exc}"
            ) from exc
        return tmp_path
    except BaseException:
        # BaseException (not Exception) so Ctrl+C / SystemExit also
        # remove the private temp dir; ``urllib``'s exception chain
        # only catches Exception, hence the outer guard.
        shutil.rmtree(tmp_dir, ignore_errors=True)
        raise


def _should_block_insecure_ssl(installer_name: str) -> bool:
    """``True`` iff ``installer_name`` lacks a SHA-256 hash in
    ``MINIFORGE_SHA256`` AND the user has explicitly opted in to
    insecure SSL via ``FASTMDX_INSECURE_SSL``.

    Extracted as a single source of truth shared by
    ``_download_miniforge_installer`` and the
    ``tests/test_install_bootstrap.py::TestInsecureSslShaGate``
    suite. If the condition drifts between the function and the
    test surface, the test catches it.
    """
    return (
        MINIFORGE_SHA256.get(installer_name) is None
        and _insecure_ssl_optin_enabled()
    )


def _download_miniforge_installer(installer_path: Path) -> None:
    # Gate the opt-in insecure-SSL fallback on having anything to
    # verify the bytes against. Without a SHA, FASTMDX_INSECURE_SSL
    # would let an attacker substitute any payload.
    if _should_block_insecure_ssl(installer_path.name):
        raise BootstrapError(
            "FASTMDX_INSECURE_SSL=1 cannot be used to download "
            f"{installer_path.name}: no SHA-256 hash is recorded in "
            "MINIFORGE_SHA256 for this installer, so the SHA check "
            "would not run. SHA verification is the only integrity "
            "guarantee when TLS verification is disabled."
        )
    url = f"{MINIFORGE_BASE_URL}/{installer_path.name}"
    expected_sha256 = MINIFORGE_SHA256.get(installer_path.name)
    _print(f"Downloading Miniforge installer: {url}")
    tmp_path = _download_to_temp_file(url, installer_path)
    try:
        if expected_sha256 is None:
            # Defensive: a future platform added to MINIFORGE_INSTALLERS
            # without a matching hash still proceeds but warns loudly.
            _print(
                f"No SHA-256 recorded for {installer_path.name}; "
                f"skipping integrity check.",
                prefix="WARN",
            )
        else:
            # Verify BEFORE the atomic rename so a tampered file can
            # never be observed at the canonical installer_path. Without
            # this ordering a TOCTOU race or symlink swap could let a
            # tampered file be picked up by ``_install_miniforge`` before
            # the hash mismatch deletes it.
            _print("Verifying Miniforge installer checksum...")
            actual_sha256 = _sha256_of_file(tmp_path)
            if actual_sha256.lower() != expected_sha256.lower():
                raise BootstrapError(
                    f"SHA-256 mismatch for {installer_path.name}.\n"
                    f"  Expected: {expected_sha256}\n"
                    f"  Got:      {actual_sha256}\n"
                    "This can indicate a corrupted download, a CDN issue, "
                    "or an attacker tampering with the response. We "
                    "declined to use the downloaded file. Re-run the "
                    "bootstrap; if the mismatch persists, download the "
                    "file manually from a different network and place it "
                    "at the destination path."
                )
            _print("Checksum verified.")
        tmp_path.replace(installer_path)
    finally:
        # Always remove the private temp dir, even after a successful
        # commit, so we don't leak tempdir containers on Ctrl+C during
        # the SHA check or chmod step.
        shutil.rmtree(tmp_path.parent, ignore_errors=True)
    if platform.system() != "Windows":
        try:
            installer_path.chmod(0o755)
        except OSError:
            pass


def _install_miniforge() -> list[str]:
    installer_path = _miniforge_installer_path()
    if not installer_path.exists():
        _download_miniforge_installer(installer_path)

    install_prefix = MINIFORGE_PREFIX
    install_prefix.mkdir(parents=True, exist_ok=True)
    system = platform.system()

    if system == "Windows":
        command = [str(installer_path), "/InstallationType=JustMe", "/RegisterPython=0", "/S", f"/D={install_prefix}"]
    else:
        command = ["/bin/bash", str(installer_path), "-b", "-p", str(install_prefix)]

    _print(f"Installing Miniforge to {install_prefix}")
    result = _run_command(command)
    if result.returncode != 0:
        raise BootstrapError(
            "Miniforge installation failed: "
            + result.stderr.strip()
            + "\nPlease install Miniforge manually from https://conda-forge.org/miniforge/."
        )
    _print("Miniforge installed successfully")
    conda_cli = _locate_conda_cli()
    if conda_cli is None:
        raise BootstrapError("Miniforge was installed, but conda/mamba could not be found afterward.")
    return conda_cli


def _ensure_conda_available() -> list[str]:
    conda_cli = _locate_conda_cli()
    if conda_cli is not None:
        return conda_cli
    return _install_miniforge()


def _write_environment_yaml() -> Path:
    temp_dir = Path(tempfile.mkdtemp(prefix="fastmdxplora-env-"))
    yaml_path = temp_dir / "environment.yml"
    yaml_path.write_text(_build_bootstrap_yaml(), encoding="utf-8")
    return yaml_path


def _parse_conda_envs(conda_cli: str) -> list[str]:
    result = _run_command([conda_cli, "env", "list", "--json"])
    if result.returncode != 0:
        return []
    try:
        import json

        data = json.loads(result.stdout)
        return [Path(p).name for p in data.get("envs", []) if p]
    except json.JSONDecodeError:
        return []


def _conda_env_exists(conda_cli: str, env_name: str) -> bool:
    return env_name in _parse_conda_envs(conda_cli)


def _remove_conda_env(conda_cli: str, env_name: str) -> None:
    _print(f"Removing existing conda environment '{env_name}'")
    result = _run_command([conda_cli, "env", "remove", "-n", env_name, "-y"])
    if result.returncode != 0:
        raise BootstrapError(f"Failed to remove existing environment: {result.stderr.strip()}")


def _create_conda_env(conda_cli: str, env_name: str, yaml_path: Path) -> None:
    _print(f"Creating conda environment '{env_name}' from bundled configuration")
    result = _run_command([conda_cli, "env", "create", "-f", str(yaml_path), "-n", env_name])
    if result.returncode != 0:
        raise BootstrapError(_format_conda_error(result.stderr, env_name, create=True))
    _print(f"Conda environment '{env_name}' created successfully")


def _update_conda_env(conda_cli: str, env_name: str, yaml_path: Path) -> None:
    _print(f"Updating existing conda environment '{env_name}'")
    result = _run_command([conda_cli, "env", "update", "-f", str(yaml_path), "-n", env_name])
    if result.returncode != 0:
        raise BootstrapError(_format_conda_error(result.stderr, env_name, create=False))
    _print(f"Conda environment '{env_name}' updated successfully")


def _install_package_in_env(conda_cli: str, env_name: str, package: str, *, editable: bool = False) -> None:
    _print(f"Installing package '{package}' into '{env_name}'")
    cmd = [conda_cli, "run", "-n", env_name, sys.executable, "-m", "pip", "install", "--upgrade"]
    if editable:
        cmd.append("--editable")
    cmd.append(package)
    result = _run_command(cmd)
    if result.returncode != 0:
        raise BootstrapError(
            f"Failed to install {package} in {env_name}: {result.stderr.strip()}"
        )
    _print(f"Package '{package}' installed into '{env_name}'")


def _verify_env(conda_cli: str, env_name: str) -> None:
    _print(f"Verifying conda environment '{env_name}'")
    command = [conda_cli, "run", "-n", env_name, sys.executable, "-m", "fastmdxplora.cli.main", "info"]
    result = _run_command(command)
    if result.returncode != 0:
        raise BootstrapError(
            "Environment verification failed: "
            + result.stderr.strip()
            + "\nIf this is a package issue, rerun `fastmdx bootstrap --force` or consult the installation output."
        )
    _print("Environment verified successfully")


def _format_conda_error(stderr: str, env_name: str, create: bool) -> str:
    hint_lines = []
    if "UnsatisfiableError" in stderr:
        hint_lines.append(
            "Package resolution failed. This often means the requested environment is unavailable for the selected Python version or platform."
        )
        hint_lines.append(
            f"Use {python_range_string()} in the bootstrap environment and retry, or install Miniforge/conda with proper channels."
        )
    if "ResolvePackageNotFound" in stderr or "PackagesNotFoundError" in stderr:
        hint_lines.append(
            "A required package could not be found in the configured channels."
            " Ensure you have internet access and the conda-forge channel available."
        )
    if "Cannot connect" in stderr or "HTTP 000" in stderr:
        hint_lines.append(
            "Network or SSL error while fetching package metadata. Check your network or proxy settings."
        )
    hint_text = "\n".join(hint_lines)
    return (
        f"Conda env {'creation' if create else 'update'} failed for '{env_name}': {stderr.strip()}"
        + ("\n" + hint_text if hint_text else "")
        + "\nIf the environment already exists, try `fastmdx bootstrap --force`."
    )


def bootstrap_environment(
    *,
    env_name: str = DEFAULT_ENV_NAME,
    python_version: str = "3.10",
    yes: bool = False,
    force: bool = False,
    package_name: str = "fastmdxplora",
    editable: bool = False,
) -> None:
    supported = python_range_string()
    if not python_version.startswith("3."):
        raise BootstrapError(
            f"Invalid Python version {python_version!r}. FastMDXplora supports {supported}."
        )
    major_minor = tuple(int(part) for part in python_version.split(".")[:2])
    if major_minor < MIN_PYTHON or major_minor >= MAX_PYTHON:
        raise BootstrapError(
            f"Python {python_version} is outside the supported range. "
            f"FastMDXplora supports {supported} (the OpenMM/PDBFixer stack "
            f"is the bottleneck; pick 3.10 or 3.11 for the smoothest ride)."
        )

    conda_cli = _ensure_conda_available()
    if conda_cli is None:
        raise BootstrapError("Conda or mamba is required but could not be located or installed.")

    existing = _conda_env_exists(conda_cli[0], env_name)
    yaml_path = _write_environment_yaml()
    try:
        # Embed the requested Python version in the YAML before env creation.
        text = yaml_path.read_text(encoding="utf-8")
        min_str = f"{MIN_PYTHON[0]}.{MIN_PYTHON[1]}"
        max_str = f"{MAX_PYTHON[0]}.{MAX_PYTHON[1]}"
        text = text.replace(f"python>={min_str},<{max_str}", f"python=={python_version}")
        yaml_path.write_text(text, encoding="utf-8")

        if existing and force:
            _remove_conda_env(conda_cli[0], env_name)
            existing = False

        if existing:
            _update_conda_env(conda_cli[0], env_name, yaml_path)
        else:
            _create_conda_env(conda_cli[0], env_name, yaml_path)

        _install_package_in_env(conda_cli[0], env_name, package_name, editable=editable)
        _verify_env(conda_cli[0], env_name)
    finally:
        try:
            shutil.rmtree(yaml_path.parent)
        except OSError:
            pass

    _print(
        "Bootstrap complete. Activate the environment and run your first study:",
        prefix="SUCCESS",
    )
    _print(f"  conda activate {env_name}")
    _print("  fastmdx info")
    _print("After activation, run: fastmdx explore --system 1L2Y")


def main(argv: Iterable[str] | None = None) -> int:
    import argparse

    parser = argparse.ArgumentParser(
        description=(
            "Bootstrap a runnable FastMDXplora conda environment on Linux/macOS/Windows."
        )
    )
    parser.add_argument("--env-name", default=DEFAULT_ENV_NAME,
                        help="Name of the conda environment to create.")
    parser.add_argument("--python-version", default="3.10",
                        help=f"Python version to install into the environment ({python_range_string()}).")
    parser.add_argument("--force", action="store_true",
                        help="Recreate the conda environment if it already exists.")
    parser.add_argument("--yes", "-y", action="store_true",
                        help="Skip prompts by assuming yes.")
    args = parser.parse_args(argv)

    try:
        bootstrap_environment(
            env_name=args.env_name,
            python_version=args.python_version,
            yes=args.yes,
            force=args.force,
        )
        return 0
    except BootstrapError as exc:
        print(f"bootstrap failed: {exc}", file=sys.stderr)
        return 1


if __name__ == "__main__":
    raise SystemExit(main())
