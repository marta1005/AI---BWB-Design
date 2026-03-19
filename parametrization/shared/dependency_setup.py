import importlib
import shutil
import subprocess
import sys
import sysconfig
from pathlib import Path
from typing import Iterable, Optional

from .pyspline_shim import patch_pyspline_for_pygeo


def project_dependency_root() -> Path:
    package_dir = Path(__file__).resolve().parent
    return package_dir.parent.parent / "project"


def ensure_local_dependency_paths() -> None:
    project_dir = project_dependency_root()
    dependency_paths = (
        project_dir,
        project_dir / "pyspline",
        project_dir / "pygeo",
    )
    for dependency_dir in dependency_paths:
        dependency_path = str(dependency_dir)
        if dependency_dir.is_dir() and dependency_path not in sys.path:
            sys.path.insert(0, dependency_path)


def _clear_modules(prefixes: Iterable[str]) -> None:
    for prefix in prefixes:
        for name in list(sys.modules):
            if name == prefix or name.startswith(prefix + "."):
                sys.modules.pop(name, None)


def _candidate_python_configs() -> list[str]:
    version_tag = f"{sys.version_info.major}.{sys.version_info.minor}"
    executable = Path(sys.executable).resolve()
    bindir = executable.parent
    sysconfig_bindir = sysconfig.get_config_var("BINDIR")
    candidates: list[str] = []
    names = (
        f"python{version_tag}-config",
        "python3-config",
        "python-config",
    )

    for base_dir in (bindir, Path(sysconfig_bindir) if sysconfig_bindir else None):
        if base_dir is None:
            continue
        for name in names:
            candidate = base_dir / name
            if candidate.exists():
                candidates.append(str(candidate))

    for name in names:
        resolved = shutil.which(name)
        if resolved:
            candidates.append(resolved)

    deduped: list[str] = []
    for item in candidates:
        if item not in deduped:
            deduped.append(item)
    return deduped


def _detect_linker_flags(python_config: Optional[str]) -> str:
    if python_config:
        for extra_args in (["--embed", "--ldflags"], ["--ldflags"]):
            try:
                result = subprocess.run(
                    [python_config, *extra_args],
                    check=True,
                    capture_output=True,
                    text=True,
                )
            except Exception:
                continue
            flags = result.stdout.strip()
            if flags:
                return flags

    flags = []
    libdir = sysconfig.get_config_var("LIBDIR")
    ldlibrary = sysconfig.get_config_var("LDLIBRARY")
    libs = sysconfig.get_config_var("LIBS")
    syslibs = sysconfig.get_config_var("SYSLIBS")
    if libdir:
        flags.append(f"-L{libdir}")
    if ldlibrary:
        libname = str(ldlibrary)
        if libname.startswith("lib") and "." in libname:
            flags.append(f"-l{libname[3:].split('.', 1)[0]}")
    if libs:
        flags.append(str(libs))
    if syslibs:
        flags.append(str(syslibs))
    return " ".join(flags).strip()


def rebuild_local_pyspline() -> None:
    ensure_local_dependency_paths()
    project_dir = project_dependency_root()
    pyspline_dir = project_dir / "pyspline"
    python_config = next(iter(_candidate_python_configs()), None)
    f2py_bin = shutil.which("f2py") or str(Path(sys.executable).resolve().parent / "f2py")
    fortran_bin = shutil.which("gfortran")
    if not fortran_bin:
        homebrew_gfortran = Path("/opt/homebrew/bin/gfortran")
        if homebrew_gfortran.exists():
            fortran_bin = str(homebrew_gfortran)

    make_cmd = [
        "make",
        "-B",
        "module",
        f"PYTHON={sys.executable}",
        f"F2PY={f2py_bin}",
    ]
    if python_config:
        make_cmd.append(f"PYTHON-CONFIG={python_config}")
    if fortran_bin:
        make_cmd.append(f"FF90={fortran_bin}")

    linker_flags = _detect_linker_flags(python_config)
    if linker_flags:
        make_cmd.append(f"LINKER_FLAGS={linker_flags}")

    clean = subprocess.run(
        ["make", "clean"],
        cwd=pyspline_dir,
        capture_output=True,
        text=True,
    )
    if clean.returncode != 0:
        raise RuntimeError(
            "Automatic pyspline rebuild failed during 'make clean'. "
            f"stdout:\n{clean.stdout}\nstderr:\n{clean.stderr}"
        )

    build = subprocess.run(
        make_cmd,
        cwd=pyspline_dir,
        capture_output=True,
        text=True,
    )
    if build.returncode != 0:
        tail = "\n".join((build.stdout + "\n" + build.stderr).splitlines()[-80:])
        raise RuntimeError(
            "Automatic pyspline rebuild failed for the current Python environment. "
            "Last build output:\n"
            f"{tail}"
        )


def load_pyspline_curve(rebuild_if_needed: bool = True):
    ensure_local_dependency_paths()
    last_error: Optional[Exception] = None

    for attempt in range(2 if rebuild_if_needed else 1):
        try:
            patch_pyspline_for_pygeo()
            try:
                from pyspline import Curve
            except Exception:
                pyspline_pkg = importlib.import_module("pyspline.pyspline")
                sys.modules["pyspline"] = pyspline_pkg
                Curve = pyspline_pkg.Curve
            return Curve
        except Exception as exc:
            last_error = exc
            if attempt == 0 and rebuild_if_needed:
                _clear_modules(("pyspline", "pygeo"))
                rebuild_local_pyspline()
                _clear_modules(("pyspline", "pygeo"))
                continue
            break

    raise RuntimeError(
        "pyspline could not be imported for the current Python environment. "
        f"Detail: {last_error}"
    ) from last_error


def load_pygeo_class(rebuild_if_needed: bool = True):
    ensure_local_dependency_paths()
    load_pyspline_curve(rebuild_if_needed=rebuild_if_needed)
    last_error: Optional[Exception] = None

    for attempt in range(2 if rebuild_if_needed else 1):
        try:
            patch_pyspline_for_pygeo()
            try:
                from pygeo import pyGeo
            except Exception:
                pygeo_pkg = importlib.import_module("pygeo.pygeo")
                sys.modules["pygeo"] = pygeo_pkg
                pyGeo = pygeo_pkg.pyGeo
            return pyGeo
        except Exception as exc:
            last_error = exc
            if attempt == 0 and rebuild_if_needed:
                _clear_modules(("pyspline", "pygeo"))
                rebuild_local_pyspline()
                _clear_modules(("pyspline", "pygeo"))
                continue
            break

    raise RuntimeError(
        "pyGeo could not be imported for the current Python environment. "
        f"Detail: {last_error}"
    ) from last_error
