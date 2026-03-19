import importlib
import sys
import types


def patch_pyspline_for_pygeo() -> None:
    try:
        surface_module = importlib.import_module("pyspline.pySurface")
        curve_module = importlib.import_module("pyspline.pyCurve")
        volume_module = importlib.import_module("pyspline.pyVolume")
        utils_module = importlib.import_module("pyspline.utils")
        shim = importlib.import_module("pyspline")
    except Exception:
        surface_module = importlib.import_module("pyspline.pyspline.pySurface")
        curve_module = importlib.import_module("pyspline.pyspline.pyCurve")
        volume_module = importlib.import_module("pyspline.pyspline.pyVolume")
        utils_module = importlib.import_module("pyspline.pyspline.utils")
        shim = importlib.import_module("pyspline.pyspline")
        sys.modules["pyspline"] = shim
        sys.modules["pyspline.pySurface"] = surface_module
        sys.modules["pyspline.pyCurve"] = curve_module
        sys.modules["pyspline.pyVolume"] = volume_module
        sys.modules["pyspline.utils"] = utils_module

    surface = surface_module.Surface
    curve = curve_module.Curve
    volume = volume_module.Volume

    if not isinstance(shim, types.ModuleType):
        shim = types.ModuleType("pyspline")
        sys.modules["pyspline"] = shim

    if not hasattr(shim, "Surface"):
        shim.Surface = surface
    if not hasattr(shim, "Curve"):
        shim.Curve = curve
    if not hasattr(shim, "Volume"):
        shim.Volume = volume
