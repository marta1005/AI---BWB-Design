# Shared Utilities

`parametrization/shared` contains utilities that are reusable across geometry
families.

This package is intentionally small and generic:

- `airfoil_io.py`: `.dat` airfoil writing helpers
- `cst.py`: CST/Kulfan evaluation utilities and cosine spacing
- `dependency_setup.py`: local `pyspline` / `pyGeo` loading and rebuild helpers
- `pyspline_shim.py`: import compatibility shim for the local `pyspline` layout

If a module is useful for both `aircraft` and `bwb`, it belongs here.
