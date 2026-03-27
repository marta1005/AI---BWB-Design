"""Reusable utilities shared across parametrization packages."""

from .airfoil_fit import (
    AirfoilSectionFitResult,
    CSTAirfoilFitOptions,
    NormalizedAirfoilSection,
    fit_airfoil_section_cst,
    load_xyz_sections_by_span,
    normalize_airfoil_section,
    split_airfoil_upper_lower,
)
from .airfoil_io import write_airfoil_dat
from .cst import CST, KulfanCSTAirfoil, bernstein_matrix, cosine_spacing, fit_kulfan_airfoil_coefficients
from .dependency_setup import (
    ensure_local_dependency_paths,
    load_pygeo_class,
    load_pyspline_curve,
    project_dependency_root,
    rebuild_local_pyspline,
)
from .pyspline_shim import patch_pyspline_for_pygeo

__all__ = [
    "AirfoilSectionFitResult",
    "CST",
    "CSTAirfoilFitOptions",
    "KulfanCSTAirfoil",
    "NormalizedAirfoilSection",
    "bernstein_matrix",
    "cosine_spacing",
    "fit_airfoil_section_cst",
    "fit_kulfan_airfoil_coefficients",
    "ensure_local_dependency_paths",
    "load_pygeo_class",
    "load_pyspline_curve",
    "load_xyz_sections_by_span",
    "normalize_airfoil_section",
    "patch_pyspline_for_pygeo",
    "project_dependency_root",
    "rebuild_local_pyspline",
    "split_airfoil_upper_lower",
    "write_airfoil_dat",
]
