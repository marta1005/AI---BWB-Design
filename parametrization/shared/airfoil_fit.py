from __future__ import annotations

from dataclasses import dataclass
from math import atan2, degrees
from pathlib import Path
from typing import Dict

import numpy as np
from scipy.interpolate import PchipInterpolator

from .cst import CST, bernstein_matrix, cosine_spacing


@dataclass(frozen=True)
class CSTAirfoilFitOptions:
    degree: int = 5
    n1: float = 0.5
    n2: float = 1.0
    shared_leading_edge: bool = False
    smoothness_weight: float = 0.0
    fit_xmin: float = 1.0e-4
    fit_xmax: float = 0.99
    sample_count: int = 301


@dataclass(frozen=True)
class NormalizedAirfoilSection:
    upper_x: np.ndarray
    upper_y: np.ndarray
    lower_x: np.ndarray
    lower_y: np.ndarray
    upper_zero_y: np.ndarray
    lower_zero_y: np.ndarray
    x_fit: np.ndarray
    y_upper: np.ndarray
    y_lower: np.ndarray
    upper_scale: float
    lower_scale: float
    upper_tail_y: float
    lower_tail_y: float
    chord: float
    twist_deg: float
    te_thickness: float
    te_mid_x: float
    te_mid_z: float
    le_x: float
    le_z: float


@dataclass(frozen=True)
class AirfoilSectionFitResult:
    x_fit: np.ndarray
    y_upper: np.ndarray
    y_lower: np.ndarray
    fit_y_upper: np.ndarray
    fit_y_lower: np.ndarray
    upper_cst: np.ndarray
    lower_cst: np.ndarray
    upper_scale: float
    lower_scale: float
    upper_tail_y: float
    lower_tail_y: float
    chord: float
    twist_deg: float
    te_thickness: float
    te_mid_x: float
    te_mid_z: float
    le_x: float
    le_z: float
    fit_xmin: float
    fit_xmax: float
    rmse_upper: float
    rmse_lower: float
    max_abs_error: float


def load_xyz_sections_by_span(path: str | Path) -> Dict[float, np.ndarray]:
    sections: Dict[float, list[list[float]]] = {}
    for raw_line in Path(path).read_text(encoding="utf-8").splitlines():
        parts = raw_line.split()
        if len(parts) != 3:
            continue
        x_val, y_val, z_val = (float(parts[0]), float(parts[1]), float(parts[2]))
        sections.setdefault(y_val, []).append([x_val, z_val])
    return {yy: np.asarray(points, dtype=float) for yy, points in sections.items()}


def split_airfoil_upper_lower(section: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
    le_idx = int(np.argmin(section[:, 0]))
    upper = section[: le_idx + 1]
    lower = section[le_idx:]
    return upper, lower


def _interp_shape_preserving(x_src: np.ndarray, y_src: np.ndarray, x_dst: np.ndarray) -> np.ndarray:
    x_src = np.asarray(x_src, dtype=float)
    y_src = np.asarray(y_src, dtype=float)
    x_dst = np.asarray(x_dst, dtype=float)

    unique_x, unique_idx = np.unique(x_src, return_index=True)
    unique_y = y_src[unique_idx]
    if unique_x.size < 2:
        return np.full_like(x_dst, float(unique_y[0]), dtype=float)

    interpolant = PchipInterpolator(unique_x, unique_y)
    return np.asarray(interpolant(x_dst), dtype=float)


def _fit_surface_cst_coefficients(
    x: np.ndarray,
    y: np.ndarray,
    degree: int,
    n1: float,
    n2: float,
    regularization: float = 1.0e-8,
    smoothness_weight: float = 0.0,
) -> np.ndarray:
    x = np.asarray(x, dtype=float)
    y = np.asarray(y, dtype=float)
    if x.ndim != 1 or y.ndim != 1 or x.size != y.size:
        raise ValueError("x and y must be 1D arrays with the same length")

    class_fun = (x**float(n1)) * ((1.0 - x) ** float(n2))
    basis = bernstein_matrix(x, degree)
    design = class_fun[:, None] * basis

    normal_matrix = design.T @ design + float(regularization) * np.eye(degree + 1)
    if smoothness_weight > 0.0 and degree >= 2:
        d2 = np.zeros((degree - 1, degree + 1), dtype=float)
        for row in range(degree - 1):
            d2[row, row] = 1.0
            d2[row, row + 1] = -2.0
            d2[row, row + 2] = 1.0
        normal_matrix = normal_matrix + float(smoothness_weight) * (d2.T @ d2)

    return np.linalg.solve(normal_matrix, design.T @ y).astype(float)


def normalize_airfoil_section(section: np.ndarray, sample_count: int = 301) -> NormalizedAirfoilSection:
    upper_raw, lower_raw = split_airfoil_upper_lower(section)

    te_upper = np.asarray(section[0], dtype=float)
    te_lower = np.asarray(section[-1], dtype=float)
    te_mid = 0.5 * (te_upper + te_lower)
    le_point = np.asarray(section[int(np.argmin(section[:, 0]))], dtype=float)
    chord_vector = te_mid - le_point
    chord = float(np.linalg.norm(chord_vector))
    if chord <= 0.0:
        raise ValueError(f"invalid chord length {chord}")
    twist_deg = float(degrees(atan2(float(chord_vector[1]), float(chord_vector[0]))))

    ex = chord_vector / chord
    ez = np.array([-ex[1], ex[0]], dtype=float)
    upper_local = np.column_stack(((upper_raw - le_point) @ ex, (upper_raw - le_point) @ ez))
    lower_local = np.column_stack(((lower_raw - le_point) @ ex, (lower_raw - le_point) @ ez))

    # Upper arrives TE->LE in the source order, so reverse it to interpolate LE->TE.
    upper_local = upper_local[::-1]

    upper_scale = float(upper_local[-1, 0])
    lower_scale = float(lower_local[-1, 0])
    if upper_scale <= 0.0 or lower_scale <= 0.0:
        raise ValueError(f"invalid surface scales {(upper_scale, lower_scale)}")

    # Match cst-modeling's normalize_foil(): each surface is scaled by its own
    # trailing-edge x-coordinate after the de-twisting rotation.
    upper_x = upper_local[:, 0] / upper_scale
    upper_z = upper_local[:, 1] / upper_scale
    lower_x = lower_local[:, 0] / lower_scale
    lower_z = lower_local[:, 1] / lower_scale

    upper_tail_y = float(upper_z[-1])
    lower_tail_y = float(lower_z[-1])
    upper_zero_y = upper_z - upper_x * upper_tail_y
    lower_zero_y = lower_z - lower_x * lower_tail_y

    x_fit = cosine_spacing(sample_count)
    # Compare and fit in the same per-surface normalized frame used by the
    # original cst-modeling code path.
    y_upper = _interp_shape_preserving(upper_x, upper_z, x_fit)
    y_lower = _interp_shape_preserving(lower_x, lower_z, x_fit)
    te_thickness = float(upper_tail_y - lower_tail_y)

    return NormalizedAirfoilSection(
        upper_x=upper_x,
        upper_y=upper_z,
        lower_x=lower_x,
        lower_y=lower_z,
        upper_zero_y=upper_zero_y,
        lower_zero_y=lower_zero_y,
        x_fit=x_fit,
        y_upper=y_upper,
        y_lower=y_lower,
        upper_scale=upper_scale,
        lower_scale=lower_scale,
        upper_tail_y=upper_tail_y,
        lower_tail_y=lower_tail_y,
        chord=chord,
        twist_deg=twist_deg,
        te_thickness=te_thickness,
        te_mid_x=float(te_mid[0]),
        te_mid_z=float(te_mid[1]),
        le_x=float(le_point[0]),
        le_z=float(le_point[1]),
    )


def fit_airfoil_section_cst(
    section: np.ndarray,
    options: CSTAirfoilFitOptions | None = None,
) -> AirfoilSectionFitResult:
    fit_options = options or CSTAirfoilFitOptions()
    normalized = normalize_airfoil_section(section, sample_count=fit_options.sample_count)

    x_fit = np.asarray(normalized.x_fit, dtype=float)
    upper_x = np.asarray(normalized.upper_x, dtype=float)
    upper_y = np.asarray(normalized.upper_y, dtype=float)
    lower_x = np.asarray(normalized.lower_x, dtype=float)
    lower_y = np.asarray(normalized.lower_y, dtype=float)
    upper_zero_y = np.asarray(normalized.upper_zero_y, dtype=float)
    lower_zero_y = np.asarray(normalized.lower_zero_y, dtype=float)
    y_upper = np.asarray(normalized.y_upper, dtype=float)
    y_lower = np.asarray(normalized.y_lower, dtype=float)
    upper_tail_y = float(normalized.upper_tail_y)
    lower_tail_y = float(normalized.lower_tail_y)
    te_thickness = float(normalized.te_thickness)

    fit_mask_upper = (upper_x >= float(fit_options.fit_xmin)) & (upper_x < float(fit_options.fit_xmax))
    fit_mask_lower = (lower_x >= float(fit_options.fit_xmin)) & (lower_x < float(fit_options.fit_xmax))
    if np.count_nonzero(fit_mask_upper) < fit_options.degree + 1:
        raise ValueError(
            f"not enough fit samples for degree {fit_options.degree}: "
            f"{np.count_nonzero(fit_mask_upper)} upper samples available"
        )
    if np.count_nonzero(fit_mask_lower) < fit_options.degree + 1:
        raise ValueError(
            f"not enough fit samples for degree {fit_options.degree}: "
            f"{np.count_nonzero(fit_mask_lower)} lower samples available"
        )

    upper_cst = _fit_surface_cst_coefficients(
        x=upper_x[fit_mask_upper],
        y=upper_zero_y[fit_mask_upper],
        degree=fit_options.degree,
        n1=fit_options.n1,
        n2=fit_options.n2,
        smoothness_weight=fit_options.smoothness_weight,
    )
    lower_cst = _fit_surface_cst_coefficients(
        x=lower_x[fit_mask_lower],
        y=lower_zero_y[fit_mask_lower],
        degree=fit_options.degree,
        n1=fit_options.n1,
        n2=fit_options.n2,
        smoothness_weight=fit_options.smoothness_weight,
    )

    fit_y_upper = CST(upper_cst, n1=fit_options.n1, n2=fit_options.n2, delta_te=0.0).evaluate(x_fit)
    fit_y_upper = fit_y_upper + x_fit * upper_tail_y
    fit_y_lower = CST(lower_cst, n1=fit_options.n1, n2=fit_options.n2, delta_te=0.0).evaluate(x_fit)
    fit_y_lower = fit_y_lower + x_fit * lower_tail_y

    error_upper = fit_y_upper - y_upper
    error_lower = fit_y_lower - y_lower

    return AirfoilSectionFitResult(
        x_fit=x_fit,
        y_upper=y_upper,
        y_lower=y_lower,
        fit_y_upper=fit_y_upper,
        fit_y_lower=fit_y_lower,
        upper_cst=upper_cst,
        lower_cst=lower_cst,
        upper_scale=float(normalized.upper_scale),
        lower_scale=float(normalized.lower_scale),
        upper_tail_y=upper_tail_y,
        lower_tail_y=lower_tail_y,
        chord=float(normalized.chord),
        twist_deg=float(normalized.twist_deg),
        te_thickness=te_thickness,
        te_mid_x=float(normalized.te_mid_x),
        te_mid_z=float(normalized.te_mid_z),
        le_x=float(normalized.le_x),
        le_z=float(normalized.le_z),
        fit_xmin=float(fit_options.fit_xmin),
        fit_xmax=float(fit_options.fit_xmax),
        rmse_upper=float(np.sqrt(np.mean(error_upper**2))),
        rmse_lower=float(np.sqrt(np.mean(error_lower**2))),
        max_abs_error=float(max(np.max(np.abs(error_upper)), np.max(np.abs(error_lower)))),
    )
